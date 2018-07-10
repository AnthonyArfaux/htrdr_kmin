/* Copyright (C) 2018 Université Paul Sabatier, |Meso|Star>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>. */

#include "htrdr.h"
#include "htrdr_c.h"
#include "htrdr_sun.h"

#include <rsys/algorithm.h>
#include <rsys/double33.h>
#include <rsys/dynamic_array_double.h>
#include <rsys/ref_count.h>
#include <rsys/math.h>

#include <star/ssp.h>

struct htrdr_sun {
  /* Short wave radiance in W.m^-2.sr^-1, for each spectral interval */
  struct darray_double radiances_sw;

  /* Short wave spectral interval boundaries, in cm^-1 */
  struct darray_double wavenumbers_sw;

  double half_angle; /* In radian */
  double cos_half_angle;
  double solid_angle; /* In sr; solid_angle = 2*PI*(1 - cos(half_angle)) */
  double frame[9];

  ref_T ref;
  struct htrdr* htrdr;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static INLINE int
cmp_dbl(const void* a, const void* b)
{
  const double d0 = *((const double*)a);
  const double d1 = *((const double*)b);
  return d0 < d1 ? -1 : (d0 > d1 ? 1 : 0);
}

static void
release_sun(ref_T* ref)
{
  struct htrdr_sun* sun;
  ASSERT(ref);
  sun = CONTAINER_OF(ref, struct htrdr_sun, ref);
  darray_double_release(&sun->radiances_sw);
  darray_double_release(&sun->wavenumbers_sw);
  MEM_RM(sun->htrdr->allocator, sun);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_sun_create(struct htrdr* htrdr, struct htrdr_sun** out_sun)
{
  const double incoming_flux_sw[] = { /* In W.m^-2 */
    12.793835026999544, 12.109561093845551, 20.365091338928245,
    23.729742422870157, 22.427697221814142, 55.626612361454150,
    102.93146523363953, 24.293596268358986, 345.73659325842243,
    218.18441435866691, 347.18437832794524, 129.49426803812202,
    50.146977730963876, 3.1197193425713365
  };
  const double wavenumbers_sw[] = { /* In cm^-1 */
    820.000, 2600.00, 3250.00, 4000.00, 4650.00,
    5150.00, 6150.00, 7700.00, 8050.00, 12850.0,
    16000.0, 22650.0, 29000.0, 38000.0, 49999.0
  };

  const size_t nspectral_intervals = sizeof(incoming_flux_sw)/sizeof(double);
  const double main_dir[3] = {0, 0, 1}; /* Default main sun direction */
  struct htrdr_sun* sun = NULL;
  size_t i;
  res_T res = RES_OK;
  ASSERT(htrdr && out_sun);
  ASSERT(sizeof(wavenumbers_sw)/sizeof(double) == nspectral_intervals+1);

  sun = MEM_CALLOC(htrdr->allocator, 1, sizeof(*sun));
  if(!sun) {
    htrdr_log_err(htrdr, "could not allocate the sun data structure.\n");
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&sun->ref);
  sun->htrdr = htrdr;
  darray_double_init(htrdr->allocator, &sun->radiances_sw);
  darray_double_init(htrdr->allocator, &sun->wavenumbers_sw);
  sun->half_angle = 4.6524e3;
  sun->cos_half_angle = cos(sun->half_angle);
  sun->solid_angle = 2*PI*(1-sun->cos_half_angle);
  d33_basis(sun->frame, main_dir);

  res = darray_double_resize(&sun->radiances_sw, nspectral_intervals);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "could not allocate the list of per spectral band radiance of the sun.\n");
    goto error;
  }
  res = darray_double_resize(&sun->wavenumbers_sw, nspectral_intervals+1);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "could not allocate the list of spectral band boundaries of the sun.\n");
    goto error;
  }

  FOR_EACH(i, 0, darray_double_size_get(&sun->radiances_sw)) {
    /* Convert the incoming flux in radiance */
    darray_double_data_get(&sun->radiances_sw)[i] = incoming_flux_sw[i] / PI;
  }
  FOR_EACH(i, 0, darray_double_size_get(&sun->wavenumbers_sw)) {
    darray_double_data_get(&sun->wavenumbers_sw)[i] = wavenumbers_sw[i];
  }

exit:
  *out_sun = sun;
  return res;
error:
  if(sun) {
    htrdr_sun_ref_put(sun);
    sun = NULL;
  }
  goto exit;
}

void
htrdr_sun_ref_get(struct htrdr_sun* sun)
{
  ASSERT(sun);
  ref_get(&sun->ref);
}

void
htrdr_sun_ref_put(struct htrdr_sun* sun)
{
  ASSERT(sun);
  ref_put(&sun->ref, release_sun);
}

void
htrdr_sun_set_direction(struct htrdr_sun* sun, const double dir[3])
{
  ASSERT(sun && dir && d3_is_normalized(dir));
  d33_basis(sun->frame, dir);
}

double*
htrdr_sun_sample_direction
  (struct htrdr_sun* sun,
   struct ssp_rng* rng,
   double dir[3])
{
  ASSERT(sun && rng && dir);
  ssp_ran_sphere_cap_uniform_local(rng, sun->cos_half_angle, dir, NULL);
  return d33_muld3(dir, sun->frame, dir);
}

double
htrdr_sun_get_solid_angle(const struct htrdr_sun* sun)
{
  ASSERT(sun);
  return sun->solid_angle;
}

double
htrdr_sun_get_radiance(const struct htrdr_sun* sun, const double wavelength)
{
  const double* wavenumbers;
  const double wnum = wavelength_to_wavenumber(wavelength);
  const double* wnum_upp;
  size_t nwavenumbers;
  size_t ispectral_band;
  ASSERT(sun && wavelength > 0);

  wavenumbers = darray_double_cdata_get(&sun->wavenumbers_sw);
  nwavenumbers = darray_double_size_get(&sun->wavenumbers_sw);
  ASSERT(nwavenumbers);

  if(wnum < wavenumbers[0] || wnum > wavenumbers[nwavenumbers-1]) {
    htrdr_log_warn(sun->htrdr,
      "the submitted wavelength is outside the sun spectrum.\n");
  }

  wnum_upp = search_lower_bound
    (&wnum, wavenumbers, nwavenumbers, sizeof(double), cmp_dbl);

  if(!wnum_upp) { /* Clamp to the upper spectral band */
    ispectral_band = nwavenumbers - 2;
    ASSERT(ispectral_band == darray_double_size_get(&sun->radiances_sw)-1);
  } else if(wnum_upp == wavenumbers) { /* Clamp to the lower spectral band */
    ispectral_band = 0;
  } else {
    ispectral_band = (size_t)(wnum_upp - wavenumbers - 1);
  }
  return darray_double_cdata_get(&sun->radiances_sw)[ispectral_band];
}

int
htrdr_sun_is_dir_in_solar_cone(const struct htrdr_sun* sun, const double dir[3])
{
  const double* main_dir;
  double dot;
  ASSERT(sun && dir && d3_is_normalized(dir));
  ASSERT(d3_is_normalized(sun->frame + 6));
  main_dir = sun->frame + 6;
  dot = d3_dot(dir, main_dir);
  return dot >= sun->cos_half_angle;
}

