/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2023 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2023 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2023 Université Paul Sabatier
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

#include "atmosphere/htrdr_atmosphere_c.h"
#include "atmosphere/htrdr_atmosphere_sun.h"

#include "core/htrdr.h"
#include "core/htrdr_log.h"

#include <rsys/algorithm.h>
#include <rsys/double33.h>
#include <rsys/ref_count.h>
#include <rsys/math.h>

#include <star/ssp.h>

struct htrdr_atmosphere_sun {
  double half_angle; /* In radian */
  double cos_half_angle;
  double solid_angle; /* In sr; solid_angle = 2*PI*(1 - cos(half_angle)) */
  double frame[9];
  double temperature; /* In K */

  ref_T ref;
  struct htrdr* htrdr;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
release_sun(ref_T* ref)
{
  struct htrdr_atmosphere_sun* sun;
  struct htrdr* htrdr;
  ASSERT(ref);
  sun = CONTAINER_OF(ref, struct htrdr_atmosphere_sun, ref);
  htrdr = sun->htrdr;
  MEM_RM(htrdr_get_allocator(htrdr), sun);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_atmosphere_sun_create
  (struct htrdr* htrdr,
   struct htrdr_atmosphere_sun** out_sun)
{
  const double main_dir[3] = {0, 0, 1}; /* Default main sun direction */
  struct htrdr_atmosphere_sun* sun = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && out_sun);

  sun = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*sun));
  if(!sun) {
    htrdr_log_err(htrdr, "could not allocate the sun data structure.\n");
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&sun->ref);

  sun->half_angle = 4.6524e-3;
  sun->temperature = 5778;
  sun->cos_half_angle = cos(sun->half_angle);
  sun->solid_angle = 2*PI*(1-sun->cos_half_angle);
  d33_basis(sun->frame, main_dir);
  htrdr_ref_get(htrdr);
  sun->htrdr = htrdr;

exit:
  *out_sun = sun;
  return res;
error:
  if(sun) {
    htrdr_atmosphere_sun_ref_put(sun);
    sun = NULL;
  }
  goto exit;
}

void
htrdr_atmosphere_sun_ref_get(struct htrdr_atmosphere_sun* sun)
{
  ASSERT(sun);
  ref_get(&sun->ref);
}

void
htrdr_atmosphere_sun_ref_put(struct htrdr_atmosphere_sun* sun)
{
  ASSERT(sun);
  ref_put(&sun->ref, release_sun);
}

void
htrdr_atmosphere_sun_set_direction
  (struct htrdr_atmosphere_sun* sun,
   const double dir[3])
{
  ASSERT(sun && dir && d3_is_normalized(dir));
  d33_basis(sun->frame, dir);
}

double
htrdr_atmosphere_sun_sample_direction
  (struct htrdr_atmosphere_sun* sun,
   struct ssp_rng* rng,
   double dir[3])
{
  ASSERT(sun && rng && dir);
  ssp_ran_sphere_cap_uniform_local(rng, sun->cos_half_angle, dir, NULL);
  d33_muld3(dir, sun->frame, dir);
  return 1.0 / htrdr_atmosphere_sun_get_solid_angle(sun);
}

double
htrdr_atmosphere_sun_get_solid_angle(const struct htrdr_atmosphere_sun* sun)
{
  ASSERT(sun);
  return sun->solid_angle;
}

double
htrdr_atmosphere_sun_get_radiance
  (const struct htrdr_atmosphere_sun* sun,
   const double wlen/*In nm*/)
{
  return htrdr_planck_monochromatic
    (wlen*1.e-9/*From nm to m*/, sun->temperature);
}

int
htrdr_atmosphere_sun_is_dir_in_solar_cone
  (const struct htrdr_atmosphere_sun* sun,
   const double dir[3])
{
  const double* main_dir;
  double dot;
  ASSERT(sun && dir && d3_is_normalized(dir));
  ASSERT(d3_is_normalized(sun->frame + 6));
  main_dir = sun->frame + 6;
  dot = d3_dot(dir, main_dir);
  return dot >= sun->cos_half_angle;
}

