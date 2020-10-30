/* Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
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
#include "htrdr_interface.h"

#include <modradurb/mrumtl.h>

#include <star/s3d.h>
#include <star/ssf.h>
#include <star/ssp.h>

#include <rsys/cstr.h>
#include <rsys/double3.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static res_T
create_bsdf_diffuse
  (struct htrdr* htrdr,
   const struct mrumtl_brdf* brdf,
   const size_t ithread,
   struct ssf_bsdf** out_bsdf)
{
  struct ssf_bsdf* bsdf = NULL;
  double reflectivity = 0;
  res_T res = RES_OK;
  ASSERT(htrdr && brdf && out_bsdf);
  ASSERT(mrumtl_brdf_get_type(brdf) == MRUMTL_BRDF_LAMBERTIAN);
  ASSERT(ithread < htrdr->nthreads);

  res = ssf_bsdf_create
    (&htrdr->lifo_allocators[ithread], &ssf_lambertian_reflection, &bsdf);
  if(res != RES_OK) goto error;

  reflectivity = mrumtl_brdf_lambertian_get_reflectivity(brdf);
  res = ssf_lambertian_reflection_setup(bsdf, reflectivity);
  if(res != RES_OK) goto error;

exit:
  *out_bsdf = bsdf;
  return res;
error:
   if(bsdf) { SSF(bsdf_ref_put(bsdf)); bsdf = NULL; }
  goto exit;
}

static res_T
create_bsdf_specular
  (struct htrdr* htrdr,
   const struct mrumtl_brdf* brdf,
   const size_t ithread,
   struct ssf_bsdf** out_bsdf)
{
  struct ssf_bsdf* bsdf = NULL;
  struct ssf_fresnel* fresnel = NULL;
  double reflectivity = 0;
  res_T res = RES_OK;
  ASSERT(htrdr && brdf && out_bsdf);
  ASSERT(mrumtl_brdf_get_type(brdf) == MRUMTL_BRDF_SPECULAR);
  ASSERT(ithread < htrdr->nthreads);

  res = ssf_bsdf_create
    (&htrdr->lifo_allocators[ithread], &ssf_specular_reflection, &bsdf);
  if(res != RES_OK) goto error;

  res = ssf_fresnel_create
    (&htrdr->lifo_allocators[ithread], &ssf_fresnel_constant, &fresnel);
  if(res != RES_OK) goto error;

  reflectivity = mrumtl_brdf_specular_get_reflectivity(brdf);
  res =  ssf_fresnel_constant_setup(fresnel, reflectivity);
  if(res != RES_OK) goto error;

  res = ssf_specular_reflection_setup(bsdf, fresnel);
  if(res != RES_OK) goto error;

exit:
  if(fresnel) SSF(fresnel_ref_put(fresnel));
  *out_bsdf = bsdf;
  return res;
error:
  if(bsdf) { SSF(bsdf_ref_put(bsdf)); bsdf = NULL; }
  goto exit;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_interface_create_bsdf
  (struct htrdr* htrdr,
   const struct htrdr_interface* interf,
   const size_t ithread,
   const double wavelength,
   const double pos[3],
   const double dir[3],
   struct ssp_rng* rng,
   struct s3d_hit* hit,
   struct ssf_bsdf** out_bsdf)
{
  enum { FRONT, BACK };
  struct ssf_bsdf* bsdf = NULL;
  const struct mrumtl_brdf* brdf = NULL;
  const struct mrumtl* mat = NULL;
  double N[3];
  double r;
  int hit_side;
  res_T res = RES_OK;
  (void)pos;
  ASSERT(htrdr && pos && hit && out_bsdf);
  ASSERT(interf && (interf->mtl_front || interf->mtl_back));

  ASSERT(htrdr && interf && pos && dir && hit && out_bsdf);
  ASSERT(d3_is_normalized(dir));

  d3_normalize(N, d3_set_f3(N, hit->normal));

  hit_side = d3_dot(N, dir) < 0 ? FRONT : BACK;

  /* Retrieve the brdf of the material on the other side of the hit side */
  switch(hit_side) {
    case BACK: mat = interf->mtl_front; break;
    case FRONT: mat = interf->mtl_back; break;
    default: FATAL("Unreachable code.\n");  break;
  }

  /* Due to numerical issue the hit side might be wrong and thus the fetched
   * material might be undefined (e.g. semi-transparent materials). Handle this
   * issue by fetching the other material. */
  if(!mat) {
    switch(hit_side) {
      case BACK: mat = interf->mtl_back; break;
      case FRONT: mat = interf->mtl_front; break;
      default: FATAL("Unreachable code.\n");  break;
    }
  }
  ASSERT(mat);

  r = ssp_rng_canonical(rng);

  res = mrumtl_fetch_brdf2(mat, wavelength, r, &brdf);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: error retreiving the MruMtl BRDF for the wavelength %g.\n",
      FUNC_NAME, wavelength);
    res = RES_BAD_ARG;
    goto error;
  }

  switch(mrumtl_brdf_get_type(brdf)) {
    case MRUMTL_BRDF_LAMBERTIAN:
      res = create_bsdf_diffuse(htrdr, brdf, ithread, &bsdf);
      break;
    case MRUMTL_BRDF_SPECULAR:
      res = create_bsdf_specular(htrdr, brdf, ithread, &bsdf);
      break;
    default: FATAL("Unreachable code.\n");  break;
  }
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "%s: could not create the BSDF -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }

exit:
  *out_bsdf = bsdf;
  return res;
error:
  if(bsdf) { SSF(bsdf_ref_put(bsdf)); bsdf = NULL; }
  goto exit;
}

