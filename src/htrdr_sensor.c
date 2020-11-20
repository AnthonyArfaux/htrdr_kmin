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
#include "htrdr_camera.h"
#include "htrdr_ground.h"
#include "htrdr_interface.h"
#include "htrdr_rectangle.h"
#include "htrdr_sensor.h"

#include <star/s3d.h>
#include <star/ssp.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static res_T
sample_camera_ray
  (struct htrdr_camera* cam,
   const size_t ipix[2],
   const double pix_sz[2],
   struct ssp_rng* rng,
   double ray_org[3],
   double ray_dir[3])
{
  double pix_samp[2];
  ASSERT(cam && ipix && pix_sz && rng && ray_org && ray_dir);

  /* Sample a position into the pixel, in the normalized image plane */
  pix_samp[0] = ((double)ipix[0] + ssp_rng_canonical(rng)) * pix_sz[0];
  pix_samp[1] = ((double)ipix[1] + ssp_rng_canonical(rng)) * pix_sz[1];

  /* Generate a ray starting from the pinhole camera and passing through the
   * pixel sample */
  htrdr_camera_ray(cam, pix_samp, ray_org, ray_dir);

  return RES_OK;
}

static res_T
sample_rectangle_ray
  (struct htrdr_rectangle* rect,
   struct htrdr* htrdr,
   const size_t ipix[2],
   const double pix_sz[2],
   struct ssp_rng* rng,
   double ray_org[3],
   double ray_dir[3])
{
  struct s3d_hit hit = S3D_HIT_NULL;
  double pix_samp[2];
  const double up_dir[3] = {0,0,1};
  const double range[2] = {0, DBL_MAX};
  double normal[3];
  ASSERT(rect && htrdr && ipix && pix_sz && rng && ray_org && ray_dir);

  /* Sample a position into the pixel, in the normalized image plane */
  pix_samp[0] = ((double)ipix[0] + ssp_rng_canonical(rng)) * pix_sz[0];
  pix_samp[1] = ((double)ipix[1] + ssp_rng_canonical(rng)) * pix_sz[1];

  /* Retrieve the world space position of pix_samp */
  htrdr_rectangle_sample_pos(rect, pix_samp, ray_org);

  /* Check that `ray_org' is not included into a geometry */
  HTRDR(ground_trace_ray(htrdr->ground, ray_org, up_dir, range, NULL, &hit));

  /* Up direction is occluded. Check if the sample must be rejected, i.e. does it
   * lies inside a geometry? */
  if(!S3D_HIT_NONE(&hit)) {
    struct htrdr_interface interf = HTRDR_INTERFACE_NULL;
    const struct htrdr_mtl* mtl = NULL;
    float N[3]; /* Normalized normal of the hit */
    float wi[3];
    float cos_wi_N;

    /* Compute the cosine between the up direction and the hit normal */
    f3_set_d3(wi, up_dir);
    f3_normalize(N, hit.normal);
    cos_wi_N = f3_dot(wi, N);

    /* Fetch the hit interface and retrieve the material into which the ray was
     * traced */
    htrdr_ground_get_interface(htrdr->ground, &hit, &interf);
    mtl = cos_wi_N < 0 ? &interf.mtl_front : &interf.mtl_back;

    /* Reject the sample if the incident direction do not travel into the sky */
    if(strcmp(mtl->name, htrdr->sky_mtl_name) != 0) return RES_BAD_OP;
  }

  /* Sample a ray direction */
  htrdr_rectangle_get_normal(rect, normal);
  ssp_ran_hemisphere_cos(rng, normal, ray_dir, NULL);

  return RES_OK;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_sensor_sample_primary_ray
  (const struct htrdr_sensor* sensor,
   struct htrdr* htrdr,
   const size_t ipix[2],
   const double pix_sz[2],
   struct ssp_rng* rng,
   double ray_org[3],
   double ray_dir[3])
{
  res_T res = RES_OK;
  switch(sensor->type) {
    case HTRDR_SENSOR_CAMERA:
      res = sample_camera_ray(sensor->camera, ipix, pix_sz, rng, ray_org, ray_dir);
      break;
    case HTRDR_SENSOR_RECTANGLE:
      res = sample_rectangle_ray(sensor->rectangle, htrdr, ipix,
        pix_sz, rng, ray_org, ray_dir);
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
  return res;
}

