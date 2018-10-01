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
#include "htrdr_camera.h"

#include <rsys/double3.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

struct htrdr_camera {
  /* Orthogonal basis of the camera */
  double axis_x[3];
  double axis_y[3];
  double axis_z[3];

  double position[3];

  ref_T ref;
  struct htrdr* htrdr;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
camera_release(ref_T* ref)
{
  struct htrdr_camera* cam;
  ASSERT(ref);
  cam = CONTAINER_OF(ref, struct htrdr_camera, ref);
  MEM_RM(cam->htrdr->allocator, cam);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_camera_create
  (struct htrdr* htrdr,
   const double position[3],
   const double target[3],
   const double up[3],
   const double proj_ratio,
   const double fov, /* In radian */
   struct htrdr_camera** out_cam)
{
  double x[3], y[3], z[3];
  double img_plane_depth;
  struct htrdr_camera* cam = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && position && target && up && out_cam);

  cam = MEM_CALLOC(htrdr->allocator, 1, sizeof(*cam));
  if(!cam) {
    htrdr_log_err(htrdr, "could not allocate the camera data structure.\n");
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&cam->ref);
  cam->htrdr = htrdr;

    if(fov <= 0) {
    htrdr_log_err(htrdr, "invalid horizontal camera field of view `%g'\n", fov);
    res = RES_BAD_ARG;
    goto error;
  }

  if(proj_ratio <= 0) {
    htrdr_log_err(htrdr, "invalid projection ratio `%g'\n", proj_ratio);
    res = RES_BAD_ARG;
    goto error;
  }

  if(d3_normalize(z, d3_sub(z, target, position)) <= 0
  || d3_normalize(x, d3_cross(x, z, up)) <= 0
  || d3_normalize(y, d3_cross(y, z, x)) <= 0) {
    htrdr_log_err(htrdr,
      "invalid camera point of view:\n"
      "  position = %g %g %g\n"
      "  target   = %g %g %g\n"
      "  up       = %g %g %g\n",
      SPLIT3(position), SPLIT3(target), SPLIT3(up));
    res = RES_BAD_ARG;
    goto error;
  }

  img_plane_depth = 1.0/tan(fov*0.5);
  d3_muld(cam->axis_x, x, proj_ratio);
  d3_set(cam->axis_y, y);
  d3_muld(cam->axis_z, z, img_plane_depth);
  d3_set(cam->position, position);

exit:
  *out_cam = cam;
  return res;
error:
  if(cam) {
    htrdr_camera_ref_put(cam);
    cam = NULL;
  }
  goto exit;
}

void
htrdr_camera_ref_get(struct htrdr_camera* cam)
{
  ASSERT(cam);
  ref_get(&cam->ref);
}

void
htrdr_camera_ref_put(struct htrdr_camera* cam)
{
  ASSERT(cam);
  ref_put(&cam->ref, camera_release);
}

void
htrdr_camera_ray
  (const struct htrdr_camera* cam,
   const double sample[2],
   double ray_org[3],
   double ray_dir[3])
{
  double x[3], y[3], len;
  (void)len;
  ASSERT(cam && sample && ray_org && ray_dir);
  ASSERT(sample[0] >= 0 || sample[0] < 1);
  ASSERT(sample[1] >= 0 || sample[1] < 1);
  d3_muld(x, cam->axis_x, sample[0]*2-1);
  d3_muld(y, cam->axis_y, sample[1]*2-1);
  d3_add(ray_dir, d3_add(ray_dir, x, y), cam->axis_z);
  len = d3_normalize(ray_dir, ray_dir);
  ASSERT(len >= 1.e-6);
  d3_set(ray_org, cam->position);
}

