/* Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
 * Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
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

#ifndef HTRDR_CAMERA_H
#define HTRDR_CAMERA_H

#include <rsys/rsys.h>

/* Forward declarations */
struct htrdr;
struct htrdr_camera;

extern LOCAL_SYM res_T
htrdr_camera_create
  (struct htrdr* htrdr,
   const double position[3],
   const double target[3],
   const double up[3],
   const double proj_ratio, /* Width / Height */
   const double fov, /* In radian */
   struct htrdr_camera** cam);

extern LOCAL_SYM void
htrdr_camera_ref_get
  (struct htrdr_camera* cam);

extern LOCAL_SYM void
htrdr_camera_ref_put
  (struct htrdr_camera* cam);

extern LOCAL_SYM void
htrdr_camera_ray
  (const struct htrdr_camera* cam,
   const double sample[2], /* In [0, 1[ */
   double ray_org[3],
   double ray_dir[3]);

#endif /* HTRDR_CAMERA_H */

