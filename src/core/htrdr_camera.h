/* Copyright (C) 2018, 2019, 2020, 2021 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019, 2021 CNRS
 * Copyright (C) 2018, 2019 Université Paul Sabatier
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

#include "core/htrdr.h"
#include <rsys/rsys.h>

/* Forward declarations */
struct htrdr;
struct htrdr_camera;

BEGIN_DECLS

HTRDR_API res_T
htrdr_camera_create
  (struct htrdr* htrdr,
   const double position[3],
   const double target[3],
   const double up[3],
   const double proj_ratio, /* Width / Height */
   const double fov, /* In radian */
   struct htrdr_camera** cam);

HTRDR_API void
htrdr_camera_ref_get
  (struct htrdr_camera* cam);

HTRDR_API void
htrdr_camera_ref_put
  (struct htrdr_camera* cam);

HTRDR_API void
htrdr_camera_ray
  (const struct htrdr_camera* cam,
   const double sample[2], /* In [0, 1[ */
   double ray_org[3],
   double ray_dir[3]);

END_DECLS

#endif /* HTRDR_CAMERA_H */

