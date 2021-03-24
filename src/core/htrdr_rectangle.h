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

#ifndef HTRDR_RECTANGLE_H
#define HTRDR_RECTANGLE_H

#include "core/htrdr.h"
#include <rsys/rsys.h>

/* Forwar declarations */
struct htrdr;
struct htrdr_rectangle; /* 2D rectangle transformed in 3D */

BEGIN_DECLS

HTRDR_CORE_API res_T
htrdr_rectangle_create
  (struct htrdr* htrdr,
   const double sz[2], /* Size of the rectangle along its local X and Y axis */
   const double pos[3], /* World space position of the rectangle center */
   const double tgt[3], /* Targeted point */
   const double up[3], /* vector orthogonal to the rectangle X axis */
   struct htrdr_rectangle** rect);

HTRDR_CORE_API void
htrdr_rectangle_ref_get
  (struct htrdr_rectangle* rect);

HTRDR_CORE_API void
htrdr_rectangle_ref_put
  (struct htrdr_rectangle* rect);

HTRDR_CORE_API void
htrdr_rectangle_sample_pos
  (const struct htrdr_rectangle* rect,
   const double sample[2], /* In [0, 1[ */
   double pos[3]);

HTRDR_CORE_API void
htrdr_rectangle_get_normal
  (const struct htrdr_rectangle* rect,
   double normal[3]);

HTRDR_CORE_API void
htrdr_rectangle_get_center
  (const struct htrdr_rectangle* rect,
   double center[3]);

HTRDR_CORE_API double*
htrdr_rectangle_get_transform
  (const struct htrdr_rectangle* rect,
   double transform[12]);

HTRDR_CORE_API double*
htrdr_rectangle_get_transform_inverse
  (const struct htrdr_rectangle* rect,
   double transform_inverse[12]);

END_DECLS

#endif /* HTRDR_RECTANGLE_H */

