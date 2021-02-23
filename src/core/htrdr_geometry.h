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

#ifndef HTRDR_GEOMETRY_H
#define HTRDR_GEOMETRY_H

#include "core/htrdr.h"

/* Forware declarations */
struct htrdr;
struct htrdr_geometry;
struct htrdr_interface;
struct htrdr_materials;
struct s3d_hit;
struct ssf_bsdf;

BEGIN_DECLS

HTRDR_CORE_API res_T
htrdr_geometry_create
  (struct htrdr* htrdr,
   const char* obj_filename,
   struct htrdr_materials* mats, /* Library of materials */
   struct htrdr_geometry** geometry);

HTRDR_CORE_API void
htrdr_geometry_ref_get
  (struct htrdr_geometry* geom);

HTRDR_CORE_API void
htrdr_geometry_ref_put
  (struct htrdr_geometry* geom);

HTRDR_CORE_API void
htrdr_geometry_get_interface
  (struct htrdr_geometry* geom,
   const struct s3d_hit* hit,
   struct htrdr_interface* interface);

HTRDR_CORE_API res_T
htrdr_geometry_create_bsdf
  (struct htrdr_geometry* geom,
   const size_t ithread,
   const double wavelength,
   const double pos[3],
   const double dir[3], /* Incoming ray */
   const struct s3d_hit* hit,
   struct htrdr_interface* interf, /* NULL <=> do not return the interface */
   struct ssf_bsdf** bsdf);

HTRDR_CORE_API res_T
htrdr_geometry_trace_ray
  (struct htrdr_geometry* geom,
   const double ray_origin[3],
   const double ray_direction[3], /* Must be normalized */
   const double ray_range[2],
   const struct s3d_hit* prev_hit,/* Previous hit. Avoid self hit. May be NULL*/
   struct s3d_hit* hit);

HTRDR_CORE_API res_T
htrdr_geometry_find_closest_point
  (struct htrdr_geometry* geom,
   const double position[3],
   const double radius,
   struct s3d_hit* hit);

HTRDR_CORE_API void
htrdr_geometry_get_aabb
  (const struct htrdr_geometry* geom,
   double lower[3],
   double upper[3]);

END_DECLS

#endif /* HTRDR_GEOMETRY_H */

