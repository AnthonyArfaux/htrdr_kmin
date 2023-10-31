/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2023 Observatoire de Paris
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

#ifndef HTRDR_GEOMETRY_H
#define HTRDR_GEOMETRY_H

#include "core/htrdr.h"

#include <star/s3d.h>

/* Forware declarations */
struct htrdr;
struct htrdr_geometry;
struct htrdr_interface;
struct htrdr_materials;
struct s3d_hit;
struct ssf_bsdf;

struct htrdr_geometry_trace_ray_args {
  double ray_org[3];
  double ray_dir[3];
  double ray_range[2];
  struct s3d_hit hit_from; /* Hit from which the ray starts */
  s3d_hit_filter_function_T filter; /* NULL <=> no user defined filter */
  void* filter_context;
};

#define HTRDR_GEOMETRY_TRACE_RAY_ARGS_NULL__ {                                 \
  {0,0,0}, /* Ray origin */                                                    \
  {0,0,1}, /* Ray direction */                                                 \
  {0,DBL_MAX}, /* Ray range */                                                 \
  S3D_HIT_NULL__, /* Hit from */                                               \
  NULL, /* User defined filter function */                                     \
  NULL /* User filter function */                                              \
}
static const struct htrdr_geometry_trace_ray_args
HTRDR_GEOMETRY_TRACE_RAY_ARGS_NULL = HTRDR_GEOMETRY_TRACE_RAY_ARGS_NULL__;

BEGIN_DECLS

HTRDR_API res_T
htrdr_geometry_create
  (struct htrdr* htrdr,
   const char* obj_filename,
   struct htrdr_materials* mats, /* Library of materials */
   struct htrdr_geometry** geometry);

HTRDR_API void
htrdr_geometry_ref_get
  (struct htrdr_geometry* geom);

HTRDR_API void
htrdr_geometry_ref_put
  (struct htrdr_geometry* geom);

HTRDR_API void
htrdr_geometry_get_interface
  (struct htrdr_geometry* geom,
   const struct s3d_hit* hit,
   struct htrdr_interface* interface);

HTRDR_API void
htrdr_geometry_get_hit_position
  (const struct htrdr_geometry* geom,
   const struct s3d_hit* hit,
   double position[3]);

HTRDR_API res_T
htrdr_geometry_trace_ray
  (struct htrdr_geometry* geom,
   const struct htrdr_geometry_trace_ray_args* args,
   struct s3d_hit* hit);

HTRDR_API res_T
htrdr_geometry_find_closest_point
  (struct htrdr_geometry* geom,
   const double position[3],
   const double radius,
   struct s3d_hit* hit);

HTRDR_API void
htrdr_geometry_get_aabb
  (const struct htrdr_geometry* geom,
   double lower[3],
   double upper[3]);

/* Empirical value relative to the extent of the geometry that represents the
 * threshold below which a numerical problem could occur. */
HTRDR_API double
htrdr_geometry_get_epsilon
  (const struct htrdr_geometry* geom);

END_DECLS

#endif /* HTRDR_GEOMETRY_H */

