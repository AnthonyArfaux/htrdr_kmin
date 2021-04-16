/* Copyright (C) 2018, 2019, 2020, 2021 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019, 2021 CNRS
 * Copyright (C) 2018, 2019, Université Paul Sabatier
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

#ifndef HTRDR_COMBUSTION_GEOMETRY_RAY_FILTER_H
#define HTRDR_COMBUSTION_GEOMETRY_RAY_FILTER_H

#include <star/s3d.h>
#include <rsys/rsys.h>

struct geometry_ray_filter_context {
  /* This is the real hit to use. The one returned by the trace ray function is
   * only used as a "early stop" of the tray tracing process */
  struct s3d_hit hit;

  struct s3d_hit hit_tmp__; /* Temporary internal data */
};

#define GEOMETRY_RAY_FILTER_CONTEXT_NULL__ {S3D_HIT_NULL__, S3D_HIT_NULL__}
static const struct geometry_ray_filter_context
GEOMETRY_RAY_FILTER_CONTEXT_NULL = GEOMETRY_RAY_FILTER_CONTEXT_NULL__;

/* Filter function used to ignore the first intersection along the ray. Note
 * that the intersection to use is output in the 'hit' member variable of the
 * 'struct geometry_filter_context' data structure pointed to by the 'ray_data'
 * input argument. The intersection returned by the ray tracing procedure
 * itself is actually just a temporary variable used to prematurely stop ray
 * tracing. */
extern LOCAL_SYM int
geometry_ray_filter_discard_first_hit
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   void* ray_data, /* Must point to a struct geometry_filter_context */
   void* filter_data);

/* Quite similar filtering function to the previous one, but this time it
 * filters the last intersection along the ray, not the first. This means that
 * if the ray intersects only one surface, it is assumed that the ray does not
 * intersect anything. */
extern LOCAL_SYM int
geometry_ray_filter_discard_last_hit
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   void* ray_data, /* Must point to a struct geometry_filter_context */
   void* filter_data);

#endif /* HTRDR_COMBUSTION_GEOMETRY_RAY_FILTER_H */
