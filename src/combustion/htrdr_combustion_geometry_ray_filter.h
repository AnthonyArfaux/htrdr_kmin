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
  struct htrdr_geometry* geom;

  /* Name of the medium attached to the interface to filter. An intersection is
   * ignorer if the side of the intersected surface points toward the medium
   * whose name is `medium_name`. */
  const char* medium_name;
};

#define GEOMETRY_RAY_FILTER_CONTEXT_NULL__ {NULL, NULL}
static const struct geometry_ray_filter_context
GEOMETRY_RAY_FILTER_CONTEXT_NULL = GEOMETRY_RAY_FILTER_CONTEXT_NULL__;

/* Filter function used to ignore the intersections with surfaces pointing
 * toward a user defined medium */
extern LOCAL_SYM int
geometry_ray_filter_discard_medium_interface
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   const float ray_range[2],
   void* ray_data, /* Must point to a struct geometry_ray_filter_context */
   void* filter_data);

#endif /* HTRDR_COMBUSTION_GEOMETRY_RAY_FILTER_H */
