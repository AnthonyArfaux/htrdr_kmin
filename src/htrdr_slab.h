/* Copyright (C) 2018 CNRS, Université Paul Sabatier, |Meso|Star>
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

#ifndef HTRDR_SLAB_H
#define HTRDR_SLAB_H

#include <rsys/rsys.h>

/* Forward declaration */
struct htrdr;

typedef res_T
(*htrdr_trace_cell_T)
  (const double org[3], /* Ray origin */
   const double dir[3], /* Ray direction. Must be normalized */
   const double range[2], /* Ray range */
   void* ctx, /* User defined data */
   int* hit); /* Hit something ? */

/* Trace a ray into a slab composed of a cell infinitely repeated in X and Y */
extern LOCAL_SYM res_T
htrdr_slab_trace_ray
  (struct htrdr* htrdr,
   const double org[3],
   const double dir[3],
   const double range[2],
   const double cell_low[2],
   const double cell_upp[2],
   htrdr_trace_cell_T trace_cell,
   const size_t max_steps, /* Max traversed cell */
   void* trace_cell_context);

#endif /* HTRDR_SLAB_H */

