/* Copyright (C) 2018-2019, 2022-2025 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2025 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2025 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2025 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2025 Observatoire de Paris
 * Copyright (C) 2022-2025 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2025 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2025 Université Paul Sabatier
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

#include "core/htrdr.h"
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

BEGIN_DECLS

/* Trace a ray into a slab composed of a cell infinitely repeated in X and Y */
HTRDR_API res_T
htrdr_slab_trace_ray
  (struct htrdr* htrdr,
   const double org[3],
   const double dir[3],
   const double range[2],
   const double cell_low[3],
   const double cell_upp[3],
   htrdr_trace_cell_T trace_cell,
   const size_t max_steps, /* Max traversed cell */
   void* trace_cell_context);

END_DECLS

#endif /* HTRDR_SLAB_H */

