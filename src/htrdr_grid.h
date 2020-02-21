/* Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
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

#ifndef HTRDR_GRID_H
#define HTRDR_GRID_H

#include <rsys/rsys.h>

/* Forwared declarations */
struct htrdr;
struct htrdr_grid;

/*******************************************************************************
 * Out of core regular grid
 ******************************************************************************/
extern LOCAL_SYM res_T
htrdr_grid_create
  (struct htrdr* htrdr,
   const size_t definition[3], /* #voxels in X, Y and Z */
   const size_t sizeof_cell, /* Size of an cell in Bytes */
   const char* filename, /* Filename where the grid data are stored */
   const int force_overwrite, /* Force the overwrite of the grid data */
   struct htrdr_grid** grid);

extern LOCAL_SYM res_T
htrdr_grid_open
  (struct htrdr* htrdr,
   const char* filename,
   struct htrdr_grid** grid);

extern LOCAL_SYM void
htrdr_grid_ref_get
  (struct htrdr_grid* grid);

extern LOCAL_SYM void
htrdr_grid_ref_put
  (struct htrdr_grid* grid);

/* Fetch the grid data from its 3D index */
extern LOCAL_SYM void*
htrdr_grid_at
  (struct htrdr_grid* grid,
   const size_t xyz[3]);

/* Retrieve the voxel data from its morton code. The morton code is computed
 * from the 3D indices following the morton_xyz_encode_u21 convention. */
extern LOCAL_SYM void*
htrdr_grid_at_mcode
  (struct htrdr_grid* grid,
   const uint64_t mcode);

extern LOCAL_SYM void
htrdr_grid_get_definition
  (struct htrdr_grid* grid,
   size_t definition[3]);

#endif /* HTRDR_GRID_H */

