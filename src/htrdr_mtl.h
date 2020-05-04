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

#ifndef HTRDR_MTL_H
#define HTRDR_MTL_H

struct htrdr_mtl;
struct mrumtl;

extern LOCAL_SYM res_T
htrdr_mtl_create
  (struct htrdr* htrdr,
   const char* filename,
   struct htrdr_mtl** mtl);

extern LOCAL_SYM void
htrdr_mtl_ref_get
  (struct htrdr_mtl* mtl);

extern LOCAL_SYM void
htrdr_mtl_ref_put
  (struct htrdr_mtl* mtl);

/* Return NULL if the material name does not exist */
extern const struct mrumtl*
htrdr_mtl_get
  (struct htrdr_mtl* mtl,
   const char* mtl_name);

#endif /* HTRDR_MTL_H */

