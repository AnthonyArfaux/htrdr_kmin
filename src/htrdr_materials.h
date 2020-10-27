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

#ifndef HTRDR_MATERIALS_H
#define HTRDR_MATERIALS_H

#include <rsys/rsys.h>

struct htrdr_materials;
struct mrumtl;

struct htrdr_mtl {
  const char* name;
  const struct mrumtl* mrumtl;
};
static const struct htrdr_mtl HTRDR_MTL_NULL;

extern LOCAL_SYM res_T
htrdr_materials_create
  (struct htrdr* htrdr,
   const char* filename,
   struct htrdr_materials** mats);

extern LOCAL_SYM void
htrdr_materials_ref_get
  (struct htrdr_materials* mats);

extern LOCAL_SYM void
htrdr_materials_ref_put
  (struct htrdr_materials* mats);

/* Return 1 if the material exist and 0 otherwise */
extern LOCAL_SYM int
htrdr_materials_find_mtl
  (struct htrdr_materials* mats,
   const char* mtl_name,
   struct htrdr_mtl* mtl);

#endif /* HTRDR_MATERIALS_H */

