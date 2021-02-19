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

#ifndef HTRDR_MATERIALS_H
#define HTRDR_MATERIALS_H

#include "core/htrdr.h"
#include <rsys/rsys.h>

/* Forward declarations */
struct htrdr_materials;
struct mrumtl;
struct s3d_hit;
struct ssf_bsdf;
struct ssp_rng;

struct htrdr_mtl {
  const char* name;
  const struct mrumtl* mrumtl;
  double temperature;
};
#define HTRDR_MTL_NULL__ {NULL, NULL, 0}
static const struct htrdr_mtl HTRDR_MTL_NULL = HTRDR_MTL_NULL__;

BEGIN_DECLS

HTRDR_CORE_API res_T
htrdr_materials_create
  (struct htrdr* htrdr,
   const char* filename,
   struct htrdr_materials** mats);

HTRDR_CORE_API void
htrdr_materials_ref_get
  (struct htrdr_materials* mats);

HTRDR_CORE_API void
htrdr_materials_ref_put
  (struct htrdr_materials* mats);

/* Return 1 if the material exist and 0 otherwise */
HTRDR_CORE_API int
htrdr_materials_find_mtl
  (struct htrdr_materials* mats,
   const char* mtl_name,
   struct htrdr_mtl* mtl);

HTRDR_CORE_API res_T
htrdr_mtl_create_bsdf
  (struct htrdr* htrdr,
   const struct htrdr_mtl* mtl,
   const size_t ithread,
   const double wavelength,
   struct ssp_rng* rng,
   struct ssf_bsdf** bsdf);

END_DECLS

#endif /* HTRDR_MATERIALS_H */

