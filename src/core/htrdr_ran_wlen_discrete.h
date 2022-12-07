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

#ifndef HTRDR_RAN_WLEN_DISCRETE_H
#define HTRDR_RAN_WLEN_DISCRETE_H

#include "core/htrdr.h"
#include <rsys/rsys.h>

struct htrdr_ran_wlen_discrete_create_args {
  void (*get)
    (const size_t i,
     double* wlen, /* In nanometer */
     double* radiance, /* In W/m²/sr/m */
     void* ctx);
  size_t nwavelengths;
  void* context; /* User defined data */
};

/* Forware declarations */
struct htrdr;
struct htrdr_ran_wlen_discrete;

BEGIN_DECLS

HTRDR_CORE_API res_T
htrdr_ran_wlen_discrete_create
  (struct htrdr* htrdr,
   const struct htrdr_ran_wlen_discrete_create_args* args,
   struct htrdr_ran_wlen_discrete** ran);

HTRDR_CORE_API void
htrdr_ran_wlen_discrete_ref_get
  (struct htrdr_ran_wlen_discrete* ran);

HTRDR_CORE_API void
htrdr_ran_wlen_discrete_ref_put
  (struct htrdr_ran_wlen_discrete* ran);

HTRDR_CORE_API double /* wavelength in nanometer */
htrdr_ran_wlen_discrete_sample
  (struct htrdr_ran_wlen_discrete* ran,
   const double r0, const double r1, /* Canonical number in [0, 1[ */
   double* pdf); /* In nm⁻¹ */

END_DECLS

#endif /* HTRDR_RAN_WLEN_DISCRETE_H */
