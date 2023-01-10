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

#ifndef HTRDR_PLANETO_H
#define HTRDR_PLANETO_H

#include "core/htrdr.h"
#include <rsys/rsys.h>

struct htrdr;
struct htrdr_planeto;
struct htrdr_planeto_args;

BEGIN_DECLS

HTRDR_API res_T
htrdr_planeto_create
  (struct htrdr* htrdr,
   const struct htrdr_planeto_args* args,
   struct htrdr_planeto** cmd);

HTRDR_API void
htrdr_planeto_ref_get
  (struct htrdr_planeto* cmd);

HTRDR_API void
htrdr_planeto_ref_put
  (struct htrdr_planeto* cmd);

HTRDR_API res_T
htrdr_planeto_run
  (struct htrdr_planeto* cmd);

HTRDR_API int
htrdr_planeto_main
  (int argc,
   char** argv);

END_DECLS

#endif /* HTRDR_PLANETO_H */
