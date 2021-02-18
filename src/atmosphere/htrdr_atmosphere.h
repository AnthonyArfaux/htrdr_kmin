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

#ifndef HTRDR_ATMOSPHERE_H
#define HTRDR_ATMOSPHERE_H

#include <rsys/rsys.h>

/* Forward declarations */
struct htrdr;
struct htrdr_atmosphere;
struct htrdr_atmosphere_args;

extern LOCAL_SYM res_T
htrdr_atmosphere_create
  (struct htrdr* htrdr,
   const struct htrdr_atmosphere_args* args,
   struct htrdr_atmosphere** cmd);

extern LOCAL_SYM void
htrdr_atmosphere_ref_get
  (struct htrdr_atmosphere* cmd);

extern LOCAL_SYM void
htrdr_atmosphere_ref_put
  (struct htrdr_atmosphere* cmd);

extern LOCAL_SYM res_T
htrdr_atmosphere_run
  (struct htrdr_atmosphere* cmd);

#endif /* HTRDR_ATMOSPHERE_H */

