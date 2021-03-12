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

#ifndef HTRDR_COMBUSTION_LASER_H
#define HTRDR_COMBUSTION_LASER_H

#include "core/htrdr_args.h"

#include <rsys/rsys.h>

struct htrdr_combustion_laser_create_args {
  struct htrdr_args_rectangle surface; /* Surface emission */
};
#define HTRDR_COMBUSTION_LASER_CREATE_ARGS_DEFAULT__ {                         \
  HTRDR_ARGS_RECTANGLE_DEFAULT__                                               \
}
static const struct htrdr_combustion_laser_create_args
HTRDR_COMBUSTION_LASER_CREATE_ARGS_DEFAULT =
  HTRDR_COMBUSTION_LASER_CREATE_ARGS_DEFAULT__;

/* Forward declaration */
struct htrdr;
struct htrdr_combustion_laser;

extern LOCAL_SYM res_T
htrdr_combustion_laser_create
  (struct htrdr* htrdr,
   const struct htrdr_combustion_laser_create_args* args,
   struct htrdr_combustion_laser** laser);

extern LOCAL_SYM void
htrdr_combustion_laser_ref_get
  (struct htrdr_combustion_laser* laser);

extern LOCAL_SYM void
htrdr_combustion_laser_ref_put
  (struct htrdr_combustion_laser* laser);

extern LOCAL_SYM void
htrdr_combustion_laser_trace_ray
  (struct htrdr_combustion_laser* laser,
   const double pos[3],
   const double dir[3],
   const double range[2],
   double distance[2]);

#endif /* HTRDR_COMBUSTION_LASER_H */
