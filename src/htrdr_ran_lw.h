/* Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
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

#ifndef HTRDR_RAN_LW_H
#define HTRDR_RAN_LW_H

#include <rsys/rsys.h>

#define HTRDR_RAN_LW_CONTINUE 0

struct htrdr;
struct htrdr_ran_lw;

extern LOCAL_SYM res_T
htrdr_ran_lw_create
  (struct htrdr* htrdr,
   const double range[2], /* Must be included in [1000, 100000] nanometers */
   /* # bands used to discretisze the LW domain. HTRDR_RAN_LW_CONTINUE <=> no
    * discretisation */
   const size_t nbands, /* Hint on #bands used to discretised the CDF */
   const double ref_temperature, /* Reference temperature */
   struct htrdr_ran_lw** ran_lw);

extern LOCAL_SYM void
htrdr_ran_lw_ref_get
  (struct htrdr_ran_lw* ran_lw);

extern LOCAL_SYM void
htrdr_ran_lw_ref_put
  (struct htrdr_ran_lw* ran_lw);

/* Return a wavelength in nanometer */
extern LOCAL_SYM double
htrdr_ran_lw_sample
  (const struct htrdr_ran_lw* ran_lw,
   const double r0, /* Canonical number in [0, 1[ */
   const double r1, /* Canonical number in [0, 1[ */
   double* pdf); /* May be NULL */

#endif /* HTRDR_RAN_LW_H */

