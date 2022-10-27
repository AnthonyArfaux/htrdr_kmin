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

#ifndef HTRDR_PLANETO_SOURCE_H
#define HTRDR_PLANETO_SOURCE_H

#include <rsys/rsys.h>

/* Forward declarations */
struct htrdr_planeto;
struct htrdr_planeto_source;
struct htrdr_planeto_source_args;
struct ssp_rng;

extern LOCAL_SYM res_T
htrdr_planeto_source_create
  (struct htrdr_planeto* cmd,
   const struct htrdr_planeto_source_args* args,
   struct htrdr_planeto_source** source);

extern LOCAL_SYM void
htrdr_planeto_source_ref_get
  (struct htrdr_planeto_source* source);

extern LOCAL_SYM void
htrdr_planeto_source_ref_put
  (struct htrdr_planeto_source* source);

/* Return the pdf of the sampled direction */
extern LOCAL_SYM double
htrdr_planeto_source_sample_direction
  (const struct htrdr_planeto_source* source,
   struct ssp_rng* rng,
   const double pos[3], /* Position from which direction is sampled */
   double dir[3]);

extern LOCAL_SYM double /* In W/m²/sr/m */
htrdr_planeto_source_get_radiance
  (const struct htrdr_planeto_source* source,
   const double wlen); /* In nanometers */

/* Return 1 if the source is targeted by the submitted ray and 0 otherwise */
extern LOCAL_SYM int
htrdr_planeto_source_is_targeted
  (const struct htrdr_planeto_source* source,
   const double pos[3], /* Ray origin */
   const double dir[3]);/* Ray direction */

#endif /* HTRDR_PLANETO_SOURCE_H */
