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

#ifndef HTRDR_ATMOSPHERE_SUN_H
#define HTRDR_ATMOSPHERE_SUN_H

#include <rsys/rsys.h>

/* Forward declaration */
struct htrdr;
struct htrdr_atmosphere_sun;
struct ssp_rng;

extern LOCAL_SYM res_T
htrdr_atmosphere_sun_create
  (struct htrdr* htrdr,
   struct htrdr_atmosphere_sun** out_sun);

extern LOCAL_SYM void
htrdr_atmosphere_sun_ref_get
  (struct htrdr_atmosphere_sun* sun);

extern LOCAL_SYM void
htrdr_atmosphere_sun_ref_put
  (struct htrdr_atmosphere_sun* sun);

/* Setup the direction *toward* the sun "center" */
extern LOCAL_SYM void
htrdr_atmosphere_sun_set_direction
  (struct htrdr_atmosphere_sun* sun,
   const double direction[3]); /* Must be normalized */

/* Return a pdf of the sampled dir */
extern LOCAL_SYM double
htrdr_atmosphere_sun_sample_direction
  (struct htrdr_atmosphere_sun* sun,
   struct ssp_rng* rng,
   double dir[3]);

extern LOCAL_SYM double
htrdr_atmosphere_sun_get_solid_angle
  (const struct htrdr_atmosphere_sun* sun);

extern LOCAL_SYM double /* W/m^2/sr/m */
htrdr_atmosphere_sun_get_radiance
  (const struct htrdr_atmosphere_sun* sun,
   const double wavelength);

extern LOCAL_SYM int
htrdr_atmosphere_sun_is_dir_in_solar_cone
  (const struct htrdr_atmosphere_sun* sun,
   const double dir[3]);

#endif /* HTRDR_ATMOSPHERE_SUN_H */
