/* Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
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

#ifndef HTRDR_SUN_H
#define HTRDR_SUN_H

#include <rsys/rsys.h>

/* Forward declaration */
struct htrdr;
struct htrdr_sun;
struct ssp_rng;

extern LOCAL_SYM res_T
htrdr_sun_create
  (struct htrdr* htrdr,
   struct htrdr_sun** out_sun);

extern LOCAL_SYM void
htrdr_sun_ref_get
  (struct htrdr_sun* sun);

extern LOCAL_SYM void
htrdr_sun_ref_put
  (struct htrdr_sun* sun);

/* Setup the direction *toward* the sun "center" */
extern LOCAL_SYM void
htrdr_sun_set_direction
  (struct htrdr_sun* sun,
   const double direction[3]); /* Must be normalized */

/* Return a direction that points *toward* the sun */
extern LOCAL_SYM double*
htrdr_sun_sample_direction
  (struct htrdr_sun* sun,
   struct ssp_rng* rng,
   double dir[3]);

extern LOCAL_SYM double
htrdr_sun_get_solid_angle
  (const struct htrdr_sun* sun);

extern LOCAL_SYM double /* W.sr^-1.m^-2 */
htrdr_sun_get_radiance
  (const struct htrdr_sun* sun,
   const double wavelength);

extern LOCAL_SYM int
htrdr_sun_is_dir_in_solar_cone
  (const struct htrdr_sun* sun,
   const double dir[3]);

extern LOCAL_SYM size_t
htrdr_sun_get_spectral_bands_count
  (const struct htrdr_sun* sun);

extern LOCAL_SYM void
htrdr_sun_get_spectral_band_bounds
  (const struct htrdr_sun* sun,
   const size_t ispectral_band,
   double bounds[2]); /* Lower and upper wavelength in nanometer */

/* Return the ranges of the spectral bands where the CIE XYZ color space is
 * defined. CIE XYZ in [band_range[0], band_range[1]] */
extern LOCAL_SYM void
htrdr_sun_get_CIE_XYZ_spectral_bands_range
  (const struct htrdr_sun* sun,
   size_t band_range[2]);

#endif /* HTRDR_SUN_H */
