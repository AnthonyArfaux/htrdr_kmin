/* Copyright (C) 2018-2019, 2022-2025 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2025 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2025 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2025 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2025 Observatoire de Paris
 * Copyright (C) 2022-2025 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2025 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2025 Université Paul Sabatier
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

#ifndef HTRDR_PLANETS_SOURCE_H
#define HTRDR_PLANETS_SOURCE_H

#include <rsys/rsys.h>

/* Forward declarations */
struct htrdr;
struct htrdr_planets_source;
struct htrdr_planets_source_args;
struct ssp_rng;

struct htrdr_planets_source_spectrum {
  const struct htrdr_planets_source* source;
  double range[2]; /* In nm. Limits are inclusive */
  size_t size; /* Number of elements representing the spectrum */
  const void* buffer; /* Pointer toward the spectrum data */
};
#define HTRDR_PLANETS_SOURCE_SPECTRUM_NULL__ {NULL, {0,0}, 0, NULL}
static const struct htrdr_planets_source_spectrum
HTRDR_PLANETS_SOURCE_SPECTRUM_NULL = HTRDR_PLANETS_SOURCE_SPECTRUM_NULL__;

extern LOCAL_SYM res_T
htrdr_planets_source_create
  (struct htrdr* htrdr,
   const struct htrdr_planets_source_args* args,
   struct htrdr_planets_source** source);

extern LOCAL_SYM void
htrdr_planets_source_ref_get
  (struct htrdr_planets_source* source);

extern LOCAL_SYM void
htrdr_planets_source_ref_put
  (struct htrdr_planets_source* source);

/* Return the pdf of the sampled direction */
extern LOCAL_SYM double
htrdr_planets_source_sample_direction
  (const struct htrdr_planets_source* source,
   struct ssp_rng* rng,
   const double pos[3], /* Position from which direction is sampled */
   double dir[3]);

extern LOCAL_SYM double /* In W/m²/sr/m */
htrdr_planets_source_get_radiance
  (const struct htrdr_planets_source* source,
   const double wlen); /* In nanometers */

/* Return the distance between the source surface and the input position. Can
 * be negative if the position is in the source */
extern LOCAL_SYM double /* In m */
htrdr_planets_source_distance_to
  (const struct htrdr_planets_source* source,
   const double pos[3]);

/* Return 1 if the source is targeted by the submitted ray and 0 otherwise */
extern LOCAL_SYM int
htrdr_planets_source_is_targeted
  (const struct htrdr_planets_source* source,
   const double pos[3], /* Ray origin */
   const double dir[3]);/* Ray direction */

extern LOCAL_SYM res_T
htrdr_planets_source_get_spectral_range
  (const struct htrdr_planets_source* source,
   double range[2]); /* In nm. Limits are inclusive */

extern LOCAL_SYM int
htrdr_planets_source_does_radiance_vary_spectrally
  (const struct htrdr_planets_source* source);

/* Get discrete spectrum data for a given range. If the boundaries of the
 * spectral range do not coincide with a discrete element, their radiance is
 * recovered from the htrdr_planets_source_get_radiance function. Note that
 * this function returns an error if the radiance from the source does not vary
 * spectrally, that is, its radiance is recovered from a constant temperature */
extern LOCAL_SYM res_T
htrdr_planets_source_get_spectrum
  (const struct htrdr_planets_source* source,
   const double range[2], /* In nm. Limits are inclusive */
   struct htrdr_planets_source_spectrum* spectrum);

/* Note that the following function profile corresponds to the type expected by
 * the discrete wavelength distribution
 * (see htrdr_ran_wlen_discrete_create_args structure) */
extern LOCAL_SYM void
htrdr_planets_source_spectrum_at
  (void* spectrum,
   size_t i, /* between [0, spectrum->size[ */
   double* wavelength, /* In nm */
   double* radiance); /* In W/m²/sr/m */

#endif /* HTRDR_PLANETS_SOURCE_H */
