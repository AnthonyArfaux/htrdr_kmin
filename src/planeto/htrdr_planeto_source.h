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
struct htrdr;
struct htrdr_planeto_source;
struct htrdr_planeto_source_args;
struct ssp_rng;

struct htrdr_planeto_source_spectrum {
  const struct htrdr_planeto_source* source;
  double range[2]; /* In nm. Limits are inclusive */
  size_t size; /* Number of elements representing the spectrum */
  const void* buffer; /* Pointer toward the spectrum data */
};

extern LOCAL_SYM res_T
htrdr_planeto_source_create
  (struct htrdr* htrdr,
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

/* Return the distance between the source surface and the input position. Can
 * be negative if the position is in the source */
extern LOCAL_SYM double /* In m */
htrdr_planeto_source_distance_to
  (const struct htrdr_planeto_source* source,
   const double pos[3]);

/* Return 1 if the source is targeted by the submitted ray and 0 otherwise */
extern LOCAL_SYM int
htrdr_planeto_source_is_targeted
  (const struct htrdr_planeto_source* source,
   const double pos[3], /* Ray origin */
   const double dir[3]);/* Ray direction */

extern LOCAL_SYM res_T
htrdr_planeto_source_get_spectral_range
  (const struct htrdr_planeto_source* source,
   double range[2]); /* In nm. Limits are inclusive */

extern LOCAL_SYM int
htrdr_planeto_source_does_radiance_vary_spectrally
  (const struct htrdr_planeto_source* source);

/* Get discrete spectrum data for a given range. If the boundaries of the
 * spectral range do not coincide with a discrete element, their radiance is
 * recovered from the htrdr_planeto_source_get_radiance function. Note that
 * this function returns an error if the radiance from the source does not vary
 * spectrally, that is, its radiance is recovered from a constant temperature */
extern LOCAL_SYM res_T
htrdr_planeto_source_get_spectrum
  (const struct htrdr_planeto_source* source,
   const double range[2], /* In nm. Limits are inclusive */
   struct htrdr_planeto_source_spectrum* spectrum);

/* Note that the following function profile corresponds to the type expected by
 * the discrete wavelength distribution
 * (see htrdr_ran_wlen_discrete_create_args structure) */
extern LOCAL_SYM void
htrdr_planeto_source_spectrum_at
  (void* spectrum,
   size_t i, /* between [0, spectrum->size[ */
   double* wavelength, /* In nm */
   double* radiance); /* In W/m²/sr/m */

#endif /* HTRDR_PLANETO_SOURCE_H */
