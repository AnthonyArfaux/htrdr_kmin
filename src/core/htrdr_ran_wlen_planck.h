/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2023 Observatoire de Paris
 * Copyright (C) 2022-2023 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2023 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2023 Université Paul Sabatier
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

#ifndef HTRDR_RAN_WLEN_PLANCK_H
#define HTRDR_RAN_WLEN_PLANCK_H

#include "core/htrdr.h"
#include <rsys/rsys.h>

#define HTRDR_WLEN_RAN_PLANCK_CONTINUE 0

/* Forward declarations */
struct htrdr;
struct htrdr_ran_wlen_planck;

BEGIN_DECLS

HTRDR_API res_T
htrdr_ran_wlen_planck_create
  (struct htrdr* htrdr,
   const double range[2], 
   /* # bands used to discretisze the spectral domain. HTRDR_WLEN_RAN_CONTINUE
    * <=> no discretisation */
   const size_t nbands, /* Hint on #bands used to discretised th CDF */
   const double ref_temperature, /* Reference temperature */
   struct htrdr_ran_wlen_planck** planck);

HTRDR_API void
htrdr_ran_wlen_planck_ref_get
  (struct htrdr_ran_wlen_planck* planck);

HTRDR_API void
htrdr_ran_wlen_planck_ref_put
  (struct htrdr_ran_wlen_planck* planck);

/* Return a wavelength in nanometer */
HTRDR_API double
htrdr_ran_wlen_planck_sample
  (const struct htrdr_ran_wlen_planck* planck,
   const double r0, /* Canonical number in [0, 1[ */
   const double r1, /* Canonical number in [0, 1[ */
   double* pdf); /* In nm^-1. May be NULL */

END_DECLS

#endif /* HTRDR_RAN_WLEN_PLANCK_H */

