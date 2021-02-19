/* Copyright (C) 2018, 2019, 2020, 2021 |Meso|Star> (contact@meso-star.com)
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

#ifndef HTRDR_CIE_XYZ_H
#define HTRDR_CIE_XYZ_H

#include "core/htrdr.h"
#include <rsys/rsys.h>

/* Wavelength boundaries of the CIE XYZ color space in nanometers */
#define HTRDR_CIE_XYZ_RANGE_DEFAULT__ {380, 780}
static const double HTRDR_CIE_XYZ_RANGE_DEFAULT[2] =
  HTRDR_CIE_XYZ_RANGE_DEFAULT__;

/* Forward declarations */
struct htrdr;
struct htrdr_cie_xyz;

BEGIN_DECLS

HTRDR_API res_T
htrdr_cie_xyz_create
  (struct htrdr* htrdr,
   const double range[2], /* Must be included in  [380, 780] nanometers */
   const size_t nbands, /* # bands used to discretisze the CIE tristimulus s*/
   struct htrdr_cie_xyz** cie);

HTRDR_API void
htrdr_cie_xyz_ref_get
  (struct htrdr_cie_xyz* cie);

HTRDR_API void
htrdr_cie_xyz_ref_put
  (struct htrdr_cie_xyz* cie);

/* Return a wavelength in nanometer */
HTRDR_API double
htrdr_cie_xyz_sample_X
  (struct htrdr_cie_xyz* cie,
   const double r0, const double r1, /* Canonical numbers in [0, 1[ */
   double* pdf); /* In nm^-1. May be NULL */

/* Return a wavelength in nanometer */
HTRDR_API double
htrdr_cie_xyz_sample_Y
  (struct htrdr_cie_xyz* cie,
   const double r0, const double r1, /* Canonical number in [0, 1[ */
   double* pdf); /* In nm^-1. May be NULL */

/* Return a wavelength in nanometer */
HTRDR_API double
htrdr_cie_xyz_sample_Z
  (struct htrdr_cie_xyz* cie,
   const double r0, const double r1, /* Canonical number in [0, 1[ */
   double* pdf); /* In nm^-1. May be NULL */

END_DECLS

#endif /* HTRDR_cie_xyz_H */

