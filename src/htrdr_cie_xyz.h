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

#ifndef HTRDR_CIE_XYZ_H
#define HTRDR_CIE_XYZ_H

#include <rsys/rsys.h>

struct htrdr;
struct htrdr_cie_xyz;

/* Wavelength boundaries of the CIE XYZ color space in nanometers */
static const double HTRDR_CIE_XYZ_RANGE_DEFAULT[2] = {380, 780};

extern LOCAL_SYM res_T
htrdr_cie_xyz_create
  (struct htrdr* htrdr,
   const double range[2], /* Must be included in  [380, 780] nanometers */
   const size_t nbands, /* # bands used to discretisze the CIE tristimulus s*/
   struct htrdr_cie_xyz** cie);

extern LOCAL_SYM void
htrdr_cie_xyz_ref_get
  (struct htrdr_cie_xyz* cie);

extern LOCAL_SYM void
htrdr_cie_xyz_ref_put
  (struct htrdr_cie_xyz* cie);

/* Return a wavelength in nanometer */
extern LOCAL_SYM double
htrdr_cie_xyz_sample_X
  (struct htrdr_cie_xyz* cie,
   const double r); /* Canonical number in [0, 1[ */

/* Return a wavelength in nanometer */
extern LOCAL_SYM double
htrdr_cie_xyz_sample_Y
  (struct htrdr_cie_xyz* cie,
   const double r); /* Canonical number in [0, 1[ */

/* Return a wavelength in nanometer */
extern LOCAL_SYM double
htrdr_cie_xyz_sample_Z
  (struct htrdr_cie_xyz* cie,
   const double r); /* Canonical number in [0, 1[ */

#endif /* HTRDR_cie_xyz_H */

