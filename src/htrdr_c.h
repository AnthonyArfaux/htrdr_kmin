/* Copyright (C) 2018 Université Paul Sabatier, |Meso|Star>
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

#ifndef HTRDR_C_H
#define HTRDR_C_H

#define SW_WAVELENGTH_MIN 380 /* In nanometer */
#define SW_WAVELENGTH_MAX 780 /* In nanometer */

/* In nanometer */
static FINLINE double
wavenumber_to_wavelength(const double nu/*In cm^-1*/)
{
  return 1.e7 / nu;
}

/* In cm^-1 */
static FINLINE double
wavelength_to_wavenumber(const double lambda/*In nanometer*/)
{
  return wavenumber_to_wavelength(lambda);
}

extern LOCAL_SYM res_T
is_file_updated
  (struct htrdr* htrdr,
   const char* filename,
   int* is_upd);

#endif /* HTRDR_C_H */

