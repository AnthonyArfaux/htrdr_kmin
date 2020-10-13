/* Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
 * Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
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

#ifndef HTRDR_INTERFACE_H
#define HTRDR_INTERFACE_H

#include <star/ssf.h>

/* Forward declaration of external data type */
struct mrumtl;
struct s3d_hit;
struct ssf_bsdf;
struct ssp_rng;

struct htrdr_interface {
  const struct mrumtl* mtl_front;
  const struct mrumtl* mtl_back;
  const struct mrumtl* mtl_thin; /* != NULL <=> thin material */
};
static const struct htrdr_interface HTRDR_INTERFACE_NULL;

extern LOCAL_SYM res_T
htrdr_interface_create_bsdf
  (struct htrdr* htrdr,
   const struct htrdr_interface* interf,
   const size_t ithread,
   const double wavelength,
   const double pos[3],
   const double dir[3], /* Normalized incoming direction */
   struct ssp_rng* rng,
   struct s3d_hit* hit,
   struct ssf_bsdf** bsdf);

#endif /* HTRDR_INTERFACE_H */

