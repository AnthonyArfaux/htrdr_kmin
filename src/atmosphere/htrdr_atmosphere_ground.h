/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
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

#ifndef HTRDR_ATMOSPHERE_GROUND_H
#define HTRDR_ATMOSPHERE_GROUND_H

#include <rsys/rsys.h>

/* Forward declarations */
struct htrdr;
struct htrdr_atmosphere_ground;
struct htrdr_interface;
struct htrdr_materials;
struct s3d_hit;
struct ssf_bsdf;

extern LOCAL_SYM res_T
htrdr_atmosphere_ground_create
  (struct htrdr* htrdr,
   const char* obj_filename, /* May be NULL <=> No ground geometry */
   struct htrdr_materials* mats, /* May be NULL if no ground geometry */
   const int repeat_ground, /* Infinitely repeat the ground in X and Y */
   struct htrdr_atmosphere_ground** ground);

extern LOCAL_SYM void
htrdr_atmosphere_ground_ref_get
  (struct htrdr_atmosphere_ground* ground);

extern LOCAL_SYM void
htrdr_atmosphere_ground_ref_put
  (struct htrdr_atmosphere_ground* ground);

extern LOCAL_SYM void
htrdr_atmosphere_ground_get_interface
  (struct htrdr_atmosphere_ground* ground,
   const struct s3d_hit* hit,
   struct htrdr_interface* interface);

extern LOCAL_SYM res_T
htrdr_atmosphere_ground_create_bsdf
  (struct htrdr_atmosphere_ground* ground,
   const size_t ithread,
   const double wavelength,
   const double pos[3],
   const double dir[3], /* Incoming ray */
   const struct s3d_hit* hit,
   struct htrdr_interface* interf, /* NULL <=> do not return the interface */
   struct ssf_bsdf** bsdf);

extern LOCAL_SYM res_T
htrdr_atmosphere_ground_trace_ray
  (struct htrdr_atmosphere_ground* ground,
   const double ray_origin[3],
   const double ray_direction[3], /* Must be normalized */
   const double ray_range[2],
   const struct s3d_hit* prev_hit,/* Previous hit. Avoid self hit. May be NULL*/
   struct s3d_hit* hit);

extern LOCAL_SYM res_T
htrdr_atmosphere_ground_find_closest_point
  (struct htrdr_atmosphere_ground* ground,
   const double position[3],
   const double radius,
   struct s3d_hit* hit);

#endif /* HTRDR_ATMOSPHERE_GROUND_H */

