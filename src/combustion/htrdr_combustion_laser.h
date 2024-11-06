/* Copyright (C) 2018-2019, 2022-2024 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2024 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2024 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2024 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2024 Observatoire de Paris
 * Copyright (C) 2022-2024 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2024 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2024 Université Paul Sabatier
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

#ifndef HTRDR_COMBUSTION_LASER_H
#define HTRDR_COMBUSTION_LASER_H

#include "core/htrdr_args.h"

#include <rsys/rsys.h>

/* Monochromatic laser */
struct htrdr_combustion_laser_create_args {
  struct htrdr_args_rectangle surface; /* Surface emission */
  double wavelength; /* In nanometers */
  double flux_density; /* In W/m^2 */
};
#define HTRDR_COMBUSTION_LASER_CREATE_ARGS_DEFAULT__ {                         \
  HTRDR_ARGS_RECTANGLE_DEFAULT__,                                              \
  -1, /* Wavelength */                                                         \
  -1, /* Flux density */                                                       \
}
static const struct htrdr_combustion_laser_create_args
HTRDR_COMBUSTION_LASER_CREATE_ARGS_DEFAULT =
  HTRDR_COMBUSTION_LASER_CREATE_ARGS_DEFAULT__;

struct htrdr_combustion_laser_mesh {
  double vertices[8/*#vertices*/*3/*#coords per vertex*/];
  /* Triangle are clock wise ordered from the inside of the laser sheet */
  unsigned triangles[10/*#triangles*/*3/*#vertices per triangle*/];
  unsigned nvertices;
  unsigned ntriangles;
};

/* Syntactic sugar to check if a laser sheet hit is valid */
#define HTRDR_COMBUSTION_LASER_HIT_NONE(Hit) ((Hit)[0] >= FLT_MAX)

/* Forward declaration */
struct htrdr;
struct htrdr_combustion_laser;

BEGIN_DECLS

extern LOCAL_SYM res_T
htrdr_combustion_laser_create
  (struct htrdr* htrdr,
   const struct htrdr_combustion_laser_create_args* args,
   struct htrdr_combustion_laser** laser);

extern LOCAL_SYM void
htrdr_combustion_laser_ref_get
  (struct htrdr_combustion_laser* laser);

extern LOCAL_SYM void
htrdr_combustion_laser_ref_put
  (struct htrdr_combustion_laser* laser);

extern LOCAL_SYM void
htrdr_combustion_laser_trace_ray
  (struct htrdr_combustion_laser* laser,
   const double pos[3],
   const double dir[3],
   const double range[2],
   double distance[2]);

extern LOCAL_SYM void
htrdr_combustion_laser_get_mesh
  (const struct htrdr_combustion_laser* laser,
   /* Max distance of the laser mesh along its infinite dimension */
   const double extent,
   struct htrdr_combustion_laser_mesh* mesh);

extern LOCAL_SYM void
htrdr_combustion_laser_get_position
  (const struct htrdr_combustion_laser* laser,
   double pos[3]);

extern LOCAL_SYM void
htrdr_combustion_laser_get_direction
  (const struct htrdr_combustion_laser* laser,
   double dir[3]); /* Normalized */

extern LOCAL_SYM double /* In W.m^2 */
htrdr_combustion_laser_get_flux_density
  (const struct htrdr_combustion_laser* laser);

extern LOCAL_SYM double /* In nm */
htrdr_combustion_laser_get_wavelength
  (const struct htrdr_combustion_laser* laser);

extern LOCAL_SYM double
htrdr_combustion_laser_compute_surface_plane_distance
  (const struct htrdr_combustion_laser* laser,
   const double pos[3]);

END_DECLS

#endif /* HTRDR_COMBUSTION_LASER_H */
