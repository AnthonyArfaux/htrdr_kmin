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

#ifndef HTRDR_COMBUSTION_ARGS_H
#define HTRDR_COMBUSTION_ARGS_H

#include "core/htrdr_args.h"
#include <rsys/rsys.h>

#include <limits.h> /* UINT_MAX support */

enum htrdr_combustion_args_output_type {
  HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE,
  HTRDR_COMBUSTION_ARGS_OUTPUT_LASER_SHEET,
  HTRDR_COMBUSTION_ARGS_OUTPUT_OCTREES,
  HTRDR_COMBUSTION_ARGS_OUTPUT_TYPES_COUNT__
};

enum htrdr_combustion_args_grid_definition_type {
  HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_AUTO,
  HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_FIXED,
  HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_TYPES_COUNT__
};

struct htrdr_combustion_args_grid_definition {
  union {
    unsigned hint; /* Hint on the grid definition to eval */
    unsigned fixed[3]; /* Fixed grid definition along the 3 axis */
  } definition;
  enum htrdr_combustion_args_grid_definition_type type;
};
#define HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_DEFAULT__ {                      \
  {256},                                                                       \
  HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_AUTO                                   \
}
static const struct htrdr_combustion_args_grid_definition
HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_DEFAULT =
 HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_DEFAULT__;

struct htrdr_combustion_args {
  struct htrdr_args_geometry geom; /* Combustion chamber geometry */

  const char* path_tetra; /* Volumetric mesh of the medium */
  const char* path_therm_props; /* Termodynamic properties of the medium */
  const char* path_refract_ids; /* Refractive indices in the medium */

  const char* path_cache; /* Path of the file to store/restore cached data */
  const char* path_output; /* Name of the output file */

  struct htrdr_args_camera camera; /* Pinhole Camera */

  struct htrdr_args_rectangle laser; /* Laser surface emission */
  double wavelength; /* Wavelength of the laser in nanometer */
  double laser_flux_density; /* In W/m^2 */

  struct htrdr_args_image image; /* Output Image */

  /* RDG-FA parameters */
  double fractal_prefactor;
  double fractal_dimension;

  struct htrdr_combustion_args_grid_definition grid;

  double optical_thickness; /* Threshold used during octree building */

  /* Miscellaneous parameters */
  unsigned nthreads; /* Hint on the number of threads to use */
  enum htrdr_combustion_args_output_type output_type;
  int precompute_normals; /* Pre-compute the tetrahedra normals */
  int force_overwriting;
  int verbose; /* Verbosity level */
  int use_simd; /* Use the SIMD instruction set if available */
  int quit; /* Stop the command */
};

#define HTRDR_COMBUSTION_ARGS_DEFAULT__ {                                      \
  HTRDR_ARGS_GEOMETRY_NULL__,                                                  \
                                                                               \
  NULL, /* Tetra path */                                                       \
  NULL, /* Therm props path */                                                 \
  NULL, /* Refractive ids path */                                              \
                                                                               \
  NULL, /* Cache path */                                                       \
  NULL, /* Output path */                                                      \
                                                                               \
  HTRDR_ARGS_CAMERA_DEFAULT__, /* Pinhole camera */                            \
                                                                               \
  HTRDR_ARGS_RECTANGLE_DEFAULT__, /* Laser surface emission */                 \
  532, /* Wavelength in nanometer */                                           \
  1, /* Flux density */                                                        \
                                                                               \
  HTRDR_ARGS_IMAGE_DEFAULT__, /* Image */                                      \
                                                                               \
  1.30, /* Gyration radius prefactor */                                        \
  1.80, /* Fractal dimension */                                                \
                                                                               \
  HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_DEFAULT__,                             \
                                                                               \
  1, /* Optical thickness */                                                   \
                                                                               \
  UINT_MAX, /* #threads */                                                     \
  HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE, /* Output type */                        \
  0, /* Precompute normals */                                                  \
  0, /* Force overwriting */                                                   \
  0, /* Verbose flag */                                                        \
  0, /* Use SIMD */                                                            \
  0  /* Stop the command */                                                    \
}
static const struct htrdr_combustion_args HTRDR_COMBUSTION_ARGS_DEFAULT =
  HTRDR_COMBUSTION_ARGS_DEFAULT__;

extern LOCAL_SYM res_T
htrdr_combustion_args_init
  (struct htrdr_combustion_args* args,
   int argc,
   char** argv);

extern LOCAL_SYM void
htrdr_combustion_args_release
  (struct htrdr_combustion_args* args);

#endif /* HTRDR_COMBUSTION_ARGS_H */
