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

struct htrdr_combustion_args {
  const char* filename_geom; /* Obj of the combustion chamber */
  const char* filename_tetra; /* Volumetric mesh of the medium */
  const char* filename_therm_props; /* Termodynamic properties of the medium */
  const char* filename_refract_ids; /* Refractive indices in the medium */

  const char* filename_cache; /* Path of the file to store/restore cached data */
  const char* filename_output; /* Name of the output file */

  struct htrdr_args_camera camera; /* Pinhole Camera */

  struct htrdr_args_rectangle laser; /* Laser surface emission */
  double wavelength; /* Wavelength of the laser in nanometer */

  struct htrdr_args_image image; /* Output Image */

  /* RDG-FA parameters */
  double gyration_radius_prefactor;
  double fractal_dimension;

  unsigned grid_max_definition[3]; /* Fixed grid definition along the 3 axes */
  unsigned auto_grid_definition_hint; /* Hint on the grid definition to eval */
  int auto_grid_definition; /* Switch between auto and fixed grid definition */

  double optical_thickness; /* Threshold used during octree building */

  /* Miscellaneous parameters */
  unsigned nthreads; /* Hint on the number of threads to use */
  int precompute_normals; /* Pre-compute the tetrahedra normals */
  int force_overwriting;
  int dump_volumetric_acceleration_structure;
  int verbose; /* Verbosity level */
  int quit; /* Stop the command */
};

#define HTRDR_COMBUSTION_ARGS_DEFAULT__ {                                      \
  NULL, /* Geom filename */                                                    \
  NULL, /* Tetra filename */                                                   \
  NULL, /* Therm props filename */                                             \
  NULL, /* Refractive ids filename */                                          \
                                                                               \
  NULL, /* Cache filename */                                                   \
  NULL, /* Output filename */                                                  \
                                                                               \
  HTRDR_ARGS_CAMERA_DEFAULT__, /* Pinhole camera */                            \
                                                                               \
  HTRDR_ARGS_RECTANGLE_DEFAULT__, /* Laser surface emission */                 \
  532, /* Wavelength in nanometer */                                           \
                                                                               \
  HTRDR_ARGS_IMAGE_DEFAULT__, /* Image */                                      \
                                                                               \
  1.30, /* Gyration radius prefactor */                                        \
  1.80, /* Fractal dimension */                                                \
                                                                               \
  {256, 256, 256}, /* Acceleration grid max definition */                      \
  256, /* Hint on grid Definition in 'auto grid definition' mode */            \
  1, /* Enable/disable 'auto grid definition' mode */                          \
                                                                               \
  1, /* Optical thickness */                                                   \
                                                                               \
  UINT_MAX, /* #threads */                                                     \
  0, /* Precompute normals */                                                  \
  0, /* Force overwriting */                                                   \
  0, /* dump volumetric acceleration structure */                              \
  0, /* Verbose flag */                                                        \
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
