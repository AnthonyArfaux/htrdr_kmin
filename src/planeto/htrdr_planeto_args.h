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

#ifndef HTRDR_PLANETO_ARGS_H
#define HTRDR_PLANETO_ARGS_H

#include "core/htrdr_args.h"

#include <rad-net/rnatm.h>
#include <rsys/rsys.h>

#include <limits.h> /* UINT_MAX */

struct htrdr_planeto_ground_args {
  char* smsh_filename; /* The Star-Mesh geometry */
  char* props_filename; /* Per triangle physical properties */
  char* mtllst_filename; /* List of used materials */
};
#define HTRDR_PLANETO_GROUND_ARGS_NULL__ {NULL,NULL,NULL}
static const struct htrdr_planeto_ground_args HTRDR_PLANETO_GROUND_ARGS_NULL =
  HTRDR_PLANETO_GROUND_ARGS_NULL__;

struct htrdr_planeto_args {
  /* System data */
  struct rnatm_gas_args gas;
  struct rnatm_aerosol_args* aerosols;
  size_t naerosols;
  struct htrdr_planeto_ground_args ground;

  /* Read/Write file where octrees are stored. May be NULL => octres are built
   * at runtime and kept in memory */
  char* octrees_storage;

  unsigned octree_definition_hint; /* Hint on octree definition */
  double optical_thickness; /* Threshold used during octree building */

  char* output; /* File where the result is written */

  /* Integration spectral domain */
  struct htrdr_args_spectral spectral_domain;

  /* Miscellaneous arguments */
  unsigned nthreads; /* Hint on the nimber of threads to use */
  int force_output_overwrite; /* Replace output if it exists */
  int verbose; /* Verbose level */
  int quit; /* Stop the command */ 
};
#define HTRDR_PLANETO_ARGS_DEFAULT__ {                                         \
  RNATM_GAS_ARGS_NULL__, /* Gas */                                             \
  NULL, /* List of aerosols */                                                 \
  0, /* Number of aerosols */                                                  \
  HTRDR_PLANETO_GROUND_ARGS_NULL__, /* Ground */                               \
                                                                               \
  NULL, /* File where to dump octrees */                                       \
                                                                               \
  512, /* Hint on octree definition */                                         \
  10, /* Optical thickness criteria */                                         \
                                                                               \
  NULL, /* Ouput file */                                                       \
                                                                               \
  HTRDR_ARGS_SPECTRAL_DEFAULT__, /* Spectral domain */                         \
                                                                               \
  UINT_MAX, /* Number of threads */                                            \
  0, /* Force output overwrite */                                              \
  0, /* Verbosity level */                                                     \
  0 /* Stop the command */                                                     \
}
static const struct htrdr_planeto_args HTRDR_PLANETO_ARGS_DEFAULT =
  HTRDR_PLANETO_ARGS_DEFAULT__;

extern LOCAL_SYM res_T
htrdr_planeto_args_init
  (struct htrdr_planeto_args* args,
   int argc,
   char** argv);

extern LOCAL_SYM void
htrdr_planeto_args_release
  (struct htrdr_planeto_args* args);

#endif /* HTRDR_PLANETO_ARGS_H */
