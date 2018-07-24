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

#ifndef HTRDR_ARGS_H
#define HTRDR_ARGS_H

#include <rsys/rsys.h>

struct htrdr_args {
  const char* filename_les; /* Path toward the HTCP file */
  const char* filename_mie; /* Path toward the Mie properties */
  const char* filename_obj; /* Path toward the 3D geometry */
  const char* output;

  struct {
    double pos[3]; /* Position */
    double tgt[3]; /* Target */
    double up[3]; /* Up vector */
    double sz[2]; /* Plane size in world space */
  } rectangle;

  struct {
    double pos[3];
    double tgt[3];
    double up[3];
    double fov_x; /* In degrees */
  } camera;

  struct {
    unsigned definition[2]; /* #pixels in X and Y */
    unsigned spp; /* #samples per pixel */
  } image;

  double main_dir[3];

  unsigned nthreads; /* Hint on the number of threads to use */
  int force_overwriting;
  int dump_vtk; /* Dump the loaded cloud properties in a VTK file */
  int verbose; /* Verbosity level */
  int quit;  /* Qui the application */
};

#define HTRDR_ARGS_DEFAULT__ {                                                 \
  NULL, /* LES filename */                                                     \
  NULL, /* Mie filename */                                                     \
  NULL, /* Obj filename */                                                     \
  NULL, /* Output filename */                                                  \
  {                                                                            \
    {0, 0, 0}, /* plane position */                                            \
    {0, 0, 1}, /* plane target */                                              \
    {0, 1, 0}, /* plane up */                                                  \
    {1, 1},    /* plane size */                                                \
  }, {                                                                         \
    {0, 0, 0}, /* Camera position */                                           \
    {0, 0, 1}, /* Camera target */                                             \
    {0, 1, 0}, /* Camera up */                                                 \
    70, /* Horizontal field of view */                                         \
  }, {                                                                         \
    {32, 32},  /* image definition */                                          \
    1          /* #samples per pixel */                                        \
  },                                                                           \
  {0, 0, -1}, /* Main direction */                                             \
  (unsigned)~0, /* #threads */                                                 \
  0, /* Force overwriting */                                                   \
  0, /* dump VTK */                                                            \
  0, /* Verbose flag */                                                        \
  0  /* Quit the application */                                                \
}
static const struct htrdr_args HTRDR_ARGS_DEFAULT = HTRDR_ARGS_DEFAULT__;

extern LOCAL_SYM res_T
htrdr_args_init
  (struct htrdr_args* args,
   int argc,
   char** argv);

extern LOCAL_SYM void
htrdr_args_release
  (struct htrdr_args* args);

#endif /* HTRDR_ARGS_H */

