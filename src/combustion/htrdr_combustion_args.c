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

#include "combustion/htrdr_combustion_args.h"

#include <getopt.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
print_help(const char* cmd)
{
  ASSERT(cmd);
  printf(
"Usage: %s [<options>] -l <laser> -r REFRACT_IDS -m TETRAHEDRA -p THERMOPROPS\n",
    cmd);
  printf(
"Render a monochromatic image within a sooting flame described according\n"
"to the RDG-FA theory and lightened by a laser source.\n\n");

  printf(
"  -C <camera>    define the rendering point of view. Refer to the\n"
"                 %s man page for the list of camera options.\n", cmd);
  printf(
"  -d             dump volumetric acceleration structures to OUTPUT\n"
"                 and exit.\n");
  printf(
"  -F FRACTAL_DIM value of the fractal dimension to use in the RDG-FA\n"
"                 model. Its default value is %g.\n",
    HTRDR_COMBUSTION_ARGS_DEFAULT.fractal_dimension);
  printf(
"  -f             overwrite the OUTPUT file if it already exists.\n");
  printf(
"  -G PREFACTOR   value gyration radius prefactor to use in the RDG-FA\n"
"                 model. Its default value is %g.\n",
    HTRDR_COMBUSTION_ARGS_DEFAULT.gyration_radius_prefactor);
  printf(
"  -g GEOMETRY    filename of the combustion chamber.\n");
  printf(
"  -h             display this help and exit.\n");
  printf(
"  -i <image>     define the image to compute. Refer to the %s man\n"
"                 page for the list of image options\n", cmd);
  printf(
"  -l <laser>     define the geometry of the laser sheet. Refer to the\n"
"                 %s man page for the list of laser options.\n", cmd);
  printf(
"  -m TETRAHEDRA  path toward the volumetric mesh.\n");
  printf(
"  -N             precompute the tetrahedra normals.\n");
  printf(
"  -O CACHE       filename of the cache file used to store/restore the\n"
"                 volumetric data. By default do not use any cache.\n");
  printf(
"  -o OUTPUT      file where data are written. If not defined, data are\n"
"                 written to standard output.\n");
  printf(
"  -p THERMOPROPS path toward the thermodynamic properties.\n");
  printf(
"  -r REFRACT_ID  path toward the per wavelength refractive\n"
"                 indices.\n");
  printf(
"  -T THRESHOLD   optical thickness used as threshold during the octree\n"
"                 building. By default its value is %g.\n",
    HTRDR_COMBUSTION_ARGS_DEFAULT.optical_thickness);
  printf(
"  -t NTHREADS    hint on the number of threads to use. By default use\n"
"                 as many threads as CPU cores.\n");
  printf(
"  -V <DEFINITION|X,Y,Z>\n"
"                 maximum definition of the volumetric acceleration\n"
"                 grids along the 3 axis. Its default value is\n"
"                 [%u, %u, %u].\n", 
    SPLIT3(HTRDR_COMBUSTION_ARGS_DEFAULT.grid_max_definition));
  printf(
"  -v             make the command verbose.\n");
  printf(
"  -w WAVELENGTH  wavelength definition of the laser in nanometer.\n"
"                 By default its value is %g.\n",
    HTRDR_COMBUSTION_ARGS_DEFAULT.wavelength);

  printf("\n");

  htrdr_fprint_copyright(cmd, stdout);
  htrdr_fprint_license(cmd, stdout);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_combustion_args_init
  (struct htrdr_combustion_args* args,
   int argc,
   char** argv)
{
  int opt;
  res_T res = RES_OK;
  ASSERT(args && argc && argv);

  *args = HTRDR_COMBUSTION_ARGS_DEFAULT;

  while((opt = getopt(argc, argv, "C:dF:fG:g:hi:l:m:NO:o:p:r:T:t:V:vw:")) != -1) {
    switch(opt) {
      /* TODO parse the options */
      case 'h':
        print_help(argv[0]);
        htrdr_combustion_args_release(args);
        args->quit = 1;
        goto exit;
      default: res = RES_BAD_ARG; break;
    }
    if(res != RES_OK) {
      if(optarg) {
        fprintf(stderr, "%s: invalid option argument '%s' -- '%c'\n",
          argv[0], optarg, opt);
      }
      goto error;
    }
  }

exit:
  return res;
error:
  htrdr_combustion_args_release(args);
  goto exit;
}

void
htrdr_combustion_args_release(struct htrdr_combustion_args* args)
{
  ASSERT(args);
  *args = HTRDR_COMBUSTION_ARGS_DEFAULT;
}

