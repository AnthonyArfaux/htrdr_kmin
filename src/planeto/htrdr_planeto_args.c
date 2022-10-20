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

#define _POSIX_C_SOURCE 200112L /* strtok_r support */

#include "planeto/htrdr_planeto_args.h"

#include <getopt.h>
#include <string.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
print_help(const char* cmd)
{
  ASSERT(cmd);
  printf(
"Usage: %s [-dfhv] [-s spectral_domain] [-t threads]\n"
"                     [-T optical_thickness] [-V octree_definition]\n"
"                     [-O octrees_storage] [-o output]\n"
"                     [-a aerosol]... -g gas -G ground\n", cmd);
  printf(
"Simulate radiative transfer in heterogeneous 3D planetary atmosphere.\n"
"See htrdr-planeto(1) man page for details\n\n");
  printf(
"  -a aerosol     define an aerosol\n");
  printf(
"  -d             write the atmospheric acceleration structures\n");
  printf(
"  -f             force overwrite the output file\n");
  printf(
"  -G ground      define the ground of the planet\n");
  printf(
"  -g gas         define the gas mixture\n");
  printf(
"  -h             display this help and exit\n");
  printf(
"  -O octrees_storage\n"
"                 file where atmospheric acceleration structures are\n"
"                 stored/loaded\n");
  printf(
"  -o output      file where the result is written. If not defined,\n"
"                 the result is written to standard output\n");
  printf(
"  -s spectral_domain\n"
"                 define the spectral domain of integration\n");
  printf(
"  -T optical_thickness\n"
"                 optical thickness criteria for octree building.\n"
"                 Default is %g\n",
    HTRDR_PLANETO_ARGS_DEFAULT.optical_thickness);
  printf(
"  -t threads     hint on the number of threads to use.\n"
"                 Default assumes as mayn threads as CPU cores\n");
  printf(
"  -V octree_definition\n"
"                 advice on the definition of the atmospheric\n"
"                 acceleration structures. Default is %u\n",
    HTRDR_PLANETO_ARGS_DEFAULT.octree_definition_hint);
  printf(
"  -v             make the command verbose\n");
  printf("\n");
  htrdr_fprint_license(cmd, stdout);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_planeto_args_init(struct htrdr_planeto_args* args, int argc, char** argv)
{
  int opt;
  res_T res = RES_OK;
  ASSERT(args && argc && argv);

  *args = HTRDR_PLANETO_ARGS_DEFAULT;

  while((opt = getopt(argc, argv, "a:dfG:g:hO:o:s:T:t:V:v")) != -1) {
    switch(opt) {
      case 'h':
        print_help(argv[0]);
        htrdr_planeto_args_release(args);
        args->quit = 1;
        goto exit;
      default: res = RES_BAD_ARG; goto error;
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
  htrdr_planeto_args_release(args);
  goto exit;
}

void
htrdr_planeto_args_release(struct htrdr_planeto_args* args)
{
  ASSERT(args);
  *args = HTRDR_PLANETO_ARGS_DEFAULT;
}
