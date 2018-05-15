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

#include "htrdr_args.h"
#include <getopt.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
print_help(const char* cmd)
{
  ASSERT(cmd);
  printf("Usage: %s -i INPUT [OPIONS]\n\n", cmd);
  printf(
"  -d               dump octree data to OUTPUT wrt the VTK ASCII file format.\n");
  printf(
"  -f               overwrite the OUTPUT file if it already exists.\n");
  printf(
"  -h               display this help and exit.\n");
  printf(
"  -i INPUT         path of the input HTCP file.\n");
  printf(
"  -o OUTPUT        file where data are written. If not defined, data are\n"
"                   written to standard output.\n");
  printf(
"  -v               make the program more verbose.\n");
  printf("\n");
  printf(
"%s (C) 2018 Université Paul Sabatier, |Meso|Star>. This is free software\n"
"released under the GNU GPL license, version 3 or later. You are free to change\n"
"or redistribute it under certain conditions <http://gnu.org/licenses/gpl.html>.\n",
    cmd);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_args_init(struct htrdr_args* args, int argc, char** argv)
{
  int opt;
  res_T res = RES_OK;
  ASSERT(args && argc && argv);

  *args = HTRDR_ARGS_DEFAULT;

  while((opt = getopt(argc, argv, "dfhi:o:v")) != -1) {
    switch(opt) {
      case 'd': args->dump_vtk = 1; break;
      case 'f': args->force_overwriting = 1; break;
      case 'h':
        print_help(argv[0]);
        htrdr_args_release(args);
        args->quit = 1;
        goto exit;
      case 'i': args->input = optarg; break;
      case 'o': args->output = optarg; break;
      case 'v': args->verbose = 1; break;
      default: res = RES_BAD_ARG; break;
    }
  }
  if(res != RES_OK) {
    if(optarg) {
      fprintf(stderr, "%s: invalid option argumet '%s' -- '%c'\n",
        argv[0], optarg, opt);
    }
    goto error;
  }
  if(!args->input) {
    fprintf(stderr, "Missing input file.\n");
    res = RES_BAD_ARG;
    goto error;
  }
exit:
  return res;
error:
  htrdr_args_release(args);
  goto exit;
}

void
htrdr_args_release(struct htrdr_args* args)
{
  ASSERT(args);
  *args = HTRDR_ARGS_DEFAULT;
}

