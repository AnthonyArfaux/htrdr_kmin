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

#include <rsys/cstr.h>
#include <getopt.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
print_help(const char* cmd)
{
  ASSERT(cmd);
  printf(
"Usage: %s [<options>] -m TETRAHEDRA -p THERMOPROPS -r REFRACT_IDS\n",
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
"  -g <geometry>  define the combustion chamber geometry. Refer to the\n"
"                 %s man page for the list of geometry options.\n", cmd);
  printf(
"  -h             display this help and exit.\n");
  printf(
"  -i <image>     define the image to compute. Refer to the %s man\n"
"                 page for the list of image options.\n", cmd);
  printf(
"  -l <laser>     define the geometry of the laser sheet. Refer to the\n"
"                 %s man page for the list of laser options.\n", cmd);
  printf(
"  -m TETRAHEDRA  path toward the volumetric mesh.\n");
  printf(
"  -N             precompute the tetrahedra normals.\n");
  printf(
"  -O CACHE       path of the cache file used to store/restore the\n"
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
"  -V <HINT|X,Y,Z>\n"
"                 definition of the volumetric acceleration grids along\n"
"                 the 3 axis. By default it is computed automatically\n"
"                 with a hint on the expected definition set to %u.\n",
    HTRDR_COMBUSTION_ARGS_DEFAULT.grid.definition.hint);
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

static res_T
parse_grid_definition
  (struct htrdr_combustion_args_grid_definition* grid_def,
   const char* str)
{
  unsigned def[3];
  size_t len;
  res_T res = RES_OK;
  ASSERT(grid_def && str);

  res = cstr_to_list_uint(str, ',', def, &len, 3);
  if(res == RES_OK && len == 2) res = RES_BAD_ARG;
  if(res != RES_OK) {
    fprintf(stderr, "Invalid grid definition `%s'.\n", str);
    goto error;
  }

  if(len == 1) {
    if(!def[0]) {
      fprintf(stderr, "Invalid null grid definition %u.\n", def[0]);
      res = RES_BAD_ARG;
      goto error;
    }
    grid_def->type = HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_AUTO;
    grid_def->definition.hint = def[0];

  } else {
    if(!def[0] || !def[1] || !def[2]) {
      fprintf(stderr,
        "Invalid null grid definition [%u, %u, %u].\n", SPLIT3(def));
      res = RES_BAD_ARG;
      goto error;
    }
    grid_def->type = HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_FIXED;
    grid_def->definition.fixed[0] = def[0];
    grid_def->definition.fixed[1] = def[1];
    grid_def->definition.fixed[2] = def[2];
  }

exit:
  return res;
error:
  goto exit;
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
      case 'C':
        res = htrdr_args_camera_parse(&args->camera, optarg);
        break;
      case 'd':
        args->dump_volumetric_acceleration_structure = 1;
        break;
      case 'F':
        res = cstr_to_double(optarg, &args->fractal_dimension);
        if(res == RES_OK && args->fractal_dimension <= 0)
          res = RES_BAD_ARG;
        break;
      case 'f': args->force_overwriting = 1; break;
      case 'G':
        res = cstr_to_double(optarg, &args->gyration_radius_prefactor);
        if(res == RES_OK && args->gyration_radius_prefactor <= 0)
          res = RES_BAD_ARG;
        break;
      case 'g':
        res = htrdr_args_geometry_parse(&args->geom, optarg);
        break;
      case 'h':
        print_help(argv[0]);
        htrdr_combustion_args_release(args);
        args->quit = 1;
        goto exit;
      case 'i':
        res = htrdr_args_image_parse(&args->image, optarg);
        break;
      case 'l':
        res = htrdr_args_rectangle_parse(&args->laser, optarg);
        break;
      case 'm': args->path_tetra = optarg; break;
      case 'N': args->precompute_normals = 1; break;
      case 'O': args->path_cache = optarg; break;
      case 'o': args->path_output = optarg; break;
      case 'p': args->path_therm_props = optarg; break;
      case 'r': args->path_refract_ids = optarg; break;
      case 'T':
        res = cstr_to_double(optarg, &args->optical_thickness);
        if(res == RES_OK && args->optical_thickness < 0) res = RES_BAD_ARG;
        break;
      case 't': /* Submit an hint on the number of threads to use */
        res = cstr_to_uint(optarg, &args->nthreads);
        if(res == RES_OK && !args->nthreads) res = RES_BAD_ARG;
        break;
      case 'V':
        res = parse_grid_definition(&args->grid, optarg);
        break;
      case 'v': args->verbose = 1; break;
      case 'w':
        res = cstr_to_double(optarg, &args->wavelength);
        break;
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

  if(!args->path_tetra) {
    fprintf(stderr, "Missing the volumetric mesh -- option '-m'\n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(!args->path_therm_props) {
    fprintf(stderr, "Missing the thermodynamic properties -- option '-p'\n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(!args->path_refract_ids) {
    fprintf(stderr, "Missing the refractive indices -- option '-r'\n");
    res = RES_BAD_ARG;
    goto error;
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
  htrdr_args_geometry_free(&args->geom);
  *args = HTRDR_COMBUSTION_ARGS_DEFAULT;
}

