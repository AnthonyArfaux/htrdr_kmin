/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2023 Observatoire de Paris
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

#include "atmosphere/htrdr_atmosphere_args.h"

#include <rsys/cstr.h>

#include <getopt.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
print_help(const char* cmd)
{
  ASSERT(cmd);
  printf("Usage: %s [option] ... -a gas\n", cmd);
  printf(
"Simulate radiative transfer in a plane-parallel atmosphere.\n"
"See htrdr-atmosphere(1) man page for details\n\n");

  printf(
"  -a gas         filename of the gas optical properties\n");
  printf(
"  -C camera      configure a perspective camera\n");
  printf(
"  -c clouds      filename of the clouds properties\n");
  printf(
"  -D azimuth,elevation\n"
"                 direction in degrees toward the sun center. By default\n"
"                 azimuth is %g and elevation is %g\n",
    HTRDR_ATMOSPHERE_ARGS_DEFAULT.sun_azimuth,
    HTRDR_ATMOSPHERE_ARGS_DEFAULT.sun_elevation);
  printf(
"  -d             dump volumetric acceleration structures to output\n"
"                 and exit\n");
  printf(
"  -f             overwrite the output file if it already exists\n");
  printf(
"  -g ground      filename of the ground geometry\n");
  printf(
"  -h             display this help and exit\n");
  printf(
"  -i image       image to compute\n");
  printf(
"  -M materials   filename of the ground materials\n");
  printf(
"  -m mie         filename of the Mie's data\n");
  printf(
"  -n sky-name    name used to identify the sky in the materials file.\n"
"                 Its default value is `%s'\n",
    HTRDR_ATMOSPHERE_ARGS_DEFAULT.sky_mtl_name);
  printf(
"  -O cache       filenaname of the cache file used to store/restore the\n"
"                 volumetric data. By default do not use any cache\n");
  printf(
"  -o output      file where data are written. If not defined, data are\n"
"                 written to standard output\n");
  printf(
"  -p rectangle   switch in flux computation by defining the rectangular\n"
"                 sensor onto which the flux is computed\n");
  printf(
"  -P camera      configure an orthoraphic camera\n");
  printf(
"  -R             infinitely repeat the ground along the X and Y axis\n");
  printf(
"  -r             infinitely repeat the clouds along the X and Y axis\n");
  printf(
"  -s spectral    define the spectral doamin of integration\n");
  printf(
"  -T optical_thickness\n"
"                 optical thickness criteria for octree building.\n"
"                 Default is %g\n",
    HTRDR_ATMOSPHERE_ARGS_DEFAULT.optical_thickness);
  printf(
"  -t threads     hint on the number of threads to use.\n"
"                 Default assumes as many threads as CPU cores\n");
  printf(
"  -V octree_definition\n"
"                 advice on the definition of the atmospheric\n"
"                 acceleration structures. By default use\n"
"                 the definition of the clouds data\n");
  printf(
"  -v             make the command verbose\n");
  printf("\n");
  htrdr_fprint_license(cmd, stdout);
}

static res_T
parse_grid_definition(struct htrdr_atmosphere_args* args, const char* str)
{
  unsigned def[3];
  size_t len;
  res_T res = RES_OK;
  ASSERT(args && str);

  res = cstr_to_list_uint(str, ',', def, &len, 3);
  if(res == RES_OK && len != 3) res = RES_BAD_ARG;
  if(res != RES_OK) {
    fprintf(stderr, "Invalid grid definition `%s'.\n", str);
    goto error;
  }

  if(!def[0] || !def[1] || !def[2]) {
    fprintf(stderr,
      "Invalid null grid definition {%u, %u, %u}.\n", SPLIT3(def));
    res = RES_BAD_ARG;
    goto error;
  }

  args->grid_max_definition[0] = def[0];
  args->grid_max_definition[1] = def[1];
  args->grid_max_definition[2] = def[2];

exit:
  return res;
error:
  goto exit;
}

static res_T
parse_sun_dir(struct htrdr_atmosphere_args* args, const char* str)
{
  double angles[2];
  size_t len;
  res_T res = RES_OK;
  ASSERT(args && str);

  res = cstr_to_list_double(str, ',', angles, &len, 2);
  if(res == RES_OK && len != 2) res = RES_BAD_ARG;
  if(res != RES_OK) {
    fprintf(stderr, "Invalid direction `%s'.\n", str);
    goto error;
  }

  if(angles[0] < 0 || angles[0] >= 360) {
    fprintf(stderr,
      "Invalid azimuth angle `%g'. Azimuth must be in [0, 360[ degrees.\n",
      angles[0]);
    res = RES_BAD_ARG;
    goto error;
  }

  if(angles[1] < 0 || angles[1] > 90) {
    fprintf(stderr,
      "Invalid elevation angle `%g'. Elevation must be in [0, 90] degrees.\n",
      angles[1]);
    res = RES_BAD_ARG;
    goto error;
  }

  args->sun_azimuth = angles[0];
  args->sun_elevation = angles[1];

exit:
  return res;
error:
  goto exit;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_atmosphere_args_init
  (struct htrdr_atmosphere_args* args,
   int argc,
   char** argv)
{
  int opt;
  res_T res = RES_OK;
  ASSERT(args && argc && argv);

  *args = HTRDR_ATMOSPHERE_ARGS_DEFAULT;

  while((opt = getopt(argc, argv, "a:C:c:D:dfg:hi:M:m:n:O:o:P:p:Rrs:T:t:V:v")) != -1) {
    switch(opt) {
      case 'a': args->filename_gas = optarg; break;
      case 'C':
        args->output_type = HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE;
        args->cam_type = HTRDR_ARGS_CAMERA_PERSPECTIVE;
        res = htrdr_args_camera_perspective_parse(&args->cam_persp, optarg);
        break;
      case 'c': args->filename_les = optarg; break;
      case 'D': res = parse_sun_dir(args, optarg); break;
      case 'd':
        args->output_type = HTRDR_ATMOSPHERE_ARGS_OUTPUT_OCTREES;
        break;
      case 'f': args->force_overwriting = 1; break;
      case 'g': args->filename_obj = optarg; break;
      case 'h':
        print_help(argv[0]);
        htrdr_atmosphere_args_release(args);
        args->quit = 1;
        goto exit;
      case 'i':
        res = htrdr_args_image_parse(&args->image, optarg);
        break;
      case 'M': args->filename_mtl = optarg; break;
      case 'm': args->filename_mie = optarg; break;
      case 'n': args->sky_mtl_name = optarg; break;
      case 'O': args->filename_cache = optarg; break;
      case 'o': args->filename_output = optarg; break;
      case 'p':
        args->output_type = HTRDR_ATMOSPHERE_ARGS_OUTPUT_FLUX_MAP;;
        res = htrdr_args_rectangle_parse(&args->flux_map, optarg);
        break;
      case 'P':
        args->output_type = HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE;
        args->cam_type = HTRDR_ARGS_CAMERA_ORTHOGRAPHIC;
        res = htrdr_args_camera_orthographic_parse(&args->cam_ortho, optarg);
        break;
      case 'r': args->repeat_clouds = 1; break;
      case 'R': args->repeat_ground = 1; break;
      case 's':
        res = htrdr_args_spectral_parse(&args->spectral, optarg);
        break;
      case 'T':
        res = cstr_to_double(optarg, &args->optical_thickness);
        if(res == RES_OK && args->optical_thickness < 0) res = RES_BAD_ARG;
        break;
      case 't': /* Submit an hint on the number of threads to use */
        res = cstr_to_uint(optarg, &args->nthreads);
        if(res == RES_OK && !args->nthreads) res = RES_BAD_ARG;
        break;
      case 'V': res = parse_grid_definition(args, optarg); break;
      case 'v': args->verbose = 1; break;
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
  if(!args->filename_gas) {
    fprintf(stderr,
      "Missing the path of the gas optical properties file -- option '-a'\n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(args->filename_obj && !args->filename_mtl) {
    fprintf(stderr,
      "Missing the path of the file listing the ground materials -- option '-M'\n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(args->filename_les && !args->filename_mie) {
    fprintf(stderr,
      "Missing the path toward the file of the Mie's data -- option '-m'\n");
    res = RES_BAD_ARG;
    goto error;
  }

  /* Setup default ref temperature if necessary */
  if(args->spectral.ref_temperature <= 0) {
    switch(args->spectral.spectral_type) {
      case HTRDR_SPECTRAL_LW:
        args->spectral.ref_temperature = HTRDR_DEFAULT_LW_REF_TEMPERATURE;
        break;
      case HTRDR_SPECTRAL_SW:
        args->spectral.ref_temperature = HTRDR_SUN_TEMPERATURE;
        break;
      case HTRDR_SPECTRAL_SW_CIE_XYZ:
        args->spectral.ref_temperature = -1; /* Unused */
        break;
      default: FATAL("Unreachable code.\n"); break;
    }
  }

exit:
  return res;
error:
  htrdr_atmosphere_args_release(args);
  goto exit;
}

void
htrdr_atmosphere_args_release(struct htrdr_atmosphere_args* args)
{
  ASSERT(args);
  *args = HTRDR_ATMOSPHERE_ARGS_DEFAULT;
}
