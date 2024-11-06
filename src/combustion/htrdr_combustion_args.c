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

#define _POSIX_C_SOURCE 200112L /* strtok_r support */

#include "combustion/htrdr_combustion_args.h"

#include <rsys/cstr.h>

#include <getopt.h>
#include <string.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
usage(void)
{
  printf("usage: htrdr-combustion [-fhINsv] [-C persp_camera_opt[:persp_camera_opt ...]]\n");
  printf("                        [-D laser_flux_density] [-d dump_type]\n");
  printf("                        [-F rdgfa_opt[:rdgfa_opt ...]]\n");
  printf("                        [-g combustion_chamber_opt[:combustion_chamber_opt...]]\n");
  printf("                        [-i image_opt[:image_opt ...]]\n");
  printf("                        [-l laser_opt[:laser_opt ...]] [-O cache] [-o output]\n");
  printf("                        [-P ortho_camera_opt[:ortho_camera_opt ...]]\n");
  printf("                        [-R flux_sensor_opt[:flux_sensor_opt ...]]\n");
  printf("                        [-T optical_thickness] [-t threads_count]\n");
  printf("                        [-V accel_struct_definition] [-w laser_wavelength]\n");
  printf("                        -m medium_geometry -p thermo_properties\n");
  printf("                        -r refractive_ids\n");
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

static res_T
parse_fractal_parameters(const char* str, void* ptr)
{
  char buf[128];
  struct htrdr_combustion_args* args = ptr;
  char* key;
  char* val;
  char* ctx;
  res_T res = RES_OK;
  ASSERT(ptr && str);

  if(strlen(str) >= sizeof(buf) -1/*NULL char*/) {
    fprintf(stderr,
      "Could not duplicate the fractal option string `%s'.\n", str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  key = strtok_r(buf, "=", &ctx);
  val = strtok_r(NULL, "",  &ctx);

  if(!strcmp(key, "prefactor")) {
    res = cstr_to_double(val, &args->fractal_prefactor);
    if(res != RES_OK) goto error;
  } else if(!strcmp(key, "dimension")) {
    res = cstr_to_double(val, &args->fractal_dimension);
    if(res != RES_OK) goto error;
  } else {
    fprintf(stderr, "Invalid fractal parameter `%s'.\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

exit:
  return res;
error:
  goto exit;
}

static res_T
parse_dump_parameter
  (const char* str,
   enum htrdr_combustion_args_output_type* output_type)
{
  res_T res = RES_OK;
  ASSERT(str && output_type);

  if(!strcmp(str, "octree")) {
    *output_type = HTRDR_COMBUSTION_ARGS_OUTPUT_OCTREES;
  } else if(!strcmp(str, "laser")) {
    *output_type = HTRDR_COMBUSTION_ARGS_OUTPUT_LASER_SHEET;
  } else {
    fprintf(stderr, "Invalid dump parameter `%s'.\n", str);
    res = RES_BAD_ARG;
    goto error;
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

  while((opt = getopt(argc, argv, "C:D:d:F:fg:hIi:l:m:NO:o:P:p:R:r:sT:t:V:vw:")) != -1) {
    switch(opt) {
      case 'C':
        args->output_type = HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE;
        args->cam_type = HTRDR_ARGS_CAMERA_PERSPECTIVE;
        res = htrdr_args_camera_perspective_parse(&args->cam_persp, optarg);
        break;
      case 'D':
        res = cstr_to_double(optarg, &args->laser_flux_density);
        if(res == RES_OK && args->laser_flux_density <= 0) res = RES_BAD_ARG;
        break;
      case 'd':
        res = parse_dump_parameter(optarg, &args->output_type);
        break;
      case 'F':
        res = cstr_parse_list(optarg, ':', parse_fractal_parameters, args);
        args->phase_func_type = HTRDR_COMBUSTION_ARGS_PHASE_FUNC_RDGFA;
        break;
      case 'f': args->force_overwriting = 1; break;
      case 'g':
        res = htrdr_args_geometry_parse(&args->geom, optarg);
        break;
      case 'h':
        usage();
        htrdr_combustion_args_release(args);
        args->quit = 1;
        goto exit;
      case 'I':
        args->phase_func_type = HTRDR_COMBUSTION_ARGS_PHASE_FUNC_ISOTROPIC;
        break;
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
      case 'P':
        args->output_type = HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE;
        args->cam_type = HTRDR_ARGS_CAMERA_ORTHOGRAPHIC;
        res = htrdr_args_camera_orthographic_parse(&args->cam_ortho, optarg);
        break;
      case 'r': args->path_refract_ids = optarg; break;
      case 'R':
        args->output_type = HTRDR_COMBUSTION_ARGS_OUTPUT_FLUX_MAP;
        res = htrdr_args_rectangle_parse(&args->flux_map, optarg);
        break;
      case 's': args->use_simd = 1; break;
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
    fprintf(stderr, "missing the volumetric mesh -- option '-m'\n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(!args->path_therm_props) {
    fprintf(stderr, "missing the thermodynamic properties -- option '-p'\n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(!args->path_refract_ids) {
    fprintf(stderr, "missing the refractive indices -- option '-r'\n");
    res = RES_BAD_ARG;
    goto error;
  }

exit:
  return res;
error:
  usage();
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

