/* Copyright (C) 2018-2019, 2022-2025 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2025 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2025 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2025 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2025 Observatoire de Paris
 * Copyright (C) 2022-2025 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2025 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2025 Université Paul Sabatier
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

#include "planets/htrdr_planets_args.h"

#include <rsys/cstr.h>
#include <rsys/stretchy_array.h>
#include <rsys/mem_allocator.h>

#include <getopt.h>
#include <string.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static INLINE res_T
check_gas_args(const struct rnatm_gas_args* args)
{
  if(!args) return RES_BAD_ARG;

  /* Filenames cannot be NULL */
  if(!args->smsh_filename
  || !args->sck_filename
  || !args->temperatures_filename)
    return RES_BAD_ARG;

  return RES_OK;
}

static INLINE res_T
check_aerosol_args(const struct rnatm_aerosol_args* args)
{
  if(!args) return RES_BAD_ARG;

  /* Filenames cannot be NULL */
  if(!args->smsh_filename
  || !args->sars_filename
  || !args->phase_fn_ids_filename
  || !args->phase_fn_lst_filename)
    return RES_BAD_ARG;

  return RES_OK;
}

static INLINE res_T
check_ground_args(const struct htrdr_planets_ground_args* args)
{
  if(!args) return RES_BAD_ARG;

  /* Filenames cannot be NULL */
  if(!args->smsh_filename
  || !args->props_filename
  || !args->mtllst_filename)
    return RES_BAD_ARG;

  return RES_OK;
}

static INLINE res_T
check_spectral_args(const struct htrdr_planets_spectral_args* args)
{
  if(!args) return RES_BAD_ARG;

  /* Invalid type */
  switch(args->type) {
    case HTRDR_SPECTRAL_LW:
    case HTRDR_SPECTRAL_SW:
    case HTRDR_SPECTRAL_SW_CIE_XYZ:
      /* Nothing to be done */
      break;
    default:
      return RES_BAD_ARG;
  }

  /* Invalid spectral range */
  if(args->wlen_range[0] < 0
  || args->wlen_range[1] < 0
  || args->wlen_range[0] > args->wlen_range[1])
    return RES_BAD_ARG;

  return RES_OK;
}

static INLINE res_T
check_volrad_budget_args(const struct htrdr_planets_volrad_budget_args* args)
{
  if(!args) return RES_BAD_ARG;

  /* Filename could not be NULL */
  if(!args->smsh_filename) return RES_BAD_ARG;

  /* Samples per tetrahedron could not be zero */
  if(!args->spt) return RES_BAD_ARG;

  return RES_OK;
}

static void
usage(void)
{
  printf("usage: htrdr-planets [-dfhNv] [-a aerosol_opt[:aerosol_opt ...]]\n");
  printf("                     [-C persp_camera_opt[:persp_camera_opt ...]]\n");
  printf("                     [-G ground_opt[:ground_opt ...]]\n");
  printf("                     [-i image_opt[:image_opt ...]] [-O accel_struct_storage]\n");
  printf("                     [-o output] [-r volrad_budget_opt[:volrad_budget_opt ...]]\n");
  printf("                     [-S source_opt[:source_opt ...]]\n");
  printf("                     [-s spectral_opt[:spectral_opt ...]] [-T optical_thickness]\n");
  printf("                     [-t threads_count] [-V accel_struct_definition]\n");
  printf("                     -g gas_opt[:gas_opt ...]\n");
}

static INLINE char*
str_dup(const char* str)
{
  size_t len = 0;
  char* dup = NULL;
  ASSERT(str);
  len = strlen(str) + 1/*NULL char*/;
  dup = mem_alloc(len);
  if(!dup) {
    return NULL;
  } else {
    return memcpy(dup, str, len);
  }
}

static res_T
parse_aerosol_parameters(const char* str, void* ptr)
{
  enum { MESH, NAME, RADPROP, PHASEFN, PHASEIDS } iparam;
  struct rnatm_aerosol_args* aerosol = NULL;
  char buf[BUFSIZ];
  struct htrdr_planets_args* args = ptr;
  char* key;
  char* val;
  char* tk_ctx;
  res_T res = RES_OK;
  ASSERT(args && str);

  if(strlen(str) >= sizeof(buf) -1/*NULL char*/) {
    fprintf(stderr, "Could not duplicate the aerosol parameter `%s'\n", str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  key = strtok_r(buf, "=", &tk_ctx);
  val = strtok_r(NULL, "", &tk_ctx);

       if(!strcmp(key, "mesh")) iparam = MESH;
  else if(!strcmp(key, "name")) iparam = NAME;
  else if(!strcmp(key, "radprop")) iparam = RADPROP;
  else if(!strcmp(key, "phasefn")) iparam = PHASEFN;
  else if(!strcmp(key, "phaseids")) iparam = PHASEIDS;
  else {
    fprintf(stderr, "Invalid aerosol parameter `%s'\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  if(!val) {
    fprintf(stderr, "Invalid null value for aerosol parameter `%s'\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  ASSERT(args->naerosols);
  aerosol = args->aerosols + (args->naerosols - 1);

  #define SET_STR(Dst) {                                                       \
    if(Dst) mem_rm(Dst);                                                       \
    if(!((Dst) = str_dup(val))) res = RES_MEM_ERR;                             \
  } (void)0
  switch(iparam) {
    case MESH: SET_STR(aerosol->smsh_filename); break;
    case NAME: SET_STR(aerosol->name); break;
    case RADPROP: SET_STR(aerosol->sars_filename); break;
    case PHASEFN: SET_STR(aerosol->phase_fn_lst_filename); break;
    case PHASEIDS: SET_STR(aerosol->phase_fn_ids_filename); break;
    default: FATAL("Unreachable code\n"); break;
  }
  #undef SET_STR
  if(res != RES_OK) {
    fprintf(stderr, "Unable to parse the aerosol parameter `%s' -- %s\n",
      str, res_to_cstr(res));
    goto error;
  }

exit:
  return res;
error:
  goto exit;
}

static res_T
parse_ground_parameters(const char* str, void* ptr)
{
  enum { BRDF, MESH, NAME, PROP } iparam;
  char buf[BUFSIZ];
  struct htrdr_planets_args* args = ptr;
  char* key;
  char* val;
  char* tk_ctx;
  res_T res = RES_OK;
  ASSERT(args && str);

  if(strlen(str) >= sizeof(buf) - 1/*NULL char*/) {
    fprintf(stderr, "Could not duplicate the ground parameter `%s'\n", str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  key = strtok_r(buf, "=", &tk_ctx);
  val = strtok_r(NULL, "", &tk_ctx);

       if(!strcmp(key, "brdf")) iparam = BRDF;
  else if(!strcmp(key, "mesh")) iparam = MESH;
  else if(!strcmp(key, "name")) iparam = NAME;
  else if(!strcmp(key, "prop")) iparam = PROP;
  else {
    fprintf(stderr, "Invalid ground parameter `%s'\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  if(!val) {
    fprintf(stderr, "Invalid null value for ground parameter `%s'\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  #define SET_STR(Dst) {                                                       \
    if(Dst) mem_rm(Dst);                                                       \
    if(!((Dst) = str_dup(val))) res = RES_MEM_ERR;                             \
  } (void)0
  switch(iparam) {
    case BRDF: SET_STR(args->ground.mtllst_filename); break;
    case MESH: SET_STR(args->ground.smsh_filename); break;
    case NAME: SET_STR(args->ground.name); break;
    case PROP: SET_STR(args->ground.props_filename); break;
    default: FATAL("Unreachable code\n"); break;
  }
  #undef SET_STR
  if(res != RES_OK) {
    fprintf(stderr, "Unable to parse the ground parameter `%s' -- %s\n",
      str, res_to_cstr(res));
    goto error;
  }

exit:
  return res;
error:
  goto exit;
}

static res_T
parse_gas_parameters(const char* str, void* ptr)
{
  enum { MESH, CK, TEMP } iparam;
  char buf[BUFSIZ];
  struct htrdr_planets_args* args = ptr;
  char* key;
  char* val;
  char* tk_ctx;
  res_T res = RES_OK;
  ASSERT(args && str);

  if(strlen(str) >= sizeof(buf) -1/*NULL char*/) {
    fprintf(stderr, "Could not duplicate the gas parameter `%s'\n", str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  key = strtok_r(buf, "=", &tk_ctx);
  val = strtok_r(NULL, "", &tk_ctx);

       if(!strcmp(key, "mesh")) iparam = MESH;
  else if(!strcmp(key, "ck")) iparam = CK;
  else if(!strcmp(key, "temp")) iparam = TEMP;
  else {
    fprintf(stderr, "Invalid gas parameter `%s'\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  if(!val) {
    fprintf(stderr, "Invalid null value for gas parameter `%s'\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  #define SET_STR(Dst) {                                                       \
    if(Dst) mem_rm(Dst);                                                       \
    if(!((Dst) = str_dup(val))) res = RES_MEM_ERR;                             \
  } (void)0
  switch(iparam) {
    case MESH: SET_STR(args->gas.smsh_filename); break;
    case CK: SET_STR(args->gas.sck_filename); break;
    case TEMP: SET_STR(args->gas.temperatures_filename); break;
    default: FATAL("Unreachable code\n"); break;
  }
  #undef SET_STR
  if(res != RES_OK) {
    fprintf(stderr, "Unable to parse the gas parameter `%s' -- %s\n",
      str, res_to_cstr(res));
    goto error;
  }

exit:
  return res;
error:
  goto exit;
}

static res_T
parse_source_parameters(const char* str, void* ptr)
{
  enum {LAT, LON, DST, RADIUS, TEMP, RAD} iparam;
  char buf[BUFSIZ];
  struct htrdr_planets_args* args = ptr;
  struct htrdr_planets_source_args* src = NULL;
  char* key;
  char* val;
  char* tk_ctx;
  res_T res = RES_OK;
  ASSERT(str && ptr);

  src = &args->source;

  if(strlen(str) >= sizeof(buf) -1/*NULL char*/) {
    fprintf(stderr, "Could not duplicate the source parameter `%s'\n", str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  key = strtok_r(buf, "=", &tk_ctx);
  val = strtok_r(NULL, "", &tk_ctx);

       if(!strcmp(key, "lat")) iparam = LAT;
  else if(!strcmp(key, "lon")) iparam = LON;
  else if(!strcmp(key, "dst")) iparam = DST;
  else if(!strcmp(key, "rad")) iparam = RAD;
  else if(!strcmp(key, "radius")) iparam = RADIUS;
  else if(!strcmp(key, "temp")) iparam = TEMP;
  else {
    fprintf(stderr, "Invalid source parameter `%s'\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  if(!val) {
    fprintf(stderr, "Invalid null value for the source parameter`%s'\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  switch(iparam) {
    case LAT:
      res = cstr_to_double(val, &src->latitude);
      if(res == RES_OK && (src->latitude < -90 || src->latitude > 90)) {
        res = RES_BAD_ARG;
      }
      break;
    case LON:
      res = cstr_to_double(val, &src->longitude);
      if(res == RES_OK && (src->longitude < -180 || src->longitude > 180)) {
        res = RES_BAD_ARG;
      }
      break;
    case DST:
      res = cstr_to_double(val, &src->distance);
      if(res == RES_OK && src->distance < 0) res = RES_BAD_ARG;
      break;
    case RAD:
      /* Use a per wavelength radiance rather than a constant temperature */
      src->temperature = -1;
      if(src->rnrl_filename) mem_rm(src->rnrl_filename);
      src->rnrl_filename = str_dup(val);
      if(!src->rnrl_filename) res = RES_MEM_ERR;
      break;
    case RADIUS:
      res = cstr_to_double(val, &src->radius);
      if(res == RES_OK && src->radius < 0) res = RES_BAD_ARG;
      break;
    case TEMP:
      /* Use a constant temperature rather than a per wavelength radiance */
      if(src->rnrl_filename) {
        mem_rm(src->rnrl_filename);
        src->rnrl_filename = NULL;
      }
      res = cstr_to_double(val, &src->temperature);
      if(res == RES_OK && src->temperature < 0) res = RES_BAD_ARG;
      break;
    default: FATAL("Unreachable code\n"); break;
  }
  if(res != RES_OK) {
    fprintf(stderr, "Unable to parse the source parameter `%s' -- %s\n",
      str, res_to_cstr(res));
    goto error;
  }

exit:
  return res;
error:
  goto exit;
}

static INLINE res_T
parse_spectral_range(const char* str, double wlen_range[2])
{
  double range[2];
  size_t len;
  res_T res = RES_OK;
  ASSERT(wlen_range && str);

  res = cstr_to_list_double(str, ',', range, &len, 2);
  if(res == RES_OK && len != 2) res = RES_BAD_ARG;
  if(res == RES_OK && range[0] > range[1]) res = RES_BAD_ARG;
  if(res == RES_OK && (range[0] < 0 || range[1] < 0)) res = RES_BAD_ARG;
  if(res != RES_OK) goto error;

  wlen_range[0] = range[0];
  wlen_range[1] = range[1];

exit:
  return res;
error:
  goto exit;
}

static res_T
parse_spectral_parameters(const char* str, void* ptr)
{
  enum {CIE_XYZ, LW, SW} iparam;
  char buf[BUFSIZ];
  struct htrdr_planets_args* args = ptr;
  struct htrdr_planets_spectral_args* spectral = NULL;
  char* key;
  char* val;
  char* tk_ctx;
  res_T res = RES_OK;
  ASSERT(str && ptr);

  spectral = &args->spectral_domain;

  if(strlen(str) >= sizeof(buf) -1/*NULL char*/) {
    fprintf(stderr, "Could not duplicate the spectral parameter `%s'\n", str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  key = strtok_r(buf, "=", &tk_ctx);
  val = strtok_r(NULL, "", &tk_ctx);

       if(!strcmp(key, "cie_xyz")) iparam = CIE_XYZ;
  else if(!strcmp(key, "lw")) iparam = LW;
  else if(!strcmp(key, "sw")) iparam = SW;
  else {
    fprintf(stderr, "Invalid spectral parameter `%s'\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  if((iparam == LW || iparam == SW) && !val) {
    fprintf(stderr,
      "Invalid null value for the spectral parameter `%s'\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  switch(iparam) {
    case CIE_XYZ:
      spectral->type = HTRDR_SPECTRAL_SW_CIE_XYZ;
      spectral->wlen_range[0] = HTRDR_RAN_WLEN_CIE_XYZ_RANGE_DEFAULT[0];
      spectral->wlen_range[1] = HTRDR_RAN_WLEN_CIE_XYZ_RANGE_DEFAULT[1];
      break;
    case LW:
      spectral->type = HTRDR_SPECTRAL_LW;
      res = parse_spectral_range(val, spectral->wlen_range);
      break;
    case SW:
      spectral->type = HTRDR_SPECTRAL_SW;
      res = parse_spectral_range(val, spectral->wlen_range);
      break;
    default: FATAL("Unreachable code\n"); break;
  }
  if(res != RES_OK) {
    fprintf(stderr, "Unable to parse the spectral parameter `%s' -- %s\n",
      str, res_to_cstr(res));
    goto error;
  }

exit:
  return res;
error:
  goto exit;
}

static res_T
parse_volrad_budget_parameters(const char* str, void* ptr)
{
  enum { MESH, SPT } iparam;
  char buf[BUFSIZ];
  struct htrdr_planets_args* args = ptr;
  char* key;
  char* val;
  char* tk_ctx;
  res_T res = RES_OK;
  ASSERT(str && ptr);

  if(strlen(str) >= sizeof(buf) -1/*NULL char*/) {
    fprintf(stderr,
      "Could not duplicate the parameters "
      "of the volumic radiative budget calculation `%s'\n",
      str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  key = strtok_r(buf, "=", &tk_ctx);
  val = strtok_r(NULL, "", &tk_ctx);

       if(!strcmp(key, "mesh")) iparam = MESH;
  else if(!strcmp(key, "spt")) iparam = SPT;
  else {
    fprintf(stderr, "Invalid volumic radiative budget parameter `%s'\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  if(!val) {
    fprintf(stderr,
      "Invalid null value for the volumic radiative budget parameter `%s'.\n",
      key);
    res = RES_BAD_ARG;
    goto error;
  }

  switch(iparam) {
    case MESH:
      if(args->volrad_budget.smsh_filename) {
        mem_rm(args->volrad_budget.smsh_filename);
      }
      if(!(args->volrad_budget.smsh_filename = str_dup(val))) {
        res = RES_MEM_ERR;
      }
      break;
    case SPT: /* Sample Per Tetrahedron */
      res = cstr_to_uint(val, &args->volrad_budget.spt);
      break;
    default: FATAL("Unreachable code\n"); break;
  }
  if(res != RES_OK) {
    fprintf(stderr,
      "Unable to parse the volumic radiative budget parameter `%s' -- %s\n",
      str, res_to_cstr(res));
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
htrdr_planets_args_init(struct htrdr_planets_args* args, int argc, char** argv)
{
  int opt;
  res_T res = RES_OK;
  ASSERT(args && argc && argv);

  *args = HTRDR_PLANETS_ARGS_DEFAULT;

  while((opt = getopt(argc, argv, "a:C:dfG:g:hi:NO:o:r:S:s:T:t:V:v")) != -1) {
    switch(opt) {
      case 'a':
        (void)sa_add(args->aerosols, 1);
        args->aerosols[args->naerosols] = RNATM_AEROSOL_ARGS_NULL;
        args->naerosols += 1;
        res = cstr_parse_list(optarg, ':', parse_aerosol_parameters, args);
        if(res == RES_OK) {
          res = check_aerosol_args(args->aerosols+args->naerosols-1);
        }
        break;
      case 'C':
        res = htrdr_args_camera_perspective_parse(&args->cam_persp, optarg);
        args->output_type = HTRDR_PLANETS_ARGS_OUTPUT_IMAGE;
        break;
      case 'd':
        args->output_type = HTRDR_PLANETS_ARGS_OUTPUT_OCTREES;
        break;
      case 'f':
        args->force_output_overwrite = 1;
        break;
      case 'G':
        res = cstr_parse_list(optarg, ':', parse_ground_parameters, args);
        if(res == RES_OK) {
          res = check_ground_args(&args->ground);
        }
        break;
      case 'g':
        res = cstr_parse_list(optarg, ':', parse_gas_parameters, args);
        if(res == RES_OK) {
          res = check_gas_args(&args->gas);
        }
        break;
      case 'h':
        usage();
        htrdr_planets_args_release(args);
        args->quit = 1;
        goto exit;
      case 'i':
        res = htrdr_args_image_parse(&args->image, optarg);
        break;
      case 'N': args->precompute_normals = 1; break;
      case 'O': args->octrees_storage = optarg; break;
      case 'o': args->output = optarg; break;
      case 'r':
        res = cstr_parse_list(optarg, ':', parse_volrad_budget_parameters, args);
        args->output_type = HTRDR_PLANETS_ARGS_OUTPUT_VOLUMIC_RADIATIVE_BUDGET;
        break;
      case 'S':
        res = cstr_parse_list(optarg, ':', parse_source_parameters, args);
        break;
      case 's':
        res = cstr_parse_list(optarg, ':', parse_spectral_parameters, args);
        break;
      case 'T':
        res = cstr_to_double(optarg, &args->optical_thickness);
        if(res == RES_OK && args->optical_thickness < 0) res = RES_BAD_ARG;
        break;
      case 't':
        res = cstr_to_uint(optarg, &args->nthreads);
        if(res == RES_OK && !args->nthreads) res = RES_BAD_ARG;
        break;
      case 'V':
        res = cstr_to_uint(optarg, &args->octree_definition_hint);
        if(res == RES_OK && !args->octree_definition_hint) res = RES_BAD_ARG;
        break;
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

  res = check_gas_args(&args->gas);
  if(res != RES_OK) {
    fprintf(stderr, "missing gas definition -- option '-g'\n");
    goto error;
  }

  if(args->output_type != HTRDR_PLANETS_ARGS_OUTPUT_OCTREES) {
    res = check_ground_args(&args->ground);
    if(res != RES_OK) {
      fprintf(stderr, "missing ground definition -- option '-G'\n");
      goto error;
    }

    /* Check the source */
    if(args->spectral_domain.type == HTRDR_SPECTRAL_SW
    || args->spectral_domain.type == HTRDR_SPECTRAL_SW_CIE_XYZ) {
      res = htrdr_planets_source_args_check(&args->source);
      if(res != RES_OK) {
        fprintf(stderr, "missing source definition -- option '-S'\n");
        goto error;
      }
    }
  }

exit:
  return res;
error:
  usage();
  htrdr_planets_args_release(args);
  goto exit;
}

void
htrdr_planets_args_release(struct htrdr_planets_args* args)
{
  size_t i;
  ASSERT(args);

  if(args->gas.smsh_filename) mem_rm(args->gas.smsh_filename);
  if(args->gas.sck_filename) mem_rm(args->gas.sck_filename);
  if(args->gas.temperatures_filename) mem_rm(args->gas.temperatures_filename);
  if(args->ground.smsh_filename) mem_rm(args->ground.smsh_filename);
  if(args->ground.props_filename) mem_rm(args->ground.props_filename);
  if(args->ground.mtllst_filename) mem_rm(args->ground.mtllst_filename);
  if(args->ground.name) mem_rm(args->ground.name);
  if(args->source.rnrl_filename) mem_rm(args->source.rnrl_filename);
  if(args->volrad_budget.smsh_filename) mem_rm(args->volrad_budget.smsh_filename);

  FOR_EACH(i, 0, args->naerosols) {
    struct rnatm_aerosol_args* aerosol = args->aerosols + i;
    if(aerosol->name) mem_rm(aerosol->name);
    if(aerosol->smsh_filename) mem_rm(aerosol->smsh_filename);
    if(aerosol->sars_filename) mem_rm(aerosol->sars_filename);
    if(aerosol->phase_fn_ids_filename) mem_rm(aerosol->phase_fn_ids_filename);
    if(aerosol->phase_fn_lst_filename) mem_rm(aerosol->phase_fn_lst_filename);
  }
  sa_release(args->aerosols);

  *args = HTRDR_PLANETS_ARGS_DEFAULT;
}

res_T
htrdr_planets_args_check(const struct htrdr_planets_args* args)
{
  size_t i;
  res_T res = RES_OK;

  if(!args) return RES_BAD_ARG;

  /* Check the gas */
  res = check_gas_args(&args->gas);
  if(res != RES_OK) return res;

  /* Check the aerosols */
  FOR_EACH(i, 0, args->naerosols) {
    res = check_aerosol_args(args->aerosols+i);
    if(res != RES_OK) return res;
  }

  /* Check the octree parameters */
  if(args->octree_definition_hint == 0
  || args->optical_thickness < 0)
    return RES_BAD_ARG;

  /* Check the spectral domain */
  res = check_spectral_args(&args->spectral_domain);
  if(res != RES_OK) return res;

  if(args->output_type != HTRDR_PLANETS_ARGS_OUTPUT_OCTREES) {
    /* Check the ground */
    res = check_ground_args(&args->ground);
    if(res != RES_OK) return res;

    /* Check the source */
    if(args->spectral_domain.type == HTRDR_SPECTRAL_SW
    || args->spectral_domain.type == HTRDR_SPECTRAL_SW_CIE_XYZ) {
      res = htrdr_planets_source_args_check(&args->source);
      if(res != RES_OK) return res;
    }
  }

  if(args->output_type == HTRDR_PLANETS_ARGS_OUTPUT_IMAGE) {
    res = htrdr_args_camera_perspective_check(&args->cam_persp);
    if(res != RES_OK) return res;

    res = htrdr_args_image_check(&args->image);
    if(res != RES_OK) return res;
  }

  if(args->output_type == HTRDR_PLANETS_ARGS_OUTPUT_VOLUMIC_RADIATIVE_BUDGET) {
    res = check_volrad_budget_args(&args->volrad_budget);
    if(res != RES_OK) return res;
  }

  /* Check miscalleneous parameters */
  if(args->nthreads == 0
  || (unsigned)args->output_type >= HTRDR_PLANETS_ARGS_OUTPUT_TYPES_COUNT__)
    return RES_BAD_ARG;

  return RES_OK;
}

res_T
htrdr_planets_source_args_check(const struct htrdr_planets_source_args* args)
{
  if(!args) return RES_BAD_ARG;

  /* Invalid position */
  if(args->latitude <-90
  || args->latitude > 90
  || args->longitude <-180
  || args->longitude > 180
  || args->distance < 0)
    return RES_BAD_ARG;

  /* Invalid radius */
  if(args->radius < 0)
    return RES_BAD_ARG;

  /* Invalid radiance */
  if((args->temperature < 0 && !args->rnrl_filename) /* Both are invalids */
  || (args->temperature >=0 &&  args->rnrl_filename)) /* Both are valids */
    return RES_BAD_ARG;

  return RES_OK;
}
