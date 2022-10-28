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
check_ground_args(const struct htrdr_planeto_ground_args* args)
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
check_spectral_args(const struct htrdr_args_spectral* spectral_domain)
{
  if(!spectral_domain) return RES_BAD_ARG;

  switch(spectral_domain->spectral_type) {
    case HTRDR_SPECTRAL_LW:
    case HTRDR_SPECTRAL_SW:

      /* Invalid reference temperature */
      if(spectral_domain->ref_temperature <= 0)
        return RES_BAD_ARG;

      /* Invalid spectral range */
      if(spectral_domain->wlen_range[0]
      >  spectral_domain->wlen_range[1])
        return RES_BAD_ARG;

      break;
    case HTRDR_SPECTRAL_SW_CIE_XYZ:
      /* Nothing to check since all parameters are implicitly defined */
      break;
    default: FATAL("Unreachable code\n"); break;
  }

  return RES_OK;
}

static void
print_help(const char* cmd)
{
  ASSERT(cmd);
  printf(
"Usage: %s [-dfhv] [-s spectral_domain] [-t threads]\n"
"                     [-T optical_thickness] [-V octree_definition]\n"
"                     [-O octrees_storage] [-o output]\n"
"                     [-a aerosol]... -g gas -G ground -S source\n", cmd);
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
"  -i image       image to compute\n");
  printf(
"  -O octrees_storage\n"
"                 file where atmospheric acceleration structures are\n"
"                 stored/loaded\n");
  printf(
"  -o output      file where the result is written. If not defined,\n"
"                 the result is written to standard output\n");
  printf(
"  -S source      define the source\n");
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
"                 Default assumes as many threads as CPU cores\n");
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
  struct htrdr_planeto_args* args = ptr;
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

  switch(iparam) {
    case MESH:
      aerosol->smsh_filename = str_dup(val);
      if(!aerosol->smsh_filename) res = RES_MEM_ERR;
      break;
    case NAME:
      aerosol->name = str_dup(val);
      if(!aerosol->name) res = RES_MEM_ERR;
      break;
    case RADPROP:
      aerosol->sars_filename = str_dup(val);
      if(!aerosol->sars_filename) res = RES_MEM_ERR;
      break;
    case PHASEFN:
      aerosol->phase_fn_lst_filename = str_dup(val);
      if(!aerosol->phase_fn_lst_filename) res = RES_MEM_ERR;
      break;
    case PHASEIDS:
      aerosol->phase_fn_ids_filename = str_dup(val);
      if(!aerosol->phase_fn_ids_filename) res = RES_MEM_ERR;
      break;
    default: FATAL("Unreachable code\n"); break;
  }
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
  struct htrdr_planeto_args* args = ptr;
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

  switch(iparam) {
    case BRDF:
      args->ground.mtllst_filename = str_dup(val);
      if(!args->ground.mtllst_filename) res = RES_MEM_ERR;
      goto error;
    case MESH:
      args->ground.smsh_filename = str_dup(val);
      if(!args->ground.smsh_filename) res = RES_MEM_ERR;
      break;
    case NAME:
      args->ground.name = str_dup(val);
      if(!args->ground.name) res = RES_MEM_ERR;
      break;
    case PROP:
      args->ground.props_filename = str_dup(val);
      if(!args->ground.props_filename) res = RES_MEM_ERR;
      break;
    default: FATAL("Unreachable code\n"); break;
  }
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
  struct htrdr_planeto_args* args = ptr;
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

  switch(iparam) {
    case MESH:
      args->gas.smsh_filename = str_dup(val);
      if(!args->gas.smsh_filename) res = RES_MEM_ERR;
      break;
    case CK:
      args->gas.sck_filename = str_dup(val);
      if(!args->gas.sck_filename) res = RES_MEM_ERR;
      break;
    case TEMP:
      args->gas.temperatures_filename = str_dup(val);
      if(!args->gas.temperatures_filename) res = RES_MEM_ERR;
      break;
    default: FATAL("Unreachable code\n"); break;
  }
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
  enum {LAT, LON, DST, RADIUS, TEMP} iparam;
  char buf[BUFSIZ];
  struct htrdr_planeto_args* args = ptr;
  struct htrdr_planeto_source_args* src = NULL;
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
    case RADIUS:
      res = cstr_to_double(val, &src->radius);
      if(res == RES_OK && src->radius < 0) res = RES_BAD_ARG;
      break;
    case TEMP:
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

  while((opt = getopt(argc, argv, "a:dfG:g:hi:O:o:S:s:T:t:V:v")) != -1) {
    switch(opt) {
      case 'a':
        sa_add(args->aerosols, 1);
        args->aerosols[args->naerosols] = RNATM_AEROSOL_ARGS_NULL;
        args->naerosols += 1;
        res = cstr_parse_list(optarg, ':', parse_aerosol_parameters, args);
        if(res == RES_OK) {
          res = check_aerosol_args(args->aerosols+args->naerosols-1);
        }
        break;
      case 'd':
        args->output_type = HTRDR_PLANETO_ARGS_OUTPUT_OCTREES;
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
        print_help(argv[0]);
        htrdr_planeto_args_release(args);
        args->quit = 1;
        goto exit;
      case 'i':
        res = htrdr_args_image_parse(&args->image, optarg);
        break;
      case 'O': args->octrees_storage = optarg; break;
      case 'o': args->output = optarg; break;
      case 'S':
        res = cstr_parse_list(optarg, ':', parse_source_parameters, args);
        break;
      case 's':
        res = htrdr_args_spectral_parse(&args->spectral_domain, optarg);
        break;
      case 'T':
        res = cstr_to_double(optarg, &args->optical_thickness);
        if(res != RES_OK && args->optical_thickness < 0) res = RES_BAD_ARG;
        break;
      case 't':
        res = cstr_to_uint(optarg, &args->nthreads);
        if(res != RES_OK && !args->nthreads) res = RES_BAD_ARG;
        break;
      case 'V':
        res = cstr_to_uint(optarg, &args->octree_definition_hint);
        if(res != RES_OK && !args->octree_definition_hint) res = RES_BAD_ARG;
        break;
      case 'v': args->verbose = 1; break;
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


  res = check_gas_args(&args->gas);
  if(res != RES_OK) {
    fprintf(stderr, "Missing gas definition -- option '-a'\n");
    goto error;
  }
  res = check_ground_args(&args->ground);
  if(res != RES_OK) {
    fprintf(stderr, "Missing ground definition -- option '-G'\n");
    goto error;
  }
  res = htrdr_planeto_source_args_check(&args->source);
  if(res != RES_OK) {
    fprintf(stderr, "Missing source definition -- option '-S'\n");
    goto error;
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
  size_t i;
  ASSERT(args);

  if(args->gas.smsh_filename) mem_rm(args->gas.smsh_filename);
  if(args->gas.sck_filename) mem_rm(args->gas.sck_filename);
  if(args->gas.temperatures_filename) mem_rm(args->gas.temperatures_filename);
  if(args->ground.smsh_filename) mem_rm(args->ground.smsh_filename);
  if(args->ground.props_filename) mem_rm(args->ground.props_filename);
  if(args->ground.mtllst_filename) mem_rm(args->ground.mtllst_filename);
  if(args->ground.name) mem_rm(args->ground.name);

  FOR_EACH(i, 0, args->naerosols) {
    struct rnatm_aerosol_args* aerosol = args->aerosols + i;
    if(aerosol->name) mem_rm(aerosol->name);
    if(aerosol->smsh_filename) mem_rm(aerosol->smsh_filename);
    if(aerosol->sars_filename) mem_rm(aerosol->sars_filename);
    if(aerosol->phase_fn_ids_filename) mem_rm(aerosol->phase_fn_ids_filename);
    if(aerosol->phase_fn_lst_filename) mem_rm(aerosol->phase_fn_lst_filename);
  }
  sa_release(args->aerosols);

  *args = HTRDR_PLANETO_ARGS_DEFAULT;
}

res_T
htrdr_planeto_args_check(const struct htrdr_planeto_args* args)
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

  /* Check the ground */
  res = check_ground_args(&args->ground);
  if(res != RES_OK) return res;

  /* Check the octree parameters */
  if(args->octree_definition_hint == 0
  || args->optical_thickness < 0)
    return RES_BAD_ARG;

  /* Check the spectral domain */
  res = check_spectral_args(&args->spectral_domain);
  if(res != RES_OK) return res;

  /* Check the source */
  res = htrdr_planeto_source_args_check(&args->source);
  if(res != RES_OK) return res;

  /* Check miscalleneous parameters */
  if(args->nthreads == 0
  || (unsigned)args->output_type >= HTRDR_PLANETO_ARGS_OUTPUT_TYPES_COUNT__)
    return RES_BAD_ARG;

  return RES_OK;
}

res_T
htrdr_planeto_source_args_check(const struct htrdr_planeto_source_args* args)
{
  if(!args) return RES_BAD_ARG;

  /* Invalid position */
  if(args->latitude <-90
  || args->latitude > 90
  || args->longitude <-180
  || args->longitude > 180
  || args->distance < 0)
    return RES_BAD_ARG;

  /* Miscellaneous parameters */
  if(args->radius < 0
  || args->temperature < 0)
    return RES_BAD_ARG;

  return RES_OK;
}
