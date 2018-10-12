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

#define _POSIX_C_SOURCE 2 /* strtok_r support */

#include "htrdr_args.h"

#include <rsys/cstr.h>
#include <rsys/double3.h>

#include <getopt.h>
#include <string.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
print_help(const char* cmd)
{
  ASSERT(cmd);
  printf("Usage: %s -i INPUT [OPIONS]\n\n", cmd);
  printf(
"  -a FILENAME      path of gas optical properties file.\n");
  printf(
"  -c FILENAME      path of the HTCP cloud properties file.\n");
  printf(
"  -C <camera>      define the rendering point of view.\n");
  printf(
"  -d               dump octree data to OUTPUT wrt the VTK ASCII file format.\n");
  printf(
"  -D AZIMUTH,ELEVATION\n"
"                   sun direction in degrees. Following the right-handed\n"
"                   convention, the azimuthal rotation is counter-clockwise\n"
"                   around the Z axis, with 0 aligned on the X axis. The\n"
"                   elevation rotation starts from 0 up to 90 at zenith.\n");
  printf(
"  -f               overwrite the OUTPUT file if it already exists.\n");
  printf(
"  -g FILENAME      path of the OBJ geometry file.\n");
  printf(
"  -h               display this help and exit.\n");
  printf(
"  -i <image>       define the image to compute.\n");
  printf(
"  -r               infinitely repeat the clouds along the X and Y axes.\n");
  printf(
"  -m FILENAME      path of the Mie data file.\n");
  printf(
"  -o OUTPUT        file where data are written. If not defined, data are\n"
"                   written to standard output.\n");
  printf(
"  -t THREADS       hint on the number of threads to use. By default use as\n"
"                   many threads as CPU cores.\n");
  printf(
"  -T THRESHOLD     optical thickness used as threshold during the octree\n"
"                   building. By default its value is `%g'.\n",
    HTRDR_ARGS_DEFAULT.optical_thickness);
  printf(
"  -v               make the program more verbose.\n");
  printf("\n");
  printf(
"%s (C) 2018 Université Paul Sabatier, |Meso|Star>. This is free software\n"
"released under the GNU GPL license, version 3 or later. You are free to change\n"
"or redistribute it under certain conditions <http://gnu.org/licenses/gpl.html>.\n",
    cmd);
}

static INLINE res_T
parse_doubleX(const char* str, double* val, const size_t sz)
{
  size_t len;
  res_T res = RES_OK;
  ASSERT(str && val);
  res = cstr_to_list_double(str, ',', val, &len, sz);
  if(res == RES_OK && len != sz) res = RES_BAD_ARG;
  return res;
}

static INLINE res_T
parse_definition(const char* str, unsigned val[2])
{
  size_t len;
  res_T res = RES_OK;
  ASSERT(str && val);
  res = cstr_to_list_uint(str, 'x', val, &len, 2);
  if(res != RES_OK) return res;
  if(len != 2) return RES_BAD_ARG;
  if(val[0] > 16384 || val[1] > 16384) return RES_BAD_ARG;
  return RES_OK;
}

static res_T
parse_fov(const char* str, double* out_fov)
{
  double fov;
  res_T res = RES_OK;
  ASSERT(str && out_fov);

  res = cstr_to_double(str, &fov);
  if(res != RES_OK) {
    fprintf(stderr, "Invalid field of view `%s'.\n", str);
    return RES_BAD_ARG;
  }
  if(fov < 30 || fov > 120) {
    fprintf(stderr, "The field of view %g is not in [30, 120].\n", fov);
    return RES_BAD_ARG;
  }
  *out_fov = fov;
  return RES_OK;
}

static res_T
parse_image_parameter(struct htrdr_args* args, const char* str)
{
  char buf[128];
  char* key;
  char* val;
  char* ctx;
  res_T res = RES_OK;
  ASSERT(str && args);

  if(strlen(str) >= sizeof(buf) -1/*NULL char*/) {
    fprintf(stderr,
      "Could not duplicate the image option string `%s'.\n", str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  key = strtok_r(buf, "=", &ctx);
  val = strtok_r(NULL, "", &ctx);

  if(!val) {
    fprintf(stderr, "Missing a value to the image option `%s'.\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  #define PARSE(Name, Func)                                                    \
    res = Func;                                                                \
    if(res != RES_OK) {                                                        \
      fprintf(stderr, "Invalid image "Name" `%s'.\n", val);                    \
      goto error;                                                              \
    } (void)0
  if(!strcmp(key, "def")) {
    PARSE("definition", parse_definition(val, args->image.definition));
  } else if(!strcmp(key, "spp")) {
    PARSE("#samples per pixel", cstr_to_uint(val, &args->image.spp));
  } else {
    fprintf(stderr, "Invalid image parameter `%s'.\n", key);
    res = RES_BAD_ARG;
    goto error;
  }
  #undef PARSE

  if(!args->image.definition[0] || !args->image.definition[1]) {
    fprintf(stderr, "The image definition cannot be null.n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(!args->image.spp) {
    fprintf(stderr, "The number of samples per pixel cannot be null.\n");
    res = RES_BAD_ARG;
    goto error;
  }

exit:
  return res;
error:
  goto exit;
}

static res_T
parse_camera_parameter(struct htrdr_args* args, const char* str)
{
  char buf[128];
  char* key;
  char* val;
  char* ctx;
  res_T res = RES_OK;

  if(strlen(str) >= sizeof(buf) -1/*NULL char*/) {
    fprintf(stderr,
      "Could not duplicate the rectangle option string `%s'.\n", str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  key = strtok_r(buf, "=", &ctx);
  val = strtok_r(NULL, "", &ctx);

  if(!val) {
    fprintf(stderr, "Missing value to the camera option `%s'.\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  #define PARSE(Name, Func) {                                                  \
    if(RES_OK != (res = Func)) {                                               \
      fprintf(stderr, "Invalid camera "Name" `%s'.\n", val);                   \
      goto error;                                                              \
    }                                                                          \
  } (void)0
  if(!strcmp(key, "pos")) {
    PARSE("position", parse_doubleX(val, args->camera.pos, 3));
  } else if(!strcmp(key, "tgt")) {
    PARSE("target", parse_doubleX(val, args->camera.tgt, 3));
  } else if(!strcmp(key, "up")) {
    PARSE("up vector", parse_doubleX(val, args->camera.up, 3));
  } else if(!strcmp(key, "fov")) {
    PARSE("field-of-view", parse_fov(val, &args->camera.fov_y));
  } else {
    fprintf(stderr, "Invalid camera parameter `%s'.\n", key);
    res = RES_BAD_ARG;
    goto error;
  }
  #undef PARSE
exit:
  return res;
error:
  goto exit;
}

static res_T
parse_multiple_parameters
  (struct htrdr_args* args,
   const char* str,
   res_T (*parse_parameter)(struct htrdr_args* args, const char* str))
{
  char buf[512];
  char* tk;
  char* ctx;
  res_T res = RES_OK;
  ASSERT(args && str);

  if(strlen(str) >= sizeof(buf) - 1/*NULL char*/) {
    fprintf(stderr, "Could not duplicate the option string `%s'.\n", str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  tk = strtok_r(buf, ":", &ctx);
  do {
    res = parse_parameter(args, tk);
    if(res != RES_OK) goto error;
    tk = strtok_r(NULL, ":", &ctx);
  } while(tk);

exit:
  return res;
error:
  goto exit;
}

static res_T
parse_sun_dir(struct htrdr_args* args, const char* str)
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
htrdr_args_init(struct htrdr_args* args, int argc, char** argv)
{
  int opt;
  res_T res = RES_OK;
  ASSERT(args && argc && argv);

  *args = HTRDR_ARGS_DEFAULT;

  while((opt = getopt(argc, argv, "a:C:c:D:dfg:hi:m:o:rT:t:v")) != -1) {
    switch(opt) {
      case 'a': args->filename_gas = optarg; break;
      case 'C':
        res = parse_multiple_parameters
          (args, optarg, parse_camera_parameter);
        break;
      case 'c': args->filename_les = optarg; break;
      case 'D': res = parse_sun_dir(args, optarg); break;
      case 'd': args->dump_vtk = 1; break;
      case 'f': args->force_overwriting = 1; break;
      case 'g': args->filename_obj = optarg; break;
      case 'h':
        print_help(argv[0]);
        htrdr_args_release(args);
        args->quit = 1;
        goto exit;
      case 'i':
        res = parse_multiple_parameters
          (args, optarg, parse_image_parameter);
        break;
      case 'm': args->filename_mie = optarg; break;
      case 'o': args->output = optarg; break;
      case 'r': args->repeat_clouds = 1; break;
      case 'T':
        res = cstr_to_double(optarg, &args->optical_thickness);
        if(res == RES_OK && args->optical_thickness < 0) res = RES_BAD_ARG;
        break;
      case 't': /* Submit an hint on the number of threads to use */
        res = cstr_to_uint(optarg, &args->nthreads);
        if(res == RES_OK && !args->nthreads) res = RES_BAD_ARG;
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
  if(!args->filename_gas) {
    fprintf(stderr,
      "Missing the path of the gas optical properties file -- option '-a'\n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(!args->filename_les) {
    fprintf(stderr,
      "Missing the path of the cloud properties file -- option '-c'\n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(!args->filename_mie) {
    fprintf(stderr,
      "Missing the path toward the file of the Mie's data -- option '-m'\n");
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

