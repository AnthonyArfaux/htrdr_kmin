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
"  -d               dump octree data to OUTPUT wrt the VTK ASCII file format.\n");
  printf(
"  -f               overwrite the OUTPUT file if it already exists.\n");
  printf(
"  -h               display this help and exit.\n");
  printf(
"  -i INPUT         path of the input HTCP file.\n");
  printf(
"  -I <image>       define the image to compute.\n");
  printf(
"  -o OUTPUT        file where data are written. If not defined, data are\n"
"                   written to standard output.\n");
  printf(
"  -r <rectangle>   define the integration plane.\n");
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
parse_rectangle_parameter(struct htrdr_args* args, const char* str)
{
  char buf[128];
  char* key;
  char* val;
  char* ctx;
  res_T res;
  ASSERT(str && args);

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
    fprintf(stderr, "Missing a value to the rectangle option `%s'.\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  #define PARSE(Name, Func)                                                    \
    res = Func;                                                                \
    if(res != RES_OK) {                                                        \
      fprintf(stderr, "Invalid rectangle "Name" `%s'.\n", val);                \
      goto error;                                                              \
    } (void)0
  if(!strcmp(key, "pos")) {
    PARSE("position", parse_doubleX(val, args->rectangle.pos, 3));
  } else if(!strcmp(key, "tgt")) {
    PARSE("target", parse_doubleX(val, args->rectangle.tgt, 3));
  } else if(!strcmp(key, "up")) {
    PARSE("up vector", parse_doubleX(val, args->rectangle.up, 3));
  } else if(!strcmp(key, "sz")) {
    PARSE("size", parse_doubleX(val, args->rectangle.sz, 2));
  } else {
    fprintf(stderr, "Invalid rectangle parameter `%s'.\n", key);
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
parse_image_parameter(struct htrdr_args* args, const char* str)
{
  char buf[128];
  char* key;
  char* val;
  char* ctx;
  res_T res;
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
  } else if(!strcmp(key, "ssp")) {
    PARSE("#samples per pixel", cstr_to_uint(val, &args->image.spp));
  } else {
    fprintf(stderr, "Invalid image parameter `%s'.\n", key);
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

  while((opt = getopt(argc, argv, "dfhI:i:o:r:v")) != -1) {
    switch(opt) {
      case 'd': args->dump_vtk = 1; break;
      case 'f': args->force_overwriting = 1; break;
      case 'h':
        print_help(argv[0]);
        htrdr_args_release(args);
        args->quit = 1;
        goto exit;
      case 'I':
        res = parse_multiple_parameters
          (args, optarg, parse_image_parameter);
        break;
      case 'i': args->input = optarg; break;
      case 'o': args->output = optarg; break;
      case 'r':
        res = parse_multiple_parameters
          (args, optarg, parse_rectangle_parameter);
        break;
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

