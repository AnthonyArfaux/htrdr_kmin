/* Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
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

#include "htrdr.h"

#include "htrdr_args.h"
#include "htrdr_version.h"

#include <rsys/cstr.h>
#include <rsys/double3.h>

#include <getopt.h>
#include <string.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
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
  if(fov <= 0 || fov >= 180) {
    fprintf(stderr, "The field of view %g is not in [30, 120].\n", fov);
    return RES_BAD_ARG;
  }
  *out_fov = fov;
  return RES_OK;
}

static res_T
parse_image_parameter(void* args, const char* str)
{
  char buf[128];
  struct htrdr_args_image* img = args;
  char* key;
  char* val;
  char* ctx;
  res_T res = RES_OK;
  ASSERT(str && img);

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
    PARSE("definition", parse_definition(val, img->definition));
  } else if(!strcmp(key, "spp")) {
    PARSE("#samples per pixel", cstr_to_uint(val, &img->spp));
  } else {
    fprintf(stderr, "Invalid image parameter `%s'.\n", key);
    res = RES_BAD_ARG;
    goto error;
  }
  #undef PARSE

  if(!img->definition[0] || !img->definition[1]) {
    fprintf(stderr, "The image definition cannot be null.n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(!img->spp) {
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
parse_camera_parameter(void* args, const char* str)
{
  char buf[128];
  struct htrdr_args_camera* cam = args;
  char* key;
  char* val;
  char* ctx;
  res_T res = RES_OK;
  ASSERT(cam && str);

  if(strlen(str) >= sizeof(buf) -1/*NULL char*/) {
    fprintf(stderr,
      "Could not duplicate the camera option string `%s'.\n", str);
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
    PARSE("position", parse_doubleX(val, cam->position, 3));
  } else if(!strcmp(key, "tgt")) {
    PARSE("target", parse_doubleX(val, cam->target, 3));
  } else if(!strcmp(key, "up")) {
    PARSE("up vector", parse_doubleX(val, cam->up, 3));
  } else if(!strcmp(key, "fov")) {
    PARSE("field-of-view", parse_fov(val, &cam->fov_y));
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
parse_rectangle_parameter(void* args, const char* str)
{
  char buf[128];
  struct htrdr_args_rectangle* rect = args;
  char* key;
  char* val;
  char* ctx;
  res_T res = RES_OK;
  ASSERT(rect && str);

  if(strlen(str) >= sizeof(buf) -1/*NULL char*/) {
    fprintf(stderr,
      "Could not duplicate the rectangle option string `%s'.\n", str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  /* pos=0,0,10.1; key <- pos, val <- 0,0,10 */
  key = strtok_r(buf, "=", &ctx);
  val = strtok_r(NULL, "", &ctx);

  if(!val) {
    fprintf(stderr, "Missing value to the rectangle option `%s'.\n", key);
    res = RES_BAD_ARG;
    goto error;
  }

  #define PARSE(Name, Func) {                                                  \
    if(RES_OK != (res = Func)) {                                               \
      fprintf(stderr, "Invalid rectangle "Name" `%s'.\n", val);                \
      goto error;                                                              \
    }                                                                          \
  } (void)0
  if(!strcmp(key, "pos")) {
    PARSE("position", parse_doubleX(val, rect->position, 3));
  } else if(!strcmp(key, "tgt")) {
    PARSE("target", parse_doubleX(val, rect->target, 3));
  } else if(!strcmp(key, "up")) {
    PARSE("up vector", parse_doubleX(val, rect->up, 3));
  } else if(!strcmp(key, "sz")) {
    PARSE("size", parse_doubleX(val, rect->size, 2));
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
parse_spectral_range(const char* str, double wlen_range[2])
{
  double range[2];
  size_t len;
  res_T res = RES_OK;
  ASSERT(wlen_range && str);

  res = cstr_to_list_double(str, ',', range, &len, 2);
  if(res == RES_OK && len != 2) res = RES_BAD_ARG;
  if(res == RES_OK && range[0] > range[1]) res = RES_BAD_ARG;
  if(res != RES_OK) {
    fprintf(stderr, "Invalid spectral range `%s'.\n", str);
    goto error;
  }
  wlen_range[0] = range[0];
  wlen_range[1] = range[1];

exit:
  return res;
error:
  goto exit;
}

static res_T
parse_spectral_parameter(void* ptr, const char* str)
{
  char buf[128];
  struct htrdr_args_spectral* args = ptr;
  char* key;
  char* val;
  char* ctx;
  res_T res = RES_OK;
  ASSERT(args && str);

  if(strlen(str) >= sizeof(buf) -1/*NULL char*/) {
    fprintf(stderr,
      "Could not duplicate the spectral option string `%s'.\n", str);
    res = RES_MEM_ERR;
    goto error;
  }
  strncpy(buf, str, sizeof(buf));

  key = strtok_r(buf, "=", &ctx);
  val = strtok_r(NULL, "",  &ctx);

  if(!strcmp(key, "cie_xyz")) {
    args->spectral_type = HTRDR_SPECTRAL_SW_CIE_XYZ;
    args->wlen_range[0] = HTRDR_CIE_XYZ_RANGE_DEFAULT[0];
    args->wlen_range[1] = HTRDR_CIE_XYZ_RANGE_DEFAULT[1];
  } else {
    if(!val) {
      fprintf(stderr, "Missing value to the spectral option `%s'.\n", key);
      res = RES_BAD_ARG;
      goto error;
    }
    if(!strcmp(key, "sw")) {
      args->spectral_type = HTRDR_SPECTRAL_SW;
      res = parse_spectral_range(val, args->wlen_range);
      if(res != RES_OK) goto error;
    } else if(!strcmp(key, "lw")) {
      args->spectral_type = HTRDR_SPECTRAL_LW;
      res = parse_spectral_range(val, args->wlen_range);
      if(res != RES_OK) goto error;
    } else if(!strcmp(key, "Tref")) {
      res = cstr_to_double(val, &args->ref_temperature);
      if(res == RES_OK && args->ref_temperature < 0) res = RES_BAD_ARG;
      if(res != RES_OK) {
        fprintf(stderr, "Invalid reference temperature Tref=%s.\n", val);
        goto error;
      }
    } else {
      fprintf(stderr, "Invalid spectral parameter `%s'.\n", key);
      res = RES_BAD_ARG;
      goto error;
    }
  }

exit:
  return res;
error:
  goto exit;
}

static res_T
parse_multiple_parameters
  (void* args,
   const char* str,
   res_T (*parse_parameter)(void* args, const char* str))
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
 * Exported functions
 ******************************************************************************/
res_T
htrdr_args_camera_parse(struct htrdr_args_camera* cam, const char* str)
{
  if(!cam || !str) return RES_BAD_ARG;
  return parse_multiple_parameters(cam, str, parse_camera_parameter);
}

res_T
htrdr_args_rectangle_parse(struct htrdr_args_rectangle* rect, const char* str)
{
  if(!rect || !str) return RES_BAD_ARG;
  return parse_multiple_parameters(rect, str, parse_rectangle_parameter);
}
 
res_T
htrdr_args_image_parse(struct htrdr_args_image* img, const char* str)
{
  if(!img || !str) return RES_BAD_ARG;
  return parse_multiple_parameters(img, str, parse_image_parameter);
}

res_T
htrdr_args_spectral_parse(struct htrdr_args_spectral* spectral, const char* str)
{
  if(!spectral || !str) return RES_BAD_ARG;
  return parse_multiple_parameters(spectral, str, parse_spectral_parameter);
}

