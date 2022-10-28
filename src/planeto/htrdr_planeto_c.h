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

#ifndef HTRDR_PLANETO_C_H
#define HTRDR_PLANETO_C_H

#include "planeto/htrdr_planeto_args.h"

#include "core/htrdr_accum.h"
#include "core/htrdr_args.h"
#include "core/htrdr_buffer.h"

#include <rsys/ref_count.h>
#include <rsys/str.h>

/* Forward declarations */
struct htrdr;
struct htrdr_pixel_format;
struct rnatm;
struct rngrd;

struct planeto_pixel_xwave {
  struct htrdr_estimate radiance; /* In W/m²/sr */
  struct htrdr_estimate radiance_temperature; /* In W/m²/sr */
  struct htrdr_accum time; /* In µs */
};
#define PLANETO_PIXEL_XWAVE_NULL__ {                                           \
  HTRDR_ESTIMATE_NULL__,                                                       \
  HTRDR_ESTIMATE_NULL__,                                                       \
  HTRDR_ACCUM_NULL__                                                           \
}
static const struct planeto_pixel_xwave PLANETO_PIXEL_XWAVE_NULL =
  PLANETO_PIXEL_XWAVE_NULL__;

struct planeto_pixel_image {
  struct htrdr_estimate X; /* In W/m²/sr */
  struct htrdr_estimate Y; /* In W/m²/sr */
  struct htrdr_estimate Z; /* In W/m²/sr */
  struct htrdr_accum time; /* In µs */
};
#define PLANETO_PIXEL_IMAGE_NULL__ {                                           \
  HTRDR_ESTIMATE_NULL__,                                                       \
  HTRDR_ESTIMATE_NULL__,                                                       \
  HTRDR_ESTIMATE_NULL__,                                                       \
  HTRDR_ACCUM_NULL__                                                           \
}

struct planeto_pixel_xwave;

struct htrdr_planeto {
  struct rnatm* atmosphere;
  struct rngrd* ground;
  struct htrdr_planeto_source* source;

  struct htrdr_args_spectral spectral_domain;

  FILE* octrees_storage;

  FILE* output;
  struct str output_name;
  enum htrdr_planeto_args_output_type output_type;

  struct htrdr_buffer_layout buf_layout;
  struct htrdr_buffer* buf; /* NULL on non master processes */
  size_t spp; /* Samples per pixel */

  ref_T ref;
  struct htrdr* htrdr;
};

extern LOCAL_SYM void
planeto_get_pixel_format
  (const struct htrdr_planeto* cmd,
   struct htrdr_pixel_format* fmt);

#endif /* HTRDR_PLANETO_C_H */
