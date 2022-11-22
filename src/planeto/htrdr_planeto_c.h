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
struct htrdr_ran_wlen_cie_xyz;
struct htrdr_ran_wlen_planck;
struct rnatm;
struct rngrd;
struct scam;

struct planeto_pixel_xwave {
  struct htrdr_accum radiance; /* In W/m²/sr */
  struct htrdr_accum time; /* In µs */
  struct htrdr_estimate radiance_temperature; /* In W/m²/sr */
};
#define PLANETO_PIXEL_XWAVE_NULL__ {                                           \
  HTRDR_ACCUM_NULL__,                                                          \
  HTRDR_ACCUM_NULL__,                                                          \
  HTRDR_ESTIMATE_NULL__                                                        \
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

struct planeto_compute_radiance_args {
  struct ssp_rng* rng;
  size_t ithread; /* Index of the thread executing the function */

  double path_org[3]; /* Origin of the path to trace */
  double path_dir[3]; /* Initial direction of the path to trace */

  double wlen; /* In nm */
  size_t iband; /* Spectral band index */
  size_t iquad; /* Quadrature point */
};
#define PLANETO_COMPUTE_RADIANCE_ARGS_NULL__ {NULL, 0, {0,0,0}, {0,0,0}, 0, 0, 0}
static const struct planeto_compute_radiance_args
PLANETO_COMPUTE_RADIANCE_ARGS_NULL = PLANETO_COMPUTE_RADIANCE_ARGS_NULL__;

struct htrdr_planeto {
  struct rnatm* atmosphere;
  struct rngrd* ground;
  struct htrdr_planeto_source* source;

  struct htrdr_args_spectral spectral_domain;
  struct htrdr_ran_wlen_cie_xyz* cie; /* HTRDR_SPECTRAL_SW_CIE_XYZ */
  struct htrdr_ran_wlen_planck* planck; /* HTRDR_SPECTRAL_<LW|SW> */

  FILE* octrees_storage;

  FILE* output;
  struct str output_name;
  enum htrdr_planeto_args_output_type output_type;

  struct scam* camera;

  struct htrdr_buffer_layout buf_layout;
  struct htrdr_buffer* buf; /* NULL on non master processes */
  size_t spp; /* Samples per pixel */

  ref_T ref;
  struct htrdr* htrdr;
};

extern LOCAL_SYM res_T
planeto_draw_map
  (struct htrdr_planeto* cmd);

extern LOCAL_SYM void
planeto_get_pixel_format
  (const struct htrdr_planeto* cmd,
   struct htrdr_pixel_format* fmt);

/* Return the radiance in W/m²/sr/m */
extern LOCAL_SYM double
planeto_compute_radiance
  (struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args);

#endif /* HTRDR_PLANETO_C_H */
