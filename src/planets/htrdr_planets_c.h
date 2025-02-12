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

#ifndef HTRDR_PLANETS_C_H
#define HTRDR_PLANETS_C_H

#include "planets/htrdr_planets_args.h"

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

struct planets_pixel_xwave {
  struct htrdr_accum radiance; /* In W/m²/sr */
  struct htrdr_accum time; /* In µs */
  struct htrdr_estimate radiance_temperature; /* In W/m²/sr */
};
#define PLANETS_PIXEL_XWAVE_NULL__ {                                           \
  HTRDR_ACCUM_NULL__,                                                          \
  HTRDR_ACCUM_NULL__,                                                          \
  HTRDR_ESTIMATE_NULL__                                                        \
}
static const struct planets_pixel_xwave PLANETS_PIXEL_XWAVE_NULL =
  PLANETS_PIXEL_XWAVE_NULL__;

struct planets_pixel_image {
  struct htrdr_estimate X; /* In W/m²/sr */
  struct htrdr_estimate Y; /* In W/m²/sr */
  struct htrdr_estimate Z; /* In W/m²/sr */
  struct htrdr_accum time; /* In µs */
};
#define PLANETS_PIXEL_IMAGE_NULL__ {                                           \
  HTRDR_ESTIMATE_NULL__,                                                       \
  HTRDR_ESTIMATE_NULL__,                                                       \
  HTRDR_ESTIMATE_NULL__,                                                       \
  HTRDR_ACCUM_NULL__                                                           \
}

struct planets_compute_radiance_args {
  struct ssp_rng* rng;
  size_t ithread; /* Index of the thread executing the function */

  double path_org[3]; /* Origin of the path to trace */
  double path_dir[3]; /* Initial direction of the path to trace */

  double wlen; /* In nm */
  size_t iband; /* Spectral band index */
  size_t iquad; /* Quadrature point */
};
#define PLANETS_COMPUTE_RADIANCE_ARGS_NULL__ {NULL, 0, {0,0,0}, {0,0,0}, 0, 0, 0}
static const struct planets_compute_radiance_args
PLANETS_COMPUTE_RADIANCE_ARGS_NULL = PLANETS_COMPUTE_RADIANCE_ARGS_NULL__;

struct htrdr_planets {
  struct rnatm* atmosphere;
  struct rngrd* ground;
  struct htrdr_planets_source* source;

  struct htrdr_planets_spectral_args spectral_domain;
  struct htrdr_ran_wlen_cie_xyz* cie; /* HTRDR_SPECTRAL_SW_CIE_XYZ */
  struct htrdr_ran_wlen_planck* planck; /* HTRDR_SPECTRAL_<LW|SW> */
  struct htrdr_ran_wlen_discrete* discrete; /* HTRDR_SPECTRAL_SW */

  FILE* octrees_storage;

  FILE* output;
  struct str output_name;
  enum htrdr_planets_args_output_type output_type;

  struct scam* camera;

  struct htrdr_buffer_layout buf_layout;
  struct htrdr_buffer* buf; /* NULL on non master processes */
  size_t spp; /* Samples per pixel */

  ref_T ref;
  struct htrdr* htrdr;
};

extern LOCAL_SYM res_T
planets_draw_map
  (struct htrdr_planets* cmd);

extern LOCAL_SYM void
planets_get_pixel_format
  (const struct htrdr_planets* cmd,
   struct htrdr_pixel_format* fmt);

/* Return the radiance in W/m²/sr/m */
extern LOCAL_SYM double
planets_compute_radiance
  (struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args);

#endif /* HTRDR_PLANETS_C_H */
