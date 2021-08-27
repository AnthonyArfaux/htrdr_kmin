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

#ifndef HTRDR_ATMOSPHERE_C_H
#define HTRDR_ATMOSPHERE_C_H

#include "atmosphere/htrdr_atmosphere_args.h"

#include "core/htrdr_accum.h"
#include "core/htrdr_buffer.h"
#include "core/htrdr_spectral.h"

#include <rsys/ref_count.h>
#include <rsys/rsys.h>
#include <rsys/str.h>

/* Define the radiance component */
enum atmosphere_radiance_cpnt_flag {
  ATMOSPHERE_RADIANCE_DIRECT = BIT(0),
  ATMOSPHERE_RADIANCE_DIFFUSE = BIT(1),
  ATMOSPHERE_RADIANCE_ALL =
    ATMOSPHERE_RADIANCE_DIRECT
  | ATMOSPHERE_RADIANCE_DIFFUSE
};

struct atmosphere_pixel_xwave {
  struct htrdr_estimate radiance; /* In W/m^2/sr */
  struct htrdr_estimate radiance_temperature; /* In K */
  struct htrdr_accum time; /* In microseconds */
};
#define ATMOSPHERE_PIXEL_XWAVE_NULL__ {                                        \
  HTRDR_ESTIMATE_NULL__, /* Radiance */                                        \
  HTRDR_ESTIMATE_NULL__, /* Radiance temperature */                            \
  HTRDR_ACCUM_NULL__ /* Time */                                                \
}
static const struct atmosphere_pixel_xwave ATMOSPHERE_PIXEL_XWAVE_NULL =
  ATMOSPHERE_PIXEL_XWAVE_NULL__;

struct atmosphere_pixel_flux {
  struct htrdr_accum flux;
  struct htrdr_accum time;
};
#define ATMOSPHERE_PIXEL_FLUX_NULL__ {                                         \
  HTRDR_ACCUM_NULL__,                                                          \
  HTRDR_ACCUM_NULL__                                                           \
}
static const struct atmosphere_pixel_flux ATMOSPHERE_PIXEL_FLUX_NULL =
  ATMOSPHERE_PIXEL_FLUX_NULL__;

struct atmosphere_pixel_image {
  struct htrdr_estimate X; /* In W/m^2/sr */
  struct htrdr_estimate Y; /* In W/m^2/sr */
  struct htrdr_estimate Z; /* In W/m^2/sr */
  struct htrdr_accum time; /* In microseconds */
};
#define ATMOSPHERE_PIXEL_IMAGE_NULL__ {                                        \
  HTRDR_ESTIMATE_NULL__, /* X */                                               \
  HTRDR_ESTIMATE_NULL__, /* Y */                                               \
  HTRDR_ESTIMATE_NULL__, /* Z */                                               \
  HTRDR_ACCUM_NULL__ /* Time */                                                \
}
static const struct atmosphere_pixel_image ATMOSPHERE_PIXEL_IMAGE_NULL =
  ATMOSPHERE_PIXEL_IMAGE_NULL__;

/* Forward declarations */
struct htsky;
struct htrdr;
struct htrdr_atmosphere_args;
struct htrdr_buffer;
struct htrdr_cie_xyz;
struct htrdr_materials;
struct htrdr_ran_wlen;
struct ssp_rng;

struct htrdr_atmosphere {
  struct htrdr_atmosphere_ground* ground;
  struct htrdr_atmosphere_sun* sun;
  struct htrdr_materials* mats;
  struct htrdr_cie_xyz* cie;
  struct htrdr_ran_wlen* ran_wlen;

  struct scam* camera; /* Camera */
  struct htrdr_rectangle* flux_map; /* Flux map */

  struct htrdr_buffer_layout buf_layout;
  struct htrdr_buffer* buf; /* NULL on non master processes */

  struct htsky* sky;
  const char* sky_mtl_name;
  enum htrdr_spectral_type spectral_type;
  double wlen_range_m[2]; /* Integration range in *meters* */
  double ref_temperature; /* Reference temperature in Kelvin */

  size_t spp; /* #samples per pixel */
  size_t width; /* Image width */
  size_t height; /* Image height */

  FILE* output;
  struct str output_name;

  unsigned grid_max_definition[3]; /* Max definition of the acceleration grids */
  unsigned nthreads; /* #threads of the process */
  enum htrdr_atmosphere_args_output_type output_type;
  int verbose; /* Verbosity level */

  ref_T ref;
  struct htrdr* htrdr;
};

extern LOCAL_SYM void
atmosphere_get_pixel_format
  (const struct htrdr_atmosphere* cmd,
   struct htrdr_pixel_format* fmt);

extern LOCAL_SYM res_T
atmosphere_draw_map
  (struct htrdr_atmosphere* cmd);

/* Return the shortwave radiance in W/m^2/sr/m */
extern LOCAL_SYM double
atmosphere_compute_radiance_sw
  (struct htrdr_atmosphere* cmd,
   const size_t ithread,
   struct ssp_rng* rng,
   const int cpnt_mask, /* Combination of enum atmosphere_radiance_cpnt_flag */
   const double pos_in[3],
   const double dir_in[3],
   const double wlen, /* In nanometer */
   const size_t iband,
   const size_t iquad);

/* Return the longwave radiance in W/m^2/sr/m */
extern LOCAL_SYM double
atmosphere_compute_radiance_lw
  (struct htrdr_atmosphere* cmd,
   const size_t ithread,
   struct ssp_rng* rng,
   const double pos_in[3],
   const double dir_in[3],
   const double wlen, /* In nanometer */
   const size_t iband,
   const size_t iquad);

#endif /* HTRDR_ATMOSPHERE_C_H */

