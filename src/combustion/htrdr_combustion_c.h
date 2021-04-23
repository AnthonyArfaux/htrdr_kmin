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

#ifndef HTRDR_COMBUSTION_C_H
#define HTRDR_COMBUSTION_C_H

#include "core/htrdr_accum.h"
#include "core/htrdr_args.h"
#include "core/htrdr_buffer.h"

#include <rsys/ref_count.h>
#include <rsys/str.h>

/* Forward declarations */
struct atrstm;
struct htrdr;
struct htrdr_camera;
struct htrdr_combustion_laser;
struct htrdr_geometry;
struct htrdr_materials;
struct htrdr_rectangle;
struct ssf_phase;
struct ssp_rng;

struct combustion_pixel {
  struct htrdr_estimate radiance; /* In W/m^2/sr */
  struct htrdr_accum time; /* In microseconds */
};
#define COMBUSTION_PIXEL_NULL__ {                                              \
  HTRDR_ESTIMATE_NULL__, /* Radiance */                                        \
  HTRDR_ACCUM_NULL__, /* Time */                                               \
}
static const struct combustion_pixel COMBUSTION_PIXEL_NULL =
  COMBUSTION_PIXEL_NULL__;

struct htrdr_combustion {
  struct htrdr_geometry* geom; /* Combustion chamber geometry */
  struct htrdr_materials* mats; /* Materials of the combustion chamber */
  struct atrstm* medium; /* Semi transparent medium */

  struct htrdr_camera* camera; /* Pinhole camera */
  struct htrdr_combustion_laser* laser; /* Laser sheet */
  double wavelength; /* Wavelength of the laser in nanometer */

  struct ssf_phase** rdgfa_phase_functions; /* Per thread RDG-FA phase func */
  struct ssf_phase** hg_phase_functions; /* Per thread Henyey-Greenstein func */

  struct htrdr_buffer_layout buf_layout;
  struct htrdr_buffer* buf; /* NULL on non master processes */
  size_t spp; /* #samples per pixel */

  FILE* output; /* Output stream */
  struct str output_name; /* Name of the output stream */
  int dump_volumetric_acceleration_structure;

  ref_T ref;
  struct htrdr* htrdr;
};

extern LOCAL_SYM void
combustion_get_pixel_format
  (const struct htrdr_combustion* cmd,
   struct htrdr_pixel_format* fmt);

extern LOCAL_SYM res_T
combustion_draw_map
  (struct htrdr_combustion* cmd);

extern LOCAL_SYM res_T
combustion_compute_radiance_sw
  (struct htrdr_combustion* cmd,
   const size_t ithread,
   struct ssp_rng* rng,
   const double pos_in[3],
   const double dir_in[3],
   double* out_weigh); /* Shortwave radiance in W/m^2/sr */

#endif /* HTRDR_COMBUSTION_C_H */
