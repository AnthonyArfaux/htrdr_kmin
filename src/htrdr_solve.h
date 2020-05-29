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

#ifndef HTRDR_SOLVE_H
#define HTRDR_SOLVE_H

#include <rsys/rsys.h>

/* Monte carlo accumulator */
struct htrdr_accum {
  double sum_weights; /* Sum of Monte-Carlo weights */
  double sum_weights_sqr; /* Sum of Monte-Carlo square weights */
  size_t nweights; /* #accumlated weights */
  size_t nfailures; /* #failures */
};
static const struct htrdr_accum HTRDR_ACCUM_NULL;

/* Monte carlo estimate */
struct htrdr_estimate {
  double E; /* Expected value */
  double SE; /* Standard error */
};
static const struct htrdr_estimate HTRDR_ESTIMATE_NULL;

struct htrdr_pixel_sw {
  struct htrdr_estimate X; /* In W/m^2/sr */
  struct htrdr_estimate Y; /* In W/m^2/sr */
  struct htrdr_estimate Z; /* In W/m^2/sr */
  struct htrdr_accum time; /* In microseconds */
};
static const struct htrdr_pixel_sw HTRDR_PIXEL_SW_NULL;

struct htrdr_pixel_lw {
  struct htrdr_estimate radiance; /* In W/m^2/sr */
  struct htrdr_estimate radiance_temperature; /* In K */
  struct htrdr_accum time; /* In microseconds */
};
static const struct htrdr_pixel_lw HTRDR_PIXEL_LW_NULL;

/* Forward declarations */
struct htrdr;
struct htrdr_camera;
struct s3d_hit;
struct ssp_rng;

/* Return the shortwave radiance in W/m^2/sr/m */
extern LOCAL_SYM double
htrdr_compute_radiance_sw
  (struct htrdr* htrdr,
   const size_t ithread,
   struct ssp_rng* rng,
   const double pos[3],
   const double dir[3],
   const double wlen, /* In nanometer */
   const size_t iband, /* Index of the spectral band */
   const size_t iquad); /* Index of the quadrature point into the band */

/* Return the longwave radiance in W/m^2/sr/m */
extern LOCAL_SYM double
htrdr_compute_radiance_lw
  (struct htrdr* htrdr,
   const size_t ithread,
   struct ssp_rng* rng,
   const double pos[3],
   const double dir[3],
   const double wlen, /* In nanometer */
   const size_t iband, /* Index of the spectral band */
   const size_t iquad); /* Index of the quadrature point into the band */

extern LOCAL_SYM res_T
htrdr_draw_radiance
  (struct htrdr* htrdr,
   const struct htrdr_camera* cam,
   const size_t width, /* Image width */
   const size_t height, /* Image height */
   const size_t spp, /* #samples per pixel, i.e. #realisations */
   /* Buffer of struct htrdr_accum[4]. May be NULL on non master processes */
   struct htrdr_buffer* buf);

extern LOCAL_SYM int
htrdr_ground_filter
  (const struct s3d_hit* hit,
   const float ray_dorg[3],
   const float ray_dir[3],
   void* ray_data,
   void* filter_data);

static FINLINE void
htrdr_accum_get_estimation
  (const struct htrdr_accum* acc,
   struct htrdr_estimate* estimate)
{
  ASSERT(acc && estimate);

  if(!acc->nweights) {
    estimate->E = 0;
    estimate->SE = 0;
  } else {
    const double N = (double)acc->nweights;
    double E, V, SE;
    E = acc->sum_weights / N;
    V = MMAX(acc->sum_weights_sqr / N - E*E, 0);
    SE = sqrt(V/N);

    estimate->E = E;
    estimate->SE = SE;
  }
}

#endif /* HTRDR_SOLVE_H */
