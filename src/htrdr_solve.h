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
#define HTRDR_ACCUM_NULL__ {0,0,0,0}
static const struct htrdr_accum HTRDR_ACCUM_NULL = HTRDR_ACCUM_NULL__;

/* Forward declarations */
struct htrdr;
struct htrdr_camera;
struct ssp_rng;

extern LOCAL_SYM double
htrdr_compute_radiance_sw
  (struct htrdr* htrdr,
   const size_t ithread,
   struct ssp_rng* rng,
   const double pos[3],
   const double dir[3],
   const size_t iband, /* Index of the spectral band */
   const size_t iquad); /* Index of the quadrature point into the band */

extern LOCAL_SYM res_T
htrdr_draw_radiance_sw
  (struct htrdr* htrdr,
   const struct htrdr_camera* cam,
   const size_t spp, /* #samples per pixel, i.e. #realisations */
   struct htrdr_buffer* buf); /* Buffer of struct htrdr_accum[3] */

#endif /* HTRDR_SOLVE_H */
