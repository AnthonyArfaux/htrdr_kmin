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

#include "planeto/htrdr_planeto_c.h"

#include "core/htrdr.h"
#include "core/htrdr_accum.h"
#include "core/htrdr_buffer.h"
#include "core/htrdr_draw_map.h"
#include "core/htrdr_log.h"
#include "core/htrdr_ran_wlen_cie_xyz.h"

#include <rad-net/rnatm.h>
#include <star/scam.h>
#include <star/ssp.h>

#include <rsys/clock_time.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
draw_pixel_xwave
  (struct htrdr* htrdr,
   const struct htrdr_draw_pixel_args* args,
   void* data)
{
  /* TODO */
  (void)htrdr, (void)args, (void)data;
}

static void
draw_pixel_image
  (struct htrdr* htrdr,
   const struct htrdr_draw_pixel_args* args,
   void* data)
{
  struct htrdr_accum XYZ[3]; /* X, Y, and Z */
  struct htrdr_accum time;
  struct htrdr_planeto* cmd;
  struct planeto_pixel_image* pixel = data;
  size_t ichannel;
  ASSERT(htrdr && htrdr_draw_pixel_args_check(args) && data);
  (void)htrdr;

  cmd = args->context;
  ASSERT(cmd);
  ASSERT(cmd->spectral_domain.spectral_type == HTRDR_SPECTRAL_SW_CIE_XYZ);
  ASSERT(cmd->output_type == HTRDR_PLANETO_ARGS_OUTPUT_IMAGE);

  /* Reset accumulators */
  XYZ[0] = HTRDR_ACCUM_NULL;
  XYZ[1] = HTRDR_ACCUM_NULL;
  XYZ[2] = HTRDR_ACCUM_NULL;
  time = HTRDR_ACCUM_NULL;

  FOR_EACH(ichannel, 0, 3) {
    size_t isamp;

    FOR_EACH(isamp, 0, args->spp) {
      struct time t0, t1;
      struct scam_sample sample = SCAM_SAMPLE_NULL;
      struct scam_ray ray = SCAM_RAY_NULL;
      double weight;
      double r0, r1, r2;
      double wlen[2]; /* Sampled wavelength into the spectral band */
      double pdf;
      size_t iband[2]; /* Sampled spectral band */
      size_t iquad; /* Sampled quadrature point into the spectral band */
      double usec;

      /* Begin the registration of the time spent to in the realisation */
      time_current(&t0);

      /* Sample a position into the pixel, in the normalized image plane */
      sample.film[0] = (double)args->pixel_coord[0]+ssp_rng_canonical(args->rng);
      sample.film[1] = (double)args->pixel_coord[1]+ssp_rng_canonical(args->rng);
      sample.film[0] *= args->pixel_normalized_size[0];
      sample.film[1] *= args->pixel_normalized_size[1];
      sample.lens[0] = ssp_rng_canonical(args->rng);
      sample.lens[1] = ssp_rng_canonical(args->rng);

      /* Generate a camera ray */
      scam_generate_ray(cmd->camera, &sample, &ray);

      r0 = ssp_rng_canonical(args->rng);
      r1 = ssp_rng_canonical(args->rng);
      r2 = ssp_rng_canonical(args->rng);

      /* Sample a wavelength according to the CIE XYZ color space */
      switch(ichannel) {
        case 0:
          wlen[0] = htrdr_ran_wlen_cie_xyz_sample_X(cmd->cie, r0, r1, &pdf);
          break;
        case 1:
          wlen[0] = htrdr_ran_wlen_cie_xyz_sample_Y(cmd->cie, r0, r1, &pdf);
          break;
        case 2:
          wlen[0] = htrdr_ran_wlen_cie_xyz_sample_Z(cmd->cie, r0, r1, &pdf);
          break;
        default: FATAL("Unreachable code.\n"); break;
      }
      pdf *= 1.e9; /* Transform the pdf from nm⁻¹ to m⁻¹ */

      /* Find the band the wavelength belongs to and sample a quadrature point */
      wlen[1] = wlen[0];
      RNATM(find_bands(cmd->atmosphere, wlen, iband));
      RNATM(band_sample_quad_pt(cmd->atmosphere, r2, iband[0], &iquad));
      ASSERT(iband[0] == iband[1]);

      /* TODO Compute the radiance in W/m²/sr/m */
#if 0
      weight = planeto_compute_radiance(cmd, args->ithread, args->rng,
        ray.org, ray.dir, wlen, iband, iquad);
#else
      weight = 0;
#endif
      ASSERT(weight >= 0);

      weight /= pdf; /* In W/m²/sr */

      /* End of time recording per realisation */
      time_sub(&t0, time_current(&t1), &t0);
      usec = (double)time_val(&t0, TIME_NSEC) * 0.001;

      /* Update pixel channel accumulators */
      XYZ[ichannel].sum_weights += weight;
      XYZ[ichannel].sum_weights_sqr += weight*weight;
      XYZ[ichannel].nweights += 1;

      /* Update the per realisation time accumulator */
      time.sum_weights += usec;
      time.sum_weights_sqr += usec*usec;
      time.nweights += 1;
    }
  }

  /* Flush pixel data */
  htrdr_accum_get_estimation(XYZ+0, &pixel->X);
  htrdr_accum_get_estimation(XYZ+1, &pixel->Y);
  htrdr_accum_get_estimation(XYZ+2, &pixel->Z);
  pixel->time = time;
}


static INLINE void
write_accum
  (const struct htrdr_accum* acc, /* Accum to dump */
   struct htrdr_accum* out_acc, /* May be NULL */
   FILE* stream)
{
  ASSERT(acc && stream);

  if(acc->nweights == 0) {
    fprintf(stream, "0 0 ");
  } else {
    struct htrdr_estimate estimate = HTRDR_ESTIMATE_NULL;

    htrdr_accum_get_estimation(acc, &estimate);
    fprintf(stream, "%g %g ", estimate.E, estimate.SE);

    if(out_acc) {
      out_acc->sum_weights += acc->sum_weights;
      out_acc->sum_weights_sqr += acc->sum_weights_sqr;
      out_acc->nweights += acc->nweights;
    }
  }
}

static INLINE void
write_pixel_image
  (const struct planeto_pixel_image* pix,
   struct htrdr_accum* time_acc, /* May be NULL */
   FILE* stream)
{
  ASSERT(pix && stream);
  fprintf(stream, "%g %g ", pix->X.E, pix->X.SE);
  fprintf(stream, "%g %g ", pix->Y.E, pix->Y.SE);
  fprintf(stream, "%g %g ", pix->Z.E, pix->Z.SE);
  write_accum(&pix->time, time_acc, stream);
  fprintf(stream, "\n");
}

static INLINE void
write_pixel_xwave
  (const struct planeto_pixel_xwave* pix,
   struct htrdr_accum* time_acc, /* May be NULL */
   FILE* stream)
{
  ASSERT(pix && stream);
  fprintf(stream, "%g %g %f %f 0 0 ",
    pix->radiance_temperature.E,
    pix->radiance_temperature.SE,
    pix->radiance.E,
    pix->radiance.SE);
  write_accum(&pix->time, time_acc, stream);
  fprintf(stream, "\n");
}

static res_T
write_buffer
  (struct htrdr_planeto* cmd,
   struct htrdr_buffer* buf,
   struct htrdr_accum* time_acc,
   FILE* stream)
{
  struct htrdr_pixel_format pixfmt;
  struct htrdr_buffer_layout layout;
  size_t x, y;
  ASSERT(cmd && buf && time_acc && stream);
  ASSERT(cmd->output_type == HTRDR_PLANETO_ARGS_OUTPUT_IMAGE);

  planeto_get_pixel_format(cmd, &pixfmt);
  htrdr_buffer_get_layout(buf, &layout);
  CHK(pixfmt.size == layout.elmt_size);

  fprintf(stream, "%lu %lu\n", layout.width, layout.height);
  *time_acc = HTRDR_ACCUM_NULL;

  FOR_EACH(y, 0, layout.height) {
    FOR_EACH(x, 0, layout.width) {
      void* pix_raw = htrdr_buffer_at(buf, x, y);
      ASSERT(IS_ALIGNED(pix_raw, pixfmt.alignment));

      switch(cmd->spectral_domain.spectral_type) {
        case HTRDR_SPECTRAL_LW:
        case HTRDR_SPECTRAL_SW:
          write_pixel_xwave(pix_raw, time_acc, stream);
          break;
        case HTRDR_SPECTRAL_SW_CIE_XYZ:
          write_pixel_image(pix_raw, time_acc, stream);
          break;
        default: FATAL("Unreachable code\n"); break;
      }
    }
  }
  return RES_OK;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
planeto_draw_map(struct htrdr_planeto* cmd)
{
  struct htrdr_draw_map_args args = HTRDR_DRAW_MAP_ARGS_NULL;
  struct htrdr_estimate path_time = HTRDR_ESTIMATE_NULL;
  struct htrdr_accum path_time_acc = HTRDR_ACCUM_NULL;
  res_T res = RES_OK;
  ASSERT(cmd && cmd->output_type == HTRDR_PLANETO_ARGS_OUTPUT_IMAGE);

  args.buffer_layout = cmd->buf_layout;
  args.spp = cmd->spp;
  args.context = cmd;
  switch(cmd->spectral_domain.spectral_type) {
    case HTRDR_SPECTRAL_LW:
    case HTRDR_SPECTRAL_SW:
      args.draw_pixel = draw_pixel_xwave;
      break;
    case HTRDR_SPECTRAL_SW_CIE_XYZ:
      args.draw_pixel = draw_pixel_image;
      break;
    default: FATAL("Unreachable code\n"); break;
  }

  res = htrdr_draw_map(cmd->htrdr, &args, cmd->buf);
  if(res != RES_OK) goto error;

  /* Nothing more to do on non master processes */
  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) goto exit;

  /* Write output image */
  res = write_buffer(cmd, cmd->buf, &path_time_acc, cmd->output);
  if(res != RES_OK) goto error;

  /* Log time per realisation */
  htrdr_accum_get_estimation(&path_time_acc, &path_time);
  htrdr_log(cmd->htrdr, "Time per radiative path (in µs): %g +/- %g\n",
    path_time.E, path_time.SE);

exit:
  return res;
error:
  goto exit;
}
