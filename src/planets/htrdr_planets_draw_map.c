/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2023 Observatoire de Paris
 * Copyright (C) 2022-2023 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2023 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2023 Université Paul Sabatier
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

#include "planets/htrdr_planets_c.h"
#include "planets/htrdr_planets_source.h"

#include "core/htrdr.h"
#include "core/htrdr_accum.h"
#include "core/htrdr_buffer.h"
#include "core/htrdr_draw_map.h"
#include "core/htrdr_log.h"
#include "core/htrdr_ran_wlen_cie_xyz.h"
#include "core/htrdr_ran_wlen_discrete.h"
#include "core/htrdr_ran_wlen_planck.h"

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
  struct planets_compute_radiance_args rad_args =
    PLANETS_COMPUTE_RADIANCE_ARGS_NULL;

  struct htrdr_accum radiance;
  struct htrdr_accum time;
  struct htrdr_planets* cmd;
  struct planets_pixel_xwave* pixel = data;
  size_t isamp = 0;
  ASSERT(htrdr && htrdr_draw_pixel_args_check(args) && data);
  (void)htrdr;

  cmd  = args->context;
  ASSERT(cmd);
  ASSERT(cmd->spectral_domain.type == HTRDR_SPECTRAL_SW
      || cmd->spectral_domain.type == HTRDR_SPECTRAL_LW);
  ASSERT(cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_IMAGE);

  /* Reset accumulators */
  radiance = HTRDR_ACCUM_NULL;
  time = HTRDR_ACCUM_NULL;

  FOR_EACH(isamp, 0, args->spp) {
    struct time t0, t1;
    struct scam_sample sample = SCAM_SAMPLE_NULL;
    struct scam_ray ray = SCAM_RAY_NULL;
    double weight;
    double r0, r1, r2;
    double wlen[2]; /* Sampled wavelength */
    double pdf;
    size_t iband[2];
    size_t iquad;
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

    /* Sample a wavelength */
    switch(cmd->spectral_domain.type) {
      case HTRDR_SPECTRAL_LW:
        wlen[0] = htrdr_ran_wlen_planck_sample(cmd->planck, r0, r1, &pdf);
        break;
      case HTRDR_SPECTRAL_SW:
        if(htrdr_planets_source_does_radiance_vary_spectrally(cmd->source)) {
          wlen[0] = htrdr_ran_wlen_discrete_sample(cmd->discrete, r0, r1, &pdf);
        } else {
          wlen[0] = htrdr_ran_wlen_planck_sample(cmd->planck, r0, r1, &pdf);
        }
        break;
      default: FATAL("Unreachable code\n"); break;

    }
    wlen[1] = wlen[0];
    pdf *= 1.e9; /* Transform the pdf from nm⁻¹ to m⁻¹ */

    /* Find the band the wavelength belongs to and sample a quadrature point */
    RNATM(find_bands(cmd->atmosphere, wlen, iband));
    RNATM(band_sample_quad_pt(cmd->atmosphere, r2, iband[0], &iquad));
    ASSERT(iband[0] == iband[1]);

    /* Compute the radiance in W/m²/sr/m */
    d3_set(rad_args.path_org, ray.org);
    d3_set(rad_args.path_dir, ray.dir);
    rad_args.rng = args->rng;
    rad_args.ithread = args->ithread;
    rad_args.wlen = wlen[0];
    rad_args.iband = iband[0];
    rad_args.iquad = iquad;
    weight = planets_compute_radiance(cmd, &rad_args);
    ASSERT(weight >= 0);

    weight /= pdf; /* In W/m²/sr */

    /* End of time recording per realisation */
    time_sub(&t0, time_current(&t1), &t0);
    usec = (double)time_val(&t0, TIME_NSEC) * 0.001;

    /* Update the radiance accumulator */
    radiance.sum_weights += weight;
    radiance.sum_weights_sqr += weight*weight;
    radiance.nweights += 1;

    /* Update the per realisation time accumulator */
    time.sum_weights += usec;
    time.sum_weights_sqr += usec*usec;
    time.nweights += 1;
  }

  /* Flush pixel data */
  pixel->radiance = radiance;
  pixel->time = time;

  if(cmd->spectral_domain.type == HTRDR_SPECTRAL_SW) {
    pixel->radiance_temperature.E = 0;
    pixel->radiance_temperature.SE = 0;
  } else {
    struct htrdr_estimate L;

    /* Wavelength in meters */
    const double wmin_m = cmd->spectral_domain.wlen_range[0] * 1.e-9;
    const double wmax_m = cmd->spectral_domain.wlen_range[1] * 1.e-9;

    /* Brightness temperatures in W/m²/sr */
    double T, Tmin, Tmax;

    htrdr_accum_get_estimation(&radiance, &L);

    T    = htrdr_radiance_temperature(cmd->htrdr, wmin_m, wmax_m, L.E);
    Tmin = htrdr_radiance_temperature(cmd->htrdr, wmin_m, wmax_m, L.E - L.SE);
    Tmax = htrdr_radiance_temperature(cmd->htrdr, wmin_m, wmax_m, L.E + L.SE);

    pixel->radiance_temperature.E = T;
    pixel->radiance_temperature.SE = Tmax - Tmin;
  }
}

static void
draw_pixel_image
  (struct htrdr* htrdr,
   const struct htrdr_draw_pixel_args* args,
   void* data)
{
  struct planets_compute_radiance_args rad_args =
    PLANETS_COMPUTE_RADIANCE_ARGS_NULL;

  struct htrdr_accum XYZ[3]; /* X, Y, and Z */
  struct htrdr_accum time;
  struct htrdr_planets* cmd;
  struct planets_pixel_image* pixel = data;
  size_t ichannel;
  ASSERT(htrdr && htrdr_draw_pixel_args_check(args) && data);
  (void)htrdr;

  cmd = args->context;
  ASSERT(cmd);
  ASSERT(cmd->spectral_domain.type == HTRDR_SPECTRAL_SW_CIE_XYZ);
  ASSERT(cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_IMAGE);

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
      SCAM(generate_ray(cmd->camera, &sample, &ray));

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

      /* Compute the radiance in W/m²/sr/m */
      d3_set(rad_args.path_org, ray.org);
      d3_set(rad_args.path_dir, ray.dir);
      rad_args.rng = args->rng;
      rad_args.ithread = args->ithread;
      rad_args.wlen = wlen[0];
      rad_args.iband = iband[0];
      rad_args.iquad = iquad;
      weight = planets_compute_radiance(cmd, &rad_args);
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
  (const struct planets_pixel_image* pix,
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
  (const struct planets_pixel_xwave* pix,
   struct htrdr_accum* radiance_acc, /* May be NULL */
   struct htrdr_accum* time_acc, /* May be NULL */
   FILE* stream)
{
  ASSERT(pix && stream);
  fprintf(stream, "%g %g ",
    pix->radiance_temperature.E,
    pix->radiance_temperature.SE);
  write_accum(&pix->radiance, radiance_acc, stream);
  fprintf(stream, "0 0 ");
  write_accum(&pix->time, time_acc, stream);
  fprintf(stream, "\n");
}

static res_T
write_buffer
  (struct htrdr_planets* cmd,
   struct htrdr_buffer* buf,
   struct htrdr_accum* radiance_acc, /* May be NULL */
   struct htrdr_accum* time_acc,
   FILE* stream)
{
  struct htrdr_pixel_format pixfmt;
  struct htrdr_buffer_layout layout;
  size_t x, y;
  ASSERT(cmd && buf && time_acc && stream);
  ASSERT(cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_IMAGE);

  planets_get_pixel_format(cmd, &pixfmt);
  htrdr_buffer_get_layout(buf, &layout);
  CHK(pixfmt.size == layout.elmt_size);

  fprintf(stream, "%lu %lu\n", layout.width, layout.height);
  *time_acc = HTRDR_ACCUM_NULL;

  FOR_EACH(y, 0, layout.height) {
    FOR_EACH(x, 0, layout.width) {
      void* pix_raw = htrdr_buffer_at(buf, x, y);
      ASSERT(IS_ALIGNED(pix_raw, pixfmt.alignment));

      switch(cmd->spectral_domain.type) {
        case HTRDR_SPECTRAL_LW:
        case HTRDR_SPECTRAL_SW:
          write_pixel_xwave(pix_raw, radiance_acc, time_acc, stream);
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
planets_draw_map(struct htrdr_planets* cmd)
{
  struct htrdr_draw_map_args args = HTRDR_DRAW_MAP_ARGS_NULL;
  struct htrdr_estimate path_time = HTRDR_ESTIMATE_NULL;
  struct htrdr_accum path_time_acc = HTRDR_ACCUM_NULL;
  struct htrdr_accum radiance_acc = HTRDR_ACCUM_NULL;
  res_T res = RES_OK;
  ASSERT(cmd && cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_IMAGE);

  args.buffer_layout = cmd->buf_layout;
  args.spp = cmd->spp;
  args.context = cmd;
  switch(cmd->spectral_domain.type) {
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
  res = write_buffer(cmd, cmd->buf, &radiance_acc, &path_time_acc, cmd->output);
  if(res != RES_OK) goto error;

  CHK(fflush(cmd->output) == 0);

  /* Log time per realisation */
  htrdr_accum_get_estimation(&path_time_acc, &path_time);
  htrdr_log(cmd->htrdr, "Time per radiative path (in µs): %g +/- %g\n",
    path_time.E, path_time.SE);

  /* Log measured radiance on the whole image */
  if(cmd->spectral_domain.type == HTRDR_SPECTRAL_LW
  || cmd->spectral_domain.type == HTRDR_SPECTRAL_SW) {
    struct htrdr_estimate L;
    double omega; /* Solid angle of the camera */

    htrdr_accum_get_estimation(&radiance_acc, &L);
    SCAM(perspective_get_solid_angle(cmd->camera, &omega));
    htrdr_log(cmd->htrdr, "Radiance in W/m²/sr: %g +/- %g\n", L.E, L.SE);
    htrdr_log(cmd->htrdr, "Radiance in W/m² (solid angle = %g sr): %g +/- %g\n",
      omega, L.E*omega, L.SE*omega);
  }

exit:
  return res;
error:
  goto exit;
}
