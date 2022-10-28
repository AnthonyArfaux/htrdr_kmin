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

#include "atmosphere/htrdr_atmosphere_c.h"
#include "atmosphere/htrdr_atmosphere_ground.h"
#include "atmosphere/htrdr_atmosphere_sun.h"

#include "core/htrdr.h"
#include "core/htrdr_buffer.h"
#include "core/htrdr_draw_map.h"
#include "core/htrdr_interface.h"
#include "core/htrdr_log.h"
#include "core/htrdr_ran_wlen_cie_xyz.h"
#include "core/htrdr_ran_wlen_planck.h"
#include "core/htrdr_rectangle.h"

#include <high_tune/htsky.h>

#include <star/s3d.h>
#include <star/scam.h>
#include <star/ssp.h>

#include <rsys/clock_time.h>
#include <rsys/str.h>

#include <string.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static res_T
sample_rectangle_ray
  (struct htrdr_atmosphere* cmd,
   struct htrdr_rectangle* rect,
   const size_t ipix[2],
   const double pix_sz[2],
   struct ssp_rng* rng,
   double ray_org[3],
   double ray_dir[3])
{
  struct s3d_hit hit = S3D_HIT_NULL;
  double pix_samp[2];
  const double up_dir[3] = {0,0,1};
  const double range[2] = {0, DBL_MAX};
  double normal[3];
  ASSERT(cmd && rect && ipix && pix_sz && rng && ray_org && ray_dir);

  /* Sample a position into the pixel, in the normalized image plane */
  pix_samp[0] = ((double)ipix[0] + ssp_rng_canonical(rng)) * pix_sz[0];
  pix_samp[1] = ((double)ipix[1] + ssp_rng_canonical(rng)) * pix_sz[1];

  /* Retrieve the world space position of pix_samp */
  htrdr_rectangle_sample_pos(rect, pix_samp, ray_org);

  /* Check that `ray_org' is not included into a geometry */
  HTRDR(atmosphere_ground_trace_ray
    (cmd->ground, ray_org, up_dir, range, NULL, &hit));

  /* Up direction is occluded. Check if the sample must be rejected, i.e. does it
   * lies inside a geometry? */
  if(!S3D_HIT_NONE(&hit)) {
    struct htrdr_interface interf = HTRDR_INTERFACE_NULL;
    const struct htrdr_mtl* mtl = NULL;
    float N[3]; /* Normalized normal of the hit */
    float wi[3];
    float cos_wi_N;

    /* Compute the cosine between the up direction and the hit normal */
    f3_set_d3(wi, up_dir);
    f3_normalize(N, hit.normal);
    cos_wi_N = f3_dot(wi, N);

    /* Fetch the hit interface and retrieve the material into which the ray was
     * traced */
    htrdr_atmosphere_ground_get_interface(cmd->ground, &hit, &interf);
    mtl = cos_wi_N < 0 ? &interf.mtl_front : &interf.mtl_back;

    /* Reject the sample if the incident direction do not travel into the sky */
    if(strcmp(mtl->name, cmd->sky_mtl_name) != 0) return RES_BAD_OP;
  }

  /* Sample a ray direction */
  htrdr_rectangle_get_normal(rect, normal);
  ssp_ran_hemisphere_cos(rng, normal, ray_dir, NULL);

  return RES_OK;
}

static void
draw_pixel_image
  (struct htrdr* htrdr,
   const struct htrdr_draw_pixel_args* args,
   void* data)
{
  struct htrdr_accum XYZ[3]; /* X, Y, and Z */
  struct htrdr_accum time;
  struct htrdr_atmosphere* cmd;
  struct atmosphere_pixel_image* pixel = data;
  size_t ichannel;
  ASSERT(htrdr && htrdr_draw_pixel_args_check(args) && data);
  (void)htrdr;

  cmd = args->context;
  ASSERT(cmd);
  ASSERT(cmd->spectral_type == HTRDR_SPECTRAL_SW_CIE_XYZ);
  ASSERT(cmd->output_type == HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE);

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
      double wlen; /* Sampled wavelength into the spectral band */
      double pdf;
      size_t iband; /* Sampled spectral band */
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

      /* Sample a spectral band and a quadrature point */
      switch(ichannel) {
        case 0: wlen = htrdr_ran_wlen_cie_xyz_sample_X(cmd->cie, r0, r1, &pdf); break;
        case 1: wlen = htrdr_ran_wlen_cie_xyz_sample_Y(cmd->cie, r0, r1, &pdf); break;
        case 2: wlen = htrdr_ran_wlen_cie_xyz_sample_Z(cmd->cie, r0, r1, &pdf); break;
        default: FATAL("Unreachable code.\n"); break;
      }
      pdf *= 1.e9; /* Transform the pdf from nm^-1 to m^-1 */

      iband = htsky_find_spectral_band(cmd->sky, wlen);
      iquad = htsky_spectral_band_sample_quadrature(cmd->sky, r2, iband);

      /* Compute the radiance in W/m^2/sr/m */
      weight = atmosphere_compute_radiance_sw(cmd, args->ithread, args->rng,
        ATMOSPHERE_RADIANCE_ALL, ray.org, ray.dir, wlen, iband, iquad);
      ASSERT(weight >= 0);

      weight /= pdf; /* In W/m^2/sr */

      /* End the registration of the per realisation time */
      time_sub(&t0, time_current(&t1), &t0);
      usec = (double)time_val(&t0, TIME_NSEC) * 0.001;

      /* Update the pixel accumulator of the current channel */
      XYZ[ichannel].sum_weights += weight;
      XYZ[ichannel].sum_weights_sqr += weight*weight;
      XYZ[ichannel].nweights += 1;

      /* Update the pixel accumulator of per realisation time */
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

static void
draw_pixel_flux
  (struct htrdr* htrdr,
   const struct htrdr_draw_pixel_args* args,
   void* data)
{
  struct htrdr_accum flux;
  struct htrdr_accum time;
  struct htrdr_atmosphere* cmd;
  struct atmosphere_pixel_flux* pixel = data;
  size_t isamp;
  ASSERT(htrdr && htrdr_draw_pixel_args_check(args) && data);
  (void)htrdr;

  cmd = args->context;
  ASSERT(cmd);
  ASSERT(cmd->output_type == HTRDR_ATMOSPHERE_ARGS_OUTPUT_FLUX_MAP);
  ASSERT(cmd->spectral_type == HTRDR_SPECTRAL_LW
      || cmd->spectral_type == HTRDR_SPECTRAL_SW);

  /* Reset the pixel accumulators */
  flux = HTRDR_ACCUM_NULL;
  time = HTRDR_ACCUM_NULL;

  FOR_EACH(isamp, 0, args->spp) {
    struct time t0, t1;
    double ray_org[3];
    double ray_dir[3];
    double weight;
    double r0, r1, r2;
    double wlen;
    size_t iband;
    size_t iquad;
    double usec;
    double band_pdf;
    res_T res = RES_OK;

    /* Begin the registration of the time spent in the realisation */
    time_current(&t0);

    res = sample_rectangle_ray(cmd, cmd->flux_map, args->pixel_coord,
      args->pixel_normalized_size, args->rng, ray_org, ray_dir);
    if(res != RES_OK) continue; /* Reject the current sample */

    r0 = ssp_rng_canonical(args->rng);
    r1 = ssp_rng_canonical(args->rng);
    r2 = ssp_rng_canonical(args->rng);

    /* Sample a wavelength */
    wlen = htrdr_ran_wlen_planck_sample(cmd->planck, r0, r1, &band_pdf);
    band_pdf *= 1.e9; /* Transform the pdf from nm^-1 to m^-1 */

    /* Select the associated band and sample a quadrature point */
    iband = htsky_find_spectral_band(cmd->sky, wlen);
    iquad = htsky_spectral_band_sample_quadrature(cmd->sky, r2, iband);

    if(cmd->spectral_type == HTRDR_SPECTRAL_LW) {
      weight = atmosphere_compute_radiance_lw(cmd, args->ithread, args->rng,
        ray_org, ray_dir, wlen, iband, iquad);
      weight *= PI / band_pdf; /* Transform weight from W/m^2/sr/m to W/m^2 */
    } else {
      double sun_dir[3];
      double N[3];
      double L_direct;
      double L_diffuse;
      double cos_N_sun_dir;
      double sun_solid_angle;
      ASSERT(cmd->spectral_type == HTRDR_SPECTRAL_SW);

      /* Compute direct contribution if necessary */
      htrdr_atmosphere_sun_sample_direction(cmd->sun, args->rng, sun_dir);
      htrdr_rectangle_get_normal(cmd->flux_map, N);
      cos_N_sun_dir = d3_dot(N, sun_dir);

      if(cos_N_sun_dir <= 0) {
        L_direct = 0;
      } else {
        L_direct = atmosphere_compute_radiance_sw(cmd, args->ithread,
          args->rng, ATMOSPHERE_RADIANCE_DIRECT, ray_org, sun_dir, wlen, iband,
          iquad);
      }

      /* Compute diffuse contribution */
      L_diffuse = atmosphere_compute_radiance_sw(cmd, args->ithread, args->rng,
        ATMOSPHERE_RADIANCE_DIFFUSE, ray_org, ray_dir, wlen, iband, iquad);

      sun_solid_angle = htrdr_atmosphere_sun_get_solid_angle(cmd->sun);

      /* Compute the weight in W/m^2/m */
      weight = cos_N_sun_dir * sun_solid_angle * L_direct + PI * L_diffuse;

      /* Importance sampling: correct weight with pdf */
      weight /= band_pdf; /* In W/m^2 */
    }

    /* End the registration of the per realisation time */
    time_sub(&t0, time_current(&t1), &t0);
    usec = (double)time_val(&t0, TIME_NSEC) * 0.001;

    /* Update the pixel accumulator of the flux */
    flux.sum_weights += weight;
    flux.sum_weights_sqr += weight*weight;
    flux.nweights += 1;

    /* Update the pixel accumulator of per realisation time */
    time.sum_weights += usec;
    time.sum_weights_sqr += usec*usec;
    time.nweights += 1;
  }

  /* Write the accumulators */
  pixel->flux = flux;
  pixel->time = time;
}

static void
draw_pixel_xwave
  (struct htrdr* htrdr,
   const struct htrdr_draw_pixel_args* args,
   void* data)
{
  struct htrdr_accum radiance;
  struct htrdr_accum time;
  struct htrdr_atmosphere* cmd;
  struct atmosphere_pixel_xwave* pixel = data;
  size_t isamp;
  double temp_min, temp_max;
  ASSERT(htrdr && htrdr_draw_pixel_args_check(args) && data);
  (void)htrdr;

  cmd = args->context;
  ASSERT(cmd->output_type == HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE);
  ASSERT(cmd->spectral_type == HTRDR_SPECTRAL_LW
      || cmd->spectral_type == HTRDR_SPECTRAL_SW);

  /* Reset the pixel accumulators */
  radiance = HTRDR_ACCUM_NULL;
  time = HTRDR_ACCUM_NULL;

  FOR_EACH(isamp, 0, args->spp) {
    struct time t0, t1;
    struct scam_sample sample = SCAM_SAMPLE_NULL;
    struct scam_ray ray = SCAM_RAY_NULL;
    double weight;
    double r0, r1, r2;
    double wlen;
    size_t iband;
    size_t iquad;
    double usec;
    double band_pdf;

    /* Begin the registration of the time spent in the realisation */
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
    wlen = htrdr_ran_wlen_planck_sample(cmd->planck, r0, r1, &band_pdf);
    band_pdf *= 1.e9; /* Transform the pdf from nm^-1 to m^-1 */

    /* Select the associated band and sample a quadrature point */
    iband = htsky_find_spectral_band(cmd->sky, wlen);
    iquad = htsky_spectral_band_sample_quadrature(cmd->sky, r2, iband);

    /* Compute the spectral radiance in W/m^2/sr/m */
    switch(cmd->spectral_type) {
      case HTRDR_SPECTRAL_LW:
        weight = atmosphere_compute_radiance_lw(cmd, args->ithread, args->rng,
          ray.org, ray.dir, wlen, iband, iquad);
        break;
      case HTRDR_SPECTRAL_SW:
        weight = atmosphere_compute_radiance_sw(cmd, args->ithread, args->rng,
          ATMOSPHERE_RADIANCE_ALL, ray.org, ray.dir, wlen, iband, iquad);
        break;
      default: FATAL("Unreachable code.\n"); break;
    }
    ASSERT(weight >= 0);
    /* Importance sampling: correct weight with pdf */
    weight /= band_pdf; /* In W/m^2/sr */

    /* End the registration of the per realisation time */
    time_sub(&t0, time_current(&t1), &t0);
    usec = (double)time_val(&t0, TIME_NSEC) * 0.001;

    /* Update the pixel accumulator */
    radiance.sum_weights += weight;
    radiance.sum_weights_sqr += weight*weight;
    radiance.nweights += 1;

    /* Update the pixel accumulator of per realisation time */
    time.sum_weights += usec;
    time.sum_weights_sqr += usec*usec;
    time.nweights += 1;
  }

  /* Compute the estimation of the pixel radiance */
  htrdr_accum_get_estimation(&radiance, &pixel->radiance);

  /* Save the per realisation integration time */
  pixel->time = time;

  /* Compute the brightness_temperature of the pixel and estimate its standard
   * error if the sources were in the medium (<=> longwave) */
  if(cmd->spectral_type == HTRDR_SPECTRAL_LW) {
    const double wlen_min = cmd->wlen_range_m[0];
    const double wlen_max = cmd->wlen_range_m[1];
    pixel->radiance_temperature.E = htrdr_radiance_temperature
      (cmd->htrdr, wlen_min, wlen_max, pixel->radiance.E);
    temp_min = htrdr_radiance_temperature
      (cmd->htrdr, wlen_min, wlen_max, pixel->radiance.E - pixel->radiance.SE);
    temp_max = htrdr_radiance_temperature
      (cmd->htrdr, wlen_min, wlen_max, pixel->radiance.E + pixel->radiance.SE);
    pixel->radiance_temperature.SE = temp_max - temp_min;
  }
}

static INLINE void
setup_draw_map_args_rectangle
  (struct htrdr_atmosphere* cmd,
   struct htrdr_draw_map_args* args)
{
  ASSERT(cmd && args);
  ASSERT(cmd->output_type == HTRDR_ATMOSPHERE_ARGS_OUTPUT_FLUX_MAP);
  *args = HTRDR_DRAW_MAP_ARGS_NULL;
  args->draw_pixel = draw_pixel_flux;
  args->buffer_layout = cmd->buf_layout;
  args->spp = cmd->spp;
  args->context = cmd;
}

static INLINE void
setup_draw_map_args_camera
  (struct htrdr_atmosphere* cmd,
   struct htrdr_draw_map_args* args)
{
  ASSERT(cmd && args);
  ASSERT(cmd->output_type == HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE);

  *args = HTRDR_DRAW_MAP_ARGS_NULL;
  args->buffer_layout = cmd->buf_layout;
  args->spp = cmd->spp;
  args->context = cmd;

  switch(cmd->spectral_type) {
    case HTRDR_SPECTRAL_LW:
    case HTRDR_SPECTRAL_SW:
      args->draw_pixel = draw_pixel_xwave;
      break;
    case HTRDR_SPECTRAL_SW_CIE_XYZ:
      args->draw_pixel = draw_pixel_image;
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
}

static INLINE void
dump_accum
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
dump_pixel_flux
  (const struct atmosphere_pixel_flux* pix,
   struct htrdr_accum* time_acc, /* May be NULL */
   struct htrdr_accum* flux_acc, /* May be NULL */
   FILE* stream)
{
  ASSERT(pix && stream);
  dump_accum(&pix->flux, flux_acc, stream);
  fprintf(stream, "0 0 0 0 ");
  dump_accum(&pix->time, time_acc, stream);
  fprintf(stream, "\n");
}

static INLINE void
dump_pixel_image
  (const struct atmosphere_pixel_image* pix,
   struct htrdr_accum* time_acc, /* May be NULL */
   FILE* stream)
{
  ASSERT(pix && stream);
  fprintf(stream, "%g %g ", pix->X.E, pix->X.SE);
  fprintf(stream, "%g %g ", pix->Y.E, pix->Y.SE);
  fprintf(stream, "%g %g ", pix->Z.E, pix->Z.SE);
  dump_accum(&pix->time, time_acc, stream);
  fprintf(stream, "\n");
}

static INLINE void
dump_pixel_xwave
  (const struct atmosphere_pixel_xwave* pix,
   struct htrdr_accum* time_acc, /* May be NULL */
   FILE* stream)
{
  ASSERT(pix && stream);
  fprintf(stream, "%g %g %f %f 0 0 ",
    pix->radiance_temperature.E,
    pix->radiance_temperature.SE,
    pix->radiance.E,
    pix->radiance.SE);
  dump_accum(&pix->time, time_acc, stream);
  fprintf(stream, "\n");
}

static res_T
dump_buffer
  (struct htrdr_atmosphere* cmd,
   struct htrdr_buffer* buf,
   struct htrdr_accum* time_acc, /* May be NULL */
   struct htrdr_accum* flux_acc, /* May be NULL */
   FILE* stream)
{
  struct htrdr_pixel_format pixfmt;
  struct htrdr_buffer_layout layout;
  size_t x, y;
  ASSERT(cmd && buf && stream);

  atmosphere_get_pixel_format(cmd, &pixfmt);
  htrdr_buffer_get_layout(buf, &layout);
  CHK(pixfmt.size == layout.elmt_size);

  fprintf(stream, "%lu %lu\n", layout.width, layout.height);

  if(time_acc) *time_acc = HTRDR_ACCUM_NULL;
  if(flux_acc) *flux_acc = HTRDR_ACCUM_NULL;

  FOR_EACH(y, 0, layout.height) {
    FOR_EACH(x, 0, layout.width) {
      void* pix_raw = htrdr_buffer_at(buf, x, y);
      ASSERT(IS_ALIGNED(pix_raw, pixfmt.alignment));

      if(cmd->output_type == HTRDR_ATMOSPHERE_ARGS_OUTPUT_FLUX_MAP) {
        const struct atmosphere_pixel_flux* pix = pix_raw;
        dump_pixel_flux(pix, time_acc, flux_acc, stream);
      } else if(cmd->spectral_type == HTRDR_SPECTRAL_SW_CIE_XYZ) {
        const struct atmosphere_pixel_image* pix = pix_raw;
        dump_pixel_image(pix, time_acc, stream);
      } else {
        const struct atmosphere_pixel_xwave* pix = pix_raw;
        dump_pixel_xwave(pix, time_acc, stream);
      }
    }
    fprintf(stream, "\n");
  }

  return RES_OK;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
atmosphere_draw_map(struct htrdr_atmosphere* cmd)
{
  struct htrdr_draw_map_args args = HTRDR_DRAW_MAP_ARGS_NULL;
  struct htrdr_accum path_time_acc = HTRDR_ACCUM_NULL;
  struct htrdr_accum flux_acc = HTRDR_ACCUM_NULL;
  struct htrdr_estimate path_time;
  struct htrdr_estimate flux;
  res_T res = RES_OK;
  ASSERT(cmd);

  args.spp = cmd->spp;

  switch(cmd->output_type) {
    case HTRDR_ATMOSPHERE_ARGS_OUTPUT_FLUX_MAP:
      setup_draw_map_args_rectangle(cmd, &args);
      break;
    case HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE:
      setup_draw_map_args_camera(cmd, &args);
      break;
    default: FATAL("Unreachable code.\n"); break;
  }

  res = htrdr_draw_map(cmd->htrdr, &args, cmd->buf);
  if(res != RES_OK) goto error;

  /* No more to do on non master processes */
  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) goto exit;

  /* Write buffer to output */
  res = dump_buffer(cmd, cmd->buf, &path_time_acc, &flux_acc, cmd->output);
  if(res != RES_OK) goto error;

  htrdr_accum_get_estimation(&path_time_acc, &path_time);
  htrdr_log(cmd->htrdr,
    "Time per radiative path (in micro seconds): %g +/- %g\n",
    path_time.E, path_time.SE);

  if(cmd->output_type == HTRDR_ATMOSPHERE_ARGS_OUTPUT_FLUX_MAP) {
    htrdr_accum_get_estimation(&flux_acc, &flux);
    htrdr_log(cmd->htrdr,
      "Radiative flux density (in W/(external m^2)): %g +/- %g\n",
      flux.E, flux.SE);
  }

exit:
  return res;
error:
  goto exit;
}

