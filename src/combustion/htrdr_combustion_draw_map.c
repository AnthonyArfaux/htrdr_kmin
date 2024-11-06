/* Copyright (C) 2018-2019, 2022-2024 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2024 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2024 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2024 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2024 Observatoire de Paris
 * Copyright (C) 2022-2024 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2024 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2024 Université Paul Sabatier
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

#include "combustion/htrdr_combustion_c.h"

#include "core/htrdr_accum.h"
#include "core/htrdr_draw_map.h"
#include "core/htrdr_log.h"
#include "core/htrdr_rectangle.h"

#include <star/scam.h>
#include <star/ssp.h>

#include <rsys/clock_time.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
draw_pixel_image
  (struct htrdr* htrdr,
   const struct htrdr_draw_pixel_args* args,
   void* data)
{
  struct htrdr_accum radiance = HTRDR_ACCUM_NULL;
  struct htrdr_accum time = HTRDR_ACCUM_NULL;
  struct htrdr_combustion* cmd = NULL;
  struct combustion_pixel_image* pixel = NULL;
  size_t isamp;
  res_T res = RES_OK;
  ASSERT(htrdr && htrdr_draw_pixel_args_check(args) && data);
  (void)htrdr;

  cmd = args->context;
  pixel = data;

  FOR_EACH(isamp, 0, args->spp) {
    struct time t0, t1;
    struct scam_sample sample = SCAM_SAMPLE_NULL;
    struct scam_ray ray = SCAM_RAY_NULL;
    double weight;
    double usec;

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

    /* Backward trace the path */
    res = combustion_compute_radiance_sw(cmd, args->ithread, args->rng,
        ray.org, ray.dir, &weight);
    if(res != RES_OK) continue; /* Reject the path */

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

  /* Save the per realisation time */
  pixel->time = time;
}

static void
draw_pixel_flux
  (struct htrdr* htrdr,
   const struct htrdr_draw_pixel_args* args,
   void* data)
{
  struct htrdr_accum flux = HTRDR_ACCUM_NULL;
  struct htrdr_accum time = HTRDR_ACCUM_NULL;
  struct htrdr_combustion* cmd = NULL;
  struct combustion_pixel_flux* pixel = NULL;
  size_t isamp;
  res_T res = RES_OK;
  ASSERT(htrdr && htrdr_draw_pixel_args_check(args) && data);
  (void)htrdr;

  cmd = args->context;
  pixel = data;

  FOR_EACH(isamp, 0, args->spp) {
    struct time t0, t1;
    double pix_samp[2];
    double ray_org[3];
    double ray_dir[3];
    double normal[3];
    double weight;
    double usec;

    /* Begin the registration of the time spent in the realisation */
    time_current(&t0);

    /* Sample a position into the pixel, in the normalized image plane */
    pix_samp[0] = (double)args->pixel_coord[0] + ssp_rng_canonical(args->rng);
    pix_samp[1] = (double)args->pixel_coord[1] + ssp_rng_canonical(args->rng);
    pix_samp[0] *= args->pixel_normalized_size[0];
    pix_samp[1] *= args->pixel_normalized_size[1];

    /* Retrieve the world space position of pix_samp */
    htrdr_rectangle_sample_pos(cmd->flux_map, pix_samp, ray_org);

    /* Sample a ray direction */
    htrdr_rectangle_get_normal(cmd->flux_map, normal);
    ssp_ran_hemisphere_cos(args->rng, normal, ray_dir, NULL);

    /* Backward trace the path */
    res = combustion_compute_radiance_sw(cmd, args->ithread,
      args->rng, ray_org, ray_dir, &weight);
    if(res != RES_OK) continue; /* Reject the path */
    weight *= PI; /* Transform form W/m^2/sr to W/m^2 */

    /* End the registration of the per realisation time */
    time_sub(&t0, time_current(&t1), &t0);
    usec = (double)time_val(&t0, TIME_NSEC) * 0.001;

    /* Update the pixel accumulator */
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

static void
dump_pixel_flux
  (const struct combustion_pixel_flux* pix,
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

static void
dump_pixel_image
  (const struct combustion_pixel_image* pix,
   struct htrdr_accum* time_acc, /* May be NULL */
   FILE* stream)
{
  ASSERT(pix && stream);
  fprintf(stream, "%g %g ", pix->radiance.E, pix->radiance.SE);
  fprintf(stream, "0 0 0 0 ");
  dump_accum(&pix->time, time_acc, stream);
  fprintf(stream, "\n");
}

static res_T
dump_buffer
  (struct htrdr_combustion* cmd,
   struct htrdr_buffer* buf,
   struct htrdr_accum* time_acc, /* May be NULL */
   struct htrdr_accum* flux_acc, /* May be NULL */
   FILE* stream)
{
  struct htrdr_pixel_format pixfmt;
  struct htrdr_buffer_layout layout;
  size_t x, y;
  ASSERT(cmd && buf && stream);

  combustion_get_pixel_format(cmd, &pixfmt);
  htrdr_buffer_get_layout(buf, &layout);
  CHK(pixfmt.size == layout.elmt_size);

  fprintf(stream, "%lu %lu\n", layout.width, layout.height);

  if(time_acc) *time_acc = HTRDR_ACCUM_NULL;

  FOR_EACH(y, 0, layout.height) {
    FOR_EACH(x, 0, layout.width) {
      void* pix_raw = htrdr_buffer_at(buf, x, y);
      ASSERT(IS_ALIGNED(pix_raw, pixfmt.alignment));
      if(cmd->output_type == HTRDR_COMBUSTION_ARGS_OUTPUT_FLUX_MAP) {
        const struct combustion_pixel_flux* pix = pix_raw;
        dump_pixel_flux(pix, time_acc, flux_acc, stream);
      } else if(cmd->output_type == HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE) {
        const struct combustion_pixel_image* pix = pix_raw;
        dump_pixel_image(pix, time_acc, stream);
      } else {
        FATAL("Unreachable code.\n");
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
combustion_draw_map(struct htrdr_combustion* cmd)
{
  struct htrdr_draw_map_args args = HTRDR_DRAW_MAP_ARGS_NULL;
  struct htrdr_accum path_time_acc = HTRDR_ACCUM_NULL;
  struct htrdr_accum flux_acc = HTRDR_ACCUM_NULL;
  struct htrdr_estimate path_time = HTRDR_ESTIMATE_NULL;
  struct htrdr_estimate flux = HTRDR_ESTIMATE_NULL;
  res_T res = RES_OK;
  ASSERT(cmd);

  /* Setup the draw map input arguments */
  switch(cmd->output_type) {
    case HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE:
      args.draw_pixel = draw_pixel_image;
      break;
    case HTRDR_COMBUSTION_ARGS_OUTPUT_FLUX_MAP:
      args.draw_pixel = draw_pixel_flux;
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
  args.buffer_layout = cmd->buf_layout;
  args.spp = cmd->spp;
  args.context = cmd;

  res = htrdr_draw_map(cmd->htrdr, &args, cmd->buf);
  if(res != RES_OK) goto error;

  /* Nothing more to do on non master processes */
  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) goto exit;

  /* Write buffer to output */
  res = dump_buffer(cmd, cmd->buf, &path_time_acc, &flux_acc, cmd->output);
  if(res != RES_OK) goto error;

  htrdr_accum_get_estimation(&path_time_acc, &path_time);
  htrdr_log(cmd->htrdr,
    "Time per radiative path (in micro seconds): %g +/- %g\n",
    path_time.E, path_time.SE);

  if(cmd->output_type == HTRDR_COMBUSTION_ARGS_OUTPUT_FLUX_MAP) {
    double map_size[2];
    double area;

    htrdr_accum_get_estimation(&flux_acc, &flux);
    htrdr_log(cmd->htrdr,
      "Radiative flux density (in W/m^2): %g +/- %g\n",
      flux.E, flux.SE);

    htrdr_rectangle_get_size(cmd->flux_map, map_size);
    area = map_size[0] * map_size[1];
    htrdr_log(cmd->htrdr,
      "Radiative flux (in W): %g +/- %g\n",
      flux.E*area, flux.SE*area);
  }

exit:
  return res;
error:
  goto exit;
}

