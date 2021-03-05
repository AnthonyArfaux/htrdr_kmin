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

#include "combustion/htrdr_combustion_c.h"

#include "core/htrdr_accum.h"
#include "core/htrdr_camera.h"
#include "core/htrdr_draw_map.h"
#include "core/htrdr_log.h"

#include <star/ssp.h>

#include <rsys/clock_time.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
draw_pixel
  (struct htrdr* htrdr,
   const struct htrdr_draw_pixel_args* args,
   void* data)
{
  struct htrdr_accum radiance = HTRDR_ACCUM_NULL;
  struct htrdr_accum time = HTRDR_ACCUM_NULL;
  struct htrdr_combustion* cmd = NULL;
  struct combustion_pixel* pixel = NULL;
  size_t isamp;
  ASSERT(htrdr && htrdr_draw_pixel_args_check(args) && data);
  (void)htrdr;

  cmd = args->context;
  pixel = data;

  FOR_EACH(isamp, 0, args->spp) {
    struct time t0, t1;
    double pix_samp[2];
    double ray_org[3];
    double ray_dir[3];
    double weight;
    double usec;

    /* Begin the registration of the time spent in the realisation */
    time_current(&t0);

    /* Sample a position into the pixel, in the normalized image plane */
    pix_samp[0] = (double)args->pixel_coord[0] + ssp_rng_canonical(args->rng);
    pix_samp[1] = (double)args->pixel_coord[1] + ssp_rng_canonical(args->rng);
    pix_samp[0] *= args->pixel_normalized_size[0];
    pix_samp[1] *= args->pixel_normalized_size[1];

    /* Generate a ray starting from the pinhole camera and passing through the
     * pixel sample */
    htrdr_camera_ray(cmd->camera, pix_samp, ray_org, ray_dir);

    /* Backward trace the path */
    weight = combustion_compute_radiance_sw(cmd, args->ithread, args->rng,
        ray_org, ray_dir, cmd->wavelength);

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
dump_pixel
  (const struct combustion_pixel* pix,
   struct htrdr_accum* time_acc, /* May be NULL */
   FILE* stream)
{
  struct htrdr_estimate time = HTRDR_ESTIMATE_NULL;
  ASSERT(pix && stream && pix->time.nweights);

  htrdr_accum_get_estimation(&pix->time, &time);

  fprintf(stream, "%g %g 0 0 0 0 %g %g\n",
    pix->radiance.E, pix->radiance.SE, time.E,time.SE);

  if(time_acc) {
    time_acc->sum_weights += pix->time.sum_weights;
    time_acc->sum_weights_sqr += pix->time.sum_weights_sqr;
    time_acc->nweights += pix->time.nweights;
  }
}

static res_T
dump_buffer
  (struct htrdr_combustion* cmd,
   struct htrdr_buffer* buf,
   struct htrdr_accum* time_acc, /* May be NULL */
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
      const struct combustion_pixel* pix = htrdr_buffer_at(buf, x, y);
      ASSERT(IS_ALIGNED(pix, pixfmt.alignment));
      dump_pixel(pix, time_acc, stream);
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
  struct htrdr_estimate path_time = HTRDR_ESTIMATE_NULL;
  res_T res = RES_OK;
  ASSERT(cmd);

  /* Setup the draw map input arguments */
  args.draw_pixel = draw_pixel;
  args.buffer_layout = cmd->buf_layout;
  args.spp = cmd->spp;
  args.context = cmd;

  res = htrdr_draw_map(cmd->htrdr, &args, cmd->buf);
  if(res != RES_OK) goto error;

  /* No more to do on non master processes */
  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) goto exit;

  /* Write buffer to output */
  res = dump_buffer(cmd, cmd->buf, &path_time_acc, cmd->output);
  if(res != RES_OK) goto error;

  htrdr_accum_get_estimation(&path_time_acc, &path_time);
  htrdr_log(cmd->htrdr,
    "Time per radiative path (in micro seconds): %g +/- %g\n",
    path_time.E, path_time.SE);

exit:
  return res;
error:
  goto exit;
}

