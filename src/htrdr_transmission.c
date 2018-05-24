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

#include "htrdr.h"
#include "htrdr_buffer.h"
#include "htrdr_rectangle.h"
#include "htrdr_solve.h"

#include <star/ssp.h>
#include <star/svx.h>
#include <rsys/math.h>

struct transmit_context {
  double tau;
};
static const struct transmit_context TRANSMIT_CONTEXT_NULL = {0};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static int
discard_hit
  (const struct svx_hit* hit,
   const double org[3],
   const double dir[3],
   const double range[2],
   void* context)
{
  struct transmit_context* ctx = context;
  double dst;
  double k;
  ASSERT(hit && ctx && !SVX_HIT_NONE(hit));
  (void)org, (void)dir, (void)range;

  k = *((double*)hit->voxel.data);
  dst = hit->distance[1] - hit->distance[0];
  ASSERT(dst >= 0);
  ctx->tau += k*dst;
  return 1;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_solve_transmission
  (struct htrdr* htrdr,
   const double pos[3],
   const double dir[3],
   double *val)
{
  struct svx_hit hit = SVX_HIT_NULL;
  struct transmit_context ctx = TRANSMIT_CONTEXT_NULL;
  const double range[2] = {0, INF};
  res_T res = RES_OK;
  ASSERT(htrdr && pos && dir && val);

  res = svx_octree_trace_ray(htrdr->clouds, pos, dir, range, NULL,
    discard_hit, &ctx, &hit);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "error computing the transmission for the position `%g %g %g' "
      "along the direction `%g %g %g'.\n",
      SPLIT3(pos), SPLIT3(dir));
    goto error;
  }
  ASSERT(SVX_HIT_NONE(&hit));
  *val = ctx.tau ? exp(-ctx.tau) : 1.0;

exit:
  return res;
error:
  goto exit;
}

res_T
htrdr_solve_transmission_buffer(struct htrdr* htrdr)
{
  struct htrdr_buffer_layout buf_layout;
  struct ssp_rng* rng = NULL;
  double pixsz[2];
  size_t ix, iy, ispp;
  res_T res = RES_OK;
  ASSERT(htrdr && htrdr->rect && htrdr->buf);

  htrdr_buffer_get_layout(htrdr->buf, &buf_layout);

  res = ssp_rng_create(htrdr->allocator, &ssp_rng_mt19937_64, &rng);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "could not allocate a RNG.\n");
    goto error;
  }

  /* Compute the pixel size in the normalized image space */
  pixsz[0] = 1.0 / (double)buf_layout.width;
  pixsz[1] = 1.0 / (double)buf_layout.height;

  FOR_EACH(iy, 0, buf_layout.height) {
    const double y = (double)iy * pixsz[1]; /* Y in normalized image space */

    FOR_EACH(ix, 0, buf_layout.width) {
      const double x = (double)ix * pixsz[0]; /* X in normalized image space */
      double* val = htrdr_buffer_at(htrdr->buf, ix, iy);
      double accum = 0;

      FOR_EACH(ispp, 0, htrdr->spp) {
        double pixsamp[2]; /* Pixel sample */
        double pos[3]; /* World space position */
        double T;

        /* Uniformaly sample the pixel in the normalized image space */
        pixsamp[0] = x + ssp_rng_canonical(rng)*pixsz[0];
        pixsamp[1] = y + ssp_rng_canonical(rng)*pixsz[1];

        /* Compute the world space position of the sample according to the
         * image plane */
        htrdr_rectangle_sample_pos(htrdr->rect, pixsamp, pos);

        /* Solve the transmission for the sample position wrt the main dir */
        res = htrdr_solve_transmission(htrdr, pos, htrdr->main_dir, &T);
        if(res != RES_OK) goto error;

        accum += T;
      }
      *val = accum / (double)htrdr->spp;
    }
  }

exit:
  if(rng) SSP(rng_ref_put(rng));
  return res;
error:
  goto exit;
}

