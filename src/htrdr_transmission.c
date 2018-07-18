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
#include "htrdr_sky.h"
#include "htrdr_solve.h"

#include <star/ssp.h>
#include <star/svx.h>
#include <rsys/math.h>

#include <omp.h>

struct transmit_context {
  const struct htrdr_sky* sky;
  double tau_sampled;
  double tau_max_min;
  double tau_min;
};
static const struct transmit_context TRANSMIT_CONTEXT_NULL = {NULL,0,0,0};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static FINLINE uint16_t
morton2D_decode(const uint32_t u32)
{
  uint32_t x = u32 & 0x55555555;
  x = (x | (x >> 1)) & 0x33333333;
  x = (x | (x >> 2)) & 0x0F0F0F0F;
  x = (x | (x >> 4)) & 0x00FF00FF;
  x = (x | (x >> 8)) & 0x0000FFFF;
  return (uint16_t)x;
}

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
  double k_ext_min;
  double k_ext_max;
  ASSERT(hit && ctx && !SVX_HIT_NONE(hit));
  (void)org, (void)dir, (void)range;

  k_ext_min = htrdr_sky_fetch_svx_voxel_property(ctx->sky, HTRDR_Kext,
    HTRDR_SVX_MIN, HTRDR_ALL_COMPONENTS, -1/*FIXME*/, &hit->voxel);
  k_ext_max = htrdr_sky_fetch_svx_voxel_property(ctx->sky, HTRDR_Kext,
    HTRDR_SVX_MAX, HTRDR_ALL_COMPONENTS, -1/*FIXME*/, &hit->voxel);

  dst = hit->distance[1] - hit->distance[0];
  ASSERT(dst >= 0);
  ctx->tau_min += k_ext_min*dst;
  ctx->tau_max_min += (k_ext_max - k_ext_min)*dst;
  return ctx->tau_max_min < ctx->tau_sampled;
}

static res_T
transmission_realisation
  (struct htrdr* htrdr,
   struct ssp_rng* rng,
   const double pos[3],
   const double dir[3],
   double *val)
{
  struct svx_hit hit = SVX_HIT_NULL;
  struct svx_tree* svx_tree = NULL;
  struct transmit_context ctx = TRANSMIT_CONTEXT_NULL;
  const double range[2] = {0, INF};
  res_T res = RES_OK;
  ASSERT(htrdr && pos && dir && val);

  ctx.tau_sampled = ssp_ran_exp(rng, 1.0);
  ctx.sky = htrdr->sky;
  svx_tree = htrdr_sky_get_svx_tree(htrdr->sky);
  res = svx_octree_trace_ray(svx_tree, pos, dir, range, NULL,
    discard_hit, &ctx, &hit);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "error computing the transmission for the position `%g %g %g' "
      "along the direction `%g %g %g'.\n",
      SPLIT3(pos), SPLIT3(dir));
    goto error;
  }
  if(!SVX_HIT_NONE(&hit)) {
    *val = 0;
  } else {
    *val = ctx.tau_min ? exp(-ctx.tau_min) : 1.0;
  }

exit:
  return res;
error:
  goto exit;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_solve_transmission_buffer(struct htrdr* htrdr)
{
  struct htrdr_buffer_layout buf_layout;
  struct ssp_rng* rng = NULL;
  double pixsz[2];
  int64_t mcode_max;
  int64_t mcode;
  ATOMIC res = RES_OK;
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

  /* Compute the maximum morton code */
  mcode_max = (int64_t)round_up_pow2(MMAX(buf_layout.height, buf_layout.width));
  mcode_max = mcode_max*mcode_max;

  /* Setup the number of threads */
  omp_set_num_threads((int)htrdr->nthreads);

  #pragma omp parallel for schedule(dynamic, 32/*chunck size*/)
  for(mcode=0; mcode<mcode_max; ++mcode) {
    size_t ispp;
    double* val;
    double x, y;
    double accum;
    uint16_t ix, iy;
    res_T res_local = RES_OK;

    if(ATOMIC_GET(&res) != RES_OK) continue;

    ix = morton2D_decode((uint32_t)(mcode>>0));
    if(ix > buf_layout.width) continue;
    iy = morton2D_decode((uint32_t)(mcode>>1));
    if(iy > buf_layout.height) continue;

    val = htrdr_buffer_at(htrdr->buf, ix, iy);

    /* Define lower left the pixel coordinate in the normalized image space */
    x = (double)ix*pixsz[0];
    y = (double)iy*pixsz[1];

    accum = 0;
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
      res_local =  transmission_realisation(htrdr, rng, pos, htrdr->main_dir, &T);
      if(res_local != RES_OK) {
        ATOMIC_SET(&res, res_local);
        continue;
      }

      accum += T;
    }
    *val = accum / (double)htrdr->spp;
  }

exit:
  if(rng) SSP(rng_ref_put(rng));
  return (res_T)res;
error:
  goto exit;
}

