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
#include "htrdr_solve.h"

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

