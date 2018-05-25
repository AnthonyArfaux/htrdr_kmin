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
#include "htrdr_rectangle.h"

#include <rsys/double3.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

struct htrdr_rectangle {
  /* Frame of the rectangle in world space */
  double axis_x[3];
  double axis_y[3];

  double position[3]; /* Center of the rectangle */
  struct htrdr* htrdr;
  ref_T ref;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
rectangle_release(ref_T* ref)
{
  struct htrdr_rectangle* rect;
  ASSERT(ref);
  rect = CONTAINER_OF(ref, struct htrdr_rectangle, ref);
  MEM_RM(rect->htrdr->allocator, rect);
}

/*******************************************************************************
 * Local fuuction
 ******************************************************************************/
res_T
htrdr_rectangle_create
  (struct htrdr* htrdr,
   const double pos[3],
   const double tgt[3],
   const double up[2],
   const double sz[2],
   struct htrdr_rectangle** out_rect)
{
  struct htrdr_rectangle* rect = NULL;
  double x[3], y[3], z[3];
  res_T res = RES_OK;
  ASSERT(htrdr && pos && tgt && up && sz && out_rect);

  rect = MEM_CALLOC(htrdr->allocator, 1, sizeof(*rect));
  if(!rect) {
    htrdr_log_err(htrdr, "could not allocate the rectangle data structure.\n");
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&rect->ref);
  rect->htrdr = htrdr;

  if(sz[0] <= 0 || sz[1] <= 0) {
    htrdr_log_err(htrdr,
      "invalid rectangle size `%g %g'. It must be strictly positive.\n",
      SPLIT2(sz));
    res = RES_BAD_ARG;
    goto error;
  }

  if(d3_normalize(z, d3_sub(z, tgt, pos)) <= 0
  || d3_normalize(x, d3_cross(x, z, up)) <= 0
  || d3_normalize(y, d3_cross(y, z, x)) <= 0) {
    htrdr_log_err(htrdr, "invalid rectangle frame:\n"
      "\tposition = %g %g %g\n"
      "\ttarget = %g %g %g\n"
      "\tup = %g %g %g\n",
      SPLIT3(pos), SPLIT3(tgt), SPLIT3(up));
    res = RES_BAD_ARG;
    goto error;
  }

  d3_muld(rect->axis_x, x, sz[0]*0.5);
  d3_muld(rect->axis_y, y, sz[1]*0.5);
  d3_set(rect->position, pos);

exit:
  *out_rect = rect;
  return res;
error:
  if(rect) {
    htrdr_rectangle_ref_put(rect);
    rect = NULL;
  }
  goto exit;
}

void
htrdr_rectangle_sample_pos
  (const struct htrdr_rectangle* rect,
   const double sample[2], /* In [0, 1[ */
   double pos[3])
{
  double x[3], y[3];
  ASSERT(rect && sample && pos);
  d3_muld(x, rect->axis_x, sample[0]*2.0 - 1.0);
  d3_muld(y, rect->axis_y, sample[1]*2.0 - 1.0);
  d3_add(pos, d3_add(pos, rect->position, x), y);
}

void
htrdr_rectangle_ref_get(struct htrdr_rectangle* rect)
{
  ASSERT(rect);
  ref_get(&rect->ref);
}

void
htrdr_rectangle_ref_put(struct htrdr_rectangle* rect)
{
  ASSERT(rect);
  ref_put(&rect->ref, rectangle_release);
}

