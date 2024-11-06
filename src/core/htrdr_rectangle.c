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

#include "core/htrdr.h"
#include "core/htrdr_log.h"
#include "core/htrdr_rectangle.h"

#include <rsys/double2.h>
#include <rsys/double3.h>
#include <rsys/double33.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

struct htrdr_rectangle {
  /* Frame of the rectangle in world space */
  double axis_x[3];
  double axis_y[3];
  double normal[3];

  double size[2];

  double local2world[12]; /* Rectangle to world transformation matrix */
  double world2local[12]; /* World to rectangle transformation matrix */

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
  struct htrdr* htrdr;
  ASSERT(ref);
  rect = CONTAINER_OF(ref, struct htrdr_rectangle, ref);
  htrdr = rect->htrdr;
  MEM_RM(htrdr_get_allocator(htrdr), rect);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Exported functions
 ******************************************************************************/
res_T
htrdr_rectangle_create
  (struct htrdr* htrdr,
   const double sz[2],
   const double pos[3],
   const double tgt[3],
   const double up[3],
   struct htrdr_rectangle** out_rect)
{
  struct htrdr_rectangle* rect = NULL;
  double x[3], y[3], z[3];
  double trans[3];
  res_T res = RES_OK;
  ASSERT(htrdr && pos && tgt && up && sz && out_rect);

  rect = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*rect));
  if(!rect) {
    htrdr_log_err(htrdr, "Could not allocate the rectangle data structure.\n");
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&rect->ref);
  htrdr_ref_get(htrdr);
  rect->htrdr = htrdr;

  if(sz[0] <= 0 || sz[1] <= 0) {
    htrdr_log_err(htrdr,
      "Invalid rectangle size `%g %g'. It must be strictly positive.\n",
      SPLIT2(sz));
    res = RES_BAD_ARG;
    goto error;
  }

  if(d3_normalize(z, d3_sub(z, tgt, pos)) <= 0
  || d3_normalize(x, d3_cross(x, z, up)) <= 0
  || d3_normalize(y, d3_cross(y, z, x)) <= 0) {
    htrdr_log_err(htrdr, "Invalid rectangle frame:\n"
      "\tposition = %g %g %g\n"
      "\ttarget = %g %g %g\n"
      "\tup = %g %g %g\n",
      SPLIT3(pos), SPLIT3(tgt), SPLIT3(up));
    res = RES_BAD_ARG;
    goto error;
  }

  /* Setup the local to world transformation matrix */
  d3_set(rect->local2world+0, x);
  d3_set(rect->local2world+3, y);
  d3_set(rect->local2world+6, z);
  d3_set(rect->local2world+9, pos);

  /* Inverse the local to world transformation matrix. Note that since the
   * represented frame is orthonormal one can transpose the original rotation
   * matrix to cheaply compute its inverse */
  d33_transpose(rect->world2local, rect->local2world);

  /* Compute the affine inverse transform */
  d3_minus(trans, pos);
  d33_muld3(rect->world2local+9, rect->world2local, trans);

  d3_muld(rect->axis_x, x, sz[0]*0.5);
  d3_muld(rect->axis_y, y, sz[1]*0.5);
  d3_set(rect->normal, z);
  d3_set(rect->position, pos);

  d2_set(rect->size, sz);

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

void
htrdr_rectangle_get_normal(const struct htrdr_rectangle* rect, double normal[3])
{
  ASSERT(rect && normal);
  d3_set(normal, rect->normal);
}

void
htrdr_rectangle_get_center(const struct htrdr_rectangle* rect, double pos[3])
{
  ASSERT(rect && pos);
  d3_set(pos, rect->position);
}

double*
htrdr_rectangle_get_transform
  (const struct htrdr_rectangle* rect,
   double transform[12])
{
  ASSERT(rect && transform);
  d3_set(transform+0, rect->local2world+0);
  d3_set(transform+3, rect->local2world+3);
  d3_set(transform+6, rect->local2world+6);
  d3_set(transform+9, rect->local2world+9);
  return transform;
}

double*
htrdr_rectangle_get_transform_inverse
  (const struct htrdr_rectangle* rect,
   double transform_inverse[12])
{
  ASSERT(rect && transform_inverse);
  d3_set(transform_inverse+0, rect->world2local+0);
  d3_set(transform_inverse+3, rect->world2local+3);
  d3_set(transform_inverse+6, rect->world2local+6);
  d3_set(transform_inverse+9, rect->world2local+9);
  return transform_inverse;
}

void
htrdr_rectangle_get_size(const struct htrdr_rectangle* rect, double size[2])
{
  ASSERT(rect && size);
  d2_set(size, rect->size);
}
