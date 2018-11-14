/* Copyright (C) 2018 CNRS, |Meso|Star>, Université Paul Sabatier
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
#include "htrdr_slab.h"

#include <rsys/cstr.h>
#include <math.h>

/*******************************************************************************
 * Local function
 ******************************************************************************/
res_T
htrdr_slab_trace_ray
  (struct htrdr* htrdr,
   const double org[3],
   const double dir[3],
   const double range[2],
   const double cell_low[2],
   const double cell_upp[2],
   htrdr_trace_cell_T trace_cell,
   const size_t max_steps,
   void* trace_cell_context)
{
  double pos[2];
  double org_cs[3]; /* Origin of the ray transformed in local cell space */
  double cell_low_ws[3]; /* Cell lower bound in world space */
  double cell_upp_ws[3]; /* Cell upper bound in world space */
  double cell_sz[3]; /* Size of a cell */
  double t_max[3], t_delta[2], t_min_z;
  size_t istep;
  int64_t xy[2]; /* 2D index of the repeated cell */
  int incr[2]; /* Index increment */
  res_T res = RES_OK;
  ASSERT(htrdr && org && dir && range && cell_low && cell_upp && trace_cell);
  ASSERT(range[0] < range[1]);

  /* Check that the ray intersects the slab */
  t_min_z  = (cell_low[2] - org[2]) / dir[2];
  t_max[2] = (cell_upp[2] - org[2]) / dir[2];
  if(t_min_z > t_max[2]) SWAP(double, t_min_z, t_max[2]);
  t_min_z  = MMAX(t_min_z,  range[0]);
  t_max[2] = MMIN(t_max[2], range[1]);
  if(t_min_z > t_max[2]) return RES_OK;

  /* Compute the size of a cell */
  cell_sz[0] = cell_upp[0] - cell_low[0];
  cell_sz[1] = cell_upp[1] - cell_low[1];
  cell_sz[2] = cell_upp[2] - cell_low[2];

  /* Define the 2D index of the current cell. (0,0) is the index of the
   * non duplicated cell */
  pos[0] = org[0] + t_min_z*dir[0];
  pos[1] = org[1] + t_min_z*dir[1];
  xy[0] = (int64_t)floor((pos[0] - cell_low[0]) / cell_sz[0]);
  xy[1] = (int64_t)floor((pos[1] - cell_low[1]) / cell_sz[1]);

  /* Define the 2D index increment wrt dir sign */
  incr[0] = dir[0] < 0 ? -1 : 1;
  incr[1] = dir[1] < 0 ? -1 : 1;

  /* Compute the world space AABB of the repeated cell currently hit */
  cell_low_ws[0] = cell_low[0] + (double)xy[0]*cell_sz[0];
  cell_low_ws[1] = cell_low[1] + (double)xy[1]*cell_sz[1];
  cell_low_ws[2] = cell_low[2];
  cell_upp_ws[0] = cell_low_ws[0] + cell_sz[0];
  cell_upp_ws[1] = cell_low_ws[1] + cell_sz[1];
  cell_upp_ws[2] = cell_upp[2];

  /* Compute the max ray intersection with the current cell */
  t_max[0] = ((dir[0]<0 ? cell_low_ws[0] : cell_upp_ws[0]) - org[0]) / dir[0];
  t_max[1] = ((dir[1]<0 ? cell_low_ws[1] : cell_upp_ws[1]) - org[1]) / dir[1];
  /*t_max[2] = ((dir[2]<0 ? cell_low_ws[2] : cell_upp_ws[2]) - org[2]) / dir[2];*/
  ASSERT(t_max[0] >= 0 && t_max[1] >= 0 && t_max[2] >= 0);

  /* Compute the distance along the ray to traverse in order to move of a
   * distance equal to the cloud size along the X and Y axis */
  t_delta[0] = (dir[0]<0 ? -cell_sz[0] : cell_sz[0]) / dir[0];
  t_delta[1] = (dir[1]<0 ? -cell_sz[1] : cell_sz[1]) / dir[1];
  ASSERT(t_delta[0] >= 0 && t_delta[1] >= 0);

  org_cs[2] = org[2];
  FOR_EACH(istep, 0, max_steps) {
    int iaxis;
    int hit;

    /* Transform the ray origin in the local cell space */
    org_cs[0] = org[0] - (double)xy[0]*cell_sz[0];
    org_cs[1] = org[1] - (double)xy[1]*cell_sz[1];

    res = trace_cell(org_cs, dir, range, trace_cell_context, &hit);
    if(res != RES_OK) {
      htrdr_log_err(htrdr,
        "%s: could not trace the ray in the repeated cells -- %s.\n",
        FUNC_NAME, res_to_cstr(res));
      goto error;
    }
    if(hit) goto exit;

    /* Define the next axis to traverse */
    iaxis = t_max[0] < t_max[1]
      ? (t_max[0] < t_max[2] ? 0 : 2)
      : (t_max[1] < t_max[2] ? 1 : 2);

    if(iaxis == 2) break; /* The ray traverse the slab */

    if(t_max[iaxis] >= range[1]) break; /* Out of bound */

    t_max[iaxis] += t_delta[iaxis];

    /* Define the 2D index of the next traversed cloud */
    xy[iaxis] += incr[iaxis];
  }

exit:
  return res;
error:
  goto exit;
}


