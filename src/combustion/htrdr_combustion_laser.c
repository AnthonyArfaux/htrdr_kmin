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

#include "combustion/htrdr_combustion_laser.h"

#include "core/htrdr.h"
#include "core/htrdr_log.h"
#include "core/htrdr_rectangle.h"

#include <rsys/cstr.h>
#include <rsys/double33.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

struct htrdr_combustion_laser {
  struct htrdr_rectangle* surface; /* Surface emission */
  double world2local[12];
  double low_ls[3], upp_ls[3]; /* Local space AABB */
  double flux_density; /* In W/m^2 */
  double wavelength; /* In nm */
  ref_T ref;
  struct htrdr* htrdr;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
laser_release(ref_T* ref)
{
  struct htrdr_combustion_laser* laser;
  struct htrdr* htrdr;
  ASSERT(ref);
  laser = CONTAINER_OF(ref, struct htrdr_combustion_laser, ref);
  if(laser->surface) htrdr_rectangle_ref_put(laser->surface);
  htrdr = laser->htrdr;
  MEM_RM(htrdr_get_allocator(htrdr), laser);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_combustion_laser_create
  (struct htrdr* htrdr,
   const struct htrdr_combustion_laser_create_args* args,
   struct htrdr_combustion_laser** out_laser)
{
  struct htrdr_combustion_laser* laser = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && args && out_laser);

  if(args->flux_density <= 0) {
    htrdr_log_err(htrdr, "Invalid flux density %g W.m^2\n", args->flux_density);
    res = RES_BAD_ARG;
    goto error;
  }

  if(args->wavelength <= 0) {
    htrdr_log_err(htrdr, "Invalid wavelength %g nm\n", args->wavelength);
    res = RES_BAD_ARG;
    goto error;
  }

  laser = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*laser));
  if(!laser) {
    htrdr_log_err(htrdr, "Could not allocate the laser data structure.\n");
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&laser->ref);
  htrdr_ref_get(htrdr);
  laser->htrdr = htrdr;
  laser->flux_density = args->flux_density;
  laser->wavelength = args->wavelength;

  res = htrdr_rectangle_create
    (laser->htrdr,
     args->surface.size,
     args->surface.position,
     args->surface.target,
     args->surface.up,
     &laser->surface);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "Could not create the laser surface emission -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  htrdr_rectangle_get_transform_inverse(laser->surface, laser->world2local);

  /* Define the laser sheet AABB in local space */
  laser->upp_ls[0] = args->surface.size[0] * 0.5;
  laser->upp_ls[1] = args->surface.size[1] * 0.5;
  laser->upp_ls[2] = 0;
  laser->low_ls[0] = -laser->upp_ls[0];
  laser->low_ls[1] = -laser->upp_ls[1];
  laser->upp_ls[2] = INF;

exit:
  *out_laser = laser;
  return res;
error:
  if(laser) { htrdr_combustion_laser_ref_put(laser); laser = NULL; }
  goto exit;
}

void
htrdr_combustion_laser_ref_get(struct htrdr_combustion_laser* laser)
{
  ASSERT(laser);
  ref_get(&laser->ref);
}

void
htrdr_combustion_laser_ref_put(struct htrdr_combustion_laser* laser)
{
  ASSERT(laser);
  ref_put(&laser->ref, laser_release);
}

void
htrdr_combustion_laser_trace_ray
  (struct htrdr_combustion_laser* laser,
   const double pos[3],
   const double dir[3],
   const double range[2],
   double distance[2])
{
  double pos_ls[3]; /* Position in local space */
  double dir_ls[3]; /* Direction in local space */
  double t_min;
  double t_max;
  int i;
  ASSERT(laser && pos && dir && range && distance);

  /* Transform the ray in local space */
  d33_muld3(pos_ls, laser->world2local, pos);
  d33_muld3(dir_ls, laser->world2local, dir);
  d3_add(pos_ls, pos_ls, laser->world2local+9);

  /* Reset the distance */
  distance[0] = INF;
  distance[1] = INF;

  /* Initialise the range */
  t_min = range[0];
  t_max = range[1];

  /* TODO one can optimise the Ray/AABB intersection test along the Z axis
   * whose AABB interval is in [0, inf) */
  FOR_EACH(i, 0, 3) {
    const double rcp_dir_ls = 1.0/dir_ls[i];
    double t0 = (laser->low_ls[i] - pos_ls[i]) * rcp_dir_ls;
    double t1 = (laser->upp_ls[i] - pos_ls[i]) * rcp_dir_ls;
    if(t0 > t1) SWAP(double, t0, t1);
    t_min = MMAX(t_min, t0);
    t_max = MMIN(t_max, t1);
    if(t_min > t_max) return; /* No intersection */
  }

  /* Save the hit distance */
  distance[0] = t_min;
  distance[1] = t_max;
}

void
htrdr_combustion_laser_get_mesh
  (const struct htrdr_combustion_laser* laser,
   const double extent,
   struct htrdr_combustion_laser_mesh* mesh)
{
  double transform[12];
  double low[3];
  double upp[3];
  int i;
  ASSERT(laser && extent > 0 && mesh);

  htrdr_rectangle_get_transform(laser->surface, transform);

  /* Define the laser sheet vertices in local space */
  low[0] = laser->low_ls[0];
  low[1] = laser->low_ls[1];
  low[2] = laser->low_ls[2];
  upp[0] = laser->upp_ls[0];
  upp[1] = laser->upp_ls[1];
  upp[2] = laser->low_ls[2] + extent;

  /* Define the mesh of the laser sheet in local space */
  d3(mesh->vertices + 0*3, low[0], low[1], low[2]);
  d3(mesh->vertices + 1*3, upp[0], low[1], low[2]);
  d3(mesh->vertices + 2*3, low[0], upp[1], low[2]);
  d3(mesh->vertices + 3*3, upp[0], upp[1], low[2]);
  d3(mesh->vertices + 4*3, low[0], low[1], upp[2]);
  d3(mesh->vertices + 5*3, upp[0], low[1], upp[2]);
  d3(mesh->vertices + 6*3, low[0], upp[1], upp[2]);
  d3(mesh->vertices + 7*3, upp[0], upp[1], upp[2]);
  mesh->nvertices = 8;

  /* Transform the laser sheet vertices in world space */
  FOR_EACH(i, 0, 8) {
    double* v = mesh->vertices + i*3;
    d33_muld3(v, transform, v);
    d3_add(v, v, transform+9);
  }

  /* Define the laser sheet triangles */
  mesh->triangles[0]  = 0; mesh->triangles[1]  = 2; mesh->triangles[2]  = 1;
  mesh->triangles[3]  = 1; mesh->triangles[4]  = 2; mesh->triangles[5]  = 3;
  mesh->triangles[6]  = 0; mesh->triangles[7]  = 4; mesh->triangles[8]  = 2;
  mesh->triangles[9]  = 2; mesh->triangles[10] = 4; mesh->triangles[11] = 6;
  mesh->triangles[12] = 1; mesh->triangles[13] = 3; mesh->triangles[14] = 5;
  mesh->triangles[15] = 5; mesh->triangles[16] = 3; mesh->triangles[17] = 7;
  mesh->triangles[18] = 1; mesh->triangles[19] = 5; mesh->triangles[20] = 0;
  mesh->triangles[21] = 0; mesh->triangles[22] = 5; mesh->triangles[23] = 4;
  mesh->triangles[24] = 3; mesh->triangles[25] = 2; mesh->triangles[26] = 7;
  mesh->triangles[27] = 7; mesh->triangles[28] = 2; mesh->triangles[29] = 6;
  mesh->ntriangles = 10;
}

void
htrdr_combustion_laser_get_position
  (const struct htrdr_combustion_laser* laser,
   double pos[3])
{
  ASSERT(laser && pos);
  htrdr_rectangle_get_center(laser->surface, pos);
}

void
htrdr_combustion_laser_get_direction
  (const struct htrdr_combustion_laser* laser,
   double dir[3])
{
  ASSERT(laser && dir);
  htrdr_rectangle_get_normal(laser->surface, dir);
}

double
htrdr_combustion_laser_get_flux_density
  (const struct htrdr_combustion_laser* laser)
{
  ASSERT(laser);
  return laser->flux_density;
}

double
htrdr_combustion_laser_get_wavelength
  (const struct htrdr_combustion_laser* laser)
{
  ASSERT(laser);
  return laser->wavelength;
}

double
htrdr_combustion_laser_compute_surface_plane_distance
  (const struct htrdr_combustion_laser* laser,
   const double pos[3])
{
  double center[3];
  double normal[3];
  double vec[3];
  ASSERT(laser && pos);
  htrdr_rectangle_get_normal(laser->surface, normal);
  htrdr_rectangle_get_center(laser->surface, center);
  d3_sub(vec, pos, center);
  return d3_dot(normal, vec);
}
