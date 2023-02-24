/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2023 Observatoire de Paris
 * Copyright (C) 2022-2023 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2023 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2023 Université Paul Sabatier
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

#include <rsys/double2.h>
#include <rsys/double3.h>

static void
dump_obj(const struct htrdr_combustion_laser_mesh* mesh, FILE* stream)
{
  unsigned i;
  ASSERT(mesh && stream);

  FOR_EACH(i, 0, mesh->nvertices) {
    fprintf(stream, "v %g %g %g\n",
      mesh->vertices[i*3+0],
      mesh->vertices[i*3+1],
      mesh->vertices[i*3+2]);
  }
  FOR_EACH(i, 0, mesh->ntriangles) {
    fprintf(stream, "f %u %u %u\n",
      mesh->triangles[i*3+0]+1,
      mesh->triangles[i*3+1]+1,
      mesh->triangles[i*3+2]+1);
  }
}

int
main(int argc, char** argv)
{
  struct htrdr_args args = HTRDR_ARGS_DEFAULT;
  struct htrdr_combustion_laser_mesh mesh;
  struct htrdr_combustion_laser_create_args laser_args =
    HTRDR_COMBUSTION_LASER_CREATE_ARGS_DEFAULT;
  struct htrdr* htrdr = NULL;
  struct htrdr_combustion_laser* laser = NULL;
  double org[3];
  double dir[3];
  double range[2];
  double hit_range[2];
  double t[2];
  double pt[3];
  double x[3], y[3], z[3];
  double plane0[4];
  double plane1[4];
  FILE* fp = NULL;

  args.verbose = 1;
  htrdr_mpi_init(argc, argv);
  CHK(htrdr_create(NULL, &args, &htrdr) == RES_OK);

  /* Setup the laser sheet */
  d3(laser_args.surface.position, 0, 0, 0);
  d3(laser_args.surface.target, 10, 10, 0);
  d3(laser_args.surface.up, 0, 1, 0);
  d2(laser_args.surface.size, 100, 50);
  laser_args.wavelength = 300;
  laser_args.flux_density = 1;
  CHK(htrdr_combustion_laser_create(htrdr, &laser_args, &laser) == RES_OK);
  htrdr_combustion_laser_get_mesh(laser, 100/*arbitrary extend*/, &mesh);

  /* Write the laser geometry */
  CHK(fp = fopen("laser.obj", "w"));
  dump_obj(&mesh, fp);
  fclose(fp);

  /* Compute the frame of the surface emission */
  d3_sub(z, laser_args.surface.target, laser_args.surface.position);
  d3_cross(x, z, laser_args.surface.up);
  d3_cross(y, x, z);
  CHK(d3_normalize(y, y) != 0);

  /* Compute the bottom plane equation of the laser sheet */
  pt[0] = laser_args.surface.position[0] - y[0]*laser_args.surface.size[1]*0.5;
  pt[1] = laser_args.surface.position[1] - y[1]*laser_args.surface.size[1]*0.5;
  pt[2] = laser_args.surface.position[2] - y[2]*laser_args.surface.size[1]*0.5;
  plane0[0] = y[0];
  plane0[1] = y[1];
  plane0[2] = y[2];
  plane0[3] = -d3_dot(y, pt);

  /* Compute the top plane equation of the laser sheet */
  pt[0] = laser_args.surface.position[0] + y[0]*laser_args.surface.size[1]*0.5;
  pt[1] = laser_args.surface.position[1] + y[1]*laser_args.surface.size[1]*0.5;
  pt[2] = laser_args.surface.position[2] + y[2]*laser_args.surface.size[1]*0.5;
  plane1[0] = y[0];
  plane1[1] = y[1];
  plane1[2] = y[2];
  plane1[3] = -d3_dot(y, pt);

  /* Trace a ray that misses the laser sheet */
  d3(org, 50, 0, 0);
  d3(dir, 0, -1, 0);
  d2(range, 0, INF);
  htrdr_combustion_laser_trace_ray(laser, org, dir, range, hit_range);
  CHK(hit_range[0] > DBL_MAX);
  CHK(hit_range[1] > DBL_MAX);

  /* Trace a ray that intersects both bottom and top laser planes */
  d3(dir, 0, 1, 0);
  htrdr_combustion_laser_trace_ray(laser, org, dir, range, hit_range);
  CHK(hit_range[0] < hit_range[1]);

  /* Compute the intersection of the ray with the bottom/top laser planes */
  t[0] = (-d3_dot(plane0, org) - plane0[3]) / d3_dot(plane0, dir);
  t[1] = (-d3_dot(plane1, org) - plane1[3]) / d3_dot(plane1, dir);

  /* Check the returned distances against the computed ones */
  CHK(eq_eps(hit_range[0], t[0], 1.e-6*hit_range[0]));
  CHK(eq_eps(hit_range[1], t[1], 1.e-6*hit_range[1]));

  /* Trace a ray that starts into the laser sheet */
  range[0] = 0.5*(hit_range[0] + hit_range[1]);
  htrdr_combustion_laser_trace_ray(laser, org, dir, range, hit_range);
  CHK(hit_range[0] < hit_range[1]);
  CHK(hit_range[0] == range[0]);
  CHK(eq_eps(hit_range[1], t[1], 1.e-6*hit_range[1]));

  htrdr_ref_put(htrdr);
  htrdr_combustion_laser_ref_put(laser);
  htrdr_mpi_finalize();

  return 0;
}

