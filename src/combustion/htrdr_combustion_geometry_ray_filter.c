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

#include "combustion/htrdr_combustion_geometry_ray_filter.h"

#include "core/htrdr_interface.h"
#include "core/htrdr_geometry.h"

#include <rsys/float3.h>

#include <string.h>

/*******************************************************************************
 * Local functions
 ******************************************************************************/
int
geometry_ray_filter_discard_medium_interface
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   const float ray_range[2],
   void* ray_data,
   void* filter_data)
{
  struct geometry_ray_filter_context* ctx = ray_data;
  struct htrdr_interface interf = HTRDR_INTERFACE_NULL;
  const struct htrdr_mtl* mtl = NULL;
  float N[3];
  int hit_front;
  int discard = 0;
  ASSERT(hit && ray_org && ray_dir && ray_range && ray_data && !S3D_HIT_NONE(hit));
  (void)ray_org, (void)ray_dir, (void)ray_range, (void)filter_data;

  /* Recover the interface of the intersected surface */
  htrdr_geometry_get_interface(ctx->geom, hit, &interf);

  /* Define if the ray intersects the front face of the surface */
  f3_normalize(N, hit->normal); /* Limit the numerical instabilities */
  hit_front = f3_dot(ray_dir, N) < 0;

  /* Recover the material in which the ray is traced */
  if(hit_front) {
    mtl = &interf.mtl_front;
  } else {
    mtl = &interf.mtl_back;
  }

  /* The material should be semi-transparent and thus could not have a BRDF */
  ASSERT(mtl->mrumtl == NULL);

  /* Discard the intersection if the ray comes from the material to filter */
  discard = !strcmp(mtl->name, ctx->medium_name);
  return discard;
}

