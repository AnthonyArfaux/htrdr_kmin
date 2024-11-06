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


#include "atmosphere/htrdr_atmosphere_ground.h"

#include "core/htrdr_log.h"
#include "core/htrdr_geometry.h"
#include "core/htrdr_slab.h"

#include <star/s3d.h>

#include <rsys/cstr.h>
#include <rsys/double2.h>
#include <rsys/double3.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

struct trace_slab_context {
  struct htrdr_geometry* geom;
  const struct s3d_hit* hit_prev;
  struct s3d_hit* hit;
};
#define TRACE_SLAB_CONTEXT_NULL__ { NULL, NULL, NULL }
static const struct trace_slab_context TRACE_SLAB_CONTEXT_NULL =
    TRACE_SLAB_CONTEXT_NULL__;

struct htrdr_atmosphere_ground {
  struct htrdr_geometry* geom;
  int repeat; /* Make the ground infinite in X and Y */
  struct htrdr* htrdr;
  ref_T ref;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static INLINE res_T
trace_slab
  (const double org[3],
   const double dir[3],
   const double range[2],
   void* context,
   int* hit)
{
  struct trace_slab_context* ctx = context;
  struct htrdr_geometry_trace_ray_args rt_args =
    HTRDR_GEOMETRY_TRACE_RAY_ARGS_NULL;
  res_T res = RES_OK;
  ASSERT(org && dir && range && context && hit);

  d3_set(rt_args.ray_org, org);
  d3_set(rt_args.ray_dir, dir);
  d2_set(rt_args.ray_range, range);
  rt_args.hit_from = ctx->hit_prev ? *ctx->hit_prev : S3D_HIT_NULL;
  res = htrdr_geometry_trace_ray(ctx->geom, &rt_args, ctx->hit);
  if(res != RES_OK) return res;

  *hit = !S3D_HIT_NONE(ctx->hit);
  return RES_OK;
}
static void
release_ground(ref_T* ref)
{
  struct htrdr_atmosphere_ground* ground;
  struct htrdr* htrdr;
  ASSERT(ref);
  ground = CONTAINER_OF(ref, struct htrdr_atmosphere_ground, ref);
  if(ground->geom) htrdr_geometry_ref_put(ground->geom);
  htrdr = ground->htrdr;
  MEM_RM(htrdr_get_allocator(ground->htrdr), ground);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_atmosphere_ground_create
  (struct htrdr* htrdr,
   const char* obj_filename, /* May be NULL */
   struct htrdr_materials* mats, /* May be NULL */
   const int repeat_ground, /* Infinitely repeat the ground in X and Y */
   struct htrdr_atmosphere_ground** out_ground)
{
  struct htrdr_atmosphere_ground* ground = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && out_ground);

  ground = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*ground));
  if(!ground) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "%s: could not allocate the ground data structure -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }
  ref_init(&ground->ref);
  ground->repeat = repeat_ground;
  htrdr_ref_get(htrdr);
  ground->htrdr = htrdr;

  if(!obj_filename) goto exit;

  if(!mats) {
    htrdr_log_err(htrdr, "%s: missing materials.\n", FUNC_NAME);
    res = RES_BAD_ARG;
    goto error;
  }

  res = htrdr_geometry_create(ground->htrdr, obj_filename, mats, &ground->geom);
  if(res != RES_OK) goto error;

exit:
  *out_ground = ground;
  return res;
error:
  if(ground) {
    htrdr_atmosphere_ground_ref_put(ground);
    ground = NULL;
  }
  goto exit;
}

void
htrdr_atmosphere_ground_ref_get(struct htrdr_atmosphere_ground* ground)
{
  ASSERT(ground);
  ref_get(&ground->ref);
}

void
htrdr_atmosphere_ground_ref_put(struct htrdr_atmosphere_ground* ground)
{
  ASSERT(ground);
  ref_put(&ground->ref, release_ground);
}

void
htrdr_atmosphere_ground_get_interface
  (struct htrdr_atmosphere_ground* ground,
   const struct s3d_hit* hit,
   struct htrdr_interface* out_interface)
{
  /* Simply wrap the geometry function */
  htrdr_geometry_get_interface(ground->geom, hit, out_interface);
}

res_T
htrdr_atmosphere_ground_trace_ray
  (struct htrdr_atmosphere_ground* ground,
   const double org[3],
   const double dir[3], /* Must be normalized */
   const double range[2],
   const struct s3d_hit* prev_hit,
   struct s3d_hit* hit)
{
  res_T res = RES_OK;
  ASSERT(ground && org && dir && range && hit);

  if(!ground->geom) { /* No ground geometry */
    *hit = S3D_HIT_NULL;
    goto exit;
  }

  if(!ground->repeat) {
    struct htrdr_geometry_trace_ray_args rt_args =
      HTRDR_GEOMETRY_TRACE_RAY_ARGS_NULL;

    d3_set(rt_args.ray_org, org);
    d3_set(rt_args.ray_dir, dir);
    d2_set(rt_args.ray_range, range);
    if(prev_hit) rt_args.hit_from = *prev_hit;
    res = htrdr_geometry_trace_ray(ground->geom, &rt_args, hit);
    if(res != RES_OK) goto error;

  } else {
    struct trace_slab_context slab_ctx = TRACE_SLAB_CONTEXT_NULL;
    double low[3], upp[3];

    htrdr_geometry_get_aabb(ground->geom, low, upp);

    *hit = S3D_HIT_NULL;
    slab_ctx.geom = ground->geom;
    slab_ctx.hit_prev = prev_hit;
    slab_ctx.hit = hit;

    res = htrdr_slab_trace_ray(ground->htrdr, org, dir, range, low, upp,
      trace_slab, 32, &slab_ctx);
    if(res != RES_OK) goto error;
  }

exit:
  return res;
error:
  goto exit;
}

res_T
htrdr_atmosphere_ground_find_closest_point
  (struct htrdr_atmosphere_ground* ground,
   const double pos[3],
   const double radius,
   struct s3d_hit* hit)
{
  double query_pos[3];
  double query_radius;
  res_T res = RES_OK;
  ASSERT(ground && pos && hit);

  if(!ground->geom) { /* No ground geometry */
    *hit = S3D_HIT_NULL;
    goto exit;
  }

  query_radius = radius;
  d3_set(query_pos, pos);

  if(ground->repeat) {
    double ground_low[3];
    double ground_upp[3];
    double ground_sz[3];
    double translation[2] = {0, 0};
    int64_t xy[2];

    /* Define the size of the ground geometry pattern */
    htrdr_geometry_get_aabb(ground->geom, ground_low, ground_upp);
    ground_sz[0] = ground_upp[0] - ground_low[0];
    ground_sz[1] = ground_upp[1] - ground_low[1];
    ground_sz[2] = ground_upp[2] - ground_low[2];

    /* Define the 2D index of the current ground instance. (0,0) is the index
     * of the non instantiated ground */
    xy[0] = (int64_t)floor((query_pos[0] - ground_low[0]) / ground_sz[0]);
    xy[1] = (int64_t)floor((query_pos[1] - ground_low[1]) / ground_sz[1]);

    /* Define the translation along the X and Y axis from world space to local
     * ground geometry space */
    translation[0] = -(double)xy[0] * ground_sz[0];
    translation[1] = -(double)xy[1] * ground_sz[1];

    /* Transform the query pos in local ground geometry space */
    query_pos[0] += translation[0];
    query_pos[1] += translation[1];
  }

  res = htrdr_geometry_find_closest_point
    (ground->geom, query_pos, query_radius, hit);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  goto exit;
}
