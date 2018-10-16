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
#include "htrdr_ground.h"

#include <rsys/cstr.h>
#include <rsys/float2.h>
#include <rsys/float3.h>

#include <star/s3d.h>
#include <star/s3daw.h>

struct ray_context {
  float range[2];
  struct s3d_hit hit_prev;
};
static const struct ray_context RAY_CONTEXT_NULL = {
  {0, INF}, S3D_HIT_NULL__
};

struct htrdr_ground {
  struct s3d_scene_view* view;
  float lower[3]; /* Ground lower bound */
  float upper[3]; /* Ground upper bound */
  int repeat; /* Make the ground infinite in X and Y */

  struct htrdr* htrdr;
  ref_T ref;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
/* Check that `hit' roughly lies on an edge. For triangular primitives, a
 * simple but approximative way is to test that its position have at least one
 * barycentric coordinate roughly equal to 0 or 1. */
static FINLINE int
hit_on_edge(const struct s3d_hit* hit)
{
  const float on_edge_eps = 1.e-4f;
  float w;
  ASSERT(hit && !S3D_HIT_NONE(hit));
  w = 1.f - hit->uv[0] - hit->uv[1];
  return eq_epsf(hit->uv[0], 0.f, on_edge_eps)
      || eq_epsf(hit->uv[0], 1.f, on_edge_eps)
      || eq_epsf(hit->uv[1], 0.f, on_edge_eps)
      || eq_epsf(hit->uv[1], 1.f, on_edge_eps)
      || eq_epsf(w, 0.f, on_edge_eps)
      || eq_epsf(w, 1.f, on_edge_eps);
}

static int
ground_filter
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   void* ray_data,
   void* filter_data)
{
  const struct ray_context* ray_ctx = ray_data;
  (void)ray_org, (void)ray_dir, (void)filter_data;

  if(!ray_ctx) return 0;

  if(S3D_PRIMITIVE_EQ(&hit->prim, &ray_ctx->hit_prev.prim)) return 1;

  if(!S3D_HIT_NONE(&ray_ctx->hit_prev) && eq_epsf(hit->distance, 0, 1.e-1f)) {
    /* If the targeted point is near of the origin, check that it lies on an
     * edge/vertex shared by the 2 primitives. */
    return hit_on_edge(&ray_ctx->hit_prev) && hit_on_edge(hit);
  }

  return hit->distance <= ray_ctx->range[0]
      || hit->distance >= ray_ctx->range[1];
}

static void
release_ground(ref_T* ref)
{
  struct htrdr_ground* ground;
  ASSERT(ref);
  ground = CONTAINER_OF(ref, struct htrdr_ground, ref);
  if(ground->view) S3D(scene_view_ref_put(ground->view));
  MEM_RM(ground->htrdr->allocator, ground);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_ground_create
  (struct htrdr* htrdr,
   const char* obj_filename,
   const int repeat_ground, /* Infinitely repeat the ground in X and Y */
   struct htrdr_ground** out_ground)
{
  struct s3d_scene* scn = NULL;
  struct s3daw* s3daw = NULL;
  struct htrdr_ground* ground = NULL;
  size_t ishape;
  size_t nshapes;
  res_T res = RES_OK;
  ASSERT(htrdr && obj_filename && out_ground);

  ground = MEM_CALLOC(htrdr->allocator, 1, sizeof(*ground));
  if(!ground) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "%s: could not allocate the ground data structure -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }
  ref_init(&ground->ref);
  ground->htrdr = htrdr;
  ground->repeat = repeat_ground;

  res = s3daw_create(&htrdr->logger, htrdr->allocator, NULL, NULL, htrdr->s3d,
    htrdr->verbose, &s3daw);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: could not create the Star-3DAW device -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }

  res = s3daw_load(s3daw, obj_filename);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "%s: could not load the obj file `%s' -- %s.\n",
      FUNC_NAME, obj_filename, res_to_cstr(res));
    goto error;
  }

  res = s3d_scene_create(htrdr->s3d, &scn);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: could not create the Star-3D scene of the ground -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }

  S3DAW(get_shapes_count(s3daw, &nshapes));
  FOR_EACH(ishape, 0, nshapes) {
    struct s3d_shape* shape;
    S3DAW(get_shape(s3daw, ishape, &shape));
    res = s3d_mesh_set_hit_filter_function(shape, ground_filter, NULL);
    if(res != RES_OK) {
      htrdr_log_err(htrdr,
        "%s: could not setup the hit filter function of the ground geometry "
        "-- %s.\n", FUNC_NAME, res_to_cstr(res));
      goto error;
    }
    res = s3d_scene_attach_shape(scn, shape);
    if(res != RES_OK) {
      htrdr_log_err(htrdr,
        "%s: could not attach the ground geometry to its Star-3D scene -- %s.\n",
        FUNC_NAME, res_to_cstr(res));
      goto error;
    }
  }

  res = s3d_scene_view_create(scn, S3D_TRACE, &ground->view);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: could not create the Star-3D scene view of the ground geometry "
      "-- %s.\n", FUNC_NAME, res_to_cstr(res));
    goto error;
  }

  res = s3d_scene_view_get_aabb(ground->view, ground->lower, ground->upper);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: could not get the ground bounding box -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }

exit:
  if(s3daw) S3DAW(ref_put(s3daw));
  if(scn) S3D(scene_ref_put(scn));
  *out_ground = ground;
  return res;
error:
  if(ground) {
    htrdr_ground_ref_put(ground);
    ground = NULL;
  }
  goto exit;
}

void
htrdr_ground_ref_get(struct htrdr_ground* ground)
{
  ASSERT(ground);
  ref_get(&ground->ref);
}

void
htrdr_ground_ref_put(struct htrdr_ground* ground)
{
  ASSERT(ground);
  ref_put(&ground->ref, release_ground);
}

res_T
htrdr_ground_trace_ray
  (struct htrdr_ground* ground,
   const double org[3],
   const double dir[3], /* Must be normalized */
   const double range[2],
   const struct s3d_hit* prev_hit,
   struct s3d_hit* hit)
{
  struct ray_context ctx = RAY_CONTEXT_NULL;
  float ray_org[3];
  float ray_dir[3];
  res_T res = RES_OK;
  ASSERT(ground && org && dir && range && hit);

  f3_set_d3(ray_org, org);
  f3_set_d3(ray_dir, dir);
  f2_set_d2(ctx.range, range);
  ctx.hit_prev = prev_hit ? *prev_hit : S3D_HIT_NULL;

  res = s3d_scene_view_trace_ray
    (ground->view, ray_org, ray_dir, ctx.range, &ctx, hit);
  if(res != RES_OK) {
    htrdr_log_err(ground->htrdr,
      "%s: could not trace the ray against the ground geometry -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }

exit:
  return res;
error:
  goto exit;
}

