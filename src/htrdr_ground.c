/* Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
 * Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
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
#include "htrdr_interface.h"
#include "htrdr_ground.h"
#include "htrdr_mtl.h"
#include "htrdr_slab.h"

#include <aw.h>
#include <rsys/clock_time.h>
#include <rsys/cstr.h>
#include <rsys/dynamic_array_double.h>
#include <rsys/dynamic_array_size_t.h>
#include <rsys/double2.h>
#include <rsys/double3.h>
#include <rsys/float2.h>
#include <rsys/float3.h>
#include <rsys/hash_table.h>

#include <star/s3d.h>
#include <star/smtl.h>

/* Define the hash table that maps an Obj vertex id to its position into the
 * vertex buffer */
#define HTABLE_NAME vertex
#define HTABLE_KEY size_t /* Obj vertex id */
#define HTABLE_DATA size_t
#include <rsys/hash_table.h>

/* Define the hash table that maps the Star-3D shape id to its interface */
#define HTABLE_NAME interface
#define HTABLE_KEY unsigned /* Star-3D shape id */
#define HTABLE_DATA struct htrdr_interface
#include <rsys/hash_table.h>

struct mesh {
  const struct darray_double* positions;
  const struct darray_size_t* indices;
};
static const struct mesh MESH_NULL;

struct ray_context {
  float range[2];
  struct s3d_hit hit_prev;
  int id[2];
};
#define RAY_CONTEXT_NULL__ {{0,INF}, S3D_HIT_NULL__, {0,0}}
static const struct ray_context RAY_CONTEXT_NULL = RAY_CONTEXT_NULL__;

struct trace_ground_context {
  struct s3d_scene_view* view;
  struct ray_context context;
  struct s3d_hit* hit;
};
static const struct trace_ground_context TRACE_GROUND_CONTEXT_NULL = {
  NULL, RAY_CONTEXT_NULL__, NULL
};

struct htrdr_ground {
  struct s3d_scene_view* view;
  float lower[3]; /* Ground lower bound */
  float upper[3]; /* Ground upper bound */
  int repeat; /* Make the ground infinite in X and Y */

  struct htable_interface interfaces; /* Map a Star3D shape to its interface */

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

static INLINE res_T
trace_ground
  (const double org[3],
   const double dir[3],
   const double range[2],
   void* context,
   int* hit)
{
  struct trace_ground_context* ctx = context;
  float ray_org[3];
  float ray_dir[3];
  float ray_range[2];
  res_T res = RES_OK;
  ASSERT(org && dir && range && context && hit);

  f3_set_d3(ray_org, org);
  f3_set_d3(ray_dir, dir);
  f2_set_d2(ray_range, range);

  res = s3d_scene_view_trace_ray
    (ctx->view, ray_org, ray_dir, ray_range, &ctx->context, ctx->hit);
  if(res != RES_OK) return res;

  *hit = !S3D_HIT_NONE(ctx->hit);
  return RES_OK;
}

static res_T
parse_shape_interface
  (struct htrdr* htrdr,
   const char* name,
   struct htrdr_interface* interf)
{
  struct str str;
  char* mtl_name_front = NULL;
  char* mtl_name_back = NULL;
  char* tk_ctx = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && name && interf);

  str_init(htrdr->allocator, &str);

  /* Locally copy the string to parse */
  res = str_set(&str, name);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "Could not locally copy the shape material string `%s' -- %s.\n",
      name, res_to_cstr(res));
    goto error;
  }

  /* Parse the name of the front/back faces */
  mtl_name_front = strtok_r(str_get(&str), ":", &tk_ctx);
  ASSERT(mtl_name_front); /* This can't be NULL */
  mtl_name_back = strtok_r(NULL, ":", &tk_ctx);
  if(!mtl_name_back) {
    htrdr_log_err(htrdr,
      "The material name of the shape back faces are missing `%s'.\n", name);
    res = RES_BAD_ARG;
    goto error;
  }

  /* Fetch the front/back materials */
  interf->mtl_front = htrdr_mtl_get(htrdr->mtl, mtl_name_front);
  interf->mtl_back = htrdr_mtl_get(htrdr->mtl, mtl_name_back);
  if(!interf->mtl_front && !interf->mtl_back) {
    htrdr_log_err(htrdr,
      "Invalid interface `%s:%s'. "
      "The front and the back materials are both uknown.\n",
      mtl_name_front, mtl_name_back);
    res = RES_BAD_ARG;
    goto error;
  }
exit:
  str_release(&str);
  return res;
error:
  interf->mtl_front = interf->mtl_back = NULL;
  goto exit;
}

static res_T
setup_mesh
  (struct htrdr* htrdr,
   const char* filename,
   struct aw_obj* obj,
   struct aw_obj_named_group* mtl,
   struct darray_double* positions,
   struct darray_size_t* indices,
   struct htable_vertex* vertices) /* Scratch data structure */
{
  size_t iface;
  res_T res = RES_OK;
  ASSERT(htrdr && filename && obj && mtl && positions && indices && vertices);

  darray_double_clear(positions);
  darray_size_t_clear(indices);
  htable_vertex_clear(vertices);

  FOR_EACH(iface, mtl->face_id, mtl->faces_count) {
    struct aw_obj_face face;
    size_t ivertex;

    AW(obj_get_face(obj, iface, &face));
    if(face.vertices_count != 3) {
      htrdr_log_err(htrdr,
        "The obj `%s' has non-triangulated polygons "
        "while currently only triangular meshes are supported.\n",
        filename);
      res = RES_BAD_ARG;
      goto error;
    }

    FOR_EACH(ivertex, 0, face.vertices_count) {
      struct aw_obj_vertex vertex;
      size_t id;
      size_t* pid;

      AW(obj_get_vertex(obj, face.vertex_id + ivertex, &vertex));
      pid = htable_vertex_find(vertices, &vertex.position_id);
      if(pid) {
        id = *pid;
      } else {
        struct aw_obj_vertex_data vdata;

        id = darray_double_size_get(positions) / 3;

        res = darray_double_resize(positions, id*3 + 3);
        if(res != RES_OK) {
          htrdr_log_err(htrdr,
            "Could not locally copy the vertex position -- %s.\n",
            res_to_cstr(res));
          goto error;
        }

        AW(obj_get_vertex_data(obj, &vertex, &vdata));
        darray_double_data_get(positions)[id*3+0] = vdata.position[0];
        darray_double_data_get(positions)[id*3+1] = vdata.position[1];
        darray_double_data_get(positions)[id*3+2] = vdata.position[2];

        res = htable_vertex_set(vertices, &vertex.position_id, &id);
        if(res != RES_OK) {
          htrdr_log_err(htrdr,
            "Could not register the vertex position -- %s.\n",
             res_to_cstr(res));
          goto error;
        }
      }

      res = darray_size_t_push_back(indices, &id);
      if(res != RES_OK) {
        htrdr_log_err(htrdr,
          "Could not locally copy the face index -- %s\n",
          res_to_cstr(res));
        goto error;
      }
    }
  }
exit:
  return res;
error:
  darray_double_clear(positions);
  darray_size_t_clear(indices);
  htable_vertex_clear(vertices);
  goto exit;
}

static void
get_position(const unsigned ivert, float position[3], void* ctx)
{
  const struct mesh* mesh = ctx;
  const double* pos = NULL;
  ASSERT(mesh);
  ASSERT(ivert < darray_double_size_get(mesh->positions) / 3);
  pos = darray_double_cdata_get(mesh->positions) + ivert*3;
  position[0] = (float)pos[0];
  position[1] = (float)pos[1];
  position[2] = (float)pos[2];
}

static void
get_indices(const unsigned itri, unsigned indices[3], void* ctx)
{
  const struct mesh* mesh = ctx;
  const size_t* ids = NULL;
  ASSERT(mesh);
  ASSERT(itri < darray_size_t_size_get(mesh->indices) / 3);
  ids = darray_size_t_cdata_get(mesh->indices) + itri*3;
  indices[0] = (unsigned)ids[0];
  indices[1] = (unsigned)ids[1];
  indices[2] = (unsigned)ids[2];
}

static res_T
create_s3d_shape
  (struct htrdr* htrdr,
   const struct darray_double* positions,
   const struct darray_size_t* indices,
   struct s3d_shape** out_shape)
{
  struct s3d_shape* shape = NULL;
  struct s3d_vertex_data vdata = S3D_VERTEX_DATA_NULL;
  struct mesh mesh = MESH_NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && positions && indices && out_shape);
  ASSERT(darray_double_size_get(positions) != 0);
  ASSERT(darray_size_t_size_get(indices) != 0);
  ASSERT(darray_double_size_get(positions)%3 == 0);
  ASSERT(darray_size_t_size_get(indices)%3 == 0);

  mesh.positions = positions;
  mesh.indices = indices;

  res = s3d_shape_create_mesh(htrdr->s3d, &shape);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "Error creating a Star-3D shape -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  vdata.usage = S3D_POSITION;
  vdata.type = S3D_FLOAT3;
  vdata.get = get_position;

  res = s3d_mesh_setup_indexed_vertices
    (shape, (unsigned int)(darray_size_t_size_get(indices)/3), get_indices,
     (unsigned int)(darray_double_size_get(positions)/3), &vdata, 1, &mesh);
  if(res != RES_OK){
    htrdr_log_err(htrdr, "Could not setup the Star-3D shape -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  res = s3d_mesh_set_hit_filter_function(shape, ground_filter, NULL);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "Could not setup the Star-3D hit filter function of the ground geometry "
       "-- %s.\n", res_to_cstr(res));
    goto error;
  }

exit:
  *out_shape = shape;
  return res;
error:
  if(shape) {
    S3D(shape_ref_put(shape));
    shape = NULL;
  }
  goto exit;
}

static res_T
setup_ground(struct htrdr_ground* ground, const char* obj_filename)
{
  struct aw_obj_desc desc;
  struct htable_vertex vertices;
  struct darray_double positions;
  struct darray_size_t indices;
  struct aw_obj* obj = NULL;
  struct s3d_shape* shape = NULL;
  struct s3d_scene* scene = NULL;
  size_t iusemtl;

  res_T res = RES_OK;
  ASSERT(obj_filename);

  htable_vertex_init(ground->htrdr->allocator, &vertices);
  darray_double_init(ground->htrdr->allocator, &positions);
  darray_size_t_init(ground->htrdr->allocator, &indices);

  res = aw_obj_create(&ground->htrdr->logger, ground->htrdr->allocator,
    ground->htrdr->verbose, &obj);
  if(res != RES_OK) {
    htrdr_log_err(ground->htrdr, "Could not create the obj loader -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  res = s3d_scene_create(ground->htrdr->s3d, &scene);
  if(res != RES_OK) {
    htrdr_log_err(ground->htrdr, "Could not create the Star-3D scene -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  /* Load the geometry data */
  res = aw_obj_load(obj, obj_filename);
  if(res != RES_OK) goto error;

  /* Fetch the descriptor of the loaded geometry */
  AW(obj_get_desc(obj, &desc));

  if(desc.usemtls_count == 0) {
    htrdr_log_err(ground->htrdr, "The obj `%s' has no material.\n", obj_filename);
    res = RES_BAD_ARG;
    goto error;
  }

  /* Setup the geometry */
  FOR_EACH(iusemtl, 0, desc.usemtls_count) {
    struct aw_obj_named_group mtl;
    struct htrdr_interface interf;
    unsigned ishape;

    AW(obj_get_mtl(obj, iusemtl , &mtl));

    res = parse_shape_interface(ground->htrdr, mtl.name, &interf);
    if(res != RES_OK) goto error;

    res = setup_mesh
      (ground->htrdr, obj_filename, obj, &mtl, &positions, &indices, &vertices);
    if(res != RES_OK) goto error;

    res = create_s3d_shape(ground->htrdr, &positions, &indices, &shape);
    if(res != RES_OK) goto error;

    S3D(shape_get_id(shape, &ishape));
    res = htable_interface_set(&ground->interfaces, &ishape, &interf);
    if(res != RES_OK) {
      htrdr_log_err(ground->htrdr,
        "Could not map the Star-3D shape to its Star-Materials -- %s.\n",
        res_to_cstr(res));
      goto error;
    }

    res = s3d_scene_attach_shape(scene, shape);
    if(res != RES_OK) {
      htrdr_log_err(ground->htrdr,
        "Could not attach a Star-3D shape to the Star-3D scene -- %s.\n",
        res_to_cstr(res));
      goto error;
    }

    S3D(shape_ref_put(shape));
    shape = NULL;
  }

  res = s3d_scene_view_create(scene, S3D_GET_PRIMITIVE|S3D_TRACE, &ground->view);
  if(res != RES_OK) {
    htrdr_log_err(ground->htrdr,
      "Could not create the Star-3D scene view -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  res = s3d_scene_view_get_aabb(ground->view, ground->lower, ground->upper);
  if(res != RES_OK) {
    htrdr_log_err(ground->htrdr,
      "Could not get the bounding box of the geometry -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

exit:
  if(obj) AW(obj_ref_put(obj));
  if(shape) S3D(shape_ref_put(shape));
  if(scene) S3D(scene_ref_put(scene));
  htable_vertex_release(&vertices);
  darray_double_release(&positions);
  darray_size_t_release(&indices);
  return res;
error:
  goto exit;
}

static void
release_ground(ref_T* ref)
{
  struct htrdr_ground* ground;
  ASSERT(ref);
  ground = CONTAINER_OF(ref, struct htrdr_ground, ref);
  if(ground->view) S3D(scene_view_ref_put(ground->view));
  htable_interface_release(&ground->interfaces);
  MEM_RM(ground->htrdr->allocator, ground);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_ground_create
  (struct htrdr* htrdr,
   const char* obj_filename, /* May be NULL */
   const int repeat_ground, /* Infinitely repeat the ground in X and Y */
   struct htrdr_ground** out_ground)
{
  char buf[128];
  struct htrdr_ground* ground = NULL;
  struct time t0, t1;
  res_T res = RES_OK;
  ASSERT(htrdr && out_ground);

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
  f3_splat(ground->lower, (float)INF);
  f3_splat(ground->upper,-(float)INF);
  htable_interface_init(ground->htrdr->allocator, &ground->interfaces);

  if(!obj_filename) goto exit;

  time_current(&t0);
  res = setup_ground(ground, obj_filename);
  if(res != RES_OK) goto error;
  time_sub(&t0, time_current(&t1), &t0);
  time_dump(&t0, TIME_ALL, NULL, buf, sizeof(buf));
  htrdr_log(ground->htrdr, "Setup ground in %s\n", buf);

exit:
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

void
htrdr_ground_get_interface
  (struct htrdr_ground* ground,
   const struct s3d_hit* hit,
   struct htrdr_interface* out_interface)
{
  struct htrdr_interface* interf = NULL;
  ASSERT(ground && hit && out_interface);

  interf = htable_interface_find(&ground->interfaces, &hit->prim.geom_id);
  ASSERT(interf);

  *out_interface = *interf;
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
  res_T res = RES_OK;
  ASSERT(ground && org && dir && range && hit);

  if(!ground->view) { /* No ground geometry */
    *hit = S3D_HIT_NULL;
    goto exit;
  }

  if(!ground->repeat) {
    struct ray_context ray_ctx = RAY_CONTEXT_NULL;
    float ray_org[3];
    float ray_dir[3];

    f3_set_d3(ray_org, org);
    f3_set_d3(ray_dir, dir);
    f2_set_d2(ray_ctx.range, range);
    ray_ctx.hit_prev = prev_hit ? *prev_hit : S3D_HIT_NULL;

    res = s3d_scene_view_trace_ray
      (ground->view, ray_org, ray_dir, ray_ctx.range, &ray_ctx, hit);
    if(res != RES_OK) {
      htrdr_log_err(ground->htrdr,
        "%s: could not trace the ray against the ground geometry -- %s.\n",
        FUNC_NAME, res_to_cstr(res));
      goto error;
    }
  } else {
    struct trace_ground_context slab_ctx = TRACE_GROUND_CONTEXT_NULL;
    double low[3], upp[3];

    *hit = S3D_HIT_NULL;
    slab_ctx.view = ground->view;
    slab_ctx.context.range[0] = (float)range[0];
    slab_ctx.context.range[1] = (float)range[1];
    slab_ctx.context.hit_prev = prev_hit ? *prev_hit : S3D_HIT_NULL;
    slab_ctx.hit = hit;

    d3_set_f3(low, ground->lower);
    d3_set_f3(upp, ground->upper);

    res = htrdr_slab_trace_ray(ground->htrdr, org, dir, range, low, upp,
      trace_ground, 32, &slab_ctx);
    if(res != RES_OK) goto error;
  }

exit:
  return res;
error:
  goto exit;
}

