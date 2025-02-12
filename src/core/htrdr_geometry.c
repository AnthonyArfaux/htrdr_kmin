/* Copyright (C) 2018-2019, 2022-2025 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2025 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2025 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2025 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2025 Observatoire de Paris
 * Copyright (C) 2022-2025 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2025 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2025 Université Paul Sabatier
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

#define _POSIX_C_SOURCE 200112L /* strtok_r support */

#include "core/htrdr.h"
#include "core/htrdr_geometry.h"
#include "core/htrdr_interface.h"
#include "core/htrdr_log.h"
#include "core/htrdr_materials.h"
#include "core/htrdr_slab.h"

#include <aw.h>
#include <star/s3d.h>

#include <rsys/clock_time.h>
#include <rsys/cstr.h>
#include <rsys/dynamic_array_double.h>
#include <rsys/dynamic_array_size_t.h>
#include <rsys/double2.h>
#include <rsys/double3.h>
#include <rsys/float2.h>
#include <rsys/float3.h>
#include <rsys/hash_table.h>
#include <rsys/ref_count.h>

#include <string.h> /* strtok_r */

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
  struct s3d_hit hit_from;

  s3d_hit_filter_function_T user_filter; /* May be NULL */
  void* user_filter_data; /* May be NULL */
};
#define RAY_CONTEXT_NULL__ {S3D_HIT_NULL__, NULL, NULL}
static const struct ray_context RAY_CONTEXT_NULL = RAY_CONTEXT_NULL__;

struct htrdr_geometry {
  struct s3d_scene_view* view;
  float lower[3]; /* Ground lower bound */
  float upper[3]; /* Ground upper bound */
  int repeat; /* Make the geom infinite in X and Y */

  /* A empirical value relative to the extent of the geometry that represents
   * the threshold below which a numerical problem could occur. */
  float epsilon;

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

/* Return 1 if the submitted hit is actually a self hit, i.e. if the ray starts
 * the primitive from which it starts */
static INLINE int
self_hit
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   const float ray_range[2],
   const struct s3d_hit* hit_from)
{
  ASSERT(hit && hit_from);
  (void)ray_org, (void)ray_dir;

  /* Internally, Star-3D relies on Embree that, due to numerically inaccuracy,
   * can find hits whose distances are not strictly included in the submitted
   * ray range. Discard these hits. */
  if(hit->distance <= ray_range[0] || hit->distance >= ray_range[1])
    return 1;

  /* The hit primitive is the one from which the the ray starts. Discard these
   * hits */
  if(S3D_PRIMITIVE_EQ(&hit->prim, &hit_from->prim))
    return 1;

  /* If the targeted point is near of the origin, check that it lies on an
   * edge/vertex shared by the 2 primitives. In such situation, we assume that
   * it is a self intersection and discard this hit */
  if(!S3D_HIT_NONE(hit_from)
  && eq_epsf(hit->distance, 0, 1.e-1f)
  && hit_on_edge(hit_from) 
  && hit_on_edge(hit)) {
    return 1;
  }

  /* No self hit */
  return 0;
}

static int
geometry_filter
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   const float ray_range[2],
   void* ray_data,
   void* filter_data)
{
  const struct ray_context* ray_ctx = ray_data;
  (void)filter_data;

  if(!ray_ctx) /* Nothing to do */
    return 0;

  if(self_hit(hit, ray_org, ray_dir, ray_range, &ray_ctx->hit_from))
    return 1; /* Discard this hit */

  if(!ray_ctx->user_filter) /* That's all */
    return 0;

  return ray_ctx->user_filter /* Invoke user filtering */
    (hit, ray_org, ray_dir, ray_range, ray_ctx->user_filter_data, filter_data);
}

static res_T
parse_shape_interface
  (struct htrdr* htrdr,
   struct htrdr_materials* mats,
   const char* name,
   struct htrdr_interface* interf)
{
  struct str str;
  char* mtl_name0 = NULL;
  char* mtl_name1 = NULL;
  char* mtl_name2 = NULL;
  char* mtl_name_front = NULL;
  char* mtl_name_thin = NULL;
  char* mtl_name_back = NULL;
  char* tk_ctx = NULL;
  int has_front = 0;
  int has_thin = 0;
  int has_back = 0;
  res_T res = RES_OK;
  ASSERT(htrdr && mats && name && interf);

  str_init(htrdr_get_allocator(htrdr), &str);

  /* Locally copy the string to parse */
  res = str_set(&str, name);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "Could not locally copy the shape material string `%s' -- %s.\n",
      name, res_to_cstr(res));
    goto error;
  }

  /* Reset the interface */
  memset(interf, 0, sizeof(*interf));

  /* Parse the name of the front/back/thin materials */
  mtl_name0 = strtok_r(str_get(&str), ":", &tk_ctx);
  mtl_name1 = strtok_r(NULL, ":", &tk_ctx);
  mtl_name2 = strtok_r(NULL, ":", &tk_ctx);
  ASSERT(mtl_name0); /* This can't be NULL */
  mtl_name_front = mtl_name0;
  if(mtl_name2) {
    mtl_name_thin = mtl_name1;
    mtl_name_back = mtl_name2;
  } else {
    mtl_name_thin = NULL;
    mtl_name_back = mtl_name1;
  }

  if(!mtl_name_back) {
    htrdr_log_err(htrdr,
      "The back material name is missing `%s'.\n", name);
    res = RES_BAD_ARG;
    goto error;
  }

  /* Fetch the interface material if any */
  if(mtl_name_thin) {
    has_thin = htrdr_materials_find_mtl(mats, mtl_name_thin, &interf->mtl_thin);
    if(!has_thin) {
      htrdr_log_err(htrdr,
        "Invalid interface `%s'. The interface material `%s' is unknown.\n",
        name, mtl_name_thin);
      res = RES_BAD_ARG;
      goto error;
    }
  }

  /* Fetch the front material */
  has_front = htrdr_materials_find_mtl(mats, mtl_name_front, &interf->mtl_front);
  if(!has_front) {
    htrdr_log_err(htrdr,
      "Invalid interface `%s'. The front material `%s' is unknown.\n",
      name, mtl_name_front);
    res = RES_BAD_ARG;
    goto error;
  }

  /* Fetch the back material */
  has_back = htrdr_materials_find_mtl(mats, mtl_name_back, &interf->mtl_back);
  if(!has_back) {
    htrdr_log_err(htrdr,
      "Invalid interface `%s'. The back material `%s' is unknown.\n",
      name, mtl_name_back);
    res = RES_BAD_ARG;
    goto error;
  }

exit:
  str_release(&str);
  return res;
error:
  *interf = HTRDR_INTERFACE_NULL;
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

  FOR_EACH(iface, mtl->face_id, mtl->face_id+mtl->faces_count) {
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

  res = s3d_shape_create_mesh(htrdr_get_s3d(htrdr), &shape);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "Error creating a Star-3D shape -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  vdata.usage = S3D_POSITION;
  vdata.type = S3D_FLOAT3;
  vdata.get = get_position;

  res = s3d_mesh_setup_indexed_vertices
    (shape,
     (unsigned int)(darray_size_t_size_get(indices)/3),
     get_indices,
     (unsigned int)(darray_double_size_get(positions)/3),
     &vdata, 1,
     &mesh);
  if(res != RES_OK){
    htrdr_log_err(htrdr, "Could not setup the Star-3D shape -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  res = s3d_mesh_set_hit_filter_function(shape, geometry_filter, NULL);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "Could not setup the Star-3D hit filter function of the geometry "
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
setup_geometry
  (struct htrdr_geometry* geom,
   struct htrdr_materials* mats,
   const char* obj_filename)
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
  ASSERT(geom && mats && obj_filename);

  htable_vertex_init(htrdr_get_allocator(geom->htrdr), &vertices);
  darray_double_init(htrdr_get_allocator(geom->htrdr), &positions);
  darray_size_t_init(htrdr_get_allocator(geom->htrdr), &indices);

  res = aw_obj_create
    (htrdr_get_logger(geom->htrdr),
     htrdr_get_allocator(geom->htrdr),
     htrdr_get_verbosity_level(geom->htrdr),
     &obj);
  if(res != RES_OK) {
    htrdr_log_err(geom->htrdr, "Could not create the obj loader -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  res = s3d_scene_create(htrdr_get_s3d(geom->htrdr), &scene);
  if(res != RES_OK) {
    htrdr_log_err(geom->htrdr, "Could not create the Star-3D scene -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  /* Load the geometry data */
  res = aw_obj_load(obj, obj_filename);
  if(res != RES_OK) goto error;

  /* Fetch the descriptor of the loaded geometry */
  AW(obj_get_desc(obj, &desc));

  if(desc.usemtls_count == 0) {
    htrdr_log_err(geom->htrdr, "The obj `%s' has no material.\n", obj_filename);
    res = RES_BAD_ARG;
    goto error;
  }

  /* Setup the geometry */
  FOR_EACH(iusemtl, 0, desc.usemtls_count) {
    struct aw_obj_named_group mtl;
    struct htrdr_interface interf;
    unsigned ishape;

    AW(obj_get_mtl(obj, iusemtl , &mtl));

    res = parse_shape_interface(geom->htrdr, mats, mtl.name, &interf);
    if(res != RES_OK) goto error;

    res = setup_mesh
      (geom->htrdr, obj_filename, obj, &mtl, &positions, &indices, &vertices);
    if(res != RES_OK) goto error;

    res = create_s3d_shape(geom->htrdr, &positions, &indices, &shape);
    if(res != RES_OK) goto error;

    S3D(shape_get_id(shape, &ishape));
    res = htable_interface_set(&geom->interfaces, &ishape, &interf);
    if(res != RES_OK) {
      htrdr_log_err(geom->htrdr,
        "Could not map the Star-3D shape to its Star-Materials -- %s.\n",
        res_to_cstr(res));
      goto error;
    }

    res = s3d_scene_attach_shape(scene, shape);
    if(res != RES_OK) {
      htrdr_log_err(geom->htrdr,
        "Could not attach a Star-3D shape to the Star-3D scene -- %s.\n",
        res_to_cstr(res));
      goto error;
    }

    S3D(shape_ref_put(shape));
    shape = NULL;
  }

  res = s3d_scene_view_create(scene, S3D_GET_PRIMITIVE|S3D_TRACE, &geom->view);
  if(res != RES_OK) {
    htrdr_log_err(geom->htrdr,
      "Could not create the Star-3D scene view -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  res = s3d_scene_view_get_aabb(geom->view, geom->lower, geom->upper);
  if(res != RES_OK) {
    htrdr_log_err(geom->htrdr,
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
  struct htrdr_geometry* geom;
  struct htrdr* htrdr;
  ASSERT(ref);
  geom = CONTAINER_OF(ref, struct htrdr_geometry, ref);
  if(geom->view) S3D(scene_view_ref_put(geom->view));
  htable_interface_release(&geom->interfaces);
  htrdr = geom->htrdr;
  MEM_RM(htrdr_get_allocator(geom->htrdr), geom);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_geometry_create
  (struct htrdr* htrdr,
   const char* obj_filename, /* May be NULL */
   struct htrdr_materials* mats,
   struct htrdr_geometry** out_ground)
{
  char buf[128];
  struct htrdr_geometry* geom = NULL;
  double low[3];
  double upp[3];
  double tmp[3];
  double extent;
  struct time t0, t1;
  res_T res = RES_OK;
  ASSERT(htrdr && obj_filename && mats && out_ground);

  geom = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*geom));
  if(!geom) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "%s: could not allocate the geom data structure -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }
  ref_init(&geom->ref);
  f3_splat(geom->lower, (float)INF);
  f3_splat(geom->upper,-(float)INF);
  htable_interface_init(htrdr_get_allocator(htrdr), &geom->interfaces);
  htrdr_ref_get(htrdr);
  geom->htrdr = htrdr;

  htrdr_log(geom->htrdr, "Loading geometry from `%s'.\n", obj_filename);
  time_current(&t0);
  res = setup_geometry(geom, mats, obj_filename);
  if(res != RES_OK) goto error;
  time_sub(&t0, time_current(&t1), &t0);
  time_dump(&t0, TIME_ALL, NULL, buf, sizeof(buf));
  htrdr_log(geom->htrdr, "Setup geom in %s\n", buf);

  htrdr_geometry_get_aabb(geom, low, upp);
  extent = d3_len(d3_sub(tmp, upp, low));
  geom->epsilon = MMAX((float)(extent * 1e-6), FLT_EPSILON);

exit:
  *out_ground = geom;
  return res;
error:
  if(geom) {
    htrdr_geometry_ref_put(geom);
    geom = NULL;
  }
  goto exit;
}

void
htrdr_geometry_ref_get(struct htrdr_geometry* geom)
{
  ASSERT(geom);
  ref_get(&geom->ref);
}

void
htrdr_geometry_ref_put(struct htrdr_geometry* geom)
{
  ASSERT(geom);
  ref_put(&geom->ref, release_ground);
}

void
htrdr_geometry_get_interface
  (struct htrdr_geometry* geom,
   const struct s3d_hit* hit,
   struct htrdr_interface* out_interface)
{
  struct htrdr_interface* interf = NULL;
  ASSERT(geom && hit && out_interface);

  interf = htable_interface_find(&geom->interfaces, &hit->prim.geom_id);
  ASSERT(interf);

  *out_interface = *interf;
}

void
htrdr_geometry_get_hit_position
  (const struct htrdr_geometry* geom,
   const struct s3d_hit* hit,
   double position[3])
{
  struct s3d_attrib attr;
  ASSERT(geom && hit && position && !S3D_HIT_NONE(hit));
  (void)geom;

  S3D(primitive_get_attrib(&hit->prim, S3D_POSITION, hit->uv, &attr));
  position[0] = attr.value[0];
  position[1] = attr.value[1];
  position[2] = attr.value[2];
}

res_T
htrdr_geometry_trace_ray
  (struct htrdr_geometry* geom,
   const struct htrdr_geometry_trace_ray_args* args,
   struct s3d_hit* hit)
{
  struct ray_context ray_ctx = RAY_CONTEXT_NULL;
  float ray_org[3];
  float ray_dir[3];
  float ray_range[2];
  res_T res = RES_OK;
  ASSERT(geom && args && hit);

  f3_set_d3(ray_org, args->ray_org);
  f3_set_d3(ray_dir, args->ray_dir);
  f2_set_d2(ray_range, args->ray_range);
  ray_ctx.hit_from = args->hit_from;
  ray_ctx.user_filter = args->filter;
  ray_ctx.user_filter_data = args->filter_context;

  res = s3d_scene_view_trace_ray
    (geom->view, ray_org, ray_dir, ray_range, &ray_ctx, hit);
  if(res != RES_OK) {
    htrdr_log_err(geom->htrdr,
      "%s: could not trace the ray against the geometry -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }

exit:
  return res;
error:
  goto exit;
}

res_T
htrdr_geometry_find_closest_point
  (struct htrdr_geometry* geom,
   const double pos[3],
   const double radius,
   struct s3d_hit* hit)
{
  float query_pos[3];
  float query_radius;
  res_T res = RES_OK;
  ASSERT(geom && pos && hit);

  query_radius = (float)radius;
  f3_set_d3(query_pos, pos);

  /* Closest point query */
  res = s3d_scene_view_closest_point
    (geom->view, query_pos, query_radius, NULL, hit);
  if(res != RES_OK) {
    htrdr_log_err(geom->htrdr,
      "%s: could not query the closest point to the geometry -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }

exit:
  return res;
error:
  goto exit;
}

void
htrdr_geometry_get_aabb
  (const struct htrdr_geometry* geom,
   double lower[3],
   double upper[3])
{
  ASSERT(geom && lower && upper);
  d3_set_f3(lower, geom->lower);
  d3_set_f3(upper, geom->upper);
}

double
htrdr_geometry_get_epsilon(const struct htrdr_geometry* geom)
{
  ASSERT(geom);
  return geom->epsilon;
}
