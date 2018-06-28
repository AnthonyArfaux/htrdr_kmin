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
#include "htrdr_sky.h"

#include <star/svx.h>
#include <high_tune/htcp.h>

#include <rsys/dynamic_array_double.h>
#include <rsys/dynamic_array_size_t.h>
#include <rsys/hash_table.h>
#include <rsys/ref_count.h>

struct split {
  size_t index; /* Index of the current htcp voxel */
  double height; /* Absolute height where the next voxel starts */
};

#define DARRAY_NAME split
#define DARRAY_DATA struct split
#include <rsys/dynamic_array.h>

struct build_octree_context {
  const struct htcp_desc* htcp_desc;
  const struct darray_split* svx2htcp_z;
  double dst_max; /* Max traversal distance */
  double tau_threshold; /* Threshold criteria for the merge process */
};

struct vertex {
  double x;
  double y;
  double z;
};

static char
vertex_eq(const struct vertex* v0, const struct vertex* v1)
{
  return eq_eps(v0->x, v1->x, 1.e-6)
      && eq_eps(v0->y, v1->y, 1.e-6)
      && eq_eps(v0->z, v1->z, 1.e-6);
}

/* Generate the htable_vertex data structure */
#define HTABLE_NAME vertex
#define HTABLE_KEY struct vertex
#define HTABLE_DATA size_t
#define HTABLE_KEY_FUNCTOR_EQ vertex_eq
#include <rsys/hash_table.h>

struct octree_data {
  struct htable_vertex vertex2id; /* Map a coordinate to its vertex id */
  struct darray_double vertices; /* Array of unique vertices */
  struct darray_double data;
  struct darray_size_t cells;
};

struct htrdr_sky {
  struct svx_tree* clouds;
  struct htcp* htcp;

  struct htcp_desc htcp_desc;

  /* Map the index in Z from the regular SVX to the irregular HTCP data */
  struct darray_split svx2htcp_z;

  ref_T ref;
  struct htrdr* htrdr;
};

/*******************************************************************************
 * Helper function
 ******************************************************************************/
static INLINE void
octree_data_init(struct octree_data* data)
{
  ASSERT(data);
  htable_vertex_init(NULL, &data->vertex2id);
  darray_double_init(NULL, &data->vertices);
  darray_double_init(NULL, &data->data);
  darray_size_t_init(NULL, &data->cells);
}

static INLINE void
octree_data_release(struct octree_data* data)
{
  ASSERT(data);
  htable_vertex_release(&data->vertex2id);
  darray_double_release(&data->vertices);
  darray_double_release(&data->data);
  darray_size_t_release(&data->cells);
}

static INLINE void
register_leaf
  (const struct svx_voxel* leaf,
   const size_t ileaf,
   void* context)
{
  struct octree_data* ctx = context;
  struct vertex v[8];
  const double* data;
  int i;
  ASSERT(leaf && ctx);
  (void)ileaf;

  data = leaf->data;

  /* Compute the leaf vertices */
  v[0].x = leaf->lower[0]; v[0].y = leaf->lower[1]; v[0].z = leaf->lower[2];
  v[1].x = leaf->upper[0]; v[1].y = leaf->lower[1]; v[1].z = leaf->lower[2];
  v[2].x = leaf->lower[0]; v[2].y = leaf->upper[1]; v[2].z = leaf->lower[2];
  v[3].x = leaf->upper[0]; v[3].y = leaf->upper[1]; v[3].z = leaf->lower[2];
  v[4].x = leaf->lower[0]; v[4].y = leaf->lower[1]; v[4].z = leaf->upper[2];
  v[5].x = leaf->upper[0]; v[5].y = leaf->lower[1]; v[5].z = leaf->upper[2];
  v[6].x = leaf->lower[0]; v[6].y = leaf->upper[1]; v[6].z = leaf->upper[2];
  v[7].x = leaf->upper[0]; v[7].y = leaf->upper[1]; v[7].z = leaf->upper[2];

  FOR_EACH(i, 0, 8) {
    size_t *pid = htable_vertex_find(&ctx->vertex2id, v+i);
    size_t id;
    if(pid) {
      id = *pid;
    } else { /* Register the leaf vertex */
      id = darray_double_size_get(&ctx->vertices)/3;
      CHK(RES_OK == htable_vertex_set(&ctx->vertex2id, v+i, &id));
      CHK(RES_OK == darray_double_push_back(&ctx->vertices, &v[i].x));
      CHK(RES_OK == darray_double_push_back(&ctx->vertices, &v[i].y));
      CHK(RES_OK == darray_double_push_back(&ctx->vertices, &v[i].z));
    }
    /* Add the vertex id to the leaf cell */
    CHK(RES_OK == darray_size_t_push_back(&ctx->cells, &id));
  }
  FOR_EACH(i, 0, HTRDR_SKY_SVX_PROPS_COUNT__) {
    CHK(RES_OK == darray_double_push_back(&ctx->data, data+i));
  }
}

static double
compute_kext(const struct htcp_desc* htcp_desc, const size_t xyz[3])
{
  const double c_ext = 6.e-10; /* Extinction cross section in m^2.particle^-1 */
  const double rho_air = 1.293; /* Volumic mass of terrestrial air in kg.m^-3 */
  const double rho_h2o = 1000; /* Volumic mass of water in kg.m^-3 */
  const double sigma = 0.1; /* Std deviation of the `r' parameter */
  const double r = 9.76617; /* Modal radius of water particles in um */
  double ql; /* Mass of the liquid water in the voxel in kg.m^-3 */
  double v; /* typical volume of a particle */
  double nv; /* Number density  */
  double k_ext; /* Extinction coefficient in m^-1 */
  ASSERT(xyz && htcp_desc);

  ql = htcp_desc_RCT_at(htcp_desc, xyz[0], xyz[1], xyz[2], 0/*time*/);
  v = 4*PI*(r*r*r)*exp(4.5*sigma*sigma) / 3.0;
  nv = 1.e18 * ql * rho_air / (v * rho_h2o);
  k_ext = c_ext*nv;

  return k_ext;
}

static void
vox_get(const size_t xyz[3], void* dst, void* context)
{
  struct build_octree_context* ctx = context;
  double kext_min;
  double kext_max;
  double* pdbl = dst;
  ASSERT(xyz && dst && context);

  if(!ctx->htcp_desc->irregular_z) {
    kext_min = kext_max = compute_kext(ctx->htcp_desc, xyz);
  } else {
    size_t ivox[3];
    size_t ivox_next;
    ASSERT(xyz[2] < darray_split_size_get(ctx->svx2htcp_z));

    ivox[0] = xyz[0];
    ivox[1] = xyz[1];
    ivox[2] = darray_split_cdata_get(ctx->svx2htcp_z)[xyz[2]].index;

    kext_min = kext_max = compute_kext(ctx->htcp_desc, ivox);

    /* Define if the SVX voxel is overlapped by 2 HTCP voxels */
    ivox_next = xyz[2] + 1 < darray_split_size_get(ctx->svx2htcp_z)
      ? darray_split_cdata_get(ctx->svx2htcp_z)[xyz[2] + 1].index
      : ivox[2];

    if(ivox_next != ivox[2]) {
      double kext;
      ASSERT(ivox[2] < ivox_next);
      ivox[2] = ivox_next;
      kext = compute_kext(ctx->htcp_desc, ivox);
      kext_min = MMIN(kext_min, kext);
      kext_max = MMAX(kext_max, kext);
    }
  }

  pdbl[HTRDR_SKY_SVX_Kext_MIN] = kext_min;
  pdbl[HTRDR_SKY_SVX_Kext_MAX] = kext_max;
}

static void
vox_merge(void* dst, const void* voxels[], const size_t nvoxs, void* ctx)
{
  double* pdbl = dst;
  double kext_min = DBL_MAX;
  double kext_max =-DBL_MAX;
  size_t ivox;
  ASSERT(dst && voxels && nvoxs);
  (void)ctx;

  FOR_EACH(ivox, 0, nvoxs) {
    const double* voxel_data = voxels[ivox];
    ASSERT(voxel_data[HTRDR_SKY_SVX_Kext_MIN]
        <= voxel_data[HTRDR_SKY_SVX_Kext_MAX]);
    kext_min = MMIN(kext_min, voxel_data[HTRDR_SKY_SVX_Kext_MIN]);
    kext_max = MMAX(kext_max, voxel_data[HTRDR_SKY_SVX_Kext_MAX]);
  }
  pdbl[HTRDR_SKY_SVX_Kext_MIN] = kext_min;
  pdbl[HTRDR_SKY_SVX_Kext_MAX] = kext_max;
}

static int
vox_challenge_merge(const void* voxels[], const size_t nvoxs, void* context)
{
  struct build_octree_context* ctx = context;
  double kext_min = DBL_MAX;
  double kext_max =-DBL_MAX;
  size_t ivox;
  ASSERT(voxels && nvoxs);

  FOR_EACH(ivox, 0, nvoxs) {
    const double* voxel_data = voxels[ivox];
    kext_min = MMIN(kext_min, voxel_data[HTRDR_SKY_SVX_Kext_MIN]);
    kext_max = MMAX(kext_max, voxel_data[HTRDR_SKY_SVX_Kext_MAX]);
  }
  return (kext_max - kext_min)*ctx->dst_max <= ctx->tau_threshold;
}

static void
release_sky(ref_T* ref)
{
  struct htrdr_sky* sky;
  ASSERT(ref);
  sky = CONTAINER_OF(ref, struct htrdr_sky, ref);
  if(sky->clouds) SVX(tree_ref_put(sky->clouds));
  if(sky->htcp) HTCP(ref_put(sky->htcp));
  darray_split_release(&sky->svx2htcp_z);
  MEM_RM(sky->htrdr->allocator, sky);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_sky_create
  (struct htrdr* htrdr,
   const char* htcp_filename,
   struct htrdr_sky** out_sky)
{
  struct htrdr_sky* sky = NULL;
  struct svx_voxel_desc vox_desc = SVX_VOXEL_DESC_NULL;
  struct build_octree_context ctx;
  size_t nvoxs[3];
  double sz[3];
  double low[3];
  double upp[3];
  res_T res = RES_OK;
  ASSERT(htrdr && htcp_filename && out_sky);

  sky = MEM_CALLOC(htrdr->allocator, 1, sizeof(*sky));
  if(!sky) {
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&sky->ref);
  sky->htrdr = htrdr;
  darray_split_init(htrdr->allocator, &sky->svx2htcp_z);

  res = htcp_create(&htrdr->logger, htrdr->allocator, htrdr->verbose, &sky->htcp);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "could not create the loader of cloud properties.\n");
    goto error;
  }

  res = htcp_load(sky->htcp, htcp_filename);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "error loading the cloud properties -- `%s'.\n",
      htcp_filename);
    goto error;
  }

  HTCP(get_desc(sky->htcp, &sky->htcp_desc));

  /* Define the number of voxels */
  nvoxs[0] = sky->htcp_desc.spatial_definition[0];
  nvoxs[1] = sky->htcp_desc.spatial_definition[1];
  nvoxs[2] = sky->htcp_desc.spatial_definition[2];

  /* Define the octree AABB */
  low[0] = sky->htcp_desc.lower[0];
  low[1] = sky->htcp_desc.lower[1];
  low[2] = sky->htcp_desc.lower[2];
  upp[0] = low[0] + (double)nvoxs[0] * sky->htcp_desc.vxsz_x;
  upp[1] = low[1] + (double)nvoxs[1] * sky->htcp_desc.vxsz_y;

  if(!sky->htcp_desc.irregular_z) { /* Refular voxel size in Z */
    upp[2] = low[2] + (double)nvoxs[2] * sky->htcp_desc.vxsz_z[0];
  } else { /* Irregular voxel size along Z */
    double min_vxsz_z;
    double len_z;
    size_t nsplits;
    size_t iz, iz2;;

    /* Find the min voxel size along Z and compute the length of a Z column */
    len_z = 0;
    min_vxsz_z = DBL_MAX;
    FOR_EACH(iz, 0, sky->htcp_desc.spatial_definition[2]) {
      len_z += sky->htcp_desc.vxsz_z[iz];
      min_vxsz_z = MMIN(min_vxsz_z, sky->htcp_desc.vxsz_z[iz]);
    }
    /* Allocate the svx2htcp LUT. The LUT is a regular table whose absolute
     * size is the size of a Z column in the htcp file. The size of its cells
     * is the minimal voxel size in Z of the htcp file */
    nsplits = (size_t)ceil(len_z / min_vxsz_z);
    res = darray_split_resize(&sky->svx2htcp_z, nsplits);
    if(res != RES_OK) {
      htrdr_log_err(htrdr,
        "could not allocate the table mapping regular to irregular Z.\n");
      goto error;
    }
    /* Setup the svx2htcp LUT. Each LUT entry stores the index of the current Z
     * voxel in the htcp file that overlaps the entry lower bound as well as the
     * lower bound in Z of the the next htcp voxel. */
    iz2 = 0;
    upp[2] = low[2] + sky->htcp_desc.vxsz_z[iz2];
    FOR_EACH(iz, 0, nsplits) {
      const double upp_z = (double)(iz + 1) * min_vxsz_z;
      darray_split_data_get(&sky->svx2htcp_z)[iz].index = iz2;
      darray_split_data_get(&sky->svx2htcp_z)[iz].height = upp[2];
      if(upp_z >= upp[2]) upp[2] += sky->htcp_desc.vxsz_z[++iz2];
    }
    ASSERT(eq_eps(upp[2] - low[2], len_z, 1.e-6));
  }

  /* Compute the size of of the AABB */
  sz[0] = upp[0] - low[0];
  sz[1] = upp[1] - low[1];
  sz[2] = upp[2] - low[2];

  /* Setup the build context */
  ctx.htcp_desc = &sky->htcp_desc;
  ctx.svx2htcp_z = &sky->svx2htcp_z;
  ctx.dst_max = sz[2];
  ctx.tau_threshold = 0.0;

  /* Setup the voxel descriptor */
  vox_desc.get = vox_get;
  vox_desc.merge = vox_merge;
  vox_desc.challenge_merge = vox_challenge_merge;
  vox_desc.context = &ctx;
  vox_desc.size = sizeof(double)*HTRDR_SKY_SVX_PROPS_COUNT__;

  /* Create the octree */
  res = svx_octree_create(htrdr->svx, low, upp, nvoxs, &vox_desc, &sky->clouds);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "could not create the octree of the cloud properties.\n");
    goto error;
  }

exit:
  *out_sky = sky;
  return res;
error:
  if(sky) {
    htrdr_sky_ref_put(sky);
    sky = NULL;
  }
  goto exit;
}

void
htrdr_sky_ref_get(struct htrdr_sky* sky)
{
  ASSERT(sky);
  ref_get(&sky->ref);
}

void
htrdr_sky_ref_put(struct htrdr_sky* sky)
{
  ASSERT(sky);
  ref_put(&sky->ref, release_sky);
}

struct svx_tree*
htrdr_sky_get_svx_tree(struct htrdr_sky* sky)
{
  ASSERT(sky);
  return sky->clouds;
}

res_T
htrdr_sky_dump_clouds_vtk(const struct htrdr_sky* sky, FILE* stream)
{
  struct svx_tree_desc desc;
  struct octree_data data;
  const double* leaf_data;
  size_t nvertices;
  size_t ncells;
  size_t i;
  ASSERT(sky && stream);

  octree_data_init(&data);
  SVX(tree_get_desc(sky->clouds, &desc));
  ASSERT(desc.type == SVX_OCTREE);

  /* Register leaf data */
  SVX(tree_for_each_leaf(sky->clouds, register_leaf, &data));
  nvertices = darray_double_size_get(&data.vertices) / 3/*#coords per vertex*/;
  ncells = darray_size_t_size_get(&data.cells)/8/*#ids per cell*/;
  ASSERT(ncells == desc.nleaves);

  /* Write headers */
  fprintf(stream, "# vtk DataFile Version 2.0\n");
  fprintf(stream, "Volume\n");
  fprintf(stream, "ASCII\n");
  fprintf(stream, "DATASET UNSTRUCTURED_GRID\n");

  /* Write vertex coordinates */
  fprintf(stream, "POINTS %lu float\n", (unsigned long)nvertices);
  FOR_EACH(i, 0, nvertices) {
    fprintf(stream, "%g %g %g\n",
      SPLIT3(darray_double_cdata_get(&data.vertices) + i*3));
  }

  /* Write the cells */
  fprintf(stream, "CELLS %lu %lu\n",
    (unsigned long)ncells,
    (unsigned long)(ncells*(8/*#verts per cell*/ + 1/*1st field of a cell*/)));
  FOR_EACH(i, 0, ncells) {
    fprintf(stream, "8 %lu %lu %lu %lu %lu %lu %lu %lu\n",
      (unsigned long)darray_size_t_cdata_get(&data.cells)[i*8+0],
      (unsigned long)darray_size_t_cdata_get(&data.cells)[i*8+1],
      (unsigned long)darray_size_t_cdata_get(&data.cells)[i*8+2],
      (unsigned long)darray_size_t_cdata_get(&data.cells)[i*8+3],
      (unsigned long)darray_size_t_cdata_get(&data.cells)[i*8+4],
      (unsigned long)darray_size_t_cdata_get(&data.cells)[i*8+5],
      (unsigned long)darray_size_t_cdata_get(&data.cells)[i*8+6],
      (unsigned long)darray_size_t_cdata_get(&data.cells)[i*8+7]);
  }

  /* Write the cell type */
  fprintf(stream, "CELL_TYPES %lu\n", (unsigned long)ncells);
  FOR_EACH(i, 0, ncells) fprintf(stream, "11\n");

  /* Write the cell data */
  leaf_data = darray_double_cdata_get(&data.data);
  fprintf(stream, "CELL_DATA %lu\n", (unsigned long)ncells);
  fprintf(stream, "SCALARS Val double %d\n", HTRDR_SKY_SVX_PROPS_COUNT__);
  fprintf(stream, "LOOKUP_TABLE default\n");
  FOR_EACH(i, 0, ncells) {
    size_t idata;
    FOR_EACH(idata, 0, HTRDR_SKY_SVX_PROPS_COUNT__) {
      fprintf(stream, "%g ", leaf_data[i*HTRDR_SKY_SVX_PROPS_COUNT__ + idata]);
    }
    fprintf(stream, "\n");
  }
  octree_data_release(&data);
  return RES_OK;
}

