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

#include <high_tune/htcp.h>
#include <star/svx.h>

#include <rsys/dynamic_array_double.h>
#include <rsys/dynamic_array_size_t.h>
#include <rsys/hash_table.h>
#include <rsys/math.h>
#include <rsys/mem_allocator.h>

#include <getopt.h>

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

struct htrdr_args {
  const char* input;
  const char* output;

  int dump_vtk;
  int verbose;
  int quit;  /* Qui the application */
};
#define HTRDR_ARGS_DEFAULT__ {NULL, NULL, 0, 0, 0}
static const struct htrdr_args HTRDR_ARGS_DEFAULT = HTRDR_ARGS_DEFAULT__;

struct octree_data {
  struct htable_vertex vertex2id; /* Map a coordinate to an id */
  struct darray_double vertices;
  struct darray_double data;
  struct darray_size_t cells;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
print_help(const char* cmd)
{
  ASSERT(cmd);
  printf("Usage: %s -i INPUT [OPIONS]\n\n", cmd);
  printf(
"  -d               dump octree data to OUTPUT wrt the VTK ASCII file format.\n");
  printf(
"  -h               display this help and exit.\n");
  printf(
"  -i INPUT         path of the input HTCP file.\n");
  printf(
"  -o OUTPUT        file where data are written. If not defined, data are\n"
"                   written to standard output.\n");
  printf(
"  -v               make the programe more verbose.\n");
  printf("\n");
  printf(
"%s (C) 2018 Université Paul Sabatier, |Meso|Star>. This is free software\n"
"released under the GNU GPL license, version 3 or later. You are free to change\n"
"or redistribute it under certain conditions <http://gnu.org/licenses/gpl.html>.\n",
    cmd);
}

static void
htrdr_args_release(struct htrdr_args* args)
{
  ASSERT(args);
  *args = HTRDR_ARGS_DEFAULT;
}

static res_T
htrdr_args_init(struct htrdr_args* args, int argc, char** argv)
{
  int opt;
  res_T res = RES_OK;
  ASSERT(args && argc && argv);

  *args = HTRDR_ARGS_DEFAULT;

  while((opt = getopt(argc, argv, "dhi:o:v")) != -1) {
    switch(opt) {
      case 'd': args->dump_vtk = 1; break;
      case 'h':
        print_help(argv[0]);
        htrdr_args_release(args);
        args->quit = 1;
        goto exit;
      case 'i': args->input = optarg; break;
      case 'o': args->output = optarg; break;
      case 'v': args->verbose = 1; break;
      default: res = RES_BAD_ARG; break;
    }
  }
  if(res != RES_OK) {
    if(optarg) {
      fprintf(stderr, "%s: invalid option argumet '%s' -- '%c'\n",
        argv[0], optarg, opt);
    }
    goto error;
  }
  if(!args->input) {
    fprintf(stderr, "Missing input file.\n");
    res = RES_BAD_ARG;
    goto error;
  }
exit:
  return res;
error:
  htrdr_args_release(args);
  goto exit;
}

static void
octree_data_init(struct octree_data* data)
{
  ASSERT(data);
  htable_vertex_init(NULL, &data->vertex2id);
  darray_double_init(NULL, &data->vertices);
  darray_double_init(NULL, &data->data);
  darray_size_t_init(NULL, &data->cells);
}

static void
octree_data_release(struct octree_data* data)
{
  ASSERT(data);
  htable_vertex_release(&data->vertex2id);
  darray_double_release(&data->vertices);
  darray_double_release(&data->data);
  darray_size_t_release(&data->cells);
}

static void
vox_get(const size_t xyz[3], void* dst, void* ctx)
{
  struct htcp_desc* desc = ctx;
  double* pdbl = dst;
  ASSERT(xyz && dst && ctx);
  *pdbl = htcp_desc_RVT_at(desc, xyz[0], xyz[1], xyz[2], 0/*time*/);
}

static void
vox_merge(void* dst, const void* voxels[], const size_t nvoxs, void* ctx)
{
  double* pdbl = dst;
  double dbl_max = -DBL_MAX;
  size_t i;
  ASSERT(dst && voxels && nvoxs);
  (void)ctx;
  FOR_EACH(i, 0, nvoxs) dbl_max = MMAX(dbl_max, *(double*)(voxels[i]));
  *pdbl = dbl_max;
}

static int
vox_challenge_merge(const void* voxels[], const size_t nvoxs, void* ctx)
{
  size_t i;
  (void)nvoxs, (void)ctx;
  FOR_EACH(i, 0, nvoxs) {
    if(*(double*)(voxels[i]) != 0) break;
  }
  return i >= nvoxs;
}

static res_T
create_octree
  (struct svx_device* svx, struct htcp* htcp, struct svx_tree** out_octree)
{
  struct svx_tree* octree = NULL;
  struct svx_voxel_desc vox_desc = SVX_VOXEL_DESC_NULL;
  struct htcp_desc desc = HTCP_DESC_NULL;
  size_t nvoxs[3];
  double low[3], upp[3];
  res_T res = RES_OK;

  HTCP(get_desc(htcp, &desc));

  /* Currently only irregular Z dimension is not supported */
  if(desc.irregular_z) {
    fprintf(stderr, "Unexpected irregular Z size.\n");
    res = RES_BAD_ARG;
    goto error;
  }

  /* Define the number of voxels */
  nvoxs[0] = desc.spatial_definition[0];
  nvoxs[1] = desc.spatial_definition[1];
  nvoxs[2] = desc.spatial_definition[2];

  /* Define the octree AABB */
  low[0] = desc.lower[0];
  low[1] = desc.lower[1];
  low[2] = desc.lower[2];
  upp[0] = low[0] + (double)nvoxs[0] * desc.vxsz_x;
  upp[1] = low[1] + (double)nvoxs[1] * desc.vxsz_y;
  upp[2] = low[2] + (double)nvoxs[2] * desc.vxsz_z[0];

  /* Setup the voxel descriptor */
  vox_desc.get = vox_get;
  vox_desc.merge = vox_merge;
  vox_desc.challenge_merge = vox_challenge_merge;
  vox_desc.context = &desc;
  vox_desc.size = sizeof(double);

  /* Create the octree */
  res = svx_octree_create(svx, low, upp, nvoxs, &vox_desc, &octree);
  if(res != RES_OK) {
    fprintf(stderr, "Could not create the octree from the HTCP data.\n");
    goto error;
  }

exit:
  *out_octree = octree;
  return res;
error:
  if(octree) {
    SVX(tree_ref_put(octree));
    octree = NULL;
  }
  goto exit;
}

static void
register_leaf_vertices
  (const struct svx_voxel* leaf,
   const size_t ileaf,
   void* context)
{
  struct octree_data* ctx = context;
  struct vertex v[8];
  int i;
  ASSERT(leaf && ctx);
  (void)ileaf;

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
  CHK(RES_OK == darray_double_push_back(&ctx->data, leaf->data));
}

static res_T
dump_vtk(struct svx_tree* octree, FILE* stream)
{
  struct svx_tree_desc desc;
  struct octree_data data;
  size_t nvertices;
  size_t ncells;
  size_t i;
  res_T res = RES_OK;

  octree_data_init(&data);
  SVX(tree_get_desc(octree, &desc));
  ASSERT(desc.type == SVX_OCTREE);

  /* Register unique vertices */
  SVX(tree_for_each_leaf(octree, register_leaf_vertices, &data));
  nvertices = darray_double_size_get(&data.vertices) / 3/*#coords per vertex*/;
  ncells = darray_size_t_size_get(&data.cells)/8/*#ides per cell*/;
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
  fprintf(stream, "CELL_DATA %lu\n", (unsigned long)ncells);
  fprintf(stream, "SCALARS Val double 1\n");
  fprintf(stream, "LOOKUP_TABLE default\n");
  FOR_EACH(i, 0, ncells) {
    fprintf(stream, "%g\n", darray_double_cdata_get(&data.data)[i]);
  }

  octree_data_release(&data);
  return res;
}

/*******************************************************************************
 * Program
 ******************************************************************************/
int
main(int argc, char** argv)
{
  struct htrdr_args args;
  struct htcp* htcp = NULL;
  struct svx_device* svx = NULL;
  struct svx_tree* octree = NULL;
  res_T res = RES_OK;
  int err = 0;

  res = htrdr_args_init(&args, argc, argv);
  if(res != RES_OK) goto error;
  if(args.quit) goto exit;

  res = htcp_create(NULL, NULL, args.verbose, &htcp);
  if(res != RES_OK) goto error;
  res = htcp_load(htcp, args.input);
  if(res != RES_OK) goto error;

  res = svx_device_create(NULL, NULL, args.verbose, &svx);
  if(res != RES_OK) goto error;
  res = create_octree(svx, htcp, &octree);
  if(res != RES_OK) goto error;

  if(args.dump_vtk) {
    res = dump_vtk(octree, stdout);
    if(res != RES_OK) goto error;
    goto exit;
  }

exit:
  if(htcp) HTCP(ref_put(htcp));
  if(svx) SVX(device_ref_put(svx));
  if(octree) SVX(tree_ref_put(octree));
  if(mem_allocated_size()) {
    fprintf(stderr, "Memory leaks: %luBytes\n",
      (unsigned long)mem_allocated_size());
    err = -1;
  }
  return err;
error:
  err = -1;
  goto exit;
}
