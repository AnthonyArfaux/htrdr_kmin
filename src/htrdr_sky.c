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

#define _POSIX_C_SOURCE 200112L /* nextafterf support */

#include "htrdr.h"
#include "htrdr_c.h"
#include "htrdr_sky.h"
#include "htrdr_sun.h"

#include <star/svx.h>
#include <high_tune/htcp.h>
#include <high_tune/htmie.h>

#include <rsys/dynamic_array_double.h>
#include <rsys/dynamic_array_size_t.h>
#include <rsys/hash_table.h>
#include <rsys/ref_count.h>

#include <math.h>

#define DRY_AIR_MOLAR_MASS 0.0289644 /* In kg.mol^-1 */
#define GAS_CONSTANT 8.3144598 /* In kg.m^2.s^-2.mol^-1.K */

struct split {
  size_t index; /* Index of the current htcp voxel */
  double height; /* Absolute height where the next voxel starts */
};

#define DARRAY_NAME split
#define DARRAY_DATA struct split
#include <rsys/dynamic_array.h>

struct build_octree_context {
  const struct htrdr_sky* sky;
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
  const struct htrdr_sky* sky;
};

struct htrdr_sky {
  struct svx_tree* clouds;
  struct svx_tree_desc cloud_desc;

  struct htrdr_sun* sun;

  struct htcp* htcp;
  struct htmie* htmie;

  struct htcp_desc htcp_desc;

  /* Map the index in Z from the regular SVX to the irregular HTCP data */
  struct darray_split svx2htcp_z;

  /* Average cross section in Short Wave, i.e. in [380, 780] nanometers */
  double Ca_avg_sw;
  double Cs_avg_sw;

  /* Average asymmetry parameter in Short Wave, in [380, 780] nanometers */
  double g_avg_sw;

  ref_T ref;
  struct htrdr* htrdr;
};

/*******************************************************************************
 * Helper function
 ******************************************************************************/
/* Compute the dry air density in the cloud */
static FINLINE double
cloud_dry_air_density
  (const struct htcp_desc* desc,
   const size_t ivox[3]) /* Index of the voxel */
{
  double P = 0; /* Pressure in Pa */
  double T = 0; /* Temperature in K */
  ASSERT(desc);
  P = htcp_desc_PABST_at(desc, ivox[0], ivox[1], ivox[2], 0);
  T = htcp_desc_T_at(desc, ivox[0], ivox[1], ivox[2], 0);
  return (P*DRY_AIR_MOLAR_MASS)/(T*GAS_CONSTANT);
}

static INLINE void
octree_data_init(const struct htrdr_sky* sky, struct octree_data* data)
{
  ASSERT(data);
  htable_vertex_init(sky->htrdr->allocator, &data->vertex2id);
  darray_double_init(sky->htrdr->allocator, &data->vertices);
  darray_double_init(sky->htrdr->allocator, &data->data);
  darray_size_t_init(sky->htrdr->allocator, &data->cells);
  data->sky = sky;
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
  double kext_min;
  double kext_max;
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

  /* Register the leaf data */
  kext_max = htrdr_sky_fetch_svx_voxel_property
    (ctx->sky, HTRDR_Kext, HTRDR_SVX_MAX, HTRDR_ALL_COMPONENTS, -1/*FIXME*/, leaf);
  kext_min = htrdr_sky_fetch_svx_voxel_property
    (ctx->sky, HTRDR_Kext, HTRDR_SVX_MIN, HTRDR_ALL_COMPONENTS, -1/*FIXME*/, leaf);
  CHK(RES_OK == darray_double_push_back(&ctx->data, &kext_min));
  CHK(RES_OK == darray_double_push_back(&ctx->data, &kext_max));
}

static void
vox_get(const size_t xyz[3], void* dst, void* context)
{
  struct build_octree_context* ctx = context;
  double rct;
  double ka, ks, kext;
  double ka_min, ka_max;
  double ks_min, ks_max;
  double kext_min, kext_max;
  double rho_da; /* Dry air density */
  float* pflt = dst;
  ASSERT(xyz && dst && context);

  if(!ctx->sky->htcp_desc.irregular_z) {
    rho_da = cloud_dry_air_density(&ctx->sky->htcp_desc, xyz);
    rct = htcp_desc_RCT_at(&ctx->sky->htcp_desc, xyz[0], xyz[1], xyz[2], 0);
    ka_min = ka_max = ka = ctx->sky->Ca_avg_sw * rho_da * rct;
    ks_min = ks_max = ks = ctx->sky->Cs_avg_sw * rho_da * rct;
    kext_min = kext_max = kext = ka + ks;
  } else {
    size_t ivox[3];
    size_t ivox_next;
    ASSERT(xyz[2] < darray_split_size_get(&ctx->sky->svx2htcp_z));

    ivox[0] = xyz[0];
    ivox[1] = xyz[1];
    ivox[2] = darray_split_cdata_get(&ctx->sky->svx2htcp_z)[xyz[2]].index;

    rho_da = cloud_dry_air_density(&ctx->sky->htcp_desc, ivox);
    rct = htcp_desc_RCT_at(&ctx->sky->htcp_desc, ivox[0], ivox[1], ivox[2], 0);

    ka_min = ka_max = ka = ctx->sky->Ca_avg_sw * rho_da * rct;
    ks_min = ks_max = ks = ctx->sky->Cs_avg_sw * rho_da * rct;
    kext_min = kext_max = kext = ka + ks;

    /* Define if the SVX voxel is overlapped by 2 HTCP voxels */
    ivox_next = xyz[2] + 1 < darray_split_size_get(&ctx->sky->svx2htcp_z)
      ? darray_split_cdata_get(&ctx->sky->svx2htcp_z)[xyz[2] + 1].index
      : ivox[2];

    if(ivox_next != ivox[2]) {
      ASSERT(ivox[2] < ivox_next);
      ivox[2] = ivox_next;

      rho_da = cloud_dry_air_density(&ctx->sky->htcp_desc, ivox);
      rct = htcp_desc_RCT_at(&ctx->sky->htcp_desc, ivox[0], ivox[1], ivox[2], 0);
      ka = ctx->sky->Ca_avg_sw * rho_da * rct;
      ks = ctx->sky->Cs_avg_sw * rho_da * rct;
      kext = ka + ks;

      ka_min = MMIN(ka_min, ka);
      ka_max = MMAX(ka_max, ka);
      ks_min = MMIN(ks_min, ks);
      ks_max = MMAX(ks_max, ks);
      kext_min = MMIN(kext_min, kext);
      kext_max = MMAX(kext_max, kext);
    }
  }

  /* Ensure that the single precision bounds include their double precision
   * version. */
  if(ka_min != (float)ka_min) ka_min = nextafterf((float)ka_min,-FLT_MAX);
  if(ka_max != (float)ka_max) ka_max = nextafterf((float)ka_max, FLT_MAX);
  if(ks_min != (float)ks_min) ks_min = nextafterf((float)ks_min,-FLT_MAX);
  if(ks_max != (float)ks_max) ks_max = nextafterf((float)ks_max, FLT_MAX);
  if(kext_min != (float)kext_min) kext_min = nextafterf((float)kext_min,-FLT_MAX);
  if(kext_max != (float)kext_max) kext_max = nextafterf((float)kext_max, FLT_MAX);

  pflt[HTRDR_Ka   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = (float)ka_min;
  pflt[HTRDR_Ka   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = (float)ka_max;
  pflt[HTRDR_Ks   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = (float)ks_min;
  pflt[HTRDR_Ks   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = (float)ks_max;
  pflt[HTRDR_Kext * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = (float)kext_min;
  pflt[HTRDR_Kext * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = (float)kext_max;
}

static void
vox_merge(void* dst, const void* voxels[], const size_t nvoxs, void* context)
{
  float* pflt = dst;
  float ka_min = FLT_MAX;
  float ka_max =-FLT_MAX;
  float ks_min = FLT_MAX;
  float ks_max =-FLT_MAX;
  float kext_min = FLT_MAX;
  float kext_max =-FLT_MAX;
  size_t ivox;
  ASSERT(dst && voxels && nvoxs);
  (void)context;

  FOR_EACH(ivox, 0, nvoxs) {
    const float* vox_data = (const float*)voxels[ivox];
    const float* ka = vox_data + HTRDR_Ka * HTRDR_SVX_OPS_COUNT__;
    const float* ks = vox_data + HTRDR_Ks * HTRDR_SVX_OPS_COUNT__;
    const float* kext = vox_data + HTRDR_Kext * HTRDR_SVX_OPS_COUNT__;
    ASSERT(ka[HTRDR_SVX_MIN] <= ka[HTRDR_SVX_MAX]);
    ASSERT(ks[HTRDR_SVX_MIN] <= ks[HTRDR_SVX_MAX]);
    ASSERT(kext[HTRDR_SVX_MIN] <= kext[HTRDR_SVX_MAX]);

    ka_min = MMIN(ka_min, ka[HTRDR_SVX_MIN]);
    ka_max = MMAX(ka_max, ka[HTRDR_SVX_MAX]);
    ks_min = MMIN(ks_min, ks[HTRDR_SVX_MIN]);
    ks_max = MMAX(ks_max, ks[HTRDR_SVX_MAX]);
    kext_min = MMIN(kext_min, kext[HTRDR_SVX_MIN]);
    kext_max = MMAX(kext_max, kext[HTRDR_SVX_MAX]);
  }

  pflt[HTRDR_Ka   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = ka_min;
  pflt[HTRDR_Ka   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = ka_max;
  pflt[HTRDR_Ks   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = ks_min;
  pflt[HTRDR_Ks   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = ks_max;
  pflt[HTRDR_Kext * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = kext_min;
  pflt[HTRDR_Kext * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = kext_max;
}

static int
vox_challenge_merge(const void* voxels[], const size_t nvoxs, void* context)
{
  struct build_octree_context* ctx = context;
  float kext_min = FLT_MAX;
  float kext_max =-FLT_MAX;
  size_t ivox;
  ASSERT(voxels && nvoxs && context);

  FOR_EACH(ivox, 0, nvoxs) {
    const float* kext = (const float*)voxels[ivox] + HTRDR_Kext * HTRDR_SVX_OPS_COUNT__;
    ASSERT(kext[HTRDR_SVX_MIN] <= kext[HTRDR_SVX_MAX]);
    kext_min = MMIN(kext_min, kext[HTRDR_SVX_MIN]);
    kext_max = MMAX(kext_max, kext[HTRDR_SVX_MAX]);
  }
  return (kext_max - kext_min)*ctx->dst_max <= ctx->tau_threshold;
}

static res_T
setup_clouds(struct htrdr_sky* sky)
{
  struct svx_voxel_desc vox_desc = SVX_VOXEL_DESC_NULL;
  struct build_octree_context ctx;
  size_t nvoxs[3];
  double low[3];
  double upp[3];
  double sz[3];
  res_T res = RES_OK;
  ASSERT(sky);

  res = htcp_get_desc(sky->htcp, &sky->htcp_desc);
  if(res != RES_OK) {
    htrdr_log_err(sky->htrdr, "could not retrieve the HTCP descriptor.\n");
    goto error;
  }

  /* Define the number of voxels */
  nvoxs[0] = sky->htcp_desc.spatial_definition[0];
  nvoxs[1] = sky->htcp_desc.spatial_definition[1];
  nvoxs[2] = sky->htcp_desc.spatial_definition[2];

  /* Define the octree AABB exepted for the Z dimension */
  low[0] = sky->htcp_desc.lower[0];
  low[1] = sky->htcp_desc.lower[1];
  low[2] = sky->htcp_desc.lower[2];
  upp[0] = low[0] + (double)nvoxs[0] * sky->htcp_desc.vxsz_x;
  upp[1] = low[1] + (double)nvoxs[1] * sky->htcp_desc.vxsz_y;

  if(!sky->htcp_desc.irregular_z) {
    /* Regular voxel size along the Z dimension: compute its upper boundary as
     * the others dimensions */
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
    /* Allocate the svx2htcp LUT. This LUT is a regular table whose absolute
     * size is the size of a Z column in the htcp file. The size of its cells
     * is the minimal voxel size in Z of the htcp file */
    nsplits = (size_t)ceil(len_z / min_vxsz_z);
    res = darray_split_resize(&sky->svx2htcp_z, nsplits);
    if(res != RES_OK) {
      htrdr_log_err(sky->htrdr,
        "could not allocate the table mapping regular to irregular Z.\n");
      goto error;
    }
    /* Setup the svx2htcp LUT. Each LUT entry stores the index of the current Z
     * voxel in the htcp file that overlaps the entry lower bound as well as the
     * lower bound in Z of the the next htcp voxel. */
    iz2 = 0;
    upp[2] = low[2] + sky->htcp_desc.vxsz_z[iz2];
    FOR_EACH(iz, 0, nsplits) {
      const double upp_z = (double)(iz + 1) * min_vxsz_z + low[2];
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
  ctx.sky = sky;
  ctx.dst_max = sz[2];
  ctx.tau_threshold = 0.7;

  /* Setup the voxel descriptor */
  vox_desc.get = vox_get;
  vox_desc.merge = vox_merge;
  vox_desc.challenge_merge = vox_challenge_merge;
  vox_desc.context = &ctx;
  vox_desc.size = sizeof(float)
    * HTRDR_SVX_OPS_COUNT__ * HTRDR_PROPERTIES_COUNT__;

  /* Create the octree */
  res = svx_octree_create
    (sky->htrdr->svx, low, upp, nvoxs, &vox_desc, &sky->clouds);
  if(res != RES_OK) {
    htrdr_log_err(sky->htrdr,
      "could not create the octree of the cloud properties.\n");
    goto error;
  }

  /* Fetch the octree descriptor for future use */
  SVX(tree_get_desc(sky->clouds, &sky->cloud_desc));

exit:
  return res;
error:
  if(sky->clouds) SVX(tree_ref_put(sky->clouds));
  darray_split_clear(&sky->svx2htcp_z);
  goto exit;
}

static void
release_sky(ref_T* ref)
{
  struct htrdr_sky* sky;
  ASSERT(ref);
  sky = CONTAINER_OF(ref, struct htrdr_sky, ref);
  if(sky->sun) htrdr_sun_ref_put(sky->sun);
  if(sky->clouds) SVX(tree_ref_put(sky->clouds));
  if(sky->htcp) HTCP(ref_put(sky->htcp));
  if(sky->htmie) HTMIE(ref_put(sky->htmie));
  darray_split_release(&sky->svx2htcp_z);
  MEM_RM(sky->htrdr->allocator, sky);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_sky_create
  (struct htrdr* htrdr,
   struct htrdr_sun* sun,
   const char* htcp_filename,
   const char* htmie_filename,
   struct htrdr_sky** out_sky)
{
  const double band_sw[2] = {SW_WAVELENGTH_MIN, SW_WAVELENGTH_MAX};
  struct htrdr_sky* sky = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && sun && htcp_filename && htmie_filename && out_sky);

  sky = MEM_CALLOC(htrdr->allocator, 1, sizeof(*sky));
  if(!sky) {
    htrdr_log_err(htrdr, "could not allocate the sky data structure.\n");
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&sky->ref);
  htrdr_sun_ref_get(sun);
  sky->htrdr = htrdr;
  sky->sun = sun;
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

  res = htmie_create(&htrdr->logger, htrdr->allocator, htrdr->verbose, &sky->htmie);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "could not create the loader of Mie's data.\n");
    goto error;
  }

  res = htmie_load(sky->htmie, htmie_filename);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "error loading the Mie's data -- `%s'.\n",
      htmie_filename);
    goto error;
  }

  sky->Ca_avg_sw = htmie_compute_xsection_absorption_average
    (sky->htmie, band_sw, HTMIE_FILTER_LINEAR);
  sky->Cs_avg_sw = htmie_compute_xsection_scattering_average
    (sky->htmie, band_sw, HTMIE_FILTER_LINEAR);
  sky->g_avg_sw = htmie_compute_asymmetry_parameter_average
    (sky->htmie, band_sw, HTMIE_FILTER_LINEAR);
  ASSERT(sky->Ca_avg_sw > 0 && sky->Cs_avg_sw > 0 && sky->g_avg_sw > 0);

  res = setup_clouds(sky);
  if(res != RES_OK) goto error;

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

double
htrdr_sky_fetch_particle_phase_function_asymmetry_parameter
  (const struct htrdr_sky* sky, const double wavelength)
{
  ASSERT(sky && wavelength > 0);
  (void)wavelength;
  return sky->g_avg_sw;
}

double
htrdr_sky_fetch_raw_property
  (const struct htrdr_sky* sky,
   const enum htrdr_sky_property prop,
   const int components_mask, /* Combination of htrdr_sky_component_flag */
   const double wavelength, /* FIXME Unused */
   const double pos[3])
{
  size_t ivox[3];
  int comp_mask = components_mask;
  double k_particle = 0;
  double k_gaz = 0;
  ASSERT(sky && pos);
  ASSERT(comp_mask & HTRDR_ALL_COMPONENTS);
  (void)wavelength; /* FIXME */

  /* Is the position outside the clouds? */
  if(pos[0] < sky->cloud_desc.lower[0]
  || pos[1] < sky->cloud_desc.lower[1]
  || pos[2] < sky->cloud_desc.lower[2]
  || pos[0] > sky->cloud_desc.upper[0]
  || pos[1] > sky->cloud_desc.upper[1]
  || pos[2] > sky->cloud_desc.upper[2]) {
    comp_mask &= ~HTRDR_PARTICLES; /* No particle */
  }

  /* Compute the index of the voxel to fetch */
  ivox[0] = (size_t)((pos[0] - sky->cloud_desc.lower[0])/sky->htcp_desc.vxsz_x);
  ivox[1] = (size_t)((pos[1] - sky->cloud_desc.lower[1])/sky->htcp_desc.vxsz_y);
  if(!sky->htcp_desc.irregular_z) {
    /* The voxels along the Z dimension have the same size */
    ivox[2] = (size_t)((pos[2] - sky->cloud_desc.lower[2])/sky->htcp_desc.vxsz_z[0]);
  } else {
    /* Irregular voxel size along the Z dimension. Compute the index of the Z
     * position in the svx2htcp_z Look Up Table and use the LUT to define the
     * voxel index into the HTCP descripptor */
    const struct split* splits = darray_split_cdata_get(&sky->svx2htcp_z);
    const size_t n = darray_split_size_get(&sky->svx2htcp_z);
    const double sz = sky->cloud_desc.upper[2] - sky->cloud_desc.lower[2];
    const double vxsz_lut = sz / (double)n;
    const size_t ilut = (size_t)((pos[2] - sky->cloud_desc.lower[2])/vxsz_lut);
    ivox[2] = splits[ilut].index + (pos[2] >= splits[ilut].height);
  }

  if(comp_mask & HTRDR_PARTICLES) {
    double rho_da = 0; /* Dry air density */
    double rct = 0; /* #droplets in kg of water per kg of dry air */
    double ql = 0; /* Droplet density In kg.m^-3 */
    double Ca = 0; /* Massic absorption cross section in m^2.kg^-1 */
    double Cs = 0; /* Massic scattering cross section in m^2.kg^-1 */

    /* Compute he dry air density */
    rho_da = cloud_dry_air_density(&sky->htcp_desc, ivox);

    /* Compute the droplet density */
    rct = htcp_desc_RCT_at(&sky->htcp_desc, ivox[0], ivox[1], ivox[2], 0);
    ql = rho_da * rct;

    /* FIXME do not use the wavelength yet. Simply use the average cross
     * section on the while short wave range */
#if 1
    if(prop == HTRDR_Ka || prop == HTRDR_Kext) Ca = sky->Ca_avg_sw;
    if(prop == HTRDR_Ks || prop == HTRDR_Kext) Cs = sky->Cs_avg_sw;
#else
    if(prop == HTRDR_Ks || prop == HTRDR_Kext) {
      Cs = htmie_fetch_xsection_scattering
        (sky->htmie, wavelength, HTMIE_FILTER_LINEAR);
    }
    if(prop == HTRDR_Ka || prop == HTRDR_Kext) {
      Ca = htmie_fetch_xsection_absorption
        (sky->htmie, wavelength, HTMIE_FILTER_LINEAR);
    }
#endif
    k_particle = ql*(Ca + Cs);
  }
  if(comp_mask & HTRDR_GAS) { /* TODO not implemented yet */ }
  return k_particle + k_gaz;
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

  octree_data_init(sky, &data);
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
  fprintf(stream, "SCALARS Kext double 2\n");
  fprintf(stream, "LOOKUP_TABLE default\n");
  FOR_EACH(i, 0, ncells) {
    fprintf(stream, "%g %g\n", leaf_data[i*2+0], leaf_data[i*2+1]);
  }
  octree_data_release(&data);
  return RES_OK;
}

double
htrdr_sky_fetch_svx_property
  (const struct htrdr_sky* sky,
   const enum htrdr_sky_property prop,
   const enum htrdr_svx_op op,
   const int components_mask, /* Combination of htrdr_sky_component_flag */
   const double wavelength,
   const double pos[3])
{
  struct svx_voxel voxel = SVX_VOXEL_NULL;
  int comp_mask = components_mask;
  ASSERT(sky && pos);
  ASSERT(comp_mask & HTRDR_ALL_COMPONENTS);

  /* Is the position outside the clouds? */
  if(pos[0] < sky->cloud_desc.lower[0]
  || pos[1] < sky->cloud_desc.lower[1]
  || pos[2] < sky->cloud_desc.lower[2]
  || pos[0] > sky->cloud_desc.upper[0]
  || pos[1] > sky->cloud_desc.upper[1]
  || pos[2] > sky->cloud_desc.upper[2]) {
    comp_mask &= ~HTRDR_PARTICLES; /* No particle */
  }

  SVX(tree_at(sky->clouds, pos, NULL, NULL, &voxel));
  return htrdr_sky_fetch_svx_voxel_property
      (sky, prop, op, comp_mask, wavelength, &voxel);
}

double
htrdr_sky_fetch_svx_voxel_property
  (const struct htrdr_sky* sky,
   const enum htrdr_sky_property prop,
   const enum htrdr_svx_op op,
   const int components_mask,
   const double wavelength, /* FIXME Unused */
   const struct svx_voxel* voxel)
{
  const float* pflt = NULL;
  int comp_mask = components_mask;
  double a, b, data;
  double gas = 0;
  double particle = 0;
  ASSERT(sky && voxel);
  ASSERT((unsigned)prop < HTRDR_PROPERTIES_COUNT__);
  ASSERT((unsigned)op < HTRDR_SVX_OPS_COUNT__);
  (void)sky, (void)wavelength;

  pflt = voxel->data;

  if(comp_mask) {
    particle = pflt[prop * HTRDR_SVX_OPS_COUNT__ + op];
  }
  if(comp_mask & HTRDR_GAS) { comp_mask &= ~HTRDR_GAS; /* TODO not implemented yet */ }

  switch(op) {
    case HTRDR_SVX_MIN:
      a = comp_mask & HTRDR_PARTICLES ? particle : DBL_MAX;
      b = comp_mask & HTRDR_GAS ? gas : DBL_MAX;
      data = MMIN(a, b);
      break;
    case HTRDR_SVX_MAX:
      a = comp_mask & HTRDR_PARTICLES ? particle : -DBL_MAX;
      b = comp_mask & HTRDR_GAS ? gas : -DBL_MAX;
      data = MMAX(a, b);
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
  return data;
}

