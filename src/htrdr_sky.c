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
#include "htrdr_grid.h"
#include "htrdr_sky.h"
#include "htrdr_sun.h"

#include <star/ssp.h>
#include <star/svx.h>
#include <high_tune/htcp.h>
#include <high_tune/htgop.h>
#include <high_tune/htmie.h>

#include <rsys/dynamic_array_double.h>
#include <rsys/dynamic_array_size_t.h>
#include <rsys/hash_table.h>
#include <rsys/ref_count.h>

#include <libgen.h>
#include <math.h>
#include <omp.h>

#define DRY_AIR_MOLAR_MASS 0.0289644 /* In kg.mol^-1 */
#define H2O_MOLAR_MASS 0.01801528 /* In kg.mol^-1 */
#define GAS_CONSTANT 8.3144598 /* In kg.m^2.s^-2.mol^-1.K */

#define NFLOATS_PER_COMPONENT (HTRDR_SVX_OPS_COUNT__ * HTRDR_PROPERTIES_COUNT__)

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
  size_t iband; /* Index of the band that overlaps the CIE XYZ color space */

  /* Precomputed voxel data of the finest level. May be NULL <=> compute the
   * voxel data at runtime. */
  struct htrdr_grid* grid;
};

#define BUILD_OCTREE_CONTEXT_NULL__ { NULL, 0, 0, 0, NULL }
static const struct build_octree_context BUILD_OCTREE_CONTEXT_NULL =
  BUILD_OCTREE_CONTEXT_NULL__;

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

/* Temporary data structure used to dump the octree data into a VTK file */
struct octree_data {
  struct htable_vertex vertex2id; /* Map a coordinate to its vertex id */
  struct darray_double vertices; /* Array of unique vertices */
  struct darray_double data; /* List of registered leaf data */
  struct darray_size_t cells; /* Ids of the cell vertices */
  size_t iband; /* Index of the band that overlaps the CIE XYZ color space */
  size_t iquad; /* Index of the quadrature point into the band */
  const struct htrdr_sky* sky;
};

/* Properties of a short wave spectral band */
struct sw_band_prop {
  /* Average cross section in the band */
  double Ca_avg; /* Absorption cross section */
  double Cs_avg; /* Scattering cross section */

  /* Average asymmetry parameter the band */
  double g_avg;
};

/* Encompass the hierarchical data structure of the cloud data and its
 * associated descriptor */
struct cloud {
  struct svx_tree* octree;
  struct svx_tree_desc octree_desc;
};

struct htrdr_sky {
  struct cloud* clouds; /* Per sw_band cloud data structure */

  struct htrdr_sun* sun; /* Sun attached to the sky */

  /* Loaders of... */
  struct htcp* htcp;   /* ... Cloud properties */
  struct htgop* htgop; /* ... Gas optical properties */
  struct htmie* htmie; /* ... Mie's data */

  struct htcp_desc htcp_desc;

  /* Map the index in Z from the regular SVX to the irregular HTCP data */
  struct darray_split svx2htcp_z;

  /* Ids and optical properties of the short wave spectral bands loaded by
   * HTGOP and that overlap the CIE XYZ color space */
  size_t sw_bands_range[2];
  struct sw_band_prop* sw_bands;

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
  ASSERT(desc && ivox);
  P = htcp_desc_PABST_at(desc, ivox[0], ivox[1], ivox[2], 0/*time*/);
  T = htcp_desc_T_at(desc, ivox[0], ivox[1], ivox[2], 0/*time*/);
  return (P*DRY_AIR_MOLAR_MASS)/(T*GAS_CONSTANT);
}

/* Compute the water molar fraction */
static FINLINE double
cloud_water_vapor_molar_fraction
  (const struct htcp_desc* desc,
   const size_t ivox[3])
{
  double rvt = 0;
  ASSERT(desc && ivox);
  rvt = htcp_desc_RVT_at(desc, ivox[0], ivox[1], ivox[2], 0/*time*/);
  return rvt / (rvt + H2O_MOLAR_MASS/DRY_AIR_MOLAR_MASS);
}

static INLINE void
octree_data_init
  (const struct htrdr_sky* sky,
   const size_t iband,
   const size_t iquad,
   struct octree_data* data)
{
  ASSERT(data);
  ASSERT(iband >= sky->sw_bands_range[0]);
  ASSERT(iband <= sky->sw_bands_range[1]);
  (void)iquad;
  htable_vertex_init(sky->htrdr->allocator, &data->vertex2id);
  darray_double_init(sky->htrdr->allocator, &data->vertices);
  darray_double_init(sky->htrdr->allocator, &data->data);
  darray_size_t_init(sky->htrdr->allocator, &data->cells);
  data->sky = sky;
  data->iband = iband;
  data->iquad = iquad;
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
  kext_max = htrdr_sky_fetch_svx_voxel_property(ctx->sky, HTRDR_Kext,
    HTRDR_SVX_MAX, HTRDR_ALL_COMPONENTS, ctx->iband, ctx->iquad, leaf);
  kext_min = htrdr_sky_fetch_svx_voxel_property(ctx->sky, HTRDR_Kext,
    HTRDR_SVX_MIN, HTRDR_ALL_COMPONENTS, ctx->iband, ctx->iquad, leaf);
  CHK(RES_OK == darray_double_push_back(&ctx->data, &kext_min));
  CHK(RES_OK == darray_double_push_back(&ctx->data, &kext_max));
}

static void
vox_get_particle
  (const size_t xyz[3],
   float particle[],
   const struct build_octree_context* ctx)
{
  double rct;
  double ka, ks, kext;
  double ka_min, ka_max;
  double ks_min, ks_max;
  double kext_min, kext_max;
  double rho_da; /* Dry air density */
  double Ca_avg;
  double Cs_avg;
  size_t i;
  ASSERT(xyz && particle && ctx);

  i = ctx->iband - ctx->sky->sw_bands_range[0];
  /* Fetch the optical properties of the spectral band */
  Ca_avg = ctx->sky->sw_bands[i].Ca_avg;
  Cs_avg = ctx->sky->sw_bands[i].Cs_avg;

  if(!ctx->sky->htcp_desc.irregular_z) {
    rho_da = cloud_dry_air_density(&ctx->sky->htcp_desc, xyz);
    rct = htcp_desc_RCT_at(&ctx->sky->htcp_desc, xyz[0], xyz[1], xyz[2], 0);
    ka_min = ka_max = ka = Ca_avg * rho_da * rct;
    ks_min = ks_max = ks = Cs_avg * rho_da * rct;
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

    ka_min = ka_max = ka = Ca_avg * rho_da * rct;
    ks_min = ks_max = ks = Cs_avg * rho_da * rct;
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
      ka = Ca_avg * rho_da * rct;
      ks = Cs_avg * rho_da * rct;
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

  particle[HTRDR_Ka  *HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = (float)ka_min;
  particle[HTRDR_Ka  *HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = (float)ka_max;
  particle[HTRDR_Ks  *HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = (float)ks_min;
  particle[HTRDR_Ks  *HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = (float)ks_max;
  particle[HTRDR_Kext*HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = (float)kext_min;
  particle[HTRDR_Kext*HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = (float)kext_max;
}

static void
vox_get_gas
  (const size_t xyz[3],
   float gas[],
   const struct build_octree_context* ctx)
{
  struct htgop_layer layer;
  struct htgop_layer_sw_spectral_interval band;
  size_t ilayer;
  size_t layer_range[2];
  size_t quad_range[2];
  double x_h2o_range[2];
  double lower[3], upper[3]; /* AABB */
  double ka[2] = {DBL_MAX, -DBL_MAX};
  double ks[2] = {DBL_MAX, -DBL_MAX};
  double kext[2] = {DBL_MAX, -DBL_MAX};
  ASSERT(xyz && gas && ctx);

  /* Define the xH2O range from the LES data */
  if(!ctx->sky->htcp_desc.irregular_z) { /* 1 LES voxel <=> 1 SVX voxel */
    double x_h2o;
    x_h2o = cloud_water_vapor_molar_fraction(&ctx->sky->htcp_desc, xyz);
    x_h2o_range[0] = x_h2o_range[1] = x_h2o;
  } else { /* A SVX voxel can be overlapped by 2 LES voxels */
    size_t ivox[3];
    size_t ivox_next;
    ASSERT(xyz[2] < darray_split_size_get(&ctx->sky->svx2htcp_z));

    ivox[0] = xyz[0];
    ivox[1] = xyz[1];
    ivox[2] = darray_split_cdata_get(&ctx->sky->svx2htcp_z)[xyz[2]].index;

    x_h2o_range[0] = cloud_water_vapor_molar_fraction(&ctx->sky->htcp_desc, ivox);

    /* Define if the SVX voxel is overlapped by 2 HTCP voxels */
    ivox_next = xyz[2] + 1 < darray_split_size_get(&ctx->sky->svx2htcp_z)
      ? darray_split_cdata_get(&ctx->sky->svx2htcp_z)[xyz[2] + 1].index
      : ivox[2];

    if(ivox_next == ivox[2]) { /* No overlap */
      x_h2o_range[1] = x_h2o_range[0];
    } else { /* Overlap */
      ASSERT(ivox[2] < ivox_next);
      ivox[2] = ivox_next;
      x_h2o_range[1] =
        cloud_water_vapor_molar_fraction(&ctx->sky->htcp_desc, ivox);
      if(x_h2o_range[0] > x_h2o_range[1])
        SWAP(double, x_h2o_range[0], x_h2o_range[1]);
    }
  }

  /* Retrieve the range of atmospheric layers that overlap the SVX voxel */
  htcp_desc_get_voxel_aabb
    (&ctx->sky->htcp_desc, xyz[0], xyz[1], xyz[2], lower, upper);
  HTGOP(position_to_layer_id(ctx->sky->htgop, lower[2], &layer_range[0]));
  HTGOP(position_to_layer_id(ctx->sky->htgop, upper[2], &layer_range[1]));

  /* For each atmospheric layer that overlaps the SVX voxel ... */
  FOR_EACH(ilayer, layer_range[0], layer_range[1]+1) {
    double k[2];

    HTGOP(get_layer(ctx->sky->htgop, ilayer, &layer));

    /* ... retrieve the considered spectral interval */
    HTGOP(layer_get_sw_spectral_interval(&layer, ctx->iband, &band));
    quad_range[0] = 0;
    quad_range[1] = band.quadrature_length-1;

    /* ... and compute the radiative properties and upd their bounds */
    HTGOP(layer_sw_spectral_interval_quadpoints_get_ka_bounds
      (&band, quad_range, x_h2o_range, k));
    ka[0] = MMIN(ka[0], k[0]);
    ka[1] = MMAX(ka[1], k[1]);
    HTGOP(layer_sw_spectral_interval_quadpoints_get_ks_bounds
      (&band, quad_range, x_h2o_range, k));
    ks[0] = MMIN(ks[0], k[0]);
    ks[1] = MMAX(ks[1], k[1]);
    HTGOP(layer_sw_spectral_interval_quadpoints_get_kext_bounds
      (&band, quad_range, x_h2o_range, k));
    kext[0] = MMIN(kext[0], k[0]);
    kext[1] = MMAX(kext[1], k[1]);
  }

  /* Ensure that the single precision bounds include their double precision
   * version. */
  if(ka[0] != (float)ka[0]) ka[0] = nextafterf((float)ka[0],-FLT_MAX);
  if(ka[1] != (float)ka[1]) ka[1] = nextafterf((float)ka[1], FLT_MAX);
  if(ks[0] != (float)ks[0]) ks[0] = nextafterf((float)ks[0],-FLT_MAX);
  if(ks[1] != (float)ks[1]) ks[1] = nextafterf((float)ks[1], FLT_MAX);
  if(kext[0] != (float)kext[0]) kext[0] = nextafterf((float)kext[0],-FLT_MAX);
  if(kext[1] != (float)kext[1]) kext[1] = nextafterf((float)kext[1], FLT_MAX);

  gas[HTRDR_Ka  *HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = (float)ka[0];
  gas[HTRDR_Ka  *HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = (float)ka[1];
  gas[HTRDR_Ks  *HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = (float)ks[0];
  gas[HTRDR_Ks  *HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = (float)ks[1];
  gas[HTRDR_Kext*HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = (float)kext[0];
  gas[HTRDR_Kext*HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = (float)kext[1];
}

static void
vox_get(const size_t xyz[3], void* dst, void* context)
{
  struct build_octree_context* ctx = context;
  ASSERT(context);

  if(ctx->grid) { /* Fetch voxel data from precomputed grid */
    const float* vox_data = htrdr_grid_at(ctx->grid, xyz);
    memcpy(dst, vox_data, NFLOATS_PER_COMPONENT*2*sizeof(float));
  } else {
    /* No precomputed grid. Compute the voxel data at runtime */
    float* par = (float*)dst + 0*NFLOATS_PER_COMPONENT; /* Particles properties */
    float* gas = (float*)dst + 1*NFLOATS_PER_COMPONENT; /* Gas properties */
    vox_get_particle(xyz, par, ctx);
    vox_get_gas(xyz, gas, ctx);
  }
}

static INLINE void
vox_merge_component
  (float* comp, const float* voxels[], const size_t off, const size_t nvoxs)
{
  float ka_min = FLT_MAX;
  float ka_max =-FLT_MAX;
  float ks_min = FLT_MAX;
  float ks_max =-FLT_MAX;
  float kext_min = FLT_MAX;
  float kext_max =-FLT_MAX;
  size_t ivox;
  ASSERT(comp && voxels && nvoxs);

  FOR_EACH(ivox, 0, nvoxs) {
    const float* ka = voxels[ivox] + off + HTRDR_Ka * HTRDR_SVX_OPS_COUNT__;
    const float* ks = voxels[ivox] + off + HTRDR_Ks * HTRDR_SVX_OPS_COUNT__;
    const float* kext = voxels[ivox] + off + HTRDR_Kext * HTRDR_SVX_OPS_COUNT__;
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

  comp[HTRDR_Ka   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = ka_min;
  comp[HTRDR_Ka   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = ka_max;
  comp[HTRDR_Ks   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = ks_min;
  comp[HTRDR_Ks   * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = ks_max;
  comp[HTRDR_Kext * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MIN] = kext_min;
  comp[HTRDR_Kext * HTRDR_SVX_OPS_COUNT__ + HTRDR_SVX_MAX] = kext_max;
}

static void
vox_merge(void* dst, const void* voxels[], const size_t nvoxs, void* context)
{
  float* par = (float*)dst + 0*NFLOATS_PER_COMPONENT; /* Particles properties */
  float* gas = (float*)dst + 1*NFLOATS_PER_COMPONENT; /* Gas properties */
  ASSERT(dst && voxels);
  (void) context;
  vox_merge_component(par, (const float**)voxels, 0*NFLOATS_PER_COMPONENT, nvoxs);
  vox_merge_component(gas, (const float**)voxels, 1*NFLOATS_PER_COMPONENT, nvoxs);
}

static INLINE int
vox_challenge_merge_component
  (const float* voxels[],
   const size_t nvoxs,
   const size_t off,
   struct build_octree_context* ctx)
{
  float kext_min = FLT_MAX;
  float kext_max =-FLT_MAX;
  size_t ivox;
  ASSERT(voxels && nvoxs && ctx);

  FOR_EACH(ivox, 0, nvoxs) {
    const float* kext = voxels[ivox] + off + HTRDR_Kext * HTRDR_SVX_OPS_COUNT__;
    ASSERT(kext[HTRDR_SVX_MIN] <= kext[HTRDR_SVX_MAX]);
    kext_min = MMIN(kext_min, kext[HTRDR_SVX_MIN]);
    kext_max = MMAX(kext_max, kext[HTRDR_SVX_MAX]);
  }
  return (kext_max - kext_min)*ctx->dst_max <= ctx->tau_threshold;
}

static int
vox_challenge_merge
  (const void* voxels[], const size_t nvoxs, void* ctx)
{
  const float** vox = (const float**)voxels;
  ASSERT(voxels);
  return vox_challenge_merge_component(vox, nvoxs, 0*NFLOATS_PER_COMPONENT, ctx)
      && vox_challenge_merge_component(vox, nvoxs, 1*NFLOATS_PER_COMPONENT, ctx);
}

/* Create/load a grid of cloud data used by SVX to build the octree. The grid
 * is saved in the directory where htrdr is run with a name generated from the
 * "htcp_filename" path. If a grid with the same name exists, the function
 * tries to load it except if the force_update flag is set. Even though the
 * grid is loaded from disk, the function will recompute and store it if the
 * definition of the loaded grid is different from the submitted definition. */
static res_T
setup_cloud_grid
  (struct htrdr_sky* sky,
   const size_t definition[3],
   const size_t iband,
   const char* htcp_filename,
   const int force_update,
   struct htrdr_grid** out_grid)
{
  struct htrdr_grid* grid = NULL;
  struct str str;
  struct build_octree_context ctx = BUILD_OCTREE_CONTEXT_NULL;
  size_t sizeof_cell;
  size_t ncells;
  uint64_t mcode;
  uint64_t mcode_max;
  char buf[16];
  size_t progress = 0;
  ATOMIC ncells_computed = 0;
  res_T res = RES_OK;
  ASSERT(sky && definition && htcp_filename && out_grid);
  ASSERT(definition[0] && definition[1] && definition[2]);
  ASSERT(iband >= sky->sw_bands_range[0] && iband <= sky->sw_bands_range[1]);

  CHK((size_t)snprintf(buf, sizeof(buf), ".%lu", iband) < sizeof(buf));

  /* Build the grid name */
  str_init(sky->htrdr->allocator, &str);
  CHK(RES_OK == str_set(&str, htcp_filename));
  CHK(RES_OK == str_set(&str, basename(str_get(&str))));
  CHK(RES_OK == str_insert(&str, 0, ".htrdr_"));
  CHK(RES_OK == str_append(&str, ".grid"));
  CHK(RES_OK == str_append(&str, buf));

  if(!force_update) {
    /* Try to open the saved grid */
    res = htrdr_grid_open(sky->htrdr, str_cget(&str), &grid);
    if(res != RES_OK) {
      htrdr_log_warn(sky->htrdr, "cannot open the grid `%s'.\n", str_cget(&str));
    } else {
      size_t grid_def[3];

      /* Check that the definition is the loaded grid is the same of the
       * submitted grid definition */
      htrdr_grid_get_definition(grid, grid_def);
      if(grid_def[0] == definition[0]
      && grid_def[1] == definition[1]
      && grid_def[2] == definition[2]) {
        htrdr_log(sky->htrdr,
          "Use the precomputed grid `%s'.\n", str_cget(&str));
        goto exit; /* No more work to do. The loaded data seems valid */
      }

      /* The grid is no more valid. Update it! */
      htrdr_grid_ref_put(grid);
      grid = NULL;
    }
  }

  sizeof_cell = NFLOATS_PER_COMPONENT * 2/*gas & particle*/ * sizeof(float);

  htrdr_log(sky->htrdr, "Compute the grid `%s'.\n", str_cget(&str));

  res = htrdr_grid_create
    (sky->htrdr, definition, sizeof_cell, str_cget(&str), 1, &grid);
  if(res != RES_OK) goto error;

  ctx.sky = sky;
  ctx.dst_max = DBL_MAX; /* Unused for grid construction */
  ctx.tau_threshold = DBL_MAX; /* Unused for grid construction */
  ctx.iband = iband;

  mcode_max = MMAX(MMAX(definition[0], definition[1]), definition[2]);
  mcode_max = round_up_pow2(mcode_max);
  mcode_max = mcode_max*mcode_max*mcode_max;

  ncells = definition[0] * definition[1] * definition[2];

  fprintf(stderr, "Generating cloud grid %lu: %3u%%", iband, 0);
  fflush(stderr);

  omp_set_num_threads((int)sky->htrdr->nthreads);

  #pragma omp parallel for
  for(mcode=0; mcode<mcode_max; ++mcode) {
    size_t xyz[3];
    size_t pcent;
    size_t n;
    if((xyz[0] = morton3D_decode_u21(mcode >> 2)) >= definition[0]) continue;
    if((xyz[1] = morton3D_decode_u21(mcode >> 1)) >= definition[1]) continue;
    if((xyz[2] = morton3D_decode_u21(mcode >> 0)) >= definition[2]) continue;
    vox_get(xyz, htrdr_grid_at_mcode(grid, mcode), &ctx);

    n = (size_t)ATOMIC_INCR(&ncells_computed);
    pcent = n * 100 / ncells;
    #pragma omp critical 
    if(pcent > progress) {
      progress = pcent;
      fprintf(stderr, "%c[2K\rGenerating cloud grid %lu: %3u%%",
        27, iband, (unsigned)pcent);
      fflush(stderr);
    }
  }
  fprintf(stderr, "\n");

exit:
  *out_grid = grid;
  str_release(&str);
  return res;
error:
  if(grid) {
    htrdr_grid_ref_put(grid);
    grid = NULL;
  }
  goto exit;
}

static res_T
setup_clouds
  (struct htrdr_sky* sky,
   const char* htcp_filename,
   const int force_cache_update)
{
  struct svx_voxel_desc vox_desc = SVX_VOXEL_DESC_NULL;
  struct build_octree_context ctx = BUILD_OCTREE_CONTEXT_NULL;
  size_t nvoxs[3];
  double low[3];
  double upp[3];
  double sz[3];
  size_t nbands;
  size_t i;
  res_T res = RES_OK;
  ASSERT(sky && sky->sw_bands);

  res = htcp_get_desc(sky->htcp, &sky->htcp_desc);
  if(res != RES_OK) {
    htrdr_log_err(sky->htrdr, "could not retrieve the HTCP descriptor.\n");
    goto error;
  }

  /* Define the number of voxels */
  nvoxs[0] = sky->htcp_desc.spatial_definition[0];
  nvoxs[1] = sky->htcp_desc.spatial_definition[1];
  nvoxs[2] = sky->htcp_desc.spatial_definition[2];

  /* Define the octree AABB excepted for the Z dimension */
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
  ctx.tau_threshold = 0.5;

  /* Setup the voxel descriptor */
  vox_desc.get = vox_get;
  vox_desc.merge = vox_merge;
  vox_desc.challenge_merge = vox_challenge_merge;
  vox_desc.context = &ctx;
  vox_desc.size = sizeof(float)
    * HTRDR_SVX_OPS_COUNT__ * HTRDR_PROPERTIES_COUNT__ * 2/*Gas & particles*/;

  /* Create as many cloud data structure than considered SW spectral bands */
  nbands = htrdr_sky_get_sw_spectral_bands_count(sky);
  sky->clouds = MEM_CALLOC(sky->htrdr->allocator, nbands, sizeof(*sky->clouds));
  if(!sky->clouds) {
    htrdr_log_err(sky->htrdr,
      "could not create the list of per SW band cloud data structure.\n");
    res = RES_MEM_ERR;
    goto error;
  }

  FOR_EACH(i, 0, nbands) {
    ctx.iband = i + sky->sw_bands_range[0];

    /* Compute grid of voxels data */
    res = setup_cloud_grid(sky, nvoxs, ctx.iband, htcp_filename,
      force_cache_update, &ctx.grid);
    if(res != RES_OK) goto error;

    /* Create the octree */
    res = svx_octree_create
      (sky->htrdr->svx, low, upp, nvoxs, &vox_desc, &sky->clouds[i].octree);
    if(res != RES_OK) {
      htrdr_log_err(sky->htrdr, "could not create the octree of the cloud "
        "properties for the band %lu.\n", (unsigned long)ctx.iband);
      goto error;
    }

    /* Fetch the octree descriptor for future use */
    SVX(tree_get_desc
      (sky->clouds[i].octree, &sky->clouds[i].octree_desc));

    if(ctx.grid) {
      htrdr_grid_ref_put(ctx.grid);
      ctx.grid = NULL;
    }
  }

exit:
  return res;
error:
  if(ctx.grid) htrdr_grid_ref_put(ctx.grid);
  if(sky->clouds) {
    FOR_EACH(i, 0, nbands) {
      if(sky->clouds[i].octree) {
        SVX(tree_ref_put(sky->clouds[i].octree));
        sky->clouds[i].octree = NULL;
      }
    }
    MEM_RM(sky->htrdr->allocator, sky->clouds);
  }
  darray_split_clear(&sky->svx2htcp_z);
  goto exit;
}

static res_T
setup_sw_bands_properties(struct htrdr_sky* sky)
{
  res_T res = RES_OK;
  size_t nbands;
  size_t i;
  ASSERT(sky);

  res = htgop_get_sw_spectral_intervals_CIE_XYZ(sky->htgop, sky->sw_bands_range);
  if(res != RES_OK) goto error;

  nbands = htrdr_sky_get_sw_spectral_bands_count(sky);
  ASSERT(nbands);
  sky->sw_bands = MEM_CALLOC
    (sky->htrdr->allocator, nbands, sizeof(*sky->sw_bands));
  if(!sky->sw_bands) {
    htrdr_log_err(sky->htrdr,
      "could not allocate the list of SW band properties.\n");
    res = RES_MEM_ERR;
    goto error;
  }

  FOR_EACH(i, 0, nbands) {
    struct htgop_spectral_interval band;
    double band_wlens[2];

    HTGOP(get_sw_spectral_interval
      (sky->htgop, i+ sky->sw_bands_range[0], &band));
    band_wlens[0] = wavenumber_to_wavelength(band.wave_numbers[1]);
    band_wlens[1] = wavenumber_to_wavelength(band.wave_numbers[0]);
    ASSERT(band_wlens[0] < band_wlens[1]);

    sky->sw_bands[i].Ca_avg = htmie_compute_xsection_absorption_average
      (sky->htmie, band_wlens, HTMIE_FILTER_LINEAR);
    sky->sw_bands[i].Cs_avg = htmie_compute_xsection_scattering_average
      (sky->htmie, band_wlens, HTMIE_FILTER_LINEAR);
    sky->sw_bands[i].g_avg = htmie_compute_asymmetry_parameter_average
      (sky->htmie, band_wlens, HTMIE_FILTER_LINEAR);
    ASSERT(sky->sw_bands[i].Ca_avg > 0);
    ASSERT(sky->sw_bands[i].Cs_avg > 0);
    ASSERT(sky->sw_bands[i].g_avg > 0);
  }

exit:
  return res;
error:
  if(sky->sw_bands) {
    MEM_RM(sky->htrdr->allocator, sky->sw_bands);
    sky->sw_bands = NULL;
  }
  goto exit;
}

static void
sample_sw_spectral_data
  (struct htgop* htgop,
   struct ssp_rng* rng,
   res_T (*sample_sw_band)(const struct htgop*, const double, size_t*),
   size_t* ispectral_band,
   size_t* iquadrature_pt)
{
  struct htgop_spectral_interval specint;
  double r1, r2;
  res_T res = RES_OK;
  ASSERT(htgop && rng && sample_sw_band && ispectral_band && iquadrature_pt);
  ASSERT(ispectral_band && iquadrature_pt);
  (void)res;
  r1 = ssp_rng_canonical(rng);
  r2 = ssp_rng_canonical(rng);
  res = sample_sw_band(htgop, r1, ispectral_band);
  ASSERT(res == RES_OK);
  HTGOP(get_sw_spectral_interval(htgop, *ispectral_band, &specint));
  HTGOP(spectral_interval_sample_quadrature(&specint, r2, iquadrature_pt));
}

static void
release_sky(ref_T* ref)
{
  struct htrdr_sky* sky;
  ASSERT(ref);
  sky = CONTAINER_OF(ref, struct htrdr_sky, ref);
  if(sky->sun) htrdr_sun_ref_put(sky->sun);
  if(sky->htcp) HTCP(ref_put(sky->htcp));
  if(sky->htgop) HTGOP(ref_put(sky->htgop));
  if(sky->htmie) HTMIE(ref_put(sky->htmie));
  if(sky->sw_bands) MEM_RM(sky->htrdr->allocator, sky->sw_bands);
  if(sky->clouds) {
    const size_t nbands = htrdr_sky_get_sw_spectral_bands_count(sky);
    size_t i;
    FOR_EACH(i, 0, nbands) {
      if(sky->clouds[i].octree) {
        SVX(tree_ref_put(sky->clouds[i].octree));
        sky->clouds[i].octree = NULL;
      }
    }
    MEM_RM(sky->htrdr->allocator, sky->clouds);
  }
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
   const char* htgop_filename,
   const char* htmie_filename,
   struct htrdr_sky** out_sky)
{
  struct htrdr_sky* sky = NULL;
  int htcp_upd = 1;
  int htmie_upd = 1;
  int htgop_upd = 1;
  int force_upd = 1;
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

  /* Load clouds properties */
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

  /* Load MIE data */
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

  /* Load the gas optical properties */
  res = htgop_create
    (&htrdr->logger, htrdr->allocator, htrdr->verbose, &sky->htgop);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "could not create the loader of the gas optical properties.\n");
    goto error;
  }
  res = htgop_load(sky->htgop, htgop_filename);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "error loading the gas optical properties -- `%s'.\n",
      htgop_filename);
    goto error;
  }

  res = setup_sw_bands_properties(sky);
  if(res != RES_OK) goto error;

  /* Define if the cached grid data must be updated */
  res = is_file_updated(sky->htrdr, htcp_filename, &htcp_upd);
  if(res != RES_OK) goto error;
  res = is_file_updated(sky->htrdr, htmie_filename, &htmie_upd);
  if(res != RES_OK) goto error;
  res = is_file_updated(sky->htrdr, htgop_filename, &htgop_upd);
  if(res != RES_OK) goto error;
  force_upd = htcp_upd || htmie_upd || htgop_upd;

  res = setup_clouds(sky, htcp_filename, force_upd);
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
  (const struct htrdr_sky* sky,
   const size_t ispectral_band,
   const size_t iquad)
{
  size_t i;
  ASSERT(sky);
  ASSERT(ispectral_band >= sky->sw_bands_range[0]);
  ASSERT(ispectral_band <= sky->sw_bands_range[1]);
  (void)iquad;
  i = ispectral_band - sky->sw_bands_range[0];
  return sky->sw_bands[i].g_avg;
}

double
htrdr_sky_fetch_raw_property
  (const struct htrdr_sky* sky,
   const enum htrdr_sky_property prop,
   const int components_mask, /* Combination of htrdr_sky_component_flag */
   const size_t iband, /* Index of the spectral band */
   const size_t iquad, /* Index of the quadrature point in the spectral band */
   const double pos[3])
{
  size_t ivox[3];
  size_t i;
  const struct svx_tree_desc* cloud_desc;
  int comp_mask = components_mask;
  double k_particle = 0;
  double k_gas = 0;
  double k = 0;
  ASSERT(sky && pos);
  ASSERT(iband >= sky->sw_bands_range[0]);
  ASSERT(iband <= sky->sw_bands_range[1]);
  ASSERT(comp_mask & HTRDR_ALL_COMPONENTS);

  i = iband - sky->sw_bands_range[0];
  cloud_desc = &sky->clouds[i].octree_desc;

  /* Is the position outside the clouds? */
  if(pos[0] < cloud_desc->lower[0]
  || pos[1] < cloud_desc->lower[1]
  || pos[2] < cloud_desc->lower[2]
  || pos[0] > cloud_desc->upper[0]
  || pos[1] > cloud_desc->upper[1]
  || pos[2] > cloud_desc->upper[2]) {
    comp_mask &= ~HTRDR_PARTICLES; /* No particle */
  }

  /* Compute the index of the voxel to fetch */
  ivox[0] = (size_t)((pos[0] - cloud_desc->lower[0])/sky->htcp_desc.vxsz_x);
  ivox[1] = (size_t)((pos[1] - cloud_desc->lower[1])/sky->htcp_desc.vxsz_y);
  if(!sky->htcp_desc.irregular_z) {
    /* The voxels along the Z dimension have the same size */
    ivox[2] = (size_t)((pos[2] - cloud_desc->lower[2])/sky->htcp_desc.vxsz_z[0]);
  } else {
    /* Irregular voxel size along the Z dimension. Compute the index of the Z
     * position in the svx2htcp_z Look Up Table and use the LUT to define the
     * voxel index into the HTCP descripptor */
    const struct split* splits = darray_split_cdata_get(&sky->svx2htcp_z);
    const size_t n = darray_split_size_get(&sky->svx2htcp_z);
    const double sz = cloud_desc->upper[2] - cloud_desc->lower[2];
    const double vxsz_lut = sz / (double)n;
    const size_t ilut = (size_t)((pos[2] - cloud_desc->lower[2])/vxsz_lut);
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

    /* Use the average cross section of the current spectral band */
    if(prop == HTRDR_Ka || prop == HTRDR_Kext) Ca = sky->sw_bands[i].Ca_avg;
    if(prop == HTRDR_Ks || prop == HTRDR_Kext) Cs = sky->sw_bands[i].Cs_avg;

    k_particle = ql*(Ca + Cs);
  }

  if(comp_mask & HTRDR_GAS) {
    struct htgop_layer layer;
    struct htgop_layer_sw_spectral_interval band;
    struct htgop_layer_sw_spectral_interval_tab tab;
    const double x_h2o = cloud_water_vapor_molar_fraction(&sky->htcp_desc, ivox);
    size_t ilayer = 0;

    /* Retrieve the quadrature point into the spectral band of the layer into
     * which `pos' lies */
    HTGOP(position_to_layer_id(sky->htgop, pos[2], &ilayer));
    HTGOP(get_layer(sky->htgop, ilayer, &layer));
    HTGOP(layer_get_sw_spectral_interval(&layer, iband, &band));
    HTGOP(layer_sw_spectral_interval_get_tab(&band, iquad, &tab));

    /* Fetch the optical properties wrt x_h2o */
    switch(prop) {
      case HTRDR_Ka:
        HTGOP(layer_sw_spectral_interval_tab_fetch_ka(&tab, x_h2o, &k_gas));
        break;
      case HTRDR_Ks:
        HTGOP(layer_sw_spectral_interval_tab_fetch_ks(&tab, x_h2o, &k_gas));
        break;
      case HTRDR_Kext:
        HTGOP(layer_sw_spectral_interval_tab_fetch_kext(&tab, x_h2o, &k_gas));
        break;
      default: FATAL("Unreachable code.\n"); break;
    }
  }

  k = k_particle + k_gas;
  return k;
}

struct svx_tree*
htrdr_sky_get_svx_tree
  (struct htrdr_sky* sky,
   const size_t ispectral_band,
   const size_t iquadrature_pt)
{
  size_t i;
  ASSERT(sky);
  ASSERT(ispectral_band >= sky->sw_bands_range[0]);
  ASSERT(ispectral_band <= sky->sw_bands_range[1]);
  (void)iquadrature_pt;
  i = ispectral_band - sky->sw_bands_range[0];
  return sky->clouds[i].octree;
}

res_T
htrdr_sky_dump_clouds_vtk
  (const struct htrdr_sky* sky,
   const size_t iband, /* Index of the spectral band */
   const size_t iquad, /* Index of the quadrature point */
   FILE* stream)
{
  struct htgop_spectral_interval specint;
  struct svx_tree_desc desc;
  struct octree_data data;
  const double* leaf_data;
  size_t nvertices;
  size_t ncells;
  size_t i;
  ASSERT(sky && stream);
  ASSERT(iband >= sky->sw_bands_range[0]);
  ASSERT(iband <= sky->sw_bands_range[1]);

  i = iband - sky->sw_bands_range[0];

  octree_data_init(sky, iband, iquad, &data);
  SVX(tree_get_desc(sky->clouds[i].octree, &desc));
  ASSERT(desc.type == SVX_OCTREE);

  /* Register leaf data */
  SVX(tree_for_each_leaf(sky->clouds[i].octree, register_leaf, &data));
  nvertices = darray_double_size_get(&data.vertices) / 3/*#coords per vertex*/;
  ncells = darray_size_t_size_get(&data.cells)/8/*#ids per cell*/;
  ASSERT(ncells == desc.nleaves);

  /* Fetch the spectral interval descriptor */
  HTGOP(get_sw_spectral_interval(sky->htgop, iband, &specint));

  /* Write headers */
  fprintf(stream, "# vtk DataFile Version 2.0\n");
  fprintf(stream, "Clouds optical properties in [%g, %g] nanometers\n",
    wavenumber_to_wavelength(specint.wave_numbers[1]),
    wavenumber_to_wavelength(specint.wave_numbers[0]));
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
   const size_t iband, /* Index of the spectral band */
   const size_t iquad, /* Index of the quadrature point in the spectral band */
   const double pos[3])
{
  struct svx_voxel voxel = SVX_VOXEL_NULL;
  size_t i;
  int comp_mask = components_mask;
  ASSERT(sky && pos);
  ASSERT(comp_mask & HTRDR_ALL_COMPONENTS);
  ASSERT(iband >= sky->sw_bands_range[0]);
  ASSERT(iband <= sky->sw_bands_range[1]);

  i = iband - sky->sw_bands_range[0];

  /* Is the position outside the clouds? */
  if(pos[0] < sky->clouds[i].octree_desc.lower[0]
  || pos[1] < sky->clouds[i].octree_desc.lower[1]
  || pos[2] < sky->clouds[i].octree_desc.lower[2]
  || pos[0] > sky->clouds[i].octree_desc.upper[0]
  || pos[1] > sky->clouds[i].octree_desc.upper[1]
  || pos[2] > sky->clouds[i].octree_desc.upper[2]) {
    comp_mask &= ~HTRDR_PARTICLES; /* No particle */
  }

  SVX(tree_at(sky->clouds[i].octree, pos, NULL, NULL, &voxel));
  return htrdr_sky_fetch_svx_voxel_property
    (sky, prop, op, comp_mask, iband, iquad, &voxel);
}

double
htrdr_sky_fetch_svx_voxel_property
  (const struct htrdr_sky* sky,
   const enum htrdr_sky_property prop,
   const enum htrdr_svx_op op,
   const int components_mask,
   const size_t ispectral_band, /* Index of the spectral band */
   const size_t iquad, /* Index of the quadrature point in the spectral band */
   const struct svx_voxel* voxel)
{
  const float* par_data = NULL;
  const float* gas_data = NULL;
  int comp_mask = components_mask;
  double gas = 0;
  double par = 0;
  ASSERT(sky && voxel);
  ASSERT((unsigned)prop < HTRDR_PROPERTIES_COUNT__);
  ASSERT((unsigned)op < HTRDR_SVX_OPS_COUNT__);
  (void)sky, (void)ispectral_band, (void)iquad;

  par_data = (const float*)voxel->data + 0*NFLOATS_PER_COMPONENT;
  gas_data = (const float*)voxel->data + 1*NFLOATS_PER_COMPONENT;

  if(comp_mask & HTRDR_PARTICLES) {
    par = par_data[prop * HTRDR_SVX_OPS_COUNT__ + op];
  }
  if(comp_mask & HTRDR_GAS) {
    gas = gas_data[prop * HTRDR_SVX_OPS_COUNT__ + op];
  }
  /* Interval arithmetic to ensure that the returned Min/Max includes the
   * Min/Max of the "gas + particles mixture" */
  return par + gas;
}

size_t
htrdr_sky_get_sw_spectral_bands_count(const struct htrdr_sky* sky)
{
  ASSERT(sky);
  return sky->sw_bands_range[1] - sky->sw_bands_range[0] + 1;
}

size_t
htrdr_sky_get_sw_spectral_band_id
  (const struct htrdr_sky* sky, const size_t i)
{
  ASSERT(sky && i < htrdr_sky_get_sw_spectral_bands_count(sky));
  return sky->sw_bands_range[0] + i;
}

res_T
htrdr_sky_get_sw_spectral_band_bounds
  (const struct htrdr_sky* sky,
   const size_t iband,
   double wavelengths[2])
{
  struct htgop_spectral_interval specint;
  res_T res = RES_OK;
  ASSERT(sky && wavelengths);

  res = htgop_get_sw_spectral_interval(sky->htgop, iband, &specint);
  if(res != RES_OK) return res;

  wavelengths[0] = wavenumber_to_wavelength(specint.wave_numbers[1]);
  wavelengths[1] = wavenumber_to_wavelength(specint.wave_numbers[0]);
  ASSERT(wavelengths[0] < wavelengths[1]);
  return RES_OK;
}

void
htrdr_sky_sample_sw_spectral_data_CIE_1931_X
  (const struct htrdr_sky* sky,
   struct ssp_rng* rng,
   size_t* ispectral_band,
   size_t* iquadrature_pt)
{
  sample_sw_spectral_data
    (sky->htgop, rng, htgop_sample_sw_spectral_interval_CIE_1931_X,
     ispectral_band, iquadrature_pt);
}

void
htrdr_sky_sample_sw_spectral_data_CIE_1931_Y
  (const struct htrdr_sky* sky,
   struct ssp_rng* rng,
   size_t* ispectral_band,
   size_t* iquadrature_pt)
{
  sample_sw_spectral_data
    (sky->htgop, rng, htgop_sample_sw_spectral_interval_CIE_1931_Y,
     ispectral_band, iquadrature_pt);
}

void
htrdr_sky_sample_sw_spectral_data_CIE_1931_Z
  (const struct htrdr_sky* sky,
   struct ssp_rng* rng,
   size_t* ispectral_band,
   size_t* iquadrature_pt)
{
  sample_sw_spectral_data
    (sky->htgop, rng, htgop_sample_sw_spectral_interval_CIE_1931_Z,
     ispectral_band, iquadrature_pt);
}

