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

#ifndef HTRDR_SKY_H
#define HTRDR_SKY_H

#include <rsys/rsys.h>
#include <star/svx.h>

/* Properties of the sky */
enum htrdr_sky_property {
  HTRDR_Ks, /* Scattering coefficient */
  HTRDR_Ka, /* Absorption coefficient */
  HTRDR_Kext, /* Extinction coefficient = Ks + Ka */
  HTRDR_PROPERTIES_COUNT__
};

/* FIXME Maybe rename this constant to avoid the confusion with the
 * htrdr_sky_component_flag enumerate */
enum htrdr_sky_component {
  HTRDR_GAS__,
  HTRDR_PARTICLES__,
  HTRDR_COMPONENTS_COUNT__
};

/* Component of the sky for which the properties are queried */
enum htrdr_sky_component_flag {
  HTRDR_GAS = BIT(HTRDR_GAS__),
  HTRDR_PARTICLES = BIT(HTRDR_PARTICLES__),
  HTRDR_ALL_COMPONENTS = HTRDR_GAS | HTRDR_PARTICLES
};

enum htrdr_svx_op {
  HTRDR_SVX_MIN,
  HTRDR_SVX_MAX,
  HTRDR_SVX_OPS_COUNT__
};

/* Forward declarations */
struct htrdr;
struct htrdr_sky;
struct htrdr_sun;
struct ssp_rng;
struct svx_tree;
struct svx_voxel;

extern LOCAL_SYM res_T
htrdr_sky_create
  (struct htrdr* htrdr,
   struct htrdr_sun* sun,
   const char* htcp_filename,
   const char* htgop_filename,
   const char* htmie_filename,
   const double optical_thickness, /* Threshold used during octree building */
   const int repeat_clouds, /* Infinitely repeat the clouds in X and Y */
   struct htrdr_sky** sky);

extern LOCAL_SYM void
htrdr_sky_ref_get
  (struct htrdr_sky* sky);

extern LOCAL_SYM void
htrdr_sky_ref_put
  (struct htrdr_sky* sky);

extern LOCAL_SYM double
htrdr_sky_fetch_particle_phase_function_asymmetry_parameter
  (const struct htrdr_sky* sky,
   const size_t ispectral_band,
   const size_t iquadrature_pt);

extern LOCAL_SYM double
htrdr_sky_fetch_raw_property
  (const struct htrdr_sky* sky,
   const enum htrdr_sky_property prop,
   const int comp_mask, /* Combination of htrdr_sky_component_flag */
   const size_t ispectral_band,
   const size_t iquadrature_pt,
   const double pos[3],
   const double k_min, /* For debug only */
   const double k_max); /* For debug only */

extern LOCAL_SYM double
htrdr_sky_fetch_svx_property
  (const struct htrdr_sky* sky,
   const enum htrdr_sky_property prop,
   const enum htrdr_svx_op op,
   const int comp_mask, /* Combination of htrdr_sky_component_flag */
   const size_t ispectral_band,
   const size_t iquadrature_pt,
   const double pos[3]);

extern LOCAL_SYM double
htrdr_sky_fetch_svx_voxel_property
  (const struct htrdr_sky* sky,
   const enum htrdr_sky_property prop,
   const enum htrdr_svx_op op,
   const int comp_mask, /* Combination of htrdr_sky_component_flag */
   const size_t ispectral_band,
   const size_t iquadrature_pt,
   const struct svx_voxel* voxel);

extern LOCAL_SYM res_T
htrdr_sky_dump_clouds_vtk
  (const struct htrdr_sky* sky,
   const size_t ispectral_band,
   const size_t iquadrature_pt,
   FILE* stream);

extern LOCAL_SYM size_t
htrdr_sky_get_sw_spectral_bands_count
  (const struct htrdr_sky* sky);

extern LOCAL_SYM size_t
htrdr_sky_get_sw_spectral_band_id
  (const struct htrdr_sky* sky,
   const size_t i); /* in [0, htrdr_sky_get_sw_spectral_bands_count[ */

extern LOCAL_SYM size_t
htrdr_sky_get_sw_spectral_band_quadrature_length
  (const struct htrdr_sky* sky,
   const size_t iband);

extern LOCAL_SYM res_T
htrdr_sky_get_sw_spectral_band_bounds
  (const struct htrdr_sky* sky,
   const size_t iband,
   double wavelengths[2]);

extern LOCAL_SYM void
htrdr_sky_sample_sw_spectral_data_CIE_1931_X
  (const struct htrdr_sky* sky,
   struct ssp_rng* rng,
   size_t* ispectral_band,
   size_t* iquadrature_pt);

extern LOCAL_SYM void
htrdr_sky_sample_sw_spectral_data_CIE_1931_Y
  (const struct htrdr_sky* sky,
   struct ssp_rng* rng,
   size_t* ispectral_band,
   size_t* iquadrature_pt);

extern LOCAL_SYM void
htrdr_sky_sample_sw_spectral_data_CIE_1931_Z
  (const struct htrdr_sky* sky,
   struct ssp_rng* rng,
   size_t* ispectral_band,
   size_t* iquadrature_pt);

extern LOCAL_SYM res_T
htrdr_sky_trace_ray
  (struct htrdr_sky* sky,
   const double ray_origin[3],
   const double ray_direction[3], /* Must be normalized */
   const double ray_range[2],
   const svx_hit_challenge_T challenge, /* NULL <=> Traversed up to the leaves */
   const svx_hit_filter_T filter, /* NULL <=> Stop RT at the 1st hit voxel */
   void* context, /* Data sent to the filter functor */
   const size_t ispectral_band,
   const size_t iquadrature_pt,
   struct svx_hit* hit);

#endif /* HTRDR_SKY_H */

