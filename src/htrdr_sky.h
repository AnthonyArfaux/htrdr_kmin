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

/* Raw sky properties */
enum htrdr_sky_property {
  HTRDR_SKY_Ks, /* Scattering coefficient */
  HTRDR_SKY_Ka /* Absorption coefficient */
};

/* Property of the sky computed by region and managed by Star-VoXel */
enum htrdr_sky_svx_property {
  HTRDR_SKY_SVX_Kext_MIN,
  HTRDR_SKY_SVX_Kext_MAX,
  HTRDR_SKY_SVX_PROPS_COUNT__
};

/* Component of the sky for which the properties are queried */
enum htrdr_sky_component_flag {
  HTRDR_SKY_GAZ = BIT(0),
  HTRDR_SKY_PARTICLE = BIT(1)
};

/* Forward declaration */
struct htrdr;
struct htrdr_sky;
struct svx_tree;
struct svx_voxel;

extern LOCAL_SYM res_T
htrdr_sky_create
  (struct htrdr* htrdr,
   const char* htcp_filename,
   struct htrdr_sky** sky);

extern LOCAL_SYM void
htrdr_sky_ref_get
  (struct htrdr_sky* sky);

extern LOCAL_SYM void
htrdr_sky_ref_put
  (struct htrdr_sky* sky);

extern LOCAL_SYM double
htrdr_sky_fetch_raw_property
  (const struct htrdr_sky* sky,
   const enum htrdr_sky_property prop,
   const int components_mask, /* Combination of htrdr_sky_component_flag */
   const double wavelength,
   const double pos[3]);

extern LOCAL_SYM struct svx_tree*
htrdr_sky_get_svx_tree
  (struct htrdr_sky* sky);

extern LOCAL_SYM double
htrdr_sky_fetch_svx_property
  (const struct htrdr_sky* sky,
   const enum htrdr_sky_svx_property prop,
   const int components_mask, /* Combination of htrdr_sky_component_flag */
   const double wavelength,
   const double pos[3]);

extern LOCAL_SYM double
htrdr_sky_fetch_svx_voxel_property
  (const struct htrdr_sky* sky,
   const enum htrdr_sky_svx_property prop,
   const int components_mask, /* Combination of htrdr_sky_component_flag */
   const double wavelength,
   const struct svx_voxel* voxel);

extern LOCAL_SYM res_T
htrdr_sky_dump_clouds_vtk
  (const struct htrdr_sky* sky,
   FILE* stream);

#endif /* HTRDR_SKY_H */

