/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2023 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2023 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2023 Université Paul Sabatier
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

#ifndef HTRDR_INTERFACE_H
#define HTRDR_INTERFACE_H

#include "core/htrdr_materials.h"
#include <star/s3d.h>
#include <rsys/double3.h>

/* Forward declaration of external data type */
struct mrumtl;
struct s3d_hit;
struct ssf_bsdf;
struct ssp_rng;

struct htrdr_interface {
  struct htrdr_mtl mtl_front;
  struct htrdr_mtl mtl_back;
  struct htrdr_mtl mtl_thin; /* != NULL <=> thin material */
};
static const struct htrdr_interface HTRDR_INTERFACE_NULL;

static INLINE const struct htrdr_mtl*
htrdr_interface_fetch_hit_mtl
  (const struct htrdr_interface* interf,
   const double dir[3], /* Incoming ray */
   const struct s3d_hit* hit)
{
  const struct htrdr_mtl* mtl = NULL;
  enum { FRONT, BACK };
  ASSERT(interf && dir && d3_is_normalized(dir) && hit && !S3D_HIT_NONE(hit));
  ASSERT(interf->mtl_front.mrumtl
    || interf->mtl_back.mrumtl
    || interf->mtl_thin.mrumtl);

  if(interf->mtl_thin.mrumtl) {
    mtl = &interf->mtl_thin;
  } else {
    double N[3];
    int hit_side;
    d3_normalize(N, d3_set_f3(N, hit->normal));
    hit_side = d3_dot(N, dir) < 0 ? FRONT : BACK;

    /* Retrieve the brdf of the material on the *other side* of the hit side */
    switch(hit_side) {
      case BACK: mtl = &interf->mtl_front; break;
      case FRONT: mtl = &interf->mtl_back; break;
      default: FATAL("Unreachable code.\n");  break;
    }

    /* Due to numerical issue the hit side might be wrong and thus the fetched
     * material might be undefined (e.g. semi-transparent materials). Handle this
     * issue by fetching the other material. */
    if(!mtl->mrumtl) {
      switch(hit_side) {
        case BACK: mtl = &interf->mtl_back; break;
        case FRONT: mtl = &interf->mtl_front; break;
        default: FATAL("Unreachable code.\n");  break;
      }
    }
    ASSERT(mtl->mrumtl);
  }

  return mtl;
}

#endif /* HTRDR_INTERFACE_H */

