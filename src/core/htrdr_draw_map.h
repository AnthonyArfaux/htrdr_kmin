/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2023 Observatoire de Paris
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

#ifndef HTRDR_DRAW_MAP_H
#define HTRDR_DRAW_MAP_H

#include "core/htrdr.h"
#include "core/htrdr_buffer.h"

#include <rsys/rsys.h>

struct htrdr_draw_pixel_args {
  size_t pixel_coord[2]; /* Image plane pixel coordinates */
  double pixel_normalized_size[2]; /* Pixel size in the normalized img plane */
  struct ssp_rng* rng; /* Random Number Generator */
  size_t spp; /* #samples per pixel */
  size_t ithread; /* Id of the thread drawing the pixel */
  void* context; /* User defined data */
};

#define HTRDR_DRAW_PIXEL_ARGS_NULL__ {                                         \
  {0, 0}, /* Image plane pixel coordinates */                                  \
  {0, 0}, /* Pixel size in the normalized img plane */                         \
  NULL, /* RNG */                                                              \
  0, /* SPP */                                                                 \
  0, /* Thread id */                                                           \
  NULL /* User data */                                                         \
}
static const struct htrdr_draw_pixel_args HTRDR_DRAW_PIXEL_ARGS_NULL =
  HTRDR_DRAW_PIXEL_ARGS_NULL__;

typedef void
(*htrdr_draw_pixel_T)
  (struct htrdr* htrdr,
   const struct htrdr_draw_pixel_args* args,
   void* pixel); /* Output data */

struct htrdr_draw_map_args {
  htrdr_draw_pixel_T draw_pixel;
  struct htrdr_buffer_layout buffer_layout;
  size_t spp; /* Samples per pixel */
  void* context; /* User defined data */
};

#define HTRDR_DRAW_MAP_ARGS_NULL__ {                                           \
  NULL, /* Draw pixel functor */                                               \
  HTRDR_BUFFER_LAYOUT_NULL__, /* Layout of the destination buffer */           \
  0, /* #Samples per pixel */                                                  \
  NULL /* User defined data */                                                 \
}
static const struct htrdr_draw_map_args HTRDR_DRAW_MAP_ARGS_NULL =
  HTRDR_DRAW_MAP_ARGS_NULL__;

static INLINE int
htrdr_draw_pixel_args_check(const struct htrdr_draw_pixel_args* args)
{
  return args
      && args->pixel_normalized_size[0] > 0
      && args->pixel_normalized_size[1] > 0
      && args->rng
      && args->spp > 0;
}

/*******************************************************************************
 * Exported symbols
 ******************************************************************************/
BEGIN_DECLS

HTRDR_API res_T
htrdr_draw_map
  (struct htrdr* htrdr,
   const struct htrdr_draw_map_args* args,
   struct htrdr_buffer* buf); /* May be NULL for non master processes */

END_DECLS

#endif /* HTRDR_DRAW_MAP_H */

