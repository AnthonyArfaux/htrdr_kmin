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

#ifndef HTRDR_BUFFER_H
#define HTRDR_BUFFER_H

#include "core/htrdr.h"

#include <rsys/math.h>
#include <rsys/rsys.h>

/*
 * Row major ordered 2D buffer
 */

struct htrdr_pixel_format {
  size_t size; /* In bytes */
  size_t alignment; /* Power of two, in Bytes */
};
#define HTRDR_PIXEL_FORMAT_NULL__ {0, 0}
static const struct htrdr_pixel_format HTRDR_PIXEL_FORMAT_NULL =
  HTRDR_PIXEL_FORMAT_NULL__;

struct htrdr_buffer_layout {
  size_t width; /* #elements in X */
  size_t height; /* #elements in Y */
  size_t pitch; /* #Bytes between 2 consecutive line */
  size_t elmt_size; /* Size of an element in the buffer */
  size_t alignment; /* Alignement of the memory */
};
#define HTRDR_BUFFER_LAYOUT_NULL__ {0,0,0,0,0}
static const struct htrdr_buffer_layout HTRDR_BUFFER_LAYOUT_NULL =
  HTRDR_BUFFER_LAYOUT_NULL__;

static INLINE int
htrdr_buffer_layout_eq
  (const struct htrdr_buffer_layout* a,
   const struct htrdr_buffer_layout* b)
{
  ASSERT(a && b);
  return a->width == b->width
      && a->height == b->height
      && a->pitch == b->pitch
      && a->elmt_size == b->elmt_size
      && a->alignment == b->alignment;
}

static INLINE int
htrdr_buffer_layout_check(const struct htrdr_buffer_layout* layout)
{
  return layout
      && layout->width
      && layout->height
      && layout->elmt_size
      && layout->width*layout->elmt_size <= layout->pitch
      && IS_POW2(layout->alignment);
}

/* Forward declarations */
struct htrdr;
struct htrdr_buffer;

BEGIN_DECLS

HTRDR_API res_T
htrdr_buffer_create
  (struct htrdr* htrdr,
   const struct htrdr_buffer_layout* layout,
   struct htrdr_buffer** buf);

HTRDR_API void
htrdr_buffer_ref_get
  (struct htrdr_buffer* buf);

HTRDR_API void
htrdr_buffer_ref_put
  (struct htrdr_buffer* buf);

HTRDR_API void
htrdr_buffer_get_layout
  (const struct htrdr_buffer* buf,
   struct htrdr_buffer_layout* layout);

HTRDR_API void*
htrdr_buffer_get_data
  (struct htrdr_buffer* buf);

HTRDR_API void*
htrdr_buffer_at
  (struct htrdr_buffer* buf,
   const size_t x,
   const size_t y);

END_DECLS

#endif /* HTRDR_BUFFER_H */
