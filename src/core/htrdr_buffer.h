/* Copyright (C) 2018, 2019, 2020, 2021 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019, 2021 CNRS
 * Copyright (C) 2018, 2019 Université Paul Sabatier
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

#include "htrdr.h"
#include <rsys/rsys.h>

/*
 * Row major ordered 2D buffer
 */

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

/* Forward declarations */
struct htrdr;
struct htrdr_buffer;

BEGIN_DECLS

HTRDR_API res_T
htrdr_buffer_create
  (struct htrdr* htrdr,
   const size_t width,
   const size_t height,
   const size_t pitch, /* #Bytes between 2 consecutive line */
   const size_t elmt_size, /* Size of an element in the buffer */
   const size_t alignment, /* Alignement of the buffer */
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
