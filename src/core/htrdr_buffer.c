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

#include "htrdr.h"
#include "htrdr_buffer.h"
#include "htrdr_log.h"

#include <rsys/math.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

struct htrdr_buffer {
  char* mem;

  size_t width;
  size_t height;
  size_t pitch;
  size_t elmtsz;
  size_t align;

  struct htrdr* htrdr;
  ref_T ref;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
buffer_release(ref_T* ref)
{
  struct htrdr_buffer* buf = NULL;
  struct htrdr* htrdr = NULL;
  ASSERT(ref);
  buf = CONTAINER_OF(ref, struct htrdr_buffer, ref);
  htrdr = buf->htrdr;
  if(buf->mem) MEM_RM(htrdr_get_allocator(htrdr), buf->mem);
  MEM_RM(htrdr_get_allocator(htrdr), buf);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_buffer_create
  (struct htrdr* htrdr,
   const size_t width,
   const size_t height,
   const size_t pitch,
   const size_t elmtsz,
   const size_t align,
   struct htrdr_buffer** out_buf)
{
  struct htrdr_buffer* buf = NULL;
  size_t memsz = 0;
  res_T res = RES_OK;
  ASSERT(htrdr && out_buf);

  if(!width || !height) {
    htrdr_log_err(htrdr, "invalid buffer definition %lux%lu.\n",
      (unsigned long)width, (unsigned long)height);
    res = RES_BAD_ARG;
    goto error;
  }
  if(pitch < width*elmtsz) {
    htrdr_log_err(htrdr,
      "invalid buffer pitch `%lu' wrt the buffer width `%lu'. "
      "The buffer pitch cannot be less than the buffer width.\n",
      (unsigned long)pitch, (unsigned long)width);
    res = RES_BAD_ARG;
    goto error;
  }
  if(!elmtsz) {
    htrdr_log_err(htrdr,
      "the size of the buffer's elements cannot be null.\n");
    res = RES_BAD_ARG;
    goto error;
  }
  if(!IS_POW2(align)) {
    htrdr_log_err(htrdr,
      "invalid buffer alignment `%lu'. It must be a power of 2.\n",
      (unsigned long)align);
    res = RES_BAD_ARG;
    goto error;
  }

  buf = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*buf));
  if(!buf) {
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&buf->ref);
  buf->width = width;
  buf->height = height;
  buf->pitch = pitch;
  buf->elmtsz = elmtsz;
  buf->align = align;
  htrdr_ref_get(htrdr);
  buf->htrdr = htrdr;

  memsz = buf->pitch * buf->height;
  buf->mem = MEM_ALLOC_ALIGNED(htrdr_get_allocator(htrdr), memsz, align);
  if(!buf->mem) {
    res = RES_MEM_ERR;
    goto error;
  }

exit:
  *out_buf = buf;
  return res;
error:
  if(buf) {
    htrdr_buffer_ref_put(buf);
    buf = NULL;
  }
  goto exit;
}

void
htrdr_buffer_ref_get(struct htrdr_buffer* buf)
{
  ASSERT(buf);
  ref_get(&buf->ref);
}

void
htrdr_buffer_ref_put(struct htrdr_buffer* buf)
{
  ASSERT(buf);
  ref_put(&buf->ref, buffer_release);
}

void
htrdr_buffer_get_layout
  (const struct htrdr_buffer* buf,
   struct htrdr_buffer_layout* layout)
{
  ASSERT(buf && layout);
  layout->width = buf->width;
  layout->height = buf->height;
  layout->pitch = buf->pitch;
  layout->elmt_size = buf->elmtsz;
  layout->alignment = buf->align;
}

void*
htrdr_buffer_get_data(struct htrdr_buffer* buf)
{
  ASSERT(buf);
  return buf->mem;
}

void*
htrdr_buffer_at(struct htrdr_buffer* buf, const size_t x, const size_t y)
{
  ASSERT(buf && x < buf->width && y < buf->height);
  return buf->mem + y*buf->pitch + x*buf->elmtsz;
}

