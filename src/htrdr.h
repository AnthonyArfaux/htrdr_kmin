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

#ifndef HTRDR_H
#define HTRDR_H

#include <rsys/logger.h>
#include <rsys/ref_count.h>

/* Forward declaration */
struct htrdr_args;
struct htrdr_buffer;
struct htrdr_rectangle;
struct mem_allocator;
struct svx_device;
struct svx_tree;

struct htrdr {
  struct svx_device* svx;
  struct svx_tree* clouds;

  struct htrdr_buffer* buf;
  struct htrdr_rectangle* rect;
  double main_dir[3]; /* Main direction */
  size_t spp; /* #samples per pixel */

  FILE* output;

  unsigned nthreads;
  int dump_vtk;
  int verbose;

  struct logger logger;
  struct mem_allocator* allocator;
};

extern LOCAL_SYM res_T
htrdr_init
  (struct mem_allocator* allocator,
   const struct htrdr_args* args,
   struct htrdr* htrdr);

extern LOCAL_SYM void
htrdr_release
  (struct htrdr* htrdr);

extern LOCAL_SYM res_T
htrdr_run
  (struct htrdr* htrdr);

extern LOCAL_SYM void
htrdr_log_err
  (struct htrdr* htrdr,
   const char* msg,
   ...)
#ifdef COMPILER_GCC
  __attribute((format(printf, 2, 3)))
#endif
  ;

extern LOCAL_SYM void
htrdr_log_warn
  (struct htrdr* htrdr,
   const char* msg,
   ...)
#ifdef COMPILER_GCC
  __attribute((format(printf, 2, 3)))
#endif
  ;

#endif /* HTRDR_H */

