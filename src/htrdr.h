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
#include <rsys/str.h>

/* Helper macro that asserts if the invocation of the htrdr function `Func'
 * returns an error. One should use this macro on htcp function calls for
 * which no explicit error checking is performed */
#ifndef NDEBUG
  #define HTRDR(Func) ASSERT(htrdr_ ## Func == RES_OK)
#else
  #define HTRDR(Func) htrdr_ ## Func
#endif

/* Forward declarations */
struct htrdr_args;
struct htrdr_buffer;
struct htrdr_sky;
struct htrdr_rectangle;
struct mem_allocator;
struct s3d_device;
struct s3d_scene;
struct ssf_bsdf;
struct ssf_phase;
struct svx_device;

struct htrdr {
  struct svx_device* svx;
  struct s3d_device* s3d;

  struct htrdr_ground* ground;
  struct htrdr_sky* sky;
  struct htrdr_sun* sun;

  struct htrdr_camera* cam;
  struct htrdr_buffer* buf;
  size_t spp; /* #samples per pixel */

  FILE* output;
  struct str output_name;

  unsigned nthreads;
  int dump_vtk;
  int verbose;

  struct logger logger;
  struct mem_allocator* allocator;
  struct mem_allocator* lifo_allocators; /* Per thread lifo allocator */
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
htrdr_log
  (struct htrdr* htrdr,
   const char* msg,
   ...)
#ifdef COMPILER_GCC
  __attribute((format(printf, 2, 3)))
#endif
  ;

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

