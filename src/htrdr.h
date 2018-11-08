/* Copyright (C) 2018 CNRS, Université Paul Sabatier, |Meso|Star>
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
struct mutext;
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
  size_t width; /* Image width */
  size_t height; /* Image height */

  FILE* output;
  struct str output_name;

  unsigned nthreads; /* #threads of the process */
  int dump_vtk; /* Dump octree VTK */
  int cache_grids; /* Use/Precompute grid caches */
  int verbose; /* Verbosity level */

  int mpi_rank; /* Rank of the process in the MPI group */
  int mpi_nprocs; /* Overall #processes in the MPI group */
  char* mpi_err_str; /* Temp buffer used to store MPI error string */
  int8_t* mpi_working_procs; /* Define the rank of active processes */
  size_t mpi_nworking_procs;

  /* Process progress percentage */
  int32_t* mpi_progress_octree;
  int32_t* mpi_progress_render;

  struct mutex* mpi_mutex; /* Protect MPI calls from concurrent threads */

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

extern LOCAL_SYM const char*
htrdr_mpi_error_string
  (struct htrdr* htrdr,
   const int mpi_err);

extern LOCAL_SYM void
htrdr_fprintf
  (struct htrdr* htrdr,
   FILE* stream,
   const char* msg,
   ...)
#ifdef COMPILER_GCC
  __attribute((format(printf, 3, 4)))
#endif
  ;

extern LOCAL_SYM void
htrdr_fflush
  (struct htrdr* htrdr,
   FILE* stream);

#endif /* HTRDR_H */

