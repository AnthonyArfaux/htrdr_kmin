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

#ifndef HTRDR_H
#define HTRDR_H

#include <rsys/rsys.h>
#include <stdio.h>

/* Library symbol management */
#if defined(HTRDR_CORE_SHARED_BUILD) /* Build shared library */
  #define HTRDR_CORE_API extern EXPORT_SYM
#elif defined(HTRDR_CORE_STATIC) /* Use/build static library */
  #define HTRDR_CORE_API extern LOCAL_SYM
#else /* Use shared library */
  #define HTRDR_CORE_API extern IMPORT_SYM
#endif

#if defined(HTRDR_SHARED_BUILD) /* Build shared library */
  #define HTRDR_API extern EXPORT_SYM
#elif defined(HTRDR_STATIC) /* Use/build static library */
  #define HTRDR_API extern LOCAL_SYM
#else /* Use shared library */
  #define HTRDR_API extern IMPORT_SYM
#endif

/* Helper macro that asserts if the invocation of the htrdr function `Func'
 * returns an error. One should use this macro on htrdr function calls for
 * which no explicit error checking is performed */
#ifndef NDEBUG
  #define HTRDR(Func) ASSERT(htrdr_ ## Func == RES_OK)
#else
  #define HTRDR(Func) htrdr_ ## Func
#endif

/* Forward declarations */
struct htrdr_buffer;
struct mem_allocator;
struct mutex;

struct htrdr_args {
  unsigned nthreads; /* #threads of the process */
  int verbose; /* Verbosity level */
};
#define HTRDR_ARGS_DEFAULT__ { (unsigned)~0, 1 }
static const struct htrdr_args HTRDR_ARGS_DEFAULT = HTRDR_ARGS_DEFAULT__;

/* Forward declaration */
struct htrdr;

static INLINE void
htrdr_fprint_copyright(const char* cmd, FILE* stream)
{
  (void)cmd;
  fprintf(stream,
"Copyright (C) 2018, 2019, 2020, 2021 |Meso|Star> <contact@meso-star.com>.\n"
"Copyright (C) 2018, 2019, 2021 CNRS.\n"
"Copyright (C) 2018, 2019 Université Paul Sabatier.\n");
}

static INLINE void
htrdr_fprint_license(const char* cmd, FILE* stream)
{
  ASSERT(cmd);
  fprintf(stream,
"%s is free software released under the GNU GPL license, version\n"
"3 or later. You are free to change or redistribute it under certain\n"
"conditions <http://gnu.org/licenses/gpl.html>.\n", cmd);
}

BEGIN_DECLS

/* Initialize the MPI execution environment. Must be called priorly to any MPI
 * invocation, e.g. at the beginning of the main function */
HTRDR_CORE_API res_T
htrdr_mpi_init
  (int argc,
   char** argv);

/* Terminate the MPI execution environment */
HTRDR_CORE_API void
htrdr_mpi_finalize
  (void);

/*******************************************************************************
 * HTRDR api
 ******************************************************************************/
HTRDR_CORE_API res_T
htrdr_create
  (struct mem_allocator* allocator,
   const struct htrdr_args* args,
   struct htrdr** htrdr);

HTRDR_CORE_API void
htrdr_ref_get
  (struct htrdr* htrdr);

HTRDR_CORE_API void
htrdr_ref_put
  (struct htrdr* htrdr);

/* Return the number of threads used by the process */
HTRDR_CORE_API size_t
htrdr_get_threads_count
  (const struct htrdr* htrdr);

/* Return the number of running processes for the current htrdr instance */
HTRDR_CORE_API size_t
htrdr_get_procs_count
  (const struct htrdr* htrdr);

HTRDR_CORE_API int
htrdr_get_mpi_rank
  (const struct htrdr* htrdr);

HTRDR_CORE_API struct mem_allocator*
htrdr_get_allocator
  (struct htrdr* htrdr);

HTRDR_CORE_API struct mem_allocator*
htrdr_get_thread_allocator
  (struct htrdr* htrdr,
   const size_t ithread);

HTRDR_CORE_API struct logger*
htrdr_get_logger
  (struct htrdr* htrdr);

HTRDR_CORE_API int
htrdr_get_verbosity_level
  (const struct htrdr* htrdr);

HTRDR_CORE_API struct s3d_device*
htrdr_get_s3d
  (struct htrdr* htrdr);

HTRDR_CORE_API res_T
htrdr_open_output_stream
  (struct htrdr* htrdr,
   const char* filename,
   const int read,
   int force_overwrite,
   FILE** out_fp);

/* TODO do not expose publicly this function(?) */
HTRDR_CORE_API const char*
htrdr_mpi_error_string
  (struct htrdr* htrdr,
   const int mpi_err);

/* TODO replace them by regular log message */
HTRDR_CORE_API void
htrdr_fprintf
  (struct htrdr* htrdr,
   FILE* stream,
   const char* msg,
   ...)
#ifdef COMPILER_GCC
  __attribute((format(printf, 3, 4)))
#endif
  ;

/* TODO remove this */
HTRDR_CORE_API void
htrdr_fflush
  (struct htrdr* htrdr,
   FILE* stream);


END_DECLS

#endif /* HTRDR_H */

