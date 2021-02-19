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

#define _POSIX_C_SOURCE 200809L /* stat.st_time support */

#include "core/htrdr.h"
#include "core/htrdr_c.h"
#include "core/htrdr_args.h"
#include "core/htrdr_log.h"
#include "core/htrdr_version.h"

#include <rsys/cstr.h>
#include <rsys/mem_allocator.h>
#include <rsys/str.h>

#include "high_tune/htsky.h"

#include <star/s3d.h>
#include <star/ssf.h>

#include <errno.h>
#include <fcntl.h> /* open */
#include <libgen.h> /* basename */
#include <stdarg.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h> /* timespec */
#include <sys/stat.h> /* S_IRUSR & S_IWUSR */

#include <omp.h>
#include <mpi.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static const char*
mpi_thread_support_string(const int val)
{
  switch(val) {
    case MPI_THREAD_SINGLE: return "MPI_THREAD_SINGLE";
    case MPI_THREAD_FUNNELED: return "MPI_THREAD_FUNNELED";
    case MPI_THREAD_SERIALIZED: return "MPI_THREAD_SERIALIZED";
    case MPI_THREAD_MULTIPLE: return "MPI_THREAD_MULTIPLE";
    default: FATAL("Unreachable code.\n"); break;
  }
}

static void
release_mpi(struct htrdr* htrdr)
{
  ASSERT(htrdr);
  if(htrdr->mpi_working_procs) {
    MEM_RM(htrdr->allocator, htrdr->mpi_working_procs);
    htrdr->mpi_working_procs = NULL;
  }
  if(htrdr->mpi_progress_octree) {
    MEM_RM(htrdr->allocator, htrdr->mpi_progress_octree);
    htrdr->mpi_progress_octree = NULL;
  }
  if(htrdr->mpi_progress_render) {
    MEM_RM(htrdr->allocator, htrdr->mpi_progress_render);
    htrdr->mpi_progress_render = NULL;
  }
  if(htrdr->mpi_err_str) {
    MEM_RM(htrdr->allocator, htrdr->mpi_err_str);
    htrdr->mpi_err_str = NULL;
  }
  if(htrdr->mpi_mutex) {
    mutex_destroy(htrdr->mpi_mutex);
    htrdr->mpi_mutex = NULL;
  }
}

static res_T
mpi_print_proc_info(struct htrdr* htrdr)
{
  char proc_name[MPI_MAX_PROCESSOR_NAME];
  int proc_name_len;
  char* proc_names = NULL;
  uint32_t* proc_nthreads = NULL;
  uint32_t nthreads = 0;
  int iproc;
  res_T res = RES_OK;
  ASSERT(htrdr);

  if(htrdr->mpi_rank == 0) {
    proc_names = MEM_CALLOC(htrdr->allocator, (size_t)htrdr->mpi_nprocs,
      MPI_MAX_PROCESSOR_NAME*sizeof(*proc_names));
    if(!proc_names) {
      res = RES_MEM_ERR;
      htrdr_log_err(htrdr,
        "could not allocate the temporary memory for MPI process names -- "
        "%s.\n", res_to_cstr(res));
      goto error;
    }

    proc_nthreads = MEM_CALLOC(htrdr->allocator, (size_t)htrdr->mpi_nprocs,
      sizeof(*proc_nthreads));
    if(!proc_nthreads) {
      res = RES_MEM_ERR;
      htrdr_log_err(htrdr,
        "could not allocate the temporary memory for the #threads of the MPI "
        "processes -- %s.\n", res_to_cstr(res));
      goto error;
    }
  }

  /* Gather process name */
  MPI(Get_processor_name(proc_name, &proc_name_len));
  MPI(Gather(proc_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, proc_names,
    MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, MPI_COMM_WORLD));

  /* Gather process #threads */
  nthreads = (uint32_t)htrdr->nthreads;
  MPI(Gather(&nthreads, 1, MPI_UINT32_T, proc_nthreads, 1, MPI_UINT32_T, 0,
    MPI_COMM_WORLD));

  if(htrdr->mpi_rank == 0) {
    FOR_EACH(iproc, 0, htrdr->mpi_nprocs) {
      htrdr_log(htrdr, "Process %d -- %s; #threads: %u\n",
        iproc, proc_names + iproc*MPI_MAX_PROCESSOR_NAME, proc_nthreads[iproc]);
    }
  }

exit:
  if(proc_names) MEM_RM(htrdr->allocator, proc_names);
  if(proc_nthreads) MEM_RM(htrdr->allocator, proc_nthreads);
  return res;
error:
  goto exit;
}

static res_T
init_mpi(struct htrdr* htrdr)
{
  size_t n;
  int err;
  res_T res = RES_OK;
  ASSERT(htrdr);

  htrdr->mpi_err_str = MEM_CALLOC
    (htrdr->allocator, htrdr->nthreads, MPI_MAX_ERROR_STRING);
  if(!htrdr->mpi_err_str) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "could not allocate the MPI error strings -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  err = MPI_Comm_rank(MPI_COMM_WORLD, &htrdr->mpi_rank);
  if(err != MPI_SUCCESS) {
    htrdr_log_err(htrdr,
      "could not determine the MPI rank of the calling process -- %s.\n",
      htrdr_mpi_error_string(htrdr, err));
    res = RES_UNKNOWN_ERR;
    goto error;
  }

  err = MPI_Comm_size(MPI_COMM_WORLD, &htrdr->mpi_nprocs);
  if(err != MPI_SUCCESS) {
    htrdr_log_err(htrdr,
      "could retrieve the size of the MPI group -- %s.\n",
      htrdr_mpi_error_string(htrdr, err));
    res = RES_UNKNOWN_ERR;
    goto error;
  }

  htrdr->mpi_working_procs = MEM_CALLOC(htrdr->allocator,
    (size_t)htrdr->mpi_nprocs, sizeof(*htrdr->mpi_working_procs));
  if(!htrdr->mpi_working_procs) {
    htrdr_log_err(htrdr,
      "could not allocate the list of working processes.\n");
    res = RES_MEM_ERR;
    goto error;
  }

  /* Initialy, all the processes are working */
  htrdr->mpi_nworking_procs = (size_t)htrdr->mpi_nprocs;
  memset(htrdr->mpi_working_procs, 0xFF,
    htrdr->mpi_nworking_procs*sizeof(*htrdr->mpi_working_procs));

  /* Allocate #processes progress statuses on the master process and only 1
   * progress status on the other ones: the master process will gather the
   * status of the other processes to report their progression. */
  n = (size_t)(htrdr->mpi_rank == 0 ? htrdr->mpi_nprocs : 1);

  htrdr->mpi_progress_octree = MEM_CALLOC
    (htrdr->allocator, n, sizeof(*htrdr->mpi_progress_octree));
  if(!htrdr->mpi_progress_octree) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "could not allocate the progress state of the octree building -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  htrdr->mpi_progress_render = MEM_CALLOC
    (htrdr->allocator, n, sizeof(*htrdr->mpi_progress_render));
  if(!htrdr->mpi_progress_render) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "could not allocate the progress state of the scene rendering -- %s.\n",
      res_to_cstr(res));
    goto error;
  }

  htrdr->mpi_mutex = mutex_create();
  if(!htrdr->mpi_mutex) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "could not create the mutex to protect MPI calls from concurrent "
      "threads -- %s.\n", res_to_cstr(res));
    goto error;
  }

  mpi_print_proc_info(htrdr);

exit:
  return res;
error:
  release_mpi(htrdr);
  goto exit;
}

static void
release_htrdr(ref_T* ref)
{
  struct htrdr* htrdr = CONTAINER_OF(ref, struct htrdr, ref);
  ASSERT(ref);

  release_mpi(htrdr);
  if(htrdr->lifo_allocators) {
    size_t i;
    FOR_EACH(i, 0, htrdr->nthreads) {
      mem_shutdown_lifo_allocator(&htrdr->lifo_allocators[i]);
    }
    MEM_RM(htrdr->allocator, htrdr->lifo_allocators);
  }
  logger_release(&htrdr->logger);

  MEM_RM(htrdr->allocator, htrdr);
}

/*******************************************************************************
 * Exported functions
 ******************************************************************************/
res_T
htrdr_mpi_init(int argc, char** argv)
{
  char str[MPI_MAX_ERROR_STRING];
  int thread_support;
  int len;
  int err = 0;
  res_T res = RES_OK;

  err = MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &thread_support);
  if(err != MPI_SUCCESS) {
    MPI_Error_string(err, str, &len);
    fprintf(stderr, "Error initializing MPI -- %s.\n", str);
    res = RES_UNKNOWN_ERR;
    goto error;
  }

  if(thread_support != MPI_THREAD_SERIALIZED) {
    fprintf(stderr, "The provided MPI implementation does not support "
      "serialized API calls from multiple threads. Provided thread support: "
      "%s.\n", mpi_thread_support_string(thread_support));
    res = RES_BAD_OP;
    goto error;
  }

exit:
  return res;
error:
  goto exit;
}

void
htrdr_mpi_finalize(void)
{
  MPI_Finalize();
}

res_T
htrdr_create
  (struct mem_allocator* mem_allocator,
   const struct htrdr_args* args,
   struct htrdr** out_htrdr)
{
  struct mem_allocator* allocator = NULL;
  struct htrdr* htrdr = NULL;
  size_t ithread;
  int nthreads_max;
  res_T res = RES_OK;
  ASSERT(args && htrdr);

  allocator = mem_allocator ? mem_allocator : &mem_default_allocator;
  htrdr = MEM_CALLOC(allocator, 0, sizeof(*htrdr));
  if(htrdr) {
    fprintf(stderr, HTRDR_LOG_ERROR_PREFIX
      "Could not allocate the htrdr data structure.\n");
    res = RES_MEM_ERR;
    goto error;
  }
  htrdr->allocator = allocator;
  ref_init(&htrdr->ref);
  nthreads_max = MMAX(omp_get_max_threads(), omp_get_num_procs());
  htrdr->verbose = args->verbose;
  htrdr->nthreads = MMIN(args->nthreads, (unsigned)nthreads_max);

  setup_logger(htrdr);

  res = init_mpi(htrdr);
  if(res != RES_OK) goto error;

  htrdr->lifo_allocators = MEM_CALLOC
    (htrdr->allocator, htrdr->nthreads, sizeof(*htrdr->lifo_allocators));
  if(!htrdr->lifo_allocators) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "%s: could not allocate the list of per thread LIFO allocator -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }

  FOR_EACH(ithread, 0, htrdr->nthreads) {
    res = mem_init_lifo_allocator
      (&htrdr->lifo_allocators[ithread], htrdr->allocator, 16384);
    if(res != RES_OK) {
      htrdr_log_err(htrdr,
        "%s: could not initialise the LIFO allocator of the thread %lu -- %s.\n",
        FUNC_NAME, (unsigned long)ithread, res_to_cstr(res));
      goto error;
    }
  }

exit:
  *out_htrdr = htrdr;
  return res;
error:
  if(htrdr) {
    htrdr_ref_put(htrdr);
    htrdr = NULL;
  }
  goto exit;
}

void
htrdr_ref_get(struct htrdr* htrdr)
{
  ASSERT(htrdr);
  ref_get(&htrdr->ref);
}

void
htrdr_ref_put(struct htrdr* htrdr)
{
  ASSERT(htrdr);
  ref_put(&htrdr->ref, release_htrdr);
}

size_t
htrdr_get_threads_count(const struct htrdr* htrdr)
{
  ASSERT(htrdr);
  return htrdr->nthreads;
}

size_t
htrdr_get_procs_count(const struct htrdr* htrdr)
{
  ASSERT(htrdr);
  return (size_t)htrdr->mpi_nprocs;
}

int
htrdr_get_mpi_rank(const struct htrdr* htrdr)
{
  ASSERT(htrdr);
  return htrdr->mpi_rank;
}

struct mem_allocator*
htrdr_get_allocator(struct htrdr* htrdr)
{
  ASSERT(htrdr);
  return htrdr->allocator;
}

struct mem_allocator*
htrdr_get_thread_allocator(struct htrdr* htrdr, const size_t ithread)
{
  ASSERT(htrdr && ithread < htrdr_get_threads_count(htrdr));
  return htrdr->lifo_allocators + ithread;
}

struct logger*
htrdr_get_logger(struct htrdr* htrdr)
{
  ASSERT(htrdr);
  return &htrdr->logger;
}

int
htrdr_get_verbosity_level(const struct htrdr* htrdr)
{
  ASSERT(htrdr);
  return htrdr->verbose;
}

const char*
htrdr_mpi_error_string(struct htrdr* htrdr, const int mpi_err)
{
  const int ithread = omp_get_thread_num();
  char* str;
  int strlen_err;
  int err;
  ASSERT(htrdr && (size_t)ithread < htrdr->nthreads);
  str = htrdr->mpi_err_str + ithread*MPI_MAX_ERROR_STRING;
  err = MPI_Error_string(mpi_err, str, &strlen_err);
  return err == MPI_SUCCESS ? str : "Invalid MPI error";
}

res_T
htrdr_open_output_stream
  (struct htrdr* htrdr,
   const char* filename,
   const int read,
   int force_overwrite,
   FILE** out_fp)
{
  FILE* fp = NULL;
  int fd = -1;
  const char* mode;
  res_T res = RES_OK;
  ASSERT(htrdr && filename && out_fp);

  mode = read ? "w+" : "w";

  if(force_overwrite) {
    fp = fopen(filename, mode);
    if(!fp) {
      htrdr_log_err(htrdr, "could not open the output file `%s'.\n", filename);
      goto error;
    }
  } else {
    const int access_flags = read ? O_RDWR : O_WRONLY;
    fd = open(filename, O_CREAT|O_EXCL|O_TRUNC|access_flags, S_IRUSR|S_IWUSR);
    if(fd >= 0) {
      fp = fdopen(fd, mode);
      if(fp == NULL) {
        htrdr_log_err(htrdr, "could not open the output file `%s'.\n", filename);
        goto error;
      }
    } else if(errno == EEXIST) {
      htrdr_log_err(htrdr, "the output file `%s' already exists. \n",
        filename);
      goto error;
    } else {
      htrdr_log_err(htrdr,
        "unexpected error while opening the output file `%s'.\n", filename);
      goto error;
    }
  }
exit:
  *out_fp = fp;
  return res;
error:
  res = RES_IO_ERR;
  if(fp) {
    CHK(fclose(fp) == 0);
    fp = NULL;
  } else if(fd >= 0) {
    CHK(close(fd) == 0);
  }
  goto exit;
}


void
htrdr_fprintf(struct htrdr* htrdr, FILE* stream, const char* msg, ...)
{
  ASSERT(htrdr && msg);
  if(htrdr->mpi_rank == 0) {
    va_list vargs_list;
    va_start(vargs_list, msg);
    vfprintf(stream, msg, vargs_list);
    va_end(vargs_list);
  }
}

void
htrdr_fflush(struct htrdr* htrdr, FILE* stream)
{
  ASSERT(htrdr);
  if(htrdr->mpi_rank == 0) {
    fflush(stream);
  }
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
void
send_mpi_progress
  (struct htrdr* htrdr, const enum htrdr_mpi_message msg, int32_t percent)
{
  ASSERT(htrdr);
  ASSERT(msg == HTRDR_MPI_PROGRESS_RENDERING);
  (void)htrdr;
  mutex_lock(htrdr->mpi_mutex);
  MPI(Send(&percent, 1, MPI_INT32_T, 0, msg, MPI_COMM_WORLD));
  mutex_unlock(htrdr->mpi_mutex);
}

void
fetch_mpi_progress(struct htrdr* htrdr, const enum htrdr_mpi_message msg)
{
  struct timespec t;
  int32_t* progress = NULL;
  int iproc;
  ASSERT(htrdr && htrdr->mpi_rank == 0);

  t.tv_sec = 0;
  t.tv_nsec = 10000000; /* 10ms */

  switch(msg) {
    case HTRDR_MPI_PROGRESS_RENDERING:
      progress = htrdr->mpi_progress_render;
     break;
    default: FATAL("Unreachable code.\n"); break;
  }

  FOR_EACH(iproc, 1, htrdr->mpi_nprocs) {
    /* Flush the last sent percentage of the process `iproc' */
    for(;;) {
      MPI_Request req;
      int flag;
      int complete;

      mutex_lock(htrdr->mpi_mutex);
      MPI(Iprobe(iproc, msg, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE));
      mutex_unlock(htrdr->mpi_mutex);

      if(flag == 0) break; /* No more message */

      mutex_lock(htrdr->mpi_mutex);
      MPI(Irecv(&progress[iproc], 1, MPI_INT32_T, iproc, msg, MPI_COMM_WORLD, &req));
      mutex_unlock(htrdr->mpi_mutex);
      for(;;) {
        mutex_lock(htrdr->mpi_mutex);
        MPI(Test(&req, &complete, MPI_STATUS_IGNORE));
        mutex_unlock(htrdr->mpi_mutex);
        if(complete) break;
        nanosleep(&t, NULL);
      }
    }
  }
}

void
print_mpi_progress(struct htrdr* htrdr, const enum htrdr_mpi_message msg)
{
  ASSERT(htrdr && htrdr->mpi_rank == 0);

  if(htrdr->mpi_nprocs == 1) {
    switch(msg) {
      case HTRDR_MPI_PROGRESS_RENDERING:
        htrdr_fprintf(htrdr, stderr, "\033[2K\rRendering: %3d%%",
          htrdr->mpi_progress_render[0]);
        break;
      default: FATAL("Unreachable code.\n"); break;
    }
    htrdr_fflush(htrdr, stderr);
  } else {
    int iproc;
    FOR_EACH(iproc, 0, htrdr->mpi_nprocs) {
      switch(msg) {
        case HTRDR_MPI_PROGRESS_RENDERING:
          htrdr_fprintf(htrdr, stderr,
            "\033[2K\rProcess %d -- rendering: %3d%%%c",
            iproc, htrdr->mpi_progress_render[iproc],
            iproc == htrdr->mpi_nprocs - 1 ? '\r' : '\n');
          break;
        default: FATAL("Unreachable code.\n"); break;
      }
    }
  }
}

void
clear_mpi_progress(struct htrdr* htrdr, const enum htrdr_mpi_message msg)
{
  ASSERT(htrdr);
  (void)msg;
  if(htrdr->mpi_nprocs > 1) {
    htrdr_fprintf(htrdr, stderr, "\033[%dA", htrdr->mpi_nprocs-1);
  }
}

int
total_mpi_progress(const struct htrdr* htrdr, const enum htrdr_mpi_message msg)
{
  const int* progress = NULL;
  int total = 0;
  int iproc;
  ASSERT(htrdr && htrdr->mpi_rank == 0);

  switch(msg) {
    case HTRDR_MPI_PROGRESS_RENDERING:
      progress = htrdr->mpi_progress_render;
      break;
    default: FATAL("Unreachable code.\n"); break;
  }

  FOR_EACH(iproc, 0, htrdr->mpi_nprocs) {
    total += progress[iproc];
  }
  total = total / htrdr->mpi_nprocs;
  return total;
}

