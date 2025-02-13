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

#define _POSIX_C_SOURCE 200112L /* nanosleep */

#include "core/htrdr_c.h"
#include "core/htrdr_proc_work.h"

#include <star/ssp.h>

#include <rsys/mutex.h>

#include <time.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
/* Return the rank of a working process */
static int
sample_working_process(struct htrdr* htrdr, struct ssp_rng* rng)
{
  int iproc, i;
  int dst_rank;
  ASSERT(htrdr && rng && htrdr->mpi_nworking_procs);

  /* Sample the index of the 1st active process */
  iproc = (int)(ssp_rng_canonical(rng) * (double)htrdr->mpi_nworking_procs);

  /* Find the rank of the sampled active process. Use a simple linear search
   * since the overall number of processes should be quite low; at most few
   * dozens.  */
  i = 0;
  FOR_EACH(dst_rank, 0, htrdr->mpi_nprocs) {
    if(htrdr->mpi_working_procs[dst_rank] == 0) continue; /* Inactive process */
    if(i == iproc) break; /* The rank of the sampled process is found */
    ++i;
  }
  ASSERT(dst_rank < htrdr->mpi_nprocs);
  return dst_rank;
}


/*******************************************************************************
 * Local functions
 ******************************************************************************/
void
proc_work_init(struct mem_allocator* allocator, struct proc_work* work)
{
  ASSERT(work);
  darray_u64_init(allocator, &work->chunks);
  work->index = 0;
  CHK(work->mutex = mutex_create());
}

void
proc_work_release(struct proc_work* work)
{
  darray_u64_release(&work->chunks);
  mutex_destroy(work->mutex);
}

void
proc_work_reset(struct proc_work* work)
{
  ASSERT(work);
  mutex_lock(work->mutex);
  darray_u64_clear(&work->chunks);
  work->index = 0;
  mutex_unlock(work->mutex);
}

void
proc_work_add_chunk(struct proc_work* work, const size_t ichunk)
{
  mutex_lock(work->mutex);
  CHK(darray_u64_push_back(&work->chunks, &ichunk) == RES_OK);
  mutex_unlock(work->mutex);
}

uint64_t
proc_work_get_chunk(struct proc_work* work)
{
  uint64_t ichunk = CHUNK_ID_NULL;
  ASSERT(work);

  mutex_lock(work->mutex);
  if(work->index >= darray_u64_size_get(&work->chunks)) {
    ichunk = CHUNK_ID_NULL;
  } else {
    ichunk = darray_u64_cdata_get(&work->chunks)[work->index];
    ++work->index;
  }
  mutex_unlock(work->mutex);
  return ichunk;
}

size_t
proc_work_get_nchunks(struct proc_work* work)
{
  size_t sz = 0;
  ASSERT(work);

  mutex_lock(work->mutex);
  sz = darray_u64_size_get(&work->chunks);
  mutex_unlock(work->mutex);
  return sz;
}

void
mpi_wait_for_request(struct htrdr* htrdr, MPI_Request* req)
{
  ASSERT(htrdr && req);

  /* Wait for process synchronisation */
  for(;;) {
    struct timespec t;
    int complete;
    t.tv_sec = 0;
    t.tv_nsec = 10000000; /* 10ms */

    mutex_lock(htrdr->mpi_mutex);
    MPI(Test(req, &complete, MPI_STATUS_IGNORE));
    mutex_unlock(htrdr->mpi_mutex);
    if(complete) break;

    nanosleep(&t, NULL);
  }
}

void
mpi_probe_thieves
  (struct htrdr* htrdr,
   struct proc_work* work,
   ATOMIC* probe_thieves)
{
  uint64_t chunks[UINT8_MAX];
  struct timespec t;
  ASSERT(htrdr && work && probe_thieves);

  if(htrdr->mpi_nprocs == 1) /* The process is alone. No thief is possible */
    return;

  t.tv_sec = 0;

  /* Protect MPI calls of multiple invocations from concurrent threads */
  #define P_MPI(Func) {                                                        \
    mutex_lock(htrdr->mpi_mutex);                                              \
    MPI(Func);                                                                 \
    mutex_unlock(htrdr->mpi_mutex);                                            \
  } (void)0

  while(ATOMIC_GET(probe_thieves)) {
    MPI_Status status;
    uint8_t i = 0;
    int msg = 0;

    /* Probe if a steal request was submitted by any processes */
    P_MPI(Iprobe(MPI_ANY_SOURCE, HTRDR_MPI_STEAL_REQUEST, MPI_COMM_WORLD, &msg,
      &status));

    if(msg) { /* A steal request was posted */
      MPI_Request req;
      uint8_t nchunks_to_steal;

      /* Asynchronously receive the steal request */
      P_MPI(Irecv(&nchunks_to_steal, 1, MPI_UINT8_T, status.MPI_SOURCE,
        HTRDR_MPI_STEAL_REQUEST, MPI_COMM_WORLD, &req));

      /* Wait for the completion of the steal request */
      mpi_wait_for_request(htrdr, &req);

      /* Thief some chunks */
      FOR_EACH(i, 0, nchunks_to_steal) {
        chunks[i] = proc_work_get_chunk(work);
      }
      P_MPI(Send(&chunks, nchunks_to_steal, MPI_UINT64_T, status.MPI_SOURCE,
        HTRDR_MPI_WORK_STEALING, MPI_COMM_WORLD));
    }

    /* Don't constantly check for thieves */
    t.tv_nsec = 500000000; /* 500ms */
    nanosleep(&t, NULL);
  }
  #undef P_MPI
}

/* Return the number of stolen tiles */
size_t
mpi_steal_work
  (struct htrdr* htrdr,
   struct ssp_rng* rng,
   struct proc_work* work)
{
  MPI_Request req;
  size_t nthieves = 0;
  uint64_t chunks[UINT8_MAX]; /* Index of the stolen chunks */
  int proc_to_steal; /* Rank of the process to steal */

  /* Empircally set the number of chunks to steal */
  const uint8_t nchunks_to_steal = MMIN((uint8_t)(htrdr->nthreads*2), 16);
  uint8_t i = 0;

  ASSERT(htrdr && rng && work && htrdr->nthreads < UINT8_MAX);

  /* Protect MPI calls of multiple invocations from concurrent threads */
  #define P_MPI(Func) {                                                        \
    mutex_lock(htrdr->mpi_mutex);                                              \
    MPI(Func);                                                                 \
    mutex_unlock(htrdr->mpi_mutex);                                            \
  } (void)0

  /* No more working process => nothing to steal */
  if(!htrdr->mpi_nworking_procs) return 0;

  /* Sample a process to steal */
  proc_to_steal = sample_working_process(htrdr, rng);

  /* Send a steal request to the sampled process and wait for a response */
  P_MPI(Send(&nchunks_to_steal, 1, MPI_UINT8_T, proc_to_steal,
    HTRDR_MPI_STEAL_REQUEST, MPI_COMM_WORLD));

  /* Receive the stolen chunks from the sampled process */
  P_MPI(Irecv(chunks, nchunks_to_steal, MPI_UINT64_T, proc_to_steal,
    HTRDR_MPI_WORK_STEALING, MPI_COMM_WORLD, &req));

  mpi_wait_for_request(htrdr, &req);

  FOR_EACH(i, 0, nchunks_to_steal) {

    if(chunks[i] != CHUNK_ID_NULL) {
      /* Save stolen chunk in job list */
      proc_work_add_chunk(work, chunks[i]);
      ++nthieves;

    } else {
      /* The process has returned at least one invalid chunk,
       * i.e. it has nothing further to do.
       * Remove it from the working process */
      ASSERT(htrdr->mpi_working_procs[proc_to_steal] != 0);
      htrdr->mpi_working_procs[proc_to_steal] = 0;
      htrdr->mpi_nworking_procs--;

      break; /* No more to steal */
    }
  }
  #undef P_MPI
  return nthieves;
}
