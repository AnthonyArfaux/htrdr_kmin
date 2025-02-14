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

#include "core/htrdr.h"
#include "core/htrdr_c.h"
#include "core/htrdr_buffer.h"
#include "core/htrdr_proc_work.h"
#include "core/htrdr_log.h"
#include "core/htrdr_solve_buffer.h"

#include <rsys/clock_time.h>
#include <rsys/cstr.h>
#include <rsys/list.h>
#include <rsys/math.h>
#include <rsys/mutex.h>
#include <rsys/ref_count.h>

#include <star/ssp.h>

#include <mpi.h>
#include <omp.h>

#define CHUNK_SIZE 32 /* Number of items in one chunk */

/* Collection of items */
struct chunk {
  struct list_node node;
  struct mem_allocator* allocator;
  ref_T ref;

  struct chunk_data {
    size_t item_sz; /* Size of an item */
    size_t item_al; /* Item alignment */
    uint64_t index; /* Chunk index in chunk space */
    /* Simulate the flexible array member of the C99 standard */
    char ALIGN(16) items[1/*Dummy*/];
  } data;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static INLINE res_T
check_solve_buffer_args(const struct htrdr_solve_buffer_args* args)
{
  if(!args) return RES_BAD_ARG;

  /* A functor must be defined */
  if(!args->solve_item) return RES_BAD_ARG;

  /* The number of realisations cannot be null */
  if(!args->nrealisations) return RES_BAD_ARG;

  /* The buffer must in one-dimensional */
  if(args->buffer_layout.height != 1) return RES_BAD_ARG;

  /* Check buffer layout consistency */
  return htrdr_buffer_layout_check(&args->buffer_layout);
}

static INLINE void
release_chunk(ref_T* ref)
{
  struct chunk* chunk = CONTAINER_OF(ref, struct chunk, ref);
  ASSERT(ref);
  MEM_RM(chunk->allocator, chunk);
}

static INLINE void
chunk_ref_get(struct chunk* chunk)
{
  ASSERT(chunk);
  ref_get(&chunk->ref);
}

static INLINE void
chunk_ref_put(struct chunk* chunk)
{
  ASSERT(chunk);
  ref_put(&chunk->ref, release_chunk);
}

static FINLINE struct chunk*
chunk_create
  (struct mem_allocator* allocator,
   const size_t item_sz, /* Size in bytes */
   const size_t item_al, /* Alignment in bytes */
   const uint64_t ichunk) /* Chunk index */
{
  struct chunk* chunk = NULL;
  const size_t header_sz = sizeof(*chunk) - 1/*dummy octet in flexible array*/;
  const size_t buf_sz = CHUNK_SIZE*item_sz;

  /* Pre conditions */
  ASSERT(allocator);
  ASSERT(IS_ALIGNED(item_sz, item_al));
  ASSERT(IS_POW2(item_al));

  chunk = MEM_ALLOC_ALIGNED(allocator, header_sz + buf_sz, 16);
  if(!chunk) goto error;
  ref_init(&chunk->ref);
  list_init(&chunk->node);
  chunk->allocator = allocator;
  chunk->data.item_sz = item_sz;
  chunk->data.item_al = item_al;
  chunk->data.index = ichunk;
  CHK(IS_ALIGNED(chunk->data.items, item_al));

exit:
  return chunk;
error:
  if(chunk) { chunk_ref_put(chunk); chunk = NULL; }
  goto exit;
}

static FINLINE void*
chunk_at(struct chunk* chunk, const size_t i)
{
  ASSERT(chunk && i < CHUNK_SIZE);
  return (void*)(chunk->data.items + i*chunk->data.item_sz);
}

static void
release_chunk_list(struct list_node* chunks)
{
  struct list_node* node = NULL;
  struct list_node* tmp = NULL;
  ASSERT(chunks);

  LIST_FOR_EACH_SAFE(node, tmp, chunks) {
    struct chunk* chunk = CONTAINER_OF(node, struct chunk, node);
    list_del(node);
    chunk_ref_put(chunk);
  }
}

static INLINE uint64_t
get_chunk
  (struct htrdr* htrdr,
   struct ssp_rng* rng,
   struct proc_work* work)
{
  uint64_t ichunk = CHUNK_ID_NULL;

  /* Make the function critical, as the entire process must be executed
   * atomically. Indeed, the first thread to query an invalid chunk steals work
   * from other processes. So, when the function exits, the other threads will
   * have valid chunks */
  #pragma omp critical
  {
    ichunk = proc_work_get_chunk(work);
    if(ichunk == CHUNK_ID_NULL) { /* No more work on this process */
      size_t nthieves = 0;

      proc_work_reset(work);
      nthieves = mpi_steal_work(htrdr, rng, work);
      if(nthieves != 0) {
        ichunk = proc_work_get_chunk(work);
      }
    }
  }

  return ichunk;
}

static INLINE void
status_update(struct htrdr* htrdr, const int32_t progress)
{
  ASSERT(htrdr);

  #pragma omp critical
  if(progress > htrdr->mpi_progress_render[0]) {
    htrdr->mpi_progress_render[0] = progress;

    /* Print update on master process */
    if(htrdr->mpi_rank == 0) {
      update_mpi_progress(htrdr, HTRDR_MPI_PROGRESS_RENDERING);

      /* Send the progress percentage to the master process */
    } else {
      send_mpi_progress(htrdr, HTRDR_MPI_PROGRESS_RENDERING, progress);
    }
  }
}

static struct ssp_rng*
rng_create
  (struct htrdr* htrdr,
   const size_t ithread,
   const size_t nchunks,
   const uint64_t ichunk)
{
  struct ssp_rng_proxy_create2_args args = SSP_RNG_PROXY_CREATE2_ARGS_NULL;
  struct ssp_rng_proxy* proxy = NULL;
  struct ssp_rng* rng = NULL;
  struct mem_allocator* allocator = NULL;
  ASSERT(htrdr);

  allocator = htrdr_get_thread_allocator(htrdr, (size_t)ithread);

  args.type = SSP_RNG_THREEFRY;
  args.sequence_offset = RNG_SEQUENCE_SIZE * (size_t)ichunk;
  args.sequence_size = RNG_SEQUENCE_SIZE;
  args.sequence_pitch = RNG_SEQUENCE_SIZE * nchunks;
  args.nbuckets = 1;
  SSP(rng_proxy_create2(allocator, &args, &proxy));
  SSP(rng_proxy_create_rng(proxy, 0, &rng));
  SSP(rng_proxy_ref_put(proxy));

  return rng;
}

static void
write_chunk_data
  (struct htrdr* htrdr,
   struct htrdr_buffer* buf,
   const struct chunk_data* chunk_data)
{
  struct htrdr_buffer_layout layout = HTRDR_BUFFER_LAYOUT_NULL;
  char* mem = NULL;
  size_t iitem = 0;
  size_t nitems = 0;

  ASSERT(htrdr && buf && chunk_data);
  ASSERT(chunk_data->index != CHUNK_ID_NULL);
  (void)htrdr;

  htrdr_buffer_get_layout(buf, &layout);
  ASSERT(layout.height == 1);
  ASSERT(layout.elmt_size == chunk_data->item_sz);

  /* Calculate the index of the first item to write */
  iitem = chunk_data->index * CHUNK_SIZE;

  /* Define the number of items to write into the buffer */
  nitems = MMIN(iitem + CHUNK_SIZE, layout.width) - iitem;

  /* Calculate destination address for writing chunk data */
  mem = htrdr_buffer_get_data(buf);
  mem = mem + iitem*layout.elmt_size;

  /* Write the chunk items into the buffer */
  memcpy(mem, chunk_data->items, nitems*layout.elmt_size);
}

static res_T
mpi_gather_chunks
  (struct htrdr* htrdr,
   const struct htrdr_buffer_layout* layout,
   struct htrdr_buffer* buf, /* Can be NULL for non master processes */
   const size_t nchunks,
   struct list_node* chunks)
{
  size_t msg_sz = 0;
  struct list_node* node = NULL;
  struct chunk* chunk = NULL; /* Temporary chunk */
  res_T res = RES_OK;

  /* Pre conditions */
  ASSERT(htrdr && layout && chunks);
  ASSERT(layout->height == 1);
  ASSERT(htrdr_buffer_layout_check(layout) == RES_OK);

  /* Compute the size of chunk data */
  msg_sz = sizeof(struct chunk_data)
    + CHUNK_SIZE*layout->elmt_size
    -1 /* Dummy octet of the flexible array */;
  ASSERT(msg_sz <= INT_MAX);

  /* Non master process: Send the computed chunk to the master process */
  if(htrdr->mpi_rank != 0) {
    LIST_FOR_EACH(node, chunks) {
      struct chunk* c = CONTAINER_OF(node, struct chunk, node);
      MPI(Send(&c->data, (int)msg_sz, MPI_CHAR, 0/*Master process*/,
        HTRDR_MPI_CHUNK_DATA, MPI_COMM_WORLD));
    }

  /* Master rrocess */
  } else {
    size_t ichunk = 0;
    ASSERT(buf);

#ifndef NDEBUG
    /* Check that the submitted buffer layout matches the buffer layout */
    {
      struct htrdr_buffer_layout buf_layout = HTRDR_BUFFER_LAYOUT_NULL;
      htrdr_buffer_get_layout(buf, &buf_layout);
      ASSERT(htrdr_buffer_layout_eq(&buf_layout, layout));
    }
#endif

    /* Write data for chunks resolved by the master process */
    LIST_FOR_EACH(node, chunks) {
      struct chunk* c = CONTAINER_OF(node, struct chunk, node);
      write_chunk_data(htrdr, buf, &c->data);
      ++ichunk;
    }

    /* There are chunks unresolved by the master process */
    if(ichunk != nchunks) {
      ASSERT(htrdr->mpi_nprocs > 1);

      /* Create a temporary chunk to receive the chunk data computed by the
       * non master processes */
      chunk = chunk_create
        (htrdr->allocator,
         layout->elmt_size,
         layout->alignment,
         CHUNK_ID_NULL); /* Dummy chunk id that is going to be overwritten */
      if(!chunk) goto error;
    }

    /* Gather chunks sent by non master processes */
    FOR_EACH(ichunk, ichunk, nchunks) {
      MPI(Recv(&chunk->data, (int)msg_sz, MPI_CHAR, MPI_ANY_SOURCE,
        HTRDR_MPI_CHUNK_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE));
      write_chunk_data(htrdr, buf, &chunk->data);
    }
  }

exit:
  if(chunk) chunk_ref_put(chunk);
  return res;
error:
  htrdr_log_err(htrdr, "Error while gathering results -- %s\n",
    res_to_cstr(res));
  goto exit;
}

static res_T
solve_chunk
  (struct htrdr* htrdr,
   const struct htrdr_solve_buffer_args* args,
   const size_t ithread,
   const uint64_t ichunk,
   struct ssp_rng* rng,
   struct chunk* chunk)
{
  struct htrdr_solve_item_args item_args = HTRDR_SOLVE_ITEM_ARGS_NULL;
  size_t i = 0;
  size_t nitems = 0;

  ASSERT(htrdr && args && rng && chunk);
  ASSERT(args->buffer_layout.height == 1);
  ASSERT(ichunk * CHUNK_SIZE < args->buffer_layout.width);

  nitems = args->buffer_layout.width;

  /* Setup the shared item arguments */
  item_args.rng = rng;
  item_args.nrealisations = args->nrealisations;
  item_args.ithread = ithread;
  item_args.context = args->context;

  FOR_EACH(i, 0, CHUNK_SIZE) {
    void* item = NULL;
    size_t item_id = ichunk*CHUNK_SIZE + i;

    if(item_id >= nitems) {
      /* The item is out of the range,
       * i.e. we've reached the total number of items to solve */
      break;
    }

    /* Fetch the item */
    item = chunk_at(chunk, i);

    /* Setup the item index */
    item_args.item_id = item_id;

    /* Solve the item */
    args->solve_item(htrdr, &item_args, item);
  }
  return RES_OK;
}

static res_T
solve_buffer
  (struct htrdr* htrdr,
   const struct htrdr_solve_buffer_args* args,
   const size_t nchunks, /* Total #chunks distributed between processes */
   struct proc_work* work,
   struct list_node* chunks)
{
  struct ssp_rng* rng_proc = NULL; /* Used to sample a working process */
  size_t nchunks_proc = 0; /* #chunks of the process */
  size_t nthreads = 0; /* #threads to use */
  ATOMIC nchunks_solved = 0; /* #chunks solved bu the process */
  ATOMIC res = RES_OK;

  /* Pre-conditions */
  ASSERT(htrdr && args && work && chunks);

  res = ssp_rng_create(htrdr->allocator, SSP_RNG_MT19937_64, &rng_proc);
  if(res != RES_OK) goto error;

  nchunks_proc = proc_work_get_nchunks(work);
  nthreads = MMIN(htrdr->nthreads, nchunks_proc);

  /* The process is not considered as a working process for himself */
  htrdr->mpi_working_procs[htrdr->mpi_rank] = 0;
  --htrdr->mpi_nworking_procs;

  omp_set_num_threads((int)nthreads);
  #pragma omp parallel
  for(;;) {
    /* Chunk */
    uint64_t ichunk = get_chunk(htrdr, rng_proc, work);
    struct chunk* chunk = NULL;

    /* Miscellaneous */
    const size_t ithread = (size_t)omp_get_thread_num();
    struct ssp_rng* rng = NULL;
    size_t n = 0;
    int32_t pcent = 0;
    res_T res_local = RES_OK;

    if(ichunk == CHUNK_ID_NULL) break; /* No more work */

    chunk = chunk_create
      (htrdr->allocator,
       args->buffer_layout.elmt_size,
       args->buffer_layout.alignment,
       ichunk);
    if(!chunk) {
      ATOMIC_SET(&res, RES_MEM_ERR);
      break;
    }

    #pragma omp critical
    list_add_tail(chunks, &chunk->node); /* Register the chunk */

    rng = rng_create(htrdr, ithread, nchunks, ichunk);
    res_local = solve_chunk(htrdr, args, ithread, ichunk, rng, chunk);

    SSP(rng_ref_put(rng));
    if(res_local != RES_OK) {
      ATOMIC_SET(&res, res_local);
      break;
    }

    /* Status update */
    n = (size_t)ATOMIC_INCR(&nchunks_solved);
    pcent = (int32_t)((double)n * 100.0 / (double)nchunks_proc + 0.5/*round*/);
    status_update(htrdr, pcent);
  }

  if(ATOMIC_GET(&res) != RES_OK) goto error;

  /* Asynchronously wait for processes completion. Use an asynchronous barrier to
   * avoid a dead lock with the `mpi_probe_thieves' thread that requires also
   * the `mpi_mutex'. */
  {
    MPI_Request req;

    mutex_lock(htrdr->mpi_mutex);
    MPI(Ibarrier(MPI_COMM_WORLD, &req));
    mutex_unlock(htrdr->mpi_mutex);

    mpi_wait_for_request(htrdr, &req);
  }

exit:
  if(rng_proc) SSP(rng_ref_put(rng_proc));
  return (res_T)res;
error:
  htrdr_log_err(htrdr, "Error while solving buffer -- %s\n",
    res_to_cstr((res_T)res));
  goto exit;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_solve_buffer
  (struct htrdr* htrdr,
   const struct htrdr_solve_buffer_args* args,
   struct htrdr_buffer* buf)
{
  /* Time registration */
  char strbuf[128];
  struct time t0, t1;

  /* Chunks */
  struct list_node chunks; /* List of solved chunks */
  size_t nchunks = 0; /* Overall number of chunks */
  size_t nchunks_proc = 0; /* #chunks for the current proc*/
  size_t nchunks_remain = 0; /* Remaining #chunks to distribute between procs */

  /* Miscellaneous */
  struct proc_work work;
  size_t nitems = 0; /* Number of Monte Carlo estimations */
  size_t i = 0;
  ATOMIC probe_thieves = 1; /* Boolean that controls thieves' polling */

  res_T res = RES_OK;
  ASSERT(htrdr && check_solve_buffer_args(args) == RES_OK);
  ASSERT(htrdr->mpi_rank != 0 || buf);

  list_init(&chunks);
  proc_work_init(htrdr->allocator, &work);

  nitems = args->buffer_layout.width;
  nchunks = (nitems + (CHUNK_SIZE-1)/*ceil*/) / CHUNK_SIZE;
  nchunks_proc = nchunks / (size_t)htrdr->mpi_nprocs;
  nchunks_remain = nchunks - nchunks_proc*(size_t)htrdr->mpi_nprocs;

  /* Distribute the remaining chunks among the processes. Each process whose
   * rank is lower than the number of remaining chunks takes an additional
   * chunk */
  if((size_t)htrdr->mpi_rank < nchunks_remain) {
    ++nchunks_proc;
  }

  /* Register the list of chunks to be processed by the current process */
  FOR_EACH(i, 0, nchunks_proc) {
    size_t ichunk = i * (size_t)htrdr->mpi_nprocs + (size_t)htrdr->mpi_rank;
    proc_work_add_chunk(&work, (uint64_t)ichunk);
  }

  /* On the master process, request and print the progress report, since the
   * other processes have been able to start the calculation */
  if(htrdr->mpi_rank == 0) {
    fetch_mpi_progress(htrdr, HTRDR_MPI_PROGRESS_RENDERING);
    print_mpi_progress(htrdr, HTRDR_MPI_PROGRESS_RENDERING);
  }

  /* Start of calculation time recording */
  time_current(&t0);

  /* Enable nested threads to enable parallelization of the solve function */
  omp_set_nested(1);

  #pragma omp parallel sections num_threads(2)
  {
    /* Polling of steal queries */
    #pragma omp section
    mpi_probe_thieves(htrdr, &work, &probe_thieves);

    #pragma omp section
    {
      solve_buffer(htrdr, args, nchunks, &work, &chunks);
      /* The process have no more work to do. Stop polling for thieves */
      ATOMIC_SET(&probe_thieves, 0);
    }
  }

  if(htrdr->mpi_rank == 0) {
    update_mpi_progress(htrdr, HTRDR_MPI_PROGRESS_RENDERING);
    htrdr_log(htrdr, "\n"); /* Add a new line after the progress statuses */
  }

  /* Stop time recording */
  time_sub(&t0, time_current(&t1), &t0);
  time_dump(&t0, TIME_ALL, NULL, strbuf, sizeof(strbuf));
  htrdr_log(htrdr, "Calculation time: %s\n", strbuf);

  /* Gather chunks on master process */
  time_current(&t0);
  res = mpi_gather_chunks(htrdr, &args->buffer_layout, buf, nchunks, &chunks);
  if(res != RES_OK) goto error;
  time_sub(&t0, time_current(&t1), &t0);
  time_dump(&t0, TIME_ALL, NULL, strbuf, sizeof(strbuf));
  htrdr_log(htrdr, "Time to gather results: %s\n", strbuf);

exit:
  release_chunk_list(&chunks);
  proc_work_release(&work);
  return res;
error:
  goto exit;
}
