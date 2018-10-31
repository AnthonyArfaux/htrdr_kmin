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

#include "htrdr.h"
#include "htrdr_c.h"
#include "htrdr_buffer.h"
#include "htrdr_camera.h"
#include "htrdr_sky.h"
#include "htrdr_solve.h"

#include <rsys/cstr.h>
#include <rsys/dynamic_array_u32.h>
#include <rsys/math.h>
#include <rsys/mutex.h>
#include <star/ssp.h>

#include <omp.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>

#define RNG_SEQUENCE_SIZE 100000

#define TILE_SIZE 32 /* Definition in X & Y of a tile */
STATIC_ASSERT(IS_POW2(TILE_SIZE), TILE_SIZE_must_be_a_power_of_2);

/* Overall work of a process */
struct proc_work {
  struct mutex* mutex;
  struct darray_u32 tiles; /* #tiles to render */
  size_t itile; /* Next tile to render in the above list of tiles */
};

#define TILE_MCODE_NULL UINT32_MAX

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static FINLINE uint16_t
morton2D_decode(const uint32_t u32)
{
  uint32_t x = u32 & 0x55555555;
  x = (x | (x >> 1)) & 0x33333333;
  x = (x | (x >> 2)) & 0x0F0F0F0F;
  x = (x | (x >> 4)) & 0x00FF00FF;
  x = (x | (x >> 8)) & 0x0000FFFF;
  return (uint16_t)x;
}

static FINLINE uint32_t
morton2D_encode(const uint16_t u16)
{
  uint32_t u32 = u16;
  u32 = (u32 | (u32 << 8)) & 0x00FF00FF;
  u32 = (u32 | (u32 << 4)) & 0X0F0F0F0F;
  u32 = (u32 | (u32 << 2)) & 0x33333333;
  u32 = (u32 | (u32 << 1)) & 0x55555555;
  return u32;
}

static INLINE void
proc_work_init(struct mem_allocator* allocator, struct proc_work* work)
{
  ASSERT(work);
  darray_u32_init(allocator, &work->tiles);
  work->itile = 0;
  CHK(work->mutex = mutex_create());
}

static INLINE void
proc_work_release(struct proc_work* work)
{
  darray_u32_release(&work->tiles);
  mutex_destroy(work->mutex);
}

static INLINE void
proc_work_reset(struct proc_work* work)
{
  ASSERT(work);
  mutex_lock(work->mutex);
  darray_u32_clear(&work->tiles);
  work->itile = 0;
  mutex_unlock(work->mutex);
}

static INLINE void
proc_work_add_tile(struct proc_work* work, const uint32_t mcode)
{
  mutex_lock(work->mutex);
  CHK(darray_u32_push_back(&work->tiles, &mcode) == RES_OK);
  mutex_unlock(work->mutex);
}

static INLINE uint32_t
proc_work_get_tile(struct proc_work* work)
{
  uint32_t mcode;
  ASSERT(work);
  mutex_lock(work->mutex);
  if(work->itile >= darray_u32_size_get(&work->tiles)) {
    mcode = TILE_MCODE_NULL;
  } else {
    mcode = darray_u32_cdata_get(&work->tiles)[work->itile];
    ++work->itile;
  }
  mutex_unlock(work->mutex);
  return mcode;
}

static INLINE size_t
proc_work_get_ntiles(struct proc_work* work)
{
  size_t sz = 0;
  ASSERT(work);
  mutex_lock(work->mutex);
  sz = darray_u32_size_get(&work->tiles);
  mutex_unlock(work->mutex);
  return sz;
}

static void
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

static void
mpi_probe_thieves
  (struct htrdr* htrdr,
   struct proc_work* work,
   ATOMIC* probe_thieves)
{
  uint32_t tiles[UINT8_MAX];
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
    size_t itile;
    int msg;

    /* Probe if a steal request was submitted by any processes */
    P_MPI(Iprobe(MPI_ANY_SOURCE, HTRDR_MPI_STEAL_REQUEST, MPI_COMM_WORLD, &msg,
      &status));

    if(msg) { /* A steal request was posted */
      MPI_Request req;
      uint8_t ntiles_to_steal;

      /* Asynchronously receive the steal request */
      P_MPI(Irecv(&ntiles_to_steal, 1, MPI_UINT8_T, status.MPI_SOURCE,
        HTRDR_MPI_STEAL_REQUEST, MPI_COMM_WORLD, &req));

      /* Wait for the completion of the steal request */
      mpi_wait_for_request(htrdr, &req);

      /* Thief some tiles */
      FOR_EACH(itile, 0, ntiles_to_steal) {
        tiles[itile] = proc_work_get_tile(work);
      }
      P_MPI(Send(&tiles, ntiles_to_steal, MPI_UINT32_T, status.MPI_SOURCE,
        HTRDR_MPI_WORK_STEALING, MPI_COMM_WORLD));
    }
    t.tv_nsec = 500000000; /* 500ms */
    nanosleep(&t, NULL);
  }
  #undef P_MPI
}

static int
mpi_sample_working_process(struct htrdr* htrdr, struct ssp_rng* rng)
{
  int iproc, i;
  int dst_rank;
  ASSERT(htrdr && rng && htrdr->mpi_nworking_procs);

  /* Sample the index of the 1st active process */
  iproc = (int)(ssp_rng_canonical(rng) * (double)htrdr->mpi_nworking_procs);

  /* Find the rank of the sampled active process. Use a simple linear search
   * since the overall number of process should be quite low; at most few
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

/* Return the number of stolen tiles */
static size_t
mpi_steal_work
  (struct htrdr* htrdr,
   struct ssp_rng* rng,
   struct proc_work* work)
{
  MPI_Request req;
  size_t itile;
  size_t nthieves = 0;
  uint32_t tiles[UINT8_MAX]; /* Morton code of the stolen tile */
  int proc_to_steal; /* Process to steal */
  uint8_t ntiles_to_steal = (uint8_t)(htrdr->nthreads*2);
  ASSERT(htrdr && rng && work && htrdr->nthreads < UINT8_MAX);

  /* Protect MPI calls of multiple invocations from concurrent threads */
  #define P_MPI(Func) {                                                        \
    mutex_lock(htrdr->mpi_mutex);                                              \
    MPI(Func);                                                                 \
    mutex_unlock(htrdr->mpi_mutex);                                            \
  } (void)0

  /* No more working process => nohting to steal */
  if(!htrdr->mpi_nworking_procs) return 0;

  /* Sample a process to steal */
  proc_to_steal = mpi_sample_working_process(htrdr, rng);

  /* Send a steal request to the sampled process and wait for a response */
  P_MPI(Send(&ntiles_to_steal, 1, MPI_UINT8_T, proc_to_steal,
    HTRDR_MPI_STEAL_REQUEST, MPI_COMM_WORLD));

  /* Receive the stolen tile from the sampled process */
  P_MPI(Irecv(tiles, ntiles_to_steal, MPI_UINT32_T, proc_to_steal,
    HTRDR_MPI_WORK_STEALING, MPI_COMM_WORLD, &req));

  mpi_wait_for_request(htrdr, &req);

  FOR_EACH(itile, 0, ntiles_to_steal) {
    if(tiles[itile] == TILE_MCODE_NULL) {
      ASSERT(htrdr->mpi_working_procs[proc_to_steal] != 0);
      htrdr->mpi_working_procs[proc_to_steal] = 0;
      htrdr->mpi_nworking_procs--;
      break;
    }
    proc_work_add_tile(work, tiles[itile]);
    ++nthieves;
  }
  #undef P_MPI
  return nthieves;
}

static res_T
mpi_gather_buffer
  (struct htrdr* htrdr,
   struct htrdr_buffer* buf)
{
  struct htrdr_buffer_layout layout;
  struct htrdr_accum* gathered_accums = NULL;
  size_t x, y;
  int iproc;
  res_T res = RES_OK;
  ASSERT(htrdr && buf);

  /* Fetch the memory layout of the submitted buffer */
  htrdr_buffer_get_layout(buf, &layout);
  ASSERT(layout.elmt_size == sizeof(struct htrdr_accum) * 3/*#channels*/);
  ASSERT(layout.alignment <= ALIGNOF(struct htrdr_accum));

  /* The process 0 allocates the memory used to store the gathered buffer lines
   * of the MPI processes */
  if(htrdr->mpi_rank == 0) {
    gathered_accums = MEM_ALLOC
      (htrdr->allocator, layout.pitch * (size_t)htrdr->mpi_nprocs);
    if(!gathered_accums) {
      res = RES_MEM_ERR;
      htrdr_log_err(htrdr,
        "could not allocate the temporary memory for MPI gathering -- %s.\n",
        res_to_cstr(res));
      goto error;
    }
  }

  FOR_EACH(y, 0, layout.height) {
    struct htrdr_accum* buf_row_accums = (struct htrdr_accum*)
      ((char*)htrdr_buffer_get_data(buf) + y * layout.pitch);
    int err; /* MPI error */

    /* Gather the buffer lines */
    mutex_lock(htrdr->mpi_mutex);
    err = MPI_Gather(buf_row_accums, (int)layout.pitch, MPI_CHAR, gathered_accums,
      (int)layout.pitch, MPI_CHAR, 0, MPI_COMM_WORLD);
    mutex_unlock(htrdr->mpi_mutex);
    if(err != MPI_SUCCESS) {
      htrdr_log_err(htrdr,
        "could not gather the buffer line `%lu' from the group of processes -- "
        "%s.\n",
        (unsigned long)y, htrdr_mpi_error_string(htrdr, err));
      res = RES_UNKNOWN_ERR;
      goto error;
    }

    /* Accumulates the gathered lines into the buffer of the process 0 */
    if(htrdr->mpi_rank == 0) {
      memset(buf_row_accums, 0, layout.pitch);
      FOR_EACH(iproc, 0, htrdr->mpi_nprocs) {
        struct htrdr_accum* proc_accums = (struct htrdr_accum*)
          ((char*)gathered_accums + (size_t)iproc * layout.pitch);
        FOR_EACH(x, 0, layout.width * 3/*#channels*/) {
          buf_row_accums[x].sum_weights += proc_accums[x].sum_weights;
          buf_row_accums[x].sum_weights_sqr += proc_accums[x].sum_weights_sqr;
          buf_row_accums[x].nweights += proc_accums[x].nweights;
          buf_row_accums[x].nfailures += proc_accums[x].nfailures;
        }
      }
    }
  }

exit:
  if(gathered_accums) MEM_RM(htrdr->allocator, gathered_accums);
  return res;
error:
  goto exit;
}

static res_T
draw_tile
  (struct htrdr* htrdr,
   const size_t ithread,
   const int64_t tile_mcode, /* For debug only */
   const size_t tile_org[2], /* Origin of the tile in pixel space */
   const size_t tile_sz[2], /* Definition of the tile */
   const double pix_sz[2], /* Size of a pixel in the normalized image plane */
   const struct htrdr_camera* cam,
   const size_t spp, /* #samples per pixel */
   struct ssp_rng* rng,
   struct htrdr_buffer* buf)
{
  size_t npixels;
  size_t mcode; /* Morton code of tile pixel */
  ASSERT(htrdr && tile_org && tile_sz && pix_sz && cam && spp && buf);
  (void)tile_mcode;
  /* Adjust the #pixels to process them wrt a morton order */
  npixels = round_up_pow2(MMAX(tile_sz[0], tile_sz[1]));
  npixels *= npixels;

  FOR_EACH(mcode, 0, npixels) {
    struct htrdr_accum* pix_accums;
    size_t ipix_tile[2]; /* Pixel coord in the tile */
    size_t ipix[2]; /* Pixel coord in the buffer */
    size_t ichannel;

    ipix_tile[0] = morton2D_decode((uint32_t)(mcode>>0));
    if(ipix_tile[0] >= tile_sz[0]) continue; /* Pixel is out of tile */
    ipix_tile[1] = morton2D_decode((uint32_t)(mcode>>1));
    if(ipix_tile[1] >= tile_sz[1]) continue; /* Pixel is out of tile */

    /* Compute the pixel coordinate */
    ipix[0] = tile_org[0] + ipix_tile[0];
    ipix[1] = tile_org[1] + ipix_tile[1];

    /* Fetch and reset the pixel accumulator */
    pix_accums = htrdr_buffer_at(buf, ipix[0], ipix[1]);

    FOR_EACH(ichannel, 0, 3) {
      size_t isamp;
      pix_accums[ichannel] = HTRDR_ACCUM_NULL;

      FOR_EACH(isamp, 0, spp) {
        double pix_samp[2];
        double ray_org[3];
        double ray_dir[3];
        double weight;
        size_t iband;
        size_t iquad;

        /* Sample a position into the pixel, in the normalized image plane */
        pix_samp[0] = ((double)ipix[0] + ssp_rng_canonical(rng)) * pix_sz[0];
        pix_samp[1] = ((double)ipix[1] + ssp_rng_canonical(rng)) * pix_sz[1];

        /* Generate a ray starting from the pinhole camera and passing through the
         * pixel sample */
        htrdr_camera_ray(cam, pix_samp, ray_org, ray_dir);

        /* Sample a spectral band and a quadrature point */
        switch(ichannel) {
          case 0:
            htrdr_sky_sample_sw_spectral_data_CIE_1931_X
              (htrdr->sky, rng, &iband, &iquad);
            break;
          case 1:
            htrdr_sky_sample_sw_spectral_data_CIE_1931_Y
              (htrdr->sky, rng, &iband, &iquad);
            break;
          case 2:
            htrdr_sky_sample_sw_spectral_data_CIE_1931_Z
              (htrdr->sky, rng, &iband, &iquad);
            break;
          default: FATAL("Unreachable code.\n"); break;
        }

        /* Compute the radiance that reach the pixel through the ray */
        weight = htrdr_compute_radiance_sw
          (htrdr, ithread, rng, ray_org, ray_dir, iband, iquad);
        ASSERT(weight >= 0);

        /* Update the pixel accumulator */
        pix_accums[ichannel].sum_weights += weight;
        pix_accums[ichannel].sum_weights_sqr += weight*weight;
        pix_accums[ichannel].nweights += 1;
      }
    }
  }
  return RES_OK;
}

static res_T
draw_image
  (struct htrdr* htrdr,
   const struct htrdr_camera* cam,
   const size_t spp,
   struct ssp_rng* rng_main,
   const size_t ntiles_x,
   const size_t ntiles_y,
   const size_t ntiles_adjusted,
   struct proc_work* work,
   struct htrdr_buffer* buf)
{
  struct htrdr_buffer_layout layout;
  MPI_Request req;
  double pix_sz[2];
  size_t nthreads = 0;
  size_t nthieves = 0;
  size_t proc_ntiles = 0;
  ATOMIC nsolved_tiles = 0;
  ATOMIC res = RES_OK;
  ASSERT(htrdr && cam && spp && rng_main && ntiles_adjusted && work && buf);
  (void)ntiles_x, (void)ntiles_y;

  /* Compute the size of a pixel in the normalized image plane */
  htrdr_buffer_get_layout(buf, &layout);
  pix_sz[0] = 1.0 / (double)layout.width;
  pix_sz[1] = 1.0 / (double)layout.height;

  nthreads = htrdr->nthreads;
  proc_ntiles = proc_work_get_ntiles(work);

  /* The process is not considered as a working process for himself */
  htrdr->mpi_working_procs[htrdr->mpi_rank] = 0;
  --htrdr->mpi_nworking_procs;

  do {
    omp_set_num_threads((int)nthreads);

    #pragma omp parallel
    for(;;) {
      const int ithread = omp_get_thread_num();
      struct ssp_rng_proxy* rng_proxy = NULL;
      struct ssp_rng* rng;
      uint32_t mcode;
      size_t tile_org[2];
      size_t tile_sz[2];
      size_t n;
      res_T res_local = RES_OK;
      int32_t pcent;

      if(ATOMIC_GET(&res) != RES_OK) continue;

      mcode = proc_work_get_tile(work);
      if(mcode == TILE_MCODE_NULL) break;

      /* Decode the morton code to retrieve the tile index  */
      tile_org[0] = morton2D_decode((uint32_t)(mcode>>0));
      tile_org[1] = morton2D_decode((uint32_t)(mcode>>1));
      ASSERT(tile_org[0] < ntiles_x && tile_org[1] < ntiles_y);

      /* Define the tile origin in pixel space */
      tile_org[0] *= TILE_SIZE;
      tile_org[1] *= TILE_SIZE;

      /* Compute the size of the tile clamped by the borders of the buffer */
      tile_sz[0] = MMIN(TILE_SIZE, layout.width - tile_org[0]);
      tile_sz[1] = MMIN(TILE_SIZE, layout.height - tile_org[1]);

      SSP(rng_proxy_create2
        (&htrdr->lifo_allocators[ithread],
         &ssp_rng_threefry,
         RNG_SEQUENCE_SIZE * (size_t)mcode, /* Offset */
         RNG_SEQUENCE_SIZE, /* Size */
         RNG_SEQUENCE_SIZE * (size_t)ntiles_adjusted, /* Pitch */
         1, &rng_proxy));
      SSP(rng_proxy_create_rng(rng_proxy, 0, &rng));

      res_local = draw_tile(htrdr, (size_t)ithread, mcode, tile_org, tile_sz,
        pix_sz, cam, spp, rng, buf);

      SSP(rng_proxy_ref_put(rng_proxy));
      SSP(rng_ref_put(rng));

      if(res_local != RES_OK) {
        ATOMIC_SET(&res, res_local);
        continue;
      }

      n = (size_t)ATOMIC_INCR(&nsolved_tiles);
      pcent = (int32_t)((double)n * 100.0 / (double)proc_ntiles + 0.5/*round*/);

      #pragma omp critical
      if(pcent > htrdr->mpi_progress_render[0]) {
        htrdr->mpi_progress_render[0] = pcent;
        if(htrdr->mpi_rank == 0) {
          update_mpi_progress(htrdr, HTRDR_MPI_PROGRESS_RENDERING);
        } else { /* Send the progress percentage to the master process */
          send_mpi_progress(htrdr, HTRDR_MPI_PROGRESS_RENDERING, pcent);
        }
      }
    }

    if(ATOMIC_GET(&res) != RES_OK) goto error;

    proc_work_reset(work);
    nthieves = mpi_steal_work(htrdr, rng_main, work);
    nthreads = MMIN(nthieves, htrdr->nthreads);
  } while(nthieves);

  /* Synchronize the process */
  mutex_lock(htrdr->mpi_mutex);
  MPI(Ibarrier(MPI_COMM_WORLD, &req));
  mutex_unlock(htrdr->mpi_mutex);

  /* Wait for synchronization */
  mpi_wait_for_request(htrdr, &req);

exit:
  return (res_T)res;
error:
  goto exit;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_draw_radiance_sw
  (struct htrdr* htrdr,
   const struct htrdr_camera* cam,
   const size_t spp,
   struct htrdr_buffer* buf)
{
  struct ssp_rng* rng_main = NULL;
  size_t ntiles_x, ntiles_y, ntiles_adjusted;
  size_t itile;
  struct proc_work work;
  struct htrdr_buffer_layout layout;
  size_t proc_ntiles_adjusted;
  ATOMIC probe_thieves = 1;
  ATOMIC res = RES_OK;
  ASSERT(htrdr && cam && buf);

  proc_work_init(htrdr->allocator, &work);

  htrdr_buffer_get_layout(buf, &layout);
  ASSERT(layout.width || layout.height || layout.elmt_size);

  if(layout.elmt_size != sizeof(struct htrdr_accum[3])/*#channels*/
  || layout.alignment < ALIGNOF(struct htrdr_accum[3])) {
    htrdr_log_err(htrdr,
      "%s: invalid buffer layout. "
      "The pixel size must be the size of 3 * accumulators.\n",
      FUNC_NAME);
    res = RES_BAD_ARG;
    goto error;
  }

  res = ssp_rng_create(htrdr->allocator, &ssp_rng_mt19937_64, &rng_main);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "%s: could not create the main RNG -- %s.\n",
      FUNC_NAME, res_to_cstr((res_T)res));
    goto error;
  }

  /* Compute the overall number of tiles */
  ntiles_x = (layout.width + (TILE_SIZE-1)/*ceil*/)/TILE_SIZE;
  ntiles_y = (layout.height+ (TILE_SIZE-1)/*ceil*/)/TILE_SIZE;

  /* Adjust the #tiles for the morton-encoding procedure */
  ntiles_adjusted = round_up_pow2(MMAX(ntiles_x, ntiles_y));
  ntiles_adjusted *= ntiles_adjusted;

  /* Define the initial number of tiles of the current process */
  proc_ntiles_adjusted = ntiles_adjusted / (size_t)htrdr->mpi_nprocs;
  if(htrdr->mpi_rank == 0) { /* Affect the remaining tiles to the master proc */
    proc_ntiles_adjusted +=
      ntiles_adjusted - proc_ntiles_adjusted*(size_t)htrdr->mpi_nprocs;
  }

  /* Define the initial list of tiles of the process */
  FOR_EACH(itile, 0, proc_ntiles_adjusted) {
    uint32_t mcode;
    uint16_t tile_org[2];

    mcode = (uint32_t)itile*(uint32_t)htrdr->mpi_nprocs
          + (uint32_t)htrdr->mpi_rank;

    tile_org[0] = morton2D_decode(mcode>>0);
    if(tile_org[0] >= ntiles_x) continue;
    tile_org[1] = morton2D_decode(mcode>>1);
    if(tile_org[1] >= ntiles_y) continue;
    proc_work_add_tile(&work, mcode);
  }

  if(htrdr->mpi_rank == 0) {
    fetch_mpi_progress(htrdr, HTRDR_MPI_PROGRESS_RENDERING);
    print_mpi_progress(htrdr, HTRDR_MPI_PROGRESS_RENDERING);
  }

  omp_set_nested(1); /* Enable nested threads for draw_image */
  #pragma omp parallel sections num_threads(2)
  {
    #pragma omp section
    mpi_probe_thieves(htrdr, &work, &probe_thieves);

    #pragma omp section
    {
      draw_image(htrdr, cam, spp, rng_main, ntiles_x, ntiles_y,
        ntiles_adjusted, &work, buf);
      /* The processes have no more work to do. Stop probing for thieves */
      ATOMIC_SET(&probe_thieves, 0);
    }
  }

  if(htrdr->mpi_rank == 0) {
    update_mpi_progress(htrdr, HTRDR_MPI_PROGRESS_RENDERING);
    fprintf(stderr, "\n"); /* Add a new line after the progress statuses */
  }

  /* Gather accum buffers from the group of processes */
  res = mpi_gather_buffer(htrdr, buf);
  if(res != RES_OK) goto error;

exit:
  proc_work_release(&work);
  if(rng_main) SSP(rng_ref_put(rng_main));
  return (res_T)res;
error:
  goto exit;
}

