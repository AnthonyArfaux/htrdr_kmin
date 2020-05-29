/* Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
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

#define _POSIX_C_SOURCE 200112L /* nanosleep && nextafter */

#include "htrdr.h"
#include "htrdr_c.h"
#include "htrdr_buffer.h"
#include "htrdr_camera.h"
#include "htrdr_cie_xyz.h"
#include "htrdr_ran_lw.h"
#include "htrdr_solve.h"

#include <high_tune/htsky.h>

#include <rsys/algorithm.h>
#include <rsys/clock_time.h>
#include <rsys/cstr.h>
#include <rsys/dynamic_array_u32.h>
#include <rsys/math.h>
#include <rsys/mutex.h>
#include <star/ssp.h>

#include <omp.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>

#define RNG_SEQUENCE_SIZE 10000

#define TILE_MCODE_NULL UINT32_MAX
#define TILE_SIZE 32 /* Definition in X & Y of a tile */
STATIC_ASSERT(IS_POW2(TILE_SIZE), TILE_SIZE_must_be_a_power_of_2);

enum pixel_format {
  PIXEL_LW,
  PIXEL_SW
};

union pixel {
  struct htrdr_pixel_lw lw;
  struct htrdr_pixel_sw sw;
};

/* Tile of row ordered image pixels */
struct tile {
  struct list_node node;
  struct mem_allocator* allocator;
  ref_T ref;

  struct tile_data {
    uint16_t x, y; /* 2D coordinates of the tile in tile space */
    enum pixel_format format;
    /* Simulate the flexible array member of the C99 standard. */
    union pixel pixels[1/*Dummy element*/];
  } data;
};

/* List of tile to compute onto the MPI process. */
struct proc_work {
  struct mutex* mutex;
  struct darray_u32 tiles; /* #tiles to render */
  size_t itile; /* Next tile to render in the above list of tiles */
};

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

static FINLINE struct tile*
tile_create(struct mem_allocator* allocator, const enum pixel_format fmt)
{
  struct tile* tile;
  const size_t tile_sz =
    sizeof(struct tile) - sizeof(union pixel)/*rm dummy pixel*/;
  const size_t buf_sz = /* Flexiblbe array element */
    TILE_SIZE*TILE_SIZE*sizeof(union pixel);
  ASSERT(allocator);

  tile = MEM_ALLOC(allocator, tile_sz+buf_sz);
  if(!tile) return NULL;

  tile->data.format = fmt;

  ref_init(&tile->ref);
  list_init(&tile->node);
  tile->allocator = allocator;
  ASSERT(IS_ALIGNED(&tile->data.pixels, ALIGNOF(union pixel)));

  return tile;
}

static INLINE void
tile_ref_get(struct tile* tile)
{
  ASSERT(tile);
  tile_ref_get(tile);
}

static INLINE void
release_tile(ref_T* ref)
{
  struct tile* tile = CONTAINER_OF(ref, struct tile, ref);
  ASSERT(ref);
  MEM_RM(tile->allocator, tile);
}

static INLINE void
tile_ref_put(struct tile* tile)
{
  ASSERT(tile);
  ref_put(&tile->ref, release_tile);
}

static FINLINE union pixel*
tile_at
  (struct tile* tile,
   const size_t x, /* In tile space */
   const size_t y) /* In tile space */
{
  ASSERT(tile && x < TILE_SIZE && y < TILE_SIZE);
  return tile->data.pixels + (y*TILE_SIZE + x);
}
static void
write_tile_data
  (struct htrdr_buffer* buf,
   const struct tile_data* tile_data)
{
  struct htrdr_buffer_layout layout = HTRDR_BUFFER_LAYOUT_NULL;
  size_t icol, irow;
  size_t irow_tile;
  size_t ncols_tile, nrows_tile;
  char* buf_mem;
  ASSERT(buf && tile_data);

  htrdr_buffer_get_layout(buf, &layout);
  buf_mem = htrdr_buffer_get_data(buf);
  ASSERT(layout.elmt_size == (tile_data->format == PIXEL_LW
    ? sizeof(struct htrdr_pixel_lw)
    : sizeof(struct htrdr_pixel_sw)));

  /* Compute the row/column of the tile origin into the buffer */
  icol = tile_data->x * (size_t)TILE_SIZE;
  irow = tile_data->y * (size_t)TILE_SIZE;

  /* Define the number of tile row/columns to write into the buffer */
  ncols_tile = MMIN(icol + TILE_SIZE, layout.width)  - icol;
  nrows_tile = MMIN(irow + TILE_SIZE, layout.height) - irow;

  /* Copy the row ordered tile data */
  FOR_EACH(irow_tile, 0, nrows_tile) {
    char* buf_row = buf_mem + (irow + irow_tile) * layout.pitch;
    char* buf_col = buf_row + icol * layout.elmt_size;
    const union pixel* tile_row = tile_data->pixels + irow_tile*TILE_SIZE;
    size_t x;

    FOR_EACH(x, 0, ncols_tile) {
      switch(tile_data->format) {
        case PIXEL_LW:
          ((struct htrdr_pixel_lw*)buf_col)[x] = tile_row[x].lw;
          break;
        case PIXEL_SW:
          ((struct htrdr_pixel_sw*)buf_col)[x] = tile_row[x].sw;
          break;
        default: FATAL("Unreachable code.\n"); break;
      }
    }
  }
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
  uint8_t ntiles_to_steal = MMIN((uint8_t)(htrdr->nthreads*2), 16);
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
mpi_gather_tiles
  (struct htrdr* htrdr,
   struct htrdr_buffer* buf,
   const size_t ntiles,
   struct list_node* tiles)
{
  /* Compute the size of the tile_data */
  const size_t msg_sz =
    sizeof(struct tile_data) - sizeof(union pixel)/*dummy*/
  + TILE_SIZE*TILE_SIZE*sizeof(union pixel);

  struct list_node* node = NULL;
  struct tile* tile = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && tiles);
  ASSERT(htrdr->mpi_rank != 0 || buf);
  (void)ntiles;

  if(htrdr->mpi_rank != 0) { /* Non master process */
    /* Send the computed tile to the master process */
    LIST_FOR_EACH(node, tiles) {
      struct tile* t = CONTAINER_OF(node, struct tile, node);
      MPI(Send(&t->data, (int)msg_sz, MPI_CHAR, 0,
        HTRDR_MPI_TILE_DATA, MPI_COMM_WORLD));
    }
  } else { /* Master process */
    size_t itile = 0;

    LIST_FOR_EACH(node, tiles) {
      struct tile* t = CONTAINER_OF(node, struct tile, node);
      write_tile_data(buf, &t->data);
      ++itile;
    }

    if(itile != ntiles) {
      ASSERT(htrdr->mpi_nprocs > 1);

      /* Create a temporary tile to receive the tile data computed by the
       * concurrent MPI processes */
      tile = tile_create
        (htrdr->allocator, htsky_is_long_wave(htrdr->sky) ? PIXEL_LW : PIXEL_SW);
      if(!tile) {
        res = RES_MEM_ERR;
        htrdr_log_err(htrdr,
          "could not allocate the temporary tile used to gather the process "
          "output data -- %s.\n", res_to_cstr(res));
        goto error;
      }

      /* Receive the tile data of the concurrent MPI processes */
      FOR_EACH(itile, itile, ntiles) {
        MPI(Recv(&tile->data, (int)msg_sz, MPI_CHAR, MPI_ANY_SOURCE,
          HTRDR_MPI_TILE_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE));
        write_tile_data(buf, &tile->data);
      }
    }
  }

exit:
  if(tile) tile_ref_put(tile);
  return res;
error:
  goto exit;
}

static void
draw_pixel_sw
  (struct htrdr* htrdr,
   const size_t ithread,
   const size_t ipix[2],
   const double pix_sz[2], /* Size of a pixel in the normalized image plane */
   const struct htrdr_camera* cam,
   const size_t spp,
   struct ssp_rng* rng,
   struct htrdr_pixel_sw* pixel)
{
  struct htrdr_accum XYZ[3]; /* X, Y, and Z */
  struct htrdr_accum time;
  size_t ichannel;
  ASSERT(ipix && ipix && pix_sz && cam && rng && pixel);
  ASSERT(!htsky_is_long_wave(htrdr->sky));

  /* Reset accumulators */
  XYZ[0] = HTRDR_ACCUM_NULL;
  XYZ[1] = HTRDR_ACCUM_NULL;
  XYZ[2] = HTRDR_ACCUM_NULL;
  time = HTRDR_ACCUM_NULL;

  FOR_EACH(ichannel, 0, 3) {
    size_t isamp;

    FOR_EACH(isamp, 0, spp) {
      struct time t0, t1;
      double pix_samp[2];
      double ray_org[3];
      double ray_dir[3];
      double weight;
      double r0, r1, r2;
      double wlen; /* Sampled wavelength into the spectral band */
      double pdf;
      size_t iband; /* Sampled spectral band */
      size_t iquad; /* Sampled quadrature point into the spectral band */
      double usec;

      /* Begin the registration of the time spent to in the realisation */
      time_current(&t0);

      /* Sample a position into the pixel, in the normalized image plane */
      pix_samp[0] = ((double)ipix[0] + ssp_rng_canonical(rng)) * pix_sz[0];
      pix_samp[1] = ((double)ipix[1] + ssp_rng_canonical(rng)) * pix_sz[1];

      /* Generate a ray starting from the pinhole camera and passing through the
       * pixel sample */
      htrdr_camera_ray(cam, pix_samp, ray_org, ray_dir);

      r0 = ssp_rng_canonical(rng);
      r1 = ssp_rng_canonical(rng);
      r2 = ssp_rng_canonical(rng);

      /* Sample a spectral band and a quadrature point */
      switch(ichannel) {
        case 0: wlen = htrdr_cie_xyz_sample_X(htrdr->cie, r0, r1, &pdf); break;
        case 1: wlen = htrdr_cie_xyz_sample_Y(htrdr->cie, r0, r1, &pdf); break;
        case 2: wlen = htrdr_cie_xyz_sample_Z(htrdr->cie, r0, r1, &pdf); break;
        default: FATAL("Unreachable code.\n"); break;
      }

      iband = htsky_find_spectral_band(htrdr->sky, wlen);
      iquad = htsky_spectral_band_sample_quadrature(htrdr->sky, r2, iband);

      /* Compute the radiance in W/m^2/sr/m */
      weight = htrdr_compute_radiance_sw
        (htrdr, ithread, rng, ray_org, ray_dir, wlen, iband, iquad);
      ASSERT(weight >= 0);

      pdf *= 1.e9; /* Transform the pdf fro nm^-1 to m^-1 */
      weight /= pdf; /* In W/m^2/sr */

      /* End the registration of the per realisation time */
      time_sub(&t0, time_current(&t1), &t0);
      usec = (double)time_val(&t0, TIME_NSEC) * 0.001;

      /* Update the pixel accumulator of the current channel */
      XYZ[ichannel].sum_weights += weight;
      XYZ[ichannel].sum_weights_sqr += weight*weight;
      XYZ[ichannel].nweights += 1;

      /* Update the pixel accumulator of per realisation time */
      time.sum_weights += usec;
      time.sum_weights_sqr += usec*usec;
      time.nweights += 1;
    }
  }

  /* Flush pixel data */
  htrdr_accum_get_estimation(XYZ+0, &pixel->X);
  htrdr_accum_get_estimation(XYZ+1, &pixel->Y);
  htrdr_accum_get_estimation(XYZ+2, &pixel->Z);
  pixel->time = time;
}

static void
draw_pixel_lw
  (struct htrdr* htrdr,
   const size_t ithread,
   const size_t ipix[2],
   const double pix_sz[2], /* Size of a pixel in the normalized image plane */
   const struct htrdr_camera* cam,
   const size_t spp,
   struct ssp_rng* rng,
   struct htrdr_pixel_lw* pixel)
{
  struct htrdr_accum radiance;
  struct htrdr_accum time;
  size_t isamp;
  double temp_min, temp_max;
  ASSERT(ipix && ipix && pix_sz && cam && rng && pixel);
  ASSERT(htsky_is_long_wave(htrdr->sky));

  /* Reset the pixel accumulators */
  radiance = HTRDR_ACCUM_NULL;
  time = HTRDR_ACCUM_NULL;

  FOR_EACH(isamp, 0, spp) {
    struct time t0, t1;
    double pix_samp[2];
    double ray_org[3];
    double ray_dir[3];
    double weight;
    double r0, r1, r2;
    double wlen;
    size_t iband;
    size_t iquad;
    double usec;
    double band_pdf;

    /* Begin the registration of the time spent in the realisation */
    time_current(&t0);

    /* Sample a position into the pixel, in the normalized image plane */
    pix_samp[0] = ((double)ipix[0] + ssp_rng_canonical(rng)) * pix_sz[0];
    pix_samp[1] = ((double)ipix[1] + ssp_rng_canonical(rng)) * pix_sz[1];

    /* Generate a ray starting from the pinhole camera and passing through the
     * pixel sample */
    htrdr_camera_ray(cam, pix_samp, ray_org, ray_dir);

    r0 = ssp_rng_canonical(rng);
    r1 = ssp_rng_canonical(rng);
    r2 = ssp_rng_canonical(rng);

    /* Sample a wavelength */
    wlen = htrdr_ran_lw_sample(htrdr->ran_lw, r0, r1, &band_pdf);

    /* Select the associated band and sample a quadrature point */
    iband = htsky_find_spectral_band(htrdr->sky, wlen);
    iquad = htsky_spectral_band_sample_quadrature(htrdr->sky, r2, iband);

    /* Compute the integrated luminance in W/m^2/sr */
    weight = htrdr_compute_radiance_lw(htrdr, ithread, rng, ray_org, ray_dir,
      wlen, iband, iquad);

    /* Importance sampling: correct weight with pdf */
    weight /= band_pdf;

    /* From integrated radiance to average radiance in W/m^2/sr/m */
    if(htrdr->wlen_range_m[0] != htrdr->wlen_range_m[1]) { 
      /* Is not monochromatic */
      weight /= (htrdr->wlen_range_m[1] - htrdr->wlen_range_m[0]) ;
    }
    ASSERT(weight >= 0);

    /* End the registration of the per realisation time */
    time_sub(&t0, time_current(&t1), &t0);
    usec = (double)time_val(&t0, TIME_NSEC) * 0.001;

    /* Update the pixel accumulator of the current channel */
    radiance.sum_weights += weight;
    radiance.sum_weights_sqr += weight*weight;
    radiance.nweights += 1;

    /* Update the pixel accumulator of per realisation time */
    time.sum_weights += usec;
    time.sum_weights_sqr += usec*usec;
    time.nweights += 1;
  }

  /* Compute the estimation of the pixel radiance */
  htrdr_accum_get_estimation(&radiance, &pixel->radiance);

  /* Save the per realisation integration time */
  pixel->time = time;

  /* Compute the brightness_temperature of the pixel and estimate its standard
   * error */
  #define BRIGHTNESS_TEMPERATURE(Radiance, Temperature) {                      \
    res_T res = brightness_temperature                                         \
      (htrdr,                                                                  \
       htrdr->wlen_range_m[0],                                                 \
       htrdr->wlen_range_m[1],                                                 \
       (Radiance),                                                             \
       &(Temperature));                                                        \
    if(res != RES_OK) {                                                        \
      htrdr_log_warn(htrdr,                                                    \
        "Could not compute the brightness temperature for the radiance %g.\n", \
        (Radiance));                                                           \
      (Temperature) = 0;                                                       \
    }                                                                          \
  } (void)0
  BRIGHTNESS_TEMPERATURE(pixel->radiance.E, pixel->radiance_temperature.E);
  BRIGHTNESS_TEMPERATURE(pixel->radiance.E - pixel->radiance.SE, temp_min);
  BRIGHTNESS_TEMPERATURE(pixel->radiance.E + pixel->radiance.SE, temp_max);
  pixel->radiance_temperature.SE = temp_max - temp_min;
  #undef BRIGHTNESS_TEMPERATURE

  /* Transform the pixel radiance from W/m^2/sr/m to W/m^/sr/nm */
  pixel->radiance.E *= 1.0e-9;
  pixel->radiance.SE *= 1.0e-9;
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
   struct tile* tile)
{
  size_t npixels;
  size_t mcode; /* Morton code of tile pixel */
  ASSERT(htrdr && tile_org && tile_sz && pix_sz && cam && spp && tile);
  (void)tile_mcode;
  /* Adjust the #pixels to process them wrt a morton order */
  npixels = round_up_pow2(MMAX(tile_sz[0], tile_sz[1]));
  npixels *= npixels;

  FOR_EACH(mcode, 0, npixels) {
    union pixel* pixel;
    size_t ipix_tile[2]; /* Pixel coord in the tile */
    size_t ipix[2]; /* Pixel coord in the buffer */

    ipix_tile[0] = morton2D_decode((uint32_t)(mcode>>0));
    if(ipix_tile[0] >= tile_sz[0]) continue; /* Pixel is out of tile */
    ipix_tile[1] = morton2D_decode((uint32_t)(mcode>>1));
    if(ipix_tile[1] >= tile_sz[1]) continue; /* Pixel is out of tile */

    /* Fetch and reset the pixel accumulator */
    pixel = tile_at(tile, ipix_tile[0], ipix_tile[1]);

    /* Compute the pixel coordinate */
    ipix[0] = tile_org[0] + ipix_tile[0];
    ipix[1] = tile_org[1] + ipix_tile[1];

    /* Draw the pixel */
    if(htsky_is_long_wave(htrdr->sky)) {
      draw_pixel_lw(htrdr, ithread, ipix, pix_sz, cam, spp, rng, &pixel->lw);
    } else {
      draw_pixel_sw(htrdr, ithread, ipix, pix_sz, cam, spp, rng, &pixel->sw);
    }
  }
  return RES_OK;
}

static res_T
draw_image
  (struct htrdr* htrdr,
   const struct htrdr_camera* cam,
   const size_t width, /* Image width */
   const size_t height, /* Image height */
   const size_t spp,
   const size_t ntiles_x,
   const size_t ntiles_y,
   const size_t ntiles_adjusted,
   const double pix_sz[2], /* Pixel size in the normalized image plane */
   struct proc_work* work,
   struct list_node* tiles)
{
  struct ssp_rng* rng_proc = NULL;
  size_t nthreads = 0;
  size_t nthieves = 0;
  size_t proc_ntiles = 0;
  enum pixel_format pixfmt;
  ATOMIC nsolved_tiles = 0;
  ATOMIC res = RES_OK;
  ASSERT(htrdr && cam && spp && ntiles_adjusted && work && tiles);
  ASSERT(pix_sz && pix_sz[0] > 0 && pix_sz[1] > 0);
  ASSERT(width && height);
  (void)ntiles_x, (void)ntiles_y;

  res = ssp_rng_create(htrdr->allocator, &ssp_rng_mt19937_64, &rng_proc);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "could not create the RNG used to sample a process "
      "to steal -- %s.\n", res_to_cstr((res_T)res));
    goto error;
  }

  proc_ntiles = proc_work_get_ntiles(work);
  nthreads = MMIN(htrdr->nthreads, proc_ntiles);

  /* The process is not considered as a working process for himself */
  htrdr->mpi_working_procs[htrdr->mpi_rank] = 0;
  --htrdr->mpi_nworking_procs;

  pixfmt = htsky_is_long_wave(htrdr->sky) ? PIXEL_LW : PIXEL_SW;

  omp_set_num_threads((int)nthreads);
  #pragma omp parallel
  for(;;) {
    const int ithread = omp_get_thread_num();
    struct ssp_rng_proxy* rng_proxy = NULL;
    struct ssp_rng* rng;
    struct tile* tile;
    uint32_t mcode = TILE_MCODE_NULL;
    size_t tile_org[2];
    size_t tile_sz[2];
    size_t n;
    res_T res_local = RES_OK;
    int32_t pcent;

    /* Get a tile to draw */
    #pragma omp critical
    {
      mcode = proc_work_get_tile(work);
      if(mcode == TILE_MCODE_NULL) { /* No more work on this process */
        /* Try to steal works to concurrent processes */
        proc_work_reset(work);
        nthieves = mpi_steal_work(htrdr, rng_proc, work);
        if(nthieves != 0) {
          mcode = proc_work_get_tile(work);
        }
      }
    }
    if(mcode == TILE_MCODE_NULL) break; /* No more work */

    /* Decode the morton code to retrieve the tile index  */
    tile_org[0] = morton2D_decode((uint32_t)(mcode>>0));
    tile_org[1] = morton2D_decode((uint32_t)(mcode>>1));
    ASSERT(tile_org[0] < ntiles_x && tile_org[1] < ntiles_y);

    /* Create the tile */
    tile = tile_create(htrdr->allocator, pixfmt);
    if(!tile) {
      ATOMIC_SET(&res, RES_MEM_ERR);
      htrdr_log_err(htrdr,
        "could not allocate the memory space of the tile (%lu, %lu) -- %s.\n",
        (unsigned long)tile_org[0], (unsigned long)tile_org[1],
         res_to_cstr((res_T)ATOMIC_GET(&res)));
      break;
    }

    /* Register the tile */
    #pragma omp critical
    list_add_tail(tiles, &tile->node);

    tile->data.x = (uint16_t)tile_org[0];
    tile->data.y = (uint16_t)tile_org[1];

    /* Define the tile origin in pixel space */
    tile_org[0] *= TILE_SIZE;
    tile_org[1] *= TILE_SIZE;

    /* Compute the size of the tile clamped by the borders of the buffer */
    tile_sz[0] = MMIN(TILE_SIZE, width - tile_org[0]);
    tile_sz[1] = MMIN(TILE_SIZE, height - tile_org[1]);

    /* Create a proxy RNG for the current tile. This proxy is used for the
     * current thread only and thus it has to manage only one RNG. This proxy
     * is initialised in order to ensure that an unique and predictable set of
     * random numbers is used for the current tile. */
    SSP(rng_proxy_create2
      (&htrdr->lifo_allocators[ithread],
       &ssp_rng_threefry,
       RNG_SEQUENCE_SIZE * (size_t)mcode, /* Offset */
       RNG_SEQUENCE_SIZE, /* Size */
       RNG_SEQUENCE_SIZE * (size_t)ntiles_adjusted, /* Pitch */
       1, &rng_proxy));
    SSP(rng_proxy_create_rng(rng_proxy, 0, &rng));

    /* Launch the tile rendering */
    res_local = draw_tile(htrdr, (size_t)ithread, mcode, tile_org, tile_sz,
      pix_sz, cam, spp, rng, tile);

    SSP(rng_proxy_ref_put(rng_proxy));
    SSP(rng_ref_put(rng));

    if(res_local != RES_OK) {
      ATOMIC_SET(&res, res_local);
      break;
    }

    /* Update the progress status */
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
  goto exit;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_draw_radiance
  (struct htrdr* htrdr,
   const struct htrdr_camera* cam,
   const size_t width,
   const size_t height,
   const size_t spp,
   struct htrdr_buffer* buf)
{
  char strbuf[128];
  struct time t0, t1;
  struct list_node tiles;
  size_t ntiles_x, ntiles_y, ntiles, ntiles_adjusted;
  size_t itile;
  struct proc_work work;
  struct htrdr_buffer_layout layout = HTRDR_BUFFER_LAYOUT_NULL;
  size_t proc_ntiles_adjusted;
  double pix_sz[2];
  ATOMIC probe_thieves = 1;
  ATOMIC res = RES_OK;
  ASSERT(htrdr && cam && width && height);
  ASSERT(htrdr->mpi_rank != 0 || buf);

  list_init(&tiles);
  proc_work_init(htrdr->allocator, &work);

  if(htrdr->mpi_rank == 0) {
    size_t pixsz, pixal;
    htrdr_buffer_get_layout(buf, &layout);
    ASSERT(layout.width || layout.height || layout.elmt_size);
    ASSERT(layout.width == width && layout.height == height);

    if(htsky_is_long_wave(htrdr->sky)) {
      pixsz = sizeof(struct htrdr_pixel_lw);
      pixal = ALIGNOF(struct htrdr_pixel_lw);
    } else {
      pixsz = sizeof(struct htrdr_pixel_sw);
      pixal = ALIGNOF(struct htrdr_pixel_sw);
    }

    if(layout.elmt_size != pixsz || layout.alignment < pixal) {
      htrdr_log_err(htrdr, "%s: invalid buffer layout.\n", FUNC_NAME);
      res = RES_BAD_ARG;
      goto error;
    }
  }

  /* Compute the overall number of tiles */
  ntiles_x = (width + (TILE_SIZE-1)/*ceil*/)/TILE_SIZE;
  ntiles_y = (height+ (TILE_SIZE-1)/*ceil*/)/TILE_SIZE;
  ntiles = ntiles_x * ntiles_y;

  /* Compute the pixel size in the normalized image plane */
  pix_sz[0] = 1.0 / (double)width;
  pix_sz[1] = 1.0 / (double)height;

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

  time_current(&t0);

  omp_set_nested(1); /* Enable nested threads for draw_image */
  #pragma omp parallel sections num_threads(2)
  {
    #pragma omp section
    mpi_probe_thieves(htrdr, &work, &probe_thieves);

    #pragma omp section
    {
      draw_image(htrdr, cam, width, height, spp, ntiles_x, ntiles_y,
        ntiles_adjusted, pix_sz, &work, &tiles);
      /* The processes have no more work to do. Stop probing for thieves */
      ATOMIC_SET(&probe_thieves, 0);
    }
  }

  if(htrdr->mpi_rank == 0) {
    update_mpi_progress(htrdr, HTRDR_MPI_PROGRESS_RENDERING);
    fprintf(stderr, "\n"); /* Add a new line after the progress statuses */
  }

  time_sub(&t0, time_current(&t1), &t0);
  time_dump(&t0, TIME_ALL, NULL, strbuf, sizeof(strbuf));
  htrdr_log(htrdr, "Rendering time: %s\n", strbuf);

  /* Gather accum buffers from the group of processes */
  time_current(&t0);
  res = mpi_gather_tiles(htrdr, buf, ntiles, &tiles);
  if(res != RES_OK) goto error;
  time_sub(&t0, time_current(&t1), &t0);
  time_dump(&t0, TIME_ALL, NULL, strbuf, sizeof(strbuf));
  htrdr_log(htrdr, "Image gathering time: %s\n", strbuf);

exit:
  { /* Free allocated tiles */
    struct list_node* node;
    struct list_node* tmp;
    LIST_FOR_EACH_SAFE(node, tmp, &tiles) {
      struct tile* tile = CONTAINER_OF(node, struct tile, node);
      list_del(node);
      tile_ref_put(tile);
    }
  }
  proc_work_release(&work);
  return (res_T)res;
error:
  goto exit;
}

