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
#include <rsys/math.h>
#include <star/ssp.h>

#include <omp.h>
#include <mpi.h>
#include <unistd.h>

#define RNG_SEQUENCE_SIZE 1000000

#define TILE_SIZE 32 /* Definition in X & Y of a tile */
STATIC_ASSERT(IS_POW2(TILE_SIZE), TILE_SIZE_must_be_a_power_of_2);

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

static res_T
gather_buffer
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
    err = MPI_Gather(buf_row_accums, (int)layout.pitch, MPI_CHAR, gathered_accums,
      (int)layout.pitch, MPI_CHAR, 0, MPI_COMM_WORLD);
    if(err != MPI_SUCCESS) {
      htrdr_log_err(htrdr,
"could not gather the buffer line `%lu' from the group of processes -- %s.\n",
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

/* Return the overall percentage */
static size_t
fetch_process_progress(struct htrdr* htrdr, size_t* progress)
{
  size_t overall_pcent;
  int iproc;
  ASSERT(htrdr && progress && htrdr->mpi_rank == 0);

  overall_pcent = progress[0];
  FOR_EACH(iproc, 1, htrdr->mpi_nprocs) {
    size_t proc_pcent = progress[iproc];

    /* Flush the last sent percentage of the process `iproc' */
    for(;;) {
      int flag;

      CHK(MPI_Iprobe(iproc, 0/*tag*/, MPI_COMM_WORLD, &flag,
        MPI_STATUS_IGNORE) == MPI_SUCCESS);
      if(flag == 0) break; /* No more message */

      CHK(MPI_Recv(&proc_pcent, sizeof(size_t), MPI_CHAR, iproc, 0/*tag*/,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE) == MPI_SUCCESS);
    }

    progress[iproc] = proc_pcent;
    overall_pcent += progress[iproc];
  }
  return overall_pcent / (size_t)htrdr->mpi_nprocs;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_draw_radiance_sw
  (struct htrdr* htrdr,
   const struct htrdr_camera* cam,
   const size_t total_spp,
   struct htrdr_buffer* buf)
{
  struct ssp_rng_proxy* rng_proxy = NULL;
  struct ssp_rng** rngs = NULL;
  size_t ntiles_x, ntiles_y, ntiles, ntiles_adjusted;
  size_t* progress = NULL;
  size_t i;
  int32_t mcode; /* Morton code of the tile */
  struct htrdr_buffer_layout layout;
  double pix_sz[2]; /* Pixel size in the normalized image plane */
  size_t spp;
  size_t overall_pcent = 0;
  ATOMIC nsolved_tiles = 0;
  ATOMIC res = RES_OK;
  ASSERT(htrdr && cam && buf);

  htrdr_buffer_get_layout(buf, &layout);
  ASSERT(layout.width || layout.height || layout.elmt_size);

  spp = total_spp / (size_t)htrdr->mpi_nprocs;

  /* Add the remaining realisations to the 1st process */
  if(htrdr->mpi_rank == 0) {
    spp += total_spp - (spp*(size_t)htrdr->mpi_nprocs);
  }

  if(layout.elmt_size != sizeof(struct htrdr_accum[3])/*#channels*/
  || layout.alignment < ALIGNOF(struct htrdr_accum[3])) {
    htrdr_log_err(htrdr,
      "%s: invalid buffer layout. "
      "The pixel size must be the size of 3 * accumulators.\n",
      FUNC_NAME);
    res = RES_BAD_ARG;
    goto error;
  }

  res = ssp_rng_proxy_create2
    (htrdr->allocator,
     &ssp_rng_mt19937_64,
     RNG_SEQUENCE_SIZE * (size_t)htrdr->mpi_rank, /* Offset */
     RNG_SEQUENCE_SIZE, /* Size */
     RNG_SEQUENCE_SIZE * (size_t)htrdr->mpi_nprocs, /* Pitch */
     htrdr->nthreads, &rng_proxy);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "%s: could not create the RNG proxy -- %s.\n",
      FUNC_NAME, res_to_cstr((res_T)res));
    goto error;
  }

  rngs = MEM_CALLOC(htrdr->allocator, htrdr->nthreads, sizeof(*rngs));
  if(!rngs) {
    htrdr_log_err(htrdr, "%s: could not allocate the RNGs list.\n", FUNC_NAME);
    res = RES_MEM_ERR;
    goto error;
  }

  FOR_EACH(i, 0, htrdr->nthreads) {
    res = ssp_rng_proxy_create_rng(rng_proxy, i, &rngs[i]);
    if(res != RES_OK) {
      htrdr_log_err(htrdr,"%s: could not create the per thread RNG -- %s.\n",
        FUNC_NAME, res_to_cstr((res_T)res));
      goto error;
    }
  }

  if(htrdr->mpi_rank == 0) {
    progress = MEM_CALLOC
      (htrdr->allocator, (size_t)htrdr->mpi_nprocs, sizeof(size_t));
  } else {
    progress = MEM_CALLOC
      (htrdr->allocator, 1, sizeof(size_t));
  }
  if(!progress) {
    htrdr_log_err(htrdr,
      "%s: could not allocate the process progress counter.\n", FUNC_NAME);
    res = RES_MEM_ERR;
    goto error;
  }

  ntiles_x = (layout.width + (TILE_SIZE-1)/*ceil*/)/TILE_SIZE;
  ntiles_y = (layout.height+ (TILE_SIZE-1)/*ceil*/)/TILE_SIZE;
  ntiles_adjusted = round_up_pow2(MMAX(ntiles_x, ntiles_y));
  ntiles_adjusted *= ntiles_adjusted;
  ntiles = ntiles_x * ntiles_y;

  pix_sz[0] = 1.0 / (double)layout.width;
  pix_sz[1] = 1.0 / (double)layout.height;

  htrdr_fprintf(htrdr, stderr, "Rendering: %3i%%", 0);
  htrdr_fflush(htrdr, stderr);

  omp_set_num_threads((int)htrdr->nthreads);
  #pragma omp parallel for schedule(static, 1/*chunck size*/)
  for(mcode=0; mcode<(int64_t)ntiles_adjusted; ++mcode) {
    const int ithread = omp_get_thread_num();
    struct ssp_rng* rng = rngs[ithread];
    size_t tile_org[2];
    size_t tile_sz[2];
    size_t pcent;
    size_t n;
    res_T res_local = RES_OK;

    /* Decode the morton code to retrieve the tile index  */
    tile_org[0] = morton2D_decode((uint32_t)(mcode>>0));
    if(tile_org[0] >= ntiles_x) continue; /* Skip border tile */
    tile_org[1] = morton2D_decode((uint32_t)(mcode>>1));
    if(tile_org[1] >= ntiles_y) continue; /* Skip border tile */

    /* Define the tile origin in pixel space */
    tile_org[0] *= TILE_SIZE;
    tile_org[1] *= TILE_SIZE;

    /* Compute the size of the tile clamped by the borders of the buffer */
    tile_sz[0] = MMIN(TILE_SIZE, layout.width - tile_org[0]);
    tile_sz[1] = MMIN(TILE_SIZE, layout.height - tile_org[1]);

    res_local = draw_tile(htrdr, (size_t)ithread, mcode, tile_org, tile_sz,
      pix_sz, cam, spp, rng, buf);
    if(res_local != RES_OK) {
      ATOMIC_SET(&res, res_local);
      continue;
    }

    n = (size_t)ATOMIC_INCR(&nsolved_tiles);
    pcent = n * 100 / ntiles;

    #pragma omp critical
    if((size_t)pcent > progress[0]) {
      progress[0] = pcent;

      if(htrdr->mpi_rank != 0) {
        /* Send the progress percentage of the process to the master process */
        CHK(MPI_Send(&pcent, sizeof(size_t), MPI_CHAR, 0/*dst*/, 0/*tag*/,
          MPI_COMM_WORLD) == MPI_SUCCESS);
      } else {
        overall_pcent = fetch_process_progress(htrdr, progress);
        htrdr_fprintf(htrdr, stderr, "%c[2K\rRendering: %3lu%%",
          27, (unsigned long)overall_pcent);
        htrdr_fflush(htrdr, stderr);
      }
    }

    if(ATOMIC_GET(&res) != RES_OK) continue;
  }

  if(htrdr->mpi_rank == 0) {
    while(overall_pcent != 100) {
      overall_pcent = fetch_process_progress(htrdr, progress);
      htrdr_fprintf(htrdr, stderr, "%c[2K\rRendering: %3lu%%",
        27, (unsigned long)overall_pcent);
      htrdr_fflush(htrdr, stderr);
      sleep(1);
    }
    htrdr_fprintf(htrdr, stderr, "\n");
  }

  /* Gather accum buffers from the group of processes */
  res = gather_buffer(htrdr, buf);
  if(res != RES_OK) goto error;

exit:
  if(rng_proxy) SSP(rng_proxy_ref_put(rng_proxy));
  if(progress) MEM_RM(htrdr->allocator, progress);
  if(rngs) {
    FOR_EACH(i, 0, htrdr->nthreads) {
      if(rngs[i]) SSP(rng_ref_put(rngs[i]));
    }
    MEM_RM(htrdr->allocator, rngs);
  }
  return (res_T)res;
error:
  goto exit;
}

