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

#define _POSIX_C_SOURCE 200809L /* stat.st_time support */

#include "htrdr.h"
#include "htrdr_c.h"
#include "htrdr_args.h"
#include "htrdr_buffer.h"
#include "htrdr_cie_xyz.h"
#include "htrdr_camera.h"
#include "htrdr_ground.h"
#include "htrdr_mtl.h"
#include "htrdr_ran_wlen.h"
#include "htrdr_rectangle.h"
#include "htrdr_sun.h"
#include "htrdr_solve.h"

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
static void
print_out(const char* msg, void* ctx)
{
  ASSERT(msg);
  (void)ctx;
#ifdef OS_UNIX
  fprintf(stderr, "\x1b[1m\x1b[32m>\x1b[0m %s", msg);
#else
  fprintf(stderr, "> %s", msg);
#endif
}

static void
print_err(const char* msg, void* ctx)
{
  ASSERT(msg);
  (void)ctx;
#ifdef OS_UNIX
  fprintf(stderr, "\x1b[31merror:\x1b[0m %s", msg);
#else
  fprintf(stderr, "error: %s", msg);
#endif
}

static void
print_warn(const char* msg, void* ctx)
{
  ASSERT(msg);
  (void)ctx;
#ifdef OS_UNIX
  fprintf(stderr, "\x1b[33mwarning:\x1b[0m %s", msg);
#else
  fprintf(stderr,"warning: %s", msg);
#endif
}

static void
log_msg
  (struct htrdr* htrdr,
   const enum log_type stream,
   const char* msg,
   va_list vargs)
{
  ASSERT(htrdr && msg);
  if(htrdr->verbose) {
    CHK(logger_vprint(&htrdr->logger, stream, msg, vargs) == RES_OK);
  }
}

static enum htsky_spectral_type
htrdr_to_sky_spectral_type(const enum htrdr_spectral_type type)
{
  enum htsky_spectral_type spectype;
  switch(type) {
    case HTRDR_SPECTRAL_LW:
      spectype = HTSKY_SPECTRAL_LW;
      break;
    case HTRDR_SPECTRAL_SW:
    case HTRDR_SPECTRAL_SW_CIE_XYZ:
      spectype = HTSKY_SPECTRAL_SW;
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
  return spectype;
}

static res_T
dump_buffer
  (struct htrdr* htrdr,
   struct htrdr_buffer* buf,
   struct htrdr_accum* time_acc, /* May be NULL */
   struct htrdr_accum* flux_acc, /* May be NULL */
   const char* stream_name,
   FILE* stream)
{
  struct htrdr_buffer_layout layout;
  size_t pixsz, pixal;
  size_t x, y;
  res_T res = RES_OK;
  ASSERT(htrdr && buf && stream_name && stream);
  (void)stream_name;

  pixsz = htrdr_spectral_type_get_pixsz
    (htrdr->spectral_type, htrdr->sensor.type);
  pixal = htrdr_spectral_type_get_pixal
    (htrdr->spectral_type, htrdr->sensor.type);

  htrdr_buffer_get_layout(buf, &layout);
  if(layout.elmt_size != pixsz || layout.alignment != pixal) {
    htrdr_log_err(htrdr, "%s: invalid buffer layout. ", FUNC_NAME);
    res = RES_BAD_ARG;
    goto error;
  }

  fprintf(stream, "%lu %lu\n", layout.width, layout.height);

  if(time_acc) *time_acc = HTRDR_ACCUM_NULL;
  if(flux_acc) *flux_acc = HTRDR_ACCUM_NULL;
  FOR_EACH(y, 0, layout.height) {
    FOR_EACH(x, 0, layout.width) {
      struct htrdr_estimate pix_time = HTRDR_ESTIMATE_NULL;
      const struct htrdr_accum* pix_time_acc = NULL;

      if(htrdr->sensor.type == HTRDR_SENSOR_RECTANGLE) {
        const struct htrdr_pixel_flux* pix = htrdr_buffer_at(buf, x, y);
        struct htrdr_estimate flux = HTRDR_ESTIMATE_NULL;

        if(pix->flux.nweights == 0) {
          fprintf(stream, "0 0 0 0 0 0 ");
        } else {
          htrdr_accum_get_estimation(&pix->flux, &flux);
          fprintf(stream, "%g %g 0 0 0 0 ", flux.E, flux.SE);

          if(flux_acc) {
            flux_acc->sum_weights += pix->flux.sum_weights;
            flux_acc->sum_weights_sqr += pix->flux.sum_weights_sqr;
            flux_acc->nweights += pix->flux.nweights;
          }
        }
        pix_time_acc = &pix->time;

      } else {
        if(htrdr->spectral_type != HTRDR_SPECTRAL_SW_CIE_XYZ){
          const struct htrdr_pixel_xwave* pix = htrdr_buffer_at(buf, x, y);
          fprintf(stream, "%g %g ",
            pix->radiance_temperature.E, pix->radiance_temperature.SE);
          fprintf(stream, "%g %g ", pix->radiance.E, pix->radiance.SE);
          fprintf(stream, "0 0 ");
          pix_time_acc = &pix->time;

        } else {
          const struct htrdr_pixel_image* pix = htrdr_buffer_at(buf, x, y);
          fprintf(stream, "%g %g ", pix->X.E, pix->X.SE);
          fprintf(stream, "%g %g ", pix->Y.E, pix->Y.SE);
          fprintf(stream, "%g %g ", pix->Z.E, pix->Z.SE);
          pix_time_acc = &pix->time;
        }
      }

      htrdr_accum_get_estimation(pix_time_acc, &pix_time);
      fprintf(stream, "%g %g\n", pix_time.E, pix_time.SE);

      if(time_acc) {
        time_acc->sum_weights += pix_time_acc->sum_weights;
        time_acc->sum_weights_sqr += pix_time_acc->sum_weights_sqr;
        time_acc->nweights += pix_time_acc->nweights;
      }
    }
    fprintf(stream, "\n");
  }

exit:
  return res;
error:
  goto exit;
}

static INLINE void
spherical_to_cartesian_dir
  (const double azimuth, /* In radians */
   const double elevation, /* In radians */
   double dir[3])
{
  double cos_azimuth;
  double sin_azimuth;
  double cos_elevation;
  double sin_elevation;
  ASSERT(azimuth >= 0 && azimuth < 2*PI);
  ASSERT(elevation >= 0 && elevation <= PI/2.0);
  ASSERT(dir);

  cos_azimuth = cos(azimuth);
  sin_azimuth = sin(azimuth);
  cos_elevation = cos(elevation);
  sin_elevation = sin(elevation);

  dir[0] = cos_elevation * cos_azimuth;
  dir[1] = cos_elevation * sin_azimuth;
  dir[2] = sin_elevation;
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

  if(htrdr->mpi_nprocs != 1)
    mpi_print_proc_info(htrdr);

exit:
  return res;
error:
  release_mpi(htrdr);
  goto exit;
}

static res_T
setup_sensor(struct htrdr* htrdr, const struct htrdr_args* args)
{
  double proj_ratio;
  res_T res = RES_OK;
  ASSERT(htrdr && args);

  htrdr->sensor.type = args->sensor_type;

  if(args->spectral_type == HTRDR_SPECTRAL_SW_CIE_XYZ
  && args->sensor_type != HTRDR_SENSOR_CAMERA) {
    htrdr_log_err(htrdr, "the CIE 1931 XYZ spectral integration can be used "
      "only with a camera sensor.\n");
    res = RES_BAD_ARG;
    goto error;
  }

  switch(args->sensor_type) {
    case HTRDR_SENSOR_CAMERA:
      proj_ratio =
        (double)args->image.definition[0]
      / (double)args->image.definition[1];
      res = htrdr_camera_create(htrdr, args->camera.pos, args->camera.tgt,
        args->camera.up, proj_ratio, MDEG2RAD(args->camera.fov_y),
        &htrdr->sensor.camera);
      break;
    case HTRDR_SENSOR_RECTANGLE:
      res = htrdr_rectangle_create(htrdr, args->rectangle.sz,
        args->rectangle.pos, args->rectangle.tgt, args->rectangle.up,
        &htrdr->sensor.rectangle);
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  goto exit;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_init
  (struct mem_allocator* mem_allocator,
   const struct htrdr_args* args,
   struct htrdr* htrdr)
{
  struct htsky_args htsky_args = HTSKY_ARGS_DEFAULT;
  double sun_dir[3];
  double spectral_range[2];
  const char* output_name = NULL;
  size_t ithread;
  int nthreads_max;
  res_T res = RES_OK;
  ASSERT(args && htrdr);

  memset(htrdr, 0, sizeof(*htrdr));

  htrdr->allocator = mem_allocator ? mem_allocator : &mem_default_allocator;

  logger_init(htrdr->allocator, &htrdr->logger);
  logger_set_stream(&htrdr->logger, LOG_OUTPUT, print_out, NULL);
  logger_set_stream(&htrdr->logger, LOG_ERROR, print_err, NULL);
  logger_set_stream(&htrdr->logger, LOG_WARNING, print_warn, NULL);
  str_init(htrdr->allocator, &htrdr->output_name);
  nthreads_max = MMAX(omp_get_max_threads(), omp_get_num_procs());
  htrdr->dump_vtk = args->dump_vtk;
  htrdr->verbose = args->verbose;
  htrdr->nthreads = MMIN(args->nthreads, (unsigned)nthreads_max);
  htrdr->spp = args->image.spp;
  htrdr->width = args->image.definition[0];
  htrdr->height = args->image.definition[1];
  htrdr->grid_max_definition[0] = args->grid_max_definition[0];
  htrdr->grid_max_definition[1] = args->grid_max_definition[1];
  htrdr->grid_max_definition[2] = args->grid_max_definition[2];
  htrdr->spectral_type = args->spectral_type;
  htrdr->ref_temperature = args->ref_temperature;

  res = init_mpi(htrdr);
  if(res != RES_OK) goto error;

  if(!args->output) {
    htrdr->output = stdout;
    output_name = "<stdout>";
  } else if(htrdr->mpi_rank != 0) {
    htrdr->output = NULL;
    output_name = "<null>";
  } else {
    res = open_output_stream
      (htrdr, args->output, 0/*read*/, args->force_overwriting, &htrdr->output);
    if(res != RES_OK) goto error;
    output_name = args->output;
  }
  res = str_set(&htrdr->output_name, output_name);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: could not store the name of the output stream `%s' -- %s.\n",
      FUNC_NAME, output_name, res_to_cstr(res));
    goto error;
  }

  /* Disable the Star-3D verbosity since the Embree backend prints some messages
   * on stdout rather than stderr. This is annoying since stdout may be used by
   * htrdr to write output data */
  res = s3d_device_create
    (&htrdr->logger, htrdr->allocator, 0, &htrdr->s3d);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: could not create the Star-3D device -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }

  /* Materials are necessary only if a ground geometry is defined */
  if(args->filename_obj) {
    res = htrdr_mtl_create(htrdr, args->filename_mtl, &htrdr->mtl);
    if(res != RES_OK) goto error;
  }

  res = htrdr_ground_create(htrdr, args->filename_obj, args->repeat_ground,
    &htrdr->ground);
  if(res != RES_OK) goto error;

  res = setup_sensor(htrdr, args);
  if(res != RES_OK) goto error;

  res = htrdr_sun_create(htrdr, &htrdr->sun);
  if(res != RES_OK) goto error;
  spherical_to_cartesian_dir
    (MDEG2RAD(args->sun_azimuth), MDEG2RAD(args->sun_elevation), sun_dir);
  htrdr_sun_set_direction(htrdr->sun, sun_dir);

  htsky_args.htcp_filename = args->filename_les;
  htsky_args.htgop_filename = args->filename_gas;
  htsky_args.htmie_filename = args->filename_mie;
  htsky_args.cache_filename = args->cache;
  htsky_args.grid_max_definition[0] = args->grid_max_definition[0];
  htsky_args.grid_max_definition[1] = args->grid_max_definition[1];
  htsky_args.grid_max_definition[2] = args->grid_max_definition[2];
  htsky_args.optical_thickness = args->optical_thickness;
  htsky_args.nthreads = htrdr->nthreads;
  htsky_args.repeat_clouds = args->repeat_clouds;
  htsky_args.verbose = htrdr->mpi_rank == 0 ? args->verbose : 0;
  htsky_args.spectral_type = htrdr_to_sky_spectral_type(args->spectral_type);
  htsky_args.wlen_range[0] = args->wlen_range[0];
  htsky_args.wlen_range[1] = args->wlen_range[1];
  res = htsky_create(&htrdr->logger, htrdr->allocator, &htsky_args, &htrdr->sky);
  if(res != RES_OK) goto error;

  HTSKY(get_raw_spectral_bounds(htrdr->sky, spectral_range));

  spectral_range[0] = MMAX(args->wlen_range[0], spectral_range[0]);
  spectral_range[1] = MMIN(args->wlen_range[1], spectral_range[1]);
  if(spectral_range[0] != args->wlen_range[0]
  || spectral_range[1] != args->wlen_range[1]) {
    htrdr_log_warn(htrdr,
      "%s: the submitted spectral range overflowed the spectral data.\n", FUNC_NAME);
  }

  htrdr->wlen_range_m[0] = spectral_range[0]*1e-9; /* Convert in meters */
  htrdr->wlen_range_m[1] = spectral_range[1]*1e-9; /* Convert in meters */

  if(htrdr->spectral_type == HTRDR_SPECTRAL_SW_CIE_XYZ) {
    size_t n;
    n = (size_t)(spectral_range[1] - spectral_range[0]);
    res = htrdr_cie_xyz_create(htrdr, spectral_range, n, &htrdr->cie);
    if(res != RES_OK) goto error;

  } else {
    size_t n;

    if(htrdr->ref_temperature <= 0) {
      htrdr_log_err(htrdr, "%s: invalid reference temperature %g K.\n",
        FUNC_NAME, htrdr->ref_temperature);
      res = RES_BAD_ARG;
      goto error;
    }

    ASSERT(htrdr->wlen_range_m[0] <= htrdr->wlen_range_m[1]);
    n = (size_t)(spectral_range[1] - spectral_range[0]);

    res = htrdr_ran_wlen_create
      (htrdr, spectral_range, n, htrdr->ref_temperature, &htrdr->ran_wlen);
    if(res != RES_OK) goto error;
  }

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

  /* Create the image buffer only on the master process; the image parts
   * rendered by the processes are gathered onto the master process. */
  if(!htrdr->dump_vtk && htrdr->mpi_rank == 0) {
    const size_t pixsz = htrdr_spectral_type_get_pixsz
      (htrdr->spectral_type, htrdr->sensor.type);
    const size_t pixal = htrdr_spectral_type_get_pixal
      (htrdr->spectral_type, htrdr->sensor.type);

    res = htrdr_buffer_create(htrdr,
      args->image.definition[0], /* Width */
      args->image.definition[1], /* Height */
      args->image.definition[0] * pixsz, /* Pitch */
      pixsz, /* Size of a pixel */
      pixal, /* Alignment of a pixel */
      &htrdr->buf);
    if(res != RES_OK) goto error;
  }

exit:
  return res;
error:
  htrdr_release(htrdr);
  goto exit;
}

void
htrdr_release(struct htrdr* htrdr)
{
  ASSERT(htrdr);
  release_mpi(htrdr);
  if(htrdr->s3d) S3D(device_ref_put(htrdr->s3d));
  if(htrdr->ground) htrdr_ground_ref_put(htrdr->ground);
  if(htrdr->sky) HTSKY(ref_put(htrdr->sky));
  if(htrdr->sun) htrdr_sun_ref_put(htrdr->sun);
  if(htrdr->sensor.camera) htrdr_camera_ref_put(htrdr->sensor.camera);
  if(htrdr->sensor.rectangle) htrdr_rectangle_ref_put(htrdr->sensor.rectangle);
  if(htrdr->buf) htrdr_buffer_ref_put(htrdr->buf);
  if(htrdr->mtl) htrdr_mtl_ref_put(htrdr->mtl);
  if(htrdr->cie) htrdr_cie_xyz_ref_put(htrdr->cie);
  if(htrdr->ran_wlen) htrdr_ran_wlen_ref_put(htrdr->ran_wlen);
  if(htrdr->output && htrdr->output != stdout) fclose(htrdr->output);
  if(htrdr->lifo_allocators) {
    size_t i;
    FOR_EACH(i, 0, htrdr->nthreads) {
      mem_shutdown_lifo_allocator(&htrdr->lifo_allocators[i]);
    }
    MEM_RM(htrdr->allocator, htrdr->lifo_allocators);
  }
  str_release(&htrdr->output_name);
  logger_release(&htrdr->logger);
}

res_T
htrdr_run(struct htrdr* htrdr)
{
  res_T res = RES_OK;
  if(htrdr->dump_vtk) {
    const size_t nbands = htsky_get_spectral_bands_count(htrdr->sky);
    size_t i;

    /* Nothing to do */
    if(htrdr->mpi_rank != 0) goto exit;

    FOR_EACH(i, 0, nbands) {
      const size_t iband = htsky_get_spectral_band_id(htrdr->sky, i);
      const size_t nquads = htsky_get_spectral_band_quadrature_length
        (htrdr->sky, iband);
      size_t iquad;
      FOR_EACH(iquad, 0, nquads) {
        res = htsky_dump_cloud_vtk(htrdr->sky, iband, iquad, htrdr->output);
        if(res != RES_OK) goto error;
        fprintf(htrdr->output, "---\n");
      }
    }
  } else {
    res = htrdr_draw_map(htrdr, &htrdr->sensor, htrdr->width, htrdr->height,
      htrdr->spp, htrdr->buf);
    if(res != RES_OK) goto error;
    if(htrdr->mpi_rank == 0) {
      struct htrdr_accum path_time_acc = HTRDR_ACCUM_NULL;
      struct htrdr_accum flux_acc = HTRDR_ACCUM_NULL;
      struct htrdr_estimate path_time;
      struct htrdr_estimate flux;

      res = dump_buffer(htrdr, htrdr->buf, &path_time_acc, &flux_acc,
        str_cget(&htrdr->output_name), htrdr->output);
      if(res != RES_OK) goto error;

      htrdr_accum_get_estimation(&path_time_acc, &path_time);
      htrdr_log(htrdr,
        "Time per radiative path (in micro seconds): %g +/- %g\n",
        path_time.E,
        path_time.SE);

      if(htrdr->sensor.type == HTRDR_SENSOR_RECTANGLE) {
        htrdr_accum_get_estimation(&flux_acc, &flux);
        htrdr_log(htrdr,
          "Radiative flux density (in W/(external m^2)): %g +/- %g\n",
          flux.E,
          flux.SE);
      }
    }
  }
exit:
  return res;
error:
  goto exit;
}

void
htrdr_log(struct htrdr* htrdr, const char* msg, ...)
{
  ASSERT(htrdr && msg);
  /* Log standard message only on master process */
  if(htrdr->mpi_rank == 0) {
    va_list vargs_list;
    va_start(vargs_list, msg);
    log_msg(htrdr, LOG_OUTPUT, msg, vargs_list);
    va_end(vargs_list);
  }
}

void
htrdr_log_err(struct htrdr* htrdr, const char* msg, ...)
{
  va_list vargs_list;
  ASSERT(htrdr && msg);
  /* Log errors on all processes */
  va_start(vargs_list, msg);
  log_msg(htrdr, LOG_ERROR, msg, vargs_list);
  va_end(vargs_list);
}

void
htrdr_log_warn(struct htrdr* htrdr, const char* msg, ...)
{
  ASSERT(htrdr && msg);
  /* Log warnings only on master process */
  if(htrdr->mpi_rank == 0) {
    va_list vargs_list;
    va_start(vargs_list, msg);
    log_msg(htrdr, LOG_WARNING, msg, vargs_list);
    va_end(vargs_list);
  }
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
double
compute_sky_min_band_len
  (struct htsky* sky,
   const double range[2])
{
  double min_band_len = DBL_MAX;
  size_t nbands;
  ASSERT(sky && range && range[0] <= range[1]);

  nbands = htsky_get_spectral_bands_count(sky);

  if(eq_eps(range[0], range[1], 1.e-6)) {
    ASSERT(nbands == 1);
    min_band_len = 0;
  } else {
    size_t i = 0;

    /* Compute the length of the current band clamped to the submitted range */
    FOR_EACH(i, 0, nbands) {
      const size_t iband = htsky_get_spectral_band_id(sky, i);
      double wlens[2];
      HTSKY(get_spectral_band_bounds(sky, iband, wlens));

      /* Adjust band boundaries to the submitted range */
      wlens[0] = MMAX(wlens[0], range[0]);
      wlens[1] = MMIN(wlens[1], range[1]);

      min_band_len = MMIN(wlens[1] - wlens[0], min_band_len);
    }
  }
  return min_band_len;
}

res_T
open_output_stream
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
send_mpi_progress
  (struct htrdr* htrdr, const enum htrdr_mpi_message msg, int32_t percent)
{
  ASSERT(htrdr);
  ASSERT(msg == HTRDR_MPI_PROGRESS_RENDERING
      || msg == HTRDR_MPI_PROGRESS_BUILD_OCTREE);
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
    case HTRDR_MPI_PROGRESS_BUILD_OCTREE:
      progress = htrdr->mpi_progress_octree;
      break;
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
      case HTRDR_MPI_PROGRESS_BUILD_OCTREE:
        htrdr_fprintf(htrdr, stderr, "\033[2K\rBuilding octree: %3d%%",
          htrdr->mpi_progress_octree[0]);
        break;
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
        case HTRDR_MPI_PROGRESS_BUILD_OCTREE:
          htrdr_fprintf(htrdr, stderr,
            "\033[2K\rProcess %d -- building octree: %3d%%%c",
            iproc, htrdr->mpi_progress_octree[iproc],
            iproc == htrdr->mpi_nprocs - 1 ? '\r' : '\n');
          break;
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
    case HTRDR_MPI_PROGRESS_BUILD_OCTREE:
      progress = htrdr->mpi_progress_octree;
      break;
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

