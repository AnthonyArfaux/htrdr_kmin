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

#define _POSIX_C_SOURCE 200809L /* stat.st_time support */

#include "htrdr.h"
#include "htrdr_c.h"
#include "htrdr_args.h"
#include "htrdr_buffer.h"
#include "htrdr_camera.h"
#include "htrdr_ground.h"
#include "htrdr_sky.h"
#include "htrdr_sun.h"
#include "htrdr_solve.h"

#include <rsys/clock_time.h>
#include <rsys/cstr.h>
#include <rsys/mem_allocator.h>
#include <rsys/str.h>

#include <star/s3d.h>
#include <star/ssf.h>
#include <star/svx.h>

#include <errno.h>
#include <fcntl.h> /* open */
#include <libgen.h> /* basename */
#include <stdarg.h>
#include <stdio.h>
#include <unistd.h>
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
  if(htrdr->verbose && htrdr->mpi_rank == 0) {
    CHK(logger_vprint(&htrdr->logger, stream, msg, vargs) == RES_OK);
  }
}

static res_T
dump_accum_buffer
  (struct htrdr* htrdr,
   struct htrdr_buffer* buf,
   const char* stream_name,
   FILE* stream)
{
  struct htrdr_buffer_layout layout;
  size_t x, y;
  res_T res = RES_OK;
  ASSERT(htrdr && buf && stream_name && stream);
  (void)stream_name;

  htrdr_buffer_get_layout(buf, &layout);
  if(layout.elmt_size != sizeof(struct htrdr_accum[3])/*#channels*/
  || layout.alignment < ALIGNOF(struct htrdr_accum[3])) {
    htrdr_log_err(htrdr,
      "%s: invalid buffer layout. "
      "The pixel size must be the size of an accumulator.\n",
      FUNC_NAME);
    res = RES_BAD_ARG;
    goto error;
  }

  fprintf(stream, "%lu %lu\n", layout.width, layout.height);
  FOR_EACH(y, 0, layout.height) {
    FOR_EACH(x, 0, layout.width) {
      const struct htrdr_accum* accums = htrdr_buffer_at(buf, x, y);
      int i;
      FOR_EACH(i, 0, 3) {
        const double N = (double)accums[i].nweights;
        double E = 0;
        double V = 0;
        double SE = 0;

        if(accums[i].nweights) {
          E = accums[i].sum_weights / N;
          V = MMAX(accums[i].sum_weights_sqr / N - E*E, 0);
          SE = sqrt(V/N);
        }
        fprintf(stream, "%g %g ", E, SE);
      }
      fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
  }

exit:
  return res;
error:
  goto exit;
}

static res_T
open_file_stamp
  (struct htrdr* htrdr,
   const char* filename,
   struct stat* out_stat, /* Stat of the submitted filename */
   int* out_fd, /* Descriptor of the opened file. Must be closed by the caller */
   struct str* stamp_filename)
{
  struct stat statbuf;
  struct str str;
  int err;
  int fd = -1;
  res_T res = RES_OK;
  ASSERT(htrdr && filename && out_fd && out_stat && stamp_filename);

  str_init(htrdr->allocator, &str);

  err = stat(filename, &statbuf);
  if(err) {
    htrdr_log_err(htrdr, "%s: could not stat the file -- %s.\n",
      filename, strerror(errno));
    res = RES_IO_ERR;
    goto error;
  }

  if(!S_ISREG(statbuf.st_mode)) {
    htrdr_log_err(htrdr, "%s: not a regular file.\n", filename);
    res = RES_IO_ERR;
    goto error;
  }

  res = create_directory(htrdr, ".htrdr/");
  if(res != RES_OK) goto error;

  #define CHK_STR(Func, ErrMsg) {                                              \
    res = str_##Func;                                                          \
    if(res != RES_OK) {                                                        \
      htrdr_log_err(htrdr, "%s: "ErrMsg"\n", filename);                        \
      goto error;                                                              \
    }                                                                          \
  } (void)0
  CHK_STR(set(&str, filename), "could not copy the filename");
  CHK_STR(set(&str, basename(str_get(&str))), "could not setup the basename");
  CHK_STR(insert(&str, 0, ".htrdr/"), "could not setup the stamp directory");
  CHK_STR(append(&str, ".stamp"), "could not setup the stamp extension");
  #undef CHK_STR

  fd = open(str_cget(&str), O_CREAT|O_RDWR, S_IRUSR|S_IWUSR);
  if(fd < 0) {
    htrdr_log_err(htrdr, "%s: could not open/create the file -- %s.\n",
      str_cget(&str), strerror(errno));
    res = RES_IO_ERR;
    goto error;
  }

  CHK(str_copy_and_clear(stamp_filename, &str) == RES_OK);

exit:
  str_release(&str);
  *out_fd = fd;
  *out_stat = statbuf;
  return res;
error:
  if(fd >= 0) {
    CHK(close(fd) == 0);
    fd = -1;
  }
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

  /* Allocate #processes progress statuses on the master process and only 1
   * progress status on the other ones: the master process will gather the
   * status of the other processes to report their progression. */
  n = (size_t)(htrdr->mpi_rank == 0 ? htrdr->mpi_nprocs : 1);

  htrdr->mpi_progress_octree = MEM_CALLOC
    (htrdr->allocator, n, sizeof(*htrdr->mpi_progress_octree));
  if(!htrdr->mpi_progress_octree) {
    htrdr_log_err(htrdr,
      "could not allocate the progress state of the octree building.\n");
    res = RES_MEM_ERR;
    goto error;
  }

  htrdr->mpi_progress_render = MEM_CALLOC
    (htrdr->allocator, n, sizeof(*htrdr->mpi_progress_render));
  if(!htrdr->mpi_progress_render) {
    htrdr_log_err(htrdr,
      "could not allocate the progress state of the scene rendering.\n");
    res = RES_MEM_ERR;
    goto error;
  }

exit:
  return res;
error:
  if(htrdr->mpi_err_str) {
    MEM_RM(htrdr->allocator, htrdr->mpi_err_str);
    htrdr->mpi_err_str = NULL;
  }
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
  double proj_ratio;
  double sun_dir[3];
  const char* output_name = NULL;
  size_t ithread;
  res_T res = RES_OK;
  ASSERT(args && htrdr);

  memset(htrdr, 0, sizeof(*htrdr));

  htrdr->allocator = mem_allocator ? mem_allocator : &mem_default_allocator;

  logger_init(htrdr->allocator, &htrdr->logger);
  logger_set_stream(&htrdr->logger, LOG_OUTPUT, print_out, NULL);
  logger_set_stream(&htrdr->logger, LOG_ERROR, print_err, NULL);
  logger_set_stream(&htrdr->logger, LOG_WARNING, print_warn, NULL);

  str_init(htrdr->allocator, &htrdr->output_name);

  htrdr->dump_vtk = args->dump_vtk;
  htrdr->cache_grids = args->cache_grids;
  htrdr->verbose = args->verbose;
  htrdr->nthreads = MMIN(args->nthreads, (unsigned)omp_get_num_procs());
  htrdr->spp = args->image.spp;

  res = init_mpi(htrdr);
  if(res != RES_OK) goto error;

  if(htrdr->cache_grids && htrdr->mpi_nprocs != 1) {
    htrdr_log_warn(htrdr, "cached grids are not supported in a MPI execution.\n");
    htrdr->cache_grids = 0;
  }

  if(!args->output) {
    htrdr->output = stdout;
    output_name = "<stdout>";
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

  res = svx_device_create
    (&htrdr->logger, htrdr->allocator, args->verbose, &htrdr->svx);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: could not create the Star-VoXel device -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
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

  res = htrdr_ground_create
    (htrdr, args->filename_obj, args->repeat_ground, &htrdr->ground);
  if(res != RES_OK) goto error;

  proj_ratio =
    (double)args->image.definition[0]
  / (double)args->image.definition[1];
  res = htrdr_camera_create(htrdr, args->camera.pos, args->camera.tgt,
    args->camera.up, proj_ratio, MDEG2RAD(args->camera.fov_y), &htrdr->cam);
  if(res != RES_OK) goto error;

  res = htrdr_buffer_create(htrdr,
    args->image.definition[0], /* Width */
    args->image.definition[1], /* Height */
    args->image.definition[0]*sizeof(struct htrdr_accum[3]), /* Pitch */
    sizeof(struct htrdr_accum[3]),
    ALIGNOF(struct htrdr_accum[3]), /* Alignment */
    &htrdr->buf);
  if(res != RES_OK) goto error;

  res = htrdr_sun_create(htrdr, &htrdr->sun);
  if(res != RES_OK) goto error;
  spherical_to_cartesian_dir
    (MDEG2RAD(args->sun_azimuth), MDEG2RAD(args->sun_elevation), sun_dir);
  htrdr_sun_set_direction(htrdr->sun, sun_dir);

  res = htrdr_sky_create(htrdr, htrdr->sun, args->filename_les,
    args->filename_gas, args->filename_mie, args->optical_thickness,
    args->repeat_clouds, &htrdr->sky);
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
      (&htrdr->lifo_allocators[ithread], htrdr->allocator, 4096);
    if(res != RES_OK) {
      htrdr_log_err(htrdr,
        "%s: could not initialise the LIFO allocator of the thread %lu -- %s.\n",
        FUNC_NAME, (unsigned long)ithread, res_to_cstr(res));
      goto error;
    }
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
  if(htrdr->s3d) S3D(device_ref_put(htrdr->s3d));
  if(htrdr->svx) SVX(device_ref_put(htrdr->svx));
  if(htrdr->ground) htrdr_ground_ref_put(htrdr->ground);
  if(htrdr->sky) htrdr_sky_ref_put(htrdr->sky);
  if(htrdr->sun) htrdr_sun_ref_put(htrdr->sun);
  if(htrdr->cam) htrdr_camera_ref_put(htrdr->cam);
  if(htrdr->buf) htrdr_buffer_ref_put(htrdr->buf);
  if(htrdr->mpi_err_str) MEM_RM(htrdr->allocator, htrdr->mpi_err_str);
  if(htrdr->mpi_progress_octree) {
    MEM_RM(htrdr->allocator, htrdr->mpi_progress_octree);
  }
  if(htrdr->mpi_progress_render) {
    MEM_RM(htrdr->allocator, htrdr->mpi_progress_render);
  }
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
    const size_t nbands = htrdr_sky_get_sw_spectral_bands_count(htrdr->sky);
    size_t i;

    /* Nothing to do */
    if(htrdr->mpi_rank != 0) goto exit;

    FOR_EACH(i, 0, nbands) {
      const size_t iband = htrdr_sky_get_sw_spectral_band_id(htrdr->sky, i);
      const size_t nquads = htrdr_sky_get_sw_spectral_band_quadrature_length
        (htrdr->sky, iband);
      size_t iquad;
      FOR_EACH(iquad, 0, nquads) {
        res = htrdr_sky_dump_clouds_vtk(htrdr->sky, iband, iquad, htrdr->output);
        if(res != RES_OK) goto error;
        fprintf(htrdr->output, "---\n");
      }
    }
  } else {
    struct time t0, t1;
    char buf[128];

    time_current(&t0);
    res = htrdr_draw_radiance_sw(htrdr, htrdr->cam, htrdr->spp, htrdr->buf);
    if(res != RES_OK) goto error;
    time_sub(&t0, time_current(&t1), &t0);
    time_dump(&t0, TIME_ALL, NULL, buf, sizeof(buf));
    htrdr_log(htrdr, "Rendering time: %s\n", buf);

    if(htrdr->mpi_rank == 0) {
      res = dump_accum_buffer
        (htrdr, htrdr->buf, str_cget(&htrdr->output_name), htrdr->output);
      if(res != RES_OK) goto error;
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
  va_list vargs_list;
  ASSERT(htrdr && msg);
  va_start(vargs_list, msg);
  log_msg(htrdr, LOG_OUTPUT, msg, vargs_list);
  va_end(vargs_list);
}

void
htrdr_log_err(struct htrdr* htrdr, const char* msg, ...)
{
  va_list vargs_list;
  ASSERT(htrdr && msg);
  va_start(vargs_list, msg);
  log_msg(htrdr, LOG_ERROR, msg, vargs_list);
  va_end(vargs_list);
}

void
htrdr_log_warn(struct htrdr* htrdr, const char* msg, ...)
{
  va_list vargs_list;
  ASSERT(htrdr && msg);
  va_start(vargs_list, msg);
  log_msg(htrdr, LOG_WARNING, msg, vargs_list);
  va_end(vargs_list);
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
extern LOCAL_SYM res_T
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

res_T
is_file_updated(struct htrdr* htrdr, const char* filename, int* out_upd)
{
  struct str stamp_filename;
  struct stat statbuf;
  ssize_t n;
  off_t size;
  struct timespec mtime;
  int fd = -1;
  int upd = 1;
  res_T res = RES_OK;
  ASSERT(htrdr && filename && out_upd);

  str_init(htrdr->allocator, &stamp_filename);

  res = open_file_stamp(htrdr, filename, &statbuf, &fd, &stamp_filename);
  if(res != RES_OK) goto error;

  n = read(fd, &mtime, sizeof(mtime));
  if(n < 0) {
    htrdr_log_err(htrdr, "%s: could not read the `mtime' data -- %s.\n",
      str_cget(&stamp_filename), strerror(errno));
    res = RES_IO_ERR;
    goto error;
  }

  upd = (size_t)n != sizeof(mtime)
      ||mtime.tv_nsec != statbuf.st_mtim.tv_nsec
      ||mtime.tv_sec  != statbuf.st_mtim.tv_sec;

  if(!upd) {
    n = read(fd, &size, sizeof(size));
    if(n < 0) {
      htrdr_log_err(htrdr, "%s: could not read the `size' data -- %s.\n",
        str_cget(&stamp_filename), strerror(errno));
      res = RES_IO_ERR;
      goto error;
    }
    upd = (size_t)n != sizeof(size) || statbuf.st_size != size;
  }

exit:
  *out_upd = upd;
  str_release(&stamp_filename);
  if(fd >= 0) CHK(close(fd) == 0);
  return res;
error:
  goto exit;
}


res_T
update_file_stamp(struct htrdr* htrdr, const char* filename)
{
  struct str stamp_filename;
  struct stat statbuf;
  int fd = -1;
  ssize_t n;
  res_T res = RES_OK;
  ASSERT(htrdr && filename);

  str_init(htrdr->allocator, &stamp_filename);

  res = open_file_stamp(htrdr, filename, &statbuf, &fd, &stamp_filename);
  if(res != RES_OK) goto error;

  #define CHK_IO(Func, ErrMsg) {                                               \
    if((Func) < 0) {                                                           \
      htrdr_log_err(htrdr, "%s: "ErrMsg" -- %s.\n",                            \
        str_cget(&stamp_filename), strerror(errno));                           \
      res = RES_IO_ERR;                                                        \
      goto error;                                                              \
    }                                                                          \
  } (void) 0

  CHK_IO(lseek(fd, 0, SEEK_SET), "could not rewind the file descriptor");

  /* NOTE: Ignore n >=0 but != sizeof(DATA). In such case stamp is currupted
   * and on the next invocation on the same filename, this function will
   * return 1 */
  n = write(fd, &statbuf.st_mtim, sizeof(statbuf.st_mtim));
  CHK_IO(n, "could not update the `mtime' data");
  n = write(fd, &statbuf.st_size, sizeof(statbuf.st_size));
  CHK_IO(n, "could not update the `size' data");

  CHK_IO(fsync(fd), "could not sync the file with storage device");

  #undef CHK_IO

exit:
  str_release(&stamp_filename);
  if(fd >= 0) CHK(close(fd) == 0);
  return res;
error:
  goto exit;
}

res_T
create_directory(struct htrdr* htrdr, const char* path)
{
  res_T res = RES_OK;
  int err;
  ASSERT(htrdr && path);

  err = mkdir(path, S_IRWXU);
  if(!err) goto exit;

  if(errno != EEXIST) {
    htrdr_log_err(htrdr, "cannot create the `%s' directory -- %s.\n",
      path, strerror(errno));
    res = RES_IO_ERR;
    goto error;
  } else {
    const int fd = open(path, O_DIRECTORY);
    if(fd < -1) {
      htrdr_log_err(htrdr, "cannot open the `%s' directory -- %s.\n",
        path, strerror(errno));
      res = RES_IO_ERR;
      goto error;
    }
    CHK(!close(fd));
  }
exit:
  return res;
error:
  goto exit;
}

void
fetch_mpi_progress(struct htrdr* htrdr, const enum htrdr_mpi_progress tag)
{
  int8_t* progress = NULL;
  int iproc;
  ASSERT(htrdr && htrdr->mpi_rank == 0);

  switch(tag) {
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
      int flag;

      CHK(MPI_Iprobe
        (iproc, tag, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE) == MPI_SUCCESS);
      if(flag == 0) break; /* No more message */

      CHK(MPI_Recv(&progress[iproc], sizeof(size_t), MPI_CHAR, iproc, tag,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE) == MPI_SUCCESS);
    }
  }
}

void
print_mpi_progress(struct htrdr* htrdr, const enum htrdr_mpi_progress tag)
{
  ASSERT(htrdr && htrdr->mpi_rank == 0);

  if(htrdr->mpi_nprocs == 1) {
    switch(tag) {
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
      switch(tag) {
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
clear_mpi_progress(struct htrdr* htrdr, const enum htrdr_mpi_progress tag)
{
  ASSERT(htrdr);
  (void)tag;
  if(htrdr->mpi_nprocs > 1) {
    htrdr_fprintf(htrdr, stderr, "\033[%dA", htrdr->mpi_nprocs-1);
  }
}

int8_t
total_mpi_progress(const struct htrdr* htrdr, const enum htrdr_mpi_progress tag)
{
  const int8_t* progress = NULL;
  int total = 0;
  int iproc;
  ASSERT(htrdr && htrdr->mpi_rank == 0);

  switch(tag) {
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
  ASSERT(total <= 100);
  return (int8_t)total;
}

