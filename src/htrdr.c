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

#define _POSIX_C_SOURCE 200112L /* open support */

#include "htrdr.h"
#include "htrdr_args.h"
#include "htrdr_buffer.h"
#include "htrdr_camera.h"
#include "htrdr_sky.h"
#include "htrdr_sun.h"
#include "htrdr_solve.h"

#include <rsys/clock_time.h>
#include <rsys/mem_allocator.h>

#include <star/s3d.h>
#include <star/s3daw.h>
#include <star/ssf.h>
#include <star/svx.h>

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h> /* open */
#include <sys/stat.h> /* S_IRUSR & S_IWUSR */

#include <omp.h>

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

static res_T
open_output_stream(struct htrdr* htrdr, const struct htrdr_args* args)
{
  FILE* fp = NULL;
  int fd = -1;
  res_T res = RES_OK;
  ASSERT(htrdr && args);

  if(args->force_overwriting) {
    fp = fopen(args->output, "w");
    if(!fp) {
      htrdr_log_err(htrdr,
        "could not open the output file `%s'.\n", args->output);
      goto error;
    }
  } else {
    fd = open(args->output, O_CREAT|O_WRONLY|O_EXCL|O_TRUNC, S_IRUSR|S_IWUSR);
    if(fd >= 0) {
      fp = fdopen(fd, "w");
      if(fp == NULL) {
        htrdr_log_err(htrdr,
          "could not open the output file `%s'.\n", args->output);
        goto error;
      }
    } else if(errno == EEXIST) {
      htrdr_log_err(htrdr,
        "the output file `%s' already exists. Use -f to overwrite it.\n",
        args->output);
      goto error;
    } else {
      htrdr_log_err(htrdr,
        "unexpected error while opening the output file `%s'.\n", args->output);
      goto error;
    }
  }
exit:
  htrdr->output = fp;
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
  if(layout.elmt_size != sizeof(struct htrdr_accum)
  || layout.alignment < ALIGNOF(struct htrdr_accum)) {
    htrdr_log_err(htrdr,
      "%s: invalid buffer layout. "
      "The pixel size must be the size of an accumulator.\n",
      FUNC_NAME);
    res = RES_BAD_ARG;
    goto error;
  }

  FOR_EACH(y, 0, layout.height) {
    FOR_EACH(x, 0, layout.width) {
      const struct htrdr_accum* accum = htrdr_buffer_at(buf, x, y);
      const double E = accum->nweights
        ? accum->sum_weights / (double)accum->nweights : 0;
      fprintf(stream, "%g ", E);
    }
    fprintf(stream, "\n");
  }

exit:
  return res;
error:
  goto exit;
}

static res_T
setup_geometry(struct htrdr* htrdr, const char* filename)
{
  struct s3d_scene* scn = NULL;
  struct s3daw* s3daw = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr);

  res = s3d_scene_create(htrdr->s3d, &scn);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "could not create the Star-3D scene.\n");
    goto error;
  }

  if(filename) {
    res = s3daw_create(&htrdr->logger, htrdr->allocator, NULL, NULL, htrdr->s3d,
      htrdr->verbose, &s3daw);
    if(res != RES_OK) {
      htrdr_log_err(htrdr, "could not create the Star-3DAW device.\n");
      goto error;
    }
    res = s3daw_load(s3daw, filename);
    if(res != RES_OK) {
      htrdr_log_err(htrdr, "could not load the obj file `%s'.\n", filename);
      goto error;
    }
    res = s3daw_attach_to_scene(s3daw, scn);
    if(res != RES_OK) {
      htrdr_log_err(htrdr,
        "could not attach the loaded geometry to the Star-3D scene.\n");
      goto error;
    }
  }

  res = s3d_scene_view_create(scn, S3D_TRACE, &htrdr->s3d_scn_view);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "could not create the Star-3D scene view.\n");
    goto error;
  }

exit:
  if(scn) S3D(scene_ref_put(scn));
  if(s3daw) S3DAW(ref_put(s3daw));
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
  double proj_ratio;
  const char* output_name = NULL;
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
  htrdr->verbose = args->verbose;
  htrdr->nthreads = MMIN(args->nthreads, (unsigned)omp_get_num_procs());
  htrdr->spp = args->image.spp;

  if(!args->output) {
    htrdr->output = stdout;
    output_name = "<stdout>";
  } else {
    res = open_output_stream(htrdr, args);
    if(res != RES_OK) goto error;
    output_name = args->output;
  }
  res = str_set(&htrdr->output_name, output_name);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "could not store the name of the output stream `%s'.\n", output_name);
    goto error;
  }

  res = svx_device_create
    (&htrdr->logger, htrdr->allocator, args->verbose, &htrdr->svx);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "could not create the Star-VoXel device.\n");
    goto error;
  }

  /* Disable the Star-3D verbosity since the Embree backend print some messages
   * on stdout rather than stderr. This is annoying since stdout may be used by
   * htrdr to write output data */
  res = s3d_device_create
    (&htrdr->logger, htrdr->allocator, 0, &htrdr->s3d);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "could not create the Star-3D device.\n");
    goto error;
  }

  proj_ratio =
    (double)args->image.definition[0]
  / (double)args->image.definition[1];
  res = htrdr_camera_create(htrdr, args->camera.pos, args->camera.tgt,
    args->camera.up, proj_ratio, MDEG2RAD(args->camera.fov_x), &htrdr->cam);
  if(res != RES_OK) goto error;

  res = htrdr_buffer_create(htrdr,
    args->image.definition[0], /* Width */
    args->image.definition[1], /* Height */
    args->image.definition[0]*sizeof(struct htrdr_accum), /* Pitch */
    sizeof(struct htrdr_accum), /* Element size */
    16, /* Alignment */
    &htrdr->buf);
  if(res != RES_OK) goto error;

  res = htrdr_sun_create(htrdr, &htrdr->sun);
  if(res != RES_OK) goto error;
  htrdr_sun_set_direction(htrdr->sun, args->main_dir);

  res = htrdr_sky_create(htrdr, htrdr->sun, args->filename_les,
    args->filename_mie, &htrdr->sky);
  if(res != RES_OK) goto error;

  res = ssf_bsdf_create
    (htrdr->allocator, &ssf_lambertian_reflection, &htrdr->bsdf);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "could not create the BSDF of the ground.\n");
    goto error;
  }
  SSF(lambertian_reflection_setup(htrdr->bsdf, 1));

  res = ssf_phase_create(htrdr->allocator, &ssf_phase_hg, &htrdr->phase_hg);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "could not create the Henyey & Greenstein phase function.\n");
    goto error;
  }
  SSF(phase_hg_setup(htrdr->phase_hg, 0));

  res = ssf_phase_create
    (htrdr->allocator, &ssf_phase_rayleigh, &htrdr->phase_rayleigh);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,"could not create the Rayleigh phase function.\n");
    goto error;
  }

  res = setup_geometry(htrdr, args->filename_obj);
  if(res != RES_OK) goto error;

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
  if(htrdr->bsdf) SSF(bsdf_ref_put(htrdr->bsdf));
  if(htrdr->phase_hg) SSF(phase_ref_put(htrdr->phase_hg));
  if(htrdr->phase_rayleigh) SSF(phase_ref_put(htrdr->phase_rayleigh));
  if(htrdr->s3d_scn_view) S3D(scene_view_ref_put(htrdr->s3d_scn_view));
  if(htrdr->s3d) S3D(device_ref_put(htrdr->s3d));
  if(htrdr->svx) SVX(device_ref_put(htrdr->svx));
  if(htrdr->sky) htrdr_sky_ref_put(htrdr->sky);
  if(htrdr->sun) htrdr_sun_ref_put(htrdr->sun);
  if(htrdr->cam) htrdr_camera_ref_put(htrdr->cam);
  if(htrdr->buf) htrdr_buffer_ref_put(htrdr->buf);
  str_release(&htrdr->output_name);
  logger_release(&htrdr->logger);
}

res_T
htrdr_run(struct htrdr* htrdr)
{
  res_T res = RES_OK;
  if(htrdr->dump_vtk) {
    res = htrdr_sky_dump_clouds_vtk(htrdr->sky, htrdr->output);
    if(res != RES_OK) goto error;
  } else {
    struct time t0, t1;
    char buf[128];

    time_current(&t0);
    res = htrdr_draw_radiance_sw(htrdr, htrdr->cam, htrdr->spp, htrdr->buf);
    if(res != RES_OK) goto error;
    time_sub(&t0, time_current(&t1), &t0);
    time_dump(&t0, TIME_ALL, NULL, buf, sizeof(buf));
    htrdr_log(htrdr, "Elapsed time: %s\n", buf);

    res = dump_accum_buffer
      (htrdr, htrdr->buf, str_cget(&htrdr->output_name), htrdr->output);
    if(res != RES_OK) goto error;
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

