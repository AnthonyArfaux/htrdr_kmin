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
#include "htrdr_rectangle.h"
#include "htrdr_sky.h"
#include "htrdr_solve.h"

#include <rsys/clock_time.h>
#include <rsys/mem_allocator.h>

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

static void
dump_buffer(struct htrdr* htrdr)
{
  struct htrdr_buffer_layout buf_layout;
  size_t x, y;
  ASSERT(htrdr);

  htrdr_buffer_get_layout(htrdr->buf, &buf_layout);

  fprintf(htrdr->output, "P3 %lu %lu\n255\n",
    (unsigned long)buf_layout.width,
    (unsigned long)buf_layout.height);
  FOR_EACH(y, 0, buf_layout.height) {
    FOR_EACH(x, 0, buf_layout.width) {
      double val = *((double*)htrdr_buffer_at(htrdr->buf, x, y));
      int i;
      i = (int)(val * 255);
      fprintf(htrdr->output, "%d %d %d\n", i, i, i);
    }
  }
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
  res_T res = RES_OK;
  ASSERT(args && htrdr);

  memset(htrdr, 0, sizeof(*htrdr));

  htrdr->allocator = mem_allocator ? mem_allocator : &mem_default_allocator;

  logger_init(htrdr->allocator, &htrdr->logger);
  logger_set_stream(&htrdr->logger, LOG_OUTPUT, print_out, NULL);
  logger_set_stream(&htrdr->logger, LOG_ERROR, print_err, NULL);
  logger_set_stream(&htrdr->logger, LOG_WARNING, print_warn, NULL);

  htrdr->dump_vtk = args->dump_vtk;
  htrdr->verbose = args->verbose;
  htrdr->nthreads = MMIN(args->nthreads, (unsigned)omp_get_num_procs());

  res = svx_device_create
    (&htrdr->logger, htrdr->allocator, args->verbose, &htrdr->svx);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "could not create the Star-VoXel device.\n");
    goto error;
  }

  if(!args->dump_vtk) { /* Legacy mode */
    const size_t elmtsz = sizeof(double);
    const size_t pitch = elmtsz * args->image.definition[0];

    /* Create the image buffer */
    res = htrdr_buffer_create(htrdr, args->image.definition[0],
      args->image.definition[1], pitch, elmtsz, 16, &htrdr->buf);
    if(res != RES_OK) goto error;

    /* TODO check the validity of the parameters */
    htrdr->main_dir[0] = args->main_dir[0];
    htrdr->main_dir[1] = args->main_dir[1];
    htrdr->main_dir[2] = args->main_dir[2];
    htrdr->spp = args->image.spp;

    /* Create the plane on which the image buffer lies */
    res = htrdr_rectangle_create(htrdr, args->rectangle.pos, args->rectangle.tgt,
      args->rectangle.up, args->rectangle.sz, &htrdr->rect);
    if(res != RES_OK) goto error;
  }

  if(!args->output) {
    htrdr->output = stdout;
  } else {
    res = open_output_stream(htrdr, args);
    if(res != RES_OK) goto error;
  }
  res = htrdr_sky_create(htrdr, args->input, &htrdr->sky);
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
  if(htrdr->svx) SVX(device_ref_put(htrdr->svx));
  if(htrdr->sky) htrdr_sky_ref_put(htrdr->sky);
  if(htrdr->buf) htrdr_buffer_ref_put(htrdr->buf);
  if(htrdr->rect) htrdr_rectangle_ref_put(htrdr->rect);
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
    res = htrdr_solve_transmission_buffer(htrdr);
    if(res != RES_OK) goto error;
    time_sub(&t0, time_current(&t1), &t0);
    time_dump(&t0, TIME_ALL, NULL, buf, sizeof(buf));
    htrdr_log(htrdr, "Elapsed time: %s\n", buf);

    dump_buffer(htrdr);
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

