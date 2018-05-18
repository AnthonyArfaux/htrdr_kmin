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
#include "htrdr_clouds.h"
#include "htrdr_rectangle.h"
#include "htrdr_solve.h"

#include <rsys/mem_allocator.h>

#include <star/svx.h>

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h> /* open */
#include <sys/stat.h> /* S_IRUSR & S_IWUSR */

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
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
      fp = fdopen(fd, "l");
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
release_htrdr(ref_T* ref)
{
  struct htrdr* htrdr = CONTAINER_OF(ref, struct htrdr, ref);
  ASSERT(htrdr);
  if(htrdr->svx) SVX(device_ref_put(htrdr->svx));
  if(htrdr->clouds) SVX(tree_ref_put(htrdr->clouds));
  if(htrdr->rect) htrdr_rectangle_ref_put(htrdr->rect);
  logger_release(&htrdr->logger);
  MEM_RM(htrdr->allocator, htrdr);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_create
  (struct mem_allocator* mem_allocator,
   const struct htrdr_args* args,
   struct htrdr** out_htrdr)
{
  struct htrdr* htrdr = NULL;
  struct mem_allocator* allocator = NULL;
  res_T res = RES_OK;
  ASSERT(args && out_htrdr);

  allocator = mem_allocator ? mem_allocator : &mem_default_allocator;
  htrdr = MEM_CALLOC(allocator, 1, sizeof(*htrdr));
  if(!htrdr) {
    fprintf(stderr, "Could not allocat the htrdr main structure.\n");
    goto error;
  }
  ref_init(&htrdr->ref);
  htrdr->allocator = allocator;

  logger_init(htrdr->allocator, &htrdr->logger);
  logger_set_stream(&htrdr->logger, LOG_ERROR, print_err, NULL);
  logger_set_stream(&htrdr->logger, LOG_WARNING, print_warn, NULL);

  htrdr->dump_vtk = args->dump_vtk;
  htrdr->verbose = args->verbose;

  /* Integration plane */
  res = htrdr_rectangle_create(htrdr, args->rectangle.pos, args->rectangle.tgt,
    args->rectangle.up, args->rectangle.sz, &htrdr->rect);
  if(res != RES_OK) goto error;

  res = svx_device_create(&htrdr->logger, allocator, args->verbose, &htrdr->svx);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "could not create the Star-VoXel device.\n");
    goto error;
  }

  if(!args->output) {
    htrdr->output = stdout;
  } else {
    res = open_output_stream(htrdr, args);
    if(res != RES_OK) goto error;
  }
  res = clouds_load(htrdr, args->verbose, args->input);
  if(res != RES_OK) goto error;

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

res_T
htrdr_run(struct htrdr* htrdr)
{
  res_T res = RES_OK;
  if(htrdr->dump_vtk) {
    res = clouds_dump_vtk(htrdr, htrdr->output);
    if(res != RES_OK) goto error;
  } else {
    double pos[3] = {3.4,2.2,0};
    double dir[3] = {0,0,1};
    double val = 0;
    res = htrdr_solve_transmission(htrdr, pos, dir, &val);
    if(res != RES_OK) goto error;
    printf(">>>> T = %g\n", val);
  }
exit:
  return res;
error:
  goto exit;
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

