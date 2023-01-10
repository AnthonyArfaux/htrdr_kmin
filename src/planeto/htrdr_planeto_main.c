/* Copyright (C) 2018, 2019, 2020, 2021 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019, 2021 CNRS
 * Copyright (C) 2018, 2019, Université Paul Sabatier
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

#include "planeto/htrdr_planeto.h"
#include "planeto/htrdr_planeto_args.h"

#include "core/htrdr_log.h"

#include <rsys/mem_allocator.h>

int
htrdr_planeto_main(int argc, char** argv)
{
  char cmd_name[] = "htrdr-planeto";
  struct htrdr_args htrdr_args = HTRDR_ARGS_DEFAULT;
  struct htrdr_planeto_args cmd_args = HTRDR_PLANETO_ARGS_DEFAULT;
  struct htrdr* htrdr = NULL;
  struct htrdr_planeto* cmd = NULL;
  const size_t memsz_begin = mem_allocated_size();
  size_t memsz_end;
  int is_mpi_init = 0;
  res_T res = RES_OK;
  int err = 0;

  /* Overwrite command name */
  argv[0] = cmd_name;

  res = htrdr_mpi_init(argc, argv);
  if(res != RES_OK) goto error;
  is_mpi_init = 1;

  res = htrdr_planeto_args_init(&cmd_args, argc, argv);
  if(res != RES_OK) goto error;
  if(cmd_args.quit) goto exit;

  htrdr_args.nthreads = cmd_args.nthreads;
  htrdr_args.verbose = cmd_args.verbose;
  res = htrdr_create(&mem_default_allocator, &htrdr_args, &htrdr);
  if(res != RES_OK) goto error;

  if(cmd_args.output_type == HTRDR_PLANETO_ARGS_OUTPUT_OCTREES
  && htrdr_get_mpi_rank(htrdr) != 0) {
    goto exit; /* Nothing to do except for the master process */
  }

  res = htrdr_planeto_create(htrdr, &cmd_args, &cmd);
  if(res != RES_OK) goto error;

  res = htrdr_planeto_run(cmd);
  if(res != RES_OK) goto error;

exit:
  htrdr_planeto_args_release(&cmd_args);
  if(is_mpi_init) htrdr_mpi_finalize();
  if(htrdr) htrdr_ref_put(htrdr);
  if(cmd) htrdr_planeto_ref_put(cmd);

  /* Check memory leaks */
  memsz_end = mem_allocated_size();
  if(memsz_begin != memsz_end) {
    ASSERT(memsz_end >= memsz_begin);
    fprintf(stderr, HTRDR_LOG_WARNING_PREFIX"Memory leaks: %lu Bytes\n",
      (unsigned long)(memsz_end - memsz_begin));
    err = -1;
  }
  return err;
error:
  err = -1;
  goto exit;
}

