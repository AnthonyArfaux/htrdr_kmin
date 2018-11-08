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
#include "htrdr_args.h"

#include <mpi.h>
#include <rsys/mem_allocator.h>

static const char*
thread_support_string(const int val)
{
  switch(val) {
    case MPI_THREAD_SINGLE: return "MPI_THREAD_SINGLE";
    case MPI_THREAD_FUNNELED: return "MPI_THREAD_FUNNELED";
    case MPI_THREAD_SERIALIZED: return "MPI_THREAD_SERIALIZED";
    case MPI_THREAD_MULTIPLE: return "MPI_THREAD_MULTIPLE";
    default: FATAL("Unreachable code.\n"); break;
  }
}

/*******************************************************************************
 * Program
 ******************************************************************************/
int
main(int argc, char** argv)
{
  struct htrdr htrdr;
  struct htrdr_args args = HTRDR_ARGS_DEFAULT;
  size_t memsz = 0;
  int err = 0;
  int is_htrdr_init = 0;
  int thread_support = 0;
  res_T res = RES_OK;

  err = MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &thread_support);
  if(err != MPI_SUCCESS) {
    fprintf(stderr, "Error initializing MPI.\n");
    goto error;
  }

  if(thread_support != MPI_THREAD_SERIALIZED) {
    fprintf(stderr, "The provided MPI implementation does not support "
      "serialized API calls from multiple threads. Provided thread support: "
      "%s.\n", thread_support_string(thread_support));
    goto error;
  }

  res = htrdr_args_init(&args, argc, argv);
  if(res != RES_OK) goto error;
  if(args.quit) goto exit;

  if(args.dump_vtk) {
    int rank;
    CHK(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);
    if(rank != 0) goto exit; /* Nothing to do except for the master process */
  }

  res = htrdr_init(NULL, &args, &htrdr);
  if(res != RES_OK) goto error;
  is_htrdr_init = 1;

  res = htrdr_run(&htrdr);
  if(res != RES_OK) goto error;

exit:
  MPI_Finalize();
  if(is_htrdr_init) htrdr_release(&htrdr);
  htrdr_args_release(&args);
  if((memsz = mem_allocated_size()) != 0) {
    fprintf(stderr, "Memory leaks: %lu Bytes\n", (unsigned long)memsz);
    err = -1;
  }
  return err;
error:
  err = -1;
  goto exit;
}
