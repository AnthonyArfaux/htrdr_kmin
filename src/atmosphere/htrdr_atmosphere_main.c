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

#include "htrdr.h"
#include "htrdr_atmosphere.h"
#include "htrdr_atmosphere_args.h"

int
main(int argc, char** argv)
{
  struct htrdr_args htrdr_args = HTRDR_ARGS_DEFAULT;
  struct htrdr_atmosphere_args cmd_args = HTRDR_ATMOSPHERE_ARGS_DEFAULT;
  struct htrdr* htrdr = NULL;
  struct htrdr_atmosphere* cmd = NULL;
  size_t memsz;
  int err = 0;
  res_T res = RES_OK;

  res = htrdr_mpi_init(argc, argv);
  if(res != RES_OK) goto error;

  res = htrdr_atmosphere_args_init(&cmd.args, argc, argv);
  if(res != RES_OK) goto error;
  if(args.quit) goto exit;

  htrdr_args.nthreads = cmd_args.nthreads;
  htrdr_args.verbose = cmd_args.verbose;
  res = htrdr_create(NULL, &htrdr_args, &htrdr);
  if(res != RES_OK) goto error;

  if(cmd_args.dump_volumetric_acceleration_structure 
  && htrdr_get_mpi_rank(htrdr) != 0) {
    goto exit; /* Nothing to do except for the master process */
  }

  res = htrdr_atmosphere_create(htrdr, &cmd_args, &cmd);
  if(res != RES_OK) goto error;

  res = htrdr_atmosphere_run(cmd);
  if(res != RES_OK) goto error;

exit:
  htrdr_mpi_finalize();
  htrdr_atmosphere_args_release(&cmd_args);
  if(htrdr) htrdr_ref_put(htrdr);
  if(cmd) htrdr_atmosphere_ref_put(cmd);

  /* Check memory leaks */
  if((memsz = mem_allocated_size()) != 0) {
    fprintf(stderr, "Memory leaks: %lu Bytes\n", (unsigned long)memsz);
    err = -1;
  }
  return err;
error:
  err = - 1;
  goto exit;
}
