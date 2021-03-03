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

#include "combustion/htrdr_combustion.h"
#include "combustion/htrdr_combustion_args.h"

#include "core/htrdr_log.h"

#include <rsys/mem_allocator.h>

int
htrdr_combustion_main(int argc, char** argv)
{
  char cmd_name[] = "htrdr-combustion";
  struct htrdr_combustion_args cmd_args = HTRDR_COMBUSTION_ARGS_DEFAULT;
  const size_t memsz_begin = MEM_ALLOCATED_SIZE(&mem_default_allocator);
  size_t memsz_end;
  res_T res = RES_OK;
  int err = 0;

  /* Overwrite command name */
  argv[0] = cmd_name;

  res = htrdr_combustion_args_init(&cmd_args, argc, argv);
  if(res != RES_OK) goto error;
  if(cmd_args.quit) goto exit;

exit:
  htrdr_combustion_args_release(&cmd_args);

  /* Check memory leaks */
  memsz_end = MEM_ALLOCATED_SIZE(&mem_default_allocator);
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

