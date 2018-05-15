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

#include <rsys/mem_allocator.h>

/*******************************************************************************
 * Program
 ******************************************************************************/
int
main(int argc, char** argv)
{
  struct htrdr htrdr;
  struct htrdr_args args;
  size_t memsz = 0;
  int htrdr_is_init = 0;
  int err = 0;
  res_T res = RES_OK;

  res = htrdr_args_init(&args, argc, argv);
  if(res != RES_OK) goto error;
  if(args.quit) goto exit;

  res = htrdr_init(NULL, &args, &htrdr);
  if(res != RES_OK) goto error;
  htrdr_is_init = 1;

  res = htrdr_run(&htrdr);
  if(res != RES_OK) goto error;

exit:
  if(htrdr_is_init) htrdr_release(&htrdr);
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
