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

#define _POSIX_C_SOURCE 200809L /* mmap support */
#define _DEFAULT_SOURCE 1 /* MAP_POPULATE support */

#include "htrdr.h"
#include "htrdr_c.h"
#include "htrdr_grid.h"

#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

#include <errno.h>
#include <sys/mman.h> /* mmap/munmap */
#include <fcntl.h>
#include <unistd.h> /* sysconf */

struct htrdr_grid {
  FILE* fp;
  char* data;
  size_t definition[3];
  size_t cell_sz;
  size_t pagesize;

  ref_T ref;
  struct htrdr* htrdr;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
grid_release(ref_T* ref)
{
  struct htrdr_grid* grid;
  ASSERT(ref);
  grid = CONTAINER_OF(ref, struct htrdr_grid, ref);
  if(grid->fp) fclose(grid->fp);
  if(grid->data) {
    size_t grid_sz;
    grid_sz =
      grid->definition[0]
    * grid->definition[1]
    * grid->definition[2]
    * grid->cell_sz;
    grid_sz = ALIGN_SIZE(grid_sz, grid->pagesize);
    if(munmap(grid->data, grid_sz)) {
      htrdr_log_err(grid->htrdr, "error unmapping the grid data -- %s.\n",
        strerror(errno));
      ASSERT(0);
    }
  }
  MEM_RM(grid->htrdr->allocator, grid);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_grid_create
  (struct htrdr* htrdr,
   const size_t definition[3],
   const size_t sizeof_cell, /* Size of an cell in Bytes */
   const char* filename,
   const int force_overwrite,
   struct htrdr_grid** out_grid)
{
  const char byte = 0;
  struct htrdr_grid* grid = NULL;
  size_t grid_sz;
  long grid_offset;
  int n;
  res_T res = RES_OK;
  ASSERT(htrdr && out_grid && filename && definition);

  if(!definition[0] || !definition[1] || !definition[2]) {
    htrdr_log_err(htrdr, "%s: invalid definition [%lu, %lu, %lu].\n", FUNC_NAME,
      (unsigned long)definition[0],
      (unsigned long)definition[1],
      (unsigned long)definition[2]);
    res = RES_BAD_ARG;
    goto error;
  }

  if(!sizeof_cell) {
    htrdr_log_err(htrdr, "%s: invalid cell size `%lu'.\n", FUNC_NAME,
      (unsigned long)sizeof_cell);
    res = RES_BAD_ARG;
    goto error;
  }

  grid = MEM_CALLOC(htrdr->allocator, 1, sizeof(*grid));
  if(!grid) {
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&grid->ref);
  grid->definition[0] = definition[0];
  grid->definition[1] = definition[1];
  grid->definition[2] = definition[2];
  grid->cell_sz = sizeof_cell;
  grid->pagesize = (size_t)sysconf(_SC_PAGESIZE);
  grid->htrdr = htrdr;

  res = open_output_stream(htrdr, filename, 1, force_overwrite, &grid->fp);
  if(res != RES_OK) goto error;

  #define WRITE(Var, N, Name) {                                                \
    if(fwrite((Var), sizeof(*(Var)), (N), grid->fp) != (N)) {                  \
      htrdr_log_err(htrdr, "%s:%s: could not write `%s' -- %s.\n",             \
        FUNC_NAME, filename, (Name), strerror(errno));                         \
      res = RES_IO_ERR;                                                        \
      goto error;                                                              \
    }                                                                          \
  } (void)0
  WRITE(&grid->pagesize, 1, "pagesize");
  WRITE(&grid->cell_sz, 1, "cell_sz");
  WRITE(grid->definition, 3, "definition");

  /* Align the grid data on pagesize */
  n = fseek
    (grid->fp, ALIGN_SIZE(ftell(grid->fp),(off_t)grid->pagesize), SEEK_SET);
  if(n < 0) {
    htrdr_log_err(htrdr,
      "%s:%s: could not align the grid data on page size -- %s.\n",
      FUNC_NAME, filename, strerror(errno));
    res = RES_IO_ERR;
    goto error;
  }

  /* Define the grid size */
  grid_sz = definition[0] * definition[1] * definition[2] * sizeof_cell;
  grid_sz = ALIGN_SIZE(grid_sz, grid->pagesize);

  /* Save the position of the grid data into the file */
  grid_offset = ftell(grid->fp);

  /* Reserve the space for the grid data */
  n = fseek(grid->fp, (long)grid_sz, SEEK_CUR);
  if(n < 0) {
    htrdr_log_err(htrdr,
      "%s:%s: could reserve the space to store the grid -- %s.\n", FUNC_NAME,
      filename, strerror(errno));
    res = RES_IO_ERR;
    goto error;
  }

  /* Write one char at the end of the file to position the EOF indicator */
  CHK(fseek(grid->fp, -1, SEEK_CUR) != -1);
  WRITE(&byte, 1, "Dummy Byte");
  #undef WRITE

  /* Avoid to be positionned on the EOF */
  rewind(grid->fp);

  grid->data = mmap(NULL, grid_sz, PROT_READ|PROT_WRITE,
    MAP_SHARED|MAP_POPULATE, fileno(grid->fp), grid_offset);
  if(grid->data == MAP_FAILED) {
    htrdr_log_err(htrdr, "%s:%s: could not map the grid data -- %s.\n",
      FUNC_NAME, filename, strerror(errno));
    res = RES_IO_ERR;
    goto error;
  }

exit:
  *out_grid = grid;
  return res;
error:
  if(grid) {
    htrdr_grid_ref_put(grid);
    grid = NULL;
  }
  goto exit;
}

res_T
htrdr_grid_open
  (struct htrdr* htrdr,
   const char* filename,
   struct htrdr_grid** out_grid)
{
  struct htrdr_grid* grid = NULL;
  size_t grid_sz;
  size_t grid_offset;
  size_t pagesize;
  int fd = -1;
  res_T res = RES_OK;
  ASSERT(htrdr && filename && out_grid);

  grid = MEM_CALLOC(htrdr->allocator, 1, sizeof(*grid));
  if(!grid) {
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&grid->ref);
  grid->pagesize = (size_t)sysconf(_SC_PAGESIZE);
  grid->htrdr = htrdr;

  fd = open(filename, O_RDWR, 0);
  if(fd < 0) {
    htrdr_log_err(htrdr, "%s: could not open `%s' -- %s.\n", FUNC_NAME,
      filename, strerror(errno));
    res = RES_IO_ERR;
    goto error;
  }
  CHK(grid->fp = fdopen(fd, "rw"));

  #define READ(Var, N, Name) {                                                 \
    if(fread((Var), sizeof(*(Var)), (N), grid->fp) != (N)) {                   \
      htrdr_log_err(htrdr, "%s:%s: could not read `%s'.\n",                    \
        FUNC_NAME, filename, Name);                                            \
      res = RES_IO_ERR;                                                        \
      goto error;                                                              \
    }                                                                          \
  } (void)0
  READ(&pagesize, 1, "pagesize");
  if(pagesize != grid->pagesize) {
    htrdr_log_err(htrdr, "%s:%s: invalid pagesize `%lu'.\n", FUNC_NAME,
      filename, (unsigned long)pagesize);
    res = RES_BAD_ARG;
    goto error;
  }

  READ(&grid->cell_sz, 1, "sizeof_cell");
  if(grid->cell_sz == 0) {
    htrdr_log_err(htrdr, "%s:%s: invalid cell size `%lu'.\n", FUNC_NAME,
      filename, (unsigned long)grid->cell_sz);
    res = RES_BAD_ARG;
    goto error;
  }

  READ(grid->definition, 3, "definition");
  if(!grid->definition[0] || !grid->definition[1] || !grid->definition[2]) {
    htrdr_log_err(htrdr, "%s:%s: invalid definition [%lu, %lu, %lu].\n",
      FUNC_NAME, filename,
      (unsigned long)grid->definition[0],
      (unsigned long)grid->definition[1],
      (unsigned long)grid->definition[2]);
    res = RES_BAD_ARG;
    goto error;
  }

  grid_offset = ALIGN_SIZE((size_t)ftell(grid->fp), grid->pagesize);
  grid_sz =
    grid->definition[0]
  * grid->definition[1]
  * grid->definition[2]
  * grid->cell_sz;
  grid_sz = ALIGN_SIZE(grid_sz, grid->pagesize);

  grid->data = mmap(NULL, grid_sz, PROT_READ|PROT_WRITE,
    MAP_SHARED|MAP_POPULATE, fileno(grid->fp), (off_t)grid_offset);

  if(grid->data == MAP_FAILED) {
    htrdr_log_err(htrdr, "%s:%s: could not map the grid data -- %s.\n",
      FUNC_NAME, filename, strerror(errno));
    res = RES_IO_ERR;
    goto error;
  }

exit:
  *out_grid = grid;
  return res;
error:
  if(grid) {
    htrdr_grid_ref_put(grid);
    grid = NULL;
  }
  goto exit;
}

void
htrdr_grid_ref_get(struct htrdr_grid* grid)
{
  ASSERT(grid);
  ref_get(&grid->ref);
}

void
htrdr_grid_ref_put(struct htrdr_grid* grid)
{
  ASSERT(grid);
  ref_put(&grid->ref, grid_release);
}

void*
htrdr_grid_at(struct htrdr_grid* grid, const size_t xyz[3])
{
  size_t slice;
  size_t pitch;
  ASSERT(grid && xyz);
  ASSERT(xyz[0] < grid->definition[0]);
  ASSERT(xyz[1] < grid->definition[1]);
  ASSERT(xyz[2] < grid->definition[2]);
  pitch = grid->definition[0] * grid->cell_sz;
  slice = grid->definition[1] * pitch;
  return grid->data + xyz[2]*slice + xyz[1]*pitch + xyz[0]*grid->cell_sz;
}

