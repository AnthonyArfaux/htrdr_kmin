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

const int32_t GRID_VERSION = 0;
const int32_t GRID_VERSION_NONE = -1;

struct htrdr_grid {
  FILE* fp;
  char* data; /* Mapped data */
  size_t definition[3]; /* Submitted definition */
  size_t def_adjusted; /* Adjusted definition along the 3 dimensions */
  size_t cell_sz; /* Size in bytes of a grid cell */
  size_t pagesize; /* Page size in bytes */
  size_t data_sz; /* Size in bytes of the overall grid data + padding */

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
  if(grid->fp) {
    rewind(grid->fp);
    CHK(fwrite(&GRID_VERSION, sizeof(int), 1, grid->fp) == 1);
    fclose(grid->fp);
  }
  if(grid->data) {
    if(munmap(grid->data, grid->data_sz)) {
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
  size_t mcode_max;
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
  WRITE(&GRID_VERSION_NONE, 1, "version");
  WRITE(&grid->pagesize, 1, "pagesize");
  WRITE(&grid->cell_sz, 1, "cell_sz");
  WRITE(grid->definition, 3, "definition");
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

  /* Adjust the grid definition in order to sort its data wrt the morton code
   * of its voxel */
  grid->def_adjusted = MMAX(MMAX(definition[0], definition[1]), definition[2]);
  grid->def_adjusted = round_up_pow2(grid->def_adjusted);
  mcode_max = grid->def_adjusted*grid->def_adjusted*grid->def_adjusted;

  /* Define the grid size */
  grid->data_sz = mcode_max * sizeof_cell;
  grid->data_sz = ALIGN_SIZE(grid->data_sz, grid->pagesize);

  /* Save the position of the grid data into the file */
  grid_offset = ftell(grid->fp);

  /* Reserve the space for the grid data */
  n = fseek(grid->fp, (long)grid->data_sz, SEEK_CUR);
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

  /* Map the grid data */
  grid->data = mmap(NULL, grid->data_sz, PROT_READ|PROT_WRITE,
    MAP_SHARED|MAP_POPULATE, fileno(grid->fp), grid_offset);

  if(grid->data == MAP_FAILED) {
    htrdr_log_err(htrdr, "%s:%s: could not map the grid data -- %s.\n",
      FUNC_NAME, filename, strerror(errno));
    grid->data = NULL;
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
  size_t grid_offset;
  size_t pagesize;
  size_t mcode_max;
  int32_t version;
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
  CHK(grid->fp = fdopen(fd, "w+"));

  #define READ(Var, N, Name) {                                                 \
    if(fread((Var), sizeof(*(Var)), (N), grid->fp) != (N)) {                   \
      htrdr_log_err(htrdr, "%s:%s: could not read `%s'.\n",                    \
        FUNC_NAME, filename, Name);                                            \
      res = RES_IO_ERR;                                                        \
      goto error;                                                              \
    }                                                                          \
  } (void)0
  READ(&version, 1, "version");
  if(version != GRID_VERSION) {
    htrdr_log_err(htrdr, "%s:%s: incompatible grid version. Loaded version is "
      "'%i' while the current version is '%i'.\n",
      FUNC_NAME, filename, version, GRID_VERSION);
    res = RES_BAD_ARG;
    goto error;
  }

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
  #undef READ

  grid_offset = ALIGN_SIZE((size_t)ftell(grid->fp), grid->pagesize);
  grid->def_adjusted = MMAX(grid->definition[0], grid->definition[1]);
  grid->def_adjusted = MMAX(grid->definition[2], grid->def_adjusted);
  grid->def_adjusted = round_up_pow2(grid->def_adjusted);
  mcode_max = grid->def_adjusted*grid->def_adjusted*grid->def_adjusted;

  /* Define the grid size */
  grid->data_sz = mcode_max * grid->cell_sz;
  grid->data_sz = ALIGN_SIZE(grid->data_sz, grid->pagesize);

  grid->data = mmap(NULL, grid->data_sz, PROT_READ|PROT_WRITE,
    MAP_SHARED|MAP_POPULATE, fileno(grid->fp), (off_t)grid_offset);

  if(grid->data == MAP_FAILED) {
    htrdr_log_err(htrdr, "%s:%s: could not map the grid data -- %s.\n",
      FUNC_NAME, filename, strerror(errno));
    grid->data = NULL;
    res = RES_IO_ERR;
    goto error;
  }

  rewind(grid->fp);
  CHK(fwrite(&GRID_VERSION_NONE, sizeof(int), 1, grid->fp) == 1);

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
  uint32_t coords[3];
  uint64_t mcode;
  ASSERT(grid && xyz);
  ASSERT(xyz[0] < grid->definition[0]);
  ASSERT(xyz[1] < grid->definition[1]);
  ASSERT(xyz[2] < grid->definition[2]);
  coords[0] = (uint32_t)xyz[0];
  coords[1] = (uint32_t)xyz[1];
  coords[2] = (uint32_t)xyz[2];
  mcode = morton_xyz_encode_u21(coords);
  return htrdr_grid_at_mcode(grid, mcode);
}

void*
htrdr_grid_at_mcode(struct htrdr_grid* grid, const uint64_t mcode)
{
  ASSERT(grid);
  ASSERT(mcode < grid->def_adjusted*grid->def_adjusted*grid->def_adjusted);
  ASSERT(morton3D_decode_u21(mcode>>2) < grid->definition[0]);
  ASSERT(morton3D_decode_u21(mcode>>1) < grid->definition[1]);
  ASSERT(morton3D_decode_u21(mcode>>0) < grid->definition[2]);
  return grid->data + mcode*grid->cell_sz;
}

void
htrdr_grid_get_definition(struct htrdr_grid* grid, size_t definition[3])
{
  ASSERT(grid && definition);
  definition[0] = grid->definition[0];
  definition[1] = grid->definition[1];
  definition[2] = grid->definition[2];
}

