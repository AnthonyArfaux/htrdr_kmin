/* Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
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

#ifndef HTRDR_C_H
#define HTRDR_C_H

#include <rsys/logger.h>
#include <rsys/ref_count.h>
#include <rsys/rsys.h>

#ifndef NDEBUG
  #define MPI(Func) ASSERT(MPI_##Func == MPI_SUCCESS)
#else
  #define MPI(Func) MPI_##Func
#endif

enum htrdr_mpi_message {
  HTRDR_MPI_PROGRESS_RENDERING,
  HTRDR_MPI_STEAL_REQUEST,
  HTRDR_MPI_WORK_STEALING,
  HTRDR_MPI_TILE_DATA
};

struct htrdr {
  unsigned nthreads; /* #threads of the process */

  int mpi_rank; /* Rank of the process in the MPI group */
  int mpi_nprocs; /* Overall #processes in the MPI group */
  char* mpi_err_str; /* Temp buffer used to store MPI error string */
  int8_t* mpi_working_procs; /* Define the rank of active processes */
  size_t mpi_nworking_procs;

  /* Process progress percentage */
  int32_t* mpi_progress_octree;
  int32_t* mpi_progress_render;

  struct mutex* mpi_mutex; /* Protect MPI calls from concurrent threads */

  int verbose;

  struct logger logger;
  struct mem_allocator* allocator;
  struct mem_allocator* lifo_allocators; /* Per thread lifo allocator */

  ref_T ref;
};

/* In nanometer */
static FINLINE double
wavenumber_to_wavelength(const double nu/*In cm^-1*/)
{
  return 1.e7 / nu;
}

/* In cm^-1 */
static FINLINE double
wavelength_to_wavenumber(const double lambda/*In nanometer*/)
{
  return wavenumber_to_wavelength(lambda);
}

static INLINE uint64_t
morton3D_encode_u21(const uint32_t u21)
{
  uint64_t u64 = u21 & ((1<<21) - 1);
  ASSERT(u21 <= ((1 << 21) - 1));
  u64 = (u64 | (u64 << 32)) & 0xFFFF00000000FFFF;
  u64 = (u64 | (u64 << 16)) & 0x00FF0000FF0000FF;
  u64 = (u64 | (u64 << 8))  & 0xF00F00F00F00F00F;
  u64 = (u64 | (u64 << 4))  & 0x30C30C30C30C30C3;
  u64 = (u64 | (u64 << 2))  & 0x9249249249249249;
  return u64;
}

static INLINE uint32_t
morton3D_decode_u21(const uint64_t u64)
{
  uint64_t tmp = (u64 & 0x9249249249249249);
  tmp = (tmp | (tmp >> 2))  & 0x30C30C30C30C30C3;
  tmp = (tmp | (tmp >> 4))  & 0xF00F00F00F00F00F;
  tmp = (tmp | (tmp >> 8))  & 0x00FF0000FF0000FF;
  tmp = (tmp | (tmp >> 16)) & 0xFFFF00000000FFFF;
  tmp = (tmp | (tmp >> 32)) & 0x00000000FFFFFFFF;
  ASSERT(tmp <= ((1<<21)-1));
  return (uint32_t)tmp;
}

static INLINE uint64_t
morton_xyz_encode_u21(const uint32_t xyz[3])
{
  return (morton3D_encode_u21(xyz[0]) << 2)
       | (morton3D_encode_u21(xyz[1]) << 1)
       | (morton3D_encode_u21(xyz[2]) << 0);
}

static INLINE void
morton_xyz_decode_u21(const uint64_t code, uint32_t xyz[3])
{
  ASSERT(xyz && code < ((1ull << 63)-1));
  xyz[0] = (uint32_t)morton3D_decode_u21(code >> 2);
  xyz[1] = (uint32_t)morton3D_decode_u21(code >> 1);
  xyz[2] = (uint32_t)morton3D_decode_u21(code >> 0);
}

/* Return the minimum length in nanometer of the sky spectral bands
 * clamped to in [range[0], range[1]]. */
extern LOCAL_SYM double
compute_sky_min_band_len
  (struct htsky* sky,
   const double range[2]);

extern LOCAL_SYM  res_T
open_output_stream
  (struct htrdr* htrdr,
   const char* filename,
   const int read, /* Enable read access */
   int force_overwrite,
   FILE** out_fp);

extern LOCAL_SYM void
send_mpi_progress
  (struct htrdr* htrdr,
   const enum htrdr_mpi_message progress,
   const int32_t percent);

extern LOCAL_SYM void
fetch_mpi_progress
  (struct htrdr* htrdr,
   const enum htrdr_mpi_message progress);

extern LOCAL_SYM void
print_mpi_progress
  (struct htrdr* htrdr,
   const enum htrdr_mpi_message progress);

extern LOCAL_SYM void
clear_mpi_progress
  (struct htrdr* htrdr,
   const enum htrdr_mpi_message progress);

extern int32_t
total_mpi_progress
  (const struct htrdr* htrdr,
   const enum htrdr_mpi_message progress);

static INLINE void
update_mpi_progress(struct htrdr* htrdr, const enum htrdr_mpi_message progress)
{
  ASSERT(htrdr);
  fetch_mpi_progress(htrdr, progress);
  clear_mpi_progress(htrdr, progress);
  print_mpi_progress(htrdr, progress);
}

static FINLINE int
cmp_dbl(const void* a, const void* b)
{
  const double d0 = *((const double*)a);
  const double d1 = *((const double*)b);
  return d0 < d1 ? -1 : (d0 > d1 ? 1 : 0);
}

#endif /* HTRDR_C_H */

