/* Copyright (C) 2018-2019, 2022-2024 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2024 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2024 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2024 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2024 Observatoire de Paris
 * Copyright (C) 2022-2024 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2024 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2024 Université Paul Sabatier
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

struct s3d_device;

struct htrdr {
  struct s3d_device* s3d;

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

extern LOCAL_SYM void
setup_logger
  (struct htrdr* htrdr);

/* Return the minimum length in nanometer of the sky spectral bands
 * clamped to in [range[0], range[1]]. */
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

