/* Copyright (C) 2018-2019, 2022-2025 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2025 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2025 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2025 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2025 Observatoire de Paris
 * Copyright (C) 2022-2025 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2025 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2025 Université Paul Sabatier
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

#ifndef HTRDR_PROC_WORK_H
#define HTRDR_PROC_WORK_H

#include <rsys/dynamic_array_u64.h>

#include <mpi.h>

#define CHUNK_ID_NULL UINT64_MAX

/* Forward declarations */
struct mutex;
struct ssp_rng;

/* List of chunks to compute onto the MPI process */
struct proc_work {
  struct mutex* mutex;
  struct darray_u64 chunks; /* #chunks to solve */
  uint64_t index; /* Next chunk to solve in the above list of chunks */
};

extern LOCAL_SYM void
proc_work_init
  (struct mem_allocator* allocator,
   struct proc_work* work);

extern LOCAL_SYM void
proc_work_release
  (struct proc_work* work);

extern LOCAL_SYM void
proc_work_reset
  (struct proc_work* work);

extern LOCAL_SYM void
proc_work_add_chunk
  (struct proc_work* work,
   const size_t ichunk);

/* Return the index of the next chunk to be processed */
extern LOCAL_SYM uint64_t
proc_work_get_chunk
  (struct proc_work* work);

extern LOCAL_SYM uint64_t
proc_work_get_nchunks
  (struct proc_work* work);

/* Wait for the completion of an MPI request */
extern LOCAL_SYM void
mpi_wait_for_request
  (struct htrdr* htrdr,
   MPI_Request* req);

/* Active polling of the "steal" request submitted by other processes to
 * relieve the current process of the chunks assigned to it.
 * The function runs until probe_thieves is set to 0 by the caller. */
extern LOCAL_SYM void
mpi_probe_thieves
  (struct htrdr* htrdr,
   struct proc_work* work,
   ATOMIC* probe_thieves);

/* Unload a working process by taking chunks from it,
 * i.e. submit a steal request and wait for it to be honored */
extern LOCAL_SYM size_t
mpi_steal_work
  (struct htrdr* htrdr,
   struct ssp_rng* rng,
   struct proc_work* work);

#endif /* HTRDR_PROC_WORK_H */
