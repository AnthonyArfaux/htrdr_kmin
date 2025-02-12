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

#include "core/htrdr.h"
#include "core/htrdr_c.h"
#include "core/htrdr_log.h"

#include <rsys/logger.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
print_out(const char* msg, void* ctx)
{
  struct htrdr* htrdr = ctx;
  ASSERT(msg);

  /* Log standard message only on master process */
  if(htrdr->verbose && htrdr->mpi_rank == 0) {
    fprintf(stderr, HTRDR_LOG_INFO_PREFIX"%s", msg);
  }
}

static void
print_err(const char* msg, void* ctx)
{
  struct htrdr* htrdr = ctx;
  ASSERT(msg);

  /* Log errors on all processes */
  if(htrdr->verbose) {
    fprintf(stderr, HTRDR_LOG_ERROR_PREFIX"%s", msg);
  }
}

static void
print_warn(const char* msg, void* ctx)
{
  struct htrdr* htrdr = ctx;
  ASSERT(msg && ctx);

  /* Log warnings only on master process */
  if(htrdr->verbose && htrdr->mpi_rank == 0) {
    fprintf(stderr, HTRDR_LOG_WARNING_PREFIX"%s", msg);
  }
}

static void
log_msg
  (struct htrdr* htrdr,
   const enum log_type stream,
   const char* msg,
   va_list vargs)
{
  ASSERT(htrdr && msg);
  CHK(logger_vprint(&htrdr->logger, stream, msg, vargs) == RES_OK);
}

/*******************************************************************************
 * Exported functions
 ******************************************************************************/
void
htrdr_log(struct htrdr* htrdr, const char* msg, ...)
{
  va_list vargs_list;
  ASSERT(htrdr && msg);
  va_start(vargs_list, msg);
  log_msg(htrdr, LOG_OUTPUT, msg, vargs_list);
  va_end(vargs_list);
}

void
htrdr_log_err(struct htrdr* htrdr, const char* msg, ...)
{
  va_list vargs_list;
  ASSERT(htrdr && msg);
  va_start(vargs_list, msg);
  log_msg(htrdr, LOG_ERROR, msg, vargs_list);
  va_end(vargs_list);
}

void
htrdr_log_warn(struct htrdr* htrdr, const char* msg, ...)
{
  va_list vargs_list;
  ASSERT(htrdr && msg);
  va_start(vargs_list, msg);
  log_msg(htrdr, LOG_WARNING, msg, vargs_list);
  va_end(vargs_list);
}

/*******************************************************************************
 * Local function
 ******************************************************************************/
void
setup_logger(struct htrdr* htrdr)
{
  logger_init(htrdr->allocator, &htrdr->logger);
  logger_set_stream(&htrdr->logger, LOG_OUTPUT, print_out, htrdr);
  logger_set_stream(&htrdr->logger, LOG_ERROR, print_err, htrdr);
  logger_set_stream(&htrdr->logger, LOG_WARNING, print_warn, htrdr);
}
