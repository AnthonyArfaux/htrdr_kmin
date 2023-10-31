/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2023 Observatoire de Paris
 * Copyright (C) 2022-2023 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2023 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2023 Université Paul Sabatier
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

#ifndef HTRDR_LOG_H
#define HTRDR_LOG_H

#include "core/htrdr.h"
#include <rsys/rsys.h>

#define HTRDR_LOG_INFO_PREFIX "\x1b[1m\x1b[32m>\x1b[0m "
#define HTRDR_LOG_ERROR_PREFIX "\x1b[31merror:\x1b[0m "
#define HTRDR_LOG_WARNING_PREFIX "\x1b[33mwarning:\x1b[0m "

struct htrdr;

BEGIN_DECLS

HTRDR_API void
htrdr_log
  (struct htrdr* htrdr,
   const char* msg,
   ...)
#ifdef COMPILER_GCC
  __attribute((format(printf, 2, 3)))
#endif
  ;

HTRDR_API void
htrdr_log_err
  (struct htrdr* htrdr,
   const char* msg,
   ...)
#ifdef COMPILER_GCC
  __attribute((format(printf, 2, 3)))
#endif
  ;

HTRDR_API void
htrdr_log_warn
  (struct htrdr* htrdr,
   const char* msg,
   ...)
#ifdef COMPILER_GCC
  __attribute((format(printf, 2, 3)))
#endif
  ;

END_DECLS

#endif /* HTRDR_LOG_H */
