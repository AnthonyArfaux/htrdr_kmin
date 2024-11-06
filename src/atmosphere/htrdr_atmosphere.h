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

#ifndef HTRDR_ATMOSPHERE_H
#define HTRDR_ATMOSPHERE_H

#include "core/htrdr.h"
#include <rsys/rsys.h>

/* Forward declarations */
struct htrdr;
struct htrdr_atmosphere;
struct htrdr_atmosphere_args;

BEGIN_DECLS

HTRDR_API res_T
htrdr_atmosphere_create
  (struct htrdr* htrdr,
   const struct htrdr_atmosphere_args* args,
   struct htrdr_atmosphere** cmd);

HTRDR_API void
htrdr_atmosphere_ref_get
  (struct htrdr_atmosphere* cmd);

HTRDR_API void
htrdr_atmosphere_ref_put
  (struct htrdr_atmosphere* cmd);

HTRDR_API res_T
htrdr_atmosphere_run
  (struct htrdr_atmosphere* cmd);

HTRDR_API int
htrdr_atmosphere_main
  (int argc,
   char** argv);

END_DECLS

#endif /* HTRDR_ATMOSPHERE_H */

