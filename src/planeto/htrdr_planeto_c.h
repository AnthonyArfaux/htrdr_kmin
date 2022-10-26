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

#ifndef HTRDR_PLANETO_C_H
#define HTRDR_PLANETO_C_H

#include "core/htrdr_args.h"

#include <rsys/ref_count.h>
#include <rsys/str.h>

/* Forward declarations */
struct htrdr;
struct rnatm;
struct rngrd;

struct htrdr_planeto {
  struct rnatm* atmosphere;
  struct rngrd* ground;

  struct htrdr_args_spectral spectral_domain;

  FILE* octrees_storage;

  FILE* output;
  struct str output_name;
  enum htrdr_planeto_args_output_type output_type;

  ref_T ref;
  struct htrdr* htrdr;
};

#endif /* HTRDR_PLANETO_C_H */
