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

#include "htrdr_combustion_c.h"

#include <rsys/double2.h>
#include <rsys/double3.h>

/*******************************************************************************
 * Local functions
 ******************************************************************************/
extern LOCAL_SYM double
combustion_compute_radiance_sw
  (struct htrdr_combustion* cmd,
   const size_t ithread,
   struct ssp_rng* rng,
   const double pos_in[3],
   const double dir_in[3],
   const double wlen) /* In nanometer */
{
  double pos[3];
  double dir[3];
  double range[2];
  double ksi; /* Throughput */
  double weight;
  ASSERT(cmd && rng && pos_in && dir_in);
  (void)cmd, (void)ithread, (void)rng, (void)pos_in, (void)dir_in, (void)wlen;
  (void)ksi;

  d3_set(pos, pos_in);
  d3_set(dir, dir_in);
  d2(range, 0, FLT_MAX);

  ksi = 1;
  weight = 0;

  /* TODO Radiative random walk */

  return weight;
}

