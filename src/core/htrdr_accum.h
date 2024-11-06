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

#ifndef HTRDR_ACCUM_H
#define HTRDR_ACCUM_H

#include <rsys/rsys.h>
#include <rsys/math.h>

/* Monte carlo accumulator */
struct htrdr_accum {
  double sum_weights; /* Sum of Monte-Carlo weights */
  double sum_weights_sqr; /* Sum of Monte-Carlo square weights */
  size_t nweights; /* #accumlated weights */
  size_t nfailures; /* #failures */
};
#define HTRDR_ACCUM_NULL__ {0, 0, 0, 0}
static const struct htrdr_accum HTRDR_ACCUM_NULL = HTRDR_ACCUM_NULL__;

/* Monte carlo estimate */
struct htrdr_estimate {
  double E; /* Expected value */
  double SE; /* Standard error */
};
#define HTRDR_ESTIMATE_NULL__ {0, 0}
static const struct htrdr_estimate HTRDR_ESTIMATE_NULL = HTRDR_ESTIMATE_NULL__;

static FINLINE void
htrdr_accum_get_estimation
  (const struct htrdr_accum* acc,
   struct htrdr_estimate* estimate)
{
  ASSERT(acc && estimate);

  if(!acc->nweights) {
    estimate->E = 0;
    estimate->SE = 0;
  } else {
    const double N = (double)acc->nweights;
    double E, V, SE;
    E = acc->sum_weights / N;
    V = MMAX(acc->sum_weights_sqr / N - E*E, 0);
    SE = sqrt(V/N);

    estimate->E = E;
    estimate->SE = SE;
  }
}

#endif /* HTRDR_SOLVE_H */
