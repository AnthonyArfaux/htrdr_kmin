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

#include "htrdr.h"
#include "htrdr_spectral.h"

res_T
brightness_temperature
  (struct htrdr* htrdr,
   const double lambda_min,
   const double lambda_max,
   const double radiance, /* In W/m2/sr/m */
   double* temperature)
{
  const size_t MAX_ITER = 100;
  const double epsilon_T = 1e-4; /* In K */
  const double epsilon_B = radiance * 1e-8;
  double T, T0, T1, T2;
  double B, B0;
  size_t i;
  res_T res = RES_OK;
  ASSERT(temperature && lambda_min <= lambda_max);

  /* Search for a brightness temperature whose radiance is greater than or
   * equal to the estimated radiance */
  T2 = 200;
  FOR_EACH(i, 0, MAX_ITER) {
    const double B2 = planck(lambda_min, lambda_max, T2);
    if(B2 >= radiance) break;
    T2 *= 2;
  }
  if(i >= MAX_ITER) { res = RES_BAD_OP; goto error; }

  B0 = T0 = T1 = 0;
  FOR_EACH(i, 0, MAX_ITER) {
    T = (T1+T2)*0.5;
    B = planck(lambda_min, lambda_max, T);

    if(B < radiance) {
      T1 = T;
    } else {
      T2 = T;
    }

    if(fabs(T-T0) < epsilon_T || fabs(B-B0) < epsilon_B)
      break;

    T0 = T;
    B0 = B;
  }
  if(i >= MAX_ITER) { res = RES_BAD_OP; goto error; }

  *temperature = T;

exit:
  return res;
error:
  htrdr_log_err(htrdr,
    "Could not compute the brightness temperature for the estimated radiance %g "
    "averaged over [%g, %g] nanometers.\n",
    radiance,
    lambda_min*1e9,
    lambda_max*1e9);
  goto exit;
}

