/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
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

#ifndef HTRDR_SPECTRAL_H
#define HTRDR_SPECTRAL_H

#include "core/htrdr.h"

#include <rsys/rsys.h>
#include <rsys/math.h> /* PI support */

#define HTRDR_SUN_TEMPERATURE 5778 /* In K */
#define HTRDR_DEFAULT_LW_REF_TEMPERATURE 290 /* Default longwave temperature in K */

enum htrdr_spectral_type {
  HTRDR_SPECTRAL_LW, /* Longwave */
  HTRDR_SPECTRAL_SW, /* Shortwave */
  HTRDR_SPECTRAL_SW_CIE_XYZ /* Shortwave wrt the CIE XYZ tristimulus */
};

/* Forwar declaration */
struct htrdr;

static FINLINE double /* In nanometer */
htrdr_wavenumber_to_wavelength(const double nu/*In cm⁻¹*/)
{
  return 1.e7 / nu;
}

static FINLINE double /* In cm⁻¹ */
wavelength_to_wavenumber(const double lambda/*In nanometer*/)
{
  return htrdr_wavenumber_to_wavelength(lambda);
}

static INLINE double
htrdr_wiebelt(const double v)
{
  int m;
  double w, v2, v4;
  /*.153989717364e+00;*/
  const double fifteen_over_pi_power_4 = 15.0/(PI*PI*PI*PI);
  const double z0 = 1.0/3.0;
  const double z1 = 1.0/8.0;
  const double z2 = 1.0/60.0;
  const double z4 = 1.0/5040.0;
  const double z6 = 1.0/272160.0;
  const double z8 = 1.0/13305600.0;

  if(v >= 2.) {
    w = 0;
    for(m=1; m<6 ;m++)
      w+=exp(-m*v)/(m*m*m*m) * (((m*v+3)*m*v+6)*m*v+6);
    w = w * fifteen_over_pi_power_4;
  } else {
    v2 = v*v;
    v4 = v2*v2;
    w = z0 - z1*v + z2*v2 - z4*v2*v2 + z6*v4*v2 - z8*v4*v4;
    w = 1. - fifteen_over_pi_power_4*v2*v*w;
  }
  ASSERT(w >= 0.0 && w <= 1.0);
  return w;
}

static INLINE double
htrdr_blackbody_fraction
  (const double lambda0, /* In meters */
   const double lambda1, /* In meters */
   const double temperature) /* In Kelvin */
{
  const double C2 = 1.43877735e-2; /* m.K */
  double x0 = C2 / lambda0;
  double x1 = C2 / lambda1;
  double v0 = x0 / temperature;
  double v1 = x1 / temperature;
  double w0 = htrdr_wiebelt(v0);
  double w1 = htrdr_wiebelt(v1);
  return w1 - w0;
}

/* Return the Planck value in W/m²/sr/m at a given wavelength */
static INLINE double
htrdr_planck_monochromatic
  (const double lambda, /* In meters */
   const double temperature) /* In Kelvin */
{
  const double c = 299792458; /* m/s */
  const double h = 6.62607015e-34; /* J.s */
  const double k = 1.380649e-23; /* J/K */
  const double lambda2 = lambda*lambda;
  const double lambda5 = lambda2*lambda2*lambda;
  const double B = ((2.0 * h * c*c) / lambda5) /* W/m²/sr/m */
                 / (exp(h*c/(lambda*k*temperature))-1.0);
  ASSERT(temperature > 0);
  return B;
}

/* Return the average Planck value in W/m²/sr/m over the
 * [lambda_min, lambda_max] interval. */
static INLINE double
htrdr_planck_interval
  (const double lambda_min, /* In meters */
   const double lambda_max, /* In meters */
   const double temperature) /* In Kelvin  */
{
  const double T2 = temperature*temperature;
  const double T4 = T2*T2;
  const double BOLTZMANN_CONSTANT = 5.6696e-8; /* W/m²/K⁴ */
  ASSERT(lambda_min < lambda_max && temperature > 0);
  return htrdr_blackbody_fraction(lambda_min, lambda_max, temperature)
       * BOLTZMANN_CONSTANT * T4 / (PI * (lambda_max-lambda_min)); /* In W/m²/sr/m */
}

/* Invoke planck_monochromatic or planck_interval whether the submitted
 * interval is null or not, respectively. The returned value is in W/m²/sr/m */
static FINLINE double
htrdr_planck
  (const double lambda_min, /* In meters */
   const double lambda_max, /* In meters */
   const double temperature) /* In Kelvin  */
{
  ASSERT(lambda_min <= lambda_max && temperature > 0);
  if(lambda_min == lambda_max) {
    return htrdr_planck_monochromatic(lambda_min, temperature);
  } else {
    return htrdr_planck_interval(lambda_min, lambda_max, temperature);
  }
}

BEGIN_DECLS

HTRDR_CORE_API res_T
htrdr_brightness_temperature
  (struct htrdr* htrdr,
   const double lambda_min, /* In meters */
   const double lambda_max, /* In meters */
   /* Averaged over [lambda_min, lambda_max]. In W/m²/sr/m */
   const double radiance,
   double* temperature);

HTRDR_CORE_API double
htrdr_radiance_temperature
  (struct htrdr* htrdr,
   const double lambda_min, /* In meters */
   const double lambda_max, /* In meters */
   const double radiance); /* In W/m²/sr */

END_DECLS

#endif /* HTRDR_SPECTRAL_H */
