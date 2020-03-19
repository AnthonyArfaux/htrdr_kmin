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

#ifndef HTRDR_C_H
#define HTRDR_C_H

#include <rsys/rsys.h>

#ifndef NDEBUG
  #define MPI(Func) ASSERT(MPI_##Func == MPI_SUCCESS)
#else
  #define MPI(Func) MPI_##Func
#endif

enum htrdr_mpi_message {
  HTRDR_MPI_PROGRESS_BUILD_OCTREE,
  HTRDR_MPI_PROGRESS_RENDERING,
  HTRDR_MPI_STEAL_REQUEST,
  HTRDR_MPI_WORK_STEALING,
  HTRDR_MPI_TILE_DATA
};

enum htrdr_estimate {
  HTRDR_ESTIMATE_X,
  HTRDR_ESTIMATE_RADIANCE = HTRDR_ESTIMATE_X,
  HTRDR_ESTIMATE_Y,
  HTRDR_ESTIMATE_Z,
  HTRDR_ESTIMATE_TIME, /* Time per realisation */
  HTRDR_ESTIMATES_COUNT__
};

struct htrdr;

#define SW_WAVELENGTH_MIN 380 /* In nanometer */
#define SW_WAVELENGTH_MAX 780 /* In nanometer */

/* In nanometer */
static FINLINE double
wavenumber_to_wavelength(const double nu/*In cm^-1*/)
{
  return 1.e7 / nu;
}

/* In cm^-1 */
static FINLINE double
wavelength_to_wavenumber(const double lambda/*In nanometer*/)
{
  return wavenumber_to_wavelength(lambda);
}

static INLINE uint64_t
morton3D_encode_u21(const uint32_t u21)
{
  uint64_t u64 = u21 & ((1<<21) - 1);
  ASSERT(u21 <= ((1 << 21) - 1));
  u64 = (u64 | (u64 << 32)) & 0xFFFF00000000FFFF;
  u64 = (u64 | (u64 << 16)) & 0x00FF0000FF0000FF;
  u64 = (u64 | (u64 << 8))  & 0xF00F00F00F00F00F;
  u64 = (u64 | (u64 << 4))  & 0x30C30C30C30C30C3;
  u64 = (u64 | (u64 << 2))  & 0x9249249249249249;
  return u64;
}

static INLINE uint32_t
morton3D_decode_u21(const uint64_t u64)
{
  uint64_t tmp = (u64 & 0x9249249249249249);
  tmp = (tmp | (tmp >> 2))  & 0x30C30C30C30C30C3;
  tmp = (tmp | (tmp >> 4))  & 0xF00F00F00F00F00F;
  tmp = (tmp | (tmp >> 8))  & 0x00FF0000FF0000FF;
  tmp = (tmp | (tmp >> 16)) & 0xFFFF00000000FFFF;
  tmp = (tmp | (tmp >> 32)) & 0x00000000FFFFFFFF;
  ASSERT(tmp <= ((1<<21)-1));
  return (uint32_t)tmp;
}

static INLINE uint64_t
morton_xyz_encode_u21(const uint32_t xyz[3])
{
  return (morton3D_encode_u21(xyz[0]) << 2)
       | (morton3D_encode_u21(xyz[1]) << 1)
       | (morton3D_encode_u21(xyz[2]) << 0);
}

static INLINE void
morton_xyz_decode_u21(const uint64_t code, uint32_t xyz[3])
{
  ASSERT(xyz && code < ((1ull << 63)-1));
  xyz[0] = (uint32_t)morton3D_decode_u21(code >> 2);
  xyz[1] = (uint32_t)morton3D_decode_u21(code >> 1);
  xyz[2] = (uint32_t)morton3D_decode_u21(code >> 0);
}

static INLINE double
wiebelt(const double v)
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
blackbody_fraction
  (const double lambda0, /* In meter */
   const double lambda1, /* In meter */
   const double temperature) /* In Kelvin */
{
  const double C2 = 1.43877735e-2; /* m.K */
  double x0 = C2 / lambda0;
  double x1 = C2 / lambda1;
  double v0 = x0 / temperature;
  double v1 = x1 / temperature;
  double w0 = wiebelt(v0);
  double w1 = wiebelt(v1);
  return w1 - w0;
}

static INLINE double
planck
  (const double lambda_min, /* In meter */
   const double lambda_max, /* In meter */
   const double temperature) /* In Kelvin  */
{
  const double T2 = temperature*temperature;
  const double T4 = T2*T2;
  const double BOLTZMANN_CONSTANT = 5.6696e-8; /* W/m^2/K^4 */
  ASSERT(lambda_min < lambda_max && temperature > 0);
  return blackbody_fraction(lambda_min, lambda_max, temperature)
       * BOLTZMANN_CONSTANT * T4 / PI; /* In W/m^2/sr */
}

extern LOCAL_SYM res_T
brightness_temperature
  (struct htrdr* htrdr,
   const double lambda_min, /* In meter */
   const double lambda_max, /* In meter */
   /* Integrated over [lambda_min, lambda_max]. In W/m^2/sr */
   const double radiance,
   double* temperature);

extern LOCAL_SYM  res_T
open_output_stream
  (struct htrdr* htrdr,
   const char* filename,
   const int read, /* Enable read access */
   int force_overwrite,
   FILE** out_fp);

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

