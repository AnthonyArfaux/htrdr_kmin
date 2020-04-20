/* Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
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

#define _POSIX_C_SOURCE 200112L /* nextafter */

#include "htrdr.h"
#include "htrdr_c.h"
#include "htrdr_cie_xyz.h"

#include <rsys/algorithm.h>
#include <rsys/cstr.h>
#include <rsys/dynamic_array_double.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

#include <math.h> /* nextafter */

struct htrdr_cie_xyz {
  struct darray_double cdf_X;
  struct darray_double cdf_Y;
  struct darray_double cdf_Z;
  double range[2]; /* Boundaries of the handled CIE XYZ color space */
  double band_len; /* Length in nanometers of a band */

  struct htrdr* htrdr;
  ref_T ref;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static INLINE double
trapezoidal_integration
  (const double lambda_lo, /* Integral lower bound. In nanometer */
   const double lambda_hi, /* Integral upper bound. In nanometer */
   double (*f_bar)(const double lambda)) /* Function to integrate */
{
  double dlambda;
  size_t i, n;
  double integral = 0;
  ASSERT(lambda_lo <= lambda_hi);
  ASSERT(lambda_lo > 0);

  n = (size_t)(lambda_hi - lambda_lo) + 1;
  dlambda = (lambda_hi - lambda_lo) / (double)n;

  FOR_EACH(i, 0, n) {
    const double lambda1 = lambda_lo + dlambda*(double)(i+0);
    const double lambda2 = lambda_hi + dlambda*(double)(i+1);
    const double f1 = f_bar(lambda1);
    const double f2 = f_bar(lambda2);
    integral += (f1 + f2)*dlambda*0.5;
  }
  return integral;
}

/* The following 3 functions are used to fit the CIE Xbar, Ybar and Zbar curved
 * has defined by the 1931 standard. These analytical fits are propsed by C.
 * Wyman, P. P. Sloan & P. Shirley in "Simple Analytic Approximations to the
 * CIE XYZ Color Matching Functions" - JCGT 2013. */
static INLINE double
fit_x_bar_1931(const double lambda)
{
  const double a = (lambda - 442.0) * (lambda < 442.0 ? 0.0624 : 0.0374);
  const double b = (lambda - 599.8) * (lambda < 599.8 ? 0.0264 : 0.0323);
  const double c = (lambda - 501.1) * (lambda < 501.1 ? 0.0490 : 0.0382);
  return 0.362*exp(-0.5*a*a) + 1.056*exp(-0.5f*b*b) - 0.065*exp(-0.5*c*c);
}

static FINLINE double
fit_y_bar_1931(const double lambda)
{
  const double a = (lambda - 568.8) * (lambda < 568.8 ? 0.0213 : 0.0247);
  const double b = (lambda - 530.9) * (lambda < 530.9 ? 0.0613 : 0.0322);
  return 0.821*exp(-0.5*a*a) + 0.286*exp(-0.5*b*b);
}

static FINLINE double
fit_z_bar_1931(const double lambda)
{
  const double a = (lambda - 437.0) * (lambda < 437.0 ? 0.0845 : 0.0278);
  const double b = (lambda - 459.0) * (lambda < 459.0 ? 0.0385 : 0.0725);
  return 1.217*exp(-0.5*a*a) + 0.681*exp(-0.5*b*b);
}

static INLINE double
sample_cie_xyz
  (const struct htrdr_cie_xyz* cie,
   const double* cdf,
   const size_t cdf_length,
   double (*f_bar)(const double lambda), /* Function to integrate */
   const double r0, /* Canonical number in [0, 1[ */
   const double r1) /* Canonical number in [0, 1[ */
{
  double r0_next = nextafter(r0, DBL_MAX);
  double* find;
  double f_min, f_max; /* CIE 1931 value for the band boundaries */
  double lambda; /* Sampled wavelength */
  double lambda_min, lambda_max; /* Boundaries of the sampled band */
  double lambda_1, lambda_2; /* Solutions if the equation to solve */
  double a, b, c, d; /* Equation parameters */
  double delta, sqrt_delta;
  size_t iband; /* Index of the sampled band */
  ASSERT(cie && cdf && cdf_length);
  ASSERT(0 <= r0 && r0 < 1);
  ASSERT(0 <= r1 && r1 < 1);

  /* Use r_next rather than r in order to find the first entry that is not less
   * than *or equal* to r */
  find = search_lower_bound(&r0_next, cdf, cdf_length, sizeof(double), cmp_dbl);
  ASSERT(find);

  /* Define and check the sampled band */
  iband = (size_t)(find - cdf);
  ASSERT(iband < cdf_length);
  ASSERT(cdf[iband] > r0 && (!iband || cdf[iband-1] <= r0));

  /* Define the boundaries of the sampled band */
  lambda_min = cie->range[0] + cie->band_len * (double)iband;
  lambda_max = lambda_min + cie->band_len;

  /* Define the value of the CIE 1931 function for the boudaries of the sampled
   * band */
  f_min = f_bar(lambda_min);
  f_max = f_bar(lambda_max);

  /* Compute the equation constants */
  a = 0.5 * (f_max - f_min) / cie->band_len;
  b = (lambda_max * f_min - lambda_min * f_max) / cie->band_len;
  c = -lambda_min * f_min + lambda_min*lambda_min * a;
  d = 0.5 * (f_max + f_min) * cie->band_len;

  delta = b*b - 4*a*(c-d*r1);
  ASSERT(delta > 0);
  sqrt_delta = sqrt(delta);

  /* Compute the roots that solve the equation */
  lambda_1 = (-b - sqrt_delta) / (2*a);
  lambda_2 = (-b + sqrt_delta) / (2*a);

  /* Select the solution */
  if(lambda_min <= lambda_1 && lambda_1 < lambda_max) {
    lambda = lambda_1;
  } else if(lambda_min <= lambda_2 && lambda_2 < lambda_max) {
    lambda = lambda_2;
  } else {
    FATAL("Unexpected error.\n");
  }

  return lambda;
}

static res_T
setup_cie_xyz
  (struct htrdr_cie_xyz* cie,
   const char* func_name,
   const double range[2],
   const size_t nbands)
{
  enum { X, Y, Z }; /* Helper constant */
  double* pdf[3] = {NULL, NULL, NULL};
  double* cdf[3] = {NULL, NULL, NULL};
  double sum[3] = {0,0,0};
  size_t i;
  res_T res = RES_OK;

  ASSERT(cie && func_name && range && nbands);
  ASSERT(range[0] >= HTRDR_CIE_XYZ_RANGE_DEFAULT[0]);
  ASSERT(range[1] <= HTRDR_CIE_XYZ_RANGE_DEFAULT[1]);
  ASSERT(range[0] < range[1]);

  /* Allocate and reset the memory space for the tristimulus CDF */
  #define SETUP_STIMULUS(Stimulus) {                                           \
    res = darray_double_resize(&cie->cdf_ ## Stimulus, nbands);                \
    if(res != RES_OK) {                                                        \
      htrdr_log_err(cie->htrdr,                                                \
        "%s: Could not reserve the memory space for the CDF "                  \
        "of the "STR(X)" stimulus -- %s.\n", func_name, res_to_cstr(res));     \
      goto error;                                                              \
    }                                                                          \
    cdf[Stimulus] = darray_double_data_get(&cie->cdf_ ## Stimulus);            \
    pdf[Stimulus] = cdf[Stimulus];                                             \
    memset(cdf[Stimulus], 0, nbands*sizeof(double));                           \
  } (void)0
  SETUP_STIMULUS(X);
  SETUP_STIMULUS(Y);
  SETUP_STIMULUS(Z);
  #undef SETUP_STIMULUS

  /* Compute the *unormalized* pdf of the tristimulus */
  cie->range[0] = range[0];
  cie->range[1] = range[1];
  cie->band_len = (range[1] - range[0]) / (double)nbands;
  FOR_EACH(i, 0, nbands) {
    const double lambda_lo = range[0] + (double)i * cie->band_len;
    const double lambda_hi = lambda_lo + cie->band_len;
    ASSERT(lambda_lo <= lambda_hi);
    pdf[X][i] = trapezoidal_integration(lambda_lo, lambda_hi, fit_x_bar_1931);
    pdf[Y][i] = trapezoidal_integration(lambda_lo, lambda_hi, fit_y_bar_1931);
    pdf[Z][i] = trapezoidal_integration(lambda_lo, lambda_hi, fit_z_bar_1931);
    sum[X] += pdf[X][i];
    sum[Y] += pdf[Y][i];
    sum[Z] += pdf[Z][i];
  }

  FOR_EACH(i, 0, nbands) {
    /* Normalize the pdf */
    pdf[X][i] /= sum[X];
    pdf[Y][i] /= sum[Y];
    pdf[Z][i] /= sum[Z];
    /* Setup the cumulative */
    if(i == 0) {
      cdf[X][i] = pdf[X][i];
      cdf[Y][i] = pdf[Y][i];
      cdf[Z][i] = pdf[Z][i];
    } else {
      cdf[X][i] = pdf[X][i] + cdf[X][i-1];
      cdf[Y][i] = pdf[Y][i] + cdf[Y][i-1];
      cdf[Z][i] = pdf[Z][i] + cdf[Z][i-1];
      ASSERT(cdf[X][i] >= cdf[X][i-1]);
      ASSERT(cdf[Y][i] >= cdf[Y][i-1]);
      ASSERT(cdf[Z][i] >= cdf[Z][i-1]);
    }
  }
  ASSERT(eq_eps(cdf[X][nbands-1], 1, 1.e-6));
  ASSERT(eq_eps(cdf[Y][nbands-1], 1, 1.e-6));
  ASSERT(eq_eps(cdf[Z][nbands-1], 1, 1.e-6));

  /* Handle numerical issue */
  cdf[X][nbands-1] = 1.0;
  cdf[Y][nbands-1] = 1.0;
  cdf[Z][nbands-1] = 1.0;

exit:
  return res;
error:
  darray_double_clear(&cie->cdf_X);
  darray_double_clear(&cie->cdf_Y);
  darray_double_clear(&cie->cdf_Z);
  goto exit;
}

static void
release_cie_xyz(ref_T* ref)
{
  struct htrdr_cie_xyz* cie = NULL;
  ASSERT(ref);
  cie = CONTAINER_OF(ref, struct htrdr_cie_xyz, ref);
  darray_double_release(&cie->cdf_X);
  darray_double_release(&cie->cdf_Y);
  darray_double_release(&cie->cdf_Z);
  MEM_RM(cie->htrdr->allocator, cie);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_cie_xyz_create
  (struct htrdr* htrdr,
   const double range[2], /* Must be included in [380, 780] nanometers */
   const size_t nbands, /* # bands used to discretisze the CIE tristimulus */
   struct htrdr_cie_xyz** out_cie)
{
  struct htrdr_cie_xyz* cie = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && range && nbands && out_cie);

  cie = MEM_CALLOC(htrdr->allocator, 1, sizeof(*cie));
  if(!cie) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "%s: could not allocate the CIE XYZ data structure -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }
  ref_init(&cie->ref);
  cie->htrdr = htrdr;
  darray_double_init(htrdr->allocator, &cie->cdf_X);
  darray_double_init(htrdr->allocator, &cie->cdf_Y);
  darray_double_init(htrdr->allocator, &cie->cdf_Z);

  res = setup_cie_xyz(cie, FUNC_NAME, range, nbands);
  if(res != RES_OK) goto error;

exit:
  *out_cie = cie;
  return res;
error:
  if(cie) htrdr_cie_xyz_ref_put(cie);
  goto exit;
}

void
htrdr_cie_xyz_ref_get(struct htrdr_cie_xyz* cie)
{
  ASSERT(cie);
  ref_get(&cie->ref);
}

void
htrdr_cie_xyz_ref_put(struct htrdr_cie_xyz* cie)
{
  ASSERT(cie);
  ref_put(&cie->ref, release_cie_xyz);
}

double
htrdr_cie_xyz_sample_X
  (struct htrdr_cie_xyz* cie,
   const double r0,
   const double r1)
{
  return sample_cie_xyz(cie, darray_double_cdata_get(&cie->cdf_X),
    darray_double_size_get(&cie->cdf_X), fit_x_bar_1931, r0, r1);
}

double
htrdr_cie_xyz_sample_Y
  (struct htrdr_cie_xyz* cie,
   const double r0,
   const double r1)
{
  return sample_cie_xyz(cie, darray_double_cdata_get(&cie->cdf_Y),
    darray_double_size_get(&cie->cdf_Y), fit_y_bar_1931, r0, r1);
}

double
htrdr_cie_xyz_sample_Z
  (struct htrdr_cie_xyz* cie,
   const double r0,
   const double r1)
{
  return sample_cie_xyz(cie, darray_double_cdata_get(&cie->cdf_Z),
    darray_double_size_get(&cie->cdf_Z), fit_z_bar_1931, r0, r1);
}

