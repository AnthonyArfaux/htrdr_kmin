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
#include "htrdr_ran_lw.h"

#include <high_tune/htsky.h>

#include <rsys/algorithm.h>
#include <rsys/cstr.h>
#include <rsys/dynamic_array_double.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

#include <math.h> /* nextafter */

struct htrdr_ran_lw {
  struct darray_double pdf;
  struct darray_double cdf;
  double range[2]; /* Boundaries of the spectral integration interval */
  double band_len; /* Length in nanometers of a band */
  double ref_temperature; /* In Kelvin */
  size_t nbands; /* # bands */

  struct htrdr* htrdr;
  ref_T ref;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static res_T
setup_ran_lw_cdf
  (struct htrdr_ran_lw* ran_lw,
   const char* func_name)
{
  double* pdf = NULL;
  double* cdf = NULL;
  double sum = 0;
  size_t i;
  res_T res = RES_OK;
  ASSERT(ran_lw && func_name && ran_lw->nbands != HTRDR_RAN_LW_CONTINUE);

  res = darray_double_resize(&ran_lw->cdf, ran_lw->nbands);
  if(res != RES_OK) {
    htrdr_log_err(ran_lw->htrdr,
      "%s: Error allocating the CDF of the long wave spectral bands -- %s.\n",
      func_name, res_to_cstr(res));
    goto error;
  }
  res = darray_double_resize(&ran_lw->pdf, ran_lw->nbands);
  if(res != RES_OK) {
    htrdr_log_err(ran_lw->htrdr,
      "%s: Error allocating the pdf of the long wave spectral bands -- %s.\n",
      func_name, res_to_cstr(res));
    goto error;
  }

  cdf = darray_double_data_get(&ran_lw->cdf);
  pdf = darray_double_data_get(&ran_lw->pdf); /* Now save pdf to correct weight */

  htrdr_log(ran_lw->htrdr,
    "Number of bands of the long wave cumulative: %lu\n",
    (unsigned long)ran_lw->nbands);

  /* Compute the *unormalized* probability to sample a long wave band */
  FOR_EACH(i, 0, ran_lw->nbands) {
    double lambda_lo = ran_lw->range[0] + (double)i * ran_lw->band_len;
    double lambda_hi = MMIN(lambda_lo + ran_lw->band_len, ran_lw->range[1]);
    ASSERT(lambda_lo<= lambda_hi);
    ASSERT(lambda_lo > ran_lw->range[0] || eq_eps(lambda_lo, ran_lw->range[0], 1.e-6));
    ASSERT(lambda_lo < ran_lw->range[1] || eq_eps(lambda_lo, ran_lw->range[1], 1.e-6));

    /* Convert from nanometer to meter */
    lambda_lo *= 1.e-9;
    lambda_hi *= 1.e-9;

    /* Compute the probability of the current band */
    pdf[i] = blackbody_fraction(lambda_lo, lambda_hi, ran_lw->ref_temperature);

    /* Update the norm */
    sum += pdf[i];
  }

  /* Compute the cumulative of the previously computed probabilities */
  FOR_EACH(i, 0, ran_lw->nbands) {
    /* Normalize the probability */
    pdf[i] /= sum;

    /* Setup the cumulative */
    if(i == 0) {
      cdf[i] = pdf[i];
    } else {
      cdf[i] = pdf[i] + cdf[i-1];
      ASSERT(cdf[i] >= cdf[i-1]);
    }
  }
  ASSERT(eq_eps(cdf[ran_lw->nbands-1], 1, 1.e-6));
  cdf[ran_lw->nbands - 1] = 1.0; /* Handle numerical issue */

exit:
  return res;
error:
  darray_double_clear(&ran_lw->cdf);
  darray_double_clear(&ran_lw->pdf);
  goto exit;
}

static double
ran_lw_sample_continue
  (const struct htrdr_ran_lw* ran_lw,
   const double r,
   const double range[2], /* In nanometer */
   const char* func_name,
   double* pdf)
{
  /* Numerical parameters of the dichotomy search */
  const size_t MAX_ITER = 100;
  const double EPSILON_LAMBDA_M = 1e-15;
  const double EPSILON_BF = 1e-6;

  /* Local variables */
  double bf = 0;
  double bf_prev = 0;
  double bf_min_max = 0;
  double lambda_m = 0;
  double lambda_m_prev = 0;
  double lambda_m_min = 0;
  double lambda_m_max = 0;
  double range_m[2] = {0, 0};
  size_t i;

  /* Check precondition */
  ASSERT(ran_lw && func_name);
  ASSERT(range && range[0] < range[1]);
  ASSERT(0 <= r && r < 1);

  /* Convert the wavelength range in meters */
  range_m[0] = range[0] * 1.0e-9;
  range_m[1] = range[1] * 1.0e-9;

  /* Setup the dichotomy search */
  lambda_m_min = range_m[0];
  lambda_m_max = range_m[1];
  bf_min_max = blackbody_fraction
    (range_m[0], range_m[1], ran_lw->ref_temperature);

  /* Numerically search the lambda corresponding to the submitted canonical
   * number */
  FOR_EACH(i, 0, MAX_ITER) {
    double r_test;
    lambda_m = (lambda_m_min + lambda_m_max) * 0.5;
    bf = blackbody_fraction(range_m[0], lambda_m, ran_lw->ref_temperature);

    r_test = bf / bf_min_max;
    if(r_test < r) {
      lambda_m_min = lambda_m;
    } else {
      lambda_m_max = lambda_m;
    }

    if(fabs(lambda_m_prev - lambda_m) < EPSILON_LAMBDA_M
    || fabs(bf_prev - bf) < EPSILON_BF)
      break;

    lambda_m_prev = lambda_m;
    bf_prev = bf;
  }
  if(i >= MAX_ITER) {
    htrdr_log_warn(ran_lw->htrdr,
      "%s: could not sample a wavelength in the range [%g, %g] nanometers "
      "for the reference temperature %g Kelvin.\n",
      func_name, SPLIT2(range), ran_lw->ref_temperature);
  }

  if(pdf) {
    const double Tref = ran_lw->ref_temperature;
    const double B_lambda = planck(lambda_m, lambda_m, Tref);
    const double B_mean = planck(range_m[0], range_m[1], Tref);
    *pdf = B_lambda / (B_mean * (range_m[1]-range_m[0]));
  }

  return lambda_m * 1.0e9; /* Convert in nanometers */
}

static double
ran_lw_sample_discrete
  (const struct htrdr_ran_lw* ran_lw,
   const double r0,
   const double r1,
   const char* func_name,
   double* pdf)
{
  const double* cdf = NULL;
  const double* find = NULL;
  double r0_next = nextafter(r0, DBL_MAX);
  double band_range[2];
  double lambda = 0;
  double pdf_continue = 0;
  double pdf_band = 0;
  size_t cdf_length = 0;
  size_t i;
  ASSERT(ran_lw && ran_lw->nbands != HTRDR_RAN_LW_CONTINUE);
  ASSERT(0 <= r0 && r0 < 1);
  ASSERT(0 <= r1 && r1 < 1);
  (void)func_name;
  (void)pdf_band;

  cdf = darray_double_cdata_get(&ran_lw->cdf);
  cdf_length = darray_double_size_get(&ran_lw->cdf);
  ASSERT(cdf_length > 0);

  /* Use r_next rather than r0 in order to find the first entry that is not less
   * than *or equal* to r0 */
  find = search_lower_bound(&r0_next, cdf, cdf_length, sizeof(double), cmp_dbl);
  ASSERT(find);

  i = (size_t)(find - cdf);
  ASSERT(i < cdf_length && cdf[i] > r0 && (!i || cdf[i-1] <= r0));

  band_range[0] = ran_lw->range[0] + (double)i*ran_lw->band_len;
  band_range[1] = band_range[0] + ran_lw->band_len;

  /* Fetch the pdf of the sampled band */
  pdf_band = darray_double_cdata_get(&ran_lw->pdf)[i];

  /* Uniformly sample a wavelength in the sampled band */
  lambda = band_range[0] + (band_range[1] - band_range[0]) * r1;
  pdf_continue = 1.0 / ((band_range[1] - band_range[0])*1.e-9);

  if(pdf) {
    *pdf = pdf_band * pdf_continue;
  }

  return lambda;
}

static void
release_ran_lw(ref_T* ref)
{
  struct htrdr_ran_lw* ran_lw = NULL;
  ASSERT(ref);
  ran_lw = CONTAINER_OF(ref, struct htrdr_ran_lw, ref);
  darray_double_release(&ran_lw->cdf);
  darray_double_release(&ran_lw->pdf);
  MEM_RM(ran_lw->htrdr->allocator, ran_lw);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_ran_lw_create
  (struct htrdr* htrdr,
   const double range[2], /* Must be included in [1000, 100000] nanometers */
   const size_t nbands, /* # bands used to discretized CDF */
   const double ref_temperature,
   struct htrdr_ran_lw** out_ran_lw)
{
  struct htrdr_ran_lw* ran_lw = NULL;
  double min_band_len = 0;
  res_T res = RES_OK;
  ASSERT(htrdr && range && out_ran_lw && ref_temperature > 0);
  ASSERT(ref_temperature > 0);
  ASSERT(range[0] >= 1000);
  ASSERT(range[1] <= 100000);
  ASSERT(range[0] <= range[1]);

  ran_lw = MEM_CALLOC(htrdr->allocator, 1, sizeof(*ran_lw));
  if(!ran_lw) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "%s: could not allocate long wave random variate data structure -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }
  ref_init(&ran_lw->ref);
  ran_lw->htrdr = htrdr;
  darray_double_init(htrdr->allocator, &ran_lw->cdf);
  darray_double_init(htrdr->allocator, &ran_lw->pdf);

  ran_lw->range[0] = range[0];
  ran_lw->range[1] = range[1];
  ran_lw->ref_temperature = ref_temperature;
  ran_lw->nbands = nbands;

  min_band_len = compute_sky_min_band_len(ran_lw->htrdr->sky, ran_lw->range);

  if(nbands == HTRDR_RAN_LW_CONTINUE) {
    ran_lw->band_len = 0;
  } else {
    ran_lw->band_len = (range[1] - range[0]) / (double)ran_lw->nbands;

    /* Adjust the band length to ensure that each sky spectral interval is
     * overlapped by at least one band */
    if(ran_lw->band_len > min_band_len) {
      ran_lw->band_len = min_band_len;
      ran_lw->nbands = (size_t)ceil((range[1] - range[0]) / ran_lw->band_len);
    }

    res = setup_ran_lw_cdf(ran_lw, FUNC_NAME);
    if(res != RES_OK) goto error;
  }

  htrdr_log(htrdr, "Long wave interval defined on [%g, %g] nanometers.\n",
    range[0], range[1]);

exit:
  *out_ran_lw = ran_lw;
  return res;
error:
  if(ran_lw) htrdr_ran_lw_ref_put(ran_lw);
  goto exit;
}

void
htrdr_ran_lw_ref_get(struct htrdr_ran_lw* ran_lw)
{
  ASSERT(ran_lw);
  ref_get(&ran_lw->ref);
}

void
htrdr_ran_lw_ref_put(struct htrdr_ran_lw* ran_lw)
{
  ASSERT(ran_lw);
  ref_put(&ran_lw->ref, release_ran_lw);
}

double
htrdr_ran_lw_sample
  (const struct htrdr_ran_lw* ran_lw,
   const double r0,
   const double r1,
   double* pdf)
{
  ASSERT(ran_lw);
  if(ran_lw->nbands != HTRDR_RAN_LW_CONTINUE) { /* Discrete */
    return ran_lw_sample_discrete(ran_lw, r0, r1, FUNC_NAME, pdf);
  } else if(eq_eps(ran_lw->range[0], ran_lw->range[1], 1.e-6)) {
    if(pdf) *pdf = 1;
    return ran_lw->range[0];
  } else { /* Continue */
    return ran_lw_sample_continue
      (ran_lw, r0, ran_lw->range, FUNC_NAME, pdf);
  }
}

