/* Copyright (C) 2018-2019, 2022-2025 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2025 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2025 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2025 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2025 Observatoire de Paris
 * Copyright (C) 2022-2025 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2025 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2025 Université Paul Sabatier
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

#include "core/htrdr.h"
#include "core/htrdr_c.h"
#include "core/htrdr_log.h"
#include "core/htrdr_ran_wlen_planck.h"
#include "core/htrdr_spectral.h"

#include <rsys/algorithm.h>
#include <rsys/cstr.h>
#include <rsys/dynamic_array_double.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

#include <math.h> /* nextafter */

struct htrdr_ran_wlen_planck {
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
setup_wlen_planck_ran_cdf
  (struct htrdr_ran_wlen_planck* planck,
   const char* func_name)
{
  double* pdf = NULL;
  double* cdf = NULL;
  double sum = 0;
  size_t i;
  res_T res = RES_OK;
  ASSERT(planck && func_name);
  ASSERT(planck->nbands != HTRDR_WLEN_RAN_PLANCK_CONTINUE);

  res = darray_double_resize(&planck->cdf, planck->nbands);
  if(res != RES_OK) {
    htrdr_log_err(planck->htrdr,
      "%s: Error allocating the CDF of the spectral bands -- %s.\n",
      func_name, res_to_cstr(res));
    goto error;
  }
  res = darray_double_resize(&planck->pdf, planck->nbands);
  if(res != RES_OK) {
    htrdr_log_err(planck->htrdr,
      "%s: Error allocating the pdf of the spectral bands -- %s.\n",
      func_name, res_to_cstr(res));
    goto error;
  }

  cdf = darray_double_data_get(&planck->cdf);
  pdf = darray_double_data_get(&planck->pdf); /* Now save pdf to correct weight */

  htrdr_log(planck->htrdr,
    "Number of bands used to speed up Planck distribution: %lu\n",
    (unsigned long)planck->nbands);

  /* Compute the *unnormalized* probability to sample a small band */
  FOR_EACH(i, 0, planck->nbands) {
    double lambda_lo = planck->range[0] + (double)i * planck->band_len;
    double lambda_hi = MMIN(lambda_lo + planck->band_len, planck->range[1]);
    ASSERT(lambda_lo<= lambda_hi);
    ASSERT(lambda_lo > planck->range[0]
        || eq_eps(lambda_lo, planck->range[0], 1.e-6));
    ASSERT(lambda_lo < planck->range[1]
        || eq_eps(lambda_lo, planck->range[1], 1.e-6));

    /* Convert from nanometer to meter */
    lambda_lo *= 1.e-9;
    lambda_hi *= 1.e-9;

    /* Compute the probability of the current band */
    pdf[i] = htrdr_blackbody_fraction
      (lambda_lo, lambda_hi, planck->ref_temperature);

    /* Update the norm */
    sum += pdf[i];
  }

  /* Compute the cumulative of the previously computed probabilities */
  FOR_EACH(i, 0, planck->nbands) {
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
  ASSERT(eq_eps(cdf[planck->nbands-1], 1, 1.e-6));
  cdf[planck->nbands - 1] = 1.0; /* Handle numerical issue */

exit:
  return res;
error:
  darray_double_clear(&planck->cdf);
  darray_double_clear(&planck->pdf);
  goto exit;
}

static double
wlen_ran_sample_continue
  (const struct htrdr_ran_wlen_planck* planck,
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
  ASSERT(planck && func_name);
  ASSERT(range && range[0] < range[1]);
  ASSERT(0 <= r && r < 1);

  /* Convert the wavelength range in meters */
  range_m[0] = range[0] * 1.0e-9;
  range_m[1] = range[1] * 1.0e-9;

  /* Setup the dichotomy search */
  lambda_m_min = range_m[0];
  lambda_m_max = range_m[1];
  bf_min_max = htrdr_blackbody_fraction
    (range_m[0], range_m[1], planck->ref_temperature);

  /* Numerically search the lambda corresponding to the submitted canonical
   * number */
  FOR_EACH(i, 0, MAX_ITER) {
    double r_test;
    lambda_m = (lambda_m_min + lambda_m_max) * 0.5;
    bf = htrdr_blackbody_fraction
      (range_m[0], lambda_m, planck->ref_temperature);

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
    htrdr_log_warn(planck->htrdr,
      "%s: could not sample a wavelength in the range [%g, %g] nanometers "
      "for the reference temperature %g Kelvin.\n",
      func_name, SPLIT2(range), planck->ref_temperature);
  }

  if(pdf) {
    const double Tref = planck->ref_temperature; /* K */

    /* W/m²/sr/m */
    const double B_lambda = htrdr_planck(lambda_m, lambda_m, Tref);
    const double B_mean = htrdr_planck(range_m[0], range_m[1], Tref);

    *pdf = B_lambda / (B_mean * (range_m[1]-range_m[0]));
    *pdf *= 1.e-9; /* Transform from m⁻¹ to nm⁻¹ */
  }

  return lambda_m * 1.e9; /* Convert in nanometers */
}

static double
wlen_ran_sample_discrete
  (const struct htrdr_ran_wlen_planck* planck,
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
  ASSERT(planck && planck->nbands != HTRDR_WLEN_RAN_PLANCK_CONTINUE);
  ASSERT(0 <= r0 && r0 < 1);
  ASSERT(0 <= r1 && r1 < 1);
  (void)func_name;
  (void)pdf_band;

  cdf = darray_double_cdata_get(&planck->cdf);
  cdf_length = darray_double_size_get(&planck->cdf);
  ASSERT(cdf_length > 0);

  /* Use r_next rather than r0 in order to find the first entry that is not less
   * than *or equal* to r0 */
  find = search_lower_bound(&r0_next, cdf, cdf_length, sizeof(double), cmp_dbl);
  ASSERT(find);

  i = (size_t)(find - cdf);
  ASSERT(i < cdf_length && cdf[i] > r0 && (!i || cdf[i-1] <= r0));

  band_range[0] = planck->range[0] + (double)i*planck->band_len;
  band_range[1] = band_range[0] + planck->band_len;

  /* Fetch the pdf of the sampled band */
  pdf_band = darray_double_cdata_get(&planck->pdf)[i];

  /* Uniformly sample a wavelength in the sampled band */
  lambda = band_range[0] + (band_range[1] - band_range[0]) * r1;
  pdf_continue = 1.0 / (band_range[1] - band_range[0]);

  if(pdf) {
    *pdf = pdf_band * pdf_continue;
  }

  return lambda;
}

static void
release_wlen_planck_ran(ref_T* ref)
{
  struct htrdr_ran_wlen_planck* planck = NULL;
  struct htrdr* htrdr = NULL;
  ASSERT(ref);
  planck = CONTAINER_OF(ref, struct htrdr_ran_wlen_planck, ref);
  darray_double_release(&planck->cdf);
  darray_double_release(&planck->pdf);
  htrdr = planck->htrdr;
  MEM_RM(htrdr_get_allocator(htrdr), planck);
  htrdr_ref_put(planck->htrdr);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_ran_wlen_planck_create
  (struct htrdr* htrdr,
   /* range must be included in [200,1000] nm for shortwave or in [1000,100000]
    * nanometers for longwave (thermal) */
   const double range[2],
   const size_t nbands, /* # bands used to discretized CDF */
   const double ref_temperature,
   struct htrdr_ran_wlen_planck** out_wlen_planck_ran)
{
  struct htrdr_ran_wlen_planck* planck = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && range && out_wlen_planck_ran && ref_temperature > 0);
  ASSERT(ref_temperature > 0);
  ASSERT(range[0] <= range[1]);

  planck = MEM_CALLOC(htrdr->allocator, 1, sizeof(*planck));
  if(!planck) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "%s: could not allocate Planck distribution -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }
  ref_init(&planck->ref);
  darray_double_init(htrdr->allocator, &planck->cdf);
  darray_double_init(htrdr->allocator, &planck->pdf);
  htrdr_ref_get(htrdr);
  planck->htrdr = htrdr;

  planck->range[0] = range[0];
  planck->range[1] = range[1];
  planck->ref_temperature = ref_temperature;
  planck->nbands = nbands;

  if(nbands == HTRDR_WLEN_RAN_PLANCK_CONTINUE) {
    planck->band_len = 0;
  } else {
    planck->band_len = (range[1] - range[0]) / (double)planck->nbands;

    res = setup_wlen_planck_ran_cdf(planck, FUNC_NAME);
    if(res != RES_OK) goto error;
  }

  htrdr_log(htrdr, "Spectral interval defined on [%g, %g] nanometers.\n",
    range[0], range[1]);

exit:
  *out_wlen_planck_ran = planck;
  return res;
error:
  if(planck) {
    htrdr_ran_wlen_planck_ref_put(planck);
    planck = NULL;
  }
  goto exit;
}

void
htrdr_ran_wlen_planck_ref_get(struct htrdr_ran_wlen_planck* planck)
{
  ASSERT(planck);
  ref_get(&planck->ref);
}

void
htrdr_ran_wlen_planck_ref_put(struct htrdr_ran_wlen_planck* planck)
{
  ASSERT(planck);
  ref_put(&planck->ref, release_wlen_planck_ran);
}

double
htrdr_ran_wlen_planck_sample
  (const struct htrdr_ran_wlen_planck* planck,
   const double r0,
   const double r1,
   double* pdf)
{
  ASSERT(planck);
  if(planck->nbands != HTRDR_WLEN_RAN_PLANCK_CONTINUE) { /* Discrete */
    return wlen_ran_sample_discrete(planck, r0, r1, FUNC_NAME, pdf);
  } else if(eq_eps(planck->range[0], planck->range[1], 1.e-6)) {
    if(pdf) *pdf = 1;
    return planck->range[0];
  } else { /* Continue */
    return wlen_ran_sample_continue
      (planck, r0, planck->range, FUNC_NAME, pdf);
  }
}

