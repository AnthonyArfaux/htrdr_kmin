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
  struct darray_double cdf;
  /* Probas to sample a sky band overlapped by the range */
  struct darray_double sky_bands_pdf;
  double range[2]; /* Boundaries of the handled CIE XYZ color space */
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

  cdf = darray_double_data_get(&ran_lw->cdf);
  pdf = cdf; /* Alias the CDF by the PDF since only one array is necessary */

  /* Compute the *unormalized* probability to sample a long wave band */
  FOR_EACH(i, 0, ran_lw->nbands) {
    double lambda_lo = ran_lw->range[0] + (double)i * ran_lw->band_len;
    double lambda_hi = lambda_lo + ran_lw->band_len;
    ASSERT(lambda_lo <= lambda_hi);

    /* Convert from nanometer to meter */
    lambda_lo *= 1.e-9;
    lambda_hi *= 1.e-9;

    /* Compute the probability of the current band */
    pdf[i] = planck(lambda_lo, lambda_hi, ran_lw->ref_temperature);

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
  goto exit;
}

static res_T
compute_sky_bands_pdf(struct htrdr_ran_lw* ran_lw, const char* func_name)
{
  double* pdf = NULL;
  double sum = 0;
  size_t i;
  size_t nbands;
  res_T res = RES_OK;
  ASSERT(ran_lw && htsky_is_long_wave(ran_lw->htrdr->sky) && func_name);

  nbands = htsky_get_spectral_bands_count(ran_lw->htrdr->sky);

  res = darray_double_resize(&ran_lw->sky_bands_pdf, nbands);
  if(res != RES_OK) {
    htrdr_log_err(ran_lw->htrdr,
      "%s: error allocating the PDF of the long wave spectral bands "
      "of the sky -- %s.\n",
      func_name, res_to_cstr(res));
    goto error;
  }

  pdf = darray_double_data_get(&ran_lw->sky_bands_pdf);

  /* Compute the *unormalized* probability to sample a long wave band */
  sum = 0;
  FOR_EACH(i, 0, nbands) {
    const size_t iband = htsky_get_spectral_band_id(ran_lw->htrdr->sky, i);
    double wlens[2];
    HTSKY(get_spectral_band_bounds(ran_lw->htrdr->sky, iband, wlens));

    /* Convert from nanometer to meter */
    wlens[0] = wlens[0] * 1.e-9;
    wlens[1] = wlens[1] * 1.e-9;

    /* Compute the probability of the current band */
    pdf[i] = planck(wlens[0], wlens[1], ran_lw->ref_temperature);

    /* Update the norm */
    sum += pdf[i];
  }

  /* Normalize the probabilities */
  FOR_EACH(i, 0, nbands)  pdf[i] /= sum;

exit:
  return res;
error:
  darray_double_clear(&ran_lw->sky_bands_pdf);
  goto exit;

}

static double
ran_lw_sample_discreet
  (const struct htrdr_ran_lw* ran_lw,
   const double r,
   const char* func_name)
{
  const double* cdf = NULL;
  const double* find = NULL;
  double r_next = nextafter(r, DBL_MAX);
  double wlen = 0;
  size_t cdf_length = 0;
  size_t i;
  ASSERT(ran_lw && ran_lw->nbands != HTRDR_RAN_LW_CONTINUE);
  ASSERT(0 <= r && r < 1);
  (void)func_name;

  cdf = darray_double_cdata_get(&ran_lw->cdf);
  cdf_length = darray_double_size_get(&ran_lw->cdf);
  ASSERT(cdf_length > 0);

  /* Use r_next rather than r in order to find the first entry that is not less
   * than *or equal* to r */
  find = search_lower_bound(&r_next, cdf, cdf_length, sizeof(double), cmp_dbl);
  ASSERT(find);

  i = (size_t)(find - cdf);
  ASSERT(i < cdf_length && cdf[i] > r && (!i || cdf[i-1] <= r));

  /* Return the wavelength at the center of the sampled band */
  wlen = ran_lw->range[0] + ran_lw->band_len * ((double)i + 0.5);
  ASSERT(ran_lw->range[0] < wlen && wlen < ran_lw->range[1]);

  return wlen;
}

static double
ran_lw_sample_continue
  (const struct htrdr_ran_lw* ran_lw,
   const double r,
   const char* func_name)
{
  /* Numerical parameters of the dichotomy search */
  const size_t MAX_ITER = 100;
  const double EPSILON_LAMBDA = 1e-6;
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
  ASSERT(ran_lw && ran_lw->nbands == HTRDR_RAN_LW_CONTINUE && func_name);
  ASSERT(0 <= r && r < 1);

  /* Convert the wavelength range in meters */
  range_m[0] = ran_lw->range[0] / 1.0e9;
  range_m[1] = ran_lw->range[1] / 1.0e9;

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

    if(fabs(lambda_m_prev - lambda_m) < EPSILON_LAMBDA
    || fabs(bf_prev - bf) < EPSILON_BF)
      break;

    lambda_m_prev = lambda_m;
    bf_prev = bf;
  }
  if(i >= MAX_ITER) {
    htrdr_log_warn(ran_lw->htrdr,
      "%s: could not sample a wavelength in the range [%g, %g] nanometers "
      "for the reference temperature %g Kelvin.\n",
      func_name, SPLIT2(ran_lw->range), ran_lw->ref_temperature);
  }

  return lambda_m * 1.0e9; /* Convert in nanometers */
}

static void
release_ran_lw(ref_T* ref)
{
  struct htrdr_ran_lw* ran_lw = NULL;
  ASSERT(ref);
  ran_lw = CONTAINER_OF(ref, struct htrdr_ran_lw, ref);
  darray_double_release(&ran_lw->cdf);
  darray_double_release(&ran_lw->sky_bands_pdf);
  MEM_RM(ran_lw->htrdr->allocator, ran_lw);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_ran_lw_create
  (struct htrdr* htrdr,
   const double range[2], /* Must be included in [1000, 100000] nanometers */
   const size_t nbands, /* # bands used to discretized the CIE tristimulus */
   const double ref_temperature,
   struct htrdr_ran_lw** out_ran_lw)
{
  struct htrdr_ran_lw* ran_lw = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && range && out_ran_lw && ref_temperature > 0);
  ASSERT(ref_temperature > 0);
  ASSERT(range[0] >= 1000);
  ASSERT(range[1] <= 100000);
  ASSERT(range[0] < range[1]);

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
  darray_double_init(htrdr->allocator, &ran_lw->sky_bands_pdf);

  ran_lw->range[0] = range[0];
  ran_lw->range[1] = range[1];
  ran_lw->ref_temperature = ref_temperature;
  ran_lw->nbands = nbands;

  if(nbands == HTRDR_RAN_LW_CONTINUE) {
    ran_lw->band_len = 0;
  } else {
    ran_lw->band_len = (range[1] - range[0]) / (double)nbands;
    res = setup_ran_lw_cdf(ran_lw, FUNC_NAME);
    if(res != RES_OK) goto error;
  }

  res = compute_sky_bands_pdf(ran_lw, FUNC_NAME);
  if(res != RES_OK) goto error;

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
   const double r)
{
  ASSERT(ran_lw);
  if(ran_lw->band_len == 0) {
    return ran_lw_sample_continue(ran_lw, r, FUNC_NAME);
  } else {
    return ran_lw_sample_discreet(ran_lw, r, FUNC_NAME);
  }
}

double
htrdr_ran_lw_get_sky_band_pdf
  (const struct htrdr_ran_lw* ran_lw,
   const size_t iband)
{
  size_t i;
  size_t nbands;
  (void)nbands;
  ASSERT(ran_lw);

  nbands = htsky_get_spectral_bands_count(ran_lw->htrdr->sky);
  ASSERT(nbands);
  ASSERT(iband >= htsky_get_spectral_band_id(ran_lw->htrdr->sky, 0));
  ASSERT(iband <= htsky_get_spectral_band_id(ran_lw->htrdr->sky, nbands-1));

  /* Compute the index of the spectral band */
  i = iband - htsky_get_spectral_band_id(ran_lw->htrdr->sky, 0);

  /* Fetch its PDF */
  ASSERT(i < darray_double_size_get(&ran_lw->sky_bands_pdf));
  return darray_double_cdata_get(&ran_lw->sky_bands_pdf)[i];
}

