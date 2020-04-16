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

#include <rsys/algorithm.h>
#include <rsys/cstr.h>
#include <rsys/dynamic_array.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

#include <math.h> /* nextafter */

struct htrdr_ran_lw {
  struct darray_double cdf;
  struct darray_double pdf;
  double range[2]; /* Boundaries of the handled CIE XYZ color space */
  double band_len; /* Length in nanometers of a band */

  struct htrdr* htrdr;
  ref_T ref;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static res_T
setup_ran_lw
  (struct htrdr_ran_lw* ran_lw,
   const char* func_name,
   const double range[2],
   const size_t nbands,
   const double ref_temperature)
{
  double* pdf = NULL;
  double* cdf = NULL;
  double sum = 0;
  size_t i;
  res_T res = RES_OK;
  ASSERT(ran_lw && func_name && range && range && ref_temperature);
  ASSERT(range[0] >= 1000);
  ASSERT(range[1] <= 100000);
  ASSERT(range[0] < range[1]);

  res = darray_double_resize(&ran_lw->cdf, nbands);
  if(res != RES_OK) {
    htrdr_log_err(ran_lw->htrdr,
      "%s: Error allocating the CDF of the long wave spectral bands -- %s.\n",
      func_name, res_to_cstr(res));
    goto error;
  }

  res = darray_double_resize(&ran_lw->pdf, nbands);
  if(res != RES_OK) {
    htrdr_log_err(ran_lw->htrdr,
      "%s: Error allocating the PDF of the long wave spectral bands -- %s.\n",
      func_name, res_to_cstr(res));
    goto error;
  }

  cdf = darray_double_data_get(&ran_lw->cdf);
  pdf = darray_double_data_get(&ran_lw->pdf);

  /* Compute the *unormalized* probability to sample a long wave band */
  ran_lw->range[0] = range[0];
  ran_lw->range[1] = range[1];
  ran_lw->band_len = (range[1] - range[0]) / (double)nbands;
  FOR_EACH(i, 0, nbands) {
    double lambda_lo = range[0] + (double)i * ran_lw->band_len;
    double lambda_hi = lambda_lo + ran_lw->band_len;
    ASSERT(lambda_lo <= lambda_hi);

    /* Convert from nanometer to meter */
    lambda_lo *= 1.e-9;
    lambda_hi *= 1.e-9;

    /* Compute the probability of the current band */
    pdf[i] = planck(lambda_lo, lambda_hi, ref_temperature);

    /* Update the norm */
    sum += pdf[i];
  }

  /* Compute the cumulative of the previously computed probabilities */
  FOR_EACH(i, 0, nbands) {
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
  ASSERT(eq_eps(cdf[nbands-1], 1, 1.e-6));
  cdf[nbands - 1] = 1.0; /* Handle numerical issue */

exit:
  return res;
error:
  darray_double_clear(&ran_lw->cdf);
  darray_double_clear(&ran_lw->pdf);
  goto exit;
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
   const size_t nbands, /* # bands used to discretisze the CIE tristimulus */
   const double ref_temperature,
   struct htrdr_ran_lw** out_ran_lw)
{
  struct htrdr_ran_lw* ran_lw = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && range && nbands && out_ran_lw);

  ran_lw = MEM_CALLOC(htrdr->allocator, 1, sizeof(*ran_lw));
  if(ran_lw) {
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

  res = setup_ran_lw(ran_lw, FUNC_NAME, range, nbands, ref_temperature);
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
   const double r,
   double* pdf)
{
  const double* cdf = NULL;
  const double* find = NULL;
  double r_next = nextafter(r, DBL_MAX);
  double wlen = 0;
  size_t cdf_length = 0;
  size_t i;
  ASSERT(ran_lw && 0 <= r && r < 1);
  ASSERT(darray_double_size_get(&ran_lw->cdf)
      == darray_double_size_get(&ran_lw->pdf));

  cdf = darray_double_cdata_get(&ran_lw->cdf);
  cdf_length = darray_double_size_get(&ran_lw->cdf);

  /* Use r_next rather than r in order to find the first entry that is not less
   * than *or equal* to r */
  find = search_lower_bound(&r_next, cdf, cdf_length, sizeof(double), cmp_dbl);
  ASSERT(find);

  i = (size_t)(find - cdf);
  ASSERT(i < cdf_length && cdf[i] > r && (!i || cdf[i-1] <= r));

  /* Return the wavelength at the center of the sampled band */
  wlen = ran_lw->range[0] + ran_lw->band_len * ((double)i + 0.5);
  ASSERT(ran_lw->range[0] < wlen && wlen < ran_lw->range[1]);

  if(pdf) {
    *pdf = darray_double_cdata_get(&ran_lw->pdf)[i];
  }

  return wlen;
}

