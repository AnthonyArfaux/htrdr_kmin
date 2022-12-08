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

#define _POSIX_C_SOURCE 200112L /* nextafter */

#include "core/htrdr_c.h"
#include "core/htrdr_log.h"
#include "core/htrdr_ran_wlen_discrete.h"

#include <rsys/algorithm.h>
#include <rsys/dynamic_array_double.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>

#include <math.h>

struct htrdr_ran_wlen_discrete {
  struct darray_double cumul;
  struct darray_double proba;
  struct darray_double radia; /* In W/m²/sr/m */
  struct darray_double wlens; /* In nm */
  double range[2]; /* Boundaries of the spectral integration interval in nm */
  size_t nbands; /* #bands */

  struct htrdr* htrdr;
  ref_T ref;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static INLINE res_T
check_htrdr_ran_wlen_discrete_create_args
  (const struct htrdr_ran_wlen_discrete_create_args* args)
{
  if(!args) return RES_BAD_ARG;

  /* Invalid number of wavelength */
  if(!args->nwavelengths)
    return RES_BAD_ARG;

  /* Invalid functor */
  if(!args->get)
    return RES_BAD_ARG;

  return RES_OK;
}

static res_T
setup_per_wlen_radiance
  (struct htrdr_ran_wlen_discrete* ran,
   const struct htrdr_ran_wlen_discrete_create_args* args)
{
  double* wlens = NULL;
  double* radia = NULL;
  size_t iwlen = 0;
  res_T res = RES_OK;
  ASSERT(ran && args);

  res = darray_double_resize(&ran->wlens, args->nwavelengths);
  if(res != RES_OK) {
    htrdr_log_err(ran->htrdr,
      "Error allocating discrete wavelength distribution wavelength list.\n");
    goto error;
  }
  res = darray_double_resize(&ran->radia, args->nwavelengths);
  if(res != RES_OK) {
    htrdr_log_err(ran->htrdr,
      "Error allocating discrete wavelength distribution radiance list.\n");
    goto error;
  }

  wlens = darray_double_data_get(&ran->wlens);
  radia = darray_double_data_get(&ran->radia);

  /* Store the discrete values */
  FOR_EACH(iwlen, 0, args->nwavelengths) {
    args->get(args->context, iwlen, wlens+iwlen, radia+iwlen);

    if(iwlen > 0 && wlens[iwlen] <= wlens[iwlen-1]) {
      htrdr_log_err(ran->htrdr,
        "Failed to calculate discrete luminance distribution probabilities. "
        "Wavelengths are not sorted in ascending order.\n");
      res = RES_BAD_ARG;
      goto error;
    }
  }

  /* Setup the spectral range */
  ran->range[0] = wlens[0];
  ran->range[1] = wlens[args->nwavelengths-1];

exit:
  return res;
error:
  darray_double_clear(&ran->wlens);
  darray_double_clear(&ran->radia);
  goto exit;
}

static res_T
setup_distribution
  (struct htrdr_ran_wlen_discrete* ran,
   const struct htrdr_ran_wlen_discrete_create_args* args)
{
  double* cumul = NULL;
  double* proba = NULL;
  double sum = 0;
  size_t iband;
  res_T res = RES_OK;
  ASSERT(ran && check_htrdr_ran_wlen_discrete_create_args(args));
  ASSERT(ran->nbands >= 1); /* At least one band */

  res = darray_double_resize(&ran->cumul, ran->nbands);
  if(res != RES_OK) {
    htrdr_log_err(ran->htrdr,
      "Error allocating the cumulative discrete wavelength distribution.\n");
    goto error;
  }
  res = darray_double_resize(&ran->proba, ran->nbands);
  if(res != RES_OK) {
    htrdr_log_err(ran->htrdr,
      "Error allocating the discrete wavelength distribution probabilities.\n");
    goto error;
  }

  cumul = darray_double_data_get(&ran->cumul);
  proba = darray_double_data_get(&ran->proba);

  /* Compute the unormalized probabilities to sample a band */
  FOR_EACH(iband, 0, ran->nbands) {
    const size_t iw0 = iband+0;
    const size_t iw1 = iband+1;
    double w0, L0;
    double w1, L1;
    double area;

    args->get(args->context, iw0, &w0, &L0);
    args->get(args->context, iw1, &w1, &L1);
    ASSERT(w0 < w1);

    area = (L0 + L1) * (w1-w0) * 0.5;

    proba[iband] = area;
    sum += area;
  }

  htrdr_log(ran->htrdr, "Discrete radiance integral = %g W/m²/sr\n", sum);

  /* Normalize the probabilities and setup the cumulative */
  FOR_EACH(iband, 0, ran->nbands) {
    proba[iband] /= sum;
    if(iband == 0) {
      cumul[iband] = proba[iband];
    } else {
      cumul[iband] = proba[iband] + cumul[iband-1];
      ASSERT(cumul[iband] >= cumul[iband-1]);
    }
  }
  ASSERT(eq_eps(cumul[ran->nbands-1], 1, 1e-6));
  cumul[ran->nbands-1] = 1.0; /* Fix numerical imprecision */

exit:
  return res;
error:
  darray_double_clear(&ran->proba);
  darray_double_clear(&ran->cumul);
  goto exit;
}

static void
release_discrete(ref_T* ref)
{
  struct htrdr_ran_wlen_discrete* ran = NULL;
  struct htrdr* htrdr = NULL;
  ASSERT(ref);
  ran = CONTAINER_OF(ref, struct htrdr_ran_wlen_discrete, ref);
  darray_double_release(&ran->cumul);
  darray_double_release(&ran->proba);
  darray_double_release(&ran->radia);
  darray_double_release(&ran->wlens);
  htrdr = ran->htrdr;
  MEM_RM(htrdr_get_allocator(htrdr), ran);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Exported functions
 ******************************************************************************/
res_T
htrdr_ran_wlen_discrete_create
  (struct htrdr* htrdr,
   const struct htrdr_ran_wlen_discrete_create_args* args,
   struct htrdr_ran_wlen_discrete** out_ran)
{
  struct htrdr_ran_wlen_discrete* ran = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && args && out_ran);
  ASSERT(check_htrdr_ran_wlen_discrete_create_args(args) == RES_OK);

  ran = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*ran));
  if(!ran) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "%s: error allocating discrete wavelength distribution\n",
      FUNC_NAME);
    goto error;
  }

  ref_init(&ran->ref);
  darray_double_init(htrdr_get_allocator(htrdr), &ran->cumul);
  darray_double_init(htrdr_get_allocator(htrdr), &ran->proba);
  darray_double_init(htrdr_get_allocator(htrdr), &ran->radia);
  darray_double_init(htrdr_get_allocator(htrdr), &ran->wlens);
  htrdr_ref_get(htrdr);
  ran->htrdr = htrdr;

  ran->nbands = args->nwavelengths - 1;

  res = setup_per_wlen_radiance(ran, args);
  if(res != RES_OK) goto error;
  res = setup_distribution(ran, args);
  if(res != RES_OK) goto error;

exit:
  *out_ran = ran;
  return res;
error:
  if(ran) { htrdr_ran_wlen_discrete_ref_put(ran); ran = NULL; }
  goto exit;
}

void
htrdr_ran_wlen_discrete_ref_get(struct htrdr_ran_wlen_discrete* ran)
{
  ASSERT(ran);
  ref_get(&ran->ref);
}

void
htrdr_ran_wlen_discrete_ref_put(struct htrdr_ran_wlen_discrete* ran)
{
  ASSERT(ran);
  ref_put(&ran->ref, release_discrete);
}

double
htrdr_ran_wlen_discrete_sample
  (struct htrdr_ran_wlen_discrete* ran,
   const double r0,
   const double r1,
   double* pdf) /* In nm⁻¹ */
{
  double* find = NULL;
  const double* proba = NULL;
  const double* cumul = NULL;
  const double* radia = NULL;
  const double* wlens = NULL;
  double w0, w1; /* Wavelengths of the sampled band in nm */
  double L0, L1; /* Radiance of the sampled band in W/m²/sr/m */
  double a, b, c, tmp;
  double lambda; /* Sampled wavelength in nm */
  size_t iband = 0;
  double r0_next = nextafter(r0, DBL_MAX);
  ASSERT(0 <= r0 && r0 < 1);
  ASSERT(0 <= r1 && r1 < 1);

  cumul = darray_double_cdata_get(&ran->cumul);
  proba = darray_double_cdata_get(&ran->proba);
  radia = darray_double_cdata_get(&ran->radia);
  wlens = darray_double_cdata_get(&ran->wlens);

  /* Sample a band. Use r0_next rather than r0 to find the first entry that is
   * not less than *or equal* to r0 */
  find = search_lower_bound
    (&r0_next, cumul, ran->nbands, sizeof(double), cmp_dbl);
  ASSERT(find);

  iband = (size_t)(find - cumul);
  ASSERT(iband < ran->nbands);
  ASSERT(cumul[iband] > r0 && (!iband || cumul[iband-1] <= r0));

  /* Retrieve the boundaries of the sampled band */
  w0 = wlens[iband+0];
  w1 = wlens[iband+1];
  L0 = radia[iband+0];
  L1 = radia[iband+1];

  /* Compute the parameters of the quadratic equation */
  tmp = (L1-L0) / (w1-w0);
  a = 0.5 * tmp;
  b = L0 - w0 * tmp;
  c = 0.5 * tmp * w0*w0 - L0*w0 - r1 * 0.5 * (L0+L1)*(w1-w0);

  /* Sample lambda in the sampled band */
  if(a == 0) {
    lambda = -c/b;
    ASSERT(w0 <= lambda && lambda <= w1);
  } else {
    const double delta = b*b - 4*a*c;
    const double sqrt_delta = sqrt(delta);
    const double root0 = (-b + sqrt_delta) / (2*a);
    const double root1 = (-b - sqrt_delta) / (2*a);

    if(w0 <= root0 && root0 <= w1) {
      ASSERT(root1 < w0 || w1 < root1);
      lambda = root0;
    } else if(w0 <= root1 && root1 <= w1) {
      ASSERT(root0 < w0 || w1 < root1);
      lambda = root1;
    } else {
      FATAL("Unreachable code\n");
    }
  }

  if(*pdf) {
    const double proba_wlen = L0+(L1-L0)/(w1-w0)*(lambda-w0);
    *pdf = proba[iband] * proba_wlen;
  }

  return lambda;
}
