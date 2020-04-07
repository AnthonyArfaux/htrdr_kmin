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
#include "htrdr_c.h"
#include "htrdr_interface.h"
#include "htrdr_ground.h"
#include "htrdr_solve.h"

#include <high_tune/htsky.h>

#include <star/s3d.h>
#include <star/ssf.h>
#include <star/ssp.h>
#include <star/svx.h>

#include <rsys/double2.h>
#include <rsys/double3.h>

enum event {
  EVENT_ABSORPTION,
  EVENT_SCATTERING,
  EVENT_NONE
};

struct filter_context {
  struct ssp_rng* rng;
  const struct htsky* htsky;
  size_t iband; /* Index of the spectral band */
  size_t iquad; /* Index of the quadrature point into the band */

  double Ts; /* Sampled optical thickness */
  double traversal_dst; /* Distance traversed along the ray */

  enum event event_type;
};
static const struct filter_context FILTER_CONTEXT_NULL = {
  NULL, NULL, 0, 0, 0.0, 0.0, EVENT_NONE
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static int
hit_filter
  (const struct svx_hit* hit,
   const double org[3],
   const double dir[3],
   const double range[2],
   void* context)
{
  struct filter_context* ctx = context;
  double kext_max;
  int pursue_traversal = 1;
  ASSERT(hit && ctx && !SVX_HIT_NONE(hit) && org && dir && range);
  (void)range;

  kext_max = htsky_fetch_svx_voxel_property(ctx->htsky, HTSKY_Kext,
    HTSKY_SVX_MAX, HTSKY_CPNT_MASK_ALL, ctx->iband, ctx->iquad, &hit->voxel);

  ctx->traversal_dst = hit->distance[0];

  for(;;) {
    double vox_dst = hit->distance[1] - ctx->traversal_dst;
    const double T = vox_dst * kext_max;

    /* A collision occurs behind `vox_dst' */
    if(ctx->Ts > T) {
      ctx->Ts -= T;
      ctx->traversal_dst = hit->distance[1];
      pursue_traversal = 1;
      break;

    /* A real/null collision occurs before `vox_dst' */
    } else {
      const double collision_dst = ctx->Ts / kext_max;
      double pos[3];
      double ks, ka;
      double r;
      double proba_abs;
      double proba_sca;

      /* Compute the traversed distance up to the challenged collision */
      ctx->traversal_dst += collision_dst;
      ASSERT(ctx->traversal_dst >= hit->distance[0]);
      ASSERT(ctx->traversal_dst <= hit->distance[1]);

      /* Compute the world space position where a collision may occur */
      pos[0] = org[0] + ctx->traversal_dst * dir[0];
      pos[1] = org[1] + ctx->traversal_dst * dir[1];
      pos[2] = org[2] + ctx->traversal_dst * dir[2];

      ka = htsky_fetch_raw_property(ctx->htsky, HTSKY_Ka, HTSKY_CPNT_MASK_ALL,
        ctx->iband, ctx->iquad, pos, -DBL_MAX, DBL_MAX);
      ks = htsky_fetch_raw_property(ctx->htsky, HTSKY_Ks, HTSKY_CPNT_MASK_ALL,
        ctx->iband, ctx->iquad, pos, -DBL_MAX, DBL_MAX);

      r = ssp_rng_canonical(ctx->rng);
      proba_abs = ka / kext_max;
      proba_sca = ks / kext_max;
      if(r < proba_abs) { /* Absorption event */
        pursue_traversal = 0;
        ctx->event_type = EVENT_ABSORPTION;
        break;
      } else if(r < proba_abs + proba_sca) { /* Scattering event */
        pursue_traversal = 0;
        ctx->event_type = EVENT_SCATTERING;
        break;
      } else { /* Null collision */
        ctx->Ts = ssp_ran_exp(ctx->rng, 1); /* Sample a  new optical thickness */
      }
    }
  }
  return pursue_traversal;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
double
htrdr_compute_radiance_lw
  (struct htrdr* htrdr,
   const size_t ithread,
   struct ssp_rng* rng,
   const double pos_in[3],
   const double dir_in[3],
   const size_t iband,
   const size_t iquad)
{
  struct s3d_hit s3d_hit = S3D_HIT_NULL;
  struct s3d_hit s3d_hit_prev = S3D_HIT_NULL;
  struct svx_hit svx_hit = SVX_HIT_NULL;
  struct ssf_phase* phase_hg = NULL;

  double wo[3];
  double pos[3];
  double dir[3];
  double range[2];
  double pos_next[3];
  double dir_next[3];
  double band_bounds[2]; /* In nanometers */
  double band_bounds_m[2]; /* In meters */
  double wlen;
  double temperature;
  double g;
  double w = 0; /* Weight */

  ASSERT(htrdr && rng && pos_in && dir_in && ithread < htrdr->nthreads);

  /* Retrieve the band boundaries */
  htsky_get_spectral_band_bounds(htrdr->sky, iband, band_bounds);

  /* Arbitrarly use the wavelength at the center of the band as the wavelength
   * to use for data that depend on wavelength rather than spectral band as the
   * BRDF */
  wlen = (band_bounds[0] + band_bounds[1]) * 0.5;
  band_bounds_m[0] = band_bounds[0] * 1e-9;
  band_bounds_m[1] = band_bounds[1] * 1e-9;

  /* Setup the phase function for this spectral band & quadrature point */
  CHK(RES_OK == ssf_phase_create
    (&htrdr->lifo_allocators[ithread], &ssf_phase_hg, &phase_hg));
  g = htsky_fetch_particle_phase_function_asymmetry_parameter
    (htrdr->sky, iband, iquad);
  SSF(phase_hg_setup(phase_hg, g));

  /* Initialise the random walk */
  d3_set(pos, pos_in);
  d3_set(dir, dir_in);
  d2(range, 0, INF);

  for(;;) {
    struct filter_context ctx = FILTER_CONTEXT_NULL;

    /* Sample an optical thickness */
    ctx.Ts = ssp_ran_exp(rng, 1);

    /* Setup the remaining fields of the hit filter context */
    ctx.rng = rng;
    ctx.htsky = htrdr->sky;
    ctx.iband = iband;
    ctx.iquad = iquad;

    /* Found the first intersection with the surface geometry */
    HTRDR(ground_trace_ray
      (htrdr->ground, pos, dir, range, &s3d_hit_prev, &s3d_hit));

    /* Fit the ray range to the surface distance along the ray */
    range[0] = 0;
    range[1] = s3d_hit.distance;

    /* Trace a ray into the participating media */
    HTSKY(trace_ray(htrdr->sky, pos, dir, range, NULL,
      hit_filter, &ctx, iband, iquad, &svx_hit));

    /* No scattering and no surface reflection.
     * Congratulation !! You are in space. */
    if(S3D_HIT_NONE(&s3d_hit) && SVX_HIT_NONE(&svx_hit)) {
      w = 0;
      break;
    }

    /* Compute the next position */
    pos_next[0] = pos[0] + dir[0]*ctx.traversal_dst;
    pos_next[1] = pos[1] + dir[1]*ctx.traversal_dst;
    pos_next[2] = pos[2] + dir[2]*ctx.traversal_dst;

    /* Absorption event. Stop the realisation */
    if(ctx.event_type == EVENT_ABSORPTION) {
      ASSERT(!SVX_HIT_NONE(&svx_hit));
      temperature = htsky_fetch_temperature(htrdr->sky, pos_next);
      w = planck(band_bounds_m[0], band_bounds_m[1], temperature);
      break;
    }

    /* Negate the incoming dir to match the convention of the SSF library */
    d3_minus(wo, dir);

    /* Scattering in the volume */
    if(ctx.event_type == EVENT_SCATTERING) {
      ASSERT(!SVX_HIT_NONE(&svx_hit));
      ssf_phase_sample(phase_hg, rng, wo, dir_next, NULL);
      s3d_hit_prev = S3D_HIT_NULL;

    /* Scattering at a surface */
    } else {
      struct htrdr_interface interf = HTRDR_INTERFACE_NULL;
      struct ssf_bsdf* bsdf = NULL;
      double bounce_reflectivity = 0;
      double N[3];
      int type;
      ASSERT(ctx.event_type == EVENT_NONE);
      ASSERT(!S3D_HIT_NONE(&s3d_hit));

      /* Fetch the hit interface and build its BSDF */
      htrdr_ground_get_interface(htrdr->ground, &s3d_hit, &interf);
      HTRDR(interface_create_bsdf
        (htrdr, &interf, ithread, wlen, pos_next, dir, rng, &s3d_hit, &bsdf));

      d3_normalize(N, d3_set_f3(N, s3d_hit.normal));
      if(d3_dot(N, wo) < 0) d3_minus(N, N);

      bounce_reflectivity = ssf_bsdf_sample
        (bsdf, rng, wo, N, dir_next, &type, NULL);
      if(!(type & SSF_REFLECTION)) { /* Handle only reflections */
        bounce_reflectivity = 0;
      }

      /* Release the BSDF */
      SSF(bsdf_ref_put(bsdf));

      if(ssp_rng_canonical(rng) >= bounce_reflectivity) { /* Absorbed at boundary */
        /* Retrieve the temperature of the interface. Anyway, since we do not
         * have this data yet, we use the temperature of the sky at the current
         * position as the temperature of the surface. */
        temperature = htsky_fetch_temperature(htrdr->sky, pos_next);
        if(temperature <= 0) {
          w = 0;
        } else {
          w = planck(band_bounds_m[0], band_bounds_m[1], temperature);
        }
        break;
      }

      s3d_hit_prev = s3d_hit;
    }
    d3_set(pos, pos_next);
    d3_set(dir, dir_next);
  }
  SSF(phase_ref_put(phase_hg));
  return w;
}

