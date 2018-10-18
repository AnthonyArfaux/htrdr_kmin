/* Copyright (C) 2018 Université Paul Sabatier, |Meso|Star>
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
#include "htrdr_ground.h"
#include "htrdr_solve.h"
#include "htrdr_sky.h"
#include "htrdr_sun.h"

#include <star/s3d.h>
#include <star/ssf.h>
#include <star/ssp.h>
#include <star/svx.h>

#include <rsys/double2.h>
#include <rsys/float2.h>
#include <rsys/float3.h>

struct scattering_context {
  struct ssp_rng* rng;
  const struct htrdr_sky* sky;
  size_t iband; /* Index of the spectrald */
  size_t iquad; /* Index of the quadrature point into the band */

  double Ts; /* Sampled optical thickness */
  double traversal_dst; /* Distance traversed along the ray */
};
static const struct scattering_context SCATTERING_CONTEXT_NULL = {
  NULL, NULL, 0, 0, 0, 0
};

struct transmissivity_context {
  struct ssp_rng* rng;
  const struct htrdr_sky* sky;
  size_t iband; /* Index of the spectrald */
  size_t iquad; /* Index of the quadrature point into the band */

  double Ts; /* Sampled optical thickness */
  double Tmin; /* Minimal optical thickness */
  double traversal_dst; /* Distance traversed along the ray */

  enum htrdr_sky_property prop;
};
static const struct transmissivity_context TRANSMISSION_CONTEXT_NULL = {
  NULL, NULL, 0, 0, 0, 0, 0, 0
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static int
scattering_hit_filter
  (const struct svx_hit* hit,
   const double org[3],
   const double dir[3],
   const double range[2],
   void* context)
{
  struct scattering_context* ctx = context;
  double vox_dst; /* Distance to traverse into the voxel */
  double ks_max;
  int pursue_traversal = 1;
  ASSERT(hit && ctx && !SVX_HIT_NONE(hit) && org && dir && range);
  (void)range;

  ks_max = htrdr_sky_fetch_svx_voxel_property(ctx->sky, HTRDR_Ks,
    HTRDR_SVX_MAX, HTRDR_ALL_COMPONENTS, ctx->iband, ctx->iquad, &hit->voxel);

  /* Compute the distance to traverse into the voxel */
  vox_dst = hit->distance[1] - hit->distance[0];

  /* Iterate until a collision occurs into the voxel or until the ray
   * does not collide the voxel */
  for(;;) {
    const double T = vox_dst * ks_max; /* Compute tau for the current leaf */

    /* A collision occurs behind `vox_dst' */
    if(ctx->Ts > T) {
      ctx->Ts -= T;
      ctx->traversal_dst = hit->distance[1];
      pursue_traversal = 1;
      break;

    /*  A real/null collision occurs before `vox_dst' */
    } else {
      double pos[3];
      double proba;
      double ks;
      const double collision_dst = ctx->Ts / ks_max;

      /* Compute the traversed distance up to the challenged collision */
      ctx->traversal_dst = hit->distance[0] + collision_dst;
      ASSERT(ctx->traversal_dst >= hit->distance[0]);
      ASSERT(ctx->traversal_dst <= hit->distance[1]);

      /* Compute the world space position where a collision may occur */
      pos[0] = org[0] + ctx->traversal_dst * dir[0];
      pos[1] = org[1] + ctx->traversal_dst * dir[1];
      pos[2] = org[2] + ctx->traversal_dst * dir[2];

      ks = htrdr_sky_fetch_raw_property(ctx->sky, HTRDR_Ks,
        HTRDR_ALL_COMPONENTS, ctx->iband, ctx->iquad, pos, -DBL_MAX, DBL_MAX);

      /* Handle the case that ks_max is not *really* the max */
      proba = ks / ks_max;

      if(ssp_rng_canonical(ctx->rng) < proba) {/* Collide <=> real scattering */
        pursue_traversal = 0;
        break;
      } else { /* Null collision */
        ctx->Ts = ssp_ran_exp(ctx->rng, 1); /* Sample a new optical thickness */
        /* Compute the remaining distance to traverse into the voxel */
        vox_dst = hit->distance[1] - ctx->traversal_dst;
      }
    }
  }
  return pursue_traversal;
}

static int
transmissivity_hit_filter
  (const struct svx_hit* hit,
   const double org[3],
   const double dir[3],
   const double range[2],
   void* context)
{
  struct transmissivity_context* ctx = context;
  int comp_mask = HTRDR_ALL_COMPONENTS;
  double vox_dst; /* Distance to traverse into the voxel */
  double k_max;
  double k_min;
  int pursue_traversal = 1;
  ASSERT(hit && ctx && !SVX_HIT_NONE(hit) && org && dir && range);
  (void)range;

  k_min = htrdr_sky_fetch_svx_voxel_property(ctx->sky, ctx->prop,
    HTRDR_SVX_MIN, comp_mask, ctx->iband, ctx->iquad, &hit->voxel);
  k_max = htrdr_sky_fetch_svx_voxel_property(ctx->sky, ctx->prop,
    HTRDR_SVX_MAX, comp_mask, ctx->iband, ctx->iquad, &hit->voxel);
  ASSERT(k_min <= k_max);

  /* Compute the distance to traverse into the voxel */
  vox_dst = hit->distance[1] - hit->distance[0];
  ctx->Tmin += vox_dst * k_min;

  /* Iterate until a collision occurs into the voxel or until the ray
   * does not collide the voxel */
  for(;;) {
    const double Tdif = vox_dst * (k_max-k_min);

    /* A collision occurs behind `vox_dst' */
    if(ctx->Ts > Tdif) {
      ctx->Ts -= Tdif;
      ctx->traversal_dst = hit->distance[1];
      pursue_traversal = 1;
      break;

    /*  A real/null collision occurs before `vox_dst' */
    } else {
      double x[3];
      double k;
      double proba;
      double collision_dst = ctx->Ts / (k_max - k_min);

      /* Compute the traversed distance up to the challenged collision */
      ctx->traversal_dst = hit->distance[0] + collision_dst;
      ASSERT(ctx->traversal_dst >= hit->distance[0]);
      ASSERT(ctx->traversal_dst <= hit->distance[1]);

      /* Compute the world space position where a collision may occur */
      x[0] = org[0] + ctx->traversal_dst * dir[0];
      x[1] = org[1] + ctx->traversal_dst * dir[1];
      x[2] = org[2] + ctx->traversal_dst * dir[2];

      k = htrdr_sky_fetch_raw_property(ctx->sky, ctx->prop,
        comp_mask, ctx->iband, ctx->iquad, x, k_min, k_max);
      ASSERT(k >= k_min && k <= k_max);

      proba = (k - k_min) / (k_max - k_min);

      if(ssp_rng_canonical(ctx->rng) < proba) { /* Collide */
        pursue_traversal = 0;
        break;
      } else { /* Null collision */
        ctx->Ts = ssp_ran_exp(ctx->rng, 1); /* Sample a new optical thickness */
        /* Compute the remaining distance to traverse into the voxel */
        vox_dst = hit->distance[1] - ctx->traversal_dst;
      }
    }
  }
  return pursue_traversal;
}

static double
transmissivity
  (struct htrdr* htrdr,
   struct ssp_rng* rng,
   const enum htrdr_sky_property prop,
   const size_t iband,
   const size_t iquad,
   const double pos[3],
   const double dir[3],
   const double range[2],
   const struct s3d_hit* hit_prev) /* May be NULL */
{
  struct s3d_hit s3d_hit;
  struct svx_hit svx_hit;
  struct transmissivity_context transmissivity_ctx = TRANSMISSION_CONTEXT_NULL;

  ASSERT(htrdr && rng && pos && dir && range);

  /* Check that the ray is not occlude along the submitted range */
  HTRDR(ground_trace_ray(htrdr->ground, pos, dir, range, hit_prev, &s3d_hit));
  if(!S3D_HIT_NONE(&s3d_hit)) return 0;

  transmissivity_ctx.rng = rng;
  transmissivity_ctx.sky = htrdr->sky;
  transmissivity_ctx.iband = iband;
  transmissivity_ctx.iquad = iquad;
  transmissivity_ctx.Ts = ssp_ran_exp(rng, 1); /* Sample an optical thickness */
  transmissivity_ctx.prop = prop;

  /* Compute the transmissivity */
  HTRDR(sky_trace_ray(htrdr->sky, pos, dir, range, NULL,
    transmissivity_hit_filter, &transmissivity_ctx, iband, iquad, &svx_hit));

  if(SVX_HIT_NONE(&svx_hit)) {
    return transmissivity_ctx.Tmin ? exp(-transmissivity_ctx.Tmin) : 1.0;
  } else {
    return 0;
  }
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
double
htrdr_compute_radiance_sw
  (struct htrdr* htrdr,
   const size_t ithread,
   struct ssp_rng* rng,
   const double pos_in[3],
   const double dir_in[3],
   const size_t iband,
   const size_t iquad)
{
  struct s3d_hit s3d_hit = S3D_HIT_NULL;
  struct svx_hit svx_hit = SVX_HIT_NULL;
  struct s3d_hit s3d_hit_prev = S3D_HIT_NULL;
  struct ssf_phase* phase_hg = NULL;
  struct ssf_phase* phase_rayleigh = NULL;
  struct ssf_bsdf* bsdf = NULL;

  double pos[3];
  double dir[3];
  double range[2];
  double pos_next[3];
  double dir_next[3];
  double band_bounds[2]; /* In nanometers */

  double wlen;
  double R;
  double r; /* Random number */
  double wo[3]; /* -dir */
  double pdf;
  double sun_solid_angle;
  double Tr; /* Overall transmissivity */
  double Tr_abs; /* Absorption transmissivity */
  double L_sun; /* Sun radiance in W.m^-2.sr^-1 */
  double sun_dir[3];
  double ksi = 1; /* Throughput */
  double w = 0; /* MC weight */
  double g = 0; /* Asymmetry parameter of the HG phase function */

  ASSERT(htrdr && rng && pos_in && dir_in && ithread < htrdr->nthreads);

  CHK(RES_OK == ssf_bsdf_create
    (&htrdr->lifo_allocators[ithread], &ssf_lambertian_reflection, &bsdf));
  CHK(RES_OK == ssf_phase_create
    (&htrdr->lifo_allocators[ithread], &ssf_phase_hg, &phase_hg));
  CHK(RES_OK == ssf_phase_create
    (&htrdr->lifo_allocators[ithread], &ssf_phase_rayleigh, &phase_rayleigh));

  SSF(lambertian_reflection_setup(bsdf, 0.5));

  /* Setup the phase function for this spectral band & quadrature point */
  g = htrdr_sky_fetch_particle_phase_function_asymmetry_parameter
    (htrdr->sky, iband, iquad);
  SSF(phase_hg_setup(phase_hg, g));

  /* Fetch sun properties. Arbitrarily use the wavelength at the center of the
   * band to retrieve the sun radiance of the current band. Note that the sun
   * spectral data are defined by bands that, actually are the same of the SW
   * spectral bands defined in the default "ecrad_opt_prot.txt" file provided
   * by the HTGOP project. */
  htrdr_sky_get_sw_spectral_band_bounds(htrdr->sky, iband, band_bounds);
  wlen = (band_bounds[0] + band_bounds[1]) * 0.5;
  sun_solid_angle = htrdr_sun_get_solid_angle(htrdr->sun);
  L_sun = htrdr_sun_get_radiance(htrdr->sun, wlen);

  d3_set(pos, pos_in);
  d3_set(dir, dir_in);

  /* Add the directly contribution of the sun  */
  if(htrdr_sun_is_dir_in_solar_cone(htrdr->sun, dir)) {
    /* Add the direct contribution of the sun */
    d2(range, 0, FLT_MAX);
    Tr = transmissivity
      (htrdr, rng, HTRDR_Kext, iband, iquad , pos, dir, range, NULL);
    w = L_sun * Tr;
  }

  /* Radiative random walk */
  for(;;) {
    struct scattering_context scattering_ctx = SCATTERING_CONTEXT_NULL;

    /* Find the first intersection with a surface */
    d2(range, 0, DBL_MAX);
    HTRDR(ground_trace_ray
      (htrdr->ground, pos, dir, range, &s3d_hit_prev, &s3d_hit));

    /* Sample an optical thickness */
    scattering_ctx.Ts = ssp_ran_exp(rng, 1);

    /* Setup the remaining scattering context fields */
    scattering_ctx.rng = rng;
    scattering_ctx.sky = htrdr->sky;
    scattering_ctx.iband = iband;
    scattering_ctx.iquad = iquad;

    /* Define if a scattering event occurs */
    d2(range, 0, s3d_hit.distance);
    HTRDR(sky_trace_ray(htrdr->sky, pos, dir, range, NULL,
      scattering_hit_filter, &scattering_ctx, iband, iquad, &svx_hit));

    /* No scattering and no surface reflection. Stop the radiative random walk */
    if(S3D_HIT_NONE(&s3d_hit) && SVX_HIT_NONE(&svx_hit)) {
      break;
    }
    ASSERT(SVX_HIT_NONE(&svx_hit)
        || (  svx_hit.distance[0] <= scattering_ctx.traversal_dst
           && svx_hit.distance[1] >= scattering_ctx.traversal_dst));

    /* Negative the incoming dir to match the convention of the SSF library */
    d3_minus(wo, dir);

    /* Compute the new position */
    pos_next[0] = pos[0] + dir[0]*scattering_ctx.traversal_dst;
    pos_next[1] = pos[1] + dir[1]*scattering_ctx.traversal_dst;
    pos_next[2] = pos[2] + dir[2]*scattering_ctx.traversal_dst;

    /* Define the absorption transmissivity from the current position to the
     * next position */
    d2(range, 0, scattering_ctx.traversal_dst);
    Tr_abs = transmissivity
      (htrdr, rng, HTRDR_Ka, iband, iquad, pos, dir, range, &s3d_hit_prev);
    if(Tr_abs <= 0) break;

    /* Sample a sun direction */
    htrdr_sun_sample_direction(htrdr->sun, rng, sun_dir);
    d2(range, 0, FLT_MAX);

    s3d_hit_prev = SVX_HIT_NONE(&svx_hit) ? s3d_hit : S3D_HIT_NULL;

    /* Check that the sun is visible from the new position */
    Tr = transmissivity(htrdr, rng, HTRDR_Kext, iband, iquad, pos_next,
      sun_dir, range, &s3d_hit_prev);

    /* Scattering at a surface */
    if(SVX_HIT_NONE(&svx_hit)) {
      double reflectivity;
      double N[3];
      int type;

      d3_normalize(N, d3_set_f3(N, s3d_hit.normal));
      if(d3_dot(N, wo) < 0) d3_minus(N, N);
      reflectivity = ssf_bsdf_sample
        (bsdf, rng, wo, N, dir_next, &type, &pdf);
      if(ssp_rng_canonical(rng) > reflectivity) break; /* Russian roulette */

      if(d3_dot(N, sun_dir) < 0) { /* Below the ground */
        R = 0;
      } else {
        R = ssf_bsdf_eval(bsdf, wo, N, sun_dir) * d3_dot(N, sun_dir);
      }

    /* Scattering in the medium */
    } else {
      struct ssf_phase* phase;
      double ks_particle; /* Scattering coefficient of the particles */
      double ks_gas; /* Scattering coefficient of the gaz */
      double ks; /* Overall scattering coefficient */

      ks_gas = htrdr_sky_fetch_raw_property(htrdr->sky, HTRDR_Ks, HTRDR_GAS,
        iband, iquad, pos_next, -DBL_MAX, DBL_MAX);
      ks_particle = htrdr_sky_fetch_raw_property(htrdr->sky, HTRDR_Ks,
        HTRDR_PARTICLES, iband, iquad, pos_next, -DBL_MAX, DBL_MAX);
      ks = ks_particle + ks_gas;

      r = ssp_rng_canonical(rng);
      if(r < ks_gas / ks) { /* Gas scattering */
        phase = phase_rayleigh;
      } else { /* Cloud scattering */
        phase = phase_hg;
      }

      ssf_phase_sample(phase, rng, wo, dir_next, NULL);
      R = ssf_phase_eval(phase, wo, sun_dir);
    }

    /* Update the MC weight */
    ksi *= Tr_abs;
    w += ksi * L_sun * sun_solid_angle * Tr * R;

    /* Setup the next random walk state */
    d3_set(pos, pos_next);
    d3_set(dir, dir_next);
  }
  SSF(bsdf_ref_put(bsdf));
  SSF(phase_ref_put(phase_hg));
  SSF(phase_ref_put(phase_rayleigh));
  return w;
}

