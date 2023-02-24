/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut Pierre-Simon Laplace
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2023 Observatoire de Paris
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

#include "atmosphere/htrdr_atmosphere_c.h"
#include "atmosphere/htrdr_atmosphere_ground.h"
#include "atmosphere/htrdr_atmosphere_sun.h"

#include "core/htrdr_interface.h"

#include <high_tune/htsky.h>

#include <star/s3d.h>
#include <star/ssf.h>
#include <star/ssp.h>
#include <star/svx.h>

#include <rsys/double2.h>
#include <rsys/float2.h>
#include <rsys/float3.h>

struct scattering_context {
  struct ssp_rng* rng;
  const struct htsky* sky;
  size_t iband; /* Index of the spectral band */
  size_t iquad; /* Index of the quadrature point into the band */

  double Ts; /* Sampled optical thickness */
  double traversal_dst; /* Distance traversed along the ray */
};
static const struct scattering_context SCATTERING_CONTEXT_NULL = {
  NULL, NULL, 0, 0, 0, 0
};

struct transmissivity_context {
  struct ssp_rng* rng;
  const struct htsky* sky;
  size_t iband; /* Index of the spectral */
  size_t iquad; /* Index of the quadrature point into the band */

  double Ts; /* Sampled optical thickness */
  double Tmin; /* Minimal optical thickness */
  double traversal_dst; /* Distance traversed along the ray */

  enum htsky_property prop;
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
  double ks_max;
  int pursue_traversal = 1;
  ASSERT(hit && ctx && !SVX_HIT_NONE(hit) && org && dir && range);
  (void)range;

  ks_max = htsky_fetch_svx_voxel_property(ctx->sky, HTSKY_Ks,
    HTSKY_SVX_MAX, HTSKY_CPNT_MASK_ALL, ctx->iband, ctx->iquad, &hit->voxel);

  ctx->traversal_dst = hit->distance[0];

  /* Iterate until a collision occurs into the voxel or until the ray
   * does not collide the voxel */
  for(;;) {
    /* Compute tau for the current leaf */
    const double vox_dst = hit->distance[1] - ctx->traversal_dst;
    const double T = vox_dst * ks_max;

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
      ctx->traversal_dst += collision_dst;
      ASSERT(ctx->traversal_dst >= hit->distance[0]);
      ASSERT(ctx->traversal_dst <= hit->distance[1]);

      /* Stop the ray whenever the traversal distance without any scattering
       * event is too high. It means the maximum scattering coefficient has a
       * very small value, and the returned radiance is null. This can only
       * happen when the voxel has a [quasi] infinite length in the propagation
       * direction. */
      if(ctx->traversal_dst > 1e9) break;

      /* Compute the world space position where a collision may occur */
      pos[0] = org[0] + ctx->traversal_dst * dir[0];
      pos[1] = org[1] + ctx->traversal_dst * dir[1];
      pos[2] = org[2] + ctx->traversal_dst * dir[2];

      ks = htsky_fetch_raw_property(ctx->sky, HTSKY_Ks,
        HTSKY_CPNT_MASK_ALL, ctx->iband, ctx->iquad, pos, -DBL_MAX, DBL_MAX);

      /* Handle the case that ks_max is not *really* the max */
      proba = ks / ks_max;

      if(ssp_rng_canonical(ctx->rng) < proba) {/* Collide <=> real scattering */
        pursue_traversal = 0;
        break;
      } else { /* Null collision */
        ctx->Ts = ssp_ran_exp(ctx->rng, 1); /* Sample a new optical thickness */
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
  int comp_mask = HTSKY_CPNT_MASK_ALL;
  double k_max;
  double k_min;
  int pursue_traversal = 1;
  ASSERT(hit && ctx && !SVX_HIT_NONE(hit) && org && dir && range);
  (void)range;

  k_min = htsky_fetch_svx_voxel_property(ctx->sky, ctx->prop,
    HTSKY_SVX_MIN, comp_mask, ctx->iband, ctx->iquad, &hit->voxel);
  k_max = htsky_fetch_svx_voxel_property(ctx->sky, ctx->prop,
    HTSKY_SVX_MAX, comp_mask, ctx->iband, ctx->iquad, &hit->voxel);
  ASSERT(k_min <= k_max);

  ctx->Tmin += (hit->distance[1] - hit->distance[0]) * k_min;
  ctx->traversal_dst = hit->distance[0];

  /* Iterate until a collision occurs into the voxel or until the ray
   * does not collide the voxel */
  for(;;) {
    const double vox_dst = hit->distance[1] - ctx->traversal_dst;
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
      ctx->traversal_dst += collision_dst;
      ASSERT(ctx->traversal_dst >= hit->distance[0]);
      ASSERT(ctx->traversal_dst <= hit->distance[1]);

      /* Compute the world space position where a collision may occur */
      x[0] = org[0] + ctx->traversal_dst * dir[0];
      x[1] = org[1] + ctx->traversal_dst * dir[1];
      x[2] = org[2] + ctx->traversal_dst * dir[2];

      k = htsky_fetch_raw_property(ctx->sky, ctx->prop,
        comp_mask, ctx->iband, ctx->iquad, x, k_min, k_max);
      ASSERT(k >= k_min && k <= k_max);

      proba = (k - k_min) / (k_max - k_min);

      if(ssp_rng_canonical(ctx->rng) < proba) { /* Collide */
        pursue_traversal = 0;
        break;
      } else { /* Null collision */
        ctx->Ts = ssp_ran_exp(ctx->rng, 1); /* Sample a new optical thickness */
      }
    }
  }
  return pursue_traversal;
}

static double
transmissivity
  (struct htrdr_atmosphere* cmd,
   struct ssp_rng* rng,
   const enum htsky_property prop,
   const size_t iband,
   const size_t iquad,
   const double pos[3],
   const double dir[3],
   const double range[2])
{
  struct svx_hit svx_hit;
  struct transmissivity_context transmissivity_ctx = TRANSMISSION_CONTEXT_NULL;

  ASSERT(cmd && rng && pos && dir && range);

  transmissivity_ctx.rng = rng;
  transmissivity_ctx.sky = cmd->sky;
  transmissivity_ctx.iband = iband;
  transmissivity_ctx.iquad = iquad;
  transmissivity_ctx.Ts = ssp_ran_exp(rng, 1); /* Sample an optical thickness */
  transmissivity_ctx.prop = prop;

  /* Compute the transmissivity */
  HTSKY(trace_ray(cmd->sky, pos, dir, range, NULL,
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
atmosphere_compute_radiance_sw
  (struct htrdr_atmosphere* cmd,
   const size_t ithread,
   struct ssp_rng* rng,
   const int cpnt_mask, /* Combination of enum htrdr_radiance_cpnt_flag */
   const double pos_in[3],
   const double dir_in[3],
   const double wlen, /* In nanometer */
   const size_t iband,
   const size_t iquad)
{
  struct s3d_hit s3d_hit = S3D_HIT_NULL;
  struct s3d_hit s3d_hit_tmp = S3D_HIT_NULL;
  struct s3d_hit s3d_hit_prev = S3D_HIT_NULL;
  struct svx_hit svx_hit = SVX_HIT_NULL;
  struct ssf_phase* phase_hg = NULL;
  struct ssf_phase* phase_rayleigh = NULL;

  double pos[3];
  double dir[3];
  double range[2];
  double pos_next[3];
  double dir_next[3];
  double band_bounds[2]; /* In nanometers */

  double R;
  double r; /* Random number */
  double wo[3]; /* -dir */
  double pdf;
  double Tr; /* Overall transmissivity */
  double Tr_abs; /* Absorption transmissivity */
  double L_sun; /* Sun radiance in W.m^-2.sr^-1 */
  double sun_dir[3];
  double ksi = 1; /* Throughput */
  double w = 0; /* MC weight */
  double g = 0; /* Asymmetry parameter of the HG phase function */

  ASSERT(cmd && rng && pos_in && dir_in);
  ASSERT(cmd->spectral_type == HTRDR_SPECTRAL_SW
      || cmd->spectral_type == HTRDR_SPECTRAL_SW_CIE_XYZ);

  /* Create the Henyey-Greenstein phase function */
  CHK(RES_OK == ssf_phase_create
    (htrdr_get_thread_allocator(cmd->htrdr, ithread),
     &ssf_phase_hg,
     &phase_hg));

  /* Create the Rayleigh phase function */
  CHK(RES_OK == ssf_phase_create
    (htrdr_get_thread_allocator(cmd->htrdr, ithread),
     &ssf_phase_rayleigh,
     &phase_rayleigh));

  /* Setup the phase function for this wavelength */
  g = htsky_fetch_per_wavelength_particle_phase_function_asymmetry_parameter
    (cmd->sky, wlen);
  SSF(phase_hg_setup(phase_hg, g));

  /* Fetch sun properties. Note that the sun spectral data are defined by bands
   * that, actually are the same of the SW spectral bands defined in the
   * default "ecrad_opt_prot.txt" file provided by the HTGOP project. */
  htsky_get_spectral_band_bounds(cmd->sky, iband, band_bounds);
  ASSERT(band_bounds[0] <= wlen  && wlen <= band_bounds[1]);
  L_sun = htrdr_atmosphere_sun_get_radiance(cmd->sun, wlen);
  d3_set(pos, pos_in);
  d3_set(dir, dir_in);

  if((cpnt_mask & ATMOSPHERE_RADIANCE_DIRECT) /* Handle direct contribution */
  && htrdr_atmosphere_sun_is_dir_in_solar_cone(cmd->sun, dir)) {
    /* Check that the ray is not occluded along the submitted range */
    d2(range, 0, FLT_MAX);
    HTRDR(atmosphere_ground_trace_ray
      (cmd->ground, pos, dir, range, NULL, &s3d_hit_tmp));
    if(!S3D_HIT_NONE(&s3d_hit_tmp)) {
      Tr = 0;
    } else {
      Tr = transmissivity
        (cmd, rng, HTSKY_Kext, iband, iquad , pos, dir, range);
      w = L_sun * Tr;
    }
  }

  if((cpnt_mask & ATMOSPHERE_RADIANCE_DIFFUSE) == 0)
    goto exit; /* Discard diffuse contribution */

  /* Radiative random walk */
  for(;;) {
    struct scattering_context scattering_ctx = SCATTERING_CONTEXT_NULL;
    struct ssf_bsdf* bsdf = NULL;
    struct ssf_phase* phase;
    double N[3];
    double bounce_reflectivity = 1;
    double sun_dir_pdf;
    int surface_scattering = 0; /* Define if hit a surface */
    int bsdf_type = 0;

    /* Find the first intersection with a surface */
    d2(range, 0, DBL_MAX);
    HTRDR(atmosphere_ground_trace_ray
      (cmd->ground, pos, dir, range, &s3d_hit_prev, &s3d_hit));

    /* Sample an optical thickness */
    scattering_ctx.Ts = ssp_ran_exp(rng, 1);

    /* Setup the remaining scattering context fields */
    scattering_ctx.rng = rng;
    scattering_ctx.sky = cmd->sky;
    scattering_ctx.iband = iband;
    scattering_ctx.iquad = iquad;

    /* Define if a scattering event occurs */
    d2(range, 0, s3d_hit.distance);
    HTSKY(trace_ray(cmd->sky, pos, dir, range, NULL,
      scattering_hit_filter, &scattering_ctx, iband, iquad, &svx_hit));

    /* No scattering and no surface reflection. Stop the radiative random walk */
    if(S3D_HIT_NONE(&s3d_hit) && SVX_HIT_NONE(&svx_hit)) {
      break;
    }
    ASSERT(SVX_HIT_NONE(&svx_hit)
        || (  svx_hit.distance[0] <= scattering_ctx.traversal_dst
           && svx_hit.distance[1] >= scattering_ctx.traversal_dst));

    /* Negate the incoming dir to match the convention of the SSF library */
    d3_minus(wo, dir);

    /* Define if the scattering occurs at a surface */
    surface_scattering = SVX_HIT_NONE(&svx_hit);

    /* Compute the new position */
    pos_next[0] = pos[0] + dir[0]*scattering_ctx.traversal_dst;
    pos_next[1] = pos[1] + dir[1]*scattering_ctx.traversal_dst;
    pos_next[2] = pos[2] + dir[2]*scattering_ctx.traversal_dst;

    /* Define the previous hit surface used to avoid self hit */
    s3d_hit_prev = surface_scattering ? s3d_hit : S3D_HIT_NULL;

    /* Define the absorption transmissivity from the current position to the
     * next position */
    d2(range, 0, scattering_ctx.traversal_dst);
    Tr_abs = transmissivity
      (cmd, rng, HTSKY_Ka, iband, iquad, pos, dir, range);
    if(Tr_abs <= 0) break;

    /* Sample the scattering direction */
    if(surface_scattering) { /* Scattering at a surface */
      struct htrdr_interface interf = HTRDR_INTERFACE_NULL;
      const struct htrdr_mtl* mtl = NULL;

      /* Fetch the hit interface material and build its BSDF */
      htrdr_atmosphere_ground_get_interface(cmd->ground, &s3d_hit, &interf);
      mtl = htrdr_interface_fetch_hit_mtl(&interf, dir, &s3d_hit);
      HTRDR(mtl_create_bsdf(cmd->htrdr, mtl, ithread, wlen, rng, &bsdf));

      /* Revert the normal if necessary to match the SSF convention */
      d3_normalize(N, d3_set_f3(N, s3d_hit.normal));
      if(d3_dot(N, wo) < 0) d3_minus(N, N);

      /* Sample scattering direction */
      bounce_reflectivity = ssf_bsdf_sample
        (bsdf, rng, wo, N, dir_next, &bsdf_type, &pdf);
      if(!(bsdf_type & SSF_REFLECTION)) { /* Handle only reflections */
        bounce_reflectivity = 0;
      }

    } else { /* Scattering in a volume */
      double ks_particle; /* Scattering coefficient of the particles */
      double ks_gas; /* Scattering coefficient of the gaz */
      double ks; /* Overall scattering coefficient */

      ks_gas = htsky_fetch_raw_property(cmd->sky, HTSKY_Ks,
        HTSKY_CPNT_FLAG_GAS, iband, iquad, pos_next, -DBL_MAX, DBL_MAX);
      ks_particle = htsky_fetch_raw_property(cmd->sky, HTSKY_Ks,
        HTSKY_CPNT_FLAG_PARTICLES, iband, iquad, pos_next, -DBL_MAX, DBL_MAX);
      ks = ks_particle + ks_gas;

      r = ssp_rng_canonical(rng);
      if(r < ks_gas / ks) { /* Gas scattering */
        phase = phase_rayleigh;
      } else { /* Cloud scattering */
        phase = phase_hg;
      }

      /* Sample scattering direction */
      ssf_phase_sample(phase, rng, wo, dir_next, NULL);
      ssf_phase_ref_get(phase);
    }

    /* Sample the direction of the direct contribution */
    if(surface_scattering && (bsdf_type & SSF_SPECULAR)) {
      if(!htrdr_atmosphere_sun_is_dir_in_solar_cone(cmd->sun, dir_next)) {
        R = 0; /* No direct lightning */
      } else {
        sun_dir[0] = dir_next[0];
        sun_dir[1] = dir_next[1];
        sun_dir[2] = dir_next[2];
        R = d3_dot(N, sun_dir)<0/* Below the ground*/ ? 0 : bounce_reflectivity;
      }
      sun_dir_pdf = 1.0;
    } else {
      /* Sample a sun direction */
      sun_dir_pdf = htrdr_atmosphere_sun_sample_direction
        (cmd->sun, rng, sun_dir);
      if(surface_scattering) {
        R = d3_dot(N, sun_dir) < 0/* Below the ground */
          ? 0 : ssf_bsdf_eval(bsdf, wo, N, sun_dir) * d3_dot(N, sun_dir);
      } else {
        R = ssf_phase_eval(phase, wo, sun_dir);
      }
    }

    /* The direct contribution to the scattering point is not null so we need
     * to compute the transmissivity from sun to scatt pt */
    if(R <= 0) {
      Tr = 0;
    } else {
      /* Check that the sun is visible from the new position */
      d2(range, 0, FLT_MAX);
      HTRDR(atmosphere_ground_trace_ray
        (cmd->ground, pos_next, sun_dir, range, &s3d_hit_prev, &s3d_hit_tmp));

      /* Compute the sun transmissivity */
      if(!S3D_HIT_NONE(&s3d_hit_tmp)) {
        Tr = 0;
      } else {
        Tr = transmissivity
          (cmd, rng, HTSKY_Kext, iband, iquad, pos_next, sun_dir, range);
      }
    }

    /* Release the scattering function */
    if(surface_scattering) {
      SSF(bsdf_ref_put(bsdf));
    } else {
      SSF(phase_ref_put(phase));
    }

    /* Update the MC weight */
    ksi *= Tr_abs;
    w += ksi * L_sun * Tr * R / sun_dir_pdf;

    /* Russian roulette wrt surface scattering */
    if(surface_scattering && ssp_rng_canonical(rng) >= bounce_reflectivity)
      break;

    /* Setup the next random walk state */
    d3_set(pos, pos_next);
    d3_set(dir, dir_next);
  }

exit:
  SSF(phase_ref_put(phase_hg));
  SSF(phase_ref_put(phase_rayleigh));
  return w;
}

