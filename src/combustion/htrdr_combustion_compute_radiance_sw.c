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

#include "combustion/htrdr_combustion_c.h"
#include "combustion/htrdr_combustion_laser.h"

#include <astoria/atrstm.h>

#include <star/ssp.h>

#include <rsys/double2.h>
#include <rsys/double3.h>

struct transmissivity_context {
  struct ssp_rng* rng;
  struct atrstm* medium;
  double wavelength;
  enum atrstm_radcoef radcoef;
};
#define TRANSMISSIVITY_CONTEXT_NULL__ {                                        \
  NULL, /* RNG */                                                              \
  NULL, /* Medium */                                                           \
  0, /* Wavelength */                                                          \
  ATRSTM_RADCOEFS_COUNT__ /* Radiative coefficient */                          \
}
static const struct transmissivity_context TRANSMISSIVITY_CONTEXT_NULL =
  TRANSMISSIVITY_CONTEXT_NULL__;

struct sample_scattering_limited_context {
  struct ssp_rng* rng;
  struct atrstm* medium;
  double wavelength;
  double traversal_dst; /* Distance traversed up to the collision */
  double ksi;
};
#define SAMPLE_SCATTERING_LIMITED_CONTEXT_NULL__ {                             \
  NULL, /* RNG */                                                              \
  NULL, /* Medium */                                                           \
  0, /* Wavelength */                                                          \
  DBL_MAX, /* Traversal dst */                                                 \
  1, /* ksi */                                                                 \
}
static const struct sample_scattering_limited_context 
SAMPLE_SCATTERING_LIMITED_CONTEXT_NULL = 
  SAMPLE_SCATTERING_LIMITED_CONTEXT_NULL__;

struct laser_sheet_scattering {
  double distance; /* Distance along the ray up to a scattering event */
  double ksi;
};
#define LASER_SHEET_SCATTERING_NULL__ {0, 1}
static const struct laser_sheet_scattering LASER_SHEET_SCATTERING_NULL =
  LASER_SHEET_SCATTERING_NULL__;

/*******************************************************************************
 * Helper function
 ******************************************************************************/
static int
transmissivity_hit_filter
  (const struct svx_hit* hit,
   const double org[3],
   const double dir[3],
   const double range[2],
   void* context)
{
  atrstm_radcoefs_svx_T radcoefs_svx;
  struct atrstm_fetch_radcoefs_args fetch_raw_args =
    ATRSTM_FETCH_RADCOEFS_ARGS_DEFAULT;
  struct atrstm_fetch_radcoefs_svx_voxel_args fetch_svx_args =
    ATRSTM_FETCH_RADCOEFS_SVX_VOXEL_ARGS_DEFAULT;
  struct transmissivity_context* ctx = context;
  double k_min = 0;
  double k_max = 0;
  double traversal_dst = 0;
  int pursue_traversal = 1;
  ASSERT(hit && org && dir && range && ctx);
  (void)range;

  /* Fetch the K<min|max> of the current traversed voxel */
  fetch_svx_args.voxel = hit->voxel;
  fetch_svx_args.radcoefs_mask = BIT(ctx->radcoef);
  fetch_svx_args.components_mask = ATRSTM_CPNTS_MASK_ALL;
  fetch_svx_args.operations_mask = ATRSTM_SVX_OP_FLAG_MIN | ATRSTM_SVX_OP_FLAG_MAX;
  ATRSTM(fetch_radcoefs_svx_voxel(ctx->medium, &fetch_svx_args, radcoefs_svx));

  k_min = radcoefs_svx[ctx->radcoef][ATRSTM_SVX_OP_MIN];
  k_max = radcoefs_svx[ctx->radcoef][ATRSTM_SVX_OP_MAX];

  /* Setup the constants of the 'fetch' function for the current voxel */
  fetch_raw_args.wavelength = ctx->wavelength;
  fetch_raw_args.radcoefs_mask = BIT(ctx->radcoef);
  fetch_raw_args.components_mask = ATRSTM_CPNTS_MASK_ALL;
  fetch_raw_args.k_min[ctx->radcoef] = k_min;
  fetch_raw_args.k_max[ctx->radcoef] = k_max;

  /* Initialised the already traversed distance to the distance from which the
   * current ray enters into the current voxel */
  traversal_dst = hit->distance[0];

  for(;;) {
    atrstm_radcoefs_T radcoefs;
    double collision_dst;
    double tau;
    double k;
    double r;

    /* Compute the distance up to the next collision to challenge */
    tau = ssp_ran_exp(ctx->rng, 1);
    collision_dst = tau / k_max;

    /* Compute the traversed distance up to the challenged collision */
    traversal_dst += collision_dst;

    /* The collision to challenge lies behind the current voxel */
    if(traversal_dst > hit->distance[1]) {
      pursue_traversal = 1;
      break;
    }
    ASSERT(traversal_dst >= hit->distance[0]);

    /* Compute the world space position where a collision may occur */
    fetch_raw_args.pos[0] = org[0] + traversal_dst * dir[0];
    fetch_raw_args.pos[1] = org[1] + traversal_dst * dir[1];
    fetch_raw_args.pos[2] = org[2] + traversal_dst * dir[2];

    /* Fetch the radiative coefficient at the collision position */
    ATRSTM(fetch_radcoefs(ctx->medium, &fetch_raw_args, radcoefs));
    k = radcoefs[ctx->radcoef];

    r = ssp_rng_canonical(ctx->rng);
    if(r < k/k_max) { /* Real collision => stop the traversal */
      pursue_traversal = 0;
      break;
    }
  }
  return pursue_traversal;
}

static double
transmissivity
  (struct htrdr_combustion* cmd,
   struct ssp_rng* rng,
   const enum atrstm_radcoef radcoef,
   const double wlen, /* In nanometer */
   const double pos[3],
   const double dir[3],
   const double range[2])
{
  struct atrstm_trace_ray_args args = ATRSTM_TRACE_RAY_ARGS_DEFAULT;
  struct transmissivity_context transmissivity_ctx = TRANSMISSIVITY_CONTEXT_NULL;
  struct svx_hit svx_hit = SVX_HIT_NULL;
  ASSERT(cmd && rng && wlen > 0 && pos && dir && range);

  /* Degenerated range => no attenuation along dir */
  if(range[1] <= range[0]) return 1;

  /* Setup the trace ray context */
  transmissivity_ctx.rng = rng;
  transmissivity_ctx.medium = cmd->medium;
  transmissivity_ctx.wavelength = wlen;
  transmissivity_ctx.radcoef = radcoef;

  /* Setup input arguments for the ray tracing into the medium */
  d3_set(args.ray_org, pos);
  d3_set(args.ray_dir, dir);
  d2_set(args.ray_range, range);
  args.filter = transmissivity_hit_filter;
  args.context = &transmissivity_ctx;

  /* Trace the ray into the medium to compute the transmissivity */
  ATRSTM(trace_ray(cmd->medium, &args, &svx_hit));

  if(SVX_HIT_NONE(&svx_hit)) {
    return 1; /* No collision with the medium */
  } else {
    return 0; /* A collision occurs */
  }
}

static int
sample_scattering_limited_hit_filter
  (const struct svx_hit* hit,
   const double org[3],
   const double dir[3],
   const double range[2],
   void* context)
{
  atrstm_radcoefs_svx_T radcoefs_svx;
  struct atrstm_fetch_radcoefs_args fetch_raw_args =
    ATRSTM_FETCH_RADCOEFS_ARGS_DEFAULT;
  struct atrstm_fetch_radcoefs_svx_voxel_args fetch_svx_args =
    ATRSTM_FETCH_RADCOEFS_SVX_VOXEL_ARGS_DEFAULT;
  struct sample_scattering_limited_context* ctx = context;
  double ks_min = 0;
  double ks_max = 0;
  double tau_max = 0;
  double traversal_dst = 0;
  int pursue_traversal = 1;
  ASSERT(hit && org && dir && range && ctx);
  (void)range;

  /* Fetch the Ks<min|max> of the current traversed voxel */
  fetch_svx_args.voxel = hit->voxel;
  fetch_svx_args.radcoefs_mask = ATRSTM_RADCOEF_FLAG_Ks;
  fetch_svx_args.components_mask = ATRSTM_CPNTS_MASK_ALL;
  fetch_svx_args.operations_mask = ATRSTM_SVX_OP_FLAG_MIN | ATRSTM_SVX_OP_FLAG_MAX;
  ATRSTM(fetch_radcoefs_svx_voxel(ctx->medium, &fetch_svx_args, radcoefs_svx));

  ks_min = radcoefs_svx[ATRSTM_RADCOEF_Ks][ATRSTM_SVX_OP_MIN];
  ks_max = radcoefs_svx[ATRSTM_RADCOEF_Ks][ATRSTM_SVX_OP_MAX];

  /* Setup the constants of the 'fetch' function for the current voxel */
  fetch_raw_args.wavelength = ctx->wavelength;
  fetch_raw_args.radcoefs_mask = ATRSTM_RADCOEF_FLAG_Ks;
  fetch_raw_args.components_mask = ATRSTM_CPNTS_MASK_ALL;
  fetch_raw_args.k_min[ATRSTM_RADCOEF_Ks] = ks_min;
  fetch_raw_args.k_max[ATRSTM_RADCOEF_Ks] = ks_max;

  /* Initialised the already traversed distance to the distance from which the
   * current ray enters into the current voxel */
  traversal_dst = hit->distance[0];

  tau_max = (hit->distance[1] - hit->distance[0]) * ks_max;

  for(;;) {
    atrstm_radcoefs_T radcoefs;
    double collision_dst;
    double tau;
    double ks;
    double r;

    r = ssp_rng_canonical(ctx->rng);

    /* Sampled optical thickness (preferential sampling) */
    tau = -log(1.0-r*(1.0-exp(-tau_max)));
    collision_dst = tau / ks_max;

    /* Compute the traversed distance up to the challenged collision */
    traversal_dst += collision_dst;

    /* Update the ksi output data */
    ctx->ksi *= 1 - exp(-tau_max);

    /* The collision to challenge lies behind the current voxel */
    if(traversal_dst > hit->distance[1]) {
      pursue_traversal = 1;
      break;
    }
    ASSERT(traversal_dst >= hit->distance[0]);

    /* Compute the world space position where a collision may occur */
    fetch_raw_args.pos[0] = org[0] + traversal_dst * dir[0];
    fetch_raw_args.pos[1] = org[1] + traversal_dst * dir[1];
    fetch_raw_args.pos[2] = org[2] + traversal_dst * dir[2];

    /* Fetch the radiative coefficient at the collision position */
    ATRSTM(fetch_radcoefs(ctx->medium, &fetch_raw_args, radcoefs));
    ks = radcoefs[ATRSTM_RADCOEF_Ks];

    r = ssp_rng_canonical(ctx->rng);
    if(r < ks/ks_max) { /* Real collision => stop the traversal */
      ctx->traversal_dst = traversal_dst;
      pursue_traversal = 0;
      break;
    }
  }
  return pursue_traversal;
}

static void
sample_scattering_limited
  (struct htrdr_combustion* cmd,
   struct ssp_rng* rng,
   const double wlen,
   const double pos[3],
   const double dir[3],
   const double range[2],
   struct laser_sheet_scattering* sample)
{
  struct atrstm_trace_ray_args args = ATRSTM_TRACE_RAY_ARGS_DEFAULT;
  struct sample_scattering_limited_context sample_scattering_limited_ctx = 
    SAMPLE_SCATTERING_LIMITED_CONTEXT_NULL;
  struct svx_hit svx_hit = SVX_HIT_NULL;
  ASSERT(cmd && rng && wlen >= 0 && pos && dir && sample);

  /* Setup the trace ray context */
  sample_scattering_limited_ctx.rng = rng;
  sample_scattering_limited_ctx.medium = cmd->medium;
  sample_scattering_limited_ctx.wavelength = wlen;
  sample_scattering_limited_ctx.traversal_dst = 0;
  sample_scattering_limited_ctx.ksi = 1;

  /* Setup the input arguments for the ray tracing into the medium */
  d3_set(args.ray_org, pos);
  d3_set(args.ray_dir, dir);
  d2_set(args.ray_range, range);
  args.filter = sample_scattering_limited_hit_filter;
  args.context = &sample_scattering_limited_ctx;

  /* Trace the ray into the medium to compute the transmissivity */
  ATRSTM(trace_ray(cmd->medium, &args, &svx_hit));

  if(SVX_HIT_NONE(&svx_hit)) { /* No collision with the medium */
    sample->distance = INF;
    sample->ksi = 0;
  } else { /* A collision occurs */
    sample->distance = sample_scattering_limited_ctx.traversal_dst;
    sample->ksi = sample_scattering_limited_ctx.ksi;
  }
}

static double
laser_once_scattered
  (struct htrdr_combustion* cmd,
   struct ssp_rng* rng,
   const double wlen, /* In nanometer */
   const double pos[3],
   const double dir[3])
{
  #define HIT_NONE(Dst) (Dst >= DBL_MAX)

  struct laser_sheet_scattering sample = LASER_SHEET_SCATTERING_NULL;
  double laser_hit_dst[2] = {0, 0};
  double range[2] = {0, DBL_MAX};
  double Tr_pos_xin = 0;
  double Tr_xin_xs = 0;
  double L = 0; /* Radiance in W/m^2/sr/m */
  ASSERT(cmd && rng && wlen >= 0 && pos && dir);

  /* Find the intersections along dir with the laser sheet */
  htrdr_combustion_laser_trace_ray(cmd->laser, pos, dir, range, laser_hit_dst);

  /* No intersection with the laser sheet => no laser contribution */
  if(HIT_NONE(laser_hit_dst[0])) return 0;

  /* Sample the scattering position into the laser sheet */
  range[0] = laser_hit_dst[0];
  range[1] = laser_hit_dst[1];
  sample_scattering_limited(cmd, rng, wlen, pos, dir, range, &sample);

  /* No scattering in the laser sheet => no laser contribution */
  if(HIT_NONE(sample.distance)) return 0;

  /* Compute the transmissivity due to absorption from 'xin' to 'xs' */
  range[0] = laser_hit_dst[0];
  range[1] = sample.distance;
  Tr_xin_xs = transmissivity(cmd, rng, ATRSTM_RADCOEF_Ka, wlen, pos, dir, range);

  /* Compute the transmissivity from 'pos' to 'xin' */
  range[0] = 0;
  range[1] = laser_hit_dst[0];
  Tr_pos_xin = transmissivity(cmd, rng, ATRSTM_RADCOEF_Kext, wlen, pos, dir, range);

  (void)Tr_xin_xs, (void)Tr_pos_xin;

  return L;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
extern LOCAL_SYM double
combustion_compute_radiance_sw
  (struct htrdr_combustion* cmd,
   const size_t ithread,
   struct ssp_rng* rng,
   const double pos_in[3],
   const double dir_in[3],
   const double wlen) /* In nanometer */
{
  double pos[3];
  double dir[3];
  double range[2];
  double ksi; /* Throughput */
  double weight;
  ASSERT(cmd && rng && pos_in && dir_in);
  (void)cmd, (void)ithread, (void)rng, (void)pos_in, (void)dir_in, (void)wlen;
  (void)ksi, (void)laser_once_scattered;

  d3_set(pos, pos_in);
  d3_set(dir, dir_in);
  d2(range, 0, FLT_MAX);

  ksi = 1;
  weight = 0;

  /* TODO Radiative random walk */

  return weight;
}

