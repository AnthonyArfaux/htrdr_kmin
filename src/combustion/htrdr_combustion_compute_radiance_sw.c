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

#include "core/htrdr.h"

#include <astoria/atrstm.h>

#include <star/ssf.h>
#include <star/ssp.h>

#include <rsys/double2.h>
#include <rsys/double3.h>
#include <rsys/double4.h>
#include <rsys/dynamic_array.h>

/* Define a position along the ray into the semi-transparent medium */
struct position {
  double distance; /* Distance up to the position  */
  struct suvm_primitive prim; /* Volumetric primitive of the position */
  double bcoords[4]; /* Local coordinate of the position into `prim' */
};
#define POSITION_NULL__ {                                                      \
  DBL_MAX, /* Distance */                                                      \
  SUVM_PRIMITIVE_NULL__, /* Primitive */                                       \
  {0, 0, 0, 0} /* Barycentric coordinates */                                   \
}
static const struct position POSITION_NULL = POSITION_NULL__;

/* Syntactic sugar to check if the position is valid */
#define POSITION_NONE(Pos) ((Pos)->distance >= DBL_MAX)

/* Common position but preferentially sampled within a limited range. Its
 * associated ksi variable defines the correction of the weight due to the
 * normalization of the sampling pdf, and the recursivity associated with the
 * null-collision technique. */
struct position_limited {
  struct position position;
  double ksi;
};
#define POSITION_LIMITED_NULL__ {POSITION_NULL__, 1}
static const struct position_limited POSITION_LIMITED_NULL =
  POSITION_LIMITED_NULL__;

struct sample_position_context {
  /* Input data */
  struct ssp_rng* rng; /* Random Number Generator */
  struct atrstm* medium; /* Semi-transparent medium */
  double wavelength; /* Wavelength to handel in nanometer */
  enum atrstm_radcoef radcoef; /* Radiative coef used to sample a position */
  double Ts; /* Sampled optical thickness */

  /* Output data */
  struct position position;
};
#define SAMPLE_POSITION_CONTEXT_NULL__ {                                       \
  NULL, /* RNG */                                                              \
  NULL, /* Medium */                                                           \
  0, /* Wavelength */                                                          \
  ATRSTM_RADCOEFS_COUNT__, /* Radiative coefficient */                         \
  0, /* Optical thickness */                                                   \
  POSITION_NULL__                                                              \
}
static const struct sample_position_context SAMPLE_POSITION_CONTEXT_NULL =
  SAMPLE_POSITION_CONTEXT_NULL__;

struct sample_scattering_limited_context {
  /* Input data */
  struct ssp_rng* rng; /* Random number generator to use */
  struct atrstm* medium; /* Semi transparent medium */
  double wavelength; /* Wavelength to handle in nanometer */

  /* Local parameters update during ray traversal */
  double ks_2hat; /* Smallest scattering upper-field over the ray range */
  double Tmax; /* Scattering optical thickness computed using ks_2hat */
  double Ume; /* Normalization of the pdf for the sampled optical thickness */
  double sampled_vox_collision_dst; /* Scattering path length */

  /* Output data */
  struct position_limited position_limited;
};
#define SAMPLE_SCATTERING_LIMITED_CONTEXT_NULL__ {                             \
  NULL, /* RNG */                                                              \
  NULL, /* Medium */                                                           \
  -1, /* Wavelength */                                                         \
  -1, /* ks_2hat */                                                            \
  -1, /* Tau max */                                                            \
  -1, /* Ume */                                                                \
  -1, /* Sampled collision dst */                                              \
  POSITION_LIMITED_NULL__                                                      \
}
static const struct sample_scattering_limited_context
SAMPLE_SCATTERING_LIMITED_CONTEXT_NULL =
  SAMPLE_SCATTERING_LIMITED_CONTEXT_NULL__;

/*******************************************************************************
 * Sample a position along a ray into the inhomogeneous medium for a given
 * radiative coefficient
 ******************************************************************************/
static int
sample_position_hit_filter
  (const struct svx_hit* hit,
   const double org[3],
   const double dir[3],
   const double range[2],
   void* context)
{
  atrstm_radcoefs_svx_T radcoefs_svx;
  struct atrstm_fetch_radcoefs_args fetch_raw_args;
  struct atrstm_fetch_radcoefs_svx_voxel_args fetch_svx_args;
  struct sample_position_context* ctx = context;
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
    /* Compute the optical thickness of the current leaf */
    const double vox_dst = hit->distance[1] - traversal_dst;
    const double T = vox_dst * k_max;

    /* A collision occurs behind the voxel */
    if(ctx->Ts > T) {
      ctx->Ts -= T;
      pursue_traversal = 1;
      break;

    /* A collision occurs _in_ the voxel */
    } else {
      const double vox_collision_dst = ctx->Ts / k_max;
      atrstm_radcoefs_T radcoefs;
      double x[3];
      double k;
      double r;

      /* Compute the traversed distance up to the challenged collision */
      traversal_dst += vox_collision_dst;

      /* Compute the world space position where a collision may occur */
      x[0] = org[0] + traversal_dst * dir[0];
      x[1] = org[1] + traversal_dst * dir[1];
      x[2] = org[2] + traversal_dst * dir[2];

      /* Fetch the radiative coefficient at the collision position */
      ATRSTM(at(ctx->medium, x, &fetch_raw_args.prim, fetch_raw_args.bcoords));
      if(SUVM_PRIMITIVE_NONE(&fetch_raw_args.prim)) {
        k = 0;
      } else {
        ATRSTM(fetch_radcoefs(ctx->medium, &fetch_raw_args, radcoefs));
        k = radcoefs[ctx->radcoef];
      }

      r = ssp_rng_canonical(ctx->rng);

      if(r > k/k_max) { /* Null collision */
        ctx->Ts = ssp_ran_exp(ctx->rng, 1);
      } else { /* Real collision */
        /* Setup output data */
        ctx->position.distance = traversal_dst;
        ctx->position.prim = fetch_raw_args.prim;
        d4_set(ctx->position.bcoords, fetch_raw_args.bcoords);

        /* Stop ray traversal */
        pursue_traversal = 0;
        break;
      }
    }
  }
  return pursue_traversal;
}

static void
sample_position
  (struct htrdr_combustion* cmd,
   struct ssp_rng* rng,
   const enum atrstm_radcoef radcoef,
   const double pos[3],
   const double dir[3],
   const double range[2],
   struct position* position)
{
  struct atrstm_trace_ray_args args = ATRSTM_TRACE_RAY_ARGS_DEFAULT;
  struct sample_position_context sample_pos_ctx = SAMPLE_POSITION_CONTEXT_NULL;
  struct svx_hit svx_hit = SVX_HIT_NULL;
  double wlen = 0;
  ASSERT(cmd && rng && pos && dir && range);

  wlen = htrdr_combustion_laser_get_wavelength(cmd->laser);

  /* Sample an optical thickness */
  sample_pos_ctx.Ts = ssp_ran_exp(rng, 1);

  /* Setup the arguments of the function invoked per voxel (i.e. the filter
   * function) */
  sample_pos_ctx.rng = rng;
  sample_pos_ctx.medium = cmd->medium;
  sample_pos_ctx.wavelength = wlen;
  sample_pos_ctx.radcoef = radcoef;

  /* Setup ray tracing arguments */
  d3_set(args.ray_org, pos);
  d3_set(args.ray_dir, dir);
  d2_set(args.ray_range, range);
  args.filter = sample_position_hit_filter;
  args.context = &sample_pos_ctx;

  /* Trace the ray into the heterogeneous medium */
  ATRSTM(trace_ray(cmd->medium, &args, &svx_hit));

  if(SVX_HIT_NONE(&svx_hit)) {
    /* No collision with the medium was found */
    *position = POSITION_NULL;
  } else {
    /* A collision occurs into the medium */
    *position = sample_pos_ctx.position;
  }
}

/*******************************************************************************
 * Preferentially sample a scattering position into an inhomogeneous medium
 * according to a limited range
 ******************************************************************************/
/* Find the constant Ks max (named ks_2hat) along the traced ray */
static int
sample_scattering_limited_find_ks_2hat
  (const struct svx_hit* hit,
   const double org[3],
   const double dir[3],
   const double range[2],
   void* context)
{
  struct sample_scattering_limited_context* ctx = context;
  struct atrstm_fetch_radcoefs_svx_voxel_args fetch_svx_args;
  atrstm_radcoefs_svx_T radcoefs_svx;
  ASSERT(hit && org && dir && range && context);
  (void)org, (void)dir;

  /* In all situations, initialise ks_2hat with the ks_max of the root node */
  if(hit->voxel.depth == 0) {
    fetch_svx_args.voxel = hit->voxel;
    fetch_svx_args.radcoefs_mask = ATRSTM_RADCOEF_FLAG_Ks;
    fetch_svx_args.components_mask = ATRSTM_CPNTS_MASK_ALL;
    fetch_svx_args.operations_mask = ATRSTM_SVX_OP_FLAG_MAX;
    ATRSTM(fetch_radcoefs_svx_voxel(ctx->medium, &fetch_svx_args, radcoefs_svx));
    ctx->ks_2hat = radcoefs_svx[ATRSTM_RADCOEF_Ks][ATRSTM_SVX_OP_MAX];

  /* Update ks_2hat with the ks_max of the current descending node until the ray
   * range was no more included into this node */
  } else if(hit->distance[0] <= range[0] && range[1] <= hit->distance[1]) {
    fetch_svx_args.voxel = hit->voxel;
    fetch_svx_args.radcoefs_mask = ATRSTM_RADCOEF_FLAG_Ks;
    fetch_svx_args.components_mask = ATRSTM_CPNTS_MASK_ALL;
    fetch_svx_args.operations_mask = ATRSTM_SVX_OP_FLAG_MAX;
    ATRSTM(fetch_radcoefs_svx_voxel(ctx->medium, &fetch_svx_args, radcoefs_svx));
    ctx->ks_2hat = radcoefs_svx[ATRSTM_RADCOEF_Ks][ATRSTM_SVX_OP_MAX];
  }

  /* Do not stop here: keep diving up to the leaves */
  return 0;
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
  struct atrstm_fetch_radcoefs_args fetch_raw_args;
  struct atrstm_fetch_radcoefs_svx_voxel_args fetch_svx_args;
  struct sample_scattering_limited_context* ctx = context;
  double ks_min = 0;
  double ks_max = 0;
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
  ASSERT(ks_max <= ctx->ks_2hat);

  /* Setup the constants of the 'fetch' function for the current voxel */
  fetch_raw_args.wavelength = ctx->wavelength;
  fetch_raw_args.radcoefs_mask = ATRSTM_RADCOEF_FLAG_Ks;
  fetch_raw_args.components_mask = ATRSTM_CPNTS_MASK_ALL;
  fetch_raw_args.k_min[ATRSTM_RADCOEF_Ks] = ks_min;
  fetch_raw_args.k_max[ATRSTM_RADCOEF_Ks] = ks_max;

  /* Initialised the already traversed distance to the distance from which the
   * current ray enters into the current voxel */
  traversal_dst = hit->distance[0];

  /* First traversed leaf */
  if(ctx->sampled_vox_collision_dst < 0) {
    ctx->Tmax = (range[1] - range[0]) * ctx->ks_2hat;
    ctx->Ume = 1 - exp(-ctx->Tmax);
  }

  for(;;) {
    atrstm_radcoefs_T radcoefs;
    double vox_dst;
    double tau;
    double ks;
    double r;

    /* Compute the remaining distance to traverse in the current voxel */
    vox_dst = hit->distance[1] - traversal_dst;

    /* A collision distance was not already sampled */
    if(ctx->sampled_vox_collision_dst < 0) {
      r = ssp_rng_canonical(ctx->rng);
      tau = -log(1.0-r*(1.0-exp(-ctx->Tmax)));
      ctx->sampled_vox_collision_dst = tau / ctx->ks_2hat;

      /* Update the ksi output data */
      ctx->position_limited.ksi *= ctx->Ume;
    }

    /* Compute the traversed distance up to the challenged collision */
    traversal_dst = traversal_dst + ctx->sampled_vox_collision_dst;

    /* The collision to challenge lies behind the current voxel */
    if(traversal_dst > hit->distance[1]) {
      ctx->sampled_vox_collision_dst -= vox_dst;
      pursue_traversal = 1;
      break;
    }
    ASSERT(traversal_dst >= hit->distance[0]);

    r = ssp_rng_canonical(ctx->rng);
    if(r >= ks_max / ctx->ks_2hat) { /* Null collision */
      ctx->sampled_vox_collision_dst = -1;
    } else { /* Collision with ks_max */
      double x[3];

      /* Compute the world space position where a collision may occur */
      x[0] = org[0] + traversal_dst * dir[0];
      x[1] = org[1] + traversal_dst * dir[1];
      x[2] = org[2] + traversal_dst * dir[2];

      /* Fetch the radiative coefficient at the collision position */
      ATRSTM(at(ctx->medium, x, &fetch_raw_args.prim, fetch_raw_args.bcoords));
      if(SUVM_PRIMITIVE_NONE(&fetch_raw_args.prim)) {
        ks = 0;
      } else {
        ATRSTM(fetch_radcoefs(ctx->medium, &fetch_raw_args, radcoefs));
        ks = radcoefs[ATRSTM_RADCOEF_Ks];
      }

      r = ssp_rng_canonical(ctx->rng);

      if(r >= ks/ks_max) { /* Null collision */
        ctx->sampled_vox_collision_dst = -1;
      } else { /* Real collision */
        /* Setup output data */
        ctx->position_limited.position.distance = traversal_dst;
        ctx->position_limited.position.prim = fetch_raw_args.prim;
        d4_set(ctx->position_limited.position.bcoords, fetch_raw_args.bcoords);

        /* Stop ray traversal */
        pursue_traversal = 0;
        break;
      }
    }
  }
  return pursue_traversal;
}

static void
sample_scattering_limited
  (struct htrdr_combustion* cmd,
   struct ssp_rng* rng,
   const double pos[3],
   const double dir[3],
   const double range[2],
   struct position_limited* position)
{
  struct atrstm_trace_ray_args args = ATRSTM_TRACE_RAY_ARGS_DEFAULT;
  struct sample_scattering_limited_context sample_scattering_limited_ctx =
    SAMPLE_SCATTERING_LIMITED_CONTEXT_NULL;
  struct svx_hit svx_hit = SVX_HIT_NULL;
  double wlen = 0;
  ASSERT(cmd && rng && pos && dir && position);

  wlen = htrdr_combustion_laser_get_wavelength(cmd->laser);

  /* Setup the trace ray context */
  sample_scattering_limited_ctx.rng = rng;
  sample_scattering_limited_ctx.medium = cmd->medium;
  sample_scattering_limited_ctx.wavelength = wlen;
  sample_scattering_limited_ctx.ks_2hat = 0;
  sample_scattering_limited_ctx.Tmax = -1;
  sample_scattering_limited_ctx.Ume = -1;
  sample_scattering_limited_ctx.sampled_vox_collision_dst = -1;
  sample_scattering_limited_ctx.position_limited = POSITION_LIMITED_NULL;

  /* Setup the input arguments for the ray tracing into the medium */
  d3_set(args.ray_org, pos);
  d3_set(args.ray_dir, dir);
  d2_set(args.ray_range, range);
  args.challenge = sample_scattering_limited_find_ks_2hat;
  args.filter = sample_scattering_limited_hit_filter;
  args.context = &sample_scattering_limited_ctx;

  /* Trace the ray into the medium to compute the transmissivity */
  ATRSTM(trace_ray(cmd->medium, &args, &svx_hit));

  if(SVX_HIT_NONE(&svx_hit)) {
    /* No scattering event was found */
    *position = POSITION_LIMITED_NULL;
  } else {
    /* A scattering event occurs into the medium */
    *position = sample_scattering_limited_ctx.position_limited;
  }
}

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static double
transmissivity
  (struct htrdr_combustion* cmd,
   struct ssp_rng* rng,
   const enum atrstm_radcoef radcoef,
   const double pos[3],
   const double dir[3],
   const double range[2])
{
  struct position position = POSITION_NULL;

  /* Degenerated range => no attenuation along dir */
  if(range[1] <= range[0]) return 1;

  sample_position(cmd, rng, radcoef, pos, dir, range, &position);

  if(POSITION_NONE(&position)) {
    return 1; /* No collision with the medium */
  } else {
    return 0; /* A collision occurs */
  }
}

/* This function computes the contribution of the in-laser scattering for the
 * current scattering position of the reverse optical path 'pos' and current
 * direction 'dir'. One contribution has to be computed for each scattering
 * position in order to compute the value of the Monte-Carlo weight for the
 * optical path (weight of the statistical realization). */
static double
laser_once_scattered
  (struct htrdr_combustion* cmd,
   const size_t ithread,
   struct ssp_rng* rng,
   const double pos[3],
   const double dir[3])
{
  /* RDG-FA phase function */
  struct atrstm_rdgfa rdgfa_param = ATRSTM_RDGFA_NULL;
  struct atrstm_fetch_rdgfa_args fetch_rdgfa_args =
    ATRSTM_FETCH_RDGFA_ARGS_DEFAULT;
  struct ssf_phase_rdgfa_setup_args setup_rdgfa_args =
    SSF_PHASE_RDGFA_SETUP_ARGS_DEFAULT;
  struct ssf_phase* phase = NULL;

  /* Scattering position into the laser sheet */
  struct position_limited sc_sample = POSITION_LIMITED_NULL;
  double xsc[3] = {0, 0, 0}; /* World space position */

  /* The transmissivities to compute */
  double Tr_ext_pos_xin = 0;
  double Tr_ext_xsc_lse = 0;
  double Tr_abs_xin_xsc = 0;

  /* Laser data */
  double laser_hit_dst[2] = {0, 0};
  double laser_flux_density = 0; /* In W/m^2 */
  double wlen = 0; /* In nm */

  /* Miscellaneous varbale */
  double wi[3] = {0, 0, 0};
  double wo[3] = {0, 0, 0};
  double range[2] = {0, DBL_MAX};
  double R = 0;

  /* Radiance to compute in W/m^2/sr */
  double L = 0;

  ASSERT(cmd && rng && pos && dir);

  /* Fetch the laser wavelength */
  wlen = htrdr_combustion_laser_get_wavelength(cmd->laser);

  /* Find the intersections along dir with the laser sheet */
  htrdr_combustion_laser_trace_ray(cmd->laser, pos, dir, range, laser_hit_dst);

  /* No intersection with the laser sheet => no laser contribution */
  if(HTRDR_COMBUSTION_LASER_HIT_NONE(laser_hit_dst)) return 0;

  /* Compute the transmissivity from 'pos' to 'xin' */
  range[0] = 0;
  range[1] = laser_hit_dst[0];
  Tr_ext_pos_xin = transmissivity(cmd, rng, ATRSTM_RADCOEF_Kext, pos, dir, range);
  if(Tr_ext_pos_xin == 0) return 0; /* no laser contribution */

  /* Sample the scattering position into the laser sheet */
  range[0] = laser_hit_dst[0];
  range[1] = laser_hit_dst[1];
  sample_scattering_limited(cmd, rng, pos, dir, range, &sc_sample);
  if(POSITION_NONE(&sc_sample.position)) return 0; /* No laser contribution */
  ASSERT(laser_hit_dst[0] <= sc_sample.position.distance);
  ASSERT(laser_hit_dst[1] >= sc_sample.position.distance);

  /* Compute the transmissivity from 'xin' to 'xsc' */
  range[0] = laser_hit_dst[0]; /* <=> xin */
  range[1] = sc_sample.position.distance; /* <=> xsc */
  Tr_abs_xin_xsc = transmissivity(cmd, rng, ATRSTM_RADCOEF_Ka, pos, dir, range);
  if(Tr_abs_xin_xsc == 0) return 0; /* No laser contribution */

  /* Compute the scattering position */
  xsc[0] = pos[0] + dir[0] * sc_sample.position.distance;
  xsc[1] = pos[1] + dir[1] * sc_sample.position.distance;
  xsc[2] = pos[2] + dir[2] * sc_sample.position.distance;

  /* Retrieve the direction toward the laser surface emission */
  htrdr_combustion_laser_get_direction(cmd->laser, wi);
  d3_minus(wi, wi);

  /* Compute the transmissivity from xsc to the laser surface emission */
  range[0] = 0;
  range[1] = htrdr_combustion_laser_compute_surface_plane_distance
    (cmd->laser, xsc);
  ASSERT(range[1] >= 0);
  Tr_ext_xsc_lse = transmissivity(cmd, rng, ATRSTM_RADCOEF_Kext, xsc, wi, range);
  if(Tr_ext_xsc_lse == 0) return 0; /* No laser contribution */

  /* Retrieve the RDG-FA phase function parameters from the semi transparent
   * medium */
  fetch_rdgfa_args.wavelength = wlen;
  fetch_rdgfa_args.prim = sc_sample.position.prim;
  fetch_rdgfa_args.bcoords[0] = sc_sample.position.bcoords[0];
  fetch_rdgfa_args.bcoords[1] = sc_sample.position.bcoords[1];
  fetch_rdgfa_args.bcoords[2] = sc_sample.position.bcoords[2];
  fetch_rdgfa_args.bcoords[3] = sc_sample.position.bcoords[3];
  ATRSTM(fetch_rdgfa(cmd->medium, &fetch_rdgfa_args, &rdgfa_param));

  /* Setup the RDG-FA phase function */
  phase = cmd->rdgfa_phase_functions[ithread];
  setup_rdgfa_args.wavelength = rdgfa_param.wavelength;
  setup_rdgfa_args.fractal_dimension = rdgfa_param.fractal_dimension;
  setup_rdgfa_args.gyration_radius = rdgfa_param.gyration_radius;
  SSF(phase_rdgfa_setup(phase, &setup_rdgfa_args));

  /* Evaluate the phase function at the scattering position */
  d3_minus(wo, dir); /* Ensure SSF convention */
  R = ssf_phase_eval(phase, wo, wi); /* In sr^-1 */

  laser_flux_density = htrdr_combustion_laser_get_flux_density(cmd->laser);

  L = Tr_ext_pos_xin
    * Tr_abs_xin_xsc
    * sc_sample.ksi * R
    * Tr_ext_xsc_lse
    * laser_flux_density;
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
   const double dir_in[3])
{
  /* RDG-FA phase function */
  struct atrstm_rdgfa rdgfa_param = ATRSTM_RDGFA_NULL;
  struct atrstm_fetch_rdgfa_args fetch_rdgfa_args =
    ATRSTM_FETCH_RDGFA_ARGS_DEFAULT;
  struct ssf_phase_rdgfa_setup_args setup_rdgfa_args =
    SSF_PHASE_RDGFA_SETUP_ARGS_DEFAULT;
  struct ssf_phase* phase = NULL;

  /* Transmissivity between the probe position (i.e. 'pos_in') and the current
   * scattering position over the reverse scattering path */
  double Tr_abs = 1;

  /* Monte carlo weight of the simulated optical path */
  double weight = 0;

  /* Miscellaneous variables */
  double pos[3] = {0, 0, 0};
  double dir[3] = {0, 0, 0};
  double range[2] = {0, DBL_MAX};
  double wlen = 0;
  ASSERT(cmd && rng && pos_in && dir_in);

  d3_set(pos, pos_in);
  d3_set(dir, dir_in);
  d2(range, 0, FLT_MAX);

  wlen = htrdr_combustion_laser_get_wavelength(cmd->laser);

  Tr_abs = 1;
  weight = 0;

  for(;;) {
    struct position scattering = POSITION_NULL;
    double laser_sheet_emissivity = 0; /* In W/m^2/sr */
    double Tr_abs_pos_xsc = 0;
    double wo[3];
    double wi[3];

    /* Handle the laser contribution */
    laser_sheet_emissivity = laser_once_scattered(cmd, ithread, rng, pos, dir);
    weight += Tr_abs * laser_sheet_emissivity;

    /* Sample a scattering position */
    sample_position(cmd, rng, ATRSTM_RADCOEF_Ks, pos, dir, range, &scattering);
    if(POSITION_NONE(&scattering)) break;

    /* Compute the absorption transmissivity */
    range[0] = 0;
    range[1] = scattering.distance;
    Tr_abs_pos_xsc = transmissivity(cmd, rng, ATRSTM_RADCOEF_Ka, pos, dir, range);
    if(Tr_abs_pos_xsc == 0) break;

    /* Update the overall absorption transmissivity of the optical path */
    Tr_abs *= Tr_abs_pos_xsc;

    /* Update the position of the optical path */
    pos[0] = pos[0] + dir[0] * scattering.distance;
    pos[1] = pos[1] + dir[1] * scattering.distance;
    pos[2] = pos[2] + dir[2] * scattering.distance;

    /* Retrieve the RDG-FA phase function parameters */
    fetch_rdgfa_args.wavelength = wlen;
    fetch_rdgfa_args.prim = scattering.prim;
    fetch_rdgfa_args.bcoords[0] = scattering.bcoords[0];
    fetch_rdgfa_args.bcoords[1] = scattering.bcoords[1];
    fetch_rdgfa_args.bcoords[2] = scattering.bcoords[2];
    fetch_rdgfa_args.bcoords[3] = scattering.bcoords[3];
    ATRSTM(fetch_rdgfa(cmd->medium, &fetch_rdgfa_args, &rdgfa_param));

    /* Setup the RDG-FA phase function corresponding to the scattering event */
    phase = cmd->rdgfa_phase_functions[ithread];
    setup_rdgfa_args.wavelength = rdgfa_param.wavelength;
    setup_rdgfa_args.fractal_dimension = rdgfa_param.fractal_dimension;
    setup_rdgfa_args.gyration_radius = rdgfa_param.gyration_radius;
    SSF(phase_rdgfa_setup(phase, &setup_rdgfa_args));

    /* Sample a new optical path direction from the phase function */
    d3_minus(wo, dir); /* Ensure SSF convention */
    ssf_phase_sample(phase, rng, wo, wi, NULL);

    /* Update the optical path direction */
    dir[0] = wi[0];
    dir[1] = wi[1];
    dir[2] = wi[2];
  }
  return weight;
}
