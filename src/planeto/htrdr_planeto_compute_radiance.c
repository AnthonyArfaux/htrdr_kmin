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

#include "planeto/htrdr_planeto_c.h"
#include "planeto/htrdr_planeto_source.h"

#include <rad-net/rnatm.h>
#include <rad-net/rngrd.h>

#include <star/s3d.h>
#include <star/ssf.h>
#include <star/ssp.h>
#include <star/suvm.h>
#include <star/svx.h>

#include <rsys/double2.h>
#include <rsys/double3.h>

/* Syntactic sugar */
#define DISTANCE_NONE(Dst) ((Dst) >= FLT_MAX)
#define SURFACE_SCATTERING(Scattering) (!S3D_HIT_NONE(&(Scattering)->hit))

struct scattering {
  /* Set to S3D_HIT_NULL if the scattering occurs in volume.*/
  struct s3d_hit hit;

  /* The surface normal defines only if scattering is on the surface. It is
   * normalized and looks towards the incoming direction */
  double normal[3];

  double distance;
};
#define SCATTERING_NULL__ {S3D_HIT_NULL__, {0,0,0}, DBL_MAX}
static const struct scattering SCATTERING_NULL = SCATTERING_NULL__;

/* Arguments of the filtering function used to sample a position */
struct sample_distance_context {
  struct ssp_rng* rng;
  struct rnatm* atmosphere;
  size_t iband;
  size_t iquad;
  double wavelength; /* In nm */
  enum rnatm_radcoef radcoef;
  double Ts; /* Sample optical thickness */

  /* Output data */
  double distance;
};
#define SAMPLE_DISTANCE_CONTEXT_NULL__ {                                       \
  NULL, NULL, 0, 0, 0, RNATM_RADCOEFS_COUNT__, 0, DBL_MAX                      \
}
static const struct sample_distance_context SAMPLE_DISTANCE_CONTEXT_NULL =
  SAMPLE_DISTANCE_CONTEXT_NULL__;

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static INLINE res_T
check_planeto_compute_radiance_args
  (const struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args)
{
  struct rnatm_band_desc band = RNATM_BAND_DESC_NULL;
  res_T res = RES_OK;

  if(!args || !args->rng)
    return RES_BAD_ARG;

  /* Invalid thread index */
  if(args->ithread >= htrdr_get_threads_count(cmd->htrdr))
    return RES_BAD_ARG;

  /* Invalid input direction */
  if(!d3_is_normalized(args->path_dir))
    return RES_BAD_ARG;

  /* Invalid wavelength */
  if(args->wlen < cmd->spectral_domain.wlen_range[0]
  || args->wlen > cmd->spectral_domain.wlen_range[1])
    return RES_BAD_ARG;

  res = rnatm_band_get_desc(cmd->atmosphere, args->iband, &band);
  if(res != RES_OK) return res;

  /* Inconsistent spectral dimension */
  if(args->wlen <  band.lower
  || args->wlen >= band.upper /* Exclusive */
  || args->iquad >= band.quad_pts_count)
    return RES_BAD_ARG;

  return RES_OK;
}

static int
sample_position_hit_filter
  (const struct svx_hit* hit,
   const double org[3],
   const double dir[3],
   const double range[2],
   void* context)
{
  struct rnatm_get_radcoef_args get_k_args = RNATM_GET_RADCOEF_ARGS_NULL;
  struct sample_distance_context* ctx = context;
  double k_min = 0;
  double k_max = 0;
  double dst_travelled = 0;
  int pursue_traversal = 1;
  ASSERT(hit && org && range && context);
  (void)range;

  dst_travelled = hit->distance[0];

  /* Get the k_min, k_max of the voxel */
  k_min = rnatm_get_k_svx_voxel
    (ctx->atmosphere, &hit->voxel, ctx->radcoef, RNATM_SVX_OP_MIN);
  k_max = rnatm_get_k_svx_voxel
    (ctx->atmosphere, &hit->voxel, ctx->radcoef, RNATM_SVX_OP_MAX);

  /* Configure the common input arguments to retrieve the radiative coefficient
   * of a given position */
  get_k_args.iband = ctx->iband;
  get_k_args.iquad = ctx->iquad;
  get_k_args.radcoef = ctx->radcoef;
  get_k_args.k_min = k_min;
  get_k_args.k_max = k_max;

  for(;;) {
    /* Compute the optical thickness of the voxel */
    const double vox_dst = hit->distance[1] - dst_travelled;
    const double T = vox_dst * k_max;

    /* A collision occurs behind the voxel */
    if(ctx->Ts > T) {
      ctx->Ts -= T;
      pursue_traversal = 1;
      break;

    /* A collision occurs in the voxel */
    } else {
      const double vox_dst_collision = ctx->Ts / k_max;
      double  k, r;

      /* Calculate the distance travelled to the collision to be queried */
      dst_travelled += vox_dst_collision;

      /* Retrieve the radiative coefficient at the collision position */
      get_k_args.pos[0] = org[0] + dst_travelled * dir[0];
      get_k_args.pos[1] = org[1] + dst_travelled * dir[1];
      get_k_args.pos[2] = org[2] + dst_travelled * dir[2];
      RNATM(get_radcoef(ctx->atmosphere, &get_k_args, &k));

      r = ssp_rng_canonical(ctx->rng);

      /* Null collision */
      if(r > k/k_max) {
        ctx->Ts = ssp_ran_exp(ctx->rng, 1);

      /* Real collision */
      } else {
        ctx->distance = dst_travelled;
        pursue_traversal = 0;
        break;
      }
    }
  }

  return pursue_traversal;
}

static double
sample_distance
  (struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args,
   const enum rnatm_radcoef radcoef,
   const double pos[3],
   const double dir[3],
   const double range[2])
{
  struct rnatm_trace_ray_args rt = RNATM_TRACE_RAY_ARGS_NULL;
  struct sample_distance_context ctx = SAMPLE_DISTANCE_CONTEXT_NULL;
  struct svx_hit hit;
  ASSERT(cmd && args && pos && dir && d3_is_normalized(dir) && range);
  ASSERT((unsigned)radcoef < RNATM_RADCOEFS_COUNT__);
  ASSERT(range[0] < range[1]);

  /* Sample an optical thickness */
  ctx.Ts = ssp_ran_exp(args->rng, 1);

  /* Setup the remaining arguments of the RT context */
  ctx.rng = args->rng;
  ctx.atmosphere = cmd->atmosphere;
  ctx.iband = args->iband;
  ctx.iquad = args->iquad;
  ctx.wavelength = args->wlen;
  ctx.radcoef = radcoef;

  /* Trace the ray into the atmosphere */
  d3_set(rt.ray_org, pos);
  d3_set(rt.ray_dir, dir);
  rt.ray_range[0] = range[0];
  rt.ray_range[1] = range[1];
  rt.filter = sample_position_hit_filter;
  rt.context = &ctx;
  rt.iband = args->iband;
  rt.iquad = args->iquad;
  RNATM(trace_ray(cmd->atmosphere, &rt, &hit));

  if(SVX_HIT_NONE(&hit)) { /* No collision found */
    return INF;
  } else { /* A (real) collision occured */
    return ctx.distance;
  }
}

static INLINE double
transmissivity
  (struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args,
   const enum rnatm_radcoef radcoef,
   const double pos[3],
   const double dir[3],
   const double range_max)
{
  double range[2];
  double dst = 0;
  ASSERT(range_max >= 0);

  range[0] = 0;
  range[1] = range_max;
  dst = sample_distance(cmd, args, radcoef, pos, dir, range);

  if(DISTANCE_NONE(dst)) {
    return 1.0; /* No collision in the medium */
  } else {
    return 0.0; /* A (real) collision occurs */
  }
}

static double
direct_contribution
  (struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args,
   const double pos[3],
   const double dir[3],
   const struct s3d_hit* hit_from)
{
  struct rngrd_trace_ray_args rt = RNGRD_TRACE_RAY_ARGS_DEFAULT;
  struct s3d_hit hit;
  double Tr;
  double Ld;
  ASSERT(cmd && args && pos && dir);

  /* Is the source hidden? */
  d3_set(rt.ray_org, pos);
  d3_set(rt.ray_dir, dir);
  if(hit_from) rt.hit_from = *hit_from;
  RNGRD(trace_ray(cmd->ground, &rt, &hit));
  if(!S3D_HIT_NONE(&hit)) return 0;

  Tr = transmissivity(cmd, args, RNATM_RADCOEF_Kext, pos, dir, INF);
  Ld = htrdr_planeto_source_get_radiance(cmd->source, args->wlen);
  return Ld * Tr;
}

static void
sample_scattering
  (struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args,
   const double pos[3],
   const double dir[3],
   struct scattering* sc)
{
  struct rngrd_trace_ray_args rt = RNGRD_TRACE_RAY_ARGS_DEFAULT;
  struct s3d_hit hit;
  double range[2];
  double dst;
  ASSERT(cmd && args && pos && dir && sc);

  *sc = SCATTERING_NULL;

  /* Look for a surface intersection */
  d3_set(rt.ray_org, pos);
  d3_set(rt.ray_dir, dir);
  d2(rt.ray_range, 0, INF);
  RNGRD(trace_ray(cmd->ground, &rt, &hit));

  /* Look for an atmospheric collision */
  range[0] = 0;
  range[1] = hit.distance;
  dst = sample_distance(cmd, args, RNATM_RADCOEF_Ks, pos, dir, range);

  /* Scattering in volume */
  if(!DISTANCE_NONE(dst)) {
    sc->distance = dst;
    sc->hit = S3D_HIT_NULL;

  /* Surface scattering */
  } else if(!S3D_HIT_NONE(&hit)) {
    /* Normalize the normal and ensure that it points to `dir' */
    d3_normalize(sc->normal, d3_set_f3(sc->normal, hit.normal));
    if(d3_dot(sc->normal, dir) > 0) d3_minus(sc->normal, sc->normal);

    sc->distance = hit.distance;
    sc->hit = hit;

  /* No scattering */
  } else {
    sc->distance = INF;
    sc->hit = S3D_HIT_NULL;
  }
}

static INLINE struct ssf_bsdf*
create_bsdf
  (struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args,
   const struct s3d_hit* hit)
{
  struct rngrd_create_bsdf_args bsdf_args = RNGRD_CREATE_BSDF_ARGS_NULL;
  struct ssf_bsdf* bsdf = NULL;
  ASSERT(!S3D_HIT_NONE(hit));

  /* Retrieve the BSDF at the intersected surface position */
  bsdf_args.prim = hit->prim;
  bsdf_args.barycentric_coords[0] = hit->uv[0];
  bsdf_args.barycentric_coords[1] = hit->uv[1];
  bsdf_args.barycentric_coords[2] = 1 - hit->uv[0] - hit->uv[1];
  bsdf_args.wavelength = args->wlen;
  bsdf_args.r = ssp_rng_canonical(args->rng);
  RNGRD(create_bsdf(cmd->ground, &bsdf_args, &bsdf));

  return bsdf;
}

static INLINE struct ssf_phase*
create_phase_fn
  (struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args,
   const double pos[3]) /* Scattering position */
{
  struct rnatm_sample_component_args sample_args =
    RNATM_SAMPLE_COMPONENT_ARGS_NULL;
  struct rnatm_cell_create_phase_fn_args phase_fn_args =
    RNATM_CELL_CREATE_PHASE_FN_ARGS_NULL;
  struct ssf_phase* phase = NULL;
  size_t cpnt;
  ASSERT(cmd && args && pos);

  /* Sample the atmospheric scattering component */
  d3_set(sample_args.pos, pos);
  sample_args.iband = args->iband;
  sample_args.iquad = args->iquad;
  sample_args.radcoef = RNATM_RADCOEF_Ks;
  sample_args.r = ssp_rng_canonical(args->rng);
  RNATM(sample_component(cmd->atmosphere, &sample_args, &cpnt));

  /* Retrieve the component cell in which the scattering position is located */
  RNATM(fetch_cell(cmd->atmosphere, pos, cpnt, &phase_fn_args.cell,
    phase_fn_args.barycentric_coords));
  ASSERT(!SUVM_PRIMITIVE_NONE(&phase_fn_args.cell.prim));

  /* Retrieve the component phase function */
  phase_fn_args.wavelength = args->wlen;
  phase_fn_args.r[0] = ssp_rng_canonical(args->rng);
  phase_fn_args.r[1] = ssp_rng_canonical(args->rng);
  RNATM(cell_create_phase_fn(cmd->atmosphere, &phase_fn_args, &phase));

  return phase;
}

/* Return the direct contribution at the scattering position */
static double
surface_scattering
  (struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args,
   const struct scattering* sc,
   const double sc_pos[3], /* Scattering position */
   const double in_dir[3], /* Incident direction */
   double sc_dir[3], /* Sampled scattering direction */
   double *out_refl) /* Surface reflectivity */
{
  struct ssf_bsdf* bsdf = NULL;
  const double* N = NULL;
  double wo[3] = {0,0,0};
  double reflectivity = 0; /* Surface reflectivity */
  double L = 0;
  int mask = 0;
  ASSERT(cmd && args && sc && sc_pos && in_dir && sc_dir && out_refl);
  ASSERT(d3_dot(sc->normal, in_dir) < 0 && d3_is_normalized(sc->normal));

  bsdf = create_bsdf(cmd, args, &sc->hit);
  N = sc->normal;
  d3_minus(wo, in_dir); /* Match StarSF convention */
  ASSERT(d3_dot(wo, N) > 0);

  /* Sample the scattering direction */
  reflectivity = ssf_bsdf_sample(bsdf, args->rng, wo, N, sc_dir, &mask, NULL);

  /* Fully absorbs transmissions */
  if(mask & SSF_TRANSMISSION) reflectivity = 0;

  /* Calculate direct contribution for specular reflection */
  if((mask & SSF_SPECULAR)
  && (mask & SSF_REFLECTION)
  && htrdr_planeto_source_is_targeted(cmd->source, sc_pos, sc_dir)) {
    const double Ld = direct_contribution(cmd, args, sc_pos, sc_dir, &sc->hit);
    L = Ld * reflectivity;

  /* Calculate direct contribution in general case */
  } else if(!(mask & SSF_SPECULAR)) {
    double wi[3], pdf;

    /* Sample a direction toward the source */
    pdf = htrdr_planeto_source_sample_direction(cmd->source, args->rng, sc_pos, wi);

    /* The source is behind the surface */
    if(d3_dot(wi, N) <= 0) {
      L = 0;

    /* The source is above the surface */
    } else {
      const double Ld = direct_contribution(cmd, args, sc_pos, wi, &sc->hit);
      const double rho = ssf_bsdf_eval(bsdf, wo, N, wi);
      const double cos_N_wi = d3_dot(N, wi);
      ASSERT(cos_N_wi > 0);

      L = Ld * rho * cos_N_wi / pdf;
    }
  }

  SSF(bsdf_ref_put(bsdf));

  *out_refl = reflectivity;
  return L;
}

/* Return the direct contribution at the scattering position */
static INLINE double
volume_scattering
  (struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args,
   const double sc_pos[3], /* Scattering position */
   const double in_dir[3], /* Incident direction */
   double sc_dir[3]) /* Sampled scattering direction */
{
  struct ssf_phase* phase = NULL;
  double wo[3] = {0,0,0};
  double wi[3] = {0,0,0};
  double L = 0;
  double pdf = 0;
  double rho = 0;
  double Ld = 0;
  ASSERT(cmd && args && sc_pos && in_dir && sc_dir);

  phase = create_phase_fn(cmd, args, sc_pos);
  d3_minus(wo, in_dir); /* Match StarSF convention */

  ssf_phase_sample(phase, args->rng, wo, sc_dir, NULL);

  /* Sample a direction toward the source */
  pdf = htrdr_planeto_source_sample_direction(cmd->source, args->rng, sc_pos, wi);

  /* Calculate the direct contribution at the scattering position */
  Ld = direct_contribution(cmd, args, sc_pos, wi, NULL);
  rho = ssf_phase_eval(phase, wo, wi);
  L = Ld * rho / pdf;

  SSF(phase_ref_put(phase));
  return L;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
double
planeto_compute_radiance
  (struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args)
{
  double pos[3];
  double dir[3];
  double L = 0; /* Radiance */
  double Tr_abs = 1; /* Absorption transmissivity */
  size_t nsc = 0; /* For debug */
  ASSERT(cmd && check_planeto_compute_radiance_args(cmd, args) == RES_OK);

  d3_set(pos, args->path_org);
  d3_set(dir, args->path_dir);

  if(htrdr_planeto_source_is_targeted(cmd->source, pos, dir)) {
    L = direct_contribution(cmd, args, pos, dir, NULL);
  }

  for(;;) {
    struct scattering sc;
    double sc_pos[3];
    double sc_dir[3];

    sample_scattering(cmd, args, pos, dir, &sc);

    /* No scattering. Stop the path */
    if(DISTANCE_NONE(sc.distance)) break;

    Tr_abs *= transmissivity(cmd, args, RNATM_RADCOEF_Ka, pos, dir, sc.distance);

    /* Full absorption. Stop the path */
    if(Tr_abs == 0) break;

    /* Compute the scattering position */
    sc_pos[0] = pos[0] + dir[0] * sc.distance;
    sc_pos[1] = pos[1] + dir[1] * sc.distance;
    sc_pos[2] = pos[2] + dir[2] * sc.distance;

    /* Scattering in volume */
    if(!SURFACE_SCATTERING(&sc)) {
      double Ls; /* Direct contribution at the scattering position */
      Ls = volume_scattering(cmd, args, sc_pos, dir, sc_dir);
      L += Tr_abs * Ls;

    /* Surface scattering */
    } else {
      double Ls; /* Direct contribution at the scattering position */
      double refl; /* Surface reflectivity */
      Ls = surface_scattering(cmd, args, &sc, sc_pos, dir, sc_dir, &refl);
      L += Tr_abs * Ls;

      /* Russian roulette wrt surface reflectivity */
      if(ssp_rng_canonical(args->rng) >= refl) break;
    }

    d3_set(pos, sc_pos);
    d3_set(dir, sc_dir);

    ++nsc;
  }

  return L;
}
