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

#include "planets/htrdr_planets_c.h"
#include "planets/htrdr_planets_source.h"

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
#define SURFACE_EVENT(Event) (!S3D_HIT_NONE(&(Event)->hit))

struct event {
  /* Set to S3D_HIT_NULL if the event occurs in volume.*/
  struct s3d_hit hit;

  /* The surface normal is defined only if event is on the surface. It is
   * normalized and looks towards the incoming direction */
  double normal[3];

  /* Cells in which the event position is located. It makes sense only for an
   * event in volume */
  struct rnatm_cell_pos cells[RNATM_MAX_COMPONENTS_COUNT];

  double distance; /* Distance from ray origin to scattering position */
};
#define EVENT_NULL__ {                                                         \
  S3D_HIT_NULL__, {0,0,0}, {RNATM_CELL_POS_NULL__}, DBL_MAX                    \
}
static const struct event EVENT_NULL = EVENT_NULL__;

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
  struct rnatm_cell_pos* cells;
  double distance;
};
#define SAMPLE_DISTANCE_CONTEXT_NULL__ {                                       \
  NULL, NULL, 0, 0, 0, RNATM_RADCOEFS_COUNT__, 0, NULL, DBL_MAX                \
}
static const struct sample_distance_context SAMPLE_DISTANCE_CONTEXT_NULL =
  SAMPLE_DISTANCE_CONTEXT_NULL__;

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static INLINE res_T
check_planets_compute_radiance_args
  (const struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args)
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
  get_k_args.cells = ctx->cells;
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
      double pos[3];
      double  k, r;

      /* Calculate the distance travelled to the collision to be queried */
      dst_travelled += vox_dst_collision;

      /* Retrieve the radiative coefficient at the collision position */
      pos[0] = org[0] + dst_travelled * dir[0];
      pos[1] = org[1] + dst_travelled * dir[1];
      pos[2] = org[2] + dst_travelled * dir[2];
      RNATM(fetch_cell_list(ctx->atmosphere, pos, get_k_args.cells, NULL));
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
  (struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args,
   struct rnatm_cell_pos* cells,
   const enum rnatm_radcoef radcoef,
   const double pos[3],
   const double dir[3],
   const double range[2])
{
  struct rnatm_trace_ray_args rt = RNATM_TRACE_RAY_ARGS_NULL;
  struct sample_distance_context ctx = SAMPLE_DISTANCE_CONTEXT_NULL;
  struct svx_hit hit;
  ASSERT(cmd && args && cells && pos && dir && d3_is_normalized(dir) && range);
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
  ctx.cells = cells;

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
  (struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args,
   const enum rnatm_radcoef radcoef,
   const double pos[3],
   const double dir[3],
   const double range_max)
{
  struct rnatm_cell_pos cells[RNATM_MAX_COMPONENTS_COUNT];
  double range[2];
  double dst = 0;
  ASSERT(range_max >= 0);

  range[0] = 0;
  range[1] = range_max;
  dst = sample_distance(cmd, args, cells, radcoef, pos, dir, range);

  if(DISTANCE_NONE(dst)) {
    return 1.0; /* No collision in the medium */
  } else {
    return 0.0; /* A (real) collision occurs */
  }
}

static double
direct_contribution
  (struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args,
   const double pos[3],
   const double dir[3],
   const struct s3d_hit* hit_from)
{
  struct rngrd_trace_ray_args rt = RNGRD_TRACE_RAY_ARGS_DEFAULT;
  struct s3d_hit hit;
  double Tr;
  double Ld;
  double src_dst;
  ASSERT(cmd && args && pos && dir);

  /* Is the source hidden? */
  d3_set(rt.ray_org, pos);
  d3_set(rt.ray_dir, dir);
  if(hit_from) rt.hit_from = *hit_from;
  RNGRD(trace_ray(cmd->ground, &rt, &hit));
  if(!S3D_HIT_NONE(&hit)) return 0;

  /* Calculate the distance between the source and `pos' */
  src_dst = htrdr_planets_source_distance_to(cmd->source, pos);
  ASSERT(src_dst >= 0);

  Tr = transmissivity(cmd, args, RNATM_RADCOEF_Kext, pos, dir, src_dst);
  Ld = htrdr_planets_source_get_radiance(cmd->source, args->wlen);
  return Ld * Tr;
}

static void
find_event
  (struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args,
   const enum rnatm_radcoef radcoef,
   const double pos[3],
   const double dir[3],
   const struct s3d_hit* hit_from,
   struct event* evt)
{
  struct rngrd_trace_ray_args rt = RNGRD_TRACE_RAY_ARGS_DEFAULT;
  struct s3d_hit hit;
  double range[2];
  double dst;
  ASSERT(cmd && args && pos && dir && hit_from && evt);

  *evt = EVENT_NULL;

  /* Look for a surface intersection */
  d3_set(rt.ray_org, pos);
  d3_set(rt.ray_dir, dir);
  d2(rt.ray_range, 0, INF);
  rt.hit_from = *hit_from;
  RNGRD(trace_ray(cmd->ground, &rt, &hit));

  /* Look for an atmospheric collision */
  range[0] = 0;
  range[1] = hit.distance;
  dst = sample_distance(cmd, args, evt->cells, radcoef, pos, dir, range);

  /* Event occurs in volume */
  if(!DISTANCE_NONE(dst)) {
    evt->distance = dst;
    evt->hit = S3D_HIT_NULL;

  /* Event is on surface */
  } else if(!S3D_HIT_NONE(&hit)) {
    /* Normalize the normal and ensure that it points to `dir' */
    d3_normalize(evt->normal, d3_set_f3(evt->normal, hit.normal));
    if(d3_dot(evt->normal, dir) > 0) d3_minus(evt->normal, evt->normal);

    evt->distance = hit.distance;
    evt->hit = hit;

  /* No event */
  } else {
    evt->distance = INF;
    evt->hit = S3D_HIT_NULL;
  }
}

static INLINE struct ssf_bsdf*
create_bsdf
  (struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args,
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
  (struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args,
   const struct rnatm_cell_pos* cells) /* Cells in which scattering occurs */
{
  struct rnatm_sample_component_args sample_args =
    RNATM_SAMPLE_COMPONENT_ARGS_NULL;
  struct rnatm_cell_create_phase_fn_args phase_fn_args =
    RNATM_CELL_CREATE_PHASE_FN_ARGS_NULL;
  struct ssf_phase* phase = NULL;
  size_t cpnt;
  ASSERT(cmd && args && cells);

  /* Sample the atmospheric scattering component */
  sample_args.cells = cells;
  sample_args.iband = args->iband;
  sample_args.iquad = args->iquad;
  sample_args.radcoef = RNATM_RADCOEF_Ks;
  sample_args.r = ssp_rng_canonical(args->rng);
  RNATM(sample_component(cmd->atmosphere, &sample_args, &cpnt));

  /* Retrieve the component cell in which the scattering position is located */
  phase_fn_args.cell = RNATM_GET_COMPONENT_CELL(cells, cpnt);
  ASSERT(!SUVM_PRIMITIVE_NONE(&phase_fn_args.cell.prim));

  /* Retrieve the component phase function */
  phase_fn_args.wavelength = args->wlen;
  phase_fn_args.r[0] = ssp_rng_canonical(args->rng);
  phase_fn_args.r[1] = ssp_rng_canonical(args->rng);
  RNATM(cell_create_phase_fn(cmd->atmosphere, &phase_fn_args, &phase));

  return phase;
}

/* In shortwave, return the contribution of the external source at the bounce
 * position. In longwave, simply return 0 */
static double
surface_bounce
  (struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args,
   const struct event* sc,
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

  /* No external source in longwave */
  if(cmd->spectral_domain.type == HTRDR_SPECTRAL_LW)
    goto exit;

  /* Calculate direct contribution for specular reflection */
  if((mask & SSF_SPECULAR)
  && (mask & SSF_REFLECTION)
  && htrdr_planets_source_is_targeted(cmd->source, sc_pos, sc_dir)) {
    const double Ld = direct_contribution(cmd, args, sc_pos, sc_dir, &sc->hit);
    L = Ld * reflectivity;

  /* Calculate direct contribution in general case */
  } else if(!(mask & SSF_SPECULAR)) {
    double wi[3], pdf;

    /* Sample a direction toward the source */
    pdf = htrdr_planets_source_sample_direction(cmd->source, args->rng, sc_pos, wi);

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

exit:
  SSF(bsdf_ref_put(bsdf));
  *out_refl = reflectivity;
  return L;
}

/* In shortwave, return the contribution at the scattering position of the
 * external source. Returns 0 in long wave */
static INLINE double
volume_scattering
  (struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args,
   const struct event* sc,
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
  ASSERT(cmd && args && sc && sc_pos && in_dir && sc_dir);

  phase = create_phase_fn(cmd, args, sc->cells);
  d3_minus(wo, in_dir); /* Match StarSF convention */

  ssf_phase_sample(phase, args->rng, wo, sc_dir, NULL);

  /* Sample a direction toward the source */
  pdf = htrdr_planets_source_sample_direction(cmd->source, args->rng, sc_pos, wi);

  /* In short wave, manage the contribution of the external source */
  switch(cmd->spectral_domain.type) {
    case HTRDR_SPECTRAL_LW:
      /* Nothing to be done */
      break;

    case HTRDR_SPECTRAL_SW:
    case HTRDR_SPECTRAL_SW_CIE_XYZ:
      /* Calculate the direct contribution at the scattering position */
      Ld = direct_contribution(cmd, args, sc_pos, wi, NULL);
      rho = ssf_phase_eval(phase, wo, wi);
      L = Ld * rho / pdf;
      break;

    default: FATAL("Unreachable code.\n"); break;
  }

  SSF(phase_ref_put(phase));
  return L;
}

static INLINE enum rnatm_radcoef
sample_volume_event_type
  (const struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args,
   struct event* evt)
{
  struct rnatm_get_radcoef_args get_k_args = RNATM_GET_RADCOEF_ARGS_NULL;
  double ka, kext;
  double r;
  ASSERT(cmd && args && evt);

  get_k_args.cells = evt->cells;
  get_k_args.iband = args->iband;
  get_k_args.iquad = args->iquad;

  /* Retrieve the absorption coefficient */
  get_k_args.radcoef = RNATM_RADCOEF_Ka;
  RNATM(get_radcoef(cmd->atmosphere, &get_k_args, &ka));

  /* Retrieve the extinction coefficient */
  get_k_args.radcoef = RNATM_RADCOEF_Kext;
  RNATM(get_radcoef(cmd->atmosphere, &get_k_args, &kext));

  r = ssp_rng_canonical(args->rng);
  if(r < ka / kext) {
    return RNATM_RADCOEF_Ka; /* Absorption */
  } else {
    return RNATM_RADCOEF_Ks; /* Scattering */
  }
}

static INLINE double
get_temperature
  (const struct htrdr_planets* cmd,
   const struct event* evt)
{
  double T = 0;
  ASSERT(cmd && evt);

  if(!SURFACE_EVENT(evt)) {
    const struct rnatm_cell_pos* cell = NULL;

    /* Get gas temperature */
    cell = &RNATM_GET_COMPONENT_CELL(evt->cells, RNATM_GAS);
    RNATM(cell_get_gas_temperature(cmd->atmosphere, cell, &T));

  } else {
    struct rngrd_get_temperature_args temp_args = RNGRD_GET_TEMPERATURE_ARGS_NULL;

    /* Get ground temperature */
    temp_args.prim = evt->hit.prim;
    temp_args.barycentric_coords[0] = evt->hit.uv[0];
    temp_args.barycentric_coords[1] = evt->hit.uv[1];
    temp_args.barycentric_coords[2] = 1 - evt->hit.uv[0] - evt->hit.uv[1];
    RNGRD(get_temperature(cmd->ground, &temp_args, &T));

  }
  return T;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
double /* Radiance in W/m²/sr/m */
planets_compute_radiance
  (struct htrdr_planets* cmd,
   const struct planets_compute_radiance_args* args)
{
  struct s3d_hit hit_from = S3D_HIT_NULL;
  struct event ev;
  double pos[3];
  double dir[3];
  double L = 0; /* Radiance in W/m²/sr/m */
  size_t nsc = 0; /* Number of surface or volume scatterings (for debug) */
  int longwave = 0;
  ASSERT(cmd && check_planets_compute_radiance_args(cmd, args) == RES_OK);

  d3_set(pos, args->path_org);
  d3_set(dir, args->path_dir);
  longwave = cmd->spectral_domain.type == HTRDR_SPECTRAL_LW;

  if(!longwave && htrdr_planets_source_is_targeted(cmd->source, pos, dir)) {
    L = direct_contribution(cmd, args, pos, dir, NULL); /* In W/m²/sr/m */
  }

  for(;;) {
    double ev_pos[3];
    double sc_dir[3];

    find_event(cmd, args, RNATM_RADCOEF_Kext, pos, dir, &hit_from, &ev);

    /* No event on surface or in volume. Stop the path */
    if(DISTANCE_NONE(ev.distance)) break;

    /* Compute the event position */
    ev_pos[0] = pos[0] + dir[0] * ev.distance;
    ev_pos[1] = pos[1] + dir[1] * ev.distance;
    ev_pos[2] = pos[2] + dir[2] * ev.distance;

    /* The event occurs on the surface */
    if(SURFACE_EVENT(&ev)) {
      double refl; /* Surface reflectivity */
      L += surface_bounce(cmd, args, &ev, ev_pos, dir, sc_dir, &refl);

      /* Check the absorption of the surface with a Russian roulette against
       * the reflectivity of the surface */
      if(ssp_rng_canonical(args->rng) >= refl) break;

      /* Save current intersected surface to avoid self-intersection when
       * sampling next event */
      hit_from = ev.hit;

    /* The event occurs in the volume */
    } else {
      enum rnatm_radcoef ev_type = sample_volume_event_type(cmd, args, &ev);
      ASSERT(ev_type == RNATM_RADCOEF_Ka || ev_type == RNATM_RADCOEF_Ks);

      /* Absorption. Stop the path */
      if(ev_type == RNATM_RADCOEF_Ka) break;

      L += volume_scattering(cmd, args, &ev, ev_pos, dir, sc_dir);
      hit_from = S3D_HIT_NULL; /* Reset surface intersection */
    }

    d3_set(pos, ev_pos);
    d3_set(dir, sc_dir);

    ++nsc;
  }

  /* From there, a valid event means that the path has stopped in surface or
   * volume absorption. In longwave, add the contribution of the internal
   * source */
  if(longwave && !DISTANCE_NONE(ev.distance)) {
    const double wlen_m = args->wlen * 1.e-9; /* wavelength in meters */
    const double temperature = get_temperature(cmd, &ev); /* In K */
    L += htrdr_planck_monochromatic(wlen_m, temperature);
  }

  return L; /* In W/m²/sr/m */
}
