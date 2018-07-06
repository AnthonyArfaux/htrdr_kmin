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

struct scattering_context {
  struct ssp_rng* rng;

  double Ts; /* Sampled optical thickness */
  double traversal_dst; /* Distance traversed along the ray */
};
static const struct scattering_context SCATTERING_CONTEXT_NULL = {
  NULL, 0, 0, 0, 0
};

struct transmissivity_context {
  struct ssp_rng* rng;
  double Ts; /* Sampled optical thickness */
  double Tmin; /* Minimal optical thickness */
  double traversal_dst; /* Distance traversed along the ray */
};
static const struct transmissivity_context TRANSMISSION_CONTEXT_NULL = {
  NULL, 0, 0;
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
  const double* vox_data = NULL;
  struct scattering_context* ctx = context;
  double dst;
  double ks_max;
  int pursue_traversal = 1;
  ASSERT(hit && ctx && !SVX_HIT_NONE(hit) && org && dir && range);
  (void)range;

  vox_data = hit->voxel.data;
  ks_max = vox_data[HTRDR_Ks_MAX];

  for(;;) {
    dst = hit->distance[1] - ctx->traversal_dst;
    T = dst * ks_max; /* Compute tau for the current leaf */

    /* A collision occurs behind `dst' */
    if(ctx->Ts > T) {
      ctx->Ts -= T;
      ctx->traversal_dst = hit.distance[1];
      pursue_traversal = 1;
      break;

    /*  A real/null collision occurs before `dst' */
    } else {
      double pos[3];
      double proba_s;
      double ks;
      const float collision_dst = ctx->Ts / ks_max;

      /* Compute the traversed distance up to the challenged collision */
      dst = ctx->traversal_dst + collision_dst;
      ctx->traversal_dst = dst;

      /* Compute the world space position where a collision may occur */
      pos[0] = org[0] + dst * dir[0];
      pos[1] = org[1] + dst * dir[1];
      pos[2] = org[2] + dst * dir[2];

      ks = fetch_ks(ctx->medium, x); /* TODO */

      /* Handle the case that ks_max is not *really* the max */
      proba_s = ks / ks_max;

      if(ssp_rng_canonical(ctx->rng) < proba) {/* Collide <=> real scattering */
        pursue_traversal = 0;
        break;
      } else { /* Null collision */
        ctx->Ts = ssp_ran_exp(rng, 1); /* Sample a new optical thickness */
      }
    }
  }
  return pursue_traversal;
}

static double
transmissivty_hit_filter
  (const struct svx_hit* hit,
   const double org[3],
   const double dir[3],
   const double range[2],
   void* context)
{
  const double* vox_data = NULL;
  struct transmissivity_context* ctx = context;
  double dst;
  double k_max;
  double k_min;
  int pursue_traversal = 1;
  ASSERT(hit && ctx && !SVX_HIT_NONE(hit) && org && dir && range);
  (void)range;

  vox_data = hit->voxel.data;
  k_max = vox_data[HTRDR_K_MAX];
  k_min = vox_data[HTRDR_K_MIN];

  dst = hit->distance[1] - ctx->traversal_dst;
  ctx->Tmin += dst * k_min;

  for(;;) {
    Tdif = dst * (k_max-k_min);

    /* A collision occurs behind `dst' */
    if(ctx->Ts > ctx->Tdif) {
      ctx->Ts -= Tdif;
      ctx->traversal_dst = hit->range[1];
      pursue_traversal = 1;
      break;

    /*  A real/null collision occurs before `dst' */
    } else {
      double x[3];
      double k;
      double dst;
      double proba_kext;
      double collision_dst = cts->Ts / (k_max - k_min);

      /* Compute the traversed distance up to the challenged collision */
      dst = ctx->traversal_dst + collision_dst;
      ctx->traversal_dst = dst;

      /* Compute the world space position where a collision may occur */
      x[0] = org[0] + dst * dir[0];
      x[1] = org[1] + dst * dir[1];
      x[2] = org[2] + dst * dir[2];

      k = fetch_k(ctx->medium, x); /* TODO */

      /* Handle the case that k_max is not *really* the max */
      proba = (k - k_min) / (k_max - k_min);

      if(ssp_rng_canonical(ctx->rng) < proba) { /* Collide */
        pursue_traversal = 0;
        break;
      } else { /* Null collision */
        ctx->Ts = ssp_ran_exp(rng, 1); /* Sample a new optical thickness */
        dst = hit->distance[1] - ctx->traversal_dst;
      }
    }
  }
  return pursue_traversal;
}

static double
transmissivity
  (struct htrdr* htrdr,
   struct ssp_rng* rng,
   const double pos[3],
   const double dir[3],
   const double range[2],
   const struct s3d_hit* hit_prev) /* May be NULL */
{
  struct s3d_hit s3d_hit;
  struct svx_hit svx_hit;

  struct transmissivity_context transmissivity_ctx = TRANSMISSION_CONTEXT_NULL;
  const struct s3d_hit s3d_hit_prev = hit_prev ? *hit_prev : S3D_HIT_NULL;
  float ray_pos[3];
  float ray_dir[3];
  float ray_range[2];

  ASSERT(htrdr && rng && pos && dir && range);

  /* Find the first intersection with a surface */
  f3_set_d3(ray_pos, pos);
  f3_set_d3(ray_dir, dir);
  f2_set_d2(ray_range, range);
  S3D(scene_view_trace(htrdr->s3d_scn, ray_pos, ray_dir, ray_range,
    &s3d_hit_prev, &s3d_hit));

  if(!S3D_HIT_NONE(&s3d_hit)) return 0;

  transmissivity_ctx->rng = rng;
  transmissivity_ctx->Ts = ssp_ran_exp(rng, 1); /* Sample an optical thickness */

  /* Compute the transmissivity */
  SVX(octree_trace_ray(htrdr->clouds, pos, dir, range, NULL,
    transmissivity_hit_filter, &transmissivity_ctx, &svx_hit));

  if(SVX_HIT_NONE(svx_hit)) {
    return transmissivity_ctx.Tmin ? exp(-ctx.Tmin) : 1.0;
  } else {
    return 0;
  }
}

/* Compute the radiance in short wave */
static double
radiance_sw
  (struct htrdr* htrdr,
   struct ssp_rng* rng,
   const double pos_in[3],
   const double dir_in[3],
   const double wavelength)
{
  struct s3d_hit s3d_hit = S3D_HIT_NULL;
  struct svx_hit svx_hit = SVX_HIT_NULL;
  struct s3d_hit s3d_hit_prev = S3D_HIT_NULL;

  double pos[3];
  double dir[3];
  double range[2];
  double pos_next[3];
  double dir_next[3];

  double R;
  double r; /* Random number */
  double wi[3]; /* -dir */
  double pdf;
  double sun_solid_angle;
  double Tr; /* Overall transmissivity */
  double Tr_abs; /* Absorption transmissivity */
  double L_sun; /* Sun radiance in W.m^-2.sr^-1 */
  double sun_dir[3];
  double ksi = 1; /* Throughput */
  double w = 0; /* MC weight */

  float ray_pos[3];
  float ray_dir[3];
  float ray_range[2];

  ASSERT(htrdr && rng && pos && dir);

  sun_solid_angle = htrdr_sun_get_solid_angle(htrdr->sun);
  d3_set(pos, pos_in);
  d3_set(dir, dir_in);

  /* Check that the sun is directly visible */
  f3_set_d3(ray_pos, pos);
  f3_set_d3(ray_dir, sun_dir);
  f2(ray_range, 0, FLT_MAX);
  S3D(scene_view_trace
    (htrdr->s3d_scn, ray_pos, ray_dir, ray_range, NULL, &s3d_hit));

  /* Add the direct contribution of the sun */
  if(!S3D_HIT_NONE(&s3d_hit)) {
    d3(range, 0, FLT_MAX);
    L_sun = htrdr_sun_get_luminance(htrdr->sun, sun_dir, wavelength);
    Tr = transmissivity(htrdr, rng, pos, sun_dir, range);
    w = L_sun * Tr;
  }

  /* Radiative random walk */
  for(;;) {
    struct scattering_context scattering_ctx = SCATTERING_CONTEXT_NULL;

    /* Setup the ray to trace */
    f3_set_d3(ray_pos, pos);
    f3_set_d3(ray_dir, dir);
    f2(ray_range, 0, FLT_MAX);

    /* Find the first intersection with a surface */
    S3D(scene_view_trace(htrdr->s3d_scn, ray_pos, ray_dir, ray_range,
      &s3d_hit_prev, &s3d_hit));

    /* Sample an optical thicknes */
    scattering_ctx->rng = rng;
    scattering_ctx->Ts = ssp_ran_exp(rng, 1);

    /* Define if a scattering event occurs */
    d2(range, 0, s3d_hit.distance);
    SVX(octree_trace_ray(htrdr->clouds, pos, dir, range, NULL,
      scattering_hit_filter, &scattering_ctx, &svx_hit));

    /* No scattering and no surface reflection. Stop the radiative random walk */
    if(S3D_HIT_NONE(&s3d_hit) && SVX_HIT_NONE(&svx_hit)) {
      w *= ksi;
      break;
    }

    /* Negative the incoming dir to match the convention of the SSF library */
    d3_minus(wi, dir);

    /* Compute the new position */
    pos_next[0] = pos[0] + dir[0]*scattering_ctx->traversal_dst;
    pos_next[1] = pos[1] + dir[1]*scattering_ctx->traversal_dst;
    pos_next[2] = pos[2] + dir[2]*scattering_ctx->traversal_dst;

    /* Define the absorption transmissivity from the current position to the
     * next position */
    d2(range, 0, scattering_ctx->traversal_dst);
    Tr_abs = transmissivity_absorption(htrdr, rng, pos, dir, range);
    if(Tr_abs <= 0) break;

    /* Define the transmissivity from the new position to the sun */
    d2(range, 0, FLT_MAX);
    htrdr_sun_sample_dir(htrdr->sun, rng, pos_next, sun_dir);
    Tr = transmissivity(htrdr, rng, pos_next, sun_dir, range);

    /* Scattering at a surface */
    if(SVX_HIT_NONE(&svx_hit)) {
      double reflectivity;
      double R;
      int type;

      reflectivity = ssf_bsdf_sample
        (htrdr->bsdf, rng, wi, s3d_hit.normal, dir_next, &type, &pdf);
      if(ssp_rng_canonical(rng) > reflectivity) break; /* Russian roulette */

      R = ssf_bsdf_eval(htrdr->bsdf, wi, s3d_hit.normal, sun_dir);

    /* Scattering in the medium */
    } else {
      struct ssf_phase* phase;
      double ks_partical; /* Scattering coefficient of the particles */
      double ks_gaz; /* Scattering coefficient of the gaz */
      double ks; /* Overall scattering coefficient */

      ks_gaz = htrdr_medium_get_ks_gaz(htrdr->medium);
      ks_cloud = hrdr_medium_get_ks_particle(htrdr->medium);
      ks = ks_particle + ks_gaz;

      r = ssp_rng_canonical(rng);
      if(r < ks_gaz / ks) { /* Gaz scattering */
        phase = htrdr->phase_rayleigh;
      } else { /* Cloud scattering */
        phase = htrdr->phase_hg;
      }

      ssf_phase_sample(phase, rng, wi, dir_next, &pdf);
      R = ssf_phase_eval(phase, wi, sun_dir);
    }

    /* Update the MC weight */
    ksi *= Tr_abs;
    w += ksi * L_sun * sun_solid_angle * Tr * R;

    /* Setup the next random walk state */
    d3_set(pos, pos_next);
    d3_set(dir, dir_next);
  }
  return w;
}


