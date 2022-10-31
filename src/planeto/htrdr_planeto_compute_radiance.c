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

#include <rsys/double2.h>
#include <rsys/double3.h>

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

/*******************************************************************************
 * Local functions
 ******************************************************************************/
double
planeto_compute_radiance
  (struct htrdr_planeto* cmd,
   const struct planeto_compute_radiance_args* args)
{
  struct rngrd_trace_ray_args rt_args = RNGRD_TRACE_RAY_ARGS_DEFAULT;
  struct rngrd_create_bsdf_args bsdf_args = RNGRD_CREATE_BSDF_ARGS_NULL;
  struct s3d_hit s3d_hit0 = S3D_HIT_NULL;
  struct s3d_hit s3d_hit1 = S3D_HIT_NULL;
  struct ssf_bsdf* bsdf = NULL;
  double pos[3];
  double N[3];
  double wo[3];
  double wi[3];
  double L = 0;
  double Li = 0;
  double cos_wi_N = 0;
  double rho = 0;
  ASSERT(cmd && check_planeto_compute_radiance_args(cmd, args) == RES_OK);

  /* Look for a surface intersection */
  d3_set(rt_args.ray_org, args->path_org);
  d3_set(rt_args.ray_dir, args->path_dir);
  d2(rt_args.ray_range, 0, INF);
  RNGRD(trace_ray(cmd->ground, &rt_args, &s3d_hit0));

  /* No surface intersection */
  if(S3D_HIT_NONE(&s3d_hit0)) goto exit; /* Nothing more to do */

  /* Calculate the position of intersection */
  pos[0] = s3d_hit0.distance * rt_args.ray_dir[0] + rt_args.ray_org[0];
  pos[1] = s3d_hit0.distance * rt_args.ray_dir[1] + rt_args.ray_org[1];
  pos[2] = s3d_hit0.distance * rt_args.ray_dir[2] + rt_args.ray_org[2];

  /* Sample a direction toward the source */
  htrdr_planeto_source_sample_direction(cmd->source, args->rng, pos, wi);
  Li = htrdr_planeto_source_get_radiance(cmd->source, args->wlen);

  /* The source is not facing the current position */
  d3_normalize(N, d3_set_f3(N, s3d_hit0.normal));
  if(d3_dot(N, rt_args.ray_dir) > 0) d3_minus(N, N);
  cos_wi_N = d3_dot(wi, N);
  if(cos_wi_N < 0) goto exit; /* Nothing more to do */

  /* Check the source visibility */
  d3_set(rt_args.ray_org, pos);
  d3_set(rt_args.ray_dir, wi);
  d2(rt_args.ray_range, 0, INF);
  rt_args.hit_from = s3d_hit0;
  RNGRD(trace_ray(cmd->ground, &rt_args, &s3d_hit1));
  if(!S3D_HIT_NONE(&s3d_hit1)) goto exit; /* Nothing more to do */

  /* Fetch the surface BSDF */
  bsdf_args.prim = s3d_hit0.prim;
  bsdf_args.barycentric_coords[0] = s3d_hit0.uv[0];
  bsdf_args.barycentric_coords[1] = s3d_hit0.uv[1];
  bsdf_args.barycentric_coords[2] = 1 - s3d_hit0.uv[0] - s3d_hit0.uv[1];
  bsdf_args.wavelength = args->wlen;
  bsdf_args.r = ssp_rng_canonical(args->rng);
  RNGRD(create_bsdf(cmd->ground, &bsdf_args, &bsdf));

  d3_minus(wo, rt_args.ray_org);
  rho = ssf_bsdf_eval(bsdf, wo, N, wi);
  SSF(bsdf_ref_put(bsdf));

  L = rho * Li * cos_wi_N;

exit:
  return L;
}
