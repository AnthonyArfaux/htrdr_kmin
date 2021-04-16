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

#include "combustion/htrdr_combustion_geometry_ray_filter.h"


/*******************************************************************************
 * Local functions
 ******************************************************************************/
int
geometry_ray_filter_discard_first_hit
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   void* ray_data,
   void* filter_data)
{
  struct geometry_ray_filter_context* ctx = ray_data;
  int discard = 0;
  ASSERT(hit && ray_org && ray_dir && ray_data && !S3D_HIT_NONE(hit));
  (void)ray_org, (void)ray_dir, (void)filter_data;

  /* Note that the surfaces are not intersected in order. In this filter
   * function, we save in ctx->hit_tmp__ the current "first intersection" and
   * in ctx->hit the second one that is the candidate intersection. */

  if(hit->distance > ctx->hit_tmp__.distance) {
    /* An intersection was already registered as the first hit. Accept the
     * input hit and update the ctx->hit if necessary */
    discard = 0;
    if(hit->distance < ctx->hit.distance) {
      ctx->hit = *hit;
    }
  } else {
    /* No intersection has been stored or the submitted intersection is closer
     * than the one recorded as "first intersection". This new intersection
     * becomes the new first hit and is therefore rejected while the previous
     * become the candidate intersection. */
    discard = 1;
    ASSERT(ctx->hit.distance >= ctx->hit_tmp__.distance);
    ctx->hit = ctx->hit_tmp__;
    ctx->hit_tmp__ = *hit;
  }
  return discard;
}

int
geometry_ray_filter_discard_last_hit
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   void* ray_data,
   void* filter_data)
{
  struct geometry_ray_filter_context* ctx = ray_data;
  int discard = 0;
  ASSERT(hit && ray_org && ray_dir && ray_data && !S3D_HIT_NONE(hit));
  (void)ray_org, (void)ray_dir, (void)filter_data;

  /* Note that the surfaces are not intersected in order. In this filter
   * function, we save in ctx->hit_tmp__ the current "last intersection" and
   * in ctx->hit the closest one exepted this last intersection. */

  if(S3D_HIT_NONE(&ctx->hit_tmp__)) {
    /* No intersection was registered. The submitted hit would be the last one.
     * So discard it and register the hit as the last intersect */
    discard = 1;
    ctx->hit_tmp__ = *hit;

  } else if(hit->distance > ctx->hit_tmp__.distance) {
    /* A "last intersection" has already been stored but is closer than the
     * submitted intersection. This new intersection could be the real “last
     * intersection”. So we save this intersection and reject it */
    discard = 1;
    ASSERT(S3D_HIT_NONE(&ctx->hit)
        || ctx->hit_tmp__.distance >= ctx->hit.distance);

    /* Also, if the nearest intersection was not already defined, we set it to
     * the last previously saved intersection; if it was already defined, it
     * would necessarily be the closest. */
    if(S3D_HIT_NONE(&ctx->hit)) {
      ctx->hit = ctx->hit_tmp__;
    }

    ctx->hit_tmp__ = *hit; /* Save last intersection */

  } else if(hit->distance < ctx->hit.distance) {
    /* A closest intersection that is not the last one is found. Keep and
     * register this intersection */
    discard = 0;
    ctx->hit = *hit;
  }

  return discard;
}
