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

#include "core/htrdr.h"
#include "core/htrdr_log.h"

#include <star/ssp.h>

#include <rsys/double3.h>
#include <rsys/ref_count.h>

struct htrdr_planeto_source {
  double position[3]; /* In m */

  double radius; /* In m */
  double temperature; /* In Kelvin */

  ref_T ref;
  struct htrdr* htrdr;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
release_source(ref_T* ref)
{
  struct htrdr_planeto_source* source;
  struct htrdr* htrdr;
  ASSERT(ref);

  source = CONTAINER_OF(ref, struct htrdr_planeto_source, ref);
  htrdr = source->htrdr;
  MEM_RM(htrdr_get_allocator(htrdr), source);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_planeto_source_create
  (struct htrdr_planeto* cmd,
   const struct htrdr_planeto_source_args* args,
   struct htrdr_planeto_source** out_source)
{
  struct htrdr_planeto_source* src = NULL;
  double dst; /* In m */
  double lat; /* In radians */
  double lon; /* In radians */
  res_T res = RES_OK;
  ASSERT(cmd && out_source);
  ASSERT(htrdr_planeto_source_args_check(args) == RES_OK);

  src = MEM_CALLOC(htrdr_get_allocator(cmd->htrdr), 1, sizeof(*src));
  if(!src) {
    htrdr_log_err(cmd->htrdr, "error allocating source\n");
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&src->ref);
  htrdr_ref_get(cmd->htrdr);
  src->htrdr = cmd->htrdr;
  src->radius = args->radius * 1e3/*From km to m*/;
  src->temperature = args->temperature;

  /* Convert latitude and longitude to radians and distance in m */
  lat = MDEG2RAD(args->latitude);
  lon = MDEG2RAD(args->longitude);
  dst = args->distance * 1e3/*From km to m*/;

  /* Compute the position of the source */
  src->position[0] = dst * cos(lat) * cos(lon);
  src->position[1] = dst * cos(lat) * sin(lon);
  src->position[2] = dst * sin(lat);

exit:
  *out_source = src;
  return res;
error:
  if(src) { htrdr_planeto_source_ref_put(src); src = NULL; }
  goto exit;
}

void
htrdr_planeto_source_ref_get(struct htrdr_planeto_source* source)
{
  ASSERT(source);
  ref_get(&source->ref);
}

void htrdr_planeto_source_ref_put(struct htrdr_planeto_source* source)
{
  ASSERT(source);
  ref_put(&source->ref, release_source);
}

double
htrdr_planeto_source_sample_direction
  (const struct htrdr_planeto_source* source,
   struct ssp_rng* rng,
   const double pos[3],
   double dir[3])
{
  double main_dir[3];
  double half_angle; /* In radians */
  double cos_half_angle;
  double dst; /* In m */
  double pdf;
  ASSERT(source && rng && pos && dir);

  /* compute the direction of `pos' toward the center of the source */
  d3_sub(main_dir, source->position, pos);

  /* Normalize the direction and keep the distance from `pos' to the center of
   * the source */
  dst = d3_normalize(main_dir, main_dir);
  CHK(dst > source->radius);

  /* Sample the source according to its solid angle,
   * i.e. 2*PI*(1 - cos(half_angle)) */
  half_angle = asin(source->radius/dst);
  cos_half_angle = cos(half_angle);
  ssp_ran_sphere_cap_uniform(rng, main_dir, cos_half_angle, dir, &pdf);

  return pdf;
}

double
htrdr_planeto_source_get_radiance
  (const struct htrdr_planeto_source* source,
   const double wlen)
{
  return htrdr_planck_monochromatic
    (wlen*1e-9/*From nm to m*/, source->temperature);
}

int
htrdr_planeto_source_is_targeted
  (const struct htrdr_planeto_source* source,
   const double pos[3],
   const double dir[3])
{
  double main_dir[3];
  double half_angle; /* In radians */
  double dst; /* In m */
  ASSERT(source && dir && d3_is_normalized(dir));

  /* compute the direction of `pos' toward the center of the source */
  d3_sub(main_dir, source->position, pos);

  /* Normalize the direction and keep the distance from `pos' to the center of
   * the source */
  dst = d3_normalize(main_dir, main_dir);
  CHK(dst > source->radius);

  /* Compute the the half angle of the source as seen from pos */
  half_angle = asin(source->radius/dst);

  return d3_dot(dir, main_dir) >= cos(half_angle);
}
