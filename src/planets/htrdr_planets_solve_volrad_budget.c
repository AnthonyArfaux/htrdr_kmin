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

#include "core/htrdr_accum.h"
#include "core/htrdr_log.h"
#include "core/htrdr_solve_buffer.h"

#include <rsys/clock_time.h>

#include <star/ssp.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
solve_volumic_radiative_budget
  (struct htrdr* htrdr,
   const struct htrdr_solve_item_args* args,
   void* data)
{
  /* Monte Carlo accumulators */
  struct htrdr_accum budget = HTRDR_ACCUM_NULL;
  struct htrdr_accum time = HTRDR_ACCUM_NULL;

  /* Miscellaneous */
  struct htrdr_planets* cmd = NULL;
  struct planets_voxel_radiative_budget* voxel = data;
  size_t i = 0;

  /* Preconditions */
  ASSERT(htrdr && args && data);
  (void)htrdr, (void)cmd;

  cmd = args->context;
  ASSERT(cmd);
  ASSERT(cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_VOLUMIC_RADIATIVE_BUDGET);

  FOR_EACH(i, 0, args->nrealisations) {
    /* Time recording */
    struct time t0, t1;
    double usec = 0;

    /* Miscellaneous */
    double w = 0;
    double r = 0;

    /* Start of realisation time recording */
    time_current(&t0);

    /* For now, the realisation statistically evaluates 4+1
     * TODO implement the actual calculation */
    r = ssp_rng_canonical(args->rng);
    if(r < 0.5) {
      w = 2;
    } else {
      w = 8;
    }

    /* Stop time recording */
    time_sub(&t0, time_current(&t1), &t0);
    usec = (double)time_val(&t0, TIME_NSEC) * 0.001;

    /* Update the volumic radiance budget accumulator */
    budget.sum_weights += w;
    budget.sum_weights_sqr += w*w;
    budget.nweights += 1;

    /* Update the per realisation time accumulator */
    time.sum_weights += usec;
    time.sum_weights_sqr += usec*usec;
    time.nweights += 1;
  }

  /* Write voxel data */
  voxel->volrad_budget = budget;
  voxel->time = time;
}

static res_T
write_buffer
  (struct htrdr_planets* cmd,
   struct htrdr_buffer* buf,
   struct htrdr_accum* budget, /* W/m^3 */
   struct htrdr_accum* time, /* us */
   FILE* stream)
{
  struct htrdr_buffer_layout layout = HTRDR_BUFFER_LAYOUT_NULL;
  size_t x = 0;

  /* Preconditions */
  ASSERT(cmd && buf && budget && time && stream);

  htrdr_buffer_get_layout(buf, &layout);
  ASSERT(layout.height == 1);
  (void)cmd;

  /* Reset global accumulators */
  *budget = HTRDR_ACCUM_NULL;
  *time = HTRDR_ACCUM_NULL;

  FOR_EACH(x, 0, layout.width) {
    struct planets_voxel_radiative_budget* voxel = htrdr_buffer_at(buf, x, 0);
    struct htrdr_estimate estim_budget = HTRDR_ESTIMATE_NULL;
    struct htrdr_estimate estim_time = HTRDR_ESTIMATE_NULL;

    htrdr_accum_get_estimation(&voxel->volrad_budget, &estim_budget);
    htrdr_accum_get_estimation(&voxel->time, &estim_time);

    /* Write the estimate of the volumic radiative budget */
    fprintf(stream, "%g %g ", estim_budget.E, estim_budget.SE);

    /* Write the accumulator of the volumic radiative budget */
    fprintf(stream, "%g %g %lu ",
      voxel->volrad_budget.sum_weights,
      voxel->volrad_budget.sum_weights_sqr,
      (unsigned long)voxel->volrad_budget.nweights);

    /* Write the estimate of the per realisation time */
    fprintf(stream, "%g %g\n", estim_time.E, estim_time.SE);

    /* Update the overall volumic radiative budget accumulator */
    budget->sum_weights += voxel->volrad_budget.sum_weights;
    budget->sum_weights_sqr += voxel->volrad_budget.sum_weights_sqr;
    budget->nweights += voxel->volrad_budget.nweights;

    /* Update the overall per realisation time accumulator */
    time->sum_weights += voxel->time.sum_weights;
    time->sum_weights_sqr += voxel->time.sum_weights_sqr;
    time->nweights += voxel->time.nweights;
  }
   return RES_OK;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
planets_solve_volrad_budget(struct htrdr_planets* cmd)
{
  struct htrdr_solve_buffer_args args = HTRDR_SOLVE_BUFFER_ARGS_NULL;
  struct htrdr_accum acc_budget = HTRDR_ACCUM_NULL;
  struct htrdr_accum acc_time = HTRDR_ACCUM_NULL;
  struct htrdr_estimate estim_budget = HTRDR_ESTIMATE_NULL;
  struct htrdr_estimate estim_time = HTRDR_ESTIMATE_NULL;
  res_T res = RES_OK;

  /* Preconditions */
  ASSERT(cmd);
  ASSERT(cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_VOLUMIC_RADIATIVE_BUDGET);

  args.solve_item = solve_volumic_radiative_budget;
  args.buffer_layout = cmd->buf_layout;
  args.nrealisations = cmd->spt;
  args.context = cmd;

  res = htrdr_solve_buffer(cmd->htrdr, &args, cmd->buf);
  if(res != RES_OK) goto error;

  /* Nothing more to do on non master processes */
  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) goto exit;

  /* Write outut data */
  res = write_buffer(cmd, cmd->buf, &acc_budget, &acc_time, cmd->output);
  if(res != RES_OK) goto error;

  htrdr_accum_get_estimation(&acc_time, &estim_time);
  htrdr_accum_get_estimation(&acc_budget, &estim_budget);

  /* Write overall results */
  htrdr_log(cmd->htrdr, "Time per radiative path (in µs): %g +/- %g\n",
    estim_time.E, estim_time.SE);
  htrdr_log(cmd->htrdr, "Volumic radiative budget (in W/m³): %g +/- %g\n",
    estim_budget.E, estim_budget.SE);

exit:
  return res;
error:
  goto exit;
}
