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

#include "core/htrdr_accum.h"
#include "core/htrdr_log.h"
#include "core/htrdr_ran_wlen_discrete.h"
#include "core/htrdr_ran_wlen_planck.h"
#include "core/htrdr_solve_buffer.h"

#include <star/smsh.h>
#include <star/ssp.h>

#include <rsys/clock_time.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
spectral_sampling
  (struct htrdr_planets* cmd,
   const struct htrdr_solve_item_args* args,
   double* out_wlen, /* Sampled wavelength [nm] */
   double* out_wlen_pdf, /* [nm^-1] */
   size_t* out_iband, /* Spectral band in which the sampled wavelength falls */
   size_t* out_iquad) /* Sampled quadrature point in the spectral band */
{
  size_t iband_range[2] = {0, 0};
  size_t iband = 0;
  size_t iquad = 0;
  double r0, r1, r2; /* Random Numbers */
  double wlen[2] = {0,0}; /* [nm] */
  double pdf = 0; /* [nm^1] */

  /* Preconditions */
  ASSERT(cmd && args && out_wlen && out_wlen_pdf && out_iband && out_iquad);

  r0 = ssp_rng_canonical(args->rng);
  r1 = ssp_rng_canonical(args->rng);
  r2 = ssp_rng_canonical(args->rng);

  /* Sample a wavelength with respect to the type of spectral integration */
  switch(cmd->spectral_domain.type) {
    /* Longwave */
    case HTRDR_SPECTRAL_LW:
      wlen[0] = htrdr_ran_wlen_planck_sample(cmd->planck, r0, r1, &pdf);
      break;
    /* Shortwave */
    case HTRDR_SPECTRAL_SW:
      if(htrdr_planets_source_does_radiance_vary_spectrally(cmd->source)) {
        wlen[0] = htrdr_ran_wlen_discrete_sample(cmd->discrete, r0, r1, &pdf);
      } else {
        wlen[0] = htrdr_ran_wlen_planck_sample(cmd->planck, r0, r1, &pdf);
      }
      break;
    default: FATAL("Unreachable code\n"); break;
  }
  wlen[1] = wlen[0];

  /* Find the band the wavelength belongs to */
  RNATM(find_bands(cmd->atmosphere, wlen, iband_range));
  ASSERT(iband_range[0] == iband_range[1]);
  iband = iband_range[0];

  /* Sample a quadrature point */
  RNATM(band_sample_quad_pt(cmd->atmosphere, r2, iband, &iquad));

  *out_wlen = wlen[0];
  *out_wlen_pdf = pdf;
  *out_iband = iband;
  *out_iquad = iquad;
}

static INLINE void
position_sampling
  (const struct htrdr_solve_item_args* args,
   const struct smsh_desc* desc,
   double pos[3])
{
  const double* v0 = NULL;
  const double* v1 = NULL;
  const double* v2 = NULL;
  const double* v3 = NULL;

  /* Preconditions */
  ASSERT(args && desc && pos);

  /* Retrieve the vertices of the tetrahedron */
  v0 = desc->nodes + desc->cells[args->item_id*4/*#vertices*/+0]*3/*#coords*/;
  v1 = desc->nodes + desc->cells[args->item_id*4/*#vertices*/+1]*3/*#coords*/;
  v2 = desc->nodes + desc->cells[args->item_id*4/*#vertices*/+2]*3/*#coords*/;
  v3 = desc->nodes + desc->cells[args->item_id*4/*#vertices*/+3]*3/*#coords*/;

  ssp_ran_tetrahedron_uniform(args->rng, v0, v1, v2, v3, pos, NULL/*pdf*/);
}

static double /* [W/m^2/sr/m] */
get_source
  (struct htrdr_planets* cmd,
   const double pos[3],
   const double wlen) /* [nm] */
{
  struct rnatm_cell_pos cell_pos = RNATM_CELL_POS_NULL;
  double temperature = 0; /* [K] */
  double source = 0; /* [W/m^2/sr/m] */
  const double wlen_m = wlen * 1.e-9; /* Wavelength [m] */
  ASSERT(cmd && pos);

  switch(cmd->spectral_domain.type) {
    case HTRDR_SPECTRAL_SW:
      /* In shortwave, the source is external to the system */
      source = 0; /* [W/m^2/sr/m] */
      break;

    case HTRDR_SPECTRAL_LW:
      RNATM(fetch_cell(cmd->atmosphere, pos, RNATM_GAS, &cell_pos));

      if(SUVM_PRIMITIVE_NONE(&cell_pos.prim)) {
        /* The position is not in the gas */
        source = 0; /* [W/m^2/sr/m] */

      } else {
        /* Fetch the source temperature */
        RNATM(cell_get_gas_temperature(cmd->atmosphere, &cell_pos, &temperature));
        source = htrdr_planck_monochromatic(wlen_m, temperature); /* [W/m^2/sr/m] */
      }
      break;

    default: FATAL("Unreachable code\n"); break;
  }

  return source; /* [W/m^2/sr/m] */
}

/* Return the total absorption coefficient,
 * i.e. the sum of the gas and aerosol ka */
static double
get_ka
  (struct htrdr_planets* cmd,
   const double pos[3],
   const size_t iband, /* Spectral band */
   const size_t iquad) /* Quadrature point */
{
  struct rnatm_cell_pos cells[RNATM_MAX_COMPONENTS_COUNT];
  struct rnatm_get_radcoef_args get_k_args = RNATM_GET_RADCOEF_ARGS_NULL;
  double ka = 0;

  ASSERT(cmd && pos);

  get_k_args.cells = cells;
  get_k_args.iband = iband;
  get_k_args.iquad = iquad;
  get_k_args.radcoef = RNATM_RADCOEF_Ka;

  /* Retrieve the list of aerosol and gas cells in which pos lies */
  RNATM(fetch_cell_list(cmd->atmosphere, pos, get_k_args.cells, NULL));

  /* Retrive the total absorption coefficient */
  RNATM(get_radcoef(cmd->atmosphere, &get_k_args, &ka));

  return ka;
}

static void
realisation
  (struct htrdr_planets* cmd,
   const struct htrdr_solve_item_args* args,
   const struct smsh_desc* volrad_mesh_desc,
   double weights[PLANETS_VOLRAD_WEIGHTS_COUNT]) /* [W/m^3] */
{
  struct planets_compute_radiance_args rad_args =
    PLANETS_COMPUTE_RADIANCE_ARGS_NULL;

  /* Spectral integration */
  double wlen = 0; /* Wavelength [nm] */
  double wlen_pdf_nm = 0; /* Wavelength pdf [nm^-1] */
  double wlen_pdf_m = 0; /* Wavelength pdf [m^-1] */
  size_t iband = 0; /* Spectral band */
  size_t iquad = 0; /* Quadrature point */

  /* Spatial & angular integration */
  double dir_src[3] = {0,0,0}; /* Direction toward the source */
  double dir[3] = {0,0,0};
  double pos[3] = {0,0,0};
  double dir_src_pdf = 0;
  double dir_pdf = 0;

  double S = 0; /* Source [W/m^2/sr/m] */
  double L_direct  = 0; /* Direct radiance [W/m^2/sr/m] */
  double L_diffuse = 0; /* Diffuse radiance [W/m^2/sr/m] */
  double ka = 0; /* Absorption coefficient */

  /* Preconditions */
  ASSERT(cmd && args && volrad_mesh_desc);
  ASSERT(cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_VOLUMIC_RADIATIVE_BUDGET);

  /* Initialise the weights */
  memset(weights, 0, sizeof(double)*PLANETS_VOLRAD_WEIGHTS_COUNT);

  spectral_sampling(cmd, args, &wlen, &wlen_pdf_nm, &iband, &iquad);
  position_sampling(args, volrad_mesh_desc, pos);
  ssp_ran_sphere_uniform(args->rng, dir, &dir_pdf);

  S = get_source(cmd, pos, wlen); /* [W/m^2/sr/m] */

  ka = get_ka(cmd, pos, iband, iquad);
  wlen_pdf_m = wlen_pdf_nm * 1.e9; /* Transform pdf from nm^-1 to m^-1 */

  /* Compute the radiance in W/m^2/sr/m */
  d3_set(rad_args.path_org, pos);
  rad_args.rng = args->rng;
  rad_args.ithread = args->ithread;
  rad_args.wlen = wlen; /* [nm] */
  rad_args.iband = iband;
  rad_args.iquad = iquad;

  if(cmd->spectral_domain.type == HTRDR_SPECTRAL_LW) {
    /* In the longwave (radiation due to the medium), simply sample a radiative
     * path for the sampled direction and position: the radiance is considered
     * as purely diffuse. */
    d3_set(rad_args.path_dir, dir);
    L_diffuse = planets_compute_radiance(cmd, &rad_args); /* [W/m^2/sr/m] */

    /* Calculate the weights [W/m^3] */
    weights[PLANETS_VOLRAD_DIRECT]  = 0.0;
    weights[PLANETS_VOLRAD_DIFFUSE] = ka * (L_diffuse - S) / (wlen_pdf_m * dir_pdf);

  } else {
    /* In the so-called shortwave region (actually, the radiation due the
     * external source) is decomposed in its direct and diffuse components */

    dir_src_pdf = htrdr_planets_source_sample_direction
      (cmd->source, args->rng, pos, dir_src);

    d3_set(rad_args.path_dir, dir_src);
    rad_args.component = PLANETS_RADIANCE_CPNT_DIRECT;
    L_direct = planets_compute_radiance(cmd, &rad_args); /* [W/m^2/sr/m] */

    d3_set(rad_args.path_dir, dir);
    rad_args.component = PLANETS_RADIANCE_CPNT_DIFFUSE;
    L_diffuse = planets_compute_radiance(cmd, &rad_args); /* [W/m^2/sr/m] */

    /* Calculate the weights [W/m^3] */
    weights[PLANETS_VOLRAD_DIRECT]  = ka * (L_direct  - S) / (wlen_pdf_m * dir_src_pdf);
    weights[PLANETS_VOLRAD_DIFFUSE] = ka * (L_diffuse - S) / (wlen_pdf_m * dir_pdf);
  }

  /* Calculate the weights [W/m^3] */
  weights[PLANETS_VOLRAD_TOTAL] =
    weights[PLANETS_VOLRAD_DIRECT]
  + weights[PLANETS_VOLRAD_DIFFUSE];
}

static void
solve_volumic_radiative_budget
  (struct htrdr* htrdr,
   const struct htrdr_solve_item_args* args,
   void* data)
{
  /* Volumic mesh on which volumic radiative budget is estimated */
  struct smsh_desc volrad_mesh_desc = SMSH_DESC_NULL;

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

  SMSH(get_desc(cmd->volrad_mesh, &volrad_mesh_desc));

  /* Initialse voxel accumulators to 0 */
  *voxel = PLANETS_VOXEL_RADIATIVE_BUDGET_NULL;

  FOR_EACH(i, 0, args->nrealisations) {
    /* Time recording */
    struct time t0, t1;

    /* Monte Carlo weights */
    double w[PLANETS_VOLRAD_WEIGHTS_COUNT] = {0}; /* [W/m^3] */
    double usec = 0; /* [us] */

    int iweight = 0;

    /* Start of realisation time recording */
    time_current(&t0);

    /* Run the realisation */
    realisation(cmd, args, &volrad_mesh_desc, w);

    /* Stop time recording */
    time_sub(&t0, time_current(&t1), &t0);
    usec = (double)time_val(&t0, TIME_NSEC) * 0.001;

    FOR_EACH(iweight, 0, PLANETS_VOLRAD_WEIGHTS_COUNT){
      /* Update the volumic radiance budget accumulator */
      voxel->volrad_budget[iweight].sum_weights += w[iweight];
      voxel->volrad_budget[iweight].sum_weights_sqr += w[iweight]*w[iweight];
      voxel->volrad_budget[iweight].nweights += 1;
    }

    /* Update the per realisation time accumulator */
    voxel->time.sum_weights += usec;
    voxel->time.sum_weights_sqr += usec*usec;
    voxel->time.nweights += 1;
  }
}

static res_T
write_buffer
  (struct htrdr_planets* cmd,
   struct htrdr_buffer* buf,
   struct htrdr_accum* out_budget, /* W/m^3 */
   struct htrdr_accum* out_time, /* us */
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
  *out_budget = HTRDR_ACCUM_NULL;
  *out_time = HTRDR_ACCUM_NULL;

  FOR_EACH(x, 0, layout.width) {
    struct planets_voxel_radiative_budget* voxel = htrdr_buffer_at(buf, x, 0);
    struct htrdr_estimate estim_time = HTRDR_ESTIMATE_NULL;
    struct htrdr_accum* budget = NULL;
    int iestim = 0;

    budget = voxel->volrad_budget;
    FOR_EACH(iestim, 0, PLANETS_VOLRAD_WEIGHTS_COUNT) {
      struct htrdr_estimate estim_budget = HTRDR_ESTIMATE_NULL;

      /* Write the estimate of the volumic radiative budget */
      htrdr_accum_get_estimation(&budget[iestim], &estim_budget);
      fprintf(stream, "%g %g ", estim_budget.E, estim_budget.SE);

      /* Write the accumulator of the volumic radiative budget */
      fprintf(stream, "%g %g ",
        budget[iestim].sum_weights,
        budget[iestim].sum_weights_sqr);
    }

    /* Write the number of realisations.
     * It must be the same for all components */
    ASSERT(budget[PLANETS_VOLRAD_TOTAL].nweights
        == budget[PLANETS_VOLRAD_DIRECT].nweights);
    ASSERT(budget[PLANETS_VOLRAD_TOTAL].nweights
        == budget[PLANETS_VOLRAD_DIFFUSE].nweights);
    fprintf(stream, "%lu ", (unsigned long)budget[PLANETS_VOLRAD_TOTAL].nweights);

    /* Write the estimate of the per realisation time */
    htrdr_accum_get_estimation(&voxel->time, &estim_time);
    fprintf(stream, "%g %g\n", estim_time.E, estim_time.SE);

    /* Update the overall (total) volumic radiative budget accumulator */
    out_budget->sum_weights += budget[PLANETS_VOLRAD_TOTAL].sum_weights;
    out_budget->sum_weights_sqr += budget[PLANETS_VOLRAD_TOTAL].sum_weights_sqr;
    out_budget->nweights += budget[PLANETS_VOLRAD_TOTAL].nweights;

    /* Update the overall per realisation time accumulator */
    out_time->sum_weights += voxel->time.sum_weights;
    out_time->sum_weights_sqr += voxel->time.sum_weights_sqr;
    out_time->nweights += voxel->time.nweights;
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
