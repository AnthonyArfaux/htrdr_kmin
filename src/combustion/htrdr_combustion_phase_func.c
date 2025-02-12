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

#include "htrdr_combustion_c.h"

#include <astoria/atrstm.h>

#include <star/ssf.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static struct ssf_phase*
combustion_fetch_phase_isotropic
  (struct htrdr_combustion* cmd,
   const double wavelength, /* In nanometer */
   const struct suvm_primitive* prim,
   const double bcoords[4],
   const size_t ithread)
{
  struct ssf_phase* phase = NULL;
  ASSERT(cmd && wavelength > 0 && prim && bcoords);
  ASSERT(ithread < htrdr_get_threads_count(cmd->htrdr));
  ASSERT(cmd->phase_func_type == HTRDR_COMBUSTION_ARGS_PHASE_FUNC_ISOTROPIC);
  (void)wavelength, (void)prim, (void)bcoords;

  /* Setup the isotropic phase function */
  phase = cmd->phase_functions[ithread];
  SSF(phase_hg_setup(phase, 0));
  return phase;
}

static struct ssf_phase*
combustion_fetch_phase_rdgfa
  (struct htrdr_combustion* cmd,
   const double wavelength, /* In nanometer */
   const struct suvm_primitive* prim,
   const double bcoords[4],
   const size_t ithread)
{
  struct atrstm_rdgfa rdgfa_param = ATRSTM_RDGFA_NULL;
  struct atrstm_fetch_rdgfa_args fetch_rdgfa_args =
    ATRSTM_FETCH_RDGFA_ARGS_DEFAULT;
  struct ssf_phase_rdgfa_setup_args setup_rdgfa_args =
    SSF_PHASE_RDGFA_SETUP_ARGS_DEFAULT;
  struct ssf_phase* phase = NULL;
  ASSERT(cmd && wavelength > 0 && prim && bcoords);
  ASSERT(ithread < htrdr_get_threads_count(cmd->htrdr));
  ASSERT(cmd->phase_func_type == HTRDR_COMBUSTION_ARGS_PHASE_FUNC_RDGFA);

  /* Retrieve the RDG-FA phase function parameters from the semi transparent
   * medium */
  fetch_rdgfa_args.wavelength = wavelength;
  fetch_rdgfa_args.prim = *prim;
  fetch_rdgfa_args.bcoords[0] = bcoords[0];
  fetch_rdgfa_args.bcoords[1] = bcoords[1];
  fetch_rdgfa_args.bcoords[2] = bcoords[2];
  fetch_rdgfa_args.bcoords[3] = bcoords[3];
  ATRSTM(fetch_rdgfa(cmd->medium, &fetch_rdgfa_args, &rdgfa_param));

  /* Setup the RDG-FA phase function */
  phase = cmd->phase_functions[ithread];
  setup_rdgfa_args.wavelength = rdgfa_param.wavelength;
  setup_rdgfa_args.fractal_dimension = rdgfa_param.fractal_dimension;
  setup_rdgfa_args.gyration_radius = rdgfa_param.gyration_radius;
  setup_rdgfa_args.simd = cmd->rdgfa_simd;
  SSF(phase_rdgfa_setup(phase, &setup_rdgfa_args));

  return phase;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
struct ssf_phase*
combustion_fetch_phase_function
  (struct htrdr_combustion* cmd,
   const double wlen, /* In nanometer */
   const struct suvm_primitive* prim,
   const double bcoords[4],
   const size_t ithread)
{
  struct ssf_phase* phase = NULL;
  switch(cmd->phase_func_type) {
    case HTRDR_COMBUSTION_ARGS_PHASE_FUNC_ISOTROPIC:
      phase = combustion_fetch_phase_isotropic(cmd, wlen, prim, bcoords, ithread);
      break;
    case HTRDR_COMBUSTION_ARGS_PHASE_FUNC_RDGFA:
      phase = combustion_fetch_phase_rdgfa(cmd, wlen, prim, bcoords, ithread);
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
  return phase;
}

