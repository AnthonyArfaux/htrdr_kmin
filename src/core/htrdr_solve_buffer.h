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

#ifndef HTRDR_SOLVE_BUFFER_H
#define HTRDR_SOLVE_BUFFER_H

struct htrdr_solve_item_args {
  struct ssp_rng* rng; /* Random Number Generator */
  size_t item_id; /* Index of the item */
  size_t nrealisations; /* #realisations to estimate the item */
  size_t ithread; /* Id of the thread solving the item */
  void* context; /* User defined data */
};
#define HTRDR_SOLVE_ITEM_ARGS_NULL__ {NULL, 0, 0, 0, NULL}
static const struct htrdr_solve_item_args HTRDR_SOLVE_ITEM_ARGS_NULL =
  HTRDR_SOLVE_ITEM_ARGS_NULL__;

typedef void
(*htrdr_solve_item_T)
  (struct htrdr* htrdr,
   const struct htrdr_solve_item_args* args,
   void* item); /* Output data */

struct htrdr_solve_buffer_args {
  htrdr_solve_item_T solve_item; /* User defined functor */
  struct htrdr_buffer_layout buffer_layout;
  size_t nrealisations; /* #realisations per item */
  void* context; /* User defined data */
};
#define HTRDR_SOLVE_BUFFER_ARGS_NULL__ {                                       \
  NULL, /* Solver item functor */                                              \
  HTRDR_BUFFER_LAYOUT_NULL__, /* Layout of the destination buffer */           \
  0, /* #realisations per item */                                              \
  NULL /* User defined data */                                                 \
}
static const struct htrdr_solve_buffer_args HTRDR_SOLVE_BUFFER_ARGS_NULL =
  HTRDR_SOLVE_BUFFER_ARGS_NULL__;

/*******************************************************************************
 * Exported symbols
 ******************************************************************************/
BEGIN_DECLS

HTRDR_API res_T
htrdr_solve_buffer
  (struct htrdr* htrdr,
   const struct htrdr_solve_buffer_args* args,
   struct htrdr_buffer* buf); /* May be NULL for non master processes */

END_DECLS

#endif /* HTRDR_SOLVE_BUFFER_H */
