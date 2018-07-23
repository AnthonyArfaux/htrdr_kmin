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

#ifndef HTRDR_H
#define HTRDR_H

#include <rsys/logger.h>
#include <rsys/ref_count.h>

/* Forward declarations */
struct htrdr_args;
struct htrdr_buffer;
struct htrdr_sky;
struct htrdr_rectangle;
struct mem_allocator;
struct s3d_device;
struct s3d_scene;
struct ssf_bsdf;
struct ssf_phase;
struct svx_device;

struct htrdr {
  struct svx_device* svx;
  struct s3d_device* s3d;

  struct s3d_scene_view* s3d_scn_view; 

  struct htrdr_sky* sky;
  struct htrdr_sun* sun;

  /* Scattering functions */
  struct ssf_bsdf* bsdf; /* BSDF of the 3D geometry */
  struct ssf_phase* phase_hg; /* Henyey & Greenstein phase function */
  struct ssf_phase* phase_rayleigh; /* Rayleigh phase function */

  struct htrdr_buffer* buf;
  struct htrdr_rectangle* rect;
  double main_dir[3]; /* Main direction */
  size_t spp; /* #samples per pixel */

  FILE* output;

  unsigned nthreads;
  int dump_vtk;
  int verbose;

  struct logger logger;
  struct mem_allocator* allocator;
};

extern LOCAL_SYM res_T
htrdr_init
  (struct mem_allocator* allocator,
   const struct htrdr_args* args,
   struct htrdr* htrdr);

extern LOCAL_SYM void
htrdr_release
  (struct htrdr* htrdr);

extern LOCAL_SYM res_T
htrdr_run
  (struct htrdr* htrdr);

extern LOCAL_SYM void
htrdr_log
  (struct htrdr* htrdr,
   const char* msg,
   ...)
#ifdef COMPILER_GCC
  __attribute((format(printf, 2, 3)))
#endif
  ;

extern LOCAL_SYM void
htrdr_log_err
  (struct htrdr* htrdr,
   const char* msg,
   ...)
#ifdef COMPILER_GCC
  __attribute((format(printf, 2, 3)))
#endif
  ;

extern LOCAL_SYM void
htrdr_log_warn
  (struct htrdr* htrdr,
   const char* msg,
   ...)
#ifdef COMPILER_GCC
  __attribute((format(printf, 2, 3)))
#endif
  ;

#endif /* HTRDR_H */

