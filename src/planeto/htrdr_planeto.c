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

#define _POSIX_C_SOURCE 200112L /* fdopen, nextafter, rint */

#include "core/htrdr.h"
#include "core/htrdr_ran_wlen_cie_xyz.h"
#include "core/htrdr_ran_wlen_planck.h"
#include "core/htrdr_log.h"

#include "planeto/htrdr_planeto.h"
#include "planeto/htrdr_planeto_args.h"
#include "planeto/htrdr_planeto_c.h"
#include "planeto/htrdr_planeto_source.h"

#include <rad-net/rnatm.h>
#include <rad-net/rngrd.h>

#include <star/scam.h>

#include <rsys/cstr.h>
#include <rsys/double3.h>
#include <rsys/mem_allocator.h>

#include <fcntl.h> /* open */
#include <math.h> /* nextafter, rint */
#include <unistd.h> /* close */
#include <sys/stat.h>

/*******************************************************************************
 * Helper function
 ******************************************************************************/
/* Calculate the minimum length of the atmospheric spectral bands for the
 * spectral domain considered */
static double
compute_min_band_len(const struct htrdr_planeto* cmd)
{
  const double* range = NULL; /* In nm */
  double len = DBL_MAX;
  size_t ibands[2];
  size_t i;
  ASSERT(cmd);

  range = cmd->spectral_domain.wlen_range;

  /* The spectral range is degenerate to a wavelength */
  if(eq_eps(range[0], range[1], 1.e-6)) {
    return 0;
  }

  RNATM(find_bands(cmd->atmosphere, cmd->spectral_domain.wlen_range, ibands));

  /* At least one band must be overlaped by the spectral domain */
  ASSERT(ibands[0]<=ibands[1]);
  FOR_EACH(i, ibands[0], ibands[1]+1) {
    double band_range[2];
    RNATM(band_get_range(cmd->atmosphere, i, band_range));

    /* Make the upper bound inclusive */
    band_range[1] = nextafter(band_range[1], 0);

    /* Clamp the band range to the spectral domain */
    band_range[0] = MMAX(band_range[0], range[0]);
    band_range[1] = MMIN(band_range[1], range[1]);
    len = MMIN(band_range[1] - band_range[0], len);
  }
  return len;
}

/* Calculate the number of fixed size spectral intervals to use for the
 * cumulative */
static size_t
compute_nintervals_for_spectral_cdf(const struct htrdr_planeto* cmd)
{
  double range_size;
  double interval_len;
  size_t nintervals;
  ASSERT(cmd);

  range_size =
    cmd->spectral_domain.wlen_range[1]
  - cmd->spectral_domain.wlen_range[0];

  /* Initially assume ~one interval per nanometer */
  nintervals = (size_t)rint(range_size);

  /* Calculate the minimum length of the atmospheric spectral bands fixed to
   * the spectral integration domain. We ensure that an interval of the
   * spectral cdf cannot be greater than this length */
  interval_len = compute_min_band_len(cmd);
  if(interval_len < (range_size / (double)nintervals)) {
    nintervals = (size_t)ceil(range_size / interval_len);
  }

  return nintervals;
}

static res_T
setup_octree_storage
  (struct htrdr_planeto* cmd,
   const struct htrdr_planeto_args* args,
   struct rnatm_create_args* rnatm_args)
{
  struct stat file_stat;
  int fd = -1;
  int err = 0;
  res_T res = RES_OK;
  ASSERT(cmd && args && rnatm_args);

  rnatm_args->octrees_storage = NULL;
  rnatm_args->load_octrees_from_storage = 0;

  if(!args->octrees_storage) goto exit;

  fd = open(args->octrees_storage, O_CREAT|O_RDWR, S_IRUSR|S_IWUSR);
  if(fd < 0) { res = RES_IO_ERR; goto error; }

  rnatm_args->octrees_storage = fdopen(fd, "w+");
  if(!rnatm_args->octrees_storage) { res = RES_IO_ERR; goto error; }

  /* From now on, manage the opened file from its pointer and not from its
   * descriptor */
  fd = -1;

  err = stat(args->octrees_storage, &file_stat);
  if(err < 0) { res = RES_IO_ERR; goto error; }

  if(file_stat.st_size != 0) {
    /* The file is not empty and therefore must contain valid octrees */
    rnatm_args->load_octrees_from_storage = 1;
  }

exit:
  cmd->octrees_storage = rnatm_args->octrees_storage;
  return res;
error:
  htrdr_log_err(cmd->htrdr, "error opening the octree storage `%s' -- %s\n",
    args->octrees_storage, res_to_cstr(res));

  if(fd >= 0) CHK(close(fd) == 0);
  if(rnatm_args->octrees_storage) CHK(fclose(rnatm_args->octrees_storage) == 0);
  rnatm_args->octrees_storage = NULL;
  rnatm_args->load_octrees_from_storage = 1;
  goto exit;
}

static res_T
setup_atmosphere
  (struct htrdr_planeto* cmd,
   const struct htrdr_planeto_args* args)
{
  struct rnatm_create_args rnatm_args = RNATM_CREATE_ARGS_DEFAULT;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  rnatm_args.gas = args->gas;
  rnatm_args.aerosols = args->aerosols;
  rnatm_args.naerosols = args->naerosols;
  rnatm_args.name = "atmosphere";
  rnatm_args.spectral_range[0] = args->spectral_domain.wlen_range[0];
  rnatm_args.spectral_range[1] = args->spectral_domain.wlen_range[1];
  rnatm_args.optical_thickness = args->optical_thickness;
  rnatm_args.grid_definition_hint = args->octree_definition_hint;
  rnatm_args.precompute_normals = 0;
  rnatm_args.logger = htrdr_get_logger(cmd->htrdr);
  rnatm_args.allocator = htrdr_get_allocator(cmd->htrdr);
  rnatm_args.nthreads = args->nthreads;
  rnatm_args.verbose = args->verbose;

  res = setup_octree_storage(cmd, args, &rnatm_args);
  if(res != RES_OK) goto error;

  res = rnatm_create(&rnatm_args, &cmd->atmosphere);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  if(cmd->atmosphere) {
    RNATM(ref_put(cmd->atmosphere));
    cmd->atmosphere = NULL;
  }
  goto exit;
}

static res_T
setup_ground
  (struct htrdr_planeto* cmd,
   const struct htrdr_planeto_args* args)
{
  struct rngrd_create_args rngrd_args = RNGRD_CREATE_ARGS_DEFAULT;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  if(cmd->output_type == HTRDR_PLANETO_ARGS_OUTPUT_OCTREES)
    goto exit;

  rngrd_args.smsh_filename = args->ground.smsh_filename;
  rngrd_args.props_filename = args->ground.props_filename;
  rngrd_args.mtllst_filename = args->ground.mtllst_filename;
  rngrd_args.name = args->ground.name;
  rngrd_args.logger = htrdr_get_logger(cmd->htrdr);
  rngrd_args.allocator = htrdr_get_allocator(cmd->htrdr);
  rngrd_args.verbose = args->verbose;

  res = rngrd_create(&rngrd_args, &cmd->ground);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  if(cmd->ground) {
    RNGRD(ref_put(cmd->ground));
    cmd->ground = NULL;
  }
  goto exit;
}

static INLINE res_T
setup_spectral_domain
  (struct htrdr_planeto* cmd,
   const struct htrdr_planeto_args* args)
{
  size_t nintervals;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  cmd->spectral_domain = args->spectral_domain;

  /* Configure the spectral distribution */
  nintervals = compute_nintervals_for_spectral_cdf(cmd);
  switch(cmd->spectral_domain.spectral_type) {
    /* Planck distribution */
    case HTRDR_SPECTRAL_LW:
    case HTRDR_SPECTRAL_SW:
      res = htrdr_ran_wlen_planck_create(cmd->htrdr,
        cmd->spectral_domain.wlen_range, nintervals,
        cmd->spectral_domain.ref_temperature, &cmd->planck);
      break;
    /* CIE XYZ distribution */
    case HTRDR_SPECTRAL_SW_CIE_XYZ:
      res = htrdr_ran_wlen_cie_xyz_create(cmd->htrdr,
        cmd->spectral_domain.wlen_range, nintervals, &cmd->cie);
      break;

    default: FATAL("Unreachable code\n"); break;
  }
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  goto exit;
}

static res_T
setup_output
  (struct htrdr_planeto* cmd,
   const struct htrdr_planeto_args* args)
{
  const char* output_name = NULL;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  /* No output stream on non master processes */
  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) {
    cmd->output = NULL;
    output_name = "<null>";

  /* Write results on stdout */
  } else if(!args->output) {
    cmd->output = stdout;
    output_name = "<stdout>";

  /* Open the output stream */
  } else {
    res = htrdr_open_output_stream(cmd->htrdr, args->output, 0/*read*/,
      args->force_output_overwrite, &cmd->output);
    if(res != RES_OK) goto error;
    output_name = args->output;
  }

  res = str_set(&cmd->output_name, output_name);
  if(res != RES_OK) {
    htrdr_log_err(cmd->htrdr, "error storing output stream name `%s' -- %s\n",
      output_name, res_to_cstr(res));
    goto error;
  }

  cmd->output_type = args->output_type;

exit:
  return res;
error:
  str_clear(&cmd->output_name);
  if(cmd->output && cmd->output != stdout) {
    CHK(fclose(cmd->output) == 0);
    cmd->output = NULL;
  }
  goto exit;
}

static INLINE res_T
setup_source
  (struct htrdr_planeto* cmd,
   const struct htrdr_planeto_args* args)
{
  res_T res = RES_OK;
  ASSERT(cmd && args);

  if(cmd->output_type == HTRDR_PLANETO_ARGS_OUTPUT_OCTREES)
    goto exit;

  res = htrdr_planeto_source_create(cmd, &args->source, &cmd->source);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  goto exit;
}

static res_T
setup_camera
  (struct htrdr_planeto* cmd,
   const struct htrdr_planeto_args* args)
{
  struct scam_perspective_args cam_args = SCAM_PERSPECTIVE_ARGS_DEFAULT;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  if(cmd->output_type != HTRDR_PLANETO_ARGS_OUTPUT_IMAGE)
    goto exit;

  ASSERT(htrdr_args_camera_perspective_check(&args->cam_persp) == RES_OK);
  ASSERT(htrdr_args_image_check(&args->image) == RES_OK);

  d3_set(cam_args.position, args->cam_persp.position);
  d3_set(cam_args.target, args->cam_persp.target);
  d3_set(cam_args.up, args->cam_persp.up);
  cam_args.field_of_view = MDEG2RAD(args->cam_persp.fov_y);
  cam_args.lens_radius = args->cam_persp.lens_radius;
  cam_args.focal_distance = args->cam_persp.focal_dst;
  cam_args.aspect_ratio =
    (double)args->image.definition[0]
  / (double)args->image.definition[1];

  res = scam_create_perspective
    (htrdr_get_logger(cmd->htrdr),
     htrdr_get_allocator(cmd->htrdr),
     htrdr_get_verbosity_level(cmd->htrdr),
     &cam_args,
     &cmd->camera);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  goto exit;
}

static res_T
setup_buffer
  (struct htrdr_planeto* cmd,
   const struct htrdr_planeto_args* args)
{
  struct htrdr_pixel_format pixfmt = HTRDR_PIXEL_FORMAT_NULL;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  if(cmd->output_type != HTRDR_PLANETO_ARGS_OUTPUT_IMAGE)
    goto exit;

  planeto_get_pixel_format(cmd, &pixfmt);

  /* Setup buffer layout */
  cmd->buf_layout.width = args->image.definition[0];
  cmd->buf_layout.height = args->image.definition[1];
  cmd->buf_layout.pitch = args->image.definition[0] * pixfmt.size;
  cmd->buf_layout.elmt_size = pixfmt.size;
  cmd->buf_layout.alignment = pixfmt.alignment;

  /* Save the number of samples per pixel */
  cmd->spp = args->image.spp;

  /* Create the image buffer only on the master process; Image parts rendered
   * by other processes are collected there */
  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) goto exit;

  res = htrdr_buffer_create(cmd->htrdr, &cmd->buf_layout, &cmd->buf);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  if(cmd->buf) { htrdr_buffer_ref_put(cmd->buf); cmd->buf = NULL; }
  goto exit;
}

static INLINE res_T
write_vtk_octrees(const struct htrdr_planeto* cmd)
{
  size_t octrees_range[2];
  res_T res = RES_OK;
  ASSERT(cmd);

  /* Nothing to do on non master process */
  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) goto exit;

  octrees_range[0] = 0;
  octrees_range[1] = rnatm_get_spectral_items_count(cmd->atmosphere) - 1;

  res = rnatm_write_vtk_octrees(cmd->atmosphere, octrees_range, cmd->output);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  goto exit;
}

static void
planeto_release(ref_T* ref)
{
  struct htrdr_planeto* cmd = CONTAINER_OF(ref, struct htrdr_planeto, ref);
  struct htrdr* htrdr = NULL;
  ASSERT(ref);

  if(cmd->atmosphere) RNATM(ref_put(cmd->atmosphere));
  if(cmd->ground) RNGRD(ref_put(cmd->ground));
  if(cmd->source) htrdr_planeto_source_ref_put(cmd->source);
  if(cmd->cie) htrdr_ran_wlen_cie_xyz_ref_put(cmd->cie);
  if(cmd->planck) htrdr_ran_wlen_planck_ref_put(cmd->planck);
  if(cmd->octrees_storage) CHK(fclose(cmd->octrees_storage) == 0);
  if(cmd->output && cmd->output != stdout) CHK(fclose(cmd->output) == 0);
  if(cmd->buf) htrdr_buffer_ref_put(cmd->buf);
  if(cmd->camera) SCAM(ref_put(cmd->camera));
  str_release(&cmd->output_name);

  htrdr = cmd->htrdr;
  MEM_RM(htrdr_get_allocator(htrdr), cmd);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Exported functions
 ******************************************************************************/
res_T
htrdr_planeto_create
  (struct htrdr* htrdr,
   const struct htrdr_planeto_args* args,
   struct htrdr_planeto** out_cmd)
{
  struct htrdr_planeto* cmd = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && out_cmd);

  res = htrdr_planeto_args_check(args);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "Invalid htrdr_planeto arguments -- %s\n",
      res_to_cstr(res));
    goto error;
  }

  cmd = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*cmd));
  if(!cmd) {
    htrdr_log_err(htrdr, "Error allocating htrdr_planeto command\n");
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&cmd->ref);
  htrdr_ref_get(htrdr);
  cmd->htrdr = htrdr;
  str_init(htrdr_get_allocator(htrdr), &cmd->output_name);

  res = setup_output(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_source(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_camera(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_ground(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_atmosphere(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_spectral_domain(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_buffer(cmd, args);
  if(res != RES_OK) goto error;

exit:
  *out_cmd = cmd;
  return res;
error:
  if(cmd) {
    htrdr_planeto_ref_put(cmd);
    cmd = NULL;
  }
  goto exit;
}

void
htrdr_planeto_ref_get(struct htrdr_planeto* cmd)
{
  ASSERT(cmd);
  ref_get(&cmd->ref);
}

void
htrdr_planeto_ref_put(struct htrdr_planeto* cmd)
{
  ASSERT(cmd);
  ref_put(&cmd->ref, planeto_release);
}

res_T
htrdr_planeto_run(struct htrdr_planeto* cmd)
{
  res_T res = RES_OK;
  ASSERT(cmd);

  switch(cmd->output_type) {
    case HTRDR_PLANETO_ARGS_OUTPUT_IMAGE:
      htrdr_log_warn(cmd->htrdr, "image rendering is not yet implemented\n");
      break;
    case HTRDR_PLANETO_ARGS_OUTPUT_OCTREES:
      res = write_vtk_octrees(cmd);
      break;
    default: FATAL("Unreachable code\n"); break;
  }
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  goto exit;
}

/*******************************************************************************
 * Local function
 ******************************************************************************/
void
planeto_get_pixel_format
  (const struct htrdr_planeto* cmd,
   struct htrdr_pixel_format* fmt)
{
  ASSERT(cmd && fmt && cmd->output_type == HTRDR_PLANETO_ARGS_OUTPUT_IMAGE);
  (void)cmd;

  switch(cmd->spectral_domain.spectral_type) {
    case HTRDR_SPECTRAL_LW:
    case HTRDR_SPECTRAL_SW:
      fmt->size = sizeof(struct planeto_pixel_xwave);
      fmt->alignment = ALIGNOF(struct planeto_pixel_xwave);
      break;
    case HTRDR_SPECTRAL_SW_CIE_XYZ:
      fmt->size = sizeof(struct planeto_pixel_image);
      fmt->alignment = ALIGNOF(struct planeto_pixel_image);
      break;
    default: FATAL("Unreachable code\n"); break;
  }
}
