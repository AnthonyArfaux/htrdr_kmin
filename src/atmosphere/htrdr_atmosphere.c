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

#define _POSIX_C_SOURCE 200112L

#include "atmosphere/htrdr_atmosphere.h"
#include "atmosphere/htrdr_atmosphere_c.h"
#include "atmosphere/htrdr_atmosphere_args.h"
#include "atmosphere/htrdr_atmosphere_ground.h"
#include "atmosphere/htrdr_atmosphere_sun.h"

#include "core/htrdr_buffer.h"
#include "core/htrdr_cie_xyz.h"
#include "core/htrdr_log.h"
#include "core/htrdr_materials.h"
#include "core/htrdr_ran_wlen.h"
#include "core/htrdr_rectangle.h"

#include <high_tune/htsky.h>

#include <star/scam.h>

#include <rsys/cstr.h>
#include <rsys/double3.h>

#include <math.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static double
compute_sky_min_band_len(struct htsky* sky, const double range[2])
{
  double min_band_len = DBL_MAX;
  size_t nbands;
  ASSERT(sky && range && range[0] <= range[1]);

  nbands = htsky_get_spectral_bands_count(sky);

  if(eq_eps(range[0], range[1], 1.e-6)) {
    ASSERT(nbands == 1);
    min_band_len = 0;
  } else {
    size_t i = 0;

    /* Compute the length of the current band clamped to the submitted range */
    FOR_EACH(i, 0, nbands) {
      const size_t iband = htsky_get_spectral_band_id(sky, i);
      double wlens[2];
      HTSKY(get_spectral_band_bounds(sky, iband, wlens));

      /* Adjust band boundaries to the submitted range */
      wlens[0] = MMAX(wlens[0], range[0]);
      wlens[1] = MMIN(wlens[1], range[1]);

      min_band_len = MMIN(wlens[1] - wlens[0], min_band_len);
    }
  }
  return min_band_len;
}

/* Compute the number of fixed size bands used to discretized the spectral
 * range */
static size_t
compute_spectral_bands_count(const struct htrdr_atmosphere* cmd)
{
  double wlen_range[2];
  double wlen_range_size;
  size_t nbands;
  double band_len;
  double band_len_max;
  ASSERT(cmd);

  /* Compute size of the spectral range in nanometers */
  wlen_range[0] = cmd->wlen_range_m[0]*1.e9;
  wlen_range[1] = cmd->wlen_range_m[1]*1.e9;
  wlen_range_size = wlen_range[1] - wlen_range[0];

  /* Define as many intervals as wavelengths count in the spectral range */
  nbands = (size_t)rint(wlen_range_size);

  /* Compute the size in nanometers of an interval */
  band_len = wlen_range_size / (double)nbands;

  /* Compute the minimum band length of the sky spectral data and define it
   * as the maximum length that the bands can have */
  band_len_max = compute_sky_min_band_len(cmd->sky, wlen_range);

  /* Adjust the bands count to ensure that each sky spectral interval is
   * overlapped by at least one band */
  if(band_len > band_len_max) {
    nbands = (size_t)ceil(wlen_range_size / band_len_max);
  }
  return nbands;
}

static enum htsky_spectral_type
htrdr_to_sky_spectral_type(const enum htrdr_spectral_type type)
{
  enum htsky_spectral_type spectype;
  switch(type) {
    case HTRDR_SPECTRAL_LW:
      spectype = HTSKY_SPECTRAL_LW;
      break;
    case HTRDR_SPECTRAL_SW:
    case HTRDR_SPECTRAL_SW_CIE_XYZ:
      spectype = HTSKY_SPECTRAL_SW;
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
  return spectype;
}

static INLINE void
spherical_to_cartesian_dir
  (const double azimuth, /* In radians */
   const double elevation, /* In radians */
   double dir[3])
{
  double cos_azimuth;
  double sin_azimuth;
  double cos_elevation;
  double sin_elevation;
  ASSERT(azimuth >= 0 && azimuth < 2*PI);
  ASSERT(elevation >= 0 && elevation <= PI/2.0);
  ASSERT(dir);

  cos_azimuth = cos(azimuth);
  sin_azimuth = sin(azimuth);
  cos_elevation = cos(elevation);
  sin_elevation = sin(elevation);

  dir[0] = cos_elevation * cos_azimuth;
  dir[1] = cos_elevation * sin_azimuth;
  dir[2] = sin_elevation;
}

static res_T
setup_camera_orthographic
  (struct htrdr_atmosphere* cmd,
   const struct htrdr_atmosphere_args* args)
{
  struct scam_orthographic_args cam_args = SCAM_ORTHOGRAPHIC_ARGS_DEFAULT;
  ASSERT(cmd && args && args->image.definition[0] && args->image.definition[1]);
  ASSERT(cmd->output_type == HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE);
  ASSERT(args->cam_type == HTRDR_ARGS_CAMERA_ORTHOGRAPHIC);

  d3_set(cam_args.position, args->cam_ortho.position);
  d3_set(cam_args.target, args->cam_ortho.target);
  d3_set(cam_args.up, args->cam_ortho.up);
  cam_args.height = args->cam_ortho.height;
  cam_args.aspect_ratio =
    (double)args->image.definition[0]
  / (double)args->image.definition[1];

  return scam_create_orthographic
    (htrdr_get_logger(cmd->htrdr),
     htrdr_get_allocator(cmd->htrdr),
     htrdr_get_verbosity_level(cmd->htrdr),
     &cam_args,
     &cmd->camera);
}

static res_T
setup_camera_perspective
  (struct htrdr_atmosphere* cmd,
   const struct htrdr_atmosphere_args* args)
{
  struct scam_perspective_args cam_args = SCAM_PERSPECTIVE_ARGS_DEFAULT;
  ASSERT(cmd && args && args->image.definition[0] && args->image.definition[1]);
  ASSERT(cmd->output_type == HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE);
  ASSERT(args->cam_type == HTRDR_ARGS_CAMERA_PERSPECTIVE);

  d3_set(cam_args.position, args->cam_persp.position);
  d3_set(cam_args.target, args->cam_persp.target);
  d3_set(cam_args.up, args->cam_persp.up);
  cam_args.aspect_ratio =
    (double)args->image.definition[0]
  / (double)args->image.definition[1];
  cam_args.field_of_view = MDEG2RAD(args->cam_persp.fov_y);
  cam_args.lens_radius = args->cam_persp.lens_radius;
  cam_args.focal_distance = args->cam_persp.focal_dst;

  return scam_create_perspective
    (htrdr_get_logger(cmd->htrdr),
     htrdr_get_allocator(cmd->htrdr),
     htrdr_get_verbosity_level(cmd->htrdr),
     &cam_args,
     &cmd->camera);
}

static res_T
setup_camera
  (struct htrdr_atmosphere* cmd,
   const struct htrdr_atmosphere_args* args)
{
  res_T res = RES_OK;
  ASSERT(cmd->output_type == HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE);
  switch(args->cam_type) {
    case HTRDR_ARGS_CAMERA_ORTHOGRAPHIC:
      res = setup_camera_orthographic(cmd, args);
      break;
    case HTRDR_ARGS_CAMERA_PERSPECTIVE:
      res = setup_camera_perspective(cmd, args);
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
  return res;
}

static res_T
setup_flux_map
  (struct htrdr_atmosphere* cmd,
   const struct htrdr_atmosphere_args* args)
{
  ASSERT(cmd && args);
  ASSERT(cmd->output_type == HTRDR_ATMOSPHERE_ARGS_OUTPUT_FLUX_MAP);

  if(args->spectral.spectral_type == HTRDR_SPECTRAL_SW_CIE_XYZ) {
    htrdr_log_err(cmd->htrdr,
      "the CIE 1931 XYZ spectral integration can be used only with a camera"
      "sensor.\n");
    return RES_BAD_ARG;
  }

  return htrdr_rectangle_create
    (cmd->htrdr,
     args->flux_map.size,
     args->flux_map.position,
     args->flux_map.target,
     args->flux_map.up,
     &cmd->flux_map);
}

static res_T
setup_sensor
  (struct htrdr_atmosphere* cmd,
   const struct htrdr_atmosphere_args* args)
{
  res_T res = RES_OK;
  switch(cmd->output_type) {
    case HTRDR_ATMOSPHERE_ARGS_OUTPUT_FLUX_MAP:
      res = setup_flux_map(cmd, args);
      break;
    case HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE:
      res = setup_camera(cmd, args);
      break;
    default: /* Nothing to do */ break;
  }
  return res;
}

static res_T
dump_volumetric_acceleration_structure(struct htrdr_atmosphere* cmd)
{
  size_t nbands;
  size_t i;
  res_T res = RES_OK;
  ASSERT(cmd);

  nbands = htsky_get_spectral_bands_count(cmd->sky);

  /* Nothing to do */
  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) goto exit;

  FOR_EACH(i, 0, nbands) {
    size_t iquad;
    const size_t iband = htsky_get_spectral_band_id(cmd->sky, i);
    const size_t nquads = htsky_get_spectral_band_quadrature_length
      (cmd->sky, iband);

    FOR_EACH(iquad, 0, nquads) {
      res = htsky_dump_cloud_vtk(cmd->sky, iband, iquad, cmd->output);
      if(res != RES_OK) goto error;
      fprintf(cmd->output, "---\n");
    }
  }

exit:
  return res;
error:
  goto exit;
}

static void
atmosphere_release(ref_T* ref)
{
  struct htrdr_atmosphere* cmd = CONTAINER_OF(ref, struct htrdr_atmosphere, ref);
  struct htrdr* htrdr = NULL;
  ASSERT(ref);

  if(cmd->ground) htrdr_atmosphere_ground_ref_put(cmd->ground);
  if(cmd->mats) htrdr_materials_ref_put(cmd->mats);
  if(cmd->sun) htrdr_atmosphere_sun_ref_put(cmd->sun);
  if(cmd->cie) htrdr_cie_xyz_ref_put(cmd->cie);
  if(cmd->ran_wlen) htrdr_ran_wlen_ref_put(cmd->ran_wlen);
  if(cmd->camera) SCAM(ref_put(cmd->camera));
  if(cmd->flux_map) htrdr_rectangle_ref_put(cmd->flux_map);
  if(cmd->buf) htrdr_buffer_ref_put(cmd->buf);
  if(cmd->sky) HTSKY(ref_put(cmd->sky));
  if(cmd->output && cmd->output != stdout) fclose(cmd->output);
  str_release(&cmd->output_name);

  htrdr = cmd->htrdr;
  MEM_RM(htrdr_get_allocator(htrdr), cmd);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Exported functions
 ******************************************************************************/
res_T
htrdr_atmosphere_create
  (struct htrdr* htrdr,
   const struct htrdr_atmosphere_args* args,
   struct htrdr_atmosphere** out_cmd)
{
  struct htrdr_atmosphere* cmd = NULL;
  struct htsky_args htsky_args = HTSKY_ARGS_DEFAULT;
  double sun_dir[3];
  double spectral_range[2];
  const char* output_name = NULL;
  size_t nintervals; /* #bands used to discretized the spectral curve */
  res_T res = RES_OK;
  ASSERT(htrdr && args && out_cmd);

  cmd = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*cmd));
  if(!cmd) {
    htrdr_log_err(htrdr,
      "%s: could not allocate the htrdr_atmosphere data.\n", FUNC_NAME);
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&cmd->ref);
  str_init(htrdr_get_allocator(htrdr), &cmd->output_name);
  cmd->output_type = args->output_type;
  cmd->verbose = args->verbose;
  cmd->spp = args->image.spp;
  cmd->width = args->image.definition[0];
  cmd->height = args->image.definition[1];
  cmd->grid_max_definition[0] = args->grid_max_definition[0];
  cmd->grid_max_definition[1] = args->grid_max_definition[1];
  cmd->grid_max_definition[2] = args->grid_max_definition[2];
  cmd->spectral_type = args->spectral.spectral_type;
  cmd->ref_temperature = args->spectral.ref_temperature;
  cmd->sky_mtl_name = args->sky_mtl_name;

  /* Get ownership on the htrdr structure */
  htrdr_ref_get(htrdr);
  cmd->htrdr = htrdr;

  if(!args->filename_output) {
    cmd->output = stdout;
    output_name = "<stdout>";
  } else if(htrdr_get_mpi_rank(htrdr) != 0) {
    cmd->output = NULL;
    output_name = "<null>";
  } else {
    res = htrdr_open_output_stream(htrdr, args->filename_output, 0/*read*/,
      args->force_overwriting, &cmd->output);
    if(res != RES_OK) goto error;
    output_name = args->filename_output;
  }
  res = str_set(&cmd->output_name, output_name);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: could not store the name of the output stream `%s' -- %s.\n",
      FUNC_NAME, output_name, res_to_cstr(res));
    goto error;
  }

  /* Materials are necessary only if a ground geometry is defined */
  if(args->filename_obj) {
    res = htrdr_materials_create(htrdr, args->filename_mtl, &cmd->mats);
    if(res != RES_OK) goto error;
  }

  res = htrdr_atmosphere_ground_create(htrdr, args->filename_obj, cmd->mats,
    args->repeat_ground, &cmd->ground);
  if(res != RES_OK) goto error;

  res = setup_sensor(cmd, args);
  if(res != RES_OK) goto error;

  res = htrdr_atmosphere_sun_create(cmd->htrdr, &cmd->sun);
  if(res != RES_OK) goto error;
  spherical_to_cartesian_dir
    (MDEG2RAD(args->sun_azimuth), MDEG2RAD(args->sun_elevation), sun_dir);
  htrdr_atmosphere_sun_set_direction(cmd->sun, sun_dir);

  htsky_args.htcp_filename = args->filename_les;
  htsky_args.htgop_filename = args->filename_gas;
  htsky_args.htmie_filename = args->filename_mie;
  htsky_args.cache_filename = args->filename_cache;
  htsky_args.grid_max_definition[0] = args->grid_max_definition[0];
  htsky_args.grid_max_definition[1] = args->grid_max_definition[1];
  htsky_args.grid_max_definition[2] = args->grid_max_definition[2];
  htsky_args.optical_thickness = args->optical_thickness;
  htsky_args.nthreads = (unsigned)htrdr_get_threads_count(htrdr);
  htsky_args.repeat_clouds = args->repeat_clouds;
  htsky_args.verbose = htrdr_get_mpi_rank(htrdr) == 0 ? args->verbose : 0;
  htsky_args.spectral_type = htrdr_to_sky_spectral_type(args->spectral.spectral_type);
  htsky_args.wlen_range[0] = args->spectral.wlen_range[0];
  htsky_args.wlen_range[1] = args->spectral.wlen_range[1];
  res = htsky_create(htrdr_get_logger(htrdr), htrdr_get_allocator(htrdr),
    &htsky_args, &cmd->sky);
  if(res != RES_OK) goto error;

  HTSKY(get_raw_spectral_bounds(cmd->sky, spectral_range));

  spectral_range[0] = MMAX(args->spectral.wlen_range[0], spectral_range[0]);
  spectral_range[1] = MMIN(args->spectral.wlen_range[1], spectral_range[1]);
  if(spectral_range[0] != args->spectral.wlen_range[0]
  || spectral_range[1] != args->spectral.wlen_range[1]) {
    htrdr_log_warn(htrdr,
      "%s: the submitted spectral range overflowed the spectral data.\n",
      FUNC_NAME);
  }

  cmd->wlen_range_m[0] = spectral_range[0]*1e-9; /* Convert in meters */
  cmd->wlen_range_m[1] = spectral_range[1]*1e-9; /* Convert in meters */

  /* Compute the number of fixed sized bands used to descrised to the spectral
   * data */
  nintervals = compute_spectral_bands_count(cmd);

  if(cmd->spectral_type == HTRDR_SPECTRAL_SW_CIE_XYZ) {
    res = htrdr_cie_xyz_create(htrdr, spectral_range, nintervals, &cmd->cie);
    if(res != RES_OK) goto error;
  } else {
    if(cmd->ref_temperature <= 0) {
      htrdr_log_err(htrdr, "%s: invalid reference temperature %g K.\n",
        FUNC_NAME, cmd->ref_temperature);
      res = RES_BAD_ARG;
      goto error;
    }
    res = htrdr_ran_wlen_create
      (htrdr, spectral_range, nintervals, cmd->ref_temperature, &cmd->ran_wlen);
    if(res != RES_OK) goto error;
  }

  if(cmd->output_type != HTRDR_ATMOSPHERE_ARGS_OUTPUT_OCTREES) {
    struct htrdr_pixel_format pixfmt = HTRDR_PIXEL_FORMAT_NULL;
    atmosphere_get_pixel_format(cmd, &pixfmt);

    /* Setup the buffer layout */
    cmd->buf_layout.width = args->image.definition[0];
    cmd->buf_layout.height = args->image.definition[1];
    cmd->buf_layout.pitch = args->image.definition[0] * pixfmt.size;
    cmd->buf_layout.elmt_size = pixfmt.size;
    cmd->buf_layout.alignment = pixfmt.alignment;

    /* Create the image buffer only on the master process; the image parts
     * rendered by the others processes are gathered onto the master process */
    if(htrdr_get_mpi_rank(htrdr) == 0) {
      res = htrdr_buffer_create(htrdr, &cmd->buf_layout, &cmd->buf);
      if(res != RES_OK) goto error;
    }
  }

exit:
  *out_cmd = cmd;
  return res;
error:
  if(cmd) {
    htrdr_atmosphere_ref_put(cmd);
    cmd = NULL;
  }
  goto exit;
}

void
htrdr_atmosphere_ref_get(struct htrdr_atmosphere* cmd)
{
  ASSERT(cmd);
  ref_get(&cmd->ref);
}

void
htrdr_atmosphere_ref_put(struct htrdr_atmosphere* cmd)
{
  ASSERT(cmd);
  ref_put(&cmd->ref, atmosphere_release);
}

res_T
htrdr_atmosphere_run(struct htrdr_atmosphere* cmd)
{
  res_T res = RES_OK;
  switch(cmd->output_type) {
    case HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE:
    case HTRDR_ATMOSPHERE_ARGS_OUTPUT_FLUX_MAP:
      res = atmosphere_draw_map(cmd);
      break;
    case HTRDR_ATMOSPHERE_ARGS_OUTPUT_OCTREES:
      res = dump_volumetric_acceleration_structure(cmd);
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
  return res;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
void
atmosphere_get_pixel_format
  (const struct htrdr_atmosphere* cmd,
   struct htrdr_pixel_format* fmt)
{
  ASSERT(cmd && fmt);
  switch(cmd->output_type) {
    case HTRDR_ATMOSPHERE_ARGS_OUTPUT_FLUX_MAP:
      fmt->size = sizeof(struct atmosphere_pixel_flux);
      fmt->alignment = ALIGNOF(struct atmosphere_pixel_flux);
      break;
    case HTRDR_ATMOSPHERE_ARGS_OUTPUT_IMAGE:
      switch(cmd->spectral_type) {
        case HTRDR_SPECTRAL_LW:
        case HTRDR_SPECTRAL_SW:
          fmt->size = sizeof(struct atmosphere_pixel_xwave);
          fmt->alignment = ALIGNOF(struct atmosphere_pixel_xwave);
          break;
        case HTRDR_SPECTRAL_SW_CIE_XYZ:
          fmt->size = sizeof(struct atmosphere_pixel_image);
          fmt->alignment = ALIGNOF(struct atmosphere_pixel_image);
          break;
        default: FATAL("Unreachable code.\n"); break;
      }
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
}

