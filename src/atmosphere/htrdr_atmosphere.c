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

#include "htrdr_atmosphere.h"
#include "htrdr_atmosphere_args.h"
#include "htrdr_atmosphere_sun.h"

#include "htrdr_buffer.h"
#include "htrdr_camera.h"
#include "htrdr_cie_xyz.h"
#include "htrdr_geometry.h"
#include "htrdr_log.h"
#include "htrdr_materials.h"

#include <star/s3d.h>

struct pixel_format {
  size_t size; /* In bytes */
  size_t alignment; /* In bytes */
};
#define PIXEL_FORMAT_NULL__ {0, 0}
static const struct pixel_format PIXEL_FORMAT_NULL = PIXEL_FORMAT_NULL__;

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
  double wlen_size;
  const size_t nbands;
  const double band_len;
  const double band_len_max;
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
  return nbands
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

static INLINE void
get_pixel_format
  (const struct htrdr_atmosphere* cmd,
   struct pixel_format* fmt)
{
  switch(cmd->sensor_type) {
    case HTRDR_SENSOR_RECTANGLE:
      fmt->size = sizeof(struct atmosphere_pixel_flux);
      fmt->alignment = ALIGNOF(struct atmosphere_pixel_flux);
      break;
    case HTRDR_SENSOR_CAMERA:
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

static res_T
setup_sensor
  (struct htrdr_atmosphere* cmd,
   const struct htrdr_atmosphere_args* args)
{
  double proj_ratio;
  res_T res = RES_OK;
  ASSERT(htrdr && args);

  cmd->sensor.type = args->sensor_type;

  if(args->spectral.type == HTRDR_SPECTRAL_SW_CIE_XYZ
  && args->sensor_type != HTRDR_SENSOR_CAMERA) {
    htrdr_log_err(cmd->htrdr, "the CIE 1931 XYZ spectral integration can be used "
      "only with a camera sensor.\n");
    res = RES_BAD_ARG;
    goto error;
  }

  switch(args->sensor_type) {
    case HTRDR_SENSOR_CAMERA:
      proj_ratio =
        (double)args->image.definition[0]
      / (double)args->image.definition[1];
      res = htrdr_camera_create(htrdr, args->camera.position,
        args->camera.target, args->camera.up, proj_ratio,
        MDEG2RAD(args->camera.fov_y), &cmd->sensor.camera);
      break;
    case HTRDR_SENSOR_RECTANGLE:
      res = htrdr_rectangle_create(htrdr, args->rectangle.size,
        args->rectangle.position, args->rectangle.target, args->rectangle.up,
        &cmd->sensor.rectangle);
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  goto exit;
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
  if(htrdr->mpi_rank != 0) goto exit;

  FOR_EACH(i, 0, nbands) {
    size_t iquad;
    const size_t iband = htsky_get_spectral_band_id(cmd->sky, i);
    const size_t nquads = htsky_get_spectral_band_quadrature_length
      (cmd->sky, iband);

    FOR_EACH(iquad, 0, nquads) {
      res = htsky_dump_cloud_vtk(htrdr->sky, iband, iquad, htrdr->output);
      if(res != RES_OK) goto error;
      fprintf(htrdr->output, "---\n");
    }
  }

exit:
  return res;
error:
  goto exit;
}

static INLINE void
dump_accum
  (const struct htrdr_accum* acc, /* Accum to dump */
   struct htrdr_accum* out_acc, /* May be NULL */
   FILE* stream)
{
  ASSERT(acc && stream);

  if(acc->nweights == 0) {
    fprintf(streamn, "0 0 ");
  } else {
    struct htrdr_estimate estimate = HTRDR_ESTIMATE_NULL;

    htrdr_accum_get_estimation(acc, &estimate);
    fprintf("%g %g\n", estimate.E, estimate.SE);

    if(out_acc) {
      out_acc->sum_weights += acc->sum_weights;
      out_acc->sum_weights_sqr += acc->sum_weights_sqr;
      out_acc->nweights += acc->nweights;
    }
  }
}

static INLINE void
dump_pixel_flux
  (const struct htrdr_pixel_flux* pix,
   struct htrdr_accum* time_acc, /* May be NULL */
   struct htrdr_accum* flux_acc, /* May be NULL */
   FILE* stream)
{
  struct htrdr_estimate pix_time = HTRDR_ESTIMATE_NULL;

  ASSERT(pix && stream);
  dump_accum(&pix->flux, flux_acc, stream_name, stream);
  fprintf(stream, "0 0 0 0 ");
  dump_accum(&pix->time, time_acc, stream_name, stream);
  fprintf(stream, "\n");
}

static INLINE void
dump_pixel_image
  (const struct htrdr_pixel_image* pix,
   struct htrdr_accum* time_acc, /* May be NULL */
   FILE* stream)
{
  ASSERT(pix && stream_name && stream);
  fprintf(stream, "%g %g ", pix->X.E, pix->X.SE);
  fprintf(stream, "%g %g ", pix->Y.E, pix->Y.SE);
  fprintf(stream, "%g %g ", pix->Z.E, pix->Z.SE);
  dump_accum(pix->time, time_acc, stream);
  fprintf(stream, "\n");
}

static INLINE void
dump_pixel_xwave
  (const struct htrdr_pixel_xwave* pix,
   struct htrdr_accum* time_acc, /* May be NULL */
   const char* stream_name,
   FILE* stream)
{
  fprintf(stream, "%g %g %f %f 0 0 ",
    pix->radiance_temperature.E,
    pix->radiance_temperature.SE,
    pix->radiance.E,
    pix->radiance.SE);
  dump_accum(pix->time, time_acc, stream);
  fprintf(stream, "\n");
}

static res_T
dump_buffer
  (struct htrdr_atmosphere* cmd,
   struct htrdr_buffer* buf,
   struct htrdr_accum* time_acc, /* May be NULL */
   struct htrdr_accum* flux_acc, /* May be NULL */
   const char* stream_name,
   FILE* stream)
{
  struct htrdr_buffer_layout layout;
  size_t pixsz, pixal;
  size_t x, y;
  res_T res = RES_OK;
  ASSERT(cmd && buf && stream_name && stream);
  (void)stream_name;

  pixsz = htrdr_spectral_type_get_pixsz(cmd->spectral_type, cmd->sensor.type);
  pixal = htrdr_spectral_type_get_pixal(cmd->spectral_type, cmd->sensor.type);

  htrdr_buffer_get_layout(buf, &layout);
  if(layout.elmt_size != pixsz || layout.alignment != pixal) {
    htrdr_log_err(cmd->htrdr, "%s: invalid buffer layout. ", FUNC_NAME);
    res = RES_BAD_ARG;
    goto error;
  }

  fprintf(stream, "%lu %lu\n", layout.width, layout.height);

  if(time_acc) *time_acc = HTRDR_ACCUM_NULL;
  if(flux_acc) *flux_acc = HTRDR_ACCUM_NULL;

  FOR_EACH(y, 0, layout.height) {
    FOR_EACH(x, 0, layout.width) {
      struct htrdr_estimate pix_time = HTRDR_ESTIMATE_NULL;
      const struct htrdr_accum* pix_time_acc = NULL;

      if(cmd->sensor.type == HTRDR_SENSOR_RECTANGLE) {
        const struct htrdr_pixel_flux* pix = htrdr_buffer_at(buf, x, y);
        dump_pixel_flux(pix, time_acc, flux_acc, stream_name, stream);
      } else if(cmd->spectral_type == HTRDR_SPECTRAL_SW_CIE_XYZ) {
        const struct htrdr_pixel_image* pix = htrdr_buffer_at(buf, x, y);
        dump_pixel_image(pix, time_acc, stream_name, stream);
      } else {
        const struct htrdr_pixel_xwave* pix =  htrdr_buffer_at(buf, x, y);
        dump_pixel_xwave(pix, time_acc, stream_name, stream);
      }
    }
    fprintf(stream, "\n");
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

  if(cmd->s3d) S3D(device_ref_put(cmd->s3d));
  if(cmd->ground) htrdr_atmosphere_ground_ref_put(cmd->ground);
  if(cmd->mats) htrdr_materials_ref_put(cmd->mats);
  if(cmd->sun) htrdr_atmosphere_sun_ref_put(cmd->sun);
  if(cmd->cie) htrdr_cie_xyz_ref_put(cmd->cie);
  if(cmd->ran_wlen) htrdr_ran_wlen_ref_put(cmd->ran_wlen);
  if(cmd->sensor.camera) htrdr_camera_ref_put(cmd->sensor.camera);
  if(cmd->sensor.rectangle) htrdr_rectangle_ref_put(cmd->sensor.rectangle);
  if(cmd->but) htrdr_buffer_ref_put(cmd->but);
  if(cmd->sky) HTSKY(ref_put(cmd->sky));
  if(cmd->output && cmd->output != stdout) fclose(cmd->output);
  str_release(cmd->output_name);

  htrdr = cmd->htrdr;
  MEM_RM(htrdr->allocator, cmd);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Local functions
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

  cmd = MEM_CALLOC(htrdr->allocator, 1, sizeof(*cmd));
  if(!cmd) {
    htrdr_log_err(htrdr,
      "%s: could not allocate the htrdr_atmosphere data.\n", FUNC_NAME);
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&cmd->ref);
  str_init(htrdr->allocator, &cmd->output_name);
  cmd->dump_vtk = args->dump_vtk;
  cmd->verbose = args->verbose;
  cmd->spp = args->image.spp;
  cmd->width = args->image.definition[0];
  cmd->height = args->image.definition[1];
  cmd->grid_max_definition[0] = args->grid_max_definition[0];

  /* Get ownership onf the htrdr structure */
  htrdr_ref_get(htrdr);
  cmd->htrdr = htrdr;

  if(!args->filename_output) {
    cmd->output = stdout;
    output_name = "<stdout>";
  } else if(htrdr->mpi_rank != 0) {
    cmd->output = NULL;
    output_name = "<null>";
  } else {
    res = htrdr_open_output_stream
      (htrdr, args->output, 0/*read*/, args->force_overwriting, &cmd->output);
    if(res != RES_OK) goto error;
    output_name = args->output;
  }
  res = str_set(&cmd->output_name, output_name);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: could not store the name of the output stream `%s' -- %s.\n",
      FUNC_NAME, output_name, res_to_cstr(res));
    goto error;
  }

  /* Disable the Star-3D verbosity since the Embree backend prints some messages
   * on stdout rather than stderr. This is annoying since stdout may be used by
   * htrdr to write output data */
  res = s3d_device_create
    (&htrdr->logger, htrdr->allocator, 0, &cmd->s3d);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: could not create the Star-3D device -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }

  /* Materials are necessary only if a ground geometry is defined */
  if(args->filename_obj) {
    res = htrdr_materials_create(htrdr, args->filename_mtl, &cmd->mats);
    if(res != RES_OK) goto error;
  }

  res = htrdr_atmosphere_ground_create(htrdr, args->filename_obj,
    args->repeat_ground, &cmd->ground);
  if(res != RES_OK) goto error;

  res = setup_sensor(cmd, args);
  if(res != RES_OK) goto error;

  res = htrdr_atmosphere_sun_create(htrdr, &cmd->sun);
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
  htsky_args.nthreads = htrdr->nthreads;
  htsky_args.repeat_clouds = args->repeat_clouds;
  htsky_args.verbose = htrdr->mpi_rank == 0 ? args->verbose : 0;
  htsky_args.spectral_type = htrdr_to_sky_spectral_type(args->spectral_type);
  htsky_args.wlen_range[0] = args->spectral.wlen_range[0];
  htsky_args.wlen_range[1] = args->spectral.wlen_range[1];
  res = htsky_create(&htrdr->logger, htrdr->allocator, &htsky_args, &cmd->sky);
  if(res != RES_OK) goto error;

  HTSKY(get_raw_spectral_bounds(cmd->sky, spectral_range));

  spectral_range[0] = MMAX(args->wlen_range[0], spectral_range[0]);
  spectral_range[1] = MMIN(args->wlen_range[1], spectral_range[1]);
  if(spectral_range[0] != args->wlen_range[0]
  || spectral_range[1] != args->wlen_range[1]) {
    htrdr_log_warn(htrdr,
      "%s: the submitted spectral range overflowed the spectral data.\n", FUNC_NAME);
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

  /* Create the image buffer only on the master process; the image parts
   * rendered by the processes are gathered onto the master process. */
  if(!cmd->dump_volumetric_acceleration_structure && htrdr->mpi_rank == 0) {
    struct pixel_format pixfmt = PIXEL_FORMAT_NULL;
    get_pixel_format(cmd, &pixfmt);

    res = htrdr_buffer_create(htrdr,
      args->image.definition[0], /* Width */
      args->image.definition[1], /* Height */
      args->image.definition[0] * pixfmt.size, /* Pitch */
      pixfmt.size, /* Size of a pixel */
      pixfmt.alignment, /* Alignment of a pixel */
      &cmd->buf);
    if(res != RES_OK) goto error;
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
  ref_get(cmd->ref);
}

void
htrdr_atmosphere_ref_put(struct htrdr_atmosphere* cmd)
{
  ASSERT(cmd);
  ref_put(cmd->ref, atmosphere_release);
}

void
htrdr_atmosphere_run(struct htrdr_atmosphere* cmd)
{
  res_T res = RES_OK;

  if(cmd->dump_volumetric_acceleration_structure) {
    res = dump_volumetric_acceleration_structure(cmd);
    if(res != RES_OK) goto error;
  } else {
    res = draw_map(cmd);
    if(res != RES_OK) goto error;
  }

exit:
  return res;
error:
  goto exit;
}

