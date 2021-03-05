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

#include "combustion/htrdr_combustion.h"
#include "combustion/htrdr_combustion_args.h"
#include "combustion/htrdr_combustion_c.h"

#include "core/htrdr.h"
#include "core/htrdr_camera.h"
#include "core/htrdr_log.h"
#include "core/htrdr_geometry.h"
#include "core/htrdr_materials.h"
#include "core/htrdr_rectangle.h"

#include <astoria/atrstm.h>

#include <rsys/cstr.h>
#include <rsys/mem_allocator.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static res_T
setup_output
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  const char* output_name = NULL;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) {
    /* No output stream on non master processes */
    cmd->output = NULL;
    output_name = "<null>";

  } else if(!args->path_output) {
    /* Write results to standard output when no destination file is defined */
    cmd->output = stdout;
    output_name = "<stdout>";

  } else {
    /* Open the output stream */
    res = htrdr_open_output_stream(cmd->htrdr, args->path_output,
      0/*read*/, args->force_overwriting, &cmd->output);
    if(res != RES_OK) goto error;
    output_name = args->path_output;
  }

  /* Setup the output name */
  str_set(&cmd->output_name, output_name);
  if(res != RES_OK) {
    htrdr_log_err(cmd->htrdr,
      "Could not store the name of the output stream `%s' -- %s.\n",
      output_name, res_to_cstr(res));
    goto error;
  }

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

static res_T
setup_geometry
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  res_T res = RES_OK;

  if(!args->geom.path_obj) goto exit;
  ASSERT(args->geom.path_mats);

  res = htrdr_materials_create(cmd->htrdr, args->geom.path_mats, &cmd->mats);
  if(res != RES_OK) goto error;

  res = htrdr_geometry_create
    (cmd->htrdr, args->geom.path_obj, cmd->mats, &cmd->geom);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  if(cmd->mats) { htrdr_materials_ref_put(cmd->mats); cmd->mats = NULL; }
  if(cmd->geom) { htrdr_geometry_ref_put(cmd->geom); cmd->geom = NULL; }
  goto exit;
}

static res_T
setup_camera
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  double proj_ratio = 0;
  ASSERT(cmd && args && args->image.definition[0] && args->image.definition[1]);

  proj_ratio =
    (double)args->image.definition[0]
  / (double)args->image.definition[1];

  return htrdr_camera_create
    (cmd->htrdr,
     args->camera.position,
     args->camera.target,
     args->camera.up,
     proj_ratio,
     MDEG2RAD(args->camera.fov_y),
     &cmd->camera);
}

static res_T
setup_laser
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  ASSERT(cmd && args);
  cmd->wavelength = args->wavelength;
  return htrdr_rectangle_create
    (cmd->htrdr,
     args->laser.size,
     args->laser.position,
     args->laser.target,
     args->laser.up,
     &cmd->laser);
}

static res_T
setup_buffer
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  struct htrdr_pixel_format pixfmt = HTRDR_PIXEL_FORMAT_NULL;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  if(cmd->dump_volumetric_acceleration_structure) goto exit;

  combustion_get_pixel_format(cmd, &pixfmt);

  /* Setup the buffer layout */
  cmd->buf_layout.width = args->image.definition[0];
  cmd->buf_layout.height = args->image.definition[1];
  cmd->buf_layout.pitch = args->image.definition[0] * pixfmt.size;
  cmd->buf_layout.elmt_size = pixfmt.size;
  cmd->buf_layout.alignment = pixfmt.alignment;

  /* Create the image buffer only on the master process; the image parts
   * rendered by the others processes are gathered onto the master process */
  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) goto exit;

  res = htrdr_buffer_create(cmd->htrdr, &cmd->buf_layout, &cmd->buf);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  if(cmd->buf) { htrdr_buffer_ref_put(cmd->buf); cmd->buf = NULL; }
  goto exit;
}

static res_T
setup_medium
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  struct atrstm_args atrstm_args = ATRSTM_ARGS_DEFAULT;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  /* Setup the semi-transaprent medium arguments */
  atrstm_args.sth_filename = args->path_tetra;
  atrstm_args.atrtp_filename = args->path_therm_props;
  atrstm_args.atrri_filename = args->path_refract_ids;
  atrstm_args.cache_filename = args->path_cache;
  atrstm_args.spectral_type = ATRSTM_SPECTRAL_SW;
  atrstm_args.wlen_range[0] = args->wavelength;
  atrstm_args.wlen_range[1] = args->wavelength;
  atrstm_args.gyration_radius_prefactor = args->gyration_radius_prefactor;
  atrstm_args.fractal_dimension = args->fractal_dimension;
  atrstm_args.optical_thickness = args->optical_thickness;
  atrstm_args.precompute_normals = args->precompute_normals;
  atrstm_args.nthreads = args->nthreads;
  atrstm_args.verbose = args->verbose;

  switch(args->grid.type) {
    case HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_AUTO:
      atrstm_args.auto_grid_definition = 1;
      atrstm_args.auto_grid_definition_hint = args->grid.definition.hint;
      break;
    case HTRDR_COMBUSTION_ARGS_GRID_DEFINITION_FIXED:
      atrstm_args.auto_grid_definition = 0;
      atrstm_args.grid_max_definition[0] = args->grid.definition.fixed[0];
      atrstm_args.grid_max_definition[1] = args->grid.definition.fixed[1];
      atrstm_args.grid_max_definition[2] = args->grid.definition.fixed[2];
      break;
    default: FATAL("Unreachable code.\n"); break;
  }

  /* Here we go! Create the semi-transparent medium */
  res = atrstm_create
    (htrdr_get_logger(cmd->htrdr),
     htrdr_get_allocator(cmd->htrdr),
     &atrstm_args,
     &cmd->medium);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  if(cmd->medium) { ATRSTM(ref_put(cmd->medium)); cmd->medium = NULL; }
  goto exit;
}

static res_T
dump_volumetric_acceleration_structure(struct htrdr_combustion* cmd)
{
  struct atrstm_dump_svx_octree_args args = ATRSTM_DUMP_SVX_OCTREE_ARGS_DEFAULT;
  res_T res = RES_OK;
  ASSERT(cmd);

  /* Nothing to do on non master process */
  if(htrdr_get_mpi_rank(cmd->htrdr) != 0) goto exit;

  htrdr_log(cmd->htrdr, "Write volumetric acceleration structure to '%s'.\n",
    str_cget(&cmd->output_name));

  res = atrstm_dump_svx_octree(cmd->medium, &args, cmd->output);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  goto exit;
}

static void
combustion_release(ref_T* ref)
{
  struct htrdr_combustion* cmd = CONTAINER_OF(ref, struct htrdr_combustion, ref);
  struct htrdr* htrdr = NULL;
  ASSERT(ref);

  if(cmd->geom) htrdr_geometry_ref_put(cmd->geom);
  if(cmd->mats) htrdr_materials_ref_put(cmd->mats);
  if(cmd->medium) ATRSTM(ref_put(cmd->medium));
  if(cmd->camera) htrdr_camera_ref_put(cmd->camera);
  if(cmd->laser) htrdr_rectangle_ref_put(cmd->laser);
  if(cmd->buf) htrdr_buffer_ref_put(cmd->buf);
  if(cmd->output && cmd->output != stdout) CHK(fclose(cmd->output) == 0);
  str_release(&cmd->output_name);

  htrdr = cmd->htrdr;
  MEM_RM(htrdr_get_allocator(htrdr), cmd);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Exported functions
 ******************************************************************************/
res_T
htrdr_combustion_create
  (struct htrdr* htrdr,
   const struct htrdr_combustion_args* args,
   struct htrdr_combustion** out_cmd)
{
  struct htrdr_combustion* cmd = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && args && out_cmd);

  cmd = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*cmd));
  if(!cmd) {
    htrdr_log_err(htrdr, "Could not allocate the htrdr_combustion data.\n");
    res = RES_BAD_ARG;
    goto error;
  }
  ref_init(&cmd->ref);
  str_init(htrdr_get_allocator(htrdr), &cmd->output_name);

  /* Get the ownership on the htrdr structure */
  htrdr_ref_get(htrdr);
  cmd->htrdr = htrdr;

  cmd->spp = args->image.spp;
  cmd->dump_volumetric_acceleration_structure =
    args->dump_volumetric_acceleration_structure;

  res = setup_output(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_geometry(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_camera(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_laser(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_buffer(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_medium(cmd, args);
  if(res != RES_OK) goto error;

exit:
  *out_cmd = cmd;
  return res;
error:
  if(cmd) {
    htrdr_combustion_ref_put(cmd);
    cmd = NULL;
  }
  goto exit;
}

void
htrdr_combustion_ref_get(struct htrdr_combustion* cmd)
{
  ASSERT(cmd);
  ref_get(&cmd->ref);
}

void
htrdr_combustion_ref_put(struct htrdr_combustion* cmd)
{
  ASSERT(cmd);
  ref_put(&cmd->ref, combustion_release);
}

res_T
htrdr_combustion_run(struct htrdr_combustion* cmd)
{
  res_T res = RES_OK;
  ASSERT(cmd);

  if(cmd->dump_volumetric_acceleration_structure) {
    res = dump_volumetric_acceleration_structure(cmd);
    if(res != RES_OK) goto error;
  } else {
    res = combustion_draw_map(cmd);
    if(res != RES_OK) goto error;
  }

exit:
  return res;
error:
  goto exit;
}

/*******************************************************************************
 * Local functions
 ******************************************************************************/
void
combustion_get_pixel_format
  (const struct htrdr_combustion* cmd,
   struct htrdr_pixel_format* fmt)
{
  ASSERT(cmd && fmt);
  (void)cmd;
  fmt->size = sizeof(struct combustion_pixel);
  fmt->alignment = ALIGNOF(struct combustion_pixel);
}
