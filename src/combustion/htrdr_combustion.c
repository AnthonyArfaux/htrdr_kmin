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

#include "combustion/htrdr_combustion.h"
#include "combustion/htrdr_combustion_args.h"
#include "combustion/htrdr_combustion_c.h"
#include "combustion/htrdr_combustion_laser.h"

#include "core/htrdr.h"
#include "core/htrdr_log.h"
#include "core/htrdr_geometry.h"
#include "core/htrdr_materials.h"
#include "core/htrdr_rectangle.h"

#include <astoria/atrstm.h>

#include <star/scam.h>
#include <star/ssf.h>

#include <rsys/cstr.h>
#include <rsys/double3.h>
#include <rsys/mem_allocator.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
release_phase_functions(struct htrdr_combustion* cmd)
{
  size_t i;
  ASSERT(cmd);

  if(!cmd->phase_functions) return; /* Nothing to release */

  FOR_EACH(i, 0, htrdr_get_threads_count(cmd->htrdr)) {
    if(cmd->phase_functions[i]) {
      SSF(phase_ref_put(cmd->phase_functions[i]));
    }
  }
  MEM_RM(htrdr_get_allocator(cmd->htrdr), cmd->phase_functions);
  cmd->phase_functions = NULL;
}

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
  res = str_set(&cmd->output_name, output_name);
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
setup_simd
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  struct ssf_info ssf_info = SSF_INFO_NULL;
  ASSERT(cmd && args);

  cmd->rdgfa_simd = SSF_SIMD_NONE;

  if(args->phase_func_type != HTRDR_COMBUSTION_ARGS_PHASE_FUNC_RDGFA
  || args->use_simd == 0)
    return RES_OK; /* Nothing to do */

  /* Check SIMD support for the RDG-FA phase function */
  ssf_get_info(&ssf_info);
  if(ssf_info.simd_256) {
    htrdr_log(cmd->htrdr,
      "Use the SIMD-256 instruction set for the RDG-FA phase function.\n");
    cmd->rdgfa_simd = SSF_SIMD_256;
  } else if(ssf_info.simd_128) {
    htrdr_log(cmd->htrdr,
      "Use the SIMD-128 instruction set for the RDG-FA phase function.\n");
    cmd->rdgfa_simd = SSF_SIMD_128;
  } else {
    htrdr_log_warn(cmd->htrdr,
      "Cannot use SIMD for the RDG-FA phase function: the "
      "Star-ScatteringFunction library was compiled without SIMD support.\n");
    cmd->rdgfa_simd = SSF_SIMD_NONE;
  }
  return RES_OK;
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
setup_camera_orthographic
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  struct scam_orthographic_args cam_args = SCAM_ORTHOGRAPHIC_ARGS_DEFAULT;
  ASSERT(cmd && args && args->image.definition[0] && args->image.definition[1]);
  ASSERT(cmd->output_type == HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE);
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
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  struct scam_perspective_args cam_args = SCAM_PERSPECTIVE_ARGS_DEFAULT;
  ASSERT(cmd && args && args->image.definition[0] && args->image.definition[1]);
  ASSERT(cmd->output_type == HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE);
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
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  res_T res = RES_OK;
  ASSERT(cmd->output_type == HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE);
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
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  ASSERT(cmd && args);
  ASSERT(cmd->output_type == HTRDR_COMBUSTION_ARGS_OUTPUT_FLUX_MAP);
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
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  res_T res = RES_OK;
  switch(cmd->output_type) {
    case HTRDR_COMBUSTION_ARGS_OUTPUT_FLUX_MAP:
      res = setup_flux_map(cmd, args);
      break;
    case HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE:
      res = setup_camera(cmd, args);
      break;
    default: /* Nothing to do */ break;
  }
  return res;
}

static res_T
setup_laser
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  struct htrdr_combustion_laser_create_args laser_args =
    HTRDR_COMBUSTION_LASER_CREATE_ARGS_DEFAULT;
  ASSERT(cmd && args);
  cmd->wavelength = args->wavelength;
  laser_args.surface = args->laser;
  laser_args.wavelength = args->wavelength;
  laser_args.flux_density = args->laser_flux_density;
  return htrdr_combustion_laser_create(cmd->htrdr, &laser_args, &cmd->laser);
}

static res_T
setup_phase_functions(struct htrdr_combustion* cmd)
{
  struct mem_allocator* allocator = NULL;
  const struct ssf_phase_type* phase_type = NULL;
  size_t nthreads;
  size_t i;
  res_T res = RES_OK;
  ASSERT(cmd);

  nthreads = htrdr_get_threads_count(cmd->htrdr);
  allocator = htrdr_get_allocator(cmd->htrdr);

  switch(cmd->phase_func_type) {
    case HTRDR_COMBUSTION_ARGS_PHASE_FUNC_ISOTROPIC:
      htrdr_log(cmd->htrdr, "Use an isotropic phase function.\n");
      phase_type = &ssf_phase_hg;
      break;
    case HTRDR_COMBUSTION_ARGS_PHASE_FUNC_RDGFA:
      htrdr_log(cmd->htrdr, "Use the RDG-FA phase function.\n");
      phase_type = &ssf_phase_rdgfa;
      break;
    default: FATAL("Unreachable code.\n"); break;
  }

  /* Allocate the list of per thread phase function */
  cmd->phase_functions = MEM_CALLOC
    (allocator, nthreads, sizeof(*cmd->phase_functions));
  if(!cmd->phase_functions) {
    htrdr_log_err(cmd->htrdr,
      "Could not allocate the per thread RDG-FA phase function.\n");
    res = RES_MEM_ERR;
    goto error;
  }

  /* Create the per thread phase function */
  FOR_EACH(i, 0, nthreads) {
    res = ssf_phase_create(allocator, phase_type, cmd->phase_functions+i);
    if(res != RES_OK) {
      htrdr_log_err(cmd->htrdr,
        "Could not create the phase function for the thread %lu -- %s.\n",
        (unsigned long)i, res_to_cstr(res));
      goto error;
    }
  }

exit:
  return res;
error:
  release_phase_functions(cmd);
  goto exit;
}

static res_T
setup_buffer
  (struct htrdr_combustion* cmd,
   const struct htrdr_combustion_args* args)
{
  struct htrdr_pixel_format pixfmt = HTRDR_PIXEL_FORMAT_NULL;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  if(cmd->output_type != HTRDR_COMBUSTION_ARGS_OUTPUT_FLUX_MAP
  && cmd->output_type != HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE)
    goto exit;

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
  atrstm_args.fractal_prefactor = args->fractal_prefactor;
  atrstm_args.fractal_dimension = args->fractal_dimension;
  atrstm_args.optical_thickness = args->optical_thickness;
  atrstm_args.precompute_normals = args->precompute_normals;
  atrstm_args.use_simd = args->use_simd;
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

static double
compute_laser_mesh_extent(const struct htrdr_combustion* cmd)
{
  double mdm_upp[3];
  double mdm_low[3];
  double laser_dir[3];
  double laser_pos[3];
  double t[2];
  int max_axis;
  ASSERT(cmd);

  /* Retrieve the medium axis aligned bounding box */
  atrstm_get_aabb(cmd->medium, mdm_low, mdm_upp);

  /* Retrieve laser parameters */
  htrdr_combustion_laser_get_position(cmd->laser, laser_pos);
  htrdr_combustion_laser_get_direction(cmd->laser, laser_dir);

  /* Compute the dominant axis of the laser direction */
  max_axis =
     fabs(laser_dir[0]) > fabs(laser_dir[1])
  ? (fabs(laser_dir[0]) > fabs(laser_dir[2]) ? 0 : 2)
  : (fabs(laser_dir[1]) > fabs(laser_dir[2]) ? 1 : 2);

  /* Define the intersection of the laser along its dominant axis with the
   * medium bounds along this axis */
  t[0] = (mdm_low[max_axis] - laser_pos[max_axis]) / laser_dir[max_axis];
  t[1] = (mdm_upp[max_axis] - laser_pos[max_axis]) / laser_dir[max_axis];
  if(t[0] > t[1]) SWAP(double, t[0], t[1]);

  /* Use the far intersection distance as the extent of the laser mesh */
  return t[1];
}

static res_T
dump_laser_sheet(const struct htrdr_combustion* cmd)
{
  struct htrdr_combustion_laser_mesh laser_mesh;
  double extent;
  unsigned i;
  res_T res = RES_OK;
  ASSERT(cmd);

  htrdr_log(cmd->htrdr, "Write laser sheet to '%s'.\n",
    str_cget(&cmd->output_name));

  /* Compute the extent of the geometry that will represent the laser sheet */
  extent = compute_laser_mesh_extent(cmd);

  /* Retreive the mesh of the laser sheet */
  htrdr_combustion_laser_get_mesh(cmd->laser, extent, &laser_mesh);

  #define FPRINTF(Fmt, Args) {                                                 \
    const int err = fprintf(cmd->output, Fmt COMMA_##Args LIST_##Args);        \
    if(err < 0) {                                                              \
      htrdr_log_err(cmd->htrdr, "Error writing data to `%s'.\n",               \
        str_cget(&cmd->output_name));                                          \
      res = RES_IO_ERR;                                                        \
      goto error;                                                              \
    }                                                                          \
  } (void)0

  /* Write header */
  FPRINTF("# vtk DataFile Version 2.0\n", ARG0());
  FPRINTF("Laser sheet\n", ARG0());
  FPRINTF("ASCII\n", ARG0());
  FPRINTF("DATASET POLYDATA\n", ARG0());

  /* Write the vertices */
  FPRINTF("POINTS %u double\n", ARG1(laser_mesh.nvertices));
  FOR_EACH(i, 0, laser_mesh.nvertices) {
    FPRINTF("%g %g %g\n", ARG3
      (laser_mesh.vertices[i*3+0],
       laser_mesh.vertices[i*3+1],
       laser_mesh.vertices[i*3+2]));
  }

  /* Write the triangles */
  FPRINTF("POLYGONS %u %u\n",ARG2
    (laser_mesh.ntriangles,
     laser_mesh.ntriangles*4));
  FOR_EACH(i, 0, laser_mesh.ntriangles) {
    FPRINTF("3 %u %u %u\n", ARG3
      (laser_mesh.triangles[i*3+0],
       laser_mesh.triangles[i*3+1],
       laser_mesh.triangles[i*3+2]));
  }

  /* Write flux density */
  FPRINTF("CELL_DATA %u\n", ARG1(laser_mesh.ntriangles));
  FPRINTF("SCALARS Flux_density double 1\n", ARG0());
  FPRINTF("LOOKUP_TABLE default\n", ARG0());
  FOR_EACH(i, 0, laser_mesh.ntriangles) {
    FPRINTF("%g\n", ARG1(htrdr_combustion_laser_get_flux_density(cmd->laser)));
  }
  #undef FPRINTF

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
  if(cmd->camera) SCAM(ref_put(cmd->camera));
  if(cmd->flux_map) htrdr_rectangle_ref_put(cmd->flux_map);
  if(cmd->laser) htrdr_combustion_laser_ref_put(cmd->laser);
  if(cmd->buf) htrdr_buffer_ref_put(cmd->buf);
  if(cmd->output && cmd->output != stdout) CHK(fclose(cmd->output) == 0);
  release_phase_functions(cmd);
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
    res = RES_MEM_ERR;
    goto error;
  }
  ref_init(&cmd->ref);
  str_init(htrdr_get_allocator(htrdr), &cmd->output_name);

  /* Get the ownership on the htrdr structure */
  htrdr_ref_get(htrdr);
  cmd->htrdr = htrdr;

  cmd->spp = args->image.spp;
  cmd->output_type = args->output_type;
  cmd->phase_func_type = args->phase_func_type;

  res = setup_output(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_phase_functions(cmd);
  if(res != RES_OK) goto error;
  res = setup_simd(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_geometry(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_sensor(cmd, args);
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

  switch(cmd->output_type) {
    case HTRDR_COMBUSTION_ARGS_OUTPUT_FLUX_MAP:
      res = combustion_draw_map(cmd);
      break;
    case HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE:
      res = combustion_draw_map(cmd);
      break;
    case HTRDR_COMBUSTION_ARGS_OUTPUT_LASER_SHEET:
      res = dump_laser_sheet(cmd);
      break;
    case HTRDR_COMBUSTION_ARGS_OUTPUT_OCTREES:
      res = dump_volumetric_acceleration_structure(cmd);
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
  if(res != RES_OK) {
    goto error;
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
  switch(cmd->output_type) {
    case HTRDR_COMBUSTION_ARGS_OUTPUT_FLUX_MAP:
      fmt->size = sizeof(struct combustion_pixel_flux);
      fmt->alignment = ALIGNOF(struct combustion_pixel_flux);
      break;
    case HTRDR_COMBUSTION_ARGS_OUTPUT_IMAGE:
      fmt->size = sizeof(struct combustion_pixel_image);
      fmt->alignment = ALIGNOF(struct combustion_pixel_image);
      break;
    default: FATAL("Unreachable code.\n"); break;
  }
}
