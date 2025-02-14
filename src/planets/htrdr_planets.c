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

#define _POSIX_C_SOURCE 200112L /* fdopen, nextafter, rint */

#include "core/htrdr.h"
#include "core/htrdr_ran_wlen_cie_xyz.h"
#include "core/htrdr_ran_wlen_discrete.h"
#include "core/htrdr_ran_wlen_planck.h"
#include "core/htrdr_log.h"

#include "planets/htrdr_planets.h"
#include "planets/htrdr_planets_args.h"
#include "planets/htrdr_planets_c.h"
#include "planets/htrdr_planets_source.h"

#include <rad-net/rnatm.h>
#include <rad-net/rngrd.h>

#include <star/scam.h>
#include <star/smsh.h>

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
/* Calculate the number of fixed size spectral intervals to use for the
 * cumulative */
static size_t
compute_nintervals_for_spectral_cdf(const struct htrdr_planets* cmd)
{
  double range_size;
  size_t nintervals;
  ASSERT(cmd);

  range_size =
    cmd->spectral_domain.wlen_range[1]
  - cmd->spectral_domain.wlen_range[0];

  /* Initially assume ~one interval per nanometer */
  nintervals = (size_t)rint(range_size);

  return nintervals;
}

static res_T
setup_octree_storage
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args,
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
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args)
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
  rnatm_args.precompute_normals = args->precompute_normals;
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
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args)
{
  struct rngrd_create_args rngrd_args = RNGRD_CREATE_ARGS_DEFAULT;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  if(cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_OCTREES)
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

static res_T
setup_spectral_domain_sw
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args)
{
  res_T res = RES_OK;
  ASSERT(cmd && args);
  ASSERT(cmd->spectral_domain.type == HTRDR_SPECTRAL_SW);

  /* Discrete distribution */
  if(args->source.rnrl_filename) {
    struct htrdr_planets_source_spectrum spectrum;
    struct htrdr_ran_wlen_discrete_create_args discrete_args;

    res = htrdr_planets_source_get_spectrum
      (cmd->source, cmd->spectral_domain.wlen_range, &spectrum);
    if(res != RES_OK) goto error;

    discrete_args.get = htrdr_planets_source_spectrum_at;
    discrete_args.nwavelengths = spectrum.size;
    discrete_args.context = &spectrum;
    res = htrdr_ran_wlen_discrete_create
      (cmd->htrdr, &discrete_args, &cmd->discrete);
    if(res != RES_OK) goto error;

  /* Planck distribution */
  } else {
    const size_t nintervals = compute_nintervals_for_spectral_cdf(cmd);

    /* Use the source temperature as the reference temperature of the Planck
     * distribution */
    res = htrdr_ran_wlen_planck_create(cmd->htrdr,
      cmd->spectral_domain.wlen_range, nintervals, args->source.temperature,
      &cmd->planck);
    if(res != RES_OK) goto error;
  }

exit:
  return res;
error:
  goto exit;
}

static INLINE res_T
setup_spectral_domain
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args)
{
  double ground_T_range[2];
  size_t nintervals;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  cmd->spectral_domain = args->spectral_domain;

  /* Configure the spectral distribution */
  switch(cmd->spectral_domain.type) {

    case HTRDR_SPECTRAL_LW:
      res = rngrd_get_temperature_range(cmd->ground, ground_T_range);
      if(res != RES_OK) goto error;

      /* Use as the reference temperature of the Planck distribution the
       * maximum scene temperature which, in fact, should be the maximum ground
       * temperature */
      nintervals = compute_nintervals_for_spectral_cdf(cmd);
      res = htrdr_ran_wlen_planck_create(cmd->htrdr,
        cmd->spectral_domain.wlen_range, nintervals, ground_T_range[1],
        &cmd->planck);
      if(res != RES_OK) goto error;
      break;

    case HTRDR_SPECTRAL_SW:
      res = setup_spectral_domain_sw(cmd, args);
      if(res != RES_OK) goto error;
      break;

    case HTRDR_SPECTRAL_SW_CIE_XYZ:
      /* CIE XYZ distribution */
      nintervals = compute_nintervals_for_spectral_cdf(cmd);
      res = htrdr_ran_wlen_cie_xyz_create(cmd->htrdr,
        cmd->spectral_domain.wlen_range, nintervals, &cmd->cie);
      if(res != RES_OK) goto error;
      break;

    default: FATAL("Unreachable code\n"); break;
  }

exit:
  return res;
error:
  goto exit;
}

static res_T
setup_output
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args)
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
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args)
{
  res_T res = RES_OK;
  ASSERT(cmd && args);

  if(cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_OCTREES)
    goto exit;

  res = htrdr_planets_source_create(cmd->htrdr, &args->source, &cmd->source);
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  goto exit;
}

static res_T
setup_camera
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args)
{
  struct scam_perspective_args cam_args = SCAM_PERSPECTIVE_ARGS_DEFAULT;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  if(cmd->output_type != HTRDR_PLANETS_ARGS_OUTPUT_IMAGE)
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
setup_buffer_image
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args)
{
  struct htrdr_pixel_format pixfmt = HTRDR_PIXEL_FORMAT_NULL;
  res_T res = RES_OK;
  ASSERT(cmd && args);
  ASSERT(cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_IMAGE);

  planets_get_pixel_format(cmd, &pixfmt);

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

static res_T
setup_buffer_raw
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args)
{
  struct smsh_desc desc = SMSH_DESC_NULL;
  size_t sz = 0; /* Size of a voxel storing volumic radiative budget */
  size_t al = 0; /* Alignment of a voxel storing volumic radiative budget */
  res_T res = RES_OK;

  ASSERT(cmd && args);
  ASSERT(cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_VOLUMIC_RADIATIVE_BUDGET);
  ASSERT(cmd->volrad_mesh != NULL); /* The volurad mesh must be defined */

  res = smsh_get_desc(cmd->volrad_mesh, &desc);
  if(res != RES_OK) goto error;

  /* Setup buffer layout for volumic radiative budget calculation */
  sz = sizeof(struct planets_voxel_radiative_budget);
  al = ALIGNOF(struct planets_voxel_radiative_budget);
  cmd->buf_layout.width = desc.ncells;
  cmd->buf_layout.height = 1;
  cmd->buf_layout.pitch = desc.ncells * sz;
  cmd->buf_layout.elmt_size = sz;
  cmd->buf_layout.alignment = al;

  /* Save the number of samples per tetrahedron */
  cmd->spt = args->volrad_budget.spt;

  /* Create the raw buffer only on master process; buffer parts calculated by
   * other processes are collected there */
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
setup_buffer
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args)
{
  res_T res = RES_OK;
  ASSERT(cmd && args);

  switch(cmd->output_type) {
    case HTRDR_PLANETS_ARGS_OUTPUT_IMAGE:
      res = setup_buffer_image(cmd, args);
      break;
    case HTRDR_PLANETS_ARGS_OUTPUT_VOLUMIC_RADIATIVE_BUDGET:
      res = setup_buffer_raw(cmd, args);
      break;
    default: /* Nothing to do */ break;
  }
  if(res != RES_OK) goto error;

exit:
  return res;
error:
  goto exit;
}

static res_T
setup_volrad_budget_mesh
  (struct htrdr_planets* cmd,
   const struct htrdr_planets_args* args)
{
  struct smsh_create_args create_args = SMSH_CREATE_ARGS_DEFAULT;
  struct smsh_load_args load_args = SMSH_LOAD_ARGS_NULL;
  struct smsh_desc desc = SMSH_DESC_NULL;
  res_T res = RES_OK;
  ASSERT(cmd && args);

  if(cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_VOLUMIC_RADIATIVE_BUDGET)
    goto exit;

  /* Store the number of samples per tetrahedron to be used */
  cmd->spt = args->volrad_budget.spt;

  create_args.logger = htrdr_get_logger(cmd->htrdr);
  create_args.allocator = htrdr_get_allocator(cmd->htrdr);
  create_args.verbose = htrdr_get_verbosity_level(cmd->htrdr);
  res = smsh_create(&create_args, &cmd->volrad_mesh);
  if(res != RES_OK) goto error;

  load_args.path = args->volrad_budget.smsh_filename;
  res = smsh_load(cmd->volrad_mesh, &load_args);
  if(res != RES_OK) goto error;

  res = smsh_get_desc(cmd->volrad_mesh, &desc);
  if(res != RES_OK) goto error;

  /* Check that the loaded mesh is effectively a volume mesh */
  if(desc.dnode != 3 || desc.dcell != 4) {
    htrdr_log_err(cmd->htrdr,
      "%s: the volumic radiative budget calculation "
      "expects a 3D tetrahedral mesh "
      "(dimension of mesh: %u; dimension of the vertices: %u)\n",
      args->volrad_budget.smsh_filename,
      desc.dnode,
      desc.dcell);
    res = RES_BAD_ARG;
    goto error;
  }

exit:
  return res;
error:
  if(cmd->volrad_mesh) {
    SMSH(ref_put(cmd->volrad_mesh));
    cmd->volrad_mesh = NULL;
  }
  goto exit;
}

static INLINE res_T
write_vtk_octrees(const struct htrdr_planets* cmd)
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
planets_release(ref_T* ref)
{
  struct htrdr_planets* cmd = CONTAINER_OF(ref, struct htrdr_planets, ref);
  struct htrdr* htrdr = NULL;
  ASSERT(ref);

  if(cmd->atmosphere) RNATM(ref_put(cmd->atmosphere));
  if(cmd->ground) RNGRD(ref_put(cmd->ground));
  if(cmd->source) htrdr_planets_source_ref_put(cmd->source);
  if(cmd->cie) htrdr_ran_wlen_cie_xyz_ref_put(cmd->cie);
  if(cmd->discrete) htrdr_ran_wlen_discrete_ref_put(cmd->discrete);
  if(cmd->planck) htrdr_ran_wlen_planck_ref_put(cmd->planck);
  if(cmd->octrees_storage) CHK(fclose(cmd->octrees_storage) == 0);
  if(cmd->output && cmd->output != stdout) CHK(fclose(cmd->output) == 0);
  if(cmd->buf) htrdr_buffer_ref_put(cmd->buf);
  if(cmd->camera) SCAM(ref_put(cmd->camera));
  if(cmd->volrad_mesh) SMSH(ref_put(cmd->volrad_mesh));
  str_release(&cmd->output_name);

  htrdr = cmd->htrdr;
  MEM_RM(htrdr_get_allocator(htrdr), cmd);
  htrdr_ref_put(htrdr);
}

/*******************************************************************************
 * Exported functions
 ******************************************************************************/
res_T
htrdr_planets_create
  (struct htrdr* htrdr,
   const struct htrdr_planets_args* args,
   struct htrdr_planets** out_cmd)
{
  struct htrdr_planets* cmd = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && out_cmd);

  res = htrdr_planets_args_check(args);
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "Invalid htrdr_planets arguments -- %s\n",
      res_to_cstr(res));
    goto error;
  }

  cmd = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*cmd));
  if(!cmd) {
    htrdr_log_err(htrdr, "Error allocating htrdr_planets command\n");
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
  res = setup_volrad_budget_mesh(cmd, args);
  if(res != RES_OK) goto error;
  res = setup_buffer(cmd, args);
  if(res != RES_OK) goto error;

exit:
  *out_cmd = cmd;
  return res;
error:
  if(cmd) {
    htrdr_planets_ref_put(cmd);
    cmd = NULL;
  }
  goto exit;
}

void
htrdr_planets_ref_get(struct htrdr_planets* cmd)
{
  ASSERT(cmd);
  ref_get(&cmd->ref);
}

void
htrdr_planets_ref_put(struct htrdr_planets* cmd)
{
  ASSERT(cmd);
  ref_put(&cmd->ref, planets_release);
}

res_T
htrdr_planets_run(struct htrdr_planets* cmd)
{
  res_T res = RES_OK;
  ASSERT(cmd);

  switch(cmd->output_type) {
    case HTRDR_PLANETS_ARGS_OUTPUT_IMAGE:
      res = planets_draw_map(cmd);
      break;
    case HTRDR_PLANETS_ARGS_OUTPUT_OCTREES:
      res = write_vtk_octrees(cmd);
      break;
    case HTRDR_PLANETS_ARGS_OUTPUT_VOLUMIC_RADIATIVE_BUDGET:
      htrdr_log_err(cmd->htrdr,
        "The calculation of volumic radiative budget "
        "has not yet been implemented.\n");
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
planets_get_pixel_format
  (const struct htrdr_planets* cmd,
   struct htrdr_pixel_format* fmt)
{
  ASSERT(cmd && fmt && cmd->output_type == HTRDR_PLANETS_ARGS_OUTPUT_IMAGE);
  (void)cmd;

  switch(cmd->spectral_domain.type) {
    case HTRDR_SPECTRAL_LW:
    case HTRDR_SPECTRAL_SW:
      fmt->size = sizeof(struct planets_pixel_xwave);
      fmt->alignment = ALIGNOF(struct planets_pixel_xwave);
      break;
    case HTRDR_SPECTRAL_SW_CIE_XYZ:
      fmt->size = sizeof(struct planets_pixel_image);
      fmt->alignment = ALIGNOF(struct planets_pixel_image);
      break;
    default: FATAL("Unreachable code\n"); break;
  }
}
