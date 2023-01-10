/* Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
 * Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
 * Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
 * Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
 * Copyright (C) 2022-2023 Université de Reims Champagne-Ardenne
 * Copyright (C) 2022-2023 Université de Versaille Saint-Quentin
 * Copyright (C) 2018-2019, 2022-2023 Université Paul Sabatier
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

#define _POSIX_C_SOURCE 200112L /* strtok_r and wordexp support */

#include "core/htrdr.h"
#include "core/htrdr_log.h"
#include "core/htrdr_materials.h"

#include <modradurb/mrumtl.h>
#include <star/ssf.h>
#include <star/ssp.h>

#include <rsys/cstr.h>
#include <rsys/double3.h>
#include <rsys/hash_table.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>
#include <rsys/str.h>
#include <rsys/text_reader.h>

#include <string.h>
#include <wordexp.h>

struct mtl {
  struct mrumtl* mrumtl;
  double temperature; /* In Kelvin */
};
static const struct mtl MTL_NULL = {NULL, -1};

/* Generate the hash table that maps a material name to its data */
#define HTABLE_NAME name2mtl
#define HTABLE_DATA struct mtl
#define HTABLE_KEY struct str
#define HTABLE_KEY_FUNCTOR_INIT str_init
#define HTABLE_KEY_FUNCTOR_RELEASE str_release
#define HTABLE_KEY_FUNCTOR_COPY str_copy
#define HTABLE_KEY_FUNCTOR_COPY_AND_RELEASE str_copy_and_release
#define HTABLE_KEY_FUNCTOR_HASH str_hash
#define HTABLE_KEY_FUNCTOR_EQ str_eq
#include <rsys/hash_table.h>

struct htrdr_materials {
  struct htable_name2mtl name2mtl;
  struct htrdr* htrdr;
  ref_T ref;
};

/*******************************************************************************
 * Local functions
 ******************************************************************************/
static res_T
parse_material
  (struct htrdr_materials* mats,
   struct txtrdr* txtrdr,
   struct str* str) /* Scratch string */
{
  struct mrumtl_create_args mrumtl_args = MRUMTL_CREATE_ARGS_DEFAULT;
  wordexp_t wexp;
  char* tk = NULL;
  char* tk_ctx = NULL;
  struct mtl mtl = MTL_NULL;
  int err = 0;
  int wexp_is_allocated = 0;
  res_T res = RES_OK;
  ASSERT(mats && txtrdr);

  tk = strtok_r(txtrdr_get_line(txtrdr), " \t", &tk_ctx);
  ASSERT(tk);

  res = str_set(str, tk);
  if(res != RES_OK) {
    htrdr_log_err(mats->htrdr,
      "%s:%lu: could not copy the material name `%s' -- %s.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr), tk,
      res_to_cstr(res));
    goto error;
  }

  tk = strtok_r(NULL, "", &tk_ctx);
  if(!tk) {
    htrdr_log_err(mats->htrdr,
      "%s:%lu: missing the MruMtl file for the material `%s'.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
      str_cget(str));
    res = RES_BAD_ARG;
    goto error;
  }

  err = wordexp(tk, &wexp, 0);
  if(err) {
    htrdr_log_err(mats->htrdr,
      "%s:%lu: error in word expension of the mrumtl path.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr));
    res = RES_BAD_ARG;
    goto error;
  }
  wexp_is_allocated = 1;

  if(wexp.we_wordc < 1) {
    htrdr_log_err(mats->htrdr,
      "%s:%lu: missing the MruMtl file for the material `%s'.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
      str_cget(str));
    res = RES_BAD_ARG;
    goto error;
  }

  /*  Parse the mrumtl file if any */
  if(strcmp(wexp.we_wordv[0], "none")) {
    mrumtl_args.logger = htrdr_get_logger(mats->htrdr);
    mrumtl_args.allocator = htrdr_get_allocator(mats->htrdr);
    mrumtl_args.verbose = htrdr_get_verbosity_level(mats->htrdr);
    res = mrumtl_create(&mrumtl_args, &mtl.mrumtl);
    if(res != RES_OK) {
      htrdr_log_err(mats->htrdr,
        "%s:%lu: error creating the MruMtl loader for the material `%s'-- %s.\n",
        txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
        str_cget(str), res_to_cstr(res));
      goto error;
    }

    res = mrumtl_load(mtl.mrumtl, wexp.we_wordv[0]);
    if(res != RES_OK) goto error;
  }

  if(wexp.we_wordc < 2) {
    if(mtl.mrumtl) {
      htrdr_log_err(mats->htrdr,
        "%s:%lu: missing temperature for the material `%s'.\n",
        txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
        str_cget(str));
      res = RES_BAD_ARG;
      goto error;
    }
  } else {
    /* Parse the temperature */
    res = cstr_to_double(wexp.we_wordv[1], &mtl.temperature);
    if(res != RES_OK) {
      htrdr_log_err(mats->htrdr,
        "%s:%lu: error parsing the temperature `%s' for the material `%s' "
        "-- %s.\n",
        txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
        wexp.we_wordv[1], str_cget(str), res_to_cstr(res));
      goto error;
    }
  }

  /* Register the material */
  res = htable_name2mtl_set(&mats->name2mtl, str, &mtl);
  if(res != RES_OK) {
    htrdr_log_err(mats->htrdr,
      "%s:%lu: could not register the material `%s' -- %s.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
      str_cget(str), res_to_cstr(res));
    goto error;
  }

  if(wexp.we_wordc > 2) {
    htrdr_log_warn(mats->htrdr, "%s:%lu: unexpected text `%s'.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
      wexp.we_wordv[2]);
  }

exit:
  if(wexp_is_allocated) wordfree(&wexp);
  return res;
error:
  if(mtl.mrumtl) MRUMTL(ref_put(mtl.mrumtl));
  goto exit;
}

static res_T
parse_materials_list
  (struct htrdr_materials* mats,
   const char* filename,
   const char* func_name)
{
  struct txtrdr* txtrdr = NULL;
  struct str str;
  res_T res = RES_OK;
  ASSERT(mats && filename && func_name);

  str_init(htrdr_get_allocator(mats->htrdr), &str);

  res = txtrdr_file(htrdr_get_allocator(mats->htrdr), filename, '#', &txtrdr);
  if(res != RES_OK) {
    htrdr_log_err(mats->htrdr,
      "%s: could not create the text reader for the material file `%s' -- %s.\n",
      func_name, filename, res_to_cstr(res));
    goto error;
  }

  for(;;) {
    res = txtrdr_read_line(txtrdr);
    if(res != RES_OK) {
      htrdr_log_err(mats->htrdr,
        "%s: error reading a line in the material file `%s' -- %s.\n",
        func_name, filename, res_to_cstr(res));
      goto error;
    }

    if(!txtrdr_get_cline(txtrdr)) break;

    res = parse_material(mats, txtrdr, &str);
    if(res != RES_OK) goto error;
  }

exit:
  str_release(&str);
  if(txtrdr) txtrdr_ref_put(txtrdr);
  return res;
error:
  goto exit;
}

static void
materials_release(ref_T* ref)
{
  struct htable_name2mtl_iterator it, it_end;
  struct htrdr_materials* mats;
  struct htrdr* htrdr;
  ASSERT(ref);
  mats = CONTAINER_OF(ref, struct htrdr_materials, ref);

  htable_name2mtl_begin(&mats->name2mtl, &it);
  htable_name2mtl_end(&mats->name2mtl, &it_end);
  while(!htable_name2mtl_iterator_eq(&it, &it_end)) {
    struct mtl* mtl = htable_name2mtl_iterator_data_get(&it);
    /* The mrumtl can be NULL for semi transparent materials */
    if(mtl->mrumtl) MRUMTL(ref_put(mtl->mrumtl));
    htable_name2mtl_iterator_next(&it);
  }
  htable_name2mtl_release(&mats->name2mtl);
  htrdr = mats->htrdr;
  MEM_RM(htrdr_get_allocator(htrdr), mats);
  htrdr_ref_put(htrdr);
}

static res_T
create_bsdf_diffuse
  (struct htrdr* htrdr,
   const struct mrumtl_brdf* brdf,
   const size_t ithread,
   struct ssf_bsdf** out_bsdf)
{
  struct mrumtl_brdf_lambertian lambert;
  struct ssf_bsdf* bsdf = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && brdf && out_bsdf);
  ASSERT(mrumtl_brdf_get_type(brdf) == MRUMTL_BRDF_LAMBERTIAN);

  res = ssf_bsdf_create(htrdr_get_thread_allocator(htrdr, ithread),
    &ssf_lambertian_reflection, &bsdf);
  if(res != RES_OK) goto error;

  MRUMTL(brdf_get_lambertian(brdf, &lambert));
  res = ssf_lambertian_reflection_setup(bsdf, lambert.reflectivity);
  if(res != RES_OK) goto error;

exit:
  *out_bsdf = bsdf;
  return res;
error:
   if(bsdf) { SSF(bsdf_ref_put(bsdf)); bsdf = NULL; }
  goto exit;
}

static res_T
create_bsdf_specular
  (struct htrdr* htrdr,
   const struct mrumtl_brdf* brdf,
   const size_t ithread,
   struct ssf_bsdf** out_bsdf)
{
  struct mrumtl_brdf_specular spec;
  struct ssf_bsdf* bsdf = NULL;
  struct ssf_fresnel* fresnel = NULL;
  struct mem_allocator* allocator = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && brdf && out_bsdf);
  ASSERT(mrumtl_brdf_get_type(brdf) == MRUMTL_BRDF_SPECULAR);

  allocator = htrdr_get_thread_allocator(htrdr, ithread);

  res = ssf_bsdf_create(allocator, &ssf_specular_reflection, &bsdf);
  if(res != RES_OK) goto error;

  res = ssf_fresnel_create(allocator, &ssf_fresnel_constant, &fresnel);
  if(res != RES_OK) goto error;

  MRUMTL(brdf_get_specular(brdf, &spec));
  res =  ssf_fresnel_constant_setup(fresnel, spec.reflectivity);
  if(res != RES_OK) goto error;

  res = ssf_specular_reflection_setup(bsdf, fresnel);
  if(res != RES_OK) goto error;

exit:
  if(fresnel) SSF(fresnel_ref_put(fresnel));
  *out_bsdf = bsdf;
  return res;
error:
  if(bsdf) { SSF(bsdf_ref_put(bsdf)); bsdf = NULL; }
  goto exit;
}

/*******************************************************************************
 * Local symbol
 ******************************************************************************/
res_T
htrdr_materials_create
  (struct htrdr* htrdr,
   const char* filename,
   struct htrdr_materials** out_mtl)
{
  struct htrdr_materials* mats = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && filename && out_mtl);

  mats = MEM_CALLOC(htrdr_get_allocator(htrdr), 1, sizeof(*mats));
  if(!mats) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "%s: could not allocate the mats data structure -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }
  ref_init(&mats->ref);
  htrdr_ref_get(htrdr);
  mats->htrdr = htrdr;
  htable_name2mtl_init(htrdr_get_allocator(htrdr), &mats->name2mtl);

  res = parse_materials_list(mats, filename, FUNC_NAME);
  if(res != RES_OK) goto error;

exit:
  if(out_mtl) *out_mtl = mats;
  return res;
error:
  if(mats) {
    htrdr_materials_ref_put(mats);
    mats = NULL;
  }
  goto exit;
}

void
htrdr_materials_ref_get(struct htrdr_materials* mats)
{
  ASSERT(mats);
  ref_get(&mats->ref);
}

void
htrdr_materials_ref_put(struct htrdr_materials* mats)
{
  ASSERT(mats);
  ref_put(&mats->ref, materials_release);
}

int
htrdr_materials_find_mtl
  (struct htrdr_materials* mats,
   const char* name,
   struct htrdr_mtl* htrdr_mtl)
{
  struct str str;
  struct htable_name2mtl_iterator it, it_end;
  int found = 0;
  ASSERT(mats && name && htrdr_mtl);

  str_init(htrdr_get_allocator(mats->htrdr), &str);
  CHK(str_set(&str, name) == RES_OK);

  htable_name2mtl_find_iterator(&mats->name2mtl, &str, &it);
  htable_name2mtl_end(&mats->name2mtl, &it_end);
  if(htable_name2mtl_iterator_eq(&it, &it_end)) { /* No material found */
    *htrdr_mtl = HTRDR_MTL_NULL;
    found = 0;
  } else {
    struct mtl* mtl = htable_name2mtl_iterator_data_get(&it);
    ASSERT(mtl != NULL);
    htrdr_mtl->name = str_cget(htable_name2mtl_iterator_key_get(&it));
    htrdr_mtl->mrumtl = mtl->mrumtl;
    htrdr_mtl->temperature = mtl->temperature;
    found = 1;
  }
  str_release(&str);

  return found;
}

res_T
htrdr_mtl_create_bsdf
  (struct htrdr* htrdr,
   const struct htrdr_mtl* mtl,
   const size_t ithread,
   const double wavelength,
   struct ssp_rng* rng,
   struct ssf_bsdf** out_bsdf)
{
  struct ssf_bsdf* bsdf = NULL;
  const struct mrumtl_brdf* brdf = NULL;
  size_t ibrdf;
  double r;
  res_T res = RES_OK;
  ASSERT(htrdr && mtl && wavelength && rng && out_bsdf);

  r = ssp_rng_canonical(rng);

  res = mrumtl_fetch_brdf(mtl->mrumtl, wavelength, r, &ibrdf);
  if(res != RES_OK) {
    htrdr_log_err(htrdr,
      "%s: error retrieving the MruMtl BRDF for the wavelength %g.\n",
      FUNC_NAME, wavelength);
    res = RES_BAD_ARG;
    goto error;
  }

  brdf = mrumtl_get_brdf(mtl->mrumtl, ibrdf);
  switch(mrumtl_brdf_get_type(brdf)) {
    case MRUMTL_BRDF_LAMBERTIAN:
      res = create_bsdf_diffuse(htrdr, brdf, ithread, &bsdf);
      break;
    case MRUMTL_BRDF_SPECULAR:
      res = create_bsdf_specular(htrdr, brdf, ithread, &bsdf);
      break;
    default: FATAL("Unreachable code.\n");  break;
  }
  if(res != RES_OK) {
    htrdr_log_err(htrdr, "%s: could not create the BSDF -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }

exit:
  *out_bsdf = bsdf;
  return res;
error:
  if(bsdf) { SSF(bsdf_ref_put(bsdf)); bsdf = NULL; }
  goto exit;
}

