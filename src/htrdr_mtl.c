/* Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
 * Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
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

#define _POSIX_C_SOURCE 200112L /* strtok_r support */

#include "htrdr.h"
#include "htrdr_mtl.h"

#include <modradurb/mrumtl.h>

#include <rsys/cstr.h>
#include <rsys/hash_table.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>
#include <rsys/str.h>
#include <rsys/text_reader.h>

#include <string.h>

/* Generate the hash table that maps a material name to its data */
#define HTABLE_NAME name2mtl
#define HTABLE_DATA struct mrumtl*
#define HTABLE_KEY struct str
#define HTABLE_KEY_FUNCTOR_INIT str_init
#define HTABLE_KEY_FUNCTOR_RELEASE str_release
#define HTABLE_KEY_FUNCTOR_COPY str_copy
#define HTABLE_KEY_FUNCTOR_COPY_AND_RELEASE str_copy_and_release
#define HTABLE_KEY_FUNCTOR_HASH str_hash
#define HTABLE_KEY_FUNCTOR_EQ str_eq
#include <rsys/hash_table.h>

struct htrdr_mtl {
  struct htable_name2mtl name2mtl;
  struct htrdr* htrdr;
  ref_T ref;
};

/*******************************************************************************
 * Local functions
 ******************************************************************************/
static res_T
parse_material
  (struct htrdr_mtl* mtl,
   struct txtrdr* txtrdr,
   struct str* str) /* Scratch string */
{
  char* tk = NULL;
  char* tk_ctx = NULL;
  struct mrumtl* mrumtl = NULL;
  res_T res = RES_OK;
  ASSERT(mtl && txtrdr);

  tk = strtok_r(txtrdr_get_line(txtrdr), " \t", &tk_ctx);
  ASSERT(tk);

  res = str_set(str, tk);
  if(res != RES_OK) {
    htrdr_log_err(mtl->htrdr,
      "%s:%lu: could not copy the material name `%s' -- %s.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr), tk,
      res_to_cstr(res));
    goto error;
  }

  tk = strtok_r(NULL, "", &tk_ctx);
  if(!tk) {
    htrdr_log_err(mtl->htrdr,
      "%s:%lu: missing the MruMtl file for the material `%s'.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
      str_cget(str));
    res = RES_BAD_ARG;
    goto error;
  }

  res = mrumtl_create
    (&mtl->htrdr->logger, mtl->htrdr->allocator, mtl->htrdr->verbose, &mrumtl);
  if(res != RES_OK) {
    htrdr_log_err(mtl->htrdr,
      "%s:%lu: error creating the MruMtl loader for the material `%s'-- %s.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
      str_cget(str), res_to_cstr(res));
    goto error;
  }

  res = mrumtl_load(mrumtl, tk);
  if(res != RES_OK) goto error;

  /* Register the material */
  res = htable_name2mtl_set(&mtl->name2mtl, str, &mrumtl);
  if(res != RES_OK) {
    htrdr_log_err(mtl->htrdr,
      "%s:%lu: could not register the material `%s' -- %s.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
      str_cget(str), res_to_cstr(res));
    goto error;
  }

exit:
  return res;
error:
  if(mrumtl) MRUMTL(ref_put(mrumtl));
  goto exit;
}

static res_T
parse_materials_list
  (struct htrdr_mtl* mtl,
   const char* filename,
   const char* func_name)
{
  struct txtrdr* txtrdr = NULL;
  struct str str;
  res_T res = RES_OK;
  ASSERT(mtl && filename && func_name);

  str_init(mtl->htrdr->allocator, &str);

  res = txtrdr_file(mtl->htrdr->allocator, filename, '#', &txtrdr);
  if(res != RES_OK) {
    htrdr_log_err(mtl->htrdr,
      "%s: could not create the text reader for the material file `%s' -- %s.\n",
      func_name, filename, res_to_cstr(res));
    goto error;
  }

  for(;;) {
    res = txtrdr_read_line(txtrdr);
    if(res != RES_OK) {
      htrdr_log_err(mtl->htrdr,
        "%s: error reading a line in the material file `%s' -- %s.\n",
        func_name, filename, res_to_cstr(res));
      goto error;
    }

    if(!txtrdr_get_cline(txtrdr)) break;

    res = parse_material(mtl, txtrdr, &str);
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
mtl_release(ref_T* ref)
{
  struct htable_name2mtl_iterator it, it_end;
  struct htrdr_mtl* mtl;
  ASSERT(ref);
  mtl = CONTAINER_OF(ref, struct htrdr_mtl, ref);

  htable_name2mtl_begin(&mtl->name2mtl, &it);
  htable_name2mtl_end(&mtl->name2mtl, &it_end);
  while(!htable_name2mtl_iterator_eq(&it, &it_end)) {
    struct mrumtl* mrumtl = *htable_name2mtl_iterator_data_get(&it);
    MRUMTL(ref_put(mrumtl));
    htable_name2mtl_iterator_next(&it);
  }
  htable_name2mtl_release(&mtl->name2mtl);
  MEM_RM(mtl->htrdr->allocator, mtl);
}

/*******************************************************************************
 * Local symbol
 ******************************************************************************/
res_T
htrdr_mtl_create
  (struct htrdr* htrdr,
   const char* filename,
   struct htrdr_mtl** out_mtl)
{
  struct htrdr_mtl* mtl = NULL;
  res_T res = RES_OK;
  ASSERT(htrdr && filename && mtl);

  mtl = MEM_CALLOC(htrdr->allocator, 1, sizeof(*mtl));
  if(!mtl) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "%s: could not allocate the mtl data structure -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }
  ref_init(&mtl->ref);
  mtl->htrdr = htrdr;
  htable_name2mtl_init(htrdr->allocator, &mtl->name2mtl);

  res = parse_materials_list(mtl, filename, FUNC_NAME);
  if(res != RES_OK) goto error;

exit:
  if(out_mtl) *out_mtl = mtl;
  return res;
error:
  if(mtl) {
    htrdr_mtl_ref_put(mtl);
    mtl = NULL;
  }
  goto exit;
}

void
htrdr_mtl_ref_get(struct htrdr_mtl* mtl)
{
  ASSERT(mtl);
  ref_get(&mtl->ref);
}

void
htrdr_mtl_ref_put(struct htrdr_mtl* mtl)
{
  ASSERT(mtl);
  ref_put(&mtl->ref, mtl_release);
}

const struct mrumtl*
htrdr_mtl_get(struct htrdr_mtl* mtl, const char* name)
{
  struct str str;
  struct mrumtl** pmrumtl = NULL;
  struct mrumtl* mrumtl = NULL;
  ASSERT(mtl && name);

  str_init(mtl->htrdr->allocator, &str);
  CHK(str_set(&str, name) == RES_OK);

  pmrumtl = htable_name2mtl_find(&mtl->name2mtl, &str);
  if(pmrumtl) mrumtl = *pmrumtl;

  str_release(&str);
  return mrumtl;
}

