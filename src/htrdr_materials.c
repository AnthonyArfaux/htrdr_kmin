/* Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
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

#include "htrdr.h"
#include "htrdr_materials.h"

#include <modradurb/mrumtl.h>

#include <rsys/cstr.h>
#include <rsys/hash_table.h>
#include <rsys/mem_allocator.h>
#include <rsys/ref_count.h>
#include <rsys/str.h>
#include <rsys/text_reader.h>

#include <string.h>
#include <wordexp.h>

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
  wordexp_t wexp;
  char* tk = NULL;
  char* tk_ctx = NULL;
  struct mrumtl* mrumtl = NULL;
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
    res = mrumtl_create(&mats->htrdr->logger, mats->htrdr->allocator,
      mats->htrdr->verbose, &mrumtl);
    if(res != RES_OK) {
      htrdr_log_err(mats->htrdr,
        "%s:%lu: error creating the MruMtl loader for the material `%s'-- %s.\n",
        txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
        str_cget(str), res_to_cstr(res));
      goto error;
    }

    res = mrumtl_load(mrumtl, wexp.we_wordv[0]);
    if(res != RES_OK) goto error;
  }

  /* Register the material */
  res = htable_name2mtl_set(&mats->name2mtl, str, &mrumtl);
  if(res != RES_OK) {
    htrdr_log_err(mats->htrdr,
      "%s:%lu: could not register the material `%s' -- %s.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
      str_cget(str), res_to_cstr(res));
    goto error;
  }

  if(wexp.we_wordc > 1) {
    htrdr_log_warn(mats->htrdr, "%s:%lu: unexpected text `%s'.\n",
      txtrdr_get_name(txtrdr), (unsigned long)txtrdr_get_line_num(txtrdr),
      wexp.we_wordv[1]);
  }

exit:
  if(wexp_is_allocated) wordfree(&wexp);
  return res;
error:
  if(mrumtl) MRUMTL(ref_put(mrumtl));
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

  str_init(mats->htrdr->allocator, &str);

  res = txtrdr_file(mats->htrdr->allocator, filename, '#', &txtrdr);
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
mtl_release(ref_T* ref)
{
  struct htable_name2mtl_iterator it, it_end;
  struct htrdr_materials* mats;
  ASSERT(ref);
  mats = CONTAINER_OF(ref, struct htrdr_materials, ref);

  htable_name2mtl_begin(&mats->name2mtl, &it);
  htable_name2mtl_end(&mats->name2mtl, &it_end);
  while(!htable_name2mtl_iterator_eq(&it, &it_end)) {
    struct mrumtl* mrumtl = *htable_name2mtl_iterator_data_get(&it);
    /* The mrumtl can be NULL for semi transparent materials */
    if(mrumtl) MRUMTL(ref_put(mrumtl));
    htable_name2mtl_iterator_next(&it);
  }
  htable_name2mtl_release(&mats->name2mtl);
  MEM_RM(mats->htrdr->allocator, mats);
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

  mats = MEM_CALLOC(htrdr->allocator, 1, sizeof(*mats));
  if(!mats) {
    res = RES_MEM_ERR;
    htrdr_log_err(htrdr,
      "%s: could not allocate the mats data structure -- %s.\n",
      FUNC_NAME, res_to_cstr(res));
    goto error;
  }
  ref_init(&mats->ref);
  mats->htrdr = htrdr;
  htable_name2mtl_init(htrdr->allocator, &mats->name2mtl);

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
  ref_put(&mats->ref, mtl_release);
}

int
htrdr_materials_find_mtl
  (struct htrdr_materials* mats,
   const char* name,
   struct htrdr_mtl* mtl)
{
  struct str str;
  struct htable_name2mtl_iterator it, it_end;
  int found = 0;
  ASSERT(mats && name && mtl);

  str_init(mats->htrdr->allocator, &str);
  CHK(str_set(&str, name) == RES_OK);

  htable_name2mtl_find_iterator(&mats->name2mtl, &str, &it);
  htable_name2mtl_end(&mats->name2mtl, &it_end);
  if(htable_name2mtl_iterator_eq(&it, &it_end)) { /* No material found */
    *mtl = HTRDR_MTL_NULL;
    found = 0;
  } else {
    mtl->name = str_cget(htable_name2mtl_iterator_key_get(&it));
    mtl->mrumtl = *htable_name2mtl_iterator_data_get(&it);
    found = 1;
  }
  str_release(&str);

  return found;
}

