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

#include "planeto/htrdr_planeto_args.h"
#include "planeto/htrdr_planeto_source.h"

#include "core/htrdr.h"

#include <stdio.h>

static void
write_per_wlen_radiances
  (FILE* fp,
   const size_t pagesize,
   const size_t size,
   const size_t szelmt,
   const size_t alelmt)
{
  const char byte = 0;
  size_t i;

  CHK(fp);

  /* Header */
  CHK(fwrite(&pagesize, sizeof(pagesize), 1, fp) == 1);
  CHK(fwrite(&size, sizeof(size), 1, fp) == 1);
  CHK(fwrite(&szelmt, sizeof(szelmt), 1, fp) == 1);
  CHK(fwrite(&alelmt, sizeof(alelmt), 1, fp) == 1);

  /* Padding */
  CHK(fseek(fp, (long)ALIGN_SIZE((size_t)ftell(fp), pagesize), SEEK_SET) == 0);

  FOR_EACH(i, 0, size) {
    const double w = (double)i;
    const double L = (double)(100 + i);

    CHK(fwrite(&w, sizeof(w), 1, fp) == 1);
    CHK(fwrite(&L, sizeof(L), 1, fp) == 1);
  }

  /* Padding. Write one char to position the EOF indicator */
  CHK(fseek(fp, (long)ALIGN_SIZE((size_t)ftell(fp), pagesize)-1, SEEK_SET) == 0);
  CHK(fwrite(&byte, sizeof(byte), 1, fp) == 1);

  CHK(fflush(fp) == 0);
}

static void
test_spectrum(struct htrdr* htrdr)
{
  struct htrdr_planeto_source_args source_args = HTRDR_PLANETO_SOURCE_ARGS_NULL;
  struct htrdr_planeto_source_spectrum spectrum;
  struct htrdr_planeto_source* source = NULL;

  FILE* fp = NULL;
  char rnrl_filename[] = "rnrl.bin";
  double range[2];
  double w, L;

  CHK(fp = fopen(rnrl_filename, "w"));
  write_per_wlen_radiances(fp, 4096, 10, 16, 16);
  CHK(fclose(fp) == 0);

  source_args.rnrl_filename = rnrl_filename;
  source_args.longitude = 0;
  source_args.latitude = 0;
  source_args.distance = 0;
  source_args.radius = 1e8;
  source_args.temperature = -1;
  CHK(htrdr_planeto_source_create(htrdr, &source_args, &source) == RES_OK);
  CHK(htrdr_planeto_source_does_radiance_vary_spectrally(source) == 1);
  CHK(htrdr_planeto_source_get_spectral_range(source, range) == RES_OK);
  CHK(range[0] == 0);
  CHK(range[1] == 9);

  range[0] = 0; range[1] = 10;
  CHK(htrdr_planeto_source_get_spectrum(source, range, &spectrum) == RES_BAD_ARG);

  range[0] = 1; range[1] = 3;
  CHK(htrdr_planeto_source_get_spectrum(source, range, &spectrum) == RES_OK);
  CHK(spectrum.source == source);
  CHK(spectrum.range[0] == 1);
  CHK(spectrum.range[1] == 3);
  CHK(spectrum.size == 3);

  htrdr_planeto_source_spectrum_at(&spectrum, 0, &w, &L);
  CHK(w == 1 && L == 101);
  htrdr_planeto_source_spectrum_at(&spectrum, 1, &w, &L);
  CHK(w == 2 && L == 102);
  htrdr_planeto_source_spectrum_at(&spectrum, 2, &w, &L);
  CHK(w == 3 && L == 103);

  range[0] = 1.7; range[1] = 1.95;
  CHK(htrdr_planeto_source_get_spectrum(source, range, &spectrum) == RES_OK);
  CHK(spectrum.source == source);
  CHK(spectrum.range[0] = 1.7);
  CHK(spectrum.range[1] = 1.95);
  CHK(spectrum.size == 2);
  htrdr_planeto_source_spectrum_at(&spectrum, 0, &w, &L);
  CHK(w == 1.7 && eq_eps(L, 101.7, 1.e-6));
  htrdr_planeto_source_spectrum_at(&spectrum, 1, &w, &L);
  CHK(w == 1.95 && eq_eps(L, 101.95, 1.e-6));

  range[0] = 2; range[1] = 2.01;
  CHK(htrdr_planeto_source_get_spectrum(source, range, &spectrum) == RES_OK);
  CHK(spectrum.size == 2);
  htrdr_planeto_source_spectrum_at(&spectrum, 0, &w, &L);
  CHK(w == 2 && L == 102);
  htrdr_planeto_source_spectrum_at(&spectrum, 1, &w, &L);
  CHK(w == 2.01 && eq_eps(L, 102.01, 1.e-6));

  range[0] = 5.1; range[1] = 6;
  CHK(htrdr_planeto_source_get_spectrum(source, range, &spectrum) == RES_OK);
  CHK(spectrum.size == 2);
  htrdr_planeto_source_spectrum_at(&spectrum, 0, &w, &L);
  CHK(w == 5.1 && eq_eps(L, 105.1, 1.e-6));
  htrdr_planeto_source_spectrum_at(&spectrum, 1, &w, &L);
  CHK(w == 6 && L == 106);

  range[0] = 7.5; range[1] = 9;
  CHK(htrdr_planeto_source_get_spectrum(source, range, &spectrum) == RES_OK);
  CHK(spectrum.size == 3);
  htrdr_planeto_source_spectrum_at(&spectrum, 0, &w, &L);
  CHK(w == 7.5 && eq_eps(L, 107.5, 1.e-6));
  htrdr_planeto_source_spectrum_at(&spectrum, 1, &w, &L);
  CHK(w == 8 && L == 108);
  htrdr_planeto_source_spectrum_at(&spectrum, 2, &w, &L);
  CHK(w == 9 && L == 109);

  range[0] = 0.9; range[1] = 7.456;
  CHK(htrdr_planeto_source_get_spectrum(source, range, &spectrum) == RES_OK);
  CHK(spectrum.size == 9);
  htrdr_planeto_source_spectrum_at(&spectrum, 0, &w, &L);
  CHK(w == 0.9 && eq_eps(L, 100.9, 1.e-6));
  htrdr_planeto_source_spectrum_at(&spectrum, 1, &w, &L);
  CHK(w == 1 && eq_eps(L, 101, 1.e-6));
  htrdr_planeto_source_spectrum_at(&spectrum, 2, &w, &L);
  CHK(w == 2 && eq_eps(L, 102, 1.e-6));
  htrdr_planeto_source_spectrum_at(&spectrum, 3, &w, &L);
  CHK(w == 3 && eq_eps(L, 103, 1.e-6));
  htrdr_planeto_source_spectrum_at(&spectrum, 4, &w, &L);
  CHK(w == 4 && eq_eps(L, 104, 1.e-6));
  htrdr_planeto_source_spectrum_at(&spectrum, 5, &w, &L);
  CHK(w == 5 && eq_eps(L, 105, 1.e-6));
  htrdr_planeto_source_spectrum_at(&spectrum, 6, &w, &L);
  CHK(w == 6 && eq_eps(L, 106, 1.e-6));
  htrdr_planeto_source_spectrum_at(&spectrum, 7, &w, &L);
  CHK(w == 7 && eq_eps(L, 107, 1.e-6));
  htrdr_planeto_source_spectrum_at(&spectrum, 8, &w, &L);
  CHK(w == 7.456 && eq_eps(L, 107.456, 1.e-6));

  htrdr_planeto_source_ref_put(source);
}

static void
test_spectrum_fail(struct htrdr* htrdr)
{
  struct htrdr_planeto_source_args source_args = HTRDR_PLANETO_SOURCE_ARGS_NULL;
  struct htrdr_planeto_source* source = NULL;
  FILE* fp = NULL;
  char rnrl_filename[] = "rnrl.bin";
  double w, L;

  source_args.rnrl_filename = rnrl_filename;
  source_args.longitude = 0;
  source_args.latitude = 0;
  source_args.distance = 0;
  source_args.radius = 1e8;
  source_args.temperature = -1;

  /* Wrong item size */
  CHK(fp = fopen(rnrl_filename, "w"));
  write_per_wlen_radiances(fp, 4096, 10, 8, 16);
  CHK(fclose(fp) == 0);
  CHK(htrdr_planeto_source_create(htrdr, &source_args, &source) == RES_BAD_ARG);

  /* Wrong item alignment */
  CHK(fp = fopen(rnrl_filename, "w"));
  write_per_wlen_radiances(fp, 4096, 10, 16, 32);
  CHK(fclose(fp) == 0);
  CHK(htrdr_planeto_source_create(htrdr, &source_args, &source) == RES_BAD_ARG);

  CHK(fp = fopen(rnrl_filename, "w"));
  write_per_wlen_radiances(fp, 4096, 4, 16, 16);

  /* Overwrite sorted items by unsorted items */
  CHK(fseek(fp, 4096, SEEK_SET) == 0);
  w = 10; L = 1;
  CHK(fwrite(&w, sizeof(w), 1, fp) == 1);
  CHK(fwrite(&L, sizeof(L), 1, fp) == 1);
  w = 11; L = 2;
  CHK(fwrite(&w, sizeof(w), 1, fp) == 1);
  CHK(fwrite(&L, sizeof(L), 1, fp) == 1);
  w = 9; L = 3;
  CHK(fwrite(&w, sizeof(w), 1, fp) == 1);
  CHK(fwrite(&L, sizeof(L), 1, fp) == 1);
  w = 12; L = 4;
  CHK(fwrite(&w, sizeof(w), 1, fp) == 1);
  CHK(fwrite(&L, sizeof(L), 1, fp) == 1);
  CHK(fclose(fp) == 0);

  /* Unsorted items */
  CHK(htrdr_planeto_source_create(htrdr, &source_args, &source) == RES_BAD_ARG);
}

int
main(int argc, char** argv)
{
  struct htrdr_args args = HTRDR_ARGS_DEFAULT;
  struct htrdr* htrdr = NULL;

  args.verbose = 1;
  htrdr_mpi_init(argc, argv);
  CHK(htrdr_create(NULL, &args, &htrdr) == RES_OK);

  test_spectrum(htrdr);
  test_spectrum_fail(htrdr);

  htrdr_ref_put(htrdr);
  htrdr_mpi_finalize();
  return 0;
}
