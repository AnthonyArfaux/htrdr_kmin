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

#ifdef HTRDR_BUILD_ATMOSPHERE
  #include "atmosphere/htrdr_atmosphere.h"
#endif
#ifdef HTRDR_BUILD_COMBUSTION
  #include "combustion/htrdr_combustion.h"
#endif
#ifdef HTRDR_BUILD_PLANETO
  #include "planeto/htrdr_planeto.h"
#endif

#include "core/htrdr_log.h"
#include "core/htrdr_version.h"

#include <rsys/rsys.h>
#include <string.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
print_usage(const char* cmd)
{
  ASSERT(cmd);
  printf("Usage: %s [-v] [-h] <mode> [option] ...\n", cmd);
}

static void
print_help(const char* cmd)
{
  ASSERT(cmd);

  print_usage(cmd);
  printf("Simulate radiative transfer.\n");
  printf("See htrdr(1) man page for details\n");
  printf("\n");

  printf(
"  -h             display this help and exit\n");
  printf(
"  -v             display version information and exit\n");
  printf("\n");

  printf("These are %s modes:\n", cmd);
  printf("\n");
  printf(
"  atmosphere     Radiative transfer in a plane-parallel atmosphere\n");
  printf(
"  combustion     Radiative transfer within a sooting flame\n");
  printf(
"  planeto        Radiative transfer in a 3D planetory atmosphere\n");
  printf("\n");

  htrdr_fprint_license(cmd, stdout);
}

/*******************************************************************************
 * Program
 ******************************************************************************/
int
main(int argc, char** argv)
{
  char cmd_name[] = "htrdr";
  int err = 0;

  if(argc < 2) {
    print_usage(argv[0]);
    err = -1;
    goto error;
  }

  /* Atmosphere mode */
  if(!strcmp(argv[1], "atmosphere")) {
#ifdef HTRDR_BUILD_ATMOSPHERE
    err = htrdr_atmosphere_main(argc-1, argv+1);
    if(err) goto error;
#else
    fprintf(stderr,
      "The atmosphere mode is not available in this htrdr build.\n");
    err = 1;
    goto error;
#endif

  /* Combustion mode */
  } else  if(!strcmp(argv[1], "combustion")) {
#ifdef HTRDR_BUILD_COMBUSTION
    err = htrdr_combustion_main(argc-1, argv+1);
    if(err) goto error;
#else
    fprintf(stderr,
      "The combustion mode is not available in this htrdr build.\n");
    err = 1;
    goto error;
#endif

  /* Planeto mode */
  } else if(!strcmp(argv[1], "planeto")) {
#ifdef HTRDR_BUILD_PLANETO
    err = htrdr_planeto_main(argc-1, argv+1);
    if(err) goto error;
#else
    fprintf(stderr,
      "The planeto mode is not available in this htrdr build.\n");
    err = 1;
    goto error;
#endif

  /* Version */
  } else if(!strcmp(argv[1], "-v")) {
    printf("%s version %d.%d.%d\n",
      argv[0],
      HTRDR_VERSION_MAJOR,
      HTRDR_VERSION_MINOR,
      HTRDR_VERSION_PATCH);
    goto exit;

  /* Help */
  } else if(!strcmp(argv[1], "-h")) {
    print_help(cmd_name);
    goto exit;

  /* Fallback */
  } else {
    fprintf(stderr, "Unknown option: %s\n", argv[1]);
    print_usage(cmd_name);
    err = -1;
    goto error;
  }

exit:
  return err;
error:
  goto exit;
}

