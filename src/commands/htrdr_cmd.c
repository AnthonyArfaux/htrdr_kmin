/* Copyright (C) 2018, 2019, 2020, 2021 |Meso|Star> (contact@meso-star.com)
 * Copyright (C) 2018, 2019, 2021 CNRS
 * Copyright (C) 2018, 2019 Université Paul Sabatier
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

#include "atmosphere/htrdr_atmosphere.h"
#include "combustion/htrdr_combustion.h"
#include "core/htrdr_version.h"

#include <string.h>

/*******************************************************************************
 * Helper functions
 ******************************************************************************/
static void
print_usage(const char* cmd)
{
  ASSERT(cmd);
  printf("Usage: %s [--version] [--help] <mode> [<args>]\n", cmd);
}

static void
print_help(const char* cmd)
{
  ASSERT(cmd);

  print_usage(cmd);
  printf("\n");

  printf(
"  --version      display version information and exit.\n");
  printf(
"  --help         display this help and exit.\n");
  printf("\n");

  printf("These are %s available modes:\n", cmd);
  printf("\n");
  printf(
"  atmosphere     Radiative transfer computations in a cloudy atmosphere.\n");
  printf(
"  combustion     Radiative transfer computations in a combustion medium.\n");
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
    err = htrdr_atmosphere_main(argc-1, argv+1);
    if(err) goto error;

  /* Combustion mode */
  } else if(!strcmp(argv[1], "combustion")) {
    err = htrdr_combustion_main(argc-1, argv+1);
    if(err) goto error;

  /* Version */
  } else if(!strcmp(argv[1], "--version")) {
    printf("%s version %d.%d.%d\n",
      argv[0],
      HTRDR_VERSION_MAJOR,
      HTRDR_VERSION_MINOR,
      HTRDR_VERSION_PATCH);
    goto exit;

  /* Help */
  } else if(!strcmp(argv[1], "--help")) {
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

