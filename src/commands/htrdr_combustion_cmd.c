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

#ifdef HTRDR_BUILD_COMBUSTION
  #include "combustion/htrdr_combustion.h"
#else
  #include <stdio.h>
#endif

int
main(int argc, char** argv)
{
#ifdef HTRDR_BUILD_COMBUSTION
  return htrdr_combustion_main(argc, argv);
#else
  (void)argc, (void)argv;
  fprintf(stderr,
    "The htrdr-combustion command is not available in this htrdr build.\n");
  return 1;
#endif
}
