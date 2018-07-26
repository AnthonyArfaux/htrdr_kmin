/* Copyright (C) 2018 Université Paul Sabatier, |Meso|Star>
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

#ifndef HTRDR_OPENEXR_H
#define HTRDR_OPENEXR_H

#include <rsys/rsys.h>
#include <stdio.h>

BEGIN_DECLS

extern res_T
htrdr_openexr_write
  (const float* pixels, 
   const size_t width,
   const size_t height,
   const size_t pitch,
   const char* output_name,
   FILE* output);

END_DECLS

#endif /* HTRDR_OPENEXR_H */

