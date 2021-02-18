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

#ifndef HTRDR_ARGS_H
#define HTRDR_ARGS_H

#include "htrdr_cie_xyz.h"

#include <rsys/rsys.h>

/* Arguments of a pinhole camera sensor */
struct htrdr_args_camera {
  double position[3]; /* Focal point */
  double target[3]; /* Targeted position */
  double up[3]; /* Up vector of the camera */
  double fov_y; /* Vertical field of view of the camera, in degrees */
};
#define HTRDR_ARGS_CAMERA_DEFAULT__ {                                          \
  {@HTRDR_ARGS_DEFAULT_CAMERA_POS@}, /* position */                            \
  {@HTRDR_ARGS_DEFAULT_CAMERA_TGT@}, /* target */                              \
  {@HTRDR_ARGS_DEFAULT_CAMERA_UP@}, /* Camera up */                            \
  @HTRDR_ARGS_DEFAULT_CAMERA_FOV@, /* Vertical field of view */                \
}
static const struct htrdr_args_camera HTRDR_ARGS_CAMERA_DEFAULT =
  HTRDR_ARGS_CAMERA_DEFAULT__;

/* Arguments of a rectangular sensor */
struct htrdr_args_rectangle {
  double position[3]; /* Center of the renctangle */
  double target[3]; /* Targeted point (rectangle.normal = target - position) */
  double up[3]; /* Up vector of the rectangle */
  double size[2]; /* Plane size */
};
#define HTRDR_ARGS_RECTANGLE_DEFAULT__ {                                       \
  {@HTRDR_ARGS_DEFAULT_RECTANGLE_POS@}, /* Rectangle center */                 \
  {@HTRDR_ARGS_DEFAULT_RECTANGLE_TGT@}, /* Rectangle target */                 \
  {@HTRDR_ARGS_DEFAULT_RECTANGLE_UP@}, /* Rectangle up */                      \
  {@HTRDR_ARGS_DEFAULT_RECTANGLE_SZ@}, /* Rectangle size */                    \
}
static const struct htrdr_args_rectangle HTRDR_ARGS_RECTANGLE_DEFAULT =
  HTRDR_ARGS_RECTANGLE_DEFAULT__;

/* Arguments of an image */
struct htrdr_args_image {
  unsigned definition[2]; /* #pixels in X and Y */
  unsigned spp; /* #samples per pixel */
};
#define HTRDR_ARGS_IMAGE_DEFAULT__ {                                           \
  {@HTRDR_ARGS_DEFAULT_IMG_WIDTH@, @HTRDR_ARGS_DEFAULT_IMG_HEIGHT@},           \
  @HTRDR_ARGS_DEFAULT_IMG_SPP@                                                 \
}
static const struct htrdr_args_image HTRDR_ARGS_IMAGE_DEFAULT =
  HTRDR_ARGS_IMAGE_DEFAULT__;

/* Arguments of the spectral domain */
struct htrdr_args_spectral {
  double wlen_range[2]; /* Spectral range of integration in nm */
  double ref_temperature; /* Planck reference temperature in Kelvin */
  enum htrdr_spectral_type spectral_type;
};
#define HTRDR_ARGS_SPECTRAL_DEFAULT__ {                                        \
  HTRDR_CIE_XYZ_RANGE_DEFAULT__, /* Spectral range */                          \
  -1, /* Reference temperature */                                              \
  HTRDR_SPECTRAL_SW_CIE_XYZ, /* Spectral type */                               \
}
static const struct htrdr_args_spectral HTRDR_ARGS_SPECTRAL_DEFAULT =
  HTRDR_ARGS_SPECTRAL_DEFAULT__;

/*******************************************************************************
 * Exported functions
 ******************************************************************************/
HTRDR_API res_T
htrdr_args_camera_parse
  (struct htrdr_args_camera* cam, 
   const char* str);

HTRDR_API res_T
htrdr_args_rectangle_parse
  (struct htrdr_args_rectangle* rect,
   const char* str);
 
HTRDR_API res_T
htrdr_args_image_parse
  (struct htrdr_args_image* img,
   const char* str);

HTRDR_API res_T
htrdr_args_spectral_parse
  (struct htrdr_args_spectral* spectral,
   const char* str);

#endif /* HTRDR_ARGS_H */
