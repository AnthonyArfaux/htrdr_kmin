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

#ifndef HTRDR_SENSOR_H
#define HTRDR_SENSOR_H

#include <rsys/rsys.h>

/* Forward declarations */
struct htrdr_ground;
struct ssp_rng;

enum htrdr_sensor_type {
  HTRDR_SENSOR_CAMERA,
  HTRDR_SENSOR_RECTANGLE
};

struct htrdr_sensor {
  struct htrdr_camera* camera;
  struct htrdr_rectangle* rectangle;
  enum htrdr_sensor_type type;
};

extern LOCAL_SYM res_T
htrdr_sensor_sample_primary_ray
  (const struct htrdr_sensor* sensor,
   struct htrdr_ground* ground,
   const size_t ipix[2],
   const double pix_sz[2],
   struct ssp_rng* rng,
   double ray_org[3],
   double ray_dir[3]);

#endif /* HTRDR_SENSOR_H */

