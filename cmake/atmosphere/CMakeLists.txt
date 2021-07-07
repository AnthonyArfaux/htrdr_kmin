# Copyright (C) 2018, 2019, 2020, 2021 |Meso|Star> (contact@meso-star.com)
# Copyright (C) 2018, 2019, 2021 CNRS
# Copyright (C) 2018, 2019 Université Paul Sabatier
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.1)
project(htrdr-atmosphere)

################################################################################
# Check dependencies
################################################################################
find_package(HTSky 0.2 REQUIRED)
find_package(RCMake 0.3 REQUIRED)
find_package(RSys 0.11 REQUIRED)
find_package(Star3D 0.8 REQUIRED)
find_package(StarSF 0.6 REQUIRED)
find_package(StarSP 0.9 REQUIRED)
find_package(StarVX 0.1 REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${RCMAKE_SOURCE_DIR})
include(rcmake)
include(rcmake_runtime)

include_directories(
  ${HTSky_INCLUDE_DIR}
  ${RSys_INCLUDE_DIR}
  ${Star3D_INCLUDE_DIR}
  ${StarSF_INCLUDE_DIR}
  ${StarSP_INCLUDE_DIR}
  ${StarVX_INCLUDE_DIR}
  ${HTRDR_BUILD_DIR}
  ${HTRDR_SOURCE_DIR})

################################################################################
# Generate files
################################################################################
set(HTRDR_ATMOSPHERE_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD "1" CACHE INTERNAL "")
set(HTRDR_ATMOSPHERE_ARGS_DEFAULT_SKY_MTL_NAME "\"air\"" CACHE INTERNAL "")

configure_file(${HTRDR_SOURCE_DIR}/atmosphere/htrdr_atmosphere_args.h.in
  ${HTRDR_BUILD_DIR}/atmosphere/htrdr_atmosphere_args.h @ONLY)

################################################################################
# Configure and define targets
################################################################################
set(HTRDR_ATMOSPHERE_FILES_SRC
  htrdr_atmosphere_args.c
  htrdr_atmosphere.c
  htrdr_atmosphere_compute_radiance_lw.c
  htrdr_atmosphere_compute_radiance_sw.c
  htrdr_atmosphere_draw_map.c
  htrdr_atmosphere_ground.c
  htrdr_atmosphere_main.c
  htrdr_atmosphere_sun.c)

set(HTRDR_ATMOSPHERE_FILES_INC
  htrdr_atmosphere_c.h
  htrdr_atmosphere_ground.h
  htrdr_atmosphere.h
  htrdr_atmosphere_sun.h)

set(HTRDR_ATMOSPHERE_FILES_INC2 # Generated files
  htrdr_atmosphere_args.h)

# Prepend each file in the `HTRDR_FILES_<SRC|INC>' list by `HTRDR_SOURCE_DIR'
rcmake_prepend_path(HTRDR_ATMOSPHERE_FILES_SRC ${HTRDR_SOURCE_DIR}/atmosphere)
rcmake_prepend_path(HTRDR_ATMOSPHERE_FILES_INC ${HTRDR_SOURCE_DIR}/atmosphere)
rcmake_prepend_path(HTRDR_ATMOSPHERE_FILES_INC2 ${HTRDR_BUILD_DIR}/atmosphere)

# Atmosphere library
add_library(htrdr-atmosphere SHARED
  ${HTRDR_ATMOSPHERE_FILES_SRC}
  ${HTRDR_ATMOSPHERE_FILES_INC}
  ${HTRDR_ATMOSPHERE_FILES_INC2})
target_link_libraries(htrdr-atmosphere htrdr-core HTSky RSys Star3D StarSF StarSP)

if(CMAKE_COMPILER_IS_GNUCC)
  target_link_libraries(htrdr-atmosphere m)
endif()

set_target_properties(htrdr-atmosphere PROPERTIES
  DEFINE_SYMBOL HTRDR_SHARED_BUILD
  VERSION ${VERSION}
  SOVERSION ${VERSION_MAJOR})

################################################################################
# Define output & install directories
################################################################################
install(TARGETS htrdr-atmosphere
  ARCHIVE DESTINATION bin
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

