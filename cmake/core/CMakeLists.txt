# Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
# Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
# Copyright (C) 2022-2023 Institut Pierre-Simon Laplace
# Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
# Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
# Copyright (C) 2022-2023 Observatoire de Paris
# Copyright (C) 2022-2023 Université de Reims Champagne-Ardenne
# Copyright (C) 2022-2023 Université de Versaille Saint-Quentin
# Copyright (C) 2018-2019, 2022-2023 Université Paul Sabatier
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
project(htrdr-core C)

################################################################################
# Check dependencies
################################################################################
find_package(AW 2.0 REQUIRED)
find_package(MruMtl 0.1 REQUIRED)
find_package(RCMake 0.3 REQUIRED)
find_package(RSys 0.11 REQUIRED)
find_package(Star3D 0.8 REQUIRED)
find_package(StarCamera 0.0 REQUIRED)
find_package(StarSF 0.6 REQUIRED)
find_package(StarSP 0.12 REQUIRED)
find_package(OpenMP 1.2 REQUIRED)
find_package(MPI 2 REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${RCMAKE_SOURCE_DIR})
include(rcmake)
include(rcmake_runtime)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

include_directories(
  ${AW_INCLUDE_DIR}
  ${MruMtl_INCLUDE_DIR}
  ${RSys_INCLUDE_DIR}
  ${Star3D_INCLUDE_DIR}
  ${StarCamera_INCLUDE_DIR}
  ${StarSF_INCLUDE_DIR}
  ${StarSP_INCLUDE_DIR}
  ${StarVX_INCLUDE_DIR}
  ${HTSky}
  ${HTRDR_BUILD_DIR}
  ${HTRDR_SOURCE_DIR}
  ${MPI_INCLUDE_PATH})

################################################################################
# Generate files
################################################################################
set(HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MIN "0.0" CACHE INTERNAL "")
set(HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MAX "180.0" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_CAMERA_POS "0,0,0" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_CAMERA_TGT "0,1,0" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_CAMERA_UP  "0,0,1" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_CAMERA_ORTHOGRAPHIC_HEIGHT "1" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOV "70" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_LENS_RADIUS "0" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOCAL_DST "1" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_RECTANGLE_POS "0,0,0" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_RECTANGLE_TGT "0,0,1" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_RECTANGLE_UP "0,1,0" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_RECTANGLE_SZ "1,1" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_IMG_WIDTH "320" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_IMG_HEIGHT "240" CACHE INTERNAL "")
set(HTRDR_ARGS_DEFAULT_IMG_SPP "1" CACHE INTERNAL "")

configure_file(${HTRDR_SOURCE_DIR}/core/htrdr_args.h.in
  ${HTRDR_BUILD_DIR}/core/htrdr_args.h @ONLY)
configure_file(${HTRDR_SOURCE_DIR}/core/htrdr_version.h.in
  ${HTRDR_BUILD_DIR}/core/htrdr_version.h @ONLY)

################################################################################
# Configure and define targets
################################################################################
set(HTRDR_CORE_FILES_SRC
  htrdr.c
  htrdr_args.c
  htrdr_buffer.c
  htrdr_draw_map.c
  htrdr_geometry.c
  htrdr_log.c
  htrdr_materials.c
  htrdr_ran_wlen_cie_xyz.c
  htrdr_ran_wlen_discrete.c
  htrdr_ran_wlen_planck.c
  htrdr_rectangle.c
  htrdr_slab.c
  htrdr_spectral.c)
set(HTRDR_CORE_FILES_INC
  htrdr.h
  htrdr_accum.h
  htrdr_buffer.h
  htrdr_c.h
  htrdr_draw_map.h
  htrdr_geometry.c
  htrdr_interface.h
  htrdr_log.h
  htrdr_materials.h
  htrdr_ran_wlen_cie_xyz.h
  htrdr_ran_wlen_discrete.h
  htrdr_ran_wlen_planck.h
  htrdr_rectangle.h
  htrdr_sensor.h
  htrdr_slab.h
  htrdr_spectral.h)
set(HTRDR_CORE_FILES_INC2 # Generated files
  htrdr_args.h
  htrdr_version.h)

# Prepend each file in the `HTRDR_FILES_<SRC|INC>' list by `HTRDR_SOURCE_DIR'
rcmake_prepend_path(HTRDR_CORE_FILES_SRC ${HTRDR_SOURCE_DIR}/core)
rcmake_prepend_path(HTRDR_CORE_FILES_INC ${HTRDR_SOURCE_DIR}/core)
rcmake_prepend_path(HTRDR_CORE_FILES_INC2 ${HTRDR_BUILD_DIR}/core)

# Core library
add_library(htrdr-core SHARED
  ${HTRDR_CORE_FILES_SRC}
  ${HTRDR_CORE_FILES_INC}
  ${HTRDR_CORE_FILES_INC2})
target_link_libraries(htrdr-core AW MruMtl RSys Star3D StarCamera StarSF StarSP)

if(CMAKE_COMPILER_IS_GNUCC)
  target_link_libraries(htrdr-core m)
  set_target_properties(htrdr-core PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")
endif()

set_target_properties(htrdr-core PROPERTIES
  COMPILE_FLAGS "${OpenMP_C_FLAGS}"
  DEFINE_SYMBOL HTRDR_CORE_SHARED_BUILD
  VERSION ${VERSION}
  SOVERSION ${VERSION_MAJOR})

################################################################################
# Define output & install directories
################################################################################
install(TARGETS htrdr-core
  ARCHIVE DESTINATION bin
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

