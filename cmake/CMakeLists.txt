# Copyright (C) 2018, 2019, 2020 |Meso|Star> (contact@meso-star.com)
# Copyright (C) 2018, 2019 CNRS, Université Paul Sabatier
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

cmake_minimum_required(VERSION 2.8)
project(htrdr C)
enable_testing()

set(HTRDR_SOURCE_DIR ${PROJECT_SOURCE_DIR}/../src)
option(NO_TEST "Do not build the tests" OFF)

################################################################################
# Check dependencies
################################################################################
find_package(AW 2.0 REQUIRED)
find_package(HTSky 0.2 REQUIRED)
find_package(MruMtl 0.0 REQUIRED)
find_package(RCMake 0.3 REQUIRED)
find_package(RSys 0.7 REQUIRED)
find_package(Star3D 0.6 REQUIRED)
find_package(StarSF 0.6 REQUIRED)
find_package(StarSP 0.8 REQUIRED)
find_package(StarVX 0.1 REQUIRED)
find_package(OpenMP 1.2 REQUIRED)
find_package(MPI 1 REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${RCMAKE_SOURCE_DIR})
include(rcmake)
include(rcmake_runtime)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

include_directories(
  ${AW_INCLUDE_DIR}
  ${MruMtl_INCLUDE_DIR}
  ${RSys_INCLUDE_DIR}
  ${Star3D_INCLUDE_DIR}
  ${StarSF_INCLUDE_DIR}
  ${StarSP_INCLUDE_DIR}
  ${StarVX_INCLUDE_DIR}
  ${HTSky}
  ${HTRDR_SOURCE_DIR}
  ${MPI_INCLUDE_PATH}
  ${CMAKE_CURRENT_BINARY_DIR})

################################################################################
# Generate files
################################################################################
set(VERSION_MAJOR 0)
set(VERSION_MINOR 4)
set(VERSION_PATCH 0)
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

set(HTRDR_ARGS_DEFAULT_CAMERA_POS "0,0,0")
set(HTRDR_ARGS_DEFAULT_CAMERA_TGT "0,1,0")
set(HTRDR_ARGS_DEFAULT_CAMERA_UP  "0,0,1")
set(HTRDR_ARGS_DEFAULT_CAMERA_FOV "70")
set(HTRDR_ARGS_DEFAULT_IMG_WIDTH "320")
set(HTRDR_ARGS_DEFAULT_IMG_HEIGHT "240")
set(HTRDR_ARGS_DEFAULT_IMG_SPP "1")
set(HTRDR_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD "1")

configure_file(${HTRDR_SOURCE_DIR}/htrdr_version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/htrdr_version.h @ONLY)

configure_file(${HTRDR_SOURCE_DIR}/htrdr_args.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/htrdr_args.h @ONLY)

################################################################################
# Add sub projects
################################################################################
add_subdirectory(doc)

################################################################################
# Configure and define targets
################################################################################
set(HTRDR_FILES_SRC
  htrdr.c
  htrdr_args.c
  htrdr_buffer.c
  htrdr_camera.c
  htrdr_cie_xyz.c
  htrdr_compute_radiance_sw.c
  htrdr_compute_radiance_lw.c
  htrdr_draw_radiance.c
  htrdr_grid.c
  htrdr_ground.c
  htrdr_interface.c
  htrdr_main.c
  htrdr_mtl.c
  htrdr_ran_wlen.c
  htrdr_rectangle.c
  htrdr_slab.c
  htrdr_sun.c)
set(HTRDR_FILES_INC
  htrdr.h
  htrdr_c.h
  htrdr_buffer.h
  htrdr_camera.h
  htrdr_cie_xyz.h
  htrdr_grid.h
  htrdr_ground.h
  htrdr_interface.h
  htrdr_mtl.h
  htrdr_ran_wlen.h
  htrdr_rectangle.h
  htrdr_slab.h
  htrdr_sun.h
  htrdr_solve.h)
set(HTRDR_FILES_INC2 # Generated files
  htrdr_args.h
  htrdr_version.h)
set(HTRDR_FILES_DOC COPYING README.md)

# Prepend each file in the `HTRDR_FILES_<SRC|INC>' list by `HTRDR_SOURCE_DIR'
rcmake_prepend_path(HTRDR_FILES_SRC ${HTRDR_SOURCE_DIR})
rcmake_prepend_path(HTRDR_FILES_INC ${HTRDR_SOURCE_DIR})
rcmake_prepend_path(HTRDR_FILES_INC2 ${CMAKE_CURRENT_BINARY_DIR})
rcmake_prepend_path(HTRDR_FILES_DOC ${PROJECT_SOURCE_DIR}/../)

add_executable(htrdr ${HTRDR_FILES_SRC} ${HTRDR_FILES_INC} ${HTRDR_FILES_INC2})
target_link_libraries(htrdr AW HTSky MruMtl RSys Star3D StarSF StarSP)

if(CMAKE_COMPILER_IS_GNUCC)
  target_link_libraries(htrdr m)
  set_target_properties(htrdr PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")
endif()

set_target_properties(htrdr PROPERTIES
  COMPILE_FLAGS "${OpenMP_C_FLAGS}"
  VERSION ${VERSION}
  SOVERSION ${VERSION_MAJOR})

################################################################################
# Define output & install directories
################################################################################
install(TARGETS htrdr
  ARCHIVE DESTINATION bin
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)
install(FILES ${HTRDR_FILES_DOC} DESTINATION share/doc/htrdr)

