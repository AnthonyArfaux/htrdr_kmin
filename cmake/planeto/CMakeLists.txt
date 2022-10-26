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
project(htrdr-planeto)

################################################################################
# Check dependencies
################################################################################
find_package(RCMake 0.4 REQUIRED)
find_package(RNATM 0.0 REQUIRED)
find_package(RNGRD 0.0 REQUIRED)
find_package(RSys 0.11 REQUIRED)
find_package(Star3D 0.8 REQUIRED)
find_package(StarCamera 0.0 REQUIRED)
find_package(StarSF 0.6 REQUIRED)
find_package(StarSP 0.12 REQUIRED)
find_package(StarVX 0.1 REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${RCMAKE_SOURCE_DIR})
include(rcmake)
include(rcmake_runtime)

include_directories(
  ${RNATM_INCLUDE_DIR}
  ${RNGRD}
  ${RSys_INCLUDE_DIR}
  ${Star3D_INCLUDE_DIR}
  ${StarCamera_INCLUDE_DIR}
  ${StarSF_INCLUDE_DIR}
  ${StarSP_INCLUDE_DIR}
  ${StarVX_INCLUDE_DIR}
  ${HTRDR_BUILD_DIR}
  ${HTRDR_SOURCE_DIR})

################################################################################
# Configure and define targets
################################################################################
set(HTRDR_PLANETO_FILES_SRC
  htrdr_planeto.c
  htrdr_planeto_args.c
  htrdr_planeto_main.c)

set(HTRDR_PLANETO_FILES_INC
  htrdr_planeto.h
  htrdr_planeto_c.h
  htrdr_planeto_args.h)

# Prepend each file in the `HTRDR_FILES_<SRC|INC>' list by `HTRDR_SOURCE_DIR'
rcmake_prepend_path(HTRDR_PLANETO_FILES_SRC ${HTRDR_SOURCE_DIR}/planeto)
rcmake_prepend_path(HTRDR_PLANETO_FILES_INC ${HTRDR_SOURCE_DIR}/planeto)

# Planeto library
add_library(htrdr-planeto SHARED
  ${HTRDR_PLANETO_FILES_SRC}
  ${HTRDR_PLANETO_FILES_INC})
target_link_libraries(htrdr-planeto
  htrdr-core RNATM RNGRD RSys Star3D StarCamera StarSF StarSP)

set_target_properties(htrdr-planeto PROPERTIES
  DEFINE_SYMBOL HTRDR_SHARED_BUILD
  VERSION ${VERSION}
  SOVERSION ${VERSION_MAJOR})

################################################################################
# Define output & install directories
################################################################################
install(TARGETS htrdr-planeto
  ARCHIVE DESTINATION bin
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

