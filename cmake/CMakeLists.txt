# Copyright (C) 2018 Université Paul Sabatier, |Meso|Star>
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
find_package(RCMake 0.3 REQUIRED)
find_package(RSys 0.6 REQUIRED)
find_package(StarSP 0.8 REQUIRED)
find_package(StarVX REQUIRED)
find_package(HTCP REQUIRED)
find_package(HTMIE REQUIRED)
find_package(OpenMP 1.2 REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${RCMAKE_SOURCE_DIR})
include(rcmake)
include(rcmake_runtime)

include_directories(
  ${RSys_INCLUDE_DIR}
  ${SVX_INCLUDE_DIR}
  ${HTCP_INCLUDE_DIR})

################################################################################
# Configure and define targets
################################################################################
set(VERSION_MAJOR 0)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

set(HTRDR_FILES_SRC
  htrdr.c
  htrdr_args.c
  htrdr_buffer.c
  htrdr_main.c
  htrdr_rectangle.c
  htrdr_sky.c
  htrdr_sun.c
  htrdr_transmission.c)
set(HTRDR_FILES_INC
  htrdr.h
  htrdr_c.h
  htrdr_args.h
  htrdr_buffer.h
  htrdr_rectangle.h
  htrdr_sky.h
  htrdr_sun.h
  htrdr_solve.h)
set(HTDRD_FILES_DOC COPYING README.md)

# Prepend each file in the `HTRDR_FILES_<SRC|INC>' list by `HTRDR_SOURCE_DIR'
rcmake_prepend_path(HTRDR_FILES_SRC ${HTRDR_SOURCE_DIR})
rcmake_prepend_path(HTRDR_FILES_INC ${HTRDR_SOURCE_DIR})
rcmake_prepend_path(HTRDR_FILES_DOC ${PROJECT_SOURCE_DIR}/../)

add_executable(htrdr ${HTRDR_FILES_SRC} ${HTRDR_FILES_INC})
target_link_libraries(htrdr HTCP HTMIE RSys SVX StarSP)

if(CMAKE_COMPILER_IS_GNUCC)
  target_link_libraries(htrdr m)
  set_target_properties(htrdr PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif()

set_target_properties(htrdr PROPERTIES
  DEFINE_SYMBOL HTRDR_SHARED_BUILD
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

