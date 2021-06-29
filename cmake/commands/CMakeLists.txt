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
project(htrdr-commands C)

include_directories(${HTRDR_BUILD_DIR} ${HTRDR_SOURCE_DIR})

set(_link_libraries "")

if(HTRDR_BUILD_ATMOSPHERE)
  list(APPEND _link_libraries htrdr-atmosphere)
endif()
if(HTRDR_BUILD_COMBUSTION)
  list(APPEND _link_libraries htrdr-combustion)
endif()

################################################################################
# Check dependencies
################################################################################
find_package(RSys 0.11 REQUIRED)
include_directories(${RSys_INCLUDE_DIR})

################################################################################
# Configure and define targets
################################################################################
add_executable(htrdr_cmd ${HTRDR_SOURCE_DIR}/commands/htrdr_cmd.c)
target_link_libraries(htrdr_cmd ${_link_libraries})
set_target_properties(htrdr_cmd PROPERTIES
  COMPILE_DEFINITIONS "${HTRDR_BUILD_COMMANDS}"
  OUTPUT_NAME htrdr)

add_executable(htrdr_atmosphere_cmd
  ${HTRDR_SOURCE_DIR}/commands/htrdr_atmosphere_cmd.c)
set_target_properties(htrdr_atmosphere_cmd PROPERTIES
  COMPILE_DEFINITIONS "${HTRDR_BUILD_COMMANDS}"
  OUTPUT_NAME htrdr-atmosphere)
if(HTRDR_BUILD_ATMOSPHERE)
  target_link_libraries(htrdr_atmosphere_cmd htrdr-atmosphere)
endif()

add_executable(htrdr_combustion_cmd
  ${HTRDR_SOURCE_DIR}/commands/htrdr_combustion_cmd.c)
set_target_properties(htrdr_combustion_cmd PROPERTIES
  COMPILE_DEFINITIONS "${HTRDR_BUILD_COMMANDS}"
  OUTPUT_NAME htrdr-combustion)
if(HTRDR_BUILD_COMBUSTION)
  target_link_libraries(htrdr_combustion_cmd htrdr-combustion)
endif()

################################################################################
# Define output & install directories
################################################################################
install(TARGETS htrdr_cmd htrdr_atmosphere_cmd htrdr_combustion_cmd
  ARCHIVE DESTINATION bin
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)
