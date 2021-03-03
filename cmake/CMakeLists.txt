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
project(htrdr C)

set(HTRDR_SOURCE_DIR ${PROJECT_SOURCE_DIR}/../src)
set(HTRDR_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(VERSION_MAJOR 0)
set(VERSION_MINOR 7)
set(VERSION_PATCH 0)
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

################################################################################
# Add sub projects
################################################################################
add_subdirectory(core)
add_subdirectory(atmosphere)
add_subdirectory(commands)
add_subdirectory(combustion)
add_subdirectory(doc)

################################################################################
# Define output & install directories
################################################################################
set(HTRDR_FILES_DOC
  ${PROJECT_SOURCE_DIR}/../COPYING
  ${PROJECT_SOURCE_DIR}/../README.md)
install(FILES ${HTRDR_FILES_DOC} DESTINATION share/doc/htrdr)

