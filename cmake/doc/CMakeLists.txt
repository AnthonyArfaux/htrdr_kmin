# Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
# Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
# Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
# Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
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

set(HTRDR_DOC_DIR ${PROJECT_SOURCE_DIR}/../doc)

################################################################################
# Look for the a2x program
################################################################################
find_program(SCDOC NAMES scdoc)
if(NOT SCDOC)
  message(WARNING
    "The `scdoc' program is missing. "
    "The htrdr man pages cannot be generated.")
  return()
endif()

################################################################################
# ROFF man pages
################################################################################
set(MAN_SOURCES 
  ${HTRDR_DOC_DIR}/htrdr.1.scd
  ${HTRDR_DOC_DIR}/htrdr-image.5.scd
  ${HTRDR_DOC_DIR}/htrdr-materials.5.scd
  ${HTRDR_DOC_DIR}/htrdr-obj.5.scd)

if(HTRDR_BUILD_ATMOSPHERE)
  configure_file(${HTRDR_DOC_DIR}/htrdr-atmosphere.1.scd.in
    ${CMAKE_CURRENT_BINARY_DIR}/htrdr-atmosphere.1.scd @ONLY)
  list(APPEND MAN_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/htrdr-atmosphere.1.scd)
endif()

if(HTRDR_BUILD_COMBUSTION)
  configure_file(${HTRDR_DOC_DIR}/htrdr-combustion.1.scd.in
    ${CMAKE_CURRENT_BINARY_DIR}/htrdr-combustion.1.scd @ONLY)
  list(APPEND MAN_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/htrdr-combustion.1.scd)
endif()

if(HTRDR_BUILD_PLANETO)
  configure_file(${HTRDR_DOC_DIR}/htrdr-planeto.1.scd.in
    ${CMAKE_CURRENT_BINARY_DIR}/htrdr-planeto.1.scd @ONLY)
  list(APPEND MAN_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/htrdr-planeto.1.scd)
endif()

set(MAN_FILES)
set(MAN5_FILES)
set(MAN1_FILES)

foreach(_src IN LISTS MAN_SOURCES)
  get_filename_component(_man ${_src} NAME)
  string(REGEX MATCH "^.*\.[1-9]" _man "${_man}")

  add_custom_command(
    OUTPUT ${_man}
    COMMAND ${SCDOC} < ${_src} > ${_man}
    DEPENDS ${_src}
    COMMENT "Build ROFF man page ${_man}"
    VERBATIM)
  list(APPEND MAN_FILES ${_man})

  string(REGEX MATCH "^.*.5$" _man5 ${_man})
  string(REGEX MATCH "^.*.1$" _man1 ${_man})
  if(_man1)
    list(APPEND MAN1_FILES ${CMAKE_CURRENT_BINARY_DIR}/${_man1})
  elseif(_man5)
    list(APPEND MAN5_FILES ${CMAKE_CURRENT_BINARY_DIR}/${_man5})
  else()
    message(FATAL_ERROR "Unexpected man type")
  endif()
endforeach()
add_custom_target(man-roff ALL DEPENDS ${MAN_FILES})

install(FILES ${MAN1_FILES} DESTINATION share/man/man1)
install(FILES ${MAN5_FILES} DESTINATION share/man/man5)
