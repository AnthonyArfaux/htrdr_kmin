# Copyright (C) 2018 CNRS, Université Paul Sabatier, |Meso|Star>
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
find_program(A2X NAMES a2x a2x.py)
if(NOT A2X)
  unset(_roff)
  message(WARNING
    "The `a2x' program is missing. "
    "The htrdr man pages cannot be generated.")
endif()

################################################################################
# Copy doc files
################################################################################
set(MAN_NAMES htrdr-image.5)

set(MAN_FILES)
foreach(_name IN LISTS MAN_NAMES)
  set(_src ${HTRDR_DOC_DIR}/${_name}.txt)
  set(_dst ${CMAKE_CURRENT_BINARY_DIR}/${_name}.txt)
  add_custom_command(
    OUTPUT ${_dst}
    COMMAND ${CMAKE_COMMAND} -E copy ${_src} ${_dst}
    DEPENDS ${_src}
    COMMENT "Copy the asciidoc ${_src}"
    VERBATIM)
  list(APPEND MAN_FILES ${_dst})
endforeach()
add_custom_target(man-copy ALL DEPENDS ${MAN_FILES})

configure_file(${HTRDR_SOURCE_DIR}/../doc/htrdr.1.txt.in
  ${CMAKE_CURRENT_BINARY_DIR}/htrdr.1.txt @ONLY)
list(APPEND MAN_NAMES htrdr.1)

################################################################################
# ROFF man pages
################################################################################
set(A2X_OPTS -dmanpage -fmanpage)
set(MAN_FILES)
set(MAN5_FILES)
set(MAN1_FILES)
foreach(_name IN LISTS MAN_NAMES)
  set(_man ${CMAKE_CURRENT_BINARY_DIR}/${_name})
  set(_txt ${CMAKE_CURRENT_BINARY_DIR}/${_name}.txt)

  add_custom_command(
    OUTPUT ${_man}
    COMMAND ${A2X} ${A2X_OPTS} ${_txt}
    DEPENDS man-copy ${_txt}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Build ROFF man page ${_man}"
    VERBATIM)
  list(APPEND MAN_FILES ${_man})

  string(REGEX MATCH "^.*.5$" _man5 ${_man})
  string(REGEX MATCH "^.*.1$" _man1 ${_man})
  if(_man1)
    list(APPEND MAN1_FILES ${_man1})
  elseif(_man5)
    list(APPEND MAN5_FILES ${_man5})
  else()
    message(FATAL_ERROR "Unexpected man type")
  endif()
endforeach()
add_custom_target(man-roff ALL DEPENDS ${MAN_FILES})

install(FILES ${MAN1_FILES} DESTINATION share/man/man1)
install(FILES ${MAN5_FILES} DESTINATION share/man/man5)

