# SPDX-FileCopyrightText: 2024 Julien Y. Dutheil
#
# SPDX-License-Identifier: GPL-3.0-or-later

# CMake script for NISCE
# Author: Julien Dutheil
# Created: 25/04/2023

set(SOURCE_FILES nisce.cpp)
add_executable(nisce ${SOURCE_FILES})

# Link (static or shared)
if (BUILD_STATIC)
  target_link_libraries (nisce ${BPP_LIBS_STATIC})
  set_target_properties (nisce PROPERTIES LINK_SEARCH_END_STATIC TRUE)
else (BUILD_STATIC)
  target_link_libraries (nisce ${BPP_LIBS_SHARED})
  set_target_properties (nisce PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif (BUILD_STATIC)

install (TARGETS nisce DESTINATION ${CMAKE_INSTALL_BINDIR})
