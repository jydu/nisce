# SPDX-FileCopyrightText: 2024 Julien Y. Dutheil
#
# SPDX-License-Identifier: GPL-3.0-or-later

# CMake script for NISCE
# Author: Julien Dutheil
# Created: 25/04/2023

set(SOURCE_FILES nisce.cpp)
add_executable(nisce ${SOURCE_FILES})

target_link_libraries (nisce ${BPP_LIBS_SHARED})
target_link_libraries (nisce ${LIBS})
install (TARGETS nisce DESTINATION ${CMAKE_INSTALL_BINDIR})
