# SPDX-FileCopyrightText: 2024 Julien Y. Dutheil
#
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.7)
project(Nisce CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Compile options
add_compile_options(-Wall -Weffc++ -Wshadow -Wconversion)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)
#static linkage?
if(NOT BUILD_STATIC)
  set(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
endif()
if(BUILD_STATIC)
  message(STATUS "Static linkage requested.")
  add_compile_options(-static -static-libgcc)
endif()


include (GNUInstallDirs)
# Find dependencies (add install directory to search)
if (CMAKE_INSTALL_PREFIX)
  set (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif (CMAKE_INSTALL_PREFIX)

find_package (bpp-core3 1.0.0 REQUIRED)

#Find boost libraries
#SET(Boost_USE_STATIC_LIBS ${BUILD_STATIC})
SET(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE( Boost 1.34.0 COMPONENTS math_c99 )

IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  SET(LIBS ${LIBS} ${Boost_LIBRARIES})
  MESSAGE("-- Boost libraries found here:")
  MESSAGE("   includes: ${Boost_INCLUDE_DIRS}")
  MESSAGE("   libraries: ${Boost_LIBRARIES}")
ENDIF()

FIND_PACKAGE( TBB REQUIRED )

IF(TBB_FOUND)
  INCLUDE_DIRECTORIES(${tbb_INCLUDE_DIRS})
  SET(LIBS ${LIBS} TBB::tbb)
  MESSAGE("-- TBB libraries found.")
ENDIF()

# Subdirectories
ADD_SUBDIRECTORY(Nisce)

