##
##  CMakeFile for FIOs-EDOs.x
##  Vigo, 2018-10-10
##
##  Copyright 2018 Nicolás Otero Martínez and Marcos Mandado Alonso
##  
##  This file and CMake configuration is based on those from 
##  astroTools package, see: http://astrotools.sf.net/
##  
##  FIOs-EDOs.x is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##  
##  FIOs-EDOs.x is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public License
##  along with FIOs-EDOs.x.  If not, see <http://www.gnu.org/licenses/>.
##
##  To install this software on your system, see the INSTALL file
##  or the manual in the root directory of this package.
##


cmake_minimum_required( VERSION 2.8 )
# Set build type. Do this *before* we set the project name:
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo Profile."
    FORCE )
endif( NOT CMAKE_BUILD_TYPE )
set( CMAKE_CONFIGURATION_TYPES "${CMAKE_BUILD_TYPE}" CACHE INTERNAL "internal" )

## >> DEFAULT COMPILER << ## Force compiler you want
#set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_COMPILER ifort)
#set(CMAKE_GENERATOR_FC gfortran)
set(CMAKE_GENERATOR_FC ifort)
message("CMAKE_Fortran_COMPILER is ${CMAKE_Fortran_COMPILER}")
message("CMAKE_GENERATOR_FC is ${CMAKE_GENERATOR_FC}")


# Project name and language:
project( FIO-POL-HYP-mkl.x Fortran )


# Search in the CMake/ directory for CMake modules:
list( APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake )


# Various compile/optimization options that we may want to enable:
include( SetCompileOptions )


# Get the compiler name (for compiler flags and to find libraries):
get_filename_component( Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME )



# Find dependencies:
#find_package( LibSUFR REQUIRED )
#set( INCLUDE_FLAGS "-I${LibSUFR_INCLUDES}" )  # will be transferred to CompilerFlags



# Set FORTRAN compiler flags:
include( CompilerFlags_Fortran )


message( STATUS "Compiling for an ${CMAKE_SYSTEM_PROCESSOR} architecture" )
message( STATUS "" )


# Set source files:
include( FileList )


add_executable( FIOs-EDOs-mkl.x ${SRC_FILES})
#target_link_libraries( EDOs.x  ${LibSUFR_LIBRARIES} )


# Install the binaries:
#install( TARGETS  EDOs.x  RUNTIME DESTINATION bin )

# Install documents:
#if( NOT EXISTS VERSION )
#  install( FILES doc/CHANGELOG doc/LICENCE doc/README doc/VERSION DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/astrotools-${PKG_VERSION} )
#else( NOT EXISTS VERSION )
#  install( FILES CHANGELOG LICENCE README VERSION DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/astrotools-${PKG_VERSION} )
#endif( NOT EXISTS VERSION )

# Install man pages:
#install( DIRECTORY man/ DESTINATION ${CMAKE_INSTALL_PREFIX}/share/man )

