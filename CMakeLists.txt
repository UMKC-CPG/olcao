# Top level CMake file for OLCAO

cmake_minimum_required(VERSION 2.8)
project(OLCAO)
enable_language (Fortran)

# Make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Establish the installation location
set (CMAKE_INSTALL_PREFIX $ENV{OLCAO_DIR})

# Get the full path of the compiler from the FC environment variable and
#   extract the compiler name from that. Then report the name and the ID (which
#   is found by cmake itself). The idea is that FC may be equal to h5fc but
#   the underlying compiler may be gfortran or intel etc. We need all of those
#   pieces of information.
get_filename_component (CMAKE_Fortran_COMPILER $ENV{FC} PROGRAM)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("Fortran_COMPILER_NAME= ${CMAKE_Fortran_COMPILER_ID}")
message("CMAKE_Fortran_COMPILER_ID= ${CMAKE_Fortran_COMPILER_ID}")

if (CMAKE_Fortran_COMPILER_ID MATCHES "gfortran")
   set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
   set (CMAKE_EXE_LINKER_FLAGS "-llapack -lblas")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
   set (CMAKE_EXE_LINKER_FLAGS "-mkl")
endif (CMAKE_Fortran_COMPILER_ID MATCHES "gfortran")
message("CMAKE_Fortran_FLAGS_RELEASE= ${CMAKE_Fortran_FLAGS_RELEASE}")

message("CMAKE_Fortran_COMPLIER= ${CMAKE_Fortran_COMPILER}")
message("Fortran_COMPILER_NAME= ${Fortran_COMPILER_NAME}")
message("CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")

# Proceed to the source directory and its subdirectories.
add_subdirectory(src)


