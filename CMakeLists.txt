# Top level CMake file for OLCAO

# Set the minimum required version, the project name, and enable Fortran.
cmake_minimum_required(VERSION 3.1.0)
project(OLCAO)
enable_language (Fortran)

# Make sure that the default build type is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Establish the installation location on the basis of the $OLCAO_DIR
#   environment variable that is defined in the .olcao file. Make sure that you
#   source this file (source .olcao) before you try to build.
set (CMAKE_INSTALL_PREFIX $ENV{OLCAO_DIR})
message("CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")

# Get the full path of the compiler from the FC environment variable and
#   extract the compiler name from that. Then report the name and the ID (which
#   is found by cmake itself). The idea is that FC may be equal to h5fc but
#   the underlying compiler that h5fc wraps may be gfortran, intel, pgf90,
#   etc. We need all of those pieces of information. (Note that FC should be
#   defined before attempting to build.)
# Another tricky note is that the Fortran compiler that was used to compile
#   the HDF5 library must also be used as the compiler for OLCAO. Otherwise the
#   module files will not be properly recognized. (At least it is the case that
#   gfortran and intel cannot read eachother's .mod files. I.e. the HDF5.mod
#   file that is created when HDF5 is compiled must be readable by the compiler
#   that you are using to compile OLCAO.)
get_filename_component (CMAKE_Fortran_COMPILER $ENV{FC} PROGRAM)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("Fortran_COMPILER_NAME= ${Fortran_COMPILER_NAME}")
message("CMAKE_Fortran_COMPLIER= ${CMAKE_Fortran_COMPILER}")
message("CMAKE_Fortran_COMPILER_ID= ${CMAKE_Fortran_COMPILER_ID}")

# Define the compiler and linker flags that are needed for the different
#   compilers. So far only gfortran and ifort are considered. You can just
#   change the appropriate line or add a new elseif block as you see fit if
#   you have a particular compiler requirement. I think that this might be
#   easier than attempting to read the user's $FCFLAGS and $LDFLAGS variables
#   because the user may change those to compile other programs and may not
#   retain a record of what they were when he/she compiles OLCAO. If the user
#   stores it here then he/she will have a record. I guess just be careful to
#   retain that information if you get a new version from github. Or something.
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fimplicit-none -Wall")
   set (CMAKE_Fortran_FLAGS_DEBUG "-g -fimplicit-none -Wall")
   set (CMAKE_EXE_LINKER_FLAGS "-llapack -lblas")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   set (CMAKE_Fortran_FLAGS_RELEASE "-ipo -O3 -no-prec-div -fp-model fast=2 -heap-arrays")
#   set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -heap-arrays")
   set (CMAKE_Fortran_FLAGS_DEBUG "-g -check all -traceback -fp-stack-check -heap-arrays")
   set (CMAKE_EXE_LINKER_FLAGS "-mkl=sequential")
endif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
message("CMAKE_Fortran_FLAGS_RELEASE= ${CMAKE_Fortran_FLAGS_RELEASE}")
message("CMAKE_Fortran_FLAGS_DEBUG= ${CMAKE_Fortran_FLAGS_DEBUG}")
message("CMAKE_EXE_LINKER_FLAGS= ${CMAKE_EXE_LINKER_FLAGS}")

# Proceed to the source directory and its subdirectories.
add_subdirectory(src)
