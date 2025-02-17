# CMake file for OLCAO application: makeKPoints_cv

# Document the current location
message("At the makeKPoints with computational visualization application.")

add_executable(makekpoints_cv ../../../../kinds.f90 ../../../../constants.f90
   ../makekpoints.F90)
add_definitions(-DCOMPVIS)

install(TARGETS makekpoints_cv DESTINATION bin)
