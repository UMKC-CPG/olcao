# CMake file for OLCAO application: makeKPoints

# Document the current location
message("At the makeKPoints application.")

add_executable(makeKPoints ../kinds.f90 ../constants.f90 makekpoints.f90)

install(TARGETS makeKPoints DESTINATION bin)
