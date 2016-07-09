# CMake file for OLCAO application: auxiliary OLCAO

# Document the current location
message("At the auxiliary OLCAO application level.")

add_executable(OLCAOkkc ../../kinds.f90 ../../constants.f90 OLCAOkkc.f90)
add_executable(OLCAOrhoV ../../kinds.f90 ../../constants.f90 OLCAOrhoV.f90)

install(TARGETS OLCAOkkc OLCAOrhoV DESTINATION bin)
