# CMake file for OLCAO application: auxiliary OLCAO

# Document the current location
message("At the auxiliary upOLCAO application level.")

add_executable(OLCAOkkc ../../kinds.f90 ../../constants.f90 OLCAOkkc.f90)
add_executable(OLCAOrhoV ../../kinds.f90 ../../constants.f90 OLCAOrhoV.f90)
add_executable(gaussIntgTest ../../kinds.f90 ../../constants.f90 gaussIntgTest.f90)
add_subdirectory(applySpaceGroup)
add_subdirectory(atomSCF)
add_subdirectory(contract)
add_subdirectory(gaussFit)
add_subdirectory(makeKPoints)

install(TARGETS OLCAOkkc OLCAOrhoV gaussIntgTest DESTINATION bin)
