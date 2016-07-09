# CMake file for OLCAO application: applySpaceGroup

# Document the current location
message("At the applySpaceGroup application.")

add_executable(applySpaceGroup ../kinds.f90 ../constants.f90
   applySpaceGroup.f90)

install(TARGETS applySpaceGroup DESTINATION bin)
