# CMake file for OLCAO application: gaussFit

# Document the current location
message("At the gaussFit application.")

add_executable(gaussFit ../kinds.f90 ../readData.f90 dataModules.f90
   gaussFit.f90)

install(TARGETS gaussFit DESTINATION bin)
