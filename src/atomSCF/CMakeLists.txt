# CMake file for OLCAO application: atomSCF

# Document the current location
message("At the atomSCF application.")

add_executable(atomSCF ../kinds.f90 ../constants.f90 ../readData.f90
   ../radialGrid.f90 dataModules.f90 splines.f90 readSCFInput.f90 initRho.f90
   mixer.f90 totalEnergy.f90 implicitSCFInput.f90 computeIonicPot.f90
   exchCorrSubs.f90 computeElectPot.f90 computeOrbitals.f90
   checkConvg.f90 atomSCF.f90)

install(TARGETS atomSCF DESTINATION bin)
