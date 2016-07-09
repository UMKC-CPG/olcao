# CMake file for OLCAO application: contract

# Document the current location
message("At the contract application.")

add_executable(contract ../kinds.f90 ../constants.f90 ../readData.f90
   ../writeData.f90 ../radialGrid.f90 dataModules.f90 matrixElements.f90
   chargeDensity.f90 matrixDiag.f90 normalization.f90 printResults.f90
   parseInput.f90 contract.f90)

install(TARGETS contract DESTINATION bin)
