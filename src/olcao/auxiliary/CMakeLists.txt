# CMake file for OLCAO application: auxiliary OLCAO

# Document the current location
message("At the auxiliary OLCAO application level.")

add_executable(OLCAOkkc ../../kinds.f90 ../../constants.f90 OLCAOkkc.f90)
add_executable(OLCAOrhoV ../../kinds.f90 ../../constants.f90 OLCAOrhoV.f90)
add_executable(num2CKE ../../kinds.f90 numerical2CKineticE.f90)
add_executable(num2CMom ../../kinds.f90 numerical2CMomentum.f90)
add_executable(num2COvlp ../../kinds.f90 numerical2COverlap.f90)
add_executable(num3CNuc ../../kinds.f90 numerical3CNuclear.f90)
add_executable(num3COvlp ../../kinds.f90 numerical3COverlap.f90)

install(TARGETS OLCAOkkc OLCAOrhoV num2CKE num2CMom num2COvlp num3CNuc
      num3COvlp DESTINATION bin)
