# CMake file for OLCAO application: real OLCAO

# Document the current location
message("At the real OLCAO application.")

# Library for gaussIntegrals.f90
add_library(gGaussIntegrals STATIC ../../kinds.f90 ../../constants.f90 
  ../../elementData.f90 ../gaussIntegrals.f90)

# Library containing common source files that all executables need
add_library(gCommon OBJECT ../../kinds.f90 ../../constants.f90 
   ../../elementData.f90 ../../readData.f90 ../stringSubs.f90 ../mathSubs.f90 
   ../sortSubs.f90 ../interfaces.F90 ../timeStamps.f90 ../commandLine.f90 
   ../potSites.f90 ../potTypes.f90 ../atomicTypes.f90 ../atomicSites.f90 
   ../potential.f90 ../lattice.f90 ../kpoints.f90 ../gaussRelations.f90 
   ../hdf5SetupIntg.F90 ../hdf5SetupElec.f90 ../hdf5SetupExco.f90
   ../hdf5Setup.f90 ../exchCorr.f90 ../input.f90)

add_executable(gOLCAOsetup $<TARGET_OBJECTS:gCommon> ../basis.f90 
   ../elecStat.f90 ../intgOrtho.F90 ../intgSaving.F90 ../integralsSCF.F90 
   ../coreCharge.f90 ../OLCAO_Setup.f90 ../OLCAO_SetupWrap.f90 )
target_link_libraries(gOLCAOsetup gGaussIntegrals)

add_executable(gOLCAOmain $<TARGET_OBJECTS:gCommon> ../hdf5MainEval.f90 
   ../hdf5MainEvec.F90 ../hdf5MainPotRho.f90 ../hdf5Main.f90 
   ../hdf5PSCFBand.F90 ../elecStat.f90 ../matrixSubs.F90 ../secularEqn.F90 
   ../populate.F90 ../valeCharge.F90 ../potentialUpdate.F90 ../dos.F90 
   ../bond.F90 ../bond3C.F90 ../OLCAO_Main.f90 ../OLCAO_MainWrap.f90 )

add_executable(gOLCAOintg $<TARGET_OBJECTS:gCommon> ../basis.f90 
   ../hdf5PSCFIntg.f90 ../intgOrtho.F90 ../intgSaving.F90 ../integralsPSCF.F90 
   ../OLCAO_Intg.f90 ../OLCAO_IntgWrap.f90 )
target_link_libraries(gOLCAOintg gGaussIntegrals)

add_executable(gOLCAOband $<TARGET_OBJECTS:gCommon> ../basis.f90 
   ../hdf5PSCFBand.F90 ../hdf5PSCFIntg.f90 ../hdf5MainEval.f90 
   ../hdf5MainEvec.F90 ../hdf5MainPotRho.f90 ../hdf5Main.f90 ../intgOrtho.F90 
   ../intgSaving.F90 ../integralsPSCF.F90 ../matrixSubs.F90 ../secularEqn.F90 
   ../populate.F90 ../OLCAO_Band.F90 ../OLCAO_BandWrap.f90 )
target_link_libraries(gOLCAOband gGaussIntegrals)

add_executable(gOLCAOdos $<TARGET_OBJECTS:gCommon> ../hdf5PSCFBand.F90 
   ../hdf5MainEval.f90 ../hdf5MainEvec.F90 ../hdf5MainPotRho.f90 
   ../hdf5Main.f90 ../matrixSubs.F90 ../secularEqn.F90 ../populate.F90 
   ../dos.F90 ../OLCAO_DOS.f90 )

add_executable(gOLCAObond $<TARGET_OBJECTS:gCommon> ../hdf5PSCFBand.F90 
   ../hdf5MainEval.f90 ../hdf5MainEvec.F90 ../hdf5MainPotRho.f90 
   ../hdf5Main.f90 ../matrixSubs.F90 ../secularEqn.F90 ../populate.F90 
   ../bond.F90 ../bond3C.F90 ../OLCAO_Bond.f90 )

add_executable(gOLCAOoptc $<TARGET_OBJECTS:gCommon> ../basis.f90 
   ../hdf5PSCFBand.F90 ../hdf5PSCFIntg.f90 ../hdf5MainEval.f90 
   ../hdf5MainEvec.F90 ../hdf5MainPotRho.f90 ../hdf5Main.f90 ../intgOrtho.F90 
   ../intgSaving.F90 ../integralsPSCF.F90 ../matrixSubs.F90 ../secularEqn.F90 
   ../populate.F90 ../optc.F90 ../optcPrint.F90 ../OLCAO_Optc.f90 
   ../OLCAO_OptcWrap.f90 )
target_link_libraries(gOLCAOoptc gGaussIntegrals)

add_executable(gOLCAOwave $<TARGET_OBJECTS:gCommon> ../basis.f90 
   ../hdf5PSCFBand.F90 ../hdf5MainEval.f90 ../hdf5MainEvec.F90 
   ../hdf5MainPotRho.f90 ../hdf5Main.f90 ../matrixSubs.F90 ../secularEqn.F90 
   ../populate.F90 ../openDX.f90 ../wave.F90 ../OLCAO_Wave.f90 )

# At this point we need to add the -DGAMMA preprocessor definition to the
#   compile command for these gamma k-point executables. This addition applies
#   only to the current directory.
add_definitions(-DGAMMA)

install(TARGETS gOLCAOsetup gOLCAOmain gOLCAOintg gOLCAOband gOLCAOdos
   gOLCAObond gOLCAOoptc gOLCAOwave DESTINATION bin)

