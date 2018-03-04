# CMake file for OLCAO application: complex OLCAO

# Document the current location
message("At the complex OLCAO application.")

# Library containing just Kinds and Constants so that the GaussIntegral library
# only depends on a library source file that almost never change.
add_library(kindsConsts STATIC ../../kinds.f90 ../../constants.f90
   ../../elementData.f90)

# Library for gaussIntegrals.f90
add_library(GaussIntegrals STATIC ../gaussIntegrals.f90)
target_link_libraries(GaussIntegrals kindsConsts)

# Library containing common source files that all executables need
add_library(Common OBJECT ../../readData.f90 ../stringSubs.f90 ../mathSubs.f90 
   ../sortSubs.f90 ../interfaces.F90 ../timeStamps.f90 ../commandLine.f90 
   ../potSites.f90 ../potTypes.f90 ../atomicTypes.f90 ../atomicSites.f90 
   ../potential.f90 ../lattice.f90 ../kpoints.f90 ../gaussRelations.f90 
   ../hdf5SetupIntg.F90 ../hdf5SetupElec.f90 ../hdf5SetupExco.f90
   ../hdf5Setup.f90 ../exchCorr.f90 ../input.f90)

add_executable(OLCAOsetup $<TARGET_OBJECTS:Common> ../basis.f90 
   ../elecStat.f90 ../intgOrtho.F90 ../intgSaving.F90 ../integralsSCF.F90 
   ../coreCharge.f90 ../OLCAO_Setup.f90 ../OLCAO_SetupWrap.f90 )
target_link_libraries(OLCAOsetup GaussIntegrals)

add_executable(OLCAOmain $<TARGET_OBJECTS:Common> ../hdf5MainEval.f90 
   ../hdf5MainEvec.F90 ../hdf5MainPotRho.f90 ../hdf5Main.f90 
   ../hdf5PSCFBand.F90 ../elecStat.f90 ../matrixSubs.F90 ../secularEqn.F90 
   ../populate.F90 ../valeCharge.F90 ../potentialUpdate.F90 ../dos.F90 
   ../bond.F90 ../bond3C.F90 ../OLCAO_Main.f90 ../OLCAO_MainWrap.f90 )
target_link_libraries(OLCAOmain kindsConsts)

add_executable(OLCAOintg $<TARGET_OBJECTS:Common> ../basis.f90 
   ../hdf5PSCFIntg.f90 ../intgOrtho.F90 ../intgSaving.F90 ../integralsPSCF.F90 
   ../OLCAO_Intg.f90 ../OLCAO_IntgWrap.f90 )
target_link_libraries(OLCAOintg GaussIntegrals)

add_executable(OLCAOband $<TARGET_OBJECTS:Common> ../basis.f90 
   ../hdf5PSCFBand.F90 ../hdf5PSCFIntg.f90 ../hdf5MainEval.f90 
   ../hdf5MainEvec.F90 ../hdf5MainPotRho.f90 ../hdf5Main.f90 ../intgOrtho.F90 
   ../intgSaving.F90 ../integralsPSCF.F90 ../matrixSubs.F90 ../secularEqn.F90 
   ../populate.F90 ../OLCAO_Band.F90 ../OLCAO_BandWrap.f90 )
target_link_libraries(OLCAOband GaussIntegrals)

add_executable(OLCAOdos $<TARGET_OBJECTS:Common> ../hdf5PSCFBand.F90 
   ../hdf5MainEval.f90 ../hdf5MainEvec.F90 ../hdf5MainPotRho.f90 
   ../hdf5Main.f90 ../matrixSubs.F90 ../secularEqn.F90 ../populate.F90 
   ../dos.F90 ../OLCAO_DOS.f90 )
target_link_libraries(OLCAOdos kindsConsts )

add_executable(OLCAObond $<TARGET_OBJECTS:Common> ../hdf5PSCFBand.F90 
   ../hdf5MainEval.f90 ../hdf5MainEvec.F90 ../hdf5MainPotRho.f90 
   ../hdf5Main.f90 ../matrixSubs.F90 ../secularEqn.F90 ../populate.F90 
   ../bond.F90 ../bond3C.F90 ../OLCAO_Bond.f90 )
target_link_libraries(OLCAObond kindsConsts)

add_executable(OLCAOoptc $<TARGET_OBJECTS:Common> ../basis.f90 
   ../hdf5PSCFBand.F90 ../hdf5PSCFIntg.f90 ../hdf5MainEval.f90 
   ../hdf5MainEvec.F90 ../hdf5MainPotRho.f90 ../hdf5Main.f90 ../intgOrtho.F90 
   ../intgSaving.F90 ../integralsPSCF.F90 ../matrixSubs.F90 ../secularEqn.F90 
   ../populate.F90 ../optc.F90 ../optcPrint.F90 ../OLCAO_Optc.f90 
   ../OLCAO_OptcWrap.f90 )
target_link_libraries(OLCAOoptc GaussIntegrals)

add_executable(OLCAOwave $<TARGET_OBJECTS:Common> ../basis.f90 
   ../hdf5PSCFBand.F90 ../hdf5MainEval.f90 ../hdf5MainEvec.F90 
   ../hdf5MainPotRho.f90 ../hdf5Main.f90 ../matrixSubs.F90 ../secularEqn.F90 
   ../populate.F90 ../openDX.f90 ../wave.F90 ../OLCAO_Wave.f90 )
target_link_libraries(OLCAOwave kindsConsts )

install(TARGETS OLCAOsetup OLCAOmain OLCAOintg OLCAOband OLCAOdos
   OLCAObond OLCAOoptc OLCAOwave DESTINATION bin)
