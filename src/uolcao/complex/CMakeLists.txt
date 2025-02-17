# CMake file for unified OLCAO application: complex uOLCAO

# Document the current location
message("At the complex uOLCAO application.")

# Library for gaussIntegrals.f90
add_library(uGaussIntegrals STATIC ../../kinds.f90 ../../constants.f90
   ../gaussIntegrals.f90)

add_executable(uOLCAO ../../kinds.f90 ../../constants.f90
   ../../elementData.f90 ../../readDataSerial.f90 ../stringSubs.f90
   ../mathSubs.f90 ../sortSubs.f90 ../interfaces.F90 ../timeStamps.f90
   ../commandLine.f90 ../potSites.f90 ../potTypes.f90 ../atomicTypes.f90
   ../atomicSites.f90 ../potential.f90 ../lattice.f90 ../kpoints.f90
   ../gaussRelations.f90 ../hdf5PSCFIntg.F90 ../hdf5PSCFEigVal.f90
   ../hdf5PSCFEigVec.F90 ../hdf5PSCFField.F90 ../hdf5PSCF.F90
   ../hdf5SCFIntg.F90 ../hdf5SCFElec.f90 ../hdf5SCFExco.f90
   ../hdf5SCFEigVal.f90 ../hdf5SCFEigVec.F90 ../hdf5SCFPotRho.f90
   ../hdf5Field.F90 ../hdf5SCF.F90 ../exchCorr.f90 ../input.f90 ../basis.f90
   ../elecStat.f90 ../intgOrtho.F90 ../matrixSubs.F90 ../secularEqn.F90
   ../populate.F90 ../intgSaving.F90 ../integrals.F90 ../integrals3Terms.F90
   ../coreCharge.f90 ../forces.F90 ../valeCharge.F90 ../potentialUpdate.F90
   ../dos.F90 ../bond.F90 ../bond3C.F90 ../optc.F90 ../optcPrint.F90
   ../loen.f90 ../openDX.f90 ../xdmf_vtk.f90 ../field.F90 ../OLCAO.F90
   ../OLCAOWrap.f90 )
target_link_libraries(uOLCAO uGaussIntegrals)

install(TARGETS uOLCAO DESTINATION bin)
