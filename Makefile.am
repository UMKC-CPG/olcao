bin_PROGRAMS = makekpoints atomSCF gaussFit contract applySpaceGroup OLCAOkkc \
   OLCAOrhoV gOLCAOsetup gOLCAOmain gOLCAOintg gOLCAOband gOLCAOdos \
   gOLCAObond gOLCAOoptc gOLCAOwave OLCAOsetup OLCAOmain OLCAOintg OLCAOband \
   OLCAOdos OLCAObond OLCAOoptc OLCAOwave

clean-local:
	rm -f *.mod

real:	gOLCAOsetup gOLCAOmain gOLCAOintg gOLCAOband gOLCAObond gOLCAOoptc \
      gOLCAOdos gOLCAOwave

complex:	OLCAOsetup OLCAOmain OLCAOintg OLCAOband OLCAObond OLCAOoptc \
         OLCAOdos OLCAOwave

aux: makekpoints atomSCF gaussFit contract applySpaceGroup OLCAOkkc OLCAOrhoV

all: real complex aux

gOLCAOsetup_CPPFLAGS = -DGAMMA
gOLCAOmain_CPPFLAGS  = -DGAMMA
gOLCAOintg_CPPFLAGS  = -DGAMMA
gOLCAOband_CPPFLAGS  = -DGAMMA
gOLCAObond_CPPFLAGS  = -DGAMMA
gOLCAOoptc_CPPFLAGS  = -DGAMMA
gOLCAOwave_CPPFLAGS  = -DGAMMA
gOLCAOdos_CPPFLAGS   = -DGAMMA

makekpoints_SOURCES = src/kinds.f90 src/constants.f90 \
   src/makeKPoints/makekpoints.f90

atomSCF_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/radialGrid.f90 src/atomSCF/dataModules.f90 src/atomSCF/splines.f90 \
   src/atomSCF/readSCFInput.f90 src/atomSCF/initRho.f90 src/atomSCF/mixer.f90 \
   src/atomSCF/totalEnergy.f90 src/atomSCF/implicitSCFInput.f90 \
   src/atomSCF/computeIonicPot.f90 src/atomSCF/exchCorrSubs.f90 \
   src/atomSCF/computeElectPot.f90 src/atomSCF/computeOrbitals.f90 \
   src/atomSCF/checkConvg.f90 src/atomSCF/atomSCF.f90

gaussFit_SOURCES = src/kinds.f90 src/readData.f90 src/gaussFit/dataModules.f90 \
   src/gaussFit/gaussFit.f90

contract_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/writeData.f90 src/radialGrid.f90 src/contract/dataModules.f90 \
   src/contract/matrixElements.f90 src/contract/chargeDensity.f90 \
   src/contract/matrixDiag.f90 src/contract/normalization.f90 \
   src/contract/printResults.f90 src/contract/parseInput.f90 \
   src/contract/contract.f90

applySpaceGroup_SOURCES = src/kinds.f90 src/constants.f90 \
   src/applySpaceGroup/applySpaceGroup.f90

OLCAOkkc_SOURCES = src/kinds.f90 src/constants.f90 \
   src/olcao/auxiliary/OLCAOkkc.f90

OLCAOrhoV_SOURCES = src/kinds.f90 src/constants.f90 \
   src/olcao/auxiliary/OLCAOrhoV.f90

OLCAOsetup_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/basis.f90 src/olcao/kpoints.f90 src/olcao/gaussRelations.f90 \
   src/olcao/gaussIntegrals.f90 src/olcao/hdf5SetupIntg.F90 \
   src/olcao/hdf5SetupElec.f90 src/olcao/hdf5SetupExco.f90 \
   src/olcao/hdf5Setup.f90 src/olcao/exchCorr.f90 src/olcao/elecStat.f90 \
   src/olcao/intgOrtho.F90 src/olcao/intgSaving.F90 src/olcao/integralsSCF.F90 \
   src/olcao/coreCharge.f90 src/olcao/input.f90 src/olcao/OLCAO_Setup.f90 \
   src/olcao/OLCAO_SetupWrap.f90
OLCAOmain_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/kpoints.f90 src/olcao/hdf5SetupIntg.F90 \
   src/olcao/hdf5SetupElec.f90 src/olcao/hdf5SetupExco.f90 \
   src/olcao/hdf5Setup.f90 src/olcao/hdf5MainEval.f90 \
   src/olcao/hdf5MainEvec.F90 src/olcao/hdf5Main.f90 \
   src/olcao/hdf5PSCFBand.F90 src/olcao/exchCorr.f90 src/olcao/elecStat.f90 \
   src/olcao/matrixSubs.F90 src/olcao/input.f90 src/olcao/secularEqn.F90 \
   src/olcao/populate.F90 src/olcao/valeCharge.F90 \
   src/olcao/potentialUpdate.f90 src/olcao/dos.F90 src/olcao/bond.F90 \
   src/olcao/bond3C.F90 src/olcao/OLCAO_Main.f90 src/olcao/OLCAO_MainWrap.f90
OLCAOintg_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 \
   src/olcao/timeStamps.f90 src/olcao/commandLine.f90 src/olcao/potSites.f90 \
   src/olcao/potTypes.f90 \
   src/olcao/atomicTypes.f90 src/olcao/atomicSites.f90 src/olcao/potential.f90 \
   src/olcao/lattice.f90 src/olcao/basis.f90 src/olcao/kpoints.f90 \
   src/olcao/gaussIntegrals.f90 src/olcao/hdf5PSCFIntg.f90 \
   src/olcao/hdf5SetupIntg.F90 src/olcao/hdf5SetupElec.f90 \
   src/olcao/hdf5SetupExco.f90 src/olcao/hdf5Setup.f90 src/olcao/exchCorr.f90 \
   src/olcao/intgOrtho.F90 src/olcao/intgSaving.F90 \
   src/olcao/integralsPSCF.F90 src/olcao/input.f90 src/olcao/OLCAO_Intg.f90 \
   src/olcao/OLCAO_IntgWrap.f90
OLCAOband_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/basis.f90 src/olcao/kpoints.f90 src/olcao/gaussIntegrals.f90 \
   src/olcao/hdf5PSCFBand.F90 src/olcao/hdf5PSCFIntg.f90 \
   src/olcao/hdf5SetupIntg.F90 src/olcao/hdf5SetupElec.f90 \
   src/olcao/hdf5SetupExco.f90 src/olcao/hdf5Setup.f90 \
   src/olcao/hdf5MainEval.f90 src/olcao/hdf5MainEvec.F90 \
   src/olcao/hdf5Main.f90 src/olcao/exchCorr.f90 src/olcao/intgOrtho.F90 \
   src/olcao/intgSaving.F90 src/olcao/integralsPSCF.F90 \
   src/olcao/matrixSubs.F90 src/olcao/input.f90 src/olcao/secularEqn.F90 \
   src/olcao/populate.F90 src/olcao/OLCAO_Band.F90 \
   src/olcao/OLCAO_BandWrap.f90
OLCAOdos_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/kpoints.f90 src/olcao/hdf5PSCFBand.F90 \
   src/olcao/hdf5SetupIntg.F90 src/olcao/hdf5SetupElec.f90 \
   src/olcao/hdf5SetupExco.f90 src/olcao/hdf5Setup.f90 \
   src/olcao/hdf5MainEval.f90 src/olcao/hdf5MainEvec.F90 \
   src/olcao/hdf5Main.f90 src/olcao/exchCorr.f90 src/olcao/matrixSubs.F90 \
   src/olcao/input.f90 src/olcao/secularEqn.F90 src/olcao/populate.F90 \
   src/olcao/dos.F90 src/olcao/OLCAO_DOS.f90
OLCAObond_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/kpoints.f90 src/olcao/hdf5PSCFBand.F90 \
   src/olcao/hdf5SetupIntg.F90 src/olcao/hdf5SetupElec.f90 \
   src/olcao/hdf5SetupExco.f90 src/olcao/hdf5Setup.f90 \
   src/olcao/hdf5MainEval.f90 src/olcao/hdf5MainEvec.F90 \
   src/olcao/hdf5Main.f90 src/olcao/exchCorr.f90 src/olcao/matrixSubs.F90 \
   src/olcao/input.f90 src/olcao/secularEqn.F90 src/olcao/populate.F90 \
   src/olcao/bond.F90 src/olcao/bond3C.F90 src/olcao/OLCAO_Bond.f90
OLCAOoptc_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/basis.f90 src/olcao/kpoints.f90 src/olcao/gaussIntegrals.f90 \
   src/olcao/hdf5PSCFBand.F90 src/olcao/hdf5PSCFIntg.f90 \
   src/olcao/hdf5SetupIntg.F90 src/olcao/hdf5SetupElec.f90 \
   src/olcao/hdf5SetupExco.f90 src/olcao/hdf5Setup.f90 \
   src/olcao/hdf5MainEval.f90 src/olcao/hdf5MainEvec.F90 \
   src/olcao/hdf5Main.f90 src/olcao/exchCorr.f90 src/olcao/intgOrtho.F90 \
   src/olcao/intgSaving.F90 src/olcao/integralsPSCF.F90 \
   src/olcao/matrixSubs.F90 src/olcao/input.f90 src/olcao/secularEqn.F90 \
   src/olcao/populate.F90 src/olcao/optc.F90 src/olcao/optcPrint.F90 \
   src/olcao/OLCAO_Optc.f90 src/olcao/OLCAO_OptcWrap.f90
OLCAOwave_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/elementData.f90 src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 \
   src/olcao/sortSubs.f90 src/olcao/interfaces.F90 src/olcao/timeStamps.f90 \
   src/olcao/commandLine.f90 src/olcao/potSites.f90 src/olcao/potTypes.f90 \
   src/olcao/atomicTypes.f90 src/olcao/atomicSites.f90 src/olcao/potential.f90 \
   src/olcao/lattice.f90 src/olcao/basis.f90 src/olcao/kpoints.f90 \
   src/olcao/hdf5PSCFBand.F90 src/olcao/hdf5SetupIntg.F90 \
   src/olcao/hdf5SetupElec.f90 src/olcao/hdf5SetupExco.f90 \
   src/olcao/hdf5Setup.f90 src/olcao/hdf5MainEval.f90 \
   src/olcao/hdf5MainEvec.F90 src/olcao/hdf5Main.f90 src/olcao/exchCorr.f90 \
   src/olcao/matrixSubs.F90 src/olcao/input.f90 src/olcao/secularEqn.F90 \
   src/olcao/populate.F90 src/olcao/openDX.f90 src/olcao/wave.F90 \
   src/olcao/OLCAO_Wave.f90

gOLCAOsetup_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/basis.f90 src/olcao/kpoints.f90 src/olcao/gaussRelations.f90 \
   src/olcao/gaussIntegrals.f90 src/olcao/hdf5SetupIntg.F90 \
   src/olcao/hdf5SetupElec.f90 src/olcao/hdf5SetupExco.f90 \
   src/olcao/hdf5Setup.f90 src/olcao/exchCorr.f90 src/olcao/elecStat.f90 \
   src/olcao/intgOrtho.F90 src/olcao/intgSaving.F90 src/olcao/integralsSCF.F90 \
   src/olcao/coreCharge.f90 src/olcao/input.f90 src/olcao/OLCAO_Setup.f90 \
   src/olcao/OLCAO_SetupWrap.f90
gOLCAOmain_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/kpoints.f90 src/olcao/hdf5SetupIntg.F90 \
   src/olcao/hdf5SetupElec.f90 src/olcao/hdf5SetupExco.f90 \
   src/olcao/hdf5Setup.f90 src/olcao/hdf5MainEval.f90 \
   src/olcao/hdf5MainEvec.F90 src/olcao/hdf5Main.f90 \
   src/olcao/hdf5PSCFBand.F90 src/olcao/exchCorr.f90 src/olcao/elecStat.f90 \
   src/olcao/matrixSubs.F90 src/olcao/input.f90 src/olcao/secularEqn.F90 \
   src/olcao/populate.F90 src/olcao/valeCharge.F90 \
   src/olcao/potentialUpdate.f90 src/olcao/dos.F90 src/olcao/bond.F90 \
   src/olcao/bond3C.F90 src/olcao/OLCAO_Main.f90 src/olcao/OLCAO_MainWrap.f90
gOLCAOintg_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/basis.f90 src/olcao/kpoints.f90 src/olcao/gaussIntegrals.f90 \
   src/olcao/hdf5PSCFIntg.f90 src/olcao/hdf5SetupIntg.F90 \
   src/olcao/hdf5SetupElec.f90 src/olcao/hdf5SetupExco.f90 \
   src/olcao/hdf5Setup.f90 src/olcao/exchCorr.f90 src/olcao/intgOrtho.F90 \
   src/olcao/intgSaving.F90 src/olcao/integralsPSCF.F90 src/olcao/input.f90 \
   src/olcao/OLCAO_Intg.f90 src/olcao/OLCAO_IntgWrap.f90
gOLCAOband_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/basis.f90 src/olcao/kpoints.f90 src/olcao/gaussIntegrals.f90 \
   src/olcao/hdf5PSCFBand.F90 src/olcao/hdf5PSCFIntg.f90 \
   src/olcao/hdf5SetupIntg.F90 src/olcao/hdf5SetupElec.f90 \
   src/olcao/hdf5SetupExco.f90 src/olcao/hdf5Setup.f90 \
   src/olcao/hdf5MainEval.f90 src/olcao/hdf5MainEvec.F90 \
   src/olcao/hdf5Main.f90 src/olcao/exchCorr.f90 src/olcao/intgOrtho.F90 \
   src/olcao/intgSaving.F90 src/olcao/integralsPSCF.F90 \
   src/olcao/matrixSubs.F90 src/olcao/input.f90 src/olcao/secularEqn.F90 \
   src/olcao/populate.F90 src/olcao/OLCAO_Band.F90 src/olcao/OLCAO_BandWrap.f90
gOLCAOdos_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/kpoints.f90 src/olcao/hdf5PSCFBand.F90 \
   src/olcao/hdf5SetupIntg.F90 src/olcao/hdf5SetupElec.f90 \
   src/olcao/hdf5SetupExco.f90 src/olcao/hdf5Setup.f90 \
   src/olcao/hdf5MainEval.f90 src/olcao/hdf5MainEvec.F90 \
   src/olcao/hdf5Main.f90 src/olcao/exchCorr.f90 src/olcao/matrixSubs.F90 \
   src/olcao/input.f90 src/olcao/secularEqn.F90 src/olcao/populate.F90 \
   src/olcao/dos.F90 src/olcao/OLCAO_DOS.f90
gOLCAObond_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/kpoints.f90 src/olcao/hdf5PSCFBand.F90 \
   src/olcao/hdf5SetupIntg.F90 src/olcao/hdf5SetupElec.f90 \
   src/olcao/hdf5SetupExco.f90 src/olcao/hdf5Setup.f90 \
   src/olcao/hdf5MainEval.f90 src/olcao/hdf5MainEvec.F90 \
   src/olcao/hdf5Main.f90 src/olcao/exchCorr.f90 src/olcao/matrixSubs.F90 \
   src/olcao/input.f90 src/olcao/secularEqn.F90 src/olcao/populate.F90 \
   src/olcao/bond.F90 src/olcao/bond3C.F90 src/olcao/OLCAO_Bond.f90
gOLCAOoptc_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 src/olcao/sortSubs.f90 \
   src/olcao/interfaces.F90 src/olcao/timeStamps.f90 src/olcao/commandLine.f90 \
   src/olcao/potSites.f90 src/olcao/potTypes.f90 src/olcao/atomicTypes.f90 \
   src/olcao/atomicSites.f90 src/olcao/potential.f90 src/olcao/lattice.f90 \
   src/olcao/basis.f90 src/olcao/kpoints.f90 src/olcao/gaussIntegrals.f90 \
   src/olcao/hdf5PSCFBand.F90 src/olcao/hdf5PSCFIntg.f90 \
   src/olcao/hdf5SetupIntg.F90 src/olcao/hdf5SetupElec.f90 \
   src/olcao/hdf5SetupExco.f90 src/olcao/hdf5Setup.f90 \
   src/olcao/hdf5MainEval.f90 src/olcao/hdf5MainEvec.F90 \
   src/olcao/hdf5Main.f90 src/olcao/exchCorr.f90 src/olcao/intgOrtho.F90 \
   src/olcao/intgSaving.F90 src/olcao/integralsPSCF.F90 \
   src/olcao/matrixSubs.F90 src/olcao/input.f90 src/olcao/secularEqn.F90 \
   src/olcao/populate.F90 \
   src/olcao/optc.F90 src/olcao/optcPrint.F90 src/olcao/OLCAO_Optc.f90 \
   src/olcao/OLCAO_OptcWrap.f90
gOLCAOwave_SOURCES = src/kinds.f90 src/constants.f90 src/readData.f90 \
   src/elementData.f90 src/olcao/stringSubs.f90 src/olcao/mathSubs.f90 \
   src/olcao/sortSubs.f90 src/olcao/interfaces.F90 src/olcao/timeStamps.f90 \
   src/olcao/commandLine.f90 src/olcao/potSites.f90 src/olcao/potTypes.f90 \
   src/olcao/atomicTypes.f90 src/olcao/atomicSites.f90 src/olcao/potential.f90 \
   src/olcao/lattice.f90 src/olcao/basis.f90 src/olcao/kpoints.f90 \
   src/olcao/hdf5PSCFBand.F90 src/olcao/hdf5SetupIntg.F90 \
   src/olcao/hdf5SetupElec.f90 src/olcao/hdf5SetupExco.f90 \
   src/olcao/hdf5Setup.f90 src/olcao/hdf5MainEval.f90 \
   src/olcao/hdf5MainEvec.F90 src/olcao/hdf5Main.f90 src/olcao/exchCorr.f90 \
   src/olcao/matrixSubs.F90 src/olcao/input.f90 src/olcao/secularEqn.F90 \
   src/olcao/populate.F90 src/olcao/openDX.f90 src/olcao/wave.F90 \
   src/olcao/OLCAO_Wave.f90

dist_data_DATA = src/data/atomicBDB.tgz src/data/atomicPDB.tgz \
   src/data/spaceDB.tgz src/data/sybdDB.tgz src/data/elements.dat \
   src/data/xc_code.dat

dist_bin_SCRIPTS = src/scripts/addAtoms src/scripts/appendCol \
   src/scripts/atomPosDX \
   src/scripts/batchChange src/scripts/batchClean src/scripts/batchSub \
   src/scripts/bondAnalysis \
   src/scripts/cif2olcao src/scripts/collectSpectra src/scripts/copyInputs \
   src/scripts/dump2skl \
   src/scripts/gaussOverlapRT \
   src/scripts/insert3cbo src/scripts/intGen \
   src/scripts/lmp2skl src/scripts/makeBDB src/scripts/makeBOND \
   src/scripts/makeFittedRhoV src/scripts/makeinput src/scripts/makePDOS \
   src/scripts/makePotDB src/scripts/makeSGDB src/scripts/makeSYBD \
   src/scripts/massGather src/scripts/massPrep src/scripts/massSubmit \
   src/scripts/mkdirhier src/scripts/modStruct \
   src/scripts/normalize \
   src/scripts/olcao src/scripts/olcaoObject.py src/scripts/olcao.py \
   src/scripts/olcaoStructures.py \
   src/scripts/pdb2skl src/scripts/plotgraph src/scripts/polcao \
   src/scripts/rpdf src/scripts/runIsoAtoms \
   src/scripts/scanSpectra src/scripts/selectSYBD src/scripts/skl2vasp \
   src/scripts/spectra2ODX \
   src/scripts/unpackOLCAODB \
   src/scripts/vasp2skl \
   src/scripts/xcMesh src/scripts/xyz2skl

dist_perllib_SCRIPTS = src/scripts/ElementData.pm \
   src/scripts/StructureControl.pm
