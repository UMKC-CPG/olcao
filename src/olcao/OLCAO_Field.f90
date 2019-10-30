program fieldDataPlot

   ! The goal of this program is to produce plottable data of the solid state
   !   wave function, data derived form it (e.g. the charge density), and/or
   !   ancillary data such as the potential function. It permits plotting of
   !   specific wave function states according to energy level by numerically
   !   evaluating the analytical wave functions, wave functions squared, etc.
   !   on a defined 3D uniform mesh. The data may be stored in either plain
   !   text OpenDX format or in XDMF+HDF5 format for Paraview.

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_ElementData,     only: initElementData
   use O_Populate,        only: populateStates
   use O_Basis,           only: renormalizeBasis
   use O_Potential,       only: spin, initPotCoeffs
   use O_CommandLine,     only: parseWaveCommandLine
   use O_Input,           only: numStates, parseInput
   use O_Field,           only: computeFieldMesh, cleanUpField
   use O_KPoints,         only: numKPoints, computePhaseFactors
   use O_FieldHDF5,       only: initFieldHDF5, closeFieldHDF5
   use O_PSCFBandHDF5,    only: accessPSCFBandHDF5, closeAccessPSCFBandHDF5
   use O_SecularEquation, only: energyEigenValues, readEnergyEigenValuesBand
   use O_Lattice,         only: initializeLattice, initializeFindVec, &
         & initialize3DMesh

   ! Make sure that there are not accidental variable declarations.
   implicit none


   ! Open the potential file that will be read from in this program.
   open (unit=8,file='fort.8',status='old',form='formatted')


   ! Initialize the logging labels.
   call initOperationLabels


   ! Initialize element data from periodic table of the elements.
   call initElementData


   ! Parse the command line parameters
   call parseWaveCommandLine


   ! Read in the input to initialize all the key data structure variables.
   call parseInput


   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file.  These values can, however, be easily determined from the input
   !   file.
   call getImplicitInfo


   ! Initialize the HDF5 data file that will be used.
   call initFieldHDF5


   ! Access the HDF5 data stored by band.
   call accessPSCFBandHDF5(numStates)


   ! Create real-space and reciprocal-space super lattices out of the primitive
   !   lattice.  These "supercells" must be big enough so as to include all the
   !   points within a sphere bounded by the negligability limit.  Points
   !   outside the sphere are considered negligable.  The '0' indicates that it
   !   is not necessary to consider the reciprocal lattice.
   call initializeLattice (0)


   ! Setup the necessary data structures so that we can easily find the lattice
   !   vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! Compute the kPoint phase factors.
   call computePhaseFactors


   ! Prepare the basis functions for easy access (and angular normalization).
   call renormalizeBasis


   ! Read the provided potential coefficients.
   call initPotCoeffs


   ! Allocate space to store the energy eigen values, and then read them in.
   allocate (energyEigenValues (numStates,numKPoints,spin))
   call readEnergyEigenValuesBand(numStates)


   ! Populate the electron states to find the highest occupied state (Fermi
   !   energ for metals).
   call populateStates


   ! Initialize certain parameters for constructing and traversing the 3D mesh.
   call initialize3DMesh


   ! Compute the requested field values for each mesh point and store in HDF5.
   call computeFieldMesh


   ! Clean up any left over arrays that need to be deallocated.
   call cleanUpField


   ! Close the HDF5 band data file.
   call closeAccessPSCFBandHDF5


   ! Close the HDF5 field data file.
   call closeFieldHDF5


   ! Open a file to signal completion of the program.
   open (unit=2,file='fort.2',status='new')

end program fieldDataPlot

subroutine getImplicitInfo

   ! Import necessary modules.
   use O_TimeStamps
   use O_ExchangeCorrelation, only: makeSampleVectors
   use O_AtomicSites,         only: getAtomicSiteImplicitInfo
   use O_AtomicTypes,         only: getAtomicTypeImplicitInfo
   use O_PotSites,            only: getPotSiteImplicitInfo
   use O_PotTypes,            only: getPotTypeImplicitInfo
   use O_Lattice,             only: getRecipCellVectors
   use O_KPoints,             only: convertKPointsToXYZ
   use O_Potential,           only: initPotStructures

   implicit none

   call timeStampStart(2)

   ! Subroutines need to be called in this order due to data dependencies.
   call makeSampleVectors

   call getAtomicTypeImplicitInfo
   call getAtomicSiteImplicitInfo
   call getPotSiteImplicitInfo
   call getPotTypeImplicitInfo

   call getRecipCellVectors
   call convertKPointsToXYZ

   call initPotStructures

   call timeStampEnd(2)

end subroutine getImplicitInfo
