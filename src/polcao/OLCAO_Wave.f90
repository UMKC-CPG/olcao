program waveFnPlotDX

   ! The goal of this program is to plot specific wave function states by
   !   evaluating the analytical wave functions or wave functions squared on
   !   a defined 3D uniform mesh.

   ! Import the necessary modules.
   use O_Kinds
   use O_Constants
   use O_ElementData
   use O_CommandLine
   use O_TimeStamps
   use O_Input
   use O_KPoints
   use O_Lattice
   use O_Basis
   use O_Populate
   use O_Potential
   use O_SecularEquation
   use O_PSCFBandHDF5
   use O_Wave
   use O_OpenDX

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


   ! Access the HDF5 data stored by band.
   call accessPSCFBandHDF5


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
   call readEnergyEigenValuesBand


   ! Populate the electron states to find the highest occupied state (Fermi
   !   energ for metals).
   call populateStates


   ! Initialize certain parameters for constructing and traversing the 3D mesh.
   call initialize3DMesh

   ! Compute the wave function squared values for each mesh point and print.
   call computeWaveFnMesh

   ! Clean up any left over arrays that need to be deallocated.
   call cleanUpWave

   ! Close the HDF data file.
   call closeAccessPSCFBandHDF5

   ! Open a file to signal completion of the program.
   open (unit=2,file='fort.2',status='new')

end program waveFnPlotDX

subroutine getImplicitInfo

   ! Import necessary modules.
   use O_ExchangeCorrelation
   use O_AtomicSites
   use O_AtomicTypes
   use O_PotSites
   use O_PotTypes
   use O_Lattice
   use O_KPoints
   use O_Potential
   use O_Populate
   use O_TimeStamps

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
