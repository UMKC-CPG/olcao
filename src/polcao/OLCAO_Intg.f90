module O_Intg

contains

subroutine intgPSCF

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Input,         only: parseInput
   use O_Potential,     only: initPotCoeffs
   use O_Basis,         only: renormalizeBasis
   use O_CommandLine,   only: doMOME, parseIntgCommandLine
   use O_IntegralsPSCF, only: intgAndOrMom
   use O_PSCFIntgHDF5,  only: initPSCFIntgHDF5, closePSCFIntgHDF5, setMOMEStatus
   use O_Lattice,       only: initializeLattice, initializeFindVec
   use MPI

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define global mpi parameters
   integer :: mpiRank
   integer :: mpiSize
   integer :: mpiErr

   ! Open the potential file that will be read from in this program.
   open (unit=8,file='fort.8',status='old',form='formatted')

   ! Initialize the MPI interface
   call MPI_INIT (mpierr)
   call MPI_COMM_RANK (MPI_COMM_WORLD,mpiRank,mpierr)
   call MPI_COMM_SIZE (MPI_COMM_WORLD,mpiSize,mpierr)

   ! Initialize the logging labels.
   call initOperationLabels


   ! Parse the command line parameters
   call parseIntgCommandLine


   ! Read in the input to initialize all the key data structure variables.
   call parseInput


   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file.  These values can, however, be easily determined from the input
   !   file.
   call getImplicitInfo


   ! Create real-space super lattice out of the primitive lattice.  These
   !   "supercells" must be big enough so as to include all the points within
   !   sphere bounded by the negligability limit.  Points outside the sphere
   !   are considered negligable and are ignored.
   call initializeLattice (1)


   ! Setup the necessary data structures so that we can easily find the lattice
   !   vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! Renormalize the basis functions
   call renormalizeBasis

!write (20,*) "Got here 0"
!call flush (20)

   ! Prepare the HDF5 files for the post SCF integrals.
   if (mpiRank == 0) then
      call initPSCFIntgHDF5
   endif

!write (20,*) "Got here 1"
!call flush (20)

   ! Read the provided potential coefficients.
   call initPotCoeffs

!write (20,*) "Got here 2"
!call flush (20)
   ! Compute the integrals and the momentum matrix elements (MME or MOME) (if
   !   requested).  The MME request is given on the command line and the
   !   subroutine gets that flag from that module.  The integrals are also
   !   optional and are controlled by the intgAndOrMom parameter.  Obviously,
   !   in this program the whole point is to compute the integrals.  However,
   !   in other programs (optc) the integrals do not need to be computed and
   !   only the MME need to be computed.  In that case the same subroutine can
   !   be used, except that only the MME are computed, not the overlap and
   !   hamiltonian.
   call intgAndOrMom(1,doMOME)


   ! Record to the HDF5 file an attribute that says whether or not the MOME
   !   were computed.
   if (mpiRank == 0) then
      call setMOMEStatus(doMOME)
   endif


   ! Close the HDF5 PSCF integrals file.
   if (mpiRank == 0) then
      call closePSCFIntgHDF5
   endif


   ! Open a file to signal completion of the program.
   if (mpiRank == 0) then
      open (unit=2,file='fort.2',status='new')
      close (2)
   endif

   ! End the MPI interface
   call MPI_FINALIZE (mpierr)

end subroutine intgPSCF


subroutine getImplicitInfo

   ! Import necessary modules.
   use O_AtomicSites, only: getAtomicSiteImplicitInfo
   use O_AtomicTypes, only: getAtomicTypeImplicitInfo
   use O_PotSites, only: getPotSiteImplicitInfo
   use O_PotTypes, only: getPotTypeImplicitInfo
   use O_Lattice, only: getRecipCellVectors
   use O_Potential, only: initPotStructures
   use O_TimeStamps

   implicit none

   call timeStampStart(2)

   call getAtomicTypeImplicitInfo
   call getAtomicSiteImplicitInfo
   call getPotSiteImplicitInfo
   call getPotTypeImplicitInfo

   ! Used to later obtain the inverse matrix of the real cell vectors.
   call getRecipCellVectors

   call initPotStructures

   call timeStampEnd(2)

end subroutine getImplicitInfo

end module O_Intg
