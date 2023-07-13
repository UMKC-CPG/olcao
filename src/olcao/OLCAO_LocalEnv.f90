module O_LocalEnv

contains

subroutine analyzeLocalEnv

   use O_Kinds

   ! Import necessary modules.
   use O_CommandLine, only: parseLoEnCommandLine
   use O_Input,       only: parseInput, loenCode
   use O_Lattice,     only: initializeLattice, initializeFindVec,&
                          & cleanUpLattice
   use O_TimeStamps,  only: initOperationLabels

   ! Import the HDF5 module.
   use HDF5

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Initialize the logging labels.
   call initOperationLabels


   ! Parse the command line parameters
   call parseLoEnCommandLine


   ! Read in the input to initialize all the key data structure variables.
   call parseInput


   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file.  These values can, however, be easily determined from the input
   !   file.
   call getImplicitInfo


   ! Create real-space and reciprocal-space super lattices out of the primitive
   !   lattice.  These "supercells" must be big enough so as to include all the
   !   points within a sphere bounded by the negligability limit.  Points
   !   outside the sphere are considered negligable.
   call initializeLattice (1) ! Zero because we don't need the reciprocal lat.


   ! Setup the necessary data structures so that we can easily find the lattice
   !   vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! Compute the requested local environment metric.
   if (loenCode == 1) then
      call bispec ! Use the bispectrum component method.
   endif

   ! Close the output file
   close (20)

   ! Open a file to signal completion of the program.
   open (unit=2,file='fort.2',status='unknown')

end subroutine analyzeLocalEnv


subroutine bispec

   ! Import necessary modules
   use O_Input, only: twoj1, twoj2
   use O_MathSubs
   use O_LOEN
   use O_TimeStamps,  only: timeStampStart, timeStampEnd

   ! Print the start tag.
   call timeStampStart(28)

   ! Precompute a list of factorial values. The largest possible factorial
   !   that will be requested will be from the calculation of the Clebsch-
   !   Gordan coefficients when we need to compute (j + j1 + j2 + 1)! where
   !   j = j1 + j2 and 2*j1 and 2*j2 are given in the olcao.dat input file.
   ! The factorial argument can be algebraically reduced to:
   !   2*j1 + 2*j2 + 1.
   ! The largest input number that will not overflow the default integer
   !   size is 12. This is because 12! = 479,001,600 while 13! =
   !   6,227,020,800 which is greater than the maximum default signed
   !   integer of 2,147,483,647.
   if (twoj1 + twoj2 + 1 > 12) then
      write (20,fmt="(a,a,i,a)") &
            & "Maximum default factorial argument is 12. ", &
            & "You requested: ",twoj1 + twoj2 + 1,"from twoj1 + twoj2 + 1:"
      write (20,fmt="(a,i)") "twoj1 = ", twoj1
      write (20,fmt="(a,i)") "twoj2 = ", twoj2
      write (20,fmt="(a)") &
            & "Program needs modification to accomodate larger arguments."
      write (20,fmt="(a)") "Stopping"
      stop
   else
      call computeFactorials((twoj1 + twoj2) + 1)
   endif

   ! As it says. Make a list of the neighbors of each atom.
   call createNeighborList

   ! Compute the set of bispectrum components for each atom.
   call computeBispectrumComponent

   ! Print the end tag.
   call timeStampEnd(28)

end subroutine bispec


subroutine getImplicitInfo

   ! Import necessary modules.
   use O_ExchangeCorrelation
   use O_AtomicTypes
   use O_AtomicSites
   use O_PotSites
   use O_PotTypes
   use O_Lattice
   use O_KPoints
   use O_Potential
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

end module O_LocalEnv
