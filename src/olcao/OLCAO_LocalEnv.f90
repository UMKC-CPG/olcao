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
   ! The largest input number that will not overflow the 8 byte integer
   !   size is 20. This is because 20! = 2,432,902,008,176,640,000 = ~2.432e18
   !   21! = 51,090,942,171,709,440,000 = ~ 5.109e19 which is greater than the
   !   maximum 8 byte signed integer of 9,223,372,036,854,775,807 = ~9.223e18.
   ! An alternative would be to use the gfortran compiler which permits a 16
   !   byte integer with a maximum value of:
   !   170,141,183,460,469,231,731,687,303,715,884,105,727 = ~170.1e36
   ! Another alternative would be to use floating point numbers for the
   !   factorial and permit approximate values. It isn't clear what impact
   !   that would have on the results though.
   if (twoj1 + twoj2 + 1 > 20) then
      write (20,fmt="(a,a,i5,a)") &
            & "Maximum default factorial argument is 20. ", &
            & "You requested: ",twoj1 + twoj2 + 1, " from twoj1 + twoj2 + 1:"
      write (20,fmt="(a,i5)") "twoj1 = ", twoj1
      write (20,fmt="(a,i5)") "twoj2 = ", twoj2
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
