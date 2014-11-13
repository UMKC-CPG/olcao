module ConvergenceSubs

public

contains

subroutine checkConvergence

   ! Use necessary modules.
   use O_Kinds

   ! Use necessary data modules.
   use ExecutionData
   use PotData

   ! Use necessary external modules.
   use O_RadialGrid

   implicit none

   ! Define local variables.
   integer :: i
   real (kind=double) :: currDeltaPot ! Current amount of potential change.
   real (kind=double) :: maxDeltaPot  ! Most significant potential change.

   ! If the system has already converged then move along.
   if (converged == 1) return

   ! Compute the current level of convergence.

   ! Initialize the record of the most significant potential change.
   maxDeltaPot = 0.0_double

   ! Check each radial point for potential changes.
   do i = 1, numRadialPoints

      ! Spin down potential first.
      currDeltaPot = (hartreePotDn(i) - hartreePotDn2(i)) / &
            & (1.0_double + hartreePotDn(i) + hartreePotUp(i))
      maxDeltaPot = max(abs(currDeltaPot),maxDeltaPot)

      ! Spin up potential second.
      currDeltaPot = (hartreePotUp(i) - hartreePotUp2(i)) / &
            & (1.0_double + hartreePotUp(i) + hartreePotDn(i))
      maxDeltaPot = max(abs(currDeltaPot),maxDeltaPot)
   enddo

   ! Check if the potential has converged to within the requested tolerance.
   if (maxDeltaPot < tolerance) then
      converged = 1
   endif

   ! If the convergence is sufficiently close then we can switch to using the
   !   alternate wave function finding algorithm that requires a good initial
   !   guess for the energy eigen values.  This is a direct integration scheme.
   if ((maxDeltaPot < 0.01_double*iteration) .and. (solverChoice == 1)) then
      solverChoice = 2
   endif

   ! If the current most significat potential change is greater than the last
   !   most significant potential change, then we reduce the mixing factor, but
   !   not below 0.25.
   if (maxDeltaPot > oldMaxDeltaPot) then
      mixingFactor = 0.8_double * mixingFactor
   endif
   if (mixingFactor < 0.25_double) then
      mixingFactor = 0.25_double
   endif

   ! Record the current largest change for comparison in the next iteration.
   oldMaxDeltaPot = maxDeltaPot


end subroutine checkConvergence


subroutine printConvgPot

   ! Use necessary modules
   use O_Kinds

   ! Use necessary data modules.
   use PotData
   use AtomData
   use ChargeDensityData

   ! Use necessary external module.
   use O_RadialGrid

   implicit none

   ! Define local variables.
   integer :: i
   real (kind=double) :: nucExpArg
   real (kind=double) :: ionicPotOverR
   real (kind=double) :: screeningPot
   real (kind=double), allocatable, dimension (:) :: totalPot
   real (kind=double), allocatable, dimension (:) :: fittablePot

   ! Allocate space to hold the different numerical forms of the potential.
   allocate (totalPot(numRadialPoints))
   allocate (fittablePot(numRadialPoints))

   ! Calculate the fittable electronic potential.  Interestingly, the old
   !   program computed the potential in an apparently non-spin polarized way.
   do i = 2, numRadialPoints

      ! Compute the nuclear potential exponential argument as a Gaussian
      !   function of r.
      nucExpArg = nuclearAlpha * radialPoints(i) * radialPoints(i)

      ! Compute the ionic potential as a function of 1/r.  Down spin only
      !   considered.  Also, only the l=1 QN is considered.  Why?
      ionicPotOverR = ionicPotDn(1,i) / radialPoints(i)

      ! Compute the screening potential.
      screeningPot = hartreePotDn(i)

      ! Compute the total potential.
      totalPot(i) = ionicPotOverR + screeningPot

      ! Compute the fittable electronic potential.
      fittablePot(i) = totalPot(i) - ionicPotOverR * exp(-nucExpArg)
   enddo

   ! Apply a factor of 0.5 to switch from Rydberg to a.u.
   totalPot(2:)    = totalPot(2:)    * 0.5_double
   fittablePot(2:) = fittablePot(2:) * 0.5_double

   ! Open the potential and charge density data files for writing.
   open (unit=44,form="formatted",file="numericalPot.dat",status="new")
   open (unit=45,form="formatted",file="numericalRho.dat",status="new")

   ! Prepare the header of the data file to be fitted with gaussians.
   write (44,fmt="(a15)") "COEFFS0_ALPHAS1"
   write (44,*) 1                 ! Flag to fit using alpha modifications.
   write (44,fmt="(a18)") "NUM_GAUSSIAN_TERMS"
   write (44,*) 30                ! Number of gaussians to use.
   write (44,fmt="(a18)") "MIN_GAUSSIAN_ALPHA"
   write (44,*) 0.12              ! Lowest term.
   write (44,fmt="(a18)") "MAX_GAUSSIAN_ALPHA"
   write (44,*) 10000000.0        ! Highest term.
   write (44,fmt="(a20)") "NUM_NUMERICAL_POINTS"
   write (44,*) numRadialPoints-1 ! Number of data points.
   write (44,fmt="(a17)") "RMS_RADIAL_WEIGHT"
   write (44,*) 0                 ! RMS Radial weight exponent.
   write (44,fmt="(a12)") "RADIAL_COEFF"
   write (44,*) 0                 ! Exponent to an implicit 1/r factor.

   ! Print the numerical potential and charge density data.
   do i = 2, numRadialPoints
      write (44,*) radialPoints(i), fittablePot(i), totalPot(i)
      write (45,*) radialPoints(i), rhoDn(i), rhoUp(i), rhoTotal(i)
   enddo


   ! Close the data file.
   close (44)
   close (45)

   ! Deallocate unnecessary space.
   deallocate (totalPot)
   deallocate (fittablePot)

end subroutine printConvgPot


end module ConvergenceSubs
