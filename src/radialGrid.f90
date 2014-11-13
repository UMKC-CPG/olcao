module O_RadialGrid

   ! Import the necessary modules.
   use O_Kinds

   public

   integer :: numRadialPoints
   real (kind=double) :: radialMaxDist ! This is the maximum distance from the
         ! atomic center that we will compute any values (charge, waveFn, etc.)
   real (kind=double) :: aaWhatever    ! This will modify the number of points
         ! on the radial mesh.  Increasing this value will increase the number
         ! of points with the vast majority of them near the nucleus.
   real (kind=double) :: bbWhatever    ! This will also modify the number of
         ! points on the radial mesh.  Increasing thie value will increase the
         ! number of points with the vast majority of them near the nucleus.

   ! The formal origin for aaWhatever and bbWhatever was left out of the
   !   documentation of the original code and I do not want to spend the time
   !   now to try and rediscover it.  This is why exhaustive and explicit
   !   documentation is so helpful and important.

   ! Points along the radial direction with a weighting for many more points
   !   (i.e. higher resolution) near the nucleus (r=0).  It may be possible
   !   that the values in radialPoints2 are the distances of separation between
   !   points in radialPoints, but I'm not sure and this must be checked.
   real (kind=double), allocatable, dimension (:) :: radialPoints
   real (kind=double), allocatable, dimension (:) :: radialPoints2

   ! Factors for numerical integration of the radial points using a trapizoidal
   !   rule (?!?).
   real (kind=double), allocatable, dimension (:) :: radialIntFactors

contains

subroutine setupRadialGrid (atomicNumber)

   ! Import necessary definitions.
   use O_Kinds

   implicit none

   ! Define passed parameters.
   real (kind=double) :: atomicNumber

   ! Define local variables.
   integer :: i ! Loop index variable.
   real (kind=double) :: expDistFactor
   real (kind=double) :: expDistFactor2

   ! Compute the exponential distance factors.
   expDistFactor  = exp(-aaWhatever)/atomicNumber
   expDistFactor2 = 1.0_double/bbWhatever

   ! Find the number of points that are needed to make the grid with the
   !   requested maxRadialDist.  In the old program the limit for the number
   !   of points was reached when a distance (r) was calculated past the max
   !   distance limit according to r = a*(exp(b*(i-1))-1) where:
   !   a = exp(-aaWhatever)/atomicNumber = expDistFactor
   !   b = 1/bbWhatever = expDistFactor2
   !   i = integers from 1 to large integer (e.g. 1000).
   !   This number is calculated by solving for i.
   numRadialPoints = floor(log(radialMaxDist/expDistFactor+1.0_double) / &
         & expDistFactor2 + 1.0_double)

   ! Allocate space for the storage of the radial distances.
   allocate (radialPoints  (numRadialPoints))
   allocate (radialPoints2 (numRadialPoints))

   ! Compute the values of the radial mesh in two forms.  I'm not sure why
   !   there are two forms and what their difference really is in terms of
   !   rigorous theory.  The number of points near the nucleus is much
   !   greater than the number of points farther away by a significant
   !   margin.  It may be possible that radialPoints2 is actually the
   !   distance between values in radialPoints but this must be checked!
   do i = 1, numRadialPoints
      radialPoints(i) = expDistFactor * (exp(expDistFactor2 * (i-1)) - 1)
      radialPoints2(i) = (radialPoints(i) + expDistFactor) * expDistFactor2
   enddo

   ! Allocate and fill an array used for crude integration over radial
   !   values.
   allocate (radialIntFactors(numRadialPoints))
   radialIntFactors(1:numRadialPoints:2) = 2.0_double
   radialIntFactors(2:numRadialPoints:2) = 4.0_double

end subroutine setupRadialGrid

end module O_RadialGrid
