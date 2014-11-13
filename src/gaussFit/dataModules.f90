module ExecutionData

   ! Import the necessary modules.
   use O_Kinds

   integer :: line              ! Current line number.
   integer :: fitByAlphas       ! Flag to fit data by varying alphas or coeffs.
   integer :: iter              ! Number of iterations.
   integer :: match             ! Flag for comparing steps of the iterations.
   integer :: patternMoveCheck  ! Flag for choosing the mode of the next step.
   integer :: stopCount         ! Counts how many terms have converged.
   integer :: lQN               ! l_QN of the major or associate major comp.
   integer :: radialWeight      ! Radial weighting factor for RMS calc.

end module ExecutionData

module GaussianData

   ! Import the necessary modules.
   use O_Kinds

   integer :: numTerms ! Number of gaussian functions to use in the fitting.
   real (kind=double) :: minTerm ! Minimum exponential alpha.
   real (kind=double) :: maxTerm ! Maximum exponential alpha.
   real (kind=double), allocatable, dimension (:)   :: alphas ! List of the
         ! exponential alphas.  Dim=numTerms
   real (kind=double), allocatable, dimension (:)   :: As ! List of the
         ! coefficients to the gaussian functions.  Dim=numTerms
   real (kind=double), allocatable, dimension (:)   :: gaussFnSum ! Summation
         ! of all numerically evaluated gaussian functions.  Dim=numPoints
   real (kind=double), allocatable, dimension (:)   :: gaussFnSumdr ! Summation
         ! of the derivative of all num eval gaussian fns.  Dim=numPoints
   real (kind=double), allocatable, dimension (:,:) :: gaussFn ! Numerical
         ! evaluation of each gaussian function's exponential component only.
         ! Dim=numPoints,numTerms

end module GaussianData

module GivenData

   ! Import the necessary modules.
   use O_Kinds

   integer :: numPoints ! Number of numerical data points.
   real (kind=double), allocatable, dimension (:) :: radialValue ! Radial
         ! values at which the numerical functions are evaluated. Dim=numPoints
   real (kind=double), allocatable, dimension (:) :: radialValueSqrd ! Radial
         ! values from above, squard for faster computation.
   real (kind=double), allocatable, dimension (:) :: dataFn ! Data values to be
         ! fit with sum of gaussians.  Dim=numPoints
   real (kind=double), allocatable, dimension (:) :: dataFndr ! Derivative of
         ! the data values to be fit.  Dim=numPoints
end module GivenData


module FittingData

   ! Import the necessary modules.
   use O_Kinds

   real (kind=double) :: bestCompRMS ! The current best RMS value.
   real (kind=double) :: testCompRMS ! RMS from test move of coefficients.
   real (kind=double) :: bestRMSdr ! The current best RMS of the derivative.
   real (kind=double) :: testRMSdr ! The RMS of the derivative from test move.
   real (kind=double) :: backupRMS ! The RMS value before a pattern move.  ???
   real (kind=double) :: savedCoeff ! Used for exploratory moves of each term.
   real (kind=double), allocatable, dimension(:) :: weight ! Weighting factors
         ! for integration.  Dim=numPoints
   real (kind=double), allocatable, dimension(:) :: weightSqrd ! Weighting
         ! factors from above squard for faster computation.
   real (kind=double), allocatable, dimension(:) :: step ! Size of the change
         ! to apply to each term to be modified when searching for the best
         ! fit.  Dim=numTerms
   real (kind=double), allocatable, dimension(:) :: stopStep ! Minimum size of
         ! the change to apply to each term to be modified when searching for
         ! the best fit.  Dim=numTerms
   real (kind=double), allocatable, dimension(:) :: coeffs ! The set of
         ! coefficients being changed (contains either the alphas or the As).
         ! Dim=numTerms
   real (kind=double), allocatable, dimension(:) :: bestCoeffs ! The current
         ! best choice for the coefficients being changed.  Dim=numTerms

end module FittingData
