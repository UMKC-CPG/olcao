program gaussfit

! The purpose of this program is to take a given numerically expressed function
!   and fit it with a series of gaussian functions of a specified range.  Then
!   the program will vary the coefficients to achieve the best possible fit.

   ! Import necessary definitions.
   use O_Kinds

   ! Import necessary subroutines.
   use O_ReadDataSubs

   ! Import necessary data modules.
   use ExecutionData
   use GaussianData
   use GivenData
   use FittingData

   implicit none

   real (kind=double) :: dummyValue   ! Dummy value for reading data.
   real (kind=double) :: expFactor ! Used to determine alphas from min,max.
   integer            :: readUnit  ! The unit number of the file from which
                                   ! we are reading.
   integer            :: writeUnit  ! The unit number of the file to which
                                    ! we are writing.

   ! Loop control variables.
   integer :: i

   ! Set the unit for reading to standard input.
   readUnit = 5
   ! Open files for writing.
   open(44,file='gauss.fit',status='new',form='formatted')
   open(20,file='gaussFit.out',status='new',form='formatted')
   writeUnit = 20
   ! Read input parameters and numerical data.
   call readData(readUnit,writeUnit,fitByAlphas,len('COEFFS0_ALPHAS1'),'COEFFS0_ALPHAS1')
   call readData(readUnit,writeUnit,numTerms,len('NUM_GAUSSIAN_TERMS'),'NUM_GAUSSIAN_TERMS')
   call readData(readUnit,writeUnit,minTerm,len('MIN_GAUSSIAN_ALPHA'),'MIN_GAUSSIAN_ALPHA')
   call readData(readUnit,writeUnit,maxTerm,len('MAX_GAUSSIAN_ALPHA'),'MAX_GAUSSIAN_ALPHA')
   call readData(readUnit,writeUnit,numPoints,len('NUM_NUMERICAL_POINTS'),'NUM_NUMERICAL_POINTS')
   call readData(readUnit,writeUnit,radialWeight,len('RMS_RADIAL_WEIGHT'),'RMS_RADIAL_WEIGHT')
   ! Read the exponent for the division of the radial wave function.  The
   !   wave function must be in the form 1/r^-l * f(r),g(r) where l is the
   !   angular quantum number.  Note that although the minor component has a
   !   different angular character from the major component they will be
   !   treated the same here.
   call readData(readUnit,writeUnit,lQN,len('RADIAL_COEFF'),'RADIAL_COEFF')


   ! Initialize parameters based on read in data.


   ! Initialize program variables.
   patternMoveCheck = 0


   ! Allocate space to hold data.
   allocate(radialValue(numPoints))
   allocate(radialValueSqrd(numPoints))
   allocate(dataFn(numPoints))
   allocate(dataFndr(numPoints))
   allocate(gaussFn(numPoints,numTerms))
   allocate(gaussFnSum(numPoints))
   allocate(gaussFnSumdr(numPoints))
   allocate(weight(numPoints))
   allocate(weightSqrd(numPoints))
   allocate(step(numTerms))
   allocate(stopStep(numTerms))
   allocate(coeffs(numTerms))
   allocate(bestCoeffs(numTerms))
   allocate(alphas(numTerms))
   allocate(As(numTerms))


   ! Read the numerical data.
   do i = 1, numPoints
      read (5,*) radialValue(i), dataFn(i), dummyValue
   enddo

   ! Compute the squard radial values.
   radialValueSqrd(:) = radialValue(:)**2



   ! Divide by the radial value ^ l_QN of the major component.
   dataFn(:) = dataFn(:)/(radialValue(:)**real(lQN))



   ! Compute the derivative of the wave function.
!   call computeDerivative(radialValue,dataFn,dataFndr,numPoints)


   ! Compute the gaussian exponential terms (alphas).
   expFactor=(maxTerm/minTerm)**((1.0_double/real(numTerms-1,double)))**1d0
   do i = 1, numTerms
      alphas(i) = minTerm*expFactor**((i-1))**1d0
!      write (44,*) "Exponential alpha ",i," = ",alphas(i)
   enddo

   ! There are two ways that the fitting can be performed.  The first will
   !   modify the gaussian coefficients (A) while the second will modify the
   !   gaussian exponential coefficients (alpha).  The form of the function is:
   !   A * exp(-alpha * r^2).

   ! The first method will be to compute a trial set of coefficients by solving
   !   the linear coefficient problem Ax=b with:
   !     A(ij) = integral(g(r,alpha_i)*g(r,alpha_j)*weight(r))dr; (g=gaussian)
   !     x(i)  = coefficients to be determined.
   !     b(i)  = integral(g(r,alpha_i)*f(r)*weight(r))dr; (f=numerical data)
   !     The weighting factors are for simpson integration.
   !   Then, we will vary the coefficients one at a time and check for
   !   improvement as in the method of Hooke and Jeeves.  (J. of Assn. Comp.
   !   Mach., 8, 212-229, 1961).

   ! The second will be to use the given trial set of exponential alphas and
   !   solve the linear coefficient problem Ax=b with the same definitions as
   !   above.  Then we will vary the exponential alphas one at a time, compute
   !   new coefficients and check for an improvement in the fit.  The method is
   !   as used above.

   ! In both schemes the variables coeff(1:numTerms) and bestCoeff(1:numTerms)
   !   are used to define the values that will be modified.  In the first case
   !   they contain the As and in the second case they contain the alphas.


   ! Compute the weighting factors for trapizoidal integration.
   do i = 1, numPoints - 1
      weight(i) = (radialValue(i+1) - radialValue(i))/2.0_double
   enddo
   weight(numPoints) = 0.0_double

   ! Compute the weights squard.
   weightSqrd(:) = weight(:)**2

   ! Begin computing the trial set of coefficients.
   call gaussCalc

   ! Record the set of coefficients (As or alphas) to be modified.
   if (fitByAlphas == 0) then
      coeffs(:)     = As(:)
      bestCoeffs(:) = As(:)
   else
      coeffs(:)     = alphas(:)
      bestCoeffs(:) = alphas(:)
   endif

   ! Compute the initial step size based on the coefficients.
   step(:) = coeffs(:)/50.0_double

   ! Compute the stopping point when the step size should not be reduced more.
   stopStep(:) = step(:)/100.0_double

   ! Finished computing the trial set of coefficients.




   ! Begin to find optimal coefficients.

   ! Calculate the initial RMS value.
   call getRMS
   bestCompRMS = testCompRMS

   ! Progress through all terms changing the coefficients to improve the result
   !   as determined by the RMS.  Repeat until there is no significant change
   !   to any coefficient between passes.
   iter = 0
   do while (.true.)

      ! Record the fact that this iteration is being done.
      iter = iter + 1
      if (mod(iter,10) == 0) then
         write (6,ADVANCE="NO",FMT="(a1)") "|"
      else
         write (6,ADVANCE="NO",FMT="(a1)") "."
      endif
      if (mod(iter,50) == 0) then
         write (6,*) " ",iter, bestCompRMS
      endif
      call flush (6)

      ! Consider improvements to each term.
      do i = 1, numTerms

         ! Save the value of the coefficient being modified.
         savedCoeff = coeffs(i)

         ! Modify the coefficient a step up.
         coeffs(i) = savedCoeff + step(i)

         ! Compute the RMS.
         call getRMS

         ! If it improved then cycle to the next term.
         if (testCompRMS < bestCompRMS) then
            bestCompRMS = testCompRMS
            cycle
         endif

         ! It did not improve so we now modify the coefficient a step down.
         coeffs(i) = savedCoeff - step(i)

         ! Compute the RMS.
         call getRMS

         ! If it improved then cycle to the next term.
         if (testCompRMS < bestCompRMS) then
            bestCompRMS = testCompRMS
            cycle
         endif

         ! It did not improve so we return to the original coefficient.
         coeffs(i) = savedCoeff
      enddo


      ! At this point we must compare the current set of coefficients with
      !   the previous set of coefficients to determine if they have changed
      !   at all.  If they have not changed then we halve the step size and
      !   try again until we reach the stopping step size.  If they have
      !   changed then we save the current set of coefficients as the best so
      !   far and prepare a new set to try and modify.


      if (patternMoveCheck == 0) then

         ! Check if the current coefficients and best-so-far match.
         match = 1  ! Assume they match
         do i = 1, numTerms
            if (abs(bestCoeffs(i)-coeffs(i)) >= 1.0d-10) then
               match = 0
               exit
            endif
         enddo

         if (match == 0) then

            ! If they do not match, then perform a pattern move.
            call patternMove

            ! Save and update the RMS.
            backupRMS = bestCompRMS
            call getRMS

            ! Record the fact that we made a pattern move.
            patternMoveCheck = 1
         else

            ! If they match, then halve the step size and check for the
            !   stopping condition.  The stopping condition must be met for
            !   every term.
            stopCount = 0
            do i = 1, numTerms
               step(i) = step(i)/2.0_double
               if (abs(step(i)) < abs(stopStep(i))) then
                  stopCount = stopCount + 1
               endif
            enddo

            ! Check if the stopping condition was met.
            if (stopCount == numTerms) then
               call outputResults
               exit
            endif
         endif
      else

         ! If the bestRMS is better than the backup we proceed with another
         !   patern move.  Otherwise go back to the bestCoeff and try again.
         if (bestCompRMS < backupRMS) then

            ! Perform a pattern move.
            call patternMove

            ! Save and update the RMS.
            backupRMS = bestCompRMS
            call getRMS
         else

            ! Revert to the previous coefficients.
            coeffs(:) = bestCoeffs(:)

            ! Revert to the old RMS
            bestCompRMS = backupRMS

            ! Record the fact that we did not do a pattern move.
            patternMoveCheck = 0
         endif
      endif
   enddo

   ! Finalize the tracking with a hard return and the last bestCompRMS
   if (mod(iter,50) /= 0) then
      write (6,*) " ",iter, bestCompRMS
   endif
   call flush (6)

   
end program gaussfit


subroutine gaussCalc

   ! Import necessary definitions.
   use O_Kinds

   ! Import necessary data modules.
   use GaussianData
   use GivenData
   use FittingData

   implicit none

   ! Define local variables.
   integer :: i,j,k ! Loop index variables.
   integer :: info  ! Used for LAPACK dspov
   real (kind=double) :: expTerm
   real (kind=double), allocatable, dimension (:,:) :: gaussME ! Matrix Element
   real (kind=double), allocatable, dimension (:,:) :: b ! For Ax=b

   interface
      subroutine dposv ( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
         use O_Kinds
         implicit none
         character :: uplo
         integer   :: info, lda, ldb, n, nrhs
         real (kind=double), dimension(lda,n) :: a
         real (kind=double), dimension(ldb,nrhs) :: b
      end subroutine dposv
   end interface

   ! Allocate space to hold temp data.
   allocate(gaussME(numTerms,numTerms))
   allocate(b(numTerms,1))

   ! Compute the numerical values of each gaussian, and the initial total
   !   function from the sum of gaussians.
   gaussFnSum(:) = 0.0_double
!   gaussFn(:,:) = 0.0_double
   do i = 1, numTerms
      do j = 1, numPoints
         expTerm = exp(-alphas(i)*radialValueSqrd(j))
!         gaussFn(j,i) = gaussFn(j,i) + expTerm
         gaussFn(j,i) = expTerm
         gaussFnSum(j)  = gaussFnSum(j) + expTerm
      enddo
   enddo

   ! Compute the matrix element between two gaussian functions.  (A from Ax=b)
   gaussME(:,:) = 0.0_double
   do i = 1, numTerms
      do j = 1, i
         do k = 1, numPoints-1
            gaussME(j,i) = gaussME(j,i) + &
                  & (gaussFn(k,i)+gaussFn(k+1,i)) * &
                  & (gaussFn(k,j)+gaussFn(k+1,j)) * weightSqrd(k)
         enddo
      enddo
   enddo

   ! Compute the right-hand-side of the linear problem (b from Ax=b).
   b(:,1) = 0.0_double
   do i = 1, numTerms
      do j = 1, numPoints-1
         b(i,1) = b(i,1) + &
               & (gaussFn(j,i)+gaussFn(j+1,i)) * &
               & (dataFn(j)+dataFn(j+1)) * weightSqrd(j)
      enddo
   enddo

   ! Solve for the coefficients (x) from Ax=b.  (Solution saved in b.)
   call dposv('U',numTerms,1,gaussME,numTerms,b,numTerms,info)

   ! Save the results of the calculation (the A coefficients).
   As(:) = b(:,1)

   ! Deallocate unnecessary arrays.
   deallocate (gaussME)
   deallocate (b)

end subroutine gaussCalc

! This function will calculate the root-mean-square error for a given
!   numerical function and a sum of analytic gaussian fitting functions.  It is
!   called to perform the same function with two slightly different parameter
!   sets.  In the first case, the gaussian functions are completely defined
!   except for the coefficients (As).  In the second case the gaussien
!   exponential alphas are new and so the gaussian functions and the associated
!   coefficients (As) need to be recalculated.  In the first case, the variable
!   'coeff' holds the As, in the second case it holds the exponential coeffs.
subroutine getRMS

   ! Import necessary definitions.
   use O_Kinds

   ! Import necessary data modeuls.
   use ExecutionData
   use GaussianData
   use GivenData
   use FittingData

   implicit none

   ! Declare local variables.
   integer :: i,j
   real (kind=double) :: RMS


   if (fitByAlphas == 0) then
      As(:) = coeffs(:)
   else
      alphas(:) = coeffs(:)
      call gaussCalc
   endif

   ! Compute the sum of gaussians with the appropriate coefficients.
   gaussFnSum(:) = 0.0_double
   do i = 1, numTerms
      do j = 1, numPoints
         gaussFnSum(j) = gaussFnSum(j) + As(i) * gaussFn(j,i)
      enddo
   enddo

!   ! Compute the derivative of the numerical sum of gaussians.
!   call computeDerivative(radialValue,gaussFnSum,gaussFnSumdr,numPoints)

   ! Compute the RMS between the given numerical data and the numerical
   !   gaussians.
   RMS = 0.0_double
   do i = 1, numPoints-1
!      RMS = RMS + sqrt(((dataFn(i)+dataFn(i+1)) - &
!                      & (gaussFnSum(i)+gaussFnSum(i+1)))**2)*weight(i)
      RMS = RMS + (sqrt(((dataFn(i)+dataFn(i+1)) - &
                      & (gaussFnSum(i)+gaussFnSum(i+1)))**2)*weight(i)) * &
                      & ((radialValue(i) + radialValue(i+1))**radialWeight) / &
                      & 2.0_double
   enddo

   ! Correct for the case where the radial value is not used to weight the
   !   fitting.  This is because each term was multiplied by 1/2.
   if (radialWeight == 0) then
      RMS = RMS * 2.0_double
   endif

!   ! Compute the RMS between the derivative of the given numerical data and the
!   !   derivative of the numerical gaussians.
!   RMSdr = 0.0_double
!   do i = 150, numPoints-1
!      RMSdr = RMSdr + (sqrt(((dataFndr(i)+dataFndr(i+1)) - &
!                      & (gaussFnSumdr(i)+gaussFnSumdr(i+1)))**2)*weight(i)) * &
!                      & ((radialValue(i) + radialValue(i+1))**1)/2.0_double
!   enddo

!   compRMS = RMS + RMSdr
   testCompRMS = RMS

end subroutine getRMS


!subroutine computeDerivative(radialValue,dataFn,dataFndr,numPoints)
!
!   use O_Kinds
!
!   implicit none
!
!   ! Define the passed dummy parameters.
!   real (kind=double), dimension(numPoints) :: radialValue
!   real (kind=double), dimension(numPoints) :: dataFn
!   real (kind=double), dimension(numPoints) :: dataFndr
!   integer :: numPoints
!
!
!   ! Define local variables.
!   integer :: i
!   real (kind=double) :: forwardRDiff
!   real (kind=double) :: backwardRDiff
!   real (kind=double) :: diffSum
!
!
!   ! The derivative at the first point and last points are approximated as the
!   !   slope of the line between the end points and their neighbors.
!   dataFndr(1) = (dataFn(2) - dataFn(1)) / (radialValue(2) - radialValue(1))
!   dataFndr(numPoints) = (dataFn(numPoints) - dataFn(numPoints-1)) / &
!         & (radialValue(numPoints) - radialValue(numPoints-1))
!
!   ! The remaining derivatives are computed as the weighted average of the
!   !   slopes between the target point and its forward and backward neighbors.
!   !   The weight is determined by the ratio of the distance of the neighbor
!   !   points from the target point.
!   do i = 2, numPoints-1
!      forwardRDiff  = radialValue(i+1) - radialValue(i)
!      backwardRDiff = radialValue(i) - radialValue(i-1)
!      diffSum = forwardRDiff + backwardRDiff
!
!      ! Simplifed (D=delta; +index = (i+1)-(i); -index = (i)-(i-1)):
!      !   df/dr = ((Df+ / Dr+)*(1/Dr+) + (Df- / Dr-)*(1/Dr-)) * (Dr+ + Dr-) / 2
!      !   The (1/Dr) terms and the (Dr+ + Dr-) are used to weight the Df/Dr
!      !   terms such that the non uniform radial values data points can be
!      !   used.
!
!      dataFndr(i) =((dataFn(i+1)-dataFn(i))/(forwardRDiff**2)   + &
!            &       (dataFn(i)-dataFn(i-1))/(backwardRDiff**2)) * &
!            &       diffSum/2.0_double
!   enddo
!end subroutine computeDerivative


subroutine patternMove

   ! Import necessary definitions.
   use O_Kinds

   ! Import necessary data modules.
   use GaussianData
   use FittingData

   implicit none

   ! Declare local variables
   integer :: i
   real (kind=double) :: tempCoeff

   ! Save the current coefficients as the best coefficients so far.  Also,
   !   make the next guess equal to double the step size between the current
   !   and the best.
   do i = 1, numTerms
      tempCoeff = coeffs(i)
      coeffs(i) = 2.0_double*coeffs(i) - bestCoeffs(i)
      bestCoeffs(i) = tempCoeff
   enddo

end subroutine patternMove

subroutine outputResults

   ! Import necessary definitions.
   use O_Kinds

   ! Import necessary data modules.
   use ExecutionData
   use GaussianData
   use GivenData
   use FittingData

   implicit none

   ! Declare local variables.
   integer :: i,j

   ! Write out the numerical data.
!   do i = 1, numPoints
!      write (12,*) radialValue(i),gaussFnSum(i)
!   enddo

   ! Write out the coefficients, exponential terms, and the final RMS.
!   write (6,*)

   if (fitByAlphas == 1) then

      ! Get the current set of As for this converged set of alphas.
      alphas(:) = coeffs(:)
      call gaussCalc

      ! Compute the sum of gaussians with the appropriate coefficients.
      gaussFnSum(:) = 0.0_double
      do i = 1, numTerms
         do j = 1, numPoints
            gaussFnSum(j) = gaussFnSum(j) + As(i) * gaussFn(j,i)
         enddo
      enddo
   else
      As(:) = coeffs(:)
   endif


   write (44,*) numTerms
   do i = 1, numTerms
      write (44,fmt="(2e18.8)") As(i), alphas(i) ! Coeffs, exp alphas
!      write (6,*) i,coeff(i),alphas(i)
   enddo
!   write (20,fmt="(a12,e16.8,a15,e16.8)") "Final RMS = ",bestRMS,&
!         & " Final RMSdr = ",bestRMSdr
!   close (20)

   do i = 1, numPoints
      write (20,*) radialValue(i),dataFn(i),gaussFnSum(i),dataFndr(i),&
            & gaussFnSumdr(i)
   enddo

   close (20)
   close (44)

   stop
end subroutine outputResults
