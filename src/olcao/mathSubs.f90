module O_MathSubs

   ! Import necessary modules.
   use O_Kinds
   use O_Constants

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

! This is a calculation of an integral over the error function.
function error_fn(x)

   ! Use precision parameters
   use O_Kinds

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double) :: x
   real (kind=double) :: error_fn

   ! Define the local variables used in this subroutine.
   real (kind=double) :: absX
   real (kind=double) :: absXSqrd
   real (kind=double) :: invAbsXSqrd
   real (kind=double) :: invSqrtPi
   real (kind=double) :: a
   real (kind=double), dimension(6) :: p, q

!-----------------------------------------------
!
!     SANDIA MATHEMATICAL PROGRAM LIBRARY
!     APPLIED MATHEMATICS DIVISION 2613
!     SANDIA LABORATORIES
!     ALBUQUERQUE, NEW MEXICO  87185
!     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                    ISSUED BY SANDIA LABORATORIES                     *
!  *                   A PRIME CONTRACTOR TO THE                       *
!  *                UNITED STATES DEPARTMENT OF ENERGY                 *
!  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * * *
!  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE   *
!  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE      *
!  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,    *
!  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES  *
!  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL      *
!  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR     *
!  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS      *
!  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE          *
!  * OWNED RIGHTS.                                                     *
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS     *
!  * PART IS SAND77-1441.                                              *
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     WRITTEN BY J.E. VOGEL FROM APPROXIMATIONS DERIVED BY W.J. CODY .
!
!     ABSTRACT
!
!          ERF(X) COMPUTES 2.0/SQRT(PI) TIMES THE INTEGRAL FROM 0 TO X
!          OF EXP(-X**2). THIS IS DONE USING RATIONAL APPROXIMATIONS.
!          ELEVEN CORRECT SIGNIFICANT FIGURES ARE PROVIDED.
!
!     DESCRIPTION OF PARAMETERS
!
!          X MAY BE ANY REAL VALUE
!
!     ERF IS DOCUMENTED COMPLETELY IN SC-M-70-275
!
! Note that this procedure has been modified to use if-then-else-endif
!   statements instead of GO TO statements.  It has also been changed so that
!   the appropriate constants will be defined only when explicitly needed
!   for each case.  For the OLCAO purposes, that great majority of the time
!   is spent in the > 6.0 case where there is no need for any of the constants
!   defined previously above.  So, those assignments have been moved into the
!   body of the function.  Finally, to improve the efficiency even more, the
!   number of arrays that need to be allocated for this routine has been
!   reduced from 6 to 2.  The p1, p2, and p3 are now just p, and the same with
!   the q1, q2, and q3.  This is because only one of each are used for any
!   given case.


   absX = abs(x)
   if (absX <= 6.0_double) then
      absXSqrd = absX * absX
      if (absX <= 4.0_double) then
         if (absX <= 0.46875_double) then
            ! Assign parameters where x is between 0.46875 and 4.
            p(1) = 242.6679552305318_double
            p(2) = 21.97926161829415_double
            p(3) = 6.996383488619136_double
            p(4) = -0.03560984370181539_double
            q(1) = 215.0588758698612_double
            q(2) = 91.16490540451490_double
            q(3) = 15.08279763040779_double
            q(4) = 1.000000000000000_double

            a = absX*(p(1) + absXSqrd * (p(2) + &
              & absXSqrd * (p(3) + absXSqrd * p(4))))
            a = a  / (q(1) + absXSqrd * (q(2) + &
              & absXSqrd * (q(3) + absXSqrd * q(4))))
            if (x < 0.0_double) a = -a
            error_fn = a
         else
            ! Assign parameters where x is between 0.46875 and 4.
            p(1) = 22.898992851659_double
            p(2) = 26.094746956075_double
            p(3) = 14.571898596926_double
            p(4) = 4.2677201070898_double
            p(5) = 0.56437160686381_double
            p(6) = -0.0000060858151959688_double
            q(1) = 22.898985749891_double
            q(2) = 51.933570687552_double
            q(3) = 50.273202863803_double
            q(4) = 26.288795758761_double
            q(5) = 7.5688482293618_double
            q(6) = 1.0000000000000_double

            a = exp((-absXSqrd)) * (p(1) + absX * (p(2) + &
              & absX * (p(3) + absX * (p(4) + absX * (p(5) + &
              & absX * p(6))))))
            a = a/(q(1) + absX * (q(2) + absX * (q(3) + absX * &
              & (q(4) + absX * (q(5) + absX * q(6))))))
            error_fn = sign(1.0_double - a,x)
         endif
      else
         ! Assign parameters where x is between 4 and 6.
         invSqrtPi = 0.564189583547756_double
         p(1) = -0.0121308276389978_double
         p(2) = -0.1199039552681460_double
         p(3) = -0.243911029488626_double
         p(4) = -0.0324319519277746_double
         q(1) = 0.0430026643452770_double
         q(2) = 0.489552441961437_double
         q(3) = 1.43771227937118_double
         q(4) = 1.00000000000000_double

         invAbsXSqrd = 1.0_double / absXSqrd
         a = invAbsXSqrd * (p(1) + invAbsXSqrd * (p(2) + invAbsXSqrd * &
           & (p(3) + invAbsXSqrd * p(4)))) / (q(1) + invAbsXSqrd * (q(2) + &
           & invAbsXSqrd * (q(3) + invAbsXSqrd * q(4))))
         a = exp(-absXSqrd) * (invSqrtPi + a) / absX
         error_fn = sign(1.0_double - a,x)
      endif
   else
      error_fn = x/absX
   endif

end function error_fn


function stepFunction (x,stepFnRange)
   use O_Kinds
   implicit none
   real (kind=double), intent(in) :: x
   real (kind=double), intent(in) :: stepFnRange
   real (kind=double) :: stepFunction

   ! Note regarding the value of the step function range. The default input
   !   value is currently 11.5. This will cause the fermi function to cut
   !   off (and not evaluate the exponential) when the fermi function would
   !   evaluate to something any closer to 1.0 than 0.99999 or any closer to
   !   0.0 than 0.00001. This can be seen by solving the stepFunction for x
   !   x=ln(1/stepFn - 1) and putting 0.99999 or 0.00001 in for stepFn.

   if     (x >  stepFnRange) then
      stepFunction = 0.0_double
   elseif (x < -stepFnRange) then
      stepFunction = 1.0_double
   else
      stepFunction = 1.0_double / (1.0_double + exp(x))
   endif
end function stepFunction


subroutine crossProduct (answer, vector1, vector2)

   ! Import the necessary kinds definitions.
   use O_Kinds

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define the passed parameters.
   real (kind=double), dimension(3) :: answer
   real (kind=double), dimension(3) :: vector1  ! ax + ay + az
   real (kind=double), dimension(3) :: vector2  ! bx + by + bz

   ! Compute the cross product.  Note the correctly reversed sign when
   !   computing answer(2).
   answer(1) = vector1(2)*vector2(3) - vector1(3)*vector2(2) ! ay*bz - az*by
   answer(2) = vector1(3)*vector2(1) - vector1(1)*vector2(3) ! az*bx - ax*bz
   answer(3) = vector1(1)*vector2(2) - vector1(2)*vector2(1) ! ax*by - ay*bx

end subroutine crossProduct


end module O_MathSubs
