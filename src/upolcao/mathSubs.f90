module O_MathSubs

   ! Import necessary modules.
   use O_Kinds
   use O_Constants

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   ! Define module variables.
   integer(8), allocatable, dimension(:) :: preCompFactorial

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


subroutine computeFactorials(maxFact)

   implicit none

   ! Define the passed parameters.
   integer :: maxFact ! Highest factorial parameter.

   ! Define local variables.
   integer(8) :: i

   allocate (preCompFactorial(0:maxFact))

   preCompFactorial(0) = 1 ! This is 0!
   preCompFactorial(1) = 1 ! This is 1!
   do i = 2, maxFact
      preCompFactorial(i) = preCompFactorial(i-1) * i
   enddo

end subroutine computeFactorials


function wignerD(twoj, twom, twomp, eulerCoords)

   ! Use necessary modules.
   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! Define function definition variable.
   complex (kind=double) :: wignerD

   ! Define passed parameters.
   integer :: twoj, twom, twomp
   real(kind=double), dimension(3) :: eulerCoords ! alpha, beta, gamma
         ! Pulled from theta and phi of the ttpCoords.

   ! Define local variables.
   real (kind=double) :: term1Exp, term3Exp
   complex (kind=double) :: term1, term3
   real (kind=double) :: term2


   ! Condition 1
!   if ((j < 0) .or. (.not. (modulo(j,1.0) == 0)) .or. &
!         & (modulo(j,1.0) == 0.5 .and. ((modulo(m,1.0) /= 0) &
!         & .or. (modulo(mp,1.0) /= 0)))) then
!      write (20,*) "Invalid input parameters:"
!      write (20,fmt="(a,3e10.3)") "j, m, mp: ", j, m, mp
!      write (20,*) "Parameter j must be non-negative integer or half-integer."
!      write (20,*) "Parameters m and mp must be between -j and j."
!      stop
!   endif

!   ! Condition 2
!   if ((eulerCoords(1) < 0) .or. (eulerCoords(1) > 2.0 * pi) &
!         & .or. (eulerCoords(2) < 0) .or. (eulerCoords(2) > pi) &
!         & .or. (eulerCoords(3) < 0) .or. (eulerCoords(3) > 2.0 * pi)) then
!      write (20,*) "Invalid input parameters:"
!      write (20,fmt="(a,3e12.3)") "phi theta phi: ", eulerCoords(1), &
!            & eulerCoords(2), eulerCoords(3)
!      write (20,fmt="(a,3e12.3)") "phi theta phi: ",2.0*pi,2.0*pi,2.0*pi
!   endif
!   if ((theta_0 < 0) .or. (theta_0 > 2.0 * pi) .or. (theta < 0) &
!         & .or. (theta > 2.0 * pi) .or. (phi < 0) .or. (phi > 2.0 * pi)) then
!      write (20,*) "Invalid input parameters:"
!      write (20,fmt="(a,3e10.3)") "j, m, mp: ", j, m, mp
!      write (20,*) "Parameter j must be non-negative integer or half-integer."
!   endif
!write (6,fmt="(a,3i3,3e14.5)") "twoj twom twomp EulerCoords",twoj,&
!   & twom,twomp,eulerCoords(:)

   ! Compute term1.
   term1Exp = real(twom) / 2.0_double * eulerCoords(1)
   term1 = cmplx(cos(term1Exp), -sin(term1Exp), double)

   ! Compute term2.
   term2 = smalld(twoj, twom, twomp, eulerCoords(2))

   ! Compute term3.
   term3Exp = real(twomp) / 2.0_double * eulerCoords(3)
   term3 = cmplx(cos(term3Exp), -sin(term3Exp), double)

!write (6,fmt="(a,2e15.4)") "WignerD term1", term1
!write (6,fmt="(a,2e15.4)") "WignerD term2", term2
!write (6,fmt="(a,2e15.4)") "WignerD term3", term3

   wignerD = term1 * term2 * term3

end function wignerD


function smalld (twoj, twom, twomp, theta)

   ! Use necessary modules
   use O_Kinds

   implicit none

   ! Define function return variable.
   real(kind=double) :: smalld

   ! Define passed parameters.
   integer :: twoj, twom, twomp
   real(kind=double) :: theta

   ! Define local variables.
   integer :: twok, twokStart, twokEnd, halfk
   real(kind=double) :: dsum
!   real(kind=double) :: coeff
   real(kind=double) :: cos_theta_2
   real(kind=double) :: sin_theta_2
   real(kind=double) :: numerator
   real(kind=double) :: denominator


   ! Following equation (4) in section 4.3.1 of "Quantum Theory of Angular
   !   Momentum: Irreducible Tensors, Spherical Harmonics, Vector Coupling
   !   Coefficients, 3nj Symbols" by Varshalovich DA, Moskalev AN, and
   !   Khersonski VK.; Singapore; Teaneck, NJ, USA: World Scientific Pub;
   !   1988. 514 p.; Equation found at the bottom of page 76.
   

   cos_theta_2 = cos(0.5_double * theta)
   sin_theta_2 = sin(0.5_double * theta)
!write (6,*) "cos sin",cos_theta_2,sin_theta_2

   twokStart = max(0, twom - twomp)
   twokEnd = min(twoj + twom, twoj - twomp)

!write (6,*) "twoj, twom, twomp = ", twoj, twom, twomp
!write (6,*) "twokStart twokEnd", twokStart, twokEnd
   smalld = sqrt( &
         &   real(preCompFactorial((twoj+twom)/2),double) &
         & * real(preCompFactorial((twoj-twom)/2),double) &
         & * real(preCompFactorial((twoj+twomp)/2),double) &
         & * real(preCompFactorial((twoj-twomp)/2),double))
!write (6,*) "smalld", smalld

   dsum = 0.0_double
   do twok = twokStart, twokEnd, 2
      halfk = int(real(twok) / 2.0_double)
      numerator = (-1.0d0)**halfk &
            & * cos_theta_2**(twoj - twok + (twom - twomp)/2) &
            & * sin_theta_2**(twok + (twomp - twom)/2)
!write (6,*) "(twom - twomp)/2",(twom - twomp)/2
!write (6,*) "(twomp - twom)/2",(twomp - twom)/2
!write (6,fmt="(a,i3,3e14.5)") "twok 1 2 3", twok, (-1.0d0)**halfk, &
!      & cos_theta_2**(twoj - twok + (twom - twomp)/2), &
!      & sin_theta_2**(twok + (twomp - twom)/2)

      denominator = real(preCompFactorial(halfk),double) &
            & * real(preCompFactorial((twoj + twom - twok)/2),double) &
            & * real(preCompFactorial((twoj - twomp - twok)/2),double) &
            & * real(preCompFactorial((twomp - twom + twok)/2),double)

!write (6,*) "numerator denom ", numerator, denominator
      dsum = dsum + numerator / denominator
   enddo

   smalld = smalld * dsum

end function smalld


function hypersphericalHarmonic4D(twoj,twom,twomp,ttpCoords)

   ! Use necessary modules
   use O_Kinds

   implicit none

   ! Define function return variable.
   complex(kind=double) :: hypersphericalHarmonic4D

   ! Define passed parameters.
   integer :: twoj, twom, twomp
   real(kind=double), dimension(3) :: ttpCoords ! theta_0, theta, phi

   ! Define local variables.
   integer :: i
   real(kind=double) :: term2Param
   real(kind=double), dimension(3) :: ttpTempCoords
   complex(kind=double) :: term1, term2, term3

   ! Initialize the accumulation variable.
   hypersphericalHarmonic4D = cmplx(0.0, 0.0, double)
   do i = -twoj, twoj, 2
      ttpTempCoords(1) = ttpCoords(3)  !  Phi
      ttpTempCoords(2) = ttpCoords(2)  !  Theta
      ttpTempCoords(3) = -ttpCoords(3) ! -Phi

      ! Compute term 1.
      term1 = wignerD(twoj, twom, i, ttpTempCoords)

      ! Compute term 2.
      term2Param = real(i)/2.0_double * ttpCoords(1)
      term2 = cmplx(cos(term2Param), -sin(term2Param), double)

      ! Compute term 3.
      ttpTempCoords(2) = -ttpCoords(2) ! -Theta
      term3 = wignerD(twoj, i, twomp, ttpTempCoords)

      ! HSH = product(term1,term2,term3).
      hypersphericalHarmonic4D = hypersphericalHarmonic4D &
            & + term1 * term2 * term3

!write (6,*) "i HSH = ", i, hypersphericalHarmonic4D
!write (6,*) "term1 = ", term1
!write (6,*) "term2 = ", term2
!write (6,*) "term3 = ", term3
   enddo
   
end function hypersphericalHarmonic4D


function clebschGordan(twoj1, twoj2, twoj, twom1, twom2, twom)

   ! Use necessary modules
   use O_Kinds

   implicit none

   ! Define function return variable.
   real (kind=double) :: clebschGordan

   ! Define passed paramters.
   integer :: twoj1, twoj2, twoj, twom1, twom2, twom

   ! Define local variables.
   integer :: twoz ! Index in the summation for computing CGC values.
   integer :: twocgcMin ! Min value of twoz.
   integer :: twocgcMax ! Max value of twoz.
   real(kind=double) :: preFactor
   real(kind=double) :: coefficient
   real(kind=double) :: numerator
   real(kind=double) :: denominator
!write (6,fmt="(a,i)") "(twoj1 + twoj2 - twoj)/2", (twoj1 + twoj2 - twoj)/2
!write (6,fmt="(a,i)") "(twoj1 - twoj2 + twoj)/2", (twoj1 - twoj2 + twoj)/2
!write (6,fmt="(a,i)") "(-twoj1 + twoj2 + twoj)/2", (-twoj1 + twoj2 + twoj)/2
!write (6,fmt="(a,i)") "(twoj + twoj1 + twoj2 + 2)/2", (twoj + twoj1 + twoj2 + 2)/2

   ! Following QTAM chapter 8, we see in sub-section 8.1.1 the constraints on
   !   parameters (j, j1, j2, m, m1, m2) to the CGCs, and in section 8.2 we
   !   see the explicit definition of the CGC values. Specifically, we use
   !   equations (1) from section 8.2 of QTAM (page 237) to define the
   !   preFactor (Delta(abc) with a = j, b = j1, c = j2).
   preFactor = sqrt( &
         & real(preCompFactorial((twoj1 + twoj2 - twoj)/2),double) &
         & * real(preCompFactorial((twoj1 - twoj2 + twoj)/2),double) &
         & * real(preCompFactorial((-twoj1 + twoj2 + twoj)/2),double) &
         & / real(preCompFactorial((twoj + twoj1 + twoj2 + 2)/2),double))
!write (6,fmt="(a,e15.4)") "preFactor", preFactor
!
!write (6,fmt="(a,i)") "(twoj + twom)/2", (twoj + twom)/2
!write (6,fmt="(a,i)") "(twoj - twom)/2", (twoj - twom)/2
!write (6,fmt="(a,i)") "(twoj + 1)", twoj + 1
!write (6,fmt="(a,i)") "(twoj1 + twom1)/2", (twoj1 + twom1)/2
!write (6,fmt="(a,i)") "(twoj1 - twom1)/2", (twoj1 - twom1)/2
!write (6,fmt="(a,i)") "(twoj2 + twom2)/2", (twoj2 + twom2)/2
!write (6,fmt="(a,i)") "(twoj2 - twom2)/2", (twoj2 - twom2)/2

   ! Then, we use equation (5) in section 8.2 of QTAM (page 238) to compute
   !   explicit values for the CGCs.
   coefficient = sqrt( &
      & real(preCompFactorial((twoj + twom)/2),double) &
      & * real(preCompFactorial((twoj - twom)/2) * (twoj + 1),double) &
      & / real(preCompFactorial((twoj1 + twom1)/2),double) &
      & * real(preCompFactorial((twoj1 - twom1)/2),double) &
      & * real(preCompFactorial((twoj2 + twom2)/2),double) &
      & * real(preCompFactorial((twoj2 - twom2)/2),double))
!write (6,fmt="(a,e15.4)") "coefficient", coefficient

   clebschGordan = 0.0_double
   twocgcMin = max(0, twom1 - twoj1, twoj2 - twoj1 + twom)
   twocgcMax = min(twoj2 + twoj + twom1, twoj - twoj1 + twoj2,&
         & twoj + twom)
!write (6,fmt="(a,2i)") "cgc Min,Max", twocgcMin, twocgcMax
!   cgcMin = (max(0, twom1 - twoj1, twoj2 - twoj1 + twom))/2
!   cgcMax = (min(twoj2 + twoj + twom1, twoj - twoj1 + twoj2,&
!         & twoj + twom))/2
!!   clebschGordan = 0.0_double
!!   cgcMin = max(0, int(m1 - j1), int(j2 - j1 + m))
!!   cgcMax = min(int(j2 + j + m1), int(j - j1 + j2), int(j + m))

   do twoz = twocgcMin, twocgcMax, 2
!write (6,*) "twoz = ", twoz
!write (6,fmt="(a,i)") "(twoj2 + twom2 + i)/2", (twoj2 + twom2 + twoz)

      numerator = real(((-1)**((twoj2 + twom2 + twoz)/2)),double) &
         & * real(preCompFactorial((twoj + twoj2 + twom1 - twoz)/2),double) &
         & * real(preCompFactorial((twoj1 - twom1 + twoz)/2),double)
!      numerator = real(((-1)**(j2 + m2 + i)) &
!            & * preCompFactorial(int(j2 + j + m1 - i)) &
!            & * preCompFactorial(int(j1 - m1 + i)),double)

      denominator = real(preCompFactorial(twoz/2),double) &
         & * real(preCompFactorial((twoj - twoj1 + twoj2 - twoz)/2),double) &
         & * real(preCompFactorial((twoj + twom - twoz)/2),double) &
         & * real(preCompFactorial((twoj1 - twoj2 - twom + twoz)/2),double)
!      denominator = real(preCompFactorial(i) &
!            & * preCompFactorial(int(j - j1 + j2 - i)) &
!            & * preCompFactorial(int(j + m - i)) &
!            & * preCompFactorial(int(j1 - j2 - m + i)),double)

      clebschGordan = clebschGordan + numerator / denominator
   enddo

   clebschGordan = preFactor * coefficient * clebschGordan

end function clebschGordan

end module O_MathSubs
