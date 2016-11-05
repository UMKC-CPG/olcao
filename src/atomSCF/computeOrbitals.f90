module ComputeOrbitalsSubs

private
public :: computeOrbitals

contains

subroutine computeOrbitals

   ! Use necessary modules.
   use ExecutionData

   implicit none

   ! On the first iteration (or any time when the energy eigen values are not
   !   sufficiently accurate for method #2) use method #1 to solve for the
   !   energy eigen values.
   if (solverChoice == 1) then
      call solveWaveFn1
   else
      call solveWaveFn2
   endif

end subroutine computeOrbitals


! Finds the non-relativistic wave function using finite differences and matrix
!   diagonalizations.  An initial guess for the eigenvalues is not needed.
subroutine solveWaveFn1

   ! Use necessary modules.
   use O_Kinds
   use O_RadialGrid
   use AtomData
   use ChargeDensityData
   use WaveFnData
   use EnergyData
   use PotData

   implicit none

   ! Define local variables.
   integer :: i,j,k,l
   integer :: llM1  ! lQN * (lQN - 1)
   integer :: ki
   integer, allocatable, dimension(:,:) :: maxQNn
   real (kind=double) :: orbRadialRho
   real (kind=double) :: a
   real (kind=double) :: b
   real (kind=double) :: c1
   real (kind=double) :: c2
   real (kind=double), allocatable, dimension (:)   :: dk

   ! Define variables specific to the LAPACK call (include eigen value and
   !   eigen vector results).
   integer :: info
   integer :: numFound ! Number of eigen values found.
   integer, allocatable, dimension(:)   :: iwork
   integer, allocatable, dimension(:)   :: ifail
   real (kind=double) :: abstol
   real (kind=double), allocatable, dimension (:)   :: work
   real (kind=double), allocatable, dimension (:)   :: d
   real (kind=double), allocatable, dimension (:)   :: sd
   real (kind=double), allocatable, dimension (:)   :: w  ! Energy Eigen Values
   real (kind=double), allocatable, dimension (:,:) :: z  ! Wave Function

   ! Define the interface to the LAPACK subroutine for obtaining the eigen
   !   values and vectors from a symmetric tridiagonal matrix.
   interface
      subroutine dstevx (JOBZ,RANGE,N,D,E,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,&
            & IWORK,IFAIL,INFO)
         use O_Kinds
         character :: JOBZ
         character :: RANGE
         integer :: N
         integer :: M
         integer :: IL
         integer :: IU
         integer :: LDZ
         integer :: INFO
         real (kind=double) :: ABSTOL
         real (kind=double) :: VL
         real (kind=double) :: VU
         integer, dimension (5*N) :: IWORK
         integer, dimension (N)   :: IFAIL
         real (kind=double), dimension(N)     :: D
         real (kind=double), dimension(N-1)   :: E
         real (kind=double), dimension(N)     :: W
         real (kind=double), dimension(LDZ,M) :: Z
         real (kind=double), dimension(5*N)   :: WORK
      end subroutine dstevx

      subroutine dstev (JOBZ,N,D,E,Z,LDZ,WORK,INFO)
         use O_Kinds
         character :: JOBZ
         integer :: N
         integer :: LDZ
         integer :: INFO
         real (kind=double), dimension(N)     :: D
         real (kind=double), dimension(N-1)   :: E
         real (kind=double), dimension(LDZ,N) :: Z
         real (kind=double), dimension(2*N-2) :: WORK
      end subroutine dstev
   end interface

   ! Allocate space to hold the maximum n QN for each l,s pair of all orbitals.
   allocate (maxQNn(2,maxQNl)) ! 2 = spin down (1), spin up (2)

   ! Allocate space to hold the diagonal elements (before and after the hartree
   !   potential effect is accounted for), off-diagonal elements, real and
   !   integer work space, convergence testing space, and the resultant eigen
   !   values.  (The eigen vectors are stored in an array allocated separately
   !   for each l,s pair later.)
   allocate (d(numRadialPoints-1))
   allocate (dk(numRadialPoints-1))
   allocate (sd(numRadialPoints-2))
   allocate (work(5*(numRadialPoints-1)))
   allocate (iwork(5*(numRadialPoints-1)))
   allocate (ifail(numRadialPoints-1))
   allocate (w(numRadialPoints-1))

   ! Initialize the energy eigen values for all orbitals.
   eigenValues(:numOrbitals) = 0.0_double

   ! Initialize the tolerance for diagonalization.  The negative value will
   !   force the use of the machine dependent value within the LAPACK library.
   abstol = -1.0_double

   ! Initialize the current values for charge density.
   rhoDn(:) = 0.0_double
   rhoUp(:) = 0.0_double

   ! Find the maximum n quantum number given l and s.  (NOTE:  Zero spin is
   !   treated as down.)
   do i = 1,2
      do j = 1, maxQNl

         ! Initialize the maximum QN for this spin-orbital angmo pair.
         maxQNn(i,j) = 0

         do k = 1, numOrbitals

            ! The maximum n QN must be non-negative.
            if (orbitalQNn(k) <= 0) cycle

            ! The maximum n QN must belong to an orbital with angular momentum
            !   QN l such that n-1=l.  Recall that the maxQNl value is actually
            !   the maximum l QN + 1 for the purpose of array access etc.
            if (orbitalQNl(k) /= j-1) cycle

            ! The maximum n QN must belong to an orbital with spin up or spin
            !   down in the spin polarized case, and spin down (zero) only in
            !   the non spin polarized case.
            if ((orbitalQNs(k)-0.1_double)*(i-1.5_double) < 0.0_double) cycle

            ! Save the n QN for this pair since only the max value could not
            !   have cycled yet.
            maxQNn(i,j) = orbitalQNn(k)
         enddo
      enddo
   enddo


   ! Set up the Hamiltonian matrix for the kinetic energy.  Only the diagonal
   !   depends on the potential.  I have no idea what most of these paramters
   !   mean or do.  This loop was modified from the original.  The index i
   !   goes from 2 to numRadialPoints instead of 3 to numRadialPoints.  The
   !   idea is to put the diagonal elements in the indices 1 to
   !   numRadialPoints-1.  Some theoretical backing here would be really really
   !   helpful to make sure that such simple modifications do not introduce
   !   errors.  Now we must be concerned with the values in dk and d.
   a  = exp(-aaWhatever) / atomicNumber
   b  = 1.0_double / bbWhatever
   c2 = -1.0_double / b**2.0_double
   c1 = -2.0_double * c2 + 0.25_double
   dk(1) = c1 / (radialPoints(2) + a)**2.0_double
   do i = 3, numRadialPoints
      dk(i-1)  = c1 / (radialPoints(i) + a)**2.0_double
      sd(i-2)  = c2 / ((radialPoints(i) + a)*(radialPoints(i-1) + a))
   enddo

   ! Loop over spin down=1 and spin up=2
   do i = 1,2

      ! Start a loop over s, p, d, ... states
      do j = 1, maxQNl

         ! Skip any orbitals with the n QN == 0.  (Non existant in this calc.)
         if (maxQNn(i,j) == 0) cycle

         ! Compute the current value of lQN * (lQN - 1)
         llM1 = j*(j-1)

         ! Consider each radial value to apply the potential.
         do k = 2, numRadialPoints

            ! Spin down i = 1; spin up otherwise
            if (i == 1) then
               d(k-1) = dk(k-1) + (ionicPotDn(j,k) + llM1/radialPoints(k)) / &
                     & radialPoints(k) + hartreePotDn2(k)
            else
               d(k-1) = dk(k-1) + (ionicPotUp(j,k) + llM1/radialPoints(k)) / &
                     & radialPoints(k) + hartreePotUp2(k)
            endif
         enddo

         ! Allocate space necessary for this computation.
         allocate (z(numRadialPoints-1,maxQNn(i,j)))

         ! Diagonalize and obtain eigen values and eigen vectors of the first
         !   few (maxQNn(i,j) terms.
         call dstevx ('V','I',numRadialPoints-1,d,sd,0.0_double,1.0_double,&
               & 1,maxQNn(i,j),abstol,numFound,w,z,numRadialPoints-1,work,&
               & iwork,ifail,info)

         ! Save the energy levels and add to the charge density.
         ki = 1
         do k = 1, numOrbitals

            ! The maximum n QN must be non-negative.
            if (orbitalQNn(k) <= 0) cycle

            ! The maximum n QN must belong to an orbital with angular momentum
            !   QN l such that n-1=l.
            if (orbitalQNl(k) /= j-1) cycle

            ! The maximum n QN must belong to an orbital with spin up or spin
            !   down in the spin polarized case, and spin down (zero) only in
            !   the non spin polarized case.
            if ((orbitalQNs(k)-0.1_double)*(i-1.5_double) < 0.0_double) cycle

            ! Save the energy eigen value for this maxQNn.
            eigenValues(k) = w(ki)

            do l = 2, numRadialPoints

               ! Compute the charge contribution from this orbital at this
               !   radial distance.
               orbRadialRho = orbitalCharge(k) * z(l-1,ki)**2.0_double / &
                     & radialPoints2(l)

               ! Increment the charge density in the radial direction for the
               !   appropriate spin.
               if (i == 1) then
                  rhoDn(l) = rhoDn(l) + orbRadialRho
               else
                  rhoUp(l) = rhoUp(l) + orbRadialRho
               endif
            enddo

            ! Increment the eigen value index.
            ki = ki + 1
         enddo

         ! Deallocate the eigen vectors for this l,s pair.
         deallocate (z)
      enddo
   enddo

   ! Deallocate unnecessary arrays.
   deallocate (maxQNn)
   deallocate (d)
   deallocate (dk)
   deallocate (sd)
   deallocate (work)
   deallocate (iwork)
   deallocate (ifail)
   deallocate (w)


end subroutine solveWaveFn1



subroutine solveWaveFn2

   ! Use necessary modules
   use O_Kinds
   use AtomData
   use WaveFnData
   use O_RadialGrid
   use PotData
   use EnergyData
   use ChargeDensityData
   use ExecutionData

   implicit none

   ! Define local variables.
   integer :: i,j
   real (kind=double) :: orbRadialRho
   real (kind=double), allocatable, dimension (:) :: v

   allocate (v(numRadialPoints))
   allocate (waveFn(numRadialPoints))
   allocate (dr1WaveFn(numRadialPoints))

   ! Initialize the current values for charge density.
   rhoDn(:)   = 0.0_double
   rhoUp(:)   = 0.0_double
   rhoCore(:) = 0.0_double

   ! Begin looping over orbitals.  Note that spin zero is treated as down.
   do i = 1, numOrbitals

      ! The n QN must be non-negative.
      if (orbitalQNn(i) <= 0) cycle

      ! There must be charge present in this orbital (occupied) unless the
      !   system has converged.
      if ((orbitalCharge(i) == 0.0_double) .and. (converged == 0)) cycle

      ! All eigen values must negative.  I'm not sure I like this.  How about
      !   a rigid shift of all values to force the highest occupied orbital to
      !   have energy = 0?  This is a bit unnecessary anyway since the
      !   intNonRel subroutine will change the eigenValue too.
      if (eigenValues(i) >= 0.0_double) then
         eigenValues(i) = -0.0001_double
      endif

      ! Set up the potential.
      do j = 2, numRadialPoints

         ! Spin down or spin non polarized.
         if (orbitalQNs(i) < 0.1_double) then
            v(j) = ionicPotDn(orbitalQNl(i)+1,j) / radialPoints(j) + &
                  & hartreePotDn2(j)
         else ! Spin up.
            v(j) = ionicPotUp(orbitalQNl(i)+1,j) / radialPoints(j) + &
                  & hartreePotUp2(j)
         endif

         ! Some non relativistic thingy.  Dunno!?!
         if (doRelativistic == 0) then
            v(j) = v(j) + (orbitalQNl(i)*(orbitalQNl(i)+1))/ radialPoints(j)**2
         endif
      enddo

      ! Call the appropriate integration subroutine.
      if (doRelativistic == 0) then
         call intNonRel(v,i)
!      else
!         call intRel(v,i)
      endif

      ! Accumulate the charge density.
      do j = 1, numRadialPoints

         ! Compute the charge contribution from this orbital at this radial
         !   distance.
         orbRadialRho = orbitalCharge(i) * waveFn(j) * waveFn(j)

         if (doRelativistic == 1) then
            orbRadialRho = orbRadialRho + orbitalCharge(i) * dr1WaveFn(j) * &
                  & dr1WaveFn(j)
         endif

         ! Increment the charge density in the radial direction for the
         !   appropriate spin.
         if (orbitalQNs(i) < 0.1_double) then
            rhoDn(j) = rhoDn(j) + orbRadialRho
         else
            rhoUp(j) = rhoUp(j) + orbRadialRho
         endif

         ! Accumulate data for the core orbitals.
         if (i <= numCoreOrb) then
            rhoCore(j) = rhoCore(j) + orbRadialRho
         endif
      enddo

      ! Compute various quantities on the last iteration.
      if (converged == 1) then
         call orbitalInfo (i)
      endif

! Not used.  Strange.
!      ! Save the radial wave functions.
!      if ((i > numCoreOrb) .and. (doRelativistic == 0)) then
!         waveFnQNl(orbitalQNl(i)+1,2:) = waveFn(2:)
!      endif
   enddo

   deallocate (v)
   deallocate (waveFn)
   deallocate (dr1WaveFn)

end subroutine solveWaveFn2



subroutine intNonRel(v,iorb)

   ! Use necessary modules
   use O_Kinds
   use O_RadialGrid
   use ExecutionData
   use AtomData
   use PotData
   use EnergyData
   use WaveFnData

   implicit none

   ! Define passed paramters.
   real (kind=double), dimension (numRadialPoints) :: v
   integer :: iorb

   ! Define local variables.
   integer :: lp
   integer :: ll
   integer :: j
   integer :: jj
   integer :: nctp
   integer :: ninf
   integer :: juflow
   integer :: nodes
   integer :: istop
   integer :: itmax
   real (kind=double) :: expzer
   real (kind=double) :: tol
   real (kind=double) :: abc1
   real (kind=double) :: abc2
   real (kind=double) :: abc3
   real (kind=double) :: abc4
   real (kind=double) :: abc5
   real (kind=double) :: amc0
   real (kind=double) :: amc1
   real (kind=double) :: amc2
   real (kind=double) :: amc3
   real (kind=double) :: amc4
   real (kind=double) :: zeff
   real (kind=double) :: aa
   real (kind=double) :: a
   real (kind=double) :: bb
   real (kind=double) :: b
   real (kind=double) :: vzero
   real (kind=double) :: var0
   real (kind=double) :: emax
   real (kind=double) :: emin
   real (kind=double) :: fa0
   real (kind=double) :: fa1
   real (kind=double) :: fa2
   real (kind=double) :: fa3
   real (kind=double) :: fa4
   real (kind=double) :: fa5
   real (kind=double) :: fb0
   real (kind=double) :: fb1
   real (kind=double) :: fb2
   real (kind=double) :: fb3
   real (kind=double) :: fb4
   real (kind=double) :: fb5
   real (kind=double) :: arp
   real (kind=double) :: arc
   real (kind=double) :: arctp
   real (kind=double) :: brp
   real (kind=double) :: brc
   real (kind=double) :: brctp
   real (kind=double) :: alf
   real (kind=double) :: factor
   real (kind=double) :: dev
   real (kind=double) :: evold

!--machine dependent parameter-
!--require exp(-2*sqrt(expzer)) to be within the range of real*8
!
   expzer = 1.D3

   tol = 1.D-10
!
!  integration coefficients
!
   itmax = 50
   abc1 =  1901.D0/720.D0
   abc2 = -1387.D0/360.D0
   abc3 =   109.D0/30.D0
   abc4 =  -637.D0/360.D0
   abc5 =   251.D0/720.D0
   amc0 =   251.D0/720.D0
   amc1 =   323.D0/360.D0
   amc2 =   -11.D0/30.D0
   amc3 =    53.D0/360.D0
   amc4 =   -19.D0/720.D0

   b = 1.0_double / bbWhatever
   a = exp(-aaWhatever) / atomicNumber
!
!  determine effective charge and vzero for startup of
!  outward integration
!
!  ar = r**(l+1) * (1 + aa r + bb r**2 + ... )
!
!  aa = -znuc / lp     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)
!        znuc = 
!
   lp = orbitalQNl(iorb)+1
   zeff = 0.D0
   if ((orbitalQNs(iorb) .lt. 0.1D0) .and. (ionicPotDn(lp,2) .lt. -0.1D0)) then
      zeff=atomicNumber
   endif
   if ((orbitalQNs(iorb) .gt. 0.1D0) .and. (ionicPotUp(lp,2) .lt. -0.1D0)) then
      zeff=atomicNumber
   endif
   aa = -zeff/lp
   vzero = -2*zeff*aa
   if ((orbitalQNs(iorb) .lt. 0.1D0) .and. (zeff .eq. 0.D0)) &
         & vzero=vzero+ionicPotDn(lp,2)/radialPoints(2)
   if ((orbitalQNs(iorb) .gt. 0.1D0) .and. (zeff .eq. 0.D0)) &
         & vzero=vzero+ionicPotUp(lp,2)/radialPoints(2)
   if (orbitalQNs(iorb) .lt. 0.1D0) vzero=vzero+hartreePotDn2(2)
   if (orbitalQNs(iorb) .gt. 0.1D0) vzero=vzero+hartreePotUp2(2)
   var0 = 0.D0
   if (orbitalQNl(iorb) .eq. 0) var0=-2*zeff
   if (orbitalQNl(iorb) .eq. 1) var0=2.D0

   emax = v(numRadialPoints)
   emin = -1.D+20
   if (eigenValues(iorb) .gt. emax) eigenValues(iorb) = emax - 0.01D0
10 if (itmax .lt. 2) then
      write(20,*) "iorb=",iorb,"iter=",iteration,"ev=",eigenValues(iorb)*0.5,&
            & "nodes=",nodes
   endif
   if (itmax .eq. 0) return
   if (eigenValues(iorb) .ge. v(numRadialPoints)) stop 'eigenV >= v(numRadPts)'
!
!  find practical infinity ninf and classical turning
!  point nctp for orbital
!
! This loop is making a comparison between the numerical radial potential
!   function and the energy eigen value of this orbital to set up some useful
!   parameters.  It also will initialize the wave function and the derivative
!   of the wave function for r>some determined value to be zero.  This value
!   is the practical infinity.
   do while (.true.)
      nctp = numRadialPoints
      ninf = numRadialPoints
      do jj=2,numRadialPoints
         j = numRadialPoints-jj+2
         waveFn(j)    = 0.D0
         dr1WaveFn(j) = 0.D0
         if (radialPoints(j)*radialPoints(j) * &
               & (v(j)-eigenValues(iorb)) .gt. expzer) ninf=j
         if (v(j) .lt. eigenValues(iorb)) then
            exit
         else
            nctp = j
         endif
      enddo
      if (nctp > 6) then
         exit
      else
         eigenValues(iorb) = 0.9D0 * eigenValues(iorb)
      endif
   enddo
!
!  outward integration from 1 to nctp
!  startup
!
   bb = (vzero-eigenValues(iorb))/(4*lp+2)
   waveFn(1)    = 0.D0
   dr1WaveFn(1) = 0.D0
   if (orbitalQNl(iorb) .eq. 0) dr1WaveFn(1) = b*a
   do j=2,5
      waveFn(j) = radialPoints(j)**lp * &
            & (1+(aa+bb*radialPoints(j))*radialPoints(j))
      dr1WaveFn(j) = radialPoints2(j)*radialPoints(j)**orbitalQNl(iorb) * &
            & (lp+(aa*(lp+1)+bb*(lp+2)*radialPoints(j))*radialPoints(j))
   enddo
   fa5 =   dr1WaveFn(1)
   fb5 = b*dr1WaveFn(1) + radialPoints2(1)*radialPoints2(1) * var0
   fa4 =   dr1WaveFn(2)
   fb4 = b*dr1WaveFn(2) + radialPoints2(2)*radialPoints2(2) * &
         & (v(2)-eigenValues(iorb))*waveFn(2)
   fa3 =   dr1WaveFn(3)
   fb3 = b*dr1WaveFn(3) + radialPoints2(3)*radialPoints2(3) * &
         & (v(3)-eigenValues(iorb))*waveFn(3)
   fa2 =   dr1WaveFn(4)
   fb2 = b*dr1WaveFn(4) + radialPoints2(4)*radialPoints2(4) * &
         & (v(4)-eigenValues(iorb))*waveFn(4)
   fa1 =   dr1WaveFn(5)
   fb1 = b*dr1WaveFn(5) + radialPoints2(5)*radialPoints2(5) * &
         & (v(5)-eigenValues(iorb))*waveFn(5)
!
!--- set underflow trap
!
   juflow = 1
   do j=2,numRadialPoints
      if(lp*abs(log(radialPoints(j))).ge.sqrt(expzer)) juflow=j
   enddo
!
!  integration loop
!
   nodes = 0
   do j=6,nctp
!
!     predictor (Adams-Bashforth)
!
      arp = waveFn(j-1) + abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5
      brp = dr1WaveFn(j-1) + abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5
      fa0 = brp
      fb0 = b*brp + radialPoints2(j)*radialPoints2(j) * &
            & (v(j)-eigenValues(iorb))*arp
!
!     corrector (Adams-Moulton)
!
      arc = waveFn(j-1)    + amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4
      brc = dr1WaveFn(j-1) + amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4
      fa5 = fa4
      fb5 = fb4
      fa4 = fa3
      fb4 = fb3
      fa3 = fa2
      fb3 = fb2
      fa2 = fa1
      fb2 = fb1
      fa1 = brc
      fb1 = b*brc + radialPoints2(j)*radialPoints2(j) * &
            & (v(j)-eigenValues(iorb))*arc
      waveFn(j)    = arc + amc0*(fa1-fa0)
      dr1WaveFn(j) = brc + amc0*(fb1-fb0)
      fa1 = dr1WaveFn(j)
      fb1 = b*dr1WaveFn(j) + radialPoints2(j)*radialPoints2(j) * &
            & (v(j)-eigenValues(iorb))*waveFn(j)
!
!     count nodes - if no underflow
!
      if (j.gt.juflow .and. waveFn(j)*waveFn(j-1).le.0) nodes=nodes+1
   enddo
!
!  end outward integration
!
!  if number of nodes correct, start inward integration
!  else modify energy stepwise and try again
!
   if (nodes .ne. orbitalQNn(iorb)-orbitalQNl(iorb)-1) goto 200
   arctp = waveFn(nctp)
   brctp = dr1WaveFn(nctp)
!
!  inward integration from ninf to nctp
!  startup
!
   do jj=1,5
      j = ninf-jj+1
      alf = v(j) - eigenValues(iorb)
      if (alf .lt. 0.D0) alf = 0.D0
      alf = sqrt(alf)
      waveFn(j) = exp(-alf*radialPoints(j))
      dr1WaveFn(j) = -radialPoints2(j)*alf*waveFn(j)
   enddo
   fa5 =   dr1WaveFn(ninf)
   fb5 = b*dr1WaveFn(ninf) + radialPoints2(ninf)*radialPoints2(ninf) * &
         & (v(ninf)-eigenValues(iorb))*waveFn(ninf)
   fa4 =   dr1WaveFn(ninf-1)
   fb4 = b*dr1WaveFn(ninf-1) + radialPoints2(ninf-1)*radialPoints2(ninf-1) * &
         & (v(ninf-1)-eigenValues(iorb))*waveFn(ninf-1)
   fa3 =   dr1WaveFn(ninf-2)
   fb3 = b*dr1WaveFn(ninf-2) + radialPoints2(ninf-2)*radialPoints2(ninf-2) * &
         & (v(ninf-2)-eigenValues(iorb))*waveFn(ninf-2)
   fa2 =   dr1WaveFn(ninf-3)
   fb2 = b*dr1WaveFn(ninf-3) + radialPoints2(ninf-3)*radialPoints2(ninf-3) * &
         & (v(ninf-3)-eigenValues(iorb))*waveFn(ninf-3)
   fa1 =   dr1WaveFn(ninf-4)
   fb1 = b*dr1WaveFn(ninf-4) + radialPoints2(ninf-4)*radialPoints2(ninf-4) * &
         & (v(ninf-4)-eigenValues(iorb))*waveFn(ninf-4)
!
!  integration loop
!
   istop = ninf-nctp
   if (istop >= 5) then
      do jj=5,istop
         j = ninf-jj
!
!        predictor (Adams-Bashforth)
!
         arp = waveFn(j+1)    - (abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5)
         brp = dr1WaveFn(j+1) - (abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5)
         fa0 = brp
         fb0 = b*brp + radialPoints2(j)*radialPoints2(j) * &
               & (v(j)-eigenValues(iorb))*arp
!
!        corrector (Adams-Moulton)
!
         arc = waveFn(j+1)    - (amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4)
         brc = dr1WaveFn(j+1) - (amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4)
         fa5 = fa4
         fb5 = fb4
         fa4 = fa3
         fb4 = fb3
         fa3 = fa2
         fb3 = fb2
         fa2 = fa1
         fb2 = fb1
         fa1 = brc
         fb1 = b*brc + radialPoints2(j)*radialPoints2(j) * &
               & (v(j)-eigenValues(iorb))*arc
         waveFn(j)    = arc - amc0*(fa1-fa0)
         dr1WaveFn(j) = brc - amc0*(fb1-fb0)
         fa1 = dr1WaveFn(j)
         fb1 = b*dr1WaveFn(j) + radialPoints2(j)*radialPoints2(j) * &
               & (v(j)-eigenValues(iorb))*waveFn(j)
      enddo
!
!     end inward integration
!
!     rescale ar and br outside nctp to match ar(nctp) from
!     outward integration
!
   endif

   factor = arctp/waveFn(nctp)
   do j=nctp,ninf
      waveFn(j)    = factor * waveFn(j)
      dr1WaveFn(j) = factor * dr1WaveFn(j)
   enddo
!
!  find normalizing factor
!
   factor = 0.D0
   ll = 4
   do j=2,ninf
      factor = factor + ll*waveFn(j)*waveFn(j)*radialPoints2(j)
      ll = 6 - ll
   enddo
!   factor = factor / 3.0_double
   factor = factor / 3
!
!  modify eigenvalue ev
!
   dev = arctp * (brctp-dr1WaveFn(nctp)) / (factor * radialPoints2(nctp))
   if (5*abs(dev) .gt. -eigenValues(iorb)) then
!      dev=sign(eigenValues(iorb),dev)/5.0_double
      dev=sign(eigenValues(iorb),dev)/5
   endif
   itmax = itmax-1
   evold = eigenValues(iorb)
   eigenValues(iorb) = eigenValues(iorb) + dev
!   if (eigenValues(iorb) .gt. emax) eigenValues(iorb) = (evold+emax)/2.0_double
!   if (eigenValues(iorb) .lt. emin) eigenValues(iorb) = (evold+emin)/2.0_double
   if (eigenValues(iorb) .gt. emax) eigenValues(iorb) = (evold+emax)/2
   if (eigenValues(iorb) .lt. emin) eigenValues(iorb) = (evold+emin)/2
   if (abs(dev) .gt. tol*(1-eigenValues(iorb))) goto 10
!
!  normalize wavefunction and change br from d(ar)/dj to d(ar)/dr
!
   factor = 1 / sqrt(factor)
   do j=1,ninf
      waveFn(j)    = factor*waveFn(j)
      dr1WaveFn(j) = factor*dr1WaveFn(j) / radialPoints2(j)
   enddo
   return
!
!  too many nodes
!  decrease ev
!
200 if (nodes .lt. orbitalQNn(iorb)-orbitalQNl(iorb)-1) goto 210
   if (eigenValues(iorb) .lt. emax) emax = eigenValues(iorb)
   eigenValues(iorb) = eigenValues(iorb) + eigenValues(iorb)/10
!   eigenValues(iorb) = eigenValues(iorb) + eigenValues(iorb)/10.0_double
   itmax = itmax-1
   goto 10
!
!  too few nodes
!  increase ev
!
210 if (eigenValues(iorb) .gt. emin) emin = eigenValues(iorb)
   eigenValues(iorb) = eigenValues(iorb) - eigenValues(iorb)/10
!   eigenValues(iorb) = eigenValues(iorb) - eigenValues(iorb)/10.0_double
   itmax = itmax-1
   goto 10

end subroutine intNonRel

!subroutine intRel(v,iorb)
!
!   ! Use necessary modules
!   use O_Kinds
!   use O_RadialGrid
!   use ExecutionData
!   use AtomData
!   use PotData
!   use EnergyData
!
!   implicit none
!
!   ! Define passed paramters.
!   real (kind=double), dimension (numRadialPoints) :: v
!   integer :: iorb
!
!
!   stop 'RELATIVISTIC INTEGRATION NOT YET IMPLEMENTED FROM OLD PROGRAM'
!
!end subroutine intRel


subroutine orbitalInfo (currentOrb)

   ! Use necessary modules.
   use O_Kinds
   use O_Constants
   use AtomData
   use O_RadialGrid
   use WaveFnData
   use PotData
   use EnergyData
   use ExecutionData

   implicit none

   ! Define passed parameters.
   integer :: currentOrb  ! The current orbital index number.

   ! Define local variables.
   integer :: i
   integer :: numZeros
   integer :: numExtrema
   integer :: intFactor
   integer :: insideZero    ! Numerical issues can make a zero 'cover' a range.
   integer :: insideExtrema ! As for the insideZero.
   integer :: lP1 ! lQN + 1
   integer :: llP1  ! lQN * (lQN + 1)
   real (kind=double) :: accumulatedCharge
   real (kind=double) :: radial90Percent
   real (kind=double) :: radial99Percent
   real (kind=double) :: fineStructFactor
   real (kind=double) :: currDr1WaveFn
   real (kind=double) :: lastDr1WaveFn
   real (kind=double), dimension(20) :: radialZero    ! Likely < 20 zeros.
   real (kind=double), dimension(20) :: radialExtrema ! Likely < 20 extrema.
   real (kind=double), dimension(20) :: extremaWaveFn ! WaveFn value at extrema
   real (kind=double), dimension(20) :: extremaDr1WaveFn ! Deriv. near extrema
   real (kind=double), allocatable, dimension(:) :: waveFnSqrd
   real (kind=double), allocatable, dimension(:) :: chargeDensity

   ! Record the orbital angular QN plus 1 and l(l+1) for easy reference.
   lP1 = orbitalQNl(currentOrb)+1
   llP1  = orbitalQNl(currentOrb) * lP1

   ! Compute a factor based on the fine structure constant.
   fineStructFactor = 2.0_double * fineStructure * 10.0_double**(-3.0_double)

   ! Compute the zeros and the extrema of this orbital.

   ! Initialize the count of the number of zeros and extrema for this orbital.
   numZeros   = 0
   numExtrema = 0

   ! Assume this value for the first radial zero point.
   radialZero(1) = 0.0_double

   ! Record the value for the derivative of the radial wave function.  In the
   !   case of a relativistic calculation, then we must also account for...
   if (doRelativistic == 0) then
      currDr1WaveFn = dr1WaveFn(2)
   else
      if (orbitalQNs(currentOrb) < 0.1_double) then
         currDr1WaveFn = lP1*waveFn(2)/radialPoints(2) + &
               & (eigenValues(currentOrb) - ionicPotDn(lP1,2) / &
               & radialPoints(2) - hartreePotDn2(2) + fineStructFactor * &
               & fineStructFactor) * dr1WaveFn(2) / fineStructFactor
      else
         currDr1WaveFn = lP1*waveFn(2)/radialPoints(2) + &
               & (eigenValues(currentOrb) - ionicPotUp(lP1,2) / &
               & radialPoints(2) - hartreePotUp2(2) + fineStructFactor * &
               & fineStructFactor) * dr1WaveFn(2) / fineStructFactor
      endif
   endif

   ! Loop over the radial points looking for extrema and the points where the
   !   wave function is zero.
   do i = 3, numRadialPoints

      ! Abort if the number of extrema is greater than the principle quantum
      !   number of this orbital minus the orbital angular momentum.  (i.e.
      !   we will only find the number of extrema equal to n-l QNs.)
      if (numExtrema >= orbitalQNn(currentOrb) - orbitalQNl(currentOrb)) then
         exit
      endif

      ! Only record a new zero if the product of the wave function at the
      !   current radial point and the previous one is negative.  (i.e. If
      !   r(x)*r(x-1) is positive then the values of r at x and x-1 are the same
      !   sign.  If it is negative then there has been a sign change and the
      !   zero is a point between the two radial points.)
      ! Note that the use of the "insideZero" will distort the true location of
      !   the zero especially if the zero range is extended.  A zero range is
      !   a sequence of two or more wave function values that are all zero.
      if (waveFn(i)*waveFn(i-1) <= 0.0_double) then

         ! Make sure that we are not already inside of a zero range.
         if (insideZero == 0) then

            ! Mark the fact that we are now 'inside' a zero range.
            insideZero = 1

            ! Increment the number of zeros.
            numZeros = numZeros + 1

            ! Save the location of the zero as a weighted average of the radial
            !   positions on either side of the zero point of the wave function.
            radialZero(numZeros) = (waveFn(i)   * radialPoints(i-1) - &
                                 &  waveFn(i-1) * radialPoints(i)) / &
                                 & (waveFn(i) - waveFn(i-1))
         endif
      else
         ! We are not in a zero range.
         insideZero = 0
      endif

      ! Shift the previous value for the derivative of the wave function to the
      !   "last" variable.  Then, compute the new value for the current radial
      !   point.  Note again that the relativistic computation is different.
      lastDr1WaveFn = currDr1WaveFn
      if (doRelativistic == 0) then
         currDr1WaveFn = dr1WaveFn(i)
      else
         if (orbitalQNs(currentOrb) < 0.1_double) then
            currDr1WaveFn = lP1*waveFn(i)/radialPoints(i) + &
                  & (eigenValues(currentOrb) - ionicPotDn(lP1,i) / &
                  & radialPoints(i) - hartreePotDn2(i) + fineStructFactor * &
                  & fineStructFactor) * dr1WaveFn(i) / fineStructFactor
         else
            currDr1WaveFn = lP1*waveFn(i)/radialPoints(i) + &
                  & (eigenValues(currentOrb) - ionicPotUp(lP1,i) / &
                  & radialPoints(i) - hartreePotUp2(i) + fineStructFactor * &
                  & fineStructFactor) * dr1WaveFn(i) / fineStructFactor
         endif
      endif

      ! Only record a new extrema if the product of the derivative of the wave
      !   function at the current radial point and the previous one is
      !   negative.
      ! Note that the use of the "insideExtrema" will distort the true location
      !   of the zero especially if the zero range is extended.  A zero range
      !   is a sequence of two or more wave function values that are all zero.
      if (currDr1WaveFn*lastDr1WaveFn <= 0.0_double) then

         ! Make sure that we are not already inside a extrema range.
         if (insideExtrema == 0) then

            ! Mark the fact that we are now inside an extrema range.
            insideExtrema = 1

            ! Increment the number of extrema.
            numExtrema = numExtrema + 1

            ! Save the location of the extrema as a weighted average of the
            !   radial positions on either side of the zero point of the wave
            !   function's derivative.
            radialExtrema(numExtrema) = (currDr1WaveFn * radialPoints(i-1) - &
                                      &  lastDr1WaveFn * radialPoints(i)) / &
                                      & (currDr1WaveFn - lastDr1WaveFn)

            ! Compute the value of the wave function at the radial extrema as
            !   the average between the two points with a correction for the
            !   curvature of the line (extrema could be a maxima or a minima).
            extremaWaveFn(numExtrema) = (waveFn(i)+waveFn(i-1))/2.0_double - &
                  & (currDr1WaveFn**2 + lastDr1WaveFn**2) * (radialPoints(i) - &
                  & radialPoints(i-1)) / (4.0_double*(currDr1WaveFn - &
                  & lastDr1WaveFn))

            ! Record the derivative of the wave function immediately after the
            !   extrema.  (Should be next to 0.)
            extremaDr1WaveFn(numExtrema) = dr1WaveFn(i)
         endif
      else
         ! We not inside an extrema range.
         insideExtrema = 0
      endif
   enddo

   ! Compute the kinetic and potential energy of this orbital.  The potential
   !   includes only the interaction with the nuclear part.

   ! Allocate space to hold the wave function squared and then compute it.
   allocate (waveFnSqrd(numRadialPoints))
   waveFnSqrd(:) = waveFn(:)*waveFn(:)

   ! Allocate space to hold the charge density and then compute it with
   !   possible relativistic effects.
   allocate (chargeDensity(numRadialPoints))
   if (doRelativistic == 0) then
      chargeDensity(:) = waveFnSqrd(:)
   else
      chargeDensity(:) = waveFnSqrd(:) + dr1WaveFn(:)
   endif

   ! Compute the radial value for when the accumulated charge density equals
   !   90% and 99% of the total (1.0).  Note that we count backwards from
   !   numRadialPoints to improve numerical accuracy because the outer values
   !   are much smaller and would be swamped by the larger values near the
   !   nucleus.  (When adding a large and a small number, the small one is
   !   ignored due to numerical rounding on computers.  However, many many
   !   small values added together can be substantial and should not be
   !   neglected by adding them one-by-one to a large value.)
   accumulatedCharge = 0.0_double
   radial90Percent = 0.0_double
   radial99Percent = 0.0_double
   do i=numRadialPoints,2,-1

      ! Accumulate the charge.
      accumulatedCharge = accumulatedCharge + waveFnSqrd(i)*radialPoints2(i)

      ! Check if the accumulated amount is sufficient and that we have not
      !   saved the radial value yet so as to avoid repeatedly saving.
      if ((accumulatedCharge >= 0.01_double) .and. &
        & (radial99Percent == 0.0_double)) then
         radial99Percent = radialPoints(i)
      endif

      ! Check if the accumulated amount is sufficient and that we have not
      !   saved the radial value yet so as to avoid repeatedly saving.
      if ((accumulatedCharge >= 0.1_double) .and. &
        & (radial90Percent == 0.0_double)) then
         radial90Percent = radialPoints(i)

         ! We got the necessary values so we can stop here.
         exit
      endif
   enddo

   ! Compute the KE and PE for this orbital from the farthest radius in to zero.

   ! Determine the radial integral factor for the last radial point index.
   if (mod(numRadialPoints,2) == 0) then
      intFactor = 4
   else
      intFactor = 2
   endif

   ! Initialize the orbitalKE and orbitalPE.
   orbitalKE(currentOrb) = 0.0_double
   orbitalPE(currentOrb) = 0.0_double

   ! Accumulate our way in from the farthest radial point.
   do i = numRadialPoints,2,-1

      ! KE first.
      orbitalKE(currentOrb) = orbitalKE(currentOrb) + intFactor * &
            & (dr1WaveFn(i)**2 + waveFnSqrd(i)*llP1/radialPoints(i)**2) * &
            & radialPoints2(i)

      ! PE next.
      if (orbitalQNs(currentOrb) < 0.1_double) then
         orbitalPE(currentOrb) = orbitalPE(currentOrb) + intFactor * &
               & waveFnSqrd(i) * ionicPotDn(lP1,i) * radialPoints2(i) / &
               & radialPoints(i)
      else
         orbitalPE(currentOrb) = orbitalPE(currentOrb) + intFactor * &
               & waveFnSqrd(i) * ionicPotUp(lP1,i) * radialPoints2(i) / &
               & radialPoints(i)
      endif

      ! Oscillate to the next integration factor.
      intFactor = 6 - intFactor
   enddo

   ! Explicitly add record the first index of the KE (PE(1) = 0).
   orbitalKE(currentOrb) = orbitalKE(currentOrb) + dr1WaveFn(1)**2 * &
         & radialPoints2(1)

!   ! Compute the KE for this orbital.
!   orbitalKE(currentOrb) = sum(radialIntFactors(:) * &
!         & (dr1WaveFn(:)*dr1WaveFn(:) + waveFnSqrd(:) * llP1 / &
!         & (radialPoints(:)*radialPoints(:))) * radialPoints2(:))

   ! Compute the PE for this orbital from the farthest radius in to zero.

   


   ! Compute the PE for this orbital.
!   if (orbitalQNs(currentOrb) < 0.1_double) then
!      orbitalPE(currentOrb) = sum(radialIntFactors(:) * waveFnSqrd(:) * &
!            & ionicPotDn(lP1,:) * radialPoints2(:) / radialPoints(:))
!   else
!      orbitalPE(currentOrb) = sum(radialIntFactors(:) * waveFnSqrd(:) * &
!            & ionicPotUp(lP1,:) * radialPoints2(:) / radialPoints(:))
!   endif

   ! Divide by three for some reason.  ?!?!
   orbitalKE(currentOrb) = orbitalKE(currentOrb) / 3.0_double
   orbitalPE(currentOrb) = orbitalPE(currentOrb) / 3.0_double

   ! Deallocate now unnecessary arrays.
   deallocate (waveFnSqrd)
   deallocate (chargeDensity)

   ! Record the results.
   write (20,*) "------------------------------------------------------------"
   write (20,10) " n =",orbitalQNn(currentOrb),"  l =",orbitalQNl(currentOrb),&
         & "  s =",orbitalQNs(currentOrb)
   write (20,20) "ev=",eigenValues(currentOrb)*0.5_double,&
         & "  ek=",orbitalKE(currentOrb)*0.5_double,&
         & "  ep=",orbitalPE(currentOrb)*0.5_double
   write (20,30) "radialZeros   =",radialZero(:numZeros)
   write (20,30) "radialExtrema =",radialExtrema(:numExtrema)
   write (20,30) "wfextrema     =",extremaWaveFn(:numExtrema)
   write (20,30) "wfdrextrema   =",extremaDr1WaveFn(:numExtrema)
   write (20,30) "r 90%, 99%    =",radial90Percent,radial99Percent
   write (20,*) "------------------------------------------------------------"
   call flush (20)

   10 format(a4,i2,a5,i2,a5,f4.1)
   20 format(8x,a3,e15.8,a5,e15.8,a5,e15.8)
   30 format(8x,a15,2x,8f8.3)

end subroutine orbitalInfo


      SUBROUTINE TRIDIB(N,EPS1,D,E,E2,LB,UB,M11,M,W,IND,IERR,RV4,RV5)
!
      DOUBLE PRECISION EPS1,D,E,E2,W,RV4,RV5,MACHEP,LB,UB, &
     & U,V,T1,T2,XU,X0,X1
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM
      DIMENSION D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
!     REAL U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP
!     REAL ABS,AMAX1,AMIN1,FLOAT
      INTEGER IND(M)
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BISECT,
!     NUM. MATH. 9, 386-393(1967) BY BARTH, MARTIN, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971).
!
!     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
!     SYMMETRIC MATRIX BETWEEN SPECIFIED BOUNDARY INDICES,
!     USING BISECTION.
!
!     ON INPUT-
!
!        N IS THE ORDER OF THE MATRIX,
!
!        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
!          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
!          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
!          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
!          PRECISION AND THE 1-NORM OF THE SUBMATRIX,
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
!
!        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
!          E2(1) IS ARBITRARY,
!
!        M11 SPECIFIES THE LOWER BOUNDARY INDEX FOR THE DESIRED
!          EIGENVALUES,
!
!        M SPECIFIES THE NUMBER OF EIGENVALUES DESIRED.  THE UPPER
!          BOUNDARY INDEX M22 IS THEN OBTAINED AS M22=M11+M-1.
!
!     ON OUTPUT-
!
!        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
!          (LAST) DEFAULT VALUE,
!
!        D AND E ARE UNALTERED,
!
!        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
!          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
!          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
!          E2(1) IS ALSO SET TO ZERO,
!
!        LB AND UB DEFINE AN INTERVAL CONTAINING EXACTLY THE DESIRED
!          EIGENVALUES,
!
!        W CONTAINS, IN ITS FIRST M POSITIONS, THE EIGENVALUES
!          BETWEEN INDICES M11 AND M22 IN ASCENDING ORDER,
!
!        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
!          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
!          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
!          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.,
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          3*N+1      IF MULTIPLE EIGENVALUES AT INDEX M11 MAKE
!                     UNIQUE SELECTION IMPOSSIBLE,
!          3*N+2      IF MULTIPLE EIGENVALUES AT INDEX M22 MAKE
!                     UNIQUE SELECTION IMPOSSIBLE,
!
!        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
!
!     NOTE THAT SUBROUTINE TQL1, IMTQL1, OR TQLRAT IS GENERALLY FASTER
!     THAN TRIDIB, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
!
!                **********
      MACHEP = 2.0D0**(-47)
!
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = 0.0D0
!     ********** LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
!                INTERVAL CONTAINING ALL THE EIGENVALUES **********
      DO 40 I = 1, N
         X1 = U
         U = 0.0D0
         IF (I .NE. N) U = ABS(E(I+1))
         XU = DMIN1(D(I)-(X1+U),XU)
         X0 = DMAX1(D(I)+(X1+U),X0)
         IF (I .EQ. 1) GO TO 20
         IF (ABS(E(I)) .GT. MACHEP * (ABS(D(I)) + ABS(D(I-1)))) &
     &      GO TO 40
   20    E2(I) = 0.0D0
   40 CONTINUE
!
      X1 = DMAX1(ABS(XU),ABS(X0)) * MACHEP * FLOAT(N)
      XU = XU - X1
      T1 = XU
      X0 = X0 + X1
      T2 = X0
!     ********** DETERMINE AN INTERVAL CONTAINING EXACTLY
!                THE DESIRED EIGENVALUES **********
      P = 1
      Q = N
      M1 = M11 - 1
      IF (M1 .EQ. 0) GO TO 75
      ISTURM = 1
   50 V = X1
      X1 = XU + (X0 - XU) * 0.5D0
      IF (X1 .EQ. V) GO TO 980
      GO TO 320
   60 IF (S - M1) 65, 73, 70
   65 XU = X1
      GO TO 50
   70 X0 = X1
      GO TO 50
   73 XU = X1
      T1 = X1
   75 M22 = M1 + M
      IF (M22 .EQ. N) GO TO 90
      X0 = T2
      ISTURM = 2
      GO TO 50
   80 IF (S - M22) 65, 85, 70
   85 T2 = X1
   90 Q = 0
      R = 0
!     ********** ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
!                INTERVAL BY THE GERSCHGORIN BOUNDS **********
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0D0
!
      DO 120 Q = P, N
         X1 = U
         U = 0.0D0
         V = 0.0D0
         IF (Q .EQ. N) GO TO 110
         U = ABS(E(Q+1))
         V = E2(Q+1)
  110    XU = DMIN1(D(Q)-(X1+U),XU)
         X0 = DMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0D0) GO TO 140
  120 CONTINUE
!
  140 X1 = DMAX1(ABS(XU),ABS(X0)) * MACHEP
      IF (EPS1 .LE. 0.0D0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
!     ********** CHECK FOR ISOLATED ROOT WITHIN INTERVAL **********
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * FLOAT(Q-P+1)
      LB = DMAX1(T1,XU-X1)
      UB = DMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
!     ********** FIND ROOTS BY BISECTION **********
      X0 = UB
      ISTURM = 5
!
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
!     ********** LOOP FOR K-TH EIGENVALUE
!                FOR K=M2 STEP -1 UNTIL M1 DO --
!                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) **********
      K = M2
  250    XU = LB
!     ********** FOR I=K STEP -1 UNTIL M1 DO -- **********
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
!
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
!     ********** NEXT BISECTION STEP **********
  300    X1 = (XU + X0) * 0.5D0
         IF ((X0 - XU) .LE. (2.0D0 * MACHEP * &
     &      (ABS(XU) + ABS(X0)) + ABS(EPS1))) GO TO 420
!     ********** IN-LINE PROCEDURE FOR STURM SEQUENCE **********
  320    S = P - 1
         U = 1.0D0
!
         DO 340 I = P, Q
            IF (U .NE. 0.0D0) GO TO 325
            V = ABS(E(I)) / MACHEP
            IF (E2(I) .EQ. 0.0D0) V = 0.0D0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0D0) S = S + 1
  340    CONTINUE
!
         GO TO (60,80,200,220,360), ISTURM
!     ********** REFINE INTERVALS **********
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
!     ********** K-TH EIGENVALUE FOUND **********
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
!     ********** ORDER EIGENVALUES TAGGED WITH THEIR
!                SUBMATRIX ASSOCIATIONS **********
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
!
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
!
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
!
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
!
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
!     ********** SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING
!                EXACTLY THE DESIRED EIGENVALUES **********
  980 IERR = 3 * N + ISTURM
 1001 LB = T1
      UB = T2
      RETURN
!     ********** LAST CARD OF TRIDIB **********
      END SUBROUTINE
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z, &
     &                  IERR,RV1,RV2,RV3,RV4,RV6)
!
      DOUBLE PRECISION D,E,E2,W,Z,RV1,RV2,RV3,RV4,RV6, &
     &       U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,MACHEP
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      DIMENSION D(N),E(N),E2(N),W(M),Z(NM,M), &
     &       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
!     REAL U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,MACHEP
!     REAL SQRT,ABS,FLOAT
      INTEGER IND(M)
! 
!     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
!     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
!
!     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
!     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
!     USING INVERSE ITERATION.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRIX,
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
!
!        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
!          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
!          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
!          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
!          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
!          0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0
!          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
!          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
!          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE,
!
!        M IS THE NUMBER OF SPECIFIED EIGENVALUES,
!
!        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER,
!
!        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
!          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
!          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
!          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
!
!     ON OUTPUT-
!
!        ALL INPUT ARRAYS ARE UNALTERED,
!
!        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
!          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO,
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
!                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS,
!
!        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
!
!                **********
      MACHEP = 2.0D0**(-47)
!
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      TAG = 0
      ORDER = 1.0D0 - E2(1)
      Q = 0
!     ********** ESTABLISH AND PROCESS NEXT SUBMATRIX **********
  100 P = Q + 1
!
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. 0.0D0) GO TO 140
  120 CONTINUE
!     ********** FIND VECTORS BY INVERSE ITERATION **********
  140 TAG = TAG + 1
      S = 0
!
      DO 920 R = 1, M
         IF (IND(R) .NE. TAG) GO TO 920
         ITS = 1
         X1 = W(R)
         IF (S .NE. 0) GO TO 510
!     ********** CHECK FOR ISOLATED ROOT **********
         XU = 1.0D0
         IF (P .NE. Q) GO TO 490
         RV6(P) = 1.0D0
         GO TO 870
  490    NORM = ABS(D(P))
         IP = P + 1
!
         DO 500 I = IP, Q
  500    NORM = NORM + ABS(D(I)) + ABS(E(I))
!     ********** EPS2 IS THE CRITERION FOR GROUPING,
!                EPS3 REPLACES ZERO PIVOTS AND EQUAL
!                ROOTS ARE MODIFIED BY EPS3,
!                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW **********
         EPS2 = 1.0D-3 * NORM
         EPS3 = MACHEP * NORM
         UK = FLOAT(Q-P+1)
         EPS4 = UK * EPS3
         UK = EPS4 / SQRT(UK)
         S = P
  505    GROUP = 0
         GO TO 520
!     ********** LOOK FOR CLOSE OR COINCIDENT ROOTS **********
  510    IF (ABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0D0) X1 = X0 + ORDER * EPS3
!     ********** ELIMINATION WITH INTERCHANGES AND
!                INITIALIZATION OF VECTOR **********
  520    V = 0.0D0
!
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (ABS(E(I)) .LT. ABS(U)) GO TO 540
!     ********** WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
!                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY **********
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0D0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0D0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
!
         IF (U .EQ. 0.0D0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0D0
         RV3(Q) = 0.0D0
!     ********** BACK SUBSTITUTION
!                FOR I=Q STEP -1 UNTIL P DO -- **********
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
!     ********** ORTHOGONALIZE WITH RESPECT TO PREVIOUS
!                MEMBERS OF GROUP **********
         IF (GROUP .EQ. 0) GO TO 700
         J = R
!
         DO 680 JJ = 1, GROUP
  630       J = J - 1
            IF (IND(J) .NE. TAG) GO TO 630
            XU = 0.0D0
!
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
!
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
!
  680    CONTINUE
!
  700    NORM = 0.0D0
!
         DO 720 I = P, Q
  720    NORM = NORM + ABS(RV6(I))
!
         IF (NORM .GE. 1.0D0) GO TO 840
!     ********** FORWARD SUBSTITUTION **********
         IF (ITS .EQ. 5) GO TO 830
         IF (NORM .NE. 0.0D0) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
!
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
!     ********** ELIMINATION OPERATIONS ON NEXT VECTOR
!                ITERATE **********
  780    DO 820 I = IP, Q
            U = RV6(I)
!     ********** IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
!                WAS PERFORMED EARLIER IN THE
!                TRIANGULARIZATION PROCESS **********
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
!
         ITS = ITS + 1
         GO TO 600
!     ********** SET ERROR -- NON-CONVERGED EIGENVECTOR **********
  830    IERR = -R
         XU = 0.0D0
         GO TO 870
!     ********** NORMALIZE SO THAT SUM OF SQUARES IS
!                1 AND EXPAND TO FULL ORDER **********
  840    U = 0.0D0
!
         DO 860 I = P, Q
  860    U = U + RV6(I)**2
!
         XU = 1.0 / SQRT(U)
!
  870    DO 880 I = 1, N
  880    Z(I,R) = 0.0D0
!
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
!
         X0 = X1
  920 CONTINUE
!
      IF (Q .LT. N) GO TO 100
 1001 RETURN
!     ********** LAST CARD OF TINVIT **********
      END subroutine



end module ComputeOrbitalsSubs
