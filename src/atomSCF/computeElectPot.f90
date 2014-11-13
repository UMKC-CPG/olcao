module ComputeElectPotSubs

private
public :: computeElectPot

contains

subroutine computeElectPot

   ! Import the necessary data and definition modules.
   use O_Kinds
   use AtomData
   use O_RadialGrid
   use ChargeDensityData

   implicit none

   ! Compute total charge density divided by the radial value from the up+down
   !   charge density distribution of the current iteration.
   rhoOverR(1) = 0.0_double
   rhoOverR(2:) = (rhoUp(2:) + rhoDn(2:)) / radialPoints(2:)

   ! Fit the data with splines and perform the cumulative integrations.
   call splineFitAndInt

   ! Renormalize the charge density.
   call renormChargeDensity

   ! Compute the Hartree potential and (if convergence was reached) the Hartree
   !   contribution to total energy.
   call hartreePotTE

   ! Compute the exchange and correlation contributions to the total energy.
   call exchCorrTE

end subroutine computeElectPot



subroutine splineFitAndInt

   ! Import necessary data and subroutine modules.
   use O_Kinds
   use O_RadialGrid
   use ChargeDensityData
   use SplineSubs

   implicit none

   ! Define local spline associated variables.
   integer :: ierr, isx
   real (kind=double) :: a1,b1,an,bn ! Spline end condition constraints.
   real (kind=double), allocatable, dimension (:) :: work

   ! Initialize values for the spline fitting.
   isx = 0
   a1  = 0.0_double
   b1  = 0.0_double
   an  = 0.0_double
   bn  = 0.0_double

   ! Allocate space for the spline fitting work.
   allocate (work(3*numRadialPoints))

   ! Fit the charge density / radial value by splines to obtain the first and
   !   second derivatives of the curves.
   call splift(radialPoints,rhoOverR,dr1RhoOverR,dr2RhoOverR,numRadialPoints,&
         & work,ierr,isx,a1,b1,an,bn)
   if (ierr /= 1) then
      write (20,*) "Error in splift = ",ierr
   endif

   ! Deallocate space for the spline fitting work.
   deallocate (work)

   ! Compute the integral of the charge density/r from r(1) = 0 to
   !   r(numRadialPoints).
   call spliq(radialPoints,rhoOverR,dr1RhoOverR,dr2RhoOverR,numRadialPoints,&
         & 0.0_double,radialPoints,numRadialPoints,cumulInt2,ierr)
   if (ierr /= 1) then
      write (20,*) "Error in spliq call #2 = ",ierr
   endif

   ! Compute the charge density and its derivatives times the radial value.
   !   (Essentially undo the over r part.)
   dr2RhoTotal(:) = radialPoints(:)*dr2RhoOverR(:) + 2.0_double*dr1RhoOverR(:)
   dr1RhoTotal(:) = radialPoints(:)*dr1RhoOverR(:) + rhoOverR(:)
      rhoTotal(:) = radialPoints(:)*rhoOverR(:)

   ! Compute the integrals of the charge density from r(1) = 0 to
   !   r(numRadialPoints).
   call spliq(radialPoints,rhoTotal,dr1RhoTotal,dr2RhoTotal,numRadialPoints,&
         & 0.0_double,radialPoints,numRadialPoints,cumulInt1,ierr)
   if (ierr /= 1) then
      write (20,*) "Error in spliq call #1 = ",ierr
   endif
end subroutine splineFitAndInt



subroutine renormChargeDensity

   ! Import necessary modules.
   use O_Kinds
   use AtomData
   use ExecutionData
   use ChargeDensityData
   use O_RadialGrid

   implicit none

   ! Check the normalization of the charge density and compute a scaling factor
   !   if necessary.  The scaling factor is the ratio of the expected charge
   !   over the total integrated charge.  Then, apply the renormalization.
   if (electronCharge /= 0.0_double) then
      renormFactor = electronCharge / cumulInt1(numRadialPoints)

      rhoUp(:) = renormFactor * rhoUp(:)
      rhoDn(:) = renormFactor * rhoDn(:)
   endif


   ! Print a note if the charge density is renormalized by any more than 1%.
   if ((iteration > 0) .and. &
         & (abs(electronCharge-cumulInt1(numRadialPoints))) > 0.01_double) then
      write (20,*) "Warning.  Charge density rescaled in iter ", iteration,&
            & " by ",renormFactor
   endif
   
end subroutine renormChargeDensity



subroutine hartreePotTE

   ! Import necessary modules.
   use O_Kinds
   use AtomData
   use ChargeDensityData
   use O_RadialGrid
   use PotData
   use EnergyData
   use ExecutionData

   implicit none

   ! Define local variables.

   ! Allocate space for the Hartree potential as a function of r.
   if (.not. allocated (hartreePotUp)) then
      allocate (hartreePotUp(numRadialPoints))
      allocate (hartreePotDn(numRadialPoints))
      allocate (hartreePotUp2(numRadialPoints))
      allocate (hartreePotDn2(numRadialPoints))
   endif

   ! Compute the Hartree potential.
   hartreePotUp(:) = 2.0_double * renormFactor * &
         & (cumulInt1(:)/radialPoints(:) + cumulInt2(numRadialPoints) - &
         &  cumulInt2(:))
   hartreePotDn(:) = hartreePotUp(:)

   ! Compute the Hartree contribution to the total energy once the system has
   !   converged or we have reached the last iteration.
   if (converged == 1) then

      ! Sum the Hartree energy.  I think that we start from index #2 because
      !   index #1 is for the r=0 point.
      hartreeE = sum((rhoUp(2:) + rhoDn(2:)) * hartreePotDn(2:) * &
            & radialPoints2(2:) * radialIntFactors(2:))

      ! Finalize the Hartree energy with a final integration factor.
      hartreeE = hartreeE / 6.0_double
   endif
end subroutine hartreePotTE



subroutine exchCorrTE

   ! Use necessary modules.
   use O_Kinds
   use O_Constants

   ! Import necessary data modules.
   use ExecutionData
   use AtomData
   use O_RadialGrid
   use ChargeDensityData
   use PotData
   use EnergyData

   ! Import necessary subroutine modules.
   use ExchCorrSubs

   implicit none

   ! Define local variables.
   integer :: i
   real (kind=double) :: oneThird
   real (kind=double) :: fourThirds
   real (kind=double) :: two43m2    ! 2^(4/3) - 2
   real (kind=double) :: alpha0
   real (kind=double) :: xAlpha
   real (kind=double) :: beta
   real (kind=double) :: sb
   real (kind=double) :: alb


   real (kind=double) :: vxc    ! ExchCorr pot for one radial iteration.
   real (kind=double) :: vc     !     Corr pot for one radial iteration.
   real (kind=double) :: exc    ! ExchCorr energy for one radial iteration.
   real (kind=double) :: ec     !     Corr energy for one radial iteration.
   real (kind=double) :: rhoSum ! Total rho.
   real (kind=double) :: rs     ! r sub s
   real (kind=double) :: z      ! (rhoDn-RhoUp) / (rhoDn+rhoUp)
   real (kind=double) :: fz     ! f(z)
   real (kind=double) :: fzp    ! 1st derivative of f(z) wrt z.
   real (kind=double) :: vxf
   real (kind=double) :: exf
   real (kind=double) :: vcf
   real (kind=double) :: ecf
   real (kind=double) :: vxcp
   real (kind=double) :: vxcf
   real (kind=double) :: vxcd
   real (kind=double) :: vxcu
   real (kind=double) :: vcp
   real (kind=double) :: vcd
   real (kind=double) :: vcu
   real (kind=double) :: vxp
   real (kind=double) :: excp
   real (kind=double) :: excf
   real (kind=double) :: exct
   real (kind=double) :: ecp
   real (kind=double) :: ect
   real (kind=double) :: exp1 ! The 1 prevents overlap with 'exp' namespace.

   ! Initialize constants for exchange and correlation evaluation.
   oneThird   = 1.0_double / 3.0_double
   fourThirds = 4.0_double * oneThird
   two43m2    = 2.0_double**fourThirds - 2.0_double
   alpha0     = (4.0_double / (9.0_double * pi))**oneThird
   if (exchCorrCode == 0) then
      xAlpha  = 1.0_double ! For exchange only.  No correlation.
   else
      xAlpha  = 2.0_double * oneThird
   endif

   ! Initialize arrays for the exchange and correlations values.
   if (.not. allocated(vxcRadial)) then
      allocate (vxcRadial(numRadialPoints))
      allocate ( vcRadial(numRadialPoints))
      allocate (excRadial(numRadialPoints))
      allocate ( ecRadial(numRadialPoints))
   endif
   vxcRadial(:) = 0.0_double
    vcRadial(:) = 0.0_double
   excRadial(:) = 0.0_double
    ecRadial(:) = 0.0_double

   ! Initiate a loop over the radial direction to compute the exchange and
   !   correlation contributions to the total energy.
   do i = 2, numRadialPoints

      ! Compute the total charge density (after it had been normalized
      !   previously).
      rhoSum = rhoUp(i) + rhoDn(i)

      ! Cycle to the next iteration if this charge density is an unphysical
      !   number or 0.
      if (rhoSum <= 0.0_double) cycle

      ! Initialize certain parameters for exchange correlation evaluation.
      rs = (3.0_double*radialPoints(i)**2.0_double/rhoSum)**oneThird
      z   = 0.0_double
      fz  = 0.0_double
      fzp = 0.0_double

      ! Initialize parameters for spin polarized exchange correlation.
      if (doSpinPol == 1) then
         z   = (rhoDn(i) - rhoUp(i)) / rhoSum
         fz  = ((1.0_double+z)**fourThirds + (1.0_double-z)**fourThirds - &
               & 2.0_double) / two43m2
         fzp = fourThirds*((1.0_double+z)**oneThird-(1.0_double-z)**oneThird)/&
               & two43m2
      endif

      ! Compute beginning values of the exchange potential and exchange energy.
      vxp = -3.0_double*xAlpha / (pi*alpha0*rs)
      exp1=  3.0_double*vxp/4.0_double

      ! Apply modifications necessary for a relativistic calculation.
      if (doRelativistic == 1) then
         beta = 0.014_double / rs
         sb   = sqrt(1+beta*beta)
         alb  = log(beta+sb)
         vxp = vxp * (-0.5_double + 1.5_double * alb / (beta*sb))
         exp1= exp1* (1.0_double - 1.5_double * &
               & ((beta*sb-alb) / beta**2.0_double)**2.0_double)
      endif

      ! Finalize the exchange potential and exchange energy.
      vxf = 2.0_double**oneThird * vxp
      exf = 2.0_double**oneThird * exp1


      ! Initialize the correlation potential and correlation energy parameters.
      vcp = 0.0_double
      ecp = 0.0_double
      vcf = 0.0_double
      ecf = 0.0_double

      ! Call the requested subroutine for evaluation of the exchange and
      !   correlation potential and energy values.
      if     (exchCorrCode == 1) then
         call wigner(rs,vcp,ecp)
      elseif (exchCorrCode == 2) then
         call hidenLundqvist(rs,vcp,ecp)
      elseif (exchCorrCode == 3) then
         call gunnarsonLundqvistWilkins(rs,vcp,ecp,vcf,ecf)
      elseif (exchCorrCode == 4) then
         call vonBarthHedin(rs,vcp,ecp,vcf,ecf)
      elseif (exchCorrCode == 5) then
         call ceperleyAlder(rs,vcp,ecp,vcf,ecf)
      elseif (exchCorrCode == 6) then
         call vonBarthHedinGC(rs,vcp,ecp)
      endif

      ! Collect the paramagnetic exchange correlation potential and energy.
      vxcp = vxp + vcp
      excp = exp1+ ecp

      ! Collect the ferromagnetic exchange correlation potential and energy.
      vxcf = vxf + vcf
      excf = exf + ecf

      ! Initialize the up and down exchange correlation potentials.
      vxcd = vxcp
      vxcu = vxcp

      ! Initialize the up and down correlation potentials.
      vcd = vcp
      vcu = vcp

      ! Initialize the exchange correlation energy and correlation energy.
      exct = excp
      ect  = ecp

      ! Apply the effects of spin polarization (if necessary).
      if (z /= 0.0_double) then
         vxcd = vxcd + fz*(vxcf-vxcp) + (1.0_double-z)*fzp*(excf-excp)
         vxcu = vxcu + fz*(vxcf-vxcp) - (1.0_double+z)*fzp*(excf-excp)
         vcd  = vcd  + fz*(vcf-vcp)   + (1.0_double-z)*fzp*(ecf-ecp)
         vcu  = vcu  + fz*(vcf-vcp)   - (1.0_double+z)*fzp*(ecf-ecp)
         exct = exct + fz*(excf-excp)
         ect  = ect  + fz*(ecf-ecp)
      endif

      ! Collect the Hartree potential.  ???
      hartreePotDn(i) = hartreePotDn(i) + vxcd
      hartreePotUp(i) = hartreePotUp(i) + vxcu

      ! Compute the final exchange and correlation potentials and energies for
      !   this radial value.
      vxcRadial(i) = radialIntFactors(i)*(rhoDn(i)*vxcd + rhoUp(i)*vxcu) * &
            & radialPoints2(i)
      vcRadial(i)  = radialIntFactors(i)*(rhoDn(i)*vcd  + rhoUp(i)*vcu)  * &
            & radialPoints2(i)
      excRadial(i) = radialIntFactors(i)*rhoSum * exct * radialPoints2(i)
      ecRadial(i)  = radialIntFactors(i)*rhoSum * ect  * radialPoints2(i)
   enddo

   ! Accumulate the exchange and correlation energy and potential for all
   !   radial values.
   vxc = sum(vxcRadial(2:))
   vc  = sum( vcRadial(2:))
   exc = sum(excRadial(2:))
   ec  = sum( ecRadial(2:))

   ! Save the computed components of the total energy.
   totalEnergy(4) = hartreeE
   totalEnergy(5) = vxc / 3.0_double
   totalEnergy(6) = (3.0_double*vc - 4.0_double*ec) / 3.0_double
   totalEnergy(7) = exc / 3.0_double

   ! Save some aspect of the Hartree potentials in the so-far unused first
   !   index.
   hartreePotDn(1) = hartreePotDn(2) - (hartreePotDn(3) - hartreePotDn(2)) * &
         & radialPoints(2) / (radialPoints(3) - radialPoints(2))
   hartreePotUp(1) = hartreePotUp(2) - (hartreePotUp(3) - hartreePotUp(2)) * &
         & radialPoints(2) / (radialPoints(3) - radialPoints(2))
end subroutine exchCorrTE


end module ComputeElectPotSubs
