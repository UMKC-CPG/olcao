module O_OptcPrint

contains

subroutine printOptcResults

   ! Import necessary data modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_CommandLine,     only: stateSet
   use O_Lattice,         only: realCellVolume
   use O_KPoints,         only: numKPoints, kPointWeight
   use O_OptcTransitions, only: maxTransEnergy, energyMin, energyScale
   use O_Input,           only: sigmaOPTC, deltaOPTC, sigmaPACS, deltaPACS
   use O_Constants,       only: dim3, pi, auTime, eCharge, hPlank, hartree

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables
   integer :: i ! Loop index variables
   integer :: numEnergyPoints
   real (kind=double) :: sigma
   real (kind=double) :: energyDelta
   real (kind=double) :: conversionFactor
   real (kind=double) :: conversionFactorEps2
   real (kind=double) :: sigmaSqrtPi
   real (kind=double), allocatable, dimension (:)     :: kPointFactor
   real (kind=double), allocatable, dimension (:,:,:) :: optcCond


   ! Initialize variables.
   if (stateSet == 0) then ! Standard optical properties calculation.
      sigma           = sigmaOPTC
      energyDelta     = deltaOPTC
      energyMin       = deltaOPTC ! Start as close to 0 as possible.
      numEnergyPoints = maxTransEnergy / energyDelta + 1
   else ! PACS calculation, Sigma(E) calculations never call this subroutine.
      sigma           = sigmaPACS
      energyDelta     = deltaPACS
      ! The energyMin was already determined for PACS calculations.
      numEnergyPoints = (maxTransEnergy - energyMin) / energyDelta + 1
   endif ! Sigma(E) calculations never call this subroutine.

   ! Initialize local conversion parameters

   ! Written to agree with Cohen and Chelikowsky, "Electronic Structure and
   !   Optical Properties of Semiconductors" section 4.1 equation 4.10.  Note
   !   that we first compute sigma = e2(w) * w / (4*pi)  Note that equation
   !   4.10 comes from Ehrenreich and Cohen, Phys. Rev., 115, 786, (1959).
   ! w=omega

   ! The equation that we use is:
   !   Sigma(omega) = 2 * pi * e^2 / (3 * m^2 * omega) * (1/(2*Pi)^3) *
   !      sum(ij) [Int(BZ) [delta(Ej - Ei - hbar*omega) * |Mij(k)|^2 d3k]]
   !   omega = angular frequency
   !   Mij = Momentum matrix element
   !   Note that the integration over BZ will produce a factor of (2*Pi)^3 / 
   !      Omega.  Where Omega equals the cell volume.
   ! 

   ! In a.u. hbar=1, e=1, and m=1 (Though they still have units)

   ! In a.u. we have E=energy, T=time, M=mass, L=length.

   ! The first step is to apply the terms in the coefficient that are not
   !   equal to one.  (Note that we do not divide by three here since we
   !   will perform the averaging later.  We also do not change the units of
   !   the cell volume since that will be accounted for later.)
   conversionFactor = (2.0_double*pi) / realCellVolume

   ! The next step is to note that our calculation uses energy instead of
   !   angular frequency (omega), and that our energy is in eV.  We must
   !   convert it back from eV to a.u., then we use the fact that hbar in a.u.
   !   equals 1 to relate energy to frequency.  Since we will divide by eV we
   !   convert it to a.u. by multiplying by the hartree factor.
!   conversionFactor = conversionFactor * hartree

   ! Along a similar vein, the delta term (delta(Ej - Ei - hbar*omega)) has
   !   units of inverse energy, but this too is in eV and must be converted
   !   back to a.u. by another multiplication by the hartree factor.
!   conversionFactor = conversionFactor * hartree

   ! The next step is to do dimensional analysis and understand the unit
   !   conversion.  The equation has the following units in a.u.:
   
   ! Sigma = EL / (M^2 * T^-1) * 1/L^3 * 1/E * (ML/T)^2
   ! Sigma = 1/T
   ! Charge^2 (Q^2) is ML^3/T^2 = EL.  (In cgs charge is sqrt(g cm^3 / s^2).)
   ! Convert 1/T in a.u. to 1/s is cgs by applying the conversion factor of
   !   2.418884326505x10-17 s = 1 a.u. of time (taken from NIST).  Note that
   !   the factor of 1e-17 is not included in the constant auTime and it must
   !   be accounted for.
   conversionFactor = conversionFactor/(auTime)

   ! The last step for the conductivity is to adjust the result to have the
   !   appropriate order of magnitude.  The result should be in units of
   !   1e15 * 1/sec.  The 1/1e-17 from auTime gives us 1e17 1/s.  To put the
   !   result in the right units we must multiply the answer by 100.
   conversionFactor = conversionFactor * 100.0_double

   ! We also want to produce a result in terms of the unitless epsilon 2.  To
   !   do this we will multiply the conductivity by 4pi and divide by the
   !   frequency.  This means that we have to multiply by hbar in eV to
   !   have the proper units since we are giving the value in eV and it
   !   needs to be a frequency in 1/s.  The 2pi from the hbar will cancel
   !   with the 4pi to leave just 2.  We also divide by eCharge to put the
   !   value of hPlank in eV s instead of Js.
   conversionFactorEps2 = 2.0_double * hPlank / eCharge

   ! Again we must adjust for the units so we will multiply the result by
   !   1x10^-15 because of the difference between 1d-34 and 1d-19 for hPlank
   !   and eCharge.  I'm not yet sure why we don't apply this factor.  CHECK!
!   conversionFactorEps2 = conversionFactorEps2 * 1.0d-15
   conversionFactorEps2 = conversionFactorEps2

   ! Used for normalization of the convoluted Gaussian
   sigmaSqrtPi = sqrt(pi) * sigma

   ! Allocate space for local arrays and matrices.
   allocate (kPointFactor (numKPoints))

   ! Begin setting up and initializing the optical conductivity parameters.

   ! Allocate space to hold the energy scale.
   allocate (energyScale (numEnergyPoints))
   allocate (optcCond    (dim3,numEnergyPoints,spin))


   ! Initialize the optical conductivity parameter
   optcCond(:,:,:) = 0.0_double

   ! Assign values to the energy range
   do i = 1, numEnergyPoints
      energyScale(i) = energyMin + energyDelta * (i-1)
   enddo

   ! Fill in factor for broadening based on kpoint weight.  The 0.5 must be
   !   included because we have already accounted for the fact that each state
   !   (in the spin non-polarized case) has two electrons and when we consider
   !   the kpoint weighting we don't want to re-multiply by 2.0.
   !   Recall that the sum(kPointWeight(:)) == 2.  In the spin polarized case
   !   the division by "spin" will divide by two because now we say that there
   !   is only one electron per state.  Also note that we must divide by
   !   hartree to have the right units for sigmaSqrtPi.
   kPointFactor(:) = kPointWeight(:) * 0.5_double / sigmaSqrtPi / hartree / &
         & real(spin,double)


   ! Compute the optical conductivity broadened appropriately
   call getOptcCond (numEnergyPoints,optcCond,kPointFactor,sigma)


   ! Output the computed results.
   if (stateSet == 1) then ! Doing PACS calculation.
      call printXAS(numEnergyPoints, optcCond, conversionFactor)
   else
      call printCond (numEnergyPoints, optcCond, conversionFactor)
      call printEps2 (numEnergyPoints, optcCond, conversionFactorEps2)
   endif

end subroutine printOptcResults


subroutine getOptcCond (numEnergyPoints, optcCond, kPointFactor, sigma)

   ! Import the necessary data modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_KPoints,         only: numKPoints
   use O_OptcTransitions, only: energyScale, energyDiff, transCounter, &
         & transitionProb

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: numEnergyPoints
   real (kind=double), dimension (:,:,:) :: optcCond
   real (kind=double), dimension (:)     :: kPointFactor
   real (kind=double) :: sigma

   ! Define local variables
   real (kind=double) :: broadenEnergyDiff
   real (kind=double) :: expAlpha
   integer :: h,i,j,k

   do h = 1, spin
      do i = 1, numKPoints
         do j = 1, transCounter(i,h)
            do k = 1, numEnergyPoints

               ! Determine the energy difference between the current transition
               !   energy and the current energy scale point.
               broadenEnergyDiff = energyDiff(j,i,h) - energyScale(k)

               ! Determine the exponential alpha for the broadening factor.
               expAlpha = broadenEnergyDiff * broadenEnergyDiff / &
                     & (sigma * sigma)

               ! If the exponential alpha is too large, then we don't have to
               !   complete the broadening for this set because it will not have
               !   a significant effect.
               if (expAlpha < 50.0_double) then
                  optcCond(:,k,h) = optcCond(:,k,h)+transitionProb(:,j,i,h) * &
                        & exp(-expAlpha) * kPointFactor(i)
               endif
            enddo
         enddo
      enddo
   enddo

end subroutine getOptcCond


subroutine printXAS (numEnergyPoints, optcCond, conversionFactor)

   ! Include the modules we need.
   use O_Kinds
   use O_Potential,       only: spin
   use O_Constants,       only: hartree
   use O_OptcTransitions, only: energyScale

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: numEnergyPoints
   real (kind=double), dimension (:,:,:) :: optcCond
   real (kind=double) :: conversionFactor

   ! Define the local variables
   integer :: h,i

   do h = 1, spin

     ! Insert the header for the XANES/ELNES data.
     write (49+h,fmt="(5a15)") "Energy","totalXANES","xXANES",&
           & "yXANES","zXANES"

      do i = 1, numEnergyPoints

         ! Adjust the XANES/ELNES spectra for the correct units.
         optcCond(:,i,h) = optcCond(:,i,h) * conversionFactor / energyScale(i)

         ! Record the spectra to disk, making sure to convert the scale to au.
         write (49+h,fmt="(5e15.7)") energyScale(i)*hartree,&
               & sum(optcCond(:,i,h))/3.0_double,optcCond(:,i,h)
      enddo
   enddo

end subroutine printXAS


subroutine printCond (numEnergyPoints, optcCond, conversionFactor)

   ! Include the modules we need
   use O_Kinds
   use O_Potential,       only: spin
   use O_Constants,       only: hartree
   use O_OptcTransitions, only: energyScale

   ! Define the dummy variables passed to this subroutine.
   integer :: numEnergyPoints
   real (kind=double), dimension (:,:,:) :: optcCond
   real (kind=double) :: conversionFactor

   ! Define the local variables
   integer :: h,i

   do h = 1, spin

      ! Insert the header for the optical conductivity data.
      write (39+h,fmt="(5a15)") "Energy","totalCond","xCond","yCond","zCond"

      do i = 1, numEnergyPoints

         ! Adjust the optical conductivity for the correct units.
         optcCond(:,i,h) = optcCond(:,i,h) * conversionFactor / energyScale(i)

         ! Record the optical conductivity to disk, making sure to use au
         !   instead of eV.
         write (39+h,fmt="(5e15.7)") energyScale(i)*hartree,&
               & sum(optcCond(:,i,h))/3.0_double,optcCond(:,i,h)
      enddo
   enddo

end subroutine printCond


subroutine printEps2 (numEnergyPoints, optcCond, conversionFactorEps2)

   ! Include the modules we need
   use O_Kinds
   use O_Potential,       only: spin
   use O_OptcTransitions, only: energyScale
   use O_Constants,       only: dim3, hartree

   ! Define the dummy variables passed to this subroutine.
   integer :: numEnergyPoints
   real (kind=double), dimension (:,:,:) :: optcCond
   real (kind=double) :: conversionFactorEps2

   ! Define the local variables
   integer :: h,i
   real (kind=double), dimension(dim3) :: epsilon2
   real (kind=double) :: totalEpsilon2

   do h = 1, spin

      ! Insert the header for the epsilon 2 data.
      write (49+h,fmt="(5a15)") "Energy","totalEps2","xEps2","yEps2","zEps2"

      do i = 1, numEnergyPoints

         ! Compute epsilon2 for x, y, z, and the total.
         epsilon2(:) = optcCond(:,i,h) * conversionFactorEps2 / energyScale(i)
         totalEpsilon2 = sum(epsilon2(:)) / 3.0_double

         ! Record the value of epsilon2 to disk, making sure to use au instead
         !   of eV for the energy scale.
         write (49+h,fmt="(5e15.7)") energyScale(i)*hartree,&
               & totalEpsilon2,epsilon2(:)

      enddo
   enddo

end subroutine printEps2

end module O_OptcPrint
