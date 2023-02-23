module O_OptcPrint

contains

subroutine printOptcResults

   ! Import necessary data modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_CommandLine,     only: stateSet
   use O_Lattice,         only: realCellVolume
   use O_KPoints,         only: numKPoints, kPointWeight
   use O_OptcTransitions, only: maxTransEnergy, energyMin, energyScale,&
                                & sumNumPartials
   use O_Input,           only: sigmaOPTC, deltaOPTC, sigmaPACS, deltaPACS,&
                                & detailCodePOPTC
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
   real (kind=double), allocatable, dimension (:,:,:,:,:) :: optcCondPOPTC


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

   ! Allocate space to hold the appropriate optical conductivity and then
   !   initialize.
   if (detailCodePOPTC == 0) then
      allocate (optcCond    (dim3,numEnergyPoints,spin))
      optcCond(:,:,:) = 0.0_double
   else
      allocate (optcCondPOPTC(sumNumPartials,sumNumPartials,dim3,&
            & numEnergyPoints,spin))
      optcCondPOPTC(:,:,:,:,:) = 0.0_double
   endif

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


   ! Output the computed results for either regular or POPTC.
   if (detailCodePOPTC == 0) then ! Regular total optical properties.

      ! Compute the optical conductivity broadened appropriately
      call getOptcCond (numEnergyPoints,optcCond,kPointFactor,sigma)

      if (stateSet == 1) then ! Doing PACS calculation.
         call printSpectrum(0,numEnergyPoints,optcCond,conversionFactor)
      else ! Do optical conductivity followed by epsilon 2.
         call printSpectrum(1,numEnergyPoints,optcCond,conversionFactor)
         call printSpectrum(2,numEnergyPoints,optcCond,conversionFactorEps2)
      endif

   else ! Partial optical properties

      ! Compute the POPTC optical conductivity broadened appropriately
      call getOptcCondPOPTC (numEnergyPoints,optcCondPOPTC,kPointFactor,sigma)

      if (stateSet == 1) then ! Doing PACS calculation.
         call printSpectrumPOPTC(0,numEnergyPoints,optcCondPOPTC,&
               & conversionFactor)
      else ! Do optical conductivity followed by epsilon 2.
         call printSpectrumPOPTC(1,numEnergyPoints,optcCondPOPTC,&
               & conversionFactor)
         call printSpectrumPOPTC(2,numEnergyPoints,optcCondPOPTC,&
               & conversionFactorEps2)
      endif
   endif

   ! Deallocate arrays.
   deallocate (energyScale)
   deallocate (kPointFactor)
   if (detailCodePOPTC == 0) then
      deallocate (optcCond)
   else
      deallocate (optcCondPOPTC)
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
               !   complete the broadening for this set because it will not
               !   have a significant effect.
               if (expAlpha < 50.0_double) then
                  optcCond(:,k,h) = optcCond(:,k,h) &
                        & + transitionProb(:,j,i,h) * exp(-expAlpha) &
                        & * kPointFactor(i)
               endif
            enddo
         enddo
      enddo
   enddo

end subroutine getOptcCond


subroutine getOptcCondPOPTC (numEnergyPoints, optcCondPOPTC, &
      & kPointFactor, sigma)

   ! Import the necessary data modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_KPoints,         only: numKPoints
   use O_OptcTransitions, only: energyScale, energyDiff, transCounter, &
         & transitionProbPOPTC

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: numEnergyPoints
   real (kind=double), dimension (:,:,:,:,:) :: optcCondPOPTC
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
               !   complete the broadening for this set because it will not
               !   have a significant effect.
               if (expAlpha < 50.0_double) then
                  optcCondPOPTC(:,:,:,k,h) = optcCondPOPTC(:,:,:,k,h) &
                        & + transitionProbPOPTC(:,:,:,j,i,h) &
                        & * exp(-expAlpha) * kPointFactor(i)
               endif
            enddo
         enddo
      enddo
   enddo

end subroutine getOptcCondPOPTC


subroutine printSpectrum (specType,numEnergyPoints,spectrum,conversionFactor)

   ! Include the modules we need.
   use O_Kinds
   use O_Potential,       only: spin
   use O_Constants,       only: hartree
   use O_OptcTransitions, only: energyScale

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: specType ! 0 = XAS; 1 = conductivity; 2 = epsilon2
   integer :: numEnergyPoints
   real (kind=double), dimension (:,:,:) :: spectrum
   real (kind=double) :: conversionFactor

   ! Define the local variables
   character*75 :: header
   integer :: h,i
   integer :: unitBase

   ! Customize the output for the current spectrum type.
   if (specType == 0) then ! XANES/ELNES
      unitBase = 49
      write (header, fmt="(5a15)") "Energy","totalXANES","xXANES",&
            & "yXANES","zXANES"
   elseif (specType == 1) then ! Optical Conductivity
      unitBase = 39
      write (header,fmt="(5a15)") "Energy","totalCond","xCond","yCond","zCond"
   elseif (specType == 2) then ! Epsilon2
      unitBase = 49
      write (header,fmt="(5a15)") "Energy","totalEps2","xEps2","yEps2","zEps2"
   endif

   ! Print the total (if spin == 1) or spin up and spin down (if spin == 2).
   do h = 1, spin

      ! Insert the appropriate header.
      write(unitBase+h,fmt="(a75)") header

      do i = 1, numEnergyPoints

         ! Adjust the spectrum for the correct units. NOTE: When this is
         !   called for printing either XANES/ELNES or the optical
         !   conductivity then the optcCond data structure is modified. For
         !   printing epsilon2, that modification side-effect is an expected
         !   prerequisite operation. 
         spectrum(:,i,h) = spectrum(:,i,h) * conversionFactor / energyScale(i)

         ! Record the spectra to disk, making sure to convert the scale to eV.
         write (unitBase+h,fmt="(5e15.7)") energyScale(i)*hartree,&
               & sum(spectrum(:,i,h))/3.0_double,spectrum(:,i,h)
      enddo
   enddo
end subroutine printSpectrum


subroutine printSpectrumPOPTC (specType,numEnergyPoints,spectrumPOPTC,&
      & conversionFactor)

   ! Include the modules we need.
   use O_Kinds
   use O_Potential,       only: spin
   use O_Constants,       only: hartree, lAngMomCount
   use O_AtomicSites,     only: numAtomSites, atomSites
   use O_AtomicTypes,     only: numAtomTypes, atomTypes
   use O_OptcTransitions, only: energyScale, sumNumPartials
   use O_Input,           only: detailCodePOPTC

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: specType ! 0 = XAS; 1 = conductivity; 2 = epsilon2
   integer :: numEnergyPoints
   real (kind=double), dimension (:,:,:,:,:) :: spectrumPOPTC
   real (kind=double) :: conversionFactor

   ! Define the local variables
   character*75 :: header
   integer :: h,i,j,k,l,m,n,o,p,q
   integer :: unitBase
   integer :: sequenceNum
   integer :: currentTypeI ! Initial
   integer :: currentTypeF ! Final
   integer :: poptcI ! pOptc initial
   integer :: poptcF ! pOptc final
   character*1, dimension (lAngMomCount) :: QN_lLetter
   character*14, dimension (4,7) :: QN_mLetter

   ! Customize the output for the current spectrum type.
   if (specType == 0) then ! XANES/ELNES
      unitBase = 149
   elseif (specType == 1) then ! Optical Conductivity
      unitBase = 139
   elseif (specType == 2) then ! Epsilon2
      unitBase = 149
   endif

   ! Define the QN_l letters.
   QN_lLetter(1) = 's'
   QN_lLetter(2) = 'p'
   QN_lLetter(3) = 'd'
   QN_lLetter(4) = 'f'

   ! Define the QN_m resolved letters.
   QN_mLetter(1,1) = 'r'
   QN_mLetter(2,1) = 'x'
   QN_mLetter(2,2) = 'y'
   QN_mLetter(2,3) = 'z'
   QN_mLetter(3,1) = 'xy'
   QN_mLetter(3,2) = 'xz'
   QN_mLetter(3,3) = 'yz'
   QN_mLetter(3,4) = 'xx~yy'
   QN_mLetter(3,5) = '2zz~xx~yy'
   QN_mLetter(4,1) = 'xyz'
   QN_mLetter(4,2) = 'xxz~yyz'
   QN_mLetter(4,3) = 'xxx~3yyx'
   QN_mLetter(4,4) = '3xxy~yyy'
   QN_mLetter(4,5) = '2zzz~3xxz~3yyz'
   QN_mLetter(4,6) = '4zzx~xxx~yyx'
   QN_mLetter(4,7) = '4zzy~xxy~yyy'


   ! Print the total (if spin == 1) or spin up and spin down (if spin == 2).
   do h = 1, spin

      ! Print the key bits of information for the pOptc output.
      write (unitBase+h,fmt="(a7)") 'STYLE 2'
      write (unitBase+h,fmt="(a10,i6)") 'NUM_UNITS ', (sumNumPartials**2)+1
      write (unitBase+h,fmt="(a11,i9)") 'NUM_POINTS ', numEnergyPoints
 
      ! Print the energy scale used by all atoms, converting to eV.
      do i = 1, numEnergyPoints
         write (unitBase+h,fmt="(f16.8)") energyScale(i) * hartree
      enddo

      ! Regardless of the decomposition style, we will always print the
      !   total spectrum first.
      sequenceNum = 1
      write (unitBase+h,fmt="(a13,i5)") 'SEQUENCE_NUM ',sequenceNum
      write (unitBase+h,fmt="(a20)") 'ELEMENT_1_NAME total'
      write (unitBase+h,fmt="(a20)") 'ELEMENT_2_NAME total'
      write (unitBase+h,fmt="(a12)") 'COL_LABELS 4'
      write (unitBase+h,fmt="(a11)") 'TOTAL x y z'

      do i = 1, numEnergyPoints

         ! Adjust the spectrum for the correct units. NOTE: When this is
         !   called for printing either XANES/ELNES or the optical
         !   conductivity then the optcCond data structure is *modified*. For
         !   printing epsilon2, that modification side-effect is an expected
         !   prerequisite operation. 
         spectrumPOPTC(:,:,:,i,h) = spectrumPOPTC(:,:,:,i,h) &
               & * conversionFactor / energyScale(i)

         ! Record the total spectrum to disk.
         write (unitBase+h,fmt="(4e15.7)") &
               & sum(spectrumPOPTC(:,:,:,i,h)) / 3.0_double,&
               & sum(spectrumPOPTC(:,:,1,i,h)), &
               & sum(spectrumPOPTC(:,:,2,i,h)), &
               & sum(spectrumPOPTC(:,:,3,i,h))
      enddo


      if (detailCodePOPTC == 1) then ! Decomposed by type.

         ! Print the Partial contributions to the optical conductivity.
         do i = 1, numAtomTypes
            do k  = 1, numAtomTypes

               sequenceNum = sequenceNum + 1

               ! Print the total for each pair
               write (unitBase+h,fmt="(a13,i5)") 'SEQUENCE_NUM ',sequenceNum
               write (unitBase+h,fmt="(a15,a3)") 'ELEMENT_1_NAME ',&
                     & atomTypes(i)%elementName
               write (unitBase+h,fmt="(a15,a3)") 'ELEMENT_2_NAME ',&
                     & atomTypes(k)%elementName
               write (unitBase+h,fmt="(a12)") 'COL_LABELS 4'
               write (unitBase+h,fmt="(a11)") 'TOTAL x y z'


               do l = 1, numEnergyPoints
                  
                  ! Record the optical conductivity to disk, making sure to
                  !   use au instead of eV.
                  write (unitBase+h,fmt="(4e15.7)") &
                        & sum(spectrumPOPTC(i,k,:,l,h)) / 3.0_double,&
                        & spectrumPOPTC(i,k,:,l,h)
               enddo
            enddo
         enddo

      elseif (detailCodePOPTC == 2) then ! Decompose by atom.

         ! Print the partial contributions.
         do i = 1, numAtomSites

            ! Obtain the type of the current initial state atom.
            currentTypeI = atomSites(i)%atomTypeAssn

            do k  = 1, numAtomSites

               ! Obtain the type of the current final state atom.
               currentTypeF = atomSites(k)%atomTypeAssn

               sequenceNum = sequenceNum + 1

               ! Print the total for each pair.
               write (unitBase+h,fmt="(a13,i5)") 'SEQUENCE_NUM ',sequenceNum
               write (unitBase+h,fmt="(a15,a3)") 'ELEMENT_1_NAME ',&
                     & atomTypes(currentTypeI)%elementName
               write (unitBase+h,fmt="(a15,a3)") 'ELEMENT_2_NAME ',&
                     & atomTypes(currentTypeF)%elementName
               write (unitBase+h,fmt="(a12)") 'COL_LABELS 4'
               write (unitBase+h,fmt="(a11)") 'TOTAL x y z'


               do l = 1, numEnergyPoints

                  ! Record the spectrum to disk.
                  write (unitBase+h,fmt="(4e15.7)") &
                        & sum(spectrumPOPTC(i,k,:,l,h)) / 3.0_double,&
                        & spectrumPOPTC(i,k,:,l,h)
               enddo
            enddo
         enddo

      elseif (detailCodePOPTC == 3) then ! Decompose by atom and QN_nl

         ! Print the partial contributions to the spectrum.
         poptci = 0
         do i = 1, numAtomTypes
            do j = 1, lAngMomCount  ! 1=s; 2=p; 3=d; 4=f
               do k = 1, atomTypes(i)%numQN_lValeRadialFns(j)
                  poptci = poptci + 1
                  poptcf = 0
                  do l = 1, numAtomTypes
                     do m = 1, lAngMomCount  ! 1=s; 2=p; 3=d; 4=f
                        do n = 1, atomTypes(l)%numQN_lValeRadialFns(m)
                           poptcf = poptcf + 1

                           sequenceNum = sequenceNum + 1
                   

                           ! Print the total for each pair
                           write (unitBase+h,fmt="(a13,i5)") 'SEQUENCE_NUM ',&
                                 & sequenceNum
                           write (unitBase+h,fmt="(a15,a,a1,i1,a1)") &
                                 & 'ELEMENT_1_NAME ', &
                                 & trim(atomTypes(i)%elementName), &
                                 & "_",atomTypes(i)%numQN_lCoreRadialFns(j) &
                                 & + k + j - 1, QN_lLetter(j)
                           write (unitBase+h,fmt="(a15,a,a1,i1,a1)") &
                                 & 'ELEMENT_2_NAME ', &
                                 & trim(atomTypes(l)%elementName), &
                                 & "_",atomTypes(l)%numQN_lCoreRadialFns(m) &
                                 & + n + m - 1, QN_lLetter(m)
                           write (unitBase+h,fmt="(a12)") 'COL_LABELS 4'
                           write (unitBase+h,fmt="(a11)") 'TOTAL x y z'

                           do o = 1, numEnergyPoints

                              ! Record the spectrum to disk.
                              write (unitBase+h,fmt="(4e15.7)") &
                                    & sum(spectrumPOPTC(poptci,poptcf,:,o,h))&
                                    & / 3.0_double, &
                                    & spectrumPOPTC(poptci,poptcf,:,o,h)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

      elseif (detailCodePOPTC == 4) then ! Decompose by atom and QN_nlm

         ! Print the Partial contributions to the optical conductivity.
         poptci = 0

         do i = 1, numAtomSites

            ! Obtain the type of the current initial state atom.
            currentTypeI = atomSites(i)%atomTypeAssn

            do j = 1, lAngMomCount  ! 1=s; 2=p; 3=d; 4=f
               do k = 1, atomTypes(currentTypeI)%numQN_lValeRadialFns(j)
                  do l = 1, (j-1)*2+1
                     poptci = poptci +1
                     poptcf = 0
                     do m = 1, numAtomSites

                        ! Obtain the type of the current final state atom.
                        currentTypeF = atomSites(m)%atomTypeAssn
                        do n = 1,lAngMomCount  ! 1=s; 2=p; 3=d; 4=f
                           do o = 1,&
                              & atomTypes(currentTypeF)%numQN_lValeRadialFns(n)
                              do p = 1, (n-1)*2+1

                                 poptcf = poptcf + 1
                                 sequenceNum = sequenceNum + 1
 
                                 ! Print the total for each pair
                                 write (unitBase+h,fmt="(a13,i5)") &
                                      & 'SEQUENCE_NUM ', sequenceNum
                                 write (unitBase+h,fmt="(a15,a,i1,a1,a1,a)") &
                                      & 'ELEMENT_1_NAME ', trim(&
                                      & atomTypes(currentTypeI)%elementName),&
                                      & atomTypes(currentTypeI)%&
                                      & numQN_lCoreRadialFns(j)+k+j-1, &
                                      & QN_lLetter(j),"_",QN_mLetter(j,l)
                                 write (unitBase+h,fmt="(a15,a,i1,a1,a1,a)")&
                                      & 'ELEMENT_2_NAME ', trim(&
                                      & atomTypes(currentTypeF)%elementName),&
                                      & atomTypes(currentTypeF)%&
                                      & numQN_lCoreRadialFns(n)+o+n-1, &
                                      & QN_lLetter(n),"_",QN_mLetter(n,p)
                                 write (unitBase+h,fmt="(a12)") 'COL_LABELS 4'
                                 write (unitBase+h,fmt="(a11)") 'TOTAL x y z'


                                 do q = 1, numEnergyPoints

                                    ! Record the spectrum to disk.
                                    write (unitBase+h,fmt="(4e15.7)") &
                                           & sum(spectrumPOPTC(poptci, &
                                           & poptcf,:,q,h)) / 3.0_double, &
                                           & spectrumPOPTC(poptci, &
                                           & poptcf,:,q,h)
                                 enddo ! q
                              enddo ! p
                           enddo ! o
                        enddo ! n
                     enddo ! m
                  enddo ! l
               enddo ! k
            enddo ! j
         enddo ! i
      endif
   enddo ! h

end subroutine printSpectrumPOPTC

!subroutine printXAS (numEnergyPoints,optcCond,conversionFactor)
!
!   ! Include the modules we need.
!   use O_Kinds
!   use O_Potential,       only: spin
!   use O_Constants,       only: hartree
!   use O_OptcTransitions, only: energyScale
!
!   ! Make sure that there are not accidental variable declarations.
!   implicit none
!
!   ! Define the dummy variables passed to this subroutine.
!   integer :: numEnergyPoints
!   real (kind=double), dimension (:,:,:) :: optcCond
!   real (kind=double) :: conversionFactor
!
!   ! Define the local variables
!   integer :: h,i
!
!   do h = 1, spin
!
!     ! Insert the header for the XANES/ELNES data.
!     write (49+h,fmt="(5a15)") "Energy","totalXANES","xXANES",&
!           & "yXANES","zXANES"
!
!      do i = 1, numEnergyPoints
!
!         ! Adjust the XANES/ELNES spectra for the correct units.
!         optcCond(:,i,h) = optcCond(:,i,h) * conversionFactor / energyScale(i)
!
!         ! Record the spectra to disk, making sure to convert the scale to au.
!         write (49+h,fmt="(5e15.7)") energyScale(i)*hartree,&
!               & sum(optcCond(:,i,h))/3.0_double,optcCond(:,i,h)
!      enddo
!   enddo
!
!end subroutine printXAS


!subroutine printCond (numEnergyPoints, optcCond, conversionFactor)
!
!   ! Include the modules we need
!   use O_Kinds
!   use O_Potential,       only: spin
!   use O_Constants,       only: hartree
!   use O_OptcTransitions, only: energyScale
!
!   ! Define the dummy variables passed to this subroutine.
!   integer :: numEnergyPoints
!   real (kind=double), dimension (:,:,:) :: optcCond
!   real (kind=double) :: conversionFactor
!
!   ! Define the local variables
!   integer :: h,i
!
!   do h = 1, spin
!
!      ! Insert the header for the optical conductivity data.
!      write (39+h,fmt="(5a15)") "Energy","totalCond","xCond","yCond","zCond"
!
!      do i = 1, numEnergyPoints
!
!         ! Adjust the optical conductivity for the correct units.
!         optcCond(:,i,h) = optcCond(:,i,h) * conversionFactor / energyScale(i)
!
!         ! Record the optical conductivity to disk, making sure to use au
!         !   instead of eV.
!         write (39+h,fmt="(5e15.7)") energyScale(i)*hartree,&
!               & sum(optcCond(:,i,h))/3.0_double,optcCond(:,i,h)
!      enddo
!   enddo
!
!end subroutine printCond


!subroutine printEps2 (numEnergyPoints, optcCond, conversionFactorEps2)
!
!   ! Include the modules we need
!   use O_Kinds
!   use O_Potential,       only: spin
!   use O_OptcTransitions, only: energyScale
!   use O_Constants,       only: dim3, hartree
!
!   ! Define the dummy variables passed to this subroutine.
!   integer :: numEnergyPoints
!   real (kind=double), dimension (:,:,:) :: optcCond
!   real (kind=double) :: conversionFactorEps2
!
!   ! Define the local variables
!   integer :: h,i
!   real (kind=double), dimension(dim3) :: epsilon2
!   real (kind=double) :: totalEpsilon2
!
!   do h = 1, spin
!
!      ! Insert the header for the epsilon 2 data.
!      write (49+h,fmt="(5a15)") "Energy","totalEps2","xEps2","yEps2","zEps2"
!
!      do i = 1, numEnergyPoints
!
!         ! Compute epsilon2 for x, y, z, and the total.
!         epsilon2(:) = optcCond(:,i,h) * conversionFactorEps2 / energyScale(i)
!         totalEpsilon2 = sum(epsilon2(:)) / 3.0_double
!
!         ! Record the value of epsilon2 to disk, making sure to use au instead
!         !   of eV for the energy scale.
!         write (49+h,fmt="(5e15.7)") energyScale(i)*hartree,&
!               & totalEpsilon2,epsilon2(:)
!
!      enddo
!   enddo
!
!end subroutine printEps2

end module O_OptcPrint
