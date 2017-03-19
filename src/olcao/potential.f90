module O_Potential

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define variables for the solid state potential.
   integer :: numAlphas ! The total number of potential terms (alphas) in the
         !   system.
   integer :: potDim ! Total potential dimension determined as the sum number
         !   of alphas for all potential types.
   integer :: plusUJForm ! A number that selects the form in which atoms with a
         !   plusUJ term will be provided on input. 0 = element; 1 = sequential
         !   type; 2 = individual atom.
   integer :: numPlusUJItems ! Number of plusUJ items for which there is a
         !   Hubbard U and Hund J term that needs to be applied.
   integer :: numPlusUJAtoms ! Number of atoms for which there is a plus UJ
         !   that needs to be applied.
   integer, allocatable, dimension (:) :: plusUJItemID ! List of the ID numbers
         !   of the items that have a plus UJ term (1 to numPlusUJItems).
   integer, allocatable, dimension (:) :: plusUJAtomID ! List of the atom
         !   numbers that have a plus UJ term (1 to numPlusUJAtoms).
   integer, allocatable, dimension (:) :: plusUJAtomSize ! List of the on-site
         !   matrix sizes (5 or 7 for d or f) for the atoms that have a plus
         !   UJ term (1 to numPlusUJAtoms).
   integer, allocatable, dimension (:) :: plusUJAtomValeIndex ! List of the
         !   index numbers within the valence dimension where the orbital with
         !   a plusUJ contribution is indexed.
   real (kind=double), allocatable, dimension (:)   :: plusUJAtomGSElectrons
         ! Number of electrons in the ground state that occupy the highest d or
         !   f orbital for each atom with a UJ term.
   real (kind=double), allocatable, dimension (:)   :: potAlphas ! An ordered
         !   list of all the alphas of all types in the system.
   real (kind=double), allocatable, dimension (:)   :: intgConsts ! A list of
         !   constants of integration for the terms of each type in the system.
   real (kind=double), allocatable, dimension (:,:) :: plusUJItemValue ! The
         !   Hubbard U and Hund J values for each item with a plus UJ
         !   contribution. The first index holds the U and J values. The second
         !   runs over numPlusUJItems.
   real (kind=double), allocatable, dimension (:,:) :: plusUJAtomValue ! The
         !   Hubbard U and Hund J values for each atom with a plus UJ
         !   contribution. The first index holds the U and J values. The second
         !   runs over numPlusUJAtoms.
   real (kind=double), allocatable, dimension (:)   :: typeSpinSplit ! The
         !   initial splitting used for each type.  This will be applied to
         !   every term in potDim as the spinSplitFactor below.
   real (kind=double), allocatable, dimension (:)   :: spinSplitFactor ! A
         !   factor -1.0 to 1.0 for each potential alpha (term) that is used
         !   to create the initial spin splitting kick.  The number (x say) is
         !   used as follows:  (up-down) = (up+down)*x.
   real (kind=double), allocatable, dimension (:,:) :: potCoeffs ! An ordered
         !   list of the potential coefficients for each spin orientation.
         !   Index1=coeff; Index2=spin(1=up,2=dn)
!#ifndef GAMMA
!   complex (kind=double), allocatable, dimension (:,:,:,:) :: plusUJ
   real (kind=double), allocatable, dimension (:,:,:,:) :: plusUJ
!#else
!   real (kind=double), allocatable, dimension (:,:,:) :: plusUJGamma
!#endif
         !   Contributions to the Hamiltonian from UJ terms. The first index
         !   holds a 1:7 array of modifiers to d or f orbital terms of specific
         !   atoms. The second index runs the range from 1 to numPlusUJAtoms.
         !   The third index covers spin up and spin down. For the non-Gamma
         !   case, there is a fourth index for the number of kpoints.
         ! See: Anisimov VI, Zaanen J, Andersen OK. Band theory and Mott
         !   insulators: Hubbard U instead of Stoner I. Physical Review B,
         !   1991;44(3):943. Available from:
         !   http://dx.doi.org/10.1103/PhysRevB.44.943
         ! See: Ching WY, Gu Z, Xu Y-N. Theoretical calculation of the optical
         !   properties of Y3Fe5O12. Journal of Applied Physics, 2001 Jun 1,
         !   89(11):68835. Available from:
         !   http://jap.aip.org/resource/1/japiau/v89/i11/p6883_s1
         ! Note that there is an error in equation (1) of the Ching paper. The
         !   sigma in the first m' summation should be a -sigma indicating that
         !   the charge being referred to is associated with a spin that is
         !   opposite to the spin of the potential being computed. The Anisimov
         !   paper has the correct formulation.

   ! Define variables for controlling convergence of the potential.
   integer :: feedbackLevel
   integer :: lastIteration
   integer :: currIteration
   integer :: xcCode
   integer :: converged
   real (kind=double) :: relaxFactor
   real (kind=double) :: convgTest

   integer :: spin
   integer :: rel
   integer :: GGA

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


! Set up data structures and values associated with the total potential
!   function.
subroutine initPotStructures

   ! Use necessary modules.
   use O_Kinds
   use O_Constants, only: pi
   use O_ElementData, only: initElementData, numUJElectrons
   use O_AtomicSites, only: numAtomSites, atomSites
   use O_AtomicTypes, only: atomTypes
   use O_PotTypes, only: numPotTypes, potTypes

   implicit none

   ! Define local variables.
   integer :: i,j
   integer :: currentElementID
   integer :: currentTypeID
   integer :: currentNucCharge
   integer :: currAtomSite
   integer :: currAtomType
   real (kind=double) :: pi32  ! pi^(3/2)

   ! Initialize variables
   pi32            = pi ** 1.5_double
   potDim          = 0
   call initElementData

   ! Compute the value for potDim.
   potDim = potTypes(numPotTypes)%cumulAlphaSum+potTypes(numPotTypes)%numAlphas

   ! Allocate space for the spin split factor, integration constants, and
   !   actual alphas in one long ordered list.
   allocate (spinSplitFactor (potDim))
   allocate (intgConsts      (potDim))
   allocate (potAlphas       (potDim))

   ! Compute the constants of integration, copy the potential alphas into the
   !   total list from the individual types, and copy the spin split factor
   !   from each type to all terms of that type.
   do i = 1, numPotTypes

      ! Calculate the integration constants where the constant (j) is the
      !   multiplicity * (pi^(3/2)) / alpha(j) * alpha(j)^1/2.
      ! Also store the alphas of all the types in one array.
      do j = 1, potTypes(i)%numAlphas

         spinSplitFactor(j + potTypes(i)%cumulAlphaSum) = typeSpinSplit(i)

         potAlphas(j + potTypes(i)%cumulAlphaSum) = potTypes(i)%alphas(j)

         intgConsts(j + potTypes(i)%cumulAlphaSum) = &
               & real(potTypes(i)%multiplicity,double) * &
               & pi32 / (potTypes(i)%alphas(j))**(3.0_double/2.0_double)
      enddo
   enddo



   ! At this point we will "fix" the representation of the plusUJ factors to be
   !   completely on a per-atom basis. (That is, the input file may specify the
   !   plusUJ factors for specific elements or specific types and here we will
   !   express the plusUJ factors for each atom. The plusUJ factors here are the
   !   U and J values.)
   select case (plusUJForm)
   case (0) ! Elements
      ! For each atom, go through the list of Items and determine if there
      !   is one with the same element ID number. If so, then increment the
      !   number of numPlusUJAtoms and go to the next atom, i.
      numPlusUJAtoms = 0
      do i = 1, numAtomSites
         currentElementID = atomTypes(atomSites(i)%atomTypeAssn)%elementID
         do j = 1, numPlusUJItems
            if (currentElementID == plusUJItemID(j)) then
               numPlusUJAtoms = numPlusUJAtoms + 1
               exit
            endif
         enddo
      enddo

      ! Once counted, allocated space to hold the ID and Value data.
      allocate(plusUJAtomID(numPlusUJAtoms))
      allocate(plusUJAtomValue(2,numPlusUJAtoms))

      ! Finally, for each atom again determine if there is an Item with a
      !   matching elementID number. If so, then assign the Value and ID
      !   numbers before exiting to the next atom, i. (The weird aspect is
      !   that we reset the numPlusUJAtoms. This has to be done to keep
      !   track of the array index positions in the newly generated Value
      !   and ID arrays.)
      numPlusUJAtoms = 0
      do i = 1, numAtomSites
         currentElementID = atomTypes(atomSites(i)%atomTypeAssn)%elementID
         do j = 1, numPlusUJItems
            if (currentElementID == plusUJItemID(j)) then
               numPlusUJAtoms = numPlusUJAtoms + 1
               plusUJAtomID(numPlusUJAtoms) = i
               plusUJAtomValue(:,numPlusUJAtoms) = plusUJItemValue(:,j)
               exit
            endif
         enddo
      enddo
   case (1) ! Type
      ! For each atom, go through the list of Items and determine if there
      !   is one with the same sequential type ID number. If so, then
      !   increment the number of numPlusUJAtoms and go to the next atom, i.
      numPlusUJAtoms = 0
      do i = 1, numAtomSites
         currentTypeID = atomSites(i)%atomTypeAssn
         do j = 1, numPlusUJItems
            if (currentTypeID == plusUJItemID(j)) then
               numPlusUJAtoms = numPlusUJAtoms + 1
               exit
            endif
         enddo
      enddo

      ! Once counted, allocated space to hold the ID and Value data.
      allocate(plusUJAtomID(numPlusUJAtoms))
      allocate(plusUJAtomValue(2,numPlusUJAtoms))

      ! Finally, for each atom again determine if there is an Item with a
      !   matching sequential type ID number. If so, then assign the Value
      !   and ID numbers before exiting to the next atom, i. (The weird
      !   aspect is that we reset the numPlusUJAtoms. This has to be done to
      !   keep track of the array index positions in the newly generated
      !   Value and ID arrays.)
      numPlusUJAtoms = 0
      do i = 1, numAtomSites
         currentTypeID = atomSites(i)%atomTypeAssn
         do j = 1, numPlusUJItems
            if (currentTypeID == plusUJItemID(j)) then
               numPlusUJAtoms = numPlusUJAtoms + 1
               plusUJAtomID(numPlusUJAtoms) = i
               plusUJAtomValue(:,numPlusUJAtoms) = plusUJItemValue(:,j)
               exit
            endif
         enddo
      enddo
   case (2) ! Atoms
      ! This is the easiest because it is just a straight copy procedure.
      numPlusUJAtoms = numPlusUJItems
      allocate(plusUJAtomID(numPlusUJAtoms))
      allocate(plusUJAtomValue(2,numPlusUJAtoms))
      plusUJAtomID(:) = plusUJItemID(:)
      plusUJAtomValue(:,:) = plusUJItemValue(:,:)
   end select

   ! In every case, now that we have the atom IDs we can compute the array
   !   sizes for each atom. Note the implicit connection between atomic types
   !   and potential types in the "select" statement.
   allocate(plusUJAtomSize(numPlusUJAtoms))
   allocate(plusUJAtomGSElectrons(numPlusUJAtoms))
   do i = 1, numPlusUJAtoms

      ! Determine the number of protons in this atom. We use this number to
      !   define the highest occupied type of orbital (d or f) and thus the
      !   m quantum number for the orbital to be treated with a UJ term. The
      !   proton number (a.k.a. Z number) is also used to determine the number
      !   of electrons in the d or f orbitals in the ground state by
      !   referencing data stored in the elements.dat file.
      currentNucCharge = &
            & int(potTypes(atomSites(plusUJAtomID(i))%atomTypeAssn)%nucCharge)

      ! This determination was assigned just by looking at the periodic table
      !   of the elements. The whole process will need to be improved to deal
      !   with multiple plusUJ orbital contributions from one atom if that
      !   should ever be implemented.
      select case (currentNucCharge)
         case (21 : 57)
            plusUJAtomSize(i) = 5
         case (58 : 70)
            plusUJAtomSize(i) = 7
         case (71 : 90)
            plusUJAtomSize(i) = 5
         case (91 : 103)
            plusUJAtomSize(i) = 7
         case default
            write (20,*) "An atom has been assigned a UJ on-site term, but ",&
                  & "its Z number is not high" 
            write (20,*) "enough to contain any d or f valence orbitals."
            stop
      end select

      ! Record the number of electrons in the ground state (GS) of the atom.
      !   This data was taken from a periodic table of the elements and should
      !   probably be improved upon in some way for more complicated scenarios.
      plusUJAtomGSElectrons(i) = real(numUJElectrons(currentNucCharge),double)
   enddo

   ! Similarly, we can now also compute the index number within the valence
   !   dimension that the orbitals each atom with a plusUJ term will occupy.
   allocate(plusUJAtomValeIndex(numPlusUJAtoms))
   do i = 1, numPlusUJAtoms
      ! Get the atom number and type number of the current plusUJ atom.
      currAtomSite = plusUJAtomID(i)
      currAtomType = atomSites(currAtomSite)%atomTypeAssn

      ! Obtain the index number by starting with the accumulated valence states
      !   of all atoms before this one and then add 1 for every s state, and 3
      !   for every p state.
      plusUJAtomValeIndex(i) = atomSites(currAtomSite)%cumulValeStates + &
            & 1 * atomTypes(currAtomType)%numQN_lValeRadialFns(1) + &
            & 3 * atomTypes(currAtomType)%numQN_lValeRadialFns(2)

      ! At this point we must determine if the plusUJ term is for a d or f
      !   orbital. Presently, we can only include a plusUJ term for the highest
      !   occupied d or f state. Therefore, if there are no f orbitals the
      !   plusUJ must be for a d orbital. The first d orbital will always be
      !   the highest occupied one, any fully occupied d orbital would be in the
      !   core. The same is true if the plusUJ is for an f orbital. Therefore,
      !   once all the s and p orbitals have been accounted for, the next index
      !   in the valence dimension will be the start of the highest occupied d
      !   orbital (if the plusUJ is for a d orbital). In the case that we are
      !   applying the approach to an f orbital, then we need to add on the d
      !   orbital count. Then, the first f orbital will be the highest occupied
      !   one.
      ! So, if numQN_lRadialFns(4) /= 0 then we add on the d orbitals, otherwise
      !   we are already at the right index.
      if (atomTypes(currAtomType)%numQN_lValeRadialFns(4) /= 0) then
         plusUJAtomValeIndex(i) = plusUJAtomValeIndex(i) + &
               & 5 * (atomTypes(currAtomType)%numQN_lValeRadialFns(3))
      endif

      ! A final tricky point. The plusUJAtomValeIndex number is one less than
      !   the actual starting index so that when it is used in a loop, the loop
      !   can run from 1 to the Size for that atom and when added to the Index
      !   number it will cover the correct range.
   enddo

   ! Deallocate arrays that will not be used further.
   deallocate (typeSpinSplit)
   deallocate (plusUJItemID)
   deallocate (plusUJItemValue)

end subroutine initPotStructures

subroutine setPotControlParameters (fbL,lastIt,corrCode,rlxFact,cTest,&
      & plusUJFormTemp,numPlusUJItemsTemp,plusUJItemIDTemp,&
      & plusUJItemValueTemp,typeSpinSplitTemp)

   ! Use necessary modules.
   use O_Kinds
   use O_PotTypes, only: numPotTypes, potTypes

   implicit none

   ! Define passed dummy variables.
   integer :: fbL
   integer :: lastIt
   integer :: corrCode
   integer :: plusUJFormTemp
   integer :: numPlusUJItemsTemp
   integer, dimension(:) :: plusUJItemIDTemp
   real (kind=double) :: rlxFact
   real (kind=double) :: cTest
   real (kind=double), dimension(:,:) :: plusUJItemValueTemp
   real (kind=double), dimension(:) :: typeSpinSplitTemp

   ! Define local variables.
   integer :: i
   integer :: info
   integer :: numCodes
   integer :: xcCodeParam
   integer :: spinParam
   integer :: relParam
   integer :: GGAParam
   character(len=25) :: functionalName
   character(len=100) :: dataDirectory
   character(len=100) :: xcCodeDataFile

   ! Set potential parameters read from input.
   feedbacklevel   = fbL
   lastIteration   = lastIt
   xcCode          = corrCode
   relaxFactor     = rlxFact
   convgTest       = cTest
   plusUJForm      = plusUJFormTemp
   numPlusUJItems  = numPlusUJItemsTemp

   ! At this stage we merely need to transfer the temporary data about the
   !   set of items with some Hubbard U and Hund J contribution into variables
   !   that are held specifically in this module. The items may represent
   !   elements, types, or individual atoms. At a later point, the data will be
   !   expanded to a form where every atom will be represented.
   ! Presently the code only allows *one* pair of UJ terms per atom, so
   !   elements with some combination of d and f orbitals cannot have multiple
   !   UJ terms yet. We only treat the highest occupied d or f orbital.
   allocate (plusUJItemValue(2,numPlusUJItems))
   allocate (plusUJItemID(numPlusUJItems))
   plusUJItemValue(:,:) = plusUJItemValueTemp(:,:)
   plusUJItemID(:) = plusUJItemIDTemp(:)
   ! The size of the data (either 1:5 or 1:7) will be allocated and computed
   !   later in the initPotStructures subroutine.

   ! Allocate space to hold the initially read in spin splitting for each type.
   allocate (typeSpinSplit(numPotTypes))
   typeSpinSplit(:) = typeSpinSplitTemp(:)

   ! Initialize control parameters.
   currIteration = 1
   converged     = 0

   ! Open the xc_code.dat file
   call get_environment_variable("OLCAO_DATA", dataDirectory, info)
   xcCodeDataFile=trim(dataDirectory)//"/xc_code.dat"
   open (unit=9, file=xcCodeDataFile)

   ! Read past the header line and then get the total number of codes.
   read (9,*)
   read (9,*) numCodes

   ! Read each line(functional type). Once the correct one is found use
   ! it to set "spin", "rel", "GGA".
   do i = 1, numCodes
      read (9,*) functionalName, xcCodeParam, spinParam, relParam, GGAParam
      if (xcCodeParam.eq.xcCode) then
         spin = spinParam
         rel = relParam
         GGA = GGAParam
      endif
   enddo
   close (9)

end subroutine setPotControlParameters

subroutine initPotCoeffs

   ! Include the necessary modules
   use O_Kinds
   use O_PotTypes, only: numPotTypes, potTypes
   use O_KPoints, only: numKPoints

   ! Define local variables
   integer :: i,j,k ! Loop index variables
   integer :: potTermCount
   integer :: tempInt
   real (kind=double) :: spaceHolder0 ! Spin down coefficients
   real (kind=double) :: spaceHolder1 ! Gaussian exponential alphas
   real (kind=double) :: spaceHolder2 ! Total (core+vale (up+down)) chrg dsty.
   real (kind=double) :: spaceHolder3 ! Valence (up+down) charge density.
   real (kind=double) :: spaceHolder4 ! Valence (up-down) charge density.

   ! Allocate space for the potential coefficients and the alphas.
   allocate (potCoeffs(potDim,spin))

   ! Read the existing potential coefficients for each term, or for the
   !   spin up and then spin down terms separately if we are doing a spin
   !   polarized calculation. Note that by the time this subroutine is called
   !   we already know many of the numbers provided in this file (e.g. the
   !   number of potential types and the spin). We will ignore those quantities
   !   from the file and keep the ones from the olcao.dat input file. (Although
   !   we do make a consistency check.)
   read (8,*) ! File header that says number of types.
   do i = 1, spin

      ! Initialize the counter of potential terms.
      potTermCount = 0

      read (8,*)  ! Read tag indicating spin up or spin down.
      do j = 1, numPotTypes
         read (8,*) tempInt ! The number of terms for this type.

         if (tempInt /= potTypes(j)%numAlphas) then
            write (20,*) "Mismatch between scfV and olcao.dat:"
            write (20,*) "scfV claims type ",j," has ",tempInt," terms"
            write (20,*) "olcao.dat claims type ",j," has ",&
                  & potTypes(j)%numAlphas," terms"
            stop
         endif

         do k = 1, potTypes(j)%numAlphas

            ! Increment the counter.
            potTermCount = potTermCount + 1

            read (8,*) potCoeffs(potTermCount,i), spaceHolder1, spaceHolder2, &
                  & spaceHolder3, spaceHolder4
         enddo
      enddo
   enddo

   ! In the event that spin==1 then we need to read past the duplication of
   !   potential information to get to the +UJ terms. If spin==2, then the
   !   loop above already read in the spin down potential terms.
   if (spin == 1) then
      read (8,*) ! SPIN_DN tag
      do j = 1, numPotTypes
         read (8,*) spaceHolder1 ! Num terms data
         do k = 1, potTypes(j)%numAlphas
            read (8,*) spaceHolder0, spaceHolder1, spaceHolder2, spaceHolder3,&
                  & spaceHolder4 ! Duplicate data that we can ignore
         enddo
      enddo
   endif

   ! At this point we have reached the section that contains +UJ data (if any).

   ! Allocate space for the +UJ terms. We allocate a 1:7 space for each atom
   !   and spin even though some atoms will only fill 1:5 of the available
   !   spaces.
!#ifndef GAMMA
   allocate (plusUJ(7,numPlusUJAtoms,spin,numKPoints))
!#else
!   allocate (plusUJGamma(7,numPlusUJAtoms,spin))
!#endif

   ! Read past the number of atoms that have +UJ terms. (We already know at
   !   this point how many there are anyway.) There is a tag line and an
   !   actual number line.
   read (8,*) ! Line with a "NUM_PLUSUJ_TERMS" tag on it.
   read (8,*) tempInt ! Number of atoms with a +UJ term.
   if (tempInt /= numPlusUJAtoms) then
      write (20,*) "Mismatch between scfV and olcao.dat:"
      write (20,*) "scfV claims that there are ",tempInt," PlusUJ terms"
      write (20,*) "olcao.dat claims that there are ",numPlusUJAtoms,"."
      stop
   endif

   ! As discussed in the potentialUpdate.f90 file, this loop must go from 1 to
   !   2 regardless of the value of spin. This will preserve consistency in the
   !   structure of the data file for all calculations. As a consequence, if
   !   the value of spin==1 then the second set of data is redundent and it will
   !   be read in, but ignored.
   do j = 1, 2
      read (8,*) ! Line with a "TOTAL__OR__SPIN_UP" tag on it for i==1 or a
                 !   "SPIN_DN" tag on it for i==2.

      ! Read the 1:7 +UJ array for each atom. Note that for d-type orbitals the
      !   array elements that are "outside" the 1:5 array will all be zero.
      do k = 1, numPlusUJAtoms
         read (8,*) tempInt ! Atom ID number
         if (tempInt /= plusUJAtomID(k)) then
            write (20,*) "Mismatch between scfV and olcao.dat:"
            write (20,*) "scfV claims atom ",tempInt," is atom number ",j
            write (20,*) "olcao.dat claims that it is ", plusUJAtomID(k)
            stop
         endif
         do i = 1, numKPoints
!#ifndef GAMMA
         if (spin == 1) then
            read (8,*) plusUJ(:,k,1,i) ! See above docs about spin index.
         else
            read (8,*) plusUJ(:,k,j,i)
         endif
         enddo
!#else
!         if (spin == 1) then
!            read (8,*) plusUJGamma(:,k,1) ! See above docs about spin index.
!         else
!            read (8,*) plusUJGamma(:,k,j)
!         endif
!#endif
      enddo
   enddo

end subroutine initPotCoeffs


subroutine cleanUpPotential

   implicit none

   deallocate (spinSplitFactor)
   deallocate (intgConsts)
   deallocate (potAlphas)

   ! This is allocated for main and intg but not other programs.
   if (allocated(potCoeffs)) then
      deallocate (potCoeffs)
   endif

end subroutine cleanUpPotential


end module O_Potential
