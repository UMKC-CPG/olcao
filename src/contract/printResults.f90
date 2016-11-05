module PrintResultsSubs

public

   ! Define the module variables.
   integer, allocatable, dimension(:) :: orbIsInverted ! This is used to
         ! keep track of the signs of the orbital basis functions so that the
         ! POVRay scenes all have consistent colors for positive and negative
         ! values of the basis function.

contains

subroutine printWaveFns

   ! Use necessary definition modules.
   use O_Kinds

   ! Use necessary subroutine modules.
   use O_WriteDataSubs

   ! Use necessary data modules.
   use AtomData
   use GaussianBasisData
   use MatrixElementData

   ! Demand all variables be declared.
   implicit none

   ! Define local variables.
   integer :: i,j
   integer :: fileUnit
   integer, dimension (2) :: orbitalIDNumbers
   integer, dimension (3) :: numOrbitalsPerBasis
   integer, dimension (5) :: componentIDNumbers

   ! Define the file unit to use for output.
   fileUnit = 70

   ! Open the wave function file for printing.
   open(unit=fileUnit,form="formatted",file="waveFn.dat",status="new")

   ! Write the number of exponential alphas to use for each of the angular
   !   momentum orbital types.
   call writeLabel("NUM_ALPHA_S_P_D_F",fileUnit)
   call writeData(size(numBasisGaussians),numBasisGaussians(:),fileUnit)

   ! Write the complete set of exponential alphas.
   call writeLabel("ALPHAS",fileUnit)
   call writeData(maxNumBasisGaussians,basisAlphas,fileUnit)

   ! Write the number of core orbitals for each basis.  (All basis types have
   !   the same core.)
   numOrbitalsPerBasis(:) = sum(numCoreOrbitals(:))
   call writeLabel("NUM_CORE_RADIAL_FNS",fileUnit)
   call writeData(3,numOrbitalsPerBasis(:),fileUnit)

   ! The values in the orbitalIDNumbers are:
   !   (1)  Number of components for this n,l orbital.  Non-relativistic
   !        orbitals defined by an n,l pair have only one component.
   !        Relativistic orbitals will have one or two components.  (e.g. s =
   !        s 1/2; p = p 1/2, p 3/2; d = d 3/2, d 5/2; etc...)
   !   (2)  Basis tag.  1=MB,FB,EB; 2=FB,EB; 3=EB

   ! The values in the componentIDNumbers are:
   !   (1)  Principle quantum number (n).
   !   (2)  Angular quantum number (l).
   !   (3)  Spin+Angular quantum number (j) * 2.
   !   (4)  Number of electron states in this orbital component (not spin
   !        degenerate).
   !   (5)  Index number for component.

   if (sum(numCoreOrbitals(:)) > 0) then
      ! Write the core orbital gaussian expansion coefficients.
      call writeLabel("NL_RADIAL_FUNCTIONS",fileUnit)
      do i = 1, 4 ! s,p,d,f   ! QN_l loop
         do j = 1, numCoreOrbitals(i)   ! QN_n loop
            orbitalIDNumbers(1) = 1     ! Always 1 component for non-rel calc.
            orbitalIDNumbers(2) = 1     ! Core orbitals are used in all basis.
            call writeData(2,orbitalIDNumbers(:),fileUnit)
            componentIDNumbers(1) = j+(i-1) ! QN_n
            componentIDNumbers(2) = i-1     ! QN_l
            componentIDNumbers(3) = 2*(i-1) ! QN_j*2=QN_l*2 for non-rel calcs.
            componentIDNumbers(4) = (2*i-1)*2 ! Num of states in this orbital.
            componentIDNumbers(5) = 1       ! Only 1 component for non-rel.
            call writeData(5,componentIDNumbers(:),fileUnit)
            call writeData(numBasisTerms(i),&
                  & hamiltonianME(:numBasisTerms(i),j,i),fileUnit)
         enddo
      enddo
   endif

   ! Write the number of valence orbitals.
   numOrbitalsPerBasis(1) = sum(numValeOrbitalsPerBasis(:,1))
   numOrbitalsPerBasis(2) = sum(numValeOrbitalsPerBasis(:,1:2))
   numOrbitalsPerBasis(3) = sum(numValeOrbitalsPerBasis(:,1:3))
   call writeLabel("NUM_VALE_RADIAL_FNS",fileUnit)
   call writeData(3,numOrbitalsPerBasis(:),fileUnit)

   ! Write the valence orbital gaussian expansion coefficients.
   call writeLabel("NL_RADIAL_FUNCTIONS",fileUnit)
   do i = 1, 4 ! s,p,d,f
      do j = 1, numValeOrbitals(i)
         orbitalIDNumbers(1) = 1     ! Always only 1 component for non-rel calc.
         orbitalIDNumbers(2) = basisID(j,i) ! Basis sets this orbital is in.
         call writeData(2,orbitalIDNumbers(:),fileUnit)
         componentIDNumbers(1) = j+numCoreOrbitals(i)+(i-1) ! QN_n
         componentIDNumbers(2) = i-1     ! QN_l
         componentIDNumbers(3) = 2*(i-1) ! QN_j * 2 = QN_l*2 for non-rel calcs.
         componentIDNumbers(4) = (2*i-1)*2 ! Number of states in this orbital.
         componentIDNumbers(5) = 1       ! First and only component for non-rel.
         call writeData(5,componentIDNumbers(:),fileUnit)
         call writeData(numBasisTerms(i),&
               & hamiltonianME(:numBasisTerms(i),j+numCoreOrbitals(i),i),&
               & fileUnit)
      enddo
   enddo


   ! Close the wave function file.
   close(70)

end subroutine printWaveFns


subroutine printWaveFnsNumerical

   ! Use necessary definition modules.
   use O_Kinds
   use O_RadialGrid
   use O_Constants

   ! Use necessary subroutine modules.
   use O_WriteDataSubs

   ! Use necessary data modules.
   use AtomData
   use GaussianBasisData
   use MatrixElementData
   use ExecutionData

   ! Demand all variables be declared.
   implicit none

   ! Define local variables.
   integer :: i,j,k
   integer :: fileUnit
   integer :: columnCount
   character*1, dimension (4) :: QN_lLabel
   character*5, allocatable, dimension (:) :: columnLabels
   real (kind=double), allocatable, dimension (:,:) :: alphasNumerical
   real (kind=double), allocatable, dimension (:,:) :: waveFnsNumerical

   ! Define the file unit to use for output.
   fileUnit = 71

   ! Open the wave function file for printing.
   open(unit=fileUnit,form="formatted",file="waveFn.plot",status="new")

   ! Initialize the radial grid with hard coded values if it has not been done
   !   already due to an optional charge plot request.
   if (doCharge /= 1) then
      radialMaxDist = 80.0_double
      aaWhatever    =  5.0_double
      bbWhatever    = 50.0_double
      call setupRadialGrid(atomicNumber)

      ! Convert the radial grid from atomic units to angstroms.
      radialPoints(:) = radialPoints(:) * bohrRad
   endif

   ! Allocate space to hold the numerical evaluation of each Gaussian (alpha)
   !   on the grid.
   allocate (alphasNumerical(numRadialPoints, numBasisGaussians(1)))

   ! Allocate space to hold the numerical radial grid and the numerical wave
   !   functions of each orbital.
   allocate (waveFnsNumerical(numRadialPoints,&
         & sum(numCoreOrbitals(:)) + sum(numValeOrbitalsPerBasis(:,1:3)) + 1))

   ! Allocate space to hold a label for each column.
   allocate (columnLabels(sum(numCoreOrbitals(:)) + &
         & sum(numValeOrbitalsPerBasis(:,1:3)) + 1))

   ! Initialize the possible orbital labels for each QN_l.
   QN_lLabel(1) = "s"
   QN_lLabel(2) = "p"
   QN_lLabel(3) = "d"
   QN_lLabel(4) = "f"

   ! Define the label for the radial position.
   columnLabels(1) = "R    "

   ! Define the labels for the orbitals.
   columnCount = 1 ! The radial positions is #1.
   do i = 1, 4
      do j = 1, numCoreOrbitals(i)
         columnCount = columnCount + 1
         write(columnLabels(columnCount),fmt="(i1,a1)") j+(i-1),&
               & QN_lLabel(i)
      enddo
      do j = 1, numValeOrbitals(i)
         columnCount = columnCount + 1
         write(columnLabels(columnCount),fmt="(i1,a1)") &
               & j+numCoreOrbitals(i)+(i-1), QN_lLabel(i)
      enddo
   enddo



   ! Compute the value of the gaussian functions for each exponential alpha on
   !   the radial grid.
   do i = 1, numBasisGaussians(1)
      alphasNumerical(:,i) = exp(-basisAlphas(i)*radialPoints(:)**2)
   enddo


   ! Copy the radial grid into waveFnsNumerical column #1 for easier printing.
   waveFnsNumerical(:,1) = radialPoints(:)


   ! Compute the sum of the Gaussian basis functions with the appropriate
   !   coefficients for each orbital.
   columnCount = 1 ! The radial positions is #1 so we start here to skip it.
   waveFnsNumerical(:,2:) = 0.0_double
   do i = 1, 4  ! Consider each orbital type s,p,d,f
      do j = 1, numCoreOrbitals(i)   ! Consider each core orbital of s,p,d,f
         columnCount = columnCount + 1
         do k = 1, numBasisTerms(i)  ! Consider each Gaussian basis fn
            waveFnsNumerical(:,columnCount) = waveFnsNumerical(:,columnCount)+&
                  & hamiltonianME(k,j,i) * alphasNumerical(:,k)
         enddo
      enddo
      do j = 1, numValeOrbitals(i)
         columnCount = columnCount + 1
         do k = 1, numBasisTerms(i)
            waveFnsNumerical(:,columnCount) = waveFnsNumerical(:,columnCount)+&
                  & hamiltonianME(k,j+numCoreOrbitals(i),i)*alphasNumerical(:,k)
         enddo
      enddo
   enddo

   ! Make space to record which orbitals needed to have their sign flipped so
   !   as to appear "positive" in the plot.  This is needed later in the POVRay
   !   scenes to make color assignments all consistent.
   allocate (orbIsInverted(columnCount-1))

   ! Apply corrections and modifications for clear plotting.
   do i = 2, columnCount

      ! Correct the sign so that they all start positive, and record it for the
      !   POVRay scenes produced later.  Note for the waveFnsNumerical that
      !   the angular aspect is not incorporated.
      if (waveFnsNumerical(1,i) < 0.0_double) then
         waveFnsNumerical(:,i) = waveFnsNumerical(:,i) * (-1.0_double)
         orbIsInverted(i-1) = 1
      else
         orbIsInverted(i-1) = -1
      endif

      ! Multiply the value by the radius to show nodal character.
!      waveFnsNumerical(:,i) = waveFnsNumerical(:,i) * waveFnsNumerical(:,1)
   enddo

   ! Print the column labels.
   write (71,fmt="(40(6x,a5,1x))") columnLabels(:)


   ! Print the results.
   do i = 1, numRadialPoints
      write (71,fmt="(40e12.4)") waveFnsNumerical(i,:)
   enddo

   close (71)


   ! Deallocate unused arrays.
   deallocate (alphasNumerical)
   deallocate (waveFnsNumerical)
   deallocate (columnLabels)

end subroutine printWaveFnsNumerical



subroutine printPOVRayScenes

   ! Use necessary definition modules.
   use O_Kinds
   use O_Constants

   ! Use necessary subroutine modules.
   use O_WriteDataSubs

   ! Use necessary data modules.
   use AtomData
   use GaussianBasisData
   use MatrixElementData
   use ExecutionData

   ! Demand that all variables be declared.
   implicit none

   ! Define local variables.
   integer :: i,j,k
   integer :: fileUnit
   integer :: orbitalCount ! Count of the number of different QN_nl orbitals.
         ! This is used to reference the orbIsInverted array in the correct
         ! order that it was created in.

   ! Initialize the orbitalCount to zero.
   orbitalCount = 0

   ! Loop over all the orbitals. Each orbital will be plotted in a separate
   !   file.  First, the range of QN_n s type orbitals (core and valence) will
   !   be plotted, then the range of QN_n p type orbitals will be plotted (core
   !   followed by valence).  When a p type or higher QN_l orbital is plotted
   !   all of its QN_m will be plotted in sequence.
   do i = 1, 4 ! s,p,d,f  (Covers the l quantum number).
      do j = 1, numCoreOrbitals(i) ! Covers the n quantum number.

         ! Increment the count of the orbitals.
         orbitalCount = orbitalCount + 1

         do k = 1, 2*i-1 ! Covers the m quantum number

            ! Create the file for this orbital with the appropriate file name.
            call prepPOVRayFile(i,j,k,fileUnit,1)

            ! Write the scene header information.
            call printPOVRayHeader(fileUnit)

            ! Write the orbital scene information.
            call createIsoUnion(i,j,k,fileUnit,orbitalCount)

            ! Close the file.
            close (fileUnit)

         enddo
      enddo
      do j = 1, numValeOrbitals(i) ! Covers the n quantum number.

         ! Increment the count of the orbitals.
         orbitalCount = orbitalCount + 1

         do k = 1, 2*i-1 ! Covers the m quantum number

            ! Create the file for this orbital with the appropriate file name.
            call prepPOVRayFile(i,j+numCoreOrbitals(i),k,fileUnit,2)

            ! Write the scene header information.
            call printPOVRayHeader(fileUnit)

            ! Write the orbital scene information.
            call createIsoUnion(i,j+numCoreOrbitals(i),k,fileUnit,orbitalCount)

            ! Close the file.
            close (fileUnit)

         enddo
      enddo
   enddo

   ! Deallocate the array used to record the sign of each basis function.
   deallocate (orbIsInverted)

end subroutine printPOVRayScenes




subroutine printPOVRayHeader(fileUnit)

   ! Demand that all variables be declared.
   implicit none

   ! Define passed parameters.
   integer :: fileUnit

   ! Print the header information for the scene requesting include files.
   write (fileUnit,*) '//Include necessary libraries for color, finish, ',&
         & 'and the ability to process functions.'
   write (fileUnit,*) '#include "colors.inc"'
   write (fileUnit,*) '#include "finish.inc"'
   write (fileUnit,*) '#include "functions.inc"'

   ! Print the camera definition.
   write (fileUnit,*) 
   write (fileUnit,*) '// Orthographic camera'
   write (fileUnit,*) 'camera {'
   write (fileUnit,*) '   orthographic'
   write (fileUnit,*) '   location <0, -15, 0>'
   write (fileUnit,*) '   look_at   <0, 0,  0>'
   write (fileUnit,*) '   right     1.33*x  // aspect'
   write (fileUnit,*) '   direction <0,0,10> // direction and zoom'
   write (fileUnit,*) '   angle 67 //field (overrides direction zoom)'
   write (fileUnit,*) '}'

   ! Print the lighting definition.
   write (fileUnit,*) 
   write (fileUnit,*) '// Define the light source.'
   write (fileUnit,*) 'light_source {<-100,-200,-100> colour rgb 1}'

   ! Choose the background color.
   write (fileUnit,*) 
   write (fileUnit,*) '// Define the background color.'
   write (fileUnit,*) 'background { color White }'

   ! Declare the parameters for visualizing an isosurface of the basis fns.
   write (fileUnit,*) 
   write (fileUnit,*) '// Declare object visualization parameters.'
   write (fileUnit,*) '#declare RADIUS = 6;'
   write (fileUnit,*) '#declare THRESH = -0.15;'

end subroutine printPOVRayHeader



subroutine prepPOVRayFile(i,j,k,fileUnit,coreValeCode)

   ! Use necessary definition modules.
   use O_Kinds
   use O_Constants

   ! Use necessary subroutine modules.
   use O_WriteDataSubs

   ! Use necessary data modules.
   use AtomData
   use GaussianBasisData
   use MatrixElementData
   use ExecutionData

   ! Demand that all variables be declared.
   implicit none

   ! Define the passed parameters.
   integer, intent (IN) :: i,j,k
   integer, intent (INOUT) :: fileUnit
   integer, intent (IN) :: coreValeCode

   ! Define local variables.
   character*6  :: coreVale
   character*15 :: fileName
   character*1  :: QN_lLabel
   character*3  :: QN_nlmLabel

   ! Initialize whether we are writing a core or a valence orbital.
   if (coreValeCode == 1) then
      coreVale = "-core-"
   elseif (coreValeCode == 2) then
      coreVale = "-vale-"
   endif

   ! Initialize the orbital label for this QN_nlm.
   select case (i)
   case (1)
      QN_lLabel = "s"
   case (2)
      QN_lLabel = "p"
   case (3)
      QN_lLabel = "d"
   case (4)
      QN_lLabel = "f"
   end select
   write (QN_nlmLabel,fmt="(i1,a1,i1)") j+(i-1),QN_lLabel,k

   ! Define the file unit to use for output.
   fileUnit = 72

   ! Create the output file name for the current basis function.  The
   !   name will look something like si-core-2p.pov or b-core-1s.pov.
   if (len_trim(elementName) == 1) then
      write(fileName,fmt="(a1,a6,a3,a4)") elementName,coreVale,QN_nlmLabel,&
            & ".pov"
   else
      write(fileName,fmt="(a2,a6,a3,a4)") elementName,coreVale,QN_nlmLabel,&
            & ".pov"
   endif

   ! Open the basis function file for printing.
   open (unit=fileUnit,form="formatted",file=fileName,status="new")

end subroutine prepPOVRayFile


subroutine createIsoUnion(i,j,k,fileUnit,orbitalCount)

   ! Use necessary definition modules.
   use O_Kinds
   use O_Constants

   ! Use necessary subroutine modules.
   use O_WriteDataSubs

   ! Use necessary data modules.
   use AtomData
   use GaussianBasisData
   use MatrixElementData
   use ExecutionData

   ! Demand that all variables be declared.
   implicit none

   ! Define the passed parameters.
   integer, intent (IN) :: i,j,k
   integer, intent (IN) :: fileUnit
   integer, intent (IN) :: orbitalCount

   ! Define local variables.
   integer :: l,m
   integer :: lastIteration ! Either 1 or 2 depending on whether the s-type or
         ! a greater orbital angular momentum orbital is being made.
   real (kind=double) :: normalization
   real (kind=double) :: red
   real (kind=double) :: green
   real (kind=double) :: blue
   character*3 :: orbSign
   character*8 :: orbSignNote
   character*25 :: lmType

   ! Determine the current lmType and the normalization value.
   select case (i)
   case (1)
      lmType = "1"
      normalization = sqrt(1.0_double / 4.0_double / pi)
   case (2)
      select case (k)
      case (1)
         lmType = "x"
         normalization = sqrt(3.0_double / 4.0_double / pi)
      case (2)
         lmType = "y"
         normalization = sqrt(3.0_double / 4.0_double / pi)
      case (3)
         lmType = "z"
         normalization = sqrt(3.0_double / 4.0_double / pi)
      end select
   case (3)
      select case (k)
      case (1)
         lmType = "x*y"
         normalization = sqrt(15.0_double / 4.0_double / pi)
      case (2)
         lmType = "x*z"
         normalization = sqrt(15.0_double / 4.0_double / pi)
      case (3)
         lmType = "y*z"
         normalization = sqrt(15.0_double / 4.0_double / pi)
      case (4)
         lmType = "x*x - y*y"
         normalization = sqrt(15.0_double / 16.0_double / pi)
      case (5)
         lmType = "2*z*z - x*x - y*y"
         normalization = sqrt(5.0_double / 16.0_double / pi)
      end select
   case (4)
      select case (k)
      case (1)
         lmType = "x*y*z"
         normalization = sqrt(105.0_double / 4.0_double / pi)
      case (2)
         lmType = "z*(x*x - y*y)"
         normalization = sqrt(105.0_double / 16.0_double / pi)
      case (3)
         lmType = "x*(x*x - 3*y*y)"
         normalization = sqrt(35.0_double / 32.0_double / pi)
      case (4)
         lmType = "y*(y*y - 3*x*x)"
         normalization = sqrt(35.0_double / 32.0_double / pi)
      case (5)
         lmType = "z*(2*z*z - 3*x*x - 3*y*y)"
         normalization = sqrt(7.0_double / 16.0_double / pi)
      case (6)
         lmType = "x*(4*z*z - x*x - y*y)"
         normalization = sqrt(21.0_double / 32.0_double / pi)
      case (7)
         lmType = "y*(4*z*z - x*x - y*y)"
         normalization = sqrt(21.0_double / 32.0_double / pi)
      end select

      case default
         write (6,*) "printResults.f90 subroutine: createIsoUnion"
         write (6,*) "Only 1,2,3,4 for s,p,d,f have been implemented"
         stop
   end select

   ! Create the union of the negative and positive components of the basis fn.
   write (fileUnit,*)
   write (fileUnit,fmt="(a)") 'union {'

   ! Make the positive and negative components.  If we are doing an s-type
   !   orbital, then we only need to do this once because there are only
   !   positive values for this function.  Otherwise, we make one isosurface
   !   for the positive and one for the negative.
   if (i == 1) then
      lastIteration = 1
   else
      lastIteration = 2
   endif

   ! Make the isosurface(s).
   do l = 1, lastIteration

      ! If the orbital is inverted, then just printing it as is should be used
      !   to produce a red lobe (i.e. the negative portion) and printing it
      !   with a -1* prefix should be used to print the positive portion (blue).
      !   We will always print the positive portion first because in some cases
      !   there is no negative portion (e.g. 1s orbitals).
      ! Note that there is great potential for confusion here because the
      !   threshold in the POVRay scene is a negative number.  Hence, we are
      !   making the terms negative to show a "positive" lobe.
      if (l == 1) then
         orbSignNote = "Positive"
         red     = 0.0_double
         green   = 0.0_double
         blue    = 1.0_double
         if (orbIsInverted(orbitalCount) == 1) then
            orbSign = " 1*"
         else
            orbSign = "-1*"
         endif
      else
         orbSignNote = "Negative"
         red     = 1.0_double
         green   = 0.0_double
         blue    = 0.0_double
         if (orbIsInverted(orbitalCount) == 1) then
            orbSign = "-1*"
         else
            orbSign = " 1*"
         endif
      endif

      ! Begin the isosurface and function definition sections.  Make sure to
      !   include a note about the sign of this isosurface.
      write (fileUnit,fmt="(a,a)") '   isosurface { // Sign = ',orbSignNote
      write (fileUnit,fmt="(a)") '      function {'

      do m = 1, numBasisTerms(i)
         write (fileUnit,advance="NO",fmt="(a8,a3,e12.5,a2,a,a7,e10.3,a19)") &
               & '        ',&
               & orbSign,hamiltonianME(m,j,i)*normalization,'*(',&
               & lmType,')*exp(-',basisAlphas(m),'*(x*x + y*y + z*z))'

         ! Add the "+" sign to all except the last term.
         if (m < numBasisTerms(i)) then
            write (fileUnit,fmt="(a1)") "+"
         else
            write (fileUnit,*)
         endif

      enddo

      write (fileUnit,fmt="(a)") '      }'
      write (fileUnit,fmt="(a)") '      contained_by{sphere{0,RADIUS}}'
      write (fileUnit,fmt="(a)") '      threshold THRESH'
      write (fileUnit,fmt="(a)") '      accuracy 0.001'
      write (fileUnit,fmt="(a)") '      max_gradient 4'
      write (fileUnit,fmt="(a,f5.3,a,f5.3,a,f5.3,a)") &
            & '      pigment {rgbt <',red,',',green,',',blue,',0>}'
      write (fileUnit,fmt="(a)") '      finish {phong 0.5 phong_size 10}'
      write (fileUnit,fmt="(a)") '   }'
   enddo
   write (fileUnit,fmt="(a)") '}'

end subroutine createIsoUnion

end module PrintResultsSubs
