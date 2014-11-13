module ChargeDensityDataModule

   ! Import the precision variables
   use O_Kinds

   character*14 :: gaussFile
   character*12, allocatable, dimension(:) :: meshFile
   character*12, allocatable, dimension(:) :: colLabel
   integer :: numCols
   integer :: maxNumSitesPerType
   integer :: maxNumTermsPerType
   integer :: numAtoms
   integer :: numTypes
   integer :: numAtomicPotTerms
   integer :: neighborCell
   integer :: noOpenDX
   integer :: numTotalMeshPoints
   integer, dimension (3) :: numMeshPoints
   integer, allocatable, dimension (:) :: numSitesPerType
   integer, allocatable, dimension (:) :: numTermsPerType
   real (kind=double), dimension (3)   :: fractStrideLength
   real (kind=double), dimension (3)   :: latticeMag
   real (kind=double), dimension (3)   :: latticeAngles  ! bc, ac, ab
   real (kind=double), dimension (3)   :: planeAngles    ! bc,a ; ac,b ; ab,c
   real (kind=double), dimension (3)   :: fractCrossArea ! bc, ac, ab
   real (kind=double), dimension (3,3) :: unitNormal     ! (xyz,bc ac ab)
   real (kind=double), dimension (3,3) :: latticeVectors ! (xyz,abc)
   real (kind=double), dimension (3,3) :: latticeUnitVectors ! (xyz,abc)
   real (kind=double), allocatable, dimension (:,:,:) :: profile ! Index1=a,b,c
         ! Index2 = Pot,potDiff,RhoV,RhoVDiff,RhoT;
         ! Index3 = 1..numMeshPoints(a,b,c)
   real (kind=double), allocatable, dimension (:,:,:) :: sitesData
   real (kind=double), allocatable, dimension (:,:,:) :: gaussCoeffs
         ! Index1=numCols; Index2=maxNumTermsPerType; Index3=numTypes
   real (kind=double), allocatable, dimension (:,:)   :: gaussAlphas

   contains

! This subroutine will set some default (hard coded) values for some array
!   sizes and other values.  (Basically, this will just be an easy place to
!   change values if the program is modified or extended in the future.)
subroutine initEnv

   ! Make sure that no variables are accidentally defined.
   implicit none

   numCols = 5

   allocate (meshFile(numCols))
   allocate (colLabel(numCols))
   meshFile(1) = 'pot.dx'
   meshFile(2) = 'potDelta.dx'
   meshFile(3) = 'rhoV.dx'
   meshFile(4) = 'rhoVDelta.dx'
   meshFile(5) = 'rhoT.dx'
   colLabel(1) = '         pot'
   colLabel(2) = '    potDelta'
   colLabel(3) = '        rhoV'
   colLabel(4) = '   rhoVDelta'
   colLabel(5) = '        rhoT'

end subroutine initEnv


! This subroutine will read the command line parameters (CLPs) that were given.
!   The first CLP is the name of the file with the potential coefficients in
!   it.  This is usually "gs_scf-pot.dat", but it could be a different edge
!   from XANES calculations too such as "1s_scf-pot.dat".  The next three CLPs
!   are the number of points in each direction of the mesh for evaluating the
!   charge density on.  The first is the a axis, then b, and finally c.
subroutine readCommandLine

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define the locally used variables.
   character*25 :: commandBuffer
   integer :: i

   open(10,file='command',position='APPEND',status='unknown',form='formatted')
   write (10,ADVANCE='NO',fmt='(a17,a1)') 'makeFittedRho.exe',' '

   ! Store the first argument on the command line.  This should be the name of
   !   the potential file that contains the charge data.
   call getarg(1,commandBuffer)

   ! Read the variable from the command buffer into the file name for the
   !   charge data.  (It is probably something like gs_scfV-fb.dat.)
   read (commandBuffer,*) gaussFile
   write (10,ADVANCE='NO',fmt='(a15,a1)') gaussFile,' '


   ! Loop through the three axes to obtain the number of mesh points in each
   !   direction consecutively from the command line.
   do i = 1, 3

      ! Store the argument from the command line in the command buffer.
      call getarg(i+1,commandBuffer)

      ! Read the argument into the numMeshPoints array.
      read (commandBuffer,*) numMeshPoints(i)
      write (10,ADVANCE='NO',fmt='(i5,a1)') numMeshPoints(i),' '
   enddo

   ! Compute the total number of mesh points.
   numTotalMeshPoints = product(numMeshPoints(:))


   ! Store the next argument on the command line.  This argument specifies how
   !   many neighboring replicated cells should be included in the charge
   !   density plot of the primary cell.
   call getarg(5,commandBuffer)

   ! Read the variable from the command buffer.  It is important to note that
   !   the value of 0 means that we check neighbors from -0 to 0 (i.e. none).
   !   A value of 1 means we check from -1 to 1, that is the 26 neighboring
   !   cells and the primary cell.  A value of 2 means we check from -2 to 2,
   !   that is the 124 neighboring cells and the primary cell.  0 should be
   !   used when computing the charge density of a biomolecule that has
   !   significant vacuum around it in its cell.  1 should be used for most
   !   other regular solids.  2 should be used for the special case of a shard
   !   shaped primary cell because the second nearest neighbors may be close
   !   enough to affect the charge density in the primary cell.  Obviously,
   !   the higher the number, the longer the calculation.
   read (commandBuffer,*) neighborCell
   write (10,ADVANCE='NO',fmt='(i5,a1)') neighborCell,' '

   ! Store the next command line argument.  This argument specifies which type
   !   of plot to create.  If the value is a 0, then an openDX style data file
   !   is created for the 3D mesh.  If the value is non-zero, then there will
   !   be no openDX style data files created.
   call getarg(6,commandBuffer)
   read (commandBuffer,*) noOpenDX
   write (10,fmt='(i4)') noOpenDX

   close (10)

end subroutine readCommandLine




! The purpose of this subroutine is to open the two input files and make a
!   preliminary parse of them.  The goal is to obtain the limits on the
!   dimensions of the data structures that hold the pot site data and charge
!   data.
subroutine openFilesPrelimRead

   ! Import the necessary kinds definitions.
   use O_Kinds

   ! Make sure that no variables are accidentally defined.
   implicit none


   ! Define the locally used variables.
   integer :: i
   integer :: ioStatus ! I/O Status (-1 when EOF (End Of File)).
   integer :: currentAtom
   integer :: currentTypeNumber
   integer :: previousTypeNumber
   real (kind=double) :: currentAlpha
   real (kind=double) :: currentMaxAlpha
   real (kind=double), dimension(numCols) :: placeHolder
   character*8 :: charPlaceHolder1

   open(4,file='structure.dat',status='old',form='formatted')
   open(8,file=gaussFile,status='old',form='formatted')



   ! Read the structure file and determine the number of sites as well as the
   !   maximum number of potential sites for any given type.

   ! Read past the lattice parameters.
   do i = 1,5
      read (4,*)
   enddo

   ! Read the number of atomic/potential sites.
   read (4,*) numAtoms

   ! Read past the header for the list of atomic sites.
   read (4,*)

   ! Count the number of types in the system.

   ! Initialize the ioStatus to the noerror state.
   ioStatus = 0

   ! Initialize the count of the number of types.
   numTypes = 0

   ! Read until we reach the end of the file.
   do i = 1, numAtoms

      ! The first value is the atom number and the second is the type number.
      read (4,*,IOSTAT=ioStatus) placeHolder(1), currentTypeNumber

      ! Record the highest type number found as the number of types.
      if (currentTypeNumber > numTypes) then
         numTypes = currentTypeNumber
      endif
   enddo

   ! Allocate space to hold the number of potential sites for each potential
   !   type.  Also allocate space to hold the number of terms for each
   !   potential type.
   allocate (numSitesPerType(numTypes))
   allocate (numTermsPerType(numTypes))

   ! Prepare to make a second pass through the structure file to count other
   !   data.
   rewind (4)

   ! Read past the lattice parameters, headers, atomic sites, and more headers.
   do i = 1,10+numAtoms
      read (4,*)
   enddo

   ! Initialize the counts for the number of sites and terms for each type to
   !   zero.
   numSitesPerType(:) = 0
   numTermsPerType(:) = 0

   ! Initialize the previous type number.
   previousTypeNumber = 0  ! The first site is always type #1 (so prev=0).

   ! Loop through the remaining potential sites.
   do i = 1,numAtoms

      ! Read the potential site information.
      read (4,*) currentAtom,currentTypeNumber,placeHolder(1:3),&
            & charPlaceHolder1

      ! Determine if the current type is different than the last type.  It
      !   should be true that the atoms are in type sorted order in this file.
      if (currentTypeNumber /= previousTypeNumber) then

         ! Initialize the number of sites for this type.
         numSitesPerType(currentTypeNumber) = 1

         ! Make the current type number become the previous type number for the
         !   next iteration.
         previousTypeNumber = currentTypeNumber
      else

         ! Increment the number of sites for this type.
         numSitesPerType(currentTypeNumber) = &
               & numSitesPerType(currentTypeNumber) + 1
      endif
   enddo

   ! Compute the maximum number of sites of all types.
   maxNumSitesPerType = maxval(numSitesPerType(:))


   ! Read the potRho file and obtain the maximum number of terms for any type,
   !   and the number of terms for each type.  Also obtain the total number
   !   of atomic potential terms in the system.
   numAtomicPotTerms  = 0
   maxNumTermsPerType = 0
   currentTypeNumber  = 1
   currentMaxAlpha    = 0.0_double
   ioStatus           = 0
   do while (ioStatus == 0)

      ! Read the current line and if there is no error or EOF, then continue.
      read (8,*,IOSTAT=ioStatus) currentAlpha,placeHolder(1:numCols)

      ! Process the line if the read was successful.
      if (ioStatus == 0) then

         if (currentAlpha < currentMaxAlpha) then

            ! We have found a new potential type beginning here.
            currentTypeNumber = currentTypeNumber + 1

            ! Initialize the number of terms for this new type.
            numTermsPerType(currentTypeNumber) = 1

            ! Record that the current maximum alpha for this term is the
            !   current alpha.  This will be incremented with each read of a
            !   line until we find a new alpha that is less than it again.
            currentMaxAlpha = currentAlpha
         else
            ! Increase the value of the currentMaxAlpha so that if a new type
            !   starts on the next line it can be detected (because the first
            !   alpha for the new type will be less than the last alpha for the
            !   current type).
            currentMaxAlpha = currentAlpha

            ! Increment the number of terms for the current type.
            numTermsPerType(currentTypeNumber) = &
                  & numTermsPerType(currentTypeNumber) + 1
         endif

         ! Always increment the number of atomic potential terms.
         numAtomicPotTerms = numAtomicPotTerms + 1
      endif
   enddo

   ! Compute the maximum number of terms in any given type.
   maxNumTermsPerType = maxval(numTermsPerType(:))


   ! Allocate space to hold the information for the sites of the Gaussians, and
   !   the definitions of those Gaussians in terms of the charge.
   allocate (sitesData (3,maxNumSitesPerType,numTypes))
   allocate (gaussCoeffs (numCols,maxNumTermsPerType,numTypes))
   allocate (gaussAlphas (maxNumTermsPerType,numTypes))



   ! Rewind both files to be re-read for the actual data.
   rewind (4)
   rewind (8)

end subroutine openFilesPrelimRead



subroutine readData

   ! Import the necessary kinds values.
   use O_Kinds

   ! Make sure that no variables are accidentally defined.
   implicit none


   ! Define the locally used variables.
   integer :: i,j
   integer :: currentTypeNumber
   integer :: currentAtomNumber
   character*8 :: charPlaceHolder1

   ! Read the charge data first.
   do i = 1, numTypes
      do j = 1, numTermsPerType(i)
         read (8,*) gaussAlphas(j,i),gaussCoeffs(1:numCols,j,i)
      enddo
   enddo

   ! Read the lattice parameters next (in a.u.)
   read (4,*)
   read (4,*) latticeVectors(:,:)

   ! Skip the atomic site data.
   read (4,*)
   read (4,*)
   read (4,*)
   do i = 1, numAtoms
      read (4,*)
   enddo


   ! Read the potential site data.
   read (4,*)
   read (4,*)
   read (4,*)
   do i = 1, numTypes
      do j = 1, numSitesPerType(i)
         read (4,*) currentAtomNumber, currentTypeNumber, sitesData(:,j,i), &
               & charPlaceHolder1
      enddo
   enddo

   ! Determine the fractional stride length for each step in the mesh for each
   !   axis.
   do i = 1, 3
      fractStrideLength(i) = 1 / real(numMeshPoints(i),double)
   enddo

   ! Compute the properties of the lattice including the a,b,c magnitudes,
   !   the angles between the lattice vectors (alpha,beta,gamma:  bc,ac,ab),
   !   and the cross sectional areas of each vector pair bc, ac, ab.
   call latticeProperties

end subroutine readData



! Compute the properties of the lattice including the a,b,c magnitudes,
!   the angles between the lattice vectors (alpha,beta,gamma:  bc,ac,ab),
!   and the cross sectional areas of each vector pair bc, ac, ab.
subroutine latticeProperties

   ! Use necessary modules
   use O_Kinds
   use O_Constants

   ! Make sure nothing funny is declared.
   implicit none

   ! Declare local variables.
   integer :: i
   integer :: index1, index2

   ! Compute the magnitudes of a, b, and c.
   do i = 1, 3
      latticeMag(i) = sqrt(sum(latticeVectors(:,i)**2))
   enddo
!write (6,*) "latticeMag(:)=",latticeMag(:)


   ! Compute the unit vectors of the lattice vectors.
   do i = 1, 3
      latticeUnitVectors(:,i) = latticeVectors(:,i) / latticeMag(i)
!write (6,*) i,"latticeUnitVectors(:,i) = ",latticeUnitVectors(:,i)
   enddo


   ! Compute the angles between the vectors (alpha,beta,gamma:  bc,ac,ab).
   do i = 1,3
      index1 = mod(i,3)+1
      index2 = mod(i+1,3)+1
      latticeAngles(i) = acos(dot_product(latticeVectors(:,index1),&
            & latticeVectors(:,index2)) / latticeMag(index1) / &
            & latticeMag(index2))
   enddo
!write (6,*) "latticeAngles bc ac ab = ",latticeAngles(:)*180.0_double/pi

   ! Compute the fractional cross sectional areas (bc, ac, ab).  This is the
   !   cross sectional area of a unit in the mesh and essentially represents
   !   the amount of area that one data point represents.  NOTE:  The cross
   !   sectional area of ab is not influenced by the c vector's orientation.
   !   If a and b were perpendicular then the fractCrossArea(3) would simply be
   !   a*b.  However, they may be non-orthogonal with an angle between them
   !   (gamma) /= 90 degrees.  In this case the area is computed as the area
   !   for a parallelogram = a*b*sin(gamma).
   do i = 1,3
      index1 = mod(i,3)+1
      index2 = mod(i+1,3)+1
      fractCrossArea(i) = latticeMag(index1) * fractStrideLength(index1) * &
                        & latticeMag(index2) * fractStrideLength(index2) * &
                        & sin(latticeAngles(i))
!write (6,*) "sin(theta)=",sin(latticeAngles(i))
!write (6,*) latticeMag(mod(i,3)+1),latticeMag(mod(i+1,3)+1)
   enddo
!write (6,*) "fractCrossArea(:)=",fractCrossArea(:)


   ! Compute the unit normal vectors that define the orientation of each
   !   plane (bc, ac, ab).  (This is used to help find the volume of space
   !   represented by one data point of the mesh.)
   do i = 1,3
      index1 = mod(i,3)+1
      index2 = mod(i+1,3)+1
      call crossProduct(unitNormal(:,i),latticeUnitVectors(:,index1),&
            latticeUnitVectors(:,index2))
      unitNormal(:,i) = unitNormal(:,i) / sin(latticeAngles(i))
!write (6,*) i,"unitNormal(:,i) = ",unitNormal(:,i)
   enddo


   ! Compute the angles between the unit normal vectors of the planes and the
   !   lattice unit vector not in the plane.  Angles are (bc,a ; ac,b ; ab,c)
   !   where the first is the plane (normal) and the second is the other unit
   !   lattice vector.  This angle is necessary because (for example) if the
   !   c lattice vector is *not* perpendicular to the ab plane then the volume
   !   represented by a mesh data point will *not* be area(ab)*mag(c).  Instead
   !   it will be area(ab)*mag(c)*sin(this angle).
   do i = 1,3
      index1 = mod(i,3)+1
      index2 = mod(i+1,3)+1
      planeAngles(i) = &
            asin(dot_product(unitNormal(:,i),latticeUnitVectors(:,i)))
   enddo
!write (6,*) "planeAngles bc,a ac,b ab,c = ",planeAngles(:)*180.0_double/pi

end subroutine latticeProperties


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


subroutine printODXHead

   ! Import the necessary kinds definitions.
   use O_Kinds

   ! Import program constants
   use O_Constants

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define local variables.
   integer :: i

   open (unit=17,file='lattice.dx',form='formatted')
   do i = 1, numCols
      open (unit=17+i,file=meshFile(i),form='formatted')
   enddo

   do i = 1,numCols
      write (17+i,*) "object 1 class gridpositions counts ",numMeshPoints(:)
      write (17+i,*) "origin 0 0 0"
      write (17+i,*) "delta ",latticeVectors(:,1)*fractStrideLength(1)*bohrRad
      write (17+i,*) "delta ",latticeVectors(:,2)*fractStrideLength(2)*bohrRad
      write (17+i,*) "delta ",latticeVectors(:,3)*fractStrideLength(3)*bohrRad
      write (17+i,*)
      write (17+i,*) "object 2 class gridconnections counts ",numMeshPoints(:)
      write (17+i,*)
      write (17+i,*) "object 3 class array type float rank 0 items ",&
            & numTotalMeshPoints," data follows"
   enddo

end subroutine printODXHead



subroutine computeMeshValues

   ! Import the necessary kinds definitions.
   use O_Kinds

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define the locally used variables
   integer :: i
   integer :: a,b,c ! Counters for the number of points along each axis.
   integer :: aLatticeShift,bLatticeShift,cLatticeShift
   integer :: currentType
   integer :: currentTerm
   integer :: currentSite
   integer :: totalPoints ! Updates user with progress.
   integer :: pointLineCounter ! Counts the number of points to be printed/line.
   integer :: pointsPerSTDOUTDot
   real (kind=double) :: currentSepSqrd
   real (kind=double) :: currentExp
   real (kind=double), dimension(numCols) :: currentPointValue


   ! Allocate space to hold the profiles and initialize the values to zero for
   !   accumulation.
   allocate(profile(3,numCols,maxVal(numMeshPoints(:))))
   profile (:,:,:) = 0.0_double

   ! Initialize a counter for the total number of points.
   totalPoints = 0

   ! Initialize the counter for the number of points to be printed per line.
   !   Every 5 points a new line will be created.
   pointLineCounter = 0

   ! Initialize the current separation to avoid compiler warnings.
   currentSepSqrd = 0

   ! Initialize the number of points that are represented by one dot on the
   !   standard output for the user to watch the progress.
   pointsPerSTDOUTDot = 1000

   ! Loop over the three axial directions for the number of points requested
   !   on the command line.  These are not indented to aid in reading the code.
   !   Also note that the indeces go from 0 to n-1 because we multiply the
   !   index by the stride length and we want to start at 0,0,0 and we want to
   !   end at a,b,c.
   do a = 0, numMeshPoints(1)-1
   do b = 0, numMeshPoints(2)-1
   do c = 0, numMeshPoints(3)-1

      ! Increment the total number of points considered.
      totalPoints = totalPoints + 1

      ! Print a record of the progress.
      if (mod(totalPoints,pointsPerSTDOUTDot*50) == 0) then
         write (6,FMT="(a2,i12)") "| ",totalPoints
      elseif (mod(totalPoints,pointsPerSTDOUTDot*10) == 0) then
         write (6,ADVANCE="NO",FMT="(a1)") "|"
      elseif (mod(totalPoints,pointsPerSTDOUTDot) == 0) then
         write (6,ADVANCE="NO",FMT="(a1)") "."
      endif
      call flush (6)

      ! Initialize the value of the current mesh point.
      currentPointValue(:) = 0.0_double

      ! Initiate a set of loops through the requested number of neighboring
      !   cells that border the primary cell.  Again, I will only indent the
      !   group as a whole and not the individual loops.
      do aLatticeShift = -neighborCell, neighborCell
      do bLatticeShift = -neighborCell, neighborCell
      do cLatticeShift = -neighborCell, neighborCell

         ! Perform another triple loop where we consider each term for each
         !   site of each type.  Between the 2nd and 3rd loops we compute the
         !   separation distance between the current mesh point and the current
         !   site.
         do currentType = 1, numTypes
         do currentSite = 1, numSitesPerType(currentType)
            call getSeparation (a,b,c,aLatticeShift,bLatticeShift,&
                  & cLatticeShift,sitesData(:,currentSite,currentType),&
                  & currentSepSqrd)
         do currentTerm = 1, numTermsPerType(currentType)

            ! Evaluate the Gaussian defined by the alpha of the current term
            !   with the coefficient given by the total charge density.
            currentExp = gaussAlphas(currentTerm,currentType) * currentSepSqrd

            ! Update only if the value is not negligable.  If it is negligable,
            !   then we can abort the loop since no others for this site will
            !   count either.
            if (currentExp < 15) then
               currentPointValue(:) = currentPointValue(:) + &
                     & gaussCoeffs(:,currentTerm,currentType) * &
                     & exp(-currentExp)
            else
               exit
            endif
         enddo
         enddo
         enddo

      enddo
      enddo
      enddo

      if (noOpenDX == 0) then

         ! Print the accumulated result for the currentPointValue
         do i = 1,numCols
            write (17+i,ADVANCE="NO",fmt="(1x,f12.8)") currentPointValue(i)
         enddo

         ! Increment the counter for the number of values printed for this line.
         pointLineCounter = pointLineCounter + 1

         ! If the number of values printed this line == 5, then print a newline
         !   and reset the pointLineCounter.
         if (pointLineCounter == 5) then
            do i = 1,numCols
               write (17+i,*)
            enddo
            pointLineCounter = 0
         endif
      endif

      ! Accumulate the data for the profiles.
      do i = 1,numCols  ! Loop over all columns of data
         profile(1,i,a+1) = profile(1,i,a+1) + currentPointValue(i)
         profile(2,i,b+1) = profile(2,i,b+1) + currentPointValue(i)
         profile(3,i,c+1) = profile(3,i,c+1) + currentPointValue(i)
      enddo

   enddo
   enddo
   enddo

   ! Finalize the record of points written to stdout.
   if (mod(totalPoints,pointsPerSTDOUTDot*50) /= 0) then
      write (6,*)
   endif

   ! Finalize printing the openDX data file.
   if (noOpenDX == 0) then
      ! Write a newline to finish out any uneven (incomplete) lines if needed.
      if (pointLineCounter /= 0) then
         do i = 1,numCols
            write (17+i,*)
         enddo
      endif
   endif

   ! Print out the profile data.
   call printProfile

end subroutine computeMeshValues



subroutine printProfile

   ! Import necessary variables.
   use O_Kinds

   ! Import program constants
   use O_Constants

   ! Make sure nothing funny is accidentally declared.
   implicit none

   ! Declare local variables.
   integer :: i,j
   integer :: index1,index2
   real (kind=double), dimension(numCols) :: integral

   ! Open the profile data files.
   open (unit=30,file='profile_a.dat',form='formatted')
   open (unit=31,file='profile_b.dat',form='formatted')
   open (unit=32,file='profile_c.dat',form='formatted')

   ! Print the header.
   write (30,fmt="(6a12)") "aPos",colLabel(:)
   write (31,fmt="(6a12)") "bPos",colLabel(:)
   write (32,fmt="(6a12)") "cPos",colLabel(:)

   ! Adjust the charge profiles to average out the cross sectional area effect
   !   and the fact that the profiles are simple accumulations of the other
   !   axes.  For the case of the potential we simply obtain the average over
   !   the plane.
   do i = 1, 3
      index1 = mod(i,3)+1
      index2 = mod(i+1,3)+1
      profile(i,1:2,:) = profile(i,1:2,:) * hartree / &
            & numMeshPoints(index1) / numMeshPoints(index2)
      profile(i,3:numCols,:numMeshPoints(i)) = &
            & profile(i,3:numCols,:numMeshPoints(i)) * fractCrossArea(i)
   enddo

   ! Print the profiles.
   do i = 1,3
      do j = 1, numMeshPoints(i)
         write (29+i,fmt="(6e12.4)") (j-1) * fractStrideLength(i) * &
               & latticeMag(i) * bohrRad, profile(i,1:2,j), &
               & profile(i,3:numCols,j) / bohrRad
      enddo
   enddo

   ! Close the profile data files.
   close (30)
   close (31)
   close (32)

   ! Integrate the charge and print the result to standard output.
   do i = 1, 3
      integral(:) = 0.0_double
      do j = 1, numMeshPoints(i)
         integral(1:2) = integral(1:2) + profile(i,1:2,j) / &
               numMeshPoints(i)
         integral(3:numCols) = integral(3:numCols) + profile(i,3:numCols,j) * &
               & fractStrideLength(i) * latticeMag(i) * sin(planeAngles(i))
      enddo
      write (6,fmt="(6a12)") colLabel(:),"Integrated"
      write (6,fmt="(5e12.4)") integral(:)
   enddo

end subroutine printProfile



subroutine getSeparation (a,b,c,aLatticeShift,bLatticeShift,cLatticeShift,&
      & currentSitePos,currentSepSqrd)

   ! Import the necessary kinds definitions.
   use O_Kinds

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define the passed variables.
   integer :: a,b,c ! Counters for the number of points along each axis.
   integer :: aLatticeShift,bLatticeShift,cLatticeShift
   real (kind=double), dimension (3) :: currentSitePos
   real (kind=double) :: currentSepSqrd

   ! Define local variables.
   integer :: i
   real (kind=double), dimension (3) :: currentPoint


   ! Compute the location in x,y,z coordinates of the current mesh point by
   !   using the fractional stride length, the lattice vectors, and the lattice
   !   shift.
   do i = 1, 3
      currentPoint(i) = (a*fractStrideLength(1)+aLatticeShift) * &
                      &  latticeVectors(i,1) + &
                      & (b*fractStrideLength(2)+bLatticeShift) * &
                      &  latticeVectors(i,2) + &
                      & (c*fractStrideLength(3)+cLatticeShift) * &
                      &  latticeVectors(i,3)
   enddo

   ! Compute the separation between the current mesh point and the current
   !   site point.
   currentSepSqrd = sum((currentPoint(:)-currentSitePos(:))**2)

end subroutine getSeparation



subroutine printODXTail

   ! Import the necessary kinds definitions.
   use O_Kinds

   ! Import program constants
   use O_Constants

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define local variables.
   integer :: i

   ! Write the information for the lattice.
   write (17,*) "object 1 class gridpositions counts 2 2 2"
   write (17,*) "origin 0 0 0"
   write (17,*) "delta ",latticeVectors(:,1) * bohrRad
   write (17,*) "delta ",latticeVectors(:,2) * bohrRad
   write (17,*) "delta ",latticeVectors(:,3) * bohrRad
   write (17,*)
   write (17,*) "object 2 class gridconnections counts 2 2 2"
   write (17,*)
   write (17,*) "object 3 class array type float rank 0 items 8 data follows"
   write (17,*) "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0"
   write (17,*) 'attribute "dep" string "positions"'
   write (17,*)
   write (17,*) 'object "lattice" class field'
   write (17,*) 'component "positions" value 1'
   write (17,*) 'component "connections" value 2'
   write (17,*) 'component "data" value 3'

   ! End the openDX lattice file.
   write (17,*) 'end'
   close (17)

   do i = 1,numCols
      ! Write the tail for the pot, rhoT, and rhoV and the define the dataField
      !   object that ties all the components for the data sets together.
      write (17+i,*) 'attribute "dep" string "positions"'
      write (17+i,*)
      write (17+i,*) 'object "dataField" class field'
      write (17+i,*) 'component "positions" value 1'
      write (17+i,*) 'component "connections" value 2'
      write (17+i,*) 'component "data" value 3'

      ! End the openDX field file.
      write (17+i,*) 'end'

      close (17+i)
   enddo

end subroutine printODXTail

end module ChargeDensityDataModule


program makeRho

   ! Import the precision variables
   use O_Kinds

   ! Import program constants
   use O_Constants

   ! Import necessary data modules
   use ChargeDensityDataModule

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Initialize environment.
   call initEnv

   ! Parse the command line to determine the mesh size
   call readCommandLine

   ! Open the data files and make a first pass through them to collect the
   !   necessary information to assign the dimensions of the matrices.
   call openFilesPrelimRead


   ! Re-read both files (structure.dat and output xx_scf-pot.dat) to collect
   !   all the necessary data.
   call readData


   ! Print out the beginning information for the openDX file.  This includes
   !   the mesh parameters, the number of points in each direction, the atomic
   !   positions and size descriptions.
   if (noOpenDX == 0) then
      call printODXHead
   endif

   ! Loop through the mesh points and evalute the charge at each point from
   !   each term of each type.
   call computeMeshValues


   ! Print out the last little tail of the openDX file.  This includes the
   !   definition of the field that contains the three objects for plotting
   !   the charge density.  This will also write the extra file that contains
   !   the information for the lattice frame.
   if (noOpenDX == 0) then
      call printODXTail
   endif


end program makeRho
