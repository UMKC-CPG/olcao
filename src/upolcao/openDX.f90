module O_OpenDX

   contains

subroutine printODXFieldHead(numFields)

   ! Import the necessary data modules.
   use O_Constants, only: bohrRad
   use O_Lattice,   only: realVectors, realFractStrideLength, numMeshPoints

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define passed parameters
   integer :: numFields

   ! Define local variables.
   integer :: i ! Loop index.
   character*7 :: fileName

   do i = 1, numFields

      ! Create the fortran file name.
      write (fileName,fmt='(a,i2)') "fort.",57+i

      ! Open the field file.
      open (unit=57+i,file=fileName,status='new',form='formatted')

      ! Write the header for this file describing the mesh step size in each
      !   of the a,b,c axes (given as x,y,z vectors), the number of mesh points
      !   in each a,b,c direction, and the total number of mesh points.
      write (57+i,*) "object 1 class gridpositions counts ",numMeshPoints(:)
      write (57+i,*) "origin 0 0 0"
      write (57+i,*) "delta ",realVectors(:,1)*realFractStrideLength(1)*bohrRad
      write (57+i,*) "delta ",realVectors(:,2)*realFractStrideLength(2)*bohrRad
      write (57+i,*) "delta ",realVectors(:,3)*realFractStrideLength(3)*bohrRad
      write (57+i,*)
      write (57+i,*) "object 2 class gridconnections counts ",numMeshPoints(:)
      write (57+i,*)
      write (57+i,*) "object 3 class array type float rank 0 items ",&
            & product(numMeshPoints(:))," data follows"
      call flush (57+i)
   enddo

end subroutine printODXFieldHead



subroutine printODXFieldTail(numFields)


   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define passed parameters
   integer :: numFields

   ! Define local variables.
   integer :: i ! Loop index.

   do i = 1,numFields 

      ! Write the tail for the charge density and the define the dataField
      !   object that ties all the components for the charge density together.
      write (57+i,*) 'attribute "dep" string "positions"'
      write (57+i,*)
      write (57+i,*) 'object "dataField" class field'
      write (57+i,*) 'component "positions" value 1'
      write (57+i,*) 'component "connections" value 2'
      write (57+i,*) 'component "data" value 3'

      ! End the openDX field file.
      write (57+i,*) 'end'

      ! Close the openDX field file.
      close (57+i)
   enddo

end subroutine printODXFieldTail


subroutine printODXLattice

   ! Import the necessary modules.
   use O_Constants, only: bohrRad
   use O_Lattice,   only: realVectors

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Open the lattice definition file for OpenDX.
   open (unit=57,file='fort.57',status='new',form='formatted') ! ODX Lattice

   ! Write the information for the lattice.
   write (57,*) "object 1 class gridpositions counts 2 2 2"
   write (57,*) "origin 0 0 0"
   write (57,*) "delta ",realVectors(:,1)*bohrRad
   write (57,*) "delta ",realVectors(:,2)*bohrRad
   write (57,*) "delta ",realVectors(:,3)*bohrRad
   write (57,*)
   write (57,*) "object 2 class gridconnections counts 2 2 2"
   write (57,*)
   write (57,*) "object 3 class array type float rank 0 items 8 data follows"
   write (57,*) "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0"
   write (57,*) 'attribute "dep" string "positions"'
   write (57,*)
   write (57,*) 'object "lattice" class field'
   write (57,*) 'component "positions" value 1'
   write (57,*) 'component "connections" value 2'
   write (57,*) 'component "data" value 3'

   ! End the openDX lattice file.
   write (57,*) 'end'

   ! Close the openDX lattice file.
   close (57)

end subroutine printODXLattice

subroutine printODXAtomPos

   ! Import the necessary modules.
   use O_Kinds
   use O_Constants,   only: bohrRad
   use O_PotTypes,    only: potTypes
   use O_ElementData, only: colorDX, covalRadii
   use O_AtomicSites, only: numAtomSites, atomSites

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define local variables.
   integer :: i
   integer :: currType ! Current atomic type in a loop.

   ! Open the atom position definition file for OpenDX.
   open (unit=56,file='fort.56',status='new',form='formatted') ! ODX Atoms

   ! Write the information for the atomic positions.
   write (56,fmt='(a,i10,a)') &
         & "object 1 class array type float rank 1 shape 3 items ", &
         & numAtomSites," data follows"
   do i = 1, numAtomSites
      write (56,fmt='(3e16.8)') atomSites(i)%cartPos(:)*bohrRad
   enddo
   write (56,*)


   ! Print the color values for each atom.
   write (56,fmt='(a,i10,a)') &
         & "object 2 class array type float rank 0 items ", &
         & numAtomSites," data follows"
   do i = 1, numAtomSites
      currType = atomSites(i)%atomTypeAssn
      write (56,fmt='(e16.8)') colorDX(int(potTypes(currType)%nucCharge))
   enddo
   write (56,fmt='(a)') 'attribute "dep" string "positions"'
   write (56,*)


   ! Add the atom sphere size information.
   write (56,fmt='(a,i10,a)') &
         & "object 3 class array type float rank 0 items ", &
         & numAtomSites, " data follows"
   do i = 1, numAtomSites
      currType = atomSites(i)%atomTypeAssn
      write (56,fmt='(e16.8)') covalRadii(int(potTypes(currType)%nucCharge))
   enddo
   write (56,fmt='(a)') 'attribute "dep" string "positions"'
   write (56,*)

   ! Combine the above three objects into the spheres field.
   write (56,fmt='(a)') 'object "atoms" class field'
   write (56,fmt='(a)') 'component "positions" value 1'
   write (56,fmt='(a)') 'component "data" value 2'
   write (56,fmt='(a)') 'component "sizes" value 3'
   write (56,fmt='(a)') 'end'
   write (56,*)

   ! End the openDX atom position file.
   write (56,fmt='(a)') 'end'

   ! Close the openDX atom position file.
   close (56)

end subroutine printODXAtomPos


end module O_OpenDX
