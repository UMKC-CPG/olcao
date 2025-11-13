module O_XDMF_VTK

   contains

subroutine printXDMFMetaFile

   ! Import the necessary data modules.
   use O_Constants, only: bohrRad
   use O_Lattice,   only: realVectors, realFractStrideLength, numMeshPoints
   use O_FieldHDF5, only: dataSetNames, groupNames, fileName, numDataSets

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define local variables.
   integer :: i,j
   integer :: dataSetIndex

   ! Open the field file.
   open (unit=78,file="fort.78",status='new',form='formatted')

   ! Write the header for this file opening the XML tags.
   write (78,fmt="(a)") "<?xml version=""1.0"" encoding=""utf-8""?>"
   write (78,fmt="(a)") "<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
   write (78,fmt="(a)") "<Xdmf Version=""3.0"" " // &
         & "xmlns:xi=""[http://www.w3.org/2001/XInclude]"">"
   write (78,fmt="(a)") "  <Domain>"
   write (78,fmt="(a)") "    <Grid Name=""Grid"" GridType=""Uniform"">"
   write (78,fmt="(a)") "      <Geometry Name=""Geometry"" " // &
         & "GeometryType=""XYZ"">"
   write (78,fmt="(a,3i6,a,a,a)") "        <DataItem NumberType=" // &
         & """Float"" Dimensions=""",numMeshPoints(:)," 3""" // &
         & " Format=""HDF"">",trim(fileName),":/meshGroup/mesh</DataItem>"
   write (78,fmt="(a)") "      </Geometry>"
   write (78,fmt="(a,3i10,a)") "      <Topology Dimensions=""", &
         & numMeshPoints(:),""" TopologyType=""3DSMesh""/>"
   dataSetIndex = 0
   do i = 1, size(groupNames)
      do j = 1, 4
         dataSetIndex = dataSetIndex + 1
         write (78,fmt="(a,a,a)") "      <Attribute Center=""Node"" Name=""", &
               & dataSetNames(dataSetIndex),""" AttributeType=""Scalar"">"
         write (78,ADVANCE="NO",fmt="(a,3i10)") "        <DataItem " // &
               & "NumberType=""Float"" Dimensions=""",numMeshPoints(:)
         write (78,fmt="(7a)") """ Format=""HDF"">",trim(fileName),":/", &
               & trim(groupNames(i)),"/",dataSetNames(dataSetIndex),"</DataItem>"
         write (78,fmt="(a)") "      </Attribute>"
      enddo
   enddo
   write (78,fmt="(a)") "    </Grid>"
   write (78,fmt="(a)") "  </Domain>"
   write (78,fmt="(a)") "</Xdmf>"

   ! Close the XDMF file.
   close (78)

end subroutine printXDMFMetaFile


subroutine printVTKHead

   ! Import the necessary data modules.
   use O_Constants, only: bohrRad
   use O_Lattice,   only: realVectors, realFractStrideLength, numMeshPoints

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define local variables.
   integer :: i,j

   ! Open the field file.
   open (unit=56,file="fort.56",status='new',form='formatted')

   ! Print the VTK Polydata header.
   write (56,fmt="(a)") "# vtk DataFile Version 4.2"
   write (56,fmt="(a)") "This is a model produced by the KC-O OLCAO program."
   write (56,fmt="(a)") "ASCII"
   write (56,fmt="(a)") "DATASET POLYDATA"
   write (56,*)


end subroutine printVTKHead


subroutine printVTKLattice

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

end subroutine printVTKLattice

subroutine printVTKAtomPos

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

end subroutine printVTKAtomPos


end module O_XDMF_VTK
