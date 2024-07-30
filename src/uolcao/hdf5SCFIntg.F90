module O_SCFIntegralsHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define the main subgroup from scf_fid that holds the multicenter
   !   Gaussian integral results.
   integer(hid_t) :: atomIntgGroup_gid ! Atom and potential overlap group.

   ! Begin group, dataspace, dataset definitions for atomIntgGroup_gid.

   ! Define the group IDs under atomIntgGroup_gid. Each group here will hold
   !   either a dataset for each kpoint, or a group for each kpoint.
   integer(hid_t) :: atomOverlap_gid
   integer(hid_t) :: atomKEOverlap_gid
   integer(hid_t) :: atomMVOverlap_gid
   integer(hid_t) :: atomNucOverlap_gid
   integer(hid_t) :: atomDMOverlap_gid
   integer(hid_t) :: atomMMOverlap_gid
   integer(hid_t) :: atomPotOverlap_gid

   ! Define the group IDs of the dynamically numbered subgroups. I.e., for
   !   the groups that hold another group, we define that collection of
   !   subgroups here.
   integer(hid_t), allocatable, dimension (:) :: atomDMxyzOL_gid
   integer(hid_t), allocatable, dimension (:) :: atomMMxyzOL_gid
   integer(hid_t), allocatable, dimension (:) :: atomPotTermOL_gid

   ! The dataspaces of each dataset in atomIntgGroup_gid are the same in
   !   all characteristics (type, dimension, etc.) and therefore can be
   !   given a static ID definition now.
   integer(hid_t) :: valeVale_dsid

   ! Each of the below datasets for atomIntgGroup will use the same property
   !   list.
   integer(hid_t) :: valeVale_plid

   ! Define array that holds the dimensions of the dataset.
   integer(hsize_t), dimension (2) :: atomDims

   ! Define array that holds the dimensions of the chunk.
   integer (hsize_t), dimension (2) :: atomDimsChunk

   ! The number of datasets under each of the atomIntgGroup groups will
   !   vary depending on the problem.  (# of kpoints, potential dimension,
   !   or xyz terms). Thus, the datasets must be given IDs dynamically.
   integer(hid_t), allocatable, dimension (:)   :: atomOverlap_did
   integer(hid_t), allocatable, dimension (:)   :: atomKEOverlap_did
   integer(hid_t), allocatable, dimension (:)   :: atomMVOverlap_did
   integer(hid_t), allocatable, dimension (:)   :: atomNucOverlap_did
   integer(hid_t), allocatable, dimension (:,:) :: atomDMOverlap_did
   integer(hid_t), allocatable, dimension (:,:) :: atomMMOverlap_did
   integer(hid_t), allocatable, dimension (:,:) :: atomPotOverlap_did

   ! Define the attribute IDs that will be used for each group. One attribute
   !   is sufficient for each integral type except the potential overlap.
   !   For that one, we will need a set of attributes, one for each potential
   !   term.
   integer(hid_t) :: atomOverlap_aid
   integer(hid_t) :: atomKEOverlap_aid
   integer(hid_t) :: atomMVOverlap_aid
   integer(hid_t) :: atomNucOverlap_aid
   integer(hid_t) :: atomDMOverlap_aid
   integer(hid_t) :: atomMMOverlap_aid
   integer(hid_t), allocatable, dimension (:) :: atomPotTermOL_aid


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initSCFIntegralHDF5 (scf_fid, attribInt_dsid, attribIntDims)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints     ! For numKPoints
   use O_Potential   ! For potDim
   use O_AtomicSites ! For valeDim

   ! Define the passed parameters.
   integer(hid_t) :: scf_fid
   integer(hid_t) :: attribInt_dsid
   integer(hsize_t), dimension (1) :: attribIntDims

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
#ifndef GAMMA
   atomDims(1) = 2 ! Real and imaginary are needed.
#else
   atomDims(1) = 1 ! Real needed only
#endif
   atomDims(2) = valeDim*(valeDim+1)/2 ! Linear storage of 1/2 matrix.

   ! Check that the chunk size is not too large. The assumption here is that
   !   the number being stored are 8 byte reals and that we should not go over
   !   2 billion bytes. Note that if x*y = >250M and we want a*b = 250M then
   !   the additional requirement x/y = a/b leads to b = sqrt(250M/>250M)*y.
   !   Thus a = 250M/b.
   if (atomDims(1) * atomDims(2) > 250000000) then
      atomDimsChunk(2) = int(250000000/atomDims(1))
      atomDimsChunk(1) = atomDims(1)
   else
      atomDimsChunk(1) = atomDims(1)
      atomDimsChunk(2) = atomDims(2)
   endif

   ! Create the Integral group within the scf HDF5 file.
   call h5gcreate_f (scf_fid,"/atomIntgGroup",atomIntgGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom intg group'

   ! The basic layout of the contents of the integral group are as follows:

   ! For integrals that have only one matrix per kpoint (overlap, kinetic
   !   energy, nuclear potential, etc.) there will be one group for that
   !   integral type that will include one data set for every kpoint.

   ! For integrals that have multiple matrices per kpoint (electron potential,
   !   momentum matrix, dipole moment, etc.) we will create one group per
   !   kpoint and each kpoint group will contain all the relevant datasets
   !   for that kpoint.

   ! Create the subgroups within the atomIntgGroup. Even if a group (type of
   !   integral calculation) will not be done with the current call to olcao
   !   scf, we create all groups anyway so that if a future olcao scf call
   !   does request that group, it will already be available to access.
   call h5gcreate_f (atomIntgGroup_gid,"atomOverlap",atomOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom overlap group'
   call h5gcreate_f (atomIntgGroup_gid,"atomKEOverlap",atomKEOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create kinetic energy overlap group'
   call h5gcreate_f (atomIntgGroup_gid,"atomMVOverlap",atomMVOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create mass velocity overlap group'
   call h5gcreate_f (atomIntgGroup_gid,"atomNucOverlap",atomNucOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create atom nuclear overlap group'
   call h5gcreate_f (atomIntgGroup_gid,"atomDMOverlap",atomDMOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create dipole moment overlap group'
   call h5gcreate_f (atomIntgGroup_gid,"atomMMOverlap",atomMMOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create momentum matrix overlap group'
   call h5gcreate_f (atomIntgGroup_gid,"atomPotOverlap",atomPotOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create atom potential overlap group'

   ! For the integrals that have multiple different types of matrices, we
   !   need to create space to hold the group IDs of each matrix.
   allocate (atomDMxyzOL_gid (3)) ! Needs x,y,z matrices all done together
   allocate (atomMMxyzOL_gid (3)) ! Needs x,y,z matrices all done together
   allocate (atomPotTermOL_gid (potDim)) ! One matrix per pot term, separate.

   ! Create the subgroups for the DM and MM terms and assign a gid name
   !   following 1=x, 2=y, and 3=z.
   do i = 1, 3
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5gcreate_f (atomDMOverlap_gid,currentName,&
            & atomDMxyzOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to create dipole moment xyz group'

      call h5gcreate_f (atomMMOverlap_gid,currentName,&
            & atomMMxyzOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to create momentum matrix xyz group'
   enddo

   ! Create the subgroups for the potential terms and assign a gid name equal
   !   to the term number.
   do i = 1, potDim
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5gcreate_f (atomPotOverlap_gid,currentName,atomPotTermOL_gid(i),&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to create potential term group'
   enddo

   ! Sufficient space must be allocated to hold the dataset IDs for all
   !   matrices that need to be computed.
   allocate (atomOverlap_did     (numKPoints))
   allocate (atomKEOverlap_did   (numKPoints))
   allocate (atomMVOverlap_did   (numKPoints))
   allocate (atomNucOverlap_did  (numKPoints))
   allocate (atomDMOverlap_did   (numKPoints,3))
   allocate (atomMMOverlap_did   (numKPoints,3))
   allocate (atomPotOverlap_did  (numKPoints,potDim))

   ! Allocate space to hold the attribute IDs for the potential matrics.
   allocate (atomPotTermOL_aid  (potDim))

   ! Create the dataspace that will be used for each dataset in atomIntgGroup
   !   and all of its subgroups.  The same dataspace definition works for all
   !   of the datasets.
   call h5screate_simple_f(2,atomDims,valeVale_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create vale vale dsid'

   ! Define the properties of the datasets to be made.

   ! Create the property list first.  Then set the properties one at a time.
   call h5pcreate_f      (H5P_DATASET_CREATE_F,valeVale_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create vale vale plid'
   call h5pset_layout_f  (valeVale_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale plid layout'
   call h5pset_chunk_f   (valeVale_plid,2,atomDimsChunk,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale plid chunk size'
!   call h5pset_shuffle_f (valeVale_plid,hdferr)
   call h5pset_deflate_f (valeVale_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale for deflation'

   ! Define the dataspace of the computed "status" attribute that will be used
   !   for each group.
   call h5screate_simple_f (1,attribIntDims(1),attribInt_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create the attribInt_dsid'

   ! Create the datasets that will be used for all subgroups of atomIntgGroup.
   do i = 1, numKPoints
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5dcreate_f (atomOverlap_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,atomOverlap_did(i),hdferr,valeVale_plid)
      if (hdferr /= 0) stop 'Failed to create atom overlap did'

      call h5dcreate_f (atomKEOverlap_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,atomKEOverlap_did(i),hdferr,valeVale_plid)
      if (hdferr /= 0) stop 'Failed to create KE overlap did'

      call h5dcreate_f (atomMVOverlap_gid,currentName,H5T_NATIVE_DOUBLE,&
         & valeVale_dsid,atomMVOverlap_did(i),hdferr,valeVale_plid)
      if (hdferr /= 0) stop 'Failed to create MV overlap did'

      call h5dcreate_f (atomNucOverlap_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,atomNucOverlap_did(i),hdferr,valeVale_plid)
      if (hdferr /= 0) stop 'Failed to create nuclear overlap did'

      do j = 1, 3
         write (currentName,fmt="(i7.7)") j
         currentName = trim (currentName)

         call h5dcreate_f (atomDMxyzOL_gid(j),currentName,&
               & H5T_NATIVE_DOUBLE,valeVale_dsid,atomDMOverlap_did(i,j),&
               & hdferr,valeVale_plid)
         if (hdferr /= 0) stop 'Failed to create DM overlap did'

         call h5dcreate_f (atomMMxyzOL_gid(j),currentName,&
               & H5T_NATIVE_DOUBLE,valeVale_dsid,atomMMOverlap_did(i,j),&
               & hdferr,valeVale_plid)
         if (hdferr /= 0) stop 'Failed to create MM overlap did'
      enddo

      do j = 1, potDim
         write (currentName,fmt="(i7.7)") j
         currentName = trim (currentName)

         call h5dcreate_f(atomPotTermOL_gid(j),currentName,&
               & H5T_NATIVE_DOUBLE,valeVale_dsid,atomPotOverlap_did(i,j),&
               & hdferr,valeVale_plid)
         if (hdferr /= 0) stop 'Failed to create potential overlap did'
      enddo
   enddo

   ! Create the attribute that indicates completion of the calculation of
   !   each integral type.
   call h5acreate_f (atomOverlap_gid,"status",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,atomOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom overlap aid'

   call h5acreate_f (atomKEOverlap_gid,"status",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,atomKEOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create KE overlap aid'

   call h5acreate_f (atomMVOverlap_gid,"status",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,atomMVOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create MV overlap aid'

   call h5acreate_f (atomNucOverlap_gid,"status",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,atomNucOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create Nuc overlap aid'

   call h5acreate_f (atomDMOverlap_gid,"status",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,atomDMOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create DM overlap aid'

   call h5acreate_f (atomMMOverlap_gid,"status",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,atomMMOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create MM overlap aid'

   do i = 1, potDim
      call h5acreate_f (atomPotTermOL_gid(i),"status",H5T_NATIVE_INTEGER,&
            & attribInt_dsid,atomPotTermOL_aid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to create potential overlap aid (i)'
   enddo

   ! Initialize all dataset attributes to the uncomputed (status = zero) state.
   call h5awrite_f(atomOverlap_aid,H5T_NATIVE_INTEGER,0,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomOverlap_aid'

   call h5awrite_f(atomKEOverlap_aid,H5T_NATIVE_INTEGER,0,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomKEOverlap_aid'

   call h5awrite_f(atomMVOverlap_aid,H5T_NATIVE_INTEGER,0,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomMVOverlap_aid'

   call h5awrite_f(atomNucOverlap_aid,H5T_NATIVE_INTEGER,0,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomNucOverlap_aid'

   call h5awrite_f(atomDMOverlap_aid,H5T_NATIVE_INTEGER,0,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomDMOverlap_aid'

   call h5awrite_f(atomMMOverlap_aid,H5T_NATIVE_INTEGER,0,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomMMOverlap_aid'

   do i = 1, potDim
      call h5awrite_f(atomPotTermOL_aid(i),H5T_NATIVE_INTEGER,0,&
            & attribIntDims,hdferr)
      if (hdferr /= 0) stop 'Failed to initialize atomPotTermOL_aid(i)'
   enddo

   ! At this point, we flush all meta data to the SCF HDF5 file. Then, the
   !   HDF5 file is primed for use.
   call h5fflush_f(scf_fid,H5F_SCOPE_GLOBAL_F,hdferr)

end subroutine initSCFIntegralHDF5


subroutine accessSCFIntegralHDF5 (scf_fid)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints,     only: numKPoints
   use O_Potential,   only: potDim
   use O_AtomicSites, only: valeDim

   ! Define the passed parameters.
   integer(hid_t) :: scf_fid

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
#ifndef GAMMA
   atomDims(1) = 2 ! Real and imaginary are needed.
#else
   atomDims(1) = 1 ! Real needed only
#endif
   atomDims(2) = valeDim*(valeDim+1)/2 ! Linear storage of 1/2 matrix.

   ! Open the Integral group within the scf HDF5 file.
   call h5gopen_f (scf_fid,"/atomIntgGroup",atomIntgGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom intg group'

   ! Open the subgroups within the atomIntgGroup. We open all the subgroups
   !   even if we don't need them so that the close process is still the
   !   same as the "init" case.
   call h5gopen_f (atomIntgGroup_gid,"atomOverlap",atomOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomKEOverlap",atomKEOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open kinetic energy overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomMVOverlap",atomMVOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open mass velocity overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomNucOverlap",atomNucOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom nuclear overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomDMOverlap",atomDMOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open dipole moment overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomMMOverlap",atomDMOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentum matrix overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomPotOverlap",atomPotOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom potential overlap group'

   ! For the integrals that have multiple matrices per kpoint, we need to
   !   create space to hold the group IDs of each kpoint group.
   allocate (atomDMxyzOL_gid (3)) ! Needs x,y,z matrices
   allocate (atomMMxyzOL_gid (3)) ! Needs x,y,z matrices
   allocate (atomPotTermOL_gid (potDim)) ! One matrix per pot term.

   ! Now, we open the xyz subgroups for the integrals that have them.
   do i = 1, 3
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5gopen_f (atomDMOverlap_gid,currentName,atomDMxyzOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open dipole moment xyz group'

      call h5gopen_f (atomMMOverlap_gid,currentName,atomMMxyzOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open momentum matrix xyz group'
   enddo

   ! Now, we can open the potential subgroups for each term.
   do i = 1, potDim
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5gopen_f (atomPotOverlap_gid,currentName,atomPotTermOL_gid(i),&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open potential term group'
   enddo

   ! Allocate space to hold the dataset IDs.
   allocate (atomOverlap_did (numKPoints))
   allocate (atomKEOverlap_did (numKPoints))
   allocate (atomMVOverlap_did (numKPoints))
   allocate (atomNucOverlap_did (numKPoints))
   allocate (atomDMOverlap_did (numKPoints,3))
   allocate (atomMMOverlap_did (numKPoints,3))
   allocate (atomPotOverlap_did (numKPoints,potDim))

   ! Allocate space to hold the attribute IDs for the potential terms.
   allocate (atomPotTermOL_aid (potDim))

   ! Open the datasets that will be used for all subgroups of atomIntgGroup.
   do i = 1, numKPoints
      write (currentName,fmt="(i7.7)") i
      call h5dopen_f (atomOverlap_gid,currentName,atomOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open atom overlap did'

      call h5dopen_f (atomKEOverlap_gid,currentName,atomKEOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open KE overlap did'

      call h5dopen_f (atomMVOverlap_gid,currentName,atomMVOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open MV overlap did'

      call h5dopen_f (atomNucOverlap_gid,currentName,atomNucOverlap_did(i),&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open nuclear overlap did'

      do j = 1, 3
         write (currentName,fmt="(i7.7)") j
         call h5dopen_f (atomDMxyzOL_gid(j),currentName,&
               & atomDMOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open DM overlap did'

         call h5dopen_f (atomMMxyzOL_gid(j),currentName,&
               & atomMMOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open MM overlap did'
      enddo

      do j = 1, potDim
         write (currentName,fmt="(i7.7)") j
         call h5dopen_f(atomPotTermOL_gid(j),currentName,&
               & atomPotOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open potential overlap did'
      enddo
   enddo

   ! Open the attributes of each integral group.
   call h5aopen_f (atomOverlap_gid,'status',atomOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom overlap aid'

   call h5aopen_f (atomKEOverlap_gid,'status',atomKEOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom KE overlap aid'

   call h5aopen_f (atomMVOverlap_gid,'status',atomMVOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom MV overlap aid'

   call h5aopen_f (atomNucOverlap_gid,'status',atomNucOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open nuclear overlap aid'

   call h5aopen_f (atomDMOverlap_gid,'status',atomDMOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open DM overlap aid'

   call h5aopen_f (atomMMOverlap_gid,'status',atomMMOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open MM overlap aid'

   do i = 1, potDim
      call h5aopen_f (atomPotTermOL_gid(i),'status',atomPotTermOL_aid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open potential overlap aid'
   enddo

   ! Obtain the properties of the datasets that were just opened.  They are all
   !   the same and so only one copy is necessary.  (Actually, this value is
   !   not really used, but in the "close" subroutine we close this id so we
   !   should make sure to have it for both the setup and main calls to the
   !   close routine.)
   call h5dget_create_plist_f (atomOverlap_did(1),valeVale_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain vale vale plid'

   ! Obtain the dataspace that is used for each dataset in atomIntgGroup
   !   and all of its subgroups.  The same dataspace definition works for all
   !   of the datasets.  (Same as for the plist above.)
   call h5dget_space_f(atomOverlap_did(1),valeVale_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain vale vale dsid'

end subroutine accessSCFIntegralHDF5


subroutine closeSCFIntegralHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints,   only: numKPoints
   use O_Potential, only: potDim

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: i,j
   integer :: hdferr

   ! Close the property list first.
   call h5pclose_f (valeVale_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeVale_plid.'

   ! Close the datasets next.
   do i = 1, numKPoints
      call h5dclose_f (atomOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomOverlap_did.'

      call h5dclose_f (atomKEOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomKEOverlap_did.'

      call h5dclose_f (atomMVOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomMVOverlap_did.'

      call h5dclose_f (atomNucOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomNucOverlap_did.'

      do j = 1, 3
         call h5dclose_f (atomDMOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomDMOverlap_did.'
      enddo

      do j = 1, 3
         call h5dclose_f (atomMMOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomMMOverlap_did.'
      enddo

      do j = 1, potDim
         call h5dclose_f (atomPotOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomPotOverlap_did.'
      enddo
   enddo

   ! Close the data spaces next.
   call h5sclose_f (valeVale_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeVale_dsid.'

   ! Close the groups.
   do i = 1, potDim
      call h5gclose_f (atomPotTermOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomPotKPointOL_gid.'
   enddo

   do i = 1, 3
      call h5gclose_f (atomDMxyzOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomDMxyzOL_gid.'
   enddo

   do i = 1, 3
      call h5gclose_f (atomMMxyzOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomMMxyzOL_gid.'
   enddo

   call h5gclose_f (atomOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomOverlap_gid.'

   call h5gclose_f (atomKEOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomKEOverlap_gid.'

   call h5gclose_f (atomMVOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomMVOverlap_gid.'

   call h5gclose_f (atomNucOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomNucOverlap_gid.'

   call h5gclose_f (atomPotOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomPotOverlap_gid.'

   call h5gclose_f (atomMMOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomDMxyzOverlap_gid.'

   call h5gclose_f (atomDMOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomDMxyzOverlap_gid.'

   call h5gclose_f (atomIntgGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomIntgGroup_gid.'

   ! The attributes tracking completion are closed when written.

   ! Deallocate the arrays that hold the id numbers.
   deallocate (atomPotTermOL_gid)
   deallocate (atomDMxyzOL_gid)
   deallocate (atomMMxyzOL_gid)
   deallocate (atomOverlap_did)
   deallocate (atomKEOverlap_did)
   deallocate (atomMVOverlap_did)
   deallocate (atomNucOverlap_did)
   deallocate (atomDMOverlap_did)
   deallocate (atomMMOverlap_did)
   deallocate (atomPotOverlap_did)

   ! Note that the attributes are closed as soon as they are finished.
   
end subroutine closeSCFIntegralHDF5

end module O_SCFIntegralsHDF5
