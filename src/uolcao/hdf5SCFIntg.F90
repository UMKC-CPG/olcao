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

   ! If a gamma k-point calculation is being done, then we only need one
   !   component. Otherwise, we need two.
   integer :: numComponents

   ! Define the main subgroup from scf_fid that holds the multicenter
   !   Gaussian integral results.
   integer(hid_t) :: atomIntgGroup_gid ! Atom and potential overlap group.

   ! Begin group, dataspace, dataset definitions for atomIntgGroup_gid.

   ! Define the group IDs under atomIntgGroup_gid. Each group here will hold
   !   either a dataset for each kpoint, or a group for each kpoint.
   integer(hid_t) :: atomOverlap_gid
   integer(hid_t) :: atomOverlapCV_gid ! For overlap coreVale
   integer(hid_t) :: atomKEOverlap_gid
   integer(hid_t) :: atomMVOverlap_gid
   integer(hid_t) :: atomNPOverlap_gid
   integer(hid_t) :: atomDMOverlap_gid
   integer(hid_t) :: atomMMOverlap_gid
   integer(hid_t) :: atomKOverlap_gid
   integer(hid_t) :: atomPotOverlap_gid

   ! Define the group IDs of the dynamically numbered subgroups. I.e., for
   !   the groups that hold another group, we define that collection of
   !   subgroups here.
   integer(hid_t), allocatable, dimension (:) :: atomDMxyzOL_gid
   integer(hid_t), allocatable, dimension (:) :: atomMMxyzOL_gid
   integer(hid_t), allocatable, dimension (:) :: atomKxyzOL_gid
   integer(hid_t), allocatable, dimension (:) :: atomPotTermOL_gid

   ! The dataspaces of each dataset in atomIntgGroup_gid are the same in
   !   all characteristics (type, dimension, etc.) and therefore can be
   !   given a static ID definition now.
   integer(hid_t) :: valeVale_dsid
   integer(hid_t) :: coreVale_dsid ! For overlap coreVale only.

   ! Each of the below datasets for atomIntgGroup will use the same property
   !   list.
   integer(hid_t) :: valeVale_plid
   integer(hid_t) :: coreVale_plid ! For overlap coreVale only.

   ! Define array that holds the dimensions of the dataset.
   integer(hsize_t), dimension (2) :: packedVVDims
   integer(hsize_t), dimension (2) :: fullCVDims

   ! Define array that holds the dimensions of the chunk.
   integer (hsize_t), dimension (2) :: packedVVDimsChunk
   integer (hsize_t), dimension (2) :: fullCVDimsChunk

   ! The number of datasets under each of the atomIntgGroup groups will
   !   vary depending on the problem.  (# of kpoints, potential dimension,
   !   or xyz terms). Thus, the datasets must be given IDs dynamically.
   integer(hid_t), allocatable, dimension (:)   :: atomOverlap_did
   integer(hid_t), allocatable, dimension (:,:) :: atomOverlapCV_did !coreVale
   integer(hid_t), allocatable, dimension (:)   :: atomKEOverlap_did
   integer(hid_t), allocatable, dimension (:)   :: atomMVOverlap_did
   integer(hid_t), allocatable, dimension (:)   :: atomNPOverlap_did
   integer(hid_t), allocatable, dimension (:,:) :: atomDMOverlap_did
   integer(hid_t), allocatable, dimension (:,:) :: atomMMOverlap_did
   integer(hid_t), allocatable, dimension (:,:) :: atomKOverlap_did
   integer(hid_t), allocatable, dimension (:,:) :: atomPotOverlap_did

   ! Define the attribute IDs that will be used for each group. One attribute
   !   is sufficient for each integral type except the potential overlap.
   !   For that one, we will need a set of attributes, one for each potential
   !   term.
   integer(hid_t) :: atomOverlap_aid ! No CV needed. Regular accounts for both.
   integer(hid_t) :: atomKEOverlap_aid
   integer(hid_t) :: atomMVOverlap_aid
   integer(hid_t) :: atomNPOverlap_aid
   integer(hid_t) :: atomDMOverlap_aid
   integer(hid_t) :: atomMMOverlap_aid
   integer(hid_t) :: atomKOverlap_aid
   integer(hid_t), allocatable, dimension (:) :: atomPotTermOL_aid


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initSCFIntegralHDF5 (scf_fid, attribInt_dsid, attribIntDims)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints,     only: numKPoints
   use O_Potential,   only: potDim
   use O_AtomicSites, only: coreDim, valeDim

   ! Define the passed parameters.
   integer(hid_t) :: scf_fid
   integer(hid_t) :: attribInt_dsid
   integer(hsize_t), dimension (1) :: attribIntDims

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Determine the number of components for some data structures.
#ifndef GAMMA
   numComponents = 2 ! Real and imaginary are needed.
#else
   numComponents = 1 ! Real only is needed.
#endif

   ! Initialize data structure dimensions.
   packedVVDims(1) = numComponents
   packedVVDims(2) = valeDim*(valeDim+1)/2 ! Linear storage of 1/2 matrix.
   if (coreDim > 0) then
      fullCVDims(1) = coreDim
      fullCVDims(2) = valeDim
   endif

   ! Check that the chunk size is not too large. The assumption here is that
   !   the number being stored are 8 byte reals and that we should not go over
   !   2 billion bytes. Note that if x*y = >250M (where x=packedVVDims(1) and
   !   y=packedVVDims(2)) and we want a*b = 250M (where a=packedVVDimsChunk(1)
   !   and b=packedVVDimsChunk(2)) and we require x=a, then b = int(250M/a).
   if (packedVVDims(1) * packedVVDims(2) > 250000000) then
      packedVVDimsChunk(2) = int(250000000/packedVVDims(1))
      packedVVDimsChunk(1) = packedVVDims(1)
   else
      packedVVDimsChunk(1) = packedVVDims(1)
      packedVVDimsChunk(2) = packedVVDims(2)
   endif

   ! Doing it for fullCV, we now have x*y = >250M and a*b for the CV case. We
   !   will ignore the "first" index which can be 1 or 2 for real only or real
   !   plus imaginary because we will always store both separately (if both
   !   real and imaginary are needed). So a*b = 250M. Here, we want to keep
   !   the x/y and a/b ratio equal. So: x*y = >250M. We require a*b = 250M so
   !   we say x*y / q = a*b = 250M. Thus, q = x*y / 250M. To keep the x/y to
   !   a/b ratio, we use q=p^2 and get: p = sqrt(x*y / 250M). Then we compute
   !   a=x/p and b=y/p. Written another way: a=x/sqrt(x*y / 250M) and
   !   b=y/sqrt(x*y / 250M), or a=x/sqrt(x*y/250M) b=y/sqrt(x*y/250M).
   if (coreDim > 0) then
      if (product(fullCVDims(:)) > 250000000) then
         fullCVDimsChunk(1) = int(fullCVDims(1) / &
               & sqrt(product(fullCVDims(:)) / 250000000.0d0))
         fullCVDimsChunk(2) = int(fullCVDims(2) / &
               & sqrt(product(fullCVDims(:)) / 250000000.0d0))
      else
         fullCVDimsChunk(1) = fullCVDims(1)
         fullCVDimsChunk(2) = fullCVDims(2)
      endif
   endif

   ! Create the Integral group within the scf HDF5 file.
   call h5gcreate_f (scf_fid,"/atomIntgGroup",atomIntgGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom intg group SCF'

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
   if (coreDim > 0) then
      call h5gcreate_f (atomIntgGroup_gid,"atomOverlapCV",atomOverlapCV_gid,&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to create atom overlapCV group'
   endif

   call h5gcreate_f (atomIntgGroup_gid,"atomOverlap",atomOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom overlap group'

   call h5gcreate_f (atomIntgGroup_gid,"atomKEOverlap",atomKEOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create kinetic energy overlap group'

   call h5gcreate_f (atomIntgGroup_gid,"atomMVOverlap",atomMVOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create mass velocity overlap group'

   call h5gcreate_f (atomIntgGroup_gid,"atomNPOverlap",atomNPOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create atom nuclear overlap group'

   call h5gcreate_f (atomIntgGroup_gid,"atomDMOverlap",atomDMOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create dipole moment overlap group'

   call h5gcreate_f (atomIntgGroup_gid,"atomMMOverlap",atomMMOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create momentum matrix overlap group'

   call h5gcreate_f (atomIntgGroup_gid,"atomKOverlap",atomKOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create Koverlap group'

   call h5gcreate_f (atomIntgGroup_gid,"atomPotOverlap",atomPotOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create atom potential overlap group'

   ! For the integrals that have multiple different types of matrices, we
   !   need to create space to hold the group IDs of each matrix.
   allocate (atomDMxyzOL_gid (3)) ! Needs x,y,z matrices all done together
   allocate (atomMMxyzOL_gid (3)) ! Needs x,y,z matrices all done together
   allocate (atomKxyzOL_gid (3)) ! Needs x,y,z matrices all done together
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

      call h5gcreate_f (atomKOverlap_gid,currentName,&
            & atomKxyzOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to create KOverlap xyz group'
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
   if (coreDim > 0) then
      allocate (atomOverlapCV_did (numComponents,numKPoints))
   else
      allocate (atomOverlapCV_did (1,1)) ! Unused, but needs to be allocated.
   endif
   allocate (atomOverlap_did    (numKPoints))
   allocate (atomKEOverlap_did  (numKPoints))
   allocate (atomMVOverlap_did  (numKPoints))
   allocate (atomNPOverlap_did  (numKPoints))
   allocate (atomDMOverlap_did  (numKPoints,3))
   allocate (atomMMOverlap_did  (numKPoints,3))
   allocate (atomKOverlap_did  (numKPoints,3))
   allocate (atomPotOverlap_did (numKPoints,potDim))

   ! Allocate space to hold the attribute IDs for the potential matrics.
   allocate (atomPotTermOL_aid  (potDim))

   ! Create the dataspaces that will be used for each dataset in atomIntgGroup
   !   and all of its subgroups.  The same dataspace definition works for all
   !   of the VV datasets.
   call h5screate_simple_f(2,packedVVDims,valeVale_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create packed vale vale dsid'
   if (coreDim > 0) then
      call h5screate_simple_f(2,fullCVDims,coreVale_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create full core vale dsid'
   endif

   ! Define the properties of the datasets to be made.

   ! Create the VV property list first.  Then set the properties one at a time.
   call h5pcreate_f      (H5P_DATASET_CREATE_F,valeVale_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create vale vale plid'
   call h5pset_layout_f  (valeVale_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale plid layout'
   call h5pset_chunk_f   (valeVale_plid,2,packedVVDimsChunk,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale plid chunk size'
!   call h5pset_shuffle_f (valeVale_plid,hdferr)
   call h5pset_deflate_f (valeVale_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale for deflation'

   ! Create the CV property list first.  Then set the properties one at a time.
   if (coreDim > 0) then
      call h5pcreate_f      (H5P_DATASET_CREATE_F,coreVale_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to create core vale plid'
      call h5pset_layout_f  (coreVale_plid,H5D_CHUNKED_F,hdferr)
      if (hdferr /= 0) stop 'Failed to set core vale plid layout'
      call h5pset_chunk_f   (coreVale_plid,2,fullCVDimsChunk,hdferr)
      if (hdferr /= 0) stop 'Failed to set core vale plid chunk size'
!      call h5pset_shuffle_f (coreVale_plid,hdferr)
      call h5pset_deflate_f (coreVale_plid,1,hdferr)
      if (hdferr /= 0) stop 'Failed to set core vale for deflation'
   endif

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

      call h5dcreate_f (atomNPOverlap_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,atomNPOverlap_did(i),hdferr,valeVale_plid)
      if (hdferr /= 0) stop 'Failed to create nuclear overlap did'

      do j = 1, 3
         call h5dcreate_f (atomDMxyzOL_gid(j),currentName,&
               & H5T_NATIVE_DOUBLE,valeVale_dsid,atomDMOverlap_did(i,j),&
               & hdferr,valeVale_plid)
         if (hdferr /= 0) stop 'Failed to create DM overlap did'

         call h5dcreate_f (atomMMxyzOL_gid(j),currentName,&
               & H5T_NATIVE_DOUBLE,valeVale_dsid,atomMMOverlap_did(i,j),&
               & hdferr,valeVale_plid)
         if (hdferr /= 0) stop 'Failed to create MM overlap did'

         call h5dcreate_f (atomKxyzOL_gid(j),currentName,&
               & H5T_NATIVE_DOUBLE,valeVale_dsid,atomKOverlap_did(i,j),&
               & hdferr,valeVale_plid)
         if (hdferr /= 0) stop 'Failed to create Koverlap did'
      enddo

      do j = 1, potDim
         call h5dcreate_f(atomPotTermOL_gid(j),currentName,&
               & H5T_NATIVE_DOUBLE,valeVale_dsid,atomPotOverlap_did(i,j),&
               & hdferr,valeVale_plid)
         if (hdferr /= 0) stop 'Failed to create potential overlap did'
      enddo

      if (coreDim > 0) then
         do j = 1, numComponents
            if (j == 1) write (currentName,fmt="(a4,i7.7)") "real",i
            if (j == 2) write (currentName,fmt="(a4,i7.7)") "imag",i
            currentName = trim (currentName)
            call h5dcreate_f (atomOverlapCV_gid,currentName,H5T_NATIVE_DOUBLE,&
                  & coreVale_dsid,atomOverlapCV_did(j,i),hdferr,coreVale_plid)
            if (hdferr /= 0) stop 'Failed to create atom overlapCV did'
         enddo
      endif
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

   call h5acreate_f (atomNPOverlap_gid,"status",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,atomNPOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create Nuc overlap aid'

   call h5acreate_f (atomDMOverlap_gid,"status",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,atomDMOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create DM overlap aid'

   call h5acreate_f (atomMMOverlap_gid,"status",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,atomMMOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create MM overlap aid'

   call h5acreate_f (atomKOverlap_gid,"status",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,atomKOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create Koverlap aid'

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

   call h5awrite_f(atomNPOverlap_aid,H5T_NATIVE_INTEGER,0,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomNPOverlap_aid'

   call h5awrite_f(atomDMOverlap_aid,H5T_NATIVE_INTEGER,0,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomDMOverlap_aid'

   call h5awrite_f(atomMMOverlap_aid,H5T_NATIVE_INTEGER,0,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomMMOverlap_aid'

   call h5awrite_f(atomKOverlap_aid,H5T_NATIVE_INTEGER,0,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomKOverlap_aid'

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
   use O_AtomicSites, only: coreDim, valeDim

   ! Define the passed parameters.
   integer(hid_t) :: scf_fid

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Determine the number of components for some data structures.
#ifndef GAMMA
   numComponents = 2 ! Real and imaginary are needed.
#else
   numComponents = 1 ! Real only is needed.
#endif

   ! Initialize data structure dimensions.
   packedVVDims(1) = numComponents
   packedVVDims(2) = valeDim*(valeDim+1)/2 ! Linear storage of 1/2 matrix.
   if (coreDim > 0) then
      fullCVDims(1) = coreDim
      fullCVDims(2) = valeDim
   endif

   ! Open the Integral group within the scf HDF5 file.
   call h5gopen_f (scf_fid,"/atomIntgGroup",atomIntgGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom intg group'

   ! Open the subgroups within the atomIntgGroup. We open all the subgroups
   !   even if we don't need them so that the close process is still the
   !   same as the "init" case.
   if (coreDim > 0) then
      call h5gopen_f (atomIntgGroup_gid,"atomOverlapCV",atomOverlapCV_gid,&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open atom overlapCV group'
   endif
   call h5gopen_f (atomIntgGroup_gid,"atomOverlap",atomOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomKEOverlap",atomKEOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open kinetic energy overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomMVOverlap",atomMVOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open mass velocity overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomNPOverlap",atomNPOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom nuclear overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomDMOverlap",atomDMOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open dipole moment overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomMMOverlap",atomMMOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentum matrix overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomKOverlap",atomKOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open Koverlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomPotOverlap",atomPotOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom potential overlap group'

   ! For the integrals that have multiple matrices per kpoint, we need to
   !   create space to hold the group IDs of each kpoint group.
   allocate (atomDMxyzOL_gid (3)) ! Needs x,y,z matrices
   allocate (atomMMxyzOL_gid (3)) ! Needs x,y,z matrices
   allocate (atomKxyzOL_gid (3)) ! Needs x,y,z matrices
   allocate (atomPotTermOL_gid (potDim)) ! One matrix per pot term.

   ! Now, we open the xyz subgroups for the integrals that have them.
   do i = 1, 3
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5gopen_f (atomDMOverlap_gid,currentName,atomDMxyzOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open dipole moment xyz group'

      call h5gopen_f (atomMMOverlap_gid,currentName,atomMMxyzOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open momentum matrix xyz group'

      call h5gopen_f (atomKOverlap_gid,currentName,atomKxyzOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open KOverlap xyz group'
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
   if (coreDim > 0) then
      allocate (atomOverlapCV_did (numComponents,numKPoints))
   else
      allocate (atomOverlapCV_did (1,1)) ! Unused, but needs to be allocated.
   endif
   allocate (atomOverlap_did (numKPoints))
   allocate (atomKEOverlap_did (numKPoints))
   allocate (atomMVOverlap_did (numKPoints))
   allocate (atomNPOverlap_did (numKPoints))
   allocate (atomDMOverlap_did (numKPoints,3))
   allocate (atomMMOverlap_did (numKPoints,3))
   allocate (atomKOverlap_did (numKPoints,3))
   allocate (atomPotOverlap_did (numKPoints,potDim))

   ! Allocate space to hold the attribute IDs for the potential terms.
   allocate (atomPotTermOL_aid (potDim))

   ! Open the datasets that will be used for all subgroups of atomIntgGroup.
   do i = 1, numKPoints
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5dopen_f (atomOverlap_gid,currentName,atomOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open atom overlap did'

      call h5dopen_f (atomKEOverlap_gid,currentName,atomKEOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open KE overlap did'

      call h5dopen_f (atomMVOverlap_gid,currentName,atomMVOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open MV overlap did'

      call h5dopen_f (atomNPOverlap_gid,currentName,atomNPOverlap_did(i),&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open nuclear overlap did'

      do j = 1, 3
         call h5dopen_f (atomDMxyzOL_gid(j),currentName,&
               & atomDMOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open DM overlap did'

         call h5dopen_f (atomMMxyzOL_gid(j),currentName,&
               & atomMMOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open MM overlap did'

         call h5dopen_f (atomKxyzOL_gid(j),currentName,&
               & atomKOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open Koverlap did'
      enddo

      do j = 1, potDim
         call h5dopen_f(atomPotTermOL_gid(j),currentName,&
               & atomPotOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open potential overlap did'
      enddo

      if (coreDim > 0) then
         do j = 1, numComponents
            if (j == 1) write (currentName,fmt="(a4,i7.7)") "real",i
            if (j == 2) write (currentName,fmt="(a4,i7.7)") "imag",i
            currentName = trim (currentName)
            call h5dopen_f (atomOverlapCV_gid,currentName,&
                  & atomOverlapCV_did(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to create atom overlapCV did'
         enddo
      endif
   enddo

   ! Open the attributes of each integral group.
   call h5aopen_f (atomOverlap_gid,'status',atomOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom overlap aid'

   call h5aopen_f (atomKEOverlap_gid,'status',atomKEOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom KE overlap aid'

   call h5aopen_f (atomMVOverlap_gid,'status',atomMVOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom MV overlap aid'

   call h5aopen_f (atomNPOverlap_gid,'status',atomNPOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open nuclear overlap aid'

   call h5aopen_f (atomDMOverlap_gid,'status',atomDMOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open DM overlap aid'

   call h5aopen_f (atomMMOverlap_gid,'status',atomMMOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open MM overlap aid'

   call h5aopen_f (atomKOverlap_gid,'status',atomKOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open Koverlap aid'

   do i = 1, potDim
      call h5aopen_f (atomPotTermOL_gid(i),'status',atomPotTermOL_aid(i),&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open potential overlap aid'
   enddo

   ! Obtain the properties of the datasets that were just opened.  All VV are
   !   the same and so only one copy is necessary.  (Actually, this value is
   !   not really used, but in the "close" subroutine we close this id so we
   !   should make sure to have it for both the setup and main calls to the
   !   close routine.) Also open one for the CV data.
   call h5dget_create_plist_f (atomOverlap_did(1),valeVale_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain vale vale plid'
   if (coreDim > 0) then
      call h5dget_create_plist_f (atomOverlapCV_did(1,1),coreVale_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to obtain core vale plid'
   endif

   ! Obtain the dataspace that is used for each VV dataset in atomIntgGroup
   !   and all of its subgroups.  The same dataspace definition works for all
   !   of the datasets.  (Same as for the plist above.)
   call h5dget_space_f(atomOverlap_did(1),valeVale_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain vale vale dsid'
   if (coreDim > 0) then
      call h5dget_space_f(atomOverlapCV_did(1,1),coreVale_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to obtain core vale dsid'
   endif

end subroutine accessSCFIntegralHDF5


subroutine closeSCFIntegralHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints
   use O_Potential, only: potDim
   use O_AtomicSites, only: coreDim

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

      call h5dclose_f (atomNPOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomNPOverlap_did.'

      do j = 1, 3
         call h5dclose_f (atomDMOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomDMOverlap_did.'
      enddo

      do j = 1, 3
         call h5dclose_f (atomMMOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomMMOverlap_did.'
      enddo

      do j = 1, 3
         call h5dclose_f (atomKOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomKOverlap_did.'
      enddo

      do j = 1, potDim
         call h5dclose_f (atomPotOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomPotOverlap_did.'
      enddo

      if (coreDim > 0) then
         do j = 1, numComponents
            call h5dclose_f (atomOverlapCV_did(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to close atomOverlapCV_did.'
         enddo
      endif

   enddo

   ! Close the data spaces next.
   call h5sclose_f (valeVale_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeVale_dsid.'
   if (coreDim > 0) then
      call h5sclose_f (coreVale_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close coreVale_dsid.'
   endif

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

   do i = 1, 3
      call h5gclose_f (atomKxyzOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomKxyzOL_gid.'
   enddo

   if (coreDim > 0) then
      call h5gclose_f (atomOverlapCV_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomOverlapCV_gid.'
   endif

   call h5gclose_f (atomOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomOverlap_gid.'

   call h5gclose_f (atomKEOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomKEOverlap_gid.'

   call h5gclose_f (atomMVOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomMVOverlap_gid.'

   call h5gclose_f (atomNPOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomNPOverlap_gid.'

   call h5gclose_f (atomPotOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomPotOverlap_gid.'

   call h5gclose_f (atomDMOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomDMxyzOverlap_gid.'

   call h5gclose_f (atomMMOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomMMxyzOverlap_gid.'

   call h5gclose_f (atomKOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomKxyzOverlap_gid.'

   call h5gclose_f (atomIntgGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomIntgGroup_gid.'

   ! The attributes tracking completion are closed when written.

   ! Deallocate the arrays that hold the id numbers.
   deallocate (atomOverlapCV_did) ! Sometimes unused. Always allocated.
   deallocate (atomOverlap_did)
   deallocate (atomKEOverlap_did)
   deallocate (atomMVOverlap_did)
   deallocate (atomNPOverlap_did)
   deallocate (atomDMOverlap_did)
   deallocate (atomMMOverlap_did)
   deallocate (atomKOverlap_did)
   deallocate (atomPotOverlap_did)
   deallocate (atomPotTermOL_aid)
   deallocate (atomPotTermOL_gid)
   deallocate (atomDMxyzOL_gid)
   deallocate (atomMMxyzOL_gid)
   deallocate (atomKxyzOL_gid)

   ! Note that the attributes are closed as soon as they are finished.
   
end subroutine closeSCFIntegralHDF5

end module O_SCFIntegralsHDF5
