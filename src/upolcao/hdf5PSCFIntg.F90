module O_PSCFIntegralsHDF5

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

   ! Define the main subgroup from pscf_fid that holds the multicenter
   !   Gaussian integral results.
   integer(hid_t) :: atomIntgGroupPSCF_gid ! Atom and potential overlap group.

   ! Begin group, dataspace, dataset definitions for atomIntgGroupPSCF_gid.

   ! Define the group IDs under atomIntgGroupPSCF_gid. Each group here will
   !   hold either a dataset for each kpoint, or a group for each kpoint. The
   !   SYBD versions are for the set of kpoints in symmetric band structure
   !   diagrams.
   integer(hid_t) :: atomOverlapPSCF_gid
   integer(hid_t) :: atomOverlapCV_PSCF_gid
   integer(hid_t) :: atomHamOverlapPSCF_gid
   integer(hid_t) :: atomDMOverlapPSCF_gid
   integer(hid_t) :: atomMMOverlapPSCF_gid
   integer(hid_t) :: atomOverlapSYBD_PSCF_gid
   integer(hid_t) :: atomOverlapCV_SYBD_PSCF_gid
   integer(hid_t) :: atomHamOverlapSYBD_PSCF_gid
   integer(hid_t) :: atomDMOverlapSYBD_PSCF_gid
   integer(hid_t) :: atomMMOverlapSYBD_PSCF_gid

   ! Define the group IDs of the dynamically numbered subgroups. I.e., for
   !   the groups that hold another group, we define that collection of
   !   subgroups here.
   integer(hid_t), allocatable, dimension (:) :: atomDMxyzOL_PSCF_gid
   integer(hid_t), allocatable, dimension (:) :: atomMMxyzOL_PSCF_gid
   integer(hid_t), allocatable, dimension (:) :: atomDMxyzOL_SYBD_PSCF_gid
   integer(hid_t), allocatable, dimension (:) :: atomMMxyzOL_SYBD_PSCF_gid

   ! The dataspaces of each dataset in atomIntgGroupPSCF_gid are the same in
   !   all characteristics (type, dimension, etc.) and therefore can be
   !   given a static ID definition now.
   integer(hid_t) :: valeValePSCF_dsid
   integer(hid_t) :: coreValePSCF_dsid ! For overlap coreVale only.

   ! Each of the below datasets for atomIntgGroup will use the same property
   !   list.
   integer(hid_t) :: valeValePSCF_plid
   integer(hid_t) :: coreValePSCF_plid ! For overlap coreVale only.

   ! Define array that holds the dimensions of the dataset.
   integer(hsize_t), dimension (2) :: packedVVDimsPSCF
   integer(hsize_t), dimension (2) :: fullCVDimsPSCF

   ! Define array that holds the dimensions of the chunk.
   integer (hsize_t), dimension (2) :: packedVVDimsPSCFChunk
   integer (hsize_t), dimension (2) :: fullCVDimsPSCFChunk

   ! The number of datasets under each of the atomIntgGroup groups will
   !   vary depending on the problem.  (# of kpoints, spin, or xyz terms).
   !   Thus, the datasets must be given IDs dynamically.
   integer(hid_t), allocatable, dimension (:)   :: atomOverlapPSCF_did
   integer(hid_t), allocatable, dimension (:,:) :: atomOverlapCV_PSCF_did
   integer(hid_t), allocatable, dimension (:,:) :: atomHamOverlapPSCF_did
   integer(hid_t), allocatable, dimension (:,:) :: atomDMOverlapPSCF_did
   integer(hid_t), allocatable, dimension (:,:) :: atomMMOverlapPSCF_did
   integer(hid_t), allocatable, dimension (:)   :: atomOverlapSYBD_PSCF_did
   integer(hid_t), allocatable, dimension (:,:) :: atomOverlapCV_SYBD_PSCF_did
   integer(hid_t), allocatable, dimension (:,:) :: atomHamOverlapSYBD_PSCF_did
   integer(hid_t), allocatable, dimension (:,:) :: atomDMOverlapSYBD_PSCF_did
   integer(hid_t), allocatable, dimension (:,:) :: atomMMOverlapSYBD_PSCF_did

   ! Define the attribute IDs that will be used for each group. One attribute
   !   is sufficient for each integral type except the potential overlap.
   !   For that one, we will need a set of attributes, one for each potential
   !   term.
   integer(hid_t) :: atomOverlapPSCF_aid ! No CV aid needed. This does both.
   integer(hid_t) :: atomHamOverlapPSCF_aid
   integer(hid_t) :: atomDMOverlapPSCF_aid
   integer(hid_t) :: atomMMOverlapPSCF_aid
   integer(hid_t) :: atomOverlapSYBD_PSCF_aid
   integer(hid_t) :: atomHamOverlapSYBD_PSCF_aid
   integer(hid_t) :: atomDMOverlapSYBD_PSCF_aid
   integer(hid_t) :: atomMMOverlapSYBD_PSCF_aid


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initPSCFIntegralHDF5 (pscf_fid, attribIntPSCF_dsid,&
     & attribIntDimsPSCF)

   ! Import any necessary definition modules.
   use HDF5
   use O_MPI

   ! Import necessary object modules.
   use O_Potential, only: spin
   use O_KPoints, only: numKPoints_HDF5, numPathKP_HDF5
   use O_AtomicSites, only: coreDim, valeDim

   ! Define the passed parameters.
   integer(hid_t) :: pscf_fid
   integer(hid_t) :: attribIntPSCF_dsid
   integer(hsize_t), dimension (1) :: attribIntDimsPSCF

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
   packedVVDimsPSCF(1) = numComponents
   packedVVDimsPSCF(2) = valeDim*(valeDim+1)/2 ! Linear storage of 1/2 matrix.
   if (coreDim > 0) then
      fullCVDimsPSCF(1) = coreDim
      fullCVDimsPSCF(2) = valeDim
   endif

   ! Check that the chunk size is not too large. The assumption here is that
   !   the number being stored are 8 byte reals and that we should not go over
   !   2 billion bytes. Note that if x*y = >250M (where x=packedVVDimsPSCF(1)
   !   and y=packedVVDimsPSCF(2)) and we want a*b = 250M (where
   !   a=packedVVDimsPSCFChunk(1) and b=packedVVDimsPSCFChunk(2)) and we
   !   require x=a, then b = int(250M/a).
   if (packedVVDimsPSCF(1) * packedVVDimsPSCF(2) > 250000000) then
      packedVVDimsPSCFChunk(2) = int(250000000/packedVVDimsPSCF(1))
      packedVVDimsPSCFChunk(1) = packedVVDimsPSCF(1)
   else
      packedVVDimsPSCFChunk(1) = packedVVDimsPSCF(1)
      packedVVDimsPSCFChunk(2) = packedVVDimsPSCF(2)
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
      if (product(fullCVDimsPSCF(:)) > 250000000) then
         fullCVDimsPSCFChunk(1) = int(fullCVDimsPSCF(1) / &
               & sqrt(product(fullCVDimsPSCF(:)) / 250000000.0d0))
         fullCVDimsPSCFChunk(2) = int(fullCVDimsPSCF(2) / &
               & sqrt(product(fullCVDimsPSCF(:)) / 250000000.0d0))
      else
         fullCVDimsPSCFChunk(1) = fullCVDimsPSCF(1)
         fullCVDimsPSCFChunk(2) = fullCVDimsPSCF(2)
      endif
   endif

   ! For the integrals that have multiple different types of matrices, we
   !   need to create space to hold the group IDs of each matrix.
   allocate (atomDMxyzOL_PSCF_gid (3)) ! x,y,z matrices all done together
   allocate (atomMMxyzOL_PSCF_gid (3)) ! x,y,z matrices all done together
   allocate (atomDMxyzOL_SYBD_PSCF_gid (3)) ! x,y,z matrices all done together
   allocate (atomMMxyzOL_SYBD_PSCF_gid (3)) ! x,y,z matrices all done together

   ! Sufficient space must be allocated to hold the dataset IDs for all
   !   matrices that need to be computed.
   if (coreDim > 0) then
      allocate (atomOverlapCV_PSCF_did  (numComponents,numKPoints_HDF5))
   else
      allocate (atomOverlapCV_PSCF_did  (1,1))
   endif
   allocate (atomOverlapPSCF_did     (numKPoints_HDF5))
   allocate (atomHamOverlapPSCF_did  (numKPoints_HDF5,spin))
   allocate (atomDMOverlapPSCF_did   (numKPoints_HDF5,3))
   allocate (atomMMOverlapPSCF_did   (numKPoints_HDF5,3))

   if (coreDim > 0) then
      allocate (atomOverlapCV_SYBD_PSCF_did  (numComponents,numPathKP_HDF5))
   else
      allocate (atomOverlapCV_SYBD_PSCF_did  (1,1))
   endif
   allocate (atomOverlapSYBD_PSCF_did     (numPathKP_HDF5))
   allocate (atomHamOverlapSYBD_PSCF_did  (numPathKP_HDF5,spin))
   allocate (atomDMOverlapSYBD_PSCF_did   (numPathKP_HDF5,3))
   allocate (atomMMOverlapSYBD_PSCF_did   (numPathKP_HDF5,3))

   ! Only process 0 opens the HDF5 structure.
   if (mpiRank /= 0) return

   ! Create the Integral group within the pscf HDF5 file.
   call h5gcreate_f (pscf_fid,"/atomIntgGroup",atomIntgGroupPSCF_gid,hdferr)
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
   !   pscf, we create all groups anyway so that if a future olcao pscf call
   !   does request that group, it will already be available to access.
   call h5gcreate_f (atomIntgGroupPSCF_gid,"atomOverlap",&
         & atomOverlapPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom overlap group'

   call h5gcreate_f (atomIntgGroupPSCF_gid,"atomHamOverlap",&
         & atomHamOverlapPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom hamiltonian overlap group'

   if (coreDim > 0) then
      call h5gcreate_f (atomIntgGroupPSCF_gid,"atomOverlapCV",&
            & atomOverlapCV_PSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to create atom overlapCV group'
   endif

   call h5gcreate_f (atomIntgGroupPSCF_gid,"atomDMOverlap",&
         & atomDMOverlapPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create dipole moment overlap group'

   call h5gcreate_f (atomIntgGroupPSCF_gid,"atomMMOverlap",&
         & atomMMOverlapPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create momentum matrix overlap group'

   call h5gcreate_f (atomIntgGroupPSCF_gid,"atomOverlapSYBD",&
         & atomOverlapSYBD_PSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom overlap SYBD group'

   if (coreDim > 0) then
      call h5gcreate_f (atomIntgGroupPSCF_gid,"atomOverlapCV_SYBD",&
            & atomOverlapCV_SYBD_PSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to create atom overlapCV SYBD group'
   endif

   call h5gcreate_f (atomIntgGroupPSCF_gid,"atomHamOverlapSYBD",&
         & atomHamOverlapSYBD_PSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom hamiltonian overlap SYBD group'

   call h5gcreate_f (atomIntgGroupPSCF_gid,"atomDMOverlapSYBD",&
         & atomDMOverlapSYBD_PSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create dipole moment overlap SYBD group'

   call h5gcreate_f (atomIntgGroupPSCF_gid,"atomMMOverlapSYBD",&
         & atomMMOverlapSYBD_PSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create momentum matrix overlap SYBD group'

   ! Create the subgroups for the DM and MM terms and assign a gid name
   !   following 1=x, 2=y, and 3=z.
   do i = 1, 3
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5gcreate_f (atomDMOverlapPSCF_gid,currentName,&
            & atomDMxyzOL_PSCF_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to create dipole moment xyz group'

      call h5gcreate_f (atomMMOverlapPSCF_gid,currentName,&
            & atomMMxyzOL_PSCF_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to create momentum matrix xyz group'

      call h5gcreate_f (atomDMOverlapSYBD_PSCF_gid,currentName,&
            & atomDMxyzOL_SYBD_PSCF_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to create dipole moment xyz SYBD group'

      call h5gcreate_f (atomMMOverlapSYBD_PSCF_gid,currentName,&
            & atomMMxyzOL_SYBD_PSCF_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to create momentum matrix xyz SYBD group'
   enddo

   ! Create the dataspace that will be used for each dataset in atomIntgGroup
   !   and all of its subgroups.  The same dataspace definition works for all
   !   of the datasets.
   call h5screate_simple_f(2,packedVVDimsPSCF,valeValePSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create vale vale dsid'
   if (coreDim > 0) then
      call h5screate_simple_f(2,fullCVDimsPSCF,coreValePSCF_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create full core vale dsid'
   endif

   ! Define the properties of the datasets to be made.

   ! Create the VV property list.  Then set the properties one at a time.
   call h5pcreate_f      (H5P_DATASET_CREATE_F,valeValePSCF_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create vale vale plid'
   call h5pset_layout_f  (valeValePSCF_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale plid layout'
   call h5pset_chunk_f   (valeValePSCF_plid,2,packedVVDimsPSCFChunk,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale plid chunk size'
!   call h5pset_shuffle_f (valeValePSCF_plid,hdferr)
   call h5pset_deflate_f (valeValePSCF_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale for deflation'

   ! Create the CV property list.  Then set the properties one at a time.
   if (coreDim > 0) then
      call h5pcreate_f      (H5P_DATASET_CREATE_F,coreValePSCF_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to create core vale plid'
      call h5pset_layout_f  (coreValePSCF_plid,H5D_CHUNKED_F,hdferr)
      if (hdferr /= 0) stop 'Failed to set core vale plid layout'
      call h5pset_chunk_f   (coreValePSCF_plid,2,fullCVDimsPSCFChunk,hdferr)
      if (hdferr /= 0) stop 'Failed to set core vale plid chunk size'
!      call h5pset_shuffle_f (coreValePSCF_plid,hdferr)
      call h5pset_deflate_f (coreValePSCF_plid,1,hdferr)
      if (hdferr /= 0) stop 'Failed to set core vale for deflation'
   endif

   ! Define the dataspace of the computed "status" attribute that will be used
   !   for each group.
   call h5screate_simple_f (1,attribIntDimsPSCF(1),attribIntPSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create the attribIntPSCF_dsid'

   ! Create the datasets that will be used for all subgroups of atomIntgGroup.
   do i = 1, numKPoints_HDF5
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5dcreate_f (atomOverlapPSCF_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeValePSCF_dsid,atomOverlapPSCF_did(i),hdferr,&
            & valeValePSCF_plid)
      if (hdferr /= 0) stop 'Failed to create atom overlap did'

      do j = 1, spin
         call h5dcreate_f (atomHamOverlapPSCF_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeValePSCF_dsid,atomHamOverlapPSCF_did(i,j),hdferr,&
               & valeValePSCF_plid)
         if (hdferr /= 0) stop 'Failed to create hamiltonian overlap did'
      enddo

      do j = 1, 3
         call h5dcreate_f (atomDMxyzOL_PSCF_gid(j),currentName,&
               & H5T_NATIVE_DOUBLE,valeValePSCF_dsid,&
               & atomDMOverlapPSCF_did(i,j),hdferr,valeValePSCF_plid)
         if (hdferr /= 0) stop 'Failed to create DM overlap did'

         call h5dcreate_f (atomMMxyzOL_PSCF_gid(j),currentName,&
               & H5T_NATIVE_DOUBLE,valeValePSCF_dsid,&
               & atomMMOverlapPSCF_did(i,j),hdferr,valeValePSCF_plid)
         if (hdferr /= 0) stop 'Failed to create MM overlap did'
      enddo

      if (coreDim > 0) then
         do j = 1, numComponents
            if (j == 1) write (currentName,fmt="(a4,i7.7)") "real",i
            if (j == 2) write (currentName,fmt="(a4,i7.7)") "imag",i
            currentName = trim (currentName)
            call h5dcreate_f (atomOverlapCV_PSCF_gid,currentName,&
                  & H5T_NATIVE_DOUBLE,coreValePSCF_dsid,&
                  & atomOverlapCV_PSCF_did(j,i),hdferr,coreValePSCF_plid)
            if (hdferr /= 0) stop 'Failed to create atom overlapCV did'
         enddo
      endif
   enddo

   ! Repeat for the SYBD Path Kpoints.
   do i = 1, numPathKP_HDF5
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5dcreate_f (atomOverlapSYBD_PSCF_gid,currentName,&
            & H5T_NATIVE_DOUBLE,valeValePSCF_dsid,&
            & atomOverlapSYBD_PSCF_did(i),hdferr,valeValePSCF_plid)
      if (hdferr /= 0) stop 'Failed to create atom overlap SYBD did'

      do j = 1, spin
         call h5dcreate_f (atomHamOverlapSYBD_PSCF_gid,currentName,&
               & H5T_NATIVE_DOUBLE,valeValePSCF_dsid,&
               & atomHamOverlapSYBD_PSCF_did(i,j),hdferr,valeValePSCF_plid)
         if (hdferr /= 0) stop 'Failed to create hamiltonian overlap SYBD did'
      enddo

      do j = 1, 3
         call h5dcreate_f (atomDMxyzOL_SYBD_PSCF_gid(j),currentName,&
               & H5T_NATIVE_DOUBLE,valeValePSCF_dsid,&
               & atomDMOverlapSYBD_PSCF_did(i,j),hdferr,valeValePSCF_plid)
         if (hdferr /= 0) stop 'Failed to create DM overlap SYBD did'

         call h5dcreate_f (atomMMxyzOL_SYBD_PSCF_gid(j),currentName,&
               & H5T_NATIVE_DOUBLE,valeValePSCF_dsid,&
               & atomMMOverlapSYBD_PSCF_did(i,j),hdferr,valeValePSCF_plid)
         if (hdferr /= 0) stop 'Failed to create MM overlap SYBD did'
      enddo

      if (coreDim > 0) then
         do j = 1, numComponents
            if (j == 1) write (currentName,fmt="(a4,i7.7)") "real",i
            if (j == 2) write (currentName,fmt="(a4,i7.7)") "imag",i
            currentName = trim (currentName)
            call h5dcreate_f (atomOverlapCV_SYBD_PSCF_gid,currentName,&
                  & H5T_NATIVE_DOUBLE,coreValePSCF_dsid,&
                  & atomOverlapCV_SYBD_PSCF_did(j,i),hdferr,coreValePSCF_plid)
            if (hdferr /= 0) stop 'Failed to create atom overlapCV SYBD did'
         enddo
      endif
   enddo

   ! Create the attribute that indicates completion of the calculation of
   !   each integral type.
   call h5acreate_f (atomOverlapPSCF_gid,"status",H5T_NATIVE_INTEGER,&
         & attribIntPSCF_dsid,atomOverlapPSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom overlap aid'

   call h5acreate_f (atomHamOverlapPSCF_gid,"status",H5T_NATIVE_INTEGER,&
         & attribIntPSCF_dsid,atomHamOverlapPSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create hamiltonian overlap aid'

   call h5acreate_f (atomDMOverlapPSCF_gid,"status",H5T_NATIVE_INTEGER,&
         & attribIntPSCF_dsid,atomDMOverlapPSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create DM overlap aid'

   call h5acreate_f (atomMMOverlapPSCF_gid,"status",H5T_NATIVE_INTEGER,&
         & attribIntPSCF_dsid,atomMMOverlapPSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create MM overlap aid'

   call h5acreate_f (atomOverlapSYBD_PSCF_gid,"status",H5T_NATIVE_INTEGER,&
         & attribIntPSCF_dsid,atomOverlapSYBD_PSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom overlap SYBD aid'

   call h5acreate_f (atomHamOverlapSYBD_PSCF_gid,"status",H5T_NATIVE_INTEGER,&
         & attribIntPSCF_dsid,atomHamOverlapSYBD_PSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create hamiltonian overlap SYBD aid'

   call h5acreate_f (atomDMOverlapSYBD_PSCF_gid,"status",H5T_NATIVE_INTEGER,&
         & attribIntPSCF_dsid,atomDMOverlapSYBD_PSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create DM overlap SYBD aid'

   call h5acreate_f (atomMMOverlapSYBD_PSCF_gid,"status",H5T_NATIVE_INTEGER,&
         & attribIntPSCF_dsid,atomMMOverlapSYBD_PSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create MM overlap SYBD aid'

   ! Initialize all dataset attributes to the uncomputed (status = zero) state.
   call h5awrite_f(atomOverlapPSCF_aid,H5T_NATIVE_INTEGER,0,attribIntDimsPSCF,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomOverlapPSCF_aid'

   call h5awrite_f(atomHamOverlapPSCF_aid,H5T_NATIVE_INTEGER,0,&
         & attribIntDimsPSCF,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomHamOverlapPSCF_aid'

   call h5awrite_f(atomDMOverlapPSCF_aid,H5T_NATIVE_INTEGER,0,&
         & attribIntDimsPSCF,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomDMOverlapPSCF_aid'

   call h5awrite_f(atomMMOverlapPSCF_aid,H5T_NATIVE_INTEGER,0,&
         & attribIntDimsPSCF,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomMMOverlapPSCF_aid'

   call h5awrite_f(atomOverlapSYBD_PSCF_aid,H5T_NATIVE_INTEGER,0,&
         & attribIntDimsPSCF,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomOverlapSYBD_PSCF_aid'

   call h5awrite_f(atomHamOverlapSYBD_PSCF_aid,H5T_NATIVE_INTEGER,0,&
         & attribIntDimsPSCF,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomHamOverlapSYBD_PSCF_aid'

   call h5awrite_f(atomDMOverlapSYBD_PSCF_aid,H5T_NATIVE_INTEGER,0,&
         & attribIntDimsPSCF,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomDMOverlapSYBD_PSCF_aid'

   call h5awrite_f(atomMMOverlapSYBD_PSCF_aid,H5T_NATIVE_INTEGER,0,&
         & attribIntDimsPSCF,hdferr)
   if (hdferr /= 0) stop 'Failed to initialize atomMMOverlapSYBD_PSCF_aid'

   ! At this point, we flush all meta data to the PSCF HDF5 file. Then, the
   !   HDF5 file is primed for use.
   call h5fflush_f(pscf_fid,H5F_SCOPE_GLOBAL_F,hdferr)

end subroutine initPSCFIntegralHDF5


subroutine accessPSCFIntegralHDF5 (pscf_fid)

   ! Import any necessary definition modules.
   use HDF5
   use O_MPI

   ! Import necessary object modules.
   use O_Potential,   only: spin
   use O_KPoints,     only: numKPoints_HDF5, numPathKP_HDF5
   use O_AtomicSites, only: coreDim, valeDim

   ! Define the passed parameters.
   integer(hid_t) :: pscf_fid

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
   packedVVDimsPSCF(1) = numComponents
   packedVVDimsPSCF(2) = valeDim*(valeDim+1)/2 ! Linear storage of 1/2 matrix.
   if (coreDim > 0) then
      fullCVDimsPSCF(1) = coreDim
      fullCVDimsPSCF(2) = valeDim
   endif

   ! For the integrals that have multiple matrices per kpoint, we need to
   !   create space to hold the group IDs of each kpoint group.
   allocate (atomDMxyzOL_PSCF_gid (3)) ! Needs x,y,z matrices
   allocate (atomMMxyzOL_PSCF_gid (3)) ! Needs x,y,z matrices
   allocate (atomDMxyzOL_SYBD_PSCF_gid (3)) ! x,y,z matrices all done together
   allocate (atomMMxyzOL_SYBD_PSCF_gid (3)) ! x,y,z matrices all done together

   ! Allocate space to hold the dataset IDs.
   if (coreDim > 0) then
      allocate (atomOverlapCV_PSCF_did (numComponents,numKPoints_HDF5))
   else
      allocate (atomOverlapCV_PSCF_did (1,1))
   endif
   allocate (atomOverlapPSCF_did (numKPoints_HDF5))
   allocate (atomHamOverlapPSCF_did (numKPoints_HDF5,spin))
   allocate (atomDMOverlapPSCF_did (numKPoints_HDF5,3))
   allocate (atomMMOverlapPSCF_did (numKPoints_HDF5,3))

   if (coreDim > 0) then
      allocate (atomOverlapCV_SYBD_PSCF_did (numComponents,numPathKP_HDF5))
   else
      allocate (atomOverlapCV_SYBD_PSCF_did (1,1))
   endif
   allocate (atomOverlapSYBD_PSCF_did (numPathKP_HDF5))
   allocate (atomHamOverlapSYBD_PSCF_did (numPathKP_HDF5,spin))
   allocate (atomDMOverlapSYBD_PSCF_did (numPathKP_HDF5,3))
   allocate (atomMMOverlapSYBD_PSCF_did (numPathKP_HDF5,3))

   ! Only process 0 opens the HDF5 structure.
   if (mpiRank /= 0) return

   ! Open the Integral group within the pscf HDF5 file.
   call h5gopen_f (pscf_fid,"/atomIntgGroup",atomIntgGroupPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom intg group'

   ! Open the subgroups within the atomIntgGroup. We open all the subgroups
   !   even if we don't need them so that the close process is still the
   !   same as the "init" case.
   call h5gopen_f (atomIntgGroupPSCF_gid,"atomOverlap",&
         & atomOverlapPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom overlap group'
   call h5gopen_f (atomIntgGroupPSCF_gid,"atomHamOverlap",&
         & atomHamOverlapPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom hamiltonian overlap group'
   if (coreDim > 0) then
      call h5gopen_f (atomIntgGroupPSCF_gid,"atomOverlapCV",&
            & atomOverlapCV_PSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open atom overlapCV group'
   endif
   call h5gopen_f (atomIntgGroupPSCF_gid,"atomDMOverlap",&
         & atomDMOverlapPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open dipole moment overlap group'
   call h5gopen_f (atomIntgGroupPSCF_gid,"atomMMOverlap",&
         & atomMMOverlapPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentum matrix overlap group'

   call h5gopen_f (atomIntgGroupPSCF_gid,"atomOverlapSYBD",&
         & atomOverlapSYBD_PSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom overlap SYBD group'
   if (coreDim > 0) then
      call h5gopen_f (atomIntgGroupPSCF_gid,"atomOverlapCV_SYBD",&
            & atomOverlapCV_SYBD_PSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open atom overlapCV SYBD group'
   endif
   call h5gopen_f (atomIntgGroupPSCF_gid,"atomHamOverlapSYBD",&
         & atomHamOverlapSYBD_PSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom hamiltonian overlap SYBD group'
   call h5gopen_f (atomIntgGroupPSCF_gid,"atomDMOverlapSYBD",&
         & atomDMOverlapSYBD_PSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open dipole moment overlap SYBD group'
   call h5gopen_f (atomIntgGroupPSCF_gid,"atomMMOverlapSYBD",&
         & atomMMOverlapSYBD_PSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentum matrix overlap SYBD group'

   ! Now, we open the xyz subgroups for the integrals that have them.
   do i = 1, 3
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5gopen_f (atomDMOverlapPSCF_gid,currentName,&
            & atomDMxyzOL_PSCF_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open dipole moment xyz group'

      call h5gopen_f (atomMMOverlapPSCF_gid,currentName,&
            & atomMMxyzOL_PSCF_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open momentum matrix xyz group'

      call h5gopen_f (atomDMOverlapSYBD_PSCF_gid,currentName,&
            & atomDMxyzOL_SYBD_PSCF_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open dipole moment xyz SYBD group'

      call h5gopen_f (atomMMOverlapSYBD_PSCF_gid,currentName,&
            & atomMMxyzOL_SYBD_PSCF_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open momentum matrix xyz SYBD group'
   enddo

   ! Open the datasets that will be used for all subgroups of atomIntgGroup.
   do i = 1, numKPoints_HDF5
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5dopen_f (atomOverlapPSCF_gid,currentName,&
            & atomOverlapPSCF_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open atom overlap did'

      do j = 1, spin
         call h5dopen_f (atomHamOverlapPSCF_gid,currentName,&
               & atomHamOverlapPSCF_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open hamiltonian overlap did'
      enddo

      do j = 1, 3
         call h5dopen_f (atomDMxyzOL_PSCF_gid(j),currentName,&
               & atomDMOverlapPSCF_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open DM overlap did'

         call h5dopen_f (atomMMxyzOL_PSCF_gid(j),currentName,&
               & atomMMOverlapPSCF_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open MM overlap did'
      enddo

      if (coreDim > 0) then
         do j = 1, numComponents
            if (j == 1) write (currentName,fmt="(a4,i7.7)") "real",i
            if (j == 2) write (currentName,fmt="(a4,i7.7)") "imag",i
            currentName = trim (currentName)
            call h5dopen_f (atomOverlapCV_PSCF_gid,currentName,&
                  & atomOverlapCV_PSCF_did(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to create atom overlapCV did'
         enddo
      endif
   enddo

   ! Repeat for the SYBD Path KPoints
   do i = 1, numPathKP_HDF5
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5dopen_f (atomOverlapSYBD_PSCF_gid,currentName,&
            & atomOverlapSYBD_PSCF_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open atom overlap SYBD did'

      do j = 1, spin
         call h5dopen_f (atomHamOverlapSYBD_PSCF_gid,currentName,&
               & atomHamOverlapSYBD_PSCF_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open hamiltonian overlap SYBD did'
      enddo

      do j = 1, 3
         call h5dopen_f (atomDMxyzOL_SYBD_PSCF_gid(j),currentName,&
               & atomDMOverlapSYBD_PSCF_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open DM overlap SYBD did'

         call h5dopen_f (atomMMxyzOL_SYBD_PSCF_gid(j),currentName,&
               & atomMMOverlapSYBD_PSCF_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open MM overlap SYBD did'
      enddo

      if (coreDim > 0) then
         do j = 1, numComponents
            if (j == 1) write (currentName,fmt="(a4,i7.7)") "real",i
            if (j == 2) write (currentName,fmt="(a4,i7.7)") "imag",i
            currentName = trim (currentName)

            call h5dopen_f (atomOverlapCV_SYBD_PSCF_gid,currentName,&
                  & atomOverlapCV_SYBD_PSCF_did(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to create atom overlapCV SYBD did'
         enddo
      endif
   enddo

   ! Open the attributes of each integral group.
   call h5aopen_f (atomOverlapPSCF_gid,'status',atomOverlapPSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom overlap aid'

   call h5aopen_f (atomHamOverlapPSCF_gid,'status',atomHamOverlapPSCF_aid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open hamiltonian overlap aid'

   call h5aopen_f (atomDMOverlapPSCF_gid,'status',atomDMOverlapPSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open DM overlap aid'

   call h5aopen_f (atomMMOverlapPSCF_gid,'status',atomMMOverlapPSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open MM overlap aid'

   call h5aopen_f (atomOverlapSYBD_PSCF_gid,'status',atomOverlapSYBD_PSCF_aid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open atom overlap SYBD aid'

   call h5aopen_f (atomHamOverlapSYBD_PSCF_gid,'status',&
         & atomHamOverlapSYBD_PSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open hamiltonian overlap SYBD aid'

   call h5aopen_f (atomDMOverlapSYBD_PSCF_gid,'status',&
         & atomDMOverlapSYBD_PSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open DM overlap SYBD aid'

   call h5aopen_f (atomMMOverlapSYBD_PSCF_gid,'status',&
         & atomMMOverlapSYBD_PSCF_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open MM overlap SYBD aid'

   ! Obtain the properties of the datasets that were just opened. They are all
   !   the same and so only one copy is necessary.  (Actually, this value is
   !   not really used, but in the "close" subroutine we close this id so we
   !   should make sure to have it for both the setup and main calls to the
   !   close routine.) Also open one for the CV data.
   call h5dget_create_plist_f (atomOverlapPSCF_did(1),valeValePSCF_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain vale vale plid'
   if (coreDim > 0) then
      call h5dget_create_plist_f (atomOverlapCV_PSCF_did(1,1),&
            & coreValePSCF_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to obtain core vale plid'
   endif

   ! Obtain the dataspace that is used for each dataset in atomIntgGroup
   !   and all of its subgroups.  The same dataspace definition works for all
   !   of the datasets.  (Same as for the plist above.)
   call h5dget_space_f(atomOverlapPSCF_did(1),valeValePSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain vale vale dsid'
   if (coreDim > 0) then
      call h5dget_space_f(atomOverlapCV_PSCF_did(1,1),coreValePSCF_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to obtain core vale dsid'
   endif

end subroutine accessPSCFIntegralHDF5


subroutine closePSCFIntegralHDF5

   ! Import any necessary definition modules.
   use HDF5
   use O_MPI

   ! Import necessary object modules.
   use O_Potential, only: spin
   use O_KPoints, only: numKPoints_HDF5, numPathKP_HDF5
   use O_AtomicSites, only: coreDim

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: i,j
   integer :: hdferr

   if (mpiRank == 0) then

      ! Close the property list first.
      call h5pclose_f (valeValePSCF_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to close valeValePSCF_plid.'

      ! Close the datasets next.
      do i = 1, numKPoints_HDF5
         call h5dclose_f (atomOverlapPSCF_did(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomOverlapPSCF_did.'

         do j = 1, spin
            call h5dclose_f (atomHamOverlapPSCF_did(i,j),hdferr)
            if (hdferr /= 0) stop 'Failed to close atomHamOverlapPSCF_did.'
         enddo

         do j = 1, 3
            call h5dclose_f (atomDMOverlapPSCF_did(i,j),hdferr)
            if (hdferr /= 0) stop 'Failed to close atomDMOverlapPSCF_did.'
         enddo

         do j = 1, 3
            call h5dclose_f (atomMMOverlapPSCF_did(i,j),hdferr)
            if (hdferr /= 0) stop 'Failed to close atomMMOverlapPSCF_did.'
         enddo

         if (coreDim > 0) then
            do j = 1, numComponents
               call h5dclose_f (atomOverlapCV_PSCF_did(j,i),hdferr)
               if (hdferr /= 0) stop 'Failed to close atomOverlapCV_PSCF_did.'
            enddo
         endif
      enddo

      ! Repeat for the SYBD Path KPoints
      do i = 1, numPathKP_HDF5
         call h5dclose_f (atomOverlapSYBD_PSCF_did(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomOverlapSYBD_PSCF_did.'

         do j = 1, spin
            call h5dclose_f (atomHamOverlapSYBD_PSCF_did(i,j),hdferr)
            if (hdferr /= 0) stop &
                  & 'Failed to close atomHamOverlapSYBD_PSCF_did.'
         enddo

         do j = 1, 3
            call h5dclose_f (atomDMOverlapSYBD_PSCF_did(i,j),hdferr)
            if (hdferr /= 0) stop 'Failed to close atomDMOverlapSYBD_PSCF_did.'
         enddo

         do j = 1, 3
            call h5dclose_f (atomMMOverlapSYBD_PSCF_did(i,j),hdferr)
            if (hdferr /= 0) stop 'Failed to close atomMMOverlapSYBD_PSCF_did.'
         enddo

         if (coreDim > 0) then
            do j = 1, numComponents
               call h5dclose_f (atomOverlapCV_SYBD_PSCF_did(j,i),hdferr)
               if (hdferr /= 0) stop &
                     & 'Failed to close atomOverlapCV_SYBD_PSCF_did.'
            enddo
         endif
      enddo

      ! Close the data spaces next.
      call h5sclose_f (valeValePSCF_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close valeValePSCF_dsid.'
      if (coreDim > 0) then
         call h5sclose_f (coreValePSCF_dsid,hdferr)
         if (hdferr /= 0) stop 'Failed to close coreValePSCF_dsid.'
      endif

      ! Close the groups.
      do i = 1, 3
         call h5gclose_f (atomMMxyzOL_PSCF_gid(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomMMxyzOL_PSCF_gid.'
         call h5gclose_f (atomMMxyzOL_SYBD_PSCF_gid(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomMMxyzOL_SYBD_PSCF_gid.'
      enddo

      do i = 1, 3
         call h5gclose_f (atomDMxyzOL_PSCF_gid(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomDMxyzOL_PSCF_gid.'
         call h5gclose_f (atomDMxyzOL_SYBD_PSCF_gid(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomDMxyzOL_SYBD_PSCF_gid.'
      enddo

      call h5gclose_f (atomOverlapPSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomOverlapPSCF_gid.'

      if (coreDim > 0) then
         call h5gclose_f (atomOverlapCV_PSCF_gid,hdferr)
         if (hdferr /= 0) stop 'Failed to close atomOverlapCV_PSCF_gid.'
      endif

      call h5gclose_f (atomHamOverlapPSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomHamOverlapPSCF_gid.'

      call h5gclose_f (atomMMOverlapPSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomMMxyzOverlapPSCF_gid.'

      call h5gclose_f (atomDMOverlapPSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomDMxyzOverlapPSCF_gid.'

      call h5gclose_f (atomOverlapSYBD_PSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomOverlapSYBD_PSCF_gid.'

      if (coreDim > 0) then
         call h5gclose_f (atomOverlapCV_SYBD_PSCF_gid,hdferr)
         if (hdferr /= 0) stop 'Failed to close atomOverlapCV_SYBD_PSCF_gid.'
      endif

      call h5gclose_f (atomHamOverlapSYBD_PSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomHamOverlapSYBD_PSCF_gid.'

      call h5gclose_f (atomMMOverlapSYBD_PSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomMMxyzOverlapSYBD_PSCF_gid.'

      call h5gclose_f (atomDMOverlapSYBD_PSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomDMxyzOverlapSYBD_PSCF_gid.'

      call h5gclose_f (atomIntgGroupPSCF_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomIntgGroupPSCF_gid.'
   endif

   ! The attributes tracking completion are closed when written.

   ! Deallocate the arrays that hold the id numbers.
   deallocate (atomOverlapCV_PSCF_did) ! Sometimes unused, always allocated.
   deallocate (atomOverlapPSCF_did)
   deallocate (atomHamOverlapPSCF_did)
   deallocate (atomDMOverlapPSCF_did)
   deallocate (atomMMOverlapPSCF_did)
   deallocate (atomDMxyzOL_PSCF_gid)
   deallocate (atomMMxyzOL_PSCF_gid)

   deallocate (atomOverlapCV_SYBD_PSCF_did) ! Sometimes unused, always alloc.
   deallocate (atomOverlapSYBD_PSCF_did)
   deallocate (atomHamOverlapSYBD_PSCF_did)
   deallocate (atomDMOverlapSYBD_PSCF_did)
   deallocate (atomMMOverlapSYBD_PSCF_did)
   deallocate (atomDMxyzOL_SYBD_PSCF_gid)
   deallocate (atomMMxyzOL_SYBD_PSCF_gid)

   ! Note that the attributes are closed as soon as they are finished.
   
end subroutine closePSCFIntegralHDF5

end module O_PSCFIntegralsHDF5
