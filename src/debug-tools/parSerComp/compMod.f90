module compMod
use dataStructs
use matSubs
use hdf5Routines

implicit none

contains
subroutine get1DDiffs(mat1, mat2, outfile)
  implicit none

  ! passed variables
  real (kind=double), intent(in), dimension(:) :: mat1, mat2
  character(len=*), intent(in) :: outfile

  ! local variables
  real (kind=double) :: diff, percent
  integer :: i

  if (size(mat1,1) == size(mat2,1)) then
    open(UNIT=101, FILE=outfile, STATUS='replace', form='formatted')
    do i=1, size(mat1,1)
      diff = abs(mat1(i) - mat2(i))
      percent = (diff/mat2(i)) * 100
      if (percent > 1.0d-7) then
        write (101,*) "i: ", i
        write (101,*) "   Mat1: ", mat1(i)
        write (101,*) "   Mat2: ", mat2(i)
        write (101,*) "   Diff: ", diff
        write (101,*) "    %:   ", percent
      endif
    enddo
    flush(101)
    close(101)
  else
    print *, "get1DDiffs recieved 2 matrices not of same size.",outfile
  endif
end subroutine get1DDiffs

subroutine get2DDiffs(mat1, mat2, outfile)
  implicit none

  ! passed variables
  real (kind=double), intent(in), dimension(:,:) :: mat1, mat2
  character(len=*), intent(in) :: outfile

  ! local variables
  real (kind=double) :: diff, percent
  integer :: i,j

  if (size(mat1,1)==size(mat2,1) .and. size(mat1,2)==size(mat2,2)) then
    open(UNIT=102, FILE=outfile, STATUS='replace', form='formatted')
    do i=1, size(mat1,1)
      do j=1, size(mat1,2)
        diff = abs(mat1(i,j) - mat2(i,j))
        percent = (diff/mat2(i,j)) * 100
        if (percent > 1.0d-7) then
          write (102,*) "i,j: ", i,j
          write (102,*) "    Mat1: ", mat1(i,j)
          write (102,*) "    Mat2: ", mat2(i,j)
          write (102,*) "    Diff: ", diff
          write (102,*) "    %:    ", percent
        endif
      enddo
    enddo
    flush(102)
    close(102)
  else
    print *, "get1DDiffs recieved 2 matrices not of same size.",outfile
  endif
end subroutine get2DDiffs

subroutine get2DDiffsCmplx(mat1, mat2, outfile)
  implicit none

  ! passed variables
  complex (kind=double), intent(in), dimension(:,:) :: mat1, mat2
  character(len=*), intent(in) :: outfile

  ! local variables
  real (kind=double) :: diff, percent
  integer :: i,j

  if (size(mat1,1)==size(mat2,1) .and. size(mat1,2)==size(mat2,2)) then
    open(UNIT=102, FILE=outfile, STATUS='replace', form='formatted')
    do i=1, size(mat1,1)
      do j=1, size(mat1,2)
        diff = abs(mat1(i,j) - mat2(i,j))
        percent = (diff/mat2(i,j)) * 100
        if (percent > 1.0d-7) then
          write (102,*) "i,j: ", i,j
          write (102,*) "    Mat1: ", mat1(i,j)
          write (102,*) "    Mat2: ", mat2(i,j)
          write (102,*) "    Diff: ", diff
          write (102,*) "    %:    ", percent
        endif
      enddo
    flush(102)
    enddo
    close(102)
  else
    print *, "get1DDiffs recieved 2 matrices not of same size.",outfile
  endif
end subroutine get2DDiffsCmplx

subroutine get3DDiffs(mat1, mat2, outfile)
  implicit none

  ! passed variables
  real (kind=double), intent(in), dimension(:,:,:) :: mat1, mat2
  character(len=*), intent(in) :: outfile

  ! local variables
  real (kind=double) :: diff,percent
  integer :: i,j,k

  if (size(mat1,1)==size(mat2,1) .and. size(mat1,2)==size(mat2,2) &
    & .and. size(mat1,3)==size(mat2,3)) then
    open(UNIT=103, FILE=outfile, STATUS='replace', form='formatted')
    do k=1, size(mat1,3)
      do i=1, size(mat1,1)
        do j=1, size(mat1,2)
          diff = abs(mat1(i,j,k) - mat2(i,j,k))
          percent = (diff/mat2(i,j,k)) * 100
          if (percent > 1.0d-7) then
            write (103,*) "i,j,k: ", i,j,k
            write (103,*) "      Mat1: ", mat1(i,j,k)
            write (103,*) "      Mat2: ", mat2(i,j,k)
            write (103,*) "      Diff: ", diff
            write (103,*) "      %:    ", percent
          endif
        enddo
      enddo
    enddo
    flush(103)
    close(103)
  else
    print *, "get1DDiffs recieved 2 matrices not of same size.",outfile
  endif
end subroutine get3DDiffs

subroutine elecStatComp()
  implicit none

  ! declare matrices for parallel components
  real (kind=double),allocatable, dimension(:,:) :: parNLResidQ
  real (kind=double),allocatable, dimension(:,:) :: parNLNeutQPot
  real (kind=double),allocatable, dimension(:)   :: parNLNucQPot
  real (kind=double),allocatable, dimension(:)   :: parlocNucQPot
  real (kind=double),allocatable, dimension(:,:) :: parlocNeutQPot
  real (kind=double),allocatable, dimension(:,:) :: parPAO

  ! declare matrices for serial components
  real (kind=double),allocatable, dimension(:,:) :: serNLResidQ
  real (kind=double),allocatable, dimension(:,:) :: serNLNeutQPot
  real (kind=double),allocatable, dimension(:)   :: serNLNucQPot
  real (kind=double),allocatable, dimension(:)   :: serlocNucQPot
  real (kind=double),allocatable, dimension(:,:) :: serlocNeutQPot
  real (kind=double),allocatable, dimension(:,:) :: serPAO

  logical, dimension(6) :: checks

  integer :: hdferr

  ! read in Parallel stuff
  call accessSetupHDF5("setup-par.hdf5")

  allocate (parNLResidQ(numPotTypes,potDim))
  allocate (parNLNeutQPot(potDim,potDim))
  allocate (parNLNucQPot(potDim))
  allocate (parlocNucQPot(potDim))
  allocate (parlocNeutQPot(potDim,potDim))
  allocate (parPAO(potDim,potDim))

  call h5dread_f(nonLocalResidualQ_did, H5T_NATIVE_DOUBLE, &
    & parNLResidQ(:,:), potTypesPot, hdferr)
  if (hdferr /= 0) stop 'Failed to read parallel nlresidq'

  call h5dread_f(nonLocalNeutQPot_did, H5T_NATIVE_DOUBLE, &
    & parNLNeutQPot, potDims2, hdferr)
  if (hdferr /= 0) stop 'Failed to read parallel nlneutqpot'
  
  call h5dread_f(nonLocalNucQPot_did, H5T_NATIVE_DOUBLE, &
    & parNLNucQPot, potDims1, hdferr)
  if (hdferr /= 0) stop 'Failed to read parallel nlnucqpot'
  
  call h5dread_f(localNucQPot_did, H5T_NATIVE_DOUBLE, &
    & parlocNucQPot, potDims1, hdferr)
  if (hdferr /= 0) stop 'Failed to read parallel locnucqpot'
  
  call h5dread_f(localNeutQPot_did, H5T_NATIVE_DOUBLE, &
    & parlocNeutQPot, potDims2, hdferr)
  if (hdferr /= 0) stop 'Failed to read parallel locneutqpot'

  call h5dread_f(potAlphaOverlap_did,H5T_NATIVE_DOUBLE, &
    & parPAO, potDims2, hdferr)
  if (hdferr /= 0) stop 'Failed to read parallel pao'

  
  call closeSetupHDF5()


  ! Read in serial stuff
  call accessSetupHDF5("setup-ser.hdf5")

  allocate (serNLResidQ(numPotTypes,potDim))
  allocate (serNLNeutQPot(potDim,potDim))
  allocate (serNLNucQPot(potDim))
  allocate (serlocNucQPot(potDim))
  allocate (serlocNeutQPot(potDim,potDim))
  allocate (serPAO(potDim,potDim))

  call h5dread_f(nonLocalResidualQ_did, H5T_NATIVE_DOUBLE, &
    & serNLResidQ, potTypesPot, hdferr)
  if (hdferr /= 0) stop 'Failed to read serial nlresidq'

  call h5dread_f(nonLocalNeutQPot_did, H5T_NATIVE_DOUBLE, &
    & serNLNeutQPot, potDims2, hdferr)
  if (hdferr /= 0) stop 'Failed to read serial nlneutqpot'
  
  call h5dread_f(nonLocalNucQPot_did, H5T_NATIVE_DOUBLE, &
    & serNLNucQPot, potDims1, hdferr)
  if (hdferr /= 0) stop 'Failed to read serial nlnucqpot'
  
  call h5dread_f(localNucQPot_did, H5T_NATIVE_DOUBLE, &
    & serlocNucQPot, potDims1, hdferr)
  if (hdferr /= 0) stop 'Failed to read serial locnucqpot'
  
  call h5dread_f(localNeutQPot_did, H5T_NATIVE_DOUBLE, &
    & serlocNeutQPot, potDims2, hdferr)
  if (hdferr /= 0) stop 'Failed to read serial locneutqpot'

  call h5dread_f(potAlphaOverlap_did,H5T_NATIVE_DOUBLE, &
    & serPAO, potDims2, hdferr)
  if (hdferr /= 0) stop 'Failed to read serial pao'

  call closeSetupHDF5()

  ! Compare the suckers
  checks(1) = ALL(parNLresidQ == serNLResidQ)
  checks(2) = ALL(parNLNeutQPot == serNLNeutQPot)
  checks(3) = ALL(parNLNucQPot == serNLNucQPot)
  checks(4) = ALL(parlocNucQPot == serlocNucQPot)
  checks(5) = ALL(parlocNeutQPot == serlocNeutQPot)
  checks(6) = ALL(parPAO == serPAO)

  print *, "nonLocalResidualQ: ", checks(1) 
  print *, "nonLocalNeutQPot:  ", checks(2) 
  print *, "nonLocalNucQPot:   ", checks(3) 
  print *, "localNucQPot:      ", checks(4) 
  print *, "localNeutQPot:     ", checks(5) 
  print *, "potAlphaOverlap:   ", checks(6) 

  if (.not. checks(1)) then
    call get2DDiffs(parNLresidQ,serNLResidQ,'nonLocalResidualQ')
  endif
  if (.not. checks(2)) then
    call get2DDiffs(parNLNeutQPot,serNLNeutQPot,'nonLocalNeutQPot')
  endif
  if (.not. checks(3)) then
    call get1DDiffs(parNLNucQPot,serNLNucQPot,'nonLocalNucQPot')
  endif
  if (.not. checks(4)) then
    call get1DDiffs(parlocNucQPot,serlocNucQPot,'localNucQPot')
  endif
  if (.not. checks(5)) then
    call get2DDiffs(parlocNeutQPot,serlocNeutQPot,'localNeutQPot')
  endif
  if (.not. checks(6)) then
    call get2DDiffs(parPAO,serPAO,'potAlphaOverlap')
  endif

  ! Deallocate parallel
  deallocate (parNLResidQ)
  deallocate (parNLNeutQPot)
  deallocate (parNLNucQPot)
  deallocate (parlocNucQPot)
  deallocate (parlocNeutQPot)
  deallocate (parPAO)
  ! Deallocate Serial
  deallocate (serNLResidQ)
  deallocate (serNLNeutQPot)
  deallocate (serNLNucQPot)
  deallocate (serlocNucQPot)
  deallocate (serlocNeutQPot)
  deallocate (serPAO)
end subroutine elecStatComp

subroutine exchCorrComp()
  implicit none

  ! Parallel
  real (kind=double), allocatable, dimension(:) :: parRW
  real (kind=double), allocatable, dimension(:,:) :: parECO
  real (kind=double), allocatable, dimension(:,:) :: parERO

  ! Serial
  real (kind=double), allocatable, dimension(:) :: serRW
  real (kind=double), allocatable, dimension(:,:) :: serECO
  real (kind=double), allocatable, dimension(:,:) :: serERO

  integer :: i
  integer :: hdferr

  logical, dimension(3) :: checks

  character*30 :: rwfile

  allocate (parRW (maxNumRayPoints)) 
  allocate (parERO(potDim,potDim))

  allocate (serRW (maxNumRayPoints)) 
  allocate (serERO(potDim,potDim))

  do i=1,numPotSites
    ! Read in parallellllll components
    call accessSetupHDF5("setup-par.hdf5")
  
    call h5dread_f (numPoints_did(i),H5T_NATIVE_INTEGER,&
          & numRayPoints,numPoints,hdferr)
    if (hdferr /= 0) stop 'Failed to read num ray points'

    call h5dread_f (radialWeight_did(i),H5T_NATIVE_DOUBLE,&
          & parRW,points,hdferr)
    if (hdferr /= 0) stop 'Failed to read radial weights'

    call h5dread_f (exchRhoOp_did(i),H5T_NATIVE_DOUBLE,&
          & parERO(:,:),potPoints,hdferr)
    if (hdferr /= 0) stop 'Failed to read exch rho operator'
    call closeSetupHDF5()
  
    ! Read in Serial Components
    call accessSetupHDF5("setup-ser.hdf5")
    call h5dread_f (radialWeight_did(i),H5T_NATIVE_DOUBLE,&
          & serRW,points,hdferr)
    if (hdferr /= 0) stop 'Failed to read radial weights'

    call h5dread_f (exchRhoOp_did(i),H5T_NATIVE_DOUBLE,&
          & serERO(:,:),potPoints,hdferr)
    if (hdferr /= 0) stop 'Failed to read exch rho operator'
    call closeSetupHDF5()
  
    ! Comp that shit
    checks(1) = ALL(parRW == serRW)
    checks(2) = ALL(parERO == serERO)

    print *, "radialWeight: ", checks(1) 
    print *, "exchRhoOp:  ", checks(2) 


    if (.not. checks(1)) then
      write (rwfile,*) i,"RW"
      rwfile = trim(rwfile)
      call get1DDiffs(parRW,serRW,rwfile)
    endif
    if (.not. checks(2)) then
      write (rwfile,*) i,"ERO"
      rwfile = trim(rwfile)
      call get2DDiffs(parERO,serERO,rwfile)
    endif

  enddo  
  deallocate (parRW)
  deallocate (parERO)
  deallocate (serRW)
  deallocate (serERO)

  ! Do the exchCorrOverlap Diff
  allocate (parECO(potDim,maxNumRayPoints))
  allocate (serECO(potDim,maxNumRayPoints))
  call accessSetupHDF5("setup-par.hdf5")
  call h5dread_f (exchCorrOverlap_did, H5T_NATIVE_DOUBLE, &
        & parECO, potDims2, hdferr)
  if (hdferr /= 0) stop 'Failed to read exch corr overlap'
  call closeSetupHDF5()

  call accessSetupHDF5("setup-ser.hdf5")
  call h5dread_f (exchCorrOverlap_did, H5T_NATIVE_DOUBLE, &
        & parECO, potDims2, hdferr)
  if (hdferr /= 0) stop 'Failed to read exch corr overlap'
  call closeSetupHDF5()
  
  checks(3) = ALL(parECO == serECO)

  print *, "exchCorrOverlap: ", checks(3)

  if (.not. checks(3)) then
    call get2DDiffs(parECO,serECO,'exchCorrOverlap')
  endif

  deallocate (parECO)
  deallocate (serECO)
end subroutine exchCorrComp

subroutine integComp()
  implicit none

  call overlapComp()
  call nucComp()
  call keComp()
  call elecComp()
end subroutine integComp

subroutine overlapComp()
  implicit none

  complex(kind=double), allocatable, dimension(:,:) :: parOVLP
  real(kind=double), allocatable, dimension(:,:) :: tempOVLP
  complex(kind=double), allocatable, dimension(:,:) :: serOVLP
  logical :: check
  character*30 :: currentname
  integer (hsize_t), dimension(2) :: matdims

  integer :: i
  integer :: valeDim
  valeDim = systemVars%valeDim
 
  allocate (parOVLP(valeDim,valeDim))
  allocate (serOVLP(valeDim,valeDim))
  do i=1,systemVars%numKPoints
    ! Parallel
    call accessSetupHDF5("setup-par.hdf5")

    allocate (tempOVLP(valeDim,valeDim))
    matdims = (/valeDim,valeDim/)
    call readPackedParMatrix(atomOverlap_did(i),tempOVLP, &
      & matdims, valeDim, valeDim)
    call unpackParMatrix(parOVLP, tempOVLP, valeDim, 0)

    deallocate(tempOVLP)
    call closeSetupHDF5()

    ! Serial
    allocate (tempOVLP(2,valeDim*(valeDim+1)/2))
    matdims = (/2,size(tempOVLP,2)/)
    call accessSetupHDF5("setup-ser.hdf5")
    call readPackedSerMatrix(atomOverlap_did(i),tempOVLP, &
      & matdims, 2, size(tempOVLP,2))

    call unpackSerMatrix(serOVLP, tempOVLP, valeDim, 0)

    call closeSetupHDF5()
    deallocate(tempOVLP)

    ! Compare
    check = ALL(parOVLP == serOVLP)

    print *, "Overlap, KPOINT:",i, check

    if (.not. check) then
      write (currentname,*) i,"kp-Overlap"
      currentname = trim(currentname)
      call get2DDiffsCmplx(parOVLP,serOVLP,currentname)
    endif
  enddo

  deallocate (parOVLP)
  deallocate (serOVLP)
end subroutine

subroutine nucComp()
  implicit none

  complex(kind=double), allocatable, dimension(:,:) :: parOVLP
  real(kind=double), allocatable, dimension(:,:) :: tempOVLP
  complex(kind=double), allocatable, dimension(:,:) :: serOVLP
  logical :: check
  character*30 :: currentname
  integer (hsize_t), dimension(2) :: matdims

  integer :: i
  integer :: valeDim
  valeDim = systemVars%valeDim
 
  allocate (parOVLP(valeDim,valeDim))
  allocate (serOVLP(valeDim,valeDim))
  do i=1,systemVars%numKPoints
    ! Parallel
    call accessSetupHDF5("setup-par.hdf5")

    allocate (tempOVLP(valeDim,valeDim))
    matdims = (/valeDim,valeDim/)
    call readPackedParMatrix(atomNucOverlap_did(i),tempOVLP, &
      & matdims, valeDim, valeDim)
    call unpackParMatrix(parOVLP, tempOVLP, valeDim, 0)

    deallocate(tempOVLP)
    call closeSetupHDF5()

    ! Serial
    allocate (tempOVLP(2,valeDim*(valeDim+1)/2))
    matdims = (/2,size(tempOVLP,2)/)
    call accessSetupHDF5("setup-ser.hdf5")
    call readPackedSerMatrix(atomNucOverlap_did(i),tempOVLP, &
      & matdims, 2, size(tempOVLP,2))

    call unpackSerMatrix(serOVLP, tempOVLP, valeDim, 0)

    call closeSetupHDF5()
    deallocate(tempOVLP)

    ! Compare
    check = ALL(parOVLP == serOVLP)

    print *, "NucOVLP, KPOINT:",i, check

    if (.not. check) then
      write (currentname,*) i,"kp-NucOvlp"
      currentname = trim(currentname)
      call get2DDiffsCmplx(parOVLP,serOVLP,currentname)
    endif
  enddo

  deallocate (parOVLP)
  deallocate (serOVLP)

end subroutine

subroutine keComp()
  implicit none

  complex(kind=double), allocatable, dimension(:,:) :: parOVLP
  real(kind=double), allocatable, dimension(:,:) :: tempOVLP
  complex(kind=double), allocatable, dimension(:,:) :: serOVLP
  logical :: check
  character*30 :: currentname
  integer (hsize_t), dimension(2) :: matdims

  integer :: i
  integer :: valeDim
  valeDim = systemVars%valeDim
 
  allocate (parOVLP(valeDim,valeDim))
  allocate (serOVLP(valeDim,valeDim))
  do i=1,systemVars%numKPoints
    ! Parallel
    call accessSetupHDF5("setup-par.hdf5")

    allocate (tempOVLP(valeDim,valeDim))
    matdims = (/valeDim,valeDim/)
    call readPackedParMatrix(atomKEOverlap_did(i),tempOVLP, &
      & matdims, valeDim, valeDim)
    call unpackParMatrix(parOVLP, tempOVLP, valeDim, 0)

    deallocate(tempOVLP)
    call closeSetupHDF5()

    ! Serial
    allocate (tempOVLP(2,valeDim*(valeDim+1)/2))
    matdims = (/2,size(tempOVLP,2)/)
    call accessSetupHDF5("setup-ser.hdf5")
    call readPackedSerMatrix(atomKEOverlap_did(i),tempOVLP, &
      & matdims, 2, size(tempOVLP,2))

    call unpackSerMatrix(serOVLP, tempOVLP, valeDim, 0)

    call closeSetupHDF5()
    deallocate(tempOVLP)

    ! Compare
    check = ALL(parOVLP == serOVLP)

    print *, "KEovlp, KPOINT:",i, check

    if (.not. check) then
      write (currentname,*) i,"kp-KEovlp"
      currentname = trim(currentname)
      call get2DDiffsCmplx(parOVLP,serOVLP,currentname)
    endif
  enddo

  deallocate (parOVLP)
  deallocate (serOVLP)

end subroutine

subroutine elecComp()
  implicit none

  complex(kind=double), allocatable, dimension(:,:) :: parOVLP
  real(kind=double), allocatable, dimension(:,:) :: tempOVLP
  complex(kind=double), allocatable, dimension(:,:) :: serOVLP
  logical :: check
  character(len=128) :: currentname
  character(len=128) :: temp
  integer (hsize_t), dimension(2) :: matdims


  integer :: i, j, k
  integer :: cntr=0
  integer :: valeDim
  valeDim = systemVars%valeDim
 
  allocate (parOVLP(valeDim,valeDim))
  allocate (serOVLP(valeDim,valeDim))
  do i=1,systemVars%numKPoints
    do j = 1, systemVars%potDim
      ! Parallel
      call accessSetupHDF5("setup-par.hdf5")

      allocate (tempOVLP(valeDim,valeDim))
      matdims = (/valeDim,valeDim/)
      call readPackedParMatrix(atomPotOverlap_did(i,j),tempOVLP, &
        & matdims, valeDim, valeDim)
      call unpackParMatrix(parOVLP, tempOVLP, valeDim, 0)

      deallocate(tempOVLP)
      call closeSetupHDF5()

      ! Serial
      allocate (tempOVLP(2,valeDim*(valeDim+1)/2))
      matdims = (/2,size(tempOVLP,2)/)
      call accessSetupHDF5("setup-ser.hdf5")
      call readPackedSerMatrix(atomPotOverlap_did(i,j),tempOVLP, &
        & matdims, 2, size(tempOVLP,2))

      call unpackSerMatrix(serOVLP, tempOVLP, valeDim, 0)

      call closeSetupHDF5()
      deallocate(tempOVLP)

      ! Compare
      check = ALL(parOVLP == serOVLP)

      print *, "elecOvlp, KPOINT,potDim:",i,j, check
      if (.not. check) then
        write (currentname,*) i,"kp-"
        write (temp,*) j
        currentname = trim(adjustl(currentname)) // trim(adjustl(temp)) &
          & //"dim-ElecOvlp"
        currentname = trim(adjustl(currentname))
        call get2DDiffsCmplx(parOVLP,serOVLP,currentname)
      endif
    enddo
  enddo

  deallocate (parOVLP)
  deallocate (serOVLP)

end subroutine

end module compMod
