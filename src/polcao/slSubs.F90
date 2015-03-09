!! This module and its subroutines were  modeled after various scalapack 
!! interfaces in the Global Arrays, and NWChem packages from 
!! Pacific Northwest National Lab (PNNL)
!! Author: James E. Currie
!! Email: jecyrd@mail.umkc.edu
!!
!! See notes above the oga_zhegvx subroutine (in this file) 
!! for additional info explaining motivation and non-elegance.
module O_SlSubs
  use MPI
  implicit none

  contains



integer*4 function slgetmxproc(n,nnodes)
      implicit none
      INTGR4 nnodes
      INTGR4 n
      INTGR4 i
      double precision fact
      INTGR4 nmax,nprocs,twoi
      double precision nprocs0
      double precision otto
      parameter(nmax=11,fact=((7108d0*7108d0)/1024d0),otto=8d0)
!new      parameter(nmax=11,fact=((7108d0*7108d0)/512d0),otto=8d0)
!     lower bound of 8 procs
      nprocs0=max((n*n)/fact,otto)

!     try to get powers of two

      do i = nmax, 0, -1 
         if(nint(nprocs0/(2d0**i)).eq.1) goto 1
      enddo
      i=4
1     twoi=2**i
      slgetmxproc=min(nnodes,twoi)
      return
end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is an interface from the Global Arrays Toolkit 
! to the PZHEGVX Scalapack subroutine. Documentation for the subroutine 
! can be found at:
! http://www.netlib.org/scalapack/explore-html/d7/dff/pzhegvx_8f.html
!
! Inspiration was drawn heavily from the Global Arrays Toolkit source code
! subroutine ga_pdsygv in the scalapack.F file from version: 5.3
! Global Arrays Website: hpc.pnl.gov/globalarrays
!
! At the time no interface existed for this particular subroutine which
! the non-expert scalapack version ZHEGV was used in the current
! (at the time of this creation) serial implementation of OLCAO (rev. 29).
! Additionally it was suspected that the ga_diag subroutine might work.
! However testing showed that this subroutine while originally looking
! very simplistic was actually sort of unweildy in the sense of compiling
! OLCAO, which led me to just implement it myself.
!
! The subroutine ASSUMES that the ga_block_cyclic, or 
! ga_set_block_cyclic_proc or any other global routines have NOT been used 
! to specify the topology of the problem. This was done at the time to
! speed the implementation of the subroutine by assuming that most who will
! encounter this subroutine are ignorant to most global array features.
! (Lets face it, most who encounter OLCAO are ignorant to programming in
! general) In that way we we can prevent bugs and make my life easier
! simultaneously.
!
! From Scalapack documentation:
! Purpose
! =======
! PZHEGVX computes all the eigenvalues, and optionally, the eigenvectors
! of a complex generalized Hermitian-definite eigenproblem, of the form
! sub( A )*x=(lambda)*sub( B )*x,  
! sub( A )*sub( B )x=(lambda)*x,  or
! sub( B )*sub( A )*x=(lambda)*x.
! Here sub( A ) denoting A( IA:IA+N-1, JA:JA+N-1 ) is assumed to be
! Hermitian, and sub( B ) denoting B( IB:IB+N-1, JB:JB+N-1 ) is assumed
! to be Hermitian positive definite.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine oga_pzhegv(ga_a, ga_b, ga_s, eigenVals)
  use O_kinds

  implicit none

#include "mafdecls.fh"
#include "global.fh"

  ! Define passed parameters
  integer, intent(inout) :: ga_a  ! Matrix A
  integer, intent(inout) :: ga_b  ! Matrix B
  integer, intent(inout) :: ga_s  ! Matrix B

  ! Eigenvalues, remember eigenvalues are real
  real (kind=double), intent(out), dimension(:) :: eigenVals

  ! Define variables for use in this subroutine
  logical status

  character*1 :: jobz, uplo, range    ! Input parameters for scalapack.
                               ! We know these ahead of time.

  integer :: ha, adra        ! Addresses for A global array
  integer :: hb, adrb        ! Addresses for B global array

  logical :: oactive         ! This is true iff this process participates

  ! Vars to hold information about the locally stored GA's
  integer :: Adim1, Adim2, typeA
  integer :: Bdim1, Bdim2, typeB

  integer :: mpA, nqA        ! rows/cols of A held by the processor
  integer :: mpB, nqB        ! rows/cols of B held by the processor
  integer   :: mps, nqs        ! node ID for this process

  integer :: lda, ldb        ! leading dimensions
  integer :: elemA, elemB    ! Global array elements
  integer :: numroc

  integer :: nb              ! block size
  integer, dimension(9) :: descA, descB   ! BLACS Array descriptors

  integer :: ngaps, hgap, adrgaps
  integer :: iclu, hclustr, adrclustr
  integer :: hfail, adrfail, lf
  integer :: liwork,hiwork, adriwork
  integer :: lcwork,hcwork, adrcwork
  integer :: liwork4, lcwork4

  real (kind=double) :: vl, vu, abstol, orfac
  integer :: nn,mq0,np
  integer :: il, ui
  integer :: m, nz
  integer :: n
  integer :: info, one4, zero4, two4, four4
  integer :: info8, dblsize
  integer :: two4n

  real (kind=double) :: pdlamch, dum
  integer :: iceil, indxg2p
  logical uses_sl_A, uses_sl_B
  integer alen, blen
  integer, dimension(2) :: block_dims_A, block_dims_B, blocks
  integer, dimension(2) :: gridA(2), gridB(2)
  logical use_direct

  integer :: i,j, me

  external pdlamch, iceil, indxg2p

  parameter(zero4=0,one4=1,two4=2,four4=4)

  ! check environment
  me = ga_nodeid()

  ! check GA info for input arrays
  call ga_check_handle(ga_a, 'In subroutine oga_zhegv, & 
    & passed ga_a is not a valid Global Arrays handle.')
  call ga_check_handle(ga_b, 'In subroutine oga_zhegv, &
    & passed ga_b is not a valid Global Arrays handle.')

  call ga_inquire(ga_a, typeA, dimA18, dimA28)
  call ga_inquire(ga_b, typeB, dimB18, dimB28)
!  call ga_inquire(ga_eVects, typeB, dimB18, dimB28)

  n = dimA1

  if (nb .lt. 1) nb = 1

  if (dimA1 .ne. dimA2) then
    call ga_error('oga_pzhegv: matrix A not square ',0)
  endif
  
  if (dimB1 .ne. dimB2) then
    call ga_error('oga_pzhegv: matrix B not square ',0)
  endif

  if (dimB1 .ne. n) then
    call ga_error('oga_pzhegv: size matrix A and B differ',0)
  endif

  ! initialize Scalapack interface
  call SLinit2(n)
  

  ! call ga_sync before to enter in BLACS and after
  call ga_sync()
 
  call blacs_pinfo(iam, nnodes)

  ! determine optimal nprocs for eigensolvers based on matrix size n

  maxproc=slgetmxproc(n,nnodes)
  call FindGrid(maxproc, nprow2, npcol2)

  call blacs_get( zero4, zero4, islctxt2 )
  call blacs_gridinit(islctxt2, 'R', nprow2, npcol2)
  
  if(iam.lt.maxproc) then
     call blacs_gridinfo(iSLctxt2, nprow2, npcol2, myrow2, mycol2)
  else
     nprow2=0
     npcol2=0
     myrow2=0
     mycol2=0
  endif
  init2=.true.

  call ga_sync()

  oactive = iam .lt. maxproc
  if (oactive) then

  ! find SBS format parameters
    mpA = numroc(dimA1, nb, myrow2, zero4, nprow2)
    nqA = numroc(dimA2, nb, mycol2, zero4, npcol2)

    mpA = numroc(dimB1, nb, myrow2, zero4, nprow2)
    nqB = numroc(dimB2, nb, mycol2, zero4, npcol2)

    mpS = numroc(dimS1, nb, myrow2, zero4, nprow2)
    nqS = numroc(dimS2, nb, mycol2, zero4, npcol2)

    lda = max(one4, mpA)
    ldb = max(one4, mpB)
    ldb = max(one4, mpS)

  ! Let scalapack check for errors
    elemA = mpA * nqA
    status = .true.
    if (elemA .ne. 0) &
      & status = ma_push_get(MT_DCPL, elemA, 'a', ha, adra)
    if (.not. status) &
      & call ga_error('oga_zhegv: mem alloc failed A ', -1)

    elemB = mpB * nqB

    if (elemB .ne. 0) &
      & status = ma_push_get(MT_DCPL, elemB, 'b', hb, adrb)
    if (.not. status) &
      & call ga_error('oga_zhegv: mem alloc failed B ', -1)

  ! copy ga_a to A using block cyclic scalapack format

    call ga_to_ZSL(ga_a, dimA1, dimA2, nb, nb, dcpl_mb(adrA), &
      & lda, mpA, nqA)
    call ga_to_ZSL(ga_b, dimB1, dimB2, nb, nb, dcpl_mb(adrB), &
      & ldb, mpB, nqB)
  
    call ga_sync()

    ! Copy Global Arrays to local using block cyclick scalapack format
    call ga_to_SL2(ga_a, Adim1, Adim2, nb, nb, dbl_mb(adrA), lda, mpA, nqA)
    call ga_to_SL2(ga_b, Bdim1, Bdim2, nb, nb, dbl_mb(adrB), ldb, mpB, nqB)
    call ga_to_SL2(ga_s, Sdim1, Sdim2, nb, nb, dbl_mb(adrS), lds, mpS, nqS)


    ! fill Scalapack matrix descriptors
    call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4, &
      & islctxt2, lda, info)
    info8=info
    if(info8.ne.0) &
      & call ga_error('oga_zhegv: descinit A failed ', -info8)

    call descinit(descB, dimB1, dimB2, nb, nb, zero4, zero4, &
      & islctxt2, ldb, info)
    info8=info
    if (info8 .ne. 0) &
      & call ga_error('oga_zhegv: descinit B failed ', -info8)

    jobz = 'V'
    uplo = 'L'

    nn = max(n,nb,two4)
    np0 = numroc(nn, nb, zero4, zero4, nprow2)
    mq0 = numroc(nn, nb, zero4, zero4, npcol2)

    ! get work parameters
    jobz = 'V'

    ! make 
    !  call pzhegvx
    call PZHEGVX(IBTYPE, JOBZ, RANGE, UPLO, N A, IA, JA, &
                & DESCA, B, IB, JB, DESCB, VL, VU, IL, UI, &
                & ABSTOL, M, NZ, W, ORFAC, Z, IZ, JZ, DESCZ, &
                & WORK LWORK, RWORK, LRWORK, IWORK, LIWORK, &
                & IFAIl, ICLUSTR, GAP, INFO)


    ! Throw error flags for lapack with ga_error
    if (nz .ne. n) then
      if (info .ne. 0) then
        if (info .gt. 0) then
          call ga_error('oga_pdzhegv: argument is illegal ', info)
        else
          call ga_error('oga_pdzhegv: eigenvectors failed to converge',info)
        endif
      endif
    endif

    ! Copy solution matrix back to g_c
    call ga_from_SL2(g_b, Adim1, Bdim2, nb, nb, dbl_mb(adrB), &
      & ldb, mpb, nqb)

    ! Deallocate work/SL arrays
    if (lcwork .ne. 0) status = ma_pop_stack(hcwork)
    if (liwork .ne. 0) status = ma_pop_stack(hiwork)
    if (if .ne. 0) status = ma_pop_stack(hfail)
    if (iclu .ne. 0) status = ma_pop_stack(hclustr)
    if (ngaps .ne. 0) status = ma_pop_stack(hgap)
    if (elemS .ne. 0) status = ma_pop_stack(hs)
    if (elemB .ne. 0) status = ma_pop_stack(hb)
    if (elemA .ne. 0) status = ma_pop_stack(ha)
  endif    

  ! Sync global arrays, and broadcast eigenvalues back
  call ga_sync()
  if (maxproc .lt. nnodes .or. dima1 .le. nb) then
    ! Broadcast evals
    dblsize = ma_sizeof(MT_DBL, 1, MT_BYTE) * dima1
    call ga_brdcst(1688, eval, dblsize, 0)
  endif

  ! Sync and close the SL interface
  call ga_sync()
  if(iam.lt.maxproc)call blacs_gridexit(islctxt2)
  init2=.false.
   
end subroutine oga_pdzhegv
end module O_SlSubs
