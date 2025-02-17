module O_MPI

   use MPI_F08

   ! Make sure that no funny variables are implicitly declared.
   implicit none

   ! Define access.
   public

   ! Define the MPI variables.
   integer :: mpiRank
   integer :: mpiSize
   integer :: mpierr
   integer :: mpiStatus(MPI_STATUS_SIZE)

   ! Consider the following global 103 x 103 sized matrix. Each "element"
   !   shown represents a 10x10 sub-matrix except for the last column (which
   !   contains a 10x3 sub-matrix (10 rows, 3 cols), and the bottom row which
   !   contains a 3x10 sub-matrix (3 rows, 10 cols). A further exception is
   !   the bottom right corner sub-matrix, which is 3x3 in size.
   !
   !    ---------------------------------------------
   !    |   |   |   |   |   |   |   |   |   |   |   |
   !    ---------------------------------------------
   !    |   |   |   |   |   |   |   |   |   |   |   |
   !    ---------------------------------------------
   !    |   |   |   |   |   |   |   |   |   |   |   |
   !    ---------------------------------------------
   !    |   |   |   |   |   |   |   |   |   |   |   |
   !    ---------------------------------------------
   !    |   |   |   |   |   |   |   |   |   |   |   |
   !    ---------------------------------------------
   !    |   |   |   |   |   |   |   |   |   |   |   |
   !    ---------------------------------------------
   !    |   |   |   |   |   |   |   |   |   |   |   |
   !    ---------------------------------------------
   !    |   |   |   |   |   |   |   |   |   |   |   |
   !    ---------------------------------------------
   !    |   |   |   |   |   |   |   |   |   |   |   |
   !    ---------------------------------------------
   !    |   |   |   |   |   |   |   |   |   |   |   |
   !    ---------------------------------------------
   !    |   |   |   |   |   |   |   |   |   |   |   |
   !    ---------------------------------------------

   ! We will operate on the matrix with a set of processes organized in a block
   !   cyclic scheme because that will ensure maximum utilization of all
   !   processes with minimal communication. So, the global matrix will be
   !   divided up by, say, six processes in a 2-row x 3-col process grid as
   !   follows (where the numbers correspond to process numbers):
   !
   !    ---------------------------------------------
   !    |1  |2  |3  |1  |2  |3  |1  |2  |3  |1  |2  |
   !    ---------------------------------------------
   !    |4  |5  |6  |4  |5  |6  |4  |5  |6  |4  |5  |
   !    ---------------------------------------------
   !    |1  |2  |3  |1  |2  |3  |1  |2  |3  |1  |2  |
   !    ---------------------------------------------
   !    |4  |5  |6  |4  |5  |6  |4  |5  |6  |4  |5  |
   !    ---------------------------------------------
   !    |1  |2  |3  |1  |2  |3  |1  |2  |3  |1  |2  |
   !    ---------------------------------------------
   !    |4  |5  |6  |4  |5  |6  |4  |5  |6  |4  |5  |
   !    ---------------------------------------------
   !    |1  |2  |3  |1  |2  |3  |1  |2  |3  |1  |2  |
   !    ---------------------------------------------
   !    |4  |5  |6  |4  |5  |6  |4  |5  |6  |4  |5  |
   !    ---------------------------------------------
   !    |1  |2  |3  |1  |2  |3  |1  |2  |3  |1  |2  |
   !    ---------------------------------------------
   !    |4  |5  |6  |4  |5  |6  |4  |5  |6  |4  |5  |
   !    ---------------------------------------------
   !    |1  |2  |3  |1  |2  |3  |1  |2  |3  |1  |2  |
   !    ---------------------------------------------
   !
   ! You will note, obviously, that we did not assign process 1 to cover all
   !   elements in the upper left and process 6 to cover all elements in the
   !   lower right (and etc. for the others) because that naive coverage
   !   scheme would yield poor efficiency for parallel ScaLAPACK and PBLAS
   !   operations. Using that naive scheme, as a Gaussian elimination (for
   !   example) progressed, only the process in the lower right would be
   !   "working" while the other processes would site idle.

   ! Each process will have its own single (locally allocated) matrix that
   !   will consist of elements from different blocks the global matrix as
   !   indicated above. Because things may not "fit" perfectly, it is likely
   !   that different processes will have different numbers of blocks (shown)
   !   and even that some blocks may be of different sizes (not shown).
   !
   !    ---------------------------------------------
   !    |1aa|2ab|3ac|1ad|2ae|3af|1ag|2ah|3ai|1aj|2ak|
   !    ---------------------------------------------
   !    |4ba|5bb|6bc|4bd|5be|6bf|4bg|5bh|6bi|4bj|5bk|
   !    ---------------------------------------------
   !    |1ca|2cb|3cc|1cd|2ce|3cf|1cg|2ch|3ci|1cj|2ck|
   !    ---------------------------------------------
   !    |4da|5db|6dc|4dd|5de|6df|4dg|5dh|6di|4dj|5dk|
   !    ---------------------------------------------
   !    |1ea|2eb|3ec|1ed|2ee|3ef|1eg|2eh|3ei|1ej|2ek|
   !    ---------------------------------------------
   !    |4fa|5fb|6fc|4fd|5fe|6ff|4fg|5fh|6fi|4fj|5fk|--------
   !    ---------------------------------------------       |
   !    |1gg|2gb|3gc|1gd|2ge|3gf|1gg|2gh|3gi|1gj|2gk|       |
   !    ---------------------------------------------       |
   !    |4ha|5hb|6hc|4hd|5he|6hf|4hg|5hh|6hi|4hj|5hk|       |
   !    ---------------------------------------------       |
   !    |1ia|2ib|3ic|1id|2ie|3if|1ig|2ih|3ii|1ij|2ik|       |
   !    ---------------------------------------------       |
   !    |4ja|5jb|6jc|4jd|5je|6jf|4jg|5jh|6ji|4jj|5jk|       |
   !    ---------------------------------------------       |
   !    |1ka|2kb|3kc|1kd|2ke|3kf|1kg|2kh|3ki|1kj|2kk|       |
   !    ---------------------------------------------       |
   !                                                        |
   !                        --------------------------------|
   !                        |
   !                        v (One locally allocated matrix / process.)
   !
   !          1                  2                3
   !  -----------------  -----------------  -------------
   !  |1aa|1ad|1ag|1aj|  |2ab|2ae|2ah|2ak|  |3ac|3af|3ai|
   !  -----------------  -----------------  -------------
   !  |1ca|1cd|1cg|1cj|  |2cb|2ce|2ch|2ck|  |3cc|3cf|3ci|
   !  -----------------  -----------------  -------------
   !  |1ea|1ed|1eg|1ej|  |2eb|2ee|2eh|2ek|  |3ec|3ef|3ei|
   !  -----------------  -----------------  -------------
   !  |1ga|1gd|1gg|1gj|  |2gb|2ge|2gh|2gk|  |3gc|3gf|3gi|
   !  -----------------  -----------------  -------------
   !  |1ia|1id|1ig|1ij|  |2ib|2ie|2ih|2ik|  |3ic|3if|3ii|
   !  -----------------  -----------------  -------------
   !  |1ka|1kd|1kg|1kj|  |2kb|2ke|2kh|2kk|  |3kc|3kf|3ki|
   !  -----------------  -----------------  -------------
   !
   !          4                  5                6
   !  -----------------  -----------------  -------------
   !  |4ba|4bd|4bg|4bj|  |5bb|5be|5bh|5bk|  |6bc|6bf|6bi|
   !  -----------------  -----------------  -------------
   !  |4da|4dd|4dg|4dj|  |5db|5de|5dh|5dk|  |6dc|6df|6di|
   !  -----------------  -----------------  -------------
   !  |4fa|4fd|4fg|4fj|  |5fb|5fe|5fh|5fk|  |6fc|6ff|6fi|
   !  -----------------  -----------------  -------------
   !  |4ha|4hd|4hg|4hj|  |5hb|5he|5hh|5hk|  |6hc|6hf|6hi|
   !  -----------------  -----------------  -------------
   !  |4ja|4jd|4jg|4jj|  |5jb|5je|5jh|5jk|  |6jc|6jf|6ji|
   !  -----------------  -----------------  -------------
   !
   ! All sub-matrices with 'k' as the second index come from the last column
   !   of the original matrix. When that matrix was divvied up to the process
   !   grid, the last colum consisted of 10x3 submatrices.
   ! The last column (3ai, 3ci, ..., 6ji) consist of 10x3 sub-matrices and the
   !   last row (4ja, 4jd, ..., 6ji) consist of 3x10 sub-matrics. (Except, of
   !   course, for element 6ji itself which is a 3x3 sub-matrix.)
   !
   ! So, each process needs to be able to allocate the space it needs for
   !   its elements and it needs to be able to make a back-and-forth mapping
   !   between the indices of the specific matrix elements within the
   !   sub-matrices it owns and the indices of the specific matrix elements
   !   in the global matrix.
   ! For example, consider sub-matrix 2ge. Since each sub-matrix of process 2
   !   from 2ab down to 2ge consists of 10x10 elements we can say that the
   !   element in the upper left corner of 2ge (with indices starting at 1,1
   !   in the upper left corner of sub-matrix 2ab) must be indexed in the
   !   local array of process 2 as 31,11 (i.e., row 31, column 11). The
   !   local matrix element in the lower right corner of sub-matrix 2ge is
   !   at local indices of 40,20. Now, the element at 31,11 in the local
   !   matrix of process 2 is at the _global_ index of 61,41 (i.e., row 61,
   !   column 41). Similarly, the element at 40,20 in the local matrix is
   !   at index 70,50 of the global matrix.

   ! The parallel process grid.
   integer :: pGridRows  ! In the example above, this is 2.
   integer :: pGridCols  ! In the example above, this is 3.

   ! The descriptors of distributed matrices. Note that a process "covers"
   !   elements in the global matrix following a block cyclic scheme. So,
   !   while a process may be identified with a specific row/col in the
   !   process grid, the same process will be identified with a disjoint
   !   set of matrix elements from the global matrix (again, following a
   !   block cyclic scheme).
   ! Everything is easy when everything fits exactly (i.e., the global matrix
   !   size and the number of processes in each dimension of the process grid
   !   divide perfectly). However, that is not usually the case and so we
   !   must include notes about extra rows or cols.
   type MatrixDescriptor
      integer :: globalRows  ! Matrix data dimension, not grid dimension.
      integer :: globalCols  ! Matrix data dimension, not grid dimension.

      integer :: myPGridRow  ! Position in the process grid.
      integer :: myPGridCol  ! Position in the process grid.

      integer :: numRowBlocks  ! mb; Number of blocks per global matrix row.
      integer :: numColBlocks  ! nb; Number of blocks per global matrix col.

      integer :: extraRows
      integer :: extraCols

      integer :: localRows  ! Local data dimension.
      integer :: localCols  ! Local data dimension.


   end type MatrixDescriptor

   contains

subroutine initMPI

   call MPI_INIT (mpierr)
   call MPI_COMM_RANK (MPI_COMM_WORLD,mpiRank,mpierr)
   call MPI_COMM_SIZE (MPI_COMM_WORLD,mpiSize,mpierr)
   call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

end subroutine initMPI


subroutine stopMPI(errorMessage)
   
   use MPI_F08

   implicit none

   ! Declare passed parameters
   character(len=*) :: errorMessage

   if (mpiRank == 0) then
      write (20,*) errorMessage
   endif
   
   call MPI_FINALIZE(mpierr)

   stop

end subroutine stopMPI


subroutine closeMPI

   use MPI_F08

   call MPI_FINALIZE (mpierr)

end subroutine closeMPI


! As discussed below, the use of PBLAS and ScaLAPACK strongly require the use
!   of a block cyclic distribution. Further, the work will be optimally load
!   balanced if the distribution blocks are square and if every process has
!   the same number of blocks. These conditions are dependent on the number
!   of MPI processes. For example, if there are 25 processes, then we can
!   form a 5x5 grid. However, if there are 26 processes, then we can only
!   form a rectangle with a 2x13 or 13x2 grid which will probably perform
!   sub-optimally. Worse, if there are 27 processes (a prime number) then
!   we can only have a 1x27 or 27x1 process grid which is not good either.
! A second consideration is the ratio between the number of processes and
!   the size of the matrices that need to be computed. Two issues to
!   consider are the memory requirements and the on-disk chunk size. The
!   memory requirement is the total memory that a single process will
!   need to allocate times the number of processes on a single compute
!   node where the physical memory resides. For example, if a node has
!   100 GB available but the matrix that needs to be worked on requires
!   150 GB, then the job will need to be spread across two nodes and the
!   number of processes per node should be close to even. Otherwise, it may
!   be possible to request say 125 GB on one node and 25 GB on the other.
!   While the total is less than the 200 GB available, it will still cause
!   one node to be oversubscribed (and thus cause a fatal error). In a
!   similar vein, it is important to make the on-disk chunk size as large as
!   possible so that the maximum compression is achieved and the greatest
!   write (and, later, read) efficiency is achieved.
! It is incumbent on the user to request a number of processes that will
!   yield a high quality (square) process grid that maximizes the chunk
!   size while also using node memory efficiently.
subroutine createProcessGrid

   implicit none

   ! Define the process grid based on mpiSize.
   call findMostSquareProcessGrid(mpiSize)

   ! Create the descriptors that will be used for the various types of
   !   interaction integral matrices.



end subroutine createProcessGrid


! Note, this just finds the most square process grid. It does not define the
!   blocks because that requires knowledge about the sizes of the matrices
!   that need to be created (coreDim and valeDim).
subroutine findMostSquareProcessGrid(numProcs)

   use O_Kinds

   implicit none

   ! Define passed parameters.
   integer, intent(in) :: numProcs

   ! Local variables.
   integer :: maxInt, productInt
   real(kind=double) :: numProcsSqrt

   ! Compute the closest integer to the square root of numProcs. It is possible
   !   that the square root of a perfect square yields a number that cannot be
   !   exactly represented in binary and thus has finite precision rounding
   !   error. Hence, we must check if the square of the floor and ceiling
   !   return the original numProcs. If one does, then we have corrected for
   !   the error for the perfect square case.
   ! If numProcs is _not_ a perfect square, then (else) we need to find an
   !   integer divisor of numProcs and partner it with the product that
   !   yields numProcs.
   ! The math we want to arrive at here is (maxInt) x (productInt) = numProcs.

   ! Get the square root of the number of processes.
   numProcsSqrt = sqrt(real(numProcs,double))

   ! The (if) and (elseif) look for perfect squares. The (else) finds the
   !   next best thing. The (if) checks for finite precision errors that
   !   rounded up, the (elseif) checks for finite precision errors that
   !   rounded down. If the square root happened to be exact, then either
   !   will get the right answer.
   if (floor(numProcsSqrt) ** 2 == numProcs) then
      maxInt = floor(numProcsSqrt)
      productInt = maxInt
   elseif (ceiling(numProcsSqrt) ** 2 == numProcs) then
      maxInt = ceiling(numProcsSqrt)
      productInt = maxInt
   else
      ! Initialize one term in the product to the maximum possible integer
      !   that could form a product that equals numProcs.
      maxInt = floor(numProcsSqrt)

      ! Then, decrement the term until we find an integer that evenly
      !   divides numProcs. At worst, maxInt will equal one.
      do while (mod(numProcs,maxInt) /= 0)
         maxInt = maxInt - 1
      enddo

      ! At this point, we have found a perfect divisor for numProcs that is
      !   as close to square as possible.
      productInt = numProcs / maxInt
   endif

   ! The product can now be used to define the process grid.
   pGridRows = maxInt
   pGridCols = productInt

end subroutine findMostSquareProcessGrid


! Collective matrix creation and manipulation operations require use of a
!   good data distribution pattern for efficient parallel execution.
!   Specifically, all the integral subroutines make large matrices and then
!   perform matrix transpositions and multiplications (PBLAS) on them to
!   orthogonalize the valence orbital portions of the basis set against the
!   core orbitals. Each matrix is the written to disk using parallel HDF5
!   writes. Those same matrices are later read back in from disk to create
!   the Hamiltonian matrix and solve the Schroedinger equation eigen value
!   problem (presently using ScaLAPACK).
! The optimal data distribution for the PBLAS and ScaLAPACK routines is
!   block cyclic. Within that constraint, the HDF5 file will require a chunk
!   size that is exactly equal to the block size so that each write activity
!   from each process will only touch one chunk. (Some testing will need to
!   be done to determine exactly how collective write commands function with
!   the HDF5 compressed chunks. Specifically, if each process writes a
!   complicated hyperslab will it work well even if the chunk size and block
!   size match? Or perhaps it would be better to have a larger chunk size
!   and have a series of collective writes called where each process
!   contributes to a portion of the chunk with one block of their block
!   cyclic portion of the whole matrix? Or something else?)
! In any case, the creation of the matrix is also a bit tricky because it
!   is assembled from a "quilt" of non-uniform rectangular sub-matrices that
!   almost certainly will not align with the ideal block cyclic data
!   distribution that the PBLAS and ScaLAPACK would prefer. Each sub-matrix
!   represents the appropriate matrix elements of interaction for the set
!   of orbital basis functions from a pair of two atoms. So, once we
!   determine the matrix elements that each process is responsible for, we
!   need to decide how to compute them. Here are some options:
! We could ask each process to compute _only_ those elements that it is
!   responsible for. This turns out to be quite ridiculous because if a sub-
!   matrix crosses a block boundary in the block cyclic distribution then
!   it would require one process to compute _some_ orbital-orbital
!   interactions of an atom pair, but not the others. Then, another process
!   (or processes) would have to compute the other orbital-orbital
!   interactions of that pair. Splitting the work in such a way would
!   require extensive additional logic that would be very expensive to
!   develop, execute, and maintain. The one benefit of this approach is
!   that each process would compute _only_ what it needed and no more.
! Alternatively, we can distribute the list of all atom pairs as well as
!   possible to ensure that each process compute _most_ of what it needs.
!   Then, the portions of each sub-matrix that a process computes but
!   doesn't need must be communicated to another process that did not
!   compute that atom pair, but which needs the results to participate in
!   the block cyclic PBLAS and ScaLAPACK calls. The benefits here are that
!   each atom pair is computed only once and that we do not have to drive
!   any logic to make processes perform partial computation of the set of
!   orbital-orbital interactions for an atom pair. The down side is that
!   each process must send the parts of the block cyclic distribution that
!   its neighbors need (and receive the parts that its neighbors have).
!   As with the first approach, the logic required to make that
!   determination will be tedious and costly.
! The next approach (and the one we adopt) is to again distribute the
!   atom pairs among the processes, but this time to do so in such as way
!   that each process is guaranteed to compute all the elements it needs
!   for its portions of the block cyclic distribution. To make it clear,
!   this will cause some atom pairs to be computed more than once by
!   different processes. Only the portions of the results that fit into the
!   blocks for this process will be kept, the others will be discarded.
!   The benefits of this approach are that each process will compute all
!   elements for its portion of the block cyclic distribution with no
!   additional communication needed and no additional logic will be needed
!   to compute specific partial subsets of anything. The expectation is
!   that the price (repeated computation of a few atom pairs) will be far
!   outweighed by the benefit (no communication and simple logic).
! To execute this approach, we need to determine the ideal block cyclic
!   distribution and then form a list for each process of which atom pairs
!   it must compute to fill all of its blocks. In this module (and this
!   subroutine) we will do the first part. We will define the block cyclic
!   distribution and construct the so-called "array descriptors" that
!   encode the block cyclic distribution. A subroutine in the atomicSites
!   module will use these descriptors to build the atom pair list for
!   each process.
subroutine createBlockCyclicDistribution

   implicit none


end subroutine createBlockCyclicDistribution


! This subroutine will balance a one-dimensional array across a set of MPI
!   processes. The input (toBalance) is the number of array indices that
!   need to be divvied up. The output (initialIdx and finalIdx) are the
!   index numbers in that array that the calling process will be
!   responsible for.
subroutine loadBalMPI(toBalance,initialIdx,finalIdx)
   implicit none
 
   ! Define passed parameters.
   integer, intent(in) :: toBalance ! Quantity to divvy up.
   integer, intent(out) :: initialIdx, finalIdx ! Array indices to do
 
   ! Define local variables.
   integer :: jobsPer, remainder
 
   ! Compute the initial number of indices from the total array that
   !   each process should be accountable for. Also get the remainder
   !   which will always be less than mpiSize.
   jobsPer = int(toBalance / mpiSize)
   remainder = mod(toBalance,mpisize)

   ! Compute the beginning and ending indices of the total array that this
   !   process will do.
   initialIdx = (jobsPer * mpiRank) + 1
   finalIdx = (jobsPer * (mpiRank+1))

   ! If the array is not perfectly distributable among the processes, then
   !   increment the start and end indices such that the higher ranked
   !   processes each take care of one extra array element.

   ! Specifically, first, make process (mpiSize - remainder) start at the
   !   same place, but now do one extra element.
   if (mpiRank == (mpiSize - remainder)) then
      finalIdx = finalIdx + 1
   endif

   ! Then, make any process with rank > (mpiSize - remainder) start at an
   !   index that is higher by the difference between the remainder and (the
   !   difference between the process rank and the total number of
   !   processes). I.e., processes with mpiRank closer to mpiSize will need
   !   to start at an index that is higher. Ultimately, the process with
   !   mpiRank == mpiSize will need to shift its start point by the
   !   magnitude of the remainder and the end point that plus one.
   !   Speaking of which, similar math applies to the end point for the
   !   other processes too.
   if (mpiRank > (mpiSize - remainder)) then
      initialIdx = initialIdx + (remainder - (mpiSize - mpiRank))
      finalIdx = finalIdx + (remainder - (mpiSize - (mpiRank+1)))
   endif

!    This is only needed for C type array indices. i.e. first index=0
!    initialIdx = initialIdx - 1
!    finalIdx = finalIdx - 1

end subroutine loadBalMPI


subroutine checkAttributeHDF5(aid,readTarget,hdf5Status)

   ! Use necessary modules.
   use HDF5

   ! Make sure no variables are implicitly decalred.
   implicit none

   ! Define passed dummy parameters.
   integer(hsize_t), intent(in) :: aid ! attribute id
   character(len=*), intent(in) :: readTarget
   integer, intent(out) :: hdf5Status

   ! Define local variables.
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim
   integer :: hdferr

   ! Assume that the HDF5 attribute indicates that the associated dataset
   !   has not yet been computed.
   hdf5Status = 0
   attribIntDims(1) = 1
   if (mpiRank == 0) then
      call h5aread_f(aid,H5T_NATIVE_INTEGER,hdf5Status,&
            & attribIntDims,hdferr)
   endif
   call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
   if (hdferr /= 0) call stopMPI("Failed to read " // readTarget)

   call MPI_BCAST(hdf5Status,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
   if (hdf5Status == 1) then
      if (mpiRank == 0) then
         write(20,*) readTarget, " already exists. Skipping."
         call h5aclose_f(aid,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) call stopMPI(readtarget)
   endif


end subroutine checkAttributeHDF5


end module O_MPI
