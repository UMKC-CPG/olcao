! Author:  Ben Walker and Paul Rulis
!
! This program is a numerical integrator that provides test numbers
! for Two-Electron Integrals, AKA Electron Repulsion Integrals.
! This version of the program is intended to output the numerical
! integration solutions for the 20x20x20x20 array of basis functions.

program numerical4CRepulsionAll

   ! Load necessary modules.
   use mpi
   use O_kinds

   ! Make sure that nothing funny is accidentally declared.
   implicit none

   ! Declare variables
   real (kind=double) :: tempMatrix
   real (kind=double) :: minRange
   real (kind=double) :: stepSize
   real (kind=double) :: currentPos
   real (kind=double) :: eps !epsilon
   real (kind=double) :: point1, point2
   real (kind=double) :: x1, x2, y1, y2, z1, z2, f, fSum, r12
   real (kind=double) :: a1, a2, a3, a4
   real (kind=double), dimension(3) :: A, B, C, D
   real (kind=double), dimension(4) :: alphas
   real (kind=double), dimension(6) :: prod
   real (kind=double), dimension(4,3) :: sites
   real (kind=double), allocatable, dimension (:) :: fSumResults
   real (kind=double), allocatable, dimension (:) :: fSumResultsMaster
   real (kind=double), allocatable, dimension (:,:) :: r12Matrix
   real (kind=double), allocatable, dimension (:,:,:,:) :: coeffs
   integer :: counter
   integer :: matrixRange
   integer :: meshPoints
   integer :: i, j, k, l, m, n, p, q, r, s, o
   integer, dimension(3) :: l1, l2, l3, l4
   integer, dimension(20,3) :: array
   integer, allocatable, dimension (:) :: oIndex

   ! Declare parallel variables
   integer :: mpirank, mpisize, mpierr, mpistatus(MPI_STATUS_SIZE)
   integer :: tagFSum, tagOIndex, workPerProcess, beginRange, endRange

   ! Call MPI initialization subroutines
   call MPI_INIT(mpierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, mpierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)
   call MPI_BARRIER(MPI_COMM_WORLD, mpierr)

   tagFSum   = 1
   tagOIndex = 2

   eps = 0.00000001_double

   ! Set up the triads
   array(1,:)  = (/0,0,0/)
   array(2,:)  = (/1,0,0/)
   array(3,:)  = (/0,1,0/)
   array(4,:)  = (/0,0,1/)
   array(5,:)  = (/2,0,0/)
   array(6,:)  = (/0,2,0/)
   array(7,:)  = (/0,0,2/)
   array(8,:)  = (/1,1,0/)
   array(9,:)  = (/1,0,1/)
   array(10,:) = (/0,1,1/)
   array(11,:) = (/1,1,1/)
   array(12,:) = (/2,1,0/)
   array(13,:) = (/2,0,1/)
   array(14,:) = (/1,2,0/)
   array(15,:) = (/0,2,1/)
   array(16,:) = (/1,0,2/)
   array(17,:) = (/0,1,2/)
   array(18,:) = (/3,0,0/)
   array(19,:) = (/0,3,0/)
   array(20,:) = (/0,0,3/)
 
   ! Read all the data.
   open (unit=50,file='input',form='formatted')
   read (50,*) sites(1,1), sites(1,2), sites(1,3)
   read (50,*) sites(2,1), sites(2,2), sites(2,3)
   read (50,*) sites(3,1), sites(3,2), sites(3,3)
   read (50,*) sites(4,1), sites(4,2), sites(4,3)
   read (50,*) alphas(1), alphas(2), alphas(3), alphas(4)
   read (50,*) stepSize, matrixRange, minRange
   close (50)

   ! Compute the number of mesh points. 
   meshPoints  = int(-minRange*2.0_double/stepSize)
   if (mpirank == 0) then
      write (6,*) "Num mesh points  = ",meshPoints
      call flush (6)
   endif
 
   allocate (coeffs(meshPoints+1,4,4,3))
   allocate (r12Matrix(meshPoints+1,meshPoints+1))

   ! Precompute all possible distance values into a symmetric r12 matrix.
   do i = 1, meshPoints+1
      point1 = (minRange + (i*stepSize))
      do j = 1, meshPoints+1
         point2 = (minRange + (j*stepSize))
         r12Matrix(j,i) = (point2 - point1)**2
      enddo
   enddo

   ! Precompute all possible terms.
   ! Note that the order of the loop indices is a bit atypical here so that
   !   they can be accessed in an optimal sequence in the later (real)
   !   computation.
   do i = 1,meshPoints+1
      currentPos = minRange + ((i-1)*stepSize)
   ! Index guide for coeffs:
   !   Index 1: meshPoints
   !   Index 2: 1-4 for the QN_l values 0-3.
   !   Index 3: 1-4 for centers A,B,C,D. Note A1 is x coordinate of center A.
   !   Index 4: 1-3 for x, y, z
      do j = 1, 3 ! x,y,z
         do k = 1,4 ! Sites A,B,C,D
             coeffs(i,1,k,j) = 1.0_double                   ! e.g. (x1 - A1)**0
             coeffs(i,2,k,j) = (currentPos - sites(k,j))    ! e.g. (x1 - A1)**1
             coeffs(i,3,k,j) = (currentPos - sites(k,j))**2 ! e.g. (x1 - A1)**2
             coeffs(i,4,k,j) = (currentPos - sites(k,j))**3 ! e.g. (x1 - A1)**3
         enddo
         do k = 1,4 ! Sites A,B,C,D (and associated alphas)
            if (alphas(k)*coeffs(i,3,k,j) > 15.0_double) then
               coeffs(i,:,k,j) = 0.0_double
            else
               coeffs(i,:,k,j) = coeffs(i,:,k,j)*exp(-alphas(k)*coeffs(i,3,k,j))
            endif
         enddo
      enddo
   enddo

   ! Determine the amount of work for each process.  This is the number of
   !   integral solutions to solve. For s,p,d,f there are 20x20x20x20 total
   !   integrals to solve for that will be evenly divided among the processes.
   workPerProcess = (matrixRange**4)/mpisize
   if ((mod(matrixRange**4,mpisize) /= 0) .and. (mpirank == 0)) then
      write (6,*) "Make sure that you ask for a number of processors that"
      write (6,*) "divides evenly into",matrixRange**4,"."
      write (6,*) "This will give you proper load balancing."
      stop
   endif

   ! Record the number of steps per process.
   if (mpirank == 0) then
      write (6,*) "Work per process = ",workPerProcess
      call flush (6)
   endif

   ! Determine the range for this process.
   beginRange = mpirank*workPerProcess + 1
   endRange   = beginRange + workPerProcess - 1

   ! Allocate space to hold the result.  For every integral we will save the
   !   actual integral value and the value of the "o" loop index associated
   !   with it.
   allocate (fSumResults(workPerProcess))
   allocate (oIndex(workPerProcess))
   if (mpirank == 0) then
      allocate (fSumResultsMaster(matrixRange**4))
   endif

   ! Start processing integrals.
   do o = beginRange, endRange

      ! Assume that all the processes will progress at about the same speed as
      !   process zero.
      if (mpirank == 0) then
         if (mod(o,10) .eq. 0) then
            write (6,ADVANCE="NO",FMT="(a1)") "|"
         else
            write (6,ADVANCE="NO",FMT="(a1)") "."
         endif
         if (mod(o,50) .eq. 0) then
            write (6,*) " ",o
         endif
         call flush (6)
      endif

      ! Compute the values of p,q,r,s. Note the ordering.
!CHECK
      s = mod(o-1,matrixRange**1) / matrixrange**0 + 1
      r = mod(o-1,matrixRange**2) / matrixRange**1 + 1
      q = mod(o-1,matrixRange**3) / matrixRange**2 + 1
      p = mod(o-1,matrixRange**4) / matrixRange**3 + 1

      ! Assign QN_lxyz values for each dimension (xyz) for each
      !   center.  Note the order of elements.
!CHECK
      l1(:) = array(p,:) + 1
      l2(:) = array(q,:) + 1
      l3(:) = array(r,:) + 1
      l4(:) = array(s,:) + 1

      fSum = 0.0d0 ! Initialize fSum to zero

      do i = 1, meshPoints+1
         prod(1) = coeffs(i,l1(1),1,1)*coeffs(i,l2(1),2,1)
         if (prod(1) == 0.0_double) cycle

      do j = 1, meshPoints+1
         prod(2) = coeffs(j,l3(1),3,1)*coeffs(j,l4(1),4,1)
         if (prod(2) == 0.0_double) cycle

      do k = 1, meshPoints+1
         prod(3) = coeffs(k,l1(2),1,2)*coeffs(k,l2(2),2,2)
         if (prod(3) == 0.0_double) cycle

      do l = 1, meshPoints+1
         prod(4) = coeffs(l,l3(2),3,2)*coeffs(l,l4(2),4,2)
         if (prod(4) == 0.0_double) cycle

         ! Precompute the sum of two distance components of r12.
         tempMatrix = r12Matrix(j,i) + r12Matrix(l,k)

      do m = 1, meshPoints+1
         prod(5) = coeffs(m,l1(3),1,3)*coeffs(m,l2(3),2,3)
         if (prod(5) == 0.0_double) cycle

      do n = 1, meshPoints+1
         ! Compute the square of r12. The root is taken later only when needed.
         r12 = tempMatrix + r12Matrix(n,m)

         if (r12 > eps) then
            prod(6) = coeffs(n,l3(3),3,3)*coeffs(n,l4(3),4,3)
            if (prod(6) == 0.0_double) cycle

            f = product(prod(:))/sqrt(r12)
            fSum = fSum + f
         end if
      end do ! n
      end do ! m
      end do ! l
      end do ! k
      end do ! j
      end do ! i

      ! Store the result from this process. It will be transmitted to process
      !   0 later.
      fSumResults(o-beginRange+1) = fSum
      oIndex(o-beginRange+1) = o

   enddo ! o

   if (mpirank == 0) then

      write (6,*)
      write (6,*) "Accumulating results..."
      call flush (6)

      ! Copy computed results into the master list.
      do i = 1, workPerProcess
         fSumResultsMaster(i) = fSumResults(i)*stepSize**6
      enddo

      do i = 1, mpisize-1

         ! Receive two messages from process i. The first contains the fSum
         !   the second contains the associated array index values.  Note that
         !   the buffer is the same array that held the rank=0 data just copied
         !   into the master list above.
         call MPI_RECV(fSumResults,workPerProcess,MPI_DOUBLE_PRECISION,&
               & i,tagFSum,MPI_COMM_WORLD,mpistatus,mpierr)
         call MPI_RECV(oIndex,workPerProcess,MPI_INTEGER,&
               & i,tagOIndex,MPI_COMM_WORLD,mpistatus,mpierr)

         ! Insert the result into the master list.
         do j = 1, workPerProcess
            fSumResultsMaster(oIndex(j)) = fSumResults(j)*stepSize**6
         enddo
      enddo
   else

      ! Send the results to process 0.
      call MPI_SEND(fSumResults,workPerProcess,MPI_DOUBLE_PRECISION,0,tagFSum,&
            & MPI_COMM_WORLD,mpierr)
      call MPI_SEND(oIndex,workPerProcess,MPI_INTEGER,0,tagOIndex,&
            & MPI_COMM_WORLD,mpierr)
   endif

   if (mpirank == 0) then

      ! Record the result.
      counter = 0
      do i = 1, matrixRange
         do j = 1, matrixRange
            do k = 1, matrixRange
               do l = 1, matrixRange
                  counter = counter + 1
                  write (6,fmt="(4i5,e17.8)") i,j,k,l,fSumResultsMaster(counter)
                  call flush (6)
               enddo
            enddo
         enddo
      enddo
   endif

   call MPI_FINALIZE(mpierr)
end program numerical4CRepulsionAll
