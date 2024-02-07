program GaussianIntegrals

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define program variables.
   real (kind=double), allocatable, dimension (:,:,:,:) :: pc
   real (kind=double), allocatable, dimension (:,:,:,:) :: sh
   real (kind=double), dimension (3) :: alphas
   real (kind=double), dimension (3,3) :: pos ! x,y,z of A,B,C
   real (kind=double), dimension (3,2) :: temp_alphas
   real (kind=double), dimension (3,3,2) :: temp_pos ! x,y,z of A,B,C
   real (kind=double), dimension (3) :: temp_alphas_step
   real (kind=double), dimension (3,3) :: temp_pos_step ! x,y,z of A,B,C
   real (kind=double) :: cell_size, step_size ! Parameters for the numerical
         ! integration process. Read carefully! The cell size is given as a
         ! positive number that defines the maximum extent of the space to
         ! integrate along any given axis. Hence the actual cell is double
         ! that number. I.e., integrate along x, y, z axes from -cell_size
         ! to +cell_size. The step size is just what it appears to be, but
         ! note that it is for numerical integration and is a very different
         ! concept and is unrelated to the "num_steps" defined below.
   integer :: num_segments, num_steps
   integer :: h, i, p, q

   ! Allocate space to hold the appropriately sized pc and sh matrices. The
   !   last index in both pc and sh is a 2 to hold analytical solutions
   !   (value = 1) and numerical solutions (value = 2). The second to last
   !   index is the max of three or (2*max_lam+1) to allow for the momentum
   !   matrix and multi-level nuclear solutions (but all elements are not
   !   always used). The pc matrix is used multiple times to produce the three
   !   momentum matrix sets (x,y,z).
   allocate (pc(1,1,3,2))
   allocate (sh(1,1,3,2))

   ! Open the control file.
   open (10, file="intgcontrol", status='old')

   ! The first bits of information in the control file will control the
   !   numerical integration. The
   read (10,*) cell_size, step_size

   ! Read the number of segments.
   read (10,*) num_segments

   ! For each segment, read in the number of steps, allocate space to hold
   !   relevant information for each step, compute the step descriptions,
   !   compute the integral matrices, store the results, deallocate for the
   !   next segment.
   do h = 1, num_segments

      ! Read in the number of steps for this segment.
      read (10,*) num_steps
      if (num_steps == 1) then
   
         ! Read the parameters for the current segment.
         read (10,*) temp_alphas(:,1), temp_pos(:,:,1)

         ! Fill in dummy values for the temp step data. (It will not be used.)
         temp_alphas_step(:) = 0.0d0
         temp_pos_step(:,:) = 0.0d0
      else
   
         ! Read the beginning and ending parameters for the current segment.
         read (10,*) temp_alphas(:,1), temp_pos(:,:,1) ! First step
         read (10,*) temp_alphas(:,2), temp_pos(:,:,2) ! Last step

         ! Compute the size of the step for each parameter.
         temp_alphas_step(:) = (temp_alphas(:,2) - temp_alphas(:,1))/num_steps
         temp_pos_step(:,:) = (temp_pos(:,:,2) - temp_pos(:,:,1))/num_steps
      endif

      ! Start iterating over the steps.
      do i = 1, num_steps
         alphas(:) = temp_alphas(:,1) + (i-1) * temp_alphas_step(:)
         pos(:,:) = temp_pos(:,:,1) + (i-1) * temp_pos_step(:,:)


         ! Compute the pc and sh integral results for the current parameters.
         call delectron3CIntgAnaBB(alphas(1),alphas(2),alphas(3),pos(:,1),&
               & pos(:,2),pos(:,3),pc(:,:,:,1),sh(:,:,:,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call delectron3CIntgNumBB(alphas(1),alphas(2),alphas(3),pos(:,1),&
               & pos(:,2),pos(:,3),pc(:,:,:,2),sh(:,:,:,2),cell_size,step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,8,alphas,pos,pc(:,:,1,:),sh(:,:,1,:),&
               & "DElexbb.dat")
         call print_pc_sh(h,i,9,alphas,pos,pc(:,:,2,:),sh(:,:,2,:),&
               & "DEleybb.dat")
         call print_pc_sh(h,i,10,alphas,pos,pc(:,:,3,:),sh(:,:,3,:),&
               & "DElezbb.dat")

      enddo
   enddo

   ! Deallocate space.
   deallocate(pc)
   deallocate(sh)


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module/program subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

   subroutine print_pc_sh(h,i,unitPlus,alphas,pos,pc,sh,filename)

      ! Use necessary modules
      use O_Kinds

      ! Make sure no funny variables are created.
      implicit none

      ! Define dummy variables.
      integer :: h, i, unitPlus
      real(kind=double), dimension(3) :: alphas
      real(kind=double), dimension(3,3) :: pos ! xyz (1st idx) of ABC (2nd)
      real(kind=double), dimension(1,1,2) :: pc ! Last idx 1 = ana; 2 = num
      real(kind=double), dimension(1,1,2) :: sh ! Last idx 1 = ana; 2 = num
      character*11 :: filename

      ! Define local variables.
      logical :: io_opened
      integer :: j,k

      ! Open the output file if it isn't already open.
      inquire(299+unitPlus, OPENED = io_opened)
      if (.not. io_opened) then
         open (299+unitPlus, file=filename, status="unknown")
      endif

      write (299+unitPlus,fmt="(2i5)",advance="NO") h, i
      write (299+unitPlus,fmt="(3e16.8)",advance="NO") alphas(:)
      write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,1)
      write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,2)
      write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,3)

      ! Print the result difference for the given pc matrix.
      do j = 1, 1
         do k = 1, 1
            write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,1)
            write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,2)
            write (299+unitPlus,fmt="(e16.8)",advance="NO") &
                  & pc(k,j,1) - pc(k,j,2)
         enddo
      enddo

      ! Print the result difference for the given sh matrix.
      do j = 1, 1
         do k = 1, 1
            write (299+unitPlus,fmt="(e16.8)",advance="NO") &
                  & sh(k,j,1) - sh(k,j,2)
         enddo
      enddo

      ! Write an endline to prepare for the next step.
      write (299+unitPlus, *) ""
   end subroutine

   subroutine delectron3CIntgAnaBB(a1,a2,a3,A,B,C,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16,:): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20,:): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2, a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   real (kind=double), dimension (1,1,3), intent(out) :: pc
   real (kind=double), dimension (1,1,3), intent(out) :: sh

   ! Define local variables.
   real (kind=double), dimension (1,1) :: pc_ol
   real (kind=double), dimension (3) :: P, PA, PB, d
   real (kind=double), dimension (3) :: G, GA, GB, GC, PC_3C
   real (kind=double) :: zeta, inv_2zeta, xi
   real (kind=double) :: zeta3C, inv_2zeta3C, zetaFactor
   real (kind=double) :: preFactorOL
   real (kind=double) :: inv_2zeta_a, inv_2zeta_b

   ! Initialize local variables.
   zeta = a1 + a2
   zeta3C = zeta + a3
   zetaFactor = (a2+a3)/zeta3C
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   inv_2zeta_a = 1.0d0 / (2.0d0 * a1)
   inv_2zeta_b = 1.0d0 / (2.0d0 * a2)
   inv_2zeta3C = 1.0d0 / (2.0d0 * zeta3C)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   G = (a1*A + a2*B + a3*C) / zeta3C
   PA = P - A
   PB = P - B
   GA = G - A
   GB = G - B
   GC = G - C
   PC_3C = P - C
   d = A - B
   preFactorOL = ((pi/zeta3C)**1.5) &
         & * exp(-xi*sum(d**2)-(zeta*a3/zeta3C)*sum(PC_3C**2))

pc_ol(1,1) = preFactorOL

pc(1,1,1) = 2.0d0*(a2+a3)*(GB(1)*pc_ol(1,1))

pc(1,1,2) = 2.0d0*(a2+a3)*(GB(2)*pc_ol(1,1))

pc(1,1,3) = 2.0d0*(a2+a3)*(GB(3)*pc_ol(1,1))

sh(1,1,1) = pc(1,1,1)

sh(1,1,2) = pc(1,1,2)

sh(1,1,3) = pc(1,1,3)


   end subroutine delectron3CIntgAnaBB

   subroutine delectron3CIntgNumBB(a1,a2,a3,A,B,C,pc,sh,cell_size,step_size)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16,3): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20,3): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2, a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   real (kind=double), dimension (1,1,3), intent(out) :: pc
   real (kind=double), dimension (1,1,3), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i, j
   integer :: num_steps
   integer, dimension (1,3) :: triads
   integer, dimension (1,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, curr_pos
   real (kind=double), dimension (3) :: xyz
   real (kind=double), dimension (3,2) :: xyz_I ! Indices=xyz, noprime||prime

   ! Before we proceed with the calculation we need to understand a bit more
   !   about exactly what is being computed. The form of the integration is:
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        d/dR_xyz2  [
   !        exp(-zeta3*( (Px-Rx3)**2 + (Py-Ry3)**2 + (Pz-Rz3)**2 )) *
   !        (Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   ! We use Pxyz for points in space; Rxyz1,2,3 for atomic sites 1,2 and 3;
   !   zeta1,2,3 for the exponential decay factors; lxyz1,2 for the angular
   !   momentum numbers for atoms 1 and 2,and lxyz3=0 (Always s-type) for the 
   !   angular momentum number for atom 3, and d/d_xyz for independent
   !   derivatives (i.e.,we have three separate triple integrals).
   !   The integral must be computed over all space for all different possible
   !   values of lxyz1,2 for some arbitrarily chosen test coordinates and
   !   decay rates.
   ! Note importantly, that in this case we are considering the solution for
   !   the scenario when the middle term is at site "B" and the last term is
   !   also at site "B". The B is in quotes to emphasize the following idea:
   !   Consider that we want to take the derivative of the Hamiltonian with
   !   respect to the coordinate of, say, atom #7 in order to get the force on
   !   atom #7. For the three center overlap contribution to the Hamiltonian,
   !   (commonly expressed as <A|C|B>) we find that atom #7 will be used in
   !   each of the A, C, and B terms. Sometimes it will be <7|C|B> with
   !   arbitrary other atoms for C and B. Other times it will be <A|7|B> or
   !   <A|C|7> or <A|7|7>. Presently, we are dealing with the <A|7|7> case.
   ! Note, importantly, that the derivative is with respect to the atomic
   !   coordinate of both the middle term and the second basis function. We
   !   are dealing with the scenario when C=B. Therefore, we will rewrite the
   !   integral now to reflect that.
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        d/dR_xyz2  [
   !        exp(-zeta3*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 )) *
   !        (Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   ! Let's now rewrite this one more time to condense terms.
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        d/dR_xyz2  [ exp((-zeta3 - zeta2)*
   !        ((Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 )) *
   !        (Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 ] dxdydz
   ! Now, focusing on just the d/dRx2 version of the set of triple integrals,
   !   we will pull all terms with Py and Pz out of the dx integral to get:
   ! SS { [(Py-Ry1)**ly1 * (Pz-Rz1)**lz1 *
   !       exp(-zeta1*((Py-Ry1)**2 + (Pz-Rz1)**2 ))]*
   !       [exp((-zeta3-zeta2)*((Py-Ry2)**2 + (Pz-Rz2)**2 ))]*
   !       [(Py-Ry2)**ly2 * (Pz-Rz2)**lz2] *
   !     S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * d/dRx2 [
   !       exp((-zeta3-zeta2)*((Px-Rx3)**2)) * (Px-Rx2)**lx2 ]
   !       ] dx}
   !     dydz
   ! Applying the derivative, the internal 1D dx integral has the form:
   !   Ix' = S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) *                
   !           [- lx2*(Px-Rx2)**(lx2-1)
   !            + 2*zeta2*zeta3*(Px-Rx2)**(lx2+1)]
   !           * exp((-zeta2-zeta3)*(Px-Rx2)**2)
   !        ]
   ! Each of the other integrals (Iy, Iz) will have the form of a simple
   !   1D overlap integral (without a derivative):
   !   Iy = S [(Py-Ry1)**ly1 * exp(-zeta1*((Py-Ry1)**2)) *
   !           (Py-Ry2)**ly2 * exp(-zeta2*((Py-Ry2)**2)) * 
   !                           exp(-zeta3*((Py-Ry3)**2))] dy
   !   Iz = S [(Pz-Rz1)**lz1 * exp(-zeta1*((Pz-Rz1)**2)) *
   !           (Pz-Rz2)**lz2 * exp(-zeta2*((Pz-Rz2)**2)) *
   !                           exp(-zeta3*((Pz-Rz3)**2))] dz
   ! The total integral !! for the d/dx version !! is thus Ix' * Iy * Iz.
   !   Similarly, the total integrals for d/dy and d/dz will have the form:
   !   Ix*Iy'*Iz and Ix*Iy*Iz' where the Iy' and Iz' are the appropriate
   !   analogs of the Ix' and the Ix is the analog of the Iy or Iz.
   ! Thus, every integral is a product of independent 1D integrals and the
   !   solution of the total integral is a set of products of integrals that
   !   can all be computed at once.
   ! With regard to the term in each integral that has the **(lx2-1) form
   !   (or **(ly2-1) or **(lz2-1)), some special care must be taken. This
   !   term represents an angular momentum and the integral will have the
   !   form of an overlap integral. Because we cannot have a negative angular
   !   momentum therefore we must discard this term when lx2-1 < 0.

   ! Initialize local variables.
   triads(1,:) = (/0,0,0/)
   conversion(1,1,:) = (/1,0,0/)
   conversion(1,2,:) = (/1,1,1/)

   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   do p = 1, 1
      do q = 1, 1

         ! Assign l1 and l2 values for each gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize sum variables.
         xyz_I(:,:) = 0.0d0

         ! Start a loop over x coordinates.
         do i = 0, num_steps
            curr_pos = start_pos + (i*step_size)

            ! Compute the no-prime integrals first.
            do j = 1, 3
               xyz_I(j,2) = xyz_I(j,2) &
                     & + noPrimeDElectronBB(step_size, curr_pos, A(j), &
                     & B(j), C(j), a1, a2, a3, l1(j), l2(j))
            enddo

            ! Compute the prime integrals second.
            do j = 1, 3
               xyz_I(j,1) = xyz_I(j,1) &
                     & + primeDElectronBB(step_size, curr_pos, A(j), &
                     & B(j), C(j), a1, a2, a3, l1(j), l2(j))
            enddo
         enddo

         pc(q,p,1) = xyz_I(1,1)*xyz_I(2,2)*xyz_I(3,2)
         pc(q,p,2) = xyz_I(1,2)*xyz_I(2,1)*xyz_I(3,2)
         pc(q,p,3) = xyz_I(1,2)*xyz_I(2,2)*xyz_I(3,1)
      enddo
    enddo

sh(1,1,1) = pc(1,1,1)

sh(1,1,2) = pc(1,1,2)

sh(1,1,3) = pc(1,1,3)


   end subroutine delectron3CIntgNumBB


   function noPrimeDElectronBB(step_size, curr_pos, A, B, C, &
         & a1, a2, a3, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos
      real (kind=double), intent(in) :: A, B, C, a1, a2, a3
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: noPrimeDElectronBB

      ! Compute the integral part.
      noPrimeDElectronBB = step_size &
            & * (curr_pos - A)**l1 * (curr_pos - B)**l2 &
            & * exp(-a1*(curr_pos - A)**2) * exp(-a2*(curr_pos - B)**2) &
            & * exp(-a3*(curr_pos - B)**2)

      return
 
   end function noPrimeDElectronBB


   function primeDElectronBB(step_size, curr_pos, A, B, C, a1, a2, a3, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos
      real (kind=double), intent(in) :: A, B, C, a1, a2, a3
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: primeDElectronBB

      ! Compute each "internal" term of the prime integral. Compare each of
      !   these lines with the 1D electron derivative matrix equation produced
      !   by the osrecurintg_makenum.py script (appropriately separated into
      !   terms) and the expression in Nuha's dissertation. FIX: Check
      !   equation number from Nuha's dissertation.

      ! Last term in the sympy script-produced equation. FIX: As above.
      primeDElectronBB = 2.0d0 * (a2 + a3) * (curr_pos - B)**(l2+1)

      ! First line in the equation produced by sympy. As mentioned above,
      !   we only compute this if the angular momentum will allow it.
      if (l2 >= 1) then
         primeDElectronBB = primeDElectronBB - l2 * (curr_pos - B)**(l2-1)
      endif

      ! Multiply prime integral by the preceeding primitive gaussian
      !   coefficient and exponential and multiply by the succeeding
      !   exponential. (We have already multiplied by the succeeding
      !   primitive gaussian coefficient in the above lines.)
      primeDElectronBB = primeDElectronBB &
            & * (curr_pos-A)**l1 * exp(-a1*(curr_pos-A)**2) &
            & * exp(-a2*(curr_pos-B)**2) * exp(-a3*(curr_pos-B)**2)

      ! Finally, multiply by the step size.
      primeDElectronBB = primeDElectronBB * step_size

      return

   end function primeDElectronBB

end program GaussianIntegrals
