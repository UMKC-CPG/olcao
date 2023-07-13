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
   allocate (pc(20,20,7,2))
   allocate (sh(16,16,3,2))

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


         ! Compute the pc and sh integral results for the current parameters
         !   using the analytic formulas.
         call overlap2CIntgAna(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,1,1),sh(:,:,1,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call overlap2CIntgNum(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,1,2),sh(:,:,1,2),cell_size,step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,1,alphas,pos,pc(:,:,1,:),sh(:,:,1,:),&
               & "overlap.dat")


         ! Compute the pc and sh integral results for the current parameters
         !   using analytical formulas.
         call kinetic2CIntgAna(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,1,1),sh(:,:,1,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call kinetic2CIntgNum(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,1,2),sh(:,:,1,2),cell_size,step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,2,alphas,pos,pc(:,:,1,:),sh(:,:,1,:),&
               & "kinetic.dat")


         ! Compute the pc and sh integral results for the current parameters
         !   using analytical formulas.
         call nuclear3CIntgAna(alphas(1),alphas(2),alphas(3),pos(:,1),&
               & pos(:,2),pos(:,3),pc(:,:,:,1),sh(:,:,1,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call nuclear3CIntgNum(alphas(1),alphas(2),alphas(3),pos(:,1),&
               & pos(:,2),pos(:,3),pc(:,:,1,2),sh(:,:,1,2),cell_size,&
               & step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,3,alphas,pos,pc(:,:,1,:),sh(:,:,1,:),&
               & "nuclear.dat")


         ! Compute the pc and sh integral results for the current parameters
         !   using analytical formulas.
         call electron3CIntgAna(alphas(1),alphas(2),alphas(3),pos(:,1),&
               & pos(:,2),pos(:,3),pc(:,:,1,1),sh(:,:,1,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call electron3CIntgNum(alphas(1),alphas(2),alphas(3),pos(:,1),&
               & pos(:,2),pos(:,3),pc(:,:,1,2),sh(:,:,1,2),cell_size,&
               & step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,4,alphas,pos,pc(:,:,1,:),sh(:,:,1,:),&
               & "electrn.dat")


         ! Compute the pc and sh integral results for the current parameters.
         call momentum2CIntgAna(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,:,1),sh(:,:,:,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call momentum2CIntgNum(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,:,2),sh(:,:,:,2),cell_size,step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,5,alphas,pos,pc(:,:,1,:),sh(:,:,1,:),&
               & "momentx.dat")
         call print_pc_sh(h,i,6,alphas,pos,pc(:,:,2,:),sh(:,:,2,:),&
               & "momenty.dat")
         call print_pc_sh(h,i,7,alphas,pos,pc(:,:,3,:),sh(:,:,3,:),&
               & "momentz.dat")

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
      real(kind=double), dimension(20,20,2) :: pc ! Last idx 1 = ana; 2 = num
      real(kind=double), dimension(16,16,2) :: sh ! Last idx 1 = ana; 2 = num
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
      do j = 1, 20
         do k = 1, 20
            write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,1)
            write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,2)
            write (299+unitPlus,fmt="(e16.8)",advance="NO") &
                  & pc(k,j,1) - pc(k,j,2)
         enddo
      enddo

      ! Print the result difference for the given sh matrix.
      do j = 1, 16
         do k = 1, 16
            write (299+unitPlus,fmt="(e16.8)",advance="NO") &
                  & sh(k,j,1) - sh(k,j,2)
         enddo
      enddo

      ! Write an endline to prepare for the next step.
      write (299+unitPlus, *) ""
   end subroutine

   subroutine overlap2CIntgAna(a1,a2,A,B,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   real (kind=double), dimension (20,20), intent(out) :: pc
   real (kind=double), dimension (16,16), intent(out) :: sh

   ! Define local variables.
   real (kind=double), dimension (3) :: P, PA, PB, d
   real (kind=double) :: zeta, inv_2zeta, xi, preFactorOL

   ! Initialize local variables.
   zeta = a1 + a2
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   d = A - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d**2))

pc(1,1) = preFactorOL

pc(2,1) = PA(1)*preFactorOL

pc(3,1) = PA(2)*preFactorOL

pc(4,1) = PA(3)*preFactorOL

pc(5,1) = PA(1)*pc(2,1) + inv_2zeta*(preFactorOL)

pc(6,1) = PA(2)*pc(3,1) + inv_2zeta*(preFactorOL)

pc(7,1) = PA(3)*pc(4,1) + inv_2zeta*(preFactorOL)

pc(8,1) = PA(2)*pc(2,1)

pc(9,1) = PA(3)*pc(2,1)

pc(10,1) = PA(3)*pc(3,1)

pc(11,1) = PA(3)*pc(8,1)

pc(12,1) = PA(2)*pc(5,1)

pc(13,1) = PA(3)*pc(5,1)

pc(14,1) = PA(1)*pc(6,1)

pc(15,1) = PA(3)*pc(6,1)

pc(16,1) = PA(1)*pc(7,1)

pc(17,1) = PA(2)*pc(7,1)

pc(18,1) = PA(1)*pc(5,1) + inv_2zeta*(2*pc(2,1))

pc(19,1) = PA(2)*pc(6,1) + inv_2zeta*(2*pc(3,1))

pc(20,1) = PA(3)*pc(7,1) + inv_2zeta*(2*pc(4,1))

pc(1,2) = PB(1)*preFactorOL

pc(2,2) = PA(1)*pc(1,2) + inv_2zeta*(preFactorOL)

pc(3,2) = PA(2)*pc(1,2)

pc(4,2) = PA(3)*pc(1,2)

pc(5,2) = PA(1)*pc(2,2) + inv_2zeta*(pc(1,2) + pc(2,1))

pc(6,2) = PA(2)*pc(3,2) + inv_2zeta*(pc(1,2))

pc(7,2) = PA(3)*pc(4,2) + inv_2zeta*(pc(1,2))

pc(8,2) = PA(2)*pc(2,2)

pc(9,2) = PA(3)*pc(2,2)

pc(10,2) = PA(3)*pc(3,2)

pc(11,2) = PA(3)*pc(8,2)

pc(12,2) = PA(2)*pc(5,2)

pc(13,2) = PA(3)*pc(5,2)

pc(14,2) = PA(1)*pc(6,2) + inv_2zeta*(pc(6,1))

pc(15,2) = PA(3)*pc(6,2)

pc(16,2) = PA(1)*pc(7,2) + inv_2zeta*(pc(7,1))

pc(17,2) = PA(2)*pc(7,2)

pc(18,2) = PA(1)*pc(5,2) + inv_2zeta*(2*pc(2,2) + pc(5,1))

pc(19,2) = PA(2)*pc(6,2) + inv_2zeta*(2*pc(3,2))

pc(20,2) = PA(3)*pc(7,2) + inv_2zeta*(2*pc(4,2))

pc(1,3) = PB(2)*preFactorOL

pc(2,3) = PA(1)*pc(1,3)

pc(3,3) = PA(2)*pc(1,3) + inv_2zeta*(preFactorOL)

pc(4,3) = PA(3)*pc(1,3)

pc(5,3) = PA(1)*pc(2,3) + inv_2zeta*(pc(1,3))

pc(6,3) = PA(2)*pc(3,3) + inv_2zeta*(pc(1,3) + pc(3,1))

pc(7,3) = PA(3)*pc(4,3) + inv_2zeta*(pc(1,3))

pc(8,3) = PA(2)*pc(2,3) + inv_2zeta*(pc(2,1))

pc(9,3) = PA(3)*pc(2,3)

pc(10,3) = PA(3)*pc(3,3)

pc(11,3) = PA(3)*pc(8,3)

pc(12,3) = PA(2)*pc(5,3) + inv_2zeta*(pc(5,1))

pc(13,3) = PA(3)*pc(5,3)

pc(14,3) = PA(1)*pc(6,3)

pc(15,3) = PA(3)*pc(6,3)

pc(16,3) = PA(1)*pc(7,3)

pc(17,3) = PA(2)*pc(7,3) + inv_2zeta*(pc(7,1))

pc(18,3) = PA(1)*pc(5,3) + inv_2zeta*(2*pc(2,3))

pc(19,3) = PA(2)*pc(6,3) + inv_2zeta*(2*pc(3,3) + pc(6,1))

pc(20,3) = PA(3)*pc(7,3) + inv_2zeta*(2*pc(4,3))

pc(1,4) = PB(3)*preFactorOL

pc(2,4) = PA(1)*pc(1,4)

pc(3,4) = PA(2)*pc(1,4)

pc(4,4) = PA(3)*pc(1,4) + inv_2zeta*(preFactorOL)

pc(5,4) = PA(1)*pc(2,4) + inv_2zeta*(pc(1,4))

pc(6,4) = PA(2)*pc(3,4) + inv_2zeta*(pc(1,4))

pc(7,4) = PA(3)*pc(4,4) + inv_2zeta*(pc(1,4) + pc(4,1))

pc(8,4) = PA(2)*pc(2,4)

pc(9,4) = PA(3)*pc(2,4) + inv_2zeta*(pc(2,1))

pc(10,4) = PA(3)*pc(3,4) + inv_2zeta*(pc(3,1))

pc(11,4) = PA(3)*pc(8,4) + inv_2zeta*(pc(8,1))

pc(12,4) = PA(2)*pc(5,4)

pc(13,4) = PA(3)*pc(5,4) + inv_2zeta*(pc(5,1))

pc(14,4) = PA(1)*pc(6,4)

pc(15,4) = PA(3)*pc(6,4) + inv_2zeta*(pc(6,1))

pc(16,4) = PA(1)*pc(7,4)

pc(17,4) = PA(2)*pc(7,4)

pc(18,4) = PA(1)*pc(5,4) + inv_2zeta*(2*pc(2,4))

pc(19,4) = PA(2)*pc(6,4) + inv_2zeta*(2*pc(3,4))

pc(20,4) = PA(3)*pc(7,4) + inv_2zeta*(2*pc(4,4) + pc(7,1))

pc(1,5) = PB(1)*pc(1,2) + inv_2zeta*(preFactorOL)

pc(2,5) = PA(1)*pc(1,5) + inv_2zeta*(2*pc(1,2))

pc(3,5) = PA(2)*pc(1,5)

pc(4,5) = PA(3)*pc(1,5)

pc(5,5) = PA(1)*pc(2,5) + inv_2zeta*(pc(1,5) + 2*pc(2,2))

pc(6,5) = PA(2)*pc(3,5) + inv_2zeta*(pc(1,5))

pc(7,5) = PA(3)*pc(4,5) + inv_2zeta*(pc(1,5))

pc(8,5) = PA(2)*pc(2,5)

pc(9,5) = PA(3)*pc(2,5)

pc(10,5) = PA(3)*pc(3,5)

pc(11,5) = PA(3)*pc(8,5)

pc(12,5) = PA(2)*pc(5,5)

pc(13,5) = PA(3)*pc(5,5)

pc(14,5) = PA(1)*pc(6,5) + inv_2zeta*(2*pc(6,2))

pc(15,5) = PA(3)*pc(6,5)

pc(16,5) = PA(1)*pc(7,5) + inv_2zeta*(2*pc(7,2))

pc(17,5) = PA(2)*pc(7,5)

pc(18,5) = PA(1)*pc(5,5) + inv_2zeta*(2*pc(2,5) + 2*pc(5,2))

pc(19,5) = PA(2)*pc(6,5) + inv_2zeta*(2*pc(3,5))

pc(20,5) = PA(3)*pc(7,5) + inv_2zeta*(2*pc(4,5))

pc(1,6) = PB(2)*pc(1,3) + inv_2zeta*(preFactorOL)

pc(2,6) = PA(1)*pc(1,6)

pc(3,6) = PA(2)*pc(1,6) + inv_2zeta*(2*pc(1,3))

pc(4,6) = PA(3)*pc(1,6)

pc(5,6) = PA(1)*pc(2,6) + inv_2zeta*(pc(1,6))

pc(6,6) = PA(2)*pc(3,6) + inv_2zeta*(pc(1,6) + 2*pc(3,3))

pc(7,6) = PA(3)*pc(4,6) + inv_2zeta*(pc(1,6))

pc(8,6) = PA(2)*pc(2,6) + inv_2zeta*(2*pc(2,3))

pc(9,6) = PA(3)*pc(2,6)

pc(10,6) = PA(3)*pc(3,6)

pc(11,6) = PA(3)*pc(8,6)

pc(12,6) = PA(2)*pc(5,6) + inv_2zeta*(2*pc(5,3))

pc(13,6) = PA(3)*pc(5,6)

pc(14,6) = PA(1)*pc(6,6)

pc(15,6) = PA(3)*pc(6,6)

pc(16,6) = PA(1)*pc(7,6)

pc(17,6) = PA(2)*pc(7,6) + inv_2zeta*(2*pc(7,3))

pc(18,6) = PA(1)*pc(5,6) + inv_2zeta*(2*pc(2,6))

pc(19,6) = PA(2)*pc(6,6) + inv_2zeta*(2*pc(3,6) + 2*pc(6,3))

pc(20,6) = PA(3)*pc(7,6) + inv_2zeta*(2*pc(4,6))

pc(1,7) = PB(3)*pc(1,4) + inv_2zeta*(preFactorOL)

pc(2,7) = PA(1)*pc(1,7)

pc(3,7) = PA(2)*pc(1,7)

pc(4,7) = PA(3)*pc(1,7) + inv_2zeta*(2*pc(1,4))

pc(5,7) = PA(1)*pc(2,7) + inv_2zeta*(pc(1,7))

pc(6,7) = PA(2)*pc(3,7) + inv_2zeta*(pc(1,7))

pc(7,7) = PA(3)*pc(4,7) + inv_2zeta*(pc(1,7) + 2*pc(4,4))

pc(8,7) = PA(2)*pc(2,7)

pc(9,7) = PA(3)*pc(2,7) + inv_2zeta*(2*pc(2,4))

pc(10,7) = PA(3)*pc(3,7) + inv_2zeta*(2*pc(3,4))

pc(11,7) = PA(3)*pc(8,7) + inv_2zeta*(2*pc(8,4))

pc(12,7) = PA(2)*pc(5,7)

pc(13,7) = PA(3)*pc(5,7) + inv_2zeta*(2*pc(5,4))

pc(14,7) = PA(1)*pc(6,7)

pc(15,7) = PA(3)*pc(6,7) + inv_2zeta*(2*pc(6,4))

pc(16,7) = PA(1)*pc(7,7)

pc(17,7) = PA(2)*pc(7,7)

pc(18,7) = PA(1)*pc(5,7) + inv_2zeta*(2*pc(2,7))

pc(19,7) = PA(2)*pc(6,7) + inv_2zeta*(2*pc(3,7))

pc(20,7) = PA(3)*pc(7,7) + inv_2zeta*(2*pc(4,7) + 2*pc(7,4))

pc(1,8) = PB(2)*pc(1,2)

pc(2,8) = PA(1)*pc(1,8) + inv_2zeta*(pc(1,3))

pc(3,8) = PA(2)*pc(1,8) + inv_2zeta*(pc(1,2))

pc(4,8) = PA(3)*pc(1,8)

pc(5,8) = PA(1)*pc(2,8) + inv_2zeta*(pc(1,8) + pc(2,3))

pc(6,8) = PA(2)*pc(3,8) + inv_2zeta*(pc(1,8) + pc(3,2))

pc(7,8) = PA(3)*pc(4,8) + inv_2zeta*(pc(1,8))

pc(8,8) = PA(2)*pc(2,8) + inv_2zeta*(pc(2,2))

pc(9,8) = PA(3)*pc(2,8)

pc(10,8) = PA(3)*pc(3,8)

pc(11,8) = PA(3)*pc(8,8)

pc(12,8) = PA(2)*pc(5,8) + inv_2zeta*(pc(5,2))

pc(13,8) = PA(3)*pc(5,8)

pc(14,8) = PA(1)*pc(6,8) + inv_2zeta*(pc(6,3))

pc(15,8) = PA(3)*pc(6,8)

pc(16,8) = PA(1)*pc(7,8) + inv_2zeta*(pc(7,3))

pc(17,8) = PA(2)*pc(7,8) + inv_2zeta*(pc(7,2))

pc(18,8) = PA(1)*pc(5,8) + inv_2zeta*(2*pc(2,8) + pc(5,3))

pc(19,8) = PA(2)*pc(6,8) + inv_2zeta*(2*pc(3,8) + pc(6,2))

pc(20,8) = PA(3)*pc(7,8) + inv_2zeta*(2*pc(4,8))

pc(1,9) = PB(3)*pc(1,2)

pc(2,9) = PA(1)*pc(1,9) + inv_2zeta*(pc(1,4))

pc(3,9) = PA(2)*pc(1,9)

pc(4,9) = PA(3)*pc(1,9) + inv_2zeta*(pc(1,2))

pc(5,9) = PA(1)*pc(2,9) + inv_2zeta*(pc(1,9) + pc(2,4))

pc(6,9) = PA(2)*pc(3,9) + inv_2zeta*(pc(1,9))

pc(7,9) = PA(3)*pc(4,9) + inv_2zeta*(pc(1,9) + pc(4,2))

pc(8,9) = PA(2)*pc(2,9)

pc(9,9) = PA(3)*pc(2,9) + inv_2zeta*(pc(2,2))

pc(10,9) = PA(3)*pc(3,9) + inv_2zeta*(pc(3,2))

pc(11,9) = PA(3)*pc(8,9) + inv_2zeta*(pc(8,2))

pc(12,9) = PA(2)*pc(5,9)

pc(13,9) = PA(3)*pc(5,9) + inv_2zeta*(pc(5,2))

pc(14,9) = PA(1)*pc(6,9) + inv_2zeta*(pc(6,4))

pc(15,9) = PA(3)*pc(6,9) + inv_2zeta*(pc(6,2))

pc(16,9) = PA(1)*pc(7,9) + inv_2zeta*(pc(7,4))

pc(17,9) = PA(2)*pc(7,9)

pc(18,9) = PA(1)*pc(5,9) + inv_2zeta*(2*pc(2,9) + pc(5,4))

pc(19,9) = PA(2)*pc(6,9) + inv_2zeta*(2*pc(3,9))

pc(20,9) = PA(3)*pc(7,9) + inv_2zeta*(2*pc(4,9) + pc(7,2))

pc(1,10) = PB(3)*pc(1,3)

pc(2,10) = PA(1)*pc(1,10)

pc(3,10) = PA(2)*pc(1,10) + inv_2zeta*(pc(1,4))

pc(4,10) = PA(3)*pc(1,10) + inv_2zeta*(pc(1,3))

pc(5,10) = PA(1)*pc(2,10) + inv_2zeta*(pc(1,10))

pc(6,10) = PA(2)*pc(3,10) + inv_2zeta*(pc(1,10) + pc(3,4))

pc(7,10) = PA(3)*pc(4,10) + inv_2zeta*(pc(1,10) + pc(4,3))

pc(8,10) = PA(2)*pc(2,10) + inv_2zeta*(pc(2,4))

pc(9,10) = PA(3)*pc(2,10) + inv_2zeta*(pc(2,3))

pc(10,10) = PA(3)*pc(3,10) + inv_2zeta*(pc(3,3))

pc(11,10) = PA(3)*pc(8,10) + inv_2zeta*(pc(8,3))

pc(12,10) = PA(2)*pc(5,10) + inv_2zeta*(pc(5,4))

pc(13,10) = PA(3)*pc(5,10) + inv_2zeta*(pc(5,3))

pc(14,10) = PA(1)*pc(6,10)

pc(15,10) = PA(3)*pc(6,10) + inv_2zeta*(pc(6,3))

pc(16,10) = PA(1)*pc(7,10)

pc(17,10) = PA(2)*pc(7,10) + inv_2zeta*(pc(7,4))

pc(18,10) = PA(1)*pc(5,10) + inv_2zeta*(2*pc(2,10))

pc(19,10) = PA(2)*pc(6,10) + inv_2zeta*(2*pc(3,10) + pc(6,4))

pc(20,10) = PA(3)*pc(7,10) + inv_2zeta*(2*pc(4,10) + pc(7,3))

pc(1,11) = PB(3)*pc(1,8)

pc(2,11) = PA(1)*pc(1,11) + inv_2zeta*(pc(1,10))

pc(3,11) = PA(2)*pc(1,11) + inv_2zeta*(pc(1,9))

pc(4,11) = PA(3)*pc(1,11) + inv_2zeta*(pc(1,8))

pc(5,11) = PA(1)*pc(2,11) + inv_2zeta*(pc(1,11) + pc(2,10))

pc(6,11) = PA(2)*pc(3,11) + inv_2zeta*(pc(1,11) + pc(3,9))

pc(7,11) = PA(3)*pc(4,11) + inv_2zeta*(pc(1,11) + pc(4,8))

pc(8,11) = PA(2)*pc(2,11) + inv_2zeta*(pc(2,9))

pc(9,11) = PA(3)*pc(2,11) + inv_2zeta*(pc(2,8))

pc(10,11) = PA(3)*pc(3,11) + inv_2zeta*(pc(3,8))

pc(11,11) = PA(3)*pc(8,11) + inv_2zeta*(pc(8,8))

pc(12,11) = PA(2)*pc(5,11) + inv_2zeta*(pc(5,9))

pc(13,11) = PA(3)*pc(5,11) + inv_2zeta*(pc(5,8))

pc(14,11) = PA(1)*pc(6,11) + inv_2zeta*(pc(6,10))

pc(15,11) = PA(3)*pc(6,11) + inv_2zeta*(pc(6,8))

pc(16,11) = PA(1)*pc(7,11) + inv_2zeta*(pc(7,10))

pc(17,11) = PA(2)*pc(7,11) + inv_2zeta*(pc(7,9))

pc(18,11) = PA(1)*pc(5,11) + inv_2zeta*(2*pc(2,11) + pc(5,10))

pc(19,11) = PA(2)*pc(6,11) + inv_2zeta*(2*pc(3,11) + pc(6,9))

pc(20,11) = PA(3)*pc(7,11) + inv_2zeta*(2*pc(4,11) + pc(7,8))

pc(1,12) = PB(2)*pc(1,5)

pc(2,12) = PA(1)*pc(1,12) + inv_2zeta*(2*pc(1,8))

pc(3,12) = PA(2)*pc(1,12) + inv_2zeta*(pc(1,5))

pc(4,12) = PA(3)*pc(1,12)

pc(5,12) = PA(1)*pc(2,12) + inv_2zeta*(pc(1,12) + 2*pc(2,8))

pc(6,12) = PA(2)*pc(3,12) + inv_2zeta*(pc(1,12) + pc(3,5))

pc(7,12) = PA(3)*pc(4,12) + inv_2zeta*(pc(1,12))

pc(8,12) = PA(2)*pc(2,12) + inv_2zeta*(pc(2,5))

pc(9,12) = PA(3)*pc(2,12)

pc(10,12) = PA(3)*pc(3,12)

pc(11,12) = PA(3)*pc(8,12)

pc(12,12) = PA(2)*pc(5,12) + inv_2zeta*(pc(5,5))

pc(13,12) = PA(3)*pc(5,12)

pc(14,12) = PA(1)*pc(6,12) + inv_2zeta*(2*pc(6,8))

pc(15,12) = PA(3)*pc(6,12)

pc(16,12) = PA(1)*pc(7,12) + inv_2zeta*(2*pc(7,8))

pc(17,12) = PA(2)*pc(7,12) + inv_2zeta*(pc(7,5))

pc(18,12) = PA(1)*pc(5,12) + inv_2zeta*(2*pc(2,12) + 2*pc(5,8))

pc(19,12) = PA(2)*pc(6,12) + inv_2zeta*(2*pc(3,12) + pc(6,5))

pc(20,12) = PA(3)*pc(7,12) + inv_2zeta*(2*pc(4,12))

pc(1,13) = PB(3)*pc(1,5)

pc(2,13) = PA(1)*pc(1,13) + inv_2zeta*(2*pc(1,9))

pc(3,13) = PA(2)*pc(1,13)

pc(4,13) = PA(3)*pc(1,13) + inv_2zeta*(pc(1,5))

pc(5,13) = PA(1)*pc(2,13) + inv_2zeta*(pc(1,13) + 2*pc(2,9))

pc(6,13) = PA(2)*pc(3,13) + inv_2zeta*(pc(1,13))

pc(7,13) = PA(3)*pc(4,13) + inv_2zeta*(pc(1,13) + pc(4,5))

pc(8,13) = PA(2)*pc(2,13)

pc(9,13) = PA(3)*pc(2,13) + inv_2zeta*(pc(2,5))

pc(10,13) = PA(3)*pc(3,13) + inv_2zeta*(pc(3,5))

pc(11,13) = PA(3)*pc(8,13) + inv_2zeta*(pc(8,5))

pc(12,13) = PA(2)*pc(5,13)

pc(13,13) = PA(3)*pc(5,13) + inv_2zeta*(pc(5,5))

pc(14,13) = PA(1)*pc(6,13) + inv_2zeta*(2*pc(6,9))

pc(15,13) = PA(3)*pc(6,13) + inv_2zeta*(pc(6,5))

pc(16,13) = PA(1)*pc(7,13) + inv_2zeta*(2*pc(7,9))

pc(17,13) = PA(2)*pc(7,13)

pc(18,13) = PA(1)*pc(5,13) + inv_2zeta*(2*pc(2,13) + 2*pc(5,9))

pc(19,13) = PA(2)*pc(6,13) + inv_2zeta*(2*pc(3,13))

pc(20,13) = PA(3)*pc(7,13) + inv_2zeta*(2*pc(4,13) + pc(7,5))

pc(1,14) = PB(1)*pc(1,6)

pc(2,14) = PA(1)*pc(1,14) + inv_2zeta*(pc(1,6))

pc(3,14) = PA(2)*pc(1,14) + inv_2zeta*(2*pc(1,8))

pc(4,14) = PA(3)*pc(1,14)

pc(5,14) = PA(1)*pc(2,14) + inv_2zeta*(pc(1,14) + pc(2,6))

pc(6,14) = PA(2)*pc(3,14) + inv_2zeta*(pc(1,14) + 2*pc(3,8))

pc(7,14) = PA(3)*pc(4,14) + inv_2zeta*(pc(1,14))

pc(8,14) = PA(2)*pc(2,14) + inv_2zeta*(2*pc(2,8))

pc(9,14) = PA(3)*pc(2,14)

pc(10,14) = PA(3)*pc(3,14)

pc(11,14) = PA(3)*pc(8,14)

pc(12,14) = PA(2)*pc(5,14) + inv_2zeta*(2*pc(5,8))

pc(13,14) = PA(3)*pc(5,14)

pc(14,14) = PA(1)*pc(6,14) + inv_2zeta*(pc(6,6))

pc(15,14) = PA(3)*pc(6,14)

pc(16,14) = PA(1)*pc(7,14) + inv_2zeta*(pc(7,6))

pc(17,14) = PA(2)*pc(7,14) + inv_2zeta*(2*pc(7,8))

pc(18,14) = PA(1)*pc(5,14) + inv_2zeta*(2*pc(2,14) + pc(5,6))

pc(19,14) = PA(2)*pc(6,14) + inv_2zeta*(2*pc(3,14) + 2*pc(6,8))

pc(20,14) = PA(3)*pc(7,14) + inv_2zeta*(2*pc(4,14))

pc(1,15) = PB(3)*pc(1,6)

pc(2,15) = PA(1)*pc(1,15)

pc(3,15) = PA(2)*pc(1,15) + inv_2zeta*(2*pc(1,10))

pc(4,15) = PA(3)*pc(1,15) + inv_2zeta*(pc(1,6))

pc(5,15) = PA(1)*pc(2,15) + inv_2zeta*(pc(1,15))

pc(6,15) = PA(2)*pc(3,15) + inv_2zeta*(pc(1,15) + 2*pc(3,10))

pc(7,15) = PA(3)*pc(4,15) + inv_2zeta*(pc(1,15) + pc(4,6))

pc(8,15) = PA(2)*pc(2,15) + inv_2zeta*(2*pc(2,10))

pc(9,15) = PA(3)*pc(2,15) + inv_2zeta*(pc(2,6))

pc(10,15) = PA(3)*pc(3,15) + inv_2zeta*(pc(3,6))

pc(11,15) = PA(3)*pc(8,15) + inv_2zeta*(pc(8,6))

pc(12,15) = PA(2)*pc(5,15) + inv_2zeta*(2*pc(5,10))

pc(13,15) = PA(3)*pc(5,15) + inv_2zeta*(pc(5,6))

pc(14,15) = PA(1)*pc(6,15)

pc(15,15) = PA(3)*pc(6,15) + inv_2zeta*(pc(6,6))

pc(16,15) = PA(1)*pc(7,15)

pc(17,15) = PA(2)*pc(7,15) + inv_2zeta*(2*pc(7,10))

pc(18,15) = PA(1)*pc(5,15) + inv_2zeta*(2*pc(2,15))

pc(19,15) = PA(2)*pc(6,15) + inv_2zeta*(2*pc(3,15) + 2*pc(6,10))

pc(20,15) = PA(3)*pc(7,15) + inv_2zeta*(2*pc(4,15) + pc(7,6))

pc(1,16) = PB(1)*pc(1,7)

pc(2,16) = PA(1)*pc(1,16) + inv_2zeta*(pc(1,7))

pc(3,16) = PA(2)*pc(1,16)

pc(4,16) = PA(3)*pc(1,16) + inv_2zeta*(2*pc(1,9))

pc(5,16) = PA(1)*pc(2,16) + inv_2zeta*(pc(1,16) + pc(2,7))

pc(6,16) = PA(2)*pc(3,16) + inv_2zeta*(pc(1,16))

pc(7,16) = PA(3)*pc(4,16) + inv_2zeta*(pc(1,16) + 2*pc(4,9))

pc(8,16) = PA(2)*pc(2,16)

pc(9,16) = PA(3)*pc(2,16) + inv_2zeta*(2*pc(2,9))

pc(10,16) = PA(3)*pc(3,16) + inv_2zeta*(2*pc(3,9))

pc(11,16) = PA(3)*pc(8,16) + inv_2zeta*(2*pc(8,9))

pc(12,16) = PA(2)*pc(5,16)

pc(13,16) = PA(3)*pc(5,16) + inv_2zeta*(2*pc(5,9))

pc(14,16) = PA(1)*pc(6,16) + inv_2zeta*(pc(6,7))

pc(15,16) = PA(3)*pc(6,16) + inv_2zeta*(2*pc(6,9))

pc(16,16) = PA(1)*pc(7,16) + inv_2zeta*(pc(7,7))

pc(17,16) = PA(2)*pc(7,16)

pc(18,16) = PA(1)*pc(5,16) + inv_2zeta*(2*pc(2,16) + pc(5,7))

pc(19,16) = PA(2)*pc(6,16) + inv_2zeta*(2*pc(3,16))

pc(20,16) = PA(3)*pc(7,16) + inv_2zeta*(2*pc(4,16) + 2*pc(7,9))

pc(1,17) = PB(2)*pc(1,7)

pc(2,17) = PA(1)*pc(1,17)

pc(3,17) = PA(2)*pc(1,17) + inv_2zeta*(pc(1,7))

pc(4,17) = PA(3)*pc(1,17) + inv_2zeta*(2*pc(1,10))

pc(5,17) = PA(1)*pc(2,17) + inv_2zeta*(pc(1,17))

pc(6,17) = PA(2)*pc(3,17) + inv_2zeta*(pc(1,17) + pc(3,7))

pc(7,17) = PA(3)*pc(4,17) + inv_2zeta*(pc(1,17) + 2*pc(4,10))

pc(8,17) = PA(2)*pc(2,17) + inv_2zeta*(pc(2,7))

pc(9,17) = PA(3)*pc(2,17) + inv_2zeta*(2*pc(2,10))

pc(10,17) = PA(3)*pc(3,17) + inv_2zeta*(2*pc(3,10))

pc(11,17) = PA(3)*pc(8,17) + inv_2zeta*(2*pc(8,10))

pc(12,17) = PA(2)*pc(5,17) + inv_2zeta*(pc(5,7))

pc(13,17) = PA(3)*pc(5,17) + inv_2zeta*(2*pc(5,10))

pc(14,17) = PA(1)*pc(6,17)

pc(15,17) = PA(3)*pc(6,17) + inv_2zeta*(2*pc(6,10))

pc(16,17) = PA(1)*pc(7,17)

pc(17,17) = PA(2)*pc(7,17) + inv_2zeta*(pc(7,7))

pc(18,17) = PA(1)*pc(5,17) + inv_2zeta*(2*pc(2,17))

pc(19,17) = PA(2)*pc(6,17) + inv_2zeta*(2*pc(3,17) + pc(6,7))

pc(20,17) = PA(3)*pc(7,17) + inv_2zeta*(2*pc(4,17) + 2*pc(7,10))

pc(1,18) = PB(1)*pc(1,5) + inv_2zeta*(2*pc(1,2))

pc(2,18) = PA(1)*pc(1,18) + inv_2zeta*(3*pc(1,5))

pc(3,18) = PA(2)*pc(1,18)

pc(4,18) = PA(3)*pc(1,18)

pc(5,18) = PA(1)*pc(2,18) + inv_2zeta*(pc(1,18) + 3*pc(2,5))

pc(6,18) = PA(2)*pc(3,18) + inv_2zeta*(pc(1,18))

pc(7,18) = PA(3)*pc(4,18) + inv_2zeta*(pc(1,18))

pc(8,18) = PA(2)*pc(2,18)

pc(9,18) = PA(3)*pc(2,18)

pc(10,18) = PA(3)*pc(3,18)

pc(11,18) = PA(3)*pc(8,18)

pc(12,18) = PA(2)*pc(5,18)

pc(13,18) = PA(3)*pc(5,18)

pc(14,18) = PA(1)*pc(6,18) + inv_2zeta*(3*pc(6,5))

pc(15,18) = PA(3)*pc(6,18)

pc(16,18) = PA(1)*pc(7,18) + inv_2zeta*(3*pc(7,5))

pc(17,18) = PA(2)*pc(7,18)

pc(18,18) = PA(1)*pc(5,18) + inv_2zeta*(2*pc(2,18) + 3*pc(5,5))

pc(19,18) = PA(2)*pc(6,18) + inv_2zeta*(2*pc(3,18))

pc(20,18) = PA(3)*pc(7,18) + inv_2zeta*(2*pc(4,18))

pc(1,19) = PB(2)*pc(1,6) + inv_2zeta*(2*pc(1,3))

pc(2,19) = PA(1)*pc(1,19)

pc(3,19) = PA(2)*pc(1,19) + inv_2zeta*(3*pc(1,6))

pc(4,19) = PA(3)*pc(1,19)

pc(5,19) = PA(1)*pc(2,19) + inv_2zeta*(pc(1,19))

pc(6,19) = PA(2)*pc(3,19) + inv_2zeta*(pc(1,19) + 3*pc(3,6))

pc(7,19) = PA(3)*pc(4,19) + inv_2zeta*(pc(1,19))

pc(8,19) = PA(2)*pc(2,19) + inv_2zeta*(3*pc(2,6))

pc(9,19) = PA(3)*pc(2,19)

pc(10,19) = PA(3)*pc(3,19)

pc(11,19) = PA(3)*pc(8,19)

pc(12,19) = PA(2)*pc(5,19) + inv_2zeta*(3*pc(5,6))

pc(13,19) = PA(3)*pc(5,19)

pc(14,19) = PA(1)*pc(6,19)

pc(15,19) = PA(3)*pc(6,19)

pc(16,19) = PA(1)*pc(7,19)

pc(17,19) = PA(2)*pc(7,19) + inv_2zeta*(3*pc(7,6))

pc(18,19) = PA(1)*pc(5,19) + inv_2zeta*(2*pc(2,19))

pc(19,19) = PA(2)*pc(6,19) + inv_2zeta*(2*pc(3,19) + 3*pc(6,6))

pc(20,19) = PA(3)*pc(7,19) + inv_2zeta*(2*pc(4,19))

pc(1,20) = PB(3)*pc(1,7) + inv_2zeta*(2*pc(1,4))

pc(2,20) = PA(1)*pc(1,20)

pc(3,20) = PA(2)*pc(1,20)

pc(4,20) = PA(3)*pc(1,20) + inv_2zeta*(3*pc(1,7))

pc(5,20) = PA(1)*pc(2,20) + inv_2zeta*(pc(1,20))

pc(6,20) = PA(2)*pc(3,20) + inv_2zeta*(pc(1,20))

pc(7,20) = PA(3)*pc(4,20) + inv_2zeta*(pc(1,20) + 3*pc(4,7))

pc(8,20) = PA(2)*pc(2,20)

pc(9,20) = PA(3)*pc(2,20) + inv_2zeta*(3*pc(2,7))

pc(10,20) = PA(3)*pc(3,20) + inv_2zeta*(3*pc(3,7))

pc(11,20) = PA(3)*pc(8,20) + inv_2zeta*(3*pc(8,7))

pc(12,20) = PA(2)*pc(5,20)

pc(13,20) = PA(3)*pc(5,20) + inv_2zeta*(3*pc(5,7))

pc(14,20) = PA(1)*pc(6,20)

pc(15,20) = PA(3)*pc(6,20) + inv_2zeta*(3*pc(6,7))

pc(16,20) = PA(1)*pc(7,20)

pc(17,20) = PA(2)*pc(7,20)

pc(18,20) = PA(1)*pc(5,20) + inv_2zeta*(2*pc(2,20))

pc(19,20) = PA(2)*pc(6,20) + inv_2zeta*(2*pc(3,20))

pc(20,20) = PA(3)*pc(7,20) + inv_2zeta*(2*pc(4,20) + 3*pc(7,7))

sh(1,1) = pc(1,1)

sh(2,1) = pc(2,1)

sh(3,1) = pc(3,1)

sh(4,1) = pc(4,1)

sh(5,1) = pc(8,1)

sh(6,1) = pc(9,1)

sh(7,1) = pc(10,1)

sh(8,1) = pc(5,1) - pc(6,1)

sh(9,1) = 2*pc(7,1) - pc(5,1) - pc(6,1)

sh(10,1) = pc(11,1)

sh(11,1) = pc(13,1) - pc(15,1)

sh(12,1) = pc(18,1) - 3*pc(14,1)

sh(13,1) = 3*pc(12,1) - pc(19,1)

sh(14,1) = 2*pc(20,1) - 3*pc(13,1) - 3*pc(15,1)

sh(15,1) = 4*pc(16,1) - pc(18,1) - pc(14,1)

sh(16,1) = 4*pc(17,1) - pc(12,1) - pc(19,1)

sh(1,2) = pc(1,2)

sh(2,2) = pc(2,2)

sh(3,2) = pc(3,2)

sh(4,2) = pc(4,2)

sh(5,2) = pc(8,2)

sh(6,2) = pc(9,2)

sh(7,2) = pc(10,2)

sh(8,2) = pc(5,2) - pc(6,2)

sh(9,2) = 2*pc(7,2) - pc(5,2) - pc(6,2)

sh(10,2) = pc(11,2)

sh(11,2) = pc(13,2) - pc(15,2)

sh(12,2) = pc(18,2) - 3*pc(14,2)

sh(13,2) = 3*pc(12,2) - pc(19,2)

sh(14,2) = 2*pc(20,2) - 3*pc(13,2) - 3*pc(15,2)

sh(15,2) = 4*pc(16,2) - pc(18,2) - pc(14,2)

sh(16,2) = 4*pc(17,2) - pc(12,2) - pc(19,2)

sh(1,3) = pc(1,3)

sh(2,3) = pc(2,3)

sh(3,3) = pc(3,3)

sh(4,3) = pc(4,3)

sh(5,3) = pc(8,3)

sh(6,3) = pc(9,3)

sh(7,3) = pc(10,3)

sh(8,3) = pc(5,3) - pc(6,3)

sh(9,3) = 2*pc(7,3) - pc(5,3) - pc(6,3)

sh(10,3) = pc(11,3)

sh(11,3) = pc(13,3) - pc(15,3)

sh(12,3) = pc(18,3) - 3*pc(14,3)

sh(13,3) = 3*pc(12,3) - pc(19,3)

sh(14,3) = 2*pc(20,3) - 3*pc(13,3) - 3*pc(15,3)

sh(15,3) = 4*pc(16,3) - pc(18,3) - pc(14,3)

sh(16,3) = 4*pc(17,3) - pc(12,3) - pc(19,3)

sh(1,4) = pc(1,4)

sh(2,4) = pc(2,4)

sh(3,4) = pc(3,4)

sh(4,4) = pc(4,4)

sh(5,4) = pc(8,4)

sh(6,4) = pc(9,4)

sh(7,4) = pc(10,4)

sh(8,4) = pc(5,4) - pc(6,4)

sh(9,4) = 2*pc(7,4) - pc(5,4) - pc(6,4)

sh(10,4) = pc(11,4)

sh(11,4) = pc(13,4) - pc(15,4)

sh(12,4) = pc(18,4) - 3*pc(14,4)

sh(13,4) = 3*pc(12,4) - pc(19,4)

sh(14,4) = 2*pc(20,4) - 3*pc(13,4) - 3*pc(15,4)

sh(15,4) = 4*pc(16,4) - pc(18,4) - pc(14,4)

sh(16,4) = 4*pc(17,4) - pc(12,4) - pc(19,4)

sh(1,5) = pc(1,8)

sh(2,5) = pc(2,8)

sh(3,5) = pc(3,8)

sh(4,5) = pc(4,8)

sh(5,5) = pc(8,8)

sh(6,5) = pc(9,8)

sh(7,5) = pc(10,8)

sh(8,5) = pc(5,8) - pc(6,8)

sh(9,5) = 2*pc(7,8) - pc(5,8) - pc(6,8)

sh(10,5) = pc(11,8)

sh(11,5) = pc(13,8) - pc(15,8)

sh(12,5) = pc(18,8) - 3*pc(14,8)

sh(13,5) = 3*pc(12,8) - pc(19,8)

sh(14,5) = 2*pc(20,8) - 3*pc(13,8) - 3*pc(15,8)

sh(15,5) = 4*pc(16,8) - pc(18,8) - pc(14,8)

sh(16,5) = 4*pc(17,8) - pc(12,8) - pc(19,8)

sh(1,6) = pc(1,9)

sh(2,6) = pc(2,9)

sh(3,6) = pc(3,9)

sh(4,6) = pc(4,9)

sh(5,6) = pc(8,9)

sh(6,6) = pc(9,9)

sh(7,6) = pc(10,9)

sh(8,6) = pc(5,9) - pc(6,9)

sh(9,6) = 2*pc(7,9) - pc(5,9) - pc(6,9)

sh(10,6) = pc(11,9)

sh(11,6) = pc(13,9) - pc(15,9)

sh(12,6) = pc(18,9) - 3*pc(14,9)

sh(13,6) = 3*pc(12,9) - pc(19,9)

sh(14,6) = 2*pc(20,9) - 3*pc(13,9) - 3*pc(15,9)

sh(15,6) = 4*pc(16,9) - pc(18,9) - pc(14,9)

sh(16,6) = 4*pc(17,9) - pc(12,9) - pc(19,9)

sh(1,7) = pc(1,10)

sh(2,7) = pc(2,10)

sh(3,7) = pc(3,10)

sh(4,7) = pc(4,10)

sh(5,7) = pc(8,10)

sh(6,7) = pc(9,10)

sh(7,7) = pc(10,10)

sh(8,7) = pc(5,10) - pc(6,10)

sh(9,7) = 2*pc(7,10) - pc(5,10) - pc(6,10)

sh(10,7) = pc(11,10)

sh(11,7) = pc(13,10) - pc(15,10)

sh(12,7) = pc(18,10) - 3*pc(14,10)

sh(13,7) = 3*pc(12,10) - pc(19,10)

sh(14,7) = 2*pc(20,10) - 3*pc(13,10) - 3*pc(15,10)

sh(15,7) = 4*pc(16,10) - pc(18,10) - pc(14,10)

sh(16,7) = 4*pc(17,10) - pc(12,10) - pc(19,10)

sh(1,8) = pc(1,5) - pc(1,6)

sh(2,8) = pc(2,5) - pc(2,6)

sh(3,8) = pc(3,5) - pc(3,6)

sh(4,8) = pc(4,5) - pc(4,6)

sh(5,8) = pc(8,5) - pc(8,6)

sh(6,8) = pc(9,5) - pc(9,6)

sh(7,8) = pc(10,5) - pc(10,6)

sh(8,8) = pc(5,5) - pc(5,6) - pc(6,5) + pc(6,6)

sh(9,8) = 2*pc(7,5) - 2*pc(7,6) - pc(5,5) + pc(5,6) - pc(6,5) + pc(6,6)

sh(10,8) = pc(11,5) - pc(11,6)

sh(11,8) = pc(13,5) - pc(13,6) - pc(15,5) + pc(15,6)

sh(12,8) = pc(18,5) - pc(18,6) - 3*pc(14,5) + 3*pc(14,6)

sh(13,8) = 3*pc(12,5) - 3*pc(12,6) - pc(19,5) + pc(19,6)

sh(14,8) = 2*pc(20,5) - 2*pc(20,6) - 3*pc(13,5) + 3*pc(13,6) - 3*pc(15,5) + 3*&
&pc(15,6)

sh(15,8) = 4*pc(16,5) - 4*pc(16,6) - pc(18,5) + pc(18,6) - pc(14,5) + pc(14,6)

sh(16,8) = 4*pc(17,5) - 4*pc(17,6) - pc(12,5) + pc(12,6) - pc(19,5) + pc(19,6)

sh(1,9) = 2*pc(1,7) - pc(1,5) - pc(1,6)

sh(2,9) = 2*pc(2,7) - pc(2,5) - pc(2,6)

sh(3,9) = 2*pc(3,7) - pc(3,5) - pc(3,6)

sh(4,9) = 2*pc(4,7) - pc(4,5) - pc(4,6)

sh(5,9) = 2*pc(8,7) - pc(8,5) - pc(8,6)

sh(6,9) = 2*pc(9,7) - pc(9,5) - pc(9,6)

sh(7,9) = 2*pc(10,7) - pc(10,5) - pc(10,6)

sh(8,9) = 2*pc(5,7) - pc(5,5) - pc(5,6) - 2*pc(6,7) + pc(6,5) + pc(6,6)

sh(9,9) = 4*pc(7,7) - 2*pc(7,5) - 2*pc(7,6) - 2*pc(5,7) + pc(5,5) + pc(5,6) - &
&2*pc(6,7) + pc(6,5) + pc(6,6)

sh(10,9) = 2*pc(11,7) - pc(11,5) - pc(11,6)

sh(11,9) = 2*pc(13,7) - pc(13,5) - pc(13,6) - 2*pc(15,7) + pc(15,5) + pc(15,6)

sh(12,9) = 2*pc(18,7) - pc(18,5) - pc(18,6) - 6*pc(14,7) + 3*pc(14,5) + 3*pc(1&
&4,6)

sh(13,9) = 6*pc(12,7) - 3*pc(12,5) - 3*pc(12,6) - 2*pc(19,7) + pc(19,5) + pc(1&
&9,6)

sh(14,9) = 4*pc(20,7) - 2*pc(20,5) - 2*pc(20,6) - 6*pc(13,7) + 3*pc(13,5) + 3*&
&pc(13,6) - 6*pc(15,7) + 3*pc(15,5) + 3*pc(15,6)

sh(15,9) = 8*pc(16,7) - 4*pc(16,5) - 4*pc(16,6) - 2*pc(18,7) + pc(18,5) + pc(1&
&8,6) - 2*pc(14,7) + pc(14,5) + pc(14,6)

sh(16,9) = 8*pc(17,7) - 4*pc(17,5) - 4*pc(17,6) - 2*pc(12,7) + pc(12,5) + pc(1&
&2,6) - 2*pc(19,7) + pc(19,5) + pc(19,6)

sh(1,10) = pc(1,11)

sh(2,10) = pc(2,11)

sh(3,10) = pc(3,11)

sh(4,10) = pc(4,11)

sh(5,10) = pc(8,11)

sh(6,10) = pc(9,11)

sh(7,10) = pc(10,11)

sh(8,10) = pc(5,11) - pc(6,11)

sh(9,10) = 2*pc(7,11) - pc(5,11) - pc(6,11)

sh(10,10) = pc(11,11)

sh(11,10) = pc(13,11) - pc(15,11)

sh(12,10) = pc(18,11) - 3*pc(14,11)

sh(13,10) = 3*pc(12,11) - pc(19,11)

sh(14,10) = 2*pc(20,11) - 3*pc(13,11) - 3*pc(15,11)

sh(15,10) = 4*pc(16,11) - pc(18,11) - pc(14,11)

sh(16,10) = 4*pc(17,11) - pc(12,11) - pc(19,11)

sh(1,11) = pc(1,13) - pc(1,15)

sh(2,11) = pc(2,13) - pc(2,15)

sh(3,11) = pc(3,13) - pc(3,15)

sh(4,11) = pc(4,13) - pc(4,15)

sh(5,11) = pc(8,13) - pc(8,15)

sh(6,11) = pc(9,13) - pc(9,15)

sh(7,11) = pc(10,13) - pc(10,15)

sh(8,11) = pc(5,13) - pc(5,15) - pc(6,13) + pc(6,15)

sh(9,11) = 2*pc(7,13) - 2*pc(7,15) - pc(5,13) + pc(5,15) - pc(6,13) + pc(6,15)

sh(10,11) = pc(11,13) - pc(11,15)

sh(11,11) = pc(13,13) - pc(13,15) - pc(15,13) + pc(15,15)

sh(12,11) = pc(18,13) - pc(18,15) - 3*pc(14,13) + 3*pc(14,15)

sh(13,11) = 3*pc(12,13) - 3*pc(12,15) - pc(19,13) + pc(19,15)

sh(14,11) = 2*pc(20,13) - 2*pc(20,15) - 3*pc(13,13) + 3*pc(13,15) - 3*pc(15,13&
&) + 3*pc(15,15)

sh(15,11) = 4*pc(16,13) - 4*pc(16,15) - pc(18,13) + pc(18,15) - pc(14,13) + pc&
&(14,15)

sh(16,11) = 4*pc(17,13) - 4*pc(17,15) - pc(12,13) + pc(12,15) - pc(19,13) + pc&
&(19,15)

sh(1,12) = pc(1,18) - 3*pc(1,14)

sh(2,12) = pc(2,18) - 3*pc(2,14)

sh(3,12) = pc(3,18) - 3*pc(3,14)

sh(4,12) = pc(4,18) - 3*pc(4,14)

sh(5,12) = pc(8,18) - 3*pc(8,14)

sh(6,12) = pc(9,18) - 3*pc(9,14)

sh(7,12) = pc(10,18) - 3*pc(10,14)

sh(8,12) = pc(5,18) - 3*pc(5,14) - pc(6,18) + 3*pc(6,14)

sh(9,12) = 2*pc(7,18) - 6*pc(7,14) - pc(5,18) + 3*pc(5,14) - pc(6,18) + 3*pc(6&
&,14)

sh(10,12) = pc(11,18) - 3*pc(11,14)

sh(11,12) = pc(13,18) - 3*pc(13,14) - pc(15,18) + 3*pc(15,14)

sh(12,12) = pc(18,18) - 3*pc(18,14) - 3*pc(14,18) + 9*pc(14,14)

sh(13,12) = 3*pc(12,18) - 9*pc(12,14) - pc(19,18) + 3*pc(19,14)

sh(14,12) = 2*pc(20,18) - 6*pc(20,14) - 3*pc(13,18) + 9*pc(13,14) - 3*pc(15,18&
&) + 9*pc(15,14)

sh(15,12) = 4*pc(16,18) - 12*pc(16,14) - pc(18,18) + 3*pc(18,14) - pc(14,18) +&
& 3*pc(14,14)

sh(16,12) = 4*pc(17,18) - 12*pc(17,14) - pc(12,18) + 3*pc(12,14) - pc(19,18) +&
& 3*pc(19,14)

sh(1,13) = 3*pc(1,12) - pc(1,19)

sh(2,13) = 3*pc(2,12) - pc(2,19)

sh(3,13) = 3*pc(3,12) - pc(3,19)

sh(4,13) = 3*pc(4,12) - pc(4,19)

sh(5,13) = 3*pc(8,12) - pc(8,19)

sh(6,13) = 3*pc(9,12) - pc(9,19)

sh(7,13) = 3*pc(10,12) - pc(10,19)

sh(8,13) = 3*pc(5,12) - pc(5,19) - 3*pc(6,12) + pc(6,19)

sh(9,13) = 6*pc(7,12) - 2*pc(7,19) - 3*pc(5,12) + pc(5,19) - 3*pc(6,12) + pc(6&
&,19)

sh(10,13) = 3*pc(11,12) - pc(11,19)

sh(11,13) = 3*pc(13,12) - pc(13,19) - 3*pc(15,12) + pc(15,19)

sh(12,13) = 3*pc(18,12) - pc(18,19) - 9*pc(14,12) + 3*pc(14,19)

sh(13,13) = 9*pc(12,12) - 3*pc(12,19) - 3*pc(19,12) + pc(19,19)

sh(14,13) = 6*pc(20,12) - 2*pc(20,19) - 9*pc(13,12) + 3*pc(13,19) - 9*pc(15,12&
&) + 3*pc(15,19)

sh(15,13) = 12*pc(16,12) - 4*pc(16,19) - 3*pc(18,12) + pc(18,19) - 3*pc(14,12)&
& + pc(14,19)

sh(16,13) = 12*pc(17,12) - 4*pc(17,19) - 3*pc(12,12) + pc(12,19) - 3*pc(19,12)&
& + pc(19,19)

sh(1,14) = 2*pc(1,20) - 3*pc(1,13) - 3*pc(1,15)

sh(2,14) = 2*pc(2,20) - 3*pc(2,13) - 3*pc(2,15)

sh(3,14) = 2*pc(3,20) - 3*pc(3,13) - 3*pc(3,15)

sh(4,14) = 2*pc(4,20) - 3*pc(4,13) - 3*pc(4,15)

sh(5,14) = 2*pc(8,20) - 3*pc(8,13) - 3*pc(8,15)

sh(6,14) = 2*pc(9,20) - 3*pc(9,13) - 3*pc(9,15)

sh(7,14) = 2*pc(10,20) - 3*pc(10,13) - 3*pc(10,15)

sh(8,14) = 2*pc(5,20) - 3*pc(5,13) - 3*pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*&
&pc(6,15)

sh(9,14) = 4*pc(7,20) - 6*pc(7,13) - 6*pc(7,15) - 2*pc(5,20) + 3*pc(5,13) + 3*&
&pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*pc(6,15)

sh(10,14) = 2*pc(11,20) - 3*pc(11,13) - 3*pc(11,15)

sh(11,14) = 2*pc(13,20) - 3*pc(13,13) - 3*pc(13,15) - 2*pc(15,20) + 3*pc(15,13&
&) + 3*pc(15,15)

sh(12,14) = 2*pc(18,20) - 3*pc(18,13) - 3*pc(18,15) - 6*pc(14,20) + 9*pc(14,13&
&) + 9*pc(14,15)

sh(13,14) = 6*pc(12,20) - 9*pc(12,13) - 9*pc(12,15) - 2*pc(19,20) + 3*pc(19,13&
&) + 3*pc(19,15)

sh(14,14) = 4*pc(20,20) - 6*pc(20,13) - 6*pc(20,15) - 6*pc(13,20) + 9*pc(13,13&
&) + 9*pc(13,15) - 6*pc(15,20) + 9*pc(15,13) + 9*pc(15,15)

sh(15,14) = 8*pc(16,20) - 12*pc(16,13) - 12*pc(16,15) - 2*pc(18,20) + 3*pc(18,&
&13) + 3*pc(18,15) - 2*pc(14,20) + 3*pc(14,13) + 3*pc(14,15)

sh(16,14) = 8*pc(17,20) - 12*pc(17,13) - 12*pc(17,15) - 2*pc(12,20) + 3*pc(12,&
&13) + 3*pc(12,15) - 2*pc(19,20) + 3*pc(19,13) + 3*pc(19,15)

sh(1,15) = 4*pc(1,16) - pc(1,18) - pc(1,14)

sh(2,15) = 4*pc(2,16) - pc(2,18) - pc(2,14)

sh(3,15) = 4*pc(3,16) - pc(3,18) - pc(3,14)

sh(4,15) = 4*pc(4,16) - pc(4,18) - pc(4,14)

sh(5,15) = 4*pc(8,16) - pc(8,18) - pc(8,14)

sh(6,15) = 4*pc(9,16) - pc(9,18) - pc(9,14)

sh(7,15) = 4*pc(10,16) - pc(10,18) - pc(10,14)

sh(8,15) = 4*pc(5,16) - pc(5,18) - pc(5,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(9,15) = 8*pc(7,16) - 2*pc(7,18) - 2*pc(7,14) - 4*pc(5,16) + pc(5,18) + pc(5&
&,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(10,15) = 4*pc(11,16) - pc(11,18) - pc(11,14)

sh(11,15) = 4*pc(13,16) - pc(13,18) - pc(13,14) - 4*pc(15,16) + pc(15,18) + pc&
&(15,14)

sh(12,15) = 4*pc(18,16) - pc(18,18) - pc(18,14) - 12*pc(14,16) + 3*pc(14,18) +&
& 3*pc(14,14)

sh(13,15) = 12*pc(12,16) - 3*pc(12,18) - 3*pc(12,14) - 4*pc(19,16) + pc(19,18)&
& + pc(19,14)

sh(14,15) = 8*pc(20,16) - 2*pc(20,18) - 2*pc(20,14) - 12*pc(13,16) + 3*pc(13,1&
&8) + 3*pc(13,14) - 12*pc(15,16) + 3*pc(15,18) + 3*pc(15,14)

sh(15,15) = 16*pc(16,16) - 4*pc(16,18) - 4*pc(16,14) - 4*pc(18,16) + pc(18,18)&
& + pc(18,14) - 4*pc(14,16) + pc(14,18) + pc(14,14)

sh(16,15) = 16*pc(17,16) - 4*pc(17,18) - 4*pc(17,14) - 4*pc(12,16) + pc(12,18)&
& + pc(12,14) - 4*pc(19,16) + pc(19,18) + pc(19,14)

sh(1,16) = 4*pc(1,17) - pc(1,12) - pc(1,19)

sh(2,16) = 4*pc(2,17) - pc(2,12) - pc(2,19)

sh(3,16) = 4*pc(3,17) - pc(3,12) - pc(3,19)

sh(4,16) = 4*pc(4,17) - pc(4,12) - pc(4,19)

sh(5,16) = 4*pc(8,17) - pc(8,12) - pc(8,19)

sh(6,16) = 4*pc(9,17) - pc(9,12) - pc(9,19)

sh(7,16) = 4*pc(10,17) - pc(10,12) - pc(10,19)

sh(8,16) = 4*pc(5,17) - pc(5,12) - pc(5,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(9,16) = 8*pc(7,17) - 2*pc(7,12) - 2*pc(7,19) - 4*pc(5,17) + pc(5,12) + pc(5&
&,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(10,16) = 4*pc(11,17) - pc(11,12) - pc(11,19)

sh(11,16) = 4*pc(13,17) - pc(13,12) - pc(13,19) - 4*pc(15,17) + pc(15,12) + pc&
&(15,19)

sh(12,16) = 4*pc(18,17) - pc(18,12) - pc(18,19) - 12*pc(14,17) + 3*pc(14,12) +&
& 3*pc(14,19)

sh(13,16) = 12*pc(12,17) - 3*pc(12,12) - 3*pc(12,19) - 4*pc(19,17) + pc(19,12)&
& + pc(19,19)

sh(14,16) = 8*pc(20,17) - 2*pc(20,12) - 2*pc(20,19) - 12*pc(13,17) + 3*pc(13,1&
&2) + 3*pc(13,19) - 12*pc(15,17) + 3*pc(15,12) + 3*pc(15,19)

sh(15,16) = 16*pc(16,17) - 4*pc(16,12) - 4*pc(16,19) - 4*pc(18,17) + pc(18,12)&
& + pc(18,19) - 4*pc(14,17) + pc(14,12) + pc(14,19)

sh(16,16) = 16*pc(17,17) - 4*pc(17,12) - 4*pc(17,19) - 4*pc(12,17) + pc(12,12)&
& + pc(12,19) - 4*pc(19,17) + pc(19,12) + pc(19,19)


   end subroutine overlap2CIntgAna

   subroutine overlap2CIntgNum(a1,a2,A,B,pc,sh,cell_size,step_size)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   real (kind=double), dimension (20,20), intent(out) :: pc
   real (kind=double), dimension (16,16), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i
   integer :: num_steps
   integer, dimension (20,3) :: triads
   integer, dimension (16,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos
   real (kind=double), dimension (3) :: xyz, xyz_sum, xyz_soln

   ! Initialize local variables.
   triads(1,:) = (/0,0,0/)
   triads(2,:) = (/1,0,0/)
   triads(3,:) = (/0,1,0/)
   triads(4,:) = (/0,0,1/)
   triads(5,:) = (/2,0,0/)
   triads(6,:) = (/0,2,0/)
   triads(7,:) = (/0,0,2/)
   triads(8,:) = (/1,1,0/)
   triads(9,:) = (/1,0,1/)
   triads(10,:) = (/0,1,1/)
   triads(11,:) = (/1,1,1/)
   triads(12,:) = (/2,1,0/)
   triads(13,:) = (/2,0,1/)
   triads(14,:) = (/1,2,0/)
   triads(15,:) = (/0,2,1/)
   triads(16,:) = (/1,0,2/)
   triads(17,:) = (/0,1,2/)
   triads(18,:) = (/3,0,0/)
   triads(19,:) = (/0,3,0/)
   triads(20,:) = (/0,0,3/)
   conversion(1,1,:) = (/1,0,0/)
   conversion(1,2,:) = (/1,1,1/)
   conversion(2,1,:) = (/1,0,0/)
   conversion(2,2,:) = (/2,1,1/)
   conversion(3,1,:) = (/1,0,0/)
   conversion(3,2,:) = (/3,1,1/)
   conversion(4,1,:) = (/1,0,0/)
   conversion(4,2,:) = (/4,1,1/)
   conversion(5,1,:) = (/1,0,0/)
   conversion(5,2,:) = (/8,1,1/)
   conversion(6,1,:) = (/1,0,0/)
   conversion(6,2,:) = (/9,1,1/)
   conversion(7,1,:) = (/1,0,0/)
   conversion(7,2,:) = (/10,1,1/)
   conversion(8,1,:) = (/1,-1,0/)
   conversion(8,2,:) = (/5,6,1/)
   conversion(9,1,:) = (/2,-1,-1/)
   conversion(9,2,:) = (/7,5,6/)
   conversion(10,1,:) = (/1,0,0/)
   conversion(10,2,:) = (/11,1,1/)
   conversion(11,1,:) = (/1,-1,0/)
   conversion(11,2,:) = (/13,15,1/)
   conversion(12,1,:) = (/1,-3,0/)
   conversion(12,2,:) = (/18,14,1/)
   conversion(13,1,:) = (/3,-1,0/)
   conversion(13,2,:) = (/12,19,1/)
   conversion(14,1,:) = (/2,-3,-3/)
   conversion(14,2,:) = (/20,13,15/)
   conversion(15,1,:) = (/4,-1,-1/)
   conversion(15,2,:) = (/16,18,14/)
   conversion(16,1,:) = (/4,-1,-1/)
   conversion(16,2,:) = (/17,12,19/)

   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   do p = 1, 20
      do q = 1, 20

         ! Assign l1 and l2 values for each gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize sum variables.
         xyz_sum(:) = 0.0d0

         do i = 0, num_steps
            xyz(:) = (start_pos + (i*step_size))
            xyz_soln(:) = step_size * &
&(xyz(:) - A(:))**l1(:)*(xyz(:) - B(:))**l2(:)*exp(-a1*(xyz(:) - A(:))**2)*exp&
&(-a2*(xyz(:) - B(:))**2)
            xyz_sum(:) = xyz_sum(:) + xyz_soln(:)
         enddo

         pc(q,p) = product(xyz_sum(:))
      enddo
    enddo

sh(1,1) = pc(1,1)

sh(2,1) = pc(2,1)

sh(3,1) = pc(3,1)

sh(4,1) = pc(4,1)

sh(5,1) = pc(8,1)

sh(6,1) = pc(9,1)

sh(7,1) = pc(10,1)

sh(8,1) = pc(5,1) - pc(6,1)

sh(9,1) = 2*pc(7,1) - pc(5,1) - pc(6,1)

sh(10,1) = pc(11,1)

sh(11,1) = pc(13,1) - pc(15,1)

sh(12,1) = pc(18,1) - 3*pc(14,1)

sh(13,1) = 3*pc(12,1) - pc(19,1)

sh(14,1) = 2*pc(20,1) - 3*pc(13,1) - 3*pc(15,1)

sh(15,1) = 4*pc(16,1) - pc(18,1) - pc(14,1)

sh(16,1) = 4*pc(17,1) - pc(12,1) - pc(19,1)

sh(1,2) = pc(1,2)

sh(2,2) = pc(2,2)

sh(3,2) = pc(3,2)

sh(4,2) = pc(4,2)

sh(5,2) = pc(8,2)

sh(6,2) = pc(9,2)

sh(7,2) = pc(10,2)

sh(8,2) = pc(5,2) - pc(6,2)

sh(9,2) = 2*pc(7,2) - pc(5,2) - pc(6,2)

sh(10,2) = pc(11,2)

sh(11,2) = pc(13,2) - pc(15,2)

sh(12,2) = pc(18,2) - 3*pc(14,2)

sh(13,2) = 3*pc(12,2) - pc(19,2)

sh(14,2) = 2*pc(20,2) - 3*pc(13,2) - 3*pc(15,2)

sh(15,2) = 4*pc(16,2) - pc(18,2) - pc(14,2)

sh(16,2) = 4*pc(17,2) - pc(12,2) - pc(19,2)

sh(1,3) = pc(1,3)

sh(2,3) = pc(2,3)

sh(3,3) = pc(3,3)

sh(4,3) = pc(4,3)

sh(5,3) = pc(8,3)

sh(6,3) = pc(9,3)

sh(7,3) = pc(10,3)

sh(8,3) = pc(5,3) - pc(6,3)

sh(9,3) = 2*pc(7,3) - pc(5,3) - pc(6,3)

sh(10,3) = pc(11,3)

sh(11,3) = pc(13,3) - pc(15,3)

sh(12,3) = pc(18,3) - 3*pc(14,3)

sh(13,3) = 3*pc(12,3) - pc(19,3)

sh(14,3) = 2*pc(20,3) - 3*pc(13,3) - 3*pc(15,3)

sh(15,3) = 4*pc(16,3) - pc(18,3) - pc(14,3)

sh(16,3) = 4*pc(17,3) - pc(12,3) - pc(19,3)

sh(1,4) = pc(1,4)

sh(2,4) = pc(2,4)

sh(3,4) = pc(3,4)

sh(4,4) = pc(4,4)

sh(5,4) = pc(8,4)

sh(6,4) = pc(9,4)

sh(7,4) = pc(10,4)

sh(8,4) = pc(5,4) - pc(6,4)

sh(9,4) = 2*pc(7,4) - pc(5,4) - pc(6,4)

sh(10,4) = pc(11,4)

sh(11,4) = pc(13,4) - pc(15,4)

sh(12,4) = pc(18,4) - 3*pc(14,4)

sh(13,4) = 3*pc(12,4) - pc(19,4)

sh(14,4) = 2*pc(20,4) - 3*pc(13,4) - 3*pc(15,4)

sh(15,4) = 4*pc(16,4) - pc(18,4) - pc(14,4)

sh(16,4) = 4*pc(17,4) - pc(12,4) - pc(19,4)

sh(1,5) = pc(1,8)

sh(2,5) = pc(2,8)

sh(3,5) = pc(3,8)

sh(4,5) = pc(4,8)

sh(5,5) = pc(8,8)

sh(6,5) = pc(9,8)

sh(7,5) = pc(10,8)

sh(8,5) = pc(5,8) - pc(6,8)

sh(9,5) = 2*pc(7,8) - pc(5,8) - pc(6,8)

sh(10,5) = pc(11,8)

sh(11,5) = pc(13,8) - pc(15,8)

sh(12,5) = pc(18,8) - 3*pc(14,8)

sh(13,5) = 3*pc(12,8) - pc(19,8)

sh(14,5) = 2*pc(20,8) - 3*pc(13,8) - 3*pc(15,8)

sh(15,5) = 4*pc(16,8) - pc(18,8) - pc(14,8)

sh(16,5) = 4*pc(17,8) - pc(12,8) - pc(19,8)

sh(1,6) = pc(1,9)

sh(2,6) = pc(2,9)

sh(3,6) = pc(3,9)

sh(4,6) = pc(4,9)

sh(5,6) = pc(8,9)

sh(6,6) = pc(9,9)

sh(7,6) = pc(10,9)

sh(8,6) = pc(5,9) - pc(6,9)

sh(9,6) = 2*pc(7,9) - pc(5,9) - pc(6,9)

sh(10,6) = pc(11,9)

sh(11,6) = pc(13,9) - pc(15,9)

sh(12,6) = pc(18,9) - 3*pc(14,9)

sh(13,6) = 3*pc(12,9) - pc(19,9)

sh(14,6) = 2*pc(20,9) - 3*pc(13,9) - 3*pc(15,9)

sh(15,6) = 4*pc(16,9) - pc(18,9) - pc(14,9)

sh(16,6) = 4*pc(17,9) - pc(12,9) - pc(19,9)

sh(1,7) = pc(1,10)

sh(2,7) = pc(2,10)

sh(3,7) = pc(3,10)

sh(4,7) = pc(4,10)

sh(5,7) = pc(8,10)

sh(6,7) = pc(9,10)

sh(7,7) = pc(10,10)

sh(8,7) = pc(5,10) - pc(6,10)

sh(9,7) = 2*pc(7,10) - pc(5,10) - pc(6,10)

sh(10,7) = pc(11,10)

sh(11,7) = pc(13,10) - pc(15,10)

sh(12,7) = pc(18,10) - 3*pc(14,10)

sh(13,7) = 3*pc(12,10) - pc(19,10)

sh(14,7) = 2*pc(20,10) - 3*pc(13,10) - 3*pc(15,10)

sh(15,7) = 4*pc(16,10) - pc(18,10) - pc(14,10)

sh(16,7) = 4*pc(17,10) - pc(12,10) - pc(19,10)

sh(1,8) = pc(1,5) - pc(1,6)

sh(2,8) = pc(2,5) - pc(2,6)

sh(3,8) = pc(3,5) - pc(3,6)

sh(4,8) = pc(4,5) - pc(4,6)

sh(5,8) = pc(8,5) - pc(8,6)

sh(6,8) = pc(9,5) - pc(9,6)

sh(7,8) = pc(10,5) - pc(10,6)

sh(8,8) = pc(5,5) - pc(5,6) - pc(6,5) + pc(6,6)

sh(9,8) = 2*pc(7,5) - 2*pc(7,6) - pc(5,5) + pc(5,6) - pc(6,5) + pc(6,6)

sh(10,8) = pc(11,5) - pc(11,6)

sh(11,8) = pc(13,5) - pc(13,6) - pc(15,5) + pc(15,6)

sh(12,8) = pc(18,5) - pc(18,6) - 3*pc(14,5) + 3*pc(14,6)

sh(13,8) = 3*pc(12,5) - 3*pc(12,6) - pc(19,5) + pc(19,6)

sh(14,8) = 2*pc(20,5) - 2*pc(20,6) - 3*pc(13,5) + 3*pc(13,6) - 3*pc(15,5) + 3*&
&pc(15,6)

sh(15,8) = 4*pc(16,5) - 4*pc(16,6) - pc(18,5) + pc(18,6) - pc(14,5) + pc(14,6)

sh(16,8) = 4*pc(17,5) - 4*pc(17,6) - pc(12,5) + pc(12,6) - pc(19,5) + pc(19,6)

sh(1,9) = 2*pc(1,7) - pc(1,5) - pc(1,6)

sh(2,9) = 2*pc(2,7) - pc(2,5) - pc(2,6)

sh(3,9) = 2*pc(3,7) - pc(3,5) - pc(3,6)

sh(4,9) = 2*pc(4,7) - pc(4,5) - pc(4,6)

sh(5,9) = 2*pc(8,7) - pc(8,5) - pc(8,6)

sh(6,9) = 2*pc(9,7) - pc(9,5) - pc(9,6)

sh(7,9) = 2*pc(10,7) - pc(10,5) - pc(10,6)

sh(8,9) = 2*pc(5,7) - pc(5,5) - pc(5,6) - 2*pc(6,7) + pc(6,5) + pc(6,6)

sh(9,9) = 4*pc(7,7) - 2*pc(7,5) - 2*pc(7,6) - 2*pc(5,7) + pc(5,5) + pc(5,6) - &
&2*pc(6,7) + pc(6,5) + pc(6,6)

sh(10,9) = 2*pc(11,7) - pc(11,5) - pc(11,6)

sh(11,9) = 2*pc(13,7) - pc(13,5) - pc(13,6) - 2*pc(15,7) + pc(15,5) + pc(15,6)

sh(12,9) = 2*pc(18,7) - pc(18,5) - pc(18,6) - 6*pc(14,7) + 3*pc(14,5) + 3*pc(1&
&4,6)

sh(13,9) = 6*pc(12,7) - 3*pc(12,5) - 3*pc(12,6) - 2*pc(19,7) + pc(19,5) + pc(1&
&9,6)

sh(14,9) = 4*pc(20,7) - 2*pc(20,5) - 2*pc(20,6) - 6*pc(13,7) + 3*pc(13,5) + 3*&
&pc(13,6) - 6*pc(15,7) + 3*pc(15,5) + 3*pc(15,6)

sh(15,9) = 8*pc(16,7) - 4*pc(16,5) - 4*pc(16,6) - 2*pc(18,7) + pc(18,5) + pc(1&
&8,6) - 2*pc(14,7) + pc(14,5) + pc(14,6)

sh(16,9) = 8*pc(17,7) - 4*pc(17,5) - 4*pc(17,6) - 2*pc(12,7) + pc(12,5) + pc(1&
&2,6) - 2*pc(19,7) + pc(19,5) + pc(19,6)

sh(1,10) = pc(1,11)

sh(2,10) = pc(2,11)

sh(3,10) = pc(3,11)

sh(4,10) = pc(4,11)

sh(5,10) = pc(8,11)

sh(6,10) = pc(9,11)

sh(7,10) = pc(10,11)

sh(8,10) = pc(5,11) - pc(6,11)

sh(9,10) = 2*pc(7,11) - pc(5,11) - pc(6,11)

sh(10,10) = pc(11,11)

sh(11,10) = pc(13,11) - pc(15,11)

sh(12,10) = pc(18,11) - 3*pc(14,11)

sh(13,10) = 3*pc(12,11) - pc(19,11)

sh(14,10) = 2*pc(20,11) - 3*pc(13,11) - 3*pc(15,11)

sh(15,10) = 4*pc(16,11) - pc(18,11) - pc(14,11)

sh(16,10) = 4*pc(17,11) - pc(12,11) - pc(19,11)

sh(1,11) = pc(1,13) - pc(1,15)

sh(2,11) = pc(2,13) - pc(2,15)

sh(3,11) = pc(3,13) - pc(3,15)

sh(4,11) = pc(4,13) - pc(4,15)

sh(5,11) = pc(8,13) - pc(8,15)

sh(6,11) = pc(9,13) - pc(9,15)

sh(7,11) = pc(10,13) - pc(10,15)

sh(8,11) = pc(5,13) - pc(5,15) - pc(6,13) + pc(6,15)

sh(9,11) = 2*pc(7,13) - 2*pc(7,15) - pc(5,13) + pc(5,15) - pc(6,13) + pc(6,15)

sh(10,11) = pc(11,13) - pc(11,15)

sh(11,11) = pc(13,13) - pc(13,15) - pc(15,13) + pc(15,15)

sh(12,11) = pc(18,13) - pc(18,15) - 3*pc(14,13) + 3*pc(14,15)

sh(13,11) = 3*pc(12,13) - 3*pc(12,15) - pc(19,13) + pc(19,15)

sh(14,11) = 2*pc(20,13) - 2*pc(20,15) - 3*pc(13,13) + 3*pc(13,15) - 3*pc(15,13&
&) + 3*pc(15,15)

sh(15,11) = 4*pc(16,13) - 4*pc(16,15) - pc(18,13) + pc(18,15) - pc(14,13) + pc&
&(14,15)

sh(16,11) = 4*pc(17,13) - 4*pc(17,15) - pc(12,13) + pc(12,15) - pc(19,13) + pc&
&(19,15)

sh(1,12) = pc(1,18) - 3*pc(1,14)

sh(2,12) = pc(2,18) - 3*pc(2,14)

sh(3,12) = pc(3,18) - 3*pc(3,14)

sh(4,12) = pc(4,18) - 3*pc(4,14)

sh(5,12) = pc(8,18) - 3*pc(8,14)

sh(6,12) = pc(9,18) - 3*pc(9,14)

sh(7,12) = pc(10,18) - 3*pc(10,14)

sh(8,12) = pc(5,18) - 3*pc(5,14) - pc(6,18) + 3*pc(6,14)

sh(9,12) = 2*pc(7,18) - 6*pc(7,14) - pc(5,18) + 3*pc(5,14) - pc(6,18) + 3*pc(6&
&,14)

sh(10,12) = pc(11,18) - 3*pc(11,14)

sh(11,12) = pc(13,18) - 3*pc(13,14) - pc(15,18) + 3*pc(15,14)

sh(12,12) = pc(18,18) - 3*pc(18,14) - 3*pc(14,18) + 9*pc(14,14)

sh(13,12) = 3*pc(12,18) - 9*pc(12,14) - pc(19,18) + 3*pc(19,14)

sh(14,12) = 2*pc(20,18) - 6*pc(20,14) - 3*pc(13,18) + 9*pc(13,14) - 3*pc(15,18&
&) + 9*pc(15,14)

sh(15,12) = 4*pc(16,18) - 12*pc(16,14) - pc(18,18) + 3*pc(18,14) - pc(14,18) +&
& 3*pc(14,14)

sh(16,12) = 4*pc(17,18) - 12*pc(17,14) - pc(12,18) + 3*pc(12,14) - pc(19,18) +&
& 3*pc(19,14)

sh(1,13) = 3*pc(1,12) - pc(1,19)

sh(2,13) = 3*pc(2,12) - pc(2,19)

sh(3,13) = 3*pc(3,12) - pc(3,19)

sh(4,13) = 3*pc(4,12) - pc(4,19)

sh(5,13) = 3*pc(8,12) - pc(8,19)

sh(6,13) = 3*pc(9,12) - pc(9,19)

sh(7,13) = 3*pc(10,12) - pc(10,19)

sh(8,13) = 3*pc(5,12) - pc(5,19) - 3*pc(6,12) + pc(6,19)

sh(9,13) = 6*pc(7,12) - 2*pc(7,19) - 3*pc(5,12) + pc(5,19) - 3*pc(6,12) + pc(6&
&,19)

sh(10,13) = 3*pc(11,12) - pc(11,19)

sh(11,13) = 3*pc(13,12) - pc(13,19) - 3*pc(15,12) + pc(15,19)

sh(12,13) = 3*pc(18,12) - pc(18,19) - 9*pc(14,12) + 3*pc(14,19)

sh(13,13) = 9*pc(12,12) - 3*pc(12,19) - 3*pc(19,12) + pc(19,19)

sh(14,13) = 6*pc(20,12) - 2*pc(20,19) - 9*pc(13,12) + 3*pc(13,19) - 9*pc(15,12&
&) + 3*pc(15,19)

sh(15,13) = 12*pc(16,12) - 4*pc(16,19) - 3*pc(18,12) + pc(18,19) - 3*pc(14,12)&
& + pc(14,19)

sh(16,13) = 12*pc(17,12) - 4*pc(17,19) - 3*pc(12,12) + pc(12,19) - 3*pc(19,12)&
& + pc(19,19)

sh(1,14) = 2*pc(1,20) - 3*pc(1,13) - 3*pc(1,15)

sh(2,14) = 2*pc(2,20) - 3*pc(2,13) - 3*pc(2,15)

sh(3,14) = 2*pc(3,20) - 3*pc(3,13) - 3*pc(3,15)

sh(4,14) = 2*pc(4,20) - 3*pc(4,13) - 3*pc(4,15)

sh(5,14) = 2*pc(8,20) - 3*pc(8,13) - 3*pc(8,15)

sh(6,14) = 2*pc(9,20) - 3*pc(9,13) - 3*pc(9,15)

sh(7,14) = 2*pc(10,20) - 3*pc(10,13) - 3*pc(10,15)

sh(8,14) = 2*pc(5,20) - 3*pc(5,13) - 3*pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*&
&pc(6,15)

sh(9,14) = 4*pc(7,20) - 6*pc(7,13) - 6*pc(7,15) - 2*pc(5,20) + 3*pc(5,13) + 3*&
&pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*pc(6,15)

sh(10,14) = 2*pc(11,20) - 3*pc(11,13) - 3*pc(11,15)

sh(11,14) = 2*pc(13,20) - 3*pc(13,13) - 3*pc(13,15) - 2*pc(15,20) + 3*pc(15,13&
&) + 3*pc(15,15)

sh(12,14) = 2*pc(18,20) - 3*pc(18,13) - 3*pc(18,15) - 6*pc(14,20) + 9*pc(14,13&
&) + 9*pc(14,15)

sh(13,14) = 6*pc(12,20) - 9*pc(12,13) - 9*pc(12,15) - 2*pc(19,20) + 3*pc(19,13&
&) + 3*pc(19,15)

sh(14,14) = 4*pc(20,20) - 6*pc(20,13) - 6*pc(20,15) - 6*pc(13,20) + 9*pc(13,13&
&) + 9*pc(13,15) - 6*pc(15,20) + 9*pc(15,13) + 9*pc(15,15)

sh(15,14) = 8*pc(16,20) - 12*pc(16,13) - 12*pc(16,15) - 2*pc(18,20) + 3*pc(18,&
&13) + 3*pc(18,15) - 2*pc(14,20) + 3*pc(14,13) + 3*pc(14,15)

sh(16,14) = 8*pc(17,20) - 12*pc(17,13) - 12*pc(17,15) - 2*pc(12,20) + 3*pc(12,&
&13) + 3*pc(12,15) - 2*pc(19,20) + 3*pc(19,13) + 3*pc(19,15)

sh(1,15) = 4*pc(1,16) - pc(1,18) - pc(1,14)

sh(2,15) = 4*pc(2,16) - pc(2,18) - pc(2,14)

sh(3,15) = 4*pc(3,16) - pc(3,18) - pc(3,14)

sh(4,15) = 4*pc(4,16) - pc(4,18) - pc(4,14)

sh(5,15) = 4*pc(8,16) - pc(8,18) - pc(8,14)

sh(6,15) = 4*pc(9,16) - pc(9,18) - pc(9,14)

sh(7,15) = 4*pc(10,16) - pc(10,18) - pc(10,14)

sh(8,15) = 4*pc(5,16) - pc(5,18) - pc(5,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(9,15) = 8*pc(7,16) - 2*pc(7,18) - 2*pc(7,14) - 4*pc(5,16) + pc(5,18) + pc(5&
&,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(10,15) = 4*pc(11,16) - pc(11,18) - pc(11,14)

sh(11,15) = 4*pc(13,16) - pc(13,18) - pc(13,14) - 4*pc(15,16) + pc(15,18) + pc&
&(15,14)

sh(12,15) = 4*pc(18,16) - pc(18,18) - pc(18,14) - 12*pc(14,16) + 3*pc(14,18) +&
& 3*pc(14,14)

sh(13,15) = 12*pc(12,16) - 3*pc(12,18) - 3*pc(12,14) - 4*pc(19,16) + pc(19,18)&
& + pc(19,14)

sh(14,15) = 8*pc(20,16) - 2*pc(20,18) - 2*pc(20,14) - 12*pc(13,16) + 3*pc(13,1&
&8) + 3*pc(13,14) - 12*pc(15,16) + 3*pc(15,18) + 3*pc(15,14)

sh(15,15) = 16*pc(16,16) - 4*pc(16,18) - 4*pc(16,14) - 4*pc(18,16) + pc(18,18)&
& + pc(18,14) - 4*pc(14,16) + pc(14,18) + pc(14,14)

sh(16,15) = 16*pc(17,16) - 4*pc(17,18) - 4*pc(17,14) - 4*pc(12,16) + pc(12,18)&
& + pc(12,14) - 4*pc(19,16) + pc(19,18) + pc(19,14)

sh(1,16) = 4*pc(1,17) - pc(1,12) - pc(1,19)

sh(2,16) = 4*pc(2,17) - pc(2,12) - pc(2,19)

sh(3,16) = 4*pc(3,17) - pc(3,12) - pc(3,19)

sh(4,16) = 4*pc(4,17) - pc(4,12) - pc(4,19)

sh(5,16) = 4*pc(8,17) - pc(8,12) - pc(8,19)

sh(6,16) = 4*pc(9,17) - pc(9,12) - pc(9,19)

sh(7,16) = 4*pc(10,17) - pc(10,12) - pc(10,19)

sh(8,16) = 4*pc(5,17) - pc(5,12) - pc(5,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(9,16) = 8*pc(7,17) - 2*pc(7,12) - 2*pc(7,19) - 4*pc(5,17) + pc(5,12) + pc(5&
&,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(10,16) = 4*pc(11,17) - pc(11,12) - pc(11,19)

sh(11,16) = 4*pc(13,17) - pc(13,12) - pc(13,19) - 4*pc(15,17) + pc(15,12) + pc&
&(15,19)

sh(12,16) = 4*pc(18,17) - pc(18,12) - pc(18,19) - 12*pc(14,17) + 3*pc(14,12) +&
& 3*pc(14,19)

sh(13,16) = 12*pc(12,17) - 3*pc(12,12) - 3*pc(12,19) - 4*pc(19,17) + pc(19,12)&
& + pc(19,19)

sh(14,16) = 8*pc(20,17) - 2*pc(20,12) - 2*pc(20,19) - 12*pc(13,17) + 3*pc(13,1&
&2) + 3*pc(13,19) - 12*pc(15,17) + 3*pc(15,12) + 3*pc(15,19)

sh(15,16) = 16*pc(16,17) - 4*pc(16,12) - 4*pc(16,19) - 4*pc(18,17) + pc(18,12)&
& + pc(18,19) - 4*pc(14,17) + pc(14,12) + pc(14,19)

sh(16,16) = 16*pc(17,17) - 4*pc(17,12) - 4*pc(17,19) - 4*pc(12,17) + pc(12,12)&
& + pc(12,19) - 4*pc(19,17) + pc(19,12) + pc(19,19)


   end subroutine overlap2CIntgNum

   subroutine kinetic2CIntgAna(a1,a2,A,B,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   real (kind=double), dimension (20,20), intent(out) :: pc
   real (kind=double), dimension (16,16), intent(out) :: sh

   ! Define local variables.
   real (kind=double), dimension (20,20) :: pc_ol
   real (kind=double), dimension (3) :: P, PA, PB, d
   real (kind=double) :: zeta, inv_2zeta, xi, preFactorOL, preFactorKE
   real (kind=double) :: inv_2zeta_a, inv_2zeta_b

   ! Initialize local variables.
   zeta = a1 + a2
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   inv_2zeta_a = 1.0d0 / (2.0d0 * a1)
   inv_2zeta_b = 1.0d0 / (2.0d0 * a2)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   d = A - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d*d))
   preFactorKE = xi*(3 - 2*xi*sum(d*d))*preFactorOL

pc_ol(1,1) = preFactorOL

pc_ol(2,1) = PA(1)*preFactorOL

pc_ol(3,1) = PA(2)*preFactorOL

pc_ol(4,1) = PA(3)*preFactorOL

pc_ol(5,1) = PA(1)*pc_ol(2,1) + inv_2zeta*(preFactorOL)

pc_ol(6,1) = PA(2)*pc_ol(3,1) + inv_2zeta*(preFactorOL)

pc_ol(7,1) = PA(3)*pc_ol(4,1) + inv_2zeta*(preFactorOL)

pc_ol(8,1) = PA(2)*pc_ol(2,1)

pc_ol(9,1) = PA(3)*pc_ol(2,1)

pc_ol(10,1) = PA(3)*pc_ol(3,1)

pc_ol(11,1) = PA(3)*pc_ol(8,1)

pc_ol(12,1) = PA(2)*pc_ol(5,1)

pc_ol(13,1) = PA(3)*pc_ol(5,1)

pc_ol(14,1) = PA(1)*pc_ol(6,1)

pc_ol(15,1) = PA(3)*pc_ol(6,1)

pc_ol(16,1) = PA(1)*pc_ol(7,1)

pc_ol(17,1) = PA(2)*pc_ol(7,1)

pc_ol(18,1) = PA(1)*pc_ol(5,1) + inv_2zeta*(2*pc_ol(2,1))

pc_ol(19,1) = PA(2)*pc_ol(6,1) + inv_2zeta*(2*pc_ol(3,1))

pc_ol(20,1) = PA(3)*pc_ol(7,1) + inv_2zeta*(2*pc_ol(4,1))

pc_ol(1,2) = PB(1)*preFactorOL

pc_ol(2,2) = PA(1)*pc_ol(1,2) + inv_2zeta*(preFactorOL)

pc_ol(3,2) = PA(2)*pc_ol(1,2)

pc_ol(4,2) = PA(3)*pc_ol(1,2)

pc_ol(5,2) = PA(1)*pc_ol(2,2) + inv_2zeta*(pc_ol(1,2) + pc_ol(2,1))

pc_ol(6,2) = PA(2)*pc_ol(3,2) + inv_2zeta*(pc_ol(1,2))

pc_ol(7,2) = PA(3)*pc_ol(4,2) + inv_2zeta*(pc_ol(1,2))

pc_ol(8,2) = PA(2)*pc_ol(2,2)

pc_ol(9,2) = PA(3)*pc_ol(2,2)

pc_ol(10,2) = PA(3)*pc_ol(3,2)

pc_ol(11,2) = PA(3)*pc_ol(8,2)

pc_ol(12,2) = PA(2)*pc_ol(5,2)

pc_ol(13,2) = PA(3)*pc_ol(5,2)

pc_ol(14,2) = PA(1)*pc_ol(6,2) + inv_2zeta*(pc_ol(6,1))

pc_ol(15,2) = PA(3)*pc_ol(6,2)

pc_ol(16,2) = PA(1)*pc_ol(7,2) + inv_2zeta*(pc_ol(7,1))

pc_ol(17,2) = PA(2)*pc_ol(7,2)

pc_ol(18,2) = PA(1)*pc_ol(5,2) + inv_2zeta*(2*pc_ol(2,2) + pc_ol(5,1))

pc_ol(19,2) = PA(2)*pc_ol(6,2) + inv_2zeta*(2*pc_ol(3,2))

pc_ol(20,2) = PA(3)*pc_ol(7,2) + inv_2zeta*(2*pc_ol(4,2))

pc_ol(1,3) = PB(2)*preFactorOL

pc_ol(2,3) = PA(1)*pc_ol(1,3)

pc_ol(3,3) = PA(2)*pc_ol(1,3) + inv_2zeta*(preFactorOL)

pc_ol(4,3) = PA(3)*pc_ol(1,3)

pc_ol(5,3) = PA(1)*pc_ol(2,3) + inv_2zeta*(pc_ol(1,3))

pc_ol(6,3) = PA(2)*pc_ol(3,3) + inv_2zeta*(pc_ol(1,3) + pc_ol(3,1))

pc_ol(7,3) = PA(3)*pc_ol(4,3) + inv_2zeta*(pc_ol(1,3))

pc_ol(8,3) = PA(2)*pc_ol(2,3) + inv_2zeta*(pc_ol(2,1))

pc_ol(9,3) = PA(3)*pc_ol(2,3)

pc_ol(10,3) = PA(3)*pc_ol(3,3)

pc_ol(11,3) = PA(3)*pc_ol(8,3)

pc_ol(12,3) = PA(2)*pc_ol(5,3) + inv_2zeta*(pc_ol(5,1))

pc_ol(13,3) = PA(3)*pc_ol(5,3)

pc_ol(14,3) = PA(1)*pc_ol(6,3)

pc_ol(15,3) = PA(3)*pc_ol(6,3)

pc_ol(16,3) = PA(1)*pc_ol(7,3)

pc_ol(17,3) = PA(2)*pc_ol(7,3) + inv_2zeta*(pc_ol(7,1))

pc_ol(18,3) = PA(1)*pc_ol(5,3) + inv_2zeta*(2*pc_ol(2,3))

pc_ol(19,3) = PA(2)*pc_ol(6,3) + inv_2zeta*(2*pc_ol(3,3) + pc_ol(6,1))

pc_ol(20,3) = PA(3)*pc_ol(7,3) + inv_2zeta*(2*pc_ol(4,3))

pc_ol(1,4) = PB(3)*preFactorOL

pc_ol(2,4) = PA(1)*pc_ol(1,4)

pc_ol(3,4) = PA(2)*pc_ol(1,4)

pc_ol(4,4) = PA(3)*pc_ol(1,4) + inv_2zeta*(preFactorOL)

pc_ol(5,4) = PA(1)*pc_ol(2,4) + inv_2zeta*(pc_ol(1,4))

pc_ol(6,4) = PA(2)*pc_ol(3,4) + inv_2zeta*(pc_ol(1,4))

pc_ol(7,4) = PA(3)*pc_ol(4,4) + inv_2zeta*(pc_ol(1,4) + pc_ol(4,1))

pc_ol(8,4) = PA(2)*pc_ol(2,4)

pc_ol(9,4) = PA(3)*pc_ol(2,4) + inv_2zeta*(pc_ol(2,1))

pc_ol(10,4) = PA(3)*pc_ol(3,4) + inv_2zeta*(pc_ol(3,1))

pc_ol(11,4) = PA(3)*pc_ol(8,4) + inv_2zeta*(pc_ol(8,1))

pc_ol(12,4) = PA(2)*pc_ol(5,4)

pc_ol(13,4) = PA(3)*pc_ol(5,4) + inv_2zeta*(pc_ol(5,1))

pc_ol(14,4) = PA(1)*pc_ol(6,4)

pc_ol(15,4) = PA(3)*pc_ol(6,4) + inv_2zeta*(pc_ol(6,1))

pc_ol(16,4) = PA(1)*pc_ol(7,4)

pc_ol(17,4) = PA(2)*pc_ol(7,4)

pc_ol(18,4) = PA(1)*pc_ol(5,4) + inv_2zeta*(2*pc_ol(2,4))

pc_ol(19,4) = PA(2)*pc_ol(6,4) + inv_2zeta*(2*pc_ol(3,4))

pc_ol(20,4) = PA(3)*pc_ol(7,4) + inv_2zeta*(2*pc_ol(4,4) + pc_ol(7,1))

pc_ol(1,5) = PB(1)*pc_ol(1,2) + inv_2zeta*(preFactorOL)

pc_ol(2,5) = PA(1)*pc_ol(1,5) + inv_2zeta*(2*pc_ol(1,2))

pc_ol(3,5) = PA(2)*pc_ol(1,5)

pc_ol(4,5) = PA(3)*pc_ol(1,5)

pc_ol(5,5) = PA(1)*pc_ol(2,5) + inv_2zeta*(pc_ol(1,5) + 2*pc_ol(2,2))

pc_ol(6,5) = PA(2)*pc_ol(3,5) + inv_2zeta*(pc_ol(1,5))

pc_ol(7,5) = PA(3)*pc_ol(4,5) + inv_2zeta*(pc_ol(1,5))

pc_ol(8,5) = PA(2)*pc_ol(2,5)

pc_ol(9,5) = PA(3)*pc_ol(2,5)

pc_ol(10,5) = PA(3)*pc_ol(3,5)

pc_ol(11,5) = PA(3)*pc_ol(8,5)

pc_ol(12,5) = PA(2)*pc_ol(5,5)

pc_ol(13,5) = PA(3)*pc_ol(5,5)

pc_ol(14,5) = PA(1)*pc_ol(6,5) + inv_2zeta*(2*pc_ol(6,2))

pc_ol(15,5) = PA(3)*pc_ol(6,5)

pc_ol(16,5) = PA(1)*pc_ol(7,5) + inv_2zeta*(2*pc_ol(7,2))

pc_ol(17,5) = PA(2)*pc_ol(7,5)

pc_ol(18,5) = PA(1)*pc_ol(5,5) + inv_2zeta*(2*pc_ol(2,5) + 2*pc_ol(5,2))

pc_ol(19,5) = PA(2)*pc_ol(6,5) + inv_2zeta*(2*pc_ol(3,5))

pc_ol(20,5) = PA(3)*pc_ol(7,5) + inv_2zeta*(2*pc_ol(4,5))

pc_ol(1,6) = PB(2)*pc_ol(1,3) + inv_2zeta*(preFactorOL)

pc_ol(2,6) = PA(1)*pc_ol(1,6)

pc_ol(3,6) = PA(2)*pc_ol(1,6) + inv_2zeta*(2*pc_ol(1,3))

pc_ol(4,6) = PA(3)*pc_ol(1,6)

pc_ol(5,6) = PA(1)*pc_ol(2,6) + inv_2zeta*(pc_ol(1,6))

pc_ol(6,6) = PA(2)*pc_ol(3,6) + inv_2zeta*(pc_ol(1,6) + 2*pc_ol(3,3))

pc_ol(7,6) = PA(3)*pc_ol(4,6) + inv_2zeta*(pc_ol(1,6))

pc_ol(8,6) = PA(2)*pc_ol(2,6) + inv_2zeta*(2*pc_ol(2,3))

pc_ol(9,6) = PA(3)*pc_ol(2,6)

pc_ol(10,6) = PA(3)*pc_ol(3,6)

pc_ol(11,6) = PA(3)*pc_ol(8,6)

pc_ol(12,6) = PA(2)*pc_ol(5,6) + inv_2zeta*(2*pc_ol(5,3))

pc_ol(13,6) = PA(3)*pc_ol(5,6)

pc_ol(14,6) = PA(1)*pc_ol(6,6)

pc_ol(15,6) = PA(3)*pc_ol(6,6)

pc_ol(16,6) = PA(1)*pc_ol(7,6)

pc_ol(17,6) = PA(2)*pc_ol(7,6) + inv_2zeta*(2*pc_ol(7,3))

pc_ol(18,6) = PA(1)*pc_ol(5,6) + inv_2zeta*(2*pc_ol(2,6))

pc_ol(19,6) = PA(2)*pc_ol(6,6) + inv_2zeta*(2*pc_ol(3,6) + 2*pc_ol(6,3))

pc_ol(20,6) = PA(3)*pc_ol(7,6) + inv_2zeta*(2*pc_ol(4,6))

pc_ol(1,7) = PB(3)*pc_ol(1,4) + inv_2zeta*(preFactorOL)

pc_ol(2,7) = PA(1)*pc_ol(1,7)

pc_ol(3,7) = PA(2)*pc_ol(1,7)

pc_ol(4,7) = PA(3)*pc_ol(1,7) + inv_2zeta*(2*pc_ol(1,4))

pc_ol(5,7) = PA(1)*pc_ol(2,7) + inv_2zeta*(pc_ol(1,7))

pc_ol(6,7) = PA(2)*pc_ol(3,7) + inv_2zeta*(pc_ol(1,7))

pc_ol(7,7) = PA(3)*pc_ol(4,7) + inv_2zeta*(pc_ol(1,7) + 2*pc_ol(4,4))

pc_ol(8,7) = PA(2)*pc_ol(2,7)

pc_ol(9,7) = PA(3)*pc_ol(2,7) + inv_2zeta*(2*pc_ol(2,4))

pc_ol(10,7) = PA(3)*pc_ol(3,7) + inv_2zeta*(2*pc_ol(3,4))

pc_ol(11,7) = PA(3)*pc_ol(8,7) + inv_2zeta*(2*pc_ol(8,4))

pc_ol(12,7) = PA(2)*pc_ol(5,7)

pc_ol(13,7) = PA(3)*pc_ol(5,7) + inv_2zeta*(2*pc_ol(5,4))

pc_ol(14,7) = PA(1)*pc_ol(6,7)

pc_ol(15,7) = PA(3)*pc_ol(6,7) + inv_2zeta*(2*pc_ol(6,4))

pc_ol(16,7) = PA(1)*pc_ol(7,7)

pc_ol(17,7) = PA(2)*pc_ol(7,7)

pc_ol(18,7) = PA(1)*pc_ol(5,7) + inv_2zeta*(2*pc_ol(2,7))

pc_ol(19,7) = PA(2)*pc_ol(6,7) + inv_2zeta*(2*pc_ol(3,7))

pc_ol(20,7) = PA(3)*pc_ol(7,7) + inv_2zeta*(2*pc_ol(4,7) + 2*pc_ol(7,4))

pc_ol(1,8) = PB(2)*pc_ol(1,2)

pc_ol(2,8) = PA(1)*pc_ol(1,8) + inv_2zeta*(pc_ol(1,3))

pc_ol(3,8) = PA(2)*pc_ol(1,8) + inv_2zeta*(pc_ol(1,2))

pc_ol(4,8) = PA(3)*pc_ol(1,8)

pc_ol(5,8) = PA(1)*pc_ol(2,8) + inv_2zeta*(pc_ol(1,8) + pc_ol(2,3))

pc_ol(6,8) = PA(2)*pc_ol(3,8) + inv_2zeta*(pc_ol(1,8) + pc_ol(3,2))

pc_ol(7,8) = PA(3)*pc_ol(4,8) + inv_2zeta*(pc_ol(1,8))

pc_ol(8,8) = PA(2)*pc_ol(2,8) + inv_2zeta*(pc_ol(2,2))

pc_ol(9,8) = PA(3)*pc_ol(2,8)

pc_ol(10,8) = PA(3)*pc_ol(3,8)

pc_ol(11,8) = PA(3)*pc_ol(8,8)

pc_ol(12,8) = PA(2)*pc_ol(5,8) + inv_2zeta*(pc_ol(5,2))

pc_ol(13,8) = PA(3)*pc_ol(5,8)

pc_ol(14,8) = PA(1)*pc_ol(6,8) + inv_2zeta*(pc_ol(6,3))

pc_ol(15,8) = PA(3)*pc_ol(6,8)

pc_ol(16,8) = PA(1)*pc_ol(7,8) + inv_2zeta*(pc_ol(7,3))

pc_ol(17,8) = PA(2)*pc_ol(7,8) + inv_2zeta*(pc_ol(7,2))

pc_ol(18,8) = PA(1)*pc_ol(5,8) + inv_2zeta*(2*pc_ol(2,8) + pc_ol(5,3))

pc_ol(19,8) = PA(2)*pc_ol(6,8) + inv_2zeta*(2*pc_ol(3,8) + pc_ol(6,2))

pc_ol(20,8) = PA(3)*pc_ol(7,8) + inv_2zeta*(2*pc_ol(4,8))

pc_ol(1,9) = PB(3)*pc_ol(1,2)

pc_ol(2,9) = PA(1)*pc_ol(1,9) + inv_2zeta*(pc_ol(1,4))

pc_ol(3,9) = PA(2)*pc_ol(1,9)

pc_ol(4,9) = PA(3)*pc_ol(1,9) + inv_2zeta*(pc_ol(1,2))

pc_ol(5,9) = PA(1)*pc_ol(2,9) + inv_2zeta*(pc_ol(1,9) + pc_ol(2,4))

pc_ol(6,9) = PA(2)*pc_ol(3,9) + inv_2zeta*(pc_ol(1,9))

pc_ol(7,9) = PA(3)*pc_ol(4,9) + inv_2zeta*(pc_ol(1,9) + pc_ol(4,2))

pc_ol(8,9) = PA(2)*pc_ol(2,9)

pc_ol(9,9) = PA(3)*pc_ol(2,9) + inv_2zeta*(pc_ol(2,2))

pc_ol(10,9) = PA(3)*pc_ol(3,9) + inv_2zeta*(pc_ol(3,2))

pc_ol(11,9) = PA(3)*pc_ol(8,9) + inv_2zeta*(pc_ol(8,2))

pc_ol(12,9) = PA(2)*pc_ol(5,9)

pc_ol(13,9) = PA(3)*pc_ol(5,9) + inv_2zeta*(pc_ol(5,2))

pc_ol(14,9) = PA(1)*pc_ol(6,9) + inv_2zeta*(pc_ol(6,4))

pc_ol(15,9) = PA(3)*pc_ol(6,9) + inv_2zeta*(pc_ol(6,2))

pc_ol(16,9) = PA(1)*pc_ol(7,9) + inv_2zeta*(pc_ol(7,4))

pc_ol(17,9) = PA(2)*pc_ol(7,9)

pc_ol(18,9) = PA(1)*pc_ol(5,9) + inv_2zeta*(2*pc_ol(2,9) + pc_ol(5,4))

pc_ol(19,9) = PA(2)*pc_ol(6,9) + inv_2zeta*(2*pc_ol(3,9))

pc_ol(20,9) = PA(3)*pc_ol(7,9) + inv_2zeta*(2*pc_ol(4,9) + pc_ol(7,2))

pc_ol(1,10) = PB(3)*pc_ol(1,3)

pc_ol(2,10) = PA(1)*pc_ol(1,10)

pc_ol(3,10) = PA(2)*pc_ol(1,10) + inv_2zeta*(pc_ol(1,4))

pc_ol(4,10) = PA(3)*pc_ol(1,10) + inv_2zeta*(pc_ol(1,3))

pc_ol(5,10) = PA(1)*pc_ol(2,10) + inv_2zeta*(pc_ol(1,10))

pc_ol(6,10) = PA(2)*pc_ol(3,10) + inv_2zeta*(pc_ol(1,10) + pc_ol(3,4))

pc_ol(7,10) = PA(3)*pc_ol(4,10) + inv_2zeta*(pc_ol(1,10) + pc_ol(4,3))

pc_ol(8,10) = PA(2)*pc_ol(2,10) + inv_2zeta*(pc_ol(2,4))

pc_ol(9,10) = PA(3)*pc_ol(2,10) + inv_2zeta*(pc_ol(2,3))

pc_ol(10,10) = PA(3)*pc_ol(3,10) + inv_2zeta*(pc_ol(3,3))

pc_ol(11,10) = PA(3)*pc_ol(8,10) + inv_2zeta*(pc_ol(8,3))

pc_ol(12,10) = PA(2)*pc_ol(5,10) + inv_2zeta*(pc_ol(5,4))

pc_ol(13,10) = PA(3)*pc_ol(5,10) + inv_2zeta*(pc_ol(5,3))

pc_ol(14,10) = PA(1)*pc_ol(6,10)

pc_ol(15,10) = PA(3)*pc_ol(6,10) + inv_2zeta*(pc_ol(6,3))

pc_ol(16,10) = PA(1)*pc_ol(7,10)

pc_ol(17,10) = PA(2)*pc_ol(7,10) + inv_2zeta*(pc_ol(7,4))

pc_ol(18,10) = PA(1)*pc_ol(5,10) + inv_2zeta*(2*pc_ol(2,10))

pc_ol(19,10) = PA(2)*pc_ol(6,10) + inv_2zeta*(2*pc_ol(3,10) + pc_ol(6,4))

pc_ol(20,10) = PA(3)*pc_ol(7,10) + inv_2zeta*(2*pc_ol(4,10) + pc_ol(7,3))

pc_ol(1,11) = PB(3)*pc_ol(1,8)

pc_ol(2,11) = PA(1)*pc_ol(1,11) + inv_2zeta*(pc_ol(1,10))

pc_ol(3,11) = PA(2)*pc_ol(1,11) + inv_2zeta*(pc_ol(1,9))

pc_ol(4,11) = PA(3)*pc_ol(1,11) + inv_2zeta*(pc_ol(1,8))

pc_ol(5,11) = PA(1)*pc_ol(2,11) + inv_2zeta*(pc_ol(1,11) + pc_ol(2,10))

pc_ol(6,11) = PA(2)*pc_ol(3,11) + inv_2zeta*(pc_ol(1,11) + pc_ol(3,9))

pc_ol(7,11) = PA(3)*pc_ol(4,11) + inv_2zeta*(pc_ol(1,11) + pc_ol(4,8))

pc_ol(8,11) = PA(2)*pc_ol(2,11) + inv_2zeta*(pc_ol(2,9))

pc_ol(9,11) = PA(3)*pc_ol(2,11) + inv_2zeta*(pc_ol(2,8))

pc_ol(10,11) = PA(3)*pc_ol(3,11) + inv_2zeta*(pc_ol(3,8))

pc_ol(11,11) = PA(3)*pc_ol(8,11) + inv_2zeta*(pc_ol(8,8))

pc_ol(12,11) = PA(2)*pc_ol(5,11) + inv_2zeta*(pc_ol(5,9))

pc_ol(13,11) = PA(3)*pc_ol(5,11) + inv_2zeta*(pc_ol(5,8))

pc_ol(14,11) = PA(1)*pc_ol(6,11) + inv_2zeta*(pc_ol(6,10))

pc_ol(15,11) = PA(3)*pc_ol(6,11) + inv_2zeta*(pc_ol(6,8))

pc_ol(16,11) = PA(1)*pc_ol(7,11) + inv_2zeta*(pc_ol(7,10))

pc_ol(17,11) = PA(2)*pc_ol(7,11) + inv_2zeta*(pc_ol(7,9))

pc_ol(18,11) = PA(1)*pc_ol(5,11) + inv_2zeta*(2*pc_ol(2,11) + pc_ol(5,10))

pc_ol(19,11) = PA(2)*pc_ol(6,11) + inv_2zeta*(2*pc_ol(3,11) + pc_ol(6,9))

pc_ol(20,11) = PA(3)*pc_ol(7,11) + inv_2zeta*(2*pc_ol(4,11) + pc_ol(7,8))

pc_ol(1,12) = PB(2)*pc_ol(1,5)

pc_ol(2,12) = PA(1)*pc_ol(1,12) + inv_2zeta*(2*pc_ol(1,8))

pc_ol(3,12) = PA(2)*pc_ol(1,12) + inv_2zeta*(pc_ol(1,5))

pc_ol(4,12) = PA(3)*pc_ol(1,12)

pc_ol(5,12) = PA(1)*pc_ol(2,12) + inv_2zeta*(pc_ol(1,12) + 2*pc_ol(2,8))

pc_ol(6,12) = PA(2)*pc_ol(3,12) + inv_2zeta*(pc_ol(1,12) + pc_ol(3,5))

pc_ol(7,12) = PA(3)*pc_ol(4,12) + inv_2zeta*(pc_ol(1,12))

pc_ol(8,12) = PA(2)*pc_ol(2,12) + inv_2zeta*(pc_ol(2,5))

pc_ol(9,12) = PA(3)*pc_ol(2,12)

pc_ol(10,12) = PA(3)*pc_ol(3,12)

pc_ol(11,12) = PA(3)*pc_ol(8,12)

pc_ol(12,12) = PA(2)*pc_ol(5,12) + inv_2zeta*(pc_ol(5,5))

pc_ol(13,12) = PA(3)*pc_ol(5,12)

pc_ol(14,12) = PA(1)*pc_ol(6,12) + inv_2zeta*(2*pc_ol(6,8))

pc_ol(15,12) = PA(3)*pc_ol(6,12)

pc_ol(16,12) = PA(1)*pc_ol(7,12) + inv_2zeta*(2*pc_ol(7,8))

pc_ol(17,12) = PA(2)*pc_ol(7,12) + inv_2zeta*(pc_ol(7,5))

pc_ol(18,12) = PA(1)*pc_ol(5,12) + inv_2zeta*(2*pc_ol(2,12) + 2*pc_ol(5,8))

pc_ol(19,12) = PA(2)*pc_ol(6,12) + inv_2zeta*(2*pc_ol(3,12) + pc_ol(6,5))

pc_ol(20,12) = PA(3)*pc_ol(7,12) + inv_2zeta*(2*pc_ol(4,12))

pc_ol(1,13) = PB(3)*pc_ol(1,5)

pc_ol(2,13) = PA(1)*pc_ol(1,13) + inv_2zeta*(2*pc_ol(1,9))

pc_ol(3,13) = PA(2)*pc_ol(1,13)

pc_ol(4,13) = PA(3)*pc_ol(1,13) + inv_2zeta*(pc_ol(1,5))

pc_ol(5,13) = PA(1)*pc_ol(2,13) + inv_2zeta*(pc_ol(1,13) + 2*pc_ol(2,9))

pc_ol(6,13) = PA(2)*pc_ol(3,13) + inv_2zeta*(pc_ol(1,13))

pc_ol(7,13) = PA(3)*pc_ol(4,13) + inv_2zeta*(pc_ol(1,13) + pc_ol(4,5))

pc_ol(8,13) = PA(2)*pc_ol(2,13)

pc_ol(9,13) = PA(3)*pc_ol(2,13) + inv_2zeta*(pc_ol(2,5))

pc_ol(10,13) = PA(3)*pc_ol(3,13) + inv_2zeta*(pc_ol(3,5))

pc_ol(11,13) = PA(3)*pc_ol(8,13) + inv_2zeta*(pc_ol(8,5))

pc_ol(12,13) = PA(2)*pc_ol(5,13)

pc_ol(13,13) = PA(3)*pc_ol(5,13) + inv_2zeta*(pc_ol(5,5))

pc_ol(14,13) = PA(1)*pc_ol(6,13) + inv_2zeta*(2*pc_ol(6,9))

pc_ol(15,13) = PA(3)*pc_ol(6,13) + inv_2zeta*(pc_ol(6,5))

pc_ol(16,13) = PA(1)*pc_ol(7,13) + inv_2zeta*(2*pc_ol(7,9))

pc_ol(17,13) = PA(2)*pc_ol(7,13)

pc_ol(18,13) = PA(1)*pc_ol(5,13) + inv_2zeta*(2*pc_ol(2,13) + 2*pc_ol(5,9))

pc_ol(19,13) = PA(2)*pc_ol(6,13) + inv_2zeta*(2*pc_ol(3,13))

pc_ol(20,13) = PA(3)*pc_ol(7,13) + inv_2zeta*(2*pc_ol(4,13) + pc_ol(7,5))

pc_ol(1,14) = PB(1)*pc_ol(1,6)

pc_ol(2,14) = PA(1)*pc_ol(1,14) + inv_2zeta*(pc_ol(1,6))

pc_ol(3,14) = PA(2)*pc_ol(1,14) + inv_2zeta*(2*pc_ol(1,8))

pc_ol(4,14) = PA(3)*pc_ol(1,14)

pc_ol(5,14) = PA(1)*pc_ol(2,14) + inv_2zeta*(pc_ol(1,14) + pc_ol(2,6))

pc_ol(6,14) = PA(2)*pc_ol(3,14) + inv_2zeta*(pc_ol(1,14) + 2*pc_ol(3,8))

pc_ol(7,14) = PA(3)*pc_ol(4,14) + inv_2zeta*(pc_ol(1,14))

pc_ol(8,14) = PA(2)*pc_ol(2,14) + inv_2zeta*(2*pc_ol(2,8))

pc_ol(9,14) = PA(3)*pc_ol(2,14)

pc_ol(10,14) = PA(3)*pc_ol(3,14)

pc_ol(11,14) = PA(3)*pc_ol(8,14)

pc_ol(12,14) = PA(2)*pc_ol(5,14) + inv_2zeta*(2*pc_ol(5,8))

pc_ol(13,14) = PA(3)*pc_ol(5,14)

pc_ol(14,14) = PA(1)*pc_ol(6,14) + inv_2zeta*(pc_ol(6,6))

pc_ol(15,14) = PA(3)*pc_ol(6,14)

pc_ol(16,14) = PA(1)*pc_ol(7,14) + inv_2zeta*(pc_ol(7,6))

pc_ol(17,14) = PA(2)*pc_ol(7,14) + inv_2zeta*(2*pc_ol(7,8))

pc_ol(18,14) = PA(1)*pc_ol(5,14) + inv_2zeta*(2*pc_ol(2,14) + pc_ol(5,6))

pc_ol(19,14) = PA(2)*pc_ol(6,14) + inv_2zeta*(2*pc_ol(3,14) + 2*pc_ol(6,8))

pc_ol(20,14) = PA(3)*pc_ol(7,14) + inv_2zeta*(2*pc_ol(4,14))

pc_ol(1,15) = PB(3)*pc_ol(1,6)

pc_ol(2,15) = PA(1)*pc_ol(1,15)

pc_ol(3,15) = PA(2)*pc_ol(1,15) + inv_2zeta*(2*pc_ol(1,10))

pc_ol(4,15) = PA(3)*pc_ol(1,15) + inv_2zeta*(pc_ol(1,6))

pc_ol(5,15) = PA(1)*pc_ol(2,15) + inv_2zeta*(pc_ol(1,15))

pc_ol(6,15) = PA(2)*pc_ol(3,15) + inv_2zeta*(pc_ol(1,15) + 2*pc_ol(3,10))

pc_ol(7,15) = PA(3)*pc_ol(4,15) + inv_2zeta*(pc_ol(1,15) + pc_ol(4,6))

pc_ol(8,15) = PA(2)*pc_ol(2,15) + inv_2zeta*(2*pc_ol(2,10))

pc_ol(9,15) = PA(3)*pc_ol(2,15) + inv_2zeta*(pc_ol(2,6))

pc_ol(10,15) = PA(3)*pc_ol(3,15) + inv_2zeta*(pc_ol(3,6))

pc_ol(11,15) = PA(3)*pc_ol(8,15) + inv_2zeta*(pc_ol(8,6))

pc_ol(12,15) = PA(2)*pc_ol(5,15) + inv_2zeta*(2*pc_ol(5,10))

pc_ol(13,15) = PA(3)*pc_ol(5,15) + inv_2zeta*(pc_ol(5,6))

pc_ol(14,15) = PA(1)*pc_ol(6,15)

pc_ol(15,15) = PA(3)*pc_ol(6,15) + inv_2zeta*(pc_ol(6,6))

pc_ol(16,15) = PA(1)*pc_ol(7,15)

pc_ol(17,15) = PA(2)*pc_ol(7,15) + inv_2zeta*(2*pc_ol(7,10))

pc_ol(18,15) = PA(1)*pc_ol(5,15) + inv_2zeta*(2*pc_ol(2,15))

pc_ol(19,15) = PA(2)*pc_ol(6,15) + inv_2zeta*(2*pc_ol(3,15) + 2*pc_ol(6,10))

pc_ol(20,15) = PA(3)*pc_ol(7,15) + inv_2zeta*(2*pc_ol(4,15) + pc_ol(7,6))

pc_ol(1,16) = PB(1)*pc_ol(1,7)

pc_ol(2,16) = PA(1)*pc_ol(1,16) + inv_2zeta*(pc_ol(1,7))

pc_ol(3,16) = PA(2)*pc_ol(1,16)

pc_ol(4,16) = PA(3)*pc_ol(1,16) + inv_2zeta*(2*pc_ol(1,9))

pc_ol(5,16) = PA(1)*pc_ol(2,16) + inv_2zeta*(pc_ol(1,16) + pc_ol(2,7))

pc_ol(6,16) = PA(2)*pc_ol(3,16) + inv_2zeta*(pc_ol(1,16))

pc_ol(7,16) = PA(3)*pc_ol(4,16) + inv_2zeta*(pc_ol(1,16) + 2*pc_ol(4,9))

pc_ol(8,16) = PA(2)*pc_ol(2,16)

pc_ol(9,16) = PA(3)*pc_ol(2,16) + inv_2zeta*(2*pc_ol(2,9))

pc_ol(10,16) = PA(3)*pc_ol(3,16) + inv_2zeta*(2*pc_ol(3,9))

pc_ol(11,16) = PA(3)*pc_ol(8,16) + inv_2zeta*(2*pc_ol(8,9))

pc_ol(12,16) = PA(2)*pc_ol(5,16)

pc_ol(13,16) = PA(3)*pc_ol(5,16) + inv_2zeta*(2*pc_ol(5,9))

pc_ol(14,16) = PA(1)*pc_ol(6,16) + inv_2zeta*(pc_ol(6,7))

pc_ol(15,16) = PA(3)*pc_ol(6,16) + inv_2zeta*(2*pc_ol(6,9))

pc_ol(16,16) = PA(1)*pc_ol(7,16) + inv_2zeta*(pc_ol(7,7))

pc_ol(17,16) = PA(2)*pc_ol(7,16)

pc_ol(18,16) = PA(1)*pc_ol(5,16) + inv_2zeta*(2*pc_ol(2,16) + pc_ol(5,7))

pc_ol(19,16) = PA(2)*pc_ol(6,16) + inv_2zeta*(2*pc_ol(3,16))

pc_ol(20,16) = PA(3)*pc_ol(7,16) + inv_2zeta*(2*pc_ol(4,16) + 2*pc_ol(7,9))

pc_ol(1,17) = PB(2)*pc_ol(1,7)

pc_ol(2,17) = PA(1)*pc_ol(1,17)

pc_ol(3,17) = PA(2)*pc_ol(1,17) + inv_2zeta*(pc_ol(1,7))

pc_ol(4,17) = PA(3)*pc_ol(1,17) + inv_2zeta*(2*pc_ol(1,10))

pc_ol(5,17) = PA(1)*pc_ol(2,17) + inv_2zeta*(pc_ol(1,17))

pc_ol(6,17) = PA(2)*pc_ol(3,17) + inv_2zeta*(pc_ol(1,17) + pc_ol(3,7))

pc_ol(7,17) = PA(3)*pc_ol(4,17) + inv_2zeta*(pc_ol(1,17) + 2*pc_ol(4,10))

pc_ol(8,17) = PA(2)*pc_ol(2,17) + inv_2zeta*(pc_ol(2,7))

pc_ol(9,17) = PA(3)*pc_ol(2,17) + inv_2zeta*(2*pc_ol(2,10))

pc_ol(10,17) = PA(3)*pc_ol(3,17) + inv_2zeta*(2*pc_ol(3,10))

pc_ol(11,17) = PA(3)*pc_ol(8,17) + inv_2zeta*(2*pc_ol(8,10))

pc_ol(12,17) = PA(2)*pc_ol(5,17) + inv_2zeta*(pc_ol(5,7))

pc_ol(13,17) = PA(3)*pc_ol(5,17) + inv_2zeta*(2*pc_ol(5,10))

pc_ol(14,17) = PA(1)*pc_ol(6,17)

pc_ol(15,17) = PA(3)*pc_ol(6,17) + inv_2zeta*(2*pc_ol(6,10))

pc_ol(16,17) = PA(1)*pc_ol(7,17)

pc_ol(17,17) = PA(2)*pc_ol(7,17) + inv_2zeta*(pc_ol(7,7))

pc_ol(18,17) = PA(1)*pc_ol(5,17) + inv_2zeta*(2*pc_ol(2,17))

pc_ol(19,17) = PA(2)*pc_ol(6,17) + inv_2zeta*(2*pc_ol(3,17) + pc_ol(6,7))

pc_ol(20,17) = PA(3)*pc_ol(7,17) + inv_2zeta*(2*pc_ol(4,17) + 2*pc_ol(7,10))

pc_ol(1,18) = PB(1)*pc_ol(1,5) + inv_2zeta*(2*pc_ol(1,2))

pc_ol(2,18) = PA(1)*pc_ol(1,18) + inv_2zeta*(3*pc_ol(1,5))

pc_ol(3,18) = PA(2)*pc_ol(1,18)

pc_ol(4,18) = PA(3)*pc_ol(1,18)

pc_ol(5,18) = PA(1)*pc_ol(2,18) + inv_2zeta*(pc_ol(1,18) + 3*pc_ol(2,5))

pc_ol(6,18) = PA(2)*pc_ol(3,18) + inv_2zeta*(pc_ol(1,18))

pc_ol(7,18) = PA(3)*pc_ol(4,18) + inv_2zeta*(pc_ol(1,18))

pc_ol(8,18) = PA(2)*pc_ol(2,18)

pc_ol(9,18) = PA(3)*pc_ol(2,18)

pc_ol(10,18) = PA(3)*pc_ol(3,18)

pc_ol(11,18) = PA(3)*pc_ol(8,18)

pc_ol(12,18) = PA(2)*pc_ol(5,18)

pc_ol(13,18) = PA(3)*pc_ol(5,18)

pc_ol(14,18) = PA(1)*pc_ol(6,18) + inv_2zeta*(3*pc_ol(6,5))

pc_ol(15,18) = PA(3)*pc_ol(6,18)

pc_ol(16,18) = PA(1)*pc_ol(7,18) + inv_2zeta*(3*pc_ol(7,5))

pc_ol(17,18) = PA(2)*pc_ol(7,18)

pc_ol(18,18) = PA(1)*pc_ol(5,18) + inv_2zeta*(2*pc_ol(2,18) + 3*pc_ol(5,5))

pc_ol(19,18) = PA(2)*pc_ol(6,18) + inv_2zeta*(2*pc_ol(3,18))

pc_ol(20,18) = PA(3)*pc_ol(7,18) + inv_2zeta*(2*pc_ol(4,18))

pc_ol(1,19) = PB(2)*pc_ol(1,6) + inv_2zeta*(2*pc_ol(1,3))

pc_ol(2,19) = PA(1)*pc_ol(1,19)

pc_ol(3,19) = PA(2)*pc_ol(1,19) + inv_2zeta*(3*pc_ol(1,6))

pc_ol(4,19) = PA(3)*pc_ol(1,19)

pc_ol(5,19) = PA(1)*pc_ol(2,19) + inv_2zeta*(pc_ol(1,19))

pc_ol(6,19) = PA(2)*pc_ol(3,19) + inv_2zeta*(pc_ol(1,19) + 3*pc_ol(3,6))

pc_ol(7,19) = PA(3)*pc_ol(4,19) + inv_2zeta*(pc_ol(1,19))

pc_ol(8,19) = PA(2)*pc_ol(2,19) + inv_2zeta*(3*pc_ol(2,6))

pc_ol(9,19) = PA(3)*pc_ol(2,19)

pc_ol(10,19) = PA(3)*pc_ol(3,19)

pc_ol(11,19) = PA(3)*pc_ol(8,19)

pc_ol(12,19) = PA(2)*pc_ol(5,19) + inv_2zeta*(3*pc_ol(5,6))

pc_ol(13,19) = PA(3)*pc_ol(5,19)

pc_ol(14,19) = PA(1)*pc_ol(6,19)

pc_ol(15,19) = PA(3)*pc_ol(6,19)

pc_ol(16,19) = PA(1)*pc_ol(7,19)

pc_ol(17,19) = PA(2)*pc_ol(7,19) + inv_2zeta*(3*pc_ol(7,6))

pc_ol(18,19) = PA(1)*pc_ol(5,19) + inv_2zeta*(2*pc_ol(2,19))

pc_ol(19,19) = PA(2)*pc_ol(6,19) + inv_2zeta*(2*pc_ol(3,19) + 3*pc_ol(6,6))

pc_ol(20,19) = PA(3)*pc_ol(7,19) + inv_2zeta*(2*pc_ol(4,19))

pc_ol(1,20) = PB(3)*pc_ol(1,7) + inv_2zeta*(2*pc_ol(1,4))

pc_ol(2,20) = PA(1)*pc_ol(1,20)

pc_ol(3,20) = PA(2)*pc_ol(1,20)

pc_ol(4,20) = PA(3)*pc_ol(1,20) + inv_2zeta*(3*pc_ol(1,7))

pc_ol(5,20) = PA(1)*pc_ol(2,20) + inv_2zeta*(pc_ol(1,20))

pc_ol(6,20) = PA(2)*pc_ol(3,20) + inv_2zeta*(pc_ol(1,20))

pc_ol(7,20) = PA(3)*pc_ol(4,20) + inv_2zeta*(pc_ol(1,20) + 3*pc_ol(4,7))

pc_ol(8,20) = PA(2)*pc_ol(2,20)

pc_ol(9,20) = PA(3)*pc_ol(2,20) + inv_2zeta*(3*pc_ol(2,7))

pc_ol(10,20) = PA(3)*pc_ol(3,20) + inv_2zeta*(3*pc_ol(3,7))

pc_ol(11,20) = PA(3)*pc_ol(8,20) + inv_2zeta*(3*pc_ol(8,7))

pc_ol(12,20) = PA(2)*pc_ol(5,20)

pc_ol(13,20) = PA(3)*pc_ol(5,20) + inv_2zeta*(3*pc_ol(5,7))

pc_ol(14,20) = PA(1)*pc_ol(6,20)

pc_ol(15,20) = PA(3)*pc_ol(6,20) + inv_2zeta*(3*pc_ol(6,7))

pc_ol(16,20) = PA(1)*pc_ol(7,20)

pc_ol(17,20) = PA(2)*pc_ol(7,20)

pc_ol(18,20) = PA(1)*pc_ol(5,20) + inv_2zeta*(2*pc_ol(2,20))

pc_ol(19,20) = PA(2)*pc_ol(6,20) + inv_2zeta*(2*pc_ol(3,20))

pc_ol(20,20) = PA(3)*pc_ol(7,20) + inv_2zeta*(2*pc_ol(4,20) + 3*pc_ol(7,7))

pc(1,1) = preFactorKE

pc(2,1) = PA(1)*preFactorKE + 2*xi*(pc_ol(2,1))

pc(3,1) = PA(2)*preFactorKE + 2*xi*(pc_ol(3,1))

pc(4,1) = PA(3)*preFactorKE + 2*xi*(pc_ol(4,1))

pc(5,1) = PA(1)*pc(2,1) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(5,1) - inv_2ze&
&ta_a*preFactorOL)

pc(6,1) = PA(2)*pc(3,1) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(6,1) - inv_2ze&
&ta_a*preFactorOL)

pc(7,1) = PA(3)*pc(4,1) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(7,1) - inv_2ze&
&ta_a*preFactorOL)

pc(8,1) = PA(2)*pc(2,1) + 2*xi*(pc_ol(8,1))

pc(9,1) = PA(3)*pc(2,1) + 2*xi*(pc_ol(9,1))

pc(10,1) = PA(3)*pc(3,1) + 2*xi*(pc_ol(10,1))

pc(11,1) = PA(3)*pc(8,1) + 2*xi*(pc_ol(11,1))

pc(12,1) = PA(2)*pc(5,1) + 2*xi*(pc_ol(12,1))

pc(13,1) = PA(3)*pc(5,1) + 2*xi*(pc_ol(13,1))

pc(14,1) = PA(1)*pc(6,1) + 2*xi*(pc_ol(14,1))

pc(15,1) = PA(3)*pc(6,1) + 2*xi*(pc_ol(15,1))

pc(16,1) = PA(1)*pc(7,1) + 2*xi*(pc_ol(16,1))

pc(17,1) = PA(2)*pc(7,1) + 2*xi*(pc_ol(17,1))

pc(18,1) = PA(1)*pc(5,1) + inv_2zeta*(2*pc(2,1)) + 2*xi*(pc_ol(18,1) - inv_2ze&
&ta_a*2*pc_ol(2,1))

pc(19,1) = PA(2)*pc(6,1) + inv_2zeta*(2*pc(3,1)) + 2*xi*(pc_ol(19,1) - inv_2ze&
&ta_a*2*pc_ol(3,1))

pc(20,1) = PA(3)*pc(7,1) + inv_2zeta*(2*pc(4,1)) + 2*xi*(pc_ol(20,1) - inv_2ze&
&ta_a*2*pc_ol(4,1))

pc(1,2) = PB(1)*preFactorKE + 2*xi*(pc_ol(1,2))

pc(2,2) = PA(1)*pc(1,2) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(2,2))

pc(3,2) = PA(2)*pc(1,2) + 2*xi*(pc_ol(3,2))

pc(4,2) = PA(3)*pc(1,2) + 2*xi*(pc_ol(4,2))

pc(5,2) = PA(1)*pc(2,2) + inv_2zeta*(pc(1,2) + pc(2,1)) + 2*xi*(pc_ol(5,2) - i&
&nv_2zeta_a*pc_ol(1,2))

pc(6,2) = PA(2)*pc(3,2) + inv_2zeta*(pc(1,2)) + 2*xi*(pc_ol(6,2) - inv_2zeta_a&
&*pc_ol(1,2))

pc(7,2) = PA(3)*pc(4,2) + inv_2zeta*(pc(1,2)) + 2*xi*(pc_ol(7,2) - inv_2zeta_a&
&*pc_ol(1,2))

pc(8,2) = PA(2)*pc(2,2) + 2*xi*(pc_ol(8,2))

pc(9,2) = PA(3)*pc(2,2) + 2*xi*(pc_ol(9,2))

pc(10,2) = PA(3)*pc(3,2) + 2*xi*(pc_ol(10,2))

pc(11,2) = PA(3)*pc(8,2) + 2*xi*(pc_ol(11,2))

pc(12,2) = PA(2)*pc(5,2) + 2*xi*(pc_ol(12,2))

pc(13,2) = PA(3)*pc(5,2) + 2*xi*(pc_ol(13,2))

pc(14,2) = PA(1)*pc(6,2) + inv_2zeta*(pc(6,1)) + 2*xi*(pc_ol(14,2))

pc(15,2) = PA(3)*pc(6,2) + 2*xi*(pc_ol(15,2))

pc(16,2) = PA(1)*pc(7,2) + inv_2zeta*(pc(7,1)) + 2*xi*(pc_ol(16,2))

pc(17,2) = PA(2)*pc(7,2) + 2*xi*(pc_ol(17,2))

pc(18,2) = PA(1)*pc(5,2) + inv_2zeta*(2*pc(2,2) + pc(5,1)) + 2*xi*(pc_ol(18,2)&
& - inv_2zeta_a*2*pc_ol(2,2))

pc(19,2) = PA(2)*pc(6,2) + inv_2zeta*(2*pc(3,2)) + 2*xi*(pc_ol(19,2) - inv_2ze&
&ta_a*2*pc_ol(3,2))

pc(20,2) = PA(3)*pc(7,2) + inv_2zeta*(2*pc(4,2)) + 2*xi*(pc_ol(20,2) - inv_2ze&
&ta_a*2*pc_ol(4,2))

pc(1,3) = PB(2)*preFactorKE + 2*xi*(pc_ol(1,3))

pc(2,3) = PA(1)*pc(1,3) + 2*xi*(pc_ol(2,3))

pc(3,3) = PA(2)*pc(1,3) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(3,3))

pc(4,3) = PA(3)*pc(1,3) + 2*xi*(pc_ol(4,3))

pc(5,3) = PA(1)*pc(2,3) + inv_2zeta*(pc(1,3)) + 2*xi*(pc_ol(5,3) - inv_2zeta_a&
&*pc_ol(1,3))

pc(6,3) = PA(2)*pc(3,3) + inv_2zeta*(pc(1,3) + pc(3,1)) + 2*xi*(pc_ol(6,3) - i&
&nv_2zeta_a*pc_ol(1,3))

pc(7,3) = PA(3)*pc(4,3) + inv_2zeta*(pc(1,3)) + 2*xi*(pc_ol(7,3) - inv_2zeta_a&
&*pc_ol(1,3))

pc(8,3) = PA(2)*pc(2,3) + inv_2zeta*(pc(2,1)) + 2*xi*(pc_ol(8,3))

pc(9,3) = PA(3)*pc(2,3) + 2*xi*(pc_ol(9,3))

pc(10,3) = PA(3)*pc(3,3) + 2*xi*(pc_ol(10,3))

pc(11,3) = PA(3)*pc(8,3) + 2*xi*(pc_ol(11,3))

pc(12,3) = PA(2)*pc(5,3) + inv_2zeta*(pc(5,1)) + 2*xi*(pc_ol(12,3))

pc(13,3) = PA(3)*pc(5,3) + 2*xi*(pc_ol(13,3))

pc(14,3) = PA(1)*pc(6,3) + 2*xi*(pc_ol(14,3))

pc(15,3) = PA(3)*pc(6,3) + 2*xi*(pc_ol(15,3))

pc(16,3) = PA(1)*pc(7,3) + 2*xi*(pc_ol(16,3))

pc(17,3) = PA(2)*pc(7,3) + inv_2zeta*(pc(7,1)) + 2*xi*(pc_ol(17,3))

pc(18,3) = PA(1)*pc(5,3) + inv_2zeta*(2*pc(2,3)) + 2*xi*(pc_ol(18,3) - inv_2ze&
&ta_a*2*pc_ol(2,3))

pc(19,3) = PA(2)*pc(6,3) + inv_2zeta*(2*pc(3,3) + pc(6,1)) + 2*xi*(pc_ol(19,3)&
& - inv_2zeta_a*2*pc_ol(3,3))

pc(20,3) = PA(3)*pc(7,3) + inv_2zeta*(2*pc(4,3)) + 2*xi*(pc_ol(20,3) - inv_2ze&
&ta_a*2*pc_ol(4,3))

pc(1,4) = PB(3)*preFactorKE + 2*xi*(pc_ol(1,4))

pc(2,4) = PA(1)*pc(1,4) + 2*xi*(pc_ol(2,4))

pc(3,4) = PA(2)*pc(1,4) + 2*xi*(pc_ol(3,4))

pc(4,4) = PA(3)*pc(1,4) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(4,4))

pc(5,4) = PA(1)*pc(2,4) + inv_2zeta*(pc(1,4)) + 2*xi*(pc_ol(5,4) - inv_2zeta_a&
&*pc_ol(1,4))

pc(6,4) = PA(2)*pc(3,4) + inv_2zeta*(pc(1,4)) + 2*xi*(pc_ol(6,4) - inv_2zeta_a&
&*pc_ol(1,4))

pc(7,4) = PA(3)*pc(4,4) + inv_2zeta*(pc(1,4) + pc(4,1)) + 2*xi*(pc_ol(7,4) - i&
&nv_2zeta_a*pc_ol(1,4))

pc(8,4) = PA(2)*pc(2,4) + 2*xi*(pc_ol(8,4))

pc(9,4) = PA(3)*pc(2,4) + inv_2zeta*(pc(2,1)) + 2*xi*(pc_ol(9,4))

pc(10,4) = PA(3)*pc(3,4) + inv_2zeta*(pc(3,1)) + 2*xi*(pc_ol(10,4))

pc(11,4) = PA(3)*pc(8,4) + inv_2zeta*(pc(8,1)) + 2*xi*(pc_ol(11,4))

pc(12,4) = PA(2)*pc(5,4) + 2*xi*(pc_ol(12,4))

pc(13,4) = PA(3)*pc(5,4) + inv_2zeta*(pc(5,1)) + 2*xi*(pc_ol(13,4))

pc(14,4) = PA(1)*pc(6,4) + 2*xi*(pc_ol(14,4))

pc(15,4) = PA(3)*pc(6,4) + inv_2zeta*(pc(6,1)) + 2*xi*(pc_ol(15,4))

pc(16,4) = PA(1)*pc(7,4) + 2*xi*(pc_ol(16,4))

pc(17,4) = PA(2)*pc(7,4) + 2*xi*(pc_ol(17,4))

pc(18,4) = PA(1)*pc(5,4) + inv_2zeta*(2*pc(2,4)) + 2*xi*(pc_ol(18,4) - inv_2ze&
&ta_a*2*pc_ol(2,4))

pc(19,4) = PA(2)*pc(6,4) + inv_2zeta*(2*pc(3,4)) + 2*xi*(pc_ol(19,4) - inv_2ze&
&ta_a*2*pc_ol(3,4))

pc(20,4) = PA(3)*pc(7,4) + inv_2zeta*(2*pc(4,4) + pc(7,1)) + 2*xi*(pc_ol(20,4)&
& - inv_2zeta_a*2*pc_ol(4,4))

pc(1,5) = PB(1)*pc(1,2) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(1,5) - inv_2ze&
&ta_b*preFactorOL)

pc(2,5) = PA(1)*pc(1,5) + inv_2zeta*(2*pc(1,2)) + 2*xi*(pc_ol(2,5))

pc(3,5) = PA(2)*pc(1,5) + 2*xi*(pc_ol(3,5))

pc(4,5) = PA(3)*pc(1,5) + 2*xi*(pc_ol(4,5))

pc(5,5) = PA(1)*pc(2,5) + inv_2zeta*(pc(1,5) + 2*pc(2,2)) + 2*xi*(pc_ol(5,5) -&
& inv_2zeta_a*pc_ol(1,5))

pc(6,5) = PA(2)*pc(3,5) + inv_2zeta*(pc(1,5)) + 2*xi*(pc_ol(6,5) - inv_2zeta_a&
&*pc_ol(1,5))

pc(7,5) = PA(3)*pc(4,5) + inv_2zeta*(pc(1,5)) + 2*xi*(pc_ol(7,5) - inv_2zeta_a&
&*pc_ol(1,5))

pc(8,5) = PA(2)*pc(2,5) + 2*xi*(pc_ol(8,5))

pc(9,5) = PA(3)*pc(2,5) + 2*xi*(pc_ol(9,5))

pc(10,5) = PA(3)*pc(3,5) + 2*xi*(pc_ol(10,5))

pc(11,5) = PA(3)*pc(8,5) + 2*xi*(pc_ol(11,5))

pc(12,5) = PA(2)*pc(5,5) + 2*xi*(pc_ol(12,5))

pc(13,5) = PA(3)*pc(5,5) + 2*xi*(pc_ol(13,5))

pc(14,5) = PA(1)*pc(6,5) + inv_2zeta*(2*pc(6,2)) + 2*xi*(pc_ol(14,5))

pc(15,5) = PA(3)*pc(6,5) + 2*xi*(pc_ol(15,5))

pc(16,5) = PA(1)*pc(7,5) + inv_2zeta*(2*pc(7,2)) + 2*xi*(pc_ol(16,5))

pc(17,5) = PA(2)*pc(7,5) + 2*xi*(pc_ol(17,5))

pc(18,5) = PA(1)*pc(5,5) + inv_2zeta*(2*pc(2,5) + 2*pc(5,2)) + 2*xi*(pc_ol(18,&
&5) - inv_2zeta_a*2*pc_ol(2,5))

pc(19,5) = PA(2)*pc(6,5) + inv_2zeta*(2*pc(3,5)) + 2*xi*(pc_ol(19,5) - inv_2ze&
&ta_a*2*pc_ol(3,5))

pc(20,5) = PA(3)*pc(7,5) + inv_2zeta*(2*pc(4,5)) + 2*xi*(pc_ol(20,5) - inv_2ze&
&ta_a*2*pc_ol(4,5))

pc(1,6) = PB(2)*pc(1,3) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(1,6) - inv_2ze&
&ta_b*preFactorOL)

pc(2,6) = PA(1)*pc(1,6) + 2*xi*(pc_ol(2,6))

pc(3,6) = PA(2)*pc(1,6) + inv_2zeta*(2*pc(1,3)) + 2*xi*(pc_ol(3,6))

pc(4,6) = PA(3)*pc(1,6) + 2*xi*(pc_ol(4,6))

pc(5,6) = PA(1)*pc(2,6) + inv_2zeta*(pc(1,6)) + 2*xi*(pc_ol(5,6) - inv_2zeta_a&
&*pc_ol(1,6))

pc(6,6) = PA(2)*pc(3,6) + inv_2zeta*(pc(1,6) + 2*pc(3,3)) + 2*xi*(pc_ol(6,6) -&
& inv_2zeta_a*pc_ol(1,6))

pc(7,6) = PA(3)*pc(4,6) + inv_2zeta*(pc(1,6)) + 2*xi*(pc_ol(7,6) - inv_2zeta_a&
&*pc_ol(1,6))

pc(8,6) = PA(2)*pc(2,6) + inv_2zeta*(2*pc(2,3)) + 2*xi*(pc_ol(8,6))

pc(9,6) = PA(3)*pc(2,6) + 2*xi*(pc_ol(9,6))

pc(10,6) = PA(3)*pc(3,6) + 2*xi*(pc_ol(10,6))

pc(11,6) = PA(3)*pc(8,6) + 2*xi*(pc_ol(11,6))

pc(12,6) = PA(2)*pc(5,6) + inv_2zeta*(2*pc(5,3)) + 2*xi*(pc_ol(12,6))

pc(13,6) = PA(3)*pc(5,6) + 2*xi*(pc_ol(13,6))

pc(14,6) = PA(1)*pc(6,6) + 2*xi*(pc_ol(14,6))

pc(15,6) = PA(3)*pc(6,6) + 2*xi*(pc_ol(15,6))

pc(16,6) = PA(1)*pc(7,6) + 2*xi*(pc_ol(16,6))

pc(17,6) = PA(2)*pc(7,6) + inv_2zeta*(2*pc(7,3)) + 2*xi*(pc_ol(17,6))

pc(18,6) = PA(1)*pc(5,6) + inv_2zeta*(2*pc(2,6)) + 2*xi*(pc_ol(18,6) - inv_2ze&
&ta_a*2*pc_ol(2,6))

pc(19,6) = PA(2)*pc(6,6) + inv_2zeta*(2*pc(3,6) + 2*pc(6,3)) + 2*xi*(pc_ol(19,&
&6) - inv_2zeta_a*2*pc_ol(3,6))

pc(20,6) = PA(3)*pc(7,6) + inv_2zeta*(2*pc(4,6)) + 2*xi*(pc_ol(20,6) - inv_2ze&
&ta_a*2*pc_ol(4,6))

pc(1,7) = PB(3)*pc(1,4) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(1,7) - inv_2ze&
&ta_b*preFactorOL)

pc(2,7) = PA(1)*pc(1,7) + 2*xi*(pc_ol(2,7))

pc(3,7) = PA(2)*pc(1,7) + 2*xi*(pc_ol(3,7))

pc(4,7) = PA(3)*pc(1,7) + inv_2zeta*(2*pc(1,4)) + 2*xi*(pc_ol(4,7))

pc(5,7) = PA(1)*pc(2,7) + inv_2zeta*(pc(1,7)) + 2*xi*(pc_ol(5,7) - inv_2zeta_a&
&*pc_ol(1,7))

pc(6,7) = PA(2)*pc(3,7) + inv_2zeta*(pc(1,7)) + 2*xi*(pc_ol(6,7) - inv_2zeta_a&
&*pc_ol(1,7))

pc(7,7) = PA(3)*pc(4,7) + inv_2zeta*(pc(1,7) + 2*pc(4,4)) + 2*xi*(pc_ol(7,7) -&
& inv_2zeta_a*pc_ol(1,7))

pc(8,7) = PA(2)*pc(2,7) + 2*xi*(pc_ol(8,7))

pc(9,7) = PA(3)*pc(2,7) + inv_2zeta*(2*pc(2,4)) + 2*xi*(pc_ol(9,7))

pc(10,7) = PA(3)*pc(3,7) + inv_2zeta*(2*pc(3,4)) + 2*xi*(pc_ol(10,7))

pc(11,7) = PA(3)*pc(8,7) + inv_2zeta*(2*pc(8,4)) + 2*xi*(pc_ol(11,7))

pc(12,7) = PA(2)*pc(5,7) + 2*xi*(pc_ol(12,7))

pc(13,7) = PA(3)*pc(5,7) + inv_2zeta*(2*pc(5,4)) + 2*xi*(pc_ol(13,7))

pc(14,7) = PA(1)*pc(6,7) + 2*xi*(pc_ol(14,7))

pc(15,7) = PA(3)*pc(6,7) + inv_2zeta*(2*pc(6,4)) + 2*xi*(pc_ol(15,7))

pc(16,7) = PA(1)*pc(7,7) + 2*xi*(pc_ol(16,7))

pc(17,7) = PA(2)*pc(7,7) + 2*xi*(pc_ol(17,7))

pc(18,7) = PA(1)*pc(5,7) + inv_2zeta*(2*pc(2,7)) + 2*xi*(pc_ol(18,7) - inv_2ze&
&ta_a*2*pc_ol(2,7))

pc(19,7) = PA(2)*pc(6,7) + inv_2zeta*(2*pc(3,7)) + 2*xi*(pc_ol(19,7) - inv_2ze&
&ta_a*2*pc_ol(3,7))

pc(20,7) = PA(3)*pc(7,7) + inv_2zeta*(2*pc(4,7) + 2*pc(7,4)) + 2*xi*(pc_ol(20,&
&7) - inv_2zeta_a*2*pc_ol(4,7))

pc(1,8) = PB(2)*pc(1,2) + 2*xi*(pc_ol(1,8))

pc(2,8) = PA(1)*pc(1,8) + inv_2zeta*(pc(1,3)) + 2*xi*(pc_ol(2,8))

pc(3,8) = PA(2)*pc(1,8) + inv_2zeta*(pc(1,2)) + 2*xi*(pc_ol(3,8))

pc(4,8) = PA(3)*pc(1,8) + 2*xi*(pc_ol(4,8))

pc(5,8) = PA(1)*pc(2,8) + inv_2zeta*(pc(1,8) + pc(2,3)) + 2*xi*(pc_ol(5,8) - i&
&nv_2zeta_a*pc_ol(1,8))

pc(6,8) = PA(2)*pc(3,8) + inv_2zeta*(pc(1,8) + pc(3,2)) + 2*xi*(pc_ol(6,8) - i&
&nv_2zeta_a*pc_ol(1,8))

pc(7,8) = PA(3)*pc(4,8) + inv_2zeta*(pc(1,8)) + 2*xi*(pc_ol(7,8) - inv_2zeta_a&
&*pc_ol(1,8))

pc(8,8) = PA(2)*pc(2,8) + inv_2zeta*(pc(2,2)) + 2*xi*(pc_ol(8,8))

pc(9,8) = PA(3)*pc(2,8) + 2*xi*(pc_ol(9,8))

pc(10,8) = PA(3)*pc(3,8) + 2*xi*(pc_ol(10,8))

pc(11,8) = PA(3)*pc(8,8) + 2*xi*(pc_ol(11,8))

pc(12,8) = PA(2)*pc(5,8) + inv_2zeta*(pc(5,2)) + 2*xi*(pc_ol(12,8))

pc(13,8) = PA(3)*pc(5,8) + 2*xi*(pc_ol(13,8))

pc(14,8) = PA(1)*pc(6,8) + inv_2zeta*(pc(6,3)) + 2*xi*(pc_ol(14,8))

pc(15,8) = PA(3)*pc(6,8) + 2*xi*(pc_ol(15,8))

pc(16,8) = PA(1)*pc(7,8) + inv_2zeta*(pc(7,3)) + 2*xi*(pc_ol(16,8))

pc(17,8) = PA(2)*pc(7,8) + inv_2zeta*(pc(7,2)) + 2*xi*(pc_ol(17,8))

pc(18,8) = PA(1)*pc(5,8) + inv_2zeta*(2*pc(2,8) + pc(5,3)) + 2*xi*(pc_ol(18,8)&
& - inv_2zeta_a*2*pc_ol(2,8))

pc(19,8) = PA(2)*pc(6,8) + inv_2zeta*(2*pc(3,8) + pc(6,2)) + 2*xi*(pc_ol(19,8)&
& - inv_2zeta_a*2*pc_ol(3,8))

pc(20,8) = PA(3)*pc(7,8) + inv_2zeta*(2*pc(4,8)) + 2*xi*(pc_ol(20,8) - inv_2ze&
&ta_a*2*pc_ol(4,8))

pc(1,9) = PB(3)*pc(1,2) + 2*xi*(pc_ol(1,9))

pc(2,9) = PA(1)*pc(1,9) + inv_2zeta*(pc(1,4)) + 2*xi*(pc_ol(2,9))

pc(3,9) = PA(2)*pc(1,9) + 2*xi*(pc_ol(3,9))

pc(4,9) = PA(3)*pc(1,9) + inv_2zeta*(pc(1,2)) + 2*xi*(pc_ol(4,9))

pc(5,9) = PA(1)*pc(2,9) + inv_2zeta*(pc(1,9) + pc(2,4)) + 2*xi*(pc_ol(5,9) - i&
&nv_2zeta_a*pc_ol(1,9))

pc(6,9) = PA(2)*pc(3,9) + inv_2zeta*(pc(1,9)) + 2*xi*(pc_ol(6,9) - inv_2zeta_a&
&*pc_ol(1,9))

pc(7,9) = PA(3)*pc(4,9) + inv_2zeta*(pc(1,9) + pc(4,2)) + 2*xi*(pc_ol(7,9) - i&
&nv_2zeta_a*pc_ol(1,9))

pc(8,9) = PA(2)*pc(2,9) + 2*xi*(pc_ol(8,9))

pc(9,9) = PA(3)*pc(2,9) + inv_2zeta*(pc(2,2)) + 2*xi*(pc_ol(9,9))

pc(10,9) = PA(3)*pc(3,9) + inv_2zeta*(pc(3,2)) + 2*xi*(pc_ol(10,9))

pc(11,9) = PA(3)*pc(8,9) + inv_2zeta*(pc(8,2)) + 2*xi*(pc_ol(11,9))

pc(12,9) = PA(2)*pc(5,9) + 2*xi*(pc_ol(12,9))

pc(13,9) = PA(3)*pc(5,9) + inv_2zeta*(pc(5,2)) + 2*xi*(pc_ol(13,9))

pc(14,9) = PA(1)*pc(6,9) + inv_2zeta*(pc(6,4)) + 2*xi*(pc_ol(14,9))

pc(15,9) = PA(3)*pc(6,9) + inv_2zeta*(pc(6,2)) + 2*xi*(pc_ol(15,9))

pc(16,9) = PA(1)*pc(7,9) + inv_2zeta*(pc(7,4)) + 2*xi*(pc_ol(16,9))

pc(17,9) = PA(2)*pc(7,9) + 2*xi*(pc_ol(17,9))

pc(18,9) = PA(1)*pc(5,9) + inv_2zeta*(2*pc(2,9) + pc(5,4)) + 2*xi*(pc_ol(18,9)&
& - inv_2zeta_a*2*pc_ol(2,9))

pc(19,9) = PA(2)*pc(6,9) + inv_2zeta*(2*pc(3,9)) + 2*xi*(pc_ol(19,9) - inv_2ze&
&ta_a*2*pc_ol(3,9))

pc(20,9) = PA(3)*pc(7,9) + inv_2zeta*(2*pc(4,9) + pc(7,2)) + 2*xi*(pc_ol(20,9)&
& - inv_2zeta_a*2*pc_ol(4,9))

pc(1,10) = PB(3)*pc(1,3) + 2*xi*(pc_ol(1,10))

pc(2,10) = PA(1)*pc(1,10) + 2*xi*(pc_ol(2,10))

pc(3,10) = PA(2)*pc(1,10) + inv_2zeta*(pc(1,4)) + 2*xi*(pc_ol(3,10))

pc(4,10) = PA(3)*pc(1,10) + inv_2zeta*(pc(1,3)) + 2*xi*(pc_ol(4,10))

pc(5,10) = PA(1)*pc(2,10) + inv_2zeta*(pc(1,10)) + 2*xi*(pc_ol(5,10) - inv_2ze&
&ta_a*pc_ol(1,10))

pc(6,10) = PA(2)*pc(3,10) + inv_2zeta*(pc(1,10) + pc(3,4)) + 2*xi*(pc_ol(6,10)&
& - inv_2zeta_a*pc_ol(1,10))

pc(7,10) = PA(3)*pc(4,10) + inv_2zeta*(pc(1,10) + pc(4,3)) + 2*xi*(pc_ol(7,10)&
& - inv_2zeta_a*pc_ol(1,10))

pc(8,10) = PA(2)*pc(2,10) + inv_2zeta*(pc(2,4)) + 2*xi*(pc_ol(8,10))

pc(9,10) = PA(3)*pc(2,10) + inv_2zeta*(pc(2,3)) + 2*xi*(pc_ol(9,10))

pc(10,10) = PA(3)*pc(3,10) + inv_2zeta*(pc(3,3)) + 2*xi*(pc_ol(10,10))

pc(11,10) = PA(3)*pc(8,10) + inv_2zeta*(pc(8,3)) + 2*xi*(pc_ol(11,10))

pc(12,10) = PA(2)*pc(5,10) + inv_2zeta*(pc(5,4)) + 2*xi*(pc_ol(12,10))

pc(13,10) = PA(3)*pc(5,10) + inv_2zeta*(pc(5,3)) + 2*xi*(pc_ol(13,10))

pc(14,10) = PA(1)*pc(6,10) + 2*xi*(pc_ol(14,10))

pc(15,10) = PA(3)*pc(6,10) + inv_2zeta*(pc(6,3)) + 2*xi*(pc_ol(15,10))

pc(16,10) = PA(1)*pc(7,10) + 2*xi*(pc_ol(16,10))

pc(17,10) = PA(2)*pc(7,10) + inv_2zeta*(pc(7,4)) + 2*xi*(pc_ol(17,10))

pc(18,10) = PA(1)*pc(5,10) + inv_2zeta*(2*pc(2,10)) + 2*xi*(pc_ol(18,10) - inv&
&_2zeta_a*2*pc_ol(2,10))

pc(19,10) = PA(2)*pc(6,10) + inv_2zeta*(2*pc(3,10) + pc(6,4)) + 2*xi*(pc_ol(19&
&,10) - inv_2zeta_a*2*pc_ol(3,10))

pc(20,10) = PA(3)*pc(7,10) + inv_2zeta*(2*pc(4,10) + pc(7,3)) + 2*xi*(pc_ol(20&
&,10) - inv_2zeta_a*2*pc_ol(4,10))

pc(1,11) = PB(3)*pc(1,8) + 2*xi*(pc_ol(1,11))

pc(2,11) = PA(1)*pc(1,11) + inv_2zeta*(pc(1,10)) + 2*xi*(pc_ol(2,11))

pc(3,11) = PA(2)*pc(1,11) + inv_2zeta*(pc(1,9)) + 2*xi*(pc_ol(3,11))

pc(4,11) = PA(3)*pc(1,11) + inv_2zeta*(pc(1,8)) + 2*xi*(pc_ol(4,11))

pc(5,11) = PA(1)*pc(2,11) + inv_2zeta*(pc(1,11) + pc(2,10)) + 2*xi*(pc_ol(5,11&
&) - inv_2zeta_a*pc_ol(1,11))

pc(6,11) = PA(2)*pc(3,11) + inv_2zeta*(pc(1,11) + pc(3,9)) + 2*xi*(pc_ol(6,11)&
& - inv_2zeta_a*pc_ol(1,11))

pc(7,11) = PA(3)*pc(4,11) + inv_2zeta*(pc(1,11) + pc(4,8)) + 2*xi*(pc_ol(7,11)&
& - inv_2zeta_a*pc_ol(1,11))

pc(8,11) = PA(2)*pc(2,11) + inv_2zeta*(pc(2,9)) + 2*xi*(pc_ol(8,11))

pc(9,11) = PA(3)*pc(2,11) + inv_2zeta*(pc(2,8)) + 2*xi*(pc_ol(9,11))

pc(10,11) = PA(3)*pc(3,11) + inv_2zeta*(pc(3,8)) + 2*xi*(pc_ol(10,11))

pc(11,11) = PA(3)*pc(8,11) + inv_2zeta*(pc(8,8)) + 2*xi*(pc_ol(11,11))

pc(12,11) = PA(2)*pc(5,11) + inv_2zeta*(pc(5,9)) + 2*xi*(pc_ol(12,11))

pc(13,11) = PA(3)*pc(5,11) + inv_2zeta*(pc(5,8)) + 2*xi*(pc_ol(13,11))

pc(14,11) = PA(1)*pc(6,11) + inv_2zeta*(pc(6,10)) + 2*xi*(pc_ol(14,11))

pc(15,11) = PA(3)*pc(6,11) + inv_2zeta*(pc(6,8)) + 2*xi*(pc_ol(15,11))

pc(16,11) = PA(1)*pc(7,11) + inv_2zeta*(pc(7,10)) + 2*xi*(pc_ol(16,11))

pc(17,11) = PA(2)*pc(7,11) + inv_2zeta*(pc(7,9)) + 2*xi*(pc_ol(17,11))

pc(18,11) = PA(1)*pc(5,11) + inv_2zeta*(2*pc(2,11) + pc(5,10)) + 2*xi*(pc_ol(1&
&8,11) - inv_2zeta_a*2*pc_ol(2,11))

pc(19,11) = PA(2)*pc(6,11) + inv_2zeta*(2*pc(3,11) + pc(6,9)) + 2*xi*(pc_ol(19&
&,11) - inv_2zeta_a*2*pc_ol(3,11))

pc(20,11) = PA(3)*pc(7,11) + inv_2zeta*(2*pc(4,11) + pc(7,8)) + 2*xi*(pc_ol(20&
&,11) - inv_2zeta_a*2*pc_ol(4,11))

pc(1,12) = PB(2)*pc(1,5) + 2*xi*(pc_ol(1,12))

pc(2,12) = PA(1)*pc(1,12) + inv_2zeta*(2*pc(1,8)) + 2*xi*(pc_ol(2,12))

pc(3,12) = PA(2)*pc(1,12) + inv_2zeta*(pc(1,5)) + 2*xi*(pc_ol(3,12))

pc(4,12) = PA(3)*pc(1,12) + 2*xi*(pc_ol(4,12))

pc(5,12) = PA(1)*pc(2,12) + inv_2zeta*(pc(1,12) + 2*pc(2,8)) + 2*xi*(pc_ol(5,1&
&2) - inv_2zeta_a*pc_ol(1,12))

pc(6,12) = PA(2)*pc(3,12) + inv_2zeta*(pc(1,12) + pc(3,5)) + 2*xi*(pc_ol(6,12)&
& - inv_2zeta_a*pc_ol(1,12))

pc(7,12) = PA(3)*pc(4,12) + inv_2zeta*(pc(1,12)) + 2*xi*(pc_ol(7,12) - inv_2ze&
&ta_a*pc_ol(1,12))

pc(8,12) = PA(2)*pc(2,12) + inv_2zeta*(pc(2,5)) + 2*xi*(pc_ol(8,12))

pc(9,12) = PA(3)*pc(2,12) + 2*xi*(pc_ol(9,12))

pc(10,12) = PA(3)*pc(3,12) + 2*xi*(pc_ol(10,12))

pc(11,12) = PA(3)*pc(8,12) + 2*xi*(pc_ol(11,12))

pc(12,12) = PA(2)*pc(5,12) + inv_2zeta*(pc(5,5)) + 2*xi*(pc_ol(12,12))

pc(13,12) = PA(3)*pc(5,12) + 2*xi*(pc_ol(13,12))

pc(14,12) = PA(1)*pc(6,12) + inv_2zeta*(2*pc(6,8)) + 2*xi*(pc_ol(14,12))

pc(15,12) = PA(3)*pc(6,12) + 2*xi*(pc_ol(15,12))

pc(16,12) = PA(1)*pc(7,12) + inv_2zeta*(2*pc(7,8)) + 2*xi*(pc_ol(16,12))

pc(17,12) = PA(2)*pc(7,12) + inv_2zeta*(pc(7,5)) + 2*xi*(pc_ol(17,12))

pc(18,12) = PA(1)*pc(5,12) + inv_2zeta*(2*pc(2,12) + 2*pc(5,8)) + 2*xi*(pc_ol(&
&18,12) - inv_2zeta_a*2*pc_ol(2,12))

pc(19,12) = PA(2)*pc(6,12) + inv_2zeta*(2*pc(3,12) + pc(6,5)) + 2*xi*(pc_ol(19&
&,12) - inv_2zeta_a*2*pc_ol(3,12))

pc(20,12) = PA(3)*pc(7,12) + inv_2zeta*(2*pc(4,12)) + 2*xi*(pc_ol(20,12) - inv&
&_2zeta_a*2*pc_ol(4,12))

pc(1,13) = PB(3)*pc(1,5) + 2*xi*(pc_ol(1,13))

pc(2,13) = PA(1)*pc(1,13) + inv_2zeta*(2*pc(1,9)) + 2*xi*(pc_ol(2,13))

pc(3,13) = PA(2)*pc(1,13) + 2*xi*(pc_ol(3,13))

pc(4,13) = PA(3)*pc(1,13) + inv_2zeta*(pc(1,5)) + 2*xi*(pc_ol(4,13))

pc(5,13) = PA(1)*pc(2,13) + inv_2zeta*(pc(1,13) + 2*pc(2,9)) + 2*xi*(pc_ol(5,1&
&3) - inv_2zeta_a*pc_ol(1,13))

pc(6,13) = PA(2)*pc(3,13) + inv_2zeta*(pc(1,13)) + 2*xi*(pc_ol(6,13) - inv_2ze&
&ta_a*pc_ol(1,13))

pc(7,13) = PA(3)*pc(4,13) + inv_2zeta*(pc(1,13) + pc(4,5)) + 2*xi*(pc_ol(7,13)&
& - inv_2zeta_a*pc_ol(1,13))

pc(8,13) = PA(2)*pc(2,13) + 2*xi*(pc_ol(8,13))

pc(9,13) = PA(3)*pc(2,13) + inv_2zeta*(pc(2,5)) + 2*xi*(pc_ol(9,13))

pc(10,13) = PA(3)*pc(3,13) + inv_2zeta*(pc(3,5)) + 2*xi*(pc_ol(10,13))

pc(11,13) = PA(3)*pc(8,13) + inv_2zeta*(pc(8,5)) + 2*xi*(pc_ol(11,13))

pc(12,13) = PA(2)*pc(5,13) + 2*xi*(pc_ol(12,13))

pc(13,13) = PA(3)*pc(5,13) + inv_2zeta*(pc(5,5)) + 2*xi*(pc_ol(13,13))

pc(14,13) = PA(1)*pc(6,13) + inv_2zeta*(2*pc(6,9)) + 2*xi*(pc_ol(14,13))

pc(15,13) = PA(3)*pc(6,13) + inv_2zeta*(pc(6,5)) + 2*xi*(pc_ol(15,13))

pc(16,13) = PA(1)*pc(7,13) + inv_2zeta*(2*pc(7,9)) + 2*xi*(pc_ol(16,13))

pc(17,13) = PA(2)*pc(7,13) + 2*xi*(pc_ol(17,13))

pc(18,13) = PA(1)*pc(5,13) + inv_2zeta*(2*pc(2,13) + 2*pc(5,9)) + 2*xi*(pc_ol(&
&18,13) - inv_2zeta_a*2*pc_ol(2,13))

pc(19,13) = PA(2)*pc(6,13) + inv_2zeta*(2*pc(3,13)) + 2*xi*(pc_ol(19,13) - inv&
&_2zeta_a*2*pc_ol(3,13))

pc(20,13) = PA(3)*pc(7,13) + inv_2zeta*(2*pc(4,13) + pc(7,5)) + 2*xi*(pc_ol(20&
&,13) - inv_2zeta_a*2*pc_ol(4,13))

pc(1,14) = PB(1)*pc(1,6) + 2*xi*(pc_ol(1,14))

pc(2,14) = PA(1)*pc(1,14) + inv_2zeta*(pc(1,6)) + 2*xi*(pc_ol(2,14))

pc(3,14) = PA(2)*pc(1,14) + inv_2zeta*(2*pc(1,8)) + 2*xi*(pc_ol(3,14))

pc(4,14) = PA(3)*pc(1,14) + 2*xi*(pc_ol(4,14))

pc(5,14) = PA(1)*pc(2,14) + inv_2zeta*(pc(1,14) + pc(2,6)) + 2*xi*(pc_ol(5,14)&
& - inv_2zeta_a*pc_ol(1,14))

pc(6,14) = PA(2)*pc(3,14) + inv_2zeta*(pc(1,14) + 2*pc(3,8)) + 2*xi*(pc_ol(6,1&
&4) - inv_2zeta_a*pc_ol(1,14))

pc(7,14) = PA(3)*pc(4,14) + inv_2zeta*(pc(1,14)) + 2*xi*(pc_ol(7,14) - inv_2ze&
&ta_a*pc_ol(1,14))

pc(8,14) = PA(2)*pc(2,14) + inv_2zeta*(2*pc(2,8)) + 2*xi*(pc_ol(8,14))

pc(9,14) = PA(3)*pc(2,14) + 2*xi*(pc_ol(9,14))

pc(10,14) = PA(3)*pc(3,14) + 2*xi*(pc_ol(10,14))

pc(11,14) = PA(3)*pc(8,14) + 2*xi*(pc_ol(11,14))

pc(12,14) = PA(2)*pc(5,14) + inv_2zeta*(2*pc(5,8)) + 2*xi*(pc_ol(12,14))

pc(13,14) = PA(3)*pc(5,14) + 2*xi*(pc_ol(13,14))

pc(14,14) = PA(1)*pc(6,14) + inv_2zeta*(pc(6,6)) + 2*xi*(pc_ol(14,14))

pc(15,14) = PA(3)*pc(6,14) + 2*xi*(pc_ol(15,14))

pc(16,14) = PA(1)*pc(7,14) + inv_2zeta*(pc(7,6)) + 2*xi*(pc_ol(16,14))

pc(17,14) = PA(2)*pc(7,14) + inv_2zeta*(2*pc(7,8)) + 2*xi*(pc_ol(17,14))

pc(18,14) = PA(1)*pc(5,14) + inv_2zeta*(2*pc(2,14) + pc(5,6)) + 2*xi*(pc_ol(18&
&,14) - inv_2zeta_a*2*pc_ol(2,14))

pc(19,14) = PA(2)*pc(6,14) + inv_2zeta*(2*pc(3,14) + 2*pc(6,8)) + 2*xi*(pc_ol(&
&19,14) - inv_2zeta_a*2*pc_ol(3,14))

pc(20,14) = PA(3)*pc(7,14) + inv_2zeta*(2*pc(4,14)) + 2*xi*(pc_ol(20,14) - inv&
&_2zeta_a*2*pc_ol(4,14))

pc(1,15) = PB(3)*pc(1,6) + 2*xi*(pc_ol(1,15))

pc(2,15) = PA(1)*pc(1,15) + 2*xi*(pc_ol(2,15))

pc(3,15) = PA(2)*pc(1,15) + inv_2zeta*(2*pc(1,10)) + 2*xi*(pc_ol(3,15))

pc(4,15) = PA(3)*pc(1,15) + inv_2zeta*(pc(1,6)) + 2*xi*(pc_ol(4,15))

pc(5,15) = PA(1)*pc(2,15) + inv_2zeta*(pc(1,15)) + 2*xi*(pc_ol(5,15) - inv_2ze&
&ta_a*pc_ol(1,15))

pc(6,15) = PA(2)*pc(3,15) + inv_2zeta*(pc(1,15) + 2*pc(3,10)) + 2*xi*(pc_ol(6,&
&15) - inv_2zeta_a*pc_ol(1,15))

pc(7,15) = PA(3)*pc(4,15) + inv_2zeta*(pc(1,15) + pc(4,6)) + 2*xi*(pc_ol(7,15)&
& - inv_2zeta_a*pc_ol(1,15))

pc(8,15) = PA(2)*pc(2,15) + inv_2zeta*(2*pc(2,10)) + 2*xi*(pc_ol(8,15))

pc(9,15) = PA(3)*pc(2,15) + inv_2zeta*(pc(2,6)) + 2*xi*(pc_ol(9,15))

pc(10,15) = PA(3)*pc(3,15) + inv_2zeta*(pc(3,6)) + 2*xi*(pc_ol(10,15))

pc(11,15) = PA(3)*pc(8,15) + inv_2zeta*(pc(8,6)) + 2*xi*(pc_ol(11,15))

pc(12,15) = PA(2)*pc(5,15) + inv_2zeta*(2*pc(5,10)) + 2*xi*(pc_ol(12,15))

pc(13,15) = PA(3)*pc(5,15) + inv_2zeta*(pc(5,6)) + 2*xi*(pc_ol(13,15))

pc(14,15) = PA(1)*pc(6,15) + 2*xi*(pc_ol(14,15))

pc(15,15) = PA(3)*pc(6,15) + inv_2zeta*(pc(6,6)) + 2*xi*(pc_ol(15,15))

pc(16,15) = PA(1)*pc(7,15) + 2*xi*(pc_ol(16,15))

pc(17,15) = PA(2)*pc(7,15) + inv_2zeta*(2*pc(7,10)) + 2*xi*(pc_ol(17,15))

pc(18,15) = PA(1)*pc(5,15) + inv_2zeta*(2*pc(2,15)) + 2*xi*(pc_ol(18,15) - inv&
&_2zeta_a*2*pc_ol(2,15))

pc(19,15) = PA(2)*pc(6,15) + inv_2zeta*(2*pc(3,15) + 2*pc(6,10)) + 2*xi*(pc_ol&
&(19,15) - inv_2zeta_a*2*pc_ol(3,15))

pc(20,15) = PA(3)*pc(7,15) + inv_2zeta*(2*pc(4,15) + pc(7,6)) + 2*xi*(pc_ol(20&
&,15) - inv_2zeta_a*2*pc_ol(4,15))

pc(1,16) = PB(1)*pc(1,7) + 2*xi*(pc_ol(1,16))

pc(2,16) = PA(1)*pc(1,16) + inv_2zeta*(pc(1,7)) + 2*xi*(pc_ol(2,16))

pc(3,16) = PA(2)*pc(1,16) + 2*xi*(pc_ol(3,16))

pc(4,16) = PA(3)*pc(1,16) + inv_2zeta*(2*pc(1,9)) + 2*xi*(pc_ol(4,16))

pc(5,16) = PA(1)*pc(2,16) + inv_2zeta*(pc(1,16) + pc(2,7)) + 2*xi*(pc_ol(5,16)&
& - inv_2zeta_a*pc_ol(1,16))

pc(6,16) = PA(2)*pc(3,16) + inv_2zeta*(pc(1,16)) + 2*xi*(pc_ol(6,16) - inv_2ze&
&ta_a*pc_ol(1,16))

pc(7,16) = PA(3)*pc(4,16) + inv_2zeta*(pc(1,16) + 2*pc(4,9)) + 2*xi*(pc_ol(7,1&
&6) - inv_2zeta_a*pc_ol(1,16))

pc(8,16) = PA(2)*pc(2,16) + 2*xi*(pc_ol(8,16))

pc(9,16) = PA(3)*pc(2,16) + inv_2zeta*(2*pc(2,9)) + 2*xi*(pc_ol(9,16))

pc(10,16) = PA(3)*pc(3,16) + inv_2zeta*(2*pc(3,9)) + 2*xi*(pc_ol(10,16))

pc(11,16) = PA(3)*pc(8,16) + inv_2zeta*(2*pc(8,9)) + 2*xi*(pc_ol(11,16))

pc(12,16) = PA(2)*pc(5,16) + 2*xi*(pc_ol(12,16))

pc(13,16) = PA(3)*pc(5,16) + inv_2zeta*(2*pc(5,9)) + 2*xi*(pc_ol(13,16))

pc(14,16) = PA(1)*pc(6,16) + inv_2zeta*(pc(6,7)) + 2*xi*(pc_ol(14,16))

pc(15,16) = PA(3)*pc(6,16) + inv_2zeta*(2*pc(6,9)) + 2*xi*(pc_ol(15,16))

pc(16,16) = PA(1)*pc(7,16) + inv_2zeta*(pc(7,7)) + 2*xi*(pc_ol(16,16))

pc(17,16) = PA(2)*pc(7,16) + 2*xi*(pc_ol(17,16))

pc(18,16) = PA(1)*pc(5,16) + inv_2zeta*(2*pc(2,16) + pc(5,7)) + 2*xi*(pc_ol(18&
&,16) - inv_2zeta_a*2*pc_ol(2,16))

pc(19,16) = PA(2)*pc(6,16) + inv_2zeta*(2*pc(3,16)) + 2*xi*(pc_ol(19,16) - inv&
&_2zeta_a*2*pc_ol(3,16))

pc(20,16) = PA(3)*pc(7,16) + inv_2zeta*(2*pc(4,16) + 2*pc(7,9)) + 2*xi*(pc_ol(&
&20,16) - inv_2zeta_a*2*pc_ol(4,16))

pc(1,17) = PB(2)*pc(1,7) + 2*xi*(pc_ol(1,17))

pc(2,17) = PA(1)*pc(1,17) + 2*xi*(pc_ol(2,17))

pc(3,17) = PA(2)*pc(1,17) + inv_2zeta*(pc(1,7)) + 2*xi*(pc_ol(3,17))

pc(4,17) = PA(3)*pc(1,17) + inv_2zeta*(2*pc(1,10)) + 2*xi*(pc_ol(4,17))

pc(5,17) = PA(1)*pc(2,17) + inv_2zeta*(pc(1,17)) + 2*xi*(pc_ol(5,17) - inv_2ze&
&ta_a*pc_ol(1,17))

pc(6,17) = PA(2)*pc(3,17) + inv_2zeta*(pc(1,17) + pc(3,7)) + 2*xi*(pc_ol(6,17)&
& - inv_2zeta_a*pc_ol(1,17))

pc(7,17) = PA(3)*pc(4,17) + inv_2zeta*(pc(1,17) + 2*pc(4,10)) + 2*xi*(pc_ol(7,&
&17) - inv_2zeta_a*pc_ol(1,17))

pc(8,17) = PA(2)*pc(2,17) + inv_2zeta*(pc(2,7)) + 2*xi*(pc_ol(8,17))

pc(9,17) = PA(3)*pc(2,17) + inv_2zeta*(2*pc(2,10)) + 2*xi*(pc_ol(9,17))

pc(10,17) = PA(3)*pc(3,17) + inv_2zeta*(2*pc(3,10)) + 2*xi*(pc_ol(10,17))

pc(11,17) = PA(3)*pc(8,17) + inv_2zeta*(2*pc(8,10)) + 2*xi*(pc_ol(11,17))

pc(12,17) = PA(2)*pc(5,17) + inv_2zeta*(pc(5,7)) + 2*xi*(pc_ol(12,17))

pc(13,17) = PA(3)*pc(5,17) + inv_2zeta*(2*pc(5,10)) + 2*xi*(pc_ol(13,17))

pc(14,17) = PA(1)*pc(6,17) + 2*xi*(pc_ol(14,17))

pc(15,17) = PA(3)*pc(6,17) + inv_2zeta*(2*pc(6,10)) + 2*xi*(pc_ol(15,17))

pc(16,17) = PA(1)*pc(7,17) + 2*xi*(pc_ol(16,17))

pc(17,17) = PA(2)*pc(7,17) + inv_2zeta*(pc(7,7)) + 2*xi*(pc_ol(17,17))

pc(18,17) = PA(1)*pc(5,17) + inv_2zeta*(2*pc(2,17)) + 2*xi*(pc_ol(18,17) - inv&
&_2zeta_a*2*pc_ol(2,17))

pc(19,17) = PA(2)*pc(6,17) + inv_2zeta*(2*pc(3,17) + pc(6,7)) + 2*xi*(pc_ol(19&
&,17) - inv_2zeta_a*2*pc_ol(3,17))

pc(20,17) = PA(3)*pc(7,17) + inv_2zeta*(2*pc(4,17) + 2*pc(7,10)) + 2*xi*(pc_ol&
&(20,17) - inv_2zeta_a*2*pc_ol(4,17))

pc(1,18) = PB(1)*pc(1,5) + inv_2zeta*(2*pc(1,2)) + 2*xi*(pc_ol(1,18) - inv_2ze&
&ta_b*2*pc_ol(1,2))

pc(2,18) = PA(1)*pc(1,18) + inv_2zeta*(3*pc(1,5)) + 2*xi*(pc_ol(2,18))

pc(3,18) = PA(2)*pc(1,18) + 2*xi*(pc_ol(3,18))

pc(4,18) = PA(3)*pc(1,18) + 2*xi*(pc_ol(4,18))

pc(5,18) = PA(1)*pc(2,18) + inv_2zeta*(pc(1,18) + 3*pc(2,5)) + 2*xi*(pc_ol(5,1&
&8) - inv_2zeta_a*pc_ol(1,18))

pc(6,18) = PA(2)*pc(3,18) + inv_2zeta*(pc(1,18)) + 2*xi*(pc_ol(6,18) - inv_2ze&
&ta_a*pc_ol(1,18))

pc(7,18) = PA(3)*pc(4,18) + inv_2zeta*(pc(1,18)) + 2*xi*(pc_ol(7,18) - inv_2ze&
&ta_a*pc_ol(1,18))

pc(8,18) = PA(2)*pc(2,18) + 2*xi*(pc_ol(8,18))

pc(9,18) = PA(3)*pc(2,18) + 2*xi*(pc_ol(9,18))

pc(10,18) = PA(3)*pc(3,18) + 2*xi*(pc_ol(10,18))

pc(11,18) = PA(3)*pc(8,18) + 2*xi*(pc_ol(11,18))

pc(12,18) = PA(2)*pc(5,18) + 2*xi*(pc_ol(12,18))

pc(13,18) = PA(3)*pc(5,18) + 2*xi*(pc_ol(13,18))

pc(14,18) = PA(1)*pc(6,18) + inv_2zeta*(3*pc(6,5)) + 2*xi*(pc_ol(14,18))

pc(15,18) = PA(3)*pc(6,18) + 2*xi*(pc_ol(15,18))

pc(16,18) = PA(1)*pc(7,18) + inv_2zeta*(3*pc(7,5)) + 2*xi*(pc_ol(16,18))

pc(17,18) = PA(2)*pc(7,18) + 2*xi*(pc_ol(17,18))

pc(18,18) = PA(1)*pc(5,18) + inv_2zeta*(2*pc(2,18) + 3*pc(5,5)) + 2*xi*(pc_ol(&
&18,18) - inv_2zeta_a*2*pc_ol(2,18))

pc(19,18) = PA(2)*pc(6,18) + inv_2zeta*(2*pc(3,18)) + 2*xi*(pc_ol(19,18) - inv&
&_2zeta_a*2*pc_ol(3,18))

pc(20,18) = PA(3)*pc(7,18) + inv_2zeta*(2*pc(4,18)) + 2*xi*(pc_ol(20,18) - inv&
&_2zeta_a*2*pc_ol(4,18))

pc(1,19) = PB(2)*pc(1,6) + inv_2zeta*(2*pc(1,3)) + 2*xi*(pc_ol(1,19) - inv_2ze&
&ta_b*2*pc_ol(1,3))

pc(2,19) = PA(1)*pc(1,19) + 2*xi*(pc_ol(2,19))

pc(3,19) = PA(2)*pc(1,19) + inv_2zeta*(3*pc(1,6)) + 2*xi*(pc_ol(3,19))

pc(4,19) = PA(3)*pc(1,19) + 2*xi*(pc_ol(4,19))

pc(5,19) = PA(1)*pc(2,19) + inv_2zeta*(pc(1,19)) + 2*xi*(pc_ol(5,19) - inv_2ze&
&ta_a*pc_ol(1,19))

pc(6,19) = PA(2)*pc(3,19) + inv_2zeta*(pc(1,19) + 3*pc(3,6)) + 2*xi*(pc_ol(6,1&
&9) - inv_2zeta_a*pc_ol(1,19))

pc(7,19) = PA(3)*pc(4,19) + inv_2zeta*(pc(1,19)) + 2*xi*(pc_ol(7,19) - inv_2ze&
&ta_a*pc_ol(1,19))

pc(8,19) = PA(2)*pc(2,19) + inv_2zeta*(3*pc(2,6)) + 2*xi*(pc_ol(8,19))

pc(9,19) = PA(3)*pc(2,19) + 2*xi*(pc_ol(9,19))

pc(10,19) = PA(3)*pc(3,19) + 2*xi*(pc_ol(10,19))

pc(11,19) = PA(3)*pc(8,19) + 2*xi*(pc_ol(11,19))

pc(12,19) = PA(2)*pc(5,19) + inv_2zeta*(3*pc(5,6)) + 2*xi*(pc_ol(12,19))

pc(13,19) = PA(3)*pc(5,19) + 2*xi*(pc_ol(13,19))

pc(14,19) = PA(1)*pc(6,19) + 2*xi*(pc_ol(14,19))

pc(15,19) = PA(3)*pc(6,19) + 2*xi*(pc_ol(15,19))

pc(16,19) = PA(1)*pc(7,19) + 2*xi*(pc_ol(16,19))

pc(17,19) = PA(2)*pc(7,19) + inv_2zeta*(3*pc(7,6)) + 2*xi*(pc_ol(17,19))

pc(18,19) = PA(1)*pc(5,19) + inv_2zeta*(2*pc(2,19)) + 2*xi*(pc_ol(18,19) - inv&
&_2zeta_a*2*pc_ol(2,19))

pc(19,19) = PA(2)*pc(6,19) + inv_2zeta*(2*pc(3,19) + 3*pc(6,6)) + 2*xi*(pc_ol(&
&19,19) - inv_2zeta_a*2*pc_ol(3,19))

pc(20,19) = PA(3)*pc(7,19) + inv_2zeta*(2*pc(4,19)) + 2*xi*(pc_ol(20,19) - inv&
&_2zeta_a*2*pc_ol(4,19))

pc(1,20) = PB(3)*pc(1,7) + inv_2zeta*(2*pc(1,4)) + 2*xi*(pc_ol(1,20) - inv_2ze&
&ta_b*2*pc_ol(1,4))

pc(2,20) = PA(1)*pc(1,20) + 2*xi*(pc_ol(2,20))

pc(3,20) = PA(2)*pc(1,20) + 2*xi*(pc_ol(3,20))

pc(4,20) = PA(3)*pc(1,20) + inv_2zeta*(3*pc(1,7)) + 2*xi*(pc_ol(4,20))

pc(5,20) = PA(1)*pc(2,20) + inv_2zeta*(pc(1,20)) + 2*xi*(pc_ol(5,20) - inv_2ze&
&ta_a*pc_ol(1,20))

pc(6,20) = PA(2)*pc(3,20) + inv_2zeta*(pc(1,20)) + 2*xi*(pc_ol(6,20) - inv_2ze&
&ta_a*pc_ol(1,20))

pc(7,20) = PA(3)*pc(4,20) + inv_2zeta*(pc(1,20) + 3*pc(4,7)) + 2*xi*(pc_ol(7,2&
&0) - inv_2zeta_a*pc_ol(1,20))

pc(8,20) = PA(2)*pc(2,20) + 2*xi*(pc_ol(8,20))

pc(9,20) = PA(3)*pc(2,20) + inv_2zeta*(3*pc(2,7)) + 2*xi*(pc_ol(9,20))

pc(10,20) = PA(3)*pc(3,20) + inv_2zeta*(3*pc(3,7)) + 2*xi*(pc_ol(10,20))

pc(11,20) = PA(3)*pc(8,20) + inv_2zeta*(3*pc(8,7)) + 2*xi*(pc_ol(11,20))

pc(12,20) = PA(2)*pc(5,20) + 2*xi*(pc_ol(12,20))

pc(13,20) = PA(3)*pc(5,20) + inv_2zeta*(3*pc(5,7)) + 2*xi*(pc_ol(13,20))

pc(14,20) = PA(1)*pc(6,20) + 2*xi*(pc_ol(14,20))

pc(15,20) = PA(3)*pc(6,20) + inv_2zeta*(3*pc(6,7)) + 2*xi*(pc_ol(15,20))

pc(16,20) = PA(1)*pc(7,20) + 2*xi*(pc_ol(16,20))

pc(17,20) = PA(2)*pc(7,20) + 2*xi*(pc_ol(17,20))

pc(18,20) = PA(1)*pc(5,20) + inv_2zeta*(2*pc(2,20)) + 2*xi*(pc_ol(18,20) - inv&
&_2zeta_a*2*pc_ol(2,20))

pc(19,20) = PA(2)*pc(6,20) + inv_2zeta*(2*pc(3,20)) + 2*xi*(pc_ol(19,20) - inv&
&_2zeta_a*2*pc_ol(3,20))

pc(20,20) = PA(3)*pc(7,20) + inv_2zeta*(2*pc(4,20) + 3*pc(7,7)) + 2*xi*(pc_ol(&
&20,20) - inv_2zeta_a*2*pc_ol(4,20))

sh(1,1) = pc(1,1)

sh(2,1) = pc(2,1)

sh(3,1) = pc(3,1)

sh(4,1) = pc(4,1)

sh(5,1) = pc(8,1)

sh(6,1) = pc(9,1)

sh(7,1) = pc(10,1)

sh(8,1) = pc(5,1) - pc(6,1)

sh(9,1) = 2*pc(7,1) - pc(5,1) - pc(6,1)

sh(10,1) = pc(11,1)

sh(11,1) = pc(13,1) - pc(15,1)

sh(12,1) = pc(18,1) - 3*pc(14,1)

sh(13,1) = 3*pc(12,1) - pc(19,1)

sh(14,1) = 2*pc(20,1) - 3*pc(13,1) - 3*pc(15,1)

sh(15,1) = 4*pc(16,1) - pc(18,1) - pc(14,1)

sh(16,1) = 4*pc(17,1) - pc(12,1) - pc(19,1)

sh(1,2) = pc(1,2)

sh(2,2) = pc(2,2)

sh(3,2) = pc(3,2)

sh(4,2) = pc(4,2)

sh(5,2) = pc(8,2)

sh(6,2) = pc(9,2)

sh(7,2) = pc(10,2)

sh(8,2) = pc(5,2) - pc(6,2)

sh(9,2) = 2*pc(7,2) - pc(5,2) - pc(6,2)

sh(10,2) = pc(11,2)

sh(11,2) = pc(13,2) - pc(15,2)

sh(12,2) = pc(18,2) - 3*pc(14,2)

sh(13,2) = 3*pc(12,2) - pc(19,2)

sh(14,2) = 2*pc(20,2) - 3*pc(13,2) - 3*pc(15,2)

sh(15,2) = 4*pc(16,2) - pc(18,2) - pc(14,2)

sh(16,2) = 4*pc(17,2) - pc(12,2) - pc(19,2)

sh(1,3) = pc(1,3)

sh(2,3) = pc(2,3)

sh(3,3) = pc(3,3)

sh(4,3) = pc(4,3)

sh(5,3) = pc(8,3)

sh(6,3) = pc(9,3)

sh(7,3) = pc(10,3)

sh(8,3) = pc(5,3) - pc(6,3)

sh(9,3) = 2*pc(7,3) - pc(5,3) - pc(6,3)

sh(10,3) = pc(11,3)

sh(11,3) = pc(13,3) - pc(15,3)

sh(12,3) = pc(18,3) - 3*pc(14,3)

sh(13,3) = 3*pc(12,3) - pc(19,3)

sh(14,3) = 2*pc(20,3) - 3*pc(13,3) - 3*pc(15,3)

sh(15,3) = 4*pc(16,3) - pc(18,3) - pc(14,3)

sh(16,3) = 4*pc(17,3) - pc(12,3) - pc(19,3)

sh(1,4) = pc(1,4)

sh(2,4) = pc(2,4)

sh(3,4) = pc(3,4)

sh(4,4) = pc(4,4)

sh(5,4) = pc(8,4)

sh(6,4) = pc(9,4)

sh(7,4) = pc(10,4)

sh(8,4) = pc(5,4) - pc(6,4)

sh(9,4) = 2*pc(7,4) - pc(5,4) - pc(6,4)

sh(10,4) = pc(11,4)

sh(11,4) = pc(13,4) - pc(15,4)

sh(12,4) = pc(18,4) - 3*pc(14,4)

sh(13,4) = 3*pc(12,4) - pc(19,4)

sh(14,4) = 2*pc(20,4) - 3*pc(13,4) - 3*pc(15,4)

sh(15,4) = 4*pc(16,4) - pc(18,4) - pc(14,4)

sh(16,4) = 4*pc(17,4) - pc(12,4) - pc(19,4)

sh(1,5) = pc(1,8)

sh(2,5) = pc(2,8)

sh(3,5) = pc(3,8)

sh(4,5) = pc(4,8)

sh(5,5) = pc(8,8)

sh(6,5) = pc(9,8)

sh(7,5) = pc(10,8)

sh(8,5) = pc(5,8) - pc(6,8)

sh(9,5) = 2*pc(7,8) - pc(5,8) - pc(6,8)

sh(10,5) = pc(11,8)

sh(11,5) = pc(13,8) - pc(15,8)

sh(12,5) = pc(18,8) - 3*pc(14,8)

sh(13,5) = 3*pc(12,8) - pc(19,8)

sh(14,5) = 2*pc(20,8) - 3*pc(13,8) - 3*pc(15,8)

sh(15,5) = 4*pc(16,8) - pc(18,8) - pc(14,8)

sh(16,5) = 4*pc(17,8) - pc(12,8) - pc(19,8)

sh(1,6) = pc(1,9)

sh(2,6) = pc(2,9)

sh(3,6) = pc(3,9)

sh(4,6) = pc(4,9)

sh(5,6) = pc(8,9)

sh(6,6) = pc(9,9)

sh(7,6) = pc(10,9)

sh(8,6) = pc(5,9) - pc(6,9)

sh(9,6) = 2*pc(7,9) - pc(5,9) - pc(6,9)

sh(10,6) = pc(11,9)

sh(11,6) = pc(13,9) - pc(15,9)

sh(12,6) = pc(18,9) - 3*pc(14,9)

sh(13,6) = 3*pc(12,9) - pc(19,9)

sh(14,6) = 2*pc(20,9) - 3*pc(13,9) - 3*pc(15,9)

sh(15,6) = 4*pc(16,9) - pc(18,9) - pc(14,9)

sh(16,6) = 4*pc(17,9) - pc(12,9) - pc(19,9)

sh(1,7) = pc(1,10)

sh(2,7) = pc(2,10)

sh(3,7) = pc(3,10)

sh(4,7) = pc(4,10)

sh(5,7) = pc(8,10)

sh(6,7) = pc(9,10)

sh(7,7) = pc(10,10)

sh(8,7) = pc(5,10) - pc(6,10)

sh(9,7) = 2*pc(7,10) - pc(5,10) - pc(6,10)

sh(10,7) = pc(11,10)

sh(11,7) = pc(13,10) - pc(15,10)

sh(12,7) = pc(18,10) - 3*pc(14,10)

sh(13,7) = 3*pc(12,10) - pc(19,10)

sh(14,7) = 2*pc(20,10) - 3*pc(13,10) - 3*pc(15,10)

sh(15,7) = 4*pc(16,10) - pc(18,10) - pc(14,10)

sh(16,7) = 4*pc(17,10) - pc(12,10) - pc(19,10)

sh(1,8) = pc(1,5) - pc(1,6)

sh(2,8) = pc(2,5) - pc(2,6)

sh(3,8) = pc(3,5) - pc(3,6)

sh(4,8) = pc(4,5) - pc(4,6)

sh(5,8) = pc(8,5) - pc(8,6)

sh(6,8) = pc(9,5) - pc(9,6)

sh(7,8) = pc(10,5) - pc(10,6)

sh(8,8) = pc(5,5) - pc(5,6) - pc(6,5) + pc(6,6)

sh(9,8) = 2*pc(7,5) - 2*pc(7,6) - pc(5,5) + pc(5,6) - pc(6,5) + pc(6,6)

sh(10,8) = pc(11,5) - pc(11,6)

sh(11,8) = pc(13,5) - pc(13,6) - pc(15,5) + pc(15,6)

sh(12,8) = pc(18,5) - pc(18,6) - 3*pc(14,5) + 3*pc(14,6)

sh(13,8) = 3*pc(12,5) - 3*pc(12,6) - pc(19,5) + pc(19,6)

sh(14,8) = 2*pc(20,5) - 2*pc(20,6) - 3*pc(13,5) + 3*pc(13,6) - 3*pc(15,5) + 3*&
&pc(15,6)

sh(15,8) = 4*pc(16,5) - 4*pc(16,6) - pc(18,5) + pc(18,6) - pc(14,5) + pc(14,6)

sh(16,8) = 4*pc(17,5) - 4*pc(17,6) - pc(12,5) + pc(12,6) - pc(19,5) + pc(19,6)

sh(1,9) = 2*pc(1,7) - pc(1,5) - pc(1,6)

sh(2,9) = 2*pc(2,7) - pc(2,5) - pc(2,6)

sh(3,9) = 2*pc(3,7) - pc(3,5) - pc(3,6)

sh(4,9) = 2*pc(4,7) - pc(4,5) - pc(4,6)

sh(5,9) = 2*pc(8,7) - pc(8,5) - pc(8,6)

sh(6,9) = 2*pc(9,7) - pc(9,5) - pc(9,6)

sh(7,9) = 2*pc(10,7) - pc(10,5) - pc(10,6)

sh(8,9) = 2*pc(5,7) - pc(5,5) - pc(5,6) - 2*pc(6,7) + pc(6,5) + pc(6,6)

sh(9,9) = 4*pc(7,7) - 2*pc(7,5) - 2*pc(7,6) - 2*pc(5,7) + pc(5,5) + pc(5,6) - &
&2*pc(6,7) + pc(6,5) + pc(6,6)

sh(10,9) = 2*pc(11,7) - pc(11,5) - pc(11,6)

sh(11,9) = 2*pc(13,7) - pc(13,5) - pc(13,6) - 2*pc(15,7) + pc(15,5) + pc(15,6)

sh(12,9) = 2*pc(18,7) - pc(18,5) - pc(18,6) - 6*pc(14,7) + 3*pc(14,5) + 3*pc(1&
&4,6)

sh(13,9) = 6*pc(12,7) - 3*pc(12,5) - 3*pc(12,6) - 2*pc(19,7) + pc(19,5) + pc(1&
&9,6)

sh(14,9) = 4*pc(20,7) - 2*pc(20,5) - 2*pc(20,6) - 6*pc(13,7) + 3*pc(13,5) + 3*&
&pc(13,6) - 6*pc(15,7) + 3*pc(15,5) + 3*pc(15,6)

sh(15,9) = 8*pc(16,7) - 4*pc(16,5) - 4*pc(16,6) - 2*pc(18,7) + pc(18,5) + pc(1&
&8,6) - 2*pc(14,7) + pc(14,5) + pc(14,6)

sh(16,9) = 8*pc(17,7) - 4*pc(17,5) - 4*pc(17,6) - 2*pc(12,7) + pc(12,5) + pc(1&
&2,6) - 2*pc(19,7) + pc(19,5) + pc(19,6)

sh(1,10) = pc(1,11)

sh(2,10) = pc(2,11)

sh(3,10) = pc(3,11)

sh(4,10) = pc(4,11)

sh(5,10) = pc(8,11)

sh(6,10) = pc(9,11)

sh(7,10) = pc(10,11)

sh(8,10) = pc(5,11) - pc(6,11)

sh(9,10) = 2*pc(7,11) - pc(5,11) - pc(6,11)

sh(10,10) = pc(11,11)

sh(11,10) = pc(13,11) - pc(15,11)

sh(12,10) = pc(18,11) - 3*pc(14,11)

sh(13,10) = 3*pc(12,11) - pc(19,11)

sh(14,10) = 2*pc(20,11) - 3*pc(13,11) - 3*pc(15,11)

sh(15,10) = 4*pc(16,11) - pc(18,11) - pc(14,11)

sh(16,10) = 4*pc(17,11) - pc(12,11) - pc(19,11)

sh(1,11) = pc(1,13) - pc(1,15)

sh(2,11) = pc(2,13) - pc(2,15)

sh(3,11) = pc(3,13) - pc(3,15)

sh(4,11) = pc(4,13) - pc(4,15)

sh(5,11) = pc(8,13) - pc(8,15)

sh(6,11) = pc(9,13) - pc(9,15)

sh(7,11) = pc(10,13) - pc(10,15)

sh(8,11) = pc(5,13) - pc(5,15) - pc(6,13) + pc(6,15)

sh(9,11) = 2*pc(7,13) - 2*pc(7,15) - pc(5,13) + pc(5,15) - pc(6,13) + pc(6,15)

sh(10,11) = pc(11,13) - pc(11,15)

sh(11,11) = pc(13,13) - pc(13,15) - pc(15,13) + pc(15,15)

sh(12,11) = pc(18,13) - pc(18,15) - 3*pc(14,13) + 3*pc(14,15)

sh(13,11) = 3*pc(12,13) - 3*pc(12,15) - pc(19,13) + pc(19,15)

sh(14,11) = 2*pc(20,13) - 2*pc(20,15) - 3*pc(13,13) + 3*pc(13,15) - 3*pc(15,13&
&) + 3*pc(15,15)

sh(15,11) = 4*pc(16,13) - 4*pc(16,15) - pc(18,13) + pc(18,15) - pc(14,13) + pc&
&(14,15)

sh(16,11) = 4*pc(17,13) - 4*pc(17,15) - pc(12,13) + pc(12,15) - pc(19,13) + pc&
&(19,15)

sh(1,12) = pc(1,18) - 3*pc(1,14)

sh(2,12) = pc(2,18) - 3*pc(2,14)

sh(3,12) = pc(3,18) - 3*pc(3,14)

sh(4,12) = pc(4,18) - 3*pc(4,14)

sh(5,12) = pc(8,18) - 3*pc(8,14)

sh(6,12) = pc(9,18) - 3*pc(9,14)

sh(7,12) = pc(10,18) - 3*pc(10,14)

sh(8,12) = pc(5,18) - 3*pc(5,14) - pc(6,18) + 3*pc(6,14)

sh(9,12) = 2*pc(7,18) - 6*pc(7,14) - pc(5,18) + 3*pc(5,14) - pc(6,18) + 3*pc(6&
&,14)

sh(10,12) = pc(11,18) - 3*pc(11,14)

sh(11,12) = pc(13,18) - 3*pc(13,14) - pc(15,18) + 3*pc(15,14)

sh(12,12) = pc(18,18) - 3*pc(18,14) - 3*pc(14,18) + 9*pc(14,14)

sh(13,12) = 3*pc(12,18) - 9*pc(12,14) - pc(19,18) + 3*pc(19,14)

sh(14,12) = 2*pc(20,18) - 6*pc(20,14) - 3*pc(13,18) + 9*pc(13,14) - 3*pc(15,18&
&) + 9*pc(15,14)

sh(15,12) = 4*pc(16,18) - 12*pc(16,14) - pc(18,18) + 3*pc(18,14) - pc(14,18) +&
& 3*pc(14,14)

sh(16,12) = 4*pc(17,18) - 12*pc(17,14) - pc(12,18) + 3*pc(12,14) - pc(19,18) +&
& 3*pc(19,14)

sh(1,13) = 3*pc(1,12) - pc(1,19)

sh(2,13) = 3*pc(2,12) - pc(2,19)

sh(3,13) = 3*pc(3,12) - pc(3,19)

sh(4,13) = 3*pc(4,12) - pc(4,19)

sh(5,13) = 3*pc(8,12) - pc(8,19)

sh(6,13) = 3*pc(9,12) - pc(9,19)

sh(7,13) = 3*pc(10,12) - pc(10,19)

sh(8,13) = 3*pc(5,12) - pc(5,19) - 3*pc(6,12) + pc(6,19)

sh(9,13) = 6*pc(7,12) - 2*pc(7,19) - 3*pc(5,12) + pc(5,19) - 3*pc(6,12) + pc(6&
&,19)

sh(10,13) = 3*pc(11,12) - pc(11,19)

sh(11,13) = 3*pc(13,12) - pc(13,19) - 3*pc(15,12) + pc(15,19)

sh(12,13) = 3*pc(18,12) - pc(18,19) - 9*pc(14,12) + 3*pc(14,19)

sh(13,13) = 9*pc(12,12) - 3*pc(12,19) - 3*pc(19,12) + pc(19,19)

sh(14,13) = 6*pc(20,12) - 2*pc(20,19) - 9*pc(13,12) + 3*pc(13,19) - 9*pc(15,12&
&) + 3*pc(15,19)

sh(15,13) = 12*pc(16,12) - 4*pc(16,19) - 3*pc(18,12) + pc(18,19) - 3*pc(14,12)&
& + pc(14,19)

sh(16,13) = 12*pc(17,12) - 4*pc(17,19) - 3*pc(12,12) + pc(12,19) - 3*pc(19,12)&
& + pc(19,19)

sh(1,14) = 2*pc(1,20) - 3*pc(1,13) - 3*pc(1,15)

sh(2,14) = 2*pc(2,20) - 3*pc(2,13) - 3*pc(2,15)

sh(3,14) = 2*pc(3,20) - 3*pc(3,13) - 3*pc(3,15)

sh(4,14) = 2*pc(4,20) - 3*pc(4,13) - 3*pc(4,15)

sh(5,14) = 2*pc(8,20) - 3*pc(8,13) - 3*pc(8,15)

sh(6,14) = 2*pc(9,20) - 3*pc(9,13) - 3*pc(9,15)

sh(7,14) = 2*pc(10,20) - 3*pc(10,13) - 3*pc(10,15)

sh(8,14) = 2*pc(5,20) - 3*pc(5,13) - 3*pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*&
&pc(6,15)

sh(9,14) = 4*pc(7,20) - 6*pc(7,13) - 6*pc(7,15) - 2*pc(5,20) + 3*pc(5,13) + 3*&
&pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*pc(6,15)

sh(10,14) = 2*pc(11,20) - 3*pc(11,13) - 3*pc(11,15)

sh(11,14) = 2*pc(13,20) - 3*pc(13,13) - 3*pc(13,15) - 2*pc(15,20) + 3*pc(15,13&
&) + 3*pc(15,15)

sh(12,14) = 2*pc(18,20) - 3*pc(18,13) - 3*pc(18,15) - 6*pc(14,20) + 9*pc(14,13&
&) + 9*pc(14,15)

sh(13,14) = 6*pc(12,20) - 9*pc(12,13) - 9*pc(12,15) - 2*pc(19,20) + 3*pc(19,13&
&) + 3*pc(19,15)

sh(14,14) = 4*pc(20,20) - 6*pc(20,13) - 6*pc(20,15) - 6*pc(13,20) + 9*pc(13,13&
&) + 9*pc(13,15) - 6*pc(15,20) + 9*pc(15,13) + 9*pc(15,15)

sh(15,14) = 8*pc(16,20) - 12*pc(16,13) - 12*pc(16,15) - 2*pc(18,20) + 3*pc(18,&
&13) + 3*pc(18,15) - 2*pc(14,20) + 3*pc(14,13) + 3*pc(14,15)

sh(16,14) = 8*pc(17,20) - 12*pc(17,13) - 12*pc(17,15) - 2*pc(12,20) + 3*pc(12,&
&13) + 3*pc(12,15) - 2*pc(19,20) + 3*pc(19,13) + 3*pc(19,15)

sh(1,15) = 4*pc(1,16) - pc(1,18) - pc(1,14)

sh(2,15) = 4*pc(2,16) - pc(2,18) - pc(2,14)

sh(3,15) = 4*pc(3,16) - pc(3,18) - pc(3,14)

sh(4,15) = 4*pc(4,16) - pc(4,18) - pc(4,14)

sh(5,15) = 4*pc(8,16) - pc(8,18) - pc(8,14)

sh(6,15) = 4*pc(9,16) - pc(9,18) - pc(9,14)

sh(7,15) = 4*pc(10,16) - pc(10,18) - pc(10,14)

sh(8,15) = 4*pc(5,16) - pc(5,18) - pc(5,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(9,15) = 8*pc(7,16) - 2*pc(7,18) - 2*pc(7,14) - 4*pc(5,16) + pc(5,18) + pc(5&
&,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(10,15) = 4*pc(11,16) - pc(11,18) - pc(11,14)

sh(11,15) = 4*pc(13,16) - pc(13,18) - pc(13,14) - 4*pc(15,16) + pc(15,18) + pc&
&(15,14)

sh(12,15) = 4*pc(18,16) - pc(18,18) - pc(18,14) - 12*pc(14,16) + 3*pc(14,18) +&
& 3*pc(14,14)

sh(13,15) = 12*pc(12,16) - 3*pc(12,18) - 3*pc(12,14) - 4*pc(19,16) + pc(19,18)&
& + pc(19,14)

sh(14,15) = 8*pc(20,16) - 2*pc(20,18) - 2*pc(20,14) - 12*pc(13,16) + 3*pc(13,1&
&8) + 3*pc(13,14) - 12*pc(15,16) + 3*pc(15,18) + 3*pc(15,14)

sh(15,15) = 16*pc(16,16) - 4*pc(16,18) - 4*pc(16,14) - 4*pc(18,16) + pc(18,18)&
& + pc(18,14) - 4*pc(14,16) + pc(14,18) + pc(14,14)

sh(16,15) = 16*pc(17,16) - 4*pc(17,18) - 4*pc(17,14) - 4*pc(12,16) + pc(12,18)&
& + pc(12,14) - 4*pc(19,16) + pc(19,18) + pc(19,14)

sh(1,16) = 4*pc(1,17) - pc(1,12) - pc(1,19)

sh(2,16) = 4*pc(2,17) - pc(2,12) - pc(2,19)

sh(3,16) = 4*pc(3,17) - pc(3,12) - pc(3,19)

sh(4,16) = 4*pc(4,17) - pc(4,12) - pc(4,19)

sh(5,16) = 4*pc(8,17) - pc(8,12) - pc(8,19)

sh(6,16) = 4*pc(9,17) - pc(9,12) - pc(9,19)

sh(7,16) = 4*pc(10,17) - pc(10,12) - pc(10,19)

sh(8,16) = 4*pc(5,17) - pc(5,12) - pc(5,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(9,16) = 8*pc(7,17) - 2*pc(7,12) - 2*pc(7,19) - 4*pc(5,17) + pc(5,12) + pc(5&
&,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(10,16) = 4*pc(11,17) - pc(11,12) - pc(11,19)

sh(11,16) = 4*pc(13,17) - pc(13,12) - pc(13,19) - 4*pc(15,17) + pc(15,12) + pc&
&(15,19)

sh(12,16) = 4*pc(18,17) - pc(18,12) - pc(18,19) - 12*pc(14,17) + 3*pc(14,12) +&
& 3*pc(14,19)

sh(13,16) = 12*pc(12,17) - 3*pc(12,12) - 3*pc(12,19) - 4*pc(19,17) + pc(19,12)&
& + pc(19,19)

sh(14,16) = 8*pc(20,17) - 2*pc(20,12) - 2*pc(20,19) - 12*pc(13,17) + 3*pc(13,1&
&2) + 3*pc(13,19) - 12*pc(15,17) + 3*pc(15,12) + 3*pc(15,19)

sh(15,16) = 16*pc(16,17) - 4*pc(16,12) - 4*pc(16,19) - 4*pc(18,17) + pc(18,12)&
& + pc(18,19) - 4*pc(14,17) + pc(14,12) + pc(14,19)

sh(16,16) = 16*pc(17,17) - 4*pc(17,12) - 4*pc(17,19) - 4*pc(12,17) + pc(12,12)&
& + pc(12,19) - 4*pc(19,17) + pc(19,12) + pc(19,19)


   end subroutine kinetic2CIntgAna

   subroutine kinetic2CIntgNum(a1,a2,A,B,pc,sh,cell_size,step_size)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   real (kind=double), dimension (20,20), intent(out) :: pc
   real (kind=double), dimension (16,16), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i, j, k
   integer :: num_steps
   integer, dimension (20,3) :: triads
   integer, dimension (16,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, curr_pos
   real (kind=double), dimension (3) :: xyz, xyz_soln
   real (kind=double), dimension (3,2) :: xyz_I ! Indices=xyz, noprime||prime

   ! Before we proceed with the calculation we need to understand a bit more
   !   about exactly what is being computed. The form of the integration is:
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        -1/2 * dell^2 [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   ! We use Pxyz for points in space; Rxyz1,2 for atomic sites 1 and 2;
   !   zeta1,2 for the exponential decay factors; and lxyz1,2 for the angular
   !   momentum numbers for atoms 1 and 2.
   ! The dell operator is d^2/dx^2 + d^2/dy^2 + d^2/dz^2.
   ! The integral must be computed over all space for all different possible
   !   values of lxyz1,2 for some arbitrarily chosen test coordinates and
   !   decay rates.
   ! Because of the plus signs in the dell operator we arrive at three
   !   identical integrals of the form (with # in d^2/d#^2 = x, y, or z):
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        -1/2 * d^2/d#^2 [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   ! Now, focusing on just the d^2/dx^2 version of the set of three integrals,
   !   we will pull all terms with Py and Pz out of the dx integral to get:
   ! SS { [(Py-Ry1)**ly1 * (Pz-Rz1)**lz1 *
   !       exp(-zeta1*((Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !      [(Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !       exp(-zeta2*((Py-Ry2)**2 + (Pz-Rz2)**2 ))] *
   !     S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * -1/2 * d^2/dx^2 [
   !      (Px-Rx2)**lx2 * exp(-zeta2*((Px-Rx2)**2))]] dx
   !    } dydz
   ! Applying the derivative, the internal 1D dx integral has the form:
   !   Ix = S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * -1/2 *
   !          [(2*lx2+1)*zeta2*(Px-Rx2)**lx2 - 2*zeta2**2*(Px-Rx2)**lx2**2
   !          -0.5*lx2*(lx2-1)*(Px-Rx2)**(lx2-2)] * exp(-zeta2*(Px-Rx2)**2)
   !        ]
   ! Each of the other integrals (Iy', Iz') will have the form of a simple
   !   1D overlap integral:
   !   Iy' = S [(Py-Ry1)**ly1 * exp(-zeta1*((Py-Ry1)**2)) *
   !            (Py-Ry2)**ly2 * exp(-zeta2*((Py-Ry2)**2))] dy
   !   Iz' = S [(Pz-Rz1)**lz1 * exp(-zeta1*((Pz-Rz1)**2)) *
   !            (Pz-Rz2)**lz2 * exp(-zeta2*((Pz-Rz2)**2))] dz
   ! The total integral !! for the d^2/dx^2 version !! is thus Ix * Iy' * Iz'.
   !   The total integral (including all terms of dell^2) will have the form:
   !   Ix*Iy'*Iz' + Ix'*Iy*Iz' + Ix'*Iy'*Iz where the Iy and Iz are the
   !   appropriate analogs of the Ix and the Ix' is the analog of the Iy'
   !   or Iz'.
   ! Thus, every integral is an independent 1D integral and the solution of
   !   the total integral is a sum of products of integrals that can all be
   !   computed at once.
   ! With regard to the term in each integral that has the **(lx2-2) form
   !   (or **(ly2-2) or **(lz2-2)), some special care must be taken. This
   !   term represents an angular momentum and the integral will have the
   !   form of an overlap integral. Because we cannot have a negative angular
   !   momentum we must discard this term when lx2-2 < 0.

   ! Initialize local variables.
   triads(1,:) = (/0,0,0/)
   triads(2,:) = (/1,0,0/)
   triads(3,:) = (/0,1,0/)
   triads(4,:) = (/0,0,1/)
   triads(5,:) = (/2,0,0/)
   triads(6,:) = (/0,2,0/)
   triads(7,:) = (/0,0,2/)
   triads(8,:) = (/1,1,0/)
   triads(9,:) = (/1,0,1/)
   triads(10,:) = (/0,1,1/)
   triads(11,:) = (/1,1,1/)
   triads(12,:) = (/2,1,0/)
   triads(13,:) = (/2,0,1/)
   triads(14,:) = (/1,2,0/)
   triads(15,:) = (/0,2,1/)
   triads(16,:) = (/1,0,2/)
   triads(17,:) = (/0,1,2/)
   triads(18,:) = (/3,0,0/)
   triads(19,:) = (/0,3,0/)
   triads(20,:) = (/0,0,3/)
   conversion(1,1,:) = (/1,0,0/)
   conversion(1,2,:) = (/1,1,1/)
   conversion(2,1,:) = (/1,0,0/)
   conversion(2,2,:) = (/2,1,1/)
   conversion(3,1,:) = (/1,0,0/)
   conversion(3,2,:) = (/3,1,1/)
   conversion(4,1,:) = (/1,0,0/)
   conversion(4,2,:) = (/4,1,1/)
   conversion(5,1,:) = (/1,0,0/)
   conversion(5,2,:) = (/8,1,1/)
   conversion(6,1,:) = (/1,0,0/)
   conversion(6,2,:) = (/9,1,1/)
   conversion(7,1,:) = (/1,0,0/)
   conversion(7,2,:) = (/10,1,1/)
   conversion(8,1,:) = (/1,-1,0/)
   conversion(8,2,:) = (/5,6,1/)
   conversion(9,1,:) = (/2,-1,-1/)
   conversion(9,2,:) = (/7,5,6/)
   conversion(10,1,:) = (/1,0,0/)
   conversion(10,2,:) = (/11,1,1/)
   conversion(11,1,:) = (/1,-1,0/)
   conversion(11,2,:) = (/13,15,1/)
   conversion(12,1,:) = (/1,-3,0/)
   conversion(12,2,:) = (/18,14,1/)
   conversion(13,1,:) = (/3,-1,0/)
   conversion(13,2,:) = (/12,19,1/)
   conversion(14,1,:) = (/2,-3,-3/)
   conversion(14,2,:) = (/20,13,15/)
   conversion(15,1,:) = (/4,-1,-1/)
   conversion(15,2,:) = (/16,18,14/)
   conversion(16,1,:) = (/4,-1,-1/)
   conversion(16,2,:) = (/17,12,19/)

   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   do p = 1, 20
      do q = 1, 20

         ! Assign l1 and l2 values for each primitive Cartesian gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize sum variables.
         xyz_I(:,:) = 0.0d0

         ! Start a loop over x coordinates.
         do i = 0, num_steps
            curr_pos = start_pos + (i*step_size)

            ! Compute the no-prime integrals first.
            do j = 1, 3
               xyz_I(j,1) = xyz_I(j,1) &
                     & + noPrimeKE(step_size, curr_pos, A(j), &
                     & B(j), a1, a2, l1(j), l2(j))
            enddo

            ! Compute the primed integrals second.
            do j = 1, 3
               xyz_I(j,2) = xyz_I(j,2) &
                     & + primeKE(step_size, curr_pos, A(j), &
                     & B(j), a1, a2, l1(j), l2(j))
            enddo
         enddo

         pc(q,p) = xyz_I(1,1)*xyz_I(2,2)*xyz_I(3,2) + &
               & xyz_I(1,2)*xyz_I(2,1)*xyz_I(3,2) + &
               & xyz_I(1,2)*xyz_I(2,2)*xyz_I(3,1)
      enddo
    enddo

sh(1,1) = pc(1,1)

sh(2,1) = pc(2,1)

sh(3,1) = pc(3,1)

sh(4,1) = pc(4,1)

sh(5,1) = pc(8,1)

sh(6,1) = pc(9,1)

sh(7,1) = pc(10,1)

sh(8,1) = pc(5,1) - pc(6,1)

sh(9,1) = 2*pc(7,1) - pc(5,1) - pc(6,1)

sh(10,1) = pc(11,1)

sh(11,1) = pc(13,1) - pc(15,1)

sh(12,1) = pc(18,1) - 3*pc(14,1)

sh(13,1) = 3*pc(12,1) - pc(19,1)

sh(14,1) = 2*pc(20,1) - 3*pc(13,1) - 3*pc(15,1)

sh(15,1) = 4*pc(16,1) - pc(18,1) - pc(14,1)

sh(16,1) = 4*pc(17,1) - pc(12,1) - pc(19,1)

sh(1,2) = pc(1,2)

sh(2,2) = pc(2,2)

sh(3,2) = pc(3,2)

sh(4,2) = pc(4,2)

sh(5,2) = pc(8,2)

sh(6,2) = pc(9,2)

sh(7,2) = pc(10,2)

sh(8,2) = pc(5,2) - pc(6,2)

sh(9,2) = 2*pc(7,2) - pc(5,2) - pc(6,2)

sh(10,2) = pc(11,2)

sh(11,2) = pc(13,2) - pc(15,2)

sh(12,2) = pc(18,2) - 3*pc(14,2)

sh(13,2) = 3*pc(12,2) - pc(19,2)

sh(14,2) = 2*pc(20,2) - 3*pc(13,2) - 3*pc(15,2)

sh(15,2) = 4*pc(16,2) - pc(18,2) - pc(14,2)

sh(16,2) = 4*pc(17,2) - pc(12,2) - pc(19,2)

sh(1,3) = pc(1,3)

sh(2,3) = pc(2,3)

sh(3,3) = pc(3,3)

sh(4,3) = pc(4,3)

sh(5,3) = pc(8,3)

sh(6,3) = pc(9,3)

sh(7,3) = pc(10,3)

sh(8,3) = pc(5,3) - pc(6,3)

sh(9,3) = 2*pc(7,3) - pc(5,3) - pc(6,3)

sh(10,3) = pc(11,3)

sh(11,3) = pc(13,3) - pc(15,3)

sh(12,3) = pc(18,3) - 3*pc(14,3)

sh(13,3) = 3*pc(12,3) - pc(19,3)

sh(14,3) = 2*pc(20,3) - 3*pc(13,3) - 3*pc(15,3)

sh(15,3) = 4*pc(16,3) - pc(18,3) - pc(14,3)

sh(16,3) = 4*pc(17,3) - pc(12,3) - pc(19,3)

sh(1,4) = pc(1,4)

sh(2,4) = pc(2,4)

sh(3,4) = pc(3,4)

sh(4,4) = pc(4,4)

sh(5,4) = pc(8,4)

sh(6,4) = pc(9,4)

sh(7,4) = pc(10,4)

sh(8,4) = pc(5,4) - pc(6,4)

sh(9,4) = 2*pc(7,4) - pc(5,4) - pc(6,4)

sh(10,4) = pc(11,4)

sh(11,4) = pc(13,4) - pc(15,4)

sh(12,4) = pc(18,4) - 3*pc(14,4)

sh(13,4) = 3*pc(12,4) - pc(19,4)

sh(14,4) = 2*pc(20,4) - 3*pc(13,4) - 3*pc(15,4)

sh(15,4) = 4*pc(16,4) - pc(18,4) - pc(14,4)

sh(16,4) = 4*pc(17,4) - pc(12,4) - pc(19,4)

sh(1,5) = pc(1,8)

sh(2,5) = pc(2,8)

sh(3,5) = pc(3,8)

sh(4,5) = pc(4,8)

sh(5,5) = pc(8,8)

sh(6,5) = pc(9,8)

sh(7,5) = pc(10,8)

sh(8,5) = pc(5,8) - pc(6,8)

sh(9,5) = 2*pc(7,8) - pc(5,8) - pc(6,8)

sh(10,5) = pc(11,8)

sh(11,5) = pc(13,8) - pc(15,8)

sh(12,5) = pc(18,8) - 3*pc(14,8)

sh(13,5) = 3*pc(12,8) - pc(19,8)

sh(14,5) = 2*pc(20,8) - 3*pc(13,8) - 3*pc(15,8)

sh(15,5) = 4*pc(16,8) - pc(18,8) - pc(14,8)

sh(16,5) = 4*pc(17,8) - pc(12,8) - pc(19,8)

sh(1,6) = pc(1,9)

sh(2,6) = pc(2,9)

sh(3,6) = pc(3,9)

sh(4,6) = pc(4,9)

sh(5,6) = pc(8,9)

sh(6,6) = pc(9,9)

sh(7,6) = pc(10,9)

sh(8,6) = pc(5,9) - pc(6,9)

sh(9,6) = 2*pc(7,9) - pc(5,9) - pc(6,9)

sh(10,6) = pc(11,9)

sh(11,6) = pc(13,9) - pc(15,9)

sh(12,6) = pc(18,9) - 3*pc(14,9)

sh(13,6) = 3*pc(12,9) - pc(19,9)

sh(14,6) = 2*pc(20,9) - 3*pc(13,9) - 3*pc(15,9)

sh(15,6) = 4*pc(16,9) - pc(18,9) - pc(14,9)

sh(16,6) = 4*pc(17,9) - pc(12,9) - pc(19,9)

sh(1,7) = pc(1,10)

sh(2,7) = pc(2,10)

sh(3,7) = pc(3,10)

sh(4,7) = pc(4,10)

sh(5,7) = pc(8,10)

sh(6,7) = pc(9,10)

sh(7,7) = pc(10,10)

sh(8,7) = pc(5,10) - pc(6,10)

sh(9,7) = 2*pc(7,10) - pc(5,10) - pc(6,10)

sh(10,7) = pc(11,10)

sh(11,7) = pc(13,10) - pc(15,10)

sh(12,7) = pc(18,10) - 3*pc(14,10)

sh(13,7) = 3*pc(12,10) - pc(19,10)

sh(14,7) = 2*pc(20,10) - 3*pc(13,10) - 3*pc(15,10)

sh(15,7) = 4*pc(16,10) - pc(18,10) - pc(14,10)

sh(16,7) = 4*pc(17,10) - pc(12,10) - pc(19,10)

sh(1,8) = pc(1,5) - pc(1,6)

sh(2,8) = pc(2,5) - pc(2,6)

sh(3,8) = pc(3,5) - pc(3,6)

sh(4,8) = pc(4,5) - pc(4,6)

sh(5,8) = pc(8,5) - pc(8,6)

sh(6,8) = pc(9,5) - pc(9,6)

sh(7,8) = pc(10,5) - pc(10,6)

sh(8,8) = pc(5,5) - pc(5,6) - pc(6,5) + pc(6,6)

sh(9,8) = 2*pc(7,5) - 2*pc(7,6) - pc(5,5) + pc(5,6) - pc(6,5) + pc(6,6)

sh(10,8) = pc(11,5) - pc(11,6)

sh(11,8) = pc(13,5) - pc(13,6) - pc(15,5) + pc(15,6)

sh(12,8) = pc(18,5) - pc(18,6) - 3*pc(14,5) + 3*pc(14,6)

sh(13,8) = 3*pc(12,5) - 3*pc(12,6) - pc(19,5) + pc(19,6)

sh(14,8) = 2*pc(20,5) - 2*pc(20,6) - 3*pc(13,5) + 3*pc(13,6) - 3*pc(15,5) + 3*&
&pc(15,6)

sh(15,8) = 4*pc(16,5) - 4*pc(16,6) - pc(18,5) + pc(18,6) - pc(14,5) + pc(14,6)

sh(16,8) = 4*pc(17,5) - 4*pc(17,6) - pc(12,5) + pc(12,6) - pc(19,5) + pc(19,6)

sh(1,9) = 2*pc(1,7) - pc(1,5) - pc(1,6)

sh(2,9) = 2*pc(2,7) - pc(2,5) - pc(2,6)

sh(3,9) = 2*pc(3,7) - pc(3,5) - pc(3,6)

sh(4,9) = 2*pc(4,7) - pc(4,5) - pc(4,6)

sh(5,9) = 2*pc(8,7) - pc(8,5) - pc(8,6)

sh(6,9) = 2*pc(9,7) - pc(9,5) - pc(9,6)

sh(7,9) = 2*pc(10,7) - pc(10,5) - pc(10,6)

sh(8,9) = 2*pc(5,7) - pc(5,5) - pc(5,6) - 2*pc(6,7) + pc(6,5) + pc(6,6)

sh(9,9) = 4*pc(7,7) - 2*pc(7,5) - 2*pc(7,6) - 2*pc(5,7) + pc(5,5) + pc(5,6) - &
&2*pc(6,7) + pc(6,5) + pc(6,6)

sh(10,9) = 2*pc(11,7) - pc(11,5) - pc(11,6)

sh(11,9) = 2*pc(13,7) - pc(13,5) - pc(13,6) - 2*pc(15,7) + pc(15,5) + pc(15,6)

sh(12,9) = 2*pc(18,7) - pc(18,5) - pc(18,6) - 6*pc(14,7) + 3*pc(14,5) + 3*pc(1&
&4,6)

sh(13,9) = 6*pc(12,7) - 3*pc(12,5) - 3*pc(12,6) - 2*pc(19,7) + pc(19,5) + pc(1&
&9,6)

sh(14,9) = 4*pc(20,7) - 2*pc(20,5) - 2*pc(20,6) - 6*pc(13,7) + 3*pc(13,5) + 3*&
&pc(13,6) - 6*pc(15,7) + 3*pc(15,5) + 3*pc(15,6)

sh(15,9) = 8*pc(16,7) - 4*pc(16,5) - 4*pc(16,6) - 2*pc(18,7) + pc(18,5) + pc(1&
&8,6) - 2*pc(14,7) + pc(14,5) + pc(14,6)

sh(16,9) = 8*pc(17,7) - 4*pc(17,5) - 4*pc(17,6) - 2*pc(12,7) + pc(12,5) + pc(1&
&2,6) - 2*pc(19,7) + pc(19,5) + pc(19,6)

sh(1,10) = pc(1,11)

sh(2,10) = pc(2,11)

sh(3,10) = pc(3,11)

sh(4,10) = pc(4,11)

sh(5,10) = pc(8,11)

sh(6,10) = pc(9,11)

sh(7,10) = pc(10,11)

sh(8,10) = pc(5,11) - pc(6,11)

sh(9,10) = 2*pc(7,11) - pc(5,11) - pc(6,11)

sh(10,10) = pc(11,11)

sh(11,10) = pc(13,11) - pc(15,11)

sh(12,10) = pc(18,11) - 3*pc(14,11)

sh(13,10) = 3*pc(12,11) - pc(19,11)

sh(14,10) = 2*pc(20,11) - 3*pc(13,11) - 3*pc(15,11)

sh(15,10) = 4*pc(16,11) - pc(18,11) - pc(14,11)

sh(16,10) = 4*pc(17,11) - pc(12,11) - pc(19,11)

sh(1,11) = pc(1,13) - pc(1,15)

sh(2,11) = pc(2,13) - pc(2,15)

sh(3,11) = pc(3,13) - pc(3,15)

sh(4,11) = pc(4,13) - pc(4,15)

sh(5,11) = pc(8,13) - pc(8,15)

sh(6,11) = pc(9,13) - pc(9,15)

sh(7,11) = pc(10,13) - pc(10,15)

sh(8,11) = pc(5,13) - pc(5,15) - pc(6,13) + pc(6,15)

sh(9,11) = 2*pc(7,13) - 2*pc(7,15) - pc(5,13) + pc(5,15) - pc(6,13) + pc(6,15)

sh(10,11) = pc(11,13) - pc(11,15)

sh(11,11) = pc(13,13) - pc(13,15) - pc(15,13) + pc(15,15)

sh(12,11) = pc(18,13) - pc(18,15) - 3*pc(14,13) + 3*pc(14,15)

sh(13,11) = 3*pc(12,13) - 3*pc(12,15) - pc(19,13) + pc(19,15)

sh(14,11) = 2*pc(20,13) - 2*pc(20,15) - 3*pc(13,13) + 3*pc(13,15) - 3*pc(15,13&
&) + 3*pc(15,15)

sh(15,11) = 4*pc(16,13) - 4*pc(16,15) - pc(18,13) + pc(18,15) - pc(14,13) + pc&
&(14,15)

sh(16,11) = 4*pc(17,13) - 4*pc(17,15) - pc(12,13) + pc(12,15) - pc(19,13) + pc&
&(19,15)

sh(1,12) = pc(1,18) - 3*pc(1,14)

sh(2,12) = pc(2,18) - 3*pc(2,14)

sh(3,12) = pc(3,18) - 3*pc(3,14)

sh(4,12) = pc(4,18) - 3*pc(4,14)

sh(5,12) = pc(8,18) - 3*pc(8,14)

sh(6,12) = pc(9,18) - 3*pc(9,14)

sh(7,12) = pc(10,18) - 3*pc(10,14)

sh(8,12) = pc(5,18) - 3*pc(5,14) - pc(6,18) + 3*pc(6,14)

sh(9,12) = 2*pc(7,18) - 6*pc(7,14) - pc(5,18) + 3*pc(5,14) - pc(6,18) + 3*pc(6&
&,14)

sh(10,12) = pc(11,18) - 3*pc(11,14)

sh(11,12) = pc(13,18) - 3*pc(13,14) - pc(15,18) + 3*pc(15,14)

sh(12,12) = pc(18,18) - 3*pc(18,14) - 3*pc(14,18) + 9*pc(14,14)

sh(13,12) = 3*pc(12,18) - 9*pc(12,14) - pc(19,18) + 3*pc(19,14)

sh(14,12) = 2*pc(20,18) - 6*pc(20,14) - 3*pc(13,18) + 9*pc(13,14) - 3*pc(15,18&
&) + 9*pc(15,14)

sh(15,12) = 4*pc(16,18) - 12*pc(16,14) - pc(18,18) + 3*pc(18,14) - pc(14,18) +&
& 3*pc(14,14)

sh(16,12) = 4*pc(17,18) - 12*pc(17,14) - pc(12,18) + 3*pc(12,14) - pc(19,18) +&
& 3*pc(19,14)

sh(1,13) = 3*pc(1,12) - pc(1,19)

sh(2,13) = 3*pc(2,12) - pc(2,19)

sh(3,13) = 3*pc(3,12) - pc(3,19)

sh(4,13) = 3*pc(4,12) - pc(4,19)

sh(5,13) = 3*pc(8,12) - pc(8,19)

sh(6,13) = 3*pc(9,12) - pc(9,19)

sh(7,13) = 3*pc(10,12) - pc(10,19)

sh(8,13) = 3*pc(5,12) - pc(5,19) - 3*pc(6,12) + pc(6,19)

sh(9,13) = 6*pc(7,12) - 2*pc(7,19) - 3*pc(5,12) + pc(5,19) - 3*pc(6,12) + pc(6&
&,19)

sh(10,13) = 3*pc(11,12) - pc(11,19)

sh(11,13) = 3*pc(13,12) - pc(13,19) - 3*pc(15,12) + pc(15,19)

sh(12,13) = 3*pc(18,12) - pc(18,19) - 9*pc(14,12) + 3*pc(14,19)

sh(13,13) = 9*pc(12,12) - 3*pc(12,19) - 3*pc(19,12) + pc(19,19)

sh(14,13) = 6*pc(20,12) - 2*pc(20,19) - 9*pc(13,12) + 3*pc(13,19) - 9*pc(15,12&
&) + 3*pc(15,19)

sh(15,13) = 12*pc(16,12) - 4*pc(16,19) - 3*pc(18,12) + pc(18,19) - 3*pc(14,12)&
& + pc(14,19)

sh(16,13) = 12*pc(17,12) - 4*pc(17,19) - 3*pc(12,12) + pc(12,19) - 3*pc(19,12)&
& + pc(19,19)

sh(1,14) = 2*pc(1,20) - 3*pc(1,13) - 3*pc(1,15)

sh(2,14) = 2*pc(2,20) - 3*pc(2,13) - 3*pc(2,15)

sh(3,14) = 2*pc(3,20) - 3*pc(3,13) - 3*pc(3,15)

sh(4,14) = 2*pc(4,20) - 3*pc(4,13) - 3*pc(4,15)

sh(5,14) = 2*pc(8,20) - 3*pc(8,13) - 3*pc(8,15)

sh(6,14) = 2*pc(9,20) - 3*pc(9,13) - 3*pc(9,15)

sh(7,14) = 2*pc(10,20) - 3*pc(10,13) - 3*pc(10,15)

sh(8,14) = 2*pc(5,20) - 3*pc(5,13) - 3*pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*&
&pc(6,15)

sh(9,14) = 4*pc(7,20) - 6*pc(7,13) - 6*pc(7,15) - 2*pc(5,20) + 3*pc(5,13) + 3*&
&pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*pc(6,15)

sh(10,14) = 2*pc(11,20) - 3*pc(11,13) - 3*pc(11,15)

sh(11,14) = 2*pc(13,20) - 3*pc(13,13) - 3*pc(13,15) - 2*pc(15,20) + 3*pc(15,13&
&) + 3*pc(15,15)

sh(12,14) = 2*pc(18,20) - 3*pc(18,13) - 3*pc(18,15) - 6*pc(14,20) + 9*pc(14,13&
&) + 9*pc(14,15)

sh(13,14) = 6*pc(12,20) - 9*pc(12,13) - 9*pc(12,15) - 2*pc(19,20) + 3*pc(19,13&
&) + 3*pc(19,15)

sh(14,14) = 4*pc(20,20) - 6*pc(20,13) - 6*pc(20,15) - 6*pc(13,20) + 9*pc(13,13&
&) + 9*pc(13,15) - 6*pc(15,20) + 9*pc(15,13) + 9*pc(15,15)

sh(15,14) = 8*pc(16,20) - 12*pc(16,13) - 12*pc(16,15) - 2*pc(18,20) + 3*pc(18,&
&13) + 3*pc(18,15) - 2*pc(14,20) + 3*pc(14,13) + 3*pc(14,15)

sh(16,14) = 8*pc(17,20) - 12*pc(17,13) - 12*pc(17,15) - 2*pc(12,20) + 3*pc(12,&
&13) + 3*pc(12,15) - 2*pc(19,20) + 3*pc(19,13) + 3*pc(19,15)

sh(1,15) = 4*pc(1,16) - pc(1,18) - pc(1,14)

sh(2,15) = 4*pc(2,16) - pc(2,18) - pc(2,14)

sh(3,15) = 4*pc(3,16) - pc(3,18) - pc(3,14)

sh(4,15) = 4*pc(4,16) - pc(4,18) - pc(4,14)

sh(5,15) = 4*pc(8,16) - pc(8,18) - pc(8,14)

sh(6,15) = 4*pc(9,16) - pc(9,18) - pc(9,14)

sh(7,15) = 4*pc(10,16) - pc(10,18) - pc(10,14)

sh(8,15) = 4*pc(5,16) - pc(5,18) - pc(5,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(9,15) = 8*pc(7,16) - 2*pc(7,18) - 2*pc(7,14) - 4*pc(5,16) + pc(5,18) + pc(5&
&,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(10,15) = 4*pc(11,16) - pc(11,18) - pc(11,14)

sh(11,15) = 4*pc(13,16) - pc(13,18) - pc(13,14) - 4*pc(15,16) + pc(15,18) + pc&
&(15,14)

sh(12,15) = 4*pc(18,16) - pc(18,18) - pc(18,14) - 12*pc(14,16) + 3*pc(14,18) +&
& 3*pc(14,14)

sh(13,15) = 12*pc(12,16) - 3*pc(12,18) - 3*pc(12,14) - 4*pc(19,16) + pc(19,18)&
& + pc(19,14)

sh(14,15) = 8*pc(20,16) - 2*pc(20,18) - 2*pc(20,14) - 12*pc(13,16) + 3*pc(13,1&
&8) + 3*pc(13,14) - 12*pc(15,16) + 3*pc(15,18) + 3*pc(15,14)

sh(15,15) = 16*pc(16,16) - 4*pc(16,18) - 4*pc(16,14) - 4*pc(18,16) + pc(18,18)&
& + pc(18,14) - 4*pc(14,16) + pc(14,18) + pc(14,14)

sh(16,15) = 16*pc(17,16) - 4*pc(17,18) - 4*pc(17,14) - 4*pc(12,16) + pc(12,18)&
& + pc(12,14) - 4*pc(19,16) + pc(19,18) + pc(19,14)

sh(1,16) = 4*pc(1,17) - pc(1,12) - pc(1,19)

sh(2,16) = 4*pc(2,17) - pc(2,12) - pc(2,19)

sh(3,16) = 4*pc(3,17) - pc(3,12) - pc(3,19)

sh(4,16) = 4*pc(4,17) - pc(4,12) - pc(4,19)

sh(5,16) = 4*pc(8,17) - pc(8,12) - pc(8,19)

sh(6,16) = 4*pc(9,17) - pc(9,12) - pc(9,19)

sh(7,16) = 4*pc(10,17) - pc(10,12) - pc(10,19)

sh(8,16) = 4*pc(5,17) - pc(5,12) - pc(5,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(9,16) = 8*pc(7,17) - 2*pc(7,12) - 2*pc(7,19) - 4*pc(5,17) + pc(5,12) + pc(5&
&,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(10,16) = 4*pc(11,17) - pc(11,12) - pc(11,19)

sh(11,16) = 4*pc(13,17) - pc(13,12) - pc(13,19) - 4*pc(15,17) + pc(15,12) + pc&
&(15,19)

sh(12,16) = 4*pc(18,17) - pc(18,12) - pc(18,19) - 12*pc(14,17) + 3*pc(14,12) +&
& 3*pc(14,19)

sh(13,16) = 12*pc(12,17) - 3*pc(12,12) - 3*pc(12,19) - 4*pc(19,17) + pc(19,12)&
& + pc(19,19)

sh(14,16) = 8*pc(20,17) - 2*pc(20,12) - 2*pc(20,19) - 12*pc(13,17) + 3*pc(13,1&
&2) + 3*pc(13,19) - 12*pc(15,17) + 3*pc(15,12) + 3*pc(15,19)

sh(15,16) = 16*pc(16,17) - 4*pc(16,12) - 4*pc(16,19) - 4*pc(18,17) + pc(18,12)&
& + pc(18,19) - 4*pc(14,17) + pc(14,12) + pc(14,19)

sh(16,16) = 16*pc(17,17) - 4*pc(17,12) - 4*pc(17,19) - 4*pc(12,17) + pc(12,12)&
& + pc(12,19) - 4*pc(19,17) + pc(19,12) + pc(19,19)


   end subroutine kinetic2CIntgNum


   function primeKE(step_size, curr_pos, A, B, a1, a2, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos, A, B, a1, a2
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: primeKE

      ! Compute the integral part.
      primeKE = step_size * (curr_pos - A)**l1 * (curr_pos - B)**l2 &
            & * exp(-a1*(curr_pos - A)**2) * exp(-a2*(curr_pos - B)**2)

      return

   end function primeKE


   function noPrimeKE(step_size, curr_pos, A, B, a1, a2, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos, A, B, a1, a2
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: noPrimeKE

      ! Compute each "internal" term of the no-prime integral. Compare each
      !   of these lines with the 1D kinetic energy equation produce by the
      !   osrecurintg_makenum.py script (appropriately separated into terms)
      !   and the expression in the last equals of equation 65 in Ben Walker's
      !   dissertation. Note that a factor of -1/2 is included in all of the
      !   below expressions.

      ! Combination of terms 1 and 2 in line 1 of last equals of eqn 65. Also
      !   visible as combination of first two terms in the main parentheses
      !   of the algebraically ordered 1D kinetic energy equation produced by
      !   osrecurintg_makenum.py.
      noPrimeKE = (a2 * (2*l2 + 1)) * (curr_pos - B)**l2

      ! Line 2 of last equals in eqn 65. Also visible as third term in the
      !   script produced equation.
      noPrimeKE = noPrimeKE - 2.0d0*a2**2 * (curr_pos - B)**(l2+2)

      ! Combination of terms 1 and 2 in line 3. Also visible as the
      !   combination of the 4th and 5th terms in the script produced eqn.
      !   As mentioned above, we only compute this if the angular momentum
      !   will allow it.
      if (l2 >= 2) then
         noPrimeKE = noPrimeKE - 0.5 * l2*(l2-1) * (curr_pos - B)**(l2-2)
      endif

      ! Multiply the no-prime integral by the preceeding primitive gaussian
      !   coefficient and exponential and multiply by the succeeding
      !   exponential. (We have already multiplied by the succeeding
      !   primitive gaussian coefficient in the above lines.)
      noPrimeKE = noPrimeKE * (curr_pos-A)**l1 * exp(-a1*(curr_pos-A)**2) &
            & * exp(-a2*(curr_pos-B)**2)

      ! Finally, multiply by the step size.
      noPrimeKE = noPrimeKE * step_size

      return
      
   end function noPrimeKE

   subroutine electron3CIntgAna(a1,a2,a3,A,B,C,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2, a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   real (kind=double), dimension (20,20), intent(out) :: pc
   real (kind=double), dimension (16,16), intent(out) :: sh

   ! Define local variables.
   real (kind=double), dimension (3) :: G, GA, GB, P, PC_3C, d
   real (kind=double) :: zeta, zeta3C, inv_2zeta3C, xi
   real (kind=double) :: preFactorOL

   ! Initialize local variables.
   zeta = a1 + a2
   zeta3C = zeta + a3
   inv_2zeta3C = 1.0d0 / (2.0d0 * zeta3C)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   G = (a1*A + a2*B + a3*C) / zeta3C
   GA = G - A
   GB = G - B
   PC_3C = P - C
   d = A - B
   prefactorOL = ((pi/zeta3C)**1.5) &
         & * exp(-xi*sum(d**2)-(zeta*a3/zeta3C)*sum(PC_3C**2))

pc(1,1) = preFactorOL

pc(2,1) = GA(1)*preFactorOL

pc(3,1) = GA(2)*preFactorOL

pc(4,1) = GA(3)*preFactorOL

pc(5,1) = GA(1)*pc(2,1) + inv_2zeta3C*(preFactorOL)

pc(6,1) = GA(2)*pc(3,1) + inv_2zeta3C*(preFactorOL)

pc(7,1) = GA(3)*pc(4,1) + inv_2zeta3C*(preFactorOL)

pc(8,1) = GA(2)*pc(2,1)

pc(9,1) = GA(3)*pc(2,1)

pc(10,1) = GA(3)*pc(3,1)

pc(11,1) = GA(3)*pc(8,1)

pc(12,1) = GA(2)*pc(5,1)

pc(13,1) = GA(3)*pc(5,1)

pc(14,1) = GA(1)*pc(6,1)

pc(15,1) = GA(3)*pc(6,1)

pc(16,1) = GA(1)*pc(7,1)

pc(17,1) = GA(2)*pc(7,1)

pc(18,1) = GA(1)*pc(5,1) + inv_2zeta3C*(2*pc(2,1))

pc(19,1) = GA(2)*pc(6,1) + inv_2zeta3C*(2*pc(3,1))

pc(20,1) = GA(3)*pc(7,1) + inv_2zeta3C*(2*pc(4,1))

pc(1,2) = GB(1)*preFactorOL

pc(2,2) = GA(1)*pc(1,2) + inv_2zeta3C*(preFactorOL)

pc(3,2) = GA(2)*pc(1,2)

pc(4,2) = GA(3)*pc(1,2)

pc(5,2) = GA(1)*pc(2,2) + inv_2zeta3C*(pc(1,2) + pc(2,1))

pc(6,2) = GA(2)*pc(3,2) + inv_2zeta3C*(pc(1,2))

pc(7,2) = GA(3)*pc(4,2) + inv_2zeta3C*(pc(1,2))

pc(8,2) = GA(2)*pc(2,2)

pc(9,2) = GA(3)*pc(2,2)

pc(10,2) = GA(3)*pc(3,2)

pc(11,2) = GA(3)*pc(8,2)

pc(12,2) = GA(2)*pc(5,2)

pc(13,2) = GA(3)*pc(5,2)

pc(14,2) = GA(1)*pc(6,2) + inv_2zeta3C*(pc(6,1))

pc(15,2) = GA(3)*pc(6,2)

pc(16,2) = GA(1)*pc(7,2) + inv_2zeta3C*(pc(7,1))

pc(17,2) = GA(2)*pc(7,2)

pc(18,2) = GA(1)*pc(5,2) + inv_2zeta3C*(2*pc(2,2) + pc(5,1))

pc(19,2) = GA(2)*pc(6,2) + inv_2zeta3C*(2*pc(3,2))

pc(20,2) = GA(3)*pc(7,2) + inv_2zeta3C*(2*pc(4,2))

pc(1,3) = GB(2)*preFactorOL

pc(2,3) = GA(1)*pc(1,3)

pc(3,3) = GA(2)*pc(1,3) + inv_2zeta3C*(preFactorOL)

pc(4,3) = GA(3)*pc(1,3)

pc(5,3) = GA(1)*pc(2,3) + inv_2zeta3C*(pc(1,3))

pc(6,3) = GA(2)*pc(3,3) + inv_2zeta3C*(pc(1,3) + pc(3,1))

pc(7,3) = GA(3)*pc(4,3) + inv_2zeta3C*(pc(1,3))

pc(8,3) = GA(2)*pc(2,3) + inv_2zeta3C*(pc(2,1))

pc(9,3) = GA(3)*pc(2,3)

pc(10,3) = GA(3)*pc(3,3)

pc(11,3) = GA(3)*pc(8,3)

pc(12,3) = GA(2)*pc(5,3) + inv_2zeta3C*(pc(5,1))

pc(13,3) = GA(3)*pc(5,3)

pc(14,3) = GA(1)*pc(6,3)

pc(15,3) = GA(3)*pc(6,3)

pc(16,3) = GA(1)*pc(7,3)

pc(17,3) = GA(2)*pc(7,3) + inv_2zeta3C*(pc(7,1))

pc(18,3) = GA(1)*pc(5,3) + inv_2zeta3C*(2*pc(2,3))

pc(19,3) = GA(2)*pc(6,3) + inv_2zeta3C*(2*pc(3,3) + pc(6,1))

pc(20,3) = GA(3)*pc(7,3) + inv_2zeta3C*(2*pc(4,3))

pc(1,4) = GB(3)*preFactorOL

pc(2,4) = GA(1)*pc(1,4)

pc(3,4) = GA(2)*pc(1,4)

pc(4,4) = GA(3)*pc(1,4) + inv_2zeta3C*(preFactorOL)

pc(5,4) = GA(1)*pc(2,4) + inv_2zeta3C*(pc(1,4))

pc(6,4) = GA(2)*pc(3,4) + inv_2zeta3C*(pc(1,4))

pc(7,4) = GA(3)*pc(4,4) + inv_2zeta3C*(pc(1,4) + pc(4,1))

pc(8,4) = GA(2)*pc(2,4)

pc(9,4) = GA(3)*pc(2,4) + inv_2zeta3C*(pc(2,1))

pc(10,4) = GA(3)*pc(3,4) + inv_2zeta3C*(pc(3,1))

pc(11,4) = GA(3)*pc(8,4) + inv_2zeta3C*(pc(8,1))

pc(12,4) = GA(2)*pc(5,4)

pc(13,4) = GA(3)*pc(5,4) + inv_2zeta3C*(pc(5,1))

pc(14,4) = GA(1)*pc(6,4)

pc(15,4) = GA(3)*pc(6,4) + inv_2zeta3C*(pc(6,1))

pc(16,4) = GA(1)*pc(7,4)

pc(17,4) = GA(2)*pc(7,4)

pc(18,4) = GA(1)*pc(5,4) + inv_2zeta3C*(2*pc(2,4))

pc(19,4) = GA(2)*pc(6,4) + inv_2zeta3C*(2*pc(3,4))

pc(20,4) = GA(3)*pc(7,4) + inv_2zeta3C*(2*pc(4,4) + pc(7,1))

pc(1,5) = GB(1)*pc(1,2) + inv_2zeta3C*(preFactorOL)

pc(2,5) = GA(1)*pc(1,5) + inv_2zeta3C*(2*pc(1,2))

pc(3,5) = GA(2)*pc(1,5)

pc(4,5) = GA(3)*pc(1,5)

pc(5,5) = GA(1)*pc(2,5) + inv_2zeta3C*(pc(1,5) + 2*pc(2,2))

pc(6,5) = GA(2)*pc(3,5) + inv_2zeta3C*(pc(1,5))

pc(7,5) = GA(3)*pc(4,5) + inv_2zeta3C*(pc(1,5))

pc(8,5) = GA(2)*pc(2,5)

pc(9,5) = GA(3)*pc(2,5)

pc(10,5) = GA(3)*pc(3,5)

pc(11,5) = GA(3)*pc(8,5)

pc(12,5) = GA(2)*pc(5,5)

pc(13,5) = GA(3)*pc(5,5)

pc(14,5) = GA(1)*pc(6,5) + inv_2zeta3C*(2*pc(6,2))

pc(15,5) = GA(3)*pc(6,5)

pc(16,5) = GA(1)*pc(7,5) + inv_2zeta3C*(2*pc(7,2))

pc(17,5) = GA(2)*pc(7,5)

pc(18,5) = GA(1)*pc(5,5) + inv_2zeta3C*(2*pc(2,5) + 2*pc(5,2))

pc(19,5) = GA(2)*pc(6,5) + inv_2zeta3C*(2*pc(3,5))

pc(20,5) = GA(3)*pc(7,5) + inv_2zeta3C*(2*pc(4,5))

pc(1,6) = GB(2)*pc(1,3) + inv_2zeta3C*(preFactorOL)

pc(2,6) = GA(1)*pc(1,6)

pc(3,6) = GA(2)*pc(1,6) + inv_2zeta3C*(2*pc(1,3))

pc(4,6) = GA(3)*pc(1,6)

pc(5,6) = GA(1)*pc(2,6) + inv_2zeta3C*(pc(1,6))

pc(6,6) = GA(2)*pc(3,6) + inv_2zeta3C*(pc(1,6) + 2*pc(3,3))

pc(7,6) = GA(3)*pc(4,6) + inv_2zeta3C*(pc(1,6))

pc(8,6) = GA(2)*pc(2,6) + inv_2zeta3C*(2*pc(2,3))

pc(9,6) = GA(3)*pc(2,6)

pc(10,6) = GA(3)*pc(3,6)

pc(11,6) = GA(3)*pc(8,6)

pc(12,6) = GA(2)*pc(5,6) + inv_2zeta3C*(2*pc(5,3))

pc(13,6) = GA(3)*pc(5,6)

pc(14,6) = GA(1)*pc(6,6)

pc(15,6) = GA(3)*pc(6,6)

pc(16,6) = GA(1)*pc(7,6)

pc(17,6) = GA(2)*pc(7,6) + inv_2zeta3C*(2*pc(7,3))

pc(18,6) = GA(1)*pc(5,6) + inv_2zeta3C*(2*pc(2,6))

pc(19,6) = GA(2)*pc(6,6) + inv_2zeta3C*(2*pc(3,6) + 2*pc(6,3))

pc(20,6) = GA(3)*pc(7,6) + inv_2zeta3C*(2*pc(4,6))

pc(1,7) = GB(3)*pc(1,4) + inv_2zeta3C*(preFactorOL)

pc(2,7) = GA(1)*pc(1,7)

pc(3,7) = GA(2)*pc(1,7)

pc(4,7) = GA(3)*pc(1,7) + inv_2zeta3C*(2*pc(1,4))

pc(5,7) = GA(1)*pc(2,7) + inv_2zeta3C*(pc(1,7))

pc(6,7) = GA(2)*pc(3,7) + inv_2zeta3C*(pc(1,7))

pc(7,7) = GA(3)*pc(4,7) + inv_2zeta3C*(pc(1,7) + 2*pc(4,4))

pc(8,7) = GA(2)*pc(2,7)

pc(9,7) = GA(3)*pc(2,7) + inv_2zeta3C*(2*pc(2,4))

pc(10,7) = GA(3)*pc(3,7) + inv_2zeta3C*(2*pc(3,4))

pc(11,7) = GA(3)*pc(8,7) + inv_2zeta3C*(2*pc(8,4))

pc(12,7) = GA(2)*pc(5,7)

pc(13,7) = GA(3)*pc(5,7) + inv_2zeta3C*(2*pc(5,4))

pc(14,7) = GA(1)*pc(6,7)

pc(15,7) = GA(3)*pc(6,7) + inv_2zeta3C*(2*pc(6,4))

pc(16,7) = GA(1)*pc(7,7)

pc(17,7) = GA(2)*pc(7,7)

pc(18,7) = GA(1)*pc(5,7) + inv_2zeta3C*(2*pc(2,7))

pc(19,7) = GA(2)*pc(6,7) + inv_2zeta3C*(2*pc(3,7))

pc(20,7) = GA(3)*pc(7,7) + inv_2zeta3C*(2*pc(4,7) + 2*pc(7,4))

pc(1,8) = GB(2)*pc(1,2)

pc(2,8) = GA(1)*pc(1,8) + inv_2zeta3C*(pc(1,3))

pc(3,8) = GA(2)*pc(1,8) + inv_2zeta3C*(pc(1,2))

pc(4,8) = GA(3)*pc(1,8)

pc(5,8) = GA(1)*pc(2,8) + inv_2zeta3C*(pc(1,8) + pc(2,3))

pc(6,8) = GA(2)*pc(3,8) + inv_2zeta3C*(pc(1,8) + pc(3,2))

pc(7,8) = GA(3)*pc(4,8) + inv_2zeta3C*(pc(1,8))

pc(8,8) = GA(2)*pc(2,8) + inv_2zeta3C*(pc(2,2))

pc(9,8) = GA(3)*pc(2,8)

pc(10,8) = GA(3)*pc(3,8)

pc(11,8) = GA(3)*pc(8,8)

pc(12,8) = GA(2)*pc(5,8) + inv_2zeta3C*(pc(5,2))

pc(13,8) = GA(3)*pc(5,8)

pc(14,8) = GA(1)*pc(6,8) + inv_2zeta3C*(pc(6,3))

pc(15,8) = GA(3)*pc(6,8)

pc(16,8) = GA(1)*pc(7,8) + inv_2zeta3C*(pc(7,3))

pc(17,8) = GA(2)*pc(7,8) + inv_2zeta3C*(pc(7,2))

pc(18,8) = GA(1)*pc(5,8) + inv_2zeta3C*(2*pc(2,8) + pc(5,3))

pc(19,8) = GA(2)*pc(6,8) + inv_2zeta3C*(2*pc(3,8) + pc(6,2))

pc(20,8) = GA(3)*pc(7,8) + inv_2zeta3C*(2*pc(4,8))

pc(1,9) = GB(3)*pc(1,2)

pc(2,9) = GA(1)*pc(1,9) + inv_2zeta3C*(pc(1,4))

pc(3,9) = GA(2)*pc(1,9)

pc(4,9) = GA(3)*pc(1,9) + inv_2zeta3C*(pc(1,2))

pc(5,9) = GA(1)*pc(2,9) + inv_2zeta3C*(pc(1,9) + pc(2,4))

pc(6,9) = GA(2)*pc(3,9) + inv_2zeta3C*(pc(1,9))

pc(7,9) = GA(3)*pc(4,9) + inv_2zeta3C*(pc(1,9) + pc(4,2))

pc(8,9) = GA(2)*pc(2,9)

pc(9,9) = GA(3)*pc(2,9) + inv_2zeta3C*(pc(2,2))

pc(10,9) = GA(3)*pc(3,9) + inv_2zeta3C*(pc(3,2))

pc(11,9) = GA(3)*pc(8,9) + inv_2zeta3C*(pc(8,2))

pc(12,9) = GA(2)*pc(5,9)

pc(13,9) = GA(3)*pc(5,9) + inv_2zeta3C*(pc(5,2))

pc(14,9) = GA(1)*pc(6,9) + inv_2zeta3C*(pc(6,4))

pc(15,9) = GA(3)*pc(6,9) + inv_2zeta3C*(pc(6,2))

pc(16,9) = GA(1)*pc(7,9) + inv_2zeta3C*(pc(7,4))

pc(17,9) = GA(2)*pc(7,9)

pc(18,9) = GA(1)*pc(5,9) + inv_2zeta3C*(2*pc(2,9) + pc(5,4))

pc(19,9) = GA(2)*pc(6,9) + inv_2zeta3C*(2*pc(3,9))

pc(20,9) = GA(3)*pc(7,9) + inv_2zeta3C*(2*pc(4,9) + pc(7,2))

pc(1,10) = GB(3)*pc(1,3)

pc(2,10) = GA(1)*pc(1,10)

pc(3,10) = GA(2)*pc(1,10) + inv_2zeta3C*(pc(1,4))

pc(4,10) = GA(3)*pc(1,10) + inv_2zeta3C*(pc(1,3))

pc(5,10) = GA(1)*pc(2,10) + inv_2zeta3C*(pc(1,10))

pc(6,10) = GA(2)*pc(3,10) + inv_2zeta3C*(pc(1,10) + pc(3,4))

pc(7,10) = GA(3)*pc(4,10) + inv_2zeta3C*(pc(1,10) + pc(4,3))

pc(8,10) = GA(2)*pc(2,10) + inv_2zeta3C*(pc(2,4))

pc(9,10) = GA(3)*pc(2,10) + inv_2zeta3C*(pc(2,3))

pc(10,10) = GA(3)*pc(3,10) + inv_2zeta3C*(pc(3,3))

pc(11,10) = GA(3)*pc(8,10) + inv_2zeta3C*(pc(8,3))

pc(12,10) = GA(2)*pc(5,10) + inv_2zeta3C*(pc(5,4))

pc(13,10) = GA(3)*pc(5,10) + inv_2zeta3C*(pc(5,3))

pc(14,10) = GA(1)*pc(6,10)

pc(15,10) = GA(3)*pc(6,10) + inv_2zeta3C*(pc(6,3))

pc(16,10) = GA(1)*pc(7,10)

pc(17,10) = GA(2)*pc(7,10) + inv_2zeta3C*(pc(7,4))

pc(18,10) = GA(1)*pc(5,10) + inv_2zeta3C*(2*pc(2,10))

pc(19,10) = GA(2)*pc(6,10) + inv_2zeta3C*(2*pc(3,10) + pc(6,4))

pc(20,10) = GA(3)*pc(7,10) + inv_2zeta3C*(2*pc(4,10) + pc(7,3))

pc(1,11) = GB(3)*pc(1,8)

pc(2,11) = GA(1)*pc(1,11) + inv_2zeta3C*(pc(1,10))

pc(3,11) = GA(2)*pc(1,11) + inv_2zeta3C*(pc(1,9))

pc(4,11) = GA(3)*pc(1,11) + inv_2zeta3C*(pc(1,8))

pc(5,11) = GA(1)*pc(2,11) + inv_2zeta3C*(pc(1,11) + pc(2,10))

pc(6,11) = GA(2)*pc(3,11) + inv_2zeta3C*(pc(1,11) + pc(3,9))

pc(7,11) = GA(3)*pc(4,11) + inv_2zeta3C*(pc(1,11) + pc(4,8))

pc(8,11) = GA(2)*pc(2,11) + inv_2zeta3C*(pc(2,9))

pc(9,11) = GA(3)*pc(2,11) + inv_2zeta3C*(pc(2,8))

pc(10,11) = GA(3)*pc(3,11) + inv_2zeta3C*(pc(3,8))

pc(11,11) = GA(3)*pc(8,11) + inv_2zeta3C*(pc(8,8))

pc(12,11) = GA(2)*pc(5,11) + inv_2zeta3C*(pc(5,9))

pc(13,11) = GA(3)*pc(5,11) + inv_2zeta3C*(pc(5,8))

pc(14,11) = GA(1)*pc(6,11) + inv_2zeta3C*(pc(6,10))

pc(15,11) = GA(3)*pc(6,11) + inv_2zeta3C*(pc(6,8))

pc(16,11) = GA(1)*pc(7,11) + inv_2zeta3C*(pc(7,10))

pc(17,11) = GA(2)*pc(7,11) + inv_2zeta3C*(pc(7,9))

pc(18,11) = GA(1)*pc(5,11) + inv_2zeta3C*(2*pc(2,11) + pc(5,10))

pc(19,11) = GA(2)*pc(6,11) + inv_2zeta3C*(2*pc(3,11) + pc(6,9))

pc(20,11) = GA(3)*pc(7,11) + inv_2zeta3C*(2*pc(4,11) + pc(7,8))

pc(1,12) = GB(2)*pc(1,5)

pc(2,12) = GA(1)*pc(1,12) + inv_2zeta3C*(2*pc(1,8))

pc(3,12) = GA(2)*pc(1,12) + inv_2zeta3C*(pc(1,5))

pc(4,12) = GA(3)*pc(1,12)

pc(5,12) = GA(1)*pc(2,12) + inv_2zeta3C*(pc(1,12) + 2*pc(2,8))

pc(6,12) = GA(2)*pc(3,12) + inv_2zeta3C*(pc(1,12) + pc(3,5))

pc(7,12) = GA(3)*pc(4,12) + inv_2zeta3C*(pc(1,12))

pc(8,12) = GA(2)*pc(2,12) + inv_2zeta3C*(pc(2,5))

pc(9,12) = GA(3)*pc(2,12)

pc(10,12) = GA(3)*pc(3,12)

pc(11,12) = GA(3)*pc(8,12)

pc(12,12) = GA(2)*pc(5,12) + inv_2zeta3C*(pc(5,5))

pc(13,12) = GA(3)*pc(5,12)

pc(14,12) = GA(1)*pc(6,12) + inv_2zeta3C*(2*pc(6,8))

pc(15,12) = GA(3)*pc(6,12)

pc(16,12) = GA(1)*pc(7,12) + inv_2zeta3C*(2*pc(7,8))

pc(17,12) = GA(2)*pc(7,12) + inv_2zeta3C*(pc(7,5))

pc(18,12) = GA(1)*pc(5,12) + inv_2zeta3C*(2*pc(2,12) + 2*pc(5,8))

pc(19,12) = GA(2)*pc(6,12) + inv_2zeta3C*(2*pc(3,12) + pc(6,5))

pc(20,12) = GA(3)*pc(7,12) + inv_2zeta3C*(2*pc(4,12))

pc(1,13) = GB(3)*pc(1,5)

pc(2,13) = GA(1)*pc(1,13) + inv_2zeta3C*(2*pc(1,9))

pc(3,13) = GA(2)*pc(1,13)

pc(4,13) = GA(3)*pc(1,13) + inv_2zeta3C*(pc(1,5))

pc(5,13) = GA(1)*pc(2,13) + inv_2zeta3C*(pc(1,13) + 2*pc(2,9))

pc(6,13) = GA(2)*pc(3,13) + inv_2zeta3C*(pc(1,13))

pc(7,13) = GA(3)*pc(4,13) + inv_2zeta3C*(pc(1,13) + pc(4,5))

pc(8,13) = GA(2)*pc(2,13)

pc(9,13) = GA(3)*pc(2,13) + inv_2zeta3C*(pc(2,5))

pc(10,13) = GA(3)*pc(3,13) + inv_2zeta3C*(pc(3,5))

pc(11,13) = GA(3)*pc(8,13) + inv_2zeta3C*(pc(8,5))

pc(12,13) = GA(2)*pc(5,13)

pc(13,13) = GA(3)*pc(5,13) + inv_2zeta3C*(pc(5,5))

pc(14,13) = GA(1)*pc(6,13) + inv_2zeta3C*(2*pc(6,9))

pc(15,13) = GA(3)*pc(6,13) + inv_2zeta3C*(pc(6,5))

pc(16,13) = GA(1)*pc(7,13) + inv_2zeta3C*(2*pc(7,9))

pc(17,13) = GA(2)*pc(7,13)

pc(18,13) = GA(1)*pc(5,13) + inv_2zeta3C*(2*pc(2,13) + 2*pc(5,9))

pc(19,13) = GA(2)*pc(6,13) + inv_2zeta3C*(2*pc(3,13))

pc(20,13) = GA(3)*pc(7,13) + inv_2zeta3C*(2*pc(4,13) + pc(7,5))

pc(1,14) = GB(1)*pc(1,6)

pc(2,14) = GA(1)*pc(1,14) + inv_2zeta3C*(pc(1,6))

pc(3,14) = GA(2)*pc(1,14) + inv_2zeta3C*(2*pc(1,8))

pc(4,14) = GA(3)*pc(1,14)

pc(5,14) = GA(1)*pc(2,14) + inv_2zeta3C*(pc(1,14) + pc(2,6))

pc(6,14) = GA(2)*pc(3,14) + inv_2zeta3C*(pc(1,14) + 2*pc(3,8))

pc(7,14) = GA(3)*pc(4,14) + inv_2zeta3C*(pc(1,14))

pc(8,14) = GA(2)*pc(2,14) + inv_2zeta3C*(2*pc(2,8))

pc(9,14) = GA(3)*pc(2,14)

pc(10,14) = GA(3)*pc(3,14)

pc(11,14) = GA(3)*pc(8,14)

pc(12,14) = GA(2)*pc(5,14) + inv_2zeta3C*(2*pc(5,8))

pc(13,14) = GA(3)*pc(5,14)

pc(14,14) = GA(1)*pc(6,14) + inv_2zeta3C*(pc(6,6))

pc(15,14) = GA(3)*pc(6,14)

pc(16,14) = GA(1)*pc(7,14) + inv_2zeta3C*(pc(7,6))

pc(17,14) = GA(2)*pc(7,14) + inv_2zeta3C*(2*pc(7,8))

pc(18,14) = GA(1)*pc(5,14) + inv_2zeta3C*(2*pc(2,14) + pc(5,6))

pc(19,14) = GA(2)*pc(6,14) + inv_2zeta3C*(2*pc(3,14) + 2*pc(6,8))

pc(20,14) = GA(3)*pc(7,14) + inv_2zeta3C*(2*pc(4,14))

pc(1,15) = GB(3)*pc(1,6)

pc(2,15) = GA(1)*pc(1,15)

pc(3,15) = GA(2)*pc(1,15) + inv_2zeta3C*(2*pc(1,10))

pc(4,15) = GA(3)*pc(1,15) + inv_2zeta3C*(pc(1,6))

pc(5,15) = GA(1)*pc(2,15) + inv_2zeta3C*(pc(1,15))

pc(6,15) = GA(2)*pc(3,15) + inv_2zeta3C*(pc(1,15) + 2*pc(3,10))

pc(7,15) = GA(3)*pc(4,15) + inv_2zeta3C*(pc(1,15) + pc(4,6))

pc(8,15) = GA(2)*pc(2,15) + inv_2zeta3C*(2*pc(2,10))

pc(9,15) = GA(3)*pc(2,15) + inv_2zeta3C*(pc(2,6))

pc(10,15) = GA(3)*pc(3,15) + inv_2zeta3C*(pc(3,6))

pc(11,15) = GA(3)*pc(8,15) + inv_2zeta3C*(pc(8,6))

pc(12,15) = GA(2)*pc(5,15) + inv_2zeta3C*(2*pc(5,10))

pc(13,15) = GA(3)*pc(5,15) + inv_2zeta3C*(pc(5,6))

pc(14,15) = GA(1)*pc(6,15)

pc(15,15) = GA(3)*pc(6,15) + inv_2zeta3C*(pc(6,6))

pc(16,15) = GA(1)*pc(7,15)

pc(17,15) = GA(2)*pc(7,15) + inv_2zeta3C*(2*pc(7,10))

pc(18,15) = GA(1)*pc(5,15) + inv_2zeta3C*(2*pc(2,15))

pc(19,15) = GA(2)*pc(6,15) + inv_2zeta3C*(2*pc(3,15) + 2*pc(6,10))

pc(20,15) = GA(3)*pc(7,15) + inv_2zeta3C*(2*pc(4,15) + pc(7,6))

pc(1,16) = GB(1)*pc(1,7)

pc(2,16) = GA(1)*pc(1,16) + inv_2zeta3C*(pc(1,7))

pc(3,16) = GA(2)*pc(1,16)

pc(4,16) = GA(3)*pc(1,16) + inv_2zeta3C*(2*pc(1,9))

pc(5,16) = GA(1)*pc(2,16) + inv_2zeta3C*(pc(1,16) + pc(2,7))

pc(6,16) = GA(2)*pc(3,16) + inv_2zeta3C*(pc(1,16))

pc(7,16) = GA(3)*pc(4,16) + inv_2zeta3C*(pc(1,16) + 2*pc(4,9))

pc(8,16) = GA(2)*pc(2,16)

pc(9,16) = GA(3)*pc(2,16) + inv_2zeta3C*(2*pc(2,9))

pc(10,16) = GA(3)*pc(3,16) + inv_2zeta3C*(2*pc(3,9))

pc(11,16) = GA(3)*pc(8,16) + inv_2zeta3C*(2*pc(8,9))

pc(12,16) = GA(2)*pc(5,16)

pc(13,16) = GA(3)*pc(5,16) + inv_2zeta3C*(2*pc(5,9))

pc(14,16) = GA(1)*pc(6,16) + inv_2zeta3C*(pc(6,7))

pc(15,16) = GA(3)*pc(6,16) + inv_2zeta3C*(2*pc(6,9))

pc(16,16) = GA(1)*pc(7,16) + inv_2zeta3C*(pc(7,7))

pc(17,16) = GA(2)*pc(7,16)

pc(18,16) = GA(1)*pc(5,16) + inv_2zeta3C*(2*pc(2,16) + pc(5,7))

pc(19,16) = GA(2)*pc(6,16) + inv_2zeta3C*(2*pc(3,16))

pc(20,16) = GA(3)*pc(7,16) + inv_2zeta3C*(2*pc(4,16) + 2*pc(7,9))

pc(1,17) = GB(2)*pc(1,7)

pc(2,17) = GA(1)*pc(1,17)

pc(3,17) = GA(2)*pc(1,17) + inv_2zeta3C*(pc(1,7))

pc(4,17) = GA(3)*pc(1,17) + inv_2zeta3C*(2*pc(1,10))

pc(5,17) = GA(1)*pc(2,17) + inv_2zeta3C*(pc(1,17))

pc(6,17) = GA(2)*pc(3,17) + inv_2zeta3C*(pc(1,17) + pc(3,7))

pc(7,17) = GA(3)*pc(4,17) + inv_2zeta3C*(pc(1,17) + 2*pc(4,10))

pc(8,17) = GA(2)*pc(2,17) + inv_2zeta3C*(pc(2,7))

pc(9,17) = GA(3)*pc(2,17) + inv_2zeta3C*(2*pc(2,10))

pc(10,17) = GA(3)*pc(3,17) + inv_2zeta3C*(2*pc(3,10))

pc(11,17) = GA(3)*pc(8,17) + inv_2zeta3C*(2*pc(8,10))

pc(12,17) = GA(2)*pc(5,17) + inv_2zeta3C*(pc(5,7))

pc(13,17) = GA(3)*pc(5,17) + inv_2zeta3C*(2*pc(5,10))

pc(14,17) = GA(1)*pc(6,17)

pc(15,17) = GA(3)*pc(6,17) + inv_2zeta3C*(2*pc(6,10))

pc(16,17) = GA(1)*pc(7,17)

pc(17,17) = GA(2)*pc(7,17) + inv_2zeta3C*(pc(7,7))

pc(18,17) = GA(1)*pc(5,17) + inv_2zeta3C*(2*pc(2,17))

pc(19,17) = GA(2)*pc(6,17) + inv_2zeta3C*(2*pc(3,17) + pc(6,7))

pc(20,17) = GA(3)*pc(7,17) + inv_2zeta3C*(2*pc(4,17) + 2*pc(7,10))

pc(1,18) = GB(1)*pc(1,5) + inv_2zeta3C*(2*pc(1,2))

pc(2,18) = GA(1)*pc(1,18) + inv_2zeta3C*(3*pc(1,5))

pc(3,18) = GA(2)*pc(1,18)

pc(4,18) = GA(3)*pc(1,18)

pc(5,18) = GA(1)*pc(2,18) + inv_2zeta3C*(pc(1,18) + 3*pc(2,5))

pc(6,18) = GA(2)*pc(3,18) + inv_2zeta3C*(pc(1,18))

pc(7,18) = GA(3)*pc(4,18) + inv_2zeta3C*(pc(1,18))

pc(8,18) = GA(2)*pc(2,18)

pc(9,18) = GA(3)*pc(2,18)

pc(10,18) = GA(3)*pc(3,18)

pc(11,18) = GA(3)*pc(8,18)

pc(12,18) = GA(2)*pc(5,18)

pc(13,18) = GA(3)*pc(5,18)

pc(14,18) = GA(1)*pc(6,18) + inv_2zeta3C*(3*pc(6,5))

pc(15,18) = GA(3)*pc(6,18)

pc(16,18) = GA(1)*pc(7,18) + inv_2zeta3C*(3*pc(7,5))

pc(17,18) = GA(2)*pc(7,18)

pc(18,18) = GA(1)*pc(5,18) + inv_2zeta3C*(2*pc(2,18) + 3*pc(5,5))

pc(19,18) = GA(2)*pc(6,18) + inv_2zeta3C*(2*pc(3,18))

pc(20,18) = GA(3)*pc(7,18) + inv_2zeta3C*(2*pc(4,18))

pc(1,19) = GB(2)*pc(1,6) + inv_2zeta3C*(2*pc(1,3))

pc(2,19) = GA(1)*pc(1,19)

pc(3,19) = GA(2)*pc(1,19) + inv_2zeta3C*(3*pc(1,6))

pc(4,19) = GA(3)*pc(1,19)

pc(5,19) = GA(1)*pc(2,19) + inv_2zeta3C*(pc(1,19))

pc(6,19) = GA(2)*pc(3,19) + inv_2zeta3C*(pc(1,19) + 3*pc(3,6))

pc(7,19) = GA(3)*pc(4,19) + inv_2zeta3C*(pc(1,19))

pc(8,19) = GA(2)*pc(2,19) + inv_2zeta3C*(3*pc(2,6))

pc(9,19) = GA(3)*pc(2,19)

pc(10,19) = GA(3)*pc(3,19)

pc(11,19) = GA(3)*pc(8,19)

pc(12,19) = GA(2)*pc(5,19) + inv_2zeta3C*(3*pc(5,6))

pc(13,19) = GA(3)*pc(5,19)

pc(14,19) = GA(1)*pc(6,19)

pc(15,19) = GA(3)*pc(6,19)

pc(16,19) = GA(1)*pc(7,19)

pc(17,19) = GA(2)*pc(7,19) + inv_2zeta3C*(3*pc(7,6))

pc(18,19) = GA(1)*pc(5,19) + inv_2zeta3C*(2*pc(2,19))

pc(19,19) = GA(2)*pc(6,19) + inv_2zeta3C*(2*pc(3,19) + 3*pc(6,6))

pc(20,19) = GA(3)*pc(7,19) + inv_2zeta3C*(2*pc(4,19))

pc(1,20) = GB(3)*pc(1,7) + inv_2zeta3C*(2*pc(1,4))

pc(2,20) = GA(1)*pc(1,20)

pc(3,20) = GA(2)*pc(1,20)

pc(4,20) = GA(3)*pc(1,20) + inv_2zeta3C*(3*pc(1,7))

pc(5,20) = GA(1)*pc(2,20) + inv_2zeta3C*(pc(1,20))

pc(6,20) = GA(2)*pc(3,20) + inv_2zeta3C*(pc(1,20))

pc(7,20) = GA(3)*pc(4,20) + inv_2zeta3C*(pc(1,20) + 3*pc(4,7))

pc(8,20) = GA(2)*pc(2,20)

pc(9,20) = GA(3)*pc(2,20) + inv_2zeta3C*(3*pc(2,7))

pc(10,20) = GA(3)*pc(3,20) + inv_2zeta3C*(3*pc(3,7))

pc(11,20) = GA(3)*pc(8,20) + inv_2zeta3C*(3*pc(8,7))

pc(12,20) = GA(2)*pc(5,20)

pc(13,20) = GA(3)*pc(5,20) + inv_2zeta3C*(3*pc(5,7))

pc(14,20) = GA(1)*pc(6,20)

pc(15,20) = GA(3)*pc(6,20) + inv_2zeta3C*(3*pc(6,7))

pc(16,20) = GA(1)*pc(7,20)

pc(17,20) = GA(2)*pc(7,20)

pc(18,20) = GA(1)*pc(5,20) + inv_2zeta3C*(2*pc(2,20))

pc(19,20) = GA(2)*pc(6,20) + inv_2zeta3C*(2*pc(3,20))

pc(20,20) = GA(3)*pc(7,20) + inv_2zeta3C*(2*pc(4,20) + 3*pc(7,7))

sh(1,1) = pc(1,1)

sh(2,1) = pc(2,1)

sh(3,1) = pc(3,1)

sh(4,1) = pc(4,1)

sh(5,1) = pc(8,1)

sh(6,1) = pc(9,1)

sh(7,1) = pc(10,1)

sh(8,1) = pc(5,1) - pc(6,1)

sh(9,1) = 2*pc(7,1) - pc(5,1) - pc(6,1)

sh(10,1) = pc(11,1)

sh(11,1) = pc(13,1) - pc(15,1)

sh(12,1) = pc(18,1) - 3*pc(14,1)

sh(13,1) = 3*pc(12,1) - pc(19,1)

sh(14,1) = 2*pc(20,1) - 3*pc(13,1) - 3*pc(15,1)

sh(15,1) = 4*pc(16,1) - pc(18,1) - pc(14,1)

sh(16,1) = 4*pc(17,1) - pc(12,1) - pc(19,1)

sh(1,2) = pc(1,2)

sh(2,2) = pc(2,2)

sh(3,2) = pc(3,2)

sh(4,2) = pc(4,2)

sh(5,2) = pc(8,2)

sh(6,2) = pc(9,2)

sh(7,2) = pc(10,2)

sh(8,2) = pc(5,2) - pc(6,2)

sh(9,2) = 2*pc(7,2) - pc(5,2) - pc(6,2)

sh(10,2) = pc(11,2)

sh(11,2) = pc(13,2) - pc(15,2)

sh(12,2) = pc(18,2) - 3*pc(14,2)

sh(13,2) = 3*pc(12,2) - pc(19,2)

sh(14,2) = 2*pc(20,2) - 3*pc(13,2) - 3*pc(15,2)

sh(15,2) = 4*pc(16,2) - pc(18,2) - pc(14,2)

sh(16,2) = 4*pc(17,2) - pc(12,2) - pc(19,2)

sh(1,3) = pc(1,3)

sh(2,3) = pc(2,3)

sh(3,3) = pc(3,3)

sh(4,3) = pc(4,3)

sh(5,3) = pc(8,3)

sh(6,3) = pc(9,3)

sh(7,3) = pc(10,3)

sh(8,3) = pc(5,3) - pc(6,3)

sh(9,3) = 2*pc(7,3) - pc(5,3) - pc(6,3)

sh(10,3) = pc(11,3)

sh(11,3) = pc(13,3) - pc(15,3)

sh(12,3) = pc(18,3) - 3*pc(14,3)

sh(13,3) = 3*pc(12,3) - pc(19,3)

sh(14,3) = 2*pc(20,3) - 3*pc(13,3) - 3*pc(15,3)

sh(15,3) = 4*pc(16,3) - pc(18,3) - pc(14,3)

sh(16,3) = 4*pc(17,3) - pc(12,3) - pc(19,3)

sh(1,4) = pc(1,4)

sh(2,4) = pc(2,4)

sh(3,4) = pc(3,4)

sh(4,4) = pc(4,4)

sh(5,4) = pc(8,4)

sh(6,4) = pc(9,4)

sh(7,4) = pc(10,4)

sh(8,4) = pc(5,4) - pc(6,4)

sh(9,4) = 2*pc(7,4) - pc(5,4) - pc(6,4)

sh(10,4) = pc(11,4)

sh(11,4) = pc(13,4) - pc(15,4)

sh(12,4) = pc(18,4) - 3*pc(14,4)

sh(13,4) = 3*pc(12,4) - pc(19,4)

sh(14,4) = 2*pc(20,4) - 3*pc(13,4) - 3*pc(15,4)

sh(15,4) = 4*pc(16,4) - pc(18,4) - pc(14,4)

sh(16,4) = 4*pc(17,4) - pc(12,4) - pc(19,4)

sh(1,5) = pc(1,8)

sh(2,5) = pc(2,8)

sh(3,5) = pc(3,8)

sh(4,5) = pc(4,8)

sh(5,5) = pc(8,8)

sh(6,5) = pc(9,8)

sh(7,5) = pc(10,8)

sh(8,5) = pc(5,8) - pc(6,8)

sh(9,5) = 2*pc(7,8) - pc(5,8) - pc(6,8)

sh(10,5) = pc(11,8)

sh(11,5) = pc(13,8) - pc(15,8)

sh(12,5) = pc(18,8) - 3*pc(14,8)

sh(13,5) = 3*pc(12,8) - pc(19,8)

sh(14,5) = 2*pc(20,8) - 3*pc(13,8) - 3*pc(15,8)

sh(15,5) = 4*pc(16,8) - pc(18,8) - pc(14,8)

sh(16,5) = 4*pc(17,8) - pc(12,8) - pc(19,8)

sh(1,6) = pc(1,9)

sh(2,6) = pc(2,9)

sh(3,6) = pc(3,9)

sh(4,6) = pc(4,9)

sh(5,6) = pc(8,9)

sh(6,6) = pc(9,9)

sh(7,6) = pc(10,9)

sh(8,6) = pc(5,9) - pc(6,9)

sh(9,6) = 2*pc(7,9) - pc(5,9) - pc(6,9)

sh(10,6) = pc(11,9)

sh(11,6) = pc(13,9) - pc(15,9)

sh(12,6) = pc(18,9) - 3*pc(14,9)

sh(13,6) = 3*pc(12,9) - pc(19,9)

sh(14,6) = 2*pc(20,9) - 3*pc(13,9) - 3*pc(15,9)

sh(15,6) = 4*pc(16,9) - pc(18,9) - pc(14,9)

sh(16,6) = 4*pc(17,9) - pc(12,9) - pc(19,9)

sh(1,7) = pc(1,10)

sh(2,7) = pc(2,10)

sh(3,7) = pc(3,10)

sh(4,7) = pc(4,10)

sh(5,7) = pc(8,10)

sh(6,7) = pc(9,10)

sh(7,7) = pc(10,10)

sh(8,7) = pc(5,10) - pc(6,10)

sh(9,7) = 2*pc(7,10) - pc(5,10) - pc(6,10)

sh(10,7) = pc(11,10)

sh(11,7) = pc(13,10) - pc(15,10)

sh(12,7) = pc(18,10) - 3*pc(14,10)

sh(13,7) = 3*pc(12,10) - pc(19,10)

sh(14,7) = 2*pc(20,10) - 3*pc(13,10) - 3*pc(15,10)

sh(15,7) = 4*pc(16,10) - pc(18,10) - pc(14,10)

sh(16,7) = 4*pc(17,10) - pc(12,10) - pc(19,10)

sh(1,8) = pc(1,5) - pc(1,6)

sh(2,8) = pc(2,5) - pc(2,6)

sh(3,8) = pc(3,5) - pc(3,6)

sh(4,8) = pc(4,5) - pc(4,6)

sh(5,8) = pc(8,5) - pc(8,6)

sh(6,8) = pc(9,5) - pc(9,6)

sh(7,8) = pc(10,5) - pc(10,6)

sh(8,8) = pc(5,5) - pc(5,6) - pc(6,5) + pc(6,6)

sh(9,8) = 2*pc(7,5) - 2*pc(7,6) - pc(5,5) + pc(5,6) - pc(6,5) + pc(6,6)

sh(10,8) = pc(11,5) - pc(11,6)

sh(11,8) = pc(13,5) - pc(13,6) - pc(15,5) + pc(15,6)

sh(12,8) = pc(18,5) - pc(18,6) - 3*pc(14,5) + 3*pc(14,6)

sh(13,8) = 3*pc(12,5) - 3*pc(12,6) - pc(19,5) + pc(19,6)

sh(14,8) = 2*pc(20,5) - 2*pc(20,6) - 3*pc(13,5) + 3*pc(13,6) - 3*pc(15,5) + 3*&
&pc(15,6)

sh(15,8) = 4*pc(16,5) - 4*pc(16,6) - pc(18,5) + pc(18,6) - pc(14,5) + pc(14,6)

sh(16,8) = 4*pc(17,5) - 4*pc(17,6) - pc(12,5) + pc(12,6) - pc(19,5) + pc(19,6)

sh(1,9) = 2*pc(1,7) - pc(1,5) - pc(1,6)

sh(2,9) = 2*pc(2,7) - pc(2,5) - pc(2,6)

sh(3,9) = 2*pc(3,7) - pc(3,5) - pc(3,6)

sh(4,9) = 2*pc(4,7) - pc(4,5) - pc(4,6)

sh(5,9) = 2*pc(8,7) - pc(8,5) - pc(8,6)

sh(6,9) = 2*pc(9,7) - pc(9,5) - pc(9,6)

sh(7,9) = 2*pc(10,7) - pc(10,5) - pc(10,6)

sh(8,9) = 2*pc(5,7) - pc(5,5) - pc(5,6) - 2*pc(6,7) + pc(6,5) + pc(6,6)

sh(9,9) = 4*pc(7,7) - 2*pc(7,5) - 2*pc(7,6) - 2*pc(5,7) + pc(5,5) + pc(5,6) - &
&2*pc(6,7) + pc(6,5) + pc(6,6)

sh(10,9) = 2*pc(11,7) - pc(11,5) - pc(11,6)

sh(11,9) = 2*pc(13,7) - pc(13,5) - pc(13,6) - 2*pc(15,7) + pc(15,5) + pc(15,6)

sh(12,9) = 2*pc(18,7) - pc(18,5) - pc(18,6) - 6*pc(14,7) + 3*pc(14,5) + 3*pc(1&
&4,6)

sh(13,9) = 6*pc(12,7) - 3*pc(12,5) - 3*pc(12,6) - 2*pc(19,7) + pc(19,5) + pc(1&
&9,6)

sh(14,9) = 4*pc(20,7) - 2*pc(20,5) - 2*pc(20,6) - 6*pc(13,7) + 3*pc(13,5) + 3*&
&pc(13,6) - 6*pc(15,7) + 3*pc(15,5) + 3*pc(15,6)

sh(15,9) = 8*pc(16,7) - 4*pc(16,5) - 4*pc(16,6) - 2*pc(18,7) + pc(18,5) + pc(1&
&8,6) - 2*pc(14,7) + pc(14,5) + pc(14,6)

sh(16,9) = 8*pc(17,7) - 4*pc(17,5) - 4*pc(17,6) - 2*pc(12,7) + pc(12,5) + pc(1&
&2,6) - 2*pc(19,7) + pc(19,5) + pc(19,6)

sh(1,10) = pc(1,11)

sh(2,10) = pc(2,11)

sh(3,10) = pc(3,11)

sh(4,10) = pc(4,11)

sh(5,10) = pc(8,11)

sh(6,10) = pc(9,11)

sh(7,10) = pc(10,11)

sh(8,10) = pc(5,11) - pc(6,11)

sh(9,10) = 2*pc(7,11) - pc(5,11) - pc(6,11)

sh(10,10) = pc(11,11)

sh(11,10) = pc(13,11) - pc(15,11)

sh(12,10) = pc(18,11) - 3*pc(14,11)

sh(13,10) = 3*pc(12,11) - pc(19,11)

sh(14,10) = 2*pc(20,11) - 3*pc(13,11) - 3*pc(15,11)

sh(15,10) = 4*pc(16,11) - pc(18,11) - pc(14,11)

sh(16,10) = 4*pc(17,11) - pc(12,11) - pc(19,11)

sh(1,11) = pc(1,13) - pc(1,15)

sh(2,11) = pc(2,13) - pc(2,15)

sh(3,11) = pc(3,13) - pc(3,15)

sh(4,11) = pc(4,13) - pc(4,15)

sh(5,11) = pc(8,13) - pc(8,15)

sh(6,11) = pc(9,13) - pc(9,15)

sh(7,11) = pc(10,13) - pc(10,15)

sh(8,11) = pc(5,13) - pc(5,15) - pc(6,13) + pc(6,15)

sh(9,11) = 2*pc(7,13) - 2*pc(7,15) - pc(5,13) + pc(5,15) - pc(6,13) + pc(6,15)

sh(10,11) = pc(11,13) - pc(11,15)

sh(11,11) = pc(13,13) - pc(13,15) - pc(15,13) + pc(15,15)

sh(12,11) = pc(18,13) - pc(18,15) - 3*pc(14,13) + 3*pc(14,15)

sh(13,11) = 3*pc(12,13) - 3*pc(12,15) - pc(19,13) + pc(19,15)

sh(14,11) = 2*pc(20,13) - 2*pc(20,15) - 3*pc(13,13) + 3*pc(13,15) - 3*pc(15,13&
&) + 3*pc(15,15)

sh(15,11) = 4*pc(16,13) - 4*pc(16,15) - pc(18,13) + pc(18,15) - pc(14,13) + pc&
&(14,15)

sh(16,11) = 4*pc(17,13) - 4*pc(17,15) - pc(12,13) + pc(12,15) - pc(19,13) + pc&
&(19,15)

sh(1,12) = pc(1,18) - 3*pc(1,14)

sh(2,12) = pc(2,18) - 3*pc(2,14)

sh(3,12) = pc(3,18) - 3*pc(3,14)

sh(4,12) = pc(4,18) - 3*pc(4,14)

sh(5,12) = pc(8,18) - 3*pc(8,14)

sh(6,12) = pc(9,18) - 3*pc(9,14)

sh(7,12) = pc(10,18) - 3*pc(10,14)

sh(8,12) = pc(5,18) - 3*pc(5,14) - pc(6,18) + 3*pc(6,14)

sh(9,12) = 2*pc(7,18) - 6*pc(7,14) - pc(5,18) + 3*pc(5,14) - pc(6,18) + 3*pc(6&
&,14)

sh(10,12) = pc(11,18) - 3*pc(11,14)

sh(11,12) = pc(13,18) - 3*pc(13,14) - pc(15,18) + 3*pc(15,14)

sh(12,12) = pc(18,18) - 3*pc(18,14) - 3*pc(14,18) + 9*pc(14,14)

sh(13,12) = 3*pc(12,18) - 9*pc(12,14) - pc(19,18) + 3*pc(19,14)

sh(14,12) = 2*pc(20,18) - 6*pc(20,14) - 3*pc(13,18) + 9*pc(13,14) - 3*pc(15,18&
&) + 9*pc(15,14)

sh(15,12) = 4*pc(16,18) - 12*pc(16,14) - pc(18,18) + 3*pc(18,14) - pc(14,18) +&
& 3*pc(14,14)

sh(16,12) = 4*pc(17,18) - 12*pc(17,14) - pc(12,18) + 3*pc(12,14) - pc(19,18) +&
& 3*pc(19,14)

sh(1,13) = 3*pc(1,12) - pc(1,19)

sh(2,13) = 3*pc(2,12) - pc(2,19)

sh(3,13) = 3*pc(3,12) - pc(3,19)

sh(4,13) = 3*pc(4,12) - pc(4,19)

sh(5,13) = 3*pc(8,12) - pc(8,19)

sh(6,13) = 3*pc(9,12) - pc(9,19)

sh(7,13) = 3*pc(10,12) - pc(10,19)

sh(8,13) = 3*pc(5,12) - pc(5,19) - 3*pc(6,12) + pc(6,19)

sh(9,13) = 6*pc(7,12) - 2*pc(7,19) - 3*pc(5,12) + pc(5,19) - 3*pc(6,12) + pc(6&
&,19)

sh(10,13) = 3*pc(11,12) - pc(11,19)

sh(11,13) = 3*pc(13,12) - pc(13,19) - 3*pc(15,12) + pc(15,19)

sh(12,13) = 3*pc(18,12) - pc(18,19) - 9*pc(14,12) + 3*pc(14,19)

sh(13,13) = 9*pc(12,12) - 3*pc(12,19) - 3*pc(19,12) + pc(19,19)

sh(14,13) = 6*pc(20,12) - 2*pc(20,19) - 9*pc(13,12) + 3*pc(13,19) - 9*pc(15,12&
&) + 3*pc(15,19)

sh(15,13) = 12*pc(16,12) - 4*pc(16,19) - 3*pc(18,12) + pc(18,19) - 3*pc(14,12)&
& + pc(14,19)

sh(16,13) = 12*pc(17,12) - 4*pc(17,19) - 3*pc(12,12) + pc(12,19) - 3*pc(19,12)&
& + pc(19,19)

sh(1,14) = 2*pc(1,20) - 3*pc(1,13) - 3*pc(1,15)

sh(2,14) = 2*pc(2,20) - 3*pc(2,13) - 3*pc(2,15)

sh(3,14) = 2*pc(3,20) - 3*pc(3,13) - 3*pc(3,15)

sh(4,14) = 2*pc(4,20) - 3*pc(4,13) - 3*pc(4,15)

sh(5,14) = 2*pc(8,20) - 3*pc(8,13) - 3*pc(8,15)

sh(6,14) = 2*pc(9,20) - 3*pc(9,13) - 3*pc(9,15)

sh(7,14) = 2*pc(10,20) - 3*pc(10,13) - 3*pc(10,15)

sh(8,14) = 2*pc(5,20) - 3*pc(5,13) - 3*pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*&
&pc(6,15)

sh(9,14) = 4*pc(7,20) - 6*pc(7,13) - 6*pc(7,15) - 2*pc(5,20) + 3*pc(5,13) + 3*&
&pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*pc(6,15)

sh(10,14) = 2*pc(11,20) - 3*pc(11,13) - 3*pc(11,15)

sh(11,14) = 2*pc(13,20) - 3*pc(13,13) - 3*pc(13,15) - 2*pc(15,20) + 3*pc(15,13&
&) + 3*pc(15,15)

sh(12,14) = 2*pc(18,20) - 3*pc(18,13) - 3*pc(18,15) - 6*pc(14,20) + 9*pc(14,13&
&) + 9*pc(14,15)

sh(13,14) = 6*pc(12,20) - 9*pc(12,13) - 9*pc(12,15) - 2*pc(19,20) + 3*pc(19,13&
&) + 3*pc(19,15)

sh(14,14) = 4*pc(20,20) - 6*pc(20,13) - 6*pc(20,15) - 6*pc(13,20) + 9*pc(13,13&
&) + 9*pc(13,15) - 6*pc(15,20) + 9*pc(15,13) + 9*pc(15,15)

sh(15,14) = 8*pc(16,20) - 12*pc(16,13) - 12*pc(16,15) - 2*pc(18,20) + 3*pc(18,&
&13) + 3*pc(18,15) - 2*pc(14,20) + 3*pc(14,13) + 3*pc(14,15)

sh(16,14) = 8*pc(17,20) - 12*pc(17,13) - 12*pc(17,15) - 2*pc(12,20) + 3*pc(12,&
&13) + 3*pc(12,15) - 2*pc(19,20) + 3*pc(19,13) + 3*pc(19,15)

sh(1,15) = 4*pc(1,16) - pc(1,18) - pc(1,14)

sh(2,15) = 4*pc(2,16) - pc(2,18) - pc(2,14)

sh(3,15) = 4*pc(3,16) - pc(3,18) - pc(3,14)

sh(4,15) = 4*pc(4,16) - pc(4,18) - pc(4,14)

sh(5,15) = 4*pc(8,16) - pc(8,18) - pc(8,14)

sh(6,15) = 4*pc(9,16) - pc(9,18) - pc(9,14)

sh(7,15) = 4*pc(10,16) - pc(10,18) - pc(10,14)

sh(8,15) = 4*pc(5,16) - pc(5,18) - pc(5,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(9,15) = 8*pc(7,16) - 2*pc(7,18) - 2*pc(7,14) - 4*pc(5,16) + pc(5,18) + pc(5&
&,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(10,15) = 4*pc(11,16) - pc(11,18) - pc(11,14)

sh(11,15) = 4*pc(13,16) - pc(13,18) - pc(13,14) - 4*pc(15,16) + pc(15,18) + pc&
&(15,14)

sh(12,15) = 4*pc(18,16) - pc(18,18) - pc(18,14) - 12*pc(14,16) + 3*pc(14,18) +&
& 3*pc(14,14)

sh(13,15) = 12*pc(12,16) - 3*pc(12,18) - 3*pc(12,14) - 4*pc(19,16) + pc(19,18)&
& + pc(19,14)

sh(14,15) = 8*pc(20,16) - 2*pc(20,18) - 2*pc(20,14) - 12*pc(13,16) + 3*pc(13,1&
&8) + 3*pc(13,14) - 12*pc(15,16) + 3*pc(15,18) + 3*pc(15,14)

sh(15,15) = 16*pc(16,16) - 4*pc(16,18) - 4*pc(16,14) - 4*pc(18,16) + pc(18,18)&
& + pc(18,14) - 4*pc(14,16) + pc(14,18) + pc(14,14)

sh(16,15) = 16*pc(17,16) - 4*pc(17,18) - 4*pc(17,14) - 4*pc(12,16) + pc(12,18)&
& + pc(12,14) - 4*pc(19,16) + pc(19,18) + pc(19,14)

sh(1,16) = 4*pc(1,17) - pc(1,12) - pc(1,19)

sh(2,16) = 4*pc(2,17) - pc(2,12) - pc(2,19)

sh(3,16) = 4*pc(3,17) - pc(3,12) - pc(3,19)

sh(4,16) = 4*pc(4,17) - pc(4,12) - pc(4,19)

sh(5,16) = 4*pc(8,17) - pc(8,12) - pc(8,19)

sh(6,16) = 4*pc(9,17) - pc(9,12) - pc(9,19)

sh(7,16) = 4*pc(10,17) - pc(10,12) - pc(10,19)

sh(8,16) = 4*pc(5,17) - pc(5,12) - pc(5,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(9,16) = 8*pc(7,17) - 2*pc(7,12) - 2*pc(7,19) - 4*pc(5,17) + pc(5,12) + pc(5&
&,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(10,16) = 4*pc(11,17) - pc(11,12) - pc(11,19)

sh(11,16) = 4*pc(13,17) - pc(13,12) - pc(13,19) - 4*pc(15,17) + pc(15,12) + pc&
&(15,19)

sh(12,16) = 4*pc(18,17) - pc(18,12) - pc(18,19) - 12*pc(14,17) + 3*pc(14,12) +&
& 3*pc(14,19)

sh(13,16) = 12*pc(12,17) - 3*pc(12,12) - 3*pc(12,19) - 4*pc(19,17) + pc(19,12)&
& + pc(19,19)

sh(14,16) = 8*pc(20,17) - 2*pc(20,12) - 2*pc(20,19) - 12*pc(13,17) + 3*pc(13,1&
&2) + 3*pc(13,19) - 12*pc(15,17) + 3*pc(15,12) + 3*pc(15,19)

sh(15,16) = 16*pc(16,17) - 4*pc(16,12) - 4*pc(16,19) - 4*pc(18,17) + pc(18,12)&
& + pc(18,19) - 4*pc(14,17) + pc(14,12) + pc(14,19)

sh(16,16) = 16*pc(17,17) - 4*pc(17,12) - 4*pc(17,19) - 4*pc(12,17) + pc(12,12)&
& + pc(12,19) - 4*pc(19,17) + pc(19,12) + pc(19,19)


   end subroutine electron3CIntgAna

   subroutine electron3CIntgNum(a1,a2,a3,A,B,C,pc,sh,cell_size,step_size)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2, a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   real (kind=double), dimension (20,20), intent(out) :: pc
   real (kind=double), dimension (16,16), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i
   integer :: num_steps
   integer, dimension (20,3) :: triads
   integer, dimension (16,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos
   real (kind=double), dimension (3) :: xyz, xyz_sum, xyz_soln

   ! Initialize local variables.
   triads(1,:) = (/0,0,0/)
   triads(2,:) = (/1,0,0/)
   triads(3,:) = (/0,1,0/)
   triads(4,:) = (/0,0,1/)
   triads(5,:) = (/2,0,0/)
   triads(6,:) = (/0,2,0/)
   triads(7,:) = (/0,0,2/)
   triads(8,:) = (/1,1,0/)
   triads(9,:) = (/1,0,1/)
   triads(10,:) = (/0,1,1/)
   triads(11,:) = (/1,1,1/)
   triads(12,:) = (/2,1,0/)
   triads(13,:) = (/2,0,1/)
   triads(14,:) = (/1,2,0/)
   triads(15,:) = (/0,2,1/)
   triads(16,:) = (/1,0,2/)
   triads(17,:) = (/0,1,2/)
   triads(18,:) = (/3,0,0/)
   triads(19,:) = (/0,3,0/)
   triads(20,:) = (/0,0,3/)
   conversion(1,1,:) = (/1,0,0/)
   conversion(1,2,:) = (/1,1,1/)
   conversion(2,1,:) = (/1,0,0/)
   conversion(2,2,:) = (/2,1,1/)
   conversion(3,1,:) = (/1,0,0/)
   conversion(3,2,:) = (/3,1,1/)
   conversion(4,1,:) = (/1,0,0/)
   conversion(4,2,:) = (/4,1,1/)
   conversion(5,1,:) = (/1,0,0/)
   conversion(5,2,:) = (/8,1,1/)
   conversion(6,1,:) = (/1,0,0/)
   conversion(6,2,:) = (/9,1,1/)
   conversion(7,1,:) = (/1,0,0/)
   conversion(7,2,:) = (/10,1,1/)
   conversion(8,1,:) = (/1,-1,0/)
   conversion(8,2,:) = (/5,6,1/)
   conversion(9,1,:) = (/2,-1,-1/)
   conversion(9,2,:) = (/7,5,6/)
   conversion(10,1,:) = (/1,0,0/)
   conversion(10,2,:) = (/11,1,1/)
   conversion(11,1,:) = (/1,-1,0/)
   conversion(11,2,:) = (/13,15,1/)
   conversion(12,1,:) = (/1,-3,0/)
   conversion(12,2,:) = (/18,14,1/)
   conversion(13,1,:) = (/3,-1,0/)
   conversion(13,2,:) = (/12,19,1/)
   conversion(14,1,:) = (/2,-3,-3/)
   conversion(14,2,:) = (/20,13,15/)
   conversion(15,1,:) = (/4,-1,-1/)
   conversion(15,2,:) = (/16,18,14/)
   conversion(16,1,:) = (/4,-1,-1/)
   conversion(16,2,:) = (/17,12,19/)

   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   do p = 1, 20
      do q = 1, 20

         ! Assign l1 and l2 values for each gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize sum variables.
         xyz_sum(:) = 0.0d0

         do i = 0, num_steps
            xyz(:) = (start_pos + (i*step_size))
            xyz_soln(:) = step_size * &
&(xyz(:)-A(:))**l1(:)*(xyz(:)-B(:))**l2(:)*exp(-a1*(xyz(:)-A(:))**2)*exp(-a2*(x&
&yz(:)-B(:))**2)*exp(-a3*(xyz(:)-C(:))**2)
            xyz_sum(:) = xyz_sum(:) + xyz_soln(:)
         enddo

         pc(q,p) = product(xyz_sum(:))
      enddo
    enddo

sh(1,1) = pc(1,1)

sh(2,1) = pc(2,1)

sh(3,1) = pc(3,1)

sh(4,1) = pc(4,1)

sh(5,1) = pc(8,1)

sh(6,1) = pc(9,1)

sh(7,1) = pc(10,1)

sh(8,1) = pc(5,1) - pc(6,1)

sh(9,1) = 2*pc(7,1) - pc(5,1) - pc(6,1)

sh(10,1) = pc(11,1)

sh(11,1) = pc(13,1) - pc(15,1)

sh(12,1) = pc(18,1) - 3*pc(14,1)

sh(13,1) = 3*pc(12,1) - pc(19,1)

sh(14,1) = 2*pc(20,1) - 3*pc(13,1) - 3*pc(15,1)

sh(15,1) = 4*pc(16,1) - pc(18,1) - pc(14,1)

sh(16,1) = 4*pc(17,1) - pc(12,1) - pc(19,1)

sh(1,2) = pc(1,2)

sh(2,2) = pc(2,2)

sh(3,2) = pc(3,2)

sh(4,2) = pc(4,2)

sh(5,2) = pc(8,2)

sh(6,2) = pc(9,2)

sh(7,2) = pc(10,2)

sh(8,2) = pc(5,2) - pc(6,2)

sh(9,2) = 2*pc(7,2) - pc(5,2) - pc(6,2)

sh(10,2) = pc(11,2)

sh(11,2) = pc(13,2) - pc(15,2)

sh(12,2) = pc(18,2) - 3*pc(14,2)

sh(13,2) = 3*pc(12,2) - pc(19,2)

sh(14,2) = 2*pc(20,2) - 3*pc(13,2) - 3*pc(15,2)

sh(15,2) = 4*pc(16,2) - pc(18,2) - pc(14,2)

sh(16,2) = 4*pc(17,2) - pc(12,2) - pc(19,2)

sh(1,3) = pc(1,3)

sh(2,3) = pc(2,3)

sh(3,3) = pc(3,3)

sh(4,3) = pc(4,3)

sh(5,3) = pc(8,3)

sh(6,3) = pc(9,3)

sh(7,3) = pc(10,3)

sh(8,3) = pc(5,3) - pc(6,3)

sh(9,3) = 2*pc(7,3) - pc(5,3) - pc(6,3)

sh(10,3) = pc(11,3)

sh(11,3) = pc(13,3) - pc(15,3)

sh(12,3) = pc(18,3) - 3*pc(14,3)

sh(13,3) = 3*pc(12,3) - pc(19,3)

sh(14,3) = 2*pc(20,3) - 3*pc(13,3) - 3*pc(15,3)

sh(15,3) = 4*pc(16,3) - pc(18,3) - pc(14,3)

sh(16,3) = 4*pc(17,3) - pc(12,3) - pc(19,3)

sh(1,4) = pc(1,4)

sh(2,4) = pc(2,4)

sh(3,4) = pc(3,4)

sh(4,4) = pc(4,4)

sh(5,4) = pc(8,4)

sh(6,4) = pc(9,4)

sh(7,4) = pc(10,4)

sh(8,4) = pc(5,4) - pc(6,4)

sh(9,4) = 2*pc(7,4) - pc(5,4) - pc(6,4)

sh(10,4) = pc(11,4)

sh(11,4) = pc(13,4) - pc(15,4)

sh(12,4) = pc(18,4) - 3*pc(14,4)

sh(13,4) = 3*pc(12,4) - pc(19,4)

sh(14,4) = 2*pc(20,4) - 3*pc(13,4) - 3*pc(15,4)

sh(15,4) = 4*pc(16,4) - pc(18,4) - pc(14,4)

sh(16,4) = 4*pc(17,4) - pc(12,4) - pc(19,4)

sh(1,5) = pc(1,8)

sh(2,5) = pc(2,8)

sh(3,5) = pc(3,8)

sh(4,5) = pc(4,8)

sh(5,5) = pc(8,8)

sh(6,5) = pc(9,8)

sh(7,5) = pc(10,8)

sh(8,5) = pc(5,8) - pc(6,8)

sh(9,5) = 2*pc(7,8) - pc(5,8) - pc(6,8)

sh(10,5) = pc(11,8)

sh(11,5) = pc(13,8) - pc(15,8)

sh(12,5) = pc(18,8) - 3*pc(14,8)

sh(13,5) = 3*pc(12,8) - pc(19,8)

sh(14,5) = 2*pc(20,8) - 3*pc(13,8) - 3*pc(15,8)

sh(15,5) = 4*pc(16,8) - pc(18,8) - pc(14,8)

sh(16,5) = 4*pc(17,8) - pc(12,8) - pc(19,8)

sh(1,6) = pc(1,9)

sh(2,6) = pc(2,9)

sh(3,6) = pc(3,9)

sh(4,6) = pc(4,9)

sh(5,6) = pc(8,9)

sh(6,6) = pc(9,9)

sh(7,6) = pc(10,9)

sh(8,6) = pc(5,9) - pc(6,9)

sh(9,6) = 2*pc(7,9) - pc(5,9) - pc(6,9)

sh(10,6) = pc(11,9)

sh(11,6) = pc(13,9) - pc(15,9)

sh(12,6) = pc(18,9) - 3*pc(14,9)

sh(13,6) = 3*pc(12,9) - pc(19,9)

sh(14,6) = 2*pc(20,9) - 3*pc(13,9) - 3*pc(15,9)

sh(15,6) = 4*pc(16,9) - pc(18,9) - pc(14,9)

sh(16,6) = 4*pc(17,9) - pc(12,9) - pc(19,9)

sh(1,7) = pc(1,10)

sh(2,7) = pc(2,10)

sh(3,7) = pc(3,10)

sh(4,7) = pc(4,10)

sh(5,7) = pc(8,10)

sh(6,7) = pc(9,10)

sh(7,7) = pc(10,10)

sh(8,7) = pc(5,10) - pc(6,10)

sh(9,7) = 2*pc(7,10) - pc(5,10) - pc(6,10)

sh(10,7) = pc(11,10)

sh(11,7) = pc(13,10) - pc(15,10)

sh(12,7) = pc(18,10) - 3*pc(14,10)

sh(13,7) = 3*pc(12,10) - pc(19,10)

sh(14,7) = 2*pc(20,10) - 3*pc(13,10) - 3*pc(15,10)

sh(15,7) = 4*pc(16,10) - pc(18,10) - pc(14,10)

sh(16,7) = 4*pc(17,10) - pc(12,10) - pc(19,10)

sh(1,8) = pc(1,5) - pc(1,6)

sh(2,8) = pc(2,5) - pc(2,6)

sh(3,8) = pc(3,5) - pc(3,6)

sh(4,8) = pc(4,5) - pc(4,6)

sh(5,8) = pc(8,5) - pc(8,6)

sh(6,8) = pc(9,5) - pc(9,6)

sh(7,8) = pc(10,5) - pc(10,6)

sh(8,8) = pc(5,5) - pc(5,6) - pc(6,5) + pc(6,6)

sh(9,8) = 2*pc(7,5) - 2*pc(7,6) - pc(5,5) + pc(5,6) - pc(6,5) + pc(6,6)

sh(10,8) = pc(11,5) - pc(11,6)

sh(11,8) = pc(13,5) - pc(13,6) - pc(15,5) + pc(15,6)

sh(12,8) = pc(18,5) - pc(18,6) - 3*pc(14,5) + 3*pc(14,6)

sh(13,8) = 3*pc(12,5) - 3*pc(12,6) - pc(19,5) + pc(19,6)

sh(14,8) = 2*pc(20,5) - 2*pc(20,6) - 3*pc(13,5) + 3*pc(13,6) - 3*pc(15,5) + 3*&
&pc(15,6)

sh(15,8) = 4*pc(16,5) - 4*pc(16,6) - pc(18,5) + pc(18,6) - pc(14,5) + pc(14,6)

sh(16,8) = 4*pc(17,5) - 4*pc(17,6) - pc(12,5) + pc(12,6) - pc(19,5) + pc(19,6)

sh(1,9) = 2*pc(1,7) - pc(1,5) - pc(1,6)

sh(2,9) = 2*pc(2,7) - pc(2,5) - pc(2,6)

sh(3,9) = 2*pc(3,7) - pc(3,5) - pc(3,6)

sh(4,9) = 2*pc(4,7) - pc(4,5) - pc(4,6)

sh(5,9) = 2*pc(8,7) - pc(8,5) - pc(8,6)

sh(6,9) = 2*pc(9,7) - pc(9,5) - pc(9,6)

sh(7,9) = 2*pc(10,7) - pc(10,5) - pc(10,6)

sh(8,9) = 2*pc(5,7) - pc(5,5) - pc(5,6) - 2*pc(6,7) + pc(6,5) + pc(6,6)

sh(9,9) = 4*pc(7,7) - 2*pc(7,5) - 2*pc(7,6) - 2*pc(5,7) + pc(5,5) + pc(5,6) - &
&2*pc(6,7) + pc(6,5) + pc(6,6)

sh(10,9) = 2*pc(11,7) - pc(11,5) - pc(11,6)

sh(11,9) = 2*pc(13,7) - pc(13,5) - pc(13,6) - 2*pc(15,7) + pc(15,5) + pc(15,6)

sh(12,9) = 2*pc(18,7) - pc(18,5) - pc(18,6) - 6*pc(14,7) + 3*pc(14,5) + 3*pc(1&
&4,6)

sh(13,9) = 6*pc(12,7) - 3*pc(12,5) - 3*pc(12,6) - 2*pc(19,7) + pc(19,5) + pc(1&
&9,6)

sh(14,9) = 4*pc(20,7) - 2*pc(20,5) - 2*pc(20,6) - 6*pc(13,7) + 3*pc(13,5) + 3*&
&pc(13,6) - 6*pc(15,7) + 3*pc(15,5) + 3*pc(15,6)

sh(15,9) = 8*pc(16,7) - 4*pc(16,5) - 4*pc(16,6) - 2*pc(18,7) + pc(18,5) + pc(1&
&8,6) - 2*pc(14,7) + pc(14,5) + pc(14,6)

sh(16,9) = 8*pc(17,7) - 4*pc(17,5) - 4*pc(17,6) - 2*pc(12,7) + pc(12,5) + pc(1&
&2,6) - 2*pc(19,7) + pc(19,5) + pc(19,6)

sh(1,10) = pc(1,11)

sh(2,10) = pc(2,11)

sh(3,10) = pc(3,11)

sh(4,10) = pc(4,11)

sh(5,10) = pc(8,11)

sh(6,10) = pc(9,11)

sh(7,10) = pc(10,11)

sh(8,10) = pc(5,11) - pc(6,11)

sh(9,10) = 2*pc(7,11) - pc(5,11) - pc(6,11)

sh(10,10) = pc(11,11)

sh(11,10) = pc(13,11) - pc(15,11)

sh(12,10) = pc(18,11) - 3*pc(14,11)

sh(13,10) = 3*pc(12,11) - pc(19,11)

sh(14,10) = 2*pc(20,11) - 3*pc(13,11) - 3*pc(15,11)

sh(15,10) = 4*pc(16,11) - pc(18,11) - pc(14,11)

sh(16,10) = 4*pc(17,11) - pc(12,11) - pc(19,11)

sh(1,11) = pc(1,13) - pc(1,15)

sh(2,11) = pc(2,13) - pc(2,15)

sh(3,11) = pc(3,13) - pc(3,15)

sh(4,11) = pc(4,13) - pc(4,15)

sh(5,11) = pc(8,13) - pc(8,15)

sh(6,11) = pc(9,13) - pc(9,15)

sh(7,11) = pc(10,13) - pc(10,15)

sh(8,11) = pc(5,13) - pc(5,15) - pc(6,13) + pc(6,15)

sh(9,11) = 2*pc(7,13) - 2*pc(7,15) - pc(5,13) + pc(5,15) - pc(6,13) + pc(6,15)

sh(10,11) = pc(11,13) - pc(11,15)

sh(11,11) = pc(13,13) - pc(13,15) - pc(15,13) + pc(15,15)

sh(12,11) = pc(18,13) - pc(18,15) - 3*pc(14,13) + 3*pc(14,15)

sh(13,11) = 3*pc(12,13) - 3*pc(12,15) - pc(19,13) + pc(19,15)

sh(14,11) = 2*pc(20,13) - 2*pc(20,15) - 3*pc(13,13) + 3*pc(13,15) - 3*pc(15,13&
&) + 3*pc(15,15)

sh(15,11) = 4*pc(16,13) - 4*pc(16,15) - pc(18,13) + pc(18,15) - pc(14,13) + pc&
&(14,15)

sh(16,11) = 4*pc(17,13) - 4*pc(17,15) - pc(12,13) + pc(12,15) - pc(19,13) + pc&
&(19,15)

sh(1,12) = pc(1,18) - 3*pc(1,14)

sh(2,12) = pc(2,18) - 3*pc(2,14)

sh(3,12) = pc(3,18) - 3*pc(3,14)

sh(4,12) = pc(4,18) - 3*pc(4,14)

sh(5,12) = pc(8,18) - 3*pc(8,14)

sh(6,12) = pc(9,18) - 3*pc(9,14)

sh(7,12) = pc(10,18) - 3*pc(10,14)

sh(8,12) = pc(5,18) - 3*pc(5,14) - pc(6,18) + 3*pc(6,14)

sh(9,12) = 2*pc(7,18) - 6*pc(7,14) - pc(5,18) + 3*pc(5,14) - pc(6,18) + 3*pc(6&
&,14)

sh(10,12) = pc(11,18) - 3*pc(11,14)

sh(11,12) = pc(13,18) - 3*pc(13,14) - pc(15,18) + 3*pc(15,14)

sh(12,12) = pc(18,18) - 3*pc(18,14) - 3*pc(14,18) + 9*pc(14,14)

sh(13,12) = 3*pc(12,18) - 9*pc(12,14) - pc(19,18) + 3*pc(19,14)

sh(14,12) = 2*pc(20,18) - 6*pc(20,14) - 3*pc(13,18) + 9*pc(13,14) - 3*pc(15,18&
&) + 9*pc(15,14)

sh(15,12) = 4*pc(16,18) - 12*pc(16,14) - pc(18,18) + 3*pc(18,14) - pc(14,18) +&
& 3*pc(14,14)

sh(16,12) = 4*pc(17,18) - 12*pc(17,14) - pc(12,18) + 3*pc(12,14) - pc(19,18) +&
& 3*pc(19,14)

sh(1,13) = 3*pc(1,12) - pc(1,19)

sh(2,13) = 3*pc(2,12) - pc(2,19)

sh(3,13) = 3*pc(3,12) - pc(3,19)

sh(4,13) = 3*pc(4,12) - pc(4,19)

sh(5,13) = 3*pc(8,12) - pc(8,19)

sh(6,13) = 3*pc(9,12) - pc(9,19)

sh(7,13) = 3*pc(10,12) - pc(10,19)

sh(8,13) = 3*pc(5,12) - pc(5,19) - 3*pc(6,12) + pc(6,19)

sh(9,13) = 6*pc(7,12) - 2*pc(7,19) - 3*pc(5,12) + pc(5,19) - 3*pc(6,12) + pc(6&
&,19)

sh(10,13) = 3*pc(11,12) - pc(11,19)

sh(11,13) = 3*pc(13,12) - pc(13,19) - 3*pc(15,12) + pc(15,19)

sh(12,13) = 3*pc(18,12) - pc(18,19) - 9*pc(14,12) + 3*pc(14,19)

sh(13,13) = 9*pc(12,12) - 3*pc(12,19) - 3*pc(19,12) + pc(19,19)

sh(14,13) = 6*pc(20,12) - 2*pc(20,19) - 9*pc(13,12) + 3*pc(13,19) - 9*pc(15,12&
&) + 3*pc(15,19)

sh(15,13) = 12*pc(16,12) - 4*pc(16,19) - 3*pc(18,12) + pc(18,19) - 3*pc(14,12)&
& + pc(14,19)

sh(16,13) = 12*pc(17,12) - 4*pc(17,19) - 3*pc(12,12) + pc(12,19) - 3*pc(19,12)&
& + pc(19,19)

sh(1,14) = 2*pc(1,20) - 3*pc(1,13) - 3*pc(1,15)

sh(2,14) = 2*pc(2,20) - 3*pc(2,13) - 3*pc(2,15)

sh(3,14) = 2*pc(3,20) - 3*pc(3,13) - 3*pc(3,15)

sh(4,14) = 2*pc(4,20) - 3*pc(4,13) - 3*pc(4,15)

sh(5,14) = 2*pc(8,20) - 3*pc(8,13) - 3*pc(8,15)

sh(6,14) = 2*pc(9,20) - 3*pc(9,13) - 3*pc(9,15)

sh(7,14) = 2*pc(10,20) - 3*pc(10,13) - 3*pc(10,15)

sh(8,14) = 2*pc(5,20) - 3*pc(5,13) - 3*pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*&
&pc(6,15)

sh(9,14) = 4*pc(7,20) - 6*pc(7,13) - 6*pc(7,15) - 2*pc(5,20) + 3*pc(5,13) + 3*&
&pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*pc(6,15)

sh(10,14) = 2*pc(11,20) - 3*pc(11,13) - 3*pc(11,15)

sh(11,14) = 2*pc(13,20) - 3*pc(13,13) - 3*pc(13,15) - 2*pc(15,20) + 3*pc(15,13&
&) + 3*pc(15,15)

sh(12,14) = 2*pc(18,20) - 3*pc(18,13) - 3*pc(18,15) - 6*pc(14,20) + 9*pc(14,13&
&) + 9*pc(14,15)

sh(13,14) = 6*pc(12,20) - 9*pc(12,13) - 9*pc(12,15) - 2*pc(19,20) + 3*pc(19,13&
&) + 3*pc(19,15)

sh(14,14) = 4*pc(20,20) - 6*pc(20,13) - 6*pc(20,15) - 6*pc(13,20) + 9*pc(13,13&
&) + 9*pc(13,15) - 6*pc(15,20) + 9*pc(15,13) + 9*pc(15,15)

sh(15,14) = 8*pc(16,20) - 12*pc(16,13) - 12*pc(16,15) - 2*pc(18,20) + 3*pc(18,&
&13) + 3*pc(18,15) - 2*pc(14,20) + 3*pc(14,13) + 3*pc(14,15)

sh(16,14) = 8*pc(17,20) - 12*pc(17,13) - 12*pc(17,15) - 2*pc(12,20) + 3*pc(12,&
&13) + 3*pc(12,15) - 2*pc(19,20) + 3*pc(19,13) + 3*pc(19,15)

sh(1,15) = 4*pc(1,16) - pc(1,18) - pc(1,14)

sh(2,15) = 4*pc(2,16) - pc(2,18) - pc(2,14)

sh(3,15) = 4*pc(3,16) - pc(3,18) - pc(3,14)

sh(4,15) = 4*pc(4,16) - pc(4,18) - pc(4,14)

sh(5,15) = 4*pc(8,16) - pc(8,18) - pc(8,14)

sh(6,15) = 4*pc(9,16) - pc(9,18) - pc(9,14)

sh(7,15) = 4*pc(10,16) - pc(10,18) - pc(10,14)

sh(8,15) = 4*pc(5,16) - pc(5,18) - pc(5,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(9,15) = 8*pc(7,16) - 2*pc(7,18) - 2*pc(7,14) - 4*pc(5,16) + pc(5,18) + pc(5&
&,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(10,15) = 4*pc(11,16) - pc(11,18) - pc(11,14)

sh(11,15) = 4*pc(13,16) - pc(13,18) - pc(13,14) - 4*pc(15,16) + pc(15,18) + pc&
&(15,14)

sh(12,15) = 4*pc(18,16) - pc(18,18) - pc(18,14) - 12*pc(14,16) + 3*pc(14,18) +&
& 3*pc(14,14)

sh(13,15) = 12*pc(12,16) - 3*pc(12,18) - 3*pc(12,14) - 4*pc(19,16) + pc(19,18)&
& + pc(19,14)

sh(14,15) = 8*pc(20,16) - 2*pc(20,18) - 2*pc(20,14) - 12*pc(13,16) + 3*pc(13,1&
&8) + 3*pc(13,14) - 12*pc(15,16) + 3*pc(15,18) + 3*pc(15,14)

sh(15,15) = 16*pc(16,16) - 4*pc(16,18) - 4*pc(16,14) - 4*pc(18,16) + pc(18,18)&
& + pc(18,14) - 4*pc(14,16) + pc(14,18) + pc(14,14)

sh(16,15) = 16*pc(17,16) - 4*pc(17,18) - 4*pc(17,14) - 4*pc(12,16) + pc(12,18)&
& + pc(12,14) - 4*pc(19,16) + pc(19,18) + pc(19,14)

sh(1,16) = 4*pc(1,17) - pc(1,12) - pc(1,19)

sh(2,16) = 4*pc(2,17) - pc(2,12) - pc(2,19)

sh(3,16) = 4*pc(3,17) - pc(3,12) - pc(3,19)

sh(4,16) = 4*pc(4,17) - pc(4,12) - pc(4,19)

sh(5,16) = 4*pc(8,17) - pc(8,12) - pc(8,19)

sh(6,16) = 4*pc(9,17) - pc(9,12) - pc(9,19)

sh(7,16) = 4*pc(10,17) - pc(10,12) - pc(10,19)

sh(8,16) = 4*pc(5,17) - pc(5,12) - pc(5,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(9,16) = 8*pc(7,17) - 2*pc(7,12) - 2*pc(7,19) - 4*pc(5,17) + pc(5,12) + pc(5&
&,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(10,16) = 4*pc(11,17) - pc(11,12) - pc(11,19)

sh(11,16) = 4*pc(13,17) - pc(13,12) - pc(13,19) - 4*pc(15,17) + pc(15,12) + pc&
&(15,19)

sh(12,16) = 4*pc(18,17) - pc(18,12) - pc(18,19) - 12*pc(14,17) + 3*pc(14,12) +&
& 3*pc(14,19)

sh(13,16) = 12*pc(12,17) - 3*pc(12,12) - 3*pc(12,19) - 4*pc(19,17) + pc(19,12)&
& + pc(19,19)

sh(14,16) = 8*pc(20,17) - 2*pc(20,12) - 2*pc(20,19) - 12*pc(13,17) + 3*pc(13,1&
&2) + 3*pc(13,19) - 12*pc(15,17) + 3*pc(15,12) + 3*pc(15,19)

sh(15,16) = 16*pc(16,17) - 4*pc(16,12) - 4*pc(16,19) - 4*pc(18,17) + pc(18,12)&
& + pc(18,19) - 4*pc(14,17) + pc(14,12) + pc(14,19)

sh(16,16) = 16*pc(17,17) - 4*pc(17,12) - 4*pc(17,19) - 4*pc(12,17) + pc(12,12)&
& + pc(12,19) - 4*pc(19,17) + pc(19,12) + pc(19,19)


   end subroutine electron3CIntgNum

   subroutine nuclear3CIntgAna(a1,a2,a3,A,B,C,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2, a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   real (kind=double), dimension (20,20,7), intent(out) :: pc
   real (kind=double), dimension (16,16), intent(out) :: sh

   ! Define local variables
   integer :: m
   real (kind=double), dimension (3) :: G, GA, GB, GC, P, PC_3C, d1, d2, d3
   real (kind=double) :: zeta, zeta3C, inv_2zeta3C, U
   real (kind=double), dimension (7) :: preFactorN
   real (kind=double), dimension (7) :: F

   ! Initialize local variables.
   zeta = a1 + a2
   zeta3C = zeta + a3
   inv_2zeta3C = 1.0d0 / (2.0d0 * zeta3C)
   P = (a1*A + a2*B) / zeta
   G = (a1*A + a2*B + a3*C) / zeta3C
   GA = G - A
   GB = G - B
   GC = G - C
   PC_3C = P - C
   d1 = A - B
   d2 = A - C
   d3 = B - C
   U = zeta3C * sum(GC**2)
   call boys(U,F)

   do m = 1, 7
      preFactorN(m) = F(m) * 2.0d0*(pi/zeta3C) * &
            & exp(-sum(a1*a2*(d1**2) + a1*a3*(d2**2) + a2*a3*(d3**2))/zeta3C)
   enddo

pc(1,2,6) = GB(1)*preFactorN(6) - GC(1)*preFactorN(7)

pc(1,3,6) = GB(2)*preFactorN(6) - GC(2)*preFactorN(7)

pc(1,4,6) = GB(3)*preFactorN(6) - GC(3)*preFactorN(7)

pc(1,2,5) = GB(1)*preFactorN(5) - GC(1)*preFactorN(6)

pc(1,3,5) = GB(2)*preFactorN(5) - GC(2)*preFactorN(6)

pc(1,4,5) = GB(3)*preFactorN(5) - GC(3)*preFactorN(6)

pc(1,8,5) = GB(2)*pc(1,2,5) - GC(2)*pc(1,2,6)

pc(1,5,5) = GB(1)*pc(1,2,5) - GC(1)*pc(1,2,6) + inv_2zeta3C*((preFactorN(5) - &
&preFactorN(6)))

pc(1,6,5) = GB(2)*pc(1,3,5) - GC(2)*pc(1,3,6) + inv_2zeta3C*((preFactorN(5) - &
&preFactorN(6)))

pc(1,7,5) = GB(3)*pc(1,4,5) - GC(3)*pc(1,4,6) + inv_2zeta3C*((preFactorN(5) - &
&preFactorN(6)))

pc(1,2,4) = GB(1)*preFactorN(4) - GC(1)*preFactorN(5)

pc(1,3,4) = GB(2)*preFactorN(4) - GC(2)*preFactorN(5)

pc(1,4,4) = GB(3)*preFactorN(4) - GC(3)*preFactorN(5)

pc(1,5,4) = GB(1)*pc(1,2,4) - GC(1)*pc(1,2,5) + inv_2zeta3C*((preFactorN(4) - &
&preFactorN(5)))

pc(1,6,4) = GB(2)*pc(1,3,4) - GC(2)*pc(1,3,5) + inv_2zeta3C*((preFactorN(4) - &
&preFactorN(5)))

pc(1,7,4) = GB(3)*pc(1,4,4) - GC(3)*pc(1,4,5) + inv_2zeta3C*((preFactorN(4) - &
&preFactorN(5)))

pc(1,8,4) = GB(2)*pc(1,2,4) - GC(2)*pc(1,2,5)

pc(1,9,4) = GB(3)*pc(1,2,4) - GC(3)*pc(1,2,5)

pc(1,10,4) = GB(3)*pc(1,3,4) - GC(3)*pc(1,3,5)

pc(1,11,4) = GB(3)*pc(1,8,4) - GC(3)*pc(1,8,5)

pc(1,12,4) = GB(2)*pc(1,5,4) - GC(2)*pc(1,5,5)

pc(1,13,4) = GB(3)*pc(1,5,4) - GC(3)*pc(1,5,5)

pc(1,14,4) = GB(1)*pc(1,6,4) - GC(1)*pc(1,6,5)

pc(1,15,4) = GB(3)*pc(1,6,4) - GC(3)*pc(1,6,5)

pc(1,16,4) = GB(1)*pc(1,7,4) - GC(1)*pc(1,7,5)

pc(1,17,4) = GB(2)*pc(1,7,4) - GC(2)*pc(1,7,5)

pc(1,18,4) = GB(1)*pc(1,5,4) - GC(1)*pc(1,5,5) + inv_2zeta3C*(2*(pc(1,2,4) - p&
&c(1,2,5)))

pc(1,19,4) = GB(2)*pc(1,6,4) - GC(2)*pc(1,6,5) + inv_2zeta3C*(2*(pc(1,3,4) - p&
&c(1,3,5)))

pc(1,20,4) = GB(3)*pc(1,7,4) - GC(3)*pc(1,7,5) + inv_2zeta3C*(2*(pc(1,4,4) - p&
&c(1,4,5)))

pc(2,1,3) = GA(1)*preFactorN(3) - GC(1)*preFactorN(4)

pc(3,1,3) = GA(2)*preFactorN(3) - GC(2)*preFactorN(4)

pc(4,1,3) = GA(3)*preFactorN(3) - GC(3)*preFactorN(4)

pc(1,2,3) = GB(1)*preFactorN(3) - GC(1)*preFactorN(4)

pc(2,2,3) = GA(1)*pc(1,2,3) - GC(1)*pc(1,2,4) + inv_2zeta3C*((preFactorN(3) - &
&preFactorN(4)))

pc(3,2,3) = GA(2)*pc(1,2,3) - GC(2)*pc(1,2,4)

pc(4,2,3) = GA(3)*pc(1,2,3) - GC(3)*pc(1,2,4)

pc(1,3,3) = GB(2)*preFactorN(3) - GC(2)*preFactorN(4)

pc(2,3,3) = GA(1)*pc(1,3,3) - GC(1)*pc(1,3,4)

pc(3,3,3) = GA(2)*pc(1,3,3) - GC(2)*pc(1,3,4) + inv_2zeta3C*((preFactorN(3) - &
&preFactorN(4)))

pc(4,3,3) = GA(3)*pc(1,3,3) - GC(3)*pc(1,3,4)

pc(1,4,3) = GB(3)*preFactorN(3) - GC(3)*preFactorN(4)

pc(2,4,3) = GA(1)*pc(1,4,3) - GC(1)*pc(1,4,4)

pc(3,4,3) = GA(2)*pc(1,4,3) - GC(2)*pc(1,4,4)

pc(4,4,3) = GA(3)*pc(1,4,3) - GC(3)*pc(1,4,4) + inv_2zeta3C*((preFactorN(3) - &
&preFactorN(4)))

pc(1,5,3) = GB(1)*pc(1,2,3) - GC(1)*pc(1,2,4) + inv_2zeta3C*((preFactorN(3) - &
&preFactorN(4)))

pc(2,5,3) = GA(1)*pc(1,5,3) - GC(1)*pc(1,5,4) + inv_2zeta3C*(2*(pc(1,2,3) - pc&
&(1,2,4)))

pc(3,5,3) = GA(2)*pc(1,5,3) - GC(2)*pc(1,5,4)

pc(4,5,3) = GA(3)*pc(1,5,3) - GC(3)*pc(1,5,4)

pc(1,6,3) = GB(2)*pc(1,3,3) - GC(2)*pc(1,3,4) + inv_2zeta3C*((preFactorN(3) - &
&preFactorN(4)))

pc(2,6,3) = GA(1)*pc(1,6,3) - GC(1)*pc(1,6,4)

pc(3,6,3) = GA(2)*pc(1,6,3) - GC(2)*pc(1,6,4) + inv_2zeta3C*(2*(pc(1,3,3) - pc&
&(1,3,4)))

pc(4,6,3) = GA(3)*pc(1,6,3) - GC(3)*pc(1,6,4)

pc(1,7,3) = GB(3)*pc(1,4,3) - GC(3)*pc(1,4,4) + inv_2zeta3C*((preFactorN(3) - &
&preFactorN(4)))

pc(2,7,3) = GA(1)*pc(1,7,3) - GC(1)*pc(1,7,4)

pc(3,7,3) = GA(2)*pc(1,7,3) - GC(2)*pc(1,7,4)

pc(4,7,3) = GA(3)*pc(1,7,3) - GC(3)*pc(1,7,4) + inv_2zeta3C*(2*(pc(1,4,3) - pc&
&(1,4,4)))

pc(1,8,3) = GB(2)*pc(1,2,3) - GC(2)*pc(1,2,4)

pc(2,8,3) = GA(1)*pc(1,8,3) - GC(1)*pc(1,8,4) + inv_2zeta3C*((pc(1,3,3) - pc(1&
&,3,4)))

pc(3,8,3) = GA(2)*pc(1,8,3) - GC(2)*pc(1,8,4) + inv_2zeta3C*((pc(1,2,3) - pc(1&
&,2,4)))

pc(4,8,3) = GA(3)*pc(1,8,3) - GC(3)*pc(1,8,4)

pc(1,9,3) = GB(3)*pc(1,2,3) - GC(3)*pc(1,2,4)

pc(2,9,3) = GA(1)*pc(1,9,3) - GC(1)*pc(1,9,4) + inv_2zeta3C*((pc(1,4,3) - pc(1&
&,4,4)))

pc(3,9,3) = GA(2)*pc(1,9,3) - GC(2)*pc(1,9,4)

pc(4,9,3) = GA(3)*pc(1,9,3) - GC(3)*pc(1,9,4) + inv_2zeta3C*((pc(1,2,3) - pc(1&
&,2,4)))

pc(1,10,3) = GB(3)*pc(1,3,3) - GC(3)*pc(1,3,4)

pc(2,10,3) = GA(1)*pc(1,10,3) - GC(1)*pc(1,10,4)

pc(3,10,3) = GA(2)*pc(1,10,3) - GC(2)*pc(1,10,4) + inv_2zeta3C*((pc(1,4,3) - p&
&c(1,4,4)))

pc(4,10,3) = GA(3)*pc(1,10,3) - GC(3)*pc(1,10,4) + inv_2zeta3C*((pc(1,3,3) - p&
&c(1,3,4)))

pc(1,11,3) = GB(3)*pc(1,8,3) - GC(3)*pc(1,8,4)

pc(2,11,3) = GA(1)*pc(1,11,3) - GC(1)*pc(1,11,4) + inv_2zeta3C*((pc(1,10,3) - &
&pc(1,10,4)))

pc(3,11,3) = GA(2)*pc(1,11,3) - GC(2)*pc(1,11,4) + inv_2zeta3C*((pc(1,9,3) - p&
&c(1,9,4)))

pc(4,11,3) = GA(3)*pc(1,11,3) - GC(3)*pc(1,11,4) + inv_2zeta3C*((pc(1,8,3) - p&
&c(1,8,4)))

pc(1,12,3) = GB(2)*pc(1,5,3) - GC(2)*pc(1,5,4)

pc(2,12,3) = GA(1)*pc(1,12,3) - GC(1)*pc(1,12,4) + inv_2zeta3C*(2*(pc(1,8,3) -&
& pc(1,8,4)))

pc(3,12,3) = GA(2)*pc(1,12,3) - GC(2)*pc(1,12,4) + inv_2zeta3C*((pc(1,5,3) - p&
&c(1,5,4)))

pc(4,12,3) = GA(3)*pc(1,12,3) - GC(3)*pc(1,12,4)

pc(1,13,3) = GB(3)*pc(1,5,3) - GC(3)*pc(1,5,4)

pc(2,13,3) = GA(1)*pc(1,13,3) - GC(1)*pc(1,13,4) + inv_2zeta3C*(2*(pc(1,9,3) -&
& pc(1,9,4)))

pc(3,13,3) = GA(2)*pc(1,13,3) - GC(2)*pc(1,13,4)

pc(4,13,3) = GA(3)*pc(1,13,3) - GC(3)*pc(1,13,4) + inv_2zeta3C*((pc(1,5,3) - p&
&c(1,5,4)))

pc(1,14,3) = GB(1)*pc(1,6,3) - GC(1)*pc(1,6,4)

pc(2,14,3) = GA(1)*pc(1,14,3) - GC(1)*pc(1,14,4) + inv_2zeta3C*((pc(1,6,3) - p&
&c(1,6,4)))

pc(3,14,3) = GA(2)*pc(1,14,3) - GC(2)*pc(1,14,4) + inv_2zeta3C*(2*(pc(1,8,3) -&
& pc(1,8,4)))

pc(4,14,3) = GA(3)*pc(1,14,3) - GC(3)*pc(1,14,4)

pc(1,15,3) = GB(3)*pc(1,6,3) - GC(3)*pc(1,6,4)

pc(2,15,3) = GA(1)*pc(1,15,3) - GC(1)*pc(1,15,4)

pc(3,15,3) = GA(2)*pc(1,15,3) - GC(2)*pc(1,15,4) + inv_2zeta3C*(2*(pc(1,10,3) &
&- pc(1,10,4)))

pc(4,15,3) = GA(3)*pc(1,15,3) - GC(3)*pc(1,15,4) + inv_2zeta3C*((pc(1,6,3) - p&
&c(1,6,4)))

pc(1,16,3) = GB(1)*pc(1,7,3) - GC(1)*pc(1,7,4)

pc(2,16,3) = GA(1)*pc(1,16,3) - GC(1)*pc(1,16,4) + inv_2zeta3C*((pc(1,7,3) - p&
&c(1,7,4)))

pc(3,16,3) = GA(2)*pc(1,16,3) - GC(2)*pc(1,16,4)

pc(4,16,3) = GA(3)*pc(1,16,3) - GC(3)*pc(1,16,4) + inv_2zeta3C*(2*(pc(1,9,3) -&
& pc(1,9,4)))

pc(1,17,3) = GB(2)*pc(1,7,3) - GC(2)*pc(1,7,4)

pc(2,17,3) = GA(1)*pc(1,17,3) - GC(1)*pc(1,17,4)

pc(3,17,3) = GA(2)*pc(1,17,3) - GC(2)*pc(1,17,4) + inv_2zeta3C*((pc(1,7,3) - p&
&c(1,7,4)))

pc(4,17,3) = GA(3)*pc(1,17,3) - GC(3)*pc(1,17,4) + inv_2zeta3C*(2*(pc(1,10,3) &
&- pc(1,10,4)))

pc(1,18,3) = GB(1)*pc(1,5,3) - GC(1)*pc(1,5,4) + inv_2zeta3C*(2*(pc(1,2,3) - p&
&c(1,2,4)))

pc(2,18,3) = GA(1)*pc(1,18,3) - GC(1)*pc(1,18,4) + inv_2zeta3C*(3*(pc(1,5,3) -&
& pc(1,5,4)))

pc(3,18,3) = GA(2)*pc(1,18,3) - GC(2)*pc(1,18,4)

pc(4,18,3) = GA(3)*pc(1,18,3) - GC(3)*pc(1,18,4)

pc(1,19,3) = GB(2)*pc(1,6,3) - GC(2)*pc(1,6,4) + inv_2zeta3C*(2*(pc(1,3,3) - p&
&c(1,3,4)))

pc(2,19,3) = GA(1)*pc(1,19,3) - GC(1)*pc(1,19,4)

pc(3,19,3) = GA(2)*pc(1,19,3) - GC(2)*pc(1,19,4) + inv_2zeta3C*(3*(pc(1,6,3) -&
& pc(1,6,4)))

pc(4,19,3) = GA(3)*pc(1,19,3) - GC(3)*pc(1,19,4)

pc(1,20,3) = GB(3)*pc(1,7,3) - GC(3)*pc(1,7,4) + inv_2zeta3C*(2*(pc(1,4,3) - p&
&c(1,4,4)))

pc(2,20,3) = GA(1)*pc(1,20,3) - GC(1)*pc(1,20,4)

pc(3,20,3) = GA(2)*pc(1,20,3) - GC(2)*pc(1,20,4)

pc(4,20,3) = GA(3)*pc(1,20,3) - GC(3)*pc(1,20,4) + inv_2zeta3C*(3*(pc(1,7,3) -&
& pc(1,7,4)))

pc(2,1,2) = GA(1)*preFactorN(2) - GC(1)*preFactorN(3)

pc(3,1,2) = GA(2)*preFactorN(2) - GC(2)*preFactorN(3)

pc(4,1,2) = GA(3)*preFactorN(2) - GC(3)*preFactorN(3)

pc(8,1,2) = GA(2)*pc(2,1,2) - GC(2)*pc(2,1,3)

pc(5,1,2) = GA(1)*pc(2,1,2) - GC(1)*pc(2,1,3) + inv_2zeta3C*((preFactorN(2) - &
&preFactorN(3)))

pc(6,1,2) = GA(2)*pc(3,1,2) - GC(2)*pc(3,1,3) + inv_2zeta3C*((preFactorN(2) - &
&preFactorN(3)))

pc(7,1,2) = GA(3)*pc(4,1,2) - GC(3)*pc(4,1,3) + inv_2zeta3C*((preFactorN(2) - &
&preFactorN(3)))

pc(1,2,2) = GB(1)*preFactorN(2) - GC(1)*preFactorN(3)

pc(2,2,2) = GA(1)*pc(1,2,2) - GC(1)*pc(1,2,3) + inv_2zeta3C*((preFactorN(2) - &
&preFactorN(3)))

pc(3,2,2) = GA(2)*pc(1,2,2) - GC(2)*pc(1,2,3)

pc(4,2,2) = GA(3)*pc(1,2,2) - GC(3)*pc(1,2,3)

pc(8,2,2) = GA(2)*pc(2,2,2) - GC(2)*pc(2,2,3)

pc(5,2,2) = GA(1)*pc(2,2,2) - GC(1)*pc(2,2,3) + inv_2zeta3C*((pc(1,2,2) - pc(1&
&,2,3)) + (pc(2,1,2) - pc(2,1,3)))

pc(6,2,2) = GA(2)*pc(3,2,2) - GC(2)*pc(3,2,3) + inv_2zeta3C*((pc(1,2,2) - pc(1&
&,2,3)))

pc(7,2,2) = GA(3)*pc(4,2,2) - GC(3)*pc(4,2,3) + inv_2zeta3C*((pc(1,2,2) - pc(1&
&,2,3)))

pc(1,3,2) = GB(2)*preFactorN(2) - GC(2)*preFactorN(3)

pc(2,3,2) = GA(1)*pc(1,3,2) - GC(1)*pc(1,3,3)

pc(3,3,2) = GA(2)*pc(1,3,2) - GC(2)*pc(1,3,3) + inv_2zeta3C*((preFactorN(2) - &
&preFactorN(3)))

pc(4,3,2) = GA(3)*pc(1,3,2) - GC(3)*pc(1,3,3)

pc(8,3,2) = GA(2)*pc(2,3,2) - GC(2)*pc(2,3,3) + inv_2zeta3C*((pc(2,1,2) - pc(2&
&,1,3)))

pc(5,3,2) = GA(1)*pc(2,3,2) - GC(1)*pc(2,3,3) + inv_2zeta3C*((pc(1,3,2) - pc(1&
&,3,3)))

pc(6,3,2) = GA(2)*pc(3,3,2) - GC(2)*pc(3,3,3) + inv_2zeta3C*((pc(1,3,2) - pc(1&
&,3,3)) + (pc(3,1,2) - pc(3,1,3)))

pc(7,3,2) = GA(3)*pc(4,3,2) - GC(3)*pc(4,3,3) + inv_2zeta3C*((pc(1,3,2) - pc(1&
&,3,3)))

pc(1,4,2) = GB(3)*preFactorN(2) - GC(3)*preFactorN(3)

pc(2,4,2) = GA(1)*pc(1,4,2) - GC(1)*pc(1,4,3)

pc(3,4,2) = GA(2)*pc(1,4,2) - GC(2)*pc(1,4,3)

pc(4,4,2) = GA(3)*pc(1,4,2) - GC(3)*pc(1,4,3) + inv_2zeta3C*((preFactorN(2) - &
&preFactorN(3)))

pc(8,4,2) = GA(2)*pc(2,4,2) - GC(2)*pc(2,4,3)

pc(5,4,2) = GA(1)*pc(2,4,2) - GC(1)*pc(2,4,3) + inv_2zeta3C*((pc(1,4,2) - pc(1&
&,4,3)))

pc(6,4,2) = GA(2)*pc(3,4,2) - GC(2)*pc(3,4,3) + inv_2zeta3C*((pc(1,4,2) - pc(1&
&,4,3)))

pc(7,4,2) = GA(3)*pc(4,4,2) - GC(3)*pc(4,4,3) + inv_2zeta3C*((pc(1,4,2) - pc(1&
&,4,3)) + (pc(4,1,2) - pc(4,1,3)))

pc(1,5,2) = GB(1)*pc(1,2,2) - GC(1)*pc(1,2,3) + inv_2zeta3C*((preFactorN(2) - &
&preFactorN(3)))

pc(2,5,2) = GA(1)*pc(1,5,2) - GC(1)*pc(1,5,3) + inv_2zeta3C*(2*(pc(1,2,2) - pc&
&(1,2,3)))

pc(3,5,2) = GA(2)*pc(1,5,2) - GC(2)*pc(1,5,3)

pc(4,5,2) = GA(3)*pc(1,5,2) - GC(3)*pc(1,5,3)

pc(8,5,2) = GA(2)*pc(2,5,2) - GC(2)*pc(2,5,3)

pc(5,5,2) = GA(1)*pc(2,5,2) - GC(1)*pc(2,5,3) + inv_2zeta3C*((pc(1,5,2) - pc(1&
&,5,3)) + 2*(pc(2,2,2) - pc(2,2,3)))

pc(6,5,2) = GA(2)*pc(3,5,2) - GC(2)*pc(3,5,3) + inv_2zeta3C*((pc(1,5,2) - pc(1&
&,5,3)))

pc(7,5,2) = GA(3)*pc(4,5,2) - GC(3)*pc(4,5,3) + inv_2zeta3C*((pc(1,5,2) - pc(1&
&,5,3)))

pc(1,6,2) = GB(2)*pc(1,3,2) - GC(2)*pc(1,3,3) + inv_2zeta3C*((preFactorN(2) - &
&preFactorN(3)))

pc(2,6,2) = GA(1)*pc(1,6,2) - GC(1)*pc(1,6,3)

pc(3,6,2) = GA(2)*pc(1,6,2) - GC(2)*pc(1,6,3) + inv_2zeta3C*(2*(pc(1,3,2) - pc&
&(1,3,3)))

pc(4,6,2) = GA(3)*pc(1,6,2) - GC(3)*pc(1,6,3)

pc(8,6,2) = GA(2)*pc(2,6,2) - GC(2)*pc(2,6,3) + inv_2zeta3C*(2*(pc(2,3,2) - pc&
&(2,3,3)))

pc(5,6,2) = GA(1)*pc(2,6,2) - GC(1)*pc(2,6,3) + inv_2zeta3C*((pc(1,6,2) - pc(1&
&,6,3)))

pc(6,6,2) = GA(2)*pc(3,6,2) - GC(2)*pc(3,6,3) + inv_2zeta3C*((pc(1,6,2) - pc(1&
&,6,3)) + 2*(pc(3,3,2) - pc(3,3,3)))

pc(7,6,2) = GA(3)*pc(4,6,2) - GC(3)*pc(4,6,3) + inv_2zeta3C*((pc(1,6,2) - pc(1&
&,6,3)))

pc(1,7,2) = GB(3)*pc(1,4,2) - GC(3)*pc(1,4,3) + inv_2zeta3C*((preFactorN(2) - &
&preFactorN(3)))

pc(2,7,2) = GA(1)*pc(1,7,2) - GC(1)*pc(1,7,3)

pc(3,7,2) = GA(2)*pc(1,7,2) - GC(2)*pc(1,7,3)

pc(4,7,2) = GA(3)*pc(1,7,2) - GC(3)*pc(1,7,3) + inv_2zeta3C*(2*(pc(1,4,2) - pc&
&(1,4,3)))

pc(8,7,2) = GA(2)*pc(2,7,2) - GC(2)*pc(2,7,3)

pc(5,7,2) = GA(1)*pc(2,7,2) - GC(1)*pc(2,7,3) + inv_2zeta3C*((pc(1,7,2) - pc(1&
&,7,3)))

pc(6,7,2) = GA(2)*pc(3,7,2) - GC(2)*pc(3,7,3) + inv_2zeta3C*((pc(1,7,2) - pc(1&
&,7,3)))

pc(7,7,2) = GA(3)*pc(4,7,2) - GC(3)*pc(4,7,3) + inv_2zeta3C*((pc(1,7,2) - pc(1&
&,7,3)) + 2*(pc(4,4,2) - pc(4,4,3)))

pc(1,8,2) = GB(2)*pc(1,2,2) - GC(2)*pc(1,2,3)

pc(2,8,2) = GA(1)*pc(1,8,2) - GC(1)*pc(1,8,3) + inv_2zeta3C*((pc(1,3,2) - pc(1&
&,3,3)))

pc(3,8,2) = GA(2)*pc(1,8,2) - GC(2)*pc(1,8,3) + inv_2zeta3C*((pc(1,2,2) - pc(1&
&,2,3)))

pc(4,8,2) = GA(3)*pc(1,8,2) - GC(3)*pc(1,8,3)

pc(8,8,2) = GA(2)*pc(2,8,2) - GC(2)*pc(2,8,3) + inv_2zeta3C*((pc(2,2,2) - pc(2&
&,2,3)))

pc(5,8,2) = GA(1)*pc(2,8,2) - GC(1)*pc(2,8,3) + inv_2zeta3C*((pc(1,8,2) - pc(1&
&,8,3)) + (pc(2,3,2) - pc(2,3,3)))

pc(6,8,2) = GA(2)*pc(3,8,2) - GC(2)*pc(3,8,3) + inv_2zeta3C*((pc(1,8,2) - pc(1&
&,8,3)) + (pc(3,2,2) - pc(3,2,3)))

pc(7,8,2) = GA(3)*pc(4,8,2) - GC(3)*pc(4,8,3) + inv_2zeta3C*((pc(1,8,2) - pc(1&
&,8,3)))

pc(1,9,2) = GB(3)*pc(1,2,2) - GC(3)*pc(1,2,3)

pc(2,9,2) = GA(1)*pc(1,9,2) - GC(1)*pc(1,9,3) + inv_2zeta3C*((pc(1,4,2) - pc(1&
&,4,3)))

pc(3,9,2) = GA(2)*pc(1,9,2) - GC(2)*pc(1,9,3)

pc(4,9,2) = GA(3)*pc(1,9,2) - GC(3)*pc(1,9,3) + inv_2zeta3C*((pc(1,2,2) - pc(1&
&,2,3)))

pc(8,9,2) = GA(2)*pc(2,9,2) - GC(2)*pc(2,9,3)

pc(5,9,2) = GA(1)*pc(2,9,2) - GC(1)*pc(2,9,3) + inv_2zeta3C*((pc(1,9,2) - pc(1&
&,9,3)) + (pc(2,4,2) - pc(2,4,3)))

pc(6,9,2) = GA(2)*pc(3,9,2) - GC(2)*pc(3,9,3) + inv_2zeta3C*((pc(1,9,2) - pc(1&
&,9,3)))

pc(7,9,2) = GA(3)*pc(4,9,2) - GC(3)*pc(4,9,3) + inv_2zeta3C*((pc(1,9,2) - pc(1&
&,9,3)) + (pc(4,2,2) - pc(4,2,3)))

pc(1,10,2) = GB(3)*pc(1,3,2) - GC(3)*pc(1,3,3)

pc(2,10,2) = GA(1)*pc(1,10,2) - GC(1)*pc(1,10,3)

pc(3,10,2) = GA(2)*pc(1,10,2) - GC(2)*pc(1,10,3) + inv_2zeta3C*((pc(1,4,2) - p&
&c(1,4,3)))

pc(4,10,2) = GA(3)*pc(1,10,2) - GC(3)*pc(1,10,3) + inv_2zeta3C*((pc(1,3,2) - p&
&c(1,3,3)))

pc(8,10,2) = GA(2)*pc(2,10,2) - GC(2)*pc(2,10,3) + inv_2zeta3C*((pc(2,4,2) - p&
&c(2,4,3)))

pc(5,10,2) = GA(1)*pc(2,10,2) - GC(1)*pc(2,10,3) + inv_2zeta3C*((pc(1,10,2) - &
&pc(1,10,3)))

pc(6,10,2) = GA(2)*pc(3,10,2) - GC(2)*pc(3,10,3) + inv_2zeta3C*((pc(1,10,2) - &
&pc(1,10,3)) + (pc(3,4,2) - pc(3,4,3)))

pc(7,10,2) = GA(3)*pc(4,10,2) - GC(3)*pc(4,10,3) + inv_2zeta3C*((pc(1,10,2) - &
&pc(1,10,3)) + (pc(4,3,2) - pc(4,3,3)))

pc(1,11,2) = GB(3)*pc(1,8,2) - GC(3)*pc(1,8,3)

pc(2,11,2) = GA(1)*pc(1,11,2) - GC(1)*pc(1,11,3) + inv_2zeta3C*((pc(1,10,2) - &
&pc(1,10,3)))

pc(3,11,2) = GA(2)*pc(1,11,2) - GC(2)*pc(1,11,3) + inv_2zeta3C*((pc(1,9,2) - p&
&c(1,9,3)))

pc(4,11,2) = GA(3)*pc(1,11,2) - GC(3)*pc(1,11,3) + inv_2zeta3C*((pc(1,8,2) - p&
&c(1,8,3)))

pc(8,11,2) = GA(2)*pc(2,11,2) - GC(2)*pc(2,11,3) + inv_2zeta3C*((pc(2,9,2) - p&
&c(2,9,3)))

pc(5,11,2) = GA(1)*pc(2,11,2) - GC(1)*pc(2,11,3) + inv_2zeta3C*((pc(1,11,2) - &
&pc(1,11,3)) + (pc(2,10,2) - pc(2,10,3)))

pc(6,11,2) = GA(2)*pc(3,11,2) - GC(2)*pc(3,11,3) + inv_2zeta3C*((pc(1,11,2) - &
&pc(1,11,3)) + (pc(3,9,2) - pc(3,9,3)))

pc(7,11,2) = GA(3)*pc(4,11,2) - GC(3)*pc(4,11,3) + inv_2zeta3C*((pc(1,11,2) - &
&pc(1,11,3)) + (pc(4,8,2) - pc(4,8,3)))

pc(1,12,2) = GB(2)*pc(1,5,2) - GC(2)*pc(1,5,3)

pc(2,12,2) = GA(1)*pc(1,12,2) - GC(1)*pc(1,12,3) + inv_2zeta3C*(2*(pc(1,8,2) -&
& pc(1,8,3)))

pc(3,12,2) = GA(2)*pc(1,12,2) - GC(2)*pc(1,12,3) + inv_2zeta3C*((pc(1,5,2) - p&
&c(1,5,3)))

pc(4,12,2) = GA(3)*pc(1,12,2) - GC(3)*pc(1,12,3)

pc(8,12,2) = GA(2)*pc(2,12,2) - GC(2)*pc(2,12,3) + inv_2zeta3C*((pc(2,5,2) - p&
&c(2,5,3)))

pc(5,12,2) = GA(1)*pc(2,12,2) - GC(1)*pc(2,12,3) + inv_2zeta3C*((pc(1,12,2) - &
&pc(1,12,3)) + 2*(pc(2,8,2) - pc(2,8,3)))

pc(6,12,2) = GA(2)*pc(3,12,2) - GC(2)*pc(3,12,3) + inv_2zeta3C*((pc(1,12,2) - &
&pc(1,12,3)) + (pc(3,5,2) - pc(3,5,3)))

pc(7,12,2) = GA(3)*pc(4,12,2) - GC(3)*pc(4,12,3) + inv_2zeta3C*((pc(1,12,2) - &
&pc(1,12,3)))

pc(1,13,2) = GB(3)*pc(1,5,2) - GC(3)*pc(1,5,3)

pc(2,13,2) = GA(1)*pc(1,13,2) - GC(1)*pc(1,13,3) + inv_2zeta3C*(2*(pc(1,9,2) -&
& pc(1,9,3)))

pc(3,13,2) = GA(2)*pc(1,13,2) - GC(2)*pc(1,13,3)

pc(4,13,2) = GA(3)*pc(1,13,2) - GC(3)*pc(1,13,3) + inv_2zeta3C*((pc(1,5,2) - p&
&c(1,5,3)))

pc(8,13,2) = GA(2)*pc(2,13,2) - GC(2)*pc(2,13,3)

pc(5,13,2) = GA(1)*pc(2,13,2) - GC(1)*pc(2,13,3) + inv_2zeta3C*((pc(1,13,2) - &
&pc(1,13,3)) + 2*(pc(2,9,2) - pc(2,9,3)))

pc(6,13,2) = GA(2)*pc(3,13,2) - GC(2)*pc(3,13,3) + inv_2zeta3C*((pc(1,13,2) - &
&pc(1,13,3)))

pc(7,13,2) = GA(3)*pc(4,13,2) - GC(3)*pc(4,13,3) + inv_2zeta3C*((pc(1,13,2) - &
&pc(1,13,3)) + (pc(4,5,2) - pc(4,5,3)))

pc(1,14,2) = GB(1)*pc(1,6,2) - GC(1)*pc(1,6,3)

pc(2,14,2) = GA(1)*pc(1,14,2) - GC(1)*pc(1,14,3) + inv_2zeta3C*((pc(1,6,2) - p&
&c(1,6,3)))

pc(3,14,2) = GA(2)*pc(1,14,2) - GC(2)*pc(1,14,3) + inv_2zeta3C*(2*(pc(1,8,2) -&
& pc(1,8,3)))

pc(4,14,2) = GA(3)*pc(1,14,2) - GC(3)*pc(1,14,3)

pc(8,14,2) = GA(2)*pc(2,14,2) - GC(2)*pc(2,14,3) + inv_2zeta3C*(2*(pc(2,8,2) -&
& pc(2,8,3)))

pc(5,14,2) = GA(1)*pc(2,14,2) - GC(1)*pc(2,14,3) + inv_2zeta3C*((pc(1,14,2) - &
&pc(1,14,3)) + (pc(2,6,2) - pc(2,6,3)))

pc(6,14,2) = GA(2)*pc(3,14,2) - GC(2)*pc(3,14,3) + inv_2zeta3C*((pc(1,14,2) - &
&pc(1,14,3)) + 2*(pc(3,8,2) - pc(3,8,3)))

pc(7,14,2) = GA(3)*pc(4,14,2) - GC(3)*pc(4,14,3) + inv_2zeta3C*((pc(1,14,2) - &
&pc(1,14,3)))

pc(1,15,2) = GB(3)*pc(1,6,2) - GC(3)*pc(1,6,3)

pc(2,15,2) = GA(1)*pc(1,15,2) - GC(1)*pc(1,15,3)

pc(3,15,2) = GA(2)*pc(1,15,2) - GC(2)*pc(1,15,3) + inv_2zeta3C*(2*(pc(1,10,2) &
&- pc(1,10,3)))

pc(4,15,2) = GA(3)*pc(1,15,2) - GC(3)*pc(1,15,3) + inv_2zeta3C*((pc(1,6,2) - p&
&c(1,6,3)))

pc(8,15,2) = GA(2)*pc(2,15,2) - GC(2)*pc(2,15,3) + inv_2zeta3C*(2*(pc(2,10,2) &
&- pc(2,10,3)))

pc(5,15,2) = GA(1)*pc(2,15,2) - GC(1)*pc(2,15,3) + inv_2zeta3C*((pc(1,15,2) - &
&pc(1,15,3)))

pc(6,15,2) = GA(2)*pc(3,15,2) - GC(2)*pc(3,15,3) + inv_2zeta3C*((pc(1,15,2) - &
&pc(1,15,3)) + 2*(pc(3,10,2) - pc(3,10,3)))

pc(7,15,2) = GA(3)*pc(4,15,2) - GC(3)*pc(4,15,3) + inv_2zeta3C*((pc(1,15,2) - &
&pc(1,15,3)) + (pc(4,6,2) - pc(4,6,3)))

pc(1,16,2) = GB(1)*pc(1,7,2) - GC(1)*pc(1,7,3)

pc(2,16,2) = GA(1)*pc(1,16,2) - GC(1)*pc(1,16,3) + inv_2zeta3C*((pc(1,7,2) - p&
&c(1,7,3)))

pc(3,16,2) = GA(2)*pc(1,16,2) - GC(2)*pc(1,16,3)

pc(4,16,2) = GA(3)*pc(1,16,2) - GC(3)*pc(1,16,3) + inv_2zeta3C*(2*(pc(1,9,2) -&
& pc(1,9,3)))

pc(8,16,2) = GA(2)*pc(2,16,2) - GC(2)*pc(2,16,3)

pc(5,16,2) = GA(1)*pc(2,16,2) - GC(1)*pc(2,16,3) + inv_2zeta3C*((pc(1,16,2) - &
&pc(1,16,3)) + (pc(2,7,2) - pc(2,7,3)))

pc(6,16,2) = GA(2)*pc(3,16,2) - GC(2)*pc(3,16,3) + inv_2zeta3C*((pc(1,16,2) - &
&pc(1,16,3)))

pc(7,16,2) = GA(3)*pc(4,16,2) - GC(3)*pc(4,16,3) + inv_2zeta3C*((pc(1,16,2) - &
&pc(1,16,3)) + 2*(pc(4,9,2) - pc(4,9,3)))

pc(1,17,2) = GB(2)*pc(1,7,2) - GC(2)*pc(1,7,3)

pc(2,17,2) = GA(1)*pc(1,17,2) - GC(1)*pc(1,17,3)

pc(3,17,2) = GA(2)*pc(1,17,2) - GC(2)*pc(1,17,3) + inv_2zeta3C*((pc(1,7,2) - p&
&c(1,7,3)))

pc(4,17,2) = GA(3)*pc(1,17,2) - GC(3)*pc(1,17,3) + inv_2zeta3C*(2*(pc(1,10,2) &
&- pc(1,10,3)))

pc(8,17,2) = GA(2)*pc(2,17,2) - GC(2)*pc(2,17,3) + inv_2zeta3C*((pc(2,7,2) - p&
&c(2,7,3)))

pc(5,17,2) = GA(1)*pc(2,17,2) - GC(1)*pc(2,17,3) + inv_2zeta3C*((pc(1,17,2) - &
&pc(1,17,3)))

pc(6,17,2) = GA(2)*pc(3,17,2) - GC(2)*pc(3,17,3) + inv_2zeta3C*((pc(1,17,2) - &
&pc(1,17,3)) + (pc(3,7,2) - pc(3,7,3)))

pc(7,17,2) = GA(3)*pc(4,17,2) - GC(3)*pc(4,17,3) + inv_2zeta3C*((pc(1,17,2) - &
&pc(1,17,3)) + 2*(pc(4,10,2) - pc(4,10,3)))

pc(1,18,2) = GB(1)*pc(1,5,2) - GC(1)*pc(1,5,3) + inv_2zeta3C*(2*(pc(1,2,2) - p&
&c(1,2,3)))

pc(2,18,2) = GA(1)*pc(1,18,2) - GC(1)*pc(1,18,3) + inv_2zeta3C*(3*(pc(1,5,2) -&
& pc(1,5,3)))

pc(3,18,2) = GA(2)*pc(1,18,2) - GC(2)*pc(1,18,3)

pc(4,18,2) = GA(3)*pc(1,18,2) - GC(3)*pc(1,18,3)

pc(8,18,2) = GA(2)*pc(2,18,2) - GC(2)*pc(2,18,3)

pc(5,18,2) = GA(1)*pc(2,18,2) - GC(1)*pc(2,18,3) + inv_2zeta3C*((pc(1,18,2) - &
&pc(1,18,3)) + 3*(pc(2,5,2) - pc(2,5,3)))

pc(6,18,2) = GA(2)*pc(3,18,2) - GC(2)*pc(3,18,3) + inv_2zeta3C*((pc(1,18,2) - &
&pc(1,18,3)))

pc(7,18,2) = GA(3)*pc(4,18,2) - GC(3)*pc(4,18,3) + inv_2zeta3C*((pc(1,18,2) - &
&pc(1,18,3)))

pc(1,19,2) = GB(2)*pc(1,6,2) - GC(2)*pc(1,6,3) + inv_2zeta3C*(2*(pc(1,3,2) - p&
&c(1,3,3)))

pc(2,19,2) = GA(1)*pc(1,19,2) - GC(1)*pc(1,19,3)

pc(3,19,2) = GA(2)*pc(1,19,2) - GC(2)*pc(1,19,3) + inv_2zeta3C*(3*(pc(1,6,2) -&
& pc(1,6,3)))

pc(4,19,2) = GA(3)*pc(1,19,2) - GC(3)*pc(1,19,3)

pc(8,19,2) = GA(2)*pc(2,19,2) - GC(2)*pc(2,19,3) + inv_2zeta3C*(3*(pc(2,6,2) -&
& pc(2,6,3)))

pc(5,19,2) = GA(1)*pc(2,19,2) - GC(1)*pc(2,19,3) + inv_2zeta3C*((pc(1,19,2) - &
&pc(1,19,3)))

pc(6,19,2) = GA(2)*pc(3,19,2) - GC(2)*pc(3,19,3) + inv_2zeta3C*((pc(1,19,2) - &
&pc(1,19,3)) + 3*(pc(3,6,2) - pc(3,6,3)))

pc(7,19,2) = GA(3)*pc(4,19,2) - GC(3)*pc(4,19,3) + inv_2zeta3C*((pc(1,19,2) - &
&pc(1,19,3)))

pc(1,20,2) = GB(3)*pc(1,7,2) - GC(3)*pc(1,7,3) + inv_2zeta3C*(2*(pc(1,4,2) - p&
&c(1,4,3)))

pc(2,20,2) = GA(1)*pc(1,20,2) - GC(1)*pc(1,20,3)

pc(3,20,2) = GA(2)*pc(1,20,2) - GC(2)*pc(1,20,3)

pc(4,20,2) = GA(3)*pc(1,20,2) - GC(3)*pc(1,20,3) + inv_2zeta3C*(3*(pc(1,7,2) -&
& pc(1,7,3)))

pc(8,20,2) = GA(2)*pc(2,20,2) - GC(2)*pc(2,20,3)

pc(5,20,2) = GA(1)*pc(2,20,2) - GC(1)*pc(2,20,3) + inv_2zeta3C*((pc(1,20,2) - &
&pc(1,20,3)))

pc(6,20,2) = GA(2)*pc(3,20,2) - GC(2)*pc(3,20,3) + inv_2zeta3C*((pc(1,20,2) - &
&pc(1,20,3)))

pc(7,20,2) = GA(3)*pc(4,20,2) - GC(3)*pc(4,20,3) + inv_2zeta3C*((pc(1,20,2) - &
&pc(1,20,3)) + 3*(pc(4,7,2) - pc(4,7,3)))

pc(1,1,1) = preFactorN(1)

pc(2,1,1) = GA(1)*preFactorN(1) - GC(1)*preFactorN(2)

pc(3,1,1) = GA(2)*preFactorN(1) - GC(2)*preFactorN(2)

pc(4,1,1) = GA(3)*preFactorN(1) - GC(3)*preFactorN(2)

pc(5,1,1) = GA(1)*pc(2,1,1) - GC(1)*pc(2,1,2) + inv_2zeta3C*((preFactorN(1) - &
&preFactorN(2)))

pc(6,1,1) = GA(2)*pc(3,1,1) - GC(2)*pc(3,1,2) + inv_2zeta3C*((preFactorN(1) - &
&preFactorN(2)))

pc(7,1,1) = GA(3)*pc(4,1,1) - GC(3)*pc(4,1,2) + inv_2zeta3C*((preFactorN(1) - &
&preFactorN(2)))

pc(8,1,1) = GA(2)*pc(2,1,1) - GC(2)*pc(2,1,2)

pc(9,1,1) = GA(3)*pc(2,1,1) - GC(3)*pc(2,1,2)

pc(10,1,1) = GA(3)*pc(3,1,1) - GC(3)*pc(3,1,2)

pc(11,1,1) = GA(3)*pc(8,1,1) - GC(3)*pc(8,1,2)

pc(12,1,1) = GA(2)*pc(5,1,1) - GC(2)*pc(5,1,2)

pc(13,1,1) = GA(3)*pc(5,1,1) - GC(3)*pc(5,1,2)

pc(14,1,1) = GA(1)*pc(6,1,1) - GC(1)*pc(6,1,2)

pc(15,1,1) = GA(3)*pc(6,1,1) - GC(3)*pc(6,1,2)

pc(16,1,1) = GA(1)*pc(7,1,1) - GC(1)*pc(7,1,2)

pc(17,1,1) = GA(2)*pc(7,1,1) - GC(2)*pc(7,1,2)

pc(18,1,1) = GA(1)*pc(5,1,1) - GC(1)*pc(5,1,2) + inv_2zeta3C*(2*(pc(2,1,1) - p&
&c(2,1,2)))

pc(19,1,1) = GA(2)*pc(6,1,1) - GC(2)*pc(6,1,2) + inv_2zeta3C*(2*(pc(3,1,1) - p&
&c(3,1,2)))

pc(20,1,1) = GA(3)*pc(7,1,1) - GC(3)*pc(7,1,2) + inv_2zeta3C*(2*(pc(4,1,1) - p&
&c(4,1,2)))

pc(1,2,1) = GB(1)*preFactorN(1) - GC(1)*preFactorN(2)

pc(2,2,1) = GA(1)*pc(1,2,1) - GC(1)*pc(1,2,2) + inv_2zeta3C*((preFactorN(1) - &
&preFactorN(2)))

pc(3,2,1) = GA(2)*pc(1,2,1) - GC(2)*pc(1,2,2)

pc(4,2,1) = GA(3)*pc(1,2,1) - GC(3)*pc(1,2,2)

pc(5,2,1) = GA(1)*pc(2,2,1) - GC(1)*pc(2,2,2) + inv_2zeta3C*((pc(1,2,1) - pc(1&
&,2,2)) + (pc(2,1,1) - pc(2,1,2)))

pc(6,2,1) = GA(2)*pc(3,2,1) - GC(2)*pc(3,2,2) + inv_2zeta3C*((pc(1,2,1) - pc(1&
&,2,2)))

pc(7,2,1) = GA(3)*pc(4,2,1) - GC(3)*pc(4,2,2) + inv_2zeta3C*((pc(1,2,1) - pc(1&
&,2,2)))

pc(8,2,1) = GA(2)*pc(2,2,1) - GC(2)*pc(2,2,2)

pc(9,2,1) = GA(3)*pc(2,2,1) - GC(3)*pc(2,2,2)

pc(10,2,1) = GA(3)*pc(3,2,1) - GC(3)*pc(3,2,2)

pc(11,2,1) = GA(3)*pc(8,2,1) - GC(3)*pc(8,2,2)

pc(12,2,1) = GA(2)*pc(5,2,1) - GC(2)*pc(5,2,2)

pc(13,2,1) = GA(3)*pc(5,2,1) - GC(3)*pc(5,2,2)

pc(14,2,1) = GA(1)*pc(6,2,1) - GC(1)*pc(6,2,2) + inv_2zeta3C*((pc(6,1,1) - pc(&
&6,1,2)))

pc(15,2,1) = GA(3)*pc(6,2,1) - GC(3)*pc(6,2,2)

pc(16,2,1) = GA(1)*pc(7,2,1) - GC(1)*pc(7,2,2) + inv_2zeta3C*((pc(7,1,1) - pc(&
&7,1,2)))

pc(17,2,1) = GA(2)*pc(7,2,1) - GC(2)*pc(7,2,2)

pc(18,2,1) = GA(1)*pc(5,2,1) - GC(1)*pc(5,2,2) + inv_2zeta3C*(2*(pc(2,2,1) - p&
&c(2,2,2)) + (pc(5,1,1) - pc(5,1,2)))

pc(19,2,1) = GA(2)*pc(6,2,1) - GC(2)*pc(6,2,2) + inv_2zeta3C*(2*(pc(3,2,1) - p&
&c(3,2,2)))

pc(20,2,1) = GA(3)*pc(7,2,1) - GC(3)*pc(7,2,2) + inv_2zeta3C*(2*(pc(4,2,1) - p&
&c(4,2,2)))

pc(1,3,1) = GB(2)*preFactorN(1) - GC(2)*preFactorN(2)

pc(2,3,1) = GA(1)*pc(1,3,1) - GC(1)*pc(1,3,2)

pc(3,3,1) = GA(2)*pc(1,3,1) - GC(2)*pc(1,3,2) + inv_2zeta3C*((preFactorN(1) - &
&preFactorN(2)))

pc(4,3,1) = GA(3)*pc(1,3,1) - GC(3)*pc(1,3,2)

pc(5,3,1) = GA(1)*pc(2,3,1) - GC(1)*pc(2,3,2) + inv_2zeta3C*((pc(1,3,1) - pc(1&
&,3,2)))

pc(6,3,1) = GA(2)*pc(3,3,1) - GC(2)*pc(3,3,2) + inv_2zeta3C*((pc(1,3,1) - pc(1&
&,3,2)) + (pc(3,1,1) - pc(3,1,2)))

pc(7,3,1) = GA(3)*pc(4,3,1) - GC(3)*pc(4,3,2) + inv_2zeta3C*((pc(1,3,1) - pc(1&
&,3,2)))

pc(8,3,1) = GA(2)*pc(2,3,1) - GC(2)*pc(2,3,2) + inv_2zeta3C*((pc(2,1,1) - pc(2&
&,1,2)))

pc(9,3,1) = GA(3)*pc(2,3,1) - GC(3)*pc(2,3,2)

pc(10,3,1) = GA(3)*pc(3,3,1) - GC(3)*pc(3,3,2)

pc(11,3,1) = GA(3)*pc(8,3,1) - GC(3)*pc(8,3,2)

pc(12,3,1) = GA(2)*pc(5,3,1) - GC(2)*pc(5,3,2) + inv_2zeta3C*((pc(5,1,1) - pc(&
&5,1,2)))

pc(13,3,1) = GA(3)*pc(5,3,1) - GC(3)*pc(5,3,2)

pc(14,3,1) = GA(1)*pc(6,3,1) - GC(1)*pc(6,3,2)

pc(15,3,1) = GA(3)*pc(6,3,1) - GC(3)*pc(6,3,2)

pc(16,3,1) = GA(1)*pc(7,3,1) - GC(1)*pc(7,3,2)

pc(17,3,1) = GA(2)*pc(7,3,1) - GC(2)*pc(7,3,2) + inv_2zeta3C*((pc(7,1,1) - pc(&
&7,1,2)))

pc(18,3,1) = GA(1)*pc(5,3,1) - GC(1)*pc(5,3,2) + inv_2zeta3C*(2*(pc(2,3,1) - p&
&c(2,3,2)))

pc(19,3,1) = GA(2)*pc(6,3,1) - GC(2)*pc(6,3,2) + inv_2zeta3C*(2*(pc(3,3,1) - p&
&c(3,3,2)) + (pc(6,1,1) - pc(6,1,2)))

pc(20,3,1) = GA(3)*pc(7,3,1) - GC(3)*pc(7,3,2) + inv_2zeta3C*(2*(pc(4,3,1) - p&
&c(4,3,2)))

pc(1,4,1) = GB(3)*preFactorN(1) - GC(3)*preFactorN(2)

pc(2,4,1) = GA(1)*pc(1,4,1) - GC(1)*pc(1,4,2)

pc(3,4,1) = GA(2)*pc(1,4,1) - GC(2)*pc(1,4,2)

pc(4,4,1) = GA(3)*pc(1,4,1) - GC(3)*pc(1,4,2) + inv_2zeta3C*((preFactorN(1) - &
&preFactorN(2)))

pc(5,4,1) = GA(1)*pc(2,4,1) - GC(1)*pc(2,4,2) + inv_2zeta3C*((pc(1,4,1) - pc(1&
&,4,2)))

pc(6,4,1) = GA(2)*pc(3,4,1) - GC(2)*pc(3,4,2) + inv_2zeta3C*((pc(1,4,1) - pc(1&
&,4,2)))

pc(7,4,1) = GA(3)*pc(4,4,1) - GC(3)*pc(4,4,2) + inv_2zeta3C*((pc(1,4,1) - pc(1&
&,4,2)) + (pc(4,1,1) - pc(4,1,2)))

pc(8,4,1) = GA(2)*pc(2,4,1) - GC(2)*pc(2,4,2)

pc(9,4,1) = GA(3)*pc(2,4,1) - GC(3)*pc(2,4,2) + inv_2zeta3C*((pc(2,1,1) - pc(2&
&,1,2)))

pc(10,4,1) = GA(3)*pc(3,4,1) - GC(3)*pc(3,4,2) + inv_2zeta3C*((pc(3,1,1) - pc(&
&3,1,2)))

pc(11,4,1) = GA(3)*pc(8,4,1) - GC(3)*pc(8,4,2) + inv_2zeta3C*((pc(8,1,1) - pc(&
&8,1,2)))

pc(12,4,1) = GA(2)*pc(5,4,1) - GC(2)*pc(5,4,2)

pc(13,4,1) = GA(3)*pc(5,4,1) - GC(3)*pc(5,4,2) + inv_2zeta3C*((pc(5,1,1) - pc(&
&5,1,2)))

pc(14,4,1) = GA(1)*pc(6,4,1) - GC(1)*pc(6,4,2)

pc(15,4,1) = GA(3)*pc(6,4,1) - GC(3)*pc(6,4,2) + inv_2zeta3C*((pc(6,1,1) - pc(&
&6,1,2)))

pc(16,4,1) = GA(1)*pc(7,4,1) - GC(1)*pc(7,4,2)

pc(17,4,1) = GA(2)*pc(7,4,1) - GC(2)*pc(7,4,2)

pc(18,4,1) = GA(1)*pc(5,4,1) - GC(1)*pc(5,4,2) + inv_2zeta3C*(2*(pc(2,4,1) - p&
&c(2,4,2)))

pc(19,4,1) = GA(2)*pc(6,4,1) - GC(2)*pc(6,4,2) + inv_2zeta3C*(2*(pc(3,4,1) - p&
&c(3,4,2)))

pc(20,4,1) = GA(3)*pc(7,4,1) - GC(3)*pc(7,4,2) + inv_2zeta3C*(2*(pc(4,4,1) - p&
&c(4,4,2)) + (pc(7,1,1) - pc(7,1,2)))

pc(1,5,1) = GB(1)*pc(1,2,1) - GC(1)*pc(1,2,2) + inv_2zeta3C*((preFactorN(1) - &
&preFactorN(2)))

pc(2,5,1) = GA(1)*pc(1,5,1) - GC(1)*pc(1,5,2) + inv_2zeta3C*(2*(pc(1,2,1) - pc&
&(1,2,2)))

pc(3,5,1) = GA(2)*pc(1,5,1) - GC(2)*pc(1,5,2)

pc(4,5,1) = GA(3)*pc(1,5,1) - GC(3)*pc(1,5,2)

pc(5,5,1) = GA(1)*pc(2,5,1) - GC(1)*pc(2,5,2) + inv_2zeta3C*((pc(1,5,1) - pc(1&
&,5,2)) + 2*(pc(2,2,1) - pc(2,2,2)))

pc(6,5,1) = GA(2)*pc(3,5,1) - GC(2)*pc(3,5,2) + inv_2zeta3C*((pc(1,5,1) - pc(1&
&,5,2)))

pc(7,5,1) = GA(3)*pc(4,5,1) - GC(3)*pc(4,5,2) + inv_2zeta3C*((pc(1,5,1) - pc(1&
&,5,2)))

pc(8,5,1) = GA(2)*pc(2,5,1) - GC(2)*pc(2,5,2)

pc(9,5,1) = GA(3)*pc(2,5,1) - GC(3)*pc(2,5,2)

pc(10,5,1) = GA(3)*pc(3,5,1) - GC(3)*pc(3,5,2)

pc(11,5,1) = GA(3)*pc(8,5,1) - GC(3)*pc(8,5,2)

pc(12,5,1) = GA(2)*pc(5,5,1) - GC(2)*pc(5,5,2)

pc(13,5,1) = GA(3)*pc(5,5,1) - GC(3)*pc(5,5,2)

pc(14,5,1) = GA(1)*pc(6,5,1) - GC(1)*pc(6,5,2) + inv_2zeta3C*(2*(pc(6,2,1) - p&
&c(6,2,2)))

pc(15,5,1) = GA(3)*pc(6,5,1) - GC(3)*pc(6,5,2)

pc(16,5,1) = GA(1)*pc(7,5,1) - GC(1)*pc(7,5,2) + inv_2zeta3C*(2*(pc(7,2,1) - p&
&c(7,2,2)))

pc(17,5,1) = GA(2)*pc(7,5,1) - GC(2)*pc(7,5,2)

pc(18,5,1) = GA(1)*pc(5,5,1) - GC(1)*pc(5,5,2) + inv_2zeta3C*(2*(pc(2,5,1) - p&
&c(2,5,2)) + 2*(pc(5,2,1) - pc(5,2,2)))

pc(19,5,1) = GA(2)*pc(6,5,1) - GC(2)*pc(6,5,2) + inv_2zeta3C*(2*(pc(3,5,1) - p&
&c(3,5,2)))

pc(20,5,1) = GA(3)*pc(7,5,1) - GC(3)*pc(7,5,2) + inv_2zeta3C*(2*(pc(4,5,1) - p&
&c(4,5,2)))

pc(1,6,1) = GB(2)*pc(1,3,1) - GC(2)*pc(1,3,2) + inv_2zeta3C*((preFactorN(1) - &
&preFactorN(2)))

pc(2,6,1) = GA(1)*pc(1,6,1) - GC(1)*pc(1,6,2)

pc(3,6,1) = GA(2)*pc(1,6,1) - GC(2)*pc(1,6,2) + inv_2zeta3C*(2*(pc(1,3,1) - pc&
&(1,3,2)))

pc(4,6,1) = GA(3)*pc(1,6,1) - GC(3)*pc(1,6,2)

pc(5,6,1) = GA(1)*pc(2,6,1) - GC(1)*pc(2,6,2) + inv_2zeta3C*((pc(1,6,1) - pc(1&
&,6,2)))

pc(6,6,1) = GA(2)*pc(3,6,1) - GC(2)*pc(3,6,2) + inv_2zeta3C*((pc(1,6,1) - pc(1&
&,6,2)) + 2*(pc(3,3,1) - pc(3,3,2)))

pc(7,6,1) = GA(3)*pc(4,6,1) - GC(3)*pc(4,6,2) + inv_2zeta3C*((pc(1,6,1) - pc(1&
&,6,2)))

pc(8,6,1) = GA(2)*pc(2,6,1) - GC(2)*pc(2,6,2) + inv_2zeta3C*(2*(pc(2,3,1) - pc&
&(2,3,2)))

pc(9,6,1) = GA(3)*pc(2,6,1) - GC(3)*pc(2,6,2)

pc(10,6,1) = GA(3)*pc(3,6,1) - GC(3)*pc(3,6,2)

pc(11,6,1) = GA(3)*pc(8,6,1) - GC(3)*pc(8,6,2)

pc(12,6,1) = GA(2)*pc(5,6,1) - GC(2)*pc(5,6,2) + inv_2zeta3C*(2*(pc(5,3,1) - p&
&c(5,3,2)))

pc(13,6,1) = GA(3)*pc(5,6,1) - GC(3)*pc(5,6,2)

pc(14,6,1) = GA(1)*pc(6,6,1) - GC(1)*pc(6,6,2)

pc(15,6,1) = GA(3)*pc(6,6,1) - GC(3)*pc(6,6,2)

pc(16,6,1) = GA(1)*pc(7,6,1) - GC(1)*pc(7,6,2)

pc(17,6,1) = GA(2)*pc(7,6,1) - GC(2)*pc(7,6,2) + inv_2zeta3C*(2*(pc(7,3,1) - p&
&c(7,3,2)))

pc(18,6,1) = GA(1)*pc(5,6,1) - GC(1)*pc(5,6,2) + inv_2zeta3C*(2*(pc(2,6,1) - p&
&c(2,6,2)))

pc(19,6,1) = GA(2)*pc(6,6,1) - GC(2)*pc(6,6,2) + inv_2zeta3C*(2*(pc(3,6,1) - p&
&c(3,6,2)) + 2*(pc(6,3,1) - pc(6,3,2)))

pc(20,6,1) = GA(3)*pc(7,6,1) - GC(3)*pc(7,6,2) + inv_2zeta3C*(2*(pc(4,6,1) - p&
&c(4,6,2)))

pc(1,7,1) = GB(3)*pc(1,4,1) - GC(3)*pc(1,4,2) + inv_2zeta3C*((preFactorN(1) - &
&preFactorN(2)))

pc(2,7,1) = GA(1)*pc(1,7,1) - GC(1)*pc(1,7,2)

pc(3,7,1) = GA(2)*pc(1,7,1) - GC(2)*pc(1,7,2)

pc(4,7,1) = GA(3)*pc(1,7,1) - GC(3)*pc(1,7,2) + inv_2zeta3C*(2*(pc(1,4,1) - pc&
&(1,4,2)))

pc(5,7,1) = GA(1)*pc(2,7,1) - GC(1)*pc(2,7,2) + inv_2zeta3C*((pc(1,7,1) - pc(1&
&,7,2)))

pc(6,7,1) = GA(2)*pc(3,7,1) - GC(2)*pc(3,7,2) + inv_2zeta3C*((pc(1,7,1) - pc(1&
&,7,2)))

pc(7,7,1) = GA(3)*pc(4,7,1) - GC(3)*pc(4,7,2) + inv_2zeta3C*((pc(1,7,1) - pc(1&
&,7,2)) + 2*(pc(4,4,1) - pc(4,4,2)))

pc(8,7,1) = GA(2)*pc(2,7,1) - GC(2)*pc(2,7,2)

pc(9,7,1) = GA(3)*pc(2,7,1) - GC(3)*pc(2,7,2) + inv_2zeta3C*(2*(pc(2,4,1) - pc&
&(2,4,2)))

pc(10,7,1) = GA(3)*pc(3,7,1) - GC(3)*pc(3,7,2) + inv_2zeta3C*(2*(pc(3,4,1) - p&
&c(3,4,2)))

pc(11,7,1) = GA(3)*pc(8,7,1) - GC(3)*pc(8,7,2) + inv_2zeta3C*(2*(pc(8,4,1) - p&
&c(8,4,2)))

pc(12,7,1) = GA(2)*pc(5,7,1) - GC(2)*pc(5,7,2)

pc(13,7,1) = GA(3)*pc(5,7,1) - GC(3)*pc(5,7,2) + inv_2zeta3C*(2*(pc(5,4,1) - p&
&c(5,4,2)))

pc(14,7,1) = GA(1)*pc(6,7,1) - GC(1)*pc(6,7,2)

pc(15,7,1) = GA(3)*pc(6,7,1) - GC(3)*pc(6,7,2) + inv_2zeta3C*(2*(pc(6,4,1) - p&
&c(6,4,2)))

pc(16,7,1) = GA(1)*pc(7,7,1) - GC(1)*pc(7,7,2)

pc(17,7,1) = GA(2)*pc(7,7,1) - GC(2)*pc(7,7,2)

pc(18,7,1) = GA(1)*pc(5,7,1) - GC(1)*pc(5,7,2) + inv_2zeta3C*(2*(pc(2,7,1) - p&
&c(2,7,2)))

pc(19,7,1) = GA(2)*pc(6,7,1) - GC(2)*pc(6,7,2) + inv_2zeta3C*(2*(pc(3,7,1) - p&
&c(3,7,2)))

pc(20,7,1) = GA(3)*pc(7,7,1) - GC(3)*pc(7,7,2) + inv_2zeta3C*(2*(pc(4,7,1) - p&
&c(4,7,2)) + 2*(pc(7,4,1) - pc(7,4,2)))

pc(1,8,1) = GB(2)*pc(1,2,1) - GC(2)*pc(1,2,2)

pc(2,8,1) = GA(1)*pc(1,8,1) - GC(1)*pc(1,8,2) + inv_2zeta3C*((pc(1,3,1) - pc(1&
&,3,2)))

pc(3,8,1) = GA(2)*pc(1,8,1) - GC(2)*pc(1,8,2) + inv_2zeta3C*((pc(1,2,1) - pc(1&
&,2,2)))

pc(4,8,1) = GA(3)*pc(1,8,1) - GC(3)*pc(1,8,2)

pc(5,8,1) = GA(1)*pc(2,8,1) - GC(1)*pc(2,8,2) + inv_2zeta3C*((pc(1,8,1) - pc(1&
&,8,2)) + (pc(2,3,1) - pc(2,3,2)))

pc(6,8,1) = GA(2)*pc(3,8,1) - GC(2)*pc(3,8,2) + inv_2zeta3C*((pc(1,8,1) - pc(1&
&,8,2)) + (pc(3,2,1) - pc(3,2,2)))

pc(7,8,1) = GA(3)*pc(4,8,1) - GC(3)*pc(4,8,2) + inv_2zeta3C*((pc(1,8,1) - pc(1&
&,8,2)))

pc(8,8,1) = GA(2)*pc(2,8,1) - GC(2)*pc(2,8,2) + inv_2zeta3C*((pc(2,2,1) - pc(2&
&,2,2)))

pc(9,8,1) = GA(3)*pc(2,8,1) - GC(3)*pc(2,8,2)

pc(10,8,1) = GA(3)*pc(3,8,1) - GC(3)*pc(3,8,2)

pc(11,8,1) = GA(3)*pc(8,8,1) - GC(3)*pc(8,8,2)

pc(12,8,1) = GA(2)*pc(5,8,1) - GC(2)*pc(5,8,2) + inv_2zeta3C*((pc(5,2,1) - pc(&
&5,2,2)))

pc(13,8,1) = GA(3)*pc(5,8,1) - GC(3)*pc(5,8,2)

pc(14,8,1) = GA(1)*pc(6,8,1) - GC(1)*pc(6,8,2) + inv_2zeta3C*((pc(6,3,1) - pc(&
&6,3,2)))

pc(15,8,1) = GA(3)*pc(6,8,1) - GC(3)*pc(6,8,2)

pc(16,8,1) = GA(1)*pc(7,8,1) - GC(1)*pc(7,8,2) + inv_2zeta3C*((pc(7,3,1) - pc(&
&7,3,2)))

pc(17,8,1) = GA(2)*pc(7,8,1) - GC(2)*pc(7,8,2) + inv_2zeta3C*((pc(7,2,1) - pc(&
&7,2,2)))

pc(18,8,1) = GA(1)*pc(5,8,1) - GC(1)*pc(5,8,2) + inv_2zeta3C*(2*(pc(2,8,1) - p&
&c(2,8,2)) + (pc(5,3,1) - pc(5,3,2)))

pc(19,8,1) = GA(2)*pc(6,8,1) - GC(2)*pc(6,8,2) + inv_2zeta3C*(2*(pc(3,8,1) - p&
&c(3,8,2)) + (pc(6,2,1) - pc(6,2,2)))

pc(20,8,1) = GA(3)*pc(7,8,1) - GC(3)*pc(7,8,2) + inv_2zeta3C*(2*(pc(4,8,1) - p&
&c(4,8,2)))

pc(1,9,1) = GB(3)*pc(1,2,1) - GC(3)*pc(1,2,2)

pc(2,9,1) = GA(1)*pc(1,9,1) - GC(1)*pc(1,9,2) + inv_2zeta3C*((pc(1,4,1) - pc(1&
&,4,2)))

pc(3,9,1) = GA(2)*pc(1,9,1) - GC(2)*pc(1,9,2)

pc(4,9,1) = GA(3)*pc(1,9,1) - GC(3)*pc(1,9,2) + inv_2zeta3C*((pc(1,2,1) - pc(1&
&,2,2)))

pc(5,9,1) = GA(1)*pc(2,9,1) - GC(1)*pc(2,9,2) + inv_2zeta3C*((pc(1,9,1) - pc(1&
&,9,2)) + (pc(2,4,1) - pc(2,4,2)))

pc(6,9,1) = GA(2)*pc(3,9,1) - GC(2)*pc(3,9,2) + inv_2zeta3C*((pc(1,9,1) - pc(1&
&,9,2)))

pc(7,9,1) = GA(3)*pc(4,9,1) - GC(3)*pc(4,9,2) + inv_2zeta3C*((pc(1,9,1) - pc(1&
&,9,2)) + (pc(4,2,1) - pc(4,2,2)))

pc(8,9,1) = GA(2)*pc(2,9,1) - GC(2)*pc(2,9,2)

pc(9,9,1) = GA(3)*pc(2,9,1) - GC(3)*pc(2,9,2) + inv_2zeta3C*((pc(2,2,1) - pc(2&
&,2,2)))

pc(10,9,1) = GA(3)*pc(3,9,1) - GC(3)*pc(3,9,2) + inv_2zeta3C*((pc(3,2,1) - pc(&
&3,2,2)))

pc(11,9,1) = GA(3)*pc(8,9,1) - GC(3)*pc(8,9,2) + inv_2zeta3C*((pc(8,2,1) - pc(&
&8,2,2)))

pc(12,9,1) = GA(2)*pc(5,9,1) - GC(2)*pc(5,9,2)

pc(13,9,1) = GA(3)*pc(5,9,1) - GC(3)*pc(5,9,2) + inv_2zeta3C*((pc(5,2,1) - pc(&
&5,2,2)))

pc(14,9,1) = GA(1)*pc(6,9,1) - GC(1)*pc(6,9,2) + inv_2zeta3C*((pc(6,4,1) - pc(&
&6,4,2)))

pc(15,9,1) = GA(3)*pc(6,9,1) - GC(3)*pc(6,9,2) + inv_2zeta3C*((pc(6,2,1) - pc(&
&6,2,2)))

pc(16,9,1) = GA(1)*pc(7,9,1) - GC(1)*pc(7,9,2) + inv_2zeta3C*((pc(7,4,1) - pc(&
&7,4,2)))

pc(17,9,1) = GA(2)*pc(7,9,1) - GC(2)*pc(7,9,2)

pc(18,9,1) = GA(1)*pc(5,9,1) - GC(1)*pc(5,9,2) + inv_2zeta3C*(2*(pc(2,9,1) - p&
&c(2,9,2)) + (pc(5,4,1) - pc(5,4,2)))

pc(19,9,1) = GA(2)*pc(6,9,1) - GC(2)*pc(6,9,2) + inv_2zeta3C*(2*(pc(3,9,1) - p&
&c(3,9,2)))

pc(20,9,1) = GA(3)*pc(7,9,1) - GC(3)*pc(7,9,2) + inv_2zeta3C*(2*(pc(4,9,1) - p&
&c(4,9,2)) + (pc(7,2,1) - pc(7,2,2)))

pc(1,10,1) = GB(3)*pc(1,3,1) - GC(3)*pc(1,3,2)

pc(2,10,1) = GA(1)*pc(1,10,1) - GC(1)*pc(1,10,2)

pc(3,10,1) = GA(2)*pc(1,10,1) - GC(2)*pc(1,10,2) + inv_2zeta3C*((pc(1,4,1) - p&
&c(1,4,2)))

pc(4,10,1) = GA(3)*pc(1,10,1) - GC(3)*pc(1,10,2) + inv_2zeta3C*((pc(1,3,1) - p&
&c(1,3,2)))

pc(5,10,1) = GA(1)*pc(2,10,1) - GC(1)*pc(2,10,2) + inv_2zeta3C*((pc(1,10,1) - &
&pc(1,10,2)))

pc(6,10,1) = GA(2)*pc(3,10,1) - GC(2)*pc(3,10,2) + inv_2zeta3C*((pc(1,10,1) - &
&pc(1,10,2)) + (pc(3,4,1) - pc(3,4,2)))

pc(7,10,1) = GA(3)*pc(4,10,1) - GC(3)*pc(4,10,2) + inv_2zeta3C*((pc(1,10,1) - &
&pc(1,10,2)) + (pc(4,3,1) - pc(4,3,2)))

pc(8,10,1) = GA(2)*pc(2,10,1) - GC(2)*pc(2,10,2) + inv_2zeta3C*((pc(2,4,1) - p&
&c(2,4,2)))

pc(9,10,1) = GA(3)*pc(2,10,1) - GC(3)*pc(2,10,2) + inv_2zeta3C*((pc(2,3,1) - p&
&c(2,3,2)))

pc(10,10,1) = GA(3)*pc(3,10,1) - GC(3)*pc(3,10,2) + inv_2zeta3C*((pc(3,3,1) - &
&pc(3,3,2)))

pc(11,10,1) = GA(3)*pc(8,10,1) - GC(3)*pc(8,10,2) + inv_2zeta3C*((pc(8,3,1) - &
&pc(8,3,2)))

pc(12,10,1) = GA(2)*pc(5,10,1) - GC(2)*pc(5,10,2) + inv_2zeta3C*((pc(5,4,1) - &
&pc(5,4,2)))

pc(13,10,1) = GA(3)*pc(5,10,1) - GC(3)*pc(5,10,2) + inv_2zeta3C*((pc(5,3,1) - &
&pc(5,3,2)))

pc(14,10,1) = GA(1)*pc(6,10,1) - GC(1)*pc(6,10,2)

pc(15,10,1) = GA(3)*pc(6,10,1) - GC(3)*pc(6,10,2) + inv_2zeta3C*((pc(6,3,1) - &
&pc(6,3,2)))

pc(16,10,1) = GA(1)*pc(7,10,1) - GC(1)*pc(7,10,2)

pc(17,10,1) = GA(2)*pc(7,10,1) - GC(2)*pc(7,10,2) + inv_2zeta3C*((pc(7,4,1) - &
&pc(7,4,2)))

pc(18,10,1) = GA(1)*pc(5,10,1) - GC(1)*pc(5,10,2) + inv_2zeta3C*(2*(pc(2,10,1)&
& - pc(2,10,2)))

pc(19,10,1) = GA(2)*pc(6,10,1) - GC(2)*pc(6,10,2) + inv_2zeta3C*(2*(pc(3,10,1)&
& - pc(3,10,2)) + (pc(6,4,1) - pc(6,4,2)))

pc(20,10,1) = GA(3)*pc(7,10,1) - GC(3)*pc(7,10,2) + inv_2zeta3C*(2*(pc(4,10,1)&
& - pc(4,10,2)) + (pc(7,3,1) - pc(7,3,2)))

pc(1,11,1) = GB(3)*pc(1,8,1) - GC(3)*pc(1,8,2)

pc(2,11,1) = GA(1)*pc(1,11,1) - GC(1)*pc(1,11,2) + inv_2zeta3C*((pc(1,10,1) - &
&pc(1,10,2)))

pc(3,11,1) = GA(2)*pc(1,11,1) - GC(2)*pc(1,11,2) + inv_2zeta3C*((pc(1,9,1) - p&
&c(1,9,2)))

pc(4,11,1) = GA(3)*pc(1,11,1) - GC(3)*pc(1,11,2) + inv_2zeta3C*((pc(1,8,1) - p&
&c(1,8,2)))

pc(5,11,1) = GA(1)*pc(2,11,1) - GC(1)*pc(2,11,2) + inv_2zeta3C*((pc(1,11,1) - &
&pc(1,11,2)) + (pc(2,10,1) - pc(2,10,2)))

pc(6,11,1) = GA(2)*pc(3,11,1) - GC(2)*pc(3,11,2) + inv_2zeta3C*((pc(1,11,1) - &
&pc(1,11,2)) + (pc(3,9,1) - pc(3,9,2)))

pc(7,11,1) = GA(3)*pc(4,11,1) - GC(3)*pc(4,11,2) + inv_2zeta3C*((pc(1,11,1) - &
&pc(1,11,2)) + (pc(4,8,1) - pc(4,8,2)))

pc(8,11,1) = GA(2)*pc(2,11,1) - GC(2)*pc(2,11,2) + inv_2zeta3C*((pc(2,9,1) - p&
&c(2,9,2)))

pc(9,11,1) = GA(3)*pc(2,11,1) - GC(3)*pc(2,11,2) + inv_2zeta3C*((pc(2,8,1) - p&
&c(2,8,2)))

pc(10,11,1) = GA(3)*pc(3,11,1) - GC(3)*pc(3,11,2) + inv_2zeta3C*((pc(3,8,1) - &
&pc(3,8,2)))

pc(11,11,1) = GA(3)*pc(8,11,1) - GC(3)*pc(8,11,2) + inv_2zeta3C*((pc(8,8,1) - &
&pc(8,8,2)))

pc(12,11,1) = GA(2)*pc(5,11,1) - GC(2)*pc(5,11,2) + inv_2zeta3C*((pc(5,9,1) - &
&pc(5,9,2)))

pc(13,11,1) = GA(3)*pc(5,11,1) - GC(3)*pc(5,11,2) + inv_2zeta3C*((pc(5,8,1) - &
&pc(5,8,2)))

pc(14,11,1) = GA(1)*pc(6,11,1) - GC(1)*pc(6,11,2) + inv_2zeta3C*((pc(6,10,1) -&
& pc(6,10,2)))

pc(15,11,1) = GA(3)*pc(6,11,1) - GC(3)*pc(6,11,2) + inv_2zeta3C*((pc(6,8,1) - &
&pc(6,8,2)))

pc(16,11,1) = GA(1)*pc(7,11,1) - GC(1)*pc(7,11,2) + inv_2zeta3C*((pc(7,10,1) -&
& pc(7,10,2)))

pc(17,11,1) = GA(2)*pc(7,11,1) - GC(2)*pc(7,11,2) + inv_2zeta3C*((pc(7,9,1) - &
&pc(7,9,2)))

pc(18,11,1) = GA(1)*pc(5,11,1) - GC(1)*pc(5,11,2) + inv_2zeta3C*(2*(pc(2,11,1)&
& - pc(2,11,2)) + (pc(5,10,1) - pc(5,10,2)))

pc(19,11,1) = GA(2)*pc(6,11,1) - GC(2)*pc(6,11,2) + inv_2zeta3C*(2*(pc(3,11,1)&
& - pc(3,11,2)) + (pc(6,9,1) - pc(6,9,2)))

pc(20,11,1) = GA(3)*pc(7,11,1) - GC(3)*pc(7,11,2) + inv_2zeta3C*(2*(pc(4,11,1)&
& - pc(4,11,2)) + (pc(7,8,1) - pc(7,8,2)))

pc(1,12,1) = GB(2)*pc(1,5,1) - GC(2)*pc(1,5,2)

pc(2,12,1) = GA(1)*pc(1,12,1) - GC(1)*pc(1,12,2) + inv_2zeta3C*(2*(pc(1,8,1) -&
& pc(1,8,2)))

pc(3,12,1) = GA(2)*pc(1,12,1) - GC(2)*pc(1,12,2) + inv_2zeta3C*((pc(1,5,1) - p&
&c(1,5,2)))

pc(4,12,1) = GA(3)*pc(1,12,1) - GC(3)*pc(1,12,2)

pc(5,12,1) = GA(1)*pc(2,12,1) - GC(1)*pc(2,12,2) + inv_2zeta3C*((pc(1,12,1) - &
&pc(1,12,2)) + 2*(pc(2,8,1) - pc(2,8,2)))

pc(6,12,1) = GA(2)*pc(3,12,1) - GC(2)*pc(3,12,2) + inv_2zeta3C*((pc(1,12,1) - &
&pc(1,12,2)) + (pc(3,5,1) - pc(3,5,2)))

pc(7,12,1) = GA(3)*pc(4,12,1) - GC(3)*pc(4,12,2) + inv_2zeta3C*((pc(1,12,1) - &
&pc(1,12,2)))

pc(8,12,1) = GA(2)*pc(2,12,1) - GC(2)*pc(2,12,2) + inv_2zeta3C*((pc(2,5,1) - p&
&c(2,5,2)))

pc(9,12,1) = GA(3)*pc(2,12,1) - GC(3)*pc(2,12,2)

pc(10,12,1) = GA(3)*pc(3,12,1) - GC(3)*pc(3,12,2)

pc(11,12,1) = GA(3)*pc(8,12,1) - GC(3)*pc(8,12,2)

pc(12,12,1) = GA(2)*pc(5,12,1) - GC(2)*pc(5,12,2) + inv_2zeta3C*((pc(5,5,1) - &
&pc(5,5,2)))

pc(13,12,1) = GA(3)*pc(5,12,1) - GC(3)*pc(5,12,2)

pc(14,12,1) = GA(1)*pc(6,12,1) - GC(1)*pc(6,12,2) + inv_2zeta3C*(2*(pc(6,8,1) &
&- pc(6,8,2)))

pc(15,12,1) = GA(3)*pc(6,12,1) - GC(3)*pc(6,12,2)

pc(16,12,1) = GA(1)*pc(7,12,1) - GC(1)*pc(7,12,2) + inv_2zeta3C*(2*(pc(7,8,1) &
&- pc(7,8,2)))

pc(17,12,1) = GA(2)*pc(7,12,1) - GC(2)*pc(7,12,2) + inv_2zeta3C*((pc(7,5,1) - &
&pc(7,5,2)))

pc(18,12,1) = GA(1)*pc(5,12,1) - GC(1)*pc(5,12,2) + inv_2zeta3C*(2*(pc(2,12,1)&
& - pc(2,12,2)) + 2*(pc(5,8,1) - pc(5,8,2)))

pc(19,12,1) = GA(2)*pc(6,12,1) - GC(2)*pc(6,12,2) + inv_2zeta3C*(2*(pc(3,12,1)&
& - pc(3,12,2)) + (pc(6,5,1) - pc(6,5,2)))

pc(20,12,1) = GA(3)*pc(7,12,1) - GC(3)*pc(7,12,2) + inv_2zeta3C*(2*(pc(4,12,1)&
& - pc(4,12,2)))

pc(1,13,1) = GB(3)*pc(1,5,1) - GC(3)*pc(1,5,2)

pc(2,13,1) = GA(1)*pc(1,13,1) - GC(1)*pc(1,13,2) + inv_2zeta3C*(2*(pc(1,9,1) -&
& pc(1,9,2)))

pc(3,13,1) = GA(2)*pc(1,13,1) - GC(2)*pc(1,13,2)

pc(4,13,1) = GA(3)*pc(1,13,1) - GC(3)*pc(1,13,2) + inv_2zeta3C*((pc(1,5,1) - p&
&c(1,5,2)))

pc(5,13,1) = GA(1)*pc(2,13,1) - GC(1)*pc(2,13,2) + inv_2zeta3C*((pc(1,13,1) - &
&pc(1,13,2)) + 2*(pc(2,9,1) - pc(2,9,2)))

pc(6,13,1) = GA(2)*pc(3,13,1) - GC(2)*pc(3,13,2) + inv_2zeta3C*((pc(1,13,1) - &
&pc(1,13,2)))

pc(7,13,1) = GA(3)*pc(4,13,1) - GC(3)*pc(4,13,2) + inv_2zeta3C*((pc(1,13,1) - &
&pc(1,13,2)) + (pc(4,5,1) - pc(4,5,2)))

pc(8,13,1) = GA(2)*pc(2,13,1) - GC(2)*pc(2,13,2)

pc(9,13,1) = GA(3)*pc(2,13,1) - GC(3)*pc(2,13,2) + inv_2zeta3C*((pc(2,5,1) - p&
&c(2,5,2)))

pc(10,13,1) = GA(3)*pc(3,13,1) - GC(3)*pc(3,13,2) + inv_2zeta3C*((pc(3,5,1) - &
&pc(3,5,2)))

pc(11,13,1) = GA(3)*pc(8,13,1) - GC(3)*pc(8,13,2) + inv_2zeta3C*((pc(8,5,1) - &
&pc(8,5,2)))

pc(12,13,1) = GA(2)*pc(5,13,1) - GC(2)*pc(5,13,2)

pc(13,13,1) = GA(3)*pc(5,13,1) - GC(3)*pc(5,13,2) + inv_2zeta3C*((pc(5,5,1) - &
&pc(5,5,2)))

pc(14,13,1) = GA(1)*pc(6,13,1) - GC(1)*pc(6,13,2) + inv_2zeta3C*(2*(pc(6,9,1) &
&- pc(6,9,2)))

pc(15,13,1) = GA(3)*pc(6,13,1) - GC(3)*pc(6,13,2) + inv_2zeta3C*((pc(6,5,1) - &
&pc(6,5,2)))

pc(16,13,1) = GA(1)*pc(7,13,1) - GC(1)*pc(7,13,2) + inv_2zeta3C*(2*(pc(7,9,1) &
&- pc(7,9,2)))

pc(17,13,1) = GA(2)*pc(7,13,1) - GC(2)*pc(7,13,2)

pc(18,13,1) = GA(1)*pc(5,13,1) - GC(1)*pc(5,13,2) + inv_2zeta3C*(2*(pc(2,13,1)&
& - pc(2,13,2)) + 2*(pc(5,9,1) - pc(5,9,2)))

pc(19,13,1) = GA(2)*pc(6,13,1) - GC(2)*pc(6,13,2) + inv_2zeta3C*(2*(pc(3,13,1)&
& - pc(3,13,2)))

pc(20,13,1) = GA(3)*pc(7,13,1) - GC(3)*pc(7,13,2) + inv_2zeta3C*(2*(pc(4,13,1)&
& - pc(4,13,2)) + (pc(7,5,1) - pc(7,5,2)))

pc(1,14,1) = GB(1)*pc(1,6,1) - GC(1)*pc(1,6,2)

pc(2,14,1) = GA(1)*pc(1,14,1) - GC(1)*pc(1,14,2) + inv_2zeta3C*((pc(1,6,1) - p&
&c(1,6,2)))

pc(3,14,1) = GA(2)*pc(1,14,1) - GC(2)*pc(1,14,2) + inv_2zeta3C*(2*(pc(1,8,1) -&
& pc(1,8,2)))

pc(4,14,1) = GA(3)*pc(1,14,1) - GC(3)*pc(1,14,2)

pc(5,14,1) = GA(1)*pc(2,14,1) - GC(1)*pc(2,14,2) + inv_2zeta3C*((pc(1,14,1) - &
&pc(1,14,2)) + (pc(2,6,1) - pc(2,6,2)))

pc(6,14,1) = GA(2)*pc(3,14,1) - GC(2)*pc(3,14,2) + inv_2zeta3C*((pc(1,14,1) - &
&pc(1,14,2)) + 2*(pc(3,8,1) - pc(3,8,2)))

pc(7,14,1) = GA(3)*pc(4,14,1) - GC(3)*pc(4,14,2) + inv_2zeta3C*((pc(1,14,1) - &
&pc(1,14,2)))

pc(8,14,1) = GA(2)*pc(2,14,1) - GC(2)*pc(2,14,2) + inv_2zeta3C*(2*(pc(2,8,1) -&
& pc(2,8,2)))

pc(9,14,1) = GA(3)*pc(2,14,1) - GC(3)*pc(2,14,2)

pc(10,14,1) = GA(3)*pc(3,14,1) - GC(3)*pc(3,14,2)

pc(11,14,1) = GA(3)*pc(8,14,1) - GC(3)*pc(8,14,2)

pc(12,14,1) = GA(2)*pc(5,14,1) - GC(2)*pc(5,14,2) + inv_2zeta3C*(2*(pc(5,8,1) &
&- pc(5,8,2)))

pc(13,14,1) = GA(3)*pc(5,14,1) - GC(3)*pc(5,14,2)

pc(14,14,1) = GA(1)*pc(6,14,1) - GC(1)*pc(6,14,2) + inv_2zeta3C*((pc(6,6,1) - &
&pc(6,6,2)))

pc(15,14,1) = GA(3)*pc(6,14,1) - GC(3)*pc(6,14,2)

pc(16,14,1) = GA(1)*pc(7,14,1) - GC(1)*pc(7,14,2) + inv_2zeta3C*((pc(7,6,1) - &
&pc(7,6,2)))

pc(17,14,1) = GA(2)*pc(7,14,1) - GC(2)*pc(7,14,2) + inv_2zeta3C*(2*(pc(7,8,1) &
&- pc(7,8,2)))

pc(18,14,1) = GA(1)*pc(5,14,1) - GC(1)*pc(5,14,2) + inv_2zeta3C*(2*(pc(2,14,1)&
& - pc(2,14,2)) + (pc(5,6,1) - pc(5,6,2)))

pc(19,14,1) = GA(2)*pc(6,14,1) - GC(2)*pc(6,14,2) + inv_2zeta3C*(2*(pc(3,14,1)&
& - pc(3,14,2)) + 2*(pc(6,8,1) - pc(6,8,2)))

pc(20,14,1) = GA(3)*pc(7,14,1) - GC(3)*pc(7,14,2) + inv_2zeta3C*(2*(pc(4,14,1)&
& - pc(4,14,2)))

pc(1,15,1) = GB(3)*pc(1,6,1) - GC(3)*pc(1,6,2)

pc(2,15,1) = GA(1)*pc(1,15,1) - GC(1)*pc(1,15,2)

pc(3,15,1) = GA(2)*pc(1,15,1) - GC(2)*pc(1,15,2) + inv_2zeta3C*(2*(pc(1,10,1) &
&- pc(1,10,2)))

pc(4,15,1) = GA(3)*pc(1,15,1) - GC(3)*pc(1,15,2) + inv_2zeta3C*((pc(1,6,1) - p&
&c(1,6,2)))

pc(5,15,1) = GA(1)*pc(2,15,1) - GC(1)*pc(2,15,2) + inv_2zeta3C*((pc(1,15,1) - &
&pc(1,15,2)))

pc(6,15,1) = GA(2)*pc(3,15,1) - GC(2)*pc(3,15,2) + inv_2zeta3C*((pc(1,15,1) - &
&pc(1,15,2)) + 2*(pc(3,10,1) - pc(3,10,2)))

pc(7,15,1) = GA(3)*pc(4,15,1) - GC(3)*pc(4,15,2) + inv_2zeta3C*((pc(1,15,1) - &
&pc(1,15,2)) + (pc(4,6,1) - pc(4,6,2)))

pc(8,15,1) = GA(2)*pc(2,15,1) - GC(2)*pc(2,15,2) + inv_2zeta3C*(2*(pc(2,10,1) &
&- pc(2,10,2)))

pc(9,15,1) = GA(3)*pc(2,15,1) - GC(3)*pc(2,15,2) + inv_2zeta3C*((pc(2,6,1) - p&
&c(2,6,2)))

pc(10,15,1) = GA(3)*pc(3,15,1) - GC(3)*pc(3,15,2) + inv_2zeta3C*((pc(3,6,1) - &
&pc(3,6,2)))

pc(11,15,1) = GA(3)*pc(8,15,1) - GC(3)*pc(8,15,2) + inv_2zeta3C*((pc(8,6,1) - &
&pc(8,6,2)))

pc(12,15,1) = GA(2)*pc(5,15,1) - GC(2)*pc(5,15,2) + inv_2zeta3C*(2*(pc(5,10,1)&
& - pc(5,10,2)))

pc(13,15,1) = GA(3)*pc(5,15,1) - GC(3)*pc(5,15,2) + inv_2zeta3C*((pc(5,6,1) - &
&pc(5,6,2)))

pc(14,15,1) = GA(1)*pc(6,15,1) - GC(1)*pc(6,15,2)

pc(15,15,1) = GA(3)*pc(6,15,1) - GC(3)*pc(6,15,2) + inv_2zeta3C*((pc(6,6,1) - &
&pc(6,6,2)))

pc(16,15,1) = GA(1)*pc(7,15,1) - GC(1)*pc(7,15,2)

pc(17,15,1) = GA(2)*pc(7,15,1) - GC(2)*pc(7,15,2) + inv_2zeta3C*(2*(pc(7,10,1)&
& - pc(7,10,2)))

pc(18,15,1) = GA(1)*pc(5,15,1) - GC(1)*pc(5,15,2) + inv_2zeta3C*(2*(pc(2,15,1)&
& - pc(2,15,2)))

pc(19,15,1) = GA(2)*pc(6,15,1) - GC(2)*pc(6,15,2) + inv_2zeta3C*(2*(pc(3,15,1)&
& - pc(3,15,2)) + 2*(pc(6,10,1) - pc(6,10,2)))

pc(20,15,1) = GA(3)*pc(7,15,1) - GC(3)*pc(7,15,2) + inv_2zeta3C*(2*(pc(4,15,1)&
& - pc(4,15,2)) + (pc(7,6,1) - pc(7,6,2)))

pc(1,16,1) = GB(1)*pc(1,7,1) - GC(1)*pc(1,7,2)

pc(2,16,1) = GA(1)*pc(1,16,1) - GC(1)*pc(1,16,2) + inv_2zeta3C*((pc(1,7,1) - p&
&c(1,7,2)))

pc(3,16,1) = GA(2)*pc(1,16,1) - GC(2)*pc(1,16,2)

pc(4,16,1) = GA(3)*pc(1,16,1) - GC(3)*pc(1,16,2) + inv_2zeta3C*(2*(pc(1,9,1) -&
& pc(1,9,2)))

pc(5,16,1) = GA(1)*pc(2,16,1) - GC(1)*pc(2,16,2) + inv_2zeta3C*((pc(1,16,1) - &
&pc(1,16,2)) + (pc(2,7,1) - pc(2,7,2)))

pc(6,16,1) = GA(2)*pc(3,16,1) - GC(2)*pc(3,16,2) + inv_2zeta3C*((pc(1,16,1) - &
&pc(1,16,2)))

pc(7,16,1) = GA(3)*pc(4,16,1) - GC(3)*pc(4,16,2) + inv_2zeta3C*((pc(1,16,1) - &
&pc(1,16,2)) + 2*(pc(4,9,1) - pc(4,9,2)))

pc(8,16,1) = GA(2)*pc(2,16,1) - GC(2)*pc(2,16,2)

pc(9,16,1) = GA(3)*pc(2,16,1) - GC(3)*pc(2,16,2) + inv_2zeta3C*(2*(pc(2,9,1) -&
& pc(2,9,2)))

pc(10,16,1) = GA(3)*pc(3,16,1) - GC(3)*pc(3,16,2) + inv_2zeta3C*(2*(pc(3,9,1) &
&- pc(3,9,2)))

pc(11,16,1) = GA(3)*pc(8,16,1) - GC(3)*pc(8,16,2) + inv_2zeta3C*(2*(pc(8,9,1) &
&- pc(8,9,2)))

pc(12,16,1) = GA(2)*pc(5,16,1) - GC(2)*pc(5,16,2)

pc(13,16,1) = GA(3)*pc(5,16,1) - GC(3)*pc(5,16,2) + inv_2zeta3C*(2*(pc(5,9,1) &
&- pc(5,9,2)))

pc(14,16,1) = GA(1)*pc(6,16,1) - GC(1)*pc(6,16,2) + inv_2zeta3C*((pc(6,7,1) - &
&pc(6,7,2)))

pc(15,16,1) = GA(3)*pc(6,16,1) - GC(3)*pc(6,16,2) + inv_2zeta3C*(2*(pc(6,9,1) &
&- pc(6,9,2)))

pc(16,16,1) = GA(1)*pc(7,16,1) - GC(1)*pc(7,16,2) + inv_2zeta3C*((pc(7,7,1) - &
&pc(7,7,2)))

pc(17,16,1) = GA(2)*pc(7,16,1) - GC(2)*pc(7,16,2)

pc(18,16,1) = GA(1)*pc(5,16,1) - GC(1)*pc(5,16,2) + inv_2zeta3C*(2*(pc(2,16,1)&
& - pc(2,16,2)) + (pc(5,7,1) - pc(5,7,2)))

pc(19,16,1) = GA(2)*pc(6,16,1) - GC(2)*pc(6,16,2) + inv_2zeta3C*(2*(pc(3,16,1)&
& - pc(3,16,2)))

pc(20,16,1) = GA(3)*pc(7,16,1) - GC(3)*pc(7,16,2) + inv_2zeta3C*(2*(pc(4,16,1)&
& - pc(4,16,2)) + 2*(pc(7,9,1) - pc(7,9,2)))

pc(1,17,1) = GB(2)*pc(1,7,1) - GC(2)*pc(1,7,2)

pc(2,17,1) = GA(1)*pc(1,17,1) - GC(1)*pc(1,17,2)

pc(3,17,1) = GA(2)*pc(1,17,1) - GC(2)*pc(1,17,2) + inv_2zeta3C*((pc(1,7,1) - p&
&c(1,7,2)))

pc(4,17,1) = GA(3)*pc(1,17,1) - GC(3)*pc(1,17,2) + inv_2zeta3C*(2*(pc(1,10,1) &
&- pc(1,10,2)))

pc(5,17,1) = GA(1)*pc(2,17,1) - GC(1)*pc(2,17,2) + inv_2zeta3C*((pc(1,17,1) - &
&pc(1,17,2)))

pc(6,17,1) = GA(2)*pc(3,17,1) - GC(2)*pc(3,17,2) + inv_2zeta3C*((pc(1,17,1) - &
&pc(1,17,2)) + (pc(3,7,1) - pc(3,7,2)))

pc(7,17,1) = GA(3)*pc(4,17,1) - GC(3)*pc(4,17,2) + inv_2zeta3C*((pc(1,17,1) - &
&pc(1,17,2)) + 2*(pc(4,10,1) - pc(4,10,2)))

pc(8,17,1) = GA(2)*pc(2,17,1) - GC(2)*pc(2,17,2) + inv_2zeta3C*((pc(2,7,1) - p&
&c(2,7,2)))

pc(9,17,1) = GA(3)*pc(2,17,1) - GC(3)*pc(2,17,2) + inv_2zeta3C*(2*(pc(2,10,1) &
&- pc(2,10,2)))

pc(10,17,1) = GA(3)*pc(3,17,1) - GC(3)*pc(3,17,2) + inv_2zeta3C*(2*(pc(3,10,1)&
& - pc(3,10,2)))

pc(11,17,1) = GA(3)*pc(8,17,1) - GC(3)*pc(8,17,2) + inv_2zeta3C*(2*(pc(8,10,1)&
& - pc(8,10,2)))

pc(12,17,1) = GA(2)*pc(5,17,1) - GC(2)*pc(5,17,2) + inv_2zeta3C*((pc(5,7,1) - &
&pc(5,7,2)))

pc(13,17,1) = GA(3)*pc(5,17,1) - GC(3)*pc(5,17,2) + inv_2zeta3C*(2*(pc(5,10,1)&
& - pc(5,10,2)))

pc(14,17,1) = GA(1)*pc(6,17,1) - GC(1)*pc(6,17,2)

pc(15,17,1) = GA(3)*pc(6,17,1) - GC(3)*pc(6,17,2) + inv_2zeta3C*(2*(pc(6,10,1)&
& - pc(6,10,2)))

pc(16,17,1) = GA(1)*pc(7,17,1) - GC(1)*pc(7,17,2)

pc(17,17,1) = GA(2)*pc(7,17,1) - GC(2)*pc(7,17,2) + inv_2zeta3C*((pc(7,7,1) - &
&pc(7,7,2)))

pc(18,17,1) = GA(1)*pc(5,17,1) - GC(1)*pc(5,17,2) + inv_2zeta3C*(2*(pc(2,17,1)&
& - pc(2,17,2)))

pc(19,17,1) = GA(2)*pc(6,17,1) - GC(2)*pc(6,17,2) + inv_2zeta3C*(2*(pc(3,17,1)&
& - pc(3,17,2)) + (pc(6,7,1) - pc(6,7,2)))

pc(20,17,1) = GA(3)*pc(7,17,1) - GC(3)*pc(7,17,2) + inv_2zeta3C*(2*(pc(4,17,1)&
& - pc(4,17,2)) + 2*(pc(7,10,1) - pc(7,10,2)))

pc(1,18,1) = GB(1)*pc(1,5,1) - GC(1)*pc(1,5,2) + inv_2zeta3C*(2*(pc(1,2,1) - p&
&c(1,2,2)))

pc(2,18,1) = GA(1)*pc(1,18,1) - GC(1)*pc(1,18,2) + inv_2zeta3C*(3*(pc(1,5,1) -&
& pc(1,5,2)))

pc(3,18,1) = GA(2)*pc(1,18,1) - GC(2)*pc(1,18,2)

pc(4,18,1) = GA(3)*pc(1,18,1) - GC(3)*pc(1,18,2)

pc(5,18,1) = GA(1)*pc(2,18,1) - GC(1)*pc(2,18,2) + inv_2zeta3C*((pc(1,18,1) - &
&pc(1,18,2)) + 3*(pc(2,5,1) - pc(2,5,2)))

pc(6,18,1) = GA(2)*pc(3,18,1) - GC(2)*pc(3,18,2) + inv_2zeta3C*((pc(1,18,1) - &
&pc(1,18,2)))

pc(7,18,1) = GA(3)*pc(4,18,1) - GC(3)*pc(4,18,2) + inv_2zeta3C*((pc(1,18,1) - &
&pc(1,18,2)))

pc(8,18,1) = GA(2)*pc(2,18,1) - GC(2)*pc(2,18,2)

pc(9,18,1) = GA(3)*pc(2,18,1) - GC(3)*pc(2,18,2)

pc(10,18,1) = GA(3)*pc(3,18,1) - GC(3)*pc(3,18,2)

pc(11,18,1) = GA(3)*pc(8,18,1) - GC(3)*pc(8,18,2)

pc(12,18,1) = GA(2)*pc(5,18,1) - GC(2)*pc(5,18,2)

pc(13,18,1) = GA(3)*pc(5,18,1) - GC(3)*pc(5,18,2)

pc(14,18,1) = GA(1)*pc(6,18,1) - GC(1)*pc(6,18,2) + inv_2zeta3C*(3*(pc(6,5,1) &
&- pc(6,5,2)))

pc(15,18,1) = GA(3)*pc(6,18,1) - GC(3)*pc(6,18,2)

pc(16,18,1) = GA(1)*pc(7,18,1) - GC(1)*pc(7,18,2) + inv_2zeta3C*(3*(pc(7,5,1) &
&- pc(7,5,2)))

pc(17,18,1) = GA(2)*pc(7,18,1) - GC(2)*pc(7,18,2)

pc(18,18,1) = GA(1)*pc(5,18,1) - GC(1)*pc(5,18,2) + inv_2zeta3C*(2*(pc(2,18,1)&
& - pc(2,18,2)) + 3*(pc(5,5,1) - pc(5,5,2)))

pc(19,18,1) = GA(2)*pc(6,18,1) - GC(2)*pc(6,18,2) + inv_2zeta3C*(2*(pc(3,18,1)&
& - pc(3,18,2)))

pc(20,18,1) = GA(3)*pc(7,18,1) - GC(3)*pc(7,18,2) + inv_2zeta3C*(2*(pc(4,18,1)&
& - pc(4,18,2)))

pc(1,19,1) = GB(2)*pc(1,6,1) - GC(2)*pc(1,6,2) + inv_2zeta3C*(2*(pc(1,3,1) - p&
&c(1,3,2)))

pc(2,19,1) = GA(1)*pc(1,19,1) - GC(1)*pc(1,19,2)

pc(3,19,1) = GA(2)*pc(1,19,1) - GC(2)*pc(1,19,2) + inv_2zeta3C*(3*(pc(1,6,1) -&
& pc(1,6,2)))

pc(4,19,1) = GA(3)*pc(1,19,1) - GC(3)*pc(1,19,2)

pc(5,19,1) = GA(1)*pc(2,19,1) - GC(1)*pc(2,19,2) + inv_2zeta3C*((pc(1,19,1) - &
&pc(1,19,2)))

pc(6,19,1) = GA(2)*pc(3,19,1) - GC(2)*pc(3,19,2) + inv_2zeta3C*((pc(1,19,1) - &
&pc(1,19,2)) + 3*(pc(3,6,1) - pc(3,6,2)))

pc(7,19,1) = GA(3)*pc(4,19,1) - GC(3)*pc(4,19,2) + inv_2zeta3C*((pc(1,19,1) - &
&pc(1,19,2)))

pc(8,19,1) = GA(2)*pc(2,19,1) - GC(2)*pc(2,19,2) + inv_2zeta3C*(3*(pc(2,6,1) -&
& pc(2,6,2)))

pc(9,19,1) = GA(3)*pc(2,19,1) - GC(3)*pc(2,19,2)

pc(10,19,1) = GA(3)*pc(3,19,1) - GC(3)*pc(3,19,2)

pc(11,19,1) = GA(3)*pc(8,19,1) - GC(3)*pc(8,19,2)

pc(12,19,1) = GA(2)*pc(5,19,1) - GC(2)*pc(5,19,2) + inv_2zeta3C*(3*(pc(5,6,1) &
&- pc(5,6,2)))

pc(13,19,1) = GA(3)*pc(5,19,1) - GC(3)*pc(5,19,2)

pc(14,19,1) = GA(1)*pc(6,19,1) - GC(1)*pc(6,19,2)

pc(15,19,1) = GA(3)*pc(6,19,1) - GC(3)*pc(6,19,2)

pc(16,19,1) = GA(1)*pc(7,19,1) - GC(1)*pc(7,19,2)

pc(17,19,1) = GA(2)*pc(7,19,1) - GC(2)*pc(7,19,2) + inv_2zeta3C*(3*(pc(7,6,1) &
&- pc(7,6,2)))

pc(18,19,1) = GA(1)*pc(5,19,1) - GC(1)*pc(5,19,2) + inv_2zeta3C*(2*(pc(2,19,1)&
& - pc(2,19,2)))

pc(19,19,1) = GA(2)*pc(6,19,1) - GC(2)*pc(6,19,2) + inv_2zeta3C*(2*(pc(3,19,1)&
& - pc(3,19,2)) + 3*(pc(6,6,1) - pc(6,6,2)))

pc(20,19,1) = GA(3)*pc(7,19,1) - GC(3)*pc(7,19,2) + inv_2zeta3C*(2*(pc(4,19,1)&
& - pc(4,19,2)))

pc(1,20,1) = GB(3)*pc(1,7,1) - GC(3)*pc(1,7,2) + inv_2zeta3C*(2*(pc(1,4,1) - p&
&c(1,4,2)))

pc(2,20,1) = GA(1)*pc(1,20,1) - GC(1)*pc(1,20,2)

pc(3,20,1) = GA(2)*pc(1,20,1) - GC(2)*pc(1,20,2)

pc(4,20,1) = GA(3)*pc(1,20,1) - GC(3)*pc(1,20,2) + inv_2zeta3C*(3*(pc(1,7,1) -&
& pc(1,7,2)))

pc(5,20,1) = GA(1)*pc(2,20,1) - GC(1)*pc(2,20,2) + inv_2zeta3C*((pc(1,20,1) - &
&pc(1,20,2)))

pc(6,20,1) = GA(2)*pc(3,20,1) - GC(2)*pc(3,20,2) + inv_2zeta3C*((pc(1,20,1) - &
&pc(1,20,2)))

pc(7,20,1) = GA(3)*pc(4,20,1) - GC(3)*pc(4,20,2) + inv_2zeta3C*((pc(1,20,1) - &
&pc(1,20,2)) + 3*(pc(4,7,1) - pc(4,7,2)))

pc(8,20,1) = GA(2)*pc(2,20,1) - GC(2)*pc(2,20,2)

pc(9,20,1) = GA(3)*pc(2,20,1) - GC(3)*pc(2,20,2) + inv_2zeta3C*(3*(pc(2,7,1) -&
& pc(2,7,2)))

pc(10,20,1) = GA(3)*pc(3,20,1) - GC(3)*pc(3,20,2) + inv_2zeta3C*(3*(pc(3,7,1) &
&- pc(3,7,2)))

pc(11,20,1) = GA(3)*pc(8,20,1) - GC(3)*pc(8,20,2) + inv_2zeta3C*(3*(pc(8,7,1) &
&- pc(8,7,2)))

pc(12,20,1) = GA(2)*pc(5,20,1) - GC(2)*pc(5,20,2)

pc(13,20,1) = GA(3)*pc(5,20,1) - GC(3)*pc(5,20,2) + inv_2zeta3C*(3*(pc(5,7,1) &
&- pc(5,7,2)))

pc(14,20,1) = GA(1)*pc(6,20,1) - GC(1)*pc(6,20,2)

pc(15,20,1) = GA(3)*pc(6,20,1) - GC(3)*pc(6,20,2) + inv_2zeta3C*(3*(pc(6,7,1) &
&- pc(6,7,2)))

pc(16,20,1) = GA(1)*pc(7,20,1) - GC(1)*pc(7,20,2)

pc(17,20,1) = GA(2)*pc(7,20,1) - GC(2)*pc(7,20,2)

pc(18,20,1) = GA(1)*pc(5,20,1) - GC(1)*pc(5,20,2) + inv_2zeta3C*(2*(pc(2,20,1)&
& - pc(2,20,2)))

pc(19,20,1) = GA(2)*pc(6,20,1) - GC(2)*pc(6,20,2) + inv_2zeta3C*(2*(pc(3,20,1)&
& - pc(3,20,2)))

pc(20,20,1) = GA(3)*pc(7,20,1) - GC(3)*pc(7,20,2) + inv_2zeta3C*(2*(pc(4,20,1)&
& - pc(4,20,2)) + 3*(pc(7,7,1) - pc(7,7,2)))

sh(1,1) = pc(1,1,1)

sh(2,1) = pc(2,1,1)

sh(3,1) = pc(3,1,1)

sh(4,1) = pc(4,1,1)

sh(5,1) = pc(8,1,1)

sh(6,1) = pc(9,1,1)

sh(7,1) = pc(10,1,1)

sh(8,1) = pc(5,1,1) - pc(6,1,1)

sh(9,1) = 2*pc(7,1,1) - pc(5,1,1) - pc(6,1,1)

sh(10,1) = pc(11,1,1)

sh(11,1) = pc(13,1,1) - pc(15,1,1)

sh(12,1) = pc(18,1,1) - 3*pc(14,1,1)

sh(13,1) = 3*pc(12,1,1) - pc(19,1,1)

sh(14,1) = 2*pc(20,1,1) - 3*pc(13,1,1) - 3*pc(15,1,1)

sh(15,1) = 4*pc(16,1,1) - pc(18,1,1) - pc(14,1,1)

sh(16,1) = 4*pc(17,1,1) - pc(12,1,1) - pc(19,1,1)

sh(1,2) = pc(1,2,1)

sh(2,2) = pc(2,2,1)

sh(3,2) = pc(3,2,1)

sh(4,2) = pc(4,2,1)

sh(5,2) = pc(8,2,1)

sh(6,2) = pc(9,2,1)

sh(7,2) = pc(10,2,1)

sh(8,2) = pc(5,2,1) - pc(6,2,1)

sh(9,2) = 2*pc(7,2,1) - pc(5,2,1) - pc(6,2,1)

sh(10,2) = pc(11,2,1)

sh(11,2) = pc(13,2,1) - pc(15,2,1)

sh(12,2) = pc(18,2,1) - 3*pc(14,2,1)

sh(13,2) = 3*pc(12,2,1) - pc(19,2,1)

sh(14,2) = 2*pc(20,2,1) - 3*pc(13,2,1) - 3*pc(15,2,1)

sh(15,2) = 4*pc(16,2,1) - pc(18,2,1) - pc(14,2,1)

sh(16,2) = 4*pc(17,2,1) - pc(12,2,1) - pc(19,2,1)

sh(1,3) = pc(1,3,1)

sh(2,3) = pc(2,3,1)

sh(3,3) = pc(3,3,1)

sh(4,3) = pc(4,3,1)

sh(5,3) = pc(8,3,1)

sh(6,3) = pc(9,3,1)

sh(7,3) = pc(10,3,1)

sh(8,3) = pc(5,3,1) - pc(6,3,1)

sh(9,3) = 2*pc(7,3,1) - pc(5,3,1) - pc(6,3,1)

sh(10,3) = pc(11,3,1)

sh(11,3) = pc(13,3,1) - pc(15,3,1)

sh(12,3) = pc(18,3,1) - 3*pc(14,3,1)

sh(13,3) = 3*pc(12,3,1) - pc(19,3,1)

sh(14,3) = 2*pc(20,3,1) - 3*pc(13,3,1) - 3*pc(15,3,1)

sh(15,3) = 4*pc(16,3,1) - pc(18,3,1) - pc(14,3,1)

sh(16,3) = 4*pc(17,3,1) - pc(12,3,1) - pc(19,3,1)

sh(1,4) = pc(1,4,1)

sh(2,4) = pc(2,4,1)

sh(3,4) = pc(3,4,1)

sh(4,4) = pc(4,4,1)

sh(5,4) = pc(8,4,1)

sh(6,4) = pc(9,4,1)

sh(7,4) = pc(10,4,1)

sh(8,4) = pc(5,4,1) - pc(6,4,1)

sh(9,4) = 2*pc(7,4,1) - pc(5,4,1) - pc(6,4,1)

sh(10,4) = pc(11,4,1)

sh(11,4) = pc(13,4,1) - pc(15,4,1)

sh(12,4) = pc(18,4,1) - 3*pc(14,4,1)

sh(13,4) = 3*pc(12,4,1) - pc(19,4,1)

sh(14,4) = 2*pc(20,4,1) - 3*pc(13,4,1) - 3*pc(15,4,1)

sh(15,4) = 4*pc(16,4,1) - pc(18,4,1) - pc(14,4,1)

sh(16,4) = 4*pc(17,4,1) - pc(12,4,1) - pc(19,4,1)

sh(1,5) = pc(1,8,1)

sh(2,5) = pc(2,8,1)

sh(3,5) = pc(3,8,1)

sh(4,5) = pc(4,8,1)

sh(5,5) = pc(8,8,1)

sh(6,5) = pc(9,8,1)

sh(7,5) = pc(10,8,1)

sh(8,5) = pc(5,8,1) - pc(6,8,1)

sh(9,5) = 2*pc(7,8,1) - pc(5,8,1) - pc(6,8,1)

sh(10,5) = pc(11,8,1)

sh(11,5) = pc(13,8,1) - pc(15,8,1)

sh(12,5) = pc(18,8,1) - 3*pc(14,8,1)

sh(13,5) = 3*pc(12,8,1) - pc(19,8,1)

sh(14,5) = 2*pc(20,8,1) - 3*pc(13,8,1) - 3*pc(15,8,1)

sh(15,5) = 4*pc(16,8,1) - pc(18,8,1) - pc(14,8,1)

sh(16,5) = 4*pc(17,8,1) - pc(12,8,1) - pc(19,8,1)

sh(1,6) = pc(1,9,1)

sh(2,6) = pc(2,9,1)

sh(3,6) = pc(3,9,1)

sh(4,6) = pc(4,9,1)

sh(5,6) = pc(8,9,1)

sh(6,6) = pc(9,9,1)

sh(7,6) = pc(10,9,1)

sh(8,6) = pc(5,9,1) - pc(6,9,1)

sh(9,6) = 2*pc(7,9,1) - pc(5,9,1) - pc(6,9,1)

sh(10,6) = pc(11,9,1)

sh(11,6) = pc(13,9,1) - pc(15,9,1)

sh(12,6) = pc(18,9,1) - 3*pc(14,9,1)

sh(13,6) = 3*pc(12,9,1) - pc(19,9,1)

sh(14,6) = 2*pc(20,9,1) - 3*pc(13,9,1) - 3*pc(15,9,1)

sh(15,6) = 4*pc(16,9,1) - pc(18,9,1) - pc(14,9,1)

sh(16,6) = 4*pc(17,9,1) - pc(12,9,1) - pc(19,9,1)

sh(1,7) = pc(1,10,1)

sh(2,7) = pc(2,10,1)

sh(3,7) = pc(3,10,1)

sh(4,7) = pc(4,10,1)

sh(5,7) = pc(8,10,1)

sh(6,7) = pc(9,10,1)

sh(7,7) = pc(10,10,1)

sh(8,7) = pc(5,10,1) - pc(6,10,1)

sh(9,7) = 2*pc(7,10,1) - pc(5,10,1) - pc(6,10,1)

sh(10,7) = pc(11,10,1)

sh(11,7) = pc(13,10,1) - pc(15,10,1)

sh(12,7) = pc(18,10,1) - 3*pc(14,10,1)

sh(13,7) = 3*pc(12,10,1) - pc(19,10,1)

sh(14,7) = 2*pc(20,10,1) - 3*pc(13,10,1) - 3*pc(15,10,1)

sh(15,7) = 4*pc(16,10,1) - pc(18,10,1) - pc(14,10,1)

sh(16,7) = 4*pc(17,10,1) - pc(12,10,1) - pc(19,10,1)

sh(1,8) = pc(1,5,1) - pc(1,6,1)

sh(2,8) = pc(2,5,1) - pc(2,6,1)

sh(3,8) = pc(3,5,1) - pc(3,6,1)

sh(4,8) = pc(4,5,1) - pc(4,6,1)

sh(5,8) = pc(8,5,1) - pc(8,6,1)

sh(6,8) = pc(9,5,1) - pc(9,6,1)

sh(7,8) = pc(10,5,1) - pc(10,6,1)

sh(8,8) = pc(5,5,1) - pc(5,6,1) - pc(6,5,1) + pc(6,6,1)

sh(9,8) = 2*pc(7,5,1) - 2*pc(7,6,1) - pc(5,5,1) + pc(5,6,1) - pc(6,5,1) + pc(6&
&,6,1)

sh(10,8) = pc(11,5,1) - pc(11,6,1)

sh(11,8) = pc(13,5,1) - pc(13,6,1) - pc(15,5,1) + pc(15,6,1)

sh(12,8) = pc(18,5,1) - pc(18,6,1) - 3*pc(14,5,1) + 3*pc(14,6,1)

sh(13,8) = 3*pc(12,5,1) - 3*pc(12,6,1) - pc(19,5,1) + pc(19,6,1)

sh(14,8) = 2*pc(20,5,1) - 2*pc(20,6,1) - 3*pc(13,5,1) + 3*pc(13,6,1) - 3*pc(15&
&,5,1) + 3*pc(15,6,1)

sh(15,8) = 4*pc(16,5,1) - 4*pc(16,6,1) - pc(18,5,1) + pc(18,6,1) - pc(14,5,1) &
&+ pc(14,6,1)

sh(16,8) = 4*pc(17,5,1) - 4*pc(17,6,1) - pc(12,5,1) + pc(12,6,1) - pc(19,5,1) &
&+ pc(19,6,1)

sh(1,9) = 2*pc(1,7,1) - pc(1,5,1) - pc(1,6,1)

sh(2,9) = 2*pc(2,7,1) - pc(2,5,1) - pc(2,6,1)

sh(3,9) = 2*pc(3,7,1) - pc(3,5,1) - pc(3,6,1)

sh(4,9) = 2*pc(4,7,1) - pc(4,5,1) - pc(4,6,1)

sh(5,9) = 2*pc(8,7,1) - pc(8,5,1) - pc(8,6,1)

sh(6,9) = 2*pc(9,7,1) - pc(9,5,1) - pc(9,6,1)

sh(7,9) = 2*pc(10,7,1) - pc(10,5,1) - pc(10,6,1)

sh(8,9) = 2*pc(5,7,1) - pc(5,5,1) - pc(5,6,1) - 2*pc(6,7,1) + pc(6,5,1) + pc(6&
&,6,1)

sh(9,9) = 4*pc(7,7,1) - 2*pc(7,5,1) - 2*pc(7,6,1) - 2*pc(5,7,1) + pc(5,5,1) + &
&pc(5,6,1) - 2*pc(6,7,1) + pc(6,5,1) + pc(6,6,1)

sh(10,9) = 2*pc(11,7,1) - pc(11,5,1) - pc(11,6,1)

sh(11,9) = 2*pc(13,7,1) - pc(13,5,1) - pc(13,6,1) - 2*pc(15,7,1) + pc(15,5,1) &
&+ pc(15,6,1)

sh(12,9) = 2*pc(18,7,1) - pc(18,5,1) - pc(18,6,1) - 6*pc(14,7,1) + 3*pc(14,5,1&
&) + 3*pc(14,6,1)

sh(13,9) = 6*pc(12,7,1) - 3*pc(12,5,1) - 3*pc(12,6,1) - 2*pc(19,7,1) + pc(19,5&
&,1) + pc(19,6,1)

sh(14,9) = 4*pc(20,7,1) - 2*pc(20,5,1) - 2*pc(20,6,1) - 6*pc(13,7,1) + 3*pc(13&
&,5,1) + 3*pc(13,6,1) - 6*pc(15,7,1) + 3*pc(15,5,1) + 3*pc(15,6,1)

sh(15,9) = 8*pc(16,7,1) - 4*pc(16,5,1) - 4*pc(16,6,1) - 2*pc(18,7,1) + pc(18,5&
&,1) + pc(18,6,1) - 2*pc(14,7,1) + pc(14,5,1) + pc(14,6,1)

sh(16,9) = 8*pc(17,7,1) - 4*pc(17,5,1) - 4*pc(17,6,1) - 2*pc(12,7,1) + pc(12,5&
&,1) + pc(12,6,1) - 2*pc(19,7,1) + pc(19,5,1) + pc(19,6,1)

sh(1,10) = pc(1,11,1)

sh(2,10) = pc(2,11,1)

sh(3,10) = pc(3,11,1)

sh(4,10) = pc(4,11,1)

sh(5,10) = pc(8,11,1)

sh(6,10) = pc(9,11,1)

sh(7,10) = pc(10,11,1)

sh(8,10) = pc(5,11,1) - pc(6,11,1)

sh(9,10) = 2*pc(7,11,1) - pc(5,11,1) - pc(6,11,1)

sh(10,10) = pc(11,11,1)

sh(11,10) = pc(13,11,1) - pc(15,11,1)

sh(12,10) = pc(18,11,1) - 3*pc(14,11,1)

sh(13,10) = 3*pc(12,11,1) - pc(19,11,1)

sh(14,10) = 2*pc(20,11,1) - 3*pc(13,11,1) - 3*pc(15,11,1)

sh(15,10) = 4*pc(16,11,1) - pc(18,11,1) - pc(14,11,1)

sh(16,10) = 4*pc(17,11,1) - pc(12,11,1) - pc(19,11,1)

sh(1,11) = pc(1,13,1) - pc(1,15,1)

sh(2,11) = pc(2,13,1) - pc(2,15,1)

sh(3,11) = pc(3,13,1) - pc(3,15,1)

sh(4,11) = pc(4,13,1) - pc(4,15,1)

sh(5,11) = pc(8,13,1) - pc(8,15,1)

sh(6,11) = pc(9,13,1) - pc(9,15,1)

sh(7,11) = pc(10,13,1) - pc(10,15,1)

sh(8,11) = pc(5,13,1) - pc(5,15,1) - pc(6,13,1) + pc(6,15,1)

sh(9,11) = 2*pc(7,13,1) - 2*pc(7,15,1) - pc(5,13,1) + pc(5,15,1) - pc(6,13,1) &
&+ pc(6,15,1)

sh(10,11) = pc(11,13,1) - pc(11,15,1)

sh(11,11) = pc(13,13,1) - pc(13,15,1) - pc(15,13,1) + pc(15,15,1)

sh(12,11) = pc(18,13,1) - pc(18,15,1) - 3*pc(14,13,1) + 3*pc(14,15,1)

sh(13,11) = 3*pc(12,13,1) - 3*pc(12,15,1) - pc(19,13,1) + pc(19,15,1)

sh(14,11) = 2*pc(20,13,1) - 2*pc(20,15,1) - 3*pc(13,13,1) + 3*pc(13,15,1) - 3*&
&pc(15,13,1) + 3*pc(15,15,1)

sh(15,11) = 4*pc(16,13,1) - 4*pc(16,15,1) - pc(18,13,1) + pc(18,15,1) - pc(14,&
&13,1) + pc(14,15,1)

sh(16,11) = 4*pc(17,13,1) - 4*pc(17,15,1) - pc(12,13,1) + pc(12,15,1) - pc(19,&
&13,1) + pc(19,15,1)

sh(1,12) = pc(1,18,1) - 3*pc(1,14,1)

sh(2,12) = pc(2,18,1) - 3*pc(2,14,1)

sh(3,12) = pc(3,18,1) - 3*pc(3,14,1)

sh(4,12) = pc(4,18,1) - 3*pc(4,14,1)

sh(5,12) = pc(8,18,1) - 3*pc(8,14,1)

sh(6,12) = pc(9,18,1) - 3*pc(9,14,1)

sh(7,12) = pc(10,18,1) - 3*pc(10,14,1)

sh(8,12) = pc(5,18,1) - 3*pc(5,14,1) - pc(6,18,1) + 3*pc(6,14,1)

sh(9,12) = 2*pc(7,18,1) - 6*pc(7,14,1) - pc(5,18,1) + 3*pc(5,14,1) - pc(6,18,1&
&) + 3*pc(6,14,1)

sh(10,12) = pc(11,18,1) - 3*pc(11,14,1)

sh(11,12) = pc(13,18,1) - 3*pc(13,14,1) - pc(15,18,1) + 3*pc(15,14,1)

sh(12,12) = pc(18,18,1) - 3*pc(18,14,1) - 3*pc(14,18,1) + 9*pc(14,14,1)

sh(13,12) = 3*pc(12,18,1) - 9*pc(12,14,1) - pc(19,18,1) + 3*pc(19,14,1)

sh(14,12) = 2*pc(20,18,1) - 6*pc(20,14,1) - 3*pc(13,18,1) + 9*pc(13,14,1) - 3*&
&pc(15,18,1) + 9*pc(15,14,1)

sh(15,12) = 4*pc(16,18,1) - 12*pc(16,14,1) - pc(18,18,1) + 3*pc(18,14,1) - pc(&
&14,18,1) + 3*pc(14,14,1)

sh(16,12) = 4*pc(17,18,1) - 12*pc(17,14,1) - pc(12,18,1) + 3*pc(12,14,1) - pc(&
&19,18,1) + 3*pc(19,14,1)

sh(1,13) = 3*pc(1,12,1) - pc(1,19,1)

sh(2,13) = 3*pc(2,12,1) - pc(2,19,1)

sh(3,13) = 3*pc(3,12,1) - pc(3,19,1)

sh(4,13) = 3*pc(4,12,1) - pc(4,19,1)

sh(5,13) = 3*pc(8,12,1) - pc(8,19,1)

sh(6,13) = 3*pc(9,12,1) - pc(9,19,1)

sh(7,13) = 3*pc(10,12,1) - pc(10,19,1)

sh(8,13) = 3*pc(5,12,1) - pc(5,19,1) - 3*pc(6,12,1) + pc(6,19,1)

sh(9,13) = 6*pc(7,12,1) - 2*pc(7,19,1) - 3*pc(5,12,1) + pc(5,19,1) - 3*pc(6,12&
&,1) + pc(6,19,1)

sh(10,13) = 3*pc(11,12,1) - pc(11,19,1)

sh(11,13) = 3*pc(13,12,1) - pc(13,19,1) - 3*pc(15,12,1) + pc(15,19,1)

sh(12,13) = 3*pc(18,12,1) - pc(18,19,1) - 9*pc(14,12,1) + 3*pc(14,19,1)

sh(13,13) = 9*pc(12,12,1) - 3*pc(12,19,1) - 3*pc(19,12,1) + pc(19,19,1)

sh(14,13) = 6*pc(20,12,1) - 2*pc(20,19,1) - 9*pc(13,12,1) + 3*pc(13,19,1) - 9*&
&pc(15,12,1) + 3*pc(15,19,1)

sh(15,13) = 12*pc(16,12,1) - 4*pc(16,19,1) - 3*pc(18,12,1) + pc(18,19,1) - 3*p&
&c(14,12,1) + pc(14,19,1)

sh(16,13) = 12*pc(17,12,1) - 4*pc(17,19,1) - 3*pc(12,12,1) + pc(12,19,1) - 3*p&
&c(19,12,1) + pc(19,19,1)

sh(1,14) = 2*pc(1,20,1) - 3*pc(1,13,1) - 3*pc(1,15,1)

sh(2,14) = 2*pc(2,20,1) - 3*pc(2,13,1) - 3*pc(2,15,1)

sh(3,14) = 2*pc(3,20,1) - 3*pc(3,13,1) - 3*pc(3,15,1)

sh(4,14) = 2*pc(4,20,1) - 3*pc(4,13,1) - 3*pc(4,15,1)

sh(5,14) = 2*pc(8,20,1) - 3*pc(8,13,1) - 3*pc(8,15,1)

sh(6,14) = 2*pc(9,20,1) - 3*pc(9,13,1) - 3*pc(9,15,1)

sh(7,14) = 2*pc(10,20,1) - 3*pc(10,13,1) - 3*pc(10,15,1)

sh(8,14) = 2*pc(5,20,1) - 3*pc(5,13,1) - 3*pc(5,15,1) - 2*pc(6,20,1) + 3*pc(6,&
&13,1) + 3*pc(6,15,1)

sh(9,14) = 4*pc(7,20,1) - 6*pc(7,13,1) - 6*pc(7,15,1) - 2*pc(5,20,1) + 3*pc(5,&
&13,1) + 3*pc(5,15,1) - 2*pc(6,20,1) + 3*pc(6,13,1) + 3*pc(6,15,1)

sh(10,14) = 2*pc(11,20,1) - 3*pc(11,13,1) - 3*pc(11,15,1)

sh(11,14) = 2*pc(13,20,1) - 3*pc(13,13,1) - 3*pc(13,15,1) - 2*pc(15,20,1) + 3*&
&pc(15,13,1) + 3*pc(15,15,1)

sh(12,14) = 2*pc(18,20,1) - 3*pc(18,13,1) - 3*pc(18,15,1) - 6*pc(14,20,1) + 9*&
&pc(14,13,1) + 9*pc(14,15,1)

sh(13,14) = 6*pc(12,20,1) - 9*pc(12,13,1) - 9*pc(12,15,1) - 2*pc(19,20,1) + 3*&
&pc(19,13,1) + 3*pc(19,15,1)

sh(14,14) = 4*pc(20,20,1) - 6*pc(20,13,1) - 6*pc(20,15,1) - 6*pc(13,20,1) + 9*&
&pc(13,13,1) + 9*pc(13,15,1) - 6*pc(15,20,1) + 9*pc(15,13,1) + 9*pc(15,15,1)

sh(15,14) = 8*pc(16,20,1) - 12*pc(16,13,1) - 12*pc(16,15,1) - 2*pc(18,20,1) + &
&3*pc(18,13,1) + 3*pc(18,15,1) - 2*pc(14,20,1) + 3*pc(14,13,1) + 3*pc(14,15,1)

sh(16,14) = 8*pc(17,20,1) - 12*pc(17,13,1) - 12*pc(17,15,1) - 2*pc(12,20,1) + &
&3*pc(12,13,1) + 3*pc(12,15,1) - 2*pc(19,20,1) + 3*pc(19,13,1) + 3*pc(19,15,1)

sh(1,15) = 4*pc(1,16,1) - pc(1,18,1) - pc(1,14,1)

sh(2,15) = 4*pc(2,16,1) - pc(2,18,1) - pc(2,14,1)

sh(3,15) = 4*pc(3,16,1) - pc(3,18,1) - pc(3,14,1)

sh(4,15) = 4*pc(4,16,1) - pc(4,18,1) - pc(4,14,1)

sh(5,15) = 4*pc(8,16,1) - pc(8,18,1) - pc(8,14,1)

sh(6,15) = 4*pc(9,16,1) - pc(9,18,1) - pc(9,14,1)

sh(7,15) = 4*pc(10,16,1) - pc(10,18,1) - pc(10,14,1)

sh(8,15) = 4*pc(5,16,1) - pc(5,18,1) - pc(5,14,1) - 4*pc(6,16,1) + pc(6,18,1) &
&+ pc(6,14,1)

sh(9,15) = 8*pc(7,16,1) - 2*pc(7,18,1) - 2*pc(7,14,1) - 4*pc(5,16,1) + pc(5,18&
&,1) + pc(5,14,1) - 4*pc(6,16,1) + pc(6,18,1) + pc(6,14,1)

sh(10,15) = 4*pc(11,16,1) - pc(11,18,1) - pc(11,14,1)

sh(11,15) = 4*pc(13,16,1) - pc(13,18,1) - pc(13,14,1) - 4*pc(15,16,1) + pc(15,&
&18,1) + pc(15,14,1)

sh(12,15) = 4*pc(18,16,1) - pc(18,18,1) - pc(18,14,1) - 12*pc(14,16,1) + 3*pc(&
&14,18,1) + 3*pc(14,14,1)

sh(13,15) = 12*pc(12,16,1) - 3*pc(12,18,1) - 3*pc(12,14,1) - 4*pc(19,16,1) + p&
&c(19,18,1) + pc(19,14,1)

sh(14,15) = 8*pc(20,16,1) - 2*pc(20,18,1) - 2*pc(20,14,1) - 12*pc(13,16,1) + 3&
&*pc(13,18,1) + 3*pc(13,14,1) - 12*pc(15,16,1) + 3*pc(15,18,1) + 3*pc(15,14,1)

sh(15,15) = 16*pc(16,16,1) - 4*pc(16,18,1) - 4*pc(16,14,1) - 4*pc(18,16,1) + p&
&c(18,18,1) + pc(18,14,1) - 4*pc(14,16,1) + pc(14,18,1) + pc(14,14,1)

sh(16,15) = 16*pc(17,16,1) - 4*pc(17,18,1) - 4*pc(17,14,1) - 4*pc(12,16,1) + p&
&c(12,18,1) + pc(12,14,1) - 4*pc(19,16,1) + pc(19,18,1) + pc(19,14,1)

sh(1,16) = 4*pc(1,17,1) - pc(1,12,1) - pc(1,19,1)

sh(2,16) = 4*pc(2,17,1) - pc(2,12,1) - pc(2,19,1)

sh(3,16) = 4*pc(3,17,1) - pc(3,12,1) - pc(3,19,1)

sh(4,16) = 4*pc(4,17,1) - pc(4,12,1) - pc(4,19,1)

sh(5,16) = 4*pc(8,17,1) - pc(8,12,1) - pc(8,19,1)

sh(6,16) = 4*pc(9,17,1) - pc(9,12,1) - pc(9,19,1)

sh(7,16) = 4*pc(10,17,1) - pc(10,12,1) - pc(10,19,1)

sh(8,16) = 4*pc(5,17,1) - pc(5,12,1) - pc(5,19,1) - 4*pc(6,17,1) + pc(6,12,1) &
&+ pc(6,19,1)

sh(9,16) = 8*pc(7,17,1) - 2*pc(7,12,1) - 2*pc(7,19,1) - 4*pc(5,17,1) + pc(5,12&
&,1) + pc(5,19,1) - 4*pc(6,17,1) + pc(6,12,1) + pc(6,19,1)

sh(10,16) = 4*pc(11,17,1) - pc(11,12,1) - pc(11,19,1)

sh(11,16) = 4*pc(13,17,1) - pc(13,12,1) - pc(13,19,1) - 4*pc(15,17,1) + pc(15,&
&12,1) + pc(15,19,1)

sh(12,16) = 4*pc(18,17,1) - pc(18,12,1) - pc(18,19,1) - 12*pc(14,17,1) + 3*pc(&
&14,12,1) + 3*pc(14,19,1)

sh(13,16) = 12*pc(12,17,1) - 3*pc(12,12,1) - 3*pc(12,19,1) - 4*pc(19,17,1) + p&
&c(19,12,1) + pc(19,19,1)

sh(14,16) = 8*pc(20,17,1) - 2*pc(20,12,1) - 2*pc(20,19,1) - 12*pc(13,17,1) + 3&
&*pc(13,12,1) + 3*pc(13,19,1) - 12*pc(15,17,1) + 3*pc(15,12,1) + 3*pc(15,19,1)

sh(15,16) = 16*pc(16,17,1) - 4*pc(16,12,1) - 4*pc(16,19,1) - 4*pc(18,17,1) + p&
&c(18,12,1) + pc(18,19,1) - 4*pc(14,17,1) + pc(14,12,1) + pc(14,19,1)

sh(16,16) = 16*pc(17,17,1) - 4*pc(17,12,1) - 4*pc(17,19,1) - 4*pc(12,17,1) + p&
&c(12,12,1) + pc(12,19,1) - 4*pc(19,17,1) + pc(19,12,1) + pc(19,19,1)


   end subroutine nuclear3CIntgAna

   ! This subroutine computes F_N(T) - which is the Incomplete
   ! Gamma Function, aka the Boys Function. This comes about as a 
   ! result of the integral transforms that are required to 
   ! account for the (1/r) term that is present in the Coulomb 
   ! integral, to make the integral separable in the Cartesian 
   ! directions (x,y,z).  The formalism is Eq. (29) on page 41 of the 
   ! following document:
   !
   ! ****************************************************************
   ! Petersson T., Hellsing B., "A Detailed Derivation of Gaussian
   ! Orbital-Based Matrix Elements In Electron Structure Calculations."
   ! European Journal of Physics, Vol. 31, No. 1, Page 37.  
   ! ****************************************************************
   !
   ! The Boys Function itself is defined as:
   ! F_N(T) = integral [from 0 to 1] {(t**(2*N))*exp(-T*(t**2))}dt
   !
   ! where T = (a1 + a2 + a3)*sum(P - C)**2
   !
   ! For our purposes here, N will be an integer >=0 that is a linear
   ! combination of index values that come from the summations that
   ! compute the main body of the Coulomb integral, be it Nuclear
   ! Attraction, or the more complicated Electron Repulsion Integral.
   !
   ! Generally, the solution to this integral is expressed numerically
   ! in terms of standard error functions (erf), but this closed form
   ! expression becomes numerically unstable if T = 0 or is very small 
   ! (~10E-2 or less).  To avoid the large roundoff errors that occur
   ! as a result of this, a check has to be performed on the value of T.  
   ! If T <= 10E-2, the following formalism shown by Cook (referenced below)
   ! needs to be applied (Cook: p. 260):
   !
   ! F_N(T) = (1/2)*exp(-T)*sum[from i=0 to infinity] of
   !                            {gamma(N+1/2)/gamma(N+i+3/2)}
   !
   ! This infinite series is truncated at i = 6 for our desired accuracy,
   ! including for the case where T = 0.
   !
   ! More detailed in formation on the Boys Function and numerical 
   ! approximations of it are detailed in the following sources:
   !
   ! ****************************************************************
   ! Cook, David B., "Handbook of Computational Quantum Chemistry."
   ! Dover Publications (August 22, 2005).
   ! ISBN-10: 0486443078
   ! ISBN-13: 978-0486443072
   ! ****************************************************************
   !
   ! as well as:
   !
   ! ****************************************************************
   ! Helgaker T., Taylor P., "Gaussian Basis Sets and Molecular 
   ! Integrals", Ch. 12 from "Modern Electronic Structure Theory."
   ! Pp. 809 - 815.
   ! ****************************************************************
   subroutine boys(T,F)
   use o_kinds
   use o_constants, only: pi
   implicit none
   
   ! Define passed parameters
   real (kind=double), intent (in) :: T
   real (kind=double), dimension(7), intent (out) :: F

   ! Define local variables
   integer :: N
   real (kind=double) :: sqrt_pi, erf_T, exp_T
   real (kind=double) :: a0, a1, a2, a3, a4, a5, a6, a7
   real (kind=double), dimension(8) :: S

   ! Initialize local variables that are used regardless of the method used
   !   to evaluate the boys function (expansion or analytic).
   exp_T = exp(-T)

   sqrt_pi = sqrt(pi) 
   erf_T = erf(sqrt(T))
   ! As stated earler, if T is small (<= 10E-2), the closed form solution 
   !   to the Boys Function can't be used because it becomes numerically
   !   unstable.  Instead, the following expansion is used to approximate it
   !   for small values (Cook: p. 260).  The loop defines the Nth degree Boys
   !   Function.

   if (T <= 10E-2) then
      ! The values (a0 through a7) are computed gamma(N + 0.5) values
      ! to be plugged into the below approximation of the boys Function:

      ! N = 0

      do N = 0, 6
         !
         ! S(N + 1) = 0.5d0*exp_T*(gamma(N + 0.5d0)/gamma(N + 1.5d0)&
         !   &+ (gamma(N + 0.5d0)*T)/gamma(N + 2.5d0)& 
         !   &+ (gamma(N + 0.5d0)*T**2)/gamma(N + 3.5d0)& 
         !   &+ (gamma(N + 0.5d0)*T**3)/gamma(N + 4.5d0)&
         !   &+ (gamma(N + 0.5d0)*T**4)/gamma(N + 5.5d0)& 
         !   &+ (gamma(N + 0.5d0)*T**5)/gamma(N + 6.5d0)& 
         !   &+ (gamma(N + 0.5d0)*T**6)/gamma(N + 7.5d0))

         ! This approximation is explained in the above documentation,
         ! and is required when T <= 10E-2..
    if (N == 0) then
      a0 = 1.7724538509055161     
      a1 = 0.88622692545275805     
      a2 = 1.3293403881791370     
      a3 = 3.3233509704478426     
      a4 = 11.631728396567450     
      a5 = 52.342777784553533     
      a6 = 287.88527781504416     
      a7 = 1871.2543057977896     
    else if (N == 1) then
      a0 = 0.88622692545275805     
      a1 = 1.3293403881791370     
      a2 = 3.3233509704478426     
      a3 = 11.631728396567450     
      a4 = 52.342777784553533     
      a5 = 287.88527781504416     
      a6 = 1871.2543057977896     
      a7 = 14034.407293483402     
    else if (N == 2) then
      a0 = 1.3293403881791370     
      a1 = 3.3233509704478426     
      a2 = 11.631728396567450     
      a3 = 52.342777784553533     
      a4 = 287.88527781504416     
      a5 = 1871.2543057977896     
      a6 = 14034.407293483402     
      a7 = 119292.46199460901     
    else if (N == 3) then
      a0 = 3.3233509704478426     
      a1 = 11.631728396567450     
      a2 = 52.342777784553533     
      a3 = 287.88527781504416     
      a4 = 1871.2543057977896     
      a5 = 14034.407293483402     
      a6 = 119292.46199460901     
      a7 = 1133278.3889487833     
    else if (N == 4) then
      a0 = 11.631728396567450     
      a1 = 52.342777784553533     
      a2 = 287.88527781504416     
      a3 = 1871.2543057977896     
      a4 = 14034.407293483402     
      a5 = 119292.46199460901     
      a6 = 1133278.3889487833     
      a7 = 11899423.083962219     
    else if (N == 5) then
      a0 = 52.342777784553533     
      a1 = 287.88527781504416     
      a2 = 1871.2543057977896     
      a3 = 14034.407293483402     
      a4 = 119292.46199460901     
      a5 = 1133278.3889487833     
      a6 = 11899423.083962219     
      a7 = 136843365.46556622     
    else if (N == 6) then  
      a0 = 287.88527781504416     
      a1 = 1871.2543057977896     
      a2 = 14034.407293483402     
      a3 = 119292.46199460901     
      a4 = 1133278.3889487833     
      a5 = 11899423.083962219     
      a6 = 136843365.46556622     
      a7 = 1710542068.3195722
    end if
    S(N + 1) = 0.5d0*exp_T*a0*((1/a1) + (T/a2) + ((T**2)/a3)&
      &+ ((T**3)/a4) + ((T**4)/a5) + ((T**5)/a6) + ((T**6)/a7))
  end do
  F(1) = S(1)
  F(2) = S(2)
  F(3) = S(3)
  F(4) = S(4)
  F(5) = S(5)
  F(6) = S(6)
  F(7) = S(7)
else
  F(1) = 0.5*((sqrt_pi/(T**(0.5d0)))*erf_T)

  F(2) = ((sqrt_pi/(4*T**(1.5d0)))*erf_T-exp_T*(1/(2d0*T)))
  
  F(3) = 6*((sqrt_pi/(16*T**(2.5d0)))*erf_T-exp_T*(2/(24d0*T)+1/((4)*2d0*T**(2&
  &))))
  
  F(4) = 60*((sqrt_pi/(64*T**(3.5d0)))*erf_T-exp_T*(6/(720d0*T)+2/((4)*24d0*T*&
  &*(2))+1/((16)*2d0*T**(3))))
  
  F(5) = 840*((sqrt_pi/(256*T**(4.5d0)))*erf_T-exp_T*(24/(40320d0*T)+6/((4)*720&
  &d0*T**(2))+2/((16)*24d0*T**(3))+1/((64)*2d0*T**(4))))
  
  F(6) = 15120*((sqrt_pi/(1024*T**(5.5d0)))*erf_T-exp_T*(120/(3628800d0*T)+24/(&
  &(4)*40320d0*T**(2))+6/((16)*720d0*T**(3))+2/((64)*24d0*T**(4))+1/((256)*2d0&
  &*T**(5))))
  
  F(7) = 332640*((sqrt_pi/(4096*T**(6.5d0)))*erf_T-exp_T*(720/(479001600d0*T)+1&
  &20/((4)*3628800d0*T**(2))+24/((16)*40320d0*T**(3))+6/((64)*720d0*T**(4))+2/&
  &((256)*24d0*T**(5))+1/((1024)*2d0*T**(6))))
end if
end subroutine boys

   subroutine nuclear3CIntgNum(a1,a2,a3,A,B,C,pc,sh,cell_size,step_size)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2, a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   real (kind=double), dimension (20,20), intent(out) :: pc
   real (kind=double), dimension (16,16), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i, j, k
   integer :: num_steps
   integer, dimension (20,3) :: triads
   integer, dimension (16,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, r, soln
   real (kind=double) :: xC_dist_sqrd, yC_dist_sqrd, zC_dist_sqrd
   real (kind=double), dimension (3) :: xyz, xyz_sum

   ! Initialize local variables.
   triads(1,:) = (/0,0,0/)
   triads(2,:) = (/1,0,0/)
   triads(3,:) = (/0,1,0/)
   triads(4,:) = (/0,0,1/)
   triads(5,:) = (/2,0,0/)
   triads(6,:) = (/0,2,0/)
   triads(7,:) = (/0,0,2/)
   triads(8,:) = (/1,1,0/)
   triads(9,:) = (/1,0,1/)
   triads(10,:) = (/0,1,1/)
   triads(11,:) = (/1,1,1/)
   triads(12,:) = (/2,1,0/)
   triads(13,:) = (/2,0,1/)
   triads(14,:) = (/1,2,0/)
   triads(15,:) = (/0,2,1/)
   triads(16,:) = (/1,0,2/)
   triads(17,:) = (/0,1,2/)
   triads(18,:) = (/3,0,0/)
   triads(19,:) = (/0,3,0/)
   triads(20,:) = (/0,0,3/)
   conversion(1,1,:) = (/1,0,0/)
   conversion(1,2,:) = (/1,1,1/)
   conversion(2,1,:) = (/1,0,0/)
   conversion(2,2,:) = (/2,1,1/)
   conversion(3,1,:) = (/1,0,0/)
   conversion(3,2,:) = (/3,1,1/)
   conversion(4,1,:) = (/1,0,0/)
   conversion(4,2,:) = (/4,1,1/)
   conversion(5,1,:) = (/1,0,0/)
   conversion(5,2,:) = (/8,1,1/)
   conversion(6,1,:) = (/1,0,0/)
   conversion(6,2,:) = (/9,1,1/)
   conversion(7,1,:) = (/1,0,0/)
   conversion(7,2,:) = (/10,1,1/)
   conversion(8,1,:) = (/1,-1,0/)
   conversion(8,2,:) = (/5,6,1/)
   conversion(9,1,:) = (/2,-1,-1/)
   conversion(9,2,:) = (/7,5,6/)
   conversion(10,1,:) = (/1,0,0/)
   conversion(10,2,:) = (/11,1,1/)
   conversion(11,1,:) = (/1,-1,0/)
   conversion(11,2,:) = (/13,15,1/)
   conversion(12,1,:) = (/1,-3,0/)
   conversion(12,2,:) = (/18,14,1/)
   conversion(13,1,:) = (/3,-1,0/)
   conversion(13,2,:) = (/12,19,1/)
   conversion(14,1,:) = (/2,-3,-3/)
   conversion(14,2,:) = (/20,13,15/)
   conversion(15,1,:) = (/4,-1,-1/)
   conversion(15,2,:) = (/16,18,14/)
   conversion(16,1,:) = (/4,-1,-1/)
   conversion(16,2,:) = (/17,12,19/)

   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   ! Original equation from sympy:
   !   (xyz(1)-A(1))**l1(1)*(xyz(1)-B(1))**l2(1)*(xyz(2)-A(2))**l1(2)*
   !   (xyz(2)-B(2))**l2(2)*(xyz(3)-A(3))**l1(3)*(xyz(3)-B(3))**l2(3)*
   !   exp(-a1*sum((xyz(:)-A(:))**2)) * &
   !   exp(-a2*sum((xyz(:)-B(:))**2)) * &
   !   exp(-a3*sum((xyz(:)-C(:))**2)) / &
   !   sqrt(sum((xyz(:)-C(:))**2))

   ! By-hand algebraic rearrangement to the expression below for faster
   !   computation in the triple nested loop environment.

   do p = 1, 20
      do q = 1, 20

         ! Assign l1 and l2 values for each gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize the solution
         soln = 0.0d0

         do i = 0, num_steps

            ! Compute the current x position, distance and sum.
            xyz(1) = (start_pos + (i*step_size))
            xC_dist_sqrd = (xyz(1)-C(1))**2
            xyz_sum(1) = ((xyz(1) - A(1))**l1(1))*((xyz(1) - B(1))**l2(1)) &
                  & * exp(-a1*(xyz(1)-A(1))**2 + -a2*(xyz(1)-B(1))**2 + &
                  &       -a3*xC_dist_sqrd)

            do j = 0, num_steps

               ! Compute the current y position, distance, and sum.
               xyz(2) = (start_pos + (j*step_size))
               yC_dist_sqrd = (xyz(2)-C(2))**2
               xyz_sum(2) = ((xyz(2) - A(2))**l1(2)) &
                     & * ((xyz(2) - B(2))**l2(2)) &
                     & * exp(-a1*(xyz(2)-A(2))**2 + -a2*(xyz(2)-B(2))**2 + &
                     &       -a3*yC_dist_sqrd)

               do k = 0, num_steps

                  ! Compute the current z position.
                  xyz(3) = (start_pos + (k*step_size))
                  zC_dist_sqrd = (xyz(3) - C(3))**2

                  ! Compute the distance between the electron and the nucleus.
                  r = xC_dist_sqrd + yC_dist_sqrd + zC_dist_sqrd

                  ! If the distance is zero, then cycle.
                  if (r == 0) cycle

                  ! Compute the current z xyz_sum.
                  xyz_sum(3) = ((xyz(3) - A(3))**l1(3)) &
                        & * ((xyz(3) - B(3))**l2(3)) &
                        & * exp(-a1*(xyz(3)-A(3))**2 + -a2*(xyz(3)-B(3))**2 + &
                        &       -a3*zC_dist_sqrd)

                  soln = soln + product(xyz_sum(:)) / sqrt(r)
               enddo
            enddo
         enddo

         ! Multiply the results by the mesh step size.
         soln = soln * step_size**3

         pc(q,p) = soln
      enddo
    enddo

sh(1,1) = pc(1,1)

sh(2,1) = pc(2,1)

sh(3,1) = pc(3,1)

sh(4,1) = pc(4,1)

sh(5,1) = pc(8,1)

sh(6,1) = pc(9,1)

sh(7,1) = pc(10,1)

sh(8,1) = pc(5,1) - pc(6,1)

sh(9,1) = 2*pc(7,1) - pc(5,1) - pc(6,1)

sh(10,1) = pc(11,1)

sh(11,1) = pc(13,1) - pc(15,1)

sh(12,1) = pc(18,1) - 3*pc(14,1)

sh(13,1) = 3*pc(12,1) - pc(19,1)

sh(14,1) = 2*pc(20,1) - 3*pc(13,1) - 3*pc(15,1)

sh(15,1) = 4*pc(16,1) - pc(18,1) - pc(14,1)

sh(16,1) = 4*pc(17,1) - pc(12,1) - pc(19,1)

sh(1,2) = pc(1,2)

sh(2,2) = pc(2,2)

sh(3,2) = pc(3,2)

sh(4,2) = pc(4,2)

sh(5,2) = pc(8,2)

sh(6,2) = pc(9,2)

sh(7,2) = pc(10,2)

sh(8,2) = pc(5,2) - pc(6,2)

sh(9,2) = 2*pc(7,2) - pc(5,2) - pc(6,2)

sh(10,2) = pc(11,2)

sh(11,2) = pc(13,2) - pc(15,2)

sh(12,2) = pc(18,2) - 3*pc(14,2)

sh(13,2) = 3*pc(12,2) - pc(19,2)

sh(14,2) = 2*pc(20,2) - 3*pc(13,2) - 3*pc(15,2)

sh(15,2) = 4*pc(16,2) - pc(18,2) - pc(14,2)

sh(16,2) = 4*pc(17,2) - pc(12,2) - pc(19,2)

sh(1,3) = pc(1,3)

sh(2,3) = pc(2,3)

sh(3,3) = pc(3,3)

sh(4,3) = pc(4,3)

sh(5,3) = pc(8,3)

sh(6,3) = pc(9,3)

sh(7,3) = pc(10,3)

sh(8,3) = pc(5,3) - pc(6,3)

sh(9,3) = 2*pc(7,3) - pc(5,3) - pc(6,3)

sh(10,3) = pc(11,3)

sh(11,3) = pc(13,3) - pc(15,3)

sh(12,3) = pc(18,3) - 3*pc(14,3)

sh(13,3) = 3*pc(12,3) - pc(19,3)

sh(14,3) = 2*pc(20,3) - 3*pc(13,3) - 3*pc(15,3)

sh(15,3) = 4*pc(16,3) - pc(18,3) - pc(14,3)

sh(16,3) = 4*pc(17,3) - pc(12,3) - pc(19,3)

sh(1,4) = pc(1,4)

sh(2,4) = pc(2,4)

sh(3,4) = pc(3,4)

sh(4,4) = pc(4,4)

sh(5,4) = pc(8,4)

sh(6,4) = pc(9,4)

sh(7,4) = pc(10,4)

sh(8,4) = pc(5,4) - pc(6,4)

sh(9,4) = 2*pc(7,4) - pc(5,4) - pc(6,4)

sh(10,4) = pc(11,4)

sh(11,4) = pc(13,4) - pc(15,4)

sh(12,4) = pc(18,4) - 3*pc(14,4)

sh(13,4) = 3*pc(12,4) - pc(19,4)

sh(14,4) = 2*pc(20,4) - 3*pc(13,4) - 3*pc(15,4)

sh(15,4) = 4*pc(16,4) - pc(18,4) - pc(14,4)

sh(16,4) = 4*pc(17,4) - pc(12,4) - pc(19,4)

sh(1,5) = pc(1,8)

sh(2,5) = pc(2,8)

sh(3,5) = pc(3,8)

sh(4,5) = pc(4,8)

sh(5,5) = pc(8,8)

sh(6,5) = pc(9,8)

sh(7,5) = pc(10,8)

sh(8,5) = pc(5,8) - pc(6,8)

sh(9,5) = 2*pc(7,8) - pc(5,8) - pc(6,8)

sh(10,5) = pc(11,8)

sh(11,5) = pc(13,8) - pc(15,8)

sh(12,5) = pc(18,8) - 3*pc(14,8)

sh(13,5) = 3*pc(12,8) - pc(19,8)

sh(14,5) = 2*pc(20,8) - 3*pc(13,8) - 3*pc(15,8)

sh(15,5) = 4*pc(16,8) - pc(18,8) - pc(14,8)

sh(16,5) = 4*pc(17,8) - pc(12,8) - pc(19,8)

sh(1,6) = pc(1,9)

sh(2,6) = pc(2,9)

sh(3,6) = pc(3,9)

sh(4,6) = pc(4,9)

sh(5,6) = pc(8,9)

sh(6,6) = pc(9,9)

sh(7,6) = pc(10,9)

sh(8,6) = pc(5,9) - pc(6,9)

sh(9,6) = 2*pc(7,9) - pc(5,9) - pc(6,9)

sh(10,6) = pc(11,9)

sh(11,6) = pc(13,9) - pc(15,9)

sh(12,6) = pc(18,9) - 3*pc(14,9)

sh(13,6) = 3*pc(12,9) - pc(19,9)

sh(14,6) = 2*pc(20,9) - 3*pc(13,9) - 3*pc(15,9)

sh(15,6) = 4*pc(16,9) - pc(18,9) - pc(14,9)

sh(16,6) = 4*pc(17,9) - pc(12,9) - pc(19,9)

sh(1,7) = pc(1,10)

sh(2,7) = pc(2,10)

sh(3,7) = pc(3,10)

sh(4,7) = pc(4,10)

sh(5,7) = pc(8,10)

sh(6,7) = pc(9,10)

sh(7,7) = pc(10,10)

sh(8,7) = pc(5,10) - pc(6,10)

sh(9,7) = 2*pc(7,10) - pc(5,10) - pc(6,10)

sh(10,7) = pc(11,10)

sh(11,7) = pc(13,10) - pc(15,10)

sh(12,7) = pc(18,10) - 3*pc(14,10)

sh(13,7) = 3*pc(12,10) - pc(19,10)

sh(14,7) = 2*pc(20,10) - 3*pc(13,10) - 3*pc(15,10)

sh(15,7) = 4*pc(16,10) - pc(18,10) - pc(14,10)

sh(16,7) = 4*pc(17,10) - pc(12,10) - pc(19,10)

sh(1,8) = pc(1,5) - pc(1,6)

sh(2,8) = pc(2,5) - pc(2,6)

sh(3,8) = pc(3,5) - pc(3,6)

sh(4,8) = pc(4,5) - pc(4,6)

sh(5,8) = pc(8,5) - pc(8,6)

sh(6,8) = pc(9,5) - pc(9,6)

sh(7,8) = pc(10,5) - pc(10,6)

sh(8,8) = pc(5,5) - pc(5,6) - pc(6,5) + pc(6,6)

sh(9,8) = 2*pc(7,5) - 2*pc(7,6) - pc(5,5) + pc(5,6) - pc(6,5) + pc(6,6)

sh(10,8) = pc(11,5) - pc(11,6)

sh(11,8) = pc(13,5) - pc(13,6) - pc(15,5) + pc(15,6)

sh(12,8) = pc(18,5) - pc(18,6) - 3*pc(14,5) + 3*pc(14,6)

sh(13,8) = 3*pc(12,5) - 3*pc(12,6) - pc(19,5) + pc(19,6)

sh(14,8) = 2*pc(20,5) - 2*pc(20,6) - 3*pc(13,5) + 3*pc(13,6) - 3*pc(15,5) + 3*&
&pc(15,6)

sh(15,8) = 4*pc(16,5) - 4*pc(16,6) - pc(18,5) + pc(18,6) - pc(14,5) + pc(14,6)

sh(16,8) = 4*pc(17,5) - 4*pc(17,6) - pc(12,5) + pc(12,6) - pc(19,5) + pc(19,6)

sh(1,9) = 2*pc(1,7) - pc(1,5) - pc(1,6)

sh(2,9) = 2*pc(2,7) - pc(2,5) - pc(2,6)

sh(3,9) = 2*pc(3,7) - pc(3,5) - pc(3,6)

sh(4,9) = 2*pc(4,7) - pc(4,5) - pc(4,6)

sh(5,9) = 2*pc(8,7) - pc(8,5) - pc(8,6)

sh(6,9) = 2*pc(9,7) - pc(9,5) - pc(9,6)

sh(7,9) = 2*pc(10,7) - pc(10,5) - pc(10,6)

sh(8,9) = 2*pc(5,7) - pc(5,5) - pc(5,6) - 2*pc(6,7) + pc(6,5) + pc(6,6)

sh(9,9) = 4*pc(7,7) - 2*pc(7,5) - 2*pc(7,6) - 2*pc(5,7) + pc(5,5) + pc(5,6) - &
&2*pc(6,7) + pc(6,5) + pc(6,6)

sh(10,9) = 2*pc(11,7) - pc(11,5) - pc(11,6)

sh(11,9) = 2*pc(13,7) - pc(13,5) - pc(13,6) - 2*pc(15,7) + pc(15,5) + pc(15,6)

sh(12,9) = 2*pc(18,7) - pc(18,5) - pc(18,6) - 6*pc(14,7) + 3*pc(14,5) + 3*pc(1&
&4,6)

sh(13,9) = 6*pc(12,7) - 3*pc(12,5) - 3*pc(12,6) - 2*pc(19,7) + pc(19,5) + pc(1&
&9,6)

sh(14,9) = 4*pc(20,7) - 2*pc(20,5) - 2*pc(20,6) - 6*pc(13,7) + 3*pc(13,5) + 3*&
&pc(13,6) - 6*pc(15,7) + 3*pc(15,5) + 3*pc(15,6)

sh(15,9) = 8*pc(16,7) - 4*pc(16,5) - 4*pc(16,6) - 2*pc(18,7) + pc(18,5) + pc(1&
&8,6) - 2*pc(14,7) + pc(14,5) + pc(14,6)

sh(16,9) = 8*pc(17,7) - 4*pc(17,5) - 4*pc(17,6) - 2*pc(12,7) + pc(12,5) + pc(1&
&2,6) - 2*pc(19,7) + pc(19,5) + pc(19,6)

sh(1,10) = pc(1,11)

sh(2,10) = pc(2,11)

sh(3,10) = pc(3,11)

sh(4,10) = pc(4,11)

sh(5,10) = pc(8,11)

sh(6,10) = pc(9,11)

sh(7,10) = pc(10,11)

sh(8,10) = pc(5,11) - pc(6,11)

sh(9,10) = 2*pc(7,11) - pc(5,11) - pc(6,11)

sh(10,10) = pc(11,11)

sh(11,10) = pc(13,11) - pc(15,11)

sh(12,10) = pc(18,11) - 3*pc(14,11)

sh(13,10) = 3*pc(12,11) - pc(19,11)

sh(14,10) = 2*pc(20,11) - 3*pc(13,11) - 3*pc(15,11)

sh(15,10) = 4*pc(16,11) - pc(18,11) - pc(14,11)

sh(16,10) = 4*pc(17,11) - pc(12,11) - pc(19,11)

sh(1,11) = pc(1,13) - pc(1,15)

sh(2,11) = pc(2,13) - pc(2,15)

sh(3,11) = pc(3,13) - pc(3,15)

sh(4,11) = pc(4,13) - pc(4,15)

sh(5,11) = pc(8,13) - pc(8,15)

sh(6,11) = pc(9,13) - pc(9,15)

sh(7,11) = pc(10,13) - pc(10,15)

sh(8,11) = pc(5,13) - pc(5,15) - pc(6,13) + pc(6,15)

sh(9,11) = 2*pc(7,13) - 2*pc(7,15) - pc(5,13) + pc(5,15) - pc(6,13) + pc(6,15)

sh(10,11) = pc(11,13) - pc(11,15)

sh(11,11) = pc(13,13) - pc(13,15) - pc(15,13) + pc(15,15)

sh(12,11) = pc(18,13) - pc(18,15) - 3*pc(14,13) + 3*pc(14,15)

sh(13,11) = 3*pc(12,13) - 3*pc(12,15) - pc(19,13) + pc(19,15)

sh(14,11) = 2*pc(20,13) - 2*pc(20,15) - 3*pc(13,13) + 3*pc(13,15) - 3*pc(15,13&
&) + 3*pc(15,15)

sh(15,11) = 4*pc(16,13) - 4*pc(16,15) - pc(18,13) + pc(18,15) - pc(14,13) + pc&
&(14,15)

sh(16,11) = 4*pc(17,13) - 4*pc(17,15) - pc(12,13) + pc(12,15) - pc(19,13) + pc&
&(19,15)

sh(1,12) = pc(1,18) - 3*pc(1,14)

sh(2,12) = pc(2,18) - 3*pc(2,14)

sh(3,12) = pc(3,18) - 3*pc(3,14)

sh(4,12) = pc(4,18) - 3*pc(4,14)

sh(5,12) = pc(8,18) - 3*pc(8,14)

sh(6,12) = pc(9,18) - 3*pc(9,14)

sh(7,12) = pc(10,18) - 3*pc(10,14)

sh(8,12) = pc(5,18) - 3*pc(5,14) - pc(6,18) + 3*pc(6,14)

sh(9,12) = 2*pc(7,18) - 6*pc(7,14) - pc(5,18) + 3*pc(5,14) - pc(6,18) + 3*pc(6&
&,14)

sh(10,12) = pc(11,18) - 3*pc(11,14)

sh(11,12) = pc(13,18) - 3*pc(13,14) - pc(15,18) + 3*pc(15,14)

sh(12,12) = pc(18,18) - 3*pc(18,14) - 3*pc(14,18) + 9*pc(14,14)

sh(13,12) = 3*pc(12,18) - 9*pc(12,14) - pc(19,18) + 3*pc(19,14)

sh(14,12) = 2*pc(20,18) - 6*pc(20,14) - 3*pc(13,18) + 9*pc(13,14) - 3*pc(15,18&
&) + 9*pc(15,14)

sh(15,12) = 4*pc(16,18) - 12*pc(16,14) - pc(18,18) + 3*pc(18,14) - pc(14,18) +&
& 3*pc(14,14)

sh(16,12) = 4*pc(17,18) - 12*pc(17,14) - pc(12,18) + 3*pc(12,14) - pc(19,18) +&
& 3*pc(19,14)

sh(1,13) = 3*pc(1,12) - pc(1,19)

sh(2,13) = 3*pc(2,12) - pc(2,19)

sh(3,13) = 3*pc(3,12) - pc(3,19)

sh(4,13) = 3*pc(4,12) - pc(4,19)

sh(5,13) = 3*pc(8,12) - pc(8,19)

sh(6,13) = 3*pc(9,12) - pc(9,19)

sh(7,13) = 3*pc(10,12) - pc(10,19)

sh(8,13) = 3*pc(5,12) - pc(5,19) - 3*pc(6,12) + pc(6,19)

sh(9,13) = 6*pc(7,12) - 2*pc(7,19) - 3*pc(5,12) + pc(5,19) - 3*pc(6,12) + pc(6&
&,19)

sh(10,13) = 3*pc(11,12) - pc(11,19)

sh(11,13) = 3*pc(13,12) - pc(13,19) - 3*pc(15,12) + pc(15,19)

sh(12,13) = 3*pc(18,12) - pc(18,19) - 9*pc(14,12) + 3*pc(14,19)

sh(13,13) = 9*pc(12,12) - 3*pc(12,19) - 3*pc(19,12) + pc(19,19)

sh(14,13) = 6*pc(20,12) - 2*pc(20,19) - 9*pc(13,12) + 3*pc(13,19) - 9*pc(15,12&
&) + 3*pc(15,19)

sh(15,13) = 12*pc(16,12) - 4*pc(16,19) - 3*pc(18,12) + pc(18,19) - 3*pc(14,12)&
& + pc(14,19)

sh(16,13) = 12*pc(17,12) - 4*pc(17,19) - 3*pc(12,12) + pc(12,19) - 3*pc(19,12)&
& + pc(19,19)

sh(1,14) = 2*pc(1,20) - 3*pc(1,13) - 3*pc(1,15)

sh(2,14) = 2*pc(2,20) - 3*pc(2,13) - 3*pc(2,15)

sh(3,14) = 2*pc(3,20) - 3*pc(3,13) - 3*pc(3,15)

sh(4,14) = 2*pc(4,20) - 3*pc(4,13) - 3*pc(4,15)

sh(5,14) = 2*pc(8,20) - 3*pc(8,13) - 3*pc(8,15)

sh(6,14) = 2*pc(9,20) - 3*pc(9,13) - 3*pc(9,15)

sh(7,14) = 2*pc(10,20) - 3*pc(10,13) - 3*pc(10,15)

sh(8,14) = 2*pc(5,20) - 3*pc(5,13) - 3*pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*&
&pc(6,15)

sh(9,14) = 4*pc(7,20) - 6*pc(7,13) - 6*pc(7,15) - 2*pc(5,20) + 3*pc(5,13) + 3*&
&pc(5,15) - 2*pc(6,20) + 3*pc(6,13) + 3*pc(6,15)

sh(10,14) = 2*pc(11,20) - 3*pc(11,13) - 3*pc(11,15)

sh(11,14) = 2*pc(13,20) - 3*pc(13,13) - 3*pc(13,15) - 2*pc(15,20) + 3*pc(15,13&
&) + 3*pc(15,15)

sh(12,14) = 2*pc(18,20) - 3*pc(18,13) - 3*pc(18,15) - 6*pc(14,20) + 9*pc(14,13&
&) + 9*pc(14,15)

sh(13,14) = 6*pc(12,20) - 9*pc(12,13) - 9*pc(12,15) - 2*pc(19,20) + 3*pc(19,13&
&) + 3*pc(19,15)

sh(14,14) = 4*pc(20,20) - 6*pc(20,13) - 6*pc(20,15) - 6*pc(13,20) + 9*pc(13,13&
&) + 9*pc(13,15) - 6*pc(15,20) + 9*pc(15,13) + 9*pc(15,15)

sh(15,14) = 8*pc(16,20) - 12*pc(16,13) - 12*pc(16,15) - 2*pc(18,20) + 3*pc(18,&
&13) + 3*pc(18,15) - 2*pc(14,20) + 3*pc(14,13) + 3*pc(14,15)

sh(16,14) = 8*pc(17,20) - 12*pc(17,13) - 12*pc(17,15) - 2*pc(12,20) + 3*pc(12,&
&13) + 3*pc(12,15) - 2*pc(19,20) + 3*pc(19,13) + 3*pc(19,15)

sh(1,15) = 4*pc(1,16) - pc(1,18) - pc(1,14)

sh(2,15) = 4*pc(2,16) - pc(2,18) - pc(2,14)

sh(3,15) = 4*pc(3,16) - pc(3,18) - pc(3,14)

sh(4,15) = 4*pc(4,16) - pc(4,18) - pc(4,14)

sh(5,15) = 4*pc(8,16) - pc(8,18) - pc(8,14)

sh(6,15) = 4*pc(9,16) - pc(9,18) - pc(9,14)

sh(7,15) = 4*pc(10,16) - pc(10,18) - pc(10,14)

sh(8,15) = 4*pc(5,16) - pc(5,18) - pc(5,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(9,15) = 8*pc(7,16) - 2*pc(7,18) - 2*pc(7,14) - 4*pc(5,16) + pc(5,18) + pc(5&
&,14) - 4*pc(6,16) + pc(6,18) + pc(6,14)

sh(10,15) = 4*pc(11,16) - pc(11,18) - pc(11,14)

sh(11,15) = 4*pc(13,16) - pc(13,18) - pc(13,14) - 4*pc(15,16) + pc(15,18) + pc&
&(15,14)

sh(12,15) = 4*pc(18,16) - pc(18,18) - pc(18,14) - 12*pc(14,16) + 3*pc(14,18) +&
& 3*pc(14,14)

sh(13,15) = 12*pc(12,16) - 3*pc(12,18) - 3*pc(12,14) - 4*pc(19,16) + pc(19,18)&
& + pc(19,14)

sh(14,15) = 8*pc(20,16) - 2*pc(20,18) - 2*pc(20,14) - 12*pc(13,16) + 3*pc(13,1&
&8) + 3*pc(13,14) - 12*pc(15,16) + 3*pc(15,18) + 3*pc(15,14)

sh(15,15) = 16*pc(16,16) - 4*pc(16,18) - 4*pc(16,14) - 4*pc(18,16) + pc(18,18)&
& + pc(18,14) - 4*pc(14,16) + pc(14,18) + pc(14,14)

sh(16,15) = 16*pc(17,16) - 4*pc(17,18) - 4*pc(17,14) - 4*pc(12,16) + pc(12,18)&
& + pc(12,14) - 4*pc(19,16) + pc(19,18) + pc(19,14)

sh(1,16) = 4*pc(1,17) - pc(1,12) - pc(1,19)

sh(2,16) = 4*pc(2,17) - pc(2,12) - pc(2,19)

sh(3,16) = 4*pc(3,17) - pc(3,12) - pc(3,19)

sh(4,16) = 4*pc(4,17) - pc(4,12) - pc(4,19)

sh(5,16) = 4*pc(8,17) - pc(8,12) - pc(8,19)

sh(6,16) = 4*pc(9,17) - pc(9,12) - pc(9,19)

sh(7,16) = 4*pc(10,17) - pc(10,12) - pc(10,19)

sh(8,16) = 4*pc(5,17) - pc(5,12) - pc(5,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(9,16) = 8*pc(7,17) - 2*pc(7,12) - 2*pc(7,19) - 4*pc(5,17) + pc(5,12) + pc(5&
&,19) - 4*pc(6,17) + pc(6,12) + pc(6,19)

sh(10,16) = 4*pc(11,17) - pc(11,12) - pc(11,19)

sh(11,16) = 4*pc(13,17) - pc(13,12) - pc(13,19) - 4*pc(15,17) + pc(15,12) + pc&
&(15,19)

sh(12,16) = 4*pc(18,17) - pc(18,12) - pc(18,19) - 12*pc(14,17) + 3*pc(14,12) +&
& 3*pc(14,19)

sh(13,16) = 12*pc(12,17) - 3*pc(12,12) - 3*pc(12,19) - 4*pc(19,17) + pc(19,12)&
& + pc(19,19)

sh(14,16) = 8*pc(20,17) - 2*pc(20,12) - 2*pc(20,19) - 12*pc(13,17) + 3*pc(13,1&
&2) + 3*pc(13,19) - 12*pc(15,17) + 3*pc(15,12) + 3*pc(15,19)

sh(15,16) = 16*pc(16,17) - 4*pc(16,12) - 4*pc(16,19) - 4*pc(18,17) + pc(18,12)&
& + pc(18,19) - 4*pc(14,17) + pc(14,12) + pc(14,19)

sh(16,16) = 16*pc(17,17) - 4*pc(17,12) - 4*pc(17,19) - 4*pc(12,17) + pc(12,12)&
& + pc(12,19) - 4*pc(19,17) + pc(19,12) + pc(19,19)


   end subroutine nuclear3CIntgNum

   subroutine momentum2CIntgAna(a1,a2,A,B,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   real (kind=double), dimension (20,20,3), intent(out) :: pc
   real (kind=double), dimension (16,16,3), intent(out) :: sh

   ! Define local variables.
   real (kind=double), dimension (20,20) :: pc_ol
   real (kind=double), dimension (3) :: P, PA, PB, d, preFactorMM
   real (kind=double) :: zeta, inv_2zeta, xi, preFactorOL
   real (kind=double) :: zeta_a_zeta, zeta_b_zeta

   ! Initialize local variables.
   zeta = a1 + a2
   zeta_a_zeta = a1/zeta
   zeta_b_zeta = a2/zeta
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   d = A - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d**2))
   preFactorMM(1) = preFactorOL*2.0d0*a1*PA(1) ! a or b arbitrary"
   preFactorMM(2) = preFactorOL*2.0d0*a1*PA(2) ! a or b arbitrary"
   preFactorMM(3) = preFactorOL*2.0d0*a1*PA(3) ! a or b arbitrary"

pc_ol(1,1) = preFactorOL

pc_ol(2,1) = PA(1)*preFactorOL

pc_ol(3,1) = PA(2)*preFactorOL

pc_ol(4,1) = PA(3)*preFactorOL

pc_ol(5,1) = PA(1)*pc_ol(2,1) + inv_2zeta*(preFactorOL)

pc_ol(6,1) = PA(2)*pc_ol(3,1) + inv_2zeta*(preFactorOL)

pc_ol(7,1) = PA(3)*pc_ol(4,1) + inv_2zeta*(preFactorOL)

pc_ol(8,1) = PA(2)*pc_ol(2,1)

pc_ol(9,1) = PA(3)*pc_ol(2,1)

pc_ol(10,1) = PA(3)*pc_ol(3,1)

pc_ol(11,1) = PA(3)*pc_ol(8,1)

pc_ol(12,1) = PA(2)*pc_ol(5,1)

pc_ol(13,1) = PA(3)*pc_ol(5,1)

pc_ol(14,1) = PA(1)*pc_ol(6,1)

pc_ol(15,1) = PA(3)*pc_ol(6,1)

pc_ol(16,1) = PA(1)*pc_ol(7,1)

pc_ol(17,1) = PA(2)*pc_ol(7,1)

pc_ol(18,1) = PA(1)*pc_ol(5,1) + inv_2zeta*(2*pc_ol(2,1))

pc_ol(19,1) = PA(2)*pc_ol(6,1) + inv_2zeta*(2*pc_ol(3,1))

pc_ol(20,1) = PA(3)*pc_ol(7,1) + inv_2zeta*(2*pc_ol(4,1))

pc_ol(1,2) = PB(1)*preFactorOL

pc_ol(2,2) = PA(1)*pc_ol(1,2) + inv_2zeta*(preFactorOL)

pc_ol(3,2) = PA(2)*pc_ol(1,2)

pc_ol(4,2) = PA(3)*pc_ol(1,2)

pc_ol(5,2) = PA(1)*pc_ol(2,2) + inv_2zeta*(pc_ol(1,2) + pc_ol(2,1))

pc_ol(6,2) = PA(2)*pc_ol(3,2) + inv_2zeta*(pc_ol(1,2))

pc_ol(7,2) = PA(3)*pc_ol(4,2) + inv_2zeta*(pc_ol(1,2))

pc_ol(8,2) = PA(2)*pc_ol(2,2)

pc_ol(9,2) = PA(3)*pc_ol(2,2)

pc_ol(10,2) = PA(3)*pc_ol(3,2)

pc_ol(11,2) = PA(3)*pc_ol(8,2)

pc_ol(12,2) = PA(2)*pc_ol(5,2)

pc_ol(13,2) = PA(3)*pc_ol(5,2)

pc_ol(14,2) = PA(1)*pc_ol(6,2) + inv_2zeta*(pc_ol(6,1))

pc_ol(15,2) = PA(3)*pc_ol(6,2)

pc_ol(16,2) = PA(1)*pc_ol(7,2) + inv_2zeta*(pc_ol(7,1))

pc_ol(17,2) = PA(2)*pc_ol(7,2)

pc_ol(18,2) = PA(1)*pc_ol(5,2) + inv_2zeta*(2*pc_ol(2,2) + pc_ol(5,1))

pc_ol(19,2) = PA(2)*pc_ol(6,2) + inv_2zeta*(2*pc_ol(3,2))

pc_ol(20,2) = PA(3)*pc_ol(7,2) + inv_2zeta*(2*pc_ol(4,2))

pc_ol(1,3) = PB(2)*preFactorOL

pc_ol(2,3) = PA(1)*pc_ol(1,3)

pc_ol(3,3) = PA(2)*pc_ol(1,3) + inv_2zeta*(preFactorOL)

pc_ol(4,3) = PA(3)*pc_ol(1,3)

pc_ol(5,3) = PA(1)*pc_ol(2,3) + inv_2zeta*(pc_ol(1,3))

pc_ol(6,3) = PA(2)*pc_ol(3,3) + inv_2zeta*(pc_ol(1,3) + pc_ol(3,1))

pc_ol(7,3) = PA(3)*pc_ol(4,3) + inv_2zeta*(pc_ol(1,3))

pc_ol(8,3) = PA(2)*pc_ol(2,3) + inv_2zeta*(pc_ol(2,1))

pc_ol(9,3) = PA(3)*pc_ol(2,3)

pc_ol(10,3) = PA(3)*pc_ol(3,3)

pc_ol(11,3) = PA(3)*pc_ol(8,3)

pc_ol(12,3) = PA(2)*pc_ol(5,3) + inv_2zeta*(pc_ol(5,1))

pc_ol(13,3) = PA(3)*pc_ol(5,3)

pc_ol(14,3) = PA(1)*pc_ol(6,3)

pc_ol(15,3) = PA(3)*pc_ol(6,3)

pc_ol(16,3) = PA(1)*pc_ol(7,3)

pc_ol(17,3) = PA(2)*pc_ol(7,3) + inv_2zeta*(pc_ol(7,1))

pc_ol(18,3) = PA(1)*pc_ol(5,3) + inv_2zeta*(2*pc_ol(2,3))

pc_ol(19,3) = PA(2)*pc_ol(6,3) + inv_2zeta*(2*pc_ol(3,3) + pc_ol(6,1))

pc_ol(20,3) = PA(3)*pc_ol(7,3) + inv_2zeta*(2*pc_ol(4,3))

pc_ol(1,4) = PB(3)*preFactorOL

pc_ol(2,4) = PA(1)*pc_ol(1,4)

pc_ol(3,4) = PA(2)*pc_ol(1,4)

pc_ol(4,4) = PA(3)*pc_ol(1,4) + inv_2zeta*(preFactorOL)

pc_ol(5,4) = PA(1)*pc_ol(2,4) + inv_2zeta*(pc_ol(1,4))

pc_ol(6,4) = PA(2)*pc_ol(3,4) + inv_2zeta*(pc_ol(1,4))

pc_ol(7,4) = PA(3)*pc_ol(4,4) + inv_2zeta*(pc_ol(1,4) + pc_ol(4,1))

pc_ol(8,4) = PA(2)*pc_ol(2,4)

pc_ol(9,4) = PA(3)*pc_ol(2,4) + inv_2zeta*(pc_ol(2,1))

pc_ol(10,4) = PA(3)*pc_ol(3,4) + inv_2zeta*(pc_ol(3,1))

pc_ol(11,4) = PA(3)*pc_ol(8,4) + inv_2zeta*(pc_ol(8,1))

pc_ol(12,4) = PA(2)*pc_ol(5,4)

pc_ol(13,4) = PA(3)*pc_ol(5,4) + inv_2zeta*(pc_ol(5,1))

pc_ol(14,4) = PA(1)*pc_ol(6,4)

pc_ol(15,4) = PA(3)*pc_ol(6,4) + inv_2zeta*(pc_ol(6,1))

pc_ol(16,4) = PA(1)*pc_ol(7,4)

pc_ol(17,4) = PA(2)*pc_ol(7,4)

pc_ol(18,4) = PA(1)*pc_ol(5,4) + inv_2zeta*(2*pc_ol(2,4))

pc_ol(19,4) = PA(2)*pc_ol(6,4) + inv_2zeta*(2*pc_ol(3,4))

pc_ol(20,4) = PA(3)*pc_ol(7,4) + inv_2zeta*(2*pc_ol(4,4) + pc_ol(7,1))

pc_ol(1,5) = PB(1)*pc_ol(1,2) + inv_2zeta*(preFactorOL)

pc_ol(2,5) = PA(1)*pc_ol(1,5) + inv_2zeta*(2*pc_ol(1,2))

pc_ol(3,5) = PA(2)*pc_ol(1,5)

pc_ol(4,5) = PA(3)*pc_ol(1,5)

pc_ol(5,5) = PA(1)*pc_ol(2,5) + inv_2zeta*(pc_ol(1,5) + 2*pc_ol(2,2))

pc_ol(6,5) = PA(2)*pc_ol(3,5) + inv_2zeta*(pc_ol(1,5))

pc_ol(7,5) = PA(3)*pc_ol(4,5) + inv_2zeta*(pc_ol(1,5))

pc_ol(8,5) = PA(2)*pc_ol(2,5)

pc_ol(9,5) = PA(3)*pc_ol(2,5)

pc_ol(10,5) = PA(3)*pc_ol(3,5)

pc_ol(11,5) = PA(3)*pc_ol(8,5)

pc_ol(12,5) = PA(2)*pc_ol(5,5)

pc_ol(13,5) = PA(3)*pc_ol(5,5)

pc_ol(14,5) = PA(1)*pc_ol(6,5) + inv_2zeta*(2*pc_ol(6,2))

pc_ol(15,5) = PA(3)*pc_ol(6,5)

pc_ol(16,5) = PA(1)*pc_ol(7,5) + inv_2zeta*(2*pc_ol(7,2))

pc_ol(17,5) = PA(2)*pc_ol(7,5)

pc_ol(18,5) = PA(1)*pc_ol(5,5) + inv_2zeta*(2*pc_ol(2,5) + 2*pc_ol(5,2))

pc_ol(19,5) = PA(2)*pc_ol(6,5) + inv_2zeta*(2*pc_ol(3,5))

pc_ol(20,5) = PA(3)*pc_ol(7,5) + inv_2zeta*(2*pc_ol(4,5))

pc_ol(1,6) = PB(2)*pc_ol(1,3) + inv_2zeta*(preFactorOL)

pc_ol(2,6) = PA(1)*pc_ol(1,6)

pc_ol(3,6) = PA(2)*pc_ol(1,6) + inv_2zeta*(2*pc_ol(1,3))

pc_ol(4,6) = PA(3)*pc_ol(1,6)

pc_ol(5,6) = PA(1)*pc_ol(2,6) + inv_2zeta*(pc_ol(1,6))

pc_ol(6,6) = PA(2)*pc_ol(3,6) + inv_2zeta*(pc_ol(1,6) + 2*pc_ol(3,3))

pc_ol(7,6) = PA(3)*pc_ol(4,6) + inv_2zeta*(pc_ol(1,6))

pc_ol(8,6) = PA(2)*pc_ol(2,6) + inv_2zeta*(2*pc_ol(2,3))

pc_ol(9,6) = PA(3)*pc_ol(2,6)

pc_ol(10,6) = PA(3)*pc_ol(3,6)

pc_ol(11,6) = PA(3)*pc_ol(8,6)

pc_ol(12,6) = PA(2)*pc_ol(5,6) + inv_2zeta*(2*pc_ol(5,3))

pc_ol(13,6) = PA(3)*pc_ol(5,6)

pc_ol(14,6) = PA(1)*pc_ol(6,6)

pc_ol(15,6) = PA(3)*pc_ol(6,6)

pc_ol(16,6) = PA(1)*pc_ol(7,6)

pc_ol(17,6) = PA(2)*pc_ol(7,6) + inv_2zeta*(2*pc_ol(7,3))

pc_ol(18,6) = PA(1)*pc_ol(5,6) + inv_2zeta*(2*pc_ol(2,6))

pc_ol(19,6) = PA(2)*pc_ol(6,6) + inv_2zeta*(2*pc_ol(3,6) + 2*pc_ol(6,3))

pc_ol(20,6) = PA(3)*pc_ol(7,6) + inv_2zeta*(2*pc_ol(4,6))

pc_ol(1,7) = PB(3)*pc_ol(1,4) + inv_2zeta*(preFactorOL)

pc_ol(2,7) = PA(1)*pc_ol(1,7)

pc_ol(3,7) = PA(2)*pc_ol(1,7)

pc_ol(4,7) = PA(3)*pc_ol(1,7) + inv_2zeta*(2*pc_ol(1,4))

pc_ol(5,7) = PA(1)*pc_ol(2,7) + inv_2zeta*(pc_ol(1,7))

pc_ol(6,7) = PA(2)*pc_ol(3,7) + inv_2zeta*(pc_ol(1,7))

pc_ol(7,7) = PA(3)*pc_ol(4,7) + inv_2zeta*(pc_ol(1,7) + 2*pc_ol(4,4))

pc_ol(8,7) = PA(2)*pc_ol(2,7)

pc_ol(9,7) = PA(3)*pc_ol(2,7) + inv_2zeta*(2*pc_ol(2,4))

pc_ol(10,7) = PA(3)*pc_ol(3,7) + inv_2zeta*(2*pc_ol(3,4))

pc_ol(11,7) = PA(3)*pc_ol(8,7) + inv_2zeta*(2*pc_ol(8,4))

pc_ol(12,7) = PA(2)*pc_ol(5,7)

pc_ol(13,7) = PA(3)*pc_ol(5,7) + inv_2zeta*(2*pc_ol(5,4))

pc_ol(14,7) = PA(1)*pc_ol(6,7)

pc_ol(15,7) = PA(3)*pc_ol(6,7) + inv_2zeta*(2*pc_ol(6,4))

pc_ol(16,7) = PA(1)*pc_ol(7,7)

pc_ol(17,7) = PA(2)*pc_ol(7,7)

pc_ol(18,7) = PA(1)*pc_ol(5,7) + inv_2zeta*(2*pc_ol(2,7))

pc_ol(19,7) = PA(2)*pc_ol(6,7) + inv_2zeta*(2*pc_ol(3,7))

pc_ol(20,7) = PA(3)*pc_ol(7,7) + inv_2zeta*(2*pc_ol(4,7) + 2*pc_ol(7,4))

pc_ol(1,8) = PB(2)*pc_ol(1,2)

pc_ol(2,8) = PA(1)*pc_ol(1,8) + inv_2zeta*(pc_ol(1,3))

pc_ol(3,8) = PA(2)*pc_ol(1,8) + inv_2zeta*(pc_ol(1,2))

pc_ol(4,8) = PA(3)*pc_ol(1,8)

pc_ol(5,8) = PA(1)*pc_ol(2,8) + inv_2zeta*(pc_ol(1,8) + pc_ol(2,3))

pc_ol(6,8) = PA(2)*pc_ol(3,8) + inv_2zeta*(pc_ol(1,8) + pc_ol(3,2))

pc_ol(7,8) = PA(3)*pc_ol(4,8) + inv_2zeta*(pc_ol(1,8))

pc_ol(8,8) = PA(2)*pc_ol(2,8) + inv_2zeta*(pc_ol(2,2))

pc_ol(9,8) = PA(3)*pc_ol(2,8)

pc_ol(10,8) = PA(3)*pc_ol(3,8)

pc_ol(11,8) = PA(3)*pc_ol(8,8)

pc_ol(12,8) = PA(2)*pc_ol(5,8) + inv_2zeta*(pc_ol(5,2))

pc_ol(13,8) = PA(3)*pc_ol(5,8)

pc_ol(14,8) = PA(1)*pc_ol(6,8) + inv_2zeta*(pc_ol(6,3))

pc_ol(15,8) = PA(3)*pc_ol(6,8)

pc_ol(16,8) = PA(1)*pc_ol(7,8) + inv_2zeta*(pc_ol(7,3))

pc_ol(17,8) = PA(2)*pc_ol(7,8) + inv_2zeta*(pc_ol(7,2))

pc_ol(18,8) = PA(1)*pc_ol(5,8) + inv_2zeta*(2*pc_ol(2,8) + pc_ol(5,3))

pc_ol(19,8) = PA(2)*pc_ol(6,8) + inv_2zeta*(2*pc_ol(3,8) + pc_ol(6,2))

pc_ol(20,8) = PA(3)*pc_ol(7,8) + inv_2zeta*(2*pc_ol(4,8))

pc_ol(1,9) = PB(3)*pc_ol(1,2)

pc_ol(2,9) = PA(1)*pc_ol(1,9) + inv_2zeta*(pc_ol(1,4))

pc_ol(3,9) = PA(2)*pc_ol(1,9)

pc_ol(4,9) = PA(3)*pc_ol(1,9) + inv_2zeta*(pc_ol(1,2))

pc_ol(5,9) = PA(1)*pc_ol(2,9) + inv_2zeta*(pc_ol(1,9) + pc_ol(2,4))

pc_ol(6,9) = PA(2)*pc_ol(3,9) + inv_2zeta*(pc_ol(1,9))

pc_ol(7,9) = PA(3)*pc_ol(4,9) + inv_2zeta*(pc_ol(1,9) + pc_ol(4,2))

pc_ol(8,9) = PA(2)*pc_ol(2,9)

pc_ol(9,9) = PA(3)*pc_ol(2,9) + inv_2zeta*(pc_ol(2,2))

pc_ol(10,9) = PA(3)*pc_ol(3,9) + inv_2zeta*(pc_ol(3,2))

pc_ol(11,9) = PA(3)*pc_ol(8,9) + inv_2zeta*(pc_ol(8,2))

pc_ol(12,9) = PA(2)*pc_ol(5,9)

pc_ol(13,9) = PA(3)*pc_ol(5,9) + inv_2zeta*(pc_ol(5,2))

pc_ol(14,9) = PA(1)*pc_ol(6,9) + inv_2zeta*(pc_ol(6,4))

pc_ol(15,9) = PA(3)*pc_ol(6,9) + inv_2zeta*(pc_ol(6,2))

pc_ol(16,9) = PA(1)*pc_ol(7,9) + inv_2zeta*(pc_ol(7,4))

pc_ol(17,9) = PA(2)*pc_ol(7,9)

pc_ol(18,9) = PA(1)*pc_ol(5,9) + inv_2zeta*(2*pc_ol(2,9) + pc_ol(5,4))

pc_ol(19,9) = PA(2)*pc_ol(6,9) + inv_2zeta*(2*pc_ol(3,9))

pc_ol(20,9) = PA(3)*pc_ol(7,9) + inv_2zeta*(2*pc_ol(4,9) + pc_ol(7,2))

pc_ol(1,10) = PB(3)*pc_ol(1,3)

pc_ol(2,10) = PA(1)*pc_ol(1,10)

pc_ol(3,10) = PA(2)*pc_ol(1,10) + inv_2zeta*(pc_ol(1,4))

pc_ol(4,10) = PA(3)*pc_ol(1,10) + inv_2zeta*(pc_ol(1,3))

pc_ol(5,10) = PA(1)*pc_ol(2,10) + inv_2zeta*(pc_ol(1,10))

pc_ol(6,10) = PA(2)*pc_ol(3,10) + inv_2zeta*(pc_ol(1,10) + pc_ol(3,4))

pc_ol(7,10) = PA(3)*pc_ol(4,10) + inv_2zeta*(pc_ol(1,10) + pc_ol(4,3))

pc_ol(8,10) = PA(2)*pc_ol(2,10) + inv_2zeta*(pc_ol(2,4))

pc_ol(9,10) = PA(3)*pc_ol(2,10) + inv_2zeta*(pc_ol(2,3))

pc_ol(10,10) = PA(3)*pc_ol(3,10) + inv_2zeta*(pc_ol(3,3))

pc_ol(11,10) = PA(3)*pc_ol(8,10) + inv_2zeta*(pc_ol(8,3))

pc_ol(12,10) = PA(2)*pc_ol(5,10) + inv_2zeta*(pc_ol(5,4))

pc_ol(13,10) = PA(3)*pc_ol(5,10) + inv_2zeta*(pc_ol(5,3))

pc_ol(14,10) = PA(1)*pc_ol(6,10)

pc_ol(15,10) = PA(3)*pc_ol(6,10) + inv_2zeta*(pc_ol(6,3))

pc_ol(16,10) = PA(1)*pc_ol(7,10)

pc_ol(17,10) = PA(2)*pc_ol(7,10) + inv_2zeta*(pc_ol(7,4))

pc_ol(18,10) = PA(1)*pc_ol(5,10) + inv_2zeta*(2*pc_ol(2,10))

pc_ol(19,10) = PA(2)*pc_ol(6,10) + inv_2zeta*(2*pc_ol(3,10) + pc_ol(6,4))

pc_ol(20,10) = PA(3)*pc_ol(7,10) + inv_2zeta*(2*pc_ol(4,10) + pc_ol(7,3))

pc_ol(1,11) = PB(3)*pc_ol(1,8)

pc_ol(2,11) = PA(1)*pc_ol(1,11) + inv_2zeta*(pc_ol(1,10))

pc_ol(3,11) = PA(2)*pc_ol(1,11) + inv_2zeta*(pc_ol(1,9))

pc_ol(4,11) = PA(3)*pc_ol(1,11) + inv_2zeta*(pc_ol(1,8))

pc_ol(5,11) = PA(1)*pc_ol(2,11) + inv_2zeta*(pc_ol(1,11) + pc_ol(2,10))

pc_ol(6,11) = PA(2)*pc_ol(3,11) + inv_2zeta*(pc_ol(1,11) + pc_ol(3,9))

pc_ol(7,11) = PA(3)*pc_ol(4,11) + inv_2zeta*(pc_ol(1,11) + pc_ol(4,8))

pc_ol(8,11) = PA(2)*pc_ol(2,11) + inv_2zeta*(pc_ol(2,9))

pc_ol(9,11) = PA(3)*pc_ol(2,11) + inv_2zeta*(pc_ol(2,8))

pc_ol(10,11) = PA(3)*pc_ol(3,11) + inv_2zeta*(pc_ol(3,8))

pc_ol(11,11) = PA(3)*pc_ol(8,11) + inv_2zeta*(pc_ol(8,8))

pc_ol(12,11) = PA(2)*pc_ol(5,11) + inv_2zeta*(pc_ol(5,9))

pc_ol(13,11) = PA(3)*pc_ol(5,11) + inv_2zeta*(pc_ol(5,8))

pc_ol(14,11) = PA(1)*pc_ol(6,11) + inv_2zeta*(pc_ol(6,10))

pc_ol(15,11) = PA(3)*pc_ol(6,11) + inv_2zeta*(pc_ol(6,8))

pc_ol(16,11) = PA(1)*pc_ol(7,11) + inv_2zeta*(pc_ol(7,10))

pc_ol(17,11) = PA(2)*pc_ol(7,11) + inv_2zeta*(pc_ol(7,9))

pc_ol(18,11) = PA(1)*pc_ol(5,11) + inv_2zeta*(2*pc_ol(2,11) + pc_ol(5,10))

pc_ol(19,11) = PA(2)*pc_ol(6,11) + inv_2zeta*(2*pc_ol(3,11) + pc_ol(6,9))

pc_ol(20,11) = PA(3)*pc_ol(7,11) + inv_2zeta*(2*pc_ol(4,11) + pc_ol(7,8))

pc_ol(1,12) = PB(2)*pc_ol(1,5)

pc_ol(2,12) = PA(1)*pc_ol(1,12) + inv_2zeta*(2*pc_ol(1,8))

pc_ol(3,12) = PA(2)*pc_ol(1,12) + inv_2zeta*(pc_ol(1,5))

pc_ol(4,12) = PA(3)*pc_ol(1,12)

pc_ol(5,12) = PA(1)*pc_ol(2,12) + inv_2zeta*(pc_ol(1,12) + 2*pc_ol(2,8))

pc_ol(6,12) = PA(2)*pc_ol(3,12) + inv_2zeta*(pc_ol(1,12) + pc_ol(3,5))

pc_ol(7,12) = PA(3)*pc_ol(4,12) + inv_2zeta*(pc_ol(1,12))

pc_ol(8,12) = PA(2)*pc_ol(2,12) + inv_2zeta*(pc_ol(2,5))

pc_ol(9,12) = PA(3)*pc_ol(2,12)

pc_ol(10,12) = PA(3)*pc_ol(3,12)

pc_ol(11,12) = PA(3)*pc_ol(8,12)

pc_ol(12,12) = PA(2)*pc_ol(5,12) + inv_2zeta*(pc_ol(5,5))

pc_ol(13,12) = PA(3)*pc_ol(5,12)

pc_ol(14,12) = PA(1)*pc_ol(6,12) + inv_2zeta*(2*pc_ol(6,8))

pc_ol(15,12) = PA(3)*pc_ol(6,12)

pc_ol(16,12) = PA(1)*pc_ol(7,12) + inv_2zeta*(2*pc_ol(7,8))

pc_ol(17,12) = PA(2)*pc_ol(7,12) + inv_2zeta*(pc_ol(7,5))

pc_ol(18,12) = PA(1)*pc_ol(5,12) + inv_2zeta*(2*pc_ol(2,12) + 2*pc_ol(5,8))

pc_ol(19,12) = PA(2)*pc_ol(6,12) + inv_2zeta*(2*pc_ol(3,12) + pc_ol(6,5))

pc_ol(20,12) = PA(3)*pc_ol(7,12) + inv_2zeta*(2*pc_ol(4,12))

pc_ol(1,13) = PB(3)*pc_ol(1,5)

pc_ol(2,13) = PA(1)*pc_ol(1,13) + inv_2zeta*(2*pc_ol(1,9))

pc_ol(3,13) = PA(2)*pc_ol(1,13)

pc_ol(4,13) = PA(3)*pc_ol(1,13) + inv_2zeta*(pc_ol(1,5))

pc_ol(5,13) = PA(1)*pc_ol(2,13) + inv_2zeta*(pc_ol(1,13) + 2*pc_ol(2,9))

pc_ol(6,13) = PA(2)*pc_ol(3,13) + inv_2zeta*(pc_ol(1,13))

pc_ol(7,13) = PA(3)*pc_ol(4,13) + inv_2zeta*(pc_ol(1,13) + pc_ol(4,5))

pc_ol(8,13) = PA(2)*pc_ol(2,13)

pc_ol(9,13) = PA(3)*pc_ol(2,13) + inv_2zeta*(pc_ol(2,5))

pc_ol(10,13) = PA(3)*pc_ol(3,13) + inv_2zeta*(pc_ol(3,5))

pc_ol(11,13) = PA(3)*pc_ol(8,13) + inv_2zeta*(pc_ol(8,5))

pc_ol(12,13) = PA(2)*pc_ol(5,13)

pc_ol(13,13) = PA(3)*pc_ol(5,13) + inv_2zeta*(pc_ol(5,5))

pc_ol(14,13) = PA(1)*pc_ol(6,13) + inv_2zeta*(2*pc_ol(6,9))

pc_ol(15,13) = PA(3)*pc_ol(6,13) + inv_2zeta*(pc_ol(6,5))

pc_ol(16,13) = PA(1)*pc_ol(7,13) + inv_2zeta*(2*pc_ol(7,9))

pc_ol(17,13) = PA(2)*pc_ol(7,13)

pc_ol(18,13) = PA(1)*pc_ol(5,13) + inv_2zeta*(2*pc_ol(2,13) + 2*pc_ol(5,9))

pc_ol(19,13) = PA(2)*pc_ol(6,13) + inv_2zeta*(2*pc_ol(3,13))

pc_ol(20,13) = PA(3)*pc_ol(7,13) + inv_2zeta*(2*pc_ol(4,13) + pc_ol(7,5))

pc_ol(1,14) = PB(1)*pc_ol(1,6)

pc_ol(2,14) = PA(1)*pc_ol(1,14) + inv_2zeta*(pc_ol(1,6))

pc_ol(3,14) = PA(2)*pc_ol(1,14) + inv_2zeta*(2*pc_ol(1,8))

pc_ol(4,14) = PA(3)*pc_ol(1,14)

pc_ol(5,14) = PA(1)*pc_ol(2,14) + inv_2zeta*(pc_ol(1,14) + pc_ol(2,6))

pc_ol(6,14) = PA(2)*pc_ol(3,14) + inv_2zeta*(pc_ol(1,14) + 2*pc_ol(3,8))

pc_ol(7,14) = PA(3)*pc_ol(4,14) + inv_2zeta*(pc_ol(1,14))

pc_ol(8,14) = PA(2)*pc_ol(2,14) + inv_2zeta*(2*pc_ol(2,8))

pc_ol(9,14) = PA(3)*pc_ol(2,14)

pc_ol(10,14) = PA(3)*pc_ol(3,14)

pc_ol(11,14) = PA(3)*pc_ol(8,14)

pc_ol(12,14) = PA(2)*pc_ol(5,14) + inv_2zeta*(2*pc_ol(5,8))

pc_ol(13,14) = PA(3)*pc_ol(5,14)

pc_ol(14,14) = PA(1)*pc_ol(6,14) + inv_2zeta*(pc_ol(6,6))

pc_ol(15,14) = PA(3)*pc_ol(6,14)

pc_ol(16,14) = PA(1)*pc_ol(7,14) + inv_2zeta*(pc_ol(7,6))

pc_ol(17,14) = PA(2)*pc_ol(7,14) + inv_2zeta*(2*pc_ol(7,8))

pc_ol(18,14) = PA(1)*pc_ol(5,14) + inv_2zeta*(2*pc_ol(2,14) + pc_ol(5,6))

pc_ol(19,14) = PA(2)*pc_ol(6,14) + inv_2zeta*(2*pc_ol(3,14) + 2*pc_ol(6,8))

pc_ol(20,14) = PA(3)*pc_ol(7,14) + inv_2zeta*(2*pc_ol(4,14))

pc_ol(1,15) = PB(3)*pc_ol(1,6)

pc_ol(2,15) = PA(1)*pc_ol(1,15)

pc_ol(3,15) = PA(2)*pc_ol(1,15) + inv_2zeta*(2*pc_ol(1,10))

pc_ol(4,15) = PA(3)*pc_ol(1,15) + inv_2zeta*(pc_ol(1,6))

pc_ol(5,15) = PA(1)*pc_ol(2,15) + inv_2zeta*(pc_ol(1,15))

pc_ol(6,15) = PA(2)*pc_ol(3,15) + inv_2zeta*(pc_ol(1,15) + 2*pc_ol(3,10))

pc_ol(7,15) = PA(3)*pc_ol(4,15) + inv_2zeta*(pc_ol(1,15) + pc_ol(4,6))

pc_ol(8,15) = PA(2)*pc_ol(2,15) + inv_2zeta*(2*pc_ol(2,10))

pc_ol(9,15) = PA(3)*pc_ol(2,15) + inv_2zeta*(pc_ol(2,6))

pc_ol(10,15) = PA(3)*pc_ol(3,15) + inv_2zeta*(pc_ol(3,6))

pc_ol(11,15) = PA(3)*pc_ol(8,15) + inv_2zeta*(pc_ol(8,6))

pc_ol(12,15) = PA(2)*pc_ol(5,15) + inv_2zeta*(2*pc_ol(5,10))

pc_ol(13,15) = PA(3)*pc_ol(5,15) + inv_2zeta*(pc_ol(5,6))

pc_ol(14,15) = PA(1)*pc_ol(6,15)

pc_ol(15,15) = PA(3)*pc_ol(6,15) + inv_2zeta*(pc_ol(6,6))

pc_ol(16,15) = PA(1)*pc_ol(7,15)

pc_ol(17,15) = PA(2)*pc_ol(7,15) + inv_2zeta*(2*pc_ol(7,10))

pc_ol(18,15) = PA(1)*pc_ol(5,15) + inv_2zeta*(2*pc_ol(2,15))

pc_ol(19,15) = PA(2)*pc_ol(6,15) + inv_2zeta*(2*pc_ol(3,15) + 2*pc_ol(6,10))

pc_ol(20,15) = PA(3)*pc_ol(7,15) + inv_2zeta*(2*pc_ol(4,15) + pc_ol(7,6))

pc_ol(1,16) = PB(1)*pc_ol(1,7)

pc_ol(2,16) = PA(1)*pc_ol(1,16) + inv_2zeta*(pc_ol(1,7))

pc_ol(3,16) = PA(2)*pc_ol(1,16)

pc_ol(4,16) = PA(3)*pc_ol(1,16) + inv_2zeta*(2*pc_ol(1,9))

pc_ol(5,16) = PA(1)*pc_ol(2,16) + inv_2zeta*(pc_ol(1,16) + pc_ol(2,7))

pc_ol(6,16) = PA(2)*pc_ol(3,16) + inv_2zeta*(pc_ol(1,16))

pc_ol(7,16) = PA(3)*pc_ol(4,16) + inv_2zeta*(pc_ol(1,16) + 2*pc_ol(4,9))

pc_ol(8,16) = PA(2)*pc_ol(2,16)

pc_ol(9,16) = PA(3)*pc_ol(2,16) + inv_2zeta*(2*pc_ol(2,9))

pc_ol(10,16) = PA(3)*pc_ol(3,16) + inv_2zeta*(2*pc_ol(3,9))

pc_ol(11,16) = PA(3)*pc_ol(8,16) + inv_2zeta*(2*pc_ol(8,9))

pc_ol(12,16) = PA(2)*pc_ol(5,16)

pc_ol(13,16) = PA(3)*pc_ol(5,16) + inv_2zeta*(2*pc_ol(5,9))

pc_ol(14,16) = PA(1)*pc_ol(6,16) + inv_2zeta*(pc_ol(6,7))

pc_ol(15,16) = PA(3)*pc_ol(6,16) + inv_2zeta*(2*pc_ol(6,9))

pc_ol(16,16) = PA(1)*pc_ol(7,16) + inv_2zeta*(pc_ol(7,7))

pc_ol(17,16) = PA(2)*pc_ol(7,16)

pc_ol(18,16) = PA(1)*pc_ol(5,16) + inv_2zeta*(2*pc_ol(2,16) + pc_ol(5,7))

pc_ol(19,16) = PA(2)*pc_ol(6,16) + inv_2zeta*(2*pc_ol(3,16))

pc_ol(20,16) = PA(3)*pc_ol(7,16) + inv_2zeta*(2*pc_ol(4,16) + 2*pc_ol(7,9))

pc_ol(1,17) = PB(2)*pc_ol(1,7)

pc_ol(2,17) = PA(1)*pc_ol(1,17)

pc_ol(3,17) = PA(2)*pc_ol(1,17) + inv_2zeta*(pc_ol(1,7))

pc_ol(4,17) = PA(3)*pc_ol(1,17) + inv_2zeta*(2*pc_ol(1,10))

pc_ol(5,17) = PA(1)*pc_ol(2,17) + inv_2zeta*(pc_ol(1,17))

pc_ol(6,17) = PA(2)*pc_ol(3,17) + inv_2zeta*(pc_ol(1,17) + pc_ol(3,7))

pc_ol(7,17) = PA(3)*pc_ol(4,17) + inv_2zeta*(pc_ol(1,17) + 2*pc_ol(4,10))

pc_ol(8,17) = PA(2)*pc_ol(2,17) + inv_2zeta*(pc_ol(2,7))

pc_ol(9,17) = PA(3)*pc_ol(2,17) + inv_2zeta*(2*pc_ol(2,10))

pc_ol(10,17) = PA(3)*pc_ol(3,17) + inv_2zeta*(2*pc_ol(3,10))

pc_ol(11,17) = PA(3)*pc_ol(8,17) + inv_2zeta*(2*pc_ol(8,10))

pc_ol(12,17) = PA(2)*pc_ol(5,17) + inv_2zeta*(pc_ol(5,7))

pc_ol(13,17) = PA(3)*pc_ol(5,17) + inv_2zeta*(2*pc_ol(5,10))

pc_ol(14,17) = PA(1)*pc_ol(6,17)

pc_ol(15,17) = PA(3)*pc_ol(6,17) + inv_2zeta*(2*pc_ol(6,10))

pc_ol(16,17) = PA(1)*pc_ol(7,17)

pc_ol(17,17) = PA(2)*pc_ol(7,17) + inv_2zeta*(pc_ol(7,7))

pc_ol(18,17) = PA(1)*pc_ol(5,17) + inv_2zeta*(2*pc_ol(2,17))

pc_ol(19,17) = PA(2)*pc_ol(6,17) + inv_2zeta*(2*pc_ol(3,17) + pc_ol(6,7))

pc_ol(20,17) = PA(3)*pc_ol(7,17) + inv_2zeta*(2*pc_ol(4,17) + 2*pc_ol(7,10))

pc_ol(1,18) = PB(1)*pc_ol(1,5) + inv_2zeta*(2*pc_ol(1,2))

pc_ol(2,18) = PA(1)*pc_ol(1,18) + inv_2zeta*(3*pc_ol(1,5))

pc_ol(3,18) = PA(2)*pc_ol(1,18)

pc_ol(4,18) = PA(3)*pc_ol(1,18)

pc_ol(5,18) = PA(1)*pc_ol(2,18) + inv_2zeta*(pc_ol(1,18) + 3*pc_ol(2,5))

pc_ol(6,18) = PA(2)*pc_ol(3,18) + inv_2zeta*(pc_ol(1,18))

pc_ol(7,18) = PA(3)*pc_ol(4,18) + inv_2zeta*(pc_ol(1,18))

pc_ol(8,18) = PA(2)*pc_ol(2,18)

pc_ol(9,18) = PA(3)*pc_ol(2,18)

pc_ol(10,18) = PA(3)*pc_ol(3,18)

pc_ol(11,18) = PA(3)*pc_ol(8,18)

pc_ol(12,18) = PA(2)*pc_ol(5,18)

pc_ol(13,18) = PA(3)*pc_ol(5,18)

pc_ol(14,18) = PA(1)*pc_ol(6,18) + inv_2zeta*(3*pc_ol(6,5))

pc_ol(15,18) = PA(3)*pc_ol(6,18)

pc_ol(16,18) = PA(1)*pc_ol(7,18) + inv_2zeta*(3*pc_ol(7,5))

pc_ol(17,18) = PA(2)*pc_ol(7,18)

pc_ol(18,18) = PA(1)*pc_ol(5,18) + inv_2zeta*(2*pc_ol(2,18) + 3*pc_ol(5,5))

pc_ol(19,18) = PA(2)*pc_ol(6,18) + inv_2zeta*(2*pc_ol(3,18))

pc_ol(20,18) = PA(3)*pc_ol(7,18) + inv_2zeta*(2*pc_ol(4,18))

pc_ol(1,19) = PB(2)*pc_ol(1,6) + inv_2zeta*(2*pc_ol(1,3))

pc_ol(2,19) = PA(1)*pc_ol(1,19)

pc_ol(3,19) = PA(2)*pc_ol(1,19) + inv_2zeta*(3*pc_ol(1,6))

pc_ol(4,19) = PA(3)*pc_ol(1,19)

pc_ol(5,19) = PA(1)*pc_ol(2,19) + inv_2zeta*(pc_ol(1,19))

pc_ol(6,19) = PA(2)*pc_ol(3,19) + inv_2zeta*(pc_ol(1,19) + 3*pc_ol(3,6))

pc_ol(7,19) = PA(3)*pc_ol(4,19) + inv_2zeta*(pc_ol(1,19))

pc_ol(8,19) = PA(2)*pc_ol(2,19) + inv_2zeta*(3*pc_ol(2,6))

pc_ol(9,19) = PA(3)*pc_ol(2,19)

pc_ol(10,19) = PA(3)*pc_ol(3,19)

pc_ol(11,19) = PA(3)*pc_ol(8,19)

pc_ol(12,19) = PA(2)*pc_ol(5,19) + inv_2zeta*(3*pc_ol(5,6))

pc_ol(13,19) = PA(3)*pc_ol(5,19)

pc_ol(14,19) = PA(1)*pc_ol(6,19)

pc_ol(15,19) = PA(3)*pc_ol(6,19)

pc_ol(16,19) = PA(1)*pc_ol(7,19)

pc_ol(17,19) = PA(2)*pc_ol(7,19) + inv_2zeta*(3*pc_ol(7,6))

pc_ol(18,19) = PA(1)*pc_ol(5,19) + inv_2zeta*(2*pc_ol(2,19))

pc_ol(19,19) = PA(2)*pc_ol(6,19) + inv_2zeta*(2*pc_ol(3,19) + 3*pc_ol(6,6))

pc_ol(20,19) = PA(3)*pc_ol(7,19) + inv_2zeta*(2*pc_ol(4,19))

pc_ol(1,20) = PB(3)*pc_ol(1,7) + inv_2zeta*(2*pc_ol(1,4))

pc_ol(2,20) = PA(1)*pc_ol(1,20)

pc_ol(3,20) = PA(2)*pc_ol(1,20)

pc_ol(4,20) = PA(3)*pc_ol(1,20) + inv_2zeta*(3*pc_ol(1,7))

pc_ol(5,20) = PA(1)*pc_ol(2,20) + inv_2zeta*(pc_ol(1,20))

pc_ol(6,20) = PA(2)*pc_ol(3,20) + inv_2zeta*(pc_ol(1,20))

pc_ol(7,20) = PA(3)*pc_ol(4,20) + inv_2zeta*(pc_ol(1,20) + 3*pc_ol(4,7))

pc_ol(8,20) = PA(2)*pc_ol(2,20)

pc_ol(9,20) = PA(3)*pc_ol(2,20) + inv_2zeta*(3*pc_ol(2,7))

pc_ol(10,20) = PA(3)*pc_ol(3,20) + inv_2zeta*(3*pc_ol(3,7))

pc_ol(11,20) = PA(3)*pc_ol(8,20) + inv_2zeta*(3*pc_ol(8,7))

pc_ol(12,20) = PA(2)*pc_ol(5,20)

pc_ol(13,20) = PA(3)*pc_ol(5,20) + inv_2zeta*(3*pc_ol(5,7))

pc_ol(14,20) = PA(1)*pc_ol(6,20)

pc_ol(15,20) = PA(3)*pc_ol(6,20) + inv_2zeta*(3*pc_ol(6,7))

pc_ol(16,20) = PA(1)*pc_ol(7,20)

pc_ol(17,20) = PA(2)*pc_ol(7,20)

pc_ol(18,20) = PA(1)*pc_ol(5,20) + inv_2zeta*(2*pc_ol(2,20))

pc_ol(19,20) = PA(2)*pc_ol(6,20) + inv_2zeta*(2*pc_ol(3,20))

pc_ol(20,20) = PA(3)*pc_ol(7,20) + inv_2zeta*(2*pc_ol(4,20) + 3*pc_ol(7,7))

pc(1,1,1) = preFactorMM(1)

pc(2,1,1) = PA(1)*preFactorMM(1) - zeta_b_zeta*preFactorOL

pc(3,1,1) = PA(2)*preFactorMM(1)

pc(4,1,1) = PA(3)*preFactorMM(1)

pc(5,1,1) = PA(1)*pc(2,1,1) + inv_2zeta*(preFactorMM(1)) - zeta_b_zeta*pc_ol(2&
&,1)

pc(6,1,1) = PA(2)*pc(3,1,1) + inv_2zeta*(preFactorMM(1))

pc(7,1,1) = PA(3)*pc(4,1,1) + inv_2zeta*(preFactorMM(1))

pc(8,1,1) = PA(2)*pc(2,1,1)

pc(9,1,1) = PA(3)*pc(2,1,1)

pc(10,1,1) = PA(3)*pc(3,1,1)

pc(11,1,1) = PA(3)*pc(8,1,1)

pc(12,1,1) = PA(2)*pc(5,1,1)

pc(13,1,1) = PA(3)*pc(5,1,1)

pc(14,1,1) = PA(1)*pc(6,1,1) - zeta_b_zeta*pc_ol(6,1)

pc(15,1,1) = PA(3)*pc(6,1,1)

pc(16,1,1) = PA(1)*pc(7,1,1) - zeta_b_zeta*pc_ol(7,1)

pc(17,1,1) = PA(2)*pc(7,1,1)

pc(18,1,1) = PA(1)*pc(5,1,1) + inv_2zeta*(2*pc(2,1,1)) - zeta_b_zeta*pc_ol(5,1&
&)

pc(19,1,1) = PA(2)*pc(6,1,1) + inv_2zeta*(2*pc(3,1,1))

pc(20,1,1) = PA(3)*pc(7,1,1) + inv_2zeta*(2*pc(4,1,1))

pc(1,2,1) = PB(1)*preFactorMM(1) + zeta_a_zeta*preFactorOL

pc(2,2,1) = PA(1)*pc(1,2,1) + inv_2zeta*(preFactorMM(1)) - zeta_b_zeta*pc_ol(1&
&,2)

pc(3,2,1) = PA(2)*pc(1,2,1)

pc(4,2,1) = PA(3)*pc(1,2,1)

pc(5,2,1) = PA(1)*pc(2,2,1) + inv_2zeta*(pc(1,2,1) + pc(2,1,1)) - zeta_b_zeta*&
&pc_ol(2,2)

pc(6,2,1) = PA(2)*pc(3,2,1) + inv_2zeta*(pc(1,2,1))

pc(7,2,1) = PA(3)*pc(4,2,1) + inv_2zeta*(pc(1,2,1))

pc(8,2,1) = PA(2)*pc(2,2,1)

pc(9,2,1) = PA(3)*pc(2,2,1)

pc(10,2,1) = PA(3)*pc(3,2,1)

pc(11,2,1) = PA(3)*pc(8,2,1)

pc(12,2,1) = PA(2)*pc(5,2,1)

pc(13,2,1) = PA(3)*pc(5,2,1)

pc(14,2,1) = PA(1)*pc(6,2,1) + inv_2zeta*(pc(6,1,1)) - zeta_b_zeta*pc_ol(6,2)

pc(15,2,1) = PA(3)*pc(6,2,1)

pc(16,2,1) = PA(1)*pc(7,2,1) + inv_2zeta*(pc(7,1,1)) - zeta_b_zeta*pc_ol(7,2)

pc(17,2,1) = PA(2)*pc(7,2,1)

pc(18,2,1) = PA(1)*pc(5,2,1) + inv_2zeta*(2*pc(2,2,1) + pc(5,1,1)) - zeta_b_ze&
&ta*pc_ol(5,2)

pc(19,2,1) = PA(2)*pc(6,2,1) + inv_2zeta*(2*pc(3,2,1))

pc(20,2,1) = PA(3)*pc(7,2,1) + inv_2zeta*(2*pc(4,2,1))

pc(1,3,1) = PB(2)*preFactorMM(1)

pc(2,3,1) = PA(1)*pc(1,3,1) - zeta_b_zeta*pc_ol(1,3)

pc(3,3,1) = PA(2)*pc(1,3,1) + inv_2zeta*(preFactorMM(1))

pc(4,3,1) = PA(3)*pc(1,3,1)

pc(5,3,1) = PA(1)*pc(2,3,1) + inv_2zeta*(pc(1,3,1)) - zeta_b_zeta*pc_ol(2,3)

pc(6,3,1) = PA(2)*pc(3,3,1) + inv_2zeta*(pc(1,3,1) + pc(3,1,1))

pc(7,3,1) = PA(3)*pc(4,3,1) + inv_2zeta*(pc(1,3,1))

pc(8,3,1) = PA(2)*pc(2,3,1) + inv_2zeta*(pc(2,1,1))

pc(9,3,1) = PA(3)*pc(2,3,1)

pc(10,3,1) = PA(3)*pc(3,3,1)

pc(11,3,1) = PA(3)*pc(8,3,1)

pc(12,3,1) = PA(2)*pc(5,3,1) + inv_2zeta*(pc(5,1,1))

pc(13,3,1) = PA(3)*pc(5,3,1)

pc(14,3,1) = PA(1)*pc(6,3,1) - zeta_b_zeta*pc_ol(6,3)

pc(15,3,1) = PA(3)*pc(6,3,1)

pc(16,3,1) = PA(1)*pc(7,3,1) - zeta_b_zeta*pc_ol(7,3)

pc(17,3,1) = PA(2)*pc(7,3,1) + inv_2zeta*(pc(7,1,1))

pc(18,3,1) = PA(1)*pc(5,3,1) + inv_2zeta*(2*pc(2,3,1)) - zeta_b_zeta*pc_ol(5,3&
&)

pc(19,3,1) = PA(2)*pc(6,3,1) + inv_2zeta*(2*pc(3,3,1) + pc(6,1,1))

pc(20,3,1) = PA(3)*pc(7,3,1) + inv_2zeta*(2*pc(4,3,1))

pc(1,4,1) = PB(3)*preFactorMM(1)

pc(2,4,1) = PA(1)*pc(1,4,1) - zeta_b_zeta*pc_ol(1,4)

pc(3,4,1) = PA(2)*pc(1,4,1)

pc(4,4,1) = PA(3)*pc(1,4,1) + inv_2zeta*(preFactorMM(1))

pc(5,4,1) = PA(1)*pc(2,4,1) + inv_2zeta*(pc(1,4,1)) - zeta_b_zeta*pc_ol(2,4)

pc(6,4,1) = PA(2)*pc(3,4,1) + inv_2zeta*(pc(1,4,1))

pc(7,4,1) = PA(3)*pc(4,4,1) + inv_2zeta*(pc(1,4,1) + pc(4,1,1))

pc(8,4,1) = PA(2)*pc(2,4,1)

pc(9,4,1) = PA(3)*pc(2,4,1) + inv_2zeta*(pc(2,1,1))

pc(10,4,1) = PA(3)*pc(3,4,1) + inv_2zeta*(pc(3,1,1))

pc(11,4,1) = PA(3)*pc(8,4,1) + inv_2zeta*(pc(8,1,1))

pc(12,4,1) = PA(2)*pc(5,4,1)

pc(13,4,1) = PA(3)*pc(5,4,1) + inv_2zeta*(pc(5,1,1))

pc(14,4,1) = PA(1)*pc(6,4,1) - zeta_b_zeta*pc_ol(6,4)

pc(15,4,1) = PA(3)*pc(6,4,1) + inv_2zeta*(pc(6,1,1))

pc(16,4,1) = PA(1)*pc(7,4,1) - zeta_b_zeta*pc_ol(7,4)

pc(17,4,1) = PA(2)*pc(7,4,1)

pc(18,4,1) = PA(1)*pc(5,4,1) + inv_2zeta*(2*pc(2,4,1)) - zeta_b_zeta*pc_ol(5,4&
&)

pc(19,4,1) = PA(2)*pc(6,4,1) + inv_2zeta*(2*pc(3,4,1))

pc(20,4,1) = PA(3)*pc(7,4,1) + inv_2zeta*(2*pc(4,4,1) + pc(7,1,1))

pc(1,5,1) = PB(1)*pc(1,2,1) + inv_2zeta*(preFactorMM(1)) + zeta_a_zeta*pc_ol(1&
&,2)

pc(2,5,1) = PA(1)*pc(1,5,1) + inv_2zeta*(2*pc(1,2,1)) - zeta_b_zeta*pc_ol(1,5)

pc(3,5,1) = PA(2)*pc(1,5,1)

pc(4,5,1) = PA(3)*pc(1,5,1)

pc(5,5,1) = PA(1)*pc(2,5,1) + inv_2zeta*(pc(1,5,1) + 2*pc(2,2,1)) - zeta_b_zet&
&a*pc_ol(2,5)

pc(6,5,1) = PA(2)*pc(3,5,1) + inv_2zeta*(pc(1,5,1))

pc(7,5,1) = PA(3)*pc(4,5,1) + inv_2zeta*(pc(1,5,1))

pc(8,5,1) = PA(2)*pc(2,5,1)

pc(9,5,1) = PA(3)*pc(2,5,1)

pc(10,5,1) = PA(3)*pc(3,5,1)

pc(11,5,1) = PA(3)*pc(8,5,1)

pc(12,5,1) = PA(2)*pc(5,5,1)

pc(13,5,1) = PA(3)*pc(5,5,1)

pc(14,5,1) = PA(1)*pc(6,5,1) + inv_2zeta*(2*pc(6,2,1)) - zeta_b_zeta*pc_ol(6,5&
&)

pc(15,5,1) = PA(3)*pc(6,5,1)

pc(16,5,1) = PA(1)*pc(7,5,1) + inv_2zeta*(2*pc(7,2,1)) - zeta_b_zeta*pc_ol(7,5&
&)

pc(17,5,1) = PA(2)*pc(7,5,1)

pc(18,5,1) = PA(1)*pc(5,5,1) + inv_2zeta*(2*pc(2,5,1) + 2*pc(5,2,1)) - zeta_b_&
&zeta*pc_ol(5,5)

pc(19,5,1) = PA(2)*pc(6,5,1) + inv_2zeta*(2*pc(3,5,1))

pc(20,5,1) = PA(3)*pc(7,5,1) + inv_2zeta*(2*pc(4,5,1))

pc(1,6,1) = PB(2)*pc(1,3,1) + inv_2zeta*(preFactorMM(1))

pc(2,6,1) = PA(1)*pc(1,6,1) - zeta_b_zeta*pc_ol(1,6)

pc(3,6,1) = PA(2)*pc(1,6,1) + inv_2zeta*(2*pc(1,3,1))

pc(4,6,1) = PA(3)*pc(1,6,1)

pc(5,6,1) = PA(1)*pc(2,6,1) + inv_2zeta*(pc(1,6,1)) - zeta_b_zeta*pc_ol(2,6)

pc(6,6,1) = PA(2)*pc(3,6,1) + inv_2zeta*(pc(1,6,1) + 2*pc(3,3,1))

pc(7,6,1) = PA(3)*pc(4,6,1) + inv_2zeta*(pc(1,6,1))

pc(8,6,1) = PA(2)*pc(2,6,1) + inv_2zeta*(2*pc(2,3,1))

pc(9,6,1) = PA(3)*pc(2,6,1)

pc(10,6,1) = PA(3)*pc(3,6,1)

pc(11,6,1) = PA(3)*pc(8,6,1)

pc(12,6,1) = PA(2)*pc(5,6,1) + inv_2zeta*(2*pc(5,3,1))

pc(13,6,1) = PA(3)*pc(5,6,1)

pc(14,6,1) = PA(1)*pc(6,6,1) - zeta_b_zeta*pc_ol(6,6)

pc(15,6,1) = PA(3)*pc(6,6,1)

pc(16,6,1) = PA(1)*pc(7,6,1) - zeta_b_zeta*pc_ol(7,6)

pc(17,6,1) = PA(2)*pc(7,6,1) + inv_2zeta*(2*pc(7,3,1))

pc(18,6,1) = PA(1)*pc(5,6,1) + inv_2zeta*(2*pc(2,6,1)) - zeta_b_zeta*pc_ol(5,6&
&)

pc(19,6,1) = PA(2)*pc(6,6,1) + inv_2zeta*(2*pc(3,6,1) + 2*pc(6,3,1))

pc(20,6,1) = PA(3)*pc(7,6,1) + inv_2zeta*(2*pc(4,6,1))

pc(1,7,1) = PB(3)*pc(1,4,1) + inv_2zeta*(preFactorMM(1))

pc(2,7,1) = PA(1)*pc(1,7,1) - zeta_b_zeta*pc_ol(1,7)

pc(3,7,1) = PA(2)*pc(1,7,1)

pc(4,7,1) = PA(3)*pc(1,7,1) + inv_2zeta*(2*pc(1,4,1))

pc(5,7,1) = PA(1)*pc(2,7,1) + inv_2zeta*(pc(1,7,1)) - zeta_b_zeta*pc_ol(2,7)

pc(6,7,1) = PA(2)*pc(3,7,1) + inv_2zeta*(pc(1,7,1))

pc(7,7,1) = PA(3)*pc(4,7,1) + inv_2zeta*(pc(1,7,1) + 2*pc(4,4,1))

pc(8,7,1) = PA(2)*pc(2,7,1)

pc(9,7,1) = PA(3)*pc(2,7,1) + inv_2zeta*(2*pc(2,4,1))

pc(10,7,1) = PA(3)*pc(3,7,1) + inv_2zeta*(2*pc(3,4,1))

pc(11,7,1) = PA(3)*pc(8,7,1) + inv_2zeta*(2*pc(8,4,1))

pc(12,7,1) = PA(2)*pc(5,7,1)

pc(13,7,1) = PA(3)*pc(5,7,1) + inv_2zeta*(2*pc(5,4,1))

pc(14,7,1) = PA(1)*pc(6,7,1) - zeta_b_zeta*pc_ol(6,7)

pc(15,7,1) = PA(3)*pc(6,7,1) + inv_2zeta*(2*pc(6,4,1))

pc(16,7,1) = PA(1)*pc(7,7,1) - zeta_b_zeta*pc_ol(7,7)

pc(17,7,1) = PA(2)*pc(7,7,1)

pc(18,7,1) = PA(1)*pc(5,7,1) + inv_2zeta*(2*pc(2,7,1)) - zeta_b_zeta*pc_ol(5,7&
&)

pc(19,7,1) = PA(2)*pc(6,7,1) + inv_2zeta*(2*pc(3,7,1))

pc(20,7,1) = PA(3)*pc(7,7,1) + inv_2zeta*(2*pc(4,7,1) + 2*pc(7,4,1))

pc(1,8,1) = PB(2)*pc(1,2,1)

pc(2,8,1) = PA(1)*pc(1,8,1) + inv_2zeta*(pc(1,3,1)) - zeta_b_zeta*pc_ol(1,8)

pc(3,8,1) = PA(2)*pc(1,8,1) + inv_2zeta*(pc(1,2,1))

pc(4,8,1) = PA(3)*pc(1,8,1)

pc(5,8,1) = PA(1)*pc(2,8,1) + inv_2zeta*(pc(1,8,1) + pc(2,3,1)) - zeta_b_zeta*&
&pc_ol(2,8)

pc(6,8,1) = PA(2)*pc(3,8,1) + inv_2zeta*(pc(1,8,1) + pc(3,2,1))

pc(7,8,1) = PA(3)*pc(4,8,1) + inv_2zeta*(pc(1,8,1))

pc(8,8,1) = PA(2)*pc(2,8,1) + inv_2zeta*(pc(2,2,1))

pc(9,8,1) = PA(3)*pc(2,8,1)

pc(10,8,1) = PA(3)*pc(3,8,1)

pc(11,8,1) = PA(3)*pc(8,8,1)

pc(12,8,1) = PA(2)*pc(5,8,1) + inv_2zeta*(pc(5,2,1))

pc(13,8,1) = PA(3)*pc(5,8,1)

pc(14,8,1) = PA(1)*pc(6,8,1) + inv_2zeta*(pc(6,3,1)) - zeta_b_zeta*pc_ol(6,8)

pc(15,8,1) = PA(3)*pc(6,8,1)

pc(16,8,1) = PA(1)*pc(7,8,1) + inv_2zeta*(pc(7,3,1)) - zeta_b_zeta*pc_ol(7,8)

pc(17,8,1) = PA(2)*pc(7,8,1) + inv_2zeta*(pc(7,2,1))

pc(18,8,1) = PA(1)*pc(5,8,1) + inv_2zeta*(2*pc(2,8,1) + pc(5,3,1)) - zeta_b_ze&
&ta*pc_ol(5,8)

pc(19,8,1) = PA(2)*pc(6,8,1) + inv_2zeta*(2*pc(3,8,1) + pc(6,2,1))

pc(20,8,1) = PA(3)*pc(7,8,1) + inv_2zeta*(2*pc(4,8,1))

pc(1,9,1) = PB(3)*pc(1,2,1)

pc(2,9,1) = PA(1)*pc(1,9,1) + inv_2zeta*(pc(1,4,1)) - zeta_b_zeta*pc_ol(1,9)

pc(3,9,1) = PA(2)*pc(1,9,1)

pc(4,9,1) = PA(3)*pc(1,9,1) + inv_2zeta*(pc(1,2,1))

pc(5,9,1) = PA(1)*pc(2,9,1) + inv_2zeta*(pc(1,9,1) + pc(2,4,1)) - zeta_b_zeta*&
&pc_ol(2,9)

pc(6,9,1) = PA(2)*pc(3,9,1) + inv_2zeta*(pc(1,9,1))

pc(7,9,1) = PA(3)*pc(4,9,1) + inv_2zeta*(pc(1,9,1) + pc(4,2,1))

pc(8,9,1) = PA(2)*pc(2,9,1)

pc(9,9,1) = PA(3)*pc(2,9,1) + inv_2zeta*(pc(2,2,1))

pc(10,9,1) = PA(3)*pc(3,9,1) + inv_2zeta*(pc(3,2,1))

pc(11,9,1) = PA(3)*pc(8,9,1) + inv_2zeta*(pc(8,2,1))

pc(12,9,1) = PA(2)*pc(5,9,1)

pc(13,9,1) = PA(3)*pc(5,9,1) + inv_2zeta*(pc(5,2,1))

pc(14,9,1) = PA(1)*pc(6,9,1) + inv_2zeta*(pc(6,4,1)) - zeta_b_zeta*pc_ol(6,9)

pc(15,9,1) = PA(3)*pc(6,9,1) + inv_2zeta*(pc(6,2,1))

pc(16,9,1) = PA(1)*pc(7,9,1) + inv_2zeta*(pc(7,4,1)) - zeta_b_zeta*pc_ol(7,9)

pc(17,9,1) = PA(2)*pc(7,9,1)

pc(18,9,1) = PA(1)*pc(5,9,1) + inv_2zeta*(2*pc(2,9,1) + pc(5,4,1)) - zeta_b_ze&
&ta*pc_ol(5,9)

pc(19,9,1) = PA(2)*pc(6,9,1) + inv_2zeta*(2*pc(3,9,1))

pc(20,9,1) = PA(3)*pc(7,9,1) + inv_2zeta*(2*pc(4,9,1) + pc(7,2,1))

pc(1,10,1) = PB(3)*pc(1,3,1)

pc(2,10,1) = PA(1)*pc(1,10,1) - zeta_b_zeta*pc_ol(1,10)

pc(3,10,1) = PA(2)*pc(1,10,1) + inv_2zeta*(pc(1,4,1))

pc(4,10,1) = PA(3)*pc(1,10,1) + inv_2zeta*(pc(1,3,1))

pc(5,10,1) = PA(1)*pc(2,10,1) + inv_2zeta*(pc(1,10,1)) - zeta_b_zeta*pc_ol(2,1&
&0)

pc(6,10,1) = PA(2)*pc(3,10,1) + inv_2zeta*(pc(1,10,1) + pc(3,4,1))

pc(7,10,1) = PA(3)*pc(4,10,1) + inv_2zeta*(pc(1,10,1) + pc(4,3,1))

pc(8,10,1) = PA(2)*pc(2,10,1) + inv_2zeta*(pc(2,4,1))

pc(9,10,1) = PA(3)*pc(2,10,1) + inv_2zeta*(pc(2,3,1))

pc(10,10,1) = PA(3)*pc(3,10,1) + inv_2zeta*(pc(3,3,1))

pc(11,10,1) = PA(3)*pc(8,10,1) + inv_2zeta*(pc(8,3,1))

pc(12,10,1) = PA(2)*pc(5,10,1) + inv_2zeta*(pc(5,4,1))

pc(13,10,1) = PA(3)*pc(5,10,1) + inv_2zeta*(pc(5,3,1))

pc(14,10,1) = PA(1)*pc(6,10,1) - zeta_b_zeta*pc_ol(6,10)

pc(15,10,1) = PA(3)*pc(6,10,1) + inv_2zeta*(pc(6,3,1))

pc(16,10,1) = PA(1)*pc(7,10,1) - zeta_b_zeta*pc_ol(7,10)

pc(17,10,1) = PA(2)*pc(7,10,1) + inv_2zeta*(pc(7,4,1))

pc(18,10,1) = PA(1)*pc(5,10,1) + inv_2zeta*(2*pc(2,10,1)) - zeta_b_zeta*pc_ol(&
&5,10)

pc(19,10,1) = PA(2)*pc(6,10,1) + inv_2zeta*(2*pc(3,10,1) + pc(6,4,1))

pc(20,10,1) = PA(3)*pc(7,10,1) + inv_2zeta*(2*pc(4,10,1) + pc(7,3,1))

pc(1,11,1) = PB(3)*pc(1,8,1)

pc(2,11,1) = PA(1)*pc(1,11,1) + inv_2zeta*(pc(1,10,1)) - zeta_b_zeta*pc_ol(1,1&
&1)

pc(3,11,1) = PA(2)*pc(1,11,1) + inv_2zeta*(pc(1,9,1))

pc(4,11,1) = PA(3)*pc(1,11,1) + inv_2zeta*(pc(1,8,1))

pc(5,11,1) = PA(1)*pc(2,11,1) + inv_2zeta*(pc(1,11,1) + pc(2,10,1)) - zeta_b_z&
&eta*pc_ol(2,11)

pc(6,11,1) = PA(2)*pc(3,11,1) + inv_2zeta*(pc(1,11,1) + pc(3,9,1))

pc(7,11,1) = PA(3)*pc(4,11,1) + inv_2zeta*(pc(1,11,1) + pc(4,8,1))

pc(8,11,1) = PA(2)*pc(2,11,1) + inv_2zeta*(pc(2,9,1))

pc(9,11,1) = PA(3)*pc(2,11,1) + inv_2zeta*(pc(2,8,1))

pc(10,11,1) = PA(3)*pc(3,11,1) + inv_2zeta*(pc(3,8,1))

pc(11,11,1) = PA(3)*pc(8,11,1) + inv_2zeta*(pc(8,8,1))

pc(12,11,1) = PA(2)*pc(5,11,1) + inv_2zeta*(pc(5,9,1))

pc(13,11,1) = PA(3)*pc(5,11,1) + inv_2zeta*(pc(5,8,1))

pc(14,11,1) = PA(1)*pc(6,11,1) + inv_2zeta*(pc(6,10,1)) - zeta_b_zeta*pc_ol(6,&
&11)

pc(15,11,1) = PA(3)*pc(6,11,1) + inv_2zeta*(pc(6,8,1))

pc(16,11,1) = PA(1)*pc(7,11,1) + inv_2zeta*(pc(7,10,1)) - zeta_b_zeta*pc_ol(7,&
&11)

pc(17,11,1) = PA(2)*pc(7,11,1) + inv_2zeta*(pc(7,9,1))

pc(18,11,1) = PA(1)*pc(5,11,1) + inv_2zeta*(2*pc(2,11,1) + pc(5,10,1)) - zeta_&
&b_zeta*pc_ol(5,11)

pc(19,11,1) = PA(2)*pc(6,11,1) + inv_2zeta*(2*pc(3,11,1) + pc(6,9,1))

pc(20,11,1) = PA(3)*pc(7,11,1) + inv_2zeta*(2*pc(4,11,1) + pc(7,8,1))

pc(1,12,1) = PB(2)*pc(1,5,1)

pc(2,12,1) = PA(1)*pc(1,12,1) + inv_2zeta*(2*pc(1,8,1)) - zeta_b_zeta*pc_ol(1,&
&12)

pc(3,12,1) = PA(2)*pc(1,12,1) + inv_2zeta*(pc(1,5,1))

pc(4,12,1) = PA(3)*pc(1,12,1)

pc(5,12,1) = PA(1)*pc(2,12,1) + inv_2zeta*(pc(1,12,1) + 2*pc(2,8,1)) - zeta_b_&
&zeta*pc_ol(2,12)

pc(6,12,1) = PA(2)*pc(3,12,1) + inv_2zeta*(pc(1,12,1) + pc(3,5,1))

pc(7,12,1) = PA(3)*pc(4,12,1) + inv_2zeta*(pc(1,12,1))

pc(8,12,1) = PA(2)*pc(2,12,1) + inv_2zeta*(pc(2,5,1))

pc(9,12,1) = PA(3)*pc(2,12,1)

pc(10,12,1) = PA(3)*pc(3,12,1)

pc(11,12,1) = PA(3)*pc(8,12,1)

pc(12,12,1) = PA(2)*pc(5,12,1) + inv_2zeta*(pc(5,5,1))

pc(13,12,1) = PA(3)*pc(5,12,1)

pc(14,12,1) = PA(1)*pc(6,12,1) + inv_2zeta*(2*pc(6,8,1)) - zeta_b_zeta*pc_ol(6&
&,12)

pc(15,12,1) = PA(3)*pc(6,12,1)

pc(16,12,1) = PA(1)*pc(7,12,1) + inv_2zeta*(2*pc(7,8,1)) - zeta_b_zeta*pc_ol(7&
&,12)

pc(17,12,1) = PA(2)*pc(7,12,1) + inv_2zeta*(pc(7,5,1))

pc(18,12,1) = PA(1)*pc(5,12,1) + inv_2zeta*(2*pc(2,12,1) + 2*pc(5,8,1)) - zeta&
&_b_zeta*pc_ol(5,12)

pc(19,12,1) = PA(2)*pc(6,12,1) + inv_2zeta*(2*pc(3,12,1) + pc(6,5,1))

pc(20,12,1) = PA(3)*pc(7,12,1) + inv_2zeta*(2*pc(4,12,1))

pc(1,13,1) = PB(3)*pc(1,5,1)

pc(2,13,1) = PA(1)*pc(1,13,1) + inv_2zeta*(2*pc(1,9,1)) - zeta_b_zeta*pc_ol(1,&
&13)

pc(3,13,1) = PA(2)*pc(1,13,1)

pc(4,13,1) = PA(3)*pc(1,13,1) + inv_2zeta*(pc(1,5,1))

pc(5,13,1) = PA(1)*pc(2,13,1) + inv_2zeta*(pc(1,13,1) + 2*pc(2,9,1)) - zeta_b_&
&zeta*pc_ol(2,13)

pc(6,13,1) = PA(2)*pc(3,13,1) + inv_2zeta*(pc(1,13,1))

pc(7,13,1) = PA(3)*pc(4,13,1) + inv_2zeta*(pc(1,13,1) + pc(4,5,1))

pc(8,13,1) = PA(2)*pc(2,13,1)

pc(9,13,1) = PA(3)*pc(2,13,1) + inv_2zeta*(pc(2,5,1))

pc(10,13,1) = PA(3)*pc(3,13,1) + inv_2zeta*(pc(3,5,1))

pc(11,13,1) = PA(3)*pc(8,13,1) + inv_2zeta*(pc(8,5,1))

pc(12,13,1) = PA(2)*pc(5,13,1)

pc(13,13,1) = PA(3)*pc(5,13,1) + inv_2zeta*(pc(5,5,1))

pc(14,13,1) = PA(1)*pc(6,13,1) + inv_2zeta*(2*pc(6,9,1)) - zeta_b_zeta*pc_ol(6&
&,13)

pc(15,13,1) = PA(3)*pc(6,13,1) + inv_2zeta*(pc(6,5,1))

pc(16,13,1) = PA(1)*pc(7,13,1) + inv_2zeta*(2*pc(7,9,1)) - zeta_b_zeta*pc_ol(7&
&,13)

pc(17,13,1) = PA(2)*pc(7,13,1)

pc(18,13,1) = PA(1)*pc(5,13,1) + inv_2zeta*(2*pc(2,13,1) + 2*pc(5,9,1)) - zeta&
&_b_zeta*pc_ol(5,13)

pc(19,13,1) = PA(2)*pc(6,13,1) + inv_2zeta*(2*pc(3,13,1))

pc(20,13,1) = PA(3)*pc(7,13,1) + inv_2zeta*(2*pc(4,13,1) + pc(7,5,1))

pc(1,14,1) = PB(1)*pc(1,6,1) + zeta_a_zeta*pc_ol(1,6)

pc(2,14,1) = PA(1)*pc(1,14,1) + inv_2zeta*(pc(1,6,1)) - zeta_b_zeta*pc_ol(1,14&
&)

pc(3,14,1) = PA(2)*pc(1,14,1) + inv_2zeta*(2*pc(1,8,1))

pc(4,14,1) = PA(3)*pc(1,14,1)

pc(5,14,1) = PA(1)*pc(2,14,1) + inv_2zeta*(pc(1,14,1) + pc(2,6,1)) - zeta_b_ze&
&ta*pc_ol(2,14)

pc(6,14,1) = PA(2)*pc(3,14,1) + inv_2zeta*(pc(1,14,1) + 2*pc(3,8,1))

pc(7,14,1) = PA(3)*pc(4,14,1) + inv_2zeta*(pc(1,14,1))

pc(8,14,1) = PA(2)*pc(2,14,1) + inv_2zeta*(2*pc(2,8,1))

pc(9,14,1) = PA(3)*pc(2,14,1)

pc(10,14,1) = PA(3)*pc(3,14,1)

pc(11,14,1) = PA(3)*pc(8,14,1)

pc(12,14,1) = PA(2)*pc(5,14,1) + inv_2zeta*(2*pc(5,8,1))

pc(13,14,1) = PA(3)*pc(5,14,1)

pc(14,14,1) = PA(1)*pc(6,14,1) + inv_2zeta*(pc(6,6,1)) - zeta_b_zeta*pc_ol(6,1&
&4)

pc(15,14,1) = PA(3)*pc(6,14,1)

pc(16,14,1) = PA(1)*pc(7,14,1) + inv_2zeta*(pc(7,6,1)) - zeta_b_zeta*pc_ol(7,1&
&4)

pc(17,14,1) = PA(2)*pc(7,14,1) + inv_2zeta*(2*pc(7,8,1))

pc(18,14,1) = PA(1)*pc(5,14,1) + inv_2zeta*(2*pc(2,14,1) + pc(5,6,1)) - zeta_b&
&_zeta*pc_ol(5,14)

pc(19,14,1) = PA(2)*pc(6,14,1) + inv_2zeta*(2*pc(3,14,1) + 2*pc(6,8,1))

pc(20,14,1) = PA(3)*pc(7,14,1) + inv_2zeta*(2*pc(4,14,1))

pc(1,15,1) = PB(3)*pc(1,6,1)

pc(2,15,1) = PA(1)*pc(1,15,1) - zeta_b_zeta*pc_ol(1,15)

pc(3,15,1) = PA(2)*pc(1,15,1) + inv_2zeta*(2*pc(1,10,1))

pc(4,15,1) = PA(3)*pc(1,15,1) + inv_2zeta*(pc(1,6,1))

pc(5,15,1) = PA(1)*pc(2,15,1) + inv_2zeta*(pc(1,15,1)) - zeta_b_zeta*pc_ol(2,1&
&5)

pc(6,15,1) = PA(2)*pc(3,15,1) + inv_2zeta*(pc(1,15,1) + 2*pc(3,10,1))

pc(7,15,1) = PA(3)*pc(4,15,1) + inv_2zeta*(pc(1,15,1) + pc(4,6,1))

pc(8,15,1) = PA(2)*pc(2,15,1) + inv_2zeta*(2*pc(2,10,1))

pc(9,15,1) = PA(3)*pc(2,15,1) + inv_2zeta*(pc(2,6,1))

pc(10,15,1) = PA(3)*pc(3,15,1) + inv_2zeta*(pc(3,6,1))

pc(11,15,1) = PA(3)*pc(8,15,1) + inv_2zeta*(pc(8,6,1))

pc(12,15,1) = PA(2)*pc(5,15,1) + inv_2zeta*(2*pc(5,10,1))

pc(13,15,1) = PA(3)*pc(5,15,1) + inv_2zeta*(pc(5,6,1))

pc(14,15,1) = PA(1)*pc(6,15,1) - zeta_b_zeta*pc_ol(6,15)

pc(15,15,1) = PA(3)*pc(6,15,1) + inv_2zeta*(pc(6,6,1))

pc(16,15,1) = PA(1)*pc(7,15,1) - zeta_b_zeta*pc_ol(7,15)

pc(17,15,1) = PA(2)*pc(7,15,1) + inv_2zeta*(2*pc(7,10,1))

pc(18,15,1) = PA(1)*pc(5,15,1) + inv_2zeta*(2*pc(2,15,1)) - zeta_b_zeta*pc_ol(&
&5,15)

pc(19,15,1) = PA(2)*pc(6,15,1) + inv_2zeta*(2*pc(3,15,1) + 2*pc(6,10,1))

pc(20,15,1) = PA(3)*pc(7,15,1) + inv_2zeta*(2*pc(4,15,1) + pc(7,6,1))

pc(1,16,1) = PB(1)*pc(1,7,1) + zeta_a_zeta*pc_ol(1,7)

pc(2,16,1) = PA(1)*pc(1,16,1) + inv_2zeta*(pc(1,7,1)) - zeta_b_zeta*pc_ol(1,16&
&)

pc(3,16,1) = PA(2)*pc(1,16,1)

pc(4,16,1) = PA(3)*pc(1,16,1) + inv_2zeta*(2*pc(1,9,1))

pc(5,16,1) = PA(1)*pc(2,16,1) + inv_2zeta*(pc(1,16,1) + pc(2,7,1)) - zeta_b_ze&
&ta*pc_ol(2,16)

pc(6,16,1) = PA(2)*pc(3,16,1) + inv_2zeta*(pc(1,16,1))

pc(7,16,1) = PA(3)*pc(4,16,1) + inv_2zeta*(pc(1,16,1) + 2*pc(4,9,1))

pc(8,16,1) = PA(2)*pc(2,16,1)

pc(9,16,1) = PA(3)*pc(2,16,1) + inv_2zeta*(2*pc(2,9,1))

pc(10,16,1) = PA(3)*pc(3,16,1) + inv_2zeta*(2*pc(3,9,1))

pc(11,16,1) = PA(3)*pc(8,16,1) + inv_2zeta*(2*pc(8,9,1))

pc(12,16,1) = PA(2)*pc(5,16,1)

pc(13,16,1) = PA(3)*pc(5,16,1) + inv_2zeta*(2*pc(5,9,1))

pc(14,16,1) = PA(1)*pc(6,16,1) + inv_2zeta*(pc(6,7,1)) - zeta_b_zeta*pc_ol(6,1&
&6)

pc(15,16,1) = PA(3)*pc(6,16,1) + inv_2zeta*(2*pc(6,9,1))

pc(16,16,1) = PA(1)*pc(7,16,1) + inv_2zeta*(pc(7,7,1)) - zeta_b_zeta*pc_ol(7,1&
&6)

pc(17,16,1) = PA(2)*pc(7,16,1)

pc(18,16,1) = PA(1)*pc(5,16,1) + inv_2zeta*(2*pc(2,16,1) + pc(5,7,1)) - zeta_b&
&_zeta*pc_ol(5,16)

pc(19,16,1) = PA(2)*pc(6,16,1) + inv_2zeta*(2*pc(3,16,1))

pc(20,16,1) = PA(3)*pc(7,16,1) + inv_2zeta*(2*pc(4,16,1) + 2*pc(7,9,1))

pc(1,17,1) = PB(2)*pc(1,7,1)

pc(2,17,1) = PA(1)*pc(1,17,1) - zeta_b_zeta*pc_ol(1,17)

pc(3,17,1) = PA(2)*pc(1,17,1) + inv_2zeta*(pc(1,7,1))

pc(4,17,1) = PA(3)*pc(1,17,1) + inv_2zeta*(2*pc(1,10,1))

pc(5,17,1) = PA(1)*pc(2,17,1) + inv_2zeta*(pc(1,17,1)) - zeta_b_zeta*pc_ol(2,1&
&7)

pc(6,17,1) = PA(2)*pc(3,17,1) + inv_2zeta*(pc(1,17,1) + pc(3,7,1))

pc(7,17,1) = PA(3)*pc(4,17,1) + inv_2zeta*(pc(1,17,1) + 2*pc(4,10,1))

pc(8,17,1) = PA(2)*pc(2,17,1) + inv_2zeta*(pc(2,7,1))

pc(9,17,1) = PA(3)*pc(2,17,1) + inv_2zeta*(2*pc(2,10,1))

pc(10,17,1) = PA(3)*pc(3,17,1) + inv_2zeta*(2*pc(3,10,1))

pc(11,17,1) = PA(3)*pc(8,17,1) + inv_2zeta*(2*pc(8,10,1))

pc(12,17,1) = PA(2)*pc(5,17,1) + inv_2zeta*(pc(5,7,1))

pc(13,17,1) = PA(3)*pc(5,17,1) + inv_2zeta*(2*pc(5,10,1))

pc(14,17,1) = PA(1)*pc(6,17,1) - zeta_b_zeta*pc_ol(6,17)

pc(15,17,1) = PA(3)*pc(6,17,1) + inv_2zeta*(2*pc(6,10,1))

pc(16,17,1) = PA(1)*pc(7,17,1) - zeta_b_zeta*pc_ol(7,17)

pc(17,17,1) = PA(2)*pc(7,17,1) + inv_2zeta*(pc(7,7,1))

pc(18,17,1) = PA(1)*pc(5,17,1) + inv_2zeta*(2*pc(2,17,1)) - zeta_b_zeta*pc_ol(&
&5,17)

pc(19,17,1) = PA(2)*pc(6,17,1) + inv_2zeta*(2*pc(3,17,1) + pc(6,7,1))

pc(20,17,1) = PA(3)*pc(7,17,1) + inv_2zeta*(2*pc(4,17,1) + 2*pc(7,10,1))

pc(1,18,1) = PB(1)*pc(1,5,1) + inv_2zeta*(2*pc(1,2,1)) + zeta_a_zeta*pc_ol(1,5&
&)

pc(2,18,1) = PA(1)*pc(1,18,1) + inv_2zeta*(3*pc(1,5,1)) - zeta_b_zeta*pc_ol(1,&
&18)

pc(3,18,1) = PA(2)*pc(1,18,1)

pc(4,18,1) = PA(3)*pc(1,18,1)

pc(5,18,1) = PA(1)*pc(2,18,1) + inv_2zeta*(pc(1,18,1) + 3*pc(2,5,1)) - zeta_b_&
&zeta*pc_ol(2,18)

pc(6,18,1) = PA(2)*pc(3,18,1) + inv_2zeta*(pc(1,18,1))

pc(7,18,1) = PA(3)*pc(4,18,1) + inv_2zeta*(pc(1,18,1))

pc(8,18,1) = PA(2)*pc(2,18,1)

pc(9,18,1) = PA(3)*pc(2,18,1)

pc(10,18,1) = PA(3)*pc(3,18,1)

pc(11,18,1) = PA(3)*pc(8,18,1)

pc(12,18,1) = PA(2)*pc(5,18,1)

pc(13,18,1) = PA(3)*pc(5,18,1)

pc(14,18,1) = PA(1)*pc(6,18,1) + inv_2zeta*(3*pc(6,5,1)) - zeta_b_zeta*pc_ol(6&
&,18)

pc(15,18,1) = PA(3)*pc(6,18,1)

pc(16,18,1) = PA(1)*pc(7,18,1) + inv_2zeta*(3*pc(7,5,1)) - zeta_b_zeta*pc_ol(7&
&,18)

pc(17,18,1) = PA(2)*pc(7,18,1)

pc(18,18,1) = PA(1)*pc(5,18,1) + inv_2zeta*(2*pc(2,18,1) + 3*pc(5,5,1)) - zeta&
&_b_zeta*pc_ol(5,18)

pc(19,18,1) = PA(2)*pc(6,18,1) + inv_2zeta*(2*pc(3,18,1))

pc(20,18,1) = PA(3)*pc(7,18,1) + inv_2zeta*(2*pc(4,18,1))

pc(1,19,1) = PB(2)*pc(1,6,1) + inv_2zeta*(2*pc(1,3,1))

pc(2,19,1) = PA(1)*pc(1,19,1) - zeta_b_zeta*pc_ol(1,19)

pc(3,19,1) = PA(2)*pc(1,19,1) + inv_2zeta*(3*pc(1,6,1))

pc(4,19,1) = PA(3)*pc(1,19,1)

pc(5,19,1) = PA(1)*pc(2,19,1) + inv_2zeta*(pc(1,19,1)) - zeta_b_zeta*pc_ol(2,1&
&9)

pc(6,19,1) = PA(2)*pc(3,19,1) + inv_2zeta*(pc(1,19,1) + 3*pc(3,6,1))

pc(7,19,1) = PA(3)*pc(4,19,1) + inv_2zeta*(pc(1,19,1))

pc(8,19,1) = PA(2)*pc(2,19,1) + inv_2zeta*(3*pc(2,6,1))

pc(9,19,1) = PA(3)*pc(2,19,1)

pc(10,19,1) = PA(3)*pc(3,19,1)

pc(11,19,1) = PA(3)*pc(8,19,1)

pc(12,19,1) = PA(2)*pc(5,19,1) + inv_2zeta*(3*pc(5,6,1))

pc(13,19,1) = PA(3)*pc(5,19,1)

pc(14,19,1) = PA(1)*pc(6,19,1) - zeta_b_zeta*pc_ol(6,19)

pc(15,19,1) = PA(3)*pc(6,19,1)

pc(16,19,1) = PA(1)*pc(7,19,1) - zeta_b_zeta*pc_ol(7,19)

pc(17,19,1) = PA(2)*pc(7,19,1) + inv_2zeta*(3*pc(7,6,1))

pc(18,19,1) = PA(1)*pc(5,19,1) + inv_2zeta*(2*pc(2,19,1)) - zeta_b_zeta*pc_ol(&
&5,19)

pc(19,19,1) = PA(2)*pc(6,19,1) + inv_2zeta*(2*pc(3,19,1) + 3*pc(6,6,1))

pc(20,19,1) = PA(3)*pc(7,19,1) + inv_2zeta*(2*pc(4,19,1))

pc(1,20,1) = PB(3)*pc(1,7,1) + inv_2zeta*(2*pc(1,4,1))

pc(2,20,1) = PA(1)*pc(1,20,1) - zeta_b_zeta*pc_ol(1,20)

pc(3,20,1) = PA(2)*pc(1,20,1)

pc(4,20,1) = PA(3)*pc(1,20,1) + inv_2zeta*(3*pc(1,7,1))

pc(5,20,1) = PA(1)*pc(2,20,1) + inv_2zeta*(pc(1,20,1)) - zeta_b_zeta*pc_ol(2,2&
&0)

pc(6,20,1) = PA(2)*pc(3,20,1) + inv_2zeta*(pc(1,20,1))

pc(7,20,1) = PA(3)*pc(4,20,1) + inv_2zeta*(pc(1,20,1) + 3*pc(4,7,1))

pc(8,20,1) = PA(2)*pc(2,20,1)

pc(9,20,1) = PA(3)*pc(2,20,1) + inv_2zeta*(3*pc(2,7,1))

pc(10,20,1) = PA(3)*pc(3,20,1) + inv_2zeta*(3*pc(3,7,1))

pc(11,20,1) = PA(3)*pc(8,20,1) + inv_2zeta*(3*pc(8,7,1))

pc(12,20,1) = PA(2)*pc(5,20,1)

pc(13,20,1) = PA(3)*pc(5,20,1) + inv_2zeta*(3*pc(5,7,1))

pc(14,20,1) = PA(1)*pc(6,20,1) - zeta_b_zeta*pc_ol(6,20)

pc(15,20,1) = PA(3)*pc(6,20,1) + inv_2zeta*(3*pc(6,7,1))

pc(16,20,1) = PA(1)*pc(7,20,1) - zeta_b_zeta*pc_ol(7,20)

pc(17,20,1) = PA(2)*pc(7,20,1)

pc(18,20,1) = PA(1)*pc(5,20,1) + inv_2zeta*(2*pc(2,20,1)) - zeta_b_zeta*pc_ol(&
&5,20)

pc(19,20,1) = PA(2)*pc(6,20,1) + inv_2zeta*(2*pc(3,20,1))

pc(20,20,1) = PA(3)*pc(7,20,1) + inv_2zeta*(2*pc(4,20,1) + 3*pc(7,7,1))

pc(1,1,2) = preFactorMM(2)

pc(2,1,2) = PA(1)*preFactorMM(2)

pc(3,1,2) = PA(2)*preFactorMM(2) - zeta_b_zeta*preFactorOL

pc(4,1,2) = PA(3)*preFactorMM(2)

pc(5,1,2) = PA(1)*pc(2,1,2) + inv_2zeta*(preFactorMM(2))

pc(6,1,2) = PA(2)*pc(3,1,2) + inv_2zeta*(preFactorMM(2)) - zeta_b_zeta*pc_ol(3&
&,1)

pc(7,1,2) = PA(3)*pc(4,1,2) + inv_2zeta*(preFactorMM(2))

pc(8,1,2) = PA(2)*pc(2,1,2) - zeta_b_zeta*pc_ol(2,1)

pc(9,1,2) = PA(3)*pc(2,1,2)

pc(10,1,2) = PA(3)*pc(3,1,2)

pc(11,1,2) = PA(3)*pc(8,1,2)

pc(12,1,2) = PA(2)*pc(5,1,2) - zeta_b_zeta*pc_ol(5,1)

pc(13,1,2) = PA(3)*pc(5,1,2)

pc(14,1,2) = PA(1)*pc(6,1,2)

pc(15,1,2) = PA(3)*pc(6,1,2)

pc(16,1,2) = PA(1)*pc(7,1,2)

pc(17,1,2) = PA(2)*pc(7,1,2) - zeta_b_zeta*pc_ol(7,1)

pc(18,1,2) = PA(1)*pc(5,1,2) + inv_2zeta*(2*pc(2,1,2))

pc(19,1,2) = PA(2)*pc(6,1,2) + inv_2zeta*(2*pc(3,1,2)) - zeta_b_zeta*pc_ol(6,1&
&)

pc(20,1,2) = PA(3)*pc(7,1,2) + inv_2zeta*(2*pc(4,1,2))

pc(1,2,2) = PB(1)*preFactorMM(2)

pc(2,2,2) = PA(1)*pc(1,2,2) + inv_2zeta*(preFactorMM(2))

pc(3,2,2) = PA(2)*pc(1,2,2) - zeta_b_zeta*pc_ol(1,2)

pc(4,2,2) = PA(3)*pc(1,2,2)

pc(5,2,2) = PA(1)*pc(2,2,2) + inv_2zeta*(pc(1,2,2) + pc(2,1,2))

pc(6,2,2) = PA(2)*pc(3,2,2) + inv_2zeta*(pc(1,2,2)) - zeta_b_zeta*pc_ol(3,2)

pc(7,2,2) = PA(3)*pc(4,2,2) + inv_2zeta*(pc(1,2,2))

pc(8,2,2) = PA(2)*pc(2,2,2) - zeta_b_zeta*pc_ol(2,2)

pc(9,2,2) = PA(3)*pc(2,2,2)

pc(10,2,2) = PA(3)*pc(3,2,2)

pc(11,2,2) = PA(3)*pc(8,2,2)

pc(12,2,2) = PA(2)*pc(5,2,2) - zeta_b_zeta*pc_ol(5,2)

pc(13,2,2) = PA(3)*pc(5,2,2)

pc(14,2,2) = PA(1)*pc(6,2,2) + inv_2zeta*(pc(6,1,2))

pc(15,2,2) = PA(3)*pc(6,2,2)

pc(16,2,2) = PA(1)*pc(7,2,2) + inv_2zeta*(pc(7,1,2))

pc(17,2,2) = PA(2)*pc(7,2,2) - zeta_b_zeta*pc_ol(7,2)

pc(18,2,2) = PA(1)*pc(5,2,2) + inv_2zeta*(2*pc(2,2,2) + pc(5,1,2))

pc(19,2,2) = PA(2)*pc(6,2,2) + inv_2zeta*(2*pc(3,2,2)) - zeta_b_zeta*pc_ol(6,2&
&)

pc(20,2,2) = PA(3)*pc(7,2,2) + inv_2zeta*(2*pc(4,2,2))

pc(1,3,2) = PB(2)*preFactorMM(2) + zeta_a_zeta*preFactorOL

pc(2,3,2) = PA(1)*pc(1,3,2)

pc(3,3,2) = PA(2)*pc(1,3,2) + inv_2zeta*(preFactorMM(2)) - zeta_b_zeta*pc_ol(1&
&,3)

pc(4,3,2) = PA(3)*pc(1,3,2)

pc(5,3,2) = PA(1)*pc(2,3,2) + inv_2zeta*(pc(1,3,2))

pc(6,3,2) = PA(2)*pc(3,3,2) + inv_2zeta*(pc(1,3,2) + pc(3,1,2)) - zeta_b_zeta*&
&pc_ol(3,3)

pc(7,3,2) = PA(3)*pc(4,3,2) + inv_2zeta*(pc(1,3,2))

pc(8,3,2) = PA(2)*pc(2,3,2) + inv_2zeta*(pc(2,1,2)) - zeta_b_zeta*pc_ol(2,3)

pc(9,3,2) = PA(3)*pc(2,3,2)

pc(10,3,2) = PA(3)*pc(3,3,2)

pc(11,3,2) = PA(3)*pc(8,3,2)

pc(12,3,2) = PA(2)*pc(5,3,2) + inv_2zeta*(pc(5,1,2)) - zeta_b_zeta*pc_ol(5,3)

pc(13,3,2) = PA(3)*pc(5,3,2)

pc(14,3,2) = PA(1)*pc(6,3,2)

pc(15,3,2) = PA(3)*pc(6,3,2)

pc(16,3,2) = PA(1)*pc(7,3,2)

pc(17,3,2) = PA(2)*pc(7,3,2) + inv_2zeta*(pc(7,1,2)) - zeta_b_zeta*pc_ol(7,3)

pc(18,3,2) = PA(1)*pc(5,3,2) + inv_2zeta*(2*pc(2,3,2))

pc(19,3,2) = PA(2)*pc(6,3,2) + inv_2zeta*(2*pc(3,3,2) + pc(6,1,2)) - zeta_b_ze&
&ta*pc_ol(6,3)

pc(20,3,2) = PA(3)*pc(7,3,2) + inv_2zeta*(2*pc(4,3,2))

pc(1,4,2) = PB(3)*preFactorMM(2)

pc(2,4,2) = PA(1)*pc(1,4,2)

pc(3,4,2) = PA(2)*pc(1,4,2) - zeta_b_zeta*pc_ol(1,4)

pc(4,4,2) = PA(3)*pc(1,4,2) + inv_2zeta*(preFactorMM(2))

pc(5,4,2) = PA(1)*pc(2,4,2) + inv_2zeta*(pc(1,4,2))

pc(6,4,2) = PA(2)*pc(3,4,2) + inv_2zeta*(pc(1,4,2)) - zeta_b_zeta*pc_ol(3,4)

pc(7,4,2) = PA(3)*pc(4,4,2) + inv_2zeta*(pc(1,4,2) + pc(4,1,2))

pc(8,4,2) = PA(2)*pc(2,4,2) - zeta_b_zeta*pc_ol(2,4)

pc(9,4,2) = PA(3)*pc(2,4,2) + inv_2zeta*(pc(2,1,2))

pc(10,4,2) = PA(3)*pc(3,4,2) + inv_2zeta*(pc(3,1,2))

pc(11,4,2) = PA(3)*pc(8,4,2) + inv_2zeta*(pc(8,1,2))

pc(12,4,2) = PA(2)*pc(5,4,2) - zeta_b_zeta*pc_ol(5,4)

pc(13,4,2) = PA(3)*pc(5,4,2) + inv_2zeta*(pc(5,1,2))

pc(14,4,2) = PA(1)*pc(6,4,2)

pc(15,4,2) = PA(3)*pc(6,4,2) + inv_2zeta*(pc(6,1,2))

pc(16,4,2) = PA(1)*pc(7,4,2)

pc(17,4,2) = PA(2)*pc(7,4,2) - zeta_b_zeta*pc_ol(7,4)

pc(18,4,2) = PA(1)*pc(5,4,2) + inv_2zeta*(2*pc(2,4,2))

pc(19,4,2) = PA(2)*pc(6,4,2) + inv_2zeta*(2*pc(3,4,2)) - zeta_b_zeta*pc_ol(6,4&
&)

pc(20,4,2) = PA(3)*pc(7,4,2) + inv_2zeta*(2*pc(4,4,2) + pc(7,1,2))

pc(1,5,2) = PB(1)*pc(1,2,2) + inv_2zeta*(preFactorMM(2))

pc(2,5,2) = PA(1)*pc(1,5,2) + inv_2zeta*(2*pc(1,2,2))

pc(3,5,2) = PA(2)*pc(1,5,2) - zeta_b_zeta*pc_ol(1,5)

pc(4,5,2) = PA(3)*pc(1,5,2)

pc(5,5,2) = PA(1)*pc(2,5,2) + inv_2zeta*(pc(1,5,2) + 2*pc(2,2,2))

pc(6,5,2) = PA(2)*pc(3,5,2) + inv_2zeta*(pc(1,5,2)) - zeta_b_zeta*pc_ol(3,5)

pc(7,5,2) = PA(3)*pc(4,5,2) + inv_2zeta*(pc(1,5,2))

pc(8,5,2) = PA(2)*pc(2,5,2) - zeta_b_zeta*pc_ol(2,5)

pc(9,5,2) = PA(3)*pc(2,5,2)

pc(10,5,2) = PA(3)*pc(3,5,2)

pc(11,5,2) = PA(3)*pc(8,5,2)

pc(12,5,2) = PA(2)*pc(5,5,2) - zeta_b_zeta*pc_ol(5,5)

pc(13,5,2) = PA(3)*pc(5,5,2)

pc(14,5,2) = PA(1)*pc(6,5,2) + inv_2zeta*(2*pc(6,2,2))

pc(15,5,2) = PA(3)*pc(6,5,2)

pc(16,5,2) = PA(1)*pc(7,5,2) + inv_2zeta*(2*pc(7,2,2))

pc(17,5,2) = PA(2)*pc(7,5,2) - zeta_b_zeta*pc_ol(7,5)

pc(18,5,2) = PA(1)*pc(5,5,2) + inv_2zeta*(2*pc(2,5,2) + 2*pc(5,2,2))

pc(19,5,2) = PA(2)*pc(6,5,2) + inv_2zeta*(2*pc(3,5,2)) - zeta_b_zeta*pc_ol(6,5&
&)

pc(20,5,2) = PA(3)*pc(7,5,2) + inv_2zeta*(2*pc(4,5,2))

pc(1,6,2) = PB(2)*pc(1,3,2) + inv_2zeta*(preFactorMM(2)) + zeta_a_zeta*pc_ol(1&
&,3)

pc(2,6,2) = PA(1)*pc(1,6,2)

pc(3,6,2) = PA(2)*pc(1,6,2) + inv_2zeta*(2*pc(1,3,2)) - zeta_b_zeta*pc_ol(1,6)

pc(4,6,2) = PA(3)*pc(1,6,2)

pc(5,6,2) = PA(1)*pc(2,6,2) + inv_2zeta*(pc(1,6,2))

pc(6,6,2) = PA(2)*pc(3,6,2) + inv_2zeta*(pc(1,6,2) + 2*pc(3,3,2)) - zeta_b_zet&
&a*pc_ol(3,6)

pc(7,6,2) = PA(3)*pc(4,6,2) + inv_2zeta*(pc(1,6,2))

pc(8,6,2) = PA(2)*pc(2,6,2) + inv_2zeta*(2*pc(2,3,2)) - zeta_b_zeta*pc_ol(2,6)

pc(9,6,2) = PA(3)*pc(2,6,2)

pc(10,6,2) = PA(3)*pc(3,6,2)

pc(11,6,2) = PA(3)*pc(8,6,2)

pc(12,6,2) = PA(2)*pc(5,6,2) + inv_2zeta*(2*pc(5,3,2)) - zeta_b_zeta*pc_ol(5,6&
&)

pc(13,6,2) = PA(3)*pc(5,6,2)

pc(14,6,2) = PA(1)*pc(6,6,2)

pc(15,6,2) = PA(3)*pc(6,6,2)

pc(16,6,2) = PA(1)*pc(7,6,2)

pc(17,6,2) = PA(2)*pc(7,6,2) + inv_2zeta*(2*pc(7,3,2)) - zeta_b_zeta*pc_ol(7,6&
&)

pc(18,6,2) = PA(1)*pc(5,6,2) + inv_2zeta*(2*pc(2,6,2))

pc(19,6,2) = PA(2)*pc(6,6,2) + inv_2zeta*(2*pc(3,6,2) + 2*pc(6,3,2)) - zeta_b_&
&zeta*pc_ol(6,6)

pc(20,6,2) = PA(3)*pc(7,6,2) + inv_2zeta*(2*pc(4,6,2))

pc(1,7,2) = PB(3)*pc(1,4,2) + inv_2zeta*(preFactorMM(2))

pc(2,7,2) = PA(1)*pc(1,7,2)

pc(3,7,2) = PA(2)*pc(1,7,2) - zeta_b_zeta*pc_ol(1,7)

pc(4,7,2) = PA(3)*pc(1,7,2) + inv_2zeta*(2*pc(1,4,2))

pc(5,7,2) = PA(1)*pc(2,7,2) + inv_2zeta*(pc(1,7,2))

pc(6,7,2) = PA(2)*pc(3,7,2) + inv_2zeta*(pc(1,7,2)) - zeta_b_zeta*pc_ol(3,7)

pc(7,7,2) = PA(3)*pc(4,7,2) + inv_2zeta*(pc(1,7,2) + 2*pc(4,4,2))

pc(8,7,2) = PA(2)*pc(2,7,2) - zeta_b_zeta*pc_ol(2,7)

pc(9,7,2) = PA(3)*pc(2,7,2) + inv_2zeta*(2*pc(2,4,2))

pc(10,7,2) = PA(3)*pc(3,7,2) + inv_2zeta*(2*pc(3,4,2))

pc(11,7,2) = PA(3)*pc(8,7,2) + inv_2zeta*(2*pc(8,4,2))

pc(12,7,2) = PA(2)*pc(5,7,2) - zeta_b_zeta*pc_ol(5,7)

pc(13,7,2) = PA(3)*pc(5,7,2) + inv_2zeta*(2*pc(5,4,2))

pc(14,7,2) = PA(1)*pc(6,7,2)

pc(15,7,2) = PA(3)*pc(6,7,2) + inv_2zeta*(2*pc(6,4,2))

pc(16,7,2) = PA(1)*pc(7,7,2)

pc(17,7,2) = PA(2)*pc(7,7,2) - zeta_b_zeta*pc_ol(7,7)

pc(18,7,2) = PA(1)*pc(5,7,2) + inv_2zeta*(2*pc(2,7,2))

pc(19,7,2) = PA(2)*pc(6,7,2) + inv_2zeta*(2*pc(3,7,2)) - zeta_b_zeta*pc_ol(6,7&
&)

pc(20,7,2) = PA(3)*pc(7,7,2) + inv_2zeta*(2*pc(4,7,2) + 2*pc(7,4,2))

pc(1,8,2) = PB(2)*pc(1,2,2) + zeta_a_zeta*pc_ol(1,2)

pc(2,8,2) = PA(1)*pc(1,8,2) + inv_2zeta*(pc(1,3,2))

pc(3,8,2) = PA(2)*pc(1,8,2) + inv_2zeta*(pc(1,2,2)) - zeta_b_zeta*pc_ol(1,8)

pc(4,8,2) = PA(3)*pc(1,8,2)

pc(5,8,2) = PA(1)*pc(2,8,2) + inv_2zeta*(pc(1,8,2) + pc(2,3,2))

pc(6,8,2) = PA(2)*pc(3,8,2) + inv_2zeta*(pc(1,8,2) + pc(3,2,2)) - zeta_b_zeta*&
&pc_ol(3,8)

pc(7,8,2) = PA(3)*pc(4,8,2) + inv_2zeta*(pc(1,8,2))

pc(8,8,2) = PA(2)*pc(2,8,2) + inv_2zeta*(pc(2,2,2)) - zeta_b_zeta*pc_ol(2,8)

pc(9,8,2) = PA(3)*pc(2,8,2)

pc(10,8,2) = PA(3)*pc(3,8,2)

pc(11,8,2) = PA(3)*pc(8,8,2)

pc(12,8,2) = PA(2)*pc(5,8,2) + inv_2zeta*(pc(5,2,2)) - zeta_b_zeta*pc_ol(5,8)

pc(13,8,2) = PA(3)*pc(5,8,2)

pc(14,8,2) = PA(1)*pc(6,8,2) + inv_2zeta*(pc(6,3,2))

pc(15,8,2) = PA(3)*pc(6,8,2)

pc(16,8,2) = PA(1)*pc(7,8,2) + inv_2zeta*(pc(7,3,2))

pc(17,8,2) = PA(2)*pc(7,8,2) + inv_2zeta*(pc(7,2,2)) - zeta_b_zeta*pc_ol(7,8)

pc(18,8,2) = PA(1)*pc(5,8,2) + inv_2zeta*(2*pc(2,8,2) + pc(5,3,2))

pc(19,8,2) = PA(2)*pc(6,8,2) + inv_2zeta*(2*pc(3,8,2) + pc(6,2,2)) - zeta_b_ze&
&ta*pc_ol(6,8)

pc(20,8,2) = PA(3)*pc(7,8,2) + inv_2zeta*(2*pc(4,8,2))

pc(1,9,2) = PB(3)*pc(1,2,2)

pc(2,9,2) = PA(1)*pc(1,9,2) + inv_2zeta*(pc(1,4,2))

pc(3,9,2) = PA(2)*pc(1,9,2) - zeta_b_zeta*pc_ol(1,9)

pc(4,9,2) = PA(3)*pc(1,9,2) + inv_2zeta*(pc(1,2,2))

pc(5,9,2) = PA(1)*pc(2,9,2) + inv_2zeta*(pc(1,9,2) + pc(2,4,2))

pc(6,9,2) = PA(2)*pc(3,9,2) + inv_2zeta*(pc(1,9,2)) - zeta_b_zeta*pc_ol(3,9)

pc(7,9,2) = PA(3)*pc(4,9,2) + inv_2zeta*(pc(1,9,2) + pc(4,2,2))

pc(8,9,2) = PA(2)*pc(2,9,2) - zeta_b_zeta*pc_ol(2,9)

pc(9,9,2) = PA(3)*pc(2,9,2) + inv_2zeta*(pc(2,2,2))

pc(10,9,2) = PA(3)*pc(3,9,2) + inv_2zeta*(pc(3,2,2))

pc(11,9,2) = PA(3)*pc(8,9,2) + inv_2zeta*(pc(8,2,2))

pc(12,9,2) = PA(2)*pc(5,9,2) - zeta_b_zeta*pc_ol(5,9)

pc(13,9,2) = PA(3)*pc(5,9,2) + inv_2zeta*(pc(5,2,2))

pc(14,9,2) = PA(1)*pc(6,9,2) + inv_2zeta*(pc(6,4,2))

pc(15,9,2) = PA(3)*pc(6,9,2) + inv_2zeta*(pc(6,2,2))

pc(16,9,2) = PA(1)*pc(7,9,2) + inv_2zeta*(pc(7,4,2))

pc(17,9,2) = PA(2)*pc(7,9,2) - zeta_b_zeta*pc_ol(7,9)

pc(18,9,2) = PA(1)*pc(5,9,2) + inv_2zeta*(2*pc(2,9,2) + pc(5,4,2))

pc(19,9,2) = PA(2)*pc(6,9,2) + inv_2zeta*(2*pc(3,9,2)) - zeta_b_zeta*pc_ol(6,9&
&)

pc(20,9,2) = PA(3)*pc(7,9,2) + inv_2zeta*(2*pc(4,9,2) + pc(7,2,2))

pc(1,10,2) = PB(3)*pc(1,3,2)

pc(2,10,2) = PA(1)*pc(1,10,2)

pc(3,10,2) = PA(2)*pc(1,10,2) + inv_2zeta*(pc(1,4,2)) - zeta_b_zeta*pc_ol(1,10&
&)

pc(4,10,2) = PA(3)*pc(1,10,2) + inv_2zeta*(pc(1,3,2))

pc(5,10,2) = PA(1)*pc(2,10,2) + inv_2zeta*(pc(1,10,2))

pc(6,10,2) = PA(2)*pc(3,10,2) + inv_2zeta*(pc(1,10,2) + pc(3,4,2)) - zeta_b_ze&
&ta*pc_ol(3,10)

pc(7,10,2) = PA(3)*pc(4,10,2) + inv_2zeta*(pc(1,10,2) + pc(4,3,2))

pc(8,10,2) = PA(2)*pc(2,10,2) + inv_2zeta*(pc(2,4,2)) - zeta_b_zeta*pc_ol(2,10&
&)

pc(9,10,2) = PA(3)*pc(2,10,2) + inv_2zeta*(pc(2,3,2))

pc(10,10,2) = PA(3)*pc(3,10,2) + inv_2zeta*(pc(3,3,2))

pc(11,10,2) = PA(3)*pc(8,10,2) + inv_2zeta*(pc(8,3,2))

pc(12,10,2) = PA(2)*pc(5,10,2) + inv_2zeta*(pc(5,4,2)) - zeta_b_zeta*pc_ol(5,1&
&0)

pc(13,10,2) = PA(3)*pc(5,10,2) + inv_2zeta*(pc(5,3,2))

pc(14,10,2) = PA(1)*pc(6,10,2)

pc(15,10,2) = PA(3)*pc(6,10,2) + inv_2zeta*(pc(6,3,2))

pc(16,10,2) = PA(1)*pc(7,10,2)

pc(17,10,2) = PA(2)*pc(7,10,2) + inv_2zeta*(pc(7,4,2)) - zeta_b_zeta*pc_ol(7,1&
&0)

pc(18,10,2) = PA(1)*pc(5,10,2) + inv_2zeta*(2*pc(2,10,2))

pc(19,10,2) = PA(2)*pc(6,10,2) + inv_2zeta*(2*pc(3,10,2) + pc(6,4,2)) - zeta_b&
&_zeta*pc_ol(6,10)

pc(20,10,2) = PA(3)*pc(7,10,2) + inv_2zeta*(2*pc(4,10,2) + pc(7,3,2))

pc(1,11,2) = PB(3)*pc(1,8,2)

pc(2,11,2) = PA(1)*pc(1,11,2) + inv_2zeta*(pc(1,10,2))

pc(3,11,2) = PA(2)*pc(1,11,2) + inv_2zeta*(pc(1,9,2)) - zeta_b_zeta*pc_ol(1,11&
&)

pc(4,11,2) = PA(3)*pc(1,11,2) + inv_2zeta*(pc(1,8,2))

pc(5,11,2) = PA(1)*pc(2,11,2) + inv_2zeta*(pc(1,11,2) + pc(2,10,2))

pc(6,11,2) = PA(2)*pc(3,11,2) + inv_2zeta*(pc(1,11,2) + pc(3,9,2)) - zeta_b_ze&
&ta*pc_ol(3,11)

pc(7,11,2) = PA(3)*pc(4,11,2) + inv_2zeta*(pc(1,11,2) + pc(4,8,2))

pc(8,11,2) = PA(2)*pc(2,11,2) + inv_2zeta*(pc(2,9,2)) - zeta_b_zeta*pc_ol(2,11&
&)

pc(9,11,2) = PA(3)*pc(2,11,2) + inv_2zeta*(pc(2,8,2))

pc(10,11,2) = PA(3)*pc(3,11,2) + inv_2zeta*(pc(3,8,2))

pc(11,11,2) = PA(3)*pc(8,11,2) + inv_2zeta*(pc(8,8,2))

pc(12,11,2) = PA(2)*pc(5,11,2) + inv_2zeta*(pc(5,9,2)) - zeta_b_zeta*pc_ol(5,1&
&1)

pc(13,11,2) = PA(3)*pc(5,11,2) + inv_2zeta*(pc(5,8,2))

pc(14,11,2) = PA(1)*pc(6,11,2) + inv_2zeta*(pc(6,10,2))

pc(15,11,2) = PA(3)*pc(6,11,2) + inv_2zeta*(pc(6,8,2))

pc(16,11,2) = PA(1)*pc(7,11,2) + inv_2zeta*(pc(7,10,2))

pc(17,11,2) = PA(2)*pc(7,11,2) + inv_2zeta*(pc(7,9,2)) - zeta_b_zeta*pc_ol(7,1&
&1)

pc(18,11,2) = PA(1)*pc(5,11,2) + inv_2zeta*(2*pc(2,11,2) + pc(5,10,2))

pc(19,11,2) = PA(2)*pc(6,11,2) + inv_2zeta*(2*pc(3,11,2) + pc(6,9,2)) - zeta_b&
&_zeta*pc_ol(6,11)

pc(20,11,2) = PA(3)*pc(7,11,2) + inv_2zeta*(2*pc(4,11,2) + pc(7,8,2))

pc(1,12,2) = PB(2)*pc(1,5,2) + zeta_a_zeta*pc_ol(1,5)

pc(2,12,2) = PA(1)*pc(1,12,2) + inv_2zeta*(2*pc(1,8,2))

pc(3,12,2) = PA(2)*pc(1,12,2) + inv_2zeta*(pc(1,5,2)) - zeta_b_zeta*pc_ol(1,12&
&)

pc(4,12,2) = PA(3)*pc(1,12,2)

pc(5,12,2) = PA(1)*pc(2,12,2) + inv_2zeta*(pc(1,12,2) + 2*pc(2,8,2))

pc(6,12,2) = PA(2)*pc(3,12,2) + inv_2zeta*(pc(1,12,2) + pc(3,5,2)) - zeta_b_ze&
&ta*pc_ol(3,12)

pc(7,12,2) = PA(3)*pc(4,12,2) + inv_2zeta*(pc(1,12,2))

pc(8,12,2) = PA(2)*pc(2,12,2) + inv_2zeta*(pc(2,5,2)) - zeta_b_zeta*pc_ol(2,12&
&)

pc(9,12,2) = PA(3)*pc(2,12,2)

pc(10,12,2) = PA(3)*pc(3,12,2)

pc(11,12,2) = PA(3)*pc(8,12,2)

pc(12,12,2) = PA(2)*pc(5,12,2) + inv_2zeta*(pc(5,5,2)) - zeta_b_zeta*pc_ol(5,1&
&2)

pc(13,12,2) = PA(3)*pc(5,12,2)

pc(14,12,2) = PA(1)*pc(6,12,2) + inv_2zeta*(2*pc(6,8,2))

pc(15,12,2) = PA(3)*pc(6,12,2)

pc(16,12,2) = PA(1)*pc(7,12,2) + inv_2zeta*(2*pc(7,8,2))

pc(17,12,2) = PA(2)*pc(7,12,2) + inv_2zeta*(pc(7,5,2)) - zeta_b_zeta*pc_ol(7,1&
&2)

pc(18,12,2) = PA(1)*pc(5,12,2) + inv_2zeta*(2*pc(2,12,2) + 2*pc(5,8,2))

pc(19,12,2) = PA(2)*pc(6,12,2) + inv_2zeta*(2*pc(3,12,2) + pc(6,5,2)) - zeta_b&
&_zeta*pc_ol(6,12)

pc(20,12,2) = PA(3)*pc(7,12,2) + inv_2zeta*(2*pc(4,12,2))

pc(1,13,2) = PB(3)*pc(1,5,2)

pc(2,13,2) = PA(1)*pc(1,13,2) + inv_2zeta*(2*pc(1,9,2))

pc(3,13,2) = PA(2)*pc(1,13,2) - zeta_b_zeta*pc_ol(1,13)

pc(4,13,2) = PA(3)*pc(1,13,2) + inv_2zeta*(pc(1,5,2))

pc(5,13,2) = PA(1)*pc(2,13,2) + inv_2zeta*(pc(1,13,2) + 2*pc(2,9,2))

pc(6,13,2) = PA(2)*pc(3,13,2) + inv_2zeta*(pc(1,13,2)) - zeta_b_zeta*pc_ol(3,1&
&3)

pc(7,13,2) = PA(3)*pc(4,13,2) + inv_2zeta*(pc(1,13,2) + pc(4,5,2))

pc(8,13,2) = PA(2)*pc(2,13,2) - zeta_b_zeta*pc_ol(2,13)

pc(9,13,2) = PA(3)*pc(2,13,2) + inv_2zeta*(pc(2,5,2))

pc(10,13,2) = PA(3)*pc(3,13,2) + inv_2zeta*(pc(3,5,2))

pc(11,13,2) = PA(3)*pc(8,13,2) + inv_2zeta*(pc(8,5,2))

pc(12,13,2) = PA(2)*pc(5,13,2) - zeta_b_zeta*pc_ol(5,13)

pc(13,13,2) = PA(3)*pc(5,13,2) + inv_2zeta*(pc(5,5,2))

pc(14,13,2) = PA(1)*pc(6,13,2) + inv_2zeta*(2*pc(6,9,2))

pc(15,13,2) = PA(3)*pc(6,13,2) + inv_2zeta*(pc(6,5,2))

pc(16,13,2) = PA(1)*pc(7,13,2) + inv_2zeta*(2*pc(7,9,2))

pc(17,13,2) = PA(2)*pc(7,13,2) - zeta_b_zeta*pc_ol(7,13)

pc(18,13,2) = PA(1)*pc(5,13,2) + inv_2zeta*(2*pc(2,13,2) + 2*pc(5,9,2))

pc(19,13,2) = PA(2)*pc(6,13,2) + inv_2zeta*(2*pc(3,13,2)) - zeta_b_zeta*pc_ol(&
&6,13)

pc(20,13,2) = PA(3)*pc(7,13,2) + inv_2zeta*(2*pc(4,13,2) + pc(7,5,2))

pc(1,14,2) = PB(1)*pc(1,6,2)

pc(2,14,2) = PA(1)*pc(1,14,2) + inv_2zeta*(pc(1,6,2))

pc(3,14,2) = PA(2)*pc(1,14,2) + inv_2zeta*(2*pc(1,8,2)) - zeta_b_zeta*pc_ol(1,&
&14)

pc(4,14,2) = PA(3)*pc(1,14,2)

pc(5,14,2) = PA(1)*pc(2,14,2) + inv_2zeta*(pc(1,14,2) + pc(2,6,2))

pc(6,14,2) = PA(2)*pc(3,14,2) + inv_2zeta*(pc(1,14,2) + 2*pc(3,8,2)) - zeta_b_&
&zeta*pc_ol(3,14)

pc(7,14,2) = PA(3)*pc(4,14,2) + inv_2zeta*(pc(1,14,2))

pc(8,14,2) = PA(2)*pc(2,14,2) + inv_2zeta*(2*pc(2,8,2)) - zeta_b_zeta*pc_ol(2,&
&14)

pc(9,14,2) = PA(3)*pc(2,14,2)

pc(10,14,2) = PA(3)*pc(3,14,2)

pc(11,14,2) = PA(3)*pc(8,14,2)

pc(12,14,2) = PA(2)*pc(5,14,2) + inv_2zeta*(2*pc(5,8,2)) - zeta_b_zeta*pc_ol(5&
&,14)

pc(13,14,2) = PA(3)*pc(5,14,2)

pc(14,14,2) = PA(1)*pc(6,14,2) + inv_2zeta*(pc(6,6,2))

pc(15,14,2) = PA(3)*pc(6,14,2)

pc(16,14,2) = PA(1)*pc(7,14,2) + inv_2zeta*(pc(7,6,2))

pc(17,14,2) = PA(2)*pc(7,14,2) + inv_2zeta*(2*pc(7,8,2)) - zeta_b_zeta*pc_ol(7&
&,14)

pc(18,14,2) = PA(1)*pc(5,14,2) + inv_2zeta*(2*pc(2,14,2) + pc(5,6,2))

pc(19,14,2) = PA(2)*pc(6,14,2) + inv_2zeta*(2*pc(3,14,2) + 2*pc(6,8,2)) - zeta&
&_b_zeta*pc_ol(6,14)

pc(20,14,2) = PA(3)*pc(7,14,2) + inv_2zeta*(2*pc(4,14,2))

pc(1,15,2) = PB(3)*pc(1,6,2)

pc(2,15,2) = PA(1)*pc(1,15,2)

pc(3,15,2) = PA(2)*pc(1,15,2) + inv_2zeta*(2*pc(1,10,2)) - zeta_b_zeta*pc_ol(1&
&,15)

pc(4,15,2) = PA(3)*pc(1,15,2) + inv_2zeta*(pc(1,6,2))

pc(5,15,2) = PA(1)*pc(2,15,2) + inv_2zeta*(pc(1,15,2))

pc(6,15,2) = PA(2)*pc(3,15,2) + inv_2zeta*(pc(1,15,2) + 2*pc(3,10,2)) - zeta_b&
&_zeta*pc_ol(3,15)

pc(7,15,2) = PA(3)*pc(4,15,2) + inv_2zeta*(pc(1,15,2) + pc(4,6,2))

pc(8,15,2) = PA(2)*pc(2,15,2) + inv_2zeta*(2*pc(2,10,2)) - zeta_b_zeta*pc_ol(2&
&,15)

pc(9,15,2) = PA(3)*pc(2,15,2) + inv_2zeta*(pc(2,6,2))

pc(10,15,2) = PA(3)*pc(3,15,2) + inv_2zeta*(pc(3,6,2))

pc(11,15,2) = PA(3)*pc(8,15,2) + inv_2zeta*(pc(8,6,2))

pc(12,15,2) = PA(2)*pc(5,15,2) + inv_2zeta*(2*pc(5,10,2)) - zeta_b_zeta*pc_ol(&
&5,15)

pc(13,15,2) = PA(3)*pc(5,15,2) + inv_2zeta*(pc(5,6,2))

pc(14,15,2) = PA(1)*pc(6,15,2)

pc(15,15,2) = PA(3)*pc(6,15,2) + inv_2zeta*(pc(6,6,2))

pc(16,15,2) = PA(1)*pc(7,15,2)

pc(17,15,2) = PA(2)*pc(7,15,2) + inv_2zeta*(2*pc(7,10,2)) - zeta_b_zeta*pc_ol(&
&7,15)

pc(18,15,2) = PA(1)*pc(5,15,2) + inv_2zeta*(2*pc(2,15,2))

pc(19,15,2) = PA(2)*pc(6,15,2) + inv_2zeta*(2*pc(3,15,2) + 2*pc(6,10,2)) - zet&
&a_b_zeta*pc_ol(6,15)

pc(20,15,2) = PA(3)*pc(7,15,2) + inv_2zeta*(2*pc(4,15,2) + pc(7,6,2))

pc(1,16,2) = PB(1)*pc(1,7,2)

pc(2,16,2) = PA(1)*pc(1,16,2) + inv_2zeta*(pc(1,7,2))

pc(3,16,2) = PA(2)*pc(1,16,2) - zeta_b_zeta*pc_ol(1,16)

pc(4,16,2) = PA(3)*pc(1,16,2) + inv_2zeta*(2*pc(1,9,2))

pc(5,16,2) = PA(1)*pc(2,16,2) + inv_2zeta*(pc(1,16,2) + pc(2,7,2))

pc(6,16,2) = PA(2)*pc(3,16,2) + inv_2zeta*(pc(1,16,2)) - zeta_b_zeta*pc_ol(3,1&
&6)

pc(7,16,2) = PA(3)*pc(4,16,2) + inv_2zeta*(pc(1,16,2) + 2*pc(4,9,2))

pc(8,16,2) = PA(2)*pc(2,16,2) - zeta_b_zeta*pc_ol(2,16)

pc(9,16,2) = PA(3)*pc(2,16,2) + inv_2zeta*(2*pc(2,9,2))

pc(10,16,2) = PA(3)*pc(3,16,2) + inv_2zeta*(2*pc(3,9,2))

pc(11,16,2) = PA(3)*pc(8,16,2) + inv_2zeta*(2*pc(8,9,2))

pc(12,16,2) = PA(2)*pc(5,16,2) - zeta_b_zeta*pc_ol(5,16)

pc(13,16,2) = PA(3)*pc(5,16,2) + inv_2zeta*(2*pc(5,9,2))

pc(14,16,2) = PA(1)*pc(6,16,2) + inv_2zeta*(pc(6,7,2))

pc(15,16,2) = PA(3)*pc(6,16,2) + inv_2zeta*(2*pc(6,9,2))

pc(16,16,2) = PA(1)*pc(7,16,2) + inv_2zeta*(pc(7,7,2))

pc(17,16,2) = PA(2)*pc(7,16,2) - zeta_b_zeta*pc_ol(7,16)

pc(18,16,2) = PA(1)*pc(5,16,2) + inv_2zeta*(2*pc(2,16,2) + pc(5,7,2))

pc(19,16,2) = PA(2)*pc(6,16,2) + inv_2zeta*(2*pc(3,16,2)) - zeta_b_zeta*pc_ol(&
&6,16)

pc(20,16,2) = PA(3)*pc(7,16,2) + inv_2zeta*(2*pc(4,16,2) + 2*pc(7,9,2))

pc(1,17,2) = PB(2)*pc(1,7,2) + zeta_a_zeta*pc_ol(1,7)

pc(2,17,2) = PA(1)*pc(1,17,2)

pc(3,17,2) = PA(2)*pc(1,17,2) + inv_2zeta*(pc(1,7,2)) - zeta_b_zeta*pc_ol(1,17&
&)

pc(4,17,2) = PA(3)*pc(1,17,2) + inv_2zeta*(2*pc(1,10,2))

pc(5,17,2) = PA(1)*pc(2,17,2) + inv_2zeta*(pc(1,17,2))

pc(6,17,2) = PA(2)*pc(3,17,2) + inv_2zeta*(pc(1,17,2) + pc(3,7,2)) - zeta_b_ze&
&ta*pc_ol(3,17)

pc(7,17,2) = PA(3)*pc(4,17,2) + inv_2zeta*(pc(1,17,2) + 2*pc(4,10,2))

pc(8,17,2) = PA(2)*pc(2,17,2) + inv_2zeta*(pc(2,7,2)) - zeta_b_zeta*pc_ol(2,17&
&)

pc(9,17,2) = PA(3)*pc(2,17,2) + inv_2zeta*(2*pc(2,10,2))

pc(10,17,2) = PA(3)*pc(3,17,2) + inv_2zeta*(2*pc(3,10,2))

pc(11,17,2) = PA(3)*pc(8,17,2) + inv_2zeta*(2*pc(8,10,2))

pc(12,17,2) = PA(2)*pc(5,17,2) + inv_2zeta*(pc(5,7,2)) - zeta_b_zeta*pc_ol(5,1&
&7)

pc(13,17,2) = PA(3)*pc(5,17,2) + inv_2zeta*(2*pc(5,10,2))

pc(14,17,2) = PA(1)*pc(6,17,2)

pc(15,17,2) = PA(3)*pc(6,17,2) + inv_2zeta*(2*pc(6,10,2))

pc(16,17,2) = PA(1)*pc(7,17,2)

pc(17,17,2) = PA(2)*pc(7,17,2) + inv_2zeta*(pc(7,7,2)) - zeta_b_zeta*pc_ol(7,1&
&7)

pc(18,17,2) = PA(1)*pc(5,17,2) + inv_2zeta*(2*pc(2,17,2))

pc(19,17,2) = PA(2)*pc(6,17,2) + inv_2zeta*(2*pc(3,17,2) + pc(6,7,2)) - zeta_b&
&_zeta*pc_ol(6,17)

pc(20,17,2) = PA(3)*pc(7,17,2) + inv_2zeta*(2*pc(4,17,2) + 2*pc(7,10,2))

pc(1,18,2) = PB(1)*pc(1,5,2) + inv_2zeta*(2*pc(1,2,2))

pc(2,18,2) = PA(1)*pc(1,18,2) + inv_2zeta*(3*pc(1,5,2))

pc(3,18,2) = PA(2)*pc(1,18,2) - zeta_b_zeta*pc_ol(1,18)

pc(4,18,2) = PA(3)*pc(1,18,2)

pc(5,18,2) = PA(1)*pc(2,18,2) + inv_2zeta*(pc(1,18,2) + 3*pc(2,5,2))

pc(6,18,2) = PA(2)*pc(3,18,2) + inv_2zeta*(pc(1,18,2)) - zeta_b_zeta*pc_ol(3,1&
&8)

pc(7,18,2) = PA(3)*pc(4,18,2) + inv_2zeta*(pc(1,18,2))

pc(8,18,2) = PA(2)*pc(2,18,2) - zeta_b_zeta*pc_ol(2,18)

pc(9,18,2) = PA(3)*pc(2,18,2)

pc(10,18,2) = PA(3)*pc(3,18,2)

pc(11,18,2) = PA(3)*pc(8,18,2)

pc(12,18,2) = PA(2)*pc(5,18,2) - zeta_b_zeta*pc_ol(5,18)

pc(13,18,2) = PA(3)*pc(5,18,2)

pc(14,18,2) = PA(1)*pc(6,18,2) + inv_2zeta*(3*pc(6,5,2))

pc(15,18,2) = PA(3)*pc(6,18,2)

pc(16,18,2) = PA(1)*pc(7,18,2) + inv_2zeta*(3*pc(7,5,2))

pc(17,18,2) = PA(2)*pc(7,18,2) - zeta_b_zeta*pc_ol(7,18)

pc(18,18,2) = PA(1)*pc(5,18,2) + inv_2zeta*(2*pc(2,18,2) + 3*pc(5,5,2))

pc(19,18,2) = PA(2)*pc(6,18,2) + inv_2zeta*(2*pc(3,18,2)) - zeta_b_zeta*pc_ol(&
&6,18)

pc(20,18,2) = PA(3)*pc(7,18,2) + inv_2zeta*(2*pc(4,18,2))

pc(1,19,2) = PB(2)*pc(1,6,2) + inv_2zeta*(2*pc(1,3,2)) + zeta_a_zeta*pc_ol(1,6&
&)

pc(2,19,2) = PA(1)*pc(1,19,2)

pc(3,19,2) = PA(2)*pc(1,19,2) + inv_2zeta*(3*pc(1,6,2)) - zeta_b_zeta*pc_ol(1,&
&19)

pc(4,19,2) = PA(3)*pc(1,19,2)

pc(5,19,2) = PA(1)*pc(2,19,2) + inv_2zeta*(pc(1,19,2))

pc(6,19,2) = PA(2)*pc(3,19,2) + inv_2zeta*(pc(1,19,2) + 3*pc(3,6,2)) - zeta_b_&
&zeta*pc_ol(3,19)

pc(7,19,2) = PA(3)*pc(4,19,2) + inv_2zeta*(pc(1,19,2))

pc(8,19,2) = PA(2)*pc(2,19,2) + inv_2zeta*(3*pc(2,6,2)) - zeta_b_zeta*pc_ol(2,&
&19)

pc(9,19,2) = PA(3)*pc(2,19,2)

pc(10,19,2) = PA(3)*pc(3,19,2)

pc(11,19,2) = PA(3)*pc(8,19,2)

pc(12,19,2) = PA(2)*pc(5,19,2) + inv_2zeta*(3*pc(5,6,2)) - zeta_b_zeta*pc_ol(5&
&,19)

pc(13,19,2) = PA(3)*pc(5,19,2)

pc(14,19,2) = PA(1)*pc(6,19,2)

pc(15,19,2) = PA(3)*pc(6,19,2)

pc(16,19,2) = PA(1)*pc(7,19,2)

pc(17,19,2) = PA(2)*pc(7,19,2) + inv_2zeta*(3*pc(7,6,2)) - zeta_b_zeta*pc_ol(7&
&,19)

pc(18,19,2) = PA(1)*pc(5,19,2) + inv_2zeta*(2*pc(2,19,2))

pc(19,19,2) = PA(2)*pc(6,19,2) + inv_2zeta*(2*pc(3,19,2) + 3*pc(6,6,2)) - zeta&
&_b_zeta*pc_ol(6,19)

pc(20,19,2) = PA(3)*pc(7,19,2) + inv_2zeta*(2*pc(4,19,2))

pc(1,20,2) = PB(3)*pc(1,7,2) + inv_2zeta*(2*pc(1,4,2))

pc(2,20,2) = PA(1)*pc(1,20,2)

pc(3,20,2) = PA(2)*pc(1,20,2) - zeta_b_zeta*pc_ol(1,20)

pc(4,20,2) = PA(3)*pc(1,20,2) + inv_2zeta*(3*pc(1,7,2))

pc(5,20,2) = PA(1)*pc(2,20,2) + inv_2zeta*(pc(1,20,2))

pc(6,20,2) = PA(2)*pc(3,20,2) + inv_2zeta*(pc(1,20,2)) - zeta_b_zeta*pc_ol(3,2&
&0)

pc(7,20,2) = PA(3)*pc(4,20,2) + inv_2zeta*(pc(1,20,2) + 3*pc(4,7,2))

pc(8,20,2) = PA(2)*pc(2,20,2) - zeta_b_zeta*pc_ol(2,20)

pc(9,20,2) = PA(3)*pc(2,20,2) + inv_2zeta*(3*pc(2,7,2))

pc(10,20,2) = PA(3)*pc(3,20,2) + inv_2zeta*(3*pc(3,7,2))

pc(11,20,2) = PA(3)*pc(8,20,2) + inv_2zeta*(3*pc(8,7,2))

pc(12,20,2) = PA(2)*pc(5,20,2) - zeta_b_zeta*pc_ol(5,20)

pc(13,20,2) = PA(3)*pc(5,20,2) + inv_2zeta*(3*pc(5,7,2))

pc(14,20,2) = PA(1)*pc(6,20,2)

pc(15,20,2) = PA(3)*pc(6,20,2) + inv_2zeta*(3*pc(6,7,2))

pc(16,20,2) = PA(1)*pc(7,20,2)

pc(17,20,2) = PA(2)*pc(7,20,2) - zeta_b_zeta*pc_ol(7,20)

pc(18,20,2) = PA(1)*pc(5,20,2) + inv_2zeta*(2*pc(2,20,2))

pc(19,20,2) = PA(2)*pc(6,20,2) + inv_2zeta*(2*pc(3,20,2)) - zeta_b_zeta*pc_ol(&
&6,20)

pc(20,20,2) = PA(3)*pc(7,20,2) + inv_2zeta*(2*pc(4,20,2) + 3*pc(7,7,2))

pc(1,1,3) = preFactorMM(3)

pc(2,1,3) = PA(1)*preFactorMM(3)

pc(3,1,3) = PA(2)*preFactorMM(3)

pc(4,1,3) = PA(3)*preFactorMM(3) - zeta_b_zeta*preFactorOL

pc(5,1,3) = PA(1)*pc(2,1,3) + inv_2zeta*(preFactorMM(3))

pc(6,1,3) = PA(2)*pc(3,1,3) + inv_2zeta*(preFactorMM(3))

pc(7,1,3) = PA(3)*pc(4,1,3) + inv_2zeta*(preFactorMM(3)) - zeta_b_zeta*pc_ol(4&
&,1)

pc(8,1,3) = PA(2)*pc(2,1,3)

pc(9,1,3) = PA(3)*pc(2,1,3) - zeta_b_zeta*pc_ol(2,1)

pc(10,1,3) = PA(3)*pc(3,1,3) - zeta_b_zeta*pc_ol(3,1)

pc(11,1,3) = PA(3)*pc(8,1,3) - zeta_b_zeta*pc_ol(8,1)

pc(12,1,3) = PA(2)*pc(5,1,3)

pc(13,1,3) = PA(3)*pc(5,1,3) - zeta_b_zeta*pc_ol(5,1)

pc(14,1,3) = PA(1)*pc(6,1,3)

pc(15,1,3) = PA(3)*pc(6,1,3) - zeta_b_zeta*pc_ol(6,1)

pc(16,1,3) = PA(1)*pc(7,1,3)

pc(17,1,3) = PA(2)*pc(7,1,3)

pc(18,1,3) = PA(1)*pc(5,1,3) + inv_2zeta*(2*pc(2,1,3))

pc(19,1,3) = PA(2)*pc(6,1,3) + inv_2zeta*(2*pc(3,1,3))

pc(20,1,3) = PA(3)*pc(7,1,3) + inv_2zeta*(2*pc(4,1,3)) - zeta_b_zeta*pc_ol(7,1&
&)

pc(1,2,3) = PB(1)*preFactorMM(3)

pc(2,2,3) = PA(1)*pc(1,2,3) + inv_2zeta*(preFactorMM(3))

pc(3,2,3) = PA(2)*pc(1,2,3)

pc(4,2,3) = PA(3)*pc(1,2,3) - zeta_b_zeta*pc_ol(1,2)

pc(5,2,3) = PA(1)*pc(2,2,3) + inv_2zeta*(pc(1,2,3) + pc(2,1,3))

pc(6,2,3) = PA(2)*pc(3,2,3) + inv_2zeta*(pc(1,2,3))

pc(7,2,3) = PA(3)*pc(4,2,3) + inv_2zeta*(pc(1,2,3)) - zeta_b_zeta*pc_ol(4,2)

pc(8,2,3) = PA(2)*pc(2,2,3)

pc(9,2,3) = PA(3)*pc(2,2,3) - zeta_b_zeta*pc_ol(2,2)

pc(10,2,3) = PA(3)*pc(3,2,3) - zeta_b_zeta*pc_ol(3,2)

pc(11,2,3) = PA(3)*pc(8,2,3) - zeta_b_zeta*pc_ol(8,2)

pc(12,2,3) = PA(2)*pc(5,2,3)

pc(13,2,3) = PA(3)*pc(5,2,3) - zeta_b_zeta*pc_ol(5,2)

pc(14,2,3) = PA(1)*pc(6,2,3) + inv_2zeta*(pc(6,1,3))

pc(15,2,3) = PA(3)*pc(6,2,3) - zeta_b_zeta*pc_ol(6,2)

pc(16,2,3) = PA(1)*pc(7,2,3) + inv_2zeta*(pc(7,1,3))

pc(17,2,3) = PA(2)*pc(7,2,3)

pc(18,2,3) = PA(1)*pc(5,2,3) + inv_2zeta*(2*pc(2,2,3) + pc(5,1,3))

pc(19,2,3) = PA(2)*pc(6,2,3) + inv_2zeta*(2*pc(3,2,3))

pc(20,2,3) = PA(3)*pc(7,2,3) + inv_2zeta*(2*pc(4,2,3)) - zeta_b_zeta*pc_ol(7,2&
&)

pc(1,3,3) = PB(2)*preFactorMM(3)

pc(2,3,3) = PA(1)*pc(1,3,3)

pc(3,3,3) = PA(2)*pc(1,3,3) + inv_2zeta*(preFactorMM(3))

pc(4,3,3) = PA(3)*pc(1,3,3) - zeta_b_zeta*pc_ol(1,3)

pc(5,3,3) = PA(1)*pc(2,3,3) + inv_2zeta*(pc(1,3,3))

pc(6,3,3) = PA(2)*pc(3,3,3) + inv_2zeta*(pc(1,3,3) + pc(3,1,3))

pc(7,3,3) = PA(3)*pc(4,3,3) + inv_2zeta*(pc(1,3,3)) - zeta_b_zeta*pc_ol(4,3)

pc(8,3,3) = PA(2)*pc(2,3,3) + inv_2zeta*(pc(2,1,3))

pc(9,3,3) = PA(3)*pc(2,3,3) - zeta_b_zeta*pc_ol(2,3)

pc(10,3,3) = PA(3)*pc(3,3,3) - zeta_b_zeta*pc_ol(3,3)

pc(11,3,3) = PA(3)*pc(8,3,3) - zeta_b_zeta*pc_ol(8,3)

pc(12,3,3) = PA(2)*pc(5,3,3) + inv_2zeta*(pc(5,1,3))

pc(13,3,3) = PA(3)*pc(5,3,3) - zeta_b_zeta*pc_ol(5,3)

pc(14,3,3) = PA(1)*pc(6,3,3)

pc(15,3,3) = PA(3)*pc(6,3,3) - zeta_b_zeta*pc_ol(6,3)

pc(16,3,3) = PA(1)*pc(7,3,3)

pc(17,3,3) = PA(2)*pc(7,3,3) + inv_2zeta*(pc(7,1,3))

pc(18,3,3) = PA(1)*pc(5,3,3) + inv_2zeta*(2*pc(2,3,3))

pc(19,3,3) = PA(2)*pc(6,3,3) + inv_2zeta*(2*pc(3,3,3) + pc(6,1,3))

pc(20,3,3) = PA(3)*pc(7,3,3) + inv_2zeta*(2*pc(4,3,3)) - zeta_b_zeta*pc_ol(7,3&
&)

pc(1,4,3) = PB(3)*preFactorMM(3) + zeta_a_zeta*preFactorOL

pc(2,4,3) = PA(1)*pc(1,4,3)

pc(3,4,3) = PA(2)*pc(1,4,3)

pc(4,4,3) = PA(3)*pc(1,4,3) + inv_2zeta*(preFactorMM(3)) - zeta_b_zeta*pc_ol(1&
&,4)

pc(5,4,3) = PA(1)*pc(2,4,3) + inv_2zeta*(pc(1,4,3))

pc(6,4,3) = PA(2)*pc(3,4,3) + inv_2zeta*(pc(1,4,3))

pc(7,4,3) = PA(3)*pc(4,4,3) + inv_2zeta*(pc(1,4,3) + pc(4,1,3)) - zeta_b_zeta*&
&pc_ol(4,4)

pc(8,4,3) = PA(2)*pc(2,4,3)

pc(9,4,3) = PA(3)*pc(2,4,3) + inv_2zeta*(pc(2,1,3)) - zeta_b_zeta*pc_ol(2,4)

pc(10,4,3) = PA(3)*pc(3,4,3) + inv_2zeta*(pc(3,1,3)) - zeta_b_zeta*pc_ol(3,4)

pc(11,4,3) = PA(3)*pc(8,4,3) + inv_2zeta*(pc(8,1,3)) - zeta_b_zeta*pc_ol(8,4)

pc(12,4,3) = PA(2)*pc(5,4,3)

pc(13,4,3) = PA(3)*pc(5,4,3) + inv_2zeta*(pc(5,1,3)) - zeta_b_zeta*pc_ol(5,4)

pc(14,4,3) = PA(1)*pc(6,4,3)

pc(15,4,3) = PA(3)*pc(6,4,3) + inv_2zeta*(pc(6,1,3)) - zeta_b_zeta*pc_ol(6,4)

pc(16,4,3) = PA(1)*pc(7,4,3)

pc(17,4,3) = PA(2)*pc(7,4,3)

pc(18,4,3) = PA(1)*pc(5,4,3) + inv_2zeta*(2*pc(2,4,3))

pc(19,4,3) = PA(2)*pc(6,4,3) + inv_2zeta*(2*pc(3,4,3))

pc(20,4,3) = PA(3)*pc(7,4,3) + inv_2zeta*(2*pc(4,4,3) + pc(7,1,3)) - zeta_b_ze&
&ta*pc_ol(7,4)

pc(1,5,3) = PB(1)*pc(1,2,3) + inv_2zeta*(preFactorMM(3))

pc(2,5,3) = PA(1)*pc(1,5,3) + inv_2zeta*(2*pc(1,2,3))

pc(3,5,3) = PA(2)*pc(1,5,3)

pc(4,5,3) = PA(3)*pc(1,5,3) - zeta_b_zeta*pc_ol(1,5)

pc(5,5,3) = PA(1)*pc(2,5,3) + inv_2zeta*(pc(1,5,3) + 2*pc(2,2,3))

pc(6,5,3) = PA(2)*pc(3,5,3) + inv_2zeta*(pc(1,5,3))

pc(7,5,3) = PA(3)*pc(4,5,3) + inv_2zeta*(pc(1,5,3)) - zeta_b_zeta*pc_ol(4,5)

pc(8,5,3) = PA(2)*pc(2,5,3)

pc(9,5,3) = PA(3)*pc(2,5,3) - zeta_b_zeta*pc_ol(2,5)

pc(10,5,3) = PA(3)*pc(3,5,3) - zeta_b_zeta*pc_ol(3,5)

pc(11,5,3) = PA(3)*pc(8,5,3) - zeta_b_zeta*pc_ol(8,5)

pc(12,5,3) = PA(2)*pc(5,5,3)

pc(13,5,3) = PA(3)*pc(5,5,3) - zeta_b_zeta*pc_ol(5,5)

pc(14,5,3) = PA(1)*pc(6,5,3) + inv_2zeta*(2*pc(6,2,3))

pc(15,5,3) = PA(3)*pc(6,5,3) - zeta_b_zeta*pc_ol(6,5)

pc(16,5,3) = PA(1)*pc(7,5,3) + inv_2zeta*(2*pc(7,2,3))

pc(17,5,3) = PA(2)*pc(7,5,3)

pc(18,5,3) = PA(1)*pc(5,5,3) + inv_2zeta*(2*pc(2,5,3) + 2*pc(5,2,3))

pc(19,5,3) = PA(2)*pc(6,5,3) + inv_2zeta*(2*pc(3,5,3))

pc(20,5,3) = PA(3)*pc(7,5,3) + inv_2zeta*(2*pc(4,5,3)) - zeta_b_zeta*pc_ol(7,5&
&)

pc(1,6,3) = PB(2)*pc(1,3,3) + inv_2zeta*(preFactorMM(3))

pc(2,6,3) = PA(1)*pc(1,6,3)

pc(3,6,3) = PA(2)*pc(1,6,3) + inv_2zeta*(2*pc(1,3,3))

pc(4,6,3) = PA(3)*pc(1,6,3) - zeta_b_zeta*pc_ol(1,6)

pc(5,6,3) = PA(1)*pc(2,6,3) + inv_2zeta*(pc(1,6,3))

pc(6,6,3) = PA(2)*pc(3,6,3) + inv_2zeta*(pc(1,6,3) + 2*pc(3,3,3))

pc(7,6,3) = PA(3)*pc(4,6,3) + inv_2zeta*(pc(1,6,3)) - zeta_b_zeta*pc_ol(4,6)

pc(8,6,3) = PA(2)*pc(2,6,3) + inv_2zeta*(2*pc(2,3,3))

pc(9,6,3) = PA(3)*pc(2,6,3) - zeta_b_zeta*pc_ol(2,6)

pc(10,6,3) = PA(3)*pc(3,6,3) - zeta_b_zeta*pc_ol(3,6)

pc(11,6,3) = PA(3)*pc(8,6,3) - zeta_b_zeta*pc_ol(8,6)

pc(12,6,3) = PA(2)*pc(5,6,3) + inv_2zeta*(2*pc(5,3,3))

pc(13,6,3) = PA(3)*pc(5,6,3) - zeta_b_zeta*pc_ol(5,6)

pc(14,6,3) = PA(1)*pc(6,6,3)

pc(15,6,3) = PA(3)*pc(6,6,3) - zeta_b_zeta*pc_ol(6,6)

pc(16,6,3) = PA(1)*pc(7,6,3)

pc(17,6,3) = PA(2)*pc(7,6,3) + inv_2zeta*(2*pc(7,3,3))

pc(18,6,3) = PA(1)*pc(5,6,3) + inv_2zeta*(2*pc(2,6,3))

pc(19,6,3) = PA(2)*pc(6,6,3) + inv_2zeta*(2*pc(3,6,3) + 2*pc(6,3,3))

pc(20,6,3) = PA(3)*pc(7,6,3) + inv_2zeta*(2*pc(4,6,3)) - zeta_b_zeta*pc_ol(7,6&
&)

pc(1,7,3) = PB(3)*pc(1,4,3) + inv_2zeta*(preFactorMM(3)) + zeta_a_zeta*pc_ol(1&
&,4)

pc(2,7,3) = PA(1)*pc(1,7,3)

pc(3,7,3) = PA(2)*pc(1,7,3)

pc(4,7,3) = PA(3)*pc(1,7,3) + inv_2zeta*(2*pc(1,4,3)) - zeta_b_zeta*pc_ol(1,7)

pc(5,7,3) = PA(1)*pc(2,7,3) + inv_2zeta*(pc(1,7,3))

pc(6,7,3) = PA(2)*pc(3,7,3) + inv_2zeta*(pc(1,7,3))

pc(7,7,3) = PA(3)*pc(4,7,3) + inv_2zeta*(pc(1,7,3) + 2*pc(4,4,3)) - zeta_b_zet&
&a*pc_ol(4,7)

pc(8,7,3) = PA(2)*pc(2,7,3)

pc(9,7,3) = PA(3)*pc(2,7,3) + inv_2zeta*(2*pc(2,4,3)) - zeta_b_zeta*pc_ol(2,7)

pc(10,7,3) = PA(3)*pc(3,7,3) + inv_2zeta*(2*pc(3,4,3)) - zeta_b_zeta*pc_ol(3,7&
&)

pc(11,7,3) = PA(3)*pc(8,7,3) + inv_2zeta*(2*pc(8,4,3)) - zeta_b_zeta*pc_ol(8,7&
&)

pc(12,7,3) = PA(2)*pc(5,7,3)

pc(13,7,3) = PA(3)*pc(5,7,3) + inv_2zeta*(2*pc(5,4,3)) - zeta_b_zeta*pc_ol(5,7&
&)

pc(14,7,3) = PA(1)*pc(6,7,3)

pc(15,7,3) = PA(3)*pc(6,7,3) + inv_2zeta*(2*pc(6,4,3)) - zeta_b_zeta*pc_ol(6,7&
&)

pc(16,7,3) = PA(1)*pc(7,7,3)

pc(17,7,3) = PA(2)*pc(7,7,3)

pc(18,7,3) = PA(1)*pc(5,7,3) + inv_2zeta*(2*pc(2,7,3))

pc(19,7,3) = PA(2)*pc(6,7,3) + inv_2zeta*(2*pc(3,7,3))

pc(20,7,3) = PA(3)*pc(7,7,3) + inv_2zeta*(2*pc(4,7,3) + 2*pc(7,4,3)) - zeta_b_&
&zeta*pc_ol(7,7)

pc(1,8,3) = PB(2)*pc(1,2,3)

pc(2,8,3) = PA(1)*pc(1,8,3) + inv_2zeta*(pc(1,3,3))

pc(3,8,3) = PA(2)*pc(1,8,3) + inv_2zeta*(pc(1,2,3))

pc(4,8,3) = PA(3)*pc(1,8,3) - zeta_b_zeta*pc_ol(1,8)

pc(5,8,3) = PA(1)*pc(2,8,3) + inv_2zeta*(pc(1,8,3) + pc(2,3,3))

pc(6,8,3) = PA(2)*pc(3,8,3) + inv_2zeta*(pc(1,8,3) + pc(3,2,3))

pc(7,8,3) = PA(3)*pc(4,8,3) + inv_2zeta*(pc(1,8,3)) - zeta_b_zeta*pc_ol(4,8)

pc(8,8,3) = PA(2)*pc(2,8,3) + inv_2zeta*(pc(2,2,3))

pc(9,8,3) = PA(3)*pc(2,8,3) - zeta_b_zeta*pc_ol(2,8)

pc(10,8,3) = PA(3)*pc(3,8,3) - zeta_b_zeta*pc_ol(3,8)

pc(11,8,3) = PA(3)*pc(8,8,3) - zeta_b_zeta*pc_ol(8,8)

pc(12,8,3) = PA(2)*pc(5,8,3) + inv_2zeta*(pc(5,2,3))

pc(13,8,3) = PA(3)*pc(5,8,3) - zeta_b_zeta*pc_ol(5,8)

pc(14,8,3) = PA(1)*pc(6,8,3) + inv_2zeta*(pc(6,3,3))

pc(15,8,3) = PA(3)*pc(6,8,3) - zeta_b_zeta*pc_ol(6,8)

pc(16,8,3) = PA(1)*pc(7,8,3) + inv_2zeta*(pc(7,3,3))

pc(17,8,3) = PA(2)*pc(7,8,3) + inv_2zeta*(pc(7,2,3))

pc(18,8,3) = PA(1)*pc(5,8,3) + inv_2zeta*(2*pc(2,8,3) + pc(5,3,3))

pc(19,8,3) = PA(2)*pc(6,8,3) + inv_2zeta*(2*pc(3,8,3) + pc(6,2,3))

pc(20,8,3) = PA(3)*pc(7,8,3) + inv_2zeta*(2*pc(4,8,3)) - zeta_b_zeta*pc_ol(7,8&
&)

pc(1,9,3) = PB(3)*pc(1,2,3) + zeta_a_zeta*pc_ol(1,2)

pc(2,9,3) = PA(1)*pc(1,9,3) + inv_2zeta*(pc(1,4,3))

pc(3,9,3) = PA(2)*pc(1,9,3)

pc(4,9,3) = PA(3)*pc(1,9,3) + inv_2zeta*(pc(1,2,3)) - zeta_b_zeta*pc_ol(1,9)

pc(5,9,3) = PA(1)*pc(2,9,3) + inv_2zeta*(pc(1,9,3) + pc(2,4,3))

pc(6,9,3) = PA(2)*pc(3,9,3) + inv_2zeta*(pc(1,9,3))

pc(7,9,3) = PA(3)*pc(4,9,3) + inv_2zeta*(pc(1,9,3) + pc(4,2,3)) - zeta_b_zeta*&
&pc_ol(4,9)

pc(8,9,3) = PA(2)*pc(2,9,3)

pc(9,9,3) = PA(3)*pc(2,9,3) + inv_2zeta*(pc(2,2,3)) - zeta_b_zeta*pc_ol(2,9)

pc(10,9,3) = PA(3)*pc(3,9,3) + inv_2zeta*(pc(3,2,3)) - zeta_b_zeta*pc_ol(3,9)

pc(11,9,3) = PA(3)*pc(8,9,3) + inv_2zeta*(pc(8,2,3)) - zeta_b_zeta*pc_ol(8,9)

pc(12,9,3) = PA(2)*pc(5,9,3)

pc(13,9,3) = PA(3)*pc(5,9,3) + inv_2zeta*(pc(5,2,3)) - zeta_b_zeta*pc_ol(5,9)

pc(14,9,3) = PA(1)*pc(6,9,3) + inv_2zeta*(pc(6,4,3))

pc(15,9,3) = PA(3)*pc(6,9,3) + inv_2zeta*(pc(6,2,3)) - zeta_b_zeta*pc_ol(6,9)

pc(16,9,3) = PA(1)*pc(7,9,3) + inv_2zeta*(pc(7,4,3))

pc(17,9,3) = PA(2)*pc(7,9,3)

pc(18,9,3) = PA(1)*pc(5,9,3) + inv_2zeta*(2*pc(2,9,3) + pc(5,4,3))

pc(19,9,3) = PA(2)*pc(6,9,3) + inv_2zeta*(2*pc(3,9,3))

pc(20,9,3) = PA(3)*pc(7,9,3) + inv_2zeta*(2*pc(4,9,3) + pc(7,2,3)) - zeta_b_ze&
&ta*pc_ol(7,9)

pc(1,10,3) = PB(3)*pc(1,3,3) + zeta_a_zeta*pc_ol(1,3)

pc(2,10,3) = PA(1)*pc(1,10,3)

pc(3,10,3) = PA(2)*pc(1,10,3) + inv_2zeta*(pc(1,4,3))

pc(4,10,3) = PA(3)*pc(1,10,3) + inv_2zeta*(pc(1,3,3)) - zeta_b_zeta*pc_ol(1,10&
&)

pc(5,10,3) = PA(1)*pc(2,10,3) + inv_2zeta*(pc(1,10,3))

pc(6,10,3) = PA(2)*pc(3,10,3) + inv_2zeta*(pc(1,10,3) + pc(3,4,3))

pc(7,10,3) = PA(3)*pc(4,10,3) + inv_2zeta*(pc(1,10,3) + pc(4,3,3)) - zeta_b_ze&
&ta*pc_ol(4,10)

pc(8,10,3) = PA(2)*pc(2,10,3) + inv_2zeta*(pc(2,4,3))

pc(9,10,3) = PA(3)*pc(2,10,3) + inv_2zeta*(pc(2,3,3)) - zeta_b_zeta*pc_ol(2,10&
&)

pc(10,10,3) = PA(3)*pc(3,10,3) + inv_2zeta*(pc(3,3,3)) - zeta_b_zeta*pc_ol(3,1&
&0)

pc(11,10,3) = PA(3)*pc(8,10,3) + inv_2zeta*(pc(8,3,3)) - zeta_b_zeta*pc_ol(8,1&
&0)

pc(12,10,3) = PA(2)*pc(5,10,3) + inv_2zeta*(pc(5,4,3))

pc(13,10,3) = PA(3)*pc(5,10,3) + inv_2zeta*(pc(5,3,3)) - zeta_b_zeta*pc_ol(5,1&
&0)

pc(14,10,3) = PA(1)*pc(6,10,3)

pc(15,10,3) = PA(3)*pc(6,10,3) + inv_2zeta*(pc(6,3,3)) - zeta_b_zeta*pc_ol(6,1&
&0)

pc(16,10,3) = PA(1)*pc(7,10,3)

pc(17,10,3) = PA(2)*pc(7,10,3) + inv_2zeta*(pc(7,4,3))

pc(18,10,3) = PA(1)*pc(5,10,3) + inv_2zeta*(2*pc(2,10,3))

pc(19,10,3) = PA(2)*pc(6,10,3) + inv_2zeta*(2*pc(3,10,3) + pc(6,4,3))

pc(20,10,3) = PA(3)*pc(7,10,3) + inv_2zeta*(2*pc(4,10,3) + pc(7,3,3)) - zeta_b&
&_zeta*pc_ol(7,10)

pc(1,11,3) = PB(3)*pc(1,8,3) + zeta_a_zeta*pc_ol(1,8)

pc(2,11,3) = PA(1)*pc(1,11,3) + inv_2zeta*(pc(1,10,3))

pc(3,11,3) = PA(2)*pc(1,11,3) + inv_2zeta*(pc(1,9,3))

pc(4,11,3) = PA(3)*pc(1,11,3) + inv_2zeta*(pc(1,8,3)) - zeta_b_zeta*pc_ol(1,11&
&)

pc(5,11,3) = PA(1)*pc(2,11,3) + inv_2zeta*(pc(1,11,3) + pc(2,10,3))

pc(6,11,3) = PA(2)*pc(3,11,3) + inv_2zeta*(pc(1,11,3) + pc(3,9,3))

pc(7,11,3) = PA(3)*pc(4,11,3) + inv_2zeta*(pc(1,11,3) + pc(4,8,3)) - zeta_b_ze&
&ta*pc_ol(4,11)

pc(8,11,3) = PA(2)*pc(2,11,3) + inv_2zeta*(pc(2,9,3))

pc(9,11,3) = PA(3)*pc(2,11,3) + inv_2zeta*(pc(2,8,3)) - zeta_b_zeta*pc_ol(2,11&
&)

pc(10,11,3) = PA(3)*pc(3,11,3) + inv_2zeta*(pc(3,8,3)) - zeta_b_zeta*pc_ol(3,1&
&1)

pc(11,11,3) = PA(3)*pc(8,11,3) + inv_2zeta*(pc(8,8,3)) - zeta_b_zeta*pc_ol(8,1&
&1)

pc(12,11,3) = PA(2)*pc(5,11,3) + inv_2zeta*(pc(5,9,3))

pc(13,11,3) = PA(3)*pc(5,11,3) + inv_2zeta*(pc(5,8,3)) - zeta_b_zeta*pc_ol(5,1&
&1)

pc(14,11,3) = PA(1)*pc(6,11,3) + inv_2zeta*(pc(6,10,3))

pc(15,11,3) = PA(3)*pc(6,11,3) + inv_2zeta*(pc(6,8,3)) - zeta_b_zeta*pc_ol(6,1&
&1)

pc(16,11,3) = PA(1)*pc(7,11,3) + inv_2zeta*(pc(7,10,3))

pc(17,11,3) = PA(2)*pc(7,11,3) + inv_2zeta*(pc(7,9,3))

pc(18,11,3) = PA(1)*pc(5,11,3) + inv_2zeta*(2*pc(2,11,3) + pc(5,10,3))

pc(19,11,3) = PA(2)*pc(6,11,3) + inv_2zeta*(2*pc(3,11,3) + pc(6,9,3))

pc(20,11,3) = PA(3)*pc(7,11,3) + inv_2zeta*(2*pc(4,11,3) + pc(7,8,3)) - zeta_b&
&_zeta*pc_ol(7,11)

pc(1,12,3) = PB(2)*pc(1,5,3)

pc(2,12,3) = PA(1)*pc(1,12,3) + inv_2zeta*(2*pc(1,8,3))

pc(3,12,3) = PA(2)*pc(1,12,3) + inv_2zeta*(pc(1,5,3))

pc(4,12,3) = PA(3)*pc(1,12,3) - zeta_b_zeta*pc_ol(1,12)

pc(5,12,3) = PA(1)*pc(2,12,3) + inv_2zeta*(pc(1,12,3) + 2*pc(2,8,3))

pc(6,12,3) = PA(2)*pc(3,12,3) + inv_2zeta*(pc(1,12,3) + pc(3,5,3))

pc(7,12,3) = PA(3)*pc(4,12,3) + inv_2zeta*(pc(1,12,3)) - zeta_b_zeta*pc_ol(4,1&
&2)

pc(8,12,3) = PA(2)*pc(2,12,3) + inv_2zeta*(pc(2,5,3))

pc(9,12,3) = PA(3)*pc(2,12,3) - zeta_b_zeta*pc_ol(2,12)

pc(10,12,3) = PA(3)*pc(3,12,3) - zeta_b_zeta*pc_ol(3,12)

pc(11,12,3) = PA(3)*pc(8,12,3) - zeta_b_zeta*pc_ol(8,12)

pc(12,12,3) = PA(2)*pc(5,12,3) + inv_2zeta*(pc(5,5,3))

pc(13,12,3) = PA(3)*pc(5,12,3) - zeta_b_zeta*pc_ol(5,12)

pc(14,12,3) = PA(1)*pc(6,12,3) + inv_2zeta*(2*pc(6,8,3))

pc(15,12,3) = PA(3)*pc(6,12,3) - zeta_b_zeta*pc_ol(6,12)

pc(16,12,3) = PA(1)*pc(7,12,3) + inv_2zeta*(2*pc(7,8,3))

pc(17,12,3) = PA(2)*pc(7,12,3) + inv_2zeta*(pc(7,5,3))

pc(18,12,3) = PA(1)*pc(5,12,3) + inv_2zeta*(2*pc(2,12,3) + 2*pc(5,8,3))

pc(19,12,3) = PA(2)*pc(6,12,3) + inv_2zeta*(2*pc(3,12,3) + pc(6,5,3))

pc(20,12,3) = PA(3)*pc(7,12,3) + inv_2zeta*(2*pc(4,12,3)) - zeta_b_zeta*pc_ol(&
&7,12)

pc(1,13,3) = PB(3)*pc(1,5,3) + zeta_a_zeta*pc_ol(1,5)

pc(2,13,3) = PA(1)*pc(1,13,3) + inv_2zeta*(2*pc(1,9,3))

pc(3,13,3) = PA(2)*pc(1,13,3)

pc(4,13,3) = PA(3)*pc(1,13,3) + inv_2zeta*(pc(1,5,3)) - zeta_b_zeta*pc_ol(1,13&
&)

pc(5,13,3) = PA(1)*pc(2,13,3) + inv_2zeta*(pc(1,13,3) + 2*pc(2,9,3))

pc(6,13,3) = PA(2)*pc(3,13,3) + inv_2zeta*(pc(1,13,3))

pc(7,13,3) = PA(3)*pc(4,13,3) + inv_2zeta*(pc(1,13,3) + pc(4,5,3)) - zeta_b_ze&
&ta*pc_ol(4,13)

pc(8,13,3) = PA(2)*pc(2,13,3)

pc(9,13,3) = PA(3)*pc(2,13,3) + inv_2zeta*(pc(2,5,3)) - zeta_b_zeta*pc_ol(2,13&
&)

pc(10,13,3) = PA(3)*pc(3,13,3) + inv_2zeta*(pc(3,5,3)) - zeta_b_zeta*pc_ol(3,1&
&3)

pc(11,13,3) = PA(3)*pc(8,13,3) + inv_2zeta*(pc(8,5,3)) - zeta_b_zeta*pc_ol(8,1&
&3)

pc(12,13,3) = PA(2)*pc(5,13,3)

pc(13,13,3) = PA(3)*pc(5,13,3) + inv_2zeta*(pc(5,5,3)) - zeta_b_zeta*pc_ol(5,1&
&3)

pc(14,13,3) = PA(1)*pc(6,13,3) + inv_2zeta*(2*pc(6,9,3))

pc(15,13,3) = PA(3)*pc(6,13,3) + inv_2zeta*(pc(6,5,3)) - zeta_b_zeta*pc_ol(6,1&
&3)

pc(16,13,3) = PA(1)*pc(7,13,3) + inv_2zeta*(2*pc(7,9,3))

pc(17,13,3) = PA(2)*pc(7,13,3)

pc(18,13,3) = PA(1)*pc(5,13,3) + inv_2zeta*(2*pc(2,13,3) + 2*pc(5,9,3))

pc(19,13,3) = PA(2)*pc(6,13,3) + inv_2zeta*(2*pc(3,13,3))

pc(20,13,3) = PA(3)*pc(7,13,3) + inv_2zeta*(2*pc(4,13,3) + pc(7,5,3)) - zeta_b&
&_zeta*pc_ol(7,13)

pc(1,14,3) = PB(1)*pc(1,6,3)

pc(2,14,3) = PA(1)*pc(1,14,3) + inv_2zeta*(pc(1,6,3))

pc(3,14,3) = PA(2)*pc(1,14,3) + inv_2zeta*(2*pc(1,8,3))

pc(4,14,3) = PA(3)*pc(1,14,3) - zeta_b_zeta*pc_ol(1,14)

pc(5,14,3) = PA(1)*pc(2,14,3) + inv_2zeta*(pc(1,14,3) + pc(2,6,3))

pc(6,14,3) = PA(2)*pc(3,14,3) + inv_2zeta*(pc(1,14,3) + 2*pc(3,8,3))

pc(7,14,3) = PA(3)*pc(4,14,3) + inv_2zeta*(pc(1,14,3)) - zeta_b_zeta*pc_ol(4,1&
&4)

pc(8,14,3) = PA(2)*pc(2,14,3) + inv_2zeta*(2*pc(2,8,3))

pc(9,14,3) = PA(3)*pc(2,14,3) - zeta_b_zeta*pc_ol(2,14)

pc(10,14,3) = PA(3)*pc(3,14,3) - zeta_b_zeta*pc_ol(3,14)

pc(11,14,3) = PA(3)*pc(8,14,3) - zeta_b_zeta*pc_ol(8,14)

pc(12,14,3) = PA(2)*pc(5,14,3) + inv_2zeta*(2*pc(5,8,3))

pc(13,14,3) = PA(3)*pc(5,14,3) - zeta_b_zeta*pc_ol(5,14)

pc(14,14,3) = PA(1)*pc(6,14,3) + inv_2zeta*(pc(6,6,3))

pc(15,14,3) = PA(3)*pc(6,14,3) - zeta_b_zeta*pc_ol(6,14)

pc(16,14,3) = PA(1)*pc(7,14,3) + inv_2zeta*(pc(7,6,3))

pc(17,14,3) = PA(2)*pc(7,14,3) + inv_2zeta*(2*pc(7,8,3))

pc(18,14,3) = PA(1)*pc(5,14,3) + inv_2zeta*(2*pc(2,14,3) + pc(5,6,3))

pc(19,14,3) = PA(2)*pc(6,14,3) + inv_2zeta*(2*pc(3,14,3) + 2*pc(6,8,3))

pc(20,14,3) = PA(3)*pc(7,14,3) + inv_2zeta*(2*pc(4,14,3)) - zeta_b_zeta*pc_ol(&
&7,14)

pc(1,15,3) = PB(3)*pc(1,6,3) + zeta_a_zeta*pc_ol(1,6)

pc(2,15,3) = PA(1)*pc(1,15,3)

pc(3,15,3) = PA(2)*pc(1,15,3) + inv_2zeta*(2*pc(1,10,3))

pc(4,15,3) = PA(3)*pc(1,15,3) + inv_2zeta*(pc(1,6,3)) - zeta_b_zeta*pc_ol(1,15&
&)

pc(5,15,3) = PA(1)*pc(2,15,3) + inv_2zeta*(pc(1,15,3))

pc(6,15,3) = PA(2)*pc(3,15,3) + inv_2zeta*(pc(1,15,3) + 2*pc(3,10,3))

pc(7,15,3) = PA(3)*pc(4,15,3) + inv_2zeta*(pc(1,15,3) + pc(4,6,3)) - zeta_b_ze&
&ta*pc_ol(4,15)

pc(8,15,3) = PA(2)*pc(2,15,3) + inv_2zeta*(2*pc(2,10,3))

pc(9,15,3) = PA(3)*pc(2,15,3) + inv_2zeta*(pc(2,6,3)) - zeta_b_zeta*pc_ol(2,15&
&)

pc(10,15,3) = PA(3)*pc(3,15,3) + inv_2zeta*(pc(3,6,3)) - zeta_b_zeta*pc_ol(3,1&
&5)

pc(11,15,3) = PA(3)*pc(8,15,3) + inv_2zeta*(pc(8,6,3)) - zeta_b_zeta*pc_ol(8,1&
&5)

pc(12,15,3) = PA(2)*pc(5,15,3) + inv_2zeta*(2*pc(5,10,3))

pc(13,15,3) = PA(3)*pc(5,15,3) + inv_2zeta*(pc(5,6,3)) - zeta_b_zeta*pc_ol(5,1&
&5)

pc(14,15,3) = PA(1)*pc(6,15,3)

pc(15,15,3) = PA(3)*pc(6,15,3) + inv_2zeta*(pc(6,6,3)) - zeta_b_zeta*pc_ol(6,1&
&5)

pc(16,15,3) = PA(1)*pc(7,15,3)

pc(17,15,3) = PA(2)*pc(7,15,3) + inv_2zeta*(2*pc(7,10,3))

pc(18,15,3) = PA(1)*pc(5,15,3) + inv_2zeta*(2*pc(2,15,3))

pc(19,15,3) = PA(2)*pc(6,15,3) + inv_2zeta*(2*pc(3,15,3) + 2*pc(6,10,3))

pc(20,15,3) = PA(3)*pc(7,15,3) + inv_2zeta*(2*pc(4,15,3) + pc(7,6,3)) - zeta_b&
&_zeta*pc_ol(7,15)

pc(1,16,3) = PB(1)*pc(1,7,3)

pc(2,16,3) = PA(1)*pc(1,16,3) + inv_2zeta*(pc(1,7,3))

pc(3,16,3) = PA(2)*pc(1,16,3)

pc(4,16,3) = PA(3)*pc(1,16,3) + inv_2zeta*(2*pc(1,9,3)) - zeta_b_zeta*pc_ol(1,&
&16)

pc(5,16,3) = PA(1)*pc(2,16,3) + inv_2zeta*(pc(1,16,3) + pc(2,7,3))

pc(6,16,3) = PA(2)*pc(3,16,3) + inv_2zeta*(pc(1,16,3))

pc(7,16,3) = PA(3)*pc(4,16,3) + inv_2zeta*(pc(1,16,3) + 2*pc(4,9,3)) - zeta_b_&
&zeta*pc_ol(4,16)

pc(8,16,3) = PA(2)*pc(2,16,3)

pc(9,16,3) = PA(3)*pc(2,16,3) + inv_2zeta*(2*pc(2,9,3)) - zeta_b_zeta*pc_ol(2,&
&16)

pc(10,16,3) = PA(3)*pc(3,16,3) + inv_2zeta*(2*pc(3,9,3)) - zeta_b_zeta*pc_ol(3&
&,16)

pc(11,16,3) = PA(3)*pc(8,16,3) + inv_2zeta*(2*pc(8,9,3)) - zeta_b_zeta*pc_ol(8&
&,16)

pc(12,16,3) = PA(2)*pc(5,16,3)

pc(13,16,3) = PA(3)*pc(5,16,3) + inv_2zeta*(2*pc(5,9,3)) - zeta_b_zeta*pc_ol(5&
&,16)

pc(14,16,3) = PA(1)*pc(6,16,3) + inv_2zeta*(pc(6,7,3))

pc(15,16,3) = PA(3)*pc(6,16,3) + inv_2zeta*(2*pc(6,9,3)) - zeta_b_zeta*pc_ol(6&
&,16)

pc(16,16,3) = PA(1)*pc(7,16,3) + inv_2zeta*(pc(7,7,3))

pc(17,16,3) = PA(2)*pc(7,16,3)

pc(18,16,3) = PA(1)*pc(5,16,3) + inv_2zeta*(2*pc(2,16,3) + pc(5,7,3))

pc(19,16,3) = PA(2)*pc(6,16,3) + inv_2zeta*(2*pc(3,16,3))

pc(20,16,3) = PA(3)*pc(7,16,3) + inv_2zeta*(2*pc(4,16,3) + 2*pc(7,9,3)) - zeta&
&_b_zeta*pc_ol(7,16)

pc(1,17,3) = PB(2)*pc(1,7,3)

pc(2,17,3) = PA(1)*pc(1,17,3)

pc(3,17,3) = PA(2)*pc(1,17,3) + inv_2zeta*(pc(1,7,3))

pc(4,17,3) = PA(3)*pc(1,17,3) + inv_2zeta*(2*pc(1,10,3)) - zeta_b_zeta*pc_ol(1&
&,17)

pc(5,17,3) = PA(1)*pc(2,17,3) + inv_2zeta*(pc(1,17,3))

pc(6,17,3) = PA(2)*pc(3,17,3) + inv_2zeta*(pc(1,17,3) + pc(3,7,3))

pc(7,17,3) = PA(3)*pc(4,17,3) + inv_2zeta*(pc(1,17,3) + 2*pc(4,10,3)) - zeta_b&
&_zeta*pc_ol(4,17)

pc(8,17,3) = PA(2)*pc(2,17,3) + inv_2zeta*(pc(2,7,3))

pc(9,17,3) = PA(3)*pc(2,17,3) + inv_2zeta*(2*pc(2,10,3)) - zeta_b_zeta*pc_ol(2&
&,17)

pc(10,17,3) = PA(3)*pc(3,17,3) + inv_2zeta*(2*pc(3,10,3)) - zeta_b_zeta*pc_ol(&
&3,17)

pc(11,17,3) = PA(3)*pc(8,17,3) + inv_2zeta*(2*pc(8,10,3)) - zeta_b_zeta*pc_ol(&
&8,17)

pc(12,17,3) = PA(2)*pc(5,17,3) + inv_2zeta*(pc(5,7,3))

pc(13,17,3) = PA(3)*pc(5,17,3) + inv_2zeta*(2*pc(5,10,3)) - zeta_b_zeta*pc_ol(&
&5,17)

pc(14,17,3) = PA(1)*pc(6,17,3)

pc(15,17,3) = PA(3)*pc(6,17,3) + inv_2zeta*(2*pc(6,10,3)) - zeta_b_zeta*pc_ol(&
&6,17)

pc(16,17,3) = PA(1)*pc(7,17,3)

pc(17,17,3) = PA(2)*pc(7,17,3) + inv_2zeta*(pc(7,7,3))

pc(18,17,3) = PA(1)*pc(5,17,3) + inv_2zeta*(2*pc(2,17,3))

pc(19,17,3) = PA(2)*pc(6,17,3) + inv_2zeta*(2*pc(3,17,3) + pc(6,7,3))

pc(20,17,3) = PA(3)*pc(7,17,3) + inv_2zeta*(2*pc(4,17,3) + 2*pc(7,10,3)) - zet&
&a_b_zeta*pc_ol(7,17)

pc(1,18,3) = PB(1)*pc(1,5,3) + inv_2zeta*(2*pc(1,2,3))

pc(2,18,3) = PA(1)*pc(1,18,3) + inv_2zeta*(3*pc(1,5,3))

pc(3,18,3) = PA(2)*pc(1,18,3)

pc(4,18,3) = PA(3)*pc(1,18,3) - zeta_b_zeta*pc_ol(1,18)

pc(5,18,3) = PA(1)*pc(2,18,3) + inv_2zeta*(pc(1,18,3) + 3*pc(2,5,3))

pc(6,18,3) = PA(2)*pc(3,18,3) + inv_2zeta*(pc(1,18,3))

pc(7,18,3) = PA(3)*pc(4,18,3) + inv_2zeta*(pc(1,18,3)) - zeta_b_zeta*pc_ol(4,1&
&8)

pc(8,18,3) = PA(2)*pc(2,18,3)

pc(9,18,3) = PA(3)*pc(2,18,3) - zeta_b_zeta*pc_ol(2,18)

pc(10,18,3) = PA(3)*pc(3,18,3) - zeta_b_zeta*pc_ol(3,18)

pc(11,18,3) = PA(3)*pc(8,18,3) - zeta_b_zeta*pc_ol(8,18)

pc(12,18,3) = PA(2)*pc(5,18,3)

pc(13,18,3) = PA(3)*pc(5,18,3) - zeta_b_zeta*pc_ol(5,18)

pc(14,18,3) = PA(1)*pc(6,18,3) + inv_2zeta*(3*pc(6,5,3))

pc(15,18,3) = PA(3)*pc(6,18,3) - zeta_b_zeta*pc_ol(6,18)

pc(16,18,3) = PA(1)*pc(7,18,3) + inv_2zeta*(3*pc(7,5,3))

pc(17,18,3) = PA(2)*pc(7,18,3)

pc(18,18,3) = PA(1)*pc(5,18,3) + inv_2zeta*(2*pc(2,18,3) + 3*pc(5,5,3))

pc(19,18,3) = PA(2)*pc(6,18,3) + inv_2zeta*(2*pc(3,18,3))

pc(20,18,3) = PA(3)*pc(7,18,3) + inv_2zeta*(2*pc(4,18,3)) - zeta_b_zeta*pc_ol(&
&7,18)

pc(1,19,3) = PB(2)*pc(1,6,3) + inv_2zeta*(2*pc(1,3,3))

pc(2,19,3) = PA(1)*pc(1,19,3)

pc(3,19,3) = PA(2)*pc(1,19,3) + inv_2zeta*(3*pc(1,6,3))

pc(4,19,3) = PA(3)*pc(1,19,3) - zeta_b_zeta*pc_ol(1,19)

pc(5,19,3) = PA(1)*pc(2,19,3) + inv_2zeta*(pc(1,19,3))

pc(6,19,3) = PA(2)*pc(3,19,3) + inv_2zeta*(pc(1,19,3) + 3*pc(3,6,3))

pc(7,19,3) = PA(3)*pc(4,19,3) + inv_2zeta*(pc(1,19,3)) - zeta_b_zeta*pc_ol(4,1&
&9)

pc(8,19,3) = PA(2)*pc(2,19,3) + inv_2zeta*(3*pc(2,6,3))

pc(9,19,3) = PA(3)*pc(2,19,3) - zeta_b_zeta*pc_ol(2,19)

pc(10,19,3) = PA(3)*pc(3,19,3) - zeta_b_zeta*pc_ol(3,19)

pc(11,19,3) = PA(3)*pc(8,19,3) - zeta_b_zeta*pc_ol(8,19)

pc(12,19,3) = PA(2)*pc(5,19,3) + inv_2zeta*(3*pc(5,6,3))

pc(13,19,3) = PA(3)*pc(5,19,3) - zeta_b_zeta*pc_ol(5,19)

pc(14,19,3) = PA(1)*pc(6,19,3)

pc(15,19,3) = PA(3)*pc(6,19,3) - zeta_b_zeta*pc_ol(6,19)

pc(16,19,3) = PA(1)*pc(7,19,3)

pc(17,19,3) = PA(2)*pc(7,19,3) + inv_2zeta*(3*pc(7,6,3))

pc(18,19,3) = PA(1)*pc(5,19,3) + inv_2zeta*(2*pc(2,19,3))

pc(19,19,3) = PA(2)*pc(6,19,3) + inv_2zeta*(2*pc(3,19,3) + 3*pc(6,6,3))

pc(20,19,3) = PA(3)*pc(7,19,3) + inv_2zeta*(2*pc(4,19,3)) - zeta_b_zeta*pc_ol(&
&7,19)

pc(1,20,3) = PB(3)*pc(1,7,3) + inv_2zeta*(2*pc(1,4,3)) + zeta_a_zeta*pc_ol(1,7&
&)

pc(2,20,3) = PA(1)*pc(1,20,3)

pc(3,20,3) = PA(2)*pc(1,20,3)

pc(4,20,3) = PA(3)*pc(1,20,3) + inv_2zeta*(3*pc(1,7,3)) - zeta_b_zeta*pc_ol(1,&
&20)

pc(5,20,3) = PA(1)*pc(2,20,3) + inv_2zeta*(pc(1,20,3))

pc(6,20,3) = PA(2)*pc(3,20,3) + inv_2zeta*(pc(1,20,3))

pc(7,20,3) = PA(3)*pc(4,20,3) + inv_2zeta*(pc(1,20,3) + 3*pc(4,7,3)) - zeta_b_&
&zeta*pc_ol(4,20)

pc(8,20,3) = PA(2)*pc(2,20,3)

pc(9,20,3) = PA(3)*pc(2,20,3) + inv_2zeta*(3*pc(2,7,3)) - zeta_b_zeta*pc_ol(2,&
&20)

pc(10,20,3) = PA(3)*pc(3,20,3) + inv_2zeta*(3*pc(3,7,3)) - zeta_b_zeta*pc_ol(3&
&,20)

pc(11,20,3) = PA(3)*pc(8,20,3) + inv_2zeta*(3*pc(8,7,3)) - zeta_b_zeta*pc_ol(8&
&,20)

pc(12,20,3) = PA(2)*pc(5,20,3)

pc(13,20,3) = PA(3)*pc(5,20,3) + inv_2zeta*(3*pc(5,7,3)) - zeta_b_zeta*pc_ol(5&
&,20)

pc(14,20,3) = PA(1)*pc(6,20,3)

pc(15,20,3) = PA(3)*pc(6,20,3) + inv_2zeta*(3*pc(6,7,3)) - zeta_b_zeta*pc_ol(6&
&,20)

pc(16,20,3) = PA(1)*pc(7,20,3)

pc(17,20,3) = PA(2)*pc(7,20,3)

pc(18,20,3) = PA(1)*pc(5,20,3) + inv_2zeta*(2*pc(2,20,3))

pc(19,20,3) = PA(2)*pc(6,20,3) + inv_2zeta*(2*pc(3,20,3))

pc(20,20,3) = PA(3)*pc(7,20,3) + inv_2zeta*(2*pc(4,20,3) + 3*pc(7,7,3)) - zeta&
&_b_zeta*pc_ol(7,20)

sh(1,1,1) = pc(1,1,1)

sh(2,1,1) = pc(2,1,1)

sh(3,1,1) = pc(3,1,1)

sh(4,1,1) = pc(4,1,1)

sh(5,1,1) = pc(8,1,1)

sh(6,1,1) = pc(9,1,1)

sh(7,1,1) = pc(10,1,1)

sh(8,1,1) = pc(5,1,1) - pc(6,1,1)

sh(9,1,1) = 2*pc(7,1,1) - pc(5,1,1) - pc(6,1,1)

sh(10,1,1) = pc(11,1,1)

sh(11,1,1) = pc(13,1,1) - pc(15,1,1)

sh(12,1,1) = pc(18,1,1) - 3*pc(14,1,1)

sh(13,1,1) = 3*pc(12,1,1) - pc(19,1,1)

sh(14,1,1) = 2*pc(20,1,1) - 3*pc(13,1,1) - 3*pc(15,1,1)

sh(15,1,1) = 4*pc(16,1,1) - pc(18,1,1) - pc(14,1,1)

sh(16,1,1) = 4*pc(17,1,1) - pc(12,1,1) - pc(19,1,1)

sh(1,2,1) = pc(1,2,1)

sh(2,2,1) = pc(2,2,1)

sh(3,2,1) = pc(3,2,1)

sh(4,2,1) = pc(4,2,1)

sh(5,2,1) = pc(8,2,1)

sh(6,2,1) = pc(9,2,1)

sh(7,2,1) = pc(10,2,1)

sh(8,2,1) = pc(5,2,1) - pc(6,2,1)

sh(9,2,1) = 2*pc(7,2,1) - pc(5,2,1) - pc(6,2,1)

sh(10,2,1) = pc(11,2,1)

sh(11,2,1) = pc(13,2,1) - pc(15,2,1)

sh(12,2,1) = pc(18,2,1) - 3*pc(14,2,1)

sh(13,2,1) = 3*pc(12,2,1) - pc(19,2,1)

sh(14,2,1) = 2*pc(20,2,1) - 3*pc(13,2,1) - 3*pc(15,2,1)

sh(15,2,1) = 4*pc(16,2,1) - pc(18,2,1) - pc(14,2,1)

sh(16,2,1) = 4*pc(17,2,1) - pc(12,2,1) - pc(19,2,1)

sh(1,3,1) = pc(1,3,1)

sh(2,3,1) = pc(2,3,1)

sh(3,3,1) = pc(3,3,1)

sh(4,3,1) = pc(4,3,1)

sh(5,3,1) = pc(8,3,1)

sh(6,3,1) = pc(9,3,1)

sh(7,3,1) = pc(10,3,1)

sh(8,3,1) = pc(5,3,1) - pc(6,3,1)

sh(9,3,1) = 2*pc(7,3,1) - pc(5,3,1) - pc(6,3,1)

sh(10,3,1) = pc(11,3,1)

sh(11,3,1) = pc(13,3,1) - pc(15,3,1)

sh(12,3,1) = pc(18,3,1) - 3*pc(14,3,1)

sh(13,3,1) = 3*pc(12,3,1) - pc(19,3,1)

sh(14,3,1) = 2*pc(20,3,1) - 3*pc(13,3,1) - 3*pc(15,3,1)

sh(15,3,1) = 4*pc(16,3,1) - pc(18,3,1) - pc(14,3,1)

sh(16,3,1) = 4*pc(17,3,1) - pc(12,3,1) - pc(19,3,1)

sh(1,4,1) = pc(1,4,1)

sh(2,4,1) = pc(2,4,1)

sh(3,4,1) = pc(3,4,1)

sh(4,4,1) = pc(4,4,1)

sh(5,4,1) = pc(8,4,1)

sh(6,4,1) = pc(9,4,1)

sh(7,4,1) = pc(10,4,1)

sh(8,4,1) = pc(5,4,1) - pc(6,4,1)

sh(9,4,1) = 2*pc(7,4,1) - pc(5,4,1) - pc(6,4,1)

sh(10,4,1) = pc(11,4,1)

sh(11,4,1) = pc(13,4,1) - pc(15,4,1)

sh(12,4,1) = pc(18,4,1) - 3*pc(14,4,1)

sh(13,4,1) = 3*pc(12,4,1) - pc(19,4,1)

sh(14,4,1) = 2*pc(20,4,1) - 3*pc(13,4,1) - 3*pc(15,4,1)

sh(15,4,1) = 4*pc(16,4,1) - pc(18,4,1) - pc(14,4,1)

sh(16,4,1) = 4*pc(17,4,1) - pc(12,4,1) - pc(19,4,1)

sh(1,5,1) = pc(1,8,1)

sh(2,5,1) = pc(2,8,1)

sh(3,5,1) = pc(3,8,1)

sh(4,5,1) = pc(4,8,1)

sh(5,5,1) = pc(8,8,1)

sh(6,5,1) = pc(9,8,1)

sh(7,5,1) = pc(10,8,1)

sh(8,5,1) = pc(5,8,1) - pc(6,8,1)

sh(9,5,1) = 2*pc(7,8,1) - pc(5,8,1) - pc(6,8,1)

sh(10,5,1) = pc(11,8,1)

sh(11,5,1) = pc(13,8,1) - pc(15,8,1)

sh(12,5,1) = pc(18,8,1) - 3*pc(14,8,1)

sh(13,5,1) = 3*pc(12,8,1) - pc(19,8,1)

sh(14,5,1) = 2*pc(20,8,1) - 3*pc(13,8,1) - 3*pc(15,8,1)

sh(15,5,1) = 4*pc(16,8,1) - pc(18,8,1) - pc(14,8,1)

sh(16,5,1) = 4*pc(17,8,1) - pc(12,8,1) - pc(19,8,1)

sh(1,6,1) = pc(1,9,1)

sh(2,6,1) = pc(2,9,1)

sh(3,6,1) = pc(3,9,1)

sh(4,6,1) = pc(4,9,1)

sh(5,6,1) = pc(8,9,1)

sh(6,6,1) = pc(9,9,1)

sh(7,6,1) = pc(10,9,1)

sh(8,6,1) = pc(5,9,1) - pc(6,9,1)

sh(9,6,1) = 2*pc(7,9,1) - pc(5,9,1) - pc(6,9,1)

sh(10,6,1) = pc(11,9,1)

sh(11,6,1) = pc(13,9,1) - pc(15,9,1)

sh(12,6,1) = pc(18,9,1) - 3*pc(14,9,1)

sh(13,6,1) = 3*pc(12,9,1) - pc(19,9,1)

sh(14,6,1) = 2*pc(20,9,1) - 3*pc(13,9,1) - 3*pc(15,9,1)

sh(15,6,1) = 4*pc(16,9,1) - pc(18,9,1) - pc(14,9,1)

sh(16,6,1) = 4*pc(17,9,1) - pc(12,9,1) - pc(19,9,1)

sh(1,7,1) = pc(1,10,1)

sh(2,7,1) = pc(2,10,1)

sh(3,7,1) = pc(3,10,1)

sh(4,7,1) = pc(4,10,1)

sh(5,7,1) = pc(8,10,1)

sh(6,7,1) = pc(9,10,1)

sh(7,7,1) = pc(10,10,1)

sh(8,7,1) = pc(5,10,1) - pc(6,10,1)

sh(9,7,1) = 2*pc(7,10,1) - pc(5,10,1) - pc(6,10,1)

sh(10,7,1) = pc(11,10,1)

sh(11,7,1) = pc(13,10,1) - pc(15,10,1)

sh(12,7,1) = pc(18,10,1) - 3*pc(14,10,1)

sh(13,7,1) = 3*pc(12,10,1) - pc(19,10,1)

sh(14,7,1) = 2*pc(20,10,1) - 3*pc(13,10,1) - 3*pc(15,10,1)

sh(15,7,1) = 4*pc(16,10,1) - pc(18,10,1) - pc(14,10,1)

sh(16,7,1) = 4*pc(17,10,1) - pc(12,10,1) - pc(19,10,1)

sh(1,8,1) = pc(1,5,1) - pc(1,6,1)

sh(2,8,1) = pc(2,5,1) - pc(2,6,1)

sh(3,8,1) = pc(3,5,1) - pc(3,6,1)

sh(4,8,1) = pc(4,5,1) - pc(4,6,1)

sh(5,8,1) = pc(8,5,1) - pc(8,6,1)

sh(6,8,1) = pc(9,5,1) - pc(9,6,1)

sh(7,8,1) = pc(10,5,1) - pc(10,6,1)

sh(8,8,1) = pc(5,5,1) - pc(5,6,1) - pc(6,5,1) + pc(6,6,1)

sh(9,8,1) = 2*pc(7,5,1) - 2*pc(7,6,1) - pc(5,5,1) + pc(5,6,1) - pc(6,5,1) + pc&
&(6,6,1)

sh(10,8,1) = pc(11,5,1) - pc(11,6,1)

sh(11,8,1) = pc(13,5,1) - pc(13,6,1) - pc(15,5,1) + pc(15,6,1)

sh(12,8,1) = pc(18,5,1) - pc(18,6,1) - 3*pc(14,5,1) + 3*pc(14,6,1)

sh(13,8,1) = 3*pc(12,5,1) - 3*pc(12,6,1) - pc(19,5,1) + pc(19,6,1)

sh(14,8,1) = 2*pc(20,5,1) - 2*pc(20,6,1) - 3*pc(13,5,1) + 3*pc(13,6,1) - 3*pc(&
&15,5,1) + 3*pc(15,6,1)

sh(15,8,1) = 4*pc(16,5,1) - 4*pc(16,6,1) - pc(18,5,1) + pc(18,6,1) - pc(14,5,1&
&) + pc(14,6,1)

sh(16,8,1) = 4*pc(17,5,1) - 4*pc(17,6,1) - pc(12,5,1) + pc(12,6,1) - pc(19,5,1&
&) + pc(19,6,1)

sh(1,9,1) = 2*pc(1,7,1) - pc(1,5,1) - pc(1,6,1)

sh(2,9,1) = 2*pc(2,7,1) - pc(2,5,1) - pc(2,6,1)

sh(3,9,1) = 2*pc(3,7,1) - pc(3,5,1) - pc(3,6,1)

sh(4,9,1) = 2*pc(4,7,1) - pc(4,5,1) - pc(4,6,1)

sh(5,9,1) = 2*pc(8,7,1) - pc(8,5,1) - pc(8,6,1)

sh(6,9,1) = 2*pc(9,7,1) - pc(9,5,1) - pc(9,6,1)

sh(7,9,1) = 2*pc(10,7,1) - pc(10,5,1) - pc(10,6,1)

sh(8,9,1) = 2*pc(5,7,1) - pc(5,5,1) - pc(5,6,1) - 2*pc(6,7,1) + pc(6,5,1) + pc&
&(6,6,1)

sh(9,9,1) = 4*pc(7,7,1) - 2*pc(7,5,1) - 2*pc(7,6,1) - 2*pc(5,7,1) + pc(5,5,1) &
&+ pc(5,6,1) - 2*pc(6,7,1) + pc(6,5,1) + pc(6,6,1)

sh(10,9,1) = 2*pc(11,7,1) - pc(11,5,1) - pc(11,6,1)

sh(11,9,1) = 2*pc(13,7,1) - pc(13,5,1) - pc(13,6,1) - 2*pc(15,7,1) + pc(15,5,1&
&) + pc(15,6,1)

sh(12,9,1) = 2*pc(18,7,1) - pc(18,5,1) - pc(18,6,1) - 6*pc(14,7,1) + 3*pc(14,5&
&,1) + 3*pc(14,6,1)

sh(13,9,1) = 6*pc(12,7,1) - 3*pc(12,5,1) - 3*pc(12,6,1) - 2*pc(19,7,1) + pc(19&
&,5,1) + pc(19,6,1)

sh(14,9,1) = 4*pc(20,7,1) - 2*pc(20,5,1) - 2*pc(20,6,1) - 6*pc(13,7,1) + 3*pc(&
&13,5,1) + 3*pc(13,6,1) - 6*pc(15,7,1) + 3*pc(15,5,1) + 3*pc(15,6,1)

sh(15,9,1) = 8*pc(16,7,1) - 4*pc(16,5,1) - 4*pc(16,6,1) - 2*pc(18,7,1) + pc(18&
&,5,1) + pc(18,6,1) - 2*pc(14,7,1) + pc(14,5,1) + pc(14,6,1)

sh(16,9,1) = 8*pc(17,7,1) - 4*pc(17,5,1) - 4*pc(17,6,1) - 2*pc(12,7,1) + pc(12&
&,5,1) + pc(12,6,1) - 2*pc(19,7,1) + pc(19,5,1) + pc(19,6,1)

sh(1,10,1) = pc(1,11,1)

sh(2,10,1) = pc(2,11,1)

sh(3,10,1) = pc(3,11,1)

sh(4,10,1) = pc(4,11,1)

sh(5,10,1) = pc(8,11,1)

sh(6,10,1) = pc(9,11,1)

sh(7,10,1) = pc(10,11,1)

sh(8,10,1) = pc(5,11,1) - pc(6,11,1)

sh(9,10,1) = 2*pc(7,11,1) - pc(5,11,1) - pc(6,11,1)

sh(10,10,1) = pc(11,11,1)

sh(11,10,1) = pc(13,11,1) - pc(15,11,1)

sh(12,10,1) = pc(18,11,1) - 3*pc(14,11,1)

sh(13,10,1) = 3*pc(12,11,1) - pc(19,11,1)

sh(14,10,1) = 2*pc(20,11,1) - 3*pc(13,11,1) - 3*pc(15,11,1)

sh(15,10,1) = 4*pc(16,11,1) - pc(18,11,1) - pc(14,11,1)

sh(16,10,1) = 4*pc(17,11,1) - pc(12,11,1) - pc(19,11,1)

sh(1,11,1) = pc(1,13,1) - pc(1,15,1)

sh(2,11,1) = pc(2,13,1) - pc(2,15,1)

sh(3,11,1) = pc(3,13,1) - pc(3,15,1)

sh(4,11,1) = pc(4,13,1) - pc(4,15,1)

sh(5,11,1) = pc(8,13,1) - pc(8,15,1)

sh(6,11,1) = pc(9,13,1) - pc(9,15,1)

sh(7,11,1) = pc(10,13,1) - pc(10,15,1)

sh(8,11,1) = pc(5,13,1) - pc(5,15,1) - pc(6,13,1) + pc(6,15,1)

sh(9,11,1) = 2*pc(7,13,1) - 2*pc(7,15,1) - pc(5,13,1) + pc(5,15,1) - pc(6,13,1&
&) + pc(6,15,1)

sh(10,11,1) = pc(11,13,1) - pc(11,15,1)

sh(11,11,1) = pc(13,13,1) - pc(13,15,1) - pc(15,13,1) + pc(15,15,1)

sh(12,11,1) = pc(18,13,1) - pc(18,15,1) - 3*pc(14,13,1) + 3*pc(14,15,1)

sh(13,11,1) = 3*pc(12,13,1) - 3*pc(12,15,1) - pc(19,13,1) + pc(19,15,1)

sh(14,11,1) = 2*pc(20,13,1) - 2*pc(20,15,1) - 3*pc(13,13,1) + 3*pc(13,15,1) - &
&3*pc(15,13,1) + 3*pc(15,15,1)

sh(15,11,1) = 4*pc(16,13,1) - 4*pc(16,15,1) - pc(18,13,1) + pc(18,15,1) - pc(1&
&4,13,1) + pc(14,15,1)

sh(16,11,1) = 4*pc(17,13,1) - 4*pc(17,15,1) - pc(12,13,1) + pc(12,15,1) - pc(1&
&9,13,1) + pc(19,15,1)

sh(1,12,1) = pc(1,18,1) - 3*pc(1,14,1)

sh(2,12,1) = pc(2,18,1) - 3*pc(2,14,1)

sh(3,12,1) = pc(3,18,1) - 3*pc(3,14,1)

sh(4,12,1) = pc(4,18,1) - 3*pc(4,14,1)

sh(5,12,1) = pc(8,18,1) - 3*pc(8,14,1)

sh(6,12,1) = pc(9,18,1) - 3*pc(9,14,1)

sh(7,12,1) = pc(10,18,1) - 3*pc(10,14,1)

sh(8,12,1) = pc(5,18,1) - 3*pc(5,14,1) - pc(6,18,1) + 3*pc(6,14,1)

sh(9,12,1) = 2*pc(7,18,1) - 6*pc(7,14,1) - pc(5,18,1) + 3*pc(5,14,1) - pc(6,18&
&,1) + 3*pc(6,14,1)

sh(10,12,1) = pc(11,18,1) - 3*pc(11,14,1)

sh(11,12,1) = pc(13,18,1) - 3*pc(13,14,1) - pc(15,18,1) + 3*pc(15,14,1)

sh(12,12,1) = pc(18,18,1) - 3*pc(18,14,1) - 3*pc(14,18,1) + 9*pc(14,14,1)

sh(13,12,1) = 3*pc(12,18,1) - 9*pc(12,14,1) - pc(19,18,1) + 3*pc(19,14,1)

sh(14,12,1) = 2*pc(20,18,1) - 6*pc(20,14,1) - 3*pc(13,18,1) + 9*pc(13,14,1) - &
&3*pc(15,18,1) + 9*pc(15,14,1)

sh(15,12,1) = 4*pc(16,18,1) - 12*pc(16,14,1) - pc(18,18,1) + 3*pc(18,14,1) - p&
&c(14,18,1) + 3*pc(14,14,1)

sh(16,12,1) = 4*pc(17,18,1) - 12*pc(17,14,1) - pc(12,18,1) + 3*pc(12,14,1) - p&
&c(19,18,1) + 3*pc(19,14,1)

sh(1,13,1) = 3*pc(1,12,1) - pc(1,19,1)

sh(2,13,1) = 3*pc(2,12,1) - pc(2,19,1)

sh(3,13,1) = 3*pc(3,12,1) - pc(3,19,1)

sh(4,13,1) = 3*pc(4,12,1) - pc(4,19,1)

sh(5,13,1) = 3*pc(8,12,1) - pc(8,19,1)

sh(6,13,1) = 3*pc(9,12,1) - pc(9,19,1)

sh(7,13,1) = 3*pc(10,12,1) - pc(10,19,1)

sh(8,13,1) = 3*pc(5,12,1) - pc(5,19,1) - 3*pc(6,12,1) + pc(6,19,1)

sh(9,13,1) = 6*pc(7,12,1) - 2*pc(7,19,1) - 3*pc(5,12,1) + pc(5,19,1) - 3*pc(6,&
&12,1) + pc(6,19,1)

sh(10,13,1) = 3*pc(11,12,1) - pc(11,19,1)

sh(11,13,1) = 3*pc(13,12,1) - pc(13,19,1) - 3*pc(15,12,1) + pc(15,19,1)

sh(12,13,1) = 3*pc(18,12,1) - pc(18,19,1) - 9*pc(14,12,1) + 3*pc(14,19,1)

sh(13,13,1) = 9*pc(12,12,1) - 3*pc(12,19,1) - 3*pc(19,12,1) + pc(19,19,1)

sh(14,13,1) = 6*pc(20,12,1) - 2*pc(20,19,1) - 9*pc(13,12,1) + 3*pc(13,19,1) - &
&9*pc(15,12,1) + 3*pc(15,19,1)

sh(15,13,1) = 12*pc(16,12,1) - 4*pc(16,19,1) - 3*pc(18,12,1) + pc(18,19,1) - 3&
&*pc(14,12,1) + pc(14,19,1)

sh(16,13,1) = 12*pc(17,12,1) - 4*pc(17,19,1) - 3*pc(12,12,1) + pc(12,19,1) - 3&
&*pc(19,12,1) + pc(19,19,1)

sh(1,14,1) = 2*pc(1,20,1) - 3*pc(1,13,1) - 3*pc(1,15,1)

sh(2,14,1) = 2*pc(2,20,1) - 3*pc(2,13,1) - 3*pc(2,15,1)

sh(3,14,1) = 2*pc(3,20,1) - 3*pc(3,13,1) - 3*pc(3,15,1)

sh(4,14,1) = 2*pc(4,20,1) - 3*pc(4,13,1) - 3*pc(4,15,1)

sh(5,14,1) = 2*pc(8,20,1) - 3*pc(8,13,1) - 3*pc(8,15,1)

sh(6,14,1) = 2*pc(9,20,1) - 3*pc(9,13,1) - 3*pc(9,15,1)

sh(7,14,1) = 2*pc(10,20,1) - 3*pc(10,13,1) - 3*pc(10,15,1)

sh(8,14,1) = 2*pc(5,20,1) - 3*pc(5,13,1) - 3*pc(5,15,1) - 2*pc(6,20,1) + 3*pc(&
&6,13,1) + 3*pc(6,15,1)

sh(9,14,1) = 4*pc(7,20,1) - 6*pc(7,13,1) - 6*pc(7,15,1) - 2*pc(5,20,1) + 3*pc(&
&5,13,1) + 3*pc(5,15,1) - 2*pc(6,20,1) + 3*pc(6,13,1) + 3*pc(6,15,1)

sh(10,14,1) = 2*pc(11,20,1) - 3*pc(11,13,1) - 3*pc(11,15,1)

sh(11,14,1) = 2*pc(13,20,1) - 3*pc(13,13,1) - 3*pc(13,15,1) - 2*pc(15,20,1) + &
&3*pc(15,13,1) + 3*pc(15,15,1)

sh(12,14,1) = 2*pc(18,20,1) - 3*pc(18,13,1) - 3*pc(18,15,1) - 6*pc(14,20,1) + &
&9*pc(14,13,1) + 9*pc(14,15,1)

sh(13,14,1) = 6*pc(12,20,1) - 9*pc(12,13,1) - 9*pc(12,15,1) - 2*pc(19,20,1) + &
&3*pc(19,13,1) + 3*pc(19,15,1)

sh(14,14,1) = 4*pc(20,20,1) - 6*pc(20,13,1) - 6*pc(20,15,1) - 6*pc(13,20,1) + &
&9*pc(13,13,1) + 9*pc(13,15,1) - 6*pc(15,20,1) + 9*pc(15,13,1) + 9*pc(15,15,1)

sh(15,14,1) = 8*pc(16,20,1) - 12*pc(16,13,1) - 12*pc(16,15,1) - 2*pc(18,20,1) &
&+ 3*pc(18,13,1) + 3*pc(18,15,1) - 2*pc(14,20,1) + 3*pc(14,13,1) + 3*pc(14,15,1&
&)

sh(16,14,1) = 8*pc(17,20,1) - 12*pc(17,13,1) - 12*pc(17,15,1) - 2*pc(12,20,1) &
&+ 3*pc(12,13,1) + 3*pc(12,15,1) - 2*pc(19,20,1) + 3*pc(19,13,1) + 3*pc(19,15,1&
&)

sh(1,15,1) = 4*pc(1,16,1) - pc(1,18,1) - pc(1,14,1)

sh(2,15,1) = 4*pc(2,16,1) - pc(2,18,1) - pc(2,14,1)

sh(3,15,1) = 4*pc(3,16,1) - pc(3,18,1) - pc(3,14,1)

sh(4,15,1) = 4*pc(4,16,1) - pc(4,18,1) - pc(4,14,1)

sh(5,15,1) = 4*pc(8,16,1) - pc(8,18,1) - pc(8,14,1)

sh(6,15,1) = 4*pc(9,16,1) - pc(9,18,1) - pc(9,14,1)

sh(7,15,1) = 4*pc(10,16,1) - pc(10,18,1) - pc(10,14,1)

sh(8,15,1) = 4*pc(5,16,1) - pc(5,18,1) - pc(5,14,1) - 4*pc(6,16,1) + pc(6,18,1&
&) + pc(6,14,1)

sh(9,15,1) = 8*pc(7,16,1) - 2*pc(7,18,1) - 2*pc(7,14,1) - 4*pc(5,16,1) + pc(5,&
&18,1) + pc(5,14,1) - 4*pc(6,16,1) + pc(6,18,1) + pc(6,14,1)

sh(10,15,1) = 4*pc(11,16,1) - pc(11,18,1) - pc(11,14,1)

sh(11,15,1) = 4*pc(13,16,1) - pc(13,18,1) - pc(13,14,1) - 4*pc(15,16,1) + pc(1&
&5,18,1) + pc(15,14,1)

sh(12,15,1) = 4*pc(18,16,1) - pc(18,18,1) - pc(18,14,1) - 12*pc(14,16,1) + 3*p&
&c(14,18,1) + 3*pc(14,14,1)

sh(13,15,1) = 12*pc(12,16,1) - 3*pc(12,18,1) - 3*pc(12,14,1) - 4*pc(19,16,1) +&
& pc(19,18,1) + pc(19,14,1)

sh(14,15,1) = 8*pc(20,16,1) - 2*pc(20,18,1) - 2*pc(20,14,1) - 12*pc(13,16,1) +&
& 3*pc(13,18,1) + 3*pc(13,14,1) - 12*pc(15,16,1) + 3*pc(15,18,1) + 3*pc(15,14,1&
&)

sh(15,15,1) = 16*pc(16,16,1) - 4*pc(16,18,1) - 4*pc(16,14,1) - 4*pc(18,16,1) +&
& pc(18,18,1) + pc(18,14,1) - 4*pc(14,16,1) + pc(14,18,1) + pc(14,14,1)

sh(16,15,1) = 16*pc(17,16,1) - 4*pc(17,18,1) - 4*pc(17,14,1) - 4*pc(12,16,1) +&
& pc(12,18,1) + pc(12,14,1) - 4*pc(19,16,1) + pc(19,18,1) + pc(19,14,1)

sh(1,16,1) = 4*pc(1,17,1) - pc(1,12,1) - pc(1,19,1)

sh(2,16,1) = 4*pc(2,17,1) - pc(2,12,1) - pc(2,19,1)

sh(3,16,1) = 4*pc(3,17,1) - pc(3,12,1) - pc(3,19,1)

sh(4,16,1) = 4*pc(4,17,1) - pc(4,12,1) - pc(4,19,1)

sh(5,16,1) = 4*pc(8,17,1) - pc(8,12,1) - pc(8,19,1)

sh(6,16,1) = 4*pc(9,17,1) - pc(9,12,1) - pc(9,19,1)

sh(7,16,1) = 4*pc(10,17,1) - pc(10,12,1) - pc(10,19,1)

sh(8,16,1) = 4*pc(5,17,1) - pc(5,12,1) - pc(5,19,1) - 4*pc(6,17,1) + pc(6,12,1&
&) + pc(6,19,1)

sh(9,16,1) = 8*pc(7,17,1) - 2*pc(7,12,1) - 2*pc(7,19,1) - 4*pc(5,17,1) + pc(5,&
&12,1) + pc(5,19,1) - 4*pc(6,17,1) + pc(6,12,1) + pc(6,19,1)

sh(10,16,1) = 4*pc(11,17,1) - pc(11,12,1) - pc(11,19,1)

sh(11,16,1) = 4*pc(13,17,1) - pc(13,12,1) - pc(13,19,1) - 4*pc(15,17,1) + pc(1&
&5,12,1) + pc(15,19,1)

sh(12,16,1) = 4*pc(18,17,1) - pc(18,12,1) - pc(18,19,1) - 12*pc(14,17,1) + 3*p&
&c(14,12,1) + 3*pc(14,19,1)

sh(13,16,1) = 12*pc(12,17,1) - 3*pc(12,12,1) - 3*pc(12,19,1) - 4*pc(19,17,1) +&
& pc(19,12,1) + pc(19,19,1)

sh(14,16,1) = 8*pc(20,17,1) - 2*pc(20,12,1) - 2*pc(20,19,1) - 12*pc(13,17,1) +&
& 3*pc(13,12,1) + 3*pc(13,19,1) - 12*pc(15,17,1) + 3*pc(15,12,1) + 3*pc(15,19,1&
&)

sh(15,16,1) = 16*pc(16,17,1) - 4*pc(16,12,1) - 4*pc(16,19,1) - 4*pc(18,17,1) +&
& pc(18,12,1) + pc(18,19,1) - 4*pc(14,17,1) + pc(14,12,1) + pc(14,19,1)

sh(16,16,1) = 16*pc(17,17,1) - 4*pc(17,12,1) - 4*pc(17,19,1) - 4*pc(12,17,1) +&
& pc(12,12,1) + pc(12,19,1) - 4*pc(19,17,1) + pc(19,12,1) + pc(19,19,1)

sh(1,1,2) = pc(1,1,2)

sh(2,1,2) = pc(2,1,2)

sh(3,1,2) = pc(3,1,2)

sh(4,1,2) = pc(4,1,2)

sh(5,1,2) = pc(8,1,2)

sh(6,1,2) = pc(9,1,2)

sh(7,1,2) = pc(10,1,2)

sh(8,1,2) = pc(5,1,2) - pc(6,1,2)

sh(9,1,2) = 2*pc(7,1,2) - pc(5,1,2) - pc(6,1,2)

sh(10,1,2) = pc(11,1,2)

sh(11,1,2) = pc(13,1,2) - pc(15,1,2)

sh(12,1,2) = pc(18,1,2) - 3*pc(14,1,2)

sh(13,1,2) = 3*pc(12,1,2) - pc(19,1,2)

sh(14,1,2) = 2*pc(20,1,2) - 3*pc(13,1,2) - 3*pc(15,1,2)

sh(15,1,2) = 4*pc(16,1,2) - pc(18,1,2) - pc(14,1,2)

sh(16,1,2) = 4*pc(17,1,2) - pc(12,1,2) - pc(19,1,2)

sh(1,2,2) = pc(1,2,2)

sh(2,2,2) = pc(2,2,2)

sh(3,2,2) = pc(3,2,2)

sh(4,2,2) = pc(4,2,2)

sh(5,2,2) = pc(8,2,2)

sh(6,2,2) = pc(9,2,2)

sh(7,2,2) = pc(10,2,2)

sh(8,2,2) = pc(5,2,2) - pc(6,2,2)

sh(9,2,2) = 2*pc(7,2,2) - pc(5,2,2) - pc(6,2,2)

sh(10,2,2) = pc(11,2,2)

sh(11,2,2) = pc(13,2,2) - pc(15,2,2)

sh(12,2,2) = pc(18,2,2) - 3*pc(14,2,2)

sh(13,2,2) = 3*pc(12,2,2) - pc(19,2,2)

sh(14,2,2) = 2*pc(20,2,2) - 3*pc(13,2,2) - 3*pc(15,2,2)

sh(15,2,2) = 4*pc(16,2,2) - pc(18,2,2) - pc(14,2,2)

sh(16,2,2) = 4*pc(17,2,2) - pc(12,2,2) - pc(19,2,2)

sh(1,3,2) = pc(1,3,2)

sh(2,3,2) = pc(2,3,2)

sh(3,3,2) = pc(3,3,2)

sh(4,3,2) = pc(4,3,2)

sh(5,3,2) = pc(8,3,2)

sh(6,3,2) = pc(9,3,2)

sh(7,3,2) = pc(10,3,2)

sh(8,3,2) = pc(5,3,2) - pc(6,3,2)

sh(9,3,2) = 2*pc(7,3,2) - pc(5,3,2) - pc(6,3,2)

sh(10,3,2) = pc(11,3,2)

sh(11,3,2) = pc(13,3,2) - pc(15,3,2)

sh(12,3,2) = pc(18,3,2) - 3*pc(14,3,2)

sh(13,3,2) = 3*pc(12,3,2) - pc(19,3,2)

sh(14,3,2) = 2*pc(20,3,2) - 3*pc(13,3,2) - 3*pc(15,3,2)

sh(15,3,2) = 4*pc(16,3,2) - pc(18,3,2) - pc(14,3,2)

sh(16,3,2) = 4*pc(17,3,2) - pc(12,3,2) - pc(19,3,2)

sh(1,4,2) = pc(1,4,2)

sh(2,4,2) = pc(2,4,2)

sh(3,4,2) = pc(3,4,2)

sh(4,4,2) = pc(4,4,2)

sh(5,4,2) = pc(8,4,2)

sh(6,4,2) = pc(9,4,2)

sh(7,4,2) = pc(10,4,2)

sh(8,4,2) = pc(5,4,2) - pc(6,4,2)

sh(9,4,2) = 2*pc(7,4,2) - pc(5,4,2) - pc(6,4,2)

sh(10,4,2) = pc(11,4,2)

sh(11,4,2) = pc(13,4,2) - pc(15,4,2)

sh(12,4,2) = pc(18,4,2) - 3*pc(14,4,2)

sh(13,4,2) = 3*pc(12,4,2) - pc(19,4,2)

sh(14,4,2) = 2*pc(20,4,2) - 3*pc(13,4,2) - 3*pc(15,4,2)

sh(15,4,2) = 4*pc(16,4,2) - pc(18,4,2) - pc(14,4,2)

sh(16,4,2) = 4*pc(17,4,2) - pc(12,4,2) - pc(19,4,2)

sh(1,5,2) = pc(1,8,2)

sh(2,5,2) = pc(2,8,2)

sh(3,5,2) = pc(3,8,2)

sh(4,5,2) = pc(4,8,2)

sh(5,5,2) = pc(8,8,2)

sh(6,5,2) = pc(9,8,2)

sh(7,5,2) = pc(10,8,2)

sh(8,5,2) = pc(5,8,2) - pc(6,8,2)

sh(9,5,2) = 2*pc(7,8,2) - pc(5,8,2) - pc(6,8,2)

sh(10,5,2) = pc(11,8,2)

sh(11,5,2) = pc(13,8,2) - pc(15,8,2)

sh(12,5,2) = pc(18,8,2) - 3*pc(14,8,2)

sh(13,5,2) = 3*pc(12,8,2) - pc(19,8,2)

sh(14,5,2) = 2*pc(20,8,2) - 3*pc(13,8,2) - 3*pc(15,8,2)

sh(15,5,2) = 4*pc(16,8,2) - pc(18,8,2) - pc(14,8,2)

sh(16,5,2) = 4*pc(17,8,2) - pc(12,8,2) - pc(19,8,2)

sh(1,6,2) = pc(1,9,2)

sh(2,6,2) = pc(2,9,2)

sh(3,6,2) = pc(3,9,2)

sh(4,6,2) = pc(4,9,2)

sh(5,6,2) = pc(8,9,2)

sh(6,6,2) = pc(9,9,2)

sh(7,6,2) = pc(10,9,2)

sh(8,6,2) = pc(5,9,2) - pc(6,9,2)

sh(9,6,2) = 2*pc(7,9,2) - pc(5,9,2) - pc(6,9,2)

sh(10,6,2) = pc(11,9,2)

sh(11,6,2) = pc(13,9,2) - pc(15,9,2)

sh(12,6,2) = pc(18,9,2) - 3*pc(14,9,2)

sh(13,6,2) = 3*pc(12,9,2) - pc(19,9,2)

sh(14,6,2) = 2*pc(20,9,2) - 3*pc(13,9,2) - 3*pc(15,9,2)

sh(15,6,2) = 4*pc(16,9,2) - pc(18,9,2) - pc(14,9,2)

sh(16,6,2) = 4*pc(17,9,2) - pc(12,9,2) - pc(19,9,2)

sh(1,7,2) = pc(1,10,2)

sh(2,7,2) = pc(2,10,2)

sh(3,7,2) = pc(3,10,2)

sh(4,7,2) = pc(4,10,2)

sh(5,7,2) = pc(8,10,2)

sh(6,7,2) = pc(9,10,2)

sh(7,7,2) = pc(10,10,2)

sh(8,7,2) = pc(5,10,2) - pc(6,10,2)

sh(9,7,2) = 2*pc(7,10,2) - pc(5,10,2) - pc(6,10,2)

sh(10,7,2) = pc(11,10,2)

sh(11,7,2) = pc(13,10,2) - pc(15,10,2)

sh(12,7,2) = pc(18,10,2) - 3*pc(14,10,2)

sh(13,7,2) = 3*pc(12,10,2) - pc(19,10,2)

sh(14,7,2) = 2*pc(20,10,2) - 3*pc(13,10,2) - 3*pc(15,10,2)

sh(15,7,2) = 4*pc(16,10,2) - pc(18,10,2) - pc(14,10,2)

sh(16,7,2) = 4*pc(17,10,2) - pc(12,10,2) - pc(19,10,2)

sh(1,8,2) = pc(1,5,2) - pc(1,6,2)

sh(2,8,2) = pc(2,5,2) - pc(2,6,2)

sh(3,8,2) = pc(3,5,2) - pc(3,6,2)

sh(4,8,2) = pc(4,5,2) - pc(4,6,2)

sh(5,8,2) = pc(8,5,2) - pc(8,6,2)

sh(6,8,2) = pc(9,5,2) - pc(9,6,2)

sh(7,8,2) = pc(10,5,2) - pc(10,6,2)

sh(8,8,2) = pc(5,5,2) - pc(5,6,2) - pc(6,5,2) + pc(6,6,2)

sh(9,8,2) = 2*pc(7,5,2) - 2*pc(7,6,2) - pc(5,5,2) + pc(5,6,2) - pc(6,5,2) + pc&
&(6,6,2)

sh(10,8,2) = pc(11,5,2) - pc(11,6,2)

sh(11,8,2) = pc(13,5,2) - pc(13,6,2) - pc(15,5,2) + pc(15,6,2)

sh(12,8,2) = pc(18,5,2) - pc(18,6,2) - 3*pc(14,5,2) + 3*pc(14,6,2)

sh(13,8,2) = 3*pc(12,5,2) - 3*pc(12,6,2) - pc(19,5,2) + pc(19,6,2)

sh(14,8,2) = 2*pc(20,5,2) - 2*pc(20,6,2) - 3*pc(13,5,2) + 3*pc(13,6,2) - 3*pc(&
&15,5,2) + 3*pc(15,6,2)

sh(15,8,2) = 4*pc(16,5,2) - 4*pc(16,6,2) - pc(18,5,2) + pc(18,6,2) - pc(14,5,2&
&) + pc(14,6,2)

sh(16,8,2) = 4*pc(17,5,2) - 4*pc(17,6,2) - pc(12,5,2) + pc(12,6,2) - pc(19,5,2&
&) + pc(19,6,2)

sh(1,9,2) = 2*pc(1,7,2) - pc(1,5,2) - pc(1,6,2)

sh(2,9,2) = 2*pc(2,7,2) - pc(2,5,2) - pc(2,6,2)

sh(3,9,2) = 2*pc(3,7,2) - pc(3,5,2) - pc(3,6,2)

sh(4,9,2) = 2*pc(4,7,2) - pc(4,5,2) - pc(4,6,2)

sh(5,9,2) = 2*pc(8,7,2) - pc(8,5,2) - pc(8,6,2)

sh(6,9,2) = 2*pc(9,7,2) - pc(9,5,2) - pc(9,6,2)

sh(7,9,2) = 2*pc(10,7,2) - pc(10,5,2) - pc(10,6,2)

sh(8,9,2) = 2*pc(5,7,2) - pc(5,5,2) - pc(5,6,2) - 2*pc(6,7,2) + pc(6,5,2) + pc&
&(6,6,2)

sh(9,9,2) = 4*pc(7,7,2) - 2*pc(7,5,2) - 2*pc(7,6,2) - 2*pc(5,7,2) + pc(5,5,2) &
&+ pc(5,6,2) - 2*pc(6,7,2) + pc(6,5,2) + pc(6,6,2)

sh(10,9,2) = 2*pc(11,7,2) - pc(11,5,2) - pc(11,6,2)

sh(11,9,2) = 2*pc(13,7,2) - pc(13,5,2) - pc(13,6,2) - 2*pc(15,7,2) + pc(15,5,2&
&) + pc(15,6,2)

sh(12,9,2) = 2*pc(18,7,2) - pc(18,5,2) - pc(18,6,2) - 6*pc(14,7,2) + 3*pc(14,5&
&,2) + 3*pc(14,6,2)

sh(13,9,2) = 6*pc(12,7,2) - 3*pc(12,5,2) - 3*pc(12,6,2) - 2*pc(19,7,2) + pc(19&
&,5,2) + pc(19,6,2)

sh(14,9,2) = 4*pc(20,7,2) - 2*pc(20,5,2) - 2*pc(20,6,2) - 6*pc(13,7,2) + 3*pc(&
&13,5,2) + 3*pc(13,6,2) - 6*pc(15,7,2) + 3*pc(15,5,2) + 3*pc(15,6,2)

sh(15,9,2) = 8*pc(16,7,2) - 4*pc(16,5,2) - 4*pc(16,6,2) - 2*pc(18,7,2) + pc(18&
&,5,2) + pc(18,6,2) - 2*pc(14,7,2) + pc(14,5,2) + pc(14,6,2)

sh(16,9,2) = 8*pc(17,7,2) - 4*pc(17,5,2) - 4*pc(17,6,2) - 2*pc(12,7,2) + pc(12&
&,5,2) + pc(12,6,2) - 2*pc(19,7,2) + pc(19,5,2) + pc(19,6,2)

sh(1,10,2) = pc(1,11,2)

sh(2,10,2) = pc(2,11,2)

sh(3,10,2) = pc(3,11,2)

sh(4,10,2) = pc(4,11,2)

sh(5,10,2) = pc(8,11,2)

sh(6,10,2) = pc(9,11,2)

sh(7,10,2) = pc(10,11,2)

sh(8,10,2) = pc(5,11,2) - pc(6,11,2)

sh(9,10,2) = 2*pc(7,11,2) - pc(5,11,2) - pc(6,11,2)

sh(10,10,2) = pc(11,11,2)

sh(11,10,2) = pc(13,11,2) - pc(15,11,2)

sh(12,10,2) = pc(18,11,2) - 3*pc(14,11,2)

sh(13,10,2) = 3*pc(12,11,2) - pc(19,11,2)

sh(14,10,2) = 2*pc(20,11,2) - 3*pc(13,11,2) - 3*pc(15,11,2)

sh(15,10,2) = 4*pc(16,11,2) - pc(18,11,2) - pc(14,11,2)

sh(16,10,2) = 4*pc(17,11,2) - pc(12,11,2) - pc(19,11,2)

sh(1,11,2) = pc(1,13,2) - pc(1,15,2)

sh(2,11,2) = pc(2,13,2) - pc(2,15,2)

sh(3,11,2) = pc(3,13,2) - pc(3,15,2)

sh(4,11,2) = pc(4,13,2) - pc(4,15,2)

sh(5,11,2) = pc(8,13,2) - pc(8,15,2)

sh(6,11,2) = pc(9,13,2) - pc(9,15,2)

sh(7,11,2) = pc(10,13,2) - pc(10,15,2)

sh(8,11,2) = pc(5,13,2) - pc(5,15,2) - pc(6,13,2) + pc(6,15,2)

sh(9,11,2) = 2*pc(7,13,2) - 2*pc(7,15,2) - pc(5,13,2) + pc(5,15,2) - pc(6,13,2&
&) + pc(6,15,2)

sh(10,11,2) = pc(11,13,2) - pc(11,15,2)

sh(11,11,2) = pc(13,13,2) - pc(13,15,2) - pc(15,13,2) + pc(15,15,2)

sh(12,11,2) = pc(18,13,2) - pc(18,15,2) - 3*pc(14,13,2) + 3*pc(14,15,2)

sh(13,11,2) = 3*pc(12,13,2) - 3*pc(12,15,2) - pc(19,13,2) + pc(19,15,2)

sh(14,11,2) = 2*pc(20,13,2) - 2*pc(20,15,2) - 3*pc(13,13,2) + 3*pc(13,15,2) - &
&3*pc(15,13,2) + 3*pc(15,15,2)

sh(15,11,2) = 4*pc(16,13,2) - 4*pc(16,15,2) - pc(18,13,2) + pc(18,15,2) - pc(1&
&4,13,2) + pc(14,15,2)

sh(16,11,2) = 4*pc(17,13,2) - 4*pc(17,15,2) - pc(12,13,2) + pc(12,15,2) - pc(1&
&9,13,2) + pc(19,15,2)

sh(1,12,2) = pc(1,18,2) - 3*pc(1,14,2)

sh(2,12,2) = pc(2,18,2) - 3*pc(2,14,2)

sh(3,12,2) = pc(3,18,2) - 3*pc(3,14,2)

sh(4,12,2) = pc(4,18,2) - 3*pc(4,14,2)

sh(5,12,2) = pc(8,18,2) - 3*pc(8,14,2)

sh(6,12,2) = pc(9,18,2) - 3*pc(9,14,2)

sh(7,12,2) = pc(10,18,2) - 3*pc(10,14,2)

sh(8,12,2) = pc(5,18,2) - 3*pc(5,14,2) - pc(6,18,2) + 3*pc(6,14,2)

sh(9,12,2) = 2*pc(7,18,2) - 6*pc(7,14,2) - pc(5,18,2) + 3*pc(5,14,2) - pc(6,18&
&,2) + 3*pc(6,14,2)

sh(10,12,2) = pc(11,18,2) - 3*pc(11,14,2)

sh(11,12,2) = pc(13,18,2) - 3*pc(13,14,2) - pc(15,18,2) + 3*pc(15,14,2)

sh(12,12,2) = pc(18,18,2) - 3*pc(18,14,2) - 3*pc(14,18,2) + 9*pc(14,14,2)

sh(13,12,2) = 3*pc(12,18,2) - 9*pc(12,14,2) - pc(19,18,2) + 3*pc(19,14,2)

sh(14,12,2) = 2*pc(20,18,2) - 6*pc(20,14,2) - 3*pc(13,18,2) + 9*pc(13,14,2) - &
&3*pc(15,18,2) + 9*pc(15,14,2)

sh(15,12,2) = 4*pc(16,18,2) - 12*pc(16,14,2) - pc(18,18,2) + 3*pc(18,14,2) - p&
&c(14,18,2) + 3*pc(14,14,2)

sh(16,12,2) = 4*pc(17,18,2) - 12*pc(17,14,2) - pc(12,18,2) + 3*pc(12,14,2) - p&
&c(19,18,2) + 3*pc(19,14,2)

sh(1,13,2) = 3*pc(1,12,2) - pc(1,19,2)

sh(2,13,2) = 3*pc(2,12,2) - pc(2,19,2)

sh(3,13,2) = 3*pc(3,12,2) - pc(3,19,2)

sh(4,13,2) = 3*pc(4,12,2) - pc(4,19,2)

sh(5,13,2) = 3*pc(8,12,2) - pc(8,19,2)

sh(6,13,2) = 3*pc(9,12,2) - pc(9,19,2)

sh(7,13,2) = 3*pc(10,12,2) - pc(10,19,2)

sh(8,13,2) = 3*pc(5,12,2) - pc(5,19,2) - 3*pc(6,12,2) + pc(6,19,2)

sh(9,13,2) = 6*pc(7,12,2) - 2*pc(7,19,2) - 3*pc(5,12,2) + pc(5,19,2) - 3*pc(6,&
&12,2) + pc(6,19,2)

sh(10,13,2) = 3*pc(11,12,2) - pc(11,19,2)

sh(11,13,2) = 3*pc(13,12,2) - pc(13,19,2) - 3*pc(15,12,2) + pc(15,19,2)

sh(12,13,2) = 3*pc(18,12,2) - pc(18,19,2) - 9*pc(14,12,2) + 3*pc(14,19,2)

sh(13,13,2) = 9*pc(12,12,2) - 3*pc(12,19,2) - 3*pc(19,12,2) + pc(19,19,2)

sh(14,13,2) = 6*pc(20,12,2) - 2*pc(20,19,2) - 9*pc(13,12,2) + 3*pc(13,19,2) - &
&9*pc(15,12,2) + 3*pc(15,19,2)

sh(15,13,2) = 12*pc(16,12,2) - 4*pc(16,19,2) - 3*pc(18,12,2) + pc(18,19,2) - 3&
&*pc(14,12,2) + pc(14,19,2)

sh(16,13,2) = 12*pc(17,12,2) - 4*pc(17,19,2) - 3*pc(12,12,2) + pc(12,19,2) - 3&
&*pc(19,12,2) + pc(19,19,2)

sh(1,14,2) = 2*pc(1,20,2) - 3*pc(1,13,2) - 3*pc(1,15,2)

sh(2,14,2) = 2*pc(2,20,2) - 3*pc(2,13,2) - 3*pc(2,15,2)

sh(3,14,2) = 2*pc(3,20,2) - 3*pc(3,13,2) - 3*pc(3,15,2)

sh(4,14,2) = 2*pc(4,20,2) - 3*pc(4,13,2) - 3*pc(4,15,2)

sh(5,14,2) = 2*pc(8,20,2) - 3*pc(8,13,2) - 3*pc(8,15,2)

sh(6,14,2) = 2*pc(9,20,2) - 3*pc(9,13,2) - 3*pc(9,15,2)

sh(7,14,2) = 2*pc(10,20,2) - 3*pc(10,13,2) - 3*pc(10,15,2)

sh(8,14,2) = 2*pc(5,20,2) - 3*pc(5,13,2) - 3*pc(5,15,2) - 2*pc(6,20,2) + 3*pc(&
&6,13,2) + 3*pc(6,15,2)

sh(9,14,2) = 4*pc(7,20,2) - 6*pc(7,13,2) - 6*pc(7,15,2) - 2*pc(5,20,2) + 3*pc(&
&5,13,2) + 3*pc(5,15,2) - 2*pc(6,20,2) + 3*pc(6,13,2) + 3*pc(6,15,2)

sh(10,14,2) = 2*pc(11,20,2) - 3*pc(11,13,2) - 3*pc(11,15,2)

sh(11,14,2) = 2*pc(13,20,2) - 3*pc(13,13,2) - 3*pc(13,15,2) - 2*pc(15,20,2) + &
&3*pc(15,13,2) + 3*pc(15,15,2)

sh(12,14,2) = 2*pc(18,20,2) - 3*pc(18,13,2) - 3*pc(18,15,2) - 6*pc(14,20,2) + &
&9*pc(14,13,2) + 9*pc(14,15,2)

sh(13,14,2) = 6*pc(12,20,2) - 9*pc(12,13,2) - 9*pc(12,15,2) - 2*pc(19,20,2) + &
&3*pc(19,13,2) + 3*pc(19,15,2)

sh(14,14,2) = 4*pc(20,20,2) - 6*pc(20,13,2) - 6*pc(20,15,2) - 6*pc(13,20,2) + &
&9*pc(13,13,2) + 9*pc(13,15,2) - 6*pc(15,20,2) + 9*pc(15,13,2) + 9*pc(15,15,2)

sh(15,14,2) = 8*pc(16,20,2) - 12*pc(16,13,2) - 12*pc(16,15,2) - 2*pc(18,20,2) &
&+ 3*pc(18,13,2) + 3*pc(18,15,2) - 2*pc(14,20,2) + 3*pc(14,13,2) + 3*pc(14,15,2&
&)

sh(16,14,2) = 8*pc(17,20,2) - 12*pc(17,13,2) - 12*pc(17,15,2) - 2*pc(12,20,2) &
&+ 3*pc(12,13,2) + 3*pc(12,15,2) - 2*pc(19,20,2) + 3*pc(19,13,2) + 3*pc(19,15,2&
&)

sh(1,15,2) = 4*pc(1,16,2) - pc(1,18,2) - pc(1,14,2)

sh(2,15,2) = 4*pc(2,16,2) - pc(2,18,2) - pc(2,14,2)

sh(3,15,2) = 4*pc(3,16,2) - pc(3,18,2) - pc(3,14,2)

sh(4,15,2) = 4*pc(4,16,2) - pc(4,18,2) - pc(4,14,2)

sh(5,15,2) = 4*pc(8,16,2) - pc(8,18,2) - pc(8,14,2)

sh(6,15,2) = 4*pc(9,16,2) - pc(9,18,2) - pc(9,14,2)

sh(7,15,2) = 4*pc(10,16,2) - pc(10,18,2) - pc(10,14,2)

sh(8,15,2) = 4*pc(5,16,2) - pc(5,18,2) - pc(5,14,2) - 4*pc(6,16,2) + pc(6,18,2&
&) + pc(6,14,2)

sh(9,15,2) = 8*pc(7,16,2) - 2*pc(7,18,2) - 2*pc(7,14,2) - 4*pc(5,16,2) + pc(5,&
&18,2) + pc(5,14,2) - 4*pc(6,16,2) + pc(6,18,2) + pc(6,14,2)

sh(10,15,2) = 4*pc(11,16,2) - pc(11,18,2) - pc(11,14,2)

sh(11,15,2) = 4*pc(13,16,2) - pc(13,18,2) - pc(13,14,2) - 4*pc(15,16,2) + pc(1&
&5,18,2) + pc(15,14,2)

sh(12,15,2) = 4*pc(18,16,2) - pc(18,18,2) - pc(18,14,2) - 12*pc(14,16,2) + 3*p&
&c(14,18,2) + 3*pc(14,14,2)

sh(13,15,2) = 12*pc(12,16,2) - 3*pc(12,18,2) - 3*pc(12,14,2) - 4*pc(19,16,2) +&
& pc(19,18,2) + pc(19,14,2)

sh(14,15,2) = 8*pc(20,16,2) - 2*pc(20,18,2) - 2*pc(20,14,2) - 12*pc(13,16,2) +&
& 3*pc(13,18,2) + 3*pc(13,14,2) - 12*pc(15,16,2) + 3*pc(15,18,2) + 3*pc(15,14,2&
&)

sh(15,15,2) = 16*pc(16,16,2) - 4*pc(16,18,2) - 4*pc(16,14,2) - 4*pc(18,16,2) +&
& pc(18,18,2) + pc(18,14,2) - 4*pc(14,16,2) + pc(14,18,2) + pc(14,14,2)

sh(16,15,2) = 16*pc(17,16,2) - 4*pc(17,18,2) - 4*pc(17,14,2) - 4*pc(12,16,2) +&
& pc(12,18,2) + pc(12,14,2) - 4*pc(19,16,2) + pc(19,18,2) + pc(19,14,2)

sh(1,16,2) = 4*pc(1,17,2) - pc(1,12,2) - pc(1,19,2)

sh(2,16,2) = 4*pc(2,17,2) - pc(2,12,2) - pc(2,19,2)

sh(3,16,2) = 4*pc(3,17,2) - pc(3,12,2) - pc(3,19,2)

sh(4,16,2) = 4*pc(4,17,2) - pc(4,12,2) - pc(4,19,2)

sh(5,16,2) = 4*pc(8,17,2) - pc(8,12,2) - pc(8,19,2)

sh(6,16,2) = 4*pc(9,17,2) - pc(9,12,2) - pc(9,19,2)

sh(7,16,2) = 4*pc(10,17,2) - pc(10,12,2) - pc(10,19,2)

sh(8,16,2) = 4*pc(5,17,2) - pc(5,12,2) - pc(5,19,2) - 4*pc(6,17,2) + pc(6,12,2&
&) + pc(6,19,2)

sh(9,16,2) = 8*pc(7,17,2) - 2*pc(7,12,2) - 2*pc(7,19,2) - 4*pc(5,17,2) + pc(5,&
&12,2) + pc(5,19,2) - 4*pc(6,17,2) + pc(6,12,2) + pc(6,19,2)

sh(10,16,2) = 4*pc(11,17,2) - pc(11,12,2) - pc(11,19,2)

sh(11,16,2) = 4*pc(13,17,2) - pc(13,12,2) - pc(13,19,2) - 4*pc(15,17,2) + pc(1&
&5,12,2) + pc(15,19,2)

sh(12,16,2) = 4*pc(18,17,2) - pc(18,12,2) - pc(18,19,2) - 12*pc(14,17,2) + 3*p&
&c(14,12,2) + 3*pc(14,19,2)

sh(13,16,2) = 12*pc(12,17,2) - 3*pc(12,12,2) - 3*pc(12,19,2) - 4*pc(19,17,2) +&
& pc(19,12,2) + pc(19,19,2)

sh(14,16,2) = 8*pc(20,17,2) - 2*pc(20,12,2) - 2*pc(20,19,2) - 12*pc(13,17,2) +&
& 3*pc(13,12,2) + 3*pc(13,19,2) - 12*pc(15,17,2) + 3*pc(15,12,2) + 3*pc(15,19,2&
&)

sh(15,16,2) = 16*pc(16,17,2) - 4*pc(16,12,2) - 4*pc(16,19,2) - 4*pc(18,17,2) +&
& pc(18,12,2) + pc(18,19,2) - 4*pc(14,17,2) + pc(14,12,2) + pc(14,19,2)

sh(16,16,2) = 16*pc(17,17,2) - 4*pc(17,12,2) - 4*pc(17,19,2) - 4*pc(12,17,2) +&
& pc(12,12,2) + pc(12,19,2) - 4*pc(19,17,2) + pc(19,12,2) + pc(19,19,2)

sh(1,1,3) = pc(1,1,3)

sh(2,1,3) = pc(2,1,3)

sh(3,1,3) = pc(3,1,3)

sh(4,1,3) = pc(4,1,3)

sh(5,1,3) = pc(8,1,3)

sh(6,1,3) = pc(9,1,3)

sh(7,1,3) = pc(10,1,3)

sh(8,1,3) = pc(5,1,3) - pc(6,1,3)

sh(9,1,3) = 2*pc(7,1,3) - pc(5,1,3) - pc(6,1,3)

sh(10,1,3) = pc(11,1,3)

sh(11,1,3) = pc(13,1,3) - pc(15,1,3)

sh(12,1,3) = pc(18,1,3) - 3*pc(14,1,3)

sh(13,1,3) = 3*pc(12,1,3) - pc(19,1,3)

sh(14,1,3) = 2*pc(20,1,3) - 3*pc(13,1,3) - 3*pc(15,1,3)

sh(15,1,3) = 4*pc(16,1,3) - pc(18,1,3) - pc(14,1,3)

sh(16,1,3) = 4*pc(17,1,3) - pc(12,1,3) - pc(19,1,3)

sh(1,2,3) = pc(1,2,3)

sh(2,2,3) = pc(2,2,3)

sh(3,2,3) = pc(3,2,3)

sh(4,2,3) = pc(4,2,3)

sh(5,2,3) = pc(8,2,3)

sh(6,2,3) = pc(9,2,3)

sh(7,2,3) = pc(10,2,3)

sh(8,2,3) = pc(5,2,3) - pc(6,2,3)

sh(9,2,3) = 2*pc(7,2,3) - pc(5,2,3) - pc(6,2,3)

sh(10,2,3) = pc(11,2,3)

sh(11,2,3) = pc(13,2,3) - pc(15,2,3)

sh(12,2,3) = pc(18,2,3) - 3*pc(14,2,3)

sh(13,2,3) = 3*pc(12,2,3) - pc(19,2,3)

sh(14,2,3) = 2*pc(20,2,3) - 3*pc(13,2,3) - 3*pc(15,2,3)

sh(15,2,3) = 4*pc(16,2,3) - pc(18,2,3) - pc(14,2,3)

sh(16,2,3) = 4*pc(17,2,3) - pc(12,2,3) - pc(19,2,3)

sh(1,3,3) = pc(1,3,3)

sh(2,3,3) = pc(2,3,3)

sh(3,3,3) = pc(3,3,3)

sh(4,3,3) = pc(4,3,3)

sh(5,3,3) = pc(8,3,3)

sh(6,3,3) = pc(9,3,3)

sh(7,3,3) = pc(10,3,3)

sh(8,3,3) = pc(5,3,3) - pc(6,3,3)

sh(9,3,3) = 2*pc(7,3,3) - pc(5,3,3) - pc(6,3,3)

sh(10,3,3) = pc(11,3,3)

sh(11,3,3) = pc(13,3,3) - pc(15,3,3)

sh(12,3,3) = pc(18,3,3) - 3*pc(14,3,3)

sh(13,3,3) = 3*pc(12,3,3) - pc(19,3,3)

sh(14,3,3) = 2*pc(20,3,3) - 3*pc(13,3,3) - 3*pc(15,3,3)

sh(15,3,3) = 4*pc(16,3,3) - pc(18,3,3) - pc(14,3,3)

sh(16,3,3) = 4*pc(17,3,3) - pc(12,3,3) - pc(19,3,3)

sh(1,4,3) = pc(1,4,3)

sh(2,4,3) = pc(2,4,3)

sh(3,4,3) = pc(3,4,3)

sh(4,4,3) = pc(4,4,3)

sh(5,4,3) = pc(8,4,3)

sh(6,4,3) = pc(9,4,3)

sh(7,4,3) = pc(10,4,3)

sh(8,4,3) = pc(5,4,3) - pc(6,4,3)

sh(9,4,3) = 2*pc(7,4,3) - pc(5,4,3) - pc(6,4,3)

sh(10,4,3) = pc(11,4,3)

sh(11,4,3) = pc(13,4,3) - pc(15,4,3)

sh(12,4,3) = pc(18,4,3) - 3*pc(14,4,3)

sh(13,4,3) = 3*pc(12,4,3) - pc(19,4,3)

sh(14,4,3) = 2*pc(20,4,3) - 3*pc(13,4,3) - 3*pc(15,4,3)

sh(15,4,3) = 4*pc(16,4,3) - pc(18,4,3) - pc(14,4,3)

sh(16,4,3) = 4*pc(17,4,3) - pc(12,4,3) - pc(19,4,3)

sh(1,5,3) = pc(1,8,3)

sh(2,5,3) = pc(2,8,3)

sh(3,5,3) = pc(3,8,3)

sh(4,5,3) = pc(4,8,3)

sh(5,5,3) = pc(8,8,3)

sh(6,5,3) = pc(9,8,3)

sh(7,5,3) = pc(10,8,3)

sh(8,5,3) = pc(5,8,3) - pc(6,8,3)

sh(9,5,3) = 2*pc(7,8,3) - pc(5,8,3) - pc(6,8,3)

sh(10,5,3) = pc(11,8,3)

sh(11,5,3) = pc(13,8,3) - pc(15,8,3)

sh(12,5,3) = pc(18,8,3) - 3*pc(14,8,3)

sh(13,5,3) = 3*pc(12,8,3) - pc(19,8,3)

sh(14,5,3) = 2*pc(20,8,3) - 3*pc(13,8,3) - 3*pc(15,8,3)

sh(15,5,3) = 4*pc(16,8,3) - pc(18,8,3) - pc(14,8,3)

sh(16,5,3) = 4*pc(17,8,3) - pc(12,8,3) - pc(19,8,3)

sh(1,6,3) = pc(1,9,3)

sh(2,6,3) = pc(2,9,3)

sh(3,6,3) = pc(3,9,3)

sh(4,6,3) = pc(4,9,3)

sh(5,6,3) = pc(8,9,3)

sh(6,6,3) = pc(9,9,3)

sh(7,6,3) = pc(10,9,3)

sh(8,6,3) = pc(5,9,3) - pc(6,9,3)

sh(9,6,3) = 2*pc(7,9,3) - pc(5,9,3) - pc(6,9,3)

sh(10,6,3) = pc(11,9,3)

sh(11,6,3) = pc(13,9,3) - pc(15,9,3)

sh(12,6,3) = pc(18,9,3) - 3*pc(14,9,3)

sh(13,6,3) = 3*pc(12,9,3) - pc(19,9,3)

sh(14,6,3) = 2*pc(20,9,3) - 3*pc(13,9,3) - 3*pc(15,9,3)

sh(15,6,3) = 4*pc(16,9,3) - pc(18,9,3) - pc(14,9,3)

sh(16,6,3) = 4*pc(17,9,3) - pc(12,9,3) - pc(19,9,3)

sh(1,7,3) = pc(1,10,3)

sh(2,7,3) = pc(2,10,3)

sh(3,7,3) = pc(3,10,3)

sh(4,7,3) = pc(4,10,3)

sh(5,7,3) = pc(8,10,3)

sh(6,7,3) = pc(9,10,3)

sh(7,7,3) = pc(10,10,3)

sh(8,7,3) = pc(5,10,3) - pc(6,10,3)

sh(9,7,3) = 2*pc(7,10,3) - pc(5,10,3) - pc(6,10,3)

sh(10,7,3) = pc(11,10,3)

sh(11,7,3) = pc(13,10,3) - pc(15,10,3)

sh(12,7,3) = pc(18,10,3) - 3*pc(14,10,3)

sh(13,7,3) = 3*pc(12,10,3) - pc(19,10,3)

sh(14,7,3) = 2*pc(20,10,3) - 3*pc(13,10,3) - 3*pc(15,10,3)

sh(15,7,3) = 4*pc(16,10,3) - pc(18,10,3) - pc(14,10,3)

sh(16,7,3) = 4*pc(17,10,3) - pc(12,10,3) - pc(19,10,3)

sh(1,8,3) = pc(1,5,3) - pc(1,6,3)

sh(2,8,3) = pc(2,5,3) - pc(2,6,3)

sh(3,8,3) = pc(3,5,3) - pc(3,6,3)

sh(4,8,3) = pc(4,5,3) - pc(4,6,3)

sh(5,8,3) = pc(8,5,3) - pc(8,6,3)

sh(6,8,3) = pc(9,5,3) - pc(9,6,3)

sh(7,8,3) = pc(10,5,3) - pc(10,6,3)

sh(8,8,3) = pc(5,5,3) - pc(5,6,3) - pc(6,5,3) + pc(6,6,3)

sh(9,8,3) = 2*pc(7,5,3) - 2*pc(7,6,3) - pc(5,5,3) + pc(5,6,3) - pc(6,5,3) + pc&
&(6,6,3)

sh(10,8,3) = pc(11,5,3) - pc(11,6,3)

sh(11,8,3) = pc(13,5,3) - pc(13,6,3) - pc(15,5,3) + pc(15,6,3)

sh(12,8,3) = pc(18,5,3) - pc(18,6,3) - 3*pc(14,5,3) + 3*pc(14,6,3)

sh(13,8,3) = 3*pc(12,5,3) - 3*pc(12,6,3) - pc(19,5,3) + pc(19,6,3)

sh(14,8,3) = 2*pc(20,5,3) - 2*pc(20,6,3) - 3*pc(13,5,3) + 3*pc(13,6,3) - 3*pc(&
&15,5,3) + 3*pc(15,6,3)

sh(15,8,3) = 4*pc(16,5,3) - 4*pc(16,6,3) - pc(18,5,3) + pc(18,6,3) - pc(14,5,3&
&) + pc(14,6,3)

sh(16,8,3) = 4*pc(17,5,3) - 4*pc(17,6,3) - pc(12,5,3) + pc(12,6,3) - pc(19,5,3&
&) + pc(19,6,3)

sh(1,9,3) = 2*pc(1,7,3) - pc(1,5,3) - pc(1,6,3)

sh(2,9,3) = 2*pc(2,7,3) - pc(2,5,3) - pc(2,6,3)

sh(3,9,3) = 2*pc(3,7,3) - pc(3,5,3) - pc(3,6,3)

sh(4,9,3) = 2*pc(4,7,3) - pc(4,5,3) - pc(4,6,3)

sh(5,9,3) = 2*pc(8,7,3) - pc(8,5,3) - pc(8,6,3)

sh(6,9,3) = 2*pc(9,7,3) - pc(9,5,3) - pc(9,6,3)

sh(7,9,3) = 2*pc(10,7,3) - pc(10,5,3) - pc(10,6,3)

sh(8,9,3) = 2*pc(5,7,3) - pc(5,5,3) - pc(5,6,3) - 2*pc(6,7,3) + pc(6,5,3) + pc&
&(6,6,3)

sh(9,9,3) = 4*pc(7,7,3) - 2*pc(7,5,3) - 2*pc(7,6,3) - 2*pc(5,7,3) + pc(5,5,3) &
&+ pc(5,6,3) - 2*pc(6,7,3) + pc(6,5,3) + pc(6,6,3)

sh(10,9,3) = 2*pc(11,7,3) - pc(11,5,3) - pc(11,6,3)

sh(11,9,3) = 2*pc(13,7,3) - pc(13,5,3) - pc(13,6,3) - 2*pc(15,7,3) + pc(15,5,3&
&) + pc(15,6,3)

sh(12,9,3) = 2*pc(18,7,3) - pc(18,5,3) - pc(18,6,3) - 6*pc(14,7,3) + 3*pc(14,5&
&,3) + 3*pc(14,6,3)

sh(13,9,3) = 6*pc(12,7,3) - 3*pc(12,5,3) - 3*pc(12,6,3) - 2*pc(19,7,3) + pc(19&
&,5,3) + pc(19,6,3)

sh(14,9,3) = 4*pc(20,7,3) - 2*pc(20,5,3) - 2*pc(20,6,3) - 6*pc(13,7,3) + 3*pc(&
&13,5,3) + 3*pc(13,6,3) - 6*pc(15,7,3) + 3*pc(15,5,3) + 3*pc(15,6,3)

sh(15,9,3) = 8*pc(16,7,3) - 4*pc(16,5,3) - 4*pc(16,6,3) - 2*pc(18,7,3) + pc(18&
&,5,3) + pc(18,6,3) - 2*pc(14,7,3) + pc(14,5,3) + pc(14,6,3)

sh(16,9,3) = 8*pc(17,7,3) - 4*pc(17,5,3) - 4*pc(17,6,3) - 2*pc(12,7,3) + pc(12&
&,5,3) + pc(12,6,3) - 2*pc(19,7,3) + pc(19,5,3) + pc(19,6,3)

sh(1,10,3) = pc(1,11,3)

sh(2,10,3) = pc(2,11,3)

sh(3,10,3) = pc(3,11,3)

sh(4,10,3) = pc(4,11,3)

sh(5,10,3) = pc(8,11,3)

sh(6,10,3) = pc(9,11,3)

sh(7,10,3) = pc(10,11,3)

sh(8,10,3) = pc(5,11,3) - pc(6,11,3)

sh(9,10,3) = 2*pc(7,11,3) - pc(5,11,3) - pc(6,11,3)

sh(10,10,3) = pc(11,11,3)

sh(11,10,3) = pc(13,11,3) - pc(15,11,3)

sh(12,10,3) = pc(18,11,3) - 3*pc(14,11,3)

sh(13,10,3) = 3*pc(12,11,3) - pc(19,11,3)

sh(14,10,3) = 2*pc(20,11,3) - 3*pc(13,11,3) - 3*pc(15,11,3)

sh(15,10,3) = 4*pc(16,11,3) - pc(18,11,3) - pc(14,11,3)

sh(16,10,3) = 4*pc(17,11,3) - pc(12,11,3) - pc(19,11,3)

sh(1,11,3) = pc(1,13,3) - pc(1,15,3)

sh(2,11,3) = pc(2,13,3) - pc(2,15,3)

sh(3,11,3) = pc(3,13,3) - pc(3,15,3)

sh(4,11,3) = pc(4,13,3) - pc(4,15,3)

sh(5,11,3) = pc(8,13,3) - pc(8,15,3)

sh(6,11,3) = pc(9,13,3) - pc(9,15,3)

sh(7,11,3) = pc(10,13,3) - pc(10,15,3)

sh(8,11,3) = pc(5,13,3) - pc(5,15,3) - pc(6,13,3) + pc(6,15,3)

sh(9,11,3) = 2*pc(7,13,3) - 2*pc(7,15,3) - pc(5,13,3) + pc(5,15,3) - pc(6,13,3&
&) + pc(6,15,3)

sh(10,11,3) = pc(11,13,3) - pc(11,15,3)

sh(11,11,3) = pc(13,13,3) - pc(13,15,3) - pc(15,13,3) + pc(15,15,3)

sh(12,11,3) = pc(18,13,3) - pc(18,15,3) - 3*pc(14,13,3) + 3*pc(14,15,3)

sh(13,11,3) = 3*pc(12,13,3) - 3*pc(12,15,3) - pc(19,13,3) + pc(19,15,3)

sh(14,11,3) = 2*pc(20,13,3) - 2*pc(20,15,3) - 3*pc(13,13,3) + 3*pc(13,15,3) - &
&3*pc(15,13,3) + 3*pc(15,15,3)

sh(15,11,3) = 4*pc(16,13,3) - 4*pc(16,15,3) - pc(18,13,3) + pc(18,15,3) - pc(1&
&4,13,3) + pc(14,15,3)

sh(16,11,3) = 4*pc(17,13,3) - 4*pc(17,15,3) - pc(12,13,3) + pc(12,15,3) - pc(1&
&9,13,3) + pc(19,15,3)

sh(1,12,3) = pc(1,18,3) - 3*pc(1,14,3)

sh(2,12,3) = pc(2,18,3) - 3*pc(2,14,3)

sh(3,12,3) = pc(3,18,3) - 3*pc(3,14,3)

sh(4,12,3) = pc(4,18,3) - 3*pc(4,14,3)

sh(5,12,3) = pc(8,18,3) - 3*pc(8,14,3)

sh(6,12,3) = pc(9,18,3) - 3*pc(9,14,3)

sh(7,12,3) = pc(10,18,3) - 3*pc(10,14,3)

sh(8,12,3) = pc(5,18,3) - 3*pc(5,14,3) - pc(6,18,3) + 3*pc(6,14,3)

sh(9,12,3) = 2*pc(7,18,3) - 6*pc(7,14,3) - pc(5,18,3) + 3*pc(5,14,3) - pc(6,18&
&,3) + 3*pc(6,14,3)

sh(10,12,3) = pc(11,18,3) - 3*pc(11,14,3)

sh(11,12,3) = pc(13,18,3) - 3*pc(13,14,3) - pc(15,18,3) + 3*pc(15,14,3)

sh(12,12,3) = pc(18,18,3) - 3*pc(18,14,3) - 3*pc(14,18,3) + 9*pc(14,14,3)

sh(13,12,3) = 3*pc(12,18,3) - 9*pc(12,14,3) - pc(19,18,3) + 3*pc(19,14,3)

sh(14,12,3) = 2*pc(20,18,3) - 6*pc(20,14,3) - 3*pc(13,18,3) + 9*pc(13,14,3) - &
&3*pc(15,18,3) + 9*pc(15,14,3)

sh(15,12,3) = 4*pc(16,18,3) - 12*pc(16,14,3) - pc(18,18,3) + 3*pc(18,14,3) - p&
&c(14,18,3) + 3*pc(14,14,3)

sh(16,12,3) = 4*pc(17,18,3) - 12*pc(17,14,3) - pc(12,18,3) + 3*pc(12,14,3) - p&
&c(19,18,3) + 3*pc(19,14,3)

sh(1,13,3) = 3*pc(1,12,3) - pc(1,19,3)

sh(2,13,3) = 3*pc(2,12,3) - pc(2,19,3)

sh(3,13,3) = 3*pc(3,12,3) - pc(3,19,3)

sh(4,13,3) = 3*pc(4,12,3) - pc(4,19,3)

sh(5,13,3) = 3*pc(8,12,3) - pc(8,19,3)

sh(6,13,3) = 3*pc(9,12,3) - pc(9,19,3)

sh(7,13,3) = 3*pc(10,12,3) - pc(10,19,3)

sh(8,13,3) = 3*pc(5,12,3) - pc(5,19,3) - 3*pc(6,12,3) + pc(6,19,3)

sh(9,13,3) = 6*pc(7,12,3) - 2*pc(7,19,3) - 3*pc(5,12,3) + pc(5,19,3) - 3*pc(6,&
&12,3) + pc(6,19,3)

sh(10,13,3) = 3*pc(11,12,3) - pc(11,19,3)

sh(11,13,3) = 3*pc(13,12,3) - pc(13,19,3) - 3*pc(15,12,3) + pc(15,19,3)

sh(12,13,3) = 3*pc(18,12,3) - pc(18,19,3) - 9*pc(14,12,3) + 3*pc(14,19,3)

sh(13,13,3) = 9*pc(12,12,3) - 3*pc(12,19,3) - 3*pc(19,12,3) + pc(19,19,3)

sh(14,13,3) = 6*pc(20,12,3) - 2*pc(20,19,3) - 9*pc(13,12,3) + 3*pc(13,19,3) - &
&9*pc(15,12,3) + 3*pc(15,19,3)

sh(15,13,3) = 12*pc(16,12,3) - 4*pc(16,19,3) - 3*pc(18,12,3) + pc(18,19,3) - 3&
&*pc(14,12,3) + pc(14,19,3)

sh(16,13,3) = 12*pc(17,12,3) - 4*pc(17,19,3) - 3*pc(12,12,3) + pc(12,19,3) - 3&
&*pc(19,12,3) + pc(19,19,3)

sh(1,14,3) = 2*pc(1,20,3) - 3*pc(1,13,3) - 3*pc(1,15,3)

sh(2,14,3) = 2*pc(2,20,3) - 3*pc(2,13,3) - 3*pc(2,15,3)

sh(3,14,3) = 2*pc(3,20,3) - 3*pc(3,13,3) - 3*pc(3,15,3)

sh(4,14,3) = 2*pc(4,20,3) - 3*pc(4,13,3) - 3*pc(4,15,3)

sh(5,14,3) = 2*pc(8,20,3) - 3*pc(8,13,3) - 3*pc(8,15,3)

sh(6,14,3) = 2*pc(9,20,3) - 3*pc(9,13,3) - 3*pc(9,15,3)

sh(7,14,3) = 2*pc(10,20,3) - 3*pc(10,13,3) - 3*pc(10,15,3)

sh(8,14,3) = 2*pc(5,20,3) - 3*pc(5,13,3) - 3*pc(5,15,3) - 2*pc(6,20,3) + 3*pc(&
&6,13,3) + 3*pc(6,15,3)

sh(9,14,3) = 4*pc(7,20,3) - 6*pc(7,13,3) - 6*pc(7,15,3) - 2*pc(5,20,3) + 3*pc(&
&5,13,3) + 3*pc(5,15,3) - 2*pc(6,20,3) + 3*pc(6,13,3) + 3*pc(6,15,3)

sh(10,14,3) = 2*pc(11,20,3) - 3*pc(11,13,3) - 3*pc(11,15,3)

sh(11,14,3) = 2*pc(13,20,3) - 3*pc(13,13,3) - 3*pc(13,15,3) - 2*pc(15,20,3) + &
&3*pc(15,13,3) + 3*pc(15,15,3)

sh(12,14,3) = 2*pc(18,20,3) - 3*pc(18,13,3) - 3*pc(18,15,3) - 6*pc(14,20,3) + &
&9*pc(14,13,3) + 9*pc(14,15,3)

sh(13,14,3) = 6*pc(12,20,3) - 9*pc(12,13,3) - 9*pc(12,15,3) - 2*pc(19,20,3) + &
&3*pc(19,13,3) + 3*pc(19,15,3)

sh(14,14,3) = 4*pc(20,20,3) - 6*pc(20,13,3) - 6*pc(20,15,3) - 6*pc(13,20,3) + &
&9*pc(13,13,3) + 9*pc(13,15,3) - 6*pc(15,20,3) + 9*pc(15,13,3) + 9*pc(15,15,3)

sh(15,14,3) = 8*pc(16,20,3) - 12*pc(16,13,3) - 12*pc(16,15,3) - 2*pc(18,20,3) &
&+ 3*pc(18,13,3) + 3*pc(18,15,3) - 2*pc(14,20,3) + 3*pc(14,13,3) + 3*pc(14,15,3&
&)

sh(16,14,3) = 8*pc(17,20,3) - 12*pc(17,13,3) - 12*pc(17,15,3) - 2*pc(12,20,3) &
&+ 3*pc(12,13,3) + 3*pc(12,15,3) - 2*pc(19,20,3) + 3*pc(19,13,3) + 3*pc(19,15,3&
&)

sh(1,15,3) = 4*pc(1,16,3) - pc(1,18,3) - pc(1,14,3)

sh(2,15,3) = 4*pc(2,16,3) - pc(2,18,3) - pc(2,14,3)

sh(3,15,3) = 4*pc(3,16,3) - pc(3,18,3) - pc(3,14,3)

sh(4,15,3) = 4*pc(4,16,3) - pc(4,18,3) - pc(4,14,3)

sh(5,15,3) = 4*pc(8,16,3) - pc(8,18,3) - pc(8,14,3)

sh(6,15,3) = 4*pc(9,16,3) - pc(9,18,3) - pc(9,14,3)

sh(7,15,3) = 4*pc(10,16,3) - pc(10,18,3) - pc(10,14,3)

sh(8,15,3) = 4*pc(5,16,3) - pc(5,18,3) - pc(5,14,3) - 4*pc(6,16,3) + pc(6,18,3&
&) + pc(6,14,3)

sh(9,15,3) = 8*pc(7,16,3) - 2*pc(7,18,3) - 2*pc(7,14,3) - 4*pc(5,16,3) + pc(5,&
&18,3) + pc(5,14,3) - 4*pc(6,16,3) + pc(6,18,3) + pc(6,14,3)

sh(10,15,3) = 4*pc(11,16,3) - pc(11,18,3) - pc(11,14,3)

sh(11,15,3) = 4*pc(13,16,3) - pc(13,18,3) - pc(13,14,3) - 4*pc(15,16,3) + pc(1&
&5,18,3) + pc(15,14,3)

sh(12,15,3) = 4*pc(18,16,3) - pc(18,18,3) - pc(18,14,3) - 12*pc(14,16,3) + 3*p&
&c(14,18,3) + 3*pc(14,14,3)

sh(13,15,3) = 12*pc(12,16,3) - 3*pc(12,18,3) - 3*pc(12,14,3) - 4*pc(19,16,3) +&
& pc(19,18,3) + pc(19,14,3)

sh(14,15,3) = 8*pc(20,16,3) - 2*pc(20,18,3) - 2*pc(20,14,3) - 12*pc(13,16,3) +&
& 3*pc(13,18,3) + 3*pc(13,14,3) - 12*pc(15,16,3) + 3*pc(15,18,3) + 3*pc(15,14,3&
&)

sh(15,15,3) = 16*pc(16,16,3) - 4*pc(16,18,3) - 4*pc(16,14,3) - 4*pc(18,16,3) +&
& pc(18,18,3) + pc(18,14,3) - 4*pc(14,16,3) + pc(14,18,3) + pc(14,14,3)

sh(16,15,3) = 16*pc(17,16,3) - 4*pc(17,18,3) - 4*pc(17,14,3) - 4*pc(12,16,3) +&
& pc(12,18,3) + pc(12,14,3) - 4*pc(19,16,3) + pc(19,18,3) + pc(19,14,3)

sh(1,16,3) = 4*pc(1,17,3) - pc(1,12,3) - pc(1,19,3)

sh(2,16,3) = 4*pc(2,17,3) - pc(2,12,3) - pc(2,19,3)

sh(3,16,3) = 4*pc(3,17,3) - pc(3,12,3) - pc(3,19,3)

sh(4,16,3) = 4*pc(4,17,3) - pc(4,12,3) - pc(4,19,3)

sh(5,16,3) = 4*pc(8,17,3) - pc(8,12,3) - pc(8,19,3)

sh(6,16,3) = 4*pc(9,17,3) - pc(9,12,3) - pc(9,19,3)

sh(7,16,3) = 4*pc(10,17,3) - pc(10,12,3) - pc(10,19,3)

sh(8,16,3) = 4*pc(5,17,3) - pc(5,12,3) - pc(5,19,3) - 4*pc(6,17,3) + pc(6,12,3&
&) + pc(6,19,3)

sh(9,16,3) = 8*pc(7,17,3) - 2*pc(7,12,3) - 2*pc(7,19,3) - 4*pc(5,17,3) + pc(5,&
&12,3) + pc(5,19,3) - 4*pc(6,17,3) + pc(6,12,3) + pc(6,19,3)

sh(10,16,3) = 4*pc(11,17,3) - pc(11,12,3) - pc(11,19,3)

sh(11,16,3) = 4*pc(13,17,3) - pc(13,12,3) - pc(13,19,3) - 4*pc(15,17,3) + pc(1&
&5,12,3) + pc(15,19,3)

sh(12,16,3) = 4*pc(18,17,3) - pc(18,12,3) - pc(18,19,3) - 12*pc(14,17,3) + 3*p&
&c(14,12,3) + 3*pc(14,19,3)

sh(13,16,3) = 12*pc(12,17,3) - 3*pc(12,12,3) - 3*pc(12,19,3) - 4*pc(19,17,3) +&
& pc(19,12,3) + pc(19,19,3)

sh(14,16,3) = 8*pc(20,17,3) - 2*pc(20,12,3) - 2*pc(20,19,3) - 12*pc(13,17,3) +&
& 3*pc(13,12,3) + 3*pc(13,19,3) - 12*pc(15,17,3) + 3*pc(15,12,3) + 3*pc(15,19,3&
&)

sh(15,16,3) = 16*pc(16,17,3) - 4*pc(16,12,3) - 4*pc(16,19,3) - 4*pc(18,17,3) +&
& pc(18,12,3) + pc(18,19,3) - 4*pc(14,17,3) + pc(14,12,3) + pc(14,19,3)

sh(16,16,3) = 16*pc(17,17,3) - 4*pc(17,12,3) - 4*pc(17,19,3) - 4*pc(12,17,3) +&
& pc(12,12,3) + pc(12,19,3) - 4*pc(19,17,3) + pc(19,12,3) + pc(19,19,3)


   end subroutine momentum2CIntgAna

   subroutine momentum2CIntgNum(a1,a2,A,B,pc,sh,cell_size,step_size)

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
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   real (kind=double), dimension (20,20,3), intent(out) :: pc
   real (kind=double), dimension (16,16,3), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i, j
   integer :: num_steps
   integer, dimension (20,3) :: triads
   integer, dimension (16,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, curr_pos
   real (kind=double), dimension (3) :: xyz
   real (kind=double), dimension (3,2) :: xyz_I ! Indices=xyz, noprime||prime

   ! Before we proceed with the calculation we need to understand a bit more
   !   about exactly what is being computed. The form of the integration is:
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        d/d_xyz  [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   ! We use Pxyz for points in space; Rxyz1,2 for atomic sites 1 and 2;
   !   zeta1,2 for the exponential decay factors; lxyz1,2 for the angular
   !   momentum numbers for atoms 1 and 2, and d/d_xyz for independent
   !   derivatives (i.e., we have three separate triple integrals).
   ! The integral must be computed over all space for all different possible
   !   values of lxyz1,2 for some arbitrarily chosen test coordinates and
   !   decay rates.
   ! Now, focusing on just the d/dx version of the set of triple integrals,
   !   we will pull all terms with Py and Pz out of the dx integral to get:
   ! SS { [(Py-Ry1)**ly1 * (Pz-Rz1)**lz1 *
   !       exp(-zeta1*((Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !      [(Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !       exp(-zeta2*((Py-Ry2)**2 + (Pz-Rz2)**2 ))] *
   !     S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * d/dx [
   !      (Px-Rx2)**lx2 * exp(-zeta2*((Px-Rx2)**2))]] dx
   !    } dydz
   ! Applying the derivative, the internal 1D dx integral has the form:
   !   Ix' = S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) *
   !           [(Px-Rx2)**lx2 * -zeta2*(Px-Rx2) + (Px-Rx2)**(lx2-1)]
   !           * exp(-zeta2*(Px-Rx2)**2)
   !         ]
   ! Each of the other integrals (Iy, Iz) will have the form of a simple
   !   1D overlap integral (without a derivative):
   !   Iy = S [(Py-Ry1)**ly1 * exp(-zeta1*((Py-Ry1)**2)) *
   !           (Py-Ry2)**ly2 * exp(-zeta2*((Py-Ry2)**2))] dy
   !   Iz = S [(Pz-Rz1)**lz1 * exp(-zeta1*((Pz-Rz1)**2)) *
   !           (Pz-Rz2)**lz2 * exp(-zeta2*((Pz-Rz2)**2))] dz
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
   !   momentum we must discard this term when lx2-1 < 0.

   ! Initialize local variables.
   triads(1,:) = (/0,0,0/)
   triads(2,:) = (/1,0,0/)
   triads(3,:) = (/0,1,0/)
   triads(4,:) = (/0,0,1/)
   triads(5,:) = (/2,0,0/)
   triads(6,:) = (/0,2,0/)
   triads(7,:) = (/0,0,2/)
   triads(8,:) = (/1,1,0/)
   triads(9,:) = (/1,0,1/)
   triads(10,:) = (/0,1,1/)
   triads(11,:) = (/1,1,1/)
   triads(12,:) = (/2,1,0/)
   triads(13,:) = (/2,0,1/)
   triads(14,:) = (/1,2,0/)
   triads(15,:) = (/0,2,1/)
   triads(16,:) = (/1,0,2/)
   triads(17,:) = (/0,1,2/)
   triads(18,:) = (/3,0,0/)
   triads(19,:) = (/0,3,0/)
   triads(20,:) = (/0,0,3/)
   conversion(1,1,:) = (/1,0,0/)
   conversion(1,2,:) = (/1,1,1/)
   conversion(2,1,:) = (/1,0,0/)
   conversion(2,2,:) = (/2,1,1/)
   conversion(3,1,:) = (/1,0,0/)
   conversion(3,2,:) = (/3,1,1/)
   conversion(4,1,:) = (/1,0,0/)
   conversion(4,2,:) = (/4,1,1/)
   conversion(5,1,:) = (/1,0,0/)
   conversion(5,2,:) = (/8,1,1/)
   conversion(6,1,:) = (/1,0,0/)
   conversion(6,2,:) = (/9,1,1/)
   conversion(7,1,:) = (/1,0,0/)
   conversion(7,2,:) = (/10,1,1/)
   conversion(8,1,:) = (/1,-1,0/)
   conversion(8,2,:) = (/5,6,1/)
   conversion(9,1,:) = (/2,-1,-1/)
   conversion(9,2,:) = (/7,5,6/)
   conversion(10,1,:) = (/1,0,0/)
   conversion(10,2,:) = (/11,1,1/)
   conversion(11,1,:) = (/1,-1,0/)
   conversion(11,2,:) = (/13,15,1/)
   conversion(12,1,:) = (/1,-3,0/)
   conversion(12,2,:) = (/18,14,1/)
   conversion(13,1,:) = (/3,-1,0/)
   conversion(13,2,:) = (/12,19,1/)
   conversion(14,1,:) = (/2,-3,-3/)
   conversion(14,2,:) = (/20,13,15/)
   conversion(15,1,:) = (/4,-1,-1/)
   conversion(15,2,:) = (/16,18,14/)
   conversion(16,1,:) = (/4,-1,-1/)
   conversion(16,2,:) = (/17,12,19/)

   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   do p = 1, 20
      do q = 1, 20

         ! Assign l1 and l2 values for each gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize sum variable.
         xyz_I(:,:) = 0.0d0

         ! Start a loop over x coordinates.
         do i = 0, num_steps
            curr_pos = start_pos + (i*step_size)

            ! Compute the no-prime integrals first.
            do j = 1, 3
               xyz_I(j,2) = xyz_I(j,2) &
                     & + noPrimeMM(step_size, curr_pos, A(j), &
                     & B(j), a1, a2, l1(j), l2(j))
            enddo

            ! Compute the prime integrals second.
            do j = 1, 3
               xyz_I(j,1) = xyz_I(j,1) &
                     & + primeMM(step_size, curr_pos, A(j), &
                     & B(j), a1, a2, l1(j), l2(j))
            enddo
         enddo

         pc(q,p,1) = xyz_I(1,1)*xyz_I(2,2)*xyz_I(3,2)
         pc(q,p,2) = xyz_I(1,2)*xyz_I(2,1)*xyz_I(3,2)
         pc(q,p,3) = xyz_I(1,2)*xyz_I(2,2)*xyz_I(3,1)
      enddo
    enddo

sh(1,1,1) = pc(1,1,1)

sh(2,1,1) = pc(2,1,1)

sh(3,1,1) = pc(3,1,1)

sh(4,1,1) = pc(4,1,1)

sh(5,1,1) = pc(8,1,1)

sh(6,1,1) = pc(9,1,1)

sh(7,1,1) = pc(10,1,1)

sh(8,1,1) = pc(5,1,1) - pc(6,1,1)

sh(9,1,1) = 2*pc(7,1,1) - pc(5,1,1) - pc(6,1,1)

sh(10,1,1) = pc(11,1,1)

sh(11,1,1) = pc(13,1,1) - pc(15,1,1)

sh(12,1,1) = pc(18,1,1) - 3*pc(14,1,1)

sh(13,1,1) = 3*pc(12,1,1) - pc(19,1,1)

sh(14,1,1) = 2*pc(20,1,1) - 3*pc(13,1,1) - 3*pc(15,1,1)

sh(15,1,1) = 4*pc(16,1,1) - pc(18,1,1) - pc(14,1,1)

sh(16,1,1) = 4*pc(17,1,1) - pc(12,1,1) - pc(19,1,1)

sh(1,2,1) = pc(1,2,1)

sh(2,2,1) = pc(2,2,1)

sh(3,2,1) = pc(3,2,1)

sh(4,2,1) = pc(4,2,1)

sh(5,2,1) = pc(8,2,1)

sh(6,2,1) = pc(9,2,1)

sh(7,2,1) = pc(10,2,1)

sh(8,2,1) = pc(5,2,1) - pc(6,2,1)

sh(9,2,1) = 2*pc(7,2,1) - pc(5,2,1) - pc(6,2,1)

sh(10,2,1) = pc(11,2,1)

sh(11,2,1) = pc(13,2,1) - pc(15,2,1)

sh(12,2,1) = pc(18,2,1) - 3*pc(14,2,1)

sh(13,2,1) = 3*pc(12,2,1) - pc(19,2,1)

sh(14,2,1) = 2*pc(20,2,1) - 3*pc(13,2,1) - 3*pc(15,2,1)

sh(15,2,1) = 4*pc(16,2,1) - pc(18,2,1) - pc(14,2,1)

sh(16,2,1) = 4*pc(17,2,1) - pc(12,2,1) - pc(19,2,1)

sh(1,3,1) = pc(1,3,1)

sh(2,3,1) = pc(2,3,1)

sh(3,3,1) = pc(3,3,1)

sh(4,3,1) = pc(4,3,1)

sh(5,3,1) = pc(8,3,1)

sh(6,3,1) = pc(9,3,1)

sh(7,3,1) = pc(10,3,1)

sh(8,3,1) = pc(5,3,1) - pc(6,3,1)

sh(9,3,1) = 2*pc(7,3,1) - pc(5,3,1) - pc(6,3,1)

sh(10,3,1) = pc(11,3,1)

sh(11,3,1) = pc(13,3,1) - pc(15,3,1)

sh(12,3,1) = pc(18,3,1) - 3*pc(14,3,1)

sh(13,3,1) = 3*pc(12,3,1) - pc(19,3,1)

sh(14,3,1) = 2*pc(20,3,1) - 3*pc(13,3,1) - 3*pc(15,3,1)

sh(15,3,1) = 4*pc(16,3,1) - pc(18,3,1) - pc(14,3,1)

sh(16,3,1) = 4*pc(17,3,1) - pc(12,3,1) - pc(19,3,1)

sh(1,4,1) = pc(1,4,1)

sh(2,4,1) = pc(2,4,1)

sh(3,4,1) = pc(3,4,1)

sh(4,4,1) = pc(4,4,1)

sh(5,4,1) = pc(8,4,1)

sh(6,4,1) = pc(9,4,1)

sh(7,4,1) = pc(10,4,1)

sh(8,4,1) = pc(5,4,1) - pc(6,4,1)

sh(9,4,1) = 2*pc(7,4,1) - pc(5,4,1) - pc(6,4,1)

sh(10,4,1) = pc(11,4,1)

sh(11,4,1) = pc(13,4,1) - pc(15,4,1)

sh(12,4,1) = pc(18,4,1) - 3*pc(14,4,1)

sh(13,4,1) = 3*pc(12,4,1) - pc(19,4,1)

sh(14,4,1) = 2*pc(20,4,1) - 3*pc(13,4,1) - 3*pc(15,4,1)

sh(15,4,1) = 4*pc(16,4,1) - pc(18,4,1) - pc(14,4,1)

sh(16,4,1) = 4*pc(17,4,1) - pc(12,4,1) - pc(19,4,1)

sh(1,5,1) = pc(1,8,1)

sh(2,5,1) = pc(2,8,1)

sh(3,5,1) = pc(3,8,1)

sh(4,5,1) = pc(4,8,1)

sh(5,5,1) = pc(8,8,1)

sh(6,5,1) = pc(9,8,1)

sh(7,5,1) = pc(10,8,1)

sh(8,5,1) = pc(5,8,1) - pc(6,8,1)

sh(9,5,1) = 2*pc(7,8,1) - pc(5,8,1) - pc(6,8,1)

sh(10,5,1) = pc(11,8,1)

sh(11,5,1) = pc(13,8,1) - pc(15,8,1)

sh(12,5,1) = pc(18,8,1) - 3*pc(14,8,1)

sh(13,5,1) = 3*pc(12,8,1) - pc(19,8,1)

sh(14,5,1) = 2*pc(20,8,1) - 3*pc(13,8,1) - 3*pc(15,8,1)

sh(15,5,1) = 4*pc(16,8,1) - pc(18,8,1) - pc(14,8,1)

sh(16,5,1) = 4*pc(17,8,1) - pc(12,8,1) - pc(19,8,1)

sh(1,6,1) = pc(1,9,1)

sh(2,6,1) = pc(2,9,1)

sh(3,6,1) = pc(3,9,1)

sh(4,6,1) = pc(4,9,1)

sh(5,6,1) = pc(8,9,1)

sh(6,6,1) = pc(9,9,1)

sh(7,6,1) = pc(10,9,1)

sh(8,6,1) = pc(5,9,1) - pc(6,9,1)

sh(9,6,1) = 2*pc(7,9,1) - pc(5,9,1) - pc(6,9,1)

sh(10,6,1) = pc(11,9,1)

sh(11,6,1) = pc(13,9,1) - pc(15,9,1)

sh(12,6,1) = pc(18,9,1) - 3*pc(14,9,1)

sh(13,6,1) = 3*pc(12,9,1) - pc(19,9,1)

sh(14,6,1) = 2*pc(20,9,1) - 3*pc(13,9,1) - 3*pc(15,9,1)

sh(15,6,1) = 4*pc(16,9,1) - pc(18,9,1) - pc(14,9,1)

sh(16,6,1) = 4*pc(17,9,1) - pc(12,9,1) - pc(19,9,1)

sh(1,7,1) = pc(1,10,1)

sh(2,7,1) = pc(2,10,1)

sh(3,7,1) = pc(3,10,1)

sh(4,7,1) = pc(4,10,1)

sh(5,7,1) = pc(8,10,1)

sh(6,7,1) = pc(9,10,1)

sh(7,7,1) = pc(10,10,1)

sh(8,7,1) = pc(5,10,1) - pc(6,10,1)

sh(9,7,1) = 2*pc(7,10,1) - pc(5,10,1) - pc(6,10,1)

sh(10,7,1) = pc(11,10,1)

sh(11,7,1) = pc(13,10,1) - pc(15,10,1)

sh(12,7,1) = pc(18,10,1) - 3*pc(14,10,1)

sh(13,7,1) = 3*pc(12,10,1) - pc(19,10,1)

sh(14,7,1) = 2*pc(20,10,1) - 3*pc(13,10,1) - 3*pc(15,10,1)

sh(15,7,1) = 4*pc(16,10,1) - pc(18,10,1) - pc(14,10,1)

sh(16,7,1) = 4*pc(17,10,1) - pc(12,10,1) - pc(19,10,1)

sh(1,8,1) = pc(1,5,1) - pc(1,6,1)

sh(2,8,1) = pc(2,5,1) - pc(2,6,1)

sh(3,8,1) = pc(3,5,1) - pc(3,6,1)

sh(4,8,1) = pc(4,5,1) - pc(4,6,1)

sh(5,8,1) = pc(8,5,1) - pc(8,6,1)

sh(6,8,1) = pc(9,5,1) - pc(9,6,1)

sh(7,8,1) = pc(10,5,1) - pc(10,6,1)

sh(8,8,1) = pc(5,5,1) - pc(5,6,1) - pc(6,5,1) + pc(6,6,1)

sh(9,8,1) = 2*pc(7,5,1) - 2*pc(7,6,1) - pc(5,5,1) + pc(5,6,1) - pc(6,5,1) + pc&
&(6,6,1)

sh(10,8,1) = pc(11,5,1) - pc(11,6,1)

sh(11,8,1) = pc(13,5,1) - pc(13,6,1) - pc(15,5,1) + pc(15,6,1)

sh(12,8,1) = pc(18,5,1) - pc(18,6,1) - 3*pc(14,5,1) + 3*pc(14,6,1)

sh(13,8,1) = 3*pc(12,5,1) - 3*pc(12,6,1) - pc(19,5,1) + pc(19,6,1)

sh(14,8,1) = 2*pc(20,5,1) - 2*pc(20,6,1) - 3*pc(13,5,1) + 3*pc(13,6,1) - 3*pc(&
&15,5,1) + 3*pc(15,6,1)

sh(15,8,1) = 4*pc(16,5,1) - 4*pc(16,6,1) - pc(18,5,1) + pc(18,6,1) - pc(14,5,1&
&) + pc(14,6,1)

sh(16,8,1) = 4*pc(17,5,1) - 4*pc(17,6,1) - pc(12,5,1) + pc(12,6,1) - pc(19,5,1&
&) + pc(19,6,1)

sh(1,9,1) = 2*pc(1,7,1) - pc(1,5,1) - pc(1,6,1)

sh(2,9,1) = 2*pc(2,7,1) - pc(2,5,1) - pc(2,6,1)

sh(3,9,1) = 2*pc(3,7,1) - pc(3,5,1) - pc(3,6,1)

sh(4,9,1) = 2*pc(4,7,1) - pc(4,5,1) - pc(4,6,1)

sh(5,9,1) = 2*pc(8,7,1) - pc(8,5,1) - pc(8,6,1)

sh(6,9,1) = 2*pc(9,7,1) - pc(9,5,1) - pc(9,6,1)

sh(7,9,1) = 2*pc(10,7,1) - pc(10,5,1) - pc(10,6,1)

sh(8,9,1) = 2*pc(5,7,1) - pc(5,5,1) - pc(5,6,1) - 2*pc(6,7,1) + pc(6,5,1) + pc&
&(6,6,1)

sh(9,9,1) = 4*pc(7,7,1) - 2*pc(7,5,1) - 2*pc(7,6,1) - 2*pc(5,7,1) + pc(5,5,1) &
&+ pc(5,6,1) - 2*pc(6,7,1) + pc(6,5,1) + pc(6,6,1)

sh(10,9,1) = 2*pc(11,7,1) - pc(11,5,1) - pc(11,6,1)

sh(11,9,1) = 2*pc(13,7,1) - pc(13,5,1) - pc(13,6,1) - 2*pc(15,7,1) + pc(15,5,1&
&) + pc(15,6,1)

sh(12,9,1) = 2*pc(18,7,1) - pc(18,5,1) - pc(18,6,1) - 6*pc(14,7,1) + 3*pc(14,5&
&,1) + 3*pc(14,6,1)

sh(13,9,1) = 6*pc(12,7,1) - 3*pc(12,5,1) - 3*pc(12,6,1) - 2*pc(19,7,1) + pc(19&
&,5,1) + pc(19,6,1)

sh(14,9,1) = 4*pc(20,7,1) - 2*pc(20,5,1) - 2*pc(20,6,1) - 6*pc(13,7,1) + 3*pc(&
&13,5,1) + 3*pc(13,6,1) - 6*pc(15,7,1) + 3*pc(15,5,1) + 3*pc(15,6,1)

sh(15,9,1) = 8*pc(16,7,1) - 4*pc(16,5,1) - 4*pc(16,6,1) - 2*pc(18,7,1) + pc(18&
&,5,1) + pc(18,6,1) - 2*pc(14,7,1) + pc(14,5,1) + pc(14,6,1)

sh(16,9,1) = 8*pc(17,7,1) - 4*pc(17,5,1) - 4*pc(17,6,1) - 2*pc(12,7,1) + pc(12&
&,5,1) + pc(12,6,1) - 2*pc(19,7,1) + pc(19,5,1) + pc(19,6,1)

sh(1,10,1) = pc(1,11,1)

sh(2,10,1) = pc(2,11,1)

sh(3,10,1) = pc(3,11,1)

sh(4,10,1) = pc(4,11,1)

sh(5,10,1) = pc(8,11,1)

sh(6,10,1) = pc(9,11,1)

sh(7,10,1) = pc(10,11,1)

sh(8,10,1) = pc(5,11,1) - pc(6,11,1)

sh(9,10,1) = 2*pc(7,11,1) - pc(5,11,1) - pc(6,11,1)

sh(10,10,1) = pc(11,11,1)

sh(11,10,1) = pc(13,11,1) - pc(15,11,1)

sh(12,10,1) = pc(18,11,1) - 3*pc(14,11,1)

sh(13,10,1) = 3*pc(12,11,1) - pc(19,11,1)

sh(14,10,1) = 2*pc(20,11,1) - 3*pc(13,11,1) - 3*pc(15,11,1)

sh(15,10,1) = 4*pc(16,11,1) - pc(18,11,1) - pc(14,11,1)

sh(16,10,1) = 4*pc(17,11,1) - pc(12,11,1) - pc(19,11,1)

sh(1,11,1) = pc(1,13,1) - pc(1,15,1)

sh(2,11,1) = pc(2,13,1) - pc(2,15,1)

sh(3,11,1) = pc(3,13,1) - pc(3,15,1)

sh(4,11,1) = pc(4,13,1) - pc(4,15,1)

sh(5,11,1) = pc(8,13,1) - pc(8,15,1)

sh(6,11,1) = pc(9,13,1) - pc(9,15,1)

sh(7,11,1) = pc(10,13,1) - pc(10,15,1)

sh(8,11,1) = pc(5,13,1) - pc(5,15,1) - pc(6,13,1) + pc(6,15,1)

sh(9,11,1) = 2*pc(7,13,1) - 2*pc(7,15,1) - pc(5,13,1) + pc(5,15,1) - pc(6,13,1&
&) + pc(6,15,1)

sh(10,11,1) = pc(11,13,1) - pc(11,15,1)

sh(11,11,1) = pc(13,13,1) - pc(13,15,1) - pc(15,13,1) + pc(15,15,1)

sh(12,11,1) = pc(18,13,1) - pc(18,15,1) - 3*pc(14,13,1) + 3*pc(14,15,1)

sh(13,11,1) = 3*pc(12,13,1) - 3*pc(12,15,1) - pc(19,13,1) + pc(19,15,1)

sh(14,11,1) = 2*pc(20,13,1) - 2*pc(20,15,1) - 3*pc(13,13,1) + 3*pc(13,15,1) - &
&3*pc(15,13,1) + 3*pc(15,15,1)

sh(15,11,1) = 4*pc(16,13,1) - 4*pc(16,15,1) - pc(18,13,1) + pc(18,15,1) - pc(1&
&4,13,1) + pc(14,15,1)

sh(16,11,1) = 4*pc(17,13,1) - 4*pc(17,15,1) - pc(12,13,1) + pc(12,15,1) - pc(1&
&9,13,1) + pc(19,15,1)

sh(1,12,1) = pc(1,18,1) - 3*pc(1,14,1)

sh(2,12,1) = pc(2,18,1) - 3*pc(2,14,1)

sh(3,12,1) = pc(3,18,1) - 3*pc(3,14,1)

sh(4,12,1) = pc(4,18,1) - 3*pc(4,14,1)

sh(5,12,1) = pc(8,18,1) - 3*pc(8,14,1)

sh(6,12,1) = pc(9,18,1) - 3*pc(9,14,1)

sh(7,12,1) = pc(10,18,1) - 3*pc(10,14,1)

sh(8,12,1) = pc(5,18,1) - 3*pc(5,14,1) - pc(6,18,1) + 3*pc(6,14,1)

sh(9,12,1) = 2*pc(7,18,1) - 6*pc(7,14,1) - pc(5,18,1) + 3*pc(5,14,1) - pc(6,18&
&,1) + 3*pc(6,14,1)

sh(10,12,1) = pc(11,18,1) - 3*pc(11,14,1)

sh(11,12,1) = pc(13,18,1) - 3*pc(13,14,1) - pc(15,18,1) + 3*pc(15,14,1)

sh(12,12,1) = pc(18,18,1) - 3*pc(18,14,1) - 3*pc(14,18,1) + 9*pc(14,14,1)

sh(13,12,1) = 3*pc(12,18,1) - 9*pc(12,14,1) - pc(19,18,1) + 3*pc(19,14,1)

sh(14,12,1) = 2*pc(20,18,1) - 6*pc(20,14,1) - 3*pc(13,18,1) + 9*pc(13,14,1) - &
&3*pc(15,18,1) + 9*pc(15,14,1)

sh(15,12,1) = 4*pc(16,18,1) - 12*pc(16,14,1) - pc(18,18,1) + 3*pc(18,14,1) - p&
&c(14,18,1) + 3*pc(14,14,1)

sh(16,12,1) = 4*pc(17,18,1) - 12*pc(17,14,1) - pc(12,18,1) + 3*pc(12,14,1) - p&
&c(19,18,1) + 3*pc(19,14,1)

sh(1,13,1) = 3*pc(1,12,1) - pc(1,19,1)

sh(2,13,1) = 3*pc(2,12,1) - pc(2,19,1)

sh(3,13,1) = 3*pc(3,12,1) - pc(3,19,1)

sh(4,13,1) = 3*pc(4,12,1) - pc(4,19,1)

sh(5,13,1) = 3*pc(8,12,1) - pc(8,19,1)

sh(6,13,1) = 3*pc(9,12,1) - pc(9,19,1)

sh(7,13,1) = 3*pc(10,12,1) - pc(10,19,1)

sh(8,13,1) = 3*pc(5,12,1) - pc(5,19,1) - 3*pc(6,12,1) + pc(6,19,1)

sh(9,13,1) = 6*pc(7,12,1) - 2*pc(7,19,1) - 3*pc(5,12,1) + pc(5,19,1) - 3*pc(6,&
&12,1) + pc(6,19,1)

sh(10,13,1) = 3*pc(11,12,1) - pc(11,19,1)

sh(11,13,1) = 3*pc(13,12,1) - pc(13,19,1) - 3*pc(15,12,1) + pc(15,19,1)

sh(12,13,1) = 3*pc(18,12,1) - pc(18,19,1) - 9*pc(14,12,1) + 3*pc(14,19,1)

sh(13,13,1) = 9*pc(12,12,1) - 3*pc(12,19,1) - 3*pc(19,12,1) + pc(19,19,1)

sh(14,13,1) = 6*pc(20,12,1) - 2*pc(20,19,1) - 9*pc(13,12,1) + 3*pc(13,19,1) - &
&9*pc(15,12,1) + 3*pc(15,19,1)

sh(15,13,1) = 12*pc(16,12,1) - 4*pc(16,19,1) - 3*pc(18,12,1) + pc(18,19,1) - 3&
&*pc(14,12,1) + pc(14,19,1)

sh(16,13,1) = 12*pc(17,12,1) - 4*pc(17,19,1) - 3*pc(12,12,1) + pc(12,19,1) - 3&
&*pc(19,12,1) + pc(19,19,1)

sh(1,14,1) = 2*pc(1,20,1) - 3*pc(1,13,1) - 3*pc(1,15,1)

sh(2,14,1) = 2*pc(2,20,1) - 3*pc(2,13,1) - 3*pc(2,15,1)

sh(3,14,1) = 2*pc(3,20,1) - 3*pc(3,13,1) - 3*pc(3,15,1)

sh(4,14,1) = 2*pc(4,20,1) - 3*pc(4,13,1) - 3*pc(4,15,1)

sh(5,14,1) = 2*pc(8,20,1) - 3*pc(8,13,1) - 3*pc(8,15,1)

sh(6,14,1) = 2*pc(9,20,1) - 3*pc(9,13,1) - 3*pc(9,15,1)

sh(7,14,1) = 2*pc(10,20,1) - 3*pc(10,13,1) - 3*pc(10,15,1)

sh(8,14,1) = 2*pc(5,20,1) - 3*pc(5,13,1) - 3*pc(5,15,1) - 2*pc(6,20,1) + 3*pc(&
&6,13,1) + 3*pc(6,15,1)

sh(9,14,1) = 4*pc(7,20,1) - 6*pc(7,13,1) - 6*pc(7,15,1) - 2*pc(5,20,1) + 3*pc(&
&5,13,1) + 3*pc(5,15,1) - 2*pc(6,20,1) + 3*pc(6,13,1) + 3*pc(6,15,1)

sh(10,14,1) = 2*pc(11,20,1) - 3*pc(11,13,1) - 3*pc(11,15,1)

sh(11,14,1) = 2*pc(13,20,1) - 3*pc(13,13,1) - 3*pc(13,15,1) - 2*pc(15,20,1) + &
&3*pc(15,13,1) + 3*pc(15,15,1)

sh(12,14,1) = 2*pc(18,20,1) - 3*pc(18,13,1) - 3*pc(18,15,1) - 6*pc(14,20,1) + &
&9*pc(14,13,1) + 9*pc(14,15,1)

sh(13,14,1) = 6*pc(12,20,1) - 9*pc(12,13,1) - 9*pc(12,15,1) - 2*pc(19,20,1) + &
&3*pc(19,13,1) + 3*pc(19,15,1)

sh(14,14,1) = 4*pc(20,20,1) - 6*pc(20,13,1) - 6*pc(20,15,1) - 6*pc(13,20,1) + &
&9*pc(13,13,1) + 9*pc(13,15,1) - 6*pc(15,20,1) + 9*pc(15,13,1) + 9*pc(15,15,1)

sh(15,14,1) = 8*pc(16,20,1) - 12*pc(16,13,1) - 12*pc(16,15,1) - 2*pc(18,20,1) &
&+ 3*pc(18,13,1) + 3*pc(18,15,1) - 2*pc(14,20,1) + 3*pc(14,13,1) + 3*pc(14,15,1&
&)

sh(16,14,1) = 8*pc(17,20,1) - 12*pc(17,13,1) - 12*pc(17,15,1) - 2*pc(12,20,1) &
&+ 3*pc(12,13,1) + 3*pc(12,15,1) - 2*pc(19,20,1) + 3*pc(19,13,1) + 3*pc(19,15,1&
&)

sh(1,15,1) = 4*pc(1,16,1) - pc(1,18,1) - pc(1,14,1)

sh(2,15,1) = 4*pc(2,16,1) - pc(2,18,1) - pc(2,14,1)

sh(3,15,1) = 4*pc(3,16,1) - pc(3,18,1) - pc(3,14,1)

sh(4,15,1) = 4*pc(4,16,1) - pc(4,18,1) - pc(4,14,1)

sh(5,15,1) = 4*pc(8,16,1) - pc(8,18,1) - pc(8,14,1)

sh(6,15,1) = 4*pc(9,16,1) - pc(9,18,1) - pc(9,14,1)

sh(7,15,1) = 4*pc(10,16,1) - pc(10,18,1) - pc(10,14,1)

sh(8,15,1) = 4*pc(5,16,1) - pc(5,18,1) - pc(5,14,1) - 4*pc(6,16,1) + pc(6,18,1&
&) + pc(6,14,1)

sh(9,15,1) = 8*pc(7,16,1) - 2*pc(7,18,1) - 2*pc(7,14,1) - 4*pc(5,16,1) + pc(5,&
&18,1) + pc(5,14,1) - 4*pc(6,16,1) + pc(6,18,1) + pc(6,14,1)

sh(10,15,1) = 4*pc(11,16,1) - pc(11,18,1) - pc(11,14,1)

sh(11,15,1) = 4*pc(13,16,1) - pc(13,18,1) - pc(13,14,1) - 4*pc(15,16,1) + pc(1&
&5,18,1) + pc(15,14,1)

sh(12,15,1) = 4*pc(18,16,1) - pc(18,18,1) - pc(18,14,1) - 12*pc(14,16,1) + 3*p&
&c(14,18,1) + 3*pc(14,14,1)

sh(13,15,1) = 12*pc(12,16,1) - 3*pc(12,18,1) - 3*pc(12,14,1) - 4*pc(19,16,1) +&
& pc(19,18,1) + pc(19,14,1)

sh(14,15,1) = 8*pc(20,16,1) - 2*pc(20,18,1) - 2*pc(20,14,1) - 12*pc(13,16,1) +&
& 3*pc(13,18,1) + 3*pc(13,14,1) - 12*pc(15,16,1) + 3*pc(15,18,1) + 3*pc(15,14,1&
&)

sh(15,15,1) = 16*pc(16,16,1) - 4*pc(16,18,1) - 4*pc(16,14,1) - 4*pc(18,16,1) +&
& pc(18,18,1) + pc(18,14,1) - 4*pc(14,16,1) + pc(14,18,1) + pc(14,14,1)

sh(16,15,1) = 16*pc(17,16,1) - 4*pc(17,18,1) - 4*pc(17,14,1) - 4*pc(12,16,1) +&
& pc(12,18,1) + pc(12,14,1) - 4*pc(19,16,1) + pc(19,18,1) + pc(19,14,1)

sh(1,16,1) = 4*pc(1,17,1) - pc(1,12,1) - pc(1,19,1)

sh(2,16,1) = 4*pc(2,17,1) - pc(2,12,1) - pc(2,19,1)

sh(3,16,1) = 4*pc(3,17,1) - pc(3,12,1) - pc(3,19,1)

sh(4,16,1) = 4*pc(4,17,1) - pc(4,12,1) - pc(4,19,1)

sh(5,16,1) = 4*pc(8,17,1) - pc(8,12,1) - pc(8,19,1)

sh(6,16,1) = 4*pc(9,17,1) - pc(9,12,1) - pc(9,19,1)

sh(7,16,1) = 4*pc(10,17,1) - pc(10,12,1) - pc(10,19,1)

sh(8,16,1) = 4*pc(5,17,1) - pc(5,12,1) - pc(5,19,1) - 4*pc(6,17,1) + pc(6,12,1&
&) + pc(6,19,1)

sh(9,16,1) = 8*pc(7,17,1) - 2*pc(7,12,1) - 2*pc(7,19,1) - 4*pc(5,17,1) + pc(5,&
&12,1) + pc(5,19,1) - 4*pc(6,17,1) + pc(6,12,1) + pc(6,19,1)

sh(10,16,1) = 4*pc(11,17,1) - pc(11,12,1) - pc(11,19,1)

sh(11,16,1) = 4*pc(13,17,1) - pc(13,12,1) - pc(13,19,1) - 4*pc(15,17,1) + pc(1&
&5,12,1) + pc(15,19,1)

sh(12,16,1) = 4*pc(18,17,1) - pc(18,12,1) - pc(18,19,1) - 12*pc(14,17,1) + 3*p&
&c(14,12,1) + 3*pc(14,19,1)

sh(13,16,1) = 12*pc(12,17,1) - 3*pc(12,12,1) - 3*pc(12,19,1) - 4*pc(19,17,1) +&
& pc(19,12,1) + pc(19,19,1)

sh(14,16,1) = 8*pc(20,17,1) - 2*pc(20,12,1) - 2*pc(20,19,1) - 12*pc(13,17,1) +&
& 3*pc(13,12,1) + 3*pc(13,19,1) - 12*pc(15,17,1) + 3*pc(15,12,1) + 3*pc(15,19,1&
&)

sh(15,16,1) = 16*pc(16,17,1) - 4*pc(16,12,1) - 4*pc(16,19,1) - 4*pc(18,17,1) +&
& pc(18,12,1) + pc(18,19,1) - 4*pc(14,17,1) + pc(14,12,1) + pc(14,19,1)

sh(16,16,1) = 16*pc(17,17,1) - 4*pc(17,12,1) - 4*pc(17,19,1) - 4*pc(12,17,1) +&
& pc(12,12,1) + pc(12,19,1) - 4*pc(19,17,1) + pc(19,12,1) + pc(19,19,1)

sh(1,1,2) = pc(1,1,2)

sh(2,1,2) = pc(2,1,2)

sh(3,1,2) = pc(3,1,2)

sh(4,1,2) = pc(4,1,2)

sh(5,1,2) = pc(8,1,2)

sh(6,1,2) = pc(9,1,2)

sh(7,1,2) = pc(10,1,2)

sh(8,1,2) = pc(5,1,2) - pc(6,1,2)

sh(9,1,2) = 2*pc(7,1,2) - pc(5,1,2) - pc(6,1,2)

sh(10,1,2) = pc(11,1,2)

sh(11,1,2) = pc(13,1,2) - pc(15,1,2)

sh(12,1,2) = pc(18,1,2) - 3*pc(14,1,2)

sh(13,1,2) = 3*pc(12,1,2) - pc(19,1,2)

sh(14,1,2) = 2*pc(20,1,2) - 3*pc(13,1,2) - 3*pc(15,1,2)

sh(15,1,2) = 4*pc(16,1,2) - pc(18,1,2) - pc(14,1,2)

sh(16,1,2) = 4*pc(17,1,2) - pc(12,1,2) - pc(19,1,2)

sh(1,2,2) = pc(1,2,2)

sh(2,2,2) = pc(2,2,2)

sh(3,2,2) = pc(3,2,2)

sh(4,2,2) = pc(4,2,2)

sh(5,2,2) = pc(8,2,2)

sh(6,2,2) = pc(9,2,2)

sh(7,2,2) = pc(10,2,2)

sh(8,2,2) = pc(5,2,2) - pc(6,2,2)

sh(9,2,2) = 2*pc(7,2,2) - pc(5,2,2) - pc(6,2,2)

sh(10,2,2) = pc(11,2,2)

sh(11,2,2) = pc(13,2,2) - pc(15,2,2)

sh(12,2,2) = pc(18,2,2) - 3*pc(14,2,2)

sh(13,2,2) = 3*pc(12,2,2) - pc(19,2,2)

sh(14,2,2) = 2*pc(20,2,2) - 3*pc(13,2,2) - 3*pc(15,2,2)

sh(15,2,2) = 4*pc(16,2,2) - pc(18,2,2) - pc(14,2,2)

sh(16,2,2) = 4*pc(17,2,2) - pc(12,2,2) - pc(19,2,2)

sh(1,3,2) = pc(1,3,2)

sh(2,3,2) = pc(2,3,2)

sh(3,3,2) = pc(3,3,2)

sh(4,3,2) = pc(4,3,2)

sh(5,3,2) = pc(8,3,2)

sh(6,3,2) = pc(9,3,2)

sh(7,3,2) = pc(10,3,2)

sh(8,3,2) = pc(5,3,2) - pc(6,3,2)

sh(9,3,2) = 2*pc(7,3,2) - pc(5,3,2) - pc(6,3,2)

sh(10,3,2) = pc(11,3,2)

sh(11,3,2) = pc(13,3,2) - pc(15,3,2)

sh(12,3,2) = pc(18,3,2) - 3*pc(14,3,2)

sh(13,3,2) = 3*pc(12,3,2) - pc(19,3,2)

sh(14,3,2) = 2*pc(20,3,2) - 3*pc(13,3,2) - 3*pc(15,3,2)

sh(15,3,2) = 4*pc(16,3,2) - pc(18,3,2) - pc(14,3,2)

sh(16,3,2) = 4*pc(17,3,2) - pc(12,3,2) - pc(19,3,2)

sh(1,4,2) = pc(1,4,2)

sh(2,4,2) = pc(2,4,2)

sh(3,4,2) = pc(3,4,2)

sh(4,4,2) = pc(4,4,2)

sh(5,4,2) = pc(8,4,2)

sh(6,4,2) = pc(9,4,2)

sh(7,4,2) = pc(10,4,2)

sh(8,4,2) = pc(5,4,2) - pc(6,4,2)

sh(9,4,2) = 2*pc(7,4,2) - pc(5,4,2) - pc(6,4,2)

sh(10,4,2) = pc(11,4,2)

sh(11,4,2) = pc(13,4,2) - pc(15,4,2)

sh(12,4,2) = pc(18,4,2) - 3*pc(14,4,2)

sh(13,4,2) = 3*pc(12,4,2) - pc(19,4,2)

sh(14,4,2) = 2*pc(20,4,2) - 3*pc(13,4,2) - 3*pc(15,4,2)

sh(15,4,2) = 4*pc(16,4,2) - pc(18,4,2) - pc(14,4,2)

sh(16,4,2) = 4*pc(17,4,2) - pc(12,4,2) - pc(19,4,2)

sh(1,5,2) = pc(1,8,2)

sh(2,5,2) = pc(2,8,2)

sh(3,5,2) = pc(3,8,2)

sh(4,5,2) = pc(4,8,2)

sh(5,5,2) = pc(8,8,2)

sh(6,5,2) = pc(9,8,2)

sh(7,5,2) = pc(10,8,2)

sh(8,5,2) = pc(5,8,2) - pc(6,8,2)

sh(9,5,2) = 2*pc(7,8,2) - pc(5,8,2) - pc(6,8,2)

sh(10,5,2) = pc(11,8,2)

sh(11,5,2) = pc(13,8,2) - pc(15,8,2)

sh(12,5,2) = pc(18,8,2) - 3*pc(14,8,2)

sh(13,5,2) = 3*pc(12,8,2) - pc(19,8,2)

sh(14,5,2) = 2*pc(20,8,2) - 3*pc(13,8,2) - 3*pc(15,8,2)

sh(15,5,2) = 4*pc(16,8,2) - pc(18,8,2) - pc(14,8,2)

sh(16,5,2) = 4*pc(17,8,2) - pc(12,8,2) - pc(19,8,2)

sh(1,6,2) = pc(1,9,2)

sh(2,6,2) = pc(2,9,2)

sh(3,6,2) = pc(3,9,2)

sh(4,6,2) = pc(4,9,2)

sh(5,6,2) = pc(8,9,2)

sh(6,6,2) = pc(9,9,2)

sh(7,6,2) = pc(10,9,2)

sh(8,6,2) = pc(5,9,2) - pc(6,9,2)

sh(9,6,2) = 2*pc(7,9,2) - pc(5,9,2) - pc(6,9,2)

sh(10,6,2) = pc(11,9,2)

sh(11,6,2) = pc(13,9,2) - pc(15,9,2)

sh(12,6,2) = pc(18,9,2) - 3*pc(14,9,2)

sh(13,6,2) = 3*pc(12,9,2) - pc(19,9,2)

sh(14,6,2) = 2*pc(20,9,2) - 3*pc(13,9,2) - 3*pc(15,9,2)

sh(15,6,2) = 4*pc(16,9,2) - pc(18,9,2) - pc(14,9,2)

sh(16,6,2) = 4*pc(17,9,2) - pc(12,9,2) - pc(19,9,2)

sh(1,7,2) = pc(1,10,2)

sh(2,7,2) = pc(2,10,2)

sh(3,7,2) = pc(3,10,2)

sh(4,7,2) = pc(4,10,2)

sh(5,7,2) = pc(8,10,2)

sh(6,7,2) = pc(9,10,2)

sh(7,7,2) = pc(10,10,2)

sh(8,7,2) = pc(5,10,2) - pc(6,10,2)

sh(9,7,2) = 2*pc(7,10,2) - pc(5,10,2) - pc(6,10,2)

sh(10,7,2) = pc(11,10,2)

sh(11,7,2) = pc(13,10,2) - pc(15,10,2)

sh(12,7,2) = pc(18,10,2) - 3*pc(14,10,2)

sh(13,7,2) = 3*pc(12,10,2) - pc(19,10,2)

sh(14,7,2) = 2*pc(20,10,2) - 3*pc(13,10,2) - 3*pc(15,10,2)

sh(15,7,2) = 4*pc(16,10,2) - pc(18,10,2) - pc(14,10,2)

sh(16,7,2) = 4*pc(17,10,2) - pc(12,10,2) - pc(19,10,2)

sh(1,8,2) = pc(1,5,2) - pc(1,6,2)

sh(2,8,2) = pc(2,5,2) - pc(2,6,2)

sh(3,8,2) = pc(3,5,2) - pc(3,6,2)

sh(4,8,2) = pc(4,5,2) - pc(4,6,2)

sh(5,8,2) = pc(8,5,2) - pc(8,6,2)

sh(6,8,2) = pc(9,5,2) - pc(9,6,2)

sh(7,8,2) = pc(10,5,2) - pc(10,6,2)

sh(8,8,2) = pc(5,5,2) - pc(5,6,2) - pc(6,5,2) + pc(6,6,2)

sh(9,8,2) = 2*pc(7,5,2) - 2*pc(7,6,2) - pc(5,5,2) + pc(5,6,2) - pc(6,5,2) + pc&
&(6,6,2)

sh(10,8,2) = pc(11,5,2) - pc(11,6,2)

sh(11,8,2) = pc(13,5,2) - pc(13,6,2) - pc(15,5,2) + pc(15,6,2)

sh(12,8,2) = pc(18,5,2) - pc(18,6,2) - 3*pc(14,5,2) + 3*pc(14,6,2)

sh(13,8,2) = 3*pc(12,5,2) - 3*pc(12,6,2) - pc(19,5,2) + pc(19,6,2)

sh(14,8,2) = 2*pc(20,5,2) - 2*pc(20,6,2) - 3*pc(13,5,2) + 3*pc(13,6,2) - 3*pc(&
&15,5,2) + 3*pc(15,6,2)

sh(15,8,2) = 4*pc(16,5,2) - 4*pc(16,6,2) - pc(18,5,2) + pc(18,6,2) - pc(14,5,2&
&) + pc(14,6,2)

sh(16,8,2) = 4*pc(17,5,2) - 4*pc(17,6,2) - pc(12,5,2) + pc(12,6,2) - pc(19,5,2&
&) + pc(19,6,2)

sh(1,9,2) = 2*pc(1,7,2) - pc(1,5,2) - pc(1,6,2)

sh(2,9,2) = 2*pc(2,7,2) - pc(2,5,2) - pc(2,6,2)

sh(3,9,2) = 2*pc(3,7,2) - pc(3,5,2) - pc(3,6,2)

sh(4,9,2) = 2*pc(4,7,2) - pc(4,5,2) - pc(4,6,2)

sh(5,9,2) = 2*pc(8,7,2) - pc(8,5,2) - pc(8,6,2)

sh(6,9,2) = 2*pc(9,7,2) - pc(9,5,2) - pc(9,6,2)

sh(7,9,2) = 2*pc(10,7,2) - pc(10,5,2) - pc(10,6,2)

sh(8,9,2) = 2*pc(5,7,2) - pc(5,5,2) - pc(5,6,2) - 2*pc(6,7,2) + pc(6,5,2) + pc&
&(6,6,2)

sh(9,9,2) = 4*pc(7,7,2) - 2*pc(7,5,2) - 2*pc(7,6,2) - 2*pc(5,7,2) + pc(5,5,2) &
&+ pc(5,6,2) - 2*pc(6,7,2) + pc(6,5,2) + pc(6,6,2)

sh(10,9,2) = 2*pc(11,7,2) - pc(11,5,2) - pc(11,6,2)

sh(11,9,2) = 2*pc(13,7,2) - pc(13,5,2) - pc(13,6,2) - 2*pc(15,7,2) + pc(15,5,2&
&) + pc(15,6,2)

sh(12,9,2) = 2*pc(18,7,2) - pc(18,5,2) - pc(18,6,2) - 6*pc(14,7,2) + 3*pc(14,5&
&,2) + 3*pc(14,6,2)

sh(13,9,2) = 6*pc(12,7,2) - 3*pc(12,5,2) - 3*pc(12,6,2) - 2*pc(19,7,2) + pc(19&
&,5,2) + pc(19,6,2)

sh(14,9,2) = 4*pc(20,7,2) - 2*pc(20,5,2) - 2*pc(20,6,2) - 6*pc(13,7,2) + 3*pc(&
&13,5,2) + 3*pc(13,6,2) - 6*pc(15,7,2) + 3*pc(15,5,2) + 3*pc(15,6,2)

sh(15,9,2) = 8*pc(16,7,2) - 4*pc(16,5,2) - 4*pc(16,6,2) - 2*pc(18,7,2) + pc(18&
&,5,2) + pc(18,6,2) - 2*pc(14,7,2) + pc(14,5,2) + pc(14,6,2)

sh(16,9,2) = 8*pc(17,7,2) - 4*pc(17,5,2) - 4*pc(17,6,2) - 2*pc(12,7,2) + pc(12&
&,5,2) + pc(12,6,2) - 2*pc(19,7,2) + pc(19,5,2) + pc(19,6,2)

sh(1,10,2) = pc(1,11,2)

sh(2,10,2) = pc(2,11,2)

sh(3,10,2) = pc(3,11,2)

sh(4,10,2) = pc(4,11,2)

sh(5,10,2) = pc(8,11,2)

sh(6,10,2) = pc(9,11,2)

sh(7,10,2) = pc(10,11,2)

sh(8,10,2) = pc(5,11,2) - pc(6,11,2)

sh(9,10,2) = 2*pc(7,11,2) - pc(5,11,2) - pc(6,11,2)

sh(10,10,2) = pc(11,11,2)

sh(11,10,2) = pc(13,11,2) - pc(15,11,2)

sh(12,10,2) = pc(18,11,2) - 3*pc(14,11,2)

sh(13,10,2) = 3*pc(12,11,2) - pc(19,11,2)

sh(14,10,2) = 2*pc(20,11,2) - 3*pc(13,11,2) - 3*pc(15,11,2)

sh(15,10,2) = 4*pc(16,11,2) - pc(18,11,2) - pc(14,11,2)

sh(16,10,2) = 4*pc(17,11,2) - pc(12,11,2) - pc(19,11,2)

sh(1,11,2) = pc(1,13,2) - pc(1,15,2)

sh(2,11,2) = pc(2,13,2) - pc(2,15,2)

sh(3,11,2) = pc(3,13,2) - pc(3,15,2)

sh(4,11,2) = pc(4,13,2) - pc(4,15,2)

sh(5,11,2) = pc(8,13,2) - pc(8,15,2)

sh(6,11,2) = pc(9,13,2) - pc(9,15,2)

sh(7,11,2) = pc(10,13,2) - pc(10,15,2)

sh(8,11,2) = pc(5,13,2) - pc(5,15,2) - pc(6,13,2) + pc(6,15,2)

sh(9,11,2) = 2*pc(7,13,2) - 2*pc(7,15,2) - pc(5,13,2) + pc(5,15,2) - pc(6,13,2&
&) + pc(6,15,2)

sh(10,11,2) = pc(11,13,2) - pc(11,15,2)

sh(11,11,2) = pc(13,13,2) - pc(13,15,2) - pc(15,13,2) + pc(15,15,2)

sh(12,11,2) = pc(18,13,2) - pc(18,15,2) - 3*pc(14,13,2) + 3*pc(14,15,2)

sh(13,11,2) = 3*pc(12,13,2) - 3*pc(12,15,2) - pc(19,13,2) + pc(19,15,2)

sh(14,11,2) = 2*pc(20,13,2) - 2*pc(20,15,2) - 3*pc(13,13,2) + 3*pc(13,15,2) - &
&3*pc(15,13,2) + 3*pc(15,15,2)

sh(15,11,2) = 4*pc(16,13,2) - 4*pc(16,15,2) - pc(18,13,2) + pc(18,15,2) - pc(1&
&4,13,2) + pc(14,15,2)

sh(16,11,2) = 4*pc(17,13,2) - 4*pc(17,15,2) - pc(12,13,2) + pc(12,15,2) - pc(1&
&9,13,2) + pc(19,15,2)

sh(1,12,2) = pc(1,18,2) - 3*pc(1,14,2)

sh(2,12,2) = pc(2,18,2) - 3*pc(2,14,2)

sh(3,12,2) = pc(3,18,2) - 3*pc(3,14,2)

sh(4,12,2) = pc(4,18,2) - 3*pc(4,14,2)

sh(5,12,2) = pc(8,18,2) - 3*pc(8,14,2)

sh(6,12,2) = pc(9,18,2) - 3*pc(9,14,2)

sh(7,12,2) = pc(10,18,2) - 3*pc(10,14,2)

sh(8,12,2) = pc(5,18,2) - 3*pc(5,14,2) - pc(6,18,2) + 3*pc(6,14,2)

sh(9,12,2) = 2*pc(7,18,2) - 6*pc(7,14,2) - pc(5,18,2) + 3*pc(5,14,2) - pc(6,18&
&,2) + 3*pc(6,14,2)

sh(10,12,2) = pc(11,18,2) - 3*pc(11,14,2)

sh(11,12,2) = pc(13,18,2) - 3*pc(13,14,2) - pc(15,18,2) + 3*pc(15,14,2)

sh(12,12,2) = pc(18,18,2) - 3*pc(18,14,2) - 3*pc(14,18,2) + 9*pc(14,14,2)

sh(13,12,2) = 3*pc(12,18,2) - 9*pc(12,14,2) - pc(19,18,2) + 3*pc(19,14,2)

sh(14,12,2) = 2*pc(20,18,2) - 6*pc(20,14,2) - 3*pc(13,18,2) + 9*pc(13,14,2) - &
&3*pc(15,18,2) + 9*pc(15,14,2)

sh(15,12,2) = 4*pc(16,18,2) - 12*pc(16,14,2) - pc(18,18,2) + 3*pc(18,14,2) - p&
&c(14,18,2) + 3*pc(14,14,2)

sh(16,12,2) = 4*pc(17,18,2) - 12*pc(17,14,2) - pc(12,18,2) + 3*pc(12,14,2) - p&
&c(19,18,2) + 3*pc(19,14,2)

sh(1,13,2) = 3*pc(1,12,2) - pc(1,19,2)

sh(2,13,2) = 3*pc(2,12,2) - pc(2,19,2)

sh(3,13,2) = 3*pc(3,12,2) - pc(3,19,2)

sh(4,13,2) = 3*pc(4,12,2) - pc(4,19,2)

sh(5,13,2) = 3*pc(8,12,2) - pc(8,19,2)

sh(6,13,2) = 3*pc(9,12,2) - pc(9,19,2)

sh(7,13,2) = 3*pc(10,12,2) - pc(10,19,2)

sh(8,13,2) = 3*pc(5,12,2) - pc(5,19,2) - 3*pc(6,12,2) + pc(6,19,2)

sh(9,13,2) = 6*pc(7,12,2) - 2*pc(7,19,2) - 3*pc(5,12,2) + pc(5,19,2) - 3*pc(6,&
&12,2) + pc(6,19,2)

sh(10,13,2) = 3*pc(11,12,2) - pc(11,19,2)

sh(11,13,2) = 3*pc(13,12,2) - pc(13,19,2) - 3*pc(15,12,2) + pc(15,19,2)

sh(12,13,2) = 3*pc(18,12,2) - pc(18,19,2) - 9*pc(14,12,2) + 3*pc(14,19,2)

sh(13,13,2) = 9*pc(12,12,2) - 3*pc(12,19,2) - 3*pc(19,12,2) + pc(19,19,2)

sh(14,13,2) = 6*pc(20,12,2) - 2*pc(20,19,2) - 9*pc(13,12,2) + 3*pc(13,19,2) - &
&9*pc(15,12,2) + 3*pc(15,19,2)

sh(15,13,2) = 12*pc(16,12,2) - 4*pc(16,19,2) - 3*pc(18,12,2) + pc(18,19,2) - 3&
&*pc(14,12,2) + pc(14,19,2)

sh(16,13,2) = 12*pc(17,12,2) - 4*pc(17,19,2) - 3*pc(12,12,2) + pc(12,19,2) - 3&
&*pc(19,12,2) + pc(19,19,2)

sh(1,14,2) = 2*pc(1,20,2) - 3*pc(1,13,2) - 3*pc(1,15,2)

sh(2,14,2) = 2*pc(2,20,2) - 3*pc(2,13,2) - 3*pc(2,15,2)

sh(3,14,2) = 2*pc(3,20,2) - 3*pc(3,13,2) - 3*pc(3,15,2)

sh(4,14,2) = 2*pc(4,20,2) - 3*pc(4,13,2) - 3*pc(4,15,2)

sh(5,14,2) = 2*pc(8,20,2) - 3*pc(8,13,2) - 3*pc(8,15,2)

sh(6,14,2) = 2*pc(9,20,2) - 3*pc(9,13,2) - 3*pc(9,15,2)

sh(7,14,2) = 2*pc(10,20,2) - 3*pc(10,13,2) - 3*pc(10,15,2)

sh(8,14,2) = 2*pc(5,20,2) - 3*pc(5,13,2) - 3*pc(5,15,2) - 2*pc(6,20,2) + 3*pc(&
&6,13,2) + 3*pc(6,15,2)

sh(9,14,2) = 4*pc(7,20,2) - 6*pc(7,13,2) - 6*pc(7,15,2) - 2*pc(5,20,2) + 3*pc(&
&5,13,2) + 3*pc(5,15,2) - 2*pc(6,20,2) + 3*pc(6,13,2) + 3*pc(6,15,2)

sh(10,14,2) = 2*pc(11,20,2) - 3*pc(11,13,2) - 3*pc(11,15,2)

sh(11,14,2) = 2*pc(13,20,2) - 3*pc(13,13,2) - 3*pc(13,15,2) - 2*pc(15,20,2) + &
&3*pc(15,13,2) + 3*pc(15,15,2)

sh(12,14,2) = 2*pc(18,20,2) - 3*pc(18,13,2) - 3*pc(18,15,2) - 6*pc(14,20,2) + &
&9*pc(14,13,2) + 9*pc(14,15,2)

sh(13,14,2) = 6*pc(12,20,2) - 9*pc(12,13,2) - 9*pc(12,15,2) - 2*pc(19,20,2) + &
&3*pc(19,13,2) + 3*pc(19,15,2)

sh(14,14,2) = 4*pc(20,20,2) - 6*pc(20,13,2) - 6*pc(20,15,2) - 6*pc(13,20,2) + &
&9*pc(13,13,2) + 9*pc(13,15,2) - 6*pc(15,20,2) + 9*pc(15,13,2) + 9*pc(15,15,2)

sh(15,14,2) = 8*pc(16,20,2) - 12*pc(16,13,2) - 12*pc(16,15,2) - 2*pc(18,20,2) &
&+ 3*pc(18,13,2) + 3*pc(18,15,2) - 2*pc(14,20,2) + 3*pc(14,13,2) + 3*pc(14,15,2&
&)

sh(16,14,2) = 8*pc(17,20,2) - 12*pc(17,13,2) - 12*pc(17,15,2) - 2*pc(12,20,2) &
&+ 3*pc(12,13,2) + 3*pc(12,15,2) - 2*pc(19,20,2) + 3*pc(19,13,2) + 3*pc(19,15,2&
&)

sh(1,15,2) = 4*pc(1,16,2) - pc(1,18,2) - pc(1,14,2)

sh(2,15,2) = 4*pc(2,16,2) - pc(2,18,2) - pc(2,14,2)

sh(3,15,2) = 4*pc(3,16,2) - pc(3,18,2) - pc(3,14,2)

sh(4,15,2) = 4*pc(4,16,2) - pc(4,18,2) - pc(4,14,2)

sh(5,15,2) = 4*pc(8,16,2) - pc(8,18,2) - pc(8,14,2)

sh(6,15,2) = 4*pc(9,16,2) - pc(9,18,2) - pc(9,14,2)

sh(7,15,2) = 4*pc(10,16,2) - pc(10,18,2) - pc(10,14,2)

sh(8,15,2) = 4*pc(5,16,2) - pc(5,18,2) - pc(5,14,2) - 4*pc(6,16,2) + pc(6,18,2&
&) + pc(6,14,2)

sh(9,15,2) = 8*pc(7,16,2) - 2*pc(7,18,2) - 2*pc(7,14,2) - 4*pc(5,16,2) + pc(5,&
&18,2) + pc(5,14,2) - 4*pc(6,16,2) + pc(6,18,2) + pc(6,14,2)

sh(10,15,2) = 4*pc(11,16,2) - pc(11,18,2) - pc(11,14,2)

sh(11,15,2) = 4*pc(13,16,2) - pc(13,18,2) - pc(13,14,2) - 4*pc(15,16,2) + pc(1&
&5,18,2) + pc(15,14,2)

sh(12,15,2) = 4*pc(18,16,2) - pc(18,18,2) - pc(18,14,2) - 12*pc(14,16,2) + 3*p&
&c(14,18,2) + 3*pc(14,14,2)

sh(13,15,2) = 12*pc(12,16,2) - 3*pc(12,18,2) - 3*pc(12,14,2) - 4*pc(19,16,2) +&
& pc(19,18,2) + pc(19,14,2)

sh(14,15,2) = 8*pc(20,16,2) - 2*pc(20,18,2) - 2*pc(20,14,2) - 12*pc(13,16,2) +&
& 3*pc(13,18,2) + 3*pc(13,14,2) - 12*pc(15,16,2) + 3*pc(15,18,2) + 3*pc(15,14,2&
&)

sh(15,15,2) = 16*pc(16,16,2) - 4*pc(16,18,2) - 4*pc(16,14,2) - 4*pc(18,16,2) +&
& pc(18,18,2) + pc(18,14,2) - 4*pc(14,16,2) + pc(14,18,2) + pc(14,14,2)

sh(16,15,2) = 16*pc(17,16,2) - 4*pc(17,18,2) - 4*pc(17,14,2) - 4*pc(12,16,2) +&
& pc(12,18,2) + pc(12,14,2) - 4*pc(19,16,2) + pc(19,18,2) + pc(19,14,2)

sh(1,16,2) = 4*pc(1,17,2) - pc(1,12,2) - pc(1,19,2)

sh(2,16,2) = 4*pc(2,17,2) - pc(2,12,2) - pc(2,19,2)

sh(3,16,2) = 4*pc(3,17,2) - pc(3,12,2) - pc(3,19,2)

sh(4,16,2) = 4*pc(4,17,2) - pc(4,12,2) - pc(4,19,2)

sh(5,16,2) = 4*pc(8,17,2) - pc(8,12,2) - pc(8,19,2)

sh(6,16,2) = 4*pc(9,17,2) - pc(9,12,2) - pc(9,19,2)

sh(7,16,2) = 4*pc(10,17,2) - pc(10,12,2) - pc(10,19,2)

sh(8,16,2) = 4*pc(5,17,2) - pc(5,12,2) - pc(5,19,2) - 4*pc(6,17,2) + pc(6,12,2&
&) + pc(6,19,2)

sh(9,16,2) = 8*pc(7,17,2) - 2*pc(7,12,2) - 2*pc(7,19,2) - 4*pc(5,17,2) + pc(5,&
&12,2) + pc(5,19,2) - 4*pc(6,17,2) + pc(6,12,2) + pc(6,19,2)

sh(10,16,2) = 4*pc(11,17,2) - pc(11,12,2) - pc(11,19,2)

sh(11,16,2) = 4*pc(13,17,2) - pc(13,12,2) - pc(13,19,2) - 4*pc(15,17,2) + pc(1&
&5,12,2) + pc(15,19,2)

sh(12,16,2) = 4*pc(18,17,2) - pc(18,12,2) - pc(18,19,2) - 12*pc(14,17,2) + 3*p&
&c(14,12,2) + 3*pc(14,19,2)

sh(13,16,2) = 12*pc(12,17,2) - 3*pc(12,12,2) - 3*pc(12,19,2) - 4*pc(19,17,2) +&
& pc(19,12,2) + pc(19,19,2)

sh(14,16,2) = 8*pc(20,17,2) - 2*pc(20,12,2) - 2*pc(20,19,2) - 12*pc(13,17,2) +&
& 3*pc(13,12,2) + 3*pc(13,19,2) - 12*pc(15,17,2) + 3*pc(15,12,2) + 3*pc(15,19,2&
&)

sh(15,16,2) = 16*pc(16,17,2) - 4*pc(16,12,2) - 4*pc(16,19,2) - 4*pc(18,17,2) +&
& pc(18,12,2) + pc(18,19,2) - 4*pc(14,17,2) + pc(14,12,2) + pc(14,19,2)

sh(16,16,2) = 16*pc(17,17,2) - 4*pc(17,12,2) - 4*pc(17,19,2) - 4*pc(12,17,2) +&
& pc(12,12,2) + pc(12,19,2) - 4*pc(19,17,2) + pc(19,12,2) + pc(19,19,2)

sh(1,1,3) = pc(1,1,3)

sh(2,1,3) = pc(2,1,3)

sh(3,1,3) = pc(3,1,3)

sh(4,1,3) = pc(4,1,3)

sh(5,1,3) = pc(8,1,3)

sh(6,1,3) = pc(9,1,3)

sh(7,1,3) = pc(10,1,3)

sh(8,1,3) = pc(5,1,3) - pc(6,1,3)

sh(9,1,3) = 2*pc(7,1,3) - pc(5,1,3) - pc(6,1,3)

sh(10,1,3) = pc(11,1,3)

sh(11,1,3) = pc(13,1,3) - pc(15,1,3)

sh(12,1,3) = pc(18,1,3) - 3*pc(14,1,3)

sh(13,1,3) = 3*pc(12,1,3) - pc(19,1,3)

sh(14,1,3) = 2*pc(20,1,3) - 3*pc(13,1,3) - 3*pc(15,1,3)

sh(15,1,3) = 4*pc(16,1,3) - pc(18,1,3) - pc(14,1,3)

sh(16,1,3) = 4*pc(17,1,3) - pc(12,1,3) - pc(19,1,3)

sh(1,2,3) = pc(1,2,3)

sh(2,2,3) = pc(2,2,3)

sh(3,2,3) = pc(3,2,3)

sh(4,2,3) = pc(4,2,3)

sh(5,2,3) = pc(8,2,3)

sh(6,2,3) = pc(9,2,3)

sh(7,2,3) = pc(10,2,3)

sh(8,2,3) = pc(5,2,3) - pc(6,2,3)

sh(9,2,3) = 2*pc(7,2,3) - pc(5,2,3) - pc(6,2,3)

sh(10,2,3) = pc(11,2,3)

sh(11,2,3) = pc(13,2,3) - pc(15,2,3)

sh(12,2,3) = pc(18,2,3) - 3*pc(14,2,3)

sh(13,2,3) = 3*pc(12,2,3) - pc(19,2,3)

sh(14,2,3) = 2*pc(20,2,3) - 3*pc(13,2,3) - 3*pc(15,2,3)

sh(15,2,3) = 4*pc(16,2,3) - pc(18,2,3) - pc(14,2,3)

sh(16,2,3) = 4*pc(17,2,3) - pc(12,2,3) - pc(19,2,3)

sh(1,3,3) = pc(1,3,3)

sh(2,3,3) = pc(2,3,3)

sh(3,3,3) = pc(3,3,3)

sh(4,3,3) = pc(4,3,3)

sh(5,3,3) = pc(8,3,3)

sh(6,3,3) = pc(9,3,3)

sh(7,3,3) = pc(10,3,3)

sh(8,3,3) = pc(5,3,3) - pc(6,3,3)

sh(9,3,3) = 2*pc(7,3,3) - pc(5,3,3) - pc(6,3,3)

sh(10,3,3) = pc(11,3,3)

sh(11,3,3) = pc(13,3,3) - pc(15,3,3)

sh(12,3,3) = pc(18,3,3) - 3*pc(14,3,3)

sh(13,3,3) = 3*pc(12,3,3) - pc(19,3,3)

sh(14,3,3) = 2*pc(20,3,3) - 3*pc(13,3,3) - 3*pc(15,3,3)

sh(15,3,3) = 4*pc(16,3,3) - pc(18,3,3) - pc(14,3,3)

sh(16,3,3) = 4*pc(17,3,3) - pc(12,3,3) - pc(19,3,3)

sh(1,4,3) = pc(1,4,3)

sh(2,4,3) = pc(2,4,3)

sh(3,4,3) = pc(3,4,3)

sh(4,4,3) = pc(4,4,3)

sh(5,4,3) = pc(8,4,3)

sh(6,4,3) = pc(9,4,3)

sh(7,4,3) = pc(10,4,3)

sh(8,4,3) = pc(5,4,3) - pc(6,4,3)

sh(9,4,3) = 2*pc(7,4,3) - pc(5,4,3) - pc(6,4,3)

sh(10,4,3) = pc(11,4,3)

sh(11,4,3) = pc(13,4,3) - pc(15,4,3)

sh(12,4,3) = pc(18,4,3) - 3*pc(14,4,3)

sh(13,4,3) = 3*pc(12,4,3) - pc(19,4,3)

sh(14,4,3) = 2*pc(20,4,3) - 3*pc(13,4,3) - 3*pc(15,4,3)

sh(15,4,3) = 4*pc(16,4,3) - pc(18,4,3) - pc(14,4,3)

sh(16,4,3) = 4*pc(17,4,3) - pc(12,4,3) - pc(19,4,3)

sh(1,5,3) = pc(1,8,3)

sh(2,5,3) = pc(2,8,3)

sh(3,5,3) = pc(3,8,3)

sh(4,5,3) = pc(4,8,3)

sh(5,5,3) = pc(8,8,3)

sh(6,5,3) = pc(9,8,3)

sh(7,5,3) = pc(10,8,3)

sh(8,5,3) = pc(5,8,3) - pc(6,8,3)

sh(9,5,3) = 2*pc(7,8,3) - pc(5,8,3) - pc(6,8,3)

sh(10,5,3) = pc(11,8,3)

sh(11,5,3) = pc(13,8,3) - pc(15,8,3)

sh(12,5,3) = pc(18,8,3) - 3*pc(14,8,3)

sh(13,5,3) = 3*pc(12,8,3) - pc(19,8,3)

sh(14,5,3) = 2*pc(20,8,3) - 3*pc(13,8,3) - 3*pc(15,8,3)

sh(15,5,3) = 4*pc(16,8,3) - pc(18,8,3) - pc(14,8,3)

sh(16,5,3) = 4*pc(17,8,3) - pc(12,8,3) - pc(19,8,3)

sh(1,6,3) = pc(1,9,3)

sh(2,6,3) = pc(2,9,3)

sh(3,6,3) = pc(3,9,3)

sh(4,6,3) = pc(4,9,3)

sh(5,6,3) = pc(8,9,3)

sh(6,6,3) = pc(9,9,3)

sh(7,6,3) = pc(10,9,3)

sh(8,6,3) = pc(5,9,3) - pc(6,9,3)

sh(9,6,3) = 2*pc(7,9,3) - pc(5,9,3) - pc(6,9,3)

sh(10,6,3) = pc(11,9,3)

sh(11,6,3) = pc(13,9,3) - pc(15,9,3)

sh(12,6,3) = pc(18,9,3) - 3*pc(14,9,3)

sh(13,6,3) = 3*pc(12,9,3) - pc(19,9,3)

sh(14,6,3) = 2*pc(20,9,3) - 3*pc(13,9,3) - 3*pc(15,9,3)

sh(15,6,3) = 4*pc(16,9,3) - pc(18,9,3) - pc(14,9,3)

sh(16,6,3) = 4*pc(17,9,3) - pc(12,9,3) - pc(19,9,3)

sh(1,7,3) = pc(1,10,3)

sh(2,7,3) = pc(2,10,3)

sh(3,7,3) = pc(3,10,3)

sh(4,7,3) = pc(4,10,3)

sh(5,7,3) = pc(8,10,3)

sh(6,7,3) = pc(9,10,3)

sh(7,7,3) = pc(10,10,3)

sh(8,7,3) = pc(5,10,3) - pc(6,10,3)

sh(9,7,3) = 2*pc(7,10,3) - pc(5,10,3) - pc(6,10,3)

sh(10,7,3) = pc(11,10,3)

sh(11,7,3) = pc(13,10,3) - pc(15,10,3)

sh(12,7,3) = pc(18,10,3) - 3*pc(14,10,3)

sh(13,7,3) = 3*pc(12,10,3) - pc(19,10,3)

sh(14,7,3) = 2*pc(20,10,3) - 3*pc(13,10,3) - 3*pc(15,10,3)

sh(15,7,3) = 4*pc(16,10,3) - pc(18,10,3) - pc(14,10,3)

sh(16,7,3) = 4*pc(17,10,3) - pc(12,10,3) - pc(19,10,3)

sh(1,8,3) = pc(1,5,3) - pc(1,6,3)

sh(2,8,3) = pc(2,5,3) - pc(2,6,3)

sh(3,8,3) = pc(3,5,3) - pc(3,6,3)

sh(4,8,3) = pc(4,5,3) - pc(4,6,3)

sh(5,8,3) = pc(8,5,3) - pc(8,6,3)

sh(6,8,3) = pc(9,5,3) - pc(9,6,3)

sh(7,8,3) = pc(10,5,3) - pc(10,6,3)

sh(8,8,3) = pc(5,5,3) - pc(5,6,3) - pc(6,5,3) + pc(6,6,3)

sh(9,8,3) = 2*pc(7,5,3) - 2*pc(7,6,3) - pc(5,5,3) + pc(5,6,3) - pc(6,5,3) + pc&
&(6,6,3)

sh(10,8,3) = pc(11,5,3) - pc(11,6,3)

sh(11,8,3) = pc(13,5,3) - pc(13,6,3) - pc(15,5,3) + pc(15,6,3)

sh(12,8,3) = pc(18,5,3) - pc(18,6,3) - 3*pc(14,5,3) + 3*pc(14,6,3)

sh(13,8,3) = 3*pc(12,5,3) - 3*pc(12,6,3) - pc(19,5,3) + pc(19,6,3)

sh(14,8,3) = 2*pc(20,5,3) - 2*pc(20,6,3) - 3*pc(13,5,3) + 3*pc(13,6,3) - 3*pc(&
&15,5,3) + 3*pc(15,6,3)

sh(15,8,3) = 4*pc(16,5,3) - 4*pc(16,6,3) - pc(18,5,3) + pc(18,6,3) - pc(14,5,3&
&) + pc(14,6,3)

sh(16,8,3) = 4*pc(17,5,3) - 4*pc(17,6,3) - pc(12,5,3) + pc(12,6,3) - pc(19,5,3&
&) + pc(19,6,3)

sh(1,9,3) = 2*pc(1,7,3) - pc(1,5,3) - pc(1,6,3)

sh(2,9,3) = 2*pc(2,7,3) - pc(2,5,3) - pc(2,6,3)

sh(3,9,3) = 2*pc(3,7,3) - pc(3,5,3) - pc(3,6,3)

sh(4,9,3) = 2*pc(4,7,3) - pc(4,5,3) - pc(4,6,3)

sh(5,9,3) = 2*pc(8,7,3) - pc(8,5,3) - pc(8,6,3)

sh(6,9,3) = 2*pc(9,7,3) - pc(9,5,3) - pc(9,6,3)

sh(7,9,3) = 2*pc(10,7,3) - pc(10,5,3) - pc(10,6,3)

sh(8,9,3) = 2*pc(5,7,3) - pc(5,5,3) - pc(5,6,3) - 2*pc(6,7,3) + pc(6,5,3) + pc&
&(6,6,3)

sh(9,9,3) = 4*pc(7,7,3) - 2*pc(7,5,3) - 2*pc(7,6,3) - 2*pc(5,7,3) + pc(5,5,3) &
&+ pc(5,6,3) - 2*pc(6,7,3) + pc(6,5,3) + pc(6,6,3)

sh(10,9,3) = 2*pc(11,7,3) - pc(11,5,3) - pc(11,6,3)

sh(11,9,3) = 2*pc(13,7,3) - pc(13,5,3) - pc(13,6,3) - 2*pc(15,7,3) + pc(15,5,3&
&) + pc(15,6,3)

sh(12,9,3) = 2*pc(18,7,3) - pc(18,5,3) - pc(18,6,3) - 6*pc(14,7,3) + 3*pc(14,5&
&,3) + 3*pc(14,6,3)

sh(13,9,3) = 6*pc(12,7,3) - 3*pc(12,5,3) - 3*pc(12,6,3) - 2*pc(19,7,3) + pc(19&
&,5,3) + pc(19,6,3)

sh(14,9,3) = 4*pc(20,7,3) - 2*pc(20,5,3) - 2*pc(20,6,3) - 6*pc(13,7,3) + 3*pc(&
&13,5,3) + 3*pc(13,6,3) - 6*pc(15,7,3) + 3*pc(15,5,3) + 3*pc(15,6,3)

sh(15,9,3) = 8*pc(16,7,3) - 4*pc(16,5,3) - 4*pc(16,6,3) - 2*pc(18,7,3) + pc(18&
&,5,3) + pc(18,6,3) - 2*pc(14,7,3) + pc(14,5,3) + pc(14,6,3)

sh(16,9,3) = 8*pc(17,7,3) - 4*pc(17,5,3) - 4*pc(17,6,3) - 2*pc(12,7,3) + pc(12&
&,5,3) + pc(12,6,3) - 2*pc(19,7,3) + pc(19,5,3) + pc(19,6,3)

sh(1,10,3) = pc(1,11,3)

sh(2,10,3) = pc(2,11,3)

sh(3,10,3) = pc(3,11,3)

sh(4,10,3) = pc(4,11,3)

sh(5,10,3) = pc(8,11,3)

sh(6,10,3) = pc(9,11,3)

sh(7,10,3) = pc(10,11,3)

sh(8,10,3) = pc(5,11,3) - pc(6,11,3)

sh(9,10,3) = 2*pc(7,11,3) - pc(5,11,3) - pc(6,11,3)

sh(10,10,3) = pc(11,11,3)

sh(11,10,3) = pc(13,11,3) - pc(15,11,3)

sh(12,10,3) = pc(18,11,3) - 3*pc(14,11,3)

sh(13,10,3) = 3*pc(12,11,3) - pc(19,11,3)

sh(14,10,3) = 2*pc(20,11,3) - 3*pc(13,11,3) - 3*pc(15,11,3)

sh(15,10,3) = 4*pc(16,11,3) - pc(18,11,3) - pc(14,11,3)

sh(16,10,3) = 4*pc(17,11,3) - pc(12,11,3) - pc(19,11,3)

sh(1,11,3) = pc(1,13,3) - pc(1,15,3)

sh(2,11,3) = pc(2,13,3) - pc(2,15,3)

sh(3,11,3) = pc(3,13,3) - pc(3,15,3)

sh(4,11,3) = pc(4,13,3) - pc(4,15,3)

sh(5,11,3) = pc(8,13,3) - pc(8,15,3)

sh(6,11,3) = pc(9,13,3) - pc(9,15,3)

sh(7,11,3) = pc(10,13,3) - pc(10,15,3)

sh(8,11,3) = pc(5,13,3) - pc(5,15,3) - pc(6,13,3) + pc(6,15,3)

sh(9,11,3) = 2*pc(7,13,3) - 2*pc(7,15,3) - pc(5,13,3) + pc(5,15,3) - pc(6,13,3&
&) + pc(6,15,3)

sh(10,11,3) = pc(11,13,3) - pc(11,15,3)

sh(11,11,3) = pc(13,13,3) - pc(13,15,3) - pc(15,13,3) + pc(15,15,3)

sh(12,11,3) = pc(18,13,3) - pc(18,15,3) - 3*pc(14,13,3) + 3*pc(14,15,3)

sh(13,11,3) = 3*pc(12,13,3) - 3*pc(12,15,3) - pc(19,13,3) + pc(19,15,3)

sh(14,11,3) = 2*pc(20,13,3) - 2*pc(20,15,3) - 3*pc(13,13,3) + 3*pc(13,15,3) - &
&3*pc(15,13,3) + 3*pc(15,15,3)

sh(15,11,3) = 4*pc(16,13,3) - 4*pc(16,15,3) - pc(18,13,3) + pc(18,15,3) - pc(1&
&4,13,3) + pc(14,15,3)

sh(16,11,3) = 4*pc(17,13,3) - 4*pc(17,15,3) - pc(12,13,3) + pc(12,15,3) - pc(1&
&9,13,3) + pc(19,15,3)

sh(1,12,3) = pc(1,18,3) - 3*pc(1,14,3)

sh(2,12,3) = pc(2,18,3) - 3*pc(2,14,3)

sh(3,12,3) = pc(3,18,3) - 3*pc(3,14,3)

sh(4,12,3) = pc(4,18,3) - 3*pc(4,14,3)

sh(5,12,3) = pc(8,18,3) - 3*pc(8,14,3)

sh(6,12,3) = pc(9,18,3) - 3*pc(9,14,3)

sh(7,12,3) = pc(10,18,3) - 3*pc(10,14,3)

sh(8,12,3) = pc(5,18,3) - 3*pc(5,14,3) - pc(6,18,3) + 3*pc(6,14,3)

sh(9,12,3) = 2*pc(7,18,3) - 6*pc(7,14,3) - pc(5,18,3) + 3*pc(5,14,3) - pc(6,18&
&,3) + 3*pc(6,14,3)

sh(10,12,3) = pc(11,18,3) - 3*pc(11,14,3)

sh(11,12,3) = pc(13,18,3) - 3*pc(13,14,3) - pc(15,18,3) + 3*pc(15,14,3)

sh(12,12,3) = pc(18,18,3) - 3*pc(18,14,3) - 3*pc(14,18,3) + 9*pc(14,14,3)

sh(13,12,3) = 3*pc(12,18,3) - 9*pc(12,14,3) - pc(19,18,3) + 3*pc(19,14,3)

sh(14,12,3) = 2*pc(20,18,3) - 6*pc(20,14,3) - 3*pc(13,18,3) + 9*pc(13,14,3) - &
&3*pc(15,18,3) + 9*pc(15,14,3)

sh(15,12,3) = 4*pc(16,18,3) - 12*pc(16,14,3) - pc(18,18,3) + 3*pc(18,14,3) - p&
&c(14,18,3) + 3*pc(14,14,3)

sh(16,12,3) = 4*pc(17,18,3) - 12*pc(17,14,3) - pc(12,18,3) + 3*pc(12,14,3) - p&
&c(19,18,3) + 3*pc(19,14,3)

sh(1,13,3) = 3*pc(1,12,3) - pc(1,19,3)

sh(2,13,3) = 3*pc(2,12,3) - pc(2,19,3)

sh(3,13,3) = 3*pc(3,12,3) - pc(3,19,3)

sh(4,13,3) = 3*pc(4,12,3) - pc(4,19,3)

sh(5,13,3) = 3*pc(8,12,3) - pc(8,19,3)

sh(6,13,3) = 3*pc(9,12,3) - pc(9,19,3)

sh(7,13,3) = 3*pc(10,12,3) - pc(10,19,3)

sh(8,13,3) = 3*pc(5,12,3) - pc(5,19,3) - 3*pc(6,12,3) + pc(6,19,3)

sh(9,13,3) = 6*pc(7,12,3) - 2*pc(7,19,3) - 3*pc(5,12,3) + pc(5,19,3) - 3*pc(6,&
&12,3) + pc(6,19,3)

sh(10,13,3) = 3*pc(11,12,3) - pc(11,19,3)

sh(11,13,3) = 3*pc(13,12,3) - pc(13,19,3) - 3*pc(15,12,3) + pc(15,19,3)

sh(12,13,3) = 3*pc(18,12,3) - pc(18,19,3) - 9*pc(14,12,3) + 3*pc(14,19,3)

sh(13,13,3) = 9*pc(12,12,3) - 3*pc(12,19,3) - 3*pc(19,12,3) + pc(19,19,3)

sh(14,13,3) = 6*pc(20,12,3) - 2*pc(20,19,3) - 9*pc(13,12,3) + 3*pc(13,19,3) - &
&9*pc(15,12,3) + 3*pc(15,19,3)

sh(15,13,3) = 12*pc(16,12,3) - 4*pc(16,19,3) - 3*pc(18,12,3) + pc(18,19,3) - 3&
&*pc(14,12,3) + pc(14,19,3)

sh(16,13,3) = 12*pc(17,12,3) - 4*pc(17,19,3) - 3*pc(12,12,3) + pc(12,19,3) - 3&
&*pc(19,12,3) + pc(19,19,3)

sh(1,14,3) = 2*pc(1,20,3) - 3*pc(1,13,3) - 3*pc(1,15,3)

sh(2,14,3) = 2*pc(2,20,3) - 3*pc(2,13,3) - 3*pc(2,15,3)

sh(3,14,3) = 2*pc(3,20,3) - 3*pc(3,13,3) - 3*pc(3,15,3)

sh(4,14,3) = 2*pc(4,20,3) - 3*pc(4,13,3) - 3*pc(4,15,3)

sh(5,14,3) = 2*pc(8,20,3) - 3*pc(8,13,3) - 3*pc(8,15,3)

sh(6,14,3) = 2*pc(9,20,3) - 3*pc(9,13,3) - 3*pc(9,15,3)

sh(7,14,3) = 2*pc(10,20,3) - 3*pc(10,13,3) - 3*pc(10,15,3)

sh(8,14,3) = 2*pc(5,20,3) - 3*pc(5,13,3) - 3*pc(5,15,3) - 2*pc(6,20,3) + 3*pc(&
&6,13,3) + 3*pc(6,15,3)

sh(9,14,3) = 4*pc(7,20,3) - 6*pc(7,13,3) - 6*pc(7,15,3) - 2*pc(5,20,3) + 3*pc(&
&5,13,3) + 3*pc(5,15,3) - 2*pc(6,20,3) + 3*pc(6,13,3) + 3*pc(6,15,3)

sh(10,14,3) = 2*pc(11,20,3) - 3*pc(11,13,3) - 3*pc(11,15,3)

sh(11,14,3) = 2*pc(13,20,3) - 3*pc(13,13,3) - 3*pc(13,15,3) - 2*pc(15,20,3) + &
&3*pc(15,13,3) + 3*pc(15,15,3)

sh(12,14,3) = 2*pc(18,20,3) - 3*pc(18,13,3) - 3*pc(18,15,3) - 6*pc(14,20,3) + &
&9*pc(14,13,3) + 9*pc(14,15,3)

sh(13,14,3) = 6*pc(12,20,3) - 9*pc(12,13,3) - 9*pc(12,15,3) - 2*pc(19,20,3) + &
&3*pc(19,13,3) + 3*pc(19,15,3)

sh(14,14,3) = 4*pc(20,20,3) - 6*pc(20,13,3) - 6*pc(20,15,3) - 6*pc(13,20,3) + &
&9*pc(13,13,3) + 9*pc(13,15,3) - 6*pc(15,20,3) + 9*pc(15,13,3) + 9*pc(15,15,3)

sh(15,14,3) = 8*pc(16,20,3) - 12*pc(16,13,3) - 12*pc(16,15,3) - 2*pc(18,20,3) &
&+ 3*pc(18,13,3) + 3*pc(18,15,3) - 2*pc(14,20,3) + 3*pc(14,13,3) + 3*pc(14,15,3&
&)

sh(16,14,3) = 8*pc(17,20,3) - 12*pc(17,13,3) - 12*pc(17,15,3) - 2*pc(12,20,3) &
&+ 3*pc(12,13,3) + 3*pc(12,15,3) - 2*pc(19,20,3) + 3*pc(19,13,3) + 3*pc(19,15,3&
&)

sh(1,15,3) = 4*pc(1,16,3) - pc(1,18,3) - pc(1,14,3)

sh(2,15,3) = 4*pc(2,16,3) - pc(2,18,3) - pc(2,14,3)

sh(3,15,3) = 4*pc(3,16,3) - pc(3,18,3) - pc(3,14,3)

sh(4,15,3) = 4*pc(4,16,3) - pc(4,18,3) - pc(4,14,3)

sh(5,15,3) = 4*pc(8,16,3) - pc(8,18,3) - pc(8,14,3)

sh(6,15,3) = 4*pc(9,16,3) - pc(9,18,3) - pc(9,14,3)

sh(7,15,3) = 4*pc(10,16,3) - pc(10,18,3) - pc(10,14,3)

sh(8,15,3) = 4*pc(5,16,3) - pc(5,18,3) - pc(5,14,3) - 4*pc(6,16,3) + pc(6,18,3&
&) + pc(6,14,3)

sh(9,15,3) = 8*pc(7,16,3) - 2*pc(7,18,3) - 2*pc(7,14,3) - 4*pc(5,16,3) + pc(5,&
&18,3) + pc(5,14,3) - 4*pc(6,16,3) + pc(6,18,3) + pc(6,14,3)

sh(10,15,3) = 4*pc(11,16,3) - pc(11,18,3) - pc(11,14,3)

sh(11,15,3) = 4*pc(13,16,3) - pc(13,18,3) - pc(13,14,3) - 4*pc(15,16,3) + pc(1&
&5,18,3) + pc(15,14,3)

sh(12,15,3) = 4*pc(18,16,3) - pc(18,18,3) - pc(18,14,3) - 12*pc(14,16,3) + 3*p&
&c(14,18,3) + 3*pc(14,14,3)

sh(13,15,3) = 12*pc(12,16,3) - 3*pc(12,18,3) - 3*pc(12,14,3) - 4*pc(19,16,3) +&
& pc(19,18,3) + pc(19,14,3)

sh(14,15,3) = 8*pc(20,16,3) - 2*pc(20,18,3) - 2*pc(20,14,3) - 12*pc(13,16,3) +&
& 3*pc(13,18,3) + 3*pc(13,14,3) - 12*pc(15,16,3) + 3*pc(15,18,3) + 3*pc(15,14,3&
&)

sh(15,15,3) = 16*pc(16,16,3) - 4*pc(16,18,3) - 4*pc(16,14,3) - 4*pc(18,16,3) +&
& pc(18,18,3) + pc(18,14,3) - 4*pc(14,16,3) + pc(14,18,3) + pc(14,14,3)

sh(16,15,3) = 16*pc(17,16,3) - 4*pc(17,18,3) - 4*pc(17,14,3) - 4*pc(12,16,3) +&
& pc(12,18,3) + pc(12,14,3) - 4*pc(19,16,3) + pc(19,18,3) + pc(19,14,3)

sh(1,16,3) = 4*pc(1,17,3) - pc(1,12,3) - pc(1,19,3)

sh(2,16,3) = 4*pc(2,17,3) - pc(2,12,3) - pc(2,19,3)

sh(3,16,3) = 4*pc(3,17,3) - pc(3,12,3) - pc(3,19,3)

sh(4,16,3) = 4*pc(4,17,3) - pc(4,12,3) - pc(4,19,3)

sh(5,16,3) = 4*pc(8,17,3) - pc(8,12,3) - pc(8,19,3)

sh(6,16,3) = 4*pc(9,17,3) - pc(9,12,3) - pc(9,19,3)

sh(7,16,3) = 4*pc(10,17,3) - pc(10,12,3) - pc(10,19,3)

sh(8,16,3) = 4*pc(5,17,3) - pc(5,12,3) - pc(5,19,3) - 4*pc(6,17,3) + pc(6,12,3&
&) + pc(6,19,3)

sh(9,16,3) = 8*pc(7,17,3) - 2*pc(7,12,3) - 2*pc(7,19,3) - 4*pc(5,17,3) + pc(5,&
&12,3) + pc(5,19,3) - 4*pc(6,17,3) + pc(6,12,3) + pc(6,19,3)

sh(10,16,3) = 4*pc(11,17,3) - pc(11,12,3) - pc(11,19,3)

sh(11,16,3) = 4*pc(13,17,3) - pc(13,12,3) - pc(13,19,3) - 4*pc(15,17,3) + pc(1&
&5,12,3) + pc(15,19,3)

sh(12,16,3) = 4*pc(18,17,3) - pc(18,12,3) - pc(18,19,3) - 12*pc(14,17,3) + 3*p&
&c(14,12,3) + 3*pc(14,19,3)

sh(13,16,3) = 12*pc(12,17,3) - 3*pc(12,12,3) - 3*pc(12,19,3) - 4*pc(19,17,3) +&
& pc(19,12,3) + pc(19,19,3)

sh(14,16,3) = 8*pc(20,17,3) - 2*pc(20,12,3) - 2*pc(20,19,3) - 12*pc(13,17,3) +&
& 3*pc(13,12,3) + 3*pc(13,19,3) - 12*pc(15,17,3) + 3*pc(15,12,3) + 3*pc(15,19,3&
&)

sh(15,16,3) = 16*pc(16,17,3) - 4*pc(16,12,3) - 4*pc(16,19,3) - 4*pc(18,17,3) +&
& pc(18,12,3) + pc(18,19,3) - 4*pc(14,17,3) + pc(14,12,3) + pc(14,19,3)

sh(16,16,3) = 16*pc(17,17,3) - 4*pc(17,12,3) - 4*pc(17,19,3) - 4*pc(12,17,3) +&
& pc(12,12,3) + pc(12,19,3) - 4*pc(19,17,3) + pc(19,12,3) + pc(19,19,3)


   end subroutine momentum2CIntgNum


   function noPrimeMM(step_size, curr_pos, A, B, a1, a2, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos, A, B, a1, a2
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: noPrimeMM

      ! Compute the integral part.
      noPrimeMM = step_size * (curr_pos - A)**l1 * (curr_pos - B)**l2 &
            & * exp(-a1*(curr_pos - A)**2) * exp(-a2*(curr_pos - B)**2)

      return

   end function noPrimeMM


   function primeMM(step_size, curr_pos, A, B, a1, a2, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos, A, B, a1, a2
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: primeMM

      ! Compute each "internal" term of the prime integral. Compare each of
      !   these lines with the 1D momentum matrix equation produce by the
      !   osrecurintg_makenum.py script (appropriately separated into terms).
      !   and the expression in equation 69 in Ben Walker's dissertation.

      ! Second line of equation 69. Also visible as second term in the
      !   script-produced equation.
      primeMM = - 2.0d0*a2 * (curr_pos - B)**(l2+1)

      ! First line in the equation produced by sympy and first line in
      !   equation 69 in Ben Walker's dissertation. As mentioned above,
      !   we only compute this if the angular momentum will allow it.
      if (l2 >= 1) then
         primeMM = primeMM + l2 * (curr_pos - B)**(l2-1)
      endif

      ! Multiply prime integral by the preceeding primitive gaussian
      !   coefficient and exponential and multiply by the succeeding
      !   exponential. (We have already multiplied by the succeeding
      !   primitive gaussian coefficient in the above lines.)
      primeMM = primeMM * (curr_pos-A)**l1 * exp(-a1*(curr_pos-A)**2) &
            & * exp(-a2*(curr_pos-B)**2)

      ! Finally, multiply by the step size.
      primeMM = primeMM * step_size

      return
      
   end function primeMM

end program GaussianIntegrals
