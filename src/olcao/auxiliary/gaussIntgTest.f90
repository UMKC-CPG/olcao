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
         !   using analytical formulas.
         call massvel2CIntgAna(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,1,1),sh(:,:,1,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call massvel2CIntgNum(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & pc(:,:,1,2),sh(:,:,1,2),cell_size,step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,2,alphas,pos,pc(:,:,1,:),sh(:,:,1,:),&
               & "massvel.dat")

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

   subroutine massvel2CIntgAna(a1,a2,A,B,pc,sh)

   use O_Kinds
   use O_Constants, only: pi, fineStructure

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
   real (kind=double), dimension (20,20) :: pc_ke
   real (kind=double), dimension (3) :: P, PA, PB, d
   real (kind=double) :: zeta, inv_2zeta, xi
   real (kind=double) :: preFactorOL, preFactorKE
   real (kind=double) :: preFactor02, preFactor04
   real (kind=double) :: preFactor22, preFactorMV
   real (kind=double) :: inv_2zeta_a, inv_2zeta_b, inv_8m3c2

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
   preFactor02 = sum(PB(:)**2) + 3.0d0/(2.0d0*zeta)
   preFactor04 = sum(PB(:)**4) + sum(PB(:)**2)*3.0d0/zeta &
         & + 9.0d0/(4.0d0*zeta**2)
   preFactor22 = (PB(1)*PB(2))**2 + (PB(1)*PB(3))**2 + (PB(2)*PB(3))**2 &
         & + sum(PB(:)**2)/zeta + 3.0d0/(4.0d0*zeta**2)
   preFactorMV = (fineStructure * 0.001d0)**2 / 8.0d0 &
         & * (16*a2**4*preFactor04 - 80*a2**3*preFactor02 + 60*a2**2 &
         & + 32*a2**4 * preFactor22) * preFactorOL
   inv_8m3c2 = (fineStructure * 0.001d0)**2 / 8.0d0

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

pc_ke(1,1) = preFactorKE

pc_ke(2,1) = PA(1)*preFactorKE + 2*xi*(pc_ol(2,1))

pc_ke(3,1) = PA(2)*preFactorKE + 2*xi*(pc_ol(3,1))

pc_ke(4,1) = PA(3)*preFactorKE + 2*xi*(pc_ol(4,1))

pc_ke(5,1) = PA(1)*pc_ke(2,1) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(5,1) - i&
&nv_2zeta_a*preFactorOL)

pc_ke(6,1) = PA(2)*pc_ke(3,1) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(6,1) - i&
&nv_2zeta_a*preFactorOL)

pc_ke(7,1) = PA(3)*pc_ke(4,1) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(7,1) - i&
&nv_2zeta_a*preFactorOL)

pc_ke(8,1) = PA(2)*pc_ke(2,1) + 2*xi*(pc_ol(8,1))

pc_ke(9,1) = PA(3)*pc_ke(2,1) + 2*xi*(pc_ol(9,1))

pc_ke(10,1) = PA(3)*pc_ke(3,1) + 2*xi*(pc_ol(10,1))

pc_ke(11,1) = PA(3)*pc_ke(8,1) + 2*xi*(pc_ol(11,1))

pc_ke(12,1) = PA(2)*pc_ke(5,1) + 2*xi*(pc_ol(12,1))

pc_ke(13,1) = PA(3)*pc_ke(5,1) + 2*xi*(pc_ol(13,1))

pc_ke(14,1) = PA(1)*pc_ke(6,1) + 2*xi*(pc_ol(14,1))

pc_ke(15,1) = PA(3)*pc_ke(6,1) + 2*xi*(pc_ol(15,1))

pc_ke(16,1) = PA(1)*pc_ke(7,1) + 2*xi*(pc_ol(16,1))

pc_ke(17,1) = PA(2)*pc_ke(7,1) + 2*xi*(pc_ol(17,1))

pc_ke(18,1) = PA(1)*pc_ke(5,1) + inv_2zeta*(2*pc_ke(2,1)) + 2*xi*(pc_ol(18,1) &
&- inv_2zeta_a*2*pc_ol(2,1))

pc_ke(19,1) = PA(2)*pc_ke(6,1) + inv_2zeta*(2*pc_ke(3,1)) + 2*xi*(pc_ol(19,1) &
&- inv_2zeta_a*2*pc_ol(3,1))

pc_ke(20,1) = PA(3)*pc_ke(7,1) + inv_2zeta*(2*pc_ke(4,1)) + 2*xi*(pc_ol(20,1) &
&- inv_2zeta_a*2*pc_ol(4,1))

pc_ke(1,2) = PB(1)*preFactorKE + 2*xi*(pc_ol(1,2))

pc_ke(2,2) = PA(1)*pc_ke(1,2) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(2,2))

pc_ke(3,2) = PA(2)*pc_ke(1,2) + 2*xi*(pc_ol(3,2))

pc_ke(4,2) = PA(3)*pc_ke(1,2) + 2*xi*(pc_ol(4,2))

pc_ke(5,2) = PA(1)*pc_ke(2,2) + inv_2zeta*(pc_ke(1,2) + pc_ke(2,1)) + 2*xi*(pc&
&_ol(5,2) - inv_2zeta_a*pc_ol(1,2))

pc_ke(6,2) = PA(2)*pc_ke(3,2) + inv_2zeta*(pc_ke(1,2)) + 2*xi*(pc_ol(6,2) - in&
&v_2zeta_a*pc_ol(1,2))

pc_ke(7,2) = PA(3)*pc_ke(4,2) + inv_2zeta*(pc_ke(1,2)) + 2*xi*(pc_ol(7,2) - in&
&v_2zeta_a*pc_ol(1,2))

pc_ke(8,2) = PA(2)*pc_ke(2,2) + 2*xi*(pc_ol(8,2))

pc_ke(9,2) = PA(3)*pc_ke(2,2) + 2*xi*(pc_ol(9,2))

pc_ke(10,2) = PA(3)*pc_ke(3,2) + 2*xi*(pc_ol(10,2))

pc_ke(11,2) = PA(3)*pc_ke(8,2) + 2*xi*(pc_ol(11,2))

pc_ke(12,2) = PA(2)*pc_ke(5,2) + 2*xi*(pc_ol(12,2))

pc_ke(13,2) = PA(3)*pc_ke(5,2) + 2*xi*(pc_ol(13,2))

pc_ke(14,2) = PA(1)*pc_ke(6,2) + inv_2zeta*(pc_ke(6,1)) + 2*xi*(pc_ol(14,2))

pc_ke(15,2) = PA(3)*pc_ke(6,2) + 2*xi*(pc_ol(15,2))

pc_ke(16,2) = PA(1)*pc_ke(7,2) + inv_2zeta*(pc_ke(7,1)) + 2*xi*(pc_ol(16,2))

pc_ke(17,2) = PA(2)*pc_ke(7,2) + 2*xi*(pc_ol(17,2))

pc_ke(18,2) = PA(1)*pc_ke(5,2) + inv_2zeta*(2*pc_ke(2,2) + pc_ke(5,1)) + 2*xi*&
&(pc_ol(18,2) - inv_2zeta_a*2*pc_ol(2,2))

pc_ke(19,2) = PA(2)*pc_ke(6,2) + inv_2zeta*(2*pc_ke(3,2)) + 2*xi*(pc_ol(19,2) &
&- inv_2zeta_a*2*pc_ol(3,2))

pc_ke(20,2) = PA(3)*pc_ke(7,2) + inv_2zeta*(2*pc_ke(4,2)) + 2*xi*(pc_ol(20,2) &
&- inv_2zeta_a*2*pc_ol(4,2))

pc_ke(1,3) = PB(2)*preFactorKE + 2*xi*(pc_ol(1,3))

pc_ke(2,3) = PA(1)*pc_ke(1,3) + 2*xi*(pc_ol(2,3))

pc_ke(3,3) = PA(2)*pc_ke(1,3) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(3,3))

pc_ke(4,3) = PA(3)*pc_ke(1,3) + 2*xi*(pc_ol(4,3))

pc_ke(5,3) = PA(1)*pc_ke(2,3) + inv_2zeta*(pc_ke(1,3)) + 2*xi*(pc_ol(5,3) - in&
&v_2zeta_a*pc_ol(1,3))

pc_ke(6,3) = PA(2)*pc_ke(3,3) + inv_2zeta*(pc_ke(1,3) + pc_ke(3,1)) + 2*xi*(pc&
&_ol(6,3) - inv_2zeta_a*pc_ol(1,3))

pc_ke(7,3) = PA(3)*pc_ke(4,3) + inv_2zeta*(pc_ke(1,3)) + 2*xi*(pc_ol(7,3) - in&
&v_2zeta_a*pc_ol(1,3))

pc_ke(8,3) = PA(2)*pc_ke(2,3) + inv_2zeta*(pc_ke(2,1)) + 2*xi*(pc_ol(8,3))

pc_ke(9,3) = PA(3)*pc_ke(2,3) + 2*xi*(pc_ol(9,3))

pc_ke(10,3) = PA(3)*pc_ke(3,3) + 2*xi*(pc_ol(10,3))

pc_ke(11,3) = PA(3)*pc_ke(8,3) + 2*xi*(pc_ol(11,3))

pc_ke(12,3) = PA(2)*pc_ke(5,3) + inv_2zeta*(pc_ke(5,1)) + 2*xi*(pc_ol(12,3))

pc_ke(13,3) = PA(3)*pc_ke(5,3) + 2*xi*(pc_ol(13,3))

pc_ke(14,3) = PA(1)*pc_ke(6,3) + 2*xi*(pc_ol(14,3))

pc_ke(15,3) = PA(3)*pc_ke(6,3) + 2*xi*(pc_ol(15,3))

pc_ke(16,3) = PA(1)*pc_ke(7,3) + 2*xi*(pc_ol(16,3))

pc_ke(17,3) = PA(2)*pc_ke(7,3) + inv_2zeta*(pc_ke(7,1)) + 2*xi*(pc_ol(17,3))

pc_ke(18,3) = PA(1)*pc_ke(5,3) + inv_2zeta*(2*pc_ke(2,3)) + 2*xi*(pc_ol(18,3) &
&- inv_2zeta_a*2*pc_ol(2,3))

pc_ke(19,3) = PA(2)*pc_ke(6,3) + inv_2zeta*(2*pc_ke(3,3) + pc_ke(6,1)) + 2*xi*&
&(pc_ol(19,3) - inv_2zeta_a*2*pc_ol(3,3))

pc_ke(20,3) = PA(3)*pc_ke(7,3) + inv_2zeta*(2*pc_ke(4,3)) + 2*xi*(pc_ol(20,3) &
&- inv_2zeta_a*2*pc_ol(4,3))

pc_ke(1,4) = PB(3)*preFactorKE + 2*xi*(pc_ol(1,4))

pc_ke(2,4) = PA(1)*pc_ke(1,4) + 2*xi*(pc_ol(2,4))

pc_ke(3,4) = PA(2)*pc_ke(1,4) + 2*xi*(pc_ol(3,4))

pc_ke(4,4) = PA(3)*pc_ke(1,4) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(4,4))

pc_ke(5,4) = PA(1)*pc_ke(2,4) + inv_2zeta*(pc_ke(1,4)) + 2*xi*(pc_ol(5,4) - in&
&v_2zeta_a*pc_ol(1,4))

pc_ke(6,4) = PA(2)*pc_ke(3,4) + inv_2zeta*(pc_ke(1,4)) + 2*xi*(pc_ol(6,4) - in&
&v_2zeta_a*pc_ol(1,4))

pc_ke(7,4) = PA(3)*pc_ke(4,4) + inv_2zeta*(pc_ke(1,4) + pc_ke(4,1)) + 2*xi*(pc&
&_ol(7,4) - inv_2zeta_a*pc_ol(1,4))

pc_ke(8,4) = PA(2)*pc_ke(2,4) + 2*xi*(pc_ol(8,4))

pc_ke(9,4) = PA(3)*pc_ke(2,4) + inv_2zeta*(pc_ke(2,1)) + 2*xi*(pc_ol(9,4))

pc_ke(10,4) = PA(3)*pc_ke(3,4) + inv_2zeta*(pc_ke(3,1)) + 2*xi*(pc_ol(10,4))

pc_ke(11,4) = PA(3)*pc_ke(8,4) + inv_2zeta*(pc_ke(8,1)) + 2*xi*(pc_ol(11,4))

pc_ke(12,4) = PA(2)*pc_ke(5,4) + 2*xi*(pc_ol(12,4))

pc_ke(13,4) = PA(3)*pc_ke(5,4) + inv_2zeta*(pc_ke(5,1)) + 2*xi*(pc_ol(13,4))

pc_ke(14,4) = PA(1)*pc_ke(6,4) + 2*xi*(pc_ol(14,4))

pc_ke(15,4) = PA(3)*pc_ke(6,4) + inv_2zeta*(pc_ke(6,1)) + 2*xi*(pc_ol(15,4))

pc_ke(16,4) = PA(1)*pc_ke(7,4) + 2*xi*(pc_ol(16,4))

pc_ke(17,4) = PA(2)*pc_ke(7,4) + 2*xi*(pc_ol(17,4))

pc_ke(18,4) = PA(1)*pc_ke(5,4) + inv_2zeta*(2*pc_ke(2,4)) + 2*xi*(pc_ol(18,4) &
&- inv_2zeta_a*2*pc_ol(2,4))

pc_ke(19,4) = PA(2)*pc_ke(6,4) + inv_2zeta*(2*pc_ke(3,4)) + 2*xi*(pc_ol(19,4) &
&- inv_2zeta_a*2*pc_ol(3,4))

pc_ke(20,4) = PA(3)*pc_ke(7,4) + inv_2zeta*(2*pc_ke(4,4) + pc_ke(7,1)) + 2*xi*&
&(pc_ol(20,4) - inv_2zeta_a*2*pc_ol(4,4))

pc_ke(1,5) = PB(1)*pc_ke(1,2) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(1,5) - i&
&nv_2zeta_b*preFactorOL)

pc_ke(2,5) = PA(1)*pc_ke(1,5) + inv_2zeta*(2*pc_ke(1,2)) + 2*xi*(pc_ol(2,5))

pc_ke(3,5) = PA(2)*pc_ke(1,5) + 2*xi*(pc_ol(3,5))

pc_ke(4,5) = PA(3)*pc_ke(1,5) + 2*xi*(pc_ol(4,5))

pc_ke(5,5) = PA(1)*pc_ke(2,5) + inv_2zeta*(pc_ke(1,5) + 2*pc_ke(2,2)) + 2*xi*(&
&pc_ol(5,5) - inv_2zeta_a*pc_ol(1,5))

pc_ke(6,5) = PA(2)*pc_ke(3,5) + inv_2zeta*(pc_ke(1,5)) + 2*xi*(pc_ol(6,5) - in&
&v_2zeta_a*pc_ol(1,5))

pc_ke(7,5) = PA(3)*pc_ke(4,5) + inv_2zeta*(pc_ke(1,5)) + 2*xi*(pc_ol(7,5) - in&
&v_2zeta_a*pc_ol(1,5))

pc_ke(8,5) = PA(2)*pc_ke(2,5) + 2*xi*(pc_ol(8,5))

pc_ke(9,5) = PA(3)*pc_ke(2,5) + 2*xi*(pc_ol(9,5))

pc_ke(10,5) = PA(3)*pc_ke(3,5) + 2*xi*(pc_ol(10,5))

pc_ke(11,5) = PA(3)*pc_ke(8,5) + 2*xi*(pc_ol(11,5))

pc_ke(12,5) = PA(2)*pc_ke(5,5) + 2*xi*(pc_ol(12,5))

pc_ke(13,5) = PA(3)*pc_ke(5,5) + 2*xi*(pc_ol(13,5))

pc_ke(14,5) = PA(1)*pc_ke(6,5) + inv_2zeta*(2*pc_ke(6,2)) + 2*xi*(pc_ol(14,5))

pc_ke(15,5) = PA(3)*pc_ke(6,5) + 2*xi*(pc_ol(15,5))

pc_ke(16,5) = PA(1)*pc_ke(7,5) + inv_2zeta*(2*pc_ke(7,2)) + 2*xi*(pc_ol(16,5))

pc_ke(17,5) = PA(2)*pc_ke(7,5) + 2*xi*(pc_ol(17,5))

pc_ke(18,5) = PA(1)*pc_ke(5,5) + inv_2zeta*(2*pc_ke(2,5) + 2*pc_ke(5,2)) + 2*x&
&i*(pc_ol(18,5) - inv_2zeta_a*2*pc_ol(2,5))

pc_ke(19,5) = PA(2)*pc_ke(6,5) + inv_2zeta*(2*pc_ke(3,5)) + 2*xi*(pc_ol(19,5) &
&- inv_2zeta_a*2*pc_ol(3,5))

pc_ke(20,5) = PA(3)*pc_ke(7,5) + inv_2zeta*(2*pc_ke(4,5)) + 2*xi*(pc_ol(20,5) &
&- inv_2zeta_a*2*pc_ol(4,5))

pc_ke(1,6) = PB(2)*pc_ke(1,3) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(1,6) - i&
&nv_2zeta_b*preFactorOL)

pc_ke(2,6) = PA(1)*pc_ke(1,6) + 2*xi*(pc_ol(2,6))

pc_ke(3,6) = PA(2)*pc_ke(1,6) + inv_2zeta*(2*pc_ke(1,3)) + 2*xi*(pc_ol(3,6))

pc_ke(4,6) = PA(3)*pc_ke(1,6) + 2*xi*(pc_ol(4,6))

pc_ke(5,6) = PA(1)*pc_ke(2,6) + inv_2zeta*(pc_ke(1,6)) + 2*xi*(pc_ol(5,6) - in&
&v_2zeta_a*pc_ol(1,6))

pc_ke(6,6) = PA(2)*pc_ke(3,6) + inv_2zeta*(pc_ke(1,6) + 2*pc_ke(3,3)) + 2*xi*(&
&pc_ol(6,6) - inv_2zeta_a*pc_ol(1,6))

pc_ke(7,6) = PA(3)*pc_ke(4,6) + inv_2zeta*(pc_ke(1,6)) + 2*xi*(pc_ol(7,6) - in&
&v_2zeta_a*pc_ol(1,6))

pc_ke(8,6) = PA(2)*pc_ke(2,6) + inv_2zeta*(2*pc_ke(2,3)) + 2*xi*(pc_ol(8,6))

pc_ke(9,6) = PA(3)*pc_ke(2,6) + 2*xi*(pc_ol(9,6))

pc_ke(10,6) = PA(3)*pc_ke(3,6) + 2*xi*(pc_ol(10,6))

pc_ke(11,6) = PA(3)*pc_ke(8,6) + 2*xi*(pc_ol(11,6))

pc_ke(12,6) = PA(2)*pc_ke(5,6) + inv_2zeta*(2*pc_ke(5,3)) + 2*xi*(pc_ol(12,6))

pc_ke(13,6) = PA(3)*pc_ke(5,6) + 2*xi*(pc_ol(13,6))

pc_ke(14,6) = PA(1)*pc_ke(6,6) + 2*xi*(pc_ol(14,6))

pc_ke(15,6) = PA(3)*pc_ke(6,6) + 2*xi*(pc_ol(15,6))

pc_ke(16,6) = PA(1)*pc_ke(7,6) + 2*xi*(pc_ol(16,6))

pc_ke(17,6) = PA(2)*pc_ke(7,6) + inv_2zeta*(2*pc_ke(7,3)) + 2*xi*(pc_ol(17,6))

pc_ke(18,6) = PA(1)*pc_ke(5,6) + inv_2zeta*(2*pc_ke(2,6)) + 2*xi*(pc_ol(18,6) &
&- inv_2zeta_a*2*pc_ol(2,6))

pc_ke(19,6) = PA(2)*pc_ke(6,6) + inv_2zeta*(2*pc_ke(3,6) + 2*pc_ke(6,3)) + 2*x&
&i*(pc_ol(19,6) - inv_2zeta_a*2*pc_ol(3,6))

pc_ke(20,6) = PA(3)*pc_ke(7,6) + inv_2zeta*(2*pc_ke(4,6)) + 2*xi*(pc_ol(20,6) &
&- inv_2zeta_a*2*pc_ol(4,6))

pc_ke(1,7) = PB(3)*pc_ke(1,4) + inv_2zeta*(preFactorKE) + 2*xi*(pc_ol(1,7) - i&
&nv_2zeta_b*preFactorOL)

pc_ke(2,7) = PA(1)*pc_ke(1,7) + 2*xi*(pc_ol(2,7))

pc_ke(3,7) = PA(2)*pc_ke(1,7) + 2*xi*(pc_ol(3,7))

pc_ke(4,7) = PA(3)*pc_ke(1,7) + inv_2zeta*(2*pc_ke(1,4)) + 2*xi*(pc_ol(4,7))

pc_ke(5,7) = PA(1)*pc_ke(2,7) + inv_2zeta*(pc_ke(1,7)) + 2*xi*(pc_ol(5,7) - in&
&v_2zeta_a*pc_ol(1,7))

pc_ke(6,7) = PA(2)*pc_ke(3,7) + inv_2zeta*(pc_ke(1,7)) + 2*xi*(pc_ol(6,7) - in&
&v_2zeta_a*pc_ol(1,7))

pc_ke(7,7) = PA(3)*pc_ke(4,7) + inv_2zeta*(pc_ke(1,7) + 2*pc_ke(4,4)) + 2*xi*(&
&pc_ol(7,7) - inv_2zeta_a*pc_ol(1,7))

pc_ke(8,7) = PA(2)*pc_ke(2,7) + 2*xi*(pc_ol(8,7))

pc_ke(9,7) = PA(3)*pc_ke(2,7) + inv_2zeta*(2*pc_ke(2,4)) + 2*xi*(pc_ol(9,7))

pc_ke(10,7) = PA(3)*pc_ke(3,7) + inv_2zeta*(2*pc_ke(3,4)) + 2*xi*(pc_ol(10,7))

pc_ke(11,7) = PA(3)*pc_ke(8,7) + inv_2zeta*(2*pc_ke(8,4)) + 2*xi*(pc_ol(11,7))

pc_ke(12,7) = PA(2)*pc_ke(5,7) + 2*xi*(pc_ol(12,7))

pc_ke(13,7) = PA(3)*pc_ke(5,7) + inv_2zeta*(2*pc_ke(5,4)) + 2*xi*(pc_ol(13,7))

pc_ke(14,7) = PA(1)*pc_ke(6,7) + 2*xi*(pc_ol(14,7))

pc_ke(15,7) = PA(3)*pc_ke(6,7) + inv_2zeta*(2*pc_ke(6,4)) + 2*xi*(pc_ol(15,7))

pc_ke(16,7) = PA(1)*pc_ke(7,7) + 2*xi*(pc_ol(16,7))

pc_ke(17,7) = PA(2)*pc_ke(7,7) + 2*xi*(pc_ol(17,7))

pc_ke(18,7) = PA(1)*pc_ke(5,7) + inv_2zeta*(2*pc_ke(2,7)) + 2*xi*(pc_ol(18,7) &
&- inv_2zeta_a*2*pc_ol(2,7))

pc_ke(19,7) = PA(2)*pc_ke(6,7) + inv_2zeta*(2*pc_ke(3,7)) + 2*xi*(pc_ol(19,7) &
&- inv_2zeta_a*2*pc_ol(3,7))

pc_ke(20,7) = PA(3)*pc_ke(7,7) + inv_2zeta*(2*pc_ke(4,7) + 2*pc_ke(7,4)) + 2*x&
&i*(pc_ol(20,7) - inv_2zeta_a*2*pc_ol(4,7))

pc_ke(1,8) = PB(2)*pc_ke(1,2) + 2*xi*(pc_ol(1,8))

pc_ke(2,8) = PA(1)*pc_ke(1,8) + inv_2zeta*(pc_ke(1,3)) + 2*xi*(pc_ol(2,8))

pc_ke(3,8) = PA(2)*pc_ke(1,8) + inv_2zeta*(pc_ke(1,2)) + 2*xi*(pc_ol(3,8))

pc_ke(4,8) = PA(3)*pc_ke(1,8) + 2*xi*(pc_ol(4,8))

pc_ke(5,8) = PA(1)*pc_ke(2,8) + inv_2zeta*(pc_ke(1,8) + pc_ke(2,3)) + 2*xi*(pc&
&_ol(5,8) - inv_2zeta_a*pc_ol(1,8))

pc_ke(6,8) = PA(2)*pc_ke(3,8) + inv_2zeta*(pc_ke(1,8) + pc_ke(3,2)) + 2*xi*(pc&
&_ol(6,8) - inv_2zeta_a*pc_ol(1,8))

pc_ke(7,8) = PA(3)*pc_ke(4,8) + inv_2zeta*(pc_ke(1,8)) + 2*xi*(pc_ol(7,8) - in&
&v_2zeta_a*pc_ol(1,8))

pc_ke(8,8) = PA(2)*pc_ke(2,8) + inv_2zeta*(pc_ke(2,2)) + 2*xi*(pc_ol(8,8))

pc_ke(9,8) = PA(3)*pc_ke(2,8) + 2*xi*(pc_ol(9,8))

pc_ke(10,8) = PA(3)*pc_ke(3,8) + 2*xi*(pc_ol(10,8))

pc_ke(11,8) = PA(3)*pc_ke(8,8) + 2*xi*(pc_ol(11,8))

pc_ke(12,8) = PA(2)*pc_ke(5,8) + inv_2zeta*(pc_ke(5,2)) + 2*xi*(pc_ol(12,8))

pc_ke(13,8) = PA(3)*pc_ke(5,8) + 2*xi*(pc_ol(13,8))

pc_ke(14,8) = PA(1)*pc_ke(6,8) + inv_2zeta*(pc_ke(6,3)) + 2*xi*(pc_ol(14,8))

pc_ke(15,8) = PA(3)*pc_ke(6,8) + 2*xi*(pc_ol(15,8))

pc_ke(16,8) = PA(1)*pc_ke(7,8) + inv_2zeta*(pc_ke(7,3)) + 2*xi*(pc_ol(16,8))

pc_ke(17,8) = PA(2)*pc_ke(7,8) + inv_2zeta*(pc_ke(7,2)) + 2*xi*(pc_ol(17,8))

pc_ke(18,8) = PA(1)*pc_ke(5,8) + inv_2zeta*(2*pc_ke(2,8) + pc_ke(5,3)) + 2*xi*&
&(pc_ol(18,8) - inv_2zeta_a*2*pc_ol(2,8))

pc_ke(19,8) = PA(2)*pc_ke(6,8) + inv_2zeta*(2*pc_ke(3,8) + pc_ke(6,2)) + 2*xi*&
&(pc_ol(19,8) - inv_2zeta_a*2*pc_ol(3,8))

pc_ke(20,8) = PA(3)*pc_ke(7,8) + inv_2zeta*(2*pc_ke(4,8)) + 2*xi*(pc_ol(20,8) &
&- inv_2zeta_a*2*pc_ol(4,8))

pc_ke(1,9) = PB(3)*pc_ke(1,2) + 2*xi*(pc_ol(1,9))

pc_ke(2,9) = PA(1)*pc_ke(1,9) + inv_2zeta*(pc_ke(1,4)) + 2*xi*(pc_ol(2,9))

pc_ke(3,9) = PA(2)*pc_ke(1,9) + 2*xi*(pc_ol(3,9))

pc_ke(4,9) = PA(3)*pc_ke(1,9) + inv_2zeta*(pc_ke(1,2)) + 2*xi*(pc_ol(4,9))

pc_ke(5,9) = PA(1)*pc_ke(2,9) + inv_2zeta*(pc_ke(1,9) + pc_ke(2,4)) + 2*xi*(pc&
&_ol(5,9) - inv_2zeta_a*pc_ol(1,9))

pc_ke(6,9) = PA(2)*pc_ke(3,9) + inv_2zeta*(pc_ke(1,9)) + 2*xi*(pc_ol(6,9) - in&
&v_2zeta_a*pc_ol(1,9))

pc_ke(7,9) = PA(3)*pc_ke(4,9) + inv_2zeta*(pc_ke(1,9) + pc_ke(4,2)) + 2*xi*(pc&
&_ol(7,9) - inv_2zeta_a*pc_ol(1,9))

pc_ke(8,9) = PA(2)*pc_ke(2,9) + 2*xi*(pc_ol(8,9))

pc_ke(9,9) = PA(3)*pc_ke(2,9) + inv_2zeta*(pc_ke(2,2)) + 2*xi*(pc_ol(9,9))

pc_ke(10,9) = PA(3)*pc_ke(3,9) + inv_2zeta*(pc_ke(3,2)) + 2*xi*(pc_ol(10,9))

pc_ke(11,9) = PA(3)*pc_ke(8,9) + inv_2zeta*(pc_ke(8,2)) + 2*xi*(pc_ol(11,9))

pc_ke(12,9) = PA(2)*pc_ke(5,9) + 2*xi*(pc_ol(12,9))

pc_ke(13,9) = PA(3)*pc_ke(5,9) + inv_2zeta*(pc_ke(5,2)) + 2*xi*(pc_ol(13,9))

pc_ke(14,9) = PA(1)*pc_ke(6,9) + inv_2zeta*(pc_ke(6,4)) + 2*xi*(pc_ol(14,9))

pc_ke(15,9) = PA(3)*pc_ke(6,9) + inv_2zeta*(pc_ke(6,2)) + 2*xi*(pc_ol(15,9))

pc_ke(16,9) = PA(1)*pc_ke(7,9) + inv_2zeta*(pc_ke(7,4)) + 2*xi*(pc_ol(16,9))

pc_ke(17,9) = PA(2)*pc_ke(7,9) + 2*xi*(pc_ol(17,9))

pc_ke(18,9) = PA(1)*pc_ke(5,9) + inv_2zeta*(2*pc_ke(2,9) + pc_ke(5,4)) + 2*xi*&
&(pc_ol(18,9) - inv_2zeta_a*2*pc_ol(2,9))

pc_ke(19,9) = PA(2)*pc_ke(6,9) + inv_2zeta*(2*pc_ke(3,9)) + 2*xi*(pc_ol(19,9) &
&- inv_2zeta_a*2*pc_ol(3,9))

pc_ke(20,9) = PA(3)*pc_ke(7,9) + inv_2zeta*(2*pc_ke(4,9) + pc_ke(7,2)) + 2*xi*&
&(pc_ol(20,9) - inv_2zeta_a*2*pc_ol(4,9))

pc_ke(1,10) = PB(3)*pc_ke(1,3) + 2*xi*(pc_ol(1,10))

pc_ke(2,10) = PA(1)*pc_ke(1,10) + 2*xi*(pc_ol(2,10))

pc_ke(3,10) = PA(2)*pc_ke(1,10) + inv_2zeta*(pc_ke(1,4)) + 2*xi*(pc_ol(3,10))

pc_ke(4,10) = PA(3)*pc_ke(1,10) + inv_2zeta*(pc_ke(1,3)) + 2*xi*(pc_ol(4,10))

pc_ke(5,10) = PA(1)*pc_ke(2,10) + inv_2zeta*(pc_ke(1,10)) + 2*xi*(pc_ol(5,10) &
&- inv_2zeta_a*pc_ol(1,10))

pc_ke(6,10) = PA(2)*pc_ke(3,10) + inv_2zeta*(pc_ke(1,10) + pc_ke(3,4)) + 2*xi*&
&(pc_ol(6,10) - inv_2zeta_a*pc_ol(1,10))

pc_ke(7,10) = PA(3)*pc_ke(4,10) + inv_2zeta*(pc_ke(1,10) + pc_ke(4,3)) + 2*xi*&
&(pc_ol(7,10) - inv_2zeta_a*pc_ol(1,10))

pc_ke(8,10) = PA(2)*pc_ke(2,10) + inv_2zeta*(pc_ke(2,4)) + 2*xi*(pc_ol(8,10))

pc_ke(9,10) = PA(3)*pc_ke(2,10) + inv_2zeta*(pc_ke(2,3)) + 2*xi*(pc_ol(9,10))

pc_ke(10,10) = PA(3)*pc_ke(3,10) + inv_2zeta*(pc_ke(3,3)) + 2*xi*(pc_ol(10,10)&
&)

pc_ke(11,10) = PA(3)*pc_ke(8,10) + inv_2zeta*(pc_ke(8,3)) + 2*xi*(pc_ol(11,10)&
&)

pc_ke(12,10) = PA(2)*pc_ke(5,10) + inv_2zeta*(pc_ke(5,4)) + 2*xi*(pc_ol(12,10)&
&)

pc_ke(13,10) = PA(3)*pc_ke(5,10) + inv_2zeta*(pc_ke(5,3)) + 2*xi*(pc_ol(13,10)&
&)

pc_ke(14,10) = PA(1)*pc_ke(6,10) + 2*xi*(pc_ol(14,10))

pc_ke(15,10) = PA(3)*pc_ke(6,10) + inv_2zeta*(pc_ke(6,3)) + 2*xi*(pc_ol(15,10)&
&)

pc_ke(16,10) = PA(1)*pc_ke(7,10) + 2*xi*(pc_ol(16,10))

pc_ke(17,10) = PA(2)*pc_ke(7,10) + inv_2zeta*(pc_ke(7,4)) + 2*xi*(pc_ol(17,10)&
&)

pc_ke(18,10) = PA(1)*pc_ke(5,10) + inv_2zeta*(2*pc_ke(2,10)) + 2*xi*(pc_ol(18,&
&10) - inv_2zeta_a*2*pc_ol(2,10))

pc_ke(19,10) = PA(2)*pc_ke(6,10) + inv_2zeta*(2*pc_ke(3,10) + pc_ke(6,4)) + 2*&
&xi*(pc_ol(19,10) - inv_2zeta_a*2*pc_ol(3,10))

pc_ke(20,10) = PA(3)*pc_ke(7,10) + inv_2zeta*(2*pc_ke(4,10) + pc_ke(7,3)) + 2*&
&xi*(pc_ol(20,10) - inv_2zeta_a*2*pc_ol(4,10))

pc_ke(1,11) = PB(3)*pc_ke(1,8) + 2*xi*(pc_ol(1,11))

pc_ke(2,11) = PA(1)*pc_ke(1,11) + inv_2zeta*(pc_ke(1,10)) + 2*xi*(pc_ol(2,11))

pc_ke(3,11) = PA(2)*pc_ke(1,11) + inv_2zeta*(pc_ke(1,9)) + 2*xi*(pc_ol(3,11))

pc_ke(4,11) = PA(3)*pc_ke(1,11) + inv_2zeta*(pc_ke(1,8)) + 2*xi*(pc_ol(4,11))

pc_ke(5,11) = PA(1)*pc_ke(2,11) + inv_2zeta*(pc_ke(1,11) + pc_ke(2,10)) + 2*xi&
&*(pc_ol(5,11) - inv_2zeta_a*pc_ol(1,11))

pc_ke(6,11) = PA(2)*pc_ke(3,11) + inv_2zeta*(pc_ke(1,11) + pc_ke(3,9)) + 2*xi*&
&(pc_ol(6,11) - inv_2zeta_a*pc_ol(1,11))

pc_ke(7,11) = PA(3)*pc_ke(4,11) + inv_2zeta*(pc_ke(1,11) + pc_ke(4,8)) + 2*xi*&
&(pc_ol(7,11) - inv_2zeta_a*pc_ol(1,11))

pc_ke(8,11) = PA(2)*pc_ke(2,11) + inv_2zeta*(pc_ke(2,9)) + 2*xi*(pc_ol(8,11))

pc_ke(9,11) = PA(3)*pc_ke(2,11) + inv_2zeta*(pc_ke(2,8)) + 2*xi*(pc_ol(9,11))

pc_ke(10,11) = PA(3)*pc_ke(3,11) + inv_2zeta*(pc_ke(3,8)) + 2*xi*(pc_ol(10,11)&
&)

pc_ke(11,11) = PA(3)*pc_ke(8,11) + inv_2zeta*(pc_ke(8,8)) + 2*xi*(pc_ol(11,11)&
&)

pc_ke(12,11) = PA(2)*pc_ke(5,11) + inv_2zeta*(pc_ke(5,9)) + 2*xi*(pc_ol(12,11)&
&)

pc_ke(13,11) = PA(3)*pc_ke(5,11) + inv_2zeta*(pc_ke(5,8)) + 2*xi*(pc_ol(13,11)&
&)

pc_ke(14,11) = PA(1)*pc_ke(6,11) + inv_2zeta*(pc_ke(6,10)) + 2*xi*(pc_ol(14,11&
&))

pc_ke(15,11) = PA(3)*pc_ke(6,11) + inv_2zeta*(pc_ke(6,8)) + 2*xi*(pc_ol(15,11)&
&)

pc_ke(16,11) = PA(1)*pc_ke(7,11) + inv_2zeta*(pc_ke(7,10)) + 2*xi*(pc_ol(16,11&
&))

pc_ke(17,11) = PA(2)*pc_ke(7,11) + inv_2zeta*(pc_ke(7,9)) + 2*xi*(pc_ol(17,11)&
&)

pc_ke(18,11) = PA(1)*pc_ke(5,11) + inv_2zeta*(2*pc_ke(2,11) + pc_ke(5,10)) + 2&
&*xi*(pc_ol(18,11) - inv_2zeta_a*2*pc_ol(2,11))

pc_ke(19,11) = PA(2)*pc_ke(6,11) + inv_2zeta*(2*pc_ke(3,11) + pc_ke(6,9)) + 2*&
&xi*(pc_ol(19,11) - inv_2zeta_a*2*pc_ol(3,11))

pc_ke(20,11) = PA(3)*pc_ke(7,11) + inv_2zeta*(2*pc_ke(4,11) + pc_ke(7,8)) + 2*&
&xi*(pc_ol(20,11) - inv_2zeta_a*2*pc_ol(4,11))

pc_ke(1,12) = PB(2)*pc_ke(1,5) + 2*xi*(pc_ol(1,12))

pc_ke(2,12) = PA(1)*pc_ke(1,12) + inv_2zeta*(2*pc_ke(1,8)) + 2*xi*(pc_ol(2,12)&
&)

pc_ke(3,12) = PA(2)*pc_ke(1,12) + inv_2zeta*(pc_ke(1,5)) + 2*xi*(pc_ol(3,12))

pc_ke(4,12) = PA(3)*pc_ke(1,12) + 2*xi*(pc_ol(4,12))

pc_ke(5,12) = PA(1)*pc_ke(2,12) + inv_2zeta*(pc_ke(1,12) + 2*pc_ke(2,8)) + 2*x&
&i*(pc_ol(5,12) - inv_2zeta_a*pc_ol(1,12))

pc_ke(6,12) = PA(2)*pc_ke(3,12) + inv_2zeta*(pc_ke(1,12) + pc_ke(3,5)) + 2*xi*&
&(pc_ol(6,12) - inv_2zeta_a*pc_ol(1,12))

pc_ke(7,12) = PA(3)*pc_ke(4,12) + inv_2zeta*(pc_ke(1,12)) + 2*xi*(pc_ol(7,12) &
&- inv_2zeta_a*pc_ol(1,12))

pc_ke(8,12) = PA(2)*pc_ke(2,12) + inv_2zeta*(pc_ke(2,5)) + 2*xi*(pc_ol(8,12))

pc_ke(9,12) = PA(3)*pc_ke(2,12) + 2*xi*(pc_ol(9,12))

pc_ke(10,12) = PA(3)*pc_ke(3,12) + 2*xi*(pc_ol(10,12))

pc_ke(11,12) = PA(3)*pc_ke(8,12) + 2*xi*(pc_ol(11,12))

pc_ke(12,12) = PA(2)*pc_ke(5,12) + inv_2zeta*(pc_ke(5,5)) + 2*xi*(pc_ol(12,12)&
&)

pc_ke(13,12) = PA(3)*pc_ke(5,12) + 2*xi*(pc_ol(13,12))

pc_ke(14,12) = PA(1)*pc_ke(6,12) + inv_2zeta*(2*pc_ke(6,8)) + 2*xi*(pc_ol(14,1&
&2))

pc_ke(15,12) = PA(3)*pc_ke(6,12) + 2*xi*(pc_ol(15,12))

pc_ke(16,12) = PA(1)*pc_ke(7,12) + inv_2zeta*(2*pc_ke(7,8)) + 2*xi*(pc_ol(16,1&
&2))

pc_ke(17,12) = PA(2)*pc_ke(7,12) + inv_2zeta*(pc_ke(7,5)) + 2*xi*(pc_ol(17,12)&
&)

pc_ke(18,12) = PA(1)*pc_ke(5,12) + inv_2zeta*(2*pc_ke(2,12) + 2*pc_ke(5,8)) + &
&2*xi*(pc_ol(18,12) - inv_2zeta_a*2*pc_ol(2,12))

pc_ke(19,12) = PA(2)*pc_ke(6,12) + inv_2zeta*(2*pc_ke(3,12) + pc_ke(6,5)) + 2*&
&xi*(pc_ol(19,12) - inv_2zeta_a*2*pc_ol(3,12))

pc_ke(20,12) = PA(3)*pc_ke(7,12) + inv_2zeta*(2*pc_ke(4,12)) + 2*xi*(pc_ol(20,&
&12) - inv_2zeta_a*2*pc_ol(4,12))

pc_ke(1,13) = PB(3)*pc_ke(1,5) + 2*xi*(pc_ol(1,13))

pc_ke(2,13) = PA(1)*pc_ke(1,13) + inv_2zeta*(2*pc_ke(1,9)) + 2*xi*(pc_ol(2,13)&
&)

pc_ke(3,13) = PA(2)*pc_ke(1,13) + 2*xi*(pc_ol(3,13))

pc_ke(4,13) = PA(3)*pc_ke(1,13) + inv_2zeta*(pc_ke(1,5)) + 2*xi*(pc_ol(4,13))

pc_ke(5,13) = PA(1)*pc_ke(2,13) + inv_2zeta*(pc_ke(1,13) + 2*pc_ke(2,9)) + 2*x&
&i*(pc_ol(5,13) - inv_2zeta_a*pc_ol(1,13))

pc_ke(6,13) = PA(2)*pc_ke(3,13) + inv_2zeta*(pc_ke(1,13)) + 2*xi*(pc_ol(6,13) &
&- inv_2zeta_a*pc_ol(1,13))

pc_ke(7,13) = PA(3)*pc_ke(4,13) + inv_2zeta*(pc_ke(1,13) + pc_ke(4,5)) + 2*xi*&
&(pc_ol(7,13) - inv_2zeta_a*pc_ol(1,13))

pc_ke(8,13) = PA(2)*pc_ke(2,13) + 2*xi*(pc_ol(8,13))

pc_ke(9,13) = PA(3)*pc_ke(2,13) + inv_2zeta*(pc_ke(2,5)) + 2*xi*(pc_ol(9,13))

pc_ke(10,13) = PA(3)*pc_ke(3,13) + inv_2zeta*(pc_ke(3,5)) + 2*xi*(pc_ol(10,13)&
&)

pc_ke(11,13) = PA(3)*pc_ke(8,13) + inv_2zeta*(pc_ke(8,5)) + 2*xi*(pc_ol(11,13)&
&)

pc_ke(12,13) = PA(2)*pc_ke(5,13) + 2*xi*(pc_ol(12,13))

pc_ke(13,13) = PA(3)*pc_ke(5,13) + inv_2zeta*(pc_ke(5,5)) + 2*xi*(pc_ol(13,13)&
&)

pc_ke(14,13) = PA(1)*pc_ke(6,13) + inv_2zeta*(2*pc_ke(6,9)) + 2*xi*(pc_ol(14,1&
&3))

pc_ke(15,13) = PA(3)*pc_ke(6,13) + inv_2zeta*(pc_ke(6,5)) + 2*xi*(pc_ol(15,13)&
&)

pc_ke(16,13) = PA(1)*pc_ke(7,13) + inv_2zeta*(2*pc_ke(7,9)) + 2*xi*(pc_ol(16,1&
&3))

pc_ke(17,13) = PA(2)*pc_ke(7,13) + 2*xi*(pc_ol(17,13))

pc_ke(18,13) = PA(1)*pc_ke(5,13) + inv_2zeta*(2*pc_ke(2,13) + 2*pc_ke(5,9)) + &
&2*xi*(pc_ol(18,13) - inv_2zeta_a*2*pc_ol(2,13))

pc_ke(19,13) = PA(2)*pc_ke(6,13) + inv_2zeta*(2*pc_ke(3,13)) + 2*xi*(pc_ol(19,&
&13) - inv_2zeta_a*2*pc_ol(3,13))

pc_ke(20,13) = PA(3)*pc_ke(7,13) + inv_2zeta*(2*pc_ke(4,13) + pc_ke(7,5)) + 2*&
&xi*(pc_ol(20,13) - inv_2zeta_a*2*pc_ol(4,13))

pc_ke(1,14) = PB(1)*pc_ke(1,6) + 2*xi*(pc_ol(1,14))

pc_ke(2,14) = PA(1)*pc_ke(1,14) + inv_2zeta*(pc_ke(1,6)) + 2*xi*(pc_ol(2,14))

pc_ke(3,14) = PA(2)*pc_ke(1,14) + inv_2zeta*(2*pc_ke(1,8)) + 2*xi*(pc_ol(3,14)&
&)

pc_ke(4,14) = PA(3)*pc_ke(1,14) + 2*xi*(pc_ol(4,14))

pc_ke(5,14) = PA(1)*pc_ke(2,14) + inv_2zeta*(pc_ke(1,14) + pc_ke(2,6)) + 2*xi*&
&(pc_ol(5,14) - inv_2zeta_a*pc_ol(1,14))

pc_ke(6,14) = PA(2)*pc_ke(3,14) + inv_2zeta*(pc_ke(1,14) + 2*pc_ke(3,8)) + 2*x&
&i*(pc_ol(6,14) - inv_2zeta_a*pc_ol(1,14))

pc_ke(7,14) = PA(3)*pc_ke(4,14) + inv_2zeta*(pc_ke(1,14)) + 2*xi*(pc_ol(7,14) &
&- inv_2zeta_a*pc_ol(1,14))

pc_ke(8,14) = PA(2)*pc_ke(2,14) + inv_2zeta*(2*pc_ke(2,8)) + 2*xi*(pc_ol(8,14)&
&)

pc_ke(9,14) = PA(3)*pc_ke(2,14) + 2*xi*(pc_ol(9,14))

pc_ke(10,14) = PA(3)*pc_ke(3,14) + 2*xi*(pc_ol(10,14))

pc_ke(11,14) = PA(3)*pc_ke(8,14) + 2*xi*(pc_ol(11,14))

pc_ke(12,14) = PA(2)*pc_ke(5,14) + inv_2zeta*(2*pc_ke(5,8)) + 2*xi*(pc_ol(12,1&
&4))

pc_ke(13,14) = PA(3)*pc_ke(5,14) + 2*xi*(pc_ol(13,14))

pc_ke(14,14) = PA(1)*pc_ke(6,14) + inv_2zeta*(pc_ke(6,6)) + 2*xi*(pc_ol(14,14)&
&)

pc_ke(15,14) = PA(3)*pc_ke(6,14) + 2*xi*(pc_ol(15,14))

pc_ke(16,14) = PA(1)*pc_ke(7,14) + inv_2zeta*(pc_ke(7,6)) + 2*xi*(pc_ol(16,14)&
&)

pc_ke(17,14) = PA(2)*pc_ke(7,14) + inv_2zeta*(2*pc_ke(7,8)) + 2*xi*(pc_ol(17,1&
&4))

pc_ke(18,14) = PA(1)*pc_ke(5,14) + inv_2zeta*(2*pc_ke(2,14) + pc_ke(5,6)) + 2*&
&xi*(pc_ol(18,14) - inv_2zeta_a*2*pc_ol(2,14))

pc_ke(19,14) = PA(2)*pc_ke(6,14) + inv_2zeta*(2*pc_ke(3,14) + 2*pc_ke(6,8)) + &
&2*xi*(pc_ol(19,14) - inv_2zeta_a*2*pc_ol(3,14))

pc_ke(20,14) = PA(3)*pc_ke(7,14) + inv_2zeta*(2*pc_ke(4,14)) + 2*xi*(pc_ol(20,&
&14) - inv_2zeta_a*2*pc_ol(4,14))

pc_ke(1,15) = PB(3)*pc_ke(1,6) + 2*xi*(pc_ol(1,15))

pc_ke(2,15) = PA(1)*pc_ke(1,15) + 2*xi*(pc_ol(2,15))

pc_ke(3,15) = PA(2)*pc_ke(1,15) + inv_2zeta*(2*pc_ke(1,10)) + 2*xi*(pc_ol(3,15&
&))

pc_ke(4,15) = PA(3)*pc_ke(1,15) + inv_2zeta*(pc_ke(1,6)) + 2*xi*(pc_ol(4,15))

pc_ke(5,15) = PA(1)*pc_ke(2,15) + inv_2zeta*(pc_ke(1,15)) + 2*xi*(pc_ol(5,15) &
&- inv_2zeta_a*pc_ol(1,15))

pc_ke(6,15) = PA(2)*pc_ke(3,15) + inv_2zeta*(pc_ke(1,15) + 2*pc_ke(3,10)) + 2*&
&xi*(pc_ol(6,15) - inv_2zeta_a*pc_ol(1,15))

pc_ke(7,15) = PA(3)*pc_ke(4,15) + inv_2zeta*(pc_ke(1,15) + pc_ke(4,6)) + 2*xi*&
&(pc_ol(7,15) - inv_2zeta_a*pc_ol(1,15))

pc_ke(8,15) = PA(2)*pc_ke(2,15) + inv_2zeta*(2*pc_ke(2,10)) + 2*xi*(pc_ol(8,15&
&))

pc_ke(9,15) = PA(3)*pc_ke(2,15) + inv_2zeta*(pc_ke(2,6)) + 2*xi*(pc_ol(9,15))

pc_ke(10,15) = PA(3)*pc_ke(3,15) + inv_2zeta*(pc_ke(3,6)) + 2*xi*(pc_ol(10,15)&
&)

pc_ke(11,15) = PA(3)*pc_ke(8,15) + inv_2zeta*(pc_ke(8,6)) + 2*xi*(pc_ol(11,15)&
&)

pc_ke(12,15) = PA(2)*pc_ke(5,15) + inv_2zeta*(2*pc_ke(5,10)) + 2*xi*(pc_ol(12,&
&15))

pc_ke(13,15) = PA(3)*pc_ke(5,15) + inv_2zeta*(pc_ke(5,6)) + 2*xi*(pc_ol(13,15)&
&)

pc_ke(14,15) = PA(1)*pc_ke(6,15) + 2*xi*(pc_ol(14,15))

pc_ke(15,15) = PA(3)*pc_ke(6,15) + inv_2zeta*(pc_ke(6,6)) + 2*xi*(pc_ol(15,15)&
&)

pc_ke(16,15) = PA(1)*pc_ke(7,15) + 2*xi*(pc_ol(16,15))

pc_ke(17,15) = PA(2)*pc_ke(7,15) + inv_2zeta*(2*pc_ke(7,10)) + 2*xi*(pc_ol(17,&
&15))

pc_ke(18,15) = PA(1)*pc_ke(5,15) + inv_2zeta*(2*pc_ke(2,15)) + 2*xi*(pc_ol(18,&
&15) - inv_2zeta_a*2*pc_ol(2,15))

pc_ke(19,15) = PA(2)*pc_ke(6,15) + inv_2zeta*(2*pc_ke(3,15) + 2*pc_ke(6,10)) +&
& 2*xi*(pc_ol(19,15) - inv_2zeta_a*2*pc_ol(3,15))

pc_ke(20,15) = PA(3)*pc_ke(7,15) + inv_2zeta*(2*pc_ke(4,15) + pc_ke(7,6)) + 2*&
&xi*(pc_ol(20,15) - inv_2zeta_a*2*pc_ol(4,15))

pc_ke(1,16) = PB(1)*pc_ke(1,7) + 2*xi*(pc_ol(1,16))

pc_ke(2,16) = PA(1)*pc_ke(1,16) + inv_2zeta*(pc_ke(1,7)) + 2*xi*(pc_ol(2,16))

pc_ke(3,16) = PA(2)*pc_ke(1,16) + 2*xi*(pc_ol(3,16))

pc_ke(4,16) = PA(3)*pc_ke(1,16) + inv_2zeta*(2*pc_ke(1,9)) + 2*xi*(pc_ol(4,16)&
&)

pc_ke(5,16) = PA(1)*pc_ke(2,16) + inv_2zeta*(pc_ke(1,16) + pc_ke(2,7)) + 2*xi*&
&(pc_ol(5,16) - inv_2zeta_a*pc_ol(1,16))

pc_ke(6,16) = PA(2)*pc_ke(3,16) + inv_2zeta*(pc_ke(1,16)) + 2*xi*(pc_ol(6,16) &
&- inv_2zeta_a*pc_ol(1,16))

pc_ke(7,16) = PA(3)*pc_ke(4,16) + inv_2zeta*(pc_ke(1,16) + 2*pc_ke(4,9)) + 2*x&
&i*(pc_ol(7,16) - inv_2zeta_a*pc_ol(1,16))

pc_ke(8,16) = PA(2)*pc_ke(2,16) + 2*xi*(pc_ol(8,16))

pc_ke(9,16) = PA(3)*pc_ke(2,16) + inv_2zeta*(2*pc_ke(2,9)) + 2*xi*(pc_ol(9,16)&
&)

pc_ke(10,16) = PA(3)*pc_ke(3,16) + inv_2zeta*(2*pc_ke(3,9)) + 2*xi*(pc_ol(10,1&
&6))

pc_ke(11,16) = PA(3)*pc_ke(8,16) + inv_2zeta*(2*pc_ke(8,9)) + 2*xi*(pc_ol(11,1&
&6))

pc_ke(12,16) = PA(2)*pc_ke(5,16) + 2*xi*(pc_ol(12,16))

pc_ke(13,16) = PA(3)*pc_ke(5,16) + inv_2zeta*(2*pc_ke(5,9)) + 2*xi*(pc_ol(13,1&
&6))

pc_ke(14,16) = PA(1)*pc_ke(6,16) + inv_2zeta*(pc_ke(6,7)) + 2*xi*(pc_ol(14,16)&
&)

pc_ke(15,16) = PA(3)*pc_ke(6,16) + inv_2zeta*(2*pc_ke(6,9)) + 2*xi*(pc_ol(15,1&
&6))

pc_ke(16,16) = PA(1)*pc_ke(7,16) + inv_2zeta*(pc_ke(7,7)) + 2*xi*(pc_ol(16,16)&
&)

pc_ke(17,16) = PA(2)*pc_ke(7,16) + 2*xi*(pc_ol(17,16))

pc_ke(18,16) = PA(1)*pc_ke(5,16) + inv_2zeta*(2*pc_ke(2,16) + pc_ke(5,7)) + 2*&
&xi*(pc_ol(18,16) - inv_2zeta_a*2*pc_ol(2,16))

pc_ke(19,16) = PA(2)*pc_ke(6,16) + inv_2zeta*(2*pc_ke(3,16)) + 2*xi*(pc_ol(19,&
&16) - inv_2zeta_a*2*pc_ol(3,16))

pc_ke(20,16) = PA(3)*pc_ke(7,16) + inv_2zeta*(2*pc_ke(4,16) + 2*pc_ke(7,9)) + &
&2*xi*(pc_ol(20,16) - inv_2zeta_a*2*pc_ol(4,16))

pc_ke(1,17) = PB(2)*pc_ke(1,7) + 2*xi*(pc_ol(1,17))

pc_ke(2,17) = PA(1)*pc_ke(1,17) + 2*xi*(pc_ol(2,17))

pc_ke(3,17) = PA(2)*pc_ke(1,17) + inv_2zeta*(pc_ke(1,7)) + 2*xi*(pc_ol(3,17))

pc_ke(4,17) = PA(3)*pc_ke(1,17) + inv_2zeta*(2*pc_ke(1,10)) + 2*xi*(pc_ol(4,17&
&))

pc_ke(5,17) = PA(1)*pc_ke(2,17) + inv_2zeta*(pc_ke(1,17)) + 2*xi*(pc_ol(5,17) &
&- inv_2zeta_a*pc_ol(1,17))

pc_ke(6,17) = PA(2)*pc_ke(3,17) + inv_2zeta*(pc_ke(1,17) + pc_ke(3,7)) + 2*xi*&
&(pc_ol(6,17) - inv_2zeta_a*pc_ol(1,17))

pc_ke(7,17) = PA(3)*pc_ke(4,17) + inv_2zeta*(pc_ke(1,17) + 2*pc_ke(4,10)) + 2*&
&xi*(pc_ol(7,17) - inv_2zeta_a*pc_ol(1,17))

pc_ke(8,17) = PA(2)*pc_ke(2,17) + inv_2zeta*(pc_ke(2,7)) + 2*xi*(pc_ol(8,17))

pc_ke(9,17) = PA(3)*pc_ke(2,17) + inv_2zeta*(2*pc_ke(2,10)) + 2*xi*(pc_ol(9,17&
&))

pc_ke(10,17) = PA(3)*pc_ke(3,17) + inv_2zeta*(2*pc_ke(3,10)) + 2*xi*(pc_ol(10,&
&17))

pc_ke(11,17) = PA(3)*pc_ke(8,17) + inv_2zeta*(2*pc_ke(8,10)) + 2*xi*(pc_ol(11,&
&17))

pc_ke(12,17) = PA(2)*pc_ke(5,17) + inv_2zeta*(pc_ke(5,7)) + 2*xi*(pc_ol(12,17)&
&)

pc_ke(13,17) = PA(3)*pc_ke(5,17) + inv_2zeta*(2*pc_ke(5,10)) + 2*xi*(pc_ol(13,&
&17))

pc_ke(14,17) = PA(1)*pc_ke(6,17) + 2*xi*(pc_ol(14,17))

pc_ke(15,17) = PA(3)*pc_ke(6,17) + inv_2zeta*(2*pc_ke(6,10)) + 2*xi*(pc_ol(15,&
&17))

pc_ke(16,17) = PA(1)*pc_ke(7,17) + 2*xi*(pc_ol(16,17))

pc_ke(17,17) = PA(2)*pc_ke(7,17) + inv_2zeta*(pc_ke(7,7)) + 2*xi*(pc_ol(17,17)&
&)

pc_ke(18,17) = PA(1)*pc_ke(5,17) + inv_2zeta*(2*pc_ke(2,17)) + 2*xi*(pc_ol(18,&
&17) - inv_2zeta_a*2*pc_ol(2,17))

pc_ke(19,17) = PA(2)*pc_ke(6,17) + inv_2zeta*(2*pc_ke(3,17) + pc_ke(6,7)) + 2*&
&xi*(pc_ol(19,17) - inv_2zeta_a*2*pc_ol(3,17))

pc_ke(20,17) = PA(3)*pc_ke(7,17) + inv_2zeta*(2*pc_ke(4,17) + 2*pc_ke(7,10)) +&
& 2*xi*(pc_ol(20,17) - inv_2zeta_a*2*pc_ol(4,17))

pc_ke(1,18) = PB(1)*pc_ke(1,5) + inv_2zeta*(2*pc_ke(1,2)) + 2*xi*(pc_ol(1,18) &
&- inv_2zeta_b*2*pc_ol(1,2))

pc_ke(2,18) = PA(1)*pc_ke(1,18) + inv_2zeta*(3*pc_ke(1,5)) + 2*xi*(pc_ol(2,18)&
&)

pc_ke(3,18) = PA(2)*pc_ke(1,18) + 2*xi*(pc_ol(3,18))

pc_ke(4,18) = PA(3)*pc_ke(1,18) + 2*xi*(pc_ol(4,18))

pc_ke(5,18) = PA(1)*pc_ke(2,18) + inv_2zeta*(pc_ke(1,18) + 3*pc_ke(2,5)) + 2*x&
&i*(pc_ol(5,18) - inv_2zeta_a*pc_ol(1,18))

pc_ke(6,18) = PA(2)*pc_ke(3,18) + inv_2zeta*(pc_ke(1,18)) + 2*xi*(pc_ol(6,18) &
&- inv_2zeta_a*pc_ol(1,18))

pc_ke(7,18) = PA(3)*pc_ke(4,18) + inv_2zeta*(pc_ke(1,18)) + 2*xi*(pc_ol(7,18) &
&- inv_2zeta_a*pc_ol(1,18))

pc_ke(8,18) = PA(2)*pc_ke(2,18) + 2*xi*(pc_ol(8,18))

pc_ke(9,18) = PA(3)*pc_ke(2,18) + 2*xi*(pc_ol(9,18))

pc_ke(10,18) = PA(3)*pc_ke(3,18) + 2*xi*(pc_ol(10,18))

pc_ke(11,18) = PA(3)*pc_ke(8,18) + 2*xi*(pc_ol(11,18))

pc_ke(12,18) = PA(2)*pc_ke(5,18) + 2*xi*(pc_ol(12,18))

pc_ke(13,18) = PA(3)*pc_ke(5,18) + 2*xi*(pc_ol(13,18))

pc_ke(14,18) = PA(1)*pc_ke(6,18) + inv_2zeta*(3*pc_ke(6,5)) + 2*xi*(pc_ol(14,1&
&8))

pc_ke(15,18) = PA(3)*pc_ke(6,18) + 2*xi*(pc_ol(15,18))

pc_ke(16,18) = PA(1)*pc_ke(7,18) + inv_2zeta*(3*pc_ke(7,5)) + 2*xi*(pc_ol(16,1&
&8))

pc_ke(17,18) = PA(2)*pc_ke(7,18) + 2*xi*(pc_ol(17,18))

pc_ke(18,18) = PA(1)*pc_ke(5,18) + inv_2zeta*(2*pc_ke(2,18) + 3*pc_ke(5,5)) + &
&2*xi*(pc_ol(18,18) - inv_2zeta_a*2*pc_ol(2,18))

pc_ke(19,18) = PA(2)*pc_ke(6,18) + inv_2zeta*(2*pc_ke(3,18)) + 2*xi*(pc_ol(19,&
&18) - inv_2zeta_a*2*pc_ol(3,18))

pc_ke(20,18) = PA(3)*pc_ke(7,18) + inv_2zeta*(2*pc_ke(4,18)) + 2*xi*(pc_ol(20,&
&18) - inv_2zeta_a*2*pc_ol(4,18))

pc_ke(1,19) = PB(2)*pc_ke(1,6) + inv_2zeta*(2*pc_ke(1,3)) + 2*xi*(pc_ol(1,19) &
&- inv_2zeta_b*2*pc_ol(1,3))

pc_ke(2,19) = PA(1)*pc_ke(1,19) + 2*xi*(pc_ol(2,19))

pc_ke(3,19) = PA(2)*pc_ke(1,19) + inv_2zeta*(3*pc_ke(1,6)) + 2*xi*(pc_ol(3,19)&
&)

pc_ke(4,19) = PA(3)*pc_ke(1,19) + 2*xi*(pc_ol(4,19))

pc_ke(5,19) = PA(1)*pc_ke(2,19) + inv_2zeta*(pc_ke(1,19)) + 2*xi*(pc_ol(5,19) &
&- inv_2zeta_a*pc_ol(1,19))

pc_ke(6,19) = PA(2)*pc_ke(3,19) + inv_2zeta*(pc_ke(1,19) + 3*pc_ke(3,6)) + 2*x&
&i*(pc_ol(6,19) - inv_2zeta_a*pc_ol(1,19))

pc_ke(7,19) = PA(3)*pc_ke(4,19) + inv_2zeta*(pc_ke(1,19)) + 2*xi*(pc_ol(7,19) &
&- inv_2zeta_a*pc_ol(1,19))

pc_ke(8,19) = PA(2)*pc_ke(2,19) + inv_2zeta*(3*pc_ke(2,6)) + 2*xi*(pc_ol(8,19)&
&)

pc_ke(9,19) = PA(3)*pc_ke(2,19) + 2*xi*(pc_ol(9,19))

pc_ke(10,19) = PA(3)*pc_ke(3,19) + 2*xi*(pc_ol(10,19))

pc_ke(11,19) = PA(3)*pc_ke(8,19) + 2*xi*(pc_ol(11,19))

pc_ke(12,19) = PA(2)*pc_ke(5,19) + inv_2zeta*(3*pc_ke(5,6)) + 2*xi*(pc_ol(12,1&
&9))

pc_ke(13,19) = PA(3)*pc_ke(5,19) + 2*xi*(pc_ol(13,19))

pc_ke(14,19) = PA(1)*pc_ke(6,19) + 2*xi*(pc_ol(14,19))

pc_ke(15,19) = PA(3)*pc_ke(6,19) + 2*xi*(pc_ol(15,19))

pc_ke(16,19) = PA(1)*pc_ke(7,19) + 2*xi*(pc_ol(16,19))

pc_ke(17,19) = PA(2)*pc_ke(7,19) + inv_2zeta*(3*pc_ke(7,6)) + 2*xi*(pc_ol(17,1&
&9))

pc_ke(18,19) = PA(1)*pc_ke(5,19) + inv_2zeta*(2*pc_ke(2,19)) + 2*xi*(pc_ol(18,&
&19) - inv_2zeta_a*2*pc_ol(2,19))

pc_ke(19,19) = PA(2)*pc_ke(6,19) + inv_2zeta*(2*pc_ke(3,19) + 3*pc_ke(6,6)) + &
&2*xi*(pc_ol(19,19) - inv_2zeta_a*2*pc_ol(3,19))

pc_ke(20,19) = PA(3)*pc_ke(7,19) + inv_2zeta*(2*pc_ke(4,19)) + 2*xi*(pc_ol(20,&
&19) - inv_2zeta_a*2*pc_ol(4,19))

pc_ke(1,20) = PB(3)*pc_ke(1,7) + inv_2zeta*(2*pc_ke(1,4)) + 2*xi*(pc_ol(1,20) &
&- inv_2zeta_b*2*pc_ol(1,4))

pc_ke(2,20) = PA(1)*pc_ke(1,20) + 2*xi*(pc_ol(2,20))

pc_ke(3,20) = PA(2)*pc_ke(1,20) + 2*xi*(pc_ol(3,20))

pc_ke(4,20) = PA(3)*pc_ke(1,20) + inv_2zeta*(3*pc_ke(1,7)) + 2*xi*(pc_ol(4,20)&
&)

pc_ke(5,20) = PA(1)*pc_ke(2,20) + inv_2zeta*(pc_ke(1,20)) + 2*xi*(pc_ol(5,20) &
&- inv_2zeta_a*pc_ol(1,20))

pc_ke(6,20) = PA(2)*pc_ke(3,20) + inv_2zeta*(pc_ke(1,20)) + 2*xi*(pc_ol(6,20) &
&- inv_2zeta_a*pc_ol(1,20))

pc_ke(7,20) = PA(3)*pc_ke(4,20) + inv_2zeta*(pc_ke(1,20) + 3*pc_ke(4,7)) + 2*x&
&i*(pc_ol(7,20) - inv_2zeta_a*pc_ol(1,20))

pc_ke(8,20) = PA(2)*pc_ke(2,20) + 2*xi*(pc_ol(8,20))

pc_ke(9,20) = PA(3)*pc_ke(2,20) + inv_2zeta*(3*pc_ke(2,7)) + 2*xi*(pc_ol(9,20)&
&)

pc_ke(10,20) = PA(3)*pc_ke(3,20) + inv_2zeta*(3*pc_ke(3,7)) + 2*xi*(pc_ol(10,2&
&0))

pc_ke(11,20) = PA(3)*pc_ke(8,20) + inv_2zeta*(3*pc_ke(8,7)) + 2*xi*(pc_ol(11,2&
&0))

pc_ke(12,20) = PA(2)*pc_ke(5,20) + 2*xi*(pc_ol(12,20))

pc_ke(13,20) = PA(3)*pc_ke(5,20) + inv_2zeta*(3*pc_ke(5,7)) + 2*xi*(pc_ol(13,2&
&0))

pc_ke(14,20) = PA(1)*pc_ke(6,20) + 2*xi*(pc_ol(14,20))

pc_ke(15,20) = PA(3)*pc_ke(6,20) + inv_2zeta*(3*pc_ke(6,7)) + 2*xi*(pc_ol(15,2&
&0))

pc_ke(16,20) = PA(1)*pc_ke(7,20) + 2*xi*(pc_ol(16,20))

pc_ke(17,20) = PA(2)*pc_ke(7,20) + 2*xi*(pc_ol(17,20))

pc_ke(18,20) = PA(1)*pc_ke(5,20) + inv_2zeta*(2*pc_ke(2,20)) + 2*xi*(pc_ol(18,&
&20) - inv_2zeta_a*2*pc_ol(2,20))

pc_ke(19,20) = PA(2)*pc_ke(6,20) + inv_2zeta*(2*pc_ke(3,20)) + 2*xi*(pc_ol(19,&
&20) - inv_2zeta_a*2*pc_ol(3,20))

pc_ke(20,20) = PA(3)*pc_ke(7,20) + inv_2zeta*(2*pc_ke(4,20) + 3*pc_ke(7,7)) + &
&2*xi*(pc_ol(20,20) - inv_2zeta_a*2*pc_ol(4,20))

pc(1,1) = preFactorMV

pc(2,1) = PA(1)*preFactorMV + 16*xi*inv_8m3c2*(pc_ke(2,1))

pc(3,1) = PA(2)*preFactorMV + 16*xi*inv_8m3c2*(pc_ke(3,1))

pc(4,1) = PA(3)*preFactorMV + 16*xi*inv_8m3c2*(pc_ke(4,1))

pc(5,1) = PA(1)*pc(2,1) + inv_2zeta*(preFactorMV) + 16*xi*inv_8m3c2*(pc_ke(5,1&
&) - inv_2zeta_a*preFactorKE)

pc(6,1) = PA(2)*pc(3,1) + inv_2zeta*(preFactorMV) + 16*xi*inv_8m3c2*(pc_ke(6,1&
&) - inv_2zeta_a*preFactorKE)

pc(7,1) = PA(3)*pc(4,1) + inv_2zeta*(preFactorMV) + 16*xi*inv_8m3c2*(pc_ke(7,1&
&) - inv_2zeta_a*preFactorKE)

pc(8,1) = PA(2)*pc(2,1) + 16*xi*inv_8m3c2*(pc_ke(8,1))

pc(9,1) = PA(3)*pc(2,1) + 16*xi*inv_8m3c2*(pc_ke(9,1))

pc(10,1) = PA(3)*pc(3,1) + 16*xi*inv_8m3c2*(pc_ke(10,1))

pc(11,1) = PA(3)*pc(8,1) + 16*xi*inv_8m3c2*(pc_ke(11,1))

pc(12,1) = PA(2)*pc(5,1) + 16*xi*inv_8m3c2*(pc_ke(12,1))

pc(13,1) = PA(3)*pc(5,1) + 16*xi*inv_8m3c2*(pc_ke(13,1))

pc(14,1) = PA(1)*pc(6,1) + 16*xi*inv_8m3c2*(pc_ke(14,1))

pc(15,1) = PA(3)*pc(6,1) + 16*xi*inv_8m3c2*(pc_ke(15,1))

pc(16,1) = PA(1)*pc(7,1) + 16*xi*inv_8m3c2*(pc_ke(16,1))

pc(17,1) = PA(2)*pc(7,1) + 16*xi*inv_8m3c2*(pc_ke(17,1))

pc(18,1) = PA(1)*pc(5,1) + inv_2zeta*(2*pc(2,1)) + 16*xi*inv_8m3c2*(pc_ke(18,1&
&) - inv_2zeta_a*2*pc_ke(2,1))

pc(19,1) = PA(2)*pc(6,1) + inv_2zeta*(2*pc(3,1)) + 16*xi*inv_8m3c2*(pc_ke(19,1&
&) - inv_2zeta_a*2*pc_ke(3,1))

pc(20,1) = PA(3)*pc(7,1) + inv_2zeta*(2*pc(4,1)) + 16*xi*inv_8m3c2*(pc_ke(20,1&
&) - inv_2zeta_a*2*pc_ke(4,1))

pc(1,2) = PB(1)*preFactorMV + 16*xi*inv_8m3c2*(pc_ke(1,2))

pc(2,2) = PA(1)*pc(1,2) + inv_2zeta*(preFactorMV) + 16*xi*inv_8m3c2*(pc_ke(2,2&
&))

pc(3,2) = PA(2)*pc(1,2) + 16*xi*inv_8m3c2*(pc_ke(3,2))

pc(4,2) = PA(3)*pc(1,2) + 16*xi*inv_8m3c2*(pc_ke(4,2))

pc(5,2) = PA(1)*pc(2,2) + inv_2zeta*(pc(1,2) + pc(2,1)) + 16*xi*inv_8m3c2*(pc_&
&ke(5,2) - inv_2zeta_a*pc_ke(1,2))

pc(6,2) = PA(2)*pc(3,2) + inv_2zeta*(pc(1,2)) + 16*xi*inv_8m3c2*(pc_ke(6,2) - &
&inv_2zeta_a*pc_ke(1,2))

pc(7,2) = PA(3)*pc(4,2) + inv_2zeta*(pc(1,2)) + 16*xi*inv_8m3c2*(pc_ke(7,2) - &
&inv_2zeta_a*pc_ke(1,2))

pc(8,2) = PA(2)*pc(2,2) + 16*xi*inv_8m3c2*(pc_ke(8,2))

pc(9,2) = PA(3)*pc(2,2) + 16*xi*inv_8m3c2*(pc_ke(9,2))

pc(10,2) = PA(3)*pc(3,2) + 16*xi*inv_8m3c2*(pc_ke(10,2))

pc(11,2) = PA(3)*pc(8,2) + 16*xi*inv_8m3c2*(pc_ke(11,2))

pc(12,2) = PA(2)*pc(5,2) + 16*xi*inv_8m3c2*(pc_ke(12,2))

pc(13,2) = PA(3)*pc(5,2) + 16*xi*inv_8m3c2*(pc_ke(13,2))

pc(14,2) = PA(1)*pc(6,2) + inv_2zeta*(pc(6,1)) + 16*xi*inv_8m3c2*(pc_ke(14,2))

pc(15,2) = PA(3)*pc(6,2) + 16*xi*inv_8m3c2*(pc_ke(15,2))

pc(16,2) = PA(1)*pc(7,2) + inv_2zeta*(pc(7,1)) + 16*xi*inv_8m3c2*(pc_ke(16,2))

pc(17,2) = PA(2)*pc(7,2) + 16*xi*inv_8m3c2*(pc_ke(17,2))

pc(18,2) = PA(1)*pc(5,2) + inv_2zeta*(2*pc(2,2) + pc(5,1)) + 16*xi*inv_8m3c2*(&
&pc_ke(18,2) - inv_2zeta_a*2*pc_ke(2,2))

pc(19,2) = PA(2)*pc(6,2) + inv_2zeta*(2*pc(3,2)) + 16*xi*inv_8m3c2*(pc_ke(19,2&
&) - inv_2zeta_a*2*pc_ke(3,2))

pc(20,2) = PA(3)*pc(7,2) + inv_2zeta*(2*pc(4,2)) + 16*xi*inv_8m3c2*(pc_ke(20,2&
&) - inv_2zeta_a*2*pc_ke(4,2))

pc(1,3) = PB(2)*preFactorMV + 16*xi*inv_8m3c2*(pc_ke(1,3))

pc(2,3) = PA(1)*pc(1,3) + 16*xi*inv_8m3c2*(pc_ke(2,3))

pc(3,3) = PA(2)*pc(1,3) + inv_2zeta*(preFactorMV) + 16*xi*inv_8m3c2*(pc_ke(3,3&
&))

pc(4,3) = PA(3)*pc(1,3) + 16*xi*inv_8m3c2*(pc_ke(4,3))

pc(5,3) = PA(1)*pc(2,3) + inv_2zeta*(pc(1,3)) + 16*xi*inv_8m3c2*(pc_ke(5,3) - &
&inv_2zeta_a*pc_ke(1,3))

pc(6,3) = PA(2)*pc(3,3) + inv_2zeta*(pc(1,3) + pc(3,1)) + 16*xi*inv_8m3c2*(pc_&
&ke(6,3) - inv_2zeta_a*pc_ke(1,3))

pc(7,3) = PA(3)*pc(4,3) + inv_2zeta*(pc(1,3)) + 16*xi*inv_8m3c2*(pc_ke(7,3) - &
&inv_2zeta_a*pc_ke(1,3))

pc(8,3) = PA(2)*pc(2,3) + inv_2zeta*(pc(2,1)) + 16*xi*inv_8m3c2*(pc_ke(8,3))

pc(9,3) = PA(3)*pc(2,3) + 16*xi*inv_8m3c2*(pc_ke(9,3))

pc(10,3) = PA(3)*pc(3,3) + 16*xi*inv_8m3c2*(pc_ke(10,3))

pc(11,3) = PA(3)*pc(8,3) + 16*xi*inv_8m3c2*(pc_ke(11,3))

pc(12,3) = PA(2)*pc(5,3) + inv_2zeta*(pc(5,1)) + 16*xi*inv_8m3c2*(pc_ke(12,3))

pc(13,3) = PA(3)*pc(5,3) + 16*xi*inv_8m3c2*(pc_ke(13,3))

pc(14,3) = PA(1)*pc(6,3) + 16*xi*inv_8m3c2*(pc_ke(14,3))

pc(15,3) = PA(3)*pc(6,3) + 16*xi*inv_8m3c2*(pc_ke(15,3))

pc(16,3) = PA(1)*pc(7,3) + 16*xi*inv_8m3c2*(pc_ke(16,3))

pc(17,3) = PA(2)*pc(7,3) + inv_2zeta*(pc(7,1)) + 16*xi*inv_8m3c2*(pc_ke(17,3))

pc(18,3) = PA(1)*pc(5,3) + inv_2zeta*(2*pc(2,3)) + 16*xi*inv_8m3c2*(pc_ke(18,3&
&) - inv_2zeta_a*2*pc_ke(2,3))

pc(19,3) = PA(2)*pc(6,3) + inv_2zeta*(2*pc(3,3) + pc(6,1)) + 16*xi*inv_8m3c2*(&
&pc_ke(19,3) - inv_2zeta_a*2*pc_ke(3,3))

pc(20,3) = PA(3)*pc(7,3) + inv_2zeta*(2*pc(4,3)) + 16*xi*inv_8m3c2*(pc_ke(20,3&
&) - inv_2zeta_a*2*pc_ke(4,3))

pc(1,4) = PB(3)*preFactorMV + 16*xi*inv_8m3c2*(pc_ke(1,4))

pc(2,4) = PA(1)*pc(1,4) + 16*xi*inv_8m3c2*(pc_ke(2,4))

pc(3,4) = PA(2)*pc(1,4) + 16*xi*inv_8m3c2*(pc_ke(3,4))

pc(4,4) = PA(3)*pc(1,4) + inv_2zeta*(preFactorMV) + 16*xi*inv_8m3c2*(pc_ke(4,4&
&))

pc(5,4) = PA(1)*pc(2,4) + inv_2zeta*(pc(1,4)) + 16*xi*inv_8m3c2*(pc_ke(5,4) - &
&inv_2zeta_a*pc_ke(1,4))

pc(6,4) = PA(2)*pc(3,4) + inv_2zeta*(pc(1,4)) + 16*xi*inv_8m3c2*(pc_ke(6,4) - &
&inv_2zeta_a*pc_ke(1,4))

pc(7,4) = PA(3)*pc(4,4) + inv_2zeta*(pc(1,4) + pc(4,1)) + 16*xi*inv_8m3c2*(pc_&
&ke(7,4) - inv_2zeta_a*pc_ke(1,4))

pc(8,4) = PA(2)*pc(2,4) + 16*xi*inv_8m3c2*(pc_ke(8,4))

pc(9,4) = PA(3)*pc(2,4) + inv_2zeta*(pc(2,1)) + 16*xi*inv_8m3c2*(pc_ke(9,4))

pc(10,4) = PA(3)*pc(3,4) + inv_2zeta*(pc(3,1)) + 16*xi*inv_8m3c2*(pc_ke(10,4))

pc(11,4) = PA(3)*pc(8,4) + inv_2zeta*(pc(8,1)) + 16*xi*inv_8m3c2*(pc_ke(11,4))

pc(12,4) = PA(2)*pc(5,4) + 16*xi*inv_8m3c2*(pc_ke(12,4))

pc(13,4) = PA(3)*pc(5,4) + inv_2zeta*(pc(5,1)) + 16*xi*inv_8m3c2*(pc_ke(13,4))

pc(14,4) = PA(1)*pc(6,4) + 16*xi*inv_8m3c2*(pc_ke(14,4))

pc(15,4) = PA(3)*pc(6,4) + inv_2zeta*(pc(6,1)) + 16*xi*inv_8m3c2*(pc_ke(15,4))

pc(16,4) = PA(1)*pc(7,4) + 16*xi*inv_8m3c2*(pc_ke(16,4))

pc(17,4) = PA(2)*pc(7,4) + 16*xi*inv_8m3c2*(pc_ke(17,4))

pc(18,4) = PA(1)*pc(5,4) + inv_2zeta*(2*pc(2,4)) + 16*xi*inv_8m3c2*(pc_ke(18,4&
&) - inv_2zeta_a*2*pc_ke(2,4))

pc(19,4) = PA(2)*pc(6,4) + inv_2zeta*(2*pc(3,4)) + 16*xi*inv_8m3c2*(pc_ke(19,4&
&) - inv_2zeta_a*2*pc_ke(3,4))

pc(20,4) = PA(3)*pc(7,4) + inv_2zeta*(2*pc(4,4) + pc(7,1)) + 16*xi*inv_8m3c2*(&
&pc_ke(20,4) - inv_2zeta_a*2*pc_ke(4,4))

pc(1,5) = PB(1)*pc(1,2) + inv_2zeta*(preFactorMV) + 16*xi*inv_8m3c2*(pc_ke(1,5&
&) - inv_2zeta_b*preFactorKE)

pc(2,5) = PA(1)*pc(1,5) + inv_2zeta*(2*pc(1,2)) + 16*xi*inv_8m3c2*(pc_ke(2,5))

pc(3,5) = PA(2)*pc(1,5) + 16*xi*inv_8m3c2*(pc_ke(3,5))

pc(4,5) = PA(3)*pc(1,5) + 16*xi*inv_8m3c2*(pc_ke(4,5))

pc(5,5) = PA(1)*pc(2,5) + inv_2zeta*(pc(1,5) + 2*pc(2,2)) + 16*xi*inv_8m3c2*(p&
&c_ke(5,5) - inv_2zeta_a*pc_ke(1,5))

pc(6,5) = PA(2)*pc(3,5) + inv_2zeta*(pc(1,5)) + 16*xi*inv_8m3c2*(pc_ke(6,5) - &
&inv_2zeta_a*pc_ke(1,5))

pc(7,5) = PA(3)*pc(4,5) + inv_2zeta*(pc(1,5)) + 16*xi*inv_8m3c2*(pc_ke(7,5) - &
&inv_2zeta_a*pc_ke(1,5))

pc(8,5) = PA(2)*pc(2,5) + 16*xi*inv_8m3c2*(pc_ke(8,5))

pc(9,5) = PA(3)*pc(2,5) + 16*xi*inv_8m3c2*(pc_ke(9,5))

pc(10,5) = PA(3)*pc(3,5) + 16*xi*inv_8m3c2*(pc_ke(10,5))

pc(11,5) = PA(3)*pc(8,5) + 16*xi*inv_8m3c2*(pc_ke(11,5))

pc(12,5) = PA(2)*pc(5,5) + 16*xi*inv_8m3c2*(pc_ke(12,5))

pc(13,5) = PA(3)*pc(5,5) + 16*xi*inv_8m3c2*(pc_ke(13,5))

pc(14,5) = PA(1)*pc(6,5) + inv_2zeta*(2*pc(6,2)) + 16*xi*inv_8m3c2*(pc_ke(14,5&
&))

pc(15,5) = PA(3)*pc(6,5) + 16*xi*inv_8m3c2*(pc_ke(15,5))

pc(16,5) = PA(1)*pc(7,5) + inv_2zeta*(2*pc(7,2)) + 16*xi*inv_8m3c2*(pc_ke(16,5&
&))

pc(17,5) = PA(2)*pc(7,5) + 16*xi*inv_8m3c2*(pc_ke(17,5))

pc(18,5) = PA(1)*pc(5,5) + inv_2zeta*(2*pc(2,5) + 2*pc(5,2)) + 16*xi*inv_8m3c2&
&*(pc_ke(18,5) - inv_2zeta_a*2*pc_ke(2,5))

pc(19,5) = PA(2)*pc(6,5) + inv_2zeta*(2*pc(3,5)) + 16*xi*inv_8m3c2*(pc_ke(19,5&
&) - inv_2zeta_a*2*pc_ke(3,5))

pc(20,5) = PA(3)*pc(7,5) + inv_2zeta*(2*pc(4,5)) + 16*xi*inv_8m3c2*(pc_ke(20,5&
&) - inv_2zeta_a*2*pc_ke(4,5))

pc(1,6) = PB(2)*pc(1,3) + inv_2zeta*(preFactorMV) + 16*xi*inv_8m3c2*(pc_ke(1,6&
&) - inv_2zeta_b*preFactorKE)

pc(2,6) = PA(1)*pc(1,6) + 16*xi*inv_8m3c2*(pc_ke(2,6))

pc(3,6) = PA(2)*pc(1,6) + inv_2zeta*(2*pc(1,3)) + 16*xi*inv_8m3c2*(pc_ke(3,6))

pc(4,6) = PA(3)*pc(1,6) + 16*xi*inv_8m3c2*(pc_ke(4,6))

pc(5,6) = PA(1)*pc(2,6) + inv_2zeta*(pc(1,6)) + 16*xi*inv_8m3c2*(pc_ke(5,6) - &
&inv_2zeta_a*pc_ke(1,6))

pc(6,6) = PA(2)*pc(3,6) + inv_2zeta*(pc(1,6) + 2*pc(3,3)) + 16*xi*inv_8m3c2*(p&
&c_ke(6,6) - inv_2zeta_a*pc_ke(1,6))

pc(7,6) = PA(3)*pc(4,6) + inv_2zeta*(pc(1,6)) + 16*xi*inv_8m3c2*(pc_ke(7,6) - &
&inv_2zeta_a*pc_ke(1,6))

pc(8,6) = PA(2)*pc(2,6) + inv_2zeta*(2*pc(2,3)) + 16*xi*inv_8m3c2*(pc_ke(8,6))

pc(9,6) = PA(3)*pc(2,6) + 16*xi*inv_8m3c2*(pc_ke(9,6))

pc(10,6) = PA(3)*pc(3,6) + 16*xi*inv_8m3c2*(pc_ke(10,6))

pc(11,6) = PA(3)*pc(8,6) + 16*xi*inv_8m3c2*(pc_ke(11,6))

pc(12,6) = PA(2)*pc(5,6) + inv_2zeta*(2*pc(5,3)) + 16*xi*inv_8m3c2*(pc_ke(12,6&
&))

pc(13,6) = PA(3)*pc(5,6) + 16*xi*inv_8m3c2*(pc_ke(13,6))

pc(14,6) = PA(1)*pc(6,6) + 16*xi*inv_8m3c2*(pc_ke(14,6))

pc(15,6) = PA(3)*pc(6,6) + 16*xi*inv_8m3c2*(pc_ke(15,6))

pc(16,6) = PA(1)*pc(7,6) + 16*xi*inv_8m3c2*(pc_ke(16,6))

pc(17,6) = PA(2)*pc(7,6) + inv_2zeta*(2*pc(7,3)) + 16*xi*inv_8m3c2*(pc_ke(17,6&
&))

pc(18,6) = PA(1)*pc(5,6) + inv_2zeta*(2*pc(2,6)) + 16*xi*inv_8m3c2*(pc_ke(18,6&
&) - inv_2zeta_a*2*pc_ke(2,6))

pc(19,6) = PA(2)*pc(6,6) + inv_2zeta*(2*pc(3,6) + 2*pc(6,3)) + 16*xi*inv_8m3c2&
&*(pc_ke(19,6) - inv_2zeta_a*2*pc_ke(3,6))

pc(20,6) = PA(3)*pc(7,6) + inv_2zeta*(2*pc(4,6)) + 16*xi*inv_8m3c2*(pc_ke(20,6&
&) - inv_2zeta_a*2*pc_ke(4,6))

pc(1,7) = PB(3)*pc(1,4) + inv_2zeta*(preFactorMV) + 16*xi*inv_8m3c2*(pc_ke(1,7&
&) - inv_2zeta_b*preFactorKE)

pc(2,7) = PA(1)*pc(1,7) + 16*xi*inv_8m3c2*(pc_ke(2,7))

pc(3,7) = PA(2)*pc(1,7) + 16*xi*inv_8m3c2*(pc_ke(3,7))

pc(4,7) = PA(3)*pc(1,7) + inv_2zeta*(2*pc(1,4)) + 16*xi*inv_8m3c2*(pc_ke(4,7))

pc(5,7) = PA(1)*pc(2,7) + inv_2zeta*(pc(1,7)) + 16*xi*inv_8m3c2*(pc_ke(5,7) - &
&inv_2zeta_a*pc_ke(1,7))

pc(6,7) = PA(2)*pc(3,7) + inv_2zeta*(pc(1,7)) + 16*xi*inv_8m3c2*(pc_ke(6,7) - &
&inv_2zeta_a*pc_ke(1,7))

pc(7,7) = PA(3)*pc(4,7) + inv_2zeta*(pc(1,7) + 2*pc(4,4)) + 16*xi*inv_8m3c2*(p&
&c_ke(7,7) - inv_2zeta_a*pc_ke(1,7))

pc(8,7) = PA(2)*pc(2,7) + 16*xi*inv_8m3c2*(pc_ke(8,7))

pc(9,7) = PA(3)*pc(2,7) + inv_2zeta*(2*pc(2,4)) + 16*xi*inv_8m3c2*(pc_ke(9,7))

pc(10,7) = PA(3)*pc(3,7) + inv_2zeta*(2*pc(3,4)) + 16*xi*inv_8m3c2*(pc_ke(10,7&
&))

pc(11,7) = PA(3)*pc(8,7) + inv_2zeta*(2*pc(8,4)) + 16*xi*inv_8m3c2*(pc_ke(11,7&
&))

pc(12,7) = PA(2)*pc(5,7) + 16*xi*inv_8m3c2*(pc_ke(12,7))

pc(13,7) = PA(3)*pc(5,7) + inv_2zeta*(2*pc(5,4)) + 16*xi*inv_8m3c2*(pc_ke(13,7&
&))

pc(14,7) = PA(1)*pc(6,7) + 16*xi*inv_8m3c2*(pc_ke(14,7))

pc(15,7) = PA(3)*pc(6,7) + inv_2zeta*(2*pc(6,4)) + 16*xi*inv_8m3c2*(pc_ke(15,7&
&))

pc(16,7) = PA(1)*pc(7,7) + 16*xi*inv_8m3c2*(pc_ke(16,7))

pc(17,7) = PA(2)*pc(7,7) + 16*xi*inv_8m3c2*(pc_ke(17,7))

pc(18,7) = PA(1)*pc(5,7) + inv_2zeta*(2*pc(2,7)) + 16*xi*inv_8m3c2*(pc_ke(18,7&
&) - inv_2zeta_a*2*pc_ke(2,7))

pc(19,7) = PA(2)*pc(6,7) + inv_2zeta*(2*pc(3,7)) + 16*xi*inv_8m3c2*(pc_ke(19,7&
&) - inv_2zeta_a*2*pc_ke(3,7))

pc(20,7) = PA(3)*pc(7,7) + inv_2zeta*(2*pc(4,7) + 2*pc(7,4)) + 16*xi*inv_8m3c2&
&*(pc_ke(20,7) - inv_2zeta_a*2*pc_ke(4,7))

pc(1,8) = PB(2)*pc(1,2) + 16*xi*inv_8m3c2*(pc_ke(1,8))

pc(2,8) = PA(1)*pc(1,8) + inv_2zeta*(pc(1,3)) + 16*xi*inv_8m3c2*(pc_ke(2,8))

pc(3,8) = PA(2)*pc(1,8) + inv_2zeta*(pc(1,2)) + 16*xi*inv_8m3c2*(pc_ke(3,8))

pc(4,8) = PA(3)*pc(1,8) + 16*xi*inv_8m3c2*(pc_ke(4,8))

pc(5,8) = PA(1)*pc(2,8) + inv_2zeta*(pc(1,8) + pc(2,3)) + 16*xi*inv_8m3c2*(pc_&
&ke(5,8) - inv_2zeta_a*pc_ke(1,8))

pc(6,8) = PA(2)*pc(3,8) + inv_2zeta*(pc(1,8) + pc(3,2)) + 16*xi*inv_8m3c2*(pc_&
&ke(6,8) - inv_2zeta_a*pc_ke(1,8))

pc(7,8) = PA(3)*pc(4,8) + inv_2zeta*(pc(1,8)) + 16*xi*inv_8m3c2*(pc_ke(7,8) - &
&inv_2zeta_a*pc_ke(1,8))

pc(8,8) = PA(2)*pc(2,8) + inv_2zeta*(pc(2,2)) + 16*xi*inv_8m3c2*(pc_ke(8,8))

pc(9,8) = PA(3)*pc(2,8) + 16*xi*inv_8m3c2*(pc_ke(9,8))

pc(10,8) = PA(3)*pc(3,8) + 16*xi*inv_8m3c2*(pc_ke(10,8))

pc(11,8) = PA(3)*pc(8,8) + 16*xi*inv_8m3c2*(pc_ke(11,8))

pc(12,8) = PA(2)*pc(5,8) + inv_2zeta*(pc(5,2)) + 16*xi*inv_8m3c2*(pc_ke(12,8))

pc(13,8) = PA(3)*pc(5,8) + 16*xi*inv_8m3c2*(pc_ke(13,8))

pc(14,8) = PA(1)*pc(6,8) + inv_2zeta*(pc(6,3)) + 16*xi*inv_8m3c2*(pc_ke(14,8))

pc(15,8) = PA(3)*pc(6,8) + 16*xi*inv_8m3c2*(pc_ke(15,8))

pc(16,8) = PA(1)*pc(7,8) + inv_2zeta*(pc(7,3)) + 16*xi*inv_8m3c2*(pc_ke(16,8))

pc(17,8) = PA(2)*pc(7,8) + inv_2zeta*(pc(7,2)) + 16*xi*inv_8m3c2*(pc_ke(17,8))

pc(18,8) = PA(1)*pc(5,8) + inv_2zeta*(2*pc(2,8) + pc(5,3)) + 16*xi*inv_8m3c2*(&
&pc_ke(18,8) - inv_2zeta_a*2*pc_ke(2,8))

pc(19,8) = PA(2)*pc(6,8) + inv_2zeta*(2*pc(3,8) + pc(6,2)) + 16*xi*inv_8m3c2*(&
&pc_ke(19,8) - inv_2zeta_a*2*pc_ke(3,8))

pc(20,8) = PA(3)*pc(7,8) + inv_2zeta*(2*pc(4,8)) + 16*xi*inv_8m3c2*(pc_ke(20,8&
&) - inv_2zeta_a*2*pc_ke(4,8))

pc(1,9) = PB(3)*pc(1,2) + 16*xi*inv_8m3c2*(pc_ke(1,9))

pc(2,9) = PA(1)*pc(1,9) + inv_2zeta*(pc(1,4)) + 16*xi*inv_8m3c2*(pc_ke(2,9))

pc(3,9) = PA(2)*pc(1,9) + 16*xi*inv_8m3c2*(pc_ke(3,9))

pc(4,9) = PA(3)*pc(1,9) + inv_2zeta*(pc(1,2)) + 16*xi*inv_8m3c2*(pc_ke(4,9))

pc(5,9) = PA(1)*pc(2,9) + inv_2zeta*(pc(1,9) + pc(2,4)) + 16*xi*inv_8m3c2*(pc_&
&ke(5,9) - inv_2zeta_a*pc_ke(1,9))

pc(6,9) = PA(2)*pc(3,9) + inv_2zeta*(pc(1,9)) + 16*xi*inv_8m3c2*(pc_ke(6,9) - &
&inv_2zeta_a*pc_ke(1,9))

pc(7,9) = PA(3)*pc(4,9) + inv_2zeta*(pc(1,9) + pc(4,2)) + 16*xi*inv_8m3c2*(pc_&
&ke(7,9) - inv_2zeta_a*pc_ke(1,9))

pc(8,9) = PA(2)*pc(2,9) + 16*xi*inv_8m3c2*(pc_ke(8,9))

pc(9,9) = PA(3)*pc(2,9) + inv_2zeta*(pc(2,2)) + 16*xi*inv_8m3c2*(pc_ke(9,9))

pc(10,9) = PA(3)*pc(3,9) + inv_2zeta*(pc(3,2)) + 16*xi*inv_8m3c2*(pc_ke(10,9))

pc(11,9) = PA(3)*pc(8,9) + inv_2zeta*(pc(8,2)) + 16*xi*inv_8m3c2*(pc_ke(11,9))

pc(12,9) = PA(2)*pc(5,9) + 16*xi*inv_8m3c2*(pc_ke(12,9))

pc(13,9) = PA(3)*pc(5,9) + inv_2zeta*(pc(5,2)) + 16*xi*inv_8m3c2*(pc_ke(13,9))

pc(14,9) = PA(1)*pc(6,9) + inv_2zeta*(pc(6,4)) + 16*xi*inv_8m3c2*(pc_ke(14,9))

pc(15,9) = PA(3)*pc(6,9) + inv_2zeta*(pc(6,2)) + 16*xi*inv_8m3c2*(pc_ke(15,9))

pc(16,9) = PA(1)*pc(7,9) + inv_2zeta*(pc(7,4)) + 16*xi*inv_8m3c2*(pc_ke(16,9))

pc(17,9) = PA(2)*pc(7,9) + 16*xi*inv_8m3c2*(pc_ke(17,9))

pc(18,9) = PA(1)*pc(5,9) + inv_2zeta*(2*pc(2,9) + pc(5,4)) + 16*xi*inv_8m3c2*(&
&pc_ke(18,9) - inv_2zeta_a*2*pc_ke(2,9))

pc(19,9) = PA(2)*pc(6,9) + inv_2zeta*(2*pc(3,9)) + 16*xi*inv_8m3c2*(pc_ke(19,9&
&) - inv_2zeta_a*2*pc_ke(3,9))

pc(20,9) = PA(3)*pc(7,9) + inv_2zeta*(2*pc(4,9) + pc(7,2)) + 16*xi*inv_8m3c2*(&
&pc_ke(20,9) - inv_2zeta_a*2*pc_ke(4,9))

pc(1,10) = PB(3)*pc(1,3) + 16*xi*inv_8m3c2*(pc_ke(1,10))

pc(2,10) = PA(1)*pc(1,10) + 16*xi*inv_8m3c2*(pc_ke(2,10))

pc(3,10) = PA(2)*pc(1,10) + inv_2zeta*(pc(1,4)) + 16*xi*inv_8m3c2*(pc_ke(3,10)&
&)

pc(4,10) = PA(3)*pc(1,10) + inv_2zeta*(pc(1,3)) + 16*xi*inv_8m3c2*(pc_ke(4,10)&
&)

pc(5,10) = PA(1)*pc(2,10) + inv_2zeta*(pc(1,10)) + 16*xi*inv_8m3c2*(pc_ke(5,10&
&) - inv_2zeta_a*pc_ke(1,10))

pc(6,10) = PA(2)*pc(3,10) + inv_2zeta*(pc(1,10) + pc(3,4)) + 16*xi*inv_8m3c2*(&
&pc_ke(6,10) - inv_2zeta_a*pc_ke(1,10))

pc(7,10) = PA(3)*pc(4,10) + inv_2zeta*(pc(1,10) + pc(4,3)) + 16*xi*inv_8m3c2*(&
&pc_ke(7,10) - inv_2zeta_a*pc_ke(1,10))

pc(8,10) = PA(2)*pc(2,10) + inv_2zeta*(pc(2,4)) + 16*xi*inv_8m3c2*(pc_ke(8,10)&
&)

pc(9,10) = PA(3)*pc(2,10) + inv_2zeta*(pc(2,3)) + 16*xi*inv_8m3c2*(pc_ke(9,10)&
&)

pc(10,10) = PA(3)*pc(3,10) + inv_2zeta*(pc(3,3)) + 16*xi*inv_8m3c2*(pc_ke(10,1&
&0))

pc(11,10) = PA(3)*pc(8,10) + inv_2zeta*(pc(8,3)) + 16*xi*inv_8m3c2*(pc_ke(11,1&
&0))

pc(12,10) = PA(2)*pc(5,10) + inv_2zeta*(pc(5,4)) + 16*xi*inv_8m3c2*(pc_ke(12,1&
&0))

pc(13,10) = PA(3)*pc(5,10) + inv_2zeta*(pc(5,3)) + 16*xi*inv_8m3c2*(pc_ke(13,1&
&0))

pc(14,10) = PA(1)*pc(6,10) + 16*xi*inv_8m3c2*(pc_ke(14,10))

pc(15,10) = PA(3)*pc(6,10) + inv_2zeta*(pc(6,3)) + 16*xi*inv_8m3c2*(pc_ke(15,1&
&0))

pc(16,10) = PA(1)*pc(7,10) + 16*xi*inv_8m3c2*(pc_ke(16,10))

pc(17,10) = PA(2)*pc(7,10) + inv_2zeta*(pc(7,4)) + 16*xi*inv_8m3c2*(pc_ke(17,1&
&0))

pc(18,10) = PA(1)*pc(5,10) + inv_2zeta*(2*pc(2,10)) + 16*xi*inv_8m3c2*(pc_ke(1&
&8,10) - inv_2zeta_a*2*pc_ke(2,10))

pc(19,10) = PA(2)*pc(6,10) + inv_2zeta*(2*pc(3,10) + pc(6,4)) + 16*xi*inv_8m3c&
&2*(pc_ke(19,10) - inv_2zeta_a*2*pc_ke(3,10))

pc(20,10) = PA(3)*pc(7,10) + inv_2zeta*(2*pc(4,10) + pc(7,3)) + 16*xi*inv_8m3c&
&2*(pc_ke(20,10) - inv_2zeta_a*2*pc_ke(4,10))

pc(1,11) = PB(3)*pc(1,8) + 16*xi*inv_8m3c2*(pc_ke(1,11))

pc(2,11) = PA(1)*pc(1,11) + inv_2zeta*(pc(1,10)) + 16*xi*inv_8m3c2*(pc_ke(2,11&
&))

pc(3,11) = PA(2)*pc(1,11) + inv_2zeta*(pc(1,9)) + 16*xi*inv_8m3c2*(pc_ke(3,11)&
&)

pc(4,11) = PA(3)*pc(1,11) + inv_2zeta*(pc(1,8)) + 16*xi*inv_8m3c2*(pc_ke(4,11)&
&)

pc(5,11) = PA(1)*pc(2,11) + inv_2zeta*(pc(1,11) + pc(2,10)) + 16*xi*inv_8m3c2*&
&(pc_ke(5,11) - inv_2zeta_a*pc_ke(1,11))

pc(6,11) = PA(2)*pc(3,11) + inv_2zeta*(pc(1,11) + pc(3,9)) + 16*xi*inv_8m3c2*(&
&pc_ke(6,11) - inv_2zeta_a*pc_ke(1,11))

pc(7,11) = PA(3)*pc(4,11) + inv_2zeta*(pc(1,11) + pc(4,8)) + 16*xi*inv_8m3c2*(&
&pc_ke(7,11) - inv_2zeta_a*pc_ke(1,11))

pc(8,11) = PA(2)*pc(2,11) + inv_2zeta*(pc(2,9)) + 16*xi*inv_8m3c2*(pc_ke(8,11)&
&)

pc(9,11) = PA(3)*pc(2,11) + inv_2zeta*(pc(2,8)) + 16*xi*inv_8m3c2*(pc_ke(9,11)&
&)

pc(10,11) = PA(3)*pc(3,11) + inv_2zeta*(pc(3,8)) + 16*xi*inv_8m3c2*(pc_ke(10,1&
&1))

pc(11,11) = PA(3)*pc(8,11) + inv_2zeta*(pc(8,8)) + 16*xi*inv_8m3c2*(pc_ke(11,1&
&1))

pc(12,11) = PA(2)*pc(5,11) + inv_2zeta*(pc(5,9)) + 16*xi*inv_8m3c2*(pc_ke(12,1&
&1))

pc(13,11) = PA(3)*pc(5,11) + inv_2zeta*(pc(5,8)) + 16*xi*inv_8m3c2*(pc_ke(13,1&
&1))

pc(14,11) = PA(1)*pc(6,11) + inv_2zeta*(pc(6,10)) + 16*xi*inv_8m3c2*(pc_ke(14,&
&11))

pc(15,11) = PA(3)*pc(6,11) + inv_2zeta*(pc(6,8)) + 16*xi*inv_8m3c2*(pc_ke(15,1&
&1))

pc(16,11) = PA(1)*pc(7,11) + inv_2zeta*(pc(7,10)) + 16*xi*inv_8m3c2*(pc_ke(16,&
&11))

pc(17,11) = PA(2)*pc(7,11) + inv_2zeta*(pc(7,9)) + 16*xi*inv_8m3c2*(pc_ke(17,1&
&1))

pc(18,11) = PA(1)*pc(5,11) + inv_2zeta*(2*pc(2,11) + pc(5,10)) + 16*xi*inv_8m3&
&c2*(pc_ke(18,11) - inv_2zeta_a*2*pc_ke(2,11))

pc(19,11) = PA(2)*pc(6,11) + inv_2zeta*(2*pc(3,11) + pc(6,9)) + 16*xi*inv_8m3c&
&2*(pc_ke(19,11) - inv_2zeta_a*2*pc_ke(3,11))

pc(20,11) = PA(3)*pc(7,11) + inv_2zeta*(2*pc(4,11) + pc(7,8)) + 16*xi*inv_8m3c&
&2*(pc_ke(20,11) - inv_2zeta_a*2*pc_ke(4,11))

pc(1,12) = PB(2)*pc(1,5) + 16*xi*inv_8m3c2*(pc_ke(1,12))

pc(2,12) = PA(1)*pc(1,12) + inv_2zeta*(2*pc(1,8)) + 16*xi*inv_8m3c2*(pc_ke(2,1&
&2))

pc(3,12) = PA(2)*pc(1,12) + inv_2zeta*(pc(1,5)) + 16*xi*inv_8m3c2*(pc_ke(3,12)&
&)

pc(4,12) = PA(3)*pc(1,12) + 16*xi*inv_8m3c2*(pc_ke(4,12))

pc(5,12) = PA(1)*pc(2,12) + inv_2zeta*(pc(1,12) + 2*pc(2,8)) + 16*xi*inv_8m3c2&
&*(pc_ke(5,12) - inv_2zeta_a*pc_ke(1,12))

pc(6,12) = PA(2)*pc(3,12) + inv_2zeta*(pc(1,12) + pc(3,5)) + 16*xi*inv_8m3c2*(&
&pc_ke(6,12) - inv_2zeta_a*pc_ke(1,12))

pc(7,12) = PA(3)*pc(4,12) + inv_2zeta*(pc(1,12)) + 16*xi*inv_8m3c2*(pc_ke(7,12&
&) - inv_2zeta_a*pc_ke(1,12))

pc(8,12) = PA(2)*pc(2,12) + inv_2zeta*(pc(2,5)) + 16*xi*inv_8m3c2*(pc_ke(8,12)&
&)

pc(9,12) = PA(3)*pc(2,12) + 16*xi*inv_8m3c2*(pc_ke(9,12))

pc(10,12) = PA(3)*pc(3,12) + 16*xi*inv_8m3c2*(pc_ke(10,12))

pc(11,12) = PA(3)*pc(8,12) + 16*xi*inv_8m3c2*(pc_ke(11,12))

pc(12,12) = PA(2)*pc(5,12) + inv_2zeta*(pc(5,5)) + 16*xi*inv_8m3c2*(pc_ke(12,1&
&2))

pc(13,12) = PA(3)*pc(5,12) + 16*xi*inv_8m3c2*(pc_ke(13,12))

pc(14,12) = PA(1)*pc(6,12) + inv_2zeta*(2*pc(6,8)) + 16*xi*inv_8m3c2*(pc_ke(14&
&,12))

pc(15,12) = PA(3)*pc(6,12) + 16*xi*inv_8m3c2*(pc_ke(15,12))

pc(16,12) = PA(1)*pc(7,12) + inv_2zeta*(2*pc(7,8)) + 16*xi*inv_8m3c2*(pc_ke(16&
&,12))

pc(17,12) = PA(2)*pc(7,12) + inv_2zeta*(pc(7,5)) + 16*xi*inv_8m3c2*(pc_ke(17,1&
&2))

pc(18,12) = PA(1)*pc(5,12) + inv_2zeta*(2*pc(2,12) + 2*pc(5,8)) + 16*xi*inv_8m&
&3c2*(pc_ke(18,12) - inv_2zeta_a*2*pc_ke(2,12))

pc(19,12) = PA(2)*pc(6,12) + inv_2zeta*(2*pc(3,12) + pc(6,5)) + 16*xi*inv_8m3c&
&2*(pc_ke(19,12) - inv_2zeta_a*2*pc_ke(3,12))

pc(20,12) = PA(3)*pc(7,12) + inv_2zeta*(2*pc(4,12)) + 16*xi*inv_8m3c2*(pc_ke(2&
&0,12) - inv_2zeta_a*2*pc_ke(4,12))

pc(1,13) = PB(3)*pc(1,5) + 16*xi*inv_8m3c2*(pc_ke(1,13))

pc(2,13) = PA(1)*pc(1,13) + inv_2zeta*(2*pc(1,9)) + 16*xi*inv_8m3c2*(pc_ke(2,1&
&3))

pc(3,13) = PA(2)*pc(1,13) + 16*xi*inv_8m3c2*(pc_ke(3,13))

pc(4,13) = PA(3)*pc(1,13) + inv_2zeta*(pc(1,5)) + 16*xi*inv_8m3c2*(pc_ke(4,13)&
&)

pc(5,13) = PA(1)*pc(2,13) + inv_2zeta*(pc(1,13) + 2*pc(2,9)) + 16*xi*inv_8m3c2&
&*(pc_ke(5,13) - inv_2zeta_a*pc_ke(1,13))

pc(6,13) = PA(2)*pc(3,13) + inv_2zeta*(pc(1,13)) + 16*xi*inv_8m3c2*(pc_ke(6,13&
&) - inv_2zeta_a*pc_ke(1,13))

pc(7,13) = PA(3)*pc(4,13) + inv_2zeta*(pc(1,13) + pc(4,5)) + 16*xi*inv_8m3c2*(&
&pc_ke(7,13) - inv_2zeta_a*pc_ke(1,13))

pc(8,13) = PA(2)*pc(2,13) + 16*xi*inv_8m3c2*(pc_ke(8,13))

pc(9,13) = PA(3)*pc(2,13) + inv_2zeta*(pc(2,5)) + 16*xi*inv_8m3c2*(pc_ke(9,13)&
&)

pc(10,13) = PA(3)*pc(3,13) + inv_2zeta*(pc(3,5)) + 16*xi*inv_8m3c2*(pc_ke(10,1&
&3))

pc(11,13) = PA(3)*pc(8,13) + inv_2zeta*(pc(8,5)) + 16*xi*inv_8m3c2*(pc_ke(11,1&
&3))

pc(12,13) = PA(2)*pc(5,13) + 16*xi*inv_8m3c2*(pc_ke(12,13))

pc(13,13) = PA(3)*pc(5,13) + inv_2zeta*(pc(5,5)) + 16*xi*inv_8m3c2*(pc_ke(13,1&
&3))

pc(14,13) = PA(1)*pc(6,13) + inv_2zeta*(2*pc(6,9)) + 16*xi*inv_8m3c2*(pc_ke(14&
&,13))

pc(15,13) = PA(3)*pc(6,13) + inv_2zeta*(pc(6,5)) + 16*xi*inv_8m3c2*(pc_ke(15,1&
&3))

pc(16,13) = PA(1)*pc(7,13) + inv_2zeta*(2*pc(7,9)) + 16*xi*inv_8m3c2*(pc_ke(16&
&,13))

pc(17,13) = PA(2)*pc(7,13) + 16*xi*inv_8m3c2*(pc_ke(17,13))

pc(18,13) = PA(1)*pc(5,13) + inv_2zeta*(2*pc(2,13) + 2*pc(5,9)) + 16*xi*inv_8m&
&3c2*(pc_ke(18,13) - inv_2zeta_a*2*pc_ke(2,13))

pc(19,13) = PA(2)*pc(6,13) + inv_2zeta*(2*pc(3,13)) + 16*xi*inv_8m3c2*(pc_ke(1&
&9,13) - inv_2zeta_a*2*pc_ke(3,13))

pc(20,13) = PA(3)*pc(7,13) + inv_2zeta*(2*pc(4,13) + pc(7,5)) + 16*xi*inv_8m3c&
&2*(pc_ke(20,13) - inv_2zeta_a*2*pc_ke(4,13))

pc(1,14) = PB(1)*pc(1,6) + 16*xi*inv_8m3c2*(pc_ke(1,14))

pc(2,14) = PA(1)*pc(1,14) + inv_2zeta*(pc(1,6)) + 16*xi*inv_8m3c2*(pc_ke(2,14)&
&)

pc(3,14) = PA(2)*pc(1,14) + inv_2zeta*(2*pc(1,8)) + 16*xi*inv_8m3c2*(pc_ke(3,1&
&4))

pc(4,14) = PA(3)*pc(1,14) + 16*xi*inv_8m3c2*(pc_ke(4,14))

pc(5,14) = PA(1)*pc(2,14) + inv_2zeta*(pc(1,14) + pc(2,6)) + 16*xi*inv_8m3c2*(&
&pc_ke(5,14) - inv_2zeta_a*pc_ke(1,14))

pc(6,14) = PA(2)*pc(3,14) + inv_2zeta*(pc(1,14) + 2*pc(3,8)) + 16*xi*inv_8m3c2&
&*(pc_ke(6,14) - inv_2zeta_a*pc_ke(1,14))

pc(7,14) = PA(3)*pc(4,14) + inv_2zeta*(pc(1,14)) + 16*xi*inv_8m3c2*(pc_ke(7,14&
&) - inv_2zeta_a*pc_ke(1,14))

pc(8,14) = PA(2)*pc(2,14) + inv_2zeta*(2*pc(2,8)) + 16*xi*inv_8m3c2*(pc_ke(8,1&
&4))

pc(9,14) = PA(3)*pc(2,14) + 16*xi*inv_8m3c2*(pc_ke(9,14))

pc(10,14) = PA(3)*pc(3,14) + 16*xi*inv_8m3c2*(pc_ke(10,14))

pc(11,14) = PA(3)*pc(8,14) + 16*xi*inv_8m3c2*(pc_ke(11,14))

pc(12,14) = PA(2)*pc(5,14) + inv_2zeta*(2*pc(5,8)) + 16*xi*inv_8m3c2*(pc_ke(12&
&,14))

pc(13,14) = PA(3)*pc(5,14) + 16*xi*inv_8m3c2*(pc_ke(13,14))

pc(14,14) = PA(1)*pc(6,14) + inv_2zeta*(pc(6,6)) + 16*xi*inv_8m3c2*(pc_ke(14,1&
&4))

pc(15,14) = PA(3)*pc(6,14) + 16*xi*inv_8m3c2*(pc_ke(15,14))

pc(16,14) = PA(1)*pc(7,14) + inv_2zeta*(pc(7,6)) + 16*xi*inv_8m3c2*(pc_ke(16,1&
&4))

pc(17,14) = PA(2)*pc(7,14) + inv_2zeta*(2*pc(7,8)) + 16*xi*inv_8m3c2*(pc_ke(17&
&,14))

pc(18,14) = PA(1)*pc(5,14) + inv_2zeta*(2*pc(2,14) + pc(5,6)) + 16*xi*inv_8m3c&
&2*(pc_ke(18,14) - inv_2zeta_a*2*pc_ke(2,14))

pc(19,14) = PA(2)*pc(6,14) + inv_2zeta*(2*pc(3,14) + 2*pc(6,8)) + 16*xi*inv_8m&
&3c2*(pc_ke(19,14) - inv_2zeta_a*2*pc_ke(3,14))

pc(20,14) = PA(3)*pc(7,14) + inv_2zeta*(2*pc(4,14)) + 16*xi*inv_8m3c2*(pc_ke(2&
&0,14) - inv_2zeta_a*2*pc_ke(4,14))

pc(1,15) = PB(3)*pc(1,6) + 16*xi*inv_8m3c2*(pc_ke(1,15))

pc(2,15) = PA(1)*pc(1,15) + 16*xi*inv_8m3c2*(pc_ke(2,15))

pc(3,15) = PA(2)*pc(1,15) + inv_2zeta*(2*pc(1,10)) + 16*xi*inv_8m3c2*(pc_ke(3,&
&15))

pc(4,15) = PA(3)*pc(1,15) + inv_2zeta*(pc(1,6)) + 16*xi*inv_8m3c2*(pc_ke(4,15)&
&)

pc(5,15) = PA(1)*pc(2,15) + inv_2zeta*(pc(1,15)) + 16*xi*inv_8m3c2*(pc_ke(5,15&
&) - inv_2zeta_a*pc_ke(1,15))

pc(6,15) = PA(2)*pc(3,15) + inv_2zeta*(pc(1,15) + 2*pc(3,10)) + 16*xi*inv_8m3c&
&2*(pc_ke(6,15) - inv_2zeta_a*pc_ke(1,15))

pc(7,15) = PA(3)*pc(4,15) + inv_2zeta*(pc(1,15) + pc(4,6)) + 16*xi*inv_8m3c2*(&
&pc_ke(7,15) - inv_2zeta_a*pc_ke(1,15))

pc(8,15) = PA(2)*pc(2,15) + inv_2zeta*(2*pc(2,10)) + 16*xi*inv_8m3c2*(pc_ke(8,&
&15))

pc(9,15) = PA(3)*pc(2,15) + inv_2zeta*(pc(2,6)) + 16*xi*inv_8m3c2*(pc_ke(9,15)&
&)

pc(10,15) = PA(3)*pc(3,15) + inv_2zeta*(pc(3,6)) + 16*xi*inv_8m3c2*(pc_ke(10,1&
&5))

pc(11,15) = PA(3)*pc(8,15) + inv_2zeta*(pc(8,6)) + 16*xi*inv_8m3c2*(pc_ke(11,1&
&5))

pc(12,15) = PA(2)*pc(5,15) + inv_2zeta*(2*pc(5,10)) + 16*xi*inv_8m3c2*(pc_ke(1&
&2,15))

pc(13,15) = PA(3)*pc(5,15) + inv_2zeta*(pc(5,6)) + 16*xi*inv_8m3c2*(pc_ke(13,1&
&5))

pc(14,15) = PA(1)*pc(6,15) + 16*xi*inv_8m3c2*(pc_ke(14,15))

pc(15,15) = PA(3)*pc(6,15) + inv_2zeta*(pc(6,6)) + 16*xi*inv_8m3c2*(pc_ke(15,1&
&5))

pc(16,15) = PA(1)*pc(7,15) + 16*xi*inv_8m3c2*(pc_ke(16,15))

pc(17,15) = PA(2)*pc(7,15) + inv_2zeta*(2*pc(7,10)) + 16*xi*inv_8m3c2*(pc_ke(1&
&7,15))

pc(18,15) = PA(1)*pc(5,15) + inv_2zeta*(2*pc(2,15)) + 16*xi*inv_8m3c2*(pc_ke(1&
&8,15) - inv_2zeta_a*2*pc_ke(2,15))

pc(19,15) = PA(2)*pc(6,15) + inv_2zeta*(2*pc(3,15) + 2*pc(6,10)) + 16*xi*inv_8&
&m3c2*(pc_ke(19,15) - inv_2zeta_a*2*pc_ke(3,15))

pc(20,15) = PA(3)*pc(7,15) + inv_2zeta*(2*pc(4,15) + pc(7,6)) + 16*xi*inv_8m3c&
&2*(pc_ke(20,15) - inv_2zeta_a*2*pc_ke(4,15))

pc(1,16) = PB(1)*pc(1,7) + 16*xi*inv_8m3c2*(pc_ke(1,16))

pc(2,16) = PA(1)*pc(1,16) + inv_2zeta*(pc(1,7)) + 16*xi*inv_8m3c2*(pc_ke(2,16)&
&)

pc(3,16) = PA(2)*pc(1,16) + 16*xi*inv_8m3c2*(pc_ke(3,16))

pc(4,16) = PA(3)*pc(1,16) + inv_2zeta*(2*pc(1,9)) + 16*xi*inv_8m3c2*(pc_ke(4,1&
&6))

pc(5,16) = PA(1)*pc(2,16) + inv_2zeta*(pc(1,16) + pc(2,7)) + 16*xi*inv_8m3c2*(&
&pc_ke(5,16) - inv_2zeta_a*pc_ke(1,16))

pc(6,16) = PA(2)*pc(3,16) + inv_2zeta*(pc(1,16)) + 16*xi*inv_8m3c2*(pc_ke(6,16&
&) - inv_2zeta_a*pc_ke(1,16))

pc(7,16) = PA(3)*pc(4,16) + inv_2zeta*(pc(1,16) + 2*pc(4,9)) + 16*xi*inv_8m3c2&
&*(pc_ke(7,16) - inv_2zeta_a*pc_ke(1,16))

pc(8,16) = PA(2)*pc(2,16) + 16*xi*inv_8m3c2*(pc_ke(8,16))

pc(9,16) = PA(3)*pc(2,16) + inv_2zeta*(2*pc(2,9)) + 16*xi*inv_8m3c2*(pc_ke(9,1&
&6))

pc(10,16) = PA(3)*pc(3,16) + inv_2zeta*(2*pc(3,9)) + 16*xi*inv_8m3c2*(pc_ke(10&
&,16))

pc(11,16) = PA(3)*pc(8,16) + inv_2zeta*(2*pc(8,9)) + 16*xi*inv_8m3c2*(pc_ke(11&
&,16))

pc(12,16) = PA(2)*pc(5,16) + 16*xi*inv_8m3c2*(pc_ke(12,16))

pc(13,16) = PA(3)*pc(5,16) + inv_2zeta*(2*pc(5,9)) + 16*xi*inv_8m3c2*(pc_ke(13&
&,16))

pc(14,16) = PA(1)*pc(6,16) + inv_2zeta*(pc(6,7)) + 16*xi*inv_8m3c2*(pc_ke(14,1&
&6))

pc(15,16) = PA(3)*pc(6,16) + inv_2zeta*(2*pc(6,9)) + 16*xi*inv_8m3c2*(pc_ke(15&
&,16))

pc(16,16) = PA(1)*pc(7,16) + inv_2zeta*(pc(7,7)) + 16*xi*inv_8m3c2*(pc_ke(16,1&
&6))

pc(17,16) = PA(2)*pc(7,16) + 16*xi*inv_8m3c2*(pc_ke(17,16))

pc(18,16) = PA(1)*pc(5,16) + inv_2zeta*(2*pc(2,16) + pc(5,7)) + 16*xi*inv_8m3c&
&2*(pc_ke(18,16) - inv_2zeta_a*2*pc_ke(2,16))

pc(19,16) = PA(2)*pc(6,16) + inv_2zeta*(2*pc(3,16)) + 16*xi*inv_8m3c2*(pc_ke(1&
&9,16) - inv_2zeta_a*2*pc_ke(3,16))

pc(20,16) = PA(3)*pc(7,16) + inv_2zeta*(2*pc(4,16) + 2*pc(7,9)) + 16*xi*inv_8m&
&3c2*(pc_ke(20,16) - inv_2zeta_a*2*pc_ke(4,16))

pc(1,17) = PB(2)*pc(1,7) + 16*xi*inv_8m3c2*(pc_ke(1,17))

pc(2,17) = PA(1)*pc(1,17) + 16*xi*inv_8m3c2*(pc_ke(2,17))

pc(3,17) = PA(2)*pc(1,17) + inv_2zeta*(pc(1,7)) + 16*xi*inv_8m3c2*(pc_ke(3,17)&
&)

pc(4,17) = PA(3)*pc(1,17) + inv_2zeta*(2*pc(1,10)) + 16*xi*inv_8m3c2*(pc_ke(4,&
&17))

pc(5,17) = PA(1)*pc(2,17) + inv_2zeta*(pc(1,17)) + 16*xi*inv_8m3c2*(pc_ke(5,17&
&) - inv_2zeta_a*pc_ke(1,17))

pc(6,17) = PA(2)*pc(3,17) + inv_2zeta*(pc(1,17) + pc(3,7)) + 16*xi*inv_8m3c2*(&
&pc_ke(6,17) - inv_2zeta_a*pc_ke(1,17))

pc(7,17) = PA(3)*pc(4,17) + inv_2zeta*(pc(1,17) + 2*pc(4,10)) + 16*xi*inv_8m3c&
&2*(pc_ke(7,17) - inv_2zeta_a*pc_ke(1,17))

pc(8,17) = PA(2)*pc(2,17) + inv_2zeta*(pc(2,7)) + 16*xi*inv_8m3c2*(pc_ke(8,17)&
&)

pc(9,17) = PA(3)*pc(2,17) + inv_2zeta*(2*pc(2,10)) + 16*xi*inv_8m3c2*(pc_ke(9,&
&17))

pc(10,17) = PA(3)*pc(3,17) + inv_2zeta*(2*pc(3,10)) + 16*xi*inv_8m3c2*(pc_ke(1&
&0,17))

pc(11,17) = PA(3)*pc(8,17) + inv_2zeta*(2*pc(8,10)) + 16*xi*inv_8m3c2*(pc_ke(1&
&1,17))

pc(12,17) = PA(2)*pc(5,17) + inv_2zeta*(pc(5,7)) + 16*xi*inv_8m3c2*(pc_ke(12,1&
&7))

pc(13,17) = PA(3)*pc(5,17) + inv_2zeta*(2*pc(5,10)) + 16*xi*inv_8m3c2*(pc_ke(1&
&3,17))

pc(14,17) = PA(1)*pc(6,17) + 16*xi*inv_8m3c2*(pc_ke(14,17))

pc(15,17) = PA(3)*pc(6,17) + inv_2zeta*(2*pc(6,10)) + 16*xi*inv_8m3c2*(pc_ke(1&
&5,17))

pc(16,17) = PA(1)*pc(7,17) + 16*xi*inv_8m3c2*(pc_ke(16,17))

pc(17,17) = PA(2)*pc(7,17) + inv_2zeta*(pc(7,7)) + 16*xi*inv_8m3c2*(pc_ke(17,1&
&7))

pc(18,17) = PA(1)*pc(5,17) + inv_2zeta*(2*pc(2,17)) + 16*xi*inv_8m3c2*(pc_ke(1&
&8,17) - inv_2zeta_a*2*pc_ke(2,17))

pc(19,17) = PA(2)*pc(6,17) + inv_2zeta*(2*pc(3,17) + pc(6,7)) + 16*xi*inv_8m3c&
&2*(pc_ke(19,17) - inv_2zeta_a*2*pc_ke(3,17))

pc(20,17) = PA(3)*pc(7,17) + inv_2zeta*(2*pc(4,17) + 2*pc(7,10)) + 16*xi*inv_8&
&m3c2*(pc_ke(20,17) - inv_2zeta_a*2*pc_ke(4,17))

pc(1,18) = PB(1)*pc(1,5) + inv_2zeta*(2*pc(1,2)) + 16*xi*inv_8m3c2*(pc_ke(1,18&
&) - inv_2zeta_b*2*pc_ke(1,2))

pc(2,18) = PA(1)*pc(1,18) + inv_2zeta*(3*pc(1,5)) + 16*xi*inv_8m3c2*(pc_ke(2,1&
&8))

pc(3,18) = PA(2)*pc(1,18) + 16*xi*inv_8m3c2*(pc_ke(3,18))

pc(4,18) = PA(3)*pc(1,18) + 16*xi*inv_8m3c2*(pc_ke(4,18))

pc(5,18) = PA(1)*pc(2,18) + inv_2zeta*(pc(1,18) + 3*pc(2,5)) + 16*xi*inv_8m3c2&
&*(pc_ke(5,18) - inv_2zeta_a*pc_ke(1,18))

pc(6,18) = PA(2)*pc(3,18) + inv_2zeta*(pc(1,18)) + 16*xi*inv_8m3c2*(pc_ke(6,18&
&) - inv_2zeta_a*pc_ke(1,18))

pc(7,18) = PA(3)*pc(4,18) + inv_2zeta*(pc(1,18)) + 16*xi*inv_8m3c2*(pc_ke(7,18&
&) - inv_2zeta_a*pc_ke(1,18))

pc(8,18) = PA(2)*pc(2,18) + 16*xi*inv_8m3c2*(pc_ke(8,18))

pc(9,18) = PA(3)*pc(2,18) + 16*xi*inv_8m3c2*(pc_ke(9,18))

pc(10,18) = PA(3)*pc(3,18) + 16*xi*inv_8m3c2*(pc_ke(10,18))

pc(11,18) = PA(3)*pc(8,18) + 16*xi*inv_8m3c2*(pc_ke(11,18))

pc(12,18) = PA(2)*pc(5,18) + 16*xi*inv_8m3c2*(pc_ke(12,18))

pc(13,18) = PA(3)*pc(5,18) + 16*xi*inv_8m3c2*(pc_ke(13,18))

pc(14,18) = PA(1)*pc(6,18) + inv_2zeta*(3*pc(6,5)) + 16*xi*inv_8m3c2*(pc_ke(14&
&,18))

pc(15,18) = PA(3)*pc(6,18) + 16*xi*inv_8m3c2*(pc_ke(15,18))

pc(16,18) = PA(1)*pc(7,18) + inv_2zeta*(3*pc(7,5)) + 16*xi*inv_8m3c2*(pc_ke(16&
&,18))

pc(17,18) = PA(2)*pc(7,18) + 16*xi*inv_8m3c2*(pc_ke(17,18))

pc(18,18) = PA(1)*pc(5,18) + inv_2zeta*(2*pc(2,18) + 3*pc(5,5)) + 16*xi*inv_8m&
&3c2*(pc_ke(18,18) - inv_2zeta_a*2*pc_ke(2,18))

pc(19,18) = PA(2)*pc(6,18) + inv_2zeta*(2*pc(3,18)) + 16*xi*inv_8m3c2*(pc_ke(1&
&9,18) - inv_2zeta_a*2*pc_ke(3,18))

pc(20,18) = PA(3)*pc(7,18) + inv_2zeta*(2*pc(4,18)) + 16*xi*inv_8m3c2*(pc_ke(2&
&0,18) - inv_2zeta_a*2*pc_ke(4,18))

pc(1,19) = PB(2)*pc(1,6) + inv_2zeta*(2*pc(1,3)) + 16*xi*inv_8m3c2*(pc_ke(1,19&
&) - inv_2zeta_b*2*pc_ke(1,3))

pc(2,19) = PA(1)*pc(1,19) + 16*xi*inv_8m3c2*(pc_ke(2,19))

pc(3,19) = PA(2)*pc(1,19) + inv_2zeta*(3*pc(1,6)) + 16*xi*inv_8m3c2*(pc_ke(3,1&
&9))

pc(4,19) = PA(3)*pc(1,19) + 16*xi*inv_8m3c2*(pc_ke(4,19))

pc(5,19) = PA(1)*pc(2,19) + inv_2zeta*(pc(1,19)) + 16*xi*inv_8m3c2*(pc_ke(5,19&
&) - inv_2zeta_a*pc_ke(1,19))

pc(6,19) = PA(2)*pc(3,19) + inv_2zeta*(pc(1,19) + 3*pc(3,6)) + 16*xi*inv_8m3c2&
&*(pc_ke(6,19) - inv_2zeta_a*pc_ke(1,19))

pc(7,19) = PA(3)*pc(4,19) + inv_2zeta*(pc(1,19)) + 16*xi*inv_8m3c2*(pc_ke(7,19&
&) - inv_2zeta_a*pc_ke(1,19))

pc(8,19) = PA(2)*pc(2,19) + inv_2zeta*(3*pc(2,6)) + 16*xi*inv_8m3c2*(pc_ke(8,1&
&9))

pc(9,19) = PA(3)*pc(2,19) + 16*xi*inv_8m3c2*(pc_ke(9,19))

pc(10,19) = PA(3)*pc(3,19) + 16*xi*inv_8m3c2*(pc_ke(10,19))

pc(11,19) = PA(3)*pc(8,19) + 16*xi*inv_8m3c2*(pc_ke(11,19))

pc(12,19) = PA(2)*pc(5,19) + inv_2zeta*(3*pc(5,6)) + 16*xi*inv_8m3c2*(pc_ke(12&
&,19))

pc(13,19) = PA(3)*pc(5,19) + 16*xi*inv_8m3c2*(pc_ke(13,19))

pc(14,19) = PA(1)*pc(6,19) + 16*xi*inv_8m3c2*(pc_ke(14,19))

pc(15,19) = PA(3)*pc(6,19) + 16*xi*inv_8m3c2*(pc_ke(15,19))

pc(16,19) = PA(1)*pc(7,19) + 16*xi*inv_8m3c2*(pc_ke(16,19))

pc(17,19) = PA(2)*pc(7,19) + inv_2zeta*(3*pc(7,6)) + 16*xi*inv_8m3c2*(pc_ke(17&
&,19))

pc(18,19) = PA(1)*pc(5,19) + inv_2zeta*(2*pc(2,19)) + 16*xi*inv_8m3c2*(pc_ke(1&
&8,19) - inv_2zeta_a*2*pc_ke(2,19))

pc(19,19) = PA(2)*pc(6,19) + inv_2zeta*(2*pc(3,19) + 3*pc(6,6)) + 16*xi*inv_8m&
&3c2*(pc_ke(19,19) - inv_2zeta_a*2*pc_ke(3,19))

pc(20,19) = PA(3)*pc(7,19) + inv_2zeta*(2*pc(4,19)) + 16*xi*inv_8m3c2*(pc_ke(2&
&0,19) - inv_2zeta_a*2*pc_ke(4,19))

pc(1,20) = PB(3)*pc(1,7) + inv_2zeta*(2*pc(1,4)) + 16*xi*inv_8m3c2*(pc_ke(1,20&
&) - inv_2zeta_b*2*pc_ke(1,4))

pc(2,20) = PA(1)*pc(1,20) + 16*xi*inv_8m3c2*(pc_ke(2,20))

pc(3,20) = PA(2)*pc(1,20) + 16*xi*inv_8m3c2*(pc_ke(3,20))

pc(4,20) = PA(3)*pc(1,20) + inv_2zeta*(3*pc(1,7)) + 16*xi*inv_8m3c2*(pc_ke(4,2&
&0))

pc(5,20) = PA(1)*pc(2,20) + inv_2zeta*(pc(1,20)) + 16*xi*inv_8m3c2*(pc_ke(5,20&
&) - inv_2zeta_a*pc_ke(1,20))

pc(6,20) = PA(2)*pc(3,20) + inv_2zeta*(pc(1,20)) + 16*xi*inv_8m3c2*(pc_ke(6,20&
&) - inv_2zeta_a*pc_ke(1,20))

pc(7,20) = PA(3)*pc(4,20) + inv_2zeta*(pc(1,20) + 3*pc(4,7)) + 16*xi*inv_8m3c2&
&*(pc_ke(7,20) - inv_2zeta_a*pc_ke(1,20))

pc(8,20) = PA(2)*pc(2,20) + 16*xi*inv_8m3c2*(pc_ke(8,20))

pc(9,20) = PA(3)*pc(2,20) + inv_2zeta*(3*pc(2,7)) + 16*xi*inv_8m3c2*(pc_ke(9,2&
&0))

pc(10,20) = PA(3)*pc(3,20) + inv_2zeta*(3*pc(3,7)) + 16*xi*inv_8m3c2*(pc_ke(10&
&,20))

pc(11,20) = PA(3)*pc(8,20) + inv_2zeta*(3*pc(8,7)) + 16*xi*inv_8m3c2*(pc_ke(11&
&,20))

pc(12,20) = PA(2)*pc(5,20) + 16*xi*inv_8m3c2*(pc_ke(12,20))

pc(13,20) = PA(3)*pc(5,20) + inv_2zeta*(3*pc(5,7)) + 16*xi*inv_8m3c2*(pc_ke(13&
&,20))

pc(14,20) = PA(1)*pc(6,20) + 16*xi*inv_8m3c2*(pc_ke(14,20))

pc(15,20) = PA(3)*pc(6,20) + inv_2zeta*(3*pc(6,7)) + 16*xi*inv_8m3c2*(pc_ke(15&
&,20))

pc(16,20) = PA(1)*pc(7,20) + 16*xi*inv_8m3c2*(pc_ke(16,20))

pc(17,20) = PA(2)*pc(7,20) + 16*xi*inv_8m3c2*(pc_ke(17,20))

pc(18,20) = PA(1)*pc(5,20) + inv_2zeta*(2*pc(2,20)) + 16*xi*inv_8m3c2*(pc_ke(1&
&8,20) - inv_2zeta_a*2*pc_ke(2,20))

pc(19,20) = PA(2)*pc(6,20) + inv_2zeta*(2*pc(3,20)) + 16*xi*inv_8m3c2*(pc_ke(1&
&9,20) - inv_2zeta_a*2*pc_ke(3,20))

pc(20,20) = PA(3)*pc(7,20) + inv_2zeta*(2*pc(4,20) + 3*pc(7,7)) + 16*xi*inv_8m&
&3c2*(pc_ke(20,20) - inv_2zeta_a*2*pc_ke(4,20))

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


   end subroutine massvel2CIntgAna

   subroutine massvel2CIntgNum(a1,a2,A,B,pc,sh,cell_size,step_size)

   use O_Kinds
   use O_Constants, only: pi, fineStructure

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
   integer :: p, q, h, i, j, k
   integer :: num_steps
   integer, dimension (20,3) :: triads
   integer, dimension (16,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, coeff
   real (kind=double), dimension (2) :: curr_pos
   real (kind=double), dimension (3) :: xyz, xyz_soln
   real (kind=double), dimension (3,3) :: xyz_I
      ! The second index is for prime, noprime, 2Dprime.
      ! The first index is xyz for prime and noprime while it is xy,xz,yz for
      !   the 2Dprime case.

   ! Before we proceed with the calculation we need to understand a bit more
   !   about exactly what is being computed. The form of the integration is:
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        -coeff * dell^2 . dell^2
   !        [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   ! We use Pxyz for points in space; Rxyz1,2 for atomic sites 1 and 2;
   !   zeta1,2 for the exponential decay factors; lxyz1,2 for the angular
   !   momentum numbers for atoms 1 and 2; and coeff for -1/8m^3c^2 in atomic
   !   units.
   ! A few notes about the coefficient: Recall that p_x-hat = -ih_bar d/dx.
   !   Thus: p_x^2 = -h_bar^2 d^2/dx^2 and p_x^4 = h_bar^4 d^4/dx^4.
   !   In the simplified Dirac equation (Eqn. 4.12 of the OLCAO book or Eqn.
   !   3 of the Zhong, Xu, Ching PRB v41 p10545 1990 paper) the coefficient
   !   of the kinetic energy is -1/2m with m = 1 in atomic units because
   !   p^2 = -h_bar^2 d^2/dx^2 and h_bar = 1 in atomic units. Similarly, the
   !   coefficient for the mass velocity term is -h_bar^4 / (8 m^3 c^2) which
   !   (using h_bar = m = 1 in atomic units and c = 1/(fine structure const.
   !   expressed in atomic units)) we get a coefficient of:
   !   - (fine structure in a.u.)^2 / 8. The negative sign will be explicit
   !   and the coefficient variable will only be the magnitude.
   !
   ! The dell^2 operator is d^2/dx^2 + d^2/dy^2 + d^2/dz^2 so that dell^2
   !   times dell^2 is [d^4/dx^4 + d^4/dy^4 + d^4/dz^4 + 2 d^2/dx^2 d^2/dy^2
   !   + 2 d^2/dx^2 d^2/dz^2 + 2 d^2/dy^2 d^2/dz^2].
   ! The integral must be computed over all space for all different possible
   !   values of lxyz1,2 for some arbitrarily chosen test coordinates and
   !   decay rates.
   !
   ! Because of the plus signs in the dell^2 times dell^2 operator we arrive
   !   at three identical integrals of the form (with # in d^4/d#^4 = x, y,
   !   or z):
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 *
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        coeff * d^4/d#^4 [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   !
   ! Also, because of the plus signs, we arrive at three other identical
   !   integrals of the form (with #,@ = x,y ; x,z ; y,z):
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        coeff * 2 * d^2/d#^2 d^2/d@^2
   !        [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   !
   ! Now, focusing on just the d^4/dx^4 version of the set of three integrals,
   !   we will pull all terms with Py and Pz out of the dx integral to get:
   ! SS { [(Py-Ry1)**ly1 * (Pz-Rz1)**lz1 *
   !       exp(-zeta1*((Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !      [(Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !       exp(-zeta2*((Py-Ry2)**2 + (Pz-Rz2)**2 ))] *
   !     S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * coeff * d^4/dx^4 [
   !      (Px-Rx2)**lx2 * exp(-zeta2*((Px-Rx2)**2))]] dx
   !    } dydz
   ! Applying the fourth derivative, the internal 1D dx integral has the form:
   !   Ix' = S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * coeff *
   !          [(lx2^4 - 6lx2^3 + 11lx2^2 - 6lx2) * (Px-Rx2)^(lx2-4) +
   !           (zeta2*lx2*(4*lx2*(-2lx2 + 3) - 4)) * (Px-Rx2)^(lx2-2) +
   !           (zeta2^2*2*(12*lx2*(lx2 + 1) + 1)) * (Px-Rx2)^(lx2) +
   !           (zeta2^3*16*(-3 - 2*lx2)) * (Px-Rx2)^(lx2+2) +
   !           (zeta2^4*16)*(Px-Rx2)^(lx2+4)] * exp(-zeta2*(Px-Rx2)**2)
   !        ]
   ! Each of the other integrals (Iy, Iz) will have the form of a simple
   !   1D overlap integral:
   !   Iy = S [(Py-Ry1)^ly1 * exp(-zeta1*((Py-Ry1)^2)) *
   !            (Py-Ry2)^ly2 * exp(-zeta2*((Py-Ry2)^2))] dy
   !   Iz = S [(Pz-Rz1)^lz1 * exp(-zeta1*((Pz-Rz1)^2)) *
   !            (Pz-Rz2)^lz2 * exp(-zeta2*((Pz-Rz2)^2))] dz
   ! The total integral !! for the d^4/dx^4 portion !! is thus Ix' * Iy * Iz.
   !   The total of the first part of the integral (including all terms of
   !   dell^4) will have the form: Ix'*Iy*Iz + Ix*Iy'*Iz + Ix*Iy*Iz' where
   !   the Iy' and Iz' are the appropriate analogs of the Ix' and the Ix is
   !   the analog of the Iy or Iz.
   !
   ! Now, focusing on just the 2 * d^2/d#^2 d^2/d@^2 where (say) #,@ = x,y,
   !   we will pull terms with Pz out of the dx dy integral to get:
   ! S { [(Pz-Rz1)^lz1 * exp(-zeta1*(Pz-Rz1)^2)] *
   !     [(Pz-Rz2)^lz2 * exp(-zeta2*(Pz-Rz2)^2)] *
   !   SS [(Px-Rx1)^lx1 * exp(-zeta1*((Px-Rx1)^2 + (Py-Ry1)^2)) * coeff *
   !      d^2/dx^2 d^2/dy^2 [(Px-Rx2)^lx2 * (Py-Ry2)^ly2 *
   !      exp(-zeta2*((Px-Rx2)^2 + (Py-Ry2)**2))]] dx dy
   !   } dz
   ! Applying the two second derivatives, the internal 2D dx dy integral has
   !   the form:
   ! Ix'y' = S [(Px-Rx1)^lx1 * (Py-Ry1)^ly1 *
   !    exp(-zeta1*((Px-Rx1)^2 + (Py-Ry1)^2)) coeff *
   !    lx2*ly2*(lx2*ly2 - lx2 - ly2 + 1) (Px-Rx2)^(lx2-2)*(Py-Ry2)^(ly2-2) *
   !    zeta2 2lx2(-2 lx2 ly2 - lx2 + 2ly2 + 1) (Px-Rx2)^(lx2-2)*(Py-Ry2)^ly2 *
   !    zeta2 2ly2(-2 lx2 ly2 + 2lx2 - ly2 + 1) (Px-Rx2)^lx2*(Py-Ry2)^(ly2-2) *
   !    zeta2^2 4(2(2 lx2 ly2 + lx2 + ly2) + 1) (Px-Rx2)^lx2 * (Py-Ry2)^ly2 *
   !    zeta2^2 4lx2(lx2 - 1) * (Px-Rx2)^(lx2-2)*(Py-Ry2)^(ly2+2) *
   !    zeta2^2 4ly2(ly2 - 1) * (Px-Rx2)^(lx2+2)*(Py-Ry2)^(ly2-2) *
   !    zeta2^3 8(-2lx2 - 1) * (Px-Rx2)^lx2*(Py-Ry2)^(ly2+2) *
   !    zeta2^3 8(-2ly2 - 1) * (Px-Rx2)^(lx2+2)*(Py-Ry2)^ly2 *
   !    zeta2^4 16 * (Px-Rx2)^(lx2+2)*(Py-Ry2)^(ly2+2) *
   !    exp(-zeta2*((Px-Rx2)^2 + (Py-Ry2)^2))]
   !
   ! Thus, every total integral is a sum of various independent integrals. We
   !   can do 1D overlap integrals for x, y, z; 1D fourth derivative integrals
   !   for x, y, z; and 2D second derivative integrals for xy, xz, yz pairs to
   !   construct the total integral.
   ! With regard to the term in each integral that has the **(lx2-2) form
   !   (or **(ly2-2) or **(lz2-2)), some special care must be taken. This
   !   term represents an angular momentum and the integral will have the
   !   form of an overlap integral. Because we cannot have a negative angular
   !   momentum we must discard this term when lx2-2 < 0, etc.

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

   coeff = (fineStructure * 0.001d0)**2 / 8.0d0
   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   ! Initialize a counter of the triad pq pairs.
   h = 0

   do p = 1, 20
      do q = 1, 20

         ! Assign l1 and l2 values for each primitive Cartesian gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize sum variables.
         xyz_I(:,:) = 0.0d0

         ! Start a loop over 1D coordinates.
         do i = 0, num_steps
            curr_pos(1) = start_pos + (i*step_size)

            ! Compute the prime integrals first. (Fourth derivative)
            do j = 1, 3
               xyz_I(j,1) = xyz_I(j,1) &
                     & + primeMV(step_size, curr_pos(1), A(j), &
                     & B(j), a1, a2, l1(j), l2(j))
            enddo

            ! Compute the no-prime integrals second. (No derivative)
            do j = 1, 3
               xyz_I(j,2) = xyz_I(j,2) &
                     & + noPrimeMV(step_size, curr_pos(1), A(j), &
                     & B(j), a1, a2, l1(j), l2(j))
            enddo

            ! Start a loop over a second set of 1D coordiantes.
            do j = 0, num_steps
                curr_pos(2) = start_pos + (j*step_size)

                ! Compute the two-prime integrals.

                ! Compute the xy term.
                xyz_I(1,3) = xyz_I(1,3) &
                      & + twoPrimeMV(step_size, curr_pos(1), curr_pos(2), &
                      & A(1), A(2), B(1), B(2), a1, a2, l1(1), l1(2), l2(1), &
                      & l2(2))

                ! Compute the xz term.
                xyz_I(2,3) = xyz_I(2,3) &
                      & + twoPrimeMV(step_size, curr_pos(1), curr_pos(2), &
                      & A(1), A(3), B(1), B(3), a1, a2, l1(1), l1(3), l2(1), &
                      & l2(3))

                ! Compute the yz term.
                xyz_I(3,3) = xyz_I(3,3) &
                      & + twoPrimeMV(step_size, curr_pos(1), curr_pos(2), &
                      & A(2), A(3), B(2), B(3), a1, a2, l1(2), l1(3), l2(2), &
                      & l2(3))
            enddo
         enddo

         ! Assemble the integrals associated with the fourth derivative.
         pc(q,p) = &
               & xyz_I(1,1)*xyz_I(2,2)*xyz_I(3,2) + &  ! Fourth D wrt x only
               & xyz_I(1,2)*xyz_I(2,1)*xyz_I(3,2) + &  ! Fourth D wrt y only
               & xyz_I(1,2)*xyz_I(2,2)*xyz_I(3,1)      ! Fourth D wrt z only

         ! Add on the integrals associated with the two second derivatives.
         pc(q,p) = pc(q,p) + 2.0 * ( &
               & xyz_I(1,3)*xyz_I(3,2) + & ! 2nd D wrt x,y; no D wrt z
               & xyz_I(2,3)*xyz_I(2,2) + & ! 2nd D wrt x,z; no D wrt y
               & xyz_I(3,3)*xyz_I(1,2))    ! 2nd D wrt y,z; no D wrt x

         ! Multiply by the mass velocity term's scaling coefficient.
         pc(q,p) = pc(q,p) * coeff

         ! Record progress thorugh the pq pairs.
         h = h + 1
         if (mod(h,10) .eq. 0) then
            write (6,ADVANCE="NO",FMT="(a1)") "|"
         else
            write (6,ADVANCE="NO",FMT="(a1)") "."
         endif
         if (mod(h,50) .eq. 0) then
            write (6,*) " ",h
         endif
         call flush (6)

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


   end subroutine massvel2CIntgNum


   function noPrimeMV(step_size, curr_pos, A, B, a1, a2, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos, A, B, a1, a2
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: noPrimeMV

      ! Compute the integral part.
      noPrimeMV = step_size * (curr_pos - A)**l1 * (curr_pos - B)**l2 &
            & * exp(-a1*(curr_pos - A)**2) * exp(-a2*(curr_pos - B)**2)

      return

   end function noPrimeMV


   function primeMV(step_size, curr_pos, A, B, a1, a2, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos, A, B, a1, a2
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: primeMV

      ! Compute each "internal" term of the prime integral. Compare each
      !   of these lines with the 1D mass velocity equation produced by the
      !   osrecurintg_makenum.py script (appropriately separated into terms).
      !   Conceptually, these are similar to the expression in the last
      !   equals of equation 65 in Ben Walker's dissertation for the regular
      !   kinetic energy. Note that a factor of XXX is included in all of the
      !   below expressions. The most succinctly refactored equations are in
      !   ...

      primeMV = a2**4 * 16 * (curr_pos - B)**(l2+4)
      
      primeMV = primeMV + a2**3 * (-32*l2 - 48) * (curr_pos - B)**(l2+2)
      
      primeMV = primeMV + a2**2 * (24*l2**2 + 24*l2 + 12) * (curr_pos - B)**l2
      
      if (l2 >= 2) then
         primeMV = primeMV + a2 * (-8*l2**3 + 12*l2**2 - 4*l2) * &
               & (curr_pos - B)**(l2-2)
      endif
      
      if (l2 >= 4) then
         primeMV = primeMV + (l2**4 - 6*l2**3 + 11*l2**2 - 6*l2) * &
               & (curr_pos - B)**(l2-4)
      endif

      !primeMV = 16 * a2**4 * (curr_pos - B)**(l2+4)

      !primeMV = primeMV + 16 * a2**3 * (-3 - 2*l2) * (curr_pos - B)**(l2+2)

      !primeMV = primeMV + 2 * a2**2 * (12*l2 * (l2 + 1) + 1) * &
      !      & (curr_pos - B)**l2

      !if (l2 >= 2) then
      !   primeMV = primeMV + a2 * l2 * (4*l2 * (-2*l2 + 3) - 4) * &
      !         & (curr_pos - B)**(l2-2)
      !endif

      !if (l2 >= 4) then
      !   primeMV = primeMV + (l2**4 - 6*l2**3 + 11*l2**2 - 6*l2) * &
      !         & (curr_pos - B)**(l2-4)
      !endif

      ! Multiply the prime integral by the preceeding primitive gaussian
      !   coefficient and exponential and multiply by the succeeding
      !   exponential. (We have already multiplied by the succeeding
      !   primitive gaussian coefficient in the above lines.)
      primeMV = primeMV * (curr_pos-A)**l1 * exp(-a1*(curr_pos-A)**2) &
            & * exp(-a2*(curr_pos-B)**2)

      ! Finally, multiply by the step size.
      primeMV = primeMV * step_size

      return
      
   end function primeMV


   ! The notation may be confusing so here is some clarification:
   ! As in the previous subroutines, the 1 and 2 for the alphas (a's)
   !   correspond to the first and second orbitals of the integral. I.e.,
   !   they are independent of x,y,z.
   ! The 1 and 2 for A corresponds to the position of site A with respect to
   !   the first xyz coordinate axis and second xyz coordinate axis
   !   respectively. So, when we are doing the xy term, then 1=x and 2=y.
   !   When we are doing the xz term, then 1=x and 2=z. Then, for the yz
   !   term 1=y and 2=z.
   ! The same concept described above for A applies to B.
   ! For the l variables the first number (1 or 2) corresponds to the first
   !   and second orbitals of the integrals. The second number (1 or 2)
   !   corresponds to the first xyz or second xyz coordinate axes as for the
   !   A1 A2 and B1 B2.

   function twoPrimeMV(step_size, curr_pos1, curr_pos2, A1, A2, B1, B2, &
         & alpha1, alpha2, l11, l12, l21, l22)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos1, curr_pos2
      real (kind=double), intent(in) :: A1, A2, B1, B2
      real (kind=double), intent(in) :: alpha1, alpha2
      integer, intent(in) :: l11, l12, l21, l22

      ! Define local and return variables.
      real (kind=double) :: twoPrimeMV

      ! Compute each internal term of the twoPrime integral. As with the prime
      !   integral, these terms are obtained from sympy derivations and simple
      !   algebraic manipulations.

      twoPrimeMV = alpha2**4*16 &
            & * (curr_pos1 - B1)**(l21 + 2) * (curr_pos2 - B2)**(l22 + 2)

      twoPrimeMV = twoPrimeMV + alpha2**3*(-8 - 16*l22) &
            & * (curr_pos1 - B1)**(l21 + 2) * (curr_pos2 - B2)**(l22)

      twoPrimeMV = twoPrimeMV + alpha2**3*(-16*l21 - 8) &
            & * (curr_pos1 - B1)**(l21) * (curr_pos2 - B2)**(l22 + 2)

      twoPrimeMV = twoPrimeMV + alpha2**2*(16*l21*l22 + 8*l21 + 8*l22 + 4) &
            & * (curr_pos1 - B1)**(l21) * (curr_pos2 - B2)**(l22)

      if (l21 >= 2) then
         twoPrimeMV = twoPrimeMV + alpha2*(-4*l21**2*l22 - 2*l21**2 &
               & + 4*l21*l22 + 2*l21) &
               & * (curr_pos1 - B1)**(l21-2) * (curr_pos2 - B2)**(l22)

         twoPrimeMV = twoPrimeMV + alpha2**2*(4*l21**2 - 4*l21) &
               & * (curr_pos1 - B1)**(l21-2) * (curr_pos2 - B2)**(l22+2)
      endif

      if (l22 >= 2) then
         twoPrimeMV = twoPrimeMV + alpha2*(-4*l21*l22**2 + 4*l21*l22 &
               & - 2*l22**2 + 2*l22) &
               & * (curr_pos1 - B1)**(l21) * (curr_pos2 - B2)**(l22-2)

         twoPrimeMV = twoPrimeMV + alpha2**2*(4*l22**2 - 4*l22) &
               & * (curr_pos1 - B1)**(l21+2) * (curr_pos2 - B2)**(l22-2)
      endif

      if ((l21 >= 2) .and. (l22 >= 2)) then
         twoPrimeMV = twoPrimeMV + (l21**2*l22**2 - l21**2*l22 &
               & - l21*l22**2 + l21*l22) &
               & * (curr_pos1 - B1)**(l21 - 2) * (curr_pos2 - B2)**(l22 - 2)
      endif


      !twoPrimeMV = 16*alpha2**4 &
      !      & * (curr_pos1 - B1)**(l21 + 2) * (curr_pos2 - B2)**(l22 + 2)

      !twoPrimeMV = twoPrimeMV + 8*alpha2**3 * (-1 - 2*l22) &
      !      & * (curr_pos1 - B1)**(l21 + 2) * (curr_pos2 - B2)**(l22)

      !twoPrimeMV = twoPrimeMV + 8*alpha2**3 * (-2*l21 - 1) &
      !      & * (curr_pos1 - B1)**(l21) * (curr_pos2 - B2)**(l22 + 2)

      !twoPrimeMV = twoPrimeMV + 4*alpha2**2 &
      !      & * (2*(2*l21*l22 + l21 + l22) + 1) &
      !      & * (curr_pos1 - B1)**(l21) * (curr_pos2 - B2)**(l22)

      !if (l22 >= 2) then
      !   twoPrimeMV = twoPrimeMV + 4*alpha2**2 * l22*(l22 - 1) &
      !         & * (curr_pos1 - B1)**(l21 + 2) * (curr_pos2 - B2)**(l22 - 2)

      !   twoPrimeMV = twoPrimeMV + 2*alpha2 &
      !         & * l22*(-2*l21*l22 + 2*l21 - l22 + 1) &
      !         & * (curr_pos1 - B1)**(l21) * (curr_pos2 - B2)**(l22 - 2)
      !endif

      !if (l21 >= 2) then
      !   twoPrimeMV = twoPrimeMV + 4*alpha2**2 * l21*(l21 - 1) &
      !         & * (curr_pos1 - B1)**(l21 - 2) * (curr_pos2 - B2)**(l22 + 2)

      !   twoPrimeMV = twoPrimeMV + 2*alpha2 &
      !         & * l21*(-2*l21*l22 - l21 + 2*l22 + 1) &
      !         & * (curr_pos1 - B1)**(l21 - 2) * (curr_pos2 - B2)**(l22)
      !endif

      !if ((l21 >= 2) .and. (l22 >= 2)) then
      !   twoPrimeMV = twoPrimeMV + l21*l22*(l21*l22 - l21 - l22 + 1) &
      !         & * (curr_pos1 - B1)**(l21 - 2) * (curr_pos2 - B2)**(l22 - 2)
      !endif

      ! Multiply the integral by the preceeding primitive gaussian
      !   coefficient and exponential and multiply by the succeeding
      !   exponential. (We have already multiplied by the succeeding
      !   primitive gaussian coefficient in the above lines.)
      twoPrimeMV = twoPrimeMV * (curr_pos1-A1)**l11 * (curr_pos2-A2)**l12 * &
            & exp(-alpha1*((curr_pos1-A1)**2 + (curr_pos2-A2)**2)) * &
            & exp(-alpha2*((curr_pos1-B1)**2 + (curr_pos2-B2)**2))

      ! Finally, multiply by the step size squared (2D integral).
      twoPrimeMV = twoPrimeMV * step_size * step_size

      return

   end function twoPrimeMV

end program GaussianIntegrals
