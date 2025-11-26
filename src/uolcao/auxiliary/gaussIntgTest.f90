program GaussianIntegrals

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define program variables.
   real (kind=double), allocatable, dimension (:,:,:,:) :: pc
   real (kind=double), allocatable, dimension (:,:,:,:) :: sh
   complex (kind=double), allocatable, dimension (:,:,:,:) :: pcCmplx
   complex (kind=double), allocatable, dimension (:,:,:,:) :: shCmplx
   real (kind=double), dimension (3) :: alphas
   real (kind=double), dimension (3) :: deltaK ! x,y,z of Kf-Ki
   real (kind=double), dimension (3,3) :: pos ! x,y,z of A,B,C
   real (kind=double), dimension (3,2) :: temp_alphas
   real (kind=double), dimension (3,2) :: temp_deltaK ! x,y,z of Kf-Ki
   real (kind=double), dimension (3,3,2) :: temp_pos ! x,y,z of A,B,C
   real (kind=double), dimension (3) :: temp_alphas_step
   real (kind=double), dimension (3) :: temp_deltaK_step ! x,y,z of Kf-Ki
   real (kind=double), dimension (3,3) :: temp_pos_step ! x,y,z of A,B,C
   real (kind=double) :: cell_size, step_size ! Parameters for the numerical
         ! integration process. Read carefully! The cell size is given as a
         ! positive number that defines the maximum extent of the space to
         ! integrate along any given axis. Hence the actual cell is double
         ! that number. I.e., integrate along x, y, z axes from -cell_size
         ! to +cell_size. The step size is just what it appears to be, but
         ! note that it is for numerical integration and is a very different
         ! concept and is unrelated to the "num_steps" defined below.
   integer :: num_segments, num_steps, curr_line
   integer :: h, i

   ! Allocate space to hold the appropriately sized pc and sh matrices. The
   !   last index in both pc and sh is a 2 to hold analytical solutions
   !   (value = 1) and numerical solutions (value = 2). The second to last
   !   index is the max of three or (2*max_lam+1) to allow for the momentum
   !   matrix and multi-level nuclear solutions (but all elements are not
   !   always used). The pc matrix is used multiple times to produce the three
   !   momentum matrix sets (x,y,z).
   allocate (pc(20,20,7,2))
   allocate (sh(16,16,3,2))
   allocate (pcCmplx(20,20,7,2))
   allocate (shCmplx(16,16,3,2))

   ! Open the control file.
   open (10, file="intgcontrol", status='old')

   ! Read the cell size and the grid step size for numerical integration.
   read (10,*) cell_size, step_size

   ! Read the number of segments.
   read (10,*) num_segments

   ! Initialize the line number for printing.
   curr_line = 0

   ! For each segment, read in the number of steps, allocate space to hold
   !   relevant information for each step, compute the step descriptions,
   !   compute the integral matrices, store the results, deallocate for the
   !   next segment.
   do h = 1, num_segments

      ! Read in the number of steps for this segment.
      read (10,*) num_steps
      if (num_steps == 1) then
   
         ! Read the parameters for the current segment.
         read (10,*) temp_alphas(:,1), temp_pos(:,:,1), temp_deltaK(:,1)

         ! Fill in dummy values for the temp step data. (It will not be used.)
         temp_alphas_step(:) = 0.0d0
         temp_deltaK_step(:) = 0.0d0
         temp_pos_step(:,:) = 0.0d0
      else
   
         ! Read the beginning and ending parameters for the current segment.
         !   First step followed by last step.
         read (10,*) temp_alphas(:,1), temp_pos(:,:,1), temp_deltaK(:,1)
         read (10,*) temp_alphas(:,2), temp_pos(:,:,2), temp_deltaK(:,2)

         ! Compute the size of the step for each parameter.
         temp_alphas_step(:) = (temp_alphas(:,2) - temp_alphas(:,1))/num_steps
         temp_deltaK_step(:) = (temp_deltaK(:,2) - temp_deltaK(:,1))/num_steps
         temp_pos_step(:,:) = (temp_pos(:,:,2) - temp_pos(:,:,1))/num_steps
      endif

      ! Start iterating over the steps.
      do i = 1, num_steps

         ! Increment the line number
         curr_line = curr_line + 1

         alphas(:) = temp_alphas(:,1) + (i-1) * temp_alphas_step(:)
         deltaK(:) = temp_deltaK(:,1) + (i-1) * temp_deltaK_step(:)
         pos(:,:) = temp_pos(:,:,1) + (i-1) * temp_pos_step(:,:)


         ! Compute the pc and sh integral results for the current parameters
         !   using the analytic formulas.
         call Koverlap2CIntgAna(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & deltaK(:),pcCmplx(:,:,1,1),shCmplx(:,:,1,1))

         ! Compute the pc and sh integral results for the current parameters
         !   using numerical integration.
         call Koverlap2CIntgNum(alphas(1),alphas(2),pos(:,1),pos(:,2),&
               & deltaK(:),pcCmplx(:,:,1,2),shCmplx(:,:,1,2),&
               & cell_size,step_size)

         ! Print the pc and sh integral result differences.
         call print_pc_sh(h,i,11,alphas,pos,deltaK,real(pcCmplx(:,:,1,:),&
               & double),real(shCmplx(:,:,1,:),double),"KO_Real.dat")
         call print_pc_sh(h,i,12,alphas,pos,deltaK,aimag(pcCmplx(:,:,1,:)),&
               & aimag(shCmplx(:,:,1,:)),"KO_Cmpx.dat")


      enddo
   enddo

   ! Deallocate space.
   deallocate(pc)
   deallocate(sh)
   deallocate(pcCmplx)
   deallocate(shCmplx)


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module/program subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

   subroutine print_pc_sh(h,i,unitPlus,alphas,pos,deltaK,pc,sh,filename)

      ! Use necessary modules
      use O_Kinds

      ! Make sure no funny variables are created.
      implicit none

      ! Define dummy variables.
      integer :: h, i, unitPlus
      real(kind=double), dimension(3) :: alphas
      real(kind=double), dimension(3,3) :: pos ! xyz (idx1) of atoms ABC (2nd)
      real(kind=double), dimension(3) :: deltaK
      real(kind=double), dimension(20,20,2) :: pc ! Last idx 1 = ana; 2 = num
      real(kind=double), dimension(16,16,2) :: sh ! Last idx 1 = ana; 2 = num
      character*11 :: filename
      character*14, dimension(16) :: shName
      character*3, dimension(20) :: pcName

      ! Define local variables.
      logical :: io_opened
      integer :: j,k

      ! Initialize the shName and pcName arrays.
      shName(1) = "s"
      shName(2) = "x"
      shName(3) = "y"
      shName(4) = "z"
      shName(5) = "xy"
      shName(6) = "xz"
      shName(7) = "yz"
      shName(8) = "xx-yy"
      shName(9) = "2zz-xx-yy"
      shName(10) = "xyz"
      shName(11) = "xxz-yyz"
      shName(12) = "xxx-3yyx"
      shName(13) = "3xxy-yyy"
      shName(14) = "2zzz-3xxz-3yyz"
      shName(15) = "4zzx-xxx-yyx"
      shName(16) = "4zzy-xxy-yyy"
      pcName(1) = "s"
      pcName(2) = "x"
      pcName(3) = "y"
      pcName(4) = "z"
      pcName(5) = "xx"
      pcName(6) = "yy"
      pcName(7) = "zz"
      pcName(8) = "xy"
      pcName(9) = "xz"
      pcName(10) = "yz"
      pcName(11) = "xyz"
      pcName(12) = "xxy"
      pcName(13) = "xxz"
      pcName(14) = "yyx"
      pcName(15) = "yyz"
      pcName(16) = "zzx"
      pcName(17) = "zzy"
      pcName(18) = "xxx"
      pcName(19) = "yyy"
      pcName(20) = "zzz"

      ! Open the output file if it isn't already open.
      inquire(299+unitPlus, OPENED = io_opened)
      if (.not. io_opened) then
         open (299+unitPlus, file=filename, status="unknown")

         ! Print a header line.

         ! Start with the parameters.
         write (299+unitPlus,fmt="(a)",advance="NO") "IDX SEG STEP "
         write (299+unitPlus,fmt="(a)",advance="NO") "alpha1 alpha2 alpha3 "
         write (299+unitPlus,fmt="(a)",advance="NO") "Ax Ay Az "
         write (299+unitPlus,fmt="(a)",advance="NO") "Bx By Bz "
         write (299+unitPlus,fmt="(a)",advance="NO") "Cx Cy Cz "
         write (299+unitPlus,fmt="(a)",advance="NO") "deltaKx deltaKy deltaKz "

         ! Now add the pc names.
         do j = 1, 20
            do k = 1, 20
               write (299+unitPlus,fmt="(4a)",advance="NO") &
                     & trim(pcName(k)),"_",trim(pcName(j))," "
            enddo
         enddo

         ! Now add the sh names.
         do j = 1, 16
            do k = 1, 16
               write (299+unitPlus,fmt="(4a)",advance="NO") &
                     & trim(shName(k)),"_",trim(shName(j))," "
            enddo
         enddo

      endif

      if ((num_segments > 1) .or. (num_steps > 1)) then

         ! Print the line number.
         write (299+unitPlus,fmt="(i5)",advance="NO") curr_line

         write (299+unitPlus,fmt="(2i5)",advance="NO") h, i
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") alphas(:)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,1)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,2)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,3)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") deltaK(:)

         ! Print the result difference for the given pc matrix.
         do j = 1, 20
            do k = 1, 20
               !write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,1)
               !write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,2)
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

      else
         write (299+unitPlus,fmt="(2i5)",advance="NO") h, i
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") alphas(:)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,1)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,2)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,3)
         write (299+unitPlus,fmt="(3e16.8)") deltaK(:)

         ! Print the result difference for the given pc matrix.
         do j = 1, 20
            do k = 1, 20
               write (299+unitPlus,fmt="(2a4)") pcName(k), pcName(j)
               write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,1)
               write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,2)
               write (299+unitPlus,fmt="(e16.8)") &
                     & pc(k,j,1) - pc(k,j,2)
            enddo
         enddo

         ! Print the result difference for the given sh matrix.
         do j = 1, 16
            do k = 1, 16
               write (299+unitPlus,fmt="(2a16)") shName(k), shName(j)
               write (299+unitPlus,fmt="(e16.8)") &
                     & sh(k,j,1) - sh(k,j,2)
            enddo
         enddo
      endif
   end subroutine

   subroutine Koverlap2CIntgAna(a1,a2,A,B,deltaK,pc,sh)

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
   real (kind=double), dimension (3) :: deltaK
   complex (kind=double), dimension (20,20), intent(out) :: pc
   complex (kind=double), dimension (16,16), intent(out) :: sh

   ! Define local variables.
   real (kind=double), dimension (3) :: P, PA, PB, d
   !real (kind=double), dimension (3) :: deltaK
   real (kind=double) :: zeta, inv_2zeta, xi
   complex (kind=double), dimension (3) :: preFactorKO
   real (kind=double), dimension (3) :: hermite_r
   real (kind=double), dimension (6,3) :: Hn
   complex (kind=double), dimension (6,3) :: hermite_term

   ! Initialize local variables.
   zeta = a1 + a2
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   xi = a1 * a2 / zeta
   P(:) = (a1*A(:) + a2*B(:)) / zeta
   PA(:) = P(:) - A(:)
   PB(:) = P(:) - B(:)
   d(:) = A(:) - B(:)
   
   ! Solving Hermite Polynomials
   hermite_r(:) = -deltaK(:)/(2*(zeta)**0.5)


   ! Calculate Hermite based on the value of compined_l
   !Hn(0,:) = 1.0
   Hn(1,:) = 2.0*hermite_r(:)
   Hn(2,:) = 4.0*hermite_r(:)**2 - 2.0
   Hn(3,:) = 8.0*hermite_r(:)**3 - 12.0*hermite_r(:)
   Hn(4,:) = 16.0*hermite_r(:)**4 - 48.0*hermite_r(:)**2 + 12.0
   Hn(5,:) = 32.0*hermite_r(:)**5 - 160.0*hermite_r(:)**3 + 120.0*hermite_r(:)
   Hn(6,:) = 64.0*hermite_r(:)**6 - 480.0*hermite_r(:)**4 &
         & + 720.0*hermite_r(:)**2 - 120.0

   !hermite_term(0,:) = Hn(1,:) ! Combined_l = 0
   hermite_term(1,:) = Hn(1,:) * &
         & (((cmplx(0.0d0,1.0d0,double))/(2.0d0*((zeta)**0.5)))**1)
   hermite_term(2,:) = Hn(2,:) * &
         & (((cmplx(0.0d0,1.0d0,double))/(2.0d0*((zeta)**0.5)))**2)
   hermite_term(3,:) = Hn(3,:) * &
         & (((cmplx(0.0d0,1.0d0,double))/(2.0d0*((zeta)**0.5)))**3)
   hermite_term(4,:) = Hn(4,:) * &
         & (((cmplx(0.0d0,1.0d0,double))/(2.0d0*((zeta)**0.5)))**4)
   hermite_term(5,:) = Hn(5,:) * &
         & (((cmplx(0.0d0,1.0d0,double))/(2.0d0*((zeta)**0.5)))**5)
   hermite_term(6,:) = Hn(6,:) * &
         & (((cmplx(0.0d0,1.0d0,double))/(2.0d0*((zeta)**0.5)))**6)

   ! This is the (s|K|s) integral
   preFactorKO(:) = ((pi/zeta)**0.5)*exp(-xi*d(:)*d(:))&
           & *exp(-(deltaK(:))**2/(4*zeta)) &
           & *exp(cmplx(0.0d0,-1.0d0,double)*deltaK(:)*P(:))

pc(1,1) = preFactorKO(1)*preFactorKO(2)*preFactorKO(3)

pc(2,1) = (PA(1) + hermite_term(1,1))*preFactorKO(1)*preFactorKO(2)*preFactorK&
&O(3)

pc(3,1) = (PA(2) + hermite_term(1,2))*preFactorKO(1)*preFactorKO(2)*preFactorK&
&O(3)

pc(4,1) = (PA(3) + hermite_term(1,3))*preFactorKO(1)*preFactorKO(2)*preFactorK&
&O(3)

pc(5,1) = (PA(1) + hermite_term(1,1))*pc(2,1) + (hermite_term(2,1) - hermite_t&
&erm(1,1)**2)*pc(1,1)

pc(6,1) = (PA(2) + hermite_term(1,2))*pc(3,1) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(1,1)

pc(7,1) = (PA(3) + hermite_term(1,3))*pc(4,1) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(1,1)

pc(8,1) = (PA(2) + hermite_term(1,2))*pc(2,1)

pc(9,1) = (PA(3) + hermite_term(1,3))*pc(2,1)

pc(10,1) = (PA(3) + hermite_term(1,3))*pc(3,1)

pc(11,1) = (PA(3) + hermite_term(1,3))*pc(8,1)

pc(12,1) = (PA(2) + hermite_term(1,2))*pc(5,1)

pc(13,1) = (PA(3) + hermite_term(1,3))*pc(5,1)

pc(14,1) = (PA(1) + hermite_term(1,1))*pc(6,1)

pc(15,1) = (PA(3) + hermite_term(1,3))*pc(6,1)

pc(16,1) = (PA(1) + hermite_term(1,1))*pc(7,1)

pc(17,1) = (PA(2) + hermite_term(1,2))*pc(7,1)

pc(18,1) = (PA(1) + hermite_term(1,1))*pc(5,1) + ((2.0*PA(1))*(hermite_term(2,&
&1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_te&
&rm(1,1)))*pc(1,1)

pc(19,1) = (PA(2) + hermite_term(1,2))*pc(6,1) + ((2.0*PA(2))*(hermite_term(2,&
&2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_te&
&rm(1,2)))*pc(1,1)

pc(20,1) = (PA(3) + hermite_term(1,3))*pc(7,1) + ((2.0*PA(3))*(hermite_term(2,&
&3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_te&
&rm(1,3)))*pc(1,1)

pc(1,2) = (PB(1) + hermite_term(1,1))*preFactorKO(1)*preFactorKO(2)*preFactorK&
&O(3)

pc(2,2) = (PA(1) + hermite_term(1,1))*pc(1,2) + (hermite_term(2,1) - hermite_t&
&erm(1,1)**2)*pc(1,1)

pc(3,2) = (PA(2) + hermite_term(1,2))*pc(1,2)

pc(4,2) = (PA(3) + hermite_term(1,3))*pc(1,2)

pc(5,2) = (PA(1) + hermite_term(1,1))*pc(2,2) + ((PA(1)+PB(1))*(hermite_term(2&
&,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_t&
&erm(1,1)))*pc(1,1)

pc(6,2) = (PA(2) + hermite_term(1,2))*pc(3,2) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(1,2)

pc(7,2) = (PA(3) + hermite_term(1,3))*pc(4,2) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(1,2)

pc(8,2) = (PA(2) + hermite_term(1,2))*pc(2,2)

pc(9,2) = (PA(3) + hermite_term(1,3))*pc(2,2)

pc(10,2) = (PA(3) + hermite_term(1,3))*pc(3,2)

pc(11,2) = (PA(3) + hermite_term(1,3))*pc(8,2)

pc(12,2) = (PA(2) + hermite_term(1,2))*pc(5,2)

pc(13,2) = (PA(3) + hermite_term(1,3))*pc(5,2)

pc(14,2) = (PA(1) + hermite_term(1,1))*pc(6,2) + (hermite_term(2,1) - hermite_&
&term(1,1)**2)*pc(6,1)

pc(15,2) = (PA(3) + hermite_term(1,3))*pc(6,2)

pc(16,2) = (PA(1) + hermite_term(1,1))*pc(7,2) + (hermite_term(2,1) - hermite_&
&term(1,1)**2)*pc(7,1)

pc(17,2) = (PA(2) + hermite_term(1,2))*pc(7,2)

pc(18,2) = (PA(1) + hermite_term(1,1))*pc(5,2) + ((PA(1)**2 + 2.0*PA(1)*PB(1))&
&*(hermite_term(2,1) - hermite_term(1,1)**2) + (2.0*PA(1) + PB(1))*(hermite_ter&
&m(3,1) - hermite_term(2,1)*hermite_term(1,1)) + (hermite_term(4,1) - hermite_t&
&erm(3,1)*hermite_term(1,1)))*pc(1,1)

pc(19,2) = (PA(2) + hermite_term(1,2))*pc(6,2) + ((2.0*PA(2))*(hermite_term(2,&
&2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_te&
&rm(1,2)))*pc(1,2)

pc(20,2) = (PA(3) + hermite_term(1,3))*pc(7,2) + ((2.0*PA(3))*(hermite_term(2,&
&3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_te&
&rm(1,3)))*pc(1,2)

pc(1,3) = (PB(2) + hermite_term(1,2))*preFactorKO(1)*preFactorKO(2)*preFactorK&
&O(3)

pc(2,3) = (PA(1) + hermite_term(1,1))*pc(1,3)

pc(3,3) = (PA(2) + hermite_term(1,2))*pc(1,3) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(1,1)

pc(4,3) = (PA(3) + hermite_term(1,3))*pc(1,3)

pc(5,3) = (PA(1) + hermite_term(1,1))*pc(2,3) + (hermite_term(2,1) - hermite_t&
&erm(1,1)**2)*pc(1,3)

pc(6,3) = (PA(2) + hermite_term(1,2))*pc(3,3) + ((PA(2)+PB(2))*(hermite_term(2&
&,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_t&
&erm(1,2)))*pc(1,1)

pc(7,3) = (PA(3) + hermite_term(1,3))*pc(4,3) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(1,3)

pc(8,3) = (PA(2) + hermite_term(1,2))*pc(2,3) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(2,1)

pc(9,3) = (PA(3) + hermite_term(1,3))*pc(2,3)

pc(10,3) = (PA(3) + hermite_term(1,3))*pc(3,3)

pc(11,3) = (PA(3) + hermite_term(1,3))*pc(8,3)

pc(12,3) = (PA(2) + hermite_term(1,2))*pc(5,3) + (hermite_term(2,2) - hermite_&
&term(1,2)**2)*pc(5,1)

pc(13,3) = (PA(3) + hermite_term(1,3))*pc(5,3)

pc(14,3) = (PA(1) + hermite_term(1,1))*pc(6,3)

pc(15,3) = (PA(3) + hermite_term(1,3))*pc(6,3)

pc(16,3) = (PA(1) + hermite_term(1,1))*pc(7,3)

pc(17,3) = (PA(2) + hermite_term(1,2))*pc(7,3) + (hermite_term(2,2) - hermite_&
&term(1,2)**2)*pc(7,1)

pc(18,3) = (PA(1) + hermite_term(1,1))*pc(5,3) + ((2.0*PA(1))*(hermite_term(2,&
&1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_te&
&rm(1,1)))*pc(1,3)

pc(19,3) = (PA(2) + hermite_term(1,2))*pc(6,3) + ((PA(2)**2 + 2.0*PA(2)*PB(2))&
&*(hermite_term(2,2) - hermite_term(1,2)**2) + (2.0*PA(2) + PB(2))*(hermite_ter&
&m(3,2) - hermite_term(2,2)*hermite_term(1,2)) + (hermite_term(4,2) - hermite_t&
&erm(3,2)*hermite_term(1,2)))*pc(1,1)

pc(20,3) = (PA(3) + hermite_term(1,3))*pc(7,3) + ((2.0*PA(3))*(hermite_term(2,&
&3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_te&
&rm(1,3)))*pc(1,3)

pc(1,4) = (PB(3) + hermite_term(1,3))*preFactorKO(1)*preFactorKO(2)*preFactorK&
&O(3)

pc(2,4) = (PA(1) + hermite_term(1,1))*pc(1,4)

pc(3,4) = (PA(2) + hermite_term(1,2))*pc(1,4)

pc(4,4) = (PA(3) + hermite_term(1,3))*pc(1,4) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(1,1)

pc(5,4) = (PA(1) + hermite_term(1,1))*pc(2,4) + (hermite_term(2,1) - hermite_t&
&erm(1,1)**2)*pc(1,4)

pc(6,4) = (PA(2) + hermite_term(1,2))*pc(3,4) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(1,4)

pc(7,4) = (PA(3) + hermite_term(1,3))*pc(4,4) + ((PA(3)+PB(3))*(hermite_term(2&
&,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_t&
&erm(1,3)))*pc(1,1)

pc(8,4) = (PA(2) + hermite_term(1,2))*pc(2,4)

pc(9,4) = (PA(3) + hermite_term(1,3))*pc(2,4) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(2,1)

pc(10,4) = (PA(3) + hermite_term(1,3))*pc(3,4) + (hermite_term(2,3) - hermite_&
&term(1,3)**2)*pc(3,1)

pc(11,4) = (PA(3) + hermite_term(1,3))*pc(8,4) + (hermite_term(2,3) - hermite_&
&term(1,3)**2)*pc(8,1)

pc(12,4) = (PA(2) + hermite_term(1,2))*pc(5,4)

pc(13,4) = (PA(3) + hermite_term(1,3))*pc(5,4) + (hermite_term(2,3) - hermite_&
&term(1,3)**2)*pc(5,1)

pc(14,4) = (PA(1) + hermite_term(1,1))*pc(6,4)

pc(15,4) = (PA(3) + hermite_term(1,3))*pc(6,4) + (hermite_term(2,3) - hermite_&
&term(1,3)**2)*pc(6,1)

pc(16,4) = (PA(1) + hermite_term(1,1))*pc(7,4)

pc(17,4) = (PA(2) + hermite_term(1,2))*pc(7,4)

pc(18,4) = (PA(1) + hermite_term(1,1))*pc(5,4) + ((2.0*PA(1))*(hermite_term(2,&
&1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_te&
&rm(1,1)))*pc(1,4)

pc(19,4) = (PA(2) + hermite_term(1,2))*pc(6,4) + ((2.0*PA(2))*(hermite_term(2,&
&2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_te&
&rm(1,2)))*pc(1,4)

pc(20,4) = (PA(3) + hermite_term(1,3))*pc(7,4) + ((PA(3)**2 + 2.0*PA(3)*PB(3))&
&*(hermite_term(2,3) - hermite_term(1,3)**2) + (2.0*PA(3) + PB(3))*(hermite_ter&
&m(3,3) - hermite_term(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite_t&
&erm(3,3)*hermite_term(1,3)))*pc(1,1)

pc(1,5) = (PB(1) + hermite_term(1,1))*pc(1,2) + (hermite_term(2,1) - hermite_t&
&erm(1,1)**2)*pc(1,1)

pc(2,5) = (PA(1) + hermite_term(1,1))*pc(1,5) + ((2.0*PB(1))*(hermite_term(2,1&
&) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_ter&
&m(1,1)))*pc(1,1)

pc(3,5) = (PA(2) + hermite_term(1,2))*pc(1,5)

pc(4,5) = (PA(3) + hermite_term(1,3))*pc(1,5)

pc(5,5) = (PA(1) + hermite_term(1,1))*pc(2,5) + ((PB(1)**2 + 2.0*PA(1)*PB(1))*&
&(hermite_term(2,1) - hermite_term(1,1)**2) + (2.0*PB(1) + PA(1))*(hermite_term&
&(3,1) - hermite_term(2,1)*hermite_term(1,1))+ (hermite_term(4,1)- hermite_term&
&(3,1)*hermite_term(1,1)))*pc(1,1)

pc(6,5) = (PA(2) + hermite_term(1,2))*pc(3,5) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(1,5)

pc(7,5) = (PA(3) + hermite_term(1,3))*pc(4,5) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(1,5)

pc(8,5) = (PA(2) + hermite_term(1,2))*pc(2,5)

pc(9,5) = (PA(3) + hermite_term(1,3))*pc(2,5)

pc(10,5) = (PA(3) + hermite_term(1,3))*pc(3,5)

pc(11,5) = (PA(3) + hermite_term(1,3))*pc(8,5)

pc(12,5) = (PA(2) + hermite_term(1,2))*pc(5,5)

pc(13,5) = (PA(3) + hermite_term(1,3))*pc(5,5)

pc(14,5) = (PA(1) + hermite_term(1,1))*pc(6,5) + ((2.0*PB(1))*(hermite_term(2,&
&1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_te&
&rm(1,1)))*pc(6,1)

pc(15,5) = (PA(3) + hermite_term(1,3))*pc(6,5)

pc(16,5) = (PA(1) + hermite_term(1,1))*pc(7,5) + ((2.0*PB(1))*(hermite_term(2,&
&1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_te&
&rm(1,1)))*pc(7,1)

pc(17,5) = (PA(2) + hermite_term(1,2))*pc(7,5)

pc(18,5) = (PA(1) + hermite_term(1,1))*pc(5,5) + ((2.0*PA(1)**2*PB(1) + 2.0*PA&
&(1)*PB(1)**2)*(hermite_term(2,1) - hermite_term(1,1)**2) + (PA(1)**2 + PB(1)**&
&2 + 4.0*PA(1)*PB(1))*(hermite_term(3,1) - hermite_term(2,1)*hermite_term(1,1))&
& + (2.0*PB(1) + 2.0*PA(1))*(hermite_term(4,1) - hermite_term(3,1)*hermite_term&
&(1,1)) + (hermite_term(5,1) - hermite_term(4,1)*hermite_term(1,1)))*pc(1,1)

pc(19,5) = (PA(2) + hermite_term(1,2))*pc(6,5) + ((2.0*PA(2))*(hermite_term(2,&
&2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_te&
&rm(1,2)))*pc(1,5)

pc(20,5) = (PA(3) + hermite_term(1,3))*pc(7,5) + ((2.0*PA(3))*(hermite_term(2,&
&3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_te&
&rm(1,3)))*pc(1,5)

pc(1,6) = (PB(2) + hermite_term(1,2))*pc(1,3) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(1,1)

pc(2,6) = (PA(1) + hermite_term(1,1))*pc(1,6)

pc(3,6) = (PA(2) + hermite_term(1,2))*pc(1,6) + ((2.0*PB(2))*(hermite_term(2,2&
&) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_ter&
&m(1,2)))*pc(1,1)

pc(4,6) = (PA(3) + hermite_term(1,3))*pc(1,6)

pc(5,6) = (PA(1) + hermite_term(1,1))*pc(2,6) + (hermite_term(2,1) - hermite_t&
&erm(1,1)**2)*pc(1,6)

pc(6,6) = (PA(2) + hermite_term(1,2))*pc(3,6) + ((PB(2)**2 + 2.0*PA(2)*PB(2))*&
&(hermite_term(2,2) - hermite_term(1,2)**2) + (2.0*PB(2) + PA(2))*(hermite_term&
&(3,2) - hermite_term(2,2)*hermite_term(1,2))+ (hermite_term(4,2)- hermite_term&
&(3,2)*hermite_term(1,2)))*pc(1,1)

pc(7,6) = (PA(3) + hermite_term(1,3))*pc(4,6) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(1,6)

pc(8,6) = (PA(2) + hermite_term(1,2))*pc(2,6) + ((2.0*PB(2))*(hermite_term(2,2&
&) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_ter&
&m(1,2)))*pc(2,1)

pc(9,6) = (PA(3) + hermite_term(1,3))*pc(2,6)

pc(10,6) = (PA(3) + hermite_term(1,3))*pc(3,6)

pc(11,6) = (PA(3) + hermite_term(1,3))*pc(8,6)

pc(12,6) = (PA(2) + hermite_term(1,2))*pc(5,6) + ((2.0*PB(2))*(hermite_term(2,&
&2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_te&
&rm(1,2)))*pc(5,1)

pc(13,6) = (PA(3) + hermite_term(1,3))*pc(5,6)

pc(14,6) = (PA(1) + hermite_term(1,1))*pc(6,6)

pc(15,6) = (PA(3) + hermite_term(1,3))*pc(6,6)

pc(16,6) = (PA(1) + hermite_term(1,1))*pc(7,6)

pc(17,6) = (PA(2) + hermite_term(1,2))*pc(7,6) + ((2.0*PB(2))*(hermite_term(2,&
&2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_te&
&rm(1,2)))*pc(7,1)

pc(18,6) = (PA(1) + hermite_term(1,1))*pc(5,6) + ((2.0*PA(1))*(hermite_term(2,&
&1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_te&
&rm(1,1)))*pc(1,6)

pc(19,6) = (PA(2) + hermite_term(1,2))*pc(6,6) + ((2.0*PA(2)**2*PB(2) + 2.0*PA&
&(2)*PB(2)**2)*(hermite_term(2,2) - hermite_term(1,2)**2) + (PA(2)**2 + PB(2)**&
&2 + 4.0*PA(2)*PB(2))*(hermite_term(3,2) - hermite_term(2,2)*hermite_term(1,2))&
& + (2.0*PB(2) + 2.0*PA(2))*(hermite_term(4,2) - hermite_term(3,2)*hermite_term&
&(1,2)) + (hermite_term(5,2) - hermite_term(4,2)*hermite_term(1,2)))*pc(1,1)

pc(20,6) = (PA(3) + hermite_term(1,3))*pc(7,6) + ((2.0*PA(3))*(hermite_term(2,&
&3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_te&
&rm(1,3)))*pc(1,6)

pc(1,7) = (PB(3) + hermite_term(1,3))*pc(1,4) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(1,1)

pc(2,7) = (PA(1) + hermite_term(1,1))*pc(1,7)

pc(3,7) = (PA(2) + hermite_term(1,2))*pc(1,7)

pc(4,7) = (PA(3) + hermite_term(1,3))*pc(1,7) + ((2.0*PB(3))*(hermite_term(2,3&
&) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_ter&
&m(1,3)))*pc(1,1)

pc(5,7) = (PA(1) + hermite_term(1,1))*pc(2,7) + (hermite_term(2,1) - hermite_t&
&erm(1,1)**2)*pc(1,7)

pc(6,7) = (PA(2) + hermite_term(1,2))*pc(3,7) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(1,7)

pc(7,7) = (PA(3) + hermite_term(1,3))*pc(4,7) + ((PB(3)**2 + 2.0*PA(3)*PB(3))*&
&(hermite_term(2,3) - hermite_term(1,3)**2) + (2.0*PB(3) + PA(3))*(hermite_term&
&(3,3) - hermite_term(2,3)*hermite_term(1,3))+ (hermite_term(4,3)- hermite_term&
&(3,3)*hermite_term(1,3)))*pc(1,1)

pc(8,7) = (PA(2) + hermite_term(1,2))*pc(2,7)

pc(9,7) = (PA(3) + hermite_term(1,3))*pc(2,7) + ((2.0*PB(3))*(hermite_term(2,3&
&) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_ter&
&m(1,3)))*pc(2,1)

pc(10,7) = (PA(3) + hermite_term(1,3))*pc(3,7) + ((2.0*PB(3))*(hermite_term(2,&
&3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_te&
&rm(1,3)))*pc(3,1)

pc(11,7) = (PA(3) + hermite_term(1,3))*pc(8,7) + ((2.0*PB(3))*(hermite_term(2,&
&3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_te&
&rm(1,3)))*pc(8,1)

pc(12,7) = (PA(2) + hermite_term(1,2))*pc(5,7)

pc(13,7) = (PA(3) + hermite_term(1,3))*pc(5,7) + ((2.0*PB(3))*(hermite_term(2,&
&3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_te&
&rm(1,3)))*pc(5,1)

pc(14,7) = (PA(1) + hermite_term(1,1))*pc(6,7)

pc(15,7) = (PA(3) + hermite_term(1,3))*pc(6,7) + ((2.0*PB(3))*(hermite_term(2,&
&3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_te&
&rm(1,3)))*pc(6,1)

pc(16,7) = (PA(1) + hermite_term(1,1))*pc(7,7)

pc(17,7) = (PA(2) + hermite_term(1,2))*pc(7,7)

pc(18,7) = (PA(1) + hermite_term(1,1))*pc(5,7) + ((2.0*PA(1))*(hermite_term(2,&
&1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_te&
&rm(1,1)))*pc(1,7)

pc(19,7) = (PA(2) + hermite_term(1,2))*pc(6,7) + ((2.0*PA(2))*(hermite_term(2,&
&2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_te&
&rm(1,2)))*pc(1,7)

pc(20,7) = (PA(3) + hermite_term(1,3))*pc(7,7) + ((2.0*PA(3)**2*PB(3) + 2.0*PA&
&(3)*PB(3)**2)*(hermite_term(2,3) - hermite_term(1,3)**2) + (PA(3)**2 + PB(3)**&
&2 + 4.0*PA(3)*PB(3))*(hermite_term(3,3) - hermite_term(2,3)*hermite_term(1,3))&
& + (2.0*PB(3) + 2.0*PA(3))*(hermite_term(4,3) - hermite_term(3,3)*hermite_term&
&(1,3)) + (hermite_term(5,3) - hermite_term(4,3)*hermite_term(1,3)))*pc(1,1)

pc(1,8) = (PB(2) + hermite_term(1,2))*pc(1,2)

pc(2,8) = (PA(1) + hermite_term(1,1))*pc(1,8) + (hermite_term(2,1) - hermite_t&
&erm(1,1)**2)*pc(1,3)

pc(3,8) = (PA(2) + hermite_term(1,2))*pc(1,8) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(1,2)

pc(4,8) = (PA(3) + hermite_term(1,3))*pc(1,8)

pc(5,8) = (PA(1) + hermite_term(1,1))*pc(2,8) + ((PA(1)+PB(1))*(hermite_term(2&
&,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_t&
&erm(1,1)))*pc(1,3)

pc(6,8) = (PA(2) + hermite_term(1,2))*pc(3,8) + ((PA(2)+PB(2))*(hermite_term(2&
&,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_t&
&erm(1,2)))*pc(1,2)

pc(7,8) = (PA(3) + hermite_term(1,3))*pc(4,8) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(1,8)

pc(8,8) = (PA(2) + hermite_term(1,2))*pc(2,8) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(2,2)

pc(9,8) = (PA(3) + hermite_term(1,3))*pc(2,8)

pc(10,8) = (PA(3) + hermite_term(1,3))*pc(3,8)

pc(11,8) = (PA(3) + hermite_term(1,3))*pc(8,8)

pc(12,8) = (PA(2) + hermite_term(1,2))*pc(5,8) + (hermite_term(2,2) - hermite_&
&term(1,2)**2)*pc(5,2)

pc(13,8) = (PA(3) + hermite_term(1,3))*pc(5,8)

pc(14,8) = (PA(1) + hermite_term(1,1))*pc(6,8) + (hermite_term(2,1) - hermite_&
&term(1,1)**2)*pc(6,3)

pc(15,8) = (PA(3) + hermite_term(1,3))*pc(6,8)

pc(16,8) = (PA(1) + hermite_term(1,1))*pc(7,8) + (hermite_term(2,1) - hermite_&
&term(1,1)**2)*pc(7,3)

pc(17,8) = (PA(2) + hermite_term(1,2))*pc(7,8) + (hermite_term(2,2) - hermite_&
&term(1,2)**2)*pc(7,2)

pc(18,8) = (PA(1) + hermite_term(1,1))*pc(5,8) + ((PA(1)**2 + 2.0*PA(1)*PB(1))&
&*(hermite_term(2,1) - hermite_term(1,1)**2) + (2.0*PA(1) + PB(1))*(hermite_ter&
&m(3,1) - hermite_term(2,1)*hermite_term(1,1)) + (hermite_term(4,1) - hermite_t&
&erm(3,1)*hermite_term(1,1)))*pc(1,3)

pc(19,8) = (PA(2) + hermite_term(1,2))*pc(6,8) + ((PA(2)**2 + 2.0*PA(2)*PB(2))&
&*(hermite_term(2,2) - hermite_term(1,2)**2) + (2.0*PA(2) + PB(2))*(hermite_ter&
&m(3,2) - hermite_term(2,2)*hermite_term(1,2)) + (hermite_term(4,2) - hermite_t&
&erm(3,2)*hermite_term(1,2)))*pc(1,2)

pc(20,8) = (PA(3) + hermite_term(1,3))*pc(7,8) + ((2.0*PA(3))*(hermite_term(2,&
&3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_te&
&rm(1,3)))*pc(1,8)

pc(1,9) = (PB(3) + hermite_term(1,3))*pc(1,2)

pc(2,9) = (PA(1) + hermite_term(1,1))*pc(1,9) + (hermite_term(2,1) - hermite_t&
&erm(1,1)**2)*pc(1,4)

pc(3,9) = (PA(2) + hermite_term(1,2))*pc(1,9)

pc(4,9) = (PA(3) + hermite_term(1,3))*pc(1,9) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(1,2)

pc(5,9) = (PA(1) + hermite_term(1,1))*pc(2,9) + ((PA(1)+PB(1))*(hermite_term(2&
&,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_t&
&erm(1,1)))*pc(1,4)

pc(6,9) = (PA(2) + hermite_term(1,2))*pc(3,9) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(1,9)

pc(7,9) = (PA(3) + hermite_term(1,3))*pc(4,9) + ((PA(3)+PB(3))*(hermite_term(2&
&,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_t&
&erm(1,3)))*pc(1,2)

pc(8,9) = (PA(2) + hermite_term(1,2))*pc(2,9)

pc(9,9) = (PA(3) + hermite_term(1,3))*pc(2,9) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(2,2)

pc(10,9) = (PA(3) + hermite_term(1,3))*pc(3,9) + (hermite_term(2,3) - hermite_&
&term(1,3)**2)*pc(3,2)

pc(11,9) = (PA(3) + hermite_term(1,3))*pc(8,9) + (hermite_term(2,3) - hermite_&
&term(1,3)**2)*pc(8,2)

pc(12,9) = (PA(2) + hermite_term(1,2))*pc(5,9)

pc(13,9) = (PA(3) + hermite_term(1,3))*pc(5,9) + (hermite_term(2,3) - hermite_&
&term(1,3)**2)*pc(5,2)

pc(14,9) = (PA(1) + hermite_term(1,1))*pc(6,9) + (hermite_term(2,1) - hermite_&
&term(1,1)**2)*pc(6,4)

pc(15,9) = (PA(3) + hermite_term(1,3))*pc(6,9) + (hermite_term(2,3) - hermite_&
&term(1,3)**2)*pc(6,2)

pc(16,9) = (PA(1) + hermite_term(1,1))*pc(7,9) + (hermite_term(2,1) - hermite_&
&term(1,1)**2)*pc(7,4)

pc(17,9) = (PA(2) + hermite_term(1,2))*pc(7,9)

pc(18,9) = (PA(1) + hermite_term(1,1))*pc(5,9) + ((PA(1)**2 + 2.0*PA(1)*PB(1))&
&*(hermite_term(2,1) - hermite_term(1,1)**2) + (2.0*PA(1) + PB(1))*(hermite_ter&
&m(3,1) - hermite_term(2,1)*hermite_term(1,1)) + (hermite_term(4,1) - hermite_t&
&erm(3,1)*hermite_term(1,1)))*pc(1,4)

pc(19,9) = (PA(2) + hermite_term(1,2))*pc(6,9) + ((2.0*PA(2))*(hermite_term(2,&
&2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_te&
&rm(1,2)))*pc(1,9)

pc(20,9) = (PA(3) + hermite_term(1,3))*pc(7,9) + ((PA(3)**2 + 2.0*PA(3)*PB(3))&
&*(hermite_term(2,3) - hermite_term(1,3)**2) + (2.0*PA(3) + PB(3))*(hermite_ter&
&m(3,3) - hermite_term(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite_t&
&erm(3,3)*hermite_term(1,3)))*pc(1,2)

pc(1,10) = (PB(3) + hermite_term(1,3))*pc(1,3)

pc(2,10) = (PA(1) + hermite_term(1,1))*pc(1,10)

pc(3,10) = (PA(2) + hermite_term(1,2))*pc(1,10) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(1,4)

pc(4,10) = (PA(3) + hermite_term(1,3))*pc(1,10) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(1,3)

pc(5,10) = (PA(1) + hermite_term(1,1))*pc(2,10) + (hermite_term(2,1) - hermite&
&_term(1,1)**2)*pc(1,10)

pc(6,10) = (PA(2) + hermite_term(1,2))*pc(3,10) + ((PA(2)+PB(2))*(hermite_term&
&(2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite&
&_term(1,2)))*pc(1,4)

pc(7,10) = (PA(3) + hermite_term(1,3))*pc(4,10) + ((PA(3)+PB(3))*(hermite_term&
&(2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite&
&_term(1,3)))*pc(1,3)

pc(8,10) = (PA(2) + hermite_term(1,2))*pc(2,10) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(2,4)

pc(9,10) = (PA(3) + hermite_term(1,3))*pc(2,10) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(2,3)

pc(10,10) = (PA(3) + hermite_term(1,3))*pc(3,10) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(3,3)

pc(11,10) = (PA(3) + hermite_term(1,3))*pc(8,10) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(8,3)

pc(12,10) = (PA(2) + hermite_term(1,2))*pc(5,10) + (hermite_term(2,2) - hermit&
&e_term(1,2)**2)*pc(5,4)

pc(13,10) = (PA(3) + hermite_term(1,3))*pc(5,10) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(5,3)

pc(14,10) = (PA(1) + hermite_term(1,1))*pc(6,10)

pc(15,10) = (PA(3) + hermite_term(1,3))*pc(6,10) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(6,3)

pc(16,10) = (PA(1) + hermite_term(1,1))*pc(7,10)

pc(17,10) = (PA(2) + hermite_term(1,2))*pc(7,10) + (hermite_term(2,2) - hermit&
&e_term(1,2)**2)*pc(7,4)

pc(18,10) = (PA(1) + hermite_term(1,1))*pc(5,10) + ((2.0*PA(1))*(hermite_term(&
&2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_&
&term(1,1)))*pc(1,10)

pc(19,10) = (PA(2) + hermite_term(1,2))*pc(6,10) + ((PA(2)**2 + 2.0*PA(2)*PB(2&
&))*(hermite_term(2,2) - hermite_term(1,2)**2) + (2.0*PA(2) + PB(2))*(hermite_t&
&erm(3,2) - hermite_term(2,2)*hermite_term(1,2)) + (hermite_term(4,2) - hermite&
&_term(3,2)*hermite_term(1,2)))*pc(1,4)

pc(20,10) = (PA(3) + hermite_term(1,3))*pc(7,10) + ((PA(3)**2 + 2.0*PA(3)*PB(3&
&))*(hermite_term(2,3) - hermite_term(1,3)**2) + (2.0*PA(3) + PB(3))*(hermite_t&
&erm(3,3) - hermite_term(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite&
&_term(3,3)*hermite_term(1,3)))*pc(1,3)

pc(1,11) = (PB(3) + hermite_term(1,3))*pc(1,8)

pc(2,11) = (PA(1) + hermite_term(1,1))*pc(1,11) + (hermite_term(2,1) - hermite&
&_term(1,1)**2)*pc(1,10)

pc(3,11) = (PA(2) + hermite_term(1,2))*pc(1,11) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(1,9)

pc(4,11) = (PA(3) + hermite_term(1,3))*pc(1,11) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(1,8)

pc(5,11) = (PA(1) + hermite_term(1,1))*pc(2,11) + ((PA(1)+PB(1))*(hermite_term&
&(2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite&
&_term(1,1)))*pc(1,10)

pc(6,11) = (PA(2) + hermite_term(1,2))*pc(3,11) + ((PA(2)+PB(2))*(hermite_term&
&(2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite&
&_term(1,2)))*pc(1,9)

pc(7,11) = (PA(3) + hermite_term(1,3))*pc(4,11) + ((PA(3)+PB(3))*(hermite_term&
&(2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite&
&_term(1,3)))*pc(1,8)

pc(8,11) = (PA(2) + hermite_term(1,2))*pc(2,11) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(2,9)

pc(9,11) = (PA(3) + hermite_term(1,3))*pc(2,11) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(2,8)

pc(10,11) = (PA(3) + hermite_term(1,3))*pc(3,11) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(3,8)

pc(11,11) = (PA(3) + hermite_term(1,3))*pc(8,11) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(8,8)

pc(12,11) = (PA(2) + hermite_term(1,2))*pc(5,11) + (hermite_term(2,2) - hermit&
&e_term(1,2)**2)*pc(5,9)

pc(13,11) = (PA(3) + hermite_term(1,3))*pc(5,11) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(5,8)

pc(14,11) = (PA(1) + hermite_term(1,1))*pc(6,11) + (hermite_term(2,1) - hermit&
&e_term(1,1)**2)*pc(6,10)

pc(15,11) = (PA(3) + hermite_term(1,3))*pc(6,11) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(6,8)

pc(16,11) = (PA(1) + hermite_term(1,1))*pc(7,11) + (hermite_term(2,1) - hermit&
&e_term(1,1)**2)*pc(7,10)

pc(17,11) = (PA(2) + hermite_term(1,2))*pc(7,11) + (hermite_term(2,2) - hermit&
&e_term(1,2)**2)*pc(7,9)

pc(18,11) = (PA(1) + hermite_term(1,1))*pc(5,11) + ((PA(1)**2 + 2.0*PA(1)*PB(1&
&))*(hermite_term(2,1) - hermite_term(1,1)**2) + (2.0*PA(1) + PB(1))*(hermite_t&
&erm(3,1) - hermite_term(2,1)*hermite_term(1,1)) + (hermite_term(4,1) - hermite&
&_term(3,1)*hermite_term(1,1)))*pc(1,10)

pc(19,11) = (PA(2) + hermite_term(1,2))*pc(6,11) + ((PA(2)**2 + 2.0*PA(2)*PB(2&
&))*(hermite_term(2,2) - hermite_term(1,2)**2) + (2.0*PA(2) + PB(2))*(hermite_t&
&erm(3,2) - hermite_term(2,2)*hermite_term(1,2)) + (hermite_term(4,2) - hermite&
&_term(3,2)*hermite_term(1,2)))*pc(1,9)

pc(20,11) = (PA(3) + hermite_term(1,3))*pc(7,11) + ((PA(3)**2 + 2.0*PA(3)*PB(3&
&))*(hermite_term(2,3) - hermite_term(1,3)**2) + (2.0*PA(3) + PB(3))*(hermite_t&
&erm(3,3) - hermite_term(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite&
&_term(3,3)*hermite_term(1,3)))*pc(1,8)

pc(1,12) = (PB(2) + hermite_term(1,2))*pc(1,5)

pc(2,12) = (PA(1) + hermite_term(1,1))*pc(1,12) + ((2.0*PB(1))*(hermite_term(2&
&,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_t&
&erm(1,1)))*pc(1,3)

pc(3,12) = (PA(2) + hermite_term(1,2))*pc(1,12) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(1,5)

pc(4,12) = (PA(3) + hermite_term(1,3))*pc(1,12)

pc(5,12) = (PA(1) + hermite_term(1,1))*pc(2,12) + ((PB(1)**2 + 2.0*PA(1)*PB(1)&
&)*(hermite_term(2,1) - hermite_term(1,1)**2) + (2.0*PB(1) + PA(1))*(hermite_te&
&rm(3,1) - hermite_term(2,1)*hermite_term(1,1))+ (hermite_term(4,1)- hermite_te&
&rm(3,1)*hermite_term(1,1)))*pc(1,3)

pc(6,12) = (PA(2) + hermite_term(1,2))*pc(3,12) + ((PA(2)+PB(2))*(hermite_term&
&(2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite&
&_term(1,2)))*pc(1,5)

pc(7,12) = (PA(3) + hermite_term(1,3))*pc(4,12) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(1,12)

pc(8,12) = (PA(2) + hermite_term(1,2))*pc(2,12) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(2,5)

pc(9,12) = (PA(3) + hermite_term(1,3))*pc(2,12)

pc(10,12) = (PA(3) + hermite_term(1,3))*pc(3,12)

pc(11,12) = (PA(3) + hermite_term(1,3))*pc(8,12)

pc(12,12) = (PA(2) + hermite_term(1,2))*pc(5,12) + (hermite_term(2,2) - hermit&
&e_term(1,2)**2)*pc(5,5)

pc(13,12) = (PA(3) + hermite_term(1,3))*pc(5,12)

pc(14,12) = (PA(1) + hermite_term(1,1))*pc(6,12) + ((2.0*PB(1))*(hermite_term(&
&2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_&
&term(1,1)))*pc(6,3)

pc(15,12) = (PA(3) + hermite_term(1,3))*pc(6,12)

pc(16,12) = (PA(1) + hermite_term(1,1))*pc(7,12) + ((2.0*PB(1))*(hermite_term(&
&2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_&
&term(1,1)))*pc(7,3)

pc(17,12) = (PA(2) + hermite_term(1,2))*pc(7,12) + (hermite_term(2,2) - hermit&
&e_term(1,2)**2)*pc(7,5)

pc(18,12) = (PA(1) + hermite_term(1,1))*pc(5,12) + ((2.0*PA(1)**2*PB(1) + 2.0*&
&PA(1)*PB(1)**2)*(hermite_term(2,1) - hermite_term(1,1)**2) + (PA(1)**2 + PB(1)&
&**2 + 4.0*PA(1)*PB(1))*(hermite_term(3,1) - hermite_term(2,1)*hermite_term(1,1&
&)) + (2.0*PB(1) + 2.0*PA(1))*(hermite_term(4,1) - hermite_term(3,1)*hermite_te&
&rm(1,1)) + (hermite_term(5,1) - hermite_term(4,1)*hermite_term(1,1)))*pc(1,3)

pc(19,12) = (PA(2) + hermite_term(1,2))*pc(6,12) + ((PA(2)**2 + 2.0*PA(2)*PB(2&
&))*(hermite_term(2,2) - hermite_term(1,2)**2) + (2.0*PA(2) + PB(2))*(hermite_t&
&erm(3,2) - hermite_term(2,2)*hermite_term(1,2)) + (hermite_term(4,2) - hermite&
&_term(3,2)*hermite_term(1,2)))*pc(1,5)

pc(20,12) = (PA(3) + hermite_term(1,3))*pc(7,12) + ((2.0*PA(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(1,12)

pc(1,13) = (PB(3) + hermite_term(1,3))*pc(1,5)

pc(2,13) = (PA(1) + hermite_term(1,1))*pc(1,13) + ((2.0*PB(1))*(hermite_term(2&
&,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_t&
&erm(1,1)))*pc(1,4)

pc(3,13) = (PA(2) + hermite_term(1,2))*pc(1,13)

pc(4,13) = (PA(3) + hermite_term(1,3))*pc(1,13) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(1,5)

pc(5,13) = (PA(1) + hermite_term(1,1))*pc(2,13) + ((PB(1)**2 + 2.0*PA(1)*PB(1)&
&)*(hermite_term(2,1) - hermite_term(1,1)**2) + (2.0*PB(1) + PA(1))*(hermite_te&
&rm(3,1) - hermite_term(2,1)*hermite_term(1,1))+ (hermite_term(4,1)- hermite_te&
&rm(3,1)*hermite_term(1,1)))*pc(1,4)

pc(6,13) = (PA(2) + hermite_term(1,2))*pc(3,13) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(1,13)

pc(7,13) = (PA(3) + hermite_term(1,3))*pc(4,13) + ((PA(3)+PB(3))*(hermite_term&
&(2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite&
&_term(1,3)))*pc(1,5)

pc(8,13) = (PA(2) + hermite_term(1,2))*pc(2,13)

pc(9,13) = (PA(3) + hermite_term(1,3))*pc(2,13) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(2,5)

pc(10,13) = (PA(3) + hermite_term(1,3))*pc(3,13) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(3,5)

pc(11,13) = (PA(3) + hermite_term(1,3))*pc(8,13) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(8,5)

pc(12,13) = (PA(2) + hermite_term(1,2))*pc(5,13)

pc(13,13) = (PA(3) + hermite_term(1,3))*pc(5,13) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(5,5)

pc(14,13) = (PA(1) + hermite_term(1,1))*pc(6,13) + ((2.0*PB(1))*(hermite_term(&
&2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_&
&term(1,1)))*pc(6,4)

pc(15,13) = (PA(3) + hermite_term(1,3))*pc(6,13) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(6,5)

pc(16,13) = (PA(1) + hermite_term(1,1))*pc(7,13) + ((2.0*PB(1))*(hermite_term(&
&2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_&
&term(1,1)))*pc(7,4)

pc(17,13) = (PA(2) + hermite_term(1,2))*pc(7,13)

pc(18,13) = (PA(1) + hermite_term(1,1))*pc(5,13) + ((2.0*PA(1)**2*PB(1) + 2.0*&
&PA(1)*PB(1)**2)*(hermite_term(2,1) - hermite_term(1,1)**2) + (PA(1)**2 + PB(1)&
&**2 + 4.0*PA(1)*PB(1))*(hermite_term(3,1) - hermite_term(2,1)*hermite_term(1,1&
&)) + (2.0*PB(1) + 2.0*PA(1))*(hermite_term(4,1) - hermite_term(3,1)*hermite_te&
&rm(1,1)) + (hermite_term(5,1) - hermite_term(4,1)*hermite_term(1,1)))*pc(1,4)

pc(19,13) = (PA(2) + hermite_term(1,2))*pc(6,13) + ((2.0*PA(2))*(hermite_term(&
&2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_&
&term(1,2)))*pc(1,13)

pc(20,13) = (PA(3) + hermite_term(1,3))*pc(7,13) + ((PA(3)**2 + 2.0*PA(3)*PB(3&
&))*(hermite_term(2,3) - hermite_term(1,3)**2) + (2.0*PA(3) + PB(3))*(hermite_t&
&erm(3,3) - hermite_term(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite&
&_term(3,3)*hermite_term(1,3)))*pc(1,5)

pc(1,14) = (PB(1) + hermite_term(1,1))*pc(1,6)

pc(2,14) = (PA(1) + hermite_term(1,1))*pc(1,14) + (hermite_term(2,1) - hermite&
&_term(1,1)**2)*pc(1,6)

pc(3,14) = (PA(2) + hermite_term(1,2))*pc(1,14) + ((2.0*PB(2))*(hermite_term(2&
&,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_t&
&erm(1,2)))*pc(1,2)

pc(4,14) = (PA(3) + hermite_term(1,3))*pc(1,14)

pc(5,14) = (PA(1) + hermite_term(1,1))*pc(2,14) + ((PA(1)+PB(1))*(hermite_term&
&(2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite&
&_term(1,1)))*pc(1,6)

pc(6,14) = (PA(2) + hermite_term(1,2))*pc(3,14) + ((PB(2)**2 + 2.0*PA(2)*PB(2)&
&)*(hermite_term(2,2) - hermite_term(1,2)**2) + (2.0*PB(2) + PA(2))*(hermite_te&
&rm(3,2) - hermite_term(2,2)*hermite_term(1,2))+ (hermite_term(4,2)- hermite_te&
&rm(3,2)*hermite_term(1,2)))*pc(1,2)

pc(7,14) = (PA(3) + hermite_term(1,3))*pc(4,14) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(1,14)

pc(8,14) = (PA(2) + hermite_term(1,2))*pc(2,14) + ((2.0*PB(2))*(hermite_term(2&
&,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_t&
&erm(1,2)))*pc(2,2)

pc(9,14) = (PA(3) + hermite_term(1,3))*pc(2,14)

pc(10,14) = (PA(3) + hermite_term(1,3))*pc(3,14)

pc(11,14) = (PA(3) + hermite_term(1,3))*pc(8,14)

pc(12,14) = (PA(2) + hermite_term(1,2))*pc(5,14) + ((2.0*PB(2))*(hermite_term(&
&2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_&
&term(1,2)))*pc(5,2)

pc(13,14) = (PA(3) + hermite_term(1,3))*pc(5,14)

pc(14,14) = (PA(1) + hermite_term(1,1))*pc(6,14) + (hermite_term(2,1) - hermit&
&e_term(1,1)**2)*pc(6,6)

pc(15,14) = (PA(3) + hermite_term(1,3))*pc(6,14)

pc(16,14) = (PA(1) + hermite_term(1,1))*pc(7,14) + (hermite_term(2,1) - hermit&
&e_term(1,1)**2)*pc(7,6)

pc(17,14) = (PA(2) + hermite_term(1,2))*pc(7,14) + ((2.0*PB(2))*(hermite_term(&
&2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_&
&term(1,2)))*pc(7,2)

pc(18,14) = (PA(1) + hermite_term(1,1))*pc(5,14) + ((PA(1)**2 + 2.0*PA(1)*PB(1&
&))*(hermite_term(2,1) - hermite_term(1,1)**2) + (2.0*PA(1) + PB(1))*(hermite_t&
&erm(3,1) - hermite_term(2,1)*hermite_term(1,1)) + (hermite_term(4,1) - hermite&
&_term(3,1)*hermite_term(1,1)))*pc(1,6)

pc(19,14) = (PA(2) + hermite_term(1,2))*pc(6,14) + ((2.0*PA(2)**2*PB(2) + 2.0*&
&PA(2)*PB(2)**2)*(hermite_term(2,2) - hermite_term(1,2)**2) + (PA(2)**2 + PB(2)&
&**2 + 4.0*PA(2)*PB(2))*(hermite_term(3,2) - hermite_term(2,2)*hermite_term(1,2&
&)) + (2.0*PB(2) + 2.0*PA(2))*(hermite_term(4,2) - hermite_term(3,2)*hermite_te&
&rm(1,2)) + (hermite_term(5,2) - hermite_term(4,2)*hermite_term(1,2)))*pc(1,2)

pc(20,14) = (PA(3) + hermite_term(1,3))*pc(7,14) + ((2.0*PA(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(1,14)

pc(1,15) = (PB(3) + hermite_term(1,3))*pc(1,6)

pc(2,15) = (PA(1) + hermite_term(1,1))*pc(1,15)

pc(3,15) = (PA(2) + hermite_term(1,2))*pc(1,15) + ((2.0*PB(2))*(hermite_term(2&
&,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_t&
&erm(1,2)))*pc(1,4)

pc(4,15) = (PA(3) + hermite_term(1,3))*pc(1,15) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(1,6)

pc(5,15) = (PA(1) + hermite_term(1,1))*pc(2,15) + (hermite_term(2,1) - hermite&
&_term(1,1)**2)*pc(1,15)

pc(6,15) = (PA(2) + hermite_term(1,2))*pc(3,15) + ((PB(2)**2 + 2.0*PA(2)*PB(2)&
&)*(hermite_term(2,2) - hermite_term(1,2)**2) + (2.0*PB(2) + PA(2))*(hermite_te&
&rm(3,2) - hermite_term(2,2)*hermite_term(1,2))+ (hermite_term(4,2)- hermite_te&
&rm(3,2)*hermite_term(1,2)))*pc(1,4)

pc(7,15) = (PA(3) + hermite_term(1,3))*pc(4,15) + ((PA(3)+PB(3))*(hermite_term&
&(2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite&
&_term(1,3)))*pc(1,6)

pc(8,15) = (PA(2) + hermite_term(1,2))*pc(2,15) + ((2.0*PB(2))*(hermite_term(2&
&,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_t&
&erm(1,2)))*pc(2,4)

pc(9,15) = (PA(3) + hermite_term(1,3))*pc(2,15) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(2,6)

pc(10,15) = (PA(3) + hermite_term(1,3))*pc(3,15) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(3,6)

pc(11,15) = (PA(3) + hermite_term(1,3))*pc(8,15) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(8,6)

pc(12,15) = (PA(2) + hermite_term(1,2))*pc(5,15) + ((2.0*PB(2))*(hermite_term(&
&2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_&
&term(1,2)))*pc(5,4)

pc(13,15) = (PA(3) + hermite_term(1,3))*pc(5,15) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(5,6)

pc(14,15) = (PA(1) + hermite_term(1,1))*pc(6,15)

pc(15,15) = (PA(3) + hermite_term(1,3))*pc(6,15) + (hermite_term(2,3) - hermit&
&e_term(1,3)**2)*pc(6,6)

pc(16,15) = (PA(1) + hermite_term(1,1))*pc(7,15)

pc(17,15) = (PA(2) + hermite_term(1,2))*pc(7,15) + ((2.0*PB(2))*(hermite_term(&
&2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_&
&term(1,2)))*pc(7,4)

pc(18,15) = (PA(1) + hermite_term(1,1))*pc(5,15) + ((2.0*PA(1))*(hermite_term(&
&2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_&
&term(1,1)))*pc(1,15)

pc(19,15) = (PA(2) + hermite_term(1,2))*pc(6,15) + ((2.0*PA(2)**2*PB(2) + 2.0*&
&PA(2)*PB(2)**2)*(hermite_term(2,2) - hermite_term(1,2)**2) + (PA(2)**2 + PB(2)&
&**2 + 4.0*PA(2)*PB(2))*(hermite_term(3,2) - hermite_term(2,2)*hermite_term(1,2&
&)) + (2.0*PB(2) + 2.0*PA(2))*(hermite_term(4,2) - hermite_term(3,2)*hermite_te&
&rm(1,2)) + (hermite_term(5,2) - hermite_term(4,2)*hermite_term(1,2)))*pc(1,4)

pc(20,15) = (PA(3) + hermite_term(1,3))*pc(7,15) + ((PA(3)**2 + 2.0*PA(3)*PB(3&
&))*(hermite_term(2,3) - hermite_term(1,3)**2) + (2.0*PA(3) + PB(3))*(hermite_t&
&erm(3,3) - hermite_term(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite&
&_term(3,3)*hermite_term(1,3)))*pc(1,6)

pc(1,16) = (PB(1) + hermite_term(1,1))*pc(1,7)

pc(2,16) = (PA(1) + hermite_term(1,1))*pc(1,16) + (hermite_term(2,1) - hermite&
&_term(1,1)**2)*pc(1,7)

pc(3,16) = (PA(2) + hermite_term(1,2))*pc(1,16)

pc(4,16) = (PA(3) + hermite_term(1,3))*pc(1,16) + ((2.0*PB(3))*(hermite_term(2&
&,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_t&
&erm(1,3)))*pc(1,2)

pc(5,16) = (PA(1) + hermite_term(1,1))*pc(2,16) + ((PA(1)+PB(1))*(hermite_term&
&(2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite&
&_term(1,1)))*pc(1,7)

pc(6,16) = (PA(2) + hermite_term(1,2))*pc(3,16) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(1,16)

pc(7,16) = (PA(3) + hermite_term(1,3))*pc(4,16) + ((PB(3)**2 + 2.0*PA(3)*PB(3)&
&)*(hermite_term(2,3) - hermite_term(1,3)**2) + (2.0*PB(3) + PA(3))*(hermite_te&
&rm(3,3) - hermite_term(2,3)*hermite_term(1,3))+ (hermite_term(4,3)- hermite_te&
&rm(3,3)*hermite_term(1,3)))*pc(1,2)

pc(8,16) = (PA(2) + hermite_term(1,2))*pc(2,16)

pc(9,16) = (PA(3) + hermite_term(1,3))*pc(2,16) + ((2.0*PB(3))*(hermite_term(2&
&,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_t&
&erm(1,3)))*pc(2,2)

pc(10,16) = (PA(3) + hermite_term(1,3))*pc(3,16) + ((2.0*PB(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(3,2)

pc(11,16) = (PA(3) + hermite_term(1,3))*pc(8,16) + ((2.0*PB(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(8,2)

pc(12,16) = (PA(2) + hermite_term(1,2))*pc(5,16)

pc(13,16) = (PA(3) + hermite_term(1,3))*pc(5,16) + ((2.0*PB(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(5,2)

pc(14,16) = (PA(1) + hermite_term(1,1))*pc(6,16) + (hermite_term(2,1) - hermit&
&e_term(1,1)**2)*pc(6,7)

pc(15,16) = (PA(3) + hermite_term(1,3))*pc(6,16) + ((2.0*PB(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(6,2)

pc(16,16) = (PA(1) + hermite_term(1,1))*pc(7,16) + (hermite_term(2,1) - hermit&
&e_term(1,1)**2)*pc(7,7)

pc(17,16) = (PA(2) + hermite_term(1,2))*pc(7,16)

pc(18,16) = (PA(1) + hermite_term(1,1))*pc(5,16) + ((PA(1)**2 + 2.0*PA(1)*PB(1&
&))*(hermite_term(2,1) - hermite_term(1,1)**2) + (2.0*PA(1) + PB(1))*(hermite_t&
&erm(3,1) - hermite_term(2,1)*hermite_term(1,1)) + (hermite_term(4,1) - hermite&
&_term(3,1)*hermite_term(1,1)))*pc(1,7)

pc(19,16) = (PA(2) + hermite_term(1,2))*pc(6,16) + ((2.0*PA(2))*(hermite_term(&
&2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_&
&term(1,2)))*pc(1,16)

pc(20,16) = (PA(3) + hermite_term(1,3))*pc(7,16) + ((2.0*PA(3)**2*PB(3) + 2.0*&
&PA(3)*PB(3)**2)*(hermite_term(2,3) - hermite_term(1,3)**2) + (PA(3)**2 + PB(3)&
&**2 + 4.0*PA(3)*PB(3))*(hermite_term(3,3) - hermite_term(2,3)*hermite_term(1,3&
&)) + (2.0*PB(3) + 2.0*PA(3))*(hermite_term(4,3) - hermite_term(3,3)*hermite_te&
&rm(1,3)) + (hermite_term(5,3) - hermite_term(4,3)*hermite_term(1,3)))*pc(1,2)

pc(1,17) = (PB(2) + hermite_term(1,2))*pc(1,7)

pc(2,17) = (PA(1) + hermite_term(1,1))*pc(1,17)

pc(3,17) = (PA(2) + hermite_term(1,2))*pc(1,17) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(1,7)

pc(4,17) = (PA(3) + hermite_term(1,3))*pc(1,17) + ((2.0*PB(3))*(hermite_term(2&
&,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_t&
&erm(1,3)))*pc(1,3)

pc(5,17) = (PA(1) + hermite_term(1,1))*pc(2,17) + (hermite_term(2,1) - hermite&
&_term(1,1)**2)*pc(1,17)

pc(6,17) = (PA(2) + hermite_term(1,2))*pc(3,17) + ((PA(2)+PB(2))*(hermite_term&
&(2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite&
&_term(1,2)))*pc(1,7)

pc(7,17) = (PA(3) + hermite_term(1,3))*pc(4,17) + ((PB(3)**2 + 2.0*PA(3)*PB(3)&
&)*(hermite_term(2,3) - hermite_term(1,3)**2) + (2.0*PB(3) + PA(3))*(hermite_te&
&rm(3,3) - hermite_term(2,3)*hermite_term(1,3))+ (hermite_term(4,3)- hermite_te&
&rm(3,3)*hermite_term(1,3)))*pc(1,3)

pc(8,17) = (PA(2) + hermite_term(1,2))*pc(2,17) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(2,7)

pc(9,17) = (PA(3) + hermite_term(1,3))*pc(2,17) + ((2.0*PB(3))*(hermite_term(2&
&,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_t&
&erm(1,3)))*pc(2,3)

pc(10,17) = (PA(3) + hermite_term(1,3))*pc(3,17) + ((2.0*PB(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(3,3)

pc(11,17) = (PA(3) + hermite_term(1,3))*pc(8,17) + ((2.0*PB(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(8,3)

pc(12,17) = (PA(2) + hermite_term(1,2))*pc(5,17) + (hermite_term(2,2) - hermit&
&e_term(1,2)**2)*pc(5,7)

pc(13,17) = (PA(3) + hermite_term(1,3))*pc(5,17) + ((2.0*PB(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(5,3)

pc(14,17) = (PA(1) + hermite_term(1,1))*pc(6,17)

pc(15,17) = (PA(3) + hermite_term(1,3))*pc(6,17) + ((2.0*PB(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(6,3)

pc(16,17) = (PA(1) + hermite_term(1,1))*pc(7,17)

pc(17,17) = (PA(2) + hermite_term(1,2))*pc(7,17) + (hermite_term(2,2) - hermit&
&e_term(1,2)**2)*pc(7,7)

pc(18,17) = (PA(1) + hermite_term(1,1))*pc(5,17) + ((2.0*PA(1))*(hermite_term(&
&2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_&
&term(1,1)))*pc(1,17)

pc(19,17) = (PA(2) + hermite_term(1,2))*pc(6,17) + ((PA(2)**2 + 2.0*PA(2)*PB(2&
&))*(hermite_term(2,2) - hermite_term(1,2)**2) + (2.0*PA(2) + PB(2))*(hermite_t&
&erm(3,2) - hermite_term(2,2)*hermite_term(1,2)) + (hermite_term(4,2) - hermite&
&_term(3,2)*hermite_term(1,2)))*pc(1,7)

pc(20,17) = (PA(3) + hermite_term(1,3))*pc(7,17) + ((2.0*PA(3)**2*PB(3) + 2.0*&
&PA(3)*PB(3)**2)*(hermite_term(2,3) - hermite_term(1,3)**2) + (PA(3)**2 + PB(3)&
&**2 + 4.0*PA(3)*PB(3))*(hermite_term(3,3) - hermite_term(2,3)*hermite_term(1,3&
&)) + (2.0*PB(3) + 2.0*PA(3))*(hermite_term(4,3) - hermite_term(3,3)*hermite_te&
&rm(1,3)) + (hermite_term(5,3) - hermite_term(4,3)*hermite_term(1,3)))*pc(1,3)

pc(1,18) = (PB(1) + hermite_term(1,1))*pc(1,5) + ((2.0*PB(1))*(hermite_term(2,&
&1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_te&
&rm(1,1)))*pc(1,1)

pc(2,18) = (PA(1) + hermite_term(1,1))*pc(1,18) + ((3.0*(PB(1)**2))*(hermite_t&
&erm(2,1) - hermite_term(1,1)**2) + (3*PB(1))*(hermite_term(3,1) - hermite_term&
&(2,1)*hermite_term(1,1)) + (hermite_term(4,1) - hermite_term(3,1)*hermite_term&
&(1,1)))*pc(1,1)

pc(3,18) = (PA(2) + hermite_term(1,2))*pc(1,18)

pc(4,18) = (PA(3) + hermite_term(1,3))*pc(1,18)

pc(5,18) = (PA(1) + hermite_term(1,1))*pc(2,18) + ((3.0*(PB(1)**2)*PA(1) + PB(&
&1)**3)*(hermite_term(2,1) - hermite_term(1,1)**2) + (3.0*PB(1)*PA(1) + 3.0*(PB&
&(1)**2))*(hermite_term(3,1) - hermite_term(2,1)* hermite_term(1,1)) + (PA(1) +&
& 3.0*PB(1))*(hermite_term(4,1) - hermite_term(3,1)*hermite_term(1,1))+ (hermit&
&e_term(5,1)- hermite_term(4,1)*hermite_term(1,1)))*pc(1,1)

pc(6,18) = (PA(2) + hermite_term(1,2))*pc(3,18) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(1,18)

pc(7,18) = (PA(3) + hermite_term(1,3))*pc(4,18) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(1,18)

pc(8,18) = (PA(2) + hermite_term(1,2))*pc(2,18)

pc(9,18) = (PA(3) + hermite_term(1,3))*pc(2,18)

pc(10,18) = (PA(3) + hermite_term(1,3))*pc(3,18)

pc(11,18) = (PA(3) + hermite_term(1,3))*pc(8,18)

pc(12,18) = (PA(2) + hermite_term(1,2))*pc(5,18)

pc(13,18) = (PA(3) + hermite_term(1,3))*pc(5,18)

pc(14,18) = (PA(1) + hermite_term(1,1))*pc(6,18) + ((3.0*(PB(1)**2))*(hermite_&
&term(2,1) - hermite_term(1,1)**2) + (3*PB(1))*(hermite_term(3,1) - hermite_ter&
&m(2,1)*hermite_term(1,1)) + (hermite_term(4,1) - hermite_term(3,1)*hermite_ter&
&m(1,1)))*pc(6,1)

pc(15,18) = (PA(3) + hermite_term(1,3))*pc(6,18)

pc(16,18) = (PA(1) + hermite_term(1,1))*pc(7,18) + ((3.0*(PB(1)**2))*(hermite_&
&term(2,1) - hermite_term(1,1)**2) + (3*PB(1))*(hermite_term(3,1) - hermite_ter&
&m(2,1)*hermite_term(1,1)) + (hermite_term(4,1) - hermite_term(3,1)*hermite_ter&
&m(1,1)))*pc(7,1)

pc(17,18) = (PA(2) + hermite_term(1,2))*pc(7,18)

pc(18,18) = (PA(1) + hermite_term(1,1))*pc(5,18) + ((2.0*(PB(1)**3)*PA(1) + 3.&
&0*(PA(1)**2)*(PB(1)**2))*(hermite_term(2,1) - hermite_term(1,1)**2) + ((PB(1)*&
&*3) + 6.0*(PB(1)**2)*PA(1) + 3.0*PB(1)*(PA(1)**2))*(hermite_term(3,1) - hermit&
&e_term(2,1)*hermite_term(1,1)) + (PA(1)**2+ 3.0*(PB(1)**2) + 6.0*PB(1)*PA(1))*&
&(hermite_term(4,1) - hermite_term(3,1)*hermite_term(1,1)) + (2.0*PA(1) + 3.0*P&
&B(1))*(hermite_term(5,1) - hermite_term(4,1)*hermite_term(1,1)) + (hermite_ter&
&m(6,1) - hermite_term(5,1)*hermite_term(1,1)))*pc(1,1)

pc(19,18) = (PA(2) + hermite_term(1,2))*pc(6,18) + ((2.0*PA(2))*(hermite_term(&
&2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_&
&term(1,2)))*pc(1,18)

pc(20,18) = (PA(3) + hermite_term(1,3))*pc(7,18) + ((2.0*PA(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(1,18)

pc(1,19) = (PB(2) + hermite_term(1,2))*pc(1,6) + ((2.0*PB(2))*(hermite_term(2,&
&2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_te&
&rm(1,2)))*pc(1,1)

pc(2,19) = (PA(1) + hermite_term(1,1))*pc(1,19)

pc(3,19) = (PA(2) + hermite_term(1,2))*pc(1,19) + ((3.0*(PB(2)**2))*(hermite_t&
&erm(2,2) - hermite_term(1,2)**2) + (3*PB(2))*(hermite_term(3,2) - hermite_term&
&(2,2)*hermite_term(1,2)) + (hermite_term(4,2) - hermite_term(3,2)*hermite_term&
&(1,2)))*pc(1,1)

pc(4,19) = (PA(3) + hermite_term(1,3))*pc(1,19)

pc(5,19) = (PA(1) + hermite_term(1,1))*pc(2,19) + (hermite_term(2,1) - hermite&
&_term(1,1)**2)*pc(1,19)

pc(6,19) = (PA(2) + hermite_term(1,2))*pc(3,19) + ((3.0*(PB(2)**2)*PA(2) + PB(&
&2)**3)*(hermite_term(2,2) - hermite_term(1,2)**2) + (3.0*PB(2)*PA(2) + 3.0*(PB&
&(2)**2))*(hermite_term(3,2) - hermite_term(2,2)* hermite_term(1,2)) + (PA(2) +&
& 3.0*PB(2))*(hermite_term(4,2) - hermite_term(3,2)*hermite_term(1,2))+ (hermit&
&e_term(5,2)- hermite_term(4,2)*hermite_term(1,2)))*pc(1,1)

pc(7,19) = (PA(3) + hermite_term(1,3))*pc(4,19) + (hermite_term(2,3) - hermite&
&_term(1,3)**2)*pc(1,19)

pc(8,19) = (PA(2) + hermite_term(1,2))*pc(2,19) + ((3.0*(PB(2)**2))*(hermite_t&
&erm(2,2) - hermite_term(1,2)**2) + (3*PB(2))*(hermite_term(3,2) - hermite_term&
&(2,2)*hermite_term(1,2)) + (hermite_term(4,2) - hermite_term(3,2)*hermite_term&
&(1,2)))*pc(2,1)

pc(9,19) = (PA(3) + hermite_term(1,3))*pc(2,19)

pc(10,19) = (PA(3) + hermite_term(1,3))*pc(3,19)

pc(11,19) = (PA(3) + hermite_term(1,3))*pc(8,19)

pc(12,19) = (PA(2) + hermite_term(1,2))*pc(5,19) + ((3.0*(PB(2)**2))*(hermite_&
&term(2,2) - hermite_term(1,2)**2) + (3*PB(2))*(hermite_term(3,2) - hermite_ter&
&m(2,2)*hermite_term(1,2)) + (hermite_term(4,2) - hermite_term(3,2)*hermite_ter&
&m(1,2)))*pc(5,1)

pc(13,19) = (PA(3) + hermite_term(1,3))*pc(5,19)

pc(14,19) = (PA(1) + hermite_term(1,1))*pc(6,19)

pc(15,19) = (PA(3) + hermite_term(1,3))*pc(6,19)

pc(16,19) = (PA(1) + hermite_term(1,1))*pc(7,19)

pc(17,19) = (PA(2) + hermite_term(1,2))*pc(7,19) + ((3.0*(PB(2)**2))*(hermite_&
&term(2,2) - hermite_term(1,2)**2) + (3*PB(2))*(hermite_term(3,2) - hermite_ter&
&m(2,2)*hermite_term(1,2)) + (hermite_term(4,2) - hermite_term(3,2)*hermite_ter&
&m(1,2)))*pc(7,1)

pc(18,19) = (PA(1) + hermite_term(1,1))*pc(5,19) + ((2.0*PA(1))*(hermite_term(&
&2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_&
&term(1,1)))*pc(1,19)

pc(19,19) = (PA(2) + hermite_term(1,2))*pc(6,19) + ((2.0*(PB(2)**3)*PA(2) + 3.&
&0*(PA(2)**2)*(PB(2)**2))*(hermite_term(2,2) - hermite_term(1,2)**2) + ((PB(2)*&
&*3) + 6.0*(PB(2)**2)*PA(2) + 3.0*PB(2)*(PA(2)**2))*(hermite_term(3,2) - hermit&
&e_term(2,2)*hermite_term(1,2)) + (PA(2)**2+ 3.0*(PB(2)**2) + 6.0*PB(2)*PA(2))*&
&(hermite_term(4,2) - hermite_term(3,2)*hermite_term(1,2)) + (2.0*PA(2) + 3.0*P&
&B(2))*(hermite_term(5,2) - hermite_term(4,2)*hermite_term(1,2)) + (hermite_ter&
&m(6,2) - hermite_term(5,2)*hermite_term(1,2)))*pc(1,1)

pc(20,19) = (PA(3) + hermite_term(1,3))*pc(7,19) + ((2.0*PA(3))*(hermite_term(&
&2,3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_&
&term(1,3)))*pc(1,19)

pc(1,20) = (PB(3) + hermite_term(1,3))*pc(1,7) + ((2.0*PB(3))*(hermite_term(2,&
&3) - hermite_term(1,3)**2) + (hermite_term(3,3) - hermite_term(2,3)*hermite_te&
&rm(1,3)))*pc(1,1)

pc(2,20) = (PA(1) + hermite_term(1,1))*pc(1,20)

pc(3,20) = (PA(2) + hermite_term(1,2))*pc(1,20)

pc(4,20) = (PA(3) + hermite_term(1,3))*pc(1,20) + ((3.0*(PB(3)**2))*(hermite_t&
&erm(2,3) - hermite_term(1,3)**2) + (3*PB(3))*(hermite_term(3,3) - hermite_term&
&(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite_term(3,3)*hermite_term&
&(1,3)))*pc(1,1)

pc(5,20) = (PA(1) + hermite_term(1,1))*pc(2,20) + (hermite_term(2,1) - hermite&
&_term(1,1)**2)*pc(1,20)

pc(6,20) = (PA(2) + hermite_term(1,2))*pc(3,20) + (hermite_term(2,2) - hermite&
&_term(1,2)**2)*pc(1,20)

pc(7,20) = (PA(3) + hermite_term(1,3))*pc(4,20) + ((3.0*(PB(3)**2)*PA(3) + PB(&
&3)**3)*(hermite_term(2,3) - hermite_term(1,3)**2) + (3.0*PB(3)*PA(3) + 3.0*(PB&
&(3)**2))*(hermite_term(3,3) - hermite_term(2,3)* hermite_term(1,3)) + (PA(3) +&
& 3.0*PB(3))*(hermite_term(4,3) - hermite_term(3,3)*hermite_term(1,3))+ (hermit&
&e_term(5,3)- hermite_term(4,3)*hermite_term(1,3)))*pc(1,1)

pc(8,20) = (PA(2) + hermite_term(1,2))*pc(2,20)

pc(9,20) = (PA(3) + hermite_term(1,3))*pc(2,20) + ((3.0*(PB(3)**2))*(hermite_t&
&erm(2,3) - hermite_term(1,3)**2) + (3*PB(3))*(hermite_term(3,3) - hermite_term&
&(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite_term(3,3)*hermite_term&
&(1,3)))*pc(2,1)

pc(10,20) = (PA(3) + hermite_term(1,3))*pc(3,20) + ((3.0*(PB(3)**2))*(hermite_&
&term(2,3) - hermite_term(1,3)**2) + (3*PB(3))*(hermite_term(3,3) - hermite_ter&
&m(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite_term(3,3)*hermite_ter&
&m(1,3)))*pc(3,1)

pc(11,20) = (PA(3) + hermite_term(1,3))*pc(8,20) + ((3.0*(PB(3)**2))*(hermite_&
&term(2,3) - hermite_term(1,3)**2) + (3*PB(3))*(hermite_term(3,3) - hermite_ter&
&m(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite_term(3,3)*hermite_ter&
&m(1,3)))*pc(8,1)

pc(12,20) = (PA(2) + hermite_term(1,2))*pc(5,20)

pc(13,20) = (PA(3) + hermite_term(1,3))*pc(5,20) + ((3.0*(PB(3)**2))*(hermite_&
&term(2,3) - hermite_term(1,3)**2) + (3*PB(3))*(hermite_term(3,3) - hermite_ter&
&m(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite_term(3,3)*hermite_ter&
&m(1,3)))*pc(5,1)

pc(14,20) = (PA(1) + hermite_term(1,1))*pc(6,20)

pc(15,20) = (PA(3) + hermite_term(1,3))*pc(6,20) + ((3.0*(PB(3)**2))*(hermite_&
&term(2,3) - hermite_term(1,3)**2) + (3*PB(3))*(hermite_term(3,3) - hermite_ter&
&m(2,3)*hermite_term(1,3)) + (hermite_term(4,3) - hermite_term(3,3)*hermite_ter&
&m(1,3)))*pc(6,1)

pc(16,20) = (PA(1) + hermite_term(1,1))*pc(7,20)

pc(17,20) = (PA(2) + hermite_term(1,2))*pc(7,20)

pc(18,20) = (PA(1) + hermite_term(1,1))*pc(5,20) + ((2.0*PA(1))*(hermite_term(&
&2,1) - hermite_term(1,1)**2) + (hermite_term(3,1) - hermite_term(2,1)*hermite_&
&term(1,1)))*pc(1,20)

pc(19,20) = (PA(2) + hermite_term(1,2))*pc(6,20) + ((2.0*PA(2))*(hermite_term(&
&2,2) - hermite_term(1,2)**2) + (hermite_term(3,2) - hermite_term(2,2)*hermite_&
&term(1,2)))*pc(1,20)

pc(20,20) = (PA(3) + hermite_term(1,3))*pc(7,20) + ((2.0*(PB(3)**3)*PA(3) + 3.&
&0*(PA(3)**2)*(PB(3)**2))*(hermite_term(2,3) - hermite_term(1,3)**2) + ((PB(3)*&
&*3) + 6.0*(PB(3)**2)*PA(3) + 3.0*PB(3)*(PA(3)**2))*(hermite_term(3,3) - hermit&
&e_term(2,3)*hermite_term(1,3)) + (PA(3)**2+ 3.0*(PB(3)**2) + 6.0*PB(3)*PA(3))*&
&(hermite_term(4,3) - hermite_term(3,3)*hermite_term(1,3)) + (2.0*PA(3) + 3.0*P&
&B(3))*(hermite_term(5,3) - hermite_term(4,3)*hermite_term(1,3)) + (hermite_ter&
&m(6,3) - hermite_term(5,3)*hermite_term(1,3)))*pc(1,1)

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


   end subroutine Koverlap2CIntgAna

   subroutine Koverlap2CIntgNum(a1,a2,A,B,deltaK,pc,sh,cell_size,step_size)

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
   real (kind=double), dimension (3), intent (in) :: deltaK
   complex (kind=double), dimension (20,20), intent(out) :: pc
   complex (kind=double), dimension (16,16), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, h, i
   integer :: num_steps
   integer, dimension (20,3) :: triads
   integer, dimension (16,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos
   real (kind=double), dimension (3) :: xyz
   complex (kind=double), dimension (3) :: xyz_sum, xyz_soln

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


   ! Initialize numerical mesh quantities.
   start_pos = -cell_size
   num_steps = int(cell_size * 2.0d0 / step_size) + 1  ! +1 accounts for xyz=0

   ! Initialize a counter of the triad pq pairs.
   h = 0

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
                  & (xyz(:) - A(:))**l1(:) * (xyz(:) - B(:))**l2(:) * &
                  & exp(-a1*(xyz(:) - A(:))**2) * exp(-a2*(xyz(:) - &
                  & B(:))**2) * (cos(deltaK(:)*xyz(:)) - &
                  & cmplx(0.0d0,1.0d0,double) * sin(deltaK(:)*xyz(:)))
            xyz_sum(:) = xyz_sum(:) + xyz_soln(:)
         enddo

         pc(q,p) = product(xyz_sum(:))

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


   end subroutine Koverlap2CIntgNum

end program GaussianIntegrals
