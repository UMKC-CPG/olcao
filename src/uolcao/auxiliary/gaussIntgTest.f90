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
   allocate (pc(4,4,3,2))
   allocate (sh(4,4,3,2))
   allocate (pcCmplx(4,4,3,2))
   allocate (shCmplx(4,4,3,2))

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
      real(kind=double), dimension(4,4,2) :: pc ! Last idx 1 = ana; 2 = num
      real(kind=double), dimension(4,4,2) :: sh ! Last idx 1 = ana; 2 = num
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
      endif

      if ((num_segments > 1) .or. (num_steps > 1)) then

         ! Print a header line.
         if (curr_line == 1) then

            ! Start with the parameters.
            write (299+unitPlus,fmt="(a)",advance="NO") "IDX SEG STEP "
            write (299+unitPlus,fmt="(a)",advance="NO") "alpha1 alpha2 alpha3 "
            write (299+unitPlus,fmt="(a)",advance="NO") "Ax Ay Az "
            write (299+unitPlus,fmt="(a)",advance="NO") "Bx By Bz "
            write (299+unitPlus,fmt="(a)",advance="NO") "Cx Cy Cz "
            write (299+unitPlus,fmt="(a)",advance="NO") "dltKx dltKy dltKz "

            ! Now add the pc names.
            do j = 1, 4
               do k = 1, 4
                  write (299+unitPlus,fmt="(4a)",advance="NO") &
                        & trim(pcName(k)),"_",trim(pcName(j))," "
               enddo
            enddo

            ! Now add the sh names.
            do j = 1, 4
               do k = 1, 4
                  write (299+unitPlus,fmt="(4a)",advance="NO") &
                        & trim(shName(k)),"_",trim(shName(j))," "
               enddo
            enddo

            ! End the header line.
            write (299+unitPlus,*)
         endif

         ! Print the data line.

         ! Print the line number.
         write (299+unitPlus,fmt="(i5)",advance="NO") curr_line

         write (299+unitPlus,fmt="(2i5)",advance="NO") h, i
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") alphas(:)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,1)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,2)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,3)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") deltaK(:)

         ! Print the result difference for the given pc matrix.
         do j = 1, 4
            do k = 1, 4
               write (299+unitPlus,fmt="(e16.8)",advance="NO") &
                     & pc(k,j,1) - pc(k,j,2)
            enddo
         enddo

         ! Print the result difference for the given sh matrix.
         do j = 1, 4
            do k = 1, 4
               write (299+unitPlus,fmt="(e16.8)",advance="NO") &
                     & sh(k,j,1) - sh(k,j,2)
            enddo
         enddo

         ! Write an endline to prepare for the next step.
         write (299+unitPlus, *) ""

      else

         write (299+unitPlus,fmt="(3e16.8)",advance="NO") alphas(:)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,1)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,2)
         write (299+unitPlus,fmt="(3e16.8)",advance="NO") pos(:,3)
         write (299+unitPlus,fmt="(3e16.8)") deltaK(:)

         ! Print the result difference for the given pc matrix.
         do j = 1, 4
            do k = 1, 4
               write (299+unitPlus,fmt="(2a4)",advance="NO") &
                     & pcName(k),pcName(j)
               write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,1)
               write (299+unitPlus,fmt="(e16.8)",advance="NO") pc(k,j,2)
               write (299+unitPlus,fmt="(e16.8)") &
                     & pc(k,j,1) - pc(k,j,2)
            enddo
         enddo

         ! Print the result difference for the given sh matrix.
         do j = 1, 4
            do k = 1, 4
               write (299+unitPlus,fmt="(2a16)",advance="NO") &
                     & shName(k),shName(j)
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
   complex (kind=double), dimension (4,4), intent(out) :: pc
   complex (kind=double), dimension (4,4), intent(out) :: sh

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

pc(1,2) = (PB(1) + hermite_term(1,1))*preFactorKO(1)*preFactorKO(2)*preFactorK&
&O(3)

pc(2,2) = (PA(1) + hermite_term(1,1))*pc(1,2) + (hermite_term(2,1) - hermite_t&
&erm(1,1)**2)*pc(1,1)

pc(3,2) = (PA(2) + hermite_term(1,2))*pc(1,2)

pc(4,2) = (PA(3) + hermite_term(1,3))*pc(1,2)

pc(1,3) = (PB(2) + hermite_term(1,2))*preFactorKO(1)*preFactorKO(2)*preFactorK&
&O(3)

pc(2,3) = (PA(1) + hermite_term(1,1))*pc(1,3)

pc(3,3) = (PA(2) + hermite_term(1,2))*pc(1,3) + (hermite_term(2,2) - hermite_t&
&erm(1,2)**2)*pc(1,1)

pc(4,3) = (PA(3) + hermite_term(1,3))*pc(1,3)

pc(1,4) = (PB(3) + hermite_term(1,3))*preFactorKO(1)*preFactorKO(2)*preFactorK&
&O(3)

pc(2,4) = (PA(1) + hermite_term(1,1))*pc(1,4)

pc(3,4) = (PA(2) + hermite_term(1,2))*pc(1,4)

pc(4,4) = (PA(3) + hermite_term(1,3))*pc(1,4) + (hermite_term(2,3) - hermite_t&
&erm(1,3)**2)*pc(1,1)

sh(1,1) = pc(1,1)

sh(2,1) = pc(2,1)

sh(3,1) = pc(3,1)

sh(4,1) = pc(4,1)

sh(1,2) = pc(1,2)

sh(2,2) = pc(2,2)

sh(3,2) = pc(3,2)

sh(4,2) = pc(4,2)

sh(1,3) = pc(1,3)

sh(2,3) = pc(2,3)

sh(3,3) = pc(3,3)

sh(4,3) = pc(4,3)

sh(1,4) = pc(1,4)

sh(2,4) = pc(2,4)

sh(3,4) = pc(3,4)

sh(4,4) = pc(4,4)


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
   complex (kind=double), dimension (4,4), intent(out) :: pc
   complex (kind=double), dimension (4,4), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, h, i
   integer :: num_steps
   integer, dimension (4,3) :: triads
   integer, dimension (4,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos
   real (kind=double), dimension (3) :: xyz
   complex (kind=double), dimension (3) :: xyz_sum, xyz_soln

   ! Initialize local variables.
   triads(1,:) = (/0,0,0/)
   triads(2,:) = (/1,0,0/)
   triads(3,:) = (/0,1,0/)
   triads(4,:) = (/0,0,1/)
   conversion(1,1,:) = (/1,0,0/)
   conversion(1,2,:) = (/1,1,1/)
   conversion(2,1,:) = (/1,0,0/)
   conversion(2,2,:) = (/2,1,1/)
   conversion(3,1,:) = (/1,0,0/)
   conversion(3,2,:) = (/3,1,1/)
   conversion(4,1,:) = (/1,0,0/)
   conversion(4,2,:) = (/4,1,1/)


   ! Initialize numerical mesh quantities.
   start_pos = -cell_size
   num_steps = int(cell_size * 2.0d0 / step_size) + 1  ! +1 accounts for xyz=0

   ! Initialize a counter of the triad pq pairs.
   h = 0

   do p = 1, 4
      do q = 1, 4

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

sh(1,2) = pc(1,2)

sh(2,2) = pc(2,2)

sh(3,2) = pc(3,2)

sh(4,2) = pc(4,2)

sh(1,3) = pc(1,3)

sh(2,3) = pc(2,3)

sh(3,3) = pc(3,3)

sh(4,3) = pc(4,3)

sh(1,4) = pc(1,4)

sh(2,4) = pc(2,4)

sh(3,4) = pc(3,4)

sh(4,4) = pc(4,4)


   end subroutine Koverlap2CIntgNum

end program GaussianIntegrals
