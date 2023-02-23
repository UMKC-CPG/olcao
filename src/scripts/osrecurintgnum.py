#!/usr/bin/env python3

"""
Collection of subroutines that are used to (1) produce numerical solutions to
multicenter gaussian type orbital integrals; (2) print Fortran code that, when
executed, will evaluate the numerical solutions for use in a testing
environment.
"""

import osrecurintglib as lib


def print_test_overlap_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
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
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos
   real (kind=double), dimension (3) :: xyz, xyz_sum, xyz_soln

   ! Initialize local variables.
"""
    for i in range(len(triads)):
        head += f"   triads({i+1},:) = (/{triads[i][0]}," + \
                f"{triads[i][1]},{triads[i][2]}/)\n"

    for i in range(len(conversion)):
        head += f"   conversion({i+1},1,:) = (/{conversion[i][0][0]}," + \
                f"{conversion[i][0][1]},{conversion[i][0][2]}/)\n"
        head += f"   conversion({i+1},2,:) = (/{conversion[i][1][0]}," + \
                f"{conversion[i][1][1]},{conversion[i][1][2]}/)\n"

    head += """
   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ f"{len(triads)}" + """

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, 0, num_conversion, num_conversion)


    # Print the subroutine foot.
    foot = """
   end subroutine overlap2CIntgNum
"""
    f.write(foot)


def print_test_kinetic_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
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
   real (kind=double), dimension (""" + \
           f"{len(triads)},{len(triads)}" + """), intent(out) :: pc
   real (kind=double), dimension (""" + \
           f"{len(conversion)},{len(conversion)}" + """), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i, j, k
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
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
   !   momentum or we must discard this term when lx2-2 > 0.

   ! Initialize local variables.
"""

    for i in range(len(triads)):
        head += f"   triads({i+1},:) = (/{triads[i][0]}," + \
                f"{triads[i][1]},{triads[i][2]}/)\n"

    for i in range(len(conversion)):
        head += f"   conversion({i+1},1,:) = (/{conversion[i][0][0]}," + \
                f"{conversion[i][0][1]},{conversion[i][0][2]}/)\n"
        head += f"   conversion({i+1},2,:) = (/{conversion[i][1][0]}," + \
                f"{conversion[i][1][1]},{conversion[i][1][2]}/)\n"

    head += """
   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, 0, num_conversion, num_conversion)


    # Print the subroutine foot and the subsequent functions.
    foot = """
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
"""
    f.write(foot)


def print_test_electron_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
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
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos
   real (kind=double), dimension (3) :: xyz, xyz_sum, xyz_soln

   ! Initialize local variables.
"""

    for i in range(len(triads)):
        head += f"   triads({i+1},:) = (/{triads[i][0]}," + \
                f"{triads[i][1]},{triads[i][2]}/)\n"

    for i in range(len(conversion)):
        head += f"   conversion({i+1},1,:) = (/{conversion[i][0][0]}," + \
                f"{conversion[i][0][1]},{conversion[i][0][2]}/)\n"
        head += f"   conversion({i+1},2,:) = (/{conversion[i][1][0]}," + \
                f"{conversion[i][1][1]},{conversion[i][1][2]}/)\n"

    head += """
   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, 0, num_conversion, num_conversion)


    # Print the subroutine foot.
    foot = """
   end subroutine electron3CIntgNum
"""
    f.write(foot)


def print_test_nuclear_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
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
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i, j, k
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)},3) :: triads" + """
   integer, dimension (""" + f"{len(conversion)},2,3) :: conversion" + """
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, r, soln
   real (kind=double) :: xC_dist_sqrd, yC_dist_sqrd, zC_dist_sqrd
   real (kind=double), dimension (3) :: xyz, xyz_sum

   ! Initialize local variables.
"""

    for i in range(len(triads)):
        head += f"   triads({i+1},:) = (/{triads[i][0]}," + \
                f"{triads[i][1]},{triads[i][2]}/)\n"

    for i in range(len(conversion)):
        head += f"   conversion({i+1},1,:) = (/{conversion[i][0][0]}," + \
                f"{conversion[i][0][1]},{conversion[i][0][2]}/)\n"
        head += f"   conversion({i+1},2,:) = (/{conversion[i][1][0]}," + \
                f"{conversion[i][1][1]},{conversion[i][1][2]}/)\n"

    head += """
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

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, 0, num_conversion, num_conversion)


    # Print the subroutine foot.
    foot = """
   end subroutine nuclear3CIntgNum
"""
    f.write(foot)


def print_test_momentum_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
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
   real (kind=double), dimension (""" \
    + f"{len(triads)},{len(triads)},3), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
    + f"{len(conversion)},{len(conversion)},3), intent(out) :: sh" + """
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i, j
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)},3) :: triads" + """
   integer, dimension (""" + f"{len(conversion)},2,3) :: conversion" + """
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
   !   momentum or we must discard this term when lx2-1 > 0.

   ! Initialize local variables.
"""

    for i in range(len(triads)):
        head += f"   triads({i+1},:) = (/{triads[i][0]}," + \
                f"{triads[i][1]},{triads[i][2]}/)\n"

    for i in range(len(conversion)):
        head += f"   conversion({i+1},1,:) = (/{conversion[i][0][0]}," + \
                f"{conversion[i][0][1]},{conversion[i][0][2]}/)\n"
        head += f"   conversion({i+1},2,:) = (/{conversion[i][1][0]}," + \
                f"{conversion[i][1][1]},{conversion[i][1][2]}/)\n"

    head += """
   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, 1, num_conversion, num_conversion)
    lib.print_pc_to_sh(conversion, f, 2, num_conversion, num_conversion)
    lib.print_pc_to_sh(conversion, f, 3, num_conversion, num_conversion)


    # Print the subroutine foot.
    foot = """
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
      !   and the expression in the last equals of equation 69 in Ben Walker's
      !   dissertation.

      ! Second line of equation 69. Also visible as second term in the
      !   script-produced equation.
      primeMM = primeMM - 2.0d0*a2 * (curr_pos - B)**(l2+1)

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
"""
    f.write(foot)


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. In this case, there is no program. This
    #   should only be imported as a module.
    main()
