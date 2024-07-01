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
   integer :: p, q, h, i
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

   ! Initialize a counter of the triad pq pairs.
   h = 0

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [1, 1, ""], 0, num_conversion,
                       num_conversion)


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
   integer :: p, q, h, i, j, k
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, curr_pos
   real (kind=double), dimension (3) :: xyz, xyz_soln
   real (kind=double), dimension (3,2) :: xyz_I ! Indices=xyz, prime||noprime

   ! Before we proceed with the calculation we need to understand a bit more
   !   about exactly what is being computed. The form of the integration is:
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        -1/2 * dell^2 [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   ! We use Pxyz for points in space; Rxyz1,2 for atomic sites 1 and 2;
   !   zeta1,2 for the exponential decay factors; and lxyz1,2 for the angular
   !   momentum numbers for atoms 1 and 2.
   ! The dell^2 operator is d^2/dx^2 + d^2/dy^2 + d^2/dz^2.
   ! The integral must be computed over all space for all different possible
   !   values of lxyz1,2 for some arbitrarily chosen test coordinates and
   !   decay rates.
   ! Because of the plus signs in the dell^2 operator we arrive at three
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
   !   Ix' = S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * -1/2 *
   !          [-2*zeta2 + 4*zeta2**2*(Px-Rx2)**2] *
   !          exp(-zeta2*(Px-Rx2)**2)
   !        ]
   ! Each of the other integrals (Iy, Iz) will have the form of a simple
   !   1D overlap integral:
   !   Iy = S [(Py-Ry1)**ly1 * exp(-zeta1*((Py-Ry1)**2)) *
   !            (Py-Ry2)**ly2 * exp(-zeta2*((Py-Ry2)**2))] dy
   !   Iz = S [(Pz-Rz1)**lz1 * exp(-zeta1*((Pz-Rz1)**2)) *
   !            (Pz-Rz2)**lz2 * exp(-zeta2*((Pz-Rz2)**2))] dz
   ! The total integral !! for the d^2/dx^2 version !! is thus Ix' * Iy * Iz.
   !   The total integral (including all terms of dell^2) will have the form:
   !   Ix'*Iy*Iz + Ix*Iy'*Iz + Ix*Iy*Iz' where the Iy' and Iz' are the
   !   appropriate analogs of the Ix' and the Ix is the analog of the Iy
   !   or Iz.
   ! Thus, every integral is an independent 1D integral and the solution of
   !   the total integral is a sum of products of integrals that can all be
   !   computed at once.
   ! With regard to the term in each integral that has the **(lx2-2) form
   !   (or **(ly2-2) or **(lz2-2)), some special care must be taken. This
   !   term represents an angular momentum and the integral will have the
   !   form of an overlap integral. Because we cannot have a negative angular
   !   momentum we must discard this term when lx2-2 < 0.

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

   ! Initialize a counter of the triad pq pairs.
   h = 0

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

         ! Assign l1 and l2 values for each primitive Cartesian gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize sum variables.
         xyz_I(:,:) = 0.0d0

         ! Start a loop over 1D coordinates.
         do i = 0, num_steps
            curr_pos = start_pos + (i*step_size)

            ! Compute the primed integrals first.
            do j = 1, 3
               xyz_I(j,1) = xyz_I(j,1) &
                     & + primeKE(step_size, curr_pos, A(j), &
                     & B(j), a1, a2, l1(j), l2(j))
            enddo

            ! Compute the no-prime integrals second.
            do j = 1, 3
               xyz_I(j,2) = xyz_I(j,2) &
                     & + noPrimeKE(step_size, curr_pos, A(j), &
                     & B(j), a1, a2, l1(j), l2(j))
            enddo
         enddo

         pc(q,p) = &
               & xyz_I(1,1)*xyz_I(2,2)*xyz_I(3,2) + &
               & xyz_I(1,2)*xyz_I(2,1)*xyz_I(3,2) + &
               & xyz_I(1,2)*xyz_I(2,2)*xyz_I(3,1)

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [2, 1, ""], 0, num_conversion,
                       num_conversion)


    # Print the subroutine foot and the subsequent functions.
    foot = """
   end subroutine kinetic2CIntgNum


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

      ! Compute the integral part.
      noPrimeKE = step_size * (curr_pos - A)**l1 * (curr_pos - B)**l2 &
            & * exp(-a1*(curr_pos - A)**2) * exp(-a2*(curr_pos - B)**2)

      return

   end function noPrimeKE


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

      ! Compute each "internal" term of the prime integral. Compare each
      !   of these lines with the 1D kinetic energy equation produce by the
      !   osrecurintg_makenum.py script (appropriately separated into terms)
      !   and the expression in the last equals of equation 65 in Ben Walker's
      !   dissertation. Note that a factor of -1/2 is already included in all
      !   the below expressions. Because only one "primeKE" term is present
      !   in each term of the KE integral, the -1/2 coefficient will be
      !   properly accounted for. (Perhaps this isn't the clearest way to
      !   do this.)

      ! Combination of terms 1 and 2 in line 1 of last equals of eqn 65. Also
      !   visible as combination of first two terms in the main parentheses
      !   of the algebraically ordered 1D kinetic energy equation produced by
      !   osrecurintg_makenum.py.
      primeKE = (a2 * (2*l2 + 1)) * (curr_pos - B)**l2

      ! Line 2 of last equals in eqn 65. Also visible as third term in the
      !   script produced equation.
      primeKE = primeKE - 2.0d0*a2**2 * (curr_pos - B)**(l2+2)

      ! Combination of terms 1 and 2 in line 3. Also visible as the
      !   combination of the 4th and 5th terms in the script produced eqn.
      !   As mentioned above, we only compute this if the angular momentum
      !   will allow it.
      if (l2 >= 2) then
         primeKE = primeKE - 0.5 * l2*(l2-1) * (curr_pos - B)**(l2-2)
      endif

      ! Multiply the prime integral by the preceeding primitive gaussian
      !   coefficient and exponential and multiply by the succeeding
      !   exponential. (We have already multiplied by the succeeding
      !   primitive gaussian coefficient in the above lines.)
      primeKE = primeKE * (curr_pos-A)**l1 * exp(-a1*(curr_pos-A)**2) &
            & * exp(-a2*(curr_pos-B)**2)

      ! Finally, multiply by the step size.
      primeKE = primeKE * step_size

      return
      
   end function primeKE
"""
    f.write(foot)


# The "full" test subroutine is intended to only be used for solving
#   integrals of the form (0|T|0) to guide the creation of the kinetic
#   energy pre-factor for use in the analytic solutions. Here, we find
#   that the equation we produce for the KE (0|T|0) expression differs
#   from that given in the Obara-Saika paper. The numerical values are
#   identical though. I have not checked yet, but probably the issue
#   is that the form here works the derivative on only the B site while
#   I'm guessing that the Obara-Saika form takes the average of both the
#   A and B site solutions to make a "unbiased" expression. Numerically,
#   all approaches lead to the same result though. It is just a matter
#   of presentation.
# Note that the numerical solution here is "full" because it does a full
#   3D numerical integral instead of the usual 3 1D integrals.
def print_test_kinetic_full_num(conversion, triads, f):

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
   integer :: p, q, h, i, j, k
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, curr_pos
   real (kind=double), dimension (3) :: xyz, xyz_soln
   real (kind=double), dimension (3,2) :: xyz_I ! Indices=xyz, prime||noprime

   ! Full solution variables.
   real (kind=double) :: soln
   real (kind=double) :: sum_xyz_A_2, sum_xyz_B_2
   real (kind=double) :: zeta, xi, preFactorOL
   real (kind=double), dimension (3) :: d, S, SB


   ! Before we proceed with the calculation we need to understand a bit more
   !   about exactly what is being computed. The form of the integration is:
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        -1/2 * dell^2 [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   ! We use Pxyz for points in space; Rxyz1,2 for atomic sites 1 and 2;
   !   zeta1,2 for the exponential decay factors; and lxyz1,2 for the angular
   !   momentum numbers for atoms 1 and 2.
   ! The dell^2 operator is d^2/dx^2 + d^2/dy^2 + d^2/dz^2.
   ! The integral must be computed over all space for all different possible
   !   values of lxyz1,2 for some arbitrarily chosen test coordinates and
   !   decay rates.
   ! Because of the plus signs in the dell^2 operator we arrive at three
   !   identical integrals of the form (with # in d^2/d#^2 = x, y, or z):
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        -1/2 * d^2/d#^2 [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   ! Now, focusing on just the d^2/dx^2 version of the set of three integrals,
   !   we will pull all terms with Py and Pz out of the dx integral and
   !   because we are focusing on the (0|T|0) integral we will set all angular
   !   momentum variables (lx1, lx2, ly1, ly2, lz1, lz2) to zero to get:
   ! SS { [exp(-zeta1*((Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !      [exp(-zeta2*((Py-Ry2)**2 + (Pz-Rz2)**2 ))] *
   !     S [exp(-zeta1*((Px-Rx1)**2)) * -1/2 * d^2/dx^2 [
   !      exp(-zeta2*((Px-Rx2)**2))]] dx
   !    } dydz
   ! Applying the derivative, the internal 1D dx integral has the form:
   !   Ix' = S [exp(-zeta1*((Px-Rx1)**2)) * -1/2 *
   !          [-2*zeta2 + 4*zeta2**2*(Px-Rx2)**2] *
   !          exp(-zeta2*(Px-Rx2)**2)
   !        ]
   ! Each of the other integrals (Iy, Iz) will have the form of a simple
   !   1D overlap integral:
   !   Iy = S [exp(-zeta1*((Py-Ry1)**2)) *
   !            exp(-zeta2*((Py-Ry2)**2))] dy
   !   Iz = S [exp(-zeta1*((Pz-Rz1)**2)) *
   !            exp(-zeta2*((Pz-Rz2)**2))] dz
   ! The total integral !! for the d^2/dx^2 version !! is thus Ix' * Iy * Iz.
   !   The total integral (including all terms of dell^2) will have the form:
   !   Ix'*Iy*Iz + Ix*Iy'*Iz + Ix*Iy*Iz' where the Iy' and Iz' are the
   !   appropriate analogs of the Ix' and the Ix is the analog of the Iy
   !   or Iz.
   ! However, for the purpose of developing an analytical form for the (0|T|0)
   !   case, we will include all parts of the integral and do a 3D mesh using
   !   nested loops.

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
   ! Initialize terms used to compute analytical factors.
   zeta = a1 + a2
   xi = a1 * a2 / zeta
   d = A - B
   S = (a1*A + a2*B) / zeta
   SB = S - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d*d))
   preFactorKE = a2 * (3.0d0 - 2*a2 * (sum(SB(:)**2) + 3/(2.0d0*zeta)))

   ! Initialize numerical mesh quantities.
   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   ! Initialize a counter of the triad pq pairs.
   h = 0

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

         ! Assign l1 and l2 values for each primitive Cartesian gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize sum variables.
         soln = 0.0d0

         ! Start a loop over 1D coordinates.
         do i = 0, num_steps
            ! Compute the current x position.
            xyz(1) = (start_pos + (i*step_size))

            do j = 0, num_steps
            
               ! Compute the current y position.
               xyz(2) = (start_pos + (j*step_size))

               do k = 0, num_steps

                  ! Compute the current z position.
                  xyz(3) = (start_pos + (k*step_size))

                  sum_xyz_A_2 = sum((xyz(:)-A(:))**2)
                  sum_xyz_B_2 = sum((xyz(:)-B(:))**2)

                  soln = soln + &
& exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) * preFactorKE
! The above solution is the compact form of an overlap (0||0)
!   integral with the mesh point independent preFactor for
!   kinetic energy integrals obtained through application of the
!   Obara-Saika scheme. Specifically, for the xx term:
!   (0||xx) = SB*(0||x) + 1/2zeta * (0||0)
!           = SB**2 * (0||0) + 1/2zeta * (0||0)
!           = (SB**2 + 1/2zeta) * (0||0)
!   The same is repeated for the (0||yy) and (0||zz) terms leading to
!   (0||xx + yy + zz) = (sum(SB(:)**2) + 3/2zeta) * (0||0)

! Slight rearrangement and condensation of the below solution.
!exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) * ( &
!& (4*a2**2*(xyz(1)-B(1))**2 - 2*a2) + &
!& (4*a2**2*(xyz(2)-B(2))**2 - 2*a2) + &
!& (4*a2**2*(xyz(3)-B(3))**2 - 2*a2))

! Original solution to integral (not including the -1/2 factor) from sympy.
!(a2**2*(2*xyz(1)-2*B(1))**2-2*a2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)+(a&
!&2**2*(2*xyz(2)-2*B(2))**2-2*a2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)+(a2*&
!&*2*(2*xyz(3)-2*B(3))**2-2*a2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)
               enddo
            enddo
         enddo

         ! Multiply the results by the mesh step size.
         soln = soln * step_size**3

         pc(q,p) = soln * -0.5  ! Apply -1/2 factor.

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [2, 1, ""], 0, num_conversion,
                       num_conversion)


    # Print the subroutine foot and the subsequent functions.
    foot = """
   end subroutine kinetic2CIntgNum
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
   integer :: p, q, h, i
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

   ! Initialize a counter of the triad pq pairs.
   h = 0

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [3, 1, ""], 0, num_conversion,
                       num_conversion)


    # Print the subroutine foot.
    foot = """
   end subroutine electron3CIntgNum
"""
    f.write(foot)


def print_test_nuclear_num1dFast(conversion, triads, f):

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
           + f"{len(triads)},{len(triads)},1), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, h, i, j, k
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)},3) :: triads" + """
   integer, dimension (""" + f"{len(conversion)},2,3) :: conversion" + """
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, r, soln, xyz_sum_coeff
   real (kind=double), allocatable, dimension (:,:) :: xyz, xyz_sum
   real (kind=double), allocatable, dimension (:,:) :: C_dist, C_dist_sqrd

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

   allocate(xyz(num_steps,3))
   allocate(xyz_sum(num_steps,3))
   allocate(C_dist(num_steps,3))
   allocate(C_dist_sqrd(num_steps,3))

   ! Original equation from sympy:
   !   (xyz(1)-A(1))**l1(1)*(xyz(1)-B(1))**l2(1)*(xyz(2)-A(2))**l1(2)*
   !   (xyz(2)-B(2))**l2(2)*(xyz(3)-A(3))**l1(3)*(xyz(3)-B(3))**l2(3)*
   !   exp(-a1*sum((xyz(:)-A(:))**2)) * &
   !   exp(-a2*sum((xyz(:)-B(:))**2)) * &
   !   exp(-a3*sum((xyz(:)-C(:))**2)) / &
   !   sqrt(sum((xyz(:)-C(:))**2))

   ! By-hand algebraic rearrangement to the expression below for faster
   !   computation in the triple nested loop environment.

   ! Initialize a counter of the triad pq pairs.
   h = 0

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

         ! Assign l1 and l2 values for each gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Precompute all one dimensional distance and exp factors.
         do i = 1, 3
             do j = 1, num_steps+1
                xyz(j,i) = (start_pos + ((j-1)*step_size))
                C_dist(j,i) = xyz(j,i) - C(i)
                C_dist_sqrd(j,i) = (xyz(j,i) - C(i))**2
                xyz_sum(j,i) = (xyz(j,i)-A(i))**l1(i)*(xyz(j,i)-B(i))**l2(i) &
                    & *exp(-a1*((xyz(j,i)-A(i))**2) - a2*((xyz(j,i)-B(i))**2) &
                    & - a3*(C_dist(j,i))**2)
             enddo
         enddo

         ! Initialize the solution
         soln = 0.0d0

         do i = 0, num_steps
            do j = 0, num_steps
               do k = 0, num_steps

                  ! Compute the distance between the electron and the nucleus.
                  r = C_dist_sqrd(i,1) + C_dist_sqrd(j,2) + C_dist_sqrd(k,3)

                  ! If the distance is zero, then cycle.
                  if (r == 0) cycle

                  ! Compute the current z xyz_sum.
                  xyz_sum_coeff = xyz_sum(i,1)*xyz_sum(j,2)*xyz_sum(k,3)

                  soln = soln + xyz_sum_coeff / sqrt(r)
               enddo
            enddo
         enddo

         ! Multiply the results by the mesh step size.
         soln = soln * step_size**3

         pc(q,p,1) = soln

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

   deallocate(xyz)
   deallocate(xyz_sum)
   deallocate(C_dist)
   deallocate(C_dist_sqrd)

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [4, 1, ""], 0, num_conversion,
                       num_conversion)


    # Print the subroutine foot.
    foot = """
   end subroutine nuclear3CIntgNum
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
           + f"{len(triads)},{len(triads)},1), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, h, i, j, k
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

   ! Initialize a counter of the triad pq pairs.
   h = 0

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

         pc(q,p,1) = soln

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [4, 1, ""], 0, num_conversion,
                       num_conversion)


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
   integer :: p, q, h, i, j
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
   !   momentum we must discard this term when lx2-1 < 0.

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

   ! Initialize a counter of the triad pq pairs
   h = 0

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [6, 3, ""], 0, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [6, 3, ""], 1, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [6, 3, ""], 2, num_conversion,
                       num_conversion)


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
"""
    f.write(foot)


def print_test_dipole_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
   subroutine dipole3CIntgNum(a1,a2,A,B,C,pc,sh,cell_size,step_size)

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
   real (kind=double), dimension (3), intent (in) :: A, B, C
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)},3), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)},3), intent(out) :: sh" + """
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables
   integer :: p, q, h, i
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)},3) :: triads" + """
   integer, dimension (""" + f"{len(conversion)},2,3) :: conversion" + """
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, xyz
   real (kind=double), dimension (3) :: soln, xyz_sol, xyz_sol_mu
   real (kind=double), dimension (3) :: mu_soln, xyz_sum, xyz_sum_mu

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
   !(xyz(1)-A(1))**l1(1)*(xyz(1)-B(1))**l2(1)
   !   *(xyz(2)-A(2))**l1(2)*(xyz(2)-B(2))**l2(2)
   !   *(xyz(3)-A(3))**l1(3)*(xyz(3)-B(3))**l2(3)
   !   *((xyz(1)-C(1))**mu_x(1)*(xyz(2)-C(2))**mu_y(2)
   !   * (xyz(3)-C(3))**mu_z(3))
   !   *exp(-a1*sum((xyz(:)-A(:))**2))
   !   *exp(-a2*sum((xyz(:)-B(:))**2))

   ! Because the expression can be separated along cartesian coordinates
   !   we can rewrite the integral as a product of three 1D integrals instead
   !   of needing to do a single (much more expensive) 3D integral.
   ! Note: Each x, y, z dipole moment solution is a product of three 1D
   !   integrals where each 1D integral is used more than once.

   ! Initialize a counter of the number of triad pairs completed.
   h = 0

   do p = 1,""" + f"{len(triads)}" + """ 
      do q = 1,""" + f"{len(triads)}" + """

         ! Assign l1 and l2 values for each gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize the solution
         soln(:) = 0.0d0
         xyz_sum(:) = 0.0d0
         xyz_sum_mu(:) = 0.0d0

         do i = 0, num_steps

            ! Compute the current x position, distance and sum.
            xyz = (start_pos + (i*step_size))

            xyz_sol_mu(:) = step_size * ((xyz - A(:))**l1(:)) &
                  & * ((xyz - B(:))**l2(:)) &
                                & * exp(-a1*(xyz-A(:))**2) * &
                                & exp(-a2*(xyz-B(:))**2) * &
                                & (xyz-C(:))

            xyz_sol(:) = step_size * ((xyz - A(:))**l1(:)) &
                  & *((xyz - B(:))**l2(:)) &
                  & * exp(-a1*(xyz-A(:))**2) * &
                  & exp(-a2*(xyz-B(:))**2)

            xyz_sum_mu(:) = xyz_sum_mu(:) + xyz_sol_mu(:)

            xyz_sum(:) = xyz_sum(:) + xyz_sol(:)

         enddo

         pc(q,p,1) = xyz_sum_mu(1)*xyz_sum(2)*xyz_sum(3)
         pc(q,p,2) = xyz_sum(1)*xyz_sum_mu(2)*xyz_sum(3)
         pc(q,p,3) = xyz_sum(1)*xyz_sum(2)*xyz_sum_mu(3)

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

"""

    f.write(head)
        
    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [7, 3, ""], 0, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [7, 3, ""], 1, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [7, 3, ""], 2, num_conversion,
                       num_conversion)


    # Print the subroutine foot.
    foot = """
   end subroutine dipole3CIntgNum
"""
    f.write(foot)


# The "full" test subroutine is intended to only be used for solving
#   integrals of the form (0|MV|0) to guide the creation of the mass
#   velocity pre-factor for use in the analytic solutions.
# Note that the numerical solution here is "full" because it does a full
#   3D numerical integral instead of the usual 3 1D integrals. The key
#   idea is that for the (0|MV|0) case, we show that a full 3D integral
#   can also be expressed as a full 3D integral of a simple (0||0) integral
#   times an appropriate preFactor and how that preFactor can be derived
#   so that it can be used in the analytic solution.
def print_test_massvel_full_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
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
   real (kind=double), dimension (""" + \
           f"{len(triads)},{len(triads)}" + """), intent(out) :: pc
   real (kind=double), dimension (""" + \
           f"{len(conversion)},{len(conversion)}" + """), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, h, i, j, k
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, coeff
   real (kind=double), dimension (2) :: curr_pos
   real (kind=double), dimension (3) :: xyz, xyz_soln
   real (kind=double), dimension (3,3) :: xyz_I
      ! The second index is for prime, noprime, 2Dprime.
      ! The first index is xyz for prime and noprime while it is xy,xz,yz for
      !   the 2Dprime case.

   ! Full solution variables.
   real (kind=double) :: soln
   real (kind=double) :: sum_xyz_A_2, sum_xyz_B_2, sum_xyz_B_4
   real (kind=double) :: product_xyz_A_l1, product_xyz_B_l2
   real (kind=double) :: zeta, xi
   real (kind=double) :: preFactorOL, preFactor02, preFactor04, preFactor22
   real (kind=double), dimension (3) :: d, S, SB

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
   !   we will pull all terms with Py and Pz out of the dx integral and we
   !   set all the angular momentum factors (lx1, lx2, ly1, ly2, lz1, lz2)
   !   to zero to get:
   ! SS { [exp(-zeta1*((Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !      [exp(-zeta2*((Py-Ry2)**2 + (Pz-Rz2)**2 ))] *
   !     S [exp(-zeta1*((Px-Rx1)**2)) * coeff * d^4/dx^4 [
   !      exp(-zeta2*((Px-Rx2)**2))]] dx
   !    } dydz
   ! Applying the fourth derivative, the internal 1D dx integral has the form:
   !   Ix' = S [exp(-zeta1*((Px-Rx1)**2)) * coeff *
   !          [zeta2^2*2 * (Px-Rx2) + zeta2^3*-48 * (Px-Rx2)^2 +
   !           (zeta2^4*16)*(Px-Rx2)^4] * exp(-zeta2*(Px-Rx2)**2)
   !        ]
   ! Each of the other integrals (Iy, Iz) will have the form of a simple
   !   1D overlap integral:
   !   Iy = S [exp(-zeta1*((Py-Ry1)^2)) *
   !            exp(-zeta2*((Py-Ry2)^2))] dy
   !   Iz = S [exp(-zeta1*((Pz-Rz1)^2)) *
   !            exp(-zeta2*((Pz-Rz2)^2))] dz
   ! The total integral !! for the d^4/dx^4 portion !! is thus Ix' * Iy * Iz.
   !   The total of the first part of the integral (including all terms of
   !   dell^4) will have the form: Ix'*Iy*Iz + Ix*Iy'*Iz + Ix*Iy*Iz' where
   !   the Iy' and Iz' are the appropriate analogs of the Ix' and the Ix is
   !   the analog of the Iy or Iz.
   !
   ! Now, focusing on just the 2 * d^2/d#^2 d^2/d@^2 where (say) #,@ = x,y,
   !   we will pull terms with Pz out of the dx dy integral and again setting
   !   all angular momentum values to zero we get:
   ! S { [exp(-zeta1*(Pz-Rz1)^2)] *
   !     [exp(-zeta2*(Pz-Rz2)^2)] *
   !   SS [exp(-zeta1*((Px-Rx1)^2 + (Py-Ry1)^2)) * coeff *
   !      d^2/dx^2 d^2/dy^2 [
   !      exp(-zeta2*((Px-Rx2)^2 + (Py-Ry2)**2))]] dx dy
   !   } dz
   ! Applying the two second derivatives, the internal 2D dx dy integral has
   !   the form:
   ! Ix'y' = S [exp(-zeta1*((Px-Rx1)^2 + (Py-Ry1)^2)) coeff *
   !    zeta2^2 4 *
   !    zeta2^3 -8 * (Py-Ry2)^2 *
   !    zeta2^3 -8 * (Px-Rx2)^2 *
   !    zeta2^4 16 * (Px-Rx2)^2*(Py-Ry2)^2 *
   !    exp(-zeta2*((Px-Rx2)^2 + (Py-Ry2)^2))]
   !
   ! The total integral is the combination of the fourth derivative terms and
   !   the product of second derivative terms.

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

   ! Initialize terms used to compute analytical factors.
   coeff = (fineStructure * 0.001d0)**2 / 8.0d0  ! 1/( 8 m^3 c^2)
   zeta = a1 + a2
   xi = a1 * a2 / zeta
   d = A - B
   S = (a1*A + a2*B) / zeta
   SB = S - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d*d))
   preFactor02 = (sum(SB(:)**2) + 3.0d0/(2.0d0*zeta))
   preFactor04 = sum(SB(:)**4) + sum(SB(:)**2)*3.0d0/zeta &
         & + 9.0d0/(4.0d0*zeta**2)
   preFactor22 = (SB(1)*SB(2))**2 + (SB(1)*SB(3))**2 + (SB(2)*SB(3))**2 + &
         & sum(SB(:)**2)/zeta + 3.0d0/(4.0d0*zeta**2)

   ! Initialize numerical mesh quantities.
   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   ! Initialize a counter of the triad pq pairs.
   h = 0

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

         ! Assign l1 and l2 values for each primitive Cartesian gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize sum solution.
         soln = 0.0d0

         ! Start a loop over 1D coordinates.
         do i = 0, num_steps

            ! Compute the current x position, distance, and sum.
            xyz(1) = (start_pos + (i*step_size))

            do j = 0, num_steps

               ! Compute the current y position, distance, and sum.
               xyz(2) = (start_pos + (j*step_size))

               do k = 0, num_steps

                  ! Compute the current z position, distance, and sum.
                  xyz(3) = (start_pos + (k*step_size))
                  sum_xyz_A_2 = sum((xyz(:)-A(:))**2)
                  sum_xyz_B_2 = sum((xyz(:)-B(:))**2)
                  product_xyz_A_l1 = product((xyz(:)-A(:))**l1(:))
                  product_xyz_B_l2 = product((xyz(:)-B(:))**l2(:))

!                  soln = soln + &
!&exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) * ( &
!&16*a2**4*preFactor04 - 80*a2**3*preFactor02 + 60*a2**2 + &
!&32*a2**4 * preFactor22)
! Above is for the <0|Q|0> case.

!&exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) * ( &
!&16*a2**4*sum_xyz_B_4 - 80*a2**3*sum_xyz_B_2 + 60*a2**2 + &
!&(32*a2**4*(xyz(1)-B(1))**2*(xyz(2)-B(2))**2) + &
!&(32*a2**4*(xyz(1)-B(1))**2*(xyz(3)-B(3))**2) + &
!&(32*a2**4*(xyz(2)-B(2))**2*(xyz(3)-B(3))**2))

!&exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) * ( &
!&(16*a2**4*(xyz(1)-B(1))**4-48*a2**3*(xyz(1)-B(1))**2) + &
!&(16*a2**4*(xyz(2)-B(2))**4-48*a2**3*(xyz(2)-B(2))**2) + &
!&(16*a2**4*(xyz(3)-B(3))**4-48*a2**3*(xyz(3)-B(3))**2) + 60*a2**2 + &
!&(32*a2**4*(xyz(1)-B(1))**2*(xyz(2)-B(2))**2 &
!&   -16*a2**3*(xyz(1)-B(1))**2 &
!&   -16*a2**3*(xyz(2)-B(2))**2) + &
!&(32*a2**4*(xyz(1)-B(1))**2*(xyz(3)-B(3))**2 &
!&   -16*a2**3*(xyz(1)-B(1))**2 &
!&   -16*a2**3*(xyz(3)-B(3))**2) + &
!&(32*a2**4*(xyz(2)-B(2))**2*(xyz(3)-B(3))**2 &
!&   -16*a2**3*(xyz(2)-B(2))**2 &
!&   -16*a2**3*(xyz(3)-B(3))**2))

!& exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) * ( &
!&(a2**4*(2*xyz(1)-2*B(1))**4-8*a2**3*(2*xyz(1)-2*B(1))**2-a2**3*(2*xyz(1)-2*B&
!&(1))*(8*xyz(1)-8*B(1))+12*a2**2) + &
!&(a2**4*(2*xyz(2)-2*B(2))**4-8*a2**3*(2*xyz(2)-2*B(2))**2-a2**3*(2*xyz(2)-2*B&
!&(2))*(8*xyz(2)-8*B(2))+12*a2**2) + &
!&(a2**4*(2*xyz(3)-2*B(3))**4-8*a2**3*(2*xyz(3)-2*B(3))**2-a2**3*(2*xyz(3)-2*B&
!&(3))*(8*xyz(3)-8*B(3))+12*a2**2) + &
!&2*(a2**4*(2*xyz(1)-2*B(1))**2*(2*xyz(2)-2*B(2))**2-2*a2**3*(2*xyz(1)-2*B(1))&
!&**2-2*a2**3*(2*xyz(2)-2*B(2))**2+4*a2**2) + &
!&2*(a2**4*(2*xyz(1)-2*B(1))**2*(2*xyz(3)-2*B(3))**2-2*a2**3*(2*xyz(1)-2*B(1))&
!&**2-2*a2**3*(2*xyz(3)-2*B(3))**2+4*a2**2) + &
!&2*(a2**4*(2*xyz(2)-2*B(2))**2*(2*xyz(3)-2*B(3))**2-2*a2**3*(2*xyz(2)-2*B(2))&
!&**2-2*a2**3*(2*xyz(3)-2*B(3))**2+4*a2**2))


!&(a2**4*(2*xyz(1)-2*B(1))**4-8*a2**3*(2*xyz(1)-2*B(1))**2-a2**3*(2*xyz(1)-2*B(1&
!&))*(8*xyz(1)-8*B(1))+12*a2**2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)+(a2**&
!&4*(2*xyz(2)-2*B(2))**4-8*a2**3*(2*xyz(2)-2*B(2))**2-a2**3*(2*xyz(2)-2*B(2))*(8&
!&*xyz(2)-8*B(2))+12*a2**2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)+(a2**4*(2*&
!&xyz(3)-2*B(3))**4-8*a2**3*(2*xyz(3)-2*B(3))**2-a2**3*(2*xyz(3)-2*B(3))*(8*xyz(&
!&3)-8*B(3))+12*a2**2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)+2*(a2**4*(2*xyz&
!&(1)-2*B(1))**2*(2*xyz(2)-2*B(2))**2-2*a2**3*(2*xyz(1)-2*B(1))**2-2*a2**3*(2*xy&
!&z(2)-2*B(2))**2+4*a2**2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)+2*(a2**4*(2&
!&*xyz(1)-2*B(1))**2*(2*xyz(3)-2*B(3))**2-2*a2**3*(2*xyz(1)-2*B(1))**2-2*a2**3*(&
!&2*xyz(3)-2*B(3))**2+4*a2**2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)+2*(a2**&
!&4*(2*xyz(2)-2*B(2))**2*(2*xyz(3)-2*B(3))**2-2*a2**3*(2*xyz(2)-2*B(2))**2-2*a2*&
!&*3*(2*xyz(3)-2*B(3))**2+4*a2**2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)

! All above is for the <0|Q|0> case.


                  soln = soln + &
&product_xyz_A_l1*(&
&l2(1)**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**4 &
&-6*l2(1)**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**4 &
&+11*l2(1)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**4 &
&-6*l2(1)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**4 &
&+6*l2(1)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(1)-B(1))**2 &
&-12*l2(1)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
&-6*l2(1)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(1)-B(1))**2 &
&+12*l2(1)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
&+4*l2(1)*a2**2*(xyz(1)-B(1))**l2(1)*(8*xyz(1)-8*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1)) &
&-4*l2(1)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**3*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1)) &
&+12*l2(1)**2*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
&-4*l2(1)**3*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
&+8*l2(1)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
&-8*l2(1)*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
&+a2**4*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**4*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
&-8*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2 &
&-a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(8*xyz(1) &
&-8*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
&+12*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
&+product_xyz_A_l1*(&
&l2(2)**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**4 &
&-6*l2(2)**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**4 &
&+12*l2(2)**2*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3+11*l2(2)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**4 &
&-6*l2(2)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**4 &
&+6*l2(2)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(2)-B(2))**2 &
&-12*l2(2)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
&-6*l2(2)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(2)-B(2))**2 &
&+12*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
&-4*l2(2)*a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**3*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2)) &
&+4*l2(2)*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(8*xyz(2)-8*B(2))*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2)) &
&+8*l2(2)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
&-4*l2(2)**3*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
&-8*l2(2)*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
&+a2**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**4*(xyz(3)-B(3))**l2(3) &
&-8*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2 &
&-a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(8*xyz(2)-8*B(2))*(xyz(3)-B(3))**l2(3) &
&+12*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
&+product_xyz_A_l1*(&
&l2(3)**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**4 &
&-6*l2(3)**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**4 &
&+11*l2(3)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**4 &
&-6*l2(3)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**4 &
&+6*l2(3)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(3)-B(3))**2 &
&-12*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2&
&-6*l2(3)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(3)-B(3))**2 &
&+12*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2 &
&+4*l2(3)*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(8*xyz(3)-8*B(3))/(xyz(3)-B(3)) &
&-4*l2(3)*a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))**3/(xyz(3)-B(3)) &
&-8*l2(3)*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
&-4*l2(3)**3*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
&+12*l2(3)**2*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
&+8*l2(3)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
&+a2**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))**4 &
&-8*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2 &
&-a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))*(8*xyz(3)-8*B(3)) &
&+12*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
&+2*product_xyz_A_l1*(l2(1)**2*l2(2)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
&-2*l2(1)**2*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))) &
&-l2(1)**2*l2(2)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
&+l2(1)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(1)-B(1))**2 &
&-2*l2(1)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
&-2*l2(1)*l2(2)**2*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(2)-B(2))**2) &
&-l2(1)*l2(2)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
&+4*l2(1)*l2(2)*a2**2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(2)-B(2))) &
&+2*l2(1)*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(2)-B(2))**2) &
&+2*l2(1)*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))) &
&+l2(1)*l2(2)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
&-2*l2(1)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**2*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1)) &
&+4*l2(1)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)-l2(1)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(1)-B(1))**2 &
&+2*l2(1)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
&+l2(2)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(2)-B(2))**2 &
&-2*l2(2)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
&-2*l2(2)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**2*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2)) &
&-l2(2)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(2)-B(2))**2 &
&+4*l2(2)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
&+2*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
&+a2**4*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**2*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**2*(xyz(3)-B(3))**l2(3) &
&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2 &
&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2 &
&+4*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
&+2*product_xyz_A_l1*(l2(1)**2*l2(3)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
&-2*l2(1)**2*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(1)-B(1))**2*(xyz(3)-B(3))) &
&-l2(1)**2*l2(3)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
&+l2(1)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(1)-B(1))**2 &
&-2*l2(1)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
&-2*l2(1)*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(3)-B(3))**2) &
&-l2(1)*l2(3)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
&+4*l2(1)*l2(3)*a2**2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(1)-B(1))*(xyz(3)-B(3))) &
&+2*l2(1)*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(3)-B(3))**2) &
&+2*l2(1)*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(1)-B(1))**2*(xyz(3)-B(3))) &
&+l2(1)*l2(3)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
&-2*l2(1)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))**2/(xyz(1)-B(1)) &
&+4*l2(1)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
&-l2(1)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(1)-B(1))**2 &
&+2*l2(1)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
&+l2(3)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(3)-B(3))**2 &
&-2*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2 &
&-2*l2(3)*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2*(2*xyz(3)-2*B(3))/(xyz(3)-B(3)) &
&-l2(3)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(3)-B(3))**2 &
&+4*l2(3)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
&+2*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2 &
&+a2**4*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2*(2*xyz(3)-2*B(3))**2 &
&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2 &
&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2 &
&+4*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
&+2*product_xyz_A_l1*(l2(2)**2*l2(3)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
&-2*l2(2)**2*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(2)-B(2))**2*(xyz(3)-B(3))) &
&-l2(2)**2*l2(3)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
&+l2(2)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(2)-B(2))**2 &
&-2*l2(2)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
&-2*l2(2)*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))*(xyz(3)-B(3))**2) &
&-l2(2)*l2(3)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
&+4*l2(2)*l2(3)*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(2)-B(2))*(xyz(3)-B(3))) &
&+2*l2(2)*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))*(xyz(3)-B(3))**2) &
&+2*l2(2)*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(2)-B(2))**2*(xyz(3)-B(3))) &
&+l2(2)*l2(3)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
&-2*l2(2)*a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))**2/(xyz(2)-B(2)) &
&+4*l2(2)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
&-l2(2)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(2)-B(2))**2 &
&+2*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
&+l2(3)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(3)-B(3))**2 &
&-2*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2 &
&-2*l2(3)*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2*(2*xyz(3)-2*B(3))/(xyz(3)-B(3)) &
&-l2(3)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(3)-B(3))**2 &
&+4*l2(3)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
&+2*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2 &
&+a2**4*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2*(2*xyz(3)-2*B(3))**2 &
&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2 &
&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2 &
&+4*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)
! Same as below, except we have grouped terms by the presence of a division.



!&product_xyz_A_l1*(l2(1)**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**4 &
!&-4*l2(1)**3*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
!&-6*l2(1)**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**4 &
!&+6*l2(1)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(1)-B(1))**2 &
!&-12*l2(1)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
!&+12*l2(1)**2*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3+11*l2(1)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**4 &
!&-4*l2(1)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**3*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1)) &
!&+8*l2(1)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&+4*l2(1)*a2**2*(xyz(1)-B(1))**l2(1)*(8*xyz(1)-8*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1)) &
!&-6*l2(1)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(1)-B(1))**2 &
!&+12*l2(1)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
!&-8*l2(1)*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
!&-6*l2(1)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**4 &
!&+a2**4*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**4*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&-8*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2 &
!&-a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(8*xyz(1) &
!&-8*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&+12*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
!&+product_xyz_A_l1*(l2(2)**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**4 &
!&-4*l2(2)**3*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
!&-6*l2(2)**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**4 &
!&+6*l2(2)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(2)-B(2))**2 &
!&-12*l2(2)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
!&+12*l2(2)**2*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3+11*l2(2)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**4 &
!&-4*l2(2)*a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**3*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2)) &
!&+8*l2(2)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&+4*l2(2)*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(8*xyz(2)-8*B(2))*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2)) &
!&-6*l2(2)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(2)-B(2))**2 &
!&+12*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
!&-8*l2(2)*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
!&-6*l2(2)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**4 &
!&+a2**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**4*(xyz(3)-B(3))**l2(3) &
!&-8*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2 &
!&-a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(8*xyz(2)-8*B(2))*(xyz(3)-B(3))**l2(3) &
!&+12*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
!&+product_xyz_A_l1*(l2(3)**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**4 &
!&-4*l2(3)**3*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
!&-6*l2(3)**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**4 &
!&+6*l2(3)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(3)-B(3))**2 &
!&-12*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2&
!&+12*l2(3)**2*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
!&+11*l2(3)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**4 &
!&-4*l2(3)*a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))**3/(xyz(3)-B(3)) &
!&+8*l2(3)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&+4*l2(3)*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(8*xyz(3)-8*B(3))/(xyz(3)-B(3)) &
!&-6*l2(3)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(3)-B(3))**2 &
!&+12*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2 &
!&-8*l2(3)*a2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)**3 &
!&-6*l2(3)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**4 &
!&+a2**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))**4 &
!&-8*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2 &
!&-a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))*(8*xyz(3)-8*B(3)) &
!&+12*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
!&+2*product_xyz_A_l1*(l2(1)**2*l2(2)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
!&-2*l2(1)**2*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))) &
!&-l2(1)**2*l2(2)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
!&+l2(1)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(1)-B(1))**2 &
!&-2*l2(1)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
!&-2*l2(1)*l2(2)**2*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(2)-B(2))**2) &
!&-l2(1)*l2(2)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
!&+4*l2(1)*l2(2)*a2**2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(2)-B(2))) &
!&+2*l2(1)*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(2)-B(2))**2) &
!&+2*l2(1)*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))) &
!&+l2(1)*l2(2)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
!&-2*l2(1)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**2*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1)) &
!&+4*l2(1)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)-l2(1)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(1)-B(1))**2 &
!&+2*l2(1)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
!&+l2(2)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(2)-B(2))**2 &
!&-2*l2(2)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
!&-2*l2(2)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**2*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2)) &
!&-l2(2)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(2)-B(2))**2 &
!&+4*l2(2)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&+2*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
!&+a2**4*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**2*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**2*(xyz(3)-B(3))**l2(3) &
!&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2 &
!&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2 &
!&+4*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
!&+2*product_xyz_A_l1*(l2(1)**2*l2(3)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
!&-2*l2(1)**2*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(1)-B(1))**2*(xyz(3)-B(3))) &
!&-l2(1)**2*l2(3)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
!&+l2(1)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(1)-B(1))**2 &
!&-2*l2(1)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
!&-2*l2(1)*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(3)-B(3))**2) &
!&-l2(1)*l2(3)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
!&+4*l2(1)*l2(3)*a2**2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(1)-B(1))*(xyz(3)-B(3))) &
!&+2*l2(1)*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(3)-B(3))**2) &
!&+2*l2(1)*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(1)-B(1))**2*(xyz(3)-B(3))) &
!&+l2(1)*l2(3)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
!&-2*l2(1)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))**2/(xyz(1)-B(1)) &
!&+4*l2(1)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&-l2(1)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(1)-B(1))**2 &
!&+2*l2(1)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))**2 &
!&+l2(3)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(3)-B(3))**2 &
!&-2*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2 &
!&-2*l2(3)*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2*(2*xyz(3)-2*B(3))/(xyz(3)-B(3)) &
!&-l2(3)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2/(xyz(3)-B(3))**2 &
!&+4*l2(3)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&+2*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2 &
!&+a2**4*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2*(2*xyz(3)-2*B(3))**2 &
!&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(1)-B(1))**2 &
!&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2 &
!&+4*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
!&+2*product_xyz_A_l1*(l2(2)**2*l2(3)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
!&-2*l2(2)**2*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(2)-B(2))**2*(xyz(3)-B(3))) &
!&-l2(2)**2*l2(3)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
!&+l2(2)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(2)-B(2))**2 &
!&-2*l2(2)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
!&-2*l2(2)*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))*(xyz(3)-B(3))**2) &
!&-l2(2)*l2(3)**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
!&+4*l2(2)*l2(3)*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(2)-B(2))*(xyz(3)-B(3))) &
!&+2*l2(2)*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))*(xyz(3)-B(3))**2) &
!&+2*l2(2)*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(2)-B(2))**2*(xyz(3)-B(3))) &
!&+l2(2)*l2(3)*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
!&-2*l2(2)*a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))**2/(xyz(2)-B(2)) &
!&+4*l2(2)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&-l2(2)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2/(xyz(2)-B(2))**2 &
!&+2*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))**2 &
!&+l2(3)**2*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(3)-B(3))**2 &
!&-2*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2 &
!&-2*l2(3)*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2*(2*xyz(3)-2*B(3))/(xyz(3)-B(3)) &
!&-l2(3)*a2**2*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2/(xyz(3)-B(3))**2 &
!&+4*l2(3)*a2**2*2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&+2*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(3)-B(3))**2 &
!&+a2**4*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2*(2*xyz(3)-2*B(3))**2 &
!&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(2)-B(2))**2 &
!&-2*a2**3*4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(xyz(3)-B(3))**2 &
!&+4*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3))*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)
! Same as below, except we have replaced all product_xyz_B_l2 terms with their expansion.


!&product_xyz_A_l1*(l2(1)**4*product_xyz_B_l2/(xyz(1)-B(1))**4 &
!&-4*l2(1)**3*a2*2*product_xyz_B_l2**3 &
!&-6*l2(1)**3*product_xyz_B_l2/(xyz(1)-B(1))**4 &
!&+6*l2(1)**2*a2**2*4*product_xyz_B_l2*(xyz(1)-B(1))**2/(xyz(1)-B(1))**2 &
!&-12*l2(1)**2*a2*product_xyz_B_l2/(xyz(1)-B(1))**2 &
!&+12*l2(1)**2*a2*2*product_xyz_B_l2**3+11*l2(1)**2*product_xyz_B_l2/(xyz(1)-B(1))**4 &
!&-4*l2(1)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**3*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1)) &
!&+8*l2(1)*a2**2*2*product_xyz_B_l2 &
!&+4*l2(1)*a2**2*(xyz(1)-B(1))**l2(1)*(8*xyz(1)-8*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1)) &
!&-6*l2(1)*a2**2*4*product_xyz_B_l2*(xyz(1)-B(1))**2/(xyz(1)-B(1))**2 &
!&+12*l2(1)*a2*product_xyz_B_l2/(xyz(1)-B(1))**2 &
!&-8*l2(1)*a2*2*product_xyz_B_l2**3 &
!&-6*l2(1)*product_xyz_B_l2/(xyz(1)-B(1))**4 &
!&+a2**4*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**4*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&-8*a2**3*4*product_xyz_B_l2*(xyz(1)-B(1))**2 &
!&-a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(8*xyz(1) &
!&-8*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3) &
!&+12*a2**2*product_xyz_B_l2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
!&+product_xyz_A_l1*(l2(2)**4*product_xyz_B_l2/(xyz(2)-B(2))**4 &
!&-4*l2(2)**3*a2*2*product_xyz_B_l2**3 &
!&-6*l2(2)**3*product_xyz_B_l2/(xyz(2)-B(2))**4 &
!&+6*l2(2)**2*a2**2*4*product_xyz_B_l2*(xyz(2)-B(2))**2/(xyz(2)-B(2))**2 &
!&-12*l2(2)**2*a2*product_xyz_B_l2/(xyz(2)-B(2))**2 &
!&+12*l2(2)**2*a2*2*product_xyz_B_l2**3+11*l2(2)**2*product_xyz_B_l2/(xyz(2)-B(2))**4 &
!&-4*l2(2)*a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**3*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2)) &
!&+8*l2(2)*a2**2*2*product_xyz_B_l2 &
!&+4*l2(2)*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(8*xyz(2)-8*B(2))*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2)) &
!&-6*l2(2)*a2**2*4*product_xyz_B_l2*(xyz(2)-B(2))**2/(xyz(2)-B(2))**2 &
!&+12*l2(2)*a2*product_xyz_B_l2/(xyz(2)-B(2))**2 &
!&-8*l2(2)*a2*2*product_xyz_B_l2**3 &
!&-6*l2(2)*product_xyz_B_l2/(xyz(2)-B(2))**4 &
!&+a2**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**4*(xyz(3)-B(3))**l2(3) &
!&-8*a2**3*4*product_xyz_B_l2*(xyz(2)-B(2))**2 &
!&-a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(8*xyz(2)-8*B(2))*(xyz(3)-B(3))**l2(3) &
!&+12*a2**2*product_xyz_B_l2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
!&+product_xyz_A_l1*(l2(3)**4*product_xyz_B_l2/(xyz(3)-B(3))**4 &
!&-4*l2(3)**3*a2*2*product_xyz_B_l2**3 &
!&-6*l2(3)**3*product_xyz_B_l2/(xyz(3)-B(3))**4 &
!&+6*l2(3)**2*a2**2*4*product_xyz_B_l2*(xyz(3)-B(3))**2/(xyz(3)-B(3))**2 &
!&-12*l2(3)**2*a2*product_xyz_B_l2/(xyz(3)-B(3))**2&
!&+12*l2(3)**2*a2*2*product_xyz_B_l2**3 &
!&+11*l2(3)**2*product_xyz_B_l2/(xyz(3)-B(3))**4 &
!&-4*l2(3)*a2**3*product_xyz_B_l2*(2*xyz(3)-2*B(3))**3/(xyz(3)-B(3)) &
!&+8*l2(3)*a2**2*2*product_xyz_B_l2 &
!&+4*l2(3)*a2**2*product_xyz_B_l2*(8*xyz(3)-8*B(3))/(xyz(3)-B(3)) &
!&-6*l2(3)*a2**2*4*product_xyz_B_l2*(xyz(3)-B(3))**2/(xyz(3)-B(3))**2 &
!&+12*l2(3)*a2*product_xyz_B_l2/(xyz(3)-B(3))**2 &
!&-8*l2(3)*a2*2*product_xyz_B_l2**3 &
!&-6*l2(3)*product_xyz_B_l2/(xyz(3)-B(3))**4 &
!&+a2**4*product_xyz_B_l2*(2*xyz(3)-2*B(3))**4 &
!&-8*a2**3*4*product_xyz_B_l2*(xyz(3)-B(3))**2 &
!&-a2**3*product_xyz_B_l2*(2*xyz(3)-2*B(3))*(8*xyz(3)-8*B(3)) &
!&+12*a2**2*product_xyz_B_l2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
!&+2*product_xyz_A_l1*(l2(1)**2*l2(2)**2*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
!&-2*l2(1)**2*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))) &
!&-l2(1)**2*l2(2)*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
!&+l2(1)**2*a2**2*4*product_xyz_B_l2*(xyz(2)-B(2))**2/(xyz(1)-B(1))**2 &
!&-2*l2(1)**2*a2*product_xyz_B_l2/(xyz(1)-B(1))**2 &
!&-2*l2(1)*l2(2)**2*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(2)-B(2))**2) &
!&-l2(1)*l2(2)**2*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
!&+4*l2(1)*l2(2)*a2**2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(2)-B(2))) &
!&+2*l2(1)*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(2)-B(2))**2) &
!&+2*l2(1)*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2))) &
!&+l2(1)*l2(2)*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2) &
!&-2*l2(1)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**2*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1)) &
!&+4*l2(1)*a2**2*2*product_xyz_B_l2-l2(1)*a2**2*4*product_xyz_B_l2*(xyz(2)-B(2))**2/(xyz(1)-B(1))**2 &
!&+2*l2(1)*a2*product_xyz_B_l2/(xyz(1)-B(1))**2 &
!&+l2(2)**2*a2**2*4*product_xyz_B_l2*(xyz(1)-B(1))**2/(xyz(2)-B(2))**2 &
!&-2*l2(2)**2*a2*product_xyz_B_l2/(xyz(2)-B(2))**2 &
!&-2*l2(2)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**2*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2)) &
!&-l2(2)*a2**2*4*product_xyz_B_l2*(xyz(1)-B(1))**2/(xyz(2)-B(2))**2 &
!&+4*l2(2)*a2**2*2*product_xyz_B_l2 &
!&+2*l2(2)*a2*product_xyz_B_l2/(xyz(2)-B(2))**2 &
!&+a2**4*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**2*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**2*(xyz(3)-B(3))**l2(3) &
!&-2*a2**3*4*product_xyz_B_l2*(xyz(1)-B(1))**2 &
!&-2*a2**3*4*product_xyz_B_l2*(xyz(2)-B(2))**2 &
!&+4*a2**2*product_xyz_B_l2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
!&+2*product_xyz_A_l1*(l2(1)**2*l2(3)**2*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
!&-2*l2(1)**2*l2(3)*a2*product_xyz_B_l2*(2*xyz(3)-2*B(3))/((xyz(1)-B(1))**2*(xyz(3)-B(3))) &
!&-l2(1)**2*l2(3)*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
!&+l2(1)**2*a2**2*4*product_xyz_B_l2*(xyz(3)-B(3))**2/(xyz(1)-B(1))**2 &
!&-2*l2(1)**2*a2*product_xyz_B_l2/(xyz(1)-B(1))**2 &
!&-2*l2(1)*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(3)-B(3))**2) &
!&-l2(1)*l2(3)**2*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
!&+4*l2(1)*l2(3)*a2**2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(1)-B(1))*(xyz(3)-B(3))) &
!&+2*l2(1)*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(3)-B(3))**2) &
!&+2*l2(1)*l2(3)*a2*product_xyz_B_l2*(2*xyz(3)-2*B(3))/((xyz(1)-B(1))**2*(xyz(3)-B(3))) &
!&+l2(1)*l2(3)*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2) &
!&-2*l2(1)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))**2/(xyz(1)-B(1)) &
!&+4*l2(1)*a2**2*2*product_xyz_B_l2 &
!&-l2(1)*a2**2*4*product_xyz_B_l2*(xyz(3)-B(3))**2/(xyz(1)-B(1))**2 &
!&+2*l2(1)*a2*product_xyz_B_l2/(xyz(1)-B(1))**2 &
!&+l2(3)**2*a2**2*4*product_xyz_B_l2*(xyz(1)-B(1))**2/(xyz(3)-B(3))**2 &
!&-2*l2(3)**2*a2*product_xyz_B_l2/(xyz(3)-B(3))**2 &
!&-2*l2(3)*a2**3*4*product_xyz_B_l2*(xyz(1)-B(1))**2*(2*xyz(3)-2*B(3))/(xyz(3)-B(3)) &
!&-l2(3)*a2**2*4*product_xyz_B_l2*(xyz(1)-B(1))**2/(xyz(3)-B(3))**2 &
!&+4*l2(3)*a2**2*2*product_xyz_B_l2 &
!&+2*l2(3)*a2*product_xyz_B_l2/(xyz(3)-B(3))**2 &
!&+a2**4*4*product_xyz_B_l2*(xyz(1)-B(1))**2*(2*xyz(3)-2*B(3))**2 &
!&-2*a2**3*4*product_xyz_B_l2*(xyz(1)-B(1))**2 &
!&-2*a2**3*4*product_xyz_B_l2*(xyz(3)-B(3))**2 &
!&+4*a2**2*product_xyz_B_l2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2) &
!&+2*product_xyz_A_l1*(l2(2)**2*l2(3)**2*product_xyz_B_l2/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
!&-2*l2(2)**2*l2(3)*a2*product_xyz_B_l2*(2*xyz(3)-2*B(3))/((xyz(2)-B(2))**2*(xyz(3)-B(3))) &
!&-l2(2)**2*l2(3)*product_xyz_B_l2/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
!&+l2(2)**2*a2**2*4*product_xyz_B_l2*(xyz(3)-B(3))**2/(xyz(2)-B(2))**2 &
!&-2*l2(2)**2*a2*product_xyz_B_l2/(xyz(2)-B(2))**2 &
!&-2*l2(2)*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))*(xyz(3)-B(3))**2) &
!&-l2(2)*l2(3)**2*product_xyz_B_l2/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
!&+4*l2(2)*l2(3)*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(2)-B(2))*(xyz(3)-B(3))) &
!&+2*l2(2)*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))*(xyz(3)-B(3))**2) &
!&+2*l2(2)*l2(3)*a2*product_xyz_B_l2*(2*xyz(3)-2*B(3))/((xyz(2)-B(2))**2*(xyz(3)-B(3))) &
!&+l2(2)*l2(3)*product_xyz_B_l2/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2) &
!&-2*l2(2)*a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))**2/(xyz(2)-B(2)) &
!&+4*l2(2)*a2**2*2*product_xyz_B_l2-l2(2)*a2**2*4*product_xyz_B_l2*(xyz(3)-B(3))**2/(xyz(2)-B(2))**2 &
!&+2*l2(2)*a2*product_xyz_B_l2/(xyz(2)-B(2))**2 &
!&+l2(3)**2*a2**2*4*product_xyz_B_l2*(xyz(2)-B(2))**2/(xyz(3)-B(3))**2 &
!&-2*l2(3)**2*a2*product_xyz_B_l2/(xyz(3)-B(3))**2 &
!&-2*l2(3)*a2**3*4*product_xyz_B_l2*(xyz(2)-B(2))**2*(2*xyz(3)-2*B(3))/(xyz(3)-B(3)) &
!&-l2(3)*a2**2*4*product_xyz_B_l2*(xyz(2)-B(2))**2/(xyz(3)-B(3))**2 &
!&+4*l2(3)*a2**2*2*product_xyz_B_l2 &
!&+2*l2(3)*a2*product_xyz_B_l2/(xyz(3)-B(3))**2 &
!&+a2**4*4*product_xyz_B_l2*(xyz(2)-B(2))**2*(2*xyz(3)-2*B(3))**2 &
!&-2*a2**3*4*product_xyz_B_l2*(xyz(2)-B(2))**2 &
!&-2*a2**3*4*product_xyz_B_l2*(xyz(3)-B(3))**2 &
!&+4*a2**2*product_xyz_B_l2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)
! Same as below. Here, we've just isolate terms.


!&product_xyz_A_l1*(l2(1)**4*product_xyz_B_l2/(xyz(1)-B(1))**4-4*l2(1)**3*a2*2*p&
!&roduct_xyz_B_l2**3-6*l2(1)**3*product_xyz_B_l2/(xyz(1)-B(1))**4+6*l2(1)**2*a2*&
!&*2*4*product_xyz_B_l2*(xyz(1)-B(1))**2/(xyz(1)-B(1))**2-12*l2(1)**2*a2*product&
!&_xyz_B_l2/(xyz(1)-B(1))**2+12*l2(1)**2*a2*2*product_xyz_B_l2**3+11*l2(1)**2*pr&
!&oduct_xyz_B_l2/(xyz(1)-B(1))**4-4*l2(1)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2&
!&*B(1))**3*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))+8*l2(1)*a2**&
!&2*2*product_xyz_B_l2+4*l2(1)*a2**2*(xyz(1)-B(1))**l2(1)*(8*xyz(1)-8*B(1))*(xyz&
!&(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/(xyz(1)-B(1))-6*l2(1)*a2**2*4*product_xy&
!&z_B_l2*(xyz(1)-B(1))**2/(xyz(1)-B(1))**2+12*l2(1)*a2*product_xyz_B_l2/(xyz(1)-&
!&B(1))**2-8*l2(1)*a2*2*product_xyz_B_l2**3-6*l2(1)*product_xyz_B_l2/(xyz(1)-B(1&
!&))**4+a2**4*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**4*(xyz(2)-B(2))**l2(2)*(xy&
!&z(3)-B(3))**l2(3)-8*a2**3*4*product_xyz_B_l2*(xyz(1)-B(1))**2-a2**3*(xyz(1)-B(&
!&1))**l2(1)*(2*xyz(1)-2*B(1))*(8*xyz(1)-8*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(&
!&3))**l2(3)+12*a2**2*product_xyz_B_l2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2&
!&)+product_xyz_A_l1*(l2(2)**4*product_xyz_B_l2/(xyz(2)-B(2))**4-4*l2(2)**3*a2*2&
!&*product_xyz_B_l2**3-6*l2(2)**3*product_xyz_B_l2/(xyz(2)-B(2))**4+6*l2(2)**2*a&
!&2**2*4*product_xyz_B_l2*(xyz(2)-B(2))**2/(xyz(2)-B(2))**2-12*l2(2)**2*a2*produ&
!&ct_xyz_B_l2/(xyz(2)-B(2))**2+12*l2(2)**2*a2*2*product_xyz_B_l2**3+11*l2(2)**2*&
!&product_xyz_B_l2/(xyz(2)-B(2))**4-4*l2(2)*a2**3*(xyz(1)-B(1))**l2(1)*(xyz(2)-B&
!&(2))**l2(2)*(2*xyz(2)-2*B(2))**3*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))+8*l2(2)*a2&
!&**2*2*product_xyz_B_l2+4*l2(2)*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)&
!&*(8*xyz(2)-8*B(2))*(xyz(3)-B(3))**l2(3)/(xyz(2)-B(2))-6*l2(2)*a2**2*4*product_&
!&xyz_B_l2*(xyz(2)-B(2))**2/(xyz(2)-B(2))**2+12*l2(2)*a2*product_xyz_B_l2/(xyz(2&
!&)-B(2))**2-8*l2(2)*a2*2*product_xyz_B_l2**3-6*l2(2)*product_xyz_B_l2/(xyz(2)-B&
!&(2))**4+a2**4*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**4*(&
!&xyz(3)-B(3))**l2(3)-8*a2**3*4*product_xyz_B_l2*(xyz(2)-B(2))**2-a2**3*(xyz(1)-&
!&B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(8*xyz(2)-8*B(2))*(xyz(3)-&
!&B(3))**l2(3)+12*a2**2*product_xyz_B_l2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B&
!&_2)+product_xyz_A_l1*(l2(3)**4*product_xyz_B_l2/(xyz(3)-B(3))**4-4*l2(3)**3*a2&
!&*2*product_xyz_B_l2**3-6*l2(3)**3*product_xyz_B_l2/(xyz(3)-B(3))**4+6*l2(3)**2&
!&*a2**2*4*product_xyz_B_l2*(xyz(3)-B(3))**2/(xyz(3)-B(3))**2-12*l2(3)**2*a2*pro&
!&duct_xyz_B_l2/(xyz(3)-B(3))**2+12*l2(3)**2*a2*2*product_xyz_B_l2**3+11*l2(3)**&
!&2*product_xyz_B_l2/(xyz(3)-B(3))**4-4*l2(3)*a2**3*product_xyz_B_l2*(2*xyz(3)-2&
!&*B(3))**3/(xyz(3)-B(3))+8*l2(3)*a2**2*2*product_xyz_B_l2+4*l2(3)*a2**2*product&
!&_xyz_B_l2*(8*xyz(3)-8*B(3))/(xyz(3)-B(3))-6*l2(3)*a2**2*4*product_xyz_B_l2*(xy&
!&z(3)-B(3))**2/(xyz(3)-B(3))**2+12*l2(3)*a2*product_xyz_B_l2/(xyz(3)-B(3))**2-8&
!&*l2(3)*a2*2*product_xyz_B_l2**3-6*l2(3)*product_xyz_B_l2/(xyz(3)-B(3))**4+a2**&
!&4*product_xyz_B_l2*(2*xyz(3)-2*B(3))**4-8*a2**3*4*product_xyz_B_l2*(xyz(3)-B(3&
!&))**2-a2**3*product_xyz_B_l2*(2*xyz(3)-2*B(3))*(8*xyz(3)-8*B(3))+12*a2**2*prod&
!&uct_xyz_B_l2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)+2*product_xyz_A_l1*(l2&
!&(1)**2*l2(2)**2*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2)-2*l2(1)**&
!&2*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)&
!&-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2)))-l2(1)**2*l2(2)*product_xyz_B_l2&
!&/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2)+l2(1)**2*a2**2*4*product_xyz_B_l2*(xyz(2)&
!&-B(2))**2/(xyz(1)-B(1))**2-2*l2(1)**2*a2*product_xyz_B_l2/(xyz(1)-B(1))**2-2*l&
!&2(1)*l2(2)**2*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(&
!&xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(2)-B(2))**2)-l2(1)*l2(2)**2*product_xy&
!&z_B_l2/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2)+4*l2(1)*l2(2)*a2**2*(xyz(1)-B(1))**&
!&l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**&
!&l2(3)/((xyz(1)-B(1))*(xyz(2)-B(2)))+2*l2(1)*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(2*x&
!&yz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(2)&
!&-B(2))**2)+2*l2(1)*l2(2)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2&
!&)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))**2*(xyz(2)-B(2)))+l2(1)*l2(2)*pr&
!&oduct_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(2)-B(2))**2)-2*l2(1)*a2**3*(xyz(1)-B(1))&
!&**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))**2*(xyz(3)-B(&
!&3))**l2(3)/(xyz(1)-B(1))+4*l2(1)*a2**2*2*product_xyz_B_l2-l2(1)*a2**2*4*produc&
!&t_xyz_B_l2*(xyz(2)-B(2))**2/(xyz(1)-B(1))**2+2*l2(1)*a2*product_xyz_B_l2/(xyz(&
!&1)-B(1))**2+l2(2)**2*a2**2*4*product_xyz_B_l2*(xyz(1)-B(1))**2/(xyz(2)-B(2))**&
!&2-2*l2(2)**2*a2*product_xyz_B_l2/(xyz(2)-B(2))**2-2*l2(2)*a2**3*(xyz(1)-B(1))*&
!&*l2(1)*(2*xyz(1)-2*B(1))**2*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3&
!&))**l2(3)/(xyz(2)-B(2))-l2(2)*a2**2*4*product_xyz_B_l2*(xyz(1)-B(1))**2/(xyz(2&
!&)-B(2))**2+4*l2(2)*a2**2*2*product_xyz_B_l2+2*l2(2)*a2*product_xyz_B_l2/(xyz(2&
!&)-B(2))**2+a2**4*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))**2*(xyz(2)-B(2))**l2(2&
!&)*(2*xyz(2)-2*B(2))**2*(xyz(3)-B(3))**l2(3)-2*a2**3*4*product_xyz_B_l2*(xyz(1)&
!&-B(1))**2-2*a2**3*4*product_xyz_B_l2*(xyz(2)-B(2))**2+4*a2**2*product_xyz_B_l2&
!&)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)+2*product_xyz_A_l1*(l2(1)**2*l2(3)&
!&**2*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2)-2*l2(1)**2*l2(3)*a2*p&
!&roduct_xyz_B_l2*(2*xyz(3)-2*B(3))/((xyz(1)-B(1))**2*(xyz(3)-B(3)))-l2(1)**2*l2&
!&(3)*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2)+l2(1)**2*a2**2*4*prod&
!&uct_xyz_B_l2*(xyz(3)-B(3))**2/(xyz(1)-B(1))**2-2*l2(1)**2*a2*product_xyz_B_l2/&
!&(xyz(1)-B(1))**2-2*l2(1)*l2(3)**2*a2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(x&
!&yz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((xyz(1)-B(1))*(xyz(3)-B(3))**2)-l2(1)&
!&*l2(3)**2*product_xyz_B_l2/((xyz(1)-B(1))**2*(xyz(3)-B(3))**2)+4*l2(1)*l2(3)*a&
!&2**2*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))&
!&**l2(3)*(2*xyz(3)-2*B(3))/((xyz(1)-B(1))*(xyz(3)-B(3)))+2*l2(1)*l2(3)*a2*(xyz(&
!&1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)/((&
!&xyz(1)-B(1))*(xyz(3)-B(3))**2)+2*l2(1)*l2(3)*a2*product_xyz_B_l2*(2*xyz(3)-2*B&
!&(3))/((xyz(1)-B(1))**2*(xyz(3)-B(3)))+l2(1)*l2(3)*product_xyz_B_l2/((xyz(1)-B(&
!&1))**2*(xyz(3)-B(3))**2)-2*l2(1)*a2**3*(xyz(1)-B(1))**l2(1)*(2*xyz(1)-2*B(1))*&
!&(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))**2/(xyz(1)-B(1))+4&
!&*l2(1)*a2**2*2*product_xyz_B_l2-l2(1)*a2**2*4*product_xyz_B_l2*(xyz(3)-B(3))**&
!&2/(xyz(1)-B(1))**2+2*l2(1)*a2*product_xyz_B_l2/(xyz(1)-B(1))**2+l2(3)**2*a2**2&
!&*4*product_xyz_B_l2*(xyz(1)-B(1))**2/(xyz(3)-B(3))**2-2*l2(3)**2*a2*product_xy&
!&z_B_l2/(xyz(3)-B(3))**2-2*l2(3)*a2**3*4*product_xyz_B_l2*(xyz(1)-B(1))**2*(2*x&
!&yz(3)-2*B(3))/(xyz(3)-B(3))-l2(3)*a2**2*4*product_xyz_B_l2*(xyz(1)-B(1))**2/(x&
!&yz(3)-B(3))**2+4*l2(3)*a2**2*2*product_xyz_B_l2+2*l2(3)*a2*product_xyz_B_l2/(x&
!&yz(3)-B(3))**2+a2**4*4*product_xyz_B_l2*(xyz(1)-B(1))**2*(2*xyz(3)-2*B(3))**2-&
!&2*a2**3*4*product_xyz_B_l2*(xyz(1)-B(1))**2-2*a2**3*4*product_xyz_B_l2*(xyz(3)&
!&-B(3))**2+4*a2**2*product_xyz_B_l2)*exp(-a1*sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)+&
!&2*product_xyz_A_l1*(l2(2)**2*l2(3)**2*product_xyz_B_l2/((xyz(2)-B(2))**2*(xyz(&
!&3)-B(3))**2)-2*l2(2)**2*l2(3)*a2*product_xyz_B_l2*(2*xyz(3)-2*B(3))/((xyz(2)-B&
!&(2))**2*(xyz(3)-B(3)))-l2(2)**2*l2(3)*product_xyz_B_l2/((xyz(2)-B(2))**2*(xyz(&
!&3)-B(3))**2)+l2(2)**2*a2**2*4*product_xyz_B_l2*(xyz(3)-B(3))**2/(xyz(2)-B(2))*&
!&*2-2*l2(2)**2*a2*product_xyz_B_l2/(xyz(2)-B(2))**2-2*l2(2)*l2(3)**2*a2*(xyz(1)&
!&-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xy&
!&z(2)-B(2))*(xyz(3)-B(3))**2)-l2(2)*l2(3)**2*product_xyz_B_l2/((xyz(2)-B(2))**2&
!&*(xyz(3)-B(3))**2)+4*l2(2)*l2(3)*a2**2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(&
!&2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)*(2*xyz(3)-2*B(3))/((xyz(2)-B(2))*(xy&
!&z(3)-B(3)))+2*l2(2)*l2(3)*a2*(xyz(1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(&
!&2)-2*B(2))*(xyz(3)-B(3))**l2(3)/((xyz(2)-B(2))*(xyz(3)-B(3))**2)+2*l2(2)*l2(3)&
!&*a2*product_xyz_B_l2*(2*xyz(3)-2*B(3))/((xyz(2)-B(2))**2*(xyz(3)-B(3)))+l2(2)*&
!&l2(3)*product_xyz_B_l2/((xyz(2)-B(2))**2*(xyz(3)-B(3))**2)-2*l2(2)*a2**3*(xyz(&
!&1)-B(1))**l2(1)*(xyz(2)-B(2))**l2(2)*(2*xyz(2)-2*B(2))*(xyz(3)-B(3))**l2(3)*(2&
!&*xyz(3)-2*B(3))**2/(xyz(2)-B(2))+4*l2(2)*a2**2*2*product_xyz_B_l2-l2(2)*a2**2*&
!&4*product_xyz_B_l2*(xyz(3)-B(3))**2/(xyz(2)-B(2))**2+2*l2(2)*a2*product_xyz_B_&
!&l2/(xyz(2)-B(2))**2+l2(3)**2*a2**2*4*product_xyz_B_l2*(xyz(2)-B(2))**2/(xyz(3)&
!&-B(3))**2-2*l2(3)**2*a2*product_xyz_B_l2/(xyz(3)-B(3))**2-2*l2(3)*a2**3*4*prod&
!&uct_xyz_B_l2*(xyz(2)-B(2))**2*(2*xyz(3)-2*B(3))/(xyz(3)-B(3))-l2(3)*a2**2*4*pr&
!&oduct_xyz_B_l2*(xyz(2)-B(2))**2/(xyz(3)-B(3))**2+4*l2(3)*a2**2*2*product_xyz_B&
!&_l2+2*l2(3)*a2*product_xyz_B_l2/(xyz(3)-B(3))**2+a2**4*4*product_xyz_B_l2*(xyz&
!&(2)-B(2))**2*(2*xyz(3)-2*B(3))**2-2*a2**3*4*product_xyz_B_l2*(xyz(2)-B(2))**2-&
!&2*a2**3*4*product_xyz_B_l2*(xyz(3)-B(3))**2+4*a2**2*product_xyz_B_l2)*exp(-a1*&
!&sum_xyz_A_2)*exp(-a2*sum_xyz_B_2)
! Above is for the full case. It does not actually work because of division by
!   zero. I.e., the single block of code has not been refined to distinguish
!   different angular momentum cases.

               enddo
            enddo
         enddo

         ! Multiply the results by the mesh step size.
         soln = soln * step_size**3

         pc(q,p) = soln * coeff

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
            
"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [5, 1, ""], 0, num_conversion,
                       num_conversion)


    # Print the subroutine foot and the subsequent functions.
    foot = """
   end subroutine massvel2CIntgNum
"""
    f.write(foot)


def print_test_massvel_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
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
   real (kind=double), dimension (""" + \
           f"{len(triads)},{len(triads)}" + """), intent(out) :: pc
   real (kind=double), dimension (""" + \
           f"{len(conversion)},{len(conversion)}" + """), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, h, i, j, k
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
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
   coeff = (fineStructure * 0.001d0)**2 / 8.0d0
   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   ! Initialize a counter of the triad pq pairs.
   h = 0

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [5, 1, ""], 0, num_conversion,
                       num_conversion)


    # Print the subroutine foot and the subsequent functions.
    foot = """
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
"""
    f.write(foot)


def print_test_dkinetic_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
   subroutine dkinetic2CIntgNum(a1,a2,A,B,pc,sh,cell_size,step_size)

   use O_Kinds
   use O_Constants, only: pi, fineStructure

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
   real (kind=double), dimension (""" + \
           f"{len(triads)},{len(triads)}" + """,3), intent(out) :: pc
   real (kind=double), dimension (""" + \
           f"{len(conversion)},{len(conversion)}" + """,3), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, h, i, j, k
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, coeff
   real (kind=double), dimension (2) :: curr_pos
   real (kind=double), dimension (3) :: xyz, xyz_soln
   real (kind=double), dimension (6,3) :: xyz_I
      ! The second index is for prime, noprime, 2Dprime.
      ! The first index is xyz for prime and noprime while it is xy, xz, yx,
      !   yz, zx, zy for the 2Dprime case.

   ! Before we proceed with the calculation we need to understand a bit more
   !   about exactly what is being computed. The form of the integration is:
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        coeff * dellR2 ( dell^2
   !        [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))]) dxdydz
   ! We use Pxyz for points in space; Rxyz1,2 for atomic sites 1 and 2;
   !   zeta1,2 for the exponential decay factors; lxyz1,2 for the angular
   !   momentum numbers for atoms 1 and 2; and coeff = -1/2. Note that "dell"
   !   implies derivatives with respect to electron positional coordinates
   !   x,y,z while "dellR2" implies derivatives with respect to the atomic
   !   coordinates of site 2. Derivatives of a Cartesian Gaussian wrt the
   !   atomic site are the negatives of the derivatives wrt the electron
   !   coordinates. I.e., dG/dx = -dG/dBx so that dell = -dellR2.
   !
   ! The dell^2 operator is d^2/dx^2 + d^2/dy^2 + d^2/dz^2 so that dellR2 of
   !   dell^2 (or -dell of dell^2) is:
   !     -[d^3/dx^3 + d/dx d^2/dy^2 + d/dx d^2/dz^2] x-hat
   !     -[d/dy d^2/dx^2 + d^3/dy^3 + d/dy d^2/dz^2] y-hat
   !     -[d/dz d^2/dx^2 + d/dz d^2/dy^2 + d^3/dz^3] z-hat
   ! The integral must be computed over all space for all different possible
   !   values of lxyz1,2 for some arbitrarily chosen test coordinates and
   !   decay rates.
   !
   ! Because of the plus signs in the dell of dell^2 operator we arrive
   !   at three identical integrals of the form (with # in d^3/d#^3 = x, y,
   !   or z):
   ! -SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 *
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        coeff * d^3/d#^3 [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   !
   ! Also, because of the plus signs, we arrive at six other identical
   !   integrals of the form (with #,@ = x,y ; x,z ; y,x ; y,z ; z,x ; z,y):
   ! -SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        coeff * d/d# d^2/d@^2
   !        [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   !
   ! Now, focusing on just the d^3/dx^3 version of the set of three integrals,
   !   we will pull all terms with Py and Pz out of the dx integral to get:
   ! -SS { [(Py-Ry1)**ly1 * (Pz-Rz1)**lz1 *
   !        exp(-zeta1*((Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !       [(Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*((Py-Ry2)**2 + (Pz-Rz2)**2 ))] *
   !      S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * coeff * d^3/dx^3 [
   !         (Px-Rx2)**lx2 * exp(-zeta2*((Px-Rx2)**2))]] dx
   !     } dydz
   ! Applying the third derivative, the internal 1D dx integral has the form:
   !   Ix' = S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * coeff *
   !          [lx2*(lx2^2 - 3lx2 + 2) * (Px-Rx2)^(lx2-3) +
   !           (zeta2*-6lx2^2) * (Px-Rx2)^(lx2-1) +
   !           (zeta2^2*(12*lx2 + 12)) * (Px-Rx2)^(lx2+1) +
   !           (zeta2^3*-8) * (Px-Rx2)^(lx2+3)] * exp(-zeta2*(Px - Rx2)**2)
   !        ]
   ! Each of the other integrals (Iy, Iz) will have the form of a simple
   !   1D overlap integral:
   !   Iy = S [(Py-Ry1)^ly1 * exp(-zeta1*((Py-Ry1)^2)) *
   !            (Py-Ry2)^ly2 * exp(-zeta2*((Py-Ry2)^2))] dy
   !   Iz = S [(Pz-Rz1)^lz1 * exp(-zeta1*((Pz-Rz1)^2)) *
   !            (Pz-Rz2)^lz2 * exp(-zeta2*((Pz-Rz2)^2))] dz
   ! The total integral !! for the d^3/dx^3 portion !! is thus Ix' * Iy * Iz
   !   in the x-hat direcetion. Similarly, the d^3/dy^3 portion in the y-hat
   !   direction would be Ix * Iy' * Iz and the d^3/dz^3 portion in the z-hat
   !   direction would be Ix * Iy * Iz'. The Iy' and Iz' are the appropriate
   !   analongs of Ix' and Ix is the analog of the Iy or Iz.
   !
   ! Now, focusing on just the d/d# d^2/d@^2 where (say) #,@ = x,y,
   !   we will pull terms with Pz out of the dx dy integral to get:
   ! -S { [(Pz-Rz1)^lz1 * exp(-zeta1*(Pz-Rz1)^2)] *
   !      [(Pz-Rz2)^lz2 * exp(-zeta2*(Pz-Rz2)^2)] *
   !     SS [(Px-Rx1)^lx1 * (Py-Ry1)^ly1 * 
   !        exp(-zeta1*((Px-Rx1)^2 + (Py-Ry1)^2)) * coeff *
   !        d/dx d^2/dy^2 [(Px-Rx2)^lx2 * (Py-Ry2)^ly2 *
   !        exp(-zeta2*((Px-Rx2)^2 + (Py-Ry2)**2))]] dx dy
   !    } dz
   ! Applying the two second derivatives, the internal 2D dx dy integral has
   !   the form:
   ! Ix'y" = S [(Px-Rx1)^lx1 * (Py-Ry1)^ly1 *
   !    exp(-zeta1*((Px-Rx1)^2 + (Py-Ry1)^2)) * coeff *
   !    lx2*ly2 (ly2 - 1) (Px-Rx2)^(lx2-1)*(Py-Ry2)^(ly2-2) *
   !    zeta2 2lx2 (-2ly2 - 1) (Px-Rx2)^(lx2-1)*(Py-Ry2)^(ly2) *
   !    zeta2^2 4lx2 (Px-Rx2)^(lx2-1)*(Py-Ry2)^(ly2+2) *
   !    zeta2 2ly2 (-ly2 + 1) (Px-Rx2)^(lx2+1)*(Py-Ry2)^(ly2-2) *
   !    zeta2^2 4 (2ly2 + 1) (Px-Rx2)^(lx2+1)*(Py-Ry2)^(ly2) *
   !    zeta2^3 (-8) (Px-Rx2)^(lx2+1)*(Py-Ry2)^(ly2+2) *
   !    exp(-zeta2*((Px-Rx2)^2 + (Py-Ry2)^2))] dx dy
   !
   ! Thus, every total integral in a particular Cartesian direction is a sum
   !   of various independent integrals. We can do 1D overlap integrals for
   !   x, y, z; 1D third derivative integrals for x, y, z; and 2D second
   !   derivative integrals for xy, xz, yx, yz, zx, zy pairs to construct the
   !   total integral.
   ! With regard to the term in each integral that has the **(lx2-2) form
   !   (or **(ly2-2) or **(lz2-2), or -1 or -3 form), some special care must
   !   be taken. This term represents an angular momentum and the integral will
   !   have the form of an overlap integral. Because we cannot have a negative
   !   angular momentum we must discard this term when lx2-2 < 0, etc.

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
   coeff = 0.5d0 ! Sign cancelled with minus from dellR2 = -dell
   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   ! Initialize a counter of the triad pq pairs.
   h = 0

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

         ! Assign l1 and l2 values for each primitive Cartesian gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize sum variables.
         xyz_I(:,:) = 0.0d0

         ! Start a loop over 1D coordinates.
         do i = 0, num_steps
            curr_pos(1) = start_pos + (i*step_size)

            ! Compute the prime integrals first. (Third derivative)
            do j = 1, 3
               xyz_I(j,1) = xyz_I(j,1) &
                     & + primeDK(step_size, curr_pos(1), A(j), &
                     & B(j), a1, a2, l1(j), l2(j))
            enddo

            ! Compute the no-prime integrals second. (No derivative)
            do j = 1, 3
               xyz_I(j,2) = xyz_I(j,2) &
                     & + noPrimeDK(step_size, curr_pos(1), A(j), &
                     & B(j), a1, a2, l1(j), l2(j))
            enddo

            ! Start a loop over a second set of 1D coordiantes.
            do j = 0, num_steps
                curr_pos(2) = start_pos + (j*step_size)

                ! Compute the two-prime integrals.

                ! Compute the xy term.
                xyz_I(1,3) = xyz_I(1,3) &
                      & + twoPrimeDK(step_size, curr_pos(1), curr_pos(2), &
                      & A(1), A(2), B(1), B(2), a1, a2, l1(1), l1(2), l2(1), &
                      & l2(2))

                ! Compute the xz term.
                xyz_I(2,3) = xyz_I(2,3) &
                      & + twoPrimeDK(step_size, curr_pos(1), curr_pos(2), &
                      & A(1), A(3), B(1), B(3), a1, a2, l1(1), l1(3), l2(1), &
                      & l2(3))

                ! Compute the yx term.
                xyz_I(3,3) = xyz_I(3,3) &
                      & + twoPrimeDK(step_size, curr_pos(1), curr_pos(2), &
                      & A(2), A(1), B(2), B(1), a1, a2, l1(2), l1(1), l2(2), &
                      & l2(1))

                ! Compute the yz term.
                xyz_I(4,3) = xyz_I(4,3) &
                      & + twoPrimeDK(step_size, curr_pos(1), curr_pos(2), &
                      & A(2), A(3), B(2), B(3), a1, a2, l1(2), l1(3), l2(2), &
                      & l2(3))

                ! Compute the zx term.
                xyz_I(5,3) = xyz_I(5,3) &
                      & + twoPrimeDK(step_size, curr_pos(1), curr_pos(2), &
                      & A(3), A(1), B(3), B(1), a1, a2, l1(3), l1(1), l2(3), &
                      & l2(1))

                ! Compute the zy term.
                xyz_I(6,3) = xyz_I(6,3) &
                      & + twoPrimeDK(step_size, curr_pos(1), curr_pos(2), &
                      & A(3), A(2), B(3), B(2), a1, a2, l1(3), l1(2), l2(3), &
                      & l2(2))
            enddo
         enddo

         ! Assemble the integrals associated with the third derivative.
         pc(q,p,1) = xyz_I(1,1)*xyz_I(2,2)*xyz_I(3,2)  ! Third D wrt x only
         pc(q,p,2) = xyz_I(1,2)*xyz_I(2,1)*xyz_I(3,2)  ! Third D wrt y only
         pc(q,p,3) = xyz_I(1,2)*xyz_I(2,2)*xyz_I(3,1)  ! Third D wrt z only

         ! Add on the integrals associated with the second derivatives.

         ! The x-axis qp term uses d/dx d^2/dy^2 with no D wrt z plus
         !   d/dx d^2/dz^2 with no D wrt y.
         pc(q,p,1) = pc(q,p,1) + (xyz_I(1,3)*xyz_I(3,2) &
               & + (xyz_I(2,3)*xyz_I(2,2)))

         ! The y-axis qp term uses d/dy d^2/dx^2 with no D wrt z plus
         !   d/dy d^2/dz^2 with no D wrt x.
         pc(q,p,2) = pc(q,p,2) + (xyz_I(3,3)*xyz_I(3,2) &
               & + (xyz_I(4,3)*xyz_I(1,2)))

         ! The z-axis qp term uses d/dz d^2/dx^2 with no D wrt y plus
         !   d/dz d^2/dy^2 with no D wrt x.
         pc(q,p,3) = pc(q,p,3) + (xyz_I(5,3)*xyz_I(2,2) &
               & + (xyz_I(6,3)*xyz_I(1,2)))

         ! Multiply by the coefficient.
         pc(q,p,:) = pc(q,p,:) * coeff

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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [8, 3, ""], 0, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [8, 3, ""], 1, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [6, 3, ""], 2, num_conversion,
                       num_conversion)


    # Print the subroutine foot and the subsequent functions.
    foot = """
   end subroutine dkinetic2CIntgNum


   function noPrimeDK(step_size, curr_pos, A, B, a1, a2, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos, A, B, a1, a2
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: noPrimeDK

      ! Compute the integral part.
      noPrimeDK = step_size * (curr_pos - A)**l1 * (curr_pos - B)**l2 &
            & * exp(-a1*(curr_pos - A)**2) * exp(-a2*(curr_pos - B)**2)

      return

   end function noPrimeDK


   function primeDK(step_size, curr_pos, A, B, a1, a2, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos, A, B, a1, a2
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: primeDK

      ! Compute each "internal" term of the prime integral. Compare each
      !   of these lines with the 1D mass velocity equation produced by the
      !   osrecurintg_makenum.py script (appropriately separated into terms).
      !   Conceptually, these are similar to the expression in the last
      !   equals of equation 65 in Ben Walker's dissertation for the regular
      !   kinetic energy. Note that a factor of XXX is included in all of the
      !   below expressions. The most succinctly refactored equations are in
      !   ...

      primeDK = a2**3 * -8 * (curr_pos - B)**(l2+3)

      primeDK = primeDK + a2**2 * (12*l2 + 12) * (curr_pos - B)**(l2+1)

      if (l2 >= 1) then
         primeDK = primeDK + a2 * (-6*l2**2) * &
               & (curr_pos - B)**(l2-1)
      endif

      if (l2 >= 3) then
         primeDK = primeDK + l2*(l2**2 - 3*l2 + 2) * &
               & (curr_pos - B)**(l2-3)
      endif

      ! Multiply the prime integral by the preceeding primitive gaussian
      !   coefficient and exponential and multiply by the succeeding
      !   exponential. (We have already multiplied by the succeeding
      !   primitive gaussian coefficient in the above lines.)
      primeDK = primeDK * (curr_pos-A)**l1 * exp(-a1*(curr_pos-A)**2) &
            & * exp(-a2*(curr_pos-B)**2)

      ! Finally, multiply by the step size.
      primeDK = primeDK * step_size

      return

   end function primeDK


   ! The notation may be confusing so here is some clarification:
   ! As in the previous subroutines, the 1 and 2 for the alphas (a's)
   !   correspond to the first and second orbitals of the integral. I.e.,
   !   they are independent of x,y,z.
   ! The 1 and 2 for A corresponds to the position of site A with respect to
   !   the first xyz coordinate axis and second xyz coordinate axis
   !   respectively. So, when we are doing the xy term, then 1=x and 2=y.
   !   When we are doing the xz term, then 1=x and 2=z. Then, for the yz
   !   term 1=y and 2=z. Etc. for the yx, zx, and zy terms.
   ! The same concept described above for A applies to B.
   ! For the l variables the first number (1 or 2) corresponds to the first
   !   and second orbitals of the integrals. The second number (1 or 2)
   !   corresponds to the first xyz or second xyz coordinate axes as for the
   !   A1 A2 and B1 B2.

   function twoPrimeDK(step_size, curr_pos1, curr_pos2, A1, A2, B1, B2, &
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
      real (kind=double) :: twoPrimeDK

      ! Compute each internal term of the twoPrime integral. As with the prime
      !   integral, these terms are obtained from sympy derivations and simple
      !   algebraic manipulations.

      twoPrimeDK = alpha2**3*(-8) &
            & * (curr_pos1 - B1)**(l21 + 1) * (curr_pos2 - B2)**(l22 + 2)

      twoPrimeDK = twoPrimeDK + alpha2**2*4*(2*l22 + 1) &
            & * (curr_pos1 - B1)**(l21 + 1) * (curr_pos2 - B2)**(l22)

      if (l21 >= 1) then
         twoPrimeDK = twoPrimeDK + alpha2**2*4*l21 &
               & * (curr_pos1 - B1)**(l21-1) * (curr_pos2 - B2)**(l22+2)

         twoPrimeDK = twoPrimeDK + alpha2*2*l21*(-2*l22 - 1) &
               & * (curr_pos1 - B1)**(l21-1) * (curr_pos2 - B2)**(l22)
      endif

      if (l22 >= 2) then
         twoPrimeDK = twoPrimeDK + alpha2*2*l22*(-l22 + 1) &
               & * (curr_pos1 - B1)**(l21+1) * (curr_pos2 - B2)**(l22-2)
      endif

      if ((l21 >= 1) .and. (l22 >= 2)) then
         twoPrimeDK = twoPrimeDK + l21*l22*(l22 - 1) &
               & * (curr_pos1 - B1)**(l21 - 1) * (curr_pos2 - B2)**(l22 - 2)
      endif

      ! Multiply the integral by the preceeding primitive gaussian
      !   coefficient and exponential and multiply by the succeeding
      !   exponential. (We have already multiplied by the succeeding
      !   primitive gaussian coefficient in the above lines.)
      twoPrimeDK = twoPrimeDK * (curr_pos1-A1)**l11 * (curr_pos2-A2)**l12 * &
            & exp(-alpha1*((curr_pos1-A1)**2 + (curr_pos2-A2)**2)) * &
            & exp(-alpha2*((curr_pos1-B1)**2 + (curr_pos2-B2)**2))

      ! Finally, multiply by the step size squared (2D integral).
      twoPrimeDK = twoPrimeDK * step_size * step_size

      return

   end function twoPrimeDK
"""
    f.write(foot)


# The "full" test subroutine is intended to only be used for solving
#   integrals of the form (0|d/dB del|0) to guide the creation of the
#   derivative of the kinetic energy s-s case.
# Note that the numerical solution here is "full" because it does a full
#   3D numerical integral instead of the usual 1D and 2D integrals.
def print_test_dkinetic_full_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
   subroutine dkinetic2CIntgNum(a1,a2,A,B,pc,sh,cell_size,step_size)

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
   real (kind=double), dimension (""" + \
           f"{len(triads)},{len(triads)}" + """,3), intent(out) :: pc
   real (kind=double), dimension (""" + \
           f"{len(conversion)},{len(conversion)}" + """,3), intent(out) :: sh
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, h, i, j, k
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, coeff
   real (kind=double), dimension (2) :: curr_pos
   real (kind=double), dimension (3) :: xyz, xyz_soln
   real (kind=double), dimension (3,3) :: xyz_I
      ! The second index is for prime, noprime, 2Dprime.
      ! The first index is xyz for prime and noprime while it is xy,xz,yz for
      !   the 2Dprime case.

   ! Full solution variables.
   real (kind=double), dimension(3) :: soln
   real (kind=double) :: sum_xyz_A_2, sum_xyz_B_2
   real (kind=double) :: zeta, xi
   real (kind=double) :: preFactorOL, preFactor02, preFactor04, preFactor22
   real (kind=double), dimension (3) :: d, S

   ! Before we proceed with the calculation we need to understand a bit more
   !   about exactly what is being computed. The form of the integration is:
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        coeff * dellR2 ( dell^2
   !        [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))]) dxdydz
   ! We use Pxyz for points in space; Rxyz1,2 for atomic sites 1 and 2;
   !   zeta1,2 for the exponential decay factors; lxyz1,2 for the angular
   !   momentum numbers for atoms 1 and 2; and coeff = -1/2. Note that "dell"
   !   implies derivatives with respect to electron positional coordinates
   !   x,y,z while "dellR2" implies derivatives with respect to the atomic
   !   coordinates of site 2. Derivatives of a Cartesian Gaussian wrt the
   !   atomic site are the negatives of the derivatives wrt the electron
   !   coordinates. I.e., dG/dx = -dG/dBx so that dell = -dellR2.
   !
   ! The dell^2 operator is d^2/dx^2 + d^2/dy^2 + d^2/dz^2 so that dellR2 of
   !   dell^2 (or -dell of dell^2) is:
   !     -[d^3/dx^3 + d/dx d^2/dy^2 + d/dx d^2/dz^2] x-hat
   !     -[d/dy d^2/dx^2 + d^3/dy^3 + d/dy d^2/dz^2] y-hat
   !     -[d/dz d^2/dx^2 + d/dz d^2/dy^2 + d^3/dz^3] z-hat
   ! The integral must be computed over all space for all different possible
   !   values of lxyz1,2 for some arbitrarily chosen test coordinates and
   !   decay rates.
   !
   ! Because of the plus signs in the dell of dell^2 operator we arrive
   !   at three identical integrals of the form (with # in d^3/d#^3 = x, y,
   !   or z):
   ! -SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 *
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        coeff * d^3/d#^3 [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   !
   ! Also, because of the plus signs, we arrive at six other identical
   !   integrals of the form (with #,@ = x,y ; x,z ; y,x ; y,z ; z,x ; z,y):
   ! -SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        coeff * d/d# d^2/d@^2
   !        [(Px-Rx2)**lx2 * (Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*( (Px-Rx2)**2 + (Py-Ry2)**2 + (Pz-Rz2)**2 ))] dxdydz
   !
   ! Now, focusing on just the d^3/dx^3 version of the set of three integrals,
   !   we will pull all terms with Py and Pz out of the dx integral to get:
   ! -SS { [(Py-Ry1)**ly1 * (Pz-Rz1)**lz1 *
   !        exp(-zeta1*((Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !       [(Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !        exp(-zeta2*((Py-Ry2)**2 + (Pz-Rz2)**2 ))] *
   !      S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * coeff * d^3/dx^3 [
   !         (Px-Rx2)**lx2 * exp(-zeta2*((Px-Rx2)**2))]] dx
   !     } dydz
   ! Applying the third derivative, the internal 1D dx integral has the form:
   !   Ix' = S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * coeff *
   !          [lx2*(lx2^2 - 3lx2 + 2) * (Px-Rx2)^(lx2-3) +
   !           (zeta2*-6lx2^2) * (Px-Rx2)^(lx2-1) +
   !           (zeta2^2*(12*lx2 + 12)) * (Px-Rx2)^(lx2+1) +
   !           (zeta2^3*-8) * (Px-Rx2)^(lx2+3)] * exp(-zeta2*(Px - Rx2)**2)
   !        ]
   ! Each of the other integrals (Iy, Iz) will have the form of a simple
   !   1D overlap integral:
   !   Iy = S [(Py-Ry1)^ly1 * exp(-zeta1*((Py-Ry1)^2)) *
   !            (Py-Ry2)^ly2 * exp(-zeta2*((Py-Ry2)^2))] dy
   !   Iz = S [(Pz-Rz1)^lz1 * exp(-zeta1*((Pz-Rz1)^2)) *
   !            (Pz-Rz2)^lz2 * exp(-zeta2*((Pz-Rz2)^2))] dz
   ! The total integral !! for the d^3/dx^3 portion !! is thus Ix' * Iy * Iz
   !   in the x-hat direcetion. Similarly, the d^3/dy^3 portion in the y-hat
   !   direction would be Ix * Iy' * Iz and the d^3/dz^3 portion in the z-hat
   !   direction would be Ix * Iy * Iz'. The Iy' and Iz' are the appropriate
   !   analongs of Ix' and Ix is the analog of the Iy or Iz.
   !
   ! Now, focusing on just the d/d# d^2/d@^2 where (say) #,@ = x,y,
   !   we will pull terms with Pz out of the dx dy integral to get:
   ! -S { [(Pz-Rz1)^lz1 * exp(-zeta1*(Pz-Rz1)^2)] *
   !      [(Pz-Rz2)^lz2 * exp(-zeta2*(Pz-Rz2)^2)] *
   !     SS [(Px-Rx1)^lx1 * (Py-Ry1)^ly1 * 
   !        exp(-zeta1*((Px-Rx1)^2 + (Py-Ry1)^2)) * coeff *
   !        d/dx d^2/dy^2 [(Px-Rx2)^lx2 * (Py-Ry2)^ly2 *
   !        exp(-zeta2*((Px-Rx2)^2 + (Py-Ry2)**2))]] dx dy
   !    } dz
   ! Applying the two second derivatives, the internal 2D dx dy integral has
   !   the form:
   ! Ix'y" = S [(Px-Rx1)^lx1 * (Py-Ry1)^ly1 *
   !    exp(-zeta1*((Px-Rx1)^2 + (Py-Ry1)^2)) * coeff *
   !    lx2*ly2 (ly2 - 1) (Px-Rx2)^(lx2-1)*(Py-Ry2)^(ly2-2) *
   !    zeta2 2lx2 (-2ly2 - 1) (Px-Rx2)^(lx2-1)*(Py-Ry2)^(ly2) *
   !    zeta2^2 4lx2 (Px-Rx2)^(lx2-1)*(Py-Ry2)^(ly2+2) *
   !    zeta2 2ly2 (-ly2 + 1) (Px-Rx2)^(lx2+1)*(Py-Ry2)^(ly2-2) *
   !    zeta2^2 4 (2ly2 + 1) (Px-Rx2)^(lx2+1)*(Py-Ry2)^(ly2) *
   !    zeta2^3 (-8) (Px-Rx2)^(lx2+1)*(Py-Ry2)^(ly2+2) *
   !    exp(-zeta2*((Px-Rx2)^2 + (Py-Ry2)^2))] dx dy
   !
   ! Thus, every total integral in a particular Cartesian direction is a sum
   !   of various independent integrals. We can do 1D overlap integrals for
   !   x, y, z; 1D third derivative integrals for x, y, z; and 2D second
   !   derivative integrals for xy, xz, yx, yz, zx, zy pairs to construct the
   !   total integral.
   ! With regard to the term in each integral that has the **(lx2-2) form
   !   (or **(ly2-2) or **(lz2-2), or -1 or -3 form), some special care must
   !   be taken. This term represents an angular momentum and the integral will
   !   have the form of an overlap integral. Because we cannot have a negative
   !   angular momentum we must discard this term when lx2-2 < 0, etc.


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

   ! Initialize terms used to compute analytical factors.
   coeff = -0.5d0 ! The - sign is because we did not convert d/dB_x into d/dx.
   zeta = a1 + a2
   xi = a1 * a2 / zeta
   d = A - B

   ! Initialize numerical mesh quantities.
   start_pos = -cell_size
   num_steps = cell_size * 2.0d0 / step_size + 1  ! +1 accounts for xyz=zero.

   ! Initialize a counter of the triad pq pairs.
   h = 0

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

         ! Assign l1 and l2 values for each primitive Cartesian gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Initialize sum solution.
         soln(:) = 0.0d0

         ! Start a loop over 1D coordinates.
         do i = 0, num_steps

            ! Compute the current x position, distance, and sum.
            xyz(1) = (start_pos + (i*step_size))

            do j = 0, num_steps

               ! Compute the current y position, distance, and sum.
               xyz(2) = (start_pos + (j*step_size))

               do k = 0, num_steps

                  ! Compute the current z position, distance, and sum.
                  xyz(3) = (start_pos + (k*step_size))
                  sum_xyz_A_2 = sum((xyz(:)-A(:))**2)
                  sum_xyz_B_2 = sum((xyz(:)-B(:))**2)

                  soln(:) = soln(:) + &
&exp(-a1*sum_xyz_A_2) * exp(-a2*sum_xyz_B_2) * (&
&8*a2**3*(xyz(:)-B(:))*(sum_xyz_B_2)&
&-20*a2**2*(xyz(:)-B(:)))
! Above is for the <0|d/dB del|0> case.


!(-a2**3*(-2*xyz(1)+2*B(1))*(2*xyz(1)-2*B(1))**2*exp(-a2*((xyz(1)-B(1))**2+(xyz&
!&(2)-B(2))**2+(xyz(3)-B(3))**2))-a2**3*(-2*xyz(1)+2*B(1))*(2*xyz(2)-2*B(2))**2*&
!&exp(-a2*((xyz(1)-B(1))**2+(xyz(2)-B(2))**2+(xyz(3)-B(3))**2))-a2**3*(-2*xyz(1)&
!&+2*B(1))*(2*xyz(3)-2*B(3))**2*exp(-a2*((xyz(1)-B(1))**2+(xyz(2)-B(2))**2+(xyz(&
!&3)-B(3))**2))+a2**2*(-8*xyz(1)+8*B(1))*exp(-a2*((xyz(1)-B(1))**2+(xyz(2)-B(2))&
!&**2+(xyz(3)-B(3))**2))+6*a2**2*(-2*xyz(1)+2*B(1))*exp(-a2*((xyz(1)-B(1))**2+(x&
!&yz(2)-B(2))**2+(xyz(3)-B(3))**2)))*exp(-a1*((xyz(1)-A(1))**2+(xyz(2)-A(2))**2+&
!&(xyz(3)-A(3))**2))

! All above is for the <0|d/dB del|0> case.

               enddo
            enddo
         enddo

         ! Multiply the results by the mesh step size.
         soln(:) = soln(:) * step_size**3

         pc(q,p,:) = soln(:) * coeff

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
            
"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [8, 3, ""], 0, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [8, 3, ""], 1, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [8, 3, ""], 2, num_conversion,
                       num_conversion)


    # Print the subroutine foot and the subsequent functions.
    foot = """
   end subroutine dkinetic2CIntgNum
"""
    f.write(foot)


def print_test_dnuclearcb_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
   subroutine dnuclear3CIntgNumCB(a1,a2,a3,A,B,C,pc,sh,cell_size,step_size)

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
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)},3), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)},3), intent(out) :: sh" + """
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, h, i, j, k
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)},3) :: triads" + """
   integer, dimension (""" + f"{len(conversion)},2,3) :: conversion" + """
   integer, dimension (3) :: l1, l2
   real (kind=double), dimension(3) :: soln
   real (kind=double) :: start_pos, r, xyz_sum_coeff
   real (kind=double), allocatable, dimension (:,:) :: xyz, xyz_sum
   real (kind=double), allocatable, dimension (:,:) :: xyz_sum_plus1
   real (kind=double), allocatable, dimension (:,:) :: xyz_sum_minus1
   real (kind=double), allocatable, dimension (:,:) :: B_dist, C_dist
   real (kind=double), allocatable, dimension (:,:) :: C_dist_sqrd

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

   ! Original equation for d/dx nuclear from sympy:
   ! (xyz(1)-A(1))**l1(1)*(xyz(2)-A(2))**l1(2)*(xyz(3)-A(3))**l1(3)
   ! *exp(-a1*((xyz(1)-A(1))**2+(xyz(2)-A(2))**2+(xyz(3)-A(3))**2))
   ! *exp(-a2*((xyz(1)-B(1))**2+(xyz(2)-B(2))**2+(xyz(3)-B(3))**2))
   ! *exp(-a3*((xyz(1)-C(1))**2+(xyz(2)-C(2))**2+(xyz(3)-C(3))**2))
   ! /sqrt((xyz(1)-C(1))**2+(xyz(2)-C(2))**2+(xyz(3)-C(3))**2)
   ! *(
   ! -l2(1)*(xyz(1)-B(1))**(l2(1)-1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)
   ! +2*a2*(xyz(1)-B(1))**(l2(1)+1)*(xyz(2)-B(2))**l2(2)*(xyz(3)-B(3))**l2(3)
   ! )

   ! Alocate space to hold precomputed quantities.
   allocate(xyz(num_steps,3))
   allocate(xyz_sum(num_steps,3))
   allocate(xyz_sum_plus1(num_steps,3))
   allocate(xyz_sum_minus1(num_steps,3))
   allocate(B_dist(num_steps,3))
   allocate(C_dist(num_steps,3))
   allocate(C_dist_sqrd(num_steps,3))

   ! By-hand algebraic rearrangement to the expression below for faster
   !   computation in the triple nested loop environment.

   ! Initialize a counter of the triad pq pairs.
   h = 0

   do p = 1, """ + f"{len(triads)}" + """
      do q = 1, """ + f"{len(triads)}" + """

         ! Assign l1 and l2 values for each gto.
         l1 = triads(q,:)
         l2 = triads(p,:)

         ! Precompute all one dimensional distance and exp factors.
         do i = 1, 3
             do j = 1, num_steps+1
                xyz(j,i) = (start_pos + ((j-1)*step_size))
                B_dist(j,i) = xyz(j,i) - B(i)
                C_dist(j,i) = xyz(j,i) - C(i)
                C_dist_sqrd(j,i) = (xyz(j,i) - C(i))**2
                xyz_sum(j,i) = &
                   & (xyz(j,i)-A(i))**l1(i)*(xyz(j,i)-B(i))**l2(i) &
                   & *exp(-a1*((xyz(j,i)-A(i))**2) - a2*((xyz(j,i)-B(i))**2) &
                   & - a3*(C_dist(j,i))**2)
                xyz_sum_plus1(j,i) = &
                   & (xyz(j,i)-A(i))**l1(i)*(xyz(j,i)-B(i))**(l2(i)+1) &
                   & *exp(-a1*((xyz(j,i)-A(i))**2) - a2*((xyz(j,i)-B(i))**2) &
                   & - a3*(C_dist(j,i))**2)
                if (l2(i) >= 1) then
                   xyz_sum_minus1(j,i) = &
                      & (xyz(j,i)-A(i))**l1(i)*(xyz(j,i)-B(i))**(l2(i)-1) &
                      & * exp(-a1*((xyz(j,i)-A(i))**2) &
                      & - a2*((B_dist(j,i))**2) - a3*(C_dist(j,i))**2)
                endif
             enddo
         enddo

         ! Initialize the solution
         soln(:) = 0.0d0

         do i = 0, num_steps
            do j = 0, num_steps
               do k = 0, num_steps

                  ! Compute the distance between the electron and the nucleus.
                  r = C_dist_sqrd(i,1) + C_dist_sqrd(j,2) + C_dist_sqrd(k,3)

                  ! If the distance is zero, then cycle.
                  if (abs(r) <= 0.00000001d0) cycle

                  soln(1) = soln(1) + xyz_sum_plus1(i,1) * xyz_sum(j,2) &
                        & * xyz_sum(k,3) * 2.0d0 * a2 / sqrt(r)

                  soln(2) = soln(2) + xyz_sum(i,1) * xyz_sum_plus1(j,2) &
                        & * xyz_sum(k,3) * 2.0d0 * a2 / sqrt(r)

                  soln(3) = soln(3) + xyz_sum(i,1) * xyz_sum(j,2) &
                        & * xyz_sum_plus1(k,3) * 2.0d0 * a2 / sqrt(r)

                  if (l2(1) >= 1) then
                     soln(1) = soln(1) + xyz_sum_minus1(i,1) * xyz_sum(j,2) &
                        & * xyz_sum(k,3) * (-l2(1)) / sqrt(r)
                  endif

                  if (l2(2) >= 1) then
                     soln(2) = soln(2) + xyz_sum(i,1) * xyz_sum_minus1(j,2) &
                        & * xyz_sum(k,3) * (-l2(2)) / sqrt(r)
                  endif

                  if (l2(3) >= 1) then
                     soln(3) = soln(3) + xyz_sum(i,1) * xyz_sum(j,2) &
                        & * xyz_sum_minus1(k,3) * (-l2(3)) / sqrt(r)
                  endif

               enddo
            enddo
         enddo

         ! Multiply the results by the mesh step size.
         soln(:) = soln(:) * step_size**3

         pc(q,p,:) = soln(:)

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

   deallocate(xyz)
   deallocate(xyz_sum)
   deallocate(xyz_sum_plus1)
   deallocate(xyz_sum_minus1)
   deallocate(B_dist)
   deallocate(C_dist)
   deallocate(C_dist_sqrd)

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [4, 3, ""], 0, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [4, 3, ""], 1, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [4, 3, ""], 2, num_conversion,
                       num_conversion)


    # Print the subroutine foot.
    foot = """
   end subroutine dnuclear3CIntgNumCB
"""
    f.write(foot)


def print_test_delectroncb_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
   subroutine delectron3CIntgNumCB(a1,a2,a3,A,B,C,pc,sh,cell_size,step_size)

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
   real (kind=double), dimension (""" \
    + f"{len(triads)},{len(triads)},3), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
    + f"{len(conversion)},{len(conversion)},3), intent(out) :: sh" + """
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i, j
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
   integer, dimension (3) :: l1, l2
   real (kind=double) :: start_pos, curr_pos
   real (kind=double), dimension (3) :: xyz
   real (kind=double), dimension (3,2) :: xyz_I ! Indices=xyz, noprime||prime

   ! Before we proceed with the calculation we need to understand a bit more
   !   about exactly what is being computed. The form of the integration is:
   ! SSS [(Px-Rx1)**lx1 * (Py-Ry1)**ly1 * (Pz-Rz1)**lz1 * 
   !        exp(-zeta1*( (Px-Rx1)**2 + (Py-Ry1)**2 + (Pz-Rz1)**2 ))] *
   !        d/dR_xyz2 [
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
   !   the scenario when the middle term is at site C and the last term is
   !   also at site "B". The B is in quotes to emphasize the following idea:
   !   Consider that we want to take the derivative of the Hamiltonian with
   !   respect to the coordinate of, say, atom #7 in order to get the force on
   !   atom #7. For the three center overlap contribution to the Hamiltonian,
   !   (commonly expressed as <A|C|B>) we find that atom #7 will be used in
   !   each of the A, C, and B terms. Sometimes it will be <7|C|B> with
   !   arbitrary other atoms for C and B. Other times it will be <A|7|B> or
   !   <A|C|7> or <A|7|7>. Presently, we are dealing with the <A|C|7> case.
   ! Note, importantly, that the derivative is with respect to the atomic
   !   coordinate of the second basis function only.
   ! Now, focusing on just the d/dRx2 version of the set of triple integrals,
   !   we will pull all terms with Py and Pz out of the dx integral to get:
   ! SS { [(Py-Ry1)**ly1 * (Pz-Rz1)**lz1 *
   !       exp(-zeta1*((Py-Ry1)**2 + (Pz-Rz1)**2 ))]*
   !       [exp(-zeta3*((Py-Ry3)**2 + (Pz-Rz3)**2 ))]*
   !       [(Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !       exp(-zeta2*((Py-Ry2)**2 + (Pz-Rz2)**2 ))] *
   !     S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * d/dRx2 [
   !     exp(-zeta3*((Px-Rx3)**2)) * (Px-Rx2)**lx2 * exp(-zeta2*((Px-Rx2)**2))
   !    ]] dx} dydz
   ! Applying the derivative, the internal 1D dx integral has the form:
   !   Ix' = S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) *                
   !           [-lx2*(Px-Rx2)**(lx2-1)
   !            + 2*zeta2*(Px-Rx2)**(lx2+1)] 
   !           * exp(-zeta2*(Px-Rx2)**2) * exp(-zeta3*(Px-Rx3)**2)
   !        ]
   ! Each of the other integrals (Iy, Iz) will have the form of a simple
   !   1D overlap integral (without a derivative):
   !   Iy = S [(Py-Ry1)**ly1 * exp(-zeta1*((Py-Ry1)**2)) *
   !           (Py-Ry2)**ly2 * exp(-zeta2*((Py-Ry2)**2)) * 
   !                           exp(-zeta3*((Py-Ry3)**2))] dy
   !   Iz = S [(Pz-Rz1)**lz1 * exp(-zeta1*((Pz-Rz1)**2)) *
   !           (Pz-Rz2)**lz2 * exp(-zeta2*((Pz-Rz2)**2)) *
   !                           exp(-zeta3*((Pz-Rz3)**2))] dz
   ! The total integral !! for the d/dRx2 version !! is thus Ix' * Iy * Iz.
   !   Similarly, the total integrals for d/dRy2 and d/dRz2 will have the
   !   form: Ix*Iy'*Iz and Ix*Iy*Iz' where the Iy' and Iz' are the appropriate
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
         xyz_I(:,:) = 0.0d0

         ! Start a loop over x coordinates.
         do i = 0, num_steps
            curr_pos = start_pos + (i*step_size)

            ! Compute the no-prime integrals first.
            do j = 1, 3
               xyz_I(j,2) = xyz_I(j,2) &
                     & + noPrimeDElectronCB(step_size, curr_pos, A(j), &
                     & B(j), C(j), a1, a2, a3, l1(j), l2(j))
            enddo

            ! Compute the prime integrals second.
            do j = 1, 3
               xyz_I(j,1) = xyz_I(j,1) &
                     & + primeDElectronCB(step_size, curr_pos, A(j), &
                     & B(j), C(j), a1, a2, a3, l1(j), l2(j))
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
    lib.print_pc_to_sh(conversion, f, [9, 3, ""], 0, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [9, 3, ""], 1, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [9, 3, ""], 2, num_conversion,
                       num_conversion)

    # Print the subroutine foot.
    foot = """
   end subroutine delectron3CIntgNumCB


   function noPrimeDElectronCB(step_size, curr_pos, A, B, C, &
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
      real (kind=double) :: noPrimeDElectronCB

      ! Compute the integral part.
      noPrimeDElectronCB = step_size &
            & * (curr_pos - A)**l1 * (curr_pos - B)**l2 &
            & * exp(-a1*(curr_pos - A)**2) * exp(-a2*(curr_pos - B)**2) &
            & * exp(-a3*(curr_pos - C)**2)

      return
 
   end function noPrimeDElectronCB


   function primeDElectronCB(step_size, curr_pos, A, B, C, a1, a2, a3, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos
      real (kind=double), intent(in) :: A, B, C, a1, a2, a3
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: primeDElectronCB

      ! Compute each "internal" term of the prime integral. Compare each of
      !   these lines with the 1D electron derivative matrix equation produced
      !   by the osrecurintg_makenum.py script (appropriately separated into
      !   terms) and the expression in Nuha's dissertation. FIX: Check
      !   equation number from Nuha's dissertation.

      ! Last term in the reduced equation produced by sympy. FIX: As above.
      primeDElectronCB = 2.0d0 * a2 * (curr_pos - B)**(l2+1)

      ! First term in the equation produced by sympy (after the exponentials
      !   of course). As mentioned above, we only compute this if the
      !   angular momentum will allow it.
      if (l2 >= 1) then
         primeDElectronCB = primeDElectronCB - l2 * (curr_pos - B)**(l2-1)
      endif

      ! Multiply prime integral by the preceeding primitive gaussian
      !   coefficient and exponential and multiply by the succeeding
      !   exponential. (We have already multiplied by the succeeding
      !   primitive gaussian coefficient in the above lines.)
      primeDElectronCB = primeDElectronCB &
            & * (curr_pos-A)**l1 * exp(-a1*(curr_pos-A)**2) &
            & * exp(-a2*(curr_pos-B)**2) * exp(-a3*(curr_pos-C)**2)

      ! Finally, multiply by the step size.
      primeDElectronCB = primeDElectronCB * step_size

      return

   end function primeDElectronCB
"""
    f.write(foot)


def print_test_delectronbb_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
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
   real (kind=double), dimension (""" \
    + f"{len(triads)},{len(triads)},3), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
    + f"{len(conversion)},{len(conversion)},3), intent(out) :: sh" + """
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i, j
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
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

"""
    f.write(head)

    # Print the section to assemble the sh matrix.
    num_conversion = len(conversion)
    lib.print_pc_to_sh(conversion, f, [9, 3, ""], 0, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [9, 3, ""], 1, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [9, 3, ""], 2, num_conversion,
                       num_conversion)

    # Print the subroutine foot.
    foot = """
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
"""
    f.write(foot)


def print_test_delectronbc_num(conversion, triads, f):

    # Print the subroutine header for the numerical portion.
    head = """
   subroutine delectron3CIntgNumBC(a1,a2,a3,A,B,C,pc,sh,cell_size,step_size)

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
   real (kind=double), dimension (""" \
    + f"{len(triads)},{len(triads)},3), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
    + f"{len(conversion)},{len(conversion)},3), intent(out) :: sh" + """
   real (kind=double), intent (in) :: cell_size, step_size

   ! Define local variables.
   integer :: p, q, i, j
   integer :: num_steps
   integer, dimension (""" + f"{len(triads)}" + """,3) :: triads
   integer, dimension (""" + f"{len(conversion)}" + """,2,3) :: conversion
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
   !   at another site C. The B is in quotes to emphasize the following idea:
   !   Consider that we want to take the derivative of the Hamiltonian with
   !   respect to the coordinate of, say, atom #7 in order to get the force on
   !   atom #7. For the three center overlap contribution to the Hamiltonian,
   !   (commonly expressed as <A|C|B>) we find that atom #7 will be used in
   !   each of the A, C, and B terms. Sometimes it will be <7|C|B> with
   !   arbitrary other atoms for C and B. Other times it will be <A|7|B> or
   !   <A|C|7> or <A|7|7>. Presently, we are dealing with the <A|7|B> case.
   ! Note, importantly, that the derivative is with respect to the atomic
   !   coordinate of the middle term only. Therefore, we change the derivative
   !   to be with respect to the R_xyz3 coordinate.
   ! Now, focusing on just the d/dRx3 version of the set of triple integrals,
   !   we will pull all terms with Py and Pz out of the dx integral to get:
   ! SS { [(Py-Ry1)**ly1 * (Pz-Rz1)**lz1 *
   !       exp(-zeta1*((Py-Ry1)**2 + (Pz-Rz1)**2 ))]*
   !       [exp(-zeta3*((Py-Ry3)**2 + (Pz-Rz3)**2 ))]*
   !       [(Py-Ry2)**ly2 * (Pz-Rz2)**lz2 *
   !       exp(-zeta2*((Py-Ry2)**2 + (Pz-Rz2)**2 ))] *
   !     S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) * d/dRx3 [
   !     exp(-zeta3*((Px-Rx3)**2)) * (Px-Rx2)**lx2 * exp(-zeta2*((Px-Rx2)**2))
   !    ]] dx} dydz
   ! Applying the derivative, the internal 1D dx integral has the form:
   !   Ix' = S [(Px-Rx1)**lx1 * exp(-zeta1*((Px-Rx1)**2)) *                
   !           [2*zeta3*(Px-Rx3)*(Px-Rx2)**lx2]
   !           * exp(-zeta2*(Px-Rx2)**2) * exp(-zeta3*(Px-Rx3)**2)
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
         xyz_I(:,:) = 0.0d0

         ! Start a loop over x coordinates.
         do i = 0, num_steps
            curr_pos = start_pos + (i*step_size)

            ! Compute the no-prime integrals first.
            do j = 1, 3
               xyz_I(j,2) = xyz_I(j,2) &
                     & + noPrimeDElectronBC(step_size, curr_pos, A(j), &
                     & B(j), C(j), a1, a2, a3, l1(j), l2(j))
            enddo

            ! Compute the prime integrals second.
            do j = 1, 3
               xyz_I(j,1) = xyz_I(j,1) &
                     & + primeDElectronBC(step_size, curr_pos, A(j), &
                     & B(j), C(j), a1, a2, a3, l1(j), l2(j))
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
    lib.print_pc_to_sh(conversion, f, [9, 3, ""], 0, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [9, 3, ""], 1, num_conversion,
                       num_conversion)
    lib.print_pc_to_sh(conversion, f, [9, 3, ""], 2, num_conversion,
                       num_conversion)

    # Print the subroutine foot.
    foot = """
   end subroutine delectron3CIntgNumBC


   function noPrimeDElectronBC(step_size, curr_pos, A, B, C, &
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
      real (kind=double) :: noPrimeDElectronBC

      ! Compute the integral part.
      noPrimeDElectronBC = step_size &
            & * (curr_pos - A)**l1 * (curr_pos - B)**l2 &
            & * exp(-a1*(curr_pos - A)**2) * exp(-a2*(curr_pos - B)**2) &
            & * exp(-a3*(curr_pos - C)**2)

      return
 
   end function noPrimeDElectronBC


   function primeDElectronBC(step_size, curr_pos, A, B, C, a1, a2, a3, l1, l2)

      ! Use necessary modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), intent(in) :: step_size, curr_pos
      real (kind=double), intent(in) :: A, B, C, a1, a2, a3
      integer, intent(in) :: l1, l2

      ! Define local and return variables.
      real (kind=double) :: primeDElectronBC

      ! Compute each "internal" term of the prime integral. Compare each of
      !   these lines with the 1D electron derivative matrix equation produced
      !   by the osrecurintg_makenum.py script (appropriately separated into
      !   terms) and the expression in Nuha's dissertation. FIX: Check
      !   equation number from Nuha's dissertation.

      ! This is the only term in the sympy script-produced equation.
      !   FIX: As above.
      primeDElectronBC = 2.0d0 * a3 * (curr_pos - C)

      ! Multiply prime integral by the preceeding primitive gaussian
      !   coefficient and exponential and multiply by the succeeding
      !   exponential. (We have already multiplied by the succeeding
      !   primitive gaussian coefficient in the above lines.)
      primeDElectronBC = primeDElectronBC &
            & * (curr_pos-A)**l1 * exp(-a1*(curr_pos-A)**2) &
            & * (curr_pos-B)**l2 * exp(-a2*(curr_pos-B)**2) &
            & * exp(-a3*(curr_pos-C)**2)

      ! Finally, multiply by the step size.
      primeDElectronBC = primeDElectronBC * step_size

      return

   end function primeDElectronBC
"""
    f.write(foot)


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. In this case, there is no program. This
    #   should only be imported as a module.
    main()
