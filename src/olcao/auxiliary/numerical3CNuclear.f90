! Author:  Ben Walker
!
! This program is a numerical integrator that provides test numbers
! for Nuclear Attraction Integrals, AKA Nuclear Potential Integrals.
! This version of the program is intended to output the numerical
! integration solutions for the 20x20 array of basis functions.

program numerical3CNuclear
use O_kinds
implicit none

! Declare variables
real (kind=double) :: a1, a2, a3
real (kind=double) :: x, y, z, soln, sSum, mesh
real (kind=double), dimension(3) :: A, B, C
integer :: i, j, k, p, q
real (kind=double) :: r
integer, dimension(3) :: l1, l2
integer, dimension(20,3) :: array 

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

! Define the constants
a1 = 0.3d0
a2 = 0.3d0
a3 = 0.3d0
A = 0d0
B = 0d0
C = 0d0
mesh = 0.2d0 ! Width of interval.

! Open output file
open(2,file='testOutput',status='new')

! Loop over the x, y, and z triad values
do p = 1, 20
  do q = 1, 20
    ! Properly assign l1 and l2 values for each dimension
    l1 = array(p,:)
    l2 = array(q,:)

    ! Initialize sum variable
    sSum = 0
  
    ! Computed f(x,y,z) for each mesh point.
    do i = 0, 40
      x = (-4d0 + i*mesh)
      do j = 0, 40
        y = (-4d0 + j*mesh)
        do k = 0, 40
          z = (-4d0 + k*mesh)
          r = ((x - C(1))**2) + ((y - C(2))**2) + ((z - C(3))**2)
          ! Avoid dividing by zero.
          if (r == 0) cycle
          ! Multiply each f(x,y,z) value by the mesh size and store the result.
          soln = (mesh**3)*((x - A(1))**l1(1))*((x - B(1))**l2(1))&
            &*exp(-a1*(x - A(1))**2)*exp(-a2*(x - B(1))**2)*exp(-a3*(x - C(1))**2)&
            &*((y - A(2))**l1(2))*((y - B(2))**l2(2))&
            &*exp(-a1*(y - A(2))**2)*exp(-a2*(y - B(2))**2)*exp(-a3*(y - C(2))**2)&
            &*((z - A(3))**l1(3))*((z - B(3))**l2(3))&
            &*exp(-a1*(z - A(3))**2)*exp(-a2*(z - B(3))**2)*exp(-a3*(z - C(3))**2)&
            &/sqrt(r)
          ! Sum up the results.
          sSum = sSum + soln
        end do !k
      end do !j
    end do !i
    write(2,*) sSum
  end do !q
end do !p
close(2)
end program numerical3CNuclear
