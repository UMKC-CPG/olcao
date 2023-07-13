! Author:  Ben Walker
!
! This program is a numerical integrator that provides test numbers
! for 3-center Electronic Potential Integrals.  This version of the program is 
! intended to output the numerical integration solutions for 
! the 20x20 array of basis functions.

program numerical3COverlap
use O_kinds
implicit none

! Declare variables
real (kind=double) :: x_sum, y_sum, z_sum
real (kind=double) :: x_solution, y_solution, z_solution, x, y, z
real (kind=double) :: a1, a2, a3
real (kind=double), dimension(3) :: A, B, C
real (kind=double) :: step, num_steps, max_size
integer :: i, p, q
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
a1 = 0.50d0
a2 = 0.95d0
a3 = 0.30d0
A(1) = 1.0d0
A(2) = 2.0d0
A(3) = 3.0d0
B(1) = 6.0d0
B(2) = 5.0d0
B(3) = 4.0d0
C(1) = 1.5d0
C(2) = 1.5d0
C(3) = 1.5d0

! Open output file
open(2,file='testOutput',status='new')

max_size = 30.0d0
step = 0.0001
num_steps = (2*max_size) / step + 1

! Loop over the x, y, and z triad values
do p = 1, 4
  do q = 1, 4
    ! Properly assign l1, l2, and l3 values for each dimension
    l1 = array(q,:)
    l2 = array(p,:)

    ! Initialize sum variables
    x_sum = 0
    y_sum = 0
    z_sum = 0
   
    do i = 0, num_steps
      x = (-max_size + (i*step))
      y = (-max_size + (i*step))
      z = (-max_size + (i*step))
      x_solution = step*((x - A(1))**l1(1))*((x - B(1))**l2(1))*&
        &exp(-a1*(x - A(1))**2)*exp(-a2*(x - B(1))**2)*exp(-a3*(x - C(1))**2)
      x_sum = x_sum + x_solution
      y_solution = step*((y - A(2))**l1(2))*((y - B(2))**l2(2))*&
        &exp(-a1*(y - A(2))**2)*exp(-a2*(y - B(2))**2)*exp(-a3*(y - C(2))**2)
      y_sum = y_sum + y_solution
      z_solution = step*((z - A(3))**l1(3))*((z - B(3))**l2(3))*&
        &exp(-a1*(z - A(3))**2)*exp(-a2*(z - B(3))**2)*exp(-a3*(z - C(3))**2)
      z_sum = z_sum + z_solution
    end do !i
    write(2,*) x_sum*y_sum*z_sum
  end do !q
end do !p
close(2)
end program numerical3COverlap
