! Author:  Ben Walker
!
! This program is a numerical integrator that provides test numbers
! for Momentum Matrix Integrals.  This version of the program is 
! intended to output the numerical integration solutions for 
! the 20x20 array of basis functions.

program numerical2CMomentum
use O_kinds
implicit none

! Declare variables
real (kind=double) :: x, y, z
real (kind=double) :: x_sum1,x_sum2,x_sum3,x_sum4,x_sum5,x_sum6,x_sum7
real (kind=double) :: x_solution1,x_solution2,x_solution3,x_solution4
real (kind=double) :: x_solution5,x_solution6,x_solution7
real (kind=double) :: y_sum1,y_sum2,y_sum3,y_sum4,y_sum5,y_sum6,y_sum7
real (kind=double) :: y_solution1,y_solution2,y_solution3,y_solution4
real (kind=double) :: y_solution5,y_solution6,y_solution7
real (kind=double) :: z_sum1,z_sum2,z_sum3,z_sum4,z_sum5,z_sum6,z_sum7
real (kind=double) :: z_solution1,z_solution2,z_solution3,z_solution4
real (kind=double) :: z_solution5,z_solution6,z_solution7
real (kind=double) :: a1, a2
real (kind=double), dimension(3) :: A, B
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
A = 1.0d0
B = 3.0d0

! Open output file
open(2,file='testOutput',status='new')

! Loop over the x, y, and z triad values
do p = 1, 20
  do q = 1, 20
    ! Properly assign l1 and l2 values for each dimension
    l1 = array(p,:)
    l2 = array(q,:)

    ! Initialize sum variables
    x_sum1 = 0
    x_sum2 = 0
    x_sum3 = 0
    x_sum4 = 0
    x_sum5 = 0
    x_sum6 = 0
    x_sum7 = 0
    y_sum1 = 0
    y_sum2 = 0
    y_sum3 = 0
    y_sum4 = 0
    y_sum5 = 0
    y_sum6 = 0
    y_sum7 = 0
    z_sum1 = 0
    z_sum2 = 0
    z_sum3 = 0
    z_sum4 = 0
    z_sum5 = 0
    z_sum6 = 0
    z_sum7 = 0

    do i = 0, 100
      x = (-10.0d0 + (i * 0.20d0))
      y = (-10.0d0 + (i * 0.20d0))
      z = (-10.0d0 + (i * 0.20d0))

      ! Since the Momentum Matrix Integral is really made up of a linear combination
      ! of two 2-center overlap integrals, some addional customization is
      ! necessary here.  There will be two segments, and the final result
      ! will add them together.

      ! ** SEGMENT 1 **
      if ((l2(1) - 1) >= 0) then
        x_solution1 = 0.2d0*l2(1)*((x-A(1))**l1(1))*((x-B(1))**(l2(1)-1))*&
          &exp(-0.5d0*(x-A(1))**2)*exp(-0.95d0*(x-B(1))**2)
        x_sum1 = x_sum1 + x_solution1
        y_solution1 = 0.2d0*((y-A(2))**l1(2))*((y-B(2))**l2(2))*&
          &exp(-0.5d0*(y-A(2))**2)*exp(-0.95d0*(y-B(2))**2)
        y_sum1 = y_sum1 + y_solution1
        z_solution1 = 0.2d0*((z-A(3))**l1(3))*((z-B(3))**l2(3))*&
          &exp(-0.5d0*(z-A(3))**2)*exp(-0.95d0*(z-B(3))**2)
        z_sum1 = z_sum1 + z_solution1
      end if
      
      ! ** SEGMENT 2 **
      x_solution2 = 0.2d0*2*a2*((x-A(1))**l1(1))*((x-B(1))**(l2(1)+1))*&
        &exp(-0.5d0*(x-A(1))**2)*exp(-0.95d0*(x-B(1))**2)
      x_sum2 = x_sum2 + x_solution2
      y_solution2 = 0.2d0*((y-A(2))**l1(2))*((y-B(2))**l2(2))*&
        &exp(-0.5d0*(y-A(2))**2)*exp(-0.95d0*(y-B(2))**2)
      y_sum2 = y_sum2 + y_solution2
      z_solution2 = 0.2d0*((z-A(3))**l1(3))*((z-B(3))**l2(3))*&
        &exp(-0.5d0*(z-A(3))**2)*exp(-0.95d0*(z-B(3))**2)
      z_sum2 = z_sum2 + z_solution2
    end do !i
    write(2,*) x_sum1*y_sum1*z_sum1 - x_sum2*y_sum2*z_sum2
  end do !q
end do !p
close(2)
end program numerical2CMomentum
