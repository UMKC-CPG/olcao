module InitRhoSubs

public

contains

subroutine initRho

   ! Import the necessary definitions
   use O_Kinds

   ! Import the necessary data modules
   use O_RadialGrid
   use ChargeDensityData
   use AtomData

   implicit none

   ! Define local variables.
   real (kind=double) :: aaSomething

   ! Allocate space for the up and down charge density arrays, the total
   !   charge density, the 1st and 2nd derivatives of the total, the total,
   !   1st and 2nd charge density divided by the radial value at each point,
   !   and the core charge density (not really used though).
   allocate (rhoDn       (numRadialPoints)) ! Up
   allocate (rhoUp       (numRadialPoints)) ! Down
   allocate (rhoTotal    (numRadialPoints)) ! Up + Down
   allocate (dr1RhoTotal (numRadialPoints)) ! 1st derivative
   allocate (dr2RhoTotal (numRadialPoints)) ! 2nd derivative
   allocate (rhoOverR    (numRadialPoints)) ! (Up + Down) / radial value
   allocate (dr1RhoOverR (numRadialPoints)) ! 1st derivative / radial value
   allocate (dr2RhoOverR (numRadialPoints)) ! 2nd derivative / radial value
   allocate (rhoCore     (numRadialPoints)) ! Core

   ! Allocate space to hold the cumulative integrals of the total charge
   !   density and the total charge density / radial value.  The integrals are
   !   cumulative in the sense that for each radial value r(i) the integral is
   !   computed over the range r=0..r(i) and saved in these arrays.
   allocate (cumulInt1   (numRadialPoints)) ! rhoTotal
   allocate (cumulInt2   (numRadialPoints)) ! rhoTotal / radial value

   ! I don't know what this is all about.
   aaSomething = 0.5_double

   ! Assign initial values for the charge density according to some crazy
   !   equation without ANY documentation.  ARG!  Well, it does say that
   !   rhoDn =  2 pi r**2 rho(r) but that doesn't help a whole lot since it
   !   does not have an obvious connection to the code in use.
   rhoDn(:) = electronCharge * aaSomething**3 * exp(-aaSomething * &
         & radialPoints(:)) * radialPoints(:)**2 / 4.0_double
   rhoUp(:) = rhoDn(:)

end subroutine initRho

end module InitRhoSubs
