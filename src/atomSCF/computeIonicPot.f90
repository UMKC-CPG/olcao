module ComputeIonicPotSubs

public

contains


! Again, there is not much in the way of theory to explain why the things done
!   in this subroutine are done.  I guess we just have to re figure it all out.
subroutine computeIonicPot

   ! Import the necessary modules.
   use O_Kinds
   use AtomData
   use O_RadialGrid
   use PotData

   ! Define local variables
   integer :: i,j ! Loop index variables.

   ! Allocate space for the up and down ionic potential.
   allocate (ionicPotDn (maxQNl,numRadialPoints))
   allocate (ionicPotUp (maxQNl,numRadialPoints))

   ! Initialize the ionic potential with 2 times the nuclear charge.  Recall
   !   that the ionicPotDn/Up is the ionic potential times r.
   ionicPotDn(:,:) = -2 * atomicNumber
   ionicPotUp(:,:) = -2 * atomicNumber

   ! Abort if the shell charge is zero indicating that atom is neutral.
   if (shellCharge == 0.0_double) return

   do i = 1, maxQNl
      do j = 1, numRadialPoints
         if (radialPoints(j) > shellRadius) then
            ionicPotDn(i,j) = ionicPotDn(i,j) - 2 * shellCharge
            ionicPotUp(i,j) = ionicPotUp(i,j) - 2 * shellCharge
         else

            ! In the original program the radialPoints(j) was actually
            !   radialPoints(i), but I am not sure how that can be right.

            ionicPotDn(i,j) = ionicPotDn(i,j) - 2 * shellCharge * &
                  & radialPoints(j) / shellRadius
            ionicPotUp(i,j) = ionicPotUp(i,j) - 2 * shellCharge * &
                  & radialPoints(j) / shellRadius
         endif
      enddo
   enddo
end subroutine computeIonicPot

end module ComputeIonicPotSubs
