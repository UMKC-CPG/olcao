program atomSCF

   ! Import the necessary data modules
   use ExecutionData
   use PotData

   ! Import the necessary subroutine modules
   use ReadSCFInputSubs
   use ImplicitSCFInputSubs
   use InitRhoSubs
   use ComputeIonicPotSubs
   use ComputeElectPotSubs
   use ComputeOrbitalsSubs
   use ConvergenceSubs
   use MixPotentialSubs
   use AtomicTotalEnergySubs

   implicit none

   ! Read the input data
   call readSCFInput

   ! Compute the implicit input data.
   call implicitSCFInput

   ! Initialize the initial charge density.
   call initRho

   ! Set up the ionic potential
   call computeIonicPot

   ! Set up the electric potential
   call computeElectPot

   ! Copy the computed Hartree potential to a secondary array.  The 'i' and 'o'
   !   used in vid and vod from the original program code have an unknown
   !   meaning to me.  In and out?  Who knows!
   hartreePotDn2(:) = hartreePotDn(:)
   hartreePotUp2(:) = hartreePotUp(:)

   ! Begin the self consistant iteration loop
   do while (.true.)

      ! Increment the iteration count.  (Recall it starts at 0.)
      iteration = iteration + 1

      ! Record the progress so far.
      if (mod(iteration,10) == 0) then
         write (6,ADVANCE="NO",FMT="(a1)") "|"
      else
         write (6,ADVANCE="NO",FMT="(a1)") "."
      endif
      if (mod(iteration,50) == 0) then
         write (6,*) " ",iteration
      endif
      call flush (6)

      ! Compute the orbitals
      call computeOrbitals

      ! Set up the output electronic potential from the charge density (rho).
      call computeElectPot

      ! Compute the convergence of the system.
      call checkConvergence

      ! Mix the new electronic potential of this iteration with the last.
      call mixPotential

      ! Check for convergence and compute the final potential if necessary.
      if ((converged == 1) .or. (iteration > maxIteration)) then
         call computeOrbitals
         call computeElectPot
         exit
      endif
   enddo

   ! Finalize the tracking with a hard return if necessary.
   if (mod(iteration,50) /= 0) then
      write (6,*)
   endif
   call flush (6)

   ! Compute the final total energy for this atom.
   call atomicTotalEnergy

   ! Produce the output.
   call printConvgPot

end program atomSCF
