module AtomicTotalEnergySubs

public

contains

subroutine atomicTotalEnergy

   ! Import necessary modules
   use O_Kinds
   use PotData
   use EnergyData
   use AtomData
   use O_RadialGrid

   implicit none

   ! Define local variables.
   character, dimension(5) :: angMomLabels
   integer :: i

   ! Initialize the labels of the angular momentum quantum numbers.
   angMomLabels(1) = 's'
   angMomLabels(2) = 'p'
   angMomLabels(3) = 'd'
   angMomLabels(4) = 'f'
   angMomLabels(5) = 'g'

   ! The total energy array has 10 indices that each represent some component
   !   of the total energy (or the total energy itself).
   ! totalEnergy(1) = sum of eigenvalues (ev)
   ! totalEnergy(2) = sum of orbital kinetic energies (orbitalKE)
   ! totalEnergy(3) = el-ion interaction from sum of orbital
   !                     potential energies (orbitalKE)
   ! totalEnergy(4) = electrostatic el-el interaction (from computeElectPot)
   ! totalEnergy(5) = vxc (exchange-correlation) correction to sum
   !                     of eigenvalues (from computeElectPot)
   ! totalEnergy(6) = 3 * vc - 4 * ec
   !                     correction term for virial theorem
   !                     when correlation is included (from computeElectPot)
   ! totalEnergy(7) = exchange and correlation energy (from computeElectPot)
   ! totalEnergy(8) = kinetic energy from eigenvalues
   ! totalEnergy(9) = potential energy
   ! totalEnergy(10)= total energy

   ! Compute totalEnergy(1,2,3)
   totalEnergy(1) = sum (orbitalCharge(:) * eigenValues(:))
   totalEnergy(2) = sum (orbitalCharge(:) * orbitalKE(:))
   totalEnergy(3) = sum (orbitalCharge(:) * orbitalPE(:))


   ! Compute the kinetic energy.
   totalEnergy(8) = totalEnergy(1) - totalEnergy(3) - &
         & 2.0_double*totalEnergy(4) - totalEnergy(5)

   ! Compute the potential energy.
   totalEnergy(9) = totalEnergy(3) + totalEnergy(4) + totalEnergy(7)

   ! Compute the total energy.
   totalEnergy(10) = totalEnergy(1) - totalEnergy(4) - totalEnergy(5) + &
         & totalEnergy(7)

   ! Print out the results

   ! Print the orbitals with their associated spins, occupaions, eigenvalues,
   !   kinetic energies, and potential energies.
   write (20,*) "Output data for ",elementName," orbitals"
   write (20,*) "---------------------------"
   write (20,*)
   write (20,10) "nl","s","occ","eigenvalue","kinetic energy","pot energy"
   write (20,*)
   do i = 1, numOrbitals
      write (20,11) orbitalQNn(i),angMomLabels(orbitalQNl(i)+1),orbitalQNs(i),&
            & orbitalCharge(i),eigenValues(i)*0.5_double,&
            & orbitalKE(i)*0.5_double, orbitalPE(i)*0.5_double
   enddo

   10 format (1x,a2,4x,a1,6x,a3,9x,a10,4x,a14,6x,a10)
   11 format (1x,i1,a1,f6.1,f10.4,3f17.8)

   ! Print out total energy data.
   write (20,*)
   write (20,*)
   write (20,*) "Total Energy"
   write (20,*) "------------"
   write (20,*)
   write (20,20) " sum of eigenvalues        =",totalEnergy(1)*0.5_double
   write (20,20) " kinetic energy from ek    =",totalEnergy(2)*0.5_double
   write (20,20) " el-ion interaction energy =",totalEnergy(3)*0.5_double
   write (20,20) " el-el  interaction energy =",totalEnergy(4)*0.5_double
   write (20,20) " vxc    correction         =",totalEnergy(5)*0.5_double
   write (20,20) " virial correction         =",totalEnergy(6)*0.5_double
   write (20,20) " exchange + corr energy    =",totalEnergy(7)*0.5_double
   write (20,20) " kinetic energy from ev    =",totalEnergy(8)*0.5_double
   write (20,20) " potential energy          =",totalEnergy(9)*0.5_double
   write (20,*)  "---------------------------------------------"
   write (20,20) " total energy              =",totalEnergy(10)*0.5_double

   20 format (a27,1x,f18.8)

   ! Print out virial theorem analysis.
   write (20,*)
   write (20,*)
   write (20,*) "virial theorem"
   write (20,*) "--------------"
   write (20,*)
   write (20,30) " kinetic energy  *  2      =",totalEnergy(8)
   write (20,30) " potential energy          =",totalEnergy(9)*0.5_double
   write (20,30) " virial correction         =",totalEnergy(6)*0.5_double
   write (20,*)  "---------------------------------------------"
   write (20,30) " virial sum                =",(2.0_double*totalEnergy(8) + &
         & totalEnergy(9) + totalEnergy(6))*0.5_double

   30 format (a27,1x,f18.8)

end subroutine atomicTotalEnergy

end module AtomicTotalEnergySubs
