program contractWaveFn

   ! Use the necessary program parameters
   use O_Kinds
   use O_Constants

   ! Import the necessary subroutines.
   use ExecutionData
   use AtomData
   use ParseInputSubs
   use NormalizationSubs
   use MatrixElementSubs
   use DiagonalizationSubs
   use PrintResultsSubs

   ! Import necessary modules.
   use ChargeDensityMod

   implicit none

   ! Define the global variables for the calculation.


   ! Obtain the input parameters.
   call parseContractInput

   ! Compute implicit input.
   call implicitInput

   ! Compute the angular normalization coefficients.
   call computeAngNorm

   ! Prepare the hamiltonian and overlap matrix elements.
   call prepMatrixElements

   ! Diagonalize the s, p, d, and f orbitals.
   call diagonalizeSPDF

   ! Renormalize the resultant wave functions.
   call renormWaveFns

   ! Compute and print the electron charge density if requested.
   if (doCharge == 1) then
      call makeRho
      call printRho
   endif

   ! Print the results.
   call printWaveFns

   ! Print numerical wave functions and a POVRay scene file set if requested.
   !   IMPORTANT NOTE:  These function calls must be done in this order so
   !   that the sign of the wave functions can be determined before the POVRay
   !   scenes are generated.
   if (doVis == 1) then
      call printWaveFnsNumerical
      call printPOVRayScenes
   endif

end program contractWaveFn
