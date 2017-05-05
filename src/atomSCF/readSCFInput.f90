module ReadSCFInputSubs

public

contains

subroutine readSCFInput

   ! Import necessary data modules.
   use AtomData
   use TempOrbitalData
   use O_RadialGrid
   use ExecutionData

   ! Import necessary procedure modules. 
   use O_ReadDataSubs

   implicit none

   ! Define local variables
   integer :: i          ! Loop index variable.
   character*8  :: date  ! Date
   character*10 :: time  ! Time
   integer      :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer      :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Open the file needed to read the input.
   open(5,file='atomSCF.dat',status='old',form='formatted')
   readUnit = 5

   ! Open the file needed to write output.
   open(20,file='atomSCF.out',status='new',form='formatted')
   writeUnit = 20

   ! Start parsing the input
   write (writeUnit,*) '***************************************************'
   write (writeUnit,*) '**************  Begin Parsing Input  **************'
   write (writeUnit,*) '***************************************************'

   ! Log the date and time we start.
   call date_and_time(date,time)
   write (writeUnit,100) date,time

   ! Read and echo the title of the calculation.
   call readLabel(readUnit, writeUnit)

   ! Read the name of the element abbreviated from the periodic table.
   call readData(readUnit,writeUnit,len(elementName),elementName,&
         & len('ELEMENT_NAME'),'ELEMENT_NAME')

   ! Read the alpha for the gaussian that defines the nuclear potential.
   call readData(readUnit,writeUnit,nuclearAlpha,len('NUCLEAR_ALPHA'),&
         & 'NUCLEAR_ALPHA')

   ! Read the maximum number of iterations to use.
   call readData(readUnit,writeUnit,maxIteration,len('MAX_ITERATION'),&
         & 'MAX_ITERATION')

   ! Read the level of tolerance to use for testing convergence.
   call readData(readUnit,writeUnit,tolerance,len('CONVG_TOLERANCE'),&
         & 'CONVG_TOLERANCE')

   ! Read the mixing factor to use to mix the input & output electric pot.
   call readData(readUnit,writeUnit,mixingFactor,len('MIXING_FACTOR'),&
         & 'MIXING_FACTOR')

   ! Read the type of exchange correlation to use.
   call readData(readUnit,writeUnit,exchCorrCode,&
         & len('EXCHANGE_CORRELATION_CODE'),'EXCHANGE_CORRELATION_CODE')

   ! Read the flag that says if this is a spin polarized calculation.
   call readData(readUnit,writeUnit,doSpinPol,len('SPIN_POLARIZATION_FLAG'),&
         & 'SPIN_POLARIZATION_FLAG')

   ! Read the flag that says if this is a relativistic calculation.
   call readData(readUnit,writeUnit,doRelativistic,len('RELATIVISTIC_FLAG'),&
         & 'RELATIVISTIC_FLAG')

   ! Read the atomic number (Z value).
   call readData(readUnit,writeUnit,atomicNumber,len('ATOMIC_NUMBER'),&
         & 'ATOMIC_NUMBER')

   ! Read the atomic shell charge.  ! NOT SURE WHAT THIS IS FOR YET.
   call readData(readUnit,writeUnit,shellCharge,len('SHELL_CHARGE'),&
         & 'SHELL_CHARGE')

   ! Read the atomic shell radius.  ! NOT SURE WHAT THIS IS FOR YET.
   call readData(readUnit,writeUnit,shellRadius,len('SHELL_RADIUS'),&
         & 'SHELL_RADIUS')

   ! Read the radial grid parameters.
   call readData(readUnit,writeUnit,radialMaxDist,aaWhatever,bbWhatever,&
         & len('RADIAL_GRID_PARAMETERS'),'RADIAL_GRID_PARAMETERS')

   ! Read the maximum l quantum number present in this system.
   call readData(readUnit,writeUnit,maxQNl,len('MAX_ORB_ANGMOM_QN'),&
         & 'MAX_ORB_ANGMOM_QN')
   ! Increase the value read in by 1 because all loops and array indices must
   !   start with 1 and so this will account for spin 0.  NOTE:  This can cause
   !   confusion and can be tricky.  Watch out!
   maxQNl = maxQNl + 1

   ! Read the number of core orbitals.  This should be the number of core
   !   orbitals in the case that a non spin-polarized calculation is being
   !   done.  It will be adjusted for spin-polarization in the implicit input
   !   subroutines.
   call readData(readUnit,writeUnit,numCoreOrb,len('NUM_CORE_ORBITALS'),&
         & 'NUM_CORE_ORBITALS')

   ! Read the number of valence orbitals.  As with the number of core orbitals
   !   above, this should be the number of valence orbitals in the case that a
   !   non spin-polarized calculation is being done.  It will be adjusted for
   !   spin-polarization in the implicit input subroutines.
   call readData(readUnit,writeUnit,numValeOrb,len('NUM_VALE_ORBITALS'),&
         & 'NUM_VALE_ORBITALS')

   ! Allocate space to hold the orbital quantum numbers.  Note that the number
   !   of orbitals that will actually be used depends on whether or not we are
   !   doing a spin-polarized calculation.  If not, then the total number of 
   !   orbitals is equal to the sum of the core and valence number of orbitals
   !   given in the input above.  If so, then the number of orbitals should be
   !   doubled.
   allocate (orbitalQNn    ((numValeOrb+numCoreOrb)*(1+doSpinPol)))
   allocate (orbitalQNl    ((numValeOrb+numCoreOrb)*(1+doSpinPol)))
   allocate (orbitalQNs    ((numValeOrb+numCoreOrb)*(1+doSpinPol)))
   allocate (orbitalCharge ((numValeOrb+numCoreOrb)*(1+doSpinPol)))

   ! Allocate space to temporarily hold the valence orbital information that is
   !   provided in the input file.  This data will be moved into the above
   !   allocated arrays in the implicit input section.  It is necessary to
   !   only allocate enough space to hold the data given in the input file.
   !   These will be deallocated at the end of the implicit input section that
   !   deals with the valence orbitals.
   allocate (tempValeQNn (numValeOrb))
   allocate (tempValeQNl (numValeOrb))
   allocate (tempValeOrbChargeDn (numValeOrb))
   allocate (tempValeOrbChargeUp (numValeOrb))

   ! Read the definitions of each valence orbital.  These will be adjusted in
   !   the implicit input subroutines.  This is the n and l quantum numbers and
   !   the number of electrons in this orbital's down and up spin states.
   do i = 1, numValeOrb
      read (5,*) tempValeQNn(i), tempValeQNl(i), tempValeOrbChargeDn(i),&
            & tempValeOrbChargeUp(i)
   enddo

   ! Log the date and time we start.
   call date_and_time(date,time)
   write (writeUnit,100) date,time

   ! End parsing the input
   write (writeUnit,*) '***************************************************'
   write (writeUnit,*) '***************  End Parsing Input  ***************'
   write (writeUnit,*) '***************************************************'
   call flush(writeUnit)

   ! Define all the formatted output styles.
   100 format ('Date is: ',a8,' Time is: ',a10) ! Date and time output

end subroutine readSCFInput

end module ReadSCFInputSubs
