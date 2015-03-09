program parSerComp
  use HDF5
  use matSubs
  use hdf5Routines
  use dataStructs
  use compMod

  implicit none

  ! for lindy: valedim = 11725, potdim = 86, numpottypes=4
  ! numpotsites=950, numkpoints=1

  call parseCmdLine()
  print *, "valeDim    : ", systemVars%valeDim
  print *, "numPotSites: ", systemVars%numPotSites
  print *, "numKPoints : ", systemVars%numKPoints
  print *, "potDim     : ", systemVars%potDim
  print *, "numPotTypes: ", systemVars%numPotTypes
  print *, "which      : ", systemVars%which
  valeDim =  systemVars%valeDim
  numPotSites = systemVars%numPotSites
  numKPoints =  systemVars%numKPoints
  potDim =  systemVars%potDim
  numPotTypes =  systemVars%numPotTypes

  if (systemVars%which==1) then
    call elecStatComp()
  else if (systemVars%which==3) then
    call exchCorrComp()
  else if (systemVars%which==5) then
    call integComp()
  else if (systemVars%which==4) then
    call elecStatComp()
    call exchCorrComp()
  else if (systemVars%which==6) then
    call elecStatComp()
    call integComp()
  else if (systemVars%which==8) then
    call exchCorrComp()
    call integComp()
  else if (systemVars%which==9) then
    call elecStatComp()
    call exchCorrComp()
    call integComp()
  endif

end program parSerComp
