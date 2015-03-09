module dataStructs
use HDF5

implicit none
integer, parameter :: single = kind(0.0)
integer, parameter :: double = kind(0.0d0)

integer :: numRayPoints
integer :: maxNumRayPoints

type :: sysData
  integer :: valeDim
  integer :: numPotSites
  integer :: numKPoints
  integer :: potDim
  integer :: numPotTypes
  integer :: which
end type sysData

type(sysData) :: systemVars

contains

subroutine parseCmdLine()
  implicit none

  character*25 :: commandBuffer

  call initSysVars()

  call getarg(1, commandBuffer)
  read (commandBuffer,*) systemVars%valeDim
  call getarg(2, commandBuffer)
  read (commandBuffer,*) systemVars%numPotSites
  call getarg(3, commandBuffer)
  read (commandBuffer,*) systemVars%numKPoints
  call getarg(4, commandBuffer)
  read (commandBuffer,*) systemVars%potDim
  call getarg(5, commandBuffer)
  read (commandBuffer,*) systemVars%numPotTypes
  call getarg(6, commandBuffer)
  read (commandBuffer,*) systemVars%which
end subroutine

subroutine initSysVars()
  systemVars%valeDim     = 0
  systemVars%numPotSites = 0
  systemVars%numKPoints  = 0
  systemVars%potDim      = 0
  systemVars%numPotTypes = 0
end subroutine

end module dataStructs
