program mainWrapper

   use O_Kinds
   use O_Main

   real (kind=double) :: totalEnergy
   integer :: fromExternal

   fromExternal = 0

   call mainSCF(totalEnergy,fromExternal)

end program mainWrapper
