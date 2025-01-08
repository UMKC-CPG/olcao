module O_MPI

   use MPI_F08

   ! Make sure that no funny variables are implicitly declared.
   implicit none

   ! Define access.
   public

   integer :: mpiRank
   integer :: mpiSize
   integer :: mpierr
   integer :: mpiStatus(MPI_STATUS_SIZE)

   contains

subroutine initMPI

   call MPI_INIT (mpierr)
   call MPI_COMM_RANK (MPI_COMM_WORLD,mpiRank,mpierr)
   call MPI_COMM_SIZE (MPI_COMM_WORLD,mpiSize,mpierr)
   call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

end subroutine initMPI


subroutine stopMPI(errorMessage)
   
   use MPI_F08

   implicit none

   ! Declare passed parameters
   character(len=*) :: errorMessage

   if (mpiRank == 0) then
      write (20,*) errorMessage
   endif
   
   call MPI_FINALIZE(mpierr)

   stop

end subroutine stopMPI


subroutine checkAttributeHDF5(aid,readTarget,hdf5Status)

   ! Use necessary modules.
   use HDF5

   ! Make sure no variables are implicitly decalred.
   implicit none

   ! Define passed dummy parameters.
   integer(hsize_t), intent(in) :: aid ! attribute id
   character(len=*), intent(in) :: readTarget
   integer, intent(out) :: hdf5Status

   ! Define local variables.
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim
   integer :: hdferr

   ! Assume that the HDF5 attribute indicates that the associated dataset
   !   has not yet been computed.
   hdf5Status = 0
   attribIntDims(1) = 1
   if (mpiRank == 0) then
      call h5aread_f(aid,H5T_NATIVE_INTEGER,hdf5Status,&
            & attribIntDims,hdferr)
   endif
   call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
   if (hdferr /= 0) call stopMPI("Failed to read " // readTarget)

   call MPI_BCAST(hdf5Status,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
   if (hdf5Status == 1) then
      if (mpiRank == 0) then
         write(20,*) readTarget, " already exists. Skipping."
         call h5aclose_f(aid,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) call stopMPI(readtarget)
   endif


end subroutine checkAttributeHDF5


subroutine closeMPI

   use MPI_F08

   call MPI_FINALIZE (mpierr)

end subroutine closeMPI

end module O_MPI
