module ExchCorrSubs

   public

   contains

   subroutine wigner(rs,vcp,ecp)

      ! Import necessary modules.
      use O_Kinds

      implicit none

      ! Define passed parameters.
      real (kind=double) :: rs
      real (kind=double) :: vcp
      real (kind=double) :: ecp

      ! Compute the potential and energy.
      vcp = -0.88_double*(4.0_double*rs + 3.0_double*7.8_double) / &
            (3.0_double*(rs + 7.8_double)**2.0_double)
      ecp = -0.88_double/(rs+7.8_double)
   end subroutine wigner



   subroutine hidenLundqvist(rs,vcp,ecp)

      ! Import necessary modules.
      use O_Kinds

      implicit none

      ! Define passed parameters.
      real (kind=double) :: rs
      real (kind=double) :: vcp
      real (kind=double) :: ecp

      ! Define local variables.
      real (kind=double) :: x
      real (kind=double) :: aln

      ! Compute parameters to the calculation.
      x   = rs / 21.0_double
      aln = log(1.0_double+1.0_double/x)

      ! Compute the potential and energy.
      vcp = -0.45_double * aln
      if (x > 500.0_double) then
         ecp = ((1.0_double/(6.0_double*x)-0.3_double)/x+0.75_double)/x
      else
         ecp = aln+(x**3.0_double * aln - x*x) + x/2.0_double - &
               & 1.0_double/3.0_double
      endif
      ecp = -0.045_double * ecp
   end subroutine hidenLundqvist



   subroutine gunnarsonLundqvistWilkins(rs,vcp,ecp,vcf,ecf)

      ! Import necessary modules.
      use O_Kinds

      implicit none

      ! Define passed parameters.
      real (kind=double) :: rs
      real (kind=double) :: vcp
      real (kind=double) :: ecp
      real (kind=double) :: vcf
      real (kind=double) :: ecf

      ! Define local variables.
      real (kind=double) :: x
      real (kind=double) :: aln

      ! Compute parameters to the paramagnetic calculation.
      x   = rs / 11.4_double
      aln = log(1.0_double+1.0_double/x)

      ! Compute the paramagnetic potential and energy.
      vcp = -0.0666_double * aln
      if (x > 500.0_double) then
         ecp = ((1.0_double/(6.0_double*x)-0.3_double)/x+0.75_double)/x
      else
         ecp = aln+(x**3.0_double * aln - x*x) + x/2.0_double - &
               & 1.0_double/3.0_double
      endif
      ecp = -0.0666_double * ecp

      ! Compute parameters to the ferromagnetic calculation.
      x   = rs / 15.9_double
      aln = log(1.0_double+1.0_double/x)

      ! Compute the ferromagnetic potential and energy.
      vcf = -0.0406_double * aln
      if (x > 500.0_double) then
         ecf = ((1.0_double/(6.0_double*x)-0.3_double)/x+0.75_double)/x
      else
         ecf = aln+(x**3.0_double * aln - x*x) + x/2.0_double - &
               & 1.0_double/3.0_double
      endif
      ecf = -0.0406_double * ecf
   end subroutine gunnarsonLundqvistWilkins



   subroutine vonBarthHedin(rs,vcp,ecp,vcf,ecf)

      ! Import necessary modules.
      use O_Kinds

      implicit none

      ! Define passed parameters.
      real (kind=double) :: rs
      real (kind=double) :: vcp
      real (kind=double) :: ecp
      real (kind=double) :: vcf
      real (kind=double) :: ecf

      ! Define local variables.
      real (kind=double) :: x
      real (kind=double) :: aln

      ! Compute parameters to the paramagnetic calculation.
      x   = rs / 30.0_double
      aln = log(1.0_double+1.0_double/x)

      ! Compute the paramagnetic potential and energy.
      vcp = -0.0504_double * aln
      if (x > 500.0_double) then
         ecp = ((1.0_double/(6.0_double*x)-0.3_double)/x+0.75_double)/x
      else
         ecp = aln+(x**3.0_double * aln - x*x) + x/2.0_double - &
               & 1.0_double/3.0_double
      endif
      ecp = -0.0504_double * ecp

      ! Compute parameters to the ferromagnetic calculation.
      x   = rs / 75.0_double
      aln = log(1.0_double+1.0_double/x)

      ! Compute the ferromagnetic potential and energy.
      vcf = -0.0254_double * aln
      if (x > 500.0_double) then
         ecf = ((1.0_double/(6.0_double*x)-0.3_double)/x+0.75_double)/x
      else
         ecf = aln+(x**3.0_double * aln - x*x) + x/2.0_double - &
               & 1.0_double/3.0_double
      endif
      ecf = -0.0254_double * ecf
   end subroutine vonBarthHedin


   ! The Perdew-Zunger parameterization is used.
   ! See Phys. Rev. B 23 5075 (1981).
   subroutine ceperleyAlder(rs,vcp,ecp,vcf,ecf)

      ! Import necessary modules.
      use O_Kinds

      implicit none

      ! Define passed parameters.
      real (kind=double) :: rs
      real (kind=double) :: vcp
      real (kind=double) :: ecp
      real (kind=double) :: vcf
      real (kind=double) :: ecf

      ! Define local variables.
      real (kind=double) :: te
      real (kind=double) :: be

      ! Compute the exchange and correlation terms depending on the size of rs.
      if (rs <= 1.0_double) then

         ! Paramagnetic case first.
         te = 1.0_double+(7.0_double/6.0_double)*1.0529_double*sqrt(rs) + &
               & (4.0_double/3.0_double)*0.3334_double*rs
         be = 1.0_double+1.0529_double*sqrt(rs)+0.3334_double*rs
         ecp = -0.2846_double / be
         vcp = -0.2846_double * te/be**2

         ! Ferromagnetic case second.
         te = 1.0_double+(7.0_double/6.0_double)*1.3981_double*sqrt(rs) + &
               & (4.0_double/3.0_double)*0.2611_double*rs
         be = 1.0_double+1.3981_double*sqrt(rs)+0.2611_double*rs
         ecp = -0.1686_double / be
         vcp = -0.1686_double * te/be**2
      else

         ! Paramagnetic case first.
         ecp = 2.0_double*((0.0311_double+0.0020_double*rs)*log(rs) - &
               & 0.048_double - 0.0116_double*rs)
         vcp = 2.0_double*((0.0311_double+2.0_double/3.0_double*0.0020_double*&
               & rs)*log(rs) - (0.048_double+0.0311_double/3.0_double) - &
               & (2.0_double/3.0_double*0.0116_double+0.0020_double / &
               & 3.0_double)*rs)

         ! Ferromagnetic case second.
         ecf = 2.0_double*((0.01555_double+0.0007_double*rs)*log(rs) - &
               & 0.0269_double-0.0048_double*rs)
         vcf = 2.0_double*((0.01555_double+2.0_double/3.0_double* &
               & 0.0007_double*rs)*log(rs) - (0.0269_double+0.01555_double/ &
               & 3.0_double) - (2.0_double/3.0_double*0.0048_double+ &
               & 0.0007_double/3.0_double)*rs)
      endif

   end subroutine ceperleyAlder


   ! GC means gradient corrected.
   subroutine vonBarthHedinGC(rs,vcp,ecp)

      ! Import necessary modules.
      use O_Kinds

      implicit none

      ! Define passed parameters.
      real (kind=double) :: rs
      real (kind=double) :: vcp
      real (kind=double) :: ecp

      ! Define local variables.
      real (kind=double) :: x
      real (kind=double) :: aln

      ! Compute parameters to the calculation.
      x   = rs / 30.0_double
      aln = log(1.0_double+1.0_double/x)

      ! Compute the potential and energy.
      vcp = -0.0504_double * aln
      if (x > 500.0_double) then
         ecp = ((1.0_double/(6.0_double*x)-0.3_double)/x+0.75_double)/x
      else
         ecp = aln+(x**3.0_double * aln - x*x) + x/2.0_double - &
               & 1.0_double/3.0_double
      endif
      ecp = -0.0504_double * ecp
   end subroutine vonBarthHedinGC
end module ExchCorrSubs
