module O_GaussianIntegrals

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   real (kind=double), dimension (3) :: F, AE, BE
!   real (kind=double) :: AT,ATS,EX,EY,EZ,E2,B2,C2,CCC,COEF,FS,XX

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

      subroutine nucPotInteg(g,l1,l2,al1,al2,al3,r1,r2,r3)

      use O_Kinds
      use O_Constants

      implicit none
! c
! c     g(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
! c     9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
! c     14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy
! c
! c     wo(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
! c     10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
! c     18,xxx; 19,yyy; 20,zzz
! c 
      ! define the dummy variables passed to this subroutine.
      real (kind=double), dimension (16,16) :: g
      integer :: l1, l2
      real (kind=double) :: al1
      real (kind=double) :: al2
      real (kind=double) :: al3
      real (kind=double), dimension (3) :: r1
      real (kind=double), dimension (3) :: r2
      real (kind=double), dimension (3) :: r3
      integer :: option

      ! define local variables
      real (kind=double) :: als, ex, a
      real (kind=double), dimension (3) :: b, c
      integer :: nff
      real (kind=double), dimension(20,20) :: dr
      real (kind=double), dimension (3) :: ae2, be2, f2, f3, f4, xf1, &
     & xf2, xf4, xf6, xf7, xf8, xf9, xf10, xf14, xf17, xf31, xf32, &
     & xf33, xf35, xf36, xf37, xf3
      real (kind=double), dimension (3,3) :: abe, aef, bef, aep, bep, ff
      integer :: ndd, i, j, k, l, m1, m2, m3
      real (kind=double) hoa, s0, s2, s4, s6, s8, s10, s12, dp, pd, dd
      real (kind=double) xf11, xf12, xf13, xf15, xf16, xf18, xf19, xf20
      real (kind=double) xf21, xf22, xf23, xf24, xf26, xf27, xf28, xf29
      real (kind=double) xf30, xf34
      real (kind=double), dimension(3) :: f,ae,be
      real (kind=double) :: at,ats,ey,ez,e2,b2,c2,ccc,coef,fs,xx

! c

      do 100 i=1,3
      b(i)=r2(i)-r1(i)
      c(i)=r3(i)-r1(i)
 100  continue
! c
      if(l1.eq.16.or.l2.eq.16) then
         nff=0
      else
         nff=1
      endif
      g(:,:)=0.0_double
      ndd=0
      dr(:,:)=0.0_double
      at=al1+al2+al3                                                       
      ats=at**2                                                         
      hoa=.5d0/at                                                       
      ex=(al2*b(1)+al3*c(1))/at                                           
      ey=(al2*b(2)+al3*c(2))/at                                           
      ez=(al2*b(3)+al3*c(3))/at                                           
      e2=ex*ex+ey*ey+ez*ez                                              
      b2=b(1)*b(1)+b(2)*b(2)+b(3)*b(3)                                  
      c2=c(1)*c(1)+c(2)*c(2)+c(3)*c(3)                                  
      ccc=at*e2-al2*b2-al3*c2                                             
      coef=2.d0*pi*dexp(ccc)/at                                         
      f(1)=(al2*b(1)-(al1+al2)*c(1))/at                                    
      f(2)=(al2*b(2)-(al1+al2)*c(2))/at                                    
      f(3)=(al2*b(3)-(al1+al2)*c(3))/at                                    
      fs=f(1)*f(1)+f(2)*f(2)+f(3)*f(3)                                  
      ae(1)=ex                                                          
      ae(2)=ey                                                          
      ae(3)=ez                                                          
      be(1)=ex-b(1)                                                     
      be(2)=ey-b(2)                                                     
      be(3)=ez-b(3)                                                     
      xx=at*fs
! c                                                          
      call fmtcf(xx,s0,s2,s4,s6,s8,s10,s12)
! c
      do 10 i=1,3                                                       
      ae2(i)=ae(i)*ae(i)                                                
      be2(i)=be(i)*be(i)                                                
      f2(i)=f(i)*f(i)
      f3(i)=f(i)*f2(i)
      do 11 j=1,3                                                       
      abe(i,j)=ae(i)*be(j)                                              
      aef(i,j)=ae(i)*f(j)                                               
      bef(i,j)=be(i)*f(j)                                               
      if (i.eq.j) go to 11                                              
      aep(i,j)=ae(i)*ae(j)                                              
      bep(i,j)=be(i)*be(j)                                              
      ff(i,j)=f(i)*f(j)                                                 
   11 continue                                                          
      f4(i)=f3(i)*f(i)                                                  
      xf1(i)=ae2(i)+hoa                                                 
      xf2(i)=be2(i)+hoa                                                
      xf3(i)=hoa*f(i)                                                   
      xf4(i)=abe(i,i)+hoa                                               
      xf6(i)=be(i)*s0-f(i)*s2                                           
      xf7(i)=be2(i)*f(i)                                                
      xf8(i)=ae2(i)*f(i)                                                
      xf9(i)=ae(i)*s0-f(i)*s2                                           
      xf10(i)=aef(i,i)+bef(i,i)+hoa                                     
      xf14(i)=2.d0*ae(i)+be(i)                                          
      xf17(i)=2.d0*be(i)+ae(i)                                          
      xf31(i)=f(i)*s4-be(i)*s2                                          
      xf32(i)=be(i)*s4-f(i)*s6                                          
      xf33(i)=f(i)*s4-ae(i)*s2                                          
      xf35(i)=f(i)*s8-be(i)*s6                                          
      xf36(i)=ae(i)*s4-f(i)*s6                                          
      xf37(i)=f(i)*s8-ae(i)*s6                                          
   10 continue                                                          
! c this ss                                                               
      dr(1,1)=coef*s0                                             
      do 105 i=1,3                                                      
      xf11=f(i)*s2                                                      
      xf12=f2(i)*s4                                                    
      xf13=f3(i)*s6                                                     
      xf15=3.d0*f(i)                                                    
      xf16=1.5d0*f(i)/at                                                
      xf18=ae2(i)*be(i)                                                 
      xf19=ae2(i)+be2(i)+4.d0*abe(i,i)                                  
      xf20=ae2(i)*f(i)                                                  
      xf21=ae(i)+be(i)                                                  
      xf22=ae(i)*be2(i)                                                 
      xf23=3.d0*f2(i)                                                  
      xf24=ae2(i)*be2(i)                                                
! c this is spx                                                           
      dr(1,i+1)=coef*xf6(i)                                       
! c this is pxs                                                           
      dr(i+1,1)=coef*xf9(i)                                       
      if(ndd.ne.0) go to 1000
! c this is sdxx                                                          
      dr(1,i+4)=coef*(xf2(i)*s0-(2.d0*bef(i,i)+hoa)*s2+xf12)      
! c this is dxxs                                                          
      dr(i+4,1)=coef*(xf1(i)*s0-(2.d0*aef(i,i)+hoa)*s2+xf12)      
 1000 continue
! c this is pxpx                                                          
      dr(i+1,i+1)=coef*(xf4(i)*s0-xf10(i)*s2+xf12)                
      if(ndd.ne.0) go to 1001
! c this is dxxpx                                                         
      dp=(xf18+hoa*xf14(i))*s0                                          
      dp=dp-(xf8(i)+hoa*(xf15+xf14(i))+2.d0*ae(i)*bef(i,i))*s2          
      dp=dp+(f2(i)*xf14(i)+xf16)*s4-xf13                            
      dr(i+4,i+1)=coef*dp                                         
! c this is pxdxx                                                         
      pd=(xf22+hoa*xf17(i))*s0                                          
      pd=pd-(xf7(i)+hoa*(xf15+xf17(i))+2.d0*abe(i,i)*f(i))*s2           
      pd=pd+(f2(i)*xf17(i)+xf16)*s4-xf13                               
      dr(i+1,i+4)=coef*pd                                         
! c this is dxxdxx                                                        
      dd=(xf24+hoa*xf19+.75d0/ats)*s0                                   
      dd=dd-((2.d0*abe(i,i)+3.d0/at)*xf21*f(i)+hoa*xf19+1.5d0/ats)*s2   
      dd=dd+(xf19*f2(i)+xf15*(xf21+f(i))/at+.75d0/ats)*s4              
      dd=dd-(2.d0*f3(i)*xf21+xf23/at)*s6+f4(i)*s8                       
      dr(i+4,i+4)=coef*dd                                         
 1001 continue
  105 continue                                                          
      if(ndd.ne.0) go to 1002
      k=7                                                               
      do 110 i=1,2                                                      
      do 110 j=2,3                                                      
      if (i.eq.j) go to 110                                             
      xf26=ff(i,j)*s4                                                   
      k=k+1                                                             
! c this is sdxy                                                          
      dr(1,k)=coef*(bep(i,j)*s0-(bef(j,i)+bef(i,j))*s2+xf26)      
! c this is dxys                                                          
      dr(k,1)=coef*(aep(i,j)*s0-(aef(j,i)+aef(i,j))*s2+xf26)      
  110 continue                                                          
 1002 continue
      do 115 i=1,3                                                      
      do 115 j=1,3                                                      
      if (i.eq.j) go to 115                                             
      xf27=ff(i,j)*s4                                                   
      xf28=hoa*(f2(i)+f2(j))                                          
      xf29=(aef(i,i)+bef(j,j))/at                                       
      xf30=bef(j,i)+aef(i,j)                                            
! c this is pxpy                                                          
      dr(i+1,j+1)=coef*(abe(i,j)*s0-(aef(i,j)+bef(j,i))*s2 &
     &+xf27)                                                            
      if(ndd.ne.0) go to 1003
! c this is pxdyy                                                         
      pd=xf2(j)*ae(i)*s0-f2(j)*f(i)*s6                                 
      pd=pd-(2.d0*bef(j,j)*ae(i)+hoa*ae(i)+be2(j)*f(i)+xf3(i))*s2       
      pd=pd+(f2(j)*ae(i)+2.d0*be(j)*ff(j,i)+xf3(i))*s4                 
      dr(i+1,j+4)=coef*pd                                         
! c this is dxxpy                                                         
      dp=xf1(i)*be(j)*s0-f2(i)*f(j)*s6                                 
      dp=dp-(2.d0*aef(i,i)*be(j)+hoa*be(j)+ae2(i)*f(j)+xf3(j))*s2       
      dp=dp+(f2(i)*be(j)+2.d0*ae(i)*ff(i,j)+xf3(j))*s4                 
      dr(i+4,j+1)=coef*dp                                         
! c this is dxxdyy                                                        
      dd=xf1(i)*xf2(j)*s0                                               
      dd=dd-(2.d0*abe(i,j)*xf30+hoa*(ae2(i)+be2(j))+xf29+.5d0/ats)*s2   
      dd=dd+(f2(i)*be2(j)+f2(j)*ae2(i)+4.d0*abe(i,j)*ff(i,j) &
     &+xf29+xf28+.25d0/ats)*s4                                          
      dd=dd-(2.d0*ff(i,j)*xf30+xf28)*s6                                 
      dd=dd+f2(i)*f2(j)*s8                                            
      dr(i+4,j+4)=coef*dd                                         
 1003 continue
  115 continue                                                          
      if(ndd.ne.0) go to 1004
      do 130 i=1,3                                                      
      xf34=1.5d0*f(i)/at                                                
      do 130 j=1,2                                                      
      do 130 k=2,3                                                      
      if (j.eq.k) go to 130                                             
      if (i.ne.j.and.i.ne.k) go to 120                                  
      if (j.eq.i) go to 118                                             
      l=j                                                               
      go to 119                                                         
  118 l=k                                                               
  119 continue                                                          
! c this is pxdxy                                                         
      pd=xf6(l)*xf4(i)                                                  
      pd=pd+xf31(l)*xf10(i)                                             
      pd=pd+xf32(l)*f2(i)                                              
      dr(i+1,j+k+5)=coef*pd                                       
! c this is dxxdxy                                                        
      dd=xf6(l)*(ae2(i)*be(i)+hoa*xf14(i))                              
      dd=dd+xf31(l)*(ae2(i)*f(i)+hoa*(3.d0*f(i)+2.d0*ae(i)+be(i))+ &
     &2.d0*abe(i,i)*f(i))                                         
      dd=dd+xf32(l)*(f2(i)*xf14(i)+xf34)                               
      dd=dd+xf35(l)*f3(i)                                               
      dr(i+4,j+k+5)=coef*dd                                       
! c this is dxypx                                                         
      dp=xf9(l)*xf4(i)                                                  
      dp=dp+xf33(l)*xf10(i)                                             
      dp=dp+xf36(l)*f2(i)                                              
      dr(j+k+5,i+1)=coef*dp                                       
! c this is dxydxx                                                        
      dd=xf9(l)*(ae(i)*be2(i)+hoa*xf17(i))                              
      dd=dd+xf33(l)*(xf7(i)+hoa*(3.d0*f(i)+xf17(i))+2.d0*abe(i,i)*f(i)) 
      dd=dd+xf36(l)*(f2(i)*xf17(i)+xf34)                               
      dd=dd+xf37(l)*f3(i)                                               
      dr(j+k+5,i+4)=coef*dd                                       
      go to 130                                                         
  120 continue                                                          
! c this is pxdyz                                                         
      pd=abe(i,j)*be(k)*s0-(ae(i)*(bef(j,k)+bef(k,j))+bef(j,i)* &
     &be(k))*s2                                                         
      pd=pd+(ae(i)*ff(j,k)+be(k)*ff(i,j)+be(j)*ff(i,k))*s4 &
     &-ff(i,j)*f(k)*s6                                                  
      dr(i+1,j+k+5)=coef*pd                                       
! c this is dxxdyz                                                        
      dd=xf6(k)*(xf1(i)*be(j))                                          
      dd=dd+xf31(k)*(2.d0*aef(i,i)*be(j)+hoa*be(j)+ae2(i)*f(j)+xf3(j))  
      dd=dd+xf32(k)*(f2(i)*be(j)+2.d0*ae(i)*ff(i,j)+xf3(j))            
      dd=dd+xf35(k)*(f2(i)*f(j))                                       
      dr(i+4,j+k+5)=coef*dd                                       
! c this is dyzpx                                                         
      dp=aep(j,k)*be(i)*s0-(f(j)*abe(k,i)+aef(j,k)*be(i)+ae(j) &
     &*aef(k,i))*s2                                                     
      dp=dp+(ff(j,k)*be(i)+aef(k,j)*f(i)+ae(j)*ff(k,i))*s4              
      dp=dp-ff(j,k)*f(i)*s6                                             
      dr(j+k+5,i+1)=coef*dp                                       
! c this is dyzdxx                                                        
      dd=xf9(k)*xf2(i)*ae(j)                                            
      dd=dd+xf33(k)*(2.d0*bef(i,i)*ae(j)+hoa*ae(j)+be2(i)*f(j)+xf3(j))  
      dd=dd+xf36(k)*(f2(i)*ae(j)+2.d0*be(i)*ff(i,j)+xf3(j))            
      dd=dd+xf37(k)*f2(i)*f(j)                                         
      dr(j+k+5,i+4)=coef*dd                                       
  130 continue                                                          
      do 146 i=1,2                                                      
      do 146 j=2,3                                                      
      if (i.eq.j) go to 146                                             
      do 145 k=1,2                                                      
      do 145 l=2,3                                                      
      if (k.eq.l) go to 145                                             
      if (i.eq.k.and.j.eq.l) go to 140                                  
      if (i.eq.k) go to 131                                             
      if (i.eq.l) go to 132                                             
      if (j.eq.k) go to 133                                             
      if  (j.eq.l) go to 134                                            
  131 m1=i                                                              
      m2=j                                                              
      m3=l                                                              
      go to 135                                                         
  132 m1=i                                                              
      m2=j                                                              
      m3=k                                                              
      go to 135                                                         
  133 m1=j                                                              
      m2=i                                                              
      m3=l                                                              
      go to 135                                                         
  134 m1=j                                                              
      m2=i                                                              
      m3=k                                                              
      go to 135                                                         
  135 continue                                                          
! c this is dxydyz                                                        
      dd=xf6(m3)*xf4(m1)*ae(m2)                                         
      dd=dd+xf31(m3)*(xf10(m1)*ae(m2)+xf4(m1)*f(m2))                    
      dd=dd+xf32(m3)*(f2(m1)*ae(m2)+xf10(m1)*f(m2))                    
      dd=dd+xf35(m3)*f(m2)*f2(m1)                                      
      dr(i+j+5,k+l+5)=coef*dd                                     
      go to 145                                                         
  140 continue                                                          
! c this is dxydxy                                                        
      dd=xf4(i)*xf4(j)*s0                                               
      dd=dd-(xf10(i)*xf4(j)+xf10(j)*xf4(i))*s2                          
      dd=dd+(xf10(i)*xf10(j)+xf4(i)*f2(j)+xf4(j)*f2(i))*s4            
      dd=dd-(xf10(i)*f2(j)+xf10(j)*f2(i))*s6                          
      dd=dd+f2(i)*f2(j)*s8                                            
      dr(i+j+5,k+l+5)=coef*dd                                     
  145 continue                                                          
  146 continue                                                          
 1004 continue
      if(nff.eq.0) then
        call pot1f2(dr,f,ae,be,at,ats,coef,xx)
        call pot1f1(dr,f,ae,be,at,ats,coef,xx)
      endif

      continue
      g(1,1)=dr(1,1)
      if(l2.eq.1) goto 2009
      g(1,2)=dr(1,2)
      g(1,3)=dr(1,3)
      g(1,4)=dr(1,4)
      if(l2.eq.4) goto 2009
      g(1,5)=dr(1,8)
      g(1,6)=dr(1,9)
      g(1,7)=dr(1,10)
      g(1,8)=dr(1,5)-dr(1,6)
      g(1,9)=2.0*dr(1,7)-dr(1,5)-dr(1,6)
      if(l2.eq.9) goto 2009
      g(1,10)=dr(1,11)
      g(1,11)=dr(1,13)-dr(1,15)
      g(1,12)=dr(1,18)-3.0*dr(1,14)
      g(1,13)=3.0*dr(1,12)-dr(1,19)
      g(1,14)=2.0*dr(1,20)-3.0*dr(1,13)-3.0*dr(1,15)
      g(1,15)=4.0*dr(1,16)-dr(1,18)-dr(1,14)
      g(1,16)=4.0*dr(1,17)-dr(1,12)-dr(1,19)
! c
 2009  continue
      if(l1.eq.1) return
      g(2,1)=dr(2,1)
      if(l2.eq.1) goto 3009
      g(2,2)=dr(2,2)
      g(2,3)=dr(2,3)
      g(2,4)=dr(2,4)
      if(l2.eq.4) goto 3009
      g(2,5)=dr(2,8)
      g(2,6)=dr(2,9)
      g(2,7)=dr(2,10)
      g(2,8)=dr(2,5)-dr(2,6)
      g(2,9)=2.0*dr(2,7)-dr(2,5)-dr(2,6)
      if(l2.eq.9) goto 3009
      g(2,10)=dr(2,11)
      g(2,11)=dr(2,13)-dr(2,15)
      g(2,12)=dr(2,18)-3.0*dr(2,14)
      g(2,13)=3.0*dr(2,12)-dr(2,19)
      g(2,14)=2.0*dr(2,20)-3.0*dr(2,13)-3.0*dr(2,15)
      g(2,15)=4.0*dr(2,16)-dr(2,18)-dr(2,14)
      g(2,16)=4.0*dr(2,17)-dr(2,12)-dr(2,19)
! c
 3009  continue
      g(3,1)=dr(3,1)
      if(l2.eq.1) goto 4009
      g(3,2)=dr(3,2)
      g(3,3)=dr(3,3)
      g(3,4)=dr(3,4)
      if(l2.eq.4) goto 4009
      g(3,5)=dr(3,8)
      g(3,6)=dr(3,9)
      g(3,7)=dr(3,10)
      g(3,8)=dr(3,5)-dr(3,6)
      g(3,9)=2.0*dr(3,7)-dr(3,5)-dr(3,6)
      if(l2.eq.9) goto 4009
      g(3,10)=dr(3,11)
      g(3,11)=dr(3,13)-dr(3,15)
      g(3,12)=dr(3,18)-3.0*dr(3,14)
      g(3,13)=3.0*dr(3,12)-dr(3,19)
      g(3,14)=2.0*dr(3,20)-3.0*dr(3,13)-3.0*dr(3,15)
      g(3,15)=4.0*dr(3,16)-dr(3,18)-dr(3,14)
      g(3,16)=4.0*dr(3,17)-dr(3,12)-dr(3,19)
! c
 4009  continue
      g(4,1)=dr(4,1)
      if(l2.eq.1) goto 5009
      g(4,2)=dr(4,2)
      g(4,3)=dr(4,3)
      g(4,4)=dr(4,4)
      if(l2.eq.4) goto 5009
      g(4,5)=dr(4,8)
      g(4,6)=dr(4,9)
      g(4,7)=dr(4,10)
      g(4,8)=dr(4,5)-dr(4,6)
      g(4,9)=2.0*dr(4,7)-dr(4,5)-dr(4,6)
      if(l2.eq.9) goto 5009
      g(4,10)=dr(4,11)
      g(4,11)=dr(4,13)-dr(4,15)
      g(4,12)=dr(4,18)-3.0*dr(4,14)
      g(4,13)=3.0*dr(4,12)-dr(4,19)
      g(4,14)=2.0*dr(4,20)-3.0*dr(4,13)-3.0*dr(4,15)
      g(4,15)=4.0*dr(4,16)-dr(4,18)-dr(4,14)
      g(4,16)=4.0*dr(4,17)-dr(4,12)-dr(4,19)
! c
 5009  continue
      if(l1.eq.4) return
      g(5,1)=dr(8,1)
      if(l2.eq.1) goto 6009
      g(5,2)=dr(8,2)
      g(5,3)=dr(8,3)
      g(5,4)=dr(8,4)
      if(l2.eq.4) goto 6009
      g(5,5)=dr(8,8)
      g(5,6)=dr(8,9)
      g(5,7)=dr(8,10)
      g(5,8)=dr(8,5)-dr(8,6)
      g(5,9)=2.0*dr(8,7)-dr(8,5)-dr(8,6)
      if(l2.eq.9) goto 6009
      g(5,10)=dr(8,11)
      g(5,11)=dr(8,13)-dr(8,15)
      g(5,12)=dr(8,18)-3.0*dr(8,14)
      g(5,13)=3.0*dr(8,12)-dr(8,19)
      g(5,14)=2.0*dr(8,20)-3.0*dr(8,13)-3.0*dr(8,15)
      g(5,15)=4.0*dr(8,16)-dr(8,18)-dr(8,14)
      g(5,16)=4.0*dr(8,17)-dr(8,12)-dr(8,19)
! c
 6009  continue
      g(6,1)=dr(9,1)
      if(l2.eq.1) goto 7009
      g(6,2)=dr(9,2)
      g(6,3)=dr(9,3)
      g(6,4)=dr(9,4)
      if(l2.eq.4) goto 7009
      g(6,5)=dr(9,8)
      g(6,6)=dr(9,9)
      g(6,7)=dr(9,10)
      g(6,8)=dr(9,5)-dr(9,6)
      g(6,9)=2.0*dr(9,7)-dr(9,5)-dr(9,6)
      if(l2.eq.9) goto 7009
      g(6,10)=dr(9,11)
      g(6,11)=dr(9,13)-dr(9,15)
      g(6,12)=dr(9,18)-3.0*dr(9,14)
      g(6,13)=3.0*dr(9,12)-dr(9,19)
      g(6,14)=2.0*dr(9,20)-3.0*dr(9,13)-3.0*dr(9,15)
      g(6,15)=4.0*dr(9,16)-dr(9,18)-dr(9,14)
      g(6,16)=4.0*dr(9,17)-dr(9,12)-dr(9,19)
! c
 7009  continue
      g(7,1)=dr(10,1)
      if(l2.eq.1) goto 8009
      g(7,2)=dr(10,2)
      g(7,3)=dr(10,3)
      g(7,4)=dr(10,4)
      if(l2.eq.4) goto 8009
      g(7,5)=dr(10,8)
      g(7,6)=dr(10,9)
      g(7,7)=dr(10,10)
      g(7,8)=dr(10,5)-dr(10,6)
      g(7,9)=2.0*dr(10,7)-dr(10,5)-dr(10,6)
      if(l2.eq.9) goto 8009
      g(7,10)=dr(10,11)
      g(7,11)=dr(10,13)-dr(10,15)
      g(7,12)=dr(10,18)-3.0*dr(10,14)
      g(7,13)=3.0*dr(10,12)-dr(10,19)
      g(7,14)=2.0*dr(10,20)-3.0*dr(10,13)-3.0*dr(10,15)
      g(7,15)=4.0*dr(10,16)-dr(10,18)-dr(10,14)
      g(7,16)=4.0*dr(10,17)-dr(10,12)-dr(10,19)
! c
 8009  continue
      g(8,1)=dr(5,1)-dr(6,1)
      if(l2.eq.1) goto 9009
      g(8,2)=dr(5,2)-dr(6,2)
      g(8,3)=dr(5,3)-dr(6,3)
      g(8,4)=dr(5,4)-dr(6,4)
      if(l2.eq.4) goto 9009
      g(8,5)=dr(5,8)-dr(6,8)
      g(8,6)=dr(5,9)-dr(6,9)
      g(8,7)=dr(5,10)-dr(6,10)
      g(8,8)=dr(5,5)-dr(6,5)-dr(5,6)+dr(6,6)
      g(8,9)=2.0*dr(5,7)-2.0*dr(6,7)-dr(5,5) &
     & +dr(6,5)-dr(5,6)+dr(6,6)
      if(l2.eq.9) goto 9009
      g(8,10)=dr(5,11)-dr(6,11)
      g(8,11)=dr(5,13)-dr(6,13)-dr(5,15)+dr(6,15)
      g(8,12)=dr(5,18)-dr(6,18)-3.0*dr(5,14) &
     & +3.0*dr(6,14)
      g(8,13)=3.0*dr(5,12)-3.0*dr(6,12)-dr(5,19)+dr(6,19)
      g(8,14)=2.0*dr(5,20)-2.0*dr(6,20)-3.0*dr(5,13) &
     & +3.0*dr(6,13)-3.0*dr(5,15)+3.0*dr(6,15)
      g(8,15)=4.0*dr(5,16)-4.0*dr(6,16)-dr(5,18)+dr(6,18) &
     & -dr(5,14)+dr(6,14)
      g(8,16)=4.0*dr(5,17)-4.0*dr(6,17)-dr(5,12)+dr(6,12) &
     & -dr(5,19)+dr(6,19)
! c
 9009  continue
      g(9,1)=2.0*dr(7,1)-dr(5,1)-dr(6,1)
      if(l2.eq.1) goto 10009
      g(9,2)=2.0*dr(7,2)-dr(5,2)-dr(6,2)
      g(9,3)=2.0*dr(7,3)-dr(5,3)-dr(6,3)
      g(9,4)=2.0*dr(7,4)-dr(5,4)-dr(6,4)
      if(l2.eq.4) goto 10009
      g(9,5)=2.0*dr(7,8)-dr(5,8)-dr(6,8)
      g(9,6)=2.0*dr(7,9)-dr(5,9)-dr(6,9)
      g(9,7)=2.0*dr(7,10)-dr(5,10)-dr(6,10)
      g(9,8)=2.0*dr(7,5)-dr(5,5)-dr(6,5)-2.0*dr(7,6)+dr(5,6)+dr(6,6)
      g(9,9)=4.0*dr(7,7)-2.0*dr(5,7)-2.0*dr(6,7)-2.0*dr(7,5)+dr(5,5) &
     & +dr(6,5)-2.0*dr(7,6)+dr(5,6)+dr(6,6)
      if(l2.eq.9) goto 10009
      g(9,10)=2.0*dr(7,11)-dr(5,11)-dr(6,11)
      g(9,11)=2.0*dr(7,13)-dr(5,13)-dr(6,13)-2.0*dr(7,15) &
     & +dr(5,15)+dr(6,15)
      g(9,12)=2.0*dr(7,18)-dr(5,18)-dr(6,18)-6.0*dr(7,14)+ &
     & 3.0*dr(5,14)+3.0*dr(6,14)
      g(9,13)=6.0*dr(7,12)-3.0*dr(5,12)-3.0*dr(6,12) &
     & -2.0*dr(7,19)+dr(5,19)+dr(6,19)
      g(9,14)=4.0*dr(7,20)-2.0*dr(5,20)-2.0*dr(6,20)-6.0*dr(7,13)+ &
     & 3.0*dr(5,13)+3.0*dr(6,13)-6.0*dr(7,15)+3.0*dr(5,15)+3.0*dr(6,15)
      g(9,15)=8.0*dr(7,16)-4.0*dr(5,16)-4.0*dr(6,16)-2.0*dr(7,18) &
     & +dr(5,18)+dr(6,18)-2.0*dr(7,14)+dr(5,14)+dr(6,14)
      g(9,16)=8.0*dr(7,17)-4.0*dr(5,17)-4.0*dr(6,17)-2.0*dr(7,12) &
     & +dr(5,12)+dr(6,12)-2.0*dr(7,19)+dr(5,19)+dr(6,19)
! c
 10009 continue
      if(l1.eq.9) return
      g(10,1)=dr(11,1)
      if(l2.eq.1) goto 11009
      g(10,2)=dr(11,2)
      g(10,3)=dr(11,3)
      g(10,4)=dr(11,4)
      if(l2.eq.4) goto 11009
      g(10,5)=dr(11,8)
      g(10,6)=dr(11,9)
      g(10,7)=dr(11,10)
      g(10,8)=dr(11,5)-dr(11,6)
      g(10,9)=2.0*dr(11,7)-dr(11,5)-dr(11,6)
      if(l2.eq.9) goto 11009
      g(10,10)=dr(11,11)
      g(10,11)=dr(11,13)-dr(11,15)
      g(10,12)=dr(11,18)-3.0*dr(11,14)
      g(10,13)=3.0*dr(11,12)-dr(11,19)
      g(10,14)=2.0*dr(11,20)-3.0*dr(11,13)-3.0*dr(11,15)
      g(10,15)=4.0*dr(11,16)-dr(11,18)-dr(11,14)
      g(10,16)=4.0*dr(11,17)-dr(11,12)-dr(11,19)
! c
 11009 continue
      g(11,1)=dr(13,1)-dr(15,1)
      if(l2.eq.1) goto 12009
      g(11,2)=dr(13,2)-dr(15,2)
      g(11,3)=dr(13,3)-dr(15,3)
      g(11,4)=dr(13,4)-dr(15,4)
      if(l2.eq.4) goto 12009
      g(11,5)=dr(13,8)-dr(15,8)
      g(11,6)=dr(13,9)-dr(15,9)
      g(11,7)=dr(13,10)-dr(15,10)
      g(11,8)=dr(13,5)-dr(15,5)-dr(13,6)+dr(15,6)
      g(11,9)=2.0*dr(13,7)-2.0*dr(15,7)-dr(13,5) &
     & +dr(15,5)-dr(13,6)+dr(15,6)
      if(l2.eq.9) goto 12009
      g(11,10)=dr(13,11)-dr(15,11)
      g(11,11)=dr(13,13)-dr(15,13)-dr(13,15)+dr(15,15)
      g(11,12)=dr(13,18)-dr(15,18)-3.0*dr(13,14) &
     & +3.0*dr(15,14)
      g(11,13)=3.0*dr(13,12)-3.0*dr(15,12)- &
     & dr(13,19)+dr(15,19)
      g(11,14)=2.0*dr(13,20)-2.0*dr(15,20)-3.0*dr(13,13) &
     & +3.0*dr(15,13)-3.0*dr(13,15)+3.0*dr(15,15)
      g(11,15)=4.0*dr(13,16)-4.0*dr(15,16)-dr(13,18)+dr(15,18) &
     & -dr(13,14)+dr(15,14)
      g(11,16)=4.0*dr(13,17)-4.0*dr(15,17)-dr(13,12)+dr(15,12) &
     & -dr(13,19)+dr(15,19)
! c
 12009 continue
      g(12,1)=dr(18,1)-3.0*dr(14,1)
      if(l2.eq.1) goto 13009
      g(12,2)=dr(18,2)-3.0*dr(14,2)
      g(12,3)=dr(18,3)-3.0*dr(14,3)
      g(12,4)=dr(18,4)-3.0*dr(14,4)
      if(l2.eq.4) goto 13009
      g(12,5)=dr(18,8)-3.0*dr(14,8)
      g(12,6)=dr(18,9)-3.0*dr(14,9)
      g(12,7)=dr(18,10)-3.0*dr(14,10)
      g(12,8)=dr(18,5)-3.0*dr(14,5) &
     & -dr(18,6)+3.0*dr(14,6)
      g(12,9)=2.0*dr(18,7)-6.0*dr(14,7) &
     & -dr(18,5)+3.0*dr(14,5) &
     & -dr(18,6)+3.0*dr(14,6)
      if(l2.eq.9) goto 13009
      g(12,10)=dr(18,11)-3.0*dr(14,11)
      g(12,11)=dr(18,13)-3.0*dr(14,13) &
     & -dr(18,15)+3.0*dr(14,15)
      g(12,12)=dr(18,18)-3.0*dr(14,18) &
     & -3.0*dr(18,14)+9.0*dr(14,14)
      g(12,13)=3.0*dr(18,12)-9.0*dr(14,12) &
     & -dr(18,19)+3.0*dr(14,19)
      g(12,14)=2.0*dr(18,20)-6.0*dr(14,20)-3.0*dr(18,13) &
     & +9.0*dr(14,13)-3.0*dr(18,15)+9.0*dr(14,15)
      g(12,15)=4.0*dr(18,16)-12.0*dr(14,16)-dr(18,18) &
     & +3.0*dr(14,18)-dr(18,14)+3.0*dr(14,14)
      g(12,16)=4.0*dr(18,17)-12.0*dr(14,17)-dr(18,12) &
     & +3.0*dr(14,12)-dr(18,19)+3.0*dr(14,19)
! c
 13009 continue
      g(13,1)=3.0*dr(12,1)-dr(19,1)
      if(l2.eq.1) goto 14009
      g(13,2)=3.0*dr(12,2)-dr(19,2)
      g(13,3)=3.0*dr(12,3)-dr(19,3)
      g(13,4)=3.0*dr(12,4)-dr(19,4)
      if(l2.eq.4) goto 14009
      g(13,5)=3.0*dr(12,8)-dr(19,8)
      g(13,6)=3.0*dr(12,9)-dr(19,9)
      g(13,7)=3.0*dr(12,10)-dr(19,10)
      g(13,8)=3.0*dr(12,5)-dr(19,5)-3.0*dr(12,6)+dr(19,6)
      g(13,9)=6.0*dr(12,7)-2.0*dr(19,7)-3.0*dr(12,5) &
     & +dr(19,5)-3.0*dr(12,6)+dr(19,6)
      if(l2.eq.9) goto 14009
      g(13,10)=3.0*dr(12,11)-dr(19,11)
      g(13,11)=3.0*dr(12,13)-dr(19,13)-3.0*dr(12,15)+dr(19,15)
      g(13,12)=3.0*dr(12,18)-dr(19,18)-9.0*dr(12,14)+3.0*dr(19,14)
      g(13,13)=9.0*dr(12,12)-3.0*dr(19,12)-3.0*dr(12,19)+dr(19,19)
      g(13,14)=6.0*dr(12,20)-2.0*dr(19,20)-9.0*dr(12,13) &
     & +3.0*dr(19,13)-9.0*dr(12,15)+3.0*dr(19,15)
      g(13,15)=12.0*dr(12,16)-4.0*dr(19,16)-3.0*dr(12,18) &
     & +dr(19,18)-3.0*dr(12,14)+dr(19,14)
      g(13,16)=12.0*dr(12,17)-4.0*dr(19,17)-3.0*dr(12,12) &
     & +dr(19,12)-3.0*dr(12,19)+dr(19,19)
! c
 14009 continue
      g(14,1)=2.0*dr(20,1)-3.0*dr(13,1)-3.0*dr(15,1)
      if(l2.eq.1) goto 15009
      g(14,2)=2.0*dr(20,2)-3.0*dr(13,2)-3.0*dr(15,2)
      g(14,3)=2.0*dr(20,3)-3.0*dr(13,3)-3.0*dr(15,3)
      g(14,4)=2.0*dr(20,4)-3.0*dr(13,4)-3.0*dr(15,4)
      if(l2.eq.4) goto 15009
      g(14,5)=2.0*dr(20,8)-3.0*dr(13,8)-3.0*dr(15,8)
      g(14,6)=2.0*dr(20,9)-3.0*dr(13,9)-3.0*dr(15,9)
      g(14,7)=2.0*dr(20,10)-3.0*dr(13,10)-3.0*dr(15,10)
      g(14,8)=2.0*dr(20,5)-3.0*dr(13,5)-3.0*dr(15,5) &
     & -2.0*dr(20,6)+3.0*dr(13,6)+3.0*dr(15,6)
      g(14,9)=4.0*dr(20,7)-6.0*dr(13,7)-6.0*dr(15,7) &
     & -2.0*dr(20,5)+3.0*dr(13,5)+3.0*dr(15,5) &
     & -2.0*dr(20,6)+3.0*dr(13,6)+3.0*dr(15,6)
      if(l2.eq.9) goto 15009
      g(14,10)=2.0*dr(20,11)-3.0*dr(13,11)-3.0*dr(15,11)
      g(14,11)=2.0*dr(20,13)-3.0*dr(13,13)-3.0*dr(15,13) &
     & -2.0*dr(20,15)+3.0*dr(13,15)+3.0*dr(15,15)
      g(14,12)=2.0*dr(20,18)-3.0*dr(13,18)-3.0*dr(15,18) &
     & -6.0*dr(20,14)+9.0*dr(13,14)+9.0*dr(15,14)
      g(14,13)=6.0*dr(20,12)-9.0*dr(13,12)-9.0*dr(15,12) &
     & -2.0*dr(20,19)+3.0*dr(13,19)+3.0*dr(15,19)
       g(14,14)=4.0*dr(20,20)-6.0*dr(13,20)-6.0*dr(15,20) &
     & -6.0*dr(20,13)+9.0*dr(13,13)+9.0*dr(15,13) &
     & -6.0*dr(20,15)+9.0*dr(13,15)+9.0*dr(15,15)
      g(14,15)=8.0*dr(20,16)-12.0*dr(13,16)-12.0*dr(15,16) &
     & -2.0*dr(20,18)+3.0*dr(13,18)+3.0*dr(15,18) &
     & -2.0*dr(20,14)+3.0*dr(13,14)+3.0*dr(15,14)
      g(14,16)=8.0*dr(20,17)-12.0*dr(13,17)-12.0*dr(15,17) &
     & -2.0*dr(20,12)+3.0*dr(13,12)+3.0*dr(15,12) &
     & -2.0*dr(20,19)+3.0*dr(13,19)+3.0*dr(15,19)
! c
 15009 continue
      g(15,1)=4.0*dr(16,1)-dr(18,1)-dr(14,1)
      if(l2.eq.1) goto 16009
      g(15,2)=4.0*dr(16,2)-dr(18,2)-dr(14,2)
      g(15,3)=4.0*dr(16,3)-dr(18,3)-dr(14,3)
      g(15,4)=4.0*dr(16,4)-dr(18,4)-dr(14,4)
      if(l2.eq.4) goto 16009
      g(15,5)=4.0*dr(16,8)-dr(18,8)-dr(14,8)
      g(15,6)=4.0*dr(16,9)-dr(18,9)-dr(14,9)
      g(15,7)=4.0*dr(16,10)-dr(18,10)-dr(14,10)
      g(15,8)=4.0*dr(16,5)-dr(18,5)-dr(14,5) &
     & -4.0*dr(16,6)+dr(18,6)+dr(14,6)
      g(15,9)=8.0*dr(16,7)-2.0*dr(18,7)-2.0*dr(14,7) &
     & -4.0*dr(16,5)+dr(18,5)+dr(14,5) &
     & -4.0*dr(16,6)+dr(18,6)+dr(14,6)
      if(l2.eq.9) goto 16009
      g(15,10)=4.0*dr(16,11)-dr(18,11)-dr(14,11)
      g(15,11)=4.0*dr(16,13)-dr(18,13)-dr(14,13) &
     & -4.0*dr(16,15)+dr(18,15)+dr(14,15)
      g(15,12)=4.0*dr(16,18)-dr(18,18)-dr(14,18) &
     & -12.0*dr(16,14)+3.0*dr(18,14)+3.0*dr(14,14)
      g(15,13)=12.0*dr(16,12)-3.0*dr(18,12)-3.0*dr(14,12) &
     & -4.0*dr(16,19)+dr(18,19)+dr(14,19)
      g(15,14)=8.0*dr(16,20)-2.0*dr(18,20)-2.0*dr(14,20) &
     & -12.0*dr(16,13)+3.0*dr(18,13)+3.0*dr(14,13) &
     & -12.0*dr(16,15)+3.0*dr(18,15)+3.0*dr(14,15)
      g(15,15)=16.0*dr(16,16)-4.0*dr(18,16)-4.0*dr(14,16) &
     & -4.0*dr(16,18)+dr(18,18)+dr(14,18) &
     & -4.0*dr(16,14)+dr(18,14)+dr(14,14)
      g(15,16)=16.0*dr(16,17)-4.0*dr(18,17)-4.0*dr(14,17) &
     & -4.0*dr(16,12)+dr(18,12)+dr(14,12) &
     & -4.0*dr(16,19)+dr(18,19)+dr(14,19)
! c
 16009 continue
      g(16,1)=4.0*dr(17,1)-dr(12,1)-dr(19,1)
      if(l2.eq.1) goto 17009
      g(16,2)=4.0*dr(17,2)-dr(12,2)-dr(19,2)
      g(16,3)=4.0*dr(17,3)-dr(12,3)-dr(19,3)
      g(16,4)=4.0*dr(17,4)-dr(12,4)-dr(19,4)
      if(l2.eq.4) goto 17009
      g(16,5)=4.0*dr(17,8)-dr(12,8)-dr(19,8)
      g(16,6)=4.0*dr(17,9)-dr(12,9)-dr(19,9)
      g(16,7)=4.0*dr(17,10)-dr(12,10)-dr(19,10)
      g(16,8)=4.0*dr(17,5)-dr(12,5)-dr(19,5) &
     & -4.0*dr(17,6)+dr(12,6)+dr(19,6)
      g(16,9)=8.0*dr(17,7)-2.0*dr(12,7)-2.0*dr(19,7) &
     & -4.0*dr(17,5)+dr(12,5)+dr(19,5) &
     & -4.0*dr(17,6)+dr(12,6)+dr(19,6)
      if(l2.eq.9) goto 17009
      g(16,10)=4.0*dr(17,11)-dr(12,11)-dr(19,11)
      g(16,11)=4.0*dr(17,13)-dr(12,13)-dr(19,13) &
     & -4.0*dr(17,15)+dr(12,15)+dr(19,15)
      g(16,12)=4.0*dr(17,18)-dr(12,18)-dr(19,18) &
     & -12.0*dr(17,14)+3.0*dr(12,14)+3.0*dr(19,14)
      g(16,13)=12.0*dr(17,12)-3.0*dr(12,12)-3.0*dr(19,12) &
     & -4.0*dr(17,19)+dr(12,19)+dr(19,19)
      g(16,14)=8.0*dr(17,20)-2.0*dr(12,20)-2.0*dr(19,20) &
     & -12.0*dr(17,13)+3.0*dr(12,13)+3.0*dr(19,13) &
     & -12.0*dr(17,15)+3.0*dr(12,15)+3.0*dr(19,15)
      g(16,15)=16.0*dr(17,16)-4.0*dr(12,16)-4.0*dr(19,16) &
     & -4.0*dr(17,18)+dr(12,18)+dr(19,18) &
     & -4.0*dr(17,14)+dr(12,14)+dr(19,14)
      g(16,16)=16.0*dr(17,17)-4.0*dr(12,17)-4.0*dr(19,17) &
     & -4.0*dr(17,12)+dr(12,12)+dr(19,12) &
     & -4.0*dr(17,19)+dr(12,19)+dr(19,19)
 17009 continue
      return                                                            
                                                          
      end subroutine nucPotInteg
! c
!*****************************************************
! c
      subroutine pot1f1(dr,f,ae,be,at,ats,coef,xx)

      use O_Kinds
      use O_Constants
!      use genmatrixdatamodule

      implicit none
      real (kind=double), dimension (20,20) :: dr

      integer :: i,j,k,l,m,ms,mt,mm,ma,mb,na,nb
      real (kind=double) :: s0, s2, s4, s6, s8, s10, s12
      real (kind=double) :: aeis, beis, aeic, beic, fis, fic
      real (kind=double) :: aejs, bejs, aejc, bejc, fjs, fjc
      real (kind=double) :: d1, d2, d3, d4, d5, d31, d32, d33, d34, d35
      real (kind=double) :: d331, d332, d333, d334, d335, d6, fifo
      real (kind=double) :: d01, d02, d03, d04, d001, d002, d003, d004
      real (kind=double) :: d11, d12, d13, d14, d21, d22, d23, d24
      real (kind=double) :: dd1, dd2, dd3, dd4, dd5, dd6, dd7
      real (kind=double) :: d41, d42, d43, d44, d45
      real (kind=double) :: d021, d022, d023, d024, d51, d52, d53, d54
      real (kind=double) :: d55, d551, d552, d553, d554, d555
      real (kind=double) :: d61, d62, d63, d64, d661, d662, d663, d664
    
      real (kind=double), dimension(3) :: f,ae,be
      real (kind=double) :: at,ats,coef,xx
! c
      call fmtcf(xx,s0,s2,s4,s6,s8,s10,s12)
! c
! c     this loop has exchange in x-y-z
      do 205 i=1, 3
      aeis=ae(i)*ae(i)
      beis=be(i)*be(i)
      aeic=aeis*ae(i)
      beic=beis*be(i)
      fis=f(i)*f(i)
      fic=fis*f(i)
! c     this fxxxs
      d1=s0*ae(i)*(aeis+1.5/at)
      d2=-s2*(3.*aeis*f(i)+ae(i)*1.5/at+1.5*f(i)/at)
      d3=s4*(3.*ae(i)*fis+1.5*f(i)/at)
      d4=-s6*fic
      dr(i+17,1)=coef*(d1+d2+d3+d4)
! c this  sfxxx
      d1=s0*be(i)*(beis+1.5/at)
      d2=-s2*(3.*beis*f(i)+be(i)*1.5/at+1.5*f(i)/at)
      d3=s4*(3.*be(i)*fis+1.5*f(i)/at)
      d4=-s6*fic
      dr(1,i+17)=coef*(d1+d2+d3+d4)
! c this fxxxpx
      d1=s0*(aeic*be(i)+1.5*ae(i)*(be(i)+ae(i))/at+0.75/ats)
      d2=-s2*(aeic*f(i)+4.5*ae(i)*f(i)/at+3.*aeis*be(i)*f(i) &
     &+1.5*ae(i)*be(i)/at+1.5*be(i)*f(i)/at &
     &+1.5*aeis/at+1.5/ats)
      d3=s4*(3.*aeis*fis+4.5*ae(i)*f(i)/at+3.*ae(i)*be(i)*fis &
     &+1.5*be(i)*f(i)/at+3.*fis/at+3./(4.*ats))
      d4=-s6*(3.*ae(i)*fic+3.*fis/at+fic*be(i))
      d5=s8*f(i)**4
      dr(i+17,i+1)=coef*(d1+d2+d3+d4+d5)
! c this pxfxxx
      d1=s0*(beic*ae(i)+1.5*be(i)*(ae(i)+be(i))/at+0.75/ats)
      d2=-s2*(beic*f(i)+4.5*be(i)*f(i)/at+3.*beis*ae(i)*f(i) &
     &+1.5*be(i)*ae(i)/at+1.5*ae(i)*f(i)/at &
     &+1.5*beis/at+1.5/ats)
      d3=s4*(3.*beis*fis+4.5*be(i)*f(i)/at+3.*be(i)*ae(i)*fis &
     &+1.5*ae(i)*f(i)/at+3.*fis/at+3./(4.*ats))
      d4=-s6*(3.*be(i)*fic+3.*fis/at+fic*ae(i))
      d5=s8*f(i)**4
      dr(i+1,i+17)=coef*(d1+d2+d3+d4+d5)
! c this fxxxdxx
      d31=aeic*(beis+0.5/at)+ae(i)*be(i)*(3.*ae(i)+1.5*be(i)) &
     &/at+9.*ae(i)/(4.*ats)+1.5*be(i)/ats
      d32=aeis*f(i)*(2.*be(i)*ae(i)+3.*beis+4.5/at)+be(i)*f(i) &
     &*(9.*ae(i)+1.5*be(i))/at+aeis*(0.5*ae(i)+3.*be(i))/at &
     &+1.5*ae(i)*beis/at+4.5*ae(i)/ats+3.*be(i)/ats+ &
     &15.*f(i)/(4.*ats)
      d33=aeis*f(i)*(ae(i)*f(i)+6.*be(i)*f(i)+4.5/at) &
     &+3.*ae(i)*be(i)*f(i)*(be(i)*f(i)+3./at)+3.*fis* &
     &(3.*ae(i)+2.*be(i))/at+1.5*beis*f(i)/at+ &
     &(9.*ae(i)/4.+1.5*be(i)+7.5*f(i))/ats
      d34=fic*(3.*aeis+6.*ae(i)*be(i)+beis+5./at) &
     &+fis*(9.*ae(i)/at+6.*be(i)/at)+15.*f(i)/(4.*ats)
      d35=fic*(3.*ae(i)*f(i)+2.*be(i)*f(i)+5./at)
      dr(i+17,i+4)=coef*(s0*d31-s2*d32+s4*d33-s6*d34 &
     &+s8*d35-s10*f(i)**5)
! c this dxxfxxx
      d331=beic*(aeis+0.5/at)+be(i)*ae(i)*(3.*be(i)+1.5*ae(i)) &
     &/at+9.*be(i)/(4.*ats)+1.5*ae(i)/ats
      d332=beis*f(i)*(2.*ae(i)*be(i)+3.*aeis+4.5/at)+ae(i)*f(i) &
     &*(9.*be(i)+1.5*ae(i))/at+beis*(0.5*be(i)+3.*ae(i))/at &
     &+1.5*be(i)*aeis/at+4.5*be(i)/ats+3.*ae(i)/ats+ &
     &15.*f(i)/(4.*ats)
      d333=beis*f(i)*(be(i)*f(i)+6.*ae(i)*f(i)+4.5/at) &
     &+3.*be(i)*ae(i)*f(i)*(ae(i)*f(i)+3./at)+3.*fis* &
     &(3.*be(i)+2.*ae(i))/at+1.5*aeis*f(i)/at+ &
     &(9.*be(i)/4.+1.5*ae(i)+7.5*f(i))/ats
      d334=fic*(3.*beis+6.*be(i)*ae(i)+aeis+5./at) &
     &+fis*(9.*be(i)/at+6.*ae(i)/at)+15.*f(i)/(4.*ats)
      d335=fic*(3.*be(i)*f(i)+2.*ae(i)*f(i)+5./at)
      dr(i+4,i+17)=coef*(s0*d331-s2*d332+s4*d333-s6*d334 &
     &+s8*d335-s10*f(i)**5)
! c this fxxxfxxx
      d1=s0*(be(i)*d31+aeis*(1.5*beis+ae(i)*be(i)+9./(4.* &
     &at))/at+be(i)*(4.5*ae(i)+3.*be(i)/4.)/ats+15./(8.*at**3))
      d2=-s2*(f(i)*d31+be(i)*d32+aeis*(3.*be(i)*f(i)+ &
     &ae(i)*f(i)+ae(i)*be(i))/at+ae(i)*(9.*f(i)+4.5*ae(i))/ats &
     &+3.*ae(i)*be(i)*(be(i)*f(i)+ae(i)*f(i)+ae(i)*be(i)/2.)/at &
     &+be(i)*(6.*f(i)+9.*ae(i)+be(i)*1.5)/ats+45./(8.*at**3))
      d3=s4*(f(i)*d32+be(i)*d33+ae(i)*f(i)*(4.5*ae(i)*f(i) &
     &+aeis+9.*be(i)*f(i)+6.*ae(i)*be(i))/at+3.*be(i)*f(i)* &
     &(be(i)*f(i)/2.+ae(i)*be(i))/at+(18.*ae(i)*f(i)+9.*aeis/4. &
     &+12.*be(i)*f(i)+4.5*ae(i)*be(i)+7.5*fis+3.*beis/4.) &
     &/ats+45./(8.*at**3))
      d4=-s6*(f(i)*d33+be(i)*d34+fis*(6.*ae(i)*f(i)+4.5*aeis &
     &+4.*be(i)*f(i)+9.*ae(i)*be(i)+1.5*beis) &
     &/at+f(i)*(15.*f(i)+9.*ae(i)+6.*be(i))/ats+15./(8.*at**3))
      d5=s8*(f(i)*d34+be(i)*d35+fic*(2.5*f(i)+6.*ae(i) &
     &+4.*be(i))/at+7.5*fis/ats)
      d6=-s10*(f(i)*d35+be(i)*f(i)**5+5.*f(i)**4/(2.*at))
      dr(i+17,i+17)=coef*(d1+d2+d3+d4+d5+d6+s12*f(i)**6)
  205 continue
      l=11
      do 210 i=1,3
      aeis=ae(i)*ae(i)
      beis=be(i)*be(i)
      aeic=aeis*ae(i)
      beic=beis*be(i)
      fis=f(i)*f(i)
      fic=fis*f(i)
      fifo=fis*fis
      do 210 j=1,3
      aejs=ae(j)*ae(j)
      bejs=be(j)*be(j)
      aejc=aejs*ae(j)
      bejc=bejs*be(j)
      fjs=f(j)*f(j)
      fjc=fjs*f(j)
      if(i.eq.j) go to 210
! c this loop has exchange in x-y-z and x-y,y-z,z-x
      l=l+1
! c this fxxys
      d1=s0*ae(j)*(aeis+0.5/at)
      d2=-s2*(2.*ae(i)*ae(j)*f(i)+0.5*ae(j)/at+aeis*f(j) &
     &+0.5*f(j)/at)
      d3=s4*(fis*ae(j)+f(j)*(2.*ae(i)*f(i)+0.5/at))
      dr(l,1)=coef*(d1+d2+d3-s6*fis*f(j))
! c this sfxxy
      d1=s0*be(j)*(beis+0.5/at)
      d2=-s2*(2.*be(i)*be(j)*f(i)+0.5*be(j)/at+beis*f(j) &
     &+0.5*f(j)/at)
      d3=s4*(fis*be(j)+f(j)*(2.*be(i)*f(i)+0.5/at))
      dr(1,l)=coef*(d1+d2+d3-s6*fis*f(j))
! c this fxxypx
      d01=aeis*be(i)*ae(j)+.5*ae(j)*be(i)/at+ae(i)*ae(j)/at
      d02=f(i)*(aeis*ae(j)+1.5*ae(j)/at+2.*be(i)*ae(i)*ae(j))+(aeis &
     &*be(i)+0.5*be(i)/at+ae(i)/at)*f(j)+(0.5*ae(j)*be(i)+ae(i)*ae(j)) &
     &/at
      d03=2.*ae(i)*ae(j)*fis+1.5*ae(j)*f(i)/at+f(i)*f(j)*(aeis+2.* &
     & ae(i)*be(i)+1.5/at)+ae(j)*be(i)*fis+0.5*be(i)*f(j)/at &
     & +ae(i)*f(j)/at
      d04=fis*(ae(j)*f(i)+2.*ae(i)*f(j)+be(i)*f(j))+1.5*f(i)*f(j)/at
      dr(l,i+1)=coef*(s0*d01-s2*d02+s4*d03-s6*d04+s8*fic*f(j))
! c this pxfxxy
      d001=beis*ae(i)*be(j)+.5*be(j)*ae(i)/at+be(i)*be(j)/at
      d002=f(i)*(beis*be(j)+1.5*be(j)/at+2.*ae(i)*be(i)*be(j))+(beis &
     &*ae(i)+0.5*ae(i)/at+be(i)/at)*f(j)+(0.5*be(j)*ae(i)+be(i)*be(j)) &
     &/at
      d003=2.*be(i)*be(j)*fis+1.5*be(j)*f(i)/at+f(i)*f(j)*(beis+2.* &
     & be(i)*ae(i)+1.5/at)+be(j)*ae(i)*fis+0.5*ae(i)*f(j)/at &
     & +be(i)*f(j)/at
      d004=fis*(be(j)*f(i)+2.*be(i)*f(j)+ae(i)*f(j))+1.5*f(i)*f(j)/at
      dr(i+1,l)=coef*(s0*d001-s2*d002+s4*d003-s6*d004+s8*fic*f(j))
! c this fxxypy
      d11=aeis*be(j)*ae(j)+0.5*ae(j)*be(j)/at+0.5*aeis/ &
     &at+0.25/ats
      d12=aeis*ae(j)*f(j)+0.5*ae(j)*f(j)/at+2.*ae(i)*ae(j)* &
     &be(j)*f(i)+0.5*ae(j)*be(j)/at+aeis*be(j)*f(j)+0.5*be(j) &
     &*f(j)/at+ae(i)*f(i)/at+0.5/ats+0.5*aeis/at
      d13=2.*ae(i)*ae(j)*f(i)*f(j)+0.5*f(j)*ae(j)/at+aeis*fjs &
     &+0.5*fjs/at+ae(j)*be(j)*fis+2.*ae(i)*be(j)* &
     &f(i)*f(j)+0.5*be(j)*f(j)/at+0.5*fis/at+ae(i)*f(i)/at+0.25/ats
      d14=ae(j)*fis*f(j)+2.*f(i)*ae(i)*fjs+0.5*fjs &
     & /at+be(j)*fis*f(j)+0.5*fis/at
      dr(l,j+1)=coef*(s0*d11-s2*d12+s4*d13-s6*d14+s8*fjs*fis)
! c this pyfxxy
      d11=beis*ae(j)*be(j)+0.5*be(j)*ae(j)/at+0.5*beis/ &
     &at+0.25/ats
      d12=beis*be(j)*f(j)+0.5*be(j)*f(j)/at+2.*be(i)*be(j)* &
     &ae(j)*f(i)+0.5*be(j)*ae(j)/at+beis*ae(j)*f(j)+0.5*ae(j) &
     & *f(j)/at+be(i)*f(i)/at+0.5/ats+0.5*beis/at
      d13=2.*be(i)*be(j)*f(i)*f(j)+0.5*f(j)*be(j)/at+beis*fjs &
     & +0.5*fjs/at+be(j)*ae(j)*fis+2.*be(i)*ae(j)* &
     & f(i)*f(j)+0.5*ae(j)*f(j)/at+0.5*fis/at+be(i)*f(i)/at+0.25/ats
      d14=be(j)*fis*f(j)+2.*f(i)*be(i)*fjs+0.5*fjs &
     & /at+ae(j)*fis*f(j)+0.5*fis/at
      dr(j+1,l)=coef*(s0*d11-s2*d12+s4*d13-s6*d14+s8*fjs*fis)
! c this fxxydxx
      d21=aeis*beis+0.5*aeis/at+2.*ae(i)*be(i)/at+0.5 &
     &*beis/at+0.75/ats
      d22=(2.*ae(i)*be(i)*f(i)+3.*f(i)/at)*(ae(i)+be(i))+0.5* &
     &(aeis+4.*be(i)*ae(i)+beis)/at+1.5/ats
      d23=fis*(aeis+4.*ae(i)*be(i)+beis)+3.*f(i)* &
     &(ae(i)+be(i)+f(i))/at+0.75/ats
      d24=2.*fic*(ae(i)+be(i))+3.*fis/at
      dr(l,i+4)=coef*(s0*ae(j)*d21-s2*(f(j)*d21+ae(j)*d22) &
     &+s4*(f(j)*d22+ae(j)*d23)-s6*(f(j)*d23+ae(j)*d24)+s8*(f(j) &
     &*d24+ae(j)*fifo)-s10*f(j)*fifo)
! c this dxxfxxy
      dr(i+4,l)=coef*(s0*be(j)*d21-s2*(f(j)*d21+be(j)*d22) &
     &+s4*(f(j)*d22+be(j)*d23)-s6*(f(j)*d23+be(j)*d24)+s8*(f(j) &
     &*d24+be(j)*fifo)-s10*f(j)*fifo)
! c this fxxydxy
      if(i.ne.3.and.j.ne.3)m=8
      if(i.ne.1.and.j.ne.1)m=10
      if(i.ne.2.and.j.ne.2)m=9
      dd1=be(j)*d01+0.5*(aeis*be(i)+0.5*be(i)/at+ae(i)/at)/at
      dd2=f(j)*d01+be(j)*d02+0.5*(aeis*f(i)+1.5*f(i)/at+2.*ae(i)/at+ &
     &2.*be(i)*ae(i)*f(i)+1.*be(i)/at+aeis*be(i))/at
      dd3=f(j)*d02+be(j)*d03+0.5*(2.*ae(i)*fis+3.0*f(i)/at &
     & +be(i)*fis+aeis*f(i)+2.*ae(i)*be(i)*f(i)+0.5*be(i)/at+ae(i) &
     & /at)/at
      dd4=f(j)*d03+be(j)*d04+0.5*fic/at+ae(i)*fis/at+0.75*f(i)/ats &
     &  +0.5*be(i)*fis/at
      dd5=f(j)*d04+be(j)*fic*f(j)+0.5*fic/at
      dd6=fic*fjs
      dr(l,m)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6)
! c this dxyfxxy
      dd1=ae(j)*d001+0.5*(beis*ae(i)+0.5*ae(i)/at+be(i)/at)/at
      dd2=f(j)*d001+ae(j)*d002+0.5*(beis*f(i)+1.5*f(i)/at+2.*be(i)/at+ &
     &2.*ae(i)*be(i)*f(i)+1.*ae(i)/at+beis*ae(i))/at
      dd3=f(j)*d002+ae(j)*d003+0.5*(2.*be(i)*fis+3.0*f(i)/at &
     & +ae(i)*fis+beis*f(i)+2.*be(i)*ae(i)*f(i)+0.5*ae(i)/at+be(i) &
     & /at)/at
      dd4=f(j)*d003+ae(j)*d004+0.5*fic/at+be(i)*fis/at+0.75*f(i)/ats &
     &  +0.5*ae(i)*fis/at
      dd5=f(j)*d004+ae(j)*fic*f(j)+0.5*fic/at
      dd6=fic*fjs
      dr(m,l)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6)
! c this fxxydyy
      dd1=(aeis+0.5/at)*(bejs+0.5/at)
      dd2=2.*ae(i)*be(j)*(be(j)*f(i)+ae(i)*f(j)) +0.5*(aeis &
     &+bejs)/at+(ae(i)*f(i)+be(j)*f(j))/at+0.5/ats
      dd3=fis*bejs+fjs*aeis+4.*ae(i)*be(j)*f(i)* &
     &f(j)+(ae(i)*f(i)+be(j)*f(j))/at+0.5*(fis+fjs)/at+0.25 &
     &/ats
      dd4=2.*f(i)*f(j)*(be(j)*f(i)+ae(i)*f(j))+0.5*(fis+fjs &
     &)/at
      d41=ae(j)*dd1+be(j)*(aeis+0.5/at)/at
      d42=f(j)*dd1+ae(j)*dd2+ae(i)*(be(j)*f(i)+ae(i)*f(j))/at+ae(i) &
     &*be(j)*f(i)/at+aeis*be(j)/at+be(j)/ats+0.5*f(j)/ats
      d43=ae(j)*dd3+f(j)*dd2+be(j)*fis/at+aeis*f(j)/at &
     &+2.*ae(i)*f(i)*f(j)/at+2.*ae(i)*be(j)*f(i)/at+0.5*be(j)/ &
     &ats+1.*f(j)/ats
      d44=ae(j)*dd4+f(j)*dd3+f(i)*(be(j)*f(i)+ae(i)*f(j))/at+fis &
     &*f(j)/at+ae(i)*f(i)*f(j)/at+0.5*f(j)/ats
      d45=ae(j)*fis*fjs+f(j)*dd4+f(j)*fis/at
      dr(l,j+4)=coef*(s0*d41-s2*d42+s4*d43-s6*d44+s8*d45-s10* &
     &  fjc*fis)
! c this dyyfxxy
      dd1=(beis+0.5/at)*(aejs+0.5/at)
      dd2=2.*be(i)*ae(j)*(ae(j)*f(i)+be(i)*f(j)) +0.5*(beis &
     &+aejs)/at+(be(i)*f(i)+ae(j)*f(j))/at+0.5/ats
      dd3=fis*aejs+fjs*beis+4.*be(i)*ae(j)*f(i)* &
     &f(j)+(be(i)*f(i)+ae(j)*f(j))/at+0.5*(fis+fjs)/at+0.25 &
     &/ats
      dd4=2.*f(i)*f(j)*(ae(j)*f(i)+be(i)*f(j))+0.5*(fis+fjs &
     &)/at
      d41=be(j)*dd1+ae(j)*(beis+0.5/at)/at
      d42=f(j)*dd1+be(j)*dd2+be(i)*(ae(j)*f(i)+be(i)*f(j))/at+be(i) &
     &*ae(j)*f(i)/at+beis*ae(j)/at+ae(j)/ats+0.5*f(j)/ats
      d43=be(j)*dd3+f(j)*dd2+ae(j)*fis/at+beis*f(j)/at &
     &+2.*be(i)*f(i)*f(j)/at+2.*be(i)*ae(j)*f(i)/at+0.5*ae(j)/ &
     &ats+1.*f(j)/ats
      d44=be(j)*dd4+f(j)*dd3+f(i)*(ae(j)*f(i)+be(i)*f(j))/at+fis &
     &*f(j)/at+be(i)*f(i)*f(j)/at+0.5*f(j)/ats
      d45=be(j)*fis*fjs+f(j)*dd4+f(j)*fis/at
      dr(j+4,l)=coef*(s0*d41-s2*d42+s4*d43-s6*d44+s8*d45-s10* &
     &fjc*fis)
! c this fxxyfxxy
      dd1=be(j)*ae(j)*d21+0.5*d21/at
      dd2=f(j)*ae(j)*d21+be(j)*(f(j)*d21+ae(j)*d22)+0.5*(d21+d22) &
     &/at
      dd3=f(j)*(f(j)*d21+ae(j)*d22)+be(j)*(f(j)*d22+ae(j)*d23) &
     &+0.5*(d22+d23)/at
      dd4=be(j)*(f(j)*d23+ae(j)*d24)+f(j)*(f(j)*d22+ae(j)*d23) &
     &+0.5*(d23+d24)/at
      dd5=f(j)*(f(j)*d23+ae(j)*d24)+be(j)*(f(j)*d24+ae(j)*fifo) &
     &+0.5*(d24+fifo)/at
      dd6=f(j)*(f(j)*d24+ae(j)*fifo)+be(j)*f(j)*fifo+0.5 &
     &*fifo/at
      dr(l,l)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6 &
     &+s12*fjs*fifo)
! c this fxxxpy
      d021=be(j)*(aeic+1.5*ae(i)/at)
      d022=aeic*f(j)+1.5*ae(i)*f(j)/at+3.*aeis*be(j)*f(i) &
     &+1.5*ae(i)*be(j)/at+1.5*be(j)*f(i)/at
      d023=3.*aeis*f(i)*f(j)+1.5*ae(i)*f(j)/at+1.5*f(j)*f(i)/at &
     &+3.*ae(i)*be(j)*fis+1.5*f(i)*be(j)/at
      d024=3.*ae(i)*fis*f(j)+1.5*f(i)*f(j)/at+fic*be(j)
      dr(i+17,j+1)=coef*(s0*d021-s2*d022+s4*d023-s6*d024+s8*f(j) &
     &*fic)
! c this pyfxxx
      d021=ae(j)*(beic+1.5*be(i)/at)
      d022=beic*f(j)+1.5*be(i)*f(j)/at+3.*beis*ae(j)*f(i) &
     &+1.5*be(i)*ae(j)/at+1.5*ae(j)*f(i)/at
      d023=3.*beis*f(i)*f(j)+1.5*be(i)*f(j)/at+1.5*f(j)*f(i)/at &
     &+3.*be(i)*ae(j)*fis+1.5*f(i)*ae(j)/at
      d024=3.*be(i)*fis*f(j)+1.5*f(i)*f(j)/at+fic*ae(j)
      dr(j+1,i+17)=coef*(s0*d021-s2*d022+s4*d023-s6*d024+s8*f(j) &
     &*fic)
! c this fxxxdyy
      d51=aeic*bejs+0.5*(ae(i)*bejs+aeic+1.5*ae(i) &
     &/at+2.*ae(i)*bejs)/at
      d52=3.*aeis*bejs*f(i)+2.*aeic*be(j)*f(j)+0.5* &
     &aeic/at+1.5*bejs*ae(i)/at+1.5*aeis*f(i)/at+3.* &
     &ae(i)*be(j)*f(j)/at+1.5*ae(i)/ats+1.5*bejs*f(i)/at+0.75* &
     &f(i)/ats
      d53=3.*ae(i)*bejs*fis+6.*aeis*be(j)*f(i)*f(j) &
     &+1.5*aeis*f(i)/at+1.5*bejs*f(i)/at+1.5*ae(i)*fis &
     &/at+3.*be(j)*f(i)*f(j)/at+1.5*f(i)/ats+aeic*fjs &
     &+3.*ae(i)*be(j)*f(j)/at+1.5*ae(i)*fjs/at+0.75*ae(i)/ats
      d54=6.*fis*f(j)*ae(i)*be(j)+3.*f(i)*fjs*aeis+1.5 &
     &*ae(i)*fis/at+1.5*ae(i)*fjs/at+bejs*fic+3.* &
     &be(j)*f(i)*f(j)/at+0.5*fic/at+1.5*f(i)*fjs/at &
     &+0.75*f(i)/ats
      d55=3.*ae(i)*fis*fjs+2.*be(j)*fic*f(j)+1.5*f(i) &
     &*fjs/at+0.5*fic/at
      dr(i+17,j+4)=coef*(s0*d51-s2*d52+s4*d53-s6*d54+s8*d55-s10 &
     &*fic*fjs)
! c this dyyfxxx
      d551=beic*aejs+0.5*(be(i)*aejs+beic+1.5*be(i) &
     &/at+2.*be(i)*aejs)/at
      d552=3.*beis*aejs*f(i)+2.*beic*ae(j)*f(j)+0.5* &
     &beic/at+1.5*aejs*be(i)/at+1.5*beis*f(i)/at+3.* &
     &be(i)*ae(j)*f(j)/at+1.5*be(i)/ats+1.5*aejs*f(i)/at+0.75* &
     &f(i)/ats
      d553=3.*be(i)*aejs*fis+6.*beis*ae(j)*f(i)*f(j) &
     &+1.5*beis*f(i)/at+1.5*aejs*f(i)/at+1.5*be(i)*fis &
     &/at+3.*ae(j)*f(i)*f(j)/at+1.5*f(i)/ats+beic*fjs &
     &+3.*be(i)*ae(j)*f(j)/at+1.5*be(i)*fjs/at+0.75*be(i)/ats
      d554=6.*fis*f(j)*be(i)*ae(j)+3.*f(i)*fjs*beis+1.5 &
     &*be(i)*fis/at+1.5*be(i)*fjs/at+aejs*fic+3.* &
     &ae(j)*f(i)*f(j)/at+0.5*fic/at+1.5*f(i)*fjs/at &
     &+0.75*f(i)/ats
      d555=3.*be(i)*fis*fjs+2.*ae(j)*fic*f(j)+1.5*f(i) &
     &*fjs/at+0.5*fic/at
      dr(j+4,i+17)=coef*(s0*d551-s2*d552+s4*d553-s6*d554+s8*d555-s10 &
     &*fic*fjs)
! c this fxxxdxy
      if(i.ne.3.and.j.ne.3) ms=8
      if(i.ne.2.and.j.ne.2) ms=9
      if(i.ne.1.and.j.ne.1) ms=10
      d61=aeic*be(i)+1.5*ae(i)*be(i)/at+1.5*aeis/at+0.75/ats
      d62=aeic*f(i)+4.5*ae(i)*f(i)/at+3.*aeis*be(i)*f(i)+1.5*ae(i) &
     &*be(i)/at+1.5*be(i)*f(i)/at+1.5*aeis/at+1.5/ats
      d63=3.*aeis*fis+4.5*ae(i)*f(i)/at+3.*ae(i)*be(i)* &
     &fis+1.5*f(i)*be(i)/at+3.*fis/at+0.75/ats
      d64=3.*ae(i)*fic+3.*fis/at+fic*be(i)
      dr(i+17,ms)=coef*(s0*d61*be(j)-s2*(f(j)*d61+be(j)*d62) &
     &+s4*(f(j)*d62+be(j)*d63)-s6*(f(j)*d63+be(j)*d64)+s8* &
     &(f(j)*d64+be(j)*fifo)-s10*f(j)*fifo)
! c this dxyfxxx
      d661=beic*ae(i)+1.5*be(i)*ae(i)/at+1.5*beis/at+0.75/ats
      d662=beic*f(i)+4.5*be(i)*f(i)/at+3.*beis*ae(i)*f(i)+1.5*be(i) &
     &*ae(i)/at+1.5*ae(i)*f(i)/at+1.5*beis/at+1.5/ats
      d663=3.*beis*fis+4.5*be(i)*f(i)/at+3.*be(i)*ae(i)* &
     &fis+1.5*f(i)*ae(i)/at+3.*fis/at+0.75/ats
      d664=3.*be(i)*fic+3.*fis/at+fic*ae(i)
      dr(ms,i+17)=coef*(s0*d661*ae(j)-s2*(f(j)*d661+ae(j)*d662) &
     &+s4*(f(j)*d662+ae(j)*d663)-s6*(f(j)*d663+ae(j)*d664)+s8* &
     &(f(j)*d664+ae(j)*fifo)-s10*f(j)*fifo)
! c this fxxxfxxy
      if(i.eq.1.and.j.eq.2) mt=12
      if(i.eq.1.and.j.eq.3) mt=13
      if(i.eq.2.and.j.eq.1) mt=14
      if(i.eq.2.and.j.eq.3) mt=15
      if(i.eq.3.and.j.eq.1) mt=16
      if(i.eq.3.and.j.eq.2) mt=17
      d31=aeic*(beis+0.5/at)+ae(i)*be(i)*(3.*ae(i)+1.5*be(i)) &
     &/at+9.*ae(i)/(4.*ats)+1.5*be(i)/ats
      d32=aeis*f(i)*(2.*be(i)*ae(i)+3.*beis+4.5/at)+be(i)*f(i) &
     &*(9.*ae(i)+1.5*be(i))/at+aeis*(0.5*ae(i)+3.*be(i))/at &
     &+1.5*ae(i)*beis/at+4.5*ae(i)/ats+3.*be(i)/ats+ &
     & 15.*f(i)/(4.*ats)
      d33=aeis*f(i)*(ae(i)*f(i)+6.*be(i)*f(i)+4.5/at) &
     &+3.*ae(i)*be(i)*f(i)*(be(i)*f(i)+3./at)+3.*fis* &
     &(3.*ae(i)+2.*be(i))/at+1.5*beis*f(i)/at+ &
     &(9.*ae(i)/4.+1.5*be(i)+7.5*f(i))/ats
      d34=fic*(3.*aeis+6.*ae(i)*be(i)+beis+5./at) &
     &+fis*(9.*ae(i)/at+6.*be(i)/at)+15.*f(i)/(4.*ats)
      d35=fic*(3.*ae(i)*f(i)+2.*be(i)*f(i)+5./at)
      dr(i+17,l)=coef*(s0*be(j)*d31-s2*(f(j)*d31+be(j)*d32) &
     &+s4*(f(j)*d32+be(j)*d33)-s6*(f(j)*d33+be(j)*d34)+s8*(f(j) &
     &*d34+be(j)*d35)-s10*(f(j)*d35+be(j)*f(i)**5)+s12*f(j)*f(i)**5)
! c this fxxyfxxx
      d31=beic*(aeis+0.5/at)+be(i)*ae(i)*(3.*be(i)+1.5*ae(i)) &
     &/at+9.*be(i)/(4.*ats)+1.5*ae(i)/ats
      d32=beis*f(i)*(2.*ae(i)*be(i)+3.*aeis+4.5/at)+ae(i)*f(i) &
     &*(9.*be(i)+1.5*ae(i))/at+beis*(0.5*be(i)+3.*ae(i))/at &
     &+1.5*be(i)*aeis/at+4.5*be(i)/ats+3.*ae(i)/ats+ &
     &15.*f(i)/(4.*ats)
      d33=beis*f(i)*(be(i)*f(i)+6.*ae(i)*f(i)+4.5/at) &
     &+3.*be(i)*ae(i)*f(i)*(ae(i)*f(i)+3./at)+3.*fis* &
     &(3.*be(i)+2.*ae(i))/at+1.5*aeis*f(i)/at+ &
     &(9.*be(i)/4.+1.5*ae(i)+7.5*f(i))/ats
      d34=fic*(3.*beis+6.*be(i)*ae(i)+aeis+5./at) &
     &+fis*(9.*be(i)/at+6.*ae(i)/at)+15.*f(i)/(4.*ats)
      d35=fic*(3.*be(i)*f(i)+2.*ae(i)*f(i)+5./at)
      dr(l,i+17)=coef*(s0*ae(j)*d31-s2*(f(j)*d31+ae(j)*d32) &
     &+s4*(f(j)*d32+ae(j)*d33)-s6*(f(j)*d33+ae(j)*d34)+s8*(f(j) &
     &*d34+ae(j)*d35)-s10*(f(j)*d35+ae(j)*f(i)**5)+s12*f(j)*f(i)**5)
! c this fxxx fxyy
      if(i.eq.1.and.j.eq.2) mm=14
      if(i.eq.1.and.j.eq.3) mm=16
      if(i.eq.2.and.j.eq.1) mm=12
      if(i.eq.2.and.j.eq.3) mm=17
      if(i.eq.3.and.j.eq.1) mm=13
      if(i.eq.3.and.j.eq.2) mm=15
      dd1=be(i)*d51+.5*(3.*aeis*bejs+.5*(bejs+3.*aeis &
     &  +1.5/at+2.*bejs)/at)/at
      dd2=be(i)*d52+f(i)*d51+.5*(6.*ae(i)*bejs*f(i)+3.*aeis &
     &*bejs+6.*aeis*be(j)*f(j)+1.5*aeis/at+1.5*bejs &
     &/at+1.5*aeis/at+3.*f(i)*ae(i)/at+3.*be(j)*f(j)/at+1.5/ats &
     & +1.5*bejs/at+.75/ats)/at
      dd3=f(i)*d52+be(i)*d53+.5*(3.*bejs*fis+6.*ae(i)*f(i) &
     & *bejs+12.*ae(i)*be(j)*f(i)*f(j)+6.*aeis*be(j)*f(j) &
     & +1.5*aeis/at+3.*ae(i)*f(i)/at+1.5*bejs/at+1.5*fis/at+3.* &
     & ae(i)*f(i)/at &
     & +3.*be(j)*f(j)/at+1.5/ats+3.*aeis*fjs+3.*be(j)/at* &
     & f(j)+1.5*fjs/at+.75/ats)/at
      dd4=f(i)*d53+be(i)*d54+.5*(6.*fis*f(j)*be(j)+12.*f(i)*f(j) &
     & *ae(i)*be(j)+3.*fjs*aeis+6.*f(i)*fjs*ae(i) &
     & +1.5*fis/at+3.*ae(i)*f(i)/at+1.5*fjs/at+3.*bejs &
     & *fis+3.*be(j)*f(j)/at+1.5*fis/at+1.5*fjs/at+.75/ &
     & ats)/at
      dd5=f(i)*d54+be(i)*d55+.5*(3.*fis*fjs+6.*ae(i)*f(i) &
     & *fjs+6.*be(j)*fis*f(j)+1.5*fjs/at+1.5*fis &
     & /at)/at
      dd6=f(i)*d55+be(i)*fic*fjs+fis*fjs*1.5/at
!****************************************************************
      dd7=fifo*fjs
      dr(i+17,mm)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10 &
     &  *dd6+s12*dd7)
! c   this fxyy fxxx
      dd1=ae(i)*d551+.5*(3.*beis*aejs+.5*(aejs+3.*beis &
     &  +1.5/at+2.*aejs)/at)/at
      dd2=ae(i)*d552+f(i)*d551+.5*(6.*be(i)*aejs*f(i)+3.*beis &
     &*aejs+6.*beis*ae(j)*f(j)+1.5*beis/at+1.5*ae(j) &
     & **2/at+1.5*beis/at+3.*f(i)*be(i)/at+3.*ae(j)*f(j)/at+1.5/ats &
     & +1.5*aejs/at+.75/ats)/at
      dd3=f(i)*d552+ae(i)*d553+.5*(3.*aejs*fis+6.*be(i)*f(i) &
     & *aejs+12.*be(i)*ae(j)*f(i)*f(j)+6.*beis*ae(j)*f(j) &
     & +1.5*beis/at+3.*be(i)*f(i)/at+1.5*aejs/at+1.5*fis/at+3.* &
     & be(i)*f(i)/at &
     & +3.*ae(j)*f(j)/at+1.5/ats+3.*beis*fjs+3.*ae(j)/at* &
     & f(j)+1.5*fjs/at+.75/ats)/at
      dd4=f(i)*d553+ae(i)*d554+.5*(6.*fis*f(j)*ae(j)+12.*f(i)*f(j) &
     & *be(i)*ae(j)+3.*fjs*beis+6.*f(i)*fjs*be(i) &
     & +1.5*fis/at+3.*be(i)*f(i)/at+1.5*fjs/at+3.*aejs &
     & *fis+3.*ae(j)*f(j)/at+1.5*fis/at+1.5*fjs/at+.75/ &
     & ats)/at
      dd5=f(i)*d554+ae(i)*d555+.5*(3.*fis*fjs+6.*be(i)*f(i) &
     & *fjs+6.*ae(j)*fis*f(j)+1.5*fjs/at+1.5*fis &
     & /at)/at
      dd6=f(i)*d555+ae(i)*fic*fjs+fis*fjs*1.5/at
!****************************************************************
      dd7=fifo*fjs
      dr(mm,i+17)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10 &
     &  *dd6+s12*dd7)
  210 continue
! c this loop only has 3-symmetry elements for x,y,z
      do 215 i=1,2
      aeis=ae(i)*ae(i)
      beis=be(i)*be(i)
      aeic=aeis*ae(i)
      beic=beis*be(i)
      fis=f(i)*f(i)
      fic=fis*f(i)
      do 215 j=2,3
      aejs=ae(j)*ae(j)
      bejs=be(j)*be(j)
      aejc=aejs*ae(j)
      bejc=bejs*be(j)
      fjs=f(j)*f(j)
      fjc=fjs*f(j)
      if(i.eq.j) go to 215
      k=i+j-2
      go to(21,22,23) ,k
   21 ma=12
      mb=14
      na=18
      nb=19
      go to 141
   22 ma=13
      mb=16
      na=18
      nb=20
      go to 141
   23 ma=15
      mb=17
      na=19
      nb=20
      go to 141
  141 continue
! c this fxyyfxxy
      dd1=(beis+0.5/at)*(aejs+0.5/at)
      dd2=2.*be(i)*ae(j)*(ae(j)*f(i)+be(i)*f(j))+0.5*(beis &
     &+aejs)/at+(be(i)*f(i)+ae(j)*f(j))/at+0.5/ats
      dd3=fis*aejs+fjs*beis+4.*be(i)*ae(j)*f(i)* &
     &f(j)+(be(i)*f(i)+ae(j)*f(j))/at+0.5*(fis+fjs)/at+0.25 &
     &/ats
      dd4=2.*f(i)*f(j)*(ae(j)*f(i)+be(i)*f(j))+0.5*(fis+fjs &
     &)/at
      d41=be(j)*dd1+ae(j)*(beis+0.5/at)/at
      d42=f(j)*dd1+be(j)*dd2+be(i)*(ae(j)*f(i)+be(i)*f(j))/at+be(i) &
     &*ae(j)*f(i)/at+beis*ae(j)/at+ae(j)/ats+0.5*f(j)/ats
      d43=be(j)*dd3+f(j)*dd2+ae(j)*fis/at+beis*f(j)/at &
     &+2.*be(i)*f(i)*f(j)/at+2.*be(i)*ae(j)*f(i)/at+0.5*ae(j)/ &
     &ats+1.*f(j)/ats
      d44=be(j)*dd4+f(j)*dd3+f(i)*(ae(j)*f(i)+be(i)*f(j))/at+f(i) &
     &**2*f(j)/at+be(i)*f(i)*f(j)/at+0.5*f(j)/ats
      d45=be(j)*fis*fjs+f(j)*dd4+f(j)*fis/at
      d2=ae(j)*(ae(j)*f(i)+be(i)*f(j))/at+be(i)/ats+be(i)*ae(j) &
     &**2/at+be(i)*ae(j)*f(j)/at+0.5*f(i)/ats
      d3=aejs*f(i)/at+be(i)*fjs/at+ae(j)*f(i)*f(j)*2./at &
     &+2.*be(i)*ae(j)*f(j)/at+f(i)/ats+0.5*be(i)/ats
      d4=f(j)*(ae(j)*f(i)+be(i)*f(j))/at+ae(j)*f(i)*f(j)/at+f(i)* &
     &fjs/at+0.5*f(i)/ats
      dd1=ae(i)*d41+be(i)*be(j)*(aejs+0.5/at)/at+be(i)*ae(j)/ats
      dd2=f(i)*d41+ae(i)*d42+be(i)*f(j)*(aejs+0.5/at)/at+be(j) &
     &*d2+(ae(j)*f(i)+be(i)*f(j))/ats+2.*be(i)*ae(j)/ats
      dd3=f(i)*d42+ae(i)*d43+f(j)*d2+be(j)*d3+2.*ae(j)*f(i)/ats &
     &+2.*be(i)*f(j)/ats+f(i)*f(j)/ats+be(i)*ae(j)/ats
      dd4=f(i)*d43+ae(i)*d44+f(j)*d3+be(j)*d4+(ae(j)*f(i)+be(i)* &
     &f(j))/ats+2.*f(i)*f(j)/ats
      dd5=f(i)*d44+ae(i)*d45+f(j)*d4+be(j)*fjs*f(i)/at+f(i) &
     &*f(j)/ats
      dd6=f(i)*d45+ae(i)*fis*fjc+f(i)*fjc/at
      dd7=fic*fjc
      dr(mb,ma)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10* &
     &dd6+s12*dd7)
! c this fxxyfxyy
      dd1=(aeis+0.5/at)*(bejs+0.5/at)
      dd2=2.*ae(i)*be(j)*(be(j)*f(i)+ae(i)*f(j))+0.5*(aeis &
     &+bejs)/at+(ae(i)*f(i)+be(j)*f(j))/at+0.5/ats
      dd3=fis*bejs+fjs*aeis+4.*ae(i)*be(j)*f(i)* &
     &f(j)+(ae(i)*f(i)+be(j)*f(j))/at+0.5*(fis+fjs)/at+0.25 &
     &/ats
      dd4=2.*f(i)*f(j)*(be(j)*f(i)+ae(i)*f(j))+0.5*(fis+fjs &
     &)/at
      d41=ae(j)*dd1+be(j)*(aeis+0.5/at)/at
      d42=f(j)*dd1+ae(j)*dd2+ae(i)*(be(j)*f(i)+ae(i)*f(j))/at+ae(i) &
     &*be(j)*f(i)/at+aeis*be(j)/at+be(j)/ats+0.5*f(j)/ats
      d43=ae(j)*dd3+f(j)*dd2+be(j)*fis/at+aeis*f(j)/at &
     &+2.*ae(i)*f(i)*f(j)/at+2.*ae(i)*be(j)*f(i)/at+0.5*be(j)/ &
     &ats+1.*f(j)/ats
      d44=ae(j)*dd4+f(j)*dd3+f(i)*(be(j)*f(i)+ae(i)*f(j))/at+f(i) &
     &**2*f(j)/at+ae(i)*f(i)*f(j)/at+0.5*f(j)/ats
      d45=ae(j)*fis*fjs+f(j)*dd4+f(j)*fis/at
      d2=be(j)*(be(j)*f(i)+ae(i)*f(j))/at+ae(i)/ats+ae(i)*be(j) &
     &**2/at+ae(i)*be(j)*f(j)/at+0.5*f(i)/ats
      d3=bejs*f(i)/at+ae(i)*fjs/at+be(j)*f(i)*f(j)*2./at &
     &+2.*ae(i)*be(j)*f(j)/at+f(i)/ats+0.5*ae(i)/ats
      d4=f(j)*(be(j)*f(i)+ae(i)*f(j))/at+be(j)*f(i)*f(j)/at+f(i)* &
     &fjs/at+0.5*f(i)/ats
      dd1=be(i)*d41+ae(i)*ae(j)*(bejs+0.5/at)/at+ae(i)*be(j)/ats
      dd2=f(i)*d41+be(i)*d42+ae(i)*f(j)*(bejs+0.5/at)/at+ae(j) &
     &*d2+(be(j)*f(i)+ae(i)*f(j))/ats+2.*ae(i)*be(j)/ats
      dd3=f(i)*d42+be(i)*d43+f(j)*d2+ae(j)*d3+2.*be(j)*f(i)/ats &
     &+2.*ae(i)*f(j)/ats+f(i)*f(j)/ats+ae(i)*be(j)/ats
      dd4=f(i)*d43+be(i)*d44+f(j)*d3+ae(j)*d4+(be(j)*f(i)+ae(i)* &
     &f(j))/ats+2.*f(i)*f(j)/ats
      dd5=f(i)*d44+be(i)*d45+f(j)*d4+ae(j)*fjs*f(i)/at+f(i) &
     &*f(j)/ats
      dd6=f(i)*d45+be(i)*fis*fjc+f(i)*fjc/at
      dd7=fic*fjc
      dr(ma,mb)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10* &
     &dd6+s12*dd7)
! c this fxxxfyyy
      d51=aeic*bejs+0.5*(ae(i)*bejs+aeic+1.5*ae(i) &
     &/at+2.*ae(i)*bejs)/at
      d52=3.*aeis*bejs*f(i)+2.*aeic*be(j)*f(j)+0.5* &
     &aeic/at+1.5*bejs*ae(i)/at+1.5*aeis*f(i)/at+3.* &
     &ae(i)*be(j)*f(j)/at+1.5*ae(i)/ats+1.5*bejs*f(i)/at+0.75* &
     &f(i)/ats
      d53=3.*ae(i)*bejs*fis+6.*aeis*be(j)*f(i)*f(j) &
     &+1.5*aeis*f(i)/at+1.5*bejs*f(i)/at+1.5*ae(i)*fis &
     &/at+3.*be(j)*f(i)*f(j)/at+1.5*f(i)/ats+aeic*fjs &
     &+3.*ae(i)*be(j)*f(j)/at+1.5*ae(i)*fjs/at+0.75*ae(i)/ats
      d54=6.*fis*f(j)*ae(i)*be(j)+3.*f(i)*fjs*aeis+1.5 &
     &*ae(i)*fis/at+1.5*ae(i)*fjs/at+bejs*fic+3.* &
     &be(j)*f(i)*f(j)/at+0.5*fic/at+1.5*f(i)*fjs/at &
     &+0.75*f(i)/ats
      d55=3.*ae(i)*fis*fjs+2.*be(j)*fic*f(j)+1.5*f(i) &
     &*fjs/at+0.5*fic/at
      dd1=be(j)*d51+aeic*be(j)/at+1.5*ae(i)*be(j)/ats
      dd2=f(j)*d51+be(j)*d52+3.*aeis*be(j)*f(i)/at+aeic* &
     &f(j)/at+aeic*be(j)/at+3.*ae(i)*be(j)/ats+1.5*ae(i) &
     &*f(j)/ats+1.5*be(j)*f(i)/ats
      dd3=f(j)*d52+be(j)*d53+3.*ae(i)*be(j)*fis/at+3.*ae(i) &
     &**2*f(i)*f(j)/at+3.*aeis*be(j)*f(i)/at+1.5*be(j)*f(i) &
     &/ats+1.5*f(i)*f(j)/ats+be(j)*f(i)*1.5/ats+aeic*f(j) &
     &/at+3.*ae(i)*f(j)/ats+1.5*ae(i)*be(j)/ats
      dd4=f(j)*d53+be(j)*d54+3.*ae(i)*f(i)*(be(j)*f(i)+f(i)*f(j) &
     &+ae(i)*f(j))/at+1.5*ae(i)*f(j)/ats+be(j)*fic/at &
     &+3.*f(i)*(f(j)+0.5*be(j))/ats
      dd5=f(j)*d54+be(j)*d55+fis*(3.*ae(i)*f(j)+f(i)*f(j) &
     &+be(j)*f(i))/at+1.5*f(i)*f(j)/ats
      dd6=f(j)*d55+be(j)*fic*fjs+fic*f(j)/at
      dd7=fjc*fic
      dr(na,nb)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5 &
     &-s10*dd6+s12*dd7)
! c this fyyyfxxx
      d51=beic*aejs+0.5*(be(i)*aejs+beic+1.5*be(i) &
     &/at+2.*be(i)*aejs)/at
      d52=3.*beis*aejs*f(i)+2.*beic*ae(j)*f(j)+0.5* &
     &beic/at+1.5*aejs*be(i)/at+1.5*beis*f(i)/at+3.* &
     &be(i)*ae(j)*f(j)/at+1.5*be(i)/ats+1.5*aejs*f(i)/at+0.75* &
     &f(i)/ats
      d53=3.*be(i)*aejs*fis+6.*beis*ae(j)*f(i)*f(j)  &
     &+1.5*beis*f(i)/at+1.5*aejs*f(i)/at+1.5*be(i)*fis &
     &/at+3.*ae(j)*f(i)*f(j)/at+1.5*f(i)/ats+beic*fjs &
     &+3.*be(i)*ae(j)*f(j)/at+1.5*be(i)*fjs/at+0.75*be(i)/ats
      d54=6.*fis*f(j)*be(i)*ae(j)+3.*f(i)*fjs*beis+1.5 &
     &*be(i)*fis/at+1.5*be(i)*fjs/at+aejs*fic+3.* &
     &ae(j)*f(i)*f(j)/at+0.5*fic/at+1.5*f(i)*fjs/at &
     &+0.75*f(i)/ats
      d55=3.*be(i)*fis*fjs+2.*ae(j)*fic*f(j)+1.5*f(i) &
     &*fjs/at+0.5*fic/at
      dd1=ae(j)*d51+beic*ae(j)/at+1.5*be(i)*ae(j)/ats
      dd2=f(j)*d51+ae(j)*d52+3.*beis*ae(j)*f(i)/at+beic* &
     &f(j)/at+beic*ae(j)/at+3.*be(i)*ae(j)/ats+1.5*be(i) &
     &*f(j)/ats+1.5*ae(j)*f(i)/ats
      dd3=f(j)*d52+ae(j)*d53+3.*be(i)*ae(j)*fis/at+3.*be(i) &
     &**2*f(i)*f(j)/at+3.*beis*ae(j)*f(i)/at+1.5*ae(j)*f(i) &
     &/ats+1.5*f(i)*f(j)/ats+ae(j)*f(i)*1.5/ats+beic*f(j) &
     &/at+3.*be(i)*f(j)/ats+1.5*be(i)*ae(j)/ats
      dd4=f(j)*d53+ae(j)*d54+3.*be(i)*f(i)*(ae(j)*f(i)+f(i)*f(j) &
     &+be(i)*f(j))/at+1.5*be(i)*f(j)/ats+ae(j)*fic/at &
     &+3.*f(i)*(f(j)+0.5*ae(j))/ats
      dd5=f(j)*d54+ae(j)*d55+fis*(3.*be(i)*f(j)+f(i)*f(j) &
     &+ae(j)*f(i))/at+1.5*f(i)*f(j)/ats
      dd6=f(j)*d55+ae(j)*fic*fjs+fic*f(j)/at
      dd7=fjc*fic
      dr(nb,na)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5 &
     &-s10*dd6+s12*dd7)
  215 continue
      continue
      return
      end subroutine pot1f1
! c
!*****************************************************
! c
      subroutine pot1f2(dr,f,ae,be,at,ats,coef,xx)

      use O_Kinds
      use O_Constants
!      use genmatrixdatamodule

      implicit none
      real (kind=double), dimension (3,3) :: dr(20,20)

      integer :: i,j,k,l,ko,kp,lo,lp,mq,mr
      real (kind=double) :: s0, s2, s4, s6, s8, s10, s12
      real (kind=double) :: aeis, beis, aeic, beic, fis, fic
      real (kind=double) :: aejs, bejs, aejc, bejc, fjs, fjc
      real (kind=double) :: aeks, beks, aekc, bekc, fks, fkc
      real (kind=double) :: fifo, fkfo
      real (kind=double) :: d1, d2, d3
      real (kind=double) :: dd1, dd2, dd3, dd4, dd5, dd6, dd7
      real (kind=double) :: ddd1, ddd2, ddd3, ddd4, ddd5, ddd6
      real (kind=double) :: d01, d02, d03, d04, d001, d002, d003, d004
      real (kind=double) :: d11, d12, d13, d14
      real (kind=double) :: d21, d22, d23, d24
      real (kind=double) :: d41, d42, d43, d44, d45
      real (kind=double) :: d51, d52, d53, d54, d55
      real (kind=double) :: d61, d62, d63, d64
      real (kind=double) :: d71, d72, d73, d74
      real (kind=double) :: d81, d82, d83, d84

      real (kind=double), dimension(3) :: f,ae,be
      real (kind=double) :: at,ats,coef,xx
! c
      call fmtcf(xx,s0,s2,s4,s6,s8,s10,s12)
! c
! c this loop for x-y-z, no symmetry elements.
      i=1
      j=2
      k=3
      fis=f(i)*f(i)
      fjs=f(j)*f(j)
      fks=f(k)*f(k)
! c this fxyzs
      dd1=ae(k)*ae(i)*ae(j)
      dd2=(ae(j)*f(i)+ae(i)*f(j))*ae(k)+ae(i)*ae(j)*f(k)
      dd3=f(k)*(ae(j)*f(i)+ae(i)*f(j))+f(i)*f(j)*ae(k)
      dd4=f(i)*f(j)*f(k)
      dr(11,1)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4)
! c this sfxyz
      dd1=be(k)*be(i)*be(j)
      dd2=(be(j)*f(i)+be(i)*f(j))*be(k)+be(i)*be(j)*f(k)
      dd3=f(k)*(be(j)*f(i)+be(i)*f(j))+f(i)*f(j)*be(k)
      dd4=f(i)*f(j)*f(k)
      dr(1,11)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4)
! c this fxyzfxyz
      d71=(ae(i)*be(i)+0.5/at)*(ae(j)*be(j)+0.5/at)
      d72=(ae(i)*f(i)+be(i)*f(i)+0.5/at)*(ae(j)*be(j)+0.5/at) &
     &+(ae(i)*be(i)+0.5/at)*(ae(j)*f(j)+be(j)*f(j)+0.5/at)
      d73=(ae(i)*f(i)+be(i)*f(i)+0.5/at)*(ae(j)*f(j)+be(j)*f(j) &
     &+0.5/at)+(ae(i)*be(i)+0.5/at)*fjs+(ae(j)*be(j)+0.5/at) &
     &*fis
      d74=(ae(i)*f(i)+be(i)*f(i)+0.5/at)*fjs+(ae(j)*f(j)+ &
     &be(j)*f(j)+0.5/at)*fis
      dd1=be(k)*ae(k)*d71+0.5*d71/at
      dd2=f(k)*ae(k)*d71+be(k)*(f(k)*d71+ae(k)*d72)+(0.5*d71/at &
     &+0.5*d72/at)
      dd3=f(k)*(f(k)*d71+ae(k)*d72)+be(k)*(f(k)*d72+ae(k)*d73) &
     &+(0.5*d72/at+0.5*d73/at)
      dd4=f(k)*(f(k)*d72+ae(k)*d73)+be(k)*(f(k)*d73+ae(k)*d74) &
     &+0.5*(d73+d74)/at
      dd5=f(k)*(f(k)*d73+ae(k)*d74)+be(k)*(f(k)*d74+ae(k)*fis*fjs &
     &)+0.5*(d74+fis*fjs)/at
      dd6=f(k)*(f(k)*d74+ae(k)*fis*fjs)+be(k)*f(k)*fis &
     &*fjs+0.5*fis*fjs/at
      dd7=fks*fjs*fis
      dr(11,11)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5 &
     &-s10*dd6+s12*dd7)
! c continue
! c this ioop has exchange in x-y-z,y-x,z-y,has three symetry elements
      do 220 i=1,3
      aeis=ae(i)*ae(i)
      beis=be(i)*be(i)
      aeic=aeis*ae(i)
      beic=beis*be(i)
      fis=f(i)*f(i)
      fic=fis*f(i)
      do 220 j=1,2
      aejs=ae(j)*ae(j)
      bejs=be(j)*be(j)
      aejc=aejs*ae(j)
      bejc=bejs*be(j)
      fjs=f(j)*f(j)
      fjc=fjs*f(j)
      do 220 k=2,3
      aeks=ae(k)*ae(k)
      beks=be(k)*be(k)
      aekc=aeks*ae(k)
      bekc=beks*be(k)
      fks=f(k)*f(k)
      fkc=fks*f(k)
      if(i.ne.j.and.j.ne.k.and.k.ne.i) go to 221
      go to 220
  221 continue
! c this fxyzpx
      dd1=ae(k)*ae(i)*ae(j)*be(i)+0.5*ae(k)*ae(j)/at
      dd2=ae(k)*ae(i)*ae(j)*f(i)+ae(j)*be(i)*f(i)*ae(k)+ae(i)*be(i) &
     &*ae(k)*f(j)+ae(i)*ae(j)*be(i)*f(k)+0.5*(ae(j)*ae(k)+f(j) &
     &*ae(k)+ae(j)*f(k))/at
      dd3=ae(j)*fis*ae(k)+ae(i)*f(i)*f(j)*ae(k)+ae(i)*ae(j) &
     &*f(i)*f(k)+ae(j)*be(i)*f(i)*f(k)+ae(i)*be(i)*f(j)*f(k)+be(i) &
     &*f(i)*f(j)*ae(k)+0.5*(ae(j)*f(k)+f(j)*f(k)+f(j)*ae(k))/at
      dd4=ae(j)*fis*f(k)+ae(i)*f(i)*f(j)*f(k)+fis*f(j)*ae(k) &
     &+be(i)*f(i)*f(j)*f(k)+0.5*f(j)*f(k)/at
      dd5=fis*f(j)*f(k)
      dr(11,i+1)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5)
! c this pxfxyz
      dd1=be(k)*be(i)*be(j)*ae(i)+0.5*be(k)*be(j)/at
      dd2=be(k)*be(i)*be(j)*f(i)+be(j)*ae(i)*f(i)*be(k)+be(i)*ae(i) &
     &*be(k)*f(j)+be(i)*be(j)*ae(i)*f(k)+0.5*(be(j)*be(k)+f(j) &
     &*be(k)+be(j)*f(k))/at
      dd3=be(j)*fis*be(k)+be(i)*f(i)*f(j)*be(k)+be(i)*be(j) &
     &*f(i)*f(k)+be(j)*ae(i)*f(i)*f(k)+be(i)*ae(i)*f(j)*f(k)+ae(i) &
     &*f(i)*f(j)*be(k)+0.5*(be(j)*f(k)+f(j)*f(k)+f(j)*be(k))/at
      dd4=be(j)*fis*f(k)+be(i)*f(i)*f(j)*f(k)+fis*f(j)*be(k) &
     &+ae(i)*f(i)*f(j)*f(k)+0.5*f(j)*f(k)/at
      dd5=fis*f(j)*f(k)
      dr(i+1,11)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5)
! c this fxyzdxx
      d1=beis*ae(i)+0.5*(2.*be(i)+ae(i))/at
      d2=beis*f(i)+0.5*(ae(i)+2.*be(i)+3.*f(i))/at+2.*be(i) &
     &*ae(i)*f(i)
      d3=fis*(2.*be(i)+ae(i))+1.5*f(i)/at
      dd1=ae(j)*d1
      dd2=f(j)*d1+ae(j)*d2
      dd3=f(j)*d2+ae(j)*d3
      dd4=f(j)*d3+ae(j)*fic
      ddd1=ae(k)*dd1
      ddd2=f(k)*dd1+ae(k)*dd2
      ddd3=f(k)*dd2+ae(k)*dd3
      ddd4=f(k)*dd3+ae(k)*dd4
      ddd5=f(k)*dd4+ae(k)*f(j)*fic
      ddd6=f(k)*f(j)*fic
      dr(11,i+4)=coef*(s0*ddd1-s2*ddd2+s4*ddd3-s6*ddd4+s8*ddd5 &
     &-s10*ddd6)
! c this dxxfxyz
      d1=aeis*be(i)+0.5*(2.*ae(i)+be(i))/at
      d2=aeis*f(i)+0.5*(be(i)+2.*ae(i)+3.*f(i))/at+2.*ae(i) &
     &*be(i)*f(i)
      d3=fis*(2.*ae(i)+be(i))+1.5*f(i)/at
      dd1=be(j)*d1
      dd2=f(j)*d1+be(j)*d2
      dd3=f(j)*d2+be(j)*d3
      dd4=f(j)*d3+be(j)*fic
      ddd1=be(k)*dd1
      ddd2=f(k)*dd1+be(k)*dd2
      ddd3=f(k)*dd2+be(k)*dd3
      ddd4=f(k)*dd3+be(k)*dd4
      ddd5=f(k)*dd4+be(k)*f(j)*fic
      ddd6=f(k)*f(j)*fic
      dr(i+4,11)=coef*(s0*ddd1-s2*ddd2+s4*ddd3-s6*ddd4+s8*ddd5 &
     &-s10*ddd6)
! c this fxyzdyz
      d71=(ae(j)*be(j)+0.5/at)*(ae(k)*be(k)+0.5/at)
      d72=(ae(j)*f(j)+be(j)*f(j)+0.5/at)*(ae(k)*be(k)+0.5/at) &
     &+(ae(j)*be(j)+0.5/at)*(ae(k)*f(k)+be(k)*f(k)+0.5/at)
      d73=(ae(j)*f(j)+be(j)*f(j)+0.5/at)*(ae(k)*f(k)+be(k)*f(k) &
     &+0.5/at)+(ae(j)*be(j)+0.5/at)*fks+(ae(k)*be(k)+0.5/at) &
     &*fjs
      d74=(ae(j)*f(j)+be(j)*f(j)+0.5/at)*fks+(ae(k)*f(k)+ &
     &be(k)*f(k)+0.5/at)*fjs
      dd1=ae(i)*d71
      dd2=f(i)*d71+ae(i)*d72
      dd3=f(i)*d72+ae(i)*d73
      dd4=f(i)*d73+ae(i)*d74
      dd5=f(i)*d74+ae(i)*fjs*fks
      dd6=f(i)*fjs*fks
      dr(11,11-i)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6)
! c thjs dyzfxyz
      d71=(be(j)*ae(j)+0.5/at)*(be(k)*ae(k)+0.5/at)
      d72=(be(j)*f(j)+ae(j)*f(j)+0.5/at)*(be(k)*ae(k)+0.5/at) &
     &+(be(j)*ae(j)+0.5/at)*(be(k)*f(k)+ae(k)*f(k)+0.5/at)
      d73=(be(j)*f(j)+ae(j)*f(j)+0.5/at)*(be(k)*f(k)+ae(k)*f(k) &
     &+0.5/at)+(be(j)*ae(j)+0.5/at)*fks+(be(k)*ae(k)+0.5/at) &
     &*fjs
      d74=(be(j)*f(j)+ae(j)*f(j)+0.5/at)*fks+(be(k)*f(k)+ &
     &ae(k)*f(k)+0.5/at)*fjs
      dd1=be(i)*d71
      dd2=f(i)*d71+be(i)*d72
      dd3=f(i)*d72+be(i)*d73
      dd4=f(i)*d73+be(i)*d74
      dd5=f(i)*d74+be(i)*fjs*fks
      dd6=f(i)*fjs*fks
      dr(11-i,11)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6)
! c this fxyzfxxx
      d61=beic*ae(i)+1.5*be(i)*ae(i)/at+1.5*beis/at+0.75/ats
      d62=beic*f(i)+4.5*be(i)*f(i)/at+3.*beis*ae(i)*f(i)+1.5*be(i) &
     &*ae(i)/at+1.5*ae(i)*f(i)/at+1.5*beis/at+1.5/ats
      d63=3.*beis*fis+4.5*be(i)*f(i)/at+3.*be(i)*ae(i)* &
     &fis+1.5*f(i)*ae(i)/at+3.*fis/at+0.75/ats
      d64=3.*be(i)*fic+3.*fis/at+fic*ae(i)
      dd1=ae(k)*ae(j)*d61
      dd2=f(k)*ae(j)*d61+ae(k)*(f(j)*d61+ae(j)*d62)
      dd3=f(k)*(f(j)*d61+ae(j)*d62)+ae(k)*(f(j)*d62+ae(j)*d63)
      dd4=f(k)*(f(j)*d62+ae(j)*d63)+ae(k)*(f(j)*d63+ae(j)*d64)
      dd5=f(k)*(f(j)*d63+ae(j)*d64)+ae(k)*(f(j)*d64+ae(j)*f(i)**4)
      dd6=f(k)*(f(j)*d64+ae(j)*f(i)**4)+ae(k)*f(j)*f(i)**4
      dd7=f(k)*f(j)*f(i)**4
      dr(11,i+17)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5 &
     &-s10*dd6+s12*dd7)
! c this fxxxfxyz
      d61=aeic*be(i)+1.5*ae(i)*be(i)/at+1.5*aeis/at+0.75/ats
      d62=aeic*f(i)+4.5*ae(i)*f(i)/at+3.*aeis*be(i)*f(i)+1.5*ae(i) &
     &*be(i)/at+1.5*be(i)*f(i)/at+1.5*aeis/at+1.5/ats
      d63=3.*aeis*fis+4.5*ae(i)*f(i)/at+3.*ae(i)*be(i)* &
     &fis+1.5*f(i)*be(i)/at+3.*fis/at+0.75/ats
      d64=3.*ae(i)*fic+3.*fis/at+fic*be(i)
      dd1=be(k)*be(j)*d61
      dd2=f(k)*be(j)*d61+be(k)*(f(j)*d61+be(j)*d62)
      dd3=f(k)*(f(j)*d61+be(j)*d62)+be(k)*(f(j)*d62+be(j)*d63)
      dd4=f(k)*(f(j)*d62+be(j)*d63)+be(k)*(f(j)*d63+be(j)*d64)
      dd5=f(k)*(f(j)*d63+be(j)*d64)+be(k)*(f(j)*d64+be(j)*f(i)**4)
      dd6=f(k)*(f(j)*d64+be(j)*f(i)**4)+be(k)*f(j)*f(i)**4
      dd7=f(k)*f(j)*f(i)**4
      dr(i+17,11)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5 &
     &-s10*dd6+s12*dd7)
! c this fxxxdyz
      d21=be(j)*(aeic+1.5*ae(i)/at)
      d22=aeic*f(j)+1.5*ae(i)*f(j)/at+3.*aeis*be(j)*f(i) &
     &+1.5*ae(i)*be(j)/at+1.5*be(j)*f(i)/at
      d23=3.*aeis*f(i)*f(j)+1.5*ae(i)*f(j)/at+1.5*f(j)*f(i)/at &
     &+3.*ae(i)*be(j)*fis+1.5*f(i)*be(j)/at
      d24=3.*ae(i)*fis*f(j)+1.5*f(i)*f(j)/at+fic*be(j)
      dd1=be(k)*d21
      dd2=f(k)*d21+be(k)*d22
      dd3=f(k)*d22+be(k)*d23
      dd4=f(k)*d23+be(k)*d24
      dd5=f(k)*d24+be(k)*f(j)*fic
      dd6=f(k)*f(j)*fic
      dr(i+17,11-i)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5 &
     &-s10*dd6)
! c this dyzfxxx
      d21=ae(j)*(beic+1.5*be(i)/at)
      d22=beic*f(j)+1.5*be(i)*f(j)/at+3.*beis*ae(j)*f(i) &
     &+1.5*be(i)*ae(j)/at+1.5*ae(j)*f(i)/at
      d23=3.*beis*f(i)*f(j)+1.5*be(i)*f(j)/at+1.5*f(j)*f(i)/at &
     &+3.*be(i)*ae(j)*fis+1.5*f(i)*ae(j)/at
      d24=3.*be(i)*fis*f(j)+1.5*f(i)*f(j)/at+fic*ae(j)
      dd1=ae(k)*d21
      dd2=f(k)*d21+ae(k)*d22
      dd3=f(k)*d22+ae(k)*d23
      dd4=f(k)*d23+ae(k)*d24
      dd5=f(k)*d24+ae(k)*f(j)*fic
      dr(11-i,i+17)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6)
  220 continue
! c this loop has exchange x-y-z and x-y, y-z, z-x, has six symmetry
! c elements
      l=11
      do 230 i=1,3
      aeis=ae(i)*ae(i)
      beis=be(i)*be(i)
      aeic=aeis*ae(i)
      beic=beis*be(i)
      fis=f(i)*f(i)
      fic=fis*f(i)
      fifo=fis*fis
      do 230 j=1,3
      aejs=ae(j)*ae(j)
      bejs=be(j)*be(j)
      aejc=aejs*ae(j)
      bejc=bejs*be(j)
      fjs=f(j)*f(j)
      fjc=fjs*f(j)
      do 230 k=1,3
      aeks=ae(k)*ae(k)
      beks=be(k)*be(k)
      aekc=aeks*ae(k)
      bekc=beks*be(k)
      fks=f(k)*f(k)
      fkc=fks*f(k)
      fkfo=fks*fks
      if(i.ne.j.and.j.ne.k.and.k.ne.i) go to 231
      go to 230
  231 continue
      l=l+1
! c this fxxypz
      dd1=be(k)*ae(j)*(aeis+0.5/at)
      dd2=aeis*ae(j)*f(k)+0.5*ae(j)*f(k)/at+2.*ae(i)*ae(j)*be(k) &
     &*f(i)+0.5*ae(j)*be(k)/at+aeis*be(k)*f(j)+0.5*f(j)*be(k)/at
      dd3=2.*ae(i)*ae(j)*f(i)*f(k)+0.5*ae(j)*f(k)/at+aeis*f(j) &
     &*f(k)+0.5*f(j)*f(k)/at+be(k)*(ae(j)*fis+2.*ae(i)*f(i)*f(j) &
     &+0.5*f(j)/at)
      dd4=f(k)*f(i)*(ae(j)*f(i)+2.*ae(i)*f(j))+0.5*f(j)*f(k)/at &
     &  +be(k)*fis*f(j)
      dd5=f(k)*fis*f(j)
      dr(l,k+1)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5)
! c this pzfxxy
      dd1=ae(k)*be(j)*(beis+0.5/at)
      dd2=beis*be(j)*f(k)+0.5*be(j)*f(k)/at+2.*be(i)*be(j)*ae(k) &
     &*f(i)+0.5*be(j)*ae(k)/at+beis*ae(k)*f(j)+0.5*f(j)*ae(k)/at
      dd3=2.*be(i)*be(j)*f(i)*f(k)+0.5*be(j)*f(k)/at+beis*f(j) &
     &*f(k)+0.5*f(j)*f(k)/at+ae(k)*(be(j)*fis+2.*be(i)*f(i)*f(j) &
     &+0.5*f(j)/at)
      dd4=f(k)*f(i)*(be(j)*f(i)+2.*be(i)*f(j))+0.5*f(j)*f(k)/at &
     &  +ae(k)*fis*f(j)
      dd5=f(k)*fis*f(j)
      dr(k+1,l)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5)
! c this fxxydyz
      d11=aeis*be(j)*ae(j)+0.5*ae(j)*be(j)/at+0.5*aeis/ &
     &at+0.25/ats
      d12=aeis*ae(j)*f(j)+0.5*ae(j)*f(j)/at+2.*ae(i)*ae(j)* &
     &be(j)*f(i)+0.5*ae(j)*be(j)/at+aeis*be(j)*f(j)+0.5*be(j) &
     &*f(j)/at+ae(i)*f(i)/at+0.5/ats+0.5*aeis/at
      d13=2.*ae(i)*ae(j)*f(i)*f(j)+0.5*f(j)*ae(j)/at+aeis*fjs &
     &+0.5*fjs/at+ae(j)*be(j)*fis+2.*ae(i)*be(j)* &
     &f(i)*f(j)+0.5*be(j)*f(j)/at+0.5*fis/at+ae(i)*f(i)/at+0.25/ats
      d14=ae(j)*fis*f(j)+2.*ae(i)*fjs*f(i)+0.5*fjs &
     & /at+be(j)*fis*f(j)+0.5*fis/at
      dd1=be(k)*d11
      dd2=f(k)*d11+be(k)*d12
      dd3=f(k)*d12+be(k)*d13
      dd4=f(k)*d13+be(k)*d14
      dd5=f(k)*d14+be(k)*fjs*fis
      dd6=f(k)*fjs*fis
      dr(l,11-i)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6)
! c this dyzfxxy
      d11=beis*ae(j)*be(j)+0.5*be(j)*ae(j)/at+0.5*beis/ &
     &at+0.25/ats
      d12=beis*be(j)*f(j)+0.5*be(j)*f(j)/at+2.*be(i)*be(j)* &
     &ae(j)*f(i)+0.5*be(j)*ae(j)/at+beis*ae(j)*f(j)+0.5*ae(j) &
     & *f(j)/at+be(i)*f(i)/at+0.5/ats+0.5*beis/at
      d13=2.*be(i)*be(j)*f(i)*f(j)+0.5*f(j)*be(j)/at+beis*fjs &
     & +0.5*fjs/at+be(j)*ae(j)*fis+2.*be(i)*ae(j)* &
     & f(i)*f(j)+0.5*ae(j)*f(j)/at+0.5*fis/at+be(i)*f(i)/at+0.25/ats
      d14=be(j)*fis*f(j)+2.*be(i)*fjs*f(i)+0.5*fjs &
     & /at+ae(j)*fis*f(j)+0.5*fis/at
      dd1=ae(k)*d11
      dd2=f(k)*d11+ae(k)*d12
      dd3=f(k)*d12+ae(k)*d13
      dd4=f(k)*d13+ae(k)*d14
      dd5=f(k)*d14+ae(k)*fjs*fis
      dd6=f(k)*fjs*fis
      dr(11-i,l)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6)
! c this fxxydzx
      if(j.eq.3) mq=8
      if(j.eq.1) mq=10
      if(j.eq.2) mq=9
      d01=aeis*be(i)*ae(j)+.5*ae(j)*be(i)/at+ae(i)*ae(j)/at
      d02=f(i)*(aeis*ae(j)+1.5*ae(j)/at+2.*be(i)*ae(i)*ae(j))+(aeis &
     &*be(i)+0.5*be(i)/at+ae(i)/at)*f(j)+(0.5*ae(j)*be(i)+ae(i)*ae(j)) &
     &/at
      d03=2.*ae(i)*ae(j)*fis+1.5*ae(j)*f(i)/at+f(i)*f(j)*(aeis+2.* &
     & ae(i)*be(i)+1.5/at)+ae(j)*be(i)*fis+0.5*be(i)*f(j)/at &
     & +ae(i)*f(j)/at
      d04=fis*(ae(j)*f(i)+2.*ae(i)*f(j)+be(i)*f(j))+1.5*f(i)*f(j)/at
      dd1=be(k)*d01
      dd2=f(k)*d01+be(k)*d02
      dd3=f(k)*d02+be(k)*d03
      dd4=f(k)*d03+be(k)*d04
      dd5=f(k)*d04+be(k)*fic*f(j)
      dd6=f(k)*fic*f(j)
      dr(l,mq)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6)
! c this dzxfxxy
      d001=beis*ae(i)*be(j)+.5*be(j)*ae(i)/at+be(i)*be(j)/at
      d002=f(i)*(beis*be(j)+1.5*be(j)/at+2.*ae(i)*be(i)*be(j))+(beis &
     &*ae(i)+0.5*ae(i)/at+be(i)/at)*f(j)+(0.5*be(j)*ae(i)+be(i)*be(j)) &
     &/at
      d003=2.*be(i)*be(j)*fis+1.5*be(j)*f(i)/at+f(i)*f(j)*(beis+2.* &
     & be(i)*ae(i)+1.5/at)+be(j)*ae(i)*fis+0.5*ae(i)*f(j)/at &
     & +be(i)*f(j)/at
      d004=fis*(be(j)*f(i)+2.*be(i)*f(j)+ae(i)*f(j))+1.5*f(i)*f(j)/at
      dd1=ae(k)*d001
      dd2=f(k)*d001+ae(k)*d002
      dd3=f(k)*d002+ae(k)*d003
      dd4=f(k)*d003+ae(k)*d004
      dd5=f(k)*d004+ae(k)*fic*f(j)
      dd6=f(k)*fic*f(j)
      dr(mq,l)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6)
! c this fxxydzz
      d81=(beks+0.5/at)*(aeis+0.5/at)
      d82=2.*be(k)*ae(i)*(ae(i)*f(k)+be(k)*f(i))+0.5*(beks+aeis &
     &)/at+(be(k)*f(k)+ae(i)*f(i))/at+0.5/ats
      d83=fks*aeis+fis*beks+4.*be(k)*ae(i)*f(k)*f(i) &
     &+(be(k)*f(k)+ae(i)*f(i))/at+0.5*(fks+fis)/at+0.25/ats
      d84=2.*f(k)*f(i)*(ae(i)*f(k)+be(k)*f(i))+0.5*(fks+fis &
     &)/at
      dd1=ae(j)*d81
      dd2=f(j)*d81+ae(j)*d82
      dd3=f(j)*d82+ae(j)*d83
      dd4=f(j)*d83+ae(j)*d84
      dd5=f(j)*d84+ae(j)*fks*fis
      dd6=fis*fks*f(j)
      dr(l,k+4)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6)
! c this dzzfxxy
      d81=(aeks+0.5/at)*(beis+0.5/at)
      d82=2.*ae(k)*be(i)*(be(i)*f(k)+ae(k)*f(i))+0.5*(aeks+beis &
     &)/at+(ae(k)*f(k)+be(i)*f(i))/at+0.5/ats
      d83=fks*beis+fis*aeks+4.*ae(k)*be(i)*f(k)*f(i) &
     &+(ae(k)*f(k)+be(i)*f(i))/at+0.5*(fks+fis)/at+0.25/ats
      d84=2.*f(k)*f(i)*(be(i)*f(k)+ae(k)*f(i))+0.5*(fks+fis &
     &)/at
      dd1=be(j)*d81
      dd2=f(j)*d81+be(j)*d82
      dd3=f(j)*d82+be(j)*d83
      dd4=f(j)*d83+be(j)*d84
      dd5=f(j)*d84+be(j)*fks*fis
      dd6=fis*fks*f(j)
      dr(k+4,l)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10*dd6)
! c this fxxyfxyz
      dd1=be(j)*d01+0.5*(aeis*be(i)+0.5*be(i)/at+ae(i)/at)/at
      dd2=f(j)*d01+be(j)*d02+0.5*(aeis*f(i)+1.5*f(i)/at+2.*ae(i)/at+ &
     &2.*be(i)*ae(i)*f(i)+1.*be(i)/at+aeis*be(i))/at
      dd3=f(j)*d02+be(j)*d03+0.5*(2.*ae(i)*fis+3.0*f(i)/at &
     & +be(i)*fis+aeis*f(i)+2.*ae(i)*be(i)*f(i)+0.5*be(i)/at+ae(i) &
     & /at)/at
      dd4=f(j)*d03+be(j)*d04+0.5*fic/at+ae(i)*fis/at+0.75*f(i)/ats &
     &  +0.5*be(i)*fis/at
      dd5=f(j)*d04+be(j)*fic*f(j)+0.5*fic/at
      ddd1=be(k)*dd1
      ddd2=f(k)*dd1+be(k)*dd2
      ddd3=f(k)*dd2+be(k)*dd3
      ddd4=f(k)*dd3+be(k)*dd4
      ddd5=f(k)*dd4+be(k)*dd5
      ddd6=f(k)*dd5+be(k)*fjs*fic
      dr(l,11)=coef*(s0*ddd1-s2*ddd2+s4*ddd3-s6*ddd4+s8*ddd5 &
     &-s10*ddd6+s12*f(k)*fjs*fic)
! c this fxyzfxxy
      dd1=ae(j)*d001+0.5*(beis*ae(i)+0.5*ae(i)/at+be(i)/at)/at
      dd2=f(j)*d001+ae(j)*d002+0.5*(beis*f(i)+1.5*f(i)/at+2.*be(i)/at+ &
     &2.*ae(i)*be(i)*f(i)+1.*ae(i)/at+beis*ae(i))/at
      dd3=f(j)*d002+ae(j)*d003+0.5*(2.*be(i)*fis+3.0*f(i)/at &
     & +ae(i)*fis+beis*f(i)+2.*be(i)*ae(i)*f(i)+0.5*ae(i)/at+be(i) &
     & /at)/at
      dd4=f(j)*d003+ae(j)*d004+0.5*fic/at+be(i)*fis/at+0.75*f(i)/ats &
     &  +0.5*ae(i)*fis/at
      dd5=f(j)*d004+ae(j)*fic*f(j)+0.5*fic/at
      ddd1=ae(k)*dd1
      ddd2=f(k)*dd1+ae(k)*dd2
      ddd3=f(k)*dd2+ae(k)*dd3
      ddd4=f(k)*dd3+ae(k)*dd4
      ddd5=f(k)*dd4+ae(k)*dd5
      ddd6=f(k)*dd5+ae(k)*fjs*fic
      dr(11,l)=coef*(s0*ddd1-s2*ddd2+s4*ddd3-s6*ddd4+s8*ddd5  &
     &-s10*ddd6+s12*f(k)*fjs*fic)
! c this fxxyfyyz
      if(i.eq.1.and.j.eq.2) mr=15
      if(i.eq.1.and.j.eq.3) mr=17
      if(i.eq.2.and.j.eq.1) mr=13
      if(i.eq.2.and.j.eq.3) mr=16
      if(i.eq.3.and.j.eq.1) mr=12
      if(i.eq.3.and.j.eq.2) mr=14
      dd1=(aeis+0.5/at)*(bejs+0.5/at)
      dd2=2.*ae(i)*be(j)*(be(j)*f(i)+ae(i)*f(j))+0.5*(aeis &
     &+bejs)/at+(ae(i)*f(i)+be(j)*f(j))/at+0.5/ats
      dd3=fis*bejs+fjs*aeis+4.*ae(i)*be(j)*f(i)* &
     &f(j)+(ae(i)*f(i)+be(j)*f(j))/at+0.5*(fis+fjs)/at+0.25 &
     &/ats
      dd4=2.*f(i)*f(j)*(be(j)*f(i)+ae(i)*f(j))+0.5*(fis+fjs &
     &)/at
      d41=ae(j)*dd1+be(j)*(aeis+0.5/at)/at
      d42=f(j)*dd1+ae(j)*dd2+ae(i)*(be(j)*f(i)+ae(i)*f(j))/at+ae(i) &
     &*be(j)*f(i)/at+aeis*be(j)/at+be(j)/ats+0.5*f(j)/ats
      d43=ae(j)*dd3+f(j)*dd2+be(j)*fis/at+aeis*f(j)/at &
     &+2.*ae(i)*f(i)*f(j)/at+2.*ae(i)*be(j)*f(i)/at+0.5*be(j)/ &
     &ats+1.*f(j)/ats
      d44=ae(j)*dd4+f(j)*dd3+f(i)*(be(j)*f(i)+ae(i)*f(j))/at+f(i) &
     &**2*f(j)/at+ae(i)*f(i)*f(j)/at+0.5*f(j)/ats
      d45=ae(j)*fis*fjs+f(j)*dd4+f(j)*fis/at
      dd1=be(k)*d41
      dd2=f(k)*d41+be(k)*d42
      dd3=f(k)*d42+be(k)*d43
      dd4=f(k)*d43+be(k)*d44
      dd5=f(k)*d44+be(k)*d45
      dd6=f(k)*d45+be(k)*fjc*fis
      dd7=f(k)*fjc*fis
      dr(l,mr)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10 &
     &*dd6+s12*dd7)
! c this fyyzfxxy
      dd1=(beis+0.5/at)*(aejs+0.5/at)
      dd2=2.*be(i)*ae(j)*(ae(j)*f(i)+be(i)*f(j))+0.5*(beis &
     &+aejs)/at+(be(i)*f(i)+ae(j)*f(j))/at+0.5/ats
      dd3=fis*aejs+fjs*beis+4.*be(i)*ae(j)*f(i)* &
     &f(j)+(be(i)*f(i)+ae(j)*f(j))/at+0.5*(fis+fjs)/at+0.25 &
     &/ats
      dd4=2.*f(i)*f(j)*(ae(j)*f(i)+be(i)*f(j))+0.5*(fis+fjs &
     &)/at
      d41=be(j)*dd1+ae(j)*(beis+0.5/at)/at
      d42=f(j)*dd1+be(j)*dd2+be(i)*(ae(j)*f(i)+be(i)*f(j))/at+be(i) &
     &*ae(j)*f(i)/at+beis*ae(j)/at+ae(j)/ats+0.5*f(j)/ats
      d43=be(j)*dd3+f(j)*dd2+ae(j)*fis/at+beis*f(j)/at &
     &+2.*be(i)*f(i)*f(j)/at+2.*be(i)*ae(j)*f(i)/at+0.5*ae(j)/ &
     &ats+1.*f(j)/ats
      d44=be(j)*dd4+f(j)*dd3+f(i)*(ae(j)*f(i)+be(i)*f(j))/at+f(i) &
     &**2*f(j)/at+be(i)*f(i)*f(j)/at+0.5*f(j)/ats
      d45=be(j)*fis*fjs+f(j)*dd4+f(j)*fis/at
      dd1=ae(k)*d41
      dd2=f(k)*d41+ae(k)*d42
      dd3=f(k)*d42+ae(k)*d43
      dd4=f(k)*d43+ae(k)*d44
      dd5=f(k)*d44+ae(k)*d45
      dd6=f(k)*d45+ae(k)*fjc*fis
      dd7=f(k)*fjc*fis
      dr(mr,l)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5-s10 &
     &*dd6+s12*dd7)
! c this fxxyfzzz
      d51=bekc*aeis+0.5*(be(k)*aeis+bekc+1.5*be(k) &
     &/at+2.*be(k)*aeis)/at
      d52=3.*beks*aeis*f(k)+2.*bekc*ae(i)*f(i)+0.5* &
     &bekc/at+1.5*aeis*be(k)/at+1.5*beks*f(k)/at+3.* &
     &be(k)*ae(i)*f(i)/at+1.5*be(k)/ats+1.5*aeis*f(k)/at+0.75* &
     &f(k)/ats
      d53=3.*be(k)*aeis*fks+6.*beks*ae(i)*f(k)*f(i) &
     &+1.5*beks*f(k)/at+1.5*aeis*f(k)/at+1.5*be(k)*fks &
     &/at+3.*ae(i)*f(k)*f(i)/at+1.5*f(k)/ats+bekc*fis &
     &+3.*be(k)*ae(i)*f(i)/at+1.5*be(k)*fis/at+0.75*be(k)/ats
      d54=6.*fks*f(i)*be(k)*ae(i)+3.*f(k)*fis*beks+1.5 &
     &*be(k)*fks/at+1.5*be(k)*fis/at+aeis*fkc+3.* &
     &ae(i)*f(k)*f(i)/at+0.5*fkc/at+1.5*f(k)*fis/at &
     &+0.75*f(k)/ats
      d55=3.*be(k)*fks*fis+2.*ae(i)*fkc*f(i)+1.5*f(k) &
     &*fis/at+0.5*fkc/at
      dd1=ae(j)*d51
      dd2=f(j)*d51+ae(j)*d52
      dd3=f(j)*d52+ae(j)*d53
      dd4=f(j)*d53+ae(j)*d54
      dd5=f(j)*d54+ae(j)*d55
      dd6=f(j)*d55+ae(j)*fkc*fis
      dd7=f(j)*fkc*fis
      dr(l,k+17)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5 &
     &-s10*dd6+s12*dd7)
! c this fzzzfxxy
      d51=aekc*beis+0.5*(ae(k)*beis+aekc+1.5*ae(k) &
     &/at+2.*ae(k)*beis)/at
      d52=3.*aeks*beis*f(k)+2.*aekc*be(i)*f(i)+0.5* &
     &aekc/at+1.5*beis*ae(k)/at+1.5*aeks*f(k)/at+3.* &
     &ae(k)*be(i)*f(i)/at+1.5*ae(k)/ats+1.5*beis*f(k)/at+0.75* &
     &f(k)/ats
      d53=3.*ae(k)*beis*fks+6.*aeks*be(i)*f(k)*f(i) &
     &+1.5*aeks*f(k)/at+1.5*beis*f(k)/at+1.5*ae(k)*fks &
     &/at+3.*be(i)*f(k)*f(i)/at+1.5*f(k)/ats+aekc*fis &
     &+3.*ae(k)*be(i)*f(i)/at+1.5*ae(k)*fis/at+0.75*ae(k)/ats
      d54=6.*fks*f(i)*ae(k)*be(i)+3.*f(k)*fis*aeks+1.5 &
     &*ae(k)*fks/at+1.5*ae(k)*fis/at+beis*fkc+3.* &
     &be(i)*f(k)*f(i)/at+0.5*fkc/at+1.5*f(k)*fis/at &
     &+0.75*f(k)/ats
      d55=3.*ae(k)*fks*fis+2.*be(i)*fkc*f(i)+1.5*f(k) &
     &*fis/at+0.5*fkc/at
      dd1=be(j)*d51
      dd2=f(j)*d51+be(j)*d52
      dd3=f(j)*d52+be(j)*d53
      dd4=f(j)*d53+be(j)*d54
      dd5=f(j)*d54+be(j)*d55
      dd6=f(j)*d55+be(j)*fkc*fis
      dd7=f(j)*fkc*fis
      dr(k+17,l)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5 &
     &-s10*dd6+s12*dd7)
  230 continue
! c this loop
      do 250 i=1,3
      aeis=ae(i)*ae(i)
      beis=be(i)*be(i)
      aeic=aeis*ae(i)
      beic=beis*be(i)
      fis=f(i)*f(i)
      fic=fis*f(i)
      go to (1,2,3) ,i
    1 j=2
      aejs=ae(j)*ae(j)
      bejs=be(j)*be(j)
      aejc=aejs*ae(j)
      bejc=bejs*be(j)
      fjs=f(j)*f(j)
      fjc=fjs*f(j)
      k=3
      aeks=ae(k)*ae(k)
      beks=be(k)*be(k)
      aekc=aeks*ae(k)
      bekc=beks*be(k)
      fks=f(k)*f(k)
      fkc=fks*f(k)
      go to 240
    2 j=3
      aejs=ae(j)*ae(j)
      bejs=be(j)*be(j)
      aejc=aejs*ae(j)
      bejc=bejs*be(j)
      fjs=f(j)*f(j)
      fjc=fjs*f(j)
      k=1
      aeks=ae(k)*ae(k)
      beks=be(k)*be(k)
      aekc=aeks*ae(k)
      bekc=beks*be(k)
      fks=f(k)*f(k)
      fkc=fks*f(k)
      go to 240
    3 j=1
      aejs=ae(j)*ae(j)
      bejs=be(j)*be(j)
      aejc=aejs*ae(j)
      bejc=bejs*be(j)
      fjs=f(j)*f(j)
      fjc=fjs*f(j)
      k=2
      aeks=ae(k)*ae(k)
      beks=be(k)*be(k)
      aekc=aeks*ae(k)
      bekc=beks*be(k)
      fks=f(k)*f(k)
      fkc=fks*f(k)
      go to 240
  240 continue
      go to (11,12,13) ,i
   11 lo=12
      ko=17
      lp=12
      kp=13
      go to 241
   12 lo=15
      ko=13
      lp=15
      kp=14
      go to 241
   13 lo=16
      ko=14
      lp=16
      kp=17
      go to 241
  241 continue
! c this fxxyfyzz
      d81=(beks+0.5/at)*(aeis+0.5/at)
      d82=2.*be(k)*ae(i)*(ae(i)*f(k)+be(k)*f(i))+0.5*(beks+aeis &
     &)/at+(be(k)*f(k)+ae(i)*f(i))/at+0.5/ats
      d83=fks*aeis+fis*beks+4.*be(k)*ae(i)*f(k)*f(i) &
     &+(be(k)*f(k)+ae(i)*f(i))/at+0.5*(fks+fis)/at+0.25/ats
      d84=2.*f(k)*f(i)*(ae(i)*f(k)+be(k)*f(i))+0.5*(fks+fis &
     &)/at
      dd1=be(j)*ae(j)*d81+0.5*d81/at
      dd2=f(j)*ae(j)*d81+be(j)*(f(j)*d81+ae(j)*d82)+0.5*(d81+ &
     &d82)/at
      dd3=f(j)*(f(j)*d81+ae(j)*d82)+be(j)*(f(j)*d82+ae(j)*d83) &
     &+0.5*(d82+d83)/at
      dd4=f(j)*(f(j)*d82+ae(j)*d83)+be(j)*(f(j)*d83+ae(j)*d84) &
     &+0.5*(d83+d84)/at
      dd5=f(j)*(f(j)*d83+ae(j)*d84)+be(j)*(f(j)*d84+ae(j)*fks &
     &*fis)+0.5*(d84+fks*fis)/at
      dd6=f(j)*(f(j)*d84+ae(j)*fks*fis)+be(j)*fis &
     &*fks*f(j)+0.5*fis*fks/at
      dd7=fjs*fis*fks
      dr(lo,ko)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5- &
     &s10*dd6+s12*dd7)
! c this fyzzfxxy
      d81=(aeks+0.5/at)*(beis+0.5/at)
      d82=2.*ae(k)*be(i)*(be(i)*f(k)+ae(k)*f(i))+0.5*(aeks+be(i) &
     &**2)/at+(ae(k)*f(k)+be(i)*f(i))/at+0.5/ats
      d83=fks*beis+fis*aeks+4.*ae(k)*be(i)*f(k)*f(i) &
     &+(ae(k)*f(k)+be(i)*f(i))/at+0.5*(fks+fis)/at+0.25/ats
      d84=2.*f(k)*f(i)*(be(i)*f(k)+ae(k)*f(i))+0.5*(fks+f(i) &
     &**2)/at
      dd1=ae(j)*be(j)*d81+0.5*d81/at
      dd2=f(j)*be(j)*d81+ae(j)*(f(j)*d81+be(j)*d82)+0.5*(d81+ &
     &d82)/at
      dd3=f(j)*(f(j)*d81+be(j)*d82)+ae(j)*(f(j)*d82+be(j)*d83) &
     &+0.5*(d82+d83)/at
      dd4=f(j)*(f(j)*d82+be(j)*d83)+ae(j)*(f(j)*d83+be(j)*d84) &
     &+0.5*(d83+d84)/at
      dd5=f(j)*(f(j)*d83+be(j)*d84)+ae(j)*(f(j)*d84+be(j)*fks &
     &*fis)+0.5*(d84+fks*fis)/at
      dd6=f(j)*(f(j)*d84+be(j)*fks*fis)+ae(j)*fis &
     &*fks*f(j)+0.5*fis*fks/at
      dd7=fjs*fis*fks
      dr(ko,lo)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5- &
     &s10*dd6+s12*dd7)
! c this fxxyfzxx
      d21=aeis*beis+0.5*aeis/at+2.*ae(i)*be(i)/at+0.5 &
     &*beis/at+0.75/ats
      d22=(2.*ae(i)*be(i)*f(i)+3.*f(i)/at)*(ae(i)+be(i))+0.5* &
     &(aeis+4.*be(i)*ae(i)+beis)/at+1.5/ats
      d23=fis*(aeis+4.*ae(i)*be(i)+beis)+3.*f(i)* &
     &(ae(i)+be(i)+f(i))/at+0.75/ats
      d24=2.*fic*(ae(i)+be(i))+3.*fis/at
      dd1=be(k)*d21*ae(j)
      dd2=f(k)*ae(j)*d21+be(k)*(f(j)*d21+ae(j)*d22)
      dd3=f(k)*(f(j)*d21+ae(j)*d22)+be(k)*(f(j)*d22+ae(j)*d23)
      dd4=f(k)*(f(j)*d22+ae(j)*d23)+be(k)*(f(j)*d23+ae(j)*d24)
      dd5=f(k)*(f(j)*d23+ae(j)*d24)+be(k)*(f(j)*d24+ae(j)*f(i)**4)
      dd6=f(k)*(f(j)*d24+ae(j)*f(i)**4)+be(k)*f(j)*f(i)**4
      dd7=f(k)*f(j)*f(i)**4
      dr(lp,kp)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5- &
     &s10*dd6+s12*dd7)
! c this fzxxfxxy
      d21=beis*aeis+0.5*beis/at+2.*be(i)*ae(i)/at+0.5 &
     &*aeis/at+0.75/ats
      d22=(2.*be(i)*ae(i)*f(i)+3.*f(i)/at)*(be(i)+ae(i))+0.5* &
     &(beis+4.*ae(i)*be(i)+aeis)/at+1.5/ats
      d23=fis*(beis+4.*be(i)*ae(i)+aeis)+3.*f(i)* &
     &(be(i)+ae(i)+f(i))/at+0.75/ats
      d24=2.*fic*(be(i)+ae(i))+3.*fis/at
      dd1=ae(k)*d21*be(j)
      dd2=f(k)*be(j)*d21+ae(k)*(f(j)*d21+be(j)*d22)
      dd3=f(k)*(f(j)*d21+be(j)*d22)+ae(k)*(f(j)*d22+be(j)*d23)
      dd4=f(k)*(f(j)*d22+be(j)*d23)+ae(k)*(f(j)*d23+be(j)*d24)
      dd5=f(k)*(f(j)*d23+be(j)*d24)+ae(k)*(f(j)*d24+be(j)*f(i)**4)
      dd6=f(k)*(f(j)*d24+be(j)*f(i)**4)+ae(k)*f(j)*f(i)**4
      dd7=f(k)*f(j)*f(i)**4
      dr(kp,lp)=coef*(s0*dd1-s2*dd2+s4*dd3-s6*dd4+s8*dd5- &
     &s10*dd6+s12*dd7)
  250 continue
      return
      end subroutine pot1f2
! c
! c
      subroutine fmtcf(xx,s0,s2,s4,s6,s8,s10,s12)

      use O_Kinds
      use O_Constants

      implicit none

      real (kind=double) :: xx,s0,s2,s4,s6,s8,s10,s12

      integer :: m
      real (kind=double) :: t
      real (kind=double) :: a0,a1,a2,a3,a4,a5
      real (kind=double) :: b1,b2,b3,b4,b5,b6
      real (kind=double) :: s14
      real (kind=double), dimension (8) :: f
      t=xx
      if(t-15.0) 10,10,100
 10   continue
      do 90 m=1,7
      go to (20,30,40,50,60,70,80),m
 20   a0=1.d0
      a1=0.213271302431420d0
      a2=0.0629344460255614d0
      a3=0.00769838037756759d0
      a4=0.000758433197127160d0
      a5=0.0000564691197633667d0
      b1=0.879937801660182d0
      b2=0.338450368470103d0
      b3=0.0738522953299624d0
      b4=0.0101431553402629d0
      b5=0.000955528842975585d0
      b6=0.0000720266520392572d0
      go to 85
 30   a0=1.d0/3.d0**(2.d0/3.d0)
      a1=0.0295195994716045d0
      a2=0.0128790985465415d0
      a3=0.000998165499553218d0
      a4=0.0000970927983276419d0
      a5=0.00000493839847029699d0
      b1=0.461403194579124d0
      b2=0.108494164372449d0
      b3=0.0171462934845042d0
      b4=0.00196918657845508d0
      b5=0.000160138863265254d0
      b6=0.00000857708713007233d0
      go to 85
 40   a0=1.d0/5.d0**0.4
      a1=-0.00575763488635418d0
      a2=0.00731474973333076d0
      a3=0.00251276149443393d0
      a4=0.0000264336244559094d0
      a5=0.d0
      b1=0.274754154712841d0
      b2=0.0425364830353043d0
      b3=0.00493902790955943d0
      b4=0.000437251500927601d0
      b5=0.0000288914662393981d0
      b6=0.d0
      go to 85
 50   a0=1.d0/7.d0**(2.d0/7.d0)
      a1=-0.0290110430424666d0
      a2=0.00561884370781462d0
      a3=0.0000301628267382713d0
      a4=0.0000110671035361856d0
      a5=0.d0
      b1=0.171637608242892d0
      b2=0.0187571417256877d0
      b3=0.00178536829675118d0
      b4=0.000137360778130936d0
      b5=0.00000791915206883054d0
      b6=0.d0
      go to 85
 60   a0=1.d0/9.d0**(2.d0/9.d0)
      a1=-0.0452693111179624d0
      a2=0.00490070062899003d0
      a3=-0.0000561789719979307d0
      a4=0.00000550814626951998d0
      a5=0.d0
      b1=0.108051989937231d0
      b2=0.00855924943430755d0
      b3=0.000724968571389473d0
      b4=0.0000502338223156067d0
      b5=0.00000249107837399141d0
      b6=0.d0
      go to 85
 70   a0=1.d0/11.d0**(2.d0/11.d0)
      a1=-0.0566143259316101d0
      a2=0.00455916894577203d0
      a3=-0.0000894152721395639d0
      a4=0.00000328096732308082d0
      a5=0.d0
      b1=0.0662932958471386d0
      b2=0.00383724443872493d0
      b3=0.000327167659811839d0
      b4=0.0000210430437682548d0
      b5=0.000000883562935089333d0
      b6=0.d0
      go to 85
 80   a0=1.d0/13.d0**(2.d0/13.d0)
      a1=-0.0503249167534352d0
      a2=0.00273135625430953d0
      a3=-0.00003107336248191d0
      a4=0.d0
      a5=0.d0
      b1=0.0586609328033371d0
      b2=0.00194044691497128d0
      b3=0.000109442742502192d0
      b4=0.00000613406236401726d0
      b5=0.d0
      b6=0.d0
 85   f(m)=((a0+a1*t+a2*t*t+a3*t*t*t+a4*t*t*t*t+a5*t*t*t*t*t)/ &
     & (1.d0+b1*t+b2*t*t+b3*t*t*t+b4*t*t*t*t+b5*t*t*t*t*t+ &
     & b6*t*t*t*t*t*t))**(real(m)-0.5d0)
 90   continue
      go to 120
 100  continue
      f(1)=dsqrt(0.5d0*pi)/dsqrt(2.d0*t)
      do 110 m=1,7
      f(m+1)=f(m)*(2*m-1)/(2.d0*t)
 110  continue
 120  s0=f(1)
      s2=f(2)
      s4=f(3)
      s6=f(4)
      s8=f(5)
      s10=f(6)
      s12=f(7)
      s14=f(8)
      return
      end subroutine fmtcf
! c
!************************************************************

      subroutine threeCentInteg(g,l1,l2,al1,al2,al3,r1,r2,r3)

      use O_Kinds
      use O_Constants

      implicit none
! c
! c     g(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
! c     9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
! c     14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy
! c
! c     wo(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
! c     10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
! c     18,xxx; 19,yyy; 20,zzz
! c 
      ! define the dummy variables passed to this subroutine.
      real (kind=double), dimension (16,16) :: g
      integer :: l1, l2
      real (kind=double) :: al1
      real (kind=double) :: al2
      real (kind=double) :: al3
      real (kind=double), dimension (3) :: r1
      real (kind=double), dimension (3) :: r2
      real (kind=double), dimension (3) :: r3

      ! define local variables
      real (kind=double) :: als, ex, a
      real (kind=double), dimension (3) :: b, c
      real (kind=double), dimension (20,20) :: ov
      integer :: ndd,nff,i,j,k,l,m1,m2,m3,m,ms,mt,mm,ma,mb,na,nb,mq,mr
      integer :: lo,ko,lp,kp

      real(kind=double), dimension(3) :: ae,be
      real(kind=double) :: at,ey,ez,b2,c2,e2,ccc,coef
! c

      do 100 i=1,3
      b(i)=r2(i)-r1(i)
      c(i)=r3(i)-r1(i)
 100  continue
! c
      if(l1.eq.16.or.l2.eq.16) then
         nff=0
      else
         nff=1
      endif
      g(:,:)=0.0_double
      ov(:,:)=0.0_double

      ndd=0
      at=al1+al2+al3
      ex=(al2*b(1)+al3*c(1))/at
      ey=(al2*b(2)+al3*c(2))/at
      ez=(al2*b(3)+al3*c(3))/at
      e2=ex**2+ey**2+ez**2
      b2=b(1)**2+b(2)**2+b(3)**2
      c2=c(1)**2+c(2)**2+c(3)**2
      ccc=at*e2-al2*b2-al3*c2
      coef=dsqrt(pi/at)**3*dexp(ccc)
      ae(1)=ex
      ae(2)=ey
      ae(3)=ez
      be(1)=(ex-b(1))
      be(2)=(ey-b(2))
      be(3)=(ez-b(3))
! c this is ss
      ov(1,1)=coef
      do 105 i=1,3
! c this is spx
      ov(1,i+1)=coef*be(i)
! c this is pxs
      ov(i+1,1)=coef*ae(i)
      if(ndd.ne.0) go to 1000
! c this is sdxx
      ov(1,i+4)=coef*(be(i)**2+.5d0/at)
! c this is dxxs
      ov(i+4,1)=coef*(ae(i)**2+.5d0/at)
 1000 continue
! c this is pxpx
      ov(i+1,i+1)=coef*(ae(i)*be(i)+.5d0/at)
      if(ndd.ne.0) go to 1001
! c this is dxxpx
      ov(i+4,i+1)=coef*(be(i)*ae(i)**2+.5d0*(2.d0*ae(i)+be(i))/at)
! c this is pxdxx
      ov(i+1,i+4)=coef*(ae(i)*be(i)**2+.5d0*(2.d0*be(i)+ae(i))/at)
! c this is dxxdxx
      ov(i+4,i+4)=coef*(ae(i)**2*be(i)**2+.5d0*(ae(i)**2+be(i)**2) &
     &/at+2.d0*be(i)*ae(i)/at+.75d0/at**2)
 1001 continue
  105 continue
      if(ndd.ne.0) go to 1002
      k=7
      do 110 i=1,2
      do 110 j=2,3
      if (i.eq.j) go to 110
      k=k+1
! c this is sdxy
      ov(1,k)=coef*be(i)*be(j)
! c this is dxys
      ov(k,1)=coef*ae(i)*ae(j)
  110 continue
 1002 continue
      do 115 i=1,3
      do 115 j=1,3
      if (i.eq.j) go to 115
! c this is pxpy
      ov(i+1,j+1)=coef*ae(i)*be(j)
      if(ndd.ne.0) go to 1003
! c this is pxdyy
      ov(i+1,j+4)=coef*ae(i)*(be(j)**2+.5d0/at)
! c this is dxxpy
      ov(i+4,j+1)=coef*be(j)*(ae(i)**2+.5d0/at)
! c this is dxxdyy
      ov(i+4,j+4)=coef*(ae(i)**2+.5d0/at)*(be(j)**2+.5d0/at)
 1003 continue
  115 continue
      if(ndd.ne.0) go to 1004
      do 130 i=1,3
      do 130 j=1,2
      do 130 k=2,3
      if (j.eq.k) go to 130
      if (i.ne.j.and.i.ne.k) go to 120
      if (j.eq.i) go to 118
      l=j
      go to 119
  118 l=k
  119 continue
! c this is pxdxy
      ov(i+1,j+k+5)=coef*be(l)*(ae(i)*be(i)+.5d0/at)
! c this is dxxdxy
      ov(i+4,j+k+5)=coef*be(l)*(ae(i)**2*be(i)+.5d0*(2.d0*ae(i)+be &
     &(i))/at)
! c this is dxypx
      ov(j+k+5,i+1)=coef*ae(l)*(ae(i)*be(i)+.5d0/at)
! c this is dxydxx
      ov(j+k+5,i+4)=coef*ae(l)*(ae(i)*be(i)**2+.5d0*(2.d0*be(i)+ae &
     &(i))/at)
      go to 130
  120 continue
! c this is pxdyz
      ov(i+1,j+k+5)=coef*be(j)*be(k)*ae(i)
! c this is dxxdyz
      ov(i+4,j+k+5)=coef*be(j)*be(k)*(ae(i)**2+.5d0/at)
! c this is dyzpx
      ov(j+k+5,i+1)=coef*ae(j)*ae(k)*be(i)
! c this is dyzdxx
      ov(j+k+5,i+4)=coef*ae(j)*ae(k)*(be(i)**2+.5d0/at)
  130 continue
      do 146 i=1,2
      do 146 j=2,3
      if (i.eq.j) go to 146
      do 145 k=1,2
      do 145 l=2,3
      if (k.eq.l) go to 145
      if (i.eq.k.and.j.eq.l) go to 140
      if (i.eq.k) go to 131
      if (i.eq.l) go to 132
      if (j.eq.k) go to 133
      if  (j.eq.l) go to 134
  131 m1=i
      m2=j
      m3=l
      go to 135
  132 m1=i
      m2=j
      m3=k
      go to 135
  133 m1=j
      m2=i
      m3=l
      go to 135
  134 m1=j
      m2=i
      m3=k
      go to 135
  135 continue
! c  this dxy dyz
      ov(i+j+5,k+l+5)=coef*ae(m2)*be(m3)*(ae(m1)*be(m1)+.5/at)
      go to 145
  140 continue
! c  this dxy dxy
      ov(i+j+5,k+l+5)=coef*(ae(i)*be(i)+.5/at)*(ae(j)*be(j)+.5/at)
  145 continue
  146 continue
 1004 continue
      if(nff.eq.1) go to 999
      do 205 i=1,3
! c   this loop has exchanging in x-y-z
! c   this fxxx s
      ov(i+17,1)=coef*ae(i)*(1.5/at+ae(i)**2)
! c   this s fxxx
      ov(1,i+17)=coef*be(i)*(1.5/at+be(i)**2)
! c   this fxxx px
      ov(i+17,i+1)=coef*(1.5*ae(i)**2/at+1.5*ae(i)*be(i) &
     &   /at+ae(i)**3*be(i)+.75/at**2)
! c   this px fxxx
      ov(i+1,i+17)=coef*(1.5*be(i)**2/at+1.5*be(i)*ae(i) &
     &   /at+be(i)**3*ae(i)+.75/at**2)
! c   this fxxx dxx
      ov(i+17,i+4)=coef*(9.*ae(i)/(4.*at**2)+3.*be(i)/(2.*at**2) &
     & +3.*ae(i)**2*be(i)/at+1.5*ae(i)*be(i)**2/at+.5*ae(i)**3 &
     & /at+be(i)**2*ae(i)**3)
! c   this dxx fxxx
      ov(i+4,i+17)=coef*(9.*be(i)/(4.*at**2)+3.*ae(i)/(2.*at**2) &
     &  +3.*be(i)**2*ae(i)/at+1.5*be(i)*ae(i)**2/at+.5*be(i) &
     &   **3/at+ae(i)**2*be(i)**3)
! c   this fxxx fxxx
      ov(i+17,i+17)=coef*((be(i)*27.*ae(i)+9.*be(i)**2)/(4.*at**2) &
     & +be(i)**2*(4.5*ae(i)**2+1.5*ae(i)*be(i))/at+ae(i)**3* &
     & be(i)*(1.5/at+be(i)**2)+9.*ae(i)**2/(4.*at**2)+15./(8.*at**3))
  205 continue
      l=11
      do 210 i=1,3
      do 210 j=1,3
      if (i.eq.j) go to 210
! c  this loop has exchanging in x-y-z,x-y,y-z,z-x,and six symmetry
! c  this fxxy s
      l=l+1
      ov(l,1)=coef*ae(j)*(ae(i)**2+.5/at)
! c  this s fxxy
      ov(1,l)=coef*be(j)*(be(i)**2+.5/at)
! c  this fxxy px
      ov(l,i+1)=coef*ae(j)*(.5*(2.*ae(i)+be(i))/at+be(i)*ae(i)**2)
! c  this px fxxy
      ov(i+1,l)=coef*be(j)*(.5*(2.*be(i)+ae(i))/at+ae(i)*be(i)**2)
! c  this fxxy py
      ov(l,j+1)=coef*(ae(i)**2+.5/at)*(ae(j)*be(j)+.5/at)
! c  this py fxxy
      ov(j+1,l)=coef*(be(i)**2+.5/at)*(be(j)*ae(j)+.5/at)
! c  this fxxy dxx
      ov(l,i+4)=coef*ae(j)*(.75/at**2+2.*ae(i)*be(i)/at &
     &  +.5*(ae(i)**2+be(i)**2)/at+ae(i)**2*be(i)**2)
! c  this dxx fxxy
      ov(i+4,l)=coef*be(j)*(.75/at**2+2.*be(i)*ae(i)/at &
     &  +.5*(be(i)**2+ae(i)**2)/at+be(i)**2*ae(i)**2)
! c  this fxxy dyy
      ov(l,j+4)=coef*(.5/at+ae(i)**2)*(.5*ae(j)/at+ &
     &  be(j)/at+ae(j)*be(j)**2)
! c  this dyy fxxy
      ov(j+4,l)=coef*(.5/at+be(i)**2)*(.5*be(j)/at+ &
     &  ae(j)/at+be(j)*ae(j)**2)
! c  this fxxy dxy
      if (i.ne.3.and.j.ne.3) m=8
      if (i.ne.1.and.j.ne.1) m=10
      if (i.ne.2.and.j.ne.2) m=9
      ov(l,m)=coef*(ae(j)*be(j)*(ae(i)/at+.5*be(i)/at+be(i)* &
     &  ae(i)**2)+(.5*ae(i)+.25*be(i))/at**2+.5*be(i)*ae(i)**2/at)
! c  this dxy fxxy
      ov(m,l)=coef*(be(j)*ae(j)*(be(i)/at+.5*ae(i)/at+ae(i)* &
     &  be(i)**2)+(.5*be(i)+.25*ae(i))/at**2+.5*ae(i)*be(i)**2/at)
! c  this fxxy fxxy
      ov(l,l)=coef*(.5/at+ae(j)*be(j))*(.75/at**2+2.*ae(i)*be(i)/ &
     &  at+.5*(be(i)**2+ae(i)**2)/at+ae(i)**2*be(i)**2)
! c  this fxxx py
      ov(i+17,j+1)=coef*ae(i)*be(j)*(ae(i)**2+1.5/at)
! c  this py fxxx
      ov(j+1,i+17)=coef*be(i)*ae(j)*(be(i)**2+1.5/at)
! c  this fxxx dyy
      ov(i+17,j+4)=coef*ae(i)*(.5/at+be(j)**2)*(1.5/at+ae(i)**2)
! c  this dyy fxxx
      ov(j+4,i+17)=coef*be(i)*(.5/at+ae(j)**2)*(1.5/at+be(i)**2)
! c  this fxxx dxy
      if (i.ne.3.and.j.ne.3) ms=8
      if (i.ne.2.and.j.ne.2) ms=9
      if (i.ne.1.and.j.ne.1) ms=10
      ov(i+17,ms)=coef*be(j)*(1.5*(ae(i)**2+ae(i)*be(i))/at+ &
     &  ae(i)**3*be(i)+.75/at**2)
! c  this dxy fxxx
      ov(ms,i+17)=coef*ae(j)*(1.5*(be(i)**2+be(i)*ae(i))/at+ &
     &  be(i)**3*ae(i)+.75/at**2)
! c  this fxxx fxxy
      if (i.eq.1.and.j.eq.2) mt=12
      if (i.eq.1.and.j.eq.3) mt=13
      if (i.eq.2.and.j.eq.1) mt=14
      if (i.eq.2.and.j.eq.3) mt=15
      if (i.eq.3.and.j.eq.1) mt=16
      if (i.eq.3.and.j.eq.2) mt=17
      ov(i+17,mt)=coef*be(j)*(9.*ae(i)/(4.*at**2)+1.5*be(i)/at**2 &
     &  +3.*ae(i)**2*be(i)/at+1.5*ae(i)*be(i)**2/at+.5*ae(i)**3/at &
     &  +be(i)**2*ae(i)**3)
! c  this fxxy fxxx
      ov(mt,i+17)=coef*ae(j)*(9.*be(i)/(4.*at**2)+1.5*ae(i)/at**2 &
     &  +3.*be(i)**2*ae(i)/at+1.5*be(i)*ae(i)**2/at+.5*be(i)**3/at &
     &  +ae(i)**2*be(i)**3)
! c  this fxxx fxyy
      if(i.eq.1.and.j.eq.2) mm=14
      if(i.eq.1.and.j.eq.3) mm=16
      if(i.eq.2.and.j.eq.1) mm=12
      if(i.eq.2.and.j.eq.3) mm=17
      if(i.eq.3.and.j.eq.1) mm=13
      if(i.eq.3.and.j.eq.2) mm=15
      ov(i+17,mm)=coef*(.5/at+be(j)**2)*((1.5/at+ae(i)**2)*(.5/at &
     &  +be(i)*ae(i))+ae(i)**2/at)
! c  this fxyy fxxx
      ov(mm,i+17)=coef*(.5/at+ae(j)**2)*((1.5/at+be(i)**2)*(.5/at &
     &  +ae(i)*be(i))+be(i)**2/at)
  210 continue
! c this loop has three symmentry elements for x,y,z
      do 215 i=1,2
      do 215 j=2,3
      if (i.eq.j) go to 215
      k=i+j-2
      go to (21,22,23), k
   21 ma=12
      mb=14
      na=18
      nb=19
      go to 141
   22 ma=13
      mb=16
      na=18
      nb=20
      go to 141
   23 ma=15
      mb=17
      na=19
      nb=20
      go to 141
  141 continue
! c  this fxxy fxyy
      ov(ma,mb)=coef*((.5*ae(j)+be(j))/at+ae(j)*be(j)**2)*(be(i)* &
     &  (.5/at+ae(i)**2)+ae(i)/at)
! c  this fxyy fxxy
      ov(mb,ma)=coef*((.5*be(j)+ae(j))/at+be(j)*ae(j)**2)*(ae(i)* &
     &  (.5/at+be(i)**2)+be(i)/at)
! c  this fxxx fyyy
      ov(na,nb)=coef*ae(i)*be(j)*(1.5/at+ae(i)**2)*(1.5/at+be(j)**2)
! c  this fyyy fxxx
      ov(nb,na)=coef*be(i)*ae(j)*(1.5/at+be(i)**2)*(1.5/at+ae(j)**2)
  215 continue
      i=1
      j=2
      k=3
! c  this loop only for xyz
! c   this fxyz s
      ov(11,1)=coef*ae(i)*ae(j)*ae(k)
! c   this s fxyz
      ov(1,11)=coef*be(i)*be(j)*be(k)
! c   this fxyz fxyz
      ov(11,11)=coef*(.5/at+ae(j)*be(j))*(.5/at+ae(k)*be(k)) &
     &  *(.5/at+ae(i)*be(i))
! c     continue
! c  this loop has exchange in x-y-z,y-x,z-x and three symmetry
      do 220 i=1,3
      do 220 j=1,2
      do 220 k=2,3
      if (i.ne.j.and.j.ne.k.and.k.ne.i) go to 221
      go to 220
  221 continue
! c   this fxyz px
      ov(11,i+1)=coef*ae(k)*ae(j)*(.5/at+ae(i)*be(i))
! c   this px fxyz
      ov(i+1,11)=coef*be(k)*be(j)*(.5/at+be(i)*ae(i))
! c   this fxyz dxx
      ov(11,i+4)=coef*ae(k)*ae(j)*(.5*(2.*be(i)+ae(i))/at+ae(i)* &
     &   be(i)**2)
! c   this dxx fxyz
      ov(i+4,11)=coef*be(k)*be(j)*(.5*(2.*ae(i)+be(i))/at+be(i)* &
     &   ae(i)**2)
! c   this fxyz dyz
      ov(11,11-i)=coef*ae(i)*(.5/at+ae(j)*be(j))*(.5/at+ae(k)*be(k))
! c   this dyz fxyz
      ov(11-i,11)=coef*be(i)*(.5/at+be(k)*ae(k))*(.5/at+be(j)*ae(j))
! c   this fxyz fxxx
      ov(11,i+17)=coef*ae(j)*ae(k)*((be(i)**2+1.5/at)*(ae(i)*be(i) &
     &  +.5/at)+be(i)**2/at)
! c   this fxxx fxyz
      ov(i+17,11)=coef*be(j)*be(k)*((ae(i)**2+1.5/at)*(be(i)*ae(i) &
     &  +.5/at)+ae(i)**2/at)
! c   this fxxx dyz
      ov(i+17,11-i)=coef*ae(i)*be(j)*be(k)*(ae(i)**2+1.5/at)
! c   this dyz fxxx
      ov(11-i,i+17)=coef*be(i)*ae(j)*ae(k)*(be(i)**2+1.5/at)
  220 continue
! c  this loop has exchange x-y-z and x-y,y-z,z-x and six symmentry
      l=11
      do 230 i=1,3
      do 230 j=1,3
      do 230 k=1,3
      if (i.ne.j.and.j.ne.k.and.k.ne.i) go to 231
      go to 230
  231 continue
      l=l+1
! c   this fxxy pz
      ov(l,k+1)=coef*ae(j)*be(k)*(.5/at+ae(i)**2)
! c   tihs pz fxxy
      ov(k+1,l)=coef*be(j)*ae(k)*(.5/at+be(i)**2)
! c   this fxxy dyz
      ov(l,11-i)=coef*be(k)*(.5/at+ae(i)**2)*(ae(j)*be(j)+.5/at)
! c   this dyz fxxy
      ov(11-i,l)=coef*ae(k)*(.5/at+be(i)**2)*(be(j)*ae(j)+.5/at)
! c   this fxxy dzx
      if (j.eq.3) mq=8
      if (j.eq.1) mq=10
      if (j.eq.2) mq=9
      ov(l,mq)=coef*ae(j)*be(k)*(.5*be(i)/at+ae(i)/at+be(i)*ae(i)**2)
! c   this dzx fxxy
      ov(mq,l)=coef*be(j)*ae(k)*(.5*ae(i)/at+be(i)/at+ae(i)*be(i)**2)
! c   this fxxy dzz
      ov(l,k+4)=coef*ae(j)*(.5/at+be(k)**2)*(.5/at+ae(i)**2)
! c   this dzz fxxy
      ov(k+4,l)=coef*be(j)*(.5/at+ae(k)**2)*(.5/at+be(i)**2)
! c   this fxyz fxxy
      ov(11,l)=coef*ae(k)*(.5/at+be(j)*ae(j))*(ae(i)*(.5/at+be(i)**2) &
     &   +be(i)/at)
! c   this fxxy fxyz
      ov(l,11)=coef*be(k)*(.5/at+ae(j)*be(j))*(be(i)*(.5/at+ae(i)**2) &
     &   +ae(i)/at)
! c   this fxxy fyyz
      if(i.eq.1.and.j.eq.2) mr=15
      if(i.eq.1.and.j.eq.3) mr=17
      if(i.eq.2.and.j.eq.1) mr=13
      if(i.eq.2.and.j.eq.3) mr=16
      if(i.eq.3.and.j.eq.1) mr=12
      if(i.eq.3.and.j.eq.2) mr=14
      ov(l,mr)=coef*be(k)*(.5/at+ae(i)**2)*(.5*ae(j)/at+be(j)/at &
     &  +ae(j)*be(j)**2)
! c   this fyyz fxxy
      ov(mr,l)=coef*ae(k)*(.5/at+be(i)**2)*(.5*be(j)/at+ae(j)/at &
     &  +be(j)*ae(j)**2)
! c   this fxxy fzzz
      ov(l,k+17)=coef*ae(j)*be(k)*(.5/at+ae(i)**2)*(1.5/at+be(k)**2)
! c   this fzzz fxxy
      ov(k+17,l)=coef*be(j)*ae(k)*(.5/at+be(i)**2)*(1.5/at+ae(k)**2)
  230 continue
! c  this loop has three symmentry elements in x,y,z
      do 250 i=1,3
      go to (1,2,3),  i
    1 j=2
      k=3
      go to 240
    2 j=3
      k=1
      go to 240
    3 j=1
      k=2
      go to 240
  240 continue
      go to (11,12,13), i
   11 lo=12
      ko=17
      lp=12
      kp=13
      go to 241
   12 lo=15
      ko=13
      lp=15
      kp=14
      go to 241
   13 lo=16
      ko=14
      lp=16
      kp=17
      go to 241
  241 continue
! c   this fxxy fyzz
      ov(lo,ko)=coef*(.5/at+be(k)**2)*(.5/at+ae(i)**2)*(.5/at+ae(j)*be &
     &     (j))
! c   this fyzz fxxy
      ov(ko,lo)=coef*(.5/at+ae(k)**2)*(.5/at+be(i)**2)*(.5/at+be(j) &
     &     *ae(j))
! c   this fxxy fzxx
      ov(lp,kp)=coef*be(k)*ae(j)*(.75/at**2+2.*ae(i)*be(i)/at+.5* &
     &  (be(i)**2+ae(i)**2)/at+ae(i)**2*be(i)**2)
! c   this fzxx fxxy
      ov(kp,lp)=coef*ae(k)*be(j)*(.75/at**2+2.*be(i)*ae(i)/at+.5* &
     &  (ae(i)**2+be(i)**2)/at+be(i)**2*ae(i)**2)
  250 continue
  999 continue
      continue
      g(1,1)=ov(1,1)
      if(l2.eq.1) goto 2009
      g(1,2)=ov(1,2)
      g(1,3)=ov(1,3)
      g(1,4)=ov(1,4)
      if(l2.eq.4) goto 2009
      g(1,5)=ov(1,8)
      g(1,6)=ov(1,9)
      g(1,7)=ov(1,10)
      g(1,8)=ov(1,5)-ov(1,6)
      g(1,9)=2.0*ov(1,7)-ov(1,5)-ov(1,6)
      if(l2.eq.9) goto 2009
      g(1,10)=ov(1,11)
      g(1,11)=ov(1,13)-ov(1,15)
      g(1,12)=ov(1,18)-3.0*ov(1,14)
      g(1,13)=3.0*ov(1,12)-ov(1,19)
      g(1,14)=2.0*ov(1,20)-3.0*ov(1,13)-3.0*ov(1,15)
      g(1,15)=4.0*ov(1,16)-ov(1,18)-ov(1,14)
      g(1,16)=4.0*ov(1,17)-ov(1,12)-ov(1,19)
! c
 2009  continue
      if(l1.eq.1) return
      g(2,1)=ov(2,1)
      if(l2.eq.1) goto 3009
      g(2,2)=ov(2,2)
      g(2,3)=ov(2,3)
      g(2,4)=ov(2,4)
      if(l2.eq.4) goto 3009
      g(2,5)=ov(2,8)
      g(2,6)=ov(2,9)
      g(2,7)=ov(2,10)
      g(2,8)=ov(2,5)-ov(2,6)
      g(2,9)=2.0*ov(2,7)-ov(2,5)-ov(2,6)
      if(l2.eq.9) goto 3009
      g(2,10)=ov(2,11)
      g(2,11)=ov(2,13)-ov(2,15)
      g(2,12)=ov(2,18)-3.0*ov(2,14)
      g(2,13)=3.0*ov(2,12)-ov(2,19)
      g(2,14)=2.0*ov(2,20)-3.0*ov(2,13)-3.0*ov(2,15)
      g(2,15)=4.0*ov(2,16)-ov(2,18)-ov(2,14)
      g(2,16)=4.0*ov(2,17)-ov(2,12)-ov(2,19)
! c
 3009  continue
      g(3,1)=ov(3,1)
      if(l2.eq.1) goto 4009
      g(3,2)=ov(3,2)
      g(3,3)=ov(3,3)
      g(3,4)=ov(3,4)
      if(l2.eq.4) goto 4009
      g(3,5)=ov(3,8)
      g(3,6)=ov(3,9)
      g(3,7)=ov(3,10)
      g(3,8)=ov(3,5)-ov(3,6)
      g(3,9)=2.0*ov(3,7)-ov(3,5)-ov(3,6)
      if(l2.eq.9) goto 4009
      g(3,10)=ov(3,11)
      g(3,11)=ov(3,13)-ov(3,15)
      g(3,12)=ov(3,18)-3.0*ov(3,14)
      g(3,13)=3.0*ov(3,12)-ov(3,19)
      g(3,14)=2.0*ov(3,20)-3.0*ov(3,13)-3.0*ov(3,15)
      g(3,15)=4.0*ov(3,16)-ov(3,18)-ov(3,14)
      g(3,16)=4.0*ov(3,17)-ov(3,12)-ov(3,19)
! c
 4009  continue
      g(4,1)=ov(4,1)
      if(l2.eq.1) goto 5009
      g(4,2)=ov(4,2)
      g(4,3)=ov(4,3)
      g(4,4)=ov(4,4)
      if(l2.eq.4) goto 5009
      g(4,5)=ov(4,8)
      g(4,6)=ov(4,9)
      g(4,7)=ov(4,10)
      g(4,8)=ov(4,5)-ov(4,6)
      g(4,9)=2.0*ov(4,7)-ov(4,5)-ov(4,6)
      if(l2.eq.9) goto 5009
      g(4,10)=ov(4,11)
      g(4,11)=ov(4,13)-ov(4,15)
      g(4,12)=ov(4,18)-3.0*ov(4,14)
      g(4,13)=3.0*ov(4,12)-ov(4,19)
      g(4,14)=2.0*ov(4,20)-3.0*ov(4,13)-3.0*ov(4,15)
      g(4,15)=4.0*ov(4,16)-ov(4,18)-ov(4,14)
      g(4,16)=4.0*ov(4,17)-ov(4,12)-ov(4,19)
! c
 5009  continue
      if(l1.eq.4) return
      g(5,1)=ov(8,1)
      if(l2.eq.1) goto 6009
      g(5,2)=ov(8,2)
      g(5,3)=ov(8,3)
      g(5,4)=ov(8,4)
      if(l2.eq.4) goto 6009
      g(5,5)=ov(8,8)
      g(5,6)=ov(8,9)
      g(5,7)=ov(8,10)
      g(5,8)=ov(8,5)-ov(8,6)
      g(5,9)=2.0*ov(8,7)-ov(8,5)-ov(8,6)
      if(l2.eq.9) goto 6009
      g(5,10)=ov(8,11)
      g(5,11)=ov(8,13)-ov(8,15)
      g(5,12)=ov(8,18)-3.0*ov(8,14)
      g(5,13)=3.0*ov(8,12)-ov(8,19)
      g(5,14)=2.0*ov(8,20)-3.0*ov(8,13)-3.0*ov(8,15)
      g(5,15)=4.0*ov(8,16)-ov(8,18)-ov(8,14)
      g(5,16)=4.0*ov(8,17)-ov(8,12)-ov(8,19)
! c
 6009  continue
      g(6,1)=ov(9,1)
      if(l2.eq.1) goto 7009
      g(6,2)=ov(9,2)
      g(6,3)=ov(9,3)
      g(6,4)=ov(9,4)
      if(l2.eq.4) goto 7009
      g(6,5)=ov(9,8)
      g(6,6)=ov(9,9)
      g(6,7)=ov(9,10)
      g(6,8)=ov(9,5)-ov(9,6)
      g(6,9)=2.0*ov(9,7)-ov(9,5)-ov(9,6)
      if(l2.eq.9) goto 7009
      g(6,10)=ov(9,11)
      g(6,11)=ov(9,13)-ov(9,15)
      g(6,12)=ov(9,18)-3.0*ov(9,14)
      g(6,13)=3.0*ov(9,12)-ov(9,19)
      g(6,14)=2.0*ov(9,20)-3.0*ov(9,13)-3.0*ov(9,15)
      g(6,15)=4.0*ov(9,16)-ov(9,18)-ov(9,14)
      g(6,16)=4.0*ov(9,17)-ov(9,12)-ov(9,19)
! c
 7009  continue
      g(7,1)=ov(10,1)
      if(l2.eq.1) goto 8009
      g(7,2)=ov(10,2)
      g(7,3)=ov(10,3)
      g(7,4)=ov(10,4)
      if(l2.eq.4) goto 8009
      g(7,5)=ov(10,8)
      g(7,6)=ov(10,9)
      g(7,7)=ov(10,10)
      g(7,8)=ov(10,5)-ov(10,6)
      g(7,9)=2.0*ov(10,7)-ov(10,5)-ov(10,6)
      if(l2.eq.9) goto 8009
      g(7,10)=ov(10,11)
      g(7,11)=ov(10,13)-ov(10,15)
      g(7,12)=ov(10,18)-3.0*ov(10,14)
      g(7,13)=3.0*ov(10,12)-ov(10,19)
      g(7,14)=2.0*ov(10,20)-3.0*ov(10,13)-3.0*ov(10,15)
      g(7,15)=4.0*ov(10,16)-ov(10,18)-ov(10,14)
      g(7,16)=4.0*ov(10,17)-ov(10,12)-ov(10,19)
! c
 8009  continue
      g(8,1)=ov(5,1)-ov(6,1)
      if(l2.eq.1) goto 9009
      g(8,2)=ov(5,2)-ov(6,2)
      g(8,3)=ov(5,3)-ov(6,3)
      g(8,4)=ov(5,4)-ov(6,4)
      if(l2.eq.4) goto 9009
      g(8,5)=ov(5,8)-ov(6,8)
      g(8,6)=ov(5,9)-ov(6,9)
      g(8,7)=ov(5,10)-ov(6,10)
      g(8,8)=ov(5,5)-ov(6,5)-ov(5,6)+ov(6,6)
      g(8,9)=2.0*ov(5,7)-2.0*ov(6,7)-ov(5,5) &
     & +ov(6,5)-ov(5,6)+ov(6,6)
      if(l2.eq.9) goto 9009
      g(8,10)=ov(5,11)-ov(6,11)
      g(8,11)=ov(5,13)-ov(6,13)-ov(5,15)+ov(6,15)
      g(8,12)=ov(5,18)-ov(6,18)-3.0*ov(5,14) &
     & +3.0*ov(6,14)
      g(8,13)=3.0*ov(5,12)-3.0*ov(6,12)-ov(5,19)+ov(6,19)
      g(8,14)=2.0*ov(5,20)-2.0*ov(6,20)-3.0*ov(5,13) &
     & +3.0*ov(6,13)-3.0*ov(5,15)+3.0*ov(6,15)
      g(8,15)=4.0*ov(5,16)-4.0*ov(6,16)-ov(5,18)+ov(6,18) &
     & -ov(5,14)+ov(6,14)
      g(8,16)=4.0*ov(5,17)-4.0*ov(6,17)-ov(5,12)+ov(6,12) &
     & -ov(5,19)+ov(6,19)
! c
 9009  continue
      g(9,1)=2.0*ov(7,1)-ov(5,1)-ov(6,1)
      if(l2.eq.1) goto 10009
      g(9,2)=2.0*ov(7,2)-ov(5,2)-ov(6,2)
      g(9,3)=2.0*ov(7,3)-ov(5,3)-ov(6,3)
      g(9,4)=2.0*ov(7,4)-ov(5,4)-ov(6,4)
      if(l2.eq.4) goto 10009
      g(9,5)=2.0*ov(7,8)-ov(5,8)-ov(6,8)
      g(9,6)=2.0*ov(7,9)-ov(5,9)-ov(6,9)
      g(9,7)=2.0*ov(7,10)-ov(5,10)-ov(6,10)
      g(9,8)=2.0*ov(7,5)-ov(5,5)-ov(6,5)-2.0*ov(7,6)+ov(5,6)+ov(6,6)
      g(9,9)=4.0*ov(7,7)-2.0*ov(5,7)-2.0*ov(6,7)-2.0*ov(7,5)+ov(5,5) &
     & +ov(6,5)-2.0*ov(7,6)+ov(5,6)+ov(6,6)
      if(l2.eq.9) goto 10009
      g(9,10)=2.0*ov(7,11)-ov(5,11)-ov(6,11)
      g(9,11)=2.0*ov(7,13)-ov(5,13)-ov(6,13)-2.0*ov(7,15) &
     & +ov(5,15)+ov(6,15)
      g(9,12)=2.0*ov(7,18)-ov(5,18)-ov(6,18)-6.0*ov(7,14)+ &
     & 3.0*ov(5,14)+3.0*ov(6,14)
      g(9,13)=6.0*ov(7,12)-3.0*ov(5,12)-3.0*ov(6,12) &
     & -2.0*ov(7,19)+ov(5,19)+ov(6,19)
      g(9,14)=4.0*ov(7,20)-2.0*ov(5,20)-2.0*ov(6,20)-6.0*ov(7,13)+ &
     & 3.0*ov(5,13)+3.0*ov(6,13)-6.0*ov(7,15)+3.0*ov(5,15)+3.0*ov(6,15)
      g(9,15)=8.0*ov(7,16)-4.0*ov(5,16)-4.0*ov(6,16)-2.0*ov(7,18) &
     & +ov(5,18)+ov(6,18)-2.0*ov(7,14)+ov(5,14)+ov(6,14)
      g(9,16)=8.0*ov(7,17)-4.0*ov(5,17)-4.0*ov(6,17)-2.0*ov(7,12) &
     & +ov(5,12)+ov(6,12)-2.0*ov(7,19)+ov(5,19)+ov(6,19)
! c
 10009 continue
      if(l1.eq.9) return
      g(10,1)=ov(11,1)
      if(l2.eq.1) goto 11009
      g(10,2)=ov(11,2)
      g(10,3)=ov(11,3)
      g(10,4)=ov(11,4)
      if(l2.eq.4) goto 11009
      g(10,5)=ov(11,8)
      g(10,6)=ov(11,9)
      g(10,7)=ov(11,10)
      g(10,8)=ov(11,5)-ov(11,6)
      g(10,9)=2.0*ov(11,7)-ov(11,5)-ov(11,6)
      if(l2.eq.9) goto 11009
      g(10,10)=ov(11,11)
      g(10,11)=ov(11,13)-ov(11,15)
      g(10,12)=ov(11,18)-3.0*ov(11,14)
      g(10,13)=3.0*ov(11,12)-ov(11,19)
      g(10,14)=2.0*ov(11,20)-3.0*ov(11,13)-3.0*ov(11,15)
      g(10,15)=4.0*ov(11,16)-ov(11,18)-ov(11,14)
      g(10,16)=4.0*ov(11,17)-ov(11,12)-ov(11,19)
! c
 11009 continue
      g(11,1)=ov(13,1)-ov(15,1)
      if(l2.eq.1) goto 12009
      g(11,2)=ov(13,2)-ov(15,2)
      g(11,3)=ov(13,3)-ov(15,3)
      g(11,4)=ov(13,4)-ov(15,4)
      if(l2.eq.4) goto 12009
      g(11,5)=ov(13,8)-ov(15,8)
      g(11,6)=ov(13,9)-ov(15,9)
      g(11,7)=ov(13,10)-ov(15,10)
      g(11,8)=ov(13,5)-ov(15,5)-ov(13,6)+ov(15,6)
      g(11,9)=2.0*ov(13,7)-2.0*ov(15,7)-ov(13,5) &
     & +ov(15,5)-ov(13,6)+ov(15,6)
      if(l2.eq.9) goto 12009
      g(11,10)=ov(13,11)-ov(15,11)
      g(11,11)=ov(13,13)-ov(15,13)-ov(13,15)+ov(15,15)
      g(11,12)=ov(13,18)-ov(15,18)-3.0*ov(13,14) &
     & +3.0*ov(15,14)
      g(11,13)=3.0*ov(13,12)-3.0*ov(15,12)- &
     & ov(13,19)+ov(15,19)
      g(11,14)=2.0*ov(13,20)-2.0*ov(15,20)-3.0*ov(13,13) &
     & +3.0*ov(15,13)-3.0*ov(13,15)+3.0*ov(15,15)
      g(11,15)=4.0*ov(13,16)-4.0*ov(15,16)-ov(13,18)+ov(15,18) &
     & -ov(13,14)+ov(15,14)
      g(11,16)=4.0*ov(13,17)-4.0*ov(15,17)-ov(13,12)+ov(15,12) &
     & -ov(13,19)+ov(15,19)
! c
 12009 continue
      g(12,1)=ov(18,1)-3.0*ov(14,1)
      if(l2.eq.1) goto 13009
      g(12,2)=ov(18,2)-3.0*ov(14,2)
      g(12,3)=ov(18,3)-3.0*ov(14,3)
      g(12,4)=ov(18,4)-3.0*ov(14,4)
      if(l2.eq.4) goto 13009
      g(12,5)=ov(18,8)-3.0*ov(14,8)
      g(12,6)=ov(18,9)-3.0*ov(14,9)
      g(12,7)=ov(18,10)-3.0*ov(14,10)
      g(12,8)=ov(18,5)-3.0*ov(14,5) &
     & -ov(18,6)+3.0*ov(14,6)
      g(12,9)=2.0*ov(18,7)-6.0*ov(14,7) &
     & -ov(18,5)+3.0*ov(14,5) &
     & -ov(18,6)+3.0*ov(14,6)
      if(l2.eq.9) goto 13009
      g(12,10)=ov(18,11)-3.0*ov(14,11)
      g(12,11)=ov(18,13)-3.0*ov(14,13) &
     & -ov(18,15)+3.0*ov(14,15)
      g(12,12)=ov(18,18)-3.0*ov(14,18) &
     & -3.0*ov(18,14)+9.0*ov(14,14)
      g(12,13)=3.0*ov(18,12)-9.0*ov(14,12) &
     & -ov(18,19)+3.0*ov(14,19)
      g(12,14)=2.0*ov(18,20)-6.0*ov(14,20)-3.0*ov(18,13) &
     & +9.0*ov(14,13)-3.0*ov(18,15)+9.0*ov(14,15)
      g(12,15)=4.0*ov(18,16)-12.0*ov(14,16)-ov(18,18) &
     & +3.0*ov(14,18)-ov(18,14)+3.0*ov(14,14)
      g(12,16)=4.0*ov(18,17)-12.0*ov(14,17)-ov(18,12) &
     & +3.0*ov(14,12)-ov(18,19)+3.0*ov(14,19)
! c
 13009 continue
      g(13,1)=3.0*ov(12,1)-ov(19,1)
      if(l2.eq.1) goto 14009
      g(13,2)=3.0*ov(12,2)-ov(19,2)
      g(13,3)=3.0*ov(12,3)-ov(19,3)
      g(13,4)=3.0*ov(12,4)-ov(19,4)
      if(l2.eq.4) goto 14009
      g(13,5)=3.0*ov(12,8)-ov(19,8)
      g(13,6)=3.0*ov(12,9)-ov(19,9)
      g(13,7)=3.0*ov(12,10)-ov(19,10)
      g(13,8)=3.0*ov(12,5)-ov(19,5)-3.0*ov(12,6)+ov(19,6)
      g(13,9)=6.0*ov(12,7)-2.0*ov(19,7)-3.0*ov(12,5) &
     & +ov(19,5)-3.0*ov(12,6)+ov(19,6)
      if(l2.eq.9) goto 14009
      g(13,10)=3.0*ov(12,11)-ov(19,11)
      g(13,11)=3.0*ov(12,13)-ov(19,13)-3.0*ov(12,15)+ov(19,15)
      g(13,12)=3.0*ov(12,18)-ov(19,18)-9.0*ov(12,14)+3.0*ov(19,14)
      g(13,13)=9.0*ov(12,12)-3.0*ov(19,12)-3.0*ov(12,19)+ov(19,19)
      g(13,14)=6.0*ov(12,20)-2.0*ov(19,20)-9.0*ov(12,13) &
     & +3.0*ov(19,13)-9.0*ov(12,15)+3.0*ov(19,15)
      g(13,15)=12.0*ov(12,16)-4.0*ov(19,16)-3.0*ov(12,18) &
     & +ov(19,18)-3.0*ov(12,14)+ov(19,14)
      g(13,16)=12.0*ov(12,17)-4.0*ov(19,17)-3.0*ov(12,12) &
     & +ov(19,12)-3.0*ov(12,19)+ov(19,19)
! c
 14009 continue
      g(14,1)=2.0*ov(20,1)-3.0*ov(13,1)-3.0*ov(15,1)
      if(l2.eq.1) goto 15009
      g(14,2)=2.0*ov(20,2)-3.0*ov(13,2)-3.0*ov(15,2)
      g(14,3)=2.0*ov(20,3)-3.0*ov(13,3)-3.0*ov(15,3)
      g(14,4)=2.0*ov(20,4)-3.0*ov(13,4)-3.0*ov(15,4)
      if(l2.eq.4) goto 15009
      g(14,5)=2.0*ov(20,8)-3.0*ov(13,8)-3.0*ov(15,8)
      g(14,6)=2.0*ov(20,9)-3.0*ov(13,9)-3.0*ov(15,9)
      g(14,7)=2.0*ov(20,10)-3.0*ov(13,10)-3.0*ov(15,10)
      g(14,8)=2.0*ov(20,5)-3.0*ov(13,5)-3.0*ov(15,5) &
     & -2.0*ov(20,6)+3.0*ov(13,6)+3.0*ov(15,6)
      g(14,9)=4.0*ov(20,7)-6.0*ov(13,7)-6.0*ov(15,7) &
     & -2.0*ov(20,5)+3.0*ov(13,5)+3.0*ov(15,5) &
     & -2.0*ov(20,6)+3.0*ov(13,6)+3.0*ov(15,6)
      if(l2.eq.9) goto 15009
      g(14,10)=2.0*ov(20,11)-3.0*ov(13,11)-3.0*ov(15,11)
      g(14,11)=2.0*ov(20,13)-3.0*ov(13,13)-3.0*ov(15,13) &
     & -2.0*ov(20,15)+3.0*ov(13,15)+3.0*ov(15,15)
      g(14,12)=2.0*ov(20,18)-3.0*ov(13,18)-3.0*ov(15,18) &
     & -6.0*ov(20,14)+9.0*ov(13,14)+9.0*ov(15,14)
      g(14,13)=6.0*ov(20,12)-9.0*ov(13,12)-9.0*ov(15,12) &
     & -2.0*ov(20,19)+3.0*ov(13,19)+3.0*ov(15,19)
       g(14,14)=4.0*ov(20,20)-6.0*ov(13,20)-6.0*ov(15,20) &
     & -6.0*ov(20,13)+9.0*ov(13,13)+9.0*ov(15,13) &
     & -6.0*ov(20,15)+9.0*ov(13,15)+9.0*ov(15,15)
      g(14,15)=8.0*ov(20,16)-12.0*ov(13,16)-12.0*ov(15,16) &
     & -2.0*ov(20,18)+3.0*ov(13,18)+3.0*ov(15,18) &
     & -2.0*ov(20,14)+3.0*ov(13,14)+3.0*ov(15,14)
      g(14,16)=8.0*ov(20,17)-12.0*ov(13,17)-12.0*ov(15,17) &
     & -2.0*ov(20,12)+3.0*ov(13,12)+3.0*ov(15,12) &
     & -2.0*ov(20,19)+3.0*ov(13,19)+3.0*ov(15,19)
! c
 15009 continue
      g(15,1)=4.0*ov(16,1)-ov(18,1)-ov(14,1)
      if(l2.eq.1) goto 16009
      g(15,2)=4.0*ov(16,2)-ov(18,2)-ov(14,2)
      g(15,3)=4.0*ov(16,3)-ov(18,3)-ov(14,3)
      g(15,4)=4.0*ov(16,4)-ov(18,4)-ov(14,4)
      if(l2.eq.4) goto 16009
      g(15,5)=4.0*ov(16,8)-ov(18,8)-ov(14,8)
      g(15,6)=4.0*ov(16,9)-ov(18,9)-ov(14,9)
      g(15,7)=4.0*ov(16,10)-ov(18,10)-ov(14,10)
      g(15,8)=4.0*ov(16,5)-ov(18,5)-ov(14,5) &
     & -4.0*ov(16,6)+ov(18,6)+ov(14,6)
      g(15,9)=8.0*ov(16,7)-2.0*ov(18,7)-2.0*ov(14,7) &
     & -4.0*ov(16,5)+ov(18,5)+ov(14,5) &
     & -4.0*ov(16,6)+ov(18,6)+ov(14,6)
      if(l2.eq.9) goto 16009
      g(15,10)=4.0*ov(16,11)-ov(18,11)-ov(14,11)
      g(15,11)=4.0*ov(16,13)-ov(18,13)-ov(14,13) &
     & -4.0*ov(16,15)+ov(18,15)+ov(14,15)
      g(15,12)=4.0*ov(16,18)-ov(18,18)-ov(14,18) &
     & -12.0*ov(16,14)+3.0*ov(18,14)+3.0*ov(14,14)
      g(15,13)=12.0*ov(16,12)-3.0*ov(18,12)-3.0*ov(14,12) &
     & -4.0*ov(16,19)+ov(18,19)+ov(14,19)
      g(15,14)=8.0*ov(16,20)-2.0*ov(18,20)-2.0*ov(14,20) &
     & -12.0*ov(16,13)+3.0*ov(18,13)+3.0*ov(14,13) &
     & -12.0*ov(16,15)+3.0*ov(18,15)+3.0*ov(14,15)
      g(15,15)=16.0*ov(16,16)-4.0*ov(18,16)-4.0*ov(14,16) &
     & -4.0*ov(16,18)+ov(18,18)+ov(14,18) &
     & -4.0*ov(16,14)+ov(18,14)+ov(14,14)
      g(15,16)=16.0*ov(16,17)-4.0*ov(18,17)-4.0*ov(14,17) &
     & -4.0*ov(16,12)+ov(18,12)+ov(14,12) &
     & -4.0*ov(16,19)+ov(18,19)+ov(14,19)
! c
 16009 continue
      g(16,1)=4.0*ov(17,1)-ov(12,1)-ov(19,1)
      if(l2.eq.1) goto 17009
      g(16,2)=4.0*ov(17,2)-ov(12,2)-ov(19,2)
      g(16,3)=4.0*ov(17,3)-ov(12,3)-ov(19,3)
      g(16,4)=4.0*ov(17,4)-ov(12,4)-ov(19,4)
      if(l2.eq.4) goto 17009
      g(16,5)=4.0*ov(17,8)-ov(12,8)-ov(19,8)
      g(16,6)=4.0*ov(17,9)-ov(12,9)-ov(19,9)
      g(16,7)=4.0*ov(17,10)-ov(12,10)-ov(19,10)
      g(16,8)=4.0*ov(17,5)-ov(12,5)-ov(19,5) &
     & -4.0*ov(17,6)+ov(12,6)+ov(19,6)
      g(16,9)=8.0*ov(17,7)-2.0*ov(12,7)-2.0*ov(19,7) &
     & -4.0*ov(17,5)+ov(12,5)+ov(19,5) &
     & -4.0*ov(17,6)+ov(12,6)+ov(19,6)
      if(l2.eq.9) goto 17009
      g(16,10)=4.0*ov(17,11)-ov(12,11)-ov(19,11)
      g(16,11)=4.0*ov(17,13)-ov(12,13)-ov(19,13) &
     & -4.0*ov(17,15)+ov(12,15)+ov(19,15)
      g(16,12)=4.0*ov(17,18)-ov(12,18)-ov(19,18) &
     & -12.0*ov(17,14)+3.0*ov(12,14)+3.0*ov(19,14)
      g(16,13)=12.0*ov(17,12)-3.0*ov(12,12)-3.0*ov(19,12) &
     & -4.0*ov(17,19)+ov(12,19)+ov(19,19)
      g(16,14)=8.0*ov(17,20)-2.0*ov(12,20)-2.0*ov(19,20) &
     & -12.0*ov(17,13)+3.0*ov(12,13)+3.0*ov(19,13) &
     & -12.0*ov(17,15)+3.0*ov(12,15)+3.0*ov(19,15)
      g(16,15)=16.0*ov(17,16)-4.0*ov(12,16)-4.0*ov(19,16) &
     & -4.0*ov(17,18)+ov(12,18)+ov(19,18) &
     & -4.0*ov(17,14)+ov(12,14)+ov(19,14)
      g(16,16)=16.0*ov(17,17)-4.0*ov(12,17)-4.0*ov(19,17) &
     & -4.0*ov(17,12)+ov(12,12)+ov(19,12) &
     & -4.0*ov(17,19)+ov(12,19)+ov(19,19)
 17009 continue
      return                                                            


      end subroutine threeCentInteg

      subroutine KEInteg(g,l1,l2,al1,al2,r,t)

      use O_Kinds
      use O_Constants

      implicit none
! c
! c     g(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
! c     9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
! c     14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy
! c
! c     wo(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
! c     10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
! c     18,xxx; 19,yyy; 20,zzz
! c 
      ! define the dummy variables passed to this subroutine.
      real (kind=double), dimension (16,16) :: g
      integer :: l1, l2
      real (kind=double) :: al1
      real (kind=double) :: al2
      real (kind=double), dimension (3) :: r
      real (kind=double), dimension (3) :: t

      ! define local variables
      integer :: nff
      real (kind=double), dimension (20,20) :: wk
      integer :: ndd,i,j,k,l,m1,m2,m3
      real (kind=double) :: a, as, ac, ap, at, a1s, a2s, a1c, a12, a2c
      real (kind=double) :: b, aa1, aa2, asa1, xs, axs, c, abc, asbc
      real (kind=double) :: acbc, xi2, xj2, wkw
      real (kind=double), dimension (3) :: x
! c

 100  continue
! c
      if(l1.eq.16.or.l2.eq.16) then
         nff=0
      else
         nff=1
      endif
      g(:,:)=0.0_double
      wk(:,:)=0.0_double

      ndd=0
      x(1)=t(1)-r(1)
      x(2)=t(2)-r(2)
      x(3)=t(3)-r(3)
      a=al1*al2/(al1+al2)
      as=a*a
      ac=as*a
      ap=ac*a
      at=ap*a
      a1s=al1**2
      a2s=al2**2
      a1c=al1**3
      a2c=al2**3
      a12=al1*al2
      b=pi/(al1+al2)
      b=dsqrt(b)
      b=b*b*b
      aa1=a*al1
      aa2=a*al2
      asa1=as*al1
      wk(:,:) = 0.0_double
      xs=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
      axs=a*xs
      if( xs.lt.1.0d-14 ) go to 150
      c=a*xs
      c=dexp(-c)
      abc=a*b*c
      asbc=as*b*c
      acbc=ac*b*c
! c this ss
      wk(1,1)=a*b*c*(3.d0-2.d0*a*xs)
      do 105 i=1,3
      xi2=x(i)**2
! c this is spx
      wk(1,i+1)=-asbc*(5.d0-2.d0*axs)*x(i)/al2
! c this is pxs
      wk(i+1,1)= asbc*(5.d0-2.d0*axs)*x(i)/al1
      if(ndd.ne.0) go to 1000
! c this is sdxx
      wk(1,i+4)=.5d0*abc*(3.d0*al2-5.d0*a+2.d0*a*(a-al2)*xs+14.d0*as* &
     &xi2-4.d0*ac*xs*xi2)/(a2s)
! c this is dxxs
      wk(i+4,1)=.5d0*abc*(3.d0*al1-5.d0*a+2.d0*a*(a-al1)*xs+14.d0*as* &
     &x(i)**2-4.d0*ac*xs*xi2)/(a1s)
 1000 continue
! c this is pxpx
      wk(i+1,i+1)=asbc*(5.d0-2.d0*axs-14.d0*a*xi2+4.d0*as*xs*x( &
     &i)**2)/(2.d0*a12)
      if(ndd.ne.0) go to 1001
! c this is dxxpx
      wk(i+4,i+1)=-asbc*(5.d0*al1-21.d0*a+18.d0*as*xi2+2.d0*a*(3.d0 &
     &*a-al1)*xs-4.d0*ac*xs*xi2)*x(i)/(2.d0*a1s*al2)
! c this is pxdxx
      wk(i+1,i+4)= asbc*(5.d0*al2-21.d0*a+18.d0*as*xi2+2.d0*a*(3.d0 &
     &*a-al2)*xs-4.d0*ac*xs*xi2)*x(i)/(2.d0*al1*a2s)
! c this is dxxdxx
      wkw=2.d0*a12-21.d0*as-(14.d0*a*a12-108.d0*ac)*xi2+6.d0*ac* &
     &xs+(4.d0*as*a12-24.d0*a**4)*xs*xi2-44.d0*a**4*x(i)**4+8.d0 &
     &*a**5*xs*x(i)**4
      wk(i+4,i+4)=-abc*wkw/(4.d0*a1s*a2s)
 1001 continue
  105 continue
      if(ndd.ne.0) go to 1002
      k=7
      do 110 i=1,2
      do 110 j=2,3
      if (i.eq.j) go to 110
      k=k+1
! c this is sdxy
      wk(1,k)=acbc*x(i)*x(j)*(7.d0-2.d0*axs)/(a2s)
! c this is dxys
      wk(k,1)=acbc*x(i)*x(j)*(7.d0-2.d0*axs)/(a1s)
  110 continue
 1002 continue
      do 115 i=1,3
      xi2=x(i)**2
      do 115 j=1,3
      xj2=x(j)**2
      if (i.eq.j) go to 115
! c this is pxpy
      wk(i+1,j+1)=-acbc*x(i)*x(j)*(7.d0-2.d0*axs)/(a12)
      if(ndd.ne.0) go to 1003
! c this is pxdyy
      wk(i+1,j+4)= asbc*x(i)*(5.d0*al2-7.d0*a-2.d0*a*(al2-a)*xs+18.d0* &
     &as*xj2-4.d0*ac*xs*xj2)/(2.d0*al1*a2s)
! c this is dxxpy
      wk(i+4,j+1)=-asbc*x(j)*(5.d0*al1-7.d0*a-2.d0*a*(al1-a)*xs+18.d0* &
     &as*xi2-4.d0*ac*xs*xi2)/(2.d0*a1s*al2)
! c this is dxxdyy
      wkw=7.d0*as-2.d0*a12+2.d0*as*(7.d0*al2-9.d0*a)*xi2+2.d0*as*( &
     &7.d0*al1-9.d0*a)*xj2-2.d0*ac*xs+4.d0*ac*(a-al2)*xs*xi2+4.d0* &
     &ac*(a-al1)*xs*xj2+44.d0*a**4*xi2*xj2-8.d0*a**5*xs &
     & *x(i)**2*xj2
      wk(i+4,j+4)=abc*wkw/(4.d0*a1s*a2s)
 1003 continue
  115 continue
      if(ndd.ne.0) go to 1004
      do 130 i=1,3
      xi2=x(i)**2
      do 130 j=1,2
      xj2=x(j)**2
      do 130 k=2,3
      if (j.eq.k) go to 130
      if (i.ne.j.and.i.ne.k) go to 120
      if (j.eq.i) go to 118
      l=j
      go to 119
  118 l=k
  119 continue
! c this is pxdxy
      wk(i+1,j+k+5)=-acbc*(7.d0-18.d0*a*xi2-2.d0*axs+4.d0*as*xs &
     &*xi2)*x(l)/(2.d0*al1*a2s)
! c this is dxypx
      wk(j+k+5,i+1)= acbc*(7.d0-18.d0*a*xi2-2.d0*axs+4.d0*as*xs &
     &*xi2)*x(l)/(2.d0*a1s*al2)
! c this is dxxdxy
      wk(i+4,j+k+5)=acbc*(7.d0*al1-27.d0*a+22.d0*as*xi2-2.d0*a*(al1 &
     &-3.d0*a)*xs-4.d0*ac*xs*xi2)*x(i)*x(l)/(2.d0*a1s*a2s)
! c this is dxydxx
      wk(j+k+5,i+4)=acbc*(7.d0*al2-27.d0*a+22.d0*as*xi2-2.d0*a*(al2 &
     &-3.d0*a)*xs-4.d0*ac*xs*xi2)*x(i)*x(l)/(2.d0*a1s*a2s)
      go to 130
  120 continue
! c this is pxdyz
      wk(i+1,j+k+5)= a**4*b*c*x(i)*x(j)*x(k)*(9.d0-2.d0*axs)/(al1 &
     &   *a2s)
! c this is dyzpx
      wk(j+k+5,i+1)=-a**4*b*c*x(i)*x(j)*x(k)*(9.d0-2.d0*axs)/(a1s*al2)
! c this is dxxdyz
      wk(i+4,j+k+5)=acbc*x(j)*x(k)*(7.d0*al1-9.d0*a+22.d0*as*xi2+ &
     &2.d0*a*(a-al1)*xs-4.d0*ac*xs*xi2)/(2.d0*a1s*a2s)
! c this is dyzdxx
      wk(j+k+5,i+4)=acbc*x(j)*x(k)*(7.d0*al2-9.d0*a+22.d0*as*xi2+ &
     &2.d0*a*(a-al2)*xs-4.d0*ac*xs*xi2)/(2.d0*a1s*a2s)
  130 continue
      do 146 i=1,2
      xi2=x(i)**2
      do 146 j=2,3
      xj2=x(j)**2
      if (i.eq.j) go to 146
      do 145 k=1,2
      do 145 l=2,3
      if (k.eq.l) go to 145
      if (i.eq.k.and.j.eq.l) go to 140
      if (i.eq.k) go to 131
      if (i.eq.l) go to 132
      if (j.eq.k) go to 133
      if  (j.eq.l) go to 134
  131 m1=i
      m2=j
      m3=l
      go to 135
  132 m1=i
      m2=j
      m3=k
      go to 135
  133 m1=j
      m2=i
      m3=l
      go to 135
  134 m1=j
      m2=i
      m3=k
      go to 135
! c this is dxydyz
  135 continue
      wk(i+j+5,k+l+5)=-a**4*b*c*(9.d0-2.d0*axs-22.d0*a*x(m1)**2+4.d0* &
     &as*xs*x(m1)**2)*x(m2)*x(m3)/(2.d0*a1s*a2s)
      go to 145
! c this is dxydxy
  140 continue
      wk(i+j+5,k+l+5)= acbc*(7.d0-2.d0*axs+(4.d0*as*xs-18.d0*a)*(x( &
     &i)**2+x(j)**2)+(44.d0*as-8.d0*ac*xs)*x(i)**2*x(j)**2)/(4.d0*a1s* &
     &a2s)
  145 continue
  146 continue
 1004 continue
      go to 200
  150 continue
! c this ss
      wk(1,1)=3.d0*a*b
      do 155 i=1,3
      xi2=x(i)**2
      if(ndd.ne.0) go to 1005
! c this is sdxx
      wk(1,i+4)=.5d0*a*b*(3.d0-5.d0*a/al2)/al2
! c this is dxxs
      wk(i+4,1)=.5d0*a*b*(3.d0-5.d0*a/al1)/al1
 1005 continue
! c this is pxpx
      wk(i+1,i+1)= 5.d0*b*as/(2.d0*a12)
      if(ndd.ne.0) go to 1006
! c this is dxxdxx
      wk(i+4,i+4)=-b*a*(2.d0-21.d0*as/(a12))/(4.d0*a12)
 1006 continue
  155 continue
      if(ndd.ne.0) go to 1007
      do 165 i=1,3
      do 165 j=1,3
      if (i.eq.j) go to 165
! c this is dxxdyy
      wk(i+4,j+4)=-b*a*(2.d0-7.d0*as/(a12))/(4.d0*a12)
  165 continue
      do 175 i=8,10
! c this is dxydxy
  175 wk(i,i)=7.d0*a*as*b/(4.d0*a1s*a2s)
 1007 continue
  200 continue
      if(nff.eq.0) then
         call kinef2(al1,al2,r,t,wk)
         call kinef1(al1,al2,r,t,wk)
      endif

      continue
      g(1,1)=wk(1,1)
      if(l2.eq.1) goto 2009
      g(1,2)=wk(1,2)
      g(1,3)=wk(1,3)
      g(1,4)=wk(1,4)
      if(l2.eq.4) goto 2009
      g(1,5)=wk(1,8)
      g(1,6)=wk(1,9)
      g(1,7)=wk(1,10)
      g(1,8)=wk(1,5)-wk(1,6)
      g(1,9)=2.0*wk(1,7)-wk(1,5)-wk(1,6)
      if(l2.eq.9) goto 2009
      g(1,10)=wk(1,11)
      g(1,11)=wk(1,13)-wk(1,15)
      g(1,12)=wk(1,18)-3.0*wk(1,14)
      g(1,13)=3.0*wk(1,12)-wk(1,19)
      g(1,14)=2.0*wk(1,20)-3.0*wk(1,13)-3.0*wk(1,15)
      g(1,15)=4.0*wk(1,16)-wk(1,18)-wk(1,14)
      g(1,16)=4.0*wk(1,17)-wk(1,12)-wk(1,19)
! c
 2009  continue
      if(l1.eq.1) return
      g(2,1)=wk(2,1)
      if(l2.eq.1) goto 3009
      g(2,2)=wk(2,2)
      g(2,3)=wk(2,3)
      g(2,4)=wk(2,4)
      if(l2.eq.4) goto 3009
      g(2,5)=wk(2,8)
      g(2,6)=wk(2,9)
      g(2,7)=wk(2,10)
      g(2,8)=wk(2,5)-wk(2,6)
      g(2,9)=2.0*wk(2,7)-wk(2,5)-wk(2,6)
      if(l2.eq.9) goto 3009
      g(2,10)=wk(2,11)
      g(2,11)=wk(2,13)-wk(2,15)
      g(2,12)=wk(2,18)-3.0*wk(2,14)
      g(2,13)=3.0*wk(2,12)-wk(2,19)
      g(2,14)=2.0*wk(2,20)-3.0*wk(2,13)-3.0*wk(2,15)
      g(2,15)=4.0*wk(2,16)-wk(2,18)-wk(2,14)
      g(2,16)=4.0*wk(2,17)-wk(2,12)-wk(2,19)
! c
 3009  continue
      g(3,1)=wk(3,1)
      if(l2.eq.1) goto 4009
      g(3,2)=wk(3,2)
      g(3,3)=wk(3,3)
      g(3,4)=wk(3,4)
      if(l2.eq.4) goto 4009
      g(3,5)=wk(3,8)
      g(3,6)=wk(3,9)
      g(3,7)=wk(3,10)
      g(3,8)=wk(3,5)-wk(3,6)
      g(3,9)=2.0*wk(3,7)-wk(3,5)-wk(3,6)
      if(l2.eq.9) goto 4009
      g(3,10)=wk(3,11)
      g(3,11)=wk(3,13)-wk(3,15)
      g(3,12)=wk(3,18)-3.0*wk(3,14)
      g(3,13)=3.0*wk(3,12)-wk(3,19)
      g(3,14)=2.0*wk(3,20)-3.0*wk(3,13)-3.0*wk(3,15)
      g(3,15)=4.0*wk(3,16)-wk(3,18)-wk(3,14)
      g(3,16)=4.0*wk(3,17)-wk(3,12)-wk(3,19)
! c
 4009  continue
      g(4,1)=wk(4,1)
      if(l2.eq.1) goto 5009
      g(4,2)=wk(4,2)
      g(4,3)=wk(4,3)
      g(4,4)=wk(4,4)
      if(l2.eq.4) goto 5009
      g(4,5)=wk(4,8)
      g(4,6)=wk(4,9)
      g(4,7)=wk(4,10)
      g(4,8)=wk(4,5)-wk(4,6)
      g(4,9)=2.0*wk(4,7)-wk(4,5)-wk(4,6)
      if(l2.eq.9) goto 5009
      g(4,10)=wk(4,11)
      g(4,11)=wk(4,13)-wk(4,15)
      g(4,12)=wk(4,18)-3.0*wk(4,14)
      g(4,13)=3.0*wk(4,12)-wk(4,19)
      g(4,14)=2.0*wk(4,20)-3.0*wk(4,13)-3.0*wk(4,15)
      g(4,15)=4.0*wk(4,16)-wk(4,18)-wk(4,14)
      g(4,16)=4.0*wk(4,17)-wk(4,12)-wk(4,19)
! c
 5009  continue
      if(l1.eq.4) return
      g(5,1)=wk(8,1)
      if(l2.eq.1) goto 6009
      g(5,2)=wk(8,2)
      g(5,3)=wk(8,3)
      g(5,4)=wk(8,4)
      if(l2.eq.4) goto 6009
      g(5,5)=wk(8,8)
      g(5,6)=wk(8,9)
      g(5,7)=wk(8,10)
      g(5,8)=wk(8,5)-wk(8,6)
      g(5,9)=2.0*wk(8,7)-wk(8,5)-wk(8,6)
      if(l2.eq.9) goto 6009
      g(5,10)=wk(8,11)
      g(5,11)=wk(8,13)-wk(8,15)
      g(5,12)=wk(8,18)-3.0*wk(8,14)
      g(5,13)=3.0*wk(8,12)-wk(8,19)
      g(5,14)=2.0*wk(8,20)-3.0*wk(8,13)-3.0*wk(8,15)
      g(5,15)=4.0*wk(8,16)-wk(8,18)-wk(8,14)
      g(5,16)=4.0*wk(8,17)-wk(8,12)-wk(8,19)
! c
 6009  continue
      g(6,1)=wk(9,1)
      if(l2.eq.1) goto 7009
      g(6,2)=wk(9,2)
      g(6,3)=wk(9,3)
      g(6,4)=wk(9,4)
      if(l2.eq.4) goto 7009
      g(6,5)=wk(9,8)
      g(6,6)=wk(9,9)
      g(6,7)=wk(9,10)
      g(6,8)=wk(9,5)-wk(9,6)
      g(6,9)=2.0*wk(9,7)-wk(9,5)-wk(9,6)
      if(l2.eq.9) goto 7009
      g(6,10)=wk(9,11)
      g(6,11)=wk(9,13)-wk(9,15)
      g(6,12)=wk(9,18)-3.0*wk(9,14)
      g(6,13)=3.0*wk(9,12)-wk(9,19)
      g(6,14)=2.0*wk(9,20)-3.0*wk(9,13)-3.0*wk(9,15)
      g(6,15)=4.0*wk(9,16)-wk(9,18)-wk(9,14)
      g(6,16)=4.0*wk(9,17)-wk(9,12)-wk(9,19)
! c
 7009  continue
      g(7,1)=wk(10,1)
      if(l2.eq.1) goto 8009
      g(7,2)=wk(10,2)
      g(7,3)=wk(10,3)
      g(7,4)=wk(10,4)
      if(l2.eq.4) goto 8009
      g(7,5)=wk(10,8)
      g(7,6)=wk(10,9)
      g(7,7)=wk(10,10)
      g(7,8)=wk(10,5)-wk(10,6)
      g(7,9)=2.0*wk(10,7)-wk(10,5)-wk(10,6)
      if(l2.eq.9) goto 8009
      g(7,10)=wk(10,11)
      g(7,11)=wk(10,13)-wk(10,15)
      g(7,12)=wk(10,18)-3.0*wk(10,14)
      g(7,13)=3.0*wk(10,12)-wk(10,19)
      g(7,14)=2.0*wk(10,20)-3.0*wk(10,13)-3.0*wk(10,15)
      g(7,15)=4.0*wk(10,16)-wk(10,18)-wk(10,14)
      g(7,16)=4.0*wk(10,17)-wk(10,12)-wk(10,19)
! c
 8009  continue
      g(8,1)=wk(5,1)-wk(6,1)
      if(l2.eq.1) goto 9009
      g(8,2)=wk(5,2)-wk(6,2)
      g(8,3)=wk(5,3)-wk(6,3)
      g(8,4)=wk(5,4)-wk(6,4)
      if(l2.eq.4) goto 9009
      g(8,5)=wk(5,8)-wk(6,8)
      g(8,6)=wk(5,9)-wk(6,9)
      g(8,7)=wk(5,10)-wk(6,10)
      g(8,8)=wk(5,5)-wk(6,5)-wk(5,6)+wk(6,6)
      g(8,9)=2.0*wk(5,7)-2.0*wk(6,7)-wk(5,5) &
     & +wk(6,5)-wk(5,6)+wk(6,6)
      if(l2.eq.9) goto 9009
      g(8,10)=wk(5,11)-wk(6,11)
      g(8,11)=wk(5,13)-wk(6,13)-wk(5,15)+wk(6,15)
      g(8,12)=wk(5,18)-wk(6,18)-3.0*wk(5,14) &
     & +3.0*wk(6,14)
      g(8,13)=3.0*wk(5,12)-3.0*wk(6,12)-wk(5,19)+wk(6,19)
      g(8,14)=2.0*wk(5,20)-2.0*wk(6,20)-3.0*wk(5,13) &
     & +3.0*wk(6,13)-3.0*wk(5,15)+3.0*wk(6,15)
      g(8,15)=4.0*wk(5,16)-4.0*wk(6,16)-wk(5,18)+wk(6,18) &
     & -wk(5,14)+wk(6,14)
      g(8,16)=4.0*wk(5,17)-4.0*wk(6,17)-wk(5,12)+wk(6,12) &
     & -wk(5,19)+wk(6,19)
! c
 9009  continue
      g(9,1)=2.0*wk(7,1)-wk(5,1)-wk(6,1)
      if(l2.eq.1) goto 10009
      g(9,2)=2.0*wk(7,2)-wk(5,2)-wk(6,2)
      g(9,3)=2.0*wk(7,3)-wk(5,3)-wk(6,3)
      g(9,4)=2.0*wk(7,4)-wk(5,4)-wk(6,4)
      if(l2.eq.4) goto 10009
      g(9,5)=2.0*wk(7,8)-wk(5,8)-wk(6,8)
      g(9,6)=2.0*wk(7,9)-wk(5,9)-wk(6,9)
      g(9,7)=2.0*wk(7,10)-wk(5,10)-wk(6,10)
      g(9,8)=2.0*wk(7,5)-wk(5,5)-wk(6,5)-2.0*wk(7,6)+wk(5,6)+wk(6,6)
      g(9,9)=4.0*wk(7,7)-2.0*wk(5,7)-2.0*wk(6,7)-2.0*wk(7,5)+wk(5,5) &
     & +wk(6,5)-2.0*wk(7,6)+wk(5,6)+wk(6,6)
      if(l2.eq.9) goto 10009
      g(9,10)=2.0*wk(7,11)-wk(5,11)-wk(6,11)
      g(9,11)=2.0*wk(7,13)-wk(5,13)-wk(6,13)-2.0*wk(7,15) &
     & +wk(5,15)+wk(6,15)
      g(9,12)=2.0*wk(7,18)-wk(5,18)-wk(6,18)-6.0*wk(7,14)+ &
     & 3.0*wk(5,14)+3.0*wk(6,14)
      g(9,13)=6.0*wk(7,12)-3.0*wk(5,12)-3.0*wk(6,12) &
     & -2.0*wk(7,19)+wk(5,19)+wk(6,19)
      g(9,14)=4.0*wk(7,20)-2.0*wk(5,20)-2.0*wk(6,20)-6.0*wk(7,13)+ &
     & 3.0*wk(5,13)+3.0*wk(6,13)-6.0*wk(7,15)+3.0*wk(5,15)+3.0*wk(6,15)
      g(9,15)=8.0*wk(7,16)-4.0*wk(5,16)-4.0*wk(6,16)-2.0*wk(7,18) &
     & +wk(5,18)+wk(6,18)-2.0*wk(7,14)+wk(5,14)+wk(6,14)
      g(9,16)=8.0*wk(7,17)-4.0*wk(5,17)-4.0*wk(6,17)-2.0*wk(7,12) &
     & +wk(5,12)+wk(6,12)-2.0*wk(7,19)+wk(5,19)+wk(6,19)
! c
 10009 continue
      if(l1.eq.9) return
      g(10,1)=wk(11,1)
      if(l2.eq.1) goto 11009
      g(10,2)=wk(11,2)
      g(10,3)=wk(11,3)
      g(10,4)=wk(11,4)
      if(l2.eq.4) goto 11009
      g(10,5)=wk(11,8)
      g(10,6)=wk(11,9)
      g(10,7)=wk(11,10)
      g(10,8)=wk(11,5)-wk(11,6)
      g(10,9)=2.0*wk(11,7)-wk(11,5)-wk(11,6)
      if(l2.eq.9) goto 11009
      g(10,10)=wk(11,11)
      g(10,11)=wk(11,13)-wk(11,15)
      g(10,12)=wk(11,18)-3.0*wk(11,14)
      g(10,13)=3.0*wk(11,12)-wk(11,19)
      g(10,14)=2.0*wk(11,20)-3.0*wk(11,13)-3.0*wk(11,15)
      g(10,15)=4.0*wk(11,16)-wk(11,18)-wk(11,14)
      g(10,16)=4.0*wk(11,17)-wk(11,12)-wk(11,19)
! c
 11009 continue
      g(11,1)=wk(13,1)-wk(15,1)
      if(l2.eq.1) goto 12009
      g(11,2)=wk(13,2)-wk(15,2)
      g(11,3)=wk(13,3)-wk(15,3)
      g(11,4)=wk(13,4)-wk(15,4)
      if(l2.eq.4) goto 12009
      g(11,5)=wk(13,8)-wk(15,8)
      g(11,6)=wk(13,9)-wk(15,9)
      g(11,7)=wk(13,10)-wk(15,10)
      g(11,8)=wk(13,5)-wk(15,5)-wk(13,6)+wk(15,6)
      g(11,9)=2.0*wk(13,7)-2.0*wk(15,7)-wk(13,5) &
     & +wk(15,5)-wk(13,6)+wk(15,6)
      if(l2.eq.9) goto 12009
      g(11,10)=wk(13,11)-wk(15,11)
      g(11,11)=wk(13,13)-wk(15,13)-wk(13,15)+wk(15,15)
      g(11,12)=wk(13,18)-wk(15,18)-3.0*wk(13,14) &
     & +3.0*wk(15,14)
      g(11,13)=3.0*wk(13,12)-3.0*wk(15,12)- &
     & wk(13,19)+wk(15,19)
      g(11,14)=2.0*wk(13,20)-2.0*wk(15,20)-3.0*wk(13,13) &
     & +3.0*wk(15,13)-3.0*wk(13,15)+3.0*wk(15,15)
      g(11,15)=4.0*wk(13,16)-4.0*wk(15,16)-wk(13,18)+wk(15,18) &
     & -wk(13,14)+wk(15,14)
      g(11,16)=4.0*wk(13,17)-4.0*wk(15,17)-wk(13,12)+wk(15,12) &
     & -wk(13,19)+wk(15,19)
! c
 12009 continue
      g(12,1)=wk(18,1)-3.0*wk(14,1)
      if(l2.eq.1) goto 13009
      g(12,2)=wk(18,2)-3.0*wk(14,2)
      g(12,3)=wk(18,3)-3.0*wk(14,3)
      g(12,4)=wk(18,4)-3.0*wk(14,4)
      if(l2.eq.4) goto 13009
      g(12,5)=wk(18,8)-3.0*wk(14,8)
      g(12,6)=wk(18,9)-3.0*wk(14,9)
      g(12,7)=wk(18,10)-3.0*wk(14,10)
      g(12,8)=wk(18,5)-3.0*wk(14,5) &
     & -wk(18,6)+3.0*wk(14,6)
      g(12,9)=2.0*wk(18,7)-6.0*wk(14,7) &
     & -wk(18,5)+3.0*wk(14,5) &
     & -wk(18,6)+3.0*wk(14,6)
      if(l2.eq.9) goto 13009
      g(12,10)=wk(18,11)-3.0*wk(14,11)
      g(12,11)=wk(18,13)-3.0*wk(14,13) &
     & -wk(18,15)+3.0*wk(14,15)
      g(12,12)=wk(18,18)-3.0*wk(14,18) &
     & -3.0*wk(18,14)+9.0*wk(14,14)
      g(12,13)=3.0*wk(18,12)-9.0*wk(14,12) &
     & -wk(18,19)+3.0*wk(14,19)
      g(12,14)=2.0*wk(18,20)-6.0*wk(14,20)-3.0*wk(18,13) &
     & +9.0*wk(14,13)-3.0*wk(18,15)+9.0*wk(14,15)
      g(12,15)=4.0*wk(18,16)-12.0*wk(14,16)-wk(18,18) &
     & +3.0*wk(14,18)-wk(18,14)+3.0*wk(14,14)
      g(12,16)=4.0*wk(18,17)-12.0*wk(14,17)-wk(18,12) &
     & +3.0*wk(14,12)-wk(18,19)+3.0*wk(14,19)
! c
 13009 continue
      g(13,1)=3.0*wk(12,1)-wk(19,1)
      if(l2.eq.1) goto 14009
      g(13,2)=3.0*wk(12,2)-wk(19,2)
      g(13,3)=3.0*wk(12,3)-wk(19,3)
      g(13,4)=3.0*wk(12,4)-wk(19,4)
      if(l2.eq.4) goto 14009
      g(13,5)=3.0*wk(12,8)-wk(19,8)
      g(13,6)=3.0*wk(12,9)-wk(19,9)
      g(13,7)=3.0*wk(12,10)-wk(19,10)
      g(13,8)=3.0*wk(12,5)-wk(19,5)-3.0*wk(12,6)+wk(19,6)
      g(13,9)=6.0*wk(12,7)-2.0*wk(19,7)-3.0*wk(12,5) &
     & +wk(19,5)-3.0*wk(12,6)+wk(19,6)
      if(l2.eq.9) goto 14009
      g(13,10)=3.0*wk(12,11)-wk(19,11)
      g(13,11)=3.0*wk(12,13)-wk(19,13)-3.0*wk(12,15)+wk(19,15)
      g(13,12)=3.0*wk(12,18)-wk(19,18)-9.0*wk(12,14)+3.0*wk(19,14)
      g(13,13)=9.0*wk(12,12)-3.0*wk(19,12)-3.0*wk(12,19)+wk(19,19)
      g(13,14)=6.0*wk(12,20)-2.0*wk(19,20)-9.0*wk(12,13) &
     & +3.0*wk(19,13)-9.0*wk(12,15)+3.0*wk(19,15)
      g(13,15)=12.0*wk(12,16)-4.0*wk(19,16)-3.0*wk(12,18) &
     & +wk(19,18)-3.0*wk(12,14)+wk(19,14)
      g(13,16)=12.0*wk(12,17)-4.0*wk(19,17)-3.0*wk(12,12) &
     & +wk(19,12)-3.0*wk(12,19)+wk(19,19)
! c
 14009 continue
      g(14,1)=2.0*wk(20,1)-3.0*wk(13,1)-3.0*wk(15,1)
      if(l2.eq.1) goto 15009
      g(14,2)=2.0*wk(20,2)-3.0*wk(13,2)-3.0*wk(15,2)
      g(14,3)=2.0*wk(20,3)-3.0*wk(13,3)-3.0*wk(15,3)
      g(14,4)=2.0*wk(20,4)-3.0*wk(13,4)-3.0*wk(15,4)
      if(l2.eq.4) goto 15009
      g(14,5)=2.0*wk(20,8)-3.0*wk(13,8)-3.0*wk(15,8)
      g(14,6)=2.0*wk(20,9)-3.0*wk(13,9)-3.0*wk(15,9)
      g(14,7)=2.0*wk(20,10)-3.0*wk(13,10)-3.0*wk(15,10)
      g(14,8)=2.0*wk(20,5)-3.0*wk(13,5)-3.0*wk(15,5) &
     & -2.0*wk(20,6)+3.0*wk(13,6)+3.0*wk(15,6)
      g(14,9)=4.0*wk(20,7)-6.0*wk(13,7)-6.0*wk(15,7) &
     & -2.0*wk(20,5)+3.0*wk(13,5)+3.0*wk(15,5) &
     & -2.0*wk(20,6)+3.0*wk(13,6)+3.0*wk(15,6)
      if(l2.eq.9) goto 15009
      g(14,10)=2.0*wk(20,11)-3.0*wk(13,11)-3.0*wk(15,11)
      g(14,11)=2.0*wk(20,13)-3.0*wk(13,13)-3.0*wk(15,13) &
     & -2.0*wk(20,15)+3.0*wk(13,15)+3.0*wk(15,15)
      g(14,12)=2.0*wk(20,18)-3.0*wk(13,18)-3.0*wk(15,18) &
     & -6.0*wk(20,14)+9.0*wk(13,14)+9.0*wk(15,14)
      g(14,13)=6.0*wk(20,12)-9.0*wk(13,12)-9.0*wk(15,12) &
     & -2.0*wk(20,19)+3.0*wk(13,19)+3.0*wk(15,19)
       g(14,14)=4.0*wk(20,20)-6.0*wk(13,20)-6.0*wk(15,20) &
     & -6.0*wk(20,13)+9.0*wk(13,13)+9.0*wk(15,13) &
     & -6.0*wk(20,15)+9.0*wk(13,15)+9.0*wk(15,15)
      g(14,15)=8.0*wk(20,16)-12.0*wk(13,16)-12.0*wk(15,16) &
     & -2.0*wk(20,18)+3.0*wk(13,18)+3.0*wk(15,18) &
     & -2.0*wk(20,14)+3.0*wk(13,14)+3.0*wk(15,14)
      g(14,16)=8.0*wk(20,17)-12.0*wk(13,17)-12.0*wk(15,17) &
     & -2.0*wk(20,12)+3.0*wk(13,12)+3.0*wk(15,12) &
     & -2.0*wk(20,19)+3.0*wk(13,19)+3.0*wk(15,19)
! c
 15009 continue
      g(15,1)=4.0*wk(16,1)-wk(18,1)-wk(14,1)
      if(l2.eq.1) goto 16009
      g(15,2)=4.0*wk(16,2)-wk(18,2)-wk(14,2)
      g(15,3)=4.0*wk(16,3)-wk(18,3)-wk(14,3)
      g(15,4)=4.0*wk(16,4)-wk(18,4)-wk(14,4)
      if(l2.eq.4) goto 16009
      g(15,5)=4.0*wk(16,8)-wk(18,8)-wk(14,8)
      g(15,6)=4.0*wk(16,9)-wk(18,9)-wk(14,9)
      g(15,7)=4.0*wk(16,10)-wk(18,10)-wk(14,10)
      g(15,8)=4.0*wk(16,5)-wk(18,5)-wk(14,5) &
     & -4.0*wk(16,6)+wk(18,6)+wk(14,6)
      g(15,9)=8.0*wk(16,7)-2.0*wk(18,7)-2.0*wk(14,7) &
     & -4.0*wk(16,5)+wk(18,5)+wk(14,5) &
     & -4.0*wk(16,6)+wk(18,6)+wk(14,6)
      if(l2.eq.9) goto 16009
      g(15,10)=4.0*wk(16,11)-wk(18,11)-wk(14,11)
      g(15,11)=4.0*wk(16,13)-wk(18,13)-wk(14,13) &
     & -4.0*wk(16,15)+wk(18,15)+wk(14,15)
      g(15,12)=4.0*wk(16,18)-wk(18,18)-wk(14,18) &
     & -12.0*wk(16,14)+3.0*wk(18,14)+3.0*wk(14,14)
      g(15,13)=12.0*wk(16,12)-3.0*wk(18,12)-3.0*wk(14,12) &
     & -4.0*wk(16,19)+wk(18,19)+wk(14,19)
      g(15,14)=8.0*wk(16,20)-2.0*wk(18,20)-2.0*wk(14,20) &
     & -12.0*wk(16,13)+3.0*wk(18,13)+3.0*wk(14,13) &
     & -12.0*wk(16,15)+3.0*wk(18,15)+3.0*wk(14,15)
      g(15,15)=16.0*wk(16,16)-4.0*wk(18,16)-4.0*wk(14,16) &
     & -4.0*wk(16,18)+wk(18,18)+wk(14,18) &
     & -4.0*wk(16,14)+wk(18,14)+wk(14,14)
      g(15,16)=16.0*wk(16,17)-4.0*wk(18,17)-4.0*wk(14,17) &
     & -4.0*wk(16,12)+wk(18,12)+wk(14,12) &
     & -4.0*wk(16,19)+wk(18,19)+wk(14,19)
! c
 16009 continue
      g(16,1)=4.0*wk(17,1)-wk(12,1)-wk(19,1)
      if(l2.eq.1) goto 17009
      g(16,2)=4.0*wk(17,2)-wk(12,2)-wk(19,2)
      g(16,3)=4.0*wk(17,3)-wk(12,3)-wk(19,3)
      g(16,4)=4.0*wk(17,4)-wk(12,4)-wk(19,4)
      if(l2.eq.4) goto 17009
      g(16,5)=4.0*wk(17,8)-wk(12,8)-wk(19,8)
      g(16,6)=4.0*wk(17,9)-wk(12,9)-wk(19,9)
      g(16,7)=4.0*wk(17,10)-wk(12,10)-wk(19,10)
      g(16,8)=4.0*wk(17,5)-wk(12,5)-wk(19,5) &
     & -4.0*wk(17,6)+wk(12,6)+wk(19,6)
      g(16,9)=8.0*wk(17,7)-2.0*wk(12,7)-2.0*wk(19,7) &
     & -4.0*wk(17,5)+wk(12,5)+wk(19,5) &
     & -4.0*wk(17,6)+wk(12,6)+wk(19,6)
      if(l2.eq.9) goto 17009
      g(16,10)=4.0*wk(17,11)-wk(12,11)-wk(19,11)
      g(16,11)=4.0*wk(17,13)-wk(12,13)-wk(19,13) &
     & -4.0*wk(17,15)+wk(12,15)+wk(19,15)
      g(16,12)=4.0*wk(17,18)-wk(12,18)-wk(19,18) &
     & -12.0*wk(17,14)+3.0*wk(12,14)+3.0*wk(19,14)
      g(16,13)=12.0*wk(17,12)-3.0*wk(12,12)-3.0*wk(19,12) &
     & -4.0*wk(17,19)+wk(12,19)+wk(19,19)
      g(16,14)=8.0*wk(17,20)-2.0*wk(12,20)-2.0*wk(19,20) &
     & -12.0*wk(17,13)+3.0*wk(12,13)+3.0*wk(19,13) &
     & -12.0*wk(17,15)+3.0*wk(12,15)+3.0*wk(19,15)
      g(16,15)=16.0*wk(17,16)-4.0*wk(12,16)-4.0*wk(19,16) &
     & -4.0*wk(17,18)+wk(12,18)+wk(19,18) &
     & -4.0*wk(17,14)+wk(12,14)+wk(19,14)
      g(16,16)=16.0*wk(17,17)-4.0*wk(12,17)-4.0*wk(19,17) &
     & -4.0*wk(17,12)+wk(12,12)+wk(19,12) &
     & -4.0*wk(17,19)+wk(12,19)+wk(19,19)
 17009 continue
      return                                                            



      end subroutine KEInteg
! c
!*****************************************************
! c
      subroutine kinef1(a1,a2,r,t,wk)

      use O_Kinds
      use O_Constants

      implicit none
      real (kind=double) :: a1,a2
      real (kind=double), dimension (3) :: r,t,x
      real (kind=double), dimension (20,20) :: wk

      integer :: i,j,k,l,m,ms,mm,mt,ma,mb,na,nb
      real (kind=double) :: a, as, ac, ap, at, a1s, a2s, a1c, a12, a2c
      real (kind=double) :: b, aa1, aa2, asa1, xs, axs, c, abc, asbc
      real (kind=double) :: acbc, xi2, xi4, xj2, xj4
      real (kind=double) :: f1, f2, f3, f4, f5, f6, f7, f8, f9, f10
      real (kind=double) :: f11, f12, f13, f14

      x(1)=t(1)-r(1)
      x(2)=t(2)-r(2)
      x(3)=t(3)-r(3)
      a=a1*a2/(a1+a2)
      as=a*a
      ac=as*a
      ap=ac*a
      at=ap*a
      a1s=a1**2
      a2s=a2**2
      a1c=a1**3
      a2c=a2**3
      a12=a1*a2
      b=pi/(a1+a2)
      b=dsqrt(b)
      b=b*b*b
      aa1=a*a1
      aa2=a*a2
      asa1=as*a1
      xs=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
      axs=a*xs
      c=a*xs
      c=dexp(-c)
      abc=a*b*c
      asbc=as*b*c
      acbc=ac*b*c
! c  this loop has exchange in x-y-z
      do 205 i=1,3
      xi2=x(i)**2
      xi4=xi2*xi2
! c   this fxxx s
      f1=asbc*x(i)
      f2=15.*a1-21.*a+(6.*as-6.*aa1)*xs
      f3=18.*as*xi2-4.*ac*xs*xi2
      wk(i+17,1)=f1*(f2+f3)/(2.*a1c)
! c   this s fxxx
      f2=15.*a2-21.*a+(6.*as-6.*aa2)*xs
      wk(1,i+17)=-f1*(f2+f3)/(2.*a2c)
! c   this fxxx px
      f1=asbc
      f2=15.*a1-21.*a
      f3=(6.*as-6.*a*a1)*xs
      f4=(108.*as-42.*a*a1)*xi2
      f5=-(24.*ac-12.*as*a1)*xs*xi2
      f6=-44.*ac*xi4
      f7=8.*ap*xs*xi4
      wk(i+17,i+1)=f1*(f2+f3+f4+f5+f6+f7)/(4.*a1c*a2)
! c   this px fxxx
      f2=15.*a2-21.*a
      f3=(6.*as-6.*aa2)*xs
      f4=(108.*as-42.*a*a2)*xi2
      f5=-(24.*ac-12.*as*a2)*xs*xi2
      wk(i+1,i+17)=f1*(f2+f3+f4+f5+f6+f7)/(4.*a2c*a1)
! c   this fxxx dxx
      f1=asbc*x(i)
      f2=15.*a12-63.*aa1-21.*aa2+135.*as
      f3=-(220.*ac-54.*as*a1-18.*as*a2)*xi2
      f4=-(30.*ac-18.*as*a1-6.*as*a2+6.*a*a12)*xs
      f5=(40.*ap-12.*ac*a1-4.*ac*a2)*xs*xi2
      f6=52.*ap*xi4
      f8=-8.*at*xs*xi4
      wk(i+17,i+4)=f1*(f2+f3+f4+f5+f6+f8)/(4.*a1c*a2s)
! c   this dxx fxxx
      f2=15.*a12-63.*aa2-21.*aa1+135.*as
      f3=-(220.*ac-54.*as*a2-18.*as*a1)*xi2
      f4=-(30.*ac-18.*as*a2-6.*as*a1+6.*a*a12)*xs
      f5=(40.*ap-12.*ac*a2-4.*ac*a1)*xs*xi2
      wk(i+4,i+17)=-f1*(f2+f3+f4+f5+f6+f8)/(4.*a1s*a2c)
! c   this fxxx fxxx
      f1=asbc
      f2=45.*a12-63.*aa1-63.*aa2+135.*as
      f3=-(990.*ac-324.*as*a1-324.*as*a2+126.*a*a12)*xi2
      f4=-(30.*ac-18.*as*a1-18.*as*a2+18.*a*a12)*xs
      f5=(180.*ap-72.*ac*a1-72.*ac*a2+36.*as*a12)*xs*xi2
      f6=(780.*ap-132.*ac*a1-132.*ac*a2)*xi4
      f7=-(120.*at-24.*ap*a1-24.*ap*a2)*xs*xi4
      f8=-120*at*x(i)**6
      f9=16.*a**6*xs*x(i)**6
      wk(i+17,i+17)=f1*(f2+f3+f4+f5+f6+f7+f8+f9)/(8.*a1c*a2c)
  205 continue
! c  this loop has exchanging in x-y-z,x-y,y-z,z-x and six symmetry
      l=11
      do 210 i=1,3
      xi2=x(i)**2
      xi4=xi2*xi2
      do 210 j=1,3
      xj2=x(j)**2
      xj4=xj2*xj2
      if(i.eq.j)go to 210
      l=l+1
! c   this fxxy s
      f1=asbc*x(j)
      f2=5.*a1-7.*a+18.*as*xi2
      f3=axs*(2.*(a-a1)-4.*as*xi2)
      wk(l,1)=f1*(f2+f3)/(2.*a1c)
! c   this s fxxy
      f1=-asbc*x(j)
      f2=5.*a2-7.*a+18.*as*xi2
      f3=axs*(2.*(a-a2)-4.*as*xi2)
      wk(1,l)=f1*(f2+f3)/(2.*a2c)
! c   this fxxy px
      f1=-asbc*x(i)*x(j)
      f2=7.*a*a1-27.*as+22.*ac*xi2
      f3=as*(6.*a-2.*a1)*xs
      f4=-4.*ap*xs*xi2
      wk(l,i+1)=f1*(f2+f3+f4)/(a1c*a2*2.)
! c   this px fxxy
      f2=7.*a*a2-27.*as+22.*ac*xi2
      f3=as*(6.*a-2.*a2)*xs
      wk(i+1,l)=f1*(f2+f3+f4)/(a2c*a1*2.)
! c   this fxxy py
      f1=-asbc
      f2=7.*a-5.*a1-18.*as*xi2
      f3=(14.*a*a1-18*as)*xj2
      f4=44.*ac*xj2*xi2
      f5=4.*as*(a-a1)*xs*xj2
      f6=4.*ac*xs*xi2
      f7=-8.*ap*xs*xi2*xj2
      f8=-2.*a*(a-a1)*xs
      wk(l,j+1)=f1*(f2+f3+f4+f5+f6+f7+f8)/(4.*a1c*a2)
! c this py fxxy
      f1=-asbc
      f2=7.*a-5.*a2-18.*as*xi2
      f3=(14.*a*a2-18*as)*xj2
      f4=44.*ac*xj2*xi2
      f5=4.*as*(a-a2)*xs*xj2
      f8=-2.*a*(a-a2)*xs
      wk(j+1,l)=f1*(f2+f3+f4+f5+f6+f7+f8)/(4.*a2c*a1)
! c   this fxxy dxx
      f1=-asbc*x(j)
      f6=8.*at*xs*xi4
      f2=7.*aa1-5.*a12+7.*aa2-27.*as
      f3=(132.*ac-18.*as*a1-18.*as*a2)*xi2
      f4=(6.*ac-2.*as*a1-2.*as*a2+2.*a*a12)*xs
      f5=-(24.*ap-4.*ac*a1-4.*ac*a2)*xs*xi2-52.*ap*xi4
      wk(l,i+4)=f1*(f2+f3+f4+f5+f6)/(4.*a1c*a2s)
! c   this dxx fxxy
      f1=asbc*x(j)
      wk(i+4,l)=f1*(f2+f3+f4+f5+f6)/(4.*a2c*a1s)
! c   this fxxy dyy
      f1=asbc*x(j)
      f2=27.*as-21.*aa1-7.*aa2+5.*a12
      f3=(18.*a1*as-22.*ac)*xj2
      f4=(18.*a2*as-66.*ac)*xi2
      f6=4.*ac*(a-a1)*xs*xj2
      f7=(12.*ap-4.*ac*a2)*xs*xi2
      f8=52.*ap*xj2*xi2
      f9=-8.*at*xs*xi2*xj2
      f5=-(6.*ac-6.*as*a1-2.*as*a2+2.*a*a12)*xs
      wk(l,j+4)=f1*(f2+f3+f4+f5+f6+f7+f8+f9)/(4.*a1c*a2s)
! c   this dyy fxxy
      f1=-asbc*x(j)
      f2=27.*as-21.*aa2-7.*aa1+5.*a12
      f3=(18.*a2*as-22.*ac)*xj2
      f4=(18.*a1*as-66.*ac)*xi2
      f6=4.*ac*(a-a2)*xs*xj2
      f7=(12.*ap-4.*ac*a1)*xs*xi2
      f5=-(6.*ac-6.*as*a2-2.*as*a1+2.*a*a12)*xs
      wk(j+4,l)=f1*(f2+f3+f4+f5+f6+f7+f8+f9)/(4.*a2c*a1s)
! c   this fxxy dxy
      if(i.ne.3.and.j.ne.3) m=8
      if(i.ne.1.and.j.ne.1) m=10
      if(i.ne.2.and.j.ne.2) m=9
      f1=asbc*x(i)
      f2=27.*as-7.*a*a1-22.*ac*xi2
      f3=-as*(6.*a-2.*a1)*xs
      f4=4.*ap*xs*xi2
      f5=2.*ac*(6.*a-2.*a1)*xs*xj2
      f6=(18.*a1*as-66.*ac)*xj2
      f7=52.*ap*xi2*xj2
      f8=-8.*at*xs*xi2*xj2
      wk(l,m)=f1*(f2+f3+f4+f5+f6+f7+f8)/(4.*a1c*a2s)
! c   this dxy fxxy
      f1=-asbc*x(i)
      f2=27.*as-7.*a*a2-22.*ac*xi2
      f3=-as*(6.*a-2.*a2)*xs
      f4=4.*ap*xs*xi2
      f5=2.*ac*(6.*a-2.*a2)*xs*xj2
      f6=(18.*a2*as-66.*ac)*xj2
      wk(m,l)=f1*(f2+f3+f4+f5+f6+f7+f8)/(4.*a2c*a1s)
! c   this fxxy fxxy
      f1=asbc
      f2=27.*as-2.*a12+(18.*a*a12-132.*ac)*xi2
      f3=-6.*ac*xs+(4.*a*a12-66.*ac)*xj2
      f4=-(44.*as*a12-312.*ap)*xi2*xj2
      f5=12.*ap*xs*xj2
      f6=-(4.*as*a12-24.*ap)*xs*xi2
      f7=2.*a*(4.*as*a12-24.*ap)*xs*xi2*xj2
      f8=52.*ap*xi4-8.*at*xs*xi4
      f9=16.*a**6*xs*xi4*xj2
      f10=-120.*at*xi4*xj2
      wk(l,l)=f1*(f2+f3+f4+f5+f6+f7+f8+f9+f10)/(8.*a1c*a2c)
! c this fxxx py
      f1=-acbc*x(i)*x(j)
      f2=21.*a1-27.*a+(6.*as-6.*a*a1)*xs
      f3=as*xi2*(22.-4.*axs)
      wk(i+17,j+1)=f1*(f2+f3)/(a1c*a2*2.)
! c   this py fxxx
      f2=21.*a2-27.*a+(6.*as-6.*a*a2)*xs
      wk(j+1,i+17)=f1*(f2+f3)/(a2c*a1*2.)
! c   this fxxx dyy
      f1=asbc*x(i)
      f2=27.*as+15.*a12-21.*a*a2-21.*a*a1
      f3=(18.*a2*as-22.*ac)*xi2
      f4=-(-54.*a1*as+66.*ac)*xj2
      f5=-3.*a*(2.*as-2.*a*a2-2.*a*a1+2.*a12)*xs
      f6=4.*ac*(a-a2)*xs*xi2
      f7=(12.*ap-12.*ac*a1)*xs*xj2
      f8=52.*ap*xi2*xj2
      f9=-8.*at*xs*xi2*xj2
      wk(i+17,j+4)=f1*(f2+f3+f4+f5+f6+f7+f8+f9)/(4.*a1c*a2s)
! c   this dyy fxxx
      f1=-asbc*x(i)
      f3=(18.*a1*as-22.*ac)*xi2
      f4=-(-54.*a2*as+66.*ac)*xj2
      f6=4.*ac*(a-a1)*xs*xi2
      f7=(12.*ap-12.*ac*a2)*xs*xj2
      wk(j+4,i+17)=f1*(f2+f3+f4+f5+f6+f7+f8+f9)/(4.*a2c*a1s)
! c this fxxx dxy
      if(i.ne.3.and.j.ne.3) ms=8
      if(i.ne.2.and.j.ne.2) ms=9
      if(i.ne.1.and.j.ne.1) ms=10
      f1=-asbc*x(j)
      f2=21.*a1*a-27.*as+(6.*ac-6.*a1*as)*xs
      f3=(132.*ac-54.*as*a1)*xi2
      f4=-4.*as*(6.*as-3.*a*a1)*xs*xi2
      f5=-52.*ap*xi4
      f6=8.*at*xs*xi4
      wk(i+17,ms)=f1*(f2+f3+f4+f5+f6)/(4.*a1c*a2s)
! c   this dxy fxxx
      f2=21.*a2*a-27.*as+(6.*ac-6.*a2*as)*xs
      f3=(132.*ac-54.*as*a2)*xi2
      f4=-4.*as*(6.*as-3.*a*a2)*xs*xi2
      wk(ms,i+17)=-f1*(f2+f3+f4+f5+f6)/(4.*a2c*a1s)
! c   this fxxx fxxy
      if(i.eq.1.and.j.eq.2) mt=12
      if (i.eq.1.and.j.eq.3) mt=13
      if (i.eq.2.and.j.eq.1) mt=14
      if (i.eq.2.and.j.eq.3) mt=15
      if (i.eq.3.and.j.eq.1) mt=16
      if (i.eq.3.and.j.eq.2) mt=17
      f1=-acbc*x(i)*x(j)
      f2=21.*a12-81.*aa1-27.*aa2+165.*as
      f3=-(260.*ac-66.*as*a1-22.*as*a2)*xi2
      f4=-(30.*ac-18.*as*a1-6.*as*a2+6.*a*a12)*xs
      f5=(40.*ap-12.*ac*a1-4.*ac*a2)*xs*xi2
      f6=60.*ap*xi4
      f7=-8.*at*xi4*xs
      wk(i+17,l)=f1*(f2+f3+f4+f5+f6+f7)/(4.*a1c*a2c)
! c   this fxxy fxxx
      f2=21.*a12-81.*aa2-27.*aa1+165.*as
      f3=-(260.*ac-66.*as*a2-22.*as*a1)*xi2
      f4=-(30.*ac-18.*as*a2-6.*as*a1+6.*a*a12)*xs
      f5=(40.*ap-12.*ac*a2-4.*ac*a1)*xs*xi2
      wk(l,i+17)=f1*(f2+f3+f4+f5+f6+f7)/(4.*a2c*a1c)
! c   this fxxx fxyy
      if(i.eq.1.and.j.eq.2) mm=14
      if(i.eq.1.and.j.eq.3) mm=16
      if(i.eq.2.and.j.eq.1) mm=12
      if(i.eq.2.and.j.eq.3) mm=17
      if(i.eq.3.and.j.eq.1) mm=13
      if(i.eq.3.and.j.eq.2) mm=15
      f1=asbc
      f2=27.*as-21.*aa2-21.*aa1+15.*a12
      f3=(54.*as*a1+108.*as*a2-132.*ac-42.*a*a12)*xi2
      f4=(54.*a1*as-66.*ac)*xj2
      f5=-3.*a*(2.*as-2.*a*a2-2.*a*a1+2.*a12)*xs
      f6=(24.*ap-24.*ac*a2-12.*ac*a1+12.*as*a12)*xs*xi2
      f7=(12.*ap-12.*ac*a1)*xs*xj2
      f8=(312.*ap-132.*ac*a1)*xi2*xj2
      f9=-(48.*at-24.*ap*a1)*xs*xi2*xj2
      f10=(52.*ap-44.*ac*a2)*xi4
      f11=-120.*at*xi4*xj2
      f12=-8.*ap*(a-a2)*xs*xi4
      f13=0.
      f14=16.*a**6*xs*xi4*xj2
      wk(i+17,mm)=f1*(f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14) &
     &  /(8.*a1c*a2c)
! c   this fxyy fxxx
      f3=(54.*as*a2+108.*as*a1-132.*ac-42.*a*a12)*xi2
      f4=(54.*a2*as-66.*ac)*xj2
      f6=(24.*ap-24.*ac*a1-12.*ac*a2+12.*as*a12)*xs*xi2
      f7=(12.*ap-12.*ac*a2)*xs*xj2
      f8=(312.*ap-132.*ac*a2)*xi2*xj2
      f9=-(48.*at-24.*ap*a2)*xs*xi2*xj2
      f10=(52.*ap-44.*ac*a1)*xi4
      f12=-8.*ap*(a-a1)*xs*xi4
      wk(mm,i+17)=f1*(f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14) &
     &  /(8.*a2c*a1c)
  210 continue
! c this loop has three symmentry elements for x,y,z
      do 215 i=1,2
      xi2=x(i)**2
      do 215 j=2,3
      xj2=x(j)**2
      if (i.eq.j) go to 215
      k=i+j-2
      go to (21,22,23), k
   21 ma=12
      mb=14
      na=18
      nb=19
      go to 141
   22 ma=13
      mb=16
      na=18
      nb=20
      go to 141
   23 ma=15
      mb=17
      na=19
      nb=20
      go to 141
  141 continue
! c  this fxxy fxyy
      f1=-acbc*x(i)*x(j)
      f2=99.*as-27.*a*a1-27.*a*a2+7.*a12
      f3=(22.*a2*as-78.*ac)*xi2
      f4=(22.*a1*as-78.*ac)*xj2
      f5=-(18.*ac-6.*as*a2-6.*as*a1+2.*a*a12)*xs
      f6=(12.*ap-4.*ac*a2)*xs*xi2
      f7=60.*ap*xi2*xj2
      f8=(12.*ap-4.*ac*a1)*xs*xj2
      f9=-8.*at*xs*xi2*xj2
      wk(ma,mb)=f1*(f2+f3+f4+f5+f6+f7+f8+f9)/(4.*a1c*a2c)
! c   this fxyy fxxy
      f3=(22.*a1*as-78.*ac)*xi2
      f4=(22.*a2*as-78.*ac)*xj2
      f6=(12.*ap-4.*ac*a1)*xs*xi2
      f8=(12.*ap-4.*ac*a2)*xs*xj2
      wk(mb,ma )=f1*(f2+f3+f4+f5+f6+f7+f8+f9)/(4.*a2c*a1c)
! c   this fxxx fyyy
      f1=-acbc*x(i)*x(j)
      f2=99.*as-81.*aa1-81.*aa2+63.*a12
      f3=(66.*as*a2-78.*ac)*xi2
      f4=(66.*as*a1-78.*ac)*xj2
      f5=-(18.*ac-18.*as*a2-18.*as*a1+18.*a*a12)*xs
      f6=(12.*ap-12.*ac*a2)*xs*xi2
      f7=(12.*ap-12.*ac*a1)*xs*xj2
      f8=60.*ap*xi2*xj2-8.*at*xs*xi2*xj2
      wk(na,nb)=f1*(f2+f3+f4+f5+f6+f7+f8)/(4.*a1c*a2c)
! c   this fyyy fxxx
      f3=(66.*as*a1-78.*ac)*xi2
      f4=(66.*as*a2-78.*ac)*xj2
      f6=(12.*ap-12.*ac*a1)*xs*xi2
      f7=(12.*ap-12.*ac*a2)*xs*xj2
      wk(nb,na)=f1*(f2+f3+f4+f5+f6+f7+f8)/(4.*a2c*a1c)
  215 continue
      return
      end subroutine kinef1
! c
!*****************************************************
! c
      subroutine kinef2 (a1,a2,r,t,wk)

      use O_Kinds
      use O_Constants

      implicit none
      real (kind=double) :: a1, a2
      real (kind=double), dimension (3) :: r,t,x
      real (kind=double), dimension (20,20) :: wk

      integer :: i,j,k,l,mq,mr,lo,ko,lp,kp
      real (kind=double) :: a, as, ac, ap, at, a1s, a2s, a1c, a12, a2c
      real (kind=double) :: b, aa1, aa2, asa1, xs, axs, c, abc, asbc
      real (kind=double) :: a1sa2s, a1ca2c, a1ca2s, xsa2, asxs, acbc
      real (kind=double) :: xi2, xj2, xk2
      real (kind=double) :: f1, f2, f3, f4, f5, f6, f7, f8, f9, f10
      real (kind=double) :: f11, f12, f13, f14, f15

      x(1)=t(1)-r(1)
      x(2)=t(2)-r(2)
      x(3)=t(3)-r(3)
      a=a1*a2/(a1+a2)
      as=a*a
      ac=as*a
      ap=ac*a
      at=ap*a
      a1s=a1**2
      a2s=a2**2
      a1c=a1**3
      a2c=a2**3
      a12=a1*a2
      b=pi/(a1+a2)
      b=dsqrt(b)
      b=b*b*b
      aa1=a*a1
      aa2=a*a2
      asa1=as*a1
      a1sa2s=a1s*a2s
      a1ca2c=a1c*a2c
      a1ca2s=a1c*a2s
      xs=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
      axs=a*xs
      xsa2=xs*a2
      asxs=as*xs
      c=a*xs
      c=dexp(-c)
      abc=a*b*c
      asbc=as*b*c
      acbc=ac*b*c
      i=1
      j=2
      k=3
      xi2=x(i)*x(i)
      xj2=x(j)*x(j)
      xk2=x(k)*x(k)
! c  this loop only for xyz
! c   this fxyz s
      f1=ap*b*c*x(i)*x(j)*x(k)
      f2=(9.-2.*axs)
      wk(11,1)=f1*f2/a1c
! c   this s fxyz
      wk(1,11)=-f1*f2/a2c
! c   this fxyz fxyz
      f1=-ap*b*c
      f2=-9.+22.*a*xi2+22.*a*xj2
      f3=2.*axs-4.*as*xi2*xs
      f4=-4.*as*xj2*xs
      f5=-52.*as*xi2*xj2
      f6=(22.*a)*xk2
      f7=-52.*as*xk2*xi2
      f8=-52.*as*xk2*xj2
      f9=-4.*as*xk2*xs
      f10=120.*ac*xk2*xi2*xj2
      f11=8.*ac*xk2*xi2*xs
      f12=8.*ac*xk2*xj2*xs
      f13=8.*ac*xi2*xj2*xs
      f14=-16.*ap*xi2*xj2*xk2*xs
      wk(11,11)=f1*(f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14) &
     &  /(8.*a1c*a2c)
! c   continue
! c  this loop has exchange in x-y-z,y-x,z-x and three symmetry
      do 220 i=1,3
      xi2=x(i)*x(i)
      do 220 j=1,2
      xj2=x(j)*x(j)
      do 220 k=2,3
      xk2=x(k)*x(k)
      if (i.ne.j.and.j.ne.k.and.k.ne.i) go to 221
      go to 220
  221 continue
! c   this fxyz px
      f1=-ap*b*c*x(j)*x(k)
      f2=22.*a*xi2+2.*axs
      f3=-4.*as*xi2*xs-9.
      wk(11,i+1)=f1*(f2+f3)/(2.*a1c*a2)
! c   this px fxyz
      wk(i+1,11)=f1*(f2+f3)/(2.*a2c*a1)
! c   this fxyz dxx
      f1=ap*b*c*x(i)*x(j)*x(k)
      f2=9.*a2-33.*a+26.*as*xi2
      f3=(6.*as-2.*a*a2)*xs
      f4=-4.*ac*xs*xi2
      wk(11,i+4)=f1*(f2+f3+f4)/(2.*a1c*a2s)
! c   this dxx fxyz
      f2=9.*a1-33.*a+26.*as*xi2
      f3=(6.*as-2.*a*a1)*xs
      wk(i+4,11)=-f1*(f2+f3+f4)/(2.*a2c*a1s)
! c   this fxyz dyz
      f1=ap*b*c*x(i)
      f2=9.-22.*a*xk2-22.*a*xj2
      f3=-2.*axs+4.*as*xj2*xs
      f4=4.*as*xk2*xs
      f5=52.*as*xk2*xj2
      f6=-8.*ac*xk2*xj2*xs
      wk(11,11-i)=f1*(f2+f3+f4+f5+f6)/(4.*a1c*a2s)
! c   this dyz fxyz
      wk(11-i,11)=-f1*(f2+f3+f4+f5+f6)/(4.*a2c*a1s)
! c   this fxyz fxxx
      f1=acbc*x(j)*x(k)
      f2=27.*a*a2-33.*as+(6.*ac-6.*as*a2)*xs
      f3=(156.*ac-66.*as*a2)*xi2
      f4=-4.*as*(6.*as-3.*a*a2)*xs*xi2
      f5=-60.*ap*x(i)**4+8.*at*xs*x(i)**4
      wk(11,i+17)=f1*(f2+f3+f4+f5)/(4.*a1c*a2c)
! c   this fxxx fxyz
      f2=27.*a*a1-33.*as+(6.*ac-6.*as*a1)*xs
      f3=(156.*ac-66.*as*a1)*xi2
      f4=-4.*as*(6.*as-3.*a*a1)*xs*xi2
      wk(i+17,11)=f1*(f2+f3+f4+f5)/(4.*a2c*a1c)
! c   this fxxx dyz
      f1=asbc*x(i)*x(j)*x(k)
      f2=27.*as*a1-33.*ac+as*(6.*as-6.*a*a1)*xs
      f3=26.*ap*xi2-4.*at*xi2*xs
      f4=0.
      wk(i+17,11-i)=f1*(f2+f3+f4)/(a1c*a2s*2.)
! c   this dyz f xxx
      f2=27.*as*a2-33.*ac+as*(6.*as-6.*a*a2)*xs
      wk(11-i,i+17)=-f1*(f2+f3+f4)/(a2c*a1s*2.)
  220 continue
! c  this loop has exchange x-y-z and x-y,y-z,z-x and six symmentry
      l=11
      do 230 i=1,3
      xi2=x(i)*x(i)
      do 230 j=1,3
      xj2=x(j)*x(j)
      do 230 k=1,3
      xk2=x(k)*x(k)
      if (i.ne.j.and.j.ne.k.and.k.ne.i) go to 231
      go to 230
  231 continue
      l=l+1
! c   this fxxy pz
      f1=-acbc*x(j)*x(k)
      f2=7.*a1-9.*a+22.*as*xi2
      f3=2.*a*(a-a1)*xs-4.*ac*xs*xi2
      wk(l,k+1)=f1*(f2+f3)/(a1c*a2*2.)
! c   this pz fxxy
      f2=7.*a2-9.*a+22.*as*xi2
      f3=2.*a*(a-a2)*xs-4.*ac*xs*xi2
      wk(k+1,l)=f1*(f2+f3)/(a2c*a1*2.)
! c   this fxxy dyz
      f1=asbc*x(k)
      f2=9.*as-7.*a*a1-22.*ac*xi2
      f3=(18.*as*a1-22.*ac)*xj2
      f4=52.*ap*xj2*xi2
      f5=4.*ac*(a-a1)*xs*xj2
      f6=4.*ap*xs*xi2-8.*at*xs*xi2*xj2
      f7=-2.*as*(a-a1)*xs
      wk(l,11-i)=f1*(f2+f3+f4+f5+f6+f7)/(4.*a1c*a2s )
! c   this dyz fxxy
      f2=9.*as-7.*a*a2-22.*ac*xi2
      f3=(18.*as*a2-22.*ac)*xj2
      f5=4.*ac*(a-a2)*xs*xj2
      f7=-2.*as*(a-a2)*xs
      wk(11-i,l)=-f1*(f2+f3+f4+f5+f6+f7)/(4.*a2c*a1s)
! c   this fxxy dzx
      if (j.eq.3) mq=8
      if (j.eq.1) mq=10
      if (j.eq.2) mq=9
      f1=asbc*x(i)*x(j)*x(k)
      f2=9.*as*a1-33.*ac+26.*ap*xi2
      f3=ac*(6.*a-2.*a1)*xs-4.*at*xs*xi2
      wk(l,mq)=f1*(f2+f3)/(a1c*a2s*2.)
! c   this dzx fxxy
      f2=9.*as*a2-33.*ac+26.*ap*xi2
      f3=ac*(6.*a-2.*a2)*xs-4.*at*xs*xi2
      wk(mq,l)=-f1*(f2+f3)/(a2c*a1s*2.)
! c   this fxxy dzz
      f1=asbc*x(j)
      f2=9.*as-7.*aa1-7.*aa2+5.*a12+(18.*as*a1-22.*ac)*xk2
      f3=(18.*as*a2-22.*ac)*xi2
      f4=(2.*as*a1+2.*as*a2-2.*a*a12-2.*ac)*xs+4.*ac*(a-a1)*xs*xk2
      f5=4.*ac*(a-a2)*xs*xi2
      f6=52.*ap*xi2*xk2
      f7=-8.*at*xs*xi2*xk2
      wk(l,k+4)=f1*(f2+f3+f4+f5+f6+f7)/(4.*a1c*a2s)
! c   this dzz fxxy
      f2=9.*as-7.*aa1-7.*aa2+5.*a12+(18.*as*a2-22.*ac)*xk2
      f3=(18.*as*a1-22.*ac)*xi2
      f4=(2.*as*a1+2.*as*a2-2.*a*a12-2.*ac)*xs+4.*ac*(a-a2)*xs*xk2
      f5=4.*ac*(a-a1)*xs*xi2
      wk(k+4,l)=-f1*(f2+f3+f4+f5+f6+f7)/(4.*a2c*a1s)
! c   this fxxy fxyz
      f1=-acbc*x(i)*x(k)
      f2=33.*as-9.*a*a1-26.*ac*xi2
      f3=(22.*a1*as-78.*ac)*xj2
      f4=-as*(6.*a-2.*a1)*xs+4.*ap*xs*xi2
      f5=2.*ac*(6.*a-2.*a1)*xs*xj2
      f6=60.*ap*xi2*xj2
      f7=-8.*at*xs*xi2*xj2
      wk(l,11)=f1*(f2+f3+f4+f5+f6+f7)/(4.*a1c*a2c)
! c   this fxyz fxxy
      f2=33.*as-9.*a*a2-26.*ac*xi2
      f3=(22.*a2*as-78.*ac)*xj2
      f4=-as*(6.*a-2.*a2)*xs+4.*ap*xs*xi2
      f5=2.*ac*(6.*a-2.*a2)*xs*xj2
      wk(11,l)=f1*(f2+f3+f4+f5+f6+f7)/(4.*a2c*a1c)
! c   this fxxy fyyz
      if(i.eq.1.and.j.eq.2) mr=15
      if(i.eq.1.and.j.eq.3) mr=17
      if(i.eq.2.and.j.eq.1) mr=13
      if(i.eq.2.and.j.eq.3) mr=16
      if(i.eq.3.and.j.eq.1) mr=12
      if(i.eq.3.and.j.eq.2) mr=14
      f1=-acbc*x(j)*x(k)
      f2=33.*as-27.*aa1-9.*aa2+7.*a12+(22.*a1*as-26.*ac)*xj2
      f3=(22.*a2*as-78.*ac)*xi2
      f4=-(6.*ac-6.*as*a1-2.*as*a2+2.*a*a12)*xs
      f5=4.*ac*(a-a1)*xs*xj2
      f6=(12.*ap-4.*ac*a2)*xs*xi2
      f7=60.*ap*xj2*xi2
      f8=-8.*at*xs*xi2*xj2
      wk(l,mr)=f1*(f2+f3+f4+f5+f6*f7+f8)/(4.*a1c*a2c)
! c   this fyyz fxxy
      f2=33.*as-27.*aa2-9.*aa1+7.*a12+(22.*a2*as-26.*ac)*xj2
      f3=(22.*a1*as-78.*ac)*xi2
      f4=-(6.*ac-6.*as*a2-2.*as*a1+2.*a*a12)*xs
      f5=4.*ac*(a-a2)*xs*xj2
      f6=(12.*ap-4.*ac*a1)*xs*xi2
      wk(mr,l)=f1*(f2+f3+f4+f5+f6*f7+f8)/(4.*a2c*a1c)
! c   thjis fzzz fxxy
      f1=-acbc*x(k)*x(j)
      f2=33.*as+21.*a12-27.*a*a2-27.*a*a1
      f3=(22.*as*a2-26.*ac)*xk2
      f4=-(-66.*as*a1+78.*ac)*xi2
      f5=-3.*a*(2.*as-2.*a*a2-2.*a*a1+2.*a12)*xs
      f6=4.*ac*(a-a2)*xs*xk2
      f7=(12.*ap-12.*ac*a1)*xs*xi2
      f8=60.*ap*xk2*xi2
      f9=-8.*at*xs*xk2*xi2
      wk(k+17,l)=f1*(f2+f3+f4+f5+f6+f7+f8+f9)/(4.*a1c*a2c)
! c   this fxxy fzzz
      f1=-acbc*x(k)*x(j)
      f3=(22.*as*a1-26.*ac)*xk2
      f4=-(-66.*as*a2+78.*ac)*xi2
      f6=4.*ac*(a-a1)*xs*xk2
      f7=(12.*ap-12.*ac*a2)*xs*xi2
      wk(l,k+17)=f1*(f2+f3+f4+f5+f6+f7+f8+f9)/(4.*a2c*a1c)
  230 continue
! c  this loop has three symmentry elements in x,y,z
      do 250 i=1,3
      xi2=x(i)*x(i)
      go to (1,2,3),  i
    1 j=2
      xj2=x(j)*x(j)
      k=3
      xk2=x(k)*x(k)
      go to 240
    2 j=3
      xj2=x(j)*x(j)
      k=1
      xk2=x(k)*x(k)
      go to 240
    3 j=1
      xj2=x(j)*x(j)
      k=2
      xk2=x(k)*x(k)
      go to 240
  240 continue
      go to (11,12,13), i
   11 lo=12
      ko=17
      lp=12
      kp=13
      go to 241
   12 lo=15
      ko=13
      lp=15
      kp=14
      go to 241
   13 lo=16
      ko=14
      lp=16
      kp=17
      go to 241
  241 continue
      xj2=x(j)*x(j)
      xk2=x(k)*x(k)
! c   this fxxy fyzz
      f1=-asbc
      f2=7.*aa1+7.*aa2-5.*a12-9.*as+(22.*ac+14.*a*a12-18.*as*a1 &
     &   -18.*as*a2)*xj2
      f3=(44.*ac*a1-52.*ap)*xj2*xk2
      f4=(44.*ac*a2-52.*ap)*xi2*xj2
      f5=-(18.*as*a1-22.*ac)*xk2
      f6=-(18.*as*a2-22.*ac)*xi2
      f7=2.*ac*xs-4.*ap*xs*xj2
      f8=-4.*ac*(a-a2)*xs*xi2
      f9=-4.*ac*(a-a1)*xs*xk2
      f10=-52.*ap*xi2*xk2
      f11=8.*ap*(a-a1)*xs*xk2*xj2
      f12=8.*ap*(a-a2)*xs*xi2*xj2
      f13=8.*at*xs*xi2*xk2
      f14=120.*at*xi2*xj2*xk2
      f15=-16.*a**6*xs*xi2*xj2*xk2
      wk(lo,ko)=f1*(f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14+f15 &
     &    )/(8.*a1c*a2c)
! c this  fyzz fxxy
      f3=(44.*ac*a2-52.*ap)*xj2*xk2
      f4=(44.*ac*a1-52.*ap)*xi2*xj2
      f5=-(18.*as*a2-22.*ac)*xk2
      f6=-(18.*as*a1-22.*ac)*xi2
      f8=-4.*ac*(a-a1)*xs*xi2
      f9=-4.*ac*(a-a2)*xs*xk2
      f11=8.*ap*(a-a2)*xs*xk2*xj2
      f12=8.*ap*(a-a1)*xs*xi2*xj2
      wk(ko,lo)=f1*(f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14+f15 &
     &    )/(8.*a2c*a1c)
! c   this fzxx fxxy
      f1=acbc*x(k)*x(j)
      f2=9.*aa1+9.*aa2-7.*a12-33.*as+(156.*ac-22.*as*a1-22.*as*a2)*xi2
      f3=(6.*ac-2.*as*a1-2.*as*a2+2.*a*a12)*xs+(4.*ac*a1+4.*ac*a2-24. &
     &  *ap)*xs*xi2
      f4=-60.*ap*x(i)**4+8.*at*xs*x(i)**4
      wk(kp,lp)=f1*(f2+f3+f4)/(4.*a1c*a2c)
! c   this fxxy fzxx
      wk(lp,kp)=f1*(f2+f3+f4)/(4.*a2c*a1c)
  250 continue
      return
      end subroutine kinef2


      subroutine overlapInteg(g,l1,l2,al1,al2,b,c)

      use O_Kinds
      use O_Constants

      implicit none
! c
! c     g(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
! c     9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
! c     14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy
! c
! c     wo(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
! c     10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
! c     18,xxx; 19,yyy; 20,zzz
! c 
      ! define the dummy variables passed to this subroutine.
      real (kind=double), dimension (16,16) :: g
      integer :: l1, l2
      real (kind=double) :: al1
      real (kind=double) :: al2
      real (kind=double), dimension (3) :: b
      real (kind=double), dimension (3) :: c
      integer :: option

      ! define local variables
      real (kind=double) :: als, ex, a
      real (kind=double), dimension (20,20) :: wo
      integer :: nff, ndd
      integer :: i,j,k,l,m1,m2,m3,m,ms,mt,mm,ma,mb,na,nb,mq,mr
      integer :: lo,ko,lp,kp

      real (kind=double), dimension (3) :: ae,be
      real (kind=double) :: at,ey,ez,ccc,coef,b2,c2,e2
! c

 100  continue
! c
      if(l1.eq.16.or.l2.eq.16) then
         nff=0
      else
         nff=1
      endif
      g(:,:)=0.0_double
      wo(:,:)=0.0_double
      ndd=0
      at=al1+al2
      ex=(al1*b(1)+al2*c(1))/at
      ey=(al1*b(2)+al2*c(2))/at
      ez=(al1*b(3)+al2*c(3))/at
      e2=ex**2+ey**2+ez**2
      b2=b(1)**2+b(2)**2+b(3)**2
      c2=c(1)**2+c(2)**2+c(3)**2
      ccc=at*e2-al1*b2-al2*c2
      coef=dsqrt(pi/at)**3*dexp(ccc)
      ae(1)=ex-b(1)
      ae(2)=ey-b(2)
      ae(3)=ez-b(3)
      be(1)=ex-c(1)
      be(2)=ey-c(2)
      be(3)=ez-c(3)
! c this ss
      wo(1,1)=coef
      do 105 i=1,3
! c this is spx
      wo(1,i+1)=coef*be(i)
! c this is pxs
      wo(i+1,1)=coef*ae(i)
      if(ndd.ne.0) go to 1000
! c this is sdxx
      wo(1,i+4)=coef*(be(i)**2+.5d0/at)
! c this is dxxs
      wo(i+4,1)=coef*(ae(i)**2+.5d0/at)
 1000 continue
! c this is pxpx
      wo(i+1,i+1)=coef*(ae(i)*be(i)+.5d0/at)
      if(ndd.ne.0) go to 1001
! c this is dxxpx
      wo(i+4,i+1)=coef*(be(i)*ae(i)**2+.5d0*(2.d0*ae(i)+be(i))/at)
! c this is pxdxx
      wo(i+1,i+4)=coef*(ae(i)*be(i)**2+.5d0*(2.d0*be(i)+ae(i))/at)
! c this is dxxdxx
      wo(i+4,i+4)=coef*(ae(i)**2*be(i)**2+.5d0*(ae(i)**2+be(i)**2) &
     &/at+2.d0*be(i)*ae(i)/at+.75d0/at**2)
 1001 continue
  105 continue
      if(ndd.ne.0) go to 1002
      k=7
      do 110 i=1,2
      do 110 j=2,3
      if (i.eq.j) go to 110
      k=k+1
! c this is sdxy
      wo(1,k)=coef*be(i)*be(j)
! c this is dxys
      wo(k,1)=coef*ae(i)*ae(j)
  110 continue
 1002 continue
      do 115 i=1,3
      do 115 j=1,3
      if (i.eq.j) go to 115
! c this is pxpy
      wo(i+1,j+1)=coef*ae(i)*be(j)
      if(ndd.ne.0) go to 1003
! c this is pxdyy
      wo(i+1,j+4)=coef*ae(i)*(be(j)**2+.5d0/at)
! c this is dxxpy
      wo(i+4,j+1)=coef*be(j)*(ae(i)**2+.5d0/at)
! c this is dxxdyy
      wo(i+4,j+4)=coef*(ae(i)**2+.5d0/at)*(be(j)**2+.5d0/at)
 1003 continue
  115 continue
      if(ndd.ne.0) go to 1004
      do 130 i=1,3
      do 130 j=1,2
      do 130 k=2,3
      if (j.eq.k) go to 130
      if (i.ne.j.and.i.ne.k) go to 120
      if (j.eq.i) go to 118
      l=j
      go to 119
  118 l=k
  119 continue
! c this is pxdxy
      wo(i+1,j+k+5)=coef*be(l)*(ae(i)*be(i)+.5d0/at)
! c this is dxxdxy
      wo(i+4,j+k+5)=coef*be(l)*(ae(i)**2*be(i)+.5d0*(2.d0*ae(i)+be &
     &(i))/at)
! c this is dxypx
      wo(j+k+5,i+1)=coef*ae(l)*(ae(i)*be(i)+.5d0/at)
! c this is dxydxx
      wo(j+k+5,i+4)=coef*ae(l)*(ae(i)*be(i)**2+.5d0*(2.d0*be(i)+ae &
     &(i))/at)
      go to 130
  120 continue
! c this is pxdyz
      wo(i+1,j+k+5)=coef*be(j)*be(k)*ae(i)
! c this is dxxdyz
      wo(i+4,j+k+5)=coef*be(j)*be(k)*(ae(i)**2+.5d0/at)
! c this is dyzpx
      wo(j+k+5,i+1)=coef*ae(j)*ae(k)*be(i)
! c this is dyzdxx
      wo(j+k+5,i+4)=coef*ae(j)*ae(k)*(be(i)**2+.5d0/at)
  130 continue
      do 146 i=1,2
      do 146 j=2,3
      if (i.eq.j) go to 146
      do 145 k=1,2
      do 145 l=2,3
      if (k.eq.l) go to 145
      if (i.eq.k.and.j.eq.l) go to 140
      if (i.eq.k) go to 131
      if (i.eq.l) go to 132
      if (j.eq.k) go to 133
      if  (j.eq.l) go to 134
  131 m1=i
      m2=j
      m3=l
      go to 135
  132 m1=i
      m2=j
      m3=k
      go to 135
  133 m1=j
      m2=i
      m3=l
      go to 135
  134 m1=j
      m2=i
      m3=k
      go to 135
  135 continue
! c  this dxy dyz
      wo(i+j+5,k+l+5)=coef*ae(m2)*be(m3)*(ae(m1)*be(m1)+.5/at)
      go to 145
  140 continue
! c  this dxy dxy
      wo(i+j+5,k+l+5)=coef*(ae(i)*be(i)+.5/at)*(ae(j)*be(j)+.5/at)
  145 continue
  146 continue
 1004 continue
      if(nff.eq.1) go to 999
      do 205 i=1,3
! c   this loop has exchanging in x-y-z
! c   this fxxx s
      wo(i+17,1)=coef*ae(i)*(1.5/at+ae(i)**2)
! c   this s fxxx
      wo(1,i+17)=coef*be(i)*(1.5/at+be(i)**2)
! c   this fxxx px
      wo(i+17,i+1)=coef*(1.5*ae(i)**2/at+1.5*ae(i)*be(i) &
     &   /at+ae(i)**3*be(i)+.75/at**2)
! c   this px fxxx
      wo(i+1,i+17)=coef*(1.5*be(i)**2/at+1.5*be(i)*ae(i) &
     &   /at+be(i)**3*ae(i)+.75/at**2)
! c   this fxxx dxx
      wo(i+17,i+4)=coef*(9.*ae(i)/(4.*at**2)+3.*be(i)/(2.*at**2) &
     & +3.*ae(i)**2*be(i)/at+1.5*ae(i)*be(i)**2/at+.5*ae(i)**3 &
     & /at+be(i)**2*ae(i)**3)
! c   this dxx fxxx
      wo(i+4,i+17)=coef*(9.*be(i)/(4.*at**2)+3.*ae(i)/(2.*at**2) &
     &  +3.*be(i)**2*ae(i)/at+1.5*be(i)*ae(i)**2/at+.5*be(i) &
     &   **3/at+ae(i)**2*be(i)**3)
! c   this fxxx fxxx
      wo(i+17,i+17)=coef*((be(i)*27.*ae(i)+9.*be(i)**2)/(4.*at**2) &
     & +be(i)**2*(4.5*ae(i)**2+1.5*ae(i)*be(i))/at+ae(i)**3* &
     & be(i)*(1.5/at+be(i)**2)+9.*ae(i)**2/(4.*at**2)+15./(8.*at**3))
  205 continue
      l=11
      do 210 i=1,3
      do 210 j=1,3
      if (i.eq.j) go to 210
! c  this loop has exchanging in x-y-z,x-y,y-z,z-x,and six symmetry
! c  this fxxy s
      l=l+1
      wo(l,1)=coef*ae(j)*(ae(i)**2+.5/at)
! c  this s fxxy
      wo(1,l)=coef*be(j)*(be(i)**2+.5/at)
! c  this fxxy px
      wo(l,i+1)=coef*ae(j)*(.5*(2.*ae(i)+be(i))/at+be(i)*ae(i)**2)
! c  this px fxxy
      wo(i+1,l)=coef*be(j)*(.5*(2.*be(i)+ae(i))/at+ae(i)*be(i)**2)
! c  this fxxy py
      wo(l,j+1)=coef*(ae(i)**2+.5/at)*(ae(j)*be(j)+.5/at)
! c  this py fxxy
      wo(j+1,l)=coef*(be(i)**2+.5/at)*(be(j)*ae(j)+.5/at)
! c  this fxxy dxx
      wo(l,i+4)=coef*ae(j)*(.75/at**2+2.*ae(i)*be(i)/at &
     &  +.5*(ae(i)**2+be(i)**2)/at+ae(i)**2*be(i)**2)
! c  this dxx fxxy
      wo(i+4,l)=coef*be(j)*(.75/at**2+2.*be(i)*ae(i)/at &
     &  +.5*(be(i)**2+ae(i)**2)/at+be(i)**2*ae(i)**2)
! c  this fxxy dyy
      wo(l,j+4)=coef*(.5/at+ae(i)**2)*(.5*ae(j)/at+ &
     &  be(j)/at+ae(j)*be(j)**2)
! c  this dyy fxxy
      wo(j+4,l)=coef*(.5/at+be(i)**2)*(.5*be(j)/at+ &
     &  ae(j)/at+be(j)*ae(j)**2)
! c  this fxxy dxy
      if (i.ne.3.and.j.ne.3) m=8
      if (i.ne.1.and.j.ne.1) m=10
      if (i.ne.2.and.j.ne.2) m=9
      wo(l,m)=coef*(ae(j)*be(j)*(ae(i)/at+.5*be(i)/at+be(i)* &
     &  ae(i)**2)+(.5*ae(i)+.25*be(i))/at**2+.5*be(i)*ae(i)**2/at)
! c  this dxy fxxy
      wo(m,l)=coef*(be(j)*ae(j)*(be(i)/at+.5*ae(i)/at+ae(i)* &
     &  be(i)**2)+(.5*be(i)+.25*ae(i))/at**2+.5*ae(i)*be(i)**2/at)
! c  this fxxy fxxy
      wo(l,l)=coef*(.5/at+ae(j)*be(j))*(.75/at**2+2.*ae(i)*be(i)/ &
     &  at+.5*(be(i)**2+ae(i)**2)/at+ae(i)**2*be(i)**2)
! c  this fxxx py
      wo(i+17,j+1)=coef*ae(i)*be(j)*(ae(i)**2+1.5/at)
! c  this py fxxx
      wo(j+1,i+17)=coef*be(i)*ae(j)*(be(i)**2+1.5/at)
! c  this fxxx dyy
      wo(i+17,j+4)=coef*ae(i)*(.5/at+be(j)**2)*(1.5/at+ae(i)**2)
! c  this dyy fxxx
      wo(j+4,i+17)=coef*be(i)*(.5/at+ae(j)**2)*(1.5/at+be(i)**2)
! c  this fxxx dxy
      if (i.ne.3.and.j.ne.3) ms=8
      if (i.ne.2.and.j.ne.2) ms=9
      if (i.ne.1.and.j.ne.1) ms=10
      wo(i+17,ms)=coef*be(j)*(1.5*(ae(i)**2+ae(i)*be(i))/at+ &
     &  ae(i)**3*be(i)+.75/at**2)
! c  this dxy fxxx
      wo(ms,i+17)=coef*ae(j)*(1.5*(be(i)**2+be(i)*ae(i))/at+ &
     &  be(i)**3*ae(i)+.75/at**2)
! c  this fxxx fxxy
      if (i.eq.1.and.j.eq.2) mt=12
      if (i.eq.1.and.j.eq.3) mt=13
      if (i.eq.2.and.j.eq.1) mt=14
      if (i.eq.2.and.j.eq.3) mt=15
      if (i.eq.3.and.j.eq.1) mt=16
      if (i.eq.3.and.j.eq.2) mt=17
      wo(i+17,mt)=coef*be(j)*(9.*ae(i)/(4.*at**2)+1.5*be(i)/at**2 &
     &  +3.*ae(i)**2*be(i)/at+1.5*ae(i)*be(i)**2/at+.5*ae(i)**3/at &
     &  +be(i)**2*ae(i)**3)
! c  this fxxy fxxx
      wo(mt,i+17)=coef*ae(j)*(9.*be(i)/(4.*at**2)+1.5*ae(i)/at**2 &
     &  +3.*be(i)**2*ae(i)/at+1.5*be(i)*ae(i)**2/at+.5*be(i)**3/at &
     &  +ae(i)**2*be(i)**3)
! c  this fxxx fxyy
      if(i.eq.1.and.j.eq.2) mm=14
      if(i.eq.1.and.j.eq.3) mm=16
      if(i.eq.2.and.j.eq.1) mm=12
      if(i.eq.2.and.j.eq.3) mm=17
      if(i.eq.3.and.j.eq.1) mm=13
      if(i.eq.3.and.j.eq.2) mm=15
      wo(i+17,mm)=coef*(.5/at+be(j)**2)*((1.5/at+ae(i)**2)*(.5/at &
     &  +be(i)*ae(i))+ae(i)**2/at)
! c  this fxyy fxxx
      wo(mm,i+17)=coef*(.5/at+ae(j)**2)*((1.5/at+be(i)**2)*(.5/at &
     &  +ae(i)*be(i))+be(i)**2/at)
  210 continue
! c this loop has three symmentry elements for x,y,z
      do 215 i=1,2
      do 215 j=2,3
      if (i.eq.j) go to 215
      k=i+j-2
      go to (21,22,23), k
   21 ma=12
      mb=14
      na=18
      nb=19
      go to 141
   22 ma=13
      mb=16
      na=18
      nb=20
      go to 141
   23 ma=15
      mb=17
      na=19
      nb=20
      go to 141
  141 continue
! c  this fxxy fxyy
      wo(ma,mb)=coef*((.5*ae(j)+be(j))/at+ae(j)*be(j)**2)*(be(i)* &
     &  (.5/at+ae(i)**2)+ae(i)/at)
! c  this fxyy fxxy
      wo(mb,ma)=coef*((.5*be(j)+ae(j))/at+be(j)*ae(j)**2)*(ae(i)* &
     &  (.5/at+be(i)**2)+be(i)/at)
! c  this fxxx fyyy
      wo(na,nb)=coef*ae(i)*be(j)*(1.5/at+ae(i)**2)*(1.5/at+be(j)**2)
! c  this fyyy fxxx
      wo(nb,na)=coef*be(i)*ae(j)*(1.5/at+be(i)**2)*(1.5/at+ae(j)**2)
  215 continue
      i=1
      j=2
      k=3
! c  this loop only for xyz
! c   this fxyz s
      wo(11,1)=coef*ae(i)*ae(j)*ae(k)
! c   this s fxyz
      wo(1,11)=coef*be(i)*be(j)*be(k)
! c   this fxyz fxyz
      wo(11,11)=coef*(.5/at+ae(j)*be(j))*(.5/at+ae(k)*be(k)) &
     &  *(.5/at+ae(i)*be(i))
! c     continue
! c  this loop has exchange in x-y-z,y-x,z-x and three symmetry
      do 220 i=1,3
      do 220 j=1,2
      do 220 k=2,3
      if (i.ne.j.and.j.ne.k.and.k.ne.i) go to 221
      go to 220
  221 continue
! c   this fxyz px
      wo(11,i+1)=coef*ae(k)*ae(j)*(.5/at+ae(i)*be(i))
! c   this px fxyz
      wo(i+1,11)=coef*be(k)*be(j)*(.5/at+be(i)*ae(i))
! c   this fxyz dxx
      wo(11,i+4)=coef*ae(k)*ae(j)*(.5*(2.*be(i)+ae(i))/at+ae(i)* &
     &   be(i)**2)
! c   this dxx fxyz
      wo(i+4,11)=coef*be(k)*be(j)*(.5*(2.*ae(i)+be(i))/at+be(i)* &
     &   ae(i)**2)
! c   this fxyz dyz
      wo(11,11-i)=coef*ae(i)*(.5/at+ae(k)*be(k))*(.5/at+ae(j)*be(j))
! c   this dyz fxyz
      wo(11-i,11)=coef*be(i)*(.5/at+be(k)*ae(k))*(.5/at+be(j)*ae(j))
! c   this fxyz fxxx
      wo(11,i+17)=coef*ae(j)*ae(k)*((be(i)**2+1.5/at)*(ae(i)*be(i) &
     &  +.5/at)+be(i)**2/at)
! c   this fxxx fxyz
      wo(i+17,11)=coef*be(j)*be(k)*((ae(i)**2+1.5/at)*(be(i)*ae(i) &
     &  +.5/at)+ae(i)**2/at)
! c   this fxxx dyz
      wo(i+17,11-i)=coef*ae(i)*be(j)*be(k)*(ae(i)**2+1.5/at)
! c   this dyz fxxx
      wo(11-i,i+17)=coef*be(i)*ae(j)*ae(k)*(be(i)**2+1.5/at)
  220 continue
! c  this loop has exchange x-y-z and x-y,y-z,z-x and six symmentry
      l=11
      do 230 i=1,3
      do 230 j=1,3
      do 230 k=1,3
      if (i.ne.j.and.j.ne.k.and.k.ne.i) go to 231
      go to 230
  231 continue
      l=l+1
! c   this fxxy pz
      wo(l,k+1)=coef*ae(j)*be(k)*(.5/at+ae(i)**2)
! c   tihs pz fxxy
      wo(k+1,l)=coef*be(j)*ae(k)*(.5/at+be(i)**2)
! c   this fxxy dyz
      wo(l,11-i)=coef*be(k)*(.5/at+ae(i)**2)*(ae(j)*be(j)+.5/at)
! c   this dyz fxxy
      wo(11-i,l)=coef*ae(k)*(.5/at+be(i)**2)*(be(j)*ae(j)+.5/at)
! c   this fxxy dzx
      if (j.eq.3) mq=8
      if (j.eq.1) mq=10
      if (j.eq.2) mq=9
      wo(l,mq)=coef*ae(j)*be(k)*(.5*be(i)/at+ae(i)/at+be(i)*ae(i)**2)
! c   this dzx fxxy
      wo(mq,l)=coef*be(j)*ae(k)*(.5*ae(i)/at+be(i)/at+ae(i)*be(i)**2)
! c   this fxxy dzz
      wo(l,k+4)=coef*ae(j)*(.5/at+be(k)**2)*(.5/at+ae(i)**2)
! c   this dzz fxxy
      wo(k+4,l)=coef*be(j)*(.5/at+ae(k)**2)*(.5/at+be(i)**2)
! c   this fxyz fxxy
      wo(11,l)=coef*ae(k)*(.5/at+be(j)*ae(j))*(ae(i)*(.5/at+be(i)**2) &
     &   +be(i)/at)
! c   this fxxy fxyz
      wo(l,11)=coef*be(k)*(.5/at+ae(j)*be(j))*(be(i)*(.5/at+ae(i)**2) &
     &   +ae(i)/at)
! c   this fxxy fyyz
      if(i.eq.1.and.j.eq.2) mr=15
      if(i.eq.1.and.j.eq.3) mr=17
      if(i.eq.2.and.j.eq.1) mr=13
      if(i.eq.2.and.j.eq.3) mr=16
      if(i.eq.3.and.j.eq.1) mr=12
      if(i.eq.3.and.j.eq.2) mr=14
      wo(l,mr)=coef*be(k)*(.5/at+ae(i)**2)*(.5*ae(j)/at+be(j)/at &
     &  +ae(j)*be(j)**2)
! c   this fyyz fxxy
      wo(mr,l)=coef*ae(k)*(.5/at+be(i)**2)*(.5*be(j)/at+ae(j)/at &
     &  +be(j)*ae(j)**2)
! c   this fxxy fzzz
      wo(l,k+17)=coef*ae(j)*be(k)*(.5/at+ae(i)**2)*(1.5/at+be(k)**2)
! c   this fzzz fxxy
      wo(k+17,l)=coef*be(j)*ae(k)*(.5/at+be(i)**2)*(1.5/at+ae(k)**2)
  230 continue
! c  this loop has three symmentry elements in x,y,z
      do 250 i=1,3
      go to (1,2,3),  i
    1 j=2
      k=3
      go to 240
    2 j=3
      k=1
      go to 240
    3 j=1
      k=2
      go to 240
  240 continue
      go to (11,12,13), i
   11 lo=12
      ko=17
      lp=12
      kp=13
      go to 241
   12 lo=15
      ko=13
      lp=15
      kp=14
      go to 241
   13 lo=16
      ko=14
      lp=16
      kp=17
      go to 241
  241 continue
! c   this fxxy fyzz
      wo(lo,ko)=coef*(.5/at+be(k)**2)*(.5/at+ae(i)**2)*(.5/at+ae(j)*be &
     &     (j))
! c   this fyzz fxxy
      wo(ko,lo)=coef*(.5/at+ae(k)**2)*(.5/at+be(i)**2)*(.5/at+be(j) &
     &     *ae(j))
! c   this fxxy fzxx
      wo(lp,kp)=coef*be(k)*ae(j)*(.75/at**2+2.*ae(i)*be(i)/at+.5* &
     &  (be(i)**2+ae(i)**2)/at+ae(i)**2*be(i)**2)
! c   this fzxx fxxy
      wo(kp,lp)=coef*ae(k)*be(j)*(.75/at**2+2.*be(i)*ae(i)/at+.5* &
     &  (ae(i)**2+be(i)**2)/at+be(i)**2*ae(i)**2)
  250 continue
  999 continue
      continue
      g(1,1)=wo(1,1)
      if(l2.eq.1) goto 2009
      g(1,2)=wo(1,2)
      g(1,3)=wo(1,3)
      g(1,4)=wo(1,4)
      if(l2.eq.4) goto 2009
      g(1,5)=wo(1,8)
      g(1,6)=wo(1,9)
      g(1,7)=wo(1,10)
      g(1,8)=wo(1,5)-wo(1,6)
      g(1,9)=2.0*wo(1,7)-wo(1,5)-wo(1,6)
      if(l2.eq.9) goto 2009
      g(1,10)=wo(1,11)
      g(1,11)=wo(1,13)-wo(1,15)
      g(1,12)=wo(1,18)-3.0*wo(1,14)
      g(1,13)=3.0*wo(1,12)-wo(1,19)
      g(1,14)=2.0*wo(1,20)-3.0*wo(1,13)-3.0*wo(1,15)
      g(1,15)=4.0*wo(1,16)-wo(1,18)-wo(1,14)
      g(1,16)=4.0*wo(1,17)-wo(1,12)-wo(1,19)
! c
 2009  continue
      if(l1.eq.1) return
      g(2,1)=wo(2,1)
      if(l2.eq.1) goto 3009
      g(2,2)=wo(2,2)
      g(2,3)=wo(2,3)
      g(2,4)=wo(2,4)
      if(l2.eq.4) goto 3009
      g(2,5)=wo(2,8)
      g(2,6)=wo(2,9)
      g(2,7)=wo(2,10)
      g(2,8)=wo(2,5)-wo(2,6)
      g(2,9)=2.0*wo(2,7)-wo(2,5)-wo(2,6)
      if(l2.eq.9) goto 3009
      g(2,10)=wo(2,11)
      g(2,11)=wo(2,13)-wo(2,15)
      g(2,12)=wo(2,18)-3.0*wo(2,14)
      g(2,13)=3.0*wo(2,12)-wo(2,19)
      g(2,14)=2.0*wo(2,20)-3.0*wo(2,13)-3.0*wo(2,15)
      g(2,15)=4.0*wo(2,16)-wo(2,18)-wo(2,14)
      g(2,16)=4.0*wo(2,17)-wo(2,12)-wo(2,19)
! c
 3009  continue
      g(3,1)=wo(3,1)
      if(l2.eq.1) goto 4009
      g(3,2)=wo(3,2)
      g(3,3)=wo(3,3)
      g(3,4)=wo(3,4)
      if(l2.eq.4) goto 4009
      g(3,5)=wo(3,8)
      g(3,6)=wo(3,9)
      g(3,7)=wo(3,10)
      g(3,8)=wo(3,5)-wo(3,6)
      g(3,9)=2.0*wo(3,7)-wo(3,5)-wo(3,6)
      if(l2.eq.9) goto 4009
      g(3,10)=wo(3,11)
      g(3,11)=wo(3,13)-wo(3,15)
      g(3,12)=wo(3,18)-3.0*wo(3,14)
      g(3,13)=3.0*wo(3,12)-wo(3,19)
      g(3,14)=2.0*wo(3,20)-3.0*wo(3,13)-3.0*wo(3,15)
      g(3,15)=4.0*wo(3,16)-wo(3,18)-wo(3,14)
      g(3,16)=4.0*wo(3,17)-wo(3,12)-wo(3,19)
! c
 4009  continue
      g(4,1)=wo(4,1)
      if(l2.eq.1) goto 5009
      g(4,2)=wo(4,2)
      g(4,3)=wo(4,3)
      g(4,4)=wo(4,4)
      if(l2.eq.4) goto 5009
      g(4,5)=wo(4,8)
      g(4,6)=wo(4,9)
      g(4,7)=wo(4,10)
      g(4,8)=wo(4,5)-wo(4,6)
      g(4,9)=2.0*wo(4,7)-wo(4,5)-wo(4,6)
      if(l2.eq.9) goto 5009
      g(4,10)=wo(4,11)
      g(4,11)=wo(4,13)-wo(4,15)
      g(4,12)=wo(4,18)-3.0*wo(4,14)
      g(4,13)=3.0*wo(4,12)-wo(4,19)
      g(4,14)=2.0*wo(4,20)-3.0*wo(4,13)-3.0*wo(4,15)
      g(4,15)=4.0*wo(4,16)-wo(4,18)-wo(4,14)
      g(4,16)=4.0*wo(4,17)-wo(4,12)-wo(4,19)
! c
 5009  continue
      if(l1.eq.4) return
      g(5,1)=wo(8,1)
      if(l2.eq.1) goto 6009
      g(5,2)=wo(8,2)
      g(5,3)=wo(8,3)
      g(5,4)=wo(8,4)
      if(l2.eq.4) goto 6009
      g(5,5)=wo(8,8)
      g(5,6)=wo(8,9)
      g(5,7)=wo(8,10)
      g(5,8)=wo(8,5)-wo(8,6)
      g(5,9)=2.0*wo(8,7)-wo(8,5)-wo(8,6)
      if(l2.eq.9) goto 6009
      g(5,10)=wo(8,11)
      g(5,11)=wo(8,13)-wo(8,15)
      g(5,12)=wo(8,18)-3.0*wo(8,14)
      g(5,13)=3.0*wo(8,12)-wo(8,19)
      g(5,14)=2.0*wo(8,20)-3.0*wo(8,13)-3.0*wo(8,15)
      g(5,15)=4.0*wo(8,16)-wo(8,18)-wo(8,14)
      g(5,16)=4.0*wo(8,17)-wo(8,12)-wo(8,19)
! c
 6009  continue
      g(6,1)=wo(9,1)
      if(l2.eq.1) goto 7009
      g(6,2)=wo(9,2)
      g(6,3)=wo(9,3)
      g(6,4)=wo(9,4)
      if(l2.eq.4) goto 7009
      g(6,5)=wo(9,8)
      g(6,6)=wo(9,9)
      g(6,7)=wo(9,10)
      g(6,8)=wo(9,5)-wo(9,6)
      g(6,9)=2.0*wo(9,7)-wo(9,5)-wo(9,6)
      if(l2.eq.9) goto 7009
      g(6,10)=wo(9,11)
      g(6,11)=wo(9,13)-wo(9,15)
      g(6,12)=wo(9,18)-3.0*wo(9,14)
      g(6,13)=3.0*wo(9,12)-wo(9,19)
      g(6,14)=2.0*wo(9,20)-3.0*wo(9,13)-3.0*wo(9,15)
      g(6,15)=4.0*wo(9,16)-wo(9,18)-wo(9,14)
      g(6,16)=4.0*wo(9,17)-wo(9,12)-wo(9,19)
! c
 7009  continue
      g(7,1)=wo(10,1)
      if(l2.eq.1) goto 8009
      g(7,2)=wo(10,2)
      g(7,3)=wo(10,3)
      g(7,4)=wo(10,4)
      if(l2.eq.4) goto 8009
      g(7,5)=wo(10,8)
      g(7,6)=wo(10,9)
      g(7,7)=wo(10,10)
      g(7,8)=wo(10,5)-wo(10,6)
      g(7,9)=2.0*wo(10,7)-wo(10,5)-wo(10,6)
      if(l2.eq.9) goto 8009
      g(7,10)=wo(10,11)
      g(7,11)=wo(10,13)-wo(10,15)
      g(7,12)=wo(10,18)-3.0*wo(10,14)
      g(7,13)=3.0*wo(10,12)-wo(10,19)
      g(7,14)=2.0*wo(10,20)-3.0*wo(10,13)-3.0*wo(10,15)
      g(7,15)=4.0*wo(10,16)-wo(10,18)-wo(10,14)
      g(7,16)=4.0*wo(10,17)-wo(10,12)-wo(10,19)
! c
 8009  continue
      g(8,1)=wo(5,1)-wo(6,1)
      if(l2.eq.1) goto 9009
      g(8,2)=wo(5,2)-wo(6,2)
      g(8,3)=wo(5,3)-wo(6,3)
      g(8,4)=wo(5,4)-wo(6,4)
      if(l2.eq.4) goto 9009
      g(8,5)=wo(5,8)-wo(6,8)
      g(8,6)=wo(5,9)-wo(6,9)
      g(8,7)=wo(5,10)-wo(6,10)
      g(8,8)=wo(5,5)-wo(6,5)-wo(5,6)+wo(6,6)
      g(8,9)=2.0*wo(5,7)-2.0*wo(6,7)-wo(5,5) &
     & +wo(6,5)-wo(5,6)+wo(6,6)
      if(l2.eq.9) goto 9009
      g(8,10)=wo(5,11)-wo(6,11)
      g(8,11)=wo(5,13)-wo(6,13)-wo(5,15)+wo(6,15)
      g(8,12)=wo(5,18)-wo(6,18)-3.0*wo(5,14) &
     & +3.0*wo(6,14)
      g(8,13)=3.0*wo(5,12)-3.0*wo(6,12)-wo(5,19)+wo(6,19)
      g(8,14)=2.0*wo(5,20)-2.0*wo(6,20)-3.0*wo(5,13) &
     & +3.0*wo(6,13)-3.0*wo(5,15)+3.0*wo(6,15)
      g(8,15)=4.0*wo(5,16)-4.0*wo(6,16)-wo(5,18)+wo(6,18) &
     & -wo(5,14)+wo(6,14)
      g(8,16)=4.0*wo(5,17)-4.0*wo(6,17)-wo(5,12)+wo(6,12) &
     & -wo(5,19)+wo(6,19)
! c
 9009  continue
      g(9,1)=2.0*wo(7,1)-wo(5,1)-wo(6,1)
      if(l2.eq.1) goto 10009
      g(9,2)=2.0*wo(7,2)-wo(5,2)-wo(6,2)
      g(9,3)=2.0*wo(7,3)-wo(5,3)-wo(6,3)
      g(9,4)=2.0*wo(7,4)-wo(5,4)-wo(6,4)
      if(l2.eq.4) goto 10009
      g(9,5)=2.0*wo(7,8)-wo(5,8)-wo(6,8)
      g(9,6)=2.0*wo(7,9)-wo(5,9)-wo(6,9)
      g(9,7)=2.0*wo(7,10)-wo(5,10)-wo(6,10)
      g(9,8)=2.0*wo(7,5)-wo(5,5)-wo(6,5)-2.0*wo(7,6)+wo(5,6)+wo(6,6)
      g(9,9)=4.0*wo(7,7)-2.0*wo(5,7)-2.0*wo(6,7)-2.0*wo(7,5)+wo(5,5) &
     & +wo(6,5)-2.0*wo(7,6)+wo(5,6)+wo(6,6)
      if(l2.eq.9) goto 10009
      g(9,10)=2.0*wo(7,11)-wo(5,11)-wo(6,11)
      g(9,11)=2.0*wo(7,13)-wo(5,13)-wo(6,13)-2.0*wo(7,15) &
     & +wo(5,15)+wo(6,15)
      g(9,12)=2.0*wo(7,18)-wo(5,18)-wo(6,18)-6.0*wo(7,14)+ &
     & 3.0*wo(5,14)+3.0*wo(6,14)
      g(9,13)=6.0*wo(7,12)-3.0*wo(5,12)-3.0*wo(6,12) &
     & -2.0*wo(7,19)+wo(5,19)+wo(6,19)
      g(9,14)=4.0*wo(7,20)-2.0*wo(5,20)-2.0*wo(6,20)-6.0*wo(7,13)+ &
     & 3.0*wo(5,13)+3.0*wo(6,13)-6.0*wo(7,15)+3.0*wo(5,15)+3.0*wo(6,15)
      g(9,15)=8.0*wo(7,16)-4.0*wo(5,16)-4.0*wo(6,16)-2.0*wo(7,18) &
     & +wo(5,18)+wo(6,18)-2.0*wo(7,14)+wo(5,14)+wo(6,14)
      g(9,16)=8.0*wo(7,17)-4.0*wo(5,17)-4.0*wo(6,17)-2.0*wo(7,12) &
     & +wo(5,12)+wo(6,12)-2.0*wo(7,19)+wo(5,19)+wo(6,19)
! c
 10009 continue
      if(l1.eq.9) return
      g(10,1)=wo(11,1)
      if(l2.eq.1) goto 11009
      g(10,2)=wo(11,2)
      g(10,3)=wo(11,3)
      g(10,4)=wo(11,4)
      if(l2.eq.4) goto 11009
      g(10,5)=wo(11,8)
      g(10,6)=wo(11,9)
      g(10,7)=wo(11,10)
      g(10,8)=wo(11,5)-wo(11,6)
      g(10,9)=2.0*wo(11,7)-wo(11,5)-wo(11,6)
      if(l2.eq.9) goto 11009
      g(10,10)=wo(11,11)
      g(10,11)=wo(11,13)-wo(11,15)
      g(10,12)=wo(11,18)-3.0*wo(11,14)
      g(10,13)=3.0*wo(11,12)-wo(11,19)
      g(10,14)=2.0*wo(11,20)-3.0*wo(11,13)-3.0*wo(11,15)
      g(10,15)=4.0*wo(11,16)-wo(11,18)-wo(11,14)
      g(10,16)=4.0*wo(11,17)-wo(11,12)-wo(11,19)
! c
 11009 continue
      g(11,1)=wo(13,1)-wo(15,1)
      if(l2.eq.1) goto 12009
      g(11,2)=wo(13,2)-wo(15,2)
      g(11,3)=wo(13,3)-wo(15,3)
      g(11,4)=wo(13,4)-wo(15,4)
      if(l2.eq.4) goto 12009
      g(11,5)=wo(13,8)-wo(15,8)
      g(11,6)=wo(13,9)-wo(15,9)
      g(11,7)=wo(13,10)-wo(15,10)
      g(11,8)=wo(13,5)-wo(15,5)-wo(13,6)+wo(15,6)
      g(11,9)=2.0*wo(13,7)-2.0*wo(15,7)-wo(13,5) &
     & +wo(15,5)-wo(13,6)+wo(15,6)
      if(l2.eq.9) goto 12009
      g(11,10)=wo(13,11)-wo(15,11)
      g(11,11)=wo(13,13)-wo(15,13)-wo(13,15)+wo(15,15)
      g(11,12)=wo(13,18)-wo(15,18)-3.0*wo(13,14) &
     & +3.0*wo(15,14)
      g(11,13)=3.0*wo(13,12)-3.0*wo(15,12)- &
     & wo(13,19)+wo(15,19)
      g(11,14)=2.0*wo(13,20)-2.0*wo(15,20)-3.0*wo(13,13) &
     & +3.0*wo(15,13)-3.0*wo(13,15)+3.0*wo(15,15)
      g(11,15)=4.0*wo(13,16)-4.0*wo(15,16)-wo(13,18)+wo(15,18) &
     & -wo(13,14)+wo(15,14)
      g(11,16)=4.0*wo(13,17)-4.0*wo(15,17)-wo(13,12)+wo(15,12) &
     & -wo(13,19)+wo(15,19)
! c
 12009 continue
      g(12,1)=wo(18,1)-3.0*wo(14,1)
      if(l2.eq.1) goto 13009
      g(12,2)=wo(18,2)-3.0*wo(14,2)
      g(12,3)=wo(18,3)-3.0*wo(14,3)
      g(12,4)=wo(18,4)-3.0*wo(14,4)
      if(l2.eq.4) goto 13009
      g(12,5)=wo(18,8)-3.0*wo(14,8)
      g(12,6)=wo(18,9)-3.0*wo(14,9)
      g(12,7)=wo(18,10)-3.0*wo(14,10)
      g(12,8)=wo(18,5)-3.0*wo(14,5) &
     & -wo(18,6)+3.0*wo(14,6)
      g(12,9)=2.0*wo(18,7)-6.0*wo(14,7) &
     & -wo(18,5)+3.0*wo(14,5) &
     & -wo(18,6)+3.0*wo(14,6)
      if(l2.eq.9) goto 13009
      g(12,10)=wo(18,11)-3.0*wo(14,11)
      g(12,11)=wo(18,13)-3.0*wo(14,13) &
     & -wo(18,15)+3.0*wo(14,15)
      g(12,12)=wo(18,18)-3.0*wo(14,18) &
     & -3.0*wo(18,14)+9.0*wo(14,14)
      g(12,13)=3.0*wo(18,12)-9.0*wo(14,12) &
     & -wo(18,19)+3.0*wo(14,19)
      g(12,14)=2.0*wo(18,20)-6.0*wo(14,20)-3.0*wo(18,13) &
     & +9.0*wo(14,13)-3.0*wo(18,15)+9.0*wo(14,15)
      g(12,15)=4.0*wo(18,16)-12.0*wo(14,16)-wo(18,18) &
     & +3.0*wo(14,18)-wo(18,14)+3.0*wo(14,14)
      g(12,16)=4.0*wo(18,17)-12.0*wo(14,17)-wo(18,12) &
     & +3.0*wo(14,12)-wo(18,19)+3.0*wo(14,19)
! c
 13009 continue
      g(13,1)=3.0*wo(12,1)-wo(19,1)
      if(l2.eq.1) goto 14009
      g(13,2)=3.0*wo(12,2)-wo(19,2)
      g(13,3)=3.0*wo(12,3)-wo(19,3)
      g(13,4)=3.0*wo(12,4)-wo(19,4)
      if(l2.eq.4) goto 14009
      g(13,5)=3.0*wo(12,8)-wo(19,8)
      g(13,6)=3.0*wo(12,9)-wo(19,9)
      g(13,7)=3.0*wo(12,10)-wo(19,10)
      g(13,8)=3.0*wo(12,5)-wo(19,5)-3.0*wo(12,6)+wo(19,6)
      g(13,9)=6.0*wo(12,7)-2.0*wo(19,7)-3.0*wo(12,5) &
     & +wo(19,5)-3.0*wo(12,6)+wo(19,6)
      if(l2.eq.9) goto 14009
      g(13,10)=3.0*wo(12,11)-wo(19,11)
      g(13,11)=3.0*wo(12,13)-wo(19,13)-3.0*wo(12,15)+wo(19,15)
      g(13,12)=3.0*wo(12,18)-wo(19,18)-9.0*wo(12,14)+3.0*wo(19,14)
      g(13,13)=9.0*wo(12,12)-3.0*wo(19,12)-3.0*wo(12,19)+wo(19,19)
      g(13,14)=6.0*wo(12,20)-2.0*wo(19,20)-9.0*wo(12,13) &
     & +3.0*wo(19,13)-9.0*wo(12,15)+3.0*wo(19,15)
      g(13,15)=12.0*wo(12,16)-4.0*wo(19,16)-3.0*wo(12,18) &
     & +wo(19,18)-3.0*wo(12,14)+wo(19,14)
      g(13,16)=12.0*wo(12,17)-4.0*wo(19,17)-3.0*wo(12,12) &
     & +wo(19,12)-3.0*wo(12,19)+wo(19,19)
! c
 14009 continue
      g(14,1)=2.0*wo(20,1)-3.0*wo(13,1)-3.0*wo(15,1)
      if(l2.eq.1) goto 15009
      g(14,2)=2.0*wo(20,2)-3.0*wo(13,2)-3.0*wo(15,2)
      g(14,3)=2.0*wo(20,3)-3.0*wo(13,3)-3.0*wo(15,3)
      g(14,4)=2.0*wo(20,4)-3.0*wo(13,4)-3.0*wo(15,4)
      if(l2.eq.4) goto 15009
      g(14,5)=2.0*wo(20,8)-3.0*wo(13,8)-3.0*wo(15,8)
      g(14,6)=2.0*wo(20,9)-3.0*wo(13,9)-3.0*wo(15,9)
      g(14,7)=2.0*wo(20,10)-3.0*wo(13,10)-3.0*wo(15,10)
      g(14,8)=2.0*wo(20,5)-3.0*wo(13,5)-3.0*wo(15,5) &
     & -2.0*wo(20,6)+3.0*wo(13,6)+3.0*wo(15,6)
      g(14,9)=4.0*wo(20,7)-6.0*wo(13,7)-6.0*wo(15,7) &
     & -2.0*wo(20,5)+3.0*wo(13,5)+3.0*wo(15,5) &
     & -2.0*wo(20,6)+3.0*wo(13,6)+3.0*wo(15,6)
      if(l2.eq.9) goto 15009
      g(14,10)=2.0*wo(20,11)-3.0*wo(13,11)-3.0*wo(15,11)
      g(14,11)=2.0*wo(20,13)-3.0*wo(13,13)-3.0*wo(15,13) &
     & -2.0*wo(20,15)+3.0*wo(13,15)+3.0*wo(15,15)
      g(14,12)=2.0*wo(20,18)-3.0*wo(13,18)-3.0*wo(15,18) &
     & -6.0*wo(20,14)+9.0*wo(13,14)+9.0*wo(15,14)
      g(14,13)=6.0*wo(20,12)-9.0*wo(13,12)-9.0*wo(15,12) &
     & -2.0*wo(20,19)+3.0*wo(13,19)+3.0*wo(15,19)
       g(14,14)=4.0*wo(20,20)-6.0*wo(13,20)-6.0*wo(15,20) &
     & -6.0*wo(20,13)+9.0*wo(13,13)+9.0*wo(15,13) &
     & -6.0*wo(20,15)+9.0*wo(13,15)+9.0*wo(15,15)
      g(14,15)=8.0*wo(20,16)-12.0*wo(13,16)-12.0*wo(15,16) &
     & -2.0*wo(20,18)+3.0*wo(13,18)+3.0*wo(15,18) &
     & -2.0*wo(20,14)+3.0*wo(13,14)+3.0*wo(15,14)
      g(14,16)=8.0*wo(20,17)-12.0*wo(13,17)-12.0*wo(15,17) &
     & -2.0*wo(20,12)+3.0*wo(13,12)+3.0*wo(15,12) &
     & -2.0*wo(20,19)+3.0*wo(13,19)+3.0*wo(15,19)
! c
 15009 continue
      g(15,1)=4.0*wo(16,1)-wo(18,1)-wo(14,1)
      if(l2.eq.1) goto 16009
      g(15,2)=4.0*wo(16,2)-wo(18,2)-wo(14,2)
      g(15,3)=4.0*wo(16,3)-wo(18,3)-wo(14,3)
      g(15,4)=4.0*wo(16,4)-wo(18,4)-wo(14,4)
      if(l2.eq.4) goto 16009
      g(15,5)=4.0*wo(16,8)-wo(18,8)-wo(14,8)
      g(15,6)=4.0*wo(16,9)-wo(18,9)-wo(14,9)
      g(15,7)=4.0*wo(16,10)-wo(18,10)-wo(14,10)
      g(15,8)=4.0*wo(16,5)-wo(18,5)-wo(14,5) &
     & -4.0*wo(16,6)+wo(18,6)+wo(14,6)
      g(15,9)=8.0*wo(16,7)-2.0*wo(18,7)-2.0*wo(14,7) &
     & -4.0*wo(16,5)+wo(18,5)+wo(14,5) &
     & -4.0*wo(16,6)+wo(18,6)+wo(14,6)
      if(l2.eq.9) goto 16009
      g(15,10)=4.0*wo(16,11)-wo(18,11)-wo(14,11)
      g(15,11)=4.0*wo(16,13)-wo(18,13)-wo(14,13) &
     & -4.0*wo(16,15)+wo(18,15)+wo(14,15)
      g(15,12)=4.0*wo(16,18)-wo(18,18)-wo(14,18) &
     & -12.0*wo(16,14)+3.0*wo(18,14)+3.0*wo(14,14)
      g(15,13)=12.0*wo(16,12)-3.0*wo(18,12)-3.0*wo(14,12) &
     & -4.0*wo(16,19)+wo(18,19)+wo(14,19)
      g(15,14)=8.0*wo(16,20)-2.0*wo(18,20)-2.0*wo(14,20) &
     & -12.0*wo(16,13)+3.0*wo(18,13)+3.0*wo(14,13) &
     & -12.0*wo(16,15)+3.0*wo(18,15)+3.0*wo(14,15)
      g(15,15)=16.0*wo(16,16)-4.0*wo(18,16)-4.0*wo(14,16) &
     & -4.0*wo(16,18)+wo(18,18)+wo(14,18) &
     & -4.0*wo(16,14)+wo(18,14)+wo(14,14)
      g(15,16)=16.0*wo(16,17)-4.0*wo(18,17)-4.0*wo(14,17) &
     & -4.0*wo(16,12)+wo(18,12)+wo(14,12) &
     & -4.0*wo(16,19)+wo(18,19)+wo(14,19)
! c
 16009 continue
      g(16,1)=4.0*wo(17,1)-wo(12,1)-wo(19,1)
      if(l2.eq.1) goto 17009
      g(16,2)=4.0*wo(17,2)-wo(12,2)-wo(19,2)
      g(16,3)=4.0*wo(17,3)-wo(12,3)-wo(19,3)
      g(16,4)=4.0*wo(17,4)-wo(12,4)-wo(19,4)
      if(l2.eq.4) goto 17009
      g(16,5)=4.0*wo(17,8)-wo(12,8)-wo(19,8)
      g(16,6)=4.0*wo(17,9)-wo(12,9)-wo(19,9)
      g(16,7)=4.0*wo(17,10)-wo(12,10)-wo(19,10)
      g(16,8)=4.0*wo(17,5)-wo(12,5)-wo(19,5) &
     & -4.0*wo(17,6)+wo(12,6)+wo(19,6)
      g(16,9)=8.0*wo(17,7)-2.0*wo(12,7)-2.0*wo(19,7) &
     & -4.0*wo(17,5)+wo(12,5)+wo(19,5) &
     & -4.0*wo(17,6)+wo(12,6)+wo(19,6)
      if(l2.eq.9) goto 17009
      g(16,10)=4.0*wo(17,11)-wo(12,11)-wo(19,11)
      g(16,11)=4.0*wo(17,13)-wo(12,13)-wo(19,13) &
     & -4.0*wo(17,15)+wo(12,15)+wo(19,15)
      g(16,12)=4.0*wo(17,18)-wo(12,18)-wo(19,18) &
     & -12.0*wo(17,14)+3.0*wo(12,14)+3.0*wo(19,14)
      g(16,13)=12.0*wo(17,12)-3.0*wo(12,12)-3.0*wo(19,12) &
     & -4.0*wo(17,19)+wo(12,19)+wo(19,19)
      g(16,14)=8.0*wo(17,20)-2.0*wo(12,20)-2.0*wo(19,20) &
     & -12.0*wo(17,13)+3.0*wo(12,13)+3.0*wo(19,13) &
     & -12.0*wo(17,15)+3.0*wo(12,15)+3.0*wo(19,15)
      g(16,15)=16.0*wo(17,16)-4.0*wo(12,16)-4.0*wo(19,16) &
     & -4.0*wo(17,18)+wo(12,18)+wo(19,18) &
     & -4.0*wo(17,14)+wo(12,14)+wo(19,14)
      g(16,16)=16.0*wo(17,17)-4.0*wo(12,17)-4.0*wo(19,17) &
     & -4.0*wo(17,12)+wo(12,12)+wo(19,12) &
     & -4.0*wo(17,19)+wo(12,19)+wo(19,19)
 17009 continue
      return      




      end subroutine overlapInteg

      subroutine momf(q,nx1,nx2,a1,a2,a,b)

      use O_Kinds
      use O_Constants

! c
      implicit none
! c
      real (kind=double) :: a1, a2
      real (kind=double), dimension (16,16,3) :: q
      real (kind=double), dimension (3,3) :: det
      real (kind=double), dimension (3)   :: a,b,r
      real (kind=double), dimension (7,3) :: fs,sf
      real (kind=double), dimension (7,3,3) :: fp
      real (kind=double), dimension (3,7,3) :: pf
      real (kind=double), dimension (7,5,3) :: fd
      real (kind=double), dimension (5,7,3) :: df
      real (kind=double), dimension (7,7,3) :: ff


      integer :: nx1,nx2,i,j,k,l,m,i1,ia,ib,ik,ic,il,im
      real (kind=double) :: rab2, u, alamda, expa, zeta, w, www, delta
      real (kind=double) :: twoal, b1, u2, u3, psum, sn, delt
      real (kind=double) :: c1, c2, c3, c4, c5, cdc, cdc1, cdc2, ci
      real (kind=double) :: abx, aby, abz, a12, a22, a23, a13
      real (kind=double) :: alp, alp2, alp4, alp6, alp8, alp10, alp12
      real (kind=double) :: alpx, alpy, alpz
      real (kind=double) :: ro, rox, roy, roz, rox2, roy2, roz2
      real (kind=double) :: roxy, roxz, royz, rox3, roy3, roz3
      real (kind=double) :: rox2y, rox2z, roxy2, roxz2, roy2z, royz2
      real (kind=double) :: roxyz, rox4, roy4, roz4
      real (kind=double) :: rox3y, rox3z, roxy3, roxz3, roy3z, royz3
      real (kind=double) :: rox2y2, rox2z2, roy2z2,rox2yz,roxy2z,roxyz2
      real (kind=double) :: rox5, roy5, roz5, rox4y, rox4z, roxy4, roxz4
      real (kind=double) :: royz4, rox3y2, rox3z2, rox2y3, rox2z3
      real (kind=double) :: roy3z2, roy2z3, rox3yz, roxy3z, roxyz3
      real (kind=double) :: rox2y2z, rox2yz2, roxy2z2
      real (kind=double) :: rox6, roy6, roz6
      real (kind=double) :: rox5y, rox5z, roxy5, roxz5, roy5z, royz5
      real (kind=double) :: rox4z2
      real (kind=double) :: roy4z, rox4y2, rox2y4, roxy4z, rox2z4, roy4z2
      real (kind=double) :: roy2z4, rox4yz, roxyz4, rox3y3, rox3z3
      real (kind=double) :: roy3z3, rox3y2z, rox3yz2, rox2y3z, rox2yz3
      real (kind=double) :: roxy3z2, roxy2z3, rox2y2z2
      real (kind=double) :: rox7, roy7, roz7
      real (kind=double) :: rox6y, rox6z, roxy6, roxz6, roy6z, royz6
      real (kind=double) :: rox5y2, rox5z2, rox2y5,rox2z5,roy5z2,roy2z5
      real (kind=double) :: roxyz5, rox5yz, roxy5z, rox4y3, rox4z3
      real (kind=double) :: rox3y4, rox3z4, roy4z3, roy3z4, rox4y2z
      real (kind=double) :: rox4yz2, rox2y4z, rox2yz4, roxy2z4, roxy4z2
      real (kind=double) :: rox3y3z, roxy3z3, rox3yz3, rox3y2z2
      real (kind=double) :: rox2y3z2, rox2y2z3
      real (kind=double) :: s1, s2, s3
      real (kind=double) :: x3xs, y3ys, z3zs, x3ys, x3zs, y3xs, y3zs
      real (kind=double) :: z3xs, z3ys
      real (kind=double) :: x2yxs, y2xys, z2xzs, y2zys, z2yzs, x2zxs
      real (kind=double) :: x2yys, x2zzs, y2xxs, z2yys, x2yzs, x2zys
      real (kind=double) :: y2xzs, y2zxs, z2xys, z2yxs, xyzxs, xyzys
      real (kind=double) :: xyzzs
      real (kind=double) :: p1, p2, p3
      real (kind=double) :: x3xx, y3yy, z3zz, x3yy, x3zz, y3xx
      real (kind=double) :: y2zzs, z2xxs, y3zz, z3xx, z3yy, x3xy
      real (kind=double) :: x3xz, y3yx, y3yz, z3zx, z3zy, x3yx, x3zx
      real (kind=double) :: y3xy, y3zy, z3xz, z3yz, x3yz, x3zy
      real (kind=double) :: y3xz, y3zx, z3xy, z3yx, x2yxx, x2zxx
      real (kind=double) :: y2xyy, y2zyy, z2xzz, z2yzz, x2yyy, x2zzz
      real (kind=double) :: z2xxx, z2yyy, x2zxz, y2zyz, z2xzx, z2yzy
      real (kind=double) :: x2yyx, x2zzx, y2xxy, y2zzy, z2xxz, z2yyz
      real (kind=double) :: x2yxz, x2zxy, y2xyz, y2xxx, y2zzz, x2yxy
      real (kind=double) :: y2xyx, y2zyx, z2xzy, z2yzx, x2yyz, x2zzy
      real (kind=double) :: y2xxz, y2zzx, z2xxy, x2yzx, x2zyx, y2xzy
      real (kind=double) :: y2zxy, z2xyz, z2yxz, x2yzy, x2zyz, y2zxz
      real (kind=double) :: z2xyx, z2yxy, z2yyx, y2xzx, x2yzz, x2zyy
      real (kind=double) :: y2xzz, y2zxx, z2xyy, z2yxx, xyzxx, xyzyy
      real (kind=double) :: xyzzz, xyzxy, xyzxz, xyzyx, xyzyz, xyzzx
      real (kind=double) :: xyzzy
      real (kind=double) :: d1, d2, d3, d4, d6, d7
      real (kind=double) :: x3xx2, y3yy2, z3zz2, x3xy2, x3xz2, y3yx2
      real (kind=double) :: y3yz2, z3zx2, z3zy2, x3xxy, x3xxz, y3yxy
      real (kind=double) :: y3yyz, z3zxz, z3zyz, x3yx2, x3zx2, y3xy2
      real (kind=double) :: y3zy2, z3xz2, z3yz2, x3yy2, x3zz2, y3xx2
      real (kind=double) :: y3zz2, z3xx2, z3yy2, x3yz2, y3zx2, z3xy2
      real (kind=double) :: x3zy2, y3xz2, z3yx2, x3yxy, x3zxz, y3xxy
      real (kind=double) :: y3zyz, z3xxz, z3yyz, x3zxy, x3yxz, y3xyz
      real (kind=double) :: y3zxy, z3xyz, z3yxz, x3xyz, y3yxz, z3zxy
      real (kind=double) :: x2yxx2,y2xyy2,z2xzz2,x2zxx2,y2zyy2,z2yzz2
      real (kind=double) :: x2yxy2,x2zxz2,y2xyx2,y2zyz2,z2xzx2,z2yzy2
      real (kind=double) :: x2yxz2,x2zxy2,y2xyz2,y2zyx2,z2xzy2,z2yzx2
      real (kind=double) :: x2yxxy,y2zyyz,z2xzxz,x2zxxz,y2xyxy,z2yzyz
      real (kind=double) :: x2yxxz,x2zxxy,y2xyyz,y2zyxy,z2xzyz,z2yzxz
      real (kind=double) :: x2yxyz,x2zxyz,y2xyxz,y2zyxz,z2xzxy,z2yzxy
      real (kind=double) :: x2yyx2,x2zzx2,y2xxy2,y2zzy2,z2xxz2,z2yyz2
      real (kind=double) :: x2yyy2,y2xxx2,y2zzz2,z2xxx2,z2yyy2,x2yyz2
      real (kind=double) :: x2zzy2,y2xxz2,y2zzx2,z2xxy2,z2yyx2,x2zzz2
      real (kind=double) :: x2yyxy,x2zzxz,y2xxxy,y2zzyz,z2xxxz,z2yyyz
      real (kind=double) :: x2yyxz,x2zzxy,y2xxyz,y2zzxy,z2xxyz,z2yyxz
      real (kind=double) :: x2yyyz,x2zzyz,y2xxxz,y2zzxz,z2xxxy,z2yyxy
      real (kind=double) :: x2yzx2,x2zyx2,y2xzy2,y2zxy2,z2xyz2,z2yxz2
      real (kind=double) :: x2yzy2,x2zyz2,y2xzx2,y2zxz2,z2xyx2,z2yxy2
      real (kind=double) :: x2yzz2,x2zyy2,y2xzz2,y2zxx2,z2xyy2,z2yxx2
      real (kind=double) :: x2yzxy,x2zyxz,y2xzyx,y2zxyz,z2xyxz,z2yxyz
      real (kind=double) :: x2yzxz,x2zyxy,y2xzyz,y2zxxy,z2xyyz,z2yxxz
      real (kind=double) :: x2yzyz,x2zyyz,y2xzxz,y2zxxz,z2xyxy
      real (kind=double) :: z2yxxy,xyzxx2,xyzyy2,xyzxy2,xyzxz2,xyzyx2
      real (kind=double) :: xyzyz2,xyzzx2,xyzzy2,xyzxxy,xyzxxz,xyzyxy
      real (kind=double) :: xyzyyz,xyzzxz,xyxzyz,xyzxyz,xyzyxz,xyzzxy
      real (kind=double) :: x3yyz,y2xzxy,xyzzz2,xyzzyz,y3xxz,z3xxy
      real (kind=double) :: x3zyz,y3zxz,z3yxy
      real (kind=double) :: f1, f2, f3, f4, f5, f6, f7, f8
      real (kind=double) :: x3xx3,y3yy3,z3zz3,x3xy3,x3xz3,y3yx3,y3yz3
      real (kind=double) :: z3zx3,z3zy3,x3xx2y,x3xx2z,y3yy2x,y3yy2z
      real (kind=double) :: z3zz2x,z3zz2y,x3xy2x,x3xz2x,y3yx2y,y3yz2y
      real (kind=double) :: z3zx2z,z3zy2z,x3xy2z,x3xz2y,y3yz2x,y3yx2z
      real (kind=double) :: z3zy2x,z3zx2y,x3xxyz,y3yxyz,z3zxyz,x3yx3
      real (kind=double) :: x3zx3,y3xy3,y3zy3,z3xz3,z3yz3,x3yy3,x3zz3
      real (kind=double) :: y3xx3,y3zz3,z3xx3,z3yy3,x3yz3,x3zy3,y3xz3
      real (kind=double) :: y3zx3,x3yx2y,x3zx2z,y3xy2x,y3zy2z,z3xz2x
      real (kind=double) :: z3yz2y,x3yx2z,x3zx2y,y3xy2z,y3zy2x,z3xz2y
      real (kind=double) :: z3xy3,z3yx3,z3yz2x,x3yy2x,x3zz2x,y3xx2y
      real (kind=double) :: y3zz2y,z3xx2z,z3yy2z,x3yy2z,x3zz2y,y3xx2z
      real (kind=double) :: y3zz2x,z3xx2y,z3yy2x,x3yz2x,x3zy2x,y3xz2y
      real (kind=double) :: y3zx2y,z3xy2z,z3yx2z,x3yz2y,x3zzy2z,y3xz2x
      real (kind=double) :: y3zx2z,x3zy2z,z3xy2x,z3yx2y,x3yxyz,x3zxyz
      real (kind=double) :: y3xxyz,y3zxyz,z3xxyz,z3yxyz,x2yxx3,x2zxx3
      real (kind=double) :: y2xyy3,y2zyy3,z2xzz3,z2yzz3,x2yxy3,x2zxz3
      real (kind=double) :: y2xyx3,y2zyz3,z2xzx3,z2yzy3,x2yxz3,x2zxy3
      real (kind=double) :: y2xyz3,y2zyx3,z2xzy3,z2yzx3,x2yxx2y,x2zxx2z
      real (kind=double) :: y2xyy2x,y2zyy2z,z2xzz2x,z2yzz2y,x2yxx2z
      real (kind=double) :: x2zxx2y,y2xyy2z,y2zyy2x,z2xzz2y,z2yzz2x
      real (kind=double) :: x2yxy2x,x2zxz2x,y2xyx2y,y2zyz2y,z2xzx2z
      real (kind=double) :: z2yzy2z,x2yxy2z,x2zxz2y,y2xyx2z,y2zyz2x
      real (kind=double) :: z2xzx2y,z2yzy2x,x2yxz2x,x2zxy2x,y2xyz2y
      real (kind=double) :: y2zyx2y,z2xzy2z,z2yzx2z,x2yxz2y,x2zxy2z
      real (kind=double) :: y2xyz2x,y2zyx2z,z2xzy2x,z2yzx2y,x2yxxyz
      real (kind=double) :: x2zxxyz,y2xyxyz,y2zyxyz,z2xzxyz,z2yzxyz
      real (kind=double) :: x2yyx3,x2zzx3,y2xxy3,y2zzy3,z2xxz3,z2yyz3
      real (kind=double) :: x2yyy3,x2zzz3,y2xxx3,y2zzz3,z2xxx3,z2yyy3
      real (kind=double) :: x2yyz3,x2zzy3,y2xxz3,y2zzx3,z2xxy3,z2yyx3
      real (kind=double) :: x2yyx2y,x2zzx2z,y2xxy2x,y2zzy2z,z2xxz2x
      real (kind=double) :: z2yyz2y,x2yyx2z,x2zzx2y,y2xxy2z,y2zzy2x
      real (kind=double) :: z2xxz2y,z2yyz2x,x2yyy2x,x2zzz2x,y2xxx2y
      real (kind=double) :: y2zzz2y,z2xxx2z,z2yyy2z,x2yyy2z,x2zzz2y
      real (kind=double) :: y2xxx2z,y2zzz2x,z2xxx2y,z2yyy2x,x2yyz2x
      real (kind=double) :: x2zzy2x,y2xxz2y,y2zzx2y,z2xxy2z,z2yyx2z
      real (kind=double) :: x2yyz2y,x2zzy2z,y2xxz2x,y2zzx2z,x2xxy2x
      real (kind=double) :: z2yyx2y,x2yyxyz,x2zzxyz,z2xxy2x,y2xxxyz
      real (kind=double) :: y2zzxyz,z2xxxyz,z2yyxyz,x2yzx3,x2zyx3
      real (kind=double) :: y2xzy3,y2zxy3,z2xyz3,z2yxz3,x2yzy3,x2zyz3
      real (kind=double) :: y2zxz3,y2xzx3,z2xyx3,z2yxy3,x2yzz3,x2zyy3
      real (kind=double) :: y2xzz3,y2zxx3,z2xyy3,z2yxx3,x2yzx2y,x2zyx2z
      real (kind=double) :: y2xzy2x,y2zxy2z,z2xyz2x,x2yzx2z
      real (kind=double) :: x2zyx2y,y2xzy2z,y2zxy2x,z2xyz2y,z2yxz2x
      real (kind=double) :: x2yzy2x,x2zyz2x,y2xzx2y,y2zxz2y,z2xyx2z
      real (kind=double) :: z2yxy2z,x2yzy2z,x2zyz2y,y2xzx2z,y2zxz2x
      real (kind=double) :: z2xyx2y,z2zxz2x,x2yzz2x,x2zyy2x,y2xzz2y
      real (kind=double) :: z2yxz2y,z2yxy2x,y2zxx2y,z2xyy2z,z2yxx2z
      real (kind=double) :: x2yzz2y,x2zyy2z,y2xzz2x,y2zxx2z,z2xyy2x
      real (kind=double) :: z2yxx2y,x2yzxyz,x2zyxyz,y2xzxyz,y2zxxyz
      real (kind=double) :: z2xyxyz,z2yxxyz,xyzxx3,xyzyy3,xyzzy3
      real (kind=double) :: xyzzz3,xyzxy3,xyzyx3,xyzzx3,xyzxz3,xyzyz3
      real (kind=double) :: xyzxx2y,xyzxx2z,xyzyy2x,xyzyy2z,xyzzz2y
      real (kind=double) :: xyzxy2x,xyzxz2x,xyzyz2y,xyzyx2y,xyzzy2z
      real (kind=double) :: xyzzx2z,xyzxy2z,xyzxz2y,xyzyz2x,xyzyx2z
      real (kind=double) :: xyzzx2y,xyzzy2x,xyzxxyz,xyzyxyz,xyzzxyz
      real (kind=double) :: xyzz2x,xyzzz2x




! c
! c     for q:
! c     s, p1=x, p2=y, p3=z, 
! c     d1=xy, d2=xz, d3=yz, d4=(x*x-y*y)/2, d5=(2*z*z-x*x-y*y)/sqrt(12),
! c     f1=sqrt(60)*x*y*z, f2=sqrt(15)*(x*x*z-y*y*z),
! c     f3=sqrt(2.5)*(x*x*x-3*y*y*x), f4=sqrt(2.5)*(3*x*x*y-y*y*y),
! c     f5=(2*z*z*z-3*x*x*z-3*y*y*z), f6=sqrt(1.5)*(4*z*z*x-x*x*x-y*y*x),
! c     f7=sqrt(1.5)*(4*z*z*y-x*x*y-y*y*y)
! c
      do 10 i=1,3
      det(i,i)=1.d0
      do 10 j=1,3
      if(i.eq.j) go to 10
      det(i,j)=0.d0
   10 continue
      do 20 i=1,3
      r(i)=b(i)-a(i)
   20 continue
      rab2=r(1)*r(1)+r(2)*r(2)+r(3)*r(3)
! c
      do 30 i=1,16
      do 30 j=1,16
      do 30 k=1,3
   30 q(i,j,k)=0.d0
! c
      u=1.d0/(a1+a2)
      alamda=a1*a2*u
      expa=alamda*rab2
!      if(expa.gt.emax) return
      zeta=dexp(-expa)
      w=pi*u
      www=w*w*w
      delta=dsqrt(www)
      twoal=2.d0*alamda
      b1=delta*zeta
      c1=2.d0*alamda*b1
      u2=a2*u
      u3=a1*u
      c2=u2*b1
      c3=u3*b1
      c4=u*alamda*b1
      c5=(alamda**3.d0)*b1/((a1*a2)**2.d0)
      cdc=c1
      do 40 k=1,3
!********************************************************
!*     s:del:s
!********************************************************
      q(1,1,k)=r(k)*cdc
   40 continue
! c
      if(nx2.ge.4) then
      cdc2=c3
      do 50 l=1,3
      do 50 k=1,3
!***********************************************************
!*     p:del:s
!***********************************************************
      q(1,l+1,k)=-(twoal*r(k)*r(l)-det(k,l))*cdc2
   50 continue
      endif
! c
      if(nx1.ge.4) then
      cdc1=c2
      do 60 l=1,3
      do 60 k=1,3
!***********************************************************
!*     p:del:s
!***********************************************************
      q(l+1,1,k)=(twoal*r(k)*r(l)-det(k,l))*cdc1
   60 continue
      endif
! c
      if(nx1.ge.4.and.nx2.ge.4) then
      cdc=c4
      do 70 l=1,3
      do 70 m=1,3
      do 70 k=1,3
!*************************************************************
!*     p:del:p
!*************************************************************
      q(l+1,m+1,k)=-(twoal*r(k)*r(l)*r(m)-det(k,l)*r(m)- &
     & det(l,m)*r(k)-det(k,m)*r(l))*cdc
   70 continue
      endif
! c
      if(nx2.gt.4) then
      cdc2=c3/a2
      do 120 l=1,5
      do 120 k=1,3
      if(l.gt.3)go to 90
      if(l.eq.k) then
      ia=k+1
      ib=l
      if(k.eq.3) ia=1
      go to 80
      endif
      ik=k-1
      if(k.eq.1) ik=3
      if(l.eq.ik) then
      ia=l
      ib=k
      go to 80
      endif
      ik=k+1
      if(k.eq.3) ik=1
      if(l.eq.ik) then
      psum=2.d0*(alamda**2)*r(1)*r(2)*r(3)
!***********************************************************
!*     d(ij):del(k):s
!***********************************************************
      q(1,l+4,k)=psum*cdc2
      go to 120
      endif
   80 psum=-alamda*r(ia)*(1.d0-2.d0*alamda*r(ib)**2)
!***********************************************************
!*     d(ij):del(i):s     or
!*     d(ij):del(j):s
!***********************************************************
      q(1,l+4,k)=psum*cdc2
   90 if(l.ne.4) go to 110
      if(k.eq.1) then
      ia=1
      ci=-1.d0
      go to 100
      endif
      if(k.eq.2) then
      ia=2
      ci=-1.d0
      go to 100
      endif
      if(k.eq.3) then
      ia=3
      ci=0.
      endif
  100 psum=2.d0*alamda*r(ia)*(ci+alamda*(r(1)**2-r(2)**2))
!************************************************************
!*     d(x2-y2):del:s
!************************************************************
      q(1,l+4,k)=psum*cdc2
      go to 120
  110 if(l.ne.5) go to 120
      if(k.ne.3) ci=-1.
      if(k.eq.3) ci=2.
      psum=-2.d0*alamda*r(k)*(alamda*(r(1)**2+r(2)**2-2.d0*r(3)**2) &
     & +ci)
!************************************************************
!*     d(3z2-r2):del:s
!************************************************************
      q(1,l+4,k)=psum*cdc2
  120 continue
      endif
! c
      if(nx1.gt.4) then
      cdc1=c2/a1
      do 170 l=1,5
      do 170 k=1,3
      if(l.gt.3)go to 140
      if(l.eq.k) then
      ia=k+1
      ib=l
      if(k.eq.3) ia=1
      go to 130
      endif
      ik=k-1
      if(k.eq.1) ik=3
      if(l.eq.ik) then
      ia=l
      ib=k
      go to 130
      endif
      ik=k+1
      if(k.eq.3) ik=1
      if(l.eq.ik) then
      psum=2.d0*(alamda**2)*r(1)*r(2)*r(3)
!***********************************************************
!*     d(ij):del(k):s
!***********************************************************
      q(l+4,1,k)=psum*cdc1
      go to 170
      endif
  130 psum=-alamda*r(ia)*(1.d0-2.d0*alamda*r(ib)**2)
!***********************************************************
!*     d(ij):del(i):s     or
!*     d(ij):del(j):s
!***********************************************************
      q(l+4,1,k)=psum*cdc1
  140 if(l.ne.4) go to 160
      if(k.eq.1) then
      ia=1
      ci=-1.d0
      go to 150
      endif
      if(k.eq.2) then
      ia=2
      ci=-1.d0
      go to 150
      endif
      if(k.eq.3) then
      ia=3
      ci=0.
      endif
  150 psum=2.d0*alamda*r(ia)*(ci+alamda*(r(1)**2-r(2)**2))
!************************************************************
!*     d(x2-y2):del:s
!************************************************************
      q(l+4,1,k)=psum*cdc1
      go to 170
  160 if(l.ne.5) go to 170
      if(k.ne.3) ci=-1.
      if(k.eq.3) ci=2.
      psum=-2.d0*alamda*r(k)*(alamda*(r(1)**2+r(2)**2-2.d0*r(3)**2) &
     & +ci)
!************************************************************
!*     d(3z2-r2):del:s
!************************************************************
      q(l+4,1,k)=psum*cdc1
  170 continue
      endif
! c
      if(nx1.ge.4.and.nx2.gt.4) then
      cdc2=c3/(a1*a2)
      do 350 l=1,5
      do 350 k=1,3
      do 350 m=1,3
      if(l.eq.4) go to 230
      if(l.eq.5) go to 300
      if(l.eq.k.and.k.eq.m)then
      ia=l
      ib=l+1
      ic=k
      if(l.eq.3) ib=1
      go to 180
      endif
      il=l+1
      if(l.eq.3) il=1
      if(il.eq.k.and.k.eq.m)then
      ia=l
      ib=il
      ic=k
      go to 180
      endif
      go to 190
  180 psum=(alamda**2.d0)*r(ia)*r(ib)*(3.d0-2.d0*alamda*r(ic)**2.d0)
!************************************************************
!*     d(ij):del(i):p(i)
!*     d(ij):del(j):p(j)
!************************************************************
      q(m+1,l+4,k)=-psum*cdc2
      go to 350
  190 continue
      im=m-1
      if(m.eq.1) im=3
      if(l.eq.k.and.k.eq.im) then
      ia=k
      ib=m
      go to 200
      endif
      ik=k-1
      if(k.eq.1) ik=3
      if(l.eq.ik.and.l.eq.m) then
      ia=m
      ib=k
      go to 200
      endif
      go to 210
  200 psum=-alamda*(1.d0-2.d0*alamda*r(ia)**2.d0)*(1.d0-2.d0*alamda &
     & *r(ib)**2.d0)/2.d0
!************************************************************
!*     d(ij):del(i):p(j)
!*     d(ij):del(j):p(i)
!************************************************************
      q(m+1,l+4,k)=-psum*cdc2
      go to 350
  210 continue
      im=m+1
      if(m.eq.3) im=1
      if(l.eq.k.and.l.eq.im) then
      ia=m
      ib=k+1
      ic=l
      if(k.eq.3) ib=1
      go to 220
      endif
      ik=k-1
      if(k.eq.1) ik=3
      if(l.eq.ik.and.l.eq.im) then
      ia=m
      ib=l
      ic=k
      go to 220
      endif
      ik=k+1
      if(k.eq.3) ik=1
      if(l.eq.ik.and.l.eq.m) then
      ia=k
      ib=m+1
      if(m.eq.3) ib=1
      ic=l
      go to 220
      endif
      im=m-1
      if(m.eq.1) im=3
      if(l.eq.ik.and.l.eq.im) then
      ia=k
      ib=l
      ic=m
      go to 220
      endif
      ik=k+1
      if(k.eq.3) ik=1
      if(l.eq.ik.and.k.eq.m) then
      ia=l
      ib=l+1
      if(l.eq.3) ib=1
      ic=k
      go to 220
      endif
      go to 230
  220 psum=(alamda**2.d0)*r(ia)*r(ib)*(1.d0-2.d0*alamda*r(ic)**2.d0)
!************************************************************
!*     d(ij):del(i):p(k)     or
!*     d(ij):del(j):p(k)     or
!*     d(ij):del(k):p(i)     or
!*     d(ij):del(k):p(j)     or
!*     d(ij):del(k):p(k)
!************************************************************
      q(m+1,l+4,k)=-psum*cdc2
      go to 350
  230 continue
      if(l.ne.4) go to 300
      if(k.eq.1.and.m.eq.1) then
      ia=1
      ib=2
      sn=1.d0
      go to 240
      endif
      if(k.eq.2.and.m.eq.2) then
      ia=2
      ib=1
      sn=-1.d0
      go to 240
      endif
      go to 250
  240 psum=sn*alamda*(2.d0*alamda*r(ia)**2.d0-(1.d0-alamda*(r(ia) &
     & **2.d0-r(ib)**2.d0))*(1.d0-2.d0*alamda*r(ia)**2.d0))
!************************************************************
!*     d(x2-y2):del(x):p(x)     or
!*     d(x2-y2):del(y):p(y)
!************************************************************
      q(m+1,l+4,k)=-psum*cdc2
      go to 350
  250 continue
      if(k.eq.1.and.m.eq.2) go to 260
      if(k.eq.2.and.m.eq.1) go to 260
      go to 270
  260 psum=-2.d0*(alamda**3.d0)*r(1)*r(2)*(r(1)**2.d0-r(2)**2.d0)
!************************************************************
!*     d(x2-y2):del(x):p(y)
!*     d(x2-y2):del(y):p(x)
!************************************************************
      q(m+1,l+4,k)=-psum*cdc2
      go to 350
  270 continue
      if(k.eq.1.and.m.eq.3) then
      ia=1
      ci=-1.d0
      go to 280
      endif
      if(k.eq.2.and.m.eq.3) then
      ia=2
      ci=1.d0
      go to 280
      endif
      if(k.eq.3.and.m.eq.1) then
      ia=1
      ci=-1.d0
      go to 280
      endif
      if(k.eq.3.and.m.eq.2) then
      ia=2
      ci=1.d0
      go to 280
      endif
      go to 290
  280 psum=-2.d0*(alamda**2.d0)*r(ia)*r(3)*(ci &
     & +alamda*(r(1)**2.d0-r(2)**2.d0))
!************************************************************
!*     d(x2-y2):del(x):p(z)     or
!*     d(x2-y2):del(y):p(z)     or
!*     d(x2-y2):del(z):p(x)     or
!*     d(x2-y2):del(z):p(y)
!************************************************************
      q(m+1,l+4,k)=-psum*cdc2
      go to 350
  290 continue
      if(k.eq.3.and.m.eq.3) then
      psum=(alamda**2.d0)*(r(1)**2.d0-r(2)**2.d0)*(1.d0-2.d0*alamda* &
     & r(3)**2.d0)
!************************************************************
!*     d(x2-y2):del(z):p(z)
!************************************************************
      q(m+1,l+4,k)=-psum*cdc2
      go to 350
      endif
  300 continue
      if(l.ne.5) go to 350
      if(k.eq.m.and.k.ne.3) then
      ia=k
      psum=alamda*(-2.d0*alamda*r(ia)**2.d0-(alamda*(r(1)**2.d0+ &
     & r(2)**2.d0-2.d0*r(3)**2.d0)-1.d0)*(1.d0-2.d0*alamda*r(ia) &
     & **2.d0))
!************************************************************
!*     d(3z2-r2):del(x):p(x)
!*     d(3z2-r2):del(y):p(y)
!************************************************************
      q(m+1,l+4,k)=-psum*cdc2
      go to 350
      endif
      if(k.eq.1.and.m.eq.2) go to 310
      if(k.eq.2.and.m.eq.1) go to 310
      go to 320
  310 psum=2.d0*(alamda**2.d0)*r(1)*r(2)*(alamda*(r(1)**2.d0+ &
     & r(2)**2.d0-2.d0*r(3)**2.d0)-1.d0)
!************************************************************
!*     d(3z2-r2):del(x):p(y)     or
!*     d(3z2-r2):del(y):p(x)
!************************************************************
      q(m+1,l+4,k)=-psum*cdc2
      go to 350
  320 continue
      if(m.eq.3.and.k.ne.3) then
      ia=m
      ib=k
      go to 330
      endif
      if(k.eq.3.and.m.ne.3) then
      ia=m
      ib=k
      go to 330
      endif
      go to 340
  330 psum=2.d0*(alamda**2.d0)*r(ia)*r(ib)*(alamda*(r(1)**2.d0+ &
     & r(2)**2.d0-2.d0*r(3)**2.d0)+1.d0)
!************************************************************
!*     d(3z2-r2):del(x):p(z)     or
!*     d(3z2-r2):del(y):p(z)     or
!*     d(3z2-r2):del(z):p(x)     or
!*     d(3z2-r2):del(z):p(y)
!************************************************************
      q(m+1,l+4,k)=-psum*cdc2
      go to 350
  340 continue
      if(k.eq.3.and.m.eq.3) then
      psum=alamda*(4.d0*alamda*r(3)**2.d0-(alamda*(r(1)**2.d0+ &
     & r(2)**2.d0-2.d0*r(3)**2.d0)+2.d0)*(1.d0-2.d0*alamda*r(3)**2.d0))
!************************************************************
!*     d(3z2-r2):del(z):p(z)
!************************************************************
      q(m+1,l+4,k)=-psum*cdc2
      endif
  350 continue
      endif
! c
      if(nx1.gt.4.and.nx2.ge.4) then
      cdc1=c2/(a1*a2)
      do 530 l=1,5
      do 530 k=1,3
      do 530 m=1,3
      if(l.eq.4) go to 410
      if(l.eq.5) go to 480
      if(l.eq.k.and.k.eq.m)then
      ia=l
      ib=l+1
      ic=k
      if(l.eq.3) ib=1
      go to 360
      endif
      il=l+1
      if(l.eq.3) il=1
      if(il.eq.k.and.k.eq.m)then
      ia=l
      ib=il
      ic=k
      go to 360
      endif
      go to 370
  360 psum=(alamda**2.d0)*r(ia)*r(ib)*(3.d0-2.d0*alamda*r(ic)**2.d0)
!************************************************************
!*     d(ij):del(i):p(i)
!*     d(ij):del(j):p(j)
!************************************************************
      q(l+4,m+1,k)=psum*cdc1
      go to 530
  370 continue
      im=m-1
      if(m.eq.1) im=3
      if(l.eq.k.and.k.eq.im) then
      ia=k
      ib=m
      go to 380
      endif
      ik=k-1
      if(k.eq.1) ik=3
      if(l.eq.ik.and.l.eq.m) then
      ia=m
      ib=k
      go to 380
      endif
      go to 390
  380 psum=-alamda*(1.d0-2.d0*alamda*r(ia)**2.d0)*(1.d0-2.d0*alamda &
     & *r(ib)**2.d0)/2.d0
!************************************************************
!*     d(ij):del(i):p(j)
!*     d(ij):del(j):p(i)
!************************************************************
      q(l+4,m+1,k)=psum*cdc1
      go to 530
  390 continue
      im=m+1
      if(m.eq.3) im=1
      if(l.eq.k.and.l.eq.im) then
      ia=m
      ib=k+1
      ic=l
      if(k.eq.3) ib=1
      go to 400
      endif
      ik=k-1
      if(k.eq.1) ik=3
      if(l.eq.ik.and.l.eq.im) then
      ia=m
      ib=l
      ic=k
      go to 400
      endif
      ik=k+1
      if(k.eq.3) ik=1
      if(l.eq.ik.and.l.eq.m) then
      ia=k
      ib=m+1
      if(m.eq.3) ib=1
      ic=l
      go to 400
      endif
      im=m-1
      if(m.eq.1) im=3
      if(l.eq.ik.and.l.eq.im) then
      ia=k
      ib=l
      ic=m
      go to 400
      endif
      ik=k+1
      if(k.eq.3) ik=1
      if(l.eq.ik.and.k.eq.m) then
      ia=l
      ib=l+1
      if(l.eq.3) ib=1
      ic=k
      go to 400
      endif
      go to 410
  400 psum=(alamda**2.d0)*r(ia)*r(ib)*(1.d0-2.d0*alamda*r(ic)**2.d0)
!************************************************************
!*     d(ij):del(i):p(k)     or
!*     d(ij):del(j):p(k)     or
!*     d(ij):del(k):p(i)     or
!*     d(ij):del(k):p(j)     or
!*     d(ij):del(k):p(k)
!************************************************************
      q(l+4,m+1,k)=psum*cdc1
      go to 530
  410 continue
      if(l.ne.4) go to 480
      if(k.eq.1.and.m.eq.1) then
      ia=1
      ib=2
      sn=1.d0
      go to 420
      endif
      if(k.eq.2.and.m.eq.2) then
      ia=2
      ib=1
      sn=-1.d0
      go to 420
      endif
      go to 430
  420 psum=sn*alamda*(2.d0*alamda*r(ia)**2.d0-(1.d0-alamda*(r(ia) &
     & **2.d0-r(ib)**2.d0))*(1.d0-2.d0*alamda*r(ia)**2.d0))
!************************************************************
!*     d(x2-y2):del(x):p(x)     or
!*     d(x2-y2):del(y):p(y)
!************************************************************
      q(l+4,m+1,k)=psum*cdc1
      go to 530
  430 continue
      if(k.eq.1.and.m.eq.2) go to 440
      if(k.eq.2.and.m.eq.1) go to 440
      go to 450
  440 psum=-2.d0*(alamda**3.d0)*r(1)*r(2)*(r(1)**2.d0-r(2)**2.d0)
!************************************************************
!*     d(x2-y2):del(x):p(y)
!*     d(x2-y2):del(y):p(x)
!************************************************************
      q(l+4,m+1,k)=psum*cdc1
      go to 530
  450 continue
      if(k.eq.1.and.m.eq.3) then
      ia=1
      ci=-1.d0
      go to 460
      endif
      if(k.eq.2.and.m.eq.3) then
      ia=2
      ci=1.d0
      go to 460
      endif
      if(k.eq.3.and.m.eq.1) then
      ia=1
      ci=-1.d0
      go to 460
      endif
      if(k.eq.3.and.m.eq.2) then
      ia=2
      ci=1.d0
      go to 460
      endif
      go to 470
  460 psum=-2.d0*(alamda**2.d0)*r(ia)*r(3)*(ci &
     & +alamda*(r(1)**2.d0-r(2)**2.d0))
!************************************************************
!*     d(x2-y2):del(x):p(z)     or
!*     d(x2-y2):del(y):p(z)     or
!*     d(x2-y2):del(z):p(x)     or
!*     d(x2-y2):del(z):p(y)
!************************************************************
      q(l+4,m+1,k)=psum*cdc1
      go to 530
  470 continue
      if(k.eq.3.and.m.eq.3) then
      psum=(alamda**2.d0)*(r(1)**2.d0-r(2)**2.d0)*(1.d0-2.d0*alamda* &
     & r(3)**2.d0)
!************************************************************
!*     d(x2-y2):del(z):p(z)
!************************************************************
      q(l+4,m+1,k)=psum*cdc1
      go to 530
      endif
  480 continue
      if(l.ne.5) go to 530
      if(k.eq.m.and.k.ne.3) then
      ia=k
      psum=alamda*(-2.d0*alamda*r(ia)**2.d0-(alamda*(r(1)**2.d0+ &
     & r(2)**2.d0-2.d0*r(3)**2.d0)-1.d0)*(1.d0-2.d0*alamda*r(ia) &
     & **2.d0))
!************************************************************
!*     d(3z2-r2):del(x):p(x)
!*     d(3z2-r2):del(y):p(y)
!************************************************************
      q(l+4,m+1,k)=psum*cdc1
      go to 530
      endif
      if(k.eq.1.and.m.eq.2) go to 490
      if(k.eq.2.and.m.eq.1) go to 490
      go to 500
  490 psum=2.d0*(alamda**2.d0)*r(1)*r(2)*(alamda*(r(1)**2.d0+ &
     & r(2)**2.d0-2.d0*r(3)**2.d0)-1.d0)
!************************************************************
!*     d(3z2-r2):del(x):p(y)     or
!*     d(3z2-r2):del(y):p(x)
!************************************************************
      q(l+4,m+1,k)=psum*cdc1
      go to 530
  500 continue
      if(m.eq.3.and.k.ne.3) then
      ia=m
      ib=k
      go to 510
      endif
      if(k.eq.3.and.m.ne.3) then
      ia=m
      ib=k
      go to 510
      endif
      go to 520
  510 psum=2.d0*(alamda**2.d0)*r(ia)*r(ib)*(alamda*(r(1)**2.d0+ &
     & r(2)**2.d0-2.d0*r(3)**2.d0)+1.d0)
!************************************************************
!*     d(3z2-r2):del(x):p(z)     or
!*     d(3z2-r2):del(y):p(z)     or
!*     d(3z2-r2):del(z):p(x)     or
!*     d(3z2-r2):del(z):p(y)
!************************************************************
      q(l+4,m+1,k)=psum*cdc1
      go to 530
  520 continue
      if(k.eq.3.and.m.eq.3) then
      psum=alamda*(4.d0*alamda*r(3)**2.d0-(alamda*(r(1)**2.d0+ &
     & r(2)**2.d0-2.d0*r(3)**2.d0)+2.d0)*(1.d0-2.d0*alamda*r(3)**2.d0))
!************************************************************
!*     d(3z2-r2):del(z):p(z)
!************************************************************
      q(l+4,m+1,k)=psum*cdc1
      endif
  530 continue
      endif
! c
      if(nx1.gt.4.and.nx2.gt.4) then
      cdc=c5
      do 810 l=1,5
      do 810 m=1,l
      do 810 k=1,3
      if(l.eq.4) go to 590
      if(m.eq.4) go to 590
      if(l.eq.5) go to 670
      if(m.eq.5) go to 670
      if(l.eq.m.and.k.eq.m) then
      ia=l
      ib=l+1
      if(l.eq.3) ib=1
      go to 540
      endif
      ik=k-1
      if(k.eq.1) ik=3
      if(l.eq.m.and.ik.eq.m) then
      ia=k
      ib=l
      go to 540
      endif
      go to 550
  540 psum=r(ia)*(1.d0-2.d0*alamda*r(ib)**2.d0)*(3.d0-2.d0*alamda* &
     & r(ia)**2.d0)/2.d0
!************************************************************
!*     d(ij):del(i):d(ij)     or
!*     d(ij):del(j):d(ij)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  550 continue
      im=m-1
      if(m.eq.1) im=3
      if(l.eq.k.and.l.eq.im) then
      ia=m+1
      if(m.eq.3) ia=1
      ib=l
      ic=m
      go to 560
      endif
      ik=k-1
      if(k.eq.1) ik=3
      im=m+1
      if(m.eq.3) im=1
      if(l.eq.ik.and.l.eq.im) then
      ia=m
      ib=l
      ic=k
      go to 560
      endif
      ik=k+1
      if(k.eq.3) ik=1
      if(l.eq.m.and.l.eq.ik) then
      ia=k
      ib=l
      ic=l+1
      if(l.eq.3) ic=1
      go to 560
      endif
      ik=k+1
      if(k.eq.3) ik=1
      im=m-1
      if(m.eq.1) im=3
      if(l.eq.im.and.l.eq.ik) then
      ia=l
      ib=k
      ic=m
      go to 560
      endif
      ik=k+1
      if(k.eq.3) ik=1
      im=m+1
      if(m.eq.3) im=1
      if(l.eq.im.and.l.eq.ik) then
      ia=l+1
      if(l.eq.3) ia=1
      ib=k
      ic=im
      go to 560
      endif
      go to 570
  560 psum=r(ia)*(1.d0-2.d0*alamda*r(ib)**2.d0)*(1.d0-2.d0*alamda &
     & *r(ic)**2.d0)/2.d0
!*************************************************************
!*     d(ij):del(i):d(jk)
!*     d(ij):del(j):d(ki)
!*     d(ij):del(k):d(ij)
!*     d(ij):del(k):d(jk)
!*     d(ij):del(k):d(ki)
!*************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  570 continue
      im=m+1
      if(m.eq.3) im=1
      if(l.eq.k.and.l.eq.im) then
      ia=k
      go to 580
      endif
      ik=k-1
      if(k.eq.1) ik=3
      im=m-1
      if(m.eq.1) im=3
      if(l.eq.im.and.l.eq.ik) then
      ia=k
      go to 580
      endif
      go to 590
  580 psum=-alamda*r(1)*r(2)*r(3)*(3.d0-2.d0*alamda*r(ia)**2.d0)
!************************************************************
!*     d(ij):del(i):d(ki)     or
!*     d(ij):del(j):d(jk)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  590 continue
      if(l.eq.4.and.m.eq.4) go to 650
      if(m.eq.1.and.k.eq.1) then
      ia=2
      ib=1
      sn=1.d0
      go to 600
      endif
      if(m.eq.1.and.k.eq.2) then
      ia=1
      ib=2
      sn=-1.d0
      go to 600
      endif
      go to 610
  600 psum=-alamda*r(ia)*((1.d0-2.d0*alamda*r(ib)**2.d0)*(r(1)**2.d0 &
     & -r(2)**2.d0)+sn*2.d0*r(ib)**2.d0)
!************************************************************
!*     d(x2-y2):del(x):d(xy)     or
!*     d(x2-y2):del(y):d(xy)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  610 continue
      if(m.eq.1.and.k.eq.3) go to 620
      if(m.eq.2.and.k.eq.1) go to 620
      if(m.eq.3.and.k.eq.2) go to 620
      go to 630
  620 psum=2.d0*(alamda**2.d0)*r(1)*r(2)*r(3)*(r(1)**2.d0-r(2)**2.d0)
!************************************************************
!*     d(x2-y2):del(z):d(xy)     or
!*     d(x2-y2):del(x):d(yz)     or
!*     d(x2-y2):del(y):d(zx)     or
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  630 continue
      if(m.eq.2.and.k.eq.2) then
      ia=3
      ib=2
      ci=1.d0
      sn=1.d0
      go to 640
      endif
      if(m.eq.2.and.k.eq.3) then
      ia=2
      ib=3
      ci=0.d0
      sn=1.d0
      go to 640
      endif
      if(m.eq.3.and.k.eq.1) then
      ia=3
      ib=1
      ci=-1.d0
      sn=-1.d0
      go to 640
      endif
      if(m.eq.3.and.k.eq.3) then
      ia=1
      ib=3
      ci=0.d0
      sn=-1.d0
      go to 640
      endif
      go to 650
  640 psum=r(ia)*(ci*2.d0*alamda*r(ib)**2.d0+(1.d0-2.d0*alamda*r(ib) &
     & **2.d0)*(-1.d0*sn-alamda*(r(1)**2.d0-r(2)**2.d0)))
!************************************************************
!*     d(x2-y2):del(y):d(yz)     or
!*     d(x2-y2):del(z):d(yz)     or
!*     d(x2-y2):del(x):d(zx)     or
!*     d(x2-y2):del(z):d(zx)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  650 continue
      if(l.ne.4) go to 670
      if(m.eq.4.and.k.eq.1) then
      ci=2.d0
      sn=-1.d0
      go to 660
      endif
      if(m.eq.4.and.k.eq.2) then
      ci=2.d0
      sn=1.d0
      go to 660
      endif
      if(m.eq.4.and.k.eq.3) then
      ci=0.d0
      sn=1.d0
      go to 660
      endif
      go to 670
  660 psum=2.d0*r(k)*(1.d0+ci-2.d0*alamda*(r(1)**2.d0+r(2)**2.d0) &
     & +(alamda*(r(1)**2.-r(2)**2.)*(ci*sn+alamda &
     & *(r(1)**2.-r(2)**2.))))
!************************************************************
!*     d(x2-y2):del(i):d(x2-y2)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  670 continue
      if(l.ne.5) go to 810
      if(l.eq.5.and.m.eq.5) go to 780
      if(m.eq.1.and.k.eq.1) then
      ia=2
      ib=1
      go to 680
      endif
      if(m.eq.1.and.k.eq.2) then
      ia=1
      ib=2
      go to 680
      endif
      go to 690
  680 psum=-r(ia)*(-(alamda*(r(1)**2.d0+r(2)**2.d0-2.d0*r(3)**2.d0) &
     & -2.d0)*(1.d0-2.d0*alamda*r(ib)**2.d0)-2.d0*alamda*r(ib)**2.d0)
!************************************************************
!*     d(3z2-r2):del(x):d(xy)     or
!*     d(3z2-r2):del(y):d(xy)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  690 continue
      if(m.eq.1.and.k.eq.3) go to 700
      if(m.eq.2.and.k.eq.1) go to 700
      if(m.eq.3.and.k.eq.2) go to 700
      go to 710
  700 psum=-2.d0*(alamda**2.d0)*r(1)*r(2)*r(3)*(r(1)**2.d0+r(2)**2.d0 &
     & -2.d0*r(3)**2.d0)
!************************************************************
!*     d(3z2-r2):del(z):d(xy)     or
!*     d(3z2-r2):del(x):d(yz)     or
!*     d(3z2-r2):del(y):d(zx)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  710 continue
      if(m.eq.2.and.k.eq.2) then
      ia=2
      sn=1.d0
      go to 720
      endif
      if(m.eq.3.and.k.eq.1) then
      ia=1
      sn=1.d0
      go to 720
      endif
      go to 730
  720 psum=r(3)*(alamda*(1.d0-2.d0*alamda*r(ia)**2.d0)*(r(1)**2.d0 &
     & +r(2)**2.d0-2.d0*r(3)**2.d0)+sn)
!************************************************************
!*     d(3z2-r2):del(y):d(yz)     or
!*     d(3z2-r2):del(x):d(zx)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  730 continue
      if(m.eq.2.and.k.eq.3) then
      ia=2
      go to 740
      endif
      if(m.eq.3.and.k.eq.3) then
      ia=1
      go to 740
      endif
      go to 750
  740 psum=r(ia)*(1.d0-6.d0*alamda*r(3)**2.d0+(r(1)**2.d0+r(2)**2.d0 &
     & -2.d0*r(3)**2.d0)*(alamda-2.d0*alamda*alamda*r(3)**2.d0))
!************************************************************
!*     d(3z2-r2):del(z):d(yz)     or
!*     d(3z2-r2):del(z):d(zx)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  750 continue
      if(m.eq.4.and.k.eq.1) then
      ia=1
      ib=2
      sn=-1.d0
      go to 760
      endif
      if(m.eq.4.and.k.eq.2) then
      ia=2
      ib=1
      sn=1.d0
      go to 760
      endif
      go to 770
  760 psum=sn*r(ia)*(-(alamda*(r(1)**2.d0+r(2)**2.d0-2.d0*r(3) &
     & **2.d0)-1.d0)*(3.d0-2.d0*alamda*r(ia)**2.d0)-(4.d0* &
     & alamda*r(ia)**2.d0-3.d0)+(alamda*(r(1)**2.d0+r(2)**2.d0 &
     & -2.d0*r(3)**2.d0)-2.d0)*(1.d0-2.d0*r(ib)**2.d0)+2.d0* &
     & alamda*r(ib)**2.d0)
!************************************************************
!*     d(3z2-r2):del(x):d(x2-y2)     or
!*     d(3z2-r2):del(y):d(x2-y2)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  770 continue
      if(m.eq.4.and.k.eq.3) then
      psum=-2.d0*(alamda**2.d0)*(r(1)**2.d0+r(2)**2.d0-2.d0* &
     & r(3)**2.d0)*(r(1)**2.d0-r(2)**2.d0)*r(3)
!************************************************************
!*     d(3z2-r2):del(z):d(x2-y2)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
      endif
  780 continue
      if(m.eq.5.and.k.eq.1) then
      ia=1
      ib=2
      go to 790
      endif
      if(m.eq.5.and.k.eq.2) then
      ia=2
      ib=1
      go to 790
      endif
      go to 800
  790 psum=r(ia)*(2.d0*((alamda*(r(1)**2.d0+r(2)**2.d0-2.d0*r(3)**2.d0)&
     & +1.d0)*(1.d0-2.d0*alamda*r(3)**2.d0)-4.d0*alamda*r(3)**2.d0) &
     & -(4.d0*alamda*r(ia)**2.d0-3.d0+(3.d0-2.d0*alamda*r(ia)**2.d0) &
     & *(alamda*(r(1)**2.d0+r(2)**2.d0-2.d0*r(3)**2.d0)-1.d0))-((1.d0 &
     & -2.d0*alamda*r(ib)**2.d0)*(alamda*(r(1)**2.d0+r(2)**2.d0-2.d0 &
     & *r(3)**2.d0)-2.d0)+2.d0*alamda*r(ib)**2.d0))
!************************************************************
!*     d(3z2-r2):del(x):d(3z2-r2)     or
!*     d(3z2-r2):del(y):d(3z2-r2)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      go to 810
  800 continue
      if(m.eq.5.and.k.eq.3) then
      psum=2.d0*r(3)*(11.d0+4.d0*alamda*(r(1)*r(1)+r(2)*r(2)- &
     & 2.d0*r(3)*r(3))-2.d0*alamda*(4.d0*r(3)*r(3)+r(1)*r(1)+ &
     & r(2)*r(2))+alamda*alamda*(2.d0*r(3)*r(3)-r(1)*r(1)- &
     & r(2)*r(2))**2)
!************************************************************
!*     d(3z2-r2):del(z):d(3z2-r2)
!************************************************************
      q(l+4,m+4,k)=psum*cdc
      endif
  810 continue
      do 830 l=1,5
      do 830 m=1,l
      do 830 k=1,3
      q(m+4,l+4,k)=q(l+4,m+4,k)
  830 continue
      endif
! c
      if(nx1.lt.16.and.nx2.lt.16) return
! c
      abx=r(1)
      aby=r(2)
      abz=r(3)
! c
      a12=a1*a1
      a22=a2*a2
      a13=a12*a1
      a23=a22*a2
      alp=a1*a2/(a1+a2)
      alp2=2.0*alp
      alp4=4.0*alp
      alp6=6.0*alp
      alp8=8.0*alp
      alp10=10.0*alp
      alp12=12.0*alp
      alpx=alp2*abx
      alpy=alp2*aby
      alpz=alp2*abz
      expa=alp*rab2
      ro=dexp(-expa)
      w=pi/(a1+a2)
      delt=w*w*w
      delt=dsqrt(delt)
!c1
      rox=-alpx*ro
      roy=-alpy*ro
      roz=-alpz*ro
!c2
      rox2=-alp2*ro-alpx*rox
      roy2=-alp2*ro-alpy*roy
      roz2=-alp2*ro-alpz*roz
      roxy=-alpy*rox
      roxz=-alpz*rox
      royz=-alpz*roy
!c3
      rox3=-alp4*rox-alpx*rox2
      roy3=-alp4*roy-alpy*roy2
      roz3=-alp4*roz-alpz*roz2
      rox2y=-alpy*rox2
      rox2z=-alpz*rox2
      roxy2=-alpx*roy2
      roxz2=-alpx*roz2
      roy2z=-alpz*roy2
      royz2=-alpy*roz2
      roxyz=-alpx*royz
!c4
      rox4=-alp6*rox2-alpx*rox3
      roy4=-alp6*roy2-alpy*roy3
      roz4=-alp6*roz2-alpz*roz3
      rox3y=-alpy*rox3
      rox3z=-alpz*rox3
      roxy3=-alpx*roy3
      roxz3=-alpx*roz3
      roy3z=-alpz*roy3
      royz3=-alpy*roz3
      rox2y2=-alp2*roy2-alpx*roxy2
      rox2z2=-alp2*roz2-alpx*roxz2
      roy2z2=-alp2*roz2-alpy*royz2
      rox2yz=-alpz*rox2y
      roxy2z=-alpz*roxy2
      roxyz2=-alpy*roxz2
!c5
      rox5=-alp8*rox3-alpx*rox4
      roy5=-alp8*roy3-alpy*roy4
      roz5=-alp8*roz3-alpz*roz4
      rox4y=-alpy*rox4
      rox4z=-alpz*rox4
      roxy4=-alpx*roy4
      roxz4=-alpx*roz4
      roy4z=-alpz*roy4
      royz4=-alpy*roz4
      rox3y2=-alp4*roxy2-alpx*rox2y2
      rox3z2=-alp4*roxz2-alpx*rox2z2
      rox2y3=-alp4*rox2y-alpy*rox2y2
      rox2z3=-alp4*rox2z-alpz*rox2z2
      roy3z2=-alp4*royz2-alpy*roy2z2
      roy2z3=-alp4*roy2z-alpz*roy2z2
      rox3yz=-alpz*rox3y
      roxy3z=-alpz*roy3z
      roxyz3=-alpx*royz3
      rox2y2z=-alpz*rox2y2
      rox2yz2=-alpy*rox2z2
      roxy2z2=-alpx*roy2z2
!c6
      rox6=-alp10*rox4-alpx*rox5
      roy6=-alp10*roy4-alpy*roy5
      roz6=-alp10*roz4-alpz*roz5
      rox5y=-alp8*rox3y-alpx*rox4y
      rox5z=-alp8*rox3z-alpx*rox4z
      roxy5=-alp8*roxy3-alpy*roxy4
      roxz5=-alp8*roxz3-alpz*roxz4
      roy5z=-alp8*roy3z-alpy*roy4z
      royz5=-alp8*royz3-alpz*royz4
      rox4y2=-alp6*rox2y2-alpx*rox3y2
      rox4z2=-alp6*rox2z2-alpx*rox3z2
      rox2y4=-alp6*rox2y2-alpy*rox2y3
      rox2z4=-alp6*rox2z2-alpz*rox2z3
      roy4z2=-alp6*roy2z2-alpy*roy3z2
      roy2z4=-alp6*roy2z2-alpz*roy2z3
      rox4yz=-alpz*rox4y
      roxy4z=-alpz*roxy4
      roxyz4=-alpy*roxz4
      rox3y3=-alp4*roxy3-alpx*rox2y3
      rox3z3=-alp4*roxz3-alpx*rox2z3
      roy3z3=-alp4*roy3z-alpz*roy3z2
      rox3y2z=-alpz*rox3y2
      rox3yz2=-alpy*rox3z2
      rox2y3z=-alpz*rox2y3
      rox2yz3=-alpy*rox2z3
      roxy3z2=-alpx*roy3z2
      roxy2z3=-alpx*roy2z3
      rox2y2z2=-alp2*roy2z2-alpx*roxy2z2
!c7
      rox7=-alp12*rox5-alpx*rox6
      roy7=-alp12*roy5-alpy*roy6
      roz7=-alp12*roz5-alpz*roz6
      rox6y=-alp10*rox4y-alpx*rox5y
      rox6z=-alp10*rox4z-alpx*rox5z
      roxy6=-alp10*roxy4-alpy*roxy5
      roxz6=-alp10*roxz4-alpz*roxz5
      roy6z=-alp10*roy4z-alpy*roy5z
      royz6=-alp10*royz4-alpz*royz5
      rox5y2=-alp8*rox3y2-alpx*rox4y2
      rox5z2=-alp8*rox3z2-alpx*rox4z2
      rox2y5=-alp8*rox2y3-alpy*rox2y4
      rox2z5=-alp8*rox2z3-alpz*rox2z4
      roy5z2=-alp8*roy3z2-alpy*roy4z2
      roy2z5=-alp8*roy2z3-alpz*roy2z4
      roxyz5=-alp8*roxyz3-alpz*roxyz4
      rox5yz=-alp8*rox3yz-alpx*rox4yz
      roxy5z=-alp8*roxy3z-alpy*roxy4z
      rox4y3=-alp6*rox2y3-alpx*rox3y3
      rox4z3=-alp6*rox2z3-alpx*rox3z3
      rox3y4=-alp6*rox3y2-alpy*rox3y3
      rox3z4=-alp6*rox3z2-alpz*rox3z3
      roy4z3=-alp6*roy2z3-alpy*roy3z3
      roy3z4=-alp6*roy3z2-alpz*roy3z3
      rox4y2z=-alp6*rox2y2z-alpx*rox3y2z
      rox4yz2=-alp6*rox2yz2-alpx*rox3yz2
      rox2y4z=-alp6*rox2y2z-alpy*rox2y3z
      roxy2z4=-alp6*roxy2z2-alpz*roxy2z3
      rox2yz4=-alp6*rox2yz2-alpz*rox2yz3
      roxy4z2=-alp6*roxy2z2-alpy*roxy3z2
      rox3y3z=-alp4*roxy3z-alpx*rox2y3z
      roxy3z3=-alp4*roxy3z-alpz*roxy3z2
      rox3yz3=-alp4*roxyz3-alpx*rox2yz3
      rox3y2z2=-alp4*roxy2z2-alpx*rox2y2z2
      rox2y3z2=-alp4*rox2yz2-alpy*rox2y2z2
      rox2y2z3=-alp4*rox2y2z-alpz*rox2y2z2
! c
      if(nx1.lt.16) go to 1020
      s1=delt/(8.d0*a13)
      s2=delt/(4.d0*a12)
      s3=3.d0*delt/(4.d0*a12)
!c1
      x3xs=s1*rox4+s3*rox2
      y3ys=s1*roy4+s3*roy2
      z3zs=s1*roz4+s3*roz2
!c2
      x3ys=s1*rox3y+s3*roxy
      x3zs=s1*rox3z+s3*roxz
      y3xs=s1*roxy3+s3*roxy
      y3zs=s1*roy3z+s3*royz
      z3xs=s1*roxz3+s3*roxz
      z3ys=s1*royz3+s3*royz
!c3
      x2yxs=s1*rox3y+s2*roxy
      y2xys=s1*roxy3+s2*roxy
      z2xzs=s1*roxz3+s2*roxz
      y2zys=s1*roy3z+s2*royz
      z2yzs=s1*royz3+s2*royz
      x2zxs=s1*rox3z+s2*roxz
!c4
      x2yys=s1*rox2y2+s2*roy2
      x2zzs=s1*rox2z2+s2*roz2
      y2xxs=s1*rox2y2+s2*rox2
      y2zzs=s1*roy2z2+s2*roz2
      z2xxs=s1*rox2z2+s2*rox2
      z2yys=s1*roy2z2+s2*roy2
!c5
      x2yzs=s1*rox2yz+s2*royz
      x2zys=s1*rox2yz+s2*royz
      y2xzs=s1*roxy2z+s2*roxz
      y2zxs=s1*roxy2z+s2*roxz
      z2xys=s1*roxyz2+s2*roxy
      z2yxs=s1*roxyz2+s2*roxy
!c6
      xyzxs=s1*rox2yz
      xyzys=s1*roxy2z
      xyzzs=s1*roxyz2
! c
      fs(1,1)=xyzxs
      fs(1,2)=xyzys
      fs(1,3)=xyzzs
      fs(2,1)=(x2zxs-y2zxs)
      fs(2,2)=(x2zys-y2zys)
      fs(2,3)=(x2zzs-y2zzs)
      fs(3,1)=(x3xs-3.d0*y2xxs)
      fs(3,2)=(x3ys-3.d0*y2xys)
      fs(3,3)=(x3zs-3.d0*y2xzs)
      fs(4,1)=(3.d0*x2yxs-y3xs)
      fs(4,2)=(3.d0*x2yys-y3ys)
      fs(4,3)=(3.d0*x2yzs-y3zs)
      fs(5,1)=2.d0*z3xs-3.d0*x2zxs-3.d0*y2zxs
      fs(5,2)=2.d0*z3ys-3.d0*x2zys-3.d0*y2zys
      fs(5,3)=2.d0*z3zs-3.d0*x2zzs-3.d0*y2zzs
      fs(6,1)=(4.d0*z2xxs-x3xs-y2xxs)
      fs(6,2)=(4.d0*z2xys-x3ys-y2xys)
      fs(6,3)=(4.d0*z2xzs-x3zs-y2xzs)
      fs(7,1)=(4.d0*z2yxs-x2yxs-y3xs)
      fs(7,2)=(4.d0*z2yys-x2yys-y3ys)
      fs(7,3)=(4.d0*z2yzs-x2yzs-y3zs)
! c
      do 1010 k=1,7
      do 1010 l=1,3
      q(9+k,1,l)=fs(k,l)
 1010 continue
 1020 continue
! c
      if(nx2.lt.16) go to 1040
      s1=-delt/(8.d0*a23)
      s2=-delt/(4.d0*a22)
      s3=-3.d0*delt/(4.d0*a22)
!c1
      x3xs=s1*rox4+s3*rox2
      y3ys=s1*roy4+s3*roy2
      z3zs=s1*roz4+s3*roz2
!c2
      x3ys=s1*rox3y+s3*roxy
      x3zs=s1*rox3z+s3*roxz
      y3xs=s1*roxy3+s3*roxy
      y3zs=s1*roy3z+s3*royz
      z3xs=s1*roxz3+s3*roxz
      z3ys=s1*royz3+s3*royz
!c3
      x2yxs=s1*rox3y+s2*roxy
      y2xys=s1*roxy3+s2*roxy
      z2xzs=s1*roxz3+s2*roxz
      y2zys=s1*roy3z+s2*royz
      z2yzs=s1*royz3+s2*royz
      x2zxs=s1*rox3z+s2*roxz
!c4
      x2yys=s1*rox2y2+s2*roy2
      x2zzs=s1*rox2z2+s2*roz2
      y2xxs=s1*rox2y2+s2*rox2
      y2zzs=s1*roy2z2+s2*roz2
      z2xxs=s1*rox2z2+s2*rox2
      z2yys=s1*roy2z2+s2*roy2
!c5
      x2yzs=s1*rox2yz+s2*royz
      x2zys=s1*rox2yz+s2*royz
      y2xzs=s1*roxy2z+s2*roxz
      y2zxs=s1*roxy2z+s2*roxz
      z2xys=s1*roxyz2+s2*roxy
      z2yxs=s1*roxyz2+s2*roxy
!c6
      xyzxs=s1*rox2yz
      xyzys=s1*roxy2z
      xyzzs=s1*roxyz2
! c
      sf(1,1)=xyzxs
      sf(1,2)=xyzys
      sf(1,3)=xyzzs
      sf(2,1)=(x2zxs-y2zxs)
      sf(2,2)=(x2zys-y2zys)
      sf(2,3)=(x2zzs-y2zzs)
      sf(3,1)=(x3xs-3.d0*y2xxs)
      sf(3,2)=(x3ys-3.d0*y2xys)
      sf(3,3)=(x3zs-3.d0*y2xzs)
      sf(4,1)=(3.d0*x2yxs-y3xs)
      sf(4,2)=(3.d0*x2yys-y3ys)
      sf(4,3)=(3.d0*x2yzs-y3zs)
      sf(5,1)=2.d0*z3xs-3.d0*x2zxs-3.d0*y2zxs
      sf(5,2)=2.d0*z3ys-3.d0*x2zys-3.d0*y2zys
      sf(5,3)=2.d0*z3zs-3.d0*x2zzs-3.d0*y2zzs
      sf(6,1)=(4.d0*z2xxs-x3xs-y2xxs)
      sf(6,2)=(4.d0*z2xys-x3ys-y2xys)
      sf(6,3)=(4.d0*z2xzs-x3zs-y2xzs)
      sf(7,1)=(4.d0*z2yxs-x2yxs-y3xs)
      sf(7,2)=(4.d0*z2yys-x2yys-y3ys)
      sf(7,3)=(4.d0*z2yzs-x2yzs-y3zs)
! c
      do 1030 k=1,7
      do 1030 l=1,3
      q(1,9+k,l)=sf(k,l)
 1030 continue
 1040 continue
! c
      if(nx1.lt.16.or.nx2.lt.4) go to 1050
      p1=delt/(16.d0*a13*a2)
      p2=delt/(8.d0*a12*a2)
      p3=3.d0*delt/(8.d0*a12*a2)
!c1
      x3xx=p1*rox5+p3*rox3
      y3yy=p1*roy5+p3*roy3
      z3zz=p1*roz5+p3*roz3
!c2
      x3yy=p1*rox3y2+p3*roxy2
      x3zz=p1*rox3z2+p3*roxz2
      y3xx=p1*rox2y3+p3*rox2y
      y3zz=p1*roy3z2+p3*royz2
      z3xx=p1*rox2z3+p3*rox2z
      z3yy=p1*roy2z3+p3*roy2z
!c3
      x3xy=p1*rox4y+p3*rox2y
      x3xz=p1*rox4z+p3*rox2z
      y3yx=p1*roxy4+p3*roxy2
      y3yz=p1*roy4z+p3*roy2z
      z3zx=p1*roxz4+p3*roxz2
      z3zy=p1*royz4+p3*royz2
!c4
      x3yx=p1*rox4y+p3*rox2y
      x3zx=p1*rox4z+p3*rox2z
      y3xy=p1*roxy4+p3*roxy2
      y3zy=p1*roy4z+p3*roy2z
      z3xz=p1*roxz4+p3*roxz2
      z3yz=p1*royz4+p3*royz2
!c5
      x3yz=p1*rox3yz+p3*roxyz
      x3zy=p1*rox3yz+p3*roxyz
      y3xz=p1*roxy3z+p3*roxyz
      y3zx=p1*roxy3z+p3*roxyz
      z3xy=p1*roxyz3+p3*roxyz
      z3yx=p1*roxyz3+p3*roxyz
!c6
      x2yxx=p1*rox4y+p2*rox2y
      x2zxx=p1*rox4z+p2*rox2z
      y2xyy=p1*roxy4+p2*roxy2
      y2zyy=p1*roy4z+p2*roy2z
      z2xzz=p1*roxz4+p2*roxz2
      z2yzz=p1*royz4+p2*royz2
!c7
      x2yyy=p1*rox2y3+p2*roy3
      x2zzz=p1*rox2z3+p2*roz3
      y2xxx=p1*rox3y2+p2*rox3
      y2zzz=p1*roy2z3+p2*roz3
      z2xxx=p1*rox3z2+p2*rox3
      z2yyy=p1*roy3z2+p2*roy3
!c8
      x2yxy=p1*rox3y2+p2*roxy2
      x2zxz=p1*rox3z2+p2*roxz2
      y2xyx=p1*rox2y3+p2*rox2y
      y2zyz=p1*roy3z2+p2*royz2
      z2xzx=p1*rox2z3+p2*rox2z
      z2yzy=p1*roy2z3+p2*roy2z
!c9
      x2yyx=p1*rox3y2+p2*roxy2
      x2zzx=p1*rox3z2+p2*roxz2
      y2xxy=p1*rox2y3+p2*rox2y
      y2zzy=p1*roy3z2+p2*royz2
      z2xxz=p1*rox2z3+p2*rox2z
      z2yyz=p1*roy2z3+p2*roy2z
!c10
      x2yxz=p1*rox3yz+p2*roxyz
      x2zxy=p1*rox3yz+p2*roxyz
      y2xyz=p1*roxy3z+p2*roxyz
      y2zyx=p1*roxy3z+p2*roxyz
      z2xzy=p1*roxyz3+p2*roxyz
      z2yzx=p1*roxyz3+p2*roxyz
!c11
      x2yyz=p1*rox2y2z+p2*roy2z
      x2zzy=p1*rox2yz2+p2*royz2
      y2xxz=p1*rox2y2z+p2*rox2z
      y2zzx=p1*roxy2z2+p2*roxz2
      z2xxy=p1*rox2yz2+p2*rox2y
      z2yyx=p1*roxy2z2+p2*roxy2
!c12
      x2yzx=p1*rox3yz+p2*roxyz
      x2zyx=p1*rox3yz+p2*roxyz
      y2xzy=p1*roxy3z+p2*roxyz
      y2zxy=p1*roxy3z+p2*roxyz
      z2xyz=p1*roxyz3+p2*roxyz
      z2yxz=p1*roxyz3+p2*roxyz
!c13
      x2yzy=p1*rox2y2z+p2*roy2z
      x2zyz=p1*rox2yz2+p2*royz2
      y2xzx=p1*rox2y2z+p2*rox2z
      y2zxz=p1*roxy2z2+p2*roxz2
      z2xyx=p1*rox2yz2+p2*rox2y
      z2yxy=p1*roxy2z2+p2*roxy2
!c14
      x2yzz=p1*rox2yz2+p2*royz2
      x2zyy=p1*rox2y2z+p2*roy2z
      y2xzz=p1*roxy2z2+p2*roxz2
      y2zxx=p1*rox2y2z+p2*rox2z
      z2xyy=p1*roxy2z2+p2*roxy2
      z2yxx=p1*rox2yz2+p2*rox2y
!c15
      xyzxx=p1*rox3yz
      xyzyy=p1*roxy3z
      xyzzz=p1*roxyz3
!c16
      xyzxy=p1*rox2y2z
      xyzxz=p1*rox2yz2
      xyzyx=p1*rox2y2z
      xyzyz=p1*roxy2z2
      xyzzx=p1*rox2yz2
      xyzzy=p1*roxy2z2
! c
      fp(1,1,1)=xyzxx
      fp(1,1,2)=xyzyx
      fp(1,1,3)=xyzzx
      fp(1,2,1)=xyzxy
      fp(1,2,2)=xyzyy
      fp(1,2,3)=xyzzy
      fp(1,3,1)=xyzxz
      fp(1,3,2)=xyzyz
      fp(1,3,3)=xyzzz
      fp(2,1,1)=(x2zxx-y2zxx)
      fp(2,1,2)=(x2zyx-y2zyx)
      fp(2,1,3)=(x2zzx-y2zzx)
      fp(2,2,1)=(x2zxy-y2zxy)
      fp(2,2,2)=(x2zyy-y2zyy)
      fp(2,2,3)=(x2zzy-y2zzy)
      fp(2,3,1)=(x2zxz-y2zxz)
      fp(2,3,2)=(x2zyz-y2zyz)
      fp(2,3,3)=(x2zzz-y2zzz)
      fp(3,1,1)=(x3xx-3.d0*y2xxx)
      fp(3,1,2)=(x3yx-3.d0*y2xyx)
      fp(3,1,3)=(x3zx-3.d0*y2xzx)
      fp(3,2,1)=(x3xy-3.d0*y2xxy)
      fp(3,2,2)=(x3yy-3.d0*y2xyy)
      fp(3,2,3)=(x3zy-3.d0*y2xzy)
      fp(3,3,1)=(x3xz-3.d0*y2xxz)
      fp(3,3,2)=(x3yz-3.d0*y2xyz)
      fp(3,3,3)=(x3zz-3.d0*y2xzz)
      fp(4,1,1)=(3.d0*x2yxx-y3xx)
      fp(4,1,2)=(3.d0*x2yyx-y3yx)
      fp(4,1,3)=(3.d0*x2yzx-y3zx)
      fp(4,2,1)=(3.d0*x2yxy-y3xy)
      fp(4,2,2)=(3.d0*x2yyy-y3yy)
      fp(4,2,3)=(3.d0*x2yzy-y3zy)
      fp(4,3,1)=(3.d0*x2yxz-y3xz)
      fp(4,3,2)=(3.d0*x2yyz-y3yz)
      fp(4,3,3)=(3.d0*x2yzz-y3zz)
      fp(5,1,1)=2.d0*z3xx-3.d0*x2zxx-3.d0*y2zxx
      fp(5,1,2)=2.d0*z3yx-3.d0*x2zyx-3.d0*y2zyx
      fp(5,1,3)=2.d0*z3zx-3.d0*x2zzx-3.d0*y2zzx
      fp(5,2,1)=2.d0*z3xy-3.d0*x2zxy-3.d0*y2zxy
      fp(5,2,2)=2.d0*z3yy-3.d0*x2zyy-3.d0*y2zyy
      fp(5,2,3)=2.d0*z3zy-3.d0*x2zzy-3.d0*y2zzy
      fp(5,3,1)=2.d0*z3xz-3.d0*x2zxz-3.d0*y2zxz
      fp(5,3,2)=2.d0*z3yz-3.d0*x2zyz-3.d0*y2zyz
      fp(5,3,3)=2.d0*z3zz-3.d0*x2zzz-3.d0*y2zzz
      fp(6,1,1)=(4.d0*z2xxx-x3xx-y2xxx)
      fp(6,1,2)=(4.d0*z2xyx-x3yx-y2xyx)
      fp(6,1,3)=(4.d0*z2xzx-x3zx-y2xzx)
      fp(6,2,1)=(4.d0*z2xxy-x3xy-y2xxy)
      fp(6,2,2)=(4.d0*z2xyy-x3yy-y2xyy)
      fp(6,2,3)=(4.d0*z2xzy-x3zy-y2xzy)
      fp(6,3,1)=(4.d0*z2xxz-x3xz-y2xxz)
      fp(6,3,2)=(4.d0*z2xyz-x3yz-y2xyz)
      fp(6,3,3)=(4.d0*z2xzz-x3zz-y2xzz)
      fp(7,1,1)=(4.d0*z2yxx-x2yxx-y3xx)
      fp(7,1,2)=(4.d0*z2yyx-x2yyx-y3yx)
      fp(7,1,3)=(4.d0*z2yzx-x2yzx-y3zx)
      fp(7,2,1)=(4.d0*z2yxy-x2yxy-y3xy)
      fp(7,2,2)=(4.d0*z2yyy-x2yyy-y3yy)
      fp(7,2,3)=(4.d0*z2yzy-x2yzy-y3zy)
      fp(7,3,1)=(4.d0*z2yxz-x2yxz-y3xz)
      fp(7,3,2)=(4.d0*z2yyz-x2yyz-y3yz)
      fp(7,3,3)=(4.d0*z2yzz-x2yzz-y3zz)
! c
      do 1045 k=1,7
      do 1045 m=1,3
      do 1045 l=1,3
      q(9+k,1+m,l)=fp(k,m,l)
 1045 continue
 1050 continue
! c
      if(nx1.lt.4.or.nx2.lt.16) go to 1060
      p1=delt/(16.d0*a23*a1)
      p2=delt/(8.d0*a22*a1)
      p3=3.d0*delt/(8.d0*a22*a1)
!c1
      x3xx=p1*rox5+p3*rox3
      y3yy=p1*roy5+p3*roy3
      z3zz=p1*roz5+p3*roz3
!c2
      x3yy=p1*rox3y2+p3*roxy2
      x3zz=p1*rox3z2+p3*roxz2
      y3xx=p1*rox2y3+p3*rox2y
      y3zz=p1*roy3z2+p3*royz2
      z3xx=p1*rox2z3+p3*rox2z
      z3yy=p1*roy2z3+p3*roy2z
!c3
      x3xy=p1*rox4y+p3*rox2y
      x3xz=p1*rox4z+p3*rox2z
      y3yx=p1*roxy4+p3*roxy2
      y3yz=p1*roy4z+p3*roy2z
      z3zx=p1*roxz4+p3*roxz2
      z3zy=p1*royz4+p3*royz2
!c4
      x3yx=p1*rox4y+p3*rox2y
      x3zx=p1*rox4z+p3*rox2z
      y3xy=p1*roxy4+p3*roxy2
      y3zy=p1*roy4z+p3*roy2z
      z3xz=p1*roxz4+p3*roxz2
      z3yz=p1*royz4+p3*royz2
!c5
      x3yz=p1*rox3yz+p3*roxyz
      x3zy=p1*rox3yz+p3*roxyz
      y3xz=p1*roxy3z+p3*roxyz
      y3zx=p1*roxy3z+p3*roxyz
      z3xy=p1*roxyz3+p3*roxyz
      z3yx=p1*roxyz3+p3*roxyz
!c6
      x2yxx=p1*rox4y+p2*rox2y
      x2zxx=p1*rox4z+p2*rox2z
      y2xyy=p1*roxy4+p2*roxy2
      y2zyy=p1*roy4z+p2*roy2z
      z2xzz=p1*roxz4+p2*roxz2
      z2yzz=p1*royz4+p2*royz2
!c7
      x2yyy=p1*rox2y3+p2*roy3
      x2zzz=p1*rox2z3+p2*roz3
      y2xxx=p1*rox3y2+p2*rox3
      y2zzz=p1*roy2z3+p2*roz3
      z2xxx=p1*rox3z2+p2*rox3
      z2yyy=p1*roy3z2+p2*roy3
!c8
      x2yxy=p1*rox3y2+p2*roxy2
      x2zxz=p1*rox3z2+p2*roxz2
      y2xyx=p1*rox2y3+p2*rox2y
      y2zyz=p1*roy3z2+p2*royz2
      z2xzx=p1*rox2z3+p2*rox2z
      z2yzy=p1*roy2z3+p2*roy2z
!c9
      x2yyx=p1*rox3y2+p2*roxy2
      x2zzx=p1*rox3z2+p2*roxz2
      y2xxy=p1*rox2y3+p2*rox2y
      y2zzy=p1*roy3z2+p2*royz2
      z2xxz=p1*rox2z3+p2*rox2z
      z2yyz=p1*roy2z3+p2*roy2z
!c10
      x2yxz=p1*rox3yz+p2*roxyz
      x2zxy=p1*rox3yz+p2*roxyz
      y2xyz=p1*roxy3z+p2*roxyz
      y2zyx=p1*roxy3z+p2*roxyz
      z2xzy=p1*roxyz3+p2*roxyz
      z2yzx=p1*roxyz3+p2*roxyz
!c11
      x2yyz=p1*rox2y2z+p2*roy2z
      x2zzy=p1*rox2yz2+p2*royz2
      y2xxz=p1*rox2y2z+p2*rox2z
      y2zzx=p1*roxy2z2+p2*roxz2
      z2xxy=p1*rox2yz2+p2*rox2y
      z2yyx=p1*roxy2z2+p2*roxy2
!c12
      x2yzx=p1*rox3yz+p2*roxyz
      x2zyx=p1*rox3yz+p2*roxyz
      y2xzy=p1*roxy3z+p2*roxyz
      y2zxy=p1*roxy3z+p2*roxyz
      z2xyz=p1*roxyz3+p2*roxyz
      z2yxz=p1*roxyz3+p2*roxyz
!c13
      x2yzy=p1*rox2y2z+p2*roy2z
      x2zyz=p1*rox2yz2+p2*royz2
      y2xzx=p1*rox2y2z+p2*rox2z
      y2zxz=p1*roxy2z2+p2*roxz2
      z2xyx=p1*rox2yz2+p2*rox2y
      z2yxy=p1*roxy2z2+p2*roxy2
!c14
      x2yzz=p1*rox2yz2+p2*royz2
      x2zyy=p1*rox2y2z+p2*roy2z
      y2xzz=p1*roxy2z2+p2*roxz2
      y2zxx=p1*rox2y2z+p2*rox2z
      z2xyy=p1*roxy2z2+p2*roxy2
      z2yxx=p1*rox2yz2+p2*rox2y
!c15
      xyzxx=p1*rox3yz
      xyzyy=p1*roxy3z
      xyzzz=p1*roxyz3
!c16
      xyzxy=p1*rox2y2z
      xyzxz=p1*rox2yz2
      xyzyx=p1*rox2y2z
      xyzyz=p1*roxy2z2
      xyzzx=p1*rox2yz2
      xyzzy=p1*roxy2z2
! c
      pf(1,1,1)=xyzxx
      pf(1,1,2)=xyzyx
      pf(1,1,3)=xyzzx
      pf(2,1,1)=xyzxy
      pf(2,1,2)=xyzyy
      pf(2,1,3)=xyzzy
      pf(3,1,1)=xyzxz
      pf(3,1,2)=xyzyz
      pf(3,1,3)=xyzzz
      pf(1,2,1)=(x2zxx-y2zxx)
      pf(1,2,2)=(x2zyx-y2zyx)
      pf(1,2,3)=(x2zzx-y2zzx)
      pf(2,2,1)=(x2zxy-y2zxy)
      pf(2,2,2)=(x2zyy-y2zyy)
      pf(2,2,3)=(x2zzy-y2zzy)
      pf(3,2,1)=(x2zxz-y2zxz)
      pf(3,2,2)=(x2zyz-y2zyz)
      pf(3,2,3)=(x2zzz-y2zzz)
      pf(1,3,1)=(x3xx-3.d0*y2xxx)
      pf(1,3,2)=(x3yx-3.d0*y2xyx)
      pf(1,3,3)=(x3zx-3.d0*y2xzx)
      pf(2,3,1)=(x3xy-3.d0*y2xxy)
      pf(2,3,2)=(x3yy-3.d0*y2xyy)
      pf(2,3,3)=(x3zy-3.d0*y2xzy)
      pf(3,3,1)=(x3xz-3.d0*y2xxz)
      pf(3,3,2)=(x3yz-3.d0*y2xyz)
      pf(3,3,3)=(x3zz-3.d0*y2xzz)
      pf(1,4,1)=(3.d0*x2yxx-y3xx)
      pf(1,4,2)=(3.d0*x2yyx-y3yx)
      pf(1,4,3)=(3.d0*x2yzx-y3zx)
      pf(2,4,1)=(3.d0*x2yxy-y3xy)
      pf(2,4,2)=(3.d0*x2yyy-y3yy)
      pf(2,4,3)=(3.d0*x2yzy-y3zy)
      pf(3,4,1)=(3.d0*x2yxz-y3xz)
      pf(3,4,2)=(3.d0*x2yyz-y3yz)
      pf(3,4,3)=(3.d0*x2yzz-y3zz)
      pf(1,5,1)=2.d0*z3xx-3.d0*x2zxx-3.d0*y2zxx
      pf(1,5,2)=2.d0*z3yx-3.d0*x2zyx-3.d0*y2zyx
      pf(1,5,3)=2.d0*z3zx-3.d0*x2zzx-3.d0*y2zzx
      pf(2,5,1)=2.d0*z3xy-3.d0*x2zxy-3.d0*y2zxy
      pf(2,5,2)=2.d0*z3yy-3.d0*x2zyy-3.d0*y2zyy
      pf(2,5,3)=2.d0*z3zy-3.d0*x2zzy-3.d0*y2zzy
      pf(3,5,1)=2.d0*z3xz-3.d0*x2zxz-3.d0*y2zxz
      pf(3,5,2)=2.d0*z3yz-3.d0*x2zyz-3.d0*y2zyz
      pf(3,5,3)=2.d0*z3zz-3.d0*x2zzz-3.d0*y2zzz
      pf(1,6,1)=(4.d0*z2xxx-x3xx-y2xxx)
      pf(1,6,2)=(4.d0*z2xyx-x3yx-y2xyx)
      pf(1,6,3)=(4.d0*z2xzx-x3zx-y2xzx)
      pf(2,6,1)=(4.d0*z2xxy-x3xy-y2xxy)
      pf(2,6,2)=(4.d0*z2xyy-x3yy-y2xyy)
      pf(2,6,3)=(4.d0*z2xzy-x3zy-y2xzy)
      pf(3,6,1)=(4.d0*z2xxz-x3xz-y2xxz)
      pf(3,6,2)=(4.d0*z2xyz-x3yz-y2xyz)
      pf(3,6,3)=(4.d0*z2xzz-x3zz-y2xzz)
      pf(1,7,1)=(4.d0*z2yxx-x2yxx-y3xx)
      pf(1,7,2)=(4.d0*z2yyx-x2yyx-y3yx)
      pf(1,7,3)=(4.d0*z2yzx-x2yzx-y3zx)
      pf(2,7,1)=(4.d0*z2yxy-x2yxy-y3xy)
      pf(2,7,2)=(4.d0*z2yyy-x2yyy-y3yy)
      pf(2,7,3)=(4.d0*z2yzy-x2yzy-y3zy)
      pf(3,7,1)=(4.d0*z2yxz-x2yxz-y3xz)
      pf(3,7,2)=(4.d0*z2yyz-x2yyz-y3yz)
      pf(3,7,3)=(4.d0*z2yzz-x2yzz-y3zz)
! c
      do 1055 k=1,7
      do 1055 m=1,3
      do 1055 l=1,3
      q(1+m,9+k,l)=pf(m,k,l)
 1055 continue
 1060 continue
! c
      if(nx1.lt.16.or.nx2.lt.9) go to 1070
      d1=delt/(32.d0*a13*a22)
      d2=delt/(16.d0*a12*a22)
      d3=3.d0*delt/(16.d0*a12*a22)
      d4=delt/(16.d0*a13*a2)
      d6=delt/(8.d0*a12*a2)
      d7=3.d0*delt/(8.d0*a12*a2)
!c1
      x3xx2=d1*rox6+(d3+d4)*rox4+d7*rox2
      y3yy2=d1*roy6+(d3+d4)*roy4+d7*roy2
      z3zz2=d1*roz6+(d3+d4)*roz4+d7*roz2
!c2
      x3xy2=d1*rox4y2+d3*rox2y2+d4*rox4+d7*rox2
      x3xz2=d1*rox4z2+d3*rox2z2+d4*rox4+d7*rox2
      y3yx2=d1*rox2y4+d3*rox2y2+d4*roy4+d7*roy2
      y3yz2=d1*roy4z2+d3*roy2z2+d4*roy4+d7*roy2
      z3zx2=d1*rox2z4+d3*rox2z2+d4*roz4+d7*roz2
      z3zy2=d1*roy2z4+d3*roy2z2+d4*roz4+d7*roz2
!c3
      x3xxy=d1*rox5y+d3*rox3y
      x3xxz=d1*rox5z+d3*rox3z
      y3yxy=d1*roxy5+d3*roxy3
      y3yyz=d1*roy5z+d3*roy3z
      z3zxz=d1*roxz5+d3*roxz3
      z3zyz=d1*royz5+d3*royz3
!c4
      x3yx2=d1*rox5y+(d3+d4)*rox3y+d7*roxy
      x3zx2=d1*rox5z+(d3+d4)*rox3z+d7*roxz
      y3xy2=d1*roxy5+(d3+d4)*roxy3+d7*roxy
      y3zy2=d1*roy5z+(d3+d4)*roy3z+d7*royz
      z3xz2=d1*roxz5+(d3+d4)*roxz3+d7*roxz
      z3yz2=d1*royz5+(d3+d4)*royz3+d7*royz
!c5
      x3yy2=d1*rox3y3+d3*roxy3+d4*rox3y+d7*roxy
      x3zz2=d1*rox3z3+d3*roxz3+d4*rox3z+d7*roxz
      y3xx2=d1*rox3y3+d3*rox3y+d4*roxy3+d7*roxy
      y3zz2=d1*roy3z3+d3*royz3+d4*roy3z+d7*royz
      z3xx2=d1*rox3z3+d3*rox3z+d4*roxz3+d7*roxz
      z3yy2=d1*roy3z3+d3*roy3z+d4*royz3+d7*royz
!c6
      x3yz2=d1*rox3yz2+d3*roxyz2+d4*rox3y+d7*roxy
      y3zx2=d1*rox2y3z+d3*rox2yz+d4*roy3z+d7*royz
      z3xy2=d1*roxy2z3+d3*roxy2z+d4*roxz3+d7*roxz
      x3zy2=d1*rox3y2z+d3*roxy2z+d4*rox3z+d7*roxz
      y3xz2=d1*roxy3z2+d3*roxyz2+d4*roxy3+d7*roxy
      z3yx2=d1*rox2yz3+d3*rox2yz+d4*royz3+d7*royz
!c7
      x3yxy=d1*rox4y2+d3*rox2y2
      x3zxz=d1*rox4z2+d3*rox2z2
      y3xxy=d1*rox2y4+d3*rox2y2
      y3zyz=d1*roy4z2+d3*roy2z2
      z3xxz=d1*rox2z4+d3*rox2z2
      z3yyz=d1*roy2z4+d3*roy2z2
!c8
      x3zxy=d1*rox4yz+d3*rox2yz
      x3yxz=d1*rox4yz+d3*rox2yz
      y3xyz=d1*roxy4z+d3*roxy2z
      y3zxy=d1*roxy4z+d3*roxy2z
      z3xyz=d1*roxyz4+d3*roxyz2
      z3yxz=d1*roxyz4+d3*roxyz2
!c9
      x3xyz=d1*rox4yz+d3*rox2yz
      y3yxz=d1*roxy4z+d3*roxy2z
      z3zxy=d1*roxyz4+d3*roxyz2
!c10
      x2yxx2=d1*rox5y+(d2+d4)*rox3y+d6*roxy
      y2xyy2=d1*roxy5+(d2+d4)*roxy3+d6*roxy
      z2xzz2=d1*roxz5+(d2+d4)*roxz3+d6*roxz
      x2zxx2=d1*rox5z+(d2+d4)*rox3z+d6*roxz
      y2zyy2=d1*roy5z+(d2+d4)*roy3z+d6*royz
      z2yzz2=d1*royz5+(d2+d4)*royz3+d6*royz
!c11
      x2yxy2=d1*rox3y3+d2*roxy3+d4*rox3y+d6*roxy
      x2zxz2=d1*rox3z3+d2*roxz3+d4*rox3z+d6*roxz
      y2xyx2=d1*rox3y3+d2*rox3y+d4*roxy3+d6*roxy
      y2zyz2=d1*roy3z3+d2*royz3+d4*roy3z+d6*royz
      z2xzx2=d1*rox3z3+d2*rox3z+d4*roxz3+d6*roxz
      z2yzy2=d1*roy3z3+d2*roy3z+d4*royz3+d6*royz
!c12
      x2yxz2=d1*rox3yz2+d2*roxyz2+d4*rox3y+d6*roxy
      x2zxy2=d1*rox3y2z+d2*roxy2z+d4*rox3z+d6*roxz
      y2xyz2=d1*roxy3z2+d2*roxyz2+d4*roxy3+d6*roxy
      y2zyx2=d1*rox2y3z+d2*rox2yz+d4*roy3z+d6*royz
      z2xzy2=d1*roxy2z3+d2*roxy2z+d4*roxz3+d6*roxz
      z2yzx2=d1*rox2yz3+d2*rox2yz+d4*royz3+d6*royz
!c13
      x2yxxy=d1*rox4y2+d2*rox2y2
      y2zyyz=d1*roy4z2+d2*roy2z2
      z2xzxz=d1*rox2z4+d2*rox2z2
      x2zxxz=d1*rox4z2+d2*rox2z2
      y2xyxy=d1*rox2y4+d2*rox2y2
      z2yzyz=d1*roy2z4+d2*roy2z2
!c14
      x2yxxz=d1*rox4yz+d2*rox2yz
      x2zxxy=d1*rox4yz+d2*rox2yz
      y2xyyz=d1*roxy4z+d2*roxy2z
      y2zyxy=d1*roxy4z+d2*roxy2z
      z2xzyz=d1*roxyz4+d2*roxyz2
      z2yzxz=d1*roxyz4+d2*roxyz2
!c15
      x2yxyz=d1*rox3y2z+d2*roxy2z
      x2zxyz=d1*rox3yz2+d2*roxyz2
      y2xyxz=d1*rox2y3z+d2*rox2yz
      y2zyxz=d1*roxy3z2+d2*roxyz2
      z2xzxy=d1*rox2yz3+d2*rox2yz
      z2yzxy=d1*roxy2z3+d2*roxy2z
!c16
      x2yyx2=d1*rox4y2+(d2+d4)*rox2y2+d6*roy2
      x2zzx2=d1*rox4z2+(d2+d4)*rox2z2+d6*roz2
      y2xxy2=d1*rox2y4+(d2+d4)*rox2y2+d6*rox2
      y2zzy2=d1*roy4z2+(d2+d4)*roy2z2+d6*roz2
      z2xxz2=d1*rox2z4+(d2+d4)*rox2z2+d6*rox2
      z2yyz2=d1*roy2z4+(d2+d4)*roy2z2+d6*roy2
!c17
      x2yyy2=d1*rox2y4+d2*roy4+d4*rox2y2+d6*roy2
      x2zzz2=d1*rox2z4+d2*roz4+d4*rox2z2+d6*roz2
      y2xxx2=d1*rox4y2+d2*rox4+d4*rox2y2+d6*rox2
      y2zzz2=d1*roy2z4+d2*roz4+d4*roy2z2+d6*roz2
      z2xxx2=d1*rox4z2+d2*rox4+d4*rox2z2+d6*rox2
      z2yyy2=d1*roy4z2+d2*roy4+d4*roy2z2+d6*roy2
!c18
      x2yyz2=d1*rox2y2z2+d2*roy2z2+d4*rox2y2+d6*roy2
      x2zzy2=d1*rox2y2z2+d2*roy2z2+d4*rox2z2+d6*roz2
      y2xxz2=d1*rox2y2z2+d2*rox2z2+d4*rox2y2+d6*rox2
      y2zzx2=d1*rox2y2z2+d2*rox2z2+d4*roy2z2+d6*roz2
      z2xxy2=d1*rox2y2z2+d2*rox2y2+d4*rox2z2+d6*rox2
      z2yyx2=d1*rox2y2z2+d2*rox2y2+d4*roy2z2+d6*roy2
!c19
      x2yyxy=d1*rox3y3+d2*roxy3
      x2zzxz=d1*rox3z3+d2*roxz3
      y2xxxy=d1*rox3y3+d2*rox3y
      y2zzyz=d1*roy3z3+d2*royz3
      z2xxxz=d1*rox3z3+d2*rox3z
      z2yyyz=d1*roy3z3+d2*roy3z
!c20
      x2yyxz=d1*rox3y2z+d2*roxy2z
      x2zzxy=d1*rox3yz2+d2*roxyz2
      y2xxyz=d1*rox2y3z+d2*rox2yz
      y2zzxy=d1*roxy3z2+d2*roxyz2
      z2xxyz=d1*rox2yz3+d2*rox2yz
      z2yyxz=d1*roxy2z3+d2*roxy2z
!c21
      x2yyyz=d1*rox2y3z+d2*roy3z
      x2zzyz=d1*rox2yz3+d2*royz3
      y2xxxz=d1*rox3y2z+d2*rox3z
      y2zzxz=d1*roxy2z3+d2*roxz3
      z2xxxy=d1*rox3yz2+d2*rox3y
      z2yyxy=d1*roxy3z2+d2*roxy3
!c22
      x2yzx2=d1*rox4yz+(d2+d4)*rox2yz+d6*royz
      x2zyx2=d1*rox4yz+(d2+d4)*rox2yz+d6*royz
      y2xzy2=d1*roxy4z+(d2+d4)*roxy2z+d6*roxz
      y2zxy2=d1*roxy4z+(d2+d4)*roxy2z+d6*roxz
      z2xyz2=d1*roxyz4+(d2+d4)*roxyz2+d6*roxy
      z2yxz2=d1*roxyz4+(d2+d4)*roxyz2+d6*roxy
!c23
      x2yzy2=d1*rox2y3z+d2*roy3z+d4*rox2yz+d6*royz
      x2zyz2=d1*rox2yz3+d2*royz3+d4*rox2yz+d6*royz
      y2xzx2=d1*rox3y2z+d2*rox3z+d4*roxy2z+d6*roxz
      y2zxz2=d1*roxy2z3+d2*roxz3+d4*roxy2z+d6*roxz
      z2xyx2=d1*rox3yz2+d2*rox3y+d4*roxyz2+d6*roxy
      z2yxy2=d1*roxy3z2+d2*roxy3+d4*roxyz2+d6*roxy
!c24
      x2yzz2=d1*rox2yz3+d2*royz3+d4*rox2yz+d6*royz
      x2zyy2=d1*rox2y3z+d2*roy3z+d4*rox2yz+d6*royz
      y2xzz2=d1*roxy2z3+d2*roxz3+d4*roxy2z+d6*roxz
      y2zxx2=d1*rox3y2z+d2*rox3z+d4*roxy2z+d6*roxz
      z2xyy2=d1*roxy3z2+d2*roxy3+d4*roxyz2+d6*roxy
      z2yxx2=d1*rox3yz2+d2*rox3y+d4*roxyz2+d6*roxy
!c25
      x2yzxy=d1*rox3y2z+d2*roxy2z
      x2zyxz=d1*rox3yz2+d2*roxyz2
      y2xzxy=d1*rox2y3z+d2*rox2yz
      y2zxyz=d1*roxy3z2+d2*roxyz2
      z2xyxz=d1*rox2yz3+d2*rox2yz
      z2yxyz=d1*roxy2z3+d2*roxy2z
!c26
      x2yzxz=d1*rox3yz2+d2*roxyz2
      x2zyxy=d1*rox3y2z+d2*roxy2z
      y2xzyz=d1*roxy3z2+d2*roxyz2
      y2zxxy=d1*rox2y3z+d2*rox2yz
      z2xyyz=d1*roxy2z3+d2*roxy2z
      z2yxxz=d1*rox2yz3+d2*rox2yz
!c27
      x2yzyz=d1*rox2y2z2+d2*roy2z2
      x2zyyz=d1*rox2y2z2+d2*roy2z2
      y2xzxz=d1*rox2y2z2+d2*rox2z2
      y2zxxz=d1*rox2y2z2+d2*rox2z2
      z2xyxy=d1*rox2y2z2+d2*rox2y2
      z2yxxy=d1*rox2y2z2+d2*rox2y2
!c28
      xyzxx2=d1*rox4yz+d4*rox2yz
      xyzyy2=d1*roxy4z+d4*roxy2z
      xyzzz2=d1*roxyz4+d4*roxyz2
!c29
      xyzxy2=d1*rox2y3z+d4*rox2yz
      xyzxz2=d1*rox2yz3+d4*rox2yz
      xyzyx2=d1*rox3y2z+d4*roxy2z
      xyzyz2=d1*roxy2z3+d4*roxy2z
      xyzzx2=d1*rox3yz2+d4*roxyz2
      xyzzy2=d1*roxy3z2+d4*roxyz2
!c30
      xyzxxy=d1*rox3y2z
      xyzxxz=d1*rox3yz2
      xyzyxy=d1*rox2y3z
      xyzyyz=d1*roxy3z2
      xyzzxz=d1*rox2yz3
      xyzzyz=d1*roxy2z3
!c31
      xyzxyz=d1*rox2y2z2
      xyzyxz=d1*rox2y2z2
      xyzzxy=d1*rox2y2z2
!c32
      x3yyz=d1*rox3y2z+d3*roxy2z
      y3xxz=d1*rox2y3z+d3*rox2yz
      z3xxy=d1*rox2yz3+d3*rox2yz
      x3zyz=d1*rox3yz2+d3*roxyz2
      y3zxz=d1*roxy3z2+d3*roxyz2
      z3yxy=d1*roxy2z3+d3*roxy2z
! c
      fd(1,1,1)=xyzxxy
      fd(1,1,2)=xyzyxy
      fd(1,1,3)=xyzzxy
      fd(1,2,1)=xyzxxz
      fd(1,2,2)=xyzyxz
      fd(1,2,3)=xyzzxz
      fd(1,3,1)=xyzxyz
      fd(1,3,2)=xyzyyz
      fd(1,3,3)=xyzzyz
      fd(1,4,1)=(xyzxx2-xyzxy2)
      fd(1,4,2)=(xyzyx2-xyzyy2)
      fd(1,4,3)=(xyzzx2-xyzzy2)
      fd(1,5,1)=(2.d0*xyzxz2-xyzxx2-xyzxy2)
      fd(1,5,2)=(2.d0*xyzyz2-xyzyx2-xyzyy2)
      fd(1,5,3)=(2.d0*xyzzz2-xyzzx2-xyzzy2)
      fd(2,1,1)=(x2zxxy-y2zxxy)
      fd(2,1,2)=(x2zyxy-y2zyxy)
      fd(2,1,3)=(x2zzxy-y2zzxy)
      fd(2,2,1)=(x2zxxz-y2zxxz)
      fd(2,2,2)=(x2zyxz-y2zyxz)
      fd(2,2,3)=(x2zzxz-y2zzxz)
      fd(2,3,1)=(x2zxyz-y2zxyz)
      fd(2,3,2)=(x2zyyz-y2zyyz)
      fd(2,3,3)=(x2zzyz-y2zzyz)
      fd(2,4,1)=(x2zxx2-y2zxx2-x2zxy2+y2zxy2)
      fd(2,4,2)=(x2zyx2-y2zyx2-x2zyy2+y2zyy2)
      fd(2,4,3)=(x2zzx2-y2zzx2-x2zzy2+y2zzy2)
      fd(2,5,1)=(2.d0*x2zxz2-2.d0*y2zxz2-x2zxx2+y2zxx2 &
     & -x2zxy2+y2zxy2)
      fd(2,5,2)=(2.d0*x2zyz2-2.d0*y2zyz2-x2zyx2+y2zyx2 &
     & -x2zyy2+y2zyy2)
      fd(2,5,3)=(2.d0*x2zzz2-2.d0*y2zzz2-x2zzx2+y2zzx2 &
     & -x2zzy2+y2zzy2)
      fd(3,1,1)=(x3xxy-3.d0*y2xxxy)
      fd(3,1,2)=(x3yxy-3.d0*y2xyxy)
      fd(3,1,3)=(x3zxy-3.d0*y2xzxy)
      fd(3,2,1)=(x3xxz-3.d0*y2xxxz)
      fd(3,2,2)=(x3yxz-3.d0*y2xyxz)
      fd(3,2,3)=(x3zxz-3.d0*y2xzxz)
      fd(3,3,1)=(x3xyz-3.d0*y2xxyz)
      fd(3,3,2)=(x3yyz-3.d0*y2xyyz)
      fd(3,3,3)=(x3zyz-3.d0*y2xzyz)
      fd(3,4,1)=(x3xx2-3.d0*y2xxx2-x3xy2+3.d0*y2xxy2)
      fd(3,4,2)=(x3yx2-3.d0*y2xyx2-x3yy2+3.d0*y2xyy2)
      fd(3,4,3)=(x3zx2-3.d0*y2xzx2-x3zy2+3.d0*y2xzy2)
      fd(3,5,1)=(2.d0*x3xz2-6.d0*y2xxz2-x3xx2 &
     & +3.d0*y2xxx2-x3xy2+3.d0*y2xxy2)
      fd(3,5,2)=(2.d0*x3yz2-6.d0*y2xyz2-x3yx2 &
     & +3.d0*y2xyx2-x3yy2+3.d0*y2xyy2)
      fd(3,5,3)=(2.d0*x3zz2-6.d0*y2xzz2-x3zx2 &
     & +3.d0*y2xzx2-x3zy2+3.d0*y2xzy2)
      fd(4,1,1)=(3.d0*x2yxxy-y3xxy)
      fd(4,1,2)=(3.d0*x2yyxy-y3yxy)
      fd(4,1,3)=(3.d0*x2yzxy-y3zxy)
      fd(4,2,1)=(3.d0*x2yxxz-y3xxz)
      fd(4,2,2)=(3.d0*x2yyxz-y3yxz)
      fd(4,2,3)=(3.d0*x2yzxz-y3zxz)
      fd(4,3,1)=(3.d0*x2yxyz-y3xyz)
      fd(4,3,2)=(3.d0*x2yyyz-y3yyz)
      fd(4,3,3)=(3.d0*x2yzyz-y3zyz)
      fd(4,4,1)=(3.d0*x2yxx2-y3xx2-3.d0*x2yxy2+y3xy2)
      fd(4,4,2)=(3.d0*x2yyx2-y3yx2-3.d0*x2yyy2+y3yy2)
      fd(4,4,3)=(3.d0*x2yzx2-y3zx2-3.d0*x2yzy2+y3zy2)
      fd(4,5,1)=(6.d0*x2yxz2-2.d0*y3xz2-3.d0*x2yxx2 &
     & +y3xx2-3.d0*x2yxy2+y3xy2)
      fd(4,5,2)=(6.d0*x2yyz2-2.d0*y3yz2-3.d0*x2yyx2 &
     & +y3yx2-3.d0*x2yyy2+y3yy2)
      fd(4,5,3)=(6.d0*x2yzz2-2.d0*y3zz2-3.d0*x2yzx2 &
     & +y3zx2-3.d0*x2yzy2+y3zy2)
      fd(5,1,1)=2.d0*z3xxy-3.d0*x2zxxy-3.d0*y2zxxy
      fd(5,1,2)=2.d0*z3yxy-3.d0*x2zyxy-3.d0*y2zyxy
      fd(5,1,3)=2.d0*z3zxy-3.d0*x2zzxy-3.d0*y2zzxy
      fd(5,2,1)=2.d0*z3xxz-3.d0*x2zxxz-3.d0*y2zxxz
      fd(5,2,2)=2.d0*z3yxz-3.d0*x2zyxz-3.d0*y2zyxz
      fd(5,2,3)=2.d0*z3zxz-3.d0*x2zzxz-3.d0*y2zzxz
      fd(5,3,1)=2.d0*z3xyz-3.d0*x2zxyz-3.d0*y2zxyz
      fd(5,3,2)=2.d0*z3yyz-3.d0*x2zyyz-3.d0*y2zyyz
      fd(5,3,3)=2.d0*z3zyz-3.d0*x2zzyz-3.d0*y2zzyz
      fd(5,4,1)=z3xx2-1.5d0*x2zxx2-1.5d0*y2zxx2 &
     & -z3xy2+1.5d0*x2zxy2+1.5d0*y2zxy2
      fd(5,4,2)=z3yx2-1.5d0*x2zyx2-1.5d0*y2zyx2 &
     & -z3yy2+1.5d0*x2zyy2+1.5d0*y2zyy2
      fd(5,4,3)=z3zx2-1.5d0*x2zzx2-1.5d0*y2zzx2 &
     & -z3zy2+1.5d0*x2zzy2+1.5d0*y2zzy2
      fd(5,5,1)=(4.d0*z3xz2-6.d0*x2zxz2-6.d0*y2zxz2 &
     & -2.d0*z3xx2+3.d0*x2zxx2+3.d0*y2zxx2 &
     & -2.d0*z3xy2+3.d0*x2zxy2+3.d0*y2zxy2)
      fd(5,5,2)=(4.d0*z3yz2-6.d0*x2zyz2-6.d0*y2zyz2 &
     & -2.d0*z3yx2+3.d0*x2zyx2+3.d0*y2zyx2 &
     & -2.d0*z3yy2+3.d0*x2zyy2+3.d0*y2zyy2)
      fd(5,5,3)=(4.d0*z3zz2-6.d0*x2zzz2-6.d0*y2zzz2 &
     & -2.d0*z3zx2+3.d0*x2zzx2+3.d0*y2zzx2 &
     & -2.d0*z3zy2+3.d0*x2zzy2+3.d0*y2zzy2)
      fd(6,1,1)=(4.d0*z2xxxy-x3xxy-y2xxxy)
      fd(6,1,2)=(4.d0*z2xyxy-x3yxy-y2xyxy)
      fd(6,1,3)=(4.d0*z2xzxy-x3zxy-y2xzxy)
      fd(6,2,1)=(4.d0*z2xxxz-x3xxz-y2xxxz)
      fd(6,2,2)=(4.d0*z2xyxz-x3yxz-y2xyxz)
      fd(6,2,3)=(4.d0*z2xzxz-x3zxz-y2xzxz)
      fd(6,3,1)=(4.d0*z2xxyz-x3xyz-y2xxyz)
      fd(6,3,2)=(4.d0*z2xyyz-x3yyz-y2xyyz)
      fd(6,3,3)=(4.d0*z2xzyz-x3zyz-y2xzyz)
      fd(6,4,1)=(4.d0*z2xxx2-x3xx2-y2xxx2 &
     & -4.d0*z2xxy2+x3xy2+y2xxy2)
      fd(6,4,2)=(4.d0*z2xyx2-x3yx2-y2xyx2 &
     & -4.d0*z2xyy2+x3yy2+y2xyy2)
      fd(6,4,3)=(4.d0*z2xzx2-x3zx2-y2xzx2 &
     & -4.d0*z2xzy2+x3zy2+y2xzy2)
      fd(6,5,1)=(8.d0*z2xxz2-2.d0*x3xz2-2.d0*y2xxz2 &
     & -4.d0*z2xxx2+x3xx2+y2xxx2 &
     & -4.d0*z2xxy2+x3xy2+y2xxy2)
      fd(6,5,2)=(8.d0*z2xyz2-2.d0*x3yz2-2.d0*y2xyz2 &
     & -4.d0*z2xyx2+x3yx2+y2xyx2 &
     & -4.d0*z2xyy2+x3yy2+y2xyy2)
      fd(6,5,3)=(8.d0*z2xzz2-2.d0*x3zz2-2.d0*y2xzz2 &
     & -4.d0*z2xzx2+x3zx2+y2xzx2 &
     & -4.d0*z2xzy2+x3zy2+y2xzy2)
      fd(7,1,1)=(4.d0*z2yxxy-x2yxxy-y3xxy)
      fd(7,1,2)=(4.d0*z2yyxy-x2yyxy-y3yxy)
      fd(7,1,3)=(4.d0*z2yzxy-x2yzxy-y3zxy)
      fd(7,2,1)=(4.d0*z2yxxz-x2yxxz-y3xxz)
      fd(7,2,2)=(4.d0*z2yyxz-x2yyxz-y3yxz)
      fd(7,2,3)=(4.d0*z2yzxz-x2yzxz-y3zxz)
      fd(7,3,1)=(4.d0*z2yxyz-x2yxyz-y3xyz)
      fd(7,3,2)=(4.d0*z2yyyz-x2yyyz-y3yyz)
      fd(7,3,3)=(4.d0*z2yzyz-x2yzyz-y3zyz)
      fd(7,4,1)=(4.d0*z2yxx2-x2yxx2-y3xx2 &
     & -4.d0*z2yxy2+x2yxy2+y3xy2)
      fd(7,4,2)=(4.d0*z2yyx2-x2yyx2-y3yx2 &
     & -4.d0*z2yyy2+x2yyy2+y3yy2)
      fd(7,4,3)=(4.d0*z2yzx2-x2yzx2-y3zx2 &
     & -4.d0*z2yzy2+x2yzy2+y3zy2)
      fd(7,5,1)=(8.d0*z2yxz2-2.d0*x2yxz2-2.d0*y3xz2 &
     & -4.d0*z2yxx2+x2yxx2+y3xx2-4.d0*z2yxy2+x2yxy2+y3xy2)
      fd(7,5,2)=(8.d0*z2yyz2-2.d0*x2yyz2-2.d0*y3yz2 &
     & -4.d0*z2yyx2+x2yyx2+y3yx2-4.d0*z2yyy2+x2yyy2+y3yy2)
      fd(7,5,3)=(8.d0*z2yzz2-2.d0*x2yzz2-2.d0*y3zz2 &
     & -4.d0*z2yzx2+x2yzx2+y3zx2-4.d0*z2yzy2+x2yzy2+y3zy2)
! c
      do 1065 k=1,7
      do 1065 m=1,5
      do 1065 l=1,3
      q(9+k,4+m,l)=fd(k,m,l)
 1065 continue
 1070 continue
! c
      if(nx1.lt.9.or.nx2.lt.16) go to 1080
      d1=-delt/(32.d0*a23*a12)
      d2=-delt/(16.d0*a12*a22)
      d3=-3.d0*delt/(16.d0*a12*a22)
      d4=-delt/(16.d0*a23*a1)
      d6=-delt/(8.d0*a22*a1)
      d7=-3.d0*delt/(8.d0*a22*a1)
!c1
      x3xx2=d1*rox6+(d3+d4)*rox4+d7*rox2
      y3yy2=d1*roy6+(d3+d4)*roy4+d7*roy2
      z3zz2=d1*roz6+(d3+d4)*roz4+d7*roz2
!c2
      x3xy2=d1*rox4y2+d3*rox2y2+d4*rox4+d7*rox2
      x3xz2=d1*rox4z2+d3*rox2z2+d4*rox4+d7*rox2
      y3yx2=d1*rox2y4+d3*rox2y2+d4*roy4+d7*roy2
      y3yz2=d1*roy4z2+d3*roy2z2+d4*roy4+d7*roy2
      z3zx2=d1*rox2z4+d3*rox2z2+d4*roz4+d7*roz2
      z3zy2=d1*roy2z4+d3*roy2z2+d4*roz4+d7*roz2
!c3
      x3xxy=d1*rox5y+d3*rox3y
      x3xxz=d1*rox5z+d3*rox3z
      y3yxy=d1*roxy5+d3*roxy3
      y3yyz=d1*roy5z+d3*roy3z
      z3zxz=d1*roxz5+d3*roxz3
      z3zyz=d1*royz5+d3*royz3
!c4
      x3yx2=d1*rox5y+(d3+d4)*rox3y+d7*roxy
      x3zx2=d1*rox5z+(d3+d4)*rox3z+d7*roxz
      y3xy2=d1*roxy5+(d3+d4)*roxy3+d7*roxy
      y3zy2=d1*roy5z+(d3+d4)*roy3z+d7*royz
      z3xz2=d1*roxz5+(d3+d4)*roxz3+d7*roxz
      z3yz2=d1*royz5+(d3+d4)*royz3+d7*royz
!c5
      x3yy2=d1*rox3y3+d3*roxy3+d4*rox3y+d7*roxy
      x3zz2=d1*rox3z3+d3*roxz3+d4*rox3z+d7*roxz
      y3xx2=d1*rox3y3+d3*rox3y+d4*roxy3+d7*roxy
      y3zz2=d1*roy3z3+d3*royz3+d4*roy3z+d7*royz
      z3xx2=d1*rox3z3+d3*rox3z+d4*roxz3+d7*roxz
      z3yy2=d1*roy3z3+d3*roy3z+d4*royz3+d7*royz
!c6
      x3yz2=d1*rox3yz2+d3*roxyz2+d4*rox3y+d7*roxy
      y3zx2=d1*rox2y3z+d3*rox2yz+d4*roy3z+d7*royz
      z3xy2=d1*roxy2z3+d3*roxy2z+d4*roxz3+d7*roxz
      x3zy2=d1*rox3y2z+d3*roxy2z+d4*rox3z+d7*roxz
      y3xz2=d1*roxy3z2+d3*roxyz2+d4*roxy3+d7*roxy
      z3yx2=d1*rox2yz3+d3*rox2yz+d4*royz3+d7*royz
!c7
      x3yxy=d1*rox4y2+d3*rox2y2
      x3zxz=d1*rox4z2+d3*rox2z2
      y3xxy=d1*rox2y4+d3*rox2y2
      y3zyz=d1*roy4z2+d3*roy2z2
      z3xxz=d1*rox2z4+d3*rox2z2
      z3yyz=d1*roy2z4+d3*roy2z2
!c8
      x3zxy=d1*rox4yz+d3*rox2yz
      x3yxz=d1*rox4yz+d3*rox2yz
      y3xyz=d1*roxy4z+d3*roxy2z
      y3zxy=d1*roxy4z+d3*roxy2z
      z3xyz=d1*roxyz4+d3*roxyz2
      z3yxz=d1*roxyz4+d3*roxyz2
!c9
      x3xyz=d1*rox4yz+d3*rox2yz
      y3yxz=d1*roxy4z+d3*roxy2z
      z3zxy=d1*roxyz4+d3*roxyz2
!c10
      x2yxx2=d1*rox5y+(d2+d4)*rox3y+d6*roxy
      y2xyy2=d1*roxy5+(d2+d4)*roxy3+d6*roxy
      z2xzz2=d1*roxz5+(d2+d4)*roxz3+d6*roxz
      x2zxx2=d1*rox5z+(d2+d4)*rox3z+d6*roxz
      y2zyy2=d1*roy5z+(d2+d4)*roy3z+d6*royz
      z2yzz2=d1*royz5+(d2+d4)*royz3+d6*royz
!c11
      x2yxy2=d1*rox3y3+d2*roxy3+d4*rox3y+d6*roxy
      x2zxz2=d1*rox3z3+d2*roxz3+d4*rox3z+d6*roxz
      y2xyx2=d1*rox3y3+d2*rox3y+d4*roxy3+d6*roxy
      y2zyz2=d1*roy3z3+d2*royz3+d4*roy3z+d6*royz
      z2xzx2=d1*rox3z3+d2*rox3z+d4*roxz3+d6*roxz
      z2yzy2=d1*roy3z3+d2*roy3z+d4*royz3+d6*royz
!c12
      x2yxz2=d1*rox3yz2+d2*roxyz2+d4*rox3y+d6*roxy
      x2zxy2=d1*rox3y2z+d2*roxy2z+d4*rox3z+d6*roxz
      y2xyz2=d1*roxy3z2+d2*roxyz2+d4*roxy3+d6*roxy
      y2zyx2=d1*rox2y3z+d2*rox2yz+d4*roy3z+d6*royz
      z2xzy2=d1*roxy2z3+d2*roxy2z+d4*roxz3+d6*roxz
      z2yzx2=d1*rox2yz3+d2*rox2yz+d4*royz3+d6*royz
!c13
      x2yxxy=d1*rox4y2+d2*rox2y2
      y2zyyz=d1*roy4z2+d2*roy2z2
      z2xzxz=d1*rox2z4+d2*rox2z2
      x2zxxz=d1*rox4z2+d2*rox2z2
      y2xyxy=d1*rox2y4+d2*rox2y2
      z2yzyz=d1*roy2z4+d2*roy2z2
!c14
      x2yxxz=d1*rox4yz+d2*rox2yz
      x2zxxy=d1*rox4yz+d2*rox2yz
      y2xyyz=d1*roxy4z+d2*roxy2z
      y2zyxy=d1*roxy4z+d2*roxy2z
      z2xzyz=d1*roxyz4+d2*roxyz2
      z2yzxz=d1*roxyz4+d2*roxyz2
!c15
      x2yxyz=d1*rox3y2z+d2*roxy2z
      x2zxyz=d1*rox3yz2+d2*roxyz2
      y2xyxz=d1*rox2y3z+d2*rox2yz
      y2zyxz=d1*roxy3z2+d2*roxyz2
      z2xzxy=d1*rox2yz3+d2*rox2yz
      z2yzxy=d1*roxy2z3+d2*roxy2z
!c16
      x2yyx2=d1*rox4y2+(d2+d4)*rox2y2+d6*roy2
      x2zzx2=d1*rox4z2+(d2+d4)*rox2z2+d6*roz2
      y2xxy2=d1*rox2y4+(d2+d4)*rox2y2+d6*rox2
      y2zzy2=d1*roy4z2+(d2+d4)*roy2z2+d6*roz2
      z2xxz2=d1*rox2z4+(d2+d4)*rox2z2+d6*rox2
      z2yyz2=d1*roy2z4+(d2+d4)*roy2z2+d6*roy2
!c17
      x2yyy2=d1*rox2y4+d2*roy4+d4*rox2y2+d6*roy2
      x2zzz2=d1*rox2z4+d2*roz4+d4*rox2z2+d6*roz2
      y2xxx2=d1*rox4y2+d2*rox4+d4*rox2y2+d6*rox2
      y2zzz2=d1*roy2z4+d2*roz4+d4*roy2z2+d6*roz2
      z2xxx2=d1*rox4z2+d2*rox4+d4*rox2z2+d6*rox2
      z2yyy2=d1*roy4z2+d2*roy4+d4*roy2z2+d6*roy2
!c18
      x2yyz2=d1*rox2y2z2+d2*roy2z2+d4*rox2y2+d6*roy2
      x2zzy2=d1*rox2y2z2+d2*roy2z2+d4*rox2z2+d6*roz2
      y2xxz2=d1*rox2y2z2+d2*rox2z2+d4*rox2y2+d6*rox2
      y2zzx2=d1*rox2y2z2+d2*rox2z2+d4*roy2z2+d6*roz2
      z2xxy2=d1*rox2y2z2+d2*rox2y2+d4*rox2z2+d6*rox2
      z2yyx2=d1*rox2y2z2+d2*rox2y2+d4*roy2z2+d6*roy2
!c19
      x2yyxy=d1*rox3y3+d2*roxy3
      x2zzxz=d1*rox3z3+d2*roxz3
      y2xxxy=d1*rox3y3+d2*rox3y
      y2zzyz=d1*roy3z3+d2*royz3
      z2xxxz=d1*rox3z3+d2*rox3z
      z2yyyz=d1*roy3z3+d2*roy3z
!c20
      x2yyxz=d1*rox3y2z+d2*roxy2z
      x2zzxy=d1*rox3yz2+d2*roxyz2
      y2xxyz=d1*rox2y3z+d2*rox2yz
      y2zzxy=d1*roxy3z2+d2*roxyz2
      z2xxyz=d1*rox2yz3+d2*rox2yz
      z2yyxz=d1*roxy2z3+d2*roxy2z
!c21
      x2yyyz=d1*rox2y3z+d2*roy3z
      x2zzyz=d1*rox2yz3+d2*royz3
      y2xxxz=d1*rox3y2z+d2*rox3z
      y2zzxz=d1*roxy2z3+d2*roxz3
      z2xxxy=d1*rox3yz2+d2*rox3y
      z2yyxy=d1*roxy3z2+d2*roxy3
!c22
      x2yzx2=d1*rox4yz+(d2+d4)*rox2yz+d6*royz
      x2zyx2=d1*rox4yz+(d2+d4)*rox2yz+d6*royz
      y2xzy2=d1*roxy4z+(d2+d4)*roxy2z+d6*roxz
      y2zxy2=d1*roxy4z+(d2+d4)*roxy2z+d6*roxz
      z2xyz2=d1*roxyz4+(d2+d4)*roxyz2+d6*roxy
      z2yxz2=d1*roxyz4+(d2+d4)*roxyz2+d6*roxy
!c23
      x2yzy2=d1*rox2y3z+d2*roy3z+d4*rox2yz+d6*royz
      x2zyz2=d1*rox2yz3+d2*royz3+d4*rox2yz+d6*royz
      y2xzx2=d1*rox3y2z+d2*rox3z+d4*roxy2z+d6*roxz
      y2zxz2=d1*roxy2z3+d2*roxz3+d4*roxy2z+d6*roxz
      z2xyx2=d1*rox3yz2+d2*rox3y+d4*roxyz2+d6*roxy
      z2yxy2=d1*roxy3z2+d2*roxy3+d4*roxyz2+d6*roxy
!c24
      x2yzz2=d1*rox2yz3+d2*royz3+d4*rox2yz+d6*royz
      x2zyy2=d1*rox2y3z+d2*roy3z+d4*rox2yz+d6*royz
      y2xzz2=d1*roxy2z3+d2*roxz3+d4*roxy2z+d6*roxz
      y2zxx2=d1*rox3y2z+d2*rox3z+d4*roxy2z+d6*roxz
      z2xyy2=d1*roxy3z2+d2*roxy3+d4*roxyz2+d6*roxy
      z2yxx2=d1*rox3yz2+d2*rox3y+d4*roxyz2+d6*roxy
!c25
      x2yzxy=d1*rox3y2z+d2*roxy2z
      x2zyxz=d1*rox3yz2+d2*roxyz2
      y2xzxy=d1*rox2y3z+d2*rox2yz
      y2zxyz=d1*roxy3z2+d2*roxyz2
      z2xyxz=d1*rox2yz3+d2*rox2yz
      z2yxyz=d1*roxy2z3+d2*roxy2z
!c26
      x2yzxz=d1*rox3yz2+d2*roxyz2
      x2zyxy=d1*rox3y2z+d2*roxy2z
      y2xzyz=d1*roxy3z2+d2*roxyz2
      y2zxxy=d1*rox2y3z+d2*rox2yz
      z2xyyz=d1*roxy2z3+d2*roxy2z
      z2yxxz=d1*rox2yz3+d2*rox2yz
!c27
      x2yzyz=d1*rox2y2z2+d2*roy2z2
      x2zyyz=d1*rox2y2z2+d2*roy2z2
      y2xzxz=d1*rox2y2z2+d2*rox2z2
      y2zxxz=d1*rox2y2z2+d2*rox2z2
      z2xyxy=d1*rox2y2z2+d2*rox2y2
      z2yxxy=d1*rox2y2z2+d2*rox2y2
!c28
      xyzxx2=d1*rox4yz+d4*rox2yz
      xyzyy2=d1*roxy4z+d4*roxy2z
      xyzzz2=d1*roxyz4+d4*roxyz2
!c29
      xyzxy2=d1*rox2y3z+d4*rox2yz
      xyzxz2=d1*rox2yz3+d4*rox2yz
      xyzyx2=d1*rox3y2z+d4*roxy2z
      xyzyz2=d1*roxy2z3+d4*roxy2z
      xyzzx2=d1*rox3yz2+d4*roxyz2
      xyzzy2=d1*roxy3z2+d4*roxyz2
!c30
      xyzxxy=d1*rox3y2z
      xyzxxz=d1*rox3yz2
      xyzyxy=d1*rox2y3z
      xyzyyz=d1*roxy3z2
      xyzzxz=d1*rox2yz3
      xyzzyz=d1*roxy2z3
!c31
      xyzxyz=d1*rox2y2z2
      xyzyxz=d1*rox2y2z2
      xyzzxy=d1*rox2y2z2
!c32
      x3yyz=d1*rox3y2z+d3*roxy2z
      y3xxz=d1*rox2y3z+d3*rox2yz
      z3xxy=d1*rox2yz3+d3*rox2yz
      x3zyz=d1*rox3yz2+d3*roxyz2
      y3zxz=d1*roxy3z2+d3*roxyz2
      z3yxy=d1*roxy2z3+d3*roxy2z
! c
      df(1,1,1)=xyzxxy
      df(1,1,2)=xyzyxy
      df(1,1,3)=xyzzxy
      df(2,1,1)=xyzxxz
      df(2,1,2)=xyzyxz
      df(2,1,3)=xyzzxz
      df(3,1,1)=xyzxyz
      df(3,1,2)=xyzyyz
      df(3,1,3)=xyzzyz
      df(4,1,1)=(xyzxx2-xyzxy2)
      df(4,1,2)=(xyzyx2-xyzyy2)
      df(4,1,3)=(xyzzx2-xyzzy2)
      df(5,1,1)=(2.d0*xyzxz2-xyzxx2-xyzxy2)
      df(5,1,2)=(2.d0*xyzyz2-xyzyx2-xyzyy2)
      df(5,1,3)=(2.d0*xyzzz2-xyzzx2-xyzzy2)
      df(1,2,1)=(x2zxxy-y2zxxy)
      df(1,2,2)=(x2zyxy-y2zyxy)
      df(1,2,3)=(x2zzxy-y2zzxy)
      df(2,2,1)=(x2zxxz-y2zxxz)
      df(2,2,2)=(x2zyxz-y2zyxz)
      df(2,2,3)=(x2zzxz-y2zzxz)
      df(3,2,1)=(x2zxyz-y2zxyz)
      df(3,2,2)=(x2zyyz-y2zyyz)
      df(3,2,3)=(x2zzyz-y2zzyz)
      df(4,2,1)=(x2zxx2-y2zxx2-x2zxy2+y2zxy2)
      df(4,2,2)=(x2zyx2-y2zyx2-x2zyy2+y2zyy2)
      df(4,2,3)=(x2zzx2-y2zzx2-x2zzy2+y2zzy2)
      df(5,2,1)=(2.d0*x2zxz2-2.d0*y2zxz2-x2zxx2+y2zxx2 &
     & -x2zxy2+y2zxy2)
      df(5,2,2)=(2.d0*x2zyz2-2.d0*y2zyz2-x2zyx2+y2zyx2 &
     & -x2zyy2+y2zyy2)
      df(5,2,3)=(2.d0*x2zzz2-2.d0*y2zzz2-x2zzx2+y2zzx2 &
     & -x2zzy2+y2zzy2)
      df(1,3,1)=(x3xxy-3.d0*y2xxxy)
      df(1,3,2)=(x3yxy-3.d0*y2xyxy)
      df(1,3,3)=(x3zxy-3.d0*y2xzxy)
      df(2,3,1)=(x3xxz-3.d0*y2xxxz)
      df(2,3,2)=(x3yxz-3.d0*y2xyxz)
      df(2,3,3)=(x3zxz-3.d0*y2xzxz)
      df(3,3,1)=(x3xyz-3.d0*y2xxyz)
      df(3,3,2)=(x3yyz-3.d0*y2xyyz)
      df(3,3,3)=(x3zyz-3.d0*y2xzyz)
      df(4,3,1)=(x3xx2-3.d0*y2xxx2-x3xy2+3.d0*y2xxy2)
      df(4,3,2)=(x3yx2-3.d0*y2xyx2-x3yy2+3.d0*y2xyy2)
      df(4,3,3)=(x3zx2-3.d0*y2xzx2-x3zy2+3.d0*y2xzy2)
      df(5,3,1)=(2.d0*x3xz2-6.d0*y2xxz2-x3xx2 &
     & +3.d0*y2xxx2-x3xy2+3.d0*y2xxy2)
      df(5,3,2)=(2.d0*x3yz2-6.d0*y2xyz2-x3yx2 &
     & +3.d0*y2xyx2-x3yy2+3.d0*y2xyy2)
      df(5,3,3)=(2.d0*x3zz2-6.d0*y2xzz2-x3zx2 &
     & +3.d0*y2xzx2-x3zy2+3.d0*y2xzy2)
      df(1,4,1)=(3.d0*x2yxxy-y3xxy)
      df(1,4,2)=(3.d0*x2yyxy-y3yxy)
      df(1,4,3)=(3.d0*x2yzxy-y3zxy)
      df(2,4,1)=(3.d0*x2yxxz-y3xxz)
      df(2,4,2)=(3.d0*x2yyxz-y3yxz)
      df(2,4,3)=(3.d0*x2yzxz-y3zxz)
      df(3,4,1)=(3.d0*x2yxyz-y3xyz)
      df(3,4,2)=(3.d0*x2yyyz-y3yyz)
      df(3,4,3)=(3.d0*x2yzyz-y3zyz)
      df(4,4,1)=(3.d0*x2yxx2-y3xx2-3.d0*x2yxy2+y3xy2)
      df(4,4,2)=(3.d0*x2yyx2-y3yx2-3.d0*x2yyy2+y3yy2)
      df(4,4,3)=(3.d0*x2yzx2-y3zx2-3.d0*x2yzy2+y3zy2)
      df(5,4,1)=(6.d0*x2yxz2-2.d0*y3xz2-3.d0*x2yxx2 &
     & +y3xx2-3.d0*x2yxy2+y3xy2)
      df(5,4,2)=(6.d0*x2yyz2-2.d0*y3yz2-3.d0*x2yyx2 &
     & +y3yx2-3.d0*x2yyy2+y3yy2)
      df(5,4,3)=(6.d0*x2yzz2-2.d0*y3zz2-3.d0*x2yzx2 &
     & +y3zx2-3.d0*x2yzy2+y3zy2)
      df(1,5,1)=2.d0*z3xxy-3.d0*x2zxxy-3.d0*y2zxxy
      df(1,5,2)=2.d0*z3yxy-3.d0*x2zyxy-3.d0*y2zyxy
      df(1,5,3)=2.d0*z3zxy-3.d0*x2zzxy-3.d0*y2zzxy
      df(2,5,1)=2.d0*z3xxz-3.d0*x2zxxz-3.d0*y2zxxz
      df(2,5,2)=2.d0*z3yxz-3.d0*x2zyxz-3.d0*y2zyxz
      df(2,5,3)=2.d0*z3zxz-3.d0*x2zzxz-3.d0*y2zzxz
      df(3,5,1)=2.d0*z3xyz-3.d0*x2zxyz-3.d0*y2zxyz
      df(3,5,2)=2.d0*z3yyz-3.d0*x2zyyz-3.d0*y2zyyz
      df(3,5,3)=2.d0*z3zyz-3.d0*x2zzyz-3.d0*y2zzyz
      df(4,5,1)=z3xx2-1.5d0*x2zxx2-1.5d0*y2zxx2 &
     & -z3xy2+1.5d0*x2zxy2+1.5d0*y2zxy2
      df(4,5,2)=z3yx2-1.5d0*x2zyx2-1.5d0*y2zyx2 &
     & -z3yy2+1.5d0*x2zyy2+1.5d0*y2zyy2
      df(4,5,3)=z3zx2-1.5d0*x2zzx2-1.5d0*y2zzx2 &
     & -z3zy2+1.5d0*x2zzy2+1.5d0*y2zzy2
      df(5,5,1)=(4.d0*z3xz2-6.d0*x2zxz2-6.d0*y2zxz2 &
     & -2.d0*z3xx2+3.d0*x2zxx2+3.d0*y2zxx2 &
     & -2.d0*z3xy2+3.d0*x2zxy2+3.d0*y2zxy2)
      df(5,5,2)=(4.d0*z3yz2-6.d0*x2zyz2-6.d0*y2zyz2 &
     & -2.d0*z3yx2+3.d0*x2zyx2+3.d0*y2zyx2 &
     & -2.d0*z3yy2+3.d0*x2zyy2+3.d0*y2zyy2)
      df(5,5,3)=(4.d0*z3zz2-6.d0*x2zzz2-6.d0*y2zzz2 &
     & -2.d0*z3zx2+3.d0*x2zzx2+3.d0*y2zzx2 &
     & -2.d0*z3zy2+3.d0*x2zzy2+3.d0*y2zzy2)
      df(1,6,1)=(4.d0*z2xxxy-x3xxy-y2xxxy)
      df(1,6,2)=(4.d0*z2xyxy-x3yxy-y2xyxy)
      df(1,6,3)=(4.d0*z2xzxy-x3zxy-y2xzxy)
      df(2,6,1)=(4.d0*z2xxxz-x3xxz-y2xxxz)
      df(2,6,2)=(4.d0*z2xyxz-x3yxz-y2xyxz)
      df(2,6,3)=(4.d0*z2xzxz-x3zxz-y2xzxz)
      df(3,6,1)=(4.d0*z2xxyz-x3xyz-y2xxyz)
      df(3,6,2)=(4.d0*z2xyyz-x3yyz-y2xyyz)
      df(3,6,3)=(4.d0*z2xzyz-x3zyz-y2xzyz)
      df(4,6,1)=(4.d0*z2xxx2-x3xx2-y2xxx2 &
     & -4.d0*z2xxy2+x3xy2+y2xxy2)
      df(4,6,2)=(4.d0*z2xyx2-x3yx2-y2xyx2 &
     & -4.d0*z2xyy2+x3yy2+y2xyy2)
      df(4,6,3)=(4.d0*z2xzx2-x3zx2-y2xzx2 &
     & -4.d0*z2xzy2+x3zy2+y2xzy2)
      df(5,6,1)=(8.d0*z2xxz2-2.d0*x3xz2-2.d0*y2xxz2 &
     & -4.d0*z2xxx2+x3xx2+y2xxx2 &
     & -4.d0*z2xxy2+x3xy2+y2xxy2)
      df(5,6,2)=(8.d0*z2xyz2-2.d0*x3yz2-2.d0*y2xyz2 &
     & -4.d0*z2xyx2+x3yx2+y2xyx2 &
     & -4.d0*z2xyy2+x3yy2+y2xyy2)
      df(5,6,3)=(8.d0*z2xzz2-2.d0*x3zz2-2.d0*y2xzz2 &
     & -4.d0*z2xzx2+x3zx2+y2xzx2 &
     & -4.d0*z2xzy2+x3zy2+y2xzy2)
      df(1,7,1)=(4.d0*z2yxxy-x2yxxy-y3xxy)
      df(1,7,2)=(4.d0*z2yyxy-x2yyxy-y3yxy)
      df(1,7,3)=(4.d0*z2yzxy-x2yzxy-y3zxy)
      df(2,7,1)=(4.d0*z2yxxz-x2yxxz-y3xxz)
      df(2,7,2)=(4.d0*z2yyxz-x2yyxz-y3yxz)
      df(2,7,3)=(4.d0*z2yzxz-x2yzxz-y3zxz)
      df(3,7,1)=(4.d0*z2yxyz-x2yxyz-y3xyz)
      df(3,7,2)=(4.d0*z2yyyz-x2yyyz-y3yyz)
      df(3,7,3)=(4.d0*z2yzyz-x2yzyz-y3zyz)
      df(4,7,1)=(4.d0*z2yxx2-x2yxx2-y3xx2 &
     & -4.d0*z2yxy2+x2yxy2+y3xy2)
      df(4,7,2)=(4.d0*z2yyx2-x2yyx2-y3yx2 &
     & -4.d0*z2yyy2+x2yyy2+y3yy2)
      df(4,7,3)=(4.d0*z2yzx2-x2yzx2-y3zx2 &
     & -4.d0*z2yzy2+x2yzy2+y3zy2)
      df(5,7,1)=(8.d0*z2yxz2-2.d0*x2yxz2-2.d0*y3xz2 &
     & -4.d0*z2yxx2+x2yxx2+y3xx2-4.d0*z2yxy2+x2yxy2+y3xy2)
      df(5,7,2)=(8.d0*z2yyz2-2.d0*x2yyz2-2.d0*y3yz2 &
     & -4.d0*z2yyx2+x2yyx2+y3yx2-4.d0*z2yyy2+x2yyy2+y3yy2)
      df(5,7,3)=(8.d0*z2yzz2-2.d0*x2yzz2-2.d0*y3zz2 &
     & -4.d0*z2yzx2+x2yzx2+y3zx2-4.d0*z2yzy2+x2yzy2+y3zy2)
! c
      do 1075 k=1,7
      do 1075 m=1,5
      do 1075 l=1,3
      q(4+m,9+k,l)=df(m,k,l)
 1075 continue
 1080 continue
! c
      if(nx1.lt.16.or.nx2.lt.16) return
      f1=delt/(64.d0*a13*a23)
      f2=delt/(32.d0*a12*a23)
      f3=3.d0*f2
      f4=delt/(32.d0*a13*a22)
      f5=3.d0*f4
      f6=delt/(16.d0*a12*a22)
      f7=3.d0*f6
      f8=9.d0*f6
!c1
      x3xx3=f1*rox7+(f3+f5)*rox5+f8*rox3
      y3yy3=f1*roy7+(f3+f5)*roy5+f8*roy3
      z3zz3=f1*roz7+(f3+f5)*roz5+f8*roz3
!c2
      x3xy3=f1*rox4y3+f3*rox2y3+f5*rox4y+f8*rox2y
      x3xz3=f1*rox4z3+f3*rox2z3+f5*rox4z+f8*rox2z
      y3yx3=f1*rox3y4+f3*rox3y2+f5*roxy4+f8*roxy2
      y3yz3=f1*roy4z3+f3*roy2z3+f5*roy4z+f8*roy2z
      z3zx3=f1*rox3z4+f3*rox3z2+f5*roxz4+f8*roxz2
      z3zy3=f1*roy3z4+f3*roy3z2+f5*royz4+f8*royz2
!c3
      x3xx2y=f1*rox6y+(f3+f4)*rox4y+f7*rox2y
      x3xx2z=f1*rox6z+(f3+f4)*rox4z+f7*rox2z
      y3yy2x=f1*roxy6+(f3+f4)*roxy4+f7*roxy2
      y3yy2z=f1*roy6z+(f3+f4)*roy4z+f7*roy2z
      z3zz2x=f1*roxz6+(f3+f4)*roxz4+f7*roxz2
      z3zz2y=f1*royz6+(f3+f4)*royz4+f7*royz2
!c4
      x3xy2x=f1*rox5y2+f3*rox3y2+f4*rox5+f7*rox3
      x3xz2x=f1*rox5z2+f3*rox3z2+f4*rox5+f7*rox3
      y3yx2y=f1*rox2y5+f3*rox2y3+f4*roy5+f7*roy3
      y3yz2y=f1*roy5z2+f3*roy3z2+f4*roy5+f7*roy3
      z3zx2z=f1*rox2z5+f3*rox2z3+f4*roz5+f7*roz3
      z3zy2z=f1*roy2z5+f3*roy2z3+f4*roz5+f7*roz3
!c5
      x3xy2z=f1*rox4y2z+f3*rox2y2z+f4*rox4z+f7*rox2z
      x3xz2y=f1*rox4yz2+f3*rox2yz2+f4*rox4y+f7*rox2y
      y3yz2x=f1*roxy4z2+f3*roxy2z2+f4*roxy4+f7*roxy2
      y3yx2z=f1*rox2y4z+f3*rox2y2z+f4*roy4z+f7*roy2z
      z3zy2x=f1*roxy2z4+f3*roxy2z2+f4*roxz4+f7*roxz2
      z3zx2y=f1*rox2yz4+f3*rox2yz2+f4*royz4+f7*royz2
!c6
      x3xxyz=f1*rox5yz+f3*rox3yz
      y3yxyz=f1*roxy5z+f3*roxy3z
      z3zxyz=f1*roxyz5+f3*roxyz3
!c7
      x3yx3=f1*rox6y+(f3+f5)*rox4y+f8*rox2y
      x3zx3=f1*rox6z+(f3+f5)*rox4z+f8*rox2z
      y3xy3=f1*roxy6+(f3+f5)*roxy4+f8*roxy2
      y3zy3=f1*roy6z+(f3+f5)*roy4z+f8*roy2z
      z3xz3=f1*roxz6+(f3+f5)*roxz4+f8*roxz2
      z3yz3=f1*royz6+(f3+f5)*royz4+f8*royz2
!c8
      x3yy3=f1*rox3y4+f3*roxy4+f5*rox3y2+f8*roxy2
      x3zz3=f1*rox3z4+f3*roxz4+f5*rox3z2+f8*roxz2
      y3xx3=f1*rox4y3+f3*rox4y+f5*rox2y3+f8*rox2y
      y3zz3=f1*roy3z4+f3*royz4+f5*roy3z2+f8*royz2
      z3xx3=f1*rox4z3+f3*rox4z+f5*rox2z3+f8*rox2z
      z3yy3=f1*roy4z3+f3*roy4z+f5*roy2z3+f8*roy2z
!c9
      x3yz3=f1*rox3yz3+f3*roxyz3+f5*rox3yz+f8*roxyz
      x3zy3=f1*rox3y3z+f3*roxy3z+f5*rox3yz+f8*roxyz
      y3xz3=f1*roxy3z3+f3*roxyz3+f5*roxy3z+f8*roxyz
      y3zx3=f1*rox3y3z+f3*rox3yz+f5*roxy3z+f8*roxyz
      z3xy3=f1*roxy3z3+f3*roxy3z+f5*roxyz3+f8*roxyz
      z3yx3=f1*rox3yz3+f3*rox3yz+f5*roxyz3+f8*roxyz
!c10
      x3yx2y=f1*rox5y2+(f3+f4)*rox3y2+f7*roxy2
      x3zx2z=f1*rox5z2+(f3+f4)*rox3z2+f7*roxz2
      y3xy2x=f1*rox2y5+(f3+f4)*rox2y3+f7*rox2y
      y3zy2z=f1*roy5z2+(f3+f4)*roy3z2+f7*royz2
      z3xz2x=f1*rox2z5+(f3+f4)*rox2z3+f7*rox2z
      z3yz2y=f1*roy2z5+(f3+f4)*roy2z3+f7*roy2z
!c11
      x3yx2z=f1*rox5yz+(f3+f4)*rox3yz+f7*roxyz
      x3zx2y=f1*rox5yz+(f3+f4)*rox3yz+f7*roxyz
      y3xy2z=f1*roxy5z+(f3+f4)*roxy3z+f7*roxyz
      y3zy2x=f1*roxy5z+(f3+f4)*roxy3z+f7*roxyz
      z3xz2y=f1*roxyz5+(f3+f4)*roxyz3+f7*roxyz
      z3yz2x=f1*roxyz5+(f3+f4)*roxyz3+f7*roxyz
!c12
      x3yy2x=f1*rox4y3+f3*rox2y3+f4*rox4y+f7*rox2y
      x3zz2x=f1*rox4z3+f3*rox2z3+f4*rox4z+f7*rox2z
      y3xx2y=f1*rox3y4+f3*rox3y2+f4*roxy4+f7*roxy2
      y3zz2y=f1*roy4z3+f3*roy2z3+f4*roy4z+f7*roy2z
      z3xx2z=f1*rox3z4+f3*rox3z2+f4*roxz4+f7*roxz2
      z3yy2z=f1*roy3z4+f3*roy3z2+f4*royz4+f7*royz2
!c13
      x3yy2z=f1*rox3y3z+f3*roxy3z+f4*rox3yz+f7*roxyz
      x3zz2y=f1*rox3yz3+f3*roxyz3+f4*rox3yz+f7*roxyz
      y3xx2z=f1*rox3y3z+f3*rox3yz+f4*roxy3z+f7*roxyz
      y3zz2x=f1*roxy3z3+f3*roxyz3+f4*roxy3z+f7*roxyz
      z3xx2y=f1*rox3yz3+f3*rox3yz+f4*roxyz3+f7*roxyz
      z3yy2x=f1*roxy3z3+f3*roxy3z+f4*roxyz3+f7*roxyz
!c14
      x3yz2x=f1*rox4yz2+f3*rox2yz2+f4*rox4y+f7*rox2y
      x3zy2x=f1*rox4y2z+f3*rox2y2z+f4*rox4z+f7*rox2z
      y3xz2y=f1*roxy4z2+f3*roxy2z2+f4*roxy4+f7*roxy2
      y3zx2y=f1*rox2y4z+f3*rox2y2z+f4*roy4z+f7*roy2z
      z3xy2z=f1*roxy2z4+f3*roxy2z2+f4*roxz4+f7*roxz2
      z3yx2z=f1*rox2yz4+f3*rox2yz2+f4*royz4+f7*royz2
!c15
      x3yz2y=f1*rox3y2z2+f3*roxy2z2+f4*rox3y2+f7*roxy2
      x3zy2z=f1*rox3y2z2+f3*roxy2z2+f4*rox3z2+f7*roxz2
      y3xz2x=f1*rox2y3z2+f3*rox2yz2+f4*rox2y3+f7*rox2y
      y3zx2z=f1*rox2y3z2+f3*rox2yz2+f4*roy3z2+f7*royz2
      z3xy2x=f1*rox2y2z3+f3*rox2y2z+f4*rox2z3+f7*rox2z
      z3yx2y=f1*rox2y2z3+f3*rox2y2z+f4*roy2z3+f7*roy2z
!c16
      x3yxyz=f1*rox4y2z+f3*rox2y2z
      x3zxyz=f1*rox4yz2+f3*rox2yz2
      y3xxyz=f1*rox2y4z+f3*rox2y2z
      y3zxyz=f1*roxy4z2+f3*roxy2z2
      z3xxyz=f1*rox2yz4+f3*rox2yz2
      z3yxyz=f1*roxy2z4+f3*roxy2z2
!c17
      x2yxx3=f1*rox6y+(f2+f5)*rox4y+f7*rox2y
      x2zxx3=f1*rox6z+(f2+f5)*rox4z+f7*rox2z
      y2xyy3=f1*roxy6+(f2+f5)*roxy4+f7*roxy2
      y2zyy3=f1*roy6z+(f2+f5)*roy4z+f7*roy2z
      z2xzz3=f1*roxz6+(f2+f5)*roxz4+f7*roxz2
      z2yzz3=f1*royz6+(f2+f5)*royz4+f7*royz2
!c18
      x2yxy3=f1*rox3y4+f2*roxy4+f5*rox3y2+f7*roxy2
      x2zxz3=f1*rox3z4+f2*roxz4+f5*rox3z2+f7*roxz2
      y2xyx3=f1*rox4y3+f2*rox4y+f5*rox2y3+f7*rox2y
      y2zyz3=f1*roy3z4+f2*royz4+f5*roy3z2+f7*royz2
      z2xzx3=f1*rox4z3+f2*rox4z+f5*rox2z3+f7*rox2z
      z2yzy3=f1*roy4z3+f2*roy4z+f5*roy2z3+f7*roy2z
!c19
      x2yxz3=f1*rox3yz3+f2*roxyz3+f5*rox3yz+f7*roxyz
      x2zxy3=f1*rox3y3z+f2*roxy3z+f5*rox3yz+f7*roxyz
      y2xyz3=f1*roxy3z3+f2*roxyz3+f5*roxy3z+f7*roxyz
      y2zyx3=f1*rox3y3z+f2*rox3yz+f5*roxy3z+f7*roxyz
      z2xzy3=f1*roxy3z3+f2*roxy3z+f5*roxyz3+f7*roxyz
      z2yzx3=f1*rox3yz3+f2*rox3yz+f5*roxyz3+f7*roxyz
!c20
      x2yxx2y=f1*rox5y2+(f2+f4)*rox3y2+f6*roxy2
      x2zxx2z=f1*rox5z2+(f2+f4)*rox3z2+f6*roxz2
      y2xyy2x=f1*rox2y5+(f2+f4)*rox2y3+f6*rox2y
      y2zyy2z=f1*roy5z2+(f2+f4)*roy3z2+f6*royz2
      z2xzz2x=f1*rox2z5+(f2+f4)*rox2z3+f6*rox2z
      z2yzz2y=f1*roy2z5+(f2+f4)*roy2z3+f6*roy2z
!c21
      x2yxx2z=f1*rox5yz+(f2+f4)*rox3yz+f6*roxyz
      x2zxx2y=f1*rox5yz+(f2+f4)*rox3yz+f6*roxyz
      y2xyy2z=f1*roxy5z+(f2+f4)*roxy3z+f6*roxyz
      y2zyy2x=f1*roxy5z+(f2+f4)*roxy3z+f6*roxyz
      z2xzz2y=f1*roxyz5+(f2+f4)*roxyz3+f6*roxyz
      z2yzz2x=f1*roxyz5+(f2+f4)*roxyz3+f6*roxyz
!c22
      x2yxy2x=f1*rox4y3+f2*rox2y3+f4*rox4y+f6*rox2y
      x2zxz2x=f1*rox4z3+f2*rox2z3+f4*rox4z+f6*rox2z
      y2xyx2y=f1*rox3y4+f2*rox3y2+f4*roxy4+f6*roxy2
      y2zyz2y=f1*roy4z3+f2*roy2z3+f4*roy4z+f6*roy2z
      z2xzx2z=f1*rox3z4+f2*rox3z2+f4*roxz4+f6*roxz2
      z2yzy2z=f1*roy3z4+f2*roy3z2+f4*royz4+f6*royz2
!c23
      x2yxy2z=f1*rox3y3z+f2*roxy3z+f4*rox3yz+f6*roxyz
      x2zxz2y=f1*rox3yz3+f2*roxyz3+f4*rox3yz+f6*roxyz
      y2xyx2z=f1*rox3y3z+f2*rox3yz+f4*roxy3z+f6*roxyz
      y2zyz2x=f1*roxy3z3+f2*roxyz3+f4*roxy3z+f6*roxyz
      z2xzx2y=f1*rox3yz3+f2*rox3yz+f4*roxyz3+f6*roxyz
      z2yzy2x=f1*roxy3z3+f2*roxy3z+f4*roxyz3+f6*roxyz
!c24
      x2yxz2x=f1*rox4yz2+f2*rox2yz2+f4*rox4y+f6*rox2y
      x2zxy2x=f1*rox4y2z+f2*rox2y2z+f4*rox4z+f6*rox2z
      y2xyz2y=f1*roxy4z2+f2*roxy2z2+f4*roxy4+f6*roxy2
      y2zyx2y=f1*rox2y4z+f2*rox2y2z+f4*roy4z+f6*roy2z
      z2xzy2z=f1*roxy2z4+f2*roxy2z2+f4*roxz4+f6*roxz2
      z2yzx2z=f1*rox2yz4+f2*rox2yz2+f4*royz4+f6*royz2
!c25
      x2yxz2y=f1*rox3y2z2+f2*roxy2z2+f4*rox3y2+f6*roxy2
      x2zxy2z=f1*rox3y2z2+f2*roxy2z2+f4*rox3z2+f6*roxz2
      y2xyz2x=f1*rox2y3z2+f2*rox2yz2+f4*rox2y3+f6*rox2y
      y2zyx2z=f1*rox2y3z2+f2*rox2yz2+f4*roy3z2+f6*royz2
      z2xzy2x=f1*rox2y2z3+f2*rox2y2z+f4*rox2z3+f6*rox2z
      z2yzx2y=f1*rox2y2z3+f2*rox2y2z+f4*roy2z3+f6*roy2z
!c26
      x2yxxyz=f1*rox4y2z+f2*rox2y2z
      x2zxxyz=f1*rox4yz2+f2*rox2yz2
      y2xyxyz=f1*rox2y4z+f2*rox2y2z
      y2zyxyz=f1*roxy4z2+f2*roxy2z2
      z2xzxyz=f1*rox2yz4+f2*rox2yz2
      z2yzxyz=f1*roxy2z4+f2*roxy2z2
!c27
      x2yyx3=f1*rox5y2+(f2+f5)*rox3y2+f7*roxy2
      x2zzx3=f1*rox5z2+(f2+f5)*rox3z2+f7*roxz2
      y2xxy3=f1*rox2y5+(f2+f5)*rox2y3+f7*rox2y
      y2zzy3=f1*roy5z2+(f2+f5)*roy3z2+f7*royz2
      z2xxz3=f1*rox2z5+(f2+f5)*rox2z3+f7*rox2z
      z2yyz3=f1*roy2z5+(f2+f5)*roy2z3+f7*roy2z
!c28
      x2yyy3=f1*rox2y5+f2*roy5+f5*rox2y3+f7*roy3
      x2zzz3=f1*rox2z5+f2*roz5+f5*rox2z3+f7*roz3
      y2xxx3=f1*rox5y2+f2*rox5+f5*rox3y2+f7*rox3
      y2zzz3=f1*roy2z5+f2*roz5+f5*roy2z3+f7*roz3
      z2xxx3=f1*rox5z2+f2*rox5+f5*rox3z2+f7*rox3
      z2yyy3=f1*roy5z2+f2*roy5+f5*roy3z2+f7*roy3
!c29
      x2yyz3=f1*rox2y2z3+f2*roy2z3+f5*rox2y2z+f7*roy2z
      x2zzy3=f1*rox2y3z2+f2*roy3z2+f5*rox2yz2+f7*royz2
      y2xxz3=f1*rox2y2z3+f2*rox2z3+f5*rox2y2z+f7*rox2z
      y2zzx3=f1*rox3y2z2+f2*rox3z2+f5*roxy2z2+f7*roxz2
      z2xxy3=f1*rox2y3z2+f2*rox2y3+f5*rox2yz2+f7*rox2y
      z2yyx3=f1*rox3y2z2+f2*rox3y2+f5*roxy2z2+f7*roxy2
!c30
      x2yyx2y=f1*rox4y3+(f2+f4)*rox2y3+f6*roy3
      x2zzx2z=f1*rox4z3+(f2+f4)*rox2z3+f6*roz3
      y2xxy2x=f1*rox3y4+(f2+f4)*rox3y2+f6*rox3
      y2zzy2z=f1*roy4z3+(f2+f4)*roy2z3+f6*roz3
      z2xxz2x=f1*rox3z4+(f2+f4)*rox3z2+f6*rox3
      z2yyz2y=f1*roy3z4+(f2+f4)*roy3z2+f6*roy3
!c31
      x2yyx2z=f1*rox4y2z+(f2+f4)*rox2y2z+f6*roy2z
      x2zzx2y=f1*rox4yz2+(f2+f4)*rox2yz2+f6*royz2
      y2xxy2z=f1*rox2y4z+(f2+f4)*rox2y2z+f6*rox2z
      y2zzy2x=f1*roxy4z2+(f2+f4)*roxy2z2+f6*roxz2
      z2xxz2y=f1*rox2yz4+(f2+f4)*rox2yz2+f6*rox2y
      z2yyz2x=f1*roxy2z4+(f2+f4)*roxy2z2+f6*roxy2
!c32
      x2yyy2x=f1*rox3y4+f2*roxy4+f4*rox3y2+f6*roxy2
      x2zzz2x=f1*rox3z4+f2*roxz4+f4*rox3z2+f6*roxz2
      y2xxx2y=f1*rox4y3+f2*rox4y+f4*rox2y3+f6*rox2y
      y2zzz2y=f1*roy3z4+f2*royz4+f4*roy3z2+f6*royz2
      z2xxx2z=f1*rox4z3+f2*rox4z+f4*rox2z3+f6*rox2z
      z2yyy2z=f1*roy4z3+f2*roy4z+f4*roy2z3+f6*roy2z
!c33
      x2yyy2z=f1*rox2y4z+f2*roy4z+f4*rox2y2z+f6*roy2z
      x2zzz2y=f1*rox2yz4+f2*royz4+f4*rox2yz2+f6*royz2
      y2xxx2z=f1*rox4y2z+f2*rox4z+f4*rox2y2z+f6*rox2z
      y2zzz2x=f1*roxy2z4+f2*roxz4+f4*roxy2z2+f6*roxz2
      z2xxx2y=f1*rox4yz2+f2*rox4y+f4*rox2yz2+f6*rox2y
      z2yyy2x=f1*roxy4z2+f2*roxy4+f4*roxy2z2+f6*roxy2
!c34
      x2yyz2x=f1*rox3y2z2+f2*roxy2z2+f4*rox3y2+f6*roxy2
      x2zzy2x=f1*rox3y2z2+f2*roxy2z2+f4*rox3z2+f6*roxz2
      y2xxz2y=f1*rox2y3z2+f2*rox2yz2+f4*rox2y3+f6*rox2y
      y2zzx2y=f1*rox2y3z2+f2*rox2yz2+f4*roy3z2+f6*royz2
      z2xxy2z=f1*rox2y2z3+f2*rox2y2z+f4*rox2z3+f6*rox2z
      z2yyx2z=f1*rox2y2z3+f2*rox2y2z+f4*roy2z3+f6*roy2z
!c35
      x2yyz2y=f1*rox2y3z2+f2*roy3z2+f4*rox2y3+f6*roy3
      x2zzy2z=f1*rox2y2z3+f2*roy2z3+f4*rox2z3+f6*roz3
      y2xxz2x=f1*rox3y2z2+f2*rox3z2+f4*rox3y2+f6*rox3
      y2zzx2z=f1*rox2y2z3+f2*rox2z3+f4*roy2z3+f6*roz3
      z2xxy2x=f1*rox3y2z2+f2*rox3y2+f4*rox3z2+f6*rox3
      z2yyx2y=f1*rox2y3z2+f2*rox2y3+f4*roy3z2+f6*roy3
!c36
      x2yyxyz=f1*rox3y3z+f2*roxy3z
      x2zzxyz=f1*rox3yz3+f2*roxyz3
      y2xxxyz=f1*rox3y3z+f2*rox3yz
      y2zzxyz=f1*roxy3z3+f2*roxyz3
      z2xxxyz=f1*rox3yz3+f2*rox3yz
      z2yyxyz=f1*roxy3z3+f2*roxy3z
!c37
      x2yzx3=f1*rox5yz+(f2+f5)*rox3yz+f7*roxyz
      x2zyx3=f1*rox5yz+(f2+f5)*rox3yz+f7*roxyz
      y2xzy3=f1*roxy5z+(f2+f5)*roxy3z+f7*roxyz
      y2zxy3=f1*roxy5z+(f2+f5)*roxy3z+f7*roxyz
      z2xyz3=f1*roxyz5+(f2+f5)*roxyz3+f7*roxyz
      z2yxz3=f1*roxyz5+(f2+f5)*roxyz3+f7*roxyz
!c38
      x2yzy3=f1*rox2y4z+f2*roy4z+f5*rox2y2z+f7*roy2z
      x2zyz3=f1*rox2yz4+f2*royz4+f5*rox2yz2+f7*royz2
      y2xzx3=f1*rox4y2z+f2*rox4z+f5*rox2y2z+f7*rox2z
      y2zxz3=f1*roxy2z4+f2*roxz4+f5*roxy2z2+f7*roxz2
      z2xyx3=f1*rox4yz2+f2*rox4y+f5*rox2yz2+f7*rox2y
      z2yxy3=f1*roxy4z2+f2*roxy4+f5*roxy2z2+f7*roxy2
!c39
      x2yzz3=f1*rox2yz4+f2*royz4+f5*rox2yz2+f7*royz2
      x2zyy3=f1*rox2y4z+f2*roy4z+f5*rox2y2z+f7*roy2z
      y2xzz3=f1*roxy2z4+f2*roxz4+f5*roxy2z2+f7*roxz2
      y2zxx3=f1*rox4y2z+f2*rox4z+f5*rox2y2z+f7*rox2z
      z2xyy3=f1*roxy4z2+f2*roxy4+f5*roxy2z2+f7*roxy2
      z2yxx3=f1*rox4yz2+f2*rox4y+f5*rox2yz2+f7*rox2y
!c40
      x2yzx2y=f1*rox4y2z+(f2+f4)*rox2y2z+f6*roy2z
      x2zyx2z=f1*rox4yz2+(f2+f4)*rox2yz2+f6*royz2
      y2xzy2x=f1*rox2y4z+(f2+f4)*rox2y2z+f6*rox2z
      y2zxy2z=f1*roxy4z2+(f2+f4)*roxy2z2+f6*roxz2
      z2xyz2x=f1*rox2yz4+(f2+f4)*rox2yz2+f6*rox2y
      z2yxz2y=f1*roxy2z4+(f2+f4)*roxy2z2+f6*roxy2
!c41
      x2yzx2z=f1*rox4yz2+(f2+f4)*rox2yz2+f6*royz2
      x2zyx2y=f1*rox4y2z+(f2+f4)*rox2y2z+f6*roy2z
      y2xzy2z=f1*roxy4z2+(f2+f4)*roxy2z2+f6*roxz2
      y2zxy2x=f1*rox2y4z+(f2+f4)*rox2y2z+f6*rox2z
      z2xyz2y=f1*roxy2z4+(f2+f4)*roxy2z2+f6*roxy2
      z2yxz2x=f1*rox2yz4+(f2+f4)*rox2yz2+f6*rox2y
!c42
      x2yzy2x=f1*rox3y3z+f2*roxy3z+f4*rox3yz+f6*roxyz
      x2zyz2x=f1*rox3yz3+f2*roxyz3+f4*rox3yz+f6*roxyz
      y2xzx2y=f1*rox3y3z+f2*rox3yz+f4*roxy3z+f6*roxyz
      y2zxz2y=f1*roxy3z3+f2*roxyz3+f4*roxy3z+f6*roxyz
      z2xyx2z=f1*rox3yz3+f2*rox3yz+f4*roxyz3+f6*roxyz
      z2yxy2z=f1*roxy3z3+f2*roxy3z+f4*roxyz3+f6*roxyz
!c43
      x2yzy2z=f1*rox2y3z2+f2*roy3z2+f4*rox2yz2+f6*royz2
      x2zyz2y=f1*rox2y2z3+f2*roy2z3+f4*rox2y2z+f6*roy2z
      y2xzx2z=f1*rox3y2z2+f2*rox3z2+f4*roxy2z2+f6*roxz2
      y2zxz2x=f1*rox2y2z3+f2*rox2z3+f4*rox2y2z+f6*rox2z
      z2xyx2y=f1*rox3y2z2+f2*rox3y2+f4*roxy2z2+f6*roxy2
      z2yxy2x=f1*rox2y3z2+f2*rox2y3+f4*rox2yz2+f6*rox2y
!c44
      x2yzz2x=f1*rox3yz3+f2*roxyz3+f4*rox3yz+f6*roxyz
      x2zyy2x=f1*rox3y3z+f2*roxy3z+f4*rox3yz+f6*roxyz
      y2xzz2y=f1*roxy3z3+f2*roxyz3+f4*roxy3z+f6*roxyz
      y2zxx2y=f1*rox3y3z+f2*rox3yz+f4*roxy3z+f6*roxyz
      z2xyy2z=f1*roxy3z3+f2*roxy3z+f4*roxyz3+f6*roxyz
      z2yxx2z=f1*rox3yz3+f2*rox3yz+f4*roxyz3+f6*roxyz
!c45
      x2yzz2y=f1*rox2y2z3+f2*roy2z3+f4*rox2y2z+f6*roy2z
      x2zyy2z=f1*rox2y3z2+f2*roy3z2+f4*rox2yz2+f6*royz2
      y2xzz2x=f1*rox2y2z3+f2*rox2z3+f4*rox2y2z+f6*rox2z
      y2zxx2z=f1*rox3y2z2+f2*rox3z2+f4*roxy2z2+f6*roxz2
      z2xyy2x=f1*rox2y3z2+f2*rox2y3+f4*rox2yz2+f6*rox2y
      z2yxx2y=f1*rox3y2z2+f2*rox3y2+f4*roxy2z2+f6*roxy2
!c46
      x2yzxyz=f1*rox3y2z2+f2*roxy2z2
      x2zyxyz=f1*rox3y2z2+f2*roxy2z2
      y2xzxyz=f1*rox2y3z2+f2*rox2yz2
      y2zxxyz=f1*rox2y3z2+f2*rox2yz2
      z2xyxyz=f1*rox2y2z3+f2*rox2y2z
      z2yxxyz=f1*rox2y2z3+f2*rox2y2z
!c47
      xyzxx3=f1*rox5yz+f5*rox3yz
      xyzyy3=f1*roxy5z+f5*roxy3z
      xyzzz3=f1*roxyz5+f5*roxyz3
!c48
      xyzxy3=f1*rox2y4z+f5*rox2y2z
      xyzyx3=f1*rox4y2z+f5*rox2y2z
      xyzzx3=f1*rox4yz2+f5*rox2yz2
      xyzxz3=f1*rox2yz4+f5*rox2yz2
      xyzyz3=f1*roxy2z4+f5*roxy2z2
      xyzzy3=f1*roxy4z2+f5*roxy2z2
!c49
      xyzxx2y=f1*rox4y2z+f4*rox2y2z
      xyzxx2z=f1*rox4yz2+f4*rox2yz2
      xyzyy2x=f1*rox2y4z+f4*rox2y2z
      xyzyy2z=f1*roxy4z2+f4*roxy2z2
      xyzzz2x=f1*rox2yz4+f4*rox2yz2
      xyzzz2y=f1*roxy2z4+f4*roxy2z2
!c50
      xyzxy2x=f1*rox3y3z+f4*rox3yz
      xyzxz2x=f1*rox3yz3+f4*rox3yz
      xyzyz2y=f1*roxy3z3+f4*roxy3z
      xyzyx2y=f1*rox3y3z+f4*roxy3z
      xyzzy2z=f1*roxy3z3+f4*roxyz3
      xyzzx2z=f1*rox3yz3+f4*roxyz3
!c51
      xyzxy2z=f1*rox2y3z2+f4*rox2yz2
      xyzxz2y=f1*rox2y2z3+f4*rox2y2z
      xyzyz2x=f1*rox2y2z3+f4*rox2y2z
      xyzyx2z=f1*rox3y2z2+f4*roxy2z2
      xyzzx2y=f1*rox3y2z2+f4*roxy2z2
      xyzzy2x=f1*rox2y3z2+f4*rox2yz2
!c52
      xyzxxyz=f1*rox3y2z2
      xyzyxyz=f1*rox2y3z2
      xyzzxyz=f1*rox2y2z3
! c
      ff(1,1,1)=xyzxxyz
      ff(1,1,2)=xyzyxyz
      ff(1,1,3)=xyzzxyz
      ff(2,1,1)=(x2zxxyz-y2zxxyz)
      ff(2,1,2)=(x2zyxyz-y2zyxyz)
      ff(2,1,3)=(x2zzxyz-y2zzxyz)
      ff(1,2,1)=(xyzxx2z-xyzxy2z)
      ff(1,2,2)=(xyzyx2z-xyzyy2z)
      ff(1,2,3)=(xyzzx2z-xyzzy2z)
      ff(2,2,1)=(x2zxx2z-y2zxx2z-x2zxy2z+y2zxy2z)
      ff(2,2,2)=(x2zyx2z-y2zyx2z-x2zyy2z+y2zyy2z)
      ff(2,2,3)=(x2zzx2z-y2zzx2z-x2zzy2z+y2zzy2z)
      ff(3,1,1)=(x3xxyz-3.d0*y2xxxyz)
      ff(3,1,2)=(x3yxyz-3.d0*y2xyxyz)
      ff(3,1,3)=(x3zxyz-3.d0*y2xzxyz)
      ff(1,3,1)=(xyzxx3-3.d0*xyzxy2x)
      ff(1,3,2)=(xyzyx3-3.d0*xyzyy2x)
      ff(1,3,3)=(xyzzx3-3.d0*xyzzy2x)
      ff(3,2,1)=(x3xx2z-3.d0*y2xxx2z-x3xy2z+3.d0*y2xxy2z)
      ff(3,2,2)=(x3yx2z-3.d0*y2xyx2z-x3yy2z+3.d0*y2xyy2z)
      ff(3,2,3)=(x3zx2z-3.d0*y2xzx2z-x3zy2z+3.d0*y2xzy2z)
      ff(2,3,1)=(x2zxx3-3.d0*x2zxy2x-y2zxx3+3.d0*y2zxy2x)
      ff(2,3,2)=(x2zyx3-3.d0*x2zyy2x-y2zyx3+3.d0*y2zyy2x)
      ff(2,3,3)=(x2zzx3-3.d0*x2zzy2x-y2zzx3+3.d0*y2zzy2x)
      ff(3,3,1)=(x3xx3-3.d0*y2xxx3-3.d0*x3xy2x+9.d0*y2xxy2x)
      ff(3,3,2)=(x3yx3-3.d0*y2xyx3-3.d0*x3yy2x+9.d0*y2xyy2x)
      ff(3,3,3)=(x3zx3-3.d0*y2xzx3-3.d0*x3zy2x+9.d0*y2xzy2x)
      ff(4,1,1)=(3.d0*x2yxxyz-y3xxyz)
      ff(4,1,2)=(3.d0*x2yyxyz-y3yxyz)
      ff(4,1,3)=(3.d0*x2yzxyz-y3zxyz)
      ff(1,4,1)=(3.d0*xyzxx2y-xyzxy3)
      ff(1,4,2)=(3.d0*xyzyx2y-xyzyy3)
      ff(1,4,3)=(3.d0*xyzzx2y-xyzzy3)
      ff(4,2,1)=(3.d0*x2yxx2z-y3xx2z-3.d0*x2yxy2z+y3xy2z)
      ff(4,2,2)=(3.d0*x2yyx2z-y3yx2z-3.d0*x2yyy2z+y3yy2z)
      ff(4,2,3)=(3.d0*x2yzx2z-y3zx2z-3.d0*x2yzy2z+y3zy2z)
      ff(2,4,1)=(3.d0*x2zxx2y-x2zxy3-3.d0*y2zxx2y+y2zxy3)
      ff(2,4,2)=(3.d0*x2zyx2y-x2zyy3-3.d0*y2zyx2y+y2zyy3)
      ff(2,4,3)=(3.d0*x2zzx2y-x2zzy3-3.d0*y2zzx2y+y2zzy3)
      ff(4,3,1)=(3.d0*x2yxx3-y3xx3-9.d0*x2yxy2x+3.d0*y3xy2x)
      ff(4,3,2)=(3.d0*x2yyx3-y3yx3-9.d0*x2yyy2x+3.d0*y3yy2x)
      ff(4,3,3)=(3.d0*x2yzx3-y3zx3-9.d0*x2yzy2x+3.d0*y3zy2x)
      ff(3,4,1)=(3.d0*x3xx2y-x3xy3-9.d0*y2xxx2y+3.d0*y2xxy3)
      ff(3,4,2)=(3.d0*x3yx2y-x3yy3-9.d0*y2xyx2y+3.d0*y2xyy3)
      ff(3,4,3)=(3.d0*x3zx2y-x3zy3-9.d0*y2xzx2y+3.d0*y2xzy3)
      ff(4,4,1)=(9.d0*x2yxx2y-3.d0*y3xx2y-3.d0*x2yxy3+y3xy3)
      ff(4,4,2)=(9.d0*x2yyx2y-3.d0*y3yx2y-3.d0*x2yyy3+y3yy3)
      ff(4,4,3)=(9.d0*x2yzx2y-3.d0*y3zx2y-3.d0*x2yzy3+y3zy3)
      ff(5,1,1)=(2.d0*z3xxyz-3.d0*x2zxxyz-3.d0*y2zxxyz)
      ff(5,1,2)=(2.d0*z3yxyz-3.d0*x2zyxyz-3.d0*y2zyxyz)
      ff(5,1,3)=(2.d0*z3zxyz-3.d0*x2zzxyz-3.d0*y2zzxyz)
      ff(1,5,1)=(2.d0*xyzxz3-3.d0*xyzxx2z-3.d0*xyzxy2z)
      ff(1,5,2)=(2.d0*xyzyz3-3.d0*xyzyx2z-3.d0*xyzyy2z)
      ff(1,5,3)=(2.d0*xyzzz3-3.d0*xyzzx2z-3.d0*xyzzy2z)
      ff(5,2,1)=(2.d0*z3xx2z-3.d0*x2zxx2z-3.d0*y2zxx2z &
     & -2.d0*z3xy2z+3.d0*x2zxy2z+3.d0*y2zxy2z)
      ff(5,2,2)=(2.d0*z3yx2z-3.d0*x2zyx2z-3.d0*y2zyx2z &
     & -2.d0*z3yy2z+3.d0*x2zyy2z+3.d0*y2zyy2z)
      ff(5,2,3)=(2.d0*z3zx2z-3.d0*x2zzx2z-3.d0*y2zzx2z &
     & -2.d0*z3zy2z+3.d0*x2zzy2z+3.d0*y2zzy2z)
      ff(2,5,1)=(2.d0*x2zxz3-3.d0*x2zxx2z-3.d0*x2zxy2z &
     & -2.d0*y2zxz3+3.d0*y2zxx2z+3.d0*y2zxy2z)
      ff(2,5,2)=(2.d0*x2zyz3-3.d0*x2zyx2z-3.d0*x2zyy2z &
     & -2.d0*y2zyz3+3.d0*y2zyx2z+3.d0*y2zyy2z)
      ff(2,5,3)=(2.d0*x2zzz3-3.d0*x2zzx2z-3.d0*x2zzy2z &
     & -2.d0*y2zzz3+3.d0*y2zzx2z+3.d0*y2zzy2z)
      ff(5,3,1)=(2.d0*z3xx3-3.d0*x2zxx3-3.d0*y2zxx3 &
     & -6.d0*z3xy2x+9.d0*x2zxy2x+9.d0*y2zxy2x)
      ff(5,3,2)=(2.d0*z3yx3-3.d0*x2zyx3-3.d0*y2zyx3 &
     & -6.d0*z3yy2x+9.d0*x2zyy2x+9.d0*y2zyy2x)
      ff(5,3,3)=(2.d0*z3zx3-3.d0*x2zzx3-3.d0*y2zzx3 &
     & -6.d0*z3zy2x+9.d0*x2zzy2x+9.d0*y2zzy2x)
      ff(3,5,1)=(2.d0*x3xz3-3.d0*x3xx2z-3.d0*x3xy2z &
     & -6.d0*y2xxz3+9.d0*y2xxx2z+9.d0*y2xxy2z)
      ff(3,5,2)=(2.d0*x3yz3-3.d0*x3yx2z-3.d0*x3yy2z &
     & -6.d0*y2xyz3+9.d0*y2xyx2z+9.d0*y2xyy2z)
      ff(3,5,3)=(2.d0*x3zz3-3.d0*x3zx2z-3.d0*x3zy2z &
     & -6.d0*y2xzz3+9.d0*y2xzx2z+9.d0*y2xzy2z)
      ff(5,4,1)=(6.d0*z3xx2y-9.d0*x2zxx2y-9.d0*y2zxx2y &
     & -2.d0*z3xy3+3.d0*x2zxy3+3.d0*y2zxy3)
      ff(5,4,2)=(6.d0*z3yx2y-9.d0*x2zyx2y-9.d0*y2zyx2y &
     & -2.d0*z3yy3+3.d0*x2zyy3+3.d0*y2zyy3)
      ff(5,4,3)=(6.d0*z3zx2y-9.d0*x2zzx2y-9.d0*y2zzx2y &
     & -2.d0*z3zy3+3.d0*x2zzy3+3.d0*y2zzy3)
      ff(4,5,1)=(6.d0*x2yxz3-9.d0*x2yxx2z-9.d0*x2yxy2z &
     & -2.d0*y3xz3+3.d0*y3xx2z+3.d0*y3xy2z)
      ff(4,5,2)=(6.d0*x2yyz3-9.d0*x2yyx2z-9.d0*x2yyy2z &
     & -2.d0*z3yx3+3.d0*y3yx2z+3.d0*y3yy2z)
      ff(4,5,3)=(6.d0*x2yzz3-9.d0*x2yzx2z-9.d0*x2yzy2z &
     & -2.d0*y3zz3+3.d0*y3zx2z+3.d0*y3zy2z)
      ff(5,5,1)=4.d0*z3xz3-6.d0*x2zxz3-6.d0*y2zxz3-6.d0*z3xx2z+ &
     & 9.d0*x2zxx2z+9.d0*y2zxx2z-6.d0*z3xy2z+9.d0*x2zxy2z+9.d0*y2zxy2z
      ff(5,5,2)=4.d0*z3yz3-6.d0*x2zyz3-6.d0*y2zyz3-6.d0*z3yx2z+ &
     & 9.d0*x2zyx2z+9.d0*y2zyx2z-6.d0*z3yy2z+9.d0*x2zyy2z+9.d0*y2zyy2z
      ff(5,5,3)=4.d0*z3zz3-6.d0*x2zzz3-6.d0*y2zzz3-6.d0*z3zx2z+ &
     & 9.d0*x2zzx2z+9.d0*y2zzx2z-6.d0*z3zy2z+9.d0*x2zzy2z+9.d0*y2zzy2z
      ff(6,1,1)=(4.d0*z2xxxyz-x3xxyz-y2xxxyz)
      ff(6,1,2)=(4.d0*z2xyxyz-x3yxyz-y2xyxyz)
      ff(6,1,3)=(4.d0*z2xzxyz-x3zxyz-y2xzxyz)
      ff(1,6,1)=(4.d0*xyzxz2x-xyzxx3-xyzxy2x)
      ff(1,6,2)=(4.d0*xyzyz2x-xyzyx3-xyzyy2x)
      ff(1,6,3)=(4.d0*xyzzz2x-xyzzx3-xyzzy2x)
      ff(6,2,1)=(4.d0*z2xxx2z-x3xx2z-y2xxx2z &
     & -4.d0*z2xxy2z+x3xy2z+y2xxy2z)
      ff(6,2,2)=(4.d0*z2xyx2z-x3yx2z-y2xyx2z &
     & -4.d0*z2xyy2z+x3yy2z+y2xyy2z)
      ff(6,2,3)=(4.d0*z2xzx2z-x3zx2z-y2xzx2z &
     & -4.d0*z2xzy2z+x3zy2z+y2xzy2z)
      ff(2,6,1)=(4.d0*x2zxz2x-x2zxx3-x2zxy2x &
     & -4.d0*y2zxz2x+y2zxx3+y2zxy2x)
      ff(2,6,2)=(4.d0*x2zyz2x-x2zyx3-x2zyy2x &
     & -4.d0*y2zyz2x+y2zyx3+y2zyy2x)
      ff(2,6,3)=(4.d0*x2zzz2x-x2zzx3-x2zzy2x &
     & -4.d0*y2zzz2x+y2zzx3+y2zzy2x)
      ff(6,3,1)=(4.d0*z2xxx3-x3xx3-y2xxx3 &
     & -12.d0*z2xxy2x+3.d0*x3xy2x+3.d0*y2xxy2x)
      ff(6,3,2)=(4.d0*z2xyx3-x3yx3-y2xyx3 &
     & -12.d0*z2xyy2x+3.d0*x3yy2x+3.d0*y2xyy2x)
      ff(6,3,3)=(4.d0*z2xzx3-x3zx3-y2xzx3 &
     & -12.d0*z2xzy2x+3.d0*x3zy2x+3.d0*y2xzy2x)
      ff(3,6,1)=(4.d0*x3xz2x-x3xx3-x3xy2x &
     & -12.d0*y2xxz2x+3.d0*y2xxx3+3.d0*y2xxy2x)
      ff(3,6,2)=(4.d0*x3yz2x-x3yx3-x3yy2x &
     & -12.d0*y2xyz2x+3.d0*y2xyx3+3.d0*y2xyy2x)
      ff(3,6,3)=(4.d0*x3zz2x-x3zx3-x3zy2x &
     & -12.d0*y2xzz2x+3.d0*y2xzx3+3.d0*y2xzy2x)
      ff(6,4,1)=(12.d0*z2xxx2y-3.d0*x3xx2y-3.d0*y2xxx2y &
     & -4.d0*z2xxy3+x3xy3+y2xxy3)
      ff(6,4,2)=(12.d0*z2xyx2y-3.d0*x3yx2y-3.d0*y2xyx2y &
     & -4.d0*z2xyy3+x3yy3+y2xyy3)
      ff(6,4,3)=(12.d0*z2xzx2y-3.d0*x3zx2y-3.d0*y2xzx2y &
     & -4.d0*z2xzy3+x3zy3+y2xzy3)
      ff(4,6,1)=(12.d0*x2yxz2x-3.d0*x2yxx3-3.d0*x2yxy2x &
     & -4.d0*y3xz2x+y3xx3+y3xy2x)
      ff(4,6,2)=(12.d0*x2yyz2x-3.d0*x2yyx3-3.d0*x2yyy2x &
     & -4.d0*y3yz2x+y3yx3+y3yy2x)
      ff(4,6,3)=(12.d0*x2yzz2x-3.d0*x2yzx3-3.d0*x2yzy2x &
     & -4.d0*y3zz2x+y3zx3+y3zy2x)
      ff(6,5,1)=(8.d0*z2xxz3-2.d0*x3xz3-2.d0*y2xxz3 &
     & -12.d0*z2xxx2z+3.d0*x3xx2z+3.d0*y2xxx2z &
     & -12.d0*z2xxy2z+3.d0*x3xy2z+3.d0*y2xxy2z)
      ff(6,5,2)=(8.d0*z2xyz3-2.d0*x3yz3-2.d0*y2xyz3 &
     & -12.d0*z2xyx2z+3.d0*x3yx2z+3.d0*y2xyx2z &
     & -12.d0*z2xyy2z+3.d0*x3yy2z+3.d0*y2xyy2z)
      ff(6,5,3)=(8.d0*z2xzz3-2.d0*x3zz3-2.d0*y2xzz3 &
     & -12.d0*z2xzx2z+3.d0*x3zx2z+3.d0*y2xzx2z &
     & -12.d0*z2xzy2z+3.d0*x3zy2z+3.d0*y2xzy2z)
      ff(5,6,1)=(8.d0*z3xz2x-2.d0*z3xx3-2.d0*z3xy2x &
     & -12.d0*x2zxz2x+3.d0*x2zxx3+3.d0*x2zxy2x &
     & -12.d0*y2zxz2x+3.d0*y2zxx3+3.d0*y2zxy2x)
      ff(5,6,2)=(8.d0*z3yz2x-2.d0*z3yx3-2.d0*z3yy2x &
     & -12.d0*x2zyz2x+3.d0*x2zyx3+3.d0*x2zyy2x &
     & -12.d0*y2zyz2x+3.d0*y2zyx3+3.d0*y2zyy2x)
      ff(5,6,3)=(8.d0*z3zz2x-2.d0*z3zx3-2.d0*z3zy2x &
     & -12.d0*x2zzz2x+3.d0*x2zzx3+3.d0*x2zzy2x &
     & -12.d0*y2zzz2x+3.d0*y2zzx3+3.d0*y2zzy2x)
      ff(6,6,1)=(16.d0*z2xxz2x-4.d0*x3xz2x-4.d0*y2xxz2x &
     & -4.d0*z2xxx3+x3xx3+y2xxx3-4.d0*z2xxy2x+x3xy2x+y2xxy2x)
      ff(6,6,2)=(16.d0*z2xyz2x-4.d0*x3yz2x-4.d0*y2xyz2x &
     & -4.d0*z2xyx3+x3yx3+y2xyx3-4.d0*z2xyy2x+x3yy2x+y2xyy2x)
      ff(6,6,3)=(16.d0*z2xzz2x-4.d0*x3zz2x-4.d0*y2xzz2x &
     & -4.d0*z2xzx3+x3zx3+y2xzx3-4.d0*z2xzy2x+x3zy2x+y2xzy2x)
      ff(7,1,1)=(4.d0*z2yxxyz-x2yxxyz-y3xxyz)
      ff(7,1,2)=(4.d0*z2yyxyz-x2yyxyz-y3yxyz)
      ff(7,1,3)=(4.d0*z2yzxyz-x2yzxyz-y3zxyz)
      ff(1,7,1)=(4.d0*xyzxz2y-xyzxx2y-xyzxy3)
      ff(1,7,2)=(4.d0*xyzyz2y-xyzyx2y-xyzyy3)
      ff(1,7,3)=(4.d0*xyzzz2y-xyzzx2y-xyzzy3)
      ff(7,2,1)=(4.d0*z2yxx2z-x2yxx2z-y3xx2z &
     & -4.d0*z2yxy2z+x2yxy2z+y3xy2z)
      ff(7,2,2)=(4.d0*z2yyx2z-x2yyx2z-y3yx2z &
     & -4.d0*z2yyy2z+x2yyy2z+y3yy2z)
      ff(7,2,3)=(4.d0*z2yzx2z-x2yzx2z-y3zx2z &
     & -4.d0*z2yzy2z+x2yzy2z+y3zy2z)
      ff(2,7,1)=(4.d0*x2zxz2y-x2zxx2y-x2zxy3 &
     & -4.d0*y2zxz2y+y2zxx2y+y2zxy3)
      ff(2,7,2)=(4.d0*x2zyz2y-x2zyx2y-x2zyy3 &
     & -4.d0*y2zyz2y+y2zyx2y+y2zyy3)
      ff(2,7,3)=(4.d0*x2zzz2y-x2zzx2y-x2zzy3 &
     & -4.d0*y2zzz2y+y2zzx2y+y2zzy3)
      ff(7,3,1)=(4.d0*z2yxx3-x2yxx3-y3xx3 &
     & -12.d0*z2yxy2x+3.d0*x2yxy2x+3.d0*y3xy2x)
      ff(7,3,2)=(4.d0*z2yyx3-x2yyx3-y3yx3 &
     & -12.d0*z2yyy2x+3.d0*x2yyy2x+3.d0*y3yy2x)
      ff(7,3,3)=(4.d0*z2yzx3-x2yzx3-y3zx3 &
     & -12.d0*z2yzy2x+3.d0*x2yzy2x+3.d0*y3zy2x)
      ff(3,7,1)=(4.d0*x3xz2y-x3xx2y-x3xy3 &
     & -12.d0*y2xxz2y+3.d0*y2xxx2y+3.d0*y2xxy3)
      ff(3,7,2)=(4.d0*x3yz2y-x3yx2y-x3yy3 &
     & -12.d0*y2xyz2y+3.d0*y2xyx2y+3.d0*y2xyy3)
      ff(3,7,3)=(4.d0*x3zz2y-x3zx2y-x3zy3 &
     & -12.d0*y2xzz2y+3.d0*y2xzx2y+3.d0*y2xzy3)
      ff(7,4,1)=(12.d0*z2yxx2y-3.d0*x2yxx2y-3.d0*y3xx2y &
     & -4.d0*z2yxy3+x2yxy3+y3xy3)
      ff(7,4,2)=(12.d0*z2yyx2y-3.d0*x2yyx2y-3.d0*y3yx2y &
     & -4.d0*z2yyy3+x2yyy3+y3yy3)
      ff(7,4,3)=(12.d0*z2yzx2y-3.d0*x2yzx2y-3.d0*y3zx2y &
     & -4.d0*z2yzy3+x2yzy3+y3zy3)
      ff(4,7,1)=(12.d0*x2yxz2y-3.d0*x2yxx2y-3.d0*x2yxy3 &
     & -4.d0*y3xz2y+y3xx2y+y3xy3)
      ff(4,7,2)=(12.d0*x2yyz2y-3.d0*x2yyx2y-3.d0*x2yyy3 &
     & -4.d0*y3yz2y+y3yx2y+y3yy3)
      ff(4,7,3)=(12.d0*x2yzz2y-3.d0*x2yzx2y-3.d0*x2yzy3 &
     & -4.d0*y3zz2y+y3zx2y+y3zy3)
      ff(7,5,1)=(8.d0*z2yxz3-2.d0*x2yxz3-2.d0*y3xz3 &
     & -12.d0*z2yxx2z+3.d0*x2yxx2z+3.d0*y3xx2z &
     & -12.d0*z2yxy2z+3.d0*x2yxy2z+3.d0*y3xy2z)
      ff(7,5,2)=(8.d0*z2yyz3-2.d0*x2yyz3-2.d0*y3yz3 &
     & -12.d0*z2yyx2z+3.d0*x2yyx2z+3.d0*y3yx2z &
     & -12.d0*z2yyy2z+3.d0*x2yyy2z+3.d0*y3yy2z)
      ff(7,5,3)=(8.d0*z2yzz3-2.d0*x2yzz3-2.d0*y3zz3 &
     & -12.d0*z2yzx2z+3.d0*x2yzx2z+3.d0*y3zx2z &
     & -12.d0*z2yzy2z+3.d0*x2yzy2z+3.d0*y3zy2z)
      ff(5,7,1)=(8.d0*z3xz2y-2.d0*z3xx2y-2.d0*z3xy3 &
     & -12.d0*x2zxz2y+3.d0*x2zxx2y+3.d0*x2zxy3 &
     & -12.d0*y2zxz2y+3.d0*y2zxx2y+3.d0*y2zxy3)
      ff(5,7,2)=(8.d0*z3yz2y-2.d0*z3yx2y-2.d0*z3yy3 &
     & -12.d0*x2zyz2y+3.d0*x2zyx2y+3.d0*x2zyy3 &
     & -12.d0*y2zyz2y+3.d0*y2zyx2y+3.d0*y2zyy3)
      ff(5,7,3)=(8.d0*z3zz2y-2.d0*z3zx2y-2.d0*z3zy3 &
     & -12.d0*x2zzz2y+3.d0*x2zzx2y+3.d0*x2zzy3 &
     & -12.d0*y2zzz2y+3.d0*y2zzx2y+3.d0*y2zzy3)
      ff(7,6,1)=(16.d0*z2yxz2x-4.d0*x2yxz2x-4.d0*y3xz2x &
     & -4.d0*z2yxx3+x2yxx3+y3xx3-4.d0*z2yxy2x+x2yxy2x+y3xy2x)
      ff(7,6,2)=(16.d0*z2yyz2x-4.d0*x2yyz2x-4.d0*y3yz2x &
     & -4.d0*z2yyx3+x2yyx3+y3yx3-4.d0*z2yyy2x+x2yyy2x+y3yy2x)
      ff(7,6,3)=(16.d0*z2yzz2x-4.d0*x2yzz2x-4.d0*y3zz2x &
     & -4.d0*z2yzx3+x2yzx3+y3zx3-4.d0*z2yzy2x+x2yzy2x+y3zy2x)
      ff(6,7,1)=(16.d0*z2xxz2y-4.d0*z2xxx2y-4.d0*z2xxy3 &
     & -4.d0*x3xz2y+x3xx2y+x3xy3-4.d0*y2xxz2y+y2xxx2y+y2xxy3)
      ff(6,7,2)=(16.d0*z2xyz2y-4.d0*z2xyx2y-4.d0*z2xyy3 &
     & -4.d0*x3yz2y+x3yx2y+x3yy3-4.d0*y2xyz2y+y2xyx2y+y2xyy3)
      ff(6,7,3)=(16.d0*z2xzz2y-4.d0*z2xzx2y-4.d0*z2xzy3 &
     & -4.d0*x3zz2y+x3zx2y+x3zy3-4.d0*y2xzz2y+y2xzx2y+y2xzy3)
      ff(7,7,1)=(16.d0*z2yxz2y-4.d0*x2yxz2y-4.d0*y3xz2y &
     & -4.d0*z2yxx2y+x2yxx2y+y3xx2y-4.d0*z2yxy3+x2yxy3+y3xy3)
      ff(7,7,2)=(16.d0*z2yyz2y-4.d0*x2yyz2y-4.d0*y3yz2y &
     & -4.d0*z2yyx2y+x2yyx2y+y3yx2y-4.d0*z2yyy3+x2yyy3+y3yy3)
      ff(7,7,3)=(16.d0*z2yzz2y-4.d0*x2yzz2y-4.d0*y3zz2y &
     & -4.d0*z2yzx2y+x2yzx2y+y3zx2y-4.d0*z2yzy3+x2yzy3+y3zy3)
! c
      do 1085 k=1,7
      do 1085 m=1,7
      do 1085 l=1,3
      q(9+k,9+m,l)=ff(k,m,l)
 1085 continue
      return
      end subroutine momf

end module O_GaussianIntegrals

