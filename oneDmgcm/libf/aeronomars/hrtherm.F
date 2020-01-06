c**********************************************************************

      subroutine hrtherm(ig,euvmod,rm,nespeuv,tx,iz,zenit,zday,jtot)


c     feb 2002        fgg           first version
c     nov 2002        fgg           second version

c**********************************************************************

      implicit none

c     common variables and constants


      include 'dimensions.h'
      include 'dimphys.h'
      include 'param.h'
      include 'param_v4.h'
      include "callkeys.h"


c    local parameters and variables

      real       xabsi(nabs,nlayermx) 			!densities
      real       jergs(ninter,nabs,nlayermx)
      
      integer    i,j,k,indexint          !indexes
      character  dn


c     input and output variables

      integer    ig  ,euvmod 
      integer    nespeuv
      real       rm(nlayermx,nespeuv)              !density matrix (cm^-3)
      real       jtot(nlayermx)                    !output: heating rate(erg/s)
      real       tx(nlayermx)                      !temperature
      real       zenit
      real       iz(nlayermx)
      real       zday

      ! tracer indexes for the EUV heating:
!!! ATTENTION. These values have to be identical to those in chemthermos.F90
!!! If the values are changed there, the same has to be done here  !!!
      integer,parameter :: i_co2=1
      integer,parameter :: i_o2=2
      integer,parameter :: i_o=3
      integer,parameter :: i_co=4
      integer,parameter :: i_h=5
      integer,parameter :: i_h2=8
      integer,parameter :: i_h2o=9
      integer,parameter :: i_h2o2=10
      integer,parameter :: i_o3=12
      integer,parameter :: i_n2=13
      integer,parameter :: i_n=14
      integer,parameter :: i_no=15
      integer,parameter :: i_no2=17

c*************************PROGRAM STARTS*******************************

      !If nighttime, photoabsorption coefficient set to 0
      if(zenit.gt.140.) then
         dn='n'
         else
         dn='d'
      end if
      if(dn.eq.'n') then
        do i=1,nlayermx                                    
	      jtot(i)=0.
        enddo       
        return
      endif 

      !initializations
      jergs(:,:,:)=0.
      xabsi(:,:)=0.
      jtot(:)=0.
      !All number densities to a single array, xabsi(species,layer)
      do i=1,nlayermx
         xabsi(1,i)  = rm(i,i_co2)
         xabsi(2,i)  = rm(i,i_o2)
         xabsi(3,i)  = rm(i,i_o)
         xabsi(4,i)  = rm(i,i_h2o)
         xabsi(5,i)  = rm(i,i_h2)
         xabsi(6,i)  = rm(i,i_h2o2)
         !Only if O3, N or ion chemistry requested
         if(euvmod.ge.1) then
            xabsi(7,i)  = rm(i,i_o3)
         endif
         !Only if N or ion chemistry requested
         if(euvmod.ge.2) then
            xabsi(8,i)  = rm(i,i_n2)
            xabsi(9,i)  = rm(i,i_n)
            xabsi(10,i) = rm(i,i_no)
            xabsi(13,i) = rm(i,i_no2)
         endif
         xabsi(11,i) = rm(i,i_co)
         xabsi(12,i) = rm(i,i_h)
      end do

      !Calculation of photoabsortion coefficient
      if(solvarmod.eq.0) then
         call jthermcalc(ig,euvmod,rm,nespeuv,tx,iz,zenit)
      else if (solvarmod.eq.1) then
         call jthermcalc_e107(ig,euvmod,rm,nespeuv,tx,iz,zenit,zday)
         do indexint=1,ninter
            fluxtop(indexint)=1.
         enddo
      endif

      !Total photoabsorption coefficient
      do i=1,nlayermx
         jtot(i)=0.
        do j=1,nabs
          do indexint=1,ninter
            jergs(indexint,j,i) = jfotsout(indexint,j,i) 
     $              * xabsi (j,i) * fluxtop(indexint)  
     $              / (0.5e9 * freccen(indexint))
            jtot(i)=jtot(i)+jergs(indexint,j,i)
          end do
        end do
      end do

      return

      end

