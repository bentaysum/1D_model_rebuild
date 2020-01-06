      subroutine sponge(ucov,vcov,h,pext,dt,mode)

! Sponge routine: Quench ucov, vcov and potential temperature near the
!                 top of the model
! Depending on 'mode' relaxation of variables is towards:
! mode = 0 : h -> h_mean , ucov -> 0 , vcov -> 0
! mode = 1 : h -> h_mean , ucov -> ucov_mean , vcov -> 0
! mode >= 2 : h -> h_mean , ucov -> ucov_mean , vcov -> vcov_mean
! Number of layer over which sponge is applied is 'nsponge' (read from def file)
! Time scale for quenching at top level is given by 'tetasponge' (read from
! def file) and doubles as level indexes decrease.

      implicit none
#include "dimensions.h"
#include "paramet.h"
#include "comdissip.h"
#include "comvert.h"
#include "comgeom2.h"
#include "sponge.h"

! Arguments:
!------------
      real,intent(inout) :: ucov(iip1,jjp1,llm) ! covariant zonal wind
      real,intent(inout) :: vcov(iip1,jjm,llm) ! covariant meridional wind
      real,intent(inout) :: h(iip1,jjp1,llm) ! potential temperature
      real,intent(in) :: pext(iip1,jjp1) ! extensive pressure
      real,intent(in) :: dt   ! time step
      integer,intent(in) :: mode  ! sponge mode

c   Local:
c   ------

      real,save :: sig_s(llm) !sigma au milieu des couches
      REAL vm,um,hm,ptot(jjp1)
      real,save :: cst(llm)

      INTEGER l,i,j
      integer,save :: l0 ! layer down to which sponge is applied

      real ssum

      real echelle,zkm
      logical,save :: firstcall=.true.



      if (firstcall) then

       ! build approximative sigma levels at midlayer
        do l=1,llm
          sig_s(l)=((ap(l)+ap(l+1))/preff+bp(l)+bp(l+1))/2.
        enddo

        l0=llm-nsponge+1
 
        PRINT*
        print*,'sponge mode',mode
        print*,'nsponge tetasponge ',nsponge,tetasponge
        print*,'Coeffs for the sponge layer'
        print*,'Z (km)     tau      cst'
        do l=llm,l0,-1
          ! double time scale with every level, starting from the top
          cst(l)=dt/(tetasponge*2**(llm-l))
        enddo

        echelle=10.
        do l=l0,llm
           zkm=-echelle*log(sig_s(l))
           print*,zkm,dt/cst(l),cst(l)
        enddo
        PRINT*

        firstcall=.false.
      endif ! of if (firstcall)

c-----------------------------------------------------------------------
c   calcul de la dissipation:
c   -------------------------

      do j=1,jjp1
         ptot(j)=ssum(iim,pext(1,j),1)
      enddo
 
c   potential temperature
      do l=l0,llm
         do j=1,jjp1
            hm=0.
            do i=1,iim
               hm=hm+h(i,j,l)*pext(i,j)
            enddo
            hm=hm/ptot(j)
            do i=1,iim
               h(i,j,l)=h(i,j,l)-cst(l)*(h(i,j,l)-hm)
            enddo
            h(iip1,j,l)=h(1,j,l)
         enddo
      enddo

c   zonal wind
      do l=l0,llm
         do j=2,jjm
            um=0.
            if(mode.ge.1) then
               do i=1,iim
                  um=um+0.5*ucov(i,j,l)*(pext(i,j)+pext(i+1,j))
     s               /cu(i,j)
               enddo
               um=um/ptot(j)
            endif
            do i=1,iim
               ucov(i,j,l)=ucov(i,j,l)-cst(l)*(ucov(i,j,l)-um*cu(i,j))
            enddo
            ucov(iip1,j,l)=ucov(1,j,l)
         enddo
      enddo

c  meridional wind
      do l=l0,llm
         do j=1,jjm
            vm=0.
            if(mode.ge.2) then
               do i=1,iim
                  vm=vm+vcov(i,j,l)*(pext(i,j)+pext(i,j+1))
     s               /cv(i,j)
               enddo
               vm=vm/(ptot(j)+ptot(j+1))
            endif
            do i=1,iim
               vcov(i,j,l)=vcov(i,j,l)-cst(l)*(vcov(i,j,l)-vm*cv(i,j))
            enddo
            vcov(iip1,j,l)=vcov(1,j,l)
         enddo
      enddo

      end
