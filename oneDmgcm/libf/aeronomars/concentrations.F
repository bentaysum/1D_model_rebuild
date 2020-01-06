      SUBROUTINE concentrations(pplay,pt,pdt,pq,pdq,ptimestep)
                                             
      implicit none

!=======================================================================
! CALCULATION OF MEAN MOLECULAR MASS, Cp, Akk and R
!
! mmean(ngridmx,nlayermx)	amu
! cpnew(ngridmx,nlayermx)	J/kg/K
! rnew(ngridmx,nlayermx)	J/kg/K
! akknew(ngridmx,nlayermx)	coefficient of thermal concduction
!
! version: April 2012 - Franck Lefevre
!=======================================================================

!     declarations
 
#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "comdiurn.h"
#include "chimiedata.h"
#include "tracer.h"
#include "conc.h"

!     input/output

      real pplay(ngridmx,nlayermx)
      real pt(ngridmx,nlayermx)
      real pdt(ngridmx,nlayermx)
      real pq(ngridmx,nlayermx,nqmx)
      real pdq(ngridmx,nlayermx,nqmx)
      real ptimestep

!     local variables

      integer       :: i, l, ig, iq
      integer, save :: nbq, niq(nqmx)
      real          :: ni(nqmx), ntot
      real          :: zq(ngridmx, nlayermx, nqmx)
      real          :: zt(ngridmx, nlayermx)
      real, save    :: aki(nqmx)
      real, save    :: cpi(nqmx)

      logical, save :: firstcall = .true.

      if (firstcall) then

!        find index of chemical tracers to use
!        initialize thermal conductivity and specific heat coefficients
!        !? values are estimated

         nbq = 0 ! to count number of tracers used in this subroutine

         if (igcm_co2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_co2
            aki(nbq) = 3.072e-4
            cpi(nbq) = 0.834e3
         end if
         if (igcm_co /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_co
            aki(nbq) = 4.87e-4
            cpi(nbq) = 1.034e3
         end if
         if (igcm_o /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o
            aki(nbq) = 7.59e-4
            cpi(nbq) = 1.3e3
         end if
         if (igcm_o1d /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o1d
            aki(nbq) = 7.59e-4  !?
            cpi(nbq) = 1.3e3    !?
         end if
         if (igcm_o2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o2
            aki(nbq) = 5.68e-4
            cpi(nbq) = 0.9194e3
         end if
         if (igcm_o3 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o3
            aki(nbq) = 3.00e-4  !?
            cpi(nbq) = 0.800e3  !?
         end if
         if (igcm_h /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_h
            aki(nbq) = 0.0
            cpi(nbq) = 20.780e3
         end if
         if (igcm_h2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_h2
            aki(nbq) = 36.314e-4
            cpi(nbq) = 14.266e3
         end if
         if (igcm_oh /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_oh
            aki(nbq)  = 7.00e-4 !?
            cpi(nbq)  = 1.045e3
         end if
         if (igcm_ho2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ho2
            aki(nbq) = 0.0
            cpi(nbq) = 1.065e3  !?
         end if
         if (igcm_n2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_n2
            aki(nbq) = 5.6e-4
            cpi(nbq) = 1.034e3
         end if
         if (igcm_ar /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ar
            aki(nbq) = 0.0      !?
            cpi(nbq) = 1.000e3  !?
         end if
         if (igcm_h2o_vap /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_h2o_vap
            aki(nbq) = 0.0
            cpi(nbq) = 1.870e3
         end if
         if (igcm_n /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_n
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_no /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_no
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_no2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_no2
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_n2d /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_n2d
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_co2plus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_co2plus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_oplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_oplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_o2plus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o2plus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_coplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_coplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_cplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_cplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_nplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_nplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_noplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_noplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_n2plus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_n2plus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_hplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_hplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_hco2plus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_hco2plus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         
            

         firstcall = .false.

      end if ! if (firstcall)

!     update temperature

      do l = 1,nlayermx
         do ig = 1,ngridmx
            zt(ig,l) = pt(ig,l) + pdt(ig,l)*ptimestep
         end do
      end do

!     update tracers

      do l = 1,nlayermx
         do ig = 1,ngridmx
            do i = 1,nbq
               iq = niq(i) 
               zq(ig,l,iq) = max(1.e-30, pq(ig,l,iq)
     $                                 + pdq(ig,l,iq)*ptimestep)
            end do
         end do
      end do

!     mmean : mean molecular mass
!     rnew  : specific gas constant

      mmean(:,:)  = 0.

      do l = 1,nlayermx
         do ig = 1,ngridmx
            do i = 1,nbq
               iq = niq(i) 
               mmean(ig,l) = mmean(ig,l) + zq(ig,l,iq)/mmol(iq)
            end do
            mmean(ig,l) = 1./mmean(ig,l)
            rnew(ig,l) = 8.314/mmean(ig,l)*1.e3     ! J/kg/K		
         end do
      end do

!     cpnew  : specicic heat
!     akknew : thermal conductivity cofficient
      
      cpnew(:,:)  = 0.
      akknew(:,:) = 0.

      do l = 1,nlayermx
         do ig = 1,ngridmx
            ntot = pplay(ig,l)/(1.381e-23*zt(ig,l))*1.e-6  ! in #/cm3
            do i = 1,nbq
               iq = niq(i) 
               ni(iq) = ntot*zq(ig,l,iq)*mmean(ig,l)/mmol(iq)
               cpnew(ig,l) = cpnew(ig,l) + ni(iq)*cpi(iq)
               akknew(ig,l) = akknew(ig,l) + ni(iq)*aki(iq)
            end do 
            cpnew(ig,l) = cpnew(ig,l)/ntot
            akknew(ig,l) = akknew(ig,l)/ntot
         end do
!        print*, l, mmean(1,l), cpnew(1,l), rnew(1,l)
      end do

      return
      end 
