      SUBROUTINE euvheat(pt,pdt,pplev,pplay,zzlay, &
           mu0,ptimestep,ptime,zday,pq,pdq,pdteuv)

      IMPLICIT NONE
!=======================================================================
!   subject:
!   --------
!   Computing heating rate due to EUV absorption
!
!   author:  MAC 2002
!   ------
!
!   input:
!   ----- 
!   mu0(ngridmx)           
!   pplay(ngrid,nlayer)   pressure at middle of layers (Pa)
!
!   output:
!   -------
!
!   pdteuv(ngrid,nlayer)      Heating rate (K/s)
!
!=======================================================================
!
!    0.  Declarations :
!    ------------------
!
#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "comdiurn.h"
#include "param.h"
#include "param_v4.h"
#include "chimiedata.h"
#include "tracer.h"
#include "conc.h"
!-----------------------------------------------------------------------
!    Input/Output
!    ------------

      real :: pt(ngridmx,nlayermx)
      real :: pdt(ngridmx,nlayermx)
      real :: pplev(ngridmx,nlayermx+1)
      real :: pplay(ngridmx,nlayermx)
      real :: zzlay(ngridmx,nlayermx)
      real :: mu0(ngridmx)
      real :: ptimestep,ptime
      real :: zday
      real :: pq(ngridmx,nlayermx,nqmx)
      real :: pdq(ngridmx,nlayermx,nqmx)

      real :: pdteuv(ngridmx,nlayermx)
!
!    Local variables :
!    -----------------
      integer,save :: nespeuv    ! Number of species considered

      INTEGER :: l,ig,n
      integer,save :: euvmod
      real, allocatable :: rm(:,:)   ! number density (cm-3)
      real :: zq(ngridmx,nlayermx,nqmx) ! local updated tracer quantity
      real :: zt(ngridmx,nlayermx)      ! local updated atmospheric temperature
      real :: zlocal(nlayermx)
      real :: zenit
      real :: jtot(nlayermx)
      real :: dens						! amu/cm-3
      real :: tx(nlayermx)
!      real euveff     !UV heating efficiency
      
! tracer indexes for the EUV heating:
!!! ATTENTION. These values have to be identical to those in chemthermos.F90
!!! If the values are changed there, the same has to be done here  !!!
      integer,parameter :: i_co2=1
      integer,parameter :: i_o2=2
      integer,parameter :: i_o=3
      integer,parameter :: i_co=4
      integer,parameter :: i_h=5
      integer,parameter :: i_oh=6
      integer,parameter :: i_ho2=7
      integer,parameter :: i_h2=8
      integer,parameter :: i_h2o=9
      integer,parameter :: i_h2o2=10
      integer,parameter :: i_o1d=11
      integer,parameter :: i_o3=12
      integer,parameter :: i_n2=13
      integer,parameter :: i_n=14
      integer,parameter :: i_no=15
      integer,parameter :: i_n2d=16
      integer,parameter :: i_no2=17

      
! Tracer indexes in the GCM:
      integer,save :: g_co2=0
      integer,save :: g_o=0
      integer,save :: g_o2=0
      integer,save :: g_h2=0
      integer,save :: g_h2o2=0
      integer,save :: g_h2o=0
      integer,save :: g_o3=0
      integer,save :: g_n2=0
      integer,save :: g_n=0
      integer,save :: g_no=0
      integer,save :: g_co=0
      integer,save :: g_h=0
      integer,save :: g_no2=0
      integer,save :: g_oh=0
      integer,save :: g_ho2=0
      integer,save :: g_o1d=0
      integer,save :: g_n2d=0


      logical,save :: firstcall=.true.

! Initializations and sanity checks:


      if (firstcall) then
         nespeuv=0
        ! identify the indexes of the tracers we'll need
         g_co2=igcm_co2
         if (g_co2.eq.0) then
            write(*,*) "euvheat: Error; no CO2 tracer !!!"
            write(*,*) "CO2 is always needed if calleuv=.true."
            stop
         else
            nespeuv=nespeuv+1
         endif
         g_o=igcm_o
         if (g_o.eq.0) then
            write(*,*) "euvheat: Error; no O tracer !!!"
            write(*,*) "O is always needed if calleuv=.true."
            stop
         else
            nespeuv=nespeuv+1
         endif
         g_o2=igcm_o2
         if (g_o2.eq.0) then
            write(*,*) "euvheat: Error; no O2 tracer !!!"
            write(*,*) "O2 is always needed if calleuv=.true."
            stop
         else
            nespeuv=nespeuv+1
         endif
         g_h2=igcm_h2
         if (g_h2.eq.0) then
            write(*,*) "euvheat: Error; no H2 tracer !!!"
            write(*,*) "H2 is always needed if calleuv=.true."
            stop
         else
            nespeuv=nespeuv+1
         endif
         g_oh=igcm_oh
         if (g_oh.eq.0) then
            write(*,*) "euvheat: Error; no OH tracer !!!"
            write(*,*) "OH must always be present if thermochem=T"
            stop
         else
            nespeuv=nespeuv+1  
         endif
         g_ho2=igcm_ho2
         if (g_ho2.eq.0) then
            write(*,*) "euvheat: Error; no HO2 tracer !!!"
            write(*,*) "HO2 must always be present if thermochem=T"
            stop
         else
            nespeuv=nespeuv+1  
         endif
         g_h2o2=igcm_h2o2
         if (g_h2o2.eq.0) then
            write(*,*) "euvheat: Error; no H2O2 tracer !!!"
            write(*,*) "H2O2 is always needed if calleuv=.true."
            stop
         else
            nespeuv=nespeuv+1
         endif
         g_h2o=igcm_h2o_vap
         if (g_h2o.eq.0) then
            write(*,*) "euvheat: Error; no water vapor tracer !!!"
            write(*,*) "H2O is always needed if calleuv=.true."
            stop
         else
            nespeuv=nespeuv+1
         endif
         g_o1d=igcm_o1d
         if (g_o1d.eq.0) then
            write(*,*) "euvheat: Error; no O1D tracer !!!"
            write(*,*) "O1D must always be present if thermochem=T"
            stop
         else
            nespeuv=nespeuv+1  
         endif
         g_co=igcm_co
         if (g_co.eq.0) then
            write(*,*) "euvheat: Error; no CO tracer !!!"
            write(*,*) "CO is always needed if calleuv=.true."
            stop
         else
            nespeuv=nespeuv+1
         endif
         g_h=igcm_h
         if (g_h.eq.0) then
            write(*,*) "euvheat: Error; no H tracer !!!"
            write(*,*) "H is always needed if calleuv=.true."
            stop
         else
            nespeuv=nespeuv+1
         endif
         
         euvmod = 0             !Default: C/O/H chemistry 
         !Check if O3 is present
         g_o3=igcm_o3
         if (g_o3.eq.0) then
            write(*,*) "euvheat: Error; no O3 tracer !!!"
            write(*,*) "O3 must be present if calleuv=.true."
            stop
         else 
            nespeuv=nespeuv+1
            euvmod=1
         endif

         !Nitrogen species
         !NO is used to determine if N chemistry is wanted
         !euvmod=2 -> N chemistry
         g_no=igcm_no
         if (g_no.eq.0) then
            write(*,*) "euvheat: no NO tracer"
            write(*,*) "No N species in UV heating"
         else if(g_no.ne.0) then
            nespeuv=nespeuv+1
            euvmod=2
         endif
         ! n2
         g_n2=igcm_n2
         if(euvmod.eq.2) then
            if (g_n2.eq.0) then
               write(*,*) "euvheat: Error; no N2 tracer !!!"
               write(*,*) "N2 needed if NO is in traceur.def"
               stop
            else
               nespeuv=nespeuv+1
            endif
         endif  ! Of if(euvmod.eq.2)
         ! N
         g_n=igcm_n
         if(euvmod == 2) then
            if (g_n.eq.0) then
               write(*,*) "euvheat: Error; no N tracer !!!"
               write(*,*) "N needed if NO is in traceur.def"
               stop
            else if(g_n.ne.0) then
               nespeuv=nespeuv+1
            endif
         else
            if(g_n /= 0) then
               write(*,*) "euvheat: Error: N present, but NO not!!!"
               write(*,*) "Both must be in traceur.def"
               stop
            endif
         endif   !Of if(euvmod==2)
         !NO2
         g_no2=igcm_no2
         if(euvmod == 2) then
            if (g_no2.eq.0) then
               write(*,*) "euvheat: Error; no NO2 tracer !!!"
               write(*,*) "NO2 needed if NO is in traceur.def"
               stop
            else if(g_no2.ne.0) then
               nespeuv=nespeuv+1
            endif
         else
            if(g_no2 /= 0) then
               write(*,*) "euvheat: Error: NO2 present, but NO not!!!"
               write(*,*) "Both must be in traceur.def"
               stop
            endif
         endif   !Of if(euvmod==2)
         !N2D
         g_n2d=igcm_n2d
         if(euvmod == 2) then
            if (g_n2d.eq.0) then
               write(*,*) "euvheat: Error; no N2D tracer !!!"
               write(*,*) "N2D needed if NO is in traceur.def"
               stop
            else
               nespeuv=nespeuv+1  
            endif
         else
            if(g_n2d /= 0) then
               write(*,*) "euvheat: Error: N2D present, but NO not!!!"
               write(*,*) "Both must be in traceur.def"
               stop
            endif
         endif   !Of if(euvmod==2)

         !Check if nespeuv is appropriate for the value of euvmod
         select case(euvmod)
         case(0)
            if(nespeuv.ne.11) then
               write(*,*)'euvheat: Wrong number of tracers!'
               stop
            else
               write(*,*)'euvheat: Computing absorption by',nespeuv, &
                    ' species'
            endif
         case(1)
            if(nespeuv.ne.12) then
               write(*,*)'euvheat: Wrong number of tracers!',nespeuv
               stop
            else
               write(*,*)'euvheat: Computing absorption by',nespeuv,  &
                    ' species'
            endif
         case(2)
            if(nespeuv.ne.17) then
               write(*,*)'euvheat: Wrong number of tracers!'
               stop
            else
               write(*,*)'euvheat: Computing absorption by',nespeuv,  &
                    ' species'
            endif
         end select
        
         firstcall= .false.
      endif                     ! of if (firstcall)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !Number of species if not firstcall


      !Allocate density vector
      allocate(rm(nlayermx,nespeuv))
      ! build local updated values of tracers and temperature
      do l=1,nlayermx
        do ig=1,ngridmx
          ! chemical species
          zq(ig,l,g_co2)=pq(ig,l,g_co2)+pdq(ig,l,g_co2)*ptimestep
          zq(ig,l,g_o2)=pq(ig,l,g_o2)+pdq(ig,l,g_o2)*ptimestep
          zq(ig,l,g_o)=pq(ig,l,g_o)+pdq(ig,l,g_o)*ptimestep
          zq(ig,l,g_h2)=pq(ig,l,g_h2)+pdq(ig,l,g_h2)*ptimestep
          zq(ig,l,g_h2o2)=pq(ig,l,g_h2o2)+pdq(ig,l,g_h2o2)*ptimestep
          zq(ig,l,g_h2o)=pq(ig,l,g_h2o)+pdq(ig,l,g_h2o)*ptimestep
          zq(ig,l,g_n2)=pq(ig,l,g_n2)+pdq(ig,l,g_n2)*ptimestep
          zq(ig,l,g_co)=pq(ig,l,g_co)+pdq(ig,l,g_co)*ptimestep
          zq(ig,l,g_h)=pq(ig,l,g_h)+pdq(ig,l,g_h)*ptimestep
          if(euvmod.ge.1)   &
               zq(ig,l,g_o3)=pq(ig,l,g_o3)+pdq(ig,l,g_o3)*ptimestep
          if(euvmod.eq.2) then
             zq(ig,l,g_n)=pq(ig,l,g_n)+pdq(ig,l,g_n)*ptimestep
             zq(ig,l,g_no)=pq(ig,l,g_no)+pdq(ig,l,g_no)*ptimestep
             zq(ig,l,g_no2)=pq(ig,l,g_no2)+pdq(ig,l,g_no2)*ptimestep
          endif
          if(euvmod.gt.2.or.euvmod.lt.0) then
             write(*,*)'euvheat: bad value for euvmod. Stop'
             stop
          endif
          ! atmospheric temperature
          zt(ig,l)=pt(ig,l)+pdt(ig,l)*ptimestep
        enddo
      enddo
      
      !Solar flux calculation
      if(solvarmod.eq.0) call flujo(solarcondate)

                                         ! Not recommended for long runs 
                                         ! (e.g. to build the MCD), as the
                                         ! solar conditions at the end will 
                                         ! be different to the ones initially
                                         ! set
      
      do ig=1,ngridmx
         zenit=acos(mu0(ig))*180./acos(-1.)
         
         do l=1,nlayermx
            !Conversion to number density
            dens=pplay(ig,l)/(rnew(ig,l)*zt(ig,l)) / 1.66e-21
            rm(l,i_co2)  = zq(ig,l,g_co2)  * dens / mmol(g_co2)
            rm(l,i_o2)   = zq(ig,l,g_o2)   * dens / mmol(g_o2)
            rm(l,i_o)    = zq(ig,l,g_o)    * dens / mmol(g_o)
            rm(l,i_h2)   = zq(ig,l,g_h2)   * dens / mmol(g_h2)
            rm(l,i_h2o)  = zq(ig,l,g_h2o)  * dens / mmol(g_h2o)
            rm(l,i_h2o2) = zq(ig,l,g_h2o2) * dens / mmol(g_h2o2)
            rm(l,i_co)   = zq(ig,l,g_co)   * dens / mmol(g_co)
            rm(l,i_h)    = zq(ig,l,g_h)    * dens / mmol(g_h)
            !Only if O3, N or ion chemistry requested
            if(euvmod.ge.1)   &
                 rm(l,i_o3)   = zq(ig,l,g_o3)   * dens / mmol(g_o3)
            !Only if N or ion chemistry requested
            if(euvmod.ge.2) then
               rm(l,i_n2)   = zq(ig,l,g_n2)   * dens / mmol(g_n2)
               rm(l,i_n)    = zq(ig,l,g_n)    * dens / mmol(g_n)
               rm(l,i_no)   = zq(ig,l,g_no)   * dens / mmol(g_no)          
               rm(l,i_no2)  = zq(ig,l,g_no2)  * dens / mmol(g_no2)
            endif
         enddo

!        zlocal(1)=-log(pplay(ig,1)/pplev(ig,1))
!     &            *Rnew(ig,1)*zt(ig,1)/g
         zlocal(1)=zzlay(ig,1)
         zlocal(1)=zlocal(1)/1000.
         tx(1)=zt(ig,1)

         do l=2,nlayermx
            tx(l)=zt(ig,l)
            zlocal(l)=zzlay(ig,l)/1000.
         enddo
        !Routine to calculate the UV heating
         call hrtherm (ig,euvmod,rm,nespeuv,tx,zlocal,zenit,zday,jtot)

!        euveff=0.16    !UV heating efficiency. Following Fox et al. ASR 1996
                       !should vary between 19% and 23%. Lower values 
                       !(i.e. 16%) can be used to compensate underestimation 
                       !of 15-um cooling (see Forget et al. JGR 2009 and 
                       !Gonzalez-Galindo et al. JGR 2009) for details
        !Calculates the UV heating from the total photoabsorption coefficient
        do l=1,nlayermx
          pdteuv(ig,l)=euveff*jtot(l)/10.                  &
               /(cpnew(ig,l)*pplay(ig,l)/(rnew(ig,l)*zt(ig,l)))
!     &                  *(1.52/dist_sol)**2  !The solar flux calculated in 
                                              !flujo.F is already corrected for
                                              !the actual Mars-Sun distance
        enddo	
      enddo  ! of do ig=1,ngridmx
      !Deallocations
      deallocate(rm)

      return
      end 
