      SUBROUTINE callsedim(ngrid,nlay, ptimestep,
     &                pplev,zlev, zlay, pt, rdust, rice,
     &                rsedcloud,rhocloud,
     &                pq, pdqfi, pdqsed,pdqs_sed,nq, 
     &                tau,tauscaling)
! to use  'getin'
      USE ioipsl_getincom
      USE updaterad
      IMPLICIT NONE

c=======================================================================
c      Sedimentation of the  Martian aerosols
c      depending on their density and radius
c
c      F.Forget 1999
c
c      Modified by J.-B. Madeleine 2010: Now includes the doubleq
c        technique in order to have only one call to callsedim in
c        physiq.F.
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------
      
#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "tracer.h"
#include "callkeys.h"

c
c   arguments:
c   ----------

      INTEGER ngrid              ! number of horizontal grid points
      INTEGER nlay               ! number of atmospheric layers
      REAL ptimestep             ! physics time step (s)
      REAL pplev(ngrid,nlay+1)   ! pressure at inter-layers (Pa)
      REAL pt(ngrid,nlay)        ! temperature at mid-layer (K)
      REAL zlev(ngrid,nlay+1)    ! altitude at layer boundaries
c    Aerosol radius provided by the water ice microphysical scheme:
      REAL rdust(ngrid,nlay)     ! Dust geometric mean radius (m)
      REAL rice(ngrid,nlay)      ! Ice geometric mean radius (m)

c    Traceurs :
      integer nq             ! number of tracers
      real pq(ngrid,nlay,nq)  ! tracers (kg/kg)
      real pdqfi(ngrid,nlay,nq)  ! tendency before sedimentation (kg/kg.s-1)
      real pdqsed(ngrid,nlay,nq) ! tendency due to sedimentation (kg/kg.s-1)
      real pdqs_sed(ngrid,nq)    ! flux at surface (kg.m-2.s-1)
      
c   local:
c   ------

      INTEGER l,ig, iq

      real zqi(ngridmx,nlayermx,nqmx) ! to locally store tracers
      real masse (ngridmx,nlayermx) ! Layer mass (kg.m-2)
      real epaisseur (ngridmx,nlayermx) ! Layer thickness (m)
      real wq(ngridmx,nlayermx+1) ! displaced tracer mass (kg.m-2)
      real r0(ngridmx,nlayermx) ! geometric mean radius used for
                                !   sedimentation (m)
      real r0dust(ngridmx,nlayermx) ! geometric mean radius used for
                                    !   dust (m)
!      real r0ccn(ngridmx,nlayermx)  ! geometric mean radius used for
!                                    !   CCNs (m)
c     Sedimentation radius of water ice
      real rsedcloud(ngridmx,nlayermx)
      real beta      ! correction for the shape of the ice particles (cf. newsedim)
      save beta
c     Cloud density (kg.m-3)
      real rhocloud(ngridmx,nlayermx)
      
c     for ice radius computation
      REAL Mo,No
      REAL tau(ngrid,nlay), tauscaling(ngrid)
      REAL zlay(ngrid,nlay)   ! altitude at the middle of the layers
      REAl ccntyp



c     Discrete size distributions (doubleq)
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       1) Parameters used to represent the changes in fall
c          velocity as a function of particle size;
      integer nr,ir
      parameter (nr=12) !(nr=7) ! number of bins
      real rd(nr),qr(ngridmx,nlayermx,nr)
      real rdi(nr+1)    ! extreme and intermediate radii
      real Sq(ngridmx,nlayermx)
      real rdmin,rdmax,rdimin,rdimax
      data rdmin/1.e-8/ !/1.e-7/
      data rdmax/30.e-6/
      data rdimin/1.e-10/
      data rdimax/1e-4/
      save rd, rdi

c       2) Second size distribution for the log-normal integration
c          (the mass mixing ratio is computed for each radius)

      integer ninter, iint
      parameter (ninter=4)  ! nombre de point entre chaque rayon rdi
      real rr(ninter,nr)
      save rr
      integer radpower
      real sigma0

c       3) Other local variables used in doubleq

      INTEGER idust_mass   ! index of tracer containing dust mass
                           !   mix. ratio
      INTEGER idust_number ! index of tracer containing dust number
                           !   mix. ratio
      INTEGER iccn_mass    ! index of tracer containing CCN mass
                           !   mix. ratio
      INTEGER iccn_number  ! index of tracer containing CCN number
                           !   mix. ratio
      SAVE idust_mass,idust_number
      SAVE iccn_mass,iccn_number
      

c     Firstcall:

      LOGICAL firstcall
      SAVE firstcall
      DATA firstcall/.true./

c    ** un petit test de coherence
c       --------------------------

      IF (firstcall) THEN
         
         IF(ngrid.NE.ngridmx) THEN
            PRINT*,'STOP dans callsedim'
            PRINT*,'probleme de dimensions :'
            PRINT*,'ngrid  =',ngrid
            PRINT*,'ngridmx  =',ngridmx
            STOP
         ENDIF

c       Doubleq: initialization
        IF (doubleq) THEN
         do ir=1,nr
             rd(ir)= rdmin*(rdmax/rdmin)**(float(ir-1)/float(nr-1))
         end do
         rdi(1)=rdimin
         do ir=2,nr
           rdi(ir)= sqrt(rd(ir-1)*rd(ir))
         end do
         rdi(nr+1)=rdimax

         do ir=1,nr
           do iint=1,ninter
             rr(iint,ir)=
     &        rdi(ir)*
     &        (rdi(ir+1)/rdi(ir))**(float(iint-1)/float(ninter-1))
c             write(*,*) rr(iint,ir)
           end do
         end do

      ! identify tracers corresponding to mass mixing ratio and
      ! number mixing ratio
        idust_mass=0      ! dummy initialization
        idust_number=0    ! dummy initialization

        do iq=1,nq
          if (noms(iq).eq."dust_mass") then
            idust_mass=iq
          endif
          if (noms(iq).eq."dust_number") then
            idust_number=iq
          endif
        enddo

        ! check that we did find the tracers
        if ((idust_mass.eq.0).or.(idust_number.eq.0)) then
          write(*,*) 'callsedim: error! could not identify'
          write(*,*) ' tracers for dust mass and number mixing'
          write(*,*) ' ratio and doubleq is activated!'
          stop
        endif
        ENDIF !of if (doubleq)

        IF (microphys) THEN
          iccn_mass=0
          iccn_number=0
          do iq=1,nq
            if (noms(iq).eq."ccn_mass") then
              iccn_mass=iq
            endif
            if (noms(iq).eq."ccn_number") then
              iccn_number=iq
            endif
          enddo
          ! check that we did find the tracers
          if ((iccn_mass.eq.0).or.(iccn_number.eq.0)) then
            write(*,*) 'callsedim: error! could not identify'
            write(*,*) ' tracers for ccn mass and number mixing'
            write(*,*) ' ratio and microphys is activated!'
            stop
          endif
        ENDIF !of if (microphys)

        IF (water) THEN
         write(*,*) "correction for the shape of the ice particles ?"
         beta=0.75 ! default value
         call getin("ice_shape",beta)
         write(*,*) " ice_shape = ",beta

          write(*,*) "water_param nueff Sedimentation:", nuice_sed
          IF (activice) THEN
            write(*,*) "water_param nueff Radiative:", nuice_ref
          ENDIF
        ENDIF
      
        firstcall=.false.
      ENDIF ! of IF (firstcall)

c-----------------------------------------------------------------------
c    1. Initialization
c    -----------------

      zqi(1:ngrid,1:nlay,1:nqmx) = 0.
c     Updating the mass mixing ratio with the tendencies coming
c       from other parameterizations:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do iq=1,nq
        do l=1,nlay
          do ig=1,ngrid
            zqi(ig,l,iq)=pq(ig,l,iq)+pdqfi(ig,l,iq)*ptimestep
          enddo
        enddo 
      enddo

c    Computing the different layer properties
c    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c    Mass (kg.m-2), thickness(m), crossing time (s)  etc.

      do  l=1,nlay
        do ig=1, ngrid
          masse(ig,l)=(pplev(ig,l) - pplev(ig,l+1)) /g 
          epaisseur(ig,l)= zlev(ig,l+1) - zlev(ig,l)
        end do
      end do

c =================================================================
c     Compute the geometric mean radius used for sedimentation

      if (doubleq) then
        do l=1,nlay
          do ig=1, ngrid
     
         call updaterdust(zqi(ig,l,igcm_dust_mass),
     &                    zqi(ig,l,igcm_dust_number),r0dust(ig,l),
     &                    tauscaling(ig))
          
          end do
        end do
      endif


c =================================================================
      do iq=1,nq
        if(radius(iq).gt.1.e-9) then   ! no sedim for gaz

c -----------------------------------------------------------------
c         DOUBLEQ CASE
c -----------------------------------------------------------------
          if ((doubleq.and.
     &        ((iq.eq.idust_mass).or.
     &         (iq.eq.idust_number)))) then
     
c           Computing size distribution:
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c            if ((iq.eq.idust_mass).or.(iq.eq.idust_number)) then
              do  l=1,nlay
                do ig=1, ngrid
                  r0(ig,l)=r0dust(ig,l)
                end do
              end do
              sigma0 = varian

c        Computing mass mixing ratio for each particle size
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF ((iq.EQ.idust_mass).or.(iq.EQ.iccn_mass)) then
            radpower = 2
          ELSE  ! number
            radpower = -1
          ENDIF
          Sq(1:ngrid,1:nlay) = 0.
          do ir=1,nr
            do l=1,nlay
              do ig=1,ngrid
c                ****************
c                Size distribution integration
c                (Trapezoid Integration Method)
                 qr(ig,l,ir)=0.5*(rr(2,ir)-rr(1,ir))*
     &             (rr(1,ir)**radpower)*
     &             exp(-(log(rr(1,ir)/r0(ig,l)))**2/(2*sigma0**2))
                 do iint=2,ninter-1
                   qr(ig,l,ir)=qr(ig,l,ir) +
     &             0.5*(rr(iint+1,ir)-rr(iint-1,ir))*
     &             (rr(iint,ir)**radpower)*
     &             exp(-(log(rr(iint,ir)/r0(ig,l)))**2/
     &             (2*sigma0**2))
                 end do
                 qr(ig,l,ir)=qr(ig,l,ir) +
     &             0.5*(rr(ninter,ir)-rr(ninter-1,ir))*
     &             (rr(ninter,ir)**radpower)*
     &             exp(-(log(rr(ninter,ir)/r0(ig,l)))**2/
     &             (2*sigma0**2))

c                **************** old method (not recommended!)
c                qr(ig,l,ir)=(rd(ir)**(5-3*iq))*
c    &           exp( -(log(rd(ir)/r0(ig,l)))**2 / (2*sigma0**2) )
c                ******************************

                 Sq(ig,l)=Sq(ig,l)+qr(ig,l,ir)
              enddo
            enddo
          enddo

          do ir=1,nr
            do l=1,nlay
              do ig=1,ngrid
                 qr(ig,l,ir) = zqi(ig,l,iq)*qr(ig,l,ir)/Sq(ig,l)
              enddo
            enddo
          enddo

c         Computing sedimentation for each tracer
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          zqi(1:ngrid,1:nlay,iq) = 0.
          pdqs_sed(1:ngrid,iq) = 0.

          do ir=1,nr
          
               call newsedim(ngrid,nlay,1,1,ptimestep,
     &         pplev,masse,epaisseur,pt,rd(ir),rho_dust,qr(1,1,ir),
     &         wq,0.5)

c            Tendencies
c            ~~~~~~~~~~
             do ig=1,ngrid
               pdqs_sed(ig,iq) = pdqs_sed(ig,iq)
     &                                + wq(ig,1)/ptimestep
             end do
             DO l = 1, nlay
               DO ig=1,ngrid
                 zqi(ig,l,iq)=zqi(ig,l,iq)+qr(ig,l,ir)
               ENDDO
             ENDDO
          enddo ! of do ir=1,nr
c -----------------------------------------------------------------
c         WATER CYCLE CASE
c -----------------------------------------------------------------
           else if ((iq .eq. iccn_mass) .or. (iq .eq. iccn_number)
     &       .or. (iq .eq. igcm_h2o_ice)) then
            if (microphys) then
              ! water ice sedimentation
              call newsedim(ngrid,nlay,ngrid*nlay,ngrid*nlay,
     &        ptimestep,pplev,masse,epaisseur,pt,rsedcloud,rhocloud,
     &        zqi(1,1,iq),wq,beta)
            else
              ! water ice sedimentation
              call newsedim(ngrid,nlay,ngrid*nlay,1,
     &        ptimestep,pplev,masse,epaisseur,pt,rsedcloud,rho_q(iq),
     &        zqi(1,1,iq),wq,beta)
            endif ! of if (microphys)
c           Tendencies
c           ~~~~~~~~~~
            do ig=1,ngrid 
              pdqs_sed(ig,iq)=wq(ig,1)/ptimestep
            end do
c -----------------------------------------------------------------
c         GENERAL CASE
c -----------------------------------------------------------------
          else
            call newsedim(ngrid,nlay,1,1,ptimestep,
     &      pplev,masse,epaisseur,pt,radius(iq),rho_q(iq),
     &      zqi(1,1,iq),wq,1.0)
c           Tendencies
c           ~~~~~~~~~~
            do ig=1,ngrid 
              pdqs_sed(ig,iq)=wq(ig,1)/ptimestep
            end do
          endif ! of if doubleq and if water
c -----------------------------------------------------------------

c         Compute the final tendency:
c         ---------------------------
          DO l = 1, nlay
            DO ig=1,ngrid
              pdqsed(ig,l,iq)=(zqi(ig,l,iq)-
     $        (pq(ig,l,iq) + pdqfi(ig,l,iq)*ptimestep))/ptimestep
            ENDDO
          ENDDO

        endif ! of if(radius(iq).gt.1.e-9)
c =================================================================
      enddo ! of do iq=1,nq
 
c     Update the dust particle size "rdust"
c     -------------------------------------
      if (doubleq) then
       DO l = 1, nlay
        DO ig=1,ngrid
        
     
         call updaterdust(zqi(ig,l,igcm_dust_mass),
     &                    zqi(ig,l,igcm_dust_number),rdust(ig,l),
     &                    tauscaling(ig))     

          
        ENDDO
       ENDDO
      endif ! of if (doubleq)
      
c     Update the ice particle size "rice"
c     -------------------------------------
      if (water) then
       IF(microphys) THEN 
       
       
        DO l = 1, nlay
          DO ig=1,ngrid

         call updaterice_micro(zqi(ig,l,igcm_h2o_ice),
     &    zqi(ig,l,igcm_ccn_mass),zqi(ig,l,igcm_ccn_number),
     &    tauscaling(ig),rice(ig,l),rhocloud(ig,l))
           
          ENDDO
        ENDDO
        
       ELSE
       
        DO l = 1, nlay
          DO ig=1,ngrid
          
            call updaterice_typ(zqi(ig,l,igcm_h2o_ice),
     &                      tau(ig,1),zlay(ig,l),rice(ig,l))

          ENDDO
        ENDDO
       ENDIF ! of IF(microphys)
      endif ! of if (water)

      RETURN
      END

