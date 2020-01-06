       SUBROUTINE watercloud(ngrid,nlay,ptimestep, 
     &                pplev,pplay,pdpsrf,pzlay,pt,pdt,
     &                pq,pdq,pdqcloud,pdtcloud,
     &                nq,tau,tauscaling,rdust,rice,nuice,
     &                rsedcloud,rhocloud)
! to use  'getin'
      USE ioipsl_getincom
      USE updaterad
      IMPLICIT NONE


c=======================================================================
c  Water-ice cloud formation
c  
c  Includes two different schemes:
c    - A simplified scheme (see simpleclouds.F)
c    - An improved microphysical scheme (see improvedclouds.F)
c
c  There is a time loop specific to cloud formation 
c  due to timescales smaller than the GCM integration timestep.
c
c  Authors: Franck Montmessin, Francois Forget, Ehouarn Millour, 
c           J.-B. Madeleine, Thomas Navarro
c
c  2004 - 2012
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "tracer.h"
#include "comgeomfi.h"
#include "dimradmars.h"

c   Inputs:
c   ------

      INTEGER ngrid,nlay
      INTEGER nq                 ! nombre de traceurs 
      REAL ptimestep             ! pas de temps physique (s)
      REAL pplev(ngrid,nlay+1)   ! pression aux inter-couches (Pa)
      REAL pplay(ngrid,nlay)     ! pression au milieu des couches (Pa)
      REAL pdpsrf(ngrid)         ! tendence surf pressure
      REAL pzlay(ngrid,nlay)     ! altitude at the middle of the layers
      REAL pt(ngrid,nlay)        ! temperature at the middle of the layers (K)
      REAL pdt(ngrid,nlay)       ! tendence temperature des autres param.

      real pq(ngrid,nlay,nq)     ! traceur (kg/kg)
      real pdq(ngrid,nlay,nq)    ! tendence avant condensation  (kg/kg.s-1)

      REAL tau(ngridmx,naerkind) ! Column dust optical depth at each point
      REAL tauscaling(ngridmx)   ! Convertion factor for dust amount
      real rdust(ngridmx,nlay)   ! Dust geometric mean radius (m)

c   Outputs:
c   -------

      real pdqcloud(ngrid,nlay,nq) ! tendence de la condensation H2O(kg/kg.s-1)
      REAL pdtcloud(ngrid,nlay)    ! tendence temperature due
                                   ! a la chaleur latente

      REAL rice(ngrid,nlay)    ! Ice mass mean radius (m)
                               ! (r_c in montmessin_2004)
      REAL nuice(ngrid,nlay)   ! Estimated effective variance
                               !   of the size distribution
      real rsedcloud(ngridmx,nlay) ! Cloud sedimentation radius
      real rhocloud(ngridmx,nlay)  ! Cloud density (kg.m-3)

c   local:
c   ------
      
      ! for ice radius computation
      REAL Mo,No
      REAl ccntyp
      
      ! for time loop
      INTEGER microstep  ! time subsampling step variable
      INTEGER imicro     ! time subsampling for coupled water microphysics & sedimentation
      SAVE imicro
      REAL microtimestep ! integration timestep for coupled water microphysics & sedimentation
      SAVE microtimestep
      
      ! tendency given by clouds (inside the micro loop)
      REAL subpdqcloud(ngrid,nlay,nq) ! cf. pdqcloud
      REAL subpdtcloud(ngrid,nlay)    ! cf. pdtcloud

      ! global tendency (clouds+physics)
      REAL subpdq(ngrid,nlay,nq)      ! cf. pdqcloud
      REAL subpdt(ngrid,nlay)         ! cf. pdtcloud


      INTEGER iq,ig,l
      LOGICAL,SAVE :: firstcall=.true.

c    ** un petit test de coherence
c       --------------------------

      IF (firstcall) THEN
        IF(ngrid.NE.ngridmx) THEN
            PRINT*,'STOP dans watercloud'
            PRINT*,'probleme de dimensions :'
            PRINT*,'ngrid  =',ngrid
            PRINT*,'ngridmx  =',ngridmx
            STOP
        ENDIF
         
        if (nq.gt.nqmx) then
           write(*,*) 'stop in watercloud (nq.gt.nqmx)!'
           write(*,*) 'nq=',nq,' nqmx=',nqmx
           stop
        endif
         
        write(*,*) "watercloud: igcm_h2o_vap=",igcm_h2o_vap
        write(*,*) "            igcm_h2o_ice=",igcm_h2o_ice
                
        write(*,*) "time subsampling for microphysic ?"
#ifdef MESOSCALE
        imicro = 2
#else
        imicro = 15
#endif
        call getin("imicro",imicro)
        write(*,*)"imicro = ",imicro
        
        microtimestep = ptimestep/real(imicro)
        write(*,*)"Physical timestep is",ptimestep 
        write(*,*)"Microphysics timestep is",microtimestep 

        firstcall=.false.
      ENDIF ! of IF (firstcall)
      
c-----Initialization
      subpdq(1:ngrid,1:nlay,1:nq) = 0
      subpdt(1:ngrid,1:nlay)      = 0
      
      ! default value if no ice
      rhocloud(1:ngrid,1:nlay) = rho_dust



c------------------------------------------------------------------
c Time subsampling for microphysics 
c------------------------------------------------------------------
      DO microstep=1,imicro 
      
c-------------------------------------------------------------------
c   1.  Tendencies: 
c------------------


c------ Temperature tendency subpdt
        ! Each microtimestep we give the cloud scheme a stepped entry subpdt instead of pdt
        ! If imicro=1 subpdt is the same as pdt
        DO l=1,nlay
          DO ig=1,ngrid
             subpdt(ig,l) = subpdt(ig,l)
     &        + pdt(ig,l) ! At each micro timestep we add pdt in order to have a stepped entry
          ENDDO
        ENDDO
c------ Tracers tendencies subpdq
c------ At each micro timestep we add pdq in order to have a stepped entry
        IF (microphys) THEN
          DO l=1,nlay
            DO ig=1,ngrid
              subpdq(ig,l,igcm_dust_mass) = 
     &            subpdq(ig,l,igcm_dust_mass)
     &          + pdq(ig,l,igcm_dust_mass)
              subpdq(ig,l,igcm_dust_number) = 
     &            subpdq(ig,l,igcm_dust_number)
     &          + pdq(ig,l,igcm_dust_number)
              subpdq(ig,l,igcm_ccn_mass) = 
     &            subpdq(ig,l,igcm_ccn_mass)
     &          + pdq(ig,l,igcm_ccn_mass)
              subpdq(ig,l,igcm_ccn_number) = 
     &            subpdq(ig,l,igcm_ccn_number)
     &          + pdq(ig,l,igcm_ccn_number)
            ENDDO
          ENDDO
        ENDIF
        DO l=1,nlay
          DO ig=1,ngrid
            subpdq(ig,l,igcm_h2o_ice) = 
     &          subpdq(ig,l,igcm_h2o_ice)
     &        + pdq(ig,l,igcm_h2o_ice)
            subpdq(ig,l,igcm_h2o_vap) = 
     &          subpdq(ig,l,igcm_h2o_vap)
     &        + pdq(ig,l,igcm_h2o_vap)
          ENDDO
        ENDDO
        
        
c-------------------------------------------------------------------
c   2.  Main call to the different cloud schemes:
c------------------------------------------------
        IF (microphys) THEN
           CALL improvedclouds(ngrid,nlay,microtimestep,
     &             pplay,pt,subpdt,
     &             pq,subpdq,subpdqcloud,subpdtcloud,
     &             nq,tauscaling)

        ELSE
           CALL simpleclouds(ngrid,nlay,microtimestep,
     &             pplay,pzlay,pt,subpdt,
     &             pq,subpdq,subpdqcloud,subpdtcloud,
     &             nq,tau,rice)
        ENDIF


c-------------------------------------------------------------------
c   3.  Updating tendencies after cloud scheme:
c-----------------------------------------------

        IF (microphys) THEN
          DO l=1,nlay
            DO ig=1,ngrid
              subpdq(ig,l,igcm_dust_mass) =
     &            subpdq(ig,l,igcm_dust_mass)
     &          + subpdqcloud(ig,l,igcm_dust_mass)
              subpdq(ig,l,igcm_dust_number) =
     &            subpdq(ig,l,igcm_dust_number)
     &          + subpdqcloud(ig,l,igcm_dust_number)
              subpdq(ig,l,igcm_ccn_mass) =
     &            subpdq(ig,l,igcm_ccn_mass)
     &          + subpdqcloud(ig,l,igcm_ccn_mass)
              subpdq(ig,l,igcm_ccn_number) =
     &            subpdq(ig,l,igcm_ccn_number)
     &          + subpdqcloud(ig,l,igcm_ccn_number)
            ENDDO
          ENDDO
        ENDIF
        DO l=1,nlay
          DO ig=1,ngrid
            subpdq(ig,l,igcm_h2o_ice) =
     &          subpdq(ig,l,igcm_h2o_ice)
     &        + subpdqcloud(ig,l,igcm_h2o_ice)
            subpdq(ig,l,igcm_h2o_vap) =
     &          subpdq(ig,l,igcm_h2o_vap)
     &        + subpdqcloud(ig,l,igcm_h2o_vap)
          ENDDO
        ENDDO
     
 
      ENDDO ! of DO microstep=1,imicro
      
c-------------------------------------------------------------------
c   6.  Compute final tendencies after time loop:
c------------------------------------------------

c------ Temperature tendency pdtcloud
       DO l=1,nlay
         DO ig=1,ngrid
             pdtcloud(ig,l) =
     &         subpdt(ig,l)/real(imicro)-pdt(ig,l)
          ENDDO
       ENDDO
       
c------ Tracers tendencies pdqcloud
       DO l=1,nlay
         DO ig=1,ngrid
            pdqcloud(ig,l,igcm_h2o_ice) = 
     &        subpdq(ig,l,igcm_h2o_ice)/real(imicro) 
     &       - pdq(ig,l,igcm_h2o_ice)
            pdqcloud(ig,l,igcm_h2o_vap) = 
     &        subpdq(ig,l,igcm_h2o_vap)/real(imicro) 
     &       - pdq(ig,l,igcm_h2o_vap)
         ENDDO
       ENDDO
       
       IF(microphys) THEN
        DO l=1,nlay
         DO ig=1,ngrid
            pdqcloud(ig,l,igcm_ccn_mass) = 
     &        subpdq(ig,l,igcm_ccn_mass)/real(imicro) 
     &       - pdq(ig,l,igcm_ccn_mass)
            pdqcloud(ig,l,igcm_ccn_number) = 
     &        subpdq(ig,l,igcm_ccn_number)/real(imicro) 
     &       - pdq(ig,l,igcm_ccn_number)
         ENDDO
        ENDDO
       ENDIF
       
       IF(scavenging) THEN
        DO l=1,nlay
         DO ig=1,ngrid
            pdqcloud(ig,l,igcm_dust_mass) = 
     &        subpdq(ig,l,igcm_dust_mass)/real(imicro) 
     &       - pdq(ig,l,igcm_dust_mass)
            pdqcloud(ig,l,igcm_dust_number) = 
     &        subpdq(ig,l,igcm_dust_number)/real(imicro)
     &       - pdq(ig,l,igcm_dust_number)
         ENDDO
        ENDDO
       ENDIF

c------- Due to stepped entry, other processes tendencies can add up to negative values
c------- Therefore, enforce positive values and conserve mass


       IF(microphys) THEN
        DO l=1,nlay
         DO ig=1,ngrid
          IF ((pq(ig,l,igcm_ccn_number) + 
     &      ptimestep* (pdq(ig,l,igcm_ccn_number) + 
     &        pdqcloud(ig,l,igcm_ccn_number)) .le. 1.)
     &   .or. (pq(ig,l,igcm_ccn_mass) + 
     &      ptimestep* (pdq(ig,l,igcm_ccn_mass) + 
     &        pdqcloud(ig,l,igcm_ccn_mass)) .le. 1.e-20)) THEN
         pdqcloud(ig,l,igcm_ccn_number) =
     &     - pq(ig,l,igcm_ccn_number)/ptimestep 
     &     - pdq(ig,l,igcm_ccn_number) + 1.
         pdqcloud(ig,l,igcm_dust_number) =  
     &     -pdqcloud(ig,l,igcm_ccn_number)
         pdqcloud(ig,l,igcm_ccn_mass) =
     &     - pq(ig,l,igcm_ccn_mass)/ptimestep
     &     - pdq(ig,l,igcm_ccn_mass) + 1.e-20
         pdqcloud(ig,l,igcm_dust_mass) = 
     &     -pdqcloud(ig,l,igcm_ccn_mass)
          ENDIF
         ENDDO
        ENDDO
       ENDIF

       IF(scavenging) THEN
        DO l=1,nlay
         DO ig=1,ngrid
          IF ((pq(ig,l,igcm_dust_number) + 
     &      ptimestep* (pdq(ig,l,igcm_dust_number) + 
     &        pdqcloud(ig,l,igcm_dust_number)) .le. 1.)
     &   .or. (pq(ig,l,igcm_dust_mass) + 
     &      ptimestep* (pdq(ig,l,igcm_dust_mass) + 
     &        pdqcloud(ig,l,igcm_dust_mass)) .le. 1.e-20)) THEN
         pdqcloud(ig,l,igcm_dust_number) =
     &     - pq(ig,l,igcm_dust_number)/ptimestep 
     &     - pdq(ig,l,igcm_dust_number) + 1.
         pdqcloud(ig,l,igcm_ccn_number) =  
     &     -pdqcloud(ig,l,igcm_dust_number)
         pdqcloud(ig,l,igcm_dust_mass) =
     &     - pq(ig,l,igcm_dust_mass)/ptimestep
     &     - pdq(ig,l,igcm_dust_mass) + 1.e-20
         pdqcloud(ig,l,igcm_ccn_mass) = 
     &     -pdqcloud(ig,l,igcm_dust_mass)
          ENDIF
         ENDDO
        ENDDO
       ENDIF

        DO l=1,nlay
         DO ig=1,ngrid
          IF (pq(ig,l,igcm_h2o_ice) + ptimestep*
     &       (pdq(ig,l,igcm_h2o_ice) + pdqcloud(ig,l,igcm_h2o_ice)) 
     &       .le. 1.e-8) THEN
           pdqcloud(ig,l,igcm_h2o_ice) = 
     &     - pq(ig,l,igcm_h2o_ice)/ptimestep - pdq(ig,l,igcm_h2o_ice)
           pdqcloud(ig,l,igcm_h2o_vap) = -pdqcloud(ig,l,igcm_h2o_ice)
          ENDIF 
          IF (pq(ig,l,igcm_h2o_vap) + ptimestep*
     &       (pdq(ig,l,igcm_h2o_vap) + pdqcloud(ig,l,igcm_h2o_vap)) 
     &       .le. 1.e-8) THEN
           pdqcloud(ig,l,igcm_h2o_vap) = 
     &     - pq(ig,l,igcm_h2o_vap)/ptimestep - pdq(ig,l,igcm_h2o_vap)
           pdqcloud(ig,l,igcm_h2o_ice) = -pdqcloud(ig,l,igcm_h2o_vap)
          ENDIF
         ENDDO
        ENDDO


c------Update the ice and dust particle size "rice" for output or photochemistry
c------Only rsedcloud is used for the water cycle

      IF(scavenging) THEN 
        DO l=1, nlay
         DO ig=1,ngrid

        call updaterdust(
     &    pq(ig,l,igcm_dust_mass) +                   ! dust mass
     &   (pdq(ig,l,igcm_dust_mass) +                  ! dust mass
     &    pdqcloud(ig,l,igcm_dust_mass))*ptimestep,   ! dust mass
     &    pq(ig,l,igcm_dust_number) +                 ! dust number
     &   (pdq(ig,l,igcm_dust_number) +                ! dust number
     &    pdqcloud(ig,l,igcm_dust_number))*ptimestep, ! dust number
     &    rdust(ig,l))

         ENDDO
        ENDDO
      ENDIF
        
        
      IF(microphys) THEN
       
       DO l=1, nlay
         DO ig=1,ngrid

        call updaterice_micro(
     &    pq(ig,l,igcm_h2o_ice) +                    ! ice mass
     &   (pdq(ig,l,igcm_h2o_ice) +                   ! ice mass
     &    pdqcloud(ig,l,igcm_h2o_ice))*ptimestep,    ! ice mass
     &    pq(ig,l,igcm_ccn_mass) +                   ! ccn mass
     &   (pdq(ig,l,igcm_ccn_mass) +                  ! ccn mass
     &    pdqcloud(ig,l,igcm_ccn_mass))*ptimestep,   ! ccn mass
     &    pq(ig,l,igcm_ccn_number) +                 ! ccn number
     &   (pdq(ig,l,igcm_ccn_number) +                ! ccn number
     &    pdqcloud(ig,l,igcm_ccn_number))*ptimestep, ! ccn number
     &    tauscaling(ig),rice(ig,l),rhocloud(ig,l))
          
         ENDDO
       ENDDO
        
      ELSE ! no microphys
        
        DO l=1,nlay
          DO ig=1,ngrid
          
        call updaterice_typ(
     &    pq(ig,l,igcm_h2o_ice) +                    ! ice mass
     &   (pdq(ig,l,igcm_h2o_ice) +                   ! ice mass
     &    pdqcloud(ig,l,igcm_h2o_ice))*ptimestep,    ! ice mass
     &    tau(ig,1),pzlay(ig,l),rice(ig,l)) 

          ENDDO
         ENDDO
       
       ENDIF ! of IF(microphys)
      
      
      
c     A correction if a lot of subliming CO2 fills the 1st layer FF04/2005
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Then that should not affect the ice particle radius
      do ig=1,ngridmx
        if(pdpsrf(ig)*ptimestep.gt.0.9*(pplev(ig,1)-pplev(ig,2)))then
          if(pdpsrf(ig)*ptimestep.gt.0.9*(pplev(ig,1)-pplev(ig,3)))
     &    rice(ig,2)=rice(ig,3) 
          rice(ig,1)=rice(ig,2)
        end if
      end do
       
       
       DO l=1,nlay
         DO ig=1,ngrid
           rsedcloud(ig,l)=max(rice(ig,l)*
     &                 (1.+nuice_sed)*(1.+nuice_sed)*(1.+nuice_sed),
     &                    rdust(ig,l))
!          rsedcloud(ig,l)=min(rsedcloud(ig,l),1.e-4)
         ENDDO
       ENDDO
       
! used for rad. transfer calculations
! nuice is constant because a lognormal distribution is prescribed
      nuice(1:ngrid,1:nlay)=nuice_ref 



c=======================================================================

      END
 
