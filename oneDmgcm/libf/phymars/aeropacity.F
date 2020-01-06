      SUBROUTINE aeropacity(ngrid,nlayer,nq,zday,pplay,pplev,ls,
     &    pq,tauscaling,tauref,tau,taucloudtes,aerosol,reffrad,nueffrad,
     &    QREFvis3d,QREFir3d,omegaREFvis3d,omegaREFir3d)
                                                   
! to use  'getin'
      USE ioipsl_getincom 
       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Computing aerosol optical depth in each gridbox.
c
c   author: F.Forget 
c   ------
c   update F. Montmessin (water ice scheme) 
c      and S. Lebonnois (12/06/2003) compatibility dust/ice/chemistry
c   update J.-B. Madeleine 2008-2009:
c       - added 3D scattering by aerosols;
c       - dustopacity transferred from physiq.F to callradite.F,
c           and renamed into aeropacity.F;
c   update E. Millour, march 2012:
c         - reference pressure is now set to 610Pa (not 700Pa)
c   
c   input:
c   ----- 
c   ngrid             Number of gridpoint of horizontal grid
c   nlayer            Number of layer
c   nq                Number of tracer
c   zday                  Date (time since Ls=0, in martian days)
c   ls                Solar longitude (Ls) , radian
c   pplay,pplev       pressure (Pa) in the middle and boundary of each layer
c   pq                Dust mixing ratio (used if tracer =T and active=T).
c   reffrad(ngrid,nlayer,naerkind)  Aerosol effective radius
c   QREFvis3d(ngridmx,nlayermx,naerkind) \ 3d extinction coefficients
c   QREFir3d(ngridmx,nlayermx,naerkind)  / at reference wavelengths;
c   omegaREFvis3d(ngridmx,nlayermx,naerkind) \ 3d single scat. albedo
c   omegaREFir3d(ngridmx,nlayermx,naerkind)  / at reference wavelengths;
c
c   output:
c   -------
c   tauref       Prescribed mean column optical depth at 610 Pa 
c   tau          Column total visible dust optical depth at each point
c   aerosol      aerosol(ig,l,1) is the dust optical
c                depth in layer l, grid point ig

c
c=======================================================================
#include "dimensions.h"
#include "dimphys.h"
#include "callkeys.h"
#include "comcstfi.h"
#include "comgeomfi.h"
#include "dimradmars.h"
#include "yomaer.h"
#include "tracer.h"
#include "planete.h"
#include "aerkind.h"

c-----------------------------------------------------------------------
c
c    Declarations :
c    --------------
c
c    Input/Output
c    ------------
      INTEGER ngrid,nlayer,nq

      REAL ls,zday,expfactor    
      REAL pplev(ngrid,nlayer+1),pplay(ngrid,nlayer)
      REAL pq(ngrid,nlayer,nq)
      REAL tauref(ngrid), tau(ngrid,naerkind)
      REAL aerosol(ngrid,nlayer,naerkind)
      REAL dsodust(ngridmx,nlayermx)
      REAL reffrad(ngrid,nlayer,naerkind)
      REAL nueffrad(ngrid,nlayer,naerkind)
      REAL QREFvis3d(ngridmx,nlayermx,naerkind)
      REAL QREFir3d(ngridmx,nlayermx,naerkind)
      REAL omegaREFvis3d(ngridmx,nlayermx,naerkind)
      REAL omegaREFir3d(ngridmx,nlayermx,naerkind)
c
c    Local variables :
c    -----------------
      INTEGER l,ig,iq,i,j
      INTEGER iaer           ! Aerosol index
      real topdust(ngridmx)
      real zlsconst, zp
      real taueq,tauS,tauN
c     Mean Qext(vis)/Qext(ir) profile
      real msolsir(nlayermx,naerkind)
c     Mean Qext(ir)/Qabs(ir) profile
      real mqextsqabs(nlayermx,naerkind)
c     Variables used when multiple particle sizes are used
c       for dust or water ice particles in the radiative transfer
c       (see callradite.F for more information).
      REAL taudusttmp(ngridmx)! Temporary dust opacity
                               !   used before scaling
      REAL tauscaling(ngridmx) ! Scaling factor for qdust and Ndust
      REAL taudustvis(ngridmx) ! Dust opacity after scaling
      REAL taudusttes(ngridmx) ! Dust opacity at IR ref. wav. as
                               !   "seen" by the GCM.
      REAL taucloudvis(ngridmx)! Cloud opacity at visible
                               !   reference wavelength
      REAL taucloudtes(ngridmx)! Cloud opacity at infrared
                               !   reference wavelength using
                               !   Qabs instead of Qext
                               !   (direct comparison with TES)

c   local saved variables
c   ---------------------

      REAL topdust0(ngridmx) 
      SAVE topdust0
c     Level under which the dust mixing ratio is held constant
c       when computing the dust opacity in each layer
c       (this applies when doubleq and active are true)
      INTEGER, PARAMETER :: cstdustlevel = 7

      LOGICAL,SAVE :: firstcall=.true.

! indexes of water ice and dust tracers:
      INTEGER,SAVE :: nqdust(nqmx) ! to store the indexes of dust tracers
      INTEGER,SAVE :: i_ice=0  ! water ice
      real,parameter :: odpref=610. ! DOD reference pressure (Pa)
      CHARACTER(LEN=20) :: txt ! to temporarly store text
      CHARACTER(LEN=1) :: txt2 ! to temporarly store text
! indexes of dust scatterers:
      INTEGER,SAVE :: iaerdust(naerkind)
      INTEGER,SAVE :: naerdust ! number of dust scatterers

      tau(1:ngrid,1:naerkind)=0

! identify tracers

      IF (firstcall) THEN
        ! identify scatterers that are dust
        naerdust=0
        DO iaer=1,naerkind
          txt=name_iaer(iaer)
          IF (txt(1:4).eq."dust") THEN
            naerdust=naerdust+1
            iaerdust(naerdust)=iaer
          ENDIF
        ENDDO
        ! identify tracers which are dust
        i=0
        DO iq=1,nq
          txt=noms(iq)
          IF (txt(1:4).eq."dust") THEN
          i=i+1
          nqdust(i)=iq
          ENDIF
        ENDDO

        IF (water.AND.activice) THEN
          i_ice=igcm_h2o_ice
          write(*,*) "aeropacity: i_ice=",i_ice
        ENDIF

c       altitude of the top of the aerosol layer (km) at Ls=2.76rad:
c       in the Viking year scenario
        DO ig=1,ngrid
            topdust0(ig)=60. -22.*SIN(lati(ig))**2
        END DO

c       typical profile of solsir and (1-w)^(-1):
        msolsir(1:nlayer,1:naerkind)=0
        mqextsqabs(1:nlayer,1:naerkind)=0
        WRITE(*,*) "Typical profiles of Qext(vis)/Qext(IR)"
        WRITE(*,*) "  and Qext(IR)/Qabs(IR):"
        DO iaer = 1, naerkind ! Loop on aerosol kind
          WRITE(*,*) "Aerosol # ",iaer
          DO l=1,nlayer
            DO ig=1,ngridmx
              msolsir(l,iaer)=msolsir(l,iaer)+
     &              QREFvis3d(ig,l,iaer)/
     &              QREFir3d(ig,l,iaer)
              mqextsqabs(l,iaer)=mqextsqabs(l,iaer)+
     &              (1.E0-omegaREFir3d(ig,l,iaer))**(-1)
            ENDDO
            msolsir(l,iaer)=msolsir(l,iaer)/REAL(ngridmx)
            mqextsqabs(l,iaer)=mqextsqabs(l,iaer)/REAL(ngridmx)
          ENDDO
          WRITE(*,*) "solsir: ",msolsir(:,iaer)
          WRITE(*,*) "Qext/Qabs(IR): ",mqextsqabs(:,iaer)
        ENDDO

!       load value of tauvis from callphys.def (if given there,
!       otherwise default value read from starfi.nc file will be used)
        call getin("tauvis",tauvis)

        firstcall=.false.

      END IF

c     Vertical column optical depth at "odpref" Pa 
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(iaervar.eq.1) THEN 
         do ig=1, ngridmx
          tauref(ig)=max(tauvis,1.e-9) ! tauvis=cste (set in callphys.def
                                       ! or read in starfi
        end do
      ELSE IF (iaervar.eq.2) THEN   ! << "Viking" Scenario>>

        tauref(1) = 0.7+.3*cos(ls+80.*pi/180.) ! like seen by VL1
        do ig=2,ngrid
          tauref(ig) = tauref(1)
        end do

      ELSE IF (iaervar.eq.3) THEN  ! << "MGS" scenario >>

        taueq= 0.2 +(0.5-0.2) *(cos(0.5*(ls-4.363)))**14
        tauS= 0.1 +(0.5-0.1)  *(cos(0.5*(ls-4.363)))**14
        tauN = 0.1
c	   if (peri_day.eq.150) then
c	     tauS=0.1
c	     tauN=0.1 +(0.5-0.1)  *(cos(0.5*(ls+pi-4.363)))**14
c	     taueq= 0.2 +(0.5-0.2) *(cos(0.5*(ls+pi-4.363)))**14
c           endif
        do ig=1,ngrid/2  ! Northern hemisphere
          tauref(ig)= tauN +
     &    (taueq-tauN)*0.5*(1+tanh((45-lati(ig)*180./pi)*6/60))
        end do
        do ig=ngrid/2+1, ngridmx  ! Southern hemisphere
          tauref(ig)= tauS +
     &    (taueq-tauS)*0.5*(1+tanh((45+lati(ig)*180./pi)*6/60))
        end do
      ELSE IF (iaervar.eq.5) THEN   ! << Escalier Scenario>>
c         tauref(1) = 0.2
c         if ((ls.ge.210.*pi/180.).and.(ls.le.330.*pi/180.))
c    &                              tauref(1) = 2.5
        tauref(1) = 2.5
        if ((ls.ge.30.*pi/180.).and.(ls.le.150.*pi/180.))
     &                              tauref(1) = .2

        do ig=2,ngrid
          tauref(ig) = tauref(1)
        end do
      ELSE IF ((iaervar.ge.6).and.(iaervar.le.7)) THEN
      ! cold or warm synthetic scenarios
        call read_dust_scenario(ngrid,nlayer,zday,pplev,tauref)
      ELSE IF ((iaervar.ge.24).and.(iaervar.le.30))
     &     THEN  ! << MY... dust scenarios >>
        call read_dust_scenario(ngrid,nlayer,zday,pplev,tauref)
      ELSE IF ((iaervar.eq.4).or.
     &         ((iaervar.ge.124).and.(iaervar.le.126))) THEN
       ! "old" TES assimation dust scenario (values at 700Pa in files!)
        call read_dust_scenario(ngrid,nlayer,zday,pplev,tauref)
      ELSE
        stop 'problem with iaervar in aeropacity.F'
      ENDIF

c -----------------------------------------------------------------
c Computing the opacity in each layer
c -----------------------------------------------------------------

      DO iaer = 1, naerkind ! Loop on aerosol kind
c     --------------------------------------------
        aerkind: SELECT CASE (name_iaer(iaer))
c==================================================================
        CASE("dust_conrath") aerkind      ! Typical dust profile
c==================================================================

c       Altitude of the top of the dust layer
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        zlsconst=SIN(ls-2.76)
        if (iddist.eq.1) then
          do ig=1,ngrid
             topdust(ig)=topdustref         ! constant dust layer top
          end do

        else if (iddist.eq.2) then          ! "Viking" scenario
          do ig=1,ngrid
            topdust(ig)=topdust0(ig)+18.*zlsconst
          end do

        else if(iddist.eq.3) then         !"MGS" scenario
          do ig=1,ngrid
            topdust(ig)=60.+18.*zlsconst
     &                -(32+18*zlsconst)*sin(lati(ig))**4
     &                 - 8*zlsconst*(sin(lati(ig)))**5
          end do
        endif

c       Optical depth in each layer :
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(iddist.ge.1) then

          expfactor=0.
          DO l=1,nlayer
            DO ig=1,ngrid
c             Typical mixing ratio profile 
              if(pplay(ig,l).gt.odpref
     $                        /(988.**(topdust(ig)/70.))) then
                zp=(odpref/pplay(ig,l))**(70./topdust(ig))
                 expfactor=max(exp(0.007*(1.-max(zp,1.))),1.e-3)
              else    
                expfactor=1.e-3
              endif
c             Vertical scaling function
              aerosol(ig,l,iaer)= (pplev(ig,l)-pplev(ig,l+1)) *
     &          expfactor *
     &          QREFvis3d(ig,l,iaer) / QREFvis3d(ig,1,iaer)
            ENDDO
          ENDDO

        else if(iddist.eq.0) then   
c         old dust vertical distribution function (pollack90)
          DO l=1,nlayer
             DO ig=1,ngrid
                zp=odpref/pplay(ig,l)
                aerosol(ig,l,1)= tauref(ig)/odpref *
     s           (pplev(ig,l)-pplev(ig,l+1))
     s           *max( exp(.03*(1.-max(zp,1.))) , 1.E-3 )
             ENDDO
          ENDDO
        end if

c==================================================================
        CASE("dust_doubleq") aerkind! Two-moment scheme for dust
c        (transport of mass and number mixing ratio)
c==================================================================

          DO l=1,nlayer
            IF (l.LE.cstdustlevel) THEN
c           Opacity in the first levels is held constant to 
c             avoid unrealistic values due to constant lifting:
              DO ig=1,ngrid
                aerosol(ig,l,iaer) = 
     &          (  0.75 * QREFvis3d(ig,cstdustlevel,iaer) /
     &          ( rho_dust * reffrad(ig,cstdustlevel,iaer) )  ) *
     &          pq(ig,cstdustlevel,igcm_dust_mass) *
     &          ( pplev(ig,l) - pplev(ig,l+1) ) / g
              ENDDO
            ELSE
              DO ig=1,ngrid
                aerosol(ig,l,iaer) =
     &          (  0.75 * QREFvis3d(ig,l,iaer) /
     &          ( rho_dust * reffrad(ig,l,iaer) )  ) *
     &          pq(ig,l,igcm_dust_mass) *
     &          ( pplev(ig,l) - pplev(ig,l+1) ) / g
              ENDDO
            ENDIF
          ENDDO

c==================================================================
        CASE("dust_submicron") aerkind   ! Small dust population
c==================================================================

          DO l=1,nlayer
            IF (l.LE.cstdustlevel) THEN
c           Opacity in the first levels is held constant to 
c             avoid unrealistic values due to constant lifting:
              DO ig=1,ngrid
                aerosol(ig,l,iaer) = 
     &          (  0.75 * QREFvis3d(ig,cstdustlevel,iaer) /
     &          ( rho_dust * reffrad(ig,cstdustlevel,iaer) )  ) *
     &          pq(ig,cstdustlevel,igcm_dust_submicron) *
     &          ( pplev(ig,l) - pplev(ig,l+1) ) / g
              ENDDO
            ELSE
              DO ig=1,ngrid
                aerosol(ig,l,iaer) = 
     &          (  0.75 * QREFvis3d(ig,l,iaer) /
     &          ( rho_dust * reffrad(ig,l,iaer) )  ) *
     &          pq(ig,l,igcm_dust_submicron) *
     &          ( pplev(ig,l) - pplev(ig,l+1) ) / g
              ENDDO
            ENDIF
          ENDDO

c==================================================================
        CASE("h2o_ice") aerkind             ! Water ice crystals
c==================================================================

c       1. Initialization
        aerosol(1:ngrid,1:nlayer,iaer) = 0.
        taucloudvis(1:ngrid) = 0.
        taucloudtes(1:ngrid) = 0.
c       2. Opacity calculation
        DO ig=1, ngrid
          DO l=1,nlayer
            aerosol(ig,l,iaer) = max(1E-20,
     &        (  0.75 * QREFvis3d(ig,l,iaer) /
     &        ( rho_ice * reffrad(ig,l,iaer) )  ) *
     &        pq(ig,l,i_ice) *
     &        ( pplev(ig,l) - pplev(ig,l+1) ) / g
     &                              )
            taucloudvis(ig) = taucloudvis(ig) + aerosol(ig,l,iaer)
            taucloudtes(ig) = taucloudtes(ig) + aerosol(ig,l,iaer)*
     &        QREFir3d(ig,l,iaer) / QREFvis3d(ig,l,iaer) *
     &        ( 1.E0 - omegaREFir3d(ig,l,iaer) )
          ENDDO
        ENDDO
c       3. Outputs -- Now done in physiq.F
!        IF (ngrid.NE.1) THEN
!          CALL WRITEDIAGFI(ngridmx,'tauVIS','tauext VIS refwvl',
!     &      ' ',2,taucloudvis)
!          CALL WRITEDIAGFI(ngridmx,'tauTES','tauabs IR refwvl',
!     &      ' ',2,taucloudtes)
!          IF (callstats) THEN
!            CALL wstats(ngridmx,'tauVIS','tauext VIS refwvl',
!     &        ' ',2,taucloudvis)
!            CALL wstats(ngridmx,'tauTES','tauabs IR refwvl',
!     &        ' ',2,taucloudtes)
!          ENDIF
!        ELSE
!         CALL writeg1d(ngrid,1,taucloudtes,'tautes','NU')
!        ENDIF
c==================================================================
        END SELECT aerkind
c     -----------------------------------
      ENDDO ! iaer (loop on aerosol kind)

c -----------------------------------------------------------------
c Rescaling each layer to reproduce the choosen (or assimilated)
c   dust extinction opacity at visible reference wavelength, which
c   is originally scaled to an equivalent odpref Pa pressure surface.
c -----------------------------------------------------------------

c     Temporary scaling factor
      taudusttmp(1:ngrid)=0.
      DO iaer=1,naerdust
        DO l=1,nlayer
          DO ig=1,ngrid
c           Scaling factor
            taudusttmp(ig) = taudusttmp(ig) + 
     &                       aerosol(ig,l,iaerdust(iaer))
          ENDDO
        ENDDO
      ENDDO

c     Saved scaling factor
      DO ig=1,ngrid
          tauscaling(ig) = tauref(ig) * 
     &                     pplev(ig,1) / odpref / taudusttmp(ig)
c          tauscaling(ig) = 1.e-4
      ENDDO

c     Opacity computation
      DO iaer=1,naerdust
        DO l=1,nlayer
          DO ig=1,ngrid
            aerosol(ig,l,iaerdust(iaer)) = max(1E-20,
     &                aerosol(ig,l,iaerdust(iaer))* tauscaling(ig))
          ENDDO
        ENDDO
      ENDDO
      
c output for debug
c        IF (ngrid.NE.1) THEN
c             CALL WRITEDIAGFI(ngridmx,'taudusttmp','virtual tau dust',
c     &      '#',2,taudusttmp)
c             CALL WRITEDIAGFI(ngridmx,'tausca','tauscaling',
c     &      '#',2,tauscaling)
c        ELSE
c             CALL WRITEDIAGFI(ngridmx,'taudusttmp','virtual tau dust',
c     &      '#',0,taudusttmp)
c             CALL WRITEDIAGFI(ngridmx,'tausca','tauscaling',
c     &      '#',0,tauscaling)
c        ENDIF
c -----------------------------------------------------------------
c Column integrated visible optical depth in each point
c -----------------------------------------------------------------
      DO iaer=1,naerkind
        do l=1,nlayer
           do ig=1,ngrid
             tau(ig,iaer) = tau(ig,iaer) + aerosol(ig,l,iaer)
           end do
        end do
      ENDDO
c -----------------------------------------------------------------
c Density scaled opacity and column opacity output
c -----------------------------------------------------------------
c      dsodust(1:ngrid,1:nlayer) = 0.
c      DO iaer=1,naerdust
c        DO l=1,nlayermx
c          DO ig=1,ngrid
c            dsodust(ig,l) = dsodust(ig,l) +
c     &                      aerosol(ig,l,iaerdust(iaer)) * g /
c     &                      (pplev(ig,l) - pplev(ig,l+1))
c          ENDDO
c        ENDDO
c        IF (ngrid.NE.1) THEN
c          write(txt2,'(i1.1)') iaer
c          call WRITEDIAGFI(ngridmx,'taudust'//txt2,
c     &                    'Dust col opacity',
c     &                    ' ',2,tau(1,iaerdust(iaer)))
c          IF (callstats) THEN
c            CALL wstats(ngridmx,'taudust'//txt2,
c     &                 'Dust col opacity',
c     &                 ' ',2,tau(1,iaerdust(iaer)))
c          ENDIF
c        ENDIF
c      ENDDO

c      IF (ngrid.NE.1) THEN
c       CALL WRITEDIAGFI(ngridmx,'dsodust','tau*g/dp',
c    &                    'm2.kg-1',3,dsodust)
c        IF (callstats) THEN
c          CALL wstats(ngridmx,'dsodust',
c     &               'tau*g/dp',
c     &               'm2.kg-1',3,dsodust)
c        ENDIF
c      ELSE
c        CALL WRITEDIAGFI(ngrid,"dsodust","dsodust","m2.kg-1",1,
c     &                   dsodust)
c      ENDIF ! of IF (ngrid.NE.1)
c -----------------------------------------------------------------
      return
      end 
