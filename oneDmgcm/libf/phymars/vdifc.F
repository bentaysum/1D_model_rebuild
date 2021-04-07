      SUBROUTINE vdifc(ngrid,nlay,nq,co2ice,ppopsk,
     $                ptimestep,pcapcal,lecrit,
     $                pplay,pplev,pzlay,pzlev,pz0,
     $                pu,pv,ph,pq,ptsrf,pemis,pqsurf,
     $                pdufi,pdvfi,pdhfi,pdqfi,pfluxsrf,
     $                pdudif,pdvdif,pdhdif,pdtsrf,pq2,
     $                pdqdif,pdqsdif,wstar,zcdv_true,zcdh_true,
     $                hfmax,sensibFlux
#ifdef MESOSCALE
     &                ,flag_LES
#endif
     &                )

      USE TLMvars

      IMPLICIT NONE

c=======================================================================
c
c   subject:
c   --------
c   Turbulent diffusion (mixing) for potential T, U, V and tracer
c
c   Shema implicite
c   On commence par rajouter au variables x la tendance physique
c   et on resoult en fait:
c      x(t+1) =  x(t) + dt * (dx/dt)phys(t)  +  dt * (dx/dt)difv(t+1)
c
c   author:
c   ------
c      Hourdin/Forget/Fournier
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "surfdat.h"
#include "comgeomfi.h"
#include "tracer.h"
#include "microphys.h"

c
c   arguments:
c   ----------

      INTEGER ngrid,nlay
      REAL ptimestep
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      REAL pzlay(ngrid,nlay),pzlev(ngrid,nlay+1)
      REAL pu(ngrid,nlay),pv(ngrid,nlay),ph(ngrid,nlay),pt(ngrid,nlay)
      REAL ptsrf(ngrid),pemis(ngrid)
      REAL pdufi(ngrid,nlay),pdvfi(ngrid,nlay),pdhfi(ngrid,nlay)
      REAL pfluxsrf(ngrid)
      REAL pdudif(ngrid,nlay),pdvdif(ngrid,nlay),pdhdif(ngrid,nlay)
      REAL pdtsrf(ngrid),pcapcal(ngrid)
      REAL pq2(ngrid,nlay+1)

c    Argument added for condensation:
      REAL co2ice (ngrid), ppopsk(ngrid,nlay)
      logical lecrit

      REAL pz0(ngridmx) ! surface roughness length (m)

c    Argument added to account for subgrid gustiness :

      REAL wstar(ngridmx), hfmax(ngridmx), zi(ngridmx)

c    Traceurs :
      integer nq 
      REAL pqsurf(ngrid,nq)
      real pq(ngrid,nlay,nq), pdqfi(ngrid,nlay,nq) 
      real pdqdif(ngrid,nlay,nq) 
      real pdqsdif(ngrid,nq) 
      
c   local:
c   ------

      REAL ust(ngridmx),tst(ngridmx)
 
      INTEGER ilev,ig,ilay,nlev

      REAL z4st,zdplanck(ngridmx)
      REAL zkv(ngridmx,nlayermx+1),zkh(ngridmx,nlayermx+1)
      REAL zkq(ngridmx,nlayermx+1)
      REAL zcdv(ngridmx),zcdh(ngridmx)
      REAL zcdv_true(ngridmx),zcdh_true(ngridmx)    ! drag coeff are used by the LES to recompute u* and hfx
      REAL zu(ngridmx,nlayermx),zv(ngridmx,nlayermx)
      REAL zh(ngridmx,nlayermx)
      REAL ztsrf2(ngridmx)
      REAL z1(ngridmx),z2(ngridmx)
      REAL za(ngridmx,nlayermx),zb(ngridmx,nlayermx)
      REAL zb0(ngridmx,nlayermx)
      REAL zc(ngridmx,nlayermx),zd(ngridmx,nlayermx)
      REAL zcst1
      REAL zu2(ngridmx)

      EXTERNAL SSUM,SCOPY
      REAL SSUM
      LOGICAL firstcall
      SAVE firstcall


c     variable added for CO2 condensation:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      REAL hh , zhcond(ngridmx,nlayermx)
      REAL latcond,tcond1mb
      REAL acond,bcond
      SAVE acond,bcond
      DATA latcond,tcond1mb/5.9e5,136.27/
      
c     For latent heat release from ground ice sublimation    
      REAL tsrf_lw(ngridmx)
      REAL alpha
      REAL T1,T2
      SAVE T1,T2
      DATA T1,T2/-604.3,1080.7/ ! zeros of latent heat equation for ice

c    Tracers :
c    ~~~~~~~ 
      INTEGER iq
      REAL zq(ngridmx,nlayermx,nqmx)
      REAL zq1temp(ngridmx)
      REAL rho(ngridmx) ! near surface air density
      REAL qsat(ngridmx)
      DATA firstcall/.true./

      REAL kmixmin

c    Mass-variation scheme :
c    ~~~~~~~

      INTEGER j,l
      REAL zcondicea(ngrid,nlayermx)
      REAL zt(ngridmx,nlayermx),ztcond(ngridmx,nlayermx+1)
      REAL betam(ngridmx,nlayermx),dmice(ngridmx,nlayermx)
      REAL pdtc(ngrid,nlayermx)
      REAL zhs(ngridmx,nlayermx)
      REAL ccond
      SAVE ccond

c     Theta_m formulation for mass-variation scheme :
c     ~~~~~~~

      INTEGER ico2,llnt(ngridmx)
      SAVE ico2
      REAL m_co2, m_noco2, A , B
      SAVE A, B, m_co2, m_noco2
      REAL vmr_co2(ngridmx,nlayermx)
      REAL qco2,mmean

      REAL sensibFlux(ngridmx)

#ifdef MESOSCALE
      LOGICAL flag_LES     !! pour LES avec isfflx!=0
#endif
c    ** un petit test de coherence
c       --------------------------

      IF (firstcall) THEN
         IF(ngrid.NE.ngridmx) THEN
            PRINT*,'STOP dans vdifc'
            PRINT*,'probleme de dimensions :'
            PRINT*,'ngrid  =',ngrid
            PRINT*,'ngridmx  =',ngridmx
            STOP
         ENDIF
c        To compute: Tcond= 1./(bcond-acond*log(.0095*p)) (p in pascal)
         bcond=1./tcond1mb
         acond=r/latcond
         ccond=cpp/(g*latcond)
         PRINT*,'In vdifc: Tcond(P=1mb)=',tcond1mb,' Lcond=',latcond
         PRINT*,'          acond,bcond,ccond',acond,bcond,ccond


         ico2=0

         if (tracer) then
c          Prepare Special treatment if one of the tracer is CO2 gas
           do iq=1,nqmx
             if (noms(iq).eq."co2") then
                ico2=iq
                m_co2 = 44.01E-3  ! CO2 molecular mass (kg/mol)
                m_noco2 = 33.37E-3  ! Non condensible mol mass (kg/mol)
c               Compute A and B coefficient use to compute
c               mean molecular mass Mair defined by
c               1/Mair = q(ico2)/m_co2 + (1-q(ico2))/m_noco2
c               1/Mair = A*q(ico2) + B
                A =(1/m_co2 - 1/m_noco2)
                B=1/m_noco2
             endif
           enddo
         end if

        firstcall=.false.
      ENDIF





c-----------------------------------------------------------------------
c    1. initialisation
c    -----------------

      nlev=nlay+1

c    ** calcul de rho*dz et dt*rho/dz=dt*rho**2 g/dp
c       avec rho=p/RT=p/ (R Theta) (p/ps)**kappa
c       ----------------------------------------

      DO ilay=1,nlay
         DO ig=1,ngrid
            za(ig,ilay)=(pplev(ig,ilay)-pplev(ig,ilay+1))/g
! Mass variation scheme:
            betam(ig,ilay)=-za(ig,ilay)*latcond/(cpp*ppopsk(ig,ilay))
         ENDDO
      ENDDO

      zcst1=4.*g*ptimestep/(r*r)
      DO ilev=2,nlev-1
         DO ig=1,ngrid
            zb0(ig,ilev)=pplev(ig,ilev)*
     s      (pplev(ig,1)/pplev(ig,ilev))**rcp /
     s      (ph(ig,ilev-1)+ph(ig,ilev))
            zb0(ig,ilev)=zcst1*zb0(ig,ilev)*zb0(ig,ilev)/
     s      (pplay(ig,ilev-1)-pplay(ig,ilev))
         ENDDO
      ENDDO
      DO ig=1,ngrid
         zb0(ig,1)=ptimestep*pplev(ig,1)/(r*ptsrf(ig))
      ENDDO

c    ** diagnostique pour l'initialisation
c       ----------------------------------

      IF(lecrit) THEN
         ig=ngrid/2+1
         PRINT*,'Pression (mbar) ,altitude (km),u,v,theta, rho dz'
         DO ilay=1,nlay
            WRITE(*,'(6f11.5)')
     s      .01*pplay(ig,ilay),.001*pzlay(ig,ilay),
     s      pu(ig,ilay),pv(ig,ilay),ph(ig,ilay),za(ig,ilay)
         ENDDO
         PRINT*,'Pression (mbar) ,altitude (km),zb'
         DO ilev=1,nlay
            WRITE(*,'(3f15.7)')
     s      .01*pplev(ig,ilev),.001*pzlev(ig,ilev),
     s      zb0(ig,ilev)
         ENDDO
      ENDIF

c     -----------------------------------
c     Potential Condensation temperature:
c     -----------------------------------

c     Compute CO2 Volume mixing ratio
c     -------------------------------
      if (callcond.and.(ico2.ne.0)) then
         DO ilev=1,nlay
            DO ig=1,ngrid
              qco2=MAX(1.E-30
     &             ,pq(ig,ilev,ico2)+pdqfi(ig,ilev,ico2)*ptimestep)
c             Mean air molecular mass = 1/(q(ico2)/m_co2 + (1-q(ico2))/m_noco2)
              mmean=1/(A*qco2 +B)
              vmr_co2(ig,ilev) = qco2*mmean/m_co2
            ENDDO
         ENDDO
      else
         DO ilev=1,nlay
            DO ig=1,ngrid
              vmr_co2(ig,ilev)=0.95
            ENDDO
         ENDDO
      end if

c     forecast of atmospheric temperature zt and frost temperature ztcond
c     --------------------------------------------------------------------

      DO ilev=1,nlay
         DO ig=1,ngrid
            ztcond(ig,ilev)=
     &         1./(bcond-acond*log(.01*vmr_co2(ig,ilev)*pplay(ig,ilev)))
            if (pplay(ig,ilev).lt.1e-4) ztcond(ig,ilev)=0.0 !mars Monica
         ENDDO
      ENDDO

      ztcond(:,nlay+1)=ztcond(:,nlay)

      if (callcond) then 
        DO ilev=1,nlay
          DO ig=1,ngrid
!            zhcond(ig,ilev) =
!     &  (1./(bcond-acond*log(.0095*pplay(ig,ilev))))/ppopsk(ig,ilev)
              zhcond(ig,ilev) = ztcond(ig,ilev)/ppopsk(ig,ilev)
          END DO
        END DO
      else
         call zerophys(ngrid*nlay,zhcond)
      end if


c-----------------------------------------------------------------------
c   2. ajout des tendances physiques
c      -----------------------------

      DO ilev=1,nlay
         DO ig=1,ngrid
            zu(ig,ilev)=pu(ig,ilev)+pdufi(ig,ilev)*ptimestep
            zv(ig,ilev)=pv(ig,ilev)+pdvfi(ig,ilev)*ptimestep
            zh(ig,ilev)=ph(ig,ilev)+pdhfi(ig,ilev)*ptimestep
!            zh(ig,ilev)=max(zh(ig,ilev),zhcond(ig,ilev))
         ENDDO
      ENDDO
      if(tracer) then
        DO iq =1, nq
         DO ilev=1,nlay
           DO ig=1,ngrid
              zq(ig,ilev,iq)=pq(ig,ilev,iq)+pdqfi(ig,ilev,iq)*ptimestep
           ENDDO
         ENDDO
        ENDDO
      end if

c-----------------------------------------------------------------------
c   3. schema de turbulence
c      --------------------

c    ** source d'energie cinetique turbulente a la surface
c       (condition aux limites du schema de diffusion turbulente
c       dans la couche limite
c       ---------------------

      CALL vdif_cd(ngrid,nlay,pz0,g,pzlay,pu,pv,wstar,ptsrf,ph
     &             ,zcdv_true,zcdh_true
#ifdef MESOSCALE
     &                ,flag_LES
#endif
     &           )

        zu2(:)=pu(:,1)*pu(:,1)+pv(:,1)*pv(:,1)

        IF (callrichsl) THEN
          zcdv(:)=zcdv_true(:)*sqrt(zu2(:)+
     &     (log(1.+0.7*wstar(:) + 2.3*wstar(:)**2))**2)
          zcdh(:)=zcdh_true(:)*sqrt(zu2(:)+
     &     (log(1.+0.7*wstar(:) + 2.3*wstar(:)**2))**2)

           ust(:)=sqrt(zcdv_true(:))*sqrt(zu2(:)+
     &     (log(1.+0.7*wstar(:) + 2.3*wstar(:)**2))**2)

           tst(:)=0.
           DO ig=1,ngrid
              IF (zcdh_true(ig) .ne. 0.) THEN      ! When Cd=Ch=0, u*=t*=0
                 tst(ig)=(ph(ig,1)-ptsrf(ig))*zcdh(ig)/ust(ig)
              ENDIF
           ENDDO

        ELSE
           zcdv(:)=zcdv_true(:)*sqrt(zu2(:))     ! 1 / bulk aerodynamic momentum conductance 
           zcdh(:)=zcdh_true(:)*sqrt(zu2(:))     ! 1 / bulk aerodynamic heat conductance
           ust(:)=sqrt(zcdv_true(:))*sqrt(zu2(:))
           tst(:)=(ph(:,1)-ptsrf(:))*zcdh_true(:)/sqrt(zcdv_true(:))
        ENDIF

! Some usefull diagnostics for the new surface layer parametrization :

!        call WRITEDIAGFI(ngridmx,'pcdv',
!     &              'momentum drag','no units',
!     &                         2,zcdv_true)
!        call WRITEDIAGFI(ngridmx,'pcdh',
!     &              'heat drag','no units',
!     &                         2,zcdh_true)
!        call WRITEDIAGFI(ngridmx,'ust',
!     &              'friction velocity','m/s',2,ust)
!       call WRITEDIAGFI(ngridmx,'tst',
!     &              'friction temperature','K',2,tst)
!        call WRITEDIAGFI(ngridmx,'rm-1',
!     &              'aerodyn momentum conductance','m/s',
!     &                         2,zcdv)
!        call WRITEDIAGFI(ngridmx,'rh-1',
!     &              'aerodyn heat conductance','m/s',
!     &                         2,zcdh)

c    ** schema de diffusion turbulente dans la couche limite
c       ---------------------------------------------------- 
       IF (.not. calltherm) THEN

       CALL vdif_kc(ptimestep,g,pzlev,pzlay
     &              ,pu,pv,ph,zcdv_true
     &              ,pq2,zkv,zkh,zq)

      ELSE

      pt(:,:)=ph(:,:)*ppopsk(:,:)
      CALL yamada4(ngrid,nlay,ptimestep,g,r,pplev,pt
     s   ,pzlev,pzlay,pu,pv,ph,pq,zcdv_true,pq2,zkv,zkh,zkq,ust
     s   ,9)

      ENDIF

      if ((doubleq).and.(ngrid.eq.1)) then
        kmixmin = 80. !80.! minimum eddy mix coeff in 1D
        do ilev=1,nlay
          do ig=1,ngrid
           zkh(ig,ilev) = max(kmixmin,zkh(ig,ilev))
           zkv(ig,ilev) = max(kmixmin,zkv(ig,ilev))
          end do
        end do
      end if

c    ** diagnostique pour le schema de turbulence
c       -----------------------------------------

      IF(lecrit) THEN
         PRINT*
         PRINT*,'Diagnostic for the vertical turbulent mixing'
         PRINT*,'Cd for momentum and potential temperature'

         PRINT*,zcdv(ngrid/2+1),zcdh(ngrid/2+1)
         PRINT*,'Mixing coefficient for momentum and pot.temp.'
         DO ilev=1,nlay
            PRINT*,zkv(ngrid/2+1,ilev),zkh(ngrid/2+1,ilev)
         ENDDO
      ENDIF




c-----------------------------------------------------------------------
c   4. inversion pour l'implicite sur u
c      --------------------------------

c    ** l'equation est 
c       u(t+1) =  u(t) + dt * {(du/dt)phys}(t)  +  dt * {(du/dt)difv}(t+1)
c       avec
c       /zu/ = u(t) + dt * {(du/dt)phys}(t)   (voir paragraphe 2.)
c       et
c       dt * {(du/dt)difv}(t+1) = dt * {(d/dz)[ Ku (du/dz) ]}(t+1)
c       donc les entrees sont /zcdv/ pour la condition a la limite sol
c       et /zkv/ = Ku
 
      CALL multipl((nlay-1)*ngrid,zkv(1,2),zb0(1,2),zb(1,2))
      CALL multipl(ngrid,zcdv,zb0,zb)

      DO ig=1,ngrid
         z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
         zc(ig,nlay)=za(ig,nlay)*zu(ig,nlay)*z1(ig)
         zd(ig,nlay)=zb(ig,nlay)*z1(ig)
      ENDDO

      DO ilay=nlay-1,1,-1
         DO ig=1,ngrid
            z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $         zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
            zc(ig,ilay)=(za(ig,ilay)*zu(ig,ilay)+
     $         zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
            zd(ig,ilay)=zb(ig,ilay)*z1(ig)
         ENDDO
      ENDDO

      DO ig=1,ngrid
         zu(ig,1)=zc(ig,1)
      ENDDO
      DO ilay=2,nlay
         DO ig=1,ngrid
            zu(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zu(ig,ilay-1)
         ENDDO
      ENDDO





c-----------------------------------------------------------------------
c   5. inversion pour l'implicite sur v
c      --------------------------------

c    ** l'equation est 
c       v(t+1) =  v(t) + dt * {(dv/dt)phys}(t)  +  dt * {(dv/dt)difv}(t+1)
c       avec
c       /zv/ = v(t) + dt * {(dv/dt)phys}(t)   (voir paragraphe 2.)
c       et
c       dt * {(dv/dt)difv}(t+1) = dt * {(d/dz)[ Kv (dv/dz) ]}(t+1)
c       donc les entrees sont /zcdv/ pour la condition a la limite sol
c       et /zkv/ = Kv

      DO ig=1,ngrid
         z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
         zc(ig,nlay)=za(ig,nlay)*zv(ig,nlay)*z1(ig)
         zd(ig,nlay)=zb(ig,nlay)*z1(ig)
      ENDDO

      DO ilay=nlay-1,1,-1
         DO ig=1,ngrid
            z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $         zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
            zc(ig,ilay)=(za(ig,ilay)*zv(ig,ilay)+
     $         zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
            zd(ig,ilay)=zb(ig,ilay)*z1(ig)
         ENDDO
      ENDDO

      DO ig=1,ngrid
         zv(ig,1)=zc(ig,1)
      ENDDO
      DO ilay=2,nlay
         DO ig=1,ngrid
            zv(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zv(ig,ilay-1)
         ENDDO
      ENDDO





c-----------------------------------------------------------------------
c   6. inversion pour l'implicite sur h sans oublier le couplage
c      avec le sol (conduction)
c      ------------------------

c    ** l'equation est 
c       h(t+1) =  h(t) + dt * {(dh/dt)phys}(t)  +  dt * {(dh/dt)difv}(t+1)
c       avec
c       /zh/ = h(t) + dt * {(dh/dt)phys}(t)   (voir paragraphe 2.)
c       et
c       dt * {(dh/dt)difv}(t+1) = dt * {(d/dz)[ Kh (dh/dz) ]}(t+1)
c       donc les entrees sont /zcdh/ pour la condition de raccord au sol
c       et /zkh/ = Kh
c       -------------

c Mass variation scheme:
      CALL multipl((nlay-1)*ngrid,zkh(1,2),zb0(1,2),zb(1,2))
      CALL multipl(ngrid,zcdh,zb0,zb)

c on initialise dm c
      
      pdtc(:,:)=0.
      zt(:,:)=0.
      dmice(:,:)=0.

c    ** calcul de (d Planck / dT) a la temperature d'interface
c       ------------------------------------------------------

      z4st=4.*5.67e-8*ptimestep
      IF (tke_heat_flux .eq. 0.) THEN
      DO ig=1,ngrid
         zdplanck(ig)=z4st*pemis(ig)*ptsrf(ig)*ptsrf(ig)*ptsrf(ig)
      ENDDO
      ELSE
         zdplanck(:)=0.
      ENDIF

! calcul de zc et zd pour la couche top en prenant en compte le terme
! de variation de masse (on fait une boucle pour que ça converge)

! Identification des points de grilles qui ont besoin de la correction

      llnt(:)=1
#ifdef MESOSCALE
      IF (.not.flag_LES) THEN
#endif
      DO ig=1,ngrid
         DO l=1,nlay
            if(zh(ig,l) .lt. zhcond(ig,l)) then
               llnt(ig)=300  
! 200 and 100 do not go beyond month 9 with normal dissipation
               goto 5
            endif
         ENDDO
5     continue
      ENDDO

#ifdef MESOSCALE
      ENDIF
#endif

      DO ig=1,ngrid

! Initialization of z1 and zd, which do not depend on dmice

      z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
      zd(ig,nlay)=zb(ig,nlay)*z1(ig)

      DO ilay=nlay-1,1,-1
          z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $        zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
          zd(ig,ilay)=zb(ig,ilay)*z1(ig)
      ENDDO

! Convergence loop

      DO j=1,llnt(ig)

            z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
            zc(ig,nlay)=za(ig,nlay)*zh(ig,nlay)
     &      -betam(ig,nlay)*dmice(ig,nlay)
            zc(ig,nlay)=zc(ig,nlay)*z1(ig)
!            zd(ig,nlay)=zb(ig,nlay)*z1(ig)

! calcul de zc et zd pour les couches du haut vers le bas

             DO ilay=nlay-1,1,-1
               z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $            zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
               zc(ig,ilay)=(za(ig,ilay)*zh(ig,ilay)+
     $            zb(ig,ilay+1)*zc(ig,ilay+1)-
     $            betam(ig,ilay)*dmice(ig,ilay))*z1(ig)
!               zd(ig,ilay)=zb(ig,ilay)*z1(ig)
            ENDDO

c    ** calcul de la temperature_d'interface et de sa tendance.
c       on ecrit que la somme des flux est nulle a l'interface
c       a t + \delta t,
c       c'est a dire le flux radiatif a {t + \delta t}
c       + le flux turbulent a {t + \delta t} 
c            qui s'ecrit K (T1-Tsurf) avec T1 = d1 Tsurf + c1
c            (notation K dt = /cpp*b/)        
c       + le flux dans le sol a t
c       + l'evolution du flux dans le sol lorsque la temperature d'interface
c       passe de sa valeur a t a sa valeur a {t + \delta t}.
c       ----------------------------------------------------

         z1(ig)=pcapcal(ig)*ptsrf(ig)+cpp*zb(ig,1)*zc(ig,1)
     s     +zdplanck(ig)*ptsrf(ig)+ pfluxsrf(ig)*ptimestep
         z2(ig)= pcapcal(ig)+cpp*zb(ig,1)*(1.-zd(ig,1))+zdplanck(ig)
         ztsrf2(ig)=z1(ig)/z2(ig)
!         pdtsrf(ig)=(ztsrf2(ig)-ptsrf(ig))/ptimestep  !incremented outside loop
            zhs(ig,1)=zc(ig,1)+zd(ig,1)*ztsrf2(ig)

c    ** et a partir de la temperature au sol on remonte 
c       -----------------------------------------------

            DO ilay=2,nlay
               zhs(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zhs(ig,ilay-1)
            ENDDO
            DO ilay=1,nlay
               zt(ig,ilay)=zhs(ig,ilay)*ppopsk(ig,ilay)
            ENDDO

c      Condensation/sublimation in the atmosphere
c      ------------------------------------------
c      (computation of zcondicea and dmice)

      zcondicea(ig,:)=0.
      pdtc(ig,:)=0.

      DO l=nlay , 1, -1
           IF(zt(ig,l).LT.ztcond(ig,l)) THEN
               pdtc(ig,l)=(ztcond(ig,l) - zt(ig,l))/ptimestep
               zcondicea(ig,l)=(pplev(ig,l)-pplev(ig,l+1))
     &                        *ccond*pdtc(ig,l)
              dmice(ig,l)= dmice(ig,l) + zcondicea(ig,l)*ptimestep
            END IF
      ENDDO

         ENDDO    !of Do j=1,XXX 

      ENDDO   !of Do ig=1,nlay

      pdtsrf(:)=(ztsrf2(:)-ptsrf(:))/ptimestep
      
      DO ig=1,ngrid  ! computing sensible heat flux (atm => surface)
         sensibFlux(ig)=cpp*zb(ig,1)/ptimestep*(zhs(ig,1)-ztsrf2(ig))
      ENDDO

c-----------------------------------------------------------------------
c   TRACERS
c   -------

      if(tracer) then

c     Using the wind modified by friction for lifting and  sublimation
c     ----------------------------------------------------------------

!     This is computed above and takes into account surface-atmosphere flux 
!     enhancement by subgrid gustiness and atmospheric-stability related
!     variations of transfer coefficients.

!        DO ig=1,ngrid
!          zu2(ig)=zu(ig,1)*zu(ig,1)+zv(ig,1)*zv(ig,1)
!          zcdv(ig)=zcdv_true(ig)*sqrt(zu2(ig))
!          zcdh(ig)=zcdh_true(ig)*sqrt(zu2(ig))
!        ENDDO

c       Calcul du flux vertical au bas de la premiere couche (dust) :
c       -----------------------------------------------------------
        do ig=1,ngridmx  
          rho(ig) = zb0(ig,1) /ptimestep
c          zb(ig,1) = 0.
        end do
c       Dust lifting:
        if (lifting) then
#ifndef MESOSCALE
           if (doubleq.AND.submicron) then
             do ig=1,ngrid
c              if(co2ice(ig).lt.1) then
                 pdqsdif(ig,igcm_dust_mass) =
     &             -alpha_lift(igcm_dust_mass)  
                 pdqsdif(ig,igcm_dust_number) = 
     &             -alpha_lift(igcm_dust_number)  
                 pdqsdif(ig,igcm_dust_submicron) =
     &             -alpha_lift(igcm_dust_submicron)
c              end if
             end do
           else if (doubleq) then
             do ig=1,ngrid
               if(co2ice(ig).lt.1) then ! soulevement pas constant
                 pdqsdif(ig,igcm_dust_mass) =
     &             -alpha_lift(igcm_dust_mass)  
                 pdqsdif(ig,igcm_dust_number) = 
     &             -alpha_lift(igcm_dust_number) 
               end if 
             end do
           else if (submicron) then
             do ig=1,ngrid
                 pdqsdif(ig,igcm_dust_submicron) =
     &             -alpha_lift(igcm_dust_submicron)
             end do
           else
            call dustlift(ngrid,nlay,nq,rho,zcdh_true,zcdh,co2ice,
     &                    pdqsdif)
           endif !doubleq.AND.submicron
#else
            call dustlift(ngrid,nlay,nq,rho,zcdh_true,zcdh,co2ice,
     &                    pdqsdif)
#endif
        else
           pdqsdif(1:ngrid,1:nq) = 0.
        end if

c       OU calcul de la valeur de q a la surface (water)  :
c       ----------------------------------------
        if (water) then 
            call watersat(ngridmx,ptsrf,pplev(1,1),qsat)
        end if

c      Inversion pour l'implicite sur q 
c       --------------------------------
        do iq=1,nq
          CALL multipl((nlay-1)*ngrid,zkh(1,2),zb0(1,2),zb(1,2))

          if ((water).and.(iq.eq.igcm_h2o_vap)) then 
c            This line is required to account for turbulent transport 
c            from surface (e.g. ice) to mid-layer of atmosphere:
             CALL multipl(ngrid,zcdv,zb0,zb(1,1))
             CALL multipl(ngrid,dryness,zb(1,1),zb(1,1)) 
          else ! (re)-initialize zb(:,1)
             zb(1:ngrid,1)=0
          end if

          DO ig=1,ngrid
               z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
               zc(ig,nlay)=za(ig,nlay)*zq(ig,nlay,iq)*z1(ig)
               zd(ig,nlay)=zb(ig,nlay)*z1(ig)
          ENDDO
  
          DO ilay=nlay-1,2,-1
               DO ig=1,ngrid
                z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $           zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
                zc(ig,ilay)=(za(ig,ilay)*zq(ig,ilay,iq)+
     $           zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
                zd(ig,ilay)=zb(ig,ilay)*z1(ig)
               ENDDO
          ENDDO

          if (water.and.(iq.eq.igcm_h2o_ice)) then
            ! special case for water ice tracer: do not include
            ! h2o ice tracer from surface (which is set when handling
           ! h2o vapour case (see further down).
            DO ig=1,ngrid
                z1(ig)=1./(za(ig,1)+zb(ig,1)+
     $           zb(ig,2)*(1.-zd(ig,2)))
                zc(ig,1)=(za(ig,1)*zq(ig,1,iq)+
     $         zb(ig,2)*zc(ig,2)) *z1(ig)
            ENDDO
          else ! general case
            DO ig=1,ngrid
                z1(ig)=1./(za(ig,1)+zb(ig,1)+
     $           zb(ig,2)*(1.-zd(ig,2)))
                zc(ig,1)=(za(ig,1)*zq(ig,1,iq)+
     $         zb(ig,2)*zc(ig,2) +
     $        (-pdqsdif(ig,iq)) *ptimestep) *z1(ig)  !tracer flux from surface
            ENDDO
          endif ! of if (water.and.(iq.eq.igcm_h2o_ice))
  
          IF ((water).and.(iq.eq.igcm_h2o_vap)) then 
c           Calculation for turbulent exchange with the surface (for ice)
            DO ig=1,ngrid
              zd(ig,1)=zb(ig,1)*z1(ig)
              zq1temp(ig)=zc(ig,1)+ zd(ig,1)*qsat(ig)

              pdqsdif(ig,igcm_h2o_ice)=rho(ig)*dryness(ig)*zcdv(ig)
     &                       *(zq1temp(ig)-qsat(ig))
c             write(*,*)'flux vers le sol=',pdqsdif(ig,nq)
            END DO

            DO ig=1,ngrid
              if(.not.watercaptag(ig)) then
                if ((-pdqsdif(ig,igcm_h2o_ice)*ptimestep)
     &             .gt.pqsurf(ig,igcm_h2o_ice)) then
c                 write(*,*)'on sublime plus que qsurf!'
                  pdqsdif(ig,igcm_h2o_ice)=
     &                         -pqsurf(ig,igcm_h2o_ice)/ptimestep
c                 write(*,*)'flux vers le sol=',pdqsdif(ig,nq)
                  z1(ig)=1./(za(ig,1)+ zb(ig,2)*(1.-zd(ig,2)))
                  zc(ig,1)=(za(ig,1)*zq(ig,1,igcm_h2o_vap)+
     $            zb(ig,2)*zc(ig,2) +
     $            (-pdqsdif(ig,igcm_h2o_ice)) *ptimestep) *z1(ig)
                  zq1temp(ig)=zc(ig,1)
                endif   
              endif ! if (.not.watercaptag(ig))
c             Starting upward calculations for water :
              zq(ig,1,igcm_h2o_vap)=zq1temp(ig)
              
!c             Take into account H2O latent heat in surface energy budget
!c             We solve dT/dt = (2834.3-0.28*(T-To)-0.004*(T-To)^2)*1e3*iceflux/cpp
!              tsrf_lw(ig) = ptsrf(ig) + pdtsrf(ig) *ptimestep
!             
!              alpha = exp(-4*abs(T1-T2)*pdqsdif(ig,igcm_h2o_ice)
!     &            *ptimestep/pcapcal(ig))
!
!              tsrf_lw(ig) = (tsrf_lw(ig)*(T2-alpha*T1)+T1*T2*(alpha-1))
!     &                      /(tsrf_lw(ig)*(1-alpha)+alpha*T2-T1)  ! surface temperature at t+1
!
!              pdtsrf(ig) = (tsrf_lw(ig)-ptsrf(ig))/ptimestep

               if(pqsurf(ig,igcm_h2o_ice)
     &           +pdqsdif(ig,igcm_h2o_ice)*ptimestep
     &           .gt.frost_albedo_threshold) ! if there is still ice, T cannot exceed To
     &           pdtsrf(ig) = min(pdtsrf(ig),(To-ptsrf(ig))/ptimestep) ! ice melt case
     
            ENDDO ! of DO ig=1,ngrid
          ELSE
c           Starting upward calculations for simple mixing of tracer (dust)
            DO ig=1,ngrid
               zq(ig,1,iq)=zc(ig,1)
            ENDDO
          END IF ! of IF ((water).and.(iq.eq.igcm_h2o_vap))

          DO ilay=2,nlay
             DO ig=1,ngrid
               zq(ig,ilay,iq)=zc(ig,ilay)+zd(ig,ilay)*zq(ig,ilay-1,iq)
             ENDDO
          ENDDO


!         -----------------------------------------------------------
! 	      Tangent Linear Model calculations for vertical diffusion of
!         trace gas species 
!         -----------------------------------------------------------
C           IF (tlm_on) THEN
C                call TLM_vdif(iq, za(1,:), zb(1,:), zc(1,:),
C      &				zd(1,:),ptimestep)
C           ENDIF



        enddo ! of do iq=1,nq
      end if ! of if(tracer)
      



c-----------------------------------------------------------------------
c   8. calcul final des tendances de la diffusion verticale
c      ----------------------------------------------------

      DO ilev = 1, nlay
         DO ig=1,ngrid
            pdudif(ig,ilev)=(    zu(ig,ilev)-
     $      (pu(ig,ilev)+pdufi(ig,ilev)*ptimestep)    )/ptimestep
            pdvdif(ig,ilev)=(    zv(ig,ilev)-
     $      (pv(ig,ilev)+pdvfi(ig,ilev)*ptimestep)    )/ptimestep
            hh = ph(ig,ilev)+pdhfi(ig,ilev)*ptimestep 
     $  + (latcond*dmice(ig,ilev)/cpp)/ppopsk(ig,ilev)
            pdhdif(ig,ilev)=( zhs(ig,ilev)- hh )/ptimestep
         ENDDO
      ENDDO

      if (tracer) then 
        DO iq = 1, nq
          DO ilev = 1, nlay
            DO ig=1,ngrid
              pdqdif(ig,ilev,iq)=(zq(ig,ilev,iq)-
     $      (pq(ig,ilev,iq) + pdqfi(ig,ilev,iq)*ptimestep))/ptimestep
            ENDDO
          ENDDO
        ENDDO
      end if

c    ** diagnostique final 
c       ------------------

      IF(lecrit) THEN
         PRINT*,'In vdif'
         PRINT*,'Ts (t) and Ts (t+st)'
         WRITE(*,'(a10,3a15)')
     s   'theta(t)','theta(t+dt)','u(t)','u(t+dt)'
         PRINT*,ptsrf(ngrid/2+1),ztsrf2(ngrid/2+1)
         DO ilev=1,nlay
            WRITE(*,'(4f15.7)')
     s      ph(ngrid/2+1,ilev),zhs(ngrid/2+1,ilev),
     s      pu(ngrid/2+1,ilev),zu(ngrid/2+1,ilev)

         ENDDO
      ENDIF

      RETURN
      END
