!
!
      SUBROUTINE thermcell_main_mars(ptimestep  &
     &                  ,pplay,pplev,pphi,zlev,zlay  &
     &                  ,pu,pv,pt,pq,pq2  &
     &                  ,pduadj,pdvadj,pdtadj,pdqadj,pdq2adj  &
     &                  ,fm,entr,detr,lmax,zmax  &
     &                  ,r_aspect &
     &                  ,zw2,fraca &
     &                  ,zpopsk,ztla,heatFlux,heatFlux_down &
     &                  ,buoyancyOut, buoyancyEst)

      IMPLICIT NONE

!=======================================================================
! Mars version of thermcell_main. Author : A Colaitis
!=======================================================================

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "tracer.h"
#include "callkeys.h"

! ============== INPUTS ==============

      REAL, INTENT(IN) :: ptimestep,r_aspect
      REAL, INTENT(IN) :: pt(ngridmx,nlayermx)
      REAL, INTENT(IN) :: pu(ngridmx,nlayermx)
      REAL, INTENT(IN) :: pv(ngridmx,nlayermx)
      REAL, INTENT(IN) :: pq(ngridmx,nlayermx,nqmx)
      REAL, INTENT(IN) :: pq2(ngridmx,nlayermx)
      REAL, INTENT(IN) :: pplay(ngridmx,nlayermx)
      REAL, INTENT(IN) :: pplev(ngridmx,nlayermx+1)
      REAL, INTENT(IN) :: pphi(ngridmx,nlayermx)
      REAL, INTENT(IN) :: zlay(ngridmx,nlayermx)
      REAL, INTENT(IN) :: zlev(ngridmx,nlayermx+1)

! ============== OUTPUTS ==============

      REAL, INTENT(OUT) :: pdtadj(ngridmx,nlayermx)
      REAL :: pduadj(ngridmx,nlayermx)
      REAL :: pdvadj(ngridmx,nlayermx)
      REAL :: pdqadj(ngridmx,nlayermx,nqmx)
!      REAL, INTENT(OUT) :: pdq2adj(ngridmx,nlayermx)
      REAL :: pdq2adj(ngridmx,nlayermx)
      REAL, INTENT(OUT) :: zw2(ngridmx,nlayermx+1)

! Diagnostics
      REAL, INTENT(OUT) :: heatFlux(ngridmx,nlayermx)   ! interface heatflux
     REAL, INTENT(OUT) :: heatFlux_down(ngridmx,nlayermx) ! interface heat flux from downdraft
!      REAL, INTENT(OUT) :: buoyancyOut(ngridmx,nlayermx)  ! interlayer buoyancy term
!      REAL, INTENT(OUT) :: buoyancyEst(ngridmx,nlayermx)  ! interlayer estimated buoyancy term

! dummy variables when output not needed :

!      REAL :: heatFlux(ngridmx,nlayermx)   ! interface heatflux
!      REAL :: heatFlux_down(ngridmx,nlayermx) ! interface heat flux from downdraft
      REAL :: buoyancyOut(ngridmx,nlayermx)  ! interlayer buoyancy term
      REAL :: buoyancyEst(ngridmx,nlayermx)  ! interlayer estimated buoyancy term


! ============== LOCAL ================

      INTEGER ig,k,l,ll,iq
      INTEGER lmax(ngridmx),lmin(ngridmx),lalim(ngridmx)
      REAL linter(ngridmx)
      REAL zmax(ngridmx)
      REAL ztva(ngridmx,nlayermx),zw_est(ngridmx,nlayermx+1),ztva_est(ngridmx,nlayermx)
      REAL zmax_sec(ngridmx)
      REAL zh(ngridmx,nlayermx)
      REAL zdthladj(ngridmx,nlayermx)
      REAL zdthladj_down(ngridmx,nlayermx)
      REAL ztvd(ngridmx,nlayermx)
      REAL ztv(ngridmx,nlayermx)
      REAL zu(ngridmx,nlayermx),zv(ngridmx,nlayermx),zo(ngridmx,nlayermx)
      REAL zva(ngridmx,nlayermx)
      REAL zua(ngridmx,nlayermx)

      REAL zta(ngridmx,nlayermx)
      REAL fraca(ngridmx,nlayermx+1)
      REAL q2(ngridmx,nlayermx)
      REAL rho(ngridmx,nlayermx),rhobarz(ngridmx,nlayermx),masse(ngridmx,nlayermx)
      REAL zpopsk(ngridmx,nlayermx)

      REAL wmax(ngridmx)
      REAL wmax_sec(ngridmx)
      REAL fm(ngridmx,nlayermx+1),entr(ngridmx,nlayermx),detr(ngridmx,nlayermx)

      REAL fm_down(ngridmx,nlayermx+1)

      REAL ztla(ngridmx,nlayermx)

      REAL f_star(ngridmx,nlayermx+1),entr_star(ngridmx,nlayermx)
      REAL detr_star(ngridmx,nlayermx)
      REAL alim_star_tot(ngridmx)
      REAL alim_star(ngridmx,nlayermx)
      REAL alim_star_clos(ngridmx,nlayermx)
      REAL f(ngridmx)

      REAL detrmod(ngridmx,nlayermx)

      REAL teta_th_int(ngridmx,nlayermx)
      REAL teta_env_int(ngridmx,nlayermx)
      REAL teta_down_int(ngridmx,nlayermx)

      CHARACTER (LEN=80) :: abort_message
      INTEGER ndt

! ============= PLUME VARIABLES ============

      REAL w_est(ngridmx,nlayermx+1)
      REAL wa_moy(ngridmx,nlayermx+1)
      REAL wmaxa(ngridmx)
      REAL zdz,zbuoy(ngridmx,nlayermx),zw2m
      LOGICAL activecell(ngridmx),activetmp(ngridmx)
      REAL a1,b1,ae,be,ad,bd,fdfu,b1inv,a1inv,omega,adalim
      INTEGER tic

! ==========================================

! ============= HEIGHT VARIABLES ===========

      REAL num(ngridmx)
      REAL denom(ngridmx)
      REAL zlevinter(ngridmx)
      INTEGER zlmax

! =========================================

! ============= CLOSURE VARIABLES =========

      REAL zdenom(ngridmx)
      REAL alim_star2(ngridmx)
      REAL alim_star_tot_clos(ngridmx)
      INTEGER llmax

! =========================================

! ============= FLUX2 VARIABLES ===========

      INTEGER ncorecfm1,ncorecfm2,ncorecfm3,ncorecalpha
      INTEGER ncorecfm4,ncorecfm5,ncorecfm6,ncorecfm7,ncorecfm8
      REAL zfm
      REAL f_old,ddd0,eee0,ddd,eee,zzz
      REAL fomass_max,alphamax

! =========================================

! ============= DTETA VARIABLES ===========

! rien : on prends la divergence du flux turbulent

! =========================================

! ============= DV2 VARIABLES =============
!               not used for now

      REAL qa(ngridmx,nlayermx),detr_dv2(ngridmx,nlayermx),zf,zf2
      REAL wvd(ngridmx,nlayermx+1),wud(ngridmx,nlayermx+1)
      REAL gamma0(ngridmx,nlayermx+1),gamma(ngridmx,nlayermx+1)
      REAL ue(ngridmx,nlayermx),ve(ngridmx,nlayermx)
      LOGICAL ltherm(ngridmx,nlayermx)
      REAL dua(ngridmx,nlayermx),dva(ngridmx,nlayermx)
      INTEGER iter
      INTEGER nlarga0

! =========================================

! ============== Theta_M Variables ========

      INTEGER ico2
      SAVE ico2
      REAL m_co2, m_noco2, A , B
      SAVE A, B
      LOGICAL firstcall
      save firstcall
      data firstcall/.true./
      REAL zhc(ngridmx,nlayermx)
      REAL ratiom(ngridmx,nlayermx)

! =========================================

!-----------------------------------------------------------------------
!   initialisation:
!   ---------------

      entr(:,:)=0.
      detr(:,:)=0.
      fm(:,:)=0.
!      zu(:,:)=pu(:,:)
!      zv(:,:)=pv(:,:)
      zhc(:,:)=pt(:,:)/zpopsk(:,:)
      ndt=1

! **********************************************************************
! Taking into account vertical molar mass gradients
! **********************************************************************

      if(firstcall) then
        ico2=0
        if (tracer) then
!     Prepare Special treatment if one of the tracers is CO2 gas
           do iq=1,nqmx
             if (noms(iq).eq."co2") then
                ico2=iq
                m_co2 = 44.01E-3  ! CO2 molecular mass (kg/mol)
                m_noco2 = 33.37E-3  ! Non condensible mol mass (kg/mol)
!               Compute A and B coefficient use to compute
!               mean molecular mass Mair defined by
!               1/Mair = q(ico2)/m_co2 + (1-q(ico2))/m_noco2
!               1/Mair = A*q(ico2) + B
                A =(1/m_co2 - 1/m_noco2)
                B=1/m_noco2
             end if
           enddo
        endif

      firstcall=.false.
      endif !of if firstcall

!.......................................................................
!  Special treatment for co2
!.......................................................................

      if (ico2.ne.0) then
!     Special case if one of the tracers is CO2 gas
         DO l=1,nlayermx
           DO ig=1,ngridmx
            ztv(ig,l) = zhc(ig,l)*(A*pq(ig,l,ico2)+B)
           ENDDO
         ENDDO
       else
          ztv(:,:)=zhc(:,:)
       end if


!------------------------------------------------------------------------
!                       --------------------
!
!
!                       + + + + + + + + + + +
!
!
!  wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
!  wh,wt,wo ...
!
!                       + + + + + + + + + + +  zh,zu,zv,zo,rho
!
!
!                       --------------------   zlev(1)
!                       \\\\\\\\\\\\\\\\\\\\
!
!

!-----------------------------------------------------------------------
!   Calcul des altitudes des couches
!-----------------------------------------------------------------------

!      do l=2,nlayermx
!         zlev(:,l)=0.5*(pphi(:,l)+pphi(:,l-1))/g
!      enddo
!         zlev(:,1)=0.
!         zlev(:,nlayermx+1)=(2.*pphi(:,nlayermx)-pphi(:,nlayermx-1))/g

!         zlay(:,:)=pphi(:,:)/g
!-----------------------------------------------------------------------
!   Calcul des densites
!-----------------------------------------------------------------------

      rho(:,:)=pplay(:,:)/(r*pt(:,:))

      rhobarz(:,1)=rho(:,1)

      do l=2,nlayermx
!         rhobarz(:,l)=0.5*(rho(:,l)+rho(:,l-1))
          rhobarz(:,l)=pplev(:,l)/(r*0.5*(pt(:,l)+pt(:,l-1)))
      enddo

!calcul de la masse
      do l=1,nlayermx
         masse(:,l)=(pplev(:,l)-pplev(:,l+1))/g
      enddo


!------------------------------------------------------------------
!
!             /|\
!    --------  |  F_k+1 -------   
!                              ----> D_k
!             /|\              <---- E_k , A_k
!    --------  |  F_k --------- 
!                              ----> D_k-1
!                              <---- E_k-1 , A_k-1
!
!
!    ---------------------------
!
!    ----- F_lmax+1=0 ----------         \
!            lmax     (zmax)              |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |  E
!    ---------------------------          |  D
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------  \       |
!            lalim                 |      |
!    ---------------------------   |      |
!                                  |      |
!    ---------------------------   |      |
!                                  | A    |
!    ---------------------------   |      |
!                                  |      |
!    ---------------------------   |      |
!    lmin  (=1 pour le moment)     |      |
!    ----- F_lmin=0 ------------  /      /
!
!    ---------------------------
!    //////////////////////////
!

!=============================================================================
!  Calculs initiaux ne faisant pas intervenir les changements de phase
!=============================================================================

!------------------------------------------------------------------
!  1. alim_star est le profil vertical de l'alimentation a la base du
!     panache thermique, calcule a partir de la flotabilite de l'air sec
!  2. lmin et lalim sont les indices inferieurs et superieurs de alim_star
!------------------------------------------------------------------
!
      entr_star=0. ; detr_star=0. ; alim_star=0. ; alim_star_tot=0.
      lmin=1

!-----------------------------------------------------------------------------
!  3. wmax_sec et zmax_sec sont les vitesses et altitudes maximum d'un
!     panache sec conservatif (e=d=0) alimente selon alim_star 
!     Il s'agit d'un calcul de type CAPE
!     zmax_sec est utilise pour determiner la geometrie du thermique.
!------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!calcul du melange et des variables dans le thermique
!--------------------------------------------------------------------------------

! ===========================================================================
! ===================== PLUME ===============================================
! ===========================================================================

! Initialisations des variables reeles
      ztva(:,:)=ztv(:,:)
      ztva_est(:,:)=ztva(:,:)
      ztla(:,:)=0.
      zdz=0.
      zbuoy(:,:)=0.
      w_est(:,:)=0.
      f_star(:,:)=0.
      wa_moy(:,:)=0.
      linter(:)=1.

! --------------------------------------------------------------------------
! -------------- MAIN PARAMETERS FOR THERMALS MODEL ------------------------
! --------------  see thermiques.pro and getfit.py -------------------------

!      a1=2.5 ; b1=0.0015 ; ae=0.045 ; be = 0.6  ! svn baseline

! Using broad downdraft selection
!      a1=1.60226 ; b1=0.0006 ; ae=0.0454 ; be = 0.57
!      ad = 0.0005114  ; bd = -0.662
!      fdfu = -1.9

! Using conditional sampling downdraft selection
      a1=1.4716 ; b1=0.0005698 ; ae=0.03683 ; be = 0.57421
      ad = 0.00048088  ; bd = -0.6697
!      fdfu = -1.3
      fdfu=-0.8
      a1inv=a1
      b1inv=b1
      omega=0.
      adalim=0.

! One good config for 34/35 levels
!      a1inv=a1
!      b1inv=b1
!      be=1.1*be

! Best configuration for 222 levels:

!      omega=0.06
!      b1=0.
!      a1=1.
!      a1inv=0.25*a1
!      b1inv=0.0002
!!
!!
!!      ae=0.9*ae

! Best config for norad 222 levels:
! and more specifically to C.large :

!       omega=0.06
!       omega=0.
       omega=-0.03
!       omega=0.
       a1=1.
!       b1=0.
       b1=0.0001
       a1inv=a1
       be=1.1*be
       ad = 0.0004
!       ad=0.0003
!       adalim=0.

!       b1inv=0.00025
       b1inv=b1

!       b1inv = 0.0003

!      omega=0.06
! Trying stuff :

!      ad=0.00035
!      ae=0.95*ae
!      b1=0.00055
!      omega=0.04
!
!      ad = 0.0003
!      ae=0.9*ae

!      omega=0.04
!!      b1=0.
!      a1=1.
!      a1inv=a1
!      b1inv=0.0005689
!!      be=1.1*be
!!      ae=0.96*ae


!       omega=0.06
!       a1=1.
!       b1=0.
!       be=be
!       a1inv=0.25*a1
!       b1inv=0.0002    
!       ad=1.1*ad
!       ae=1.*ae
! --------------------------------------------------------------------------
! --------------------------------------------------------------------------
! --------------------------------------------------------------------------

! Initialisation des variables entieres
      wmaxa(:)=0.
      lalim(:)=1

!-------------------------------------------------------------------------
! On ne considere comme actif que les colonnes dont les deux premieres
! couches sont instables.
!-------------------------------------------------------------------------
      activecell(:)=ztv(:,1)>ztv(:,2)
          do ig=1,ngridmx
            if (ztv(ig,1)>=(ztv(ig,2))) then
               alim_star(ig,1)=MAX((ztv(ig,1)-ztv(ig,2)),0.)  &
!     &                       *log(1.+zlev(ig,2))
     &                       *sqrt(zlev(ig,2))
!     &                       *sqrt(sqrt(zlev(ig,2)))
!     &                       /sqrt(zlev(ig,2))
!      &                       *zlev(ig,2)
!      &                     *exp(-zlev(ig,2)/1000.)
               lalim(ig)=2
               alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,1)
            endif
         enddo

       do l=2,nlayermx-1
!        do l=2,4
         do ig=1,ngridmx
           if (ztv(ig,l)>(ztv(ig,l+1)) .and. ztv(ig,1)>=ztv(ig,l) .and. (alim_star(ig,l-1).ne. 0.)) then ! .and. (zlev(ig,l+1) .lt. 1000.)) then
               alim_star(ig,l)=MAX((ztv(ig,l)-ztv(ig,l+1)),0.)  &
!     &                       *log(1.+zlev(ig,l+1))
     &                       *sqrt(zlev(ig,l+1))
!     &                       *sqrt(sqrt(zlev(ig,l+1)))
!     &                       /sqrt(zlev(ig,l+1))
!      &                       *zlev(ig,l+1)
!      &                     *exp(-zlev(ig,l+1)/1000.)
                lalim(ig)=l+1
               alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
           endif
        enddo
      enddo
      do l=1,nlayermx
         do ig=1,ngridmx
            if (alim_star_tot(ig) > 1.e-10 ) then
               alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
            endif
         enddo
      enddo

      alim_star_tot(:)=1.
!      if(alim_star(1,1) .ne. 0.) then
!      print*, 'lalim star' ,lalim(1)
!      endif

!------------------------------------------------------------------------------
! Calcul dans la premiere couche
! On decide dans cette version que le thermique n'est actif que si la premiere
! couche est instable.
! Pourrait etre change si on veut que le thermiques puisse se dÃ©clencher
! dans une couche l>1
!------------------------------------------------------------------------------

      do ig=1,ngridmx
! Le panache va prendre au debut les caracteristiques de l'air contenu
! dans cette couche.
          if (activecell(ig)) then
          ztla(ig,1)=ztv(ig,1)
!cr: attention, prise en compte de f*(1)=1 => AC : what ? f*(1) =0. ! (d'ou f*(2)=a*(1)
! dans un panache conservatif
          f_star(ig,1)=0.
          
          f_star(ig,2)=alim_star(ig,1)

          zw2(ig,2)=2.*g*(ztv(ig,1)-ztv(ig,2))/ztv(ig,2)  &
      &                     *(zlev(ig,2)-zlev(ig,1))  &
      &                    *0.4*pphi(ig,1)/(pphi(ig,2)-pphi(ig,1))
          w_est(ig,2)=zw2(ig,2)

          endif
      enddo


!==============================================================================
!boucle de calcul de la vitesse verticale dans le thermique
!==============================================================================
      do l=2,nlayermx-1
!==============================================================================


! On decide si le thermique est encore actif ou non
! AFaire : Il faut sans doute ajouter entr_star a alim_star dans ce test
          do ig=1,ngridmx
             activecell(ig)=activecell(ig) &
      &                 .and. zw2(ig,l)>1.e-10 &
      &                 .and. f_star(ig,l)+alim_star(ig,l)>1.e-10
          enddo

!---------------------------------------------------------------------------
! calcul des proprietes thermodynamiques et de la vitesse de la couche l
! sans tenir compte du detrainement et de l'entrainement dans cette
! couche
! C'est a dire qu'on suppose
! ztla(l)=ztla(l-1)
! Ici encore, on doit pouvoir ajouter entr_star (qui peut etre calculer
! avant) a l'alimentation pour avoir un calcul plus propre
!---------------------------------------------------------------------------

          do ig=1,ngridmx
             if(activecell(ig)) then
!                if(l .lt. lalim(ig)) then
!               ztva_est(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+  &
!     &            alim_star(ig,l)*ztv(ig,l))  &
!     &            /(f_star(ig,l)+alim_star(ig,l))
!                else
                ztva_est(ig,l)=ztla(ig,l-1)
!                endif

                zdz=zlev(ig,l+1)-zlev(ig,l)
                zbuoy(ig,l)=g*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)

                if (((a1*zbuoy(ig,l)/w_est(ig,l)-b1) .gt. 0.) .and. (w_est(ig,l) .ne. 0.)) then
                w_est(ig,l+1)=Max(0.0001,w_est(ig,l)+2.*zdz*a1*zbuoy(ig,l)-2.*zdz*w_est(ig,l)*b1 &
     & -2.*(1.-omega)*zdz*w_est(ig,l)*ae*(a1*zbuoy(ig,l)/w_est(ig,l)-b1)**be)
                else
                w_est(ig,l+1)=Max(0.0001,w_est(ig,l)+2.*zdz*a1inv*zbuoy(ig,l)-2.*zdz*w_est(ig,l)*b1inv)
                endif
                if (w_est(ig,l+1).lt.0.) then
                w_est(ig,l+1)=zw2(ig,l)
                endif
             endif
          enddo

!-------------------------------------------------
!calcul des taux d'entrainement et de detrainement
!-------------------------------------------------

      do ig=1,ngridmx
        if (activecell(ig)) then

          zw2m=w_est(ig,l+1)

          if((a1*(zbuoy(ig,l)/zw2m)-b1) .gt. 0.) then
          entr_star(ig,l)=f_star(ig,l)*zdz*  &
        &   MAX(0.,ae*(a1*(zbuoy(ig,l)/zw2m)-b1)**be)
!         &  MAX(0.,log(1. + 0.03*sqrt(a1*(zbuoy(ig,l)/zw2m)-b1)))
          else
          entr_star(ig,l)=0.
          endif

          if(zbuoy(ig,l) .gt. 0.) then
             if(l .lt. lalim(ig)) then

!                detr_star(ig,l)=0.
                 detr_star(ig,l) = f_star(ig,l)*zdz*              &
            &  adalim
             else

!                 detr_star(ig,l) = f_star(ig,l)*zdz*              &
!              &     0.0105*((zbuoy(ig,l)/zw2m)/0.048)**(1./1.7)
!                 detr_star(ig,l) = f_star(ig,l)*zdz*              &
!              &     0.0085*((zbuoy(ig,l)/zw2m)/0.05)**(1./1.55)

! last baseline from direct les
!                 detr_star(ig,l) = f_star(ig,l)*zdz*              &
!               &     0.065*(2.5*(zbuoy(ig,l)/zw2m))**0.75

! new param from continuity eq with a fit on dfdz
                 detr_star(ig,l) = f_star(ig,l)*zdz*              &
            &  ad
!             &  Max(0., 0.0005 - 0.55*zbuoy(ig,l)/zw2m)

!               &     MAX(0.,-0.38*zbuoy(ig,l)/zw2m+0.0005)   !svn baseline
!                &     MAX(0.,-0.38*zbuoy(ig,l)/zw2m+0.0008)     

!              &     0.014*((zbuoy(ig,l)/zw2m)/0.05)**(1./1.35)
!                 detr_star(ig,l) = f_star(ig,l)*zdz*              &
!              &     ((zbuoy(ig,l)/zw2m)/2.222)! + 0.0002)

             endif
          else
          detr_star(ig,l)=f_star(ig,l)*zdz*                        &
                &     MAX(ad,bd*zbuoy(ig,l)/zw2m)
!             &  Max(0., 0.001 - 0.45*zbuoy(ig,l)/zw2m)
!             &  Max(0., Min(0.001,0.0005 - 0.55*zbuoy(ig,l)/zw2m))


!               &     MAX(0.,-0.38*zbuoy(ig,l)/zw2m+0.0005)      !svn baseline

!              &  *5.*(-afact*zbetalpha*zbuoy(ig,l)/zw2m)
!              &  *5.*(-afact*zbuoy(ig,l)/zw2m)

! last baseline from direct les
!               &     0.065*(-2.5*(zbuoy(ig,l)/zw2m))**0.75

! new param from continuity eq with a fit on dfdz


          endif

! En dessous de lalim, on prend le max de alim_star et entr_star pour
! alim_star et 0 sinon

       if (l.lt.lalim(ig)) then
          alim_star(ig,l)=max(alim_star(ig,l),entr_star(ig,l))
          entr_star(ig,l)=0. 
       endif

! Calcul du flux montant normalise

      f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)  &
     &              -detr_star(ig,l)

      endif
      enddo


!----------------------------------------------------------------------------
!calcul de la vitesse verticale en melangeant Tl et qt du thermique
!---------------------------------------------------------------------------

      DO tic=0,5  ! internal convergence loop
      activetmp(:)=activecell(:) .and. f_star(:,l+1)>1.e-10
      do ig=1,ngridmx
       if (activetmp(ig)) then

           ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+  &
     &            (alim_star(ig,l)+entr_star(ig,l))*ztv(ig,l))  &
     &            /(f_star(ig,l+1)+detr_star(ig,l))

        endif
      enddo

      activetmp(:)=activetmp(:).and.(abs(ztla(:,l)-ztva(:,l)).gt.0.01)

      do ig=1,ngridmx
      if (activetmp(ig)) then
           ztva(ig,l) = ztla(ig,l)
           zbuoy(ig,l)=g*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)

           if (((a1*zbuoy(ig,l)/zw2(ig,l)-b1) .gt. 0.) .and. (zw2(ig,l) .ne. 0.) ) then
           zw2(ig,l+1)=Max(0.,zw2(ig,l)+2.*zdz*a1*zbuoy(ig,l)-         &
     & 2.*zdz*zw2(ig,l)*b1-2.*(1.-omega)*zdz*zw2(ig,l)*ae*(a1*zbuoy(ig,l)/zw2(ig,l)-b1)**be)
           else
           zw2(ig,l+1)=Max(0.,zw2(ig,l)+2.*zdz*a1inv*zbuoy(ig,l)-2.*zdz*zw2(ig,l)*b1inv)
           endif
      endif
      enddo

! ================ RECOMPUTE ENTR, DETR, and F FROM NEW W2 ===================

      do ig=1,ngridmx
        if (activetmp(ig)) then

          zw2m=zw2(ig,l+1)
          if(zw2m .gt. 0) then
          if((a1*(zbuoy(ig,l)/zw2m)-b1) .gt. 0.) then
          entr_star(ig,l)=f_star(ig,l)*zdz*  &
        &   MAX(0.,ae*(a1*(zbuoy(ig,l)/zw2m)-b1)**be)
!         &  MAX(0.,log(1. + 0.03*sqrt(a1*(zbuoy(ig,l)/zw2m)-b1)))
          else
          entr_star(ig,l)=0.
          endif

          if(zbuoy(ig,l) .gt. 0.) then
             if(l .lt. lalim(ig)) then

!                detr_star(ig,l)=0.
                 detr_star(ig,l) = f_star(ig,l)*zdz*              &
            &  adalim

             else
                 detr_star(ig,l) = f_star(ig,l)*zdz*              &
            &  ad
!             &  Max(0., 0.0005 - 0.55*zbuoy(ig,l)/zw2m)

             endif
          else
          detr_star(ig,l)=f_star(ig,l)*zdz*                        &
                &     MAX(ad,bd*zbuoy(ig,l)/zw2m)
!             &  Max(0.,Min(0.001,0.0005 - 0.55*zbuoy(ig,l)/zw2m))

          endif
          else
          entr_star(ig,l)=0.
          detr_star(ig,l)=0.
          endif

! En dessous de lalim, on prend le max de alim_star et entr_star pour
! alim_star et 0 sinon

        if (l.lt.lalim(ig)) then
          alim_star(ig,l)=max(alim_star(ig,l),entr_star(ig,l))
          entr_star(ig,l)=0.
        endif

! Calcul du flux montant normalise

      f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)  &
     &              -detr_star(ig,l)

      endif
      enddo
 
      ENDDO   ! of tic

!---------------------------------------------------------------------------
!initialisations pour le calcul de la hauteur du thermique, de l'inversion et de la vitesse verticale max
!---------------------------------------------------------------------------

      do ig=1,ngridmx
            if (zw2(ig,l+1)>0. .and. zw2(ig,l+1).lt.1.e-10) then
      IF (lwrite) THEN
             print*,'On tombe sur le cas particulier de thermcell_plume'
      ENDIF
                zw2(ig,l+1)=0.
                linter(ig)=l+1
            endif

        if (zw2(ig,l+1).lt.0.) then
           linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))  &
     &               -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
           zw2(ig,l+1)=0.
        endif
           wa_moy(ig,l+1)=sqrt(zw2(ig,l+1))

        if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
!   lmix est le niveau de la couche ou w (wa_moy) est maximum
            wmaxa(ig)=wa_moy(ig,l+1)
        endif
      enddo

!=========================================================================
! FIN DE LA BOUCLE VERTICALE
      enddo
!=========================================================================

!on recalcule alim_star_tot
       do ig=1,ngridmx
          alim_star_tot(ig)=0.
       enddo
       do ig=1,ngridmx
          do l=1,lalim(ig)-1
          alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
          enddo
       enddo

      do l=1,nlayermx
         do ig=1,ngridmx
            if (alim_star_tot(ig) > 1.e-10 ) then
               alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
            endif
         enddo
      enddo

! ===========================================================================
! ================= FIN PLUME ===============================================
! ===========================================================================

! ===========================================================================
! ================= HEIGHT ==================================================
! ===========================================================================

! Attention, w2 est transforme en sa racine carree dans cette routine

!-------------------------------------------------------------------------------
! Calcul des caracteristiques du thermique:zmax,wmax
!-------------------------------------------------------------------------------

!calcul de la hauteur max du thermique
      do ig=1,ngridmx
         lmax(ig)=lalim(ig)
      enddo
      do ig=1,ngridmx
         do l=nlayermx,lalim(ig)+1,-1
            if (zw2(ig,l).le.1.e-10) then
               lmax(ig)=l-1 
            endif
         enddo
      enddo

! On traite le cas particulier qu'il faudrait éviter ou le thermique
! atteind le haut du modele ...
      do ig=1,ngridmx
      if ( zw2(ig,nlayermx) > 1.e-10 ) then
          print*,'WARNING !!!!! W2 thermiques non nul derniere couche '
          lmax(ig)=nlayermx
      endif
      enddo

! pas de thermique si couche 1 stable
!      do ig=1,ngridmx
!         if (lmin(ig).gt.1) then
!             lmax(ig)=1
!             lmin(ig)=1
!             lalim(ig)=1
!         endif
!      enddo
!
! Determination de zw2 max
      do ig=1,ngridmx
         wmax(ig)=0.
      enddo

      do l=1,nlayermx
         do ig=1,ngridmx
            if (l.le.lmax(ig)) then
                if (zw2(ig,l).lt.0.)then
!                  print*,'pb2 zw2<0',zw2(ig,l)
                  zw2(ig,l)=0.
                endif
                zw2(ig,l)=sqrt(zw2(ig,l))
                wmax(ig)=max(wmax(ig),zw2(ig,l))
            else
                 zw2(ig,l)=0.
            endif
          enddo
      enddo
!   Longueur caracteristique correspondant a la hauteur des thermiques.
      do  ig=1,ngridmx
         zmax(ig)=0.
         zlevinter(ig)=zlev(ig,1)
      enddo

         num(:)=0.
         denom(:)=0.
         do ig=1,ngridmx
          do l=1,nlayermx
             num(ig)=num(ig)+zw2(ig,l)*zlev(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
             denom(ig)=denom(ig)+zw2(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
          enddo
       enddo
       do ig=1,ngridmx
       if (denom(ig).gt.1.e-10) then
          zmax(ig)=2.*num(ig)/denom(ig)
       endif
       enddo

! Attention, w2 est transforme en sa racine carree dans cette routine

! ===========================================================================
! ================= FIN HEIGHT ==============================================
! ===========================================================================

      zlmax=MAXVAL(lmax(:))+2
      if (zlmax .ge. nlayermx) then
        print*,'thermals have reached last layer of the model'
        print*,'this is not good !'
      endif

! Choix de la fonction d'alimentation utilisee pour la fermeture.

      alim_star_clos(:,:)=entr_star(:,:)+alim_star(:,:)

! ===========================================================================
! ============= CLOSURE =====================================================
! ===========================================================================

!-------------------------------------------------------------------------------
! Fermeture,determination de f
!-------------------------------------------------------------------------------
! Appel avec la version seche

      alim_star2(:)=0.
      alim_star_tot_clos(:)=0.
      f(:)=0.

! Indice vertical max (max de lalim) atteint par les thermiques sur le domaine
      llmax=1
      do ig=1,ngridmx
         if (lalim(ig)>llmax) llmax=lalim(ig)
      enddo


! Calcul des integrales sur la verticale de alim_star et de
!   alim_star^2/(rho dz)
      do k=1,llmax-1
         do ig=1,ngridmx
            if (k<lalim(ig)) then
         alim_star2(ig)=alim_star2(ig)+alim_star_clos(ig,k)*alim_star_clos(ig,k)  &
      &                    /(rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k)))
         alim_star_tot_clos(ig)=alim_star_tot_clos(ig)+alim_star_clos(ig,k)
      endif
         enddo
      enddo
 
! WARNING : MARS MODIF : we have added 2. : ratio of wmax/vmoy
! True ratio is 3.5 but wetake into account the vmoy is the one alimentating
! the thermal, so there are vs=0 into the vmoy... the true vmoy is lower. (a la louche)
! And r_aspect has been changed from 2 to 1.5 from observations
      do ig=1,ngridmx
         if (alim_star2(ig)>1.e-10) then
!            f(ig)=wmax_sec(ig)*alim_star_tot_clos(ig)/  &
!      &     (max(500.,zmax_sec(ig))*r_aspect*alim_star2(ig))
             f(ig)=wmax(ig)*alim_star_tot_clos(ig)/  &
      &     (max(500.,zmax(ig))*r_aspect*alim_star2(ig))

         endif
      enddo

! ===========================================================================
! ============= FIN CLOSURE =================================================
! ===========================================================================


! ===========================================================================
! ============= FLUX2 =======================================================
! ===========================================================================

!-------------------------------------------------------------------------------
!deduction des flux
!-------------------------------------------------------------------------------

      fomass_max=0.8
      alphamax=0.5

      ncorecfm1=0
      ncorecfm2=0
      ncorecfm3=0
      ncorecfm4=0
      ncorecfm5=0
      ncorecfm6=0
      ncorecfm7=0
      ncorecfm8=0
      ncorecalpha=0

!-------------------------------------------------------------------------
! Multiplication par le flux de masse issu de la femreture
!-------------------------------------------------------------------------

      do l=1,zlmax
         entr(:,l)=f(:)*(entr_star(:,l)+alim_star(:,l))
         detr(:,l)=f(:)*detr_star(:,l)
      enddo

      do l=1,zlmax
         do ig=1,ngridmx
            if (l.lt.lmax(ig)) then
               fm(ig,l+1)=fm(ig,l)+entr(ig,l)-detr(ig,l)
            elseif(l.eq.lmax(ig)) then
               fm(ig,l+1)=0.
               detr(ig,l)=fm(ig,l)+entr(ig,l)
            else
               fm(ig,l+1)=0.
            endif
         enddo
      enddo

! Test provisoire : pour comprendre pourquoi on corrige plein de fois
! le cas fm6, on commence par regarder une premiere fois avant les
! autres corrections.

!      do l=1,nlayermx
!         do ig=1,ngridmx
!            if (detr(ig,l).gt.fm(ig,l)) then
!               ncorecfm8=ncorecfm8+1
!            endif
!         enddo
!      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH Version en cours de test;
! par rapport a thermcell_flux, on fait une grande boucle sur "l"
! et on modifie le flux avec tous les contrï¿½les appliques d'affilee
! pour la meme couche
! Momentanement, on duplique le calcule du flux pour pouvoir comparer
! les flux avant et apres modif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do l=1,zlmax

         do ig=1,ngridmx
            if (l.lt.lmax(ig)) then
               fm(ig,l+1)=fm(ig,l)+entr(ig,l)-detr(ig,l)
            elseif(l.eq.lmax(ig)) then
               fm(ig,l+1)=0.
               detr(ig,l)=fm(ig,l)+entr(ig,l)
            else
               fm(ig,l+1)=0.
            endif
         enddo


!-------------------------------------------------------------------------
! Verification de la positivite des flux de masse
!-------------------------------------------------------------------------

         do ig=1,ngridmx

            if (fm(ig,l+1).lt.0.) then
               if((l+1) .eq. lmax(ig)) then
               detr(ig,l)=detr(ig,l)+fm(ig,l+1)
               fm(ig,l+1)=0.
               entr(ig,l+1)=0.
               ncorecfm2=ncorecfm2+1
               else
               IF (lwrite) THEN
          print*,'fm(l+1)<0 : ig, l+1,lmax :',ig,l+1,lmax(ig),fm(ig,l+1)
               ENDIF
               ncorecfm1=ncorecfm1+1
               fm(ig,l+1)=fm(ig,l)
               detr(ig,l)=entr(ig,l)
               endif
            endif

         enddo

!  Les "optimisations" du flux sont desactivecelles : moins de bidouilles
!  je considere que celles ci ne sont pas justifiees ou trop delicates
!  pour MARS, d'apres des observations LES. 
!-------------------------------------------------------------------------
!Test sur fraca croissant
!-------------------------------------------------------------------------
!      if (iflag_thermals_optflux==0) then
!         do ig=1,ngridmx
!          if (l.ge.lalim(ig).and.l.le.lmax(ig) &
!     &    .and.(zw2(ig,l+1).gt.1.e-10).and.(zw2(ig,l).gt.1.e-10) ) then
!!  zzz est le flux en l+1 a frac constant
!             zzz=fm(ig,l)*rhobarz(ig,l+1)*zw2(ig,l+1)  &
!     &                          /(rhobarz(ig,l)*zw2(ig,l))
!             if (fm(ig,l+1).gt.zzz) then
!                detr(ig,l)=detr(ig,l)+fm(ig,l+1)-zzz
!                fm(ig,l+1)=zzz
!                ncorecfm4=ncorecfm4+1
!             endif
!          endif
!        enddo
!      endif
!
!-------------------------------------------------------------------------
!test sur flux de masse croissant
!-------------------------------------------------------------------------
!      if (iflag_thermals_optflux==0) then
!         do ig=1,ngridmx
!            if ((fm(ig,l+1).gt.fm(ig,l)).and.(l.gt.lalim(ig))) then
!              f_old=fm(ig,l+1)
!              fm(ig,l+1)=fm(ig,l)
!              detr(ig,l)=detr(ig,l)+f_old-fm(ig,l+1)
!               ncorecfm5=ncorecfm5+1
!            endif
!         enddo
!      endif
!
!-------------------------------------------------------------------------
!detr ne peut pas etre superieur a fm
!-------------------------------------------------------------------------

         do ig=1,ngridmx
            if (detr(ig,l).gt.fm(ig,l)) then
               ncorecfm6=ncorecfm6+1
               detr(ig,l)=fm(ig,l)
               entr(ig,l)=fm(ig,l+1)

! Dans le cas ou on est au dessus de la couche d'alimentation et que le
! detrainement est plus fort que le flux de masse, on stope le thermique.
!            endif

            if(l.gt.lmax(ig)) then
!            if(l.gt.lalim(ig)) then
               detr(ig,l)=0.
               fm(ig,l+1)=0.
               entr(ig,l)=0.
            endif
            
            endif

         enddo

!-------------------------------------------------------------------------
!fm ne peut pas etre negatif
!-------------------------------------------------------------------------

         do ig=1,ngridmx
            if (fm(ig,l+1).lt.0.) then
               detr(ig,l)=detr(ig,l)+fm(ig,l+1)
               fm(ig,l+1)=0.
               ncorecfm2=ncorecfm2+1
            endif
         enddo

!-----------------------------------------------------------------------
!la fraction couverte ne peut pas etre superieure a 1
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH Partie a revisiter.
! Il semble qu'etaient codees ici deux optiques dans le cas
! F/ (rho *w) > 1
! soit limiter la hauteur du thermique en considerant que c'est
! la derniere chouche, soit limiter F a rho w.
! Dans le second cas, il faut en fait limiter a un peu moins
! que ca parce qu'on a des 1 / ( 1 -alpha) un peu plus loin
! dans thermcell_main et qu'il semble de toutes facons deraisonable
! d'avoir des fractions de 1..
! Ci dessous, et dans l'etat actuel, le premier des  deux if est
! sans doute inutile.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do ig=1,ngridmx
           if (zw2(ig,l+1).gt.1.e-10) then
           zfm=rhobarz(ig,l+1)*zw2(ig,l+1)*alphamax
           if ( fm(ig,l+1) .gt. zfm) then
              f_old=fm(ig,l+1)
              fm(ig,l+1)=zfm
              detr(ig,l)=detr(ig,l)+f_old-fm(ig,l+1)
              ncorecalpha=ncorecalpha+1
           endif
           endif

        enddo

! Fin de la grande boucle sur les niveaux verticaux
      enddo

!-----------------------------------------------------------------------
! On fait en sorte que la quantite totale d'air entraine dans le
! panache ne soit pas trop grande comparee a la masse de la maille
!-----------------------------------------------------------------------

      do l=1,zlmax
         do ig=1,ngridmx
            eee0=entr(ig,l)
            ddd0=detr(ig,l)
            eee=entr(ig,l)-masse(ig,l)*fomass_max/ptimestep
            ddd=detr(ig,l)-eee
            if (eee.gt.0.) then
                ncorecfm3=ncorecfm3+1
                entr(ig,l)=entr(ig,l)-eee
                if ( ddd.gt.0.) then
!   l'entrainement est trop fort mais l'exces peut etre compense par une
!   diminution du detrainement)
                   detr(ig,l)=ddd
                else
!   l'entrainement est trop fort mais l'exces doit etre compense en partie
!   par un entrainement plus fort dans la couche superieure
                   if(l.eq.lmax(ig)) then
                      detr(ig,l)=fm(ig,l)+entr(ig,l)
                   else
                      entr(ig,l+1)=entr(ig,l+1)-ddd
                      detr(ig,l)=0.
                      fm(ig,l+1)=fm(ig,l)+entr(ig,l)
                      detr(ig,l)=0.
                   endif
                endif
            endif
         enddo
      enddo
!
!              ddd=detr(ig)-entre
!on s assure que tout s annule bien en zmax
      do ig=1,ngridmx
         fm(ig,lmax(ig)+1)=0.
         entr(ig,lmax(ig))=0.
         detr(ig,lmax(ig))=fm(ig,lmax(ig))+entr(ig,lmax(ig))
      enddo

!-----------------------------------------------------------------------
! Impression du nombre de bidouilles qui ont ete necessaires
!-----------------------------------------------------------------------

!IM 090508 beg
      IF (lwrite) THEN
      if (ncorecfm1+ncorecfm2+ncorecfm3+ncorecfm4+ncorecfm5+ncorecalpha > ngridmx/4. ) then
         print*,'thermcell warning : large number of corrections'
         print*,'PB thermcell : on a du coriger ',ncorecfm1,'x fm1',&
     &     ncorecfm2,'x fm2',ncorecfm3,'x fm3 et', &
     &     ncorecfm4,'x fm4',ncorecfm5,'x fm5 et', &
     &     ncorecfm6,'x fm6', &
     &     ncorecfm7,'x fm7', &
     &     ncorecfm8,'x fm8', &
     &     ncorecalpha,'x alpha'
      endif
      ENDIF

! ===========================================================================
! ============= FIN FLUX2 ===================================================
! ===========================================================================


! ===========================================================================
! ============= TRANSPORT ===================================================
! ===========================================================================

!------------------------------------------------------------------
!   calcul du transport vertical
!------------------------------------------------------------------

! ------------------------------------------------------------------
! Transport de teta dans l'updraft : (remplace thermcell_dq)
! ------------------------------------------------------------------

      zdthladj(:,:)=0.

      if (1 .eq. 0) then
!      call thermcell_dqup(ngridmx,nlayermx,ptimestep     &
!     &     ,fm,entr,  &
!     &    masse,ztv,zdthladj)
      else


      do ig=1,ngridmx
         if(lmax(ig) .gt. 1) then
         do k=1,lmax(ig)
            zdthladj(ig,k)=(1./masse(ig,k))*(fm(ig,k+1)*ztv(ig,k+1)-    &
     &   fm(ig,k)*ztv(ig,k)+fm(ig,k)*ztva(ig,k)-fm(ig,k+1)*ztva(ig,k+1))
            if (ztv(ig,k) + ptimestep*zdthladj(ig,k) .le. 0.) then
      IF (lwrite) THEN
              print*,'Teta<0 in thermcell_dTeta up: qenv .. dq : ', ztv(ig,k),ptimestep*zdthladj(ig,k)
      ENDIF
              if(ztv(ig,k) .gt. 0.) then
              zdthladj(ig,k)=0.
              endif
            endif
         enddo
         endif
      enddo

      endif

! ------------------------------------------------------------------
! Prescription des proprietes du downdraft
! ------------------------------------------------------------------

      ztvd(:,:)=ztv(:,:)
      fm_down(:,:)=0.
      do ig=1,ngridmx
         if (lmax(ig) .gt. 1) then
         do l=1,lmax(ig)
!              if(zlay(ig,l) .le. 0.8*zmax(ig)) then
              if(zlay(ig,l) .le. zmax(ig)) then
              fm_down(ig,l) =fm(ig,l)* &
     &      max(fdfu,-4*max(0.,(zlay(ig,l)/zmax(ig)))-0.6)
              endif

!             if(zlay(ig,l) .le. 0.06*zmax(ig)) then
!          ztvd(ig,l)=ztv(ig,l)*max(0.,(1.+(sqrt((zlay(ig,l)/zmax(ig))/0.122449) - 1.)*(ztva(ig,l)/ztv(ig,l) - 1.)))
!             elseif(zlay(ig,l) .le. 0.4*zmax(ig)) then
!          ztvd(ig,l)=ztv(ig,l)*max(0.,1.-0.25*(ztva(ig,l)/ztv(ig,l) - 1.))
!             elseif(zlay(ig,l) .le. 0.7*zmax(ig)) then
!          ztvd(ig,l)=ztv(ig,l)*max(0.,(1.+(((zlay(ig,l)/zmax(ig))-0.7)/1.)*(ztva(ig,l)/ztv(ig,l) - 1.)))
!             else
!          ztvd(ig,l)=ztv(ig,l)
!            endif

!            if(zlay(ig,l) .le. 0.6*zmax(ig)) then
!          ztvd(ig,l)=ztv(ig,l)*((zlay(ig,l)/zmax(ig))/179.848 + 0.997832)
!            elseif(zlay(ig,l) .le. 0.8*zmax(ig)) then
!          ztvd(ig,l)=-ztv(ig,l)*(((zlay(ig,l)/zmax(ig))-171.74)/170.94)
!             else
!          ztvd(ig,l)=ztv(ig,l)
!            endif


!            if(zlay(ig,l) .le. 0.8*zmax(ig)) then
!          ztvd(ig,l)=ztv(ig,l)*((zlay(ig,l)/zmax(ig))/224.81 + 0.997832)
!            elseif(zlay(ig,l) .le. zmax(ig)) then
!          ztvd(ig,l)=-ztv(ig,l)*(((zlay(ig,l)/zmax(ig))-144.685)/143.885)
!             else
!          ztvd(ig,l)=ztv(ig,l)
!            endif


!             if (zbuoy(ig,l) .gt. 0.) then
!               ztvd(ig,l)=ztva(ig,l)*0.9998
!!               ztvd(ig,l)=ztv(ig,l)*0.997832
!!             else
!!               if(zlay(ig,l) .le. zmax(ig)) then            
!!               ztvd(ig,l)=ztv(ig,l)*((zlay(ig,l)/zmax(ig))/299.7 + 0.997832)
!!               endif
!             endif

            if(zlay(ig,l) .le. zmax(ig)) then            
          ztvd(ig,l)=min(ztv(ig,l),ztv(ig,l)*((zlay(ig,l)/zmax(ig))/400. + 0.997832))
!          ztvd(ig,l)=min(ztv(ig,l),ztv(ig,l)*((zlay(ig,l)/zmax(ig))/299.7 + 0.997832))
             else
          ztvd(ig,l)=ztv(ig,l)
            endif

         enddo
         endif
      enddo

! ------------------------------------------------------------------
! Transport de teta dans le downdraft : (remplace thermcell_dq)
! ------------------------------------------------------------------

       zdthladj_down(:,:)=0.

      do ig=1,ngridmx
         if(lmax(ig) .gt. 1) then
! No downdraft in the very-near surface layer, we begin at k=3
 
         do k=3,lmax(ig)
            zdthladj_down(ig,k)=(1./masse(ig,k))*(fm_down(ig,k+1)*ztv(ig,k+1)- &
     & fm_down(ig,k)*ztv(ig,k)+fm_down(ig,k)*ztvd(ig,k)-fm_down(ig,k+1)*ztvd(ig,k+1))
            if (ztv(ig,k) + ptimestep*zdthladj_down(ig,k) .le. 0.) then
      IF (lwrite) THEN
              print*,'q<0 in thermcell_dTeta down: qenv .. dq : ', ztv(ig,k),ptimestep*zdthladj_down(ig,k)
      ENDIF
              if(ztv(ig,k) .gt. 0.) then
              zdthladj(ig,k)=0.
              endif
            endif
         enddo
         endif
      enddo

!------------------------------------------------------------------
! Calcul de la fraction de l'ascendance
!------------------------------------------------------------------
      fraca(:,:)=0.
      do l=2,zlmax
         do ig=1,ngridmx
            if (zw2(ig,l).gt.1.e-10) then
            fraca(ig,l)=fm(ig,l)/(rhobarz(ig,l)*zw2(ig,l))
            else
            fraca(ig,l)=0.
            endif
         enddo
      enddo



! ===========================================================================
! ============= DV2 =========================================================
! ===========================================================================
! ROUTINE OVERIDE : ne prends pas en compte le downdraft
! de plus, le gradient de pression horizontal semble tout deregler... A VOIR

      if (0 .eq. 1) then

!------------------------------------------------------------------
!  calcul du transport vertical du moment horizontal
!------------------------------------------------------------------

! Calcul du transport de V tenant compte d'echange par gradient
! de pression horizontal avec l'environnement

!   calcul du detrainement
!---------------------------

      nlarga0=0.

      do k=1,nlayermx
         do ig=1,ngridmx
            detr_dv2(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
         enddo
      enddo

!   calcul de la valeur dans les ascendances
      do ig=1,ngridmx
         zua(ig,1)=zu(ig,1)
         zva(ig,1)=zv(ig,1)
         ue(ig,1)=zu(ig,1)
         ve(ig,1)=zv(ig,1)
      enddo

      gamma(1:ngridmx,1)=0.
      do k=2,nlayermx
         do ig=1,ngridmx
            ltherm(ig,k)=(fm(ig,k+1)+detr_dv2(ig,k))*ptimestep > 1.e-5*masse(ig,k)
            if(ltherm(ig,k).and.zmax(ig)>0.) then
               gamma0(ig,k)=masse(ig,k)  &
     &         *sqrt( 0.5*(fraca(ig,k+1)+fraca(ig,k)) )  &
     &         *0.5/zmax(ig)  &
     &         *1.
            else
               gamma0(ig,k)=0.
            endif
            if (ltherm(ig,k).and.zmax(ig)<=0.) nlarga0=nlarga0+1
          enddo
      enddo

      gamma(:,:)=0.

      do k=2,nlayermx

         do ig=1,ngridmx

            if (ltherm(ig,k)) then
               dua(ig,k)=zua(ig,k-1)-zu(ig,k-1)
               dva(ig,k)=zva(ig,k-1)-zv(ig,k-1)
            else
               zua(ig,k)=zu(ig,k)
               zva(ig,k)=zv(ig,k)
               ue(ig,k)=zu(ig,k)
               ve(ig,k)=zv(ig,k)
            endif
         enddo


! Debut des iterations
!----------------------

! AC WARNING : SALE !

      do iter=1,5
         do ig=1,ngridmx
! Pour memoire : calcul prenant en compte la fraction reelle
!              zf=0.5*(fraca(ig,k)+fraca(ig,k+1))
!              zf2=1./(1.-zf)
! Calcul avec fraction infiniement petite
               zf=0.
               zf2=1.

!  la première fois on multiplie le coefficient de freinage
!  par le module du vent dans la couche en dessous.
!  Mais pourquoi donc ???
               if (ltherm(ig,k)) then
!   On choisit une relaxation lineaire.
!                 gamma(ig,k)=gamma0(ig,k)
!   On choisit une relaxation quadratique.
                gamma(ig,k)=gamma0(ig,k)*sqrt(dua(ig,k)**2+dva(ig,k)**2)
                  zua(ig,k)=(fm(ig,k)*zua(ig,k-1)  &
     &               +(zf2*entr(ig,k)+gamma(ig,k))*zu(ig,k))  &
     &               /(fm(ig,k+1)+detr_dv2(ig,k)+entr(ig,k)*zf*zf2  &
     &                 +gamma(ig,k))
                  zva(ig,k)=(fm(ig,k)*zva(ig,k-1)  &
     &               +(zf2*entr(ig,k)+gamma(ig,k))*zv(ig,k))  &
     &               /(fm(ig,k+1)+detr_dv2(ig,k)+entr(ig,k)*zf*zf2  &
     &                 +gamma(ig,k))

!                  print*,' OUTPUT DV2 !!!!!!!!!!!!!!!!!!!!',k,zua(ig,k),zva(ig,k),zu(ig,k),zv(ig,k),dua(ig,k),dva(ig,k)
                  dua(ig,k)=zua(ig,k)-zu(ig,k)
                  dva(ig,k)=zva(ig,k)-zv(ig,k)
                  ue(ig,k)=(zu(ig,k)-zf*zua(ig,k))*zf2
                  ve(ig,k)=(zv(ig,k)-zf*zva(ig,k))*zf2
               endif
         enddo
! Fin des iterations
!--------------------
      enddo

      enddo ! k=2,nlayermx

! Calcul du flux vertical de moment dans l'environnement.
!---------------------------------------------------------
      do k=2,nlayermx
         do ig=1,ngridmx
            wud(ig,k)=fm(ig,k)*ue(ig,k)
            wvd(ig,k)=fm(ig,k)*ve(ig,k)
         enddo
      enddo
      do ig=1,ngridmx
         wud(ig,1)=0.
         wud(ig,nlayermx+1)=0.
         wvd(ig,1)=0.
         wvd(ig,nlayermx+1)=0.
      enddo

! calcul des tendances.
!-----------------------
      do k=1,nlayermx
         do ig=1,ngridmx
            pduadj(ig,k)=((detr_dv2(ig,k)+gamma(ig,k))*zua(ig,k)  &
     &               -(entr(ig,k)+gamma(ig,k))*ue(ig,k)  &
     &               -wud(ig,k)+wud(ig,k+1))  &
     &               /masse(ig,k)
            pdvadj(ig,k)=((detr_dv2(ig,k)+gamma(ig,k))*zva(ig,k)  &
     &               -(entr(ig,k)+gamma(ig,k))*ve(ig,k)  &
     &               -wvd(ig,k)+wvd(ig,k+1))  &
     &               /masse(ig,k)
         enddo
      enddo


! Sorties eventuelles.
!----------------------

!      if (nlarga0>0) then
!          print*,'WARNING !!!!!! DANS THERMCELL_DV2 '
!      print*,nlarga0,' points pour lesquels laraga=0. dans un thermique'
!          print*,'Il faudrait decortiquer ces points'
!      endif

! ===========================================================================
! ============= FIN DV2 =====================================================
! ===========================================================================

      else
!      detrmod(:,:)=0.
!      do k=1,zlmax
!         do ig=1,ngridmx
!            detrmod(ig,k)=fm(ig,k)-fm(ig,k+1) &
!     &      +entr(ig,k)
!            if (detrmod(ig,k).lt.0.) then
!               entr(ig,k)=entr(ig,k)-detrmod(ig,k)
!               detrmod(ig,k)=0.
!            endif
!         enddo
!      enddo
!
!
!      call thermcell_dqup(ngridmx,nlayermx,ptimestep                &
!     &      ,fm,entr,detrmod,  &
!     &     masse,zu,pduadj,ndt,zlmax)
!
!      call thermcell_dqup(ngridmx,nlayermx,ptimestep    &
!     &       ,fm,entr,detrmod,  &
!     &     masse,zv,pdvadj,ndt,zlmax)

      endif

!------------------------------------------------------------------
!  calcul du transport vertical de traceurs
!------------------------------------------------------------------

! We only transport co2 tracer because it is coupled to the scheme through theta_m
! The rest is transported outside the sub-timestep loop

      ratiom(:,:)=1.

      if (ico2.ne.0) then
      detrmod(:,:)=0.
      do k=1,zlmax
         do ig=1,ngridmx
            detrmod(ig,k)=fm(ig,k)-fm(ig,k+1) &
     &      +entr(ig,k)
            if (detrmod(ig,k).lt.0.) then
               entr(ig,k)=entr(ig,k)-detrmod(ig,k)
               detrmod(ig,k)=0.
            endif
         enddo
      enddo

      call thermcell_dqup(ngridmx,nlayermx,ptimestep     &
     &     ,fm,entr,detrmod,  &
     &    masse,pq(:,:,ico2),pdqadj(:,:,ico2),ndt,zlmax)

! Compute the ratio between theta and theta_m
      
       do l=1,zlmax
          do ig=1,ngridmx
             ratiom(ig,l)=1./(A*(pq(ig,l,ico2)+pdqadj(ig,l,ico2)*ptimestep)+B)
          enddo
       enddo

       endif

!------------------------------------------------------------------
!  incrementation dt
!------------------------------------------------------------------

      pdtadj(:,:)=0.
      do l=1,zlmax
         do ig=1,ngridmx
         pdtadj(ig,l)=(zdthladj(ig,l)+zdthladj_down(ig,l))*zpopsk(ig,l)*ratiom(ig,l)
!         pdtadj(ig,l)=zdthladj(ig,l)*zpopsk(ig,l)*ratiom(ig,l)
         enddo
      enddo

!------------------------------------------------------------------
!  calcul du transport vertical de la tke
!------------------------------------------------------------------

!      modname='tke'
!      call thermcell_dqupdown(ngridmx,nlayermx,ptimestep,fm,entr,detr,  &
!      &      masse,pq2,pdq2adj,ztvd,fm_down,ztv,modname,lmax)

! ===========================================================================
! ============= FIN TRANSPORT ===============================================
! ===========================================================================


!------------------------------------------------------------------
!   Calculs de diagnostiques pour les sorties
!------------------------------------------------------------------
! DIAGNOSTIQUE
! We compute interface values for teta env and th. The last interface
! value does not matter, as the mass flux is 0 there.

    
      do l=1,nlayermx-1
       do ig=1,ngridmx
        teta_th_int(ig,l)=0.5*(ztva(ig,l+1)+ztva(ig,l))*ratiom(ig,l)
        teta_down_int(ig,l) = 0.5*(ztvd(ig,l+1)+ztvd(ig,l))*ratiom(ig,l)
        teta_env_int(ig,l)=0.5*(ztv(ig,l+1)+ztv(ig,l))*ratiom(ig,l)
       enddo
      enddo
      do ig=1,ngridmx
       teta_th_int(ig,nlayermx)=teta_th_int(ig,nlayermx-1)
       teta_env_int(ig,nlayermx)=teta_env_int(ig,nlayermx-1)
       teta_down_int(ig,nlayermx)=teta_down_int(ig,nlayermx-1)
      enddo
        heatFlux(:,:)=0.
        buoyancyOut(:,:)=0.
        buoyancyEst(:,:)=0.
        heatFlux_down(:,:)=0.
      do l=1,zlmax
       do ig=1,ngridmx
        heatFlux(ig,l)=fm(ig,l)*(teta_th_int(ig,l)-teta_env_int(ig,l))/(rhobarz(ig,l))
        buoyancyOut(ig,l)=g*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
        buoyancyEst(ig,l)=g*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)
        heatFlux_down(ig,l)=fm_down(ig,l)*(teta_down_int(ig,l)-teta_env_int(ig,l))/rhobarz(ig,l)
       enddo
      enddo

      return
      end
