!
! AC 2011-01-05
!
      SUBROUTINE calltherm_interface (firstcall, & 
     & zzlev,zzlay, &
     & ptimestep,pu,pv,pt,pq,pdu,pdv,pdt,pdq,q2, &
     & pplay,pplev,pphi,zpopsk, &
     & pdu_th,pdv_th,pdt_th,pdq_th,lmax,zmaxth,pbl_dtke, &
     & pdhdif,hfmax,wstar,sensibFlux)

       USE ioipsl_getincom

      implicit none
#include "callkeys.h"
#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "tracer.h"

!--------------------------------------------------------
! Input Variables
!--------------------------------------------------------

!      REAL, INTENT(IN) :: long(ngridmx),lati(ngridmx)
      REAL, INTENT(IN) :: ptimestep
      REAL, INTENT(IN) :: pplev(ngridmx,nlayermx+1)
      REAL, INTENT(IN) :: pplay(ngridmx,nlayermx)
      REAL, INTENT(IN) :: pphi(ngridmx,nlayermx)
      REAL, INTENT(IN) :: pu(ngridmx,nlayermx),pv(ngridmx,nlayermx)
      REAL, INTENT(IN) :: pt(ngridmx,nlayermx),pq(ngridmx,nlayermx,nqmx)
      REAL, INTENT(IN) :: zzlay(ngridmx,nlayermx)
      REAL, INTENT(IN) :: zzlev(ngridmx,nlayermx+1) 
      LOGICAL, INTENT(IN) :: firstcall
      REAL, INTENT(IN) :: pdu(ngridmx,nlayermx),pdv(ngridmx,nlayermx)
      REAL, INTENT(IN) :: pdq(ngridmx,nlayermx,nqmx)
      REAL, INTENT(IN) :: pdt(ngridmx,nlayermx)
      REAL, INTENT(IN) :: q2(ngridmx,nlayermx+1)
      REAL, INTENT(IN) :: zpopsk(ngridmx,nlayermx)
      REAL, INTENT(IN) :: pdhdif(ngridmx,nlayermx)
      REAL, INTENT(IN) :: sensibFlux(ngridmx)

!--------------------------------------------------------
! Output Variables
!--------------------------------------------------------

      REAL, INTENT(OUT) :: pdu_th(ngridmx,nlayermx)
      REAL, INTENT(OUT) :: pdv_th(ngridmx,nlayermx)
      REAL, INTENT(OUT) :: pdt_th(ngridmx,nlayermx)
      REAL, INTENT(OUT) :: pdq_th(ngridmx,nlayermx,nqmx)
      INTEGER, INTENT(OUT) :: lmax(ngridmx)
      REAL, INTENT(OUT) :: zmaxth(ngridmx)
      REAL, INTENT(OUT) :: pbl_dtke(ngridmx,nlayermx+1)
      REAL, INTENT(OUT) :: wstar(ngridmx)

!--------------------------------------------------------
! Thermals local variables
!--------------------------------------------------------
      REAL zu(ngridmx,nlayermx), zv(ngridmx,nlayermx)
      REAL zt(ngridmx,nlayermx)
      REAL d_t_ajs(ngridmx,nlayermx)
      REAL d_u_ajs(ngridmx,nlayermx), d_q_ajs(ngridmx,nlayermx,nqmx)
      REAL d_v_ajs(ngridmx,nlayermx) 
      REAL fm_therm(ngridmx,nlayermx+1), entr_therm(ngridmx,nlayermx)
      REAL detr_therm(ngridmx,nlayermx),detrmod(ngridmx,nlayermx)
      REAL zw2(ngridmx,nlayermx+1)
      REAL fraca(ngridmx,nlayermx+1),zfraca(ngridmx,nlayermx+1)
      REAL ztla(ngridmx,nlayermx)
      REAL q_therm(ngridmx,nlayermx), pq_therm(ngridmx,nlayermx,nqmx)
      REAL q2_therm(ngridmx,nlayermx), dq2_therm(ngridmx,nlayermx)
      REAL lmax_real(ngridmx)
      REAL masse(ngridmx,nlayermx)
      LOGICAL qtransport_thermals,dtke_thermals
      INTEGER l,ig,iq,ii(1),k
      CHARACTER (LEN=20) modname

!--------------------------------------------------------
! Local variables for sub-timestep
!--------------------------------------------------------

      REAL d_t_the(ngridmx,nlayermx), d_q_the(ngridmx,nlayermx,nqmx)
      REAL d_u_the(ngridmx,nlayermx),d_v_the(ngridmx,nlayermx)
      REAL dq2_the(ngridmx,nlayermx)
      INTEGER isplit
      INTEGER,SAVE :: nsplit_thermals
      REAL, SAVE :: r_aspect_thermals
      REAL fact
      REAL zfm_therm(ngridmx,nlayermx+1),zdt
      REAL zentr_therm(ngridmx,nlayermx),zdetr_therm(ngridmx,nlayermx)
      REAL zheatFlux(ngridmx,nlayermx)
      REAL zheatFlux_down(ngridmx,nlayermx)
      REAL zbuoyancyOut(ngridmx,nlayermx)
      REAL zbuoyancyEst(ngridmx,nlayermx)
      REAL zzw2(ngridmx,nlayermx+1)
      REAL zmax(ngridmx)
      INTEGER ndt,zlmax

!--------------------------------------------------------
! Diagnostics
!--------------------------------------------------------

      REAL heatFlux(ngridmx,nlayermx)
      REAL heatFlux_down(ngridmx,nlayermx)
      REAL buoyancyOut(ngridmx,nlayermx)
      REAL buoyancyEst(ngridmx,nlayermx)
      REAL hfmax(ngridmx),wmax(ngridmx)
      REAL pbl_teta(ngridmx),dteta(ngridmx,nlayermx)
      REAL rpdhd(ngridmx,nlayermx)
      REAL wtdif(ngridmx,nlayermx),rho(ngridmx,nlayermx)
      REAL wtth(ngridmx,nlayermx)

!--------------------------------------------------------
! Theta_m
!--------------------------------------------------------

      INTEGER ico2
      SAVE ico2

! **********************************************************************
! Initialization
! **********************************************************************

      lmax(:)=0
      pdu_th(:,:)=0.
      pdv_th(:,:)=0.
      pdt_th(:,:)=0.
      entr_therm(:,:)=0.
      detr_therm(:,:)=0.
      q2_therm(:,:)=0.
      dq2_therm(:,:)=0.
      ztla(:,:)=0.
      pbl_dtke(:,:)=0.
      fm_therm(:,:)=0.
      zw2(:,:)=0.
      fraca(:,:)=0.
      zfraca(:,:)=0.
      if (tracer) then
         pdq_th(:,:,:)=0.
      end if
      d_t_ajs(:,:)=0.
      d_u_ajs(:,:)=0.
      d_v_ajs(:,:)=0.
      d_q_ajs(:,:,:)=0.
      heatFlux(:,:)=0.
      heatFlux_down(:,:)=0.
      buoyancyOut(:,:)=0.
      buoyancyEst(:,:)=0.
      zmaxth(:)=0.
      lmax_real(:)=0.


! **********************************************************************
! Preparing inputs for the thermals
! **********************************************************************

       zu(:,:)=pu(:,:)+pdu(:,:)*ptimestep
       zv(:,:)=pv(:,:)+pdv(:,:)*ptimestep
       zt(:,:)=pt(:,:)+pdt(:,:)*ptimestep

       pq_therm(:,:,:)=0.
       qtransport_thermals=.true. !! default setting
       !call getin("qtransport_thermals",qtransport_thermals)

       if(qtransport_thermals) then
          if(tracer) then
                pq_therm(:,:,:)=pq(:,:,:)+pdq(:,:,:)*ptimestep
          endif
       endif

       dtke_thermals=.false. !! default setting
       call getin("dtke_thermals",dtke_thermals)
       IF(dtke_thermals) THEN
          DO l=1,nlayermx
              q2_therm(:,l)=0.5*(q2(:,l)+q2(:,l+1))
          ENDDO
       ENDIF

! **********************************************************************
! Polar night mixing : theta_m
! **********************************************************************

      if(firstcall) then
        ico2=0
        if (tracer) then
!     Prepare Special treatment if one of the tracers is CO2 gas
           do iq=1,nqmx
             if (noms(iq).eq."co2") then
                ico2=iq
             end if
           enddo
        endif
      endif !of if firstcall


! **********************************************************************
! **********************************************************************
! **********************************************************************
! CALLTHERM
! **********************************************************************
! **********************************************************************
! **********************************************************************

!         r_aspect_thermals     ! Mainly control the shape of the temperature profile
                                ! in the surface layer. Decreasing it goes toward
                                ! a convective-adjustment like profile.
!         nsplit_thermals       ! Sub-timestep for the thermals. Very dependant on the
                                ! chosen timestep for the radiative transfer.
                                ! It is recommended to run with 96 timestep per day and 
                                ! iradia = 1., configuration in which thermals can run
                                ! very well with a sub-timestep of 10.
         IF (firstcall) THEN
            r_aspect_thermals=1.  ! same value is OK for GCM and mesoscale
#ifdef MESOSCALE
            !! valid for timesteps < 200s
            nsplit_thermals=4
#else
            IF ((ptimestep .le. 3699.*24./96.) .and. (iradia .eq. 1)) THEN
               nsplit_thermals=10
            ELSE
               nsplit_thermals=35
            ENDIF
#endif
            call getin("nsplit_thermals",nsplit_thermals)
            call getin("r_aspect_thermals",r_aspect_thermals)
         ENDIF

! **********************************************************************
! SUB-TIMESTEP LOOP
! **********************************************************************

         zdt=ptimestep/REAL(nsplit_thermals)

         DO isplit=1,nsplit_thermals

! Initialization of intermediary variables

!         zfm_therm(:,:)=0. !init is done inside
!         zentr_therm(:,:)=0.
!         zdetr_therm(:,:)=0.
!         zheatFlux(:,:)=0.  
!         zheatFlux_down(:,:)=0.
!         zbuoyancyOut(:,:)=0.
!         zbuoyancyEst(:,:)=0.
         zzw2(:,:)=0.
         zmax(:)=0.
         lmax(:)=0
!         d_t_the(:,:)=0. !init is done inside

!         d_u_the(:,:)=0. !transported outside
!         d_v_the(:,:)=0.
         dq2_the(:,:)=0.

         if (nqmx .ne. 0 .and. ico2 .ne. 0) then
            d_q_the(:,:,ico2)=0.
         endif

             CALL thermcell_main_mars(zdt  &
     &      ,pplay,pplev,pphi,zzlev,zzlay  &
     &      ,zu,zv,zt,pq_therm,q2_therm  &
     &      ,d_u_the,d_v_the,d_t_the,d_q_the,dq2_the  &
     &      ,zfm_therm,zentr_therm,zdetr_therm,lmax,zmax  &
     &      ,r_aspect_thermals &
     &      ,zzw2,fraca,zpopsk &
     &      ,ztla,zheatFlux,zheatFlux_down &
     &      ,zbuoyancyOut,zbuoyancyEst)

      fact=1./REAL(nsplit_thermals)

            d_t_the(:,:)=d_t_the(:,:)*ptimestep*fact
!            d_u_the(:,:)=d_u_the(:,:)*ptimestep*fact
!            d_v_the(:,:)=d_v_the(:,:)*ptimestep*fact
            dq2_the(:,:)=dq2_the(:,:)*fact
            if (ico2 .ne. 0) then
               d_q_the(:,:,ico2)=d_q_the(:,:,ico2)*ptimestep*fact
            endif

            zmaxth(:)=zmaxth(:)+zmax(:)*fact
            lmax_real(:)=lmax_real(:)+float(lmax(:))*fact
            fm_therm(:,:)=fm_therm(:,:)  &
     &      +zfm_therm(:,:)*fact
            entr_therm(:,:)=entr_therm(:,:)  &
     &       +zentr_therm(:,:)*fact
            detr_therm(:,:)=detr_therm(:,:)  &
     &       +zdetr_therm(:,:)*fact
            zfraca(:,:)=zfraca(:,:) + fraca(:,:)*fact

            heatFlux(:,:)=heatFlux(:,:) &
     &       +zheatFlux(:,:)*fact
            heatFlux_down(:,:)=heatFlux_down(:,:) &
     &       +zheatFlux_down(:,:)*fact
            buoyancyOut(:,:)=buoyancyOut(:,:) &
     &       +zbuoyancyOut(:,:)*fact
            buoyancyEst(:,:)=buoyancyEst(:,:) &
     &       +zbuoyancyEst(:,:)*fact
  

            zw2(:,:)=zw2(:,:) + zzw2(:,:)*fact

!  accumulation de la tendance

           d_t_ajs(:,:)=d_t_ajs(:,:)+d_t_the(:,:)
!           d_u_ajs(:,:)=d_u_ajs(:,:)+d_u_the(:,:)
!           d_v_ajs(:,:)=d_v_ajs(:,:)+d_v_the(:,:)
            if (ico2 .ne. 0) then
               d_q_ajs(:,:,ico2)=d_q_ajs(:,:,ico2)+d_q_the(:,:,ico2)
            endif
!            dq2_therm(:,:)=dq2_therm(:,:)+dq2_the(:,:)
!  incrementation des variables meteo

            zt(:,:) = zt(:,:) + d_t_the(:,:)
!            zu(:,:) = zu(:,:) + d_u_the(:,:)
!            zv(:,:) = zv(:,:) + d_v_the(:,:)
            if (ico2 .ne. 0) then
             pq_therm(:,:,ico2) = &
     &          pq_therm(:,:,ico2) + d_q_the(:,:,ico2)
            endif
!            q2_therm(:,:) = q2_therm(:,:) + dq2_therm(:,:)


         ENDDO ! isplit
!****************************************************************

      lmax(:)=nint(lmax_real(:))
      zlmax=MAXVAL(lmax(:))+2
      if (zlmax .ge. nlayermx) then
        print*,'thermals have reached last layer of the model'
        print*,'this is not good !'
      endif


! Now that we have computed total entrainment and detrainment, we can
! advect u, v, and q in thermals. (theta already advected). We can do
! that separatly because u,v,and q are not used in thermcell_main for
! any thermals-related computation : they are purely passive.

! mass of cells
      do l=1,nlayermx
         masse(:,l)=(pplev(:,l)-pplev(:,l+1))/g
      enddo

      detrmod(:,:)=0.
      do l=1,zlmax
         do ig=1,ngridmx
            detrmod(ig,l)=fm_therm(ig,l)-fm_therm(ig,l+1) &
     &      +entr_therm(ig,l)
            if (detrmod(ig,l).lt.0.) then
               entr_therm(ig,l)=entr_therm(ig,l)-detrmod(ig,l)
               detrmod(ig,l)=0.
            endif
         enddo
      enddo
      ndt=10
      call thermcell_dqup(ngridmx,nlayermx,ptimestep                &
     &      ,fm_therm,entr_therm,detrmod,  &
     &     masse,zu,d_u_ajs,ndt,zlmax)

      call thermcell_dqup(ngridmx,nlayermx,ptimestep    &
     &       ,fm_therm,entr_therm,detrmod,  &
     &     masse,zv,d_v_ajs,ndt,zlmax)

      if (nqmx .ne. 0.) then
      DO iq=1,nqmx
      if (iq .ne. ico2) then
      call thermcell_dqup(ngridmx,nlayermx,ptimestep     &
     &     ,fm_therm,entr_therm,detrmod,  &
     &    masse,pq_therm(:,:,iq),d_q_ajs(:,:,iq),ndt,zlmax)
      endif
      ENDDO
      endif

      if (dtke_thermals) then
      detrmod(:,:)=0.
      ndt=10
      do l=1,zlmax
         do ig=1,ngridmx
            detrmod(ig,l)=fm_therm(ig,l)-fm_therm(ig,l+1) &
     &      +entr_therm(ig,l)
            if (detrmod(ig,l).lt.0.) then
               entr_therm(ig,l)=entr_therm(ig,l)-detrmod(ig,l)
               detrmod(ig,l)=0.
            endif
         enddo
      enddo
      call thermcell_dqup(ngridmx,nlayermx,ptimestep     &
     &     ,fm_therm,entr_therm,detrmod,  &
     &    masse,q2_therm,dq2_therm,ndt,zlmax)
      endif

      DO ig=1,ngridmx
         wmax(ig)=MAXVAL(zw2(ig,:))
      ENDDO

! **********************************************************************
! **********************************************************************
! **********************************************************************
! CALLTHERM END
! **********************************************************************
! **********************************************************************
! **********************************************************************


! **********************************************************************
! Preparing outputs
! **********************************************************************

      do l=1,zlmax
        pdu_th(:,l)=d_u_ajs(:,l)
        pdv_th(:,l)=d_v_ajs(:,l)
      enddo

           if(qtransport_thermals) then
              if(tracer) then
               do iq=1,nqmx
                if (iq .ne. ico2) then
                  do l=1,zlmax
                     pdq_th(:,l,iq)=d_q_ajs(:,l,iq)
                  enddo
                else
                  do l=1,zlmax
                     pdq_th(:,l,iq)=d_q_ajs(:,l,iq)/ptimestep
                  enddo
                endif
               enddo
              endif
           endif

           IF(dtke_thermals) THEN
              DO l=2,nlayermx
                 pbl_dtke(:,l)=0.5*(dq2_therm(:,l-1)+dq2_therm(:,l))
              ENDDO
  
              pbl_dtke(:,1)=0.5*dq2_therm(:,1)
              pbl_dtke(:,nlayermx+1)=0.
           ENDIF

           do l=1,zlmax
              pdt_th(:,l)=d_t_ajs(:,l)/ptimestep
           enddo


! **********************************************************************
! Compute the free convection velocity scale for vdifc
! **********************************************************************


! Potential temperature gradient

      dteta(:,nlayermx)=0.
      DO l=1,nlayermx-1
         DO ig=1, ngridmx
            dteta(ig,l) = ((zt(ig,l+1)-zt(ig,l))/zpopsk(ig,l))          &
     &              /(zzlay(ig,l+1)-zzlay(ig,l))
         ENDDO
      ENDDO

! Computation of the pbl mixed layer temperature

      DO ig=1, ngridmx
         ii=MINLOC(abs(dteta(ig,1:lmax(ig))))
         pbl_teta(ig) = zt(ig,ii(1))/zpopsk(ig,ii(1))
      ENDDO

! we must add the heat flux from the diffusion scheme to hfmax

! compute rho as it is after the diffusion

      rho(:,:)=pplay(:,:)                                               &
     & /(r*(pt(:,:)+pdhdif(:,:)*zpopsk(:,:)*ptimestep))

! integrate -rho*pdhdif

      rpdhd(:,:)=0.

      DO ig=1,ngridmx
       DO l=1,lmax(ig)
        rpdhd(ig,l)=0.
        DO k=1,l
         rpdhd(ig,l)=rpdhd(ig,l)-rho(ig,k)*pdhdif(ig,k)*                &
     & (zzlev(ig,k+1)-zzlev(ig,k))
        ENDDO
        rpdhd(ig,l)=rpdhd(ig,l)-sensibFlux(ig)/cpp
       ENDDO
      ENDDO

! compute w'teta' from diffusion

      wtdif(:,:)=rpdhd(:,:)/rho(:,:)

! compute rho as it is after the thermals

      rho(:,:)=pplay(:,:)/(r*(zt(:,:)))
! integrate -rho*pdhdif

      DO ig=1,ngridmx
       DO l=1,lmax(ig)
        rpdhd(ig,l)=0.
        DO k=1,l
         rpdhd(ig,l)=rpdhd(ig,l)-rho(ig,k)*(pdt_th(ig,k)/zpopsk(ig,k))* &
     & (zzlev(ig,k+1)-zzlev(ig,k))
        ENDDO
        rpdhd(ig,l)=rpdhd(ig,l)+                                        &
     &    rho(ig,1)*(heatFlux(ig,1)+heatFlux_down(ig,1))
       ENDDO
      ENDDO
      rpdhd(:,nlayermx)=0.

! compute w'teta' from thermals

      wtth(:,:)=rpdhd(:,:)/rho(:,:)

! We get the max heat flux from thermals and add the contribution from the diffusion

      DO ig=1,ngridmx
        hfmax(ig)=MAXVAL(wtth(ig,:)+wtdif(ig,:))
      ENDDO
! We follow Spiga et. al 2010 (QJRMS)
! ------------

      DO ig=1, ngridmx
         IF (zmax(ig) .gt. 0.) THEN
            wstar(ig)=(g*zmaxth(ig)*hfmax(ig)/pbl_teta(ig))**(1./3.)
         ELSE
            wstar(ig)=0.
         ENDIF
      ENDDO



! **********************************************************************
! Diagnostics
! **********************************************************************
        
        if(outptherm) then
        if (ngridmx .eq. 1) then
        call WRITEDIAGFI(ngridmx,'entr_therm','entrainement thermique',&
     &                       'kg/m-2',1,entr_therm)
        call WRITEDIAGFI(ngridmx,'detr_therm','detrainement thermique',&
     &                       'kg/m-2',1,detr_therm)
        call WRITEDIAGFI(ngridmx,'fm_therm','flux masse thermique',&
     &                       'kg/m-2',1,fm_therm)
        call WRITEDIAGFI(ngridmx,'zw2','vitesse verticale thermique',&
     &                       'm/s',1,zw2)
        call WRITEDIAGFI(ngridmx,'heatFlux_up','heatFlux_updraft',&
     &                       'SI',1,heatFlux)
       call WRITEDIAGFI(ngridmx,'heatFlux_down','heatFlux_downdraft',&
     &                       'SI',1,heatFlux_down)
        call WRITEDIAGFI(ngridmx,'fraca','fraction coverage',&
     &                       'percent',1,fraca)
        call WRITEDIAGFI(ngridmx,'buoyancyOut','buoyancyOut',&
     &                       'm.s-2',1,buoyancyOut)
        call WRITEDIAGFI(ngridmx,'buoyancyEst','buoyancyEst',&
     &                       'm.s-2',1,buoyancyEst)
        call WRITEDIAGFI(ngridmx,'d_t_th',  &
     &         'tendance temp TH','K',1,d_t_ajs)
        call WRITEDIAGFI(ngridmx,'d_q_th',  &
     &         'tendance traceur TH','kg/kg',1,d_q_ajs)
        call WRITEDIAGFI(ngridmx,'zmax',  &
     &         'pbl height','m',0,zmaxth)
        call WRITEDIAGFI(ngridmx,'d_u_th',  &
     &         'tendance moment','m/s',1,pdu_th)
        call WRITEDIAGFI(ngridmx,'wtdif',  &
     &         'heat flux from diffusion','K.m/s',1,wtdif)
        call WRITEDIAGFI(ngridmx,'wtth',  &
     &         'heat flux from thermals','K.m/s',1,wtth)
        call WRITEDIAGFI(ngridmx,'wttot',  &
     &         'heat flux PBL','K.m/s',1,wtdif(:,:)+wtth(:,:))

      else

        call WRITEDIAGFI(ngridmx,'entr_therm','entrainement thermique',&
     &                       'kg/m-2',3,entr_therm)
        call WRITEDIAGFI(ngridmx,'detr_therm','detrainement thermique',&
     &                       'kg/m-2',3,detr_therm)
        call WRITEDIAGFI(ngridmx,'fm_therm','flux masse thermique',&
     &                       'kg/m-2',3,fm_therm)
        call WRITEDIAGFI(ngridmx,'zw2','vitesse verticale thermique',&
     &                       'm/s',3,zw2)
        call WRITEDIAGFI(ngridmx,'heatFlux','heatFlux',&
     &                       'SI',3,heatFlux)
        call WRITEDIAGFI(ngridmx,'buoyancyOut','buoyancyOut',&
     &                       'SI',3,buoyancyOut)
        call WRITEDIAGFI(ngridmx,'d_t_th',  &
     &         'tendance temp TH','K',3,d_t_ajs)

      endif
      endif

       END
