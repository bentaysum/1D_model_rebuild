         if (.not.callrichsl) then
       
             DO ig=1,ngrid
!! sensible heat flux in W/m2

             sensheat(ig) = zflubid(ig)-capcal(ig)*zdtsdif(ig)

!! u star in similarity theory in m/s
             ustar(ig) = 0.4
     .               * sqrt( pu(ig,1)*pu(ig,1) + pv(ig,1)*pv(ig,1) )
     .               / log( 1.E+0 + zzlay(ig,1)/z0_default )
             ENDDO

         else

            DO ig=1,ngrid

! New SL parametrization, correct formulation for sensheat :

            sensheat(ig) = (pplay(ig,1)/(r*pt(ig,1)))*cpp
     &        *sqrt(pu(ig,1)*pu(ig,1) + pv(ig,1)*pv(ig,1)
     &        + (log(1.+0.7*wstar(ig) + 2.3*wstar(ig)**2))**2)
     &        *zcdh(ig)*(tsurf(ig)-zh(ig,1))

! New SL parametrization, ustar is more accurately computed in vdif_cd :

            ustar(ig) = sqrt(zcdv(ig)*
     &   (pu(ig,1)*pu(ig,1) + pv(ig,1)*pv(ig,1) 
     &     + (log(1.+0.7*wstar(ig) + 2.3*wstar(ig)**2))**2)
     &                     )
            ENDDO   

         endif  !of if callrichsl

!         write (*,*) 'PHYS HFX cp zdts', sensheat(100), zflubid(100), 
!     .       capcal(100), 
!     .       zdtsdif(100)
!         write (*,*) 'PHYS UST', ustar(100) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! LES LES 
       IF (flag_LES) THEN        

         write (*,*) '************************************************' 
         write (*,*) '** LES mode: the difv part is only used to'
         write (*,*) '**  provide HFX and UST to the dynamics'
         write (*,*) '** NB: - dudif, dvdif, dhdif, dqdif are set to 0'
         write (*,*) '**     - tsurf is updated'     
         write (*,*) '************************************************'

         DO ig=1,ngrid
          wstar(ig)=0.  !! no additional gustiness needed in surface layer (see vdifc.F)
          DO l=1,nlayer
            zdvdif(ig,l) = 0.
            zdudif(ig,l) = 0.
            zdhdif(ig,l) = 0.
          ENDDO
         ENDDO
         IF (lifting .and. doubleq) THEN 
         !! lifted dust is injected in the first layer. 
         !! Sedimentation must be called after turbulent mixing, i.e. on next step, after WRF. 
         !! => lifted dust is not incremented before the sedimentation step.
            zdqdif(:,1,:)=0.
            zdqdif(:,1,igcm_dust_number)=-zdqsdif(:,igcm_dust_number)
            zdqdif(:,1,igcm_dust_mass)=-zdqsdif(:,igcm_dust_mass)
            zdqdif(:,2:nlayermx,:)=0.
            DO iq=1, nq
              IF ((iq .ne. igcm_dust_mass)
     &        .and. (iq .ne. igcm_dust_number)) THEN
                zdqsdif(:,iq)=0.
              ENDIF
            ENDDO
         ELSE
            DO ig=1,ngrid
             DO l=1,nlayer
              DO iq=1, nq
                zdqdif(ig,l,iq) = 0.
                zdqsdif(ig,iq) = 0. !! sortir de la boucle
              ENDDO 
             ENDDO
            ENDDO
         ENDIF
      ENDIF
!!! LES LES        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
