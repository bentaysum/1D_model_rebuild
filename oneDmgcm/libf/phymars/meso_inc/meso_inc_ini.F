ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c MESOSCALE MODEL
c 
c       No need to use startfi.nc
c               > part of the job of phyetat0 is done in inifis
c               > remaining initializations are passed here from the WRF variables
c               > beware, some operations were done by phyetat0 (ex: tracers)
c                       > if any problems, look in phyetat0
c
c    LMD_MM_MARS
c       day_ini,tsurf,tsoil,emis,q2,qsurf,co2ice are inputs
c               and locally saved variables
c                       (no need to call phyetat0)
c
      tsurf(:)=wtsurf(:)
      PRINT*,'check: tsurf ',tsurf(1),tsurf(ngridmx)
      tsoil(:,:)=wtsoil(:,:)
      PRINT*,'check: tsoil ',tsoil(1,1),tsoil(ngridmx,nsoilmx)
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     !!!new physics
c      PRINT*,'no variable z0 field for the moment. use z0_default.'
c      z0(:) = z0_default
c      PRINT*, 'check z0 ', z0(1)
      inertiedat(:,:)=wisoil(:,:)
      PRINT*,'check: inert ',inertiedat(1,1),inertiedat(ngridmx,nsoilmx)
      mlayer(0:nsoilmx-1)=wdsoil(1,:)
      PRINT*,'check: midlayer ', mlayer(:)
            !!!!!!!!!!!!!!!!! DONE in soil_setting.F 
            ! 1.5 Build layer(); following the same law as mlayer()
            ! Assuming layer distribution follows mid-layer law:
            ! layer(k)=lay1*alpha**(k-1)
            lay1=sqrt(mlayer(0)*mlayer(1))
            alpha=mlayer(1)/mlayer(0)
            do iloop=1,nsoilmx
              layer(iloop)=lay1*(alpha**(iloop-1))
            enddo

      PRINT*,'check: layer ', layer(:)

            !!!!!!!!!!!!!!!!! DONE in soil_setting.F
      tnom(:)=wtnom(:)   !! est rempli dans advtrac.h
      PRINT*,'check: tracernames ', tnom
     !!!new physics
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      emis(:)=wemis(:)
      PRINT*,'check: emis ',emis(1),emis(ngridmx)
      q2(:,:)=wq2(:,:)
      PRINT*,'check: q2 ',q2(1,1),q2(ngridmx,nlayermx+1)
      qsurf(:,:)=wqsurf(:,:)
      PRINT*,'check: qsurf ',qsurf(1,1),qsurf(ngridmx,nqmx)
      co2ice(:)=wco2ice(:)
      PRINT*,'check: co2 ',co2ice(1),co2ice(ngridmx)
      day_ini=wday_ini

c       artificially filling dyn3d/control.h is also required
c       > iphysiq is put in WRF to be set easily (cf ptimestep)
c       > day_step is simply deduced:
c
      day_step=daysec/ptimestep
      PRINT*,'Call to LMD physics:',day_step,' per Martian day'
c
      iphysiq=ptimestep
c
      !DO iq=1, nq 
      !  PRINT*, tnom(iq), pq(:,:,iq)
      !ENDDO

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !!!! see meso_inc_caps
      !!!! this is a test to change outliers' albedo and thermal inertia
      !alb_lim = 0.26
      !lat_lim = 70.
      !inertie_lim = 800. 
      !PRINT *, 'lat_lim ',lat_lim
      !PRINT *, 'alb_lim ',alb_lim
      !PRINT *, 'inertie_lim ',inertie_lim 
      !!!!
      !!!!
      !IF ( caps .and. water ) THEN
      !    do ig=1,ngridmx
      !      if ( lati(ig)*180./pi .gt. lat_lim ) then
      !         if ( albedodat(ig) .ge. alb_lim ) then
      !              albedodat(ig) = alb_surfice
      !              inertiedat(ig,1) = inertie_lim  
      !         endif
      !         if (inertiedat(ig,1) .ge. inertie_lim ) then
      !              inertiedat(ig,1) = inertie_lim
      !         endif
      !      endif  ! (lati, albedodat)
      !    end do ! (ngridmx)
      !ENDIF ! (caps)
