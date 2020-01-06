c------------------------------------------------------
c  Fill some parameters in the 'include' files
c  >> Do part of the job previously done by phyetat0.F
c  >> Complete list of parameters is found in tabfi.F
c------------------------------------------------------
c
c Values are defined in the module_model_constants.F WRF routine
c      
      ! in 'comcstfi.h'
      omeg=womeg                
      mugaz=wmugaz  
      print*,"check: rad,cpp,g,r,rcp,daysec,omeg,mugaz"
      print*,rad,cpp,g,r,rcp,daysec,omeg,mugaz
    
      ! in 'planet.h' 
      year_day=wyear_day
      periheli=wperiheli
      aphelie=waphelie
      peri_day=wperi_day
      obliquit=wobliquit
      emin_turb=wemin_turb
      lmixmin=wlmixmin
      print*,"check: year_day,periheli,aphelie,peri_day,obliquit"
      print*,year_day,periheli,aphelie,peri_day,obliquit
      print*,"check: emin_turb,lmixmin"
      print*,emin_turb,lmixmin

      ! in 'surfdat.h'
      emissiv=wemissiv
      emisice(1)=wemissiceN
      emisice(2)=wemissiceS
      albedice(1)=walbediceN
      albedice(2)=walbediceS
      iceradius(1)=wiceradiusN
      iceradius(2)=wiceradiusS
      dtemisice(1)=wdtemisiceN
      dtemisice(2)=wdtemisiceS
      z0_default=wz0
      print*,"check: z0def,emissiv,emisice,albedice,iceradius,dtemisice"
      print*,z0_default,emissiv,emisice,albedice,iceradius,dtemisice

c
c Values are defined in the WPS processing
c  
        albedodat(:)=walbedodat(:)
        !!!!! ***WRF inertiedat was moved, new physics !!
        !inertiedat(:)=winertiedat(:)
        phisfi(:)=wphisfi(:)
        print*,"check: albedodat(1),phisfi(1)"
        print*,albedodat(1),phisfi(1)
        print*,"check: albedodat(end),phisfi(end)"
        print*,albedodat(ngrid),phisfi(ngrid)

        ! NB: usually, gravity wave scheme is useless in mesoscale modeling
        ! NB: we however keep the option for coarse grid case ... 	
        zmea(:)=wzmea(:)
        zstd(:)=wzstd(:)
        zsig(:)=wzsig(:)
        zgam(:)=wzgam(:)
        zthe(:)=wzthe(:)
        print*,"check: gw param"
        print*,zmea(1),zmea(ngrid)
        print*,zstd(1),zstd(ngrid)
        print*,zsig(1),zsig(ngrid)
        print*,zgam(1),zgam(ngrid)
        print*,zthe(1),zthe(ngrid)

        !
        ! in slope.h
        !
        theta_sl(:)=wtheta(:)
        psi_sl(:)=wpsi(:)
        print*,"check: theta_sl(1),psi_sl(1)"
        print*,theta_sl(1),psi_sl(1)
        print*,"check: theta_sl(end),psi_sl(end)"
        print*,theta_sl(ngrid),psi_sl(ngrid)

        !
        ! in comsoil.h
        !
        volcapa=wvolcapa
        print*,"check: volcapa"
        print*,volcapa 

        !
        ! in surfdat.h
        !
        z0(:) = wz0tab(:)
        print*,"check: z0"
        print*,z0(1),z0(ngrid)

c*****************************************************
c Since it comes from WRF settings, we have to
c fill dtphys in the include file 
c It must be set now, because it is used afterwards
c*****************************************************
        dtphys=wdt*ptimestep 
        print*,'Physical timestep (s) ',dtphys
