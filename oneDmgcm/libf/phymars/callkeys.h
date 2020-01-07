!
! For Fortran 77/Fortran 90 compliance always use line continuation
! symbols '&' in columns 73 and 6
!
! NB: to keep commons aligned, it is better to split them in groups
!     of given types (logical, integer, real, ...)

      COMMON/callkeys_l/callrad,calldifv,calladj,callcond,callsoil      &
     &   ,season,diurnal,lwrite,calllott,callstats,calleofdump          &
     &   ,callnirco2,callnlte,callthermos,callconduct,calleuv           &
     &   ,callmolvis,callmoldiff,thermochem,thermoswater,callemis       &
     &   ,callg2d,linear,rayleigh,tracer,active,doubleq,submicron       &
     &   ,lifting,callddevil,scavenging,sedimentation,activice,water    &
     &   ,tifeedback,microphys,caps,photochem,calltherm,outptherm       &
     &   ,callrichsl,callslope,tituscap,long_mean
     
      COMMON/callkeys_i/iradia,iaervar,iddist,ilwd,ilwb,ilwn,ncouche    &
     &   ,dustbin,nltemodel,nircorr,solvarmod,solvaryear
     
      COMMON/callkeys_r/topdustref,solarcondate,semi,alphan,euveff,     &
     &   tke_heat_flux
     
      LOGICAL callrad,calldifv,calladj,callcond,callsoil,               &
     &   season,diurnal,lwrite,calllott                                 &
     &   ,callstats,calleofdump                                         &
     &   ,callnirco2,callnlte,callthermos,callconduct,                  &
     &    calleuv,callmolvis,callmoldiff,thermochem,thermoswater        &
     &   ,calltherm,outptherm,callrichsl,callslope,tituscap


      logical callemis
      logical callg2d
      logical linear

      logical long_mean
	  
      real topdustref
      real semi
      real alphan
      real solarcondate
      real euveff
      real tke_heat_flux

      integer iddist
      integer iaervar
      integer iradia
      integer ilwd
      integer ilwb
      integer ilwn
      integer ncouche
      integer solvarmod   ! model for solar EUV variation
      integer solvaryear  ! mars year for realisticly varying solar EUV 

      logical rayleigh
      logical tracer
      integer dustbin
      logical active,doubleq,submicron,lifting,callddevil,scavenging
      logical sedimentation
      logical water,activice,tifeedback,microphys,caps
      logical photochem
      integer nltemodel
      integer nircorr

      integer swrtype ! type of short wave (solar wavelength) radiative
      ! transfer to use 1: Fouquart 2: Toon.
      parameter (swrtype=2)
!      parameter (swrtype=2)
