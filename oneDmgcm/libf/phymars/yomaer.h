!     Radiative characteristics of the aerosols
! 
!   Shortwave
!   ~~~~~~~~~
! 
! tauvis: dust optical depth at reference wavelength  ("longrefvis" set
! in dimradmars.h : typically longrefvis = 0.67E-6 m, as measured by Viking )

! For the "naerkind" kind of aerosol radiative properties : 
! QVISsQREF  :  Qext / Qext("longrefvis")   <--- For both solar bands
! omegavis   :  sinle scattering albedo     <--- For both solar bands
! gvis       :  assymetry factor            <--- For both solar bands
! 
!   Longwave
!   ~~~~~~~~
! 
! For the "naerkind" kind of aerosol radiative properties : 
! QIRsQREF :  Qext / Qext("longrefvis")     <--- For the nir bandes IR
! omegaIR  :  mean single scattering albedo <--- For the nir bandes IR
! gIR      :  mean assymetry factor         <--- For the nir bandes IR
! 
      real tauvis
      real QVISsQREF(nsun,naerkind,nsizemax)
      real omegavis(nsun,naerkind,nsizemax)
      real gvis(nsun,naerkind,nsizemax)
      real QIRsQREF(nir,naerkind,nsizemax)
      real omegaIR(nir,naerkind,nsizemax)
      real gIR(nir,naerkind,nsizemax)

! Actual number of grain size classes in each domain for a
!   given aerosol:

      INTEGER          :: nsize(naerkind,2)

! Particle size axis (depend on the kind of aerosol and the
!   radiation domain)

      REAL :: radiustab(naerkind,2,nsizemax)

! Extinction coefficient at reference wavelengths;
!   These wavelengths are defined in dimradmars.h, and called
!   longrefvis and longrefir.

      REAL             :: QREFvis(naerkind,nsizemax)
      REAL             :: QREFir(naerkind,nsizemax)
      REAL             :: omegaREFvis(naerkind,nsizemax)
      REAL             :: omegaREFir(naerkind,nsizemax)

! For Fortran 77/Fortran 90 compliance always use line continuation
! symbols '&' in columns 73 and 6

      COMMON/YOMAER/tauvis,QVISsQREF,omegavis,gvis,                     &
     &              QIRsQREF,omegaIR,gIR,                               &
     &              QREFvis,QREFir,omegaREFvis,omegaREFir,              &
     &              nsize,radiustab
