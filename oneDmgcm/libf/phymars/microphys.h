!-----------------------------------------------------------------------
! INCLUDE 'microphys.h'
! Parameters and physical constants used by the microphysal scheme;
!-----------------------------------------------------------------------

!     Number of bins
      INTEGER, PARAMETER :: nbin_cld = 5

!     Reference temperature, T=273.15 K
      REAL, PARAMETER :: To = 273.15
!     Avogadro number
      DOUBLE PRECISION, PARAMETER :: nav = 6.023d23
!     Perfect gas constant
      DOUBLE PRECISION, PARAMETER :: rgp = 8.3143
!     Boltzman constant
      DOUBLE PRECISION, PARAMETER :: kbz = 1.381d-23
!     Molecular weight of H2O (kg.mol-1)
      DOUBLE PRECISION, PARAMETER :: mh2o = 18.d-3
!     Molecular weight of CO2 (kg.mol-1)
      DOUBLE PRECISION, PARAMETER :: mco2 = 44.d-3
!     Effective CO2 gas molecular radius (m)
      DOUBLE PRECISION, PARAMETER :: molco2 = 2.2d-10
!     Effective H2O gas molecular radius (m)
      DOUBLE PRECISION, PARAMETER :: molh2o = 1.2d-10
!     Surface tension of ice/vapor (N.m)
      DOUBLE PRECISION, PARAMETER :: sigh2o = 0.12
!     Activation energy for desorption of
!       water on a dust-like substrate
!       (J/molecule)
      DOUBLE PRECISION, PARAMETER :: desorp = 0.288e-19
!     Jump frequency of a water molecule (s-1)
      DOUBLE PRECISION, PARAMETER :: nus = 1.e+13
!     Estimated activation energy for
!       surface diffusion of water molecules
!       (J/molecule)
      DOUBLE PRECISION, PARAMETER :: surfdif = desorp / 10.
!     Weight of a water molecule (kg)
      DOUBLE PRECISION, PARAMETER :: m0 = mh2o / nav

!     Contact parameter ( m=cos(theta) )
!       (initialized in improvedclouds.F)
      REAL mteta

!     Volume of a water molecule (m3)
      DOUBLE PRECISION vo1
!     Radius used by the microphysical scheme (m)
      DOUBLE PRECISION rad_cld(nbin_cld)


! NB: to keep commons aligned: 
!     split them in groups (reals, integers and characters)

      COMMON/microphys/rad_cld,vo1
      
      COMMON/microphys_2/mteta
      
!     EXAMPLE:
!     COMMON/tracer/radius,rho_q,alpha_lift,alpha_devil,mmol,           &
!    & varian,r3n_q,rho_dust,rho_ice,nuice_ref,nuice_sed,               &
!    & ref_r0,dryness
!-----------------------------------------------------------------------
