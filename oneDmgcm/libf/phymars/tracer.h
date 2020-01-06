!-----------------------------------------------------------------------
! INCLUDE 'tracer.h'

      character*20  noms(nqmx)  ! name of the tracer
      real mmol(nqmx)           ! mole mass of tracer (g/mol-1) 
      real radius(nqmx)   ! dust and ice particle radius (m)
      real rho_q(nqmx)    ! tracer densities (kg.m-3)
      real alpha_lift(nqmx) ! saltation vertical flux/horiz flux ratio (m-1)
      real alpha_devil(nqmx) ! lifting coeeficient by dust devil

      real varian      ! Characteristic variance of log-normal distribution
      real r3n_q     ! used to compute r0 from number and mass mixing ratio
      real rho_dust     ! Mars dust density (kg.m-3)
      real rho_ice     ! Water ice density (kg.m-3)
      real nuice_ref   ! Effective variance of the water ice dist.
      real nuice_sed   ! Sedimentation effective variance of the water ice dist.
      real ref_r0        ! for computing reff=ref_r0*r0 (in log.n. distribution)
      
      real ccn_factor  ! ratio of nuclei for water ice particles

      real dryness(ngridmx)!"Dryness coefficient" for grnd water ice sublimation
      
! tracer indexes: these are initialized in initracer and should be 0 if the
!                 corresponding tracer does not exist
      ! dust
      integer :: igcm_dustbin(nqmx) ! for dustbin 'dust' tracers
      ! dust, special doubleq case
      integer :: igcm_dust_mass   ! dust mass mixing ratio
                                  !   (for transported dust)
      integer :: igcm_dust_number ! dust number mixing ratio
                                  !   (transported dust)
      integer :: igcm_ccn_mass   ! CCN mass mixing ratio
      integer :: igcm_ccn_number ! CCN number mixing ratio
      integer :: igcm_dust_submicron ! submicron dust mixing ratio
                                     !   (transported dust)
      ! water
      integer :: igcm_h2o_vap ! water vapour
      integer :: igcm_h2o_ice ! water ice
      ! chemistry:
      integer :: igcm_co2
      integer :: igcm_co
      integer :: igcm_o
      integer :: igcm_o1d
      integer :: igcm_o2
      integer :: igcm_o3
      integer :: igcm_h
      integer :: igcm_h2
      integer :: igcm_oh
      integer :: igcm_ho2
      integer :: igcm_h2o2
      integer :: igcm_n2
      integer :: igcm_ar
      integer :: igcm_n
      integer :: igcm_no
      integer :: igcm_no2
      integer :: igcm_n2d
      integer :: igcm_ch4
      ! Ions
      integer :: igcm_co2plus
      integer :: igcm_oplus
      integer :: igcm_o2plus
      integer :: igcm_coplus
      integer :: igcm_cplus
      integer :: igcm_nplus
      integer :: igcm_noplus
      integer :: igcm_n2plus 
      integer :: igcm_hplus
      integer :: igcm_hco2plus
      integer :: igcm_elec
      ! other tracers
      integer :: igcm_ar_n2 ! for simulations using co2 +neutral gas


! NB: to keep commons aligned: 
!     split them in groups (reals, integers and characters)
      COMMON/tracer/radius,rho_q,alpha_lift,alpha_devil,mmol,           &
     & varian,r3n_q,rho_dust,rho_ice,nuice_ref,nuice_sed,               &
     & ref_r0,ccn_factor,dryness
      COMMON/tracer2/                                                   &
     & igcm_dustbin,igcm_dust_mass,igcm_dust_number,                    &
     & igcm_ccn_mass,igcm_ccn_number,igcm_dust_submicron,               &
     & igcm_h2o_vap,igcm_h2o_ice,igcm_co2,igcm_co,igcm_o,igcm_o1d,      &
     & igcm_o2,igcm_o3,igcm_h,igcm_h2,igcm_oh,igcm_ho2,igcm_h2o2,       &
     & igcm_n2,igcm_ar,igcm_n,igcm_no,igcm_no2,igcm_n2d,                &
     & igcm_ch4,                                                        &
     & igcm_co2plus,igcm_oplus,igcm_o2plus,igcm_coplus,igcm_cplus,      &
     & igcm_nplus,igcm_noplus,igcm_n2plus,igcm_hplus,igcm_elec,         &
     & igcm_hco2plus,igcm_ar_n2!,nbqchem,niqchem
      COMMON/tracer3/noms
!-----------------------------------------------------------------------
