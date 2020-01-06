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
      integer :: igcm_13ch4
      integer :: igcm_ch3o2
      integer :: igcm_ch3
      integer :: igcm_ch3oh
      integer :: igcm_hcho
      integer :: igcm_ch3ooh
      integer :: igcm_c2h6
      integer :: igcm_ch2choh
      integer :: igcm_hoch2ch2o2
      integer :: igcm_hoch2ch2o
      integer :: igcm_ethgly
      integer :: igcm_hyetho2h
      integer :: igcm_ch3choho2
      integer :: igcm_ch3chohooh
      integer :: igcm_hcoch2o2
      integer :: igcm_glyox
      integer :: igcm_hcoco
      integer :: igcm_hooch2cho
      integer :: igcm_hoch2cho
      integer :: igcm_hochcho
      integer :: igcm_hoch2co
      integer :: igcm_hoch2co3
      integer :: igcm_hoch2co2h
      integer :: igcm_hcoco2h
      integer :: igcm_hoch2co3h
      integer :: igcm_hcoco3h
      integer :: igcm_hcoco3	  
      integer :: igcm_c2h4
      integer :: igcm_c2h5o2
      integer :: igcm_c2h2
      integer :: igcm_hcooh
      integer :: igcm_1ch2
      integer :: igcm_3ch2
	  integer :: igcm_ch2
      integer :: igcm_ch
	  integer :: igcm_ch3o
	  integer :: igcm_c2h5ooh
      integer :: igcm_ch3cho
      integer :: igcm_ch3cooo
      integer :: igcm_ch3cooh
      integer :: igcm_ch3coooh
	  integer :: igcm_hco
      integer :: igcm_c2h5
      integer :: igcm_c3h8
      integer :: igcm_ic3h7o2
      integer :: igcm_ic3h7ooh
      integer :: igcm_ch3coch3
      integer :: igcm_ch3coch2o2
      integer :: igcm_hoch2o2
      integer :: igcm_hoch2ooh
      integer :: igcm_hoch2oh
      integer :: igcm_c2h5o
      integer :: igcm_c2h5oh 
      integer :: igcm_ch3co 
      integer :: igcm_ch2oo_e
      integer :: igcm_ch2oo	  
	  
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
     & igcm_ch4,igcm_ch3o2,igcm_ch3,igcm_ch3oh,igcm_hcho,igcm_ch3ooh,   &
     & igcm_c2h6,igcm_c2h4,igcm_c2h5o2,igcm_c2h2,igcm_hcooh,igcm_ch2,   &
     & igcm_ch3o,igcm_co2plus,igcm_oplus,igcm_o2plus,igcm_coplus,       &
     & igcm_cplus,igcm_nplus,igcm_noplus,igcm_n2plus,igcm_hplus,        &
     & igcm_elec,igcm_hco2plus,igcm_ar_n2,igcm_c2h5ooh,igcm_ch3cho,     &
     & igcm_ch3cooo,igcm_ch3cooh,igcm_ch3coooh,igcm_hco,igcm_c2h5,      &
     & igcm_c3h8,igcm_ic3h7o2,igcm_ic3h7ooh,igcm_ch3coch3,igcm_13ch4,   &
     & igcm_ch3coch2o2,igcm_hoch2o2,igcm_hoch2ooh,igcm_hoch2oh,         &
     & igcm_c2h5o,igcm_c2h5oh,igcm_ch3co,igcm_ch3choho2,igcm_ch2oo_e,   &
     & igcm_ch2oo,igcm_hoch2ch2o2,igcm_hoch2ch2o,igcm_ethgly,           &
     & igcm_hyetho2h,igcm_ch3chohooh,igcm_hcoch2o2,                     &
     & igcm_glyox,igcm_hcoco,igcm_hooch2cho,igcm_hoch2cho,igcm_hochcho, &
     & igcm_hoch2co,igcm_hoch2co3,igcm_hoch2co2h,igcm_hcoco2h,          &
     & igcm_hoch2co3h,igcm_hcoco3h,igcm_hcoco3,igcm_ch2choh!,nbqchem,niqchem               
      COMMON/tracer3/noms
!-----------------------------------------------------------------------
