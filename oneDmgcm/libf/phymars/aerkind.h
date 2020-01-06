!------------------------------------------------------------------
! INCLUDE 'aerkind.h'
! Contains the names of the different scatterers
!------------------------------------------------------------------

!     Don't forget to set up the right number of scatterer
!       (naerkind) in dimradmars.h!
      character*20  name_iaer(naerkind)  ! name of the scatterers
!     Scatterer: DUST
      integer iaer_dust_conrath ! Typical dust profiles using a
                                ! Conrath type analytical equation
      integer iaer_dust_doubleq ! Dust profile is given by the
                                ! mass mixing ratio of the two-
                                ! moment scheme method (doubleq)
      integer iaer_dust_submicron ! Dust profile is given by a
                                  ! submicron population of dust
                                  ! particles
!     Scatterer: WATER
      integer iaer_h2o_ice ! Water ice particles

! NB: to keep commons aligned: 
!     split them in groups (reals, integers and characters)
      COMMON/aerkind/                                                   &
     & iaer_dust_conrath,iaer_dust_doubleq,iaer_dust_submicron,         &
     & iaer_h2o_ice
      COMMON/aerkind2/name_iaer
!------------------------------------------------------------------
