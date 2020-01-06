! ----------------------------------------------------------------------
! -----------------------------------------------------------------------
!   INCLUDE 'dimradmars.h'

!   Declaration and initialisation or radiative transfer calculations
!------------------------------------------------------------------------
!------------------------------------------------------------------------

! nflev: number of vertical layer
! ndlon,ndlo2: number of horizontal points
! Splitting of horizontal grid
! NDLO2 et ndomainsz pour le decoupage de l'appel a la physique
! ATTENTION:  Il faut  1 < ndomainsz =< ngridmx

      INTEGER  NFLEV,NDLON,NDLO2,ndomainsz

!      parameter (ndomainsz=ngridmx)
      parameter (ndomainsz=(ngridmx-1)/20 + 1)
!      parameter (ndomainsz=(ngridmx-1)/5 + 1) 

      parameter (NFLEV=nlayermx,NDLON=ndomainsz) ! avec decoupage
      parameter (NDLO2=NDLON)

! Number of kind of tracer radiative properties
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! naerkind is set in scatterers.h (built when compiling with makegcm -s #)

#include"scatterers.h"
! NB: May have to change value of nsizemax below when changing scatterers

! Reference wavelengths used to compute reference optical depth (m)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      REAL longrefir(naerkind),longrefvis(naerkind)

      REAL long1vis,long2vis,long3vis, long1ir,long2ir
      REAL long1co2,long2co2
      REAL sunfr(2)
      integer nir, nuco2
      INTEGER npademx,nabsmx,nt_pademx, NSUN


! Definition of spectral intervals at thermal infrared wavelengths (LW)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nir=4) ! Total number of thermal IR bands
      parameter (nuco2=2) ! number of bands in CO2 bands
      PARAMETER (long1ir=5.E-6 , long2ir=200.E-6)
      PARAMETER (long1co2=1.E+0 / 865.E+2 , long2co2=1.E+0 / 500.E+2)

!  Warning : the "nir" thermal IR bands are not ordered by wavelength:
!      iir=1 : central 15um CO2 bands     
!      iir=2 : CO2 band wings    [long1co2-long2co2] MINUS central band
!      iir=3 : 9 um band [long1ir - long1co2]
!      iir=4 : Far IR    [long2co2 - long2ir]
    
!  Definition of spectral interval at solar wavelengths (SW)
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      PARAMETER (NSUN=2)   ! do not change that !
!  Boundaries of spectral intervals (m) : 
      PARAMETER (long1vis=0.1E-6 , long2vis=0.5E-6 , long3vis=5.E-6)
!  Fraction of solar energy in solar band #1 [long1vis-long2vis]
      DATA sunfr(1) / 0.274490 /  
!  Fraction of solar energy in solar band #2 [long2vis-long3vis]
      DATA sunfr(2) / 0.725509 /

! Maximum number of grain size classes
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This parameter has to be set to the maximum number of particle
!   sizes contained in the optical parameter database; For example,
!   if only one grain size is used to describe dust, and 30 are used
!   to describe water-ice crystals in the visible and 15 in the IR,
!   nsizemax has to be set to 30.
! If only one grain size is considered for all the aerosols, set
!   this parameter to 1 and convolution will be turned off during
!   the radiative calculations.

      INTEGER, PARAMETER :: nsizemax = 60
!      INTEGER, PARAMETER :: nsizemax = 1

! Various initialisation for LW radiative code
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! npademx : nombre de coef de pade 
! nabsmx : ?
! nt_pademx : nombre d'intervalles de temperature pour pade

      PARAMETER (npademx=4,nabsmx=2,nt_pademx=19)
