SUBROUTINE cl_dust(cc, nesp, dens, rdust, temp, press, &
                   dust001, dust002)

IMPLICIT NONE


#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "conc.h"

! Input 
! -----
REAL cc(nlayermx,nesp) ! Species number densities 
INTEGER nesp           ! Number of tracers 
REAL dens(nlayermx)    ! Atmospheric Number Density [cm-3]
REAL rdust(nlayermx)   ! Mean Geometric radius of dust particles [m]
REAL temp(nlayermx)    ! temperature [K]
REAL press(nlayermx)   ! Pressure [Pa]
REAL dt                ! timestep

! Output 
! ------
REAL dust001(nlayermx) ! Rate of OH + dust -> Cl + Products
REAL dust002(nlayermx) ! Rate of HCl adsorption onto dust 

! Local 
! =====
INTEGER l              ! Layer iterator 
REAL particle_volume   ! Volume of single dust grain 
REAL mass_conc         ! Dust mass concentration 
REAL vol_conc          ! Dust volume concentration
REAL S                 ! Dust Surface Area per Unit volume
REAL pv                ! Water Vapour Partial Pressure
REAL pvs_a, pvs_b, pvs ! Water Vapour Saturation Pressure 
REAL RH                ! Rel. Humidity 
REAL uptake_oh(nlayermx) ! OH uptake coefficient 
REAL uptake_hcl(nlayermx) ! HCl uptake coefficient
REAL v_oh, v_hcl       ! Mean thermal velocities

! Dust Parameters
! ---------------
REAL, PARAMETER :: dust_density = 2.5e3
REAL, PARAMETER :: NA = 6.022e23 ! Avogadro's constant 

! Tracer Indexes in Chemistry
! ===========================
integer, parameter :: i_co2  =  1
integer, parameter :: i_co   =  2
integer, parameter :: i_o    =  3
integer, parameter :: i_o1d  =  4
integer, parameter :: i_o2   =  5
integer, parameter :: i_o3   =  6
integer, parameter :: i_h    =  7
integer, parameter :: i_h2   =  8
integer, parameter :: i_oh   =  9
integer, parameter :: i_ho2  = 10
integer, parameter :: i_h2o2 = 11
!      Methane Oxidation
integer, parameter :: i_ch4  = 12
integer, parameter :: i_ch3  = 13 
integer, parameter :: i_ch3o2 = 14
integer, parameter :: i_ch3ooh = 15
integer, parameter :: i_ch3oh= 16
integer, parameter :: i_ch3o  = 17
integer, parameter :: i_hcho = 18
integer, parameter :: i_hcooh = 19
integer, parameter :: i_hoch2o2 = 20
integer, parameter :: i_hoch2oh = 21 
integer, parameter :: i_hoch2ooh = 22 
integer, parameter :: i_hco = 23
!      Alkane Oxidation
integer, parameter :: i_c2h6 = 24
integer, parameter :: i_c2h5 = 25
integer, parameter :: i_c2h5o2 = 26
integer, parameter :: i_c2h5ooh = 27
integer, parameter :: i_c2h5oh = 28
integer, parameter :: i_hoch2ch2o2 = 29
integer, parameter :: i_hoch2ch2o = 30
integer, parameter :: i_ethgly = 31
integer, parameter :: i_hyetho2h = 32
integer, parameter :: i_ch3cho = 33
integer, parameter :: i_ch3choho2 = 34
integer, parameter :: i_ch3cooh = 35
integer, parameter :: i_ch3chohooh = 36
integer, parameter :: i_ch3co = 37
integer, parameter :: i_ch3cooo = 38
integer, parameter :: i_ch3coooh = 39
integer, parameter :: i_hcoch2o2 = 40
integer, parameter :: i_glyox = 41
integer, parameter :: i_hcoco = 42
integer, parameter :: i_hooch2cho = 43
integer, parameter :: i_hoch2cho = 44
integer, parameter :: i_hochcho = 45
integer, parameter :: i_hoch2co = 46
integer, parameter :: i_hoch2co3 = 47
integer, parameter :: i_hoch2co2h = 48
integer, parameter :: i_hcoco2h = 49
integer, parameter :: i_hoch2co3h = 50
integer, parameter :: i_hcoco3h = 51
integer, parameter :: i_hcoco3 = 52
integer, parameter :: i_ch2choh = 53
!      Water, nitrogen, and Families
integer, parameter :: i_h2o  = 54
integer, parameter :: i_n2   = 55
integer, parameter :: i_hox  = 56
integer, parameter :: i_ox   = 57
integer, parameter :: i_RO2  = 58
integer, parameter :: i_dust = 59
!      Chlorine Compounds 
integer, parameter :: i_cl = 60
integer, parameter :: i_clo = 61
integer, parameter :: i_cl2 = 62
integer, parameter :: i_oclo = 63
integer, parameter :: i_cl2o2 = 64
integer, parameter :: i_hcl = 65
integer, parameter :: i_hocl = 66 
integer, parameter :: i_cloo = 67 
integer, parameter :: i_ch3ocl = 68
integer, parameter :: i_clco = 69
integer, parameter :: i_clo3 = 70
integer, parameter :: i_hclo4 = 71
integer, parameter :: i_clo4 = 72
integer, parameter :: i_clox = 73





DO l = 1, nlayermx 

    ! =========================
    ! 1. : Dust Number Density 
    ! =========================

    ! 1.1 : Particle Volume 
    ! ---------------------
    particle_volume = (4./3.)*pi*(rdust(l)**3) ! m^3

    ! 1.2 : Mass Concentration 
    ! ------------------------
    !   - cc(:,i_dust) is in MMR units, remainder are number densities 
    mass_conc = cc(l,i_dust)*press(l)/(8.314*temp(l)) 

    ! 1.3 : Volume Concentration
    ! --------------------------
    vol_conc = mass_conc/dust_density

    ! 1.4 : Dust Number Density 
    ! -------------------------
    cc(l,i_dust) = 1.e-6*vol_conc/particle_volume

    ! ======================
    ! 2. : Dust Surface Area 
    ! ======================
    s = 4.*pi*(100.*rdust(l))**2
    s = s*cc(l,i_dust)

    ! =========================
    ! 3 : OH Uptake Coefficient
    ! =========================
    !   - coefficient from https://doi.org/10.1021/jp311235h
    !
    ! 3.1 : Partial Pressure of Water Vapour
    ! --------------------------------------
    pv = cc(l,i_h2o)*temp(l)*kb*1.e6 

    ! 3.2 : Saturation Pressure
    ! -------------------------
    pvs_a =  6816.*( (1./273.15) - (1./temp(l)) ) 
    pvs_b =  5.1309*LOG(273.15/temp(l))

    pvs = 6.112*EXP(pvs_a + pvs_b)

    ! 3.3 : OH Uptake 
    ! ---------------
    uptake_oh(l) = 0.2/( 1. + RH**0.36 )

    ! ============================================
    ! 4. : HCl Uptake Coeffient [constant for now]
    ! ============================================
    uptake_hcl(l) = 1.e-2 

    ! =======================
    ! 5. : Thermal Velocities 
    ! =======================
    v_oh = (8./pi)*kb*temp(l)*NA/mmol(igcm_oh)
    v_oh = SQRT(v_oh)

    v_hcl = (8./pi)*kb*temp(l)*NA/mmol(igcm_hcl)
    v_hcl = SQRT(v_hcl) 

    ! =========================
    ! 6. Loss/Production Fluxes 
    ! =========================

    ! 6.1: OH Reactions with dust 
    ! ---------------------------
    dust001(l) = 0.25*uptake_oh(l)*s*v_oh

    ! 6.2: HCl Uptake 
    ! ----------------
    dust002(l) = 0.25*uptake_hcl(l)*s*v_hcl

ENDDO ! l 

RETURN



END SUBROUTINE cl_dust


! ! INPUT 
! REAL q(nlayermx,nqmx) ! Tracer Mixing Ratios 
! REAL reff(nlayermx) ! Dust Effective Radius (m)
! REAL temp(nlayermx) ! Atmospheric temperature (k)
! REAL press(nlayermx) ! Atmospheric Pressure (hPa)
! ! LOCAL 
! INTEGER iq, l ! tracer and layer iterators 
! REAL dust_numdens(nlayermx) ! Dust Number Density 
! REAL mass_conc, vol_conc ! Mass and Volume Concentrations
! REAL particle_volume(nlayermx) 
! REAL numdens(nlayermx) ! Atmospheric Number Density

! REAL dustmass ! Mass (g) of dust in layer
! REAL cl_dustmass ! Mass (g) of Cl within dust 
! REAL N_cldust(nlayermx) ! Number of Cl atoms in dust in layer 





! ! ---------------------------------------------------------------------------------------!

! SUBROUTINE cl_dust_tendency(N_cldust,q, dq, reff, &
!                             press, temp, dt) 

! IMPLICIT NONE 

! #include "dimensions.h"
! #include "dimphys.h"
! #include "chimiedata.h"
! #include "tracer.h"
! #include "comcstfi.h"
! #include "callkeys.h"
! #include "conc.h"

! REAL N_cldust(nlayermx) ! Chlorine atoms in dust
! REAL q(nlayermx,nqmx), dq(nlayermx,nqmx) ! Tracer Mixing Ratios and tendencies 
! REAL reff(nlayermx) ! Dust Geometric radius (m)
! REAL press(nlayermx), temp(nlayermx) ! Pressure and temperature 
! REAL dt ! physical timestep

! ! LOCAL 
! REAL, PARAMETER :: release_rate = 1.E1 ! Release rate of Cl atoms from dust [molecules/s]
! REAL uptake_coeff(nlayermx) ! Uptake coefficient on HCl on dust 
! REAL S(nlayermx) ! surface area of layer
! REAL v_hcl, v_oh ! Mean thermal speed of HCl
! REAL uptake_hcl ! HCl dust uptake coefficient
! REAL hcl_nd(nlayermx), oh_nd(nlayermx), cl_nd(nlayermx) ! Number densities
! real hcl_nd0(nlayermx), oh_nd0(nlayermx), cl_nd0(nlayermx), N_cldust0(nlayermx) ! Initial Number Densities
! REAL zdens(nlayermx) ! Atmospheric Number Density  

! REAL d_cl(nlayermx), d_hcl(nlayermx), d_oh(nlayermx), d_cldust(nlayermx)

! REAL, PARAMETER :: NA = 6.022e23 ! Avogadro's constant 
! REAL, PARAMETER :: cl_molar = 35.4 ! Cl Molar Mass

! REAL, PARAMETER :: uptake_oh = 1.e0

! INTEGER l

! LOGICAL, SAVE :: firstcall = .TRUE. 


! DO l =  1, nlayermx

!     !  Atmospheric Number Density 
!     zdens(l) = press(l)/(temp(l)*kb*1.e6)

!     ! Tracer Number densities
!     ! -----------------------
!     ! HCl
!     hcl_nd0(l) = (q(l,igcm_hcl) + dq(l,igcm_hcl)*dt)*zdens(l)*mmean(1,l)/mmol(igcm_hcl)
!     ! OH 
!     oh_nd0(l) = (q(l,igcm_oh) + dq(l,igcm_oh)*dt)*zdens(l)*mmean(1,l)/mmol(igcm_oh) 
!     ! Cl
!     cl_nd0(l) =  (q(l,igcm_cl) + dq(l,igcm_cl)*dt)*zdens(l)*mmean(1,l)/mmol(igcm_cl) 
!     ! Cl in dust
!     N_cldust0(l) = N_cldust(l)


!     ! Surface Area of Dust (cm2 per cm-3)
!     ! -----------------------------------
!     s(l) = 4.*pi*(reff(l)/1.e-2)**2

!     ! Mean Molecule Thermal Velocities
!     ! --------------------------------
!     v_hcl = SQRT( kb*temp(l)*NA/mmol(igcm_hcl) ) 
!     v_oh = SQRT( kb*temp(l)*NA/mmol(igcm_oh) ) 

!     ! HCl Uptake Coefficient 
!     ! ----------------------
!     !  - methane equation used from https://doi.org/10.1016/j.icarus.2009.11.030
!     !     Gough et al. 2010 
!     uptake_hcl = exp( -45. + (18080.)/(8.3145*temp(l)))




!     ! Update the Chlorine Reservoir 
!     ! -----------------------------
!     N_cldust(l) = N_cldust0(l) + d_cldust(l)*dt 

!     ! Update Cl 
!     cl_nd(l) = cl_nd0(l) + d_cl(l)*dt 
!     ! Update OH
!     oh_nd(l) = oh_nd0(l) + d_oh(l)*dt 
!     ! Update HCl 
!     hcl_nd(l) = hcl_nd0(l) + d_hcl(l)*dt 

!     ! Calculation of tendencies for tracer mixing ratios 

!     ! Cl 
!     ! --
!     dq(l,igcm_cl) = dq(l,igcm_cl) + (  cl_nd(l) - cl_nd0(l)  )*mmol(igcm_cl) &
!                     /(mmean(1,l)*zdens(l)*dt)

!     ! HCl 
!     ! ---
!     dq(l,igcm_hcl) = dq(l,igcm_hcl) + (  hcl_nd(l) - hcl_nd0(l)  )*mmol(igcm_hcl) &
!                     /(mmean(1,l)*zdens(l)*dt)

!     ! OH 
!     ! --
!     dq(l,igcm_oh) = dq(l,igcm_oh) + (  oh_nd(l) - oh_nd0(l)  )*mmol(igcm_oh) &
!                     /(mmean(1,l)*zdens(l)*dt)


!     write(*,*)

! ENDDO 


! stop

! END SUBROUTINE cl_dust_tendency