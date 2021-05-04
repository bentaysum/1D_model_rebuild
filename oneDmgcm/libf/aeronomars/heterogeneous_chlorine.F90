SUBROUTINE heterogeneous_chlorine(cc, nesp, dens, & 
                                 temp, press, &  
                                 ice_mmr, &              
                                 rdust, rice, &
                                 dust001, dust002, ice001)

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
REAL ice_mmr(nlayermx) ! Water Ice Mass Mixing Ratio (kg/kg)
REAL rdust(nlayermx), rice(nlayermx) ! Mean Geometric radius of dust/ice particles [m]
REAL temp(nlayermx)    ! temperature [K]
REAL press(nlayermx)   ! Pressure [Pa]
! Output 
! ------
REAL dust001(nlayermx) ! Rate of OH + dust -> Cl + Products
REAL dust002(nlayermx) ! Rate of HCl adsorption onto dust 
REAL ice001(nlayermx)  ! Rate of HCl adsorption onto ice 

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
REAL uptake_hcl(nlayermx,2) ! HCl uptake coefficient on dust (1) and ice(2)
REAL v_oh, v_hcl       ! Mean thermal velocities
REAL v_h, v_ho2
REAL ice_nd            ! Ice Particle Number Density
REAL airdens           ! Air Density (g/cm3)

! Dust Parameters
! ---------------
REAL, PARAMETER :: dust_density = 1.5 ! g/cm^3
REAL, PARAMETER :: ice_density = 0.9167 ! g/cm^3
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

! Scalar Maximum Uptake on ice  
! -----------------------------
REAL, PARAMETER :: gamma_oh = 9.e-1


DO l = 1, nlayermx 

    airdens = 1.e-3*press(l)/(rnew(1,l)*temp(l)) 

    ! =========================
    ! 1. : Dust Number Density 
    ! =========================

    ! 1.1.0 : Prevents division by 0
    ! ----------------------------
    IF ( rdust(l) < 1.e-12 ) THEN 
        cc(l,i_dust) = 0. 
    ELSE

    ! 1.1.1 : Particle Volume 
    ! ---------------------
    particle_volume = 1.e6*(4./3.)*pi*(rdust(l)**3) ! cm^3

    ! 1.1.2 : Mass Concentration 
    ! ------------------------
    !   - cc(:,i_dust) is in MMR units, remainder are number densities 
    mass_conc = cc(l,i_dust)*airdens

    ! 1.1.3 : Volume Concentration
    ! --------------------------
    vol_conc = mass_conc/dust_density

    ! 1.1.4 : Dust Number Density 
    ! -------------------------
    cc(l,i_dust) = vol_conc/particle_volume

    ! 1.1.5 : Small values suppressed to 0
    ! ----------------------------------
    IF ( cc(l,i_dust) < 1.e-30*dens(l) ) cc(l,i_dust) = 0. 

    ENDIF 

    ! ======================
    ! 1.2.0 : Dust Surface Area 
    ! ======================
    s = 4.*pi*(100.*rdust(l))**2
    s = s*cc(l,i_dust)

    ! =========================
    ! 1.3.0 : OH Uptake Coefficient
    ! =========================
    !   - coefficient from https://doi.org/10.1021/jp311235h
    !
    ! 1.3.1 : Partial Pressure of Water Vapour
    ! --------------------------------------
    pv = cc(l,i_h2o)*temp(l)*kb*1.e6 

    ! 1.3.2 : Saturation Pressure
    ! -------------------------
    pvs_a =  6816.*( (1./273.15) - (1./temp(l)) ) 
    pvs_b =  5.1309*LOG(273.15/temp(l))

    pvs = 6.112*EXP(pvs_a + pvs_b)

    ! 1.3.3 : OH Uptake 
    ! ---------------
    uptake_oh(l) = gamma_oh/( 1. + RH**0.36 )

    ! ============================================
    ! 1.4.0 : HCl Uptake Coeffient [constant for now]
    ! ============================================
    uptake_hcl(l,1) = 0.

    ! =======================
    ! 1.5.0 : Thermal Velocities 
    ! =======================
    v_oh = (8./pi)*kb*temp(l)*NA/mmol(igcm_oh)
    v_oh = SQRT(v_oh)

    v_h = (8./pi)*kb*temp(l)*NA/mmol(igcm_h)
    v_h = SQRT(v_h)

    v_ho2 = (8./pi)*kb*temp(l)*NA/mmol(igcm_ho2)
    v_ho2 = SQRT(v_ho2)

    v_hcl = (8./pi)*kb*temp(l)*NA/mmol(igcm_hcl)
    v_hcl = SQRT(v_hcl) 

    ! =========================
    ! 1.6.0 Loss/Production Fluxes 
    ! =========================

    ! 1.6.1: OH Reactions with dust 
    ! ---------------------------
    dust001(l) = 0.25*uptake_oh(l)*s*( v_oh*cc(l,i_oh) &
                                  + v_h*cc(l,i_h) &
                                  + v_ho2*cc(l,i_ho2) )

    ! 1.6.2: HCl Uptake 
    ! ----------------
    dust002(l) = 0.25*uptake_hcl(l,1)*s*v_hcl



    ! =========================
    ! 2.0 : Water Ice Particles 
    ! =========================
    IF ( rice(l) < 1.e-12 ) THEN 
        ice_nd = 0. 
    ELSE

    ! Volume of Ice Particle cm^3
    ! ---------------------------
    particle_volume = (4./3.)*pi*(rice(l)**3)*1.e6 

    ! 2.0.1 : Mass Concentration
    ! --------------------------
    mass_conc = airdens*ice_mmr(l)*mmol(igcm_h2o_ice)/mmean(1,l)

    ! 2.0.2 : Volume Concentration
    ! ----------------------------
    vol_conc = mass_conc/ice_density

    ! 2.0.3 : Number Density 
    ! ----------------------
    ice_nd = vol_conc/particle_volume

    ENDIF 

    ! 2.1 : Available Ice Surface Area (cm^2)
    ! ---------------------------------------
    s = 4.*pi*( rice(l)*100. )**2 

    ! 2.2 : Uptake Coeffient 
    ! ----------------------
    uptake_hcl(l,2) = 0.1

    ! 2.3 : Rate Coeffient
    ! --------------------
    ice001(l) = 0.25*uptake_hcl(l,2)*s*v_hcl


ENDDO ! l 



RETURN



END SUBROUTINE heterogeneous_chlorine