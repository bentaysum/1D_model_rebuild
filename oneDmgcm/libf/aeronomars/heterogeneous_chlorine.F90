SUBROUTINE heterogeneous_chlorine(cc, nesp, dens, & 
                                 temp, press, &  
                                 ice_mmr, &              
                                 rdust, rice, &
                                 dust001, dust002, dust003, dust004,&
                                 ice001, ice002, ice003, ice004)

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

! Output [2nd order rate coefficients]
! ------------------------------------
real dust001(nlayermx) ! OH + dust -> ClOx + product 
real, parameter :: d1_up_gamma =8.E-1
real d1_up
real dust002(nlayermx) ! HO2 + dust -> ClOx + product 
real, parameter :: d2_up = 0.02
real dust003(nlayermx) ! H2O + dust -> ClOx + product
real, parameter :: d3_up = 0.
real dust004(nlayermx) ! Cl + dust -> Products 
real, parameter :: d4_up = 1.e-2

real ice001(nlayermx)  ! HCl + ice -> Products
real, parameter :: i1_up = 0.8
real ice002(nlayermx)  ! Cl2 + ice -> products
real, parameter :: i2_up = 0. !1.e-4
real ice003(nlayermx)  ! OClO + ice -> Products
real, parameter :: i3_up = 0.! 0.8
real ice004(nlayermx)  ! ClO + ice -> Products 
real, parameter :: i4_up = 0. !1.e-4

! Local 
! =====
INTEGER l              ! Layer iterator 
REAL particle_volume   ! Volume of single dust grain 
REAL mass_conc         ! Dust/ice mass concentration 
REAL vol_conc          ! Dust/ice volume concentration
REAL dustsurf, icesurf ! Dust/ice Surface Area per Unit volume
REAL pv                ! Water Vapour Partial Pressure
REAL pvs_a, pvs_b, pvs ! Water Vapour Saturation Pressure 
REAL RH                ! Rel. Humidity 
REAL ice_nd            ! Ice Particle Number Density
REAL airdens           ! Air Density (g/cm3)
REAL v_therm           ! Thermal velocity of molecule 

! Uptake Coefficients
! ===================

! 

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


DO l = 1, nlayermx ! Altitude step loop 
! ===================================================== !
! Stage 1 : Calculation of Surface Area per unit volume !
!           of dust and ice 
! ===================================================== !
    
    ! 1.0 : Air Density (g/cm3)
    ! -------------------------
    airdens = 1.e-3*press(l)/(rnew(1,l)*temp(l)) 

    ! +++++++++++++++++++++++++++++++++++++
    ! 1.1.0 : Dust Surface Area (cm^2/cm^3)
    ! +++++++++++++++++++++++++++++++++++++

    ! 1.1.1 : Suppress division by 0
    ! ------------------------------
    IF ( rdust(l) < 1.e-12 ) THEN 
        cc(l,i_dust) = 0. 
    ELSE

    ! 1.1.2 : Particle Volume 
    ! ---------------------
    particle_volume = 1.e6*(4./3.)*pi*(rdust(l)**3) ! cm^3

    ! 1.1.3 : Mass Concentration 
    ! ------------------------
    !   - cc(:,i_dust) is in MMR units on input,
    !     remainder are trace gas number densities 
    mass_conc = cc(l,i_dust)*airdens

    ! 1.1.4 : Volume Concentration
    ! --------------------------
    vol_conc = mass_conc/dust_density

    ! 1.1.5 : Dust Number Density 
    ! -------------------------
    cc(l,i_dust) = vol_conc/particle_volume

    ! 1.1.6 : Small values suppressed to 0
    ! ----------------------------------
    IF ( cc(l,i_dust) < 1.e-30*dens(l) ) cc(l,i_dust) = 0. 

    ENDIF ! if rdust < 1.e-12

    ! 1.1.7 : Dust Surface Area per cm^3 
    ! ----------------------------------
    dustsurf = cc(l,i_dust)*4.*pi*(100.*rdust(l))**2


    ! +++++++++++++++++++++++++++++++++++++
    ! 1.2.0 : Ice Surface Area (cm^2/cm^3)
    ! +++++++++++++++++++++++++++++++++++++

    ! 1.2.1 : Suppress division by 0
    ! ------------------------------
    IF ( rice(l) < 1.e-12 ) THEN 
        ice_nd = 0. 
    ELSE

    ! 1.2.2 : Particle Volume
    ! -----------------------
    particle_volume = (4./3.)*pi*(rice(l)**3)*1.e6 

    ! 1.2.3 : Mass Concentration
    ! --------------------------
    mass_conc = airdens*ice_mmr(l)*mmol(igcm_h2o_ice)/mmean(1,l)

    ! 1.2.4 : Volume Concentration
    ! ----------------------------
    vol_conc = mass_conc/ice_density

    ! 1.2.5 : Number Density 
    ! ----------------------
    ice_nd = vol_conc/particle_volume

    ! 1.2.6 : Small values suppressed to 0
    ! ----------------------------------
    IF ( ice_nd < 1.e-30*dens(l) ) ice_nd = 0. 

    ENDIF ! if rdust < 1.e-12


    ! 1.2.6 : Available Ice Surface Area (cm^2)
    ! ---------------------------------------
    icesurf = ice_nd*4.*pi*( rice(l)*100. )**2 


    ! ===============================================
    ! 2.0 : Rate Coefficients 
    ! ===============================================

    ! ------------------------------
    ! 2.1 : OH + Dust -> ClOx + Dust 
    ! ------------------------------

    ! 2.1.0 : OH Uptake Coefficient dependancy on Partial Pressure
    ! ------------------------------------------------------------ 

    ! 2.1.1 : Partial Pressure of Water Vapour
    !       - Buck Eq.
    ! --------------------------------------
    pv = 18.678 - ((temp(l) - 273.15)/234.5) 
    pv = pv*( (temp(l)-273.15)/( 257.14 - 273.15 + temp(l) ) )
    pv = 0.61121*EXP( pv )*1.e-3 ! Pa

    ! 2.1.2 : Saturation Pressure
    ! ---------------------------
    pvs_a = 2.07023-0.00320991*temp(l)-2484.896/temp(l)+3.56654*alog10(temp(l))
    pvs = (10.**pvs_a)*100. ! Pa
    
    ! 2.1.3 : Relative Humidity
    ! -------------------------
    RH = 100.*pv/pvs

    ! 2.1.4 : Uptake Coefficient
    ! --------------------------
    d1_up = d1_up_gamma/(1. + RH**0.36)

    ! 2.1.5 : Mean Thermal Velocity of OH
    ! -----------------------------------
    v_therm = 1.e7*(8./pi)*kb*temp(l)*NA/mmol(igcm_oh)
    v_therm = SQRT(v_therm)

    ! 2.1.6 : Rate
    ! ------------
    dust001(l) = 0.25*d1_up*dustsurf*v_therm

    ! -------------------------------
    ! 2.2 : HO2 + Dust -> ClOx + Dust 
    ! -------------------------------

    ! 2.2.1 : Mean Thermal Vel. of HO2
    ! --------------------------------
    v_therm = 1.e7*(8./pi)*kb*temp(l)*NA/mmol(igcm_ho2)
    v_therm = SQRT(v_therm)

    ! 2.2.2 : Rate
    ! ------------
    dust002(l) = 0.25*d2_up*dustsurf*v_therm

    ! -------------------------------
    ! 2.3 : H2O + Dust -> ClOx + Dust 
    ! -------------------------------

    ! 2.3.1 : Mean Thermal Vel. of HO2
    ! --------------------------------
    v_therm = 1.e7*(8./pi)*kb*temp(l)*NA/mmol(igcm_h2o_vap)
    v_therm = SQRT(v_therm)

    ! 2.3.2 : Rate
    ! ------------
    dust003(l) = 0.25*d3_up*dustsurf*v_therm

    ! -------------------------------
    ! 2.4 : Cl + Dust -> Products  
    ! -------------------------------

    ! 2.4.1 : Mean Thermal Vel. of Cl
    ! --------------------------------
    v_therm = 1.e7*(8./pi)*kb*temp(l)*NA/mmol(igcm_cl)
    v_therm = SQRT(v_therm)

    ! 2.4.2 : Rate
    ! ------------
    dust004(l) = 0.25*d4_up*dustsurf*v_therm

    ! ----------------------------
    ! 2.5 : HCl + Ice -> Products
    ! ----------------------------

    ! 2.5.1 : Thermal Vel.
    ! -------------------
    v_therm = 1.e7*(8./pi)*kb*temp(l)*NA/mmol(igcm_hcl)
    v_therm = SQRT(v_therm)

    ! 2.5.2 : Rate
    ! ------------
    ice001(l) = 0.25*i1_up*icesurf*v_therm

    ! ----------------------------
    ! 2.6 : Cl2 + Ice -> Products
    ! ----------------------------

    ! 2.6.1 : Thermal Vel.
    ! -------------------
    v_therm = 1.e7*(8./pi)*kb*temp(l)*NA/mmol(igcm_cl2)
    v_therm = SQRT(v_therm)

    ! 2.6.2 : Rate
    ! ------------
    ice002(l) = 0.25*i2_up*icesurf*v_therm

    ! ----------------------------
    ! 2.7 : OClO + Ice -> Products
    ! ----------------------------

    ! 2.7.1 : Thermal Vel.
    ! -------------------
    v_therm = 1.e7*(8./pi)*kb*temp(l)*NA/mmol(igcm_oclo)
    v_therm = SQRT(v_therm)

    ! 2.7.2 : Rate
    ! ------------
    ice003(l) = 0.25*i3_up*icesurf*v_therm

    ! ----------------------------
    ! 2.6 : ClO + Ice -> Products
    ! ----------------------------

    ! 2.6.1 : Thermal Vel.
    ! -------------------
    v_therm = 1.e7*(8./pi)*kb*temp(l)*NA/mmol(igcm_clo)
    v_therm = SQRT(v_therm)

    ! 2.6.2 : Rate
    ! ------------
    ice004(l) = 0.25*i4_up*icesurf*v_therm


ENDDO 



RETURN



END SUBROUTINE heterogeneous_chlorine