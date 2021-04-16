SUBROUTINE cl_dust(vmr, zdens, reff, temp, press, dt)

IMPLICIT NONE


#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "conc.h"

! ===============
! Input Variables 
! ===============
REAL vmr(nlayermx,nqmx)     ! Tracer Number Density
REAL zdens(nlayermx) 
REAL reff(nlayermx)         ! Dust geometric radius (m)
REAL dustsurf(nlayermx)
REAL temp(nlayermx)
REAL press(nlayermx)
REAL dt

! ===============
! Local Variables 
! ===============
INTEGER l, iq

REAL particle_volume(nlayermx) ! Mean volume of a dust particle
REAL mass_conc                 ! Dust Mass Concentration
REAL vol_conc                   ! Dust Volume Concentration
REAL dust_numdens(nlayermx)    ! Dust number density (cm-3)
REAL dustmass, cl_dustmass

REAL S(nlayermx) ! Dust Surface are (cm^2 per cm^-3)

REAL v_hcl(nlayermx), v_oh(nlayermx), &
    v_ho2(nlayermx), v_h2o(nlayermx), &
    v_cl(nlayermx), v_cl2(nlayermx) 

REAL uptake(nlayermx,6) 
REAL alpha(6), gamma_rxn(6)

INTEGER, PARAMETER :: g_oh = 1
INTEGER, PARAMETER :: g_ho2 = 2
INTEGER, PARAMETER :: g_h2o = 3 
INTEGER, PARAMETER :: g_hcl = 4
INTEGER, PARAMETER :: g_cl = 5 
INTEGER, PARAMETER :: g_cl2 = 6

REAL d_cc(nlayermx,6)


REAL, SAVE :: dust_cl(nlayermx), dust_cl0(nlayermx)

LOGICAL, SAVE :: firstcall = .True.

REAL, PARAMETER :: dustdens = 2.5e3 ! Dust Density (kg/m-3)
REAL, PARAMETER :: cl_wt = 0.5 ! Weight Percentage of Cl of dust 
REAL, PARAMETER :: NA = 6.022e23 ! Avogadro's constant 

REAL vmr0(nlayermx,nqmx)
REAL cc(nlayermx,nqmx) 

REAL delta

REAL pv ! partial pressure of water 
REAL pvs ! saturation partial pressure of water vapour 
REAL pvs_a, pvs_b 
REAL RH ! Relative Humidity 


! ============================================
! Step 1 : Dust Mixing Ratio -> Number Density
! ============================================
DO iq = 1, nqmx
    IF ( trim(noms(iq)) == "dust_mass" ) THEN 
        
        DO l = 1, nlayermx  
            ! Particle Geometric Volume 
            particle_volume(l) = (4./3.)*pi*(reff(l)*100.)**3
            ! Mass Concentration
            mass_conc =  vmr(l,iq)*press(l) & ! vmr for dust is -actually- mmr 
                      /(8.314*temp(l))
            ! Volume Concentration
            vol_conc = mass_conc/dustdens 
            ! Particle Number Density 
            dust_numdens(l) = vol_conc/particle_volume(l)

        ENDDO  ! l 

        EXIT 
    
    ENDIF ! "dust_mass"

ENDDO ! iq 


! ===============================================
! Step 2: Total Chlorine Available in Dust (cm-3)
! ===============================================
IF ( firstcall ) THEN 
    DO l = 1, nlayermx
        ! Mass of dust in layer (g) = particle volume * density of grain * number of particles 
        dustmass = dust_numdens(l)*dustdens*1.e-2*particle_volume(l)
        ! Mass of Cl in dust in layer 
        cl_dustmass = dustmass*cl_wt/100.
        ! Number of Cl atoms in dust in layer = mass of cl in dust * Avogadro's constant / molar mass of Cl
        ! Save first value as the maximum 
            dust_cl(l) = cl_dustmass*NA/mmol(igcm_cl)
            dust_cl0(l) = dust_cl(l)

    ENDDO
ENDIF 




DO l = 1, nlayermx

    ! VMR -> Number Density
    cc(l,:) = vmr(l,:)*zdens(l)
    vmr0(l,:) = vmr(l,:)

    ! HCl Uptake Coefficient on Dust
    uptake(l,g_hcl) = 5.E-1

    ! ------------------------
    ! OH and Relative Humidity 
    ! ------------------------
    pv = cc(l,igcm_h2o_vap)*temp(l)*kb*1.e6

    pvs_a = 6816.*( (1./273.15) - (1./temp(l)) ) 
    pvs_b = 5.1309*LOG(273.15/temp(l))

    pvs = 6.112*EXP(pvs_a + pvs_b)

    RH = 100.*(pv/pvs)

    uptake(l,g_oh) = 0.2/(1. + RH**0.36 )

    ! Calculate available dust surface area (cm^2 cm^-3)
    s(l) =  dust_numdens(l)*4.*pi*(reff(l)*100.)**2

    ! Thermal Velocity of Adsorbed Species OH and HCl (cm s-1)
    v_hcl(l) = SQRT( kb*temp(l)*NA/mmol(igcm_hcl) )
    v_oh(l) = SQRT( kb*temp(l)*NA/mmol(igcm_oh) )
    v_ho2(l) = SQRT( kb*temp(l)*NA/mmol(igcm_ho2) )
    v_h2o(l) = SQRT( kb*temp(l)*NA/mmol(igcm_h2o_vap) )
    v_cl(l) = SQRT( kb*temp(l)*NA/mmol(igcm_cl) )
    v_cl2(l) = SQRT( kb*temp(l)*NA/mmol(igcm_cl2) )

    ! Tendency Calculations 
    ! =====================
    ! 
    ! OH 
    ! --
    d_cc(l,g_oh) = -0.25*uptake(l,g_oh)*s(l)*v_oh(l)!*cc(l,igcm_oh)*dt 

    cc(l,igcm_oh) = cc(l,igcm_oh)/( 1. + d_cc(l,g_oh)*dt )

    ! HCl 
    ! ---
    d_cc(l,g_hcl) = -0.25*uptake(l,g_hcl)*s(l)*v_hcl(l)!*cc(l,igcm_hcl)*dt 

    cc(l,igcm_hcl) = cc(l,igcm_hcl)/( 1. + d_cc(l,g_hcl)*dt )

    ! Cl 
    ! --
    d_cc(l,g_cl) = -d_cc(l,g_oh) 

    cc(l,igcm_cl) = cc(l,igcm_cl) + d_cc(l,g_cl)*dt 


    dust_cl(l) = dust_cl(l) &
               -d_cc(l,g_hcl) + d_cc(l,g_cl) 



    ! Suppression of chlorine build up in dust
    !   - assume that initial value is the saturation point 
    !     after which cl is not capable of developing.
    ! if ( dust_cl(l) > dust_cl0(l) ) then 
        
    !     delta = dust_cl(l) - dust_cl0(l) 


    !     dust_cl(l) = dust_cl0(l) 

    ! endif 

    ! Update Tracer VMR's

    ! OH 
    ! --
    vmr(l,igcm_oh) =  cc(l,igcm_oh)/zdens(l)
   
    ! HCl 
    ! --
    vmr(l,igcm_hcl) = cc(l,igcm_hcl)/zdens(l)

    ! Cl 
    ! --
    vmr(l,igcm_cl) = cc(l,igcm_cl)/zdens(l)




ENDDO 



! Optional Saved Output 

call WRITEDIAGFI(1,"cl_in_dust","cl_in_dust", "s", &
                  1,dust_cl)

call WRITEDIAGFI(1,"dcl_dust","dcl_dust", "s", &
                  1,d_cc(:,g_cl))

call WRITEDIAGFI(1,"dhcl_dust","dhcl_dust", "s", &
                  1,d_cc(:,g_hcl))

call WRITEDIAGFI(1,"doh_dust","doh_dust", "s", &
                  1,d_cc(:,g_oh))

call WRITEDIAGFI(1,"oh_uptake","oh_uptake", "s", &
                  1,uptake(:,g_oh))

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