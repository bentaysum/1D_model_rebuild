SUBROUTINE cl_dust(vmr, zdens, reff, dustsurf, temp, press, dt)

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

REAL uptake(6) 
REAL alpha(6), gamma_rxn(6)

INTEGER, PARAMETER :: g_oh = 1
INTEGER, PARAMETER :: g_ho2 = 2
INTEGER, PARAMETER :: g_h2o = 3 
INTEGER, PARAMETER :: g_hcl = 4
INTEGER, PARAMETER :: g_cl = 5 
INTEGER, PARAMETER :: g_cl2 = 6

REAL d_cc(6)


REAL, SAVE :: dust_cl(nlayermx), dust_cl0(nlayermx)

LOGICAL, SAVE :: firstcall = .True.

REAL, PARAMETER :: dustdens = 2.5e3 ! Dust Density (kg/m-3)
REAL, PARAMETER :: cl_wt = 0.5 ! Weight Percentage of Cl of dust 
REAL, PARAMETER :: NA = 6.022e23 ! Avogadro's constant 

REAL vmr0(nlayermx,nqmx)
REAL cc(nlayermx,nqmx) 

REAL delta

IF ( firstcall ) THEN 
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
    DO l = 1, nlayermx
        ! Mass of dust in layer (g) = particle volume * density of grain * number of particles 
        dustmass = dust_numdens(l)*dustdens*1.e-2*particle_volume(l)
        ! Mass of Cl in dust in layer 
        cl_dustmass = dustmass*cl_wt/100.
        ! Number of Cl atoms in dust in layer = mass of cl in dust * Avogadro's constant / molar mass of Cl
        dust_cl(l) = cl_dustmass*NA/mmol(igcm_cl)
        ! Save first value as the maximum 
        dust_cl0(l) = dust_cl(l)
    ENDDO

    firstcall = .False.

ENDIF ! firstcall

! =====================
! Reactions of Interest 
! =====================
! 1. OH + dust -> Cl + products 
!       alpha = [0.01,0.1]
!       gamma = [0.01,0.2]
alpha(g_oh) = 0.1
gamma_rxn(g_oh) = 0.5 

! 2. HO2 + dust -> Cl + products
!       alpha = [0.02,0.4]
!       gamma = [0.025,0.1] 
alpha(g_ho2) = 0.4 
gamma_rxn(g_ho2) = 0.1 

! 3. H2O + dust -> Cl + products 
!       alpha = 0.85
alpha(g_h2o) = 0.!1.e-1 

! 4. HCl + dust -> dust
!       alpha = [0.05,0.3]
alpha(g_hcl) = 0.95

! 5. Cl + dust -> dust 
!       gamma = 2e-4
gamma_rxn(g_cl) = 2.e-4


! 6. Cl2 + dust -> dust 
!       gamma = [1.e-3,1.e-1]
gamma_rxn(g_cl2) = 0.2

! where uptake^-1 = (1/alpha) + (1/gamma)
!
! Coefficients used from Chemical Kinetics and Photochemical Data
!   for Use in Atmospheric Studies - Issue 19

DO l = 1, nlayermx

    ! Uptake Coefficients 
    uptake(g_oh) = (1./alpha(g_oh)) + (1./gamma_rxn(g_oh))
        uptake(g_oh) = 1./uptake(g_oh)

    uptake(g_ho2) = (1./alpha(g_ho2)) + (1./gamma_rxn(g_ho2))
        uptake(g_ho2) = 1./uptake(g_ho2)

    uptake(g_h2o) = alpha(g_h2o)

    uptake(g_hcl) = alpha(g_hcl) 

    uptake(g_cl) = gamma_rxn(g_cl) 

    uptake(g_cl2) = gamma_rxn(g_cl2) 

    ! VMR -> Number Density
    cc(l,:) = vmr(l,:)*zdens(l)
    vmr0(l,:) = vmr(l,:)

    ! Calculate available dust surface area (cm^2 cm^-3)
    s(l) =  dust_numdens(l)*4.*pi*(reff(l)*100.)**2

    ! Thermal Velocity of Adsorbed Species OH and HCl (cm s-1)
    v_hcl(l) = SQRT( kb*temp(l)*NA/mmol(igcm_hcl) )*100. 
    v_oh(l) = SQRT( kb*temp(l)*NA/mmol(igcm_oh) )*100.
    v_ho2(l) = SQRT( kb*temp(l)*NA/mmol(igcm_ho2) )*100.
    v_h2o(l) = SQRT( kb*temp(l)*NA/mmol(igcm_h2o_vap) )*100.
    v_cl(l) = SQRT( kb*temp(l)*NA/mmol(igcm_cl) )*100.
    v_cl2(l) = SQRT( kb*temp(l)*NA/mmol(igcm_cl2) )*100.

    ! Tendency 
    d_cc(g_oh) = -0.25*uptake(g_oh)*s(l)*v_oh(l)*cc(l,igcm_oh)*dt 

    d_cc(g_ho2) = -0.25*uptake(g_ho2)*s(l)*v_ho2(l)*cc(l,igcm_ho2)*dt 

    d_cc(g_h2o) = -0.25*uptake(g_h2o)*s(l)*v_h2o(l)*cc(l,igcm_h2o_vap)*dt 

    d_cc(g_cl2) = -0.25*uptake(g_cl2)*s(l)*v_cl2(l)*cc(l,igcm_cl2)*dt 

    d_cc(g_hcl) = -0.25*uptake(g_hcl)*s(l)*v_hcl(l)*cc(l,igcm_hcl)*dt 

    d_cc(g_cl) = -d_cc(g_oh) - d_cc(g_ho2) - d_cc(g_h2o)

    dust_cl(l) = dust_cl(l) &
               - ( d_cc(g_cl2) + d_cc(g_hcl) + d_cc(g_cl) ) &
               + ( d_cc(g_oh) + d_cc(g_h2o) + d_cc(g_ho2) )

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
    vmr(l,igcm_oh) = ( cc(l,igcm_oh) + d_cc(g_oh) )/zdens(l)

    ! HO2 
    ! --
    vmr(l,igcm_ho2) = ( cc(l,igcm_ho2) + d_cc(g_ho2) )/zdens(l)

    ! H2O 
    ! ---
    vmr(l,igcm_h2o_vap) = ( cc(l,igcm_h2o_vap) + d_cc(g_h2o) )/zdens(l)

    ! HCl 
    ! --
    vmr(l,igcm_hcl) = ( cc(l,igcm_hcl) + d_cc(g_hcl) )/zdens(l)

    ! Cl2
    ! ---
    vmr(l,igcm_cl2) = ( cc(l,igcm_cl2) + d_cc(g_cl2) )/zdens(l)

    ! OH 
    ! --
    vmr(l,igcm_cl) = ( cc(l,igcm_cl) + d_cc(g_cl) )/zdens(l)



    !  ! Tendency of Cl released from dust [molec. s-1 cm-3]
    !  d_cl =  0.25*gamma_oh(l)*s(l)*v_oh(l)*cc(l,igcm_oh)*dt

    !  ! Tendency of HCl [molec. s-1 cm-3] 
    !  d_hcl = - 0.25*gamma_hcl(l)*s(l)*v_hcl(l)*cc(l,igcm_hcl)*dt


    ! ! Tendency of OH [molec. s-1 cm-3]
    !  d_oh =  -d_cl*dt

    !  ! Chlorine in Dust 
    !  dust_cl(l) = dust_cl(l) - d_cl - d_hcl
    !  ! Check if the dust is saturated 
    !  IF ( dust_cl0(l) - dust_cl(l) < 0. ) THEN 
    !     d_hcl = d_hcl - ( dust_cl0(l) - dust_cl(l) )
    !     dust_cl(l) = dust_cl0(l)
    ! ENDIF 

    !  vmr(l,igcm_oh) = (cc(l,igcm_oh) + d_oh)/zdens(l)
    !  vmr(l,igcm_hcl) = (cc(l,igcm_hcl) + d_hcl)/zdens(l)
    !  vmr(l,igcm_cl) = (cc(l,igcm_cl) + d_cl)/zdens(l)


ENDDO 


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