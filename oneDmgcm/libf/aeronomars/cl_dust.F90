SUBROUTINE cl_dust_initialise(q,reff, &
                            temp,press, &
                            N_cldust)

IMPLICIT NONE


#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "conc.h"

! INPUT 
REAL q(nlayermx,nqmx) ! Tracer Mixing Ratios 
REAL reff(nlayermx) ! Dust Effective Radius (m)
REAL temp(nlayermx) ! Atmospheric temperature (k)
REAL press(nlayermx) ! Atmospheric Pressure (hPa)
! LOCAL 
INTEGER iq, l ! tracer and layer iterators 
REAL dust_numdens(nlayermx) ! Dust Number Density 
REAL mass_conc, vol_conc ! Mass and Volume Concentrations
REAL particle_volume(nlayermx) 
REAL numdens(nlayermx) ! Atmospheric Number Density

REAL dustmass ! Mass (g) of dust in layer
REAL cl_dustmass ! Mass (g) of Cl within dust 
REAL N_cldust(nlayermx) ! Number of Cl atoms in dust in layer 

REAL, PARAMETER :: dustdens = 2.5e3 ! Dust Density (kg/m-3)
REAL, PARAMETER :: cl_wt = 0.5 ! Weight Percentage of Cl of dust 
REAL, PARAMETER :: NA = 6.022e23 ! Avogadro's constant 
REAL, PARAMETER :: cl_molar = 35.4 ! Cl Molar Mass

! ============================================
! Step 1 : Dust Mixing Ratio -> Number Density
! ============================================
DO iq = 1, nqmx
    IF ( trim(noms(iq)) == "dust_mass" ) THEN 
        
        DO l = 1, nlayermx  
            ! Particle Geometric Volume 
            particle_volume(l) = (4./3.)*pi*reff(l)**3
            ! Mass Concentration
            mass_conc =  q(l,iq)*press(l) &
                      /(8.314*temp(l))
            ! Volume Concentration
            vol_conc = mass_conc/dustdens 
            ! Particle Number Density 
            dust_numdens(l) = 1.e-6*vol_conc/particle_volume(l)

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
    N_cldust(l) = cl_dustmass*NA/cl_molar

ENDDO


END SUBROUTINE cl_dust_initialise


! ---------------------------------------------------------------------------------------!

SUBROUTINE cl_dust_tendency(N_cldust,q, dq, reff, &
                            press, temp, dt) 

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "conc.h"

REAL N_cldust(nlayermx) ! Chlorine atoms in dust
REAL q(nlayermx,nqmx), dq(nlayermx,nqmx) ! Tracer Mixing Ratios and tendencies 
REAL reff(nlayermx) ! Dust Geometric radius (m)
REAL press(nlayermx), temp(nlayermx) ! Pressure and temperature 
REAL dt ! physical timestep
! LOCAL 
REAL, PARAMETER :: release_rate = 1.E1 ! Release rate of Cl atoms from dust [molecules/s]

REAL uptake_coeff(nlayermx) ! Uptake coefficient on HCl on dust 
REAL S(nlayermx) ! surface area of layer
REAL v_hcl, v_oh ! Mean thermal speed of HCl
REAL uptake_hcl ! HCl dust uptake coefficient
REAL hcl_nd(nlayermx), oh_nd(nlayermx) ! Number densities
REAL zdens(nlayermx) ! Atmospheric Number Density 

REAL d_cl(nlayermx), d_hcl(nlayermx), d_oh(nlayermx), d_cldust(nlayermx)

REAL, PARAMETER :: NA = 6.022e23 ! Avogadro's constant 
REAL, PARAMETER :: cl_molar = 35.4 ! Cl Molar Mass

REAL, PARAMETER :: uptake_oh = 1.e0

INTEGER l

LOGICAL, SAVE :: firstcall = .TRUE. 
REAL N_cldust0(nlayermx) ! Initial Reservoir (maximum permittable values)


IF ( firstcall ) THEN 
    N_cldust0 = N_cldust
ENDIF 

! Initialise by updating the mixing ratio with their prior calculated tendencies 
q = q + dq*dt

DO l =  1, nlayermx

    !  Atmospheric Number Density 
    zdens(l) = press(l)/(temp(l)*kb*1.e6)

    ! Tracer Number densities
    ! -----------------------
    ! HCl
    hcl_nd(l) = zdens(l)*q(l,igcm_hcl)*43.34/mmol(igcm_hcl)
    ! OH 
    oh_nd(l) = zdens(l)*q(l,igcm_oh)*43.34/mmol(igcm_oh) 

    ! Surface Area of Dust (cm2 per cm-3)
    s(l) = 4.*pi*(reff(l)/1.e-2)**2

    ! Mean Molecule Thermal Velocities
    v_hcl = SQRT( kb*temp(l)*NA/mmol(igcm_hcl) ) 
    v_oh = SQRT( kb*temp(l)*NA/mmol(igcm_oh) ) 

    ! HCl Uptake Coefficient 
    ! ----------------------
    !  - methane equation used from https://doi.org/10.1016/j.icarus.2009.11.030
    !     Gough et al. 2010 
    uptake_hcl = exp( -45. + (18080.)/(8.3145*temp(l)))

    ! Tendency of Cl released from dust [molec. s-1 cm-3]
    d_cl(l) = 0.25*uptake_oh*s(l)*v_oh*oh_nd(l)

    ! Tendency of HCl [molec. s-1 cm-3] 
    d_hcl(l) = -0.25*uptake_hcl*s(l)*v_hcl*hcl_nd(l)

    ! Tendency of Cl trapped in dust [molec. s-1 cm-3] 
    d_cldust(l) = d_hcl(l) - d_cl(l)

    ! Tendency of OH [molec. s-1 cm-3]
    d_oh(l) = -d_cl(l)


    ! Update the Chlorine Reservoir 
    N_cldust(l) = N_cldust(l) + d_cldust(l)*1800.

    write(*,*) N_cldust0(l), N_cldust(l), 100.*(N_cldust(l)-N_cldust0(l))/N_cldust0(l)
ENDDO 


stop

END SUBROUTINE cl_dust_tendency