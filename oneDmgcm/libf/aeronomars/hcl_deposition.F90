SUBROUTINE hcl_deposition(ice, rice, &
                          dens,  temp, press, &
                          ice001)

! ==================================================================
! Calculates the removal of atmospheric Hydrogen Chloride 
! due to deposition onto water ice particles.
!   Heterogeneous interactions (i.e. HCl removal by dust deposition)
! is handled in cl_dust.F90.
! ==================================================================

IMPLICIT NONE 


#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "conc.h"

! =====
! INPUT 
! =====
REAL ice(nlayermx) ! Mixing Ratio  of Ice 
REAL rice(nlayermx)               ! Geometric Radius of Ice Particles 
REAL dens(nlayermx)               ! Atmospheric Number Density
REAL temp(nlayermx)               ! Atmospheric Temperature
REAL press(nlayermx)              ! Pressure (Pa on input)

! ======
! OUTPUT
! ======
REAL ice001(nlayermx)

! =====
! LOCAL
! =====
INTEGER l ! Layer iterator
REAL icedens(nlayermx) ! Atmospheric Number density of water ice particles

REAL particle_volume
REAL mass_conc
REAL vol_conc
REAL s
REAL nuhcl

REAL, PARAMETER :: HCL_UPTAKE = 0.02 

REAL, PARAMETER :: ice_particle_density = 916.7 ! kg/m^3
REAL, PARAMETER :: NA = 6.022e23 ! Avogadro's constant 

! ==================================================================
! 1.0: Ice Particle VMR to Ice Particle Number Density
! ==================================================================
DO l = 1, nlayermx

    ! 1.1: Single particle vol. 
    ! ------------------------
    particle_volume = (4./3.)*pi*(rice(l))**3 ! m^3

    ! 1.2: Mass Concentration of Ice
    ! ------------------------------
    mass_conc = (ice(l)*mmol(igcm_h2o_ice)/mmean(1,l)) & ! VMR -> MMR conversion
            *press(l)/(8.314*temp(l))

    ! 1.3: Volume Concentration of Ice 
    ! --------------------------------
    vol_conc = mass_conc/ice_particle_density

    ! 1.4: Ice Particle Number Density
    ! --------------------------------
    IF (rice(l) < 1.e-20) THEN
        ice001(l) = 0. ! Prevents issues with divisions by 0
        continue 
    ELSE
        icedens(l) = 1.e-6*vol_conc/particle_volume
    ENDIF

! ====================================
! 2.0: Water Ice Particle Surface Area 
! ====================================
    s = 4.*pi*(rice(l)*100.)**2 ! cm^2 
    s = s*icedens(l)

! ==========================
! 3.0 : HCl Thermal Velocity 
! ==========================
    nuhcl = (8./pi)*kb*temp(l)*NA/mmol(igcm_hcl)
    nuhcl = SQRT(nuhcl)

! =====================
! 4.0 : Deposition Rate  
! =====================
    ice001(l) = 0.25*HCL_UPTAKE*s*nuhcl

ENDDO 






END SUBROUTINE