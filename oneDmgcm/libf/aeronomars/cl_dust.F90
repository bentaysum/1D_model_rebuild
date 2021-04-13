SUBROUTINE cl_dust_initialise(q,reff, &
                            temp,press,numdens)

IMPLICIT NONE


#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "conc.h"

! INPUT 
REAL q(nqmx,nlayermx) ! Tracer Mixing Ratios 
REAL reff(nlayermx) ! Dust Effective Radius (m)
REAL temp(nlayermx) ! Atmospheric temperature (k)
REAL press(nlayermx) ! Atmospheric Pressure (hPa)
REAL numdens(nlayermx) ! Atmospheric Number Density
! LOCAL 
INTEGER iq, l ! tracer and layer iterators 
REAL dustdens(nlayermx) ! Dust Number Density 
REAL particle_volume 

! Step 1 : Dust Mixing Ratio -> Number Density
!   - Geometric Radius of dust ~ 0.36*r_eff 
DO iq = 1, nqmx
    IF ( trim(noms(iq)) == "dust_mass" ) THEN 
        
        DO l = 1, nlayermx
            ! Mixing Ratio -> Mass Concentration
            dustdens(l) = q(iq,l)*press(l)*100./(Rnew(1,l)*temp(l))
            ! Mass Concentration -> Volume Concentration
            dustdens(l) = dustdens(l)/numdens(l)
            ! Volume Concentration -> Number Density
            particle_volume = 1.e-2*(4./3.)*pi*0.36*reff(l)**3
            dustdens(l) = dustdens(l)/particle_volume 
        ENDDO  ! l 

        EXIT 
    
    ENDIF ! "dust_mass"
ENDDO ! iq 

DO l = 1, nlayermx
    WRITE(*,*) l, dustdens(l)
ENDDO 

STOP


END SUBROUTINE cl_dust_initialise


! ---------------------------------------------------------------------------------------!

SUBROUTINE cl_dust_tendency()


IMPLICIT NONE 

END SUBROUTINE cl_dust_tendency