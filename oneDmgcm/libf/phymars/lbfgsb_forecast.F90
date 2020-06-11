SUBROUTINE lbfgsb_forecast(PQi, x)  

USE lbfgsb_module 

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "conc.h" 
#include "tracer.h"
! ====================================================
! Takes the passed mixing ratio vector PQ_{i} that
! comes from the optimization routines, and 
! calculate the value of the Cost Function
! COST. The gradient of COST is then found 
! via the 1-D model adjoint matrix; these are
! both passed to the L-BFGS-B routine which 
! will construct a new guess of vector PQ_{i+1}, 
! the 1-D model will rewind, and PQ_{i+1} will be
! used in the forward model. The new COST and 
! gradient will be calculated, and the process
! repeats until an ideal PQ_{I} is found to 
! minimize the value of COST.
!
! COST = ABS( PQ( O2 )_i - [Curiosity O2] ) 
! 
! GRAD = d(PQ(O2)_i)/d(PQ) (from Adjoint calculations)
!
! ====================================================

! Input Variables
! ===============
REAL, INTENT(IN) :: PQi( nlayermx, nqmx ) ! End state of 1-D forward model 
REAL*8, INTENT(INOUT) :: x( nmax ) ! Input vector of this model run 

! Local Variables
! ===============
INTEGER iq, lyr ! tracer and layer iterators 

DO iq = 1, nqmx

     IF ( trim(noms(iq)) == "o2" ) THEN 
          call LBFGSB_COST( PQi(1,iq)*1.D0 )
     ENDIF 
     
ENDDO 

END SUBROUTINE 








SUBROUTINE LBFGSB_COST(PQi_O2)

use lbfgsb_module

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "conc.h" 


! Input 
! =====
REAL*8, INTENT(IN) :: PQi_O2
! Local 
! =====
REAL*8, PARAMETER :: mmol_o2 = 32.D0
! Ouput 
! =====
REAL*8 COST 




COST = ABS( PQi_O2 - (J_O2*mmol_o2/(mmean(1,1)*1.D0)) ) 

write(*,*) COST 
WRITE(*,*) PQi_O2*(mmean(1,1)*1.D0)/mmol_o2!, J_O2*mmol_o2/(mmean(1,1)*1.D0)
STOP


RETURN 

END 
