SUBROUTINE scalar_perturbation(q)

USE ioipsl_getincom 

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "tracer.h"
! Input 
! -----
REAL q(nlayermx,nqmx)

! Perturbation variables from callphys.def 
! ----------------------------------------
CHARACTER(len=20) tracer_name ! Perurbed tracer 
REAL scalar ! Perturbation magnitude
INTEGER l_1, l_2 ! Model layer limits 

! Local loop variables 
! --------------------
INTEGER iq ! tracer indexes
INTEGER l ! model layers 






! Extracer info on the perturbation 
! ---------------------------------
CALL getin("TRACER",tracer_name)
write(*,*) TRIM(tracer_name)

CALL getin("SCALAR",scalar)
write(*,*) scalar

CALL getin("L_1",l_1)
write(*,*) L_1 

CALL getin("L_2",l_2)
write(*,*) L_2 

! Perturbation loop 
! -----------------
DO iq = 1,nqmx

    IF ( trim(noms(iq)) == trim(tracer_name) ) THEN 

        DO l = l_1, l_2 
            q(l,iq) = q(l,iq) + q(l,iq)*scalar
        ENDDO 
        
        RETURN 
    ENDIF 

ENDDO 

WRITE(*,*) "CHECK LOOP SCRIPT FOR ERRORS."
STOP 

END 