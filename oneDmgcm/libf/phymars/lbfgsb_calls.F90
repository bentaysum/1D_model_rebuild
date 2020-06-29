! ==============================================================
! L-BFGS-B related subroutines are stored here with the aims of
! finding an optimum mixing ratio state at the backtrace time-
! step, t_0, to produce the required state at the forecast time-
! step, t_N. The values, t_0 and t_N, are established in the ex-
! ternal routines held in the L-BFGS-B parent directory inside
! of the 1-D model main directory [one above the oneDmgcm fol-
! der] 
! ==============================================================



! ---------------------------------------------------------------
! Backtrace and forecast time-step extraction 
! ---------------------------------------------------------------
SUBROUTINE lbfgsb_timesteps(ndt) 

USE ioipsl_getincom 

IMPLICIT NONE 

#include "callkeys.h"

! INPUT 
! =====
integer, INTENT(INOUT) :: ndt 

write(*,*) "Aqcuiring backtrace time-step"
t_backtrace = 0 ! default value
call getin("t_backtrace",t_backtrace)
write(*,*) " t_backtrace = ", t_backtrace

write(*,*) "Aqcuiring forecast time-step"
t_forecast = 1 ! default value
call getin("t_forecast",t_forecast )
write(*,*) " t_forecast = ",t_forecast

ndt = t_forecast

RETURN

END SUBROUTINE lbfgsb_timesteps 


! ---------------------------------------------------------------
! Gradient vector calculation 
! ---------------------------------------------------------------
SUBROUTINE lbfgsb_gradient(idt)

USE TLMvars 

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "tracer.h" 
#include "callkeys.h" 

! Input Variables 
! ---------------
integer, intent(in) :: idt ! 1-D model time-step on call 

! Local Variables 
! ---------------
REAL*8, SAVE :: TRANS_MATRIX(nqmx*nlayermx,nqmx*nlayermx) ! Transition Matrix  
REAL*8 grad(nqmx*nlayermx) ! Gradient Vector 

INTEGER iq ! Tracer loop iterator 
INTEGER l ! Layer loop iterator 

! .dat File Variables 
! -------------------
CHARACTER(len=*), PARAMETER :: FILENAME = "grad.dat" 



!  idt == Backtrace Timestep 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
IF ( idt == t_backtrace ) THEN 
     TRANS_MATRIX = TRANSPOSE(TLM)
     RETURN 
ENDIF 

! Backtrace Timestep < idt < Forecast Timestep
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IF ( idt .ne. t_forecast ) THEN 
     TLM = TRANSPOSE(TLM)
     TRANS_MATRIX = MATMUL( TRANS_MATRIX, TLM )
     RETURN 
ENDIF 

! idt == Forecast Timestep
! ~~~~~~~~~~~~~~~~~~~~~~~~
IF ( idt == t_forecast ) THEN 
! 1. Initialise gradient vector 
     grad(:) = 0.D0 
     
     ! STUDYING CURIOSITY SURFACE O2 FOR NOW : MAKE INTERACTIVE AT LATER DATE 
     DO iq = 1, nqmx
          IF ( trim(noms(iq)) == "o2" ) THEN 
               grad( (iq-1)*nlayermx + 1 ) = 1.D0
               grad = MATMUL( TRANS_MATRIX, GRAD )
               GOTO 100 
          ENDIF 
     ENDDO 

! 2. Write the gradient vector to a .dat file in the oneDmgcm directory
100 OPEN( unit = 20, file = FILENAME, action = "WRITE", status = "REPLACE")
WRITE(20,"(A15,A10,A23)") "TRACER", "LAYER", "GRADIENT"  
WRITE(20,"(A48)") '************************************************'
DO iq = 1, nqmx
     DO l = 1,nlayermx
          WRITE(20,"(A15,I10, E23.15)") TRIM(NOMS(IQ)), l, &
                    grad( (iq-1)*nlayermx + l ) 
     ENDDO 
ENDDO


RETURN 

ENDIF 


END SUBROUTINE lbfgsb_gradient 


