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
     WRITE(*,"(A46)") '**********************************************'
     WRITE(*,"(2E23.15)") MAXVAL(TLM), MAXVAL(TRANS_MATRIX)
     WRITE(*,"(2E23.15)") MINVAL(TLM), MINVAL(TRANS_MATRIX) 
     WRITE(*,"(A46)") '**********************************************'
     
     IF ( ABS(MAXVAL(TLM)) .GE. 1D11 ) THEN 
          WRITE(*,*) "ALARMINGLY BIG TLM VALUE"
          CALL tlm_error
     ENDIF 
     
     RETURN 
ENDIF 


! Backtrace Timestep < idt < Forecast Timestep
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TLM = TRANSPOSE(TLM)
TRANS_MATRIX = MATMUL( TRANS_MATRIX, TLM )

WRITE(*,"(A46)") '**********************************************'
WRITE(*,"(2E23.15)") MAXVAL(TLM), MAXVAL(TRANS_MATRIX)
WRITE(*,"(2E23.15)") MINVAL(TLM), MINVAL(TRANS_MATRIX) 
WRITE(*,"(A46)") '**********************************************'

IF ( ABS(MAXVAL(TLM)) .GE. 1D11 ) THEN 
     WRITE(*,*) "ALARMINGLY BIG TLM VALUE"
     CALL tlm_error
ELSEIF ( ABS(MAXVAL(TRANS_MATRIX)) .GE. 1D11 ) THEN 
     WRITE(*,*) "ALARMINGLY BIG TRANSITION MATRIX"
     CALL trans_error(TRANS_MATRIX)
ENDIF 

DO iq = 1, nlayermx*nqmx
     DO l = 1, nlayermx*nqmx 
          IF ( TLM( iq, l ) .NE. TLM( iq, l ) ) call tlm_error 
          IF ( TRANS_MATRIX( iq, l ) .NE. TRANS_MATRIX( iq, l ) ) call trans_error( TRANS_MATRIX )
     ENDDO 
ENDDO 



IF ( idt .ne. t_forecast) RETURN 


! idt == Forecast Timestep
! ~~~~~~~~~~~~~~~~~~~~~~~~
IF ( idt == t_forecast ) THEN 
! 1. Initialise gradient vector 
     grad(:) = 0.D0 
     
     ! STUDYING CURIOSITY SURFACE O2 FOR NOW : MAKE INTERACTIVE AT LATER DATE 
     ! DO iq = 1, nqmx
          ! IF ( trim(noms(iq)) == "o2" ) THEN 
               ! grad( (iq-1)*nlayermx + 1 ) = 1.D0
               ! grad = MATMUL( TRANS_MATRIX, GRAD )
               ! GOTO 100 
          ! ENDIF 
     ! ENDDO 
     
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

CLOSE(20)
RETURN 

ENDIF 


END SUBROUTINE lbfgsb_gradient 


! ---------------------------------------------------------------
! Final mixing ratio state from 1-D model saved to a .dat
! ---------------------------------------------------------------
SUBROUTINE lbfgsb_stateoutput(q)

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "tracer.h" 

! Input Variables
! ---------------
REAL q(nlayermx,nqmx)

! Local Variables 
! ---------------
REAL, SAVE :: q_initial(nlayermx,nqmx) ! Initial state 
CHARACTER(LEN=*), PARAMETER :: FILENAME = "finalstate.dat"
INTEGER iq, l ! Loop iterators 
CHARACTER(LEN=20) tracer
LOGICAL, SAVE :: FIRSTCALL = .True.


IF ( firstcall ) THEN 
     q_initial = q 
     firstcall = .FALSE. 
     RETURN 
ENDIF 

OPEN(unit = 30, file = FILENAME, status = "REPLACE", action = "WRITE")
WRITE(30,"(A15,A10,2A15)") "TRACER", "LAYER", "INITIAL MMR", "FINAL MMR"  

DO iq = 1, nqmx     
     DO l = 1, nlayermx
          WRITE(30,"(A15,I10, 2E15.7)") TRIM(NOMS(IQ)), l, &
                      q_initial(l,iq), q(l,iq)
     ENDDO 
ENDDO 

RETURN

END SUBROUTINE


SUBROUTINE lbfgsb_stateinput( q ) 

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "tracer.h" 

! Input Variables
! ===============
REAL, INTENT(INOUT) :: q(nlayermx,nqmx)
! Local Variables 
! ===============
CHARACTER(len=*), PARAMETER :: filename = "inputstate.dat" 
INTEGER iq, lyr ! Tracer and layer loop iterators 
LOGICAL exists
CHARACTER(len=15) tracer 
INTEGER layer 

INQUIRE( FILE = FILENAME, EXIST = EXISTS) 

IF ( .NOT. EXISTS ) THEN 
     WRITE(*,*) FILENAME , " DOES NOT EXIST."
     STOP
ENDIF 

OPEN( UNIT = 70 , FILE = FILENAME, ACTION = "READ" )

DO iq = 1, nqmx
     DO lyr = 1, nlayermx
          READ(70,"(A15,I10,E15.7)") tracer, layer, q(lyr,iq)
          
          IF ( ( TRIM(TRACER) .NE. TRIM(NOMS(IQ)) ) .OR. &
               ( layer .NE. lyr ) ) THEN 
               WRITE(*,*) "ERROR IN INPUTSTATE.DAT FORMAT" 
               WRITE(*,*) "¦",TRIM(TRACER), "¦", " ¦",TRIM(NOMS(IQ)), "¦"
               WRITE(*,*) "¦",LAYER, "¦", " ¦",lyr, "¦"
               STOP
          ENDIF 
          
          q(lyr,iq) = MAX( q(lyr,iq), 1.E-31)
          
     ENDDO 
ENDDO

CLOSE(70)

return


END SUBROUTINE lbfgsb_stateinput


SUBROUTINE trans_error(trans)

USE TLMvars 

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "tracer.h" 

INTEGER iq, l 
CHARACTER(len=20) nlayermx_string, nqmx_string, nq_nl_string 
CHARACTER(len=100) FMT_HEADER, FMT_MATRIX 
REAL*8 trans(nqmx*nlayermx,nqmx*nlayermx)

WRITE( nlayermx_string, * ) nlayermx
WRITE( nqmx_string, * ) nqmx
WRITE( nq_nl_string, * ) nqmx*nlayermx

FMT_HEADER = "(A23," // ADJUSTL(nqmx_string) // "A" // ADJUSTL(nq_nl_string) // ")"
FMT_MATRIX = "(A23," // ADJUSTL(nq_nl_string) // "E23.15)" 

OPEN( UNIT = 100, FILE = "TRANSITION_ERROR.dat", ACTION = "WRITE", STATUS = "REPLACE") 

WRITE(100,ADJUSTL(FMT_HEADER)) (TRIM(noms(iq)), iq = 1, nqmx )

DO iq = 1, nqmx
     DO l = 1, nlayermx 
          WRITE(100,ADJUSTL(FMT_MATRIX)) TRIM(NOMS(IQ)), trans( (iq-1)*nlayermx + l , : ) 
     ENDDO 
ENDDO 

CLOSE(100)
WRITE(*,*) "ERROR IN TRANSITION MATRIX" 
STOP
END SUBROUTINE trans_error


SUBROUTINE tlm_error 

USE TLMvars 

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "tracer.h" 

INTEGER iq, l 
CHARACTER(len=20) nlayermx_string, nqmx_string, nq_nl_string 
CHARACTER(len=100) FMT_HEADER, FMT_MATRIX 

WRITE( nlayermx_string, * ) nlayermx
WRITE( nqmx_string, * ) nqmx
WRITE( nq_nl_string, * ) nqmx*nlayermx

FMT_HEADER = "(A23," // ADJUSTL(nqmx_string) // "A" // ADJUSTL(nq_nl_string) // ")"
FMT_MATRIX = "(A23," // ADJUSTL(nq_nl_string) // "E23.15)" 

OPEN( UNIT = 100, FILE = "TLM_ERROR.dat", ACTION = "WRITE", STATUS = "REPLACE") 

WRITE(100,ADJUSTL(FMT_HEADER)) (TRIM(noms(iq)), iq = 1, nqmx )

DO iq = 1, nqmx
     DO l = 1, nlayermx 
          WRITE(100,ADJUSTL(FMT_MATRIX)) TRIM(NOMS(IQ)), TLM( (iq-1)*nlayermx + l , : ) 
     ENDDO 
ENDDO 

CLOSE(100)
WRITE(*,*) "ERROR IN TLM" 
STOP
END SUBROUTINE tlm_error