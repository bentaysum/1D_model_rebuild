! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Subroutines for the L-BFGS-B Optimisation loop that occur cc
!  during the forecast time-step, t_N. This routines are :   cc
!                                                            cc
!  1) lbfgsb_cost   : Calculates the cost function value     cc 
!  2) lbfgsb_grad   : Calculates the gradient                cc 
!                                                            cc       
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

REAL*8 FUNCTION lbfgsb_cost( PQi )

USE lbfgsb_module

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "tracer.h"
#include "conc.h" 

! Input 
! =====
REAL, intent(in) :: PQi(nlayermx,nqmx)

! Local 
! =====
INTEGER iq ! Loop iterator
INTEGER, SAVE :: call_number = 1
REAL*8 frac 
! Output 
! ======
REAL*8 COST 
REAL*8 vmr_mmr ! VMR -> MMR conversion factor 

! ----------------------------------------------------------------! 
DO iq = 1, nqmx 

     IF ( trim(noms(iq)) == "o2" ) THEN 
          vmr_mmr = (1.D0*mmol(iq))/(1.0D0*mmean(1,iq)) 
          COST =  ABS(1.D0*PQi(1,iq) - J_o2*vmr_mmr)*1.D6
          FRAC = 1.D0*PQi(1,iq)/(J_o2*vmr_mmr)
          exit
     ENDIF 

ENDDO 

lbfgsb_cost = COST 

! ------------------------------------------------------------------- ! 
!                       COST.dat 
! ------------------------------------------------------------------- ! 
IF ( call_number == 1 ) THEN ! On first call, wipe out existing file if there is one 
     OPEN( 12, FILE = "COST.dat", ACTION = "WRITE", STATUS = "REPLACE")
     
     ! WRITE THE HEADER 
     WRITE( 12 , "(A15,4A23,27A23)" ) "Iteration" , "COST", "1D/Cur.","MAX GRAD", "MIN GRAD", (TRIM(NOMS(IQ)), iq = 1,nqmx)
     WRITE( 12 , "(I15, 4E23.14, 27E23.14)" ) call_number, lbfgsb_cost, frac, MAXVAL(g_lbfgsb(:nqmx*nlayermx)), &
                                   MINVAL(g_lbfgsb(:nqmx*nlayermx)), &
                                   ( MAXVAL( X( (iq-1)*nlayermx + 1 : iq*nlayermx  ) &
                                    - LBFGSB_FIRSTGUESS( (iq-1)*nlayermx + 1 : iq*nlayermx  )), iq = 1, nqmx)
     
     CLOSE(12) 
     
     call_number = call_number + 1 
     RETURN 

ELSE ! Consequent calls, insert data 
     OPEN( 12, FILE = "COST.dat", ACTION = "write", STATUS = "old", POSITION = "append")
     WRITE( 12 , "(I15, 4E23.14, 27E23.7)" ) call_number, lbfgsb_cost, frac, MAXVAL(g_lbfgsb(:nqmx*nlayermx)), &
                                   MINVAL(g_lbfgsb(:nqmx*nlayermx)), &
                                   ( MAXVAL( X( (iq-1)*nlayermx + 1 : iq*nlayermx  ) &
                                   - LBFGSB_FIRSTGUESS( (iq-1)*nlayermx + 1 : iq*nlayermx  )), iq = 1, nqmx)
     
     CLOSE(12)
     call_number = call_number + 1      
     RETURN 

ENDIF   





RETURN 

END FUNCTION lbfgsb_cost


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE lbfgsb_grad( i ) 

USE TLMvars 
USE lbfgsb_module

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "tracer.h"
#include "conc.h" 

! Use the Adjoint equation to backtrace the sensitivity of the 1-D model's
! forecast element at the forecast time-step t_N to the previously defined
! backtrace time-step t_0. 

! Input Vars 
! ----------
INTEGER i ! 1-D model time-step 

! Local Vars 
! ----------
REAL*8, ALLOCATABLE, SAVE :: TLM_stash(:,:,:) ! Store for the TLM matrices produced  
LOGICAL, SAVE :: FIRSTCALL = .TRUE. ! Prevents the re-allocation of TLM_stash
INTEGER iter, iq, t ! Loop iterators 
REAL*8 hatJ(nqmx*nlayermx) ! Sensitivity vector
REAL*8 adj(nqmx*nlayermx,nqmx*nlayermx) ! Adjoint matrix
REAL*8 vmr_mmr ! VMR -> MMR conversion factor 




! On the first call to the gradient routine, allocate the size of TLM_stash 
IF ( FIRSTCALL ) THEN 
     ALLOCATE( TLM_stash( nqmx*nlayermx, nqmx*nlayermx, t_N - t_0 ) ) 
     
     FIRSTCALL = .FALSE. 
     
     TLM_stash( : , : , 1) = TLM 
     
     RETURN 
ENDIF 

! When i == t_0 after the firstcall here occurs 
IF ( i == t_0 ) THEN 
     
     ! Clear the stash from the previous iteration 
     TLM_stash(:,:,:) = 0.D0
     
     ! Prescribe first values 
     TLM_stash(:,:,1) = TLM 
     
     RETURN
     
! When we are between the backtrace and the forecast timesteps 
ELSEIF ( i < t_N ) THEN
     
     iter = ( i + 1 - t_0)
     
     TLM_stash(:,:,iter) = TLM 
     
     RETURN 
ENDIF 

! ==================
! Forecast Time-step 
! ==================

! Error if we reach this point and it is not the forecast timestep 
! ----------------------------------------------------------------
IF ( i .ne. t_N ) THEN
     WRITE(*,*) "TIMESTEP MISMANEGEMENT IN lbfgsb_forecast.F90"
     STOP 
ENDIF 

TLM_stash(:,:, SIZE(TLM_stash(1,1,:)) ) = TLM 

! Initialise the sensitivity vector 
! ---------------------------------
hatJ(:) = 0.D0 
! Trivial at forecast time 
! ------------------------
DO iq = 1, nqmx 
     IF ( trim(noms(iq)) == "o2" ) THEN 
          vmr_mmr = (1.D0*mmol(iq))/(1.0D0*mmean(1,iq)) 
          hatJ( (iq-1)*nlayermx + 1 ) = 1.D0 
          GOTO 100 
     ENDIF 
ENDDO 

100 DO t = (t_N-t_0) - 1 , 1, - 1
     ADJ = Transpose( TLM_stash(:,:,t) )
     hatJ = MATMUL( ADJ, hatJ )
ENDDO 

g_lbfgsb(:nqmx*nlayermx) = hatJ*1.D6


RETURN

END SUBROUTINE lbfgsb_grad
