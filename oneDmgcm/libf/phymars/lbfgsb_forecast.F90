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

! Output 
! ======
REAL*8 COST 
REAL*8 vmr_mmr ! VMR -> MMR conversion factor 

! ----------------------------------------------------------------! 
DO iq = 1, nqmx 

     IF ( trim(noms(iq)) == "o2" ) THEN 
          vmr_mmr = (1.D0*mmol(iq))/(1.0D0*mmean(1,iq)) 
          COST = ABS( (1.D0*PQi(1,iq)) - (J_o2*vmr_mmr) ) 
          exit
     ENDIF 

ENDDO 


lbfgsb_cost = COST 

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

! INPUT VARIABLES 
! ---------------
INTEGER, INTENT(IN) :: i ! 1-D model time-step 

! LOCAL VARIABLES 
! ---------------
REAL*8, SAVE :: P(nqmx*nlayermx,nqmx*nlayermx) 
LOGICAL, SAVE :: FIRSTCALL = .TRUE. 
INTEGER iq ! Loop iterator 

! ---------------------------------------------------------------------

! On first entry ( i == t_0 ) we initialise the transition matrix P 
! as the 1st TLM and return ( gradient only calculated at forecast 
! time-step t_N 
IF ( FIRSTCALL ) THEN 
     FIRSTCALL = .False. 
     
     P = TLM 
   
     RETURN 
ENDIF 

! On following entries, we calculate this time-step's (i) transition 
! matrix P^i via the formula:
!         P^i = TLM^i x P^{i-1} 
IF ( i < t_N ) THEN 
     
     P = MATMUL( TLM, P ) 

     RETURN
! On the forecast time-step we calculate the  gradient of the cost 
! function, i.e. the sensitivity vector from the adjoint model. 
ELSEIF ( i == t_N ) THEN

     P = MATMUL( TLM, P ) 

     ! Initialise the gradient vector 
     g_lbfgsb( : nqmx*nlayermx ) = 0.D0 
     ! At the forecast time-step the sensitivity vector is trivial
     DO iq = 1, nqmx
          IF ( trim(noms(iq)) == "o2" ) THEN
               g_lbfgsb( (iq-1)*nlayermx + 1 ) = 1.D0 
               EXIT 
          ENDIF 
     ENDDO 
     
     ! Transpose the transition matrix 
     P = Transpose(P) 
     
     ! Calculation of the gradient vector
     g_lbfgsb(:nqmx*nlayermx) = MATMUL( P , g_lbfgsb(:nqmx*nlayermx) ) 
     
     ! Reset FIRSTCALL for next L-BFGS-B loop iterations 
     FIRSTCALL = .True.
     
     ! Clear P for next L-BFGS-B loop 
     P(:,:) = 0.D0 
     
     RETURN 
     
ENDIF 

! If we get to here, something's gone askew with the timesteps 
WRITE(*,*) "TIME-STEP PROBLEMS"
WRITE(*,*) "i = ", i 
WRITE(*,*) "t_0 = ", t_0 
WRITE(*,*) "t_N = ", t_N  

STOP






END SUBROUTINE lbfgsb_grad
