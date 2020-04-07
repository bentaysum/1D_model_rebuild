MODULE TLMvars
! Matrixes
! TLM = TLMtrans + TLMphoto
REAL, DIMENSION(:,:), allocatable,save :: tlm, tlm_photo, dhox_dpq, dox_dpq
REAL, DIMENSION(:,:), allocatable,save :: dhox0_dpq, dox0_dpq
! Perturbation Vectors
!REAL, DIMENSION(:), allocatable,save :: perts
! TLM semi-identity matrix
REAL, DIMENSION(:,:), allocatable, save :: tlm_ident
! Callkeys for TLM runs
LOGICAL,SAVE :: TLM_on = .False.
LOGICAL,SAVE :: TLM_read = .True.
! Conversion factors for mmr -> vmr 
REAL, DIMENSION(:,:), allocatable, save :: Avmr
! Indexes inside the tangent linear model 
INTEGER t_co2, t_co, t_o, t_o1d, t_o2, t_o3, t_h, t_h2, t_oh, &
	  t_ho2, t_h2o2, t_ch4, t_ch3, t_ch3o2, t_ch3ooh, t_ch3oh, &
	  t_ch3o, t_hcho, t_hcooh, t_hoch2o2, t_hoch2oh, t_hoch2ooh, &
	  t_hco, t_c2h6, t_c2h5, t_c2h5o2, t_c2h5ooh, t_c2h5oh, &
	  t_hoch2ch2o2, t_hoch2ch2o, t_ethgly, t_hyetho2h, t_ch3cho, &
	  t_ch2choh, t_ch3choho2, t_ch3cooh, t_ch3chohooh, t_ch3co, &
	  t_ch3cooo, t_ch3coooh, t_hcoch2o2, t_glyox, t_hcoco, t_hooch2cho, &
	  t_hoch2cho, t_hochcho, t_hoch2co, t_hoch2co3, t_hoch2co2h, &
	  t_hcoco2h, t_hcoco3h, t_hcoco3, t_hoch2co3h, t_h2ovap, t_h2oice
	  
CONTAINS

! -----------------------------------
! Returns required TLM values for the
! linearised photochemistry routine  
! -----------------------------------

SUBROUTINE tlm_extract(iter,pqj, pqi, lyrj, lyri, tlm_val)

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "conc.h"

! Input 
INTEGER, INTENT(IN) :: iter ! subtimestep iteration
INTEGER, INTENT(IN) :: pqj, pqi ! Index values for species j and i in the 
				 ! equation d[PQ_j]m/d[PQ_i]n at layers m and n2
INTEGER, INTENT(IN) :: lyrj, lyri ! Layer m and j in the above comment's differential

! Output 
REAL, INTENT(INOUT) :: tlm_val ! Value from the Tangent Linear Model transport equations 

IF (iter > 1) THEN 
	tlm_val = 0.
	RETURN
ELSE 	
	tlm_val = TLM( (pqj-1)*nlayermx + lyrj, (pqi-1)*nlayermx + lyri ) 
ENDIF 

END 





END MODULE 