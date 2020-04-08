      SUBROUTINE TLM_initialise(q)
      
      USE TLMvars
      
      IMPLICIT NONE


#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "conc.h"
      
! ---------------
! Input
! ---------------
      REAL q(nlayermx,nqmx)
! ---------------
! Local Variables
! ---------------
      INTEGER :: a, b, iq ! Loop iterators
      
          ALLOCATE(tlm(nlayermx*nqmx,nlayermx*nqmx))
          ALLOCATE(tlm_photo(nlayermx*nqmx,nlayermx*nqmx))

          ALLOCATE(Avmr(nlayermx,nqmx))

          ALLOCATE(tlm_ident(nlayermx*nqmx,nlayermx*nqmx))


          ! Initialise to 0
          tlm(:,:) = 0.
          tlm_photo(:,:) = 0.
          
          ! Establish the TLM identity matrix
          DO a = 1, nlayermx*nqmx 
          DO b = 1, nlayermx*nqmx
          IF (a==b) THEN
          TLM_ident(a,b) = 1.
          ELSE
          TLM_ident(a,b) = 0.
          ENDIF
          ENDDO
          ENDDO

          do iq = 1,nqmx
          if (TRIM(noms(iq)) == "co2") then
          t_co2 = iq

          elseif (TRIM(noms(iq)) == "co") then
          t_co = iq

          elseif (TRIM(noms(iq)) == "o") then
          t_o = iq

          elseif (TRIM(noms(iq)) == "o1d") then
          t_o1d = iq

          elseif (TRIM(noms(iq)) == "o2") then
          t_o2 = iq

          elseif (TRIM(noms(iq)) == "o3") then
          t_o3 = iq

          elseif (TRIM(noms(iq)) == "h") then
          t_h = iq

          elseif (TRIM(noms(iq)) == "h2") then
          t_h2 = iq

          elseif (TRIM(noms(iq)) == "oh") then
          t_oh = iq

          elseif (TRIM(noms(iq)) == "ho2") then
          t_ho2 = iq

          elseif (TRIM(noms(iq)) == "h2o2") then
          t_h2o2 = iq	

          elseif (TRIM(noms(iq)) == "h2o_vap") then
          t_h2ovap = iq		

          elseif (TRIM(noms(iq)) == "h2o_ice") then
          t_h2oice = iq			

          elseif (TRIM(noms(iq)) == "ch4") then
          t_ch4 = iq

          elseif (TRIM(noms(iq)) == "ch3") then
		  t_ch3 = iq 

          elseif (TRIM(noms(iq)) == "ch3o2") then
		  t_ch3o2 = iq 

          elseif (TRIM(noms(iq)) == "ch3ooh") then
		  t_ch3ooh = iq 
		  
          elseif (TRIM(noms(iq)) == "ch3oh") then
		  t_ch3oh = iq 
		  
          elseif (TRIM(noms(iq)) == "ch3o") then
		  t_ch3o = iq 
		  
          elseif (TRIM(noms(iq)) == "hcho") then
		  t_hcho = iq 
		 
          elseif (TRIM(noms(iq)) == "hcooh") then
		  t_hcooh = iq 
		  
          elseif (TRIM(noms(iq)) == "hoch2o2") then
		  t_hoch2o2 = iq 
		  
          elseif (TRIM(noms(iq)) == "hoch2oh") then
		  t_hoch2oh = iq 
		  
          elseif (TRIM(noms(iq)) == "hoch2ooh") then
		  t_hoch2ooh = iq 
		  
          elseif (TRIM(noms(iq)) == "hco") then
		  t_hco = iq 
		 
          elseif (TRIM(noms(iq)) == "c2h6") then
		  t_c2h6 = iq 
		  
          elseif (TRIM(noms(iq)) == "c2h5") then
		  t_c2h5 = iq 
		
          elseif (TRIM(noms(iq)) == "c2h5o2") then
		  t_c2h5o2 = iq 
		
          elseif (TRIM(noms(iq)) == "c2h5ooh") then
		  t_c2h5ooh = iq 
		
          elseif (TRIM(noms(iq)) == "c2h5oh") then
		  t_c2h5oh = iq 
		
          elseif (TRIM(noms(iq)) == "hoch2ch2o2") then
		  t_hoch2ch2o2 = iq 
		
          elseif (TRIM(noms(iq)) == "hoch2ch2o") then
		  t_hoch2ch2o = iq 
		
          elseif (TRIM(noms(iq)) == "ethgly") then
		  t_ethgly = iq 
		  
          elseif (TRIM(noms(iq)) == "hyetho2h") then
		  t_hyetho2h = iq 
		  
          elseif (TRIM(noms(iq)) == "ch3cho") then
		  t_ch3cho = iq 
		  
          elseif (TRIM(noms(iq)) == "ch2choh") then
		  t_ch2choh = iq 
		  
          elseif (TRIM(noms(iq)) == "ch3choho2") then
		  t_ch3choho2 = iq 
		  
          elseif (TRIM(noms(iq)) == "ch3cooh") then
		  t_ch3cooh = iq 
		  
          elseif (TRIM(noms(iq)) == "ch3chohooh") then
		  t_ch3chohooh = iq 
		  
          elseif (TRIM(noms(iq)) == "ch3c(o)") then
		  t_ch3co = iq 
		  
          elseif (TRIM(noms(iq)) == "ch3c(o)oo") then
		  t_ch3cooo = iq 
		  
          elseif (TRIM(noms(iq)) == "ch3c(o)ooh") then
		  t_ch3coooh = iq 
		  
          elseif (TRIM(noms(iq)) == "hcoch2o2") then
		  t_hcoch2o2 = iq 
		  
          elseif (TRIM(noms(iq)) == "glyox") then
		  t_glyox = iq 
		  
          elseif (TRIM(noms(iq)) == "hcoco") then
		  t_hcoco = iq 
		  
          elseif (TRIM(noms(iq)) == "hooch2cho") then
		  t_hooch2cho = iq 
		  
          elseif (TRIM(noms(iq)) == "hoch2cho") then
		  t_hoch2cho = iq 
		  
          elseif (TRIM(noms(iq)) == "hochcho") then
		  t_hochcho = iq 
		  
          elseif (TRIM(noms(iq)) == "hoch2co") then
		  t_hoch2co = iq 
		  
          elseif (TRIM(noms(iq)) == "hoch2co3") then
		  t_hoch2co3 = iq 
		  
          elseif (TRIM(noms(iq)) == "hoch2co2h") then
		  t_hoch2co2h = iq 
		  
          elseif (TRIM(noms(iq)) == "hcoco2h") then
		  t_hcoco2h = iq 
		  
          elseif (TRIM(noms(iq)) == "hcoco3h") then
		  t_hcoco3h = iq 
		  
          elseif (TRIM(noms(iq)) == "hcoco3") then
		  t_hcoco3 = iq 
		  
          elseif (TRIM(noms(iq)) == "hoch2co3h") then
		  t_hoch2co3h = iq 

          endif
		  
		  
		  Avmr(:,iq) =mmol(iq)/mmean(1,:)
		  
          enddo

c     Initialisation of perturbation vector 
          allocate(pertvector(nqmx*nlayermx))

          pertvector(:) = 0. 

c     Manually change for respective test.
c     Possibly make automatic? 
c     Comment lines > 218 to turn this clunky routine off.
          pertvector( (t_o2-1)*nlayermx + 10 )
     $              = q(10,5)*1.e-1




          RETURN 

          END SUBROUTINE      
