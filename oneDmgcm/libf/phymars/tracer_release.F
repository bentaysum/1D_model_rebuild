      SUBROUTINE tracer_release(dq)

      IMPLICIT NONE

#include "dimensions.h"
#include "dimphys.h"
#include "comgeomfi.h"
#include "surfdat.h"
#include "comsoil.h"
#include "comdiurn.h"
#include "callkeys.h"
#include "comcstfi.h"
#include "planete.h"
#include "comsaison.h"
#include "control.h"
#include "dimradmars.h"
#include "comg1d.h"
#include "tracer.h"
#include "nlteparams.h"

#include "chimiedata.h"
#include "param.h"
#include "param_v4.h"
#include "conc.h"

!     INPUT
      real dq(nlayermx,nqmx) ! Tracer tendency from 1-D model

!     LOCAL       
      real dq_plume(nlayermx,nqmx) ! Tracer tendencies to be added due to plume 
      
      character(len=*), parameter :: filename = "plume.def" 
      character(len=50) dum
      	  
      character(len=50) fileline      
      
      
      integer ERROR      
      integer eq_index, un_index     
      integer iq      
      integer lyr
            
            
      character(len=10) gasnom 
      real gasconc   
      real unit_mag
            
            
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
            
            
            
C     Initialise the tendencies to 0.
      dq_plume(:,:) = 0.
      
C     Open the "plume.def" file held in main directory
      OPEN(unit=100,file = filename)

C     Locate the line where the tracer name and concentration list begins.
C     Format specifics of the plume.def file:
C         Header lines are started with ##
C         Tracer names and concentrations are a list such as:
C                   NAME1 - CONCENTRATION1 - ppt/ppb/ppm 
C                   NAME2 - CONCENTRATION2 - ppt/ppb/ppm 
C                   NAME3 - CONCENTRATION3 - ppt/ppb/ppm 
      DO 
          READ(100,'(A)') dum
C     After the header, there should be an empty line seperating the header 
C     from the list of gases. If not, stop the program and inform user of
C     current format issue.
          IF ( dum(1:2) .ne. "##" ) GOTO 10 
      ENDDO
	
   10 IF ( len(trim(dum)) .ne. 0 ) THEN
          WRITE(*,*) "ERROR IN FORMAT OF ", filename
          WRITE(*,*) "EMPTY LINE REQUIRED BETWEEN HEADER AND TRACER ",
     $               "INFORMATION."
          STOP
      ENDIF 
      
C     Now the routine is pointing at the top of the list. We iterate through
C     the remainder of the file lines.
      DO 
          READ(100, '(A)' , END=20) fileline 
          eq_index = INDEX(fileline, "=") - 1
          
c         Tracer names and concentrations will be seperated by an = symbol"
          IF( eq_index == -1 ) THEN
               WRITE(*,*) "tracer_release.F : TRACER NAMES AND ",
     $                    "CONCENTRATIONS SEPERATED BY SPACE THEN ",
     $                    "= SYMBOL IN ", filename
          STOP
          ENDIF 
          
C         Acquire the units of the perturbation
          un_index = INDEX(fileline,"pp")
          
          IF ( TRIM(fileline(un_index:)) == "ppt" ) THEN 
               unit_mag = 1.e-12
          ELSEIF  ( TRIM(fileline(un_index:)) == "ppb" ) THEN 
               unit_mag = 1.e-9
          ELSEIF  ( TRIM(fileline(un_index:)) == "ppm" )  THEN
               unit_mag = 1.e-6
          ELSE 
               write(*,*) "tracer_release.F : NO UNIT SPECIFIED IN "
     $                    , filename
               STOP
          ENDIF
          
C     Read the name of the gas, the floating point value of the pertrubation,
C     and multiply said perturbation by the relevant magnitude.     
          gasnom = TRIM( fileline(1:eq_index) )
          READ( fileline(eq_index+2:),*) gasconc
          gasconc = gasconc*unit_mag
          
          DO iq = 1,nqmx
               IF ( trim(noms(iq)) == gasnom ) THEN 
                    DO lyr = p_low, p_high 
                      dq_plume(lyr,iq) = gasconc*mmol(iq)/mmean(1,lyr)
                    ENDDO 
                    write(*,*) dq_plume(:,iq)
                    CYCLE
               ENDIF 
          ENDDO 
          
          
          
      ENDDO 
      
   20  write(*,*) "Plume parameters successfully read..."
   
      dq = dq + dq_plume/dtphys
  
      
      
      
      
      END 