SUBROUTINE TLM_chem(dcc0_dpq, dccn_dpq, dHOX_dPQ, dOX_dPQ, &
                    dHOX0_dPQ, dOX0_dPQ, &
                   ptimestep, istep, phychemrat, &
                   dens, sza, initialisation)

USE TLMvars

IMPLICIT NONE 

#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "callkeys.h"
#include "conc.h"

! =====
! Input
! =====
real dcc0_dpq(nqmx*nlayermx,nqmx*nlayermx) ! Initial Linearised Number Denisty Array 
real dccn_dpq(nqmx*nlayermx,nqmx*nlayermx) ! Sub-timestep n Linearised Number Density Array 
real dHOX_dPQ(nlayermx,nqmx*nlayermx) ! Linearised HOx Number Density
real dOX_dPQ(nlayermx,nqmx*nlayermx) ! Linearised Ox Number Density
real dHOX0_dPQ(nlayermx,nqmx*nlayermx) ! Initial Linearised HOx Number Density
real dOX0_dPQ(nlayermx,nqmx*nlayermx) ! Initial Linearised Ox Number Density

REAL ptimestep ! Physical Time-step 
INTEGER istep ! Chemistry Sub-timestep iteration 
INTEGER phychemrat ! Number of Chemistry Sub-timestep iterations 

REAL dens(nlayermx) ! Atmospheric Number Density 
REAL sza ! Solar Zenith Angle (degrees)
LOGICAL initialisation ! Logical for initialisation of arrays 

! =====
! Local
! =====
real, SAVE :: dcc0_dpq_firstcall(nqmx*nlayermx,nqmx*nlayermx) 

INTEGER iq ! Tracer Iterator
INTEGER l ! Layer iterator

INTEGER h_j, oh_j, ho2_j ! TLM Indices
INTEGER o_j, o3_j ! TLM Indices 
INTEGER cl_j, clo_j !TLM Indices
INTEGER x_j ! TLM Indices



! ========================
! STAGE 0 : INITIALISATION 
! ========================
IF ( initialisation ) THEN 

    ! Initialise K_Pseudo
    dKpseudo_dPQ(:,:) = 0.0

    ! VMR --> MMR Conversion Factors 
    ! ------------------------------
    DO iq = 1, nqmx
        Avmr(:,iq) = mmean(1,:)/mmol(iq) 
    ENDDO ! of iq

    ! Linearised HOx 
    ! --------------
    DO l = 1, nlayermx 
        h_j = (t_h-1)*nlayermx + l
        oh_j = (t_oh-1)*nlayermx + l
        ho2_j = (t_ho2-1)*nlayermx + l

        dHOX_dPQ(l,:) = (TLM_ident( h_j, : )+TLM_trans( h_j, :)*ptimestep) &
               *Avmr(l,t_h)*dens(l)  &
               + (TLM_ident( oh_j, : ) + TLM_trans( oh_j, :)*ptimestep) &
               *Avmr(l,t_oh)*dens(l)  &
               + (TLM_ident( ho2_j, : ) + TLM_trans( ho2_j, :)*ptimestep) &
               *Avmr(l,t_ho2)*dens(l) 

        dHOX0_dPQ(l,:) = dHOX_dPQ(l,:)
    ENDDO ! l 

    ! Linearised Ox [Daylight only]
    ! -----------------------------
    IF ( sza .le. 95. ) THEN 
        DO l = 1, nlayermx
            o_j = (t_o-1)*nlayermx + l
            o3_j = (t_o3-1)*nlayermx + l

            dOX_dPQ(l,:) = 0.
            dOX0_dPQ(l,:) = 0.

            dOX_dPQ(l,:) =(TLM_ident( o3_j, : )+TLM_trans( o3_j, :)*ptimestep) &
                        *Avmr(l,t_o3)*dens(l) &
                        + (TLM_ident( o_j, : )+ TLM_trans( o_j, :)*ptimestep) &
                        *Avmr(l,t_o)*dens(l)

            dOX0_dPQ(l,:) = dOX_dPQ(l,:)    
        ENDDO ! l 
    ENDIF ! of sza 

    ! Linearised ClOx [if chlorine is active]
    IF (igcm_cl .ne. 0) THEN 

      DO l = 1, nlayermx 

        cl_j = (t_cl-1)*nlayermx + l 
        clo_j = (t_clo-1)*nlayermx + l 

        dClOx_dPQ(l,:) = (TLM_ident( cl_j, : )+TLM_trans( cl_j, :)*ptimestep) &
                        *Avmr(l,t_cl)*dens(l) &
                        + (TLM_ident( clo_j, : )+ TLM_trans( clo_j, :)*ptimestep) &
                        *Avmr(l,t_clo)*dens(l)

        dClOx0_dPQ(l,:) = dClOx_dPQ(l,:)

      ENDDO

    ENDIF 



    ! All Trace Gas Species
    ! ---------------------
    DO iq = 1, nqmx 
      DO l = 1,nlayermx

        x_j = ( iq - 1 )*nlayermx + l 

        ! Constant CH4 floor requires us to ignore 
        ! transport tendency for the surface
        IF ( (iq == t_ch4) .and. (l==1) ) THEN          
          dccn_dpq( x_j, : ) =  TLM_ident( x_j, :)*Avmr(l,iq)*dens(l)
        ELSE
          dccn_dpq( x_j, : ) = ( TLM_ident( x_j, :)+TLM_trans( x_j, :)*ptimestep) &
                           *Avmr(l,iq)*dens(l)
        ENDIF 

        dcc0_dpq( x_j, : ) = dccn_dpq( x_j, : )
      ENDDO 
    ENDDO ! nqmx

    !  Save the very first Linearised state 
    ! -------------------------------------
    dcc0_dpq_firstcall = dcc0_dpq 

RETURN 

ENDIF ! Initialisation

! ================================== 
! STAGE 1 : SUB-TIMESTEP ADVANCEMENT
! ==================================
dcc0_dpq = dccn_dpq 
dHOX0_dPQ = dHOX_dPQ 

IF ( sza .le. 95. ) dOX0_dPQ = dOX_dPQ 
IF ( igcm_cl .ne. 0 ) dClOx0_dPQ = dClOx_dPQ

IF ( istep < phychemrat ) RETURN 

! ===================================================
! STAGE 2 : PHOTOCHEMICAL TLM CALCULATION AND STORAGE 
! ===================================================
DO iq = 1,nqmx 
    DO l = 1, nlayermx
        x_j = (iq-1)*nlayermx + l
        TLM_photo( X_J, : ) = (dccn_dpq(X_J,:) - dcc0_dpq_firstcall(X_J,:)) &
                            /(Avmr(l,iq)*dens(L)*ptimestep)
    ENDDO 
ENDDO

TLM = TLM_trans + TLM_photo
TLM = tlm_ident + TLM*ptimestep

! ===================
! STAGE 3 : NAN CHECK 
! ===================
do iq = 1, nqmx*nlayermx 
  do l = 1, nlayermx*nqmx
       IF ( TLM(iq,l) .ne. TLM(iq,l) ) THEN 
            WRITE(*,*) "NAN AT ", iq, l 
            STOP 
       ENDIF 
  enddo 
enddo 


END SUBROUTINE