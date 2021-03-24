SUBROUTINE tlm_nox(iter, lyr_m, dens, nesp, &
                    nox, no, no2, rno2_no, &
                    cc, d001, d002, d003, &
                    d004, d005, & 
                    j_no2, &
                    dccn_dpq)

USE TLMvars

IMPLICIT NONE

#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "callkeys.h"
#include "conc.h"

! ====================================================================================
! 									Input Variables 
! ====================================================================================
integer iter ! iteration in the chimie routine
integer lyr_m ! layer we are differentiating 
real dens ! atmospheric number density 
integer  nesp ! number of tracers 

real nox ! NOx number density
real no, no2 ! NO, NO2 number density
real rno2_no ! Partition Function 

real cc(nesp) ! tracer number density
real d001, d002, d003, d004, d005 ! Chemical rate coefficients
real j_no2 ! Photolysis rate of NO2

real dccn_dpq(nqmx*nlayermx,nqmx*nlayermx)


! ====================================================================================
!                                   Output Variables 
! ====================================================================================
real drno2_no(nqmx*nlayermx) ! Linearised Partition function 
real A, B(2), C(5) ! coefficients

! Tracer indexing as in photochemistry.F 
! ======================================
integer, parameter :: i_co2  =  1
integer, parameter :: i_co   =  2
integer, parameter :: i_o    =  3
integer, parameter :: i_o1d  =  4
integer, parameter :: i_o2   =  5
integer, parameter :: i_o3   =  6
integer, parameter :: i_h    =  7
integer, parameter :: i_h2   =  8
integer, parameter :: i_oh   =  9
integer, parameter :: i_ho2  = 10
integer, parameter :: i_h2o2 = 11
!      Methane Oxidation
integer, parameter :: i_ch4  = 12
integer, parameter :: i_ch3  = 13 
integer, parameter :: i_ch3o2 = 14
integer, parameter :: i_ch3ooh = 15
integer, parameter :: i_ch3oh= 16
integer, parameter :: i_ch3o  = 17
integer, parameter :: i_hcho = 18
integer, parameter :: i_hcooh = 19
integer, parameter :: i_hoch2o2 = 20
integer, parameter :: i_hoch2oh = 21 
integer, parameter :: i_hoch2ooh = 22 
integer, parameter :: i_hco = 23
!      Alkane Oxidation
integer, parameter :: i_c2h6 = 24
integer, parameter :: i_c2h5 = 25
integer, parameter :: i_c2h5o2 = 26
integer, parameter :: i_c2h5ooh = 27
integer, parameter :: i_c2h5oh = 28
integer, parameter :: i_hoch2ch2o2 = 29
integer, parameter :: i_hoch2ch2o = 30
integer, parameter :: i_ethgly = 31
integer, parameter :: i_hyetho2h = 32
integer, parameter :: i_ch3cho = 33
integer, parameter :: i_ch3choho2 = 34
integer, parameter :: i_ch3cooh = 35
integer, parameter :: i_ch3chohooh = 36
integer, parameter :: i_ch3co = 37
integer, parameter :: i_ch3cooo = 38
integer, parameter :: i_ch3coooh = 39
integer, parameter :: i_hcoch2o2 = 40
integer, parameter :: i_glyox = 41
integer, parameter :: i_hcoco = 42
integer, parameter :: i_hooch2cho = 43
integer, parameter :: i_hoch2cho = 44
integer, parameter :: i_hochcho = 45
integer, parameter :: i_hoch2co = 46
integer, parameter :: i_hoch2co3 = 47
integer, parameter :: i_hoch2co2h = 48
integer, parameter :: i_hcoco2h = 49
integer, parameter :: i_hoch2co3h = 50
integer, parameter :: i_hcoco3h = 51
integer, parameter :: i_hcoco3 = 52
integer, parameter :: i_ch2choh = 53
!      Water, nitrogen, and Families
integer, parameter :: i_h2o  = 54
integer, parameter :: i_n2   = 55
integer, parameter :: i_hox  = 56
integer, parameter :: i_ox   = 57
integer, parameter :: i_RO2  = 58
integer, parameter :: i_dust = 59
!      Chlorine Compounds 
integer, parameter :: i_cl = 60
integer, parameter :: i_clo = 61
integer, parameter :: i_cl2 = 62
integer, parameter :: i_oclo = 63
integer, parameter :: i_cl2o2 = 64
integer, parameter :: i_hcl = 65
integer, parameter :: i_hocl = 66 
integer, parameter :: i_cloo = 67 
integer, parameter :: i_ch3ocl = 68
integer, parameter :: i_clco = 69
integer, parameter :: i_clo3 = 70
integer, parameter :: i_hclo4 = 71
integer, parameter :: i_clo4 = 72
integer, parameter :: i_clox = 73


! =======================================================
! Linearised Equations 
! =======================================================

! Equations:
!
! d[NO]/d[PQ] = -A1 * d[rno2_no]/d[PQ]
!
! d[NO2]/d[PQ] =  B1 * d[NO]/d[PQ]
!              +  B2 * d[rno2_no]/d[PQ]

A = nox/( (rno2_no+1.)**2 )

B(1) = rno2_no
B(2) = no

! 1.0 Linearised Partition Function
! ---------------------------------
!
! d[rno2_no]/d[PQ] = C1 * [O3]' 
!                  + C2 * [HO2]'
!				   + C3 * [ClO]'
!				   + C4 * [OClO]'
!                  - C3 * [O] ' 

C(1) =  d002/(j_no2 + d001*cc(i_o))
C(2) =  d003/(j_no2 + d001*cc(i_o))
C(3) = d004/(j_no2 + d001*cc(i_o))
C(4) = d005/(j_no2 + d001*cc(i_o))

C(5) = d001*rno2_no/(j_no2 + d001*cc(i_o))

drno2_no = C(1)*dccn_dpq( (t_o3-1)*nlayermx + lyr_m, : ) &
         + C(2)*dccn_dpq( (t_ho2-1)*nlayermx + lyr_m, : ) &
         + C(3)*dccn_dpq( (t_clo-1)*nlayermx + lyr_m, : ) &
         + C(4)*dccn_dpq( (t_oclo-1)*nlayermx + lyr_m, : ) &
         - C(5)*dccn_dpq( (t_o-1)*nlayermx + lyr_m, : ) 

! 1.1 Linearised NO 
! -----------------
dNO_dPQ(lyr_m,:) = -A*drno2_no

! 1.2 Linearised NO2
! ------------------
dNO2_dPQ(lyr_m,:) = B(1)*dNO_dPQ(lyr_m,:) &
                  + B(2)*drno2_no


  

RETURN


END