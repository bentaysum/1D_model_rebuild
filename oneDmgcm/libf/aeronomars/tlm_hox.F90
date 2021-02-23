SUBROUTINE tlm_hox(iter, lyr_m, dens,sza,&
					rh_ho2, roh_ho2, &
					dt_c, dt_p, &
					nesp, cc, cc_prev, ro2, &
					j, &
					a001, a002, a003, &
					b001, b002, b003, b004, b005, b006, &
					b007, b008, b009, &
					c001, c002, c003, c004, c005, c006, &
					c007, c008, c009, c010, c011, c012, &
					c013, c014, c015, c016, c017, c018, &
					e001, e002, e003, &
					cab001, cab002, cab003, cab004, cab005, &
					cab006, cab007, cab008, cab009, cab010, &
					cab011, cab012, cab013, cab014, cab015, &
					cab016, cab017, cab018, cab019, cab020, &
					cab021, cab022, cab023, cab024, cab025, &
					cab026, cab027, cab028, cab029, cab030, &
					cab031, cab032, cab033, cab034, cab035, &
					cab036, cab037, cab038, cab039, cab040, &
					cab041, cab042, cab043, cab044, cab045, &
					cab046, cab047, cab048, cab049, cab050, &
					cab051, cab052, cab053, cab054, cab055, &
					cab056, cab057, cab058, cab059, cab060, &
					cab061, cab062, cab063, cab064, cab065, &
					cab066, cab067, cab068, cab069, cab070, &
					cab071, cab072, cab073, cab074, cab075, &
					cab076, cab077, cab078, cab079, cab080, &
					cab081, cab082, cab083, cab084, cab085, &
					cab086, cab087, cab088, cab089, cab090, &
					cab091, cab092, cab093, cab094, cab095, &
					cab096, cab097, cab098, cab099, cab100, & 
					cab101, cab102, cab103, cab104, cab105, &
					cab106, cab107, &
                    cl001, cl002, cl003 &
                    ,cl004, cl005, cl006 &
                    ,cl007, cl008, cl009 &
                    ,cl010, cl011, cl012 &
                    ,cl013, cl014, cl015 &
                    ,cl016, cl017, cl018 &
                    ,cl019, cl020, cl021 &
                    ,cl022, cl023, cl024 &
                    ,cl025, cl026, cl027 &
                    ,cl028, cl029, cl030 &
                    ,cl031, cl032, cl033 &
                    ,cl034, cl035, cl036 &
                    ,cl037, cl038, cl039 &
                    ,cl040, cl041, cl042 &
                    ,cl043, cl044, cl045 &
                    ,cl046, cl047, cl048 &
                    ,cl049, cl050, cl051 &
                    ,cl052, cl053, cl054, &
                    dccn_dpq, dcc0_dpq, &
                    dHOX_dPQ, dHOX0_dPQ, &
                    k_pseudo)

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
real sza 
real rh_ho2, roh_ho2 ! partition functions for H, OH, HO2 
REAL rh_ho2_denominator, roh_ho2_denominator ! Partition function denominators
REAL rh_ho2_cabN, rh_ho2_cabD ! Additional terms from organic 
							  ! CAABA reactions to the numerator (N)
							  ! and denominator (D) of the rh_ho2 
							  ! function, which require division
							  ! via number density to be dimensionless.
REAL roh_ho2_cabN, roh_ho2_cabD ! Additional terms from organic 
							    ! CAABA reactions to the numerator (N)
							    ! and denominator (D) of the roh_ho2 
							    ! function, which require division
							    ! via number density to be dimensionless.
real dt_c, dt_p ! chemical and physical timestep
integer nesp ! number of species in the chemistry routines
real cc(nesp), cc_prev(nesp) ! number density of species after and before the
						 ! odd-hydrogen calculations (only H, OH and HO2 are effected)
real ro2 ! Peroxy-Radical Number Density
REAL dccn_dpq(nqmx*nlayermx,nqmx*nlayermx)
REAL, INTENT(IN) :: dcc0_dpq(nqmx*nlayermx,nqmx*nlayermx)
REAL dHOX_dPQ(nlayermx, nqmx*nlayermx), dHOX0_dPQ(nlayermx,nqmx*nlayermx)
REAL k_pseudo

real j(nd) ! photolysis values 

! Oxygen Reaction Rates 
real a001, a002, a003
! O(1D) Reaction Rates 
real b001, b002, b003, &
     b004, b005, b006, &
     b007, b008, b009
! Hydrogen Reaction rates 
real c001, c002, c003, &
     c004, c005, c006, &
     c007, c008, c009, &
     c010, c011, c012, &
     c013, c014, c015, &
     c016, c017, c018
! Carbon Reaction Rates 
real e001, e002, e003
! Organic Reaction Rates 
real cab001, cab002, cab003, &
     cab004, cab005, cab006, &
     cab007, cab008, cab009, &
     cab010, cab011, cab012, &
     cab013, cab014, cab015, &
     cab016, cab017, cab018, &
     cab019, cab020, cab021, &
     cab022, cab023, cab024, &
     cab025, cab026, cab027, &
     cab028, cab029, cab030, &
     cab031, cab032, cab033, &
     cab034, cab035, cab036, &
     cab037, cab038, cab039, &
     cab040, cab041, cab042, &
     cab043, cab044, cab045, &
     cab046, cab047, cab048, &
     cab049, cab050, cab051, &
     cab052, cab053, cab054, &
     cab055, cab056, cab057, &
     cab058, cab059, cab060, &
     cab061, cab062, cab063, &
     cab064, cab065, cab066, &
     cab067, cab068, cab069, &
     cab070, cab071, cab072, &
     cab073, cab074, cab075, &
     cab076, cab077, cab078, &
     cab079, cab080, cab081, &
     cab082, cab083, cab084, &
     cab085, cab086, cab087, &
     cab088, cab089, cab090, &
     cab091, cab092, cab093, &
     cab094, cab095, cab096, &
     cab097, cab098, cab099, &
     cab100, cab101, cab102, &
     cab103, cab104, cab105, &
     cab106, cab107
! Chlorine Reaction Rates 
real cl001, cl002, cl003 &
 ,cl004, cl005, cl006 &
 ,cl007, cl008, cl009 &
 ,cl010, cl011, cl012 &
 ,cl013, cl014, cl015 &
 ,cl016, cl017, cl018 &
 ,cl019, cl020, cl021 &
 ,cl022, cl023, cl024 &
 ,cl025, cl026, cl027 &
 ,cl028, cl029, cl030 &
 ,cl031, cl032, cl033 &
 ,cl034, cl035, cl036 &
 ,cl037, cl038, cl039 &
 ,cl040, cl041, cl042 &
 ,cl043, cl044, cl045 &
 ,cl046, cl047, cl048 &
 ,cl049, cl050, cl051 &
 ,cl052, cl053, cl054

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

! Photolysis indexes as in photochemistry.F
! =========================================
integer j_o2_o, j_o2_o1d, j_co2_o, j_co2_o1d, &
          j_o3_o1d, j_o3_o, j_h2o, j_hdo, j_h2o2, &
          j_ho2, j_no2, j_ch4_ch3_h, j_ch4_1ch2_h2, &
          j_ch4_3ch2_h_h, j_ch4_ch_h2_h, j_ch3o2h, &
          j_ch2o_co, j_ch2o_hco, j_ch3oh, j_c2h6, j_hcl, &
          j_hocl, j_clo, j_so2, j_so, j_h2s, j_so3, &
          j_hno3, j_hno4, &
          j_ch3cho_ch3, j_ch3cho_ch4, j_ch3cho_h, & 
          j_hoch2ooh, j_hoch2cho_hco, j_hoch2cho_co, &
          j_hoch2cho_oh, j_glyox_hco, j_glyox_hcho, &
          j_glyox_h2, j_ch3cooh, j_ch3coooh, j_ch3cocooh, &
          j_cl2, j_cloo, j_cl2o2, j_oclo

! Local Variables 
! ===============
integer h_j, oh_j, ho2_j ! Indices for H, HO2, and OH locations for initialisation
integer x_j ! Tracer indice for initialisation
integer iq ! Tracer iterations

! Stage 1 
! -------
REAL A_RH(6), A_ROH(6)
! REAL PHn, LHn, PHd, LHd
REAL Px, Lx, Pho2, Lho2 
REAL dPx_dPQ(nqmx*nlayermx), dPho2_dPQ(nqmx*nlayermx)
REAL dLx_dPQ(nqmx*nlayermx), dLho2_dPQ(nqmx*nlayermx)

REAL POHn, LOHn, POHd, LOHd

REAL B_H(4,nqmx), B_OH(4,nqmx)

REAL dPn_dPQ(nqmx*nlayermx), dPd_dPQ(nqmx*nlayermx)
REAL dLn_dPQ(nqmx*nlayermx), dLd_dPQ(nqmx*nlayermx)


REAL dN_dPQ(nqmx*nlayermx), &
        dD_dPQ(nqmx*nlayermx) 
REAL drh_ho2(nqmx*nlayermx),droh_ho2(nqmx*nlayermx) ! Linearised partition functions

! Stage 2
! -------
REAL GAMMA_H(3) ! Coefficients 
REAL dH_dPQ(nqmx*nlayermx) ! d( CC[H] )/d(PQ) 
! Stage 3
! -------
REAL GAMMA_HO2(2) ! Coefficients
REAL dHO2_dPQ(nqmx*nlayermx) ! d( CC[HO2] )/d(PQ) 
! Stage 4
! -------
REAL GAMMA_OH(2) ! Coefficients
REAL dOH_dPQ(nqmx*nlayermx) ! d( CC[HO2] )/d(PQ) 

! =================================================== !
!     numbering of photolysis rates
! =================================================== !
j_o2_o         =  1      ! o2 + hv     -> o + o
j_o2_o1d       =  2      ! o2 + hv     -> o + o(1d)
j_co2_o        =  3      ! co2 + hv    -> co + o
j_co2_o1d      =  4      ! co2 + hv    -> co + o(1d)
j_o3_o1d       =  5      ! o3 + hv     -> o2 + o(1d)
j_o3_o         =  6      ! o3 + hv     -> o2 + o
j_h2o          =  7      ! h2o + hv    -> h + oh
j_hdo          =  8      ! hdo + hv    -> d + oh
j_h2o2         =  9      ! h2o2 + hv   -> oh + ohOLV
j_ho2          =  10     ! ho2 + hv    -> oh + o
j_no2          =  11     ! no2 + hv    -> no + o
j_ch4_ch3_h    =  12     ! ch4 + hv    -> ch3 + h
j_ch4_1ch2_h2  =  13     ! ch4 + hv    -> 1ch2 + h2
j_ch4_3ch2_h_h =  14     ! ch4 + hv    -> 3ch2 + h + h
j_ch4_ch_h2_h  =  15     ! ch4 + hv    -> ch + h2 + h
j_ch3o2h       =  16     ! ch3o2h + hv -> ch3o + oh
j_ch2o_hco     =  17     ! ch2o + hv   -> h + hco
j_ch2o_co      =  18     ! ch2o + hv   -> h2 + co
j_ch3oh        =  19     ! ch3oh + hv  -> ch3o + h
j_c2h6         =  20     ! c2h6 + hv   -> products
j_hcl          =  21     ! hcl + hv    -> h + cl
j_hocl         =  22     ! hocl + hv   -> oh + cl
j_clo          =  23     ! clo + hv    -> cl + o
j_so2          =  24     ! so2 + hv    -> so + o
j_so           =  25     ! so + hv     -> s + o
j_h2s          =  26     ! h2s + hv    -> hs + s
j_so3          =  27     ! so2 + hv    -> so2 + o
j_hno3         =  28     ! hno3 + hv   -> oh + no2
j_hno4         =  29     ! hno4 + hv   -> ho2 + no2

j_ch3cho_ch3   =  30     ! ch3cho + hv -> ch3 + hco 
j_ch3cho_ch4   =  31     ! ch3cho + hv -> ch4 + co 
j_ch3cho_h     =  32     ! ch3cho + hv -> ch3co + h 
j_hoch2ooh     =  33     ! hoch2ooh + hv -> hoch2o + oh 
                     ! hoch2o + o2 -> hcooh + ho2 
                     ! = hoch2ooh + o2 + hv -> hcooh + ho2 + oh
                     ! handle via 
j_hoch2cho_hco =  34     ! hoch2cho + hv -> ch2oh + hco 
j_hoch2cho_co  =  35     ! hoch2cho + hv -> ch3oh + co
j_hoch2cho_oh  =  36     ! hoch2cho + hv -> ch2cho + oh 
                     ! ch2cho + o2   -> hcoch2o2 
                     ! hoch2cho + o2 + hv -> hcoch2o2 + oh
j_glyox_hco    =  37     ! glyoxal  + hv -> hco + hco 
j_glyox_h2     =  38     ! glyoxal  + hv -> h2 + 2co 
j_glyox_hcho   =  39     ! glyoxal  + hv -> hcho + co 
j_ch3cooh      =  40     ! ch3cooh  + hv -> ch3 + cooh 
j_ch3coooh     =  41     ! ch3c(o)ooh + hv -> ch3 + oh + co2 
j_ch3cocooh    =  42     ! ch3coco(oh) + hv -> products 

j_cl2          = 43      ! cl2 + hv -> 2.*cl
j_cloo         = 44      ! cloo + hv -> Products 
j_oclo         = 45      ! oclo + hv -> clo + o
j_cl2o2        = 46      ! cl2o2 + hv -> cl + cloo 

! ============================================ ! 
! STAGE 1: LINEARISING PARTITION FUNCTIONS 
! ============================================ ! 
!
! rh_ho2 and roh_ho2 have the form:
!
! rX_ho2 = (Px/[HO2] + Lho2)/(Pho2/[X] + Lx)
!
! where X = [H] or [OH].
!
! Linearised with respect to tracer mixing ratios:
!
! rX_ho2 ' = A1 * Lho2 '
!          + A2 * Px '
!          + A3 * [HO2] '
!          - A4 * Lx ' 
!          - A5 * Pho2 '
!          + A6 * [X] '
!
! Where:
!        A1 = ( Pho2/[X] + Lx )^-1 
!        A2 = A1/[HO2]
!        A3 = A2 * Px/[HO2]
!        
!        A4 = (Px/[HO2] + Lx) * (Pho2/[X] + Lx)^-2
!        A5 = A4/[X]
!        A6 = A5 * Pho2/[X]
! -------------------------------------------------

! ============================================ ! 
! 1.0: RH_HO2  
! ============================================ ! 

! ----------------------------------
! 1.1: Calculate Coefficients A1-A6
! ----------------------------------
!
!    Production of [H] 
!    -----------------
 Px =           b003*cc(i_o1d)*cc(i_h2) &
            +   b008*cc(i_ch4)*cc(i_o1d) &
            +   c002*cc(i_o)*cc_prev(i_oh) &
            +   c010*cc_prev(i_OH)*cc(i_h2) &
            +   e001*cc(i_co)*cc_prev(i_oh) &
            +   cab002*cc(i_o)*cc(i_ch4)*0.49 &
            +   cab004*cc(i_ch3)*cc(i_o3)*0.956 &
            +   cab005*cc(i_ch3)*cc(i_o)*0.83  &
            +   k_pseudo*cc(i_ch4) &
            +   j(j_h2o)*cc(i_h2o) &
            +   j(j_ch2o_hco)*cc(i_hcho) &
            +   j(j_ch4_ch3_h)*cc(i_ch4) &
            +   2.*j(j_ch4_3ch2_h_h)*cc(i_ch4) &
            +   j(j_ch4_ch_h2_h)*cc(i_ch4) &
            +   j(j_ch2o_hco)*cc(i_hcho) &
            +   j(j_ch3oh)*cc(i_ch3oh) &
            +   j(j_hcl)*cc(i_hcl) &
            +   cl008*cc(i_cl)*cc(i_h2) &
            +   cl039*cc(i_hcl)*cc(i_o1d)*0.22 

!   Loss of [H]
!   -----------
    Lx =  c003*cc(i_o3) &
          + c005*cc_prev(i_ho2) &
          + c006*cc_prev(i_ho2) &
          + c011*cc(i_o2) &
          + c018*cc_prev(i_h) &
          + cab027*cc(i_hco) &
          + cl037*cc(i_cl2) &
          + cl041*cc(i_hcl)


!   Production of [HO2]
!   -------------------
    Pho2 = c011*cc(i_o2)*cc_prev(i_h) &
            + c009*cc(i_h2o2)*cc_prev(i_oh)&
              + c012*cc(i_o)*cc(i_h2o2)&
              + c014*cc_prev(i_oh)*cc(i_o3)&
              + cab011*cc(i_ch3o2)*cc_prev(i_oh)&
              + cab013*cc(i_ch3oh)*cc_prev(i_oh)*0.85&
              + cab015*cc(i_ch3o)*cc(i_o2)&
              + cab026*cc(i_hco)*cc(i_o2)&
              + cab028*cc(i_hoch2o2)&
              + cab030*cc(i_hoch2o2)*ro2&
              + cab032*cc(i_hcooh)*cc_prev(i_oh)&
              + cab035*cc(i_hoch2oh)*cc_prev(i_oh)&
              + j(j_hoch2ooh)*cc(i_hoch2ooh) & 
              + cl011*cc(i_cl)*cc(i_h2o2) &
              + cl012*cc(i_clo)*cc_prev(i_oh)*0.94 &
              + cl050*cc(i_clo3)*cc_prev(i_oh) 


!    Loss of [HO2]
!    -------------
     Lho2 =     c001*cc(i_o)&
             + c004*cc_prev(i_h)&
             + c005*cc_prev(i_h)&
             + c006*cc_prev(i_h)&
             + c007*cc_prev(i_oh)&
             + c008*cc_prev(i_ho2)&
             + c015*cc(i_o3)&
             + c016*cc_prev(i_ho2)&
             + cab006*cc(i_ch3o2)&
             + cab007*cc(i_ch3o2)&
             + cab019*cc(i_hcho)&
             + cab029*cc(i_hoch2o2)*0.8 &
             + cl009*cc(i_cl) &
             + cl010*cc(i_cl) &
             + cl013*cc(i_clo)

!   Partition Function denominator
!   ------------------------------
    rh_ho2_denominator = Pho2/MAX( 1.e-30*dens, cc_prev(i_h) ) &
                       + Lx 

!   ============
!   Coefficients
!   ============
    A_RH(1) = 1./rh_ho2_denominator 
    A_RH(2) = A_RH(1)/cc_prev(i_ho2) 
    A_RH(3) = A_RH(2)*Px/cc_prev(i_ho2) 

    A_RH(4) = A_RH(1)*RH_HO2
    A_RH(5) = A_RH(4)/cc_prev(i_h)
    A_RH(6) = A_RH(5)*Pho2/cc_prev(i_H)


! -------------------------------------
! 1.2: Linearise Loss of HO2 
! -------------------------------------
        B_h(:,:) = 0.

        B_H(1,t_o) = c001  
        B_H(1,t_h) = c004 + c005 + c006 
        B_H(1,t_oh) = c007   
        B_H(1,t_o3) = c015  
        B_H(1,t_ho2) = c008 + c016 
        if ( t_ch3o2 .ne. 0 ) B_H(1,t_ch3o2) = cab006 + cab007
        if ( t_hcho .ne. 0 )  B_H(1,t_hcho) = cab019  
        if ( t_hoch2o2 .ne. 0 ) B_H(1,t_hoch2o2) = 0.8*cab029 

        IF ( igcm_cl .ne. 0 ) THEN 
            B_H(1,t_cl) = cl009 + cl010 
            B_H(1,t_clo) = cl013
        ENDIF 

! -------------------------------------
! 1.2: Linearise Production of H  
! -------------------------------------
        B_H(2,t_co) = e001*cc_prev(i_oh)

        B_H(2,t_o1d) = b003*cc(i_h2) &
                     + b008*cc(i_ch4) &
                     + cl039*cc(i_hcl)*0.22 

        B_H(2,t_o) = c002*cc_prev(i_oh) &
                   + cab002*cc(i_ch4)*0.49 &
                   + cab005*cc(i_ch3)*0.83

        B_H(2,t_o3) = cab004*cc(i_ch3)*0.956

        B_H(2,t_oh) = c002*cc(i_o) &
                    + c010*cc(i_h2) &
                    + e001*cc(i_co)

        B_H(2,t_h2ovap) = j(j_h2o)  

        B_H(2,t_h2) = b003*cc(i_o1d) &
                    + c010*cc_prev(i_oh) &
                    + cl008*cc(i_cl) 

        B_H(2,t_ch4) = b008*cc(i_o1d) & 
                     + cab002*cc(i_o)*0.49 &
                     + k_pseudo &
                     + j(j_ch4_ch3_h) &
                     + j(j_ch4_3ch2_h_h)*2. &
                     + j(j_ch4_ch_h2_h) 

        IF ( igcm_ch3 .ne. 0 ) B_H(2,t_ch3) = cab004*cc(i_o3)*0.956 &
                                            + cab005*cc(i_o)*0.83
        IF ( igcm_hcho .ne. 0 ) B_H(2,t_hcho) = j(j_ch2o_hco)  
        IF ( igcm_ch3oh .ne. 0 ) B_H(2,t_ch3oh) = j(j_ch3oh)  

        IF ( igcm_cl .ne. 0 ) THEN 
            B_H(2,t_hcl) = j(j_hcl) &
                         + cl039*cc(i_o1d)*0.22 
            B_H(2,t_cl) = cl008*cc(i_h2)
        ENDIF 


! -------------------------------------
! 1.3: Linearise Loss of H  
! -------------------------------------
        B_H(3,t_o3) = c003  
        B_H(3,t_ho2) = c005 + c006
        B_H(3,t_o2) = c011  
        B_H(3,t_h) = c018  
        IF ( igcm_hco .ne. 0 )  B_H(3,t_hco) = cab027   
        IF ( igcm_cl .ne. 0 ) THEN 
            B_H(3,t_cl2) = cl037 
            B_H(3,t_hcl) = cl041 
        ENDIF 


! -------------------------------------
! 1.3: Linearise Production of HO2   
! -------------------------------------
        B_H(4,t_o) = c012*cc(i_h2o2)

        B_H(4,t_o2) = c011*cc_prev(i_h) &
                    + cab015*cc(i_ch3o) &
                    + cab026*cc(i_hco) 

        B_H(4,t_o3) = c014*cc_prev(i_oh)

        B_H(4,t_h) = c011*cc(i_o2)

        B_H(4,t_oh) = c009*cc(i_h2o2) &
                    + c014*cc(i_o3) &
                    + cab011*cc(i_ch3o2) &
                    + cab013*cc(i_ch3oh)*0.85 &
                    + cab032*cc(i_hcooh) &
                    + cab035*cc(i_hoch2o2) &
                    + 0.94*cl012*cc(i_clo) &
                    + cl050*cc(i_clo3)

        B_H(4,t_h2o2) = c009*cc_prev(i_oh) &
                      + c012*cc(i_o) &
                      + cl011*cc(i_cl)

        IF ( igcm_ch3o2 .ne. 0 ) B_H(4,t_ch3o2) = cab011*cc_prev(i_oh) &
                                                + cab030*cc(i_hoch2o2)

        IF ( igcm_ch3oh .ne. 0 ) B_H(4,t_ch3oh) = cab013*cc_prev(i_oh)*0.85


        IF ( igcm_ch3o .ne. 0 ) B_H(4,t_ch3o) = cab015*cc(i_o2)

        IF ( igcm_hco .ne. 0 ) B_H(4,t_hco) = cab026*cc(i_o2)

        IF ( igcm_hoch2o2 .ne. 0 ) B_H(4,t_hoch2o2) = cab028 &
                                                    + cab030*ro2 &
                                                    + cab030*cc(i_hoch2o2)

        IF ( igcm_hcooh .ne. 0 ) B_H(4,t_hcooh) = cab032*cc_prev(i_oh)

        IF ( igcm_hoch2oh .ne. 0 ) B_H(4,t_hoch2oh) = cab038*cc_prev(i_oh)

        IF ( igcm_hoch2ooh .ne. 0 ) B_H(4,t_hoch2ooh) = j(j_hoch2ooh)

        IF ( igcm_cl .ne. 0 ) THEN 
            B_H(4,t_cl) = cl011*cc(i_h2o2)
            B_H(4,t_clo) = cl012*cc_prev(i_oh)*0.94 
            B_H(4,t_clo3) = cl050*cc_prev(i_oh)
        ENDIF 
! -------------------------------------
! 1.4: Construct the Linearised Partition 
!      Function d[rh_ho2]/d[PQ]
! -------------------------------------
        dPho2_dPQ(:) = 0.
        dLx_dPQ(:) = 0. 
        dPx_dPQ(:) = 0. 
        dLho2_dPQ(:) = 0.

        DO iq = 1,nqmx

            x_j = (iq-1)*nlayermx + lyr_m

            dLho2_dPQ = dLho2_dPQ + B_H(1,iq)*dccn_dpq( x_j,: )
            dPx_dPQ = dPx_dPQ + B_H(2,iq)*dccn_dpq( x_j,: )
            dLx_dPQ = dLx_dPQ + B_H(3,iq)*dccn_dpq( x_j,: )
            dPho2_dPQ = dPho2_dPQ + B_H(4,iq)*dccn_dpq( x_j,: )

        ENDDO

        drh_ho2 = A_RH(1)*dLho2_dPQ + A_RH(2)*dPx_dPQ &
                - A_RH(3)*dccn_dpq( (t_ho2-1)*nlayermx + lyr_m, : ) &
                - A_RH(4)*dLx_dPQ - A_RH(5)*dPho2_dPQ &
                + A_RH(6)*dccn_dpq( (t_h-1)*nlayermx + lyr_m, : ) 


! ============================================ ! 
!  ROH_HO2  
! ============================================ ! 

! ----------------------------------
! 1.5: Calculate Coefficients A1-A6
! ----------------------------------

!   Production of [OH]
!   ------------------
     Px =  c001*cc(i_o)*cc_prev(i_ho2) &
             +   c004*cc_prev(i_h)*2.*cc_prev(i_ho2) &
             +   c015*cc(i_o3)*cc_prev(i_ho2) &
             +   0.2*cab029*cc(i_hoch2o2)*cc_prev(i_ho2) &
             +   2.*j(j_ho2)*cc_prev(i_ho2) &
             +   b002*cc(i_o1d)*cc(i_h2o)*2. &
                 + b003*cc(i_o1d)*cc(i_h2) &
                 + b007*cc(i_o1d)*cc(i_ch4) &
                 + c003*cc(i_o3)*cc_prev(i_h) &
                 + c012*cc(i_h2o2)*cc(i_o) &
                 + cab020*cc(i_hcho)*cc(i_o) &
                 + cab021*cc(i_hco)*cc(i_o) &
                 + j(j_h2o)*cc(i_h2o) &
                 + j(j_ch3o2h)*cc(i_ch3ooh) &
                 + 2.*j(j_h2o2)*cc(i_h2o2) &
                 + j(j_hoch2ooh)*cc(i_hoch2ooh) &
                 + cl010*cc(i_cl)*cc_prev(i_ho2) &
                 + cl018*cc(i_cl)*cc(i_ch3ooh) &
                 + j(j_hocl)*cc(i_hocl) &
                 + cl039*cc(i_hcl)*cc(i_o1d)*0.66 &
                 + cl040*cc(i_hcl)*cc(i_o) &
                 + cl042*cc(i_hocl)*cc(i_o) &
                 + cl054*cc(i_clo4)*cc(i_hocl)

!   Loss of [OH]
!   -----------
    Lx =       c002*cc(i_o) &
            + c007*cc_prev(i_ho2) &
            + c009*cc(i_h2o2) &
            + c010*cc(i_h2) &
            + c013*cc_prev(i_oh) & 
            + c014*cc(i_o3) & 
            + c017*cc_prev(i_oh) &
            + e001*cc(i_co) &
            + cab001*cc(i_ch4) &
            + cab011*cc(i_ch3o2) &
            + cab013*cc(i_ch3oh) &
            + cab014*cc(i_ch3ooh) &
            + cab018*cc(i_hcho) &
            + cab025*cc(i_hco) &
            + cab032*cc(i_hcooh) &
            + cab033*cc(i_hoch2ooh) &
            + cab035*cc(i_hoch2oh) &
            + cl012*cc(i_clo) &
            + cl014*cc(i_hcl) &
            + cl015*cc(i_hocl) &
            + cl032*cc(i_ch3ocl) &
            + cl036*cc(i_cl2) &
            + cl048*cc(i_clo3) &
            + cl049*cc(i_clo3) &
            + cl050*cc(i_clo3) &
            + cl027*cc(i_hclo4) 


    roh_ho2_denominator = Pho2/MAX( 1.e-30*dens, cc_prev(i_oh) ) &
                       + Lx 

    A_ROH(1) = 1./roh_ho2_denominator 
    A_ROH(2) = A_ROH(1)/cc_prev(i_ho2) 
    A_ROH(3) = A_ROH(2)*Px/cc_prev(i_ho2) 

    A_ROH(4) = A_ROH(1)*ROH_HO2
    A_ROH(5) = A_ROH(4)/cc_prev(i_oh)
    A_ROH(6) = A_ROH(5)*Pho2/cc_prev(i_OH)


        B_OH(:,:) = 0. 

! ----------------------------------
! 1.7: Linearised Production of OH 
! ----------------------------------
!
        B_OH(2,t_o) = c001*cc_prev(i_ho2) &
                    + c012*cc(i_h2o2) &
                    + cab020*cc(i_hcho) &
                    + cab021*cc(i_hco) &
                    + cl040*cc(i_hcl) &
                    + cl042*cc(i_hocl)

        B_OH(2,t_o1d) = b002*cc(i_h2o)*2. &
                      + b003*cc(i_h2) &
                      + b007*cc(i_ch4) &
                      + 0.66*cl039*cc(i_hcl)

        B_OH(2,t_o3) = c015*cc_prev(i_ho2) &
                     + c003*cc_prev(i_h)

        B_OH(2,t_h2ovap) = b002*cc(i_o1d)*2. &
                         + j(j_h2o)  

        B_OH(2,t_h2) = b003*cc(i_o1d) 

        B_OH(2,t_h2o2) = c012*cc(i_o) &
                       + 2.*j(j_h2o2)  

        B_OH(2,t_h) = c004*cc_prev(i_ho2)*2. &
                    + c003*cc(i_o3)

        B_OH(2,t_ho2) = c001*cc(i_o) &
                      + c004*cc_prev(i_h)*2. &
                      + c015*cc(i_o3) &
                      + 0.2*cab029*cc(i_hoch2o2) &
                      + 2.*j(j_ho2)  &
                      + cl010*cc(i_cl)

        B_OH(2,t_ch4) = b007*cc(i_o1d)

        IF (t_hoch2o2 .ne. 0 ) B_OH(2,t_hoch2o2) = 0.2*cab029*cc_prev(i_ho2)

        IF (t_hcho .ne. 0 ) B_OH(2,t_hcho) = cab020*cc(i_o) 

        IF (t_hco .ne. 0 ) B_OH(2,t_hco) = cab021*cc(i_o)

        IF (t_ch3ooh .ne. 0 ) B_OH(2,t_ch3ooh) = j(j_ch3o2h)  &
                                               + cl018*cc(i_cl)

        IF (t_hoch2ooh .ne. 0 ) B_OH(2,t_hoch2ooh) = j(j_hoch2ooh)  

        IF ( igcm_cl .ne. 0 ) THEN 
            B_OH(2,t_cl) = cl010*cc_prev(i_ho2) &
                         + cl018*cc(i_ch3ooh) 

            B_OH(2,t_hocl) = j(j_hocl) &
                           + cl042*cc(i_o) &
                           + cl054*cc(i_clo4)

            B_OH(2,t_hcl) = 0.66*cl039*cc(i_o1d) &
                          + cl040*cc(i_o)

            B_OH(2,t_clo4) = cl054*cc(i_hocl)

        ENDIF 

! ----------------------------------
! 1.8: Linearised Loss of OH 
! ----------------------------------
        B_OH(3,t_co) = e001  

        B_OH(3,t_o) = c002     

        B_OH(3,t_o3) = c014  

        B_OH(3,t_ho2) = c007  

        B_OH(3,t_h2o2) = c009   

        B_OH(3,t_h2) = c010  

        B_OH(3,t_oh) = c013 + c017

        B_OH(3,t_ch4) = cab001   

        IF ( igcm_ch3o2 .ne. 0 ) B_OH(3,t_ch3o2) = cab011  
        IF ( igcm_ch3oh .ne. 0 ) B_OH(3,t_ch3oh) = cab013  
        IF ( igcm_ch3ooh .ne. 0 ) B_OH(3,t_ch3ooh) = cab014  
        IF ( igcm_hcho .ne. 0 ) B_OH(3,t_hcho) = cab018  
        IF ( igcm_hco .ne. 0 ) B_OH(3,t_hco) = cab025  
        IF ( igcm_hcooh .ne. 0 ) B_OH(3,t_hcooh) = cab032  
        IF ( igcm_hoch2ooh .ne. 0 ) B_OH(3,t_hoch2ooh) = cab033  
        IF ( igcm_hoch2oh .ne. 0 ) B_OH(3,t_hoch2oh) = cab035  

        IF ( igcm_cl .ne. 0 ) THEN
            B_OH(3,t_clo) = cl012 
            B_OH(3,t_hcl) = cl014 
            B_OH(3,t_hocl) = cl015 
            B_OH(3,t_ch3ocl) = cl032 
            B_OH(3,t_cl2) = cl036
            B_OH(3,t_clo3) = cl048 + cl049 + cl050 
            B_OH(3,t_hclo4) = cl027 
        ENDIF 
! -------------------------------------
! 1.9: Construct the Linearised Partition 
!      Function d[roh_ho2]/d[PQ]
! -------------------------------------
        dLx_dPQ(:) = 0. 
        dPx_dPQ(:) = 0. 

        DO iq = 1,nqmx

            x_j = (iq-1)*nlayermx + lyr_m

            dPx_dPQ = dPx_dPQ + B_OH(2,iq)*dccn_dpq( x_j,: )
            dLx_dPQ = dLx_dPQ + B_OH(3,iq)*dccn_dpq( x_j,: )

        ENDDO

        droh_ho2 = A_ROH(1)*dLho2_dPQ + A_ROH(2)*dPx_dPQ &
                - A_ROH(3)*dccn_dpq( (t_ho2-1)*nlayermx + lyr_m, : ) &
                - A_ROH(4)*dLx_dPQ - A_ROH(5)*dPho2_dPQ &
                + A_ROH(6)*dccn_dpq( (t_oh-1)*nlayermx + lyr_m, : ) 

! ============================================ ! 
! STAGE 2: LINEARISING HYDROGEN ATOMS  
! ============================================ ! 
!
! cc_H(PQ) = cc_HOx(PQ)
!          / (1. + (1.+roh_ho2(PQ)/rh_ho2(PQ) ) ) 
!
! d( cc_H(PQ) )/d( PQ ) 
!   = ( 1 + (1+roh_ho2)/rh_ho2)  )^-1 * d( cc_HOx(PQ) )/d( PQ )
!   - ( cc_HOx * (1 + (1 + roh_ho2)/rh_ho2 )^-2 )/rh_ho2 * d( roh_ho2(PQ) )/d( PQ )
!   + ( cc_HOx * (1 + (1 + roh_ho2)/rh_ho2 )^-2 )*(1 + roh_ho2)/(rh_ho2^2) * d( rh_ho2(PQ) )/d( PQ )

!   2.1: Equation Coefficients
!   --------------------------
    GAMMA_H(1) = 1./(1. + ((1.+roh_ho2)/rh_ho2))
    GAMMA_H(2) = (GAMMA_H(1)**2.)*cc(i_hox)/rh_ho2
    GAMMA_H(3) = GAMMA_H(2)*(1.+roh_ho2)/rh_ho2

!   2.2: Calculation of d( cc_H )/d( PQ )
!   -------------------------------------
    dH_dPQ = GAMMA_H(1)*dHOX_dPQ(lyr_m,:) &
            - GAMMA_H(2)*droh_ho2 &
            + GAMMA_H(3)*drh_ho2

! ============================================ ! 
! STAGE 3: LINEARISING 	HO2  
! ============================================ ! 
!
! cc_HO2(PQ) = cc_H(PQ)/rh_ho2 
!              where we note that cc_H(PQ) has
!              been advanced forward one sub-
!              chemistry timestep 
!
! d( cc_HO2(PQ) )/d( PQ )
!   = rh_ho2^-1 * d( cc_H(PQ)^t+1 )/d( PQ )
!   - (cc_H^t+1 * rh_ho2^-2) * d( rh_ho2 )/d( PQ )

!   3.1: Equation Coefficients
!   --------------------------
    GAMMA_HO2(1) = 1./rh_ho2
    GAMMA_HO2(2) = cc(i_h)/(rh_ho2**2)

!   3.2: Calculation of d( cc_H )/d( PQ )
!   -------------------------------------
    dHO2_dPQ = GAMMA_HO2(1)*dH_dPQ - GAMMA_HO2(2)*drh_ho2

! ============================================ ! 
! STAGE 4: LINEARISING OH
! ============================================ ! 
!
!   cc_OH(PQ) = cc_HO2(PQ) * ROH_HO2(PQ) 
!             noting that cc_HO2 has been advanced
!             forward one sub-chemistry timestep
!
!  d( cc_OH(PQ) )/d( PQ ) 
!   = ROH_HO2 * d( cc_HO2(PQ)^t+1 )/d( PQ ) 
!   + cc_HO2^t+1 * d( roh_ho2(PQ) )/d( PQ )

!   4.1: Equation Coefficients
!   --------------------------
    GAMMA_OH(1) = roh_ho2
    GAMMA_OH(2) = cc(i_ho2)

!   4.2: Calculation of d( cc_OH )/d( PQ )
!   --------------------------
    dOH_dPQ = GAMMA_OH(1)*dHO2_dPQ + GAMMA_OH(2)*droh_ho2

! ============================================ ! 
! STAGE 5: DEPOSIT VALUES
! ============================================ ! 

!   5.1: Indices in the linearised number density
!        vector 
!   -----------------------------------------
    h_j = (t_h-1)*nlayermx + lyr_m
    oh_j = (t_oh-1)*nlayermx + lyr_m
    ho2_j = (t_ho2-1)*nlayermx + lyr_m

!   5.2: Deposition in the vector
!   -----------------------------
    dccn_dpq( h_j, : ) = dH_dPQ
    dccn_dpq( oh_j, : ) = dOH_dPQ
    dccn_dpq( ho2_j, : ) = dHO2_dPQ



RETURN


END