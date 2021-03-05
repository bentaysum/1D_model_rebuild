SUBROUTINE tlm_sibem(iter, niter, lyr_m, dens, sza, &
					dt_c, dt_p, &
					nesp, cc, cc0, cc_hox_next, &
					j, production, loss, &
                         methane_enhancement, &
					a001, a002, a003, &
					b001, b002, b003, b004, b005, b006, &
					b007, b008, b009, &
					c001, c002, c003, c004, c005, c006, &
					c007, c008, c009, c010, c011, c012, &
					c013, c014, c015, c016, c017, c018, &
                    d001, d002, d003, &
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
                    ,cl052, cl053, cl054 &
                    ,cl055, cl056, cl057 &
                    ,cl058, cl059, &
					dccn_dpq, dcc0_dpq, &
                    ro2, no, no2, &
					dHOX_dPQ, dHOX0_dPQ, k_pseudo)

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
integer iter,niter ! iteration and max number of iters in the chimie routine
integer lyr_m ! layer we are differentiating 
real dens ! atmospheric number density 
real sza ! Solar zenith angle 
real dt_c, dt_p ! chemical and physical timestep
integer nesp ! number of species in the chemistry routines
real cc(nesp), cc0(nesp) ! number density of species after and before the
						 ! odd-hydrogen calculations (only H, OH and HO2 are effected)
real cc_hox_next ! cc(hox)^{t+1}
REAL dccn_dpq(nqmx*nlayermx,nqmx*nlayermx)
real ro2 ! Number density of peroxy radical compounds
real no, no2 ! NO, NO2 number density
REAL, INTENT(IN) :: dcc0_dpq(nqmx*nlayermx,nqmx*nlayermx)
REAL dHOX_dPQ(nlayermx,nqmx*nlayermx), dHOX0_dPQ(nlayermx,nlayermx*nqmx)
REAL k_pseudo 

real j(nd), production(nesp), loss(nesp) ! photolysis, production and loss values 
real methane_enhancement ! Methane enhancement factor 
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
real d001, d002, d003
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
 ,cl052, cl053, cl054 &
 ,cl055, cl056, cl057 &
 ,cl058, cl059

! ====================================================================================
!									Local Variables
! ====================================================================================
REAL dP_coeff(nqmx,nqmx)
REAL dL_coeff(nqmx,nqmx)
REAL dP_dPQ(nqmx,nqmx*nlayermx), dL_dPQ(nqmx,nqmx*nlayermx)

REAL A(nqmx,3) ! Coefficients
REAL A_hox(4) 

REAL dPhox_coeff(nqmx), dLhox_coeff(nqmx) 
REAL dPhox_dPQ(nqmx*nlayermx), dLhox_dPQ(nqmx*nlayermx)

REAL A_clox(3)
REAL dPclox_coeff(nqmx), dLclox_coeff(nqmx)
REAL dPclox_dPQ(nqmx*nlayermx), dLclox_dPQ(nqmx*nlayermx) 

integer iq_j,iq_i, iq
integer x_i, x_j

! Optional Outputs
! ================


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


! ====================================================
! 1.0 : Linearised Loss Rates 
! ====================================================
    dL_coeff(:,:) = 0. 
    dLhox_coeff(:) = 0.
! 1.1 : Inorganic Species [and methane]
! -------------------------------------
    ! 1.1.1: CO 
    ! ------------------------
    dL_coeff(t_co,t_oh) = e001  
    dL_coeff(t_co,t_o) = e002   

    ! 1.1.2: O2 
    ! ------------------------
    dL_coeff(t_o2,t_o) = a001  
    dL_coeff(t_o2,t_h) = c011   

    ! 1.1.3: H2 
    ! ------------------------
    dL_coeff(t_h2,t_o1d) =  b003  
    dL_coeff(t_h2,t_oh) = c010   

    ! 1.1.4: H2O 
    ! ------------------------
    dL_coeff(t_h2ovap,t_o1d) =  b002  

    ! 1.1.5: H2O2
    ! ------------------------
    dL_coeff(t_h2o2,t_oh) = c009  
    dL_coeff(t_h2o2,t_o) = c012  

    ! 1.1.6: CH4 
    ! ------------------------
    dL_coeff(t_ch4,t_o1d) = b007 + b008 + b009 
    dL_coeff(t_ch4,t_oh) = cab001   
    dL_coeff(t_ch4,t_o) = cab002   

    ! 1.1.7: HOx
    ! ------------------------
    dLhox_coeff(t_h) = 2.*c005*cc(i_ho2) &
                     + 2.*c006*cc(i_ho2) &
                     + 4.*c018*cc(i_h) &
                     + cab027*cc(i_hco) &
                     + cl037*cc(i_cl2) &
                     + cl041*cc(i_hcl)

    dLhox_coeff(t_oh) = 2.*c007*cc(i_ho2) &
                      + 4.*c013*cc(i_oh) &
                      + 4.*c017*cc(i_oh) &
                      + cab001*cc(i_ch4) &
                      + cab013*cc(i_ch3oh)*0.15 &
                      + cab014*0.6*cc(i_ch3ooh) &
                      + cab018*cc(i_hcho) &
                      + cab025*cc(i_hco) &
                      + cab034*cc(i_hoch2ooh) &
                      + 0.06*cl012*cc(i_clo) &
                      + cl014*cc(i_hcl) &
                      + cl015*cc(i_hocl) &
                      + cl032*cc(i_ch3ocl) &
                      + cl036*cc(i_cl2) &
                      + cl048*cc(i_clo3) &
                      + cl049*cc(i_clo3) &
                      + cl027*cc(i_hclo4) &
                      + cl058*cc(i_oclo) &
                      + cab107*cc(i_ch3)

    dLhox_coeff(t_ho2) = 2.*c005*cc(i_h) &
                       + 2.*c006*cc(i_h) &
                       + 2.*c007*cc(i_oh) &
                       + 4.*c008*cc(i_ho2) &
                       + 4.*c016*cc(i_ho2) &
                       + cab006*cc(i_ch3o2) &
                       + cab007*cc(i_ch3o2) &
                       + cab019*cc(i_hcho) &
                       + 0.6*cab029*cc(i_hoch2o2) &
                       + cl009*cc(i_cl) &
                       + cl013*cc(i_clo) 

    dLhox_coeff(t_ch4) = cab001*cc(i_oh) 


    ! 1.1.8: Organic Chemistry Section
    ! --------------------------------
    IF ( igcm_ch3 .ne. 0 ) THEN ! only check for ch3 to limit number of conditionals,
                                ! will cause errors if only part of the organic routine
                                ! is used.
        ! O2 
        dL_coeff(t_o2,t_ch3) = cab003    
        dL_coeff(t_o2,t_ch3o) = cab015
        dL_coeff(t_o2,t_hco) = cab026  

        ! HOx 
        dLhox_coeff(t_ch3) = cab107*cc(i_oh)
        dLhox_coeff(t_hco) = cab025*cc(i_oh) &
                           + cab027*cc(i_h)
        dLhox_coeff(t_ch3o2) = cab006*cc(i_ho2) &
                             + cab007*cc(i_ho2)
        dLhox_coeff(t_ch3oh) = 0.15*cab013*cc(i_oh)
        dLhox_coeff(t_hcho) = cab018*cc(i_oh) &
                            + cab019*cc(i_ho2)
        dLhox_coeff(t_ch3ooh) = 0.6*cab014*cc(i_oh)
        dLhox_coeff(t_hoch2o2) = 0.6*cab029*cc(i_ho2)
        dLhox_coeff(t_hoch2ooh) = cab034*cc(i_oh)

        ! CH3 
        dL_coeff(t_ch3,t_o2) = cab003   
        dL_coeff(t_ch3,t_o3) = cab004   
        dL_coeff(t_ch3,t_o) = cab005   
        dL_coeff(t_ch3,t_oh) = cab107  
        dL_coeff(t_ch3,t_ch3) = cab038  
        dL_coeff(t_ch3,t_hco) = cab022 + cab023 
        ! CH3O2 
        dL_coeff(t_ch3o2,t_ho2) = cab006 + cab007 
        dL_coeff(t_ch3o2,t_ch3o2) = cab008 + cab009 
        dL_coeff(t_ch3o2,t_o3) = cab010   
        dL_coeff(t_ch3o2,t_oh) = cab011   
        dL_coeff(t_ch3o2,t_o) = cab012   
        dL_coeff(t_ch3o2,t_hoch2o2) = cab008 + cab009 
        ! CH3OOH 
        dL_coeff(t_ch3ooh,t_oh) = cab014   
        ! CH3OH 
        dL_coeff(t_ch3oh,t_oh) = cab013   
        ! CH3O
        dL_coeff(t_ch3o,t_o2) = cab015  
        dL_coeff(t_ch3o,t_o) = cab017  
        dL_coeff(t_ch3o,t_o3) = cab016  
        ! HCHO 
        dL_coeff(t_hcho,t_oh) = cab018  
        dL_coeff(t_hcho,t_ho2) = cab019  
        dL_coeff(t_hcho,t_o) = cab020   
        ! HCOOH 
        dL_coeff(t_hcooh,t_oh) = cab032  
        ! HOCH2O2 
        dL_coeff(t_hoch2o2,t_ho2) = cab029   
        dL_coeff(t_hoch2o2,t_hoch2o2) = cab030 + cab031
        dL_coeff(t_hoch2o2,t_ch3o2) = cab030 + cab031 
        ! HOCH2OH 
        dL_coeff(t_hoch2oh,t_oh) = cab035  
        ! HOCH2OOH
        dL_coeff(t_hoch2ooh,t_oh) = cab033 + cab034
        ! HCO
        dL_coeff(t_hco,t_o) = cab021   
        dL_coeff(t_hco,t_hco) = cab024   
        dL_coeff(t_hco,t_oh) = cab025   
        dL_coeff(t_hco,t_o2) = cab026   
        dL_coeff(t_hco,t_h) = cab027   
        dL_coeff(t_hco,t_ch3) = cab022 + cab023 

        ! Chlorine and Organic interactions 
        IF ( igcm_cl .ne. 0 ) THEN 
            ! CH3 
            dL_coeff(t_ch3,t_cl2) = cl038 
            ! CH3O2 
            dL_coeff(t_ch3o2,t_clo) = cl019 + cl020 + cl021
            dL_coeff(t_ch3o2,t_cl) = cl022 
            ! CH3OOH 
            dL_coeff(t_ch3ooh,t_cl) = cl018
            ! HCHO
            dL_coeff(t_hcho,t_cl) = cl017 
            ! HCOOH
            dL_coeff(t_hcooh,t_cl) = cl043 
        ENDIF 

    ENDIF 

    ! 1.1.9: Chlorine Chemistry 
    ! -------------------------
    IF ( igcm_cl .ne. 0 ) THEN 
        ! CO 
        dL_coeff(t_co,t_cl) = cl023 
        ! O2 
        dL_coeff(t_o2,t_cl) = cl028 
        ! H2 
        dL_coeff(t_h2,t_cl) = cl008
        ! H2O2 
        dL_coeff(t_h2o2,t_cl) = cl011
        ! CH4 
        dL_coeff(t_ch4,t_cl) = cl016 
        ! HOx
        dLhox_coeff(t_cl) = cl009*cc(i_ho2) 
        dLhox_coeff(t_clo) = 0.06*cl012*cc(i_oh) &
                           + cl013*cc(i_ho2) 
        dLhox_coeff(t_clo3) = cl048*cc(i_oh) &
                            + cl049*cc(i_oh)
        dLhox_coeff(t_oclo) = cl058*cc(i_oh)
        dLhox_coeff(t_hcl) = cl014*cc(i_oh) &
                           + cl041*cc(i_hcl)
        dLhox_coeff(t_hocl) = cl015*cc(i_oh) 
        dLhox_coeff(t_hclo4) = cl027*cc(i_oh)
        dLhox_coeff(t_ch3ocl) = cl032*cc(i_oh)
        dLhox_coeff(t_cl2) = cl036*cc(i_oh) &
                           + cl037*cc(i_h) 

        ! Cl 
        dL_coeff(t_cl,t_o2) = cl028 
        dL_coeff(t_cl,t_o3) = cl001 + cl044
        dL_coeff(t_cl,t_h2) = cl008 
        dL_coeff(t_cl,t_ho2) = cl009 + cl010 
        dL_coeff(t_cl,t_h2o2) = cl011 
        dL_coeff(t_cl,t_ch4) = cl016 
        dL_coeff(t_cl,t_co) = cl023
        dL_coeff(t_cl,t_cloo) = cl024 + cl025 
        dL_coeff(t_cl,t_ch3ocl) = cl029 + cl031
        dL_coeff(t_cl,t_cl2o2) = cl030
        dL_coeff(t_cl,t_clo4) = cl053 
        dL_coeff(t_cl,t_hocl) = cl055 + cl056 
        dL_coeff(t_cl,t_oclo) = cl059 
        ! ClO 
        dL_coeff(t_clo,t_o) = cl002 
        dL_coeff(t_clo,t_clo) = cl003 + cl004 + cl005 + cl006 
        dL_coeff(t_clo,t_oh) = cl012 
        dL_coeff(t_clo,t_ho2) = cl013 
        dL_coeff(t_clo,t_clo3) = cl045 + cl046 + cl047 
        ! ClOx 
        dLclox_coeff(:) = 0.
        dLclox_coeff(t_co) = cl023*cc(i_cl)
        dLclox_coeff(t_o2) = cl028*cc(i_cl) 
        dLclox_coeff(t_o3) = cl044*cc(i_cl)
        dLclox_coeff(t_oh) = 0.06*cl012*cc(i_cl)
        dLclox_coeff(t_ho2) = cl009*cc(i_cl) &
                            + cl013*cc(i_clo) 
        dLclox_coeff(t_h2) = cl008*cc(i_cl) 
        dLclox_coeff(t_h2o2) = cl011*cc(i_cl)
        dLclox_coeff(t_ch4) = cl016*cc(i_cl) 

        dLclox_coeff(t_clo) = 4.*cl003*cc(i_clo) &
                            + 2.*cl005*cc(i_clo) &
                            + 4.*cl006*cc(i_clo) &
                            + cl013*cc(i_ho2) &
                            + cl019*cc(i_ch3o2) &
                            + cl020*cc(i_ch3o2) &
                            + cl021*cc(i_ch3o2) &
                            + cl045*cc(i_clo3) &
                            + cl046*cc(i_clo3) &
                            + cl047*cc(i_clo3)

        dLclox_coeff(t_clo3) = cl045*cc(i_clo) &
                             + cl046*cc(i_clo) &
                             + cl047*cc(i_clo)


        dLclox_coeff(t_cl) = cl008*cc(i_h2) &
                           + cl009*cc(i_ho2) &
                           + cl011*cc(i_h2o2) &
                           + 0.06*cl012*cc(i_oh) &
                           + cl016*cc(i_ch4) &
                           + cl017*cc(i_hcho) &
                           + cl018*cc(i_ch3ooh) &
                           + cl022*cc(i_ch3o2)*0.5 &
                           + cl023*cc(i_co) &
                           + cl025*cc(i_cloo) &
                           + cl028*cc(i_o2) &
                           + cl029*cc(i_ch3ocl) &
                           + cl030*cc(i_cl2o2) &
                           + cl031*cc(i_ch3ocl) &
                           + cl043*cc(i_hcooh) &
                           + cl044*cc(i_o3) &
                           + cl055*cc(i_hocl)

        dLclox_coeff(t_hocl) = cl055*cc(i_cl)

        dLclox_coeff(t_cl2o2) = cl030*cc(i_cl)

        dLclox_coeff(t_ch3ocl) = cl029*cc(i_cl) &
                               + cl031*cc(i_ch3ocl)

        dLclox_coeff(t_cloo) = cl025*cc(i_cl)


        ! Chlorine and Organic Interactions 
        IF ( igcm_ch3 .ne. 0 ) THEN
            ! Cl  
            dL_coeff(t_cl,t_hcho) = cl017 
            dL_coeff(t_cl,t_ch3ooh) = cl018
            dL_coeff(t_cl,t_ch3o2) = cl022 
            dL_coeff(t_cl,t_hcooh) = cl043 
            ! ClO 
            dL_coeff(t_clo,t_ch3o2) = cl019 + cl020 + cl021 
            ! ClOX 
            dLclox_coeff(t_hcho) = cl017*cc(i_cl)
            dLclox_coeff(t_ch3ooh) = cl018*cc(i_cl) 
            dLclox_coeff(t_ch3o2) = cl019*cc(i_clo) &
                                  + cl020*cc(i_clo) &
                                  + cl021*cc(i_clo) &
                                  + cl022*cc(i_cl)*0.5
            dLclox_coeff(t_hcooh) = cl043*cc(i_cl) 
        ENDIF 

    ENDIF 


! 1.2 : Organic Species
! ---------------------
    





    
RETURN 

END 