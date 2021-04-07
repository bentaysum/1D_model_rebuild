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
					dHOX_dPQ, dHOX0_dPQ)

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

REAL A_clox(4)
REAL dPclox_coeff(nqmx), dLclox_coeff(nqmx)
REAL dPclox_dPQ(nqmx*nlayermx), dLclox_dPQ(nqmx*nlayermx) 
REAL cc_clox_next 

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
dLclox_coeff(:) = 0. 

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

        ! HCl 
        dL_coeff(t_hcl,t_oh) = cl014 
        dL_coeff(t_hcl,t_o1d) = 0.88*cl039 
        dL_coeff(t_hcl,t_o) = cl040 
        dL_coeff(t_hcl,t_h) = cl041 
        ! Cl2 
        dL_coeff(t_cl2,t_o1d) = cl035 
        dL_coeff(t_cl2,t_oh) = cl036 
        dL_coeff(t_cl2,t_h) = cl037 
        ! ClOO 
        dL_coeff(t_cloo,t_cl) = cl024 + cl025 
        ! HOCl 
        dL_coeff(t_hocl,t_oh) = cl015 
        dL_coeff(t_hocl,t_o) = cl042 
        dL_coeff(t_hocl,t_clo4) = cl054 
        dL_coeff(t_hocl,t_cl) = cl055 + cl056
        ! Cl2O2 
        dL_coeff(t_cl2o2,t_cl) = cl030 
        ! OClO 
        dL_coeff(t_oclo,t_o) = cl051 + cl057 
        dL_coeff(t_oclo,t_o3) = cl052
        dL_coeff(t_oclo,t_oh) = cl058 
        dL_coeff(t_oclo,t_cl) = cl059
        ! CH3OCl 
        dL_coeff(t_ch3ocl,t_cl) = cl029 + cl031 
        dL_coeff(t_ch3ocl,t_oh) = cl032 
        ! ClO3 
        dL_coeff(t_clo3,t_clo) = cl045 + cl046 + cl047 
        dL_coeff(t_clo3,t_oh) = cl048 + cl049 + cl050
        ! HClO4 
        dL_coeff(t_hclo4,t_oh) = cl027 
        ! ClO4 
        dL_coeff(t_clo4,t_cl) = cl053 
        dL_coeff(t_clo4,t_hocl) = cl054

        ! ClOx Summed 
        dLclox_coeff(t_clo) = 4.*cl003*cc(i_clo) &
                           + 2.*cl005*cc(i_clo) &
                           + 4.*cl006*cc(i_clo) &
                         + 0.06*cl012*cc(i_oh) &
                              + cl013*cc(i_ho2) &
                              + cl019*cc(i_ch3o2) &
                              + cl020*cc(i_ch3o2) &
                              + cl021*cc(i_ch3o2) &
                              + 0.5*cl022*cc(i_ch3o2) &
                              + cl045*cc(i_clo3) &
                              + cl046*cc(i_clo3) &
                              + cl047*cc(i_clo3)


        dLclox_coeff(t_cl) = cl008*cc(i_h2) &
                          + cl009*cc(i_ho2) &
                          + cl011*cc(i_h2o2) &
                          + cl016*cc(i_ch4) &
                          + cl017*cc(i_hcho) &
                          + cl018*cc(i_ch3ooh) &
                          + cl023*cc(i_co) &
                          + cl025*cc(i_cloo) &
                          + cl028*cc(i_o2) &
                          + cl029*cc(i_ch3ocl) &
                          + cl030*cc(i_cl2o2) &
                          + cl031*cc(i_ch3ocl) &
                          + cl043*cc(i_hcooh) &
                          + cl044*cc(i_o3) &
                          + cl055*cc(i_hocl) 

        dLclox_coeff(t_clo3) = cl045*cc(i_clo) &
                            + cl046*cc(i_clo) &
                            + cl047*cc(i_clo)

        dLclox_coeff(t_cl2o2) = cl030*cc(i_cl)

        dLclox_coeff(t_cloo) = cl025*cc(i_cl)

        dLclox_coeff(t_ch3ocl) = cl029*cc(i_cl) &
                              + cl031*cc(i_cl)

        dLclox_coeff(t_hocl) = cl055*cc(i_cl)

        dLclox_coeff(t_co) = cl023*cc(i_cl)

        dLclox_coeff(t_o2) = cl028*cc(i_cl)

        dLclox_coeff(t_o3) = cl044*cc(i_cl)

        dLclox_coeff(t_oh) = 0.06*cl012*cc(i_clo)

        dLclox_coeff(t_h2) = cl008*cc(i_cl) 

        dLclox_coeff(t_ho2) = cl009*cc(i_cl) &
                           + cl013*cc(i_clo)

        dLclox_coeff(t_h2o2) = cl011*cc(i_cl)

        dLclox_coeff(t_ch4) = cl016*cc(i_cl)

        ! Chlorine and Organic Interactions 
        IF ( igcm_ch3 .ne. 0 ) THEN
            ! Cl2 
            dL_coeff(t_cl2,t_ch3) = cl038 

            ! ClOx 
            dLclox_coeff(t_hcho) = cl017*cc(i_cl)

            dLclox_coeff(t_ch3ooh) = cl018*cc(i_cl)

            dLclox_coeff(t_ch3o2) = cl019*cc(i_clo) &
                               + cl020*cc(i_clo) &
                               + cl021*cc(i_clo) &
                               + 0.5*cl022*cc(i_clo)

            dLclox_coeff(t_hcooh) = cl043*cc(i_cl)

         ENDIF 


    ENDIF 


    ! 1.1.10: Night Time Ox 
    ! ---------------------
    IF ( sza > 95. ) THEN 
        ! O 
        dL_coeff(t_o,t_o) = 2.*a002   
        dL_coeff(t_o,t_o2) = a001   
        dL_coeff(t_o,t_o3) = a003   
        dL_coeff(t_o,t_ho2) = c001   
        dL_coeff(t_o,t_oh) = c002   
        dL_coeff(t_o,t_h2o2) = c012   
        dL_coeff(t_o,t_co) = e002   
        dL_coeff(t_o,t_ch4) = cab002   
        ! O3 
        dL_coeff(t_o3,t_o) = a003   
        dL_coeff(t_o3,t_h) = c003   
        dL_coeff(t_o3,t_oh) = c014   
        dL_coeff(t_o3,t_ho2) = c015   

        ! Organic Reactions 
        IF ( igcm_ch3 .ne. 0 ) THEN 
            ! O 
            dL_coeff(t_o,t_ch3) = cab005   
            dL_coeff(t_o,t_ch3o2) = cab012   
            dL_coeff(t_o,t_ch3o) = cab017  
            dL_coeff(t_o,t_hcho) = cab020   
            dL_coeff(t_o,t_hco) = cab021  
            ! O3 
            dL_coeff(t_o3,t_ch3) = cab004  
            dL_coeff(t_o3,t_ch3o2) = cab010   
            dL_coeff(t_o3,t_ch3o) = cab016  
        ENDIF 

        ! Chlorine Reactions 
        IF ( igcm_cl .ne. 0 ) THEN 
            ! O 
            dL_coeff(t_o,t_clo) = cl002 
            dL_coeff(t_o,t_hcl) = cl040 
            dL_coeff(t_o,t_hocl) = cl042 
            dL_coeff(t_o,t_oclo) = cl051 + cl057 
            ! O3 
            dL_coeff(t_o3,t_cl) = cl001 + cl044 
            dL_coeff(t_o3,t_oclo) = cl051
        ENDIF 

    ENDIF 


! ====================================================
! 2.0 : Linearised Production Rates 
! ====================================================
dP_coeff(:,:) = 0.
dPhox_coeff(:) = 0. 
dPclox_coeff(:) = 0.
! -----------------------
! 2.1 : Inorganic Species 
! -----------------------
!   2.1.1: CO2 
    dP_coeff(t_co2,t_co) = e001*cc(i_oh) + e002*cc(i_o) 
    dP_coeff(t_co2,t_o) = e002*cc(i_co)  
    dP_coeff(t_co2,t_o2) = cab084*cc(i_hcoco) + cab092*cc(i_hoch2co)
    dP_coeff(t_co2,t_oh) = e001*cc(i_co) + cab032*cc(i_hcooh) &
                         + cab067*cc(i_ch3cooh) + 0.91*cab098*cc(i_hoch2co2h) &
                         + cab099*cc(i_hcoco2h) + cab103*cc(i_hcoco3h)
    dP_coeff(t_co2,t_ho2) = cab072*cc(i_ch3cooo) + cab095*cc(i_hoch2co3) &  
                          + cab106*cc(i_hcoco3)
!   2.1.2: CO
    dP_coeff(t_co,t_co2) = j(j_co2_o1d) + j(j_co2_o) 

    dP_coeff(t_co,t_o) = 0.17*cab005*cc(i_ch3) &
                       + cab021*cc(i_hco) 

    dP_coeff(t_co,t_o2) = cab026*cc(i_hco)

    dP_coeff(t_co,t_h) = cab027*cc(i_hco)

    dP_coeff(t_co,t_oh) = cab025*cc(i_hco)

!   2.1.3: O2 
    dP_coeff(t_o2,t_o3) = j(j_o3_o) + j(j_o3_o1d) &
                        + 2.*a003*cc(i_o) &
                        + 2.*b005*cc(i_o1d) &
                        + b006*cc(i_o1d) &
                        + c003*cc(i_h) &
                        + c014*cc(i_oh) &
                        + 2.*c015*cc(i_ho2) &
                        + cab004*cc(i_ch3) &
                        + 2.*cab010*cc(i_ch3o2) &
                        + cab016*cc(i_ch3o) &
                        + cl001*cc(i_cl) &
                        + cl052*cc(i_oclo)

    dP_coeff(t_o2,t_o) = 2.*a002*cc(i_o) &
                       + 2.*a003*cc(i_o3) & 
                       + c001*cc(i_ho2) &
                       + c002*cc(i_oh) &
                       + d001*no2 &
                       + cab012*cc(i_ch3o2) &
                       + 0.75*cab017*cc(i_ch3o) &
                       + cl002*cc(i_clo) &
                       + cl057*cc(i_oclo) 

    dP_coeff(t_o2,t_h) = c003*cc(i_o3) &
                        + c005*cc(i_ho2)

    dP_coeff(t_o2,t_oh) = c002*cc(i_o) &
                        + c007*cc(i_ho2) &
                        + c014*cc(i_o3) &
                        + 0.06*cl012*cc(i_clo) &
                        + cl058*cc(i_oclo)

    dP_coeff(t_o2,t_ho2) = c001*cc(i_o) &
                         + c005*cc(i_h) &
                         + c007*cc(i_oh) &
                         + 2.*c008*cc(i_ho2) &
                         + 2.*c015*cc(i_o3) &
                         + 2.*c016*cc(i_ho2) &
                         + cab006*cc(i_ch3o2) &
                         + cab007*cc(i_ch3o2) &
                         + 0.8*cab029*cc(i_ho2) &
                         + cl009*cc(i_cl) &
                         + cl013*cc(i_clo)

    dP_coeff(t_o2,t_o1d) = 2.*b005*cc(i_o3) &
                         + b006*cc(i_o3) 

!   2.1.4: H2 
    dP_coeff(t_h2,t_o) = 0.17*cab005*cc(i_ch3)

    dP_coeff(t_h2,t_h) = c005*cc(i_ho2) & 
                       + 2.*c008*cc(i_h) &
                       + cab027*cc(i_hco) &
                       + cl041*cc(i_hcl)

    dP_coeff(t_h2,t_ho2) = c005*cc(i_h) 

    dP_coeff(t_h2,t_o1d) = b009*cc(i_ch4)

    dP_coeff(t_h2,t_ch4) = b009*cc(i_o1d)  &
                         + j(j_ch4_1ch2_h2) + j(j_ch4_ch_h2_h) 
!   2.1.5: H2O 
    dP_coeff(t_h2ovap,t_h2) = c010*cc(i_oh)

    dP_coeff(t_h2ovap,t_h) = c006*cc(i_ho2) 

    dP_coeff(t_h2ovap,t_oh) = c007*cc(i_ho2) &
                            + c009*cc(i_h2o2) &
                            + c010*cc(i_h2) &
                            + 2.*c013*cc(i_oh) &
                            + cab001*cc(i_ch4) &
                            + cab013*cc(i_ch3oh) &
                            + cab014*cc(i_ch3ooh) &
                            + cab018*cc(i_hcho) &
                            + cab025*cc(i_hco) &
                            + cab032*cc(i_hcooh) &
                            + cab034*cc(i_hoch2ooh) &
                            + cl014*cc(i_hcl) &
                            + cl015*cc(i_hocl) &
                            + cl027*cc(i_hclo4) 

    dP_coeff(t_h2ovap,t_ho2) = c006*cc(i_h) &
                             + c007*cc(i_oh) &
                             + cab007*cc(i_ch3o2) &
                             + 0.3*cab009*cc(i_hoch2o2)

    dP_coeff(t_h2ovap,t_h2o2) = c009*cc(i_oh)

    dP_coeff(t_h2ovap,t_ch4) = cab001*cc(i_oh)
!   2.1.6: H2O2 
    dP_coeff(t_h2o2,t_oh) = c017*cc(i_oh)

    dP_coeff(t_h2o2,t_ho2) = c008*cc(i_ho2) + c016*cc(i_ho2) 
!   2.1.7: HOx 
    dPhox_coeff(t_o1d) = 2.*b002*cc(i_h2o) &
                       + 2.*b003*cc(i_h2) &
                       + b007*cc(I_ch4) &
                       + b008*cc(i_ch4) &
                       + 0.88*cl039*cc(i_hcl)

    dPhox_coeff(t_o) = 2.*c012*cc(i_h2o2) &
                     + cab002*cc(i_ch4) &
                     + cab005*cc(i_ch3) &
                     + 0.25*cab017*cc(i_ch3o) &
                     + cab020*cc(i_hcho) &
                     + cab021*cc(i_hco) &
                     + cl040*cc(i_hcl) &
                     + cl042*cc(i_hocl)

    dPhox_coeff(t_o2) = cab015*cc(i_ch3o) &
                      + cab026*cc(i_hco)

    dPhox_coeff(t_o3) = 0.956*cab004*cc(i_ch3)

    dPhox_coeff(t_h2ovap) = 2.*j(j_h2o) &
                          + 2.*b002*cc(i_o1d)

    dPhox_coeff(t_h2) = 2.*b003*cc(i_o1d) &
                      + cl008*cc(i_cl)

    dPhox_coeff(t_h2o2) = 2.*j(j_h2o2) &
                        + 2.*c012*cc(i_o) &
                        + cl011*cc(i_cl)

    dPhox_coeff(t_ch4) = cab002*cc(i_o) &
                       + b007*cc(i_o1d) &
                       + b008*cc(i_o1d) &
                       + j(j_ch4_ch3_h) &
                       + j(j_ch4_3ch2_h_h)*2. &
                       + j(j_ch4_ch_h2_h)

! ---------------------
! 2.2 : Organic Species 
! ---------------------
    IF ( igcm_ch3 .ne. 0 ) THEN 
        ! CO2 
        dP_coeff(t_co2,t_hcooh) = cab032*cc(i_oh)
        ! CO 
        dP_coeff(t_co,t_ch3) = 0.17*cab005*cc(i_o) &
                                               + 2.*cab022*cc(i_hco) 

        dP_coeff(t_co,t_hcho) = j(j_ch2o_co)

        dP_coeff(t_co,t_hco) = cab021*cc(i_o) &
                             + cab022*cc(i_ch3) &
                             + cab024*cc(i_hco) & 
                             + cab025*cc(i_oh) &
                             + cab026*cc(i_o2) &
                             + cab027*cc(i_h) 
        ! O2 
        dP_coeff(t_o2,t_ch3) = cab004*cc(i_o3)

        dP_coeff(t_o2,t_ch3o2) = cab006*cc(i_ho2) &
                                                   + cab007*cc(i_ho2) &
                                                   + 0.5*cab008*ro2 &
                                                   + 0.5*cab008*cc(i_ch3o2) &
                                                   + 0.5*cab009*ro2 &
                                                   + 0.5*cab009*cc(i_ch3o2) &
                                                   + 2.*cab010*cc(i_o3) &
                                                   + cab012*cc(i_o) &
                                                   + 0.5*cab031*cc(i_hoch2o2) &
                                                   + cl020*cc(i_clo)

        dP_coeff(t_o2,t_ch3o) = cab016*cc(i_o3) &
                                                 + 0.75*cab017*cc(i_o)

        dP_coeff(t_o2,t_hoch2o2) = 0.8*cab029*cc(i_ho2) &
                                                       + 0.5*cab008*cc(i_ch3o2) &
                                                       + 0.5*cab009*cc(i_ch3o2) &
                                                       + 0.5*cab031*cc(i_hoch2o2) &
                                                       + 0.5*cab031*ro2
        ! H2 
        dP_coeff(t_h2,t_ch3) = 0.17*cab005*cc(i_o)

        dP_coeff(t_h2,t_hcho) = j(j_ch2o_co)
        ! H2O 
        dP_coeff(t_h2ovap,t_ch3o2) = cab007*cc(i_oh)

        dP_coeff(t_h2ovap,t_ch3oh) = cab013*cc(i_oh)

        dP_coeff(t_h2ovap,t_ch3ooh) = cab015*cc(i_oh)

        dP_coeff(t_h2ovap,t_hcho) = cab018*cc(i_oh)

        dP_coeff(t_h2ovap,t_hco) = cab025*cc(i_oh)

        dP_coeff(t_h2ovap,t_hoch2o2) = 0.3*cab009*cc(i_ho2)

        dP_coeff(t_h2ovap,t_hcooh) = cab032*cc(i_oh)

        dP_coeff(t_h2ovap,t_hoch2ooh) = cab034*cc(i_oh) 

        dP_coeff(t_h2ovap,t_hoch2oh) = cab035*cc(i_oh) 
        ! CH4 
        dP_coeff(t_ch4,t_ch3) = cab022*cc(i_hco)
        dP_coeff(t_ch4,t_hco) = cab022*cc(i_ch3)

        ! HOx 
        dPhox_coeff(t_ch3) = 0.956*cab004*cc(i_o3) &
                                          + cab005*cc(i_o)     

        dPhox_coeff(t_ch3o) = cab015*cc(i_o2) &
                                            + 0.25*cab017*cc(i_o)  

        dPhox_coeff(t_hcho) = cab020*cc(i_o) &
                                            + 2.*j(j_ch2o_hco)

        dPhox_coeff(t_ch3ooh) = j(j_ch3o2h) &
                                                + cl018*cc(i_cl)

        dPhox_coeff(t_ch3oh) = j(j_ch3oh)

        dPhox_coeff(t_hco) = cab021*cc(i_o) &
                           + cab026*cc(i_o2)

        dPhox_coeff(t_hoch2o2) = cab028 &
                               + cab030*( cc(i_hoch2o2) + ro2 )

        dPhox_coeff(t_ch3o2) = cab030*cc(i_hoch2o2)

        ! 2.2.1: CH3 
        dP_coeff( t_ch3, t_o1d ) = b007*cc(i_ch4) 

        dP_coeff( t_ch3, t_o ) = 0.51*cab002*cc(i_ch4) &
                               + 0.75*cab017*cc(i_ch3o)

        dP_coeff( t_ch3, t_oh ) = cab001*cc(i_ch4)

        dP_coeff( t_ch3, t_ch4 ) = b007*cc(i_o1d) &
                                 + cab001*cc(i_oh) &
                                 + 0.51*cab002*cc(i_o) &
                                 + j(j_ch4_ch3_h) &
                                 + cl016*cc(i_cl)

        dP_coeff(t_ch3,t_ch3o) = 0.75*cab017*cc(i_o)

        ! 2.2.2: CH3O2 
        dP_coeff( t_ch3o2, t_ch3) = cab003*cc(i_o2)
        dP_coeff( t_ch3o2, t_o2) = cab003*cc(i_ch3)
        dP_coeff( t_ch3o2, t_o3) = cab016*cc(i_ch3o)
        dP_coeff( t_ch3o2, t_oh) = 0.6*cab014*cc(i_ch3ooh)
        dP_coeff( t_ch3o2, t_ch3ooh) = 0.6*cab014*cc(i_oh)
        dP_coeff( t_ch3o2, t_ch3o) = cab016*cc(i_o3) 
        ! 2.2.3: CH3OOH 
        dP_coeff(t_ch3ooh,t_ho2) = cab006*cc(i_ch3o2)
        dP_coeff(t_ch3ooh,t_ch3o2) = cab006*cc(i_ho2)
        ! 2.2.4: CH3OH 
        dP_coeff(t_ch3oh,t_ch3o2) = 0.5*cab009*(cc(i_ch3o2) + ro2) 
        dP_coeff(t_ch3oh,t_hoch2o2) = 0.5*cab009*cc(i_ch3o2) 
        dP_coeff(t_ch3oh,t_oh) =  cab107*cc(i_ch3)
        dP_coeff(t_ch3oh,t_ch3) = cab107*cc(i_oh)
        ! 2.2.5 CH3O 
        dP_coeff(t_ch3o,t_ch4) = 0.49*cab002*cc(i_o) &
                               + b008*cc(i_o1d) 

        dP_coeff(t_ch3o,t_oh) = cab011*cc(i_ch3o2) &
                              + 0.15*cab013*cc(i_ch3oh) &
                              + cl032*cc(i_ch3ocl)

        dP_coeff(t_ch3o,t_o1d) = b008*cc(i_ch4) 

        dP_coeff(t_ch3o,t_o) = 0.49*cab002*cc(i_ch4) &
                             + cab012*cc(i_ch3o2)

        dP_coeff(t_ch3o,t_o3) = 0.044*cab004*cc(i_ch3) &
                              + cab010*cc(i_ch3o2) 

        dP_coeff(t_ch3o,t_ch3) = 0.044*cab004*cc(i_o3) 

        dP_coeff(t_ch3o,t_ch3oh) = 0.15*cab013*cc(i_oh) &
                                 + j(j_ch3oh) 

        dP_coeff(t_ch3o,t_ch3o2) = cab008*(cc(i_ch3o2)+ro2) &
                                 + cab010*cc(i_o3) &
                                 + cab011*cc(i_oh) &
                                 + cab012*cc(i_o) &
                                 + cl019*cc(i_clo) &
                                 + cl021*cc(i_clo) &
                                 + 0.5*cl022*cc(i_cl)

        dP_coeff(t_ch3o,t_ch3ooh) = j(j_ch3o2h)   

        dP_coeff(t_ch3o,t_hoch2o2) = cab008*cc(i_ch3o2) 

        ! 2.2.6: HCHO 
        dP_coeff(t_hcho, t_o1d) = b009*cc(i_ch4)

        dP_coeff(t_hcho, t_o) = 0.83*cab005*cc(i_ch3) &
                              + 0.25*cab017*cc(i_ch3o)

        dP_coeff(t_hcho, t_o2) = cab015*cc(i_ch3o)

        dP_coeff(t_hcho, t_o3) = 0.956*cab005*cc(i_ch3) 

        dP_coeff(t_hcho, t_oh) = 0.85*cab013*cc(i_ch3oh) &
                               + 0.4*cab014*cc(i_ch3ooh)

        dP_coeff(t_hcho, t_ho2) = cab007*cc(i_ch3o2)

        dP_coeff(t_hcho, t_ch4) = b009*cc(i_o1d)

        dP_coeff(t_hcho, t_ch3) = 0.956*cab004*cc(i_o3) &
                                + 0.83*cab005*cc(i_o) 

        dP_coeff(t_hcho, t_ch3o2) = cab007*cc(i_ho2) &
                                  + 0.5*cab009*(cc(i_ch3o2)+ro2)

        dP_coeff(t_hcho, t_hoch2o2) = 0.5*cab009*cc(i_ch3o2) &
                                    + cab028

        dP_coeff(t_hcho, t_ch3oh) = 0.85*cab013*cc(i_oh)

        dP_coeff(t_hcho, t_ch3ooh) = 0.4*cab014*cc(i_oh) &
                                   + cl018*cc(i_cl)

        dP_coeff(t_hcho, t_ch3o) = cab015*cc(i_o2) &
                                 + 0.25*cab017*cc(i_o)

        dP_coeff(t_hcho, t_hco) = cab024*cc(i_hco)*2. 
        ! 2.2.7: HCO 
        dP_coeff(t_hco,t_oh) = cab018*cc(i_hcho) 
        dP_coeff(t_hco,t_o) = cab020*cc(i_hcho) 
        dP_coeff(t_hco,t_hcho) = cab018*cc(i_oh) + cab020*cc(i_o) &
                               + j(j_ch2o_co) + cl017*cc(i_cl)
        ! 2.2.8: HOCH2O2 
        dP_coeff(t_hoch2o2,t_hcho) = cab019*cc(i_ho2)
        dP_coeff(t_hoch2o2,t_ho2) = cab019*cc(i_hcho) 
        dP_coeff(t_hoch2o2,t_hoch2ooh) = cab033*cc(i_oh)
        dP_coeff(t_hoch2o2,t_oh) = cab033*cc(i_hoch2ooh) 
        ! 2.2.9: HCOOH 
        dP_coeff(t_hcooh,t_oh) = cab034*cc(i_hoch2ooh) &
                               + cab035*cc(i_hoch2oh) 

        dP_coeff(t_hcooh,t_ho2) = 0.5*cab029*cc(i_hoch2o2)
 
        dP_coeff(t_hcooh,t_ch3o2) = cab030*cc(i_hoch2o2) &
                                  + cab031*0.5*cc(i_hoch2o2) 

        dP_coeff(t_hcooh,t_hoch2o2) = 0.5*cab029*cc(i_ho2) &
                                    + cab030*(cc(i_hoch2o2) + ro2 ) &
                                    + cab031*0.5*(cc(i_hoch2o2) + ro2 ) 

        dP_coeff(t_hcooh,t_hoch2ooh) = cab034*cc(i_oh)

        dP_coeff(t_hcooh,t_hoch2oh) = cab035*cc(i_oh)
        ! 2.2.10: HOCH2OH 
        dP_coeff(t_hoch2oh,t_ch3o2) = 0.5*cab031*cc(i_hoch2o2) 
        dP_coeff(t_hoch2oh,t_hoch2o2) = 0.5*cab031*( cc(i_hoch2o2) + ro2 )
        ! 2.2.11 : HOCH2OOH 
        dP_coeff(t_hoch2ooh,t_ho2) = cab029*0.5*cc(i_hoch2o2)
        dP_coeff(t_hoch2ooh,t_hoch2o2) = cab029*0.5*cc(i_ho2)

        ! Chlorine and Organic Chemistry 
        IF ( igcm_cl .ne. 0 ) THEN 
            ! CH3 
            dP_coeff(t_ch3,t_cl) = cl016*cc(i_ch4) 
            ! CH3O 
            dP_coeff(t_ch3o,t_clo) = cl019*cc(i_ch3o2) &
                                   + cl021*cc(i_ch3o2) 
            dP_coeff(t_ch3o,t_cl) = 0.5*cl022*cc(i_ch3o2) &
                                  + cl029*cc(i_ch3ocl)
            dP_coeff(t_ch3o,t_ch3ocl) = cl029*cc(i_cl) &
                                      + cl032*cc(i_oh)
            ! HCHO
            dP_coeff(t_hcho,t_cl) = cl018*cc(i_ch3ooh)
            ! HCO 
            dP_coeff(t_hco,t_cl) = cl017*cc(i_hcho)
        ENDIF 


    ENDIF 


! 2.3 : Chlorine Chemistry 
    IF ( igcm_cl .ne. 0 ) THEN 
        ! CO 
        dP_coeff(t_co,t_clco) = cl033
        ! O2 
        dP_coeff(t_o2,t_cl) = cl001*cc(i_o3) &
                            + cl009*cc(i_ho2) &
                            + cl025*cc(i_cloo) 


        dP_coeff(t_o2,t_clo) = cl002*cc(i_o) &
                             + 2.*cl003*cc(i_clo) &
                             + 2.*cl004*cc(i_clo) & 
                             + 0.06*cl012*cc(i_oh) &
                             + cl013*cc(i_ho2) &
                             + cl020*cc(i_ch3o2) 
                             
        dP_coeff(t_o2,t_cloo) = cl025*cc(i_cl) &
                              + cl026

        dP_coeff(t_o2,t_oclo) = cl052*cc(i_o3) &
                              + cl057*cc(i_o) &
                              + cl058*cc(i_oh)
        ! H2 
        dP_coeff(t_h2,t_hcl) = cl041*cc(i_h)
        ! H2O 
        dP_coeff(t_h2ovap,t_hcl) = cl014*cc(i_oh) 
        dP_coeff(t_h2ovap,t_hocl) = cl015*cc(i_oh)
        dP_coeff(t_h2ovap,t_hclo4) = cl027*cc(i_oh) 
        ! HOx 
        dPhox_coeff(t_cl) = cl008*cc(i_h2) &
                          + cl011*cc(i_h2o2) &
                          + cl018*cc(i_ch3ooh) &
                          + cl055*cc(i_hocl)  
        dPhox_coeff(t_hcl) = j(j_hcl) &
                           + 0.88*cl039*cc(i_o1d) &
                           + cl040*cc(i_o) 
        dPhox_coeff(t_hocl) = j(j_hocl) &
                            + cl042*cc(i_o) &
                            + cl054*cc(i_clo4)
        dPhox_coeff(t_clo4) = cl054*cc(i_hocl) 

        dPhox_coeff(t_hocl) = cl055*cc(i_cl)

        ! 2.3.1 : HCL 
        dP_coeff(t_hcl,t_h) = cl037*cc(i_cl2)

        dP_coeff(t_hcl,t_oh) = 0.06*cl012*cc(i_clo)

        dP_coeff(t_hcl,t_ho2) = cl009*cc(i_cl)

        dP_coeff(t_hcl,t_h2) = cl008*cc(i_cl)

        dP_coeff(t_hcl,t_h2o2) = cl011*cc(i_cl)

        dP_coeff(t_hcl,t_ch4) = cl016*cc(i_cl)

        dP_coeff(t_hcl,t_cl) = cl008*cc(i_h2) &
                             + cl009*cc(i_ho2) &
                             + cl011*cc(i_h2o2) &
                             + cl016*cc(i_ch4) &
                             + cl017*cc(i_hcho) &
                             + cl018*cc(i_ch3ooh) &
                             + cl022*cc(i_ch3o2)*0.5 &
                             + cl031*cc(i_ch3ocl) &
                             + cl043*cc(i_hcooh) &
                             + cl056*cc(i_hocl)

        dP_coeff(t_hcl,t_clo) = 0.06*cl012*cc(i_oh)

        dP_coeff(t_hcl,t_hocl) = cl056*cc(i_cl)

        dP_coeff(t_hcl,t_cl2) = cl037*cc(i_h)

        dP_coeff(t_hcl,t_ch3ocl) = cl031*cc(i_cl)
        ! 2.3.2 : Cl2 
        dP_coeff(t_cl2,t_clo) = 2.*cl003*cc(i_clo)

        dP_coeff(t_cl2,t_cl) = cl025*cc(i_cloo) &
                             + cl029*cc(i_ch3ocl) &
                             + cl030*cc(i_cl2o2) &
                             + cl055*cc(i_hocl) 

        dP_coeff(t_cl2,t_hocl) = cl055*cc(i_cl)

        dP_coeff(t_cl2,t_cl2o2) = cl030*cc(i_cl)

        dP_coeff(t_cl2,t_cloo) = cl025*cc(i_cl)

        dP_coeff(t_cl2,t_ch3ocl) = cl029*cc(i_cl)
        ! 2.3.3: ClOO 
        dP_coeff(t_cloo,t_o2) = cl028*cc(i_cl) 

        dP_coeff(t_cloo,t_cl) = cl028*cc(i_o2) &
                              + cl030*cc(i_cl2o2) 

        dP_coeff(t_cloo,t_clo) = cl019*cc(i_ch3o2) &
                               + cl045*cc(i_clo3) 

        dP_coeff(t_cloo,t_clo3) = cl045*cc(i_clo)

        dP_coeff(t_cloo,t_cl2o2) = cl030*cc(i_cl) &
                                 + j(j_cl2o2)
        ! 2.3.4: HOCl 
        dP_coeff(t_hocl,t_oh) = cl032*cc(i_ch3ocl) &
                              + cl036*cc(i_cl2) &
                              + cl058*cc(i_oclo)

        dP_coeff(t_hocl,t_ho2) = cl013*cc(i_clo)

        dP_coeff(t_hocl,t_clo) = cl013*cc(i_ho2) 

        dP_coeff(t_hocl,t_cl2) = cl036*cc(i_oh)

        dP_coeff(t_hocl,t_oclo) = cl058*cc(i_oh)

        dP_coeff(t_hocl,t_ch3ocl) = cl032*cc(i_oh)
        ! 2.3.5: Cl2O2 
        dP_coeff(t_cl2o2,t_clo) = 2.*cl006*cc(i_clo) 
        ! 2.3.6: OClO 
        dP_coeff(t_oclo,t_oh) = cl050*cc(i_clo3) 

        dP_coeff(t_oclo,t_clo) = 2.*cl005*cc(i_clo) &
                               + cl021*cc(i_ch3o2) &
                               + cl045*cc(i_clo3) &
                               + 2.*cl046*cc(i_clo3) 

        dP_coeff(t_oclo,t_clo3) = cl045*cc(i_clo) &
                                + 2.*cl046*cc(i_clo) &
                                + cl050*cc(i_oh)
        ! 2.3.7: CH3OCl 
        dP_coeff(t_ch3ocl,t_clo) = cl020*cc(i_ch3o2)
        ! 2.3.8: ClCO 
        dP_coeff(t_clco,t_co) = cl023*cc(i_cl)
        dP_coeff(t_clco,t_cl) = cl023*cc(i_co)
        ! 2.3.9: ClO3 
        dP_coeff(t_clo3,t_o) = cl051*cc(i_oclo)
        dP_coeff(t_clo3,t_o3) = cl044*cc(i_cl) &
                              + cl052*cc(i_oclo)

        dP_coeff(t_clo3,t_cl) = cl044*cc(i_o3) & 
                              + cl053*cc(i_clo4) 

        dP_coeff(t_clo3,t_oclo) = cl051*cc(i_o) &
                                + cl052*cc(i_o3)

        dP_coeff(t_clo3,t_clo4) = cl053*cc(i_cl)
        ! 2.3.10: HClO4 
        dP_coeff(t_hclo4,t_oh) = cl048*cc(i_clo3) &
                               + cl049*cc(i_clo3)

        dP_coeff(t_hclo4,t_clo3) = cl048*cc(i_oh) &
                                 + cl049*cc(i_oh) 

        dP_coeff(t_hclo4,t_hocl) = cl054*cc(i_clo4) 
        dP_coeff(t_hclo4,t_clo4) = cl054*cc(i_hocl) 
        ! 2.3.11: ClO4 
        dP_coeff(t_clo4,t_oh) = cl027*cc(i_hclo4)
        dP_coeff(t_clo4,t_hclo4) = cl027*cc(i_oh)
        
        ! 2.3.12: ClOx 
        dPclox_coeff(t_cl2o2) = cl007 &
                              + j(j_cl2o2)

        dPclox_coeff(t_cl) = cl024*cc(i_cloo) &
                           + cl059*cc(i_oclo)

        dPclox_coeff(t_cloo) = cl024*cc(i_cl) &
                             + cl026 &
                             + j(j_cloo)

        dPclox_coeff(t_clco) = cl033 

        dPclox_coeff(t_oclo) = cl057*cc(i_o) &
                             + cl059*cc(i_cl) &
                             + j(j_oclo)

        dPclox_coeff(t_hcl) = cl014*cc(i_oh) &
                            + 0.88*cl039*cc(i_o1d) &
                            + cl040*cc(i_o) &
                            + cl041*cc(i_h) &
                            + j(j_hcl)

        dPclox_coeff(t_hocl) = cl015*cc(i_oh) &
                             + cl042*cc(i_o) &
                             + j(j_hocl)

        dPclox_coeff(t_cl2) = 2.*cl035*cc(i_o1d) &
                            + cl036*cc(i_oh) &
                            + cl037*cc(i_h) &
                            + cl038*cc(i_ch3) &
                            + 2.*j(j_cl2)

        dPclox_coeff(t_o) = cl040*cc(i_hcl) &
                          + cl042*cc(i_hocl) &
                          + cl057*cc(i_oclo)

        dPclox_coeff(t_o1d) = 2.*cl035*cc(i_cl2) &
                            + 0.88*cl039*cc(i_hcl)

        dPclox_coeff(t_h) = cl037*cc(I_cl2) &
                          + cl041*cc(i_hcl)

        dPclox_coeff(t_oh) = cl014*cc(i_hcl) &
                           + cl015*cc(i_hocl) &
                           + cl036*cc(i_cl2)

        ! Chlorine and Organic Chemistry 
        IF ( igcm_ch3 .ne. 0 ) THEN 
            ! HCL 
            dP_coeff(t_hcl,t_hcho) = cl017*cc(i_cl) 
            dP_coeff(t_hcl,t_ch3ooh) = cl018*cc(i_cl)
            dP_coeff(t_hcl,t_ch3o2) = cl022*cc(i_cl)*0.5
            dP_coeff(t_hcl,t_hcooh) = cl043*cc(i_hcooh)
            ! ClOO 
            dP_coeff(t_cloo,t_ch3o2) = cl019*cc(i_clo)
            ! OClO 
            dP_coeff(t_oclo,t_ch3o2) = cl021*cc(i_clo) 
            ! CH3OCl 
            dP_coeff(t_ch3ocl,t_ch3o2) = cl020*cc(i_clo)
            ! ClOX 
            dPclox_coeff(t_ch3) = cl038*cc(i_cl2)
        ENDIF 



    ENDIF 

    ! 2.4 : Night Time Ox 
    ! -------------------
    IF ( sza > 95. ) THEN 
        ! O 
        dP_coeff(t_o,t_h) = c006*cc(i_ho2)
        dP_coeff(t_o,t_ho2) = c006*cc(i_h)
        dP_coeff(t_o,t_oh) = 2.*c013*cc(i_oh)
        ! O3 
        dP_coeff(t_o3,t_o2) = a001*cc(i_o)
        dP_coeff(t_o3,t_o) = a001*cc(i_o2)
    ENDIF 


! ==============================
! 3.0 : Dot Product Calculations 
! ==============================
dP_dPQ(:,:) = 0. 
dL_dPQ(:,:) = 0. 

dPhox_dPQ(:) = 0. 
dLhox_dPQ(:) = 0. 

dLclox_dPQ(:) = 0.
dPclox_dPQ(:) = 0. 

DO iq_j = 1, nqmx
    
    ! SIBEM and Steady-State Species 
    ! ------------------------------
    DO iq_i = 1, nqmx

        x_i = (iq_i-1)*nlayermx + lyr_m 
        
        dP_dPQ(iq_j,:) = dP_dPQ(iq_j,:) + &
                        dP_coeff(iq_j,iq_i)*dccn_dpq( x_i, : )

        dL_dPQ(iq_j,:) = dL_dPQ(iq_j,:) + &
                        dL_coeff(iq_j,iq_i)*dccn_dpq( x_i, : )

    ENDDO 
    
    ! Odd-Hydrogen [HOx]
    ! ------------------------------- 
    x_j = (iq_j-1)*nlayermx + lyr_m
    dphox_dPQ = dphox_dPQ + dphox_coeff(iq_j)*dccn_dpq( x_j, : )
    dlhox_dPQ = dLhox_dPQ + dLhox_coeff(iq_j)*dccn_dpq( x_j, : )

    if ( igcm_cl .ne. 0 ) then 
      dPclox_dPQ = dPclox_dPQ + dPclox_coeff(iq_j)*dccn_dpq(x_j,:)
      dLclox_dPQ = dLclox_dPQ + dLclox_coeff(iq_j)*dccn_dpq(x_j,:)
    endif 

ENDDO

! ===========================
! 4.0 : Equation Coefficients
! ===========================
!
! SIBEM : cc^t+1 = ( cc0 + P*dt )/( 1 + L*dt )
!
!  cc^t+1 ' = A1 * cc0 '
!           + A2 * P'
!           - A3 * L'
!
!  A1 = 1/(1 + L*dt)
!  A2 = dt/(1 + L*dt)
!  A3 = dt*(cc0 + P*dt)/(1 + L*dt)^2 
!   
! STEADY-STATE : cc^t+1 = P/L 
!
! cc^t+1 = A1 * P'
!        - A2 * L' 
!
! A1 = 1/L 
! A2 = P/L^2 
! ----------------------------
A(:,:) = 0. 

! 4.1 : Inorganics 
! ----------------
A(t_co2,1) = 1./( 1. + loss(i_co2)*dt_c)
A(t_co,1) = 1./( 1. + loss(i_co)*dt_c)
A(t_o2,1) = 1./( 1. + loss(i_o2)*dt_c)
A(t_h2ovap,1) = 1./( 1. + loss(i_h2o)*dt_c)
A(t_h2,1) = 1./( 1. + loss(i_h2)*dt_c)
A(t_h2o2,1) = 1./( 1. + loss(i_h2o2)*dt_c)
A(t_ch4,1) = 1./( 1. + loss(i_ch4)*dt_c)

A(t_co2,2) = A(t_co2,1)*dt_c 
A(t_co,2) = A(t_co,1)*dt_c 
A(t_o2,2) = A(t_o2,1)*dt_c 
A(t_h2ovap,2) = A(t_h2ovap,1)*dt_c 
A(t_h2,2) = A(t_h2,1)*dt_c 
A(t_h2o2,2) = A(t_h2o2,1)*dt_c 
A(t_ch4,2) = A(t_ch4,1)*dt_c

A(t_co2,3) = (cc0(i_co2) + production(i_co2)*dt_c)*(A(t_co2,1)**2)*dt_c
A(t_co,3) = (cc0(i_co) + production(i_co)*dt_c)*(A(t_co,1)**2)*dt_c
A(t_o2,3) = (cc0(i_o2) + production(i_o2)*dt_c)*(A(t_o2,1)**2)*dt_c
A(t_h2,3) = (cc0(i_h2) + production(i_h2)*dt_c)*(A(t_h2,1)**2)*dt_c
A(t_h2o2,3) = (cc0(i_h2o2) + production(i_h2o2)*dt_c)*(A(t_h2o2,1)**2)*dt_c
A(t_h2ovap,3) = (cc0(i_h2o) + production(i_h2o)*dt_c)*(A(t_h2ovap,1)**2)*dt_c
A(t_ch4,3) = (cc0(i_ch4) + production(i_ch4)*dt_c)*(A(t_ch4,1)**2)*dt_c

IF ( sza > 95. ) THEN 
    ! 4.1.1: Night Time O and O3
    ! --------------------------
    A(t_o,1) = 1./( 1. + loss(i_o)*dt_c)
    A(t_o3,1) = 1./( 1. + loss(i_o3)*dt_c)

    A(t_o,2) = A(t_o,1)*dt_c 
    A(t_o3,2) = A(t_o3,1)*dt_c 

    A(t_o,3) = (cc0(i_o) + production(i_o)*dt_c)*(A(t_o,1)**2)*dt_c
    A(t_o3,3) = (cc0(i_o3) + production(i_o3)*dt_c)*(A(t_o3,1)**2)*dt_c
ENDIF 

! 4.2 : Organics
! ======================
! 4.2.1 : Steady-States
! ---------------------
IF ( igcm_ch3 .ne. 0 ) THEN 
    A(t_ch3,1) = 1./loss(i_ch3)
    A(t_ch3o2,1) = 1./loss(i_ch3o2)
    A(t_ch3o,1) = 1./loss(i_ch3o)
    A(t_hoch2o2,1) = 1./loss(i_hoch2o2)
    A(t_hco,1) = 1./loss(i_hco)

    A(t_ch3,2) = production(i_ch3)*(A(t_ch3,1)**2.)
    A(t_ch3o2,2) = production(i_ch3o2)*(A(t_ch3o2,1)**2.)
    A(t_ch3o,2) = production(i_ch3o)*(A(t_ch3o,1)**2.)
    A(t_hoch2o2,2) = production(i_hoch2o2)*(A(t_hoch2o2,1)**2.)
    A(t_hco,2) = production(i_hco)*(A(t_hco,1)**2.)
! 4.2.2 : SIBEM Coefficients
! --------------------------
    A(t_ch3ooh,1) = 1./( 1. + loss(t_ch3ooh)*dt_c)
    A(t_ch3oh,1) = 1./( 1. + loss(t_ch3oh)*dt_c)
    A(t_hcho,1) = 1./( 1. + loss(t_hcho)*dt_c)
    A(t_hcooh,1) = 1./( 1. + loss(t_hcooh)*dt_c)
    A(t_hoch2oh,1) = 1./( 1. + loss(t_hoch2oh)*dt_c)
    A(t_hoch2ooh,1) = 1./( 1. + loss(t_hoch2ooh)*dt_c)

    A(t_ch3ooh,2) = A(t_ch3ooh,1)*dt_c
    A(t_ch3oh,2) = A(t_ch3oh,1)*dt_c
    A(t_hcho,2) = A(t_hcho,1)*dt_c
    A(t_hcooh,2) = A(t_hcooh,1)*dt_c
    A(t_hoch2oh,2) = A(t_hoch2oh,1)*dt_c
    A(t_hoch2ooh,2) = A(t_hoch2ooh,1)*dt_c

    A(t_ch3ooh,3) = (cc0(i_ch3ooh) + production(i_ch3ooh)*dt_c)*(A(t_ch3ooh,1)**2)*dt_c
    A(t_ch3oh,3) = (cc0(i_ch3oh) + production(i_ch3oh)*dt_c)*(A(t_ch3oh,1)**2)*dt_c
    A(t_hcho,3) = (cc0(i_hcho) + production(i_hcho)*dt_c)*(A(t_hcho,1)**2)*dt_c
    A(t_hcooh,3) = (cc0(i_hcooh) + production(i_hcooh)*dt_c)*(A(t_hcooh,1)**2)*dt_c
    A(t_hoch2oh,3) = (cc0(i_hoch2oh) + production(i_hoch2oh)*dt_c)*(A(t_hoch2oh,1)**2)*dt_c
    A(t_hoch2ooh,3) = (cc0(i_hoch2ooh) + production(i_hoch2ooh)*dt_c)*(A(t_hoch2ooh,1)**2)*dt_c
ENDIF 

! 4.4 : Odd-Hydrogen Coefficients
! ===============================
A_hox(1) = 1./( 1. + loss(i_hox)*dt_c )
A_hox(2) = A_hox(1)*dt_c 
A_hox(3) = cc_hox_next*A_hox(1)*dt_c/cc(i_Hox)
A_hox(4) = cc_hox_next*A_hox(1)*dt_c*loss(i_Hox)/(cc(i_Hox))

! 4.5 : Odd-Chlorine Coefficients
! ===============================

IF (igcm_cl.ne.0) THEN 
  cc_clox_next = (cc0(i_clox) + production(i_clox)*dt_c)/(1. + loss(i_clox)*dt_c)
  A_clox(1) = 1./( 1. + loss(i_clox)*dt_c )
  A_clox(2) = A_clox(1)*dt_c
  A_clox(3) = cc_clox_next*A_clox(2)/cc(i_clox)
  A_clox(4) = A_clox(3)*loss(i_clox)
ENDIF


! =====================================================
! 5.0 : Calculations of the Linearised Number Densities 
! =====================================================

! 5.1 : Inorganic Chemistry 
! -------------------------
! CO2
x_j = (t_co2-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_co2,1)*dcc0_dpq( x_j , : ) + A(t_co2,2)*dP_dPQ(t_co2,:) &
                    - A(t_co2,3)*dL_dPQ( t_co2, :)
! CO
x_j = (t_co-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_co,1)*dcc0_dpq( x_j , : ) + A(t_co,2)*dP_dPQ(t_co,:) &
                    - A(t_co,3)*dL_dPQ( t_co, :)
! O2
x_j = (t_o2-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_o2,1)*dcc0_dpq( x_j , : ) + A(t_o2,2)*dP_dPQ(t_o2,:) &
                    - A(t_o2,3)*dL_dPQ( t_o2, :)
! H2
x_j = (t_h2-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_h2,1)*dcc0_dpq( x_j , : ) + A(t_h2,2)*dP_dPQ(t_h2,:) &
                    - A(t_h2,3)*dL_dPQ( t_h2, :)
! H2O
x_j = (t_h2ovap-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_h2ovap,1)*dcc0_dpq( x_j , : ) + A(t_h2ovap,2)*dP_dPQ(t_h2ovap,:) &
                    - A(t_h2ovap,3)*dL_dPQ( t_h2ovap, :)
! H2O2
x_j = (t_h2o2-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_h2o2,1)*dcc0_dpq( x_j , : ) + A(t_h2o2,2)*dP_dPQ(t_h2o2,:) &
                    - A(t_h2o2,3)*dL_dPQ( t_h2o2, :)
! CH4
x_j = (t_ch4-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_ch4,1)*dcc0_dpq( x_j , : ) + A(t_ch4,2)*dP_dPQ(t_ch4,:) &
                    - A(t_ch4,3)*dL_dPQ( t_ch4, :)
IF ( sza > 95. ) THEN 
    ! O 
    x_j = (t_o-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_o,1)*dcc0_dpq( x_j , : ) + A(t_o,2)*dP_dPQ(t_o,:) &
                        - A(t_o,3)*dL_dPQ( t_o, :)
    ! O3 
    x_j = (t_o3-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_o3,1)*dcc0_dpq( x_j , : ) + A(t_o3,2)*dP_dPQ(t_o3,:) &
                        - A(t_o3,3)*dL_dPQ( t_o3, :)
    ! O1D forced to 0
    x_j = (t_o1d-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = 0.
ENDIF 

! HOx 
dHOX_dPQ(lyr_m,:) = A_hox(1)*dHOX0_dPQ(lyr_m,:) &
                  + A_hox(2)*dPhox_dPQ &
                  - A_hox(3)*dLhox_dPQ &
                  + A_hox(4)*dHOX_dPQ(lyr_m,:)

  

! 5.2 : Organic Chemistry 
! -----------------------
IF ( igcm_ch3 .ne. 0 ) THEN 
    ! CH3 
    x_j = (t_ch3-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_ch3,1)*dP_dPQ(t_ch3,:) - A(t_ch3,2)*dL_dPQ(t_ch3,:)

    ! CH3O2 
    x_j = (t_ch3o2-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_ch3o2,1)*dP_dPQ(t_ch3o2,:) - A(t_ch3o2,2)*dL_dPQ(t_ch3o2,:)

    ! CH3OOH 
    x_j = (t_ch3ooh-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_ch3ooh,1)*dcc0_dpq( x_j , : ) + A(t_ch3ooh,2)*dP_dPQ(t_ch3ooh,:) &
                    - A(t_ch3ooh,3)*dL_dPQ( t_ch3ooh, :)

    ! CH3OH 
    x_j = (t_ch3oh-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_ch3oh,1)*dcc0_dpq( x_j , : ) + A(t_ch3oh,2)*dP_dPQ(t_ch3oh,:) &
                        - A(t_ch3oh,3)*dL_dPQ( t_ch3oh, :)

    ! CH3O 
    x_j = (t_ch3o-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_ch3o,1)*dP_dPQ(t_ch3o,:) - A(t_ch3o,2)*dL_dPQ(t_ch3o,:)

    ! HCHO 
    x_j = (t_hcho-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_hcho,1)*dcc0_dpq( x_j , : ) + A(t_hcho,2)*dP_dPQ(t_hcho,:) &
                        - A(t_hcho,3)*dL_dPQ( t_hcho, :)

    ! HCOOH 
    x_j = (t_hcooh-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_hcooh,1)*dcc0_dpq( x_j , : ) + A(t_hcooh,2)*dP_dPQ(t_hcooh,:) &
                        - A(t_hcooh,3)*dL_dPQ( t_hcooh, :)

    ! HOCH2O2
    x_j = (t_hoch2o2-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_hoch2o2,1)*dP_dPQ(t_hoch2o2,:) - A(t_hoch2o2,2)*dL_dPQ(t_hoch2o2,:)

    ! HOCH2OH
    x_j = (t_hoch2oh-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_hoch2oh,1)*dcc0_dpq( x_j , : ) + A(t_hoch2oh,2)*dP_dPQ(t_hoch2oh,:) &
                        - A(t_hoch2oh,3)*dL_dPQ( t_hoch2oh, :)

    ! HOCH2OOH 
    x_j = (t_hoch2ooh-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_hoch2ooh,1)*dcc0_dpq( x_j , : ) + A(t_hoch2ooh,2)*dP_dPQ(t_hoch2ooh,:) &
                        - A(t_hoch2ooh,3)*dL_dPQ( t_hoch2ooh, :)

    ! HCO 
    x_j = (t_hco-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_hco,1)*dP_dPQ(t_hco,:) - A(t_hco,2)*dL_dPQ(t_hco,:)

ENDIF ! igcm_ch3 .ne. 0 

! 5.3 : Chlorine Chemistry 
! ------------------------
IF ( igcm_cl .ne. 0 ) THEN 

    ! ClOx 
    ! -----
    dClOx_dPQ(lyr_m,:) = A_clox(1)*dClOx0_dPQ(lyr_m,:) &
                       + A_clox(2)*dPclox_dPQ &
                       - A_clox(3)*dLclox_dPQ &
                       + A_clox(4)*dClOx_dPQ(lyr_m,:)

    ! ! 5.3.3: HCl 
    x_j = (t_hcl-1)*nlayermx + lyr_m 
    dccn_dpq(x_j,:) = linearised_qssa(dt_c, production(i_hcl), loss(i_hcl), cc0(i_hcl), &
                                  dP_dPQ(t_hcl,:), dL_dPQ(t_hcl,:), dcc0_dpq( x_j, : ) )

    ! ! 5.3.4: HOCl 
    x_j = (t_hocl-1)*nlayermx + lyr_m 
    dccn_dpq(x_j,:) = linearised_qssa(dt_c, production(i_hocl), loss(i_hocl), cc0(i_hocl), &
                                  dP_dPQ(t_hocl,:), dL_dPQ(t_hocl,:), dcc0_dpq( x_j, : ))

    ! 5.3.5: ClO3 
    x_j = (t_clo3-1)*nlayermx + lyr_m 
    dccn_dpq(x_j,:) = linearised_qssa(dt_c, production(i_clo3), loss(i_clo3), cc0(i_clo3), &
                                  dP_dPQ(t_clo3,:), dL_dPQ(t_clo3,:), dcc0_dpq( x_j, : ))

    ! ! 5.3.6: ClO4 
    x_j = (t_clo4-1)*nlayermx + lyr_m 
    dccn_dpq(x_j,:) = linearised_qssa(dt_c, production(i_clo4), loss(i_clo4), cc0(i_clo4), &
                                  dP_dPQ(t_clo4,:), dL_dPQ(t_clo4,:), dcc0_dpq( x_j, : ))

    ! 5.3.7: ClCO 
    x_j = (t_clco-1)*nlayermx + lyr_m 
    dccn_dpq(x_j,:) = linearised_qssa(dt_c, production(i_clco), loss(i_clco), cc0(i_clco), &
                                  dP_dPQ(t_clco,:), dL_dPQ(t_clco,:), dcc0_dpq( x_j, : ))

    ! ! 5.3.8: Cl2 
    x_j = (t_cl2-1)*nlayermx + lyr_m 
    dccn_dpq(x_j,:) = linearised_qssa(dt_c, production(i_cl2), loss(i_cl2), cc0(i_cl2), &
                                  dP_dPQ(t_cl2,:), dL_dPQ(t_cl2,:), dcc0_dpq( x_j, : ))

    ! 5.3.9: Cl2O2 
    x_j = (t_cl2o2-1)*nlayermx + lyr_m 
    dccn_dpq(x_j,:) = linearised_qssa(dt_c, production(i_cl2o2), loss(i_cl2o2), cc0(i_cl2o2), &
                                  dP_dPQ(t_cl2o2,:), dL_dPQ(t_cl2o2,:), dcc0_dpq( x_j, : ))

    ! 5.3.10: OClO 
    x_j = (t_oclo-1)*nlayermx + lyr_m 
    dccn_dpq(x_j,:) = linearised_qssa(dt_c, production(i_oclo), loss(i_oclo), cc0(i_oclo), &
                                  dP_dPQ(t_oclo,:), dL_dPQ(t_oclo,:), dcc0_dpq( x_j, : ))

    ! 5.3.11 : CH3OCl
    x_j = (t_ch3ocl-1)*nlayermx + lyr_m 
    dccn_dpq(x_j,:) = linearised_qssa(dt_c, production(i_ch3ocl), loss(i_ch3ocl), cc0(i_ch3ocl), &
                                  dP_dPQ(t_ch3ocl,:), dL_dPQ(t_ch3ocl,:), dcc0_dpq( x_j, : ))

    ! 5.3.12 : HClO4
    x_j = (t_hclo4-1)*nlayermx + lyr_m 
    dccn_dpq(x_j,:) = linearised_qssa(dt_c, production(i_hclo4), loss(i_hclo4), cc0(i_hclo4), &
                                  dP_dPQ(t_hclo4,:), dL_dPQ(t_hclo4,:), dcc0_dpq( x_j, : ))

    ! 5.3.13 : ClOO
    x_j = (t_cloo-1)*nlayermx + lyr_m 
    dccn_dpq(x_j,:) = linearised_qssa(dt_c, production(i_cloo), loss(i_cloo), cc0(i_cloo), &
                                  dP_dPQ(t_cloo,:), dL_dPQ(t_cloo,:), dcc0_dpq( x_j, : ))



ENDIF



RETURN 

END 