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
                    ,cl052, cl053, cl054, &
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
 ,cl052, cl053, cl054

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


! Notes on Chlorine:
! 
! SIBEM : 
!       ClOx 
!       HCl
!       HOCl 
!       CH3OCl 
!       HClO4 
!       ClO3 
!       ClO4
!       Cl2 [SZA > 95.]
!       Cl2O2 [SZA > 95.]
!       OClO [SZA > 95.]
!
! Steady State : 
!       Cl2 [SZA .le. 95.]
!       Cl2O2 [SZA .le. 95.]
!       OClO [SZA .le. 95.]
!  
    
! 0.0 : Initialisation of Array 
! -----------------------------
dP_coeff(:,:) = 0.
dL_coeff(:,:) = 0. 
dP_dPQ(:,:) = 0.
dL_dPQ(:,:) = 0.

! 0.1 : If k_pseudo > 0., we need to linearise this TERMS
! -----------------------------
if ( ( lyr_m .le. 5 ) .and. ( k_pseudo .gt. 1.e10 ) ) then
    dKpseudo_dPQ(lyr_m,t_o1d) = -(b007 + b008 + b009)
    dKpseudo_dPQ(lyr_m,t_o) = -cab002 
    dKpseudo_dPQ(lyr_m,t_oh) = -cab001 
else
    dKpseudo_dPQ(lyr_m,:) = 0.
endif

! Remainder of Chemistry linearisation occurs in this file.
! Two methods exist for calculating number density, Steady
! State [SS], and the Semi-Implicit Backwards Euler Method [SIBEM]
!
! SS: cc = P(PQ)/L(PQ)
!
! SIBEM : cc^t+1 = ( cc0 + P(PQ)*dt_c )
!				 / ( 1 + L(PQ)*dt_c )
! 
! Steady-State assumptions applied to:
!	CH3, CH3O2, CH3O, HCO and HOCH2O2
!
! Remainder via SIBEM. 
!
! ' used to denote d( )/d(PQ) 
!
! Steady-State
! ============
!
!  cc(PQ)' = L^-1 * P(PQ)'
! 		   - P*L^-1 * L(PQ)'
!
! SIBEM
! =====
! cc^t+1 ' = 1/(1 + L*dt)  * [ cc0' + dt*P' ]
!		   - (cc0 + P*dt)/(1 + L*dt)^2 * dt * L' 
!
! ============================================
! STAGE ONE : LINEARISING THE PRODUCTION AND
!             LOSS TERMS 
! ============================================
! 
! For tracer i, the production rate is:
!
! P_i = k_ab * cc_a * cc_b 
!     + k_ac * cc_a * cc_c 
! 	  + ... 
!
! P_i = SUM_{j=1}^{nqmx}( 
!			SUM_{k=1}^{nqmx}(
!							 k_jk * cc_j * cc_k
!							)
!						)
!
! If tracers cc_j and cc_k do not react to produce
! tracer i, then k_jk = 0.
!
! d( P_i )/d(PQ) = SUM_{j=1}^{nqmx}( 
!				       SUM_{k=1}^{nqmx}(
!							 k_jk * d(cc_j)/d(PQ) * cc_k 
!						+    k_jk * d(cc_k)/d(PQ) * cc_j
!									   )
!									)
!
! This results in coefficients appearing constructed
! purely from control model run values, example being:
!
! P_i' = (k_ab*cc_b + k_ad*cc_d + k_ah*cc_h) * cc_a'
!	   + (k_ab*cc_a + k_bc*cc_c + k_bh*cc_h) * cc_b'
!	   + (k_bc*cc_b + k_cd) * cc_c' 
!	   + ... 

! 	1.1: Initialise to prevent issues
! 	---------------------------------
!   1.1.1: CO2 [SIBEM]
!   ------------------
    dP_coeff(t_co2,t_co) = e001*cc(i_oh) + e002*cc(i_o) 
    dP_coeff(t_co2,t_o) = e002*cc(i_co)  
    dP_coeff(t_co2,t_o2) = cab084*cc(i_hcoco) + cab092*cc(i_hoch2co)
    dP_coeff(t_co2,t_oh) = e001*cc(i_co) + cab032*cc(i_hcooh) &
                         + cab067*cc(i_ch3cooh) + 0.91*cab098*cc(i_hoch2co2h) &
                         + cab099*cc(i_hcoco2h) + cab103*cc(i_hcoco3h)
    dP_coeff(t_co2,t_ho2) = cab072*cc(i_ch3cooo) + cab095*cc(i_hoch2co3) &  
                          + cab106*cc(i_hcoco3)

    IF ( igcm_hcooh .ne. 0 ) dP_coeff(t_co2,t_hcooh) = cab032*cc(i_oh)

!   1.1.2: CO [SIBEM] 
!   -----------------
    dP_coeff(t_co,t_co2) = j(j_co2_o1d) + j(j_co2_o) 

    dP_coeff(t_co,t_o) = 0.17*cab005*cc(i_ch3) &
                       + cab021*cc(i_hco) 

    dP_coeff(t_co,t_o2) = cab026*cc(i_hco)

    dP_coeff(t_co,t_h) = cab027*cc(i_hco)

    dP_coeff(t_co,t_oh) = cab025*cc(i_hco)

    IF (igcm_ch3 .ne. 0 ) dP_coeff(t_co,t_ch3) = 0.17*cab005*cc(i_o) &
                                               + 2.*cab022*cc(i_hco) 

    IF (igcm_hcho .ne. 0 ) dP_coeff(t_co,t_hcho) = j(j_ch2o_co)

    IF (igcm_hco .ne. 0 ) dP_coeff(t_co,t_hco) = cab021*cc(i_o) &
                                               + cab022*cc(i_ch3) &
                                               + cab024*cc(i_hco) & 
                                               + cab025*cc(i_oh) &
                                               + cab026*cc(i_o2) &
                                               + cab027*cc(i_h) 

!   1.1.3: O2 [SIBEM]
!   -----------------
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
                       + cl002*cc(i_clo)

    dP_coeff(t_o2,t_h) = c003*cc(i_o3) &
                        + c005*cc(i_ho2)

    dP_coeff(t_o2,t_oh) = c002*cc(i_o) &
                        + c007*cc(i_ho2) &
                        + c014*cc(i_o3) &
                        + 0.06*cl012*cc(i_clo)

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

    IF (igcm_ch3 .ne. 0 ) dP_coeff(t_o2,t_ch3) = cab004*cc(i_o3)

    IF (igcm_ch3o2 .ne. 0 ) dP_coeff(t_o2,t_ch3o2) = cab006*cc(i_ho2) &
                                                   + cab007*cc(i_ho2) &
                                                   + 0.5*cab008*ro2 &
                                                   + 0.5*cab008*cc(i_ch3o2) &
                                                   + 0.5*cab009*ro2 &
                                                   + 0.5*cab009*cc(i_ch3o2) &
                                                   + 2.*cab010*cc(i_o3) &
                                                   + cab012*cc(i_o) &
                                                   + 0.5*cab031*cc(i_hoch2o2) &
                                                   + cl020*cc(i_clo)

    IF (igcm_ch3o .ne. 0 ) dP_coeff(t_o2,t_ch3o) = cab016*cc(i_o3) &
                                                 + 0.75*cab017*cc(i_o)

    IF (igcm_Hoch2o2 .ne. 0 ) dP_coeff(t_o2,t_hoch2o2) = 0.8*cab029*cc(i_ho2) &
                                                       + 0.5*cab008*cc(i_ch3o2) &
                                                       + 0.5*cab009*cc(i_ch3o2) &
                                                       + 0.5*cab031*cc(i_hoch2o2) &
                                                       + 0.5*cab031*ro2


!   1.1.4: H2 [SIBEM]
!   -----------------
    dP_coeff(t_h2,t_o) = 0.17*cab005*cc(i_ch3)

    dP_coeff(t_h2,t_h) = c005*cc(i_ho2) & 
                       + 2.*c008*cc(i_h) &
                       + cab027*cc(i_hco) &
                       + cl041*cc(i_hcl)

    dP_coeff(t_h2,t_ho2) = c005*cc(i_h) 

    dP_coeff(t_h2,t_o1d) = b009*cc(i_ch4)

    dP_coeff(t_h2,t_ch4) = b009*cc(i_o1d)  &
                         + j(j_ch4_1ch2_h2) + j(j_ch4_ch_h2_h) 

    IF (igcm_ch3 .ne. 0 ) dP_coeff(t_h2,t_ch3) = 0.17*cab005*cc(i_o)

    IF (igcm_hcho .ne. 0 ) dP_coeff(t_h2,t_hcho) = j(j_ch2o_co)

!   1.1.5: H2O [SIBEM]
!   ------------------
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

    IF (igcm_ch3o2 .ne. 0 ) dP_coeff(t_h2ovap,t_ch3o2) = cab007*cc(i_oh)

    IF (igcm_ch3oh .ne. 0 ) dP_coeff(t_h2ovap,t_ch3oh) = cab013*cc(i_oh)

    IF (igcm_ch3ooh .ne. 0 ) dP_coeff(t_h2ovap,t_ch3ooh) = cab015*cc(i_oh)

    IF (igcm_hcho .ne. 0 ) dP_coeff(t_h2ovap,t_hcho) = cab018*cc(i_oh)

    IF (igcm_hco .ne. 0 ) dP_coeff(t_h2ovap,t_hco) = cab025*cc(i_oh)

    IF (igcm_hoch2o2 .ne. 0 ) dP_coeff(t_h2ovap,t_hoch2o2) = 0.3*cab009*cc(i_ho2)

    IF (igcm_hcooh .ne. 0  )  dP_coeff(t_h2ovap,t_hcooh) = cab032*cc(i_oh)

    IF (igcm_hoch2ooh .ne. 0 ) dP_coeff(t_h2ovap,t_hoch2ooh) = cab034*cc(i_oh) 

    IF (igcm_hoch2oh .ne. 0 ) dP_coeff(t_h2ovap,t_hoch2oh) = cab035*cc(i_oh) 

!   1.1.6: H2O2 [SIBEM]
!   -------------------
    dP_coeff(t_h2o2,t_oh) = 2.D0*c017*cc(i_oh)

    dP_coeff(t_h2o2,t_ho2) = 2.D0*c008*cc(i_ho2) + 2.D0*c016*cc(i_ho2) 


!   1.1.7: CH4 [SIBEM]
!   ------------------ 
    IF ( igcm_ch3 .ne. 0 ) dP_coeff(t_ch4,t_ch3) = cab022*cc(i_hco)
    IF ( igcm_hco .ne. 0 ) dP_coeff(t_ch4,t_hco) = cab022*cc(i_ch3)


!   1.1.8: CH3 [Steady-State]
!   ------------------------- 
    IF ( igcm_ch3 .ne. 0 ) THEN 
        dP_coeff( t_ch3, t_o1d ) = b007*cc(i_ch4) 

        dP_coeff( t_ch3, t_o ) = 0.51*cab002*cc(i_ch4) &
                               + 0.75*cab017*cc(i_ch3o)

        dP_coeff( t_ch3, t_oh ) = cab001*cc(i_ch4)

        dP_coeff( t_ch3, t_ch4 ) = b007*cc(i_o1d) &
                                 + cab001*cc(i_oh) &
                                 + 0.51*cab002*cc(i_o) &
                                 + j(j_ch4_ch3_h) &
                                 + cl016*cc(i_cl)

        IF ( igcm_ch3o .ne. 0 ) dP_coeff(t_ch3,t_ch3o) = 0.75*cab017*cc(i_o)
    ENDIF 

!   1.1.9: CH3O2 [Steady-State]
!   --------------------------- 
    IF ( igcm_ch3o2 .ne. 0 ) THEN 
          dP_coeff( t_ch3o2, t_ch3) = cab003*cc(i_o2)
          dP_coeff( t_ch3o2, t_o2) = cab003*cc(i_ch3)
          dP_coeff( t_ch3o2, t_o3) = cab016*cc(i_ch3o)
          dP_coeff( t_ch3o2, t_oh) = 0.6*cab014*cc(i_ch3ooh)
          dP_coeff( t_ch3o2, t_ch3ooh) = 0.6*cab014*cc(i_oh)
          dP_coeff( t_ch3o2, t_ch3o) = cab016*cc(i_o3) 
    ENDIF 

!   1.1.10: CH3OOH [SIBEM]
!   ----------------------
    IF ( igcm_ch3ooh .ne. 0 ) THEN 
       dP_coeff(t_ch3ooh,t_ho2) = cab006*cc(i_ch3o2)
       dP_coeff(t_ch3ooh,t_ch3o2) = cab006*cc(i_ho2)
    ENDIF 

!   1.1.10: CH3OH [SIBEM]
!   ----------------------
    IF ( igcm_ch3oh .ne. 0 ) THEN 
       dP_coeff(t_ch3oh,t_ch3o2) = 0.5*cab009*(cc(i_ch3o2) + ro2) 
       dP_coeff(t_ch3oh,t_hoch2o2) = 0.5*cab009*cc(i_ch3o2) 

       dP_coeff(t_ch3oh,t_oh) =  cab107*cc(i_ch3)
       dP_coeff(t_ch3oh,t_ch3) = cab107*cc(i_oh)
    ENDIF 

!   1.1.11: CH3O [Steady-State]
!   ---------------------------
    IF ( igcm_ch3o .ne. 0 ) THEN 
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
                                 + cl020*cc(i_clo) &
                                 + 0.5*cl022*cc(i_cl)

        dP_coeff(t_ch3o,t_ch3ooh) = j(j_ch3o2h)   

        dP_coeff(t_ch3o,t_hoch2o2) = cab008*cc(i_ch3o2) 

    ENDIF 

!   1.1.12: HCHO [SIBEM]
!   --------------------
    IF ( igcm_hcho .ne. 0 ) THEN 
        
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

        dP_coeff(t_hcho, t_hco) = cab024*cc(i_hco)*2.d0 

    ENDIF 

!   1.1.13: HCOOH [SIBEM]
!   --------------------- 
    IF ( igcm_hcooh .ne. 0 ) THEN 

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


    ENDIF   

!   1.1.14: HOCH2O2 [SIBEM]
!   -----------------------
    IF ( igcm_hoch2o2 .ne. 0 ) THEN 
        dP_coeff(t_hoch2o2,t_hcho) = cab019*cc(i_ho2)
        dP_coeff(t_hoch2o2,t_ho2) = cab019*cc(i_hcho) 
        dP_coeff(t_hoch2o2,t_hoch2ooh) = cab033*cc(i_oh)
        dP_coeff(t_hoch2o2,t_oh) = cab033*cc(i_hoch2ooh) 

    ENDIF 

!   1.1.15: HOCH2OH [SIBEM]
!   ----------------------- 
    IF ( igcm_hoch2oh .ne. 0 ) THEN 
            dP_coeff(t_hoch2oh,t_ch3o2) = 0.5*cab031*cc(i_hoch2o2) 
            dP_coeff(t_hoch2oh,t_hoch2o2) = 0.5*cab031*( cc(i_hoch2o2) + ro2 )
    ENDIF 

!   1.1.16: HOCH2OOH [SIBEM]
!   ------------------------
    IF ( igcm_hoch2ooh .ne. 0 ) THEN 
        dP_coeff(t_hoch2ooh,t_ho2) = cab029*0.5*cc(i_hoch2o2)
        dP_coeff(t_hoch2ooh,t_hoch2o2) = cab029*0.5*cc(i_ho2)
    ENDIF 

!   1.1.17: HCO [Steady-State]
!   --------------------------
    IF ( igcm_hoch2ooh .ne. 0 ) THEN 
        dP_coeff(t_hco,t_oh) = cab018*cc(i_hcho) 
        dP_coeff(t_hco,t_o) = cab020*cc(i_hcho) 
        dP_coeff(t_hco,t_hcho) = cab018*cc(i_oh) + cab020*cc(i_o) &
                               + j(j_ch2o_co) + cl017*cc(i_cl)
    ENDIF 


!   1.2: NIGHT O3 AND O [SIBEM] 
!   ---------------------------
    IF ( sza > 95. ) THEN 
!       1.2.1: O3 
        dP_coeff(t_o3,t_o) = a001*cc(i_o2)

        dP_coeff(t_o3,t_o2) = a001*cc(i_o)

!       1.2.2: O
        dP_coeff(t_o,t_h) = c006*cc(i_ho2)

        dP_coeff(t_o,t_ho2) = c006*cc(i_h) 

        dP_coeff(t_o,t_oh) = 2.*c013*cc(i_oh)
    ENDIF 

     !   1.3 : Chlorine Chemistry 
!         Contributions 
!   ------------------------
    IF ( igcm_cl .ne. 0) THEN 
        ! CO 
        dP_coeff(t_co,t_clco) =cl033 
        !O2 
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

         dP_coeff(t_o2,t_oclo) = cl052*cc(i_o3) 

         ! H2
         dP_coeff(t_h2,t_hcl) = cl041*cc(i_h)

         ! H2O
         dP_coeff(t_h2ovap,t_hcl) = cl014*cc(i_oh)
         dP_coeff(t_h2ovap,t_hocl) = cl015*cc(i_oh)
         dP_coeff(t_h2ovap,t_hclo4) = cl027*cc(i_oh)

         ! CH3 
         IF ( igcm_ch3 .ne. 0 ) dP_coeff(t_ch3,t_cl) = cl016*cc(i_ch4)

         ! CH3O 
         IF ( igcm_ch3o .ne. 0 ) THEN 
         dP_coeff(t_ch3o,t_clo) = cl019*cc(i_ch3o2) &
                                + cl020*cc(i_ch3o2) 
          dP_coeff(t_ch3o,t_cl) = 0.5*cl022*cc(i_ch3o2) &
                                + cl029*cc(i_ch3ocl) 
          dP_coeff(t_ch3o,t_ch3ocl) = cl029*cc(i_cl) &
                                    + cl032*cc(i_oh)
          ENDIF 

          ! HCHO
          IF ( igcm_hcho .ne. 0) dP_coeff(t_hcho,t_cl) =cl018*cc(i_ch3ooh)
          ! HCO 
          IF ( igcm_hco .ne. 0 ) dP_coeff(t_hco,t_cl) = cl017*cc(i_hcho)

    ENDIF 

! Chlorine Compounds : Production  
! ===============================
IF (igcm_cl .ne.0) THEN 

    ! HCl 
    dP_coeff(t_hcl,t_h) = cl037*cc(i_cl2)

    dP_coeff(t_hcl,t_oh) = 0.06*cl012*cc(i_clo)

    dP_coeff(t_hcl,t_ho2) = cl009*cc(i_cl)

    dP_coeff(t_hcl,t_h2) = cl008*cc(i_cl) 

    dP_coeff(t_hcl,t_h2o2) = cl011*cc(i_cl)

    dP_coeff(t_hcl,t_ch4) = cl016*cc(i_cl)

    IF (igcm_ch3 .ne. 0) THEN 
        dP_coeff(t_hcl,t_hcho) = cl017*cc(i_cl)
        dP_coeff(t_hcl,t_ch3ooh) = cl018*cc(i_cl) 
        dP_coeff(t_hcl,t_ch3o2) = 0.5*cl022*cc(i_cl)
        dP_coeff(t_hcl,t_hcooh) = cl043*cc(i_cl)
    ENDIF 

    dP_coeff(t_hcl,t_cl) = cl008*cc(i_h2) &
                         + cl009*cc(i_ho2) &
                         + cl011*cc(i_h2o2) &
                         + cl016*cc(i_ch4) &
                         + cl017*cc(i_hcho) &
                         + cl018*cc(i_ch3ooh) &
                         + 0.5*cl022*cc(i_ch3o2) &
                         + cl031*cc(i_ch3ocl) &
                         + cl043*cc(i_hcooh)

    dP_coeff(t_hcl,t_clo) = 0.06*cl012*cc(i_oh)

    dP_coeff(t_hcl,t_ch3ocl) = cl031*cc(i_cl)

    dP_coeff(t_hcl,t_cl2) = cl037*cc(i_h) 

    ! Cl2 
    dP_coeff(t_cl2,t_cl) = cl025*cc(i_cloo) &
                         + cl029*cc(i_ch3ocl) &
                         + cl030*cc(i_cl2o2)
    dP_coeff(t_cl2,t_clo) = 2.*cl003*cc(i_clo)
    dP_coeff(t_cl2,t_cloo) = cl025*cc(i_cl)
    dP_coeff(t_cl2,t_ch3ocl) = cl029*cc(i_cl)
    dP_coeff(t_cl2,t_cl2o2) = cl030*cc(i_cl)

    ! ClOO 
    dP_coeff(t_cloo,t_o2) = cl028*cc(i_cl) 

    IF (igcm_ch3o2 .ne. 0 ) dP_coeff(t_cloo,t_ch3o2) = cl019*cc(i_clo)
    
    dP_coeff(t_cloo,t_clo) = cl019*cc(i_ch3o2) &
                           + cl045*cc(i_clo3) 
    dP_coeff(t_cloo,t_clo3) = cl045*cc(i_clo) 
    dP_coeff(t_cloo,t_cl) = cl028*cc(i_o2) &
                          + cl030*cc(i_cl2o2)
    dP_coeff(t_cloo,t_cl2o2) = cl030*cc(i_cl) &
                             + j(j_cl2o2)
    ! HOCl
    dP_coeff(t_hocl,t_oh) = cl032*cc(i_ch3ocl) &
                          + cl036*cc(i_cl2)
    dP_coeff(t_hocl,t_ho2) = cl013*cc(i_clo)

    dP_coeff(t_hocl,t_cl2) = cl036*cc(i_oh)
    dP_coeff(t_hocl,t_clo) = cl013*cc(i_ho2)
    dP_coeff(t_hocl,t_ch3ocl) = cl032*cc(i_oh)
    ! Cl2O2 
    dP_coeff(t_cl2o2,t_clo) = 2.*cl006*cc(i_clo)
    ! OClO
    dP_coeff(t_oclo,t_oh) = cl050*cc(i_clo3)

    IF ( igcm_ch3o2 .ne. 0 ) dP_coeff(t_oclo,t_ch3o2) = cl021*cc(i_clo) 

    dP_coeff(t_oclo,t_clo) = 2.*cl005*cc(i_clo) &
                           + cl021*cc(i_ch3o2) &
                           + cl045*cc(i_clo3) &
                           + 2.*cl046*cc(i_clo3)

    dP_coeff(t_oclo,t_clo3) = cl045*cc(i_clo) &
                            + 2.*cl046*cc(i_clo) &
                            + cl050*cc(i_oh)
    ! CH3OCl
    dP_coeff(t_ch3ocl,t_clo) = cl020*cc(i_ch3o2)
    IF (igcm_ch3o2 .ne. 0 ) dP_coeff(t_ch3ocl,t_ch3o2) = cl020*cc(i_clo)
    ! ClCO
    dP_coeff(t_clco,t_cl) = cl023*cc(i_co)
    dP_coeff(t_clco,t_co) = cl023*cc(i_cl)
    ! ClO3 
    dP_coeff(t_clo3,t_o) = cl051*cc(i_oclo)
    dP_coeff(t_clo3,t_o3) = cl044*cc(i_cl) &
                          + cl052*cc(i_oclo)


    dP_coeff(t_clo3,t_cl) = cl044*cc(i_o3) &
                          + cl053*cc(i_clo4) 
    dP_coeff(t_clo3,t_oclo) = cl051*cc(i_o) &
                            + cl052*cc(i_o3) 
    dP_coeff(t_clo3,t_clo4) = cl053*cc(i_cl) 
    ! HClO4
    dP_coeff(t_hclo4,t_clo3) = cl048*cc(i_oh) &
                             + cl049*cc(i_oh) 
    dP_coeff(t_hclo4,t_oh) = cl048*cc(i_clo3) &
                           + cl049*cc(i_clo3)
    dP_coeff(t_hclo4,t_hocl) = cl054*cc(i_clo4) 
    dP_coeff(t_hclo4,t_clo4) = cl054*cc(i_hocl) 
    ! ClO4 
    dP_coeff(t_clo4,t_oh) = cl027*cc(i_hclo4) 
    dP_dPQ(t_clo4,t_hclo4) = cl027*cc(i_oh)

ENDIF 




!   2.0: Construct the Linearised Production Array
!   ----------------------------------------------
    DO iq_j = 1, nqmx
        
        DO iq_i = 1, nqmx
            x_i = (iq_i-1)*nlayermx + lyr_m 
            
            dP_dPQ(iq_j,:) = dP_dPQ(iq_j,:) + &
                            dP_coeff(iq_j,iq_i)*dccn_dpq( x_i, : )
                                
        ENDDO 
        
    ENDDO

! Additional Nitrogen Chemistry
! -----------------------------
    dP_dPQ(t_o2,:) = dP_dPQ(t_o2,:) &
                   + d001*cc(i_o)*dNO2_dPQ(lyr_m,:)

! ===============================================================
! ADDITIONAL CHEMISTRY 
! ===============================================================
!
! 1 - Atmospheric lifetime of CH4 600x faster than standard chem.
! 2 - Surface lifetime of CH4 forced to 1 hour
!
! ===============================================================

! -------------------------------------------------------
! 1 ) new P[CH3] = P[CH3] + [CH4]( (tau_ch4 - 1)*L[CH4] )  
! -------------------------------------------------------
! IF ( sza < 95. ) THEN 
     ! dP_coeff(t_ch3,t_ch4) = dP_coeff(t_ch3,t_ch4) + (methane_enhancement-1) &
                                                    ! *(loss(i_ch4) - j(j_ch4_1ch2_h2) &
                                                    ! - j(j_ch4_3ch2_h_h) - j(j_ch4_ch_h2_h) &
                                                    ! - j(j_ch4_ch3_h))
     ! dP_coeff(t_ch3,t_o1d) = dP_coeff(t_ch3,t_o1d) &
                           ! + cc(i_ch4)*( (methane_enhancement-1)*(b007 + b008 + b009) )   
                           
     ! dP_coeff(t_ch3,t_o) = dP_coeff(t_ch3,t_o) &
                           ! + cc(i_ch4)*( (methane_enhancement-1)*(cab002) ) 
     ! dP_coeff(t_ch3,t_oh) = dP_coeff(t_ch3,t_oh) &
                           ! + cc(i_ch4)*( (methane_enhancement-1)*(cab001) )
                           
! ENDIF 

! ! ! -------------------------------------------------------
! ! ! 2 ) new P[CH3] = P[CH3] + [CH4]( (1/tau_ch4) - L[CH4] ) 
! ! ! -------------------------------------------------------
! IF ( (lyr_m == 1 ) .and. (sza .le. 95.) ) THEN 
     ! ! IF ( sza .le. 95. ) THEN 
          ! dP_coeff(t_ch3,t_ch4) = dP_coeff(t_ch3,t_ch4) + &
                                 ! (1./methane_enhancement) - &
                                 ! ( (b007 + b008 + b009)*cc(i_o1d) &
                                  ! +(cab001)*cc(i_oh) &
                                  ! +(cab002)*cc(i_o) ) 
                                                    
          ! dP_coeff(t_ch3,t_o1d) = dP_coeff(t_ch3,t_o1d) - &
                                 ! cc(i_ch4)*( b007 + b008 + b009 )
                                       
          ! dP_coeff(t_ch3,t_o) = dP_coeff(t_ch3,t_o) - &
                                ! cc(i_ch4)*( cab002 ) 
                                
          ! dP_coeff(t_ch3,t_oh) = dP_coeff(t_ch3,t_oh) - &
                                ! cc(i_ch4)*( cab001 ) 
     ! ! ENDIF     
! ENDIF 


! ============================================
! STAGE TWO : LINEARISING THE LOSS TERMS 
! ============================================
    ! 2.1.1: CO 
    ! ------------------------
    dL_coeff(t_co,t_oh) = e001  
    dL_coeff(t_co,t_o) = e002   

    ! 2.1.2: O2 
    ! ------------------------
    dL_coeff(t_o2,t_o) = a001  
    dL_coeff(t_o2,t_h) = c011   
    IF (igcm_ch3 .ne. 0 ) dL_coeff(t_o2,t_ch3) = cab003  
    IF (igcm_ch3o .ne. 0 ) dL_coeff(t_o2,t_ch3o) = cab015  
    IF (igcm_hco .ne. 0 ) dL_coeff(t_o2,t_hco) = cab026  

    ! 2.1.3: H2 
    ! ------------------------
    dL_coeff(t_h2,t_o1d) =  b003  
    dL_coeff(t_h2,t_oh) = c010   

    ! 2.1.4: H2O 
    ! ------------------------
    dL_coeff(t_h2ovap,t_o1d) =  b002  

    ! 2.1.5: H2O2
    ! ------------------------
    dL_coeff(t_h2o2,t_oh) = c009  
    dL_coeff(t_h2o2,t_o) = c012  

    ! 2.1.6: CH4 
    ! ------------------------
    dL_coeff(t_ch4,t_o1d) = b007 + b008 + b009 
    dL_coeff(t_ch4,t_oh) = cab001   
    dL_coeff(t_ch4,t_o) = cab002   

!   -------------------------------------
!   Forced Lifetime via k_pseudo reaction 
!   -------------------------------------
    if ( (lyr_m .le. 5) .and. (k_pseudo .ge. 1e10) ) &
    dL_coeff(t_ch4,:) = 0. 

    ! 2.1.7: CH3 
    ! -----------
    IF ( igcm_ch3 .ne. 0 ) THEN 
        dL_coeff(t_ch3,t_o2) = cab003   
        dL_coeff(t_ch3,t_o3) = cab004   
        dL_coeff(t_ch3,t_o) = cab005   
        dL_coeff(t_ch3,t_oh) = cab107  
        dL_coeff(t_ch3,t_ch3) = cab038  
        dL_coeff(t_ch3,t_hco) = cab022 + cab023
    ENDIF 

    ! 2.1.8: CH3O2 
    ! ------------ 
    IF ( igcm_ch3o2 .ne. 0 ) THEN 
        dL_coeff(t_ch3o2,t_ho2) = cab006 + cab007 
        dL_coeff(t_ch3o2,t_ch3o2) = cab008 + cab009 
        dL_coeff(t_ch3o2,t_o3) = cab010   
        dL_coeff(t_ch3o2,t_oh) = cab011   
        dL_coeff(t_ch3o2,t_o) = cab012   
        dL_coeff(t_ch3o2,t_hoch2o2) = cab008 + cab009 
    ENDIF 

    ! 2.1.9: CH3OOH 
    ! --------------
    IF ( igcm_ch3ooh .ne. 0 ) THEN 
        dL_coeff(t_ch3ooh,t_oh) = cab014   
    ENDIF 

    ! 2.1.10: CH3OH 
    ! --------------
    IF ( igcm_ch3oh .ne. 0 ) THEN 
        dL_coeff(t_ch3oh,t_oh) = cab013   
    ENDIF 

    ! 2.1.11: CH3O 
    ! -----------
    IF ( igcm_ch3o .ne. 0 ) THEN 
        dL_coeff(t_ch3o,t_o2) = cab015  
        dL_coeff(t_ch3o,t_o) = cab017  
        dL_coeff(t_ch3o,t_o3) = cab016  
    ENDIF 
    ! 2.1.12: HCHO 
    ! ------------ 
    IF ( igcm_hcho .ne. 0 ) THEN 
        dL_coeff(t_hcho,t_oh) = cab018  
        dL_coeff(t_hcho,t_ho2) = cab019  
        dL_coeff(t_hcho,t_o) = cab020   
    ENDIF 

    ! 2.1.13: HCOOH
    ! -------------
    IF ( igcm_hcooh .ne. 0 ) THEN 
        dL_coeff(t_hcooh,t_oh) = cab032  
    ENDIF 

    ! 2.1.14: HOCH2O2 
    ! ---------------
    IF ( igcm_hoch2o2 .ne. 0 ) THEN 
        dL_coeff(t_hoch2o2,t_ho2) = cab029   
        dL_coeff(t_hoch2o2,t_hoch2o2) = cab030 + cab031
        dL_coeff(t_hoch2o2,t_ch3o2) = cab030 + cab031 
    ENDIF 

    ! 2.1.15: HOCH2OH
    ! ---------------
    IF ( igcm_hoch2oh .ne. 0 ) THEN 
        dL_coeff(t_hoch2oh,t_oh) = cab035  
    ENDIF 

    ! 2.1.16: HOCH2OOH 
    ! ---------------
    IF ( igcm_hoch2ooh .ne. 0 ) THEN 
        dL_coeff(t_hoch2ooh,t_oh) = cab033 + cab034
    ENDIF 

    ! 2.1.17: HCO 
    ! -----------
    IF ( igcm_hco .ne. 0 ) THEN 
        dL_coeff(t_hco,t_o) = cab021   
        dL_coeff(t_hco,t_hco) = cab024   
        dL_coeff(t_hco,t_oh) = cab025   
        dL_coeff(t_hco,t_o2) = cab026   
        dL_coeff(t_hco,t_h) = cab027   
        dL_coeff(t_hco,t_ch3) = cab022 + cab023 
    ENDIF 

    ! 2.2: NIGHT O3 AND O 
    ! ------------------------
    IF ( sza > 95. ) THEN 
    ! 2.2.1: O  
    ! --------
        dL_coeff(t_o,t_o) = 2.*a002   
        dL_coeff(t_o,t_o2) = a001   
        dL_coeff(t_o,t_o3) = a003   
        dL_coeff(t_o,t_ho2) = c001   
        dL_coeff(t_o,t_oh) = c002   
        dL_coeff(t_o,t_h2o2) = c012   
        dL_coeff(t_o,t_co) = e002   
        dL_coeff(t_o,t_ch4) = cab002   
        if ( igcm_ch3 .ne. 0 ) dL_coeff(t_o,t_ch3) = cab005   
        if ( igcm_ch3o2 .ne. 0 ) dL_coeff(t_o,t_ch3o2) = cab012   
        if ( igcm_ch3o .ne. 0 ) dL_coeff(t_o,t_ch3o) = cab017  
        if ( igcm_hcho .ne. 0 )  dL_coeff(t_o,t_hcho) = cab020   
        if ( igcm_hco .ne. 0 ) dL_coeff(t_o,t_hco) = cab021  

!    2.2.2: O3 
!    ---------
        dL_coeff(t_o3,t_o) = a003   
        dL_coeff(t_o3,t_h) = c003   
        dL_coeff(t_o3,t_oh) = c014   
        dL_coeff(t_o3,t_ho2) = c015   
        IF ( igcm_ch3 .ne. 0 ) dL_coeff(t_o3,t_ch3) = cab004  
        IF ( igcm_ch3o2 .ne. 0 ) dL_coeff(t_o3,t_ch3o2) = cab010   
        IF ( igcm_ch3o .ne. 0 ) dL_coeff(t_o3,t_ch3o) = cab016  
    ENDIF

!    2.2.3: Chlorine Contributions 
!   ------------------------------
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
     ! CH3 
     dL_coeff(t_ch3,t_cl2) = cl038 
     ! CH3O2
     IF ( igcm_ch3o2 .ne. 0 ) THEN 
          dL_coeff(t_ch3o2,t_clo) = cl019 + cl020 + cl021 
          dL_coeff(t_ch3o2,t_cl) = cl022 
     ENDIF 
     ! CH3OOH
     IF ( igcm_ch3ooh .ne. 0 ) dL_coeff(t_ch3ooh,t_cl) = cl018
     ! HCHO
     IF ( igcm_hcho .ne. 0 ) dL_coeff(t_hcho, t_cl) = cl017 
     ! HCOOH 
     IF ( igcm_hcooh .ne. 0 ) dL_coeff(t_hcooh,t_cl) = cl043

     ! Daylight 
     IF ( sza > 95. ) THEN 
        ! O3 
        dL_coeff(t_o3,t_cl) = cl001 + cl044 
        dL_coeff(t_o3,t_oclo) = cl052 
        ! O 
        dL_coeff(t_o,t_clo) = cl002 
        dL_coeff(t_o,t_hcl) = cl040 
        dL_coeff(t_o,t_hocl) = cl042 
        dL_coeff(t_o,t_oclo) = cl051

     ENDIF 

ENDIF 

! Chlorine Chemistry : Loss Rates 
! ===============================
IF ( igcm_cl .ne. 0 ) THEN 
    ! HCl
    dL_coeff(t_hcl,t_oh) = cl014 
    dL_coeff(t_hcl,t_o1d) = 0.88*cl039 
    dL_coeff(t_hcl,t_o) = cl040 
    dL_coeff(t_hcl,t_h) = cl041 
    ! Cl2 
    dL_coeff(t_cl2,t_o1d) = cl035
    dL_coeff(t_cl2,t_oh) = cl036 
    dL_coeff(t_cl2,t_h) = cl037 
    dL_coeff(t_cl2,t_ch3) = cl038 
    ! ClOO
    dL_coeff(t_cloo,t_cl) = cl024 + cl025 
    ! HOCl
    dL_coeff(t_hocl,t_oh) = cl015
    dL_coeff(t_hocl,t_o) = cl042 
    dL_coeff(t_hocl,t_clo4) = cl054
    ! Cl2O2 
    dL_coeff(t_cl2o2,t_cl) = cl030 
    ! OClO
    dL_coeff(t_oclo,t_o) = cl051
    dL_coeff(t_oclo,t_o3) = cl052 
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
ENDIF 

!   2.3: Create Linearised loss array 
!   ---------------------------------
    DO iq_j = 1,nqmx

        DO iq_i = 1, nqmx
            x_i = (iq_i-1)*nlayermx + lyr_m 
            dL_dPQ(iq_j,:) = dL_dPQ(iq_j,:) &
                        + dL_coeff(iq_j,iq_i)*dccn_dpq( x_i,:)  
        ENDDO 

    ENDDO 

! ================
! Optional Outputs
! ================
!
! 1: Steady-State O2 VMR
! ----------------------
! At surface, we want to assess
! how sensitive the steady-state O2 
! [P/L] is to tracer species in the 
! model.
!
! Equation breaks down to:
!
! d[O2_steady]/d[PQ] =
!       *( 1/L * dP/dPQ - P/L^2 * dL/dPQ) * d(PQ^t)/d(PQ^t0)
!
! where d(PQ^t)/d(PQ^t0) is the sensitivity array of lyr_m 
! for all tracers at forecast timestep t with respect to 
! tracers at backtrace timestep t0, calculated via the
! adjoint externally. 
if ( lyr_m == 1 ) then 
o2_coefficient_array(:) = 0.
    do iq_i = 1,nqmx
         o2_coefficient_array( iq_i ) = (1./loss(i_o2))*dP_coeff(t_o2,iq_i) &
                   - (production(i_o2)/(loss(i_o2)**2.))*dL_coeff(t_o2,iq_i)
    enddo

endif 



! ===============================================================
! ADDITIONAL CHEMISTRY 
! ===============================================================
!
! 1 - Atmospheric lifetime of CH4 600x faster than standard chem.
! 2 - Surface lifetime of CH4 forced to 1 hour
!
! ===============================================================

! ----------------------------------------------------
! 1 ) new L[CH4] = SUM( J(l,j_ch4_XXX) ) + Tau_ch4*L[CH4] 
! -------------------------------------------------
! IF ( sza < 95. ) THEN 
          ! ! CH4 Loss 
          ! ! --------
          ! dL_coeff(t_ch4,t_oh) = dL_coeff(t_ch4,t_oh)*methane_enhancement
          ! dL_coeff(t_ch4,t_o) = dL_coeff(t_ch4,t_o)*methane_enhancement
          ! dL_coeff(t_ch4,t_o1d) = dL_coeff(t_ch4,t_o1d)*methane_enhancement
! ENDIF 

! ----------------------------------------------------
! 2 ) new L[CH4] = SUM( J(l,j_ch4_XXX) ) + (1/Tau_ch4)
! ----------------------------------------------------
! IF ( (lyr_m == 1) .and. (sza .le. 95.) ) THEN 
     ! ! IF ( sza < 95. ) THEN 
          ! dL_coeff(t_ch4,:) = 0. 
     ! ! ENDIF     
! ENDIF 


! ============================================
! STAGE THREE : EQUATION COEFFICIENTS
! ============================================
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

    A(:,:) = 0.

    A(t_co2,1) = 1./( 1. + loss(i_co2)*dt_c)
    A(t_co,1) = 1./( 1. + loss(i_co)*dt_c)
    A(t_o2,1) = 1./( 1. + loss(i_o2)*dt_c)
    A(t_h2ovap,1) = 1./( 1. + loss(i_h2o)*dt_c)
    A(t_h2,1) = 1./( 1. + loss(i_h2)*dt_c)
    A(t_h2o2,1) = 1./( 1. + loss(i_h2o2)*dt_c)
    A(t_ch4,1) = 1./( 1. + loss(i_ch4)*dt_c)

    ! Methane Oxidation
    ! ------------------------
    IF ( igcm_ch3ooh .ne. 0 ) THEN 
    A(t_ch3ooh,1) = 1./( 1. + loss(t_ch3ooh)*dt_c)
    ENDIF 
    IF ( igcm_ch3oh .ne. 0 ) THEN 
    A(t_ch3oh,1) = 1./( 1. + loss(t_ch3oh)*dt_c)
    ENDIF 
    IF ( igcm_hcho .ne. 0 ) THEN 
    A(t_hcho,1) = 1./( 1. + loss(t_hcho)*dt_c)
    ENDIF 
    IF ( igcm_hcho .ne. 0 ) THEN 
    A(t_hcooh,1) = 1./( 1. + loss(t_hcooh)*dt_c)
    ENDIF
    IF ( igcm_hoch2oh .ne. 0 ) THEN 
    A(t_hoch2oh,1) = 1./( 1. + loss(t_hoch2oh)*dt_c)
    ENDIF
    IF ( igcm_hoch2ooh .ne. 0 ) THEN 
    A(t_hoch2ooh,1) = 1./( 1. + loss(t_hoch2ooh)*dt_c)
    ENDIF

    IF ( sza > 95. ) THEN
    A(t_o,1) = 1./( 1. + loss(i_o)*dt_c)
    A(t_o3,1) = 1./( 1. + loss(i_o3)*dt_c)
    ENDIF 

    ! Chlorine 
    ! --------
    IF ( igcm_cl .ne. 0 ) THEN 
        A(t_hcl,1) = 1./( 1. + loss(t_hcl)*dt_c)
        A(t_hocl,1) = 1./( 1. + loss(t_hocl)*dt_c)
        A(t_ch3ocl,1) = 1./( 1. + loss(t_ch3ocl)*dt_c)
        A(t_hclo4,1) = 1./( 1. + loss(t_hclo4)*dt_c)
        A(t_clo3,1) = 1./( 1. + loss(t_clo3)*dt_c)
        A(t_clo4,1) = 1./( 1. + loss(t_clo4)*dt_c)

        IF ( sza > 95. ) THEN 
            A(t_cl2,1) = 1./( 1. + loss(t_cl2)*dt_c)
            A(t_cl2o2,1) = 1./( 1. + loss(t_cl2o2)*dt_c)
            A(t_oclo,1) = 1./( 1. + loss(t_oclo)*dt_c)
        ENDIF 

    ENDIF 

    A(:,2) = A(:,1)*dt_c 

    A(t_co2,3) = (cc0(i_co2) + production(i_co2)*dt_c)*(A(t_co2,1)**2)*dt_c
    A(t_co,3) = (cc0(i_co) + production(i_co)*dt_c)*(A(t_co,1)**2)*dt_c
    A(t_o2,3) = (cc0(i_o2) + production(i_o2)*dt_c)*(A(t_o2,1)**2)*dt_c
    A(t_h2,3) = (cc0(i_h2) + production(i_h2)*dt_c)*(A(t_h2,1)**2)*dt_c
    A(t_h2o2,3) = (cc0(i_h2o2) + production(i_h2o2)*dt_c)*(A(t_h2o2,1)**2)*dt_c
    A(t_h2ovap,3) = (cc0(i_h2o) + production(i_h2o)*dt_c)*(A(t_h2ovap,1)**2)*dt_c
    A(t_ch4,3) = (cc0(i_ch4) + production(i_ch4)*dt_c)*(A(t_ch4,1)**2)*dt_c

    ! Methane Oxidation
    ! ------------------------
    IF ( igcm_ch3ooh .ne. 0 ) THEN 
    	A(t_ch3ooh,3) = (cc0(i_ch3ooh) + production(i_ch3ooh)*dt_c)*(A(t_ch3ooh,1)**2)*dt_c
    ENDIF 		
    IF ( igcm_ch3oh .ne. 0 ) THEN 
    	A(t_ch3oh,3) = (cc0(i_ch3oh) + production(i_ch3oh)*dt_c)*(A(t_ch3oh,1)**2)*dt_c
    ENDIF 		
    IF ( igcm_hcho .ne. 0 ) THEN 
    	A(t_hcho,3) = (cc0(i_hcho) + production(i_hcho)*dt_c)*(A(t_hcho,1)**2)*dt_c
    ENDIF 		
    IF ( igcm_hcooh .ne. 0 ) THEN 
    	A(t_hcooh,3) = (cc0(i_hcooh) + production(i_hcooh)*dt_c)*(A(t_hcooh,1)**2)*dt_c
    ENDIF 		
    IF ( igcm_hoch2oh .ne. 0 ) THEN 
    	A(t_hoch2oh,3) = (cc0(i_hoch2oh) + production(i_hoch2oh)*dt_c)*(A(t_hoch2oh,1)**2)*dt_c
    ENDIF 		
    IF ( igcm_hoch2ooh .ne. 0 ) THEN 
    	A(t_hoch2ooh,3) = (cc0(i_hoch2ooh) + production(i_hoch2ooh)*dt_c)*(A(t_hoch2ooh,1)**2)*dt_c
    ENDIF 		

    IF ( sza > 95. ) THEN 
    	A(t_o,3) = (cc0(i_o) + production(i_o)*dt_c)*(A(t_o,1)**2)*dt_c
    	A(t_o3,3) = (cc0(i_o3) + production(i_o3)*dt_c)*(A(t_o3,1)**2)*dt_c
    ENDIF 

    ! Chlorine Chemistry 
    ! ------------------------
    IF ( igcm_cl .ne. 0 ) THEN 
        A(t_hcl,3) = (cc0(i_hcl) + production(i_hcl)*dt_c)*(A(t_hcl,1)**2)*dt_c
        A(t_hocl,3) = (cc0(i_hocl) + production(i_hocl)*dt_c)*(A(t_hocl,1)**2)*dt_c
        A(t_ch3ocl,3) = (cc0(i_ch3ocl) + production(i_ch3ocl)*dt_c)*(A(t_ch3ocl,1)**2)*dt_c
        A(t_hclo4,3) = (cc0(i_hclo4) + production(i_hclo4)*dt_c)*(A(t_hclo4,1)**2)*dt_c
        A(t_clo3,3) = (cc0(i_clo3) + production(i_clo3)*dt_c)*(A(t_clo3,1)**2)*dt_c
        A(t_clo4,3) = (cc0(i_clo4) + production(i_clo4)*dt_c)*(A(t_clo4,1)**2)*dt_c

        IF ( sza > 95. ) THEN 
            A(t_cl2,3) = (cc0(i_cl2) + production(i_cl2)*dt_c)*(A(t_cl2,1)**2)*dt_c
            A(t_cl2o2,3) = (cc0(i_cl2o2) + production(i_cl2o2)*dt_c)*(A(t_cl2o2,1)**2)*dt_c
            A(t_oclo,3) = (cc0(i_oclo) + production(i_oclo)*dt_c)*(A(t_oclo,1)**2)*dt_c

        ELSE 
            A(t_cl2,1) = 1./loss(i_cl2)
            A(t_cl2,2) = production(i_cl2)*(A(t_cl2,1)**2.)

            A(t_cl2o2,1) = 1./loss(i_cl2o2)
            A(t_cl2o2,2) = production(i_cl2o2)*(A(t_cl2o2,1)**2.)

            A(t_oclo,1) = 1./loss(i_oclo)
            A(t_oclo,2) = production(i_oclo)*(A(t_oclo,1)**2.)

        ENDIF 

    ! Steady States
        A(t_clco,1) = 1./loss(i_clco)
        A(t_clco,2) = production(i_clco)*(A(t_clco,1)**2.)

        A(t_cloo,1) = 1./loss(i_cloo)
        A(t_cloo,2) = production(i_cloo)*(A(t_cloo,1)**2.)


    ENDIF 

    ! Steady State Assumptions 
    ! ------------------------
    IF ( igcm_ch3 .ne. 0 ) THEN 
    	A(t_ch3,1) = 1./loss(i_ch3)
    	A(t_ch3,2) = production(i_ch3)*(A(t_ch3,1)**2.)
    ENDIF
    IF ( igcm_ch3o2 .ne. 0 ) THEN 
    	A(t_ch3o2,1) = 1./loss(i_ch3o2)
    	A(t_ch3o2,2) = production(i_ch3o2)*(A(t_ch3o2,1)**2.)
    ENDIF 
    IF ( igcm_ch3o .ne. 0 ) THEN 
    	A(t_ch3o,1) = 1./loss(i_ch3o)
    	A(t_ch3o,2) = production(i_ch3o)*(A(t_ch3o,1)**2.)
    ENDIF 
    IF ( igcm_hoch2o2 .ne. 0 ) THEN 
    	A(t_hoch2o2,1) = 1./loss(i_hoch2o2)
    	A(t_hoch2o2,2) = production(i_hoch2o2)*(A(t_hoch2o2,1)**2.)
    ENDIF 
    IF ( igcm_hco .ne. 0 ) THEN 
    	A(t_hco,1) = 1./loss(i_hco)
    	A(t_hco,2) = production(i_hco)*(A(t_hco,1)**2.)
    ENDIF 


! ============================================
! STAGE 5 : ODD-HYDROGEN FAMILY 
! ============================================
!
! cc_HOx^t+1 = (cc0_HOx + P*dt)/(1 + L/cc_HOx^t *dt )
!
! Linearised Production 
dPhox_coeff(:) = 0. 

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

IF ( igcm_ch3 .ne. 0 ) dPhox_coeff(t_ch3) = 0.956*cab004*cc(i_o3) &
                                          + cab005*cc(i_o)     

IF ( igcm_ch3o .ne. 0 ) dPhox_coeff(t_ch3o) = cab015*cc(i_o2) &
                                            + 0.25*cab017*cc(i_o)  

IF ( igcm_hcho .ne. 0 ) dPhox_coeff(t_hcho) = cab020*cc(i_o) &
                                            + 2.d0*j(j_ch2o_hco)

IF ( igcm_ch3ooh .ne. 0 ) dPhox_coeff(t_ch3ooh) = j(j_ch3o2h) &
                                                + cl018*cc(i_cl)

IF ( igcm_ch3oh .ne. 0 ) dPhox_coeff(t_ch3oh) = j(j_ch3oh)

IF ( igcm_hco .ne. 0 ) dPhox_coeff(t_hco) = cab021*cc(i_o) &
                                          + cab026*cc(i_o2)

IF ( igcm_hoch2o2 .ne. 0 ) dPhox_coeff(t_hoch2o2) = cab028 &
                                                  + cab030*( cc(i_hoch2o2) + ro2 )

IF ( igcm_ch3o2 .ne. 0 ) dPhox_coeff(t_ch3o2) = cab030*cc(i_hoch2o2)

IF ( igcm_cl .ne. 0 ) THEN 
    dPhox_coeff(t_cl) = cl008*cc(i_h2) &
                      + cl011*cc(i_h2o2) &
                      + cl018*cc(i_ch3ooh) 
    dPhox_coeff(t_hcl) = j(j_hcl) &
                       + 0.88*cl039*cc(i_o1d) &
                       + cl040*cc(i_o) 
    dPhox_coeff(t_hocl) = j(j_hocl) &
                        + cl042*cc(i_o) &
                        + cl054*cc(i_clo4)
    dPhox_coeff(t_clo4) = cl054*cc(i_hocl) 
ENDIF 


dPhox_dPQ(:) = 0.
DO iq = 1,nqmx
        x_j = (iq-1)*nlayermx + lyr_m
        dPhox_dPQ = dPhox_dPQ + dPhox_coeff(iq)*dccn_dpq( x_j, : )
ENDDO 


! Linearised Loss 
dLhox_coeff(:) = 0.

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
                  + cl027*cc(i_hclo4)

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

IF (igcm_ch3o2 .ne. 0) dLhox_coeff(t_ch3o2) = cab006*cc(i_ho2) &
                                            + cab007*cc(i_ho2)

IF (igcm_hoch2o2 .ne. 0 ) dLhox_coeff(t_hoch2o2) = 0.6*cab029*cc(i_ho2)

IF (igcm_ch3oh .ne. 0) dLhox_coeff(t_ch3oh) = 0.15*cab013*cc(i_oh)

IF (igcm_ch3ooh .ne. 0) dLhox_coeff(t_ch3ooh) = 0.6*cab014*cc(i_oh)

IF (igcm_hcho .ne. 0) dLhox_coeff(t_hcho) = cab018*cc(i_oh) &
                                          + cab019*cc(i_ho2)

IF (igcm_hco .ne. 0 ) dLhox_coeff(t_hco) = cab025*cc(i_oh) &
                                         + cab027*cc(i_h)

IF (igcm_hoch2ooh .ne. 0 ) dLhox_coeff(t_hoch2ooh) = cab034*cc(i_oh)

IF (igcm_cl .ne. 0 ) THEN 
    dLhox_coeff(t_cl) = cl009*cc(i_ho2) 
    dLhox_coeff(t_clo) = 0.06*cl012*cc(i_oh) &
                       + cl013*cc(i_ho2)
    dLhox_coeff(t_clo3) = cl048*cc(i_oh) &  
                        + cl049*cc(i_oh)
    dLhox_coeff(t_hcl) = cl014*cc(i_oh) &
                       + cl041*cc(i_h)
    dLhox_coeff(t_hocl) = cl015*cc(i_oh)
    dLhox_coeff(t_ch3ocl) = cl032*cc(i_oh)
    dLhox_coeff(t_cl2) = cl036*cc(i_oh) &
                       + cl037*cc(i_h)
    dLhox_coeff(t_hclo4) = cl027*cc(i_oh)
ENDIF 

dLhox_dPQ(:) = 0.
DO iq = 1,nqmx
	x_j = (iq-1)*nlayermx + lyr_m
	dlhox_dPQ = dLhox_dPQ + dLhox_coeff(iq)*dccn_dpq( x_j, : )
ENDDO 					


! Coefficients for HOx 
A_hox(1) = 1./( 1. + loss(i_hox)*dt_c )
A_hox(2) = A_hox(1)*dt_c 
A_hox(3) = cc_hox_next*A_hox(1)*dt_c/cc(i_Hox)
A_hox(4) = cc_hox_next*A_hox(1)*dt_c*loss(i_Hox)/(cc(i_Hox))

dHOX_dPQ(lyr_m,:) = A_hox(1)*dHOX0_dPQ(lyr_m,:) &
                  + A_hox(2)*dPhox_dPQ &
                  - A_hox(3)*dLhox_dPQ &
                  + A_hox(4)*dHOX_dPQ(lyr_m,:)

! =======================================
! Stage 5.1 : OPTIONAL CHLORINE CHEMISTRY 
! =======================================
IF ( igcm_cl .ne. 0 ) THEN 
    A_clox(1) = 1./( 1. + loss(i_clox)*dt_c)
    A_clox(2) = A_clox(1)*dt_c 
    A_clox(3) = (cc0(i_clox) + production(i_clox)*dt_c)*(A_clox(1)**2)*dt_c

    ! Linearised Production of ClOx
    dPclox_coeff(:) = 0.
    dPclox_coeff(t_h) = cl037*cc(i_cl2) & 
                      + cl041*cc(i_hcl) 
    dPclox_coeff(t_o1d) = 2.*cl035*cc(I_cl2) &
                        + 0.88*cl039*cc(i_hcl)
    dPclox_coeff(t_o) = cl040*cc(i_hcl) &
                      + cl042*cc(i_hocl)

    dPclox_coeff(t_oh) = cl014*cc(i_hcl) &
                       + cl015*cc(i_oh) &
                       + cl036*cc(i_cl2)
    dPclox_coeff(t_cl2o2) = cl007*2. &
                          + j(j_cl2o2) 
    dPclox_coeff(t_hcl) = cl014*cc(i_oh) &
                        + 0.88*cl039*cc(i_o1d) &
                        + cl040*cc(i_o) &
                        + cl041*cc(i_h) &
                        + j(j_hcl)

    dPclox_coeff(t_hocl) = cl015*cc(i_oh) &
                         + cl042*cc(i_o) &
                         + j(j_hocl) 

    dPclox_coeff(t_cl) = cl024*cc(i_cloo) 
    dPclox_coeff(t_cloo) = cl024*cc(i_cl) &
                         + cl026  &
                         + j(j_cloo) 
    dPclox_coeff(t_oclo) = j(j_oclo) 
    dPclox_coeff(t_clco) = cl033 
    dPclox_coeff(t_cl2) = 2.*cl035*cc(i_o1d) &
                        + cl036*cc(i_oh) &
                        + cl037*cc(I_h) &
                        + cl038*cc(i_ch3) &
                        + j(j_cl2)*2. 

    IF (igcm_ch3 .ne.0) dPclox_coeff(t_ch3) = cl038*cc(i_cl2)

    ! Linearised Loss of ClOX 
    dLclox_coeff(:) = 0.
    dLclox_coeff(t_co) = cl023 
    dLclox_coeff(t_o3) = cl044
    dLclox_coeff(t_o2) = cl028
    dLclox_coeff(t_h2) = cl008 
    dLclox_coeff(t_oh) = 0.06*cl012
    dLclox_coeff(t_ho2) = cl009 + cl013
    dLclox_coeff(t_h2o2) = cl011
    dLclox_coeff(t_ch4) = cl016 

    IF (igcm_ch3.ne.0) THEN 
        dLclox_coeff(t_hcho) = cl017
        dLclox_coeff(t_ch3ooh) = cl018
        dLclox_coeff(t_ch3o2) = cl019 + cl020 + cl021 + 0.5*cl022
        dLclox_coeff(t_hcooh) = cl043
    ENDIF 

    dLclox_coeff(t_clo) = 2.*cl003 + cl005 + 2.*cl006 
    dLclox_coeff(t_cloo) = cl025
    dLclox_coeff(t_ch3ocl) = cl029 + cl031
    dLclox_coeff(t_cl2o2) = cl030
    dLclox_coeff(t_clo3) = cl045 + cl046 + cl047 

    ! Linearised Production and Loss calculations
    dPclox_dPQ(:) = 0.
    dLclox_dPQ(:) = 0. 
    DO iq = 1, nqmx
        x_j = (iq-1)*nlayermx + lyr_m 
        dLclox_dPQ = dLclox_dPQ&
                    + dLclox_coeff(iq)*dccn_dpq( x_j,:)  
        dPclox_dPQ = dPclox_dPQ &
                    + dPclox_coeff(iq)*dccn_dpq( x_j,:)  
    ENDDO 

    ! ClOx Number Densities Linearised 
    ! dClOx_dPQ(lyr_m,:) = A_clox(1)*dClOx0_dPQ(lyr_m,:) &
    !                    + A_clox(2)*dPclox_dPQ &
    !                    - A_clox(3)*dLclox_dPQ
    dClOx_dPQ(lyr_m,:) = dClOx0_dPQ(lyr_m,:)

ENDIF 

! ============================================
! STAGE 6 : CALCULATIONS AND STORAGE IN ARRAYS
! ============================================

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
					
! Organics 
! ========
IF ( igcm_ch3 .ne. 0 ) THEN 
	x_j = (t_ch3-1)*nlayermx + lyr_m 
	dccn_dpq( x_j, : ) = A(t_ch3,1)*dP_dPQ(t_ch3,:) - A(t_ch3,2)*dL_dPQ(t_ch3,:)
ENDIF
IF ( igcm_ch3o2 .ne. 0 ) THEN 
	x_j = (t_ch3o2-1)*nlayermx + lyr_m 
	dccn_dpq( x_j, : ) = A(t_ch3o2,1)*dP_dPQ(t_ch3o2,:) - A(t_ch3o2,2)*dL_dPQ(t_ch3o2,:)
ENDIF
IF ( igcm_ch3ooh .ne. 0 ) THEN 
x_j = (t_ch3ooh-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_ch3ooh,1)*dcc0_dpq( x_j , : ) + A(t_ch3ooh,2)*dP_dPQ(t_ch3ooh,:) &
					- A(t_ch3ooh,3)*dL_dPQ( t_ch3ooh, :)
ENDIF 
IF ( igcm_ch3oh .ne. 0 ) THEN 
x_j = (t_ch3oh-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_ch3oh,1)*dcc0_dpq( x_j , : ) + A(t_ch3oh,2)*dP_dPQ(t_ch3oh,:) &
					- A(t_ch3oh,3)*dL_dPQ( t_ch3oh, :)
ENDIF 
IF ( igcm_ch3o .ne. 0 ) THEN 
	x_j = (t_ch3o-1)*nlayermx + lyr_m 
	dccn_dpq( x_j, : ) = A(t_ch3o,1)*dP_dPQ(t_ch3o,:) - A(t_ch3o,2)*dL_dPQ(t_ch3o,:)
ENDIF
IF ( igcm_hcho .ne. 0 ) THEN 
x_j = (t_hcho-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_hcho,1)*dcc0_dpq( x_j , : ) + A(t_hcho,2)*dP_dPQ(t_hcho,:) &
					- A(t_hcho,3)*dL_dPQ( t_hcho, :)
ENDIF 
IF ( igcm_hcooh .ne. 0 ) THEN 
x_j = (t_hcooh-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_hcooh,1)*dcc0_dpq( x_j , : ) + A(t_hcooh,2)*dP_dPQ(t_hcooh,:) &
					- A(t_hcooh,3)*dL_dPQ( t_hcooh, :)
ENDIF 
IF ( igcm_hoch2o2 .ne. 0 ) THEN 
	x_j = (t_hoch2o2-1)*nlayermx + lyr_m 
	dccn_dpq( x_j, : ) = A(t_hoch2o2,1)*dP_dPQ(t_hoch2o2,:) - A(t_hoch2o2,2)*dL_dPQ(t_hoch2o2,:)
ENDIF
IF ( igcm_hoch2oh .ne. 0 ) THEN 
x_j = (t_hoch2oh-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_hoch2oh,1)*dcc0_dpq( x_j , : ) + A(t_hoch2oh,2)*dP_dPQ(t_hoch2oh,:) &
					- A(t_hoch2oh,3)*dL_dPQ( t_hoch2oh, :)
ENDIF 
IF ( igcm_hoch2ooh .ne. 0 ) THEN 
x_j = (t_hoch2ooh-1)*nlayermx + lyr_m 
dccn_dpq( x_j, : ) = A(t_hoch2ooh,1)*dcc0_dpq( x_j , : ) + A(t_hoch2ooh,2)*dP_dPQ(t_hoch2ooh,:) &
					- A(t_hoch2ooh,3)*dL_dPQ( t_hoch2ooh, :)
ENDIF 
IF ( igcm_hco .ne. 0 ) THEN 
	x_j = (t_hco-1)*nlayermx + lyr_m 
	dccn_dpq( x_j, : ) = A(t_hco,1)*dP_dPQ(t_hco,:) - A(t_hco,2)*dL_dPQ(t_hco,:)
ENDIF

! ==================
! Chlorine Chemistry 
! ==================
IF (igcm_cl.ne.0) THEN 
    ! HCl
    x_j = (t_hcl-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_hcl,1)*dcc0_dpq( x_j , : ) + A(t_hcl,2)*dP_dPQ(t_hcl,:) &
                    - A(t_hcl,3)*dL_dPQ( t_hcl, :)
    ! HOCl
    x_j = (t_hocl-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_hocl,1)*dcc0_dpq( x_j , : ) + A(t_hocl,2)*dP_dPQ(t_hocl,:) &
                    - A(t_hocl,3)*dL_dPQ( t_hocl, :)
    ! CH3OCl
    x_j = (t_ch3ocl-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_ch3ocl,1)*dcc0_dpq( x_j , : ) + A(t_ch3ocl,2)*dP_dPQ(t_ch3ocl,:) &
                    - A(t_ch3ocl,3)*dL_dPQ( t_ch3ocl, :)
    ! HClO4
    x_j = (t_hclo4-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_hclo4,1)*dcc0_dpq( x_j , : ) + A(t_hclo4,2)*dP_dPQ(t_hclo4,:) &
                    - A(t_hclo4,3)*dL_dPQ( t_hclo4, :)
    ! ClO4
    x_j = (t_clo4-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_clo4,1)*dcc0_dpq( x_j , : ) + A(t_clo4,2)*dP_dPQ(t_clo4,:) &
                    - A(t_clo4,3)*dL_dPQ( t_clo4, :)
    ! ClO3
    x_j = (t_clo3-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_clo3,1)*dcc0_dpq( x_j , : ) + A(t_clo3,2)*dP_dPQ(t_clo3,:) &
                    - A(t_clo3,3)*dL_dPQ( t_clo3, :)
    ! ClCO 
    x_j = (t_clco-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_clco,1)*dP_dPQ(t_clco,:) - A(t_clco,2)*dL_dPQ(t_clco,:)
    ! ClOO 
    x_j = (t_cloo-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = A(t_cloo,1)*dP_dPQ(t_cloo,:) - A(t_cloo,2)*dL_dPQ(t_cloo,:)


    ! Daylight 
    ! --------
    IF ( sza > 95. ) THEN 
        ! Cl2
        x_j = (t_cl2-1)*nlayermx + lyr_m 
        dccn_dpq( x_j, : ) = A(t_cl2,1)*dcc0_dpq( x_j , : ) + A(t_cl2,2)*dP_dPQ(t_cl2,:) &
                        - A(t_cl2,3)*dL_dPQ( t_cl2, :)
        ! Cl2O2
        x_j = (t_cl2o2-1)*nlayermx + lyr_m 
        dccn_dpq( x_j, : ) = A(t_cl2o2,1)*dcc0_dpq( x_j , : ) + A(t_cl2o2,2)*dP_dPQ(t_cl2o2,:) &
                        - A(t_cl2o2,3)*dL_dPQ( t_cl2o2, :)
        ! OClO
        x_j = (t_oclo-1)*nlayermx + lyr_m 
        dccn_dpq( x_j, : ) = A(t_oclo,1)*dcc0_dpq( x_j , : ) + A(t_oclo,2)*dP_dPQ(t_oclo,:) &
                      - A(t_oclo,3)*dL_dPQ( t_oclo, :)
    ELSE 
        ! Cl2
        x_j = (t_cl2-1)*nlayermx + lyr_m 
        dccn_dpq( x_j, : ) = A(t_cl2,1)*dP_dPQ(t_cl2,:) - A(t_cl2,2)*dL_dPQ(t_cl2,:)
        ! Cl2O2
        x_j = (t_cl2o2-1)*nlayermx + lyr_m 
        dccn_dpq( x_j, : ) = A(t_cl2o2,1)*dP_dPQ(t_cl2o2,:) - A(t_cl2o2,2)*dL_dPQ(t_cl2o2,:)
        ! OClO
        x_j = (t_oclo-1)*nlayermx + lyr_m 
        dccn_dpq( x_j, : ) = A(t_oclo,1)*dP_dPQ(t_oclo,:) - A(t_oclo,2)*dL_dPQ(t_oclo,:)

    ENDIF 

ENDIF 




END 