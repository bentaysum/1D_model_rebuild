SUBROUTINE tlm_ox(iter, lyr_m, dens,&
					ro_o3, ro_o3_denominator,&
					dt_c, dt_p, &
					nesp, cc, cc0,&
					j, loss_ox, prod_ox, ccOX_tplus1,&
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
					cab106, cab107)

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
real ro_o3 ! partition function of O and O3 
real ro_o3_denominator ! denominator value of ro_o3 
real dt_c, dt_p ! chemical and physical timestep
integer nesp ! number of species in the chemistry routines
real cc(nesp), cc0(nesp) ! number density of species after and before the
						 ! odd-hydrogen calculations (only H, OH and HO2 are effected)
real j(nd) ! photolysis values 
real loss_ox, prod_ox ! loss and production terms for the Odd Oxygen summed family
real ccOX_tplus1 ! next time-steps value for CC(OX)
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

! ====================================================================================
!									Local Variables
! ====================================================================================

! Initialisation variables
integer o_j, o3_j, o1d_j ! Locations for arrays
integer x_j 
integer iq ! Tracer iterator
! O1D 
real dPo1d_coeff(nqmx), dLo1d_coeff(nqmx) 
real dPo1d_dPQ(nqmx*nlayermx), dLo1d_dPQ(nqmx*nlayermx)
REAL loss_o1d, production_o1d
real A_O1D(2)
real dO1D_dPQ(nqmx*nlayermx)
! O3 and O 
real dN_dPQ(nqmx*nlayermx), dD_dPQ(nlayermx*nqmx), dDd_dPQ(nlayermx*nqmx)
real dN_coef(nqmx), dD_coef(nqmx), dDd_coef(nqmx)
real dro_o3_gamma(4) ! Coefficients for linearising partition function ro_o3
real dro_o3_dPQ(nqmx*nlayermx) ! d(ro_o3)/d(PQ)
real A_O(2), A_O3(2) ! Coefficients for linearising O and O3 
REAL dO_dPQ(nqmx*nlayermx), dO3_dPQ(nqmx*nlayermx) 
! Odd-Oxygen Family
real ox_gamma(4) ! Coefficients
real dPox_dPQ(nqmx*nlayermx), dLox_dPQ(nqmx*nlayermx) ! Linearised loss and production 
real dPox_coef(nqmx), dLox_coef(nqmx) ! Coefficients for linearised loss and prod.

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
		j_glyox_h2, j_ch3cooh, j_ch3coooh, j_ch3cocooh

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

! ============================================ ! 
! STAGE 0: INITIALISATION OF THE ARRAYS 
! ============================================ ! 
!
! 0.1 : Linearised Odd-Oxygen ( cc(i_ox) )
! -------------------------------------------
IF ( iter == 1 ) THEN 

	o_j = (t_o-1)*nlayermx + lyr_m
	o3_j = (t_o3-1)*nlayermx + lyr_m
	
	dOX_dPQ(lyr_m,:) = 0.
	dOX0_dPQ(lyr_m,:) = 0.
	
	dOX_dPQ(lyr_m,:) = (TLM_ident( o3_j, : ) + TLM( o3_j, :)*dt_p) &
					  *Avmr(lyr_m,t_o3)*dens  &
					  + (TLM_ident( o_j, : ) + TLM( o_j, :)*dt_p) &
					  *Avmr(lyr_m,t_o)*dens

	dOX0_dPQ(lyr_m,:) = dOX_dPQ(lyr_m,:)		
			
ENDIF 

! ============================================ ! 
! STAGE 1: Excited Atomic Oxygen ( O(1D) )
! ============================================ ! 
!
! Linearise the production and loss terms
dPo1d_coeff(:) = 0.
dPo1d_coeff(t_co2) = j(j_co2_o1d)
dPo1d_coeff(t_o2) = j(j_o2_o1d)
dPo1d_coeff(t_o3) = j(j_o3_o1d)

dLo1d_coeff(:) = 0.
dLo1d_coeff(t_co2) = b001 
dLo1d_coeff(t_h2ovap) = b002
dLo1d_coeff(t_h2) = b003 
dLo1d_coeff(t_o2) = b004 
dLo1d_coeff(t_o3) = b005 + b006 
dLo1d_coeff(t_ch4) = b007 + b008 + b009 

dPo1d_dPQ(:) = 0.
dLo1d_dPQ(:) = 0.
DO iq = 1, nqmx
	x_j = (iq-1)*nlayermx + lyr_m
	dPo1d_dPQ = dPo1d_dPQ + dPo1d_coeff(iq)*dccn_dpq( x_j, : )
	dLo1d_dPQ = dLo1d_dPQ + dLo1d_coeff(iq)*dccn_dpq( x_j, : )
ENDDO

! Values in Photochemistry.F 
loss_o1d = (b001*cc(i_co2) &
			+ b002*cc(i_h2o) &
			+ b003*cc(i_h2) &
			+ b004*cc(i_o2) &
			+ b005*cc0(i_o3) &
			+ b006*cc0(i_o3) &
			+ b007*cc(i_ch4) &
			+ b008*cc(i_ch4) &
			+ b009*cc(i_ch4))
			
production_o1d = j(j_co2_o1d)*cc(I_co2) &
                        + j(j_o2_o1d)*cc(i_o2) &
                        + j(j_o3_o1d)*cc0(i_o3)
	 

A_O1D(1) =  1./loss_o1d
A_O1D(2) = production_o1d*(A_O1D(1)**2)

dO1D_dPQ = A_O1D(1)*dPo1d_dPQ - A_O1D(2)*dLo1d_dPQ

dccn_dpq( (t_o1d-1)*nlayermx + lyr_m, : ) = dO1D_dPQ


		
		

! ============================================ ! 
! STAGE 2: LINEARISED PARTITION FUNCTION
! ============================================ ! 
!
! ro_o3 = ( N ) / ( D + d/(cc[o]) )
!
! Similar to hox routine...
! dN_dPQ  
! ------
dN_coef(:) = 0. 
dN_coef(t_o) = a003 
dN_coef(t_h) = c003 
dN_coef(t_oh) = c014 
dN_coef(t_ho2) = c015 

IF ( igcm_ch3 .ne. 0 ) dN_coef(t_ch3) = cab004 
IF ( igcm_ch3o2 .ne. 0 ) dN_coef(t_ch3o2) = cab010 
IF ( igcm_ch3o .ne. 0 ) dN_coef(t_ch3o) = cab016

! dD_dPQ 
! ------
dD_coef(:) = 0. 
dD_coef(t_o2) = a001 
dD_coef(t_ch4) = cab002

IF ( igcm_ch3 .ne. 0 ) dD_coef(t_ch3) = cab005 
IF ( igcm_ch3o2 .ne. 0 ) dD_coef(t_ch3o2) = cab012
IF ( igcm_ch3o .ne. 0 ) dD_coef(t_ch3o) = cab017
IF ( igcm_hcho .ne. 0 ) dD_coef(t_hcho) = cab020 
IF ( igcm_hco .ne. 0 ) dD_coef(t_hco) = cab021 
! dDd_dPQ 
! -------
dDd_coef(:) = 0. 
dDd_coef(t_ho2) = cab074*cc(i_ch3cooo) + cab097*cc(i_hoch2co3)

dN_dPQ(:) = 0.
dD_dPQ(:) = 0. 
dDd_dPQ(:) = 0.
DO iq = 1, nqmx 
	x_j = (iq-1)*nlayermx + lyr_m
	dN_dPQ = dN_dPQ + dN_coef(iq)*dccn_dpq( x_j, : )
	dD_dPQ = dD_dPQ + dD_coef(iq)*dccn_dpq( x_j, : )
	dDd_dPQ = dDd_dPQ + dDd_coef(iq)*dccn_dpq( x_j, : )
ENDDO

! Coefficients 
! ------------
dro_o3_gamma(1) = 1./ro_o3_denominator
dro_o3_gamma(2) = ro_o3*dro_o3_gamma(1)
dro_o3_gamma(3) = dro_o3_gamma(2)/cc0(i_o)
dro_o3_gamma(4) = dro_o3_gamma(3)*(cab074*cc(i_ch3cooo)*cc(i_ho2) &
                     + cab097*cc(i_hoch2co3)*cc(i_ho2))/cc0(i_o)

dro_o3_dPQ = dro_o3_gamma(1)*dN_dPQ &
		   - dro_o3_gamma(2)*dD_dPQ &
		   - dro_o3_gamma(3)*dDd_dPQ &
		   + dro_o3_gamma(4)*dccn_dpq( (t_o-1)*nlayermx+lyr_m,:)

! ============================================ ! 
! STAGE 3: LINEARISED O3 
! ============================================ ! 
A_O3(1) =  1./(1.+ro_o3)
A_O3(2) = cc(i_ox)*(A_O3(1)**2)

dO3_dPQ = A_O3(1)*dOX_dPQ( lyr_m, : ) - A_O3(2)*dro_o3_dPQ
! ============================================ ! 
! STAGE 4: LINEARISED O 
! ============================================ ! 
A_O(1) = ro_o3
A_O(2) = cc(i_o3) 

dO_dPQ = A_O(1)*dO3_dPQ + A_O(2)*dro_o3_dPQ 

dccn_dpq( (t_o-1)*nlayermx + lyr_m, : ) = dO_dPQ
dccn_dpq( (t_o3-1)*nlayermx + lyr_m, : ) = dO3_dPQ

! ============================================ ! 
! STAGE 5: LINEARISED OX 
! ============================================ ! 
dPox_coef(:) = 0. 
dPox_coef(t_co2) = j(j_co2_o) + j(j_co2_o1d)
dPox_coef(t_o2) = 2.*j(j_o2_o) + 2.*j(j_o2_o1d)
dPox_coef(t_ho2) = j(j_ho2) + c006*cc(i_h) &
					+ cab074*cc(i_ch3cooo) + cab097*cc(i_hoch2co3)
dPox_coef(t_h) = c006*cc(i_ho2)
dPox_coef(t_oh) = 2.*c013*cc(i_oh)

dLox_coef(:) = 0.
dLox_coef(t_co) = e002*cc(i_o)
dLox_coef(t_o) = 4.*a002*cc(i_o) + 2.*a003*cc(i_o3) &
				+ c001*cc(i_ho2) + c002*cc(i_oh) &
				+ c012*cc(i_h2o2) + e002*cc(i_co) &
				+ cab002*cc(i_ch4) + cab012*cc(i_ch3o2) &
				+ cab017*cc(i_ch3o) + cab020*cc(i_hcho) &
				+ cab021*cc(i_hco) + cab037*cc(i_c2h6) 
dLox_coef(t_o3) = 2.*a003*cc(i_o) + c003*cc(i_h) &
				+ c014*cc(i_oh) + c015*cc(i_ho2) &
				+ cab004*cc(i_ch3) + cab010*cc(i_ch3o2) &
				+ cab016*cc(i_ch3o) 
dLox_coef(t_h) = c003*cc(i_o3) 
dLox_coef(t_oh) = c002*cc(i_o) + c014*cc(i_o3) 
dLox_coef(t_ho2) = c001*cc(i_o) + c015*cc(i_o3) 
dLox_coef(t_h2o2) = c012*cc(i_o) 
dLox_coef(t_ch4) = cab002*cc(i_o)

IF ( igcm_ch3 .ne. 0 ) dLox_coef(t_ch3) = cab004*cc(i_o3)
IF ( igcm_ch3o2 .ne. 0 ) dLox_coef(t_ch3o2) = cab010*cc(i_o3) + cab012*cc(i_o)
IF ( igcm_ch3o .ne. 0 ) dLox_coef(t_ch3o) = cab016*cc(i_o3) + cab017*cc(i_o)
IF ( igcm_hcho .ne. 0 ) dLox_coef(t_hcho) = cab020*cc(i_o)
IF ( igcm_hco .ne. 0 ) dLox_coef(t_hco) = cab021*cc(i_o)

dPox_dPQ(:) = 0.
dLox_dPQ(:) = 0.
DO iq = 1,nqmx
	x_j = (iq-1)*nlayermx + lyr_m 
	dPox_dPQ = dPox_dPQ + dPox_coef(iq)*dccn_dpq(x_j, : )
	dLox_dPQ = dLox_dPQ + dLox_coef(iq)*dccn_dpq(x_j, : )
ENDDO 
	
! Coefficients
ox_gamma(1) = 1./(1. + loss_ox*dt_c) 
ox_gamma(2) = ox_gamma(1)*dt_c
ox_gamma(3) = ccOX_tplus1*ox_gamma(1)*dt_c/cc(i_ox)
ox_gamma(4) = dt_c*ccOX_tplus1*ox_gamma(1)*loss_ox/cc(i_ox)

! Linearisation:
!---------------
dOX_dPQ(lyr_m,:) = ox_gamma(1)*dOX0_dPQ(lyr_m,:) &
				+ ox_gamma(2)*dPox_dPQ &
				- ox_gamma(3)*dLox_dPQ &
				+ ox_gamma(4)*dOX_dPQ(lyr_m,:)

				
RETURN 


		
END 



