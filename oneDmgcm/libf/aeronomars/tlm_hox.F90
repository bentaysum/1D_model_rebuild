SUBROUTINE tlm_hox(iter, lyr_m, dens,sza,&
					rh_ho2, roh_ho2, &
					dt_c, dt_p, &
					nesp, cc, cc_prev,&
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
real sza 
real rh_ho2, roh_ho2 ! partition functions for H, OH, HO2 
real rh_ho2_denominator, roh_ho2_denominator ! Partition function denominators
real rh_ho2_cabN, rh_ho2_cabD ! Additional terms from organic 
							  ! CAABA reactions to the numerator (N)
							  ! and denominator (D) of the rh_ho2 
							  ! function, which require division
							  ! via number density to be dimensionless.
real roh_ho2_cabN, roh_ho2_cabD ! Additional terms from organic 
							    ! CAABA reactions to the numerator (N)
							    ! and denominator (D) of the roh_ho2 
							    ! function, which require division
							    ! via number density to be dimensionless.
real dt_c, dt_p ! chemical and physical timestep
integer nesp ! number of species in the chemistry routines
real cc(nesp), cc_prev(nesp) ! number density of species after and before the
						 ! odd-hydrogen calculations (only H, OH and HO2 are effected)
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

! Local Variables 
! ===============
integer h_j, oh_j, ho2_j ! Indices for H, HO2, and OH locations for initialisation
integer x_j ! Tracer indice for initialisation
integer	iq ! Tracer iterations
		
! Stage 1 
! -------
REAL A_RH(6), A_ROH(6)

REAL B_H(4,nqmx), B_OH(4,nqmx)
REAL dN_dPQ(nqmx*nlayermx), dNn_dPQ(nqmx*nlayermx), &
		dD_dPQ(nqmx*nlayermx), dDd_dPQ(nqmx*nlayermx)
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
! 0.1 : Linearised Odd-Hydrogen ( cc(i_hox) )
! -------------------------------------------
!
! cc(i_hox) is equal to:
!			( VMR(H) + VMR(OH) + VMR(HO2) )*density 
! 		=	( PQ(H)/mmol(H) + PQ(OH)/mmol(OH) + PQ(HO2)/mmol(HO2) )*mmean*density
! at very first time-step.

IF ( iter == 1 ) THEN 
	h_j = (t_h-1)*nlayermx + lyr_m
	oh_j = (t_oh-1)*nlayermx + lyr_m
 	ho2_j = (t_ho2-1)*nlayermx + lyr_m
	
	dHOX_dPQ(lyr_m,:) = 0. 

	dHOX_dPQ(lyr_m,:) = (TLM_ident( h_j, : ) + TLM( h_j, :)*dt_p) &
					  *Avmr(lyr_m,t_h)*dens  &
					  + (TLM_ident( oh_j, : ) + TLM( oh_j, :)*dt_p) &
					  *Avmr(lyr_m,t_oh)*dens & 
					  + (TLM_ident( ho2_j, : ) + TLM( ho2_j, :)*dt_p) &
					  *Avmr(lyr_m,t_ho2)*dens 
					  
	dHOX0_dPQ(lyr_m,:) = dHOX_dPQ(lyr_m,:)
					  
! 0.2 : Linearised Number Densities ( cc(X) )
! -------------------------------------------
!
! cc(X) is equal to:
!	(mmean/mmol(X)) * ( PQ(X) ) * density 
! at very first time-step.
	DO iq = 1, nqmx
			
		Avmr(lyr_m,iq) = mmol(iq)/mmean(1,lyr_m)
			
		x_j = ( iq - 1 )*nlayermx + lyr_m 
			
		dccn_dpq( x_j, : ) = ( TLM_ident( x_j, :) + TLM( x_j, :)*dt_p) &
							*Avmr(lyr_m,iq)*dens
		
		dcc0_dpq( x_j, : ) = dccn_dpq( x_j, : )
		
	ENDDO	

ENDIF 

! Forcing O(1D) to zero at night
IF ( sza > 95. ) THEN 
	dccn_dpq( (t_o1d-1)*nlayermx + lyr_m, : ) = 0.
	dcc0_dpq( (t_o1d-1)*nlayermx + lyr_m, : ) = 0.
	dccn_dpq( :, (t_o1d-1)*nlayermx + 1 : t_o1d*nlayermx) = 0.
	dcc0_dpq( :, (t_o1d-1)*nlayermx + 1 : t_o1d*nlayermx) = 0.
ENDIF 

! =================================================!
! Stage 0.1: Creating the segments of the partition
! 	 		 functions 
!
! - rh_ho2 = (A + rh_ho2_cabN)/rh_ho2_denominator
! 		   == (A + rh_ho2_cabN)/(B + rh_ho2_cabD)

! - roh_ho2 = (B + roh_ho2_cabN)/roh_ho2_denominator
!		    == (B + roh_ho2_cabN)/(B + roh_ho2_cabD)
! =================================================! 
rh_ho2_cabN = cab002*cc(i_o)*cc(i_ch4)*0.49 &
                + cab005*cc(i_ch3)*cc(i_o3)*0.956 &
                + cab006*cc(i_ch3)*cc(i_o) &
                + j(j_ch4_ch3_h)*cc(i_ch4) &
                + j(j_ch4_3ch2_h_h)*cc(i_ch4)*2. &
                + j(j_ch4_ch_h2_h)*cc(i_ch4) &
                + j(j_ch2o_hco)*cc(i_hcho) &
                + j(j_ch3oh)*cc(i_ch3oh) &
                + j(j_ch3cho_h)*cc(i_ch3cho)

rh_ho2_cabD = cab011*cc(i_ch3o2)*cc(i_oh) &
                + cab013*cc(i_ch3oh)*cc(i_oh) &
                + cab015*cc(i_ch3o)*cc(i_o2) &
                + cab028*cc(i_hoch2o2) &
                + cab030*cc(i_hoch2o2)*cc(i_ro2) &
                + cab032*cc(i_hcooh)*cc(i_oh) &
                + cab035*cc(i_hoch2oh)*cc(i_oh) &
                + cab041*cc(i_c2h5)*cc(i_o2) &
                + 0.6*cab044*cc(i_c2h5o2)*cc(i_ro2) &
                + 0.95*cab047*cc(i_c2h5oh)*cc(i_oh) &
                + cab051*cc(i_hoch2ch2o)*cc(i_o2) &
                + cab052*cc(i_hoch2ch2o) &
                + cab053*cc(i_ethgly)*cc(i_oh) &
                + cab056*cc(i_hyetho2h)*cc(i_oh) &
                + cab062*cc(i_ch2choh)*cc(i_oh) &
                + cab064*cc(i_ch3choho2) &
                + cab078*cc(i_hcoch2o2)*cc(i_ro2)*0.6 &
                + cab090*cc(i_hochcho) &
                + cab093*cc(i_hoch2co3)*cc(i_ro2) &
                + cab098*cc(i_hoch2co2h)*cc(i_oh) &
                + cab099*cc(i_hcoco2h)*cc(i_oh) & 
                + cab101*cc(i_hoch2co3h)*cc(i_oh) &
                + cab104*cc(i_hcoco3)*cc(i_ro2) &
                + j(j_hoch2ooh)*cc(i_hoch2ooh) &
                + j(j_ch3o2h)*cc(i_hoch2ooh) &
                + j(j_ch3o2h)*cc(i_hoch2co3h) &
                + 2.*3.95*j(j_ch2o_co)*cc(i_hcoco2h) &
                + (j(j_ch3o2h) + j(j_hoch2cho_co) &
                   + j(j_hoch2cho_hco) + j(j_hoch2cho_oh)) &
                   *cc(i_hcoco3h) &
                + (j(j_ch3o2h) + j(j_hoch2cho_co) &
                   + j(j_hoch2cho_hco) + j(j_hoch2cho_oh)) &
                   *cc(i_hooch2cho) 

rh_ho2_denominator = c011*cc(i_o2) &
           + cab027*cc(i_hco) &
           + cab042*cc(i_c2h5) &
           + (rh_ho2_cabD)/MAX(cc_prev(i_h),dens*1.e-30)

roh_ho2_cabN = cab002*cc(i_o)*cc(i_ch4)*0.51 &
                 + cab017*cc(i_ch3o)*cc(i_o)*0.25 &
                 + cab020*cc(i_hcho)*cc(i_o) &
                 + cab021*cc(i_hco)*cc(i_o) &
                 + cab037*cc(i_c2h6)*cc(i_o) &
                 + cab066*cc(i_ch3choho2)*cc(i_ro2) &
                 + cab071*cc(i_ch3co)*cc(i_o2) &
                 + cab084*cc(i_hcoco)*cc(i_o2) &
                 + cab092*cc(i_hoch2co)*cc(i_o2) &
                 + j(j_ch3coooh)*cc(i_ch3coooh) &
                 + j(j_ch3o2h)*cc(i_ch3ooh) &
                 + j(j_hoch2cho_oh)*cc(i_hoch2cho) & 
                 + j(j_ch3o2h)*cc(i_ch3chohooh) &
                 + j(j_ch3o2h)*cc(i_hyetho2h) &
                + 2.*b002*cc(i_o1d)*cc(i_h2o) &
                    + b003*cc(i_o1d)*cc(i_h2) &
                    + j(j_h2o)*cc(i_h2o) &
               + b007*cc(i_ch4)*cc(i_o1d)


roh_ho2_cabD = cab015*cc(i_ch3o)*cc(i_o2) &
                 + cab026*cc(i_hco)*cc(i_o2) &
                 + cab028*cc(i_hoch2o2) &
                 + cab030*cc(i_hoch2o2)*cc(i_ro2) &
                 + cab041*cc(i_c2h5)*cc(i_o2) &
                 + 0.6*cab044*cc(i_c2h5o2)*cc(i_ro2) &
                 + cab051*cc(i_hoch2ch2o)*cc(i_o2) &
                 + cab052*cc(i_hoch2ch2o) &
                 + cab064*cc(i_ch3choho2) &
                 + cab078*cc(i_hcoch2o2)*cc(i_ro2)*0.6 &
                 + cab090*cc(i_hochcho) &
                 + cab093*cc(i_hoch2co3)*cc(i_ro2) &
                 + cab104*cc(i_hcoco3)*cc(i_ro2) &
                 + 2.*3.95*j(j_ch2o_co)*cc(i_hcoco2h) 

roh_ho2_denominator = c002*cc(i_o) &
            + c007*cc_prev(i_ho2) &
            + c009*cc(i_h2o2) &        ! ajout 20090401
            + 2.*c013*cc_prev(i_oh) &        ! ajout 20090401
            + 2.*c017*cc_prev(i_oh) &        ! ajout 20090401
            + e001*cc(i_co) &
            + cab001*cc(i_ch4) &
            + 2.*cab011*cc(i_ch3o2) &
            + 0.15*cab013*cc(i_ch3oh) &
            + 0.6*cab014*cc(i_ch3ooh) &
            + cab018*cc(i_hcho) &
            + cab025*cc(i_hco) &
            + 2.*cab032*cc(i_hcooh) &
            + cab033*cc(i_hoch2ooh) &
            + 2.*cab035*cc(i_hoch2oh) &
            + cab036*cc(i_c2h6) &
            + cab045*cc(i_c2h5ooh) &
            + cab047*cc(i_c2h5oh)*1.95 &
            + 2.*cab053*cc(i_ethgly) &
            + cab054*cc(i_hyetho2h) &
            + 2.*cab056*cc(i_hyetho2h) &
            + cab057*cc(i_ch3cho) &
            + cab058*cc(i_ch3cho) &
            + 2.*cab062*cc(i_ch2choh) &
            + cab067*cc(i_ch3cooh) &
            + cab069*cc(i_ch3chohooh) &
            + cab077*cc(i_ch3coooh) &
            + cab081*cc(i_glyox) &
            + cab085*cc(i_hooch2cho) &
            + cab088*cc(i_hoch2cho) &
            + cab089*cc(i_hoch2cho) &
            + cab098*cc(i_hoch2co2h)*2. &
            + cab099*cc(i_hcoco2h)*.2 &
            + cab100*cc(i_hoch2co3h) &
            + cab101*cc(i_hoch2co3h)*2. &
            + cab102*cc(i_hcoco3h) &
            + cab107*cc(i_ch3) &
            + (roh_ho2_cabD/MAX(cc_prev(i_oh),dens*1.e-30))

! ============================================ ! 
! STAGE 1: LINEARISING PARTITION FUNCTIONS 
! ============================================ ! 
!
! rh_ho2 and roh_ho2 have the form:
!
! rX_ho2 = ( N + n/cc(i_ho2) ) / ( D + d/cc(i_X) ) 
!
! Linearisng with respect to mixing ratios:
! 	d(rX_ho2)/d(PQ) = A_RX(1) * dN/dPQ 
!					+ A_RX(2) * dn/dPQ 
!					- A_RX(3) * d(cc(i_ho2))/dPQ 
!					- A_RX(4) * dD/dPQ 
!					- A_RX(5) * dd/dPQ 
!					+ A_RX(6) * d(cc(i_x))/dPQ 

! ----------
! 1.1 RH_HO2
! ----------
A_RH(1) = 1./rh_ho2_denominator 
A_RH(2) = A_RH(1)/cc_prev(i_ho2) 
A_RH(3) = A_RH(2)*rh_ho2_cabN/cc_prev(i_ho2)
A_RH(4) = rh_ho2*A_RH(1) 
A_RH(5) = A_RH(4)/cc_prev(i_h)
A_RH(6) =  A_RH(5)*rh_ho2_cabD/cc_prev(i_h)

! N and D = a001 * cc(X) + a002 * cc(Y) ... 
! so d( N )/d(PQ) = SUM_{i=1}^{nqmx} B_i * d(cc(i))/d(PQ)

! RH_HO2 : dN_dPQ 
! ---------------
B_h(:,:) = 0.
B_h(1,t_h) = c004 + c005 + c006 
B_h(1,t_oh) = c007 
B_h(1,t_ho2) = 2.*c008 + 2.*c016
B_h(1,t_o) = c001 
B_h(1,t_o3) = c015
IF ( igcm_ch3o2 .ne. 0 ) B_H(1,t_ch3o2) = cab006 + cab007
IF ( igcm_hcho .ne. 0 ) B_H(1,t_hcho) = cab019 
IF ( igcm_hoch2o2 .ne. 0 ) B_H(1,t_hoch2o2) = 0.8*cab029 

! RH_HO2 : dn_dPQ 
! ---------------
B_h(2,t_o) = 0.49*cab002*cc(i_ch4) + cab006*cc(i_ch3) 
B_h(2,t_o3) = 0.956*cab005*cc(i_ch3)
B_h(2,t_ch4) = 0.49*cab002*cc(i_o) + j(j_ch4_ch3_h) &
			+ j(j_ch4_3ch2_h_h)*2. + j(j_ch4_ch_h2_h)
			
IF ( igcm_ch3 .ne. 0 ) B_h(2,t_ch3) = cab006*cc(i_o) + 0.956*cab005*cc(i_o3)
IF ( igcm_ch3oh .ne. 0 ) B_h(2,t_ch3oh) = j(j_ch3oh)
IF ( igcm_hcho .ne. 0 ) B_h(2,t_hcho) = j(j_ch2o_hco)
			
! RH_HO2 : dD_dPQ 
! ---------------
B_h(3,t_o2) = c011 
IF ( igcm_hco .ne. 0 ) B_h(3,t_hco) = cab027
! RH_HO2 : dDd_dPQ 
! ---------------
B_h(4,t_oh) = cab011*cc(i_ch3o2) + cab013*cc(i_ch3oh) &
			+ cab032*cc(i_hcooh) + cab035*cc(i_hoch2oh) &
			+ 0.95*cab047*cc(i_c2h5oh) + cab053*cc(i_ethgly) &
			+ cab056*cc(i_hyetho2h) + cab062*cc(i_ch2choh) &
			+ cab098*cc(i_hoch2co2h) + cab099*cc(i_hcoco2h) &
			+ cab101*cc(i_hoch2co3h)
B_h(4,t_o2) = cab015*cc(i_ch3o) + cab041*cc(i_c2h5) &
			+ cab051*cc(i_hoch2ch2o)
IF ( igcm_ch3o2 .ne. 0 ) B_H(4,t_ch3o2) = cab011*cc(i_oh) + cab030*cc(i_hoch2o2) &
										+ 0.6*cab044*cc(i_c2h5o2) + cab078*cc(i_hcoch2o2)*0.6 &
										+ cab093*cc(i_hoch2co3) + cab104*cc(i_hcoco3)
IF ( igcm_ch3oh .ne. 0 ) B_H(4,t_ch3oh) = cab013*cc(i_oh)
IF ( igcm_ch3o .ne. 0 ) B_H(4,t_ch3o) = cab015*cc(i_o2)
IF ( igcm_hcooh .ne. 0 ) B_H(4,t_hcooh) = cab032*cc(i_oh) 
IF ( igcm_hoch2o2 .ne. 0) B_H(4,t_hoch2o2) = cab028 + cab030*(cc(i_hoch2o2) + cc(i_ro2)) &
											+ 0.6*cab044*cc(i_c2h5o2) + 0.6*cab078*cc(i_hcoch2o2) &
											+ cab093*cc(i_hoch2co3) + cab104*cc(i_hcoco3) 
IF ( igcm_hoch2oh .ne. 0 ) B_H(4,t_hoch2oh) = cab035*cc(i_oh) 
IF ( igcm_hoch2ooh .ne. 0 ) B_H(4,t_hoch2ooh) = j(j_hoch2ooh) + j(j_ch3o2h)	
			
dN_dPQ(:) = 0.
dD_dPQ(:) = 0.
dNn_dPQ(:) = 0.
dDd_dPQ(:) = 0.
DO iq = 1, nqmx
	x_j = (iq-1)*nlayermx + lyr_m
	
	dN_dPQ = dN_dPQ + B_h(1,iq)*dccn_dpq( x_j, : )
	dNn_dPQ = dNn_dPQ + B_h(2,iq)*dccn_dpq( x_j, : )
	
	dD_dPQ = dD_dPQ + B_h(3,iq)*dccn_dpq( x_j, : )
	dDd_dPQ = dDd_dPQ + B_h(4,iq)*dccn_dpq( x_j, : )
ENDDO

! d( rh_ho2 )/d( PQ )
! -------------------
drh_ho2 = A_RH(1)*dN_dPQ + A_RH(2)*dNn_dPQ &
		- A_RH(3)*dccn_dpq( (t_ho2-1)*nlayermx + lyr_m, :) &
		- A_RH(4)*dD_dPQ - A_RH(5)*dDd_dPQ &
		+ A_RH(6)*dccn_dpq( (t_h-1)*nlayermx + lyr_m, :)
		
		
! -----------
! 1.2 ROH_HO2 
! -----------
A_roh(1) = 1./roh_ho2_denominator 
A_roh(2) = A_roh(1)/cc_prev(i_ho2) 
A_roh(3) = A_roh(2)*roh_ho2_cabN/cc_prev(i_ho2)
A_roh(4) = roh_ho2*A_roh(1) 
A_roh(5) = A_roh(4)/cc_prev(i_oh)
A_roh(6) =  A_roh(5)*roh_ho2_cabD/cc_prev(i_oh)

! ROH_HO2 : dN_dPQ
B_OH(:,:) = 0. 
B_OH(1,t_o) = c001
B_OH(1,t_o3) = c003*rh_ho2 + c015
B_OH(1,t_h) = 2.*c004 
B_OH(1,t_ho2) = 2.*c008
IF ( igcm_ch3o2 .ne. 0 ) B_OH(1,t_ch3o2) = cab006 + cab007
IF ( igcm_hcho .ne. 0 ) B_OH(1,t_hcho) = cab019  
IF ( igcm_hoch2o2 .ne. 0 ) B_OH(1,t_hoch2o2) = 0.25*cab029

! ROH_HO2 : dNn_dPQ
B_OH(2,t_o) = 0.51*cab002*cc(i_ch4) + cab017*cc(i_ch3o)*0.25 &
			+ cab020*cc(i_hcho) + cab021*cc(i_hco) &
			+ cab037*cc(i_c2h6) 
B_OH(2,t_o1d) = 2.*b002*cc(i_h2o) + b003*cc(i_h2) &
				+ b007*cc(i_ch4)
B_OH(2,t_o2) = cab071*cc(i_ch3co) + cab084*cc(i_hcoco) &
				+ cab092*cc(i_hoch2co)
B_OH(2,t_h2ovap) = 2.*b002*cc(i_o1d) + j(j_h2o)
B_OH(2,t_h2) = b003*cc(i_o1d)			
B_OH(2,t_ch4) = 0.51*cab002*cc(i_o) + b007*cc(i_o1d) 
IF ( igcm_ch3o2 .ne. 0 ) B_OH(2,t_ch3o2) = cab066*cc(i_ch3choho2)
IF ( igcm_ch3ooh .ne. 0 ) B_OH(2,t_ch3ooh) = j(j_ch3o2h)
IF ( igcm_ch3o .ne. 0 ) B_OH(2,t_ch3o) = 0.25*cab017*cc(i_o)
IF ( igcm_hcho .ne. 0 ) B_OH(2,t_hcho) = cab020*cc(i_o)
IF ( igcm_hoch2o2 .ne. 0 ) B_OH(2,t_hoch2o2) = cab066*cc(i_ch3choho2) 
IF ( igcm_hco .ne. 0 ) B_OH(2,t_hco) = cab021*cc(i_o)
! ROH_HO2 : dD_dPQ
B_OH(3,t_co) = e001 
B_OH(3,t_o) = c002 
B_OH(3,t_ho2) = c007
B_OH(3,t_h2o2) = c009 
B_OH(3,t_oh) = 2.*(c013 + c017)
B_OH(3,t_ch4) = cab001 
IF ( igcm_ch3 .ne. 0 ) B_OH(3,t_ch3) = cab107
IF ( igcm_ch3o2 .ne. 0 ) B_OH(3,t_ch3o2) = 2.*cab011
IF ( igcm_ch3ooh .ne. 0 ) B_OH(3,t_ch3ooh) = 0.6*cab014
IF ( igcm_ch3oh .ne. 0 ) B_OH(3,t_ch3oh) = 0.15*cab013 
IF ( igcm_hcho .ne. 0 ) B_OH(3,t_hcho) = cab018 
IF ( igcm_ch3ooh .ne. 0 ) B_OH(3,t_hcooh) = 2.*cab032
IF ( igcm_hoch2oh .ne. 0 ) B_OH(3,t_hoch2oh) = 2.*cab035
IF ( igcm_hoch2ooh .ne. 0 ) B_OH(3,t_hoch2ooh) = cab033 
IF ( igcm_hco .ne. 0 ) B_OH(3,t_hco) = cab025
! ROH_HO2 : dDd_dPQ
B_OH(4,t_o2) = cab015*cc(i_ch3o) + cab026*cc(i_hco) &
			+ cab041*cc(i_c2h5) + cab051*cc(i_hoch2ch2o)
IF ( igcm_ch3o2 .ne. 0 ) B_OH(4,t_ch3o2) = cab030*cc(i_hoch2o2) + 0.6*cab044*cc(i_c2h5o2) &
										+ 0.6*cab078*cc(i_hcoch2o2) + cab093*cc(i_hoch2co3) &
										+ cab104*cc(i_hcoco3)
IF ( igcm_ch3o .ne. 0 ) B_OH(4,t_ch3o) = cab015*cc(i_o2)
IF ( igcm_hoch2o2 .ne. 0 ) B_OH(4,t_hoch2o2) = cab028 + cab030*(cc(i_hoch2o2) + cc(i_ro2)) &
											+ cab030*cc(i_c2h5o2) + 0.6*cab078*cc(i_hcoch2o2) &
											+ cab093*cc(i_hoch2co3) + cab104*cc(i_hcoco3)
IF ( igcm_hco .ne. 0 ) B_OH(4,t_hco) = cab026*cc(i_o2)

dN_dPQ(:) = 0.
dD_dPQ(:) = 0.
dNn_dPQ(:) = 0.
dDd_dPQ(:) = 0.
DO iq = 1, nqmx
	x_j = (iq-1)*nlayermx + lyr_m
	
	dN_dPQ = dN_dPQ + B_OH(1,iq)*dccn_dpq( x_j, : )
	dNn_dPQ = dNn_dPQ + B_OH(2,iq)*dccn_dpq( x_j, : )
	
	dD_dPQ = dD_dPQ + B_OH(3,iq)*dccn_dpq( x_j, : )
	dDd_dPQ = dDd_dPQ + B_OH(4,iq)*dccn_dpq( x_j, : )
ENDDO

! d( roh_ho2 )/d( PQ )
! -------------------
droh_ho2 = A_ROH(1)*dN_dPQ + A_ROH(2)*dNn_dPQ &
		- A_ROH(3)*dccn_dpq( (t_ho2-1)*nlayermx + lyr_m, :) &
		- A_ROH(4)*dD_dPQ - A_ROH(5)*dDd_dPQ &
		+ A_ROH(6)*dccn_dpq( (t_oh-1)*nlayermx + lyr_m, :)
			
droh_ho2 = droh_ho2 + A_ROH(1)*c003*cc(i_o3)*drh_ho2
			

			
			
! ============================================ ! 
! STAGE 2: LINEARISING HYDROGEN ATOMS  
! ============================================ ! 
!
! Equation Coefficients:
GAMMA_H(1) = 1./(1. + (1.+roh_ho2)/rh_ho2) 
GAMMA_H(2) = (GAMMA_H(1)**2.)*cc(i_hox)/rh_ho2
GAMMA_H(3) = (GAMMA_H(1)**2.)*cc(i_hox)*(1.+roh_ho2)/(rh_ho2**2.)

dH_dPQ = GAMMA_H(1)*dHOX_dPQ(lyr_m,:) &
		- GAMMA_H(2)*droh_ho2 &
		+ GAMMA_H(3)*drh_ho2

! ============================================ ! 
! STAGE 3: LINEARISING 	HO2  
! ============================================ ! 
GAMMA_HO2(1) = 1./rh_ho2
GAMMA_HO2(2) = cc(i_h)/(rh_ho2**2)

dHO2_dPQ = GAMMA_HO2(1)*dH_dPQ - GAMMA_HO2(2)*drh_ho2

! ============================================ ! 
! STAGE 4: LINEARISING OH
! ============================================ ! 
GAMMA_OH(1) = roh_ho2
GAMMA_OH(2) = cc(i_ho2)

dOH_dPQ = GAMMA_OH(1)*dHO2_dPQ + GAMMA_OH(2)*droh_ho2

h_j = (t_h-1)*nlayermx + lyr_m
oh_j = (t_oh-1)*nlayermx + lyr_m
ho2_j = (t_ho2-1)*nlayermx + lyr_m

! ============================================ ! 
! STAGE 5: DEPOSIT VALUES
! ============================================ ! 


dccn_dpq( h_j, : ) = dH_dPQ
dccn_dpq( oh_j, : ) = dOH_dPQ
dccn_dpq( ho2_j, : ) = dHO2_dPQ


RETURN


END