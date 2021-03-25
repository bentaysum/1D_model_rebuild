SUBROUTINE tlm_clox(lyr_m, rclo_cl, iter, dt, dens, &
				pcl, lcl, pclo, lclo, &
				cc, cc_prev, &
				nesp, &
                    dccn_dpq, dcc0_dpq, &
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
                    ,cl052, cl053, cl054 &
                    ,cl055, cl056, cl057 &
                    ,cl058, cl059)

USE TLMvars

IMPLICIT NONE

#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "callkeys.h"
#include "conc.h"


! Input Variables
! ===============
INTEGER lyr_m
REAL rclo_cl ! Chlorine Partition function 
integer iter
REAL dt ! chemistry timestep
REAL dens ! Atmospheric Number Density 
REAL pcl, lcl, pclo, lclo ! Cl and ClO Production and loss
REAL cc(nesp), cc_prev(nesp)
INTEGER nesp
REAL dccn_dpq(nqmx*nlayermx,nqmx*nlayermx)
REAL, INTENT(IN) :: dcc0_dpq(nqmx*nlayermx,nqmx*nlayermx)

! Local Variables 
! ===============
REAL A(6) ! Partition Function Linearising Coefficients

REAL Acl(2), Aclo(2) 

REAL dPcl_dPQ(nqmx*nlayermx), dPclo_dPQ(nqmx*nlayermx) ! Linearised Production 
REAL dLcl_dPQ(nqmx*nlayermx), dLclo_dPQ(nqmx*nlayermx) ! Linearised Loss
 
REAL pcl_coeff(nqmx), pclo_coeff(nqmx) ! Linearising Coefficients for Production
REAL lcl_coeff(nqmx), lclo_coeff(nqmx) ! Linearising Coefficients for Loss 

REAL drclo_dpq(nqmx*nlayermx) ! Linearised partition function 

INTEGER iq, l ! Tracer iterator
INTEGER x_j, X_I

REAL dCl_dPQ(nqmx*nlayermx), dClO_dPQ(nqmx*nlayermx) ! Linearised Number Densities 

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


! Photolysis indexing 
! -------------------
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

! Tracer numbering in model 
! -------------------------
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

logical :: exist

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


! ===================================================
! Tangent Linear Model : ClOx species
! ---------------------------------------------------
! 
! Main Equations:
!
!   [Cl] = [ClOx]/(1. + rclo_cl)
!
!   [ClO] = [Cl] x rclo_cl
!
! where ClOx, the sum total of Cl and ClO, is handled
! via the SIBEM equation, calculated in the TLM_sibem
! file.
!
!=====================================================
! 1.0 : Linearising Partition Function rclo_cl
! ----------------------------------------------------
!
! rclo_cl = (pclo + lcl)/(pcl + lclo)
!
! where p and l terms are the loss and production rates,
! both in units of s^-1, of Cl and ClO within the rapid
! Cl <--> ClO conversion network.
!
! Linearised equation:
!
! rclo_cl ' =  (pcl + lclo)^-1 x ( pclo' + lcl' )
!           -  (pclo + lcl)x(pcl + lclo)^-2 x ( pcl' + lclo')
! 
! ------------------------------------------------------
! 1.1 : Coefficients
! ------------------------------------------------------
    A(1) = 1./(pcl + lclo)
    A(2) = A(1)*rclo_cl

! ------------------------------------------------------
! 1.2: Linearised Production Rates 
! ------------------------------------------------------
    pcl_coeff(:) = 0.
    pclo_coeff(:) = 0. 

    ! 1.2.1 : Cl 
    ! ----------
    pcl_coeff(t_o) = cl002 
    pcl_coeff(t_clo) = 2.*cl004 + cl005 
    pcl_coeff(t_oh) = 0.94*cl012 

    ! 1.2.2 : ClO 
    ! -----------
    pclo_coeff(t_o3) = cl001 
    pclo_coeff(t_ho2) = cl010
    pclo_coeff(t_cloo) = cl024*2. 
    pclo_coeff(t_clo4) = cl053
    pclo_coeff(t_hocl) = cl056
    pclo_coeff(t_oclo) = 2.*cl059 

    IF ( igcm_ch3 .ne. 0 ) THEN 
        pclo_coeff(t_ch3o2) = cl022*0.5
    ENDIF

! ------------------------------------------------------
! 1.3: Linearised Loss Rate Coefficients
! ------------------------------------------------------
    lcl_coeff(:) = 0.
    lclo_coeff(:) = 0. 

    ! 1.3.1 : Cl  
    ! ----------
    lcl_coeff(t_o3) = cl001 
    lcl_coeff(t_ho2) = cl010 
    lcl_coeff(t_cloo) = cl024
    lcl_coeff(t_clo4) = cl053 
    lcl_coeff(t_hocl) = cl056
    lcl_coeff(t_oclo) = cl059

    ! 1.3.2 : ClO
    lclo_coeff(t_o) = cl002 
    lclo_coeff(t_clo) = cl004 + cl005 
    lclo_coeff(t_oh) = cl012 


    IF ( igcm_ch3 .ne. 0 ) THEN 
        lcl_coeff(t_ch3o2) = 0.5*cl022 
    ENDIF

! ------------------------------------------------------
! 1.4 : Calculating the Linearised P and L terms 
! ------------------------------------------------------
    dPcl_dPQ(:) = 0.
    dLcl_dPQ(:) = 0.
    dPclo_dPQ(:) = 0. 
    dLclo_dPQ(:) = 0. 

    DO iq = 1,nqmx

          x_j = (iq-1)*nlayermx + lyr_m

          ! Cl 
          dPcl_dPQ = dPcl_dPQ + pcl_coeff(iq)*dccn_dpq(x_j,:)
          dLcl_dPQ = dLcl_dPQ + lcl_coeff(iq)*dccn_dpq(x_j,:)
          ! ClO
          dPclo_dPQ = dPclo_dPQ + pclo_coeff(iq)*dccn_dpq(x_j,:) 
          dLclo_dPQ = dLclo_dPQ + lclo_coeff(iq)*dccn_dpq(x_j,:) 

    ENDDO

! ------------------------------------------------------
! 1.5 : Calculating the Linearised Partition Function 
! ------------------------------------------------------
    drclo_dpq = A(1)*( dPclo_dPQ + dLcl_dPQ ) &
              - A(2)*( dLclo_dPQ + dPcl_dPQ )

! =======================================================
! 2.0 : Linearised Cl and ClO Number Densities
! -------------------------------------------------------
!
! [Cl] ' = (1. + rclo_cl)^-1 x [ClOx]' 
!        - [ClOx] x (1. + rclo_cl)^-2 x [rclo_cl]'
!
! [ClO]' = rclo_cl x [Cl]' 
!        + [Cl] x [rclo_cl]'
! --------------------------------------------------------
! 2.1 : Coefficients
! --------------------------------------------------------
    Acl(1) = 1./(1.+rclo_cl)
    Acl(2) = cc(i_clox)*(Acl(1)**2.) 

    Aclo(1) = rclo_cl
    Aclo(2) = cc(i_cl)

! --------------------------------------------------------
! 2.2 : Calculations
! --------------------------------------------------------
    
    ! 2.2.1 : Cl 
    ! ----------
    dCl_dPQ = Acl(1)*dClOx_dPQ( lyr_m, : ) &
            - Acl(2)*drclo_dpq 

    ! 2.2.2 : ClO 
    ! -----------
    dClO_dPQ = Aclo(1)*dCl_dPQ &
             + Aclo(2)*drclo_dpq

! --------------------------------------------------------
! 2.3 : Insertion into Linearised NUmber Density Array
! --------------------------------------------------------
    ! 2.3.1: Cl 
    x_j = (t_cl-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = dCl_dPQ 

    ! 2.3.2: ClO 
    x_j = (t_clo-1)*nlayermx + lyr_m 
    dccn_dpq( x_j, : ) = dClO_dPQ 


    RETURN 

END SUBROUTINE
