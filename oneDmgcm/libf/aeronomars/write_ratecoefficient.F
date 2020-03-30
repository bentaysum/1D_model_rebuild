         subroutine write_ratecoefficient(nesp,rm,dens,j,
     $               a001, a002, a003,
     $               b001, b002, b003, b004, b005, b006,
     $               b007, b008, b009,
     $               c001, c002, c003, c004, c005, c006,
     $               c007, c008, c009, c010, c011, c012,
     $               c013, c014, c015, c016, c017, c018,
     $               d001, d002, d003,
     $               e001, e002, e003,
     $               h001, h002, h003, h004, h005,
     $               t001, t002, t003,
     $               cab001, cab002, cab003, cab004, cab005,
     $               cab006, cab007, cab008, cab009, cab010,
     $               cab011, cab012, cab013, cab014, cab015,
     $               cab016, cab017, cab018, cab019, cab020,
     $               cab021, cab022, cab023, cab024, cab025,
     $               cab026, cab027, cab028, cab029, cab030,
     $               cab031, cab032, cab033, cab034, cab035,
     $               cab036, cab037, cab038, cab039, cab040,
     $               cab041, cab042, cab043, cab044, cab045,
     $               cab046, cab047, cab048, cab049, cab050,
     $               cab051, cab052, cab053, cab054, cab055,
     $               cab056, cab057, cab058, cab059, cab060,
     $               cab061, cab062, cab063, cab064, cab065,
     $               cab066, cab067, cab068, cab069, cab070,
     $               cab071, cab072, cab073, cab074, cab075,
     $               cab076, cab077, cab078, cab079, cab080,
     $               cab081, cab082, cab083, cab084, cab085,
     $               cab086, cab087, cab088, cab089, cab090,
     $               cab091, cab092, cab093, cab094, cab095,
     $               cab096, cab097, cab098, cab099, cab100,
     $               cab101, cab102, cab103, cab104, cab105,
     $               cab106, cab107)
	  
       IMPLICIT NONE

#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h" 
#include "callkeys.h" 
       
cccccc Input variables
      integer nesp   ! Number of species 
      real rm(nlayermx,nesp) ! Mixing Ratios
      real dens(nlayermx) ! Number density
      real j(nlayermx,nd) ! Photolytic frequencies
      
      real a001(nlayermx), a002(nlayermx), a003(nlayermx)
      real b001(nlayermx), b002(nlayermx), b003(nlayermx), 
     $     b004(nlayermx), b005(nlayermx), b006(nlayermx),
     $     b007(nlayermx), b008(nlayermx), b009(nlayermx)
      real c001(nlayermx), c002(nlayermx), c003(nlayermx), 
     $     c004(nlayermx), c005(nlayermx), c006(nlayermx), 
     $     c007(nlayermx), c008(nlayermx), c009(nlayermx), 
     $     c010(nlayermx), c011(nlayermx), c012(nlayermx),
     $     c013(nlayermx), c014(nlayermx), c015(nlayermx), 
     $     c016(nlayermx), c017(nlayermx), c018(nlayermx)
      real d001(nlayermx), d002(nlayermx), d003(nlayermx)
      real e001(nlayermx), e002(nlayermx), e003(nlayermx)
      real h001(nlayermx), h002(nlayermx), h003(nlayermx),
     $     h004(nlayermx), h005(nlayermx)
      real t001(nlayermx), t002(nlayermx), t003(nlayermx)

      real cab001(nlayermx), cab002(nlayermx), cab003(nlayermx),
     $     cab004(nlayermx), cab005(nlayermx), cab006(nlayermx),
     $     cab007(nlayermx), cab008(nlayermx), cab009(nlayermx),
     $     cab010(nlayermx), cab011(nlayermx), cab012(nlayermx),
     $     cab013(nlayermx), cab014(nlayermx), cab015(nlayermx),
     $     cab016(nlayermx), cab017(nlayermx), cab018(nlayermx), 
     $     cab019(nlayermx), cab020(nlayermx), cab021(nlayermx), 
     $     cab022(nlayermx), cab023(nlayermx), cab024(nlayermx), 
     $     cab025(nlayermx), cab026(nlayermx), cab027(nlayermx), 
     $     cab028(nlayermx), cab029(nlayermx), cab030(nlayermx),
     $     cab031(nlayermx), cab032(nlayermx), cab033(nlayermx),
     $     cab034(nlayermx), cab035(nlayermx), cab036(nlayermx),
     $     cab037(nlayermx), cab038(nlayermx), cab039(nlayermx),
     $     cab040(nlayermx), cab041(nlayermx), cab042(nlayermx),
     $     cab043(nlayermx), cab044(nlayermx), cab045(nlayermx),
     $     cab046(nlayermx), cab047(nlayermx), cab048(nlayermx),
     $     cab049(nlayermx), cab050(nlayermx), cab051(nlayermx),
     $     cab052(nlayermx), cab053(nlayermx), cab054(nlayermx),
     $     cab055(nlayermx), cab056(nlayermx), cab057(nlayermx),
     $     cab058(nlayermx), cab059(nlayermx), cab060(nlayermx),
     $     cab061(nlayermx), cab062(nlayermx), cab063(nlayermx),
     $     cab064(nlayermx), cab065(nlayermx), cab066(nlayermx),
     $     cab067(nlayermx), cab068(nlayermx), cab069(nlayermx),
     $     cab070(nlayermx), cab071(nlayermx), cab072(nlayermx),
     $     cab073(nlayermx), cab074(nlayermx), cab075(nlayermx),
     $     cab076(nlayermx), cab077(nlayermx), cab078(nlayermx),
     $     cab079(nlayermx), cab080(nlayermx), cab081(nlayermx),
     $     cab082(nlayermx), cab083(nlayermx), cab084(nlayermx),
     $     cab085(nlayermx), cab086(nlayermx), cab087(nlayermx),
     $     cab088(nlayermx), cab089(nlayermx), cab090(nlayermx),
     $     cab091(nlayermx), cab092(nlayermx), cab093(nlayermx),
     $     cab094(nlayermx), cab095(nlayermx), cab096(nlayermx),
     $     cab097(nlayermx), cab098(nlayermx), cab099(nlayermx),
     $     cab100(nlayermx), cab101(nlayermx), cab102(nlayermx),
     $     cab103(nlayermx), cab104(nlayermx), cab105(nlayermx),
     $     cab106(nlayermx), cab107(nlayermx)

     
ccccc Local variables    

c     tracer indexes in the chemistry:

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


      integer j_o2_o, j_o2_o1d, j_co2_o, j_co2_o1d,
     $        j_o3_o1d, j_o3_o, j_h2o, j_hdo, j_h2o2,
     $        j_ho2, j_no2, j_ch4_ch3_h, j_ch4_1ch2_h2,
     $        j_ch4_3ch2_h_h, j_ch4_ch_h2_h, j_ch3o2h,
     $        j_ch2o_co, j_ch2o_hco, j_ch3oh, j_c2h6, j_hcl,
     $        j_hocl, j_clo, j_so2, j_so, j_h2s, j_so3,
     $        j_hno3, j_hno4,
     $        j_ch3cho_ch3, j_ch3cho_ch4, j_ch3cho_h, 
     $        j_hoch2ooh, j_hoch2cho_hco, j_hoch2cho_co,
     $        j_hoch2cho_oh, j_glyox_hco, j_glyox_hcho,
     $        j_glyox_h2, j_ch3cooh, j_ch3coooh, j_ch3cocooh
c
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     numbering of photolysis rates
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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
      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c New rates : Tropospheric Photochemistry of Ozone, its
c  			  Precursors and the Hydroxyl Radical: A 3D
c 			  modeling Study Considering Non-methane 
C 			  Hydrocarbon - Rolf von Kuhlmann [2001]
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc =
C c                              ! ch3coooh + hv -> ch3 + oh + co2 [ = (hcho + hv -> hco + h) x 0.025] (1)
C                                ! hoch2ooh + hv -> hcooh + oh + ho2 [ = (ch3ooh + hv) ] (2)
C                                ! hoch2co3h + hv -> hcho + ho2 + oh + co2 [ = (ch3ooh + hv) ] (2)
C                                ! hcoco2h + hv -> 2.*ho2 + co + co2 [ = (hcho + hv -> co + h2)*3.95 ] (1)
C                                ! ch3chohooh + hv -> ch3 + hcooh + oh [ = (ch3ooh + hv) ] (2)
C                                ! hyetho2h + hv -> hoch2ch2o + oh [ = (ch3ooh + hv) ] (2)
                                 ! hcoco3h + hv -> ho2 + co + oh + co2 [ = j_ch3o2h + SUM(j_hoch2cho)] (2)
                                 ! hooch2cho + hv -> oh + hcho + co + ho2 [= j_ch3o2h + SUM(j_hoch2cho)] (2)
                                 
                                 
!=================!
! Non-Organics    !
C !=================! 

c              CALL WRITEDIAGFI(1, 'a001',
c     &     "o + o2 + co2 -> o3 + co2", 'molecules/s',1,
c     &      a001)
     

C              CALL WRITEDIAGFI(1, 'a002', 
C      &     "o + o + co2 -> o2 + co2", 'molecules/s',1,
C      &     a002) 
    
C              CALL WRITEDIAGFI(1, 'a003',
C      &     "o + o3 -> o2 + o2", 'molecules/s',1,
C      &     a003)
     
C              CALL WRITEDIAGFI(1,'b001',
C      &     "o(1d) + co2 -> o + co2", 'molecules/s',1,
C      &     b001) 
             
C              CALL WRITEDIAGFI(1,'b002',
C      &     "o(1d) + h2o -> oh + oh", 'molecules/s',1,
C      &     b002)
     
C              CALL WRITEDIAGFI(1,'b003',
C      &     "o(1d) + h2 -> oh + h", 'molecules/s',1,
C      &     b003)
            
C              CALL WRITEDIAGFI(1,'b004',
C      &     "o(1d) + o2 -> o + o2", 'molecules/s',1,
C      &     b004)
            
C              CALL WRITEDIAGFI(1,'b005',
C      &     "o(1d) + o3 -> o2 + o2", 'molecules/s',1,
C      &     b005) 
     
C              CALL WRITEDIAGFI(1,'b006',
C      &     "o(1d) + o3 -> o2 + o + o", 'molecules/s',1,
C      &     b006*rm(:,i_o1d)*rm(:,i_o3)*dens*dens)      
 
             CALL WRITEDIAGFI(1, 'b009', 
     &     "o(1d) + ch4 -> hcho + h2", 'molecules/s',1,
     &      b009)

             CALL WRITEDIAGFI(1, 'b008', 
     &     "o(1d) + ch4 -> ch3o + h", 'molecules/s',1,
     &      b008)

             CALL WRITEDIAGFI(1, 'b007', 
     &     "o(1d) + ch4 -> ch3 + oh", 'molecules/s',1,
     &      b007)

C              CALL WRITEDIAGFI(1, 'c001',
C      &     "o + ho2 -> oh + o2", 'molecules/s',1,
C      &     c001*rm(:,i_o)*rm(:,i_ho2)*dens*dens)
         
C              CALL WRITEDIAGFI(1, 'c002',
C      &     "o + oh -> o2 + h", 'molecules/s',1,
C      &     c002*rm(:,i_o)*rm(:,i_oh)*dens*dens)

C              CALL WRITEDIAGFI(1, 'c003',
C      &     "h + o3 -> oh + o2", 'molecules/s',1,
C      &     c003*rm(:,i_h)*rm(:,i_o3)*dens*dens)

C              CALL WRITEDIAGFI(1, 'c004',
C      &     "h + ho2 -> oh + oh", 'molecules/s',1,
C      &     c004*rm(:,i_h)*rm(:,i_ho2)*dens*dens)

C              CALL WRITEDIAGFI(1, 'c005',
C      &     "h + ho2 -> h2 + o2", 'molecules/s',1,
C      &     c005*rm(:,i_h)*rm(:,i_ho2)*dens*dens)
    
C              CALL WRITEDIAGFI(1, 'c006',
C      &     "h + ho2 -> h2o + o", 'molecules/s',1,
C      &     c006*rm(:,i_h)*rm(:,i_ho2)*dens*dens)

C              CALL WRITEDIAGFI(1, 'c007',
C      &     "oh + ho2 -> h2o + o2", 'molecules/s',1,
C      &     c007*rm(:,i_oh)*rm(:,i_ho2)*dens*dens)     
     
C              CALL WRITEDIAGFI(1, 'c008',
C      &     "ho2 + ho2 -> h2o2 + o2", 'molecules/s',1,
C      &     c008*rm(:,i_ho2)*rm(:,i_ho2)*dens*dens)          

C              CALL WRITEDIAGFI(1, 'c009',
C      &     "oh + h2o2 -> h2o + ho2", 'molecules/s',1,
C      &     c009*rm(:,i_oh)*rm(:,i_h2o2)*dens*dens)        

C              CALL WRITEDIAGFI(1, 'c010',
C      &     "oh + h2 -> h2o + h", 'molecules/s',1,
C      &     c010*rm(:,i_oh)*rm(:,i_h2)*dens*dens)        
     
C              CALL WRITEDIAGFI(1, 'c011',
C      &     "h + o2 + M -> ho2 + M", 'molecules/s',1,
C      &     c011*rm(:,i_h)*rm(:,i_o2)*dens*dens)           

C              CALL WRITEDIAGFI(1, 'c012',
C      &     "o + h2o2 -> oh + ho2", 'molecules/s',1,
C      &     c012*rm(:,i_o)*rm(:,i_h2o2)*dens*dens)     
     
C              CALL WRITEDIAGFI(1, 'c013',
C      &     "oh + oh -> h2o + o", 'molecules/s',1,
C      &     c013*rm(:,i_oh)*rm(:,i_oh)*dens*dens)     
     
C              CALL WRITEDIAGFI(1, 'c014',
C      &     "oh + o3 -> ho2 + o2", 'molecules/s',1,
C      &     c014*rm(:,i_oh)*rm(:,i_o3)*dens*dens)      
     
C              CALL WRITEDIAGFI(1, 'c015',
C      &     "ho2 + o3 -> oh + o2 + o2", 'molecules/s',1,
C      &     c015*rm(:,i_ho2)*rm(:,i_o3)*dens*dens)      

C              CALL WRITEDIAGFI(1, 'c016',
C      &     "ho2 + ho2 + M -> h2o2 + o2 + M", 'molecules/s',1,
C      &     c016*rm(:,i_ho2)*rm(:,i_ho2)*dens*dens)      

C              CALL WRITEDIAGFI(1, 'c017',
C      &     "oh + oh + M -> h2o2 + M", 'molecules/s',1,
C      &     c017*rm(:,i_oh)*rm(:,i_oh)*dens*dens)      
     
C              CALL WRITEDIAGFI(1, 'c018',
C      &     "h + h + M  -> h2 + M", 'molecules/s',1,
C      &     c018*rm(:,i_h)*rm(:,i_h)*dens*dens)      
     
C              CALL WRITEDIAGFI(1, 'e001',
C      &     "oh + co  -> co2 + h", 'molecules/s',1,
C      &     e001*rm(:,i_oh)*rm(:,i_co)*dens*dens)     

C              CALL WRITEDIAGFI(1, 'e002',
C      &     "o + co + M -> co2 + M", 'molecules/s',1,
C      &     e002*rm(:,i_o)*rm(:,i_co)*dens*dens)          

!========================!
! Organic Rates       
!========================!

             CALL WRITEDIAGFI(1, 'cab001',
     &     "ch4 + oh -> ch3 + oh", 'molecules/s',1,
     &     cab001)          

             CALL WRITEDIAGFI(1, 'cab002',
     &     "ch4 + o -> 0.51*ch3 + 0.51*oh + 0.49*ch3o + 0.49*h ", 
     &     'molecules/s',1,
     &      cab002)       
     
             CALL WRITEDIAGFI(1, 'cab003',
     &     "ch3 + o2 + M -> ch3o2 + M ", 'molecules/s',1,
     &      cab003)          
     
             CALL WRITEDIAGFI(1, 'cab004',
     &     "ch3 + o3 -> 0.956*hcho + 0.956*h + 0.044*ch3o + o2 ",
     &     'molecules/s',1,
     &      cab004)          
     
             CALL WRITEDIAGFI(1, 'cab005',
     &     "ch3 + o -> 0.83*hcho + 0.83*h + 0.17*co + 0.17*h2 + 0.17*h", 
     &     'molecules/s',1,
     &      cab005)   
     
             CALL WRITEDIAGFI(1, 'cab006',
     &     " ch3o2 + ho2 -> ch3ooh + o2", 'molecules/s',1,
     &      cab006)
     
             CALL WRITEDIAGFI(1, 'cab007',
     &     " ch3o2 + ho2 -> hcho + h2o + o2 ", 'molecules/s',1,
     &      cab007)  

             CALL WRITEDIAGFI(1, 'cab008',
     &     "ch3o2 + ro2 -> ch3o + 0.5*o2", 'molecules/s',1,
     &      cab008)  
     
             CALL WRITEDIAGFI(1, 'cab009',
     &     " ch3o2 + ro2 -> 0.5*hcho + 0.5ch3oh + 0.5*o2  ",
     &     'molecules/s',1,
     &      cab009)  
     
             CALL WRITEDIAGFI(1, 'cab010',
     &     "ch3o2 + o3 -> ch3o + 2.*o2 ", 'molecules/s',1,
     &      cab010)  
     
             CALL WRITEDIAGFI(1, 'cab011',
     &     "ch3o2 + oh -> ch3o + ho2", 'molecules/s',1,
     &     cab011)  
     
             CALL WRITEDIAGFI(1, 'cab012',
     &     "ch3o2 + o -> ch3o + o2 ", 'molecules/s',1,
     &     cab012)  
     
             CALL WRITEDIAGFI(1, 'cab013',
     &     " ch3oh + oh -> 0.85*hcho + 0.85*ho2 + 0.15*ch3o + h2o ",
     &     'molecules/s',1,
     &     cab013)  
     
             CALL WRITEDIAGFI(1, 'cab014',
     &     "ch3ooh + oh -> 0.6*ch3o2 + 0.4*hcho + 0.4*oh + h2o ",
     &     'molecules/s',1,
     &     cab014)  
     
             CALL WRITEDIAGFI(1, 'cab015',
     &     "ch3o + o2 -> ho2 + hcho", 'molecules/s',1,
     &     cab015)  
     
             CALL WRITEDIAGFI(1, 'cab016',
     &     "ch3o + o3 -> ch3o2 + o2", 'molecules/s',1,
     &     cab016)          

             CALL WRITEDIAGFI(1, 'cab017',
     &     "ch3o + o -> 0.75*ch3 + 0.75*o2 + 0.25*hcho + 0.25*oh ", 
     &     'molecules/s',1,
     &      cab017)       
     
             CALL WRITEDIAGFI(1, 'cab018',
     &     "hcho + oh -> hco + h2o", 'molecules/s',1,
     &      cab018)          
     
             CALL WRITEDIAGFI(1, 'cab019',
     &     "hcho + ho2 -> hoch2o2 ",
     &     'molecules/s',1,
     &      cab019)          
     
             CALL WRITEDIAGFI(1, 'cab020',
     &     "hcho + o -> hco + oh ", 
     &     'molecules/s',1,
     &      cab020)   
     
             CALL WRITEDIAGFI(1, 'cab021',
     &     "hco + o -> co + oh ", 'molecules/s',1,
     &      cab021)
     
             CALL WRITEDIAGFI(1, 'cab022',
     &     "hco + ch3 -> ch4 + co", 'molecules/s',1,
     &      cab022)  

             CALL WRITEDIAGFI(1, 'cab023',
     &     "hco + ch3 -> ch3cho", 'molecules/s',1,
     &      cab023)  
     
             CALL WRITEDIAGFI(1, 'cab024',
     &     " 2.*hco -> hcho + co",
     &     'molecules/s',1,
     &      cab024)  
     
             CALL WRITEDIAGFI(1, 'cab025',
     &     " hco + oh -> co + h2o", 'molecules/s',1,
     &      cab025)  
     
             CALL WRITEDIAGFI(1, 'cab026',
     &     "hco + o2 -> co + ho2", 'molecules/s',1,
     &     cab026)  
     
             CALL WRITEDIAGFI(1, 'cab027',
     &     "hco + h -> co + h2", 'molecules/s',1,
     &     cab027)  
     
             CALL WRITEDIAGFI(1, 'cab028',
     &     "hoch2o2 -> hcho + ho2 ",
     &     'molecules/s',1,
     &     cab028)  
     
             CALL WRITEDIAGFI(1, 'cab029',
     &     " hoch2o2 + ho2 -> 0.5*hoch2ooh + 0.5*hcooh + 0.2*oh +"
     &     //" 0.2*ho2 + 0.3*h2o + 0.8*o2 ",
     &     'molecules/s',1,
     &     cab029)  
     
             CALL WRITEDIAGFI(1, 'cab030',
     &     "hoch2o2 + ro2 -> hcooh + ho2", 'molecules/s',1,
     &     cab030)       
     
             CALL WRITEDIAGFI(1, 'cab031',
     &     "hoch2o2 + ro2 -> 0.5*hcooh + 0.5*hoch2oh + 0.5*o2",
     &     'molecules/s',1,
     &     cab031)       
     
             CALL WRITEDIAGFI(1, 'cab032',
     &     "hcooh + oh  -> co2 + ho2 + h2o ", 'molecules/s',1,
     &     cab032)       
     
             CALL WRITEDIAGFI(1, 'cab033',
     &     "hoch2ooh + oh -> hoch2o2 ", 'molecules/s',1,
     &     cab033)       
     
             CALL WRITEDIAGFI(1, 'cab034',
     &     "hoch2ooh + oh -> hcooh + h2o + oh", 'molecules/s',1,
     &     cab034)       
     
             CALL WRITEDIAGFI(1, 'cab035',
     &     "hoch2oh + oh -> ho2 + hcooh + h2o", 'molecules/s',1,
     &     cab035)       
     
             CALL WRITEDIAGFI(1, 'cab036',
     &     "c2h6 + oh -> c2h5 + h2o", 'molecules/s',1,
     &     cab036)       
     
             CALL WRITEDIAGFI(1, 'cab037',
     &     "c2h6 + o -> c2h5 + oh ", 'molecules/s',1,
     &     cab037)       
     
             CALL WRITEDIAGFI(1, 'cab038',
     &     " ch3 + ch3 + M -> c2h6", 'molecules/s',1,
     &     cab038)       
     
             CALL WRITEDIAGFI(1, 'cab039',
     &     "c2h5 + o2 + M -> c2h5o2 + M ", 'molecules/s',1,
     &     cab039)       
     
             CALL WRITEDIAGFI(1, 'cab040',
     &     "c2h5 + c2h5 -> c2h6 + c2h4", 'molecules/s',1,
     &     cab040)       
     
             CALL WRITEDIAGFI(1, 'cab041',
     &     "c2h5 + o2 -> c2h4 + ho2", 'molecules/s',1,
     &     cab041)       
     
             CALL WRITEDIAGFI(1, 'cab042',
     &     "c2h5 + h -> ch3 + ch3", 'molecules/s',1,
     &     cab042)       
     
             CALL WRITEDIAGFI(1, 'cab043',
     &     "c2h5o2 + ho2 -> c2h5ooh + o2", 'molecules/s',1,
     &     cab043)       
     
             CALL WRITEDIAGFI(1, 'cab044',
     &     "c2h5o2 + ro2 -> 0.8*ch3cho + 0.6*ho2 + 0.2*c2h5oh + o2",
     &     'molecules/s',1,
     &     cab044)       
     
             CALL WRITEDIAGFI(1, 'cab045',
     &     "c2h5ooh + oh -> c2h5o2 + h2o", 'molecules/s',1,
     &     cab045)       
     
             CALL WRITEDIAGFI(1, 'cab046',
     &     "c2h5ooh + oh -> ch3cho + oh", 'molecules/s',1,
     &     cab046)       
     
             CALL WRITEDIAGFI(1, 'cab047',
     &     "c2h5oh + oh -> 095*c2h5o2 + 0.95*ho2 + "
     &     //"0.05*hoch2ch2o2 + h2o", 'molecules/s',1,
     &     cab047)       
     
             CALL WRITEDIAGFI(1, 'cab048',
     &     "hoch2ch2o2 + ro2 -> 0.6*hoch2ch2o + 0.2*hoch2cho +"
     &     //" 0.2*ethgly", 'molecules/s',1,
     &     cab048)       
     
             CALL WRITEDIAGFI(1, 'cab049',
     &     "hoch2ch2o2 + ho2 -> hyetho2h", 'molecules/s',1,
     &     cab049)       
     
             CALL WRITEDIAGFI(1, 'cab050',
     &     "hoch2ch2o2 + ho2 -> hoch2ch2o + oh ", 'molecules/s',1,
     &     cab050)       

             CALL WRITEDIAGFI(1, 'cab051',
     &     "hoch2ch2o + o2 -> ho2 + hoch2cho ", 'molecules/s',1,
     &     cab051)       
     
             CALL WRITEDIAGFI(1, 'cab052',
     &     "hoch2ch2o -> ho2 + hcho + hcho ", 'molecules/s',1,
     &     cab052)       
     
             CALL WRITEDIAGFI(1, 'cab053',
     &     "ethgly + oh -> hoch2cho + ho2 + h2o ", 'molecules/s',1,
     &     cab053)       
     
             CALL WRITEDIAGFI(1, 'cab054',
     &     "hyetho2h + oh -> hoch2ch2o2 + h2o ", 'molecules/s',1,
     &     cab054)       
     
             CALL WRITEDIAGFI(1, 'cab055',
     &     "hyetho2h + oh -> hoch2cho + oh + h2o", 'molecules/s',1,
     &     cab055)       
     
             CALL WRITEDIAGFI(1, 'cab056',
     &     "hyetho2h + oh -> hooch2cho + ho2 + h2o ", 'molecules/s',1,
     &     cab056)       
     
             CALL WRITEDIAGFI(1, 'cab057',
     &     "ch3cho + oh -> ch3co + h2o ", 'molecules/s',1,
     &     cab057)       
     
             CALL WRITEDIAGFI(1, 'cab058',
     &     "ch3cho + oh -> hcoch2o2 + h2o ", 'molecules/s',1,
     &     cab058)       
     
             CALL WRITEDIAGFI(1, 'cab059',
     &     "ch3cho + ho2 -> ch3choho2 ", 'molecules/s',1,
     &     cab059)       
     
             CALL WRITEDIAGFI(1, 'cab060',
     &     "ch3cho + hcooh -> ch2choh + hcooh  ", 'molecules/s',1,
     &     cab060)       
     
             CALL WRITEDIAGFI(1, 'cab061',
     &     " ch2choh + oh -> hcooh + oh + hcho ", 'molecules/s',1,
     &     cab061)       
     
             CALL WRITEDIAGFI(1, 'cab062',
     &     "ch2choh + oh -> hoch2cho + ho2 ", 'molecules/s',1,
     &     cab062)       
     
             CALL WRITEDIAGFI(1, 'cab063',
     &     "ch2choh + hcooh -> ch3cho + hcooh ", 'molecules/s',1,
     &     cab063)       
     
             CALL WRITEDIAGFI(1, 'cab064',
     &     "ch3choho2 -> ch3cho + ho2 ", 'molecules/s',1,
     &     cab064)       
     
             CALL WRITEDIAGFI(1, 'cab065',
     &     "ch3choho2 + ho2 -> 0.5*ch3chohooh + 0.3*ch3cooh + 0.2*ch3"
     &     // " + 0.2*hcooh + 0.2*oh + o2", 'molecules/s',1,
     &     cab065)       
     
             CALL WRITEDIAGFI(1, 'cab066',
     &     " ch3choho2 + ro2 -> ch3 + hcooh + oh ", 'molecules/s',1,
     &     cab066)       
     
             CALL WRITEDIAGFI(1, 'cab067',
     &     "ch3cooh + oh -> ch3 + co2 + h2o", 'molecules/s',1,
     &     cab067)       
     
             CALL WRITEDIAGFI(1, 'cab068',
     &     "ch3chohooh + oh -> ch3cooh + oh ", 'molecules/s',1,
     &     cab068)       
     
             CALL WRITEDIAGFI(1, 'cab069',
     &     "ch3chohooh + oh -> ch3choho2 ", 'molecules/s',1,
     &     cab069)       

             CALL WRITEDIAGFI(1, 'cab070',
     &     "ch3c(o) + o2 + M -> ch3c(o)oo + M ", 'molecules/s',1,
     &     cab070)       

             CALL WRITEDIAGFI(1, 'cab071',
     &     "ch3c(o) + o2 + M -> oh + hcho + co + M ", 'molecules/s',1,
     &     cab071)       

             CALL WRITEDIAGFI(1, 'cab072',
     &     "ch3c(o)oo + ho2 -> oh + ch3 + co2 ", 'molecules/s',1,
     &     cab072)       

             CALL WRITEDIAGFI(1, 'cab073',
     &     "ch3c(o)oo + ho2 -> ch3c(o)ooh", 'molecules/s',1,
     &     cab073)       

             CALL WRITEDIAGFI(1, 'cab074',
     &     "ch3c(o)oo + ho2 -> ch3cooh + o3 ", 'molecules/s',1,
     &     cab074)       

             CALL WRITEDIAGFI(1, 'cab075',
     &     "ch3c(o)oo + ro2 -> ch3 + co2 ", 'molecules/s',1,
     &     cab075)       

             CALL WRITEDIAGFI(1, 'cab076',
     &     "ch3c(o)oo + ro2 -> ch3cooh ", 'molecules/s',1,
     &     cab076)       

             CALL WRITEDIAGFI(1, 'cab077',
     &     "ch3c(o)ooh + oh -> ch3c(o)oo + h2o ", 'molecules/s',1,
     &     cab077)       

             CALL WRITEDIAGFI(1, 'cab078',
     &     "hcoch2o2 + ro2 -> 0.6*hcho + 0.6*co + 0.6*ho2 + 0.*glyox "
     &     // "+ 0.2*hoch2cho ", 'molecules/s',1,
     &     cab078)       

             CALL WRITEDIAGFI(1, 'cab079',
     &     "hcoch2o2 + ho2 -> hooch2cho ", 'molecules/s',1,
     &     cab079)       

             CALL WRITEDIAGFI(1, 'cab080',
     &     "hcoch2o2 + ho2 -> hcho + co + ho2 + oh", 'molecules/s',1,
     &     cab080)       

             CALL WRITEDIAGFI(1, 'cab081',
     &     "glyox + oh -> hcoco + h2o ", 'molecules/s',1,
     &     cab081)       

             CALL WRITEDIAGFI(1, 'cab082',
     &     "hcoco -> hco + co", 'molecules/s',1,
     &     cab082)       

             CALL WRITEDIAGFI(1, 'cab083',
     &     "hcoco + o2 -> hcoco3 ", 'molecules/s',1,
     &     cab083)       

             CALL WRITEDIAGFI(1, 'cab084',
     &     "hcoco + o2 -> oh + co + co2 ", 'molecules/s',1,
     &     cab084)       

             CALL WRITEDIAGFI(1, 'cab085',
     &     "hooch2cho + oh -> hcoch2o2 + h2o", 'molecules/s',1,
     &     cab085)       

             CALL WRITEDIAGFI(1, 'cab086',
     &     "hooch2cho + oh -> hcho + co + oh + h2o ", 'molecules/s',1,
     &     cab086)       

             CALL WRITEDIAGFI(1, 'cab087',
     &     " hooch2cho + oh -> glyox + oh + h2o ", 'molecules/s',1,
     &     cab087)       
     
             CALL WRITEDIAGFI(1, 'cab088',
     &     "hoch2cho + oh -> hoch2co + h2o ", 'molecules/s',1,
     &     cab088)       

             CALL WRITEDIAGFI(1, 'cab089',
     &     " hoch2cho + oh -> hochcho + h2o ", 'molecules/s',1,
     &     cab089)       

             CALL WRITEDIAGFI(1, 'cab090',
     &     "hochcho -> glyox + ho2  ", 'molecules/s',1,
     &     cab090)       

             CALL WRITEDIAGFI(1, 'cab091',
     &     " hoch2co + o2 + M -> hoch2co3  ", 'molecules/s',1,
     &     cab091)       

             CALL WRITEDIAGFI(1, 'cab092',
     &     "hoch2co + o2 + M -> oh + hcho + co2 + M ", 'molecules/s',1,
     &     cab092)       

             CALL WRITEDIAGFI(1, 'cab093',
     &     "hoch2co3 + ro2 -> hcho + co2 + ho2 ", 'molecules/s',1,
     &     cab093)       

             CALL WRITEDIAGFI(1, 'cab094',
     &     "hoch2co3 + ro2 -> hoch2co2h ", 'molecules/s',1,
     &     cab094)       

             CALL WRITEDIAGFI(1, 'cab095',
     &     "hoch2co3 + ho2 -> hcho + oh + ho2 + co2", 'molecules/s',1,
     &     cab095)       

             CALL WRITEDIAGFI(1, 'cab096',
     &     "hoch2co3 + ho2 -> hoch2co3h", 'molecules/s',1,
     &     cab096)       

             CALL WRITEDIAGFI(1, 'cab097',
     &     "hoch2co3 + ho2 -> hoch2co2h + o3", 'molecules/s',1,
     &     cab097)       

             CALL WRITEDIAGFI(1, 'cab098',
     &     "hoch2co2h + oh -> 0.09*hcho + 0.9*co2 + 0.91*hcoco2h + ho2"
     &     // " + h2o ", 'molecules/s',1,
     &     cab098)       

             CALL WRITEDIAGFI(1, 'cab099',
     &     "hcoco2h + oh -> co + ho2 + co2 + h2o ", 'molecules/s',1,
     &     cab099)       

             CALL WRITEDIAGFI(1, 'cab100',
     &     "hoch2co3h + oh -> hoch2co3 + h2o ", 'molecules/s',1,
     &     cab100)       

             CALL WRITEDIAGFI(1, 'cab101',
     &     "hoch2co3h + oh -> hcoco3h + ho2 ", 'molecules/s',1,
     &     cab101)       

            CALL WRITEDIAGFI(1, 'cab102',
     &     "hcoco3h + oh -> hcoco3 + h2o ", 'molecules/s',1,
     &     cab102)       

             CALL WRITEDIAGFI(1, 'cab103',
     &     "hcoco3h + oh -> co + co2 + h2o + oh ", 'molecules/s',1,
     &     cab103)       

             CALL WRITEDIAGFI(1, 'cab104',
     &     "hcoco3 + ro2 -> co + ho2 + co2 ", 'molecules/s',1,
     &     cab104)       

             CALL WRITEDIAGFI(1, 'cab105',
     &     "hcoco3 + ro2 -> hcoco2h + o2 ", 'molecules/s',1,
     &     cab105)       

             CALL WRITEDIAGFI(1, 'cab106',
     &     "hcoco3 + ho2 -> ho2 + co + co2 + oh ", 'molecules/s',1,
     &     cab106)       

              CALL WRITEDIAGFI(1, 'cab107',
     &     "ch3 + oh + M -> ch3oh + M ", 'molecules/s',1,
     &     cab107)       
ccccccccccccc 
cc Photolysis 
C ccccccccccccc
   
C           CALL WRITEDIAGFI(1, 'j_o2_o',
C      &     "o2 + hv -> 2*o ", 'molecules/s',1,
C      &     j(:,j_o2_o))      

C           CALL WRITEDIAGFI(1, 'j_o2_o1d',
C      &     "o2 + hv -> o + o(1d) ", 'molecules/s',1,
C      &     j(:,j_o2_o1d)*rm(:,i_o2)*dens)
	 
C           CALL WRITEDIAGFI(1, 'j_co2_o',
C      &     "co2 + hv -> co + o", 'molecules/s',1,
c     &     j(:,j_co2_o)*rm(:,i_co2)*dens)
    
c          CALL WRITEDIAGFI(1, 'j_co2_o1d',
c     &     "co2 + hv -> co + o(1d) ", 'molecules/s',1,
c     &     j(:,j_co2_o1d)*rm(:,i_co2)*dens)

c         CALL WRITEDIAGFI(1, 'j_o3_o',
c     &     "o3 + hv -> o2 + o ", 'molecules/s',1,
c     &     j(:,j_o3_o)*rm(:,i_o3)*dens)
c    
c          CALL WRITEDIAGFI(1, 'j_o3_o1d',
c     &     "o3 + hv -> o2 + o(1d) ", 'molecules/s',1,
c     &     j(:,j_o3_o1d)*rm(:,i_o3)*dens)
c   
c          CALL WRITEDIAGFI(1, 'j_h2o',
c     &     "h2o + hv -> h + oh ", 'molecules/s',1,
c     &     j(:,j_h2o)*rm(:,i_h2o)*dens)
c   
c          CALL WRITEDIAGFI(1, 'j_h2o2',
c     &     "h2o2 + hv -> 2*oh ", 'molecules/s',1,
c     &     j(:,j_h2o2)*rm(:,i_h2o2)*dens)
c          CALL WRITEDIAGFI(1, 'j_ho2',
c     &     "ho2 + hv -> oh + o", 'molecules/s',1,
c     &     j(:,j_ho2)*rm(:,i_ho2)*dens)
    
             CALL WRITEDIAGFI(1, 'j_ch3cho_ch3',
     &     "ch3cho + hv -> ch3 + hco ", 'molecules/s',1,
     &     j(:,j_ch3cho_ch3))       
   
             CALL WRITEDIAGFI(1, 'j_ch3cho_ch4',
     &     "ch3cho + hv -> ch4 + co ", 'molecules/s',1,
     &     j(:,j_ch3cho_ch4))       
   
             CALL WRITEDIAGFI(1, 'j_ch3cho_h',
     &     "ch3cho + hv -> ch3co + h ", 'molecules/s',1,
     &     j(:,j_ch3cho_h))       
     
             CALL WRITEDIAGFI(1, 'j_hoch2ooh',
     &     "hoch2ooh + o2 + hv -> hcooh + ho2 + oh", 'molecules/s',1,
     &     j(:,j_hoch2ooh)*rm(:,i_o2)*dens*3.5e-14)     

             CALL WRITEDIAGFI(1, 'j_hoch2cho_hco',
     &     "hoch2cho + hv -> ch3o + hco ", 'molecules/s',1,
     &     j(:,j_hoch2cho_hco))     
     
             CALL WRITEDIAGFI(1, 'j_hoch2cho_co',
     &     "hoch2cho + hv -> ch3oh + co", 'molecules/s',1,
     &     j(:,j_hoch2cho_co))     

             CALL WRITEDIAGFI(1, 'j_hoch2cho_oh',
     &     "hoch2cho + hv -> hcoch2o2 + oh", 'molecules/s',1,
     &     j(:,j_hoch2cho_oh))     
     
             CALL WRITEDIAGFI(1, 'j_glyox_hco',
     &     "glyox + hv -> hco + hco ", 'molecules/s',1,
     &     j(:,j_glyox_hco))     
     
              CALL WRITEDIAGFI(1, 'j_glyox_h2',
     &     "glyox + hv -> h2 + 2co", 'molecules/s',1,
     &     j(:,j_glyox_h2))     
     
              CALL WRITEDIAGFI(1, 'j_glyox_hcho',
     &     "glyox + hv -> hcho + co", 'molecules/s',1,
     &     j(:,j_glyox_hcho))     
     
             CALL WRITEDIAGFI(1, 'j_ch3cooh',
     &     "ch3cooh + hv -> ch3 + cooh", 'molecules/s',1,
     &     j(:,j_ch3cooh))     
     
             CALL WRITEDIAGFI(1, 'j_ch3coooh',
     &     "ch3c(o)ooh + hv -> ch3 + oh + co2", 'molecules/s',1,
     &     j(:,j_ch3coooh))     
     
             CALL WRITEDIAGFI(1, 'j_ch3coooh',
     &     "ch3coco(oh) + hv -> products ", 'molecules/s',1,
     &     j(:,j_ch3cocooh))
           
             CALL WRITEDIAGFI(1, 'j_ch4_ch3_h',
     &     "ch4 + hv -> ch3 + h", 'molecules/s',1,
     &     j(:,j_ch4_ch3_h))
           
             CALL WRITEDIAGFI(1, 'j_ch4_1ch2_h2',
     &     "ch4 + hv -> 1ch2 + h2", 'molecules/s',1,
     &     j(:,j_ch4_1ch2_h2))
           
             CALL WRITEDIAGFI(1, 'j_ch4_3ch2_h_h',
     &     "ch4 + hv -> 3ch2 + h + h", 'molecules/s',1,
     &     j(:,j_ch4_3ch2_h_h))
           
             CALL WRITEDIAGFI(1, 'j_ch4_ch_h2_h',
     &     "ch4 + hv -> ch + h2 + h", 'molecules/s',1,
     &     j(:,j_ch4_ch_h2_h))
           
             CALL WRITEDIAGFI(1, 'j_ch3o2h',
     &     "ch3ooh + hv -> ch3o + oh", 'molecules/s',1,
     &     j(:,j_ch3o2h))
           
             CALL WRITEDIAGFI(1, 'j_ch2o_hco',
     &     "hcho + hv -> h + hco", 'molecules/s',1,
     &     j(:,j_ch2o_hco))
           
             CALL WRITEDIAGFI(1, 'j_ch2o_co',
     &     "hcho + hv -> h2 + co", 'molecules/s',1,
     &     j(:,j_ch2o_co))
           
             CALL WRITEDIAGFI(1, 'j_ch3oh',
     &     "ch3oh + hv -> ch3o + h", 'molecules/s',1,
     &     j(:,j_ch3oh))
     
             CALL WRITEDIAGFI(1, 'j_c2h6',
     &     "c2h6 + hv -> products", 'molecules/s',1,
     &     j(:,j_c2h6))
C 	 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c New rates : Tropospheric Photochemistry of Ozone, its
c  			  Precursors and the Hydroxyl Radical: A 3D
c 			  modeling Study Considering Non-methane 
C 			  Hydrocarbon - Rolf von Kuhlmann [2001]
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc =
C                                ! hoch2ooh + hv -> hcooh + oh + ho2 [ = (ch3ooh + hv) ] (2)
C                                ! hoch2co3h + hv -> hcho + ho2 + oh + co2 [ = (ch3ooh + hv) ] (2)
C                                ! hcoco2h + hv -> 2.*ho2 + co + co2 [ = (hcho + hv -> co + h2)*3.95 ] (1)
C                                ! ch3chohooh + hv -> ch3 + hcooh + oh [ = (ch3ooh + hv) ] (2)
C                                ! hyetho2h + hv -> hoch2ch2o + oh [ = (ch3ooh + hv) ] (2)
                                 ! hcoco3h + hv -> ho2 + co + oh + co2 [ = j_ch3o2h + SUM(j_hoch2cho)] (2)
                                 ! hooch2cho + hv -> oh + hcho + co + ho2 [= j_ch3o2h + SUM(j_hoch2cho)] (2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c New rates : Tropospheric Photochemistry of Ozone, its
c  			  Precursors and the Hydroxyl Radical: A 3D
c 			  modeling Study Considering Non-methane 
C 			  Hydrocarbon - Rolf von Kuhlmann [2001]
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc =
C                                ! hoch2co3h + hv -> hcho + ho2 + oh + co2 [ = (ch3ooh + hv) ] (2)
C                                ! hcoco2h + hv -> 2.*ho2 + co + co2 [ = (hcho + hv -> co + h2)*3.95 ] (1)
C                                ! ch3chohooh + hv -> ch3 + hcooh + oh [ = (ch3ooh + hv) ] (2)
C                                ! hyetho2h + hv -> hoch2ch2o + oh [ = (ch3ooh + hv) ] (2)
                                 ! hcoco3h + hv -> ho2 + co + oh + co2 [ = j_ch3o2h + SUM(j_hoch2cho)] (2)
                                 ! hooch2cho + hv -> oh + hcho + co + ho2 [= j_ch3o2h + SUM(j_hoch2cho)] (2)
    
      return
      end