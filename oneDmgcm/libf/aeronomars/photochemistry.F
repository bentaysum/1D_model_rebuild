c*****************************************************************
c
c     Photochemical routine 
c
c     Author: Franck Lefevre
c     ------
c
c     Version: 17/03/2011
c
c*****************************************************************
c
      subroutine photochemistry(lswitch, zycol, sza, ptimestep, 
     $                  press,temp, dens, dist_sol, surfdust1d,
     $                  surfice1d, jo3, tau,pq,pdq)
          
      USE TLMvars 
      implicit none
c
#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h" 
#include "callkeys.h"
#include "tracer.h"
#include "conc.h"

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     input/output:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      real zycol(nlayermx,nqmx)  ! chemical species volume mixing ratio
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      real sza                   ! solar zenith angle (deg)
      real ptimestep             ! physics timestep (s)
      real press(nlayermx)       ! pressure (hpa)
      real temp(nlayermx)        ! temperature (k)
      real dens(nlayermx)        ! density (cm-3)
      real dist_sol              ! sun distance (au) 
      real surfdust1d(nlayermx)  ! dust surface area (cm2/cm3)
      real surfice1d(nlayermx)   ! ice surface area (cm2/cm3)
      real tau                   ! optical depth at 7 hpa
      real pq(nlayermx,nqmx)   ! mixing ratios (mmrs)
      real pdq(nlayermx,nqmx)  ! tracer tendency prior to photochem routine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     output:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      real jo3(nlayermx)  ! photodissociation rate o3 -> o1d
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer phychemrat         ! ratio physics timestep/chemistry timestep
      integer istep
      integer i_co2,i_o3,j_o3_o1d,lev
      integer nesp
      integer lswitch
c
      integer l, iq
      parameter (nesp = 73)      ! number of species in the chemistry (+3 for RO2, e-, and dust)
c
      real stimestep             ! standard timestep for the chemistry (s)
      real ctimestep             ! real timestep for the chemistry (s) 
      real rm(nlayermx,nesp)     ! species volume mixing ratio 
      real j(nlayermx,nd)        ! interpolated photolysis rates (s-1)
      real rmco2(nlayermx)       ! co2 mixing ratio
      real rmo3(nlayermx)        ! ozone mixing ratio
c
c     reaction rates
c
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
 
cc    reaction rates from CAABA/MECCA 
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

cc    Chlorine Reaction Scheme
      real cl001(nlayermx), cl002(nlayermx), cl003(nlayermx) 
     $ ,cl004(nlayermx), cl005(nlayermx), cl006(nlayermx) 
     $ ,cl007(nlayermx), cl008(nlayermx), cl009(nlayermx) 
     $ ,cl010(nlayermx), cl011(nlayermx), cl012(nlayermx) 
     $ ,cl013(nlayermx), cl014(nlayermx), cl015(nlayermx) 
     $ ,cl016(nlayermx), cl017(nlayermx), cl018(nlayermx) 
     $ ,cl019(nlayermx), cl020(nlayermx), cl021(nlayermx) 
     $ ,cl022(nlayermx), cl023(nlayermx), cl024(nlayermx) 
     $ ,cl025(nlayermx), cl026(nlayermx), cl027(nlayermx) 
     $ ,cl028(nlayermx), cl029(nlayermx), cl030(nlayermx) 
     $ ,cl031(nlayermx), cl032(nlayermx), cl033(nlayermx) 
     $ ,cl034(nlayermx), cl035(nlayermx), cl036(nlayermx)
     $ ,cl037(nlayermx), cl038(nlayermx), cl039(nlayermx)
     $ ,cl040(nlayermx), cl041(nlayermx), cl042(nlayermx)
     $ ,cl043(nlayermx), cl044(nlayermx), cl045(nlayermx)
     $ ,cl046(nlayermx), cl047(nlayermx), cl048(nlayermx)
     $ ,cl049(nlayermx), cl050(nlayermx), cl051(nlayermx)
     $ ,cl052(nlayermx), cl053(nlayermx), cl054(nlayermx)
     $ ,cl055(nlayermx), cl056(nlayermx), cl057(nlayermx)
     $ ,cl058(nlayermx), cl059(nlayermx)


      REAL Ncl(nlayermx)

c     Tangent Linear Model variables
      real dccn_dpq(nqmx*nlayermx,nqmx*nlayermx)
      real, save :: dcc0_dpq(nqmx*nlayermx,nqmx*nlayermx)
      real dcc0_dpq_firstcall(nqmx*nlayermx,nlayermx*nqmx)
      real dHOX_dPQ(nlayermx,nlayermx*nqmx)
      real dHOX0_dPQ(nlayermx,nlayermx*nqmx)
      real dOX_dPQ(nlayermx,nlayermx*nqmx)
      real dOX0_dPQ(nlayermx,nlayermx*nqmx)

      REAL lifetime(nlayermx,nesp)

      integer h_j, oh_j, ho2_j
      integer o_j, o3_j 
      integer x_j 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     stimestep  : standard timestep for the chemistry (s)       c
c     ctimestep  : real timestep for the chemistry (s)           c
c     phychemrat : ratio physical/chemical timestep              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      stimestep = 600.


      phychemrat = nint(ptimestep/stimestep)
c
      ctimestep = ptimestep/real(phychemrat)

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialisation of chemical species and families            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c


      call gcmtochim(zycol, lswitch, nesp, rm) 


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute reaction rates                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  


      call chemrates(lswitch, dens, press, temp,
     $               surfdust1d, surfice1d, 
     $               a001, a002, a003,
     $               b001, b002, b003, b004, b005, b006,
     $               b007, b008, b009,
     $               c001, c002, c003, c004, c005, c006,
     $               c007, c008, c009, c010, c011, c012,
     $               c013, c014, c015, c016, c017, c018,
     $               d001, d002, d003,
     $               e001, e002, e003,
     $               h001, h002, h003, h004, h005,
     $               t001, t002, t003, tau,
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
     $               cab106, cab107,
     $               cl001, cl002, cl003,
     $               cl004, cl005, cl006,
     $               cl007, cl008, cl009,
     $               cl010, cl011, cl012,
     $               cl013, cl014, cl015,
     $               cl016, cl017, cl018,
     $               cl019, cl020, cl021,
     $               cl022, cl023, cl024,
     $               cl025, cl026, cl027,
     $               cl028, cl029, cl030,
     $               cl031, cl032, cl033,
     $               cl034, cl035, cl036,
     $               cl037, cl038, cl039,
     $               cl040, cl041, cl042,
     $               cl043, cl044, cl045, cl046, cl047, cl048,
     $               cl049, cl050, cl051, cl052, cl053, cl054,
     $               cl055, cl056, cl057, cl058, cl059)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     interpolation of photolysis rates in the lookup table      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      i_co2 = 1
      i_o3  = 6 
c

      do lev = 1,lswitch-1
         rmco2(lev) = rm(lev,i_co2)
         rmo3(lev)  = rm(lev, i_o3)
      end do
c


      call phot(lswitch, press, temp, sza, tau, dist_sol, 
     $          rmco2, rmo3, j)
      
      j_o3_o1d = 5
      do lev = 1,lswitch-1
         jo3(lev) = j(lev,j_o3_o1d)
      end do
c

c     -----------------------------------
c     Tangent Linear Model Initialisation 
c     -----------------------------------
C       IF ( tlm_on ) THEN 

C         CALL TLM_chem(dcc0_dpq, dccn_dpq, dHOX_dPQ, dOX_dPQ,
C      $                dHOX0_dPQ, dOX0_dPQ,
C      $                ptimestep, istep, phychemrat,
C      $                dens, sza, .True.)

C       ENDIF ! tlm_on 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over time within the physical timestep                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    


      do istep = 1,phychemrat

c =================================================================

         call chimie(lswitch,nesp, rm, j, dens, ctimestep, ptimestep,
     $               press, temp, sza, dist_sol, surfdust1d,
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
     $               lifetime,
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
     $               cab106, cab107,
     $               cl001, cl002, cl003,
     $               cl004, cl005, cl006,
     $               cl007, cl008, cl009,
     $               cl010, cl011, cl012,
     $               cl013, cl014, cl015,
     $               cl016, cl017, cl018,
     $               cl019, cl020, cl021,
     $               cl022, cl023, cl024,
     $               cl025, cl026, cl027,
     $               cl028, cl029, cl030,
     $               cl031, cl032, cl033,
     $               cl034, cl035, cl036,
     $               cl037, cl038, cl039,
     $               cl040, cl041, cl042,
     $               cl043, cl044, cl045, cl046, cl047, cl048,
     $               cl049, cl050, cl051, cl052, cl053, cl054,
     $               cl055, cl056, cl057, cl058, cl059,
     $               dccn_dpq, dcc0_dpq,
     $               dHOX_dPQ, dHOX0_dPQ,
     $               dOX_dPQ, dOX0_dPQ)

c     --------------------
c     Tangent Linear Model
c     --------------------
C       IF ( tlm_on ) THEN 

C         CALL TLM_chem(dcc0_dpq, dccn_dpq, dHOX_dPQ, dOX_dPQ,
C      $                dHOX0_dPQ, dOX0_dPQ,
C      $                ptimestep, istep, phychemrat,
C      $                dens, sza, .False.)

C       ENDIF ! tlm_on 
      
      end do
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     save chemical species for the gcm                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 

      call chimtogcm(zycol, lswitch, nesp, rm) 


      ! Chlorine Species Lifetimes 
      ! ==========================
       call WRITEDIAGFI(1,"cl_tau","cl_tau", "s",
     $                  1,lifetime(:,60))
       call WRITEDIAGFI(1,"clo_tau","clo_tau", "s",
     $                  1,lifetime(:,61))
       call WRITEDIAGFI(1,"cl2_tau","cl2_tau", "s",
     $                  1,lifetime(:,62))
       call WRITEDIAGFI(1,"oclo_tau","oclo_tau", "s",
     $                  1,lifetime(:,63))
       call WRITEDIAGFI(1,"cl2o2_tau","cl2o2_tau", "s",
     $                  1,lifetime(:,64))
       call WRITEDIAGFI(1,"hcl_tau","hcl_tau", "s",
     $                  1,lifetime(:,65))
       call WRITEDIAGFI(1,"hocl_tau","hocl_tau", "s",
     $                  1,lifetime(:,66))
       call WRITEDIAGFI(1,"cloo_tau","cloo_tau", "s",
     $                  1,lifetime(:,67))
       call WRITEDIAGFI(1,"ch3ocl_tau","ch3ocl_tau", "s",
     $                  1,lifetime(:,68))
       call WRITEDIAGFI(1,"clco_tau","clco_tau", "s",
     $                  1,lifetime(:,69))
       call WRITEDIAGFI(1,"clo3_tau","clo3_tau", "s",
     $                  1,lifetime(:,70))
       call WRITEDIAGFI(1,"hclo4_tau","hclo4_tau", "s",
     $                  1,lifetime(:,71))
       call WRITEDIAGFI(1,"clo4_tau","clo4_tau", "s",
     $                  1,lifetime(:,72))


      ! HCl Chemical Reactions
      ! ----------------------
      ! Production Rates 
       call WRITEDIAGFI(1,"hcl_P","hcl_P", "cm-3 s-1",
     $                  1,lifetime(:,1))

       call WRITEDIAGFI(1, "hcl_p1",
     $                  "$Cl + H_2 \longrightarrow HCl + H$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,2))

       call WRITEDIAGFI(1, "hcl_p2",
     $                  "$Cl + HO_2 \longrightarrow HCl + O_2$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,3))

       call WRITEDIAGFI(1, "hcl_p3",
     $                  "$Cl + H_2O_2 \longrightarrow HCl + HO_2$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,4))

       call WRITEDIAGFI(1, "hcl_p4",
     $                  "$ClO + OH \longrightarrow HCl + O_2$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,5))

       call WRITEDIAGFI(1, "hcl_p5",
     $                  "$Cl + CH_4 \longrightarrow HCl + CH_3$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,6))

       call WRITEDIAGFI(1, "hcl_p6",
     $                  "$Cl + HCHO \longrightarrow HCl + HCO$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,7))

       call WRITEDIAGFI(1, "hcl_p7",
     $                  "$Cl + CH_3OOH \longrightarrow HCl + HCHO + OH$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,8))

       call WRITEDIAGFI(1, "hcl_p8",
     $                  "$Cl + CH_3O_2 \longrightarrow HCl + CH_2OO"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,9))

       call WRITEDIAGFI(1, "hcl_p9",
     $                  "$Cl + CH_3OCl \longrightarrow HCl + CH_2OCl$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,10))

       call WRITEDIAGFI(1, "hcl_p10",
     $                  "$Cl_2 + H \longrightarrow HCl + Cl$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,11))

       call WRITEDIAGFI(1, "hcl_p11",
     $            "$Cl + HCOOH \longrightarrow HCl + \textrm{Products}$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,12))

       call WRITEDIAGFI(1, "hcl_p12",
     $                  "$Cl + HOCl \longrightarrow HCl + ClO$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,13))

C        ! Loss Rates
       call WRITEDIAGFI(1, "hcl_l1",
     $                  "$HCl + h\mu \longrightarrow H + Cl$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,14))

       call WRITEDIAGFI(1, "hcl_l2",
     $                  "$HCl + OH \longrightarrow H_2O + Cl$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,15))

       call WRITEDIAGFI(1, "hcl_l3",
     $  "$HCl + O(^1D) \longrightarrow 0.33*H + 0.33*ClO +"
     $  // " 0.66*Cl + 0.66*OH$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,16))

       call WRITEDIAGFI(1, "hcl_l4",
     $                  "$HCl + O(^3P) \longrightarrow OH + Cl$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,17))

       call WRITEDIAGFI(1, "hcl_l5",
     $                  "$HCl + H \longrightarrow H_2 + Cl$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,18))

      ! -----------------------------
      ! Organic Compound Lifeitmes 
      ! -----------------------------
       call WRITEDIAGFI(1,"ch4_tau","ch4_tau", "s",
     $                  1,lifetime(:,19))

       call WRITEDIAGFI(1,"hcho_tau","hcho_tau", "s",
     $                  1,lifetime(:,20))

       call WRITEDIAGFI(1,"hcooh_tau","hcooh_tau", "s",
     $                  1,lifetime(:,21))

       call WRITEDIAGFI(1,"ch3ooh_tau","ch3ooh_tau", "s",
     $                  1,lifetime(:,22))


      ! Perchlorate Production
       call WRITEDIAGFI(1,"hclo4_P","hclo4_P", "cm-3 s-1",
     $                  1,lifetime(:,23))

       call WRITEDIAGFI(1, "hclo4_p1",
     $                  "$ClO_3 + OH \longrightarrow HClO_4$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,24))

       call WRITEDIAGFI(1, "hclo4_p2",
     $                  "$ClO_3 + OH + M\longrightarrow HClO_4 + M$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,25))

       call WRITEDIAGFI(1, "hclo4_p3",
     $                  "$ClO_4 + HOCl \longrightarrow HClO_4 + OH$"
     $                  , "cm-3 s-1",
     $                  1, lifetime(:,26))

       ! O=O Lifetime
       call WRITEDIAGFI(1,"o=o_tau","o=o_tau", "s",
     $                  1,lifetime(:,60))


       ! ClOX 
C        call WRITEDIAGFI(1,"clox_P","clox_P", "s",
C      $                  1,lifetime(:,27))
C        call WRITEDIAGFI(1,"clox_L","clox_L", "s",
C      $                  1,lifetime(:,28))
C        call WRITEDIAGFI(1,"pcl","pcl", "s",
C      $                  1,lifetime(:,29))
C        call WRITEDIAGFI(1,"lcl","lcl", "s",
C      $                  1,lifetime(:,30))
C        call WRITEDIAGFI(1,"pclo","pclo", "s",
C      $                  1,lifetime(:,31))
C        call WRITEDIAGFI(1,"lclo","lclo", "s",
C      $                  1,lifetime(:,32))

       ! PLACE IN WRITEDIAGFI FUNCTION 
C        call WRITEDIAGFI(1,"clox_l1","cl003", "s",
C      $                  1,lifetime(:,33))
C        call WRITEDIAGFI(1,"clox_l2","cl005", "s",
C      $                  1,lifetime(:,34))
C        call WRITEDIAGFI(1,"clox_l3","cl006", "s",
C      $                  1,lifetime(:,35))
C        call WRITEDIAGFI(1,"clox_l4","cl008", "s",
C      $                  1,lifetime(:,36))
C        call WRITEDIAGFI(1,"clox_l5","cl009", "s",
C      $                  1,lifetime(:,37))
C        call WRITEDIAGFI(1,"clox_l6","cl011", "s",
C      $                  1,lifetime(:,38))
C        call WRITEDIAGFI(1,"clox_l7","cl012", "s",
C      $                  1,lifetime(:,39))
C        call WRITEDIAGFI(1,"clox_l8","cl013", "s",
C      $                  1,lifetime(:,40))
C        call WRITEDIAGFI(1,"clox_l9","cl016", "s",
C      $                  1,lifetime(:,41))
C        call WRITEDIAGFI(1,"clox_l10","cl017", "s",
C      $                  1,lifetime(:,42))
C        call WRITEDIAGFI(1,"clox_l11","cl018", "s",
C      $                  1,lifetime(:,43))
C        call WRITEDIAGFI(1,"clox_l12","cl019", "s",
C      $                  1,lifetime(:,44))
C        call WRITEDIAGFI(1,"clox_l13","cl020", "s",
C      $                  1,lifetime(:,45))
C        call WRITEDIAGFI(1,"clox_l14","cl021", "s",
C      $                  1,lifetime(:,46))
C        call WRITEDIAGFI(1,"clox_l15","cl022", "s",
C      $                  1,lifetime(:,47))
C        call WRITEDIAGFI(1,"clox_l16","cl023", "s",
C      $                  1,lifetime(:,48))
C        call WRITEDIAGFI(1,"clox_l17","cl025", "s",
C      $                  1,lifetime(:,49))
C        call WRITEDIAGFI(1,"clox_l18","cl028", "s",
C      $                  1,lifetime(:,50))
C        call WRITEDIAGFI(1,"clox_l19","cl029", "s",
C      $                  1,lifetime(:,51))
C        call WRITEDIAGFI(1,"clox_l20","cl030", "s",
C      $                  1,lifetime(:,52))
C        call WRITEDIAGFI(1,"clox_l21","cl031", "s",
C      $                  1,lifetime(:,53))
C        call WRITEDIAGFI(1,"clox_l22","cl043", "s",
C      $                  1,lifetime(:,54))
C        call WRITEDIAGFI(1,"clox_l23","cl044", "s",
C      $                  1,lifetime(:,55))
C        call WRITEDIAGFI(1,"clox_l24","cl045", "s",
C      $                  1,lifetime(:,56))
C        call WRITEDIAGFI(1,"clox_l25","cl046", "s",
C      $                  1,lifetime(:,57))
C        call WRITEDIAGFI(1,"clox_l26","cl047", "s",
C      $                  1,lifetime(:,58))
C        call WRITEDIAGFI(1,"clox_l27","cl055", "s",
C      $                  1,lifetime(:,59))







      return
      end
c
c*****************************************************************
c
      subroutine chimie(lswitch, nesp, rm, j, dens, dt, ptimestep,
     $                  press, t, sza, dist_sol, surfdust,
     $                  a001, a002, a003,
     $                  b001, b002, b003, b004, b005, b006,
     $                  b007, b008, b009,
     $                  c001, c002, c003, c004, c005, c006,
     $                  c007, c008, c009, c010, c011, c012,
     $                  c013, c014, c015, c016, c017, c018,
     $                  d001, d002, d003,
     $                  e001, e002, e003,
     $                  h001, h002, h003, h004, h005,
     $                  t001, t002, t003,
     $                   lifetime,
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
     $               cab106, cab107,
     $               cl001, cl002, cl003,
     $               cl004, cl005, cl006,
     $               cl007, cl008, cl009,
     $               cl010, cl011, cl012,
     $               cl013, cl014, cl015,
     $               cl016, cl017, cl018,
     $               cl019, cl020, cl021,
     $               cl022, cl023, cl024,
     $               cl025, cl026, cl027,
     $               cl028, cl029, cl030,
     $               cl031, cl032, cl033,
     $               cl034, cl035, cl036,
     $               cl037, cl038, cl039,
     $               cl040, cl041, cl042,
     $               cl043, cl044, cl045, cl046, cl047, cl048,
     $               cl049, cl050, cl051, cl052, cl053, cl054,
     $               cl055, cl056, cl057, cl058, cl059,
     $               dccn_dpq, dcc0_dpq,
     $               dHOX_dPQ, dHOX0_dPQ,
     $               dOX_dPQ, dOX0_dPQ)
c
c*****************************************************************
c
      USE TLMvars 
      
      implicit none 
c
#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h" 
#include "callkeys.h" 
#include "tracer.h" 	
#include "conc.h"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     input/output:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      integer lswitch         ! interface level between chemistries
      integer nesp            ! number of species
c0.!
      real rm(nlayermx,nesp)  ! volume mixing ratios 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real dens(nlayermx)     ! density (cm-3)
      real dt                 ! chemistry timestep (s) 
      real ptimestep          ! physical timestep (s)
      real j(nlayermx,nd)     ! interpolated photolysis rates (s-1) 
      real press(nlayermx)    ! pressure (hpa)
      real t(nlayermx)        ! temperature (k)
      real sza                ! solar zenith angle (deg)
      real dist_sol           ! sun distance (au) 
      real surfdust(nlayermx) ! Dust surface area cm2 / cm3 of air
      real pq(nlayermx,nqmx)  ! mixing ratios (mmr)
      real pdq(nlayermx,nqmx) ! tracer tendency prior to photochem routine
c     reaction rates
c
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

      real cl001(nlayermx), cl002(nlayermx), cl003(nlayermx) 
     $ ,cl004(nlayermx), cl005(nlayermx), cl006(nlayermx) 
     $ ,cl007(nlayermx), cl008(nlayermx), cl009(nlayermx) 
     $ ,cl010(nlayermx), cl011(nlayermx), cl012(nlayermx) 
     $ ,cl013(nlayermx), cl014(nlayermx), cl015(nlayermx) 
     $ ,cl016(nlayermx), cl017(nlayermx), cl018(nlayermx) 
     $ ,cl019(nlayermx), cl020(nlayermx), cl021(nlayermx) 
     $ ,cl022(nlayermx), cl023(nlayermx), cl024(nlayermx) 
     $ ,cl025(nlayermx), cl026(nlayermx), cl027(nlayermx) 
     $ ,cl028(nlayermx), cl029(nlayermx), cl030(nlayermx) 
     $ ,cl031(nlayermx), cl032(nlayermx), cl033(nlayermx)
     $ ,cl034(nlayermx), cl035(nlayermx), cl036(nlayermx)
     $ ,cl037(nlayermx), cl038(nlayermx), cl039(nlayermx)
     $ ,cl040(nlayermx), cl041(nlayermx), cl042(nlayermx)
     $ ,cl043(nlayermx), cl044(nlayermx), cl045(nlayermx)
     $ ,cl046(nlayermx), cl047(nlayermx), cl048(nlayermx)
     $ ,cl049(nlayermx), cl050(nlayermx), cl051(nlayermx)
     $ ,cl052(nlayermx), cl053(nlayermx), cl054(nlayermx)
     $ ,cl055(nlayermx), cl056(nlayermx), cl057(nlayermx)
     $ ,cl058(nlayermx), cl059(nlayermx)


c     Tangent Linear Model 
      real dccn_dpq(nqmx*nlayermx,nqmx*nlayermx)
      real, INTENT(IN) :: dcc0_dpq(nqmx*nlayermx,nqmx*nlayermx)
      real dHOX_dPQ(nlayermx,nqmx*nlayermx)
      real dHOX0_dPQ(nlayermx,nqmx*nlayermx)
      real dOX_dPQ(nlayermx,nqmx*nlayermx)
      real dOX0_dPQ(nlayermx,nqmx*nlayermx)

      REAL lifetime(nlayermx,nesp)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real hetero_ice, hetero_dust
c
      integer iesp, iter, l, niter, i, k

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
     $        j_glyox_h2, j_ch3cooh, j_ch3coooh, j_ch3cocooh,
     $        j_cl2, j_cloo, j_cl2o2, j_oclo


c
      parameter (hetero_ice  = 0.)!1.)   ! switch for het. chem. on ice clouds
      parameter (hetero_dust = 0.)   ! switch for het. chem. on dust
                                     ! hetero_dust = 0. advised for the moment
c
      parameter (niter = 5)          ! iterations in the chemical scheme [originally 5]
c
      real cc0(nlayermx,nesp)        ! initial number density (cm-3)
      real cc(nlayermx,nesp)         ! final number density (cm-3)
      real nox(nlayermx)             ! nox number density (cm-3)
      real no(nlayermx)              ! no  number density (cm-3)
      real no2(nlayermx)             ! no2 number density (cm-3)
      real ro2(nlayermx)             ! R-O2 number density (cm-3)
      real production(nlayermx,nesp) ! production rate
      real loss(nlayermx,nesp)       ! loss rate

      real ro_o3, rh_ho2, roh_ho2, rno2_no

      integer iq
      real cc_prev(nesp)

      real cc_ox_next    ! Next timestep cc(i_ox) value for TLM 
      real cc_hox_next   ! " " cc(i_hox) " "
       
c     Rate for CH4 + X -> CH3 + ... proposed reaction 
c     ===============================================
      real methane_enhancement
      parameter( methane_enhancement = 1.E6 )
      real k_pseudo(nlayermx)
      
      real o2_loss_rate(nlayermx)

c     Chlorine Routine
c     ================
      real rclo_cl(nlayermx)
      real pcl(nlayermx), pclo(nlayermx)
      real lcl(nlayermx), lclo(nlayermx)

c     Dust Parameters
c     ===============
      real dustradius
      parameter ( dustradius = 1.5e-6 ) ! meters
      real dustdens 
      parameter ( dustdens = 2.5e3) ! kg/m3
      real dust_clwt
      parameter ( dust_clwt = 0.5 ) ! Chlorine weight % in soil

c     Chlorine Production Parameters
c     ==============================
      REAL Chlorine_QSSA


      real clox_prod
      parameter( clox_prod = 1.E9)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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
c (1) : A model for studies of tropospheric ozone and nonmethane hydrocarbons: Model description and ozone results   
c (2) : Sander et al. (2014) [CAABA 4.0]



cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     volume mixing ratio -> density conversion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iesp = 1,nesp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Dust Number Density [BMT 02/2021]
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if ( iesp == i_dust ) then 
            do l = 1, lswitch-1
c     Mass Mixing Ratio -> Mass Concentration
            cc0(l,i_dust) = rm(l,i_dust)*(press(l)*100.)
     $                 /(Rnew(1,l)*t(l))
c     Mass Concentration -> Volume Concentration
            cc0(l,i_dust) = cc0(l,i_dust)/dustdens
c     Volume Concentration -> Number Density (cm-3)
            cc0(l,i_dust) = cc0(l,i_dust)/
     $                ( 1.333333*dustradius**3)*1.e-6
            cc(l,i_dust) = cc0(l,i_dust)
            enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Gas-phase number density
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         else 
             do l = 1,lswitch-1
                cc0(l,iesp) = rm(l,iesp)*dens(l)
                cc(l,iesp)  = cc0(l,iesp)
             end do
         endif 
      end do
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     co2 and nox number densities (cm-3)    
c      
c     nox mixing ratio: 6.e-10 (yung and demore, 1999)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
         nox(l) = 6.e-10*dens(l)
      end do 
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Peroxy-Radical Number Densities
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do l = 1, lswitch-1

        ro2(l) = 0.

      IF (igcm_ch3o2 .ne. 0) ro2(l) = ro2(l) + cc0(l,i_ch3o2)
      IF (igcm_hoch2o2 .ne. 0 ) ro2(l) = ro2(l) + cc0(l,i_hoch2o2)
      IF (igcm_c2h5o2 .ne. 0 ) ro2(l) = ro2(l) + cc0(l,i_c2h5o2)
      IF (igcm_hoch2ch2o2 .ne. 0 ) ro2(l) = ro2(l) + cc0(l,i_hoch2ch2o2)
      IF (igcm_ch3cooo .ne. 0 ) ro2(l) = ro2(l) + cc0(l,i_ch3cooo)
      IF (igcm_hcoch2o2 .ne. 0 ) ro2(l) = ro2(l) + cc0(l,i_hcoch2o2)

      IF ( ro2(l) .le. dens(l)*1.e-29 ) ro2(l) = 0.

      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over iterations in the chemical scheme
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      k_pseudo(:)=0. 

      do iter = 1,niter
       

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        nox species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        no2/no
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
c
         rno2_no = (d002(l)*cc(l,i_o3) + d003(l)*cc(l,i_ho2))
     $            /(j(l,j_no2) +
     $              d001(l)*max(cc(l,i_o),1.e-30*dens(l)))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        no
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         no(l) = nox(l)/(1. + rno2_no)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        no2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         no2(l) = no(l)*rno2_no



cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Tangent Linear
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C       IF ( tlm_on ) THEN 

C           call tlm_nox(iter, l, dens(l), nesp, 
C      $               nox(l), no(l), no2(l), rno2_no, 
C      $               cc(l,:), d001(l), d002(l), d003(l), 
C      $               j(l,j_no2), 
C      $              dccn_dpq)

C       ENDIF 



      end do


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hox species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        photochemical equilibrium for oh and ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h/ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c_

      do l = 1,lswitch-1
      
C     BMT v3.0
c            PRODUCTION(H)/CC[HO2]
C            CCCCCCCCCCCCCCCCCCCCC
      rh_ho2 = (
     $         ( b003(l)*cc(l,i_o1d)*cc(l,i_h2)
     $         + b008(l)*cc(l,i_o1d)*cc(l,i_ch4)
     $         + c002(l)*cc(l,i_o)*cc(l,i_oh)
     $         + c010(l)*cc(l,i_oh)*cc(l,i_h2)
     $         + e001(l)*cc(l,i_co)*cc(l,i_oh)
     $         + 0.49*cab002(l)*cc(l,i_ch4)*cc(l,i_o)
     $         + 0.956*cab004(l)*cc(l,i_o3)*cc(l,i_ch3)
     $         + 0.83*cab005(l)*cc(l,i_o)*cc(l,i_ch3)
     $         + j(l,j_h2o)*cc(l,i_h2o)
     $         + j(l,j_ch4_ch3_h)*cc(l,i_ch4)
     $         + 2.*j(l,j_ch4_3ch2_h_h)*cc(l,i_ch4)
     $         + j(l,j_ch4_ch_h2_h)*cc(l,i_ch4) 
     $         + j(l,j_ch2o_hco)*cc(l,i_hcho)
     $         + j(l,j_ch3oh)*cc(l,i_ch3oh)
     $         + cl008(l)*cc(l,i_cl)*cc(l,i_h2)
     $         + cl039(l)*cc(l,i_hcl)*cc(l,i_o1d)*0.22
     $         + j(l,j_hcl)*cc(l,i_hcl)
     $         )/MAX(cc(l,i_ho2),dens(l)*1.e-30)
c            LOSS[HO2]
C            CCCCCCCCC
     $        + c001(l)*cc(l,i_o)
     $        + c004(l)*cc(l,i_h)
     $        + c005(l)*cc(l,i_h)
     $        + c006(l)*cc(l,i_h)
     $        + c007(l)*cc(l,i_oh)
     $        + c008(l)*cc(l,i_ho2)
     $        + c015(l)*cc(l,i_o3)
     $        + c016(L)*cc(l,i_ho2)
     $        + cab006(l)*cc(l,i_ch3o2)
     $        + cab007(l)*cc(l,i_ch3o2)
     $        + cab019(l)*cc(l,i_hcho)
     $        + cab029(l)*cc(l,i_hoch2o2)*0.8
     $        + cl009(l)*cc(l,i_cl)
     $        + cl010(l)*cc(l,i_cl)
     $        + cl013(l)*cc(l,i_clo))
C            PRODUCTION(HO2)/CC[H]
C            CCCCCCCCCCCCCCCCCCCCC
     $      /( c011(l)*cc(l,i_o2) 
     $       + ( c009(l)*cc(l,i_h2o2)*cc(l,i_oh)
     $         + c012(l)*cc(l,i_o)*cc(l,i_h2o2)
     $         + c014(l)*cc(l,i_oh)*cc(l,i_o3)
     $         + cab011(l)*cc(l,i_ch3o2)*cc(l,i_oh)
     $         + cab013(l)*cc(l,i_ch3oh)*cc(l,i_oh)*0.85
     $         + cab015(l)*cc(l,i_ch3o)*cc(l,i_o2)
     $         + cab026(l)*cc(l,i_hco)*cc(l,i_o2)
     $         + cab028(l)*cc(l,i_hoch2o2)
     $         + cab030(l)*cc(l,i_hoch2o2)*ro2(l)
     $         + cab032(l)*cc(l,i_hcooh)*cc(l,i_oh)
     $         + cab035(l)*cc(l,i_hoch2oh)*cc(l,i_oh)
     $         + j(l,j_hoch2ooh)*cc(l,i_hoch2ooh)
     $         + cl011(l)*cc(l,i_cl)*cc(l,i_h2o2)
     $         + cl012(l)*cc(l,i_clo)*cc(l,i_oh)*0.94
     $         + cl050(l)*cc(l,i_clo3)*cc(l,i_oh)
     $         )/MAX(cc(l,i_h),dens(l)*1.e-30)
c            LOSS[H]
c            ccccccc
     $     + c003(l)*cc(l,i_o3)
     $     + c004(l)*cc(l,i_ho2)
     $     + c005(l)*cc(l,i_ho2)
     $     + c006(l)*cc(l,i_ho2)
     $     + c011(l)*cc(l,i_o2)
     $     + c018(l)*cc(l,i_h)
     $     + cab027(l)*cc(l,i_hco)
     $     + cl037(l)*cc(l,i_cl2)
     $     + cl041(l)*cc(l,i_hcl)
     $    )

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        oh/ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          
c     BMT v3.0
c     cccccccc
C                 PRODUCTION(OH)/CC[HO2]
C                 CCCCCCCCCCCCCCCCCCCCCC
      roh_ho2 = ( c001(l)*cc(l,i_o)
     $        +   c004(l)*cc(l,i_h)*2.
     $        +   c015(l)*cc(l,i_o3)
     $        +   0.2*cab029(l)*cc(l,i_hoch2o2)
     $        +   2.*j(l,j_ho2)
     $        +   cl010(l)*cc(l,i_cl)
     $        +   ( b002(l)*cc(l,i_o1d)*cc(l,i_h2o)*2.
     $            + b003(l)*cc(l,i_o1d)*cc(l,i_h2)
     $            + b007(l)*cc(l,i_o1d)*cc(l,i_ch4)
     $            + c003(l)*cc(l,i_o3)*cc(l,i_h)
     $            + c012(l)*cc(l,i_h2o2)*cc(l,i_o)
     $            + cab020(l)*cc(l,i_hcho)*cc(l,i_o)
     $            + cab021(l)*cc(l,i_hco)*cc(l,i_o)
     $            + j(l,j_h2o)*cc(l,i_h2o)
     $            + j(l,j_ch3o2h)*cc(l,i_ch3ooh)
     $            + 2.*j(l,j_h2o2)*cc(l,i_h2o2)
     $            + j(l,j_hoch2ooh)*cc(l,i_hoch2ooh)
     $            + cl018(l)*cc(l,i_cl)*cc(l,i_ch3ooh)
     $            + cl039(l)*cc(l,i_hcl)*cc(l,i_o1d)*0.66
     $            + cl040(l)*cc(l,i_hcl)*cc(l,i_o)
     $            + cl041(l)*cc(l,i_hocl)*cc(l,i_o)
     $            + cl054(l)*cc(l,i_hocl)*cc(l,i_clo4)
     $            + cl055(l)*cc(l,i_hocl)*cc(l,i_cl)
     $            + j(l,j_hocl)*cc(l,i_hocl)
     $            )/MAX(cc(l,i_ho2),dens(l)*1.e-30)
C                 LOSS(HO2)
C                 CCCCCCCCC
     $        + c001(l)*cc(l,i_o)
     $        + c004(l)*cc(l,i_h)
     $        + c005(l)*cc(l,i_h)
     $        + c006(l)*cc(l,i_h)
     $        + c007(l)*cc(l,i_oh)
     $        + c008(l)*cc(l,i_ho2)
     $        + c015(l)*cc(l,i_o3)
     $        + c016(L)*cc(l,i_ho2)
     $        + cab006(l)*cc(l,i_ch3o2)
     $        + cab007(l)*cc(l,i_ch3o2)
     $        + cab019(l)*cc(l,i_hcho)
     $        + cab029(l)*cc(l,i_hoch2o2)*0.8
     $        + cl009(l)*cc(l,i_cl)
     $        + cl010(l)*cc(l,i_cl)
     $        + cl013(l)*cc(l,i_clo))  
     $        /
c             PRODUCTION(HO2)/CC[OH]
c             CCCCCCCCCCCCCCCCCCCCCC
     $        ( c009(l)*cc(l,i_h2o2)
     $        + c014(l)*cc(l,i_o3)
     $        + cab011(l)*cc(l,i_ch3o2)
     $        + cab013(l)*cc(l,i_ch3oh)*0.85
     $        + cab032(l)*cc(l,i_hcooh)
     $        + cab035(L)*cc(l,i_hoch2oh)
     $        + cl012(l)*cc(l,i_clo)*0.94
     $        + cl050(l)*cc(l,i_clo3)
     $        + ( c011(l)*cc(l,i_h)*cc(l,i_o2)
     $         + c012(l)*cc(l,i_o)*cc(l,i_h2o2)
     $         + cab015(l)*cc(l,i_ch3o)*cc(l,i_o2)
     $         + cab026(l)*cc(l,i_hco)*cc(l,i_o2)
     $         + cab028(l)*cc(l,i_hoch2o2)
     $         + cab030(l)*cc(l,i_hoch2o2)*ro2(l)
     $         + cl011(l)*cc(L,i_cl)*cc(l,i_h2o2)
     $         + j(l,j_hoch2ooh)*cc(l,i_hoch2ooh)
     $          )/MAX(cc(l,i_oh),dens(l)*1.e-30)
c             LOSS(OH)
C             CCCCCCCC
     $       + c002(l)*cc(l,i_o)
     $       + c007(l)*cc(l,i_ho2)
     $       + c009(l)*cc(l,i_h2o2)
     $       + c010(l)*cc(l,i_h2)
     $       + c013(l)*cc(l,i_oh)
     $       + c014(l)*cc(l,i_o3)
     $       + c017(l)*cc(l,i_oh)
     $       + e001(l)*cc(l,i_co)
     $       + cab001(l)*cc(l,i_ch4)
     $       + cab011(l)*cc(l,i_ch3o2)
     $       + cab013(l)*cc(l,i_ch3oh)
     $       + cab014(l)*cc(l,i_ch3ooh)
     $       + cab018(l)*cc(l,i_hcho)
     $       + cab025(l)*cc(l,i_hco)
     $       + cab032(l)*cc(l,i_hcooh)
     $       + cab033(l)*cc(l,i_hoch2ooh)
     $       + cab035(l)*cc(l,i_hoch2oh)
     $       + cl012(l)*cc(l,i_clo)
     $       + cl014(l)*cc(l,i_clo)
     $       + cl027(l)*cc(l,i_hclo4)
     $       + cl036(l)*Cc(l,i_oh)
     $       + cl048(l)*cc(l,i_clo3)
     $       + cl049(l)*cc(l,i_clo3)
     $       + cl050(l)*cc(l,i_clo3)
     $       + cl058(l)*cc(l,i_oclo)
     $       )

C           IF ( tlm_on ) THEN 

C          cc_prev = cc(l,:)

     
C          ENDIF
         
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         cc(l,i_h) = cc(l,i_hox)
     $                 /(1. + (1. + roh_ho2)/rh_ho2 ) 
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         cc(l,i_ho2) = cc(l,i_h)/rh_ho2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        oh
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         cc(l,i_oh) = cc(l,i_ho2)*roh_ho2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c      Tangent Linear Model - Odd-Hydrogen
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

C       IF ( tlm_on ) THEN
      
C       call tlm_hox(iter, l, dens(l), sza,
C      $  rh_ho2, roh_ho2,
C      $  dt, ptimestep,
C      $  nesp, cc(l,:), cc_prev(:), ro2(l),
C      $  j(l,:),
C      $  a001(l), a002(l), a003(l),
C      $  b001(l), b002(l), b003(l), b004(l), b005(l), b006(l),
C      $  b007(l), b008(l), b009(l),
C      $  c001(l), c002(l), c003(l), c004(l), c005(l), c006(l), 
C      $  c007(l), c008(l), c009(l), c010(l), c011(l), c012(l),
C      $  c013(l), c014(l), c015(l), c016(l), c017(l), c018(l),
C      $  e001(l), e002(l), e003(l),
C      $  cab001(l), cab002(l), cab003(l),
C      $  cab004(l), cab005(l), cab006(l),
C      $  cab007(l), cab008(l), cab009(l),
C      $  cab010(l), cab011(l), cab012(l),
C      $  cab013(l), cab014(l), cab015(l),
C      $  cab016(l), cab017(l), cab018(l), 
C      $  cab019(l), cab020(l), cab021(l), 
C      $  cab022(l), cab023(l), cab024(l), 
C      $  cab025(l), cab026(l), cab027(l), 
C      $  cab028(l), cab029(l), cab030(l),
C      $  cab031(l), cab032(l), cab033(l),
C      $  cab034(l), cab035(l), cab036(l),
C      $  cab037(l), cab038(l), cab039(l),
C      $  cab040(l), cab041(l), cab042(l),
C      $  cab043(l), cab044(l), cab045(l),
C      $  cab046(l), cab047(l), cab048(l),
C      $  cab049(l), cab050(l), cab051(l),
C      $  cab052(l), cab053(l), cab054(l),
C      $  cab055(l), cab056(l), cab057(l),
C      $  cab058(l), cab059(l), cab060(l),
C      $  cab061(l), cab062(l), cab063(l),
C      $  cab064(l), cab065(l), cab066(l),
C      $  cab067(l), cab068(l), cab069(l),
C      $  cab070(l), cab071(l), cab072(l),
C      $  cab073(l), cab074(l), cab075(l),
C      $  cab076(l), cab077(l), cab078(l),
C      $  cab079(l), cab080(l), cab081(l),
C      $  cab082(l), cab083(l), cab084(l),
C      $  cab085(l), cab086(l), cab087(l),
C      $  cab088(l), cab089(l), cab090(l),
C      $  cab091(l), cab092(l), cab093(l),
C      $  cab094(l), cab095(l), cab096(l),
C      $  cab097(l), cab098(l), cab099(l),
C      $  cab100(l), cab101(l), cab102(l),
C      $  cab103(l), cab104(l), cab105(l),
C      $  cab106(l), cab107(l),
C      $  cl001(l), cl002(l), cl003(l), 
C      $  cl004(l), cl005(l), cl006(l), 
C      $  cl007(l), cl008(l), cl009(l), 
C      $  cl010(l), cl011(l), cl012(l),
C      $  cl013(l), cl014(l), cl015(l),
C      $  cl016(l), cl017(l), cl018(l), 
C      $  cl019(l), cl020(l), cl021(l), 
C      $  cl022(l), cl023(l), cl024(l), 
C      $  cl025(l), cl026(l), cl027(l), 
C      $  cl028(l), cl029(l), cl030(l), 
C      $  cl031(l), cl032(l), cl033(l),
C      $  cl034(l), cl035(l), cl036(l),
C      $  cl037(l), cl038(l), cl039(l),
C      $  cl040(l), cl041(l), cl042(l),
C      $  cl043(l), cl044(l), cl045(l),
C      $  cl046(l), cl047(l), cl048(l),
C      $  cl049(l), cl050(l), cl051(l),
C      $  cl052(l), cl053(l), cl054(l), 
C      $  cl055(l), cl056(l), cl057(l),
C      $  cl058(l), cl059(l),
C      $  dccn_dpq, dcc0_dpq,
C      $  dHOX_dPQ, dHOX0_dPQ,
C      $  k_pseudo(l))
C        ENDIF

      end do
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ox species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        day:
c        - o1d at photochemical equilibrium 
c        - o3 at photochemical equilibrium
c        - continuity equation for ox
c        night:
c        - o1d = 0
c        - continuity equation for o3
c        - continuity equation for o
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
      if (sza .le. 95.) then
c
         
         do l = 1,lswitch-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o(1d)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            cc(l,i_o1d) = (j(l,j_co2_o1d)*cc(l,i_co2)
     $                   + j(l,j_o2_o1d)*cc(l,i_o2)
     $                   + j(l,j_o3_o1d)*cc(l,i_o3))
     $                   /(b001(l)*cc(l,i_co2)
     $                   + b002(l)*cc(l,i_h2o)
     $                   + b003(l)*cc(l,i_h2)
     $                   + b004(l)*cc(l,i_o2)
     $                   + b005(l)*cc(l,i_o3)
     $                   + b006(l)*cc(l,i_o3)
     $                   + b007(l)*cc(l,i_ch4)
     $                   + b008(l)*cc(l,i_ch4)
     $                   + b009(l)*cc(l,i_ch4)
     $                   + cl034(l)*Cc(l,i_cl2)
     $                   + cl035(l)*cc(l,i_cl2)
     $                   + cl039(l)*cc(l,i_hcl))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o/o3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccc Original Routine
            ro_o3 = (j(l,j_o3_o1d) + j(l,j_o3_o)
     $              + a003(l)*cc(l,i_o)
     $              + c003(l)*cc(l,i_h)
     $              + c014(l)*cc(l,i_oh)
     $              + c015(l)*cc(l,i_ho2)
     $              + b005(l)*cc(l,i_o1d) ! BMT Additions below
     $              + b006(l)*cc(l,i_o1d)
     $              + d002(l)*no(l)
     $              + cab004(l)*cc(l,i_ch3)
     $              + cab010(l)*cc(l,i_ch3o2)
     $              + cab016(l)*cc(l,i_ch3o)
     $              + cl052(l)*cc(l,i_oclo)
     $              )
     $              /(a001(l)*cc(l,i_o2)
     $              + cl051(l)*cc(l,i_oclo))

     
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           cc(l,i_o3) = cc(l,i_ox)/(1. + ro_o3)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            cc(l,i_o) = cc(l,i_o3)*ro_o3
     
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ox = o + o3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            production(l,i_ox) =
     $                   + j(l,j_co2_o)*cc(l,i_co2)
     $                   + j(l,j_co2_o1d)*cc(l,i_co2)
     $                   + j(l,j_ho2)*cc(l,i_ho2)
     $                   + 2.*j(l,j_o2_o)*cc(l,i_o2)
     $                   + 2.*j(l,j_o2_o1d)*cc(l,i_o2)
     $                   + c006(l)*cc(l,i_h)*cc(l,i_ho2)
     $                   + c013(l)*cc(l,i_oh)*cc(l,i_oh)
     $                   + d003(l)*cc(l,i_ho2)*no(l)
     $                   + cab074(l)*cc(l,i_ch3cooo)*cc(l,i_ho2)
     $                   + cab097(l)*cc(l,i_hoch2co3)*cc(l,i_ho2)
     $                   + cl034(l)*cc(l,i_cl2)*cc(l,i_o1d)
     $                   + cl039(l)*cc(l,i_hcl)*cc(l,i_o1d)*0.12
     $                   + j(l,j_cloo)*cc(l,i_cloo)
     $                   + j(l,j_oclo)*cc(l,i_oclo)
     $                   + j(l,j_clo)*cc(l,i_clo)


            loss(l,i_ox) = 2.*a002(l)*cc(l,i_o)*cc(l,i_o)
     $                   + 2.*a003(l)*cc(l,i_o)*cc(l,i_o3)
     $                   + c001(l)*cc(l,i_ho2)*cc(l,i_o)
     $                   + c002(l)*cc(l,i_oh)*cc(l,i_o)
     $                   + c003(l)*cc(l,i_h)*cc(l,i_o3)
     $                   + c012(l)*cc(l,i_o)*cc(l,i_h2o2)
     $                   + c014(l)*cc(l,i_o3)*cc(l,i_oh)
     $                   + c015(l)*cc(l,i_o3)*cc(l,i_ho2)
     $                   + d001(l)*cc(l,i_o)*no2(l)
     $                   + e002(l)*cc(l,i_o)*cc(l,i_co)
     $                   + cab002(l)*cc(l,i_o)*cc(l,i_ch4)
     $                   + cab004(l)*cc(l,i_o3)*cc(l,i_ch3)
     $                   + cab010(l)*cc(l,i_o3)*cc(l,i_ch3o2)
     $                   + cab012(l)*cc(l,i_o)*cc(l,i_ch3o2)
     $                   + cab016(l)*cc(l,i_o3)*cc(l,i_ch3o)
     $                   + cab017(l)*cc(l,i_o)*cc(l,i_ch3o)
     $                   + cab020(l)*cc(l,i_o)*cc(l,i_hcho)
     $                   + cab021(l)*cc(l,i_hco)*cc(l,i_o)
     $                   + cab037(l)*cc(l,i_c2h6)*cc(l,i_o)
     $                   + cl001(l)*cc(l,i_cl)*cc(l,i_o3)
     $                   + cl002(l)*cc(l,i_clo)*cc(l,i_o)
     $                   + cl040(l)*cc(l,i_hcl)*cc(l,i_o)
     $                   + cl041(l)*cc(l,i_hocl)*cc(l,i_o)
     $                   + cl044(l)*cc(l,i_cl)*cc(l,i_o3)
     $                   + cl051(l)*cc(l,i_oclo)*cc(l,i_o)
     $                   + cl052(l)*cc(l,i_oclo)*cc(l,i_o3)
     $                   + cl057(l)*cc(l,i_oclo)*cc(l,i_o)

            loss(l,i_ox) = loss(l,i_ox)/cc(l,i_ox)

            
            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c      Tangent Linear Model - Odd-Oxygen (daytime)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
            
C       IF ( tlm_on ) THEN
     
      
C        cc_ox_next = (cc0(l,i_ox) + production(l,i_ox)*dt)
C      $            /(1. + loss(l,i_ox)*dt)  
       
C       call tlm_ox(iter, l, dens(l),
C      $  ro_o3,
C      $  dt, ptimestep,
C      $  nesp, cc(l,:), cc_prev(:),
C      $  j(l,:), loss(l,i_ox), production(l,i_ox), cc_ox_next,
C      $  a001(l), a002(l), a003(l),
C      $  b001(l), b002(l), b003(l), b004(l), b005(l), b006(l),
C      $  b007(l), b008(l), b009(l),
C      $  c001(l), c002(l), c003(l), c004(l), c005(l), c006(l), 
C      $  c007(l), c008(l), c009(l), c010(l), c011(l), c012(l),
C      $  c013(l), c014(l), c015(l), c016(l), c017(l), c018(l),
C      $  d001(l), d002(l), d003(l),
C      $  e001(l), e002(l), e003(l),
C      $  cab001(l), cab002(l), cab003(l),
C      $  cab004(l), cab005(l), cab006(l),
C      $  cab007(l), cab008(l), cab009(l),
C      $  cab010(l), cab011(l), cab012(l),
C      $  cab013(l), cab014(l), cab015(l),
C      $  cab016(l), cab017(l), cab018(l), 
C      $  cab019(l), cab020(l), cab021(l), 
C      $  cab022(l), cab023(l), cab024(l), 
C      $  cab025(l), cab026(l), cab027(l), 
C      $  cab028(l), cab029(l), cab030(l),
C      $  cab031(l), cab032(l), cab033(l),
C      $  cab034(l), cab035(l), cab036(l),
C      $  cab037(l), cab038(l), cab039(l),
C      $  cab040(l), cab041(l), cab042(l),
C      $  cab043(l), cab044(l), cab045(l),
C      $  cab046(l), cab047(l), cab048(l),
C      $  cab049(l), cab050(l), cab051(l),
C      $  cab052(l), cab053(l), cab054(l),
C      $  cab055(l), cab056(l), cab057(l),
C      $  cab058(l), cab059(l), cab060(l),
C      $  cab061(l), cab062(l), cab063(l),
C      $  cab064(l), cab065(l), cab066(l),
C      $  cab067(l), cab068(l), cab069(l),
C      $  cab070(l), cab071(l), cab072(l),
C      $  cab073(l), cab074(l), cab075(l),
C      $  cab076(l), cab077(l), cab078(l),
C      $  cab079(l), cab080(l), cab081(l),
C      $  cab082(l), cab083(l), cab084(l),
C      $  cab085(l), cab086(l), cab087(l),
C      $  cab088(l), cab089(l), cab090(l),
C      $  cab091(l), cab092(l), cab093(l),
C      $  cab094(l), cab095(l), cab096(l),
C      $  cab097(l), cab098(l), cab099(l),
C      $  cab100(l), cab101(l), cab102(l),
C      $  cab103(l), cab104(l), cab105(l),
C      $  cab106(l), cab107(L),
C      $  cl001(l), cl002(l), cl003(l), 
C      $  cl004(l), cl005(l), cl006(l), 
C      $  cl007(l), cl008(l), cl009(l), 
C      $  cl010(l), cl011(l), cl012(l),
C      $  cl013(l), cl014(l), cl015(l),
C      $  cl016(l), cl017(l), cl018(l), 
C      $  cl019(l), cl020(l), cl021(l), 
C      $  cl022(l), cl023(l), cl024(l), 
C      $  cl025(l), cl026(l), cl027(l), 
C      $  cl028(l), cl029(l), cl030(l), 
C      $  cl031(l), cl032(l), cl033(l),
C      $  cl034(l), cl035(l), cl036(l),
C      $  cl037(l), cl038(l), cl039(l),
C      $  cl040(l), cl041(l), cl042(l),
C      $  cl043(l), cl044(l), cl045(l),
C      $  cl046(l), cl047(l), cl048(l),
C      $  cl049(l), cl050(l), cl051(l),
C      $  cl052(l), cl053(l), cl054(l), 
C      $  cl055(l), cl056(l), cl057(l),
C      $  cl058(l), cl059(l),
C      $  no(l), no2(l),
C      $  dccn_dpq, dcc0_dpq,
C      $  dOX_dPQ, dOX0_dPQ)



C         ENDIF
       
       
         end do
     
       
      else
c
         do l = 1,lswitch-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o(1d)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            cc(l,i_o1d) = 0.


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            production(l,i_o3) = a001(l)*cc(l,i_o2)*cc(l,i_o)
     $                        + cab074(l)*cc(l,i_ch3cooo)*cc(l,i_ho2)
     $                        + cab097(l)*cc(l,i_hoch2co3)*cc(l,i_ho2)
C c

            loss(l,i_o3) = a003(l)*cc(l,i_o)
     $                   + c003(l)*cc(l,i_h)
     $                   + c014(l)*cc(l,i_oh)
     $                   + c015(l)*cc(l,i_ho2)
     $                   + cab004(l)*cc(l,i_ch3)
     $                   + cab010(l)*cc(l,i_ch3o2)
     $                   + cab016(l)*cc(l,i_ch3o)
     $                   + cl001(l)*cc(L,i_cl)
     $                   + cl044(l)*cc(l,i_cl)
     $                   + cl051(l)*cc(l,i_oclo)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            production(l,i_o) = c006(l)*cc(l,i_h)*cc(l,i_ho2)
     $                        + c013(l)*cc(l,i_oh)*cc(l,i_oh)
c

            loss(l,i_o)  =  a001(l)*cc(l,i_o2)
     $                   + 2.*a002(l)*cc(l,i_o)
     $                   + a003(l)*cc(l,i_o3)
     $                   + c001(l)*cc(l,i_ho2)
     $                   + c002(l)*cc(l,i_oh)
     $                   + c012(l)*cc(l,i_h2o2)
     $                   + e002(l)*cc(l,i_co)
     $                   + cab002(l)*cc(l,i_ch4)
     $                   + cab005(l)*cc(l,i_ch3)
     $                   + cab012(l)*cc(l,i_ch3o2)
     $                   + cab017(l)*Cc(L,i_ch3o)
     $                   + cab020(l)*cc(l,i_hcho)
     $                   + cab021(l)*cc(l,i_hco)
     $                   + cab037(l)*cc(l,i_c2h6)
     $                   + cl002(l)*cc(l,i_clo)
     $                   + cl040(l)*cc(l,i_hcl)
     $                   + cl042(l)*cc(l,i_hocl)
     $                   + cl051(l)*cc(l,i_oclo)
     $                   + cl057(l)*cc(l,i_oclo)

       end do 

      end if


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Odd-Chlorine Chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( igcm_cl .ne. 0 ) THEN 

      DO l = 1, lswitch-1

      ! Production of Cl via ClO 
      pcl(l) = j(l,j_clo) 
     $    + cl002(l)*cc(l,i_o)
     $    + cl004(l)*cc(l,i_clo)*2.
     $    + cl005(l)*cc(l,i_clo)
     $    + cl012(l)*cc(l,i_oh)*0.94
      ! Loss of Cl 
      lcl(l) = cl001(l)*cc(l,i_o3)
     $    + cl010(l)*cc(l,i_ho2)
     $    + cl022(l)*cc(l,i_ch3o2)*0.5
     $    + cl024(l)*cc(l,i_cloo) 
     $    + cl053(l)*cc(l,i_clo4)
     $    + cl056(l)*cc(l,i_hocl)
     $    + cl059(l)*cc(l,i_oclo)

      ! Production of ClO via Cl
      pclo(l) = cl001(l)*cc(l,i_o3)
     $     + cl010(l)*cc(l,i_ho2)
     $     + cl022(l)*cc(l,i_ch3o2)*0.5
     $     + cl024(l)*cc(l,i_cloo)*2.
     $     + cl053(l)*cc(l,i_clo4)
     $     + cl056(l)*cc(l,i_hocl)
     $     + cl059(l)*cc(l,i_oclo)*2.

      ! Loss of ClO
      lclo(l) = j(l,j_clo)
     $     + cl002(l)*cc(l,i_o)
     $     + cl004(l)*cc(l,i_clo)
     $     + cl005(l)*cc(l,i_clo)
     $     + cl012(l)*cc(l,i_oh)

      rclo_cl(l) = (pclo(l) + lcl(l))/(pcl(l) + lclo(l))

c     -----------------------------------------------
C     Atmospheric state prior to timestep advancement
c     neeeded for TLM
c     ----------------------------------------------- 
C       IF ( tlm_on ) THEN 
C             cc_prev = cc(l,:)
C       ENDIF ! tlm_on 

        ! Cl
        cc(l,i_cl) = cc(l,i_clox)/(1. + rclo_cl(l))

        ! ClO
        cc(l,i_clo) = cc(l,i_cl)*rclo_cl(l)
  
c     ==============
C     Tangent Linear 
c     ==============
C       IF ( tlm_on ) THEN 

C       call tlm_clox(l, rclo_cl(l), iter,  dt, dens(l),
C      $             pcl(l), lcl(l), pclo(l), lclo(l),
C      $             cc(l,:), cc_prev(:),
C      $             nesp,
C      $             dccn_dpq, dcc0_dpq,
C      $             j(l,:),
C      $             a001(l), a002(l), a003(l),
C      $             b001(l), b002(l), b003(l), b004(l), b005(l), b006(l),
C      $             b007(l), b008(l), b009(l),
C      $             c001(l), c002(l), c003(l), c004(l), c005(l), c006(l), 
C      $             c007(l), c008(l), c009(l), c010(l), c011(l), c012(l),
C      $             c013(l), c014(l), c015(l), c016(l), c017(l), c018(l),
C      $             e001(l), e002(l), e003(l),
C      $             cab001(l), cab002(l), cab003(l),
C      $             cab004(l), cab005(l), cab006(l),
C      $             cab007(l), cab008(l), cab009(l),
C      $             cab010(l), cab011(l), cab012(l),
C      $             cab013(l), cab014(l), cab015(l),
C      $             cab016(l), cab017(l), cab018(l), 
C      $             cab019(l), cab020(l), cab021(l), 
C      $             cab022(l), cab023(l), cab024(l), 
C      $             cab025(l), cab026(l), cab027(l), 
C      $             cab028(l), cab029(l), cab030(l),
C      $             cab031(l), cab032(l), cab033(l),
C      $             cab034(l), cab035(l), cab036(l),
C      $             cab037(l), cab038(l), cab039(l),
C      $             cab040(l), cab041(l), cab042(l),
C      $             cab043(l), cab044(l), cab045(l),
C      $             cab046(l), cab047(l), cab048(l),
C      $             cab049(l), cab050(l), cab051(l),
C      $             cab052(l), cab053(l), cab054(l),
C      $             cab055(l), cab056(l), cab057(l),
C      $             cab058(l), cab059(l), cab060(l),
C      $             cab061(l), cab062(l), cab063(l),
C      $             cab064(l), cab065(l), cab066(l),
C      $             cab067(l), cab068(l), cab069(l),
C      $             cab070(l), cab071(l), cab072(l),
C      $             cab073(l), cab074(l), cab075(l),
C      $             cab076(l), cab077(l), cab078(l),
C      $             cab079(l), cab080(l), cab081(l),
C      $             cab082(l), cab083(l), cab084(l),
C      $             cab085(l), cab086(l), cab087(l),
C      $             cab088(l), cab089(l), cab090(l),
C      $             cab091(l), cab092(l), cab093(l),
C      $             cab094(l), cab095(l), cab096(l),
C      $             cab097(l), cab098(l), cab099(l),
C      $             cab100(l), cab101(l), cab102(l),
C      $             cab103(l), cab104(l), cab105(l),
C      $             cab106(l), cab107(l),
C      $             cl001(l), cl002(l), cl003(l) 
C      $            ,cl004(l), cl005(l), cl006(l) 
C      $            ,cl007(l), cl008(l), cl009(l) 
C      $            ,cl010(l), cl011(l), cl012(l) 
C      $            ,cl013(l), cl014(l), cl015(l) 
C      $            ,cl016(l), cl017(l), cl018(l) 
C      $            ,cl019(l), cl020(l), cl021(l) 
C      $            ,cl022(l), cl023(l), cl024(l) 
C      $            ,cl025(l), cl026(l), cl027(l) 
C      $            ,cl028(l), cl029(l), cl030(l) 
C      $            ,cl031(l), cl032(l), cl033(l)
C      $            ,cl034(l), cl035(l), cl036(l)
C      $            ,cl037(l), cl038(l), cl039(l)
C      $            ,cl040(l), cl041(l), cl042(l)
C      $            ,cl043(l), cl044(l), cl045(l)
C      $            ,cl046(l), cl047(l), cl048(l)
C      $            ,cl049(l), cl050(l), cl051(l)
C      $            ,cl052(l), cl053(l), cl054(l)
C      $            ,cl055(l), cl056(l), cl057(l)
C      $            ,cl058(l), cl059(l) )

C       ENDIF ! tlm_on 

      ENDDO ! l = 1, lswitch-1

      ENDIF ! igcm_cl .ne. 0 



cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     other species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        co2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         production(l,i_co2) = e001(l)*cc(l,i_oh)*cc(l,i_co)
     $                       + e002(l)*cc(l,i_o)*cc(l,i_co)
     $                       + t002(l)*cc(l,i_ch4)*16./44. ! ajout 20090401
     $                       + t003(l)*cc(l,i_co2)*8./44.  ! ajout 20090401
     $                 + cab032(l)*cc(l,i_oh)*cc(l,i_hcooh)
     $                 + cab067(l)*cc(l,i_ch3cooh)*cc(l,i_oh)
     $                 + cab072(l)*cc(l,i_ch3cooo)*cc(l,i_ho2)
     $                 + cab075(l)*cc(l,i_ch3cooo)*ro2(l) 
     $                 + cab084(l)*cc(l,i_o2)*cc(l,i_hcoco)
     $                 + j(l,j_ch3o2h)*cc(l,i_hoch2co3h)
     $                 + j(l,j_ch2o_co)*3.95*cc(l,i_hcoco2h)
     $                 + cab092(l)*cc(l,i_hoch2co)*cc(l,i_o2)
     $                 + cab093(l)*cc(l,i_hoch2co3)*ro2(l) 
     $                 + cab095(l)*cc(l,i_hoch2co3)*cc(l,i_ho2)
     $                 + 0.91*cab098(l)*cc(l,i_hoch2co2h)*cc(l,i_oh)
     $                 + cab099(l)*cc(l,i_oh)*cc(l,i_hcoco2h)
     $                 + cab103(l)*cc(l,i_oh)*cc(l,i_hcoco3h)
     $                 + cab104(l)*cc(l,i_hcoco3)*ro2(l) 
     $                 + cab106(l)*cc(l,i_hcoco3)*cc(l,i_ho2)
     $                 + j(l,j_ch3coooh)*cc(l,i_ch3coooh)
     $                 + 2.*( j(l,j_ch3o2h)
     $                 + j(l,j_hoch2cho_co)
     $                 + j(l,j_hoch2cho_hco)
     $                 + j(l,j_hoch2cho_oh) )*cc(l,i_hcoco3h)
     
         loss(l,i_co2) = j(l,j_co2_o)
     $                 + j(l,j_co2_o1d)
     $                 + t003(l)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        co
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_co) = j(l,j_co2_o)*cc(l,i_co2)
     $                      + j(l,j_co2_o1d)*cc(l,i_co2)
     $                      + t003(l)*cc(l,i_co2)
     $                   + 0.17*cab005(l)*cc(l,i_ch3)*cc(l,i_o)     
     $                   + j(l,j_ch2o_co)*cc(l,i_hcho)
     $                   + cab021(l)*cc(l,i_hco)*cc(l,i_o)
     $                   + cab022(l)*cc(l,i_hco)*cc(l,i_ch3)
     $                   + cab024(l)*cc(l,i_hco)*cc(l,i_hco)
     $                   + cab025(l)*cc(l,i_hco)*cc(l,i_oh)
     $                   + cab026(l)*cc(l,i_hco)*cc(l,i_o2)
     $                   + cab027(l)*cc(l,i_hco)*cc(l,i_h)
     $                   + cab071(l)*cc(l,i_ch3co)*cc(l,i_o2)
     $                   + 0.6*cab078(l)*cc(l,i_hcoch2o2)*ro2(l) 
     $                   + cab080(l)*cc(l,i_hcoch2o2)*cc(l,i_ho2)
     $                   + cab082(l)*cc(l,i_hcoco)
     $                   + cab084(l)*cc(l,i_o2)*cc(l,i_hcoco)
     $                   + cab086(l)*cc(l,i_oh)*cc(l,i_hooch2cho)
     $                   + j(l,j_ch2o_co)*3.95*cc(l,i_hcoco2h)
     $                   + cab099(l)*cc(l,i_oh)*cc(l,i_hcoco2h)
     $                   + cab103(l)*cc(l,i_oh)*cc(l,i_hcoco3h)
     $                   + cab104(l)*cc(l,i_hcoco3)*ro2(l) 
     $                   + cab106(l)*cc(l,i_hcoco3)*cc(l,i_ho2)
     $                   + j(l,j_ch3cho_ch4)*cc(l,i_ch3cho)
     $                   + 2.*j(l,j_glyox_h2)*cc(l,i_glyox)
     $                   + j(l,j_glyox_hcho)*cc(l,i_glyox)
     $                 + 2.*( j(l,j_ch3o2h)
     $                 + j(l,j_hoch2cho_co)
     $                 + j(l,j_hoch2cho_hco)
     $                 + j(l,j_hoch2cho_oh) )*cc(l,i_hcoco3h)
     $                 + ( j(l,j_ch3o2h)
     $                 + j(l,j_hoch2cho_co)
     $                 + j(l,j_hoch2cho_hco)
     $                 + j(l,j_hoch2cho_oh) )*cc(l,i_hoch2cho)
     $                 + cl033(l)*cc(l,i_clco)

         loss(l,i_co) = e001(l)*cc(l,i_oh)
     $                + e002(l)*cc(l,i_o)
     $                + cl023(l)*cc(l,i_cl)
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_o2) = 
     $                  j(l,j_o3_o)*cc(l,i_o3)
     $                + j(l,j_o3_o1d)*cc(l,i_o3)
     $                + a002(l)*cc(l,i_o)*cc(l,i_o)
     $                + 2.*a003(l)*cc(l,i_o)*cc(l,i_o3)
     $                + 2.*b005(l)*cc(l,i_o1d)*cc(l,i_o3)
     $                + b006(l)*cc(l,i_o1d)*cc(l,i_o3)
     $                + c001(l)*cc(l,i_o)*cc(l,i_ho2)
     $                + c002(l)*cc(l,i_o)*cc(l,i_oh)
     $                + c003(l)*cc(l,i_h)*cc(l,i_o3)
     $                + c005(l)*cc(l,i_h)*cc(l,i_ho2)
     $                + c007(l)*cc(l,i_oh)*cc(l,i_ho2)
     $                + c008(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                + c014(l)*cc(l,i_o3)*cc(l,i_oh)
     $                + 2.*c015(l)*cc(l,i_o3)*cc(l,i_ho2)
     $                + c016(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                + d001(l)*cc(l,i_o)*no2(l)
     $                + cab004(l)*cc(l,i_ch3)*cc(l,i_o3)
     $                + cab006(l)*cc(l,i_ch3o2)*cc(l,i_ho2)
     $                + cab007(L)*cc(l,i_ch3o2)*cc(l,i_ho2)
     $                + 0.5*cab008(l)*cc(l,i_ch3o2)*ro2(l) 
     $                + 0.5*cab009(l)*cc(l,i_ch3o2)*ro2(l) 
     $                + 2.*cab010(l)*cc(l,i_ch3o2)*cc(l,i_o3)
     $                + cab012(l)*cc(l,i_ch3o2)*cc(l,i_o)
     $                + cab016(l)*cc(l,i_ch3o)*cc(l,i_o3)
     $                + 0.75*cab017(l)*cc(l,i_ch3o)*cc(l,i_o)
     $                + 0.8*cab029(l)*cc(l,i_hoch2o2)*cc(l,i_ho2)
     $                + 0.5*cab031(l)*cc(l,i_hoch2o2)*ro2(l) 
     $                + cab043(l)*cc(l,i_c2h5o2)*cc(l,i_ho2)
     $                + cab044(l)*cc(l,i_c2h5o2)*ro2(l) 
     $                + cab065(l)*cc(l,i_ch3choho2)*cc(L,i_ho2)
     $                + cab105(l)*cc(l,i_hcoco3)*ro2(l)*0.1
     $                + cl001(l)*cc(l,i_cl)*cc(l,i_o3)
     $                + cl002(l)*cc(l,i_clo)*cc(l,i_o)
     $                + cl003(l)*cc(l,i_clO)*cc(l,i_clo)
     $                + cl004(l)*cc(l,i_clo)*cc(l,i_clo)
     $                + cl009(l)*cc(l,i_cl)*cc(l,i_ho2)
     $                + cl012(l)*cc(l,i_clo)*cc(l,i_oh)*0.06
     $                + cl013(l)*cc(l,i_clo)*cc(l,i_ho2)
     $                + cl020(l)*cc(l,i_clo)*cc(l,i_ch3o2)
     $                + cl025(l)*cc(l,i_cloo)*cc(l,i_cl)
     $                + cl026(l)*cc(l,i_cloo) 
     $                + cl052(l)*cc(l,i_oclo)*Cc(l,i_o3)
     $                + cl057(l)*cc(l,i_oclo)*cc(l,i_o)
     $                + cl058(l)*cc(l,i_oclo)*Cc(l,i_oh)

         loss(l,i_o2) = j(l,j_o2_o)
     $                + j(l,j_o2_o1d)
     $                + a001(l)*cc(l,i_o)
     $                + c011(l)*cc(l,i_h)
     $                + cab003(l)*cc(l,i_ch3)
     $                + cab015(l)*cc(l,i_ch3o)
     $                + cab026(l)*cc(l,i_hco)
     $                + cab039(l)*cc(l,i_c2h5)
     $                + cab041(l)*cc(l,i_c2h5)
     $                + cab051(L)*cc(l,i_hoch2ch2o)
     $                + cab070(l)*cc(l,i_ch3co)
     $                + cab071(l)*cc(l,i_ch3co)
     $                + cab083(l)*cc(l,i_hcoco)
     $                + cab091(l)*cc(l,i_hoch2co)
     $                + cab092(l)*cc(l,i_hoch2co)
     $                + cl028(l)*cc(l,i_cl)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_h2) = c005(l)*cc(l,i_h)*cc(l,i_ho2)
     $                      + c018(l)*cc(l,i_h)*cc(l,i_h)
     $                      + b009(l)*cc(l,i_o1d)*cc(l,i_ch4)
     $                   + 0.17*cab005(l)*cc(l,i_ch3)*cc(l,i_o)    
     $                   + j(l,j_ch4_1ch2_h2)*cc(l,i_ch4)
     $                   + j(l,j_ch4_ch_h2_h)*cc(l,i_ch4)
     $                   + j(l,j_ch2o_co)*cc(l,i_hcho)
     $                   + cab027(l)*cc(l,i_hco)*cc(l,i_h)
     $                   + j(l,j_glyox_h2)*cc(l,i_glyox)
     $                   + cl041(l)*cc(l,i_hcl)*cc(l,i_h)

         loss(l,i_h2) = b003(l)*cc(l,i_o1d)
     $                + c010(l)*cc(l,i_oh)
     $                + cl008(l)*cc(l,i_cl)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h2o
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_h2o) = 
     $                   c006(l)*cc(l,i_h)*cc(l,i_ho2)
     $                 + c007(l)*cc(l,i_oh)*cc(l,i_ho2)
     $                 + c009(l)*cc(l,i_oh)*cc(l,i_h2o2)
     $                 + c010(l)*cc(l,i_oh)*cc(l,i_h2)
     $                 + c013(l)*cc(l,i_oh)*cc(l,i_oh)
     $                 + h004(l)*cc(l,i_h2o2)*hetero_ice
     $                 + cab001(l)*cc(l,i_ch4)*cc(l,i_oh)
     $                 + cab007(L)*cc(l,i_ch3o2)*cc(l,i_ho2)
     $                 + cab013(l)*cc(l,i_ch3oh)*cc(l,i_oh)
     $                 + cab014(l)*cc(l,i_ch3ooh)*cc(l,i_oh)
     $                 + cab018(l)*cc(l,i_hcho)*cc(l,i_oh)
     $                 + cab025(l)*cc(l,i_hco)*cc(l,i_oh)
     $                 + 0.3*cab029(l)*cc(l,i_hoch2o2)*cc(l,i_ho2)
     $                 + cab032(l)*cc(l,i_oh)*cc(l,i_hcooh)
     $                 + cab034(l)*cc(l,i_hoch2ooh)*cc(l,i_oh)
     $                 + cab035(l)*cc(l,i_hoch2oh)*cc(l,i_oh)
     $                 + cab036(l)*cc(l,i_c2h6)*cc(l,i_oh)
     $                 + cab045(l)*cc(l,i_c2h5ooh)*cc(l,i_oh)
     $                 + cab047(l)*cc(l,i_c2h5oh)*cc(l,i_oh)
     $                 + cab053(l)*cc(l,i_ethgly)*cc(l,i_oh)
     $                 + cab054(l)*cc(l,i_hyetho2h)*cc(l,i_oh)
     $                 + cab055(l)*cc(l,i_hyetho2h)*cc(l,i_oh)
     $                 + cab056(l)*cc(l,i_oh)*cc(l,i_hyetho2h)
     $                 + cab057(l)*cc(l,i_oh)*cc(l,i_ch3cho)
     $                 + cab058(l)*cc(l,i_oh)*cc(l,i_ch3cho)
     $                 + cab067(l)*cc(l,i_ch3cooh)*cc(l,i_oh)
     $                 + cab077(l)*cc(l,i_oh)*cc(l,i_ch3coooh)
     $                 + cab081(l)*cc(l,i_oh)*cc(l,i_glyox)
     $                 + cab085(l)*cc(l,i_hooch2cho)*cc(l,i_oh)
     $                 + cab086(l)*cc(l,i_hooch2cho)*cc(l,i_oh)
     $                 + cab087(l)*cc(l,i_hooch2cho)*cc(l,i_oh)
     $                 + cab088(l)*cc(l,i_oh)*cc(l,i_hoch2cho)
     $                 + cab089(l)*cc(l,i_oh)*cc(l,i_hoch2cho)
     $                 + cab099(l)*cc(l,i_oh)*cc(l,i_hcoco2h)
     $                 + cab100(l)*cc(l,i_oh)*Cc(l,i_hoch2co3h)
     $                 + cab102(l)*cc(l,i_oh)*cc(l,i_hcoco3h)
     $                 + cab103(l)*cc(l,i_oh)*cc(l,i_hcoco3h)
     $                 + cl014(l)*cc(l,i_hcl)*cc(l,i_oh)
     $                 + cl015(l)*cc(l,i_hocl)*cc(l,i_oh)
     $                 + cl027(l)*cc(l,i_hclo4)*cc(l,i_oh)

         loss(l,i_h2o) = j(l,j_h2o)
     $                 + b002(l)*cc(l,i_o1d)
     $                 + t001(l)
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h2o2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_h2o2) =
     $                    c008(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                  + c016(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                  + c017(l)*cc(l,i_oh)*cc(l,i_oh)
c    $                  + 0.5*h001(l)*cc(l,i_ho2)*hetero_ice
c    $                  + 0.5*h002(l)*cc(l,i_oh)*hetero_ice
c
         loss(l,i_h2o2) = j(l,j_h2o2)
     $                  + c009(l)*cc(l,i_oh)
     $                  + c012(l)*cc(l,i_o)
     $                  + h004(l)*hetero_ice
     $                  + h005(l)*hetero_dust
     $                  + cl011(l)*cc(L,i_cl)
c



C cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C c        ch4
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_ch4) = cab022(l)*cc(l,i_hco)*cc(l,i_ch3)
     $                       + j(l,j_ch3cho_ch4)*cc(l,i_ch3cho)

    
         loss(l,i_ch4) = b007(l)*cc(l,i_o1d)
     $                 + b008(l)*cc(l,i_o1d)
     $                 + b009(l)*cc(l,i_o1d)
     $                 + cab001(l)*cc(l,i_oh)
     $                 + cab002(l)*cc(l,i_o)
     $                 + j(l,j_ch4_1ch2_h2)
     $                 + j(l,j_ch4_ch3_h)
     $                 + j(l,j_ch4_3ch2_h_h)
     $                 + j(l,j_ch4_ch_h2_h)
     $                 + cl016(l)*cc(l,i_cl)
     
C         IF ( ( l <= 1 ) .and. (sza <= 95. ) ) THEN 

C           k_pseudo(l) =(1./( 31540000.))
C      $             -  loss(l,i_ch4)

C         loss(l,i_ch4) = loss(l,i_ch4) + k_pseudo(l)

C         else

        k_pseudo(l) = 0.



C         ENDIF 


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hox = h + oh + ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_hox) = 
     $                   2.*j(l,j_h2o)*cc(l,i_h2o)
     $                 + 2.*j(l,j_h2o2)*cc(l,i_h2o2)
     $                 + 2.*b002(l)*cc(l,i_o1d)*cc(l,i_h2o)
     $                 + 2.*b003(l)*cc(l,i_o1d)*cc(l,i_h2)
     $                 + 2.*c012(l)*cc(l,i_o)*cc(l,i_h2o2)
     $                 + 2.*t001(l)*cc(l,i_h2o)
     $                 + cab002(l)*cc(l,i_ch4)*cc(l,i_o)
     $                 + b007(l)*cc(l,i_ch4)*cc(l,i_o1d)
     $                 + b008(l)*cc(l,i_ch4)*cc(l,i_o1d)
     $                 + 0.956*cab004(l)*cc(l,i_ch3)*cc(l,i_o3)           
     $                 + cab005(l)*cc(l,i_ch3)*cc(l,i_o)
     $                 + cab015(l)*cc(l,i_ch3o)*cc(l,i_o2)
     $                 + 0.25*cab017(l)*cc(l,i_ch3o)*cc(l,i_o)
     $                 + cab020(l)*cc(l,i_hcho)*cc(l,i_o)
     $                 + j(l,j_ch4_ch3_h)*cc(l,i_ch4)
     $                 + 2.*j(l,j_ch4_3ch2_h_h)*cc(l,i_ch4)
     $                 + j(l,j_ch4_ch_h2_h)*cc(l,i_ch4)
     $                 + j(l,j_ch3o2h)*cc(l,i_ch3ooh)
     $                 + j(l,j_ch2o_hco)*cc(l,i_hcho)*2.
     $                 + j(l,j_ch3oh)*cc(l,i_ch3oh)
     $                 + cab021(l)*cc(l,i_hco)*cc(l,i_o)
     $                 + cab026(l)*cc(l,i_hco)*cc(l,i_o2)
     $                 + cab028(l)*cc(l,i_hoch2o2)
     $                 + cab030(l)*cc(l,i_hoch2o2)*ro2(l) 
     $                 + cab037(l)*cc(l,i_c2h6)*cc(l,i_o)
     $                 + cab041(l)*cc(l,i_c2h5)*cc(l,i_o2)
     $                 + cab044(l)*cc(l,i_c2h5o2)*ro2(l) *0.6
     $                 + cab051(L)*cc(l,i_hoch2ch2o)*cc(l,i_o2)
     $                 + cab052(l)*cc(l,i_hoch2ch2o)
     $                 + cab064(l)*cc(l,i_ch3choho2)
     $                 + cab066(L)*cc(l,i_ch3choho2)*ro2(l) 
     $                 + cab071(l)*cc(l,i_ch3co)*cc(l,i_o2)
     $                 + 0.6*cab078(l)*cc(l,i_hcoch2o2)*ro2(l) 
     $                 + cab080(l)*cc(l,i_hcoch2o2)*cc(l,i_ho2)
     $                 + cab084(l)*cc(l,i_o2)*cc(l,i_hcoco)
     $                 + 0.025*j(l,j_ch2o_hco)*cc(l,i_ch3coooh)
     $                 + 2.*j(l,j_ch3o2h)*cc(l,i_hoch2co3h)
     $                 + j(l,j_ch2o_co)*3.95*cc(l,i_hcoco2h)*2.
     $                 + j(l,j_ch3o2h)*cc(l,i_ch3chohooh)
     $                 + cab090(l)*cc(l,i_hochcho)
     $                 + cab092(l)*cc(l,i_hoch2co)*cc(l,i_o2)
     $                 + cab093(l)*cc(l,i_hoch2co3)*ro2(l) 
     $                 + cab095(l)*cc(l,i_ho2)*cc(l,i_hoch2co3)
     $                 + cab104(l)*cc(l,i_hcoco3)*ro2(l) 
     $                 + cab106(l)*cc(l,i_ho2)*cc(l,i_hcoco3)
     $                 + j(l,j_ch3o2h)*cc(l,i_hyetho2h)
     $                 + j(l,j_ch3cho_h)*cc(l,i_ch3cho)
     $                 + 2.*j(l,j_hoch2ooh)*cc(l,j_hoch2ooh)
     $                 + j(l,j_hoch2cho_oh)*cc(l,i_hoch2cho)
     $                 + j(l,j_ch3coooh)*cc(l,i_ch3coooh)
     $                 + 2.*( j(l,j_ch3o2h)
     $                 + j(l,j_hoch2cho_co)
     $                 + j(l,j_hoch2cho_hco)
     $                 + j(l,j_hoch2cho_oh) )*cc(l,i_hcoco3h)
     $                 + 2.*( j(l,j_ch3o2h)
     $                 + j(l,j_hoch2cho_co)
     $                 + j(l,j_hoch2cho_hco)
     $                 + j(l,j_hoch2cho_oh) )*cc(l,i_hoch2cho)
     $                 + cl008(l)*cc(l,i_h2)*cc(L,i_cl)     
     $                 + cl011(l)*cc(l,i_h2o2)*cc(l,i_cl)
     $                 + cl018(l)*cc(l,i_cl)*cc(l,i_ch3ooh)
     $                 + 0.88*cl039(l)*cc(l,i_hcl)*cc(l,i_o1d)
     $                 + cl040(l)*cc(l,i_hcl)*cc(l,i_o)
     $                 + cl041(l)*cc(l,i_hocl)*cc(l,i_o)
     $                 + cl054(l)*cc(l,i_clo4)*cc(l,i_hocl)
     $                 + cl055(l)*cc(l,i_hocl)*cc(l,i_cl)
     $                 + j(l,j_hcl)*cc(l,i_hcl)
     $                 + j(l,j_hocl)*cc(l,i_hocl)

          loss(l,i_hox) = 2.*c005(l)*cc(l,i_h)*cc(l,i_ho2)
     $                 + 2.*c006(l)*cc(l,i_h)*cc(l,i_ho2)
     $                 + 2.*c007(l)*cc(l,i_oh)*cc(l,i_ho2)
     $                 + 2.*c008(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                 + 2.*c013(l)*cc(l,i_oh)*cc(l,i_oh)
     $                 + 2.*c016(l)*cc(l,i_ho2)*cc(l,i_ho2)
     $                 + 2.*c017(l)*cc(l,i_oh)*cc(l,i_oh)
     $                 + 2.*c018(l)*cc(l,i_h)*cc(l,i_h)
     $                 + h001(l)*cc(l,i_ho2)*hetero_ice
     $                 + h002(l)*cc(l,i_oh)*hetero_ice
     $                 + h003(l)*cc(l,i_ho2)*hetero_dust
     $                 + cab001(l)*cc(l,i_ch4)*cc(l,i_oh)
     $                 + cab006(l)*cc(l,i_ch3o2)*cc(l,i_ho2)
     $                 + cab007(L)*cc(l,i_ch3o2)*cc(l,i_ho2)
     $                 + 0.15*cab013(l)*cc(l,i_ch3oh)*cc(l,i_oh)
     $                 + 0.6*cab014(l)*cc(l,i_ch3ooh)*cc(l,i_oh)
     $                 + cab018(l)*cc(l,i_hcho)*cc(l,i_oh)
     $                 + cab019(l)*cc(l,i_hcho)*cc(l,i_ho2)
     $                 + cab025(l)*cc(l,i_hco)*cc(l,i_oh)
     $                 + cab027(l)*cc(l,i_hco)*cc(l,i_h)
     $                 + 0.6*cab029(l)*cc(l,i_hoch2o2)*cc(l,i_ho2)
     $                 + cab034(l)*cc(l,i_hoch2ooh)*cc(l,i_oh)
     $                 + cab036(l)*cc(l,i_c2h6)*cc(l,i_oh)
     $                 + cab042(l)*cc(l,i_c2h5)*cc(l,i_h)
     $                 + cab043(l)*cc(l,i_c2h5o2)*cc(l,i_ho2)
     $                 + cab045(l)*cc(l,i_c2h5ooh)*cc(l,i_oh)
     $                 + 0.05*cab047(l)*cc(l,i_c2h5oh)*cc(l,i_oh)
     $                 + cab049(l)*cc(l,i_hoch2ch2o2)*cc(l,i_ho2)
     $                 + cab054(l)*cc(l,i_hyetho2h)*cc(l,i_oh)
     $                 + cab057(l)*cc(l,i_oh)*cc(l,i_ch3cho)
     $                 + cab058(l)*cc(l,i_oh)*cc(l,i_ch3cho)
     $                 + cab059(l)*cc(l,i_ho2)*cc(l,i_ch3cho)
     $                 + 0.8*cab065(l)*cc(l,i_ch3choho2)*cc(l,i_ho2)
     $                 + cab067(l)*cc(l,i_ch3cooh)*cc(l,i_oh)
     $                 + cab069(l)*cc(l,i_ch3chohooh)*cc(l,i_oh)
     $                 + cab073(l)*cc(l,i_ch3cooo)*cc(l,i_ho2)
     $                 + cab074(l)*cc(l,i_ch3cooo)*cc(l,i_ho2)
     $                 + cab077(l)*cc(l,i_oh)*cc(l,i_ch3coooh)
     $                 + cab079(l)*cc(l,i_hcoch2o2)*cc(l,i_ho2)
     $                 + cab080(l)*cc(l,i_hcoch2o2)*cc(l,i_ho2)
     $                 + cab081(l)*cc(l,i_oh)*cc(l,i_glyox)
     $                 + cab085(l)*cc(l,i_oh)*cc(l,i_hooch2cho)
     $                 + cab088(l)*cc(l,i_oh)*cc(l,i_hoch2cho)
     $                 + cab089(l)*cc(l,i_oh)*cc(l,i_hoch2cho)
     $                 + cab096(l)*cc(l,i_ho2)*cc(l,i_hoch2co3)
     $                 + cab097(l)*cc(l,i_ho2)*cc(l,i_hoch2co3)
     $                 + cab100(l)*cc(l,i_oh)*Cc(l,i_hoch2co3h)
     $                 + cab102(l)*cc(l,i_oh)*cc(l,i_hcoco3h)
     $                 + cab107(l)*cc(l,i_oh)*cc(l,i_ch3)
     $                 + cl009(l)*cc(l,i_ho2)*cc(l,i_cl)
     $                 + cl012(l)*cc(l,i_oh)*cc(l,i_clo)*0.06
     $                 + cl013(l)*cc(l,i_clo)*cc(l,i_ho2)
     $                 + cl014(l)*cc(l,i_hcl)*cc(l,i_oh)
     $                 + cl015(l)*cc(l,i_hocl)*cc(l,i_oh)
     $                 + cl027(l)*cc(l,i_hclo4)*cc(l,i_oh)
     $                 + cl032(l)*Cc(l,i_ch3ocl)*cc(l,i_oh)
     $                 + cl036(l)*cc(l,i_cl2)*cc(l,i_oh)
     $                 + cl037(l)*cc(l,i_cl2)*cc(l,i_h)
     $                 + cl041(l)*cc(l,i_hcl)*cc(l,i_h)
     $                 + cl048(l)*cc(l,i_clo3)*cc(l,i_oh)
     $                 + cl049(l)*cc(l,i_clo3)*cc(l,i_oh)
     $                 + cl058(l)*cc(l,i_oclo)*cc(l,i_oh)


      loss(l,i_hox) = loss(l,i_hox)/cc(l,i_hox)
c


      
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C c        ch3
C C cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


          production(l,i_ch3) = b007(l)*cc(l,i_ch4)*cc(l,i_o1d)
     $                        + cab001(l)*cc(l,i_ch4)*cc(l,i_oh)
     $                        + 0.51*cab002(l)*cc(l,i_ch4)*cc(l,i_o)
     $                        + j(l,j_ch4_ch3_h)*cc(l,i_ch4)
     $                        + cab042(l)*cc(l,i_c2h5)*cc(l,i_h)*2.
     $                        + 0.2*cab065(l)*cc(l,i_ch3choho2)
     $                                       *cc(l,i_ho2)
     $                        + cab066(l)*cc(l,i_ch3choho2)*ro2(l) 
     $                        + cab067(l)*cc(l,i_ch3cooh)*cc(l,i_oh)
     $                        + cab072(l)*cc(l,i_ch3cooo)*cc(l,i_ho2)
     $                        + cab075(l)*cc(l,i_ch3cooo)*ro2(l) 
     $                        + j(l,j_ch3o2h)*cc(l,i_ch3chohooh)
     $                        + j(l,j_ch3cho_ch3)*cc(l,i_ch3cho)
     $                        + j(l,j_ch3cooh)*cc(l,i_ch3cooh)
     $                        + j(l,j_ch3coooh)*cc(l,i_ch3coooh)
     $                        + 0.75*cab017(l)*cc(l,i_o)*cc(l,i_ch3o) 
     $                        + cl016(l)*cc(l,i_cl)*cc(l,i_ch4)

          loss(l,i_ch3) = cab003(l)*cc(l,i_o2)
     $                  + cab004(l)*cc(l,i_o3)
     $                  + cab005(l)*cc(l,i_o)
     $                  + cab022(L)*cc(l,i_hco)
     $                  + cab023(L)*cc(l,i_hco)
     $                  + cab038(l)*cc(l,i_ch3)
     $                  + cab107(l)*cc(l,i_oh)
     $                  + cl038(l)*cc(l,i_cl2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch3o2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          production(l,i_ch3o2) =  cab003(l)*cc(l,i_ch3)*cc(l,i_o2)
     $                         + 0.6*cab014(l)*cc(L,i_ch3ooh)*cc(l,i_oh)
     $                         + cab016(l)*cc(l,i_ch3o)*cc(l,i_o3)
     
	 
	 
          loss(l,i_ch3o2) = cab006(l)*cc(l,i_ho2)
     $                    + cab007(l)*cc(l,i_ho2)
     $                    + cab008(l)*ro2(l) 
     $                    + cab009(l)*ro2(l) 
     $                    + cab010(l)*cc(l,i_o3)
     $                    + cab011(l)*cc(l,i_oh)
     $                    + cab012(l)*cc(l,i_o)
     $                    + cl019(l)*cc(l,i_clo)
     $                    + cl020(l)*cc(l,i_clo)
     $                    + cl021(l)*cc(l,i_clo)
     $                    + cl022(l)*cc(l,i_cl)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch3ooh
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          
          production(l,i_ch3ooh) = cab006(l)*cc(l,i_ch3o2)*cc(l,i_ho2)

          loss(l,i_ch3ooh) = cab014(l)*cc(l,i_oh)
     $                     + j(l,j_ch3o2h) 
     $                     + cl018(l)*cc(l,i_cl)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch3oh
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         production(l,i_ch3oh) = 0.5*cab009(l)*cc(l,i_ch3o2)*ro2(l) 
     $                         + cab107(l)*cc(l,i_ch3)*cc(l,i_oh)
     $                         + j(l,j_hoch2cho_co)*cc(l,i_hoch2cho)
         
         
        loss(l,i_ch3oh) = cab013(l)*cc(l,i_oh)
     $                  + j(l,j_ch3oh)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch3o
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         production(l,i_ch3o) =  0.49*cab002(l)*cc(l,i_ch4)*cc(l,i_o)
     $                         + b008(l)*cc(l,i_ch4)*cc(l,i_o1d)
     $                         + 0.044*cab004(l)*cc(l,i_ch3)*cc(l,i_o3)
     $                         + cab008(l)*cc(l,i_ch3o2)*ro2(l) 
     $                         + cab010(l)*cc(l,i_ch3o2)*cc(l,i_o3)
     $                         + cab011(l)*cc(l,i_ch3o2)*cc(l,i_oh)
     $                         + cab012(l)*cc(l,i_ch3o2)*cc(l,i_o)
     $                         + 0.15*cab013(l)*cc(l,i_ch3oh)*cc(l,i_oh)
     $                         + j(l,j_ch3o2h)*cc(l,i_ch3ooh)
     $                         + j(l,j_ch3oh)*cc(l,i_ch3oh)
     $                         + j(l,j_hoch2cho_hco)*cc(l,i_hoch2cho)
     $                         + cl019(l)*cc(l,i_clo)*cc(l,i_ch3o2)
     $                         + cl021(l)*cc(l,i_clo)*cc(L,i_ch3o2)
     $                         + cl022(l)*cc(l,i_cl)*cc(l,i_ch3o2)*0.5
     $                         + cl029(l)*cc(l,i_cl)*cc(l,i_ch3ocl) 
     $                         + cl032(l)*Cc(l,i_ch3ocl)*cc(l,i_oh)

         loss(l,i_ch3o) = cab015(l)*cc(l,i_o2)
     $                  + cab016(l)*cc(l,i_o3)
     $                  + cab017(l)*cc(l,i_o)
     
    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hcho (formaldehyde)
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          production(l,i_hcho) = b009(l)*cc(l,i_ch4)*cc(l,i_o1d)
     $                   + 0.956*cab004(l)*cc(l,i_ch3)*cc(l,i_o3)
     $                   + 0.83*cab005(l)*cc(l,i_ch3)*cc(l,i_o)     
     $                   + cab007(L)*cc(l,i_ch3o2)*cc(l,i_ho2)
     $                   + 0.5*cab009(l)*cc(l,i_ch3o2)*ro2(l) 
     $                   + 0.85*cab013(l)*cc(l,i_ch3oh)*cc(l,i_oh)
     $                   + 0.4*cab014(l)*cc(L,i_ch3ooh)*cc(l,i_oh)
     $                   + cab015(l)*cc(l,i_ch3o)*cc(l,i_o2)
     $                   + 0.25*cab017(l)*cc(l,i_ch3o)*cc(l,i_o)
     $                   + cab024(l)*cc(l,i_hco)*cc(l,i_hco)
     $                   + cab028(l)*cc(l,i_hoch2o2)
     $                   + 2.*cab052(l)*cc(l,i_hoch2ch2o)
     $                   + cab061(l)*cc(l,i_ch2choh)*cc(l,i_oh)
     $                   + cab071(l)*cc(l,i_ch3co)*cc(l,i_o2)
     $                   + 0.6*cab078(l)*cc(l,i_hcoch2o2)*ro2(l) 
     $                   + cab080(l)*cc(l,i_hcoch2o2)*cc(l,i_ho2)
     $                   + cab086(l)*cc(l,i_oh)*cc(l,i_hooch2cho)
     $                   + j(l,j_ch3o2h)*cc(l,i_hoch2co3h)
     $                   + cab092(l)*cc(l,i_hoch2co)*cc(l,i_o2)
     $                   + cab093(l)*cc(l,i_hoch2co3)*ro2(l) 
     $                   + cab095(l)*cc(l,i_hoch2co3)*cc(l,i_ho2)
     $                   + 0.09*cab098(l)*cc(l,i_hoch2co2h)*cc(l,i_oh)
     $                   + j(l,j_glyox_hcho)*cc(l,i_glyox)
     $                   + ( j(l,j_ch3o2h)
     $                   + j(l,j_hoch2cho_co)
     $                   + j(l,j_hoch2cho_hco)
     $                   + j(l,j_hoch2cho_oh) )*cc(l,i_hoch2cho)
     $                   + cl018(l)*cc(l,i_cl)*cc(l,i_ch3ooh)

          loss(l,i_hcho) = cab018(l)*cc(l,i_oh)
     $                   + cab019(l)*cc(l,i_ho2)
     $                   + cab020(l)*cc(l,i_o)
     $                   + j(l,j_ch2o_co)
     $                   + j(l,j_ch2o_hco)
     $                   + cl017(l)*cc(l,i_cl)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hco 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        production(l,i_hco) = cab018(l)*cc(l,i_hcho)*cc(l,i_oh)
     $                      + cab020(l)*cc(l,i_hcho)*cc(l,i_o)
     $                      + j(l,j_ch2o_co)*cc(l,i_hcho)
     $                      + cab082(l)*cc(l,i_hcoco)
     $                      + j(l,j_ch3cho_ch3)*cc(l,i_ch3cho)
     $                      + j(l,j_hoch2cho_hco)*cc(l,i_hoch2cho)
     $                      + 2.*j(l,j_glyox_hco)*cc(l,i_glyox)
     $                      + cl017(l)*cc(l,i_cl)*cc(l,i_hcho)

        loss(l,i_hco) = cab021(l)*cc(l,i_o)
     $                + cab022(l)*cc(l,i_ch3)
     $                + cab023(l)*cc(l,i_ch3)
     $                + cab024(l)*cc(l,i_hco)
     $                + cab025(l)*cc(l,i_oh)
     $                + cab026(l)*cc(l,i_o2)
     $                + cab027(l)*cc(l,i_h)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hoch2o2
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        production(l,i_hoch2o2) = cab019(l)*cc(l,i_hcho)*cc(l,i_ho2)
     $                          + cab033(l)*cc(l,i_hoch2ooh)*cc(l,i_oh)
     
        loss(l,i_hoch2o2) = cab028(l)
     $                    + cab029(l)*cc(l,i_ho2)
     $                    + cab030(l)*ro2(l) 
     $                    + cab031(l)*ro2(l) 
        

ccccccccccccccccccccccccccccccccc
c        hcooh
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           production(l,i_hcooh) = 
     $              0.5*cab029(l)*cc(l,i_hoch2o2)*cc(l,i_ho2)
     $            + cab030(l)*cc(l,i_hoch2o2)*ro2(l) 
     $            + 0.5*cab031(l)*cc(l,i_hoch2o2)*ro2(l) 
     $            + cab034(l)*cc(l,i_hoch2ooh)*cc(l,i_oh)
     $            + cab035(l)*cc(l,i_hoch2oh)*cc(l,i_oh)
     $            + cab061(l)*cc(L,i_ch2choh)*cc(l,i_oh)
     $            + 0.2*cab065(l)*cc(l,i_ch3choho2)*cc(L,i_ho2)
     $            + cab066(L)*cc(l,i_ch3choho2)*ro2(l) 

          loss(l,i_hcooh) = cab032(l)*cc(l,i_oh)
     $                    + cl043(l)*cc(l,i_hcooh)
          


     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hoch2oh
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        production(l,i_hoch2oh) = 
     $                   0.5*cab031(l)*cc(l,i_hoch2o2)*ro2(l) 
        
        loss(l,i_hoch2oh) = cab035(l)*cc(l,i_oh)
        

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hoch2ooh
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        production(l,i_hoch2ooh) = 
     $                0.5*cab029(l)*cc(l,i_hoch2o2)*cc(l,i_ho2)
        
        loss(l,i_hoch2ooh) = cab033(l)*cc(l,i_oh)
     $                     + cab034(l)*cc(l,i_oh)
     $                     + j(l,j_hoch2ooh)

          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        c2h6
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          production(l, i_c2h6) = cab038(l)*cc(l,i_ch3)*cc(l,i_ch3)
     $                          + cab040(l)*cc(l,i_c2h5)*cc(l,i_c2h5)
          
          loss(l, i_c2h6) = j(l,j_c2h6)
     $                    + cab036(l)*cc(l,i_oh)
     $                    + cab037(l)*cc(l,i_o)
         
          

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        c2h5
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        production(l,i_c2h5) = cab036(l)*cc(l,i_c2h6)*cc(l,i_oh)
     $                       + cab037(l)*cc(l,i_c2h6)*cc(l,i_o)
     
        loss(l,i_c2h5) = cab039(l)*cc(l,i_o2)
     $                 + cab040(l)*cc(l,i_c2h5)
     $                 + cab041(l)*cc(l,i_o2)
     $                 + cab042(l)*cc(l,i_h)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        c2h5o2
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l, i_c2h5o2) = 
     $              cab039(l)*cc(l,i_c2h5)*cc(l,i_o2)
     $            + cab045(l)*cc(l,i_c2h5ooh)*cc(l,i_oh)
     $            + 0.95*cab047(l)*cc(l,i_c2h5oh)*cc(l,i_oh)
     
          loss(l,i_c2h5o2) = cab043(l)*cc(l,i_ho2)
     $                     + cab044(l)*ro2(l) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        c2h5ooh
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         production(l, i_c2h5ooh) =
     $                   cab043(l)*cc(l,i_ho2)*cc(l,i_c2h5o2)
         
         
         loss(l,i_c2h5ooh) = cab045(l)*cc(l,i_oh)
     $                     + cab046(l)*cc(l,i_oh)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        c2h5oh
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l, i_c2h5oh) =    
     $                  + 0.2*cab044(l)*cc(l,i_c2h5o2)*ro2(l) 


          loss(l,i_c2h5oh) = cab047(l)*cc(l,i_oh)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hoch2ch2o2
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hoch2ch2o2) = 
     $              0.05*cab047(l)*cc(l,i_c2h5oh)*cc(l,i_oh)
     $            + cab054(l)*cc(l,i_hyetho2h)*cc(l,i_oh)

          loss(l,i_hoch2ch2o2) = cab048(l)*ro2(l) 
     $                         + cab049(l)*cc(l,i_ho2)
     $                         + cab050(l)*cc(l,i_ho2)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hoch2ch2o
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hoch2ch2o) = 
     $              0.6*cab048(l)*cc(l,i_hoch2ch2o2)*ro2(l) 
     $            + cab050(l)*cc(l,i_hoch2ch2o2)*cc(l,i_ho2)
     $            + j(l,j_ch3o2h)*cc(l,i_hyetho2h)
     
          loss(l,i_hoch2ch2o) = cab051(l)*cc(l,i_o2)
     $                        + cab052(l)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ethgly
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_ethgly) = 
     $              0.2*cab048(l)*cc(l,i_hoch2ch2o2)*ro2(l) 


          loss(l,i_ethgly) = cab053(l)*cc(l,i_oh)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hyetho2h
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hyetho2h) = 
     $              cab049(l)*cc(l,i_hoch2ch2o2)*cc(l,i_ho2)

          loss(l,i_hyetho2h) = cab054(l)*cc(l,i_oh)
     $                       + cab055(l)*cc(l,i_oh)
     $                       + cab056(l)*cc(l,i_oh)
     $                       + j(l,j_ch3o2h)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch3cho
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_ch3cho) = 
     $                    cab023(l)*cc(l,i_hco)*cc(l,i_ch3)
     $                  + 0.8*cab044(l)*cc(l,i_c2h5o2)*ro2(l) 
     $                  + cab046(l)*cc(l,i_c2h5ooh)*cc(l,i_oh)
     $                  + cab063(l)*cc(l,i_ch2choh)*cc(l,i_hcooh)
     $                  + cab064(l)*cc(l,i_ch3choho2)
     
          loss(l,i_ch3cho) = cab057(l)*cc(l,i_oh)
     $                     + cab058(l)*cc(l,i_oh)
     $                     + cab059(l)*cc(l,i_ho2)
     $                     + cab060(l)*cc(l,i_hcooh)
     $                     + j(l,j_ch3cho_ch3)
     $                     + j(l,j_ch3cho_ch4)
     $                     + j(l,j_ch3cho_h)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch2choh
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_ch2choh) = 
     $              cab060(l)*cc(l,i_ch3cho)*cc(l,i_hcooh)

          loss(l,i_ch2choh) = cab061(l)*cc(l,i_oh)
     $                      + cab062(l)*cc(l,i_oh)
     $                      + cab063(l)*cc(l,i_hcooh)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch3choho2
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_ch3choho2) = 
     $       cab059(l)*cc(l,i_ch3cho)*cc(l,i_ho2)
     $     + cab069(l)*cc(l,i_ch3chohooh)*cc(l,i_oh)

          loss(l,i_ch3choho2) = cab064(l)
     $                        + cab065(l)*cc(l,i_ho2)
     $                        + cab066(l)*ro2(l) 
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch3cooh
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_ch3cooh) = 
     $         0.3*cab065(l)*cc(l,i_ch3choho2)*cc(l,i_ho2)
     $       + cab068(l)*cc(l,i_ch3chohooh)*cc(l,i_oh)
     $       + cab074(l)*cc(l,i_ch3cooo)*cc(l,i_ho2)
     $       + cab076(l)*cc(l,i_ch3cooo)*ro2(l) 

          loss(l,i_ch3cooh) = cab067(l)*cc(l,i_oh)
     $                      + j(l,j_ch3cooh)
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch3chohooh
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_ch3chohooh) = 
     $              0.5*cab065(l)*cc(l,i_ho2)*cc(l,i_ch3choho2)

          loss(l,i_ch3chohooh) = cab068(l)*cc(l,i_oh)
     $                         + cab069(l)*cc(l,i_oh)
     $                         + j(l,j_ch3o2h)
    
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch3co
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_ch3co) = 
     $         cab057(l)*cc(l,i_oh)*cc(l,i_ch3cho)
     $       + j(l,j_ch3cho_h)*cc(l,i_ch3cho)

          loss(l,i_ch3co) =  cab070(l)*cc(l,i_o2)
     $                    +  cab071(l)*cc(l,i_o2)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch3cooo
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_ch3cooo) = 
     $              cab070(l)*cc(l,i_o2)*cc(l,i_ch3co)
     $            + cab077(l)*cc(l,i_oh)*cc(l,i_ch3coooh)

          loss(l,i_ch3cooo) = cab072(l)*cc(l,i_ho2)
     $                      + cab073(l)*cc(l,i_ho2)
     $                      + cab074(l)*cc(l,i_ho2)
     $                      + cab075(l)*ro2(l) 
     $                      + cab076(l)*ro2(l) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch3coooh
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_ch3coooh) = 
     $              cab073(l)*cc(l,i_ch3cooo)*cc(l,i_Ho2)

          loss(l,i_ch3coooh) = cab077(l)*cc(l,i_oh)
     $                       + j(l,j_ch3coooh)

          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        glyox
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_glyox) =    
     $              0.2*cab078(l)*cc(l,i_hcoch2o2)*ro2(l) 
     $            + cab087(l)*cc(l,i_oh)*cc(l,i_hooch2cho)
     $            + cab090(l)*cc(l,i_hochcho)
     
          loss(l,i_glyox) = cab081(l)*cc(l,i_oh)
     $                    + j(l,j_glyox_h2)
     $                    + j(l,j_glyox_hcho)
     $                    + j(l,j_glyox_hco)
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hcoco
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hcoco) = cab081(l)*cc(l,i_oh)*cc(l,i_glyox)

          loss(l,i_hcoco) = cab082(l)
     $                    + cab083(l)*cc(l,i_o2)
     $                    + cab084(l)*cc(l,i_o2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hooch2cho
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hooch2cho) = 
     $              cab056(l)*cc(l,i_oh)*cc(l,i_hyetho2h)
     $            + cab079(l)*cc(l,i_hcoch2o2)*cc(l,i_ho2)
     
          loss(l,i_hooch2cho) = cab085(l)*cc(l,i_oh)
     $                        + cab086(l)*cc(l,i_oh)
     $                        + cab087(l)*cc(l,i_oh)
     $                        + j(l,j_ch3o2h)
     $                        + j(l,j_hoch2cho_co)
     $                        + j(l,j_hoch2cho_hco)
     $                        + j(l,j_hoch2cho_oh) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hoch2cho
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hoch2cho) = 
     $              0.2*cab048(l)*cc(l,i_hoch2ch2o2)*ro2(l) 
     $            + cab051(L)*cc(l,i_hoch2ch2o)*cc(l,i_o2)
     $            + cab053(l)*cc(l,i_ethgly)*cc(l,i_oh)
     $            + cab055(l)*cc(l,i_hyetho2h)*cc(l,i_oh)
     $            + cab062(l)*cc(l,i_ch2choh)*cc(l,i_oh)
     $            + 0.2*cab078(l)*cc(l,i_hcoch2o2)*ro2(l) 
     
          loss(l,i_hoch2cho) = cab088(l)*cc(l,i_oh)
     $                       + cab089(l)*cc(l,i_oh)
     $                       + j(l,j_hoch2cho_co)
     $                       + j(l,j_hoch2cho_hco)
     $                       + j(l,j_hoch2cho_oh)
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hochcho
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hochcho) =
     $              cab089(l)*cc(l,i_oh)*cc(l,i_hoch2cho)

          loss(l,i_hochcho) = cab090(l)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hoch2co
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hoch2co) = 
     $              cab088(l)*cc(l,i_oh)*cc(l,i_hoch2cho)

          loss(l,i_hoch2co) = cab091(l)*cc(l,i_o2)
     $                      + cab092(l)*cc(l,i_o2)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hoch2co3
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hoch2co3) = 
     $              cab091(l)*cc(l,i_o2)*cc(l,i_hoch2co)
     $            + cab100(l)*cc(l,i_oh)*Cc(l,i_hoch2co3h)

          loss(l,i_hoch2co3) = 
     $              cab093(l)*ro2(l) 
     $            + cab094(l)*ro2(l) 
     $            + cab095(l)*cc(l,i_ho2)
     $            + cab096(l)*cc(l,i_ho2)
     $            + cab097(l)*cc(l,i_ho2)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hoch2co2h 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hoch2co2h) = 
     $              cab094(l)*cc(l,i_hoch2co3)*ro2(l) 
     $            + cab097(l)*cc(l,i_hoch2co3)*cc(l,i_ho2)
     
          loss(l,i_hoch2co2h) = cab098(l)*cc(l,i_oh)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hcoco2h
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hcoco2h) = 
     $         0.91*cab098(l)*cc(l,i_hoch2co2h)*cc(l,i_oh)
     $       + cab105(l)*cc(l,i_hcoco3)*ro2(l) 

          loss(l,i_hcoco2h) = j(l,j_ch2o_co)*3.95
     $                      + cab099(l)*cc(l,i_oh)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hoch2co3h 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hoch2co3h) =
     $         cab096(l)*cc(l,i_hoch2co3)*cc(l,i_ho2)

          loss(l,i_hoch2co3h) = j(l,j_ch3o2h)
     $                        + cab100(l)*cc(l,i_oh)
     $                        + cab101(l)*cc(l,i_oh)
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hcoco3h
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hcoco3h) =
     $         cab101(l)*cc(l,i_hoch2co3h)*cc(l,i_oh)
     
          loss(l,i_hcoco3h) = cab102(l)*cc(l,i_oh)
     $                      + cab103(l)*cc(l,i_oh)
     $                      + j(l,j_ch3o2h)
     $                      + j(l,j_hoch2cho_co)
     $                      + j(l,j_hoch2cho_hco)
     $                      + j(l,j_hoch2cho_oh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hcoco3
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          production(l,i_hcoco3) = 
     $              cab083(l)*cc(l,i_hcoco)*cc(l,i_o2)
     $            + cab102(l)*cc(l,i_hcoco3h)*cc(l,i_oh)
     
          loss(l,i_hcoco3) = cab104(L)*ro2(l) 
     $                     + cab105(l)*ro2(l) 
     $                     + cab106(l)*cc(l,i_ho2)
     
  
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Chlorine Chemistry
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ClOx = Cl + ClO 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      production(l,i_clox) = cl007(l)*cc(l,i_cl2o2)
     $                     + cl014(l)*cc(l,i_hcl)*cc(l,i_oh)
     $                     + cl015(l)*cc(l,i_hocl)*cc(l,i_oh)
     $                     + cl024(l)*cc(l,i_cloo)*cc(l,i_cl)
     $                     + cl026(l)*cc(l,i_cloo)
     $                     + cl033(l)*cc(l,i_clco)
     $                     + cl035(l)*cc(l,i_o1d)*cc(l,i_cl2)*2.
     $                     + cl036(l)*cc(l,i_cl2)*cc(L,i_oh)
     $                     + cl037(l)*cc(l,i_cl2)*cc(L,i_h)
     $                     + cl038(l)*cc(l,i_cl2)*cc(l,i_ch3)
     $                     + 0.88*cl039(l)*cc(l,i_hcl)*cc(l,i_o1d)
     $                     + cl040(l)*cc(l,i_hcl)*cc(l,i_o)
     $                     + cl041(l)*cc(l,i_hcl)*cc(l,i_h)
     $                     + cl042(l)*cc(l,i_hocl)*cc(l,i_o)
     $                     + cl057(l)*cc(l,i_oclo)*cc(l,i_o)
     $                     + cl059(l)*cc(l,i_oclo)*cc(l,i_cl)
     $                     + j(l,j_cl2)*cc(l,i_cl2)*2.
     $                     + j(l,j_cloo)*cc(l,i_cloo)
     $                     + j(l,j_oclo)*cc(l,i_oclo)
     $                     + j(l,j_cl2o2)*cc(l,i_cl2o2)
     $                     + j(l,j_hcl)*cc(l,i_hcl)
     $                     + j(l,j_hocl)*cc(l,i_hocl)

C       loss(l,i_clox) = 2.*cl003(l)*cc(l,i_clo) 
C      $               + cl005(l)*cc(l,i_clo)
C      $               + 2.*cl006(l)*cc(l,i_clo)
C      $               + cl008(l)*cc(l,i_h2)
C      $               + cl009(l)*cc(l,i_ho2)
C      $               + cl011(l)*cc(l,i_h2o2)
C      $               + 0.06*cl012(l)*cc(l,i_oh)
C      $               + cl013(l)*cc(l,i_ho2)
C      $               + cl016(l)*cc(l,i_ch4)
C      $               + cl017(l)*cc(l,i_hcho)
C      $               + cl018(l)*cc(l,i_ch3ooh)
C      $               + cl019(l)*cc(l,i_ch3o2)
C      $               + cl020(l)*cc(l,i_ch3o2)
C      $               + cl021(l)*cc(l,i_ch3o2)
C      $               + 0.5*cl022(l)*cc(l,i_ch3o2)
C      $               + cl023(l)*Cc(l,i_co)
C      $               + cl025(l)*cc(l,i_cloo)
C      $               + cl028(l)*Cc(l,i_o2)
C      $               + cl029(l)*cc(l,i_ch3ocl)
C      $               + cl030(l)*cc(l,i_cl2o2)
C      $               + cl031(l)*cc(l,i_ch3ocl)
C      $               + cl043(l)*cc(l,i_hcooh)
C      $               + cl044(l)*cc(l,i_o3)
C      $               + cl045(l)*cc(l,i_clo3)
C      $               + cl046(l)*cc(l,i_clo3)
C      $               + cl047(l)*cc(l,i_clo3)
C      $               + cl055(l)*cc(l,i_hocl)


      loss(l,i_clox) = 2.*cl003(l)*cc(l,i_clo)*cc(l,i_clo) 
     $               + cl005(l)*cc(l,i_clo)*cc(l,i_clo)
     $               + 2.*cl006(l)*cc(l,i_clo)*cc(l,i_clo)
     $               + cl008(l)*cc(l,i_h2)*cc(l,i_cl)
     $               + cl009(l)*cc(l,i_ho2)*cc(l,i_cl)
     $               + cl011(l)*cc(l,i_h2o2)*cc(l,i_cl)
     $               + 0.06*cl012(l)*cc(l,i_oh)*cc(l,i_clo)
     $               + cl013(l)*cc(l,i_ho2)*cc(l,i_clo)
     $               + cl016(l)*cc(l,i_ch4)*cc(l,i_cl)
     $               + cl017(l)*cc(l,i_hcho)*cc(l,i_cl)
     $               + cl018(l)*cc(l,i_ch3ooh)*cc(l,i_cl)
     $               + cl019(l)*cc(l,i_ch3o2)*cc(l,i_clo)
     $               + cl020(l)*cc(l,i_ch3o2)*cc(l,i_clo)
     $               + cl021(l)*cc(l,i_ch3o2)*cc(l,i_clo)
     $               + 0.5*cl022(l)*cc(l,i_ch3o2)*cc(l,i_clo)
     $               + cl023(l)*Cc(l,i_co)*cc(l,i_cl)
     $               + cl025(l)*cc(l,i_cloo)*cc(l,i_cl)
     $               + cl028(l)*Cc(l,i_o2)*cc(l,i_cl)
     $               + cl029(l)*cc(l,i_ch3ocl)*cc(l,i_cl)
     $               + cl030(l)*cc(l,i_cl2o2)*cc(l,i_cl)
     $               + cl031(l)*cc(l,i_ch3ocl)*cc(l,i_cl)
     $               + cl043(l)*cc(l,i_hcooh)*cc(l,i_cl)
     $               + cl044(l)*cc(l,i_o3)*cc(l,i_cl)
     $               + cl045(l)*cc(l,i_clo3)*cc(l,i_clo)
     $               + cl046(l)*cc(l,i_clo3)*cc(l,i_clO)
     $               + cl047(l)*cc(l,i_clo3)*cc(l,i_clo)
     $               + cl055(l)*cc(l,i_hocl)*cc(l,i_cl)

      loss(l,i_clox) = loss(l,i_clox)/MAX( cc(l,i_clox), 1.e-30*dens(L))

      ! Suppression Attempt
C         IF ( cc(l,i_cl) <= 1.e-30*dens(l) ) cc(l,i_cl) = 0.
C         IF ( cc(l,i_clo) <= 1.e-30*dens(l) ) cc(l,i_clo) = 0.

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     HCl 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        production(l,i_hcl) = cl008(l)*cc(l,i_cl)*cc(l,i_h2)
     $                      + cl009(l)*cc(l,i_cl)*cc(l,i_ho2)
     $                      + cl011(l)*cc(l,i_cl)*cc(l,i_h2o2)
     $                      + cl012(l)*cc(l,i_clo)*cc(l,i_oh)*0.06
     $                      + cl016(l)*cc(l,i_cl)*cc(l,i_ch4)
     $                      + cl017(l)*cc(l,i_cl)*cc(l,i_hcho)
     $                      + cl018(l)*cc(l,i_cl)*cc(l,i_ch3ooh)
     $                      + cl022(l)*cc(l,i_cl)*cc(l,i_ch3o2)*0.5
     $                      + cl031(l)*cc(l,i_cl)*cc(l,i_ch3ocl)
     $                      + cl037(l)*cc(l,i_cl2)*cc(l,i_h)
     $                      + cl043(l)*cc(l,i_cl)*cc(l,i_hcooh)
     $                      + cl056(l)*cc(l,i_cl)*Cc(l,i_hocl)

        loss(l,i_hcl) = j(l,j_hcl) 
     $                + cl014(l)*cc(l,i_oh)
     $                + cl039(l)*cc(l,i_o1d)*0.88
     $                + cl040(l)*cc(l,i_o) 
     $                + cl041(l)*cc(l,i_h) 


Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Cl2 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      production(l,i_cl2) = cl003(l)*cc(l,i_clo)*cc(l,i_clo)
     $                    + cl025(l)*cc(l,i_cloo)*cc(l,i_cl)
     $                    + cl029(l)*cc(l,i_cl)*cc(l,i_ch3ocl) 
     $                    + cl030(l)*cc(l,i_cl2o2)*cc(l,i_cl)
     $                    + cl055(l)*cc(l,i_hocl)*cc(l,i_cl)

      loss(l,i_cl2) = cl035(l)*cc(l,i_o1d)
     $              + cl036(l)*cc(l,i_oh)
     $              + cl037(l)*cc(l,i_h)
     $              + cl038(l)*cc(l,i_ch3)
     $              + j(l,j_cl2)

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ClOO 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      production(l,i_cloo) = cl019(l)*cc(l,i_clo)*cc(l,i_ch3o2)
     $                     + cl028(l)*cc(l,i_cl)*cc(l,i_o2)
     $                     + cl030(l)*cc(l,i_cl2o2)*cc(l,i_cl)
     $                     + cl045(l)*cc(l,i_clo)*cc(l,i_clo3)
     $                     + j(l,j_cl2o2)*cc(l,i_cl2o2)

      loss(l,i_cloo) = cl024(l)*cc(l,i_cl)
     $               + cl025(l)*cc(l,i_cl)
     $               + cl026(l)
     $               + j(l,j_cloo) 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     HOCl 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      production(l,i_hocl) = cl013(l)*cc(l,i_clo)*cc(l,i_ho2)
     $                     + cl032(l)*Cc(l,i_ch3ocl)*cc(l,i_oh)
     $                     + cl036(l)*cc(l,i_cl2)*cc(l,i_oh)
     $                     + cl058(l)*cc(l,i_oclo)*Cc(l,i_oh)

      loss(l,i_hocl) = cl015(l)*cc(l,i_oh)
     $               + cl042(l)*cc(l,i_o)
     $               + cl054(l)*cc(l,i_clo4) 
     $               + cl055(l)*cc(l,i_cl)
     $               + cl056(l)*cc(l,i_cl)
     $              + j(l,j_hocl) 

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Cl2O2 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      production(l,i_cl2o2) = cl006(l)*cc(l,i_clo)*cc(l,i_clo)
     $                      

      loss(l,i_cl2o2) = cl007(l)
     $                + cl030(l)*cc(l,i_cl)
     $                + j(l,j_cl2o2)

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OClO
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      production(l,i_oclo) = cl005(l)*cc(l,i_clo)*Cc(l,i_clo)
     $                     + cl021(l)*cc(l,i_clo)*Cc(l,i_ch3o2)
     $                     + cl045(l)*cc(l,i_clo)*cc(l,i_clo3)
     $                     + 2.*cl046(l)*cc(l,i_clo)*cc(l,i_clo3)
     $                     + cl050(l)*cc(l,i_clo3)*cc(l,i_oh)


      loss(l,i_oclo) = cl051(l)*cc(l,i_o)
     $               + cl052(l)*cc(l,i_o3)
     $               + cl057(l)*cc(l,i_o)
     $               + cl058(l)*cc(l,i_oh)
     $               + cl059(l)*cc(l,i_cl)
     $               + j(l,j_oclo)

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CH3OCl
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      production(l,i_ch3ocl) = cl020(l)*cc(l,i_clo)*cc(l,i_ch3o2)

      loss(l,i_ch3ocl) = cl029(l)*cc(l,i_cl)
     $                 + cl031(l)*cc(l,i_cl)
     $                 + cl032(l)*cc(l,i_oh)

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ClCO
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      production(l,i_clco) = cl023(l)*cc(l,i_co)*cc(l,i_cl)

      loss(l,i_clco) = cl033(l)

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ClO3
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      production(l,i_clo3) = cl044(l)*cc(l,i_cl)*cc(l,i_o3)
     $                     + cl051(l)*cc(l,i_oclo)*cc(l,i_o)
     $                     + cl052(l)*cc(l,i_oclo)*cc(l,i_o3) 
     $                     + cl053(l)*cc(l,i_clo4)*cc(l,i_cl)

      loss(l,i_clo3) = cl045(l)*cc(l,i_clo)
     $               + cl046(l)*cc(l,i_clo)
     $               + cl047(l)*cc(l,i_clo)
     $               + cl048(l)*cc(l,i_oh)
     $               + cl049(l)*cc(l,i_oh)
     $               + cl050(l)*cc(l,i_oh)


Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     HClO4
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      production(l,i_hclo4) = cl048(l)*cc(l,i_clo3)*cc(l,i_oh)
     $                      + cl049(l)*cc(l,i_clo3)*cc(l,i_oh)
     $                      + cl054(l)*cc(l,i_clo4)*cc(l,i_hocl)

      loss(l,i_hclo4) = cl027(l)*cc(l,i_oh)


Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ClO4
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      production(l,i_clo4) = cl027(l)*cc(l,i_hclo4)*cc(l,i_oh)

      loss(l,i_clo4) = cl053(l)*cc(l,i_cl)
     $               + cl054(l)*cc(l,i_hocl)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c      Tangent Linear Model - Inorganic/Long-lived species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
C           IF (tlm_on) THEN 
          
          
C       cc_hox_next = (cc0(l,i_hox) + production(l,i_hox)*dt)
C      $            /(1. + loss(l,i_hox)*dt)

C                call tlm_sibem(iter,niter, l, dens(l), sza,
C      $         dt, ptimestep,  
C      $         nesp, cc(l,:), cc0(l,:), cc_hox_next,
C      $         j(l,:), production(l,:), loss(l,:),  
C      $         0.,!methane_enhancement,
C      $         a001(l), a002(l), a003(l),  
C      $         b001(l), b002(l), b003(l), b004(l), b005(l), b006(l),  
C      $         b007(l), b008(l), b009(l),  
C      $         c001(l), c002(l), c003(l), c004(l), c005(l), c006(l),  
C      $         c007(l), c008(l), c009(l), c010(l), c011(l), c012(l),  
C      $         c013(l), c014(l), c015(l), c016(l), c017(l), c018(l),  
C      $         d001(l), d002(l), d003(l),
C      $         e001(l), e002(l), e003(l),  
C      $         cab001(l), cab002(l), cab003(l), cab004(l), cab005(l),  
C      $         cab006(l), cab007(l), cab008(l), cab009(l), cab010(l),  
C      $         cab011(l), cab012(l), cab013(l), cab014(l), cab015(l),  
C      $         cab016(l), cab017(l), cab018(l), cab019(l), cab020(l),  
C      $         cab021(l), cab022(l), cab023(l), cab024(l), cab025(l),  
C      $         cab026(l), cab027(l), cab028(l), cab029(l), cab030(l),  
C      $         cab031(l), cab032(l), cab033(l), cab034(l), cab035(l),  
C      $         cab036(l), cab037(l), cab038(l), cab039(l), cab040(l),  
C      $         cab041(l), cab042(l), cab043(l), cab044(l), cab045(l),  
C      $         cab046(l), cab047(l), cab048(l), cab049(l), cab050(l),  
C      $         cab051(l), cab052(l), cab053(l), cab054(l), cab055(l),  
C      $         cab056(l), cab057(l), cab058(l), cab059(l), cab060(l),  
C      $         cab061(l), cab062(l), cab063(l), cab064(l), cab065(l),  
C      $         cab066(l), cab067(l), cab068(l), cab069(l), cab070(l),  
C      $         cab071(l), cab072(l), cab073(l), cab074(l), cab075(l),  
C      $         cab076(l), cab077(l), cab078(l), cab079(l), cab080(l),  
C      $         cab081(l), cab082(l), cab083(l), cab084(l), cab085(l),  
C      $         cab086(l), cab087(l), cab088(l), cab089(l), cab090(l),  
C      $         cab091(l), cab092(l), cab093(l), cab094(l), cab095(l),  
C      $         cab096(l), cab097(l), cab098(l), cab099(l), cab100(l),   
C      $         cab101(l), cab102(l), cab103(l), cab104(l), cab105(l),  
C      $         cab106(l), cab107(l),
C      $         cl001(l), cl002(l), cl003(l),
C      $         cl004(l), cl005(l), cl006(l),
C      $         cl007(l), cl008(l), cl009(l),
C      $         cl010(l), cl011(l), cl012(l), 
C      $         cl013(l), cl014(l), cl015(l),
C      $         cl016(l), cl017(l), cl018(l),
C      $         cl019(l), cl020(l), cl021(l),
C      $         cl022(l), cl023(l), cl024(l),
C      $         cl025(l), cl026(l), cl027(l),
C      $         cl028(l), cl029(l), cl030(l),
C      $         cl031(l), cl032(l), cl033(l),
C      $         cl034(l), cl035(l), cl036(l),
C      $         cl037(l), cl038(l), cl039(l),
C      $         cl040(l), cl041(l), cl042(l),
C      $         cl043(l), cl044(l), cl045(l),
C      $         cl046(l), cl047(l), cl048(l),
C      $         cl049(l), cl050(l), cl051(l),
C      $         cl052(l), cl053(l), cl054(l), 
C      $         cl055(l), cl056(l), cl057(l),
C      $         cl058(l), cl059(l),
C      $         dccn_dpq, dcc0_dpq,
C      $         ro2(l), no(l), no2(l),
C      $         dHOX_dPQ, dHOX0_dPQ, k_pseudo)          

C           ENDIF 
  
      end do
                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     update number densities
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     long-lived species
c
      do l = 1,lswitch-1
                  
         cc(l,i_co2) = (cc0(l,i_co2) + production(l,i_co2)*dt)
     $              /(1. + loss(l,i_co2)*dt)
         cc(l,i_co)  = (cc0(l,i_co) + production(l,i_co)*dt)
     $              /(1. + loss(l,i_co)*dt)

         cc(l,i_o2)  = (cc0(l,i_o2) + production(l,i_o2)*dt)
     $              /(1. + loss(l,i_o2)*dt)
  
         cc(l,i_h2)  = (cc0(l,i_h2) + production(l,i_h2)*dt)
     $              /(1. + loss(l,i_h2)*dt)
         cc(l,i_h2o2)= (cc0(l,i_h2o2) + production(l,i_h2o2)*dt)
     $              /(1. + loss(l,i_h2o2)*dt)
         cc(l,i_h2o) = (cc0(l,i_h2o) + production(l,i_h2o)*dt)
     $              /(1. + loss(l,i_h2o)*dt)
    
         cc(l,i_hox) = (cc0(l,i_hox) + production(l,i_hox)*dt)
     $              /(1. + loss(l,i_hox)*dt)
     
C !      Methane Oxidation
         cc(l,i_ch4) = (cc0(l,i_ch4) + production(l,i_ch4)*dt)
     $              /(1. + loss(l,i_ch4)*dt)

         if ( igcm_ch3 .ne. 0 )
     $   cc(l,i_ch3) = production(l,i_ch3)/loss(l,i_ch3)
         
         if ( igcm_ch3o2 .ne. 0 )
     $   cc(l,i_ch3o2) = production(l,i_ch3o2)/loss(l,i_ch3o2)
         
         if ( igcm_ch3ooh .ne. 0 )
     $    cc(l,i_ch3ooh) = (cc0(l,i_ch3ooh) + production(l,i_ch3ooh)*dt)
     $                  /(1. + loss(l,i_ch3ooh)*dt)
     
         if ( igcm_ch3oh .ne. 0 )
     $    cc(l,i_ch3oh) = (cc0(l,i_ch3oh) + production(l,i_ch3oh)*dt)
     $              /(1. + loss(l,i_ch3oh)*dt)

         if ( igcm_ch3o .ne. 0 )
     $   cc(l,i_ch3o) =  production(l,i_ch3o)/loss(l,i_ch3o)
         
         if ( igcm_hcho .ne. 0 )
     $   cc(l,i_hcho) = (cc0(l,i_hcho) + production(l,i_hcho)*dt)
     $              /(1. + loss(l,i_hcho)*dt)

         if ( igcm_hco .ne. 0 )
     $   cc(l,i_hco) = production(l,i_hco)/loss(l,i_hco)   

         if ( igcm_hoch2o2 .ne. 0 )
     $    cc(l,i_hoch2o2) = production(l,i_hoch2o2)/loss(l,i_hoch2o2)
     
         if ( igcm_hcooh .ne. 0 )
     $    cc(l,i_hcooh) = (cc0(l,i_hcooh) + production(l,i_hcooh)*dt)
     $              /(1. + loss(l,i_hcooh)*dt)

         if ( igcm_hoch2oh .ne. 0 )
     $    cc(l,i_hoch2oh) = (cc0(l,i_hoch2oh) 
     $                   + production(l,i_hoch2oh)*dt)
     $              /(1. + loss(l,i_hoch2oh)*dt)
     
         if ( igcm_hoch2ooh .ne. 0 )
     $    cc(l,i_hoch2ooh) = (cc0(l,i_hoch2ooh) 
     $                    + production(l,i_hoch2ooh)*dt)
     $              /(1. + loss(l,i_hoch2ooh)*dt)
     

C !     Alkane Oxidation
         if ( igcm_c2h6 .ne. 0 ) then
         cc(l,i_c2h6) = (cc0(l,i_c2h6) + production(l,i_c2h6)*dt)
     $              /(1. + loss(l,i_c2h6)*dt)

         cc(l,i_c2h5) = production(l,i_c2h5)/loss(l,i_c2h5)   

         cc(l,i_c2h5ooh) = (cc0(l,i_c2h5ooh) 
     $                   + production(l,i_c2h5ooh)*dt)
     $              /(1. + loss(l,i_c2h5ooh)*dt)
     
         cc(l,i_c2h5oh) = (cc0(l,i_c2h5oh) + production(l,i_c2h5oh)*dt)
     $              /(1. + loss(l,i_c2h5oh)*dt)

         cc(l,i_c2h5o2) = (cc0(l,i_c2h5o2) + production(l,i_c2h5o2)*dt)
     $              /(1. + loss(l,i_c2h5o2)*dt)
     
         cc(l,i_hoch2ch2o)=production(l,i_hoch2ch2o)/loss(l,i_hoch2ch2o)
         
         cc(l,i_hoch2ch2o2)=production(l,i_hoch2ch2o2)/
     $                     loss(l,i_hoch2ch2o2)

         cc(l,i_ethgly) = (cc0(l,i_ethgly) 
     $                    + production(l,i_ethgly)*dt)
     $              /(1. + loss(l,i_ethgly)*dt)
     
         cc(l,i_hyetho2h) = (cc0(l,i_hyetho2h) 
     $                    + production(l,i_hyetho2h)*dt)
     $              /(1. + loss(l,i_hyetho2h)*dt)
     
         cc(l,i_ch3cho) = (cc0(l,i_ch3cho) + production(l,i_ch3cho)*dt)
     $              /(1. + loss(l,i_ch3cho)*dt)
     
        cc(l,i_ch2choh) =(cc0(l,i_ch2choh) + production(l,i_ch2choh)*dt)
     $              /(1. + loss(l,i_ch2choh)*dt)
     
         cc(l,i_ch3choho2)=production(l,i_ch3choho2)/
     $                     loss(l,i_ch3choho2)     
     
         cc(l,i_ch3cooh) = (cc0(l,i_ch3cooh) 
     $                   + production(l,i_ch3cooh)*dt)
     $              /(1. + loss(l,i_ch3cooh)*dt)
     
         cc(l,i_ch3chohooh) = (cc0(l,i_ch3chohooh) 
     $                      + production(l,i_ch3chohooh)*dt)
     $              /(1. + loss(l,i_ch3chohooh)*dt)
     
         cc(l,i_ch3co)=production(l,i_ch3co)/
     $                     loss(l,i_ch3co)
     
         cc(l,i_ch3cooo)=production(l,i_ch3cooo)/
     $                     loss(l,i_ch3cooo)
     
         cc(l,i_ch3coooh) = (cc0(l,i_ch3coooh) 
     $                    + production(l,i_ch3coooh)*dt)
     $              /(1. + loss(l,i_ch3coooh)*dt)
     
         cc(l,i_hcoch2o2)=production(l,i_hcoch2o2)/
     $                     loss(l,i_hcoch2o2)
     
         cc(l,i_glyox) = (cc0(l,i_glyox) + production(l,i_glyox)*dt)
     $              /(1. + loss(l,i_glyox)*dt)
     
         cc(l,i_hcoco)=production(l,i_hcoco)/
     $                     loss(l,i_hcoco)
     
         cc(l,i_hooch2cho) = (cc0(l,i_hooch2cho) 
     $                     + production(l,i_hooch2cho)*dt)
     $              /(1. + loss(l,i_hooch2cho)*dt)

         cc(l,i_hoch2cho) = (cc0(l,i_hoch2cho) 
     $                    + production(l,i_hoch2cho)*dt)
     $              /(1. + loss(l,i_hoch2cho)*dt)
     
         cc(l,i_hochcho)=production(l,i_hochcho)/
     $                     loss(l,i_hochcho)

         cc(l,i_hoch2co)=production(l,i_hoch2co)/
     $                     loss(l,i_hoch2co)

         cc(l,i_hoch2co3)=production(l,i_hoch2co3)/
     $                     loss(l,i_hoch2co3)
     
         cc(l,i_hoch2co2h) = (cc0(l,i_hoch2co2h) 
     $                     + production(l,i_hoch2co2h)*dt)
     $              /(1. + loss(l,i_hoch2co2h)*dt)
     
         cc(l,i_hcoco2h) = (cc0(l,i_hcoco2h) 
     $                   + production(l,i_hcoco2h)*dt)
     $              /(1. + loss(l,i_hcoco2h)*dt)
     
         cc(l,i_hoch2co3h) = (cc0(l,i_hoch2co3h) 
     $                     + production(l,i_hoch2co3h)*dt)
     $              /(1. + loss(l,i_hoch2co3h)*dt)
     
         cc(l,i_hcoco3h) = (cc0(l,i_hcoco3h) 
     $                   + production(l,i_hcoco3h)*dt)
     $              /(1. + loss(l,i_hcoco3h)*dt)
    
         cc(l,i_hcoco3)=production(l,i_hcoco3)/
     $                     loss(l,i_hcoco3)
     
      ENDIF 

      IF (igcm_ch3o2 .ne. 0) ro2(l) = ro2(l) + cc(l,i_ch3o2)
      IF (igcm_hoch2o2 .ne. 0 ) ro2(l) = ro2(l) + cc(l,i_hoch2o2)
      IF (igcm_c2h5o2 .ne. 0 ) ro2(l) = ro2(l) + cc(l,i_c2h5o2)
      IF (igcm_hoch2ch2o2 .ne. 0 ) ro2(l) = ro2(l) + cc(l,i_hoch2ch2o2)
      IF (igcm_ch3cooo .ne. 0 ) ro2(l) = ro2(l) + cc(l,i_ch3cooo)
      IF (igcm_hcoch2o2 .ne. 0 ) ro2(l) = ro2(l) + cc(l,i_hcoch2o2)

      IF ( ro2(l) .le. dens(l)*1.e-29 ) ro2(l) = 0.

ccccccccccccccccccccccccc
cccccc Chlorine Chemistry
ccccccccccccccccccccccccc
      IF ( igcm_cl .ne. 0 ) THEN 
        ! HCl
         cc(l,i_hcl) = Chlorine_QSSA(cc0(l,i_hcl),
     $                               production(l,i_hcl),
     $                               loss(l,i_hcl),
     $                               dt)

        ! HOCl
         cc(l,i_hocl) = Chlorine_QSSA(cc0(l,i_hocl),
     $                               production(l,i_hocl),
     $                               loss(l,i_hocl),
     $                               dt)

        ! ClO3 
         cc(l,i_clo3) = Chlorine_QSSA(cc0(l,i_clo3),
     $                               production(l,i_clo3),
     $                               loss(l,i_clo3),
     $                               dt)

        ! ClO4
         cc(l,i_clo4) = Chlorine_QSSA(cc0(l,i_clo4),
     $                               production(l,i_clo4),
     $                               loss(l,i_clo4),
     $                               dt)

        ! ClCO 
         cc(l,i_clco) = Chlorine_QSSA(cc0(l,i_clco),
     $                               production(l,i_clco),
     $                               loss(l,i_clco),
     $                               dt)

        ! Cl2
         cc(l,i_cl2) = Chlorine_QSSA(cc0(l,i_cl2),
     $                               production(l,i_cl2),
     $                               loss(l,i_cl2),
     $                               dt)

        ! Cl2O2 
         cc(l,i_cl2o2) = Chlorine_QSSA(cc0(l,i_cl2o2),
     $                               production(l,i_cl2o2),
     $                               loss(l,i_cl2o2),
     $                               dt)

        ! ClOO
         cc(l,i_cloo) = Chlorine_QSSA(cc0(l,i_cloo),
     $                               production(l,i_cloo),
     $                               loss(l,i_cloo),
     $                               dt)

        ! OClO
         cc(l,i_oclo) = Chlorine_QSSA(cc0(l,i_oclo),
     $                               production(l,i_oclo),
     $                               loss(l,i_oclo),
     $                               dt)

        ! CH3OCl
         cc(l,i_ch3ocl) = Chlorine_QSSA(cc0(l,i_ch3ocl),
     $                               production(l,i_ch3ocl),
     $                               loss(l,i_ch3ocl),
     $                               dt)

        ! HClO4
         cc(l,i_hclo4) = Chlorine_QSSA(cc0(l,i_hclo4),
     $                               production(l,i_hclo4),
     $                               loss(l,i_hclo4),
     $                               dt)
C         ! ClOX
C          cc(l,i_clox) = Chlorine_QSSA(cc0(l,i_clox),
C      $                               production(l,i_clox),
C      $                               loss(l,i_clox),
C      $                               dt)


C         cc(l,i_ch3ocl) = (cc0(l,i_ch3ocl) + production(l,i_ch3ocl)*dt)
C      $              /(1. + loss(l,i_ch3ocl)*dt)

C         cc(l,i_hclo4) = (cc0(l,i_hclo4) + production(l,i_hclo4)*dt)
C      $              /(1. + loss(l,i_hclo4)*dt)


        cc(l,i_clox) = (cc0(l,i_clox) + production(l,i_clox)*dt)
     $                /(1. + loss(l,i_clox)*dt)


ccccccccccccccccccccccccccc
c     Daylight Conditionals
ccccccccccccccccccccccccccc
C       IF (sza .le. 95.) THEN 
C         cc(l,i_oclo) = production(l,i_oclo)/loss(l,i_oclo)

C       ELSE
C         cc(l,i_oclo) = (cc0(l,i_oclo) + production(l,i_oclo)*dt)
C      $              /(1. + loss(l,i_oclo)*dt)

C       ENDIF 


      ENDIF 


      end do

       
c     ox species
      if (sza .le. 95.) then
         do l = 1,lswitch-1
            cc(l,i_ox) = (cc0(l,i_ox) + production(l,i_ox)*dt)
     $                   /(1. + loss(l,i_ox)*dt)         
         end do
      else
         do l = 1,lswitch-1
            cc(l,i_o)  = (cc0(l,i_o) + production(l,i_o)*dt)
     $                   /(1. + loss(l,i_o)*dt)
            cc(l,i_o3) = (cc0(l,i_o3) + production(l,i_o3)*dt)
     $                   /(1. + loss(l,i_o3)*dt)
            cc(l,i_ox) = cc(l,i_o) + cc(l,i_o3)          
         end do
      end if 







cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end of loop over iterations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      end do
  

      ! ================
      ! Loss Rate of O2:
      ! ================
      o2_loss_rate =  j(:,j_o2_o)*cc(:,i_o2)
     &             +  j(:,j_o2_o1d)*cc(:,i_o2)
     &             +  j(:,j_ho2)*cc(:,i_ho2)
     &             +  j(:,j_h2o2)*cc(:,i_h2o2) 
     &             +  j(:,j_ch3o2h)*cc(:,i_ch3ooh)     
     &             +  c004(:)*cc(:,i_h)*cc(:,i_ho2)
     &             +  c006(:)*cc(:,i_h)*cc(:,i_ho2)
     &             +  d003(:)*cc(:,i_ho2)*no(:)
     &             +  cab007*cc(:,i_ch3o2)*cc(:,i_ho2) 
     &             +  cab008*cc(:,i_ch3o2)*ro2(:)*1.5
     &             +  cab009*cc(:,i_ch3o2)*ro2(:)*1.5
     &             +  cab014*cc(:,i_ch3ooh)*cc(:,i_oh)*0.4 
     &             +  0.5*cab029*cc(:,i_hoch2o2)*cc(:,i_ho2)
     &             +  cab034*cc(:,i_hoch2ooh)*cc(:,i_hoch2ooh)
     $             + cl010*cc(:,i_cl)*cc(:,i_ho2)
     $             + cl018*cc(:,i_cl)*cc(:,i_ch3ooh)
     $             + cl021*cc(:,i_clo)*cc(:,i_ch3o2)
     $             + cl022*cc(:,i_cl)*Cc(:,i_ch3o2)*0.5
     $             + cl024*cc(:,i_cloo)*cc(:,i_cl)
     $             + cl046*cc(:,i_clo)*cc(:,i_clo3)
     $             + cl047*cc(:,i_clo)*cc(:,i_clo3)
     $             + cl048*cc(:,i_clo3)*cc(:,i_oh)
     $             + cl049*cc(:,i_clo3)*cc(:,i_oh)
     $             + cl054*cc(:,i_clo4)*cc(:,i_hocl)
     $             + j(:,j_cloo)*cc(:,i_cloo) 

      ! ---------------------------------------------
      ! Atmospheric Lifetimes [Inverse loss rate = s]
      ! ---------------------------------------------
      lifetime = 1./loss

      ! -----------------------------
      ! HCl Production and Loss Rates
      ! -----------------------------
      lifetime(:,1) = production(:,i_hcl)
      lifetime(:,2) = cl008*cc(:,i_cl)*cc(:,i_h2)
      lifetime(:,3) = cl009(:)*cc(:,i_cl)*cc(:,i_ho2)
      lifetime(:,4) = cl011(:)*cc(:,i_cl)*cc(:,i_h2o2)
      lifetime(:,5) = cl012(:)*cc(:,i_clo)*cc(:,i_oh)*0.06
      lifetime(:,6) = cl016(:)*cc(:,i_cl)*cc(:,i_ch4)
      lifetime(:,7) = cl017(:)*cc(:,i_cl)*cc(:,i_hcho) 
      lifetime(:,8) = cl018(:)*cc(:,i_cl)*cc(:,i_ch3ooh)
      lifetime(:,9) = cl022(:)*cc(:,i_cl)*cc(:,i_ch3o2)*0.5
      lifetime(:,10) = cl031(:)*cc(:,i_cl)*cc(:,i_ch3ocl)
      lifetime(:,11) = cl037(:)*cc(:,i_cl2)*cc(:,i_h)
      lifetime(:,12) = cl043(:)*cc(:,i_cl)*cc(:,i_hcooh)
      lifetime(:,13) = cl056(:)*cc(:,i_cl)*Cc(:,i_hocl)

      lifetime(:,14) = j(:,j_hcl) 
      lifetime(:,15) = cl014(:)*cc(:,i_oh)
      lifetime(:,16) = cl039(:)*cc(:,i_o1d)*0.88
      lifetime(:,17) = cl040(:)*cc(:,i_o) 
      lifetime(:,18) = cl041(:)*cc(:,i_h) 

      ! Primary Organic Lifetimes
      lifetime(:,19) = 1./loss(:,i_ch4)
      lifetime(:,20) = 1./loss(:,i_hcho)
      lifetime(:,21) = 1./loss(:,i_hcooh)
      lifetime(:,22) = 1./loss(:,i_ch3ooh)

      ! HClO4 Production Rates
      lifetime(:,23) = production(:,i_hclo4) 
      lifetime(:,24) = cl048(:)*cc(:,i_clo3)*cc(:,i_oh)
      lifetime(:,25) = cl049(:)*cc(:,i_clo3)*cc(:,i_oh)
      lifetime(:,26) = cl054(:)*cc(:,i_clo4)*cc(:,i_hocl)

      ! ClOx
      lifetime(:,27) = production(:,i_clox)
      lifetime(:,28) = loss(:,i_clox)
      lifetime(:,29) = pcl
      lifetime(:,30) = lcl 
      lifetime(:,31) = pclo
      lifetime(:,32) = lclo 

      ! Cl2 
      lifetime(:,33) =  2.*cl003(:)*cc(:,i_clo) 
      lifetime(:,34) =  cl005(:)*cc(:,i_clo)
      lifetime(:,35) =  2.*cl006(:)*cc(:,i_clo)
      lifetime(:,36) =  cl008(:)*cc(:,i_h2)
      lifetime(:,37) =  cl009(:)*cc(:,i_ho2)
      lifetime(:,38) =  cl011(:)*cc(:,i_h2o2)
      lifetime(:,39) =  0.06*cl012(:)*cc(:,i_oh)
      lifetime(:,40) =  cl013(:)*cc(:,i_ho2)
      lifetime(:,41) =  cl016(:)*cc(:,i_ch4)
      lifetime(:,42) =  cl017(:)*cc(:,i_hcho)
      lifetime(:,43) =  cl018(:)*cc(:,i_ch3ooh)
      lifetime(:,44) =  cl019(:)*cc(:,i_ch3o2)
      lifetime(:,45) = cl020(:)*cc(:,i_ch3o2)
      lifetime(:,46) = cl021(:)*cc(:,i_ch3o2)
      lifetime(:,47) = 0.5*cl022(:)*cc(:,i_ch3o2)
      lifetime(:,48) = cl023(:)*Cc(:,i_co)
      lifetime(:,49) = cl025(:)*cc(:,i_cloo)
      lifetime(:,50) = cl028(:)*Cc(:,i_o2)
      lifetime(:,51) = cl029(:)*cc(:,i_ch3ocl)
      lifetime(:,52) =  cl030(:)*cc(:,i_cl2o2)
      lifetime(:,53) =  cl031(:)*cc(:,i_ch3ocl)
      lifetime(:,54) =  cl043(:)*cc(:,i_hcooh)
      lifetime(:,55) =  cl044(:)*cc(:,i_o3)
      lifetime(:,56) =  cl045(:)*cc(:,i_clo3)
      lifetime(:,57) =  cl046(:)*cc(:,i_clo3)
      lifetime(:,58) =  cl047(:)*cc(:,i_clo3)
      lifetime(:,59) =  cl055(:)*cc(:,i_hocl)

      ! O=O Loss Rate
      lifetime(:,60) = o2_loss_rate


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     density -> volume mixing ratio conversion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iesp = 1,nesp
         if ( iesp == i_dust ) cycle 
         do l = 1,lswitch-1
            rm(l,iesp) = max(cc(l,iesp)/dens(l), 1.e-30)
         end do
      end do

c
      return
      end
c
c*****************************************************************
c
      subroutine phot(lswitch, press, temp, sza, tauref, dist_sol, 
     $                rmco2, rmo3, j)
c
c*****************************************************************
c
      implicit none
c
#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h" 
#include "comcstfi.h"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
      integer lswitch         ! interface level between chemistries
      real press(nlayermx)             ! pressure (hPa)
      real temp(nlayermx)              ! temperature (K)
      real sza                         ! solar zenith angle (deg)
      real tauref                      ! optical depth at 7 hpa
      real dist_sol                    ! sun distance (AU) 
      real rmco2(nlayermx)             ! co2 mixing ratio
      real rmo3(nlayermx)              ! ozone mixing ratio
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     output:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real j(nlayermx,nd)              ! interpolated photolysis rates (s-1) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer icol, ij, indsza, indtau, indcol, indozo, indtemp, 
     $        iozo, isza, itau, it, l
c
      real col(nlayermx)               ! overhead air column   (molecule cm-2)
      real colo3(nlayermx)             ! overhead ozone column (molecule cm-2)
      real poids(2,2,2,2,2)            ! 5D interpolation weights 
      real tref                        ! temperature  at 1.9 hPa in the gcm (K)
      real table_temp(ntemp)           ! temperatures at 1.9 hPa in jmars   (K)
      real cinf, csup, cicol, 
     $     ciozo, cisza, citemp, citau
      real colo3min, dp, coef
      real ratio_o3(nlayermx)
      real tau
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     day/night criterion
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (sza .le. 95.) then 

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     temperatures at 1.9 hPa in lookup table
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      table_temp(1) = 226.2
      table_temp(2) = 206.2
      table_temp(3) = 186.2
      table_temp(4) = 169.8
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     interpolation in solar zenith angle
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      indsza = nsza - 1
      do isza = 1,nsza
         if (szatab(isza) .ge. sza) then
            indsza = min(indsza,isza - 1)
            indsza = max(indsza, 1)
         end if
      end do
      cisza = (sza - szatab(indsza))
     $       /(szatab(indsza + 1) - szatab(indsza))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     interpolation in dust (tau)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      tau = min(tauref, tautab(ntau))
      tau = max(tau, tautab(1))
c
      indtau = ntau - 1
      do itau = 1,ntau
         if (tautab(itau) .ge. tau) then
            indtau = min(indtau,itau - 1)
            indtau = max(indtau, 1)
         end if
      end do
      citau = (tau - tautab(indtau))
     $       /(tautab(indtau + 1) - tautab(indtau))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     co2 and ozone columns
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     co2 column at model top (molecule.cm-2)
c
      col(lswitch-1) = 6.022e22*rmco2(lswitch-1)*press(lswitch-1)*100.
     $                 /(mugaz*g)
c
c     ozone column at model top
c
      colo3(lswitch-1) = 0.
c
c     co2 and ozone columns for other levels (molecule.cm-2)
c
      do l = lswitch-2,1,-1
         dp = (press(l) - press(l+1))*100.
         col(l) = col(l+1) 
     $            + (rmco2(l+1) + rmco2(l))*0.5
     $              *6.022e22*dp/(mugaz*g)
         col(l) = min(col(l), colairtab(0))
         colo3(l) = colo3(l+1)
     $            + (rmo3(l+1) + rmo3(l))*0.5
     $              *6.022e22*dp/(mugaz*g)
      end do
c
c     ratio ozone column/minimal theoretical column (0.1 micron-atm) 
c
c     ro3 = 7.171e-10 is the o3 mixing ratio for a uniform
c     profile giving a column 0.1 micron-atmosphere at
c     a surface pressure of 10 hpa.
c
      do l = 1,lswitch-1
         colo3min    = col(l)*7.171e-10
         ratio_o3(l) = colo3(l)/colo3min
         ratio_o3(l) = min(ratio_o3(l), table_ozo(nozo)*10.)
         ratio_o3(l) = max(ratio_o3(l), 1.)
      end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     temperature dependence
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     1) search for temperature at 1.9 hPa (tref): vertical interpolation
c
      tref = temp(1)
      do l = (lswitch-1)-1,1,-1
         if (press(l) .gt. 1.9) then
            cinf = (press(l) - 1.9)
     $            /(press(l) - press(l+1))
            csup = 1. - cinf
            tref = cinf*temp(l+1) + csup*temp(l)
            go to 10
         end if
      end do 
 10   continue
c
c     2) interpolation in temperature
c
      tref = min(tref,table_temp(1))
      tref = max(tref,table_temp(ntemp))
c
      do it = 2, ntemp
         if (table_temp(it) .le. tref) then
            citemp = (log(tref) - log(table_temp(it)))
     $              /(log(table_temp(it-1)) - log(table_temp(it)))
            indtemp = it - 1
            goto 20
         end if
      end do
  20  continue
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over vertical levels
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
c
c     interpolation in air column
c
         do icol = 0,200
            if (colairtab(icol) .lt. col(l)) then
               cicol = (log(col(l)) - log(colairtab(icol)))
     $                /(log(colairtab(icol-1)) - log(colairtab(icol)))
               indcol = icol - 1
               goto 30
            end if
         end do
 30      continue
c
cc    interpolation in ozone column
c
         indozo = nozo - 1
         do iozo = 1,nozo
            if (table_ozo(iozo)*10. .ge. ratio_o3(l)) then
               indozo = min(indozo, iozo - 1)
               indozo = max(indozo, 1)
            end if
         end do
         ciozo = (ratio_o3(l) - table_ozo(indozo)*10.)
     $          /(table_ozo(indozo + 1)*10. - table_ozo(indozo)*10.)
c



cc    4-dimensional interpolation weights
c
c     poids(temp,sza,co2,o3,tau)
c
         poids(1,1,1,1,1) = citemp
     $                     *(1.-cisza)*cicol*(1.-ciozo)*(1.-citau)
         poids(1,1,1,2,1) = citemp
     $                     *(1.-cisza)*cicol*ciozo*(1.-citau)
         poids(1,1,2,1,1) = citemp
     $                     *(1.-cisza)*(1.-cicol)*(1.-ciozo)*(1.-citau)
         poids(1,1,2,2,1) = citemp
     $                     *(1.-cisza)*(1.-cicol)*ciozo*(1.-citau)
         poids(1,2,1,1,1) = citemp
     $                     *cisza*cicol*(1.-ciozo)*(1.-citau)
         poids(1,2,1,2,1) = citemp
     $                     *cisza*cicol*ciozo*(1.-citau)
         poids(1,2,2,1,1) = citemp
     $                     *cisza*(1.-cicol)*(1.-ciozo)*(1.-citau)
         poids(1,2,2,2,1) = citemp
     $                     *cisza*(1.-cicol)*ciozo*(1.-citau)
         poids(2,1,1,1,1) = (1.-citemp)
     $                     *(1.-cisza)*cicol*(1.-ciozo)*(1.-citau)
         poids(2,1,1,2,1) = (1.-citemp)
     $                     *(1.-cisza)*cicol*ciozo*(1.-citau)
         poids(2,1,2,1,1) = (1.-citemp)
     $                     *(1.-cisza)*(1.-cicol)*(1.-ciozo)*(1.-citau)
         poids(2,1,2,2,1) = (1.-citemp)
     $                     *(1.-cisza)*(1.-cicol)*ciozo*(1.-citau)
         poids(2,2,1,1,1) = (1.-citemp)
     $                     *cisza*cicol*(1.-ciozo)*(1.-citau)
         poids(2,2,1,2,1) = (1.-citemp)
     $                     *cisza*cicol*ciozo*(1.-citau)
         poids(2,2,2,1,1) = (1.-citemp)
     $                     *cisza*(1.-cicol)*(1.-ciozo)*(1.-citau)
         poids(2,2,2,2,1) = (1.-citemp)
     $                     *cisza*(1.-cicol)*ciozo*(1.-citau)
c
         poids(1,1,1,1,2) = citemp
     $                     *(1.-cisza)*cicol*(1.-ciozo)*citau
         poids(1,1,1,2,2) = citemp
     $                     *(1.-cisza)*cicol*ciozo*citau
         poids(1,1,2,1,2) = citemp
     $                     *(1.-cisza)*(1.-cicol)*(1.-ciozo)*citau
         poids(1,1,2,2,2) = citemp
     $                     *(1.-cisza)*(1.-cicol)*ciozo*citau
         poids(1,2,1,1,2) = citemp
     $                     *cisza*cicol*(1.-ciozo)*citau
         poids(1,2,1,2,2) = citemp
     $                     *cisza*cicol*ciozo*citau
         poids(1,2,2,1,2) = citemp
     $                     *cisza*(1.-cicol)*(1.-ciozo)*citau
         poids(1,2,2,2,2) = citemp
     $                     *cisza*(1.-cicol)*ciozo*citau
         poids(2,1,1,1,2) = (1.-citemp)
     $                     *(1.-cisza)*cicol*(1.-ciozo)*citau
         poids(2,1,1,2,2) = (1.-citemp)
     $                     *(1.-cisza)*cicol*ciozo*citau
         poids(2,1,2,1,2) = (1.-citemp)
     $                     *(1.-cisza)*(1.-cicol)*(1.-ciozo)*citau
         poids(2,1,2,2,2) = (1.-citemp)
     $                     *(1.-cisza)*(1.-cicol)*ciozo*citau
         poids(2,2,1,1,2) = (1.-citemp)
     $                     *cisza*cicol*(1.-ciozo)*citau
         poids(2,2,1,2,2) = (1.-citemp)
     $                     *cisza*cicol*ciozo*citau
         poids(2,2,2,1,2) = (1.-citemp)
     $                     *cisza*(1.-cicol)*(1.-ciozo)*citau
         poids(2,2,2,2,2) = (1.-citemp)
     $                     *cisza*(1.-cicol)*ciozo*citau
c
cc    4-dimensional interpolation in the lookup table
c
         do ij = 1,nd
           j(l,ij) =
     $     poids(1,1,1,1,1)
     $     *jphot(indtemp,indsza,indcol,indozo,indtau,ij)
     $   + poids(1,1,1,2,1)
     $     *jphot(indtemp,indsza,indcol,indozo+1,indtau,ij)
     $   + poids(1,1,2,1,1)
     $     *jphot(indtemp,indsza,indcol+1,indozo,indtau,ij)
     $   + poids(1,1,2,2,1)
     $     *jphot(indtemp,indsza,indcol+1,indozo+1,indtau,ij)
     $   + poids(1,2,1,1,1)
     $     *jphot(indtemp,indsza+1,indcol,indozo,indtau,ij)
     $   + poids(1,2,1,2,1)
     $     *jphot(indtemp,indsza+1,indcol,indozo+1,indtau,ij)
     $   + poids(1,2,2,1,1)
     $     *jphot(indtemp,indsza+1,indcol+1,indozo,indtau,ij)
     $   + poids(1,2,2,2,1)
     $     *jphot(indtemp,indsza+1,indcol+1,indozo+1,indtau,ij)
     $   + poids(2,1,1,1,1)
     $     *jphot(indtemp+1,indsza,indcol,indozo,indtau,ij)
     $   + poids(2,1,1,2,1)
     $     *jphot(indtemp+1,indsza,indcol,indozo+1,indtau,ij)
     $   + poids(2,1,2,1,1)
     $     *jphot(indtemp+1,indsza,indcol+1,indozo,indtau,ij)
     $   + poids(2,1,2,2,1)
     $     *jphot(indtemp+1,indsza,indcol+1,indozo+1,indtau,ij)
     $   + poids(2,2,1,1,1)
     $     *jphot(indtemp+1,indsza+1,indcol,indozo,indtau,ij)
     $   + poids(2,2,1,2,1)
     $     *jphot(indtemp+1,indsza+1,indcol,indozo+1,indtau,ij)
     $   + poids(2,2,2,1,1)
     $     *jphot(indtemp+1,indsza+1,indcol+1,indozo,indtau,ij)
     $   + poids(2,2,2,2,1)
     $     *jphot(indtemp+1,indsza+1,indcol+1,indozo+1,indtau,ij)
c
     $   + poids(1,1,1,1,2)
     $     *jphot(indtemp,indsza,indcol,indozo,indtau+1,ij)
     $   + poids(1,1,1,2,2)
     $     *jphot(indtemp,indsza,indcol,indozo+1,indtau+1,ij)
     $   + poids(1,1,2,1,2)
     $     *jphot(indtemp,indsza,indcol+1,indozo,indtau+1,ij)
     $   + poids(1,1,2,2,2)
     $     *jphot(indtemp,indsza,indcol+1,indozo+1,indtau+1,ij)
     $   + poids(1,2,1,1,2)
     $     *jphot(indtemp,indsza+1,indcol,indozo,indtau+1,ij)
     $   + poids(1,2,1,2,2)
     $     *jphot(indtemp,indsza+1,indcol,indozo+1,indtau+1,ij)
     $   + poids(1,2,2,1,2)
     $     *jphot(indtemp,indsza+1,indcol+1,indozo,indtau+1,ij)
     $   + poids(1,2,2,2,2)
     $     *jphot(indtemp,indsza+1,indcol+1,indozo+1,indtau+1,ij)
     $   + poids(2,1,1,1,2)
     $     *jphot(indtemp+1,indsza,indcol,indozo,indtau+1,ij)
     $   + poids(2,1,1,2,2)
     $     *jphot(indtemp+1,indsza,indcol,indozo+1,indtau+1,ij)
     $   + poids(2,1,2,1,2)
     $     *jphot(indtemp+1,indsza,indcol+1,indozo,indtau+1,ij)
     $   + poids(2,1,2,2,2)
     $     *jphot(indtemp+1,indsza,indcol+1,indozo+1,indtau+1,ij)
     $   + poids(2,2,1,1,2)
     $     *jphot(indtemp+1,indsza+1,indcol,indozo,indtau+1,ij)
     $   + poids(2,2,1,2,2)
     $     *jphot(indtemp+1,indsza+1,indcol,indozo+1,indtau+1,ij)
     $   + poids(2,2,2,1,2)
     $     *jphot(indtemp+1,indsza+1,indcol+1,indozo,indtau+1,ij)
     $   + poids(2,2,2,2,2)
     $     *jphot(indtemp+1,indsza+1,indcol+1,indozo+1,indtau+1,ij)
         end do
		 

c
cc    correction for sun distance
c
         do ij = 1,nd
            j(l,ij) = j(l,ij)*(1.52/dist_sol)**2.
         end do 
              
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end of loop over vertical levels
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      end do
c

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     night
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      else

      do ij = 1,nd
         do l = 1,lswitch-1
            j(l,ij) = 0.
         end do
      end do
c
      end if
c
      return
      end
c
c*****************************************************************
c
      subroutine gcmtochim(zycol, lswitch, nesp, rm)
c
c*****************************************************************
c
      implicit none
c
#include "dimensions.h"
#include "dimphys.h"
#include "callkeys.h"
#include "tracer.h"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      real zycol(nlayermx,nqmx)! species volume mixing ratio in the gcm
c
      integer nesp             ! number of species in the chemistry
      integer lswitch          ! interface level between chemistries
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     outputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real rm(nlayermx,nesp)   ! species volume mixing ratio 
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer l,iq
      
c     tracer indexes in the chemistry:

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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

       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialise chemical species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
         rm(l,i_co2)  = max(zycol(l, igcm_co2),     1.e-30)
         rm(l,i_co)   = max(zycol(l, igcm_co),      1.e-30)
         rm(l,i_o)    = max(zycol(l, igcm_o),       1.e-30)
         rm(l,i_o1d)  = max(zycol(l, igcm_o1d),     1.e-30)
         rm(l,i_o2)   = max(zycol(l, igcm_o2),      1.e-30)
         rm(l,i_o3)   = max(zycol(l, igcm_o3),      1.e-30)
         rm(l,i_h)    = max(zycol(l, igcm_h),       1.e-30)
         rm(l,i_h2)   = max(zycol(l, igcm_h2),      1.e-30)
         rm(l,i_oh)   = max(zycol(l, igcm_oh),      1.e-30)
         rm(l,i_ho2)  = max(zycol(l, igcm_ho2),     1.e-30)
         rm(l,i_h2o2) = max(zycol(l, igcm_h2o2),    1.e-30)
         rm(l,i_n2)   = max(zycol(l, igcm_n2),      1.e-30)
         rm(l,i_h2o)  = max(zycol(l, igcm_h2o_vap), 1.e-30)
      end do 

      if (igcm_ch4 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch4) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch4) = max(zycol(l,igcm_ch4), 1.e-30)
         end do
      end if
      if (igcm_ch3 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3) = max(zycol(l,igcm_ch3), 1.e-30)
         end do
      end if
      if (igcm_ch3o2 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3o2) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3o2) = max(zycol(l,igcm_ch3o2), 1.e-30)
         end do
      end if
      if (igcm_ch3ooh .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3ooh) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3ooh) = max(zycol(l,igcm_ch3ooh), 1.e-30)
         end do
      end if
      if (igcm_ch3oh .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3oh) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3oh) = max(zycol(l,igcm_ch3oh), 1.e-30)
         end do
      end if
      if (igcm_ch3o .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3o) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3o) = max(zycol(l,igcm_ch3o), 1.e-30)
         end do
      end if
      if (igcm_hcho .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hcho) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hcho) = max(zycol(l,igcm_hcho), 1.e-30)
         end do
      end if
      if (igcm_hcooh .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hcooh) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hcooh) = max(zycol(l,igcm_hcooh), 1.e-30)
         end do
      end if
      if (igcm_hoch2o2 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hoch2o2) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hoch2o2) = max(zycol(l,igcm_hoch2o2), 1.e-30)
         end do
      end if
      if (igcm_hoch2oh .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hoch2oh) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hoch2oh) = max(zycol(l,igcm_hoch2oh), 1.e-30)
         end do
      end if
      if (igcm_hoch2ooh .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hoch2ooh) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hoch2ooh) = max(zycol(l,igcm_hoch2ooh), 1.e-30)
         end do
      end if
      if (igcm_hco .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hco) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hco) = max(zycol(l,igcm_hco), 1.e-30)
         end do
      end if
      if (igcm_c2h6 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_c2h6) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_c2h6) = max(zycol(l,igcm_c2h6), 1.e-30)
         end do
      end if
      if (igcm_c2h5 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_c2h5) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_c2h5) = max(zycol(l,igcm_c2h5), 1.e-30)
         end do
      end if
      if (igcm_c2h5o2 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_c2h5o2) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_c2h5o2) = max(zycol(l,igcm_c2h5o2), 1.e-30)
         end do
      end if
      if (igcm_c2h5ooh .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_c2h5ooh) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_c2h5ooh) = max(zycol(l,igcm_c2h5ooh), 1.e-30)
         end do
      end if
      if (igcm_c2h5oh .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_c2h5oh) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_c2h5oh) = max(zycol(l,igcm_c2h5oh), 1.e-30)
         end do
      end if
      if (igcm_hoch2ch2o2 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hoch2ch2o2) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hoch2ch2o2) = max(zycol(l,igcm_hoch2ch2o2), 1.e-30)
         end do
      end if
      if (igcm_hoch2ch2o .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hoch2ch2o) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hoch2ch2o) = max(zycol(l,igcm_hoch2ch2o), 1.e-30)
         end do
      end if
      if (igcm_ethgly .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ethgly) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ethgly) = max(zycol(l,igcm_ethgly), 1.e-30)
         end do
      end if
      if (igcm_hyetho2h .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hyetho2h) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hyetho2h) = max(zycol(l,igcm_hyetho2h), 1.e-30)
         end do
      end if
      if (igcm_ch3cho .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3cho) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3cho) = max(zycol(l,igcm_ch3cho), 1.e-30)
         end do
      end if
      if (igcm_ch2choh .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch2choh) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch2choh) = max(zycol(l,igcm_ch2choh), 1.e-30)
         end do
      end if  
      if (igcm_ch3choho2 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3choho2) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3choho2) = max(zycol(l,igcm_ch3choho2), 1.e-30)
         end do
      end if
      if (igcm_ch3cooh .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3cooh) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3cooh) = max(zycol(l,igcm_ch3cooh), 1.e-30)
         end do
      end if
      if (igcm_ch3chohooh .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3chohooh) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3chohooh) = max(zycol(l,igcm_ch3chohooh), 1.e-30)
         end do
      end if
      if (igcm_ch3co .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3co) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3co) = max(zycol(l,igcm_ch3co), 1.e-30)
         end do
      end if
      if (igcm_ch3cooo .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3cooo) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3cooo) = max(zycol(l,igcm_ch3cooo), 1.e-30)
         end do
      end if
      if (igcm_ch3coooh .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3coooh) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3coooh) = max(zycol(l,igcm_ch3coooh), 1.e-30)
         end do
      end if
      if (igcm_hcoch2o2 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hcoch2o2) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hcoch2o2) = max(zycol(l,igcm_hcoch2o2), 1.e-30)
         end do
      end if
      if (igcm_glyox .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_glyox) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_glyox) = max(zycol(l,igcm_glyox), 1.e-30)
         end do
      end if
      if (igcm_hcoco .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hcoco) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hcoco) = max(zycol(l,igcm_hcoco), 1.e-30)
         end do
      end if
      if (igcm_hooch2cho .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hooch2cho) = 0.
         end do
      else
         do l = 1,lswitch-1      
            rm(l,i_hooch2cho) = max(zycol(l,igcm_hooch2cho), 1.e-30)
         end do
      end if
      if (igcm_hoch2cho .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hoch2cho) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hoch2cho) = max(zycol(l,igcm_hoch2cho), 1.e-30)
         end do
      end if
      if (igcm_hochcho .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hochcho) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hochcho) = max(zycol(l,igcm_hochcho), 1.e-30)
         end do
      end if
      if (igcm_hoch2co .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hoch2co) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hoch2co) = max(zycol(l,igcm_hoch2co), 1.e-30)
         end do
      end if
      if (igcm_hoch2co3 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hoch2co3) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hoch2co3) = max(zycol(l,igcm_hoch2co3), 1.e-30)
         end do
      end if
      if (igcm_hoch2co2h .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hoch2co2h) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hoch2co2h) = max(zycol(l,igcm_hoch2co2h), 1.e-30)
         end do
      end if
      if (igcm_hcoco2h .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hcoco2h) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hcoco2h) = max(zycol(l,igcm_hcoco2h), 1.e-30)
         end do
      end if
      if (igcm_hcoco3h .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hcoco3h) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hcoco3h) = max(zycol(l,igcm_hcoco3h), 1.e-30)
         end do
      end if      
      if (igcm_hcoco3 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hcoco3) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hcoco3) = max(zycol(l,igcm_hcoco3), 1.e-30)
         end do
      end if
      if (igcm_hoch2co3h .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hoch2co3h) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hoch2co3h) = max(zycol(l,igcm_hoch2co3h), 1.e-30)
         end do
      end if

      if (igcm_dust_mass .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_dust) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_dust) = max(zycol(l,igcm_dust_mass), 1.e-30)
         end do
      end if

!     Chlorine Chemistry

      if (igcm_cl .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_cl) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_cl) = max(zycol(l,igcm_cl), 0.)
         end do
      end if
      if (igcm_clo .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_clo) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_clo) = max(zycol(l,igcm_clo), 0.)
         end do
      end if
      if (igcm_cl2 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_cl2) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_cl2) = max(zycol(l,igcm_cl2), 0.)
         end do
      end if
      if (igcm_oclo .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_oclo) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_oclo) = max(zycol(l,igcm_oclo), 0.)
         end do
      end if
      if (igcm_cl2o2 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_cl2o2) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_cl2o2) = max(zycol(l,igcm_cl2o2), 0.)
         end do
      end if
      if (igcm_hcl .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hcl) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hcl) = max(zycol(l,igcm_hcl), 0.)
         end do
      end if
      if (igcm_hocl .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hocl) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hocl) = max(zycol(l,igcm_hocl),0.)
         end do
      end if
      if (igcm_cloo .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_cloo) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_cloo) = max(zycol(l,igcm_cloo), 0.)
         end do
      end if
      if (igcm_ch3ocl .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_ch3ocl) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_ch3ocl) = max(zycol(l,igcm_ch3ocl), 0.)
         end do
      end if
      if (igcm_clco .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_clco) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_clco) = max(zycol(l,igcm_clco), 0.)
         end do
      end if
      if (igcm_clo3 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_clo3) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_clo3) = max(zycol(l,igcm_clo3), 0.)
         end do
      end if      
      if (igcm_hclo4 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_hclo4) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_hclo4) = max(zycol(l,igcm_hclo4), 0.)
         end do
      end if
      if (igcm_clo4 .eq. 0) then
         do l = 1,lswitch-1
            rm(l,i_clo4) = 0.
         end do
      else
         do l = 1,lswitch-1
            rm(l,i_clo4) = max(zycol(l,igcm_clo4), 0.)
         end do
      end if


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialise chemical families                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
         rm(l,i_hox) = rm(l,i_h) 
     $               + rm(l,i_oh)
     $               + rm(l,i_ho2)
         rm(l,i_ox)  = rm(l,i_o) 
     $               + rm(l,i_o3)

         rm(l,i_clox) = rm(l,i_cl) + rm(l,i_clo)

c        Suppressing the production of Cl compounds
c        when no chlorine is present
         IF ( rm(l,i_clox) .le. 1.e-28 ) THEN 
            rm(l,i_clox) = 0.! 
         ENDIF 

      end do
c
      return
      end
c
c*****************************************************************
c
      subroutine chimtogcm(zycol, lswitch, nesp, rm)
c
c*****************************************************************
c
      implicit none
c
#include "dimensions.h"
#include "dimphys.h"
#include "callkeys.h"
#include "tracer.h"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer nesp               ! number of species in the chemistry
      integer lswitch            ! interface level between chemistries
c      
      real rm(nlayermx,nesp)     ! species volume mixing ratio 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     output:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      real zycol(nlayermx,nqmx)  ! species volume mixing ratio in the gcm
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer l,iq
      
c     tracer indexes in the chemistry:

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     save mixing ratios for the gcm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
         zycol(l, igcm_co2)     = rm(l,i_co2) 
         zycol(l, igcm_co)      = rm(l,i_co) 
         zycol(l, igcm_o)       = rm(l,i_o) 
         zycol(l, igcm_o1d)     = rm(l,i_o1d)
         zycol(l, igcm_o2)      = rm(l,i_o2) 
         zycol(l, igcm_o3)      = rm(l,i_o3) 
         zycol(l, igcm_h)       = rm(l,i_h)  
         zycol(l, igcm_h2)      = rm(l,i_h2) 
         zycol(l, igcm_oh)      = rm(l,i_oh) 
         zycol(l, igcm_ho2)     = rm(l,i_ho2) 
         zycol(l, igcm_h2o2)    = rm(l,i_h2o2)
         zycol(l, igcm_n2)      = rm(l,i_n2)
         zycol(l, igcm_h2o_vap) = rm(l,i_h2o)
      end do 

      if (igcm_ch4 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch4) = rm(l,i_ch4)
         end do
      end if
      if (igcm_ch3 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3) = rm(l,i_ch3)
         end do
      end if
      if (igcm_ch3o2 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3o2) = rm(l,i_ch3o2)
         end do
      end if
      if (igcm_ch3ooh .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3ooh) = rm(l,i_ch3ooh)
         end do
      end if
      if (igcm_ch3oh .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3oh) = rm(l,i_ch3oh)
         end do
      end if
      if (igcm_ch3o .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3o) = rm(l,i_ch3o)
         end do
      end if
      if (igcm_hcho .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hcho) = rm(l,i_hcho)
         end do
      end if
      if (igcm_hcooh .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hcooh) = rm(l,i_hcooh)
         end do
      end if
      if (igcm_hoch2o2 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hoch2o2) = rm(l,i_hoch2o2)
         end do
      end if
      if (igcm_hoch2oh .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hoch2oh) = rm(l,i_hoch2oh)
         end do
      end if
      if (igcm_hoch2ooh .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hoch2ooh) = rm(l,i_hoch2ooh)
         end do
      end if
      if (igcm_hco .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hco) = rm(l,i_hco)
         end do
      end if
      if (igcm_c2h6 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_c2h6) = rm(l,i_c2h6)
         end do
      end if
      if (igcm_c2h5 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_c2h5) = rm(l,i_c2h5)
         end do
      end if
      if (igcm_c2h5o2 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_c2h5o2) = rm(l,i_c2h5o2)
         end do
      end if
      if (igcm_c2h5ooh .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_c2h5ooh) = rm(l,i_c2h5ooh)
         end do
      end if
      if (igcm_c2h5oh .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_c2h5oh) = rm(l,i_c2h5oh)
         end do
      end if
      if (igcm_hoch2ch2o2 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hoch2ch2o2) = rm(l,i_hoch2ch2o2)
         end do
      end if
      if (igcm_hoch2ch2o .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hoch2ch2o) = rm(l,i_hoch2ch2o)
         end do
      end if
      if (igcm_ethgly .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ethgly) = rm(l,i_ethgly)
         end do
      end if
      if (igcm_hyetho2h .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hyetho2h) = rm(l,i_hyetho2h)
         end do
      end if
      if (igcm_ch3cho .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3cho) = rm(l,i_ch3cho)
         end do
      end if
      if (igcm_ch3choho2 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3choho2) = rm(l,i_ch3choho2)
         end do
      end if
      if (igcm_ch3cooh .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3cooh) = rm(l,i_ch3cooh)
         end do
      end if
      if (igcm_ch3chohooh .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3chohooh) = rm(l,i_ch3chohooh)
         end do
      end if
      if (igcm_ch3co .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3co) = rm(l,i_ch3co)
         end do
      end if
      if (igcm_ch3cooo .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3cooo) = rm(l,i_ch3cooo)
         end do
      end if
      if (igcm_ch3coooh .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3coooh) = rm(l,i_ch3coooh)
         end do
      end if
      if (igcm_hcoch2o2 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hcoch2o2) = rm(l,i_hcoch2o2)
         end do
      end if
      if (igcm_glyox .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_glyox) = rm(l,i_glyox)
         end do
      end if
      if (igcm_hcoco .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hcoco) = rm(l,i_hcoco)
         end do
      end if
      if (igcm_hooch2cho .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hooch2cho) = rm(l,i_hooch2cho)
         end do
      end if
      if (igcm_hoch2cho .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hoch2cho) = rm(l,i_hoch2cho)
         end do
      end if
      if (igcm_hochcho .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hochcho) = rm(l,i_hochcho)
         end do
      end if
      if (igcm_hoch2co .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hoch2co) = rm(l,i_hoch2co)
         end do
      end if
      if (igcm_hoch2co3 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hoch2co3) = rm(l,i_hoch2co3)
         end do
      end if
      if (igcm_hoch2co2h .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hoch2co2h) = rm(l,i_hoch2co2h)
         end do
      end if
      if (igcm_hcoco2h .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hcoco2h) = rm(l,i_hcoco2h)
         end do
      end if
      if (igcm_hoch2co3h .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hoch2co3h) = rm(l,i_hoch2co3h)
         end do
      end if
      if (igcm_hcoco3h .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hcoco3h) = rm(l,i_hcoco3h)
         end do
      end if
      if (igcm_hcoco3 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hcoco3) = rm(l,i_hcoco3)
         end do
      end if
      if (igcm_ch2choh .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch2choh) = rm(l,i_ch2choh)
         end do
      end if
      
      if (igcm_dust_mass .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_dust_mass) = rm(l,i_dust)
         end do
      end if
      
      if (igcm_cl .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_cl) = rm(l,i_cl)
         end do
      end if
      if (igcm_clo .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_clo) = rm(l,i_clo)
         end do
      end if
      if (igcm_cl2 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_cl2) = rm(l,i_cl2)
         end do
      end if
      if (igcm_oclo .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_oclo) = rm(l,i_oclo)
         end do
      end if
      if (igcm_cl2o2 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_cl2o2) = rm(l,i_cl2o2)
         end do
      end if
      if (igcm_hcl .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hcl) = rm(l,i_hcl)
         end do
      end if
      if (igcm_hocl .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hocl) = rm(l,i_hocl)
         end do
      end if
      if (igcm_cloo .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_cloo) = rm(l,i_cloo)
         end do
      end if
      if (igcm_ch3ocl .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_ch3ocl) = rm(l,i_ch3ocl)
         end do
      end if
      if (igcm_clco .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_clco) = rm(l,i_clco)
         end do
      end if
      if (igcm_clo3 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_clo3) = rm(l,i_clo3)
         end do
      end if
      if (igcm_hclo4 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_hclo4) = rm(l,i_hclo4)
         end do
      end if
      if (igcm_clo4 .ne. 0) then
         do l = 1,lswitch-1
            zycol(l,igcm_clo4) = rm(l,i_clo4)
         end do
      end if



      return
      
      end
c
c*****************************************************************
c
      subroutine chemrates(lswitch, dens, press, t,
     $                     surfdust1d, surfice1d, 
     $                     a001, a002, a003,
     $                     b001, b002, b003, b004, b005, b006,
     $                     b007, b008, b009,
     $                     c001, c002, c003, c004, c005, c006,
     $                     c007, c008, c009, c010, c011, c012,
     $                     c013, c014, c015, c016, c017, c018,
     $                     d001, d002, d003,
     $                     e001, e002, e003,
     $                     h001, h002, h003, h004, h005,
     $                     t001, t002, t003, tau,
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
     $               cab106, cab107,
     $               cl001, cl002, cl003,
     $               cl004, cl005, cl006,
     $               cl007, cl008, cl009,
     $               cl010, cl011, cl012,
     $               cl013, cl014, cl015,
     $               cl016, cl017, cl018,
     $               cl019, cl020, cl021,
     $               cl022, cl023, cl024,
     $               cl025, cl026, cl027,
     $               cl028, cl029, cl030,
     $               cl031, cl032, cl033,
     $               cl034, cl035, cl036,
     $               cl037, cl038, cl039,
     $               cl040, cl041, cl042,
     $               cl043, cl044, cl045, cl046, cl047, cl048,
     $               cl049, cl050, cl051, cl052, cl053, cl054,
     $               cl055, cl056, cl057, cl058, cl059)
c
c*****************************************************************
c
      implicit none
c
#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer lswitch                  ! interface level between chemistries

      real dens(nlayermx)              ! density (cm-3)
      real press(nlayermx)             ! pressure (hpa)
      real t(nlayermx)                 ! temperature (k)
      real surfdust1d(nlayermx)        ! dust surface area (cm^2/cm^3)
      real surfice1d(nlayermx)         ! ice surface area (cm^2/cm^3)
      real tribo                       ! switch for triboelectricity
      real tau                         ! dust opacity at 7 hpa
c
      parameter (tribo   = 0.)         ! switch for triboelectricity
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     outputs:                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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

      real cl001(nlayermx), cl002(nlayermx), cl003(nlayermx) 
     $ ,cl004(nlayermx), cl005(nlayermx), cl006(nlayermx) 
     $ ,cl007(nlayermx), cl008(nlayermx), cl009(nlayermx) 
     $ ,cl010(nlayermx), cl011(nlayermx), cl012(nlayermx) 
     $ ,cl013(nlayermx), cl014(nlayermx), cl015(nlayermx) 
     $ ,cl016(nlayermx), cl017(nlayermx), cl018(nlayermx) 
     $ ,cl019(nlayermx), cl020(nlayermx), cl021(nlayermx) 
     $ ,cl022(nlayermx), cl023(nlayermx), cl024(nlayermx) 
     $ ,cl025(nlayermx), cl026(nlayermx), cl027(nlayermx) 
     $ ,cl028(nlayermx), cl029(nlayermx), cl030(nlayermx) 
     $ ,cl031(nlayermx), cl032(nlayermx), cl033(nlayermx)
     $ ,cl034(nlayermx), cl035(nlayermx), cl036(nlayermx)
     $ ,cl037(nlayermx), cl038(nlayermx), cl039(nlayermx)
     $ ,cl040(nlayermx), cl041(nlayermx), cl042(nlayermx)
     $ ,cl043(nlayermx), cl044(nlayermx), cl045(nlayermx)
     $ ,cl046(nlayermx), cl047(nlayermx), cl048(nlayermx)
     $ ,cl049(nlayermx), cl050(nlayermx), cl051(nlayermx)
     $ ,cl052(nlayermx), cl053(nlayermx), cl054(nlayermx)
     $ ,cl055(nlayermx), cl056(nlayermx), cl057(nlayermx)
     $ ,cl058(nlayermx), cl059(nlayermx)

      real k_clo_clo 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:                                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real ak0, ak1, rate, rate1, rate2, xpo, xpo1, xpo2
      real ef, efmax, lossh2o, lossch4, lossco2
c
      integer l
      real k1a, k1b, k1a0, k1b0, k1ainf
      real x, y, fc, fx
c
c     KIE calculations for 13CH4 + OH / 13CH4 + O1D  
      real dum(2)
      real r_coef, theta
      real kie_oh, kie_o1d
      real, parameter :: KIE_OH_U = 1.00562 ! Mean KIE of 13CH4 + OH 
      real, parameter :: KIE_OH_STD = 0.00042 ! STD of KIE_OH_U 
                                              ! Values from https://doi.org/10.1029/JD095iD13p22455
      real, parameter :: KIE_O1D_U = 1.013 ! Mean KIE of 13CH4 + O(1D)
      real, parameter :: KIE_O1D_STD = 0.0006 ! STD of KIE_O1D_STD

cccc Three-body calculator 
      real k_3body  ! the function 
cccc Organic parameters [CAABA v4.0]
      real k_ch3o2, k_ch3ooh_oh, k_roohro, k_rohro,
     $     k_s, k_RO2sORO2,
     $     k_t, k_ro2rco3, k_ro2poro2,k_ro2ho2_1,
     $     k_ro2ho2_2,  kapho2,
     $     k_co2h, k_ch3co2h
      real, parameter :: f_sooh = 8. 
      real, parameter :: f_soh = 3.44
      real, parameter :: k_rchohch2o2_oh = 0.1
      real, parameter :: f_pch2oh = 1.29 
      real, parameter :: f_tooh = 8.
      real, parameter :: f_toh = 2.68
      real, parameter :: rcoch2o2_oh = 0.15
      real, parameter :: rcoch2o2_ooh = 0.85
      real, parameter :: kdec = 1.e6
      real, parameter :: rco3_o3 = 0.1
      real, parameter :: rco3_oh = 0.69 
      real, parameter :: rco3_ooh = 0.21 
      real, parameter :: f_o = 8.15
      real, parameter :: f_cho = 0.55
      real, parameter :: f_co2h = 1.67
     
      real j_5_2, k0, kinf, A, B 
     


     
     
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute reaction rates
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
      do l = 1,lswitch-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        oxygen compounds
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      a001: o + o2 + co2 -> o3 + co2
c
c        jpl 2003
c
c        co2 efficiency as a third body (2.075)
c        from sehested et al., j. geophys. res., 100, 1995.
C c
         a001(l) = 2.075
     $             *6.0e-34*(t(l)/300.)**(-2.4)*dens(l)
c
c        mulcahy and williams, 1968
c
c        a001(l) = 2.68e-33*(t(l)/298.)**(-2.4)*dens(l)
c
c        nair et al., 1994
c
C        a001(l) = 1.3e-34*exp(724./t(l))*dens(l)
c
ccc      a002: o + o + co2 -> o2 + co2
c
c        Tsang and Hampson, J. Chem. Phys. Ref. Data, 15, 1087, 1986
c
c        a002(l) = 2.5*5.2e-35*exp(900./t(l))*dens(l)
c
c        Campbell and Gray, Chem. Phys. Lett., 18, 607, 1973
c
c        a002(l) = 1.2e-32*(300./t(l))**(2.0)*dens(l)  ! yung expression
c
         a002(l) = 2.5*9.46e-34*exp(485./t(l))*dens(l) ! nist expression
c
c        baulch et al., 1976 confirmed by smith and robertson, 2008
c
c        a002(l) = 2.5*2.76e-34*exp(720./t(l))*dens(l)
c
ccc      a003: o + o3 -> o2 + o2
c
c        jpl 2003
c
         a003(l) = 8.0e-12*exp(-2060./t(l))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        reactions with o(1d)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      b001: o(1d) + co2  -> o + co2
c
c        jpl 2003
c
c        b001(l) = 7.4e-11*exp(120./t(l))
c
c        jpl 2006
c
         b001(l) = 7.5e-11*exp(115./t(l))
c
ccc      b002: o(1d) + h2o  -> oh + oh
c
c        jpl 2003
c
c        b002(l) = 2.2e-10
c
c        jpl 2006
c 
         b002(l) = 1.63e-10*exp(60./t(l))
c    
ccc      b003: o(1d) + h2  -> oh + h
c
c        jpl 2011
c
         b003(l) = 1.2e-10
c    
ccc      b004: o(1d) + o2  -> o + o2
c
c        jpl 2003
c
c        b004(l) = 3.2e-11*exp(70./t(l))
c
c        jpl 2006
c
         b004(l) = 3.3e-11*exp(55./t(l))
c    
ccc      b005: o(1d) + o3  -> o2 + o2
c
c        jpl 2003
c
         b005(l) = 1.2e-10
c    
ccc      b006: o(1d) + o3  -> o2 + o + o
c
c        jpl 2003
c
         b006(l) = 1.2e-10
c    
C ccc      b007: o(1d) + ch4 -> ch3 + oh
C c
C c        jpl 2003
C c
          b007(l) = 1.75e-10*0.75
C c
C ccc      b008: o(1d) + ch4 -> ch3o + h
C c
C c        jpl 2003
C c
          b008(l) = 1.75e-10*0.20
C c
C ccc      b009: o(1d) + ch4 -> ch2o + h2
C c
C c        jpl 2003
C c
          b009(l) = 1.75e-10*0.05
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hydrogen compounds    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      c001: o + ho2 -> oh + o2
c
c        jpl 2003

c 
         c001(l) = 3.0e-11*exp(200./t(l))
c
ccc      c002: o + oh -> o2 + h
c
c        jpl 2011
c 
         c002(l) = 1.8e-11*exp(180./t(l))
c
c        robertson and smith, j. chem. phys. a 110, 6673, 2006
c
c        c002(l) = 11.2e-11*t(l)**(-0.32)*exp(177./t(l))
c
ccc      c003: h + o3 -> oh + o2
c
c        jpl 2003
c 
         c003(l) = 1.4e-10*exp(-470./t(l))
c
ccc      c004: h + ho2 -> oh + oh
c
c        jpl 2003
c 
c        c004(l) = 8.1e-11*0.90
c
c        jpl 2006
c
         c004(l) = 7.2e-11
c
ccc      c005: h + ho2 -> h2 + o2
c
c        jpl 2003
c 
c        c005(l) = 8.1e-11*0.08
c
c        jpl 2006
c
         c005(l) = 6.9e-12
c
ccc      c006: h + ho2 -> h2o + o
c
c        jpl 2003
c 
c        c006(l) = 8.1e-11*0.02
c
c        jpl 2006
c
         c006(l) = 1.6e-12
c
ccc      c007: oh + ho2 -> h2o + o2
c
c        jpl 2003
c 
         c007(l) = 4.8e-11*exp(250./t(l))
c
c        jpl 2003 +20% d'apres canty et al., grl, 2006
c 
c        c007(l) = 4.8e-11*exp(250./t(l))*1.2
c
ccc      c008: ho2 + ho2 -> h2o2 + o2
c
c        jpl 2003
c 
c        c008(l) = 2.3e-13*exp(600./t(l))
c
c        christensen et al., grl, 13, 2002
c
         c008(l) = 1.5e-12*exp(19./t(l))
         
cccccc Water Vapor enhancement: 
cccccc "The Effects of Water Vapor on the CH3O2 Self-Reaction 
cccccc  and Reaction with HO2" - English, A. M. et al. [2008]
cccccc DOI : 10.1021/jp800727a
cccccc c008(l) = 3.4e-13*EXP(594/T(l))*(1. + 2.4e-21*exp(2200/T(l))*[H2O] )
         
c
ccc      c009: oh + h2o2 -> h2o + ho2
c
c        jpl 2003
c 
c        c009(l) = 2.9e-12*exp(-160./t(l))
c
c        jpl 2006
c
         c009(l) = 1.8e-12
c
ccc      c010: oh + h2 -> h2o + h
c
c        jpl 2003
c 
c        c010(l) = 5.5e-12*exp(-2000./t(l))
c
c        jpl 2006
c
         c010(l) = 2.8e-12*exp(-1800./t(l))

ccc      c011: h + o2 + co2 -> ho2 + co2
c
c        jpl 2011
c
         ak0 = 2.5*4.4e-32*(t(l)/300.)**(-1.3)
         ak1 = 7.5e-11*(t(l)/300.)**(0.2)
c
         rate = (ak0*dens(l))/(1. + ak0*dens(l)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(l))/ak1)**2)
         c011(l) = rate*0.6**xpo
         
         
c
ccc      c012: o + h2o2 -> oh + ho2
c
c        jpl 2003
c 
         c012(l) = 1.4e-12*exp(-2000./t(l))
c
ccc      c013: oh + oh -> h2o + o
c
c        jpl 2003
c 
c        c013(l) = 4.2e-12*exp(-240./t(l))
c
c        jpl 2006
c
         c013(l) = 1.8e-12
c
ccc      c014: oh + o3 -> ho2 + o2
c
c        jpl 2003
c
         c014(l) = 1.7e-12*exp(-940./t(l))
c
c        jpl 2000
c
c        c014(l) = 1.5e-12*exp(-880./t(l))
c
c        nair et al., 1994 (jpl 1997)
c
c        c014(l) = 1.6e-12*exp(-940./t(l))
c
ccc      c015: ho2 + o3 -> oh + o2 + o2
c
c        jpl 2003
c
         c015(l) = 1.0e-14*exp(-490./t(l))
c
c        jpl 2000
c
c        c015(l) = 2.0e-14*exp(-680./t(l))
c
c        nair et al., 1994 (jpl 1997)
c
c        c015(l) = 1.1e-14*exp(-500./t(l))
c
ccc      c016: ho2 + ho2 + co2 -> h2o2 + o2 + co2
c
c        jpl 2011
c
         c016(l) = 2.5*2.1e-33
     $              *exp(920./t(l))*dens(l)
c
ccc      c017: oh + oh + co2 -> h2o2 + co2
c
c        jpl 2003
c
         ak0 = 2.5*6.9e-31*(t(l)/300.)**(-1.0)
         ak1 = 2.6e-11*(t(l)/300.)**(0.0)
c
c        jpl 1997
c
c        ak0 = 2.5*6.2e-31*(t(l)/300.)**(-1.0)
c        ak1 = 2.6e-11*(t(l)/300.)**(0.0)
c
c        nair et al., 1994
c
c        ak0 = 2.5*7.1e-31*(t(l)/300.)**(-0.8)
c        ak1 = 1.5e-11*(t(l)/300.)**(0.0)
c
         rate = (ak0*dens(l))/(1. + ak0*dens(l)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(l))/ak1)**2)
         c017(l) = rate*0.6**xpo
c


ccc      c018: h + h + co2 -> h2 + co2
c
c        baulch et al., 1992
c
c        c018(l) = 2.5*8.85e-33*(t(l)/298.)**(-0.6)*dens(l)
c
c        baulch et al., 2005
c
         c018(l) = 2.5*1.8e-30*(t(l)**(-1.0))*dens(l)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        nitrogen compounds
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      d001: no2 + o -> no + o2
c
c        jpl 2003
c
c        d001(l) = 5.6e-12*exp(180./t(l))
c
ccc      jpl 2006
c
         d001(l) = 5.1e-12*exp(210./t(l))
c
ccc      d002: no + o3 -> no2 + o2
c
c        jpl 2003
c
         d002(l) = 3.0e-12*exp(-1500./t(l))
c
ccc      d003: no + ho2 -> no2 + oh
c
c        jpl 2011
c
         d003(l) = 3.3e-12*exp(270./t(l))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        carbon compounds
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      e001: oh + co -> co2 + h
c
c        jpl 2003
c
c        e001(l) = 1.5e-13*(1 + 0.6*press(l)/1013.)
c
c        mccabe et al., grl, 28, 3135, 2001
c
c        e001(l) = 1.57e-13 + 3.54e-33*dens(l)
c
c        atkinson et al. 2006
c
c        e001(l) = 1.44e-13 + 3.43e-33*dens(l)
c
c        joshi et al., 2006
c         
         k1a0 = 1.34*2.5*dens(l)
     $         *1/(1/(3.62e-26*t(l)**(-2.739)*exp(-20./t(l)))
     $         + 1/(6.48e-33*t(l)**(0.14)*exp(-57./t(l))))   ! corrige de l'erreur publi
         k1b0 = 1.17e-19*t(l)**(2.053)*exp(139./t(l))
     $        + 9.56e-12*t(l)**(-0.664)*exp(-167./t(l))
         k1ainf = 1.52e-17*t(l)**(1.858)*exp(28.8/t(l))
     $          + 4.78e-8*t(l)**(-1.851)*exp(-318./t(l))
         x = k1a0/(k1ainf - k1b0)
         y = k1b0/(k1ainf - k1b0)
         fc = 0.628*exp(-1223./t(l)) + (1. - 0.628)*exp(-39./t(l))
     $      + exp(-t(l)/255.)
         fx = fc**(1./(1. + (alog(x))**2))                   ! corrige de l'erreur publi
         k1a = k1a0*((1. + y)/(1. + x))*fx
         k1b = k1b0*(1./(1.+x))*fx
c
         e001(l) = k1a + k1b
c
ccc      e002: o + co + m -> co2 + m
c
c        tsang and hampson, 1986.
c
         e002(l) = 2.5*6.5e-33*exp(-2184./t(l))*dens(l)
c
c        baulch et al., butterworths, 1976.
c
c        e002(l) = 1.6e-32*exp(-2184./t(l))*dens(l)
c
ccc      e003: ch4 + oh -> ch3 + h2o
c
c        jpl 2003
c
c        e003(l) = 2.45e-12*exp(-1775./t(l))
c
c        jpl 2003, three-parameter expression
c
c          e003(l) = 2.80e-14*(t(l)**0.667)*exp(-1575./t(l))

	  
C ! cccccccccccccccccccccc
C ! ccccc 13CH4 Reactions 
C ! cccccccccccccccccccccc
       
C ! c     13ch4 + oh -> ch3 + h2o
C       ! CALL RANDOM_NUMBER(dum)

C       ! r_coef = (-2.*LOG(dum(1)))**0.5

C       ! theta = 2.*3.14159265359*dum(2)

C       ! kie_oh = kie_oh_u + kie_oh_std*r_coef*sin(theta)

     
C ! c     13ch4 + o(1d) -> 0.75*ch3 + 0.75*oh + 0.2*ch3o + 0.2*h + 0.1*ch2o + 0.1*h2
C       ! kie_o1d = kie_o1d_u + kie_o1d_std*r_coef*sin(theta)

C       ! org037(l) = (b007(l) + b008(l) + b009(l))/kie_o1d 

       
       
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Methane Oxidation - CAABA/MECCA  Version 4.0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccc Parameter calculations
         k_ch3o2 = 1.03e-13*exp(365./t(l))
         k_ch3ooh_oh = 5.3e-12*exp(190./t(l))
         k_ch3co2h = 4.e-14*exp(850./t(l))
         k_roohro = 0.6*k_ch3ooh_oh
         k_rohro = 2.1e-18*t(l)*t(l)*exp(-85./t(l))
         k_s = 4.5e-18*t(l)*t(l)*exp(253./t(l))
         k_RO2sORO2 = 2.*(7.7E-15*exp(1330./t(l))*k_ch3o2)**0.5
         k_t = 2.12e-18*t(l)*t(l)*exp(696./t(l))
         k_ro2rco3 = 2.*2.e-12*exp(500./t(l))
         k_ro2poro2 = 2.*7.5e-13*exp(500./t(l))
         k_ro2ho2_1 = 2.91e-13*exp(1300./t(l))*(1.-exp(-0.245*1.))
         k_ro2ho2_2 = 2.91e-13*exp(1300./t(l))*(1.-exp(-0.245*2.))
         kapho2 = 5.e-13*exp(980./t(l))
         k_co2h = 0.7*k_ch3co2h

         
cccccccc
cc ch4 
cccccccc
c        cab001: ch4 + oh -> ch3 + h2o
         cab001(l) = 1.85e-20*EXP(2.82*LOG(t(l)) - 987./t(l))

c        cab002: ch4 + o -> 0.51*ch3 + 0.51*oh + 0.49*ch3o + 0.49*h 
         cab002(l) = 6.03e-18*(t(l)**2.17)*EXP(-3619./t(l))
         
cccccccc 
cc ch3 
cccccccc 
c        cab003: ch3 + o2 + M -> ch3o2 + M 
         cab003(l) = k_3body(t(l),dens(l), 7.e-31, 3.,
     $              1.8e-12, -1.1, 0.33, .True.)
          
c        cab004: ch3 + o3 -> 0.956*hcho + 0.956*h + 0.044*ch3o + o2 
         cab004(l) = 5.1e-12*EXP(-210./t(l))
         
c        cab005: ch3 + o -> 0.83*hcho + 0.83*h + 0.17*co + 0.17*h2 + 0.17*h 
         cab005(l) = 1.3e-10
         
ccccccccc 
cc ch3o2 
ccccccccc 
c        cab006: ch3o2 + ho2 -> ch3ooh + o2
         cab006(l) = 3.8e-13*exp(780./t(l))
     $             /( 1. + 1./(498.*exp(1160./t(l))))

c        cab007: ch3o2 + ho2 -> hcho + h2o + o2 
         cab007(l) = 3.8e-13*exp(780./t(l))
     $             /( 1. + (498.*exp(-1160./t(l))))
         
c        cab008: ch3o2 + ro2 -> ch3o + 0.5*o2
         cab008(l) = 7.4e-13*exp(-520./t(l))*2. ! *[RO2]
         
C        cab009: ch3o2 + ro2 -> 0.5*hcho + 0.5ch3oh + 0.5*o2  
         cab009(l) = (k_ch3o2 - 7.4e-13*exp(-520./t(l)))*2. ! *[RO2]
         
c        cab010: ch3o2 + o3 -> ch3o + 2.*o2 
         cab010(l) = 2.9e-16*exp(-1000./t(l))
         
c        cab011: ch3o2 + oh -> ch3o + ho2 
         cab011(l) = 1.4e-10 
         
c        cab012: ch3o2 + o -> ch3o + o2 
         cab012(l) = 4.3e-11 
         
ccccccccc
cc ch3oh 
ccccccccc
c       cab013: ch3oh + oh -> 0.85*hcho + 0.85*ho2 + 0.15*ch3o h2o 
        cab013(l) = 6.38e-18*(t(l)**2.)*exp(144./t(l))
         
ccccccccc
cc ch3ooh 
ccccccccc
c       cab014: ch3ooh + oh -> 0.6*ch3o2 + 0.4*hcho + 0.4*oh + h2o 
        cab014(l) = k_ch3ooh_oh
         
ccccccc
cc ch3o
ccccccc
c       cab015: ch3o + o2 -> ho2 + hcho 
        cab015(l) = 1.3e-14*exp(-663./t(l))

c       cab016: ch3o + o3 -> ch3o2 + o2 
        cab016(l) = 2.53e-14
        
c       cab017: ch3o + o -> 0.75*ch3 + 0.75*o2 + 0.25*hcho + 0.25*oh 
        cab017(l) = 2.5e-11
c        
ccccccc
cc hcho
ccccccc
c       cab018: hcho + oh -> hco + h2o
        cab018(l) = 5.5e-12*exp(-125/t(L)) ! [JPL Atmospheric Handbook Evaluation 18]

c       cab019: hcho + ho2 -> hoch2o2 
        cab019(l) = 9.7e-15*exp(625./t(l)) 
        
c       cab020: hcho + o -> hco + oh 
        cab020(l) = 2.99e-11*exp( -1.28e4/(8.314*t(l))) ! [Evaluated chemical kinetic data for the reactions of atomic oxygen O(3P) with saturated organic compounds in the gas phase]
        
ccccccc 
cc hco 
ccccccc
c      cab021: hco + o -> co + oh 
       cab021(l) = 5.e-11 ! [Evaluated kinetic data for combustion modelling]
     
c      cab022: hco + ch3 -> ch4 + co [The application of an intracavity laser spectroscopy method for elementary processes study in gas-phase reactions]
       cab022(l) = 4.4e-11 
    
c      cab023: hco + ch3 -> ch3cho [ Investigation of the Recombination of the HCO Radical in an Atmosphere of Argon and Helium by the Method of Internal Resonator Laser Spectroscopy]
       cab023(l) = 4.42e-11 
       
c      cab024: 2.*hco -> hcho + co [Quantitative Detection of HCO Behind Shock Waves: The Thermal Decomposition of HCO]
       cab024(l) =  4.48e-11 
       
c      cab025: hco + oh -> co + h2o [Evaluated kinetic data for combustion modelling]
       cab025(l) = 1.69e-10
       
c      cab026: hco + o2 -> co + ho2 [JPL Handbook Evaluation 18]
       cab026(l) = 5.2e-12 
       
c      cab027: hco + h -> co + h2 [Quantitative Detection of HCO Behind Shock Waves: The Thermal Decomposition of HCO]
       cab027(l) = 1.83e-10
       
cccccccccc 
cc hoch2o2 
cccccccccc 
c     cab028: hoch2o2 -> hcho + ho2 
      cab028(l) = 2.4e12*exp(-7000./t(l))
      
c     cab029: hoch2o2 + ho2 -> 0.5*hoch2ooh + 0.5*hcooh + 0.2*oh + 0.2*ho2 + 0.3*h2o + 0.8*o2 
      cab029(l) = 5.6e-15*exp(2300./t(l))

c     cab030: hoch2o2 + ro2 -> hcooh + ho2 
      cab030(l) = 2.*(k_ch3o2*5.5e-12)**0.5
      
c     cab031: hoch2o2 + ro2 -> 0.5*hcooh + 0.5*hoch2oh + 0.5*o2 
      cab031(l) = 2.*(k_ch3o2*5.7e-14*exp(750./t(l)))**0.5
      
cccccccccc
cc hcooh 
cccccccccc
c     cab032: hcooh + oh  -> co2 + ho2 + h2o 
      cab032(l) = 2.94e-14*exp(786./t(l)) 
     $            + 9.85e-13*exp(-1036./t(l))
     
cccccccccc
cc hoch2ooh
cccccccccc
c     cab033: hoch2ooh + oh -> hoch2o2 
      cab033(l) = k_roohro
      
c     cab034: hoch2ooh + oh -> hcooh + h2o + oh     
      cab034(l) = k_rohro + k_s*f_sooh*f_soh
     
cccccccccc
cc hoch2oh
cccccccccc
c     cab035: hoch2oh + oh -> ho2 + hcooh + h2o     
      cab035(l) = 2.*k_rohro + k_s*f_soh*f_soh
      
cccccccccc 
cc c2h6 
cccccccccc 
c     cab036: c2h6 + oh -> c2h5 + h2o [JPL Data Evaluation 18]
      cab036(l) = 7.66e-12*exp(-1020./t(l))
      
c     cab037: c2h6 + o -> c2h5 + oh [Chemical kinetic data sheets for high-temperature reactions. Part II] 
      cab037(l) = (2.21e-15*(t(l)/298.)**6.5)*exp(-132./t(l))

c     cab038: ch3 + ch3 + M -> c2h6 [Rate constant for the reaction CH3 + CH3 → C2H6 at T = 155 K and model calculation of the CH3 abundance in the atmospheres of Saturn and Neptune]
      kinf = 4.504e-11*exp(70.12/t(l))
      k0   =  9.682e-28*exp(569.4/t(l)) 
      j_5_2 = 163.*exp(-1001./t(l))
      
      if ( j_5_2 <= 1.1) then 
          j_5_2 = 1.1 
      endif 
      
      cab038(l) = (-(kinf + k0*dens(l)*2.5) +
     $          sqrt(
     $               (kinf + k0*dens(l)*2.5)**2. 
     $             + 4.*(j_5_2 -1.)*kinf*k0*dens(l)*2.5
     $              ) )
     $            /(2.*(j_5_2-1.))
     
 
cccccccccc 
cc c2h5 
cccccccccc 
c     cab039: c2h5 + o2 + M -> c2h5o2 + M [JPL Data Evaluation 18]
C        kinf = 8e-12*(t(l)/300.)
C        k0 = 1.5e-28*(300./t(l))**3.
       
C        cab039(l) = (k0*dens(l)*2.5/(1. + (k0*dens(l)*2.5/kinf))) 
C      $           *0.6
C      $           **(1./(1. + log10(k0*dens(l)*2.5/kinf)**2.) )
      cab039(l) = 1.5e-28*((298./t(l))**3.)*dens(l)*2.5



c     cab040: c2h5 + c2h5 -> c2h6 + c2h4 [Temperature coefficients of the rates of Cl atom reactions with C2H6, C2H5, and C2H4. The rates of disproportionation and recombination of ethyl radicals]
      cab040(l) = 2.01e-12 
      
c     cab041: c2h5 + o2 -> c2h4 + ho2 [JPL Data Evaluation 18 (UPPER LIMIT) ]
      cab041(l) = 1.9e-14 
      
c     cab042: c2h5 + h -> ch3 + ch3 [ Kinetics of the reaction of methyl radicals with oxygen]
      cab042(l) = 7.95e-11*exp(-132./t(l))
      
ccccccccc 
cc c2h5o2 
ccccccccc
c     cab043: c2h5o2 + ho2 -> c2h5ooh + o2 
      cab043(l) = 7.5e-13*exp(700./t(l))
      
c     cab044: c2h5o2 + ro2 -> 0.8*ch3cho + 0.6*ho2 + 0.2*c2h5oh + o2
      cab044(l) = 2.*(7.6e-14*k_ch3o2)**0.5
      
ccccccccc
cc c2h5ooh
ccccccccc 
c     cab045: c2h5ooh + oh -> c2h5o2 + h2o 
      cab045(l) = k_roohro

c     cab046: c2h5ooh + oh -> ch3cho + oh 
      cab046(l) = k_s*f_sooh
      
ccccccccc
cc c2h5oh
ccccccccc 
c     cab047: c2h5oh + oh -> 0.95*c2h5o2 + 0.95*ho2 + 0.05*hoch2ch2o2 + h2o 
      cab047(l) = 3.e-12*exp(20./t(l))
 
ccccccccc
cc hoch2ch2o2
ccccccccc 
c     cab048: hoch2ch2o2 + ro2 -> 0.6*hoch2ch2o + 0.2*hoch2cho + 0.2*ethgly
      cab048(l) = 2.*(7.8e-14*exp(1000./t(l))*k_ch3o2)**(0.5)
      
c     cab049: hoch2ch2o2 + ho2 -> HYETHO2H
      cab049(l) = 1.53e-13*exp(1300./t(l))*(1.-k_rchohch2o2_oh)
 
c     cab050: hoch2ch2o2 + ho2 -> hoch2ch2o + oh 
      cab050(l) = 1.53e-13*exp(1300./t(l))*k_rchohch2o2_oh

ccccccccc
cc hoch2ch2o
ccccccccc 
c     cab051: hoch2ch2o + o2 -> ho2 + hoch2cho 
      cab051(l) = 6e-14*exp(-550./t(l))
      
c     cab052: hoch2ch2o -> ho2 + hcho + hcho 
      cab052(l) = 9.5e13*exp(-5988./t(l))

ccccccccc
cc ethgly
ccccccccc 
c     cab053: ethgly + oh -> hoch2cho + ho2 + h2o 
      cab053(l) = 2.*k_s*f_soh*f_pch2oh + 2.*k_rohro

ccccccccc
cc hyetho2h
ccccccccc 
c     cab054: hyetho2h + oh -> hoch2ch2o2 + h2o 
      cab054(l) = k_roohro
      
c     cab055: hyetho2h + oh -> hoch2cho + oh + h2o
      cab055(l) = k_s*f_sooh*f_pch2oh
      
c     cab056: hyetho2h + oh -> hooch2cho + ho2 + h2o 
      cab056(l) = k_s*f_soh*f_pch2oh + k_rohro
      
ccccccccc
cc ch3cho
ccccccccc 
c     cab057: ch3cho + oh -> ch3co + h2o 
      cab057(l) = 4.4e-12*exp(365./t(l))*0.95

c     cab058: ch3cho + oh -> hcoch2o2 + h2o 
      cab058(l) = 4.4e-12*exp(365./t(l))*0.05
      
c     cab059: ch3cho + ho2 -> ch3choho2 
      cab059(l) = 3.46e12*exp(-12500./(t(l)*1.98))
     $          /(6.34e26*exp(-1470./(t(l)*1.98)))

c     cab060: ch3cho + hcooh -> ch2choh + hcooh  
      cab060(l) = (1.17e-19*t(l)**2.209)*exp(-556./(1.987*t(l)))
      
ccccccccc 
cc ch2choh 
ccccccccc 
c     cab061: ch2choh + oh -> hcooh + oh + hcho 
      cab061(l) = 4.3e-11 
      
c     cab062: ch2choh + oh -> hoch2cho + ho2 
      cab062(l) = 2.4e-11 
      
c     cab063: ch2choh + hcooh -> ch3cho + hcooh 
      cab063(l) = (4.67e-26*t(l)**3.286)*exp(-556./(1.987*t(l)))

ccccccccc 
cc ch3choho2
ccccccccc 
c     cab064: ch3choho2 -> ch3cho + ho2 
      cab064(l) = 3.46e12*exp(-12500./(t(l)*1.98))
     
c     cab065: ch3choho2 + ho2 -> 0.5*ch3chohooh + 0.3*ch3cooh + 0.2*ch3 + 0.2*hcooh + 0.2*oh + o2
      cab065(l) = 5.6e-15*exp(2300./t(l))
     
c     cab066: ch3choho2 + ro2 -> ch3 + hcooh + oh 
      cab066(L) = k_RO2sORO2

ccccccccc 
cc ch3cooh
ccccccccc 
c     cab067: ch3cooh + oh -> ch3 + co2 + h2o
      cab067(l) =  4.e-14*exp(850./t(l))
     
ccccccccc 
cc ch3chohooh
ccccccccc 
c     cab068: ch3chohooh + oh -> ch3cooh + oh 
      cab068(l) = k_t*f_tooh*f_toh + k_rohro
      
c     cab069: ch3chohooh + oh -> ch3choho2 
      cab069(l) = k_roohro 
      
ccccccccc 
cc ch3c(o)
ccccccccc 
c     cab070: ch3c(o) + o2 + M -> ch3c(o)oo + M 
      cab070(l) = 5.1e-12*(1. - 1./(1. + 9.4e-18*dens(l)*2.5))
      
c     cab070: ch3c(o) + o2 + M -> oh + hcho + co + M 
      cab071(l) = 5.1e-12/(1.+9.4e-18*dens(l)*2.5)
 
ccccccccc 
cc ch3c(o)oo
ccccccccc 
c     cab072: ch3c(o)oo + ho2 -> oh + ch3 + co2 
      cab072(L) = 5.2e-13*exp(980./t(l))*1.507*0.61

c     cab073: ch3c(o)oo + ho2 -> ch3c(o)ooh
      cab073(L) = cab072(l)*0.23/0.61

c     cab074: ch3c(o)oo + ho2 -> ch3cooh + o3 
      cab074(L) = cab074(l)*0.16/0.61

c     cab075: ch3c(o)oo + ro2 -> ch3 + co2 
      cab075(l) = k_ro2rco3*0.9
     
c     cab076: ch3c(o)oo + ro2 -> ch3cooh 
      cab076(l) = k_ro2rco3*0.1 
      
ccccccccc 
cc ch3c(o)ooh
ccccccccc 
c     cab077: ch3c(o)ooh + oh -> ch3c(o)oo + h2o 
      cab077(l) = k_roohro 
     
ccccccccc 
cc hcoch2o2
ccccccccc 
c     cab078: hcoch2o2 + ro2 -> 0.6*hcho + 0.6*co + 0.6*ho2 + 0.2*glyox + 0.2*hoch2cho 
      cab078(l) = k_ro2poro2

c     cab079: hcoch2o2 + ho2 -> hooch2cho 
      cab079(l) = k_ro2ho2_2*rcoch2o2_ooh 
     
c     cab080: hcoch2o2 + ho2 -> hcho + co + ho2 + oh 
      cab080(l) = k_ro2ho2_2*rcoch2o2_oh 
      
ccccccccc 
cc glyox
ccccccccc 
c     cab081: glyox + oh -> hcoco + h2o 
      cab081(l) = 3.1e-12*exp(340./t(l))
      
ccccccccc 
cc hcoco 
ccccccccc 
c     cab082: hcoco -> co + co + ho2 [CAABA 4.0]
c      cab082(l) = 7.e11*exp(-3160./t(l)) ! + 5.e-12*cc(l,i_o2)  
c             hcoco -> hco + co [The Atmospheric Chemistry of the HC(O)CO Radical]
       cab082(l) = 1.4e12*exp(-3159./(l))

c     cab083: hcoco + o2 -> hcoco3 
      cab083(l) = 5.e-12*3.2*exp(-550./t(l)) 
 
c     cab084: hcoco + o2 -> oh + co + co2 
      cab084(l) = 5.e-12*(1. - 3.2*exp(-550./t(l)))

ccccccccccc
cc hooch2cho
ccccccccccc      
c     cab085: hooch2cho + oh -> hcoch2o2 ! + h2o
      cab085(l) = k_roohro 

c     cab086: hooch2cho + oh -> hcho + co + oh ! + h2o 
      cab086(l) = 0.8*8.e-12 
     
c     cab087: hooch2cho + oh -> glyox + oh ! + h2o 
      cab087(l) = k_s*f_sooh*f_cho*0.55 
       
ccccccccccc
cc hoch2cho
ccccccccccc 
c     cab088: hoch2cho + oh -> hoch2co + h2o 
      cab088(l) = 8.e-12*0.8
     
c     cab089: hoch2cho + oh -> hochcho + h2o 
      cab089(l) = 8.e-12*0.2

ccccccccccc
cc hochcho
ccccccccccc 
c     cab090: hochcho -> glyox + ho2 
      cab090(l) = kdec 
      
ccccccccccc
cc hoch2co
ccccccccccc
c     cab091: hoch2co + o2 + M -> hoch2co3 
      cab091(l) = 5.1e-12*(1. - 1./(1. + 1.85e-18*dens(l)*2.5))
      
c     cab092: hoch2co + o2 + M -> oh + hcho + co2 + M 
      cab092(l) = 5.1e-12/(1. + 1.85e-18*dens(l)*2.5)

ccccccccccc
cc hoch2co3
ccccccccccc
c     cab093: hoch2co3 + ro2 -> hcho + co2 + ho2 
      cab093(l) = k_ro2rco3*0.9

c     cab094: hoch2co3 + ro2 -> hoch2co2h 
      cab094(l) = k_ro2rco3*0.1

c     cab095: hoch2co3 + ho2 -> hcho + oh + ho2 + co2
      cab095(l) = kapho2*rco3_oh

c     cab096: hoch2co3 + ho2 -> hoch2co3h
      cab096(l) = kapho2*rco3_ooh 
      
c     cab097: hoch2co3 + ho2 -> hoch2co2h + o3
      cab097(l) = kapho2*rco3_o3 
      
ccccccccccc
cc hoch2co2h
ccccccccccc
c     cab098: hoch2co2h + oh -> 0.09*hcho + 0.9*co2 + 0.91*hcoco2h + ho2 + h2o 
      cab098(l) = k_co2h + k_s*f_soh*f_co2h 

ccccccccccc
cc hcoco2h
ccccccccccc
c     cab099: hcoco2h + oh -> co + ho2 + co2 + h2o 
      cab099(l) = k_co2h + k_t*f_o*f_co2h 

cccccccccccc
cc hoch2co3h
cccccccccccc
c     cab100: hoch2co3h + oh -> hoch2co3 + h2o 
      cab100(l) = k_roohro
    
c     cab101: hoch2co3h + oh -> hcoco3h + ho2 
      cab101(l) = k_s*f_soh*f_co2h 
  
cccccccccccc  
cc hcoco3h
cccccccccccc
c     cab102: hcoco3h + oh -> hcoco3 + h2o 
      cab102(l) = k_roohro 
      
c     cab103: hcoco3h + oh -> co + co2 + h2o + oh 
      cab103(l) = k_t*f_o*f_co2h 
      
cccccccccccc 
cc hcoco3  
cccccccccccc 
c     cab104: hcoco3 + ro2 -> co + ho2 + co2 
      cab104(l) = k_ro2rco3*0.9
c     cab105: hcoco3 + ro2 -> hcoco2h + o2 
      cab105(l) = k_ro2rco3*0.1 
c     cab106: hcoco3 + ho2 -> ho2 + co + co2 + oh 
      cab106(l) = kapho2
       
ccccccccccccccccccccccccccccccc 
cc Additional Rate Coefficients
ccccccccccccccccccccccccccccccc 
c    cab107: CH3 + OH + M -> CH3OH 
      cab107(l) = 2.5*dens(l)*2.48E-27*(298./t(l))**(3.8)
      
ccccccccccccc 
cc c2h4 
ccccccccccccc 
c     cab107: c2h4 + o3 -> hcho + ch2oo*
c     cab107(L) = 9.1e-15*exp(-2580./t(l)) 
c     cab108: c2h4 + oh + M(o2) -> hoch2ch2o2
c     cab108(l) = k_3body( t(l), cc(l,i_o2), 8.6e-29, 
c                         3.1, 9.e-1, 0.85, 0.48, .true.)

ccccccccccccc 
cc ch2oo*
ccccccccccccc 
c    cab109: ch2oo* -> 0.37*ch2oo + 0.47*co + 0.47*h2o + 0.16*ho2 + 0.16*co + 0.16*oh 
c    cab109(l) = kdec 
 
ccccccccccccc 
cc ch2oo
ccccccccccccc 
c    cab110: ch2oo -> co + ho2 + oh 
c    cab110(l) = 1.124e14*exp(-10000./t(l))

c    cab111: ch2oo + h2o -> hoch2ooh 
c    cab111(l) = k_ch2oo_no2*3.6e-6 

c    cab112: ch2oo + co -> hcho + co2 
c    cab112(l) = 3.6e-14 

c    cab113: ch2oo + hcooh -> 2.*hcooh 
c    cab113(l) = 1.e-10 

c    cab114: ch2oo + hcho -> 2.*LCARBON
C    cab114(l) = 1.7e-12 

C    cab115: ch2oo + ch3oh -> 2.*LCARBON 
c    cab115(l) = 5.e-12 

c    cab116: ch2oo + ch3o2 -> 2.*LCARBON
c    cab116(L) = 5.e-12 

c    cab117: ch2oo + ho2 -> LCARBON 
c    cab117(l) = 5.e-12

c    cab118: ch2oo + o3 -> hcho + 2.*o2 
c    cab118(l) = 1.e-12 

c    cab119: ch2oo + ch2oo -> 2.*hcho + o2 
c    cab119(l) = 6.e-11 

       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        heterogenous chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        k = (surface*v*gamma)/4 (s-1)
c        v = 100*sqrt(8rt/(pi*m))  (cm s-1)
c
ccc      h001: ho2 + ice -> products
c
c        cooper and abbatt, 1996: gamma = 0.025
c
         h001(l) = surfice1d(l)
     $            *100.*sqrt(8.*8.31*t(l)/(33.e-3*pi))*0.025/4.
c
c        h002: oh + ice -> products
c
c        cooper and abbatt, 1996: gamma = 0.03
c
         h002(l) = surfice1d(l)
     $            *100.*sqrt(8.*8.31*t(l)/(17.e-3*pi))*0.03/4.
c
c        h003: ho2 + dust -> products
c
c        jacob, 2000: gamma = 0.2
c        see dereus et al., atm. chem. phys., 2005
c
c        h003(l) = surfdust1d(l)
c    $            *100.*sqrt(8.*8.31*t(l)/(33.e-3*pi))*0.2/4.
         h003(l) = 0.     ! advised
c
ccc      h004: h2o2 + ice -> products
c
c        gamma = 1.e-3    test value
c
c        h004(l) = surfice1d(l)
c    $            *100.*sqrt(8.*8.31*t(l)/(34.e-3*pi))*0.001/4.
         h004(l) = 0.     ! advised
c
c        h005: h2o2 + dust -> products
c
c        gamma = 5.e-4
c        see dereus et al., atm. chem. phys., 2005
c
         h005(l) = surfdust1d(l)
     $            *100.*sqrt(8.*8.31*t(l)/(34.e-3*pi))*5.e-4/4.
         h005(l) = 0.     ! advised


c     Chlorine Scheme 22/04/2020
c     ==========================

      k_clo_clo = 2.5*dens(l)*2.0E-32*(300./t(l))**(4.)

c     cl001 : cl + o3 -> clo + o2 
c     MECCA 4.0
      cl001(l) = 2.8E-11*EXP(-250./T(L))

c     cl002 : clo + o -> cl + o2 
c     MECCA 4.0
      cl002(l) = 2.5E-11*EXP(110./T(L))

c     cl003 : clo + clo -> cl2 + o2 
c     MECCA 4.0
      cl003(l) = 1.0E-12*EXP(-1590./T(L))

c     cl004 : clo + clo -> 2*cl + o2 
c     MECCA 4.0
      cl004(l) = 3.0E-11*EXP(-2450./T(L))

c     cl005 : clo + clo -> cl + oclo  
c     MECCA 4.0
      cl005(l) = 3.5E-13*EXP(-1370./T(L))

c     cl006 : clo + clo -> cl2o2 
c     MECCA 4.0
      cl006(l) = K_CLO_CLO

c     cl007 : cl2o2 -> clo + clo  
c     MECCA 4.0
      cl007(l) = K_CLO_CLO/(2.16E-27*EXP(8537./T(L)))

c     cl008 : cl + h2 -> hcl + h  
c     MECCA 4.0
      cl008(l) = 3.9E-11*EXP(-2310./T(L))

c     cl009 : cl + ho2 -> hcl + o2 
c     MECCA 4.0
      cl009(l) = 4.4E-11 - 7.5E-11*EXP(-620./T(L))

c     cl010 : cl + ho2 -> clo + oh 
c     MECCA 4.0
      cl010(l) = 7.5E-11*EXP(-620./T(L))

c     cl011 : cl + h2o2 -> hcl + ho2 
c     MECCA 4.0
      cl011(l) = 1.1E-11*EXP(-980./T(L))

c     cl012 : clo + oh -> 0.94*cl + 0.94*ho2 + 0.06*hcl + 0.06*o2 
c     MECCA 4.0
      cl012(l) = 7.3E-12*EXP(300./T(L))

c     cl013 : clo + ho2 -> hocl + o2 
c     MECCA 4.0
      cl013(l) = 2.2E-12*EXP(340./T(L))

c     cl014 : hcl + oh -> cl + h2o 
c     MECCA 4.0
      cl014(l) = 1.7E-12*EXP(-230./T(L))

c     cl015 : hocl + oh -> clo + h2o 
c     MECCA 4.0
      cl015(l) = 3.0E-12*EXP(-500./T(L))

c     cl016 : cl + ch4 -> hcl + ch3 
c     MECCA 4.0
      cl016(l) = 6.6E-12*EXP(-1240./T(L))

c     cl017 : cl + hcho -> hcl + hco 
c     Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry
      cl017(l) = 8.2E-11*EXP(-34./T(L))

c     cl018 : cl + ch3ooh -> hcho + hcl + oh 
c     MECCA 4.0
      cl018(l) = 5.9E-11

c     cl019 : clo + ch3o2 -> ch3o + cloo 
c     Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry
      cl019(l) = 4.92E-12*EXP(-261.5/T(L))

c     cl020 : clo + ch3o2 -> ch3ocl + o2 
c     Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry
      cl020(l) = 2.6E-13*EXP(260./T(L))

c     cl021 : clo + ch3o2 -> ch3o + oclo  
c     Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry
      cl021(l) = 1.E-15

c     cl022 : cl + ch3o2 -> 0.5*clo + 0.5*ch3o + 0.5*hcl + 0.5*ch2oo 
c     MECCA 4.0
      cl022(l) = 1.6E-10

c     cl023 : co + cl + M -> ClCO 
c     Evaluated kinetic and photochemical data for atmospheric chemistry: Volume III - gas phase reactions of inorganic halogens
      cl023(l) = dens(l)*2.5*1.33E-33*(298./T(L))**(3.8)

c     cl024 : cloo + cl -> 2*clo 
c     Chemical kinetics and photochemical data for use in stratospheric modeling. Evaluation number 12
      cl024(l) = 1.2E-11

c     cl025 : cloo + cl -> o2 + cl2 
c     Chemical kinetics and photochemical data for use in stratospheric modeling. Evaluation number 12
      cl025(l) = 2.31E-10

c     cl026 : cloo -> o2 + cl 
c     Evaluated kinetic and photochemical data for atmospheric chemistry: Volume III - gas phase reactions of inorganic halogens
      cl026(l) = 2.8E-10*EXP(-1820./T(L))*dens(l)

c     cl027 : HClO4 + OH -> ClO4 + H2O
c     An Ab Initio Chemical Kinetic Study on the Reactions of H, OH, and Cl with HOClO3
      cl027(l) = exp(-1664./t(l))*1.24e-8/(t(l)**2.9)

C c     cl028 : cl + o2 + M -> ClOO + M  
c       Evaluated kinetic and photochemical data for atmospheric chemistry: Volume III - gas phase reactions of inorganic halogens
      cl028(l) = 2.5*dens(l)*1.44e-33*(298./t(l))**3.9

c     cl029 : cl + ch3ocl -> cl2 + ch3o 
c     Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry
      cl029(l) = 5.2E-11

c     cl030 : cl + cl2o2 -> cl2 + cloo 
c     Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry
      cl030(l) = 1.E-10

c     cl031 : cl + ch3ocl -> hcl + ch2ocl  
c     Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry
      cl031(l) = 9.15E-12

c     cl032 : oh + ch3ocl -> ch3o + hocl
c     Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry
      cl032(l) = 2.4E-12*EXP(-360./T(L))

c     cl033 : clco -> co + cl 
c     Evaluated kinetic and photochemical data for atmospheric chemistry: Volume III - gas phase reactions of inorganic halogens
      cl033(l) = 4.1E-10*EXP(-2960./T(L))*dens(l)

c     cl034 : cl2 + o1d -> o + cl2
c     JPL 2015
      cl034(l) = 0.25*2.7e-10

c     cl035 : cl2 + o1d -> cl + clo
c     JPL 2015
      cl035(l) = 0.75*2.7e-10 

c     cl036 : cl2 + oh -> hocl + cl 
c     JPL205 
      cl036(l) = 2.6e-12*exp(-1100./t(l))

c     cl037 : cl2 + h -> hcl + cl 
c     UV Absorption Spectrum and Self-Reaction Kinetics of the Cyclohexadienyl Radical, and Stability of a Series of Cyclohexadienyl-Type Radicals
      cl037(l) = 8.e-11*exp(-3460./(8.314*t(l)))

c     cl038 : cl2 + ch3 -> ch3cl + cl 
c     Rate constants and hydrogen isotope substitution effects in the CH3+HCl and CH3+Cl2 reaction
      cl038(l) = exp(5520./(8.314*t(l)))*1.67e-13*(T(l)/298.)**(2.52)

c     cl039 : hcl + o1d -> 0.12*o + 0.12*hcl + 0.22*h + 0.22*clo + 0.66*cl + 0.66*oh 
c     JPL 2015
      cl039(l) = 1.5e-10

c     cl040 : hcl + o -> oh + cl 
c     JPL 2015
      cl040(l) = 1.1e-11*exp(-3300./(t(l)))

c     cl041 : hcl + h -> h2 + cl 
c     Kinetic studies of reactions of hydrogen atoms with HCl, Cl2 and NOCl, and chlorine atoms with H2 and NOCl
      cl041(l) = 2.01e-11*exp(-14880./(8.314*t(l)))

c     cl042 : hocl + o -> clo + oh
c     Evaluated kinetic and photochemical data for atmospheric chemistry: Volume III - gas phase reactions of inorganic halogens
      cl042(l) = 1.7e-13

c     cl043 : hcooh + cl -> hcl + products
c     Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry. Supplement IV, IUPAC Subcommittee on Gas Kinetic Data Evaluation for Atmospheric Chemistry
      cl043(l) = 2.01e-13

C     cl044 : Cl + O3 + M -> ClO3 + M 
      cl044(l) = 2.5e-31*dens(l)

c     cl045 : ClO + ClO3 -> ClOO + OClO 
      cl045(l) = 1.85e-18*exp(-2417./t(l))*(t(l)**2)*0.28 

c     cl046 : ClO + ClO3 -> OClO + OClO 
      cl046(l) = 1.42e-18*exp(-2870./t(l))*(t(l)**2)*0.11
 
c     cl407 : ClO + ClO3 + M -> Cl2O4 + M
      k0 = 3.57e-1*exp(10079./t(l))/(t(l)**10.19)
      kinf = 1.43e-10*exp(-82./t(l))*t(l)**0.09
      
      A = k0*dens(l)/( 1. + (k0*dens(l))/kinf )
      B = 1. + log10( k0*dens(l)/kinf )**2.

      cl047(l) = A*0.6**(1./B)

c     cl048 : ClO3 + OH -> HClO4 
      cl048(l) = 6.67e-13

c     cl049 : ClO3 + OH + M -> HClO4
      k0 = 1.9e36*exp(-5542./t(l)) 
      k0 = k0/(t(l)**15.3)

      kinf = 3.2e-10*t(l)**(0.07)*exp(-25./t(l))
      A = k0*dens(l)/( 1. + (k0*dens(l))/kinf )
      B = 1. + log10( k0*dens(l)/kinf )**2.

      cl049(l) = A*0.6**(1./B)

c     cl050 : ClO3 + OH -> OClO + HO2 
      cl050(l) = 2.1e-10*exp(-18./t(l))*t(l)**0.09

c     cl051 : OClO + O + M -> ClO3 
      k0 = 7.25e-31*(300./t(l))**(3.1) 
      kinf = 8.3e-12

      A = k0*dens(l)/( 1. + (k0*dens(l))/kinf )
      B = 1. + log10( k0*dens(l)/kinf )**2.

      cl051(l) = A*0.6**(1./B)

c     cl052 : OClO + O3 -> ClO3 + O2 
      cl052(l) = 2.1e-12*exp(-4700./t(l))

c     cl053 : clo4 + cl -> clo3 + clo 
      cl053(l) = 1.98e-10*exp(-49./t(l))*( t(l)/298. )**0.16

c     cl054 : clo4 + hocl -> oh + hclo4 
      cl054(l) = 2.57e-14*exp(-1017./t(l))*(t(l)/298.)**1.73

c     cl055 : hocl + cl -> cl2 + oh
      cl055(l) = 3.196e-12*exp(-130./t(l)) 

c     cl056 : hocl + cl -> hcl + clo 
      cl056(l) = 2.04e-13*exp(-130./t(l))

c     cl057 : oclo + o -> clo + o2 
      cl057(l) = 2.4e-12*exp(-960./t(l)) 

c     cl058 : oclo + oh -> hocl + o2 
      cl058(l) = 1.4e-12*exp(600./t(l)) 

c     cl059: oclo + cl -> clo + clo 
      cl059(l) = 3.4e-11*exp(160./t(l))

      end do
     
      
     

      

      if (tribo .eq. 1.) then
c
c        electrochemical reactions
c
c        efmax: maximum electric field (kv.m-1)
c
         efmax = 23.3
c
c        ef: actual electric field, scaled by tau.
c
c        if (tau .ge. 1.) then
c           ef = efmax
c        else
c           ef = 0.
c        end if
c        ef = min(efmax,efmax*tau/1.0)
c
         ef = (efmax/0.5)*tau - (efmax/0.5)*0.5
c
         ef = max(ef, 0.)
         ef = min(ef, efmax)
c
ccc      t001: h2o + e -> oh + h-
c
c        lossh2o: fit of oh/h- production rates
c        given by delory et al., astrobiology, 6, 451, 2006
c
         if (ef .eq. 0.) then
            lossh2o = 0.
         else if (ef .lt. 10.) then
            lossh2o = 0.054136*exp(1.0978*ef)
         else if (ef .lt. 16.) then
            lossh2o = 64.85*exp(0.38894*ef)
         else if (ef .le. 20.) then
            lossh2o = 0.2466*exp(0.73719*ef)
         else
            lossh2o = 2.3269e-8*exp(1.546*ef)
         end if
c
c        production rates are given for h2o = 20 prec. microns.
c        t001 is converted to first-order reaction rate
c        assuming h2o number density at the surface =  5e13 mol cm-3
c
         do l = 1,21                     ! 70 km
            t001(l) = lossh2o/5.e13      ! s-1
         end do
         do l = 22,lswitch-1
            t001(l) = 0.
         end do
c
ccc      t002: ch4 + e -> products
c
c        lossch4: fit of ch4 loss rates
c        given by farrell et al., grl, 33, 2006
c
         if (ef .eq. 0.) then
            lossch4 = 0.
         else if (ef .gt. 20.) then
            lossch4 = 1.113e-21*exp(1.6065*ef)
         else if (ef .gt. 17.5) then
            lossch4 = 1.e-15*exp(0.92103*ef)
         else if (ef .gt. 14.) then
            lossch4 = 1.e-13*exp(0.65788*ef)
         else
            lossch4 = 8.9238e-15*exp(0.835*ef)
         end if
c
         do l = 1,21               ! 70 km
            t002(l) = lossch4      ! s-1
         end do
         do l = 22,lswitch-1
            t002(l) = 0.
         end do
c
ccc      t003: co2 + e -> co + o-
c
c        lossco2: fit of co/o- production rates
c        given by delory et al., astrobiology, 6, 451, 2006
c
         if (ef .eq. 0.) then
            lossco2 = 0.
         else if (ef .lt. 10.) then
            lossco2 = 22.437*exp(1.045*ef)
         else if (ef .lt. 16.) then
            lossco2 = 17518.*exp(0.37896*ef)
         else if (ef .lt. 20.) then
            lossco2 = 54.765*exp(0.73946*ef)
         else
            lossco2 = 4.911e-6*exp(1.5508*ef)
         end if
c
c        production rates are assumed to be given for p = 6 hPa
c        lossco2 is converted to first-order reaction rate
c        assuming co2 number density at the surface =  2e17 mol cm-3
c
         do l = 1,21                     ! 70 km
            t003(l) = lossco2/2.e17      ! s-1
         end do
         do l = 22,lswitch-1
            t003(l) = 0.
         end do
      else
         do l = 1,lswitch-1
            t001(l) = 0.
            t002(l) = 0.
            t003(l) = 0.
         end do
      end if
c

      
      return


      end

ccccccccccccccccccccccccccccccccccccc
cc 2-body equivalent rate coefficient
cc for 3 body reactions 
ccccccccccccccccccccccccccccccccccccc
      REAL FUNCTION k_3body(t_atm, M_dens, k0_300, n_3bod,
     $              kinf_300, m_3bod, fc_3bod, iupac)
     
          IMPLICIT NONE 
          
          REAL t_atm, M_dens, k0_300, n_3bod, kinf_300, m_3bod, fc_3bod 
          LOGICAL iupac
          REAL ratio, k_holding 

          REAL k0, kinf, N_cap 
          
          k0 = k0_300*(300./t_atm)**n_3bod
          kinf = kinf_300*(300./t_atm)**m_3bod
          ratio = k0*2.5*M_dens/kinf 
          
          IF ( iupac ) THEN 
               N_cap = 0.75 - 1.27*log10(fc_3bod)
               k_holding = 1./( 1. + ( log10(ratio)/N_cap)**2.)
               k_3body = (2.5*M_dens*k0/(1. + ratio))*fc_3bod**k_holding
          ELSE 
               k_holding = 1./( 1. + ( log10(ratio))**2.)
               k_3body = (2.5*M_dens*k0/(1. + ratio))*fc_3bod**k_holding
          ENDIF 
         
          RETURN 
      END 

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC Quasi-Steady-State Approximation 
cc for Cl- species in the 1-D model 
cccccccccccccccccccccccccccccccccccccc
      REAL FUNCTION Chlorine_QSSA(cc0,P,L,h)


      IMPLICIT NONE 

      REAL cc0 ! Initial number density 
      REAL P, L ! Production and Loss Rates
      REAL h ! Timestep 

      REAL cc ! Final number density

      IF ( h*L < 0.01 ) THEN 
!     1.0 : SIBEM 
        Chlorine_QSSA = (cc0 + P*h)/(1. + L*h)

      ELSEIF ( h*L .le. 10. ) THEN 
!     2.0 : qssa
        Chlorine_QSSA = cc0*exp(-h*L) + (P/L)*(1. - exp(-h*L) )

      ELSE
!     3.0 : Steady-State
        Chlorine_QSSA = P/L 
      ENDIF 


      RETURN
      END 