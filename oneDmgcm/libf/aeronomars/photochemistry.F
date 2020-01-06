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
     $                  surfice1d, jo3, tau)
c
      implicit none
c
#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h" 
#include "callkeys.h"
#include "tracer.h"
c
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
c
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
      parameter (nesp = 16)      ! number of species in the chemistry
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
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     stimestep  : standard timestep for the chemistry (s)       c
c     ctimestep  : real timestep for the chemistry (s)           c
c     phychemrat : ratio physical/chemical timestep              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      stimestep = 600. ! standard value : 10 mn
c
      phychemrat = nint(ptimestep/stimestep)
c
      ctimestep = ptimestep/real(phychemrat)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialisation of chemical species and families            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call gcmtochim(zycol, lswitch, nesp, rm) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute reaction rates                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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
     $               t001, t002, t003, tau)
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
c
      j_o3_o1d = 5
c
      do lev = 1,lswitch-1
         jo3(lev) = j(lev,j_o3_o1d)
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over time within the physical timestep                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do istep = 1,phychemrat
         call chimie(lswitch,nesp, rm, j, dens, ctimestep, 
     $               press, temp, sza, dist_sol,
     $               a001, a002, a003,
     $               b001, b002, b003, b004, b005, b006,
     $               b007, b008, b009,
     $               c001, c002, c003, c004, c005, c006,
     $               c007, c008, c009, c010, c011, c012,
     $               c013, c014, c015, c016, c017, c018,
     $               d001, d002, d003,
     $               e001, e002, e003,
     $               h001, h002, h003, h004, h005,
     $               t001, t002, t003)
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     save chemical species for the gcm                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call chimtogcm(zycol, lswitch, nesp, rm) 
c
      return
      end
c
c*****************************************************************
c
      subroutine chimie(lswitch, nesp, rm, j, dens, dt, 
     $                  press, t, sza, dist_sol,
     $                  a001, a002, a003,
     $                  b001, b002, b003, b004, b005, b006,
     $                  b007, b008, b009,
     $                  c001, c002, c003, c004, c005, c006,
     $                  c007, c008, c009, c010, c011, c012,
     $                  c013, c014, c015, c016, c017, c018,
     $                  d001, d002, d003,
     $                  e001, e002, e003,
     $                  h001, h002, h003, h004, h005,
     $                  t001, t002, t003)
c
c*****************************************************************
c
      implicit none 
c
#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h" 
#include "callkeys.h" 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     input/output:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      integer lswitch         ! interface level between chemistries
      integer nesp            ! number of species
c
      real rm(nlayermx,nesp)  ! volume mixing ratios 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inputs:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real dens(nlayermx)     ! density (cm-3)
      real dt                 ! chemistry timestep (s) 
      real j(nlayermx,nd)     ! interpolated photolysis rates (s-1) 
      real press(nlayermx)    ! pressure (hpa)
      real t(nlayermx)        ! temperature (k)
      real sza                ! solar zenith angle (deg)
      real dist_sol           ! sun distance (au) 
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
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real hetero_ice, hetero_dust
c
      integer iesp, iter, l, niter
      integer i_co2, i_co, i_o2, i_h2, i_h2o, i_h2o2, i_hox, i_ox,
     $        i_o1d, i_o, i_o3, i_h, i_oh, i_ho2, i_ch4, i_n2
      integer j_o2_o, j_o2_o1d, j_co2_o, j_co2_o1d,
     $        j_o3_o1d, j_o3_o, j_h2o, j_hdo, j_h2o2,
     $        j_ho2, j_no2, j_ch4_ch3_h, j_ch4_1ch2_h2,
     $        j_ch4_3ch2_h_h, j_ch4_ch_h2_h, j_ch3o2h,
     $        j_ch2o_co, j_ch2o_hco, j_ch3oh, j_c2h6, j_hcl,
     $        j_hocl, j_clo, j_so2, j_so, j_h2s, j_so3,
     $        j_hno3, j_hno4 
c
      parameter (hetero_ice  = 1.)   ! switch for het. chem. on ice clouds
      parameter (hetero_dust = 0.)   ! switch for het. chem. on dust
                                     ! hetero_dust = 0. advised for the moment
c
      parameter (niter = 5)          ! iterations in the chemical scheme
c
      real cc0(nlayermx,nesp)        ! initial number density (cm-3)
      real cc(nlayermx,nesp)         ! final number density (cm-3)
      real nox(nlayermx)             ! nox number density (cm-3)
      real no(nlayermx)              ! no  number density (cm-3)
      real no2(nlayermx)             ! no2 number density (cm-3)
      real production(nlayermx,nesp) ! production rate
      real loss(nlayermx,nesp)       ! loss rate
c
      real ro_o3, rh_ho2, roh_ho2, rno2_no
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      i_co2    =  1
      i_co     =  2
      i_o      =  3
      i_o1d    =  4
      i_o2     =  5
      i_o3     =  6
      i_h      =  7
      i_h2     =  8
      i_oh     =  9
      i_ho2    = 10 
      i_h2o2   = 11
      i_ch4    = 12
      i_h2o    = 13
      i_n2     = 14
      i_hox    = 15
      i_ox     = 16
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
      j_h2o2         =  9      ! h2o2 + hv   -> oh + oh
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
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     volume mixing ratio -> density conversion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iesp = 1,nesp
         do l = 1,lswitch-1
            cc0(l,iesp) = rm(l,iesp)*dens(l)
            cc(l,iesp)  = cc0(l,iesp)
         end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c     loop over iterations in the chemical scheme
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iter = 1,niter
c
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
c
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        hox species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        photochemical equilibrium for oh and ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h/ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
c
         rh_ho2 = (c001(l)*cc(l,i_o)
     $           + c004(l)*cc(l,i_h)
     $           + c005(l)*cc(l,i_h)
     $           + c006(l)*cc(l,i_h)
     $           + c007(l)*cc(l,i_oh)
     $           + 2.*c008(l)*cc(l,i_ho2)
     $           + c015(l)*cc(l,i_o3)
     $           + 2.*c016(l)*cc(l,i_ho2)
     $           + d003(l)*no(l)             ! ajout 20090401
     $           + j(l,j_ho2)
     $           + h001(l)*hetero_ice 
     $           + h003(l)*hetero_dust)
     $          /( c011(l)*cc(l,i_o2)
     $           + t001(l)*cc(l,i_h2o)
     $            /max(cc(l,i_h),dens(l)*1.e-30)   ! ajout 20090401
     $           )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        oh/ho2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         roh_ho2 = (c001(l)*cc(l,i_o)
     $            + c003(l)*cc(l,i_o3)*rh_ho2
     $            + 2.*c004(l)*cc(l,i_h)
     $            + 2.*c008(l)*cc(l,i_ho2)
     $            + c015(l)*cc(l,i_o3)
     $            + d003(l)*no(l)
     $            + j(l,j_ho2)
     $            + 2.*b002(l)*cc(l,i_o1d)*cc(l,i_h2o)         ! ajout 20101210
     $             /max(cc(l,i_ho2),dens(l)*1.e-30)            ! ajout 20101210
     $            + b003(l)*cc(l,i_o1d)*cc(l,i_h2)             ! ajout 20101210
     $             /max(cc(l,i_ho2),dens(l)*1.e-30)            ! ajout 20101210
     $            + j(l,j_h2o)*cc(l,i_h2o)
     $             /max(cc(l,i_ho2),dens(l)*1.e-30)
     $            + t001(l)*cc(l,i_h2o)               ! suppression 20090401
     $             /max(cc(l,i_ho2),dens(l)*1.e-30)   ! suppression 20090401
     $            )
     $            /(c002(l)*cc(l,i_o)
     $            + c007(l)*cc(l,i_ho2)
     $            + c009(l)*cc(l,i_h2o2)         ! ajout 20090401
     $            + 2.*c013(l)*cc(l,i_oh)        ! ajout 20090401
     $            + 2.*c017(l)*cc(l,i_oh)        ! ajout 20090401
     $            + e001(l)*cc(l,i_co)
     $            + h002(l)*hetero_ice)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         cc(l,i_h) = cc(l,i_hox)
     $                 /(1. + (1. + roh_ho2)/rh_ho2)
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
      end do
c
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
     $                   + b006(l)*cc(l,i_o3))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o/o3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            ro_o3 = (j(l,j_o3_o1d) + j(l,j_o3_o)
     $              + a003(l)*cc(l,i_o)
     $              + c003(l)*cc(l,i_h)
     $              + c014(l)*cc(l,i_oh)
     $              + c015(l)*cc(l,i_ho2)
     $              )
     $              /(a001(l)*cc(l,i_o2))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            cc(l,i_o3) = cc(l,i_ox)/(1. + ro_o3)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            cc(l,i_o) = cc(l,i_o3)*ro_o3
c
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
c
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
c
            loss(l,i_ox) = loss(l,i_ox)/cc(l,i_ox)
c
         end do
c
      else
c
         do l = 1,lswitch-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o(1d)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            cc(l,i_o1d) = 0.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        o3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            production(l,i_o3) = a001(l)*cc(l,i_o2)*cc(l,i_o)
c
            loss(l,i_o3) = a003(l)*cc(l,i_o)
     $                   + c003(l)*cc(l,i_h)
     $                   + c014(l)*cc(l,i_oh)
     $                   + c015(l)*cc(l,i_ho2)
c

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
c
         end do
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     other species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do l = 1,lswitch-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        co2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_co2) = e001(l)*cc(l,i_oh)*cc(l,i_co)
     $                       + e002(l)*cc(l,i_o)*cc(l,i_co)
     $                       + t002(l)*cc(l,i_ch4)*16./44. ! ajout 20090401
     $                       + t003(l)*cc(l,i_co2)*8./44.  ! ajout 20090401
c
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
c
         loss(l,i_co) = e001(l)*cc(l,i_oh)
     $                + e002(l)*cc(l,i_o)
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
c
         loss(l,i_o2) = j(l,j_o2_o)
     $                + j(l,j_o2_o1d)
     $                + a001(l)*cc(l,i_o)
     $                + c011(l)*cc(l,i_h)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        h2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_h2) = c005(l)*cc(l,i_h)*cc(l,i_ho2)
     $                      + c018(l)*cc(l,i_h)*cc(l,i_h)
c
         loss(l,i_h2) = b003(l)*cc(l,i_o1d)
     $                + c010(l)*cc(l,i_oh)
c
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
c
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
c
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
c
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
c
         loss(l,i_hox) = loss(l,i_hox)/cc(l,i_hox)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ch4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         production(l,i_ch4) = 0.
c
         loss(l,i_ch4) = j(l,j_ch4_ch3_h)
     $                 + j(l,j_ch4_1ch2_h2)
     $                 + j(l,j_ch4_3ch2_h_h)
     $                 + j(l,j_ch4_ch_h2_h)
     $                 + b007(l)*cc(l,i_o1d) 
     $                 + b008(l)*cc(l,i_o1d)
     $                 + b009(l)*cc(l,i_o1d)
     $                 + e003(l)*cc(l,i_oh)
c
      end do
c
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
         cc(l,i_ch4) = (cc0(l,i_ch4) + production(l,i_ch4)*dt)
     $              /(1. + loss(l,i_ch4)*dt)
      end do
c
c     ox species
c
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
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end of loop over iterations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     density -> volume mixing ratio conversion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iesp = 1,nesp
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end of loop over vertical levels
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      end do
c
      else
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     night
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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
c
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

      integer,parameter :: i_co2  =  1
      integer,parameter :: i_co   =  2
      integer,parameter :: i_o    =  3
      integer,parameter :: i_o1d  =  4
      integer,parameter :: i_o2   =  5
      integer,parameter :: i_o3   =  6
      integer,parameter :: i_h    =  7
      integer,parameter :: i_h2   =  8
      integer,parameter :: i_oh   =  9
      integer,parameter :: i_ho2  = 10
      integer,parameter :: i_h2o2 = 11
      integer,parameter :: i_ch4  = 12
      integer,parameter :: i_h2o  = 13
      integer,parameter :: i_n2   = 14
      integer,parameter :: i_hox  = 15
      integer,parameter :: i_ox   = 16
c      
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
c
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

      integer,parameter :: i_co2  =  1
      integer,parameter :: i_co   =  2
      integer,parameter :: i_o    =  3
      integer,parameter :: i_o1d  =  4
      integer,parameter :: i_o2   =  5
      integer,parameter :: i_o3   =  6
      integer,parameter :: i_h    =  7
      integer,parameter :: i_h2   =  8
      integer,parameter :: i_oh   =  9
      integer,parameter :: i_ho2  = 10
      integer,parameter :: i_h2o2 = 11
      integer,parameter :: i_ch4  = 12
      integer,parameter :: i_h2o  = 13
      integer,parameter :: i_n2   = 14
      integer,parameter :: i_hox  = 15
      integer,parameter :: i_ox   = 16
c      
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
c
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
     $                     t001, t002, t003, tau)
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
c
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute reaction rates
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c
         a001(l) = 2.075
     $             *6.0e-34*(t(l)/300.)**(-2.4)*dens(l)
c
c        mulcahy and williams, 1968
c
c        a001(l) = 2.68e-33*(t(l)/298.)**(-2.4)*dens(l)
c
c        nair et al., 1994
c
c        a001(l) = 1.3e-34*exp(724./t(l))*dens(l)
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
ccc      b007: o(1d) + ch4 -> ch3 + oh
c
c        jpl 2003
c
         b007(l) = 1.5e-10*0.75
c
ccc      b008: o(1d) + ch4 -> ch3o + h
c
c        jpl 2003
c
         b008(l) = 1.5e-10*0.20
c
ccc      b009: o(1d) + ch4 -> ch2o + h2
c
c        jpl 2003
c
         b009(l) = 1.5e-10*0.05
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
c
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
         e003(l) = 2.80e-14*(t(l)**0.667)*exp(-1575./t(l))
c
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
c
      end do
c
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
