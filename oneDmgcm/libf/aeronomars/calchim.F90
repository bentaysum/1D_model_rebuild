      subroutine calchim(ptimestep,pplay,pplev,pt,pdt,dist_sol,mu0,         &
                         zzlev,zzlay,zday,pq,pdq,dqchim,dqschim,dqcloud,    &
                         dqscloud,tauref,co2ice,                            &
                         pu,pdu,pv,pdv,surfdust,surfice)

      implicit none

!=======================================================================
!
!   subject:
!   --------
!
!  Prepare the call for the photochemical module, and send back the
!  tendencies from photochemistry in the chemical species mass mixing ratios
!
!   Author:   Sebastien Lebonnois (08/11/2002)
!   -------
!    update 12/06/2003 for water ice clouds and compatibility with dust
!    update 07/2003 for coupling with thermosphere (Monica Angelats-i-Coll)
!    update 03/05/2005 cosmetic changes (Franck Lefevre)
!    update sept. 2008 identify tracers by their names (Ehouarn Millour)
!    update 05/12/2011 synchronize with latest version of chemistry (Franck Lefevre)
!    update 16/03/2012 optimization (Franck Lefevre)
!
!   Arguments:
!   ----------
!
!  Input:
!
!    ptimestep                  timestep (s)
!    pplay(ngridmx,nlayermx)    Pressure at the middle of the layers (Pa)
!    pplev(ngridmx,nlayermx+1)  Intermediate pressure levels (Pa)
!    pt(ngridmx,nlayermx)       Temperature (K)
!    pdt(ngridmx,nlayermx)      Temperature tendency (K)
!    pu(ngridmx,nlayermx)       u component of the wind (ms-1)
!    pdu(ngridmx,nlayermx)      u component tendency (K)
!    pv(ngridmx,nlayermx)       v component of the wind (ms-1)
!    pdv(ngridmx,nlayermx)      v component tendency (K)
!    dist_sol                   distance of the sun (AU)
!    mu0(ngridmx)               cos of solar zenith angle (=1 when sun at zenith)
!    pq(ngridmx,nlayermx,nqmx)  Advected fields, ie chemical species here
!    pdq(ngridmx,nlayermx,nqmx) Previous tendencies on pq
!    tauref(ngridmx)            Optical depth at 7 hPa
!    co2ice(ngridmx)            co2 ice surface layer (kg.m-2)
!    surfdust(ngridmx,nlayermx) dust surface area (m2/m3)
!    surfice(ngridmx,nlayermx)  ice surface area (m2/m3)
!
!  Output:
!
!    dqchim(ngridmx,nlayermx,nqmx) ! tendencies on pq due to chemistry
!    dqschim(ngridmx,nqmx)         ! tendencies on qsurf 
!
!=======================================================================

#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "conc.h"

!     input:

      real :: ptimestep
      real :: pplay(ngridmx,nlayermx)    ! pressure at the middle of the layers
      real :: zzlay(ngridmx,nlayermx)    ! pressure at the middle of the layers
      real :: pplev(ngridmx,nlayermx+1)  ! intermediate pressure levels
      real :: zzlev(ngridmx,nlayermx+1)  ! altitude at layer boundaries
      real :: pt(ngridmx,nlayermx)       ! temperature
      real :: pdt(ngridmx,nlayermx)      ! temperature tendency
      real :: pu(ngridmx,nlayermx)       ! u component of the wind (m.s-1)
      real :: pdu(ngridmx,nlayermx)      ! u component tendency
      real :: pv(ngridmx,nlayermx)       ! v component of the wind (m.s-1)
      real :: pdv(ngridmx,nlayermx)      ! v component tendency
      real :: dist_sol                   ! distance of the sun (AU)
      real :: mu0(ngridmx)               ! cos of solar zenith angle (=1 when sun at zenith)
      real :: pq(ngridmx,nlayermx,nqmx)  ! tracers mass mixing ratio
      real :: pdq(ngridmx,nlayermx,nqmx) ! previous tendencies
      real :: zday                       ! date (time since Ls=0, in martian days)
      real :: tauref(ngridmx)            ! optical depth at 7 hPa
      real :: co2ice(ngridmx)            ! co2 ice surface layer (kg.m-2)
      real :: surfdust(ngridmx,nlayermx) ! dust surface area (m2/m3)
      real :: surfice(ngridmx,nlayermx)  !  ice surface area (m2/m3)

!     output:

      real :: dqchim(ngridmx,nlayermx,nqmx) ! tendencies on pq due to chemistry
      real :: dqschim(ngridmx,nqmx)         ! tendencies on qsurf 
      real :: dqcloud(ngridmx,nlayermx,nqmx)! tendencies on pq due to condensation
      real :: dqscloud(ngridmx,nqmx)        ! tendencies on qsurf 

!     local variables:

      integer,save :: nbq                   ! number of tracers used in the chemistry
      integer,save :: niq(nqmx)             ! array storing the indexes of the tracers
      integer :: iloc(1)            ! index of major species
      integer :: ig,l,i,iq,iqmax
      integer :: foundswitch, lswitch
      integer,save :: chemthermod

      integer,save :: i_co2  = 0
      integer,save :: i_co   = 0
      integer,save :: i_o    = 0
      integer,save :: i_o1d  = 0
      integer,save :: i_o2   = 0
      integer,save :: i_o3   = 0
      integer,save :: i_h    = 0
      integer,save :: i_h2   = 0
      integer,save :: i_oh   = 0
      integer,save :: i_ho2  = 0
      integer,save :: i_h2o2 = 0
      integer,save :: i_ch4  = 0
      integer,save :: i_n2   = 0
      integer,save :: i_h2o  = 0
      integer,save :: i_n    = 0
      integer,save :: i_no   = 0
      integer,save :: i_no2  = 0
      integer,save :: i_n2d  = 0
      integer,save :: i_co2plus=0
      integer,save :: i_oplus=0
      integer,save :: i_o2plus=0
      integer,save :: i_coplus=0
      integer,save :: i_cplus=0
      integer,save :: i_nplus=0
      integer,save :: i_noplus=0
      integer,save :: i_n2plus=0
      integer,save :: i_hplus=0
      integer,save :: i_hco2plus=0
      integer,save :: i_elec=0

      integer :: ig_vl1

      real    :: latvl1, lonvl1
      real    :: zq(ngridmx,nlayermx,nqmx) ! pq+pdq*ptimestep before chemistry
                                           ! new mole fraction after
      real    :: zt(ngridmx,nlayermx)      ! temperature
      real    :: zu(ngridmx,nlayermx)      ! u component of the wind
      real    :: zv(ngridmx,nlayermx)      ! v component of the wind
      real    :: taucol                    ! optical depth at 7 hPa

      logical,save :: firstcall = .true.
      logical,save :: depos = .false.      ! switch for dry deposition

!     for each column of atmosphere:

      real :: zpress(nlayermx)       !  Pressure (mbar)
      real :: zdens(nlayermx)        !  Density  (cm-3)
      real :: ztemp(nlayermx)        !  Temperature (K)
      real :: zlocal(nlayermx)       !  Altitude (km)
      real :: zycol(nlayermx,nqmx)   !  Composition (mole fractions)
      real :: szacol                 !  Solar zenith angle
      real :: surfice1d(nlayermx)    !  Ice surface area (cm2/cm3)
      real :: surfdust1d(nlayermx)   !  Dust surface area (cm2/cm3)
      real :: jo3(nlayermx)          !  Photodissociation rate O3->O1D (s-1)

!     for output:

      logical :: output                 ! to issue calls to writediagfi and stats
      parameter (output = .true.)
      real :: jo3_3d(ngridmx,nlayermx)  ! Photodissociation rate O3->O1D (s-1)

!=======================================================================
!     initialization of the chemistry (first call only)
!=======================================================================

      if (firstcall) then

         if (photochem) then
            print*,'calchim: Read photolysis lookup table'
            call read_phototable
         end if
         ! find index of chemical tracers to use
         ! Listed here are all tracers that can go into photochemistry
         nbq = 0 ! to count number of tracers
         ! Species ALWAYS present if photochem=.T. or thermochem=.T.
         i_co2 = igcm_co2
         if (i_co2 == 0) then
            write(*,*) "calchim: Error; no CO2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_co2
         end if
         i_co = igcm_co
         if (i_co == 0) then
            write(*,*) "calchim: Error; no CO tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_co
         end if
         i_o = igcm_o
         if (i_o == 0) then
            write(*,*) "calchim: Error; no O tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o
         end if
         i_o1d = igcm_o1d
         if (i_o1d == 0) then
            write(*,*) "calchim: Error; no O1D tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o1d
         end if
         i_o2 = igcm_o2
         if (i_o2 == 0) then
            write(*,*) "calchim: Error; no O2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o2
         end if
         i_o3 = igcm_o3
         if (i_o3 == 0) then
            write(*,*) "calchim: Error; no O3 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o3
         end if
         i_h = igcm_h
         if (i_h == 0) then
            write(*,*) "calchim: Error; no H tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h
         end if
         i_h2 = igcm_h2
         if (i_h2 == 0) then
            write(*,*) "calchim: Error; no H2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h2
         end if
         i_oh = igcm_oh
         if (i_oh == 0) then
            write(*,*) "calchim: Error; no OH tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_oh
         end if
         i_ho2 = igcm_ho2
         if (i_ho2 == 0) then
            write(*,*) "calchim: Error; no HO2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_ho2
         end if
         i_h2o2 = igcm_h2o2
         if (i_h2o2 == 0) then
            write(*,*) "calchim: Error; no H2O2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h2o2
         end if
         i_ch4 = igcm_ch4
         if (i_ch4 == 0) then
            write(*,*) "calchim: Error; no CH4 tracer !!!"
            write(*,*) "CH4 will be ignored in the chemistry"
         else
            nbq = nbq + 1
            niq(nbq) = i_ch4
         end if
         i_n2 = igcm_n2
         if (i_n2 == 0) then
            write(*,*) "calchim: Error; no N2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_n2
         end if
         i_h2o = igcm_h2o_vap
         if (i_h2o == 0) then
            write(*,*) "calchim: Error; no water vapor tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h2o
         end if
         !Check tracers needed for thermospheric chemistry
         if(thermochem) then
            chemthermod=0  !Default: C/O/H chemistry
            !Nitrogen chemistry
            !NO is used to determine if N chemistry is wanted
            !chemthermod=2 -> N chemistry
            i_no = igcm_no
            if (i_no == 0) then
               write(*,*) "calchim: no NO tracer"
               write(*,*) "C/O/H themosp chemistry only "
            else
               nbq = nbq + 1
               niq(nbq) = i_no
               chemthermod=2
               write(*,*) "calchim: NO in traceur.def"
               write(*,*) "Nitrogen chemistry included"
            end if
            ! N
            i_n = igcm_n
            if(chemthermod == 2) then
               if (i_n == 0) then
                  write(*,*) "calchim: Error; no N tracer !!!"
                  write(*,*) "N is needed if NO is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_n
               end if
            else 
               if (i_n /= 0) then
                  write(*,*) "calchim: Error: N is present, but NO is not!!!"
                  write(*,*) "Both must be in traceur.def if N chemistry is wanted"
                  stop
               endif
            endif    !Of if(chemthermod == 2) 
            ! NO2
            i_no2 = igcm_no2
            if(chemthermod == 2) then
               if (i_no2 == 0) then
                  write(*,*) "calchim: Error; no NO2 tracer !!!"
                  write(*,*) "NO2 is needed if NO is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_no2
               end if
            else
               if (i_no2 /= 0) then
                  write(*,*) "calchim: Error: N is present, but NO is not!!!"
                  write(*,*) "Both must be in traceur.def if N chemistry is wanted"
                  stop
               endif
            endif     !Of if(chemthermod == 2)
            ! N(2D)
            if(chemthermod == 2) then
               i_n2d = igcm_n2d
               if (i_n2d == 0) then
                  write(*,*) "calchim: Error; no N2D !!!"
                  write(*,*) "N2D is needed if NO is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_n2d
               end if
            else
               if (i_n2d /= 0) then
                  write(*,*) "calchim: Error: N2D is present, but NO is not!!!"
                  write(*,*) "Both must be in traceur.def if N chemistry wanted"
                  stop
               endif
            endif    !Of if(chemthermod == 2)
            ! Ions
            ! O2+ is used to determine if ion chemistry is needed
            ! chemthermod=3 -> ion chemistry
            i_o2plus = igcm_o2plus
            if(chemthermod == 2) then
               if (i_o2plus == 0) then
                  write(*,*) "calchim: no O2+ tracer; no ion chemistry"
               else
                  nbq = nbq + 1
                  niq(nbq) = i_o2plus
                  chemthermod = 3
                  write(*,*) "calchim: O2+ in traceur.def"
                  write(*,*) "Ion chemistry included"
               end if
            else
               if (i_o2plus /= 0) then
                  write(*,*) "calchim: O2+ is present, but NO is not!!!"
                  write(*,*) "Both must be in traceur.def if ionosphere wanted"
                  stop
               endif
            endif   !Of if(chemthermod == 2)
            ! CO2+
            i_co2plus = igcm_co2plus
            if(chemthermod == 3) then
               if (i_co2plus == 0) then
                  write(*,*) "calchim: Error; no CO2+ tracer !!!"
                  write(*,*) "CO2+ is needed if O2+ is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_co2plus
               end if
            else
               if (i_co2plus /= 0) then
                  write(*,*) "calchim: Error: CO2+ is present, but O2+ is not!!!"
                  write(*,*) "Both must be in traceur.def if ionosphere wanted"
                  stop
               endif
            endif    !Of if(chemthermod == 3)
            ! O+
            i_oplus = igcm_oplus
            if(chemthermod == 3) then
               if (i_oplus == 0) then
                  write(*,*) "calchim: Error; no O+ tracer !!!"
                  write(*,*) "O+ is needed if O2+ is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_oplus
               end if
            else
               if (i_oplus /= 0) then
                  write(*,*) "calchim: Error: O+ is present, but O2+ is not!!!"
                  write(*,*) "Both must be in traceur.def if ionosphere wanted"
                  stop
               endif
            endif   !Of if (chemthermod == 3)
            ! CO+
            i_coplus = igcm_coplus
            if(chemthermod == 3) then
               if (i_coplus == 0) then
                  write(*,*) "calchim: Error; no CO+ tracer !!!"
                  write(*,*) "CO+ is needed if O2+ is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_coplus
               end if
            else
               if (i_coplus /= 0) then
                  write(*,*) "calchim: Error: CO+ is present, but O2+ is not!!!"
                  write(*,*) " Both must be in traceur.def if ionosphere wanted"
                  stop
               endif
            endif   ! Of if (chemthermod == 3)
            ! C+
            i_cplus = igcm_cplus
            if(chemthermod == 3) then
               if (i_cplus == 0) then
                  write(*,*) "calchim: Error; no C+ tracer !!!"
                  write(*,*) "C+ is needed if O2+ is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_cplus
               end if
            else
               if (i_cplus /= 0) then
                  write(*,*) "calchim: Error; C+ is present, but O2+ is not!!!"
                  write(*,*) "Both must be in traceur.def if ionosphere wanted"
                  stop
               endif
            endif   ! Of if (chemthermod == 3)
            ! N+
            i_nplus = igcm_nplus
            if(chemthermod == 3) then
               if (i_nplus == 0) then
                  write(*,*) "calchim: Error; no N+ tracer !!!"
                  write(*,*) "N+ is needed if O2+ is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_nplus
               end if
            else
               if (i_nplus /= 0) then
                  write(*,*) "calchim: Error: N+ is present, but O2+ is not!!!"
                  write(*,*) "Both must be in traceur.def if ionosphere wanted"
                  stop
               endif
            endif   !Of if (chemthermod == 3)
            ! NO+
            i_noplus = igcm_noplus
            if(chemthermod == 3) then
               if (i_noplus == 0) then
                  write(*,*) "calchim: Error; no NO+ tracer !!!"
                  write(*,*) "NO+ is needed if O2+ is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_noplus
               end if
            else
               if (i_noplus /= 0) then
                  write(*,*) "calchim: Error: NO+ is present, but O2+ is not!!!"
                  write(*,*) "Both must be in traceur.def if ionosphere wanted"
                  stop
               endif
            endif   !Of if (chemthermod == 3)
            ! N2+
            i_n2plus = igcm_n2plus
            if (chemthermod == 3) then
               if (i_n2plus == 0) then
                  write(*,*) "calchim: Error; no N2+ tracer !!!"
                  write(*,*) "N2+ is needed if O2+ is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_n2plus
               end if
            else
               if (i_n2plus /= 0) then
                  write(*,*) "calchim: Error: N2+ is present, but O2+ is not!!!"
                  write(*,*) "Both must be in traceur.def if ionosphere wanted"
                  stop
               endif
            endif   !Of if (chemthermod == 3)
            !H+
            i_hplus = igcm_hplus
            if (chemthermod == 3) then
               if (i_hplus == 0) then
                  write(*,*) "calchim: Error; no H+ tracer !!!"
                  write(*,*) "H+ is needed if O2+ is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_hplus
               end if
            else
               if (i_hplus /= 0) then
                  write(*,*) "calchim: Error: H+ is present, but O2+ is not!!!"
                  write(*,*) "Both must be in traceur.def if ionosphere wanted"
                  stop
               endif
            endif   !Of if (chemthermod == 3)
            ! HCO2+
            i_hco2plus = igcm_hco2plus
            if(chemthermod == 3) then
               if (i_hco2plus == 0) then
                  write(*,*) "calchim: Error; no HCO2+ tracer !!!"
                  write(*,*) "HCO2+ is needed if O2+ is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_hco2plus
               end if
            else
               if (i_hco2plus /= 0) then
                  write(*,*) "calchim: Error: HCO2+ is present, but O2+ is not!!!"
                  write(*,*) "Both must be in traceur.def if ionosphere wanted"
                  stop
               endif
            endif    !Of if(chemthermod == 3)
            !e-
            i_elec = igcm_elec
            if(chemthermod == 3) then
               if (i_elec == 0) then
                  write(*,*) "calchim: Error; no e- tracer !!!"
                  write(*,*) "e- is needed if O2+ is in traceur.def"
                  stop
               else
                  nbq = nbq + 1
                  niq(nbq) = i_elec
               end if
            else
               if(i_elec /= 0) then
                  write(*,*) "calchim: Error: e- is present, but O2+ is not!!!"
                  write(*,*) "Both must be in traceur.def if ionosphere wanted"
                  stop
               endif
            endif   !Of if (chemthermod == 3)
         endif      !Of thermochem

         write(*,*) 'calchim: found nbq    = ',nbq,' tracers'
               
         firstcall = .false.
      end if ! if (firstcall)

! Initializations

      zycol(:,:)    = 0.
      dqchim(:,:,:) = 0.
      dqschim(:,:)  = 0.

!     latvl1= 22.27
!     lonvl1= -47.94
!     ig_vl1= 1+ int( (1.5-(latvl1-90.)*jjm/180.)  -2 )*iim +    &
!             int(1.5+(lonvl1+180)*iim/360.)

!=======================================================================
!     loop over grid
!=======================================================================

      do ig = 1,ngridmx
         
         foundswitch = 0
         do l = 1,nlayermx
            do i = 1,nbq
               iq = niq(i) ! get tracer index
               zq(ig,l,iq) = pq(ig,l,iq) + pdq(ig,l,iq)*ptimestep
               zycol(l,iq) = zq(ig,l,iq)*mmean(ig,l)/mmol(iq)
            end do
            zt(ig,l)  = pt(ig,l) + pdt(ig,l)*ptimestep
            zu(ig,l)  = pu(ig,l) + pdu(ig,l)*ptimestep
            zv(ig,l)  = pv(ig,l) + pdv(ig,l)*ptimestep
            zpress(l) = pplay(ig,l)/100.
            ztemp(l)  = zt(ig,l)
            zdens(l)  = zpress(l)/(kb*1.e4*ztemp(l))
            zlocal(l) = zzlay(ig,l)/1000.

!           surfdust1d and surfice1d: conversion from m2/m3 to cm2/cm3

            surfdust1d(l) = surfdust(ig,l)*1.e-2
            surfice1d(l)  = surfice(ig,l)*1.e-2

!           search for switch index between regions

            if (photochem .and. thermochem) then
               if (foundswitch == 0 .and. pplay(ig,l) < 1.e-1) then
                  lswitch = l
                  foundswitch = 1
               end if
            end if
            if (.not. photochem) then
               lswitch = 22
            end if
            if (.not. thermochem) then
               lswitch = min(50,nlayermx+1)
            end if

         end do ! of do l=1,nlayermx

         szacol = acos(mu0(ig))*180./pi
         taucol = tauref(ig)*(700./610.)  ! provisoire en attente de nouveau jmars

!=======================================================================
!     call chemical subroutines
!=======================================================================

!        chemistry in lower atmosphere

         if (photochem) then
            call photochemistry(lswitch,zycol,szacol,ptimestep,    &
                                zpress,ztemp,zdens,dist_sol,       &
                                surfdust1d,surfice1d,jo3,taucol)

!        ozone photolysis, for output

            do l = 1,nlayermx
               jo3_3d(ig,l) = jo3(l)
            end do

!        condensation of h2o2

            call perosat(ig,ptimestep,pplev,pplay,                 &
                         ztemp,zycol,dqcloud,dqscloud)
         end if

!        chemistry in upper atmosphere
        
         if (thermochem) then
            call chemthermos(ig,lswitch,chemthermod,zycol,ztemp,zdens,  &
                             zpress,zlocal,szacol,ptimestep,zday)
         end if

!        dry deposition

         if (depos) then
            call deposition(ig, ig_vl1, pplay, pplev, zzlay, zzlev,& 
                            zu, zv, zt, zycol, ptimestep, co2ice)
         end if
!=======================================================================
!     tendencies
!=======================================================================

!     index of the most abundant species at each level

!         major(:) = maxloc(zycol, dim = 2)

!     tendency for the most abundant species = - sum of others
         do l = 1,nlayermx
            iloc=maxloc(zycol(l,:))
            iqmax=iloc(1)
            do i = 1,nbq
               iq = niq(i) ! get tracer index
               if (iq /= iqmax) then
                  dqchim(ig,l,iq) = (zycol(l,iq)*mmol(iq)/mmean(ig,l)  &
                                   - zq(ig,l,iq))/ptimestep
                  dqchim(ig,l,iqmax) = dqchim(ig,l,iqmax)              &
                                     - dqchim(ig,l,iq) 
               end if
            end do
         end do ! of do l = 1,nlayermx

!=======================================================================
!     end of loop over grid
!=======================================================================

      end do ! of do ig=1,ngridmx

!=======================================================================
!     write outputs
!=======================================================================

! value of parameter 'output' to trigger writting of outputs
! is set above at the declaration of the variable.

      if (photochem .and. output) then
         if (ngridmx > 1) then
            call writediagfi(ngridmx,'jo3','j o3->o1d',    &
                             's-1',3,jo3_3d(1,1))
           if (callstats) then
              call wstats(ngridmx,'jo3','j o3->o1d',       &
                          's-1',3,jo3_3d(1,1))
           endif
         end if ! of if (ngridmx.gt.1)
      end if ! of if (output)

      return
      end

