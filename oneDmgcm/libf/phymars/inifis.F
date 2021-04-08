      SUBROUTINE inifis(
     $           ngrid,nlayer
     $           ,day_ini,pdaysec,ptimestep
     $           ,plat,plon,parea
     $           ,prad,pg,pr,pcpp
#ifdef MESOSCALE
#include "meso_inc/meso_inc_inifisinvar.F"
#endif
     $           )
!
!=======================================================================
!
!   purpose:
!   -------
!
!   Initialisation for the physical parametrisations of the LMD 
!   martian atmospheric general circulation modele.
!
!   author: Frederic Hourdin 15 / 10 /93
!   -------
!   modified: Sebastien Lebonnois 11/06/2003 (new callphys.def)
!             Ehouarn Millour (oct. 2008) tracers are now identified
!              by their names and may not be contiguously
!              stored in the q(:,:,:,:) array
!             E.M. (june 2009) use getin routine to load parameters
!             adapted to the mesoscale use - Aymeric Spiga - 01/2007-07/2011
!
!
!   arguments:
!   ----------
!
!   input:
!   ------
!
!    ngrid                 Size of the horizontal grid.
!                          All internal loops are performed on that grid.
!    nlayer                Number of vertical layers.
!    pdayref               Day of reference for the simulation
!    pday                  Number of days counted from the North. Spring
!                          equinoxe.
!
!=======================================================================
!
!-----------------------------------------------------------------------
!   declarations:
!   -------------
! to use  'getin'
      USE ioipsl_getincom 
      IMPLICIT NONE
#include "dimensions.h"
#include "dimphys.h"
#include "planete.h"
#include "comcstfi.h"
#include "comsaison.h"
#include "comdiurn.h"
#include "comgeomfi.h"
#include "callkeys.h"
#include "surfdat.h"
#include "dimradmars.h"
#include "yomaer.h"
#include "datafile.h"
#include "slope.h"
#include "microphys.h"
#include "tracer.h"
#ifdef MESOSCALE
#include "comsoil.h"     !!MESOSCALE -- needed to fill volcapa
#include "meso_inc/meso_inc_inifisvar.F"
#endif
      REAL prad,pg,pr,pcpp,pdaysec

      REAL ptimestep 
      INTEGER day_ini

      INTEGER ngrid,nlayer
      REAL plat(ngrid),plon(ngrid),parea(ngridmx)
      INTEGER ig,ierr
 
!      EXTERNAL iniorbit,orbite
      EXTERNAL SSUM
      REAL SSUM
 
      CHARACTER ch1*12
      CHARACTER ch80*80
      CHARACTER tester*100 

!      logical chem, h2o

!      chem = .false.
!      h2o = .false.

      rad=prad
      cpp=pcpp
      g=pg
      r=pr
      rcp=r/cpp
      daysec=pdaysec
      dtphys=ptimestep
#ifdef MESOSCALE
#include "meso_inc/meso_inc_inifisini.F"
#endif

! --------------------------------------------------------
!     The usual Tests
!     --------------
      IF (nlayer.NE.nlayermx) THEN
         PRINT*,'STOP in inifis'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'nlayer     = ',nlayer
         PRINT*,'nlayermx   = ',nlayermx
         STOP
      ENDIF

      IF (ngrid.NE.ngridmx) THEN
         PRINT*,'STOP in inifis'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'ngrid     = ',ngrid
         PRINT*,'ngridmx   = ',ngridmx
         STOP
      ENDIF

! --------------------------------------------------------------
!  Reading the "callphys.def" file controlling some key options
! --------------------------------------------------------------
     
      ! check that 'callphys.def' file is around
      OPEN(99,file='callphys.def',status='old',form='formatted'
     &     ,iostat=ierr)
      CLOSE(99)
      
      IF(ierr.EQ.0) THEN
         PRINT*
         PRINT*
         PRINT*,'--------------------------------------------'
         PRINT*,' inifis: Parameters for the physics (callphys.def)'
         PRINT*,'--------------------------------------------'

         write(*,*) "Directory where external input files are:"
! PIP
         !datafile="/u/forget/WWW/datagcm/datafile"
         datafile="/Users/ppalmer/Desktop/oneDmodel/mgcm-datafile"
         call getin("datadir",datafile) ! default path
         write(*,*) " datafile = ",trim(datafile)

         write(*,*) "Run with or without tracer transport ?"
         tracer=.false. ! default value
         call getin("tracer",tracer)
         write(*,*) " tracer = ",tracer

         write(*,*) "Diurnal cycle ?"
         write(*,*) "(if diurnal=False, diurnal averaged solar heating)"
         diurnal=.true. ! default value
         call getin("diurnal",diurnal)
         write(*,*) " diurnal = ",diurnal

         write(*,*) "Seasonal cycle ?"
         write(*,*) "(if season=False, Ls stays constant, to value ",
     &   "set in 'start'"
         season=.true. ! default value
         call getin("season",season)
         write(*,*) " season = ",season

!        MCD Interpolation keys
         write(*,*) "MCD v5.3 tracer profiles of CO2, CO, O2, H2 and ",
     &              "H2O vapour?"
         mcd_profiles=.True. ! default value
         call getin("mcd_profiles",mcd_profiles)
         write(*,*) " mcd_profiles = ",mcd_profiles

         IF ( mcd_profiles ) THEN 
            mcd_co2 = .False. 
            mcd_co = .False.
            mcd_h2 = .False.
            mcd_h2o = .False.
            mcd_o2 = .False.

         call getin("mcd_co2",mcd_co2)
         call getin("mcd_co",mcd_co)
         call getin("mcd_h2",mcd_h2)
         call getin("mcd_h2o",mcd_h2o)
         call getin("mcd_o2",mcd_o2)

         write(*,*) "mcd_co2 = ", mcd_co2
         write(*,*) "mcd_co = ", mcd_co
         write(*,*) "mcd_h2 = ", mcd_h2
         write(*,*) "mcd_h2o = ", mcd_h2o
         write(*,*) "mcd_o2 = ", mcd_o2

         call getin("mcd_dir",mcd_dir)
         call getin("mcd_file",mcd_file)

         ENDIF 

         write(*,*) "Longitudinal mean values of MCD v5.3 variables?"
         long_mean=.true. ! default value
         call getin("long_mean",long_mean)
         write(*,*) " long_mean = ",long_mean

!       NOMAD/ACS Profiles Used?
        write(*,*) "Use NOMAD/ACS Tracer Profiles?"
        NOMAD_ACS = .false. ! default
        call getin("NOMAD_ACS", NOMAD_ACS)
        write(*,*)"NOMAD_ACS = ", NOMAD_ACS


         
!        Tangent Linear Model Parameters 
         call_tlm = .False. ! default
         call getin("call_tlm",call_tlm)
         write(*,*) "call_tlm : ", call_tlm 
         
        tlm_day = 10 ! default 
        call getin("tlm_day",tlm_day)
        write(*,*) "tlm_day : ", tlm_day 

        tlm_lt = 0.  ! default 
        call getin("tlm_lt",tlm_lt)
        write(*,*) "tlm_lt : ", tlm_lt 
            

!        Plume parameters and keys 
C          write(*,*) "Introducing an instantaneous plume?"
C          plume=.False. ! default value
C          call getin("plume",plume)
C          write(*,*) " plume = ",plume 
         
C          IF ( plume ) THEN 
         
C          write(*,*) "Days after day0 to release plume?"
C          plume_day= 0 ! default value
C          call getin("plume_day",plume_day)
C          write(*,*) " plume_day = ",plume_day 

C          write(*,*) "Local time to release plume?"
C          plume_lt= 0. ! default value
C          call getin("plume_lt",plume_lt)
C          write(*,*) " plume_lt = ",plume_lt 
         
C          write(*,*) "Plume lower cell"
C          p_low= 1 ! default value
C          call getin("p_low",p_low)
C          write(*,*) " p_low = ",p_low 
         
C          write(*,*) "Plume upper cell"
C          p_high= 1 ! default value
C          call getin("p_high",p_high)
C          write(*,*) " p_high = ",p_high 
         
C          ENDIF 

         
         
         
		 
         write(*,*) "Write some extra output to the screen ?"
         lwrite=.false. ! default value
         call getin("lwrite",lwrite)
         write(*,*) " lwrite = ",lwrite

         write(*,*) "Save statistics in file stats.nc ?"
#ifdef MESOSCALE
         callstats=.false. ! default value
#else
         callstats=.true. ! default value
#endif
         call getin("callstats",callstats)
         write(*,*) " callstats = ",callstats

         write(*,*) "Save EOF profiles in file 'profiles' for ",
     &              "Climate Database?"
         calleofdump=.false. ! default value
         call getin("calleofdump",calleofdump)
         write(*,*) " calleofdump = ",calleofdump

         write(*,*) "Dust scenario: 1=constant dust (read from startfi",
     &   " or set as tauvis); 2=Viking scenario; =3 MGS scenario,",
     &   "=6 cold (low dust) scenario; =7 warm (high dust) scenario ",
     &   "=24,25 ... 30 :Mars Year 24, ... or 30 from TES assimilation"
         iaervar=3 ! default value
         call getin("iaervar",iaervar)
         write(*,*) " iaervar = ",iaervar

         write(*,*) "Reference (visible) dust opacity at 610 Pa ",
     &   "(matters only if iaervar=1)"
         ! NB: default value of tauvis is set/read in startfi.nc file
         call getin("tauvis",tauvis)
         write(*,*) " tauvis = ",tauvis

         write(*,*) "Dust vertical distribution:"
         write(*,*) "(=1 top set by topdustref parameter;",
     & " =2 Viking scenario; =3 MGS scenario)"
         iddist=3 ! default value
         call getin("iddist",iddist)
         write(*,*) " iddist = ",iddist

         write(*,*) "Dust top altitude (km). (Matters only if iddist=1)"
         topdustref= 90.0 ! default value
         call getin("topdustref",topdustref)
         write(*,*) " topdustref = ",topdustref

         write(*,*) "Prescribed surface thermal flux (H/(rho*cp),K m/s)"
         tke_heat_flux=0. ! default value
         call getin("tke_heat_flux",tke_heat_flux)
         write(*,*) " tke_heat_flux = ",tke_heat_flux
         write(*,*) " 0 means the usual schemes are computing"

         write(*,*) "call radiative transfer ?"
         callrad=.true. ! default value
         call getin("callrad",callrad)
         write(*,*) " callrad = ",callrad

         write(*,*) "call slope insolation scheme ?",
     &              "(matters only if callrad=T)"
#ifdef MESOSCALE
         callslope=.true. ! default value
#else
         callslope=.false. ! default value (not supported yet)
#endif
         call getin("callslope",callslope)
         write(*,*) " callslope = ",callslope

         write(*,*) "call NLTE radiative schemes ?",
     &              "(matters only if callrad=T)"
         callnlte=.false. ! default value
         call getin("callnlte",callnlte)
         write(*,*) " callnlte = ",callnlte
         
         nltemodel=0    !default value
         write(*,*) "NLTE model?"
         write(*,*) "0 -> old model, static O"
         write(*,*) "1 -> old model, dynamic O"
         write(*,*) "2 -> new model"
         write(*,*) "(matters only if callnlte=T)"
         call getin("nltemodel",nltemodel)
         write(*,*) " nltemodel = ",nltemodel

         write(*,*) "call CO2 NIR absorption ?",
     &              "(matters only if callrad=T)"
         callnirco2=.false. ! default value
         call getin("callnirco2",callnirco2)
         write(*,*) " callnirco2 = ",callnirco2

         write(*,*) "New NIR NLTE correction ?",
     $              "0-> old model (no correction)",
     $              "1-> new correction",
     $              "(matters only if callnirco2=T)"
#ifdef MESOSCALE
         nircorr=0      !default value. this is OK below 60 km.
#else
         nircorr=0      !default value
#endif
         call getin("nircorr",nircorr)
         write(*,*) " nircorr = ",nircorr

         write(*,*) "call turbulent vertical diffusion ?"
         calldifv=.true. ! default value
         call getin("calldifv",calldifv)
         write(*,*) " calldifv = ",calldifv

         write(*,*) "call thermals ?"
         calltherm=.false. ! default value
         call getin("calltherm",calltherm)
         write(*,*) " calltherm = ",calltherm

         write(*,*) "output thermal diagnostics ?"
         outptherm=.false. ! default value
         call getin("outptherm",outptherm)
         write(*,*) " outptherm = ",outptherm

         write(*,*) "call convective adjustment ?"
         calladj=.true. ! default value
         call getin("calladj",calladj)
         write(*,*) " calladj = ",calladj
         
         if (calltherm .and. (.not. calladj)) then
          print*,'Convadj has to be activated when using thermals'
          stop
         endif

         write(*,*) "call Richardson-based surface layer ?"
         callrichsl=.false. ! default value
         call getin("callrichsl",callrichsl)
         write(*,*) " callrichsl = ",callrichsl

         if (calltherm .and. .not.callrichsl) then
          print*,'WARNING WARNING WARNING'
          print*,'if calltherm=T we strongly advise that '
          print*,'you use the new surface layer scheme '
          print*,'by setting callrichsl=T '
         endif

         if (calladj .and. callrichsl .and. (.not. calltherm)) then
          print*,'You should not be calling the convective adjustment
     & scheme with the Richardson surface-layer and without the thermals
     &. This approach is not
     & physically consistent and can lead to unrealistic friction
     & values.'
          print*,'If you want to use the Ri. surface-layer, either
     & activate thermals OR de-activate the convective adjustment.'
          stop
         endif

         write(*,*) "call CO2 condensation ?"
         callcond=.true. ! default value
         call getin("callcond",callcond)
         write(*,*) " callcond = ",callcond

         write(*,*)"call thermal conduction in the soil ?"
         callsoil=.true. ! default value
         call getin("callsoil",callsoil)
         write(*,*) " callsoil = ",callsoil
         

         write(*,*)"call Lott's gravity wave/subgrid topography ",
     &             "scheme ?"
         calllott=.true. ! default value
         call getin("calllott",calllott)
         write(*,*)" calllott = ",calllott


         write(*,*)"rad.transfer is computed every iradia",
     &             " physical timestep"
         iradia=1 ! default value
         call getin("iradia",iradia)
         write(*,*)" iradia = ",iradia
         

         write(*,*)"Output of the exchange coefficient mattrix ?",
     &             "(for diagnostics only)"
         callg2d=.false. ! default value
         call getin("callg2d",callg2d)
         write(*,*)" callg2d = ",callg2d

         write(*,*)"Rayleigh scattering : (should be .false. for now)"
         rayleigh=.false.
         call getin("rayleigh",rayleigh)
         write(*,*)" rayleigh = ",rayleigh


! TRACERS:

! dustbin
         write(*,*)"Transported dust ? (if >0, use 'dustbin' dust bins)"
         dustbin=0 ! default value
         call getin("dustbin",dustbin)
         write(*,*)" dustbin = ",dustbin
! active
         write(*,*)"Radiatively active dust ? (matters if dustbin>0)"
         active=.false. ! default value
         call getin("active",active)
         write(*,*)" active = ",active

! Test of incompatibility:
! if active is used, then dustbin should be > 0

         if (active.and.(dustbin.lt.1)) then
           print*,'if active is used, then dustbin should > 0'
           stop
         endif
! doubleq
         write(*,*)"use mass and number mixing ratios to predict",
     &             " dust size ?"
         doubleq=.false. ! default value
         call getin("doubleq",doubleq)
         write(*,*)" doubleq = ",doubleq
! submicron
         submicron=.false. ! default value
         call getin("submicron",submicron)
         write(*,*)" submicron = ",submicron

! Test of incompatibility:
! if doubleq is used, then dustbin should be 2

         if (doubleq.and.(dustbin.ne.2)) then
           print*,'if doubleq is used, then dustbin should be 2'
           stop
         endif
         if (doubleq.and.submicron.and.(nqmx.LT.3)) then
           print*,'If doubleq is used with a submicron tracer,'
           print*,' then the number of tracers has to be'
           print*,' larger than 3.'
           stop
         endif
! lifting
         write(*,*)"dust lifted by GCM surface winds ?"
         lifting=.false. ! default value
         call getin("lifting",lifting)
         write(*,*)" lifting = ",lifting

! Test of incompatibility:
! if lifting is used, then dustbin should be > 0

         if (lifting.and.(dustbin.lt.1)) then
           print*,'if lifting is used, then dustbin should > 0'
           stop
         endif
! callddevil
         write(*,*)" dust lifted by dust devils ?"
         callddevil=.false. !default value
         call getin("callddevil",callddevil)
         write(*,*)" callddevil = ",callddevil

! Test of incompatibility:
! if dustdevil is used, then dustbin should be > 0

         if (callddevil.and.(dustbin.lt.1)) then
           print*,'if dustdevil is used, then dustbin should > 0'
           stop
         endif
! sedimentation
         write(*,*) "Gravitationnal sedimentation ?"
         sedimentation=.true. ! default value
         call getin("sedimentation",sedimentation)
         write(*,*) " sedimentation = ",sedimentation
! activice
         write(*,*) "Radiatively active transported atmospheric ",
     &              "water ice ?"
         activice=.false. ! default value
         call getin("activice",activice)
         write(*,*) " activice = ",activice
! water
         write(*,*) "Compute water cycle ?"
         water=.false. ! default value
         call getin("water",water)
         write(*,*) " water = ",water

! thermal inertia feedback
         write(*,*) "Activate the thermal inertia feedback ?"
         tifeedback=.false. ! default value
         call getin("tifeedback",tifeedback)
         write(*,*) " tifeedback = ",tifeedback

! Test of incompatibility:

         if (tifeedback.and..not.water) then
           print*,'if tifeedback is used,'
           print*,'water should be used too'
           stop
         endif

         if (tifeedback.and..not.callsoil) then
           print*,'if tifeedback is used,'
           print*,'callsoil should be used too'
           stop
         endif

         if (activice.and..not.water) then
           print*,'if activice is used, water should be used too'
           stop
         endif

         if (water.and..not.tracer) then
           print*,'if water is used, tracer should be used too'
           stop
         endif
         
! water ice clouds effective variance distribution for sedimentaion       
        write(*,*) "effective variance for water ice clouds ?"
        nuice_sed=0.45 
        call getin("nuice_sed",nuice_sed)
        write(*,*) "water_param nueff Sedimentation:", nuice_sed
         
! ccn factor if no scavenging         
        write(*,*) "water param CCN reduc. factor ?", ccn_factor
        ccn_factor = 4.5
        call getin("ccn_factor",ccn_factor)
        write(*,*)" ccn_factor = ",ccn_factor
        write(*,*)"Careful: only used when microphys=F, otherwise"
        write(*,*)"the contact parameter is used instead;"

! microphys
         write(*,*)"Microphysical scheme for water-ice clouds?"
         microphys=.false. ! default value
         call getin("microphys",microphys)
         write(*,*)" microphys = ",microphys

! microphysical parameter contact       
         write(*,*) "water contact parameter ?"
         mteta  = 0.95
         call getin("mteta",mteta)
         write(*,*) "mteta = ", mteta

! scavenging
         write(*,*)"Dust scavenging by H2O/CO2 snowfall ?"
         scavenging=.false. ! default value
         call getin("scavenging",scavenging)
         write(*,*)" scavenging = ",scavenging
         

! Test of incompatibility:
! if scavenging is used, then dustbin should be > 0

         if ((microphys.and..not.doubleq).or.
     &       (microphys.and..not.water)) then
             print*,'if microphys is used, then doubleq,'
             print*,'and water must be used!'
             stop
         endif
         if (microphys.and..not.scavenging) then
             print*,''
             print*,'----------------WARNING-----------------'
             print*,'microphys is used without scavenging !!!'
             print*,'----------------WARNING-----------------'
             print*,''
         endif

         if ((scavenging.and..not.microphys).or.
     &       (scavenging.and.(dustbin.lt.1))) then
             print*,'if scavenging is used, then microphys'
             print*,'must be used!'
             stop
         endif

! Test of incompatibility:

         write(*,*) "Permanent water caps at poles ?",
     &               " .true. is RECOMMENDED"
         write(*,*) "(with .true., North cap is a source of water ",
     &   "and South pole is a cold trap)"
         caps=.true. ! default value
         call getin("caps",caps)
         write(*,*) " caps = ",caps

! albedo_h2o_ice
         write(*,*) "water ice albedo ?"
         albedo_h2o_ice=0.45
         call getin("albedo_h2o_ice",albedo_h2o_ice)
         write(*,*) " albedo_h2o_ice = ",albedo_h2o_ice
! inert_h2o_ice
         write(*,*) "water ice thermal inertia ?"
         inert_h2o_ice=2400 ! (J.m^-2.K^-1.s^-1/2)
         call getin("inert_h2o_ice",inert_h2o_ice)
         write(*,*) " inert_h2o_ice = ",inert_h2o_ice
! frost_albedo_threshold
         write(*,*) "frost thickness threshold for albedo ?"
         frost_albedo_threshold=0.005 ! 5.4 mic (i.e 0.005 kg.m-2)
         call getin("frost_albedo_threshold",
     &    frost_albedo_threshold)
         write(*,*) " frost_albedo_threshold = ",
     &            frost_albedo_threshold

! call Titus crocus line -- DEFAULT IS NONE
         write(*,*) "Titus crocus line ?"
         tituscap=.false.  ! default value
         call getin("tituscap",tituscap)
         write(*,*) "tituscap",tituscap
                     

         write(*,*) "photochemistry: include chemical species"
         photochem=.false. ! default value
         call getin("photochem",photochem)
         write(*,*) " photochem = ",photochem


! THERMOSPHERE

         write(*,*) "call thermosphere ?"
         callthermos=.false. ! default value
         call getin("callthermos",callthermos)
         write(*,*) " callthermos = ",callthermos
         

         write(*,*) " water included without cycle ",
     &              "(only if water=.false.)"
         thermoswater=.false. ! default value
         call getin("thermoswater",thermoswater)
         write(*,*) " thermoswater = ",thermoswater

         write(*,*) "call thermal conduction ?",
     &    " (only if callthermos=.true.)"
         callconduct=.false. ! default value
         call getin("callconduct",callconduct)
         write(*,*) " callconduct = ",callconduct

         write(*,*) "call EUV heating ?",
     &   " (only if callthermos=.true.)"
         calleuv=.false.  ! default value
         call getin("calleuv",calleuv)
         write(*,*) " calleuv = ",calleuv

         write(*,*) "call molecular viscosity ?",
     &   " (only if callthermos=.true.)"
         callmolvis=.false. ! default value
         call getin("callmolvis",callmolvis)
         write(*,*) " callmolvis = ",callmolvis

         write(*,*) "call molecular diffusion ?",
     &   " (only if callthermos=.true.)"
         callmoldiff=.false. ! default value
         call getin("callmoldiff",callmoldiff)
         write(*,*) " callmoldiff = ",callmoldiff
         

         write(*,*) "call thermospheric photochemistry ?",
     &   " (only if callthermos=.true.)"
         thermochem=.false. ! default value
         call getin("thermochem",thermochem)
         write(*,*) " thermochem = ",thermochem

         write(*,*) "Method to include solar variability"
         write(*,*) "0-> old method (using solarcondate); ",
     &                  "1-> variability wit E10.7"
         solvarmod=1
         call getin("solvarmod",solvarmod)
         write(*,*) " solvarmod = ",solvarmod

         write(*,*) "date for solar flux calculation:",
     &   " (1985 < date < 2002)",
     $   " (Only used if solvarmod=0)"
         write(*,*) "(Solar min=1996.4 ave=1993.4 max=1990.6)"
         solarcondate=1993.4 ! default value
         call getin("solarcondate",solarcondate)
         write(*,*) " solarcondate = ",solarcondate
         
         write(*,*) "Solar variability as observed for MY: "
         write(*,*) "Only if solvarmod=1"
         solvaryear=24
         call getin("solvaryear",solvaryear)
         write(*,*) " solvaryear = ",solvaryear

         write(*,*) "UV heating efficiency:",
     &   "measured values between 0.19 and 0.23 (Fox et al. 1996)",
     &   "lower values may be used to compensate low 15 um cooling"
         euveff=0.21 !default value
         call getin("euveff",euveff)
         write(*,*) " euveff = ", euveff

         if (.not.callthermos) then
           if (thermoswater) then
             print*,'if thermoswater is set, callthermos must be true'
             stop
           endif          
           if (callconduct) then
             print*,'if callconduct is set, callthermos must be true'
             stop
           endif        
           if (calleuv) then
             print*,'if calleuv is set, callthermos must be true'
             stop
           endif         
           if (callmolvis) then
             print*,'if callmolvis is set, callthermos must be true'
             stop
           endif        
           if (callmoldiff) then
             print*,'if callmoldiff is set, callthermos must be true'
             stop
           endif          
           if (thermochem) then
             print*,'if thermochem is set, callthermos must be true'
             stop
           endif          
        endif

! Test of incompatibility:
! if photochem is used, then water should be used too

         if (photochem.and..not.water) then
           print*,'if photochem is used, water should be used too'
           stop
         endif

! if callthermos is used, then thermoswater should be used too 
! (if water not used already)

         if (callthermos .and. .not.water) then
           if (callthermos .and. .not.thermoswater) then
             print*,'if callthermos is used, water or thermoswater 
     &               should be used too'
             stop
           endif
         endif

         PRINT*,'--------------------------------------------'
         PRINT*
         PRINT*
      ELSE
         write(*,*)
         write(*,*) 'Cannot read file callphys.def. Is it here ?'
         stop
      ENDIF

8000  FORMAT(t5,a12,l8)
8001  FORMAT(t5,a12,i8)

      PRINT*
      PRINT*,'inifis: daysec',daysec
      PRINT*
      PRINT*,'inifis: The radiative transfer is computed:'
      PRINT*,'           each ',iradia,' physical time-step'
      PRINT*,'        or each ',iradia*dtphys,' seconds'
      PRINT*
! --------------------------------------------------------------
!  Managing the Longwave radiative transfer
! --------------------------------------------------------------

!     In most cases, the run just use the following values :
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      callemis=.true.     
!     ilwd=10*int(daysec/dtphys) ! bug before 22/10/01       
      ilwd=1
      ilwn=1 !2
      ilwb=1 !2
      linear=.true.        
      ncouche=3
      alphan=0.4
      semi=0

!     BUT people working hard on the LW may want to read them in 'radia.def' 
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      OPEN(99,file='radia.def',status='old',form='formatted'
     .     ,iostat=ierr)
      IF(ierr.EQ.0) THEN
         write(*,*) 'inifis: Reading radia.def !!!'
         READ(99,fmt='(a)') ch1
         READ(99,*) callemis
         WRITE(*,8000) ch1,callemis

         READ(99,fmt='(a)') ch1
         READ(99,*) iradia
         WRITE(*,8001) ch1,iradia

         READ(99,fmt='(a)') ch1
         READ(99,*) ilwd
         WRITE(*,8001) ch1,ilwd

         READ(99,fmt='(a)') ch1
         READ(99,*) ilwn
         WRITE(*,8001) ch1,ilwn

         READ(99,fmt='(a)') ch1
         READ(99,*) linear
         WRITE(*,8000) ch1,linear

         READ(99,fmt='(a)') ch1
         READ(99,*) ncouche
         WRITE(*,8001) ch1,ncouche

         READ(99,fmt='(a)') ch1
         READ(99,*) alphan
         WRITE(*,*) ch1,alphan

         READ(99,fmt='(a)') ch1
         READ(99,*) ilwb
         WRITE(*,8001) ch1,ilwb


         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callg2d
         WRITE(*,8000) ch1,callg2d

         READ(99,fmt='(a)') ch1
         READ(99,*) semi
         WRITE(*,*) ch1,semi
      end if
      CLOSE(99)

!-----------------------------------------------------------------------
!     Some more initialization:
!     ------------------------

      ! in 'comgeomfi.h'
      CALL SCOPY(ngrid,plon,1,long,1)
      CALL SCOPY(ngrid,plat,1,lati,1)
      CALL SCOPY(ngrid,parea,1,area,1)
      totarea=SSUM(ngridmx,area,1)

      ! in 'comdiurn.h'
      DO ig=1,ngrid
         sinlat(ig)=sin(plat(ig))
         coslat(ig)=cos(plat(ig))
         sinlon(ig)=sin(plon(ig))
         coslon(ig)=cos(plon(ig))
      ENDDO

      pi=2.*asin(1.) ! NB: pi is a common in comcstfi.h

!     managing the tracers, and tests:
!     -------------------------------
!     Ehouarn: removed; as these tests are now done in initracer.F
!      if(tracer) then
!
!!          when photochem is used, nqchem_min is the rank
!!          of the first chemical species
!
!! Ehouarn: nqchem_min is now meaningless and no longer used
!!       nqchem_min = 1
!       if (photochem .or. callthermos) then
!         chem = .true.
!       end if
!
!       if (water .or. thermoswater) h2o = .true.
!
!!          TESTS
!
!       print*,'inifis: TRACERS:'
!       write(*,*) "    chem=",chem,"    h2o=",h2o
!!       write(*,*) "   doubleq=",doubleq
!!       write(*,*) "   dustbin=",dustbin
!
!       if ((doubleq).and.(h2o).and.
!     $     (chem)) then
!         print*,' 2 dust tracers (doubleq)'
!         print*,' 1 water vapour tracer'
!         print*,' 1 water ice tracer'
!         print*,nqmx-4,' chemistry tracers'
!       endif
!
!       if ((doubleq).and.(h2o).and.
!     $     .not.(chem)) then
!         print*,' 2 dust tracers (doubleq)'
!         print*,' 1 water vapour tracer'
!         print*,' 1 water ice tracer'
!         if (nqmx.LT.4) then
!           print*,'nqmx should be at least equal to'
!           print*,'4 with these options.'
!           stop
!         endif
!       endif
!
!       if (.not.(doubleq).and.(h2o).and.
!     $     (chem)) then
!         if (dustbin.gt.0) then
!           print*,dustbin,' dust bins'
!         endif
!         print*,nqmx-2-dustbin,' chemistry tracers'
!         print*,' 1 water vapour tracer'
!         print*,' 1 water ice tracer'
!       endif
!
!       if (.not.(doubleq).and.(h2o).and.
!     $     .not.(chem)) then
!         if (dustbin.gt.0) then
!           print*,dustbin,' dust bins'
!         endif
!         print*,' 1 water vapour tracer'
!         print*,' 1 water ice tracer'
!         if (nqmx.gt.(dustbin+2)) then
!           print*,'nqmx should be ',(dustbin+2),
!     $            ' with these options...'
!		   print*,'(or check callphys.def)'
!         endif
!         if (nqmx.lt.(dustbin+2)) then
!           write(*,*) "inifis: nqmx.lt.(dustbin+2)"
!           stop
!         endif
!       endif
!
!      endif ! of if (tracer)
!
!      RETURN
      END
