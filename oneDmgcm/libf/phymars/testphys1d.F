
      PROGRAM testphys1d
! to use  'getin'
      USE TLMvars
      USE ioipsl_getincom 
      IMPLICIT NONE

c=======================================================================
c   subject:
c   --------
c   PROGRAM useful to run physical part of the martian GCM in a 1D column
c       
c Can be compiled with a command like (e.g. for 25 layers)
c  "makegcm -p mars -d 25 testphys1d"
c It requires the files "testphys1d.def" "callphys.def"
c   and a 'run.def' file (containing a "INCLUDEDEF=callphys.def" line)
c      and a file describing the sigma layers (e.g. "z2sig.def")
c
c   author: Frederic Hourdin, R.Fournier,F.Forget
c   -------
c   
c   update: 12/06/2003 including chemistry (S. Lebonnois) 
c                            and water ice (F. Montmessin)
c 
c=======================================================================

#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "comgeomfi.h"
#include "surfdat.h"
#include "slope.h"
#include "comsoil.h"
#include "comdiurn.h"
#include "callkeys.h"
#include "comcstfi.h"
#include "planete.h"
#include "comsaison.h"
#include "yomaer.h"
#include "control.h"
#include "comvert.h"
#include "netcdf.inc"
#include "comg1d.h"
#include "logic.h"
#include "advtrac.h"
#include "tracer.h"

c --------------------------------------------------------------
c  Declarations
c --------------------------------------------------------------
c
      INTEGER unitstart      ! unite d'ecriture de "startfi"
      INTEGER nlayer,nlevel,nsoil,ndt
      INTEGER ilayer,ilevel,isoil,idt,iq
      LOGICAl firstcall,lastcall
c
      real,parameter :: odpref=610. ! DOD reference pressure (Pa)
c
      INTEGER day0          ! date initial (sol ; =0 a Ls=0)
      REAL day           ! date durant le run
      REAL date          ! number of days (and fractions of) that have occured
      REAL time             ! time (0<time<1 ; time=0.5 a midi)
      REAL play(nlayermx)   ! Pressure at the middle of the layers (Pa)
      REAL plev(nlayermx+1) ! intermediate pressure levels (pa)
      REAL psurf,tsurf      
      REAL u(nlayermx),v(nlayermx)  ! zonal, meridional wind
      REAL gru,grv   ! prescribed "geostrophic" background wind
      REAL temp(nlayermx)   ! temperature at the middle of the layers
      REAL q(nlayermx,nqmx) ! tracer mixing ratio (e.g. kg/kg)
      REAL qsurf(nqmx)      ! tracer surface budget (e.g. kg.m-2)
      REAL tsoil(nsoilmx)   ! subsurface soik temperature (K)
      REAL co2ice           ! co2ice layer (kg.m-2)
      REAL emis             ! surface layer
      REAL q2(nlayermx+1)   ! Turbulent Kinetic Energy
      REAL zlay(nlayermx)   ! altitude estimee dans les couches (km)

c    Physical and dynamical tandencies (e.g.  m.s-2, K/s, Pa/s)
      REAL du(nlayermx),dv(nlayermx),dtemp(nlayermx)
      REAL dudyn(nlayermx),dvdyn(nlayermx),dtempdyn(nlayermx)
      REAL dpsurf    
      REAL dq(nlayermx,nqmx)
      REAL dqdyn(nlayermx,nqmx)

c   Various intermediate variables
      INTEGER thermo
      REAL zls
      REAL phi(nlayermx),h(nlayermx),s(nlayermx)
      REAL pks, ptif, w(nlayermx)
      REAL qtotinit, mqtot(nqmx),qtot
      INTEGER ierr, aslun
      REAL tmp1(0:nlayermx),tmp2(0:nlayermx)
      Logical  tracerdyn
      integer :: nq=1 ! number of tracers

      character*2 str2
      character (len=7) :: str7
      character(len=44) :: txt
      !pip
      character(len=100) :: infile      

      ! Plume Parameters
      integer ndt_plume ! Time-step for tracer insertion 

      ! TLM Parameters
      integer ndt_tlm ! Time-step for tlm activation

      
c=======================================================================

c=======================================================================
c INITIALISATION
c=======================================================================

c ------------------------------------------------------
c  Prescribed constants to be set here
c ------------------------------------------------------

      pi=2.E+0*asin(1.E+0)

c     Mars planetary constants
c     ----------------------------
      rad=3397200.               ! mars radius (m)  ~3397200 m
      daysec=88775.              ! length of a sol (s)  ~88775 s
      omeg=4.*asin(1.)/(daysec)  ! rotation rate (rad.s-1)
      g=3.72                     ! gravity (m.s-2) ~3.72  
      mugaz=43.49                ! atmosphere mola mass (g.mol-1) ~43.49
      rcp=.256793                ! = r/cp  ~0.256793
      r= 8.314511E+0 *1000.E+0/mugaz
      cpp= r/rcp
      year_day = 669             ! lenght of year (sols) ~668.6
      periheli = 206.66          ! minimum sun-mars distance (Mkm) ~206.66
      aphelie = 249.22           ! maximum sun-mars distance (Mkm) ~249.22
      peri_day =  485.           ! perihelion date (sols since N. Spring)
      obliquit = 25.2            ! Obliquity (deg) ~25.2         
 
c     Planetary Boundary Layer and Turbulence parameters 
c     --------------------------------------------------
      z0_default =  1.e-2        ! surface roughness (m) ~0.01 
      emin_turb = 1.e-6          ! minimal turbulent energy ~1.e-8
      lmixmin = 30               ! mixing length ~100
 
c     cap properties and surface emissivities
c     ----------------------------------------------------
      emissiv= 0.95              ! Bare ground emissivity ~.95
      emisice(1)=0.95            ! Northern cap emissivity
      emisice(2)=0.95            ! Southern cap emisssivity
      albedice(1)=0.5            ! Northern cap albedo
      albedice(2)=0.5            ! Southern cap albedo
      iceradius(1) = 100.e-6     ! mean scat radius of CO2 snow (north)
      iceradius(2) = 100.e-6     ! mean scat radius of CO2 snow (south)
      dtemisice(1) = 2.          ! time scale for snow metamorphism (north)
      dtemisice(2) = 2.          ! time scale for snow metamorphism (south


c ------------------------------------------------------
c  Loading run parameters from "run.def" file
c ------------------------------------------------------


! check if 'run.def' file is around (otherwise reading parameters
! from callphys.def via getin() routine won't work.
!      open(99,file='run.def',status='old',form='formatted',
!     &     iostat=ierr)
!      if (ierr.ne.0) then
!        write(*,*) 'Cannot find required file "run.def"'
!        write(*,*) '  (which should contain some input parameters'
!        write(*,*) '   along with the following line:'
!        write(*,*) '   INCLUDEDEF=callphys.def'
!        write(*,*) '   )'
!        write(*,*) ' ... might as well stop here ...'
!        stop
!      else
!        close(99)
!      endif


      ! pip
      infile = 'run.def'
      call checkfile(infile,99)

      
! check if we are going to run with or without tracers
      write(*,*) "Run with or without tracer transport ?"
      tracer=.false. ! default value
      call getin("tracer",tracer)
      write(*,*) " tracer = ",tracer

! while we're at it, check if there is a 'traceur.def' file
! and preocess it, if necessary. Otherwise initialize tracer names
      if (tracer) then
      ! load tracer names from file 'traceur.def'
        open(90,file='traceur.def',status='old',form='formatted',
     &       iostat=ierr)
        if (ierr.ne.0) then
          write(*,*) 'Cannot find required file "traceur.def"'
          write(*,*) ' If you want to run with tracers, I need it'
          write(*,*) ' ... might as well stop here ...'
          stop
        else
          write(*,*) "testphys1d: Reading file traceur.def"
          ! read number of tracers:
          read(90,*,iostat=ierr) nq
          if (ierr.ne.0) then
            write(*,*) "testphys1d: error reading number of tracers"
            write(*,*) "   (first line of traceur.def) "
            stop
          else
            ! check that the number of tracers is indeed nqmx
            if (nq.ne.nqmx) then
              write(*,*) "testphys1d: error, wrong number of tracers:"
              write(*,*) "nq=",nq," whereas nqmx=",nqmx
              stop
            endif
          endif
        endif
        ! read tracer names from file traceur.def
        do iq=1,nqmx
          read(90,*,iostat=ierr) tnom(iq)
          if (ierr.ne.0) then
            write(*,*) 'testphys1d: error reading tracer names...'
            stop
          endif
        enddo
        close(90)

		! Added in so that noms and tnom are both filled JAH 6/16
        do iq=1,nqmx
            noms(iq) = tnom(iq)
        enddo
		
        ! initialize tracers here:
        write(*,*) "testphys1d: initializing tracers"
        q(:,:)=0 ! default, set everything to zero
        qsurf(:)=0
        ! "smarter" initialization of some tracers
        ! (get values from "profile_*" files, if these are available)
        do iq=1,nqmx
          txt=""
          write(txt,"(a)") tnom(iq)
          write(*,*)"  tracer:",trim(txt)
          ! CO2
          if (txt.eq."co2") then
            q(:,iq)=0.95   ! kg /kg of atmosphere
            qsurf(iq)=0. ! kg/m2 (not used for CO2)
            ! even better, look for a "profile_co2" input file
            open(91,file='profile_co2',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayermx
                read(91,*) q(ilayer,iq)
              enddo
            endif
            close(91)
          endif ! of if (txt.eq."co2")
          ! Allow for an initial profile of argon
          ! Can also be used to introduce a decaying tracer
          ! in the 1D (TBD) to study thermals
          if (txt.eq."ar") then
            !look for a "profile_ar" input file
            open(91,file='profile_ar',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayermx
                read(91,*) q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_ar file!"
            endif
            close(91)
          endif ! of if (txt.eq."ar")

          ! WATER VAPOUR
          if (txt.eq."h2o_vap") then
            !look for a "profile_h2o_vap" input file
            open(91,file='profile_h2o_vap',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayermx
                read(91,*) q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_h2o_vap file!"
            endif
            close(91)
          endif ! of if (txt.eq."h2o_ice")
          ! WATER ICE
          if (txt.eq."h2o_ice") then
            !look for a "profile_h2o_vap" input file
            open(91,file='profile_h2o_ice',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayermx
                read(91,*) q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_h2o_ice file!"
            endif
            close(91)
          endif ! of if (txt.eq."h2o_ice")
          ! DUST
          !if (txt(1:4).eq."dust") then
          !  q(:,iq)=0.4    ! kg/kg of atmosphere
          !  qsurf(iq)=100 ! kg/m2
          !endif
          ! DUST MMR
          if (txt.eq."dust_mass") then
            !look for a "profile_dust_mass" input file
            open(91,file='profile_dust_mass',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayermx
                read(91,*) q(ilayer,iq)
!                write(*,*) "l=",ilayer," q(ilayer,iq)=",q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_dust_mass file!"
            endif
            close(91)
          endif ! of if (txt.eq."dust_mass")
          ! DUST NUMBER
          if (txt.eq."dust_number") then
            !look for a "profile_dust_number" input file
            open(91,file='profile_dust_number',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayermx
                read(91,*) q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_dust_number file!"
            endif
            close(91)
          endif ! of if (txt.eq."dust_number")
          ! NB: some more initializations (chemistry) is done later
          ! CCN MASS
          if (txt.eq."ccn_mass") then
            !look for a "profile_ccn_mass" input file
            open(91,file='profile_ccn_mass',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayermx
                read(91,*) q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_ccn_mass file!"
            endif
            close(91)
          endif ! of if (txt.eq."ccn_mass")
          ! CCN NUMBER
          if (txt.eq."ccn_number") then
            !look for a "profile_ccn_number" input file
            open(91,file='profile_ccn_number',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayermx
                read(91,*) q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_ccn_number file!"
            endif
            close(91)
          endif ! of if (txt.eq."ccn_number")
        enddo ! of do iq=1,nqmx

      else
      ! we still need to set (dummy) tracer names for physdem1
        nq=nqmx
        do iq=1,nq
          write(str7,'(a1,i2.2)')'q',iq
          tnom(iq)=str7
        enddo
      ! and just to be clean, also initialize tracers to zero for physdem1
        q(:,:)=0 
        qsurf(:)=0      
      endif ! of if (tracer)
      
      !write(*,*) "testphys1d q", q(1,:)
      !write(*,*) "testphys1d qsurf", qsurf

c  Date and local time at beginning of run
c  ---------------------------------------
c    Date (in sols since spring solstice) at beginning of run
      day0 = 0 ! default value for day0
      write(*,*) 'Initial date (in martian sols ; =0 at Ls=0)?'
      call getin("day0",day0)
      day=float(day0)
      write(*,*) " day0 = ",day0
c  Local time at beginning of run 
      time=0 ! default value for time
      write(*,*)'Initial local time (in hours, between 0 and 24)?'
      call getin("time",time)
      write(*,*)" time = ",time
      time=time/24.E+0 ! convert time (hours) to fraction of sol

c  Discretization (Definition of grid and time steps)
c  --------------
c
      nlayer=nlayermx
      nlevel=nlayer+1
      nsoil=nsoilmx

      day_step=48 ! default value for day_step
      PRINT *,'Number of time steps per sol ?'
      call getin("day_step",day_step)
      write(*,*) " day_step = ",day_step

      ndt=10 ! default value for ndt
      PRINT *,'Number of sols to run ?'
      call getin("ndt",ndt)
      write(*,*) " ndt = ",ndt

      ndt=ndt*day_step     
      dtphys=daysec/day_step  

c Imposed surface pressure [MOVED BELOW]
c ------------------------------------
c
      psurf=610. ! default value for psurf
      PRINT *,'Surface pressure (Pa) ?'
      call getin("psurf",psurf)
      write(*,*) " psurf = ",psurf
      
      
      
      
      
c Reference pressures
      pa=20.   ! transition pressure (for hybrid coord.)
      preff=610.      ! reference surface pressure
 
c Aerosol properties
c --------------------------------
      tauvis=0.2 ! default value for tauvis (dust opacity)
      write(*,'("Reference dust opacity at ",f4.0," Pa ?")')odpref
      call getin("tauvis",tauvis)
      write(*,*) " tauvis = ",tauvis

c Orbital parameters
c ------------------
      print *,'Min. distance Sun-Mars (Mkm)?'
      call getin("periheli",periheli)
      write(*,*) " periheli = ",periheli

      print *,'Max. distance Sun-Mars (Mkm)?'
      call getin("aphelie",aphelie)
      write(*,*) " aphelie = ",aphelie

      print *,'Day of perihelion?'
      call getin("periday",peri_day)
      write(*,*) " periday = ",peri_day

      print *,'Obliquity?'
      call getin("obliquit",obliquit)
      write(*,*) " obliquit = ",obliquit
 
c  latitude/longitude
c  ------------------
      lati(1)=0 ! default value for lati(1)
      PRINT *,'latitude (in degrees) ?'
      call getin("latitude",lati(1))
      write(*,*) " latitude = ",lati(1)
      lati(1)=lati(1)*pi/180.E+0
      long(1)=0.E+0
      long(1)=long(1)*pi/180.E+0

c  Initialize albedo / soil thermal inertia
c  ----------------------------------------
c
      albedodat(1)=0.2 ! default value for albedodat
      PRINT *,'Albedo of bare ground ?'
      call getin("albedo",albedodat(1))
      write(*,*) " albedo = ",albedodat(1)

      inertiedat(1,1)=400 ! default value for inertiedat
      PRINT *,'Soil thermal inertia (SI) ?'
      call getin("inertia",inertiedat(1,1))
      write(*,*) " inertia = ",inertiedat(1,1)

      z0(1)=z0_default ! default value for roughness
      write(*,*) 'Surface roughness length z0 (m)?'
      call getin("z0",z0(1))
      write(*,*) " z0 = ",z0(1)

! Initialize local slope parameters (only matters if "callslope"
! is .true. in callphys.def)
      ! slope inclination angle (deg) 0: horizontal, 90: vertical
      theta_sl(1)=0.0 ! default: no inclination
      call getin("slope_inclination",theta_sl(1))
      ! slope orientation (deg)
      ! 0 == Northward, 90 == Eastward, 180 == Southward, 270 == Westward
      psi_sl(1)=0.0 ! default value
      call getin("slope_orientation",psi_sl(1))
      
c
c  for the gravity wave scheme
c  ---------------------------------
c
      zmea(1)=0.E+0
      zstd(1)=0.E+0
      zsig(1)=0.E+0
      zgam(1)=0.E+0
      zthe(1)=0.E+0


c   Specific initializations for "physiq"
c   -------------------------------------
c   mesh surface (not a very usefull quantity in 1D)
      area(1)=1.E+0

c   surface geopotential is not used (or useful) since in 1D
c   everything is controled by surface pressure
      phisfi(1)=0.E+0

c  "inifis" does some initializations (some of which have already been
c  done above!) and loads parameters set in callphys.def

!Mars possible matter with dtphys in input and include!!!
      CALL inifis(1,llm,day0,daysec,dtphys,
     .            lati,long,area,rad,g,r,cpp)

c   Initialization to take into account prescribed winds
c   ------------------------------------------------------
      ptif=2.E+0*omeg*sinlat(1)
 
c    geostrophic wind
      gru=10. ! default value for gru
      PRINT *,'zonal eastward component of the geostrophic wind (m/s) ?'
      call getin("u",gru)
      write(*,*) " u = ",gru
      grv=0. !default value for grv
      PRINT *,'meridional northward component of the geostrophic',
     &' wind (m/s) ?'
      call getin("v",grv)
      write(*,*) " v = ",grv

c     Initialize winds  for first time step
      DO ilayer=1,nlayer
         u(ilayer)=gru
         v(ilayer)=grv
      ENDDO

c     Initialize turbulente kinetic energy
      DO ilevel=1,nlevel
         q2(ilevel)=0.E+0
      ENDDO

c  CO2 ice on the surface
c  -------------------
      co2ice=0.E+0 ! default value for co2ice
      PRINT *,'Initial CO2 ice on the surface (kg.m-2)'
      call getin("co2ice",co2ice)
      write(*,*) " co2ice = ",co2ice

c
c  emissivity
c  ----------
      emis=emissiv
      IF (co2ice.eq.1.E+0) THEN
         emis=emisice(1) ! northern hemisphere
         IF(lati(1).LT.0) emis=emisice(2) ! southern hemisphere
      ENDIF

 

c  Compute pressures and altitudes of atmospheric levels 
c  ----------------------------------------------------------------
!     Initialise the surface pressure from the MCDv5.3 look-up table
      CALL iniorbit(aphelie,periheli,year_day,peri_day,obliquit)
      call solarlong(day+time,zls)
      call mcd_phys(time*24.,zlay,zls*180./pi,
     &              temp, psurf, u, v, 
     &              q2, day_step)

c    Vertical Coordinates
c    """"""""""""""""""""
      hybrid=.true.
      PRINT *,'Hybrid coordinates ?'
      call getin("hybrid",hybrid)
      write(*,*) " hybrid = ", hybrid

      CALL  disvert

      DO ilevel=1,nlevel
        plev(ilevel)=ap(ilevel)+psurf*bp(ilevel)
      ENDDO

      DO ilayer=1,nlayer
        play(ilayer)=aps(ilayer)+psurf*bps(ilayer)
      ENDDO

      DO ilayer=1,nlayer
        zlay(ilayer)=-200.E+0 *r*log(play(ilayer)/plev(1))
     &   /g
      ENDDO


c  Initialize temperature profile
c  --------------------------------------
      pks=psurf**rcp

C altitude in km in profile: divide zlay by 1000
      tmp1(0)=0.E+0
      DO ilayer=1,nlayer
        tmp1(ilayer)=zlay(ilayer)/1000.E+0
      ENDDO

      call profile(nlayer+1,tmp1,tmp2)

      tsurf=tmp2(0)
      DO ilayer=1,nlayer
        temp(ilayer)=tmp2(ilayer)
      ENDDO
      
c        
C     Second call to the MCD v5.3 routine will initialise the
C     temperature, wind, and q2 profiles.
      call mcd_phys(time*24.,zlay,zls*180./pi,
     &              temp, psurf, u, v, 
     &              q2, day_step)
     


C       tsurf = temp(1)
      
! Initialize soil properties and temperature
! ------------------------------------------
      volcapa=1.e6 ! volumetric heat capacity
      DO isoil=1,nsoil
         inertiedat(1,isoil)=inertiedat(1,1) ! soil thermal inertia
         tsoil(isoil)=tsurf  ! soil temperature
      ENDDO

! Initialize depths
! -----------------
      do isoil=0,nsoil-1
        mlayer(isoil)=2.e-4*(2.**(isoil-0.5)) ! mid-layer depth
      enddo
      do isoil=1,nsoil
        layer(isoil)=2.e-4*(2.**(isoil-1)) ! layer depth
      enddo

c    Initialize traceurs
c    ---------------------------

      
      if (photochem.or.callthermos) then
         write(*,*) 'Initializing chemical species'
         ! thermo=0: initialize over all atmospheric layers
         thermo=0
         call inichim_newstart(q,qsurf,psurf,1,0)
      endif



      
c Check if the surface is a water ice reservoir 
c --------------------------------------------------
      watercaptag(ngridmx)=.false. ! Default: no water ice reservoir
      print *,'Water ice cap on ground ?'
      call getin("watercaptag",watercaptag)
      write(*,*) " watercaptag = ",watercaptag
      

c    Initialization for GRADS outputs in "g1d.dat" and "g1d.ctl"
c    ----------------------------------------------------------------
c    (output done in "writeg1d", typically called by "physiq.F")

        g1d_nlayer=nlayer
        g1d_nomfich='g1d.dat'
        g1d_unitfich=40
        g1d_nomctl='g1d.ctl'
        g1d_unitctl=41
        g1d_premier=.true.
        g2d_premier=.true.

c  Write a "startfi" file
c  --------------------
c  This file will be read during the first call to "physiq".
c  It is needed to transfert physics variables to "physiq"...

      call physdem1("startfi.nc",long,lati,nsoilmx,nqmx,
     .              dtphys,float(day0),time,tsurf,
     .              tsoil,co2ice,emis,q2,qsurf,area,albedodat,
     .              inertiedat,zmea,zstd,zsig,zgam,zthe)

     
      ndt_plume =  (plume_day*day_step) + NINT((plume_lt/24.)*day_step)

c=======================================================================
c  1D MODEL TIME STEPPING LOOP
c=======================================================================
c
      firstcall=.true.
      lastcall=.false.

      if ( call_tlm ) then
c     Plume release activates after this many model steps 
      ndt_plume =  (plume_day*day_step) + NINT((plume_lt/24.)*day_step)
c     TLM activates after this many model steps 
      ndt_tlm = (tlm_day*day_step) + NINT((tlm_lt/24.)*day_step)
      endif 

      DO idt=1,ndt
c        IF (idt.eq.ndt) lastcall=.true.
        IF (idt.eq.ndt-day_step-1) then       !test
         lastcall=.true.
         call solarlong(day*1.0,zls)
         write(103,*) 'Ls=',zls*180./pi
         write(103,*) 'Lat=', lati(1)*180./pi
         write(103,*) 'Tau=', tauvis/odpref*psurf
         write(103,*) 'RunEnd - Atmos. Temp. File'
         write(103,*) 'RunEnd - Atmos. Temp. File'
         write(104,*) 'Ls=',zls*180./pi
         write(104,*) 'Lat=', lati(1)
         write(104,*) 'Tau=', tauvis/odpref*psurf
         write(104,*) 'RunEnd - Atmos. Temp. File'
        ENDIF

c     compute geopotential
c     ~~~~~~~~~~~~~~~~~~~~~
      DO ilayer=1,nlayer
        s(ilayer)=(aps(ilayer)/psurf+bps(ilayer))**rcp
        h(ilayer)=cpp*temp(ilayer)/(pks*s(ilayer))
      ENDDO
      phi(1)=pks*h(1)*(1.E+0-s(1))
      DO ilayer=2,nlayer
         phi(ilayer)=phi(ilayer-1)+
     &               pks*(h(ilayer-1)+h(ilayer))*.5E+0
     &                  *(s(ilayer-1)-s(ilayer))

      ENDDO

c  Begin the tracer interpolation routine for the MCD v5.3 look-up
c  table for tracers CO2, CO, H2, H2O vapour, and O2.
          IF ( mcd_profiles ) THEN 
          call mcd_trac(time,zlay,zls*180./pi, 
     &                     q(:,1), q(:,2), 
     &                     q(:,5), q(:,8),
     &                     q(:,15))
          WRITE(*,*) "mcd [CO2, CO, O2, H2, H2O_vap]"
          ENDIF 
      

c     Tangent Linear Model Activation
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        IF ( call_tlm ) THEN 
          IF ( idt == ndt_tlm ) THEN 
               
               call_tlm = .False.     ! Prevents multiple initialisations.
               
               tlm_on = .True.        ! Routines will now engage.
               
               mcd_profiles = .False. ! Interpolation routines will 
                                      ! invalidate TLM results.
               
               calladj = .False. ! Convective adjustments require a 
                                  ! conditional that the TLM will not
                                  ! be capable of performing.
               callcond = .False. ! Again, conditionals activate this 
                                  ! module.
               call TLM_initialise

          ENDIF
          
        ENDIF 


c       call physics
c       --------------------
!      write(*,*) "testphys1d avant q", q(1,:)
      CALL physiq (1,llm,nqmx,
     ,     firstcall,lastcall,
     ,     day,time,dtphys,
     ,     plev,play,phi,
     ,     u, v,temp, q,  
     ,     w,
C - outputs
     s     du, dv, dtemp, dq,dpsurf,tracerdyn)
!      write(*,*) "testphys1d apres q", q(1,:)

c     Interpolate physical variables (temperature, wind,
c     q2 and surface pressure) 
      call mcd_phys(time*24.,zlay,zls*180./pi,
     &              temp, psurf, u, v, 
     &              q2, day_step)


      IF (tlm_on) call TLM_ADJ_out(idt,time*24.)


c       wind increment : specific for 1D
c       --------------------------------
 
c       The physics compute the tendencies on u and v,
c       here we just add Coriolos effect
c
C       DO ilayer=1,nlayer
C          du(ilayer)=du(ilayer)+ptif*(v(ilayer)-grv)
C          dv(ilayer)=dv(ilayer)+ptif*(-u(ilayer)+gru)
C       ENDDO

c       For some tests : No coriolis force at equator
c       if(lati(1).eq.0.) then
C           DO ilayer=1,nlayer
C              du(ilayer)=du(ilayer)+ (gru-u(ilayer))/1.e4
C              dv(ilayer)=dv(ilayer)+ (grv-v(ilayer))/1.e4
C           ENDDO
c       end if
c      
c
c       Compute time for next time step
c       ---------------------------------------
        firstcall=.false.
        time=time+dtphys/daysec
        IF (time.gt.1.E+0) then
            time=time-1.E+0
            day=day+1
        ENDIF

c       compute winds and temperature for next time step
c       ----------------------------------------------------------

C         DO ilayer=1,nlayer
C            u(ilayer)=u(ilayer)+dtphys*du(ilayer)
C            v(ilayer)=v(ilayer)+dtphys*dv(ilayer)
C            temp(ilayer)=temp(ilayer)+dtphys*dtemp(ilayer)
C         ENDDO

c       compute pressure for next time step 
c       ----------------------------------------------------------

C            psurf=psurf+dtphys*dpsurf   ! surface pressure change
           DO ilevel=1,nlevel
             plev(ilevel)=ap(ilevel)+psurf*bp(ilevel)
           ENDDO
           DO ilayer=1,nlayer
             play(ilayer)=aps(ilayer)+psurf*bps(ilayer)
           ENDDO

          date = idt*(dtphys/daysec) 

!      Plume/perturbation scenarios engaged here
        IF ( plume ) THEN 
           
           IF ( idt == ndt_plume ) THEN

           
           call tracer_release(dq)
           
           plume = .False.
           mcd_profiles = .False. 
           
           ENDIF
        ENDIF 
        
!       increment tracers
        DO iq = 1, nqmx
          DO ilayer=1,nlayer
             q(ilayer,iq)=q(ilayer,iq)+dtphys*dq(ilayer,iq)
          ENDDO
        ENDDO

!     Check if the MCD v5.3 tracer interpolation routine has reached
!     the desired shutdown point 
      IF ( mcd_profiles ) THEN 
          
          IF (  ( date >= mcd_dayfin ) .and. 
     &      ( ABS(time*24. - mcd_ltfin) < 1.e-1 ) )
     &    mcd_profiles = .False.
         
      ENDIF 
      
      
        
        
      ENDDO   ! of idt=1,ndt ! end of time stepping loop

c    ========================================================
c    OUTPUTS
c    ========================================================

c    finalize and close grads files "g1d.dat" and "g1d.ctl"

c        CALL endg1d(1,nlayer,zphi/(g*1000.),ndt)
        CALL endg1d(1,nlayer,zlay/1000.,ndt)

c    ========================================================
      END
 
c***********************************************************************
c***********************************************************************
c     Dummy subroutines used only in 3D, but required to
c     compile testphys1d (to cleanly use writediagfi)

      subroutine gr_fi_dyn(nfield,ngrid,im,jm,pfi,pdyn)

      IMPLICIT NONE

      INTEGER im,jm,ngrid,nfield
      REAL pdyn(im,jm,nfield)
      REAL pfi(ngrid,nfield)
      
      if (ngrid.ne.1) then
        write(*,*) "gr_fi_dyn error: in 1D ngrid should be 1!!!"
        stop
      endif
      
      pdyn(1,1,1:nfield)=pfi(1,1:nfield)
      
      end
 
c***********************************************************************
c***********************************************************************

#include "../dyn3d/disvert.F"
