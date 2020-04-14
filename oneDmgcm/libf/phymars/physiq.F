      SUBROUTINE physiq(
     $            ngrid,nlayer,nq
     $            ,firstcall,lastcall
     $            ,pday,ptime,ptimestep
     $            ,idt,ndt_tlm
     $            ,pplev,pplay,pphi
     $            ,pu,pv,pt,pq
     $            ,pw
     $            ,pdu,pdv,pdt,pdq,pdpsrf,tracerdyn
#ifdef MESOSCALE
#include "meso_inc/meso_inc_invar.F"
#endif
     $            )

      USE TLMvars
      
      IMPLICIT NONE
c=======================================================================
c
c   subject:
c   -------
c
c   Organisation of the physical parametrisations of the LMD 
c   martian atmospheric general circulation model.
c
c   The GCM can be run without or with tracer transport
c   depending on the value of Logical "tracer" in file  "callphys.def"
c   Tracers may be water vapor, ice OR chemical species OR dust particles
c
c   SEE comments in initracer.F about numbering of tracer species...
c
c   It includes:
c
c      1. Initialization:
c      1.1 First call initializations
c      1.2 Initialization for every call to physiq
c      1.2.5 Compute mean mass and cp, R and thermal conduction coeff.
c      2. Compute radiative transfer tendencies
c         (longwave and shortwave) for CO2 and aerosols.
c      3. Gravity wave and subgrid scale topography drag :
c      4. Vertical diffusion (turbulent mixing):
c      5. Convective adjustment
c      6. Condensation and sublimation of carbon dioxide.
c      7.  TRACERS :
c       7a. water and water ice
c       7b. call for photochemistry when tracers are chemical species
c       7c. other scheme for tracer (dust) transport (lifting, sedimentation)
c       7d. updates (CO2 pressure variations, surface budget)
c      8. Contribution to tendencies due to thermosphere
c      9. Surface and sub-surface temperature calculations
c     10. Write outputs :
c           - "startfi", "histfi" (if it's time)
c           - Saving statistics (if "callstats = .true.")
c           - Dumping eof (if "calleofdump = .true.")
c           - Output any needed variables in "diagfi" 
c     11. Diagnostic: mass conservation of tracers
c 
c   author: 
c   ------- 
c           Frederic Hourdin	15/10/93
c           Francois Forget		1994
c           Christophe Hourdin	02/1997 
c           Subroutine completly rewritten by F.Forget (01/2000)
c           Introduction of the photochemical module: S. Lebonnois (11/2002)
c           Introduction of the thermosphere module: M. Angelats i Coll (2002)
c           Water ice clouds: Franck Montmessin (update 06/2003)
c           Radiatively active tracers: J.-B. Madeleine (10/2008-06/2009)
c             Nb: See callradite.F for more information.
c           Mesoscale lines: Aymeric Spiga (2007 - 2011) -- check MESOSCALE flags
c           jul 2011 malv+fgg: Modified calls to NIR heating routine and 15 um cooling parameterization
c           
c   arguments:
c   ----------
c
c   input:
c   ------
c    ecri                  period (in dynamical timestep) to write output
c    ngrid                 Size of the horizontal grid.
c                          All internal loops are performed on that grid.
c    nlayer                Number of vertical layers.
c    nq                    Number of advected fields
c    firstcall             True at the first call
c    lastcall              True at the last call
c    pday                  Number of days counted from the North. Spring
c                          equinoxe.
c    ptime                 Universal time (0<ptime<1): ptime=0.5 at 12:00 UT
c    ptimestep             timestep (s)
c  
c    idt                   Model iteration 
c    ndt_tlm               iteration to engage the TLM routines    
c
c    pplay(ngrid,nlayer)   Pressure at the middle of the layers (Pa)
c    pplev(ngrid,nlayer+1) intermediate pressure levels (pa)
c    pphi(ngrid,nlayer)    Geopotential at the middle of the layers (m2s-2)
c    pu(ngrid,nlayer)      u component of the wind (ms-1)
c    pv(ngrid,nlayer)      v component of the wind (ms-1)
c    pt(ngrid,nlayer)      Temperature (K)
c    pq(ngrid,nlayer,nq)   Advected fields
c    pudyn(ngrid,nlayer)    | 
c    pvdyn(ngrid,nlayer)    | Dynamical temporal derivative for the
c    ptdyn(ngrid,nlayer)    | corresponding variables
c    pqdyn(ngrid,nlayer,nq) |
c    pw(ngrid,?)           vertical velocity
c
c   output:
c   -------
c
c    pdu(ngrid,nlayermx)       |
c    pdv(ngrid,nlayermx)       |  Temporal derivative of the corresponding
c    pdt(ngrid,nlayermx)       |  variables due to physical processes.
c    pdq(ngrid,nlayermx,nqmx)  |
c    pdpsrf(ngrid)             |
c    tracerdyn                 call tracer in dynamical part of GCM ?

c
c=======================================================================
c
c    0.  Declarations :
c    ------------------

#include "dimensions.h"
#include "dimphys.h"
#include "comgeomfi.h"
#include "surfdat.h"
#include "comsoil.h"
#include "comdiurn.h"
#include "callkeys.h"
#include "comcstfi.h"
#include "planete.h"
#include "comsaison.h"
#include "control.h"
#include "dimradmars.h"
#include "comg1d.h"
#include "tracer.h"
#include "nlteparams.h"

#include "chimiedata.h"
#include "param.h"
#include "param_v4.h"
#include "conc.h"

#include "netcdf.inc"

#include "slope.h"

#ifdef MESOSCALE
#include "wrf_output_2d.h"
#include "wrf_output_3d.h"
#include "advtrac.h"   !!! this is necessary for tracers (in dyn3d)
#include "meso_inc/meso_inc_var.F"
#endif

c Arguments :
c -----------

c   inputs:
c   -------
      INTEGER ngrid,nlayer,nq
      REAL ptimestep
      REAL pplev(ngridmx,nlayer+1),pplay(ngridmx,nlayer)
      REAL pphi(ngridmx,nlayer)
      REAL pu(ngridmx,nlayer),pv(ngridmx,nlayer)
      REAL pt(ngridmx,nlayer),pq(ngridmx,nlayer,nq)
      REAL pw(ngridmx,nlayer) !Mars pvervel transmit par dyn3d
      REAL zh(ngridmx,nlayermx)      ! potential temperature (K)
      LOGICAL firstcall,lastcall

      REAL pday
      REAL ptime 

      INTEGER idt, ndt_tlm 

      logical tracerdyn

c   outputs:
c   --------
c     physical tendencies
      REAL pdu(ngridmx,nlayer),pdv(ngridmx,nlayer)
      REAL pdt(ngridmx,nlayer),pdq(ngridmx,nlayer,nq)
      REAL pdpsrf(ngridmx) ! surface pressure tendency


c Local saved variables:
c ----------------------
c     aerosol (dust or ice) extinction optical depth  at reference wavelength 
c     "longrefvis" set in dimradmars.h , for one of the "naerkind"  kind of
c      aerosol optical properties  :
      REAL aerosol(ngridmx,nlayermx,naerkind)

      INTEGER day_ini  ! Initial date of the run (sol since Ls=0) 
      INTEGER icount     ! counter of calls to physiq during the run.
      REAL tsurf(ngridmx)            ! Surface temperature (K)
      REAL tsoil(ngridmx,nsoilmx)    ! sub-surface temperatures (K)
      REAL co2ice(ngridmx)           ! co2 ice surface layer (kg.m-2)  
      REAL albedo(ngridmx,2)         ! Surface albedo in each solar band
      REAL emis(ngridmx)             ! Thermal IR surface emissivity
      REAL dtrad(ngridmx,nlayermx)   ! Net atm. radiative heating rate (K.s-1)
      REAL fluxrad_sky(ngridmx)      ! rad. flux from sky absorbed by surface (W.m-2)
      REAL fluxrad(ngridmx)          ! Net radiative surface flux (W.m-2)
      REAL capcal(ngridmx)           ! surface heat capacity (J m-2 K-1)
      REAL fluxgrd(ngridmx)          ! surface conduction flux (W.m-2)
      REAL qsurf(ngridmx,nqmx)       ! tracer on surface (e.g. kg.m-2)
      REAL q2(ngridmx,nlayermx+1)    ! Turbulent Kinetic Energy
      
c     Variables used by the water ice microphysical scheme:
      REAL rice(ngridmx,nlayermx)    ! Water ice geometric mean radius (m)
      REAL nuice(ngridmx,nlayermx)   ! Estimated effective variance
                                     !   of the size distribution
      real rsedcloud(ngridmx,nlayermx) ! Cloud sedimentation radius
      real rhocloud(ngridmx,nlayermx)  ! Cloud density (kg.m-3)
      REAL surfdust(ngridmx,nlayermx) ! dust surface area (m2/m3, if photochemistry)
      REAL surfice(ngridmx,nlayermx)  !  ice surface area (m2/m3, if photochemistry)
      REAL inertiesoil(ngridmx,nsoilmx)! Time varying subsurface
                                       ! thermal inertia (J.s-1/2.m-2.K-1)
                                       ! (used only when tifeedback=.true.)

c     Variables used by the slope model
      REAL sl_ls, sl_lct, sl_lat
      REAL sl_tau, sl_alb, sl_the, sl_psi
      REAL sl_fl0, sl_flu
      REAL sl_ra, sl_di0
      REAL sky

      SAVE day_ini, icount
      SAVE aerosol, tsurf,tsoil
      SAVE co2ice,albedo,emis, q2
      SAVE capcal,fluxgrd,dtrad,fluxrad,fluxrad_sky,qsurf

      REAL stephan   
      DATA stephan/5.67e-08/  ! Stephan Boltzman constant
      SAVE stephan

c Local variables :
c -----------------

      REAL CBRT
      EXTERNAL CBRT

      CHARACTER*80 fichier 
      INTEGER l,ig,ierr,igout,iq,i, tapphys

      REAL fluxsurf_lw(ngridmx)      !incident LW (IR) surface flux (W.m-2)
      REAL fluxsurf_sw(ngridmx,2)    !incident SW (solar) surface flux (W.m-2)
      REAL fluxtop_lw(ngridmx)       !Outgoing LW (IR) flux to space (W.m-2)
      REAL fluxtop_sw(ngridmx,2)     !Outgoing SW (solar) flux to space (W.m-2)
      REAL tauref(ngridmx)           ! Reference column optical depth at odpref
      real,parameter :: odpref=610. ! DOD reference pressure (Pa)
      REAL tau(ngridmx,naerkind)     ! Column dust optical depth at each point
      REAL zls                       !  solar longitude (rad)
      REAL zday                      ! date (time since Ls=0, in martian days)
      REAL zzlay(ngridmx,nlayermx)   ! altitude at the middle of the layers
      REAL zzlev(ngridmx,nlayermx+1) ! altitude at layer boundaries
      REAL latvl1,lonvl1             ! Viking Lander 1 point (for diagnostic)

c     Tendancies due to various processes:
      REAL dqsurf(ngridmx,nqmx)
      REAL zdtlw(ngridmx,nlayermx)     ! (K/s)
      REAL zdtsw(ngridmx,nlayermx)     ! (K/s)
      REAL cldtlw(ngridmx,nlayermx)     ! (K/s) LW heating rate for clear area
      REAL cldtsw(ngridmx,nlayermx)     ! (K/s) SW heating rate for clear area
      REAL zdtnirco2(ngridmx,nlayermx) ! (K/s)
      REAL zdtnlte(ngridmx,nlayermx)   ! (K/s)
      REAL zdtsurf(ngridmx)            ! (K/s)
      REAL zdtcloud(ngridmx,nlayermx)
      REAL zdvdif(ngridmx,nlayermx),zdudif(ngridmx,nlayermx)  ! (m.s-2)
      REAL zdhdif(ngridmx,nlayermx), zdtsdif(ngridmx)         ! (K/s)
      REAL zdvadj(ngridmx,nlayermx),zduadj(ngridmx,nlayermx)  ! (m.s-2)
      REAL zdhadj(ngridmx,nlayermx)                           ! (K/s)
      REAL zdtgw(ngridmx,nlayermx)                            ! (K/s)
      REAL zdugw(ngridmx,nlayermx),zdvgw(ngridmx,nlayermx)    ! (m.s-2)
      REAL zdtc(ngridmx,nlayermx),zdtsurfc(ngridmx)
      REAL zdvc(ngridmx,nlayermx),zduc(ngridmx,nlayermx)

      REAL zdqdif(ngridmx,nlayermx,nqmx), zdqsdif(ngridmx,nqmx)
      REAL zdqsed(ngridmx,nlayermx,nqmx), zdqssed(ngridmx,nqmx)
      REAL zdqdev(ngridmx,nlayermx,nqmx), zdqsdev(ngridmx,nqmx)
      REAL zdqadj(ngridmx,nlayermx,nqmx)
      REAL zdqc(ngridmx,nlayermx,nqmx)
      REAL zdqcloud(ngridmx,nlayermx,nqmx)
      REAL zdqscloud(ngridmx,nqmx)
      REAL zdqchim(ngridmx,nlayermx,nqmx)
      REAL zdqschim(ngridmx,nqmx)

      REAL zdteuv(ngridmx,nlayermx)    ! (K/s)
      REAL zdtconduc(ngridmx,nlayermx) ! (K/s)
      REAL zdumolvis(ngridmx,nlayermx)
      REAL zdvmolvis(ngridmx,nlayermx)
      real zdqmoldiff(ngridmx,nlayermx,nqmx)

c     Local variable for local intermediate calcul:
      REAL zflubid(ngridmx)
      REAL zplanck(ngridmx),zpopsk(ngridmx,nlayermx)
      REAL zdum1(ngridmx,nlayermx)
      REAL zdum2(ngridmx,nlayermx)
      REAL ztim1,ztim2,ztim3, z1,z2
      REAL ztime_fin
      REAL zdh(ngridmx,nlayermx)
      INTEGER length
      PARAMETER (length=100)

c local variables only used for diagnostic (output in file "diagfi" or "stats")
c -----------------------------------------------------------------------------
      REAL ps(ngridmx), zt(ngridmx,nlayermx)
      REAL zu(ngridmx,nlayermx),zv(ngridmx,nlayermx)
      REAL zq(ngridmx,nlayermx,nqmx)
      REAL fluxtop_sw_tot(ngridmx), fluxsurf_sw_tot(ngridmx)
      character*2 str2
      character*5 str5
      real zdtdif(ngridmx,nlayermx), zdtadj(ngridmx,nlayermx)
      REAL tauscaling(ngridmx)   ! Convertion factor for qdust and Ndust
      SAVE tauscaling            ! in case iradia NE 1
      real rdust(ngridmx,nlayermx) ! dust geometric mean radius (m)
      integer igmin, lmin
      logical tdiag

      real co2col(ngridmx)        ! CO2 column
      REAL zplev(ngrid,nlayermx+1),zplay(ngrid,nlayermx)
      REAL zstress(ngrid), cd
      real hco2(nqmx),tmean, zlocal(nlayermx)
      real rho(ngridmx,nlayermx)  ! density
      real vmr(ngridmx,nlayermx)  ! volume mixing ratio
      real rhopart(ngridmx,nlayermx) ! number density of a given species
      real colden(ngridmx,nqmx)   ! vertical column of tracers
      REAL mtot(ngridmx)          ! Total mass of water vapor (kg/m2)
      REAL icetot(ngridmx)        ! Total mass of water ice (kg/m2)
      REAL Nccntot(ngridmx)       ! Total number of ccn (nbr/m2)
      REAL Mccntot(ngridmx)       ! Total mass of ccn (kg/m2)
      REAL rave(ngridmx)          ! Mean water ice effective radius (m)
      REAL opTES(ngridmx,nlayermx)! abs optical depth at 825 cm-1
      REAL tauTES(ngridmx)        ! column optical depth at 825 cm-1
      REAL Qabsice                ! Water ice absorption coefficient
      REAL taucloudtes(ngridmx)! Cloud opacity at infrared
                               !   reference wavelength using
                               !   Qabs instead of Qext
                               !   (direct comparison with TES)
                               
      REAL dqdustsurf(ngridmx) ! surface q dust flux (kg/m2/s)
      REAL dndustsurf(ngridmx) ! surface n dust flux (number/m2/s)
      REAL ndust(ngridmx,nlayermx) ! true n dust (kg/kg)
      REAL qdust(ngridmx,nlayermx) ! true q dust (kg/kg)
      REAL nccn(ngridmx,nlayermx)  ! true n ccn (kg/kg)
      REAL qccn(ngridmx,nlayermx)  ! true q ccn (kg/kg)

c Test 1d/3d scavenging
      real h2otot(ngridmx)
      REAL satu(ngridmx,nlayermx)  ! satu ratio for output
      REAL zqsat(ngridmx,nlayermx) ! saturation

      REAL time_phys

! Added for new NLTE scheme

      real co2vmr_gcm(ngridmx,nlayermx)
      real n2vmr_gcm(ngridmx,nlayermx)
      real ovmr_gcm(ngridmx,nlayermx)
      real covmr_gcm(ngridmx,nlayermx)


c Variables for PBL
      REAL zz1(ngridmx)
      REAL lmax_th_out(ngridmx),zmax_th(ngridmx)
      REAL, SAVE :: wstar(ngridmx)
      REAL, SAVE :: hfmax_th(ngridmx)
      REAL pdu_th(ngridmx,nlayermx),pdv_th(ngridmx,nlayermx)
      REAL pdt_th(ngridmx,nlayermx),pdq_th(ngridmx,nlayermx,nqmx)
      INTEGER lmax_th(ngridmx),dimout,n_out,n
      CHARACTER(50) zstring
      REAL dtke_th(ngridmx,nlayermx+1)
      REAL zcdv(ngridmx), zcdh(ngridmx)
      REAL, ALLOCATABLE, DIMENSION(:,:) :: T_out
      REAL, ALLOCATABLE, DIMENSION(:,:) :: u_out ! Interpolated teta and u at z_out
      REAL u_out1(ngridmx), T_out1(ngridmx)
      REAL, ALLOCATABLE, DIMENSION(:) :: z_out     ! height of interpolation between z0 and z1 [meters]
      REAL ustar(ngridmx),tstar(ngridmx)  ! friction velocity and friction potential temp
      REAL L_mo(ngridmx),vhf(ngridmx),vvv(ngridmx)
      REAL zu2(ngridmx),sensibFlux(ngridmx)

c=======================================================================

       IF ( tlm_on ) THEN 
          TLM(:,:) = 0.D0
       ENDIF 
c 1. Initialisation:
c -----------------
c  1.1   Initialisation only at first call
c  ---------------------------------------
      IF (firstcall) THEN

c        variables set to 0
c        ~~~~~~~~~~~~~~~~~~
         aerosol(:,:,:)=0
         dtrad(:,:)=0

#ifndef MESOSCALE
         fluxrad(:)=0
         wstar(:)=0.
#else
#include "meso_inc/meso_inc_ini_restart.F"
#endif

c        read startfi 
c        ~~~~~~~~~~~~
#ifndef MESOSCALE
! Read netcdf initial physical parameters.
         CALL phyetat0 ("startfi.nc",0,0,
     &         nsoilmx,nq,
     &         day_ini,time_phys,
     &         tsurf,tsoil,emis,q2,qsurf,co2ice)
#else
#include "meso_inc/meso_inc_ini.F"
#endif

         if (pday.ne.day_ini) then
           write(*,*) "PHYSIQ: ERROR: bad synchronization between ",
     &                "physics and dynamics"
           write(*,*) "dynamics day: ",pday
           write(*,*) "physics day:  ",day_ini
           stop
         endif

         write (*,*) 'In physiq day_ini =', day_ini

c        initialize tracers
c        ~~~~~~~~~~~~~~~~~~
         tracerdyn=tracer
         IF (tracer) THEN
            CALL initracer(qsurf,co2ice)
         ENDIF  ! end tracer

c        Initialize albedo and orbital calculation
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         CALL surfini(ngrid,co2ice,qsurf,albedo)
         CALL iniorbit(aphelie,periheli,year_day,peri_day,obliquit)

c        initialize soil 
c        ~~~~~~~~~~~~~~~
         IF (callsoil) THEN
c           Thermal inertia feedback:
            IF (tifeedback) THEN
                CALL soil_tifeedback(ngrid,nsoilmx,qsurf,inertiesoil)
                CALL soil(ngrid,nsoilmx,firstcall,inertiesoil,
     s             ptimestep,tsurf,tsoil,capcal,fluxgrd)
            ELSE
                CALL soil(ngrid,nsoilmx,firstcall,inertiedat,
     s             ptimestep,tsurf,tsoil,capcal,fluxgrd)
            ENDIF ! of IF (tifeedback)
         ELSE
            PRINT*,
     &     'PHYSIQ WARNING! Thermal conduction in the soil turned off'
            DO ig=1,ngrid
               capcal(ig)=1.e5
               fluxgrd(ig)=0.
            ENDDO
         ENDIF
         icount=1

#ifndef MESOSCALE
c        Initialize thermospheric parameters
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if (callthermos) then
            if(solvarmod.eq.0) call param_read
            if(solvarmod.eq.1) call param_read_e107
         endif
#endif
c        Initialize R and Cp as constant

         if (.not.callthermos .and. .not.photochem) then
                 do l=1,nlayermx
                  do ig=1,ngridmx
                   rnew(ig,l)=r
                   cpnew(ig,l)=cpp
                   mmean(ig,l)=mugaz
                   enddo
                  enddo  
         endif         

         if(callnlte.and.nltemodel.eq.2) call nlte_setup
         if(callnirco2.and.nircorr.eq.1) call NIR_leedat

        IF (tracer.AND.water.AND.(ngridmx.NE.1)) THEN
          write(*,*)"physiq: water_param Surface water ice albedo:", 
     .                  albedo_h2o_ice
        ENDIF
                   
      ENDIF        !  (end of "if firstcall")


c    Initialise the TLM
c     -------------------
       if (tlm_on) then
          TLM(:,:) = 0.
          TLM_photo(:,:) = 0.
       endif   

c ---------------------------------------------------
c 1.2   Initializations done at every physical timestep:
c ---------------------------------------------------
c
      IF (ngrid.NE.ngridmx) THEN
         PRINT*,'STOP in PHYSIQ'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'ngrid     = ',ngrid
         PRINT*,'ngridmx   = ',ngridmx
         STOP
      ENDIF

c     Initialize various variables
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      pdv(:,:)=0
      pdu(:,:)=0
      pdt(:,:)=0
      pdq(:,:,:)=0
      pdpsrf(:)=0
      zflubid(:)=0
      zdtsurf(:)=0
      dqsurf(:,:)=0
      igout=ngrid/2+1 


      zday=pday+ptime ! compute time, in sols (and fraction thereof)

c     Compute Solar Longitude (Ls) :
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (season) then
         call solarlong(zday,zls)
      else
         call solarlong(float(day_ini),zls)
      end if

c     Compute geopotential at interlayers
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     ponderation des altitudes au niveau des couches en dp/p

      DO l=1,nlayer
         DO ig=1,ngrid
            zzlay(ig,l)=pphi(ig,l)/g
         ENDDO
      ENDDO
      DO ig=1,ngrid
         zzlev(ig,1)=0.
         zzlev(ig,nlayer+1)=1.e7    ! dummy top of last layer above 10000 km...
      ENDDO
      DO l=2,nlayer
         DO ig=1,ngrid
            z1=(pplay(ig,l-1)+pplev(ig,l))/(pplay(ig,l-1)-pplev(ig,l))
            z2=(pplev(ig,l)+pplay(ig,l))/(pplev(ig,l)-pplay(ig,l))
            zzlev(ig,l)=(z1*zzlay(ig,l-1)+z2*zzlay(ig,l))/(z1+z2)
         ENDDO
      ENDDO


!     Potential temperature calculation not the same in physiq and dynamic

c     Compute potential temperature
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO l=1,nlayer
         DO ig=1,ngrid 
            zpopsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**rcp
            zh(ig,l)=pt(ig,l)/zpopsk(ig,l)
         ENDDO
      ENDDO

#ifndef MESOSCALE
c-----------------------------------------------------------------------
c    1.2.5 Compute mean mass, cp, and R
c    --------------------------------

      if(photochem.or.callthermos) then
         call concentrations(pplay,pt,pdt,pq,pdq,ptimestep)
      endif
#endif
c----------------------------------------------------------------------
c    2. Compute radiative tendencies :
c------------------------------------


      IF (callrad) THEN
         IF( MOD(icount-1,iradia).EQ.0) THEN

c          Local Solar zenith angle
c          ~~~~~~~~~~~~~~~~~~~~~~~~
           CALL orbite(zls,dist_sol,declin)

           IF(diurnal) THEN
               ztim1=SIN(declin)
               ztim2=COS(declin)*COS(2.*pi*(zday-.5))
               ztim3=-COS(declin)*SIN(2.*pi*(zday-.5))

               CALL solang(ngrid,sinlon,coslon,sinlat,coslat,
     s         ztim1,ztim2,ztim3, mu0,fract)

           ELSE
               CALL mucorr(ngrid,declin, lati, mu0, fract,10000.,rad)
           ENDIF

c          NLTE cooling from CO2 emission
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           IF(callnlte) then
              if(nltemodel.eq.0.or.nltemodel.eq.1) then
                 CALL nltecool(ngrid,nlayer,nq,pplay,pt,pq,zdtnlte)
              else if(nltemodel.eq.2) then
                co2vmr_gcm(1:ngrid,1:nlayer)=
     &                      pq(1:ngrid,1:nlayer,igcm_co2)*
     &                      mmean(1:ngrid,1:nlayer)/mmol(igcm_co2)
                n2vmr_gcm(1:ngrid,1:nlayer)=
     &                      pq(1:ngrid,1:nlayer,igcm_n2)*
     &                      mmean(1:ngrid,1:nlayer)/mmol(igcm_n2)
                covmr_gcm(1:ngrid,1:nlayer)=
     &                      pq(1:ngrid,1:nlayer,igcm_co)*
     &                      mmean(1:ngrid,1:nlayer)/mmol(igcm_co)
                ovmr_gcm(1:ngrid,1:nlayer)=
     &                      pq(1:ngrid,1:nlayer,igcm_o)*
     &                      mmean(1:ngrid,1:nlayer)/mmol(igcm_o)
                 
                 CALL nlte_tcool(ngrid,nlayer,pplay*9.869e-6,
     $                pt,zzlay,co2vmr_gcm, n2vmr_gcm, covmr_gcm, 
     $                ovmr_gcm,  zdtnlte )

                 zdtnlte(1:ngrid,1:nlayer)=
     &                             zdtnlte(1:ngrid,1:nlayer)/86400.
              endif
           else
              zdtnlte(:,:)=0.
           endif

c          Find number of layers for LTE radiation calculations
           IF(MOD(iphysiq*(icount-1),day_step).EQ.0)
     &          CALL nlthermeq(ngrid,nlayer,pplev,pplay)

c          Note: Dustopacity.F has been transferred to callradite.F
         
c          Call main radiative transfer scheme
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c          Transfer through CO2 (except NIR CO2 absorption)
c            and aerosols (dust and water ice)

c          Radiative transfer
c          ------------------

           CALL callradite(icount,ngrid,nlayer,nq,zday,zls,pq,albedo,
     $     emis,mu0,pplev,pplay,pt,tsurf,fract,dist_sol,igout,
     $     zdtlw,zdtsw,fluxsurf_lw,fluxsurf_sw,fluxtop_lw,fluxtop_sw,
     $     tauref,tau,aerosol,tauscaling,taucloudtes,rdust,rice,
     $     nuice,co2ice)

c          Outputs for basic check (middle of domain)
c          ------------------------------------------
           write(*,'("Ls =",f11.6," check lat =",f10.6,
     &               " lon =",f11.6)')
     &           zls*180./pi,lati(igout)*180/pi,long(igout)*180/pi
           write(*,'(" tauref(",f4.0," Pa) =",f9.6,
     &             " tau(",f4.0," Pa) =",f9.6)')
     &            odpref,tauref(igout),
     &            odpref,tau(igout,1)*odpref/pplev(igout,1)
c          ---------------------------------------------------------
c          Call slope parameterization for direct and scattered flux
c          ---------------------------------------------------------
           IF(callslope) THEN
            DO ig=1,ngrid  
              sl_the = theta_sl(ig)
              IF (sl_the .ne. 0.) THEN
                ztim1=fluxsurf_sw(ig,1)+fluxsurf_sw(ig,2)
                DO l=1,2
                 sl_lct = ptime*24. + 180.*long(ig)/pi/15.
                 sl_ra  = pi*(1.0-sl_lct/12.)
                 sl_lat = 180.*lati(ig)/pi
                 sl_tau = tau(ig,1) !il faudrait iaerdust(iaer)
                 sl_alb = albedo(ig,l)
                 sl_psi = psi_sl(ig)
                 sl_fl0 = fluxsurf_sw(ig,l)
                 sl_di0 = 0.
                 if (mu0(ig) .gt. 0.) then
                  sl_di0 = mu0(ig)*(exp(-sl_tau/mu0(ig)))
                  sl_di0 = sl_di0*1370./dist_sol/dist_sol
                  sl_di0 = sl_di0/ztim1
                  sl_di0 = fluxsurf_sw(ig,l)*sl_di0
                 endif
                 ! you never know (roundup concern...)
                 if (sl_fl0 .lt. sl_di0) sl_di0=sl_fl0
                 !!!!!!!!!!!!!!!!!!!!!!!!!!
                 CALL param_slope( mu0(ig), declin, sl_ra, sl_lat, 
     &                             sl_tau, sl_alb, sl_the, sl_psi,
     &                             sl_di0, sl_fl0, sl_flu )
                 !!!!!!!!!!!!!!!!!!!!!!!!!!
                 fluxsurf_sw(ig,l) = sl_flu
                ENDDO
              !!! compute correction on IR flux as well
              sky= (1.+cos(pi*theta_sl(ig)/180.))/2.
              fluxsurf_lw(ig)= fluxsurf_lw(ig)*sky
              ENDIF
            ENDDO
           ENDIF

c          CO2 near infrared absorption
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           zdtnirco2(:,:)=0
           if (callnirco2) then
              call nirco2abs (ngrid,nlayer,pplay,dist_sol,nq,pq,
     .                       mu0,fract,declin, zdtnirco2)
           endif

c          Radiative flux from the sky absorbed by the surface (W.m-2)
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           DO ig=1,ngrid
               fluxrad_sky(ig)=emis(ig)*fluxsurf_lw(ig)
     $         +fluxsurf_sw(ig,1)*(1.-albedo(ig,1))
     $         +fluxsurf_sw(ig,2)*(1.-albedo(ig,2))
           ENDDO


c          Net atmospheric radiative heating rate (K.s-1)
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           IF(callnlte) THEN
              CALL blendrad(ngrid, nlayer, pplay,
     &             zdtsw, zdtlw, zdtnirco2, zdtnlte, dtrad)
           ELSE
              DO l=1,nlayer
                 DO ig=1,ngrid
                    dtrad(ig,l)=zdtsw(ig,l)+zdtlw(ig,l)
     &                          +zdtnirco2(ig,l)
                  ENDDO
              ENDDO
           ENDIF

        ENDIF ! of if(mod(icount-1,iradia).eq.0)

c    Transformation of the radiative tendencies:
c    -------------------------------------------

c          Net radiative surface flux (W.m-2)
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~
c
           DO ig=1,ngrid
               zplanck(ig)=tsurf(ig)*tsurf(ig)
               zplanck(ig)=emis(ig)*
     $         stephan*zplanck(ig)*zplanck(ig)
               fluxrad(ig)=fluxrad_sky(ig)-zplanck(ig)
               IF(callslope) THEN
                 sky= (1.+cos(pi*theta_sl(ig)/180.))/2.
                 fluxrad(ig)=fluxrad(ig)+(1.-sky)*zplanck(ig)
               ENDIF
           ENDDO

         DO l=1,nlayer
            DO ig=1,ngrid
               pdt(ig,l)=pdt(ig,l)+dtrad(ig,l)
            ENDDO
         ENDDO

      ENDIF ! of IF (callrad)
      
c-----------------------------------------------------------------------
c    3. Gravity wave and subgrid scale topography drag :
c    -------------------------------------------------


      IF(calllott)THEN

        CALL calldrag_noro(ngrid,nlayer,ptimestep,
     &                 pplay,pplev,pt,pu,pv,zdtgw,zdugw,zdvgw)
 
        DO l=1,nlayer
          DO ig=1,ngrid
            pdv(ig,l)=pdv(ig,l)+zdvgw(ig,l)
            pdu(ig,l)=pdu(ig,l)+zdugw(ig,l)
            pdt(ig,l)=pdt(ig,l)+zdtgw(ig,l)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c    4. Vertical diffusion (turbulent mixing):
c    -----------------------------------------

      IF (calldifv) THEN

         DO ig=1,ngrid
            zflubid(ig)=fluxrad(ig)+fluxgrd(ig)
         ENDDO

         zdum1(:,:)=0
         zdum2(:,:)=0
         do l=1,nlayer
            do ig=1,ngrid
               zdh(ig,l)=pdt(ig,l)/zpopsk(ig,l)
            enddo
         enddo


#ifdef MESOSCALE
      IF (.not.flag_LES) THEN
#endif
c ----------------------
c Treatment of a special case : using new surface layer (Richardson based)
c without using the thermals in gcm and mesoscale can yield problems in
c weakly unstable situations when winds are near to 0. For those cases, we add
c a unit subgrid gustiness. Remember that thermals should be used we using the 
c Richardson based surface layer model.
        IF ( .not.calltherm .and. callrichsl ) THEN
          DO ig=1, ngridmx
             IF (zh(ig,1) .lt. tsurf(ig)) THEN
               wstar(ig)=1.
               hfmax_th(ig)=0.2
             ELSE
               wstar(ig)=0.
               hfmax_th(ig)=0.
             ENDIF      
          ENDDO
        ENDIF
c ----------------------
#ifdef MESOSCALE
      ENDIF
#endif

         IF (tke_heat_flux .ne. 0.) THEN
             zz1(:)=(pt(:,1)+pdt(:,1)*ptimestep)*(r/g)*
     &                 (-alog(pplay(:,1)/pplev(:,1)))
             pdt(:,1)=pdt(:,1) + (tke_heat_flux/zz1(:))*zpopsk(:,1)
         ENDIF

c        Calling vdif (Martian version WITH CO2 condensation)
         CALL vdifc(ngrid,nlayer,nq,co2ice,zpopsk,
     $        ptimestep,capcal,lwrite,
     $        pplay,pplev,zzlay,zzlev,z0,
     $        pu,pv,zh,pq,tsurf,emis,qsurf,
     $        zdum1,zdum2,zdh,pdq,zflubid,
     $        zdudif,zdvdif,zdhdif,zdtsdif,q2,
     &        zdqdif,zdqsdif,wstar,zcdv,zcdh,hfmax_th,sensibFlux
#ifdef MESOSCALE
     &        ,flag_LES
#endif
     &        )


#ifdef MESOSCALE
#include "meso_inc/meso_inc_les.F"
#endif
         DO l=1,nlayer
            DO ig=1,ngrid
               pdv(ig,l)=pdv(ig,l)+zdvdif(ig,l)
               pdu(ig,l)=pdu(ig,l)+zdudif(ig,l)
               pdt(ig,l)=pdt(ig,l)+zdhdif(ig,l)*zpopsk(ig,l)

               zdtdif(ig,l)=zdhdif(ig,l)*zpopsk(ig,l) ! for diagnostic only

            ENDDO
         ENDDO

          DO ig=1,ngrid
             zdtsurf(ig)=zdtsurf(ig)+zdtsdif(ig)
          ENDDO

         if (tracer) then
#ifdef MESOSCALE
         IF (.not.flag_LES) THEN
#endif 
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                 pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqdif(ig,l,iq) 
              ENDDO
            ENDDO
           ENDDO
           DO iq=1, nq
              DO ig=1,ngrid
                 dqsurf(ig,iq)=dqsurf(ig,iq) + zdqsdif(ig,iq)
              ENDDO
           ENDDO
#ifdef MESOSCALE
         ENDIF
#endif 
         end if ! of if (tracer)

      ELSE    
         DO ig=1,ngrid
            zdtsurf(ig)=zdtsurf(ig)+
     s      (fluxrad(ig)+fluxgrd(ig))/capcal(ig)
         ENDDO
#ifdef MESOSCALE
         IF (flag_LES) THEN
            write(*,*) 'LES mode !' 
            write(*,*) 'Please set calldifv to T in callphys.def'
            STOP
         ENDIF
#endif
      ENDIF ! of IF (calldifv)
      

      

c-----------------------------------------------------------------------
c   5. Thermals :
c   -----------------------------

      if(calltherm) then
 
        call calltherm_interface(firstcall,
     $ zzlev,zzlay,
     $ ptimestep,pu,pv,pt,pq,pdu,pdv,pdt,pdq,q2,
     $ pplay,pplev,pphi,zpopsk,
     $ pdu_th,pdv_th,pdt_th,pdq_th,lmax_th,zmax_th,
     $ dtke_th,zdhdif,hfmax_th,wstar,sensibFlux)
     
         DO l=1,nlayer
           DO ig=1,ngrid
              pdu(ig,l)=pdu(ig,l)+pdu_th(ig,l)
              pdv(ig,l)=pdv(ig,l)+pdv_th(ig,l)
              pdt(ig,l)=pdt(ig,l)+pdt_th(ig,l)
              q2(ig,l)=q2(ig,l)+dtke_th(ig,l)*ptimestep
           ENDDO
        ENDDO
 
        DO ig=1,ngrid
          q2(ig,nlayer+1)=q2(ig,nlayer+1)+dtke_th(ig,nlayer+1)*ptimestep
        ENDDO      
    
        if (tracer) then
        DO iq=1,nq
         DO l=1,nlayer
           DO ig=1,ngrid
             pdq(ig,l,iq)=pdq(ig,l,iq)+pdq_th(ig,l,iq)
           ENDDO
         ENDDO
        ENDDO
        endif

        lmax_th_out(:)=real(lmax_th(:))

        else   !of if calltherm
        lmax_th(:)=0
        wstar(:)=0.
        hfmax_th(:)=0.
        lmax_th_out(:)=0.
        end if

c-----------------------------------------------------------------------
c   5. Dry convective adjustment:
c   -----------------------------

      IF(calladj) THEN

         DO l=1,nlayer
            DO ig=1,ngrid
               zdh(ig,l)=pdt(ig,l)/zpopsk(ig,l)
            ENDDO
         ENDDO
         zduadj(:,:)=0
         zdvadj(:,:)=0
         zdhadj(:,:)=0

         CALL convadj(ngrid,nlayer,nq,ptimestep,
     $                pplay,pplev,zpopsk,lmax_th,
     $                pu,pv,zh,pq,
     $                pdu,pdv,zdh,pdq,
     $                zduadj,zdvadj,zdhadj,
     $                zdqadj)


         DO l=1,nlayer
            DO ig=1,ngrid
               pdu(ig,l)=pdu(ig,l)+zduadj(ig,l)
               pdv(ig,l)=pdv(ig,l)+zdvadj(ig,l)
               pdt(ig,l)=pdt(ig,l)+zdhadj(ig,l)*zpopsk(ig,l)

               zdtadj(ig,l)=zdhadj(ig,l)*zpopsk(ig,l) ! for diagnostic only
            ENDDO
         ENDDO

         if(tracer) then 
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                 pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqadj(ig,l,iq) 
              ENDDO
            ENDDO
           ENDDO
         end if
      ENDIF ! of IF(calladj)


      
c-----------------------------------------------------------------------
c   6. Carbon dioxide condensation-sublimation:
c   -------------------------------------------

      IF (tituscap) THEN
         !!! get the actual co2 seasonal cap from Titus observations
         CALL geticecover( ngrid, 180.*zls/pi,
     .                  180.*long/pi, 180.*lati/pi, co2ice )
         co2ice = co2ice * 10000.
      ENDIF

      IF (callcond) THEN
         CALL newcondens(ngrid,nlayer,nq,ptimestep,
     $              capcal,pplay,pplev,tsurf,pt,
     $              pphi,pdt,pdu,pdv,zdtsurf,pu,pv,pq,pdq,
     $              co2ice,albedo,emis,
     $              zdtc,zdtsurfc,pdpsrf,zduc,zdvc,zdqc,
     $              fluxsurf_sw,zls) 

         DO l=1,nlayer
           DO ig=1,ngrid
             pdt(ig,l)=pdt(ig,l)+zdtc(ig,l)
             pdv(ig,l)=pdv(ig,l)+zdvc(ig,l)
             pdu(ig,l)=pdu(ig,l)+zduc(ig,l)
           ENDDO
         ENDDO
         DO ig=1,ngrid
            zdtsurf(ig) = zdtsurf(ig) + zdtsurfc(ig)
         ENDDO

         IF (tracer) THEN
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqc(ig,l,iq) 
              ENDDO
            ENDDO
           ENDDO
         ENDIF ! of IF (tracer)

      ENDIF  ! of IF (callcond)
      
c-----------------------------------------------------------------------
c   7. Specific parameterizations for tracers 
c:   -----------------------------------------

      if (tracer) then 

c   7a. Water and ice
c     ---------------

c        ---------------------------------------
c        Water ice condensation in the atmosphere
c        ----------------------------------------
         IF (water) THEN

           call watercloud(ngrid,nlayer,ptimestep,
     &                pplev,pplay,pdpsrf,zzlay, pt,pdt,
     &                pq,pdq,zdqcloud,zdtcloud,
     &                nq,tau,tauscaling,rdust,rice,nuice,
     &                rsedcloud,rhocloud)

c Temperature variation due to latent heat release
           if (activice) then
               pdt(1:ngrid,1:nlayer) = 
     &          pdt(1:ngrid,1:nlayer) + 
     &          zdtcloud(1:ngrid,1:nlayer)
           endif
           
! increment water vapour and ice atmospheric tracers tendencies
           pdq(1:ngrid,1:nlayer,igcm_h2o_vap) = 
     &       pdq(1:ngrid,1:nlayer,igcm_h2o_vap) + 
     &       zdqcloud(1:ngrid,1:nlayer,igcm_h2o_vap)
           pdq(1:ngrid,1:nlayer,igcm_h2o_ice) = 
     &       pdq(1:ngrid,1:nlayer,igcm_h2o_ice) + 
     &       zdqcloud(1:ngrid,1:nlayer,igcm_h2o_ice)

! increment dust and ccn masses and numbers
           if (microphys) then
             pdq(1:ngrid,1:nlayer,igcm_dust_mass) = 
     &         pdq(1:ngrid,1:nlayer,igcm_dust_mass) + 
     &         zdqcloud(1:ngrid,1:nlayer,igcm_dust_mass)
             pdq(1:ngrid,1:nlayer,igcm_dust_number) = 
     &         pdq(1:ngrid,1:nlayer,igcm_dust_number) + 
     &         zdqcloud(1:ngrid,1:nlayer,igcm_dust_number)
             pdq(1:ngrid,1:nlayer,igcm_ccn_mass) = 
     &         pdq(1:ngrid,1:nlayer,igcm_ccn_mass) + 
     &         zdqcloud(1:ngrid,1:nlayer,igcm_ccn_mass)
             pdq(1:ngrid,1:nlayer,igcm_ccn_number) = 
     &         pdq(1:ngrid,1:nlayer,igcm_ccn_number) + 
     &         zdqcloud(1:ngrid,1:nlayer,igcm_ccn_number)
           endif

! We need to check that we have Nccn & Ndust > 0
! This is due to single precision rounding problems
         if (scavenging) then

           where (pq(:,:,igcm_ccn_mass) + 
     &      ptimestep*pdq(:,:,igcm_ccn_mass) < 0.)
              pdq(:,:,igcm_ccn_mass) = 
     &        - pq(:,:,igcm_ccn_mass)/ptimestep + 1.e-30
              pdq(:,:,igcm_ccn_number) = 
     &        - pq(:,:,igcm_ccn_number)/ptimestep + 1.e-30
           end where
           where (pq(:,:,igcm_ccn_number) + 
     &      ptimestep*pdq(:,:,igcm_ccn_number) < 0.)
              pdq(:,:,igcm_ccn_mass) = 
     &        - pq(:,:,igcm_ccn_mass)/ptimestep + 1.e-30
              pdq(:,:,igcm_ccn_number) = 
     &        - pq(:,:,igcm_ccn_number)/ptimestep + 1.e-30
           end where
           where (pq(:,:,igcm_dust_mass) + 
     &      ptimestep*pdq(:,:,igcm_dust_mass) < 0.)
              pdq(:,:,igcm_dust_mass) = 
     &        - pq(:,:,igcm_dust_mass)/ptimestep + 1.e-30
              pdq(:,:,igcm_dust_number) = 
     &        - pq(:,:,igcm_dust_number)/ptimestep + 1.e-30
           end where
           where (pq(:,:,igcm_dust_number) + 
     &      ptimestep*pdq(:,:,igcm_dust_number) < 0.)
              pdq(:,:,igcm_dust_mass) = 
     &        - pq(:,:,igcm_dust_mass)/ptimestep + 1.e-30
              pdq(:,:,igcm_dust_number) = 
     &        - pq(:,:,igcm_dust_number)/ptimestep + 1.e-30
           end where
           
         endif ! of if scavenging
                      

         END IF  ! of IF (water)

c   7b. Aerosol particles
c     -------------------

c        ----------
c        Dust devil :
c        ----------
         IF(callddevil) then 
           call dustdevil(ngrid,nlayer,nq, pplev,pu,pv,pt, tsurf,q2,
     &                zdqdev,zdqsdev)
 
           if (dustbin.ge.1) then
              do iq=1,nq
                 DO l=1,nlayer
                    DO ig=1,ngrid
                       pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqdev(ig,l,iq)
                    ENDDO
                 ENDDO
              enddo
              do iq=1,nq
                 DO ig=1,ngrid
                    dqsurf(ig,iq)= dqsurf(ig,iq) + zdqsdev(ig,iq)
                 ENDDO
              enddo
           endif  ! of if (dustbin.ge.1)

         END IF ! of IF (callddevil)

c        ------------- 
c        Sedimentation :   acts also on water ice
c        ------------- 
         IF (sedimentation) THEN 
           !call zerophys(ngrid*nlayer*nq, zdqsed)
           zdqsed(1:ngrid,1:nlayer,1:nq)=0
           !call zerophys(ngrid*nq, zdqssed)
           zdqssed(1:ngrid,1:nq)=0

           call callsedim(ngrid,nlayer, ptimestep,
     &                pplev,zzlev, zzlay, pt, rdust, rice,
     &                rsedcloud,rhocloud,
     &                pq, pdq, zdqsed, zdqssed,nq, 
     &                tau,tauscaling)

           DO iq=1, nq
             DO l=1,nlayer
               DO ig=1,ngrid
                    pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqsed(ig,l,iq)
               ENDDO
             ENDDO
           ENDDO
           DO iq=1, nq
             DO ig=1,ngrid
                dqsurf(ig,iq)= dqsurf(ig,iq) + zdqssed(ig,iq)
             ENDDO
           ENDDO
         END IF   ! of IF (sedimentation)
        
c Add lifted dust to tendancies after sedimentation in the LES
#ifdef MESOSCALE
#include "meso_inc/meso_inc_lift_les.F"
#endif

c
c   7c. Chemical species
c     ------------------

#ifndef MESOSCALE
c        --------------
c        photochemistry :
c        --------------
         IF (photochem .or. thermochem) then

!           dust and ice surface area
            call surfacearea(ngrid, nlayer, ptimestep, pplay, zzlay, 
     $                       pt, pq, pdq, nq, 
     $                       rdust, rice, tau, tauscaling, 
     $                       surfdust, surfice)
!           call photochemistry
            call calchim(ptimestep,pplay,pplev,pt,pdt,dist_sol,mu0,
     $                   zzlev,zzlay,zday,pq,pdq,zdqchim,zdqschim,
     $                   zdqcloud,zdqscloud,tauref,co2ice,
     $                   pu,pdu,pv,pdv,surfdust,surfice)

           ! increment values of tracers:
           DO iq=1,nq ! loop on all tracers; tendencies for non-chemistry
                      ! tracers is zero anyways
             DO l=1,nlayer
               DO ig=1,ngrid
                 pdq(ig,l,iq)=pdq(ig,l,iq)+zdqchim(ig,l,iq)
               ENDDO
             ENDDO
           ENDDO ! of DO iq=1,nq
           
           ! add condensation tendency for H2O2
           if (igcm_h2o2.ne.0) then
             DO l=1,nlayer
               DO ig=1,ngrid
                 pdq(ig,l,igcm_h2o2)=pdq(ig,l,igcm_h2o2)
     &                                +zdqcloud(ig,l,igcm_h2o2)
               ENDDO
             ENDDO
           endif

           ! increment surface values of tracers:
           DO iq=1,nq ! loop on all tracers; tendencies for non-chemistry
                      ! tracers is zero anyways
             DO ig=1,ngrid
               dqsurf(ig,iq)=dqsurf(ig,iq)+zdqschim(ig,iq)
             ENDDO
           ENDDO ! of DO iq=1,nq

           ! add condensation tendency for H2O2
           if (igcm_h2o2.ne.0) then
             DO ig=1,ngrid
               dqsurf(ig,igcm_h2o2)=dqsurf(ig,igcm_h2o2)
     &                                +zdqscloud(ig,igcm_h2o2)
             ENDDO
           endif

         END IF  ! of IF (photochem.or.thermochem)
#endif

c   7d. Updates
c     ---------

        DO iq=1, nq
          DO ig=1,ngrid

c       ---------------------------------
c       Updating tracer budget on surface
c       ---------------------------------
            qsurf(ig,iq)=qsurf(ig,iq)+ptimestep*dqsurf(ig,iq)

          ENDDO  ! (ig)
        ENDDO    ! (iq)

      endif !  of if (tracer) 

#ifndef MESOSCALE

c-----------------------------------------------------------------------
c   8. THERMOSPHERE CALCULATION
c-----------------------------------------------------------------------

      if (callthermos) then
        call thermosphere(pplev,pplay,dist_sol,
     $     mu0,ptimestep,ptime,zday,tsurf,zzlev,zzlay,
     &     pt,pq,pu,pv,pdt,pdq,
     $     zdteuv,zdtconduc,zdumolvis,zdvmolvis,zdqmoldiff)

        DO l=1,nlayer
          DO ig=1,ngrid
            dtrad(ig,l)=dtrad(ig,l)+zdteuv(ig,l)
            pdt(ig,l)=pdt(ig,l)+zdtconduc(ig,l)
     &                         +zdteuv(ig,l)
            pdv(ig,l)=pdv(ig,l)+zdvmolvis(ig,l)
            pdu(ig,l)=pdu(ig,l)+zdumolvis(ig,l)
            DO iq=1, nq
              pdq(ig,l,iq)=pdq(ig,l,iq)+zdqmoldiff(ig,l,iq)
            ENDDO
          ENDDO
        ENDDO

      endif ! of if (callthermos)
#endif
c-----------------------------------------------------------------------
c   9. Surface  and sub-surface soil temperature
c-----------------------------------------------------------------------
c
c
c   9.1 Increment Surface temperature:
c   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DO ig=1,ngrid
         tsurf(ig)=tsurf(ig)+ptimestep*zdtsurf(ig) 
      ENDDO

c  Prescribe a cold trap at south pole (except at high obliquity !!)
c  Temperature at the surface is set there to be the temperature
c  corresponding to equilibrium temperature between phases of CO2


      IF (tracer.AND.water.AND.(ngridmx.NE.1)) THEN
#ifndef MESOSCALE
         if (caps.and.(obliquit.lt.27.)) then
           ! NB: Updated surface pressure, at grid point 'ngrid', is
           !     ps(ngrid)=pplev(ngrid,1)+pdpsrf(ngrid)*ptimestep
           tsurf(ngrid)=1./(1./136.27-r/5.9e+5*alog(0.0095*
     &                     (pplev(ngrid,1)+pdpsrf(ngrid)*ptimestep)))
         endif
#endif
c       -------------------------------------------------------------
c       Change of surface albedo in case of ground frost
c       everywhere except on the north permanent cap and in regions
c       covered by dry ice. 
c              ALWAYS PLACE these lines after newcondens !!!
c       -------------------------------------------------------------
         do ig=1,ngrid
           if ((co2ice(ig).eq.0).and.
     &        (qsurf(ig,igcm_h2o_ice).gt.frost_albedo_threshold)) then
              albedo(ig,1) = albedo_h2o_ice
              albedo(ig,2) = albedo_h2o_ice
c              write(*,*) "frost thickness", qsurf(ig,igcm_h2o_ice)
c              write(*,*) "physiq.F frost :"
c     &        ,lati(ig)*180./pi, long(ig)*180./pi
           endif
         enddo  ! of do ig=1,ngrid
      ENDIF  ! of IF (tracer.AND.water.AND.(ngridmx.NE.1))

c
c   9.2 Compute soil temperatures and subsurface heat flux:
c   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (callsoil) THEN
c       Thermal inertia feedback
        IF (tifeedback) THEN
         CALL soil_tifeedback(ngrid,nsoilmx,qsurf,inertiesoil)
         CALL soil(ngrid,nsoilmx,.false.,inertiesoil,
     s          ptimestep,tsurf,tsoil,capcal,fluxgrd)
        ELSE
         CALL soil(ngrid,nsoilmx,.false.,inertiedat,
     s          ptimestep,tsurf,tsoil,capcal,fluxgrd)
        ENDIF
      ENDIF
      
c ======================================================================
c 10.0: Construct the 1-D model Jacobian [TLM] and advance perturbation
c       vector [PERTS] to next time-step then call the netCDF save routine
c ======================================================================

          IF (tlm_on) THEN

C                TLM = TLM_trans + TLM_photo
C                TLM = tlm_ident + TLM*ptimestep

               TLM( (t_co2-1)*nlayermx + 1 :
     $              t_co2*nlayermx, : ) = 0.D0


              call TLM_ADJ_out(idt,ptime*24.)
              

          ENDIF
          

c 10.1: Conditional activation of TLM routines
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       IF ( idt == ndt_tlm
     $      .and. call_tlm ) THEN 
            
            call_tlm = .False.     ! Prevents multiple initialisations.
            
            tlm_on = .True.        ! Routines will now engage.
            
            mcd_profiles = .False. ! Interpolation routines will 
                                   ! invalidate TLM results.
            
            calladj = .False. ! Convective adjustments require a 
                               ! conditional that the TLM will not
                               ! be capable of performing.
            callcond = .False. ! Again, conditionals activate this 
                               ! module.
            call TLM_initialise(pq(1,:,:))

       

       ELSEIF ( idt == ndt_tlm 
     $      .and. pertloop ) THEN 

           CALL scalar_perturbation(pq(1,:,:))

            mcd_profiles = .False. ! Interpolation routines will 
                                   ! invalidate TLM results.
            
            calladj = .False. ! Convective adjustments require a 
                               ! conditional that the TLM will not
                               ! be capable of performing.
            callcond = .False. ! Again, conditionals activate this 
                               ! module.
            pertloop = .False.

       ENDIF 












c-----------------------------------------------------------------------
c  10. Write output files
c  ----------------------

c Save variables for eventual restart in MMM and LES
#ifdef MESOSCALE
#include "meso_inc/meso_inc_save_restart.F"
#endif

c    -------------------------------
c    Dynamical fields incrementation
c    -------------------------------
c (FOR OUTPUT ONLY : the actual model integration is performed in the dynamics)
      ! temperature, zonal and meridional wind
      DO l=1,nlayer
        DO ig=1,ngrid
          zt(ig,l)=pt(ig,l)  + pdt(ig,l)*ptimestep
          zu(ig,l)=pu(ig,l)  + pdu(ig,l)*ptimestep
          zv(ig,l)=pv(ig,l)  + pdv(ig,l)*ptimestep
        ENDDO
      ENDDO

      ! tracers
      DO iq=1, nq
        DO l=1,nlayer
          DO ig=1,ngrid
            zq(ig,l,iq)=pq(ig,l,iq) +pdq(ig,l,iq)*ptimestep
          ENDDO
        ENDDO
      ENDDO

      ! surface pressure
      DO ig=1,ngrid
          ps(ig)=pplev(ig,1) + pdpsrf(ig)*ptimestep
      ENDDO

      ! pressure
      DO l=1,nlayer
        DO ig=1,ngrid
             zplev(ig,l)=pplev(ig,l)/pplev(ig,1)*ps(ig)
             zplay(ig,l)=pplay(ig,l)/pplev(ig,1)*ps(ig)
        ENDDO
      ENDDO

      ! Density 
      DO l=1,nlayer
         DO ig=1,ngrid
            rho(ig,l) = zplay(ig,l)/(rnew(ig,l)*zt(ig,l))
         ENDDO
      ENDDO

      ! Potential Temperature

       DO ig=1,ngridmx
          DO l=1,nlayermx
              zh(ig,l) = zt(ig,l)*(zplev(ig,1)/zplay(ig,l))**rcp
          ENDDO
       ENDDO


c    Compute surface stress : (NB: z0 is a common in surfdat.h)
c     DO ig=1,ngrid
c        cd = (0.4/log(zzlay(ig,1)/z0(ig)))**2
c        zstress(ig) = rho(ig,1)*cd*(zu(ig,1)**2 + zv(ig,1)**2)
c     ENDDO

c     Sum of fluxes in solar spectral bands (for output only)
      DO ig=1,ngrid
             fluxtop_sw_tot(ig)=fluxtop_sw(ig,1) + fluxtop_sw(ig,2)
             fluxsurf_sw_tot(ig)=fluxsurf_sw(ig,1) + fluxsurf_sw(ig,2)
      ENDDO
c ******* TEST ******************************************************
      ztim1 = 999
      DO l=1,nlayer
        DO ig=1,ngrid
           if (pt(ig,l).lt.ztim1) then
               ztim1 = pt(ig,l)
               igmin = ig
               lmin = l 
           end if
        ENDDO
      ENDDO
      if(min(pt(igmin,lmin),zt(igmin,lmin)).lt.70.) then
        write(*,*) 'PHYSIQ: stability WARNING :'
        write(*,*) 'pt, zt Tmin = ', pt(igmin,lmin), zt(igmin,lmin),
     &              'ig l =', igmin, lmin
      end if
c *******************************************************************
      
c     ---------------------
c     Outputs to the screen 
c     ---------------------

      IF (lwrite) THEN
         PRINT*,'Global diagnostics for the physics'
         PRINT*,'Variables and their increments x and dx/dt * dt'
         WRITE(*,'(a6,a10,2a15)') 'Ts','dTs','ps','dps'
         WRITE(*,'(2f10.5,2f15.5)')
     s   tsurf(igout),zdtsurf(igout)*ptimestep,
     s   pplev(igout,1),pdpsrf(igout)*ptimestep
         WRITE(*,'(a4,a6,5a10)') 'l','u','du','v','dv','T','dT'
         WRITE(*,'(i4,6f10.5)') (l,
     s   pu(igout,l),pdu(igout,l)*ptimestep,
     s   pv(igout,l),pdv(igout,l)*ptimestep,
     s   pt(igout,l),pdt(igout,l)*ptimestep,
     s   l=1,nlayer)
      ENDIF ! of IF (lwrite)

c        ----------------------------------------------------------
c        ----------------------------------------------------------
c        INTERPOLATIONS IN THE SURFACE-LAYER
c        ----------------------------------------------------------
c        ----------------------------------------------------------

           n_out=0 ! number of elements in the z_out array.
                   ! for z_out=[3.,2.,1.,0.5,0.1], n_out must be set
                   ! to 5
           IF (n_out .ne. 0) THEN

           IF(.NOT. ALLOCATED(z_out)) ALLOCATE(z_out(n_out))
           IF(.NOT. ALLOCATED(T_out)) ALLOCATE(T_out(ngrid,n_out))
           IF(.NOT. ALLOCATED(u_out)) ALLOCATE(u_out(ngrid,n_out))

           z_out(:)=[3.,2.,1.,0.5,0.1]
           u_out(:,:)=0.
           T_out(:,:)=0.

           call pbl_parameters(ngrid,nlayer,ps,zplay,z0,
     & g,zzlay,zzlev,zu,zv,wstar,hfmax_th,zmax_th,tsurf,zh,z_out,n_out,
     & T_out,u_out,ustar,tstar,L_mo,vhf,vvv)

#ifndef MESOSCALE
            IF (ngrid .eq. 1) THEN
               dimout=0
            ELSE
               dimout=2
            ENDIF
            DO n=1,n_out
               write(zstring, '(F8.6)') z_out(n)
               call WRITEDIAGFI(ngrid,'T_out_'//trim(zstring),
     &   'potential temperature at z_out','K',dimout,T_out(:,n))
               call WRITEDIAGFI(ngrid,'u_out_'//trim(zstring),
     &   'horizontal velocity norm at z_out','m/s',dimout,u_out(:,n))
            ENDDO
            call WRITEDIAGFI(ngrid,'u_star',
     &   'friction velocity','m/s',dimout,ustar)
            call WRITEDIAGFI(ngrid,'teta_star',
     &   'friction potential temperature','K',dimout,tstar)
            call WRITEDIAGFI(ngrid,'L',
     &   'Monin Obukhov length','m',dimout,L_mo)
            call WRITEDIAGFI(ngrid,'vvv',
     &   'Vertical velocity variance at zout','m',dimout,vvv)
            call WRITEDIAGFI(ngrid,'vhf',
     &   'Vertical heat flux at zout','m',dimout,vhf)
#else
         T_out1(:)=T_out(:,1)
         u_out1(:)=u_out(:,1)
#endif

         ENDIF

c        ----------------------------------------------------------
c        ----------------------------------------------------------
c        END OF SURFACE LAYER INTERPOLATIONS
c        ----------------------------------------------------------
c        ----------------------------------------------------------

      IF (ngrid.NE.1) THEN

#ifndef MESOSCALE
c        -------------------------------------------------------------------
c        Writing NetCDF file  "RESTARTFI" at the end of the run
c        -------------------------------------------------------------------
c        Note: 'restartfi' is stored just before dynamics are stored
c              in 'restart'. Between now and the writting of 'restart',
c              there will have been the itau=itau+1 instruction and
c              a reset of 'time' (lastacll = .true. when itau+1= itaufin)
c              thus we store for time=time+dtvr

         IF(lastcall) THEN
            ztime_fin = ptime + ptimestep/(float(iphysiq)*daysec) 
            write(*,*)'PHYSIQ: for physdem ztime_fin =',ztime_fin
            call physdem1("restartfi.nc",long,lati,nsoilmx,nq,
     .              ptimestep,pday,
     .              ztime_fin,tsurf,tsoil,co2ice,emis,q2,qsurf,
     .              area,albedodat,inertiedat,zmea,zstd,zsig,
     .              zgam,zthe)
         ENDIF
#endif

c        -------------------------------------------------------------------
c        Calculation of diagnostic variables written in both stats and
c          diagfi files
c        -------------------------------------------------------------------

         if (tracer) then
         
           if(doubleq) then
              do ig=1,ngrid 
                dqdustsurf(ig) = 
     &                zdqssed(ig,igcm_dust_mass)*tauscaling(ig)
                dndustsurf(ig) = 
     &                zdqssed(ig,igcm_dust_number)*tauscaling(ig)
                ndust(ig,:) =
     &                pq(ig,:,igcm_dust_number)*tauscaling(ig)
                qdust(ig,:) =
     &                pq(ig,:,igcm_dust_mass)*tauscaling(ig)
              enddo
              if (scavenging) then
                do ig=1,ngrid 
                  dqdustsurf(ig) = dqdustsurf(ig) + 
     &                     zdqssed(ig,igcm_ccn_mass)*tauscaling(ig)
                  dndustsurf(ig) = dndustsurf(ig) + 
     &                     zdqssed(ig,igcm_ccn_number)*tauscaling(ig)
                  nccn(ig,:) =
     &                     pq(ig,:,igcm_ccn_number)*tauscaling(ig)
                  qccn(ig,:) =
     &                     pq(ig,:,igcm_ccn_mass)*tauscaling(ig)
                enddo
              endif
           endif
                    
           if (water) then
             mtot(:)=0
             icetot(:)=0
             rave(:)=0
             tauTES(:)=0
             do ig=1,ngrid 
               do l=1,nlayermx
                 mtot(ig) = mtot(ig) + 
     &                      zq(ig,l,igcm_h2o_vap) * 
     &                      (pplev(ig,l) - pplev(ig,l+1)) / g
                 icetot(ig) = icetot(ig) + 
     &                        zq(ig,l,igcm_h2o_ice) * 
     &                        (pplev(ig,l) - pplev(ig,l+1)) / g
c                Computing abs optical depth at 825 cm-1 in each
c                  layer to simulate NEW TES retrieval
                 Qabsice = min(
     &             max(0.4e6*rice(ig,l)*(1.+nuice_ref)-0.05 ,0.),1.2
     &                        )
                 opTES(ig,l)= 0.75 * Qabsice * 
     &             zq(ig,l,igcm_h2o_ice) *
     &             (pplev(ig,l) - pplev(ig,l+1)) / g
     &             / (rho_ice * rice(ig,l) * (1.+nuice_ref))
                 tauTES(ig)=tauTES(ig)+ opTES(ig,l) 
               enddo
c              rave(ig)=rave(ig)/max(icetot(ig),1.e-30)       ! mass weight
c               if (icetot(ig)*1e3.lt.0.01) rave(ig)=0.
             enddo
             call watersat(ngridmx*nlayermx,zt,pplay,zqsat)
             satu(:,:) = zq(:,:,igcm_h2o_vap)/zqsat(:,:)

             if (scavenging) then
               Nccntot(:)= 0
               Mccntot(:)= 0
               rave(:)=0
               do ig=1,ngrid 
                 do l=1,nlayermx
                    Nccntot(ig) = Nccntot(ig) + 
     &              zq(ig,l,igcm_ccn_number)*tauscaling(ig)
     &              *(pplev(ig,l) - pplev(ig,l+1)) / g
                    Mccntot(ig) = Mccntot(ig) + 
     &              zq(ig,l,igcm_ccn_mass)*tauscaling(ig)
     &              *(pplev(ig,l) - pplev(ig,l+1)) / g
cccc Column integrated effective ice radius 
cccc is weighted by total ice surface area (BETTER than total ice mass) 
                    rave(ig) = rave(ig) + 
     &                      tauscaling(ig) *
     &                      zq(ig,l,igcm_ccn_number) *
     &                      (pplev(ig,l) - pplev(ig,l+1)) / g * 
     &                      rice(ig,l) * rice(ig,l)*  (1.+nuice_ref)
                 enddo
               rave(ig)=(icetot(ig)/rho_ice+Mccntot(ig)/rho_dust)*0.75
     &                  /max(pi*rave(ig),1.e-30) ! surface weight
               if (icetot(ig)*1e3.lt.0.01) rave(ig)=0.
               enddo
             else ! of if (scavenging)
               rave(:)=0
               do ig=1,ngrid 
                 do l=1,nlayermx
                 rave(ig) = rave(ig) + 
     &                      zq(ig,l,igcm_h2o_ice) *
     &                      (pplev(ig,l) - pplev(ig,l+1)) / g * 
     &                      rice(ig,l) * (1.+nuice_ref)
                 enddo
                 rave(ig) = max(rave(ig) / 
     &             max(icetot(ig),1.e-30),1.e-30) ! mass weight
               enddo
             endif ! of if (scavenging)

           endif ! of if (water)
         endif ! of if (tracer)

c        -----------------------------------------------------------------
c        WSTATS: Saving statistics
c        -----------------------------------------------------------------
c        ("stats" stores and accumulates 8 key variables in file "stats.nc"
c        which can later be used to make the statistic files of the run:
c        "stats")          only possible in 3D runs !
         
       IF (callstats) THEN

        call wstats(ngrid,"ps","Surface pressure","Pa",2,ps)
        call wstats(ngrid,"tsurf","Surface temperature","K",2,tsurf)
        call wstats(ngrid,"co2ice","CO2 ice cover",
     &                "kg.m-2",2,co2ice)
        call wstats(ngrid,"tauref","reference dod at 610 Pa","NU",
     &                2,tauref)
        call wstats(ngrid,"fluxsurf_lw",
     &                "Thermal IR radiative flux to surface","W.m-2",2,
     &                fluxsurf_lw)
        call wstats(ngrid,"fluxsurf_sw",
     &                "Solar radiative flux to surface","W.m-2",2,
     &                fluxsurf_sw_tot)
        call wstats(ngrid,"fluxtop_lw",
     &                "Thermal IR radiative flux to space","W.m-2",2,
     &                fluxtop_lw)
        call wstats(ngrid,"fluxtop_sw",
     &                "Solar radiative flux to space","W.m-2",2,
     &                fluxtop_sw_tot)
        call wstats(ngrid,"temp","Atmospheric temperature","K",3,zt)
        call wstats(ngrid,"u","Zonal (East-West) wind","m.s-1",3,zu)
        call wstats(ngrid,"v","Meridional (North-South) wind",
     &                "m.s-1",3,zv)
        call wstats(ngrid,"w","Vertical (down-up) wind",
     &                "m.s-1",3,pw)
        call wstats(ngrid,"rho","Atmospheric density","kg/m3",3,rho)
        call wstats(ngrid,"pressure","Pressure","Pa",3,pplay)
          call wstats(ngrid,"q2",
     &                "Boundary layer eddy kinetic energy",
     &                "m2.s-2",3,q2)
          call wstats(ngrid,"emis","Surface emissivity","w.m-1",2,
     &                emis)
c          call wstats(ngrid,"ssurf","Surface stress","N.m-2",
c    &                2,zstress)
c          call wstats(ngrid,"sw_htrt","sw heat.rate",
c    &                 "W.m-2",3,zdtsw)
c          call wstats(ngrid,"lw_htrt","lw heat.rate",
c    &                 "W.m-2",3,zdtlw)

          if (calltherm) then
            call wstats(ngrid,"zmax_th","Height of thermals",
     &                "m",2,zmax_th)
            call wstats(ngrid,"hfmax_th","Max thermals heat flux",
     &                "K.m/s",2,hfmax_th)
            call wstats(ngrid,"wstar",
     &                "Max vertical velocity in thermals",
     &                "m/s",2,wstar)
          endif

           if (tracer) then
             if (water) then
               vmr=zq(1:ngridmx,1:nlayermx,igcm_h2o_vap)
     &      *mmean(1:ngridmx,1:nlayermx)/mmol(igcm_h2o_vap)
               call wstats(ngrid,"vmr_h2ovap",
     &                    "H2O vapor volume mixing ratio","mol/mol",
     &                    3,vmr)
               vmr=zq(1:ngridmx,1:nlayermx,igcm_h2o_ice)
     &      *mmean(1:ngridmx,1:nlayermx)/mmol(igcm_h2o_ice)
               call wstats(ngrid,"vmr_h2oice",
     &                    "H2O ice volume mixing ratio","mol/mol",
     &                    3,vmr)
               vmr=zqsat(1:ngridmx,1:nlayermx)
     &      *mmean(1:ngridmx,1:nlayermx)/mmol(igcm_h2o_vap)
               call wstats(ngrid,"vmr_h2osat",
     &                    "saturation volume mixing ratio","mol/mol",
     &                    3,vmr)
               call wstats(ngrid,"h2o_ice_s",
     &                    "surface h2o_ice","kg/m2",
     &                    2,qsurf(1,igcm_h2o_ice))
               call wstats(ngrid,'albedo',
     &                         'albedo',
     &                         '',2,albedo(1:ngridmx,1))
               call wstats(ngrid,"mtot",
     &                    "total mass of water vapor","kg/m2",
     &                    2,mtot)
               call wstats(ngrid,"icetot",
     &                    "total mass of water ice","kg/m2",
     &                    2,icetot)
               call wstats(ngrid,"reffice",
     &                    "Mean reff","m",
     &                    2,rave)
              call wstats(ngrid,"Nccntot",
     &                    "condensation nuclei","Nbr/m2",
     &                    2,Nccntot)
              call wstats(ngrid,"Mccntot",
     &                    "condensation nuclei mass","kg/m2",
     &                    2,Mccntot)
              call wstats(ngrid,"rice",
     &                    "Ice particle size","m",
     &                    3,rice)
               if (.not.activice) then
                 call wstats(ngrid,"tauTESap",
     &                      "tau abs 825 cm-1","",
     &                      2,tauTES)
               else
                 call wstats(ngridmx,'tauTES',
     &                   'tau abs 825 cm-1',
     &                  '',2,taucloudtes)
               endif

             endif ! of if (water)
             
             
           if (dustbin.ne.0) then
          
             call wstats(ngridmx,'tau','taudust','SI',2,tau(1,1))
             
             if (doubleq) then
c            call wstats(ngridmx,'qsurf','qsurf',
c     &                       'kg.m-2',2,qsurf(1,igcm_dust_mass))
c            call wstats(ngridmx,'Nsurf','N particles',
c     &                       'N.m-2',2,qsurf(1,igcm_dust_number))
c            call wstats(ngridmx,'dqsdev','ddevil lift',
c    &                       'kg.m-2.s-1',2,zdqsdev(1,1))
c            call wstats(ngridmx,'dqssed','sedimentation',
c     &                       'kg.m-2.s-1',2,zdqssed(1,1))
c            call wstats(ngridmx,'dqsdif','diffusion',
c     &                       'kg.m-2.s-1',2,zdqsdif(1,1))
               call wstats(ngridmx,'dqsdust',
     &                        'deposited surface dust mass',
     &                        'kg.m-2.s-1',2,dqdustsurf)
               call wstats(ngridmx,'dqndust',
     &                        'deposited surface dust number',
     &                        'number.m-2.s-1',2,dndustsurf)
               call wstats(ngridmx,'reffdust','reffdust',
     &                        'm',3,rdust*ref_r0)
               call wstats(ngridmx,'dustq','Dust mass mr',
     &                        'kg/kg',3,qdust)
               call wstats(ngridmx,'dustN','Dust number',
     &                        'part/kg',3,ndust)
             else
               do iq=1,dustbin
                 write(str2(1:2),'(i2.2)') iq
                 call wstats(ngridmx,'q'//str2,'mix. ratio',
     &                         'kg/kg',3,zq(1,1,iq))
                 call wstats(ngridmx,'qsurf'//str2,'qsurf',
     &                         'kg.m-2',2,qsurf(1,iq))
               end do
             endif ! (doubleq)

             if (scavenging) then
               call wstats(ngridmx,'ccnq','CCN mass mr',
     &                        'kg/kg',3,qccn)
               call wstats(ngridmx,'ccnN','CCN number',
     &                        'part/kg',3,nccn)
             endif ! (scavenging)
          
          endif ! (dustbin.ne.0)


             if (thermochem.or.photochem) then
                do iq=1,nq
                  if (noms(iq) .ne. "dust_mass" .and.
     $                 noms(iq) .ne. "dust_number" .and.
     $                 noms(iq) .ne. "ccn_mass" .and.
     $                 noms(iq) .ne. "ccn_number" .and.
     $                 noms(iq) .ne. "h2o_vap" .and.
     $                 noms(iq) .ne. "h2o_ice") then
                    vmr(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,iq)
     &                          *mmean(1:ngrid,1:nlayer)/mmol(iq)
                    rhopart(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,iq)
     &                          *rho(1:ngrid,1:nlayer)*n_avog/
     &                           (1000*mmol(iq))
                   call wstats(ngrid,"vmr_"//trim(noms(iq)),
     $                         "Volume mixing ratio","mol/mol",3,vmr)
!                   call wstats(ngrid,"rho_"//trim(noms(iq)),
!     $                     "Number density","cm-3",3,rhopart)
!                   call writediagfi(ngrid,"rho_"//trim(noms(iq)),
!     $                     "Number density","cm-3",3,rhopart)
                   if ((noms(iq).eq."o") .or. (noms(iq).eq."co2").or.
     $                 (noms(iq).eq."o3")) then                     
                      call writediagfi(ngrid,"vmr_"//trim(noms(iq)), 
     $                         "Volume mixing ratio","mol/mol",3,vmr)
                   end if
                   do ig = 1,ngrid
                      colden(ig,iq) = 0.                           
                   end do
                   do l=1,nlayer                                    
                      do ig=1,ngrid                                  
                         colden(ig,iq) = colden(ig,iq) + zq(ig,l,iq)
     $                                  *(pplev(ig,l)-pplev(ig,l+1)) 
     $                                  *6.022e22/(mmol(iq)*g)       
                      end do                                        
                   end do                                          
                   call wstats(ngrid,"c_"//trim(noms(iq)),           
     $                         "column","mol cm-2",2,colden(1,iq))  
                   call writediagfi(ngrid,"c_"//trim(noms(iq)),  
     $                             "column","mol cm-2",2,colden(1,iq))
                  end if ! of if (noms(iq) .ne. "dust_mass" ...)
                end do ! of do iq=1,nq
             end if ! of if (thermochem.or.photochem)

           end if ! of if (tracer)

           IF(lastcall) THEN
             write (*,*) "Writing stats..."
             call mkstats(ierr)
           ENDIF

         ENDIF !if callstats

c        (Store EOF for Mars Climate database software)
         IF (calleofdump) THEN
            CALL eofdump(ngrid, nlayer, zu, zv, zt, rho, ps)
         ENDIF


#ifdef MESOSCALE
      !!!
      !!! OUTPUT FIELDS
      !!!
      wtsurf(1:ngrid) = tsurf(1:ngrid)    !! surface temperature
      wco2ice(1:ngrid) = co2ice(1:ngrid)  !! co2 ice
      TAU_lay(:)=tau(:,1)!!true opacity (not a reference like tauref)
      IF (tracer) THEN
      mtot(1:ngrid) = mtot(1:ngrid) * 1.e6 / rho_ice
      icetot(1:ngrid) = icetot(1:ngrid) * 1.e6 / rho_ice
      !! JF 
      IF (igcm_dust_mass .ne. 0) THEN
        qsurfdust(1:ngrid) = qsurf(1:ngrid,igcm_dust_mass)
      ENDIF
      IF (igcm_h2o_ice .ne. 0) THEN      
        qsurfice(1:ngrid) = qsurf(1:ngrid,igcm_h2o_ice)
        vmr=1.e6 * zq(1:ngrid,1:nlayer,igcm_h2o_ice)
     .           *mmean(1:ngridmx,1:nlayermx) / mmol(igcm_h2o_ice)
      ENDIF
      !! Dust quantity integration along the vertical axe
      dustot(:)=0
      IF (igcm_dust_mass .ne. 0) THEN
      do ig=1,ngrid
       do l=1,nlayermx
        dustot(ig) = dustot(ig) +
     &               zq(ig,l,igcm_dust_mass)
     &               *  (pplev(ig,l) - pplev(ig,l+1)) / g
       enddo
      enddo
      ENDIF
      ENDIF
      !! TAU water ice as seen by TES
      if (activice) tauTES = taucloudtes
c AUTOMATICALLY GENERATED FROM REGISTRY
#include "fill_save.inc"
#else
#ifndef MESOINI

c        ==========================================================
c        WRITEDIAGFI: Outputs in netcdf file "DIAGFI", containing
c          any variable for diagnostic (output with period
c          "ecritphy", set in "run.def")
c        ==========================================================
c        WRITEDIAGFI can ALSO be called from any other subroutines
c        for any variables !!
c        call WRITEDIAGFI(ngrid,"emis","Surface emissivity","w.m-1",2,
c    &                  emis)
c        call WRITEDIAGFI(ngrid,"pplay","Pressure","Pa",3,zplay)
c        call WRITEDIAGFI(ngrid,"pplev","Pressure","Pa",3,zplev)
         call WRITEDIAGFI(ngrid,"tsurf","Surface temperature","K",2,
     &                  tsurf)
         call WRITEDIAGFI(ngrid,"ps","surface pressure","Pa",2,ps)
         call WRITEDIAGFI(ngrid,"co2ice","co2 ice thickness"
     &                                         ,"kg.m-2",2,co2ice)

         call WRITEDIAGFI(ngrid,"temp7","temperature in layer 7",
     &                  "K",2,zt(1,7))
         call WRITEDIAGFI(ngrid,"fluxsurf_lw","fluxsurf_lw","W.m-2",2,
     &                  fluxsurf_lw)
         call WRITEDIAGFI(ngrid,"fluxsurf_sw","fluxsurf_sw","W.m-2",2,
     &                  fluxsurf_sw_tot)
         call WRITEDIAGFI(ngrid,"fluxtop_lw","fluxtop_lw","W.m-2",2,
     &                  fluxtop_lw)
         call WRITEDIAGFI(ngrid,"fluxtop_sw","fluxtop_sw","W.m-2",2,
     &                  fluxtop_sw_tot)
         call WRITEDIAGFI(ngrid,"temp","temperature","K",3,zt)
         call WRITEDIAGFI(ngrid,"u","Zonal wind","m.s-1",3,zu)
         call WRITEDIAGFI(ngrid,"v","Meridional wind","m.s-1",3,zv)
         call WRITEDIAGFI(ngrid,"w","Vertical wind","m.s-1",3,pw)
         call WRITEDIAGFI(ngrid,"rho","density","none",3,rho)
c        call WRITEDIAGFI(ngrid,"q2","q2","kg.m-3",3,q2)
c        call WRITEDIAGFI(ngrid,'Teta','T potentielle','K',3,zh)
         call WRITEDIAGFI(ngrid,"pressure","Pressure","Pa",3,zplay)
c        call WRITEDIAGFI(ngrid,"ssurf","Surface stress","N.m-2",2,
c    &                  zstress)
c        call WRITEDIAGFI(ngridmx,'sw_htrt','sw heat. rate',
c    &                   'w.m-2',3,zdtsw)
c        call WRITEDIAGFI(ngridmx,'lw_htrt','lw heat. rate',
c    &                   'w.m-2',3,zdtlw)
            if (.not.activice) then
               CALL WRITEDIAGFI(ngridmx,'tauTESap',
     &                         'tau abs 825 cm-1',
     &                         '',2,tauTES)
             else
               CALL WRITEDIAGFI(ngridmx,'tauTES',
     &                         'tau abs 825 cm-1',
     &                         '',2,taucloudtes)
             endif

#else
     !!! this is to ensure correct initialisation of mesoscale model
        call WRITEDIAGFI(ngrid,"tsurf","Surface temperature","K",2,
     &                  tsurf)
        call WRITEDIAGFI(ngrid,"ps","surface pressure","Pa",2,ps) 
        call WRITEDIAGFI(ngrid,"co2ice","co2 ice thickness","kg.m-2",2,
     &                  co2ice)
        call WRITEDIAGFI(ngrid,"temp","temperature","K",3,zt)
        call WRITEDIAGFI(ngrid,"u","Zonal wind","m.s-1",3,zu)
        call WRITEDIAGFI(ngrid,"v","Meridional wind","m.s-1",3,zv)
        call WRITEDIAGFI(ngrid,"emis","Surface emissivity","w.m-1",2,
     &                  emis)
        call WRITEDIAGFI(ngrid,"tsoil","Soil temperature",
     &                       "K",3,tsoil)
        call WRITEDIAGFI(ngrid,"inertiedat","Soil inertia",
     &                       "K",3,inertiedat)
#endif


c        ----------------------------------------------------------
c        Outputs of the CO2 cycle
c        ----------------------------------------------------------

         if (tracer.and.(igcm_co2.ne.0)) then
!          call WRITEDIAGFI(ngrid,"co2_l1","co2 mix. ratio in 1st layer",
!    &                     "kg/kg",2,zq(1,1,igcm_co2))
          call WRITEDIAGFI(ngrid,"co2","co2 mass mixing ratio",
     &                     "kg/kg",3,zq(1,1,igcm_co2))
        
         ! Compute co2 column
         co2col(:)=0
         do l=1,nlayermx
           do ig=1,ngrid
             co2col(ig)=co2col(ig)+
     &                  zq(ig,l,igcm_co2)*(pplev(ig,l)-pplev(ig,l+1))/g
           enddo
         enddo
         call WRITEDIAGFI(ngrid,"co2col","CO2 column","kg.m-2",2,
     &                  co2col)
         endif ! of if (tracer.and.(igcm_co2.ne.0))

c        ----------------------------------------------------------
c        Outputs of the water cycle
c        ----------------------------------------------------------
         if (tracer) then
           if (water) then

#ifdef MESOINI
            !!!! waterice = q01, voir readmeteo.F90
            call WRITEDIAGFI(ngridmx,'q01',noms(igcm_h2o_ice),
     &                      'kg/kg',3,
     &                       zq(1:ngridmx,1:nlayermx,igcm_h2o_ice))
            !!!! watervapor = q02, voir readmeteo.F90
            call WRITEDIAGFI(ngridmx,'q02',noms(igcm_h2o_vap),
     &                      'kg/kg',3,
     &                       zq(1:ngridmx,1:nlayermx,igcm_h2o_vap))
            !!!! surface waterice qsurf02 (voir readmeteo)
            call WRITEDIAGFI(ngridmx,'qsurf02','surface tracer',
     &                      'kg.m-2',2,
     &                       qsurf(1:ngridmx,igcm_h2o_ice))
#endif

            CALL WRITEDIAGFI(ngridmx,'mtot',
     &                       'total mass of water vapor',
     &                       'kg/m2',2,mtot)
            CALL WRITEDIAGFI(ngridmx,'icetot',
     &                       'total mass of water ice',
     &                       'kg/m2',2,icetot)
            vmr=zq(1:ngridmx,1:nlayermx,igcm_h2o_ice)
     &      *mmean(1:ngridmx,1:nlayermx)/mmol(igcm_h2o_ice)
            call WRITEDIAGFI(ngridmx,'vmr_h2oice','h2o ice vmr',
     &                       'mol/mol',3,vmr)
            vmr=zq(1:ngridmx,1:nlayermx,igcm_h2o_vap)
     &      *mmean(1:ngridmx,1:nlayermx)/mmol(igcm_h2o_vap)
            call WRITEDIAGFI(ngridmx,'vmr_h2ovap','h2o vap vmr',
     &                       'mol/mol',3,vmr)
            CALL WRITEDIAGFI(ngridmx,'reffice',
     &                       'Mean reff',
     &                       'm',2,rave)
            CALL WRITEDIAGFI(ngrid,"Nccntot",
     &                    "condensation nuclei","Nbr/m2",
     &                    2,Nccntot)
            CALL WRITEDIAGFI(ngrid,"Mccntot",
     &                    "mass condensation nuclei","kg/m2",
     &                    2,Mccntot)
            call WRITEDIAGFI(ngridmx,'rice','Ice particle size',
     &                       'm',3,rice)
            call WRITEDIAGFI(ngridmx,'h2o_ice_s',
     &                       'surface h2o_ice',
     &                       'kg.m-2',2,qsurf(1,igcm_h2o_ice))
            CALL WRITEDIAGFI(ngridmx,'albedo',
     &                         'albedo',
     &                         '',2,albedo(1:ngridmx,1))
              if (tifeedback) then
                 call WRITEDIAGSOIL(ngridmx,"soiltemp",
     &                              "Soil temperature","K",
     &                              3,tsoil)
                 call WRITEDIAGSOIL(ngridmx,'soilti',
     &                       'Soil Thermal Inertia',
     &                       'J.s-1/2.m-2.K-1',3,inertiesoil)
              endif
           endif !(water)


           if (water.and..not.photochem) then
             iq=nq
c            write(str2(1:2),'(i2.2)') iq
c            call WRITEDIAGFI(ngridmx,'dqs'//str2,'dqscloud',
c    &                       'kg.m-2',2,zdqscloud(1,iq))
c            call WRITEDIAGFI(ngridmx,'dqch'//str2,'var chim',
c    &                       'kg/kg',3,zdqchim(1,1,iq))
c            call WRITEDIAGFI(ngridmx,'dqd'//str2,'var dif',
c    &                       'kg/kg',3,zdqdif(1,1,iq))
c            call WRITEDIAGFI(ngridmx,'dqa'//str2,'var adj',
c    &                       'kg/kg',3,zdqadj(1,1,iq))
c            call WRITEDIAGFI(ngridmx,'dqc'//str2,'var c',
c    &                       'kg/kg',3,zdqc(1,1,iq))
           endif  !(water.and..not.photochem)
         endif

c        ----------------------------------------------------------
c        Outputs of the dust cycle
c        ----------------------------------------------------------

        call WRITEDIAGFI(ngridmx,'tauref',
     &                    'Dust ref opt depth','NU',2,tauref)

         if (tracer.and.(dustbin.ne.0)) then
          call WRITEDIAGFI(ngridmx,'tau','taudust','SI',2,tau(1,1))
           if (doubleq) then
c            call WRITEDIAGFI(ngridmx,'qsurf','qsurf',
c     &                       'kg.m-2',2,qsurf(1,igcm_dust_mass))
c            call WRITEDIAGFI(ngridmx,'Nsurf','N particles',
c     &                       'N.m-2',2,qsurf(1,igcm_dust_number))
c            call WRITEDIAGFI(ngridmx,'dqsdev','ddevil lift',
c    &                       'kg.m-2.s-1',2,zdqsdev(1,1))
c            call WRITEDIAGFI(ngridmx,'dqssed','sedimentation',
c     &                       'kg.m-2.s-1',2,zdqssed(1,1))
c            call WRITEDIAGFI(ngridmx,'dqsdif','diffusion',
c     &                       'kg.m-2.s-1',2,zdqsdif(1,1))
             call WRITEDIAGFI(ngridmx,'dqsdust',
     &                        'deposited surface dust mass',
     &                        'kg.m-2.s-1',2,dqdustsurf)
             call WRITEDIAGFI(ngridmx,'dqndust',
     &                        'deposited surface dust number',
     &                        'number.m-2.s-1',2,dndustsurf)
             call WRITEDIAGFI(ngridmx,'reffdust','reffdust',
     &                        'm',3,rdust*ref_r0)
             call WRITEDIAGFI(ngridmx,'dustq','Dust mass mr',
     &                        'kg/kg',3,qdust)
             call WRITEDIAGFI(ngridmx,'dustN','Dust number',
     &                        'part/kg',3,ndust)
#ifdef MESOINI
             call WRITEDIAGFI(ngridmx,'dustq','Dust mass mr',
     &                        'kg/kg',3,qdust)
             call WRITEDIAGFI(ngridmx,'dustN','Dust number',
     &                        'part/kg',3,ndust)
             call WRITEDIAGFI(ngridmx,'ccn','Nuclei mass mr',
     &                        'kg/kg',3,qccn)
             call WRITEDIAGFI(ngridmx,'ccnN','Nuclei number',
     &                        'part/kg',3,nccn)
#endif
           else
             do iq=1,dustbin
               write(str2(1:2),'(i2.2)') iq
               call WRITEDIAGFI(ngridmx,'q'//str2,'mix. ratio',
     &                         'kg/kg',3,zq(1,1,iq))
               call WRITEDIAGFI(ngridmx,'qsurf'//str2,'qsurf',
     &                         'kg.m-2',2,qsurf(1,iq))
             end do
           endif ! (doubleq)

           if (scavenging) then
             call WRITEDIAGFI(ngridmx,'ccnq','CCN mass mr',
     &                        'kg/kg',3,qccn)
             call WRITEDIAGFI(ngridmx,'ccnN','CCN number',
     &                        'part/kg',3,nccn)
           endif ! (scavenging)

c          if (submicron) then
c            call WRITEDIAGFI(ngridmx,'dustsubm','subm mass mr',
c    &                        'kg/kg',3,pq(1,1,igcm_dust_submicron))
c          endif ! (submicron)
         end if  ! (tracer.and.(dustbin.ne.0))


c        ----------------------------------------------------------
c        Thermospheric outputs
c        ----------------------------------------------------------

         if(callthermos) then 

            call WRITEDIAGFI(ngridmx,"q15um","15 um cooling","K/s",
     $           3,zdtnlte)
            call WRITEDIAGFI(ngridmx,"quv","UV heating","K/s",
     $           3,zdteuv)
            call WRITEDIAGFI(ngridmx,"cond","Thermal conduction","K/s",
     $           3,zdtconduc)
            call WRITEDIAGFI(ngridmx,"qnir","NIR heating","K/s",
     $           3,zdtnirco2)

         endif  !(callthermos)

c        ----------------------------------------------------------
c        ----------------------------------------------------------
c        PBL OUTPUS
c        ----------------------------------------------------------
c        ----------------------------------------------------------

c        ----------------------------------------------------------
c        Outputs of thermals
c        ----------------------------------------------------------
         if (calltherm) then

!        call WRITEDIAGFI(ngrid,'dtke',
!     &              'tendance tke thermiques','m**2/s**2',
!     &                         3,dtke_th)
!        call WRITEDIAGFI(ngrid,'d_u_ajs',
!     &              'tendance u thermiques','m/s',
!     &                         3,pdu_th*ptimestep)
!        call WRITEDIAGFI(ngrid,'d_v_ajs',
!     &              'tendance v thermiques','m/s',
!     &                         3,pdv_th*ptimestep)
!        if (tracer) then
!        if (nq .eq. 2) then
!        call WRITEDIAGFI(ngrid,'deltaq_th',
!     &              'delta q thermiques','kg/kg',
!     &                         3,ptimestep*pdq_th(:,:,2))
!        endif
!        endif

        call WRITEDIAGFI(ngridmx,'zmax_th',
     &              'hauteur du thermique','m',
     &                         2,zmax_th)
        call WRITEDIAGFI(ngridmx,'hfmax_th',
     &              'maximum TH heat flux','K.m/s',
     &                         2,hfmax_th)
        call WRITEDIAGFI(ngridmx,'wstar',
     &              'maximum TH vertical velocity','m/s',
     &                         2,wstar)

         endif

c        ----------------------------------------------------------
c        ----------------------------------------------------------
c        END OF PBL OUTPUS
c        ----------------------------------------------------------
c        ----------------------------------------------------------


c        ----------------------------------------------------------
c        Output in netcdf file "diagsoil.nc" for subterranean
c          variables (output every "ecritphy", as for writediagfi)
c        ----------------------------------------------------------

         ! Write soil temperature
!        call writediagsoil(ngrid,"soiltemp","Soil temperature","K",
!    &                     3,tsoil)
         ! Write surface temperature
!        call writediagsoil(ngrid,"tsurf","Surface temperature","K",
!    &                     2,tsurf)

c        ==========================================================
c        END OF WRITEDIAGFI
c        ==========================================================
#endif

      ELSE     ! if(ngrid.eq.1)

         write(*,'("Ls =",f11.6," tauref(",f4.0," Pa) =",f9.6)')
     &    zls*180./pi,odpref,tauref
        CALL WRITEDIAGFI(ngridmx, "zls", 
     &                   "Solar Zenith Angle", 'degs', 
     &                   0, zls*180./pi)
c      ----------------------------------------------------------------------
c      Output in grads file "g1d" (ONLY when using testphys1d)
c      (output at every X physical timestep)
c      ----------------------------------------------------------------------
c
c        CALL writeg1d(ngrid,1,fluxsurf_lw,'Fs_ir','W.m-2')
c         CALL writeg1d(ngrid,1,tsurf,'tsurf','K')
c         CALL writeg1d(ngrid,1,ps,'ps','Pa')
         
c         CALL writeg1d(ngrid,nlayer,zt,'T','K')
c        CALL writeg1d(ngrid,nlayer,pu,'u','m.s-1')
c        CALL writeg1d(ngrid,nlayer,pv,'v','m.s-1')
c        CALL writeg1d(ngrid,nlayer,pw,'w','m.s-1')

! THERMALS STUFF 1D
         if(calltherm) then

        call WRITEDIAGFI(ngridmx,'lmax_th',
     &              'hauteur du thermique','point',
     &                         0,lmax_th_out)
        call WRITEDIAGFI(ngridmx,'zmax_th',
     &              'hauteur du thermique','m',
     &                         0,zmax_th)
        call WRITEDIAGFI(ngridmx,'hfmax_th',
     &              'maximum TH heat flux','K.m/s',
     &                         0,hfmax_th)
        call WRITEDIAGFI(ngridmx,'wstar',
     &              'maximum TH vertical velocity','m/s',
     &                         0,wstar)

         co2col(:)=0.
         if (tracer) then
         do l=1,nlayermx
           do ig=1,ngrid
             co2col(ig)=co2col(ig)+
     &                  zq(ig,l,1)*(pplev(ig,l)-pplev(ig,l+1))/g
         enddo
         enddo

         end if
         call WRITEDIAGFI(ngrid,'co2col','integrated co2 mass'          &
     &                                      ,'kg/m-2',0,co2col)
         endif ! of if (calltherm)

C          call WRITEDIAGFI(ngrid,'w','vertical velocity'                 &
C      &                              ,'m/s',1,pw)
C          call WRITEDIAGFI(ngrid,"q2","q2","kg.m-3",1,q2)
C          call WRITEDIAGFI(ngrid,"tsurf","Surface temperature","K",0,
C      &                  tsurf)
C          call WRITEDIAGFI(ngrid,"u","u wind","m/s",1,zu)
C          call WRITEDIAGFI(ngrid,"v","v wind","m/s",1,zv)

         call WRITEDIAGFI(ngrid,"pplay","Pressure","Pa",1,zplay)
         call WRITEDIAGFI(ngrid,"pplev","Pressure","Pa",1,zplev)
C          call WRITEDIAGFI(ngrid,"rho","rho","kg.m-3",1,rho)
C !         call WRITEDIAGFI(ngrid,"dtrad","rad. heat. rate",              &
C !     &              "K.s-1",1,dtrad/zpopsk)
C !        call WRITEDIAGFI(ngridmx,'sw_htrt','sw heat. rate',
C !     &                   'w.m-2',1,zdtsw/zpopsk)
C !        call WRITEDIAGFI(ngridmx,'lw_htrt','lw heat. rate',
C !     &                   'w.m-2',1,zdtlw/zpopsk)
C          call WRITEDIAGFI(ngrid,"co2ice","co2 ice thickness"
C      &                                   ,"kg.m-2",0,co2ice)

! or output in diagfi.nc (for testphys1d)
         call WRITEDIAGFI(ngridmx,'ps','Surface pressure','Pa',0,ps)
         call WRITEDIAGFI(ngridmx,'temp','Temperature',
     &                       'K',1,zt)
     
         colden(:,:) = 0.
         if(tracer) then
         
c           CALL writeg1d(ngrid,1,tau,'tau','SI')
            do iq=1,nq
c              CALL writeg1d(ngrid,nlayer,zq(1,1,iq),noms(iq),'kg/kg') 

          do l=1,nlayer                                    
               do ig=1,ngrid                                  
                    colden(ig,iq) = colden(ig,iq) + zq(ig,l,iq)
     $                                  *(pplev(ig,l)-pplev(ig,l+1)) 
     $                                 *6.022e22/(mmol(iq)*g)       
               end do                                        
          end do       
               


               call WRITEDIAGFI(ngridmx,trim(noms(iq)) // "_col",
     &              trim(noms(iq))//" column", 'cm-2',0,colden(1,iq))
     
               call WRITEDIAGFI(ngridmx,trim(noms(iq)),
     &              trim(noms(iq)),'kg/kg',1,zq(1,1,iq))

c            ===========================

            end do
           if (doubleq) then
             call WRITEDIAGFI(ngridmx,'rdust','rdust',
     &                        'm',1,rdust)
           endif
           if (water.AND.tifeedback) then
             call WRITEDIAGFI(ngridmx,"soiltemp",
     &                              "Soil temperature","K",
     &                              1,tsoil)
             call WRITEDIAGFI(ngridmx,'soilti',
     &                       'Soil Thermal Inertia',
     &                       'J.s-1/2.m-2.K-1',1,inertiesoil)
           endif
         end if
         
cccccccccccccccccc scavenging & water outputs 1D TN ccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF (water) THEN

           if (.not.activice) then

             tauTES=0
             do l=1,nlayermx
               Qabsice = min(
     &             max(0.4e6*rice(1,l)*(1.+nuice_ref)-0.05 ,0.),1.2
     &                      )
               opTES(1,l)= 0.75 * Qabsice *
     &             zq(1,l,igcm_h2o_ice) *
     &             (pplev(1,l) - pplev(1,l+1)) / g
     &             / (rho_ice * rice(1,l) * (1.+nuice_ref))
               tauTES=tauTES+ opTES(1,l)
             enddo
             CALL WRITEDIAGFI(ngridmx,'tauTESap',
     &                         'tau abs 825 cm-1',
     &                         '',0,tauTES)
           else

             CALL WRITEDIAGFI(ngridmx,'tauTES',
     &                         'tau abs 825 cm-1',
     &                         '',0,taucloudtes)
           endif
     
           mtot = 0
           icetot = 0
           h2otot = qsurf(1,igcm_h2o_ice)

           do l=1,nlayer
             mtot = mtot +  zq(1,l,igcm_h2o_vap) 
     &                 * (pplev(1,l) - pplev(1,l+1)) / g
             icetot = icetot +  zq(1,l,igcm_h2o_ice) 
     &                 * (pplev(1,l) - pplev(1,l+1)) / g
           end do
           h2otot = h2otot+mtot+icetot

             CALL WRITEDIAGFI(ngridmx,'h2otot',
     &                         'h2otot',
     &                         'kg/m2',0,h2otot)
             CALL WRITEDIAGFI(ngridmx,'mtot',
     &                         'mtot',
     &                         'kg/m2',0,mtot)
             CALL WRITEDIAGFI(ngridmx,'icetot',
     &                         'icetot',
     &                         'kg/m2',0,icetot)

           if (scavenging) then

             rave = 0
             do l=1,nlayer
cccc Column integrated effective ice radius
cccc is weighted by total ice surface area (BETTER) 
             rave = rave + tauscaling(1) *
     &              zq(1,l,igcm_ccn_number) *
     &              (pplev(1,l) - pplev(1,l+1)) / g * 
     &              rice(1,l) * rice(1,l)*  (1.+nuice_ref)
             enddo
             rave=icetot*0.75/max(rave*pi*rho_ice,1.e-30) ! surface weight

              Nccntot= 0
              call watersat(ngridmx*nlayermx,zt,pplay,zqsat)
               do l=1,nlayermx
                   Nccntot = Nccntot + 
     &              zq(1,l,igcm_ccn_number)*tauscaling(1)
     &              *(pplev(1,l) - pplev(1,l+1)) / g
                   satu(1,l) = zq(1,l,igcm_h2o_vap)/zqsat(1,l)
                   satu(1,l) = (max(satu(1,l),float(1))-1)
!     &                      * zq(1,l,igcm_h2o_vap) *
!     &                      (pplev(1,l) - pplev(1,l+1)) / g
               enddo
               call WRITEDIAGFI(ngrid,"satu","vap in satu","kg/kg",1,
     &                    satu)
               CALL WRITEDIAGFI(ngridmx,'Nccntot',
     &                         'Nccntot',
     &                         'nbr/m2',0,Nccntot)

             call WRITEDIAGFI(ngridmx,'zdqsed_dustq'
     & ,'sedimentation q','kg.m-2.s-1',1,zdqsed(1,:,igcm_dust_mass))
             call WRITEDIAGFI(ngridmx,'zdqsed_dustN'
     &,'sedimentation N','Nbr.m-2.s-1',1,
     &                                 zdqsed(1,:,igcm_dust_number))

           else ! of if (scavenging)

cccc Column integrated effective ice radius 
cccc is weighted by total ice mass         (LESS GOOD)
             rave = 0
             do l=1,nlayer
               rave = rave + zq(1,l,igcm_h2o_ice)
     &              * (pplev(1,l) - pplev(1,l+1)) / g
     &              *  rice(1,l) * (1.+nuice_ref)
             enddo
             rave=max(rave/max(icetot,1.e-30),1.e-30) ! mass weight
           endif ! of if (scavenging)

C            CALL WRITEDIAGFI(ngridmx,'reffice',
C      &                      'reffice',
C      &                      'm',0,rave)

C            do iq=1,nq
C              call WRITEDIAGFI(ngridmx,trim(noms(iq))//'_s',
C      &            trim(noms(iq))//'_s','kg/kg',0,qsurf(1,iq))
C            end do
     
C         call WRITEDIAGFI(ngridmx,'zdqcloud_ice','cloud ice',
C      &            'kg.m-2.s-1',1,zdqcloud(1,:,igcm_h2o_ice))
C         call WRITEDIAGFI(ngridmx,'zdqcloud_vap','cloud vap',
C      &            'kg.m-2.s-1',1,zdqcloud(1,:,igcm_h2o_vap))
C         call WRITEDIAGFI(ngridmx,'zdqcloud','cloud ice',
C      &            'kg.m-2.s-1',1,zdqcloud(1,:,igcm_h2o_ice)
C      &                          +zdqcloud(1,:,igcm_h2o_vap))
     
C           call WRITEDIAGFI(ngrid,"rice","ice radius","m",1,
C      &                    rice)
         ENDIF ! of IF (water)
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


         zlocal(1)=-log(pplay(1,1)/pplev(1,1))* Rnew(1,1)*zt(1,1)/g

         do l=2,nlayer-1
            tmean=zt(1,l)
            if(zt(1,l).ne.zt(1,l-1))
     &        tmean=(zt(1,l)-zt(1,l-1))/log(zt(1,l)/zt(1,l-1))
              zlocal(l)= zlocal(l-1)
     &        -log(pplay(1,l)/pplay(1,l-1))*rnew(1,l)*tmean/g
         enddo
         zlocal(nlayer)= zlocal(nlayer-1)-
     &                   log(pplay(1,nlayer)/pplay(1,nlayer-1))*
     &                   rnew(1,nlayer)*tmean/g

      END IF       ! if(ngrid.ne.1)

      icount=icount+1
      RETURN
      END
