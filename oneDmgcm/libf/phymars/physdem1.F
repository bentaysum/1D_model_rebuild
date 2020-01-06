      subroutine physdem1(filename,lonfi,latfi,nsoil,nq,
     .                   phystep,day_ini,
     .                   time,tsurf,tsoil,co2ice,emis,q2,qsurf,
     .                   airefi,alb,ith,pzmea,pzstd,pzsig,pzgam,pzthe)

      implicit none
c-------------------------------------------------------------
c
c create physics (re-)start data file "restartfi.nc"
c
c
c
#include "dimensions.h"
#include "paramet.h"
c-----------------------------------------------------------------------
#include "comvert.h"
#include "comgeom2.h"
#include "control.h"
#include "comdissnew.h"
#include "logic.h"
#include "ener.h"
#include "netcdf.inc"
#include "dimphys.h"
#include"advtrac.h"
#include"callkeys.h"
c
      INTEGER nid,iq
      INTEGER, parameter :: ivap=1
      REAL, parameter :: qsolmax= 150.0
      character (len=*) :: filename
      character (len=7) :: str7

      REAL day_ini
      INTEGER nsoil,nq
      integer ierr,idim1,idim2,idim3,idim4,idim5,nvarid

c
      REAL phystep,time
      REAL latfi(ngridmx), lonfi(ngridmx)
!      REAL champhys(ngridmx)
      REAL tsurf(ngridmx)
      INTEGER length
      PARAMETER (length=100)
      REAL tab_cntrl(length)

c

!      EXTERNAL defrun_new,iniconst,geopot,inigeom,massdair,pression
!      EXTERNAL exner_hyb , SSUM
c
#include "serre.h"
#include "clesph0.h"
#include "fxyprim.h"
#include "comgeomfi.h"
#include "surfdat.h"
#include "comsoil.h"
#include "planete.h"
#include "dimradmars.h"
#include "yomaer.h"
#include "comcstfi.h"

      real co2ice(ngridmx),tsoil(ngridmx,nsoil),emis(ngridmx)
      real q2(ngridmx, llm+1),qsurf(ngridmx,nq)
      real airefi(ngridmx)
      real alb(ngridmx),ith(ngridmx,nsoil)
      real pzmea(ngridmx),pzstd(ngridmx)
      real pzsig(ngridmx),pzgam(ngridmx),pzthe(ngridmx)
      integer ig

! flag which identifies if we are using old tracer names (qsurf01,...)
      logical :: oldtracernames=.false.
      integer :: count
      character(len=30) :: txt ! to store some text
! indexes of water vapour & water ice tracers (if any):
      integer :: i_h2o_vap=0
      integer :: i_h2o_ice=0
c-----------------------------------------------------------------------

      ! copy airefi(:) to area(:)
      CALL SCOPY(ngridmx,airefi,1,area,1)
      ! note: area() is defined in comgeomfi.h

      DO ig=1,ngridmx
         albedodat(ig)=alb(ig) ! note: albedodat() is defined in surfdat.h
         zmea(ig)=pzmea(ig) ! note: zmea() is defined in surfdat.h
         zstd(ig)=pzstd(ig) ! note: zstd() is defined in surfdat.h
         zsig(ig)=pzsig(ig) ! note: zsig() is defined in surfdat.h
         zgam(ig)=pzgam(ig) ! note: zgam() is defined in surfdat.h
         zthe(ig)=pzthe(ig) ! note: zthe() is defined in surfdat.h
      ENDDO

      inertiedat(:,:)=ith(:,:) ! note inertiedat() is defined in comsoil.h
c
c  things to store in the physics start file:
c
      ierr = NF_CREATE(adjustl(filename),  
     &      IOR(NF_CLOBBER,NF_64BIT_OFFSET), nid)
      IF (ierr.NE.NF_NOERR) THEN
        WRITE(6,*)'physdem1: Problem creating file ',adjustl(filename)
        write(6,*) NF_STRERROR(ierr)
        CALL ABORT
      ENDIF
c
      ierr = NF_PUT_ATT_TEXT (nid, NF_GLOBAL, "title", 18,
     .                       "Physics start file")
c
      ierr = NF_DEF_DIM (nid,"index",length,idim1)
      if (ierr.ne.NF_NOERR) then
        WRITE(6,*)'physdem1: Problem defining index dimension'
	write(6,*) NF_STRERROR(ierr)
        call abort
      endif
c
      ierr = NF_DEF_DIM (nid,"physical_points",ngridmx,idim2)
      if (ierr.ne.NF_NOERR) then
        WRITE(6,*)'physdem1: Problem defining physical_points dimension'
	write(6,*) NF_STRERROR(ierr)
        call abort
      endif
c
      ierr = NF_DEF_DIM (nid,"subsurface_layers",nsoil,idim3)
      if (ierr.ne.NF_NOERR) then
      WRITE(6,*)'physdem1: Problem defining subsurface_layers dimension'
	write(6,*) NF_STRERROR(ierr)
        call abort
      endif
c
!      ierr = NF_DEF_DIM (nid,"nlayer+1",llm+1,idim4)
      ierr = NF_DEF_DIM (nid,"nlayer_plus_1",llm+1,idim4)
      if (ierr.ne.NF_NOERR) then
        WRITE(6,*)'physdem1: Problem defining nlayer+1 dimension'
        write(6,*) NF_STRERROR(ierr)
	call abort
      endif
c
      ierr = NF_DEF_DIM (nid,"number_of_advected_fields",nq,idim5)
      if (ierr.ne.NF_NOERR) then
        WRITE(6,*)'physdem1: Problem defining advected fields dimension'
        WRITE(6,*)' nq = ',nq,' and ierr = ', ierr
        write(6,*) NF_STRERROR(ierr)
      endif

      ierr = NF_ENDDEF(nid) ! exit NetCDF define mode

c clear tab_cntrl(:) array
      DO ierr = 1, length
         tab_cntrl(ierr) = 0.0
      ENDDO

      write(*,*) "physdem1: ngridmx: ",ngridmx

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Fill control array tab_cntrl(:) with paramleters for this run
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Informations on the physics grid
      tab_cntrl(1) = float(ngridmx)  ! number of nodes on physics grid
      tab_cntrl(2) = float(nlayermx) ! number of atmospheric layers
      tab_cntrl(3) = day_ini + int(time)         ! initial day 
      tab_cntrl(4) = time -int(time)            ! initiale time of day

c Informations about Mars, used by dynamics and physics
      tab_cntrl(5) = rad      ! radius of Mars (m) ~3397200
      tab_cntrl(6) = omeg     ! rotation rate (rad.s-1)
      tab_cntrl(7) = g        ! gravity (m.s-2) ~3.72
      tab_cntrl(8) = mugaz    ! Molar mass of the atmosphere (g.mol-1) ~43.49
      tab_cntrl(9) = rcp      !  = r/cp  ~0.256793 (=kappa dans dynamique)
      tab_cntrl(10) = daysec  ! length of a sol (s)  ~88775

      tab_cntrl(11) = phystep  ! time step in the physics
      tab_cntrl(12) = 0.
      tab_cntrl(13) = 0.

c Informations about Mars, only for physics
      tab_cntrl(14) = year_day  ! length of year (sols) ~668.6
      tab_cntrl(15) = periheli  ! min. Sun-Mars distance (Mkm) ~206.66
      tab_cntrl(16) = aphelie   ! max. SUn-Mars distance (Mkm) ~249.22
      tab_cntrl(17) = peri_day  ! date of perihelion (sols since N. spring)
      tab_cntrl(18) = obliquit  ! Obliquity of the planet (deg) ~23.98

c Boundary layer and turbulence
      tab_cntrl(19) = z0_default   ! default surface roughness (m) ~0.01
      tab_cntrl(20) = lmixmin   ! mixing length ~100
      tab_cntrl(21) = emin_turb ! minimal energy ~1.e-8

c Optical properties of polar caps and ground emissivity
      tab_cntrl(22) = albedice(1)  ! Albedo of northern cap ~0.5
      tab_cntrl(23) = albedice(2)  ! Albedo of southern cap ~0.5
      tab_cntrl(24) = emisice(1)   ! Emissivity of northern cap ~0.95
      tab_cntrl(25) = emisice(2)   ! Emissivity of southern cap ~0.95
      tab_cntrl(26) = emissiv      ! Emissivity of martian soil ~.95
      tab_cntrl(31) = iceradius(1) ! mean scat radius of CO2 snow (north)
      tab_cntrl(32) = iceradius(2) ! mean scat radius of CO2 snow (south)
      tab_cntrl(33) = dtemisice(1) ! time scale for snow metamorphism (north)
      tab_cntrl(34) = dtemisice(2) ! time scale for snow metamorphism (south)

c dust aerosol properties
      tab_cntrl(27) = tauvis      ! mean visible optical depth

      tab_cntrl(28) = 0. 
      tab_cntrl(29) = 0.
      tab_cntrl(30) = 0.

! Soil properties:
      tab_cntrl(35) = volcapa ! soil volumetric heat capacity
      
c
!      write(*,*) "physdem1: tab_cntrl():",tab_cntrl
      
      ierr = NF_REDEF (nid) ! Enter NetCDF (re-)define mode
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'physdem1: Failed to swich to NetCDF define mode'
         CALL abort
      ENDIF
      ! define variable
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "controle", NF_DOUBLE, 1, idim1,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "controle", NF_FLOAT, 1, idim1,nvarid)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'physdem1: Failed to define controle'
         CALL abort
      ENDIF
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 18,
     .                        "Control parameters")
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'physdem1: Failed to define controle title attribute'
         CALL abort
      ENDIF
      ierr = NF_ENDDEF(nid) ! Leave NetCDF define mode
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'physdem1: Failed to swich out of NetCDF define mode'
         CALL abort
      ENDIF
      ! write variable
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tab_cntrl)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
#endif
      IF (ierr.NE.NF_NOERR) THEN
         PRINT*, 'physdem1: Failed to write controle data'
         CALL abort
      ENDIF

! write mid-layer depths mlayer() !known from comsoil.h

      ierr = NF_REDEF (nid) ! Enter NetCDF (re-)define mode
      ! define variable
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"soildepth",NF_DOUBLE,1,idim3,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"soildepth",NF_FLOAT,1,idim3,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 20,
     .                        "Soil mid-layer depth")
      ierr = NF_ENDDEF(nid) ! Leave NetCDF define mode
      ! write variable
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,mlayer)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,mlayer)
#endif

c

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "longitude", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "longitude", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 26,
     .                        "Longitudes of physics grid")
      ierr = NF_ENDDEF(nid)

#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,lonfi)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,lonfi)
#endif

c

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "latitude", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "latitude", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 25,
     .                        "Latitudes of physics grid")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,latfi)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,latfi)
#endif

c

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "area", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "area", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9,
     .                        "Mesh area")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,area)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,area)
#endif

c

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "phisfi", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "phisfi", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 27,
     .                        "Geopotential at the surface")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,phisfi)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,phisfi)
#endif

c

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "albedodat", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "albedodat", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 21,
     .                        "Albedo of bare ground")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,albedodat)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,albedodat)
#endif

c
c   some data for Francois Lott's programs
c

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZMEA", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZMEA", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19,
     .                        "Relief: mean relief")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zmea)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zmea)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZSTD", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZSTD", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 26,
     .                        "Relief: standard deviation")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zstd)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zstd)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZSIG", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZSIG", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                        "Relief: sigma parameter")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zsig)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zsig)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZGAM", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZGAM", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                        "Relief: gamma parameter")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zgam)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zgam)
#endif
c
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "ZTHE", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "ZTHE", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                        "Relief: theta parameter")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,zthe)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,zthe)
#endif

c Write the physical fields

! CO2 Ice Cover

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "co2ice", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "co2ice", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 13,
     .                        "CO2 ice cover")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,co2ice)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,co2ice)
#endif

! Soil Thermal inertia

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"inertiedat",NF_DOUBLE,
     &                   2,(/idim2,idim3/),nvarid)
#else
      ierr = NF_DEF_VAR (nid,"inertiedat",NF_FLOAT,
     &                   2,(/idim2,idim3/),nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 20,
     .                        "Soil thermal inertia")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,inertiedat)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,inertiedat)
#endif

!  Surface temperature

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "tsurf", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "tsurf", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19,
     .                        "Surface temperature")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tsurf)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tsurf)
#endif

! Soil temperature

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"tsoil",NF_DOUBLE,2,(/idim2,idim3/),nvarid)
#else
!      ierr = NF_DEF_VAR (nid, "tsoil", NF_FLOAT, 2, idim2,nvarid)
      ierr = NF_DEF_VAR (nid,"tsoil",NF_FLOAT,2,(/idim2,idim3/),nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 16,
     .                        "Soil temperature")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tsoil)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tsoil)
#endif

c emissivity

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "emis", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "emis", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 18,
     .                        "Surface emissivity")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,emis)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,emis)
#endif

! surface roughness length (z0 is a common in surfdat.h)

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "z0", NF_DOUBLE, 1, idim2,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "z0", NF_FLOAT, 1, idim2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 24,
     .                        "Surface roughness length")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,z0)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,z0)
#endif


c planetary boundary layer

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"q2",NF_DOUBLE,2,(/idim2,idim4/),nvarid)
#else
      ierr = NF_DEF_VAR (nid,"q2",NF_FLOAT, 2,(/idim2,idim4/),nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 17,
     .                        "pbl wind variance")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,q2)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,q2)
#endif
      IF (ierr.NE.NF_NOERR) THEN
        PRINT*, 'physdem1: Failed to write q2'
        CALL abort
      ENDIF

c tracers

! Preliminary stuff: check if tracers follow old naming convention (qsurf01,
!                    qsurf02, ...)
      count=0
      do iq=1,nqmx
        txt= " "
        write(txt,'(a1,i2.2)')'q',iq
        if (txt.ne.tnom(iq)) then ! use tracer names stored in dynamics
          ! did not find old tracer name
          exit ! might as well stop here
        else
          ! found old tracer name
          count=count+1
        endif
      enddo
      if (count.eq.nqmx) then
        write(*,*) "physdem1:tracers seem to follow old naming ",
     &             "convention (qsurf01,qsurf02,...)"
        write(*,*) "   => will work for now ... "
        write(*,*) "      but you should run newstart to rename them"
        oldtracernames=.true.
        ! Moreover, if computing water cycle with ice, move surface ice
        ! back to qsurf(nqmx)
        IF (water) THEN
          !"loop" to avoid potential out-of-bounds on arrays
          write(*,*)'physdem1: moving surface water ice to index ',nqmx
          do iq=nqmx,nqmx
          qsurf(1:ngridmx,iq)=qsurf(1:ngridmx,iq-1)
          qsurf(1:ngridmx,iq-1)=0
          enddo
        ENDIF
      endif ! of if (count.eq.nqmx)

      IF(nq.GE.1) THEN
! preliminary stuff: look for water vapour & water ice tracers (if any)
        if (.not.oldtracernames) then
         do iq=1,nq
           if (tnom(iq).eq."h2o_vap") then
             i_h2o_vap=iq
           endif
           if (tnom(iq).eq."h2o_ice") then
             i_h2o_ice=iq
           endif
         enddo ! of iq=1,nq
         ! handle special case of only water vapour tracer (no ice)
         if ((i_h2o_vap.ne.0).and.(i_h2o_ice.eq.0)) then
          ! then the index of (surface) ice is i_h2o_vap
          i_h2o_ice=i_h2o_vap
         endif
        endif ! of if (.not.oldtracernames) 

         DO iq=1,nq
           IF (oldtracernames) THEN
             txt=" "
             write(txt,'(a5,i2.2)')'qsurf',iq
           ELSE
             txt=tnom(iq)
             ! Exception: there is no water vapour surface tracer
             if (txt.eq."h2o_vap") then
               write(*,*)"physdem1: skipping water vapour tracer"
               if (i_h2o_ice.eq.i_h2o_vap) then
               ! then there is no "water ice" tracer; but still
               ! there is some water ice on the surface
                 write(*,*)"          writting water ice instead"
                 txt="h2o_ice"
               else
               ! there is a "water ice" tracer which has been / will be
               ! delt with in due time
                 cycle
               endif ! of if (igcm_h2o_ice.eq.igcm_h2o_vap)
             endif ! of if (txt.eq."h2o_vap")
           ENDIF ! of IF (oldtracernames)

           ierr=NF_REDEF(nid)
           IF (ierr.NE.NF_NOERR) THEN
             PRINT*, 'physdem1: Failed to swich to NetCDF define mode'
             CALL abort
           ENDIF
#ifdef NC_DOUBLE
           ierr=NF_DEF_VAR(nid,txt,NF_DOUBLE,1,idim2,nvarid)
#else
           ierr=NF_DEF_VAR(nid,txt,NF_FLOAT,1,idim2,nvarid)
#endif
           IF (ierr.NE.NF_NOERR) THEN
             PRINT*, 'physdem1: Failed to define ',trim(txt)
             CALL abort
           ENDIF
           ierr=NF_PUT_ATT_TEXT (nid, nvarid, "title", 17,
     &                        "tracer on surface")
           IF (ierr.NE.NF_NOERR) THEN
             PRINT*, 'physdem1: Failed to define ',trim(txt),
     &               ' title attribute'
             CALL abort
           ENDIF
           ierr=NF_ENDDEF(nid)
           IF (ierr.NE.NF_NOERR) THEN
             PRINT*, 'physdem1: Failed to swich out of define mode'
             CALL abort
           ENDIF
           
#ifdef NC_DOUBLE
            ierr=NF_PUT_VAR_DOUBLE (nid,nvarid,qsurf(1,iq))
#else
            ierr=NF_PUT_VAR_REAL (nid,nvarid,qsurf(1,iq))
#endif
           IF (ierr.NE.NF_NOERR) THEN
             PRINT*, 'physdem1: Failed to write ',trim(txt)
             CALL abort
           ENDIF
         ENDDO ! of DO iq=1,nq
      ENDIF ! of IF(nq.GE.1)

c close file
      ierr = NF_CLOSE(nid)

      RETURN

      END
