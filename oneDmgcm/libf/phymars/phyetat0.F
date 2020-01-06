      SUBROUTINE phyetat0 (fichnom,tab0,Lmodif,nsoil,nq,
     .           day_ini,time,
     .           tsurf,tsoil,emis,q2,qsurf,co2ice)

      use netcdf

      implicit none
c======================================================================
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c  Adaptation à Mars : Yann Wanherdrick 
c Objet: Lecture de l etat initial pour la physique
! Modifs: Aug.2010 EM : use NetCDF90 to load variables (enables using
!                      r4 or r8 restarts independently of having compiled
!                      the GCM in r4 or r8)
c======================================================================
!#include "netcdf.inc"
#include "dimensions.h"
#include "dimphys.h"
!#include "comgeomfi.h"
#include "surfdat.h"
#include "planete.h"
#include "dimradmars.h"
#include "yomaer.h"
#include "comcstfi.h"
!#include "tracer.h"
#include"advtrac.h"
c======================================================================
      INTEGER nbsrf !Mars nbsrf a 1 au lieu de 4
      PARAMETER (nbsrf=1) ! nombre de sous-fractions pour une maille
!======================================================================
!  Arguments:
!  ---------
!  inputs:
      character*(*) fichnom ! "startfi.nc" file
      integer tab0
      integer Lmodif
      integer nsoil ! # of soil layers
      integer nq
      integer day_ini
      real time

!  outputs:
      real tsurf(ngridmx) ! surface temperature
      real tsoil(ngridmx,nsoil) ! soil temperature
      real emis(ngridmx) ! surface emissivity
      real q2(ngridmx, llm+1) ! 
      real qsurf(ngridmx,nq) ! tracers on surface
      real co2ice(ngridmx) ! co2 ice cover

!======================================================================
!  Local variables:

      real surffield(ngridmx) ! to temporarily store a surface field
      real xmin,xmax ! to display min and max of a field
c
      INTEGER ig,iq,lmax
      INTEGER nid, nvarid
      INTEGER ierr, i, nsrf
!      integer isoil 
!      INTEGER length
!      PARAMETER (length=100)
      CHARACTER*7 str7
      CHARACTER*2 str2
      CHARACTER*1 yes
c
      REAL p_rad,p_omeg,p_g,p_mugaz,p_daysec
      INTEGER nqold

! flag which identifies if 'startfi.nc' file is using old names (qsurf01,...)
      logical :: oldtracernames=.false.
      integer :: count
      character(len=30) :: txt ! to store some text

 
c
c Ouvrir le fichier contenant l etat initial:
c

      ierr=nf90_open(fichnom,NF90_NOWRITE,nid)
      IF (ierr.NE.nf90_noerr) THEN
        write(6,*)' Pb d''ouverture du fichier ',trim(fichnom)
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF

! Preliminary stuff: check if tracers follow old naming convention (qsurf01,
!                    qsurf02, ...)
      count=0
      do iq=1,nqmx
        txt= " "
        write(txt,'(a5,i2.2)')'qsurf',iq
        ierr=nf90_inq_varid(nid,txt,nvarid)
        if (ierr.ne.nf90_noerr) then
          ! did not find old tracer name
          exit ! might as well stop here
        else
          ! found old tracer name
          count=count+1
        endif
      enddo
      if (count.eq.nqmx) then
        write(*,*) "phyetat0:tracers seem to follow old naming ",
     &             "convention (qsurf01,qsurf02,...)"
        write(*,*) "   => will work for now ... "
        write(*,*) "      but you should run newstart to rename them"
        oldtracernames=.true.
      endif

c modifications possibles des variables de tab_cntrl
      PRINT*
      write(*,*) 'TABFI de phyeta0',Lmodif,tab0
      call tabfi (nid,Lmodif,tab0,day_ini,lmax,p_rad,
     .              p_omeg,p_g,p_mugaz,p_daysec,time)
c
c Read latitudes (coordinates): No need, these are provided by the dynamics
c
!      ierr=nf90_inq_varid(nid,"latitude",nvarid)
!      IF (ierr.NE.nf90_noerr) THEN
!         PRINT*, 'phyetat0: Le champ <latitude> est absent'
!         write(*,*)trim(nf90_strerror(ierr))
!         CALL abort
!      ENDIF
!      ierr=nf90_get_var(nid,nvarid,lati)
!      IF (ierr.NE.nf90_noerr) THEN
!         PRINT*, 'phyetat0: Lecture echouee pour <latitude>'
!         write(*,*)trim(nf90_strerror(ierr))
!         CALL abort
!      ENDIF
c
c read longitudes (coordinates): No need, these are provided by the dynamics
c
!      ierr=nf90_inq_varid(nid,"longitude",nvarid)
!      IF (ierr.NE.nf90_noerr) THEN
!         PRINT*, 'phyetat0: Le champ <longitude> est absent'
!         write(*,*)trim(nf90_strerror(ierr))
!         CALL abort
!      ENDIF
!      ierr=nf90_get_var(nid,nvarid,long)
!      IF (ierr.NE.nf90_noerr) THEN
!         PRINT*, 'phyetat0: Lecture echouee pour <longitude>'
!         write(*,*)trim(nf90_strerror(ierr))
!         CALL abort
!      ENDIF
c
c Read areas of meshes: No need, these are provided by the dynamics
c
!      ierr=nf90_inq_varid(nid,"area",nvarid)
!      IF (ierr.NE.nf90_noerr) THEN
!         PRINT*, 'phyetat0: Le champ <area> est absent'
!         write(*,*)trim(nf90_strerror(ierr))
!         CALL abort
!      ENDIF
!      ierr=nf90_get_var(nid,nvarid,area)
!      IF (ierr.NE.nf90_noerr) THEN
!         PRINT*, 'phyetat0: Lecture echouee pour <area>'
!         write(*,*)trim(nf90_strerror(ierr))
!         CALL abort
!      ENDIF
!      xmin = 1.0E+20
!      xmax = -1.0E+20
!      xmin = MINVAL(area)
!      xmax = MAXVAL(area)
!      PRINT*,'Aires des mailles <area>:', xmin, xmax
c
c Lecture du geopotentiel au sol:
c
      ierr=nf90_inq_varid(nid,"phisfi",nvarid)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Le champ <phisfi> est absent'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,phisfi)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <phisfi>'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(phisfi)
      xmax = MAXVAL(phisfi)
      PRINT*,'Geopotentiel au sol <phisfi>:', xmin, xmax
c
c Lecture de l''albedo du sol nu:
c
      ierr=nf90_inq_varid(nid,"albedodat",nvarid)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Le champ <albedodat> est absent'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,albedodat)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <albedodat>'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(albedodat)
      xmax = MAXVAL(albedodat)
      PRINT*,'Albedo du sol nu <albedodat>:', xmin, xmax
c
c Lecture de l''inertie thermique du sol:
c
!      ierr = NF_INQ_VARID (nid, "inertiedat", nvarid)
!      IF (ierr.NE.nf90_noerr) THEN
!         PRINT*, 'phyetat0: Le champ <inertiedat> est absent'
!         CALL abort
!      ENDIF
!#ifdef NC_DOUBLE
!      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, inertiedat)
!#else
!      ierr = NF_GET_VAR_REAL(nid, nvarid, inertiedat)
!#endif
!      IF (ierr.NE.nf90_noerr) THEN
!         PRINT*, 'phyetat0: Lecture echouee pour <inertiedat>'
!         CALL abort
!      ENDIF
!      xmin = 1.0E+20
!      xmax = -1.0E+20
!      xmin = MINVAL(inertiedat)
!      xmax = MAXVAL(inertiedat)
!      PRINT*,'Inertie thermique du sol <inertiedat>:', xmin, xmax
c
c ZMEA
c
      ierr=nf90_inq_varid(nid,"ZMEA",nvarid)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Le champ <ZMEA> est absent'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,zmea)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZMEA>'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zmea(i),xmin)
         xmax = MAX(zmea(i),xmax)
      ENDDO
      PRINT*,'<zmea>:', xmin, xmax
c
c ZSTD
c
      ierr=nf90_inq_varid(nid,"ZSTD",nvarid)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Le champ <ZSTD> est absent'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,zstd)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZSTD>'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zstd(i),xmin)
         xmax = MAX(zstd(i),xmax)
      ENDDO
      PRINT*,'<zstd>:', xmin, xmax
c
c ZSIG
c
      ierr=nf90_inq_varid(nid,"ZSIG",nvarid)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Le champ <ZSIG> est absent'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,zsig)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZSIG>'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zsig(i),xmin)
         xmax = MAX(zsig(i),xmax)
      ENDDO
      PRINT*,'<zsig>:', xmin, xmax
c
c ZGAM
c
      ierr=nf90_inq_varid(nid,"ZGAM",nvarid)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Le champ <ZGAM> est absent'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,zgam)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZGAM>'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zgam(i),xmin)
         xmax = MAX(zgam(i),xmax)
      ENDDO
      PRINT*,'<zgam>:', xmin, xmax
c
c ZTHE
c
      ierr=nf90_inq_varid(nid,"ZTHE",nvarid)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Le champ <ZTHE> est absent'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,zthe)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <ZTHE>'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, ngridmx
         xmin = MIN(zthe(i),xmin)
         xmax = MAX(zthe(i),xmax)
      ENDDO
      PRINT*,'<zthe>:', xmin, xmax
c
c CO2 ice cover
c
      ierr=nf90_inq_varid(nid,"co2ice",nvarid)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Le champ <co2ice> est absent'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,co2ice)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <co2ice>'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(co2ice)
      xmax = MAXVAL(co2ice)
      PRINT*,'CO2 ice cover <co2ice>:', xmin, xmax
c
c Lecture des temperatures du sol:
c
      ierr=nf90_inq_varid(nid,"tsurf",nvarid)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Le champ <tsurf> est absent'
! Ehouarn: dump this part (causes problems with netcdf90 interface)
!         PRINT*, '          Mais je vais essayer de lire TS**'
!         IF (nbsrf.GT.99) THEN
!            PRINT*, "Trop de sous-mailles"
!            CALL abort
!         ENDIF
!         DO nsrf = 1, nbsrf
!           WRITE(str2,'(i2.2)') nsrf
!           ierr=nf90_inq_varid(nid,"TS"//str2,nvarid)
!           IF (ierr.NE.nf90_noerr) THEN
!              PRINT*, "phyetat0: Le champ <TS"//str2//"> est absent"
!              write(*,*)trim(nf90_strerror(ierr))
!              CALL abort
!           ENDIF
!           ierr=nf90_get_var(nid,nvarid,tsurf(1,nsrf))
!           IF (ierr.NE.nf90_noerr) THEN
!             PRINT*, "phyetat0: Lecture echouee pour <TS"//str2//">"
!             write(*,*)trim(nf90_strerror(ierr))
!             CALL abort
!           ENDIF
!           xmin = 1.0E+20
!           xmax = -1.0E+20
!           xmin = MINVAL(tsurf)
!           xmax = MAXVAL(tsurf)
!           PRINT*,'Temperature du sol TS**:', nsrf, xmin, xmax
!         ENDDO
      ELSE
         PRINT*, 'phyetat0: Le champ <tsurf> est present'
         PRINT*, '          J ignore donc les autres temperatures TS**'
         !ierr=nf90_get_var(nid,nvarid,tsurf(1,1))
         ierr=nf90_get_var(nid,nvarid,tsurf)
         IF (ierr.NE.nf90_noerr) THEN
            PRINT*, "phyetat0: Lecture echouee pour <TSURF>"
            write(*,*)trim(nf90_strerror(ierr))
            CALL abort
         ENDIF
         xmin = 1.0E+20
         xmax = -1.0E+20
         xmin = MINVAL(tsurf)
         xmax = MAXVAL(tsurf)
         PRINT*,'Temperature du sol <tsurf>', xmin, xmax
! Ehouarn: Mars => nbsurf=1
!         IF (nbsrf >= 2) THEN
!            DO nsrf = 2, nbsrf
!               DO i = 1, ngridmx
!                  tsurf(i,nsrf) = tsurf(i,1)
!               ENDDO
!            ENDDO
!         ENDIF
      ENDIF
c
c Lecture des temperatures du sol profond:
c
!      IF (nsoil.GT.99 .OR. nbsrf.GT.99) THEN
!         PRINT*, "Trop de couches ou sous-mailles"
!         CALL abort
!      ENDIF
!      DO nsrf = 1, nbsrf
!         DO isoil=1, nsoil
!            WRITE(str7,'(i2.2,"srf",i2.2)') isoil, nsrf
!            ierr = NF_INQ_VARID (nid, 'tsoil', nvarid)
!            IF (ierr.NE.nf90_noerr) THEN
!               PRINT*, "phyetat0: Le champ <tsoil> est absent"
!               PRINT*, "          Il prend donc la valeur de surface"
!               DO i=1, ngridmx
!                  tsoil(i,isoil,nsrf)=tsurf(i,nsrf)
!               ENDDO
!            ELSE
!#ifdef NC_DOUBLE
!              ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tsoil(1,1,nsrf))
!#else
!              ierr = NF_GET_VAR_REAL(nid, nvarid, tsoil(1,1,nsrf))
!#endif
!               IF (ierr.NE.nf90_noerr) THEN
!                  PRINT*, "Lecture echouee pour <tsoil>"
!                  CALL abort
!               ENDIF
!            ENDIF
!         ENDDO
!      ENDDO
!      xmin = 1.0E+20
!      xmax = -1.0E+20
!      xmin = MINVAL(tsoil)
!      xmax = MAXVAL(tsoil)
!      PRINT*,'Temperatures du sol profond <tsoil>', xmin, xmax
c
c Surface emissivity
c
      ierr=nf90_inq_varid(nid,"emis",nvarid)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Le champ <emis> est absent'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,emis)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <emis>'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(emis)
      xmax = MAXVAL(emis)
      PRINT*,'Surface emissivity <emis>:', xmin, xmax

!
! surface roughness length (NB: z0 is a common in surfdat.h)
!
      ierr=nf90_inq_varid(nid,"z0",nvarid)
      if (ierr.ne.nf90_noerr) then
        write(*,*) 'phyetat0: did not find z0 field!'
        write(*,*) 'will use constant value of z0_default instead'
        z0(:)=z0_default
      else
        ierr=nf90_get_var(nid,nvarid,z0)
        if (ierr.ne.nf90_noerr) then
          write(*,*) 'phyetat0: Failed loading <z0>'
          write(*,*)trim(nf90_strerror(ierr))
          call abort
        endif
        xmin=1.0E+20
        xmax=-1.0E+20
        xmin=minval(z0)
        xmax=maxval(z0)
        write(*,*)'Surface roughness <z0>:',xmin,xmax
      endif


c
c pbl wind variance
c
      ierr=nf90_inq_varid(nid,"q2",nvarid)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Le champ <q2> est absent'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,q2)
      IF (ierr.NE.nf90_noerr) THEN
         PRINT*, 'phyetat0: Lecture echouee pour <q2>'
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(q2)
      xmax = MAXVAL(q2)
      PRINT*,'pbl wind variance <q2>:', xmin, xmax
c
c tracer on surface
c

      IF(nq.GE.1) THEN
         nqold=nq
         DO iq=1,nq
!            str7(1:5)='qsurf'
!            WRITE(str7(6:7),'(i2.2)') iq
!            ierr = NF_INQ_VARID (nid,str7,nvarid)
           IF (oldtracernames) THEN
             txt=" "
             write(txt,'(a5,i2.2)')'qsurf',iq
           ELSE
             txt=tnom(iq)
             if (txt.eq."h2o_vap") then
               ! There is no surface tracer for h2o_vap;
               ! "h2o_ice" should be loaded instead
               txt="h2o_ice"
               write(*,*) 'phyetat0: loading surface tracer',
     &                     ' h2o_ice instead of h2o_vap'
             endif
           ENDIF ! of IF (oldtracernames) THEN
           ierr=nf90_inq_varid(nid,txt,nvarid)
           IF (ierr.NE.nf90_noerr) THEN
             write(*,*) 'PHYETAT0: WARNING : surface tracer',trim(txt),
     &                  ' not found in file'
             write(*,*) trim(txt), ' set to 0'
             do ig=1,ngridmx
               qsurf(ig,iq)=0.
             end do
             nqold=min(iq-1,nqold)
           ELSE
           !ierr=nf90_get_var(nid,nvarid,qsurf(1,iq))
           ierr=nf90_get_var(nid,nvarid,surffield)
           qsurf(1:ngridmx,iq)=surffield(1:ngridmx)
             IF (ierr.NE.nf90_noerr) THEN
               PRINT*, 'phyetat0: Lecture echouee pour <',trim(txt),'>'
               write(*,*)trim(nf90_strerror(ierr))
               CALL abort
             ENDIF
           ENDIF
           xmin = 1.0E+20
           xmax = -1.0E+20
           xmin = MINVAL(qsurf(1:ngridmx,iq))
           xmax = MAXVAL(qsurf(1:ngridmx,iq))
           PRINT*,'tracer on surface <',trim(txt),'>:',xmin,xmax
         ENDDO
         if ((nqold.lt.nq).and.(nqold.ge.1)) then
c        case when new tracer are added in addition to old ones
             write(*,*)'qsurf 1 to ', nqold,'were already present'
             write(*,*)'qsurf ', nqold+1,' to ', nqmx,'are new'
!            yes=' '
!            do while ((yes.ne.'y').and.(yes.ne.'n'))
!             write(*,*) 'Would you like to reindex qsurf # 1 ->',nqold
!             write(*,*) 'to #',nqmx-nqold+1,'->', nqmx,'   (y or n) ?'
!             read(*,fmt='(a)') yes
!            end do
!            if (yes.eq.'y') then
!              write(*,*) 'OK, let s reindex qsurf'
!                 do ig=1,ngridmx
!                    do iq=nqmx,nqmx-nqold+1,-1
!                       qsurf(ig,iq)=qsurf(ig,iq-nqmx+nqold)
!                    end do
!                    do iq=nqmx-nqold,1,-1
!                       qsurf(ig,iq)= 0.
!                    end do
!                 end do
!            end if
         end if
      ENDIF

! Call to soil_settings, in order to read soil temperatures,
! as well as thermal inertia and volumetric heat capacity

      call soil_settings(nid,ngridmx,nsoil,tsurf,tsoil)
c
c Fermer le fichier:
c
      ierr=nf90_close(nid)
c
      RETURN
      END
