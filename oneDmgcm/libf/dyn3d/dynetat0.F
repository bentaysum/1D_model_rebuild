      SUBROUTINE dynetat0(fichnom,nq,vcov,ucov,
     .                    teta,q,masse,ps,phis,time)
      
      use netcdf
      
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van / L.Fairhead
c   -------
c
c   objet:
c   ------
c
c   Lecture de l'etat initial
c
c   Modifs: Oct.2008 read in tracers by name. Ehouarn Millour
c           Aug.2010 use NetCDF90 to load variables (enables using
c                    r4 or r8 restarts independently of having compiled
c                    the GCM in r4 or r8)
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "temps.h"
#include "comconst.h"
#include "comvert.h"
#include "comgeom2.h"
#include "ener.h"
#include "description.h"
#include "serre.h"
#include "logic.h"
#include"advtrac.h"

c   Arguments:
c   ----------

      CHARACTER*(*) fichnom
      INTEGER nq
      REAL vcov(iip1,jjm,llm),ucov(iip1,jjp1,llm),teta(iip1,jjp1,llm)
      REAL q(iip1,jjp1,llm,nq)
      REAL masse(iip1,jjp1,llm)
      REAL ps(iip1,jjp1),phis(iip1,jjp1)

      REAL time ! fraction of day the fields correspond to

c   Variables 
c
      REAL traceur(iip1,jjp1,llm) ! to temporarily store a tracer
      INTEGER length,iq,i,j,l
      PARAMETER (length = 100)
      REAL tab_cntrl(length) ! tableau des parametres du run
      INTEGER ierr, nid, nvarid, nqold
      CHARACTER  str3*3,yes*1

c-----------------------------------------------------------------------

c  Ouverture NetCDF du fichier etat initial

      ierr=nf90_open(fichnom,NF90_NOWRITE,nid)
      IF (ierr.NE.nf90_noerr) THEN
        write(6,*)' Pb d''ouverture du fichier ',trim(fichnom)
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF

c
      ierr=nf90_inq_varid(nid,"controle",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <controle> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,tab_cntrl)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echoue pour <controle>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      im         = tab_cntrl(1)
      jm         = tab_cntrl(2)
      lllm       = tab_cntrl(3)
      day_ini    = tab_cntrl(4)
      rad        = tab_cntrl(5)
      omeg       = tab_cntrl(6)
      g          = tab_cntrl(7)
      cpp        = tab_cntrl(8)
      kappa      = tab_cntrl(9)
      daysec     = tab_cntrl(10)
      dtvr       = tab_cntrl(11)
      etot0      = tab_cntrl(12)
      ptot0      = tab_cntrl(13)
      ztot0      = tab_cntrl(14)
      stot0      = tab_cntrl(15)
      ang0       = tab_cntrl(16)
      pa         = tab_cntrl(17)
      preff      = tab_cntrl(18)
c
      clon       = tab_cntrl(19)
      clat       = tab_cntrl(20)
      grossismx  = tab_cntrl(21)
      grossismy  = tab_cntrl(22)
c
      IF ( tab_cntrl(23).EQ.1. )  THEN
        fxyhypb  = . TRUE .
        dzoomx   = tab_cntrl(24)
        dzoomy   = tab_cntrl(25)
        taux     = tab_cntrl(27)
        tauy     = tab_cntrl(28)
      ELSE
        fxyhypb = . FALSE .
        ysinus  = . FALSE .
        IF( tab_cntrl(26).EQ.1. ) ysinus = . TRUE. 
      ENDIF
c   .................................................................
c
c
      PRINT*,'dynetat0: rad,omeg,g,cpp,kappa',rad,omeg,g,cpp,kappa
 
      IF(   im.ne.iim           )  THEN
          PRINT 1,im,iim
          STOP
      ELSE  IF( jm.ne.jjm       )  THEN
          PRINT 2,jm,jjm
          STOP
      ELSE  IF( lllm.ne.llm     )  THEN
          PRINT 3,lllm,llm
          STOP
      ENDIF

      ierr=nf90_inq_varid(nid,"rlonu",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <rlonu> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,rlonu)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <rlonu>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_inq_varid(nid,"rlatu",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <rlatu> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,rlatu)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <rlatu>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_inq_varid(nid,"rlonv",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <rlonv> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,rlonv)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <rlonv>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_inq_varid(nid,"rlatv",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <rlatv> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,rlatv)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour rlatv"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_inq_varid(nid,"cu",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <cu> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,cu)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <cu>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_inq_varid(nid,"cv",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <cv> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,cv)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <cv>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_inq_varid(nid,"aire",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <aire> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,aire)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <aire>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_inq_varid(nid,"phisinit",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <phisinit> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,phis)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <phisinit>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_inq_varid(nid,"Time",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
	     ierr=nf90_inq_varid(nid,"temps", nvarid)
		 IF (ierr .NE. nf90_noerr) THEN
           PRINT*, "dynetat0: <Time> or <temps> absent"
           write(*,*)trim(nf90_strerror(ierr))
           CALL abort
         ENDIF
      ENDIF
      ierr=nf90_get_var(nid,nvarid,time)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee <Time>/<temps>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_inq_varid(nid,"ucov",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <ucov> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,ucov)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <ucov>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
 
      ierr=nf90_inq_varid(nid,"vcov",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <vcov> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,vcov)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <vcov>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_inq_varid(nid,"teta",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <teta> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,teta)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <teta>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF


      IF(nq.GE.1) THEN
        write(*,*) 'dynetat0: loading tracers'
         IF(nq.GT.99) THEN
            PRINT*, "Trop de traceurs"
            CALL abort
         ENDIF
         nqold=nq
         DO iq=1,nq
!           str3(1:1)='q'
!           WRITE(str3(2:3),'(i2.2)') iq
!           ierr =  NF_INQ_VARID (nid, str3, nvarid)
! NB: tracers are now read in using their name ('tnom' from advtrac.h)
!           write(*,*) "  loading tracer:",trim(tnom(iq))
           ierr=nf90_inq_varid(nid,tnom(iq),nvarid)
           IF (ierr .NE. nf90_noerr) THEN
!              PRINT*, "dynetat0: Le champ <"//str3//"> est absent"
              PRINT*, "dynetat0: Le champ <"//trim(tnom(iq))//
     &                "> est absent"
              PRINT*, "          Il est donc initialise a zero"
              CALL initial0(ijp1llm,q(1,1,1,iq))
              nqold=min(iq-1,nqold)
           ELSE
           !ierr=nf90_get_var(nid,nvarid,q(1,1,1,iq))
           ierr=nf90_get_var(nid,nvarid,traceur)
             IF (ierr .NE. nf90_noerr) THEN
!                 PRINT*, "dynetat0: Lecture echouee pour "//str3
               PRINT*, "dynetat0: Lecture echouee pour "//trim(tnom(iq))
               CALL abort
             ENDIF
             q(:,:,:,iq)=traceur(:,:,:)
           ENDIF
         ENDDO
         if ((nqold.lt.nq).and.(nqold.ge.1)) then   
c        case when new tracer are added in addition to old ones
             write(*,*)'tracers 1 to ', nqold,'were already present'
             write(*,*)'tracers ', nqold+1,' to ', nqmx,'are new'
!             yes=' '
!            do while ((yes.ne.'y').and.(yes.ne.'n'))
!             write(*,*) 'Would you like to reindex tracer # 1 ->',nqold
!             write(*,*) 'to #',nqmx-nqold+1,'->', nqmx,'   (y or n) ?'
!             read(*,fmt='(a)') yes
!            end do
!            if (yes.eq.'y') then
!              write(*,*) 'OK, let s reindex the tracers'
!              do l=1,llm
!                do j=1,jjp1
!                  do i=1,iip1
!                    do iq=nqmx,nqmx-nqold+1,-1
!                       q(i,j,l,iq)=q(i,j,l,iq-nqmx+nqold)   
!                    end do
!                    do iq=nqmx-nqold,1,-1
!                       q(i,j,l,iq)= 0.
!                    end do
!                  end do
!                end do
!              end do
!            end if 
         end if
      ENDIF

      ierr=nf90_inq_varid(nid,"masse",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <masse> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,masse)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <masse>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_inq_varid(nid,"ps",nvarid)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Le champ <ps> est absent"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF
      ierr=nf90_get_var(nid,nvarid,ps)
      IF (ierr .NE. nf90_noerr) THEN
         PRINT*, "dynetat0: Lecture echouee pour <ps>"
         write(*,*)trim(nf90_strerror(ierr))
         CALL abort
      ENDIF

      ierr=nf90_close(nid)

      ! day_ini=day_ini+INT(time) ! obsolete stuff ; 0<time<1 anyways
      ! time=time-INT(time)

  1   FORMAT(//10x,'la valeur de im =',i4,2x,'lue sur le fichier de dem
     *arrage est differente de la valeur parametree iim =',i4//)
   2  FORMAT(//10x,'la valeur de jm =',i4,2x,'lue sur le fichier de dem
     *arrage est differente de la valeur parametree jjm =',i4//)
   3  FORMAT(//10x,'la valeur de lmax =',i4,2x,'lue sur le fichier dema
     *rrage est differente de la valeur parametree llm =',i4//)
   4  FORMAT(//10x,'la valeur de dtrv =',i4,2x,'lue sur le fichier dema
     *rrage est differente de la valeur  dtinteg =',i4//)

      RETURN
      END
