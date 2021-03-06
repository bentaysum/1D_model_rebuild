      SUBROUTINE disvert 

c    Auteur :  F. Forget Y. Wanherdrick, P. Levan
c    Nouvelle version 100% Mars !!

      IMPLICIT NONE

#include "dimensions.h"
#include "paramet.h"
#include "comvert.h"
#include "comconst.h"
#include "logic.h"
c
c=======================================================================
c    Discretisation verticale en coordonn�e hybride OU sigma
c
c=======================================================================
c
c   declarations:
c   -------------
c
c
      INTEGER l,ll
      REAL snorm
      REAL alpha,beta,gama,delta,deltaz,h,quoi,quand
      REAL zsig(llm)
      INTEGER np,ierr
      integer :: ierr1,ierr2,ierr3,ierr4
      REAL x

      REAL SSUM
      EXTERNAL SSUM
      real newsig 
      REAL dz0,dz1,nhaut,sig1,esig,csig,zz
      real tt,rr,gg, prevz
      real s(llm),dsig(llm) 

      integer iz 
      real z, ps,p

c
c-----------------------------------------------------------------------
c
      pi=2.*ASIN(1.)

! Ouverture possible de fichiers typiquement E.T.

         open(99,file="esasig.def",status='old',form='formatted',
     s   iostat=ierr2)
         if(ierr2.ne.0) then
              close(99)
              open(99,file="z2sig.def",status='old',form='formatted',
     s        iostat=ierr4)
         endif

c-----------------------------------------------------------------------
c   cas 1 on lit les options dans esasig.def:
c   ----------------------------------------

      IF(ierr2.eq.0) then

c        Lecture de esasig.def :
c        Systeme peu souple, mais qui respecte en theorie
c        La conservation de l'energie (conversion Energie potentielle
c        <-> energie cinetique, d'apres la note de Frederic Hourdin...

         PRINT*,'*****************************'
         PRINT*,'WARNING Lecture de esasig.def'
         PRINT*,'*****************************'
         READ(99,*) h
         READ(99,*) dz0
         READ(99,*) dz1
         READ(99,*) nhaut
         CLOSE(99)

         dz0=dz0/h
         dz1=dz1/h

         sig1=(1.-dz1)/tanh(.5*(llm-1)/nhaut)

         esig=1.

         do l=1,20
            esig=-log((1./sig1-1.)*exp(-dz0)/esig)/(llm-1.)
         enddo
         csig=(1./sig1-1.)/(exp(esig)-1.)

         DO L = 2, llm
            zz=csig*(exp(esig*(l-1.))-1.)
            sig(l) =1./(1.+zz)
     &      * tanh(.5*(llm+1-l)/nhaut)
         ENDDO
         sig(1)=1.
         sig(llm+1)=0.
         quoi      = 1. + 2.* kappa
         s( llm )  = 1.
         s(llm-1) = quoi
         IF( llm.gt.2 )  THEN
            DO  ll = 2, llm-1
               l         = llm+1 - ll
               quand     = sig(l+1)/ sig(l)
               s(l-1)    = quoi * (1.-quand) * s(l)  + quand * s(l+1)
            ENDDO
         END IF
c
         snorm=(1.-.5*sig(2)+kappa*(1.-sig(2)))*s(1)+.5*sig(2)*s(2)
         DO l = 1, llm
            s(l)    = s(l)/ snorm
         ENDDO

c-----------------------------------------------------------------------
c   cas 2 on lit les options dans z2sig.def:
c   ----------------------------------------

      ELSE IF(ierr4.eq.0) then
         PRINT*,'****************************'
         PRINT*,'Lecture de z2sig.def'
         PRINT*,'****************************'

         READ(99,*) h
         do l=1,llm
            read(99,*) zsig(l)
         end do
         CLOSE(99)

         sig(1) =1
         do l=2,llm
           sig(l) = 0.5 * ( exp(-zsig(l)/h) + exp(-zsig(l-1)/h) )
         end do
         sig(llm+1) =0

c-----------------------------------------------------------------------
      ELSE
         write(*,*) 'didn t you forget something ??? '
         write(*,*) 'We need the file  z2sig.def ! (OR esasig.def) '
         stop
      ENDIF
c-----------------------------------------------------------------------

      DO l=1,llm
        nivsigs(l) = REAL(l)
      ENDDO

      DO l=1,llmp1
        nivsig(l)= REAL(l)
      ENDDO

 
c-----------------------------------------------------------------------
c    ....  Calculs  de ap(l) et de bp(l)  ....
c    .........................................
c
c   .....  pa et preff sont lus  sur les fichiers start par dynetat0 .....
c-----------------------------------------------------------------------
c

      if (hybrid) then
         write(*,*) "*******************************"
         write(*,*) "Systeme en coordonnees hybrides"
         write(*,*) 
c        Coordonnees hybrides avec mod
         DO l = 1, llm

         call sig_hybrid(sig(l),pa,preff,newsig)
            bp(l) = EXP( 1. - 1./(newsig**2)  )
            ap(l) = pa * (newsig - bp(l) )
         enddo
         bp(llmp1) = 0.
         ap(llmp1) = 0.
      else
         write(*,*) "****************************"
         write(*,*) "Systeme en coordonnees sigma"
         write(*,*) 
c        Pour ne pas passer en coordonnees hybrides
         DO l = 1, llm
            ap(l) = 0.
            bp(l) = sig(l)
         ENDDO
         ap(llmp1) = 0.
      endif

      bp(llmp1) =   0.

      PRINT *,' BP '
      PRINT *,  bp
      PRINT *,' AP '
      PRINT *,  ap

c     Calcul au milieu des couches :
c     WARNING : le choix de placer le milieu des couches au niveau de
c     pression interm�diaire est arbitraire et pourrait etre modifi�.
c     Le calcul du niveau pour la derniere couche 
c     (on met la meme distance (en log pression)  entre P(llm)
c     et P(llm -1) qu'entre P(llm-1) et P(llm-2) ) est
c     Specifique.  Ce choix est sp�cifi� ici ET dans exner_hyb.F

      DO l = 1, llm-1
       aps(l) =  0.5 *( ap(l) +ap(l+1)) 
       bps(l) =  0.5 *( bp(l) +bp(l+1)) 
      ENDDO
     
      if (hybrid) then
         aps(llm) = aps(llm-1)**2 / aps(llm-2) 
         bps(llm) = 0.5*(bp(llm) + bp(llm+1))
      else
         bps(llm) = bps(llm-1)**2 / bps(llm-2) 
         aps(llm) = 0.
      end if

      PRINT *,' BPs '
      PRINT *,  bps
      PRINT *,' APs'
      PRINT *,  aps

      DO l = 1, llm
       presnivs(l) = aps(l)+bps(l)*preff
       pseudoalt(l) = -10.*log(presnivs(l)/preff)
      ENDDO

      PRINT *,' PRESNIVS' 
      PRINT *,presnivs 
      PRINT *,'Pseudo altitude des Presnivs : ' 
      PRINT *,pseudoalt

c     --------------------------------------------------
c     This can be used to plot the vertical discretization
c     (> xmgrace -nxy testhybrid.tab       (z = H*log(p(l)/pref))
c     --------------------------------------------------
c     open (53,file='testhybrid.tab')
c     do iz=0,60
c       z = -10 + min(iz,60-iz)
c       ps = preff*exp(-z/10)
c       do l=1,llm
c          zsig(l)= -10.*log((aps(l) + bps(l)*ps)/preff)
c       end do
c       write(53,*)iz, (zsig(l),l=1,llm,1)
c      end do
c      close(53)
c     --------------------------------------------------


      RETURN
      END

c ************************************************************
      subroutine sig_hybrid(sig,pa,preff,newsig)
c     ----------------------------------------------
c     Subroutine utilisee pour calculer des valeurs de sigma modifie
c     pour conserver les coordonnees verticales decrites dans
c     esasig.def/z2sig.def lors du passage en coordonnees hybrides
c     F. Forget 2002
c     Connaissant sig (niveaux "sigma" ou on veut mettre les couches)
c     L'objectif est de calculer newsig telle que
c       (1 -pa/preff)*exp(1-1./newsig**2)+(pa/preff)*newsig = sig
c     Cela ne se r�soud pas analytiquement: 
c     => on r�soud par iterration bourrine 
c     ----------------------------------------------
c     Information  : where exp(1-1./x**2) become << x
c           x      exp(1-1./x**2) /x
c           1           1
c           0.68       0.5
c           0.5        1.E-1
c           0.391      1.E-2
c           0.333      1.E-3
c           0.295      1.E-4
c           0.269      1.E-5
c           0.248      1.E-6
c        => on peut utiliser newsig = sig*preff/pa si sig*preff/pa < 0.25


      implicit none
      real x1, x2, sig,pa,preff, newsig, F
      integer j

      newsig = sig
      x1=0
      x2=1
	  if (sig.ge.1) then
	       newsig= sig
      else if (sig*preff/pa.ge.0.25) then
        DO J=1,9999  ! nombre d''iteration max
          F=((1 -pa/preff)*exp(1-1./newsig**2)+(pa/preff)*newsig)/sig
c         write(0,*) J, ' newsig =', newsig, ' F= ', F
          if (F.gt.1) then
              X2 = newsig
              newsig=(X1+newsig)*0.5
          else
              X1 = newsig
              newsig=(X2+newsig)*0.5
          end if
c         Test : on arete lorsque on approxime sig � moins de 0.01 m pr�s 
c         (en pseudo altiude) :
          IF(abs(10.*log(F)).LT.1.E-5) goto 999
        END DO
       else   !    if (sig*preff/pa.le.0.25) then
	        newsig= sig*preff/pa
       end if
 999   continue
       Return
      END
