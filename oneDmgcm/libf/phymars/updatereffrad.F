      SUBROUTINE updatereffrad(ngrid,nlayer,
     &                rdust,rice,nuice,
     &                reffrad,nueffrad,
     &                pq,tauscaling,tau,pplay)
       USE updaterad
       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Subroutine designed to update the aerosol size distribution used by
c     the radiative transfer scheme. This size distribution is assumed
c     to be a log-normal distribution, with effective radius "reffrad" and
c     variance "nueffrad".
c   At firstcall, "rice" and "nuice" are not known, because
c     the H2O ice microphysical scheme is called after the radiative
c     transfer in physiq.F. That's why we assess the size of the 
c     water-ice particles at firstcall (see part 1.2 below).
c
c   author:   
c   ------
c   J.-B. Madeleine (2009-2010)
c
c=======================================================================
c
c    Declarations :
c    -------------
c
#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "dimradmars.h"
#include "tracer.h"
#include "aerkind.h"
#include "yomaer.h"

c-----------------------------------------------------------------------
c     Inputs:
c     ------

      INTEGER ngrid,nlayer
c     Ice geometric mean radius (m)
      REAL :: rice(ngridmx,nlayermx)
c     Estimated effective variance of the size distribution (n.u.)
      REAL :: nuice(ngridmx,nlayermx)
c     Tracer mass mixing ratio (kg/kg)
      REAL pq(ngrid,nlayer,nqmx)
      REAL rdust(ngridmx,nlayermx) ! Dust geometric mean radius (m)
      
      REAL pplay(ngrid,nlayer) ! altitude at the middle of the layers
      REAL tau(ngrid,naerkind)


c     Outputs:
c     -------

c     Aerosol effective radius used for radiative transfer (meter)
      REAL :: reffrad(ngridmx,nlayermx,naerkind)
c     Aerosol effective variance used for radiative transfer (n.u.)
      REAL :: nueffrad(ngridmx,nlayermx,naerkind)

c     Local variables:
c     ---------------

      INTEGER :: ig,l          ! 3D grid indices
      INTEGER :: iaer          ! Aerosol index

c     Number of cloud condensation nuclei near the surface
c     (only used at firstcall). This value is taken from 
c     Montmessin et al. 2004 JGR 109 E10004 p5 (2E6 part m-3), and
c     converted to part kg-1 using a typical atmospheric density.

      REAL, PARAMETER :: ccn0 = 1.3E8
      
c     For microphysics only:      
      REAL Mo,No                       ! Mass and number of ccn
      REAL rhocloud(ngridmx,nlayermx)  ! Cloud density (kg.m-3)
      REAL tauscaling(ngridmx)         ! Convertion factor for qccn and Nccn

      LOGICAL firstcall
      DATA firstcall/.true./
      SAVE firstcall

      REAL CBRT
      EXTERNAL CBRT

      REAL,SAVE :: nueffdust(ngridmx,nlayermx) ! Dust effective variance

c     Local saved variables:
c     ---------------------


c==================================================================
c 1. Update radius from fields from dynamics or initial state
c==================================================================

c       1.1 Dust particles
c       ------------------
        IF (doubleq.AND.active) THEN
          DO l=1,nlayer
            DO ig=1, ngrid
              call updaterdust(pq(ig,l,igcm_dust_mass),
     &                         pq(ig,l,igcm_dust_number),rdust(ig,l))
              nueffdust(ig,l) = exp(varian**2.)-1.
             ENDDO
           ENDDO
        ELSE
          DO l=1,nlayer
            DO ig=1, ngrid
              rdust(ig,l) = 0.8E-6
              nueffdust(ig,l) = 0.3
            ENDDO
          ENDDO
        ENDIF
        
c       1.2 Water-ice particles
c       -----------------------

        IF (water.AND.activice) THEN 
         IF (microphys) THEN
         
          IF (firstcall) THEN
            DO l=1,nlayer
              DO ig=1,ngrid
                call updaterice_micro(pq(ig,l,igcm_h2o_ice),
     &                                pq(ig,l,igcm_ccn_mass),
     &                                pq(ig,l,igcm_ccn_number),
     &                                1.e-3,rice(ig,l),
     &                                rhocloud(ig,l))
                nuice(ig,l) = nuice_ref
              ENDDO
            ENDDO
          firstcall = .false.
c    At firstcall, the true number and true mass of cloud condensation nuclei are not known.
c    Indeed it is scaled on the prescribed dust opacity via a 'tauscaling' coefficient
c    computed after radiative transfer.
          ELSE
            DO l=1,nlayer
              DO ig=1,ngrid
                call updaterice_micro(pq(ig,l,igcm_h2o_ice),
     &                                pq(ig,l,igcm_ccn_mass),
     &                                pq(ig,l,igcm_ccn_number),
     &                                tauscaling(ig),rice(ig,l),
     &                                rhocloud(ig,l))
                nuice(ig,l) = nuice_ref
               ENDDO
             ENDDO
          ENDIF ! of if firstcall
          
        ELSE ! if not microphys
         
          DO l=1,nlayer
            DO ig=1,ngrid    
              call updaterice_typ(pq(ig,l,igcm_h2o_ice),
     &                          tau(ig,1),pplay(ig,l),rice(ig,l)) 
              nuice(ig,l) = nuice_ref
            ENDDO
          ENDDO
 
        ENDIF ! of if microphys
       ENDIF ! of if (water.AND.activice)

c==================================================================
c 2. Radius used in the radiative transfer code (reffrad)
c==================================================================

      DO iaer = 1, naerkind ! Loop on aerosol kind
        aerkind: SELECT CASE (name_iaer(iaer))
c==================================================================
        CASE("dust_conrath") aerkind         ! Typical dust profile
c==================================================================
          DO l=1,nlayer
            DO ig=1,ngrid
              reffrad(ig,l,iaer) = rdust(ig,l) * 
     &          (1.e0 + nueffdust(ig,l))**2.5
              nueffrad(ig,l,iaer) = nueffdust(ig,l)
            ENDDO
          ENDDO
c==================================================================
        CASE("dust_doubleq") aerkind! Two-moment scheme for dust
c==================================================================
          DO l=1,nlayer
            DO ig=1,ngrid
              reffrad(ig,l,iaer) = rdust(ig,l) * ref_r0
              nueffrad(ig,l,iaer) = nueffdust(ig,l)
            ENDDO
          ENDDO
c==================================================================
        CASE("dust_submicron") aerkind   ! Small dust population
c==================================================================
          DO l=1,nlayer
            DO ig=1,ngrid
              reffrad(ig,l,iaer)=radius(igcm_dust_submicron)
              nueffrad(ig,l,iaer)=0.03
            ENDDO
          ENDDO     
c==================================================================
        CASE("h2o_ice") aerkind             ! Water ice crystals
c==================================================================
          DO l=1,nlayer
            DO ig=1,ngrid
c             About reffice, do not confuse the mass mean radius
c             (rayon moyen massique) and the number median radius
c             (or geometric mean radius, rayon moyen géométrique).
c             rice is a mass mean radius, whereas rdust
c             is a geometric mean radius:
c             number median rad = mass mean rad x exp(-1.5 sigma0^2)
c             (Montmessin et al. 2004 paragraph 30). Therefore:
              reffrad(ig,l,iaer)=rice(ig,l)*(1.+nuice_ref)
              nueffrad(ig,l,iaer)=nuice_ref
            ENDDO
          ENDDO
c==================================================================
        END SELECT aerkind
      ENDDO ! iaer (loop on aerosol kind)

      RETURN
      END 
