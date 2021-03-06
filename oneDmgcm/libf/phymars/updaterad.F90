module updaterad



! This module intents to group together all ice and dust radius computations done in the GCM,
! so that it stays coherent through the code and make (numerous) radius bugs easier to debug.
! All different thresholds values are defined below.

! Radius computation do not always occur on the whole grid (cf. improvedcloud).
! So, subroutines are designed for scalar values instead of tables


! T. Navarro, June 2012



! For instance, if R^3 is lower than r3icemin, then R is set to ricemin.
! So, ricemin should be the cubic root of r3icemin, but not necessarily ...

real, parameter :: r3icemin = 1.e-30   ! ie ricemin  = 0.0001 microns
real, parameter :: ricemin  = 1.e-10

real, parameter :: r3icemax = 125.e-12 ! ie ricemax  = 500 microns
real, parameter :: ricemax  = 500.e-6 

real, parameter :: qice_threshold  = 1.e-15 ! 1.e-10



real, parameter :: nccn_threshold  =  1.
real, parameter :: qccn_threshold  =  1.e-20

real, parameter :: r3ccnmin = 1.e-21    ! ie rccnmin = 0.1 microns
real, parameter :: rccnmin  = 0.1e-6

real, parameter :: r3ccnmax = 125.e-12  ! ie rccnmax  = 500 microns
real, parameter :: rccnmax  = 500.e-6  



real, parameter :: ndust_threshold  = 1.
real, parameter :: qdust_threshold  = 1.e-20

real, parameter :: r3dustmin = 1.e-21  ! ie rdustmin = 0.1 microns
real, parameter :: rdustmin  = 0.1e-6

real, parameter :: r3dustmax = 125.e-12 ! ie rdustmax  = 500 microns
real, parameter :: rdustmax  = 500.e-6  


real, parameter :: rdust0 = 0.8e-6


      

contains


!============================================================================
!============================================================================
!============================================================================
! Update ice radius if microphys == true 
subroutine updaterice_micro(qice,qccn,nccn,coeff,rice,rhocloud)
implicit none

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "tracer.h"

real, intent(in)  :: qice,qccn,nccn
real, intent(in)  :: coeff         ! this coeff is tauscaling if microphy = T (possibly ccn_factor^-1 otherwise)
real, intent(out) :: rice,rhocloud ! rhocloud is needed for sedimentation and is also a good diagnostic variable
real nccn_true,qccn_true ! radius of the ice crystal core to the power of 3

    
nccn_true = max(nccn * coeff, 1.e-30)
qccn_true = max(qccn * coeff, 1.e-30)


!! Nota: It is very dangerous to apply a threshold on qccn or nccn to force rice to be ricemin.
!! Indeed, one can obtain ricemin for small but non negligible qice values, and therefore hugely opaque clouds.
 
if (qice .ge. qice_threshold) then

  rhocloud = (qice*rho_ice + qccn_true*rho_dust) / (qice + qccn_true)
  rhocloud = min(max(rhocloud,rho_ice),rho_dust)

  rice = (qice + qccn_true) * 0.75 / pi / rhocloud / nccn_true
  
  if (rice .le. r3icemin) then
    rice = ricemin
  else if (rice .ge. r3icemax) then
    rice = ricemax
  else
    rice = rice**(1./3.) ! here rice is always positive
  endif

else

  rice = ricemin
  rhocloud = rho_dust

endif

end subroutine updaterice_micro
!============================================================================
!============================================================================
!============================================================================





!============================================================================
!============================================================================
!============================================================================
! Update ice radius from a typical profile if microphys == false
subroutine updaterice_typ(qice,tau,pzlay,rice)
implicit none

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "tracer.h"

real, intent(in)  :: qice
real, intent(in)  :: tau   ! tau for dust
real, intent(in)  :: pzlay ! altitude at the middle of the layers
real, intent(out) :: rice
real rccn,nccn ! radius  and number of ice crystals

!   Typical CCN profile following Montmessin et al. 2004
!   (N0=2e6 m-3 has been converted into N0=1.3e8 kg-1, otherwise the equation for rice is not homogeneous...)
     nccn  = 1.3e+8*max(tau,0.001)/0.1*exp(-pzlay/10000.)
!   The previously used profile was not correct:
!   Nccn=( epaisseur/masse ) * 2.e+6/0.1*max(tau,0.001)*exp(-pzlay/10000.)
     
if (nccn .le. 1) then

  rice = ricemin

else

! Typical dust radius profile:
  rccn = max(rdust0*exp(-pzlay/18000.),1.e-9)
  rice  = qice * 0.75 / pi / rho_ice / nccn + rccn*rccn*rccn
       
  if (rice .le. r3icemin) then
    rice = ricemin
  else if (rice .ge. r3icemax) then
    rice = ricemax
  else
    rice = rice**(1./3.) ! here rice is always positive
  endif

endif
 
end subroutine updaterice_typ
!============================================================================
!============================================================================
!============================================================================






!============================================================================
!============================================================================
!============================================================================
! This subroutine computes the geometric mean radius(or number median radius)
! For a lognormal distribution :
! geometric mean radius = mass mean radius x exp(-1.5 sigma0^2)
! To be used with doubleq == true. otherwise, rdust is constant !!!
subroutine updaterdust(qdust,ndust,rdust,tauscaling)
implicit none

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "tracer.h"

real, intent(in) :: qdust,ndust ! needed if doubleq
real, intent(in), optional :: tauscaling ! useful for realistic thresholds
real, intent(out) :: rdust
real coeff

if (present(tauscaling)) then
  coeff = tauscaling ! thresholds on realistic values
else
  coeff = 1. ! thresholds on virtual values
endif

if ((ndust .le. ndust_threshold/coeff) .or. (qdust .le. qdust_threshold/coeff)) then

   rdust = rdustmin

else

   rdust  = r3n_q * qdust / ndust

   if (rdust .le. r3dustmin) then
       rdust = rdustmin
   else if (rdust .ge. r3dustmax) then
       rdust = rdustmax
   else
      rdust = rdust**(1./3.)
   endif

endif
     

end subroutine updaterdust
!============================================================================
!============================================================================
!============================================================================






!============================================================================
!============================================================================
!============================================================================
! This subroutine computes the mass mean radius, 
! used for heterogenous nucleation on CCNs in microphysics.
! For a lognormal distribution :
! geometric mean radius = mass mean radius x exp(-1.5 sigma0^2)
subroutine updaterccn(qccn,nccn,rccn,tauscaling)
implicit none

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "tracer.h"

real, intent(in) :: qccn,nccn ! needed if doubleq
real, intent(in), optional :: tauscaling ! useful for realistic thresholds

real, intent(out) :: rccn

real coeff


if (present(tauscaling)) then
  coeff = tauscaling ! threshold on realistic values
else
  coeff = 1. ! threshold on virtual values
endif

if ((nccn .le. nccn_threshold/coeff) .or. (qccn .le. qccn_threshold/coeff)) then

   rccn = rccnmin

else

   rccn  =  qccn * 0.75 / pi / rho_dust / nccn

   if (rccn .le. r3ccnmin) then
       rccn = rccnmin
   else if (rccn .ge. r3ccnmax) then
       rccn = rccnmax
   else
       rccn = rccn**(1./3.)
   endif

endif
     
end subroutine updaterccn
!============================================================================
!============================================================================
!============================================================================



end module updaterad

