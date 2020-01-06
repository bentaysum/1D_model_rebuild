      subroutine thermcell_dqup(ngrid,nlayer,ptimestep,fm,entr,detr,  &
     &    masse0,q_therm,dq_therm,ndt,zlmax)
      implicit none

!=======================================================================
!
!   Calcul du transport verticale dans la couche limite en presence
!   de "thermiques" explicitement representes
!   calcul du dq/dt une fois qu'on connait les ascendances
!   Version modifiee pour prendre les downdrafts a la place de la 
!   subsidence compensatoire
!
!   Version with sub-timestep for Martian thin layers
!
!=======================================================================

#include "dimensions.h"
#include "dimphys.h"

! ============================ INPUTS ============================

      INTEGER, INTENT(IN) :: ngrid,nlayer
      REAL, INTENT(IN) :: ptimestep
      REAL, INTENT(IN) :: fm(ngridmx,nlayermx+1)
      REAL, INTENT(IN) :: entr(ngridmx,nlayermx)
      REAL, INTENT(IN) :: detr(ngridmx,nlayermx)
      REAL, INTENT(IN) :: q_therm(ngridmx,nlayermx)
      REAL, INTENT(IN) :: masse0(ngridmx,nlayermx)
      INTEGER, INTENT(IN) :: ndt
      INTEGER, INTENT(IN) :: zlmax

! ============================ OUTPUTS ===========================

      REAL, INTENT(OUT) :: dq_therm(ngridmx,nlayermx)  ! dq/dt -> derivative

! ============================ LOCAL =============================

      REAL q(ngridmx,nlayermx)
      REAL qa(ngridmx,nlayermx)
      INTEGER ig,k,i
      REAL invflux0(ngridmx,nlayermx)
      REAL ztimestep

! =========== Init ==============================================

      qa(:,:)=q_therm(:,:)
      q(:,:)=q_therm(:,:)

! ====== Computing q ============================================

      dq_therm(:,:)=0.
      ztimestep=ptimestep/real(ndt)
      invflux0(:,:)=ztimestep/masse0(:,:)      

      do i=1,ndt

      do ig=1,ngrid
         qa(ig,1)=q(ig,1)
      enddo

      do k=2,zlmax
         do ig=1,ngridmx
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.  &
     &         1.e-5*masse0(ig,k)) then
         qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+entr(ig,k)*q(ig,k))  &
     &     /(fm(ig,k+1)+detr(ig,k))
            else
               qa(ig,k)=q(ig,k)
            endif
         enddo
      enddo

      do k=1,zlmax
           q(:,k)=q(:,k)+         &
     & (detr(:,k)*qa(:,k)-entr(:,k)*q(:,k) &
     &    -fm(:,k)*q(:,k)+fm(:,k+1)*q(:,k+1))  &
     &               *invflux0(:,k)
      enddo

      enddo !of do i=1,ndt

! ====== Derivative ==============================================


         do k=1,zlmax
          dq_therm(:,k)=(q(:,k)-q_therm(:,k))/ptimestep
         enddo

! ==============

      return
      end

