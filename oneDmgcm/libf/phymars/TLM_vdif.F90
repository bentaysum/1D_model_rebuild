subroutine TLM_vdif(iq, za, zb, zc, zd, ptimestep)

! 20/08/2018 - Ben Taysum
! Updated - 03/09/2018 - Ben Taysum
! Re-made - 25/09/2018 - Ben Taysum
!
! Matrix structures:
!
! perts = [nlayermx*nqmx]
! TLM = [nqmx*nlayermx,nqmx*nlayermx]
!
! Here we calculate the indices of the matrix M in the
! equation:
!			[zc]'' = M x [pq]''
!
! and the indices of M will be added to the relevant
! TLM indices.

USE TLMvars

IMPLICIT NONE

#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "conc.h"

! =================================================

! Input variables
INTEGER :: iq 
REAL :: za(nlayermx), zb(nlayermx), zc(nlayermx), zd(nlayermx)
REAL :: ptimestep

! Local variables
REAL :: z1(nlayermx) ! Z1 from vdifc
REAL :: M(nlayermx,nlayermx) 
REAL :: zcdiff(nlayermx,nlayermx)! The matrix used to compute zc''
INTEGER :: x, y ! Loop iterators
INTEGER :: i, n, a, b, c, d, e ! Loop iterators matching the LaTeX document
REAL    :: coeff ! Holds the summation values
! =================================================

!==========================
! Step One : Initialise M and the TLM identity matrix
!==========================

M(:,:) = 0.
zcdiff(:,:) = 0.
!===============================
! Step Two : Compute z1 values =
!===============================

!==== Top layer (nlayermx)
z1(nlayermx) = 1./(za(nlayermx) + zb(nlayermx))
!==== Mid layers (1 < n < nlayermx)
DO y = 2,nlayermx-1
	z1(y) = 1./(za(y)+zb(y) + zb(y+1)*(1.-zd(y+1)))
ENDDO
!==== Bottom layers (n=1)
z1(1) = 1./(za(1) + zb(1) + zb(2)*(1.-zd(2)))    

!======================================
! Step Three : ZCDIFF    calculations =
!======================================
DO n = 1,nlayermx
	DO i = 1, nlayermx
		IF (i .lt. n) THEN
			zcdiff(n,i) = 0.
		ELSE IF (i == n) THEN
			zcdiff(n,i) = za(i)*z1(i)
		ELSE IF (i .gt. n) THEN
			zcdiff(n,i) = za(i)*PRODUCT(z1(n:i))*PRODUCT(zb(n+1:i))
		ENDIF
	ENDDO
ENDDO

!======================================
! Step Four : M matrix   calculations =
!======================================
DO n = 1,nlayermx
	IF ( n == 1) THEN
		M(n,:) = zcdiff(n,:)		
	ELSE
		M(n,:) = zcdiff(n,:) + zd(n)*M(n-1,:)
	ENDIF
ENDDO

DO n = 1,nlayermx
	M(n,n) = M(n,n) - 1.
ENDDO

!============================================== 
! Step Four : Allocating the correct index of =
! 			  M to the correct index of TLM   =
!==============================================
TlM_trans( (iq-1)*nlayermx + 1 : iq*nlayermx,  (iq-1)*nlayermx + 1 : iq*nlayermx) = M/ptimestep
end subroutine
