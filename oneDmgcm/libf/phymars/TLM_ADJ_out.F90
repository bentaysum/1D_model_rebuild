SUBROUTINE TLM_ADJ_out(idt,ptime)

USE ioipsl_getincom 

USE TLMvars
USE netcdf

IMPLICIT NONE

! Commons
#include "callkeys.h"
#include "dimensions.h"
#include "dimphys.h"
#include "paramet.h"
#include "control.h"
#include "comvert.h"
#include "comgeom.h"
#include "description.h"
#include "surfdat.h"
#include "chimiedata.h"
#include "tracer.h"
#include "conc.h"

! Aiming to now produce a binary datafile holding the TLM (or adjoint) matrices 
! to save computational burden experienced with netcdf routines.

! ================================================================================== !
! ============================= Variable Declaration =============================== !
! ================================================================================== !
! INPUT 
! -----
integer :: idt 
real :: ptime ! local time 00:00 - 23:59 hrs
! LOCAL 
! -----
logical, save :: firstcall = .True.
integer,save :: t_idx = 1
integer, save :: ndt = 0
character(len=20) filename
character(len=*), parameter :: directory = "/exports/csce/datastore/geos/users/s1215319/paper2/tlmfiles/"
integer, parameter :: nx = nlayermx*nqmx
integer, parameter :: ny = nlayermx*nqmx
integer, save :: nt 

real :: vmr(nlayermx)
integer :: iq

! STORAGE 
! -------
real*8, allocatable, save :: tangent_matrix(:,:,:)
integer x, y, t ! Loop iterators
integer, parameter :: wl = 8
integer :: length
integer :: l
integer :: iostat
character(len=15) :: HEADER_FMT

! WRITING TO TERMINAL 
integer :: loc_max(2)
integer :: loc_min(2)
integer :: idx_j(2), idx_i(2) 

! STAGE 1: INITIALISATION
! 		   - Construct a 3D structure that can hold the TLM sets
! 			 prior to inserting into the binary data file.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IF ( firstcall ) THEN
	
	firstcall = .False.
	
	call getin("ndt",ndt)
	
	ndt = ndt*day_step
	
	nt = ndt - idt + 1
	ALLOCATE ( tangent_matrix(nx,ny,nt) )
! Save the initial model state (in VMR units)
	length = nlayermx*wl
	
	open (unit = 15, file = directory // "tracer_index.txt")
	
	! Save idt, so we can jump directly to relevant time index for control and perturbed files
	! in post-analysis.
	HEADER_FMT = "(A12,A3,I4)"
	write(15,"(A12,A3,F5.2)") "initial lt", " : ", ptime
	write(15,HEADER_FMT) "initial idt", " : ", idt
	write(15,HEADER_FMT) "N steps", " : ", nt
	write(15,HEADER_FMT) "nqmx", " : ", nqmx 
	write(15,HEADER_FMT) "nlayermx", " : ", nlayermx
	write(15,"(A16)") "################"
	
	do iq = 1, nqmx
		
		write(15,'(A12,I2)') ADJUSTL(NOMS(IQ)), iq
		
	enddo 
	
	CLOSE(15)

ENDIF
! write(*,*) "--------------------------------"
! write(*,*) "TLM MAX: ", MAXVAL(tangent_matrix)
! write(*,*) "TLM MIN: ", MINVAL(tangent_matrix)
! write(*,*) "--------------------------------"

! STAGE 2: ADD TO THE STORAGE STRUCTURE
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	tangent_matrix(:,:,t_idx) = TLM 


loc_max = MAXLOC( tangent_matrix(:,:,t_idx) )
loc_min = MINLOC( tangent_matrix(:,:,t_idx) ) 

idx_j(1) = CEILING( loc_max(1)/FLOAT(nlayermx) ) 
idx_j(2) = CEILING( loc_min(1)/FLOAT(nlayermx) ) 

idx_i(1) = CEILING( loc_max(2)/FLOAT(nlayermx) ) 
idx_i(2) = CEILING( loc_min(2)/FLOAT(nlayermx) ) 



write(*,*) "===================================="
write(*,"(E15.7,2A10)")MAXVAL(tangent_matrix(:,:,t_idx)), trim(noms(idx_j(1))), trim(noms(idx_i(1)))
write(*,"(E15.7,2A10)")MINVAL(tangent_matrix(:,:,t_idx)), trim(noms(idx_j(2))), trim(noms(idx_i(2)))
write(*,*) "===================================="

! STAGE 3: CONSTRUCT THE BINARY FILE AND SAVE THE STRUCTURE
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IF ( idt == ndt ) THEN 
				
		filename = "tlm.bin"
			
	length = nx*ny*wl
		
	open (unit=11, file= directory // filename, access='direct', recl=length, iostat=iostat)
	
	! DO k = 1, nt 
		! WRITE(11,rec=k) ( (tangent_matrix(x,y,k),x=1,nx ), y=1,ny)
	! ENDDO 
	do t = 1, nt		
		write (11, rec=t) ((tangent_matrix(x,y,t),x=1,nx), y=1,ny)
	end do
	
	CLOSE(11)
	
	write(*,*) SHAPE(tangent_matrix), SIZE(tangent_matrix)

	
ENDIF 

t_idx = t_idx + 1
 
 
END SUBROUTINE
