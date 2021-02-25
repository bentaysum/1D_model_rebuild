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
! character(len=*), parameter :: directory = "/exports/csce/datastore/geos/users/s1215319/paper2/tlmfiles/"
character(len=*), parameter :: directory = "/scratch/local/s1215319/"
! character(len=*), parameter :: directory = "/exports/csce/datastore/geos/users/" &
!                                 // "s1215319/paper3/version2/tlmfiles/o2sensitivity/standard/"

integer, parameter :: nx = nlayermx*nqmx
integer, parameter :: ny = nlayermx*nqmx
integer, save :: nt 

real :: vmr(nlayermx)
integer :: iq

! STORAGE 
! -------
integer x, y ! Loop iterators
integer,save :: t 
integer, parameter :: wl = 4
integer :: length
integer :: l
integer :: iostat
character(len=15) :: HEADER_FMT
character(len=3), SAVE :: day0_string
integer :: day0 

! WRITING TO TERMINAL 
integer :: loc_max(2)
integer :: loc_min(2)
integer :: idx_j(2), idx_i(2) 




! STAGE 1: INITIALISATION
! 		   - Construct a 3D structure that can hold the TLM sets
! 			 prior to inserting into the binary data file.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IF ( firstcall ) THEN
    
    t = 0 

    firstcall = .False.

     call getin("day0",day0)
     
     write(day0_string,"(I3)") day0 
     
    call getin("ndt",ndt)

    ndt = ndt*day_step

    nt = ndt - idt + 1

    ! Save the initial model state (in VMR units)
    length = nlayermx*wl

    ! open (unit = 15, file = directory // ADJUSTL(trim(day0_string)) // "_tracer_index.txt")
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

    ! Length of file to be created
    length = nx*ny*wl

    ! If files exist remove them
    ! --------------------------
    open(unit=11, iostat=iostat, file=directory // "tlm.bin", status='old')
    if (iostat == 0) close(11, status='delete')
    close(11)

    ! open(unit=111, iostat=iostat, file=directory // "steady_o2coefficients.bin", status='old')
    ! if (iostat == 0) close(111, status='delete')
    ! close(111)

    ! Adjoint/Tangent Linear Matrix 
    ! =============================
    ! Create the file in scratch space
    ! open (unit=11, file= directory // "tlm.bin", access='direct', recl=length, &
    !  STATUS = "REPLACE", iostat=iostat)
    open(11,file= directory // "tlm.bin",status="new",action="write",access='stream',form='unformatted')

    ! Coefficients for O2 Steady-State
    ! ================================
    ! open(111,file= directory // "steady_o2coefficients.bin",status="new",action="write",access='stream',form='unformatted')


    CLOSE(11)
    ! CLOSE(111)

ENDIF
! write(*,*) "--------------------------------"
! write(*,*) "TLM MAX: ", MAXVAL(tangent_matrix)
! write(*,*) "TLM MIN: ", MINVAL(tangent_matrix)
! write(*,*) "--------------------------------"

! STAGE 2: ADD TO THE STORAGE STRUCTURE
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! OPEN(11,file= directory // "tlm.bin" ,action='write',position='append')
open(11,file= directory // "tlm.bin", status="old",action="write",&
    access='stream',form='unformatted',position='append')
! open(111,file= directory // "steady_o2coefficients.bin", &
!     status="old",action="write",access='stream',form='unformatted',position='append')


loc_max = MAXLOC( TLM )
loc_min = MINLOC( TLM ) 

idx_j(1) = CEILING( loc_max(1)/FLOAT(nlayermx) ) 
idx_j(2) = CEILING( loc_min(1)/FLOAT(nlayermx) ) 

idx_i(1) = CEILING( loc_max(2)/FLOAT(nlayermx) ) 
idx_i(2) = CEILING( loc_min(2)/FLOAT(nlayermx) ) 


write(*,*) "===================================="
write(*,"(E15.7,2A10)") MAXVAL(TLM), trim(noms(idx_j(1))), trim(noms(idx_i(1)))
write(*,"(E15.7,2A10)") MINVAL(TLM), trim(noms(idx_j(2))), trim(noms(idx_i(2)))
write(*,*) "===================================="

t = t + 1
! write(*,"(E15.7,A)")  (t/nt)*1.e2 , "%"
WRITE(*,*) T, NT
! WRITE(11,*) ((TLM(x,y), x=1,nx), y=1, ny)
WRITE(11) TLM
CLOSE(11) 

! WRITE(111) o2_coefficient_array
! CLOSE(111)

RETURN 



 
 
END SUBROUTINE
