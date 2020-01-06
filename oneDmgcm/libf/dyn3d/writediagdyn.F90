subroutine writediagdyn(name,title,units,dimpx,px)

! Write variable 'name' to NetCDF file 'diagdyn.nc'.
! The variable may be 3D (iip1,jjp1,llm) dynamical field,
! a 2D (iip1,jjp1) surface field, or a simple scalar (0D variable).
!
! Calls to 'writediagdyn' can originate from anywhere in the program;
! An initialisation of variable 'name' is done if it is the first time
! that this routine is called with given 'name'; otherwise data is appended
! (yielding the sought time series of the variable)
!
! NB: the rate a which outputs are made can be changed (see parameter isample)
!
implicit none

#include"dimensions.h"
#include"paramet.h"
#include"control.h"
#include"netcdf.inc"

! Arguments:
character(len=*),intent(in) :: name ! 'name' of the variable
character(len=*),intent(in) :: title ! 'long_name' attribute of the variable
character(len=*),intent(in) :: units ! 'units' attribute of the variable
integer,intent(in) :: dimpx ! dimension of the variable (3,2 or 0)
real,dimension(iip1,jjp1,llm),intent(in) :: px ! variable
! Note: px might actually be dimensio(iip1,jjp1,1); no problem

! Local variables:
!real,dimension(iip1,jjp1,llm) :: data3 ! to store 3D data
! Note iip1,jjp1 known from paramet.h; nsoilmx known from dimphys.h
!real,dimension(iip1,jjp1) :: data2 ! to store 2D data
!real :: data0 ! to store 0D data
integer :: i,j,l ! for loops
integer :: ig0

real,save :: date ! time counter (in elapsed days)
! sample rate at which data is to be written to output
!integer,parameter :: isample=1
integer,save :: isample ! initialized during Initialization step
integer,save :: ntime=0 ! counter to internally store time steps
character(len=20),save :: firstname="1234567890"
integer,save :: zitau=0

character(len=30) :: filename="diagdyn.nc"

! NetCDF stuff:
integer :: nid ! NetCDF output file ID
integer :: varid ! NetCDF ID of a variable
integer :: ierr ! NetCDF routines return code
integer,dimension(4) :: id ! NetCDF IDs of the dimensions of the variable
integer,dimension(4) :: edges,corners

! 1. Initialization step
if (firstname.eq."1234567890") then
  ! Initialize isample
  !isample=1
  isample=iphysiq
  
  ! Store 'name' as 'firstname'
  firstname=name
  ! From now on, if 'name'.eq.'firstname', then it is a new time cycle
  
  ! Create output NetCDF file
  ierr=NF_CREATE(filename,IOR(NF_CLOBBER,NF_64BIT_OFFSET),nid)
  if (ierr.ne.NF_NOERR) then
    write(*,*)'writediagdyn: Error, failed creating file '//trim(filename)
    stop
  endif
  
  ! Define dimensions and axis attributes
  call iniwritediagdyn(nid)
  
  ! set zitau to -1 to be compatible with zitau incrementation step below
  zitau=-1
  
else
  ! If not an initialization call, simply open the NetCDF file
  ierr=NF_OPEN(filename,NF_WRITE,nid)
  if (ierr.ne.NF_NOERR) then
    write(*,*)'writediagdyn: Error, failed opening file '//trim(filename)
    stop
  endif
endif ! of if (firstname.eq."1234567890")

! 2. Increment local time counter, if necessary
if (name.eq.firstname) then
  ! if we run across 'firstname', then it is a new dynamical time step
  zitau=zitau+1
endif

! 3. Write data, if the time index matches the sample rate
if (mod(zitau+1,isample).eq.0) then

! 3.1 If first call at this date, update 'time' variable
  if (name.eq.firstname) then
    ntime=ntime+1
    date=float(zitau+1)/float(day_step)
    ! Note: day_step is known from control.h
    
    ! Get NetCDF ID for "time"
    ierr=NF_INQ_VARID(nid,"time",varid)
    ! Add the current value of date to the "time" array
#ifdef NC_DOUBLE
    ierr=NF_PUT_VARA_DOUBLE(nid,varid,ntime,1,date)
#else
    ierr=NF_PUT_VARA_REAL(nid,varid,ntime,1,date)
#endif
    if (ierr.ne.NF_NOERR) then
      write(*,*)"writediagdyn: Failed writing date to time variable"
      stop 
    endif
  endif ! of if (name.eq.firstname)

! 3.2 Write the variable to the NetCDF file
if (dimpx.eq.3) then ! 3.2.1. Case of a 3D variable
  ! Write (append) the variable to the NetCDF file
  ! 3.2.1.A. Get the ID of the variable
  ierr=NF_INQ_VARID(nid,name,varid)
  if (ierr.ne.NF_NOERR) then
    ! If we failed geting the variable's ID, we assume it is because
    ! the variable doesn't exist yet and must be created.
    ! Start by obtaining corresponding dimensions IDs
    ierr=NF_INQ_DIMID(nid,"longitude",id(1))
    ierr=NF_INQ_DIMID(nid,"latitude",id(2))
    ierr=NF_INQ_DIMID(nid,"altitude",id(3))
    ierr=NF_INQ_DIMID(nid,"time",id(4))
    ! Tell the world about it
    write(*,*) "====================="
    write(*,*) "writediagdyn: creating variable "//trim(name)
    call def_var_diagdyn(nid,name,title,units,4,id,varid,ierr)
  endif ! of if (ierr.ne.NF_NOERR)
  
  ! 3.2.1.B. Prepare things to be able to write/append the variable
  corners(1)=1
  corners(2)=1
  corners(3)=1
  corners(4)=ntime
  
  edges(1)=iip1
  edges(2)=jjp1
  edges(3)=llm
  edges(4)=1
  
  ! 3.2.1.C. Write the slab of data
#ifdef NC_DOUBLE
  ierr=NF_PUT_VARA_DOUBLE(nid,varid,corners,edges,px)
#else
  ierr=NF_PUT_VARA_REAL(nid,varid,corners,edges,px)
#endif
  if (ierr.ne.NF_NOERR) then
    write(*,*) "writediagdyn: Error: Failed writing "//trim(name)//&
               " to file "//trim(filename)//" at time",date
  endif

elseif (dimpx.eq.2) then ! 3.2.2 Case of a 2D variable
  ! Write (append) the variable to the NetCDF file
  ! 3.2.2.A. Get the ID of the variable
  ierr=NF_INQ_VARID(nid,name,varid)
  if (ierr.ne.NF_NOERR) then
    ! If we failed geting the variable's ID, we assume it is because
    ! the variable doesn't exist yet and must be created.
    ! Start by obtaining corresponding dimensions IDs
    ierr=NF_INQ_DIMID(nid,"longitude",id(1))
    ierr=NF_INQ_DIMID(nid,"latitude",id(2))
    ierr=NF_INQ_DIMID(nid,"time",id(3))
    ! Tell the world about it
    write(*,*) "====================="
    write(*,*) "writediagdyn: creating variable "//trim(name)
    call def_var(nid,name,title,units,3,id,varid,ierr)
  endif ! of if (ierr.ne.NF_NOERR)

  ! 3.2.2.B. Prepare things to be able to write/append the variable
  corners(1)=1
  corners(2)=1
  corners(3)=ntime
  
  edges(1)=iip1
  edges(2)=jjp1
  edges(3)=1
  
  ! 3.2.2.C. Write the slab of data
#ifdef NC_DOUBLE
  ierr=NF_PUT_VARA_DOUBLE(nid,varid,corners,edges,px)
#else
  ierr=NF_PUT_VARA_REAL(nid,varid,corners,edges,px)
#endif
  if (ierr.ne.NF_NOERR) then
    write(*,*) "writediagdyn: Error: Failed writing "//trim(name)//&
               " to file "//trim(filename)//" at time",date
  endif

elseif (dimpx.eq.0) then ! 3.2.3. Case of a 0D variable
  ! Write (append) the variable to the NetCDF file
  ! 3.2.3.A. Get the ID of the variable
  ierr=NF_INQ_VARID(nid,name,varid)
  if (ierr.ne.NF_NOERR) then
    ! If we failed geting the variable's ID, we assume it is because
    ! the variable doesn't exist yet and must be created.
    ! Start by obtaining corresponding dimensions IDs
    ierr=NF_INQ_DIMID(nid,"time",id(1))
    ! Tell the world about it
    write(*,*) "====================="
    write(*,*) "writediagdyn: creating variable "//trim(name)
    call def_var(nid,name,title,units,1,id,varid,ierr)
  endif ! of if (ierr.ne.NF_NOERR)

  ! B.2. Prepare things to be able to write/append the variable
  corners(1)=ntime
  
  edges(1)=1

  ! B.3. Write the data
#ifdef NC_DOUBLE
  ierr=NF_PUT_VARA_DOUBLE(nid,varid,corners,edges,px)
#else
  ierr=NF_PUT_VARA_REAL(nid,varid,corners,edges,px)
#endif
  if (ierr.ne.NF_NOERR) then
    write(*,*) "writediagdyn: Error: Failed writing "//trim(name)//&
               " to file "//trim(filename)//" at time",date
  endif

endif ! of if (dimpx.eq.3) elseif (dimpx.eq.2) ...


endif ! of if (mod(zitau+1,isample).eq.0)

! 4. Close the NetCDF file
ierr=NF_CLOSE(nid)

end subroutine writediagdyn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine iniwritediagdyn(nid)

! Initialization routine for 'writediagdyn'. Here we create/define
! dimensions (longitude, latitude, altitude and time) and eventually
! other fixed (time-independent) parameters.

implicit none

#include"dimensions.h"
#include"paramet.h"
#include"comgeom.h"
#include"comconst.h"
#include"comvert.h"
#include"netcdf.inc"

! Arguments:
integer,intent(in) :: nid ! NetCDF output file ID

! Local variables:

! NetCDF stuff:
integer :: ierr ! NetCDF routines return code
integer :: idim_rlatu ! ID of the 'latitude' dimension
integer :: idim_rlonv ! ID of the 'longitude' dimension
integer :: idim_alt ! ID of the 'altitude' dimension
integer :: idim_time  ! ID of the 'time' dimension
integer :: varid ! to store NetCDF ID of a variable
integer,dimension(2) :: dimids ! to store IDs of dimensions of a variable
character(len=60) :: text ! to store some text


! 1. Define the dimensions
! Switch to NetCDF define mode
ierr=NF_REDEF(nid)

! Define the dimensions
ierr=NF_DEF_DIM(nid,"longitude",iip1,idim_rlonv)
! iip1 known from paramet.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not define longitude dimension"
endif
ierr=NF_DEF_DIM(nid,"latitude",jjp1,idim_rlatu)
! jjp1 known from paramet.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not define latitude dimension"
endif
ierr=NF_DEF_DIM(nid,"altitude",llm,idim_alt)
! llm known from dimensions.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not define depth dimension"
endif
ierr=NF_DEF_DIM(nid,"time",NF_UNLIMITED,idim_time)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not define time dimension"
endif

! Switch out of NetCDF define mode
ierr=NF_ENDDEF(nid)

! 2. Define (as variables) and write dimensions, as well as their attributes
! 2.1. Longitude
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,"longitude",NF_DOUBLE,1,idim_rlonv,varid)
#else
ierr=NF_DEF_VAR(nid,"longitude",NF_FLOAT,1,idim_rlonv,varid)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not define longitude variable"
endif

! Longitude attributes
text="East longitude"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="degrees_east"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

! Write longitude to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
#ifdef NC_DOUBLE
ierr=NF_PUT_VAR_DOUBLE(nid,varid,rlonv*(180./pi))
#else
ierr=NF_PUT_VAR_REAL(nid,varid,rlonv*(180./pi))
#endif
! Note: rlonv is known from comgeom.h and pi from comconst.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not write longitude variable"
endif

! 2.2. Latitude
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,"latitude",NF_DOUBLE,1,idim_rlatu,varid)
#else
ierr=NF_DEF_VAR(nid,"latitude",NF_FLOAT,1,idim_rlatu,varid)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not define latitude variable"
endif

! Latitude attributes
text="North latitude"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="degrees_north"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

! Write latitude to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
#ifdef NC_DOUBLE
ierr=NF_PUT_VAR_DOUBLE(nid,varid,rlatu*(180./pi))
#else
ierr=NF_PUT_VAR_REAL(nid,varid,rlatu*(180./pi))
#endif
! Note: rlatu is known from comgeom.h and pi from comconst.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not write latitude variable"
endif

! 2.3. Altitude
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,"altitude",NF_DOUBLE,1,idim_alt,varid)
#else
ierr=NF_DEF_VAR(nid,"altitude",NF_FLOAT,1,idim_alt,varid)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not define altitude variable"
endif

! Depth attributes
text="Pseudo-altitude"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="km"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)
text="up"
ierr=NF_PUT_ATT_TEXT(nid,varid,"positive",len_trim(text),text)

! Write depth to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
#ifdef NC_DOUBLE
ierr=NF_PUT_VAR_DOUBLE(nid,varid,pseudoalt)
#else
ierr=NF_PUT_VAR_REAL(nid,varid,pseudoalt)
#endif
! Note pseudoalt is known from comvert.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not write altitude variable"
endif

! 2.4. Time
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,"time",NF_DOUBLE,1,idim_time,varid)
#else
ierr=NF_DEF_VAR(nid,"time",NF_FLOAT,1,idim_time,varid)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not define time variable"
endif

! time attributes
text="Time"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="days since 0000-01-01 00:00:00"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Note no need to write time variable here; it is done in writediagsoil.

! 3. Other variables to be included

! 3.1 mesh area surrounding each horizontal point
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
dimids(1)=idim_rlonv ! ID of the 'longitude' dimension
dimids(2)=idim_rlatu ! ID of the 'latitude' dimension
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,"area",NF_DOUBLE,2,dimids,varid)
#else
ierr=NF_DEF_VAR(nid,"area",NF_FLOAT,2,dimids,varid)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not define area variable"
endif

! Area attributes
text="Mesh area"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="m2"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

! Write area to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
#ifdef NC_DOUBLE
ierr=NF_PUT_VAR_DOUBLE(nid,varid,aire)
#else
ierr=NF_PUT_VAR_REAL(nid,varid,aire)
#endif
! Note: aire is known from comgeom.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritediagdyn: Error, could not write area variable"
endif

end subroutine iniwritediagdyn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine def_var_diagdyn(nid,name,title,units,nbdim,dimids,nvarid,ierr)

! This subroutine defines variable 'name' in a (pre-existing and opened)
! NetCDF file (known from its NetCDF ID 'nid').
! The number of dimensions 'nbdim' of the variable, as well as the IDs of
! corresponding dimensions must be set (in array 'dimids').
! Upon successfull definition of the variable, 'nvarid' contains the
! NetCDF ID of the variable.
! The variables' attributes 'title' (Note that 'long_name' would be more
! appropriate) and 'units' are also set. 
! Modifs: Aug2010 Ehouarn: enforce outputs to be real*4

implicit none

#include "netcdf.inc"

integer,intent(in) :: nid ! NetCDF file ID
character(len=*),intent(in) :: name ! the variable's name
character(len=*),intent(in) :: title ! 'title' attribute of variable
character(len=*),intent(in) :: units ! 'units' attribute of variable
integer,intent(in) :: nbdim ! number of dimensions of the variable
integer,dimension(nbdim),intent(in) :: dimids ! NetCDF IDs of the dimensions
                                              ! the variable is defined along
integer,intent(out) :: nvarid ! NetCDF ID of the variable
integer,intent(out) :: ierr ! returned NetCDF staus code

! 1. Switch to NetCDF define mode 
ierr=NF_REDEF(nid)

! 2. Define the variable
#ifdef NC_DOUBLE
ierr = NF_DEF_VAR (nid,adjustl(name),NF_DOUBLE,nbdim,dimids,nvarid)
#else
ierr = NF_DEF_VAR (nid,adjustl(name),NF_FLOAT,nbdim,dimids,nvarid)
#endif
if(ierr/=NF_NOERR) then
   write(*,*) "def_var_diagdyn: Failed defining variable "//trim(name)
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

! 3. Write attributes
ierr=NF_PUT_ATT_TEXT(nid,nvarid,"title",&
                     len_trim(adjustl(title)),adjustl(title))
if(ierr/=NF_NOERR) then
   write(*,*) "def_var_diagdyn: Failed writing title attribute for "//trim(name)
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

ierr=NF_PUT_ATT_TEXT(nid,nvarid,"units",&
                     len_trim(adjustl(units)),adjustl(units))
if(ierr/=NF_NOERR) then
   write(*,*) "def_var_diagdyn: Failed writing units attribute for "//trim(name)
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

! 4. Switch out of NetCDF define mode
ierr = NF_ENDDEF(nid)

end subroutine def_var_diagdyn
