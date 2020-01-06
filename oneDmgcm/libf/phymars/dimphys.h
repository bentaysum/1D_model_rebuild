!-----------------------------------------------------------------------
!   INCLUDE 'dimphys.h'

! ngridmx : number of horizontal grid points
! note: the -1/jjm term will be 0; unless jj=1
      integer, parameter :: ngridmx = (2+(jjm-1)*iim - 1/jjm)   
! nlayermx : number of atmospheric layers
      integer, parameter :: nlayermx = llm 
! nsoilmx : number of subterranean layers
!EM: old soil routine:      integer, parameter :: nsoilmx = 10
      integer, parameter :: nsoilmx = 18 
!-----------------------------------------------------------------------
