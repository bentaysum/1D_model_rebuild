!-----------------------------------------------------------------------
! INCLUDE datafile.h

!  Address of the directory containing tables of data needed by the GCM    
      COMMON/datadirectory/datafile
      character (len=150) :: datafile
! NB: default value for 'datafile' is set in inifis.F
!-----------------------------------------------------------------------
