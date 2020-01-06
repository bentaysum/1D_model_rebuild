!-----------------------------------------------------------------------
! INCLUDE 'control.h'
! For Fortran 77/Fortran 90 compliance always use line continuation
! symbols '&' in columns 73 and 6
!

      COMMON/control_i/ndynstep,day_step,                               &
     &              iperiod,iconser,idissip,iphysiq ,                   &
     &              anneeref
      COMMON/control_r/periodav,ecritphy,nday_r

      INTEGER ndynstep ! # of dynamical time steps to run (if negative or not specified in run.def, nday_r is used instead)
      INTEGER day_step ! # of dynamical time steps per day
      INTEGER iperiod  ! make a Matsuno step before avery iperiod-1 LF steps
      INTEGER iconser !
      INTEGER idissip ! apply dissipation every idissip dynamical step
      INTEGER iphysiq ! call physics every iphysiq dynamical steps
      INTEGER anneeref ! reference year # ! not used
      REAL periodav
      REAL ecritphy ! output data in "diagfi.nc" every ecritphy dynamical steps 
      real nday_r ! number of days to run (possibly including a fraction of day)

!-----------------------------------------------------------------------
