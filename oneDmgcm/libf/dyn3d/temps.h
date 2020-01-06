c-----------------------------------------------------------------------
c INCLUDE 'temps.h'

      COMMON/temps_i/day_ini,day_end,anne_ini,itaufin
      COMMON/temps_r/dt

      INTEGER  itaufin  ! total number of dynamical steps for the run
      INTEGER*4 day_ini ! initial day # of simulation sequence
      INTEGER*4 day_end ! final day # ; i.e. day # when this simulation ends
      INTEGER*4 anne_ini ! initial year # of simulation sequence ? Not used.
      REAL dt ! (dynamics) time step (changes if doing Matsuno or LF step) 

c-----------------------------------------------------------------------
