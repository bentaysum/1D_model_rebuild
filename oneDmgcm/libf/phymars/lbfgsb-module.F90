MODULE lbfgsb_module 
! Contains global parameters relevant to the optimization routines
! that emply the L-BFGS-B optimization schemes.

INTEGER t_0, t_N ! Backtrace and forecast time-steps [temporal range of optimization scheme] 
REAL*8, ALLOCATABLE :: LBFGSB_FIRSTGUESS(:) ! Mixing ratios of the first guess [i.e. the 1-D model control at backtrace timestep]

END 