MODULE lbfgsb_module 
! Contains global parameters relevant to the optimization routines
! that emply the L-BFGS-B optimization schemes.

INTEGER t_0, t_N ! Backtrace and forecast time-steps [temporal range of optimization scheme] 
REAL*8, ALLOCATABLE :: LBFGSB_FIRSTGUESS(:) ! Mixing ratios of the first guess [i.e. the 1-D model control at backtrace timestep]

INTEGER, PARAMETER :: nmax = 1024 ! Maximum sizes for the L-BFGS-B routines
INTEGER, PARAMETER :: mmax = 17 ! "  "  "  "

REAL*8 J_o2 ! Curiosity O2 mixing ratio

REAL*8, ALLOCATABLE :: Adjoint_Transition(:,:) 
REAL*8 g_lbfgsb(nmax) ! Gradient of the cost function
END 