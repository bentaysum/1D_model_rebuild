c    --------------------------------------------------------------
c    Optimization routine that calls the 1-D photochemsitry model 
c    to create an ideal initial mixing ratio concentration for to
c    minimize a specified cost function. The cost function here
c    is defined as:
c 
c    f = ABS( [1-D Surface O2] - [Curiosity Measured Surface O2] ) 
c    
c    which has a gradient:
c
c    g = d[1-D Surface O2]/d[ 1-D Input MMRS ]
c
c    which is the sensitivity vector that the Adjoint to the 1-D
c    model code calculates. 

      PROGRAM main_optimze
      
      USE lbfgsb_mod 
      
      IMPLICIT NONE 
      
c   =================================================================
c                          L-BFGS-B Variables 
c   =================================================================
c     Declare the variables needed by the code.
c       A description of all these variables is given at the end of 
c       the driver.
 
      character*60     task, csave
      logical          lsave(4)
      integer          n, m,
     +                 iwa(3*nmax), isave(44)
      double precision f, 
     +                 l(nmax), u(nmax), dsave(29), 
     +                 wa(2*mmax*nmax + 5*nmax + 11*mmax*mmax + 8*mmax)
     
    
c   =================================================================
c                   Optimization Routine Variables 
c   =================================================================

      INTEGER, PARAMETER :: bash_unit = 10 ! Unit for the bash script  
      CHARACTER(len=10) day0 ! Read at run time 
      INTEGER ndt ! Number of time-steps for the 1-D model to run for 
      INTEGER, PARAMETER :: spinup_sols = 10 ! Number of sols for the 1-D model to spin
      INTEGER, PARAMETER :: day_step = 48 ! Time-steps per sol 
      INTEGER sol_run ! Number of sols to run for -after- forecast time-step 
      CHARACTER(len=100) ONED_HOME 


      ONED_HOME = "/home/s1215319/mgcm/oneDmgcm/"


c     ----------------------------------------------------------------
c     Acquire time-steps t_0 (the backtrace timestep) and t_N (the 
c     forecast timestep), and proceed to generate a bash file that 
c     will overwrite the run.def input file of the 1-D model with the 
c     required variables :
c     day0 - Establishes the solar longitude 
c     ndt  - number of sols to run for (INCLUDING the 10 sol spin-up)
c     -----------------------------------------------------------------
      call lbfgsb_init
      

c     Open the new bash file 
      OPEN( bash_unit , FILE = "APPEND_INPUT", ACTION = "WRITE",
     $                  STATUS = "REPLACE") 
      
C    Choose day0 
C    -----------
      WRITE(*,"(A42,F6.2)") "SOLAR LONGITUDE OF CURIOSITY DATA POINT = "
     $ , J_ls
      WRITE(*,*) "CORRESPONDING DAY0 [CHECK USERMANUAL.PDF] : "
      READ(*,*) day0

C     Write overwrite commands
      WRITE(bash_unit,*) "#!/bin/bash"
      WRITE(bash_unit,"(A21,A3,A2,A29,A7)") 
     $  "sed -i '/day0/c\day0=", TRIM(day0), "' ", TRIM(ONED_HOME), 
     $  "run.def"
      WRITE(bash_unit,"(A19,I3,A2,A29,A7)") 
     $  "sed -i '/ndt/c\ndt=", t_N, "' ", TRIM(ONED_HOME), 
     $  "run.def"
      CLOSE(bash_unit)
c     Execute the bash script 
      CALL execute_command_line("chmod a+x APPEND_INPUT &&" 
     $ // "./APPEND_INPUT") 
    





      
      
      
      END PROGRAM main_optimze





c     --------------------------------------------------------------
c             DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     n is an INTEGER variable that must be set by the user to the
c       number of variables.  It is not altered by the routine.
c
c     m is an INTEGER variable that must be set by the user to the
c       number of corrections used in the limited memory matrix.
c       It is not altered by the routine.  Values of m < 3  are
c       not recommended, and large values of m can result in excessive
c       computing time. The range  3 <= m <= 20 is recommended. 
c
c     x is a DOUBLE PRECISION array of length n.  On initial entry
c       it must be set by the user to the values of the initial
c       estimate of the solution vector.  Upon successful exit, it
c       contains the values of the variables at the best point
c       found (usually an approximate solution).
c
c     l is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the lower bounds on the variables. If
c       the i-th variable has no lower bound, l(i) need not be defined.
c
c     u is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the upper bounds on the variables. If
c       the i-th variable has no upper bound, u(i) need not be defined.
c
c     nbd is an INTEGER array of dimension n that must be set by the
c       user to the type of bounds imposed on the variables:
c       nbd(i)=0 if x(i) is unbounded,
c              1 if x(i) has only a lower bound,
c              2 if x(i) has both lower and upper bounds, 
c              3 if x(i) has only an upper bound.
c
c     f is a DOUBLE PRECISION variable.  If the routine setulb returns
c       with task(1:2)= 'FG', then f must be set by the user to
c       contain the value of the function at the point x.
c
c     g is a DOUBLE PRECISION array of length n.  If the routine setulb
c       returns with taskb(1:2)= 'FG', then g must be set by the user to
c       contain the components of the gradient at the point x.
c
c     factr is a DOUBLE PRECISION variable that must be set by the user.
c       It is a tolerance in the termination test for the algorithm.
c       The iteration will stop when
c
c        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c       where epsmch is the machine precision which is automatically
c       generated by the code. Typical values for factr on a computer
c       with 15 digits of accuracy in double precision are:
c       factr=1.d+12 for low accuracy;
c             1.d+7  for moderate accuracy; 
c             1.d+1  for extremely high accuracy.
c       The user can suppress this termination test by setting factr=0.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       The user can suppress this termination test by setting pgtol=0.
c
c     wa is a DOUBLE PRECISION  array of length 
c       (2mmax + 5)nmax + 11mmax^2 + 8mmax used as workspace.
c       This array must not be altered by the user.
c
c     iwa is an INTEGER  array of length 3nmax used as
c       workspace. This array must not be altered by the user.
c
c     task is a CHARACTER string of length 60.
c       On first entry, it must be set to 'START'.
c       On a return with task(1:2)='FG', the user must evaluate the
c         function f and gradient g at the returned value of x.
c       On a return with task(1:5)='NEW_X', an iteration of the
c         algorithm has concluded, and f and g contain f(x) and g(x)
c         respectively.  The user can decide whether to continue or stop
c         the iteration. 
c       When
c         task(1:4)='CONV', the termination test in L-BFGS-B has been 
c           satisfied;
c         task(1:4)='ABNO', the routine has terminated abnormally
c           without being able to satisfy the termination conditions,
c           x contains the best approximation found,
c           f and g contain f(x) and g(x) respectively;
c         task(1:5)='ERROR', the routine has detected an error in the
c           input parameters;
c       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
c         contains additional information that the user can print.
c       This array should not be altered unless the user wants to
c          stop the run for some reason.  See driver2 or driver3
c          for a detailed explanation on how to stop the run 
c          by assigning task(1:4)='STOP' in the driver.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave  is a CHARACTER working array of length 60.
c
c     lsave is a LOGICAL working array of dimension 4.
c       On exit with task = 'NEW_X', the following information is
c         available:
c       lsave(1) = .true.  the initial x did not satisfy the bounds;
c       lsave(2) = .true.  the problem contains bounds;
c       lsave(3) = .true.  each variable has upper and lower bounds.
c
c     isave is an INTEGER working array of dimension 44.
c       On exit with task = 'NEW_X', it contains information that
c       the user may want to access:
c         isave(30) = the current iteration number;
c         isave(34) = the total number of function and gradient
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints at the current
c                         iteration;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in isave
c
c     dsave is a DOUBLE PRECISION working array of dimension 29.
c       On exit with task = 'NEW_X', it contains information that
c         the user may want to access:
c         dsave(2) = the value of f at the previous iteration;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(13) = the infinity norm of the projected gradient;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in dsave
c
c     --------------------------------------------------------------
c           END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------