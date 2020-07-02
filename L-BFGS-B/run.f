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

      PROGRAM main_optimize
      
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
      character(len=30) dummy_1, dummy_2 ! Dummy variables for file parsing 
      integer line ! Loop iterator 
      real finalstate(nqmx*nlayermx) ! Final state from 1-D model run with X as input 
      real initialstate(nqmx*nlayermx) ! Initial state in REAL*4 format
      
      logical lbfgsb_firstcall
      logical Existance
      character(len=100) filecheck
      
      integer iq ! Iterator 
      integer lyr ! Iterator
      integer iql ! Iterator
      
      CHARACTER(LEN=*), PARAMETER :: BENSOUTPUT = "bensoutput.txt"
      INTEGER activate
      INTEGER errorstat
      INTEGER filenumber 
      CHARACTER(len=3) filenumber_string
c     ----------------------------------------------------------------
c     Acquire time-steps t_0 (the backtrace timestep) and t_N (the 
c     forecast timestep), and proceed to generate a bash file that 
c     will overwrite the run.def input file of the 1-D model with the 
c     required variables :
c     day0 - Establishes the solar longitude 
c     ndt  - number of sols to run for (INCLUDING the 10 sol spin-up)
c     -----------------------------------------------------------------


c     Acquisition of tracer names from traceur.def 
      OPEN(200,FILE = ONED_HOME // "traceur.def", ACTION = "READ") 
      READ(200,*) 
      DO iq = 1, nqmx 
          READ(200,*) noms(iq) 
          write(*,*) TRIM(NOMS(IQ)), " | ", iq 
      ENDDO 
      
      lbfgsb_firstcall = .TRUE.

c     Purge the 1-D model directory of temporarily stored .dat files:
c     - grad.dat 
c     - finalstate.dat 
c     - inputstate.dat 
c     prior to routine starting. These files are searched in the 1-D 
c     model to activate certain conditionals and must not be present
c     on the very first call. 
      filecheck = "/home/s1215319/mgcm/oneDmgcm/finalstate.dat"
      INQUIRE(file = filecheck, EXIST = existance) 
      IF ( existance ) call system("rm " // filecheck)

      filecheck = "/home/s1215319/mgcm/oneDmgcm/grad.dat"
      INQUIRE(file = filecheck, EXIST = existance) 
      IF ( existance ) call system("rm " // filecheck)
      
      filecheck = "/home/s1215319/mgcm/oneDmgcm/inputstate.dat"
      INQUIRE(file = filecheck, EXIST = existance) 
      IF ( existance ) call system("rm " // filecheck)

      call lbfgsb_init
      
      call makeinput 
      
C     ===============
C     Very First Call 
C     ===============
C     Need to extract the 1-D atmospheric state at the specified backtrace 
c     time-step to pose as the first guess of X for the L-BFGS-B routines. 
      call system("cd /home/s1215319/mgcm/oneDmgcm && "
     $ // "./testphys1d.e ", errorstat)
      IF ( errorstat .ne. 0 ) CALL testphys1d_error
      
      filenumber = 0
      write(filenumber_string,"(I1)") filenumber
      call system("mv /scratch/local/s1215319/organics/diagfi.nc " 
     $ // "/exports/csce/datastore/geos/users/s1215319/L-BFGS-B/"
     $ // "temphold/diagfi_" // TRIM(filenumber_string) // ".nc"
     $ ,errorstat)
      IF ( errorstat .ne. 0 ) CALL testphys1d_error
     
      call system("cp /home/s1215319/mgcm/oneDmgcm/finalstate.dat " 
     $ // "/exports/csce/datastore/geos/users/s1215319/L-BFGS-B/"
     $ // "temphold/final_" // TRIM(filenumber_string) // ".dat"
     $ ,errorstat)
      IF ( errorstat .ne. 0 ) CALL testphys1d_error
      
c     =====================
c     Define the dimensions
c     =====================
      n = nqmx*nlayermx
      m = 15

c     =========================
c     Extract the first X state
c     =========================
c     File Layout:
c     / tracer name / model layer / mmr at t_0 / mmr at t_N /  
      OPEN(unit = 50, 
     $     file = "/home/s1215319/mgcm/oneDmgcm/finalstate.dat",
     $     action = "READ")
      READ(50,*) dummy_1 ! First line = Header 
      iq = 1 
      DO line = 1, nqmx*nlayermx 
          READ(50,"(A15,A10,2E15.7)") DUMMY_1, DUMMY_2, 
     $            initialstate(line), finalstate(line)
     
c     =============================
c     Assign upper and lower bounds  
c     =============================
          X(line) = MAX(1.e-31,DBLE(initialstate(line))) ! Double conversion for L-BFGS-B
          nbd(line) = 2

          l(line) = MAX(X(line) - X(line)*0.9D0,1.D-31) 
          u(line) = MIN(X(line) + X(line)*10.D0,0.99) 
          
      ENDDO 
      CLOSE(50)
      
c     =========================
c     Optimization routine loop 
c     =========================
      task = 'START' 
      f = 0.d0 
      OPEN( UNIT = 222, FILE = BENSOUTPUT, ACTION = "WRITE",
     $      STATUS = "REPLACE" )
      activate = 0 
      
      write( 222 , "(A9,2A23, 16A23)" ) "TASK",
     $     "f", "MAX[g]",  ( noms(iq) , iq = 1, nqmx ) 
      write( 222,"(A368)") 
     $            "----------------------------------------------------"
     $          //"----------------------------------------------------"
     $          //"----------------------------------------------------"
     $          //"----------------------------------------------------"
     $          //"----------------------------------------------------"
     $          //"----------------------------------------------------"
     $          //"----------------------------------------------------"
     $          //"----"

      
111   CONTINUE 
      
      call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,
     +            csave,lsave,isave,dsave)
     
              activate = activate + 1
              WRITE(222,"(A9,2E23.15,16E23.15)") TASK(1:7),f, maxval(g),
     $       (X( (iq-1)*nlayermx + 1 ) , iq = 1,nqmx )
     
Cccccccccccccccccc
cc Incomplete
cccccccccccccccccc     

      if ( task(1:2) .eq. 'FG' ) then
              
          ! On the first ask for F and G values, all we need to do is
          ! read files made by the call to the 1-D model above. 
          if ( lbfgsb_firstcall ) then 
          
               lbfgsb_firstcall = .FALSE. 
             
               call COST(f)
               call GRAD 
               

               
               
               goto 111 
               
          else 
          
          OPEN( UNIT = 50, 
     $          FILE = "/home/s1215319/mgcm/oneDmgcm/inputstate.dat",
     $          ACTION = "WRITE",
     $          STATUS = "REPLACE")
          
          DO iq = 1, nqmx 
          
               DO lyr = 1, nlayermx
               
               write(50,"(A15,I10,E15.7)" ) adjustl(noms(iq)), 
     $                     lyr,
     $                     real(X((iq-1)*nlayermx + lyr))
               ENDDO 
               
          ENDDO 
          
          close(50)
          
          call system("cd /home/s1215319/mgcm/oneDmgcm &&"
     $                // " ./testphys1d.e ", errorstat)
          
          if ( errorstat .ne. 0 ) call testphys1d_error
          
      filenumber = filenumber + 1
      write(filenumber_string,"(I3)") filenumber
      write(*,*) "mv /scratch/local/s1215319/organics/diagfi.nc " 
     $ // "/exports/csce/datastore/geos/users/s1215319/L-BFGS-B/"
     $ // "temphold/diagfi_" // TRIM(ADJUSTL(filenumber_string))//".nc"
      call system("mv /scratch/local/s1215319/organics/diagfi.nc " 
     $ // "/exports/csce/datastore/geos/users/s1215319/L-BFGS-B/"
     $ // "temphold/diagfi_" // TRIM(ADJUSTL(filenumber_string)) //".nc"
     $ ,errorstat)
     
      IF ( errorstat .ne. 0 ) CALL testphys1d_error          
      
          call COST(f)
          call GRAD 
          
               
          
          goto 111 
                   
          endif 
         
      endif 
      
ccccccccccccccccccc
cc Incomplete
ccccccccccccccccccc
      if ( task(1:5) .eq. 'NEW_X' ) then 
          OPEN( UNIT = 50, 
     $          FILE = "/home/s1215319/mgcm/oneDmgcm/inputstate.dat",
     $          ACTION = "WRITE",
     $          STATUS = "REPLACE")
          
          DO iq = 1, nqmx 
          
               DO lyr = 1, nlayermx
               
               write(50,"(A15,I10,E15.7)" ) adjustl(noms(iq)), 
     $                     lyr,
     $                     real(X((iq-1)*nlayermx + lyr))
               ENDDO 
               
          ENDDO 
          
          close(50)
          
          GOTO 111
          
      endif 

     
      END PROGRAM main_optimize





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

      SUBROUTINE COST(f)
      
      USE lbfgsb_mod
      
      IMPLICIT NONE 
      
      ! INPUT 
      ! =====
      REAL*8, INTENT(OUT) :: f 
      
      ! LOCAL 
      ! =====
      INTEGER i ! loop iterator 
      CHARACTER(len = 100) dummy_1
      CHARACTER(len=15) tracer 
      INTEGER layer 
      REAL*4 initial, final 
      
      ! ------------------------------------------------------------
      OPEN(UNIT=20,FILE="/home/s1215319/mgcm/oneDmgcm/finalstate.dat"
     $    ,ACTION="READ")
      READ(20,*) dummy_1 
      DO i = 1, nqmx*nlayermx
         READ(20,"(A15,I10,2E15.7)") tracer, layer, initial, final
      
      IF ( ADJUSTL(tracer) == "o2" ) THEN 
          
          IF ( layer == 1 ) THEN 
               
               f = ABS(final*1.D0-J_o2)
               close(20)
               RETURN 
               
          ENDIF 
      
      ENDIF 
          
      ENDDO 
      
      WRITE(*,*) "ERROR IN COST FUNCTION : SHOULDN'T REACH HERE" 
      
      STOP 
       
      END SUBROUTINE COST 
      
      
      SUBROUTINE GRAD
      
      USE lbfgsb_mod
      
      IMPLICIT NONE 

      ! Variables
      ! ---------
      CHARACTER(len=100) dummy 
      CHARACTER(len=15) dummy_1
      CHARACTER(len=10) dummy_2
      INTEGER i ! Loop iterator 
      
      OPEN(UNIT=20,FILE="/home/s1215319/mgcm/oneDmgcm/grad.dat")
      
      READ(20,*) dummy 
      READ(20,*) dummy
      DO i = 1 ,nqmx*nlayermx
          READ(20,"(A15,A10,E23.15)") dummy_1, dummy_2, g(i)
          IF ( ABS(g(i)) < 1.e-6 ) g(i) = 0.D0
          
          IF ( adjustl(dummy_1) == "o2" ) THEN
               g(i) = 0.D0
          ELSE  
               g(i) = g(i)
          ENDIF 
          
      ENDDO 
      
      RETURN 
      
      END SUBROUTINE GRAD 
      
      
      
      
      
      
      
      SUBROUTINE testphys1d_error 
      
      IMPLICIT NONE 

      WRITE(*,*) " ERROR IN TESTPHYS1D " 
      STOP

      END SUBROUTINE testphys1d_error
     