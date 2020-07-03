SUBROUTINE makeinput

! Edits the run.def input file for the 1-D photochemistry submodule
! by creating a bash script and executing it.

USE lbfgsb_mod 

IMPLICIT NONE 

INTEGER, PARAMETER :: bash_unit = 10 ! Unit for the bash script  
CHARACTER(len=10) day0 ! Read at run time 
INTEGER ndt ! Number of time-steps for the 1-D model to run for 
INTEGER, PARAMETER :: spinup_sols = 10 ! Number of sols for the 1-D model to spin
INTEGER, PARAMETER :: day_step = 48 ! Time-steps per sol 
INTEGER sol_run ! Number of sols to run for -after- forecast time-step 






! 1.0) Open the bash script 
OPEN( bash_unit , FILE = "APPEND_INPUT", ACTION = "WRITE", STATUS = "REPLACE") 

! 1.1) Select day0 for the 1-D model 
WRITE(*,"(A42,F6.2)") "SOLAR LONGITUDE OF CURIOSITY DATA POINT = " , J_ls
WRITE(*,*) "CORRESPONDING DAY0 [CHECK USERMANUAL.PDF] : "
READ(*,*) day0

! 1.2) Construct the bash script 
     WRITE(bash_unit,*) "#!/bin/bash"
     ! 1.2.1) Set day0 
     WRITE(bash_unit,"(A21,A3,A2,A29,A7)") "sed -i '/day0/c\day0=", TRIM(day0), &
                    "' ", TRIM(ONED_HOME), "run.def"
     ! 1.2.2) Set ndt to an arbitrarily large value
     WRITE(bash_unit,"(A19,I3,A2,A29,A7)") "sed -i '/ndt/c\ndt=", 100, &
                    "' ", TRIM(ONED_HOME), "run.def"
     ! 1.2.3) Set latitude to 4.5 South 
     WRITE(bash_unit,"(A35,A29,A7)") "sed -i '/latitude/c\latitude=-4.5' ", &
                         TRIM(ONED_HOME), "run.def"
     ! 1.2.4) Set longitudinal mean off for MCD interpolations 
     WRITE(bash_unit,"(A40,A29,A12)") "sed -i '/long_mean/c\long_mean=.false.' ", &
                         TRIM(ONED_HOME), "callphys.def"
     ! 1.2.5) Set longitude to 137.4 
     WRITE(bash_unit,"(A38,A29,A12)") "sed -i '/longitude/c\longitude=137.4' ", &
                         TRIM(ONED_HOME), "callphys.def"
     ! 1.2.6) Define backtrace time in callphys.def 
     WRITE(bash_unit,"(A35,I3,A2,A29,A12)") "sed -i '/t_backtrace/c\t_backtrace=", t_0,"' ", &
                         TRIM(ONED_HOME), "callphys.def"
     ! 1.2.7) Define forecast time in callphys.def 
     WRITE(bash_unit,"(A35,I6,A2,A29,A12)") "sed -i '/t_forecast/c\t_forecast=", t_N,"' ", &
                         TRIM(ONED_HOME), "callphys.def"
                         
CLOSE(bash_unit)
CALL execute_command_line("chmod a+x APPEND_INPUT && ./APPEND_INPUT") 
! CALL execute_command_line("rm APPEND_INPUT") 


END SUBROUTINE