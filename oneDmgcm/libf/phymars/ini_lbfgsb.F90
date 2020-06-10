SUBROUTINE ini_lbfgsb(day0, day_step, t_0, t_N) 

! BMT - 09/06/2020
! 
! Initialisation of the L-BFGS-B optimization routines. 
! -----------------------------------------------------
! 
! NOTES : 1-D model should be allowed its spin-up period
! 		  of 10 sols to enable sufficient mixing and
!		  calculations of HOx and Ox abundances.
!
! 1.0) Select Curiosity data point to steady
! 2.0) Acquire day0 point for 1-D model with regard to Curiosity datapoint selected 
! 3.0) Choose forecast time-step (time-steps AFTER 10 sol spin-up) and backtrace 
!	   time-step (time-steps AFTER 10 sol spin-up) 

IMPLICIT NONE 

! Input Variables
! ===============
INTEGER, INTENT(INOUT) :: day0  
INTEGER, INTENT(IN) :: day_step ! Time-steps per sol 
INTEGER, INTENT(OUT) :: t_0, t_N ! Backtrace and forecast time-steps
! Local Variables 
! ===============
REAL, PARAMETER :: spin_up = 10. ! Number of sols allowed for model spin-up 
INTEGER sol_backtrace 
REAL lt_backtrace

INTEGER i ! Loop iterator 

REAL day0_float 

CHARACTER(len=200), PARAMETER :: CURIOSITY_DIRECT = &
				"/exports/csce/datastore/geos/users/s1215319/paper2/curiousity_oxygen/curiosity_SAM_data/"


INTEGER array_days(668) ! Use to correctly assign day0 for the
REAL array_zls(668)     ! user selected Curiosity data point. 

CHARACTER(len=18) DUMMY 

INTEGER, PARAMETER :: N_curiosity = 19 ! Number of data points collected in trainer_tables_S1.txt 
REAL curiosity_zls_array(N_curiosity) 	   ! Solar Longitude of curiosity collected points 
REAL curiosity_lt_array(N_curiosity) 		   ! local Time of curiosity collected points 
REAL curiosity_o2_array(N_curiosity) 		   ! Curiosity collected value of surface O2 VMR
REAL curiosity_co2_array(N_curiosity) 	   ! Curiosity collected value of surface CO2 VMR [Constraint for L-BFGS-B] 
REAL curiosity_co_array(N_curiosity) 		   ! Curiosity collected value of surface CO VMR [Constraint for L-BFGS-B] 

REAL J_ls, J_lt, J_o2, J_co2, J_co ! Selected individual values of the forecast from the above arrays

INTEGER tid
REAL fp, co2_er, ar, ar_er, n2, n2_er, co_er, o2_er(2)  ! Other data not needed at present
REAL sol 

INTEGER data_point
INTEGER iostat 
CHARACTER(len=1) confirm

REAL day0s(668)
REAL zls_grid(668) ! Arrays to interpolate our day0 from 

INTEGER fortran_is_a_fussy_diva(1)

! =============================================================
! Stage One: Reading the Curiosity Data Table 
! =============================================================
111 OPEN(UNIT=10,FILE=TRIM(CURIOSITY_DIRECT)//"trainer_table_S1.txt",STATUS="OLD", ACTION="READ") 
DO i = 1, 6
	READ(10,*) DUMMY
ENDDO 
write(*,"(A54)") "===========================================================" 
write(*,*) "Curiosity data point choices:"
write(*,"(A54)") "===========================================================" 
WRITE(*,"(A4,5A10)") "", "ZLS", "CO2 VMR", "O2 VMR", "CO VMR", "LT (hrs)"
write(*,"(A54)") "-------------------------------------------------------------"
DO i = 1, 19 
	READ(10,*) curiosity_zls_array(i), sol, tid, fp, curiosity_co2_array(i), co2_er, &
														ar, ar_er, n2, n2_er, curiosity_o2_array(i), o2_er(1), o2_er(2), &
														curiosity_co_array(i), co_er
																	
	curiosity_lt_array(i) = (sol - floor(sol))*24. 
	
	write(*,"(I4,F10.2,F10.4,F10.5,F10.5, F10.2)") i, curiosity_zls_array(i), curiosity_co2_array(i),&
											curiosity_o2_array(i), curiosity_co_array(i), curiosity_lt_array(i)
ENDDO 
CLOSE(10)

! 1.1 : User select data point in presented table above to work with
 WRITE(*,*) "Point to optimize for [INTEGER 1 <= x <= 19] :"
	READ(*,"(I10)",iostat=iostat) data_point 
	IF ( data_point < 1 .or. data_point > 19 ) THEN 
		WRITE(*,*) "BETWEEN 1 AND 19" 
		GOTO 111
	ENDIF 
	IF ( iostat .ne. 0 ) GOTO 111
! 1.1.2 : Display to the user their choices and ask if they want to proceed
write(*,"(A54)") "--------------------------------------------------------"
WRITE(*,*) "User choices:" 
WRITE(*,"(A4,5A10)") "", "ZLS", "CO2 VMR", "O2 VMR", "CO VMR", "LT (hrs)"
write(*,"(I4,F10.2,F10.4,F10.5,F10.5, F10.2)") data_point, curiosity_zls_array(data_point), curiosity_co2_array(data_point),&
										curiosity_o2_array(data_point), curiosity_co_array(data_point), curiosity_lt_array(data_point)
222 write(*,*) "Correct [ y/n ] ? : "
READ(*,"(A1)") confirm
IF ( confirm == "n" .or. confirm == "N" ) GOTO 111 
IF ( confirm == "y" .or. confirm == "Y" ) GOTO 333
GOTO 222
! 1.1.3 : Set these values as the forecast values 
333 J_co = curiosity_co_array(data_point)
    J_co2 = curiosity_co2_array(data_point)
    J_o2 = curiosity_o2_array(data_point)
	J_lt = curiosity_lt_array(data_point)
	J_ls = curiosity_zls_array(data_point)

! ==========================================================================
! Stage 2 : Establishing the spatial and temporal variables of the 1-D model
! 		    which [obviously] will be fixed throughout the L-BFGS-B optimis-
!			ation process.
! ==========================================================================

! 2.1 : Create a grid of day0 vs. zls in the 1-D model to interpolate the app-
!	  ropriate day0 to begin with for this analysis.

day0s(1) = 0.
zls_grid(1) = 0.
DO i = 2, 668 
	day0s(i) = (i-1)*1.e0 
	call solarlong(day0s(i), zls_grid(i) )
	zls_grid(i) = zls_grid(i)*180./3.14
ENDDO
! 2.2 : Interpolate the value of day0 using the Curiosity data point's solar longitude.
call interp_line(zls_grid,day0s,668,J_ls,day0_float,1)

! 2.3 : Set day0 as integer (floor) of the interpolated day0_float 
call solarlong(day0_float, sol )
day0 = INT(day0_float) 
write(*,*) "day0 = ", day0,  " | Ls = ", sol*180./3.14
write(*,"(A54)") "--------------------------------------------------------"


! ==========================================================================
! Stage 3 : What time-steps in 1-D model space are we forecasting from and 
! 			backtracing towards?
! ==========================================================================

! 3.1 : Forecast time-step ( t_N ) 
! 		Model should always begin at LT = 00:00 hrs. The model will exit spin-up,
! 		run for day_step time-steps (1 sol), then continue to the temporal location
!		(lt index) of the Curiosity data point.
t_N = (spin_up + 1)*day_step + INT( (24./day_step )*J_lt )*day_step

! 3.2 : Ask for number of sols after spin-up and local time to backtrace model to 
444 write(*,*) "Sols (INT) and LT (FLOAT) after spin-up to backtrace model towards:" 
	READ(*,"(I2,F5.2)", iostat = iostat) sol_backtrace, lt_backtrace 
	IF ( iostat .ne. 0 ) GOTO 444
	IF ( lt_backtrace < 0. .or. lt_backtrace > 24. ) THEN 
		WRITE(*,*) "0. <= LT <= 24." 
		GOTO 444
	ENDIF
	
	t_0 = (spin_up)*day_step + day_step*sol_backtrace + INT( (24./day_step )*lt_backtrace )*day_step
	
	IF ( t_0 >= t_N ) THEN 
		write(*,"(A19,I4,A9, I4)") "t_0 >= t_N | t_0 = ", t_0, " , t_N = ", t_N 
		GOTO 444
	ENDIF 
write(*,"(A63)") "--------------------------------------------------------"
WRITE(*,"(A30, A3, A30)")  "BACKTRACE", " | ", "FORECAST"
WRITE(*,"(2A15, A3, 2A15)") "SOL", "LT" , " | ", "SOL", "LT"
WRITE(*,"(I15, F15.2, A3, I15, F15.2)") 1, J_lt , " | ", sol_backtrace, lt_backtrace 
write(*,*) "Correct [y/n] ? : "
READ(*,*) confirm 
IF ( confirm == "y" .or. confirm == "Y" ) GOTO 555
GOTO 444


	
555 stop 

END SUBROUTINE 