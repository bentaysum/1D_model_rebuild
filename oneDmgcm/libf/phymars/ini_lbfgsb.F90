SUBROUTINE ini_lbfgsb() 

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

IMPLICIT NONE 

! Local Variables 
! ===============
INTEGER i ! Loop iterator 

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
333 write(*,"(A54)") "--------------------------------------------------------"
! 1.1.3 : Set these values as the forecast values 
J_co = curiosity_co_array(data_point)
J_co2 = curiosity_co2_array(data_point)
J_o2 = curiosity_o2_array(data_point)
J_lt = curiosity_lt_array(data_point)
J_ls = curiosity_zls_array(data_point)



stop 

END SUBROUTINE 