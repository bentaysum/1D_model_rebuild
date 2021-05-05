SUBROUTINE tgo_profiles(ig, z, press, temp, &
						dust_mmr, &
						h2oice_mmr, &
						h2ovap_mmr, &
						rdust, &
						rice )

IMPLICIT NONE



#include "dimensions.h"
#include "dimphys.h"
#include "tracer.h"
#include "conc.h"
#include "comcstfi.h"
#include "callkeys.h"


! ===========================================
! Uses profiles of:
!	- Water Ice 
!	- Dust 
!	- Water Vapour
!	
! from measurements made by the ExoMars TGO.
! ===========================================

! ==========
! Input Vars
! ==========
INTEGER ig ! Grid point
REAL z(nlayermx) ! Altitude of 1-D Model
REAL press(nlayermx) ! Pressure (Pa)
REAL temp(nlayermx)  ! Temperature (K)
REAL dust_mmr(nlayermx) ! Dust Mass Mixing Ratio (kg/kg)
REAL h2oice_mmr(nlayermx) ! Water Ice Mass Mixing Ratio (kg/kg)
REAL h2ovap_mmr(nlayermx) ! Water Vapor Mass Mixing Ratio (kg/kg)
REAL rdust(nlayermx) ! Geometric Radius of dust (m)
REAL rice(nlayermx)  ! Geometric Radius of ice (m)

! ==========
! Local Vars
! ==========

! File ID's
INTEGER, PARAMETER :: dust_ice_ID = 110
INTEGER, PARAMETER :: h2o_ID = 120 

! Number of Altitude Steps in NOMAD Data Files 
INTEGER, PARAMETER :: Nalt = 81

! Altitude grid of NOMAD Data (0 -> 80 km , 1 km resolution)
REAL alt_grid(Nalt) !  

! Array Holding all NOMAD Data extracted from files
REAL data_grid(Nalt,5) 
! Indices for data_grid
INTEGER,PARAMETER  :: i_dustnd=1
INTEGER,PARAMETER  :: i_rdust=2
INTEGER,PARAMETER  :: i_icend=3
INTEGER,PARAMETER  :: i_rice=4
INTEGER,PARAMETER  :: i_h2ovap=5

! Ice Number density (cm-3)
! 	- 1D model requires a conversion to MMR (kg/kg) 
REAL ice_1D(nlayermx) 
! Dust Number density (cm-3)
!	- 1D model requires a conversion to MMR (kg/kg) 
REAL dust_nd(nlayermx) 

! Loop iterator
INTEGER i, l

! Holds File headers
CHARACTER(len=100) DUMMY

! Holds Unwanted Columns from H2O Files
!	- H2O files have 4 data columns made 
!	  via a 1st, 2nd, 3rd order interpolation/
!	  extrapolation routine, and "constant"
!	  that assigns constant values outside of the 
!     TGO's measured altitude range.
REAL h2ovmr_NOMAD(Nalt)
! Used to disregard unwanted columns upon 
! Read Command 
REAL dumh2o_1, dumh2o_2, dumh2o_3 

CHARACTER(len=200) NOMADH2O, NOMADICEDUST

! TGO Datafiles
! 	- pre-interpolated onto a 0->80 km altitude grid 
! --------------------------------------------------


! CHARACTER(len=*), PARAMETER :: TGO_dir = "NOMAD_ACS/" ! Main Directory of TGO Data
! CHARACTER(len=*), PARAMETER :: TGO_dust_ice = "DUST_ICE/Gridded_Data/" ! TGO Dust 
! CHARACTER(len=*), PARAMETER :: TGO_h2o  = "H2O/Gridded/" 

! CHARACTER(len=*), PARAMETER :: DUST_ICE_FILE = "gridded_aerosol_3482_I.txt"
! CHARACTER(len=*), PARAMETER :: H2O_FILE = "60Lat_240-260Ls_H2ONOMAD_MY34.txt"

! Particle Number Density -> Mass Mixing Ratio
! --------------------------------------------
REAL dust_volume, ice_volume ! Single particle volumes (m-3)
REAL dust_volcon, ice_volcon ! Volume Concentrations (m3/m3)
REAL dust_masscon, ice_masscon ! Mass Concentrations 
 
REAL airdens ! Density of Air parcel 
 
REAL,PARAMETER :: dustdens = 1.5 ! Dust Density (g cm-3) 
REAL,PARAMETER :: icedens = 0.9168  ! Ice density (g m-3)


! ==========================================

! ------------------------------
! Stage 1 : Reading TGO Profiles
! ------------------------------
OPEN( UNIT = dust_ice_ID , FILE = TRIM(NOMAD_ACS_DIR) // "/DUST_ICE/Gridded_Data/" // &
								  TRIM(NOMAD_DUST_FILE) // ".txt", ACTION = "READ" )

OPEN( UNIT = h2o_ID, FILE = TRIM(NOMAD_ACS_DIR) // "/H2O/Gridded/" // &
								  TRIM(NOMAD_H2O_FILE) // ".txt", ACTION = "READ" )
								  


READ(  dust_ice_ID, * ) DUMMY ! Skip Header 
READ( h2o_ID, * ) DUMMY ! Skip Header 

DO i = 1, Nalt

	READ( dust_ice_ID, *  ) alt_grid(i), data_grid(i,i_rdust), data_grid(i,i_rice), &
								data_grid(i,i_dustnd), data_grid(i,i_icend) 
								
	READ ( h2o_ID, * ) alt_grid(i), h2ovmr_NOMAD(i), dumh2o_1, dumh2o_2, dumh2o_3 
	
	
ENDDO



! -----------------------------------------------
! Stage 2 : Interpolating onto the 1-D Model Grid
! -----------------------------------------------

! -------------------
! Dust Number Density
! -------------------
call interp_line( alt_grid, data_grid(:,i_dustnd), Nalt, &
				  z, dust_nd, nlayermx )

! --------------------------------
! Water Vapour Volume Mixing ratio
! --------------------------------
call interp_line( alt_grid, h2ovmr_NOMAD, Nalt, &
				  z, h2ovap_mmr, nlayermx )
! VMR -> MMR 
	h2ovap_mmr = h2ovap_mmr*mmol(igcm_h2o_vap)/mmean(ig,:)

! -----------
! Dust Radius 
! -----------
call interp_line( alt_grid, data_grid(:,i_rdust), Nalt, &
				  z, rdust, nlayermx )
! um -> m 
rdust = rdust*1.e-6
 
! ----------
! Ice Radius 
! ----------
call interp_line( alt_grid, data_grid(:,i_rice), Nalt, &
				  z, rice, nlayermx )
! um -> m
rice = rice*1.e-6 

! ------------------
! Ice Number Density 
! ------------------
call interp_line( alt_grid, data_grid(:,i_icend), Nalt, &
				  z, ice_1D, nlayermx )

				  
! ! -------------------------------------------------------------------
! ! Stage 3 : Dust and Ice Particle Number Density -> Mass Mixing Ratio 
! ! -------------------------------------------------------------------
 
DO l = 1, nlayermx
	

	! Air Density at layer l (g cm-3)
	airdens = 1.e-3*press(l)/(rnew(ig,l)*temp(l))
	
	! ========
	! 3.1 Dust 
	! ========
	
	! Volume of particle [cm3] 
	dust_volume = 1.e6*(4./3.)*pi*(rdust(l)**3)
	
	! Volume Concentration [cm3/cm3]
	dust_volcon = dust_nd(l)*dust_volume
	
	! Mass Concentration 
	dust_masscon = dust_volcon*dustdens
	
	! Mass Mixing Ratio 
	dust_mmr(l) = dust_masscon/airdens 

	! =============
	! 3.2 Water Ice 
	! =============
	
	! Particle volume (cm3)
	ice_volume = 1.e6*(4./3.)*pi*(rice(l)**3)
	
	! Volume concentration
	ice_volcon = ice_1D(l)*ice_volume 
	
	! Mass Concentration
	ice_masscon = ice_volcon*icedens 
	
	! Mass Mixing Ratio 
	h2oice_mmr(l) = ice_masscon/airdens	

ENDDO

! ===========
! Close Units 
! ===========
CLOSE(dust_ice_ID)
CLOSE(h2o_ID)

				  
				  
				 
				 
END SUBROUTINE tgo_profiles