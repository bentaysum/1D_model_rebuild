      SUBROUTINE mcd_dust(lt,alt,zls,
     & 			      q)

      ! To acquire longitude...
      USE ioipsl_getincom 


      IMPLICIT NONE



#include "dimensions.h"
#include "dimphys.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "conc.h"
#include "netcdf.inc"
#include "comgeomfi.h"

!     INPUT PARAMETERS
      real lt            ! 1D Model local time
      real alt(nlayermx) ! 1D Model altitude layers (km)
      real zls           ! 1D Model zls
      real q(nlayermx,nqmx)
!     1D Longitude 
      real longitude
!     VALUES FOR THE 1D SUBMODULE
      real dust(nlayermx) ! Mixing ratio (kg/kg)
      real elec(nlayermx) ! Number density (e-/cm^-3)
      real localtime 

! 	  LOCAL VARIABLES
!      ---------------
      LOGICAL, save :: firscall = .TRUE.
      INTEGER, parameter :: N_LT =25   ! Local time
      INTEGER, parameter :: N_ALT = 30 ! Altitude 
      INTEGER, parameter :: N_ZLS = 37 ! Solar Longitude
      INTEGER, parameter ::N_LAT = 41 ! Latitude
      INTEGER, parameter ::N_LONG = 11! Longitude (East) 
      
      ! Required for interpolation with time 
      REAL, SAVE :: LT_mcd(N_LT), LS_mcd(N_ZLS)
      ! Static: lat, long and altitudes are invariant with time
      REAL LAT_mcd(N_LAT), LONG_mcd(N_LONG), ALT_mcd(N_ALT)
     
      ! Straight from MCD dataset
      REAL DUST_MCD(N_LT, N_ALT, N_ZLS, 
     $              N_LAT, N_LONG) 
      REAL ELEC_MCD(N_LT, N_ALT, N_ZLS, 
     $              N_LAT, N_LONG) 

      ! Interpolated with longitude/zonal mean taken 
      REAL DUST_long(N_LT, N_ALT, N_ZLS, N_LAT)
      REAL ELEC_long(N_LT, N_ALT, N_ZLS, N_LAT)
      ! Interpolated with latitude 
      REAL DUST_lat(N_LT, N_ALT, N_ZLS)
      REAL ELEC_lat(N_LT, N_ALT, N_ZLS)
      ! Interpolated with altitude 
      REAL,SAVE :: DUST_alt(N_LT, NLAYERMX, N_ZLS)
      REAL,SAVE :: ELEC_alt(N_LT, NLAYERMX, N_ZLS)

      ! Interpolated with solar longitude 
      REAL DUST_zls(N_LT,NLAYERMX)
      REAL ELEC_zls(N_LT,NLAYERMX)

      
      ! Loop iterators
      INTEGER w, x, y, z, t

!      NETCDF VARIABLES
      INTEGER retval ! Error handle
      CHARACTER(LEN=100), PARAMETER :: FILE_NAME = 
     $                   "/home/s1215319/mgcm/mcd-datafile/"
     $                 //"dustprofiles_MY24.nc"
      INTEGER NCID ! File ID 
      INTEGER VARID(7) ! Variable ID's 



      localtime = lt
!     0: IF THE MODEL IS NOT OPERATING VIA DIURNAL CHEMISTRY, LT = 12:00
      if ( .not. diurnal ) localtime = 12.
      
!=======================================================================
!     STAGE 1
!           - read in the netCDF file if this is the first time the 
!             routine is being called.
!=======================================================================
      IF ( firscall ) THEN 
          firscall = .FALSE. 
          
          ! Prevent issues with solar longitude lying outwith interpolation
          ! bounds [0., 360.]
          IF ( ZLS > 359.9 ) zls = 0. 
          
          retval = nf_open(TRIM(FILE_NAME), NF_NOWRITE, ncid)
          if (retval .ne. nf_noerr) call handle_err(retval)

          ! Variable ID Extractions 
           retval = nf_inq_varid(ncid, 'Dust MMR', varid(1))
           if (retval .ne. nf_noerr) call handle_err(retval)
           retval = nf_inq_varid(ncid, 'Local Time', varid(2))
           if (retval .ne. nf_noerr) call handle_err(retval)
           retval = nf_inq_varid(ncid, 'Altitude', varid(3))
           if (retval .ne. nf_noerr) call handle_err(retval)
           retval = nf_inq_varid(ncid, 'Solar Longitude', varid(4))
           if (retval .ne. nf_noerr) call handle_err(retval)           
           retval = nf_inq_varid(ncid, 'Latitude', varid(5))
           if (retval .ne. nf_noerr) call handle_err(retval)
           retval = nf_inq_varid(ncid, 'Longitude', varid(6))
           if (retval .ne. nf_noerr) call handle_err(retval)
           retval = nf_inq_varid(ncid, 'Electron Number Density'
     $                               , varid(7))
           if (retval .ne. nf_noerr) call handle_err(retval)
          ! Variable Value Extractions 
           retval = nf_get_var_real(ncid, varid(1), DUST_MCD)
           if (retval .ne. nf_noerr) call handle_err(retval)
           retval = nf_get_var_real(ncid, varid(2), LT_mcd)
           if (retval .ne. nf_noerr) call handle_err(retval)
           retval = nf_get_var_real(ncid, varid(3), ALT_mcd)
           if (retval .ne. nf_noerr) call handle_err(retval)
           retval = nf_get_var_real(ncid, varid(4), LS_mcd)
           if (retval .ne. nf_noerr) call handle_err(retval)
           retval = nf_get_var_real(ncid, varid(5), LAT_mcd)
           if (retval .ne. nf_noerr) call handle_err(retval)
           retval = nf_get_var_real(ncid, varid(6), LONG_mcd)
           if (retval .ne. nf_noerr) call handle_err(retval)
           
           ! Optional electron densities
           retval = nf_get_var_real(ncid, varid(7), ELEC_mcd)
           if (retval .ne. nf_noerr) call handle_err(retval)
           
           retval = nf_close(ncid)
           if (retval .ne. nf_noerr) call handle_err(retval)
     

!=======================================================================
!     STAGE 2: Interpolation with LONGITUDE 
!=======================================================================
           ! 1D Model working with specific longitude
           IF ( .NOT. LONG_MEAN ) THEN 
           
          ! Read longitude for MCD parameters
           longitude=0. ! default value for ndt
           call getin("longitude",longitude)
           
               DO x = 1, N_LT
                    DO y = 1, N_ALT
                         DO z = 1, N_ZLS 
                              DO t = 1, N_LAT 
               call interp_line(LONG_mcd,DUST_MCD(x,y,z,t,:),N_LONG
     &                             ,longitude, DUST_long(x,y,z,t)
     &                             ,1)                              
               call interp_line(LONG_mcd,ELEC_MCD(x,y,z,t,:),N_LONG
     &                             ,longitude, ELEC_long(x,y,z,t)
     &                             ,1)       
                              ENDDO
                         ENDDO 
                    ENDDO
               ENDDO
           
           ! 1D Model using zonal mean values
           ELSE 
               DO x = 1, N_LT
                    DO y = 1, N_ALT
                         DO z = 1, N_ZLS 
                              DO t = 1, N_LAT 
               DUST_long(x,y,z,t) = SUM(DUST_MCD(x,y,z,t,:))/N_LONG        
               ELEC_long(x,y,z,t) = SUM(ELEC_MCD(x,y,z,t,:))/N_LONG        
                              ENDDO
                         ENDDO 
                    ENDDO
               ENDDO           
           ENDIF 
!=======================================================================
!     STAGE 3: Interpolation with LATITUDE 
!=======================================================================
          DO x = 1, N_LT
               DO y = 1, N_ALT
                    DO z = 1, N_ZLS 
               call interp_line(LAT_mcd,DUST_long(x,y,z,:),N_lat
     &                             ,lati(1)*180./pi, DUST_lat(x,y,z)
     &                             ,1)             
               IF (igcm_elec .ne. 0 ) 
     &         call interp_line(LAT_mcd,ELEC_long(x,y,z,:),N_lat
     &                             ,lati(1)*180./pi, ELEC_lat(x,y,z)
     &                             ,1)           
                    ENDDO 
               ENDDO
          ENDDO        
!=======================================================================
!     STAGE 4: Interpolation with ALTITUDE 
!=======================================================================
          DO x = 1, N_LT
                    DO z = 1, N_ZLS 
               call interp_line(ALT_mcd,DUST_lat(x,:,z),N_alt
     &                             ,alt(:), DUST_alt(x,:,z)
     &                             ,nlayermx)   
               IF ( igcm_elec .ne. 0 ) 
     &         call interp_line(ALT_mcd,ELEC_lat(x,:,z),N_alt
     &                             ,alt(:), ELEC_alt(x,:,z)
     &                             ,nlayermx)    
                    ENDDO 
          ENDDO           
      ENDIF 


!=======================================================================
!     STAGE 5: Interpolation with SOLAR LONGITUDE 
!=======================================================================
      DO x = 1, N_LT 
          DO y = 1, nlayermx
               call interp_line(LS_mcd,DUST_alt(x,y,:),N_zls
     &                             ,zls, DUST_zls(x,y)
     &                             ,1)    
               IF ( igcm_elec .ne. 0 ) 
     &         call interp_line(LS_mcd,ELEC_alt(x,y,:),N_zls
     &                             ,zls, ELEC_zls(x,y)
     &                             ,1)    
     
     
          ENDDO 
      ENDDO 

!=======================================================================
!     STAGE 5: Interpolation with LOCAL TIME [final output]
!=======================================================================
      DO y = 1, nlayermx
               call interp_line(LT_mcd,DUST_zls(:,y),N_lt
     &                             ,localtime, dust(y)
     &                             ,1)    
               IF ( igcm_elec .ne. 0 )
     &         call interp_line(LT_mcd,ELEC_zls(:,y),N_lt
     &                             ,localtime, elec(y)
     &                             ,1)    
     
      ENDDO 

      q(:,igcm_dust_mass) = dust
      IF ( igcm_elec .ne. 0 ) q(:,igcm_elec) = elec
      
      RETURN 

      END







