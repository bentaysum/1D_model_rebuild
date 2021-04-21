      SUBROUTINE mcd_trac(lt,alt,zls,
     & 			          co2, co, 
     &                        o2, h2,
     &                        h2ovap, h2oice)

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

!     VALUES FOR THE 1D SUBMODULE
      real co2(nlayermx), co(nlayermx), o2(nlayermx)
      real h2(nlayermx), h2ovap(nlayermx), h2oice(nlayermx)
! 	  LOCAL VARIABLES
!     === longitude ===
      real, save :: mcd_longitude 
!     === Local Time ===
      real :: local_time
!     
!     === netCDF parameters ===
      integer N_LT, N_ALT, N_ZLS, N_LAT, N_LONG
      parameter (N_LT =25, N_ALT = 30, N_ZLS = 37, 
     $           N_LAT = 41, N_LONG = 11)
      integer w, x, y, z, t ! loop iterators     
      integer retval            
      integer NCID        

      character(len=300) FILE_NAME
      integer co2id, coid, h2id, o2id, h2oid, h2oice_id
      integer altid, latid, ltid, zlsid, longid

!     === Mars Climate Database arrays ===
      real mcdv5_co2(N_LT, N_ALT, N_ZLS, N_LAT, N_LONG),  
     $     mcdv5_co(N_LT, N_ALT, N_ZLS, N_LAT, N_LONG),
     $     mcdv5_h2(N_LT, N_ALT, N_ZLS, N_LAT, N_LONG), 
     $     mcdv5_o2(N_LT, N_ALT, N_ZLS, N_LAT, N_LONG),
     $     mcdv5_h2o(N_LT, N_ALT, N_ZLS, N_LAT, N_LONG),
     $     mcdv5_h2oice(N_LT, N_ALT, N_ZLS, N_LAT, N_LONG)

      real mcdv5_lt(N_LT), mcdv5_alt(N_ALT), 
     $     mcdv5_zls(N_ZLS), mcdv5_lat(N_LAT),
     $     mcdv5_long(N_LONG)
!     === Arrays interpolated via longitude ===
      real int_long_co2(N_LT, N_ALT, N_ZLS, N_LAT), 
     $     int_long_co(N_LT, N_ALT, N_ZLS, N_LAT),
     $     int_long_h2(N_LT, N_ALT, N_ZLS, N_LAT), 
     $     int_long_o2(N_LT, N_ALT, N_ZLS, N_LAT),
     $     int_long_h2o(N_LT, N_ALT, N_ZLS, N_LAT),
     $     int_long_h2oice(N_LT, N_ALT, N_ZLS, N_LAT)

!     === Arrays interpolated via latitude ===
      real int_lat_co2(N_LT, N_ALT, N_ZLS),
     $     int_lat_co(N_LT, N_ALT, N_ZLS),
     $     int_lat_h2(N_LT, N_ALT, N_ZLS),
     $     int_lat_o2(N_LT, N_ALT, N_ZLS),
     $     int_lat_h2o(N_LT, N_ALT, N_ZLS),
     $     int_lat_h2oice(N_LT, N_ALT, N_ZLS)
!     === Arrays interpolated via altitude ===
      real, save :: fin_co2(N_LT,nlayermx,N_ZLS), 
     $              fin_co(N_LT,nlayermx,N_ZLS),
     $              fin_h2(N_LT,nlayermx,N_ZLS),
     $              fin_o2(N_LT,nlayermx,N_ZLS),
     $              fin_h2o(N_LT,nlayermx,N_ZLS),
     $              fin_h2oice(N_LT,nlayermx,N_ZLS)


!     === Arrays interpolated via zls ===
      real zco2(N_LT,nlayermx), zco(N_LT,nlayermx)
      real zo2(N_LT,nlayermx), zh2(N_LT,nlayermx),
     $     zh2o(N_LT,nlayermx), zh2oice(N_LT,nlayermx)
!     OTHER
      logical, save :: first_call = .True.

      local_time = lt*24.
!     0: IF THE MODEL IS NOT OPERATING VIA DIURNAL CHEMISTRY, LT = 12:00
      if ( .not. diurnal ) local_time = 12.
     

!=======================================================================
!     STAGE 1
!           - read in the netCDF file if this is the first time the 
!             routine is being called. The idea is: remove one dimension
!             from the T, ps, and wind arrays from the netCDF4 file by
!             interpolating with latitude. This allows the resultant 3D
!			  arrays to be interpolated repeatedly with less computatio-
! 			  nal demand.
!=======================================================================
      

      FILE_NAME = TRIM(mcd_dir) // TRIM(mcd_file)

	
      IF (first_call) THEN
      first_call = .false.
	  
		  if (zls > 359.) then
		  zls = 0.0
		  endif 
	  !     Open the netcdf file
      retval = nf_open(TRIM(FILE_NAME), NF_NOWRITE, ncid)
      if (retval .ne. nf_noerr) call get_error_hand(retval)

      
      ! Are we using longitudinal mean values?
      if ( long_mean ) then 
          ! Get IDs
           retval = nf_inq_varid(ncid, 'co2_vmr_longmean', co2id)
           if (retval .ne. nf_noerr) call get_error_hand(retval)

           retval = nf_inq_varid(ncid, 'co_vmr_longmean', coid)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
            
           retval = nf_inq_varid(ncid, 'h2_vmr_longmean', h2id)
           if (retval .ne. nf_noerr) call get_error_hand(retval)

           retval = nf_inq_varid(ncid, 'o2_vmr_longmean', o2id)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
            
           retval = nf_inq_varid(ncid, 'h2ovap_vmr_longmean', h2oid)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
           
           retval = nf_inq_varid(ncid, 'h2oice_vmr_longmean', h2oice_id)
           if (retval .ne. nf_noerr) call get_error_hand(retval)



          ! Place into the relevant sized arrays for interpolation 
          ! further down
           retval = nf_get_var_real(ncid, h2oice_id, int_long_h2oice)
           if (retval .ne. nf_noerr) call get_error_hand(retval)

           retval = nf_get_var_real(ncid, h2oid, int_long_h2o)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
           
           retval = nf_get_var_real(ncid, co2id, int_long_co2)
           if (retval .ne. nf_noerr) call get_error_hand(retval)

           retval = nf_get_var_real(ncid, h2id, int_long_h2)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
           
           retval = nf_get_var_real(ncid, o2id, int_long_o2)
           if (retval .ne. nf_noerr) call get_error_hand(retval)

           retval = nf_get_var_real(ncid, coid, int_long_co)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
           
      else 
          ! Read longitude for MCD parameters
           mcd_longitude=0. ! default value for ndt
           call getin("longitude",mcd_longitude) 
     
           ! Get IDs
           retval = nf_inq_varid(ncid, 'co2_vmr', co2id)
           if (retval .ne. nf_noerr) call get_error_hand(retval)

           retval = nf_inq_varid(ncid, 'co_vmr', coid)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
            
           retval = nf_inq_varid(ncid, 'h2_vmr', h2id)
           if (retval .ne. nf_noerr) call get_error_hand(retval)

           retval = nf_inq_varid(ncid, 'o2_vmr', o2id)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
            
           retval = nf_inq_varid(ncid, 'h2ovap_vmr', h2oid)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
           
           retval = nf_inq_varid(ncid, 'h2oice_vmr', h2oice_id)
           if (retval .ne. nf_noerr) call get_error_hand(retval)


           retval = nf_inq_varid(ncid, 'Longitude', longid)
           if (retval .ne. nf_noerr) call get_error_hand(retval)       

           ! Extract and place into relevant sized array
           retval = nf_get_var_real(ncid, co2id, mcdv5_co2)
           if (retval .ne. nf_noerr) call get_error_hand(retval)

           retval = nf_get_var_real(ncid, h2id, mcdv5_h2)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
           
           retval = nf_get_var_real(ncid, o2id, mcdv5_o2)
           if (retval .ne. nf_noerr) call get_error_hand(retval)

           retval = nf_get_var_real(ncid, coid, mcdv5_co)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
           
           retval = nf_get_var_real(ncid, h2oid, mcdv5_h2o)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
           
           retval = nf_get_var_real(ncid, h2oice_id, mcdv5_h2oice)
           if (retval .ne. nf_noerr) call get_error_hand(retval)

           retval = nf_get_var_real(ncid, longid, mcdv5_long)
           if (retval .ne. nf_noerr) call get_error_hand(retval)
           
	
      endif ! long_mean 
       
      retval = nf_inq_varid(ncid, 'Altitude', altid)
      if (retval .ne. nf_noerr) call get_error_hand(retval)
     
      retval = nf_inq_varid(ncid, 'Solar Longitude', zlsid)
      if (retval .ne. nf_noerr) call get_error_hand(retval)

      retval = nf_inq_varid(ncid, 'Local Time', ltid)
      if (retval .ne. nf_noerr) call get_error_hand(retval)

      retval = nf_inq_varid(ncid, 'Latitude', latid)
      if (retval .ne. nf_noerr) call get_error_hand(retval)

      retval = nf_get_var_real(ncid, altid, mcdv5_ALT)
      if (retval .ne. nf_noerr) call get_error_hand(retval)

      retval = nf_get_var_real(ncid, zlsid, mcdv5_ZLS)
      if (retval .ne. nf_noerr) call get_error_hand(retval)

      retval = nf_get_var_real(ncid, ltid, mcdv5_LT)
      if (retval .ne. nf_noerr) call get_error_hand(retval)

      retval = nf_get_var_real(ncid, latid, mcdv5_LAT)
      if (retval .ne. nf_noerr) call get_error_hand(retval)

!     Close the netcdf4 file
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call get_error_hand(retval)
          
      
C !=======================================================================
C !     STAGE 2
C !			- Interpolate the arrays with the submodule longitude if 
c                required
C !=======================================================================
      if ( .NOT. long_mean ) then 
      
          do x = 1, N_LT 
               do y = 1, N_ALT
                    do z = 1, N_ZLS
                         do w = 1, N_LAT  
                         
                       CALL interp_line(mcdv5_long,mcdv5_co2(x,y,z,w,:),
     $                                    N_LONG, mcd_longitude, 
     $                                    int_long_co2(x,y,z,w), 1)
                       CALL interp_line(mcdv5_long,mcdv5_co(x,y,z,w,:),
     $                                    N_LONG, mcd_longitude, 
     $                                    int_long_co(x,y,z,w), 1)
                       CALL interp_line(mcdv5_long,mcdv5_o2(x,y,z,w,:),
     $                                    N_LONG, mcd_longitude, 
     $                                    int_long_o2(x,y,z,w), 1)
                       CALL interp_line(mcdv5_long,mcdv5_h2(x,y,z,w,:),
     $                                    N_LONG, mcd_longitude, 
     $                                    int_long_h2(x,y,z,w), 1)     
                       CALL interp_line(mcdv5_long,mcdv5_h2o(x,y,z,w,:),
     $                                    N_LONG, mcd_longitude, 
     $                                    int_long_h2o(x,y,z,w), 1)
     
                    CALL interp_line(mcdv5_long,mcdv5_h2oice(x,y,z,w,:),
     $                                    N_LONG, mcd_longitude, 
     $                                    int_long_h2oice(x,y,z,w), 1)     
                         enddo 
                    enddo
               enddo
          enddo 

      endif

C !=======================================================================
C !     STAGE 3
C !			- Interpolate the arrays with the submodule latitude
C !=======================================================================

      do x = 1, N_LT 
          do y = 1, N_ALT
              do z = 1, N_ZLS 

                  call interp_line(mcdv5_lat,int_long_co2(x,y,z,:),N_LAT
     &                            ,lati(1)*180./pi,int_lat_co2(x,y,z),1)

                  call interp_line(mcdv5_lat,int_long_co(x,y,z,:),N_LAT
     &                            ,lati(1)*180./pi,int_lat_co(x,y,z),1)
	 
                  call interp_line(mcdv5_lat,int_long_h2(x,y,z,:),N_LAT
     &                            ,lati(1)*180./pi,int_lat_h2(x,y,z),1)

                  call interp_line(mcdv5_lat,int_long_o2(x,y,z,:),N_LAT
     &                             ,lati(1)*180./pi,int_lat_o2(x,y,z),1)
	 
                  call interp_line(mcdv5_lat,int_long_h2o(x,y,z,:),N_LAT
     &                            ,lati(1)*180./pi,int_lat_h2o(x,y,z),1)

              call interp_line(mcdv5_lat,int_long_h2oice(x,y,z,:),N_LAT
     &                         ,lati(1)*180./pi,int_lat_h2oice(x,y,z),1)


              end do ! of N_ZLS
          end do ! of N_ALT
      end do ! of N_LT
   
C !=======================================================================
C !     STAGE 4
C !			- Interpolate the arrays with the submodule altitude cells
C !=======================================================================

      do x = 1,N_LT
          do z = 1,N_ZLS
c               do y = 1,nlayermx

               call interp_line(mcdv5_ALT*1000.,int_lat_co2(x,:,z),N_ALT
     &                             ,alt,fin_co2(x,:,z),nlayermx)
               call interp_line(mcdv5_ALT*1000.,int_lat_co(x,:,z),N_ALT
     &                             ,alt,fin_co(x,:,z),nlayermx)
               call interp_line(mcdv5_ALT*1000.,int_lat_o2(x,:,z),N_ALT
     &                             ,alt,fin_o2(x,:,z),nlayermx)
               call interp_line(mcdv5_ALT*1000.,int_lat_h2(x,:,z),N_ALT
     &                             ,alt,fin_h2(x,:,z),nlayermx)
               call interp_line(mcdv5_ALT*1000.,int_lat_h2o(x,:,z),N_ALT
     &                             ,alt,fin_h2o(x,:,z),nlayermx)
            call interp_line(mcdv5_ALT*1000.,int_lat_h2oice(x,:,z),N_ALT
     &                             ,alt,fin_h2oice(x,:,z),nlayermx)


C               end do ! of N_ALT
          end do ! of N_ZLS
      end do ! of N_LT
      
      endif ! first_call 

      
!=======================================================================
!     STAGE 4
!			- This section and the following ones run inside the time 
!             loop of the 1D submodule.
! 			  First, we interpolate with the submodule zls
!=======================================================================

      do x = 1,N_LT
          do y = 1, nlayermx
                  call interp_line(mcdv5_zls,fin_co2(x,y,:),N_ZLS
     &                             ,zls ,zco2(x,y),1)
                  call interp_line(mcdv5_zls,fin_co(x,y,:),N_ZLS
     &                             ,zls ,zco(x,y),1)
                  call interp_line(mcdv5_zls,fin_o2(x,y,:),N_ZLS
     &                             ,zls ,zo2(x,y),1)
                  call interp_line(mcdv5_zls,fin_h2(x,y,:),N_ZLS
     &                             ,zls ,zh2(x,y),1)
                  call interp_line(mcdv5_zls,fin_h2o(x,y,:),N_ZLS
     &                             ,zls ,zh2o(x,y),1)     
                  call interp_line(mcdv5_zls,fin_h2oice(x,y,:),N_ZLS
     &                             ,zls ,zh2oice(x,y),1)     

     
         end do ! of nlayermx
      end do ! of N_LT


      
C !=======================================================================
C !     STAGE 5
C !			- Finally, interpolate with local time in the 1D submodule
C ! 			  to extract the final values
C !=======================================================================
      


      do y = 1, nlayermx
      
c         CO2 
C         ===
          IF ( mcd_co2 ) THEN

          call interp_line(mcdv5_lt, zco2(:,y), N_LT, 
     &                     local_time, co2(y), 1)
     
          co2(y) = co2(y)*mmol(igcm_co2)/43.34!mmean(1,y)
          
          ENDIF

c        H2O 
C        ===
          IF ( mcd_h2o ) THEN

          call interp_line(mcdv5_lt, zh2o(:,y), N_LT,
     &                     local_time, h2ovap(y), 1)

          h2ovap(y)  = h2ovap(y)*mmol(igcm_h2o_vap)/mmean(1,y) 

          call interp_line(mcdv5_lt, zh2oice(:,y), N_LT,
     &                     local_time, h2oice(y), 1)

          h2oice(y)  = h2oice(y)*mmol(igcm_h2o_ice)/mmean(1,y) 


          ENDIF

C         O2 
C         ===
          IF ( mcd_o2 ) THEN 
            call interp_line(mcdv5_lt, zo2(:,y), N_LT,
     &                    local_time, o2(y), 1)

            o2(y) = o2(y)*mmol(5)/mmean(1,y) 

          ENDIF 

C         H2 
C         ===
          IF ( mcd_h2 ) THEN 
            call interp_line(mcdv5_lt, zh2(:,y), N_LT,
     &                     local_time, h2(y), 1)
           h2(y)  = h2(y)*mmol(8)/mmean(1,y) 
          ENDIF

c         CO 
c         ===
          IF ( mcd_co ) THEN 
           call interp_line(mcdv5_lt, zco(:,y), N_LT,
     &                     local_time, co(y), 1)
           co(y)  = co(y)*mmol(2)/mmean(1,y) 
          ENDIF 


      end do
      
  
         
      END







      subroutine get_error_hand(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2

      END