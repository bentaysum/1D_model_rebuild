      SUBROUTINE mcd_phys(lt,alt,zls,
     & 			      temp, ps, zonwind, merwind,
     &                    q2, ndt)

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
      integer ndt        ! number of timesteps per sol

!     VALUES FOR THE 1D SUBMODULE
      real temp(nlayermx)                         ! Temperature
      real ps                                   ! Surface Pressure
      real zonwind(nlayermx), merwind(nlayermx) ! Winds
      real q2(nlayermx)

! 	  LOCAL VARIABLES
!     === netCDF parameters ===
      integer N_LT, N_ALT, N_ZLS, N_LAT, N_LONG
      integer x, y, z, t ! loop iterators     
      integer retval            
      integer NCID        
      parameter (N_LT =25, N_ALT = 30, N_ZLS = 37, N_LAT = 41, 
     $           N_LONG = 11)
      character(len=100) FILE_NAME
      integer Tid, psid, merid, zonid, altid, latid, ltid, zlsid, Q2id,
     $        longid
!     === Mars Climate Database arrays ===
      real, save :: v5mcd_T(N_LT, N_ALT, N_ZLS, N_LAT, N_LONG), 
     $     v5mcd_ps(N_LT, N_ZLS, N_LAT, N_LONG),
     $     v5mcd_u(N_LT, N_ALT, N_ZLS, N_LAT, N_LONG), 
     $     v5mcd_v(N_LT, N_ALT, N_ZLS, N_LAT, N_LONG),
     $     v5mcd_Q2(N_LT, N_ALT, N_ZLS, N_LAT, N_LONG)
     
      real, save :: v5mcd_lt(N_LT), v5mcd_alt(N_ALT), v5mcd_zls(N_ZLS)
      real v5mcd_long(N_LONG), v5mcd_lat(N_LAT)
!     === Arrays interpolated via longitude ===
      real int_long_T(N_LT, N_ALT, N_ZLS, N_LAT), 
     $     int_long_ps(N_LT, N_ZLS, N_LAT),
     $     int_long_u(N_LT, N_ALT, N_ZLS, N_LAT), 
     $     int_long_v(N_LT, N_ALT, N_ZLS, N_LAT),
     $     int_long_Q2(N_LT, N_ALT, N_ZLS, N_LAT) 
!     === Arrays interpolated via latitude ===
      real int_lat_T(N_LT, N_ALT, N_ZLS)
      real int_lat_u(N_LT, N_ALT, N_ZLS), int_lat_v(N_LT, N_ALT, N_ZLS)
      real int_lat_Q2(N_LT, N_ALT, N_ZLS)
      real, save :: int_lat_ps(N_LT, N_ZLS)
!     === Arrays interpolated via altitude ===
      real, save ::  fin_T(N_LT, nlayermx, N_ZLS), 
     $               fin_u(N_LT, nlayermx, N_ZLS), 
     $               fin_v(N_LT, nlayermx, N_ZLS), 
     $               fin_q2(N_LT, nlayermx, N_ZLS)

!     === Arrays interpolated via zls ===
      real zT(N_LT,nlayermx), zps(N_LT)
      real zu(N_LT,nlayermx), zv(N_LT,nlayermx), zQ2(N_LT,nlayermx)
!     OTHER
      logical, save :: first_call = .True.
      logical, save :: secondcall = .False.

	  
!=======================================================================
!     STAGE 1
!           - read in the netCDF file if this is the first time the 
!             routine is being called. The idea is: remove one dimension
!             from the T, ps, and wind arrays from the netCDF4 file by
!             interpolating with latitude. This allows the resultant 3D
!			  arrays to be interpolated repeatedly with less computatio-
! 			  nal demand.
!=======================================================================
      
      FILE_NAME = "/home/s1215319/mgcm_v0.1/mcd-datafile/"
     $                 // "mcd_T_PS_U_V_Q2.nc"


      
      IF (first_call) THEN
      first_call = .FALSE.
       
      if (zls > 359.) then 
           zls = 0.0
      endif
		  
!     Open the netcdf file
      retval = nf_open(TRIM(FILE_NAME), NF_NOWRITE, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     Extract the variable ids
      if (long_mean) then ! Longitudinal mean values?
           ! ID's
           retval = nf_inq_varid(ncid, 'ps_mean', psid)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_inq_varid(ncid, 'T_mean', Tid)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_inq_varid(ncid, 'Q2_mean', Q2id)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_inq_varid(ncid, 'zon_mean', zonid)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_inq_varid(ncid, 'mer_mean', merid)
           if (retval .ne. nf_noerr) call handle_err(retval)

          ! Extraction
           retval = nf_get_var_real(ncid, psid, int_long_ps)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_get_var_real(ncid, Tid, int_long_T)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_get_var_real(ncid, q2id, int_long_q2)
           if (retval .ne. nf_noerr) call handle_err(retval)
           
           retval = nf_get_var_real(ncid, merid, int_long_v)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_get_var_real(ncid, zonid, int_long_u)
           if (retval .ne. nf_noerr) call handle_err(retval)
      
           
      else 
          ! ID's
           retval = nf_inq_varid(ncid, 'ps', psid)
           if (retval .ne. nf_noerr) call handle_err(retval)
          
           retval = nf_inq_varid(ncid, 'T', Tid)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_inq_varid(ncid, 'Q2', Q2id)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_inq_varid(ncid, 'zon_wind', zonid)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_inq_varid(ncid, 'mer_wind', merid)
           if (retval .ne. nf_noerr) call handle_err(retval)
      
          ! Extraction
           retval = nf_get_var_real(ncid, psid, v5mcd_ps)
           if (retval .ne. nf_noerr) call handle_err(retval)
           
           retval = nf_get_var_real(ncid, Tid, v5mcd_T)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_get_var_real(ncid, q2id, v5mcd_q2)
           if (retval .ne. nf_noerr) call handle_err(retval)
           
           retval = nf_get_var_real(ncid, merid, v5mcd_v)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_get_var_real(ncid, zonid, v5mcd_u)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_inq_varid(ncid, 'longitude', longid)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_get_var_real(ncid, longid, v5mcd_long)
           if (retval .ne. nf_noerr) call handle_err(retval)

      endif ! if long_mean
      
           retval = nf_inq_varid(ncid, 'solar longitude', zlsid)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_inq_varid(ncid, 'local time', ltid)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_inq_varid(ncid, 'latitude', latid)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_get_var_real(ncid, zlsid, v5mcd_ZLS)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_get_var_real(ncid, ltid, v5mcd_LT)
           if (retval .ne. nf_noerr) call handle_err(retval)

           retval = nf_get_var_real(ncid, latid, v5mcd_LAT)
           if (retval .ne. nf_noerr) call handle_err(retval)

      
      retval = nf_inq_varid(ncid, 'altitude', altid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_get_var_real(ncid, altid, v5mcd_ALT)
      if (retval .ne. nf_noerr) call handle_err(retval)
!     Close the netcdf4 file
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!=======================================================================
!     STAGE 2
!			- Interpolate the arrays with the submodule longitude if 
!                required
!=======================================================================

      if ( .not. long_mean ) then
      do x = 1, N_LT 
          do y = 1, N_ALT
              do z = 1, N_ZLS 
                   do t = 1, N_LAT

                  call interp_line(v5mcd_long,v5mcd_q2(x,y,z,t,:),N_LONG
     &                             ,long(1)*180./pi,int_long_q2(x,y,z,t)
     &                             ,1)

                   call interp_line(v5mcd_long,v5mcd_T(x,y,z,t,:),N_LONG
     &                             ,long(1)*180./pi,int_long_T(x,y,z,t)
     &                             ,1)

                   call interp_line(v5mcd_long,v5mcd_u(x,y,z,t,:),N_LONG
     &                             ,long(1)*180./pi,int_long_u(x,y,z,t)
     &                             ,1)
     
                   call interp_line(v5mcd_long,v5mcd_v(x,y,z,t,:),N_LONG
     &                             ,long(1)*180./pi,int_long_v(x,y,z,t)
     &                             ,1)
                
                   call interp_line(v5mcd_long,v5mcd_ps(x,z,t,:),N_LONG
     &                             ,long(1)*180./pi,int_long_ps(x,z,t)
     &                             ,1)
               
                  end do ! of N_LAT
              end do ! of N_ZLS
          end do ! of N_ALT
      end do ! of N_LT
      endif
      
!=======================================================================
!     STAGE 3
!			- Interpolate the arrays with the submodule latitude
!=======================================================================
      do x = 1, N_LT 
          do y = 1, N_ALT
              do z = 1, N_ZLS 

                   call interp_line(v5mcd_lat,int_long_q2(x,y,z,:),N_LAT
     &                             ,lati(1)*180./pi,int_lat_q2(x,y,z),1)

                  call interp_line(v5mcd_lat,int_long_T(x,y,z,:),N_LAT
     &                             ,lati(1)*180./pi,int_lat_T(x,y,z),1)
                  
                  call interp_line(v5mcd_lat,int_long_u(x,y,z,:),N_LAT
     &                             ,lati(1)*180./pi,int_lat_u(x,y,z),1) 

                  call interp_line(v5mcd_lat,int_long_v(x,y,z,:),N_LAT
     &                             ,lati(1)*180./pi,int_lat_v(x,y,z),1) 

                  call interp_line(v5mcd_lat,int_long_ps(x,z,:),N_LAT
     &                             ,lati(1)*180./pi,int_lat_ps(x,z),1) 
               
              end do ! of N_ZLS
          end do ! of N_ALT
      end do ! of N_LT

!=======================================================================
!     STAGE 4
!			- Acquire surface pressure then return to testphys1d to
!                enable vertical grid construction
!=======================================================================
      
      do x = 1,N_LT
           
          call interp_line(v5mcd_zls,int_lat_ps(x,:),N_ZLS
     &                       ,zls,zps(x),1) 
      enddo
          
      call interp_line(v5mcd_lt, zps(:), N_LT,
     &                     lt, ps, 1)

       secondcall = .True. 
           
       RETURN 
      
      ENDIF ! of .firstcall. 
      
      
C !=======================================================================
C !     STAGE 5
C !			- Interpolate the arrays with the submodule altitude cells
C !=======================================================================
      
      IF ( secondcall ) THEN 
      
      secondcall = .False. 

      do x = 1,N_LT
          do z = 1,N_ZLS

                call interp_line(v5mcd_ALT*1000.,int_lat_T(x,:,z),N_ALT
     &                             ,alt,fin_T(x,:,z),nlayermx)

                call interp_line(v5mcd_ALT*1000.,int_lat_q2(x,:,z),N_ALT
     &                             ,alt,fin_q2(x,:,z),nlayermx)

                call interp_line(v5mcd_ALT*1000.,int_lat_u(x,:,z),N_ALT
     &                             ,alt,fin_u(x,:,z),nlayermx) 

                call interp_line(v5mcd_ALT*1000.,int_lat_v(x,:,z),N_ALT
     &                             ,alt,fin_v(x,:,z),nlayermx) 

          end do ! of N_ZLS
      end do ! of N_LT
      
      endif  ! end of .secondcall.
       
C !=======================================================================
C !     STAGE 5
C !			- This section and the following ones run inside the time 
C !             loop of the 1D submodule.
C ! 			  First, we interpolate with the submodule zls
C !=======================================================================


      do x = 1,N_LT
      
C           call interp_line(v5mcd_zls,int_lat_ps(x,:),N_ZLS
C      &                       ,zls,zps(x),1) 
      
          do y = 1, nlayermx
                  call interp_line(v5mcd_zls,fin_T(x,y,:),N_ZLS
     &                             ,zls ,zT(x,y),1)

                   call interp_line(v5mcd_zls,fin_q2(x,y,:),N_ZLS
     &                             ,zls ,zq2(x,y),1)
   
                  call interp_line(v5mcd_zls,fin_u(x,y,:),N_ZLS
     &                             ,zls,zu(x,y),1) 

                  call interp_line(v5mcd_zls,fin_v(x,y,:),N_ZLS
     &                             ,zls,zv(x,y),1) 

         end do ! of nlayermx
      end do ! of N_LT

C !=======================================================================
C !     STAGE 5
C !			- Finally, interpolate with local time in the 1D submodule
C ! 			  to extract the final values
C !=======================================================================
      

C       call interp_line(v5mcd_lt, zps(:), N_LT,
C      &                     lt, ps, 1)

      do y = 1, nlayermx
		 
          call interp_line(v5mcd_lt, zT(:,y), N_LT, lt, temp(y), 1)
          
          call interp_line(v5mcd_lt, zq2(:,y), N_LT,
     &                     lt, q2(y), 1)
          
          call interp_line(v5mcd_lt, zu(:,y), N_LT,
     &                     lt, zonwind(y), 1)
          
          call interp_line(v5mcd_lt, zv(:,y), N_LT,
     &                     lt, merwind(y), 1)

      end do
      


      END













      subroutine handle_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2

      END
      