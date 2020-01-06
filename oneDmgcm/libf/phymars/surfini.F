      SUBROUTINE surfini(ngrid,piceco2,qsurf,psolaralb)
   ! to use  'getin'
      USE ioipsl_getincom
      use netcdf
      IMPLICIT NONE
c=======================================================================
c
c   creation des calottes pour l'etat initial
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------
#include "dimensions.h"
#include "dimphys.h"
#include "surfdat.h"
#include "callkeys.h"
#include "tracer.h"
#include "comgeomfi.h"
#include "comcstfi.h"

#include "datafile.h"

      INTEGER ngrid,ig,icap,iq,alternate
      REAL  piceco2(ngrid),psolaralb(ngrid,2)
      REAL qsurf(ngrid,nqmx) !tracer on surface (kg/m2)
      REAL icedryness ! ice dryness
      
      ! longwatercaptag is watercaptag. Trick for some compilers
      LOGICAL, DIMENSION(100000) :: longwatercaptag

      EXTERNAL ISMIN,ISMAX
      INTEGER ISMIN,ISMAX
      
! There are 3 different modes for ice distribution:
! icelocationmode = 1 ---> based on data from surface.nc
! icelocationmode = 2 ---> directly predefined for GCM resolutions 32x24 or 64x48
! icelocationmode = 3 ---> based on logical relations for latitude and longitude
! For visualisation : > /u/tnalmd/bin/watercaps gcm_txt_output_file
      INTEGER,SAVE :: icelocationmode = 2
       
       
      !in case icelocationmode == 1
      INTEGER i,j
      INTEGER     imd,jmd
      PARAMETER   (imd=360,jmd=180)
      REAL        zdata(imd,jmd)
      REAL        zelat,zelon 

      INTEGER nb_ice(ngrid,2)              ! number of counts | detected ice for GCM grid
      INTEGER latice(jjm,2),lonice (iim,2) ! number of counts | detected ice along lat & lon axis

      REAL step,count,ratiolat

      INTEGER   ierr,nid,nvarid
      
      REAL,SAVE :: min_icevalue = 500.
      character(len=50) :: string = 'thermal'
      
      character (len=100) :: zedatafile
c
c=======================================================================

c     water ice outliers
c     ------------------------------------------

      IF ((water) .and. (caps)) THEN
     
c Perennial H20 north cap defined by watercaptag=true (allows surface to be
c hollowed by sublimation in vdifc).

c We might not want albedodat to be modified because it is used to write 
c restart files. Instead, albedo is directly modified when needed (i.e. 
c if we have watercaptag and no co2 ice), below and in albedocaps.F90

c       "Dryness coefficient" controlling the evaporation and
c        sublimation from the ground water ice (close to 1)
c        HERE, the goal is to correct for the fact
c        that the simulated permanent water ice polar caps
c        is larger than the actual cap and the atmospheric
c        opacity not always realistic.

         alternate = 0
         
         if (ngridmx .ne. 1) then
           watercaptag(:) = .false.
           longwatercaptag(:) = .false.
         endif
         
         write(*,*) "Ice dryness ?"
         icedryness=1. ! default value
         call getin("icedryness",icedryness)
         write(*,*) " icedryness = ",icedryness
         dryness (:) = icedryness
         
        
#ifdef MESOSCALE

      do ig=1,ngridmx

         !write(*,*) "all qsurf to zero. dirty."
         do iq=1,nqmx
         qsurf(ig,iq)=0.  !! on jette les inputs GCM
                          !! on regle juste watercaptag
                          !! il faudrait garder les inputs GCM 
                          !! si elles sont consequentes
         enddo
         if ( ( lati(ig)*180./pi .gt. 70. ) .and.
     .        ( albedodat(ig) .ge. 0.26   ) )  then
                 write(*,*)"outlier ",ig
                 watercaptag(ig)  = .true.
                 dryness(ig)      = 1.
                 albedodat(ig)    = albedo_h2o_ice  !! pour output 
         else
                 watercaptag(ig)  = .false.
                 dryness(ig)      = 1.
         endif 
         
      enddo
#else



      IF (ngridmx .eq. 1) THEN ! special case for 1d --> do nothing
      
         print*, 'ngridmx = 1, do no put ice caps in surfini.F'

      ELSE IF (icelocationmode .eq. 1) THEN
      
         print*,'Surfini: ice caps defined from surface.nc'
            
! This method detects ice as gridded value above min_icevalue in the field "string" from surface.nc
! Typically, it is for thermal inertia above 500 tiu.
! Two conditions are verified:
! 1. GCM ice caps are defined such as area is conserved for a given latitude
! (the approximation is that all points within the GCM latitude resolution have the same area).
! 2. caps are placed to fill the GCM points with the most detected ice first.
      

           
        zedatafile = trim(datafile)
 
        
        ierr=nf90_open(trim(zedatafile)//'/surface.nc',
     &  NF90_NOWRITE,nid)
     
      IF (ierr.NE.nf90_noerr) THEN
       write(*,*)'Error : cannot open file surface.nc '
       write(*,*)'(in phymars/surfini.F)'
       write(*,*)'It should be in :',trim(zedatafile),'/'
       write(*,*)'1) You can set this path in the callphys.def file:'
       write(*,*)'   datadir=/path/to/the/datafiles'
       write(*,*)'2) If necessary, surface.nc (and other datafiles)'
       write(*,*)'   can be obtained online on:'
       write(*,*)' http://www.lmd.jussieu.fr/~forget/datagcm/datafile'
       CALL ABORT
      ENDIF
      
      
      ierr=nf90_inq_varid(nid, string, nvarid)
      if (ierr.ne.nf90_noerr) then
        write(*,*) 'surfini error, cannot find ',trim(string)
        write(*,*) ' in file ',trim(zedatafile),'/surface.nc'
        write(*,*)trim(nf90_strerror(ierr))
        stop
      endif

      ierr=nf90_get_var(nid, nvarid, zdata)

      if (ierr.ne.nf90_noerr) then
        write(*,*) 'surfini: error failed loading ',trim(string)
        write(*,*)trim(nf90_strerror(ierr))
        stop
      endif
 
                     
      ierr=nf90_close(nid)
 

      nb_ice(:,1) = 1 ! default: there is no ice
      latice(:,1) = 1
      lonice(:,1) = 1
      nb_ice(:,2) = 0
      latice(:,2) = 0
      lonice(:,2) = 0
      !print*,'jjm,iim',jjm,iim ! jjm =  nb lati , iim = nb longi

      ! loop over the GCM grid - except for poles (ig=1 and ngridmx)
      do ig=2,ngridmx-1
      
        ! loop over the surface file grid      
        do i=1,imd
          do j=1,jmd
            zelon = i - 180.
            zelat = 90. - j 
            if ((abs(lati(ig)*180./pi-zelat) .le. 90./real(jjm)) .and.
     &        (abs(long(ig)*180./pi-zelon) .le. 180./real(iim))) then
              ! count all points in that GCM grid point
              nb_ice(ig,1) = nb_ice(ig,1) + 1
              if (zdata(i,j) > min_icevalue)
                 ! count all detected points in that GCM grid point
     &           nb_ice(ig,2) = nb_ice(ig,2) + 1
             endif
          enddo
        enddo  

        ! projection of nb_ice on GCM lat and lon axes
        latice(1+(ig-2)/iim,:) =
     &     latice(1+(ig-2)/iim,:) + nb_ice(ig,:)
        lonice(1+mod(ig-2,iim),:) = 
     &     lonice(1+mod(ig-2,iim),:) + nb_ice(ig,:) ! lonice is USELESS ...

      enddo ! of do ig=2,ngridmx-1
     

     
      ! special case for poles
      nb_ice(1,2)   = 1  ! ice prescribed on north pole
      latice(1,:)   = nb_ice(1,:)
      lonice(1,:)   = nb_ice(1,:)
      latice(jjm,:) = nb_ice(ngridmx,:)
      lonice(iim,:) = nb_ice(ngridmx,:)
      
     
!      print*, 'latice TOT', latice(:,1)
!      print*, 'latice FOUND', latice(:,2)
!      print*, 'lonice TOT', lonice(:,1)
!      print*, 'lonice FOUND', lonice(:,2)
      
!      print*, 'lat ratio', int(real(latice(:,2))/real(latice(:,1))*iim)
!      print*, 'lon ratio', int(real(lonice(:,2))/real(lonice(:,1))*jjm)
      
!      print*,''
!      print*,'sum lat', sum(latice(:,1)), sum(lonice(:,1))
!      print*,'sum lon', sum(latice(:,2)), sum(lonice(:,2))
      
    
      ! loop over GCM latitudes. CONSIDER ONLY NORTHERN HEMISPHERE
      do i=1,jjm/2
      step  = 1. ! threshold to add ice cap
      count = 0. ! number of ice GCM caps at this latitude
      ! ratiolat is the ratio of area covered by ice within this GCM latitude range
      ratiolat  = real(latice(i,2))/real(latice(i,1))
      !print*,'i',i,(i-1)*iim+2,i*iim+1
     
        ! put ice caps while there is not enough ice,
        ! as long as the threshold is above 20%
        do while ( (count .le. ratiolat*iim ) .and. (step .ge. 0.2))
          count = 0.
          ! loop over GCM longitudes
          do j=1,iim
            ! if the detected ice ratio in the GCM grid point 
            ! is more than 'step', then add ice
            if (real(nb_ice((i-1)*iim+1+j,2)) 
     &        / real(nb_ice((i-1)*iim+1+j,1)) .ge. step) then
                  watercaptag((i-1)*iim+1+j) = .true.
                  count = count + 1
            endif
          enddo ! of do j=1,iim
          !print*, 'step',step,count,ratiolat*iim
          step = step - 0.01
        enddo ! of do while
      !print*, 'step',step,count,ratiolat*iim

      enddo ! of do i=1,jjm/2
            

      ELSE IF (icelocationmode .eq. 2) THEN
      
        print*,'Surfini: predefined ice caps'
      
        if ((iim .eq. 32) .and. (jjm .eq. 24)) then ! 32x24
           
          print*,'water ice caps distribution for 32x24 resolution'
          longwatercaptag(1:9)    = .true. ! central cap - core
          longwatercaptag(26:33)  = .true. ! central cap
          longwatercaptag(1:33)  = .true. ! central cap
          longwatercaptag(56)  = .true. ! central cap
          longwatercaptag(58)  = .true. ! central cap
          longwatercaptag(60)  = .true. ! central cap
          longwatercaptag(62)  = .true. ! central cap
          longwatercaptag(64)  = .true. ! central cap
!---------------------   OUTLIERS  ----------------------------

        else if ((iim .eq. 64) .and. (jjm .eq. 48)) then ! 64x48

          print*,'water ice caps distribution for 64x48 resolution'
          longwatercaptag(1:65)   = .true. ! central cap - core
          longwatercaptag(75:85)  = .true. ! central cap 
          longwatercaptag(93:114) = .true. ! central cap
!---------------------   OUTLIERS  ----------------------------
          if (.true.) then
          longwatercaptag(136)    = .true. ! outlier, lat = 78.75
          longwatercaptag(138)    = .true. ! outlier, lat = 78.75
          longwatercaptag(140)    = .true. ! outlier, lat = 78.75
          longwatercaptag(142)    = .true. ! outlier, lat = 78.75
          longwatercaptag(161)    = .true. ! outlier, lat = 78.75
          longwatercaptag(163)    = .true. ! outlier, lat = 78.75
          longwatercaptag(165)    = .true. ! outlier, lat = 78.75
          longwatercaptag(183)    = .true. ! outlier, lat = 78.75
          longwatercaptag(185)    = .true. ! outlier, lat = 78.75
          longwatercaptag(187)    = .true. ! outlier, lat = 78.75
          longwatercaptag(189)    = .true. ! outlier, lat = 78.75
          longwatercaptag(191)    = .true. ! outlier, lat = 78.75
          longwatercaptag(193)    = .true. ! outlier, lat = 78.75
          longwatercaptag(194)    = .true. ! outlier, lat = 75
          longwatercaptag(203)    = .true. ! outlier, lat = 75
          longwatercaptag(207)    = .true. ! outlier, lat = 75
          longwatercaptag(244)    = .true. ! outlier, lat = 75
          longwatercaptag(246)    = .true. ! outlier, lat = 75
          longwatercaptag(250)    = .true. ! outlier, lat = 75
          longwatercaptag(252)    = .true. ! outlier, lat = 75
          longwatercaptag(254)    = .true. ! outlier, lat = 75
          longwatercaptag(256)    = .true. ! outlier, lat = 75
          endif
!--------------------------------------------------------------       

        !Added as a test for UK version 
        else if ((iim .eq. 72) .and. (jjm .eq. 36)) then
        
        longwatercaptag(1:72)   = .true.
        longwatercaptag(73:100) = .true.
            
        else if (ngridmx .ne. 1) then
        
      print*,'No predefined ice location for this resolution :',iim,jjm
      print*,'Please change icelocationmode in surfini.F'
      print*,'Or add some new definitions ...'
      call abort
          
        endif

        do ig=1,ngridmx
          if (longwatercaptag(ig)) watercaptag(ig) = .true.
        enddo


      ELSE IF (icelocationmode .eq. 3) THEN
      
        print*,'Surfini: ice caps defined by lat and lon values'

         do ig=1,ngridmx
         
c-------- Towards olympia planitia water caps -----------
c-------------------------------------------------------- 

       if ( ( ( lati(ig)*180./pi .ge. 77.  ) .and. ! cap #2
     .           ( lati(ig)*180./pi .le. 80.  ) .and.
     .           ( long(ig)*180./pi .ge. 110. ) .and.
     .           ( long(ig)*180./pi .le. 181. ) )
     .         .or.

     .         ( ( lati(ig)*180./pi .ge. 75.  ) .and. ! cap #4 (Korolev crater)
     .           ( lati(ig)*180./pi .le. 76.  ) .and.
     .           ( long(ig)*180./pi .ge. 150. ) .and.
     .           ( long(ig)*180./pi .le. 168. ) )
     .         .or.
     .         ( ( lati(ig)*180./pi .ge. 77 ) .and. ! cap #5
     .           ( lati(ig)*180./pi .le. 80.  ) .and.
     .           ( long(ig)*180./pi .ge. -150.) .and.
     .           ( long(ig)*180./pi .le. -110.) ) )
     .         then
             
               if ((alternate .eq. 0)) then  ! 1/2 en 64x48 sinon trop large en lat
              !    watercaptag(ig)=.true.
                  alternate = 1
               else
                  alternate = 0
               endif !end if alternate = 0
               
       endif

c----------- Opposite olympia planitia water cap --------
c-------------------------------------------------------- 

        if ( ( ( lati(ig)*180./pi     .ge.  80 ) .and.
     .         ( lati(ig)*180./pi     .le.  84 ) )
     .         .and.
     .       ( ( long(ig)*180./pi .lt. -95. ) .or.       !!! 32x24
     .         ( long(ig)*180./pi .gt.  85. ) ) ) then   !!! 32x24
!     .     ( ( ( long(ig)*180./pi .ge. -29. ) .and.       !!! 64x48
!     .         ( long(ig)*180./pi .le.  90. ) ) .or.      !!! 64x48
!     .       ( ( long(ig)*180./pi .ge. -77. ) .and.       !!! 64x48
!     .         ( long(ig)*180./pi .le. -70. ) ) ) ) then  !!! 64x48
        !   watercaptag(ig)=.true.
        endif


c -------------------- Central cap ----------------------
c-------------------------------------------------------- 

      if (abs(lati(ig)*180./pi).gt.80)
     .       watercaptag(ig)=.true.
           
c--------------------------------------------------------
c--------------------------------------------------------
      end do ! of (ngridmx)


       ELSE
      
         print*, 'In surfini.F, icelocationmode is ', icelocationmode
         print*, 'It should be 1, 2 or 3.'
         call abort 

       ENDIF ! of if (icelocation)
       
       
        ! print caps locations - useful for plots too 
        print*,'latitude | longitude | ig'
        do ig=1,ngridmx
          dryness (ig) = icedryness

          if (watercaptag(ig)) then
             print*,'ice water cap', lati(ig)*180./pi,
     .              long(ig)*180./pi, ig
          endif
        enddo
       
#endif      

       ENDIF ! (caps & water)
       

c ===============================================================
c      INITIAL ALBEDO
c ===============================================================

         write(*,*)"surfini: water frost thickness",
     s     frost_albedo_threshold
         write(*,*)"surfini: water ice albedo:", albedo_h2o_ice
         write(*,*)"surfini: water ice TI:", inert_h2o_ice

c        To start with : Initial albedo = observed dataset 
c        -------------------------------------------------
         DO ig=1,ngrid
              psolaralb(ig,1)=albedodat(ig)
              psolaralb(ig,2)=albedodat(ig)
         END DO
         PRINT*,'minimum albedo sans water caps',
     s     albedodat(ISMIN(ngrid,albedodat,1))
         PRINT*,'maximum albedo sans water caps',
     s     albedodat(ISMAX(ngrid,albedodat,1))

c        initial albedo if permanent H2O ice is present 
c        ------------------------------------------------
         IF ((water) .and. (caps)) THEN
           DO ig=1,ngrid
            IF (watercaptag(ig)) THEN
              psolaralb(ig,1) = albedo_h2o_ice
              psolaralb(ig,2) = albedo_h2o_ice
            ENDIF
           END DO
           PRINT*,'minimum albedo avec water caps',
     s     psolaralb(ISMIN(ngrid,psolaralb,1),1)
           PRINT*,'maximum albedo avec water caps',
     s     psolaralb(ISMAX(ngrid,psolaralb,1),1)

         ENDIF

c      changing initial albedo if CO2 ice is present 
c      -------------------------------------------

       DO ig=1,ngrid
         IF (piceco2(ig) .GT. 0.) THEN
             IF(ig.GT.ngrid/2+1) THEN
                icap=2
             ELSE
                icap=1
             ENDIF
             psolaralb(ig,1) = albedice(icap)
             psolaralb(ig,2) = albedice(icap)    
         END IF
       END DO

c      changing initial albedo if water ice frost is present
c      -------------------------------------------
       IF (water) THEN
          do iq=1,nqmx
c          if there is frost and surface albedo is set to albedo_h2o_ice
           if(noms(iq).eq."h2o_ice") then
             do ig=1,ngrid
             
              if ((watercaptag(ig).eqv..false.) 
     &     .and. (qsurf(ig,iq).lt.-frost_albedo_threshold)) then
              print*, ''
              print*, '!!! PROBLEM in SURFINI !!!!'
              print*, 'FOUND NEGATIVE SURFACE ICE VALUE WHERE
     & WATERCAPTAG IS FALSE'
              print*, ''
              print*, 'ig,qsurf,threshold' , 
     &         ig, qsurf(ig,iq), -frost_albedo_threshold
              print*, ''
              print*, '1) Check h2o_ice in startfi and ice
     & distribution in surfini'
              print*, '2) Use ini_h2osurf option in newstart'
              print*, ''
              CALL ABORT
            endif
             
              if ((piceco2(ig) .eq. 0.).and.
     &          (qsurf(ig,iq).gt.frost_albedo_threshold)) then
                     psolaralb(ig,1) = albedo_h2o_ice
                     psolaralb(ig,2) = albedo_h2o_ice
c                     PRINT*,'surfini.F frost',
c     &                  lati(ig)*180./pi, long(ig)*180./pi
               endif 
              enddo
           endif 
          end do
          PRINT*,'minimum albedo avec givre et co2',
     s     psolaralb(ISMIN(ngrid,psolaralb,1),1)
          PRINT*,'maximum albedo avec givre et co2',
     s     psolaralb(ISMAX(ngrid,psolaralb,1),1)
       END IF
         
      RETURN
      END
