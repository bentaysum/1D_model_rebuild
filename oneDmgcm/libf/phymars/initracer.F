      SUBROUTINE initracer(qsurf,co2ice)

       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Initialization related to tracer 
c   (transported dust, water, chemical species, ice...)
c
c   Name of the tracer
c
c   Test of dimension :
c   Initialize COMMON tracer in tracer.h, using tracer names provided
c   by the dynamics in "advtrac.h"
c
c   Old conventions: (not used any more)
c
c   If water=T : q(iq=nqmx) is the water mass mixing ratio
c     and q(iq=nqmx-1) is the ice mass mixing ratio

c   If there is transported dust, it uses iq=1 to iq=dustbin
c   If there is no transported dust : dustbin=0
c   If doubleq=T : q(iq=1) is the dust mass mixing ratio
c                  q(iq=2) is the dust number mixing ratio 

c
c   author: F.Forget
c   ------
c    Modifs: Franck Montmessin, Sebastien Lebonnois (june 2003)
c            Ehouarn Millour (oct. 2008) identify tracers by their names
c=======================================================================


#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "tracer.h"
#include "advtrac.h"
#include "comgeomfi.h"

#include "surfdat.h"

      real qsurf(ngridmx,nqmx)       ! tracer on surface (e.g.  kg.m-2)
      real co2ice(ngridmx)           ! co2 ice mass on surface (e.g.  kg.m-2)
      integer iq,ig,count
      real r0_lift , reff_lift, nueff_lift
c     Ratio of small over large dust particles (used when both 
c       doubleq and the submicron mode are active); In Montmessin
c       et al. (2002), a value of 25 has been deduced;
      real, parameter :: popratio = 25.
      character(len=20) :: txt ! to store some text

c-----------------------------------------------------------------------
c  radius(nqmx)      ! aerosol particle radius (m)
c  rho_q(nqmx)       ! tracer densities (kg.m-3)
c  alpha_lift(nqmx)  ! saltation vertical flux/horiz flux ratio (m-1)
c  alpha_devil(nqmx) ! lifting coeeficient by dust devil
c  rho_dust          ! Mars dust density
c  rho_ice           ! Water ice density
c  nuice_ref         ! Effective variance nueff of the
c                    !   water-ice size distributions
c  doubleq           ! if method with mass (iq=1) and number(iq=2) mixing ratio
c  varian            ! Characteristic variance of log-normal distribution
c-----------------------------------------------------------------------

! Initialization: get tracer names from the dynamics and check if we are
!                 using 'old' tracer convention ('q01',q02',...)
!                 or new convention (full tracer names)
      ! check if tracers have 'old' names
      count=0
      do iq=1,nqmx
        txt=" "
        write(txt,'(a1,i2.2)') 'q',iq
        if (txt.eq.tnom(iq)) then
          count=count+1
        endif
      enddo ! of do iq=1,nqmx
      
      if (count.eq.nqmx) then
        write(*,*) "initracer: tracers seem to follow old naming ",
     &             "convention (q01,q02,...)"
        write(*,*) "you should run newstart to rename them"
        stop
      endif

      ! copy tracer names from dynamics
      do iq=1,nqmx
        noms(iq)=tnom(iq)
      enddo
      
c------------------------------------------------------------
c         NAME and molar mass of the tracer 
c------------------------------------------------------------
   
! Identify tracers by their names: (and set corresponding values of mmol)
      ! 0. initialize tracer indexes to zero:
      do iq=1,nqmx
        igcm_dustbin(iq)=0
      enddo
      igcm_dust_mass=0
      igcm_dust_number=0
      igcm_ccn_mass=0
      igcm_ccn_number=0
      igcm_dust_submicron=0
      igcm_h2o_vap=0
      igcm_h2o_ice=0
      igcm_co2=0
      igcm_co=0
      igcm_o=0
      igcm_o1d=0
      igcm_o2=0
      igcm_o3=0
      igcm_h=0
      igcm_h2=0
      igcm_oh=0
      igcm_ho2=0
      igcm_h2o2=0
      igcm_ch4=0
      igcm_hcooh=0
      igcm_13ch4=0
      igcm_ch3o2=0
      igcm_ch3oh=0
      igcm_ch3=0
      igcm_hcho=0
      igcm_ch3ooh=0
      igcm_c2h6=0
      igcm_c2h5=0
      igcm_c2h5o2=0
      igcm_c2h5ooh=0
      igcm_c2h5oh=0
      igcm_ch3cho=0
      igcm_ch3cooo=0
      igcm_ch3cooh=0
      igcm_ch3coooh=0
      igcm_ch3co=0
      igcm_hco=0
      igcm_ch3o=0
      igcm_hoch2ch2o2=0
      igcm_hoch2ch2o=0
      igcm_ethgly=0
      igcm_hyetho2h=0
      igcm_ch3choho2=0
      igcm_ch3chohooh=0
      igcm_hcoch2o2=0
      igcm_glyox=0
      igcm_hcoco=0
      igcm_hooch2cho=0
      igcm_hoch2cho=0
      igcm_hochcho=0
      igcm_hoch2co=0
      igcm_hoch2co3=0
      igcm_hoch2co2h=0
      igcm_hcoco2h=0
      igcm_hoch2co3h=0
      igcm_hcoco3h=0
      igcm_hcoco3=0
      igcm_ch2choh=0
      igcm_ch2oo_e=0
      igcm_ch2oo=0  
      igcm_c2h4=0
      igcm_c2h2=0
      igcm_1ch2=0
      igcm_3ch2=0
      igcm_ch2=0
      igcm_ch=0

      igcm_c3h8=0
      igcm_ic3h7o2=0
      igcm_ic3h7ooh=0
      igcm_ch3coch3=0
      igcm_ch3coch2o2=0
      igcm_hoch2o2=0
      igcm_hoch2ooh=0
      igcm_hoch2oh=0
      igcm_c2h5o=0

!     Chlorine Scheme 22/04/2020
      igcm_cl=0
      igcm_cl2=0
      igcm_clo=0
      igcm_oclo=0
      igcm_cl2o2=0
      igcm_hcl=0
      igcm_hocl=0
      igcm_cloo=0
      igcm_ch3ocl=0
      igcm_clco=0
      igcm_clo3=0
      igcm_hclo4=0

      igcm_n2=0
      igcm_ar=0
      igcm_ar_n2=0
      igcm_n=0
      igcm_no=0
      igcm_no2=0
      igcm_n2d=0
      igcm_co2plus=0
      igcm_oplus=0
      igcm_o2plus=0
      igcm_coplus=0
      igcm_cplus=0
      igcm_nplus=0
      igcm_noplus=0
      igcm_n2plus=0
      igcm_hplus=0
      igcm_hco2plus=0
      igcm_elec=0

      ! 1. find dust tracers
      count=0
      if (dustbin.gt.0) then
        do iq=1,nqmx
          txt=" "
          write(txt,'(a4,i2.2)')'dust',count+1
          if (noms(iq).eq.txt) then
            count=count+1
            igcm_dustbin(count)=iq
            mmol(iq)=100.
          endif
        enddo !do iq=1,nqmx
      endif ! of if (dustbin.gt.0)
      if (doubleq) then
        do iq=1,nqmx
          if (noms(iq).eq."dust_mass") then
            igcm_dust_mass=iq
            count=count+1
          endif
          if (noms(iq).eq."dust_number") then
            igcm_dust_number=iq
            count=count+1
          endif
        enddo
      endif ! of if (doubleq)
      if (microphys) then
        do iq=1,nqmx
          if (noms(iq).eq."ccn_mass") then
            igcm_ccn_mass=iq
            count=count+1
          endif
          if (noms(iq).eq."ccn_number") then
            igcm_ccn_number=iq
            count=count+1
          endif
        enddo
      endif ! of if (microphys)
      if (submicron) then
        do iq=1,nqmx
          if (noms(iq).eq."dust_submicron") then
            igcm_dust_submicron=iq
            mmol(iq)=100.
            count=count+1
          endif
        enddo
      endif ! of if (submicron)
      ! 2. find chemistry and water tracers
      do iq=1,nqmx
        if (noms(iq).eq."co2") then
          igcm_co2=iq
          mmol(igcm_co2)=44.
          count=count+1
        endif
        if (noms(iq).eq."co") then
          igcm_co=iq
          mmol(igcm_co)=28.
          count=count+1
        endif
        if (noms(iq).eq."o") then
          igcm_o=iq
          mmol(igcm_o)=16.
          count=count+1
        endif
        if (noms(iq).eq."o1d") then
          igcm_o1d=iq
          mmol(igcm_o1d)=16.
          count=count+1
        endif
        if (noms(iq).eq."o2") then
          igcm_o2=iq
          mmol(igcm_o2)=32.
          count=count+1
        endif
        if (noms(iq).eq."o3") then
          igcm_o3=iq
          mmol(igcm_o3)=48.
          count=count+1
        endif
        if (noms(iq).eq."h") then
          igcm_h=iq
          mmol(igcm_h)=1.
          count=count+1
        endif
        if (noms(iq).eq."h2") then
          igcm_h2=iq
          mmol(igcm_h2)=2.
          count=count+1
        endif
        if (noms(iq).eq."oh") then
          igcm_oh=iq
          mmol(igcm_oh)=17.
          count=count+1
        endif
        if (noms(iq).eq."ho2") then
          igcm_ho2=iq
          mmol(igcm_ho2)=33.
          count=count+1
        endif
        if (noms(iq).eq."h2o2") then
          igcm_h2o2=iq
          mmol(igcm_h2o2)=34.
          count=count+1
        endif
        if (noms(iq).eq."n2") then
          igcm_n2=iq
          mmol(igcm_n2)=28.
          count=count+1
        endif
        if (noms(iq).eq."ch4") then
          igcm_ch4=iq
          mmol(igcm_ch4)=16.
          count=count+1
        endif
        if (noms(iq).eq."ch3") then
          igcm_ch3=iq
          mmol(igcm_ch3)=15.
          count=count+1
        endif
        if (noms(iq).eq."ch3o2") then
          igcm_ch3o2=iq
          mmol(igcm_ch3o2)=47.03
          count=count+1
        endif
        if (noms(iq).eq."ch3ooh") then
          igcm_ch3ooh=iq
          mmol(igcm_ch3ooh)=60.
          count=count+1
        endif
        if (noms(iq).eq."ch3oh") then
          igcm_ch3oh=iq
          mmol(igcm_ch3oh)=32.
          count=count+1
        endif
        if (noms(iq).eq."ch3o") then
          igcm_ch3o=iq
          mmol(igcm_ch3o)=31.
          count=count+1
        endif
        if (noms(iq).eq."hcho") then
          igcm_hcho=iq
          mmol(igcm_hcho)=30.
          count=count+1
        endif
        if (noms(iq).eq."hcooh") then
          igcm_hcooh=iq
          mmol(igcm_hcooh)=46.
          count=count+1
        endif
        if (noms(iq).eq."hoch2o2") then
          igcm_hoch2o2=iq
          mmol(igcm_hoch2o2)=63.
          count=count+1
        endif   
        if (noms(iq).eq."hoch2oh") then
          igcm_hoch2oh=iq
          mmol(igcm_hoch2oh)=48.
          count=count+1
        endif 

        if (noms(iq).eq."hoch2ooh") then
          igcm_hoch2ooh=iq
          mmol(igcm_hoch2ooh)=64.
          count=count+1
        endif  
        if (noms(iq).eq."hco") then
          igcm_hco=iq
          mmol(igcm_hco)=29.
          count=count+1
        endif   
        if (noms(iq).eq."c2h6") then
          igcm_c2h6=iq
          mmol(igcm_c2h6)=30.
          count=count+1
        endif
        if (noms(iq).eq."c2h5") then
          igcm_c2h5=iq
          mmol(igcm_c2h5)=29.
          count=count+1
        endif   
        if (noms(iq).eq."c2h5o2") then
          igcm_c2h5o2=iq
          mmol(igcm_c2h5o2)=61.
          count=count+1
        endif
        if (noms(iq).eq."c2h5ooh") then
          igcm_c2h5ooh=iq
          mmol(igcm_c2h5ooh)=62.
          count=count+1
        endif
        if (noms(iq).eq."c2h5oh") then
          igcm_c2h5oh=iq
          mmol(igcm_c2h5oh)=46.
          count=count+1
        endif  
        if (noms(iq).eq."hoch2ch2o2") then
          igcm_hoch2ch2o2=iq
          mmol(igcm_hoch2ch2o2)=77.
          count=count+1
        endif  
        if (noms(iq).eq."hoch2ch2o") then
          igcm_hoch2ch2o=iq
          mmol(igcm_hoch2ch2o)=61.
          count=count+1
        endif
        if (noms(iq).eq."ethgly") then
          igcm_ethgly=iq
          mmol(igcm_ethgly)=78.
          count=count+1
        endif
        if (noms(iq).eq."hyetho2h") then
          igcm_hyetho2h=iq
          mmol(igcm_hyetho2h)=78.
          count=count+1
        endif
        if (noms(iq).eq."ch3cho") then
          igcm_ch3cho=iq
          mmol(igcm_ch3cho)=44.
          count=count+1
        endif
        if (noms(iq).eq."ch2choh") then
          igcm_ch2choh=iq
          mmol(igcm_ch2choh)=44
          count=count+1
        endif
        if (noms(iq).eq."ch3choho2") then
          igcm_ch3choho2=iq
          mmol(igcm_ch3choho2)=77.
          count=count+1
        endif  
        if (noms(iq).eq."ch3cooh") then
          igcm_ch3cooh=iq
          mmol(igcm_ch3cooh)=60.
          count=count+1
        endif
        if (noms(iq).eq."ch3chohooh") then
          igcm_ch3chohooh=iq
          mmol(igcm_ch3chohooh)=78.
          count=count+1
        endif
        if (noms(iq).eq."ch3c(o)") then
          igcm_ch3co=iq
          mmol(igcm_ch3co)=43.
          count=count+1
        endif  
        if (noms(iq).eq."ch3c(o)oo") then
          igcm_ch3cooo=iq
          mmol(igcm_ch3cooo)=75.
          count=count+1
        endif
        if (noms(iq).eq."ch3c(o)ooh") then
          igcm_ch3coooh=iq
          mmol(igcm_ch3coooh)=76.
          count=count+1
        endif
        if (noms(iq).eq."hcoch2o2") then
          igcm_hcoch2o2=iq
          mmol(igcm_hcoch2o2)=75.
          count=count+1
        endif
        if (noms(iq).eq."glyox") then
          igcm_glyox=iq
          mmol(igcm_glyox)=58.
          count=count+1
        endif
        if (noms(iq).eq."hcoco") then
          igcm_hcoco=iq
          mmol(igcm_hcoco)=57.
          count=count+1
        endif
        if (noms(iq).eq."hooch2cho") then
          igcm_hooch2cho=iq
          mmol(igcm_hooch2cho)=76.
          count=count+1
        endif
        if (noms(iq).eq."hoch2cho") then
          igcm_hoch2cho=iq
          mmol(igcm_hoch2cho)=60.
          count=count+1
        endif
        if (noms(iq).eq."hochcho") then
          igcm_hochcho=iq
          mmol(igcm_hochcho)=59.
          count=count+1
        endif
        if (noms(iq).eq."hoch2co") then
          igcm_hoch2co=iq
          mmol(igcm_hoch2co)=59.
          count=count+1
        endif
        if (noms(iq).eq."hoch2co3") then
          igcm_hoch2co3=iq
          mmol(igcm_hoch2co3)=91.
          count=count+1
        endif
        if (noms(iq).eq."hoch2co2h") then
          igcm_hoch2co2h=iq
          mmol(igcm_hoch2co2h)=91.
          count=count+1
        endif
        if (noms(iq).eq."hcoco2h") then
          igcm_hcoco2h=iq
          mmol(igcm_hcoco2h)=74.
          count=count+1
        endif
        if (noms(iq).eq."hcoco3h") then
          igcm_hcoco3h=iq
          mmol(igcm_hcoco3h)=90.
          count=count+1
        endif
        if (noms(iq).eq."hcoco3") then
          igcm_hcoco3=iq
          mmol(igcm_hcoco3)=89.
          count=count+1
        endif        
        if (noms(iq).eq."hoch2co3h") then
          igcm_hoch2co3h=iq
          mmol(igcm_hoch2co3h)=92.
          count=count+1
        endif        

c     Chlorine Scheme 22/04/2020 
        if (noms(iq).eq."cl") then
          igcm_cl=iq
          mmol(igcm_cl)=35.45
          count=count+1
        endif
        if (noms(iq).eq."clo") then
          igcm_clo=iq
          mmol(igcm_clo)=51.54
          count=count+1
        endif
        if (noms(iq).eq."cl2") then
          igcm_cl2=iq
          mmol(igcm_cl2)=70.9
          count=count+1
        endif
        if (noms(iq).eq."oclo") then
          igcm_oclo=iq
          mmol(igcm_oclo)=67.45
          count=count+1
        endif
        if (noms(iq).eq."cl2o2") then
          igcm_cl2o2=iq
          mmol(igcm_cl2o2)=102.9
          count=count+1
        endif
        if (noms(iq).eq."hcl") then
          igcm_hcl=iq
          mmol(igcm_hcl)=36.46
          count=count+1
        endif
        if (noms(iq).eq."hocl") then
          igcm_hocl=iq
          mmol(igcm_hocl)=52.46
          count=count+1
        endif
        if (noms(iq).eq."cloo") then
          igcm_cloo=iq
          mmol(igcm_cloo)=67.45
          count=count+1
        endif
        if (noms(iq).eq."ch3ocl") then
          igcm_ch3ocl=iq
          mmol(igcm_ch3ocl)=66.49
          count=count+1
        endif
        if (noms(iq).eq."clco") then
          igcm_clco=iq
          mmol(igcm_clco)=63.46
          count=count+1
        endif
        if (noms(iq).eq."clo3") then
          igcm_clo3=iq
          mmol(igcm_clo3)=83.45
          count=count+1
        endif
        if (noms(iq).eq."hclo4") then
          igcm_hclo4=iq
          mmol(igcm_hclo4)=100.46
          count=count+1
        endif


        if (noms(iq).eq."ar") then
          igcm_ar=iq
          mmol(igcm_ar)=40.
          count=count+1
        endif
        if (noms(iq).eq."n") then
          igcm_n=iq
          mmol(igcm_n)=14.
          count=count+1
        endif
        if (noms(iq).eq."no") then
          igcm_no=iq
          mmol(igcm_no)=30.
          count=count+1
        endif
        if (noms(iq).eq."no2") then
          igcm_no2=iq
          mmol(igcm_no2)=46.
          count=count+1
        endif
        if (noms(iq).eq."n2d") then
          igcm_n2d=iq
          mmol(igcm_n2d)=28.
          count=count+1
        endif
        if (noms(iq).eq."co2plus") then
          igcm_co2plus=iq
          mmol(igcm_co2plus)=44.
          count=count+1
        endif
        if (noms(iq).eq."oplus") then
          igcm_oplus=iq
          mmol(igcm_oplus)=16.
          count=count+1
        endif
        if (noms(iq).eq."o2plus") then
          igcm_o2plus=iq
          mmol(igcm_o2plus)=32.
          count=count+1
        endif
        if (noms(iq).eq."coplus") then
          igcm_coplus=iq
          mmol(igcm_coplus)=28.
          count=count+1
        endif
        if (noms(iq).eq."cplus") then
          igcm_cplus=iq
          mmol(igcm_cplus)=12.
          count=count+1
        endif
        if (noms(iq).eq."nplus") then
          igcm_nplus=iq
          mmol(igcm_nplus)=14.
          count=count+1
        endif
        if (noms(iq).eq."noplus") then
          igcm_noplus=iq
          mmol(igcm_noplus)=30.
          count=count+1
        endif
        if (noms(iq).eq."n2plus") then
          igcm_n2plus=iq
          mmol(igcm_n2plus)=28.
          count=count+1
        endif
        if (noms(iq).eq."hplus") then
          igcm_hplus=iq
          mmol(igcm_hplus)=1.
          count=count+1
        endif
        if (noms(iq).eq."hco2plus") then
          igcm_hco2plus=iq
          mmol(igcm_hco2plus)=45.
          count=count+1
        endif
        if (noms(iq).eq."elec") then
          igcm_elec=iq
          mmol(igcm_elec)=1./1822.89
          count=count+1
        endif
        if (noms(iq).eq."h2o_vap") then
          igcm_h2o_vap=iq
          mmol(igcm_h2o_vap)=18.
          count=count+1
        endif
        if (noms(iq).eq."h2o_ice") then
          igcm_h2o_ice=iq
          mmol(igcm_h2o_ice)=18.
          count=count+1
        endif
        ! Other stuff: e.g. for simulations using co2 + neutral gaz
        if (noms(iq).eq."Ar_N2") then
          igcm_ar_n2=iq
          mmol(igcm_ar_n2)=30.
          count=count+1
        endif

      enddo ! of do iq=1,nqmx
      
      ! check that we identified all tracers:
      if (count.ne.nqmx) then
        write(*,*) "initracer: found only ",count," tracers"
        write(*,*) "               expected ",nqmx
        do iq=1,count
          write(*,*)'      ',iq,' ',trim(noms(iq))
        enddo
        stop
      else
        write(*,*) "initracer: found all expected tracers, namely:"
        do iq=1,nqmx
          write(*,*)'      ',iq,' ',trim(noms(iq))
        enddo
      endif

      ! if water cycle but iceparty=.false., there will nevertheless be
      ! water ice at the surface (iceparty is not used anymore, but this
      ! part is still relevant, as we want to stay compatible with the
      ! older versions).
      if (water.and.(igcm_h2o_ice.eq.0)) then
        igcm_h2o_ice=igcm_h2o_vap ! so that qsurf(i_h2o_ice) is identified
                                  ! even though there is no q(i_h2o_ice)
      else
       ! surface ice qsurf(i_h2o_ice) was loaded twice by phyetat0,
       ! as qsurf(i_h2o_vap) & as qsurf(i_h2o_ice), so to be clean:
       if (igcm_h2o_vap.ne.0) then
         qsurf(1:ngridmx,igcm_h2o_vap)=0
       endif
      endif

c------------------------------------------------------------
c     Initialisation tracers ....
c------------------------------------------------------------
      call zerophys(nqmx,rho_q)

      rho_dust=2500.  ! Mars dust density (kg.m-3)
      rho_ice=920.    ! Water ice density (kg.m-3)
      nuice_ref=0.1   ! Effective variance nueff of the
                      ! water-ice size distribution
      !!!nuice_sed=0.45   ! Sedimentation effective variance
                      ! of the water-ice size distribution

      if (doubleq) then
c       "doubleq" technique 
c       -------------------
c      (transport of mass and number mixing ratio)
c       iq=1: Q mass mixing ratio, iq=2: N number mixing ratio

        if( (nqmx.lt.2).or.(water.and.(nqmx.lt.4)) ) then
          write(*,*)'initracer: nqmx is too low : nqmx=', nqmx
          write(*,*)'water= ',water,' doubleq= ',doubleq   
        end if

        nueff_lift = 0.5
        varian=sqrt(log(1.+nueff_lift))

        rho_q(igcm_dust_mass)=rho_dust
        rho_q(igcm_dust_number)=rho_dust

c       Intermediate calcul for computing geometric mean radius r0
c       as a function of mass and number mixing ratio Q and N
c       (r0 = (r3n_q * Q/ N)^(1/3))
        r3n_q = exp(-4.5*varian**2)*(3./4.)/(pi*rho_dust)

c       Intermediate calcul for computing effective radius reff
c       from geometric mean radius r0
c       (reff = ref_r0 * r0)
        ref_r0 = exp(2.5*varian**2)
        
c       lifted dust :
c       '''''''''''
        reff_lift = 3.0e-6 !3.e-6 !Effective radius of lifted dust (m)
        alpha_devil(igcm_dust_mass)=9.e-9   !  dust devil lift mass coeff
c       alpha_lift(igcm_dust_mass)=3.0e-15  !  Lifted mass coeff
        alpha_lift(igcm_dust_mass)=1.e-6 !1.e-6 !Lifted mass coeff

        r0_lift = reff_lift/ref_r0
        alpha_devil(igcm_dust_number)=r3n_q*
     &                        alpha_devil(igcm_dust_mass)/r0_lift**3
        alpha_lift(igcm_dust_number)=r3n_q*
     &                        alpha_lift(igcm_dust_mass)/r0_lift**3

        radius(igcm_dust_mass) = reff_lift
        radius(igcm_dust_number) = reff_lift

        write(*,*) "initracer: doubleq_param reff_lift:", reff_lift
        write(*,*) "initracer: doubleq_param nueff_lift:", nueff_lift
        write(*,*) "initracer: doubleq_param alpha_lift:",
     &    alpha_lift(igcm_dust_mass)
      else

       ! initialize varian, which may be used (e.g. by surfacearea)
       ! even with conrath dust
       nueff_lift = 0.5
       varian=sqrt(log(1.+nueff_lift))

       if (dustbin.gt.1) then
        print*,'initracer: STOP!',
     $   ' properties of dust need to be set in initracer !!!'
        stop

       else if (dustbin.eq.1) then

c       This will be used for 1 dust particle size:
c       ------------------------------------------
        radius(igcm_dustbin(1))=3.e-6
        alpha_lift(igcm_dustbin(1))=0.0e-6
        alpha_devil(igcm_dustbin(1))=7.65e-9
        rho_q(igcm_dustbin(1))=rho_dust

       endif
      end if    ! (doubleq)


c     Scavenging of dust particles by H2O clouds:
c     ------------------------------------------
c     Initialize the two tracers used for the CCNs
      if (water.AND.doubleq.AND.scavenging) then
        radius(igcm_ccn_mass) = radius(igcm_dust_mass)
        alpha_lift(igcm_ccn_mass) = 1e-30
        alpha_devil(igcm_ccn_mass) = 1e-30
        rho_q(igcm_ccn_mass) = rho_dust

        radius(igcm_ccn_number) = radius(igcm_ccn_mass)
        alpha_lift(igcm_ccn_number) = alpha_lift(igcm_ccn_mass)
        alpha_devil(igcm_ccn_number) = alpha_devil(igcm_ccn_mass)
        rho_q(igcm_ccn_number) = rho_q(igcm_ccn_mass)
      endif ! of if (water.AND.doubleq.AND.scavenging)

c     Submicron dust mode:
c     --------------------

      if (submicron) then
        radius(igcm_dust_submicron)=0.1e-6
        rho_q(igcm_dust_submicron)=rho_dust
        if (doubleq) then
c         If doubleq is also active, we use the population ratio:
          alpha_lift(igcm_dust_submicron) = 
     &      alpha_lift(igcm_dust_number)*popratio*
     &      rho_q(igcm_dust_submicron)*4./3.*pi*
     &      radius(igcm_dust_submicron)**3.
          alpha_devil(igcm_dust_submicron)=1.e-30
        else
          alpha_lift(igcm_dust_submicron)=1e-6
          alpha_devil(igcm_dust_submicron)=1.e-30
        endif ! (doubleq)
      end if  ! (submicron)

c     Initialization for photochemistry:
c     ---------------------------------
      if (photochem) then
      ! initialize chemistry+water (water will be correctly initialized below)
      ! by initializing everything which is not dust ...
        do iq=1,nqmx
          txt=noms(iq)
          if (txt(1:4).ne."dust") then
            radius(iq)=0.
            alpha_lift(iq) =0.
            alpha_devil(iq)=0.
          endif
        enddo ! do iq=1,nqmx
      endif

c     Initialization for water vapor
c     ------------------------------
      if(water) then
         radius(igcm_h2o_vap)=0.
         alpha_lift(igcm_h2o_vap) =0.
         alpha_devil(igcm_h2o_vap)=0.
         if(water.and.(nqmx.ge.2)) then
           radius(igcm_h2o_ice)=3.e-6
           rho_q(igcm_h2o_ice)=rho_ice
           alpha_lift(igcm_h2o_ice) =0.
           alpha_devil(igcm_h2o_ice)=0.
         elseif(water.and.(nqmx.lt.2)) then
            write(*,*) 'nqmx is too low : nqmx=', nqmx
            write(*,*) 'water= ',water
         endif

      end if  ! (water)

c     Output for records:
c     ~~~~~~~~~~~~~~~~~~
      write(*,*)
      Write(*,*) '******** initracer : dust transport parameters :'
      write(*,*) 'alpha_lift = ', alpha_lift
      write(*,*) 'alpha_devil = ', alpha_devil
      write(*,*) 'radius  = ', radius
      if(doubleq) then
        write(*,*) 'reff_lift (um) =  ', reff_lift
        write(*,*) 'size distribution variance  = ', varian
        write(*,*) 'r3n_q , ref_r0 : ', r3n_q , ref_r0
      end if

!
!     some extra (possibly redundant) sanity checks for tracers: 
!     ---------------------------------------------------------

       if (doubleq) then 
       ! verify that we indeed have dust_mass and dust_number tracers 
         if (igcm_dust_mass.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use doubleq option without ",
     &                "a dust_mass tracer !"
           stop
         endif
         if (igcm_dust_number.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use doubleq option without ",
     &                "a dust_number tracer !"
           stop
         endif
       endif

       if ((.not.doubleq).and.(dustbin.gt.0)) then
       ! verify that we indeed have 'dustbin' dust tracers
         count=0
         do iq=1,dustbin
           if (igcm_dustbin(iq).ne.0) then
             count=count+1
           endif
         enddo
         if (count.ne.dustbin) then
           write(*,*) "initracer: error !!"
           write(*,*) "  dusbin is set to ",dustbin,
     &                " but we only have the following dust tracers:"
           do iq=1,count
             write(*,*)"   ",trim(noms(igcm_dustbin(iq)))
           enddo
           stop
         endif
       endif

       if (water) then
       ! verify that we indeed have h2o_vap and h2o_ice tracers
         if (igcm_h2o_vap.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use water option without ",
     &                "an h2o_vap tracer !"
           stop
         endif
         if (igcm_h2o_ice.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use water option without ",
     &                "an h2o_ice tracer !"
           stop
         endif
       endif

       if (scavenging) then
       ! verify that we indeed have ccn_mass and ccn_number tracers
         if (igcm_ccn_mass.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use scavenging option without ",
     &                "a ccn_mass tracer !"
           stop
         endif
         if (igcm_ccn_number.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use scavenging option without ",
     &                "a ccn_number tracer !"
           stop
         endif
       endif ! of if (scavenging)

       if (photochem .or. callthermos) then
       ! verify that we indeed have the chemistry tracers
         if (igcm_co2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a co2 tracer !"
         stop
         endif
         if (igcm_co.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a co tracer !"
         stop
         endif
         if (igcm_o.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a o tracer !"
         stop
         endif
         if (igcm_o1d.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a o1d tracer !"
         stop
         endif
         if (igcm_o2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "an o2 tracer !"
         stop
         endif
         if (igcm_o3.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "an o3 tracer !"
         stop
         endif
         if (igcm_h.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a h tracer !"
         stop
         endif
         if (igcm_h2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a h2 tracer !"
         stop
         endif
         if (igcm_oh.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "an oh tracer !"
         stop
         endif
         if (igcm_ho2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a ho2 tracer !"
         stop
         endif
         if (igcm_h2o2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a h2o2 tracer !"
         stop
         endif
         if (igcm_n2.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "a n2 tracer !"
         stop
         endif
         if (igcm_ar.eq.0) then
           write(*,*) "initracer: error !!"
           write(*,*) "  cannot use chemistry option without ",
     &                "an ar tracer !"
         stop
         endif
       endif ! of if (photochem .or. callthermos)

      end
