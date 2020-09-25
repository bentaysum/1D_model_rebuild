 subroutine inichim_newstart(pq, qsurf, ps, flagh2o, flagthermo, zls)

      USE ioipsl_getincom 

      implicit none

!=======================================================================
!
!  Purpose:
!  --------
!
!  Initialization of the chemistry for use in the building of a new start file
!  used by program newstart.F
!  also used by program testphys1d.F
!
!  Authors: 
!  -------
!  Initial version 11/2002 by Sebastien Lebonnois
!  Revised 07/2003 by Monica Angelats-i-Coll to use input files
!  Modified 10/2008 Identify tracers by their names Ehouarn Millour
!  Modified 11/2011 Addition of methane Franck Lefevre
!  Rewritten 04/2012 Franck Lefevre
!
!  Arguments:
!  ----------
!
!  pq(iip1,jjp1,llm,nqmx)  Advected fields, ie chemical species here
!  qsurf(ngridmx,nqmx)     Amount of tracer on the surface (kg/m2)
!  ps(iip1,jjp1)           Surface pressure (Pa)
!  flagh2o                 flag for initialisation of h2o (1: yes / 0: no)
!  flagthermo              flag for initialisation of thermosphere only (1: yes / 0: no)
!
!=======================================================================

#include "dimensions.h"
#include "dimphys.h"
#include "paramet.h"
#include "tracer.h"
#include "comvert.h"
#include "callkeys.h"
#include "datafile.h"

! inputs :

      real,intent(in) :: ps(iip1,jjp1)            ! surface pressure in the gcm (Pa)   
      integer,intent(in) :: flagh2o               ! flag for h2o initialisation
      integer,intent(in) :: flagthermo            ! flag for thermosphere initialisation only
      real, intent(in) :: zls                     ! solar longitude (degrees)
! outputs :

      real,intent(out) :: pq(iip1,jjp1,llm,nqmx)  ! advected fields, ie chemical species
      !PIP
      real,intent(out) :: qsurf(ngridmx,nqmx)     ! surface values (kg/m2) of tracers

! local :

      integer :: iq, i, j, l, n, nbqchem
      integer :: count, ierr, dummy
      real    :: mmean(iip1,jjp1,llm)             ! mean molecular mass (g)
      real    :: pgcm                             ! pressure at each layer in the gcm (Pa)
      real    :: real_dummy

      integer, parameter         :: nalt = 252    ! number of altitudes in the initialization files
      integer                    :: nspe          ! number of species in the initialization files
      integer, allocatable       :: niq(:)        ! local index of species in initialization files
      real, dimension(nalt)      :: tinit, zzfile ! temperature in initialization files
      real, dimension(nalt)      :: pinit         ! pressure in initialization files
      real, dimension(nalt)      :: densinit      ! total number density in initialization files
      real, allocatable          :: vmrinit(:,:)  ! mixing ratios in initialization files
      real, allocatable          :: vmrint(:)     ! mixing ratio interpolated onto the gcm vertical grid
      real                       :: vmr

      character(len=20)          :: txt           ! to store some text
      logical                    :: flagnitro     ! checks if N species present

      real :: ch4vmr

  
! 1. identify tracers by their names: (and set corresponding values of mmol)

! 1.1 initialize tracer indexes to zero:

      do iq = 1,nqmx
        igcm_dustbin(iq) = 0
      end do

      igcm_dust_mass   = 0
      igcm_dust_number = 0
      igcm_ccn_mass    = 0
      igcm_ccn_number  = 01
      igcm_h2o_vap     = 0
      igcm_h2o_ice     = 0
      igcm_co2         = 0
      igcm_co          = 0
      igcm_o           = 0
      igcm_o1d         = 0
      igcm_o2          = 0
      igcm_o3          = 0
      igcm_h           = 0
      igcm_h2          = 0
      igcm_oh          = 0
      igcm_ho2         = 0
      igcm_h2o2        = 0
      igcm_ch4         = 0
      igcm_ch3o2       = 0 ! Added 24/04/2018
      igcm_ch3         = 0 ! Added 25/04/2018
      igcm_ch3oh       = 0 ! Added 27/04/2018
      igcm_hcho        = 0 ! Added 27/04/2018
      igcm_ch3ooh      = 0 ! Added 27/04/2018
      igcm_c2h6        = 0 ! Added 12/06/2018
	  igcm_hoch2ch2o2  = 0
	  igcm_hoch2ch2o   = 0
	  igcm_ethgly 	   = 0
	  igcm_hyetho2h    = 0
	  igcm_ch3chohooh  = 0
	  igcm_hcoch2o2    = 0 
	  igcm_glyox       = 0
	  igcm_hcoco       = 0
	  igcm_hooch2cho   = 0
	  igcm_hoch2cho    = 0
	  igcm_hochcho 	   = 0
	  igcm_hoch2co     = 0
	  igcm_hoch2co3    = 0
	  igcm_hoch2co2h   = 0
	  igcm_hcoco2h     = 0
	  igcm_hoch2co3h   = 0
	  igcm_hcoco3h     = 0
	  igcm_hcoco3 	   = 0
      igcm_ch2choh     = 0
      igcm_c2h4        = 0 ! Added 12/06/2018
      igcm_c2h5o2      = 0 ! Added 12/06/2018
      igcm_c2h2        = 0 ! Added 13/06/2018
      igcm_hcooh       = 0 ! Added 08/07/2018
      igcm_1ch2        = 0 ! Added 28/05/2018
      igcm_3ch2        = 0 ! Added 28/05/2018
      igcm_ch          = 0 ! Added 28/05/2018
      igcm_ch2          = 0 ! Added 28/05/2018
	  igcm_ch3o        = 0
	  igcm_c2h5ooh     = 0
	  igcm_ch3cho      = 0
	  igcm_ch3cooo     = 0
	  igcm_ch3cooh     = 0
	  igcm_ch3coooh    = 0 
	  igcm_c3h8        = 0 
	  igcm_ic3h7o2     = 0
	  igcm_ic3h7ooh    = 0
	  igcm_ch3coch3    = 0
	  igcm_ch3coch2o2  = 0
      igcm_hoch2o2	   = 0
      igcm_hoch2ooh	   = 0
      igcm_hoch2oh	   = 0
      igcm_c2h5o       = 0
      igcm_c2h5oh      = 0
      igcm_ch3co       = 0
      igcm_ch3choho2   = 0 
      igcm_ch2oo_e     = 0
      igcm_ch2oo	   = 0
	  igcm_hco  	   = 0
      igcm_c2h5        = 0
      igcm_n2          = 0
      igcm_ar          = 0
      igcm_ar_n2       = 0
      igcm_co2plus     = 0
      igcm_oplus       = 0
      igcm_o2plus      = 0
      igcm_coplus      = 0
      igcm_cplus       = 0
      igcm_nplus       = 0
      igcm_noplus      = 0
      igcm_n2plus      = 0
      igcm_hplus       = 0
      igcm_hco2plus    = 0
      igcm_elec        = 0

!   Chlorine species 
     igcm_cl = 0
     igcm_clo = 0
     igcm_cl2 = 0
     igcm_oclo = 0
     igcm_cl2o2 = 0
     igcm_hcl = 0
     igcm_hocl = 0
     igcm_cloo = 0
     igcm_ch3ocl = 0
     igcm_clco = 0


! 1.2 find dust tracers

      count = 0

      if (dustbin > 0) then
         do iq = 1,nqmx
            txt = " "
            write(txt,'(a4,i2.2)') 'dust', count + 1
            if (noms(iq) == txt) then
               count = count + 1
               igcm_dustbin(count) = iq
               mmol(iq) = 100.
            end if
         end do !do iq=1,nqmx
      end if ! of if (dustbin.gt.0)

      if (doubleq) then
         do iq = 1,nqmx
            if (noms(iq) == "dust_mass") then
               igcm_dust_mass = iq
               count = count + 1
            end if
            if (noms(iq) == "dust_number") then
               igcm_dust_number = iq
               count = count + 1
            end if
         end do
      end if ! of if (doubleq)
      if (scavenging) then
         do iq = 1,nqmx
            if (noms(iq) == "ccn_mass") then
               igcm_ccn_mass = iq
               count = count + 1
            end if
            if (noms(iq) == "ccn_number") then
               igcm_ccn_number = iq
               count = count + 1
            end if
         end do
      end if ! of if (scavenging)

      if (submicron) then
         do iq=1,nqmx
            if (noms(iq) == "dust_submicron") then
               igcm_dust_submicron = iq
               mmol(iq) = 100.
               count = count + 1
            end if
         end do
      end if ! of if (submicron)

! 1.3 find chemistry and water tracers

      print*,noms

      
      nbqchem = 0
      do iq = 1,nqmx
         if (noms(iq) == "co2") then
            igcm_co2 = iq
            mmol(igcm_co2) = 44.
            count = count + 1
            nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "co") then
           igcm_co = iq
           mmol(igcm_co) = 28.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "o") then
           igcm_o = iq
           mmol(igcm_o) = 16.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "o1d") then
           igcm_o1d = iq
           mmol(igcm_o1d) = 16.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "o2") then
           igcm_o2 = iq
           mmol(igcm_o2) = 32.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "o3") then
           igcm_o3 = iq
           mmol(igcm_o3) = 48.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "h") then
           igcm_h = iq
           mmol(igcm_h) = 1.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "h2") then
           igcm_h2 = iq
           mmol(igcm_h2) = 2.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "oh") then
           igcm_oh = iq
           mmol(igcm_oh) = 17.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ho2") then
           igcm_ho2 = iq
           mmol(igcm_ho2) = 33.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "h2o2") then
           igcm_h2o2 = iq
           mmol(igcm_h2o2) = 34.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch4") then
           igcm_ch4 = iq
           mmol(igcm_ch4) = 16.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch3") then
           igcm_ch3 = iq
           mmol(igcm_ch3) = 15.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch3o2") then
           igcm_ch3o2 = iq
           mmol(igcm_ch3o2) = 47.03
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch3ooh") then
           igcm_ch3ooh = iq
           mmol(igcm_ch3ooh) = 48.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch3oh") then
           igcm_ch3oh = iq
           mmol(igcm_ch3oh) = 32.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch3o") then
           igcm_ch3o = iq
           mmol(igcm_ch3o) = 31.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "hcho") then
           igcm_hcho = iq
           mmol(igcm_hcho) = 30.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "hcooh") then
           igcm_hcooh = iq
           mmol(igcm_hcooh) = 46.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "hoch2o2") then
           igcm_hoch2o2 = iq
           mmol(igcm_hoch2o2) = 63.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hoch2ooh") then
           igcm_hoch2ooh = iq
           mmol(igcm_hoch2ooh) = 64.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hoch2oh") then
           igcm_hoch2oh = iq
           mmol(igcm_hoch2oh) = 48.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hco") then
           igcm_hco = iq
           mmol(igcm_hco) = 29.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "c2h6") then
           igcm_c2h6 = iq
           mmol(igcm_c2h6) = 30.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "c2h5") then
           igcm_c2h5 = iq
           mmol(igcm_c2h5) = 29.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "c2h5o2") then
           igcm_c2h5o2 = iq
           mmol(igcm_c2h5o2) = 61.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "c2h5ooh") then
           igcm_c2h5ooh = iq
           mmol(igcm_c2h5ooh) = 62.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "c2h5oh") then
           igcm_c2h5oh = iq
           mmol(igcm_c2h5oh) = 46.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "hoch2ch2o2") then
           igcm_hoch2ch2o2 = iq
           mmol(igcm_hoch2ch2o2) = 77.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "hoch2ch2o") then
           igcm_hoch2ch2o = iq
           mmol(igcm_hoch2ch2o) = 61.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ethgly") then
           igcm_ethgly = iq
           mmol(igcm_ethgly) = 62.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hyetho2h") then
           igcm_hyetho2h = iq
           mmol(igcm_hyetho2h) = 78.
           count = count + 1
           nbqchem = nbqchem + 1
        end if		
        if (noms(iq) == "ch3cho") then
           igcm_ch3cho = iq
           mmol(igcm_ch3cho) = 44.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch2choh") then
           igcm_ch2choh = iq
           mmol(igcm_ch2choh) = 44.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch3choho2") then
           igcm_ch3choho2 = iq
           mmol(igcm_ch3choho2) = 77.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch3cooh") then
           igcm_ch3cooh = iq
           mmol(igcm_ch3cooh) = 60.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch3chohooh") then
           igcm_ch3chohooh = iq
           mmol(igcm_ch3chohooh) = 78.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "ch3c(o)") then
           igcm_ch3co = iq
           mmol(igcm_ch3co) = 43.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch3c(o)oo") then
           igcm_ch3cooo = iq
           mmol(igcm_ch3cooo) = 75.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ch3c(o)ooh") then
           igcm_ch3coooh = iq
           mmol(igcm_ch3coooh) = 76.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "hcoch2o2") then
           igcm_hcoch2o2 = iq
           mmol(igcm_hcoch2o2) = 75.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "glyox") then
           igcm_glyox = iq
           mmol(igcm_glyox) = 58.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "hcoco") then
           igcm_hcoco = iq
           mmol(igcm_hcoco) = 57.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hooch2cho") then
           igcm_hooch2cho = iq
           mmol(igcm_hooch2cho) = 76.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hoch2cho") then
           igcm_hoch2cho = iq
           mmol(igcm_hoch2cho) = 60.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "hochcho") then
           igcm_hochcho = iq
           mmol(igcm_hochcho) = 75.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hoch2co") then
           igcm_hoch2co = iq
           mmol(igcm_hoch2co) = 59.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "hoch2co3") then
           igcm_hoch2co3 = iq
           mmol(igcm_hoch2co3) = 91.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hoch2co2h") then
           igcm_hoch2co2h = iq
           mmol(igcm_hoch2co2h) = 76.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hcoco2h") then
           igcm_hcoco2h = iq
           mmol(igcm_hcoco2h) = 74.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hoch2co3h") then
           igcm_hoch2co3h = iq
           mmol(igcm_hoch2co3h) = 92.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hcoco3h") then
           igcm_hcoco3h = iq
           mmol(igcm_hcoco3h) = 90.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
        if (noms(iq) == "hcoco3") then
           igcm_hcoco3 = iq
           mmol(igcm_hcoco3) = 89.
           count = count + 1
           nbqchem = nbqchem + 1
        end if	
				
        if (noms(iq) == "n2") then
           igcm_n2 = iq
           mmol(igcm_n2) = 28.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "n") then
           igcm_n = iq
           mmol(igcm_n) = 14.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "n2d") then
           igcm_n2d = iq
           mmol(igcm_n2d) = 14.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "no") then
           igcm_no = iq
           mmol(igcm_no) = 30.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "no2") then
           igcm_no2 = iq
           mmol(igcm_no2) = 46.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "ar") then
           igcm_ar = iq
           mmol(igcm_ar) = 40.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "h2o_vap") then
           igcm_h2o_vap = iq
           mmol(igcm_h2o_vap) = 18.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (trim(noms(iq)) == "h2o_ice") then
           igcm_h2o_ice = iq
           mmol(igcm_h2o_ice) = 18.
           count = count + 1
           nbqchem = nbqchem + 1
        end if


!     Chlorine Scheme 22/04/2020 
        if (noms(iq).eq."cl") then
          igcm_cl=iq
          mmol(igcm_cl)=35.45
          count=count+1
          nbqchem = nbqchem + 1
        endif
        if (noms(iq).eq."clo") then
          igcm_clo=iq
          mmol(igcm_clo)=51.54
          count=count+1
          nbqchem = nbqchem + 1
        endif
        if (noms(iq).eq."cl2") then
          igcm_cl2=iq
          mmol(igcm_cl2)=70.9
          count=count+1
          nbqchem = nbqchem + 1
        endif
        if (noms(iq).eq."oclo") then
          igcm_oclo=iq
          mmol(igcm_oclo)=67.45
          count=count+1
          nbqchem = nbqchem + 1
        endif
        if (noms(iq).eq."cl2o2") then
          igcm_cl2o2=iq
          mmol(igcm_cl2o2)=102.9
          count=count+1
          nbqchem = nbqchem + 1
        endif
        if (noms(iq).eq."hcl") then
          igcm_hcl=iq
          mmol(igcm_hcl)=36.46
          count=count+1
          nbqchem = nbqchem + 1
        endif
        if (noms(iq).eq."hocl") then
          igcm_hocl=iq
          mmol(igcm_hocl)=52.46
          count=count+1
          nbqchem = nbqchem + 1
        endif
        if (noms(iq).eq."cloo") then
          igcm_cloo=iq
          mmol(igcm_cloo)=67.45
          count=count+1
          nbqchem = nbqchem + 1
        endif
        if (noms(iq).eq."ch3ocl") then
          igcm_ch3ocl=iq
          mmol(igcm_ch3ocl)=66.49
          count=count+1
          nbqchem = nbqchem + 1
        endif
        if (noms(iq).eq."clco") then
          igcm_clco=iq
          mmol(igcm_clco)=63.46
          count=count+1
          nbqchem = nbqchem + 1
        endif
        
        
        
        
        
        
        
! 1.4 find ions

        if (noms(iq) == "co2plus") then
           igcm_co2plus = iq
           mmol(igcm_co2plus) = 44.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "oplus") then
           igcm_oplus = iq
           mmol(igcm_oplus) = 16.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "o2plus") then
           igcm_o2plus = iq
           mmol(igcm_o2plus) = 32.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "coplus") then
           igcm_coplus = iq
           mmol(igcm_coplus) = 28.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "cplus") then
           igcm_cplus = iq
           mmol(igcm_cplus) = 12.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "nplus") then
           igcm_nplus = iq
           mmol(igcm_nplus) = 14.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "noplus") then
           igcm_noplus = iq
           mmol(igcm_noplus) = 30.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "n2plus") then
           igcm_n2plus = iq
           mmol(igcm_n2plus) = 28.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "hplus") then
           igcm_hplus = iq
           mmol(igcm_hplus) = 1.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "hco2plus") then
           igcm_hco2plus = iq
           mmol(igcm_hco2plus) = 45.
           count = count + 1
           nbqchem = nbqchem + 1
        end if
        if (noms(iq) == "elec") then
           igcm_elec = iq
           mmol(igcm_elec) = 1./1822.89
           count = count + 1
        end if

! 1.5 find idealized non-condensible tracer

        if (noms(iq) == "Ar_N2") then
           igcm_ar_n2 = iq
           mmol(igcm_ar_n2) = 30.
           count = count + 1
        end if
		
      end do ! of do iq=1,nqmx

	  
! 1.6 check that we identified all tracers:

      if (count /= nqmx) then
         write(*,*) "inichim_newstart: found only ",count," tracers"
         write(*,*) "                  expected ",nqmx
         do iq = 1,count
            write(*,*) '      ', iq, ' ', trim(noms(iq))
         end do
         stop
      else
         write(*,*) "inichim_newstart: found all expected tracers"
         do iq = 1,nqmx
            write(*,*) '      ', iq, ' ', trim(noms(iq))
         end do
      end if

! 1.7 check if nitrogen species are present:

      if(igcm_no == 0) then
         !check that no N species is in traceur.def
         if(igcm_n /= 0 .or. igcm_no2 /= 0 .or. igcm_n2d /= 0) then
            write(*,*)'inichim_newstart error:'
            write(*,*)'N, NO2 and/or N2D are in traceur.def, but not NO'
            write(*,*)'stop'
            stop
         endif
         flagnitro = .false.
         nspe = 14
      else
         !check that all N species are in traceur.def
         if(igcm_n == 0 .or. igcm_no2 == 0 .or. igcm_n2d == 0) then
            write(*,*)'inichim_newstart error:'
            write(*,*)'if NO is in traceur.def, N, NO2 and N2D must also be'
            write(*,*)'stop'
            stop
         endif
         flagnitro = .true.
         nspe = 18
      endif

! 1.8 allocate arrays

      allocate(niq(nspe))
      allocate(vmrinit(nalt,nspe))
      allocate(vmrint(nspe))

! 2. load in chemistry data for initialization:

! order of major species in initialization file:
!
!    1: co2 
!    2: ar
!    3: n2  
!    4: o2  
!    5: co  
!    6: o   
!    7: h2
!
! order of minor species in initialization file:
!
!    1: h  
!    2: oh 
!    3: ho2 
!    4: h2o 
!    5: h2o2 
!    6: o1d
!    7: o3
!
! order of nitrogen species in initialization file:
!
!    1: n
!    2: no
!    3: no2
!    4: n2d

! major species:

      niq(1) = igcm_co2
      niq(2) = igcm_ar
      niq(3) = igcm_n2
      niq(4) = igcm_o2
      niq(5) = igcm_co
      niq(6) = igcm_o
      niq(7) = igcm_h2

! minor species:

      niq(8)  = igcm_h
      niq(9)  = igcm_oh
      niq(10) = igcm_ho2
      niq(11) = igcm_h2o_vap
      niq(12) = igcm_h2o2
      niq(13) = igcm_o1d 
      niq(14) = igcm_o3

! nitrogen species:

      if (flagnitro) then
         niq(15) = igcm_n
         niq(16) = igcm_no
         niq(17) = igcm_no2
         niq(18) = igcm_n2d         
      end if

! 2.1 open initialization files

      open(210, iostat=ierr,file=trim(datafile)//'/atmosfera_LMD_may.dat')
      if (ierr /= 0) then
         write(*,*)'Error : cannot open file atmosfera_LMD_may.dat '
         write(*,*)'(in aeronomars/inichim_newstart.F)'
         write(*,*)'It should be in :', trim(datafile),'/'
         write(*,*)'1) You can change this path in callphys.def with'
         write(*,*)'   datadir=/path/to/datafiles/'
         write(*,*)'2) If necessary atmosfera_LMD_may.dat (and others)'
         write(*,*)'   can be obtained online on:'
         write(*,*)' http://www.lmd.jussieu.fr/~forget/datagcm/datafile'
         stop
      end if
      open(220, iostat=ierr,file=trim(datafile)//'/atmosfera_LMD_min.dat')
      if (ierr /= 0) then
         write(*,*)'Error : cannot open file atmosfera_LMD_min.dat '
         write(*,*)'(in aeronomars/inichim_newstart.F)'
         write(*,*)'It should be in :', trim(datafile),'/'
         write(*,*)'1) You can change this path in callphys.def with'
         write(*,*)'   datadir=/path/to/datafiles/'
         write(*,*)'2) If necessary atmosfera_LMD_min.dat (and others)'
         write(*,*)'   can be obtained online on:'
         write(*,*)' http://www.lmd.jussieu.fr/~forget/datagcm/datafile'
         stop
      end if
      if(flagnitro) then
         open(230, iostat=ierr,file=trim(datafile)//'/atmosfera_LMD_nitr.dat')
         if (ierr.ne.0) then
            write(*,*)'Error : cannot open file atmosfera_LMD_nitr.dat '
            write(*,*)'(in aeronomars/inichim_newstart.F)'
            write(*,*)'It should be in :', datafile
            write(*,*)'1) You can change this directory address in '
            write(*,*)'   file phymars/datafile.h'
            write(*,*)'2) If necessary atmosfera_LMD_nitr.dat (and others)'
            write(*,*)'   can be obtained online on:'
            write(*,*)' http://www.lmd.jussieu.fr/~forget/datagcm/datafile'
            STOP
         endif
      endif   ! Of if(flagnitro)

! 2.2 read initialization files

! major species

      read(210,*)
      do l = 1,nalt
         read(210,*) dummy, tinit(l), pinit(l), densinit(l), &
                     (vmrinit(l,n), n = 1,7)
         pinit(l) = pinit(l)*100.              ! conversion in Pa
         pinit(l) = log(pinit(l))              ! for the vertical interpolation
      end do
      close(210)

! minor species

      read(220,*)
      do l = 1,nalt
         read(220,*) real_dummy, (vmrinit(l,n), n = 8,14)
      end do 
      close(220)

! nitrogen species

      if (flagnitro) then
         read(230,*)
         do l = 1,nalt
            read(230,*) dummy, (vmrinit(l,n), n = 15,18)
         end do
         close(230)
      end if
      
! 3. initialization of tracers

      do i = 1,iip1
         do j = 1,jjp1
            do l = 1,llm

               pgcm = aps(l) + bps(l)*ps(i,j)  ! gcm pressure
               pgcm = log(pgcm)                ! for the vertical interpolation
               mmean(i,j,l) = 0.

! 3.1 vertical interpolation



               do n = 1,nspe
          
                       call intrplf(pgcm,vmr,pinit,vmrinit(:,n),nalt)
                 
                      vmrint(n) = vmr
                      iq = niq(n)
                      mmean(i,j,l) = mmean(i,j,l) + vmrint(n)*mmol(iq)
               end do


! 3.2 attribute mixing ratio: - all layers or only thermosphere
!                             - with our without h2o 

               if (flagthermo == 0 .or. (flagthermo == 1 .and. exp(pgcm) < 0.1)) then
                  do n = 1,nspe
                     iq = niq(n)
                     if (iq /= igcm_h2o_vap .or. flagh2o == 1) then
                        pq(i,j,l,iq) = vmrint(n)*mmol(iq)/mmean(i,j,l)
                     end if
                  end do
               end if

            end do
         end do
      end do

! set surface values of chemistry tracers to zero

      if (flagthermo == 0) then
         ! NB: no problem for "surface water vapour" tracer which is always 0
         do n = 1,nspe
            iq = niq(n)
            qsurf(1:ngridmx,iq) = 0.
         end do
      end if

! 3.3 initialization of tracers not contained in the initialization files

! methane : 0 ppbv

      if (igcm_ch4 /= 0) then


        ! Polynomial line of best fit of Webster 2018
         ! vmr = (-8.1565118E-14*zls**6 &
          ! +   8.3522676E-11*zls**5 &
          ! -   3.0856216E-08*zls**4 &
          ! +   4.8741198E-06*zls**3 &
          ! -   2.9426485E-04*zls**2 &
          ! +   4.6374896E-03*zls &
          ! +   2.8769755E-01)*1.e-9
!ch4vmr = 410.e-12
  
ch4vmr=410.e-12

          do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch4) = ch4vmr*mmol(igcm_ch4)/mmean(i,j,l)
                  ! if ( l .le. 7 ) pq(i,j,l,igcm_ch4) = 10.e-9*mmol(igcm_ch4)/mmean(i,j,l)
               end do
            end do
         end do
         
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch4) = 0.
      end if

! 13CH4
      if (igcm_13ch4 /= 0) then
         vmr =  0.!14.e-16
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_13ch4) = vmr*mmol(igcm_13ch4)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_13ch4) = 0.
      end if
	  
! C3H8
      if (igcm_c3h8 /= 0) then
         vmr = 0.!2.25e-9       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_c3h8) = vmr*mmol(igcm_c3h8)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_c3h8) = 0.
      end if  	
	  
! c2h6 : 0 ppbv

      if (igcm_c2h6 /= 0) then
         vmr =  0.!14.e-9
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_c2h6) = vmr*mmol(igcm_c2h6)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_c2h6) = 0.
      end if

! c2h4 : 0 ppbv

      if (igcm_c2h4 /= 0) then
         vmr = 0.!.e-12 
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_c2h4) = vmr*mmol(igcm_c2h4)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_c2h4) = 0.
      end if

! c2h2 : 0 ppbv

      if (igcm_c2h2 /= 0) then
         vmr = 0.!1.e-12
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_c2h2) = vmr*mmol(igcm_c2h2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_c2h2) = 0.
      end if

! ch3o2 : 0 ppbv

      if (igcm_ch3o2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3o2) = vmr*mmol(igcm_ch3o2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3o2) = 0.
      end if

! ch3 : 0 ppbv

      if (igcm_ch3 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3) = vmr*mmol(igcm_ch3)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3) = 0.
      end if

! ch3oh : 0 ppbv

      if (igcm_ch3oh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3oh) = vmr*mmol(igcm_ch3oh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3oh) = 0.
      end if

! formaldehyde (hcho) : 0 ppbv

      if (igcm_hcho /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hcho) = vmr*mmol(igcm_hcho)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hcho) = 0.
      end if

! CH3OOH
      if (igcm_ch3ooh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3ooh) = vmr*mmol(igcm_ch3ooh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3ooh) = 0.
      end if   

! C2H5O2
      if (igcm_c2h5o2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_c2h5o2) = vmr*mmol(igcm_c2h5o2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_c2h5o2) = 0.
      end if   
! 1CH2
      if (igcm_1ch2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_1ch2) = vmr*mmol(igcm_1ch2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_1ch2) = 0.
      end if  

! 3CH2
      if (igcm_3ch2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_3ch2) = vmr*mmol(igcm_3ch2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_3ch2) = 0.
      end if 

! CH2	  
      if (igcm_ch2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch2) = vmr*mmol(igcm_ch2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch2) = 0.
      end if  


! CH
      if (igcm_ch /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch) = vmr*mmol(igcm_ch)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch) = 0.
      end if  

! CH3O
      if (igcm_ch3o /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3o) = vmr*mmol(igcm_ch3o)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3o) = 0.
      end if  

	  
! HCOOH
      if (igcm_hcooh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hcooh) = vmr*mmol(igcm_hcooh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hcooh) = 0.
      end if  

! C2H5OOH
      if (igcm_c2h5ooh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_c2h5ooh) = vmr*mmol(igcm_c2h5ooh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_c2h5ooh) = 0.
      end if  
	  
! CH3CHO
      if (igcm_ch3cho /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3cho) = vmr*mmol(igcm_ch3cho)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3cho) = 0.
      end if  

! CH3COOO
      if (igcm_ch3cooo /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3cooo) = vmr*mmol(igcm_ch3cooo)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3cooo) = 0.
      end if  
	  
! CH3COOH
      if (igcm_ch3cooh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3cooh) = vmr*mmol(igcm_ch3cooh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3cooh) = 0.
      end if  
	
! CH3C(O)OOH
      if (igcm_ch3coooh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3coooh) = vmr*mmol(igcm_ch3coooh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3coooh) = 0.
      end if  	
	  
! HCO
      if (igcm_hco /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hco) = vmr*mmol(igcm_hco)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hco) = 0.
      end if  	
	  
! C2H5
      if (igcm_c2h5 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_c2h5) = vmr*mmol(igcm_c2h5)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_c2h5) = 0.
      end if  	
	  
! HOCH2O2
      if (igcm_hoch2o2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hoch2o2) = vmr*mmol(igcm_hoch2o2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hoch2o2) = 0.
      end if 
	  
! HOCH2OOH
      if (igcm_hoch2ooh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hoch2ooh) = vmr*mmol(igcm_hoch2ooh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hoch2ooh) = 0.
      end if 
	  
! HOCH2OH
      if (igcm_hoch2oh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hoch2oh) = vmr*mmol(igcm_hoch2oh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hoch2oh) = 0.
      end if 
	  
! C2H5O
      if (igcm_c2h5o /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_c2h5o) = vmr*mmol(igcm_c2h5o)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_c2h5o) = 0.
      end if 
	  
! C2H5OH
      if (igcm_c2h5oh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_c2h5oh) = vmr*mmol(igcm_c2h5oh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_c2h5oh) = 0.
      end if 
	  
! CH3CO
      if (igcm_ch3co /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3co) = vmr*mmol(igcm_ch3co)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3co) = 0.
      end if 
	  
! CH3CHOHO2
      if (igcm_ch3choho2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3choho2) = vmr*mmol(igcm_ch3choho2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3choho2) = 0.
      end if  
	  
! HOCH2CH2O2
      if (igcm_hoch2ch2o2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hoch2ch2o2) = vmr*mmol(igcm_hoch2ch2o2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hoch2ch2o2) = 0.
      end if  
	  
! HOCH2CH2O
      if (igcm_hoch2ch2o /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hoch2ch2o) = vmr*mmol(igcm_hoch2ch2o2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hoch2ch2o) = 0.
      end if  
	  
! ETHGLY
      if (igcm_ethgly /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ethgly) = vmr*mmol(igcm_ethgly)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ethgly) = 0.
      end if  
	  
! HYETHO2H
      if (igcm_hyetho2h /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hyetho2h) = vmr*mmol(igcm_hyetho2h)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hyetho2h) = 0.
      end if  
	  
! CH3CHOHOOH
      if (igcm_ch3chohooh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3chohooh) = vmr*mmol(igcm_ch3chohooh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3chohooh) = 0.
      end if  
	  
! HCOCH2O2
      if (igcm_hcoch2o2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hcoch2o2) = vmr*mmol(igcm_hcoch2o2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hcoch2o2) = 0.
      end if  
	  
! GLYOX
      if (igcm_glyox /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_glyox) = vmr*mmol(igcm_glyox)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_glyox) = 0.
      end if  
	  
! HCOCO
      if (igcm_hcoco /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hcoco) = vmr*mmol(igcm_hcoco)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hcoco) = 0.
      end if  
	 
! HOOCH2CHO
      if (igcm_hooch2cho /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hooch2cho) = vmr*mmol(igcm_hooch2cho)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hooch2cho) = 0.
      end if  
	  
! HOCH2CHO
      if (igcm_hoch2cho /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hoch2cho) = vmr*mmol(igcm_hoch2cho)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hoch2cho) = 0.
      end if  
	  
! HOCHCHO
      if (igcm_hochcho /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hochcho) = vmr*mmol(igcm_hochcho)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hochcho) = 0.
      end if 

! HOCH2CO
      if (igcm_hoch2co /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hoch2co) = vmr*mmol(igcm_hoch2co)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hoch2co) = 0.
      end if  
	  
! HOCH2CO3
      if (igcm_hoch2co3 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hoch2co3) = vmr*mmol(igcm_hoch2co3)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hoch2co3) = 0.
      end if  
	  
! HOCH2CO2H
      if (igcm_hoch2co2h /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hoch2co2h) = vmr*mmol(igcm_hoch2co2h)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hoch2co2h) = 0.
      end if  
	  
! HCOCO2H
      if (igcm_hcoco2h /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hcoco2h) = vmr*mmol(igcm_hcoco2h)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hcoco2h) = 0.
      end if  

! HOCH2CO3H
      if (igcm_hoch2co3h /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hoch2co3h) = vmr*mmol(igcm_hoch2co3h)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hoch2co3h) = 0.
      end if 
	  
! HCOCO3H
      if (igcm_hcoco3h /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hcoco3h) = vmr*mmol(igcm_hcoco3h)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hcoco3h) = 0.
      end if 
	  
! HCOCO3 
      if (igcm_hcoco3 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hcoco3) = vmr*mmol(igcm_hcoco3)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hcoco3) = 0.
      end if 
	 
! Chlrorine species 22/04/2020
      ! Cl
      if (igcm_cl /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_cl) = vmr*mmol(igcm_cl)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_cl) = 0.
      end if 
      ! ClO
      if (igcm_clo /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_clo) = vmr*mmol(igcm_clo)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_clo) = 0.
      end if 
      ! OClO
      if (igcm_oclo /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_oclo) = vmr*mmol(igcm_oclo)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_oclo) = 0.
      end if 
      ! Cl2
      if (igcm_cl2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_cl2) = vmr*mmol(igcm_cl2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_cl2) = 0.
      end if 
      ! ClO2
      if (igcm_cl2o2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_cl2o2) = vmr*mmol(igcm_cl2o2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_cl2o2) = 0.
      end if 
      ! HCl
      if (igcm_hcl /= 0) then
         vmr = 0.
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hcl) = vmr*mmol(igcm_hcl)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hcl) = 0.
      end if 
      ! HOCl
      if (igcm_hocl /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_hocl) = vmr*mmol(igcm_hocl)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_hocl) = 0.
      end if 
      ! ClOO
      if (igcm_cloo /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_cloo) = vmr*mmol(igcm_cloo)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_cloo) = 0.
      end if 
      ! CH3COCl
      if (igcm_ch3ocl /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3ocl) = vmr*mmol(igcm_ch3ocl)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3ocl) = 0.
      end if 
      ! ClCO
      if (igcm_clco /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_clco) = vmr*mmol(igcm_clco)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_clco) = 0.
      end if 
      
      
      
      
      
      
      
	  
! CH2OO*
      if (igcm_ch2oo_e /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch2oo_e) = vmr*mmol(igcm_ch2oo_e)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch2oo_e) = 0.
      end if  
! CH2OO
      if (igcm_ch2oo /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch2oo) = vmr*mmol(igcm_ch2oo)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch2oo) = 0.
      end if  
	  
! i-C3H7O2
      if (igcm_ic3h7o2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ic3h7o2) = vmr*mmol(igcm_ic3h7o2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ic3h7o2) = 0.
      end if 
	  
! i-C3H7OOH
      if (igcm_ic3h7ooh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ic3h7ooh) = vmr*mmol(igcm_ic3h7ooh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ic3h7ooh) = 0.
      end if 
	  
! CH3COCH3
      if (igcm_ch3coch3 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3coch3) = vmr*mmol(igcm_ch3coch3)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3coch3) = 0.
      end if 

! CH3COCH2O2
      if (igcm_ch3coch2o2 /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch3coch2o2) = vmr*mmol(igcm_ch3coch2o2)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch3coch2o2) = 0.
      end if 
	  
! CH2CHOH
      if (igcm_ch2choh /= 0) then
         vmr = 0.       
         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  pq(i,j,l,igcm_ch2choh) = vmr*mmol(igcm_ch2choh)/mmean(i,j,l)
               end do
            end do
         end do
         ! set surface value to zero
         qsurf(1:ngridmx,igcm_ch2choh) = 0.
      end if 
	  
! ions: 0

      if (igcm_co2plus /= 0) then
         !check that all required ions are in traceur.def
         if (igcm_o2plus == 0 .or. igcm_oplus == 0 .or. igcm_coplus == 0          &
              .or. igcm_cplus == 0 .or. igcm_nplus == 0 .or. igcm_noplus == 0    & 
              .or. igcm_n2plus == 0 .or. igcm_hplus == 0 .or. igcm_hco2plus == 0) then 
            write(*,*)'inichim_newstart error:'
            write(*,*)'if co2plus is in traceur.def, all other ions must also be'
            write(*,*)'o2plus, oplus, coplus, cplus, nplus, noplus, n2plus'
            write(*,*)'hplus, hco2plus and elec'
            write(*,*)'stop'
            stop
         end if

         do i = 1,iip1
            do j = 1,jjp1
               do l = 1,llm
                  ! all ions to 0     
                  pq(i,j,l,igcm_co2plus)  = 0.
                  pq(i,j,l,igcm_o2plus)   = 0.
                  pq(i,j,l,igcm_oplus)    = 0.
                  pq(i,j,l,igcm_coplus)   = 0.
                  pq(i,j,l,igcm_cplus)    = 0.
                  pq(i,j,l,igcm_nplus)    = 0.
                  pq(i,j,l,igcm_noplus)   = 0.
                  pq(i,j,l,igcm_n2plus)   = 0.
                  pq(i,j,l,igcm_hplus)    = 0.
                  pq(i,j,l,igcm_hco2plus) = 0.
                  pq(i,j,l,igcm_elec)     = 0.
               end do
            end do
         end do

         ! surface value to 0

         qsurf(1:ngridmx,igcm_co2plus)  = 0.
         qsurf(1:ngridmx,igcm_o2plus)   = 0.
         qsurf(1:ngridmx,igcm_oplus)    = 0.
         qsurf(1:ngridmx,igcm_coplus)   = 0.
         qsurf(1:ngridmx,igcm_cplus)    = 0.
         qsurf(1:ngridmx,igcm_nplus)    = 0.
         qsurf(1:ngridmx,igcm_noplus)   = 0.
         qsurf(1:ngridmx,igcm_n2plus)   = 0.
         qsurf(1:ngridmx,igcm_hplus)    = 0.
         qsurf(1:ngridmx,igcm_hco2plus) = 0.
         qsurf(1:ngridmx,igcm_elec)     = 0.

      else

         if (igcm_o2plus /= 0 .or. igcm_oplus /= 0 .or. igcm_coplus /= 0          &
              .or. igcm_cplus /= 0 .or. igcm_nplus /= 0 .or. igcm_noplus /= 0    & 
              .or. igcm_n2plus /= 0 .or. igcm_hplus /= 0 .or. igcm_hco2plus /= 0) then 
            write(*,*)'inichim_newstart error:'
            write(*,*)'some ions are in traceur.def, but not co2plus'
            write(*,*)'stop'
            stop
         end if
      end if    ! of if(igcm_co2 /= 0)
      
      ! deallocations

      deallocate(niq)
      deallocate(vmrinit)
      deallocate(vmrint)


      end
