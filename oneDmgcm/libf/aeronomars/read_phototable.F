      subroutine read_phototable
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   subject:
c   --------
c
c   read photolysis lookup table
c
c   VERSION: 17/03/2011
c
c   Author:   Franck Lefevre
c
c   Arguments:
c   ----------
c
c   The output variable is jphot and is put in common chimiedata.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
! to use  'getin'
      USE netcdf
      USE ioipsl_getincom 

      implicit none
c
#include "netcdf.inc"
#include "dimensions.h"
#include "chimiedata.h"
#include "datafile.h"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer fic, ij, iozo, isza, itemp, iz, itau, ierr
c
      real xsza
      character(len=128) :: phototable ! photolysis table file name
	  
ccc Organic netCDF File 
      character(len=128) :: netcdf_phototable ! photolysis for new organics
      integer ncid    ! file ID
      integer retval  ! Error code
      real org_j(ntemp, nsza, 201, nozo, ntau)
      integer orgvarid 
      integer numorgs 
      parameter(numorgs = 13)
      character(len=30) organic_reactions(numorgs)
ccc Chlorine compounds inside organic netCDF
      real chlor_j(ntemp,nsza, 201, nozo,ntau)
      integer clvarid 
      integer numchlor
      parameter (numchlor = 4)
      character(len=30) chlor_reactions(numchlor)

      integer i 
c
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
! set photolysis table input file name
      phototable="jmars.20111014" ! default
      ! look for a " phototable= ..." option in def files
      call getin("phototable",phototable)

      fic = 81

      open(fic, form = 'formatted', status = 'old',
     &file =trim(datafile)//"/"//trim(phototable),iostat=ierr)

      IF (ierr.NE.0) THEN
        write(*,*)'Error : cannot open photolysis lookup table ',
     &             trim(phototable)
        write(*,*)'(in aeronomars/init_chemistry.F)'
        write(*,*)'It should be in :',trim(datafile),'/'
        write(*,*)'1) You can change this directory path in physiq.def'
        write(*,*)'   with:'
        write(*,*)'   datadir=/path/to/the/directory'
        write(*,*)'2) You can change the input phototable file name in'
        write(*,*)'   physiq.def with:'
        write(*,*)'   phototable=filename'
        STOP
      ENDIF

c
      print*, 'read photolysis lookup table ',trim(phototable)
c


      do itau = 1,ntau
         do itemp = 1,ntemp
            do iozo = 1,nozo
               do isza = 1,nsza
                  do iz = 200,0,-1
                     read(fic,*) colairtab(iz), xsza, table_ozo(iozo)
c                    print*, colairtab(iz), xsza, table_ozo(iozo)
                     read(fic,'(7e11.4)')
     $               (jphot(itemp,isza,iz,iozo,itau,ij),
     $                                       ij= 1, nd-numorgs-numchlor)
                     do ij = 1,nd-numorgs-numchlor
                        if (jphot(itemp,isza,iz,iozo,itau,ij)
     $                     .eq. 1.e-30) then
                           jphot(itemp,isza,iz,iozo,itau,ij) = 0.
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      print*, 'lookup table...ok'
      close(fic)
	  
       
       
       
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  New organic rate coefficients read in from a 
c  netCDF file located in datafile      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! NetCDF file name 
C        netcdf_phototable = "/jmars_organics.nc"
       netcdf_phototable = "/jmars_organics_cl.nc"

       write(*,*) "NEW ORGANIC PHOTOLYSIS:"
       WRITE(*,*) trim(datafile) // trim(netcdf_phototable)
       retval = nf_open(trim(datafile) // trim(netcdf_phototable), 
     $                        NF_NOWRITE, ncid)
     
       if (retval .ne. nf_noerr) call net_err(retval)


C Create array with reaction name strings and iterate through to 
c save file size
       organic_reactions(1) = "CH3CHO -> CH3 + HCO"
       organic_reactions(2) = "CH3CHO -> CH4 + CO"   
       organic_reactions(3) = "CH3CHO -> CH3CO + H"
       organic_reactions(4) = "HOCH2OOH -> HOCH2O. + OH"
       organic_reactions(5) = "HOCH2CHO -> CH2OH + HCO"
       organic_reactions(6) = "HOCH2CHO -> CH3OH + CO"
       organic_reactions(7) = "HOCH2CHO -> CH2CHO + OH"
       organic_reactions(8) = "CHOCHO -> HCO + HCO"
       organic_reactions(9) = "CHOCHO -> H2 + 2CO"
       organic_reactions(10) = "CHOCHO -> CH2O + CO"
       organic_reactions(11) = "CH3COOH -> CH3 + COOH"
       organic_reactions(12) = "CH3CO(OOH) -> Products"
       organic_reactions(13) = "CH3COCO(OH) -> Products"
       
c     Chlorine species 
       chlor_reactions(1) = "Cl2 -> Cl + Cl"
       chlor_reactions(2) = "ClOO -> Products"
       chlor_reactions(3) = "OClO -> Products"
       chlor_reactions(4) = "ClOOCl -> Cl + ClOO "




       ! Iterate through the reaction strings and store
       ! the values extracted from the netCDF table and
       ! place in the jphot array

       do i = 1, numorgs
       
       retval = nf_inq_varid(ncid,TRIM(organic_reactions(i)),orgvarid)
       if (retval .ne. nf_noerr) call net_err(retval)
       retval = nf_get_var_real(ncid, orgvarid, org_j)
       if (retval .ne. nf_noerr) call net_err(retval)       

       jphot(:,:,:,:,:,ij+i-1) = org_j
	   

       end do 


       do i = 1, numchlor
       
       retval = nf_inq_varid(ncid,TRIM(chlor_reactions(i)),clvarid)
       if (retval .ne. nf_noerr) call net_err(retval)
       retval = nf_get_var_real(ncid, clvarid, chlor_j)
       if (retval .ne. nf_noerr) call net_err(retval)       

       jphot(:,:,:,:,:,ij+numorgs + i - 1) = chlor_j
     

       end do 

       
       retval = nf_close(ncid)
       if (retval .ne. nf_noerr) call net_err(retval)       


      return
      end

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine net_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end