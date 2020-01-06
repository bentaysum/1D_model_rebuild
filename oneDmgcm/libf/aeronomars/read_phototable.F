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
      USE ioipsl_getincom 

      implicit none
c
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
     $               (jphot(itemp,isza,iz,iozo,itau,ij), ij= 1,nd)
                     do ij = 1,nd
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
c
      return
      end
