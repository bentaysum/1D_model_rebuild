      subroutine checkfile(infile,iunit)

      implicit none
      
      character*100 :: infile
      character*100 :: errmsg

      integer       :: iunit, ierr
      
      
      open(iunit,file=trim(infile),
     &     status='old',form='formatted',
     &     iostat=ierr)
      close(99)

      write(errmsg,'(A)') trim(infile)//' does not exist; STOP'
      print*, errmsg
      stop
      
      end  
