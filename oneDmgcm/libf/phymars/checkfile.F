      subroutine checkfile(infile,iunit)

      implicit none
      
      character (len=100) :: infile
      character (len=100) :: errmsg

      integer       :: iunit, ierr
      
      
      open(iunit,file=trim(infile),
     &     status='old',form='formatted',
     &     iostat=ierr)
      close(99)

      if (ierr.ne.0) then
         write(errmsg,'(A)') trim(infile)//' does not exist; STOP'
         print*, errmsg
         stop
      endif
      
      end  
