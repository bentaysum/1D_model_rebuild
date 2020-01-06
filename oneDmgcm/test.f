      program test

      implicit none

      character*100 :: infile

      infile = 'runpop'
      
      call checkfile(infile,99)
      
      stop
      end
      
