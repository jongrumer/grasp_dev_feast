      subroutine copy7t9
      integer fil_1,utfil
      parameter (fil_1=7, utfil=9)
      character rad11*1000
      open(fil_1,file='rcsf.out',status='unknown')
      open(unit=utfil,file='fil1.dat',status='unknown')
      do 
      read(utfil,999,end=100) rad11
      write(fil_1,999) trim(rad11)
      enddo
  100 close  (utfil,status='delete')
      close  (fil_1)
      return
  999 format(A)
      end
