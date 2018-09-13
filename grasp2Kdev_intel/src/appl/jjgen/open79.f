      subroutine open79(i)
      integer fil_1,utfil
      parameter (fil_1=7, utfil=9)
      i1=mod(i,2)
      if(i1.eq.0) then
      open(fil_1,file='fil1.dat',status='unknown')
      open(unit=utfil,file='rcsf.out',status='unknown')
      else
      open(fil_1,file='rcsf.out',status='unknown')
      open(unit=utfil,file='fil1.dat',status='unknown')
      endif
  100 close  (utfil)
      rewind (fil_1)
      return
      end
