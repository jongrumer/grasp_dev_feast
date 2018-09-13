*     last edited July 31, 1996
      subroutine Matcin(lock,closed,med,varmax,cfmax,nmax,
     &                  minJ,maxJ,lim)
      integer logfil
      parameter (logfil=31)
      integer varmax,i,nmax,lmax,minJ,maxJ,lim(15)
      integer cfmax,nmaks,mshell,lmaks
      character X*1,orb(0:10)
      logical lock(1:15,0:10),closed(1:15,0:10),all,med(1:15,0:10),lima
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
      nmaks = 1
      lmaks = 0
      do 10 i=1,15
         do 10 j=0,min(10,i-1)
            if (med(i,j)) then
               nmaks = i
               lmaks = max(j,lmaks)
            endif
   10 continue 
   60 continue
      if (nmaks.LE.9) then
         write(*,201) 'Highest n-number? (',nmaks,'..15)'
      else
         write(*,202) 'Highest n-number? (',nmaks,'..15)'
      endif
      read(*,*,err=60) nmax
      nmax = max(nmax,nmaks)
      nmax = min(nmax,15)
      write(logfil,*) nmax,' Highest principal quantum number.'
   70 write(*,400) 'Highest l-number? (',orb(lmaks),'..',
     :             orb(min(10,nmax-1)),')'
      read(*,1000) X
      lmax = -1
      do 71 i=0,min(10,nmax-1)
   71    if (X.EQ.orb(i)) lmax=i
      lmax = max(lmaks,lmax)
*     if (lmax.EQ.-1) goto 70
      write(logfil,*) lmax,' Highest orbital angular momentum.'
      write(*,200) 'Are all these nl-subshells active? (n/*)'
      read(*,1000) X
      all = .NOT.(X.EQ.'n' .OR. X.EQ.'N')
      write(logfil,*) all,' all subshells active.'
      do 72 i=1,15
   72    lim(i) = 0
      if (nmax.GE.2) then
c******************* modified by yu zou, 3/6/00
c this option cannot run correctly. It is not provided at present.
        go to 80
        write(*,200) 'Limitations on population of n-subshells? (y/*)'
         read(*,1000) X
         lima = X.EQ.'y' .OR. X.EQ.'Y'
   80    lima=.false.
c******************* modified by yu zou, 3/6/00
         write(logfil,*) lima,
     :                      ' limitations on population of n-subshells.'
         if (lima) then
            mshell = 0
            do 85 i=1,nmax-1
               mshell = mshell + 2*i*i
   83          continue
               if (i.EQ.1) then
                  write(*,200)
     :                    'Minimum number of electrons with n=1? (0..2)'
               elseif (i.LT.10) then
                  if (mshell.LT.100) then
                     write(*,208) 'Minimum number of electrons with n<='
     :                            ,i,'? (0..',mshell,')'
                  else
                     write(*,208) 'Minimum number of electrons with n<='
     :                            ,i,'? (0..)'
                  endif
               else
                  write(*,202) 'Minimum number of electrons with n<=',i,
     :                         '? (0..)'
               endif
               read(*,*,err=83) lim(i)
               if (lim(i).GT.mshell) lim(i) = mshell
               write(logfil,*) lim(i),
     :                     ' is minimum number of electrons with n =',i
   85       continue
         endif
      endif
      do 150 i=1,15
         do 150 j=0,min(10,i-1)
            if (nmax.GE.i .AND. lmax.GE.j .AND. .NOT.closed(i,j)) then
               if (all) then
                  lock(i,j) = .FALSE.
               else
                  if (i.LT.10) then
                     write(*,204) i,orb(j),' inactive or active? ',
     :                           '(i/*)'
                  else
                     write(*,205) i,orb(j),' inactive or active? ',
     :                           '(i/*)'
                  endif
                  read(*,1000) X
                  write(logfil,*) X,i,orb(j),' inactive, active, etc...'
                  lock(i,j)   = X.EQ.'i' .OR. X.EQ.'I'
               endif
            else
               lock(i,j) = .TRUE.
               if (closed(i,j)) then
                  if (i.LT.10) then
                     write(*,204) i,orb(j),' is a closed shell.'
                  else
                     write(*,205) i,orb(j),' is a closed shell.'
                  endif
               endif
            endif
  150 continue
 1100 write(*,400) 'Resulting 2*J-number? lower, higher ',
     :              '(J=1 -> 2*J=2 etc.)'
      read(*,*,ERR=1100) minJ,maxJ
      write(logfil,*) minJ,' to',maxJ,' is the resulting term.'
  160 write(*,200) 'Number of excitations = ? (0..)'
      read(*,*,err=160) varmax
      write(logfil,*) varmax,' number of excitations.'
  170 write(*,400) 'Maximum number of uncoupled configuration',
     :             ' states? (0..)'
      read(*,*,err=170) cfmax
      write(logfil,*) cfmax,' maximum number '
      write(*,*)
      
  200 format(' ',A,I1,A,A,I1,A)
  201 format(' ',A,I1,A,A,I2,A)
  202 format(' ',A,I2,A,A,I1,A)
  203 format(' ',A,I2,A,A,I2,A)
  204 format(' ',I1,3A)
  205 format(' ',I2,3A)
  208 format(' ',A,I1,A,I2,A)
  300 format(' ',A,I1,A)
  301 format(' ',A,I2,A)
  400 format(' ',7A)
  401 format(' ',2A,I1,A)
  402 format(' ',2A,I2,A)
 1000 format(A,A,A)
 2000 format(I1,A)
      return
      end
