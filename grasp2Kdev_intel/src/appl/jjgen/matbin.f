*     last edited August 2, 1996
      subroutine Matbin(org,lock,closed,varmax,skal,second,anel0,
     :                  par0,low,nmax,lim,dubbel,minJ,maxJ)    
     
      integer logfil
      parameter (logfil=31)
      integer org(1:15,0:10),anel,par,anel0,par0,low(1:15,0:10)
      integer varmax,i,j,skal,nmax,lmax,em,nenter,anela,anelb
      integer lim(15),block,mshell,minJ,maxJ,tmp
      character X*1,orb(0:10),L(0:20),Y*2
      logical lock(1:15,0:10),closed(1:15,0:10),second,all,
     :        dubbel(1:15,0:10),lima
      data (L(i),i=0,20)/'S','P','D','F','G','H','I','K','L','M','N'
     :                         ,'O','Q','R','T','U','V','W','X','Y','Z'/
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/

   40 if (.NOT. second) then
         write(*,*)
         write(*,*) 'You have the possibility to generate another list'
         write(*,*) 'This list must have the same 2*J values as       '
         write(*,*) 'previous lists of the same parity                '
         write(*,200) 'Generate another list? (y/*)'
         read(*,1000) X
         second = X.EQ.'y' .OR. X.EQ.'Y'
         write(logfil,*) second,' Generate another list.'
         if (.NOT.second) return
      endif
      anel  = 0
      anela = 0
      anelb = 0
      par   = 0
      skal = 20
   60 write(*,200) 'Highest n-number? (1..15)'
      read(*,*,err=60) nmax
      nmax = max(nmax,1)
      nmax = min(nmax,15)
      write(logfil,*) nmax,' Highest principal quantum number.'
   70 write(*,400) 'Highest l-number? (s..',orb(min(10,nmax-1)),')'
      read(*,1000) X
      lmax = -1
      do 71 i=0,min(10,nmax-1)
   71    if (X.EQ.orb(i)) lmax=i
      if (lmax.EQ.-1) goto 70
      write(logfil,*) lmax,' Highest orbital angular momentum.'
      write(*,200) 'Are all these nl-subshells active? (n/*)'
      read(*,1000) X
      all = .NOT.(X.EQ.'n' .OR. X.EQ.'N')
      write(logfil,*) all,' all subshells active.'
      do 72 i=1,15
   72    lim(i) = 0
      if (nmax.GE.2) then
        write(*,200) 'Limitations on population of n-subshells? (y/*)'
         read(*,1000) X
         lima = X.EQ.'y' .OR. X.EQ.'Y'
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
   95 continue
      if (nmax.LT.10) then
         write(*,200)
     :      'Highest n-number in reference configuration? (1..',nmax,')'
      else
         write(*,202)
     :      'Highest n-number in reference configuration? (1..',nmax,')'
      endif
      read(*,*,err=95) nenter
      nenter = max(nenter,1)
      nenter = min(nenter,nmax)
      write(logfil,*) nenter,' highest n-number.'
      block  = 0
      do 160 i=1,15
         do 150 j=0,min(10,i-1)
            low(i,j)    = 0
            dubbel(i,j) = .FALSE.
            if (nmax.GE.i .AND. lmax.GE.j .AND. .NOT.closed(i,j)) then
               if (nenter.GE.i) then
                  em = 2 + 4*j
                  if (em.LT.10) then
  100                continue
                     if (i.LE.9) then
                        write(*,200) 'Number of electrons in ',i,orb(j),
     :                                  '? (0..',em,')'
                     else
                        write(*,202) 'Number of electrons in ',i,orb(j),
     :                                  '? (0..',em,')'
                     endif
                     read(*,*,err=100) org(i,j)
                     if (org(i,j).LT.0 .OR. org(i,j).GT.em) goto 100
                  else
  101                continue
                     if (i.LT.10) then
                        write(*,201) 'Number of electrons in ',i,orb(j),
     :                                  '? (0..',em,')'
                     else
                        write(*,203) 'Number of electrons in ',i,orb(j),
     :                                  '? (0..',em,')'
                     endif
                     read(*,*,err=101) org(i,j)
                     if (org(i,j).LT.0 .OR. org(i,j).GT.em) goto 101
                  endif
                  write(logfil,*) org(i,j),' number of electrons in',i,
     :                        orb(j)
                  if (all) then
                     lock(i,j) = .FALSE.
                  else
                     if (org(i,j).GT.1) then
                        if (org(i,j).LE.10) then
                           write(*,201)
     :          'Inactive, active or minimum? (i/*/0..',org(i,j)-1,')'
                        else 
                           write(*,202)
     :          'Inactive, active or minimum? (i/*/0..',org(i,j)-1,')'
                        endif
                        read(*,1000) Y
                     elseif (org(i,j).EQ.1) then
                        write(*,201) 'Inactive or active? (i/*)'
                        read(*,1000) Y
                     else
                        write(*,400) 'Inactive, active or doubled  ',
     :                               'excited? (i/*/d)'
                        read(*,1000) Y
                        dubbel(i,j) =  Y(1:1).EQ.'d' .OR.  Y(1:1).EQ.'D'
                     endif
                     lock(i,j) =  Y(1:1).EQ.'i' .OR.  Y(1:1).EQ.'I' 
                     if (Y(1:1).GE.'0' .AND. Y(1:1).LE.'9') then
                        if (org(i,j).GT.0) then
                           tmp = ICHAR(Y(1:1))-ICHAR('0')
                           if (Y(2:2).GE.'1' .AND. Y(2:2).LE.'9')
     :                        tmp = tmp*10 + ICHAR(Y(2:2))-ICHAR('0')
                           low(i,j) = min(org(i,j),tmp)
                        endif
                     endif
                     write(logfil,1000) Y,' inactive, active, etc...'
                  endif
                  if (.NOT.lock(i,j)) anela = anela + org(i,j)
                  anel = anel + org(i,j)
                  par  = mod(par+j*org(i,j),2)
               elseif (all) then
                  org(i,j) = 0
                  lock(i,j) = .FALSE.
               else
                  org(i,j) = 0
                  closed(i,j) = .FALSE.
                  if (i.LT.10) then
                     write(*,204) i,orb(j),' inactive, active or ',
     :                         'doubled excited? (i/*/d)'
                  else
                     write(*,205) i,orb(j),' inactive, active or ',
     :                         'doubled excited? (i/*/d)'
                  endif
                  read(*,1000) X
                  dubbel(i,j) = X.EQ.'d' .OR. X.EQ.'D'
                  lock(i,j)   = X.EQ.'i' .OR. X.EQ.'I'
                  write(logfil,*) X,i,orb(j),' inactive, active, etc...'
               endif
            else
               org(i,j)  = 0
               lock(i,j) = .TRUE.
               if (closed(i,j)) then
                  if (i.LT.10) then
                     write(*,204) i,orb(j),' is a closed shell.'
                  else
                     write(*,205) i,orb(j),' is a closed shell.'
                  endif
                  em    = 2 + 4*j
                  anel  = anel + em
                  block = block + em
               endif
            endif
            anelb = anelb + low(i,j)
  150    continue
         lim(i) = lim(i) - block
         if (lim(i).LT.0) lim(i) = 0
  160 continue
      if (anel.NE.anel0) then
         if (anel0.LT.10) write(*,300)
     :             'Wrong number of electrons. The first list had ',
     :                 anel0,' electrons.'
         if (anel0.GE.10) write(*,301)
     :             'Wrong number of electrons. The first list had ',
     :                 anel0,' electrons.'
         if (anel.LT.10) write(*,300) 'This list has ',anel,
     :                                 ' electrons.'
         if (anel.GE.10) write(*,301) 'This list has ',anel,
     :                                 ' electrons.'
         second = .FALSE.
         goto 40
      endif
 1100 write(*,400) 'Resulting 2*J-number? lower, higher ',
     :              '(J=1 -> 2*J=2 etc.)'
      read(*,*,ERR=1100) minJ,maxJ
      if (anel .EQ. 2*(anel/2)) then
         if (minJ .NE. 2*(minJ/2) .OR. maxJ .NE. 2*(maxJ/2)) then
            write(*,*) 'The resulting 2*J-numbers should be even'
            goto 1100
         endif
      else
         if (minJ .EQ. 2*(minJ/2) .OR. maxJ .EQ. 2*(maxJ/2)) then
            write(*,*) 'The resulting 2*J-numbers should be odd'
            goto 1100
         endif
      endif
      write(logfil,*) minJ,' to',maxJ,' is the resulting term.'
*     if (par.NE.par0) then
*        write(*,200) 'Wrong parity.'
*        if (par0.EQ.0) write(*,*)
*    :           'The first list had even parity and this list has odd.'
*        if (par0.EQ.1) write(*,*)
*    :           'The first list had odd parity and this list has even.'
*        second = .FALSE.
*        goto 40
*     endif
      par0 = par
      anelb = anela - anelb
 1200 continue
      if (anelb.LT.10) then
         write(*,200) 'Number of excitations = ? (0..',anelb,')'
         read(*,*,err=1200) varmax
      else
         write(*,202) 'Number of excitations = ? (0..',anelb,')'
         read(*,*,err=1200) varmax
      endif
      write(logfil,*) varmax,' number of excitations.'
  200 format(' ',A,I1,A,A,I1,A)
  201 format(' ',A,I1,A,A,I2,A)
  202 format(' ',A,I2,A,A,I1,A)
  203 format(' ',A,I2,A,A,I2,A)
  204 format(' ',I1,3A)
  205 format(' ',I2,3A)
  208 format(' ',A,I1,A,I2,A)
  300 format(' ',A,I1,A)
  301 format(' ',A,I2,A)
  400 format(' ',3A)
  401 format(' ',2A,I1,A)
  402 format(' ',2A,I2,A)
 1000 format(A,A,A)
 2000 format(I1,A)
      return
      end
