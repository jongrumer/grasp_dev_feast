*     last edited August 1, 1996
      subroutine Blandc(varmax,cfmax,lock,med,minJ,maxJ,nmax,
     &                  posn,posl,lim)
      integer fil_1,fil_2
      parameter (fil_1=7, fil_2=8)
      integer varmax,cfmax,minJ,maxJ,nmax,posn(110),posl(110)
      integer cf,org(1:15,0:10),i,j,n,k,l,l1,antal,tal,antalc
      integer low(1:15,0:10),tot,lista(15000,15,0:10)
      integer ansats(1:15,0:10,0:1),par,lim(15)
      logical lock(1:15,0:10),finns
      logical lik(15000),med(1:15,0:10)
      integer start,stopp,skal,duplet,kvar
      character rad*500,orb(0:10)
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
      cf     = 0
      antalc = 0
      tot    = 0
      skal   = 20
      finns  = .FALSE.
      do 1 i=1,15
         do 1 j=0,min(10,i-1)
            org(i,j)  = 0
    1       low(i,j)  = 0
      open(unit=7,status='scratch')
      read(fil_2,1000) rad
      write(fil_1,1000) rad
      read(fil_2,1000) rad
      write(fil_1,1000) rad
      read(fil_2,1000) rad
      write(fil_1,1000) rad  
      read(fil_2,1000) 
      do 930 i=1,500
  930    rad(i:i) = ' '
      start =-2
      stopp = 0
      do 940 k=1,110
         i=posn(k)
         j=posl(k)
         if (med(i,j) .OR. .NOT.lock(i,j)) then 
            start                = start + 5
            stopp                = stopp + 5
            rad(start:start)     = CHAR(ICHAR('0')+i)
            rad(start+1:start+1) = orb(j)
            if (j.GE.1) then
               rad(start+2:start+2) = '-'
               start                = start + 5
               stopp                = stopp + 5
               rad(start:start)     = CHAR(ICHAR('0')+i)
               rad(start+1:start+1) = orb(j)
            endif
         endif
  940 continue
      write(fil_1,1000) rad
      read(fil_2,1000) rad
      write(fil_1,1000) rad
    2 read(fil_2,1000,end=200) rad
      read(fil_2,1000,end=200)
      read(fil_2,1000,end=200)
      tot = tot+1
      do 10 i=1,15
         do 10 j=0,min(10,i-1)
   10       lista(tot,i,j) = 0
      do 20 i=1,skal
         n = 9*(i-1)+3
         l = n+1
         tal = ichar(rad(n:n)) - ichar('0')
         if (tal.GE.1 .AND. tal.LE.15) then
            l1 = -1
            do 15 j=0,tal-1
   15          if (orb(j).EQ.rad(l:l)) l1=j
            if (l1 .EQ. -1) goto 30
         else
            goto 30
         endif
         antal = ichar(rad(l+3:l+3)) - ichar('0')
         if (antal.GE.0 .AND. antal.LE.9) then
            antal = antal*10
         else
            antal = 0
         endif                    
         antal = antal + ichar(rad(l+4:l+4)) - ichar('0')
   20    lista(tot,tal,l1) = lista(tot,tal,l1) + antal
   30 continue
      goto 2
  200 if (tot.EQ.0) then
         write(*,1005) 'Nothing in inputfile!'
         stop
      endif
      if (tot.LT.10) then
         write(*,2001) 'Number of csf:s in inputfile = ',tot
      elseif (tot.LT.100) then
         write(*,2002) 'Number of csf:s in inputfile = ',tot 
      elseif (tot.LT.1000) then
         write(*,2003) 'Number of csf:s in inputfile = ',tot 
      elseif (tot.LT.10000) then
         write(*,2004) 'Number of csf:s in inputfile = ',tot 
      else
         write(*,2005) 'Number of csf:s in inputfile = ',tot 
      endif  
      duplet = 0
      do 250 i=1,tot
  250    lik(i) = .FALSE.
      if (tot.GE.2) then
         do 305 i=1,tot-1
            if (.NOT.lik(i)) then
               do 302 j=i+1,tot
                  if (.NOT.lik(j)) then
                     do 300 k=1,nmax
                        do 300 l=0,min(10,k-1)
  300                      if (.NOT.(lista(i,k,l).EQ. 
     :                               lista(j,k,l))) goto 302
                     lik(j) = .TRUE.
                   endif
  302          continue 
            endif 
  305    continue
      endif
      do 310 i=1,tot
  310    if (lik(i)) duplet = duplet + 1
      if (duplet.LT.10) then
         write(*,2001) 'Number of duplicat csf:s in file = ',duplet
      elseif (duplet.LT.100) then
         write(*,2002) 'Number of duplicat csf:s in file = ',duplet 
      elseif (duplet.LT.1000) then
         write(*,2003) 'Number of duplicat csf:s in file = ',duplet 
      elseif (duplet.LT.10000) then
         write(*,2004) 'Number of duplicat csf:s in file = ',duplet 
      else
         write(*,2005) 'Number of duplicat csf:s in file = ',duplet 
      endif
      kvar = tot - duplet
      write(*,*)
      do 320 i=1,tot
         if (.NOT.lik(i)) then
            if (kvar.GT.1) then
               if (kvar.LT.10) then
                  write(*,1001) kvar,' csf:s still to be expanded.'
               elseif (kvar.LT.100) then
                  write(*,1002) kvar,' csf:s still to be expanded.' 
               elseif (kvar.LT.1000) then
                  write(*,1003) kvar,' csf:s still to be expanded.' 
               elseif (kvar.LT.10000) then
                  write(*,1004) kvar,' csf:s still to be expanded.' 
               else
                  write(*,1006) kvar,' csf:s still to be expanded.'
               endif
            else
               write(*,1005) 'The last csf is still to be expanded.'
            endif
            kvar = kvar-1
            par = 0
            do 315 k=1,nmax
               do 315 l=0,min(k-1,10)
  315             org(k,l) = lista(i,k,l)
            if (finns) then
               open(unit=21,status='scratch')
               call Blandb(org,nmax,varmax,lock,21,low,lim,
     :                     posn,posl,minJ,maxJ)
               rewind(21)
               call Mergeb(antalc)
               if (antalc.LT.10) then
                  write(*,2001) 'Number of uncoupled csf:s = ',antalc
               elseif (antalc.LT.100) then
                  write(*,2002) 'Number of uncoupled csf:s = ',antalc 
               elseif (antalc.LT.1000) then
                  write(*,2003) 'Number of uncoupled csf:s = ',antalc
               elseif (antalc.LT.10000) then 
                  write(*,2004) 'Number of uncoupled csf:s = ',antalc
               else 
                  write(*,2005) 'Number of uncoupled csf:s = ',antalc
               endif 
            else
               open(unit=20,status='scratch')
               call Blandb(org,nmax,varmax,lock,20,low,lim,
     :                     posn,posl,minJ,maxJ)
               rewind(20)
               finns = .TRUE.
               antalc = 0
               write(*,1005) 
     :             'The first configuration has been expanded.'
            endif
            if (antalc .GE. cfmax) then
	       write(*,1005) 
     :             'Maximum number of uncoupled csf:s exceeded'
       	       goto 345
            endif
         else
         endif
  320 continue
  345 write(*,*)
      write(*,1005) 'Preparing the couplings of the csf:s.'

  350 if (nmax.LT.15) then
         do 391 i=nmax+1,15
            do 391 j=0,min(10,i-1)
               ansats(i,j,0) = 0
  391          ansats(i,j,1) = 0
      endif
      cf   = 0
  490 do 491 i=1,15
         read(20,5000,end=492) (ansats(i,j,0),j=0,min(10,i-1))
  491    read(20,5000,end=492) (ansats(i,j,1),j=0,min(10,i-1))
       par = 0
       do 496 i=1,15
          do 496 j=0,min(10,i-1)
             do 496 k=0,min(j,1)
  496           par = mod(par+j*ansats(i,j,k),2)
      call Gen(ansats,posn,posl,skal,cf,.TRUE.,minJ,maxJ,par)              
      goto 490
  492 continue
      rewind(fil_1)
      if (cf.EQ.0) then
         write(*,1005) 'No configuration state has been generated.'
      elseif (cf.EQ.1) then
         write(*,1005) 'One configuration state has been generated.'
      elseif (cf.LT.10) then
         write(*,1001) cf,' configuration states have been generated.'
      elseif (cf.LT.100) then
         write(*,1002) cf,' configuration states have been generated.'
      elseif (cf.LT.1000) then
         write(*,1003) cf,' configuration states have been generated.'
      elseif (cf.LT.10000) then
         write(*,1004) cf,' configuration states have been generated.'
      elseif (cf.LT.100000) then
         write(*,1006) cf,' configuration states have been generated.'
      else
         write(*,*) cf,' configuration states have been generated.'
      endif
 1000 format(A)
 1001 format(' ',I1,A)
 1002 format(' ',I2,A)
 1003 format(' ',I3,A)
 1004 format(' ',I4,A)
 1005 format(' ',A)
 1006 format(' ',I5,A)
 2001 format(' ',A,I1,'.')   
 2002 format(' ',A,I2,'.')                                                
 2003 format(' ',A,I3,'.')
 2004 format(' ',A,I4,'.')                 
 2005 format(' ',A,I5,'.')                 

 5000 format(11I2)
      return
      end
