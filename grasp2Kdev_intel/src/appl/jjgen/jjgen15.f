************************************************************************
*     ------------------------------------------------------------------
*     JJGEN  -- program to generate
*
*     Written by: Lennart Sturesson
*                 
*     2:nd version
*
*     Last edited Januar 2, 1997
*
*     ------------------------------------------------------------------
*
      program Jjgen15
      integer logfil
      parameter (logfil=31)
      integer org(1:15,0:10),varmax,skal,anel,par,
     :        low(1:15,0:10),posn(110),posl(110),nmax,
     :        lim(15),minJ,maxJ,cfmax
      logical lock(1:15,0:10),closed(1:15,0:10),slut,med(1:15,0:10),
     :        dubbel(1:15,0:10),advexp,second  
      character X
      open(unit=logfil,file='rcsf.log',status='unknown')
      write(*,*)
      write(*,*) 'RCSFGENERATE'
      write(*,*) 'This program generates lists of CSFs'
      write(*,*) 'Outputfiles: rcsf.out, rcsf.log'
      write(*,*) 
   10 write(*,200) '    *  : new list'
      write(*,200) '    a  : add to existing list'
      write(*,200) '    e  : expand existing list'
      write(*,200) '    q  : quit'
      read(*,100) X
      write(logfil,200) ' Option : ',X
      call Reffa(posn,posl)
      advexp = .FALSE.
      if (X.EQ.'a' .OR. X.EQ.'A') then
         call Adder(closed,med,slut,anel,par,.FALSE.)
         if (slut) then
            write(*,200)
            write(*,200) 'The clist.inp-file is not readable! '
            stop
         endif
         write(logfil,200) ' New reference set.'
         second = .TRUE.
      elseif (X.EQ.'e' .OR. X.EQ.'E') then
         advexp = .TRUE.
         call Adder(closed,med,slut,anel,par,.TRUE.)
         if (slut) then
            write(*,200)
            write(*,200) 'The clist.inp-file is not readable! '
            stop
         endif
         write(logfil,200) ' File as reference sets.'
         call 
     :      Matcin(lock,closed,med,varmax,cfmax,nmax,minJ,maxJ,lim)
         call 
     :      Blandc(varmax,cfmax,lock,med,minJ,maxJ,nmax,posn,posl,lim)
         second = .FALSE.
      else
         call Matain(org,lock,closed,varmax,skal,nmax,anel,par,
     :                              low,minJ,maxJ,lim,dubbel)
         call Fivelines(org,lock,closed,.TRUE.,posn,posl)
         call Blanda(org,varmax,lock,minJ,maxJ,skal,nmax,low,
     :                            posn,posl,lim,dubbel,.TRUE.) 
         second = .FALSE.
      endif
      if (.NOT.advexp)
     :   call Matbin(org,lock,closed,varmax,skal,second,anel,
     :                      par,low,nmax,lim,dubbel,minJ,maxJ) 
      if (second) then
         call Fivelines(org,lock,closed,.FALSE.,posn,posl)
         call Blanda(org,varmax,lock,minJ,maxJ,skal,nmax,low,
     :                          posn,posl,lim,dubbel,.FALSE.) 
         call Merge(.FALSE.,posn,posl)
C         write(*,200) 'The merged file is called rcsf.out.'
      else
         call Merge(.TRUE.,posn,posl)
C         write(*,200) 'The generated file is called rcsf.out.'
      endif


      call rcsfblock


      stop
  100 format(A)
  200 format(' ',10A)
  300 format(' ',A,I2,A)
  400 format(' ',A,I3,A)
      end
