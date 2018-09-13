         call Merge(.TRUE.,posn,posl)
         call copy9t7 
         do  
         call Matbin(org,lock,closed,varmax,skal,second,anel,
     :                      par,low,nmax,lim,dubbel,minJ,maxJ) 
         if(.not.second) exit
         call Fivelines(org,lock,closed,.FALSE.,posn,posl)
         call Blanda(org,varmax,lock,minJ,maxJ,skal,nmax,low,
     :                          posn,posl,lim,dubbel,.FALSE.) 
         call Merge(.F.,posn,posl)
         call copy9t7 
         second = .f.
         enddo
C         write(*,200) 'The merged file is called rcsf.out.'
