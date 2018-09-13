*     last edited July 30, 1996
      subroutine Kopp2(pos,rad3,J,Jprim,par,antko)
      integer pos,i,k,J(20),par,tal,Jprim(20),antko(20)
      logical first
      character rad3*200

      first = .TRUE.
      do 5 i=1,200
    5    rad3(i:i) = ' '
      i = max(1,pos-1)
      k = 9*pos-2 
      if (J(i).EQ.2*(J(i)/2)) then
         tal = J(i)/20
         if (tal.NE.0) then
            rad3(k+2:k+2) = char(48+tal)
         endif
         tal = J(i)/2-tal*10
         rad3(k+3:k+3) = char(48+tal)
      else
         tal = J(i)/10
         if (tal.NE.0) then
            rad3(k:k) = char(48+tal)
         endif
         tal = J(i)-tal*10
         rad3(k+1:k+1) = char(48+tal)
         rad3(k+2:k+3) = '/2'
      endif
      if (par.EQ.0) then
         rad3(k+4:k+4) = '+'
      else
         rad3(k+4:k+4) = '-'
      endif
      if (pos.GT.2) then
         if (Jprim(1).NE.0 .OR. 
     :       .NOT.(Jprim(1).EQ.0 .AND. antko(1).EQ.1)) first = .FALSE. 
         do 10 i=1,pos-2
            if (first .AND. (Jprim(i+1).NE.0 .OR. 
     :                 (Jprim(i+1).EQ.0 .AND. antko(i+1).NE.1))) then
               first = .FALSE. 
            elseif (.NOT.first .AND. Jprim(i+1).NE.0) then 
               k = 9*(i+1)
               if (J(i).EQ.2*(J(i)/2)) then
                  tal = J(i)/20
	          if (tal.NE.0) then
                     rad3(k+2:k+2) = char(48+tal)
                  endif
                  tal = J(i)/2-tal*10
                  rad3(k+3:k+3) = char(48+tal)
               else
                  tal = J(i)/10
                  if (tal.NE.0) then
                     rad3(k:k) = char(48+tal)
                  endif
                  tal = J(i)-tal*10
                  rad3(k+1:k+1) = char(48+tal)
                  rad3(k+2:k+3) = '/2'
               endif
            endif
   10    continue
      endif
      return
      end
