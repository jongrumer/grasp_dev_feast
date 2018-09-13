*     last edited July 30, 1996
      subroutine Kopp1(pos,rad2,J,S,antko)
      integer pos,i,k,J(20),tal,S(20),antko(20)
      character rad2*200
      do 5 i=1,200    
    5    rad2(i:i) = ' ' 
      do 10 i=1,pos
         k = 9*i
         if (J(i).EQ.2*(J(i)/2)) then
            if (.NOT.(J(i).EQ.0 .AND. antko(i).EQ.1)) then
               if (S(i).NE.-1) then
                  rad2(k-4:k-4) = ';'   
                  rad2(k-5:k-5) = char(48+S(i))
               endif
               tal = J(i)/20
               if (tal.NE.0) then
C Bug                  rad2(k:k) = char(48+tal)
                  rad2(k-1:k-1) = char(48+tal)
               endif
               tal = J(i)/2-tal*10
               rad2(k:k) = char(48+tal)
            endif
         else
            if (S(i).NE.-1) then
               rad2(k-4:k-4) = ';'
               rad2(k-5:k-5) = char(48+S(i))
            endif
            tal = J(i)/10
            if (tal.NE.0) then
               rad2(k-3:k-3) = char(48+tal)
            endif
            tal = J(i)-tal*10
            rad2(k-2:k-2) = char(48+tal)
            rad2(k-1:k) = '/2'
         endif
   10    continue

      return
      end
