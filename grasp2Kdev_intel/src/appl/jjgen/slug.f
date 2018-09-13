
*     last edited November 2, 1995
      subroutine Slug(i,j,varmax,varupp,varned,ansats,org,lock,dubbel,
     :                                                 low,start,stopp)
      integer i,j,varmax,varupp(1:15,0:10),varned(1:15,0:10),minmax
      integer ansats(1:15,0:10,0:1),org(1:15,0:10),start,stopp,iold,jold
      integer low(1:15,0:10)
      logical lock,dubbel(1:15,0:10)
      if (i.EQ.1) then
         varupp(1,0) = 0
         varned(1,0) = 0
      else
         if (j.EQ.0) then
            iold = i - 1
            jold = min(10,iold-1)
         else
            iold = i
            jold = j - 1
         endif
         if (jold.EQ.0) then
            varupp(i,j) = varupp(iold,jold) +
     :                 max(0,ansats(iold,jold,0)-org(iold,jold))
            varned(i,j) = varned(iold,jold) +
     :                 max(0,org(iold,jold)-ansats(iold,jold,0))
         else
            varupp(i,j) = varupp(iold,jold) + max(0,ansats(iold,jold,0)
     :                 +ansats(iold,jold,1)-org(iold,jold))
            varned(i,j) = varned(iold,jold) + max(0,org(iold,jold)
     :                 -ansats(iold,jold,0)-ansats(iold,jold,1))
         endif
      endif
      if (lock) then
         start = org(i,j)
         stopp = org(i,j)
         return
      endif
      if (j.GE.5) then
         minmax = 4
      else
         minmax = 4*j + 2
      endif
      start = min(minmax,org(i,j)+(varmax-varupp(i,j)))
      if (dubbel(i,j)) start = 2*(start/2)
      stopp = max(low(i,j),org(i,j)-(varmax-varned(i,j)))
      return
      end
