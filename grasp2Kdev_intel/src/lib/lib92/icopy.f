      subroutine  icopy(n,ix,incx,iy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      integer ix(*),iy(*)
      integer i,incx,incy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        iy(i) = ix(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        iy(i) = ix(i)
        iy(i + 1) = ix(i + 1)
        iy(i + 2) = ix(i + 2)
        iy(i + 3) = ix(i + 3)
        iy(i + 4) = ix(i + 4)
        iy(i + 5) = ix(i + 5)
        iy(i + 6) = ix(i + 6)
   50 continue
      return
      end
