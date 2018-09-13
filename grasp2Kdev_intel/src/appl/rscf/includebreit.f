      program includebreit

! Written by Per Jonsson, September 2014

      implicit none
      integer :: i, j, k, n, nlayer, norb, norblayer, nsymmetrymatch
      integer ::  norbcomp, ncsf, nwrite, ncount, nblock, ncsfblock(50)
      integer, allocatable :: nbreit(:,:)
      character(len=300) :: string1, string2, string3, name
      character(len=1500) :: orbitalstring
      character(len=3) :: orb(300),orblayer(25),orbcomp(300)
      character(len=4) :: orbrel(300)
      character(len=100) :: orbitallayer

      ncsfblock = 0

      write(*,*) 'Name of state'
      read(*,'(a)') name
      write(*,*) 'Give orbital set'
      read(*,'(a)') orbitallayer

      open(unit=36,file=trim(name)//'.c',status='old')
      open(unit=521,status='scratch')


! Read five first line save line four containing the string of orbitals

      read(36,'(a)') string1
      read(36,'(a)') string1
      read(36,'(a)') string1
      read(36,'(a)') orbitalstring
      read(36,'(a)') string1

! Number of orbitals

      norb = (len_trim(orbitalstring)+1)/5

! Get individual orbital strings non-relativistic notation

      read(unit=orbitalstring,fmt='(300(x,a3,x))') orb(1:300)

! Number of orbitals in orbital layer

      norblayer = (len_trim(orbitallayer)+1)/4

! Get individual orbital layer strings non-relativistic notation

      read(unit=orbitallayer,fmt='(25(a3,x))') orblayer(1:25)

! For current orbital layer find out the compliment orbitals 

      norbcomp = 0
      do i = 1,norb
         nsymmetrymatch = 0
         do j = 1,norblayer
            if (orb(i)(3:3).eq.orblayer(j)(3:3)) then
               nsymmetrymatch = 1
               if (lgt(orb(i)(1:2),orblayer(j)(1:2))) then
                  norbcomp = norbcomp + 1
                  orbcomp(norbcomp) =  orb(i)
               end if
            end if
         end do
         if (nsymmetrymatch.eq.0) then
            norbcomp = norbcomp + 1
            orbcomp(norbcomp) =  orb(i)
         end if 
      end do 
         
      write(*,*)

      rewind(36)
      read(36,'(a)') string1
      read(36,'(a)') string1
      read(36,'(a)') string1
      read(36,'(a)') orbitalstring

!  Find out and write the orbitals for this layer in relativistic notation

      read(unit=orbitalstring,fmt='(300(x,a4))') orbrel(1:300)
      orbitalstring = ''
      ncount = 0
      do i = 1,norb
         nwrite = 0
         do j = 1,norbcomp
            if (orbrel(i)(1:3).eq.orbcomp(j)) then
               nwrite = 1
            end if
         end do

!  Write out the ones not in the complement orbital set

         if (nwrite.eq.0) then
            orbitalstring = orbitalstring(1:ncount*5)//' '//orbrel(i)
            ncount = ncount + 1
         end if
      end do
      read(36,'(a)') string1

      nblock = 1
      ncsfblock(nblock) = 0
      do 
         read(36,'(a)',end=99) string1
         if (string1(2:2).eq.'*') then
            read(36,'(a)') string1
            nblock = nblock + 1
            ncsfblock(nblock) = 0
         end if
         read(36,'(a)') string2
         read(36,'(a)') string3
         ncsfblock(nblock) = ncsfblock(nblock) + 1 

!  A CSF should be kept if there are no complementary orbitals in the string

         n = 0
         do i = 1,norbcomp
            n = index(trim(string1),orbcomp(i))
            if (n.ne.0) exit
         end do
         if (n.eq.0) then
            write(*,*) '   Breit',ncsfblock(nblock),trim(string1)
            write(521,*) 1
         else 
            write(*,*) 'No Breit',ncsfblock(nblock),trim(string1)
            write(521,*) 0
         end if
      end do

   99 continue

      rewind(unit=521)

      do i = 1,nblock
         write(*,*) i,ncsfblock(i)
      end do

      allocate(nbreit(maxval(ncsfblock),nblock))

      do i = 1,nblock
         do j = 1,ncsfblock(i)
            read(521,*) nbreit(j,i)
            write(139,*) i,j,'xxx',nbreit(j,i)
         end do
      end do


      end program includebreit
   
