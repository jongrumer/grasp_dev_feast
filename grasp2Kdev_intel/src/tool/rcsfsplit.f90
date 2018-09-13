program rcsfsplit

! Written by Per Jonsson, August 2014

implicit none
integer :: i, j, k, n, nlayer, norb, norblayer, nsymmetrymatch, norbcomp, ncsf, nwrite, ncount
character(len=100) :: string1, string2, string3, name
character(len=1500) :: orbitalstring
character(len=3) :: orb(300),orblayer(25),orbcomp(300)
character(len=4) :: orbrel(300)
character(len=100) :: orbitallayer(15),label(15)

write(*,*) 'RCSFSPLIT'
write(*,*) 'Splits a list name.c of CSFs into a number of lists with CSFs that '
write(*,*) 'can be formed from different sets of active orbitals. '
write(*,*) 'Orbital sets are specified by giving the orbitals with highest n for '
write(*,*) 'each symmetry. Orbitals occupy three positions and are right justified e.g.'
write(*,*) ' 4s, 4p, 4d '
write(*,*) 'Input file: name.c'
write(*,*) 'Output files: namelabel1.c, namelabel2.c, ...'
write(*,*) 

write(*,*) 'Name of state'
read(*,'(a)') name
write(*,*) 'Number of orbital sets'
read(*,*) nlayer
do i = 1,nlayer
   write(*,*) 'Give orbital set',i
   read(*,'(a)') orbitallayer(i)
   write(*,*) 'Give file label'
   read(*,'(a)') label(i)
end do

open(unit=36,file=trim(name)//'.c',status='old')

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

do k = 1,nlayer
   open(unit=48+k,file=trim(name)//trim(label(k))//'.c',status='unknown')

! Number of orbitals in orbital layer

   norblayer = (len_trim(orbitallayer(k))+1)/4

! Get individual orbital layer strings non-relativistic notation

   read(unit=orbitallayer(k),fmt='(25(a3,x))') orblayer(1:25)

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
   write(48+k,'(a)') trim(string1)
   read(36,'(a)') string1
   write(48+k,'(a)') trim(string1)
   read(36,'(a)') string1
   write(48+k,'(a)') trim(string1)
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
   write(48+k,'(a)') trim(orbitalstring)
   read(36,'(a)') string1
   write(48+k,'(a)') trim(string1)

   ncsf = 0
   do 
      read(36,'(a)',end=99) string1
      if (string1(2:2).eq.'*') then
         write(48+k,'(a)') trim(string1) 
         read(36,'(a)') string1
      end if
      read(36,'(a)') string2
      read(36,'(a)') string3

!  A CSF should be kept if there are no complementary orbitals in the string

      n = 0
      do i = 1,norbcomp
         n = index(trim(string1),orbcomp(i))
         if (n.ne.0) exit
      end do
      if (n.eq.0) then
         write(48+k,'(a)') trim(string1) 
         write(48+k,'(a)') trim(string2) 
         write(48+k,'(a)') trim(string3) 
         ncsf = ncsf + 1
      end if
   end do

99 continue

   write(*,*) 'List ',trim(name)//trim(label(k))//'.c',' based on orbital set',k
   write(*,*) ncsf,'CSFs'
end do

end program rcsfsplit
   
