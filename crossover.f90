!Program to calculate crossover temperatures for a given transition state
!on a classical minimum energy path (i.e. the top of the barrier)
program crossover
use instantonmod
use mcmod_mass
implicit none

integer::                        i,j, i1, i2, j1, j2, idof1, idof2
integer::                        lwork, liwork, info, dummyint
character::                      dummylabel
double precision, allocatable::  tstate(:,:), hess(:,:,:,:), hessmat(:,:)
double precision, allocatable::  work(:), etasquared(:)
integer,allocatable::            iwork(:)
namelist /CROSSDATA/ ndim, natom,xunit

  ndim=3
  natom=1
  xunit=2

  read(5, nml=CROSSDATA)

  ndof= ndim*natom

  allocate(mass(natom),label(natom))
  open(18, file="masses.dat", status="old")
  do j=1,natom
     read(18,*)label(j), mass(j)
  end do
  close(18)


allocate(tstate(ndim,natom), hess(ndim, natom, ndim, natom))
allocate(etasquared(ndim*natom), hessmat(ndim*natom, ndim*natom))
open(20, file="transition.xyz")

read(20,*) dummyint
read(20,*) dummylabel

do i=1, natom
      read(20,*) label(i), (tstate(j,i), j=1,ndim)
end do

if (xunit .eq. 2) then
   tstate(:,:)= tstate(:,:)/0.529177d0
end if

hess(:,:,:,:)=0.0d0

call Vdoubleprime(tstate, hess)

do i1=1, ndim
   do j1=1,natom
      do i2=1, ndim
         do j2=1,ndim
            idof1= (j1-1)*ndim + i1
            idof2= (j2-1)*ndim + i2
            hessmat(idof1,idof2)= hess(i1,j1,i2,j2)/sqrt(mass(j1)*mass(j2))
         end do
      end do
   end do
end do

lwork= 2*ndof +1
liwork= 1

allocate(work(lwork), iwork(liwork))

call dsyevd('N', 'U', ndim*natom, hessmat, ndim*natom, etasquared,work,lwork,iwork,liwork,info)
write(*,*) etasquared(:)

write(*,*) "Crossover beta=", 2.0*pi/sqrt(-etasquared(1))

end program
