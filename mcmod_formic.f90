module mcmod_mass
  use pes,wp=>pes_wp
  implicit none
  double precision::               V0
  integer,parameter::              atom1=1, atom2=2, atom3=3
  integer::                        n, ndim, ndof, natom, xunit, totdof

  !PES variables
  integer(kind=4), parameter :: nki(0:2)=(/4,4,2/)
  integer, parameter :: nk=10
  integer, parameter :: iord(0:nk-1)=(/1,2,3,4,5,6,7,8,9,10/)
  real(kind=wp) :: x1_cf, y1_cf, z1_cf
  type(cx_t) :: x4y4z2_pc = cx_null
  real(kind=wp), allocatable :: x4y4z2_cf(:)

contains

  subroutine V_init()
    character(len=80) :: dirname = "coef"
    integer (kind=4):: iun
    integer::  nb
    logical :: b0
    character(len=255) :: chd

    call pes_getiun(iun)
    b0 = dirname(len_trim(dirname):len_trim(dirname)).eq.'/'
    if (b0) then
       chd = dirname
    else
       chd = trim(dirname)//'/'
    endif

    write (*,*) 'Principal data directory: ', chd(1:len_trim(chd))
    write (*,*) 'Reading pcf-x4y4z2.dat'
    open (iun, status='old', file=trim(chd)//'pcf-x4y4z2.dat')
    read (iun,*) x4y4z2_pc
    read (iun,*) nb
    if (nb.ne.pes_x4y4z2_nb(x4y4z2_pc%dg)) then
       stop 'pes_init: x4y4z2 dimension error'
    endif

    allocate (x4y4z2_cf(0:nb-1))
    if (1.le.nb) then
       read (iun,*) x4y4z2_cf
    endif
    close (iun)

    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x)
    implicit none
    real(kind=wp),dimension(:,:),intent(in)::x
    double precision:: V
    real (kind=wp) :: r(0:nk-1,0:nk-1)

    call pes_dists(x, r)
    V = cx_f442(nki, r, x4y4z2_pc, x4y4z2_cf) +379.193815017316d0

    return
  end function V

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision::     grad(:,:), x(:,:)
    double precision::     potplus, potminus
    double precision, parameter::  eps=1.0d-5

    do i= 1,ndim
       do j=1,natom
          x(i,j)= x(i,j) + eps
          potplus= V(x)
          x(i,j)= x(i,j) - 2.0d0*eps
          potminus= V(x)
          x(i,j)= x(i,j) + eps
          grad(i,j)= (potplus-potminus)/(2.0d0*eps)
       end do
    end do
    return
  end subroutine Vprime
  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1, eps
    integer::              i, j
    double precision, allocatable::     gradplus(:, :), gradminus(:, :)

    eps=1d-4
    allocate(gradplus(ndim, natom), gradminus(ndim, natom))
    do i= 1, ndim
       do j= 1, natom
          x(i,j)= x(i,j) + eps
          call Vprime(x, gradplus)
          x(i,j)= x(i,j) - 2.0d0*eps
          call Vprime(x, gradminus)
          x(i,j)= x(i,j) + eps
          hess(i,j,:,:)= (gradplus(:,:)-gradminus(:,:))/(2.0d0*eps)          
       end do
    end do
    deallocate(gradplus, gradminus)
    return
  end subroutine Vdoubleprime

end module mcmod_mass
