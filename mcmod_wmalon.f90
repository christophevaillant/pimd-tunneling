module mcmod_mass
  use pes,wp=>pes_wp
  implicit none
  double precision::               V0
  integer,parameter::              atom1=1, atom2=2, atom3=4
  integer::                        n, ndim, ndof, natom, xunit, totdof
  character, allocatable::         label(:)

contains
  subroutine V_init()
    character (len=*), parameter :: dname='./coefs/'
    call pes0_init (dname)
    call pes1_init (pes_h4c3o2s_sysall)

    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x)
    implicit none
    real(kind=wp),dimension(:,:),intent(in)::x
    real(kind=wp)::pot
    ! ::::::::::::::::::::
    real(kind=wp),dimension(3,size(x)/3)::xn
    double precision:: V

    ! xn=reshape(x,(/3,size(x)/3/))

    V=pes_h4c3o2s_pot(x)

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
