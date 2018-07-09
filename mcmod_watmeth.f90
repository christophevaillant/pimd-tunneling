module mcmod_mass
  use watermethane_mod
  implicit none
  double precision::               V0
  integer,parameter::              atom1=1, atom2=2, atom3=3
  integer::                        n, ndim, ndof, natom, xunit, totdof

contains

  subroutine V_init()
    V0=0.0d0
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x)
    implicit none
    double precision::     v, x(:,:)
    double precision, allocatable:: dummy1(:)
    integer::              i,j

    !the input x for this routine is in cartesian coordinates
    !but the coordinates *have* to be varied in a rigid way
    !otherwise the potential gives you nonsense.
    call wmrb(x, dummy1, V, .false.)

    return
  end function V

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision::     grad(:,:), x(:,:), dummy1
    double precision, allocatable:: gradtemp(:), xtemp(:)
    
    call wmrb_grad(x, grad)
    
    return
  end subroutine Vprime
  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1, eps
    integer::              i, j
    double precision, allocatable::     gradplus(:, :), gradminus(:, :)

    !could probably work out an analytic hessian
    !but not necessarily worth it yet
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
