module mcmod_mass
  implicit none
  double precision::               V0
  integer,parameter::              atom1=1, atom2=2, atom3=3
  double precision, allocatable::  x0(:), y0(:)
  double precision::               a0, b0, rho0
  integer::                        n, ndim, ndof, natom, xunit, totdof
  integer::                        m, potcount
contains

  subroutine V_init()
    implicit none
    integer:: k
    double precision, parameter::    pi=3.14159265358979d0

    m=6
    a0=2.0d0
    b0=0.2d0
    rho0=3.0d0
    allocate(x0(m), y0(m))
    do k=1, m
       x0(k)= rho0*cos(dble(k)*2.0d0*pi/dble(m))
       y0(k)= rho0*sin(dble(k)*2.0d0*pi/dble(m))
    end do

    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x)
    implicit none
    double precision::     v, x(:,:), answer
    integer::              k, i

    answer=0.0d0
    do k=1,m
       answer= answer - 0.5*exp(-a0*((x(1,1)-x0(k))**2 + (x(2,1) - y0(k))**2))
       answer= answer - 0.5*exp(-b0*((x(1,1)-x0(k))**2 + (x(2,1) - y0(k))**2))
    end do
    V=answer - V0
    return
  end function V

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j, k
    double precision::     grad(:,:), x(:,:), dummy1
    double precision, allocatable:: gradtemp(:), dummy2(:)
    potcount=potcount+1
    grad(1,1)=0.0d0
    grad(2,1)=0.0d0
    do k=1,m
       grad(1,1)= grad(1,1) + a0*(x(1,1)-x0(k))*exp(-a0*((x(1,1)-x0(k))**2 + (x(2,1) - y0(k))**2))
       grad(1,1)= grad(1,1) + b0*(x(1,1)-x0(k))*exp(-b0*((x(1,1)-x0(k))**2 + (x(2,1) - y0(k))**2))
       grad(2,1)= grad(2,1) + a0*(x(2,1)-y0(k))*exp(-a0*((x(1,1)-x0(k))**2 + (x(2,1) - y0(k))**2))
       grad(2,1)= grad(2,1) + b0*(x(2,1)-y0(k))*exp(-b0*((x(1,1)-x0(k))**2 + (x(2,1) - y0(k))**2))
       ! write(*,*) k, grad(1,1), grad(2,1)
    end do

    return
  end subroutine Vprime

  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1
    double precision, allocatable::   dummy2(:), hesstemp(:)
    double precision, allocatable::     gradplus(:, :), gradminus(:, :)
    integer::              k, i, j
    double precision:: u, dvdu, d2vdu2, dudx, dudy, eps

    hess(:,:,:,:)=0.0d0
    do k=1, m
       u= (x(1,1)- x0(k))**2 + (x(2,1)-y0(k))**2
       dvdu= a0*exp(-a0*u) + b0*exp(-b0*u)
       d2vdu2= -a0**2*exp(-a0*u) - b0**2*exp(-b0*u)
       dudx= x(1,1) - x0(k)
       dudy= x(2,1) - y0(k)
       !d2Vdx2:
       hess(1,1,1,1)= (d2vdu2*dudx + dvdu)*dudx
       !d2Vdxdy:
       hess(2,1,1,1)= (d2vdu2*dudy + dvdu)*dudx
       hess(1,1,2,1)= (d2vdu2*dudy + dvdu)*dudx
       !d2Vdy2:
       hess(2,1,2,1)= (d2vdu2*dudy + dvdu)*dudy
    end do
    ! eps=1d-4
    ! allocate(gradplus(ndim, natom), gradminus(ndim, natom))
    ! do i= 1, ndim
    !    do j= 1, natom
    !       x(i,j)= x(i,j) + eps
    !       call Vprime(x, gradplus)
    !       x(i,j)= x(i,j) - 2.0d0*eps
    !       call Vprime(x, gradminus)
    !       x(i,j)= x(i,j) + eps
    !       hess(i,j,:,:)= (gradplus(:,:)-gradminus(:,:))/(2.0d0*eps)          
    !    end do
    ! end do
    ! deallocate(gradplus, gradminus)

    return
  end subroutine Vdoubleprime

end module mcmod_mass
