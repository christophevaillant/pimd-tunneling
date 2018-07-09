module mcmod_mass
  implicit none
  double precision::               V0
  integer,parameter::              atom1=1, atom2=2, atom3=3
  integer::                        n, ndim, ndof, natom, xunit, totdof

contains

  subroutine V_init()

    call init_ccpol(3, 1,1,0)

    write(*,*) "Potential initializaton complete"
    V0=0.0d0
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x)
    double precision:: V, x(:,:), xtemp(ndim,natom), Etot
    double precision, parameter::  ang=0.529177d0
    integer::  i,j
    xtemp(:,:)= x(:,:)*ang
    call ccpol(Etot, xtemp(:,1), xtemp(:,2), xtemp(:,3), xtemp(:,4), xtemp(:,5), xtemp(:,6))
    ! if (abs(Etot) .lt. 5.d3) then
       V=(Etot/627.510d0)  -V0
    ! else
    !    V=0.0d0
    !    write(*,*) "cut:", Etot
    !    do i=1, natom
    !       write(*,*) (x(j,i)*ang, j=1,ndim)
    !    end do
    !    stop
    ! end if

    return

  end function V

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision::     grad(:,:), x(:,:)
    double precision::     potplus, potminus, eps
    eps=1d-4
    do i= 1,ndim
       do j=1,natom
          x(i,j)= x(i,j) + eps
          potplus= V(x)
          x(i,j)= x(i,j) - 2.0d0*eps
          potminus= V(x)
          x(i,j)= x(i,j) + eps
          grad(i,j)= (potplus-potminus)/(2.0d0*eps)
          ! write(*,*) i,j, x(i,j),grad(i,j)
       end do
    end do
    return
  end subroutine Vprime
  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1, eps
    integer::              i, j
    double precision::     gradplus(ndim, natom), gradminus(ndim, natom)
    eps=1d-5
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
    return
  end subroutine Vdoubleprime


end module mcmod_mass
