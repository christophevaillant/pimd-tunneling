module mcmod_mass
  implicit none
  double precision::               V0
  integer,parameter::              atom1=1, atom2=4, atom3=7 !pick the big atoms, O, O and Ha
  integer::                        n, ndim, ndof, natom, xunit, totdof
  character(len=5), allocatable::               at_name(:)
  character, allocatable::         label(:)

contains

  subroutine V_init()
    V0=0.0d0
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x)
    implicit none
    double precision::     v, x(:,:)
    double precision, allocatable:: dummy1(:),dummy2(:), xtemp(:)
    integer::              i,j

    if (.not. allocated(at_name)) then
       allocate(at_name(natom))
       do i=1,natom
          at_name(i)= trim(label(i))//char(0)
       end do
    end if

    allocate(xtemp(natom*ndim),dummy1(natom*ndim))
    do i=1,ndim
       do j=1,natom
          xtemp((j-1)*ndim +i)= x(i,j)*0.529177d0
       end do
    end do
    call calc_pot_link2f90(xtemp, at_name, natom, V)
    V= V*1.59362d-3! - V0
    deallocate(xtemp, dummy1)
    return
  end function V

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision::     grad(:,:), x(:,:), dummy1
    double precision, allocatable:: gradtemp(:), xtemp(:)

    if (.not. allocated(at_name)) then
       allocate(at_name(natom))
       do i=1,natom
          at_name(i)= trim(label(i))//char(0)
       end do
    end if

    allocate(gradtemp(ndof), xtemp(natom*ndim))
    do i=1,ndim
       do j=1,natom
          xtemp((j-1)*ndim +i)= x(i,j)*0.529177d0
       end do
    end do
    call calc_pot_link2f90_g(xtemp, at_name, natom, dummy1, gradtemp)
    do i= 1,ndim
       do j=1,natom
          grad(i,j)= gradtemp(ndim*(j-1)+i)*1.59362d-3*0.529177d0
       end do
    end do
    deallocate(gradtemp, xtemp)
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
