module mcmod_mass
  use malonaldehyde
  implicit none
  double precision::               V0, eps2=1.0d-5
  integer,parameter::              atom1=1, atom2=2, atom3=4
  integer::                        n, ndim, ndof, natom, xunit, totdof
  character, allocatable::         label(:)

contains
  !---------------------------------------------------------------------
  subroutine V_init()
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x)
    implicit none
    double precision::     v, x(:,:)
    double precision, allocatable:: dummy1(:),dummy2(:)
    
    call pes(x,0,V,dummy1,dummy2)
    V= V- V0
    return
  end function V

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision::     grad(:,:), x(:,:), dummy1
    double precision, allocatable:: gradtemp(:), dummy2(:)

    allocate(gradtemp(ndof))
    call pes(x,1,dummy1,gradtemp,dummy2)
    do i= 1,ndim
       do j=1,natom
          grad(i,j)= gradtemp(ndim*(j-1)+i)
       end do
    end do
    deallocate(gradtemp)
    return
  end subroutine Vprime

  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1
    double precision, allocatable::   dummy2(:), hesstemp(:)
    integer::              idof1, idof2, i1,j1,i2,j2, ij

    allocate(hesstemp(ndof*(ndof +1)/2),dummy2(ndof))
    hesstemp=0.0d0
    hess=0.0d0
    call pes(x,2,dummy1,dummy2,hesstemp)
    ij=0
    do idof1=1,ndof
       do idof2=1,idof1
          ij=ij+1
          i1= mod(idof1, ndim)
          if (i1.eq.0) i1=ndim
          j1= 1 + (idof1-i1)/ndim
          i2= mod(idof2, ndim)
          if (i2.eq.0) i2=ndim
          j2= 1 + (idof2-i2)/ndim

          hess(i1,j1,i2,j2)= hesstemp(ij)
          hess(i2,j2,i1,j1)= hesstemp(ij)
       end do
    end do
    deallocate(hesstemp, dummy2)
    return
  end subroutine Vdoubleprime

end module mcmod_mass
