module mcmod_mass
  implicit none
  double precision::               V0
  integer,parameter::              atom1=1, atom2=2, atom3=3
  double precision::               omegaforce, r0
  integer::                        n, ndim, ndof, natom, xunit, totdof
  character, allocatable::         label(:)

contains
  subroutine V_init()
    omegaforce=10000.0d0
    r0=20.0d0
    V0=0.0
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x)
    implicit none
    double precision::     v, x(:,:), r, answer
    integer::              i,j

    r= sqrt(x(1,1)**2 + x(2,1)**2)
    answer= 0.5d0*omegaforce**2*(r-r0)**2
    ! answer= 1.0d0 -0.5*(exp(-2.0d0*(r-3.0d0)**2) + exp(-0.2d0*(r-3.0d0)**2))
    if (answer.ne.answer) then
       write(*,*) "NaN in pot!"
       write(*,*) x
       stop
    end if
    V=answer-V0
    return
  end function V

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision::     grad(:,:), x(:,:), r, u, dvdr
    double precision::     eplus, Eminus, xtest(2), eps
    double precision, allocatable:: gradtest(:,:)

    eps=1d-5
    r= sqrt(x(1,1)**2 + x(2,1)**2)
    grad(1,1)= omegaforce**2*x(1,1)*(1.0d0 - r0/r)
    grad(2,1)= omegaforce**2*x(2,1)*(1.0d0 - r0/r)
    if (grad(1,1).ne.grad(1,1) .or.grad(2,1).ne.grad(2,1) ) then
       write(*,*) "NaN in grad!"
       write(*,*) x
       write(*,*) grad
       stop
    end if

    ! allocate(gradtest(ndim, natom))
    ! eps=1d-6
    ! do i= 1, ndim
    !    do j= 1, natom
    !       x(i,j)= x(i,j) + eps
    !       eplus= V(x)
    !       x(i,j)= x(i,j) - 2.0d0*eps
    !       eminus= V(x)
    !       x(i,j)= x(i,j) + eps
    !       gradtest(i,j)= (Eplus-Eminus)/(2.0d0*eps)
    !       write(*,*)"-----------"
    !       write(*,*) i, j, gradtest(i,j)
    !       write(*,*) i,j,grad(i,j)
    !    end do
    ! end do
    ! deallocate(gradtest)

    return
  end subroutine Vprime
  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), r,u, dvdr, eps
    double precision, allocatable::     gradplus(:, :), gradminus(:, :)
    double precision, allocatable::   hesstest(:,:,:,:)
    integer::              i, j

    r= sqrt(x(1,1)**2 + x(2,1)**2)
    dvdr= omegaforce**2*(r-r0)
    hess(1,1,1,1)= x(1,1)**2*omegaforce**2*r0/r**3
    hess(1,1,2,1)= x(1,1)*x(2,1)*omegaforce**2*r0/r**3
    hess(2,1,1,1)= x(1,1)*x(2,1)*omegaforce**2*r0/r**3
    hess(2,1,2,1)= x(2,1)**2*omegaforce**2*r0/r**3

    ! u= (r-3.0d0)**2
    ! dvdr=(2.0d0 - 4.0d0*u)*exp(-2.0d0*u) + (0.2d0 - 4.0d-2*u)*exp(-0.2d0*u)
    ! hess(1,1,1,1)= dvdr*x(2,1)**2/r**3
    ! hess(2,1,2,1)= dvdr*x(1,1)**2/r**3
    ! hess(1,1,2,1)= dvdr*x(1,1)*x(2,1)/r**3
    ! hess(2,1,1,1)= hess(1,1,2,1)

    ! eps=1d-4
    ! allocate(hesstest(ndim, natom, ndim, natom))
    ! allocate(gradplus(ndim, natom), gradminus(ndim, natom))
    ! do i= 1, ndim
    !    do j= 1, natom
    !       x(i,j)= x(i,j) + eps
    !       call Vprime(x, gradplus)
    !       x(i,j)= x(i,j) - 2.0d0*eps
    !       call Vprime(x, gradminus)
    !       x(i,j)= x(i,j) + eps
    !       hesstest(i,j,:,:)= (gradplus(:,:)-gradminus(:,:))/(2.0d0*eps)
    !       write(*,*)"-----------"
    !       write(*,*) i, j, hesstest(i,j,:,1)
    !       write(*,*) i,j,hess(i,j,:,1)
    !    end do
    ! end do
    ! deallocate(gradplus, gradminus, hesstest)
    ! stop

    return
  end subroutine Vdoubleprime


end module mcmod_mass
