module instantonmod
  use mcmod_mass
  implicit none
  double precision, parameter::    pi=3.14159265358979d0
  double precision::               beta, betan, UMtilde
  double precision, allocatable::  xtilde(:,:,:),well1(:,:), well2(:,:), mass(:)
  logical::                        fixedends

  public :: QsortC
  private :: Partition
contains


  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !linear polymer potential
  function UM(x,a,b)
    implicit none
    integer::            i,j,k
    double precision,intent(in)::   x(:,:,:)
    double precision,intent(in),optional:: a(:,:),b(:,:)
    double precision:: UM, pot

    UM=0.0d0
    do i=1, N-1, 1
       pot= V(x(i,:,:),i)
       write(*,*) i, pot
       UM=UM + pot
       do j=1, ndim
          do k=1, natom
             UM=UM+ (0.5d0*mass(k)/betan**2)*(x(i+1,j,k)-x(i,j,k))**2
          end do
       end do
    end do
    UM=UM+ V(x(N,:,:),i)
    if (fixedends) then
       !UM=UM + V(a(:,:))+ V(b(:,:)) !already set these to 0!
       do j=1, ndim
          do k=1, natom
             UM=UM+ (0.5d0*mass(k)/betan**2)*(x(1,j,k)-a(j,k))**2
             UM=UM+ (0.5d0*mass(k)/betan**2)*(b(j,k)-x(N,j,k))**2
          end do
       end do
    end if

    return
  end function UM

  !---------------------------------------------------------------------
  !linear polymer force
  subroutine UMprime(x, answer,a,b)
    implicit none
    integer::            i,j,k
    double precision, intent(in)::   x(:,:,:)
    double precision, intent(out)::  answer(:,:,:)
    double precision, intent(in), optional:: a(:,:),b(:,:)
    double precision,allocatable:: grad(:,:)

    allocate(grad(ndim,natom))
    do i=1, N
       do j=1,ndim
          do k=1,natom
             if (i.eq.1) then
                if (fixedends) then
                   answer(1,j,k)=mass(k)*(2.0*x(1,j,k) - a(j,k) - x(2,j,k))/betan**2
                else
                   answer(1,j,k)=mass(k)*(x(1,j,k) - x(2,j,k))/betan**2                   
                end if
             else if (i.eq.N) then
                if (fixedends) then
                   answer(N,j,k)=mass(k)*(2.0*x(N,j,k) - x(N-1,j,k) - b(j,k))/betan**2
                else
                   answer(N,j,k)=mass(k)*(x(N,j,k) - x(N-1,j,k))/betan**2
                end if
             else
                answer(i,j,k)= mass(k)*(2.0*x(i,j,k) - x(i-1,j,k) - x(i+1,j,k))/betan**2
             end if
          end do
       end do
       call Vprime(x(i,:,:),grad(:,:),i)
       do j=1, ndim
          do k=1, natom
             answer(i,j,k)= answer(i,j,k)+ grad(j,k)
          end do
       end do
    end do
    deallocate(grad)

    return
  end subroutine UMprime
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !linear polymer force and energy
  subroutine UMforceenergy(x, answer,UM,a,b)
    implicit none
    integer::            i,j,k
    double precision,intent(in)::   x(:,:,:)
    double precision,intent(in),optional:: a(:,:),b(:,:)
    double precision, intent(out)::  answer(:,:,:), UM
    double precision,allocatable:: grad(:,:)
    double precision:: energy

    allocate(grad(ndim,natom))
    UM=0.0d0
    do i=1, N, 1
       call potforce(x(i,:,:),grad,energy,i)
       UM=UM+ energy
       do j=1, ndim
          do k=1, natom
             if (i .lt. N) UM=UM+ (0.5d0*mass(k)/betan**2)*(x(i+1,j,k)-x(i,j,k))**2
          end do
       end do
       do j=1,ndim
          do k=1,natom
             if (i.eq.1) then
                if (fixedends) then
                   answer(1,j,k)=mass(k)*(2.0*x(1,j,k) - a(j,k) - x(2,j,k))/betan**2
                else
                   answer(1,j,k)=mass(k)*(x(1,j,k) - x(2,j,k))/betan**2                   
                end if
             else if (i.eq.N) then
                if (fixedends) then
                   answer(N,j,k)=mass(k)*(2.0*x(N,j,k) - x(N-1,j,k) - b(j,k))/betan**2
                else
                   answer(N,j,k)=mass(k)*(x(N,j,k) - x(N-1,j,k))/betan**2
                end if
             else
                answer(i,j,k)= mass(k)*(2.0*x(i,j,k) - x(i-1,j,k) - x(i+1,j,k))/betan**2
             end if
          end do
       end do
       do j=1, ndim
          do k=1, natom
             answer(i,j,k)= answer(i,j,k)+ grad(j,k)
          end do
       end do
    end do

    if (fixedends) then
       ! UM=UM + V(a(:,:))+ V(b(:,:))
       do j=1, ndim
          do k=1, natom
             UM=UM+ (0.5d0*mass(k)/betan**2)*(x(1,j,k)-a(j,k))**2
             UM=UM+ (0.5d0*mass(k)/betan**2)*(b(j,k)-x(N,j,k))**2
          end do
       end do
    end if

    deallocate(grad)

    return
  end subroutine UMforceenergy

  !---------------------------------------------------------------------
  !linear polymer hessian
  subroutine UMhessian(x, singlewell,answer)
    implicit none
    integer::            i, j1, k1, j2, k2, idof1, idof2
    integer::            fulldof1, fulldof2, index
    double precision, intent(in)::   x(:,:,:)
    logical, intent(in)::   singlewell
    double precision, intent(out):: answer(:,:)
    double precision, allocatable:: hess(:,:,:,:)

    allocate(hess(ndim, natom, ndim, natom))

    answer=0.0d0
    hess=0.0d0
    do i=1, n, 1
       if ((i .eq. 1 .and. singlewell) .or. .not. singlewell) then
          call Vdoubleprime(x(i,:,:), hess,i)
       end if
       do j1=1,ndim
          do k1=1,natom
             do j2=1,ndim
                do k2=1,natom
                   !The DoF label for the individual bead
                   idof1= (k1-1)*ndim + j1
                   idof2= (k2-1)*ndim + j2
                   !The DoF label for the whole matrix
                   fulldof1=ndof*(i-1) + idof1
                   fulldof2=ndof*(i-1) + idof2

                   if (fulldof2 .lt. fulldof1) cycle

                   if (idof1.eq.idof2) then
                      answer(1,fulldof1)= 2.0d0/betan**2&
                           +hess(j2,k2,j1,k1)/sqrt(mass(k1)*mass(k2)) 
                      if (i.gt.1) answer(ndof+1,fulldof1)=-1.0d0/betan**2
                   else
                      index=1 + fulldof2 -fulldof1
                      if (index.lt.0) cycle
                      answer(index,fulldof1)= &
                        +hess(j2,k2,j1,k1)/sqrt(mass(k1)*mass(k2)) 
                   end if
                end do
             end do
          end do
       end do
    end do
    deallocate(hess)
    return
  end subroutine UMhessian

  !---------------------------------------------------------------------
  !Align a vector of atoms

  subroutine get_align(atomsin,theta1, theta2, theta3, origin)
    implicit none
    double precision,intent(in)::     atomsin(:,:)
    double precision,intent(out)::     theta1, theta2, theta3, origin(ndim)
    double precision::     workvec(3), atoms(ndim,natom)
    integer::              i,j,k


    atoms(:,:)=0.0d0 
    theta1=0.0d0
    theta2=0.0d0
    theta3=0.0d0
    if (ndim .ne. 3) then
       write(*,*) "Wrong number of dimensions; change align_atoms subroutine!"
       stop
    end if
    !-----------------------------------------
    !Put atom1 at origin
    origin(:)= atomsin(:,atom1)

    do i=1, natom
       atoms(:,i)= atomsin(:,i) - origin(:)
    end do
    !-----------------------------------------
    !Align vector between atom1 and atom2 to x axis
    !first rotate about z-axis to align with zx plane
    workvec(:)= atoms(:,atom2) - atoms(:,atom1)
    theta1= atan2(workvec(2),workvec(1))
    call rotate_atoms(atoms, 3, theta1)

    !rotate about y-axis to align with z-axis
    workvec(:)= atoms(:,atom2) - atoms(:,atom1)
    theta2= atan2(workvec(3), workvec(1))
    call rotate_atoms(atoms, 2, theta2)

    !-----------------------------------------
    !Align vector between atom1 and atom3 to xz plane
    !rotate about x-axis
    workvec(:)= atoms(:,atom3) - atoms(:,atom1)
    theta3= -atan2(workvec(2),workvec(3))
    
    return
  end subroutine get_align

  !---------------------------------------------------------------------
  !Align a vector of atoms

  subroutine align_atoms(atomsin, theta1,theta2,theta3, origin, atomsout)
    implicit none
    double precision,intent(in)::     atomsin(:,:),theta1, theta2, theta3 , origin(ndim)
    double precision, intent(out):: atomsout(:,:)
    double precision::     workvec(3)
    integer::              i,j,k

    atomsout(:,:)=0.0d0
    workvec(:)=0.0d0
    if (ndim .ne. 3) then
       write(*,*) "Wrong number of dimensions; change align_atoms subroutine!"
       stop
    end if

    if (any((/abs(theta1).gt.1D-10,abs(theta2).gt.1D-10,abs(theta3).gt.1D-10/))) then
       !-----------------------------------------
       !Put atom1 at origin
       do i=1, natom
          atomsout(:,i)= atomsin(:,i) - atomsin(:,atom1)
       end do
       !-----------------------------------------
       !Align vector between atom1 and atom2 to x axis
       !first rotate about z-axis to align with zx plane
       workvec(:)= atomsout(:,atom2) - atomsout(:,atom1)
       call rotate_atoms(atomsout, 3, theta1)

       !rotate about y-axis to align with z-axis
       workvec(:)= atomsout(:,atom2) - atomsout(:,atom1)
       call rotate_atoms(atomsout, 2, theta2)

       !-----------------------------------------
       !Align vector between atom1 and atom3 to xz plane
       !rotate about x-axis
       workvec(:)= atomsout(:,atom3) - atomsout(:,atom1)
       call rotate_atoms(atomsout, 1, theta3)
    else
       atomsout(:,:)= atomsin(:,:)
    end if
    ! call centreofmass(atomsout, workvec)
    ! do i=1,natom
    !    atomsout(:,i)=atomsout(:,i) - workvec(:)
    ! end do
    return
  end subroutine align_atoms

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !general rotation of a vector of angle theta about axis
  subroutine rotate_vec(vec,axis,theta)
    implicit none
    double precision::     rotmatrix(3,3), vec(:)
    double precision::     theta
    integer::              i, axis, j,k

    if (axis.eq. 1) then
       j=2
       k=3
    else if(axis.eq.2) then
       j=1
       k=3
    else
       j=1
       k=2
    end if
    rotmatrix(:,:)=0.0d0
    rotmatrix(axis,axis)=1.0d0
    rotmatrix(j,j)= cos(theta)
    rotmatrix(k,k)= cos(theta)
    rotmatrix(k,j)= -sin(theta)
    rotmatrix(j,k)= sin(theta)
    vec(:)= matmul(rotmatrix(:,:), vec(:))
    return
  end subroutine rotate_vec

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !rotate the atoms by an angle theta about an axis
  subroutine rotate_atoms(atoms,axis,theta)
    implicit none
    double precision,intent(in):: theta
    integer, intent(in)::  axis
    double precision, intent(inout):: atoms(:,:)
    double precision::     rotmatrix(3,3)
    integer::              i,j,k

    if (abs(theta) .gt. 1D-10) then
       if (axis.eq. 1) then
          j=2
          k=3
       else if(axis.eq.2) then
          j=1
          k=3
       else
          j=1
          k=2
       end if
       rotmatrix(:,:)=0.0d0
       rotmatrix(axis,axis)=1.0d0
       rotmatrix(j,j)= cos(theta)
       rotmatrix(k,k)= cos(theta)
       rotmatrix(k,j)= -sin(theta)
       rotmatrix(j,k)= sin(theta)
       do i=1, natom
          atoms(:,i)= matmul(rotmatrix(:,:), atoms(:,i))
       end do
    end if
    return
  end subroutine rotate_atoms

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !spline algorithms from numerical recipes
  FUNCTION assert_eq(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq
    if (n1 == n2 .and. n2 == n3) then
       assert_eq=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq'
    end if
  END FUNCTION assert_eq

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  SUBROUTINE spline(x,y,yp1,ypn,y2)
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(IN) :: x,y
    double precision, INTENT(IN) :: yp1,ypn
    double precision, DIMENSION(:), INTENT(OUT) :: y2
    ! Given arrays x and y of length N containing a tabulated function, i.e., yi = f (xi ), with x1 <
    ! x2 < . . . < xN , and given values yp1 and ypn for the first derivative of the interpolating
    ! function at points 1 and N , respectively, this routine returns an array y2 of length N
    ! that contains the second derivatives of the interpolating function at the tabulated points
    ! xi . If yp1 and/or ypn are equal to 1 × 10 30 or larger, the routine is signaled to set the
    ! corresponding boundary condition for a natural spline, with zero second derivative on that
    ! boundary.
    INTEGER :: n
    double precision, DIMENSION(size(x)) :: a,b,c,r
    n=assert_eq(size(x),size(y),size(y2),"spline")
    c(1:n-1)=x(2:n)-x(1:n-1)
    ! Set up the tridiagonal equations.
    r(1:n-1)=6.0d0*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0d0*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    if (yp1 > 0.99d30) then
       ! The lower boundary condition is set either to be “natural”
       r(1)=0.0
       c(1)=0.0
    else
       ! or else to have a specified first derivative.
       r(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       c(1)=0.5
    end if
    ! The upper boundary condition is set either to be
    ! or else to have a specified first derivative.

    if (ypn > 0.99d30) then
       ! “natural”
       r(n)=0.0
       a(n)=0.0
    else
       r(n)=(-3.0d0/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
       a(n)=0.5
    end if
    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
    return
  END SUBROUTINE spline
  
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  SUBROUTINE tridag(a,b,c,r,u)
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(IN) :: a,b,c,r
    double precision, DIMENSION(:), INTENT(OUT) :: u
    ! Solves for a vector u of size N the tridiagonal linear set given by equation (2.4.1) using a
    ! serial algorithm. Input vectors b (diagonal elements) and r (right-hand sides) have size N ,
    ! while a and c (off-diagonal elements) are size N − 1.
    double precision, DIMENSION(size(b)) :: gam
    ! One vector of workspace, gam is needed.
    INTEGER:: n,j
    double precision :: bet
    n=assert_eqn((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),"tridag_ser")
    bet=b(1)
    if (bet == 0.0) then
       write(*,*) "tridag_ser: Error at code stage 1"
       stop
    end if
    ! If this happens then you should rewrite your equations as a set of order N − 1, with u2
    ! trivially eliminated.
    u(1)=r(1)/bet
    do j=2,n
       ! Decomposition and forward substitution.
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j-1)*gam(j)
       if (bet == 0.0) then
          write(*,*) "tridag_ser: Error at code stage 2"
          stop
       end if
       u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1
       ! Backsubstitution.
       u(j)=u(j)-gam(j+1)*u(j+1)
    end do
  END SUBROUTINE tridag

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) "nrerror: an assert_eq failed with this tag:", &
            string
       STOP "program terminated by assert_eqn"
    end if
  END FUNCTION assert_eqn

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  FUNCTION splint(xa,ya,y2a,x)
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    double precision, INTENT(IN) :: x
    double precision :: splint
    ! Given the arrays xa and ya, which tabulate a function (with the xai "s in increasing or
    ! decreasing order), and given the array y2a, which is the output from spline above, and
    ! given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
    ! and y2a are all of the same size.
    INTEGER :: khi,klo,n
    double precision :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),"splint")
    klo=max(min(locate(xa,x),n-1),1)
    ! We will find the right place in the table by means of locate"s bisection algorithm. This is
    ! optimal if sequential calls to this routine are at random values of x. If sequential calls are in
    ! order, and closely spaced, one would do better to store previous values of klo and khi and
    ! test if they remain appropriate on the next call.
    khi=klo+1
    ! klo and khi now bracket the input value of x.
    h=xa(khi)-xa(klo)
    if (h == 0.0) then
       write(*,*) "bad xa input in splint"
       stop
    end if
    ! The xa"s must be distinct.
    a=(xa(khi)-x)/h
    ! Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0
  END FUNCTION splint

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  FUNCTION splin_grad(xa,ya,y2a,x)
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    double precision, INTENT(IN) :: x
    double precision :: splin_grad
    ! Given the arrays xa and ya, which tabulate a function (with the xai "s in increasing or
    ! decreasing order), and given the array y2a, which is the output from spline above, and
    ! given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
    ! and y2a are all of the same size.
    INTEGER :: khi,klo,n
    double precision :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),"splin_grad")
    klo=max(min(locate(xa,x),n-1),1)
    ! We will find the right place in the table by means of locate"s bisection algorithm. This is
    ! optimal if sequential calls to this routine are at random values of x. If sequential calls are in
    ! order, and closely spaced, one would do better to store previous values of klo and khi and
    ! test if they remain appropriate on the next call.
    khi=klo+1
    ! klo and khi now bracket the input value of x.
    h=xa(khi)-xa(klo)
    if (h == 0.0) then
       write(*,*) "bad xa input in splin_grad"
       stop
    end if
    ! The xa"s must be distinct.
    a=(xa(khi)-x)/h
    ! Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h
    splin_grad=((ya(khi) - ya(klo))/h) +((1.0d0-3.0d0*a**2)*y2a(klo)+(3.0d0*b**2-1.0d0)*y2a(khi))*h/6.0d0
  END FUNCTION splin_grad

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  FUNCTION locate(xx,x)
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(IN) :: xx
    double precision, INTENT(IN) :: x
    INTEGER :: locate
    ! Given an array xx(1:N ), and given a value x, returns a value j such that x is between
    ! xx(j ) and xx(j + 1) . xx must be monotonic, either increasing or decreasing. j = 0 or
    ! j = N is returned to indicate that x is out of range.
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    ! True if ascending order of table, false otherwise.
    jl=0
    ! Initialize lower
    ju=n+1
    ! and upper limits.
    do
       if (ju-jl <= 1) exit
       ! Repeat until this condition is satisfied.
       jm=(ju+jl)/2
       ! Compute a midpoint,
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
          ! and replace either the lower limit
       else
          ju=jm
          ! or the upper limit, as appropriate.
       end if
    end do
    if (x == xx(1)) then
       ! Then set the output, being careful with the endpoints.
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if
  END FUNCTION locate

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !sort algorithms
recursive subroutine QsortC(A)
  double precision, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
subroutine Partition(A, marker)
  double precision, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  double precision :: temp
  double precision :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

  !---------------------------------------------------------------------
!---------------------------------------------------------------------
!returns the centre of mass
subroutine centreofmass(x, com)
  implicit none
  double precision::    x(:,:), com(:)
  integer::             i,j

  com(:)=0.0d0
  do i=1,ndim
     do j=1,natom
        com(i)= com(i) + mass(j)*x(i,j)
     end do
  end do
  com(:)=com(:)/sum(mass(:))
  return
end subroutine centreofmass

  !---------------------------------------------------------------------
!---------------------------------------------------------------------
!optimizes and returns the instanton
  subroutine instanton(xtilde,a,b)
    implicit none
    double precision, intent(inout)::xtilde(:,:,:)
    double precision, intent(in), optional:: a(:,:),b(:,:)
    integer::                        iprint, m, iflag, mp,idof, maxiter
    integer::                        i, lp, count, iw, j,k
    double precision::               xtol, gtol, stpmin, stpmax
    double precision::               f, xtemp(ndim,natom)

    double precision::               factr, com(ndim)
    double precision, allocatable::  fprime(:,:,:), work(:), fprimework(:)
    double precision, allocatable::  lb(:), ub(:), dsave(:), xwork(:)
    integer, allocatable::           nbd(:), iwork(:), isave(:)
    logical::                        lsave(4)
    character(len=60)::              task, csave

    allocate(lb(totdof), ub(totdof),fprime(n,ndim,natom), nbd(totdof))
    allocate(fprimework(totdof), xwork(totdof))
    do i=1, n, 1
       do j=1,ndim
          do k=1,natom
             idof= ((k-1)*ndim + j -1)*n +i
             if (fixedends) then
                lb(idof)= a(j,k)
                ub(idof)= a(j,k)
             else
                lb(idof)= 0.0d0
                ub(idof)= 0.0d0
             end if
             nbd(idof)=0
          end do
       end do
    end do
      
    !------------------------
    !perform minimization
    task='START'
    m=8
    iprint=-1
    xtol= 1d-5
    iw=totdof*(2*m+5) + 11*m**2 + 8*m
    allocate(work(iw), iwork(3*totdof), isave(44), dsave(29))
    iflag=0
    ! eps2= 1.0d-5 !gradient convergence
    factr=1.0d6
    maxiter=40
    if (fixedends) then
       if (potforcepresent) then
          call UMforceenergy(xtilde,fprime,f,a,b)
       else
          f= UM(xtilde,a,b)
          call UMprime(xtilde,fprime,a,b)
       end if
    else
       if (potforcepresent) then
          call UMforceenergy(xtilde,fprime,f)
       else
          f= UM(xtilde)
          call UMprime(xtilde,fprime)
       end if
    end if
    count=0
    do while( task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
         task.eq.'START')
       count=count+1
       xwork=reshape(xtilde,(/totdof/))
       fprimework= reshape(fprime,(/totdof/))
       call setulb(totdof,m,xwork,lb,ub,nbd,f,fprimework,factr,eps2,work&
            ,iwork,task,iprint, csave,lsave,isave,dsave,maxiter)
       if (task(1:2) .eq. 'FG') then
          write(*,*) "iteration", count
          xtilde= reshape(xwork,(/n,ndim,natom/))
          if (fixedends) then
             if (potforcepresent) then
                call UMforceenergy(xtilde,fprime,f,a,b)
             else
                f= UM(xtilde,a,b)
                call UMprime(xtilde,fprime,a,b)
             end if
          else
             if (potforcepresent) then
                call UMforceenergy(xtilde,fprime,f)
             else
                f= UM(xtilde)
                call UMprime(xtilde,fprime)
             end if
          end if
       end if
    end do
    if (task(1:5) .eq. "ERROR" .or. task(1:4) .eq. "ABNO") then
       write(*,*) "Error:"
       write(*,*) task
    end if

    deallocate(work, lb, ub, fprime, fprimework,xwork)
    deallocate(iwork, nbd, isave, dsave)

    return
  end subroutine instanton

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !returns the fluctuation prefactor
  subroutine detJ(x, etasquared, singlewell)
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol
  double precision,intent(in)::    x(:,:,:)
  logical, intent(in)::             singlewell
  double precision,intent(inout)::  etasquared(:)
  integer::                        nout, ldz, lwork, liwork, info,i
  integer,allocatable::            isuppz(:), iwork(:)
  double precision, allocatable::  work(:), z(:,:), H(:,:)
  ! !get diagonal hessian
  jobz='N'
  range='A'
  uplo='U'
  abstol=1.0d-8
  lwork= 2*totdof
  liwork= 1
  info=0
  ldz=totdof
  vl=0.0d0
  vu=0.0d0
  nout=0
  allocate(work(lwork), iwork(liwork), H(ndof+1,totdof))
  H=0.0d0
  etasquared=0.0d0
  call UMhessian(x,singlewell,H)
  write(*,*) "Hessian is cooked."

  call DSBEVD('N', 'L', totdof, ndof, H, ndof+1, etasquared, z, 1, work, lwork, iwork, liwork, info)
  deallocate(work, iwork, H)
  return
  end subroutine detJ

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !returns the centre of a path
  function findmiddle(x1,x2,lampath,path, splinepath)
    integer::          jmax,j
    double precision:: findmiddle,x1,x2,xacc
    double precision:: lampath(:), path(:), splinepath(:)
    parameter (jmax=40)
    parameter (xacc= 1d-6)
    double precision:: dx,f,fmid,xmid

    fmid=splin_grad(lampath, path, splinepath, x2)
    f=splin_grad(lampath, path, splinepath, x1)
    if(f*fmid.ge.0.) then
       write(*,*) 'root must be bracketed in findmiddle'
       write(*,*) 'you probably did a terrible job at finding the right path'
       write(*,*) x1,f
       write(*,*) x2,fmid
       stop
    end if
    if(f.lt.0.)then
       findmiddle=x1
       dx=x2-x1
    else
       findmiddle=x2
       dx=x1-x2
    endif
    do j=1,JMAX
       dx=dx*.5
       xmid=findmiddle+dx
       fmid=splin_grad(lampath, path, splinepath, xmid)
       if(fmid.le.0.)findmiddle=xmid
       if(abs(dx).lt.xacc .or. fmid.eq.0.) return
    end do
    write(*,*)'too many bisections in findmiddle'
    stop
  end function findmiddle

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !calculates the euclidean distance between two configurations
  function eucliddist(x1, x2)
    implicit none
    double precision,intent(in)::   x1(:,:), x2(:,:)
    double precision::  eucliddist
    integer::            i,j
    
    eucliddist=0.0d0
    do i=1, ndim
       do j=1, natom
          eucliddist= eucliddist+ (x1(i,j)- x2(i,j))**2
       end do
    end do
    eucliddist= sqrt(eucliddist)
    return
  end function eucliddist

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine read_path(instapath,centre,npath,path,Vpath, splinepath, lampath)
    implicit none
    integer, intent(inout)::   npath
    logical, intent(in)::   instapath, centre
    double precision,allocatable, intent(out):: path(:,:,:), Vpath(:), splinepath(:), lampath(:)
    double precision, allocatable::   initpath(:,:)
    double precision:: a,b,xmiddle, origin(3), theta1, theta2, theta3
    integer::   i,j,k, dummy
    character:: dummylabel, dummystr(28)    

    !----------------------------------
    !allocate and initialize
    allocate(lampath(npath),initpath(ndim, natom))
    allocate(path(npath, ndim, natom), Vpath(npath),splinepath(npath))


    path(:,:,:)= 0.0d0
    initpath(:,:)=0.0d0

    !----------------------------------
    !open file and read in
    open(15, file="path.xyz")
    do i=1, npath
       read(15,*) dummy
       read(15,'(28A)') dummystr
       do j=1, natom
          read(15,*) dummylabel, (initpath(k,j), k=1,ndim)
       end do
       if (xunit .eq. 2) then
          initpath(:,:)= initpath(:,:)/0.529177d0
       end if
       if (i.eq.1) then
          lampath(1)=0.0d0
          call get_align(initpath,theta1, theta2, theta3, origin)
          call align_atoms(initpath,theta1, theta2, theta3, origin, path(i,:,:))
       else
          call align_atoms(initpath,theta1, theta2, theta3, origin, path(i,:,:))
          lampath(i)= lampath(i-1) + eucliddist(path(i-1,:,:), path(i,:,:))
       end if
       Vpath(i)= V(path(i,:,:))
    end do
    lampath(:)= lampath(:)/lampath(npath)
    deallocate(initpath)
    close(15)

    !----------------------------------
    !calculate splines for interpolation
    do i=1,ndim
       do j=1,natom
          splinepath(:)=0.0d0
          call spline(lampath(:), path(:,i,j), 1.0d31, 1.0d31, splinepath(:))
          do k=1,n
             xtilde(k,i,j)= splint(lampath, path(:,i,j), splinepath(:), dble(k-1)/dble(n-1))
          end do
       end do
    end do
    !----------------------------------
    !find instanton from initial guess
    if (instapath) then
       write(*,*) "Finding instanton"
       if (fixedends) then
          call instanton(xtilde,well1,well2)
       else
          call instanton(xtilde)
       end if
       deallocate(path,lampath,vpath,splinepath)
       if (fixedends) then
          npath=n+2
          allocate(path(npath, ndim, natom), Vpath(npath),splinepath(npath),lampath(npath))
          path(2:n+1,:,:)=xtilde(:,:,:)
          path(1,:,:)= well1(:,:)
          path(npath,:,:)= well2(:,:)
       else
          npath=n
          allocate(path(npath, ndim, natom), Vpath(npath),splinepath(npath),lampath(npath))
          path(:,:,:)= xtilde(:,:,:)
       end if
       write(*,*) "Found instanton."
       open(19, file="instanton.xyz")
       do i=1,n
          write(19,*) natom
          write(19,*) "Energy of minimum",i
          do j=1, natom
             write(19,*)  label(j), (xtilde(i,k,j)*0.529177d0, k=1,ndim)
          end do
       end do

       !work out reaction coordinate
       do i=1, npath
          if (i.eq.1) then
             lampath(1)=0.0d0
          else
             lampath(i)= lampath(i-1) + eucliddist(path(i-1,:,:), path(i,:,:))
          end if
          Vpath(i)= V(path(i,:,:))
       end do
       lampath(:)= lampath(:)/lampath(npath)
       close(19)
    end if

     !-------------------------
     !Find the centre to make sure this is symmetric
     if (centre) then
        write(*,*) size(lampath), size(Vpath), size(splinepath)
        xmiddle= findmiddle(0.3d0, 0.7d0, lampath, Vpath, splinepath)
        a= 2.0d0  - 4.0d0*xmiddle
        b= 4.0d0*xmiddle - 1.0d0
        do i=1, npath
           if (a.ge.0) then
              lampath(i)= -0.5d0*b/a + sqrt((lampath(i)/a) + (0.5*b/a)**2)
           else
              lampath(i)= -0.5d0*b/a - sqrt((lampath(i)/a) + (0.5*b/a)**2)
           end if
        end do
        write(*,*)"Centred,", a, b, xmiddle
     end if

     !----------------------------------
     !calculate splines for interpolation
     write(*,*) size(lampath), size(Vpath), size(splinepath)
     do i=1,ndim
        do j=1,natom
           splinepath(:)=0.0d0
           call spline(lampath(:), path(:,i,j), 1.0d31, 1.0d31, splinepath(:))
        end do
     end do
     open(20, file="aligned.xyz")
     do i=1,n
        write(20,*) natom
        write(20,*) "rotation angles:", theta1, theta2, theta3
        do j=1, natom
           write(20,*)  label(j), (xtilde(i,k,j)*0.529177d0, k=1,ndim)
        end do
     end do
     close(20)

    return
  end subroutine read_path

end module instantonmod
