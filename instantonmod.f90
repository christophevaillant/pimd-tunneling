module instantonmod
  use mcmod_mass
  implicit none
  double precision, parameter::    pi=3.14159265358979d0
  double precision::               beta, betan, UMtilde
  double precision, allocatable::  well1(:,:), well2(:,:), mass(:)
  character, allocatable::         label(:)
  logical::                        fixedends

  public :: QsortC
  private :: Partition
contains


  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  function UM(x,a,b)
    implicit none
    integer::            i,j,k
    double precision::   x(:,:,:), UM,a(:,:),b(:,:)

    UM=0.0d0
    do i=1, N-1, 1
       UM=UM+ V(x(i,:,:))
       do j=1, ndim
          do k=1, natom
             UM=UM+ (0.5d0*mass(k)/betan**2)*(x(i+1,j,k)-x(i,j,k))**2
          end do
       end do
    end do
    UM=UM+ V(a(:,:))+ V(b(:,:))+ V(x(N,:,:))
    if (fixedends) then
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
  subroutine UMprime(x,a,b, answer)
    implicit none
    integer::            i,j,k
    double precision::   x(:,:,:), answer(:,:,:),a(:,:),b(:,:)
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
       call Vprime(x(i,:,:),grad(:,:))
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
  subroutine UMhessian(x, answer)
    implicit none
    integer::            i, j1, k1, j2, k2, idof1, idof2
    integer::            fulldof1, fulldof2, index
    double precision::   x(:,:,:), answer(:,:)
    double precision, allocatable:: hess(:,:,:,:)

    allocate(hess(ndim, natom, ndim, natom))

    answer=0.0d0
    do i=1, n, 1
       hess=0.0d0
       call Vdoubleprime(x(i,:,:), hess)
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
    double precision::     atomsin(:,:), workvec(3), atoms(ndim,natom)
    double precision::     theta1, theta2, theta3, origin(ndim)
    integer::              i,j,k, atom1, atom2, atom3


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
    write(*,*) "Alignment angles:", theta1, theta2, theta3
    
    return
  end subroutine get_align

  !---------------------------------------------------------------------
  !Align a vector of atoms

  subroutine align_atoms(atomsin, theta1,theta2,theta3, origin, atomsout)
    implicit none
    double precision,intent(in)::     atomsin(:,:),theta1, theta2, theta3 , origin(ndim)
    double precision, intent(out):: atomsout(:,:)
    double precision::     workvec(3)
    integer::              i,j,k, atom1, atom2, atom3

    atomsout(:,:)=0.0d0
    workvec(:)=0.0d0
    if (ndim .ne. 3) then
       write(*,*) "Wrong number of dimensions; change align_atoms subroutine!"
       stop
    end if

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
    
    return
  end subroutine align_atoms

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

  subroutine rotate_atoms(atoms,axis,theta)
    implicit none
    double precision::     rotmatrix(3,3), atoms(:,:)
    double precision::     theta
    integer::              i, axis, j,k

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
  !spline algorithms
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

  subroutine instanton(xtilde,a,b)
    implicit none
    integer::                        iprint, m, iflag, mp,idof, maxiter
    integer::                        i, lp, count, iw, j,k, dof
    double precision::               eps2, xtol, gtol, stpmin, stpmax
    double precision::               f, xtilde(:,:,:), xtemp(ndim,natom)
    double precision::               factr, a(:,:),b(:,:), com(ndim)
    double precision, allocatable::  fprime(:,:,:), work(:), fprimework(:)
    double precision, allocatable::  lb(:), ub(:), dsave(:), xwork(:)
    integer, allocatable::           nbd(:), iwork(:), isave(:)
    logical::                        lsave(4)
    character(len=60)::              task, csave

    dof= n*ndim*natom
    allocate(lb(dof), ub(dof),fprime(n,ndim,natom), nbd(dof))
    allocate(fprimework(dof), xwork(dof))
    ! call centreofmass(a, com)
    ! do i=1,ndim
    !    do j=1,natom
    !       a(i,j)= a(i,j)- com(i)
    !    end do
    ! end do
    ! call centreofmass(b, com)
    ! do i=1,ndim
    !    do j=1,natom
    !       b(i,j)= b(i,j)- com(i)
    !    end do
    ! end do 
    do i=1, n, 1
       ! call centreofmass(xtilde(i,:,:), com)
       do j=1,ndim
          do k=1,natom
             ! xtilde(i,j,k)= xtilde(i,j,k) - com(j)
             idof= ((k-1)*ndim + j -1)*n +i
                lb(idof)= a(j,k)
                ub(idof)= a(j,k)
                nbd(idof)=0
          end do
       end do
    end do
      
    !------------------------
    !perform minimization
    task='START'
    m=8
    iprint=-1
    xtol= 1d-8
    iw=dof*(2*m+5) + 11*m**2 + 8*m
    allocate(work(iw), iwork(3*dof), isave(44), dsave(29))
    iflag=0
    eps2= 1.0d-5 !gradient convergence
    factr=1.0d5
    maxiter=40
    f= UM(xtilde,a,b)
    call UMprime(xtilde,a,b,fprime)
    count=0
    do while( task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
         task.eq.'START')
       count=count+1
       xwork=reshape(xtilde,(/dof/))
       fprimework= reshape(fprime,(/dof/))
       call setulb(dof,m,xwork,lb,ub,nbd,f,fprimework,factr,eps2,work&
            ,iwork,task,iprint, csave,lsave,isave,dsave,maxiter)
       if (task(1:2) .eq. 'FG') then
          xtilde= reshape(xwork,(/n,ndim,natom/))
          f= UM(xtilde,a,b)
          call UMprime(xtilde,a,b,fprime)
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
  subroutine detJ(x, etasquared)
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol
  double precision::               x(:,:,:), etasquared(:)
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
  call UMhessian(x,H)
  write(*,*) "Hessian is cooked."

  call DSBEVD('N', 'L', totdof, ndof, H, ndof+1, etasquared, z, 1, work, lwork, iwork, liwork, info)
  deallocate(work, iwork, H)
  return
  end subroutine detJ

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

end module instantonmod
