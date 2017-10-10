module mcmod_mass
  use malonaldehyde
  implicit none
  double precision, parameter::    pi=3.14159265358979d0
  double precision::               beta, betan, UMtilde
  integer::                        n, ndim, ndof, natom, xunit, totdof
  double precision, allocatable::  well1(:,:), well2(:,:), mass(:)

  public :: QsortC
  private :: Partition
contains
  !---------------------------------------------------------------------
  subroutine V_init()
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine detJ(xtilde, etasquared)
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol
  double precision::               xtilde(:,:,:), etasquared(:)
  integer::                        nout, ldz, lwork, liwork, info,i
  integer,allocatable::            isuppz(:), iwork(:)
  double precision, allocatable::  work(:), z(:,:), H(:,:)
  ! !get diagonal hessian
  jobz='N'
  range='A'
  uplo='U'
  abstol=1.0d-8
  lwork= 2*totdof +1!26*totdof !
  liwork= 1!10*totdof !
  info=0
  ldz=totdof
  vl=0.0d0
  vu=0.0d0
  nout=0
  allocate(isuppz(2*totdof), work(lwork), iwork(liwork), z(totdof,totdof), H(totdof,totdof))
  H=0.0d0
  etasquared=0.0d0
  call UMhessian(xtilde,H)
  ! call dsyevr(jobz, range, uplo, totdof, H, totdof, vl, vu, totdof, totdof, abstol, nout, etasquared,&
  !      z, ldz, isuppz, work, lwork, iwork, liwork, info)
  call dsyevd(jobz, uplo, totdof, H, totdof, etasquared,work,lwork,iwork,liwork,info)
  ! write(*,*) info
  ! do i=1,totdof
  ! write(*,*) i,etasquared(i)
  ! end do

  deallocate(isuppz, work, iwork, z,H)
  return
  end subroutine detJ
  !---------------------------------------------------------------------
  function V(x)
    implicit none
    double precision::     v, x(:,:)
    double precision, allocatable:: dummy1(:),dummy2(:)
    
    call pes(x,0,V,dummy1,dummy2)
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
          ! write(*,*) i1,j1,i2,j2,idof1, idof2,ij,hesstemp(ij)
       end do
    end do
    deallocate(hesstemp, dummy2)
    return
  end subroutine Vdoubleprime
  !---------------------------------------------------------------------
  function UM(x,a,b)
    implicit none
    integer::            i,j,k
    double precision::   x(:,:,:), UM,a(:,:),b(:,:)

    UM=0.0d0
    do i=2, N-1, 1
       UM=UM+ V(x(i,:,:))
       do j=1, ndim
          do k=1, natom
             UM=UM+ (0.5d0*mass(k)/betan**2)*(x(i+1,j,k)-x(i,j,k))**2
          end do
       end do
    end do
    UM=UM+ V(x(n,:,:))+ V(x(1,:,:))
    do j=1, ndim
       do k=1, natom
          UM=UM+ (0.5d0*mass(k)/betan**2)*(x(1,j,k)-a(j,k))**2
          UM=UM+ (0.5d0*mass(k)/betan**2)*(b(j,k)-x(N,j,k))**2
       end do
    end do

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
                answer(1,j,k)=mass(k)*(2.0*x(1,j,k) - a(j,k) - x(2,j,k))/betan**2
             else if (i.eq.N) then
                answer(N,j,k)=mass(k)*(2.0*x(N,j,k) - x(N-1,j,k) - b(j,k))/betan**2
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
                   idof1= (k1-1)*ndim + j1
                   idof2= (k2-1)*ndim + j2
                   answer(n*(idof1-1) + i,n*(idof2-1) + i)= hess(j1,k1,j2,k2)/sqrt(mass(k1)*mass(k2))
                   if (idof1.eq.idof2) then
                      answer(n*(idof1-1) + i,n*(idof1-1) + i)= &
                           answer(n*(idof1-1) + i,n*(idof1-1) + i) +2.0d0/betan**2
                      if (i.gt.1) answer(n*(idof1-1) + i,n*(idof1-1) + i-1)=&
                           answer(n*(idof1-1) + i,n*(idof1-1) + i-1) -1.0d0/betan**2
                      if (i.lt.n) answer(n*(idof1-1) + i,n*(idof1-1) + i+1)=&
                           answer(n*(idof1-1) + i,n*(idof1-1) + i+1)-1.0d0/betan**2
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
    double precision::   x1(:,:), x2(:,:), eucliddist
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

  subroutine instanton(xtilde,a,b)
    implicit none
    integer::                        iprint, m, iflag, mp,idof
    integer::                        i, lp, count, iw, j,k, dof
    double precision::               eps, xtol, gtol, stpmin, stpmax
    double precision::               f, xtilde(:,:,:)
    double precision::               factr, a(:,:),b(:,:)
    double precision, allocatable::  fprime(:,:,:), work(:), fprimework(:)
    double precision, allocatable::  lb(:), ub(:), dsave(:), xwork(:)
    integer, allocatable::           nbd(:), iwork(:), isave(:)
    logical::                        lsave(4)
    character(len=60)::              task, csave

    dof= n*ndim*natom
    allocate(lb(dof), ub(dof),fprime(n,ndim,natom), nbd(dof))
    allocate(fprimework(dof), xwork(dof))
    do i=1, n, 1
       do j=1,ndim
          do k=1,natom
             idof= ((k-1)*ndim + j -1)*n +i
             lb(idof)= a(j,k)
             ub(idof)= b(j,k)
             nbd(idof)=0
             xtilde(i,j,k)= a(j,k) + (b(j,k)-a(j,k))*dble(i-1)/dble(n-1)
          end do
       end do
    end do
      
    !------------------------
    !perform minimization

    task='START'
    m=10
    iprint=-1
    xtol= 1d-8
    iw=dof*(2*m+5) + 11*m**2 + 8*m
    allocate(work(iw), iwork(3*dof), isave(44), dsave(29))
    iflag=0
    eps= 1.0d-8
    factr=1.0d4
    f= UM(xtilde,a,b)
    call UMprime(xtilde,a,b,fprime)
    count=0
    do while( task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
         task.eq.'START')
       ! write(*,*) f
       count=count+1
       xwork=reshape(xtilde,(/dof/))
       fprimework= reshape(fprime,(/dof/))
       call setulb(dof,m,xwork,lb,ub,nbd,f,fprimework,factr,eps,work&
            ,iwork,task,iprint, csave,lsave,isave,dsave)
       if (task(1:2) .eq. 'FG') then
          xtilde= reshape(xwork,(/n,ndim,natom/))
          f= UM(xtilde,a,b)
          call UMprime(xtilde,a,b,fprime)
          ! write(*,*) count, f, dot_product(fprime, fprime)
       end if
    end do
    if (task(1:5) .eq. "ERROR" .or. task(1:4) .eq. "ABNO") then
         write(*,*) "Error:"
         write(*,*) task
         ! stop
      end if
    deallocate(work, lb, ub, fprime, fprimework,xwork)
    deallocate(iwork, nbd, isave, dsave)

    return
  end subroutine instanton

end module mcmod_mass
