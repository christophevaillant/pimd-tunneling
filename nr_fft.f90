module nr_fft
  implicit none
  !Module of subroutines from Numerical Recipes
INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
! Symbolic names for kind types of single- and double-precision reals:
INTEGER, PARAMETER :: SP = KIND(1.0)
INTEGER, PARAMETER :: DP = KIND(1.0D0)
! Symbolic names for kind types of single- and double-precision complex:
INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
! Symbolic name for kind type of default logical:
INTEGER, PARAMETER :: LGT = KIND(.true.)
! Frequently used mathematical constants (with precision to spare):
REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp

INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8

  INTERFACE assert
     MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
  END INTERFACE assert

  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE assert_eq

  INTERFACE cumsum
     MODULE PROCEDURE cumsum_r,cumsum_i
  END INTERFACE cumsum

  INTERFACE arth
     MODULE PROCEDURE arth_i
  END INTERFACE arth

  INTERFACE swap
     MODULE PROCEDURE swap_i,swap_r,swap_rv, &
          swap_z,swap_zv,swap_zm, &
          masked_swap_rs,masked_swap_rv,masked_swap_rm
  END INTERFACE swap

contains

  SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rv

  SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rm

  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(SP), INTENT(INOUT) :: a,b
    LOGICAL(LGT), INTENT(IN) :: mask
    REAL(SP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_rs

  SUBROUTINE swap_r(a,b)
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r
  SUBROUTINE swap_rv(a,b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv
  SUBROUTINE swap_z(a,b)
    COMPLEX(DPC), INTENT(INOUT) :: a,b
    COMPLEX(DPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_z
  SUBROUTINE swap_zv(a,b)
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zv
  SUBROUTINE swap_zm(a,b)
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zm

  SUBROUTINE swap_i(a,b)
    INTEGER(I4B), INTENT(INOUT) :: a,b
    INTEGER(I4B) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i

  SUBROUTINE assert1(n1,string)
    ! Report and die if any logical is false (used for arg range checking).
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1
    if (.not. n1) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert1'
    end if
  END SUBROUTINE assert1
  SUBROUTINE assert2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2
    if (.not. (n1 .and. n2)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert2'
    end if
  END SUBROUTINE assert2
  SUBROUTINE assert3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3
    if (.not. (n1 .and. n2 .and. n3)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert3'
    end if
  END SUBROUTINE assert3
  SUBROUTINE assert4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3,n4
    if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert4'
    end if
  END SUBROUTINE assert4
  SUBROUTINE assert_v(n,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, DIMENSION(:), INTENT(IN) :: n
    if (.not. all(n)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert_v'
    end if
  END SUBROUTINE assert_v


  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2
  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3
  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4
  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn

  subroutine sinft(y)
    implicit none
    double precision, dimension(:), intent(inout) :: y
    double precision, dimension(size(y)/2+1) :: wi
    double precision, dimension(size(y)/2) :: y1,y2
    integer(i4b) :: n,nh
    n=size(y)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in sinft')
    nh=n/2
    wi=aimag(zroots_unity(n+n,nh+1))
    y(1)=0.0d0
    y1=wi(2:nh+1)*(y(2:nh+1)+y(n:nh+1:-1))
    y2=0.5d0*(y(2:nh+1)-y(n:nh+1:-1))
    y(2:nh+1)=y1+y2
    y(n:nh+1:-1)=y1-y2
    call realft(y,+1)
    y(1)=0.5d0*y(1)
    y(2)=0.0d0
    y1=cumsum(y(1:n-1:2))
    y(1:n-1:2)=y(2:n:2)
    y(2:n:2)=y1
  end subroutine sinft

  SUBROUTINE realft(data,isign,zdata)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(DPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
    INTEGER(I4B) :: n,ndum,nh,nq
    COMPLEX(DPC), DIMENSION(size(data)/4) :: w
    COMPLEX(DPC), DIMENSION(size(data)/4-1) :: h1,h2
    COMPLEX(DPC), DIMENSION(:), POINTER :: cdata
    COMPLEX(DPC) :: z
    REAL(DP) :: d1=0.5_dp,d2
    n=size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_dp')
    nh=n/2
    nq=n/4
    if (present(zdata)) then
       ndum=assert_eq(n/2,size(zdata),'realft_dp')
       cdata=>zdata
       if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=dpc)
    else
       allocate(cdata(n/2))
       cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=dpc)
    end if
    if (isign == 1) then
       d2=-0.5_dp
       call four1(cdata,+1)
    else
       d2=0.5_dp
    end if
    w=zroots_unity(sign(n,isign),n/4)
    w=cmplx(-aimag(w),real(w),kind=dpc)
    h1=d1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
    h2=d2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
    cdata(2:nq)=h1+w(2:nq)*h2
    cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
    z=cdata(1)
    if (isign == 1) then
       cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=dpc)
    else
       cdata(1)=cmplx(d1*(real(z)+aimag(z)),d1*(real(z)-aimag(z)),kind=dpc)
       call four1(cdata,-1)
    end if
    if (present(zdata)) then
       if (isign /= 1) then
          data(1:n-1:2)=real(cdata)
          data(2:n:2)=aimag(cdata)
       end if
    else
       data(1:n-1:2)=real(cdata)
       data(2:n:2)=aimag(cdata)
       deallocate(cdata)
    end if
  END SUBROUTINE realft

  SUBROUTINE four1(data,isign)
    IMPLICIT NONE
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
    REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
    INTEGER(I4B) :: n,m1,m2,j
    n=size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_dp')
    m1=2**ceiling(0.5d0*log(real(n,sp))/0.693147d0)
    m2=n/m1
    allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat=reshape(data,shape(dat))
    call fourrow(dat,isign)
    theta=arth_i(0,isign,m1)*TWOPI_D/n
    wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
    w=cmplx(1.0_dp,0.0_dp,kind=dpc)
    do j=2,m2
       w=w*wp+w
       dat(:,j)=dat(:,j)*w
    end do
    temp=transpose(dat)
    call fourrow(temp,isign)
    data=reshape(temp,shape(data))
    deallocate(dat,w,wp,theta,temp)
  END SUBROUTINE four1

  SUBROUTINE fourrow(data,isign)
    IMPLICIT NONE
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
    REAL(DP) :: theta
    COMPLEX(DPC), DIMENSION(size(data,1)) :: temp
    COMPLEX(DPC) :: w,wp
    COMPLEX(DPC) :: ws
    n=size(data,2)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_dp')
    n2=n/2
    j=n2
    do i=1,n-2
       if (j > i) call swap_zv(data(:,j+1),data(:,i+1))
       m=n2
       do
          if (m < 2 .or. j < m) exit
          j=j-m
          m=m/2
       end do
       j=j+m
    end do
    mmax=1
    do
       if (n <= mmax) exit
       istep=2*mmax
       theta=PI_D/(isign*mmax)
       wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
       w=cmplx(1.0_dp,0.0_dp,kind=dpc)
       do m=1,mmax
          ws=w
          do i=m,n,istep
             j=i+mmax
             temp=ws*data(:,j)
             data(:,j)=data(:,i)-temp
             data(:,i)=data(:,i)+temp
          end do
          w=w*wp+w
       end do
       mmax=istep
    end do
  END SUBROUTINE fourrow

  RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    REAL(DP), OPTIONAL, INTENT(IN) :: seed
    REAL(DP), DIMENSION(size(arr)) :: ans
    INTEGER(I4B) :: n,j
    REAL(DP) :: sd
    n=size(arr)
    if (n == 0_i4b) RETURN
    sd=0.0d0
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  END FUNCTION cumsum_r

  RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: seed
    INTEGER(I4B), DIMENSION(size(arr)) :: ans
    INTEGER(I4B) :: n,j,sd
    n=size(arr)
    if (n == 0_i4b) RETURN
    sd=0_i4b
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  END FUNCTION cumsum_i

  FUNCTION zroots_unity(n,nn)
    INTEGER(I4B), INTENT(IN) :: n,nn
    COMPLEX(DPC), DIMENSION(nn) :: zroots_unity
    INTEGER(I4B) :: k
    REAL(DP) :: theta
    zroots_unity(1)=1.0
    theta=TWOPI/n
    k=1
    do
       if (k >= nn) exit
       zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),DPC)
       zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
            zroots_unity(2:min(k,nn-k))
       k=2*k
    end do
  END FUNCTION zroots_unity


    FUNCTION arth_i(first,increment,n)
      INTEGER(I4B), INTENT(IN) :: first,increment,n
      INTEGER(I4B), DIMENSION(n) :: arth_i
      INTEGER(I4B) :: k,k2,temp
      if (n > 0) arth_i(1)=first
      if (n <= NPAR_ARTH) then
         do k=2,n
            arth_i(k)=arth_i(k-1)+increment
         end do
      else
         do k=2,NPAR2_ARTH
            arth_i(k)=arth_i(k-1)+increment
         end do
         temp=increment*NPAR2_ARTH
         k=NPAR2_ARTH
         do
            if (k >= n) exit
            k2=k+k
            arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
            temp=temp+temp
            k=k2
         end do
      end if
    END FUNCTION arth_i

end module nr_fft
