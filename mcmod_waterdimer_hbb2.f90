module mcmod_mass
  implicit none
  double precision::               V0
  integer,parameter::              atom1=1, atom2=2, atom3=3
  integer::                        nw
  integer,dimension(:,:),allocatable::conn
  integer::                        n, ndim, ndof, natom, xunit, totdof

contains

  subroutine V_init()
    use pes, only: pes0_init,pes1_init
    integer::i
    integer,dimension(nw)::idx_o
    character (len=*), parameter:: pes_x6y3_sysall=&
         'x1 y1 x2 x1y1 y2 x3 x2y1 x1y2 y3 x4 x3y1 x2y2 x1y3 '// &
         'x5 x4y1 x3y2 x2y3 x6 x5y1 x4y2 x3y3 x6y1 x5y2 x4y3 x5y2 x4y3 '// &
         'x6y2 x5y3 x6y3'
    
    nw=2

    ! 3-body pot
    call pes0_init (dir='./coef-3b')
    call pes1_init (pes_x6y3_sysall)

    ! dimer pot
    call prepot()

    ! H-bonded 3b
    allocate(conn(nw*2,nw-1))
    idx_o=(/(i,i=1,nw)/)
    do i=1,nw
       idx_o(i)=0
       conn(i*2-1,:)=pack(idx_o,mask=idx_o.ne.0)
       conn(i*2,:)=conn(i*2-1,:)
       idx_o(i)=i
    end do

    write(*,*) "Potential initializaton complete"

    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x)
  use pes,only:pes_x6y3_pot
    implicit none
    double precision::     v, x(:,:)
    double precision,dimension(6,3)::x1,x2
    double precision,dimension(3,9)::x3
    double precision::e3

    double precision::e1,e2,vect(3)
    integer::              i,j, k,fo

    V=0.0d0
    fo=natom/3*2
    do i=1,natom/3-1
       x2(1,:)=x(:,fo+i)        ! O1
       x2(2,:)=x(:,i*2-1)       ! H1
       x2(3,:)=x(:,i*2)         ! H2
       do j=i+1,natom/3
          x2(4,:)=x(:,fo+j)     ! O2
          x2(5,:)=x(:,j*2-1)    ! H1
          x2(6,:)=x(:,j*2)      ! H2
          vect(:)=x2(4,:)-x2(1,:)
          x1=x2
          vect(:)=vect(:)*200.d0
          x1(4,:)=x2(4,:)+vect(:)
          x1(5,:)=x2(5,:)+vect(:)
          x1(6,:)=x2(6,:)+vect(:)
          call calcpot(e2,x2) 
          call calcpot(e1,x1) 
          V=V+e2-e1*(1.d0-1.d0/dble(natom/3-1))
       end do
    end do

    e3=0.d0
    do i=1,natom/3-2
       x3(:,7)=x(:,fo+i)          ! O1
       x3(:,1)=x(:,i*2-1)         ! H1
       x3(:,2)=x(:,i*2)           ! H1'
       do j=i+1,natom/3-1
          x3(:,8)=x(:,fo+j)       ! O2
          x3(:,3)=x(:,j*2-1)      ! H2
          x3(:,4)=x(:,j*2)        ! H2'
          do k=j+1,natom/3
             x3(:,9)=x(:,fo+k)       ! O3
             x3(:,5)=x(:,k*2-1)      ! H3
             x3(:,6)=x(:,k*2)        ! H3'
             e3=pes_x6y3_pot(x3)
             V=V+e3
          end do
       end do
    end do
    ! V= V
    return
  end function V

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision::     grad(:,:), x(:,:)
    double precision::     potplus, potminus

    eps=1d-4
    do i= 1,ndim
       do j=1,natom
          x(i,j)= x(i,j) + eps
          potplus= V(x)
          x(i,j)= x(i,j) - 2.0d0*eps
          potminus= V(x)
          x(i,j)= x(i,j) + eps
          grad(i,j)= (potplus-potminus)/(2.0d0*eps)
       end do
    end do
    return
  end subroutine Vprime
  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1
    integer::              i, j
    double precision::     gradplus(ndim, natom), gradminus(ndim, natom)

    eps=1d-4
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
