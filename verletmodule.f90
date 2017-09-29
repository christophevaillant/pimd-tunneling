include 'mkl_vsl.fi'
module verletint
  use mcmod_mass
  use MKL_VSL_TYPE
  use MKL_VSL
  use nr_fft, only : sinft

  implicit none

  integer::                         NMC, Noutput, imin, iproc
  integer, allocatable::            ipar(:)
  double precision::                dt
  double precision,allocatable::    transmatrix(:,:),dpar(:), beadvec(:,:)
  double precision,allocatable::    beadmass(:,:), lam(:)
  double precision, allocatable::   c1(:,:), c2(:,:)
  double precision::                alpha1, alpha2, alpha3, alpha4
  double precision::                beta1, beta2, beta3, tau, gamma
  logical::                         iprint, use_mkl
  integer::                         errcode_poisson, rmethod_poisson
  integer::                         brng_poisson, seed_poisson, stat
  type (vsl_stream_state)::         stream_poisson,stream_normal
  integer::                         errcode_normal, rmethod_normal, brng_normal
  integer::                         seed_normal, rk

contains

  subroutine gauleg(x1,x2,x,w,nintegral)
    implicit none
    integer::             i,j,m,nintegral
    double precision::    x1,x2,x(:),w(:)
    double precision::    EPS
    parameter (EPS=3.d-14)
    double precision      p1,p2,p3,pp,xl,xm,z,z1
    ! EPS is the relative precision.
    ! Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
    ! arrays x(1:n) and w(1:n) of length n , containing the abscissas and weights of the Gauss-
    ! Legendre n -point quadrature formula.
    ! High precision is a good idea for this routine.
    m=(nintegral+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    do i=1,m
       ! Loop over the desired roots.
       z=cos(pi*(i-0.25d0)/(nintegral+0.5d0))
1       continue
       p1=1.d0
       p2=0.d0
       do j=1,nintegral
          p3=p2
          p2=p1
          p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
       enddo
       pp=nintegral*(z*p1-p2)/(z*z-1.d0)
       z1=z
       z=z1-p1/pp
       if(abs(z-z1).gt.EPS) goto 1
       x(i)=xm-xl*z
       x(nintegral+1-i)=xm+xl*z
       w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
       w(nintegral+1-i)=w(i)
    enddo
    return
  end subroutine gauleg

  !-----------------------------------------------------
  !-----------------------------------------------------
  !main function for propagating a MD trajectory
  subroutine propagate_pimd_nm(xprop,vprop,a,b,dbdl,dHdr)
    implicit none
    double precision::     xprop(:,:,:), vprop(:,:,:),randno, sigma
    double precision::     a(:,:),b(:,:), dHdr, totenergy, kin,dbdl(:,:)
    double precision, allocatable:: pprop(:), tempp(:), tempv(:)
    integer::              i,j,k,l,count, time1,time2,imax, irate, dofi
    integer (kind=4)::     rkick(1)

    allocate(pprop(n), tempp(n),tempv(n))
    count=0
    dHdr=0.0d0
    kin=0.0d0
    errcode_poisson = virngpoisson( rmethod_poisson, stream_poisson, 1, rkick(1), dble(Noutput))
    ! open(20, file="andersen_m.dat")
    do i=1, NMC, 1
       count=count+1
       if (count .ge. rkick(1)) then
          count=0
          if (iprint) write(*,*) 100*dble(i)/dble(NMC),0.5d0*dHdr/(betan*dble(i-imin))!dHdr/dble(i)
          do j=1,ndim
             do k=1,natom
                dofi=(k-1)*ndim +j
                sigma=sqrt(1.0d0/betan)
                errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,pprop,0.0d0,sigma)
                do l=1,n
                   ! tempp(l)= pprop((dofi-1)*n +l)*sqrt(beadmass(k,l))
                   tempp(l)= pprop(l)*sqrt(beadmass(k,l))
                end do
                if (.not. use_mkl) then
                   call nmtransform_backward(tempp, tempv, 0)
                else
                   call nmtransform_backward_nr(tempp, tempv,0)
                end if
                do l=1,n
                   vprop(l, j, k) = tempv(l)
                end do
             end do
          end do
          errcode_poisson = virngpoisson( rmethod_poisson, stream_poisson, 1, rkick(1), dble(Noutput))
          ! write(*,*)"-----------------"
          count=0
          ! write(20,*)i, dHdr/dble(i)
       end if
       call time_step_test(xprop, vprop)
       if (i .gt. imin) then
       do j=1,ndim
          do k=1,natom
             dHdr= dHdr+mass(k)*(-xprop(n,j,k))*dbdl(j,k)
          end do
       end do
    end if
    end do
    ! close(20)
    ! stop
    dHdr=dHdr/dble(NMC-imin)
    deallocate(pprop, tempv, tempp)
    return
  end subroutine propagate_pimd_nm
  !-----------------------------------------------------
  !-----------------------------------------------------
  !normal mode forward transformation for linear polymer
  subroutine nmtransform_forward(xprop, qprop, bead)
    implicit none
    double precision::        xprop(:), qprop(:)
    integer::                 i,j,bead
    qprop(:)=xprop(:)
    call dsymv('U', n, 1.0d0,transmatrix, n, xprop,1,0.0d0, qprop,1)
    if (bead .gt. 0) then
       do i=1,n
          qprop(i)= qprop(i) - beadvec(i,bead)
       end do
    end if
    return
  end subroutine nmtransform_forward
  !-----------------------------------------------------
  !-----------------------------------------------------
  !normal mode backward transformation for linear polymer
  subroutine nmtransform_backward(qprop, xprop,bead)
    implicit none
    double precision::        xprop(:), qprop(:)
    integer::                 i,j,bead
    if (bead .gt. 0) then
       do i=1,n
          qprop(i)= qprop(i) + beadvec(i,bead)
       end do
    end if
    call dsymv('U', n, 1.0d0, transmatrix, n, qprop,1,0.0d0, xprop,1)
    if (bead .gt.0) then
       do i=1,n
          qprop(i)= qprop(i) - beadvec(i,bead)
       end do
    end if
    return
  end subroutine nmtransform_backward
  !-----------------------------------------------------
  !-----------------------------------------------------
  !routine taking a step in time using normal mode verlet algorithm
  !from ceriotti et al 2010 paper.
  subroutine time_step_nm(x, v)
    implicit none
    double precision::    x(:,:,:), v(:,:,:), omegak
    double precision, allocatable::  newv(:,:,:), newx(:,:,:),force(:,:)
    double precision, allocatable::  q(:,:,:), p(:,:,:)
    integer::             i,j,k, dofi

    allocate(newv(n,ndim,natom), newx(n,ndim,natom),force(ndim,natom))
    allocate(q(n,ndim,natom),p(n,ndim,natom))
    newv(:,:,:)=0.0d0
    newx(:,:,:)=0.0d0
    force(:,:)=0.0d0
    q(:,:,:)=0.0d0
    p(:,:,:)=0.0d0
    !-------------
    !step 2
    do i=1,ndim
       do j=1,natom
          dofi=(j-1)*ndim + i
          if (.not. use_mkl) then
             call nmtransform_forward(v(:,i,j), p(:,i,j), 0)
             call nmtransform_forward(x(:,i,j), q(:,i,j), dofi)
          else
             call nmtransform_forward_nr(v(:,i,j), p(:,i,j), 0)
             call nmtransform_forward_nr(x(:,i,j), q(:,i,j), dofi)
          end if
       end do
    end do
    !-------------
    !step 3
    do i=1, n
       do j=1,ndim
          do k=1,natom
             omegak= sqrt(mass(k)/beadmass(k,i))*lam(i)
             newv(i,j,k)= p(i,j,k)*cos(0.5d0*dt*omegak) - &
                  q(i,j,k)*omegak*beadmass(k,i)*sin(0.5d0*omegak*dt)
             newx(i,j,k)= q(i,j,k)*cos(0.5d0*dt*omegak) + &
                  p(i,j,k)*sin(0.5d0*omegak*dt)/(omegak*beadmass(k,i))
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in 1st NM propagation"
                write(*,*) omegak, mass(k), beadmass(k,i), lam(i)
                write(*,*) p(i,j,k), v(i,j,k)
                stop
             end if
          end do
       end do
    end do
    !-------------
    !step 4
    do i=1,ndim
       do j=1,natom
          dofi=(j-1)*ndim + i
          if (.not. use_mkl) then
             call nmtransform_backward(newv(:,i,j), v(:,i,j), 0)
             call nmtransform_backward(newx(:,i,j), x(:,i,j),dofi)
          else
             call nmtransform_backward_nr(newv(:,i,j), v(:,i,j),0)
             call nmtransform_backward_nr(newx(:,i,j), x(:,i,j),dofi)
          end if
       end do
    end do
    !-------------
    !step 1
    do i=1,n
       call Vprime(x(i,:,:),force)
       do j=1,ndim
          do k=1,natom
             newv(i,j,k)= v(i,j,k) - force(j,k)*dt
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in pot propagation"
                stop
             end if
          end do
       end do
    end do
    !-------------
    !step 2
    do i=1,ndim
       do j=1,natom
          if (.not. use_mkl) then
             call nmtransform_forward(newv(:,i,j), p(:,i,j), 0)
          else
             call nmtransform_forward_nr(newv(:,i,j), p(:,i,j), 0)
          end if
          q(:,i,j)=newx(:,i,j)
       end do
    end do
    !-------------
    !step 3
    do i=1, n
       do j=1,ndim
          do k=1,natom
             omegak= sqrt(mass(k)/beadmass(k,i))*lam(i)
             newv(i,j,k)= p(i,j,k)*cos(0.5d0*dt*omegak) - &
                  q(i,j,k)*omegak*beadmass(k,i)*sin(0.5d0*omegak*dt)
             newx(i,j,k)= q(i,j,k)*cos(0.5d0*dt*omegak) + &
                  p(i,j,k)*sin(0.5d0*omegak*dt)/(omegak*beadmass(k,i))
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in 1st NM propagation"
                stop
             end if
          end do
       end do
    end do
    !-------------
    !step 4
    do i=1,ndim
       do j=1,natom
          dofi=(j-1)*ndim + i
          if (.not. use_mkl) then
             call nmtransform_backward(newv(:,i,j), v(:,i,j), 0)
             call nmtransform_backward(newx(:,i,j), x(:,i,j),dofi)
          else
             call nmtransform_backward_nr(newv(:,i,j), v(:,i,j),0)
             call nmtransform_backward_nr(newx(:,i,j), x(:,i,j),dofi)
          end if
       end do
    end do
    deallocate(newv, newx,force)
    deallocate(q,p)
    return
  end subroutine time_step_nm

  !-----------------------------------------------------
  !-----------------------------------------------------
  !routine taking a step in time using normal mode verlet algorithm
  !from ceriotti et al 2010 paper.
  subroutine time_step_test(x, v)
    implicit none
    double precision::    x(:,:,:), v(:,:,:), omegak
    double precision, allocatable::  newv(:,:,:), newx(:,:,:),force(:,:)
    double precision, allocatable::  q(:,:,:), p(:,:,:)
    integer::             i,j,k,dofi

    allocate(newv(n,ndim,natom), newx(n,ndim,natom),force(ndim,natom))
    allocate(q(n,ndim,natom),p(n,ndim,natom))
    !-------------
    !step 1
    do i=1,n
       call Vprime(x(i,:,:),force)
       do j=1,ndim
          do k=1,natom
             newv(i,j,k)= v(i,j,k) - force(j,k)*dt
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in pot propagation"
                stop
             end if
          end do
       end do
    end do
    !-------------
    !step 2
    do i=1,ndim
       do j=1,natom
          if (.not. use_mkl) then
             dofi= (j-1)*ndim+i
             call nmtransform_forward(newv(:,i,j), p(:,i,j), 0)
             call nmtransform_forward(x(:,i,j), q(:,i,j), dofi)
          else
             call nmtransform_forward_nr(newv(:,i,j), p(:,i,j), 0)
             call nmtransform_forward_nr(x(:,i,j), q(:,i,j), dofi)
          end if
       end do
    end do
    !-------------
    !step 3
    do i=1, n
       do j=1,ndim
          do k=1,natom
             omegak= sqrt(mass(k)/beadmass(k,i))*lam(i)
             newv(i,j,k)= p(i,j,k)*cos(dt*omegak) - &
                  q(i,j,k)*omegak*beadmass(k,i)*sin(omegak*dt)
             newx(i,j,k)= q(i,j,k)*cos(dt*omegak) + &
                  p(i,j,k)*sin(omegak*dt)/(omegak*beadmass(k,i))
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in 1st NM propagation"
                stop
             end if
          end do
       end do
    end do
    !-------------
    !step 4
    do i=1,ndim
       do j=1,natom
          dofi=(j-1)*ndim + i
          if (.not. use_mkl) then
             call nmtransform_backward(newv(:,i,j), v(:,i,j), 0)
             call nmtransform_backward(newx(:,i,j), x(:,i,j),dofi)
          else
             call nmtransform_backward_nr(newv(:,i,j), v(:,i,j),0)
             call nmtransform_backward_nr(newx(:,i,j), x(:,i,j),dofi)
          end if
       end do
    end do
    !-------------
    !step 5
    ! do i=1,n
    !    call Vprime(x(i,:,:),force)
    !    do j=1,ndim
    !       do k=1,natom
    !          newv(i,j,k)= v(i,j,k) - 0.5d0*force(j,k)*dt
    !          if (newv(i,j,k) .ne. newv(i,j,k)) then
    !             write(*,*) "NaN in pot propagation"
    !             stop
    !          end if
    !       end do
    !    end do
    ! end do

    deallocate(newv, newx,force)
    deallocate(q,p)
    return
  end subroutine time_step_test
  !-----------------------------------------------------
  !-----------------------------------------------------
  !initialize normal mode routines
  subroutine init_nm(a,b)
    double precision::    a(:,:),b(:,:)
    integer::             i,j,k,l, dofi

    do i=1, n
       lam(i)= 2.0d0*sin(dble(i)*pi/dble(2*n +2))/betan
       do j=1, natom
          beadmass(j,i)=mass(j)*(lam(i)*tau)**2
       end do
       do l=i,n
          transmatrix(i,l)= sin(dble(i*l)*pi/dble(n+1))*sqrt(2.0d0/dble(n+1))
          transmatrix(l,i)= transmatrix(i,l)
          if (transmatrix(i,l) .ne. transmatrix(i,l)) then
             write(*,*)"Nan!"
             stop
          end if
       end do
       do j=1,ndim
          do k=1,natom
             dofi= (k-1)*ndim +j
             beadvec(i, dofi)= a(j,k)*sin(dble(i)*pi/dble(n+1)) &
                  + b(j,k)*sin(dble(n*i)*pi/dble(n+1))
             beadvec(i,dofi)= beadvec(i,dofi)*sqrt(2.0d0/dble(n+1))
             beadvec(i,dofi)= beadvec(i,dofi)/(lam(i)*betan)**2
          end do
       end do
    end do

  end subroutine init_nm
  !-----------------------------------------------------
  !-----------------------------------------------------
  !free normal mode arrays

  subroutine free_nm()
    deallocate(transmatrix,beadmass,beadvec,lam)
  end subroutine free_nm
  !-----------------------------------------------------
  !-----------------------------------------------------
  !allocate normal mode arrays

  subroutine alloc_nm(iproc)
    implicit none
    integer::    itime, irate, imax
    integer, optional:: iproc

    call system_clock(itime,irate,imax)
    seed_normal= mod(itime,1000)
    call system_clock(itime,irate,imax)
    seed_poisson= mod(itime,1000)
    brng_normal = VSL_BRNG_MT19937
    brng_poisson = VSL_BRNG_MT19937
    rmethod_normal = VSL_RNG_METHOD_GAUSSIAN_ICDF
    rmethod_poisson = VSL_RNG_METHOD_POISSON_POISNORM
    errcode_normal = vslnewstream( stream_normal, brng_normal, seed_normal )
    errcode_poisson = vslnewstream( stream_poisson, brng_poisson, seed_poisson )
    
    if (present(iproc)) then
       write(*,*)"proc", iproc, "running with seeds:", seed_normal, seed_poisson
    else
       write(*,*)"running with seeds:", seed_normal, seed_poisson
    end if

  end subroutine alloc_nm
  !-----------------------------------------------------
  !-----------------------------------------------------
  !main function for propagating a MD trajectory with Langevin thermostat
  subroutine propagate_pimd_pile(xprop,vprop,a,b,dbdl,dHdr)
    implicit none
    double precision::     xprop(:,:,:), vprop(:,:,:),randno
    double precision::     a(:,:),b(:,:), dHdr, dbdl(:,:)
    double precision, allocatable:: force(:,:,:)
    integer::              i,j,k,count, time1,time2,imax, irate
    integer (kind=4)::     rkick(1)

    allocate(c1(natom,n), c2(natom,n), force(n, ndim, natom))
    do i=1,n
       do j=1,natom
          c1(j,i)= exp(-gamma*dt*lam(i)*sqrt(mass(j)/beadmass(j,i)))!
          c2(j,i)= sqrt(1.0d0- c1(j,i)**2)
       end do
    end do
    count=0
    dHdr=0.0d0
    force(:,:,:)=0.0d0
    do i=1, NMC, 1
       count=count+1
       if (count .ge. Noutput .and. i .gt. imin) then
          count=0
          if (iprint) write(*,*) 100*dble(i)/dble(NMC-imin), 0.5d0*dHdr/(betan*dble(i-imin))
       end if
       call time_step_test_pile(xprop, vprop)
       if (i.gt.imin) then
          do j=1,ndim
             do k=1,natom
                dHdr= dHdr+mass(k)*(-xprop(n,j,k))*dbdl(j,k)
             end do
          end do
       end if
    end do
    dHdr= dHdr/dble(NMC-imin)
    deallocate(c1, c2)
    return
  end subroutine propagate_pimd_pile
  !-----------------------------------------------------
  !-----------------------------------------------------
  !main function for propagating a MD trajectory with Langevin thermostat
  subroutine propagate_pimd_higher(xprop,vprop,a,b,dbdl,dHdr)
    implicit none
    double precision::     xprop(:,:,:), vprop(:,:,:),randno
    double precision::     a(:,:),b(:,:), dHdr, dbdl(:,:)
    double precision, allocatable:: force(:,:,:)
    integer::              i,j,k,count, time1,time2,imax, irate
    integer (kind=4)::     rkick(1)

    allocate(c1(natom,n), c2(natom,n), force(n, ndim, natom))
    do i=1,n
       do j=1,natom
          c1(j,i)= exp(-dt*lam(i)*sqrt(mass(j)/beadmass(j,i)))!
          c2(j,i)= sqrt(1.0d0- c1(j,i)**2)
       end do
    end do
    beta1= 1.0d0/(2.0d0- 2.0d0**(1.0d0/3.0d0))
    alpha1=0.5d0*beta1
    alpha2=0.5d0*(1.0d0- 2.0d0**(1.0d0/3.0d0))*beta1
    beta2=-2.0d0**(1.0d0/3.0d0)*beta1
    alpha3=alpha2
    beta3=beta1
    alpha4=alpha1
    count=0
    dHdr=0.0d0
    force(:,:,:)=0.0d0
    do i=1, NMC, 1
       count=count+1
       if (count .ge. Noutput .and. i .gt. imin) then
          count=0
          if (iprint) write(*,*) 100*dble(i)/dble(NMC-imin), 0.5d0*dHdr/(betan*dble(i-imin))
       end if
       call time_step_higher(xprop, vprop)
       if (i.gt.imin) then
          do j=1,ndim
             do k=1,natom
                dHdr= dHdr+mass(k)*(-xprop(n,j,k))*dbdl(j,k)
             end do
          end do
       end if
    end do
    dHdr= dHdr/dble(NMC-imin)
    deallocate(c1, c2)
    return
  end subroutine propagate_pimd_higher
  !-----------------------------------------------------
  !-----------------------------------------------------
  !routine taking a step in time using normal mode verlet algorithm
  !modified from ceriotti et al 2010 paper to call pot once
  subroutine time_step_pile(x, vprop)
    implicit none
    double precision::    x(:,:,:), vprop(:,:,:), omegak
    double precision, allocatable::  newv(:,:,:), newx(:,:,:),force(:,:,:)
    double precision, allocatable::  q(:,:,:), p(:,:,:)
    double precision, allocatable::  pprop(:)
    integer::             i,j,k,dofi

    allocate(pprop(n*ndof), force(n,ndim, natom))
    allocate(p(n,ndim,natom),q(n,ndim,natom))
    allocate(newv(n,ndim,natom), newx(n,ndim,natom))
    !-------------
    !Thermostat step
    errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n*ndof,pprop,0.0d0,1.0d0)
    do i=1,ndim
       do j=1,natom
          dofi= (j-1)*ndim+i
          if (.not. use_mkl) then
             call nmtransform_forward(vprop(:,i,j), p(:,i,j), 0)
             call nmtransform_forward(x(:,i,j), q(:,i,j), dofi)
          else
             call nmtransform_forward_nr(vprop(:,i,j), p(:,i,j), 0)
             call nmtransform_forward_nr(x(:,i,j), q(:,i,j), dofi)
          end if
          do k=1,n
             p(k,i,j)= (c1(j,k)**2)*p(k,i,j) + &
                  sqrt(beadmass(j,k)/betan)*c2(j,k)*sqrt(1.0+c1(j,k)**2)*pprop((dofi-1)*n +k)
          end do
       end do
    end do
    
    do i=1, n
       do j=1,ndim
          do k=1,natom
             omegak= sqrt(mass(k)/beadmass(k,i))*lam(i)
             newv(i,j,k)= p(i,j,k)*cos(0.5d0*dt*omegak) - &
                  q(i,j,k)*(omegak*beadmass(k,i))*sin(0.5d0*omegak*dt)
             newx(i,j,k)= q(i,j,k)*cos(0.5d0*dt*omegak) + &
                  p(i,j,k)*sin(0.5d0*omegak*dt)/(omegak*beadmass(k,i))
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in 1st PILE propagation."
                stop
             end if
          end do
       end do
    end do

    do i=1,ndim
       do j=1,natom
          dofi= (j-1)*ndim +i
          if (.not. use_mkl) then
             call nmtransform_backward(newv(:,i,j),vprop(:,i,j), 0)
             call nmtransform_backward(newx(:,i,j), x(:,i,j),dofi)
          else
             call nmtransform_backward_nr(newv(:,i,j),vprop(:,i,j), 0)
             call nmtransform_backward_nr(newx(:,i,j), x(:,i,j),dofi)
          end if
       end do
    end do
    !-------------
    !step 1
    do i=1,n
       call Vprime(x(i,:,:),force(i,:,:))
       ! write(*,*) i, x(i, 1, 1), force(i, 1, 1), x(i,2,1),force(i, 2, 1)
       do j=1,ndim
          do k=1,natom
             newv(i,j,k)= vprop(i,j,k) - force(i,j,k)*dt 
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in pot propagation"
                ! stop
             end if
          end do
       end do
    end do
    !-------------
    !step 2
    do i=1,ndim
       do j=1,natom
          if (.not. use_mkl) then
             call nmtransform_forward(newv(:,i,j), p(:,i,j), 0)
          else
             call nmtransform_forward_nr(newv(:,i,j), p(:,i,j), 0)
          end if
          q(:,i,j)=newx(:,i,j)
       end do
    end do
    !-------------
    !step 3
    do i=1, n
       do j=1,ndim
          do k=1,natom
             omegak= sqrt(mass(k)/beadmass(k,i))*lam(i)
             newv(i,j,k)= p(i,j,k)*cos(0.5d0*dt*omegak) - &
                  q(i,j,k)*omegak*beadmass(k,i)*sin(0.5d0*omegak*dt)
             newx(i,j,k)= q(i,j,k)*cos(0.5d0*dt*omegak) + &
                  p(i,j,k)*sin(0.5d0*omegak*dt)/(omegak*beadmass(k,i))
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in 1st NM propagation"
                write(*,*), i,j,k, p(i,j,k), q(i,j,k)
                stop
             end if
          end do
       end do
    end do
    !-------------
    !step 4
    do i=1,ndim
       do j=1,natom
          dofi=(j-1)*ndim + i
          if (.not. use_mkl) then
             call nmtransform_backward(newv(:,i,j), vprop(:,i,j), 0)
             call nmtransform_backward(newx(:,i,j), x(:,i,j),dofi)
          else
             call nmtransform_backward_nr(newv(:,i,j), vprop(:,i,j),0)
             call nmtransform_backward_nr(newx(:,i,j), x(:,i,j),dofi)
          end if
       end do
    end do

    deallocate(pprop)
    deallocate(newv, newx,force)
    deallocate(q,p)
    return
  end subroutine time_step_pile
  !-----------------------------------------------------
  !-----------------------------------------------------
  !routine taking a step in time using normal mode verlet algorithm
  !from ceriotti et al 2010 paper.
  subroutine time_step_test_pile(x, vprop)
    implicit none
    double precision::    x(:,:,:), vprop(:,:,:), omegak
    double precision, allocatable::  newv(:,:,:), newx(:,:,:),force(:,:,:)
    double precision, allocatable::  q(:,:,:), p(:,:,:)
    double precision, allocatable::  pprop(:)
    integer::             i,j,k,dofi

    allocate(pprop(n*ndof), force(n,ndim, natom))
    allocate(p(n,ndim,natom),q(n,ndim,natom))
    allocate(newv(n,ndim,natom), newx(n,ndim,natom))
    do i=1,n
       call Vprime(x(i,:,:),force(i,:,:))
    end do
    do i=1,n
       do j=1,ndim
          do k=1,natom
             vprop(i,j,k)= vprop(i,j,k) - 0.5d0*force(i,j,k)*dt 
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in pot propagation"
                stop
             end if
          end do
       end do
    end do
    !-------------
    !Thermostat step
    errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n*ndof,pprop,0.0d0,1.0d0)
    do i=1,ndim
       do j=1,natom
          dofi= (j-1)*ndim+i
          if (.not. use_mkl) then
             call nmtransform_forward(vprop(:,i,j), p(:,i,j), 0)
             ! call nmtransform_forward(x(:,i,j), q(:,i,j), dofi)
          else
             call nmtransform_forward_nr(vprop(:,i,j), p(:,i,j), 0)
             ! call nmtransform_forward_nr(x(:,i,j), q(:,i,j), dofi)
          end if
          do k=1,n
             p(k,i,j)= (c1(j,k)**2)*p(k,i,j) + &
                  sqrt(beadmass(j,k)/betan)*c2(j,k)*sqrt(1.0+c1(j,k)**2)*pprop((dofi-1)*n +k)
          end do
       end do
    end do
    
    do i=1,ndim
       do j=1,natom
          dofi= (j-1)*ndim +i
          if (.not. use_mkl) then
             call nmtransform_backward(p(:,i,j),vprop(:,i,j), 0)
             ! call nmtransform_backward(newx(:,i,j), x(:,i,j),dofi)
          else
             call nmtransform_backward_nr(p(:,i,j),vprop(:,i,j), 0)
             ! call nmtransform_backward_nr(newx(:,i,j), x(:,i,j),dofi)
          end if
       end do
    end do
    !-------------
    !step 1
       ! write(*,*) i, x(i, 1, 1), force(i, 1, 1), x(i,2,1),force(i, 2, 1)
    do i=1,n
       do j=1,ndim
          do k=1,natom
             newv(i,j,k)= vprop(i,j,k) - 0.5d0*force(i,j,k)*dt 
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in pot propagation"
                stop
             end if
          end do
       end do
    end do
    !-------------
    !step 2
    do i=1,ndim
       do j=1,natom
          dofi= (j-1)*ndim +i
          if (.not. use_mkl) then
             call nmtransform_forward(newv(:,i,j), p(:,i,j), 0)
             call nmtransform_forward(x(:,i,j), q(:,i,j), dofi)
          else
             call nmtransform_forward_nr(newv(:,i,j), p(:,i,j), 0)
             call nmtransform_forward_nr(x(:,i,j), q(:,i,j), dofi)
          end if
          ! q(:,i,j)=newx(:,i,j)
       end do
    end do
    !-------------
    !step 3
    do i=1, n
       do j=1,ndim
          do k=1,natom
             omegak= sqrt(mass(k)/beadmass(k,i))*lam(i)
             newv(i,j,k)= p(i,j,k)*cos(dt*omegak) - &
                  q(i,j,k)*omegak*beadmass(k,i)*sin(omegak*dt)
             newx(i,j,k)= q(i,j,k)*cos(dt*omegak) + &
                  p(i,j,k)*sin(omegak*dt)/(omegak*beadmass(k,i))
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in 1st NM propagation"
                write(*,*), i,j,k, p(i,j,k), q(i,j,k)
                stop
             end if
          end do
       end do
    end do
    !-------------
    !step 4
    do i=1,ndim
       do j=1,natom
          dofi=(j-1)*ndim + i
          if (.not. use_mkl) then
             call nmtransform_backward(newv(:,i,j), vprop(:,i,j), 0)
             call nmtransform_backward(newx(:,i,j), x(:,i,j),dofi)
          else
             call nmtransform_backward_nr(newv(:,i,j), vprop(:,i,j),0)
             call nmtransform_backward_nr(newx(:,i,j), x(:,i,j),dofi)
          end if
       end do
    end do

    ! do i=1,n
    !    call Vprime(x(i,:,:),force(i,:,:))
    !    do j=1,ndim
    !       do k=1,natom
    !          newv(i,j,k)= vprop(i,j,k) - 0.5d0*force(i,j,k)*dt 
    !          if (newv(i,j,k) .ne. newv(i,j,k)) then
    !             write(*,*) "NaN in pot propagation"
    !             stop
    !          end if
    !       end do
    !    end do
    ! end do

    deallocate(pprop)
    deallocate(newv, newx,force)
    deallocate(q,p)
    return
  end subroutine time_step_test_pile

  !-----------------------------------------------------
  !-----------------------------------------------------
  !normal mode forward transformation for linear polymer
  subroutine nmtransform_forward_nr(xprop, qprop, bead)
    implicit none
    double precision::        xprop(:), qprop(:)
    double precision, allocatable:: qprop_nr(:)
    integer::                 i,j,bead
    allocate(qprop_nr(1:n+1))
    qprop_nr(1)=0.0d0
    qprop_nr(2:n+1)=xprop(1:n)
    call sinft(qprop_nr)
    qprop(1:n)= qprop_nr(2:n+1)*sqrt(2.0d0/(dble(n+1)))
    if (bead.gt.0) qprop(:)= qprop(:) - beadvec(:,bead)
    deallocate(qprop_nr)
    return
  end subroutine nmtransform_forward_nr
  !-----------------------------------------------------
  !-----------------------------------------------------
  !normal mode backward transformation for linear polymer
  subroutine nmtransform_backward_nr(qprop, xprop,bead)
    implicit none
    double precision::        xprop(:), qprop(:)
    double precision, allocatable:: xprop_nr(:)
    integer::                 i,j,bead
    allocate(xprop_nr(1:n+1))
    if (bead .gt. 0) then
       xprop_nr(2:n+1)=qprop(1:n) + beadvec(1:n, bead)
    else 
       xprop_nr(2:n+1)=qprop(1:n)
    end if
    call sinft(xprop_nr)
    xprop(1:n)= xprop_nr(2:n+1)*sqrt(2.0d0/(dble(n+1)))
    deallocate(xprop_nr)
    return
  end subroutine nmtransform_backward_nr

  !-----------------------------------------------------
  !-----------------------------------------------------
  subroutine time_step_higher(xprop, vprop)
    implicit none
    double precision::    xprop(:,:,:), vprop(:,:,:), omegak
    double precision::    pprop(n*ndof)
    double precision, allocatable::  force(:,:,:)
    double precision, allocatable::  q(:,:,:), p(:,:,:)
    integer::             i,j,k,dofi

    allocate(p(n,ndim,natom),q(n,ndim,natom))

    q(:,:,:)= xprop(:,:,:)
    p(:,:,:)= vprop(:,:,:)
    call step_nm(alpha4*dt, q, p, .true.)
    call step_v(beta3*dt,q, p)
    call step_nm(alpha3*dt, q, p, .true.)
    call step_v(beta2*dt,q, p)
    call step_nm(alpha2*dt, q, p, .true.)
    call step_v(beta1*dt,q, p)
    call step_nm(alpha1*dt, q, p,.false.)
    
    errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n*ndof,pprop,0.0d0,1.0d0)
    do i=1,ndim
       do j=1,natom
          dofi= (j-1)*ndim+i
          do k=1,n
             p(k,i,j)= (c1(j,k)**2)*p(k,i,j) + &
                  sqrt(beadmass(j,k)/betan)*c2(j,k)*sqrt(1.0+c1(j,k)**2)*pprop((dofi-1)*n +k)
          end do
          if (.not. use_mkl) then
             call nmtransform_backward(p(:,i,j), vprop(:,i,j), 0)
             call nmtransform_backward(q(:,i,j), xprop(:,i,j), dofi)
          else
             call nmtransform_backward_nr(p(:,i,j), vprop(:,i,j), 0)
             call nmtransform_backward_nr(q(:,i,j), xprop(:,i,j), dofi)
          end if
       end do
    end do

    deallocate(q,p)
    return
  end subroutine time_step_higher

  !-----------------------------------------------------
  !-----------------------------------------------------
  subroutine step_nm(time, x, p, transform)
    implicit none
    integer::          i, j, k, dofi
    double precision:: omegak, newpi(n, ndim, natom)
    double precision:: q(n, ndim, natom), pi(n, ndim, natom)
    double precision:: x(:,:,:), p(:,:,:) , time
    logical::          transform

       do i=1,ndim
          do j=1,natom
             dofi= (j-1)*ndim+i
             if (.not. use_mkl) then
                call nmtransform_forward(p(:,i,j), pi(:,i,j), 0)
                call nmtransform_forward(x(:,i,j), q(:,i,j), dofi)
             else
                call nmtransform_forward_nr(p(:,i,j), pi(:,i,j), 0)
                call nmtransform_forward_nr(x(:,i,j), q(:,i,j), dofi)
             end if
          end do
       end do

    do i=1, n
       do j=1,ndim
          do k=1,natom
             omegak= sqrt(mass(k)/beadmass(k,i))*lam(i)
             newpi(i,j,k)= pi(i,j,k)*cos(time*omegak) - &
                  q(i,j,k)*omegak*beadmass(k,i)*sin(omegak*time)
             q(i,j,k)= q(i,j,k)*cos(time*omegak) + &
                  pi(i,j,k)*sin(omegak*time)/(omegak*beadmass(k,i))
             if (newpi(i,j,k) .ne. newpi(i,j,k)) then
                write(*,*) "NaN in 1st NM propagation"
                stop
             end if
          end do
       end do
    end do
    pi(:,:,:)=newpi(:,:,:)
    
    if (transform) then
    do i=1,ndim
       do j=1,natom
          dofi=(j-1)*ndim + i
          if (.not. use_mkl) then
             call nmtransform_backward(pi(:,i,j), p(:,i,j), 0)
             call nmtransform_backward(q(:,i,j), x(:,i,j),dofi)
          else
             call nmtransform_backward_nr(pi(:,i,j), p(:,i,j),0)
             call nmtransform_backward_nr(q(:,i,j), x(:,i,j),dofi)
          end if
       end do
    end do
    else
       x(:,:,:)= q(:,:,:)
       p(:,:,:)= pi(:,:,:)
    end if

       
  end subroutine step_nm
  !-----------------------------------------------------
  !-----------------------------------------------------

  subroutine step_v(time,x,p)
    implicit none
    integer::         i,j,k
    double precision:: time, p(:,:,:), x(:,:,:)
    double precision:: force(ndim, natom)

    do i=1,n
       call Vprime(x(i,:,:),force(:,:))
       do j=1,ndim
          do k=1,natom
             p(i,j,k)= p(i,j,k) - force(j,k)*time
             if (p(i,j,k) .ne. p(i,j,k)) then
                write(*,*) "NaN in pot propagation"
                stop
             end if
          end do
       end do
    end do

  end subroutine step_v
  !-----------------------------------------------------
  !-----------------------------------------------------

end module verletint