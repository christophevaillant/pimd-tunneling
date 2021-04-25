include 'mkl_vsl.fi'
module verletint
  use mcmod_mass
  use utilsmod
  use MKL_VSL_TYPE
  use MKL_VSL
  use nr_fft, only : sinft
  use instantonmod
  implicit none

  integer::                         NMC, Noutput, imin, restart, iproc, restartnmc,restartunit
  integer, allocatable::            ipar(:)
  double precision::                dt, dHdrlimit
  double precision,allocatable::    transmatrix(:,:),dpar(:), beadvec(:,:)
  double precision,allocatable::    beadmass(:,:), lam(:)
  double precision, allocatable::  path(:,:,:), lampath(:), splinepath(:,:,:), Vpath(:)
  double precision, allocatable::  hesspath(:,:,:), splinehess(:,:,:)  
  double precision, allocatable::   c1(:,:), c2(:,:)
  double precision::                alpha1, alpha2, alpha3, alpha4
  double precision::                beta1, beta2, beta3, tau, gamma
  logical::                         iprint, use_mkl, cayley, readhess
  integer::                         errcode_poisson, rmethod_poisson
  integer::                         brng_poisson, seed_poisson, stat
  type (vsl_stream_state)::         stream_poisson,stream_normal
  integer::                         errcode_normal, rmethod_normal, brng_normal
  integer::                         seed_normal, rk
  character::                       procstr

contains
  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(xi, x, p)
    double precision::  x(:,:,:), p(:,:,:), xi
    double precision, allocatable:: vel(:), tempp(:), dists(:), interphess(:,:,:)
    double precision, allocatable:: etasquared(:),eigvecs(:,:),tempx(:)
    double precision::  stdev, xieff
    integer::           i,j,k, dofi, imin(1),i2,j2,k2,idof1,idof2

    allocate(vel(n),tempp(n), dists(n))
    do i=1,ndim
       do j=1,natom
          do k=1,n
             xieff= dble(k-1)*xi/dble(n-1)

             x(k,i,j)= splint(lampath, path(:,i,j), splinepath(:,i,j), xieff)
          end do
       end do
    end do
    if (readhess) then
       allocate(interphess(n,ndof,ndof), etasquared(totdof), eigvecs(totdof,totdof))
       allocate(tempx(totdof))
       !------------------
       !generate random numbers
       stdev=sqrt(1.0d0/beta)
       errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,totdof,tempx,0.0d0,stdev)
       !------------------
       !interpolate hessian
       do i=1,ndof
          do j=1,ndof
             do k=1,n
                xieff= dble(k-1)*xi/dble(n-1)
                interphess(k,i,j)= splint(lampath, hesspath(:,i,j), splinehess(:,i,j), xieff)
             end do
          end do
       end do

       !------------------
       !calculate "instanton" e.vectors and e.values
       call detJ(x,etasquared,.false.,interphess,eigvecs)
       ! gamma=tau*sqrt(abs(etasquared(totdof)))
       ! write(*,*) iproc, i, 2.0d0*pi/sqrt(abs(etasquared(totdof))), gamma
       !------------------
       !add random displacements to path
       do idof1=2,totdof
          if (etasquared(idof1) .lt. 0.0) cycle
          do i2= 1,n
             do j2=1,ndim
                do k2=1,natom
                   idof2= natom*(j2-1 + ndim*(i2-1)) +k2
                   x(i2,j2,k2)= x(i2,j2,k2) + sqrt(1.0/(etasquared(idof1)*mass(k2)))&
                        *tempx(idof1)*eigvecs(idof2,idof1)
                end do
             end do
          end do
       end do
             
       deallocate(interphess, etasquared,eigvecs,tempx)
    end if

       ! write(procstr,"(I1)") iproc
       ! open(400,file="initial"//trim(procstr)//".xyz")
       ! do i=1,n
       !    write(400,*) natom
       !    write(400,*) i
       !    do k2=1,natom
       !       write(400,*) label(k2), x(i,1,k2)*0.529177d0, x(i,2,k2)*0.529177d0, x(i,3,k2)*0.529177d0
       !    end do
       ! end do
       ! close(400)
    ! x(:,:,:)= xtilde(:,:,:)

    do i=1,ndim
       do k=1,natom
          dofi=(k-1)*ndim +i
          stdev=sqrt(1.0d0/betan)
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,vel,0.0d0,stdev)
          do j=1,n
             vel(j)= vel(j)*sqrt(beadmass(k,j))
          end do
          call nmtransform_backward(vel, tempp, 0)
          do j=1,n
             p(j,i,k)= tempp(j)
          end do
       end do
    end do

    deallocate(vel, tempp)
    return
  end subroutine init_path

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for calculating gauss-legendre points
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

  subroutine write_restart(xprop,pprop,ii,dHdr)
    implicit none
    double precision, intent(in)::   xprop(:,:,:), pprop(:,:,:), dHdr
    integer, intent(in)::            ii
    integer::                        i,j,k

    rewind(restartunit)
    do i=1,n
       write(restartunit,*) natom
       write(restartunit,*) dHdr
       do j=1,natom
          write(restartunit,*) label(j), xprop(i,1,j), xprop(i,2,j), xprop(i,3,j)
       end do
    end do
    do i=1,n
       write(restartunit,*) natom
       write(restartunit,*) ii
       do j=1,natom
          write(restartunit,*) label(j), pprop(i,1,j), pprop(i,2,j), pprop(i,3,j)
       end do
    end do

    return
  end subroutine write_restart

  !-----------------------------------------------------
  !-----------------------------------------------------
  !main function for propagating a MD trajectory
  subroutine propagate_pimd_nm(xprop,vprop,a,b,dbdl,dHdr)
    implicit none
    double precision::     xprop(:,:,:), vprop(:,:,:),randno, sigma
    double precision::     a(:,:),b(:,:), dHdr, totenergy, kin,dbdl(:,:), contr
    double precision, allocatable:: pprop(:), tempp(:), tempv(:)
    integer::              i,j,k,l,count, time1,time2,imax, irate, dofi,ii
    integer (kind=4)::     rkick(1)

    allocate(pprop(n), tempp(n),tempv(n))
    count=0
    if (restart .lt. 2) dHdr=0.0d0
    kin=0.0d0
    errcode_poisson = virngpoisson( rmethod_poisson, stream_poisson, 1, rkick(1), dble(Noutput))
    do ii=1, NMC, 1
       count=count+1
       if ((mod(ii,Noutput) .eq. 0) .and. (restart .gt. 0)) then
          call write_restart(xprop,vprop,ii+restartnmc,dHdr)
       end if
       if (count .ge. rkick(1)) then
          count=0
          if (iprint) write(*,*) 100*dble(ii)/dble(NMC),dHdr/(betan**2*dble(ii-imin))!dHdr/dble(i)
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
       call time_step_nm(xprop, vprop)
       if (ii.gt.imin) then
          contr=0.0d0
          do j=1,ndim
             do k=1,natom
                contr=contr+mass(k)*(-xprop(n,j,k))*dbdl(j,k)
             end do
          end do
          dHdr= dHdr+contr
       end if
    end do
    if ((restart .gt. 0)) call write_restart(xprop,vprop,nmc+restartnmc,dHdr)
    dHdr=dHdr/dble(NMC+restartnmc-imin)
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
  subroutine time_step_nm(xprop, pprop)
    implicit none
    double precision, intent(inout)::    xprop(:,:,:), pprop(:,:,:)
    double precision, allocatable::  force(:,:,:), pip(:,:,:)

    allocate(force(n, ndim, natom), pip(n, ndim, natom))
    call step_nm(0.5d0*dt,xprop,pprop ,.true.)
    call step_v(dt, xprop, pprop, force, .true.)
    call step_nm(0.5d0*dt,xprop,pprop ,.true.)
    deallocate(force, pip)
    return
  end subroutine time_step_nm
  !-----------------------------------------------------
  !-----------------------------------------------------
  !initialize normal mode routines
  subroutine init_nm(a,b)
    double precision::    a(:,:),b(:,:)
    integer::             i,j,k,l, dofi,tdofi

    do i=1, n
       do k=1,ndim
          do j=1, natom
             dofi= natom*(k-1 + ndim*(i-1)) +j
             lam(i)= 2.0d0*sin(dble(i)*pi/dble(2*n +2))/betan
             beadmass(j,i)=mass(j)*(lam(i)*tau)**2
          end do
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
             tdofi= natom*(j-1 + ndim*(i-1)) +k
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
    integer:: iproc

    call system_clock(itime,irate,imax)
    seed_normal= mod(itime+5*iproc,1000)
    call system_clock(itime,irate,imax)
    seed_poisson= mod(itime+5*iproc,1000)
    brng_normal = VSL_BRNG_MT19937
    brng_poisson = VSL_BRNG_MT19937
    rmethod_normal = VSL_RNG_METHOD_GAUSSIAN_ICDF
    rmethod_poisson = VSL_RNG_METHOD_POISSON_POISNORM
    errcode_normal = vslnewstream( stream_normal, brng_normal, seed_normal )
    errcode_poisson = vslnewstream( stream_poisson, brng_poisson, seed_poisson )
    
    write(*,*)"proc", iproc, "running with seeds:", seed_normal, seed_poisson

  end subroutine alloc_nm
  !-----------------------------------------------------
  !-----------------------------------------------------
  !main function for propagating a MD trajectory with Langevin thermostat
  subroutine propagate_pimd_pile(xprop,vprop,a,b,dbdl,xi,dHdr)
    implicit none
    double precision::     xprop(:,:,:), vprop(:,:,:),randno, Eold
    double precision::     potenergy, springenergy, kinenergy, totenergy,xi
    double precision::     a(:,:),b(:,:), dHdr, dbdl(:,:), contr
    integer::              i,j,k,count, time1,time2,imax, irate, jj,idof,ii
    integer (kind=4)::     rkick(1)

    allocate(c1(natom,n), c2(natom,n))
    do i=1,n
       do j=1,natom
          c1(j,i)= exp(-gamma*dt*lam(i)*sqrt(mass(j)/beadmass(j,i)))!
          c2(j,i)= sqrt(1.0d0- c1(j,i)**2)
       end do
    end do
    count=0
    if (restart .lt. 2)dHdr=0.0d0
    do ii=1, NMC, 1
       count=count+1
       if (count .ge. Noutput .and. ii .gt. imin) then
          count=0
          if (iprint) write(*,*) 100*dble(ii)/dble(NMC-imin), dHdr/(betan**2*dble(ii-imin))
          if (restart .gt. 0) call write_restart(xprop,vprop,ii+restartnmc,dHdr)
       end if
       call time_step_pile(xprop, vprop)
       if (ii.gt.imin) then
          contr=0.0d0
          do j=1,ndim
             do k=1,natom
                contr=contr+mass(k)*(-xprop(n,j,k))*dbdl(j,k)
             end do
          end do
          if (abs(contr) .lt. dHdrlimit .or. dHdrlimit .lt. 0.0) then
             dHdr= dHdr+contr
          else
             write(*,*) "Over limit", contr, ", reinitialize path"
             call init_path(xi, xprop, vprop)
          end if
       end if
    end do
    if ((restart .gt. 0)) call write_restart(xprop,vprop,nmc+restartnmc,dHdr)
    dHdr= dHdr/dble(NMC+restartnmc-imin)
    deallocate(c1, c2)
    return
  end subroutine propagate_pimd_pile


  !-----------------------------------------------------
  !-----------------------------------------------------
  !routine taking a step in time using normal mode verlet algorithm
  !from ceriotti et al 2010 paper.
  subroutine time_step_pile(xprop, pprop)
    implicit none
    double precision, intent(inout)::    xprop(:,:,:), pprop(:,:,:)
    double precision, allocatable::  force(:,:,:), pip(:,:,:)

    allocate(force(n, ndim, natom), pip(n, ndim, natom))
    call step_v(dt, xprop, pprop, force, .true.)
    call step_nm(0.5d0*dt,xprop,pprop ,.true.)
    call step_langevin(pprop)
    call step_nm(0.5d0*dt,xprop,pprop ,.true.)

    return
  end subroutine time_step_pile
  !-----------------------------------------------------
  !-----------------------------------------------------
  !routine taking a step in time using normal mode verlet algorithm
  !from ceriotti et al 2010 paper.
  subroutine time_step_ffpile(xprop, pprop)
    implicit none
    double precision, intent(inout)::    xprop(:,:,:), pprop(:,:,:)
    double precision, allocatable::  force(:,:,:), pip(:,:,:)

    allocate(force(n, ndim, natom), pip(n, ndim, natom))
    call step_v(dt, xprop, pprop, force, .true.)
    call step_nm(0.5d0*dt,xprop,pprop ,.true.)
    call step_fflangevin(pprop)
    call step_nm(0.5d0*dt,xprop,pprop ,.true.)

    deallocate(force, pip)
    return
  end subroutine time_step_ffpile

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
  subroutine step_nm(time, x, p, transform)
    implicit none
    integer::          i, j, k, dofi
    double precision:: omegak, newpi(n, ndim, natom)
    double precision:: q(n, ndim, natom), pip(n, ndim, natom)
    double precision:: x(:,:,:), p(:,:,:) , time
    logical::          transform

       do i=1,ndim
          do j=1,natom
             dofi= (j-1)*ndim+i
             if (.not. use_mkl) then
                call nmtransform_forward(p(:,i,j), pip(:,i,j), 0)
                call nmtransform_forward(x(:,i,j), q(:,i,j), dofi)
             else
                call nmtransform_forward_nr(p(:,i,j), pip(:,i,j), 0)
                call nmtransform_forward_nr(x(:,i,j), q(:,i,j), dofi)
             end if
          end do
       end do

    do i=1, n
       do j=1,ndim
          do k=1,natom
             dofi= natom*(j-1 + ndim*(i-1)) +k
             omegak= sqrt(mass(k)/beadmass(k,i))*lam(i)
             if (cayley) then
                newpi(i,j,k)= pip(i,j,k)*(4.0d0- omegak**2*time**2) &
                     - 4.0d0*q(i,j,k)*beadmass(k,i)*omegak**2*time
                newpi(i,j,k)=newpi(i,j,k)/(4.0d0+ omegak**2*time**2)
                q(i,j,k)= q(i,j,k)*(4.0d0- omegak**2*time**2) &
                     + 4.0d0*pip(i,j,k)*time/beadmass(k,i)
                q(i,j,k)=q(i,j,k)/(4.0d0+ omegak**2*time**2)
             else
                newpi(i,j,k)= pip(i,j,k)*cos(time*omegak) - &
                     q(i,j,k)*omegak*beadmass(k,i)*sin(omegak*time)
                q(i,j,k)= q(i,j,k)*cos(time*omegak) + &
                     pip(i,j,k)*sin(omegak*time)/(omegak*beadmass(k,i))
             end if
             if (newpi(i,j,k) .ne. newpi(i,j,k)) then
                write(*,*) "NaN in 1st NM propagation"
                stop
             end if
          end do
       end do
    end do
    pip(:,:,:)=newpi(:,:,:)
    
    if (transform) then
    do i=1,ndim
       do j=1,natom
          dofi=(j-1)*ndim + i
          if (.not. use_mkl) then
             call nmtransform_backward(pip(:,i,j), p(:,i,j), 0)
             call nmtransform_backward(q(:,i,j), x(:,i,j),dofi)
          else
             call nmtransform_backward_nr(pip(:,i,j), p(:,i,j),0)
             call nmtransform_backward_nr(q(:,i,j), x(:,i,j),dofi)
          end if
       end do
    end do
    else
       x(:,:,:)= q(:,:,:)
       p(:,:,:)= pip(:,:,:)
    end if

       
  end subroutine step_nm
  !-----------------------------------------------------
  !-----------------------------------------------------

  subroutine step_v(time,x,p, force, recalculate)
    implicit none
    double precision, intent(in)::  time
    double precision, intent(inout):: force(:, :,:),p(:,:,:), x(:,:,:)
    integer::         i,j,k
    logical, intent(in):: recalculate

    do i=1,n
       if (recalculate) call Vprime(x(i,:,:),force(i,:,:))
       do j=1,ndim
          do k=1,natom
             p(i,j,k)= p(i,j,k) - force(i,j,k)*time
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

  !-----------------------------------------------------
  !-----------------------------------------------------
  subroutine step_fflangevin(pprop)
    implicit none
    double precision, intent(inout):: pprop(:,:,:)
    double precision, allocatable:: p(:,:,:), pk(:)
    integer:: i,j,k, dofi

    allocate(p(n,ndim,natom), pk(n*ndof))
    errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n*ndof,pk,0.0d0,1.0d0)
    do i=1,ndim
       do j=1,natom
          dofi= calcidof(i,j)
          if (.not. use_mkl) then
             call nmtransform_forward(pprop(:,i,j), p(:,i,j), 0)
          else
             call nmtransform_forward_nr(pprop(:,i,j), p(:,i,j), 0)
          end if
          do k=1,n
             pprop(k,i,j)= p(k,i,j) !seems silly, but it isn't!
             p(k,i,j)= (c1(j,k)**2)*p(k,i,j) + &
                  sqrt(beadmass(j,k)/betan)*c2(j,k)*sqrt(1.0+c1(j,k)**2)*pk((dofi-1)*n +k)
          end do
       end do
    end do
    do k=1,n
       do j=1,natom
          p(k,:,j)= norm2(p(k,:,j))*pprop(k,:,j)/norm2(pprop(k,:,j)) !see, not so silly!
       end do
    end do
    do i=1,ndim
       do j=1,natom
          dofi= calcidof(i,j)
          if (.not. use_mkl) then
             call nmtransform_backward(p(:,i,j),pprop(:,i,j), 0)
          else
             call nmtransform_backward_nr(p(:,i,j),pprop(:,i,j), 0)
          end if
       end do
    end do
    deallocate(p,pk)
  end subroutine step_fflangevin

    !-----------------------------------------------------
  !-----------------------------------------------------
  subroutine step_langevin(pprop)
    implicit none
    double precision, intent(inout):: pprop(:,:,:)
    double precision, allocatable:: p(:,:,:), pk(:)
    integer:: i,j,k, dofi

    allocate(p(n,ndim,natom), pk(n*ndof))
    errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n*ndof,pk,0.0d0,1.0d0)
    do i=1,ndim
       do j=1,natom
          dofi= calcidof(i,j)
          if (.not. use_mkl) then
             call nmtransform_forward(pprop(:,i,j), p(:,i,j), 0)
          else
             call nmtransform_forward_nr(pprop(:,i,j), p(:,i,j), 0)
          end if
          do k=1,n
             p(k,i,j)= (c1(j,k)**2)*p(k,i,j) + &
                  sqrt(beadmass(j,k)/betan)*c2(j,k)*sqrt(1.0+c1(j,k)**2)*pk((dofi-1)*n +k)
          end do
       end do
    end do
    do i=1,ndim
       do j=1,natom
          dofi= calcidof(i,j)
          if (.not. use_mkl) then
             call nmtransform_backward(p(:,i,j),pprop(:,i,j), 0)
          else
             call nmtransform_backward_nr(p(:,i,j),pprop(:,i,j), 0)
          end if
       end do
    end do
    deallocate(p,pk)
  end subroutine step_langevin

    !TODO: replace idof calculations with this function
    function calcidof(i,j) !i for dimensions, j for atom
      implicit none
      integer:: calcidof
      integer, intent(in):: i,j
      calcidof= (j-1)*ndim + i
    end function calcidof

end module verletint
