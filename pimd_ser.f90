program pimd
  use mcmod_mass
  use verletint
  use MKL_VSL_TYPE
  use MKL_VSL

  implicit none
  !general variables
  integer::                        i, j, irej, randind,count,k,jj, nrep
  integer::                        thermostat,dofi, npath
  logical::                        potential_test, instapath
  double precision::               Isum, Isum0,Isqprev, xcutoff, UHnew, UHold
  double precision::               diagonals, offdiags, Pacc, delta, stdev
  double precision::               UHinitial, f, randno,answer, com, totm
  double precision::               lndetj, lndetj0, skink, theta, phi
  double precision::               dHdr, dx, finalI,sigmaA, pot, grad, rmspath
  double precision::               a,b,c, xmiddle
  double precision, allocatable::  x(:,:,:), vel(:), tempv(:), tempp(:)
  double precision, allocatable::  p(:,:,:), xi(:), dbdxi(:,:,:), Vpath(:)
  double precision, allocatable::  path(:,:,:), lampath(:), splinepath(:), y(:,:,:)
  integer::                        ndofrb, dummy
  character::                      dummylabel, dummystr(28)
  !gauss-legendre variables
  integer::                        nintegral,ii, time1, time2, imax, irate
  double precision, allocatable::  weights(:), xint(:,:,:), integrand(:), sigma(:)
  !lapack variables
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol
  integer::                        nout, ldz, lwork, liwork, info
  integer,allocatable::            isuppz(:), iwork(:)
  double precision, allocatable::  diag(:), work(:), z(:,:)
  namelist /MCDATA/ n, beta, NMC, Noutput,dt, iprint,imin,tau,npath,&
       nintegral,nrep, use_mkl, thermostat, ndim, natom, xunit,&
       potential_test, instapath, gamma
  !-------------------------
  !Set default system parameters then read in namelist
  iprint=.false.
  n= 100
  beta= 100.0d0
  NMC= (5e6)
  Noutput= 1e5
  dt= 1d-3
  nintegral= 5
  nrep=1
  thermostat=1
  ndim=3
  natom=1
  xunit=1
  use_mkl=.false.
  imin=0
  tau=1.0d0
  npath=0
  potential_test=.false.
  instapath=.false.
  gamma=1.0d0

  read(5, nml=MCDATA)
  betan= beta/dble(n+1)
  call system_clock(time1, irate, imax)
  write(*,*)"Running with parameters (in a.u.):"
  write(*,*) "beta, betan, n=", beta, betan, n
  write(*,*) "NMC, Noutput, dt=", NMC, Noutput, dt
  write(*,*) "with integration points", nintegral
  write(*,*) "tau=", tau
  if (thermostat .eq. 1) then
     write(*,*) "Running with Andersen thermostat"
  else if (thermostat .eq. 2) then
     write(*,*) "Running with Langevin thermostat"
  end if
  call V_init()
  !-------------------------
  !-------------------------
  !Read in initial wells, and masses
  allocate(mass(natom),label(natom))
  open(18, file="masses.dat", status="old")
  do j=1,natom
     read(18,*)label(j), mass(j)
  end do
  close(18)

  ndof=ndim*natom

  if (.not. instapath) then
     allocate(path(npath, ndim, natom), lampath(npath), Vpath(npath))
     open(15, file="path.xyz")
     do i=1, npath
        read(15,*) dummy
        read(15,'(28A)') dummystr!, dummyE !'(A,G25.15)'
        do j=1, natom
           read(15,*) dummylabel, path(i, 1, j),path(i, 2, j),path(i, 3, j) !'(A2,4X,3G20.10)'
        end do
        if (i.eq.1) then
           lampath(1)=0.0d0
        else
           lampath(i)= lampath(i-1) + eucliddist(path(i-1,:,:), path(i,:,:))!dble(i-1)/dble(npath-1)
        end if
        Vpath(i)= V(path(i,:,:))
     end do
     lampath(:)= lampath(:)/lampath(npath)
     close(15)

     !xunit=1 means bohr
     !xunit=2 means angstroms
     if (xunit .eq. 2) then
        path(:,:,:) = path(:,:,:)/0.529177d0
     end if
  else
     npath=n
     allocate(well1(ndim,natom), well2(ndim,natom))
     open(15, file="well1.dat", status="old")
     open(16, file="well2.dat", status="old")
     do j=1,natom
        read(15,*) well1(1,j),well1(2,j),well1(3,j)
        read(16,*) well2(1,j),well2(2,j),well2(3,j)
     end do
     close(15)
     close(16)
     if (xunit .eq. 2) then
        well1(:,:)= well1(:,:)/0.529177d0
        well2(:,:)= well2(:,:)/0.529177d0
     end if
     V0=V(well1)
     write(*,*) "Potential at wells:", V(well1), V(well2)
     !xunit=1 means bohr
     !xunit=2 means angstroms
     allocate(path(npath, ndim, natom), lampath(npath), Vpath(npath))
     call instanton(path,well1,well2)
     write(*,*) "Found instanton."
     open(19, file="instanton.xyz")
     do i=1,npath
        write(19,*) "9"
        write(19,*) "Energy of minimum",i
        do j=1, natom
           write(19,*)  label(j), (path(i,k,j), k=1,ndim)
        end do
        if (i.eq.1) then
           lampath(1)=0.0d0
        else
           lampath(i)= lampath(i-1) + eucliddist(path(i-1,:,:), path(i,:,:))!dble(i-1)/dble(npath-1)
        end if
        Vpath(i)= V(path(i,:,:))
     end do
     lampath(:)= lampath(:)/lampath(npath)
     close(19)
  end if
  allocate(xtilde(n, ndim,natom))
  xtilde(:,:,:)= path(:,:,:)
  !-------------------------
  !Find the centre to make sure this is symmetric
  ! allocate(splinepath(npath))
  ! write(*,*) size(lampath), size(Vpath), size(splinepath)
  ! call spline(lampath, Vpath, 1.0d31, 1.0d31, splinepath)
  ! xmiddle= findmiddle(0.4d0, 0.6d0, lampath, Vpath, splinepath)
  ! deallocate(splinepath)
  ! a= 2.0d0  - 4.0d0*xmiddle
  ! b= 4.0d0*xmiddle - 1.0d0
  ! c= 0.0d0
  ! do i=1, npath
  !    if (a.ge.0) then
  !       lampath(i)= -0.5d0*b/a + sqrt((lampath(i)/a) + (0.5*b/a)**2)
  !    else
  !       lampath(i)= -0.5d0*b/a - sqrt((lampath(i)/a) + (0.5*b/a)**2)
  !    end if
  !    write(*,*) i, dble(i-1)/dble(npath-1), lampath(i)
  ! end do
  ! write(*,*)"Centred,", a, b, xmiddle

  !-------------------------
  !-------------------------
  !set up gauss-legendre quadrature arrays
  allocate(weights(nintegral), xi(nintegral),xint(nintegral,ndim,natom), integrand(nintegral))
  allocate(dbdxi(nintegral,ndim,natom), sigma(nintegral),splinepath(npath))
  if (potential_test) allocate(y(1000,ndim,natom))
  call gauleg(0d0,1.0d0,xi, weights,nintegral)
  ! do i=1,ndim
  !    do j=1,natom
  !       do k=1,1000
  !          y(k,i,j)=(well2(i,j)-well1(i,j))*dble(k-1)/dble(999) + well1(i,j)
  !       end do
  !       do k=1, nintegral
  !          xint(k,i,j)= (well2(i,j)-well1(i,j))*xi(k) + well1(i,j)
  !          dbdxi(k,i,j)= (well2(i,j)-well1(i,j))
  !       end do
  !    end do
  ! end do
  do i=1,ndim
     do j=1,natom
        splinepath(:)=0.0d0
        call spline(lampath(:), path(:,i,j), 1.0d31, 1.0d31, splinepath(:))
        do k=1, nintegral
           xint(k,i,j)= splint(lampath, path(:,i,j), splinepath(:), xi(k))
           dbdxi(k,i,j)= splin_grad(lampath, path(:,i,j), splinepath(:), xi(k))
        end do
        if (potential_test) then
           do k=1,1000
              y(k,i,j)= splint(lampath, path(:,i,j), splinepath(:), dble(k)*1d-3)
           end do
        end if
     end do
  end do

  do i=1, nintegral
     write(*,*) xi(i), V(xint(i,:,:))
  end do

  if (potential_test) then
     open(20, file="pottest.dat")
     open(30, file="interp_path.xyz")
     open(40, file="path_ref.xyz")
     open(50, file="pot_ref.dat")
     do i=1,1000
        pot= V(y(i,:,:))
        write(20,*) dble(i)*1d-3, pot!, grad
        write(30,*) "9"
        write(30,*) "Energy of minimum",i
        do j=1, natom
           write(30,*)  label(j), (y(i,k,j), k=1,ndim)
        end do
     end do
     ! do i=1,npath
     !    pot= V(path(i,:,:))
     !    write(50,*) lampath(i),pot
     !    write(40,*) "9"
     !    write(40,*) lampath(i)
     !    do j=1, natom
     !       write(40,*)  label(j), (path(i,k,j), k=1,ndim)
     !    end do
     ! end do
     close(20)
     close(30)
     close(40)
     deallocate(y)
  end if
  deallocate(splinepath,lampath)

  !-------------------------
  !-------------------------
  !Start loop over integration points
  integrand(:)=0.0d0
  call random_seed()
  allocate(x(n,ndim,natom), vel(n),p(n,ndim,natom))
  allocate(tempp(n), tempv(n))
  call alloc_nm()

  do ii=1, nintegral
     sigma(ii)=0.0d0
     integrand(ii)=0.0d0
     do jj=1, nrep
        call init_nm(path(1,:,:),xint(ii,:,:))
        call init_path(path(1,:,:), xint(ii,:,:), x, p)
        if (thermostat .eq. 1) then
           call propagate_pimd_nm(x,p, path(1,:,:),xint(ii,:,:),dbdxi(ii,:,:),dHdr)
        else if (thermostat .eq. 2) then
           call propagate_pimd_pile(x,p, path(1,:,:),xint(ii,:,:), dbdxi(ii,:,:),dHdr)
        else if (thermostat .eq. 3) then
           call propagate_pimd_higher(x,p, path(1,:,:),xint(ii,:,:),dbdxi(ii,:,:),dHdr)
        else
           write(*,*) "Incorrect thermostat option."
           stop
        end if
        integrand(ii)=integrand(ii) + dHdr/(dble(nrep)*betan**2)
        sigma(ii)= sigma(ii)+ (dHdr/betan**2)**2
     end do
     sigma(ii)= sigma(ii)/dble(nrep)
     sigma(ii)= sigma(ii) - integrand(ii)**2
     ! write(*,*)"-----------------"
     write(*,*) ii, xi(ii), integrand(ii), sqrt(sigma(ii))
     ! write(*,*)"-----------------"
  end do
  call free_nm()
  deallocate(x, vel,p, tempv, tempp)
  answer= 0.0d0
  sigmaA=0.0d0

  do i=1, nintegral
     if (integrand(i) .eq. integrand(i)) &
          answer= answer+ weights(i)*integrand(i)
     sigmaA= sigmaA + sigma(i)*weights(i)**2
  end do

  finalI= exp(-answer*betan)
  write(*,*) "beta, betan, n=", beta, betan, n
  write(*,*) "final Delta A=", answer , "+/-", sqrt(sigmaA)/sqrt(dble(nrep))
  write(*,*) "q/q0=", finalI, betan*sqrt(sigmaA)*finalI/sqrt(dble(nrep))
  write(*,*) "q0/q=", 1.0d0/finalI, betan*sqrt(sigmaA)/(finalI*sqrt(dble(nrep)))
  call system_clock(time2, irate, imax)
  write(*,*) "Ran in", dble(time2-time1)/dble(irate), "s"
  deallocate(path, xtilde)
end program pimd
