program pimd
  use mcmod_mass
  use verletint
  use MKL_VSL_TYPE
  use MKL_VSL
  use instantonmod

  implicit none
  !general variables
  integer::                        i, j, irej, randind,count,k,jj, nrep
  integer::                        thermostat,dofi, npath
  logical::                        potential_test, centre, instapath, readpath, alignwell
  double precision::               Isum, Isum0,Isqprev, xcutoff, UHnew, UHold
  double precision::               diagonals, offdiags, Pacc, delta, stdev
  double precision::               UHinitial, f, randno,answer, com, totm
  double precision::               lndetj, lndetj0, skink, theta, phi
  double precision::               dHdr, dx, finalI,sigmaA, pot, grad, rmspath
  double precision::               a,b,c, xmiddle
  double precision::               theta1, theta2, theta3
  double precision, allocatable::  origin(:), wellinit(:,:), vgrad(:,:)
  double precision, allocatable::  x(:,:,:), vel(:), tempv(:), tempp(:), initpath(:,:)
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
       potential_test, instapath, gamma,centre, alignwell, fixedends
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
  centre=.false.
  gamma=1.0d0
  dHdrlimit=-1.0
  alignwell=.false.
  readpath=.true.
  fixedends=.true.

  read(5, nml=MCDATA)
  betan= beta/dble(n+1)
  call system_clock(time1, irate, imax)
  write(*,*)"Running with parameters (in a.u.):"
  write(*,*) "beta, betan, n=", beta, betan, n
  write(*,*) "NMC, Noutput, dt=", NMC, Noutput, dt
  write(*,*) "with integration points", nintegral
  write(*,*) "tau=", tau
  write(*,*) "gamma=", gamma
  if (thermostat .eq. 1) then
     write(*,*) "Running with Andersen thermostat"
  else if (thermostat .eq. 2) then
     write(*,*) "Running with Langevin thermostat"
  end if
  call V_init(0)
  !-------------------------
  !-------------------------
  !Read in initial wells, and masses
  allocate(origin(ndim))
  allocate(mass(natom),label(natom),xtilde(n, ndim, natom))
  open(18, file="masses.dat", status="old")
  do j=1,natom
     read(18,*)label(j), mass(j)
  end do
  close(18)

  ndof=ndim*natom
  totdof=n*ndim*natom
  allocate(well1(ndim,natom), well2(ndim,natom), wellinit(ndim,natom))
  open(15, file="well1.dat", status="old")
  open(16, file="well2.dat", status="old")
  do j=1,natom
     read(15,*) (well1(i,j),i=1,ndim)
     read(16,*) (well2(i,j),i=1,ndim)
  end do
  close(15)
  close(16)
  ! xunit=1 means bohr
  ! xunit=2 means angstroms
  if (xunit .eq. 2) then
     well1(:,:)= well1(:,:)/0.529177d0
     well2(:,:)= well2(:,:)/0.529177d0
  end if
  call get_align(well1,theta1, theta2, theta3, origin)
  wellinit(:,:)= well1(:,:)
  call align_atoms(wellinit, theta1, theta2, theta3, origin, well1)
  if (alignwell) call get_align(well2,theta1, theta2, theta3, origin)
  wellinit(:,:)= well2(:,:)
  call align_atoms(wellinit, theta1, theta2, theta3, origin, well2)
  V0=V(well1)
  write(*,*) "Potential at wells:", V(well1), V(well2)
  allocate(vgrad(ndim,natom))
  call Vprime(well1, vgrad)
  write(*,*) "With norm of grad:", norm2(reshape(vgrad, (/ndim*natom/)))

  !-------------------------
  !-------------------------
  !Read in initial wells, and masses
  if (readpath) then
        call read_path(instapath,centre,npath,path,Vpath,splinepath,lampath)
  end if
  !-------------------------
  !-------------------------
  !set up gauss-legendre quadrature arrays
  allocate(weights(nintegral), xi(nintegral),xint(nintegral,ndim,natom))
  allocate(dbdxi(nintegral,ndim,natom))
  call gauleg(0d0,1.0d0,xi, weights,nintegral)

  do i=1,ndim
     do j=1,natom
        do k=1, nintegral
           xint(k,i,j)= splint(lampath, path(:,i,j), splinepath(:), xi(k))
           dbdxi(k,i,j)= splin_grad(lampath, path(:,i,j), splinepath(:), xi(k))
        end do
     end do
  end do
  deallocate(splinepath)

  do i=1, nintegral
     write(*,*) xi(i), V(xint(i,:,:)), norm2(reshape(dbdxi(i,:,:), (/ndim*natom/)))
  end do

  !-------------------------
  !-------------------------
  !Start loop over integration points
  allocate(integrand(nintegral), sigma(nintegral))
  integrand(:)=0.0d0
  call random_seed()
  allocate(x(n,ndim,natom), vel(n),p(n,ndim,natom))
  allocate(tempp(n), tempv(n))
  allocate(transmatrix(n,n),beadmass(natom,n),beadvec(n,ndof), lam(n))
  call alloc_nm(1)
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
  deallocate(path, xtilde,lampath)
end program pimd
