program pimd
  use mcmod_mass
  use verletint
  use MKL_VSL_TYPE
  use MKL_VSL
  implicit none

  include 'mpif.h'
  ! include 'mkl_service.fi'
  !general variables
  integer::                        i, j, irej, randind,count,k,l,jj, nrep
  integer::                        thermostat,dofi, npath
  double precision::               Isum, Isum0,Isqprev, xcutoff, UHnew, UHold
  double precision::               diagonals, offdiags, Pacc, delta, stdev
  double precision::               UHinitial, f, randno,answer, com, totm
  double precision::               lndetj, lndetj0, skink, theta, phi, dummyE
  double precision::               dHdr, dx, finalI,sigmaA, a, b, xmiddle
  double precision::               theta1, theta2, theta3
  double precision, allocatable::  x(:,:,:), initpath(:,:)
  double precision, allocatable::  origin(:), wellinit(:,:)
  double precision, allocatable::  xi(:), dbdxi(:,:,:), y(:,:,:), pinit(:,:,:)
  double precision, allocatable::  path(:,:,:), lampath(:), splinepath(:), Vpath(:)
  integer::                        ndofrb, dummy
  logical::                        instapath, centre
  character::                      dummylabel, dummystr(28)
  !gauss-legendre variables
  integer::                        nintegral,ii, time1, time2, imax, irate
  double precision, allocatable::  weights(:), xint(:,:,:), integrand(:), sigma(:)
  !lapack variables
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol
  integer::                        nout, ldz, lwork, liwork, info
  integer,allocatable::            isuppz(:), iwork(:)
  !MPI variables
  integer::                        ierr, nproc, ncalcs
  integer, allocatable::           mpi_int_send(:)
  double precision, allocatable::  mpi_double_send(:), endpoints(:,:,:), allendpoints(:,:,:)
  double precision, allocatable::  gradpoints(:,:,:), allgradpoints(:,:,:), startpoint(:,:)
  double precision, allocatable::  finalintegrand(:), allintegrands(:)
  integer (kind=4)::               ierrmkl
  integer, dimension(MPI_STATUS_SIZE) :: rstatus
  namelist /MCDATA/ n, beta, NMC, Noutput,dt, iprint,imin,tau,npath, gamma,&
       nintegral,nrep, use_mkl, thermostat, ndim, natom, xunit, instapath, centre,&
       dHdrlimit

  !initialize MPI
  nproc=0
  iproc=0
  ncalcs=0
  ierr=0
  call MPI_Init(ierr)
  !get the processor ID
  call MPI_Comm_rank ( MPI_COMM_WORLD, iproc, ierr )
  !get the total number of processors
  call MPI_Comm_size ( MPI_COMM_WORLD, nproc, ierr )

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
  gamma=1.0d0
  instapath=.false.
  dHdrlimit=-1.0

  !Read in namelist variables for root proc, and spread
  !it to other procs
  allocate(mpi_int_send(8), mpi_double_send(5))
  
  call V_init()
  if (iproc .eq. 0) then
     read(*, nml=MCDATA)
     betan= beta/dble(n+1)
     call system_clock(time1, irate, imax)
     write(*,*)"Running with parameters (in a.u.):"
     write(*,*) "beta, betan, n=", beta, betan, n
     write(*,*) "NMC, Noutput, dt=", NMC, Noutput, dt
     write(*,*) "with integration points", nintegral
     write(*,*) "and", nrep, "repetitions"
     write(*,*) "tau=", tau
     write(*,*) "gamma=", gamma
     if (thermostat .eq. 1) then
        write(*,*) "Running with Andersen thermostat"
     else if (thermostat .eq. 2) then
        write(*,*) "Running with Langevin thermostat"
     end if

     ncalcs= nintegral*nrep/nproc
     if (mod(nintegral*nrep, nproc) .ne. 0) ncalcs=ncalcs+1

  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  ierr=0
  ! call MPI_Bcast(mpi_int_send, 8, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ncalcs, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(N, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(NMC, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Noutput, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(thermostat, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(imin, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ndim, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(natom, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)

  allocate(xtilde(n, ndim, natom),mass(natom), label(natom))
  if (iproc.eq.0) then
     open(18, file="masses.dat", status="old")
     do j=1,natom
        read(18,*)label(j), mass(j)
     end do
     close(18)
  end if

  ! call MPI_Bcast(mpi_double_send, 5, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(beta, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(betan, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(dt, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(tau, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(gamma, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mass, natom, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(dHdrlimit, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(label, natom, MPI_CHARACTER, 0,MPI_COMM_WORLD, ierr)
  ndof=ndim*natom

  !Set up integration path and integration points/weights
  if (iproc .eq. 0) then
     !-------------------------
     !-------------------------
     !Read in initial wells, and masses
     allocate(initpath(ndim, natom),path(npath, ndim, natom), lampath(npath))
     allocate(origin(ndim))
     open(15, file="path.xyz")
     do i=1, npath
        read(15,*) dummy
        read(15,'(28A)') dummystr!, dummyE !'(A,G25.15)'
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
     end do
     lampath(:)= lampath(:)/lampath(npath)
     deallocate(initpath)
     close(15)
     open(20, file="aligned.xyz")
     do i=1,npath
        write(20,*) natom
        write(20,*) "Energy of minimum",i
        do j=1, natom
           write(20,*)  label(j), (path(i,k,j)*0.529177d0, k=1,ndim)
        end do
     end do
     close(20)
     !xunit=1 means bohr
     !xunit=2 means angstroms
     if (instapath) then
        allocate(well1(ndim,natom), well2(ndim,natom), wellinit(ndim,natom))
        open(15, file="well1.dat", status="old")
        open(16, file="well2.dat", status="old")
        do j=1,natom
           read(15,*) well1(1,j),well1(2,j),well1(3,j)
           read(16,*) well2(1,j),well2(2,j),well2(3,j)
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
        wellinit(:,:)= well2(:,:)
        call align_atoms(wellinit, theta1, theta2, theta3, origin, well2)

        write(*,*) "Potential at wells:", V(well1), V(well2)
        allocate(splinepath(npath))
        do i=1,ndim
           do j=1,natom
              splinepath(:)=0.0d0
              call spline(lampath(:), path(:,i,j), 1.0d31, 1.0d31, splinepath(:))
              do k=1,n
                 xtilde(k,i,j)= splint(lampath, path(:,i,j), splinepath(:), dble(k-1)/dble(n-1))
              end do
           end do
        end do
        deallocate(lampath,path, splinepath)
        call instanton(xtilde,well1,well2)
        npath=n
        allocate(lampath(npath), Vpath(npath), path(npath,ndim,natom))
        path(:,:,:)=xtilde(:,:,:)
        write(*,*) "Found instanton."
        open(19, file="instanton.xyz")
        do i=1,npath
           write(19,*) natom
           write(19,*) "Energy of minimum",i
           do j=1, natom
              write(19,*)  label(j), (xtilde(i,k,j)*0.529177d0, k=1,ndim)
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
     !-------------------------
     !Find the centre to make sure this is symmetric
     if (centre) then
        allocate(splinepath(npath))
        write(*,*) size(lampath), size(Vpath), size(splinepath)
        call spline(lampath, Vpath, 1.0d31, 1.0d31, splinepath)
        xmiddle= findmiddle(0.3d0, 0.7d0, lampath, Vpath, splinepath)
        deallocate(splinepath)
        a= 2.0d0  - 4.0d0*xmiddle
        b= 4.0d0*xmiddle - 1.0d0
        do i=1, npath
           ! lampath(i)= a*lampath(i)**2 + lampath(i)*b
           if (a.ge.0) then
              lampath(i)= -0.5d0*b/a + sqrt((lampath(i)/a) + (0.5*b/a)**2)
           else
              lampath(i)= -0.5d0*b/a - sqrt((lampath(i)/a) + (0.5*b/a)**2)
           end if
           ! write(*,*) i, dble(i-1)/dble(npath-1), lampath(i)
        end do
        write(*,*)"Centred,", a, b, xmiddle
     end if
     !-------------------------
     !-------------------------
     !set up gauss-legendre quadrature arrays
     allocate(weights(nintegral), xi(nintegral),xint(nintegral,ndim,natom))
     allocate(dbdxi(nintegral,ndim,natom))
     call gauleg(0d0,1.0d0,xi, weights,nintegral)

     allocate(y(n,ndim,natom), splinepath(npath))
     do i=1,ndim
        do j=1,natom
           splinepath(:)=0.0d0
           call spline(lampath(:), path(:,i,j), 1.0d31, 1.0d31, splinepath(:))
           do k=1,n
              ! y(k,i,j)=a*(dble(k)*1d-3)**2 + b*(dble(k)*1d-3) + c
              y(k,i,j)= splint(lampath, path(:,i,j), splinepath(:), dble(k-1)/dble(n-1))
           end do
           do k=1, nintegral
              xint(k,i,j)= splint(lampath, path(:,i,j), splinepath(:), xi(k))
              !a*xi(k)**2 + b*xi(k) + c
              dbdxi(k,i,j)= splin_grad(lampath, path(:,i,j), splinepath(:), xi(k))
           end do
        end do
     end do
     deallocate(splinepath,lampath)

     do i=1, nintegral
        write(*,*) xi(i), V(xint(i,:,:)), norm2(reshape(dbdxi(i,:,:), (/ndim*natom/)))
     end do

  end if

  !put a barrier to make sure the procs are synced
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  !-------------------------------------------
  !Set up end points array to distribute to the processors
  !ncalcs is number of calculations for each processor
  allocate(endpoints(ncalcs, ndim, natom), allendpoints(ncalcs*nproc, ndim, natom))
  allocate(gradpoints(ncalcs, ndim, natom), allgradpoints(ncalcs*nproc, ndim, natom))
  allocate(startpoint(ndim, natom))
  allendpoints= 0.0d0
  allgradpoints=0.0d0
  endpoints=0.0d0
  gradpoints=0.0d0
  startpoint=0.0d0
  if (iproc .eq. 0) then
     startpoint(:,:)= path(1,:,:)
     do i=1, nintegral
        do j=1, nrep
           do k=1,ndim
              do l=1, natom
                 allendpoints(nrep*(i-1) + j, k, l)= xint(i, k, l)
                 allgradpoints(nrep*(i-1) + j, k, l)= dbdxi(i, k, l)
              end do
           end do
        end do
     end do
     deallocate(path, xint, dbdxi)
  end if

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Bcast(startpoint, ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(xtilde, n*ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)

  if (iproc.eq.0) then
     do ii=1,nproc-1
        do i=1,ncalcs
           call MPI_Send(allendpoints(ncalcs*(ii) + i,:,:), ndof, MPI_DOUBLE_PRECISION, ii, 1, MPI_COMM_WORLD, ierr)
        end do
     end do
     endpoints(1:ncalcs, :, :) = allendpoints(1:ncalcs, :, :)
  else if (iproc.ne.0) then
     do i=1,ncalcs
        call MPI_Recv(endpoints(i,:,:), ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus,ierr)
     end do
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (iproc.eq.0) then
     do ii=1,nproc-1
        do i=1,ncalcs
           call MPI_Send(allgradpoints(ncalcs*(ii) + i,:,:), ndof, MPI_DOUBLE_PRECISION, ii, 1, MPI_COMM_WORLD, ierr)
        end do
     end do
     gradpoints(1:ncalcs, :, :) = allgradpoints(1:ncalcs, :, :)
  else if (iproc.ne.0) then
     do i=1,ncalcs
        call MPI_Recv(gradpoints(i,:,:), ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus,ierr)
     end do
  end if
  deallocate(allgradpoints, allendpoints)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  !-------------------------
  !-------------------------
  !Start loop over integration points
  allocate(integrand(ncalcs))
  integrand(:)=0.0d0
  ! call random_seed()
  allocate(x(n,ndim,natom), pinit(n,ndim,natom))
  allocate(transmatrix(n,n),beadmass(natom,n),beadvec(n,ndof), lam(n))
  call alloc_nm(iproc)
  do ii=1, ncalcs
     if (all(abs(endpoints(ii,:,:)) .lt. 1d-10)) then
        integrand(ii)=0.0d0
        cycle
     end if
     call init_nm(startpoint,endpoints(ii,:,:))
     call init_path(startpoint, endpoints(ii,:,:), x, pinit)
     if (thermostat .eq. 1) then
        call propagate_pimd_nm(x,pinit, startpoint,endpoints(ii,:,:), gradpoints(ii,:,:),dHdr)
     else if (thermostat .eq. 2) then
        call propagate_pimd_pile(x,pinit, startpoint,endpoints(ii,:,:), gradpoints(ii,:,:),dHdr)
     else if (thermostat .eq. 3) then
        call propagate_pimd_higher(x,pinit, startpoint,endpoints(ii,:,:), gradpoints(ii,:,:),dHdr)
     else
        write(*,*) "Incorrect thermostat option."
        stop
     end if
     integrand(ii)=dHdr/(betan**2)
     write(*,*) ii,iproc, integrand(ii)
  end do
  call free_nm()
  deallocate(x, pinit)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  !Collate data from all processors and have the root finalize the calculations
  allocate(allintegrands(ncalcs*nproc))
  
  call MPI_Gather(integrand, ncalcs, MPI_DOUBLE_PRECISION, allintegrands, ncalcs, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  if (iproc .eq. 0) then
     write(*,*) "final answers from all processors:"
     do i=1, ncalcs*nproc
        write(*,*) i, allintegrands(i)
     end do
     
     allocate(finalintegrand(nintegral), sigma(nintegral))
     finalintegrand=0.0d0
     sigma=0.0d0

     do ii=1,nintegral
        do jj=1,nrep
           finalintegrand(ii) = finalintegrand(ii) + allintegrands(nrep*(ii-1) + jj)
           sigma(ii)= sigma(ii)+ (allintegrands(nrep*(ii-1) + jj))**2
        end do
        finalintegrand(ii) = finalintegrand(ii)/dble(nrep)
        sigma(ii)= sigma(ii)/dble(nrep)
        sigma(ii)= sigma(ii) - finalintegrand(ii)**2
     end do

     answer= 0.0d0
     sigmaA=0.0d0

     do i=1, nintegral
        if (integrand(i) .eq. integrand(i)) &
             answer= answer+ weights(i)*finalintegrand(i)
        sigmaA= sigmaA + sigma(i)*weights(i)**2
     end do

     finalI= exp(-answer*betan)
     write(*,*) "beta, betan, n=", beta, betan, n
     write(*,*) "final Delta A=", answer , "+/-", sqrt(sigmaA)
     write(*,*) "q/q0=", finalI, betan*sqrt(sigmaA)*finalI
     write(*,*) "q0/q=", 1.0d0/finalI, betan*sqrt(sigmaA)/finalI
     call system_clock(time2, irate, imax)
     write(*,*) "Ran in", dble(time2-time1)/dble(irate), "s"
  end if
  deallocate(xtilde, label)

  call MPI_FINALIZE(ierr)

end program pimd
