program pimd
  use mcmod_mass
  use verletint
  use MKL_VSL_TYPE
  use MKL_VSL
  use instantonmod
  implicit none

  ! include 'mpif.h'
  ! include 'mkl_service.fi'
  !general variables
  integer::                        i, j, irej, randind,count,k,l,jj, nrep
  integer::                        thermostat,dofi, npath, dummyint
  double precision::               Isum, Isum0,Isqprev, xcutoff, UHnew, UHold
  double precision::               diagonals, offdiags, Pacc, delta, stdev
  double precision::               UHinitial, f, randno,answer, com, totm
  double precision::               lndetj, lndetj0, skink, theta, phi
  double precision::               dHdr, dx, finalI,sigmaA, a, b, xmiddle
  double precision::               theta1, theta2, theta3
  double precision, allocatable::  x(:,:,:), initpath(:,:)
  double precision, allocatable::  origin(:), wellinit(:,:)
  double precision, allocatable::  xi(:), dbdxi(:,:,:), pinit(:,:,:)
  character(len=30)::              format_string, restartstr
  logical::                        instapath, centre, readpath, alignwell
  character(len=3)::               dummychar, iprocchar, iichar
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
  double precision, allocatable::  endpoints(:,:,:), allendpoints(:,:,:)
  double precision, allocatable::  gradpoints(:,:,:), allgradpoints(:,:,:), startpoint(:,:)
  double precision, allocatable::  finalintegrand(:), allintegrands(:)
  double precision, allocatable::  allxipoints(:), xipoints(:)
  integer (kind=4)::               ierrmkl
  integer, dimension(MPI_STATUS_SIZE) :: rstatus
  namelist /MCDATA/ n, beta, NMC, Noutput,dt, iprint,imin,tau,npath, gamma,&
       nintegral,nrep, use_mkl, thermostat, ndim, natom, xunit, instapath, centre,&
       dHdrlimit, readpath, alignwell, fixedends, cayley, readhess,basename, &
       atom1,atom2,atom3, restart

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
  centre=.false.
  dHdrlimit=-1.0
  readpath=.true.
  alignwell=.false.
  fixedends=.true.
  cayley=.false.
  readhess=.false.
  atom1=1
  atom2=2
  atom3=3
  basename=""
  restart=0
  
  !Read in namelist variables for root proc, and spread
  !it to other procs
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
  call MPI_Bcast(atom1, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(atom2, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(atom3, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(basename, 20, MPI_CHARACTER, 0,MPI_COMM_WORLD, ierr)
  call V_init(iproc)

  if (iproc .eq. 0) then
     !Read in wells and work out the instanton path in order to interpolate
     allocate(origin(ndim))
     allocate(well1(ndim,natom), well2(ndim,natom), wellinit(ndim,natom))
     open(15, file="well1.dat", status="old")
     open(16, file="well2.dat", status="old")
     do j=1,natom
        read(15,*) (well1(i,j),i=1,ndim) !well1(2,j),well1(3,j)
        read(16,*) (well2(i,j),i=1,ndim) !well2(1,j),well2(2,j),well2(3,j)
     end do
     close(15)
     close(16)
     ! xunit=1 means bohr
     ! xunit=2 means angstroms
     if (xunit .eq. 2) then
        well1(:,:)= well1(:,:)/0.529177d0
        well2(:,:)= well2(:,:)/0.529177d0
     end if
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  ierr=0
  call MPI_Bcast(ncalcs, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(N, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(NMC, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Noutput, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(thermostat, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(imin, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ndim, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(natom, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(restart, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fixedends, 1, MPI_LOGICAL, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(readhess, 1, MPI_LOGICAL, 0,MPI_COMM_WORLD, ierr)

  allocate(xtilde(n, ndim, natom),mass(natom), label(natom))
  if (iproc.eq.0) then
     open(18, file="masses.dat", status="old")
     do j=1,natom
        read(18,*)label(j), mass(j)
     end do
     close(18)
     call get_align(well1,theta1, theta2, theta3, origin)

     wellinit(:,:)= well1(:,:)
     call align_atoms(wellinit, theta1, theta2, theta3, origin, well1)
     if (alignwell) call get_align(well2,theta1, theta2, theta3, origin)
     wellinit(:,:)= well2(:,:)
     call align_atoms(wellinit, theta1, theta2, theta3, origin, well2)
     V0=V(well1)
     write(*,*) "Potential at wells:", V(well1), V(well2)
     deallocate(wellinit)

     if (alignwell) then
        open(25, file="aligned_wells.xyz")
        write(25,*) natom
        write(25,*) "rotation angles:", theta1, theta2, theta3
        do j=1, natom
           write(25,*)  label(j), (well1(k,j)*0.529177d0, k=1,ndim)
        end do
        write(25,*) natom
        write(25,*) "rotation angles:", theta1, theta2, theta3
        do j=1, natom
           write(25,*)  label(j), (well2(k,j)*0.529177d0, k=1,ndim)
        end do
        close(25)
     end if
  end if

  call MPI_Bcast(beta, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(betan, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(dt, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(tau, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(gamma, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mass, natom, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(dHdrlimit, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(label, natom, MPI_CHARACTER, 0,MPI_COMM_WORLD, ierr)
  ndof=ndim*natom
  totdof=n*ndim*natom
  !Set up integration path and integration points/weights
  if (iproc .eq. 0) then
     !-------------------------
     !-------------------------
     !Read in initial wells, and masses
     if (readpath) then
        call read_path(instapath,centre,npath,path,Vpath,splinepath,lampath)
        if (readhess) then
           call read_hess(npath,path,lampath,splinehess,hesspath)
        end if
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
              xint(k,i,j)= splint(lampath, path(:,i,j), splinepath(:,i,j), xi(k))
              dbdxi(k,i,j)= splin_grad(lampath, path(:,i,j), splinepath(:,i,j), xi(k))
           end do
        end do
     end do

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
  allocate(startpoint(ndim, natom), allxipoints(ncalcs*nproc), xipoints(ncalcs))
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
           allxipoints(nrep*(i-1)+j)=xi(i)
        end do
     end do
     deallocate(xint, dbdxi)
  end if

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Bcast(startpoint, ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(xtilde, n*ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(V0, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(npath, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  if (iproc .ne. 0) allocate(path(npath,ndim,natom), lampath(npath), splinepath(npath,ndim,natom))
  call MPI_Bcast(path, npath*ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(lampath, npath, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(splinepath, npath*ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  if (readhess) then
     if (iproc .ne. 0) allocate(hesspath(npath,ndof,ndof), splinehess(npath,ndof,ndof))
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     !need to split up data transfer because of bullshit mpi reasons
     do i=1,ndof
        do j=1,ndof
           call MPI_Bcast(splinehess(:,i,j), npath, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
           call MPI_Bcast(hesspath(:,i,j), npath, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
        end do
     end do
  end if

  if (iproc.eq.0) then
     do ii=1,nproc-1
        do i=1,ncalcs
           call MPI_Send(allendpoints(ncalcs*(ii) + i,:,:), ndof, MPI_DOUBLE_PRECISION, ii, 1, MPI_COMM_WORLD, ierr)
           call MPI_Send(allxipoints(ncalcs*(ii) + i), 1, MPI_DOUBLE_PRECISION, ii, 1, MPI_COMM_WORLD, ierr)
        end do
     end do
     endpoints(1:ncalcs, :, :) = allendpoints(1:ncalcs, :, :)
     xipoints(:)= allxipoints(1:ncalcs)
  else if (iproc.ne.0) then
     do i=1,ncalcs
        call MPI_Recv(endpoints(i,:,:), ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus,ierr)
        call MPI_Recv(xipoints(i), 1, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus,ierr)
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
     !need to check status of calculation for restart
     if (restart .lt. 2) then
        call init_path(xipoints(ii), x, pinit)
        restartnmc=0
     end if
     if (restart .ge. 1) then
        restartunit= 600+ iproc
        !-----------------------
        !format statement for reading in
        if (iproc .lt. 10) then
           format_string = "(A12,I1,A1"
        else if (iproc .ge. 10 .and. iproc .lt. 100) then
           format_string = "(A12,I2,A1"
        else if (iproc .ge. 100) then
           format_string = "(A12,I3,A1"
        endif
        
        if (ii .lt. 10) then
           format_string= trim(format_string) // ",I1,A4)"
        else if (ii .ge. 10 .and. ii .lt. 100) then
           format_string= trim(format_string) // ",I2,A4)"
        else if (ii .ge. 100 .and. ii .lt. 1000) then
           format_string= trim(format_string) // ",I3,A4)"
        end if
        
        write(restartstr,trim(format_string)) "restart_proc", iproc, "_", ii,".xyz"
        open(restartunit,file=trim(restartstr))
     end if
     if (restart .eq. 2) then
        do i=1,n
           read(restartunit,*) dummyint
           read(restartunit,*) dHdr
           do j=1,natom
              read(restartunit,*) dummychar, x(i,1,j), x(i,2,j), x(i,3,j)
           end do
        end do
        do i=1,n
           read(restartunit,*) dummyint
           read(restartunit,*) restartnmc
           do j=1,natom
              read(restartunit,*) dummychar, pinit(i,1,j), pinit(i,2,j), pinit(i,3,j)
           end do
        end do
     end if
     if (thermostat .eq. 1) then
        call propagate_pimd_nm(x,pinit, startpoint,endpoints(ii,:,:), gradpoints(ii,:,:),dHdr)
     else if (thermostat .eq. 2) then
        call propagate_pimd_pile(x,pinit, startpoint,endpoints(ii,:,:), gradpoints(ii,:,:),xipoints(ii),dHdr)
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
     write(*,*) "final Delta A=", answer , "+/-", sqrt(sigmaA),sqrt(sigmaA/nrep)
     write(*,*) "q/q0=", finalI, betan*sqrt(sigmaA)*finalI,betan*sqrt(sigmaA/nrep)*finalI
     write(*,*) "q0/q=", 1.0d0/finalI, betan*sqrt(sigmaA)/finalI, betan*sqrt(sigmaA/nrep)/finalI
     call system_clock(time2, irate, imax)
     write(*,*) "Ran in", dble(time2-time1)/dble(irate), "s"
  end if
  deallocate(xtilde, label)

  call MPI_FINALIZE(ierr)

end program pimd
