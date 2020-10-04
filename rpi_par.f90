program rpi
  use mcmod_mass
  use verletint
  use instantonmod
  use parallelmod
  implicit none


  double precision, allocatable::   theta(:),phi(:),xharm(:,:,:)
  double precision, allocatable::   weightstheta(:),weightsphi(:), origin(:)
  double precision, allocatable::   eta(:),weightseta(:)
  double precision, allocatable::   etasquared(:),wellinit(:,:)
  double precision, allocatable::   initpath(:,:), xtilderot(:,:,:), wellrot(:,:)
  double precision::                lndetj, lndetj0, skink, psi, cutofftheta, cutoffphi
  double precision::                delta, omega, gammetilde, Ibeta
  double precision::                theta1, theta2, theta3
  integer::                         i, j,k,npath, dummy, zerocount, npoints
  integer::                         ii,jj,kk
  character::                      dummylabel, dummystr(28)
  !MPI variables
  integer::                        ierr, nproc, ncalcs
  integer, dimension(MPI_STATUS_SIZE) :: rstatus
  double precision, allocatable::  endpoints(:,:,:), allendpoints(:,:,:)
  double precision, allocatable::  results(:), allresults(:)
  logical::                        alignwell, angular

  namelist /RPIDATA/ n, beta, ndim, natom,npath,xunit, npoints, cutofftheta,cutoffphi, alignwell,&
       fixedends, angular

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
  n= 100
  beta= 30.0d0
  ndim=3
  natom=1
  xunit=1
  npath=0
  npoints=10
  cutofftheta=6.5d0
  cutoffphi=3.2d0
  alignwell=.false.
  fixedends=.true.

  call V_init(iproc)

  if (iproc .eq. 0) then
     read(5, nml=RPIDATA)
     betan= beta/dble(n)
     write(*,*)"Running with parameters (in a.u.):"
     write(*,*) "beta, betan, n=", beta, betan, n
     if (angular) then
        ncalcs= (npoints**3)/nproc
        if (mod(npoints**3, nproc) .ne. 0) ncalcs=ncalcs+1
     end if
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  ierr=0
  if (angular) then
     call MPI_Bcast(ncalcs, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
     call MPI_Bcast(npoints, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  end if
  call MPI_Bcast(N, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ndim, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(natom, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fixedends, 1, MPI_LOGICAL, 0,MPI_COMM_WORLD, ierr)
  ndof=ndim*natom
  totdof= n*ndof
  if ((iproc.eq. nproc-1 .and. mod(npoints**3,nproc) > 0) .and. angular) then
     ncalcs= ncalcs-mod(npoints**3,nproc)
  end if
  
    if (.not. angular) then
     ncalcs= N/nproc
     if (iproc .lt. mod(N, nproc)) ncalcs=ncalcs+1
  end if


  write(*,*) "ncalcs=", ncalcs, "on iproc", iproc

  allocate(mass(natom), label(natom), xtilde(n, ndim, natom))
  allocate(well1(ndim,natom),well2(ndim,natom), wellinit(ndim,natom))

  if (iproc.eq.0) then
     open(18, file="masses.dat", status="old")
     do j=1,natom
        read(18,*)label(j), mass(j)
     end do
     close(18)
     open(30, file="well1.dat", status='old')
     open(40, file="well2.dat", status='old')
     do j=1, natom
        read(30,*)(well1(i,j), i=1,ndim)
        read(40,*)(well2(i,j), i=1,ndim)
     end do
     close(30)
     close(40)

     if (xunit .eq. 2) then
        well1(:,:)= well1(:,:)/0.529177d0
        well2(:,:)= well2(:,:)/0.529177d0
     end if
     allocate(origin(ndim))

     call get_align(well1,theta1, theta2, theta3, origin)
     wellinit(:,:)= well1(:,:)
     call align_atoms(wellinit, theta1, theta2, theta3, origin, well1)
     if (alignwell) call get_align(well2,theta1, theta2, theta3, origin)
     wellinit(:,:)= well2(:,:)
     call align_atoms(wellinit, theta1, theta2, theta3, origin, well2)
     V0=V(well1)
     write(*,*) "Potential at wells:", V(well1), V(well2)

     !-------------------------
     !obtain instanton solution, x_tilde
     call read_path(.false.,.false.,npath,path,Vpath,splinepath,lampath)
     xtilde=0.0d0
     do i=1,ndim
        do j=1,natom
           do k=1,n
              xtilde(k,i,j)= splint(lampath, path(:,i,j), splinepath(:,i,j), dble(k-1)/dble(n-1))
           end do
        end do
     end do
     deallocate(Vpath, path, splinepath,lampath)
  end if


  call MPI_Bcast(beta, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(betan, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(V0, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mass, natom, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(well1, ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(well2, ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(label, natom, MPI_CHARACTER, 0,MPI_COMM_WORLD, ierr)

  write(*,*) "Starting instanton search..."
  call parallel_instanton(xtilde,iproc, nproc,well1,well2)
  if (iproc .eq. 0) then
     write(*,*) "Found instanton."
     open(19, file="instanton.xyz")
     do i=1,n
        write(19,*) natom
        write(19,*) "Energy of minimum",i
        do j=1, natom
           write(19,*)  label(j), (xtilde(i,k,j)*0.529177d0, k=1,ndim)
        end do
     end do
     close(19)
  end if
  call MPI_Bcast(xtilde, totdof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)

  !------------------------------------
  !work out Q_0
  allocate( etasquared(totdof))
  if (iproc .eq. 0) then
     allocate(xharm(n,ndim, natom))
     do i=1,n
        xharm(i,:,:)= well1(:,:)
     end do
     etasquared(:)=0.0d0
     call detJ(xharm,etasquared, .true.)
     lndetj0= 0.0d0
     zerocount=0
     do i=1,totdof
        if (etasquared(i) .gt. 0.0d0) then
           ! write(*,*) i,etasquared(i)
           lndetj0= lndetj0+ log(etasquared(i))
        else
           zerocount=zerocount+1
        end if
     end do
     write(*,*) "Skipped ", zerocount, "states"
     write(*,*) "lndetJ0", lndetj0
     deallocate(xharm)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Bcast(lndetj0, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  !put a barrier to make sure the procs are synced

  !------------------------------------
  !loop over solid angle points
  !-------------------------------------------
  !Set up end points array to distribute to the processors
  !ncalcs is number of calculations for each processor
  if (angular) then
     allocate(wellrot(ndim,natom))
     allocate(endpoints(ncalcs, ndim, natom), allendpoints(npoints**3, ndim, natom))
     allendpoints=0.0d0
     if (iproc .eq. 0) then
        allocate(theta(npoints), weightstheta(npoints))
        allocate(phi(npoints), weightsphi(npoints))
        allocate(eta(npoints), weightseta(npoints))
        call gauleg(0d0,min(2.0d0*pi, cutofftheta),theta, weightstheta,npoints)
        call gauleg(0d0,min(pi, cutoffphi),phi, weightsphi,npoints)
        call gauleg(0d0,min(2.0d0*pi, cutofftheta),eta, weightseta,npoints)
        i=1
        do ii=1,npoints
           do jj=1,npoints
              do kk=1,npoints
                 wellrot(:,:)= well2(:,:)
                 call rotate_atoms(wellrot,1,eta(kk))
                 call rotate_atoms(wellrot,3,phi(jj))
                 call rotate_atoms(wellrot,1,theta(ii))
                 allendpoints(i, :,:)= wellrot(:,:)
                 i=i+1
              end do
           end do
        end do
     end if
     endpoints=0.0d0

     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     write(*,*) iproc, "passed barrier"

     if (iproc.eq.0) then
        do ii=1,nproc-1
           do i=1,ncalcs
              call MPI_Send(allendpoints(ncalcs*ii + i,:,:), ndof, MPI_DOUBLE_PRECISION, ii, 1, MPI_COMM_WORLD, ierr)
           end do
        end do
        endpoints(1:ncalcs, :, :) = allendpoints(1:ncalcs, :, :)
     else
        do i=1,ncalcs
           call MPI_Recv(endpoints(i,:,:), ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus,ierr)
        end do
     end if
     deallocate(allendpoints)
     allocate(xtilderot(n,ndim,natom), results(ncalcs))
     do ii=1, ncalcs
        xtilderot(:,:,:)= xtilde(:,:,:)
        call instanton(xtilderot,well1,endpoints(ii,:,:))
        call detJ(xtilderot, etasquared, .false.)

        lndetj= 0.0d0
        zerocount=0
        do i=2,totdof
           if (etasquared(i) .gt. 0.0d0) then
              lndetj= lndetj+ log(etasquared(i))
           else
              zerocount=zerocount+1
           end if
        end do
        write(*,*) "Skipped ", zerocount, "states"
        write(*,*) "lndetJ", lndetj
        gammetilde= exp(0.5d0*(lndetJ-lndetJ0))
        !-------------------------
        !Put it all together and output.
        Skink= betan*UM(xtilderot, well1, endpoints(ii,:,:))
        omega= betan*exp(-skink)*sqrt(skink/(2.0d0*pi))/gammetilde
        if (omega*N > 1.0d0) then
           Ibeta=1.0
        else
           Ibeta= tanh(omega*N)
        end if
        results(ii)= Ibeta
        write(*,*) iproc, ii, Ibeta, omega, skink
     end do
     write(*,*) "Processor", iproc, "is done."
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     !Collate data from all processors and have the root finalize the calculations
     allocate(allresults(ncalcs*nproc))

     call MPI_Gather(results, ncalcs, MPI_DOUBLE_PRECISION, allresults, ncalcs, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

     if (iproc .eq. 0) then
        open(90, file="angularI.dat")
        i=1
        do ii=1,npoints
           do jj=1,npoints
              do kk=1,npoints
                 write(90,*) theta(ii), phi(jj), eta(kk),weightstheta(ii)*weightsphi(jj)*weightseta(kk), allresults(i)
                 i=i+1
              end do
           end do
        end do
        close(90)
        deallocate(theta, phi, eta)
        deallocate(weightstheta,weightsphi, weightseta)
        deallocate(well2)
     end if
     deallocate(xtilderot, endpoints, well1)
  else !if not looping over angles
       !-------------------------
     !work out Hessian at the instantons
     !-------------------------
     !Put it all together and output.
     call parallel_detJ(xtilde, iproc, nproc, etasquared)
     if (iproc .eq. 0) then
        lndetj= 0.0d0
        zerocount=0
        write(*,*) 1,etasquared(1)
        do i=2,totdof
           ! write(*,*) i,etasquared(i)
           if (etasquared(i) .gt. 0.0d0) then
              lndetj= lndetj+ log(etasquared(i))
           else
              zerocount=zerocount+1
           end if
        end do
        write(*,*) "Skipped ", zerocount, "states"
        write(*,*) "lndetJ", lndetj, lndetj0
        gammetilde= exp(0.5d0*(lndetJ-lndetJ0))
     end if
     if (fixedends) then
        Skink= betan*parallel_UM(xtilde, iproc,nproc,well1, well2)
     else
        Skink= betan*parallel_UM(xtilde,iproc,nproc)
     end if
     if (iproc .eq. 0) then
        omega= betan*exp(-skink)*sqrt(skink/(2.0d0*pi))/gammetilde
        delta= 2.0d0*omega/betan
        write(*,*) "phi=", gammetilde
        write(*,*) "s_kink=", skink
        write(*,*) "theta, N*theta=",omega,omega*N
        write(*,*) "h_ij=", omega/betan, 219475.0d0*omega/betan
        write(*,*) "RPI delta=", delta, delta*219475.0d0
     end if
  end if

  deallocate(etasquared, xtilde)
  call MPI_FINALIZE(ierr)
end program
