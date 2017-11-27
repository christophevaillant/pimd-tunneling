program rpi
  use mcmod_mass
  use verletint
  implicit none

  include 'mpif.h'

  double precision, allocatable::   theta(:),phi(:), xtilde(:,:,:), xharm(:,:,:), H(:,:)
  double precision, allocatable::   weightstheta(:),weightsphi(:), origin(:)
  double precision, allocatable::   eta(:),weightseta(:)
  double precision, allocatable::   HHarm(:,:), etasquared(:),Vpath(:), wellinit(:,:)
  double precision, allocatable::  path(:,:,:), lampath(:), splinepath(:)
  double precision, allocatable::   initpath(:,:), xtilderot(:,:,:), wellrot(:,:)
  double precision::                lndetj, lndetj0, skink, psi, cutofftheta, cutoffphi
  double precision::                delta, omega, gammetilde, Ibeta
  double precision::                theta1, theta2, theta3
  integer::                         i, j,k,npath, dummy, zerocount, npoints
  integer::                         ii,jj,kk
  character, allocatable::         label(:)
  character::                      dummylabel, dummystr(28)
  !MPI variables
  integer::                        ierr, nproc, ncalcs
  integer, dimension(MPI_STATUS_SIZE) :: rstatus
  double precision, allocatable::  endpoints(:,:,:), allendpoints(:,:,:)
  double precision, allocatable::  results(:), allresults(:)

  namelist /RPIDATA/ n, beta, ndim, natom,npath,xunit, eps, npoints, cutofftheta,cutoffphi

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
  eps=1d-3
  npoints=10
  call V_init()

  if (iproc .eq. 0) then
     read(5, nml=RPIDATA)
     betan= beta/dble(n)
     write(*,*)"Running with parameters (in a.u.):"
     write(*,*) "beta, betan, n=", beta, betan, n
     ncalcs= (npoints**3)/nproc
     if (mod(npoints**3, nproc) .ne. 0) ncalcs=ncalcs+1

  end if

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  ierr=0
  call MPI_Bcast(ncalcs, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(N, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ndim, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(natom, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  ndof=ndim*natom
  totdof= n*ndof

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
     wellinit(:,:)= well2(:,:)
     call align_atoms(wellinit, theta1, theta2, theta3, origin, well2)

     write(*,*) "Potential at wells:", V(well1), V(well2)

     !-------------------------
     !obtain instanton solution, x_tilde
     allocate(initpath(ndim, natom), lampath(npath), Vpath(npath))
     allocate(path(npath, ndim, natom),splinepath(npath))
     path=0.0d0
     open(15, file="path.xyz")
     do i=1, npath
        read(15,*) dummy
        read(15,'(28A)') dummystr!, dummyE !'(A,G25.15)'
        do j=1, natom
           read(15,*) dummylabel, (initpath(k,j), k=1,ndim)
        end do
        if (i.eq.1) then
           ! call align_atoms(
           lampath(1)=0.0d0
           call get_align(initpath,theta1, theta2, theta3, origin)
        else
           lampath(i)= lampath(i-1) + eucliddist(path(i-1,:,:), path(i,:,:))!dble(i-1)/dble(npath-1)
        end if
        call align_atoms(initpath,theta1, theta2, theta3, origin, path(i,:,:))
        ! path(i,:,:)= initpath(:,:)
        Vpath(i)= V(path(i,:,:))
     end do
     lampath(:)= lampath(:)/lampath(npath)
     close(15)
     open(20, file="aligned.xyz")
     do i=1,npath
        write(20,*) natom
        write(20,*) "Energy of minimum",i
        do j=1, natom
           if (xunit .eq. 1) then
              write(20,*)  label(j), (path(i,k,j)*0.529177d0, k=1,ndim)
           else
              write(20,*)  label(j), (path(i,k,j), k=1,ndim)
           end if
        end do
     end do
     close(20)

     deallocate(initpath, origin)

     if (xunit .eq. 2) then
        path(:,:,:) = path(:,:,:)/0.529177d0
     end if

     xtilde=0.0d0
     splinepath=0.0d0
     do i=1,ndim
        do j=1,natom
           splinepath(:)=0.0d0
           call spline(lampath(:), path(:,i,j), 1.0d31, 1.0d31, splinepath(:))
           do k=1,n
              xtilde(k,i,j)= splint(lampath, path(:,i,j), splinepath(:), dble(k-1)/dble(n-1))
           end do
        end do
     end do
     deallocate(lampath,Vpath, path, splinepath)
     call instanton(xtilde,well1,well2)
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

  call MPI_Bcast(eps, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(beta, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(betan, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mass, natom, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(well1, ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(well2, ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
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
     call detJ(xharm, etasquared)
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
     deallocate(xharm)
  end if
  call MPI_Bcast(lndetj0, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  !put a barrier to make sure the procs are synced
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  !------------------------------------
  !loop over solid angle points
  !-------------------------------------------
  !Set up end points array to distribute to the processors
  !ncalcs is number of calculations for each processor
  allocate(xtilderot(n,ndim,natom), wellrot(ndim,natom))
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
              call rotate_atoms(wellrot,2,phi(jj))
              call rotate_atoms(wellrot,1,theta(ii))
              allendpoints(i, :,:)= wellrot(:,:)
              i=i+1
           end do
        end do
     end do
  end if

  endpoints=0.0d0
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  if (iproc.eq.0) then
     do ii=1,nproc-1
        do i=1,ncalcs
           call MPI_Send(allendpoints(ncalcs*ii + i,:,:), ndof, MPI_DOUBLE_PRECISION, ii, 1, MPI_COMM_WORLD, ierr)
        end do
     end do
     endpoints(1:ncalcs, :, :) = allendpoints(1:ncalcs, :, :)
  else if (iproc.ne.0) then
     do i=1,ncalcs
        call MPI_Recv(endpoints(i,:,:), ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus,ierr)
     end do
  end if

  allocate(results(ncalcs))
  do ii=1, ncalcs
     xtilderot(:,:,:)= xtilde(:,:,:)
     call instanton(xtilderot,well1,endpoints(ii,:,:))
     call detJ(xtilderot, etasquared)

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
     Ibeta= tanh(omega*N)
     results(ii)= Ibeta
     write(*,*) iproc, ii, Ibeta
  end do
  deallocate(xtilderot)
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
  end if

  deallocate(etasquared, xtilde)
  call MPI_FINALIZE(ierr)
end program rpi