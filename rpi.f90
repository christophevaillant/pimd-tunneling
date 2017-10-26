program rpi
  use mcmod_mass
  implicit none
  double precision, allocatable::   xtilde(:,:,:), xharm(:,:,:), H(:,:)
  double precision, allocatable::   HHarm(:,:), etasquared(:),Vpath(:)
  double precision, allocatable::  path(:,:,:), lampath(:), splinepath(:)
  double precision, allocatable::   initpath(:,:), xtilderot(:,:,:), wellrot(:,:)
  double precision::                lndetj, lndetj0, skink, psi, cutoff
  double precision::                theta,delta, phi, omega, gamma, Ibeta
  integer::                         i, j,k,npath, dummy, zerocount, npoints
  integer::                         ii,jj,kk
  character, allocatable::         label(:)
  character::                      dummylabel, dummystr(28)
  logical::                        angular
  namelist /RPIDATA/ n, beta, ndim, natom,npath,xunit, eps, angular, npoints, cutoff

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
  angular=.false.

  read(5, nml=RPIDATA)
  betan= beta/dble(n)
  allocate(mass(natom),label(natom))
  open(18, file="masses.dat", status="old")
  do j=1,natom
     read(18,*)label(j), mass(j)
  end do
  close(18)

  ndof=ndim*natom
  totdof= n*ndof
  call V_init()

  allocate(well1(ndim,natom),well2(ndim,natom))
  open(30, file="well1.dat", status='old')
  open(40, file="well2.dat", status='old')
  do j=1, natom
     read(30,*)(well1(i,j), i=1,ndim)
     read(40,*)(well2(i,j), i=1,ndim)
  end do
  close(30)
  close(40)
  write(*,*) "Potential at wells:", V(well1), V(well2)
  !-------------------------
  !obtain instanton solution, x_tilde
  allocate(initpath(ndim, natom), lampath(npath), Vpath(npath))
  allocate(path(npath, ndim, natom))
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
     else
        lampath(i)= lampath(i-1) + eucliddist(path(i-1,:,:), path(i,:,:))!dble(i-1)/dble(npath-1)
     end if
     ! call align_atoms(initpath, path(i,:,:))
     path(i,:,:)= initpath(:,:)
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

  deallocate(initpath)

  if (xunit .eq. 2) then
     path(:,:,:) = path(:,:,:)/0.529177d0
  end if

  allocate(xtilde(n, ndim, natom),splinepath(npath))
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

  !------------------------------------
  !work out Q_0
  allocate(xharm(n,ndim, natom), etasquared(totdof))
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

  call detJ(xtilde, etasquared)
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
  write(*,*) "lndetJ", lndetj, lndetj0
  gamma= exp(0.5d0*(lndetJ-lndetJ0))

  !------------------------------------
  !loop over solid angle points
  if (angular) then
     allocate(xtilderot(n,ndim,natom), wellrot(ndim,natom))
     open(90, file="angularI.dat")
     xtilderot(:,:,:)= xtilde(:,:,:)
     do ii=1,npoints
        theta= cutoff*dble(ii-1)/dble(npoints-1)
        do jj=1,npoints
           ! do k=1,10
           phi= 0.5d0*cutoff*dble(jj-1)/dble(npoints-1)
              ! do i=1,n
              !    dtheta= dble(i-1)*theta/dble(n-1)
              !    ! dphi= dble(i-1)*phi/dble(n-1)
              !    call rotate_atoms(xtilderot(i,:,:),1,dtheta)
              !    ! call rotate_atoms(xtilderot(i,:,:),1,dphi)
              ! end do
              wellrot(:,:)= well2(:,:)
              call rotate_atoms(wellrot,1,theta)
              call rotate_atoms(wellrot,2,phi)
              ! call rotate_atoms(wellrot,3,psi)
              call instanton(xtilderot,well1,wellrot)
              ! call detJ(xtilderot, etasquared)

              ! lndetj= 0.0d0
              ! zerocount=0
              ! do i=2,totdof
              !    if (etasquared(i) .gt. 0.0d0) then
              !       lndetj= lndetj+ log(etasquared(i))
              !    else
              !       zerocount=zerocount+1
              !    end if
              ! end do
              ! write(*,*) "Skipped ", zerocount, "states"
              ! write(*,*) "lndetJ", lndetj, lndetj0
              !-------------------------
              !Put it all together and output.
              Skink= betan*UM(xtilderot, well1, wellrot)
              omega= betan*exp(-skink)*sqrt(skink/(2.0d0*pi))/gamma
              Ibeta= tanh(omega*N)
              write(90,*) theta, phi, Ibeta, 219475.0d0*omega/betan
              write(*,*) theta,phi, Ibeta, 219475.0d0*omega/betan
              ! end do
           ! end do
        end do
     end do
     deallocate(xtilderot, wellrot)
     close(90)
  else
  !-------------------------
  !work out Hessian at the instantons
     !-------------------------
     !Put it all together and output.
     gamma= exp(0.5d0*(lndetJ-lndetJ0))
     Skink= betan*UM(xtilde, well1, well2)
     omega= betan*exp(-skink)*sqrt(skink/(2.0d0*pi))/gamma
     delta= 2.0d0*omega/betan
     write(*,*) "beta=", beta, "n=", n
     write(*,*) "beta_n=", betan
     write(*,*) "phi=", gamma
     write(*,*) "s_kink=", skink
     write(*,*) "theta, N*theta=",omega,omega*N
     write(*,*) "h_ij=", omega/betan, 219475.0d0*omega/betan
     write(*,*) "RPI delta=", delta, delta*219475.0d0
  end if

  deallocate(xharm, etasquared, xtilde)
end program rpi
