program rpi
  use mcmod_mass
  use verletint
  use instantonmod
  implicit none
  double precision, allocatable::   theta(:),phi(:), xharm(:,:,:), H(:,:,:,:)
  double precision, allocatable::   weightstheta(:),weightsphi(:), origin(:)
  double precision, allocatable::   eta(:),weightseta(:)
  double precision, allocatable::   HHarm(:,:), etasquared(:), wellinit(:,:), grad(:,:)
  double precision, allocatable::   initpath(:,:), xtilderot(:,:,:), wellrot(:,:)
  double precision::                lndetj, lndetj0, skink, psi, cutofftheta, cutoffphi
  double precision::                delta, omega, gammetilde, Ibeta, Vwell1, Vwell2
  double precision::                theta1, theta2, theta3
  integer::                         i, j,k,npath, dummy, zerocount, npoints
  integer::                         ii,jj,kk,l
  character::                      dummylabel, dummystr(28)
  logical::                        angular, output_instanton, readpath, alignwell, alignpath, checkhess
    integer::                        nout, ldz, lwork, liwork, info
  integer,allocatable::            isuppz(:), iwork(:)
  double precision, allocatable::  work(:)

  namelist /RPIDATA/ n, beta, ndim, natom,npath,xunit, angular, npoints, cutofftheta,cutoffphi,&
       output_instanton,readpath, alignwell, fixedends, alignpath, checkhess

  !-------------------------
  !Set default system parameters then read in namelist
  n= 100
  beta= 30.0d0
  ndim=3
  natom=1
  xunit=1
  npath=0
  npoints=10
  angular=.false.
  output_instanton=.false.
  readpath=.true.
  alignwell=.false.
  fixedends=.true.
  alignpath=.true.
  potforcepresent=.false.
  checkhess=.false.
  
  read(5, nml=RPIDATA)
  betan= beta/dble(n)
  allocate(mass(natom),label(natom))
  open(18, file="masses.dat", status="old")
  do j=1,natom
     read(18,*)label(j), mass(j)
  end do
  close(18)

  if (angular) fixedends=.true.
  ndof=ndim*natom
  totdof= n*ndof

  call V_init(0)
  allocate(well1(ndim,natom),well2(ndim,natom), wellinit(ndim,natom))
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

  if (alignwell) then
     call get_align(well1,theta1, theta2, theta3, origin)
     wellinit(:,:)= well1(:,:)
     call align_atoms(wellinit, theta1, theta2, theta3, origin, well1)
     call get_align(well2,theta1, theta2, theta3, origin)
     wellinit(:,:)= well2(:,:)
     call align_atoms(wellinit, theta1, theta2, theta3, origin, well2)
     open(21, file="aligned_ends.xyz")
     write(21,*) natom
     write(21,*) "Well1"
     do i=1, natom
        write(21,*)  label(i), (well1(k,i)*0.529177d0, k=1,ndim)
     end do
     write(21,*) natom
     write(21,*) "Well2"
     do i=1, natom
        write(21,*)  label(i), (well2(k,i)*0.529177d0, k=1,ndim)
     end do
     close(21)
  end if
  ! write(*,*) "Potential at wells:", V(well1), V(well2), V0
  V0=V(well1)
  Vwell1= 0.0d0
  Vwell2= V(well2)
  write(*,*) "Potential at wells:", Vwell1, Vwell2, V0
  allocate(grad(ndim,natom))
  call Vprime(well1, grad)
  write(*,*) "With norm of grad:", norm2(reshape(grad, (/ndim*natom/)))
  write(*,*) "beta=", beta, "n=", n
  write(*,*) "beta_n=", betan

  if (checkhess) then
     write(*,*)"Calculating Hessian as test"
     allocate(H(ndim,natom,ndim,natom), etasquared(ndof))
     call Vdoubleprime(well1,H)
     open(22, file="hessian_well1.dat")
     do j=1,natom
        do i=1,ndim
           write(22,*) ((H(i,j,k,l), k=1,ndim),l=1,natom)           
        end do
        do l=1,natom
           H(:,j,:,l) = H(:,j,:,l)/sqrt(mass(j)*mass(l))
        end do
     end do
     ! !get diagonal hessian
     lwork= 2*ndof+1
     liwork= 1
     info=0
     allocate(work(lwork), iwork(liwork))
     call DSYEVD('N', 'L', ndof, reshape(H, (/ndof,ndof/)), ndof, etasquared, work, lwork, iwork, liwork, info)
     write(*,*)"normal modes:"
     do j=1, ndof
        write(*,*) j,etasquared(j)*219475.0d0
     end do
     deallocate(H,etasquared)
  end if
  !-------------------------
  !obtain instanton solution, x_tilde
  allocate(xtilde(n, ndim, natom))
  if (readpath) then
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

  else
     !do a quick and dirty linear interpolation
     write(*,*) "well1"
     do j=1, natom
        write(*,*)  label(j), (well1(k,j)*0.529177d0, k=1,ndim)
     end do
     write(*,*) "well2"
     do j=1, natom
        write(*,*)  label(j), (well2(k,j)*0.529177d0, k=1,ndim)
     end do

     do i=1,ndim
        do j=1,natom
           do k=1,n
              xtilde(k,i,j)= (dble(n-k)*well1(i,j) + dble(k-1)*well2(i,j))/dble(n-1)
           end do
        end do
     end do
     open(20, file="aligned.xyz")
     do i=1,n
        write(20,*) natom
        write(20,*) "Energy of minimum",i
        do j=1, natom
           write(20,*)  label(j), (xtilde(i,k,j)*0.529177d0, k=1,ndim)
        end do
     end do
     close(20)

  end if

  write(*,*) "Locating instanton"
  call instanton(xtilde,well1,well2)
  write(*,*) "Found instanton."
  open(19, file="instanton.xyz")

  allocate(Vpath(n))
  do i=1, n
     Vpath(i)= V(xtilde(i,:,:))
  end do
  if (minval(Vpath) .lt. 0.0d0) then
     i= minloc(Vpath, dim=1)
     V0=0
     V0= V(xtilde(i,:,:))
     write(*,*) "New V0=", V0
     write(*,*) "New wells=", V(well1), V(well2)
  end if

  do i=1,n
     ! write(*,*) i, V(xtilde(i,:,:))
     write(19,*) natom
     write(19,*) "Energy of minimum", Vpath(i)
     do j=1, natom
        write(19,*)  label(j), (xtilde(i,k,j)*0.529177d0, k=1,ndim), (0.0d0, k=ndim+1,3)
     end do
  end do
  close(19)


  allocate(lampath(n))
  open(21, file="instantonpathenergy.csv")
  do i=1,n
     if (i.eq.1) then
        lampath(1)=0.0d0
     else
        lampath(i)= lampath(i-1) + eucliddist(xtilde(i-1,:,:), xtilde(i,:,:))
     end if
  end do
  lampath(:)= lampath(:)/lampath(n)
  do i=1,n
     write(21,*) lampath(i),Vpath(i)
  end do
  close(21)
  deallocate(lampath)
  !------------------------------------
  !work out Q_0
  allocate(xharm(n,ndim, natom), etasquared(totdof))
  do i=1,n
     xharm(i,:,:)= well1(:,:)
  end do
  etasquared(:)=0.0d0
  call detJ(xharm, etasquared,.true.)
  lndetj0= 0.0d0
  zerocount=0
  do i=1,totdof
     ! write(*,*) i,etasquared(i)
     if (etasquared(i) .gt. 0.0d0) then
        lndetj0= lndetj0+ log(etasquared(i))
     else
        zerocount=zerocount+1
     end if
  end do
  write(*,*) "Skipped ", zerocount, "states"

  !------------------------------------
  !loop over solid angle points
  if (angular) then
     allocate(theta(npoints), weightstheta(npoints))
     allocate(phi(npoints), weightsphi(npoints))
     allocate(eta(npoints), weightseta(npoints))
     allocate(xtilderot(n,ndim,natom), wellrot(ndim,natom))
     if (ndim.eq.3) then
     call gauleg(0d0,min(2.0d0*pi, cutofftheta),theta, weightstheta,npoints)
     call gauleg(0d0,min(pi, cutoffphi),phi, weightsphi,npoints)
     call gauleg(0d0,min(2.0d0*pi, cutofftheta),eta, weightseta,npoints)
     open(90, file="angularI.dat")
     xtilderot(:,:,:)= xtilde(:,:,:)
     do ii=1,npoints
        do jj=1,npoints
           do kk=1,npoints
              wellrot(:,:)= well2(:,:)
              xtilderot(:,:,:)= xtilde(:,:,:)
              call rotate_atoms(wellrot,1,eta(kk))
              call rotate_atoms(wellrot,3,phi(jj))
              call rotate_atoms(wellrot,1,theta(ii))
              call instanton(xtilderot,well1,wellrot)
              if (output_instanton) then
                 open(19, file="instanton.xyz")
                 do i=1,n
                    write(19,*) natom
                    write(19,*) "Energy of minimum",i
                    do j=1, natom
                       write(19,*)  label(j), (xtilderot(i,k,j)*0.529177d0, k=1,ndim)
                    end do
                 end do
                 close(19)
              end if

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
              write(*,*) "lndetJ", lndetj, lndetj0
              gammetilde= exp(0.5d0*(lndetJ-lndetJ0))
              !-------------------------
              !Put it all together and output.
              Skink= betan*UM(xtilderot, well1, wellrot)
              omega= betan*exp(-skink)*sqrt(skink/(2.0d0*pi))/gammetilde
              if (omega*N > 1.0d0) then
                 Ibeta=1.0
              else
                 Ibeta= tanh(omega*N)
              end if
              write(90,*) theta(ii), phi(jj), eta(kk),weightstheta(ii)*weightsphi(jj)*weightseta(kk), Ibeta
              write(*,*) theta(ii),phi(jj),eta(kk), weightstheta(ii)*weightsphi(jj)*weightseta(kk), Ibeta, lndetj
              ! end do
           end do
        end do
     end do
     else if (ndim.eq.2) then
     open(90, file="angularI.dat")
        call gauleg(0d0,min(2.0d0*pi, cutofftheta),eta, weightseta,npoints)
           do kk=1,npoints
              ! eta(kk)= pi*dble(kk)/dble(npoints)
              wellrot(:,:)= well2(:,:)
              xtilderot(:,:,:)= xtilde(:,:,:)
              call rotate_atoms(wellrot,1,eta(kk))
              call instanton(xtilderot,well1,wellrot)
              if (output_instanton) then
                 open(19, file="instanton.xyz")
                 do i=1,n
                    write(19,*) natom
                    write(19,*) "Energy of minimum",i
                    do j=1, natom
                       write(19,*)  label(j), (xtilderot(i,k,j), k=1,ndim)
                    end do
                 end do
                 close(19)
              end if

              call detJ(xtilderot, etasquared,.false.)

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
              gammetilde= exp(0.5d0*(lndetJ-lndetJ0))
              !-------------------------
              !Put it all together and output.
              Skink= betan*UM(xtilderot, well1, wellrot)
              omega= betan*exp(-skink)*sqrt(skink/(2.0d0*pi))/gammetilde
              Ibeta= tanh(omega*N)
              write(90,*) eta(kk),weightseta(kk), Ibeta
              write(*,*) eta(kk), weightseta(kk), Ibeta, lndetj, skink
              ! end do
           end do
     end if
     deallocate(xtilderot, wellrot)
     close(90)
  else
  !-------------------------
  !work out Hessian at the instantons
     !-------------------------
     !Put it all together and output.
     call detJ(xtilde, etasquared,.false.)
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

     if (fixedends) then
        Skink= betan*UM(xtilde, well1, well2)
     else
        Skink= betan*UM(xtilde)
     end if
     omega= betan*exp(-skink)*sqrt(skink/(2.0d0*pi))/gammetilde
     delta= 2.0d0*omega/betan
     write(*,*) "phi=", gammetilde
     write(*,*) "s_kink=", skink
     write(*,*) "theta, N*theta=",omega,omega*N
     write(*,*) "h_ij=", omega/betan, 219475.0d0*omega/betan
     write(*,*) "RPI delta=", delta, delta*219475.0d0
  end if

  deallocate(Vpath)

  deallocate(xharm, etasquared, xtilde)
end program rpi
