module parallelmod
  use instantonmod
  implicit none
  include 'mpif.h'

contains

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !optimizes and returns the instanton
  subroutine parallel_instanton(xtilde,iproc,nproc,a,b)
    implicit none
    double precision, intent(inout)::xtilde(:,:,:)
    double precision, intent(in), optional:: a(:,:),b(:,:)
    integer, intent(in)::            iproc, nproc
    integer::                        iprint, m, iflag, mp,idof, maxiter
    integer::                        i, lp, count, iw, j,k, ierr
    double precision::               xtol, gtol, stpmin, stpmax
    double precision::               f, xtemp(ndim,natom)

    double precision::               factr, com(ndim)
    double precision, allocatable::  fprime(:,:,:), work(:), fprimework(:)
    double precision, allocatable::  lb(:), ub(:), dsave(:), xwork(:)
    integer, allocatable::           nbd(:), iwork(:), isave(:)
    logical::                        lsave(4)
    character(len=60)::              task, csave


    if (iproc .eq.0) then
       allocate(lb(totdof), ub(totdof),fprime(n,ndim,natom), nbd(totdof))
       allocate(fprimework(totdof), xwork(totdof))
       do i=1, n, 1
          do j=1,ndim
             do k=1,natom
                idof= ((k-1)*ndim + j -1)*n +i
                if (fixedends) then
                   lb(idof)= a(j,k)
                   ub(idof)= a(j,k)
                else
                   lb(idof)= 0.0d0
                   ub(idof)= 0.0d0
                end if
                nbd(idof)=0
             end do
          end do
       end do

       !------------------------
       !perform minimization
       m=8
       iprint=-1
       xtol= 1d-5
       iw=totdof*(2*m+5) + 11*m**2 + 8*m
       allocate(work(iw), iwork(3*totdof), isave(44), dsave(29))
       iflag=0
       ! eps2= 1.0d-5 !gradient convergence
       factr=1.0d6
       maxiter=40
    end if
    task='START'

    if (fixedends) then
       if (potforcepresent) then
          call parallel_UMforceenergy(xtilde,iproc, nproc, fprime,f,a,b)
       else
          f= parallel_UM(xtilde,iproc, nproc, a,b)
          call parallel_UMprime(xtilde,iproc, nproc, fprime,a,b)
       end if
    else
       if (potforcepresent) then
          call parallel_UMforceenergy(xtilde,iproc, nproc, fprime,f)
       else
          write(*,*) "iproc ", iproc, "starting potential evaluation: l73"
          f= parallel_UM(xtilde,iproc, nproc)
          call parallel_UMprime(xtilde,iproc, nproc, fprime)
       end if
    end if

    count=0
    call MPI_Bcast(task, 5, MPI_CHARACTER, 0,MPI_COMM_WORLD, ierr)
    do while( task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
         task.eq.'START')
       if (iproc .eq. 0) then
          count=count+1
          write(*,*) "Iteration ", count, f
          xwork=reshape(xtilde,(/totdof/))
          fprimework= reshape(fprime,(/totdof/))
          call setulb(totdof,m,xwork,lb,ub,nbd,f,fprimework,factr,eps2,work&
               ,iwork,task,iprint, csave,lsave,isave,dsave,maxiter)
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
          call MPI_Bcast(task, 5, MPI_CHARACTER, 0,MPI_COMM_WORLD, ierr)
       end if
       if (task(1:2) .eq. 'FG') then
          if (iproc .eq. 0) then
             write(*,*) "iteration", count
             xtilde= reshape(xwork,(/n,ndim,natom/))
          end if
          if (fixedends) then
             if (potforcepresent) then
                call parallel_UMforceenergy(xtilde,iproc, nproc, fprime,f,a,b)
             else
                f= parallel_UM(xtilde,iproc, nproc, a,b)
                call parallel_UMprime(xtilde,iproc, nproc, fprime,a,b)
             end if
          else
             if (potforcepresent) then
                call parallel_UMforceenergy(xtilde,iproc, nproc, fprime,f)
             else
                f= parallel_UM(xtilde,iproc, nproc)
                call parallel_UMprime(xtilde,iproc, nproc, fprime)
             end if
          end if
       end if
    end do
    if (task(1:5) .eq. "ERROR" .or. task(1:4) .eq. "ABNO") then
       write(*,*) "Error:"
       write(*,*) task
    end if

    if (iproc .eq. 0) then
       deallocate(work, lb, ub, fprime, fprimework,xwork)
       deallocate(iwork, nbd, isave, dsave)
    end if
    return
  end subroutine parallel_instanton
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !returns the fluctuation prefactor
  subroutine parallel_detJ(x, iproc, nproc,etasquared)
    double precision,intent(in)::    x(:,:,:)
    double precision,intent(inout)::  etasquared(:)
    integer, intent(in)::            iproc, nproc
    character::                      jobz, range, uplo
    double precision::               vl, vu, abstol
    integer::                        nout, ldz, lwork, liwork, info,i
    integer,allocatable::            isuppz(:), iwork(:)
    double precision, allocatable::  work(:), z(:,:), H(:,:)


    if (iproc .eq. 0) then
       jobz='N'
       lwork= 2*totdof
       liwork= 1
       range='A'
       uplo='U'
       abstol=1.0d-8
       info=0
       ldz=totdof
       vl=0.0d0
       vu=0.0d0
       nout=0
       allocate(work(lwork), iwork(liwork), H(ndof+1,totdof))
       H=0.0d0
       etasquared=0.0d0
    end if

    call parallel_UMhessian(x,iproc, nproc, H)
    write(*,*) "Hessian is cooked."

    if (iproc .eq. 0) then
       call DSBEVD(jobz, 'L', totdof, ndof, H, ndof+1, etasquared, z, 1, work, lwork, iwork,&
            liwork, info)
       deallocate(work, iwork, H)
    end if
    return
  end subroutine parallel_detJ


  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !linear polymer potential
  function parallel_UM(x,iproc, nproc,a,b)
    implicit none
    double precision,intent(in)::   x(:,:,:)
    double precision,intent(in),optional:: a(:,:),b(:,:)
    integer, intent(in)::          iproc, nproc
    double precision, allocatable:: xpart(:,:,:), Vpart(:), Vall(:)
    double precision:: parallel_UM, pot
    integer::            i,j,k, ncalcs, ierr, startind, ncalcproc
    integer, dimension(MPI_STATUS_SIZE) :: rstatus

    parallel_UM=0.0d0

    !Begin Parallel parts!
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ncalcs= N/nproc
    if (iproc .lt. mod(N, nproc)) ncalcs=ncalcs+1
    write(*,*) "iproc ", iproc, "starting potential evaluation with ncalcs=", ncalcs
    allocate(xpart(ncalcs,ndim,natom),Vpart(ncalcs))
    if (iproc .eq. 0) then
       !need to send x to all the procs
       startind=1+ncalcs
       do i=1,nproc-1
          ncalcproc= N/nproc
          if (i .lt. mod(N, nproc)) ncalcproc=ncalcproc+1
          write(*,*) "iproc ", iproc, "sending to ", i, ncalcproc, "slices of x"
          write(*,*) "starting from ", startind
          do j=1, ncalcproc
             call MPI_Send(x(startind+j-1,:,:), ndof, MPI_DOUBLE_PRECISION, i, 1,&
                  MPI_COMM_WORLD, ierr)
          end do
          startind= startind+ ncalcproc
       end do
       xpart(1:ncalcs,:,:) = x(1:ncalcs,:,:)
    else
       !need to receive x to all procs
       write(*,*) "iproc ", iproc, "receiving ", ncalcs, "slices of x"
       do i=1,ncalcs
          call MPI_Recv(xpart(i,:,:),ndof, MPI_DOUBLE_PRECISION, 0, 1,&
               MPI_COMM_WORLD, rstatus, ierr)
       end do
    end if
    write(*,*) "iproc ", iproc, "received the x: l204"

    do i=1, ncalcs
       !need to calculate the potential
       Vpart(i)= V(xpart(i,:,:))
    end do
    write(*,*) "iproc ", iproc, "finished calculating"
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !gather all the results
    allocate(Vall(n))
    call MPI_Gather(Vpart,ncalcs,MPI_DOUBLE_PRECISION, Vall, ncalcs, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD, ierr)
    write(*,*) "iproc ", iproc, "finished gathering"
    deallocate(xpart,Vpart)

    !Do easy bit
    if (iproc .eq. 0) then
       do i=1, N-1, 1
          parallel_UM=parallel_UM + Vall(i)
          do j=1, ndim
             do k=1, natom
                parallel_UM=parallel_UM+ (0.5d0*mass(k)/betan**2)*(x(i+1,j,k)-x(i,j,k))**2
             end do
          end do
       end do
       parallel_UM=parallel_UM+ Vall(n)
       if (fixedends) then
          do j=1, ndim
             do k=1, natom
                parallel_UM=parallel_UM+ (0.5d0*mass(k)/betan**2)*(x(1,j,k)-a(j,k))**2
                parallel_UM=parallel_UM+ (0.5d0*mass(k)/betan**2)*(b(j,k)-x(N,j,k))**2
             end do
          end do
       end if
    end if
    deallocate(Vall)
    write(*,*) iproc, "reached the final barrier"
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    return
  end function PARALLEL_UM

  !---------------------------------------------------------------------
  !linear polymer force
  subroutine parallel_UMprime(x,iproc,nproc, answer,a,b)
    implicit none
    double precision, intent(in)::   x(:,:,:)
    double precision, intent(out)::  answer(:,:,:)
    double precision, intent(in), optional:: a(:,:),b(:,:)
    integer, intent(in)::          iproc, nproc
    double precision,allocatable:: gradpart(:,:,:), xpart(:,:,:), gradall(:,:,:)
    integer::            i,j,k, ncalcs, ierr, startind, ncalcproc
    integer, dimension(MPI_STATUS_SIZE) :: rstatus


    !Begin Parallel parts!
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ncalcs= N/nproc
    if (iproc .lt. mod(N, nproc)) ncalcs=ncalcs+1
    allocate(xpart(ncalcs,ndim,natom),gradpart(ncalcs,ndim,natom))
    if (iproc .eq. 0) then
       !need to send x to all the procs
       startind=1+ncalcs
       do i=1,nproc-1
          ncalcproc= N/nproc
          if (i .lt. mod(N, nproc)) ncalcproc=ncalcproc+1
          do j=1, ncalcs
             call MPI_Send(x(startind+j-1,:,:), ndof, MPI_DOUBLE_PRECISION, i, 1, MPI_COMM_WORLD, ierr)
          end do
          startind= startind+ ncalcproc
       end do
       xpart(1:ncalcs,:,:) = x(1:ncalcs,:,:)
    else
       !need to receive x to all procs
       do i=1,ncalcs
          call MPI_Recv(xpart(i,:,:),ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus, ierr)
       end do
    end if
    do i=1, ncalcs
       !need to calculate the potential
       call Vprime(xpart(i,:,:),gradpart(i,:,:))
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !gather all the results
    allocate(gradall(n,ndim,natom))
    call MPI_Gather(gradpart,ncalcs*ndim*natom,MPI_DOUBLE_PRECISION, gradall, ncalcs*ndim*natom,&
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    deallocate(xpart,gradpart)

    do i=1, N
       do j=1,ndim
          do k=1,natom
             if (i.eq.1) then
                if (fixedends) then
                   answer(1,j,k)=mass(k)*(2.0*x(1,j,k) - a(j,k) - x(2,j,k))/betan**2
                else
                   answer(1,j,k)=mass(k)*(x(1,j,k) - x(2,j,k))/betan**2                   
                end if
             else if (i.eq.N) then
                if (fixedends) then
                   answer(N,j,k)=mass(k)*(2.0*x(N,j,k) - x(N-1,j,k) - b(j,k))/betan**2
                else
                   answer(N,j,k)=mass(k)*(x(N,j,k) - x(N-1,j,k))/betan**2
                end if
             else
                answer(i,j,k)= mass(k)*(2.0*x(i,j,k) - x(i-1,j,k) - x(i+1,j,k))/betan**2
             end if
          end do
       end do
       do j=1, ndim
          do k=1, natom
             answer(i,j,k)= answer(i,j,k)+ gradall(i,j,k)
          end do
       end do
    end do
    deallocate(gradall)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    return
  end subroutine Parallel_UMprime
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !linear polymer force and energy
  subroutine parallel_UMforceenergy(x, iproc, nproc,answer,energy,a,b)
    implicit none
    double precision,intent(in)::   x(:,:,:)
    double precision,intent(in),optional:: a(:,:),b(:,:)
    double precision, intent(out)::  answer(:,:,:), energy
    integer, intent(in)::          iproc, nproc
    double precision, allocatable:: xpart(:,:,:), Vpart(:), Vall(:)
    double precision,allocatable:: gradpart(:,:,:), gradall(:,:,:)
    integer::            i,j,k, ncalcs, ierr, startind, ncalcproc
    integer, dimension(MPI_STATUS_SIZE) :: rstatus

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    energy=0.0d0
    !Begin Parallel parts!
    ncalcs= N/nproc
    if (iproc .lt. mod(N, nproc)) ncalcs=ncalcs+1
    allocate(xpart(ncalcs,ndim,natom),Vpart(ncalcs))
    allocate(gradpart(ncalcs,ndim,natom))
    if (iproc .eq. 0) then
       !need to send x to all the procs
       startind=1+ncalcs
       do i=1,nproc-1
          ncalcproc= N/nproc
          if (i .lt. mod(N, nproc)) ncalcproc=ncalcproc+1
          do j=1, ncalcs
             call MPI_Send(x(startind+j-1,:,:), ndof, MPI_DOUBLE_PRECISION, i, 1, MPI_COMM_WORLD, ierr)
          end do
          startind= startind+ ncalcproc
       end do
       xpart(1:ncalcs,:,:) = x(1:ncalcs,:,:)
    else
       !need to receive x to all procs
       do i=1,ncalcs
          call MPI_Recv(xpart(i,:,:),ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus, ierr)
       end do
    end if
    do i=1, ncalcs
       !need to calculate the potential
       call potforce(x(i,:,:),gradpart(i,:,:),Vpart(i))
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !gather all the results
    allocate(Vall(n),gradall(n,ndim,natom))
    call MPI_Gather(Vpart,ncalcs,MPI_DOUBLE_PRECISION, Vall, ncalcs, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD, ierr)
    call MPI_Gather(gradpart,ncalcs*ndim*natom,MPI_DOUBLE_PRECISION, gradall, ncalcs*ndim*natom,&
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    deallocate(xpart,Vpart, gradpart)

    do i=1, N, 1
       energy=energy+ Vall(i)
       do j=1, ndim
          do k=1, natom
             if (i .lt. N) energy=energy+ (0.5d0*mass(k)/betan**2)*(x(i+1,j,k)-x(i,j,k))**2
          end do
       end do
       do j=1,ndim
          do k=1,natom
             if (i.eq.1) then
                if (fixedends) then
                   answer(1,j,k)=mass(k)*(2.0*x(1,j,k) - a(j,k) - x(2,j,k))/betan**2
                else
                   answer(1,j,k)=mass(k)*(x(1,j,k) - x(2,j,k))/betan**2                   
                end if
             else if (i.eq.N) then
                if (fixedends) then
                   answer(N,j,k)=mass(k)*(2.0*x(N,j,k) - x(N-1,j,k) - b(j,k))/betan**2
                else
                   answer(N,j,k)=mass(k)*(x(N,j,k) - x(N-1,j,k))/betan**2
                end if
             else
                answer(i,j,k)= mass(k)*(2.0*x(i,j,k) - x(i-1,j,k) - x(i+1,j,k))/betan**2
             end if
          end do
       end do
       do j=1, ndim
          do k=1, natom
             answer(i,j,k)= answer(i,j,k)+ gradall(i,j,k)
          end do
       end do
    end do

    if (fixedends) then
       ! energy=energy + V(a(:,:))+ V(b(:,:))
       do j=1, ndim
          do k=1, natom
             energy=energy+ (0.5d0*mass(k)/betan**2)*(x(1,j,k)-a(j,k))**2
             energy=energy+ (0.5d0*mass(k)/betan**2)*(b(j,k)-x(N,j,k))**2
          end do
       end do
    end if

    deallocate(gradall, Vall)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    return
  end subroutine Parallel_UMforceenergy

  !---------------------------------------------------------------------
  !linear polymer hessian
  subroutine parallel_UMhessian(x, iproc, nproc,answer)
    implicit none
    double precision, intent(in)::   x(:,:,:)
    integer, intent(in)::          iproc, nproc
    double precision, intent(out):: answer(:,:)
    double precision, allocatable:: hess(:,:,:,:)
    double precision,allocatable:: hesspart(:,:,:,:,:), xpart(:,:,:), hessall(:,:,:,:,:)
    integer::            i, j1, k1, j2, k2, idof1, idof2, startind, ncalcproc
    integer::            fulldof1, fulldof2, index,j,k, ncalcs, ierr
    integer, dimension(MPI_STATUS_SIZE) :: rstatus

    answer=0.0d0
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !Begin Parallel parts!
    ncalcs= N/nproc
    if (iproc .lt. mod(N, nproc)) ncalcs=ncalcs+1
    allocate(xpart(ncalcs,ndim,natom),hesspart(ncalcs,ndim,natom,ndim,natom))
    if (iproc .eq. 0) then
       !need to send x to all the procs
       startind=1+ncalcs
       do i=1,nproc-1
          ncalcproc= N/nproc
          if (i .lt. mod(N, nproc)) ncalcproc=ncalcproc+1
          do j=1, ncalcs
             call MPI_Send(x(startind+j-1,:,:), ndof, MPI_DOUBLE_PRECISION, i, 1, MPI_COMM_WORLD, ierr)
          end do
          startind= startind+ ncalcproc
       end do
       xpart(1:ncalcs,:,:) = x(1:ncalcs,:,:)
    else
       !need to receive x to all procs
       do i=1,ncalcs
          call MPI_Recv(xpart(i,:,:),ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus, ierr)
       end do
    end if
    do i=1, ncalcs
       !need to calculate the potential
       call Vdoubleprime(xpart(i,:,:), hesspart(i,:,:,:,:))
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !gather all the results
    allocate(hessall(n,ndim,natom,ndim,natom))
    call MPI_Gather(hesspart,ncalcs*ndof*ndof,MPI_DOUBLE_PRECISION, hessall, ncalcs*ndof*ndof,&
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    deallocate(xpart,hesspart)

    if (iproc .eq. 0) then
       do i=1, n, 1
          do j1=1,ndim
             do k1=1,natom
                do j2=1,ndim
                   do k2=1,natom
                      !The DoF label for the individual bead
                      idof1= (k1-1)*ndim + j1
                      idof2= (k2-1)*ndim + j2
                      !The DoF label for the whole matrix
                      fulldof1=ndof*(i-1) + idof1
                      fulldof2=ndof*(i-1) + idof2

                      if (fulldof2 .lt. fulldof1) cycle

                      if (idof1.eq.idof2) then
                         answer(1,fulldof1)= 2.0d0/betan**2&
                              +hessall(i,j2,k2,j1,k1)/sqrt(mass(k1)*mass(k2)) 
                         if (i.gt.1) answer(ndof+1,fulldof1)=-1.0d0/betan**2
                      else
                         index=1 + fulldof2 -fulldof1
                         if (index.lt.0) cycle
                         answer(index,fulldof1)= &
                              +hessall(i,j2,k2,j1,k1)/sqrt(mass(k1)*mass(k2)) 
                      end if
                   end do
                end do
             end do
          end do
       end do
    end if
    deallocate(hessall)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    return
  end subroutine Parallel_UMhessian

end module parallelmod
