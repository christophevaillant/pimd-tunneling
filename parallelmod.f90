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
    integer::                        i, lp, count, iw, j,k
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
       task='START'
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
          f= parallel_UM(xtildeiproc, nproc)
          call UMprime(xtilde,iproc, nproc, fprime)
       end if
    end if
    count=0
    do while( task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
         task.eq.'START')
       count=count+1
       xwork=reshape(xtilde,(/totdof/))
       fprimework= reshape(fprime,(/totdof/))
       call setulb(totdof,m,xwork,lb,ub,nbd,f,fprimework,factr,eps2,work&
            ,iwork,task,iprint, csave,lsave,isave,dsave,maxiter)
       if (task(1:2) .eq. 'FG') then
          write(*,*) "iteration", count
          xtilde= reshape(xwork,(/n,ndim,natom/))
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

    deallocate(work, lb, ub, fprime, fprimework,xwork)
    deallocate(iwork, nbd, isave, dsave)

    return
  end subroutine parallel_instanton

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !linear polymer potential
  function parallel_UM(x,iproc, nproc,a,b)
    implicit none
    double precision,intent(in)::   x(:,:,:)
    double precision,intent(in),optional:: a(:,:),b(:,:)
    integer, intent(in)::          iproc, nproc
    double precision, allocatable:: xpart(:,:,:), Vpart(:), Vall(:)
    double precision:: UM, pot
    integer::            i,j,k, ncalcs, ierr
    integer, dimension(MPI_STATUS_SIZE) :: rstatus

    UM=0.0d0

    !Begin Parallel parts!
    ncalcs= N/nproc
    if (iproc .lt. mod(N, nproc)) ncalcs=ncalcs+1
    allocate(xpart(ncalcs,ndim,natom),Vpart(ncalcs))
    if (iproc .eq. 0) then
       !need to send x to all the procs
       do i=1,nproc-1
          do j=1, ncalcs
             call MPI_Send(x(ncalcs*i +j,:,:), ndof, MPI_DOUBLE_PRECISION, j, 1, MPI_COMM_WORLD, ierr)
          end do
       end do
       xpart(1:ncalcs,:,:) = x(1:ncalcs,:,:)
    else
       !need to receive x to all procs
       do i=1,ncalcs
          call MPI_Recv(xpart(i,:,:),ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus, ierr))
       end do
    end if
    do i=1, ncalcs
       !need to calculate the potential
       Vpart(i)= V(xpart(i,:,:))
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !gather all the results
    allocate(Vall(n))
    call MPI_Gather(Vpart,nalcs,MPI_DOUBLE_PRECISION, Vall, ncalcs, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD, ierr)
    deallocate(xpart,Vpart)

    !Do easy bit
    if (iproc .eq. 0) then
       do i=1, N-1, 1
          UM=UM + Vall(i)
          do j=1, ndim
             do k=1, natom
                UM=UM+ (0.5d0*mass(k)/betan**2)*(x(i+1,j,k)-x(i,j,k))**2
             end do
          end do
       end do
       UM=UM+ V(x(N,:,:))
       if (fixedends) then
          do j=1, ndim
             do k=1, natom
                UM=UM+ (0.5d0*mass(k)/betan**2)*(x(1,j,k)-a(j,k))**2
                UM=UM+ (0.5d0*mass(k)/betan**2)*(b(j,k)-x(N,j,k))**2
             end do
          end do
       end if
    end if

    return
  end function PARALLEL_UM

  !---------------------------------------------------------------------
  !linear polymer force
  subroutine parallel_UMprime(x, answer,a,b)
    implicit none
    integer::            i,j,k
    double precision, intent(in)::   x(:,:,:)
    double precision, intent(out)::  answer(:,:,:)
    double precision, intent(in), optional:: a(:,:),b(:,:)
    double precision,allocatable:: grad(:,:)

    allocate(grad(ndim,natom))

    !Begin Parallel parts!
    ncalcs= N/nproc
    if (iproc .lt. mod(N, nproc)) ncalcs=ncalcs+1
    allocate(xpart(ncalcs,ndim,natom),Vpart(ncalcs))
    if (iproc .eq. 0) then
       !need to send x to all the procs
       do i=1,nproc-1
          do j=1, ncalcs
             call MPI_Send(x(ncalcs*i +j,:,:), ndof, MPI_DOUBLE_PRECISION, j, 1, MPI_COMM_WORLD, ierr)
          end do
       end do
       xpart(1:ncalcs,:,:) = x(1:ncalcs,:,:)
    else
       !need to receive x to all procs
       do i=1,ncalcs
          call MPI_Recv(xpart(i,:,:),ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus, ierr))
       end do
    end if
    do i=1, ncalcs
       !need to calculate the potential
       Vpart(i)= V(xpart(i,:,:))
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !gather all the results
    allocate(Vall(n))
    call MPI_Gather(Vpart,nalcs,MPI_DOUBLE_PRECISION, Vall, ncalcs, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD, ierr)
    deallocate(xpart,Vpart)

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
       call Vprime(x(i,:,:),grad(:,:))
       do j=1, ndim
          do k=1, natom
             answer(i,j,k)= answer(i,j,k)+ grad(j,k)
          end do
       end do
    end do
    deallocate(grad)

    return
  end subroutine Parallel_UMprime
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !linear polymer force and energy
  subroutine parallel_UMforceenergy(x, answer,UM,a,b)
    implicit none
    integer::            i,j,k
    double precision,intent(in)::   x(:,:,:)
    double precision,intent(in),optional:: a(:,:),b(:,:)
    double precision, intent(out)::  answer(:,:,:), UM
    double precision,allocatable:: grad(:,:)
    double precision:: energy

    allocate(grad(ndim,natom))
    UM=0.0d0
    !Begin Parallel parts!
    ncalcs= N/nproc
    if (iproc .lt. mod(N, nproc)) ncalcs=ncalcs+1
    allocate(xpart(ncalcs,ndim,natom),Vpart(ncalcs))
    if (iproc .eq. 0) then
       !need to send x to all the procs
       do i=1,nproc-1
          do j=1, ncalcs
             call MPI_Send(x(ncalcs*i +j,:,:), ndof, MPI_DOUBLE_PRECISION, j, 1, MPI_COMM_WORLD, ierr)
          end do
       end do
       xpart(1:ncalcs,:,:) = x(1:ncalcs,:,:)
    else
       !need to receive x to all procs
       do i=1,ncalcs
          call MPI_Recv(xpart(i,:,:),ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus, ierr))
       end do
    end if
    do i=1, ncalcs
       !need to calculate the potential
       Vpart(i)= V(xpart(i,:,:))
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !gather all the results
    allocate(Vall(n))
    call MPI_Gather(Vpart,nalcs,MPI_DOUBLE_PRECISION, Vall, ncalcs, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD, ierr)
    deallocate(xpart,Vpart)

    do i=1, N, 1
       call potforce(x(i,:,:),grad,energy)
       UM=UM+ energy
       do j=1, ndim
          do k=1, natom
             if (i .lt. N) UM=UM+ (0.5d0*mass(k)/betan**2)*(x(i+1,j,k)-x(i,j,k))**2
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
             answer(i,j,k)= answer(i,j,k)+ grad(j,k)
          end do
       end do
    end do

    if (fixedends) then
       ! UM=UM + V(a(:,:))+ V(b(:,:))
       do j=1, ndim
          do k=1, natom
             UM=UM+ (0.5d0*mass(k)/betan**2)*(x(1,j,k)-a(j,k))**2
             UM=UM+ (0.5d0*mass(k)/betan**2)*(b(j,k)-x(N,j,k))**2
          end do
       end do
    end if

    deallocate(grad)

    return
  end subroutine Parallel_UMforceenergy

  !---------------------------------------------------------------------
  !linear polymer hessian
  subroutine parallel_UMhessian(x, singlewell,answer,inithess)
    implicit none
    integer::            i, j1, k1, j2, k2, idof1, idof2
    integer::            fulldof1, fulldof2, index
    double precision, intent(in)::   x(:,:,:)
    double precision, intent(in),optional::  inithess(:,:,:)
    logical, intent(in)::   singlewell
    double precision, intent(out):: answer(:,:)
    double precision, allocatable:: hess(:,:,:,:)

    allocate(hess(ndim, natom, ndim, natom))

    answer=0.0d0
    hess=0.0d0

    !Begin Parallel parts!
    ncalcs= N/nproc
    if (iproc .lt. mod(N, nproc)) ncalcs=ncalcs+1
    allocate(xpart(ncalcs,ndim,natom),Vpart(ncalcs))
    if (iproc .eq. 0) then
       !need to send x to all the procs
       do i=1,nproc-1
          do j=1, ncalcs
             call MPI_Send(x(ncalcs*i +j,:,:), ndof, MPI_DOUBLE_PRECISION, j, 1, MPI_COMM_WORLD, ierr)
          end do
       end do
       xpart(1:ncalcs,:,:) = x(1:ncalcs,:,:)
    else
       !need to receive x to all procs
       do i=1,ncalcs
          call MPI_Recv(xpart(i,:,:),ndof, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, rstatus, ierr))
       end do
    end if
    do i=1, ncalcs
       !need to calculate the potential
       Vpart(i)= V(xpart(i,:,:))
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !gather all the results
    allocate(Vall(n))
    call MPI_Gather(Vpart,nalcs,MPI_DOUBLE_PRECISION, Vall, ncalcs, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD, ierr)
    deallocate(xpart,Vpart)

    do i=1, n, 1
       if (present(inithess)) then
          do j1=1,ndim
             do k1=1,natom
                idof1= ndim*(k1-1) + j1
                do j2=1,j1
                   do k2=1,k1
                      idof2= ndim*(k2-1) + j2
                      hess(j1,k1,j2,k2)= inithess(i,idof1,idof2)
                   end do
                end do
             end do
          end do
       else
          if ((i .eq. 1 .and. singlewell) .or. .not. singlewell) then
             call Vdoubleprime(x(i,:,:), hess)
          end if
       end if
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
                           +hess(j2,k2,j1,k1)/sqrt(mass(k1)*mass(k2)) 
                      if (i.gt.1) answer(ndof+1,fulldof1)=-1.0d0/betan**2
                   else
                      index=1 + fulldof2 -fulldof1
                      if (index.lt.0) cycle
                      answer(index,fulldof1)= &
                        +hess(j2,k2,j1,k1)/sqrt(mass(k1)*mass(k2)) 
                   end if
                end do
             end do
          end do
       end do
    end do
    deallocate(hess)
    return
  end subroutine Parallel_UMhessian

end module parallelmod
