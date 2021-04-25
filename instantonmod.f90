module instantonmod
  use mcmod_mass
  implicit none
  double precision, parameter::    pi=3.14159265358979d0
  double precision::               beta, betan, UMtilde
  double precision, allocatable::  xtilde(:,:,:),well1(:,:), well2(:,:), mass(:)
  logical::                        fixedends

contains


  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !linear polymer potential
  function UM(x,a,b)
    implicit none
    integer::            i,j,k
    double precision,intent(in)::   x(:,:,:)
    double precision,intent(in),optional:: a(:,:),b(:,:)
    double precision:: UM, pot

    UM=0.0d0
    do i=1, N-1, 1
       pot= V(x(i,:,:))
       UM=UM + pot
       do j=1, ndim
          do k=1, natom
             UM=UM+ (0.5d0*mass(k)/betan**2)*(x(i+1,j,k)-x(i,j,k))**2
          end do
       end do
    end do
    UM=UM+ V(x(N,:,:))
    if (fixedends) then
       !UM=UM + V(a(:,:))+ V(b(:,:)) !already set these to 0!
       do j=1, ndim
          do k=1, natom
             UM=UM+ (0.5d0*mass(k)/betan**2)*(x(1,j,k)-a(j,k))**2
             UM=UM+ (0.5d0*mass(k)/betan**2)*(b(j,k)-x(N,j,k))**2
          end do
       end do
    end if

    return
  end function UM

  !---------------------------------------------------------------------
  !linear polymer force
  subroutine UMprime(x, answer,a,b)
    implicit none
    integer::            i,j,k
    double precision, intent(in)::   x(:,:,:)
    double precision, intent(out)::  answer(:,:,:)
    double precision, intent(in), optional:: a(:,:),b(:,:)
    double precision,allocatable:: grad(:,:)

    allocate(grad(ndim,natom))
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
  end subroutine UMprime
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !linear polymer force and energy
  subroutine UMforceenergy(x, answer,UM,a,b)
    implicit none
    integer::            i,j,k
    double precision,intent(in)::   x(:,:,:)
    double precision,intent(in),optional:: a(:,:),b(:,:)
    double precision, intent(out)::  answer(:,:,:), UM
    double precision,allocatable:: grad(:,:)
    double precision:: energy

    allocate(grad(ndim,natom))
    UM=0.0d0
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
  end subroutine UMforceenergy

  !---------------------------------------------------------------------
  !linear polymer hessian
  subroutine UMhessian(x, singlewell,answer,inithess)
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
  end subroutine UMhessian


  !---------------------------------------------------------------------
!---------------------------------------------------------------------
!optimizes and returns the instanton
  subroutine instanton(xtilde,a,b)
    implicit none
    double precision, intent(inout)::xtilde(:,:,:)
    double precision, intent(in), optional:: a(:,:),b(:,:)
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
    if (fixedends) then
       if (potforcepresent) then
          call UMforceenergy(xtilde,fprime,f,a,b)
       else
          f= UM(xtilde,a,b)
          call UMprime(xtilde,fprime,a,b)
       end if
    else
       if (potforcepresent) then
          call UMforceenergy(xtilde,fprime,f)
       else
          f= UM(xtilde)
          call UMprime(xtilde,fprime)
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
                call UMforceenergy(xtilde,fprime,f,a,b)
             else
                f= UM(xtilde,a,b)
                call UMprime(xtilde,fprime,a,b)
             end if
          else
             if (potforcepresent) then
                call UMforceenergy(xtilde,fprime,f)
             else
                f= UM(xtilde)
                call UMprime(xtilde,fprime)
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
  end subroutine instanton

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !returns the fluctuation prefactor
  subroutine detJ(x, etasquared, singlewell,interphess,eigvecs)
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol
  double precision, intent(in),optional::  interphess(:,:,:)
  double precision, intent(out), optional:: eigvecs(:,:)
  double precision,intent(in)::    x(:,:,:)
  logical, intent(in)::             singlewell
  double precision,intent(inout)::  etasquared(:)
  integer::                        nout, ldz, lwork, liwork, info,i
  integer,allocatable::            isuppz(:), iwork(:)
  double precision, allocatable::  work(:), z(:,:), H(:,:)
  ! !get diagonal hessian
  if (present(eigvecs)) then
     jobz='V'
     allocate(z(totdof,totdof))
     lwork= 1 + 5*totdof + 2*totdof**2
     liwork= 3 + 5*totdof
  else
     jobz='N'
     lwork= 2*totdof
     liwork= 1
  end if
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
  call UMhessian(x,singlewell,H)
  write(*,*) "Hessian is cooked."

  if (present(eigvecs)) then
     call DSBEVD(jobz, 'L', totdof, ndof, H, ndof+1, etasquared, z, totdof, work, lwork, iwork, liwork, info)
     eigvecs=z
     deallocate(z)
  else
     call DSBEVD(jobz, 'L', totdof, ndof, H, ndof+1, etasquared, z, 1, work, lwork, iwork, liwork, info)
  end if
  deallocate(work, iwork, H)
  return
  end subroutine detJ

end module instantonmod
