module mcmod_mass
  implicit none
  double precision::               V0, eps2=1.0d-5
  integer::                        atom1, atom2, atom3
  character(len=:), allocatable::  procdir, beaddir(:)
  character, allocatable::         label(:)
  character(len=20)::              basename
  integer::                        n, ndim, ndof, natom, xunit, totdof, potprocs
  logical::                        potforcepresent
  
contains

  ! subroutine file_init(force,i)
  !   implicit none
  !   integer, intent(in), optional:: i
  !   logical, intent(in):: force
  !   character(len=25)::  wfufilename, fulldir

  !   if (present(i)) then
  !      fulldir= trim(procdir) // "/" // trim(beaddir(i))
  !   else
  !      fulldir= trim(procdir)
  !   end if
    
  !   if (force) then
  !      call EXECUTE_COMMAND_LINE("cp " //trim(basename)// "_force.com "//trim(fulldir))
  !      open(50, file=trim(fulldir)//"/"//trim(basename)//"_force.com", status="old", position="append")
  !   else
  !      call EXECUTE_COMMAND_LINE("cp " //trim(basename)// ".com "//trim(fulldir))
  !      open(50, file=trim(fulldir)//"/"//trim(basename)//".com", status="old", position="append")
  !   end if

  !   close(50)
  !   return
  ! end subroutine file_init
  
  subroutine V_init(iproc)
    integer, intent(in):: iproc
    character(len=7)::  format_string
    character:: status
    logical::    ex
    integer::  i,readstat, ierr

    potforcepresent= .true.

    !-----------------------
    !format statement for reading in
    if (iproc .lt. 10) then
       format_string = "(A4,I1)"
       allocate(character(len=5)::procdir)
    else if (iproc .ge. 10 .and. iproc .lt. 100) then
       format_string = "(A4,I2)"
       allocate(character(len=6)::procdir)
    else if (iproc .ge. 100) then
       format_string = "(A4,I3)"
       allocate(character(len=7)::procdir)
    endif
    write(procdir,format_string) "proc", iproc
    V0=0.0d0
    !-----------------------
    !make directory structure
    inquire(file=procdir//"/geometry", exist=ex)
    if (.not. ex) then
       call EXECUTE_COMMAND_LINE("mkdir "//procdir)
       call EXECUTE_COMMAND_LINE("mkfifo "//procdir//"/geometry")
       call EXECUTE_COMMAND_LINE("mkfifo "//procdir//"/gradient")
       call EXECUTE_COMMAND_LINE("mkfifo "//procdir//"/statusmol")
       call EXECUTE_COMMAND_LINE("mkfifo "//procdir//"/statuspi")
       call EXECUTE_COMMAND_LINE("cp " //trim(basename)// ".com "//trim(procdir))
       call EXECUTE_COMMAND_LINE("cp geometry.xyz "//trim(procdir))
    end if
    ! call CHDIR(procdir)
    open(1000, file=trim(procdir)//"/geometry",form="formatted", status="old")
    open(2000, file=trim(procdir)//"/gradient",form="formatted", status="old")
    open(3000, file=trim(procdir)//"/statusmol",form="formatted", status="old")
    open(4000, file=trim(procdir)//"/statuspi",form="formatted", status="old")
    call EXECUTE_COMMAND_LINE("molpro -n 1 --no-xml-output --nouse-logfile --no-flush6 -d ./" // trim(procdir) &
         //" -s ./"//trim(procdir) //"/" // trim(basename)//".com", wait=.false.) !
    status=" "
    readstat=1
    
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x)
    implicit none
    double precision, intent(in)::     x(:,:)
    double precision::  v,throwaway
    double precision, allocatable:: dummy1(:,:)
    character(len=200)::   intext
    character(len=25)::    inword,fulldir
    integer::              i,j,k, ierr, linecount

    allocate(dummy1(ndim,natom))
    call potforce(x,dummy1,V)
    deallocate(dummy1)
    ! write(*,*)V
    return
  end function V

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    double precision, intent(in)::  x(:,:)
    double precision, intent(out):: grad(:,:)
    integer::              i,j,k,linecount
    character::            throw0
    character(len=25)::    fulldir
    double precision::     dummy1, eps, potplus, potminus
    double precision::     throw1,throw2,throw3,throw4
    double precision, allocatable:: xtemp(:,:)

    call potforce(x,grad,throw1)
    return
  end subroutine Vprime

    !---------------------------------------------------------------------
  subroutine potforce(x,grad,energy)
    implicit none
    double precision,intent(in)::     x(:,:)
    double precision, intent(out)::  grad(:,:), energy
    double precision::               throwaway
    double precision, allocatable::  dummy1(:),dummy2(:), xtemp(:)
    character(len=200)::   intext
    character(len=25)::    inword1, inword2, inword3,fulldir
    integer::              i,j,k, ierr, linecount, readstat
    character::            status
    double precision::     throw1,throw2,throw3,throw4

    ! call CHDIR(procdir)
    do i=1,natom
       do j=1, 3
          write(1000,*) x(j,i)
       end do
    end do
    flush(1000)
    write(4000,*)'C'
    flush(4000)
    readstat=1
    do
       read(3000,*,iostat=readstat) status
       if (readstat.eq.0 .and. status.eq."G") exit
    end do
    read(2000,*) energy
    energy=energy-V0
    do i=1,natom
       do j=1,3
          read(2000,*)grad(j,i)
       end do
    end do
    ! call CHDIR("..")
    

    return
  end subroutine Potforce
!---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1, eps
    integer::              i, j,k,l
    double precision:: pot1, pot2, pot3, pot4
    ! double precision, allocatable:: xtemp(:,:)
    double precision, allocatable::     gradplus(:, :), gradminus(:, :)


    allocate(gradplus(ndim,natom), gradminus(ndim,natom))
    eps=1d-2
    do i= 1, ndim
       do j= 1, natom
          x(i,j)= x(i,j) + eps
          call Vprime(x, gradplus)
          x(i,j)= x(i,j) - 2.0d0*eps
          call Vprime(x, gradminus)
          x(i,j)= x(i,j) + eps
          hess(i,j,:,:)= (gradplus(:,:)-gradminus(:,:))/(2.0d0*eps)          
       end do
    end do

    deallocate(gradplus,gradminus)
    return
  end subroutine Vdoubleprime

  subroutine V_finalize()
    call CHDIR(procdir)
    write(4000,*)'F'
    close(1000)
    close(2000)
    close(3000)
    close(4000)
    call EXECUTE_COMMAND_LINE("rm geometry gradient statusmol statuspi")
    call CHDIR("..")
    deallocate(procdir)
  end subroutine V_finalize
  
end module mcmod_mass
