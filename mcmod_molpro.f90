module mcmod_mass
  implicit none
  double precision::               V0, eps2=1.0d-3
  integer::                        atom1, atom2, atom3
  character(len=:), allocatable::  procdir
  character, allocatable::         label(:)
  character(len=20)::              basename
  character(len=15)::              informat
  character(len=2)::               procstring
  integer::                        n, ndim, ndof, natom, xunit, totdof, potprocs
  logical::                        potforcepresent
  
contains

  subroutine V_init(iproc)
    integer, intent(in):: iproc
    character(len=7)::  format_string
    logical::    ex
    namelist /POTDATA/ basename, atom1,atom2,atom3, informat, potprocs

    atom1=1
    atom2=2
    atom3=3
    potforcepresent= .true.
    informat="(A22,3X,F18.12)"
    potprocs=1
    
    read(5, nml=POTDATA)

    if (iproc < 10) then
       format_string = "(A4,I1)"
       allocate(character(len=5)::procdir)
    else
       format_string = "(A4,I2)"
       allocate(character(len=6)::procdir)
    endif
    write(procstring, "(I2)") potprocs
    write(procdir,format_string) "proc", iproc
    V0=0.0d0
    inquire(file=procdir, exist=ex)
    if (.not. ex) then
       call SYSTEM("mkdir "//procdir)
    end if
    call SYSTEM("cp " //trim(basename)// ".com "//procdir)
    call SYSTEM("cp " //trim(basename)// "_force.com "//procdir)
    
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x)
    implicit none
    double precision::     v, x(:,:),throwaway
    double precision, allocatable:: dummy1(:),dummy2(:), xtemp(:)
    character(len=200)::   intext
    character(len=20)::    inword
    integer::              i,j,k, ierr, linecount

    !--------------------------------
    !Output geometry to geometry file
    call CHDIR(procdir)
    open(1000, file="geometry.xyz")
    write(1000,*) natom
    write(1000,*) "Geometry of current point"
    do j=1, natom
       write(1000,*)  label(j), (x(k,j), k=1,ndim)
    end do
    close(1000)
    !--------------------------------
    !Run Molpro
    call SYSTEM("molpro -s -n "//procstring//" "//trim(basename)//".com")
    !--------------------------------
    !Search through output and find the energy
    open(unit=2000, file=trim(basename)//".out", status="OLD", access="SEQUENTIAL")
    linecount=0
    do
       read(2000,*,END=10)
       linecount=linecount+1
    end do
10  rewind(2000)

    ! do i=1,linecount-3
    !    read(2000,*)
    ! end do
    ! read(2000,*) V, throwaway
    ! close(2000)
    V=0.0d0
    do i=1,linecount-6
       read(2000,*)
    end do
    read(2000,informat) intext,V
    close(2000)

    V=V-V0
    write(*,*)V
    call CHDIR("..")
    
    return
  end function V

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j,k,linecount
    character::            throw0
    double precision::     grad(:,:), x(:,:), dummy1, eps, potplus, potminus
    double precision::     throw1,throw2,throw3,throw4

    !--------------------------------
    !Output geometry to geometry file
    call CHDIR(procdir)
    open(1000, file="geometry.xyz")
    write(1000,*) natom
    write(1000,*) "Geometry of current point"
    do j=1, natom
       write(1000,*)  label(j), (x(k,j), k=1,ndim)
    end do
    close(1000)
    !--------------------------------
    !Run Molpro
    call SYSTEM("molpro -s -n "//procstring//" "//trim(basename)//"_force.com")
    !--------------------------------
    !Search through output and find the energy
    open(unit=3000, file="forces.xyz", status="OLD", access="SEQUENTIAL")
    read(3000,*)
    read(3000,*)
    do i=1,natom
       read(3000,*) throw0,throw1,throw2,throw3,throw4, grad(1,i), grad(2,i), grad(3,i)
       grad(1,i)=grad(1,i)*0.529177d0*3.6749d-2
       grad(2,i)=grad(2,i)*0.529177d0*3.6749d-2
       grad(3,i)=grad(3,i)*0.529177d0*3.6749d-2
    end do
    close(3000)
    call CHDIR("..")
    
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
    character(len=20)::    inword1, inword2, inword3
    integer::              i,j,k, ierr, linecount
    character::            throw0
    double precision::     throw1,throw2,throw3,throw4

    !--------------------------------
    !Output geometry to geometry file
    call CHDIR(procdir)
    open(1000, file="geometry.xyz")
    write(1000,*) natom
    write(1000,*) "Geometry of current point"
    do j=1, natom
       write(1000,*)  label(j), (x(k,j), k=1,ndim)
    end do
    close(1000)
    !--------------------------------
    !Run Molpro
    call SYSTEM("molpro -s -n "//procstring//" "//trim(basename)//"_force.com")
    !--------------------------------
    !Search through output and find the energy
    open(unit=2000, file=trim(basename)//"_force.out", status="OLD", access="SEQUENTIAL")
    linecount=0
    do
       read(2000,*,END=10)
       linecount=linecount+1
    end do
10  rewind(2000)

    ! do i=1,linecount-3
    !    read(2000,*)
    ! end do
    ! read(2000,*) energy, throwaway
    energy=0.0d0
    grad(:,:)=0.0d0
    do i=1,linecount-6
       read(2000,*)
    end do
    read(2000,informat) intext,energy
    close(2000)
    energy=energy-V0
    write(*,*)energy
    !--------------------------------
    !Search through output and find the energy
    open(unit=3000, file="forces.xyz", status="OLD", access="SEQUENTIAL")
    read(3000,*)
    read(3000,*)
    do i=1,natom
       read(3000,*) throw0,throw1,throw2,throw3,throw4, grad(1,i), grad(2,i), grad(3,i)
       grad(1,i)=grad(1,i)*0.529177d0*3.6749d-2
       grad(2,i)=grad(2,i)*0.529177d0*3.6749d-2
       grad(3,i)=grad(3,i)*0.529177d0*3.6749d-2
    end do
    close(3000)
    call CHDIR("..")
    
    return
  end subroutine Potforce
!---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1, eps
    integer::              i, j
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
    deallocate(procdir)
    call SYSTEM("cd ..")
  end subroutine V_finalize
  
end module mcmod_mass
