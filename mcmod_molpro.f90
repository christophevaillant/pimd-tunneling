module mcmod_mass
  implicit none
  double precision::               V0, eps2=1.0d-3
  integer::                        atom1, atom2, atom3
  character(len=:), allocatable::  procdir, beaddir(:)
  character, allocatable::         label(:)
  character(len=20)::              basename
  character(len=15)::              informat
  character(len=2)::               procstring
  integer::                        n, ndim, ndof, natom, xunit, totdof, potprocs
  logical::                        potforcepresent
  
contains

  subroutine file_init(force,i)
    implicit none
    integer, intent(in), optional:: i
    logical, intent(in):: force
    character(len=25)::  wfufilename, fulldir

    if (present(i)) then
       fulldir= trim(procdir) // "/" // trim(beaddir(i))
    else
       fulldir= trim(procdir)
    end if
    
    if (force) then
       call SYSTEM("cp " //trim(basename)// "_force.com "//trim(fulldir))
       open(50, file=trim(fulldir)//"/"//trim(basename)//"_force.com", status="old", position="append")
    else
       call SYSTEM("cp " //trim(basename)// ".com "//trim(fulldir))
       open(50, file=trim(fulldir)//"/"//trim(basename)//".com", status="old", position="append")
    end if
    ! if (present(i)) then
    !    write(wfufilename,"(A14,I3,A4)") "acetylcyanide_", i, ".wfu"
    !    write(50,*) "FILE,2," // trim(wfufilename) // ",UNKNOWN"
       ! write(50,"(A6,I4,A2)")"start,", 7000+i,".2"
       ! write(50,*)
       ! write(50,"(A8,I4,A2)")"orbital,",7000+i,".2"
    ! else
    !    wfufilename="acetylcyanide.wfu"
    !    write(50,*) "FILE,2," // trim(wfufilename) // ",UNKNOWN"
    !    ! write(50,"(A6,I4,A2)")"start,", 7000,".2"
    !    ! write(50,*)
    !    ! write(50,"(A8,I4,A2)")"orbital,",7000,".2"
    ! end if
    ! close(50)
    return
  end subroutine file_init
  
  subroutine V_init(iproc)
    integer, intent(in):: iproc
    character(len=7)::  format_string
    logical::    ex
    integer::  i
    namelist /POTDATA/ basename, atom1,atom2,atom3, informat, potprocs

    !-----------------------
    !read in default parameters
    atom1=1
    atom2=2
    atom3=3
    potforcepresent= .false.
    informat="(A22,3X,F18.12)"
    potprocs=1
    
    read(5, nml=POTDATA)

    !-----------------------
    !format statement for reading in
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
    !-----------------------
    !make directory structure
    inquire(file=procdir, exist=ex)
    if (.not. ex) then
       call SYSTEM("mkdir "//procdir)
    end if

    allocate(character(len=8)::beaddir(n))
    do i=1, n
       if (i .lt. 10) then
          format_string = "(A4,I1)"
       else if (i .lt. 100) then
          format_string = "(A4,I2)"
       else if (i .lt. 1000) then
          format_string = "(A4,I3)"
       end if
       write(beaddir(i),format_string) "bead", i
       call SYSTEM("mkdir "//procdir//"/"//beaddir(i))
    end do
    
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function V(x, bead)
    implicit none
    integer, intent(in), optional:: bead
    double precision, intent(in)::     x(:,:)
    double precision::  v,throwaway
    double precision, allocatable:: dummy1(:),dummy2(:), xtemp(:)
    character(len=200)::   intext
    character(len=25)::    inword,fulldir
    integer::              i,j,k, ierr, linecount

    !--------------------------------
    !Output geometry to geometry file
    if (present(bead)) then
       call file_init(.false., bead)
    else
       call file_init(.false.)
    end if
    
    if (present(bead)) then
       fulldir= trim(procdir) // "/" // trim(beaddir(bead))
    else
       fulldir= trim(procdir)
    end if
    
    call CHDIR(trim(fulldir))
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
  subroutine Vprime(x, grad, bead)
    implicit none
    integer, intent(in), optional:: bead
    double precision, intent(in)::  x(:,:)
    double precision, intent(out):: grad(:,:)
    integer::              i,j,k,linecount
    character::            throw0
    double precision::     dummy1, eps, potplus, potminus
    ! double precision::     throw1,throw2,throw3,throw4
    double precision, allocatable:: xtemp(:,:)

    allocate(xtemp(ndim,natom))
    eps=1d-2
    do i= 1, ndim
       do j= 1, natom
          xtemp(:,:)= x(:,:)
          xtemp(i,j)= xtemp(i,j) + eps
          potplus=V(xtemp,bead)
          xtemp(i,j)= xtemp(i,j) - 2.0d0*eps
          potminus=V(xtemp,bead)
          grad(i,j)= (potplus-potminus)/(2.0d0*eps)          
       end do
    end do
    deallocate(xtemp)
    ! if (present(bead)) then
    !    call file_init(.true., bead)
    ! else
    !    call file_init(.true.)
    ! end if
    ! !--------------------------------
    ! !Output geometry to geometry file
    ! call CHDIR(procdir)
    ! open(1000, file="geometry.xyz")
    ! write(1000,*) natom
    ! write(1000,*) "Geometry of current point"
    ! do j=1, natom
    !    write(1000,*)  label(j), (x(k,j), k=1,ndim)
    ! end do
    ! close(1000)
    ! !--------------------------------
    ! !Run Molpro
    ! call SYSTEM("molpro -s -n "//procstring//" "//trim(basename)//"_force.com")
    ! !--------------------------------
    ! !Search through output and find the energy
    ! open(unit=3000, file="forces.xyz", status="OLD", access="SEQUENTIAL")
    ! read(3000,*)
    ! read(3000,*)
    ! do i=1,natom
    !    read(3000,*) throw0,throw1,throw2,throw3,throw4, grad(1,i), grad(2,i), grad(3,i)
    !    grad(1,i)=-grad(1,i)!*0.529177d0*3.6749d-2
    !    grad(2,i)=-grad(2,i)!*0.529177d0*3.6749d-2
    !    grad(3,i)=-grad(3,i)!*0.529177d0*3.6749d-2
    ! end do
    ! close(3000)
    ! call CHDIR("..")
    
    return
  end subroutine Vprime

    !---------------------------------------------------------------------
  subroutine potforce(x,grad,energy, bead)
    implicit none
    double precision,intent(in)::     x(:,:)
    integer, intent(in), optional:: bead
    double precision, intent(out)::  grad(:,:), energy
    double precision::               throwaway
    double precision, allocatable::  dummy1(:),dummy2(:), xtemp(:)
    character(len=200)::   intext
    character(len=20)::    inword1, inword2, inword3
    integer::              i,j,k, ierr, linecount
    character::            throw0
    double precision::     throw1,throw2,throw3,throw4

    if (present(bead)) then
       call file_init(.true., bead)
    else
       call file_init(.true.)
    end if

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
  subroutine  Vdoubleprime(x,hess, bead)
    implicit none
    integer, intent(in), optional:: bead
    double precision::     hess(:,:,:,:), x(:,:), dummy1, eps
    integer::              i, j,k,l
    double precision:: pot1, pot2, pot3, pot4
    double precision, allocatable:: xtemp(:,:)
    ! double precision, allocatable::     gradplus(:, :), gradminus(:, :)

    allocate(xtemp(ndim,natom))
    eps=1d-2
    do i= 1, ndim
       do j= 1, natom
          do k=i, ndim
             do l=j,natom
                xtemp(:,:)= x(:,:)
                xtemp(i,j)= xtemp(i,j) + eps
                xtemp(k,l)= xtemp(k,l)+ eps
                pot1=V(xtemp,bead)
                xtemp(k,l)= xtemp(k,l) - 2.0d0*eps
                pot2=V(xtemp,bead)
                xtemp(i,j)= xtemp(i,j) - 2.0d0*eps
                xtemp(k,l)= xtemp(k,l) + 2.0d0*eps
                pot3=V(xtemp,bead)                
                xtemp(k,l)= xtemp(k,l) - 2.0d0*eps
                pot4=V(xtemp,bead)                
                hess(i,j,k,l)= (pot1 - pot2- pot3 + pot4)/(4.0d0*eps)
                hess(k,l,i,j)= hess(i,j,k,l)
             end do
          end do
       end do
    end do

    deallocate(xtemp)
    ! allocate(gradplus(ndim,natom), gradminus(ndim,natom))
    ! eps=1d-2
    ! do i= 1, ndim
    !    do j= 1, natom
    !       x(i,j)= x(i,j) + eps
    !       call Vprime(x, gradplus)
    !       x(i,j)= x(i,j) - 2.0d0*eps
    !       call Vprime(x, gradminus)
    !       x(i,j)= x(i,j) + eps
    !       hess(i,j,:,:)= (gradplus(:,:)-gradminus(:,:))/(2.0d0*eps)          
    !    end do
    ! end do

    ! deallocate(gradplus,gradminus)
    return
  end subroutine Vdoubleprime

  subroutine V_finalize()
    deallocate(procdir)
    call SYSTEM("cd ..")
  end subroutine V_finalize
  
end module mcmod_mass
