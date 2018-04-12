!Module containing the LCBOP potential for C-C interactions
!with periodic boundary conditions set up for graphene
!Potential from: Los and Fasolino, PRB 68, 024107 (2003)
module graphenemod
  implicit none

  !general purpose variables
  integer::              natoms
  double precision, parameter::    pi=3.14159265358979d0
  !-------------------------------------------------------------
  !Parameters for the potential
  
  double precision:: r1, r2, gamma, A, B1, B2, alpha, beta1, beta2
  double precision:: dH, C1, C4, C6, LH, kappa, P0, P1 !P0 and P1 are R0 and R1 in the paper
  double precision:: r0, r1lr, r2lr, epsilon1, epsilon2, lambda1, lambda2

  !Energy params (eV)
  data A /35652.94452d0/, B1 /18614.83652d0/, B2 /32.01993977d0/
  data epsilon1 /6.093133d0/, epsilon2 /2.617755d0/
  !Distance params (angstroms)
  data r1 /1.7d0/, r2 /2.3d0/, gamma /1.8d0/
  data alpha /6.26781252d0/, beta1 /5.83045680/, beta2 /1.16864228d0/
  data dH /0.14d0/, C1 /3.3d0/, C4 /220.0d0/, C6 /-5434.715d0/
  data LH /0.688316d0/, kappa /1.619070d0/, P0 /1.612316d0/, P1 /5.485568d0/
  data r0 /3.716163d0/, r1lr /5.5d0/, r2lr /6.d0/
  data lambda1 /1.359381d0/, lambda2 /2.073944d0/

  !Interpolation of G function parameters
  double precision::     gspline(6),gfunc(6),costheta(6)
  data gfunc(:)/5.48948d-3,8.188859d-2,1.5709129d-1,7.772d-1,6.780d0,24.40d0/
  data costheta(:)/-1.0d0, -0.5d0,-0.3333333333333333333333d0,0.0d0,0.5d0,1.0d0/

  !Interpolation of Fconj matrix
  double precision::    Fconj0(4,4),Fconj1(4,4)
  data Fconj0(1,:) /0.0d0, 0.034993d0, -0.009085d0, -0.229403d0/
  data Fconj0(2,:) /0.034993d0, 0.0d0, -0.058546d0, -0.147668d0/
  data Fconj0(3,:) /-0.009085d0, -0.058546d0, 0.0d0, -0.083991d0/
  data Fconj0(4,:) /-0.229403d0, -0.147667d0, -0.083991d0, 0.0d0/

  data Fconj1(1,:) /0.0d0, 0.100921d0, 0.071525d0, -0.229403d0/
  data Fconj1(2,:) /0.100921d0, 0.239564d0, 0.010324d0, -0.147667d0/
  data Fconj1(3,:) /0.071525d0, 0.010324d0, 0.161180d0, -0.083991d0/
  data Fconj1(4,:) /-0.229403d0, -0.147667d0, -0.083991d0, 0.0d0/

  !-------------------------------------------------------------
  private
  public::               graphenepot, carbon_init

contains
  !-------------------------------------------------------------
  !Initialization of the potential
  subroutine carbon_init()
    implicit none

    !Set up the spline for the G function
    call spline(costheta, gfunc, 0.0d0, 0.0d0, gspline)
    
  end subroutine carbon_init

  !-------------------------------------------------------------
  !Potential with periodic boundary conditions
  !Assumes 6 atoms per unit cell
  subroutine graphenepot(x, rcut, V)
    double precision, intent(inout):: x(18)
    double precision, intent(in)::  rcut
    double precision, intent(out):: V
    double precision, allocatable:: superx(:,:), finalx(:,:), tempx(:,:)
    double precision::   rij, a1(3), a2(3), a, yl,yu, xcopy(18)
    data a /2.4595d0/
    data a1(:) /3.6893d0, -2.13d0, 0.0d0/
    data a2(:) /3.6893d0, 2.13d0, 0.0d0/
    integer:: ncells,cells, i, j,k, m
    character(len=100)::  stringnumber, filename

    xcopy(:)= x(:)
    V=0.0d0
    !Work out if coordinates are in the unit cell
    do i=1,6
       if (x(3*(i-1)+1) .le. 0) then
          yl= x(3*(i-1)+1)*(a1(2)/a1(1))+a1(2)
          yu= x(3*(i-1)+1)*(a2(2)/a2(1))+a2(2)
          ! write(*,*) "Left:"
          ! write(*,*) x(3*(i-1)+1), x(3*(i-1)+2), yl, yu
          !if upper left of unit cell
          if (x(3*(i-1)+2) .gt. yu) x(3*(i-1)+1:3*(i-1)+3)= &
               x(3*(i-1)+1:3*(i-1)+3)+ a1(1:3)
          !if lower left of unit cell
          if (x(3*(i-1)+2) .lt. yl) x(3*(i-1)+1:3*(i-1)+3)= &
               x(3*(i-1)+1:3*(i-1)+3)+ a2(1:3)
       else if (x(3*(i-1)+1) .gt. 0) then
          yl= (x(3*(i-1)+1) - a2(1))*(a2(2)/a2(1)) + a1(2) + a2(2)
          yu= (x(3*(i-1)+1) - a1(1))*(a1(2)/a1(1)) + a1(2) + a2(2)
          ! write(*,*) "Right:"
          ! write(*,*) x(3*(i-1)+1), x(3*(i-1)+2), yl, yu
          !if upper right of unit cell
          if (x(3*(i-1)+2) .gt. yu) x(3*(i-1)+1:3*(i-1)+3)= &
               x(3*(i-1)+1:3*(i-1)+3)- a2(1:3)
          !if lower right of unit cell
          if (x(3*(i-1)+2) .lt. yl) x(3*(i-1)+1:3*(i-1)+3)= &
               x(3*(i-1)+1:3*(i-1)+3)- a1(1:3)
       end if
    end do

    ! open(121, file="unitcell.xyz")
    ! write(121,*) 6
    ! write(121,*) "unit cell"
    ! do k=1, 6
    !    write(121,*) "C", (x(3*(k-1)+j), j=1,3)
    ! end do
    ! close(121)

    !Need to work out how many unit cells
    !ncells is total number of cells, cells is length of grid
    cells= (nint(rcut/a) + 1)
    if (mod(cells,2).eq.0.0) cells=cells+1
    ncells= cells**2

    !Now work out coordinates of all atoms and work out which ones are relevent
    allocate(superx(ncells*6,3))
    do i=1, cells
       do j=1,cells
          do k=1,6
             superx(6*(cells*(j-1) + i-1)+k, 1)= x(3*(k-1)+1) + &
                  (i-1-(cells/2))*a1(1) + (j-1-(cells/2))*a2(1)
             superx(6*(cells*(j-1) + i-1)+k, 2)= x(3*(k-1)+2) + &
                  (i-1-(cells/2))*a1(2) + (j-1-(cells/2))*a2(2)
             superx(6*(cells*(j-1) + i-1)+k, 3)= x(3*(k-1) + 3) !considering only planar
          end do
       end do
    end do

    ! open(120, file="supercell.xyz")
    ! write(120,*) ncells*6
    ! write(120,*) "unit cell"
    ! do k=1, ncells*6
    !    write(120,*) "C", (superx(k,j), j=1,3)
    ! end do
    ! close(120)
    ! stop

    !Loop over cells in the unit cell
    allocate(tempx(ncells*6,3))
    do i=1, 6
       tempx(:,:)=0.0d0
       natoms=0
       !loop over final atoms (starts at i+1 to avoid double counting)
       do j=1, 6*ncells
             rij= distance(x(3*(i-1)+1:3*(i-1)+3), superx(j,1:3))
             ! write(*,*)i,j, "---------"
             ! write(*,*) x(3*(i-1)+1:3*(i-1)+3)
             ! write(*,*) superx(j,1:3)
             if (rij .lt. rcut) then
                natoms=natoms+1
                tempx(natoms,:)= superx(j,:)
             !    write(*,*) "Accepted:", rij
             ! else
             !    write(*,*) "Rejected:", rij
             end if
       end do
       ! write(*,*) i,natoms, rcut
       allocate(finalx(natoms, 3))
       finalx(1:natoms,:)= tempx(1:natoms,:)

       ! write(stringnumber,'(I1)') i
       ! filename= trim("cell_" // trim(stringnumber) // ".xyz")
       ! open(120, file=filename)
       ! write(120,*) natoms
       ! write(120,*) "supercell"
       ! do k=1, natoms
       !    write(120,*) "C", (finalx(k,j), j=1,3)
       ! end do
       ! close(120)

       V=V+ pot(finalx)*23.061D0

       deallocate(finalx)
       x(:)=xcopy(:)
    end do
    ! stop
    deallocate(superx,tempx)
    return
  end subroutine graphenepot
    
  !-------------------------------------------------------------
  !Potential without periodic boundary conditions
  function pot(x)
    implicit none
    double precision, intent(in)::   x(:,:)
    double precision::   r, pot
    integer::            i,j
    
    pot=0.0d0
    !Ewald sum
    do i=1,natoms
       do j=i+1, natoms
          r= distance(x(i,:), x(j,:))
          pot= pot + 0.5*potij(x, i,j)
       end do
    end do
    return
  end function pot

  !-------------------------------------------------------------
  !Two-body potential
  function potij(x,i,j)
    implicit none
    double precision, intent(in):: x(:,:)
    integer, intent(in)::          i,j
    integer::           k
    double precision::    cut, r, potij, Elr, cut2, rlr

    r= distance(x(i,:), x(j,:))
    cut=cutoff(r)
    potij= cut*potsr(x,i,j) + (1.0d0-cut)*potlr(x,i,j)
    return
  end function potij

  !-------------------------------------------------------------
  !Cutoff function
  function cutoff(r)
    implicit none
    double precision, intent(in):: r
    double precision:: x, cutoff
    
    x= (r-r1)/(r2-r1)
    ! if (abs(x**3 -1.0d0) .gt. 5D-3) then
    cutoff= heaviside(-x) + &
         heaviside(x)*heaviside(1.0d0-x)*exp(gamma*x**3/(x**3 - 1.0d0))
    ! else
    !    cutoff=1.0d10*x
    ! end if
    ! write(*,*) x, (x**3 - 1.0d0),gamma*x**3
    return
  end function cutoff

  !-------------------------------------------------------------
  !Heaviside step function
  function heaviside(r)
    implicit none
    double precision, intent(in):: r
    double precision:: heaviside

    heaviside= 0.0d0
    if (r .ge.0.0d0) heaviside=1.0d0
    return
  end function heaviside

  !-------------------------------------------------------------
  !Short range potential
  function potsr(x,i,j)
    implicit none
    double precision, intent(in):: x(:,:)
    integer, intent(in)::          i,j
    double precision::             r, potsr

    r= distance(x(i,:), x(j,:))
    potsr= A*exp(-alpha*r) - &
         0.5d0*(Bij(x,i,j)+ Bij(x,j,i) + Fconj(x,i,j))*&
         (B1*exp(-beta1*r) + B2*exp(-beta2*r))
    return
  end function potsr

  !-------------------------------------------------------------
  !Bond order

  function Bij(x,i,j)
    implicit none
    double precision, intent(in):: x(:,:)
    integer, intent(in)::          i,j
    double precision::             r, rij, rik, rjk,thetaijk
    double precision::             ksum, Bij
    integer::                      k

    ksum=0.0d0
    rij= distance(x(i,:), x(j,:))
    do k=1, natoms
       if ((k.eq.i) .or. (k.eq.j)) cycle
       rik= distance(x(i,:), x(k,:))
       rjk= distance(x(j,:), x(k,:))
       thetaijk= 0.5d0*(rjk**2 + rij**2 - rik**2)/(rjk*rij) !cosine rule
       ksum=ksum+ cutoff(rik)*splint(costheta, gfunc,gspline, thetaijk)&
            *Hfunc(rij-rik)
    end do
    Bij= 1.0d0/sqrt(1.0d0+ ksum)

    return
  end function Bij

  !-------------------------------------------------------------
  !H functions
  function Hfunc(x)
    implicit none
    double precision, intent(in):: x
    double precision:: Hfunc
    
    Hfunc=0.0d0

    if (x .lt. -dH) then
       Hfunc= (1.0/(1.0 + (kappa*(x+dH))**10.0d0))**0.1d0
       Hfunc= LH*(1.0d0 + kappa*(x+dH)*Hfunc)
    else if ((-dH .le. x) .and. (x.le.dH)) then
       Hfunc= 1.0d0 + C1*x + 0.5d0*(C1*x)**2 + C4*x**4 + C6*x**6
    else
       Hfunc= P0 + P1*(x-dH)
    end if
    
    return
  end function Hfunc

  !-------------------------------------------------------------
  !Utility to work out a radial distance
  function distance(ri,rj)
    implicit none
    double precision, intent(in):: ri(:), rj(:)
    double precision:: distance
    integer:: i
    
    distance=0.0d0
    do i=1,3
       distance=distance + (ri(i)- rj(i))**2
    end do
    if (distance .ne. 0.0d0) distance=sqrt(distance)
    ! if (distance .ne. distance) then
    !    write(*,*) "Nan in distance"
    !    write(*,*) ri(:)
    !    write(*,*) rj(:)
    !    do i=1,3
    !       write(*,*) (ri(i)- rj(i))**2
    !    end do
    ! end if
    return
  end function distance

  !-------------------------------------------------------------
  !Fconj function which is a horrible faff
  function Fconj(x, i,j)
    implicit none
    double precision, intent(in):: x(:,:)
    integer, intent(in)::  i,j
    double precision::   Ni, Nj, Nconj, Nij, Nji
    double precision::   Nelij, Nelji,rij,rik
    double precision:: Mij, Mji, Fconj
    integer::  k,l

    rij= distance(x(i,:), x(j,:))
    Ni= coordination(x,i)
    Nj= coordination(x,j)
    Nij= min(3.0d0, Ni-cutoff(rij))
    Nji= min(3.0d0, Nj-cutoff(rij))

    Mij= 0.0d0
    do k=1,natoms
       if ((k.eq.i) .or. (k.eq.j)) cycle
       rik= distance(x(i,:), x(k,:))
       Mij= Mij+ cutoff(rik)*Fij(rik)
    end do
    Mij= min(3.0d0, Mij)

    Mji= 0.0d0
    do k=1,natoms
       if ((k.eq.i) .or. (k.eq.j)) cycle
       rik= distance(x(j,:), x(k,:))
       Mji= Mji+ cutoff(rik)*Fij(rik)
    end do
    Mji= min(3.0d0, Mji)

    Nelij= (4.0d0-Mij)/(Nij+1.0-Mij)
    Nelji= (4.0d0-Mji)/(Nji+1.0-Mji)

    Nconj= (Nij+1.0d0)*(Nji+1.0d0)*(Nelij+Nelji) - 4.0d0*(Nij+Nji+2.0d0)
    Nconj= Nconj/(Nij*(3.0d0-Nij)*(Nji+1.0d0) + Nji*(3.0d0-Nji)*(Nij+1.0d0) + 1d-10)


    Fconj= (1.0d0- Nconj)*Fconjeval(Nij,Nji,0) + Nconj*Fconjeval(Nij,Nji,1)
  end function Fconj
  !-------------------------------------------------------------
  function Fconjeval(Nij,Nji,m)
    implicit none
    integer, intent(in):: m
    double precision, intent(in):: Nij, Nji
    double precision:: x, y, Fconjeval
    integer:: k,l
    
    x= Nij - int(Nij)
    y= Nji - int(Nji)
    k= int(Nij)
    l= int(Nji)

    if ((x.gt.0 .and. x.lt.1) .or. (y.gt.0 .and. y.lt.1)) then
       if (m.eq.0) then
          Fconjeval= (1.0d0-y)*(1.0d0-x)*(Fconj0(max(1,min(4,k+1)),max(1,min(4,l+1))) + x**2*Fconjtilde(k,l,1,0,m) +&
               y**2*Fconjtilde(k,l,0,1,m))
          Fconjeval= Fconjeval + (1.0d0-y)*x*(Fconj0(max(1,min(4,k+2)),max(1,min(4,l+1))) + (1.0d0-x**2)*Fconjtilde(k,l,0,0,m) +&
               y**2*Fconjtilde(k,l,1,1,m))
          Fconjeval= Fconjeval + (1.0d0-x)*y*(Fconj0(max(1,min(4,k+1)),max(1,min(4,l+2))) + (1.0d0-y**2)*Fconjtilde(k,l,0,0,m) +&
               x**2*Fconjtilde(k,l,1,1,m))
          Fconjeval= Fconjeval + x*y*(Fconj0(max(1,min(4,k+2)),max(1,min(4,l+2))) + (1.0d0-y**2)*Fconjtilde(k,l,1,0,m) +&
               (1.0d0-x**2)*Fconjtilde(k,l,0,1,m))
       else
          Fconjeval= (1.0d0-y)*(1.0d0-x)*(Fconj1(max(1,min(4,k+1)),max(1,min(4,l+1))) + x**2*Fconjtilde(k,l,1,0,m) +&
               y**2*Fconjtilde(k,l,0,1,m))
          Fconjeval= Fconjeval + (1.0d0-y)*x*(Fconj1(max(1,min(4,k+2)),max(1,min(4,l+1))) + (1.0d0-x**2)*Fconjtilde(k,l,0,0,m) +&
               y**2*Fconjtilde(k,l,1,1,m))
          Fconjeval= Fconjeval + (1.0d0-x)*y*(Fconj1(max(1,min(4,k+1)),max(1,min(4,l+2))) + (1.0d0-y**2)*Fconjtilde(k,l,0,0,m) +&
               x**2*Fconjtilde(k,l,1,1,m))
          Fconjeval= Fconjeval + x*y*(Fconj1(max(1,min(4,k+2)),max(1,min(4,l+2))) + (1.0d0-y**2)*Fconjtilde(k,l,1,0,m) +&
               (1.0d0-x**2)*Fconjtilde(k,l,0,1,m))
       end if
    else if (m.eq.0) then
       Fconjeval= Fconj0(max(1,min(4,k+1)),max(1,min(4,l+1)))
    else if (m.eq.1) then
       Fconjeval= Fconj1(max(1,min(4,k+1)),max(1,min(4,l+1)))
    end if
    return
  end function Fconjeval
  !-------------------------------------------------------------
  function Fconjtilde(k,l,n,m,p)
    implicit none
    integer, intent(in):: k,l,n,m,p
    double precision:: Fconjtilde

    if ((k.eq.0) .or. (k.eq.3)) then
       Fconjtilde=0.0d0
    else if (p .eq. 0) then
       Fconjtilde= 0.5d0*(Fconj0(max(1,min(4,k+2+n)),max(1,min(4,l+m))) - &
            Fconj0(max(1,min(4,k-2+n)),max(1,min(4,l+m))))
    else if (p.eq.1) then
       Fconjtilde= 0.5d0*(Fconj1(max(1,min(4,k+2+n)),max(1,min(4,l+m))) - &
            Fconj1(max(1,min(4,k-2+n)),max(1,min(4,l+m))))
    end if
    if (p .eq. 0) then
       Fconjtilde= Fconjtilde - Fconj0(max(1,min(4,k+2)), max(1,min(4,l+m+1))) + Fconj0(max(1,min(4,k+1)), max(1,min(4,l+m+1)))
       Fconjtilde= Fconjtilde*(-1.0d0)**(k+n)
    else if (p.eq.1) then
       Fconjtilde= Fconjtilde - Fconj1(max(1,min(4,k+2)), max(1,min(4,l+m+1))) + Fconj1(max(1,min(4,k+1)), max(1,min(4,l+m+1)))
       Fconjtilde= Fconjtilde*(-1.0d0)**(l+m)
    end if

    return
  end function Fconjtilde
  !-------------------------------------------------------------
  function coordination(x,i)
    implicit none
    double precision, intent(in):: x(:,:)
    integer, intent(in):: i
    integer:: k
    double precision::  rik, coordination

    coordination=0.0d0
    do k=1,natoms
       rik= distance(x(i,:), x(k,:))
       coordination= coordination+ cutoff(rik)
    end do
    return
  end function coordination

  !-------------------------------------------------------------
  function Fij(r)
    implicit none
    double precision, intent(in):: r
    double precision:: Fij

    Fij= heaviside(r-3.0d0) + 0.5d0*heaviside(r-2.0d0)*&
         heaviside(3.0d0-r)*(1.0d0-cos(pi*r-2.0d0))
    return
  end function Fij
  !-------------------------------------------------------------
  function potlr(x,i,j)
    implicit none
    double precision, intent(in):: x(:,:)
    integer, intent(in):: i,j
    double precision:: rij, potlr
    
    rij= distance(x(i,:), x(j,:))
    potlr= heaviside(r0-rij)*Vmorse(rij,1) + &
         heaviside(rij-r0)*Vmorse(rij,2)
    return
  end function potlr
  !-------------------------------------------------------------
  function Vmorse(r, p)
    implicit none
    double precision, intent(in):: r
    integer, intent(in):: p
    double precision:: Vmorse
    
    if (p.eq.1) then
       Vmorse= epsilon1*(exp(-2.0d0*lambda1*(r-r0)) &
            - 2.0d0*exp(-lambda1*(r-r0))) + epsilon1-epsilon2
    else
       Vmorse= epsilon2*(exp(-2.0d0*lambda2*(r-r0)) &
            - 2.0d0*exp(-lambda2*(r-r0)))
    end if
    return
  end function Vmorse

  !---------------------------------------------------------------------
  !spline algorithms
  FUNCTION assert_eq(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq
    if (n1 == n2 .and. n2 == n3) then
       assert_eq=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq'
    end if
  END FUNCTION assert_eq

  SUBROUTINE spline(x,y,yp1,ypn,y2)
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(IN) :: x,y
    double precision, INTENT(IN) :: yp1,ypn
    double precision, DIMENSION(:), INTENT(OUT) :: y2
    ! Given arrays x and y of length N containing a tabulated function, i.e., yi = f (xi ), with x1 <
    ! x2 < . . . < xN , and given values yp1 and ypn for the first derivative of the interpolating
    ! function at points 1 and N , respectively, this routine returns an array y2 of length N
    ! that contains the second derivatives of the interpolating function at the tabulated points
    ! xi . If yp1 and/or ypn are equal to 1 × 10 30 or larger, the routine is signaled to set the
    ! corresponding boundary condition for a natural spline, with zero second derivative on that
    ! boundary.
    INTEGER :: n
    double precision, DIMENSION(size(x)) :: a,b,c,r
    n=assert_eq(size(x),size(y),size(y2),"spline")
    c(1:n-1)=x(2:n)-x(1:n-1)
    ! Set up the tridiagonal equations.
    r(1:n-1)=6.0d0*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0d0*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    if (yp1 > 0.99d30) then
       ! The lower boundary condition is set either to be “natural”
       r(1)=0.0
       c(1)=0.0
    else
       ! or else to have a specified first derivative.
       r(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       c(1)=0.5
    end if
    ! The upper boundary condition is set either to be
    ! or else to have a specified first derivative.

    if (ypn > 0.99d30) then
       ! “natural”
       r(n)=0.0
       a(n)=0.0
    else
       r(n)=(-3.0d0/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
       a(n)=0.5
    end if
    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
    return
  END SUBROUTINE spline
  
  SUBROUTINE tridag(a,b,c,r,u)
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(IN) :: a,b,c,r
    double precision, DIMENSION(:), INTENT(OUT) :: u
    ! Solves for a vector u of size N the tridiagonal linear set given by equation (2.4.1) using a
    ! serial algorithm. Input vectors b (diagonal elements) and r (right-hand sides) have size N ,
    ! while a and c (off-diagonal elements) are size N − 1.
    double precision, DIMENSION(size(b)) :: gam
    ! One vector of workspace, gam is needed.
    INTEGER:: n,j
    double precision :: bet
    n=assert_eqn((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),"tridag_ser")
    bet=b(1)
    if (bet == 0.0) then
       write(*,*) "tridag_ser: Error at code stage 1"
       stop
    end if
    ! If this happens then you should rewrite your equations as a set of order N − 1, with u2
    ! trivially eliminated.
    u(1)=r(1)/bet
    do j=2,n
       ! Decomposition and forward substitution.
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j-1)*gam(j)
       if (bet == 0.0) then
          write(*,*) "tridag_ser: Error at code stage 2"
          stop
       end if
       u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1
       ! Backsubstitution.
       u(j)=u(j)-gam(j+1)*u(j+1)
    end do
  END SUBROUTINE tridag

  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) "nrerror: an assert_eq failed with this tag:", &
            string
       STOP "program terminated by assert_eqn"
    end if
  END FUNCTION assert_eqn

  FUNCTION splint(xa,ya,y2a,x)
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    double precision, INTENT(IN) :: x
    double precision :: splint
    ! Given the arrays xa and ya, which tabulate a function (with the xai "s in increasing or
    ! decreasing order), and given the array y2a, which is the output from spline above, and
    ! given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
    ! and y2a are all of the same size.
    INTEGER :: khi,klo,n
    double precision :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),"splint")
    klo=max(min(locate(xa,x),n-1),1)
    ! We will find the right place in the table by means of locate"s bisection algorithm. This is
    ! optimal if sequential calls to this routine are at random values of x. If sequential calls are in
    ! order, and closely spaced, one would do better to store previous values of klo and khi and
    ! test if they remain appropriate on the next call.
    khi=klo+1
    ! klo and khi now bracket the input value of x.
    h=xa(khi)-xa(klo)
    if (h == 0.0) then
       write(*,*) "bad xa input in splint"
       stop
    end if
    ! The xa"s must be distinct.
    a=(xa(khi)-x)/h
    ! Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0
  END FUNCTION splint

  FUNCTION locate(xx,x)
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(IN) :: xx
    double precision, INTENT(IN) :: x
    INTEGER :: locate
    ! Given an array xx(1:N ), and given a value x, returns a value j such that x is between
    ! xx(j ) and xx(j + 1) . xx must be monotonic, either increasing or decreasing. j = 0 or
    ! j = N is returned to indicate that x is out of range.
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    ! True if ascending order of table, false otherwise.
    jl=0
    ! Initialize lower
    ju=n+1
    ! and upper limits.
    do
       if (ju-jl <= 1) exit
       ! Repeat until this condition is satisfied.
       jm=(ju+jl)/2
       ! Compute a midpoint,
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
          ! and replace either the lower limit
       else
          ju=jm
          ! or the upper limit, as appropriate.
       end if
    end do
    if (x == xx(1)) then
       ! Then set the output, being careful with the endpoints.
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if
  END FUNCTION locate

end module graphenemod
