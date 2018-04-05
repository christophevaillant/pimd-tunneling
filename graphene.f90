!Module containing the LCBOP potential for C-C interactions
!with periodic boundary conditions set up for graphene
!Potential from: Los and Fasolino, PRB 68, 024107 (2003)
module graphenemod
  implicit none

  !all the functions
  double precision::     pot, potij, potsr, potlr, cutoff, switch
  double precision::     heaviside, Bij, periodpot

  !general purpose variables
  integer::              i, j, k, natoms

  !-------------------------------------------------------------
  !Parameters for the potential

  !-------------------------------------------------------------
  private
  public::               pot, periodpot

contains

  !-------------------------------------------------------------
  !Potential with periodic boundary conditions
  subroutine periodpot(x, rcut, a, V)
    double precision, intent(in):: x(:,:), rcut, a
    double precision, intent(out):: V

    V=pot(x)

    return
  end subroutine periodpot
    
  !-------------------------------------------------------------
  !Potential without periodic boundary conditions
  function pot(x)
    implicit none
    double precision, intent(in)::   x(:,:)
    double precision::   r
    
    pot=0.0d0
    !Ewald sum
    do i=1,natoms
       do j=1, i
          pot= pot + 0.5*potij(x, i,j)
       end do
    end do

    return
  end function pot

  !-------------------------------------------------------------
  !Two-body potential
  function potij(r)
    implicit none
    double precision, intent(in):: x(:,:)
    integer, intent(in)::          i,j
    double precision::    cut, r

    r= sqrt(sum(x(i,:)**2- x(j,:)**2))
    cut=cutoff(r)
    potij= cut*potsr(x,i,j) + (1.0d0-cut)*potlr(r)

    return
  end function potij

  !-------------------------------------------------------------
  !Cutoff function
  function cutoff(r)
    implicit none
    double precision, intent(in):: r

    cutoff= heaviside(-r) + &
         heaviside(r)*heaviside(1.0d0-r)*exp(gamma*r**3/(r**3 - 1.0d0))
    return
  end function cutoff

  !-------------------------------------------------------------
  !Heaviside step function
  function heaviside(r)
    implicit none
    double precision, intent(in):: r

    heaviside= 0.0d0
    if (x .ge.0.0d0) heaviside=1.0d0
    return
  end function heaviside

  !-------------------------------------------------------------
  !Short range potential
  function potsr(x,i,j)
    implicit none
    double precision, intent(in):: x
    integer, intent(in)::          i,j
    double precision::             r

    r= sqrt(sum(x(i,:)**2- x(j,:)**2))
    potsr= A*exp(-alpha*r) - Bij(x,i,j)*(B1*exp(-beta1*r) + B2*exp(-beta2*r))
  end function potsr

  !-------------------------------------------------------------
  !Bond order

  function Bij(x,i,j)
    implicit none
    double precision, intent(in):: x
    integer, intent(in)::          i,j
    double precision::             r, rij, rik, thetaijk

  end function Bij
    


end module graphenemod
