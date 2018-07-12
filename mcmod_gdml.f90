module mcmod_mass
  implicit none
  double precision::               V0
  integer,parameter::              atom1=1, atom2=2, atom3=3
  integer::                        n, ndim, ndof, natom, xunit, totdof
! Written by John Morgan, modified by Christophe Vaillant
! Wrapper for Gradient Domain Machine Learning.
! Calls out to the GDML libary.
!
! GDML_SETUP: called from KEYWORDS, sets state information
! GMIN_GDML_WRAPPER: called from POTENTIAL, calculates potential and gradients

INTERFACE
    SUBROUTINE GDML_SETUP() BIND (C)
        USE ISO_C_BINDING
        IMPLICIT NONE
    END SUBROUTINE

    SUBROUTINE GDML_PREDICT(GDMLDOF,CRDS, GRAD, ENRG) BIND (C)
        USE ISO_C_BINDING
        IMPLICIT NONE
        INTEGER (C_INT)::  GDMLDOF
        REAL (C_DOUBLE), DIMENSION(GDMLDOF) :: CRDS
        REAL (C_DOUBLE), DIMENSION(GDMLDOF) :: GRAD
        REAL (C_DOUBLE) :: ENRG
    END SUBROUTINE

    SUBROUTINE GDML_FINALIZE() BIND (C)
        USE ISO_C_BINDING
        IMPLICIT NONE
    END SUBROUTINE
END INTERFACE

CONTAINS
    !---------------------------------------------------------------------------
    !
    ! Calls the setup function in the GDML library. The configuration file must
    ! be present in the current working directory. This also actas as a check
    ! that the file is present and that GDML has been compiled in.
    !
    !---------------------------------------------------------------------------

  subroutine V_init()
    V0=0.0d0
    call GDML_SETUP()
    return
  end subroutine V_init


    !---------------------------------------------------------------------------
    !
    ! Calls the GDML potential to calculate the energy and gradients.
    !
    ! CRDS: coordinates provided
    ! GRAD: return for the gradients
    ! ENRG: return for the energy
    !
    !---------------------------------------------------------------------------
    function V(x)
      implicit none
      double precision::     v, x(:,:)
      double precision, allocatable:: dummy1(:),xtemp(:)
      integer::              i,j

      allocate(xtemp(natom*ndim),dummy1(natom*ndim))
      do i=1,ndim
         do j=1,natom
            xtemp((j-1)*ndim +i)= x(i,j)*0.529177d0
         end do
      end do


        CALL GDML_PREDICT(ndof, xtemp, dummy1, V)

      deallocate(xtemp,dummy1)
      return
    end function V

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision::     grad(:,:), x(:,:), dummy1
    double precision, allocatable:: gradtemp(:), xtemp(:)

    allocate(gradtemp(ndof), xtemp(natom*ndim))
    do i=1,ndim
       do j=1,natom
          xtemp((j-1)*ndim +i)= x(i,j)*0.529177d0
       end do
    end do
        CALL GDML_PREDICT(ndof,xtemp, gradtemp, dummy1)
    do i= 1,ndim
       do j=1,natom
          grad(i,j)= gradtemp(ndim*(j-1)+i)*1.59362d-3*0.529177d0
       end do
    end do
    deallocate(gradtemp, xtemp)
    return
  end subroutine Vprime

    !---------------------------------------------------------------------------
    !
    ! Called after we don't need to call the potential wrapper anymore. Cleans
    ! up the Cython interface.
    !
    !---------------------------------------------------------------------------
    SUBROUTINE GDML_FINALISE_WRAPPER()

        IMPLICIT NONE

        CALL GDML_FINALIZE()

    END SUBROUTINE GDML_FINALISE_WRAPPER

  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1, eps
    integer::              i, j
    double precision, allocatable::     gradplus(:, :), gradminus(:, :)

    eps=1d-4
    allocate(gradplus(ndim, natom), gradminus(ndim, natom))
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
    deallocate(gradplus, gradminus)
    return
  end subroutine Vdoubleprime


  END MODULE MCMOD_MASS
