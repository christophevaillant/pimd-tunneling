module mcmod_mass
  implicit none
  double precision::               V0
  integer,parameter::              atom1=1, atom2=2, atom3=3
  integer::                        n, ndim, ndof, natom, xunit, totdof
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

    SUBROUTINE GDML_PREDICT(CRDS, GRAD, ENRG) BIND (C)
        USE COMMONS, ONLY: NATOMS
        USE ISO_C_BINDING
        IMPLICIT NONE
        REAL (C_DOUBLE), DIMENSION(NATOMS) :: CRDS
        REAL (C_DOUBLE), DIMENSION(NATOMS) :: GRAD
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
    SUBROUTINE GDML_SETUP_WRAPPER ()

        IMPLICIT NONE

        CALL GDML_SETUP()

    END SUBROUTINE GDML_SETUP_WRAPPER

    !---------------------------------------------------------------------------
    !
    ! Calls the GDML potential to calculate the energy and gradients.
    !
    ! CRDS: coordinates provided
    ! GRAD: return for the gradients
    ! ENRG: return for the energy
    !
    !---------------------------------------------------------------------------
    SUBROUTINE GDML_PREDICT_WRAPPER(CRDS, GRAD, ENRG)

        USE COMMONS, ONLY : NATOMS

        IMPLICIT NONE

        DOUBLE PRECISION, DIMENSION(3*NATOMS), INTENT(IN) :: CRDS
        DOUBLE PRECISION, DIMENSION(3*NATOMS), INTENT(OUT) :: GRAD
        DOUBLE PRECISION, INTENT(OUT) :: ENRG

        CALL GDML_PREDICT(CRDS, GRAD, ENRG)

    END SUBROUTINE GDML_PREDICT_WRAPPER

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

END MODULE GDML_MOD
