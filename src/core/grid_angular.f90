MODULE GRID_ANGULAR
!> Angular (theta,phi) grid: Gauss-Legendre quadrature in theta (built on
!> FORTVSH's SHGLQ), uniform quadrature in phi. Truncation is triangular
!> (Mmax=Lmax) throughout this codebase.
USE KINDS,   ONLY: dp, i4
USE GLOBALS, ONLY: pi
USE VSH,     ONLY: SHGLQ
IMPLICIT NONE
PRIVATE
PUBLIC :: ANGULAR_GRID_T, BUILD_ANGULAR_GRID

TYPE :: ANGULAR_GRID_T
  INTEGER(i4) :: LMAX      = 0
  INTEGER(i4) :: NLM       = 0
  INTEGER(i4) :: LMAX_QUAD = 0
  INTEGER(i4) :: NTHETA    = 0
  INTEGER(i4) :: NPHI      = 0
  LOGICAL     :: DEALIAS   = .FALSE.
  REAL(KIND=dp), ALLOCATABLE :: COSTHETA(:), THETA(:), WTHETA(:)
  REAL(KIND=dp), ALLOCATABLE :: PHI(:)
END TYPE ANGULAR_GRID_T

CONTAINS

!> @param AGRID Output grid.
!> @param LMAX Spherical-harmonic truncation degree.
!> @param DEALIAS Optional, default .FALSE.. When .TRUE., sizes the
!>   quadrature for the 3/2 rule (exact integration of quadratic products
!>   of Lmax-truncated fields) rather than the minimal exact sizing for a
!>   purely linear right-hand side.
SUBROUTINE BUILD_ANGULAR_GRID(AGRID, LMAX, DEALIAS)
  TYPE(ANGULAR_GRID_T), INTENT(OUT) :: AGRID
  INTEGER(KIND=i4),     INTENT(IN)  :: LMAX
  LOGICAL, INTENT(IN), OPTIONAL     :: DEALIAS
  INTEGER(KIND=i4) :: IQ, IP
  REAL(KIND=dp)    :: DPHI

  AGRID%LMAX    = LMAX
  AGRID%NLM     = (LMAX+1)**2
  AGRID%DEALIAS = .FALSE.
  IF (PRESENT(DEALIAS)) AGRID%DEALIAS = DEALIAS

  IF (AGRID%DEALIAS) THEN
    AGRID%LMAX_QUAD = CEILING(1.5_dp*LMAX) + 2
    AGRID%NPHI      = 3*LMAX + 2
  ELSE
    AGRID%LMAX_QUAD = LMAX
    AGRID%NPHI      = 2*LMAX + 1
  END IF
  AGRID%NTHETA = AGRID%LMAX_QUAD + 1

  ALLOCATE(AGRID%COSTHETA(AGRID%NTHETA))
  ALLOCATE(AGRID%THETA(AGRID%NTHETA))
  ALLOCATE(AGRID%WTHETA(AGRID%NTHETA))
  CALL SHGLQ(AGRID%COSTHETA, AGRID%WTHETA, AGRID%LMAX_QUAD)
  DO IQ = 1, AGRID%NTHETA
    AGRID%THETA(IQ) = ACOS(AGRID%COSTHETA(IQ))
  END DO

  ALLOCATE(AGRID%PHI(AGRID%NPHI))
  DPHI = 2.0_dp*pi/REAL(AGRID%NPHI, KIND=dp)
  DO IP = 1, AGRID%NPHI
    AGRID%PHI(IP) = REAL(IP-1, KIND=dp)*DPHI
  END DO
END SUBROUTINE BUILD_ANGULAR_GRID

END MODULE GRID_ANGULAR
