MODULE GRID_RADIAL
!> Radial grid: finite-difference nodes, either a spherical shell
!> (r_min..r_max, genuine boundary conditions at both ends) or a full
!> sphere (r=0..r_max, half-step-staggered so no node sits exactly at the
!> origin, where l(l+1)/r**2 is singular for l>=1).
USE KINDS, ONLY: dp, i4
IMPLICIT NONE
PRIVATE
PUBLIC :: RADIAL_GRID_T, BUILD_RADIAL_GRID, RADIAL_UNIFORM

INTEGER(KIND=i4), PARAMETER :: RADIAL_UNIFORM = 0

TYPE :: RADIAL_GRID_T
  INTEGER(KIND=i4) :: N            = 0
  LOGICAL          :: FULL_SPHERE  = .FALSE.
  REAL(KIND=dp)    :: R_MIN        = 0.0_dp
  REAL(KIND=dp)    :: R_MAX        = 0.0_dp
  INTEGER(KIND=i4) :: STRETCH_TYPE = RADIAL_UNIFORM
  REAL(KIND=dp)    :: STRETCH_PARAM = 0.0_dp
  INTEGER(KIND=i4) :: FD_ORDER     = 4
  REAL(KIND=dp), ALLOCATABLE :: R(:)
END TYPE RADIAL_GRID_T

CONTAINS

!> @param RGRID Output grid.
!> @param N Number of radial nodes.
!> @param R_MIN Inner boundary (ignored, forced to 0, when FULL_SPHERE).
!> @param R_MAX Outer boundary.
!> @param FULL_SPHERE .TRUE. for a full sphere (half-step-staggered grid,
!>   origin regularity handled by RADIAL_OPERATORS), .FALSE. for a shell
!>   (uniform grid spanning [R_MIN,R_MAX] with genuine BCs at both ends).
!> @param FD_ORDER Optional nominal interior finite-difference order,
!>   default 4.
SUBROUTINE BUILD_RADIAL_GRID(RGRID, N, R_MIN, R_MAX, FULL_SPHERE, FD_ORDER)
  TYPE(RADIAL_GRID_T), INTENT(OUT) :: RGRID
  INTEGER(KIND=i4),    INTENT(IN)  :: N
  REAL(KIND=dp),       INTENT(IN)  :: R_MIN, R_MAX
  LOGICAL,              INTENT(IN)  :: FULL_SPHERE
  INTEGER(KIND=i4), INTENT(IN), OPTIONAL :: FD_ORDER
  INTEGER(KIND=i4) :: I
  REAL(KIND=dp) :: DR

  RGRID%N            = N
  RGRID%FULL_SPHERE   = FULL_SPHERE
  RGRID%STRETCH_TYPE  = RADIAL_UNIFORM
  RGRID%STRETCH_PARAM = 0.0_dp
  RGRID%FD_ORDER      = 4
  IF (PRESENT(FD_ORDER)) RGRID%FD_ORDER = FD_ORDER

  ALLOCATE(RGRID%R(N))

  IF (FULL_SPHERE) THEN
    RGRID%R_MIN = 0.0_dp
    RGRID%R_MAX = R_MAX
    DR = R_MAX / (REAL(N, KIND=dp) - 0.5_dp)
    DO I = 1, N
      RGRID%R(I) = (REAL(I, KIND=dp) - 0.5_dp) * DR
    END DO
  ELSE
    RGRID%R_MIN = R_MIN
    RGRID%R_MAX = R_MAX
    DR = (R_MAX - R_MIN) / REAL(N-1, KIND=dp)
    DO I = 1, N
      RGRID%R(I) = R_MIN + REAL(I-1, KIND=dp)*DR
    END DO
  END IF
END SUBROUTINE BUILD_RADIAL_GRID

END MODULE GRID_RADIAL
