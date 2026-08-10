!> Checks the vacuum outer-boundary row mutations in BOUNDARY_CONDITIONS:
!> APPLY_VACUUM_BC_POLOIDAL (Robin row on Phi) and APPLY_VACUUM_BC_TOROIDAL
!> (Dirichlet row on Psi). Part A/B are direct structural checks (only row
!> N changes, by exactly the intended amount). Part C is a physics check:
!> Phi(r)=r**(-l) satisfies d(Phi)/dr + (l/r)*Phi = 0 at every r, not just
!> r_out, so applying the augmented row to it should return ~0 (up to the
!> D1 stencil's own FD truncation error) regardless of l.
PROGRAM TEST_BOUNDARY_CONDITIONS
USE KINDS,               ONLY: dp, i4
USE GRID_RADIAL,         ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,    ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE BOUNDARY_CONDITIONS, ONLY: APPLY_VACUUM_BC_POLOIDAL, APPLY_VACUUM_BC_TOROIDAL
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R   = 40
REAL(KIND=dp),    PARAMETER :: R_MIN = 0.5_dp
REAL(KIND=dp),    PARAMETER :: R_MAX = 1.0_dp

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
INTEGER(KIND=i4) :: N_FAIL

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)

N_FAIL = 0
CALL RUN_POLOIDAL_ROW_CHECK(N_FAIL)
CALL RUN_TOROIDAL_ROW_CHECK(N_FAIL)
CALL RUN_POLOIDAL_PHYSICS_CHECK(N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_boundary_conditions"
END IF

CONTAINS

!> Part A: APPLY_VACUUM_BC_POLOIDAL must change only D1(N,N), by exactly
!> l/r_out, and leave every other entry (including the rest of row N)
!> untouched.
SUBROUTINE RUN_POLOIDAL_ROW_CHECK(N_FAIL)
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  INTEGER(KIND=i4), PARAMETER :: L = 3
  REAL(KIND=dp), ALLOCATABLE :: D1_TEST(:,:)
  REAL(KIND=dp) :: EXPECTED_DIAG, MAX_OFFDIAG_ERR
  INTEGER(KIND=i4) :: J

  ALLOCATE(D1_TEST(N_R,N_R))
  D1_TEST = OPS%D1
  CALL APPLY_VACUUM_BC_POLOIDAL(D1_TEST, L, RGRID)

  EXPECTED_DIAG = OPS%D1(N_R,N_R) + REAL(L,KIND=dp)/RGRID%R(N_R)
  MAX_OFFDIAG_ERR = 0.0_dp
  DO J = 1, N_R
    IF (J /= N_R) THEN
      MAX_OFFDIAG_ERR = MAX(MAX_OFFDIAG_ERR, ABS(D1_TEST(N_R,J)-OPS%D1(N_R,J)))
    END IF
  END DO
  MAX_OFFDIAG_ERR = MAX(MAX_OFFDIAG_ERR, MAXVAL(ABS(D1_TEST(1:N_R-1,:)-OPS%D1(1:N_R-1,:))))

  IF (ABS(D1_TEST(N_R,N_R)-EXPECTED_DIAG) > 1.0E-14_dp .OR. MAX_OFFDIAG_ERR > 1.0E-14_dp) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3,A,ES10.3)') "FAIL  poloidal_row  diag_err=", &
      ABS(D1_TEST(N_R,N_R)-EXPECTED_DIAG), "  offdiag_err=", MAX_OFFDIAG_ERR
  ELSE
    WRITE(*,'(A)') "PASS  poloidal_row"
  END IF
  DEALLOCATE(D1_TEST)
END SUBROUTINE RUN_POLOIDAL_ROW_CHECK

!> Part B: APPLY_VACUUM_BC_TOROIDAL must set row N to the identity row
!> and leave every other row untouched.
SUBROUTINE RUN_TOROIDAL_ROW_CHECK(N_FAIL)
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp), ALLOCATABLE :: OP_TEST(:,:)
  REAL(KIND=dp) :: ROW_ERR, REST_ERR

  ALLOCATE(OP_TEST(N_R,N_R))
  OP_TEST = OPS%D1
  CALL APPLY_VACUUM_BC_TOROIDAL(OP_TEST, RGRID)

  ROW_ERR  = MAXVAL(ABS(OP_TEST(N_R,1:N_R-1))) + ABS(OP_TEST(N_R,N_R)-1.0_dp)
  REST_ERR = MAXVAL(ABS(OP_TEST(1:N_R-1,:)-OPS%D1(1:N_R-1,:)))

  IF (ROW_ERR > 1.0E-14_dp .OR. REST_ERR > 1.0E-14_dp) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3,A,ES10.3)') "FAIL  toroidal_row  row_err=", ROW_ERR, &
      "  rest_err=", REST_ERR
  ELSE
    WRITE(*,'(A)') "PASS  toroidal_row"
  END IF
  DEALLOCATE(OP_TEST)
END SUBROUTINE RUN_TOROIDAL_ROW_CHECK

!> Part C: physics check. Phi(r)=r**(-l) exactly satisfies the vacuum
!> Robin condition at every radius (d(Phi)/dr = -l*r**(-l-1) =
!> -(l/r)*Phi identically), so dotting the augmented row N against this
!> profile sampled on the grid should return ~0, bounded only by D1's
!> own FD truncation error -- a check that the row assembly (not just
!> its diagonal entry) is correct, for two different l.
SUBROUTINE RUN_POLOIDAL_PHYSICS_CHECK(N_FAIL)
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp),    PARAMETER :: TOL = 1.0E-4_dp
  INTEGER(KIND=i4), PARAMETER :: N_CASES = 2
  INTEGER(KIND=i4), PARAMETER :: L_CASE(N_CASES) = (/1, 4/)
  REAL(KIND=dp), ALLOCATABLE :: D1_TEST(:,:), PHI(:)
  REAL(KIND=dp) :: RESIDUAL
  INTEGER(KIND=i4) :: K

  ALLOCATE(D1_TEST(N_R,N_R), PHI(N_R))
  DO K = 1, N_CASES
    D1_TEST = OPS%D1
    CALL APPLY_VACUUM_BC_POLOIDAL(D1_TEST, L_CASE(K), RGRID)
    PHI = RGRID%R ** (-L_CASE(K))
    RESIDUAL = DOT_PRODUCT(D1_TEST(N_R,:), PHI)
    IF (ABS(RESIDUAL) > TOL) THEN
      N_FAIL = N_FAIL + 1
      WRITE(*,'(A,I0,A,ES10.3)') "FAIL  poloidal_physics  l=", L_CASE(K), &
        "  residual=", RESIDUAL
    ELSE
      WRITE(*,'(A,I0,A,ES10.3)') "PASS  poloidal_physics  l=", L_CASE(K), &
        "  residual=", RESIDUAL
    END IF
  END DO
  DEALLOCATE(D1_TEST, PHI)
END SUBROUTINE RUN_POLOIDAL_PHYSICS_CHECK

END PROGRAM TEST_BOUNDARY_CONDITIONS
