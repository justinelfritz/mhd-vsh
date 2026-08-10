!> Checks LINEAR_SOLVE's FACTORIZE_DENSE/SOLVE_FACTORED against a known
!> solution: build a well-conditioned tridiagonal A and a chosen X_TRUE,
!> form B=A*X_TRUE, then recover X_TRUE from (A,B) alone. Also checks
!> that a single factorization can be reused for a second, independent
!> right-hand side (the whole point of separating FACTORIZE from SOLVE).
PROGRAM TEST_LINEAR_SOLVE
USE KINDS,        ONLY: dp, i4
USE LINEAR_SOLVE, ONLY: DENSE_FACTORS_T, FACTORIZE_DENSE, SOLVE_FACTORED
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N   = 5
REAL(KIND=dp),    PARAMETER :: TOL = 1.0E-11_dp

REAL(KIND=dp) :: A(N,N)
TYPE(DENSE_FACTORS_T) :: FACTORS
INTEGER(KIND=i4) :: I, N_FAIL

A = 0.0_dp
DO I = 1, N
  A(I,I) = 4.0_dp
  IF (I > 1) A(I,I-1) = 1.0_dp
  IF (I < N) A(I,I+1) = 1.0_dp
END DO

CALL FACTORIZE_DENSE(FACTORS, A)

N_FAIL = 0
IF (FACTORS%INFO /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,I0)') "FAIL  factorize_info  info=", FACTORS%INFO
ELSE
  WRITE(*,'(A)') "PASS  factorize_info"
END IF

CALL RUN_MULTI_RHS_CHECK(A, FACTORS, N_FAIL)
CALL RUN_REUSE_CHECK(A, FACTORS, N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_linear_solve"
END IF

CONTAINS

!> Part A: a single SOLVE_FACTORED call with two right-hand-side columns
!> at once must recover both columns of X_TRUE.
SUBROUTINE RUN_MULTI_RHS_CHECK(A, FACTORS, N_FAIL)
  REAL(KIND=dp),          INTENT(IN)    :: A(N,N)
  TYPE(DENSE_FACTORS_T),  INTENT(IN)    :: FACTORS
  INTEGER(KIND=i4),       INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp) :: X_TRUE(N,2), B(N,2), X(N,2)
  REAL(KIND=dp) :: MAX_ERR
  INTEGER(KIND=i4) :: J

  DO J = 1, N
    X_TRUE(J,1) = REAL(J, KIND=dp)
    X_TRUE(J,2) = REAL(N-J+1, KIND=dp)
  END DO
  B = MATMUL(A, X_TRUE)

  CALL SOLVE_FACTORED(X, FACTORS, B)
  MAX_ERR = MAXVAL(ABS(X-X_TRUE))
  IF (MAX_ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  multi_rhs  max_err=", MAX_ERR
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  multi_rhs  max_err=", MAX_ERR
  END IF
END SUBROUTINE RUN_MULTI_RHS_CHECK

!> Part B: the same FACTORS, solved a second time against an unrelated
!> right-hand side, must still recover the correct answer -- confirms
!> SOLVE_FACTORED doesn't mutate/consume the stored factorization.
SUBROUTINE RUN_REUSE_CHECK(A, FACTORS, N_FAIL)
  REAL(KIND=dp),          INTENT(IN)    :: A(N,N)
  TYPE(DENSE_FACTORS_T),  INTENT(IN)    :: FACTORS
  INTEGER(KIND=i4),       INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp) :: X_TRUE(N,1), B(N,1), X(N,1)
  REAL(KIND=dp) :: MAX_ERR
  INTEGER(KIND=i4) :: J

  DO J = 1, N
    X_TRUE(J,1) = SIN(REAL(J, KIND=dp))
  END DO
  B = MATMUL(A, X_TRUE)

  CALL SOLVE_FACTORED(X, FACTORS, B)
  MAX_ERR = MAXVAL(ABS(X-X_TRUE))
  IF (MAX_ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  reuse  max_err=", MAX_ERR
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  reuse  max_err=", MAX_ERR
  END IF
END SUBROUTINE RUN_REUSE_CHECK

END PROGRAM TEST_LINEAR_SOLVE
