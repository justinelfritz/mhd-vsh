!> Regression test for SPHERICAL_BESSEL_JN_YN. Two independent checks,
!> neither a self-consistency triviality:
!>  (1) closed-form cross-checks at n=0,1,2 (hand-derivable from
!>      elementary trig identities, catches sign/indexing bugs);
!>  (2) the Wronskian identity j_n(x)*y_{n-1}(x)-j_{n-1}(x)*y_n(x)=1/x**2,
!>      an EXACT mathematical identity for every n>=1 and every x>0 --
!>      verified by hand at n=1 while writing the design plan for this
!>      module (expands to cos**2(x)/x**2+sin**2(x)/x**2=1/x**2 exactly).
!> Both are swept across x=0.5..20 and n up to 20 to stress Miller's
!> algorithm's stability margin (N_START heuristic in JN_MILLER) at the
!> actual scale this project needs, not just at one convenient point.
!> Derivatives (DJ, DY) are checked separately against a central finite
!> difference, an independent check of the derivative identity itself.
PROGRAM TEST_SPHERICAL_BESSEL
USE KINDS,             ONLY: dp, i4
USE SPHERICAL_BESSEL,  ONLY: SPHERICAL_BESSEL_JN_YN
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: NMAX = 20
REAL(KIND=dp),    PARAMETER :: TOL = 1.0E-9_dp
REAL(KIND=dp),    PARAMETER :: TOL_FD = 1.0E-6_dp
REAL(KIND=dp),    PARAMETER :: H_FD = 1.0E-5_dp
REAL(KIND=dp),    PARAMETER :: X_VALUES(6) = (/0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 20.0_dp/)

REAL(KIND=dp) :: J(0:NMAX), Y(0:NMAX), DJ(0:NMAX), DY(0:NMAX)
REAL(KIND=dp) :: J_LO(0:NMAX), Y_LO(0:NMAX), J_HI(0:NMAX), Y_HI(0:NMAX)
REAL(KIND=dp) :: X, J0_EXACT, J1_EXACT, J2_EXACT, Y0_EXACT, Y1_EXACT, Y2_EXACT
REAL(KIND=dp) :: WRONSKIAN, FD_DERIV
INTEGER(KIND=i4) :: N_FAIL, IX, N
CHARACTER(LEN=64) :: LABEL

N_FAIL = 0

DO IX = 1, SIZE(X_VALUES)
  X = X_VALUES(IX)
  CALL SPHERICAL_BESSEL_JN_YN(X, NMAX, J, Y, DJ, DY)

  ! --- Closed-form cross-checks at n=0,1,2 ---
  J0_EXACT = SIN(X)/X
  J1_EXACT = SIN(X)/X**2 - COS(X)/X
  J2_EXACT = (3.0_dp/X**3 - 1.0_dp/X)*SIN(X) - (3.0_dp/X**2)*COS(X)
  Y0_EXACT = -COS(X)/X
  Y1_EXACT = -COS(X)/X**2 - SIN(X)/X
  Y2_EXACT = (-3.0_dp/X**3 + 1.0_dp/X)*COS(X) - (3.0_dp/X**2)*SIN(X)

  WRITE(LABEL,'(A,F5.2)') "j0_x=", X; CALL CHECK_CLOSE(LABEL, J(0), J0_EXACT, TOL, N_FAIL)
  WRITE(LABEL,'(A,F5.2)') "j1_x=", X; CALL CHECK_CLOSE(LABEL, J(1), J1_EXACT, TOL, N_FAIL)
  WRITE(LABEL,'(A,F5.2)') "j2_x=", X; CALL CHECK_CLOSE(LABEL, J(2), J2_EXACT, TOL, N_FAIL)
  WRITE(LABEL,'(A,F5.2)') "y0_x=", X; CALL CHECK_CLOSE(LABEL, Y(0), Y0_EXACT, TOL, N_FAIL)
  WRITE(LABEL,'(A,F5.2)') "y1_x=", X; CALL CHECK_CLOSE(LABEL, Y(1), Y1_EXACT, TOL, N_FAIL)
  WRITE(LABEL,'(A,F5.2)') "y2_x=", X; CALL CHECK_CLOSE(LABEL, Y(2), Y2_EXACT, TOL, N_FAIL)

  ! --- Wronskian identity, n=1..NMAX: exact for every n, every x>0 ---
  DO N = 1, NMAX
    WRONSKIAN = J(N)*Y(N-1) - J(N-1)*Y(N)
    WRITE(LABEL,'(A,I0,A,F5.2)') "wronskian_n=", N, "_x=", X
    CALL CHECK_CLOSE(LABEL, WRONSKIAN, 1.0_dp/X**2, TOL, N_FAIL)
  END DO

  ! --- Derivative check against central finite difference, n=0 and n=NMAX ---
  CALL SPHERICAL_BESSEL_JN_YN(X-H_FD, NMAX, J_LO, Y_LO)
  CALL SPHERICAL_BESSEL_JN_YN(X+H_FD, NMAX, J_HI, Y_HI)

  FD_DERIV = (J_HI(0)-J_LO(0))/(2.0_dp*H_FD)
  WRITE(LABEL,'(A,F5.2)') "dj0_fd_x=", X; CALL CHECK_CLOSE(LABEL, DJ(0), FD_DERIV, TOL_FD, N_FAIL)

  FD_DERIV = (Y_HI(0)-Y_LO(0))/(2.0_dp*H_FD)
  WRITE(LABEL,'(A,F5.2)') "dy0_fd_x=", X; CALL CHECK_CLOSE(LABEL, DY(0), FD_DERIV, TOL_FD, N_FAIL)

  FD_DERIV = (J_HI(NMAX)-J_LO(NMAX))/(2.0_dp*H_FD)
  WRITE(LABEL,'(A,F5.2)') "djN_fd_x=", X; CALL CHECK_CLOSE(LABEL, DJ(NMAX), FD_DERIV, TOL_FD, N_FAIL)

  FD_DERIV = (Y_HI(NMAX)-Y_LO(NMAX))/(2.0_dp*H_FD)
  WRITE(LABEL,'(A,F5.2)') "dyN_fd_x=", X; CALL CHECK_CLOSE(LABEL, DY(NMAX), FD_DERIV, TOL_FD, N_FAIL)
END DO

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_spherical_bessel"
END IF

CONTAINS

SUBROUTINE CHECK_CLOSE(NAME, VAL, EXPECTED, TOL, N_FAIL)
  CHARACTER(LEN=*), INTENT(IN)    :: NAME
  REAL(KIND=dp),    INTENT(IN)    :: VAL, EXPECTED, TOL
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp) :: ERR
  ERR = ABS(VAL-EXPECTED)/MAX(ABS(EXPECTED), 1.0E-300_dp)
  IF (ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,ES10.3)') "FAIL  ", NAME, "  relerr=", ERR
  ELSE
    WRITE(*,'(A,A,A,ES10.3)') "PASS  ", NAME, "  relerr=", ERR
  END IF
END SUBROUTINE CHECK_CLOSE

END PROGRAM TEST_SPHERICAL_BESSEL
