!> Validates the REGIME_INTERFACE/TIMESTEPPER plumbing end-to-end using
!> a toy regime (explicit-Euler exponential decay, dY/dt=-K*Y) defined
!> entirely within this test -- no real regime exists yet, so this
!> exercises exactly the two pieces this milestone adds: that a
!> concrete ADVANCE procedure matching REGIME_ADVANCE_I's CLASS(*) state
!> can be correctly unpacked via SELECT TYPE, and that RUN calls it the
!> right number of times with the right T at each call.
MODULE TEST_TOY_REGIME
!> A minimal regime for testing: state is a single real value Y evolved
!> by explicit Euler for dY/dt=-K*Y, i.e. Y <- Y*(1-K*DT) each step.
!> Also records the T seen on its most recent call, to check
!> TIMESTEPPER::RUN's bookkeeping independently of the toy physics.
USE KINDS, ONLY: dp
IMPLICIT NONE
PRIVATE
PUBLIC :: TOY_STATE_T, TOY_ADVANCE, TOY_K

REAL(KIND=dp), PARAMETER :: TOY_K = 0.3_dp

TYPE :: TOY_STATE_T
  REAL(KIND=dp) :: Y      = 0.0_dp
  REAL(KIND=dp) :: T_LAST = -1.0_dp
END TYPE TOY_STATE_T

CONTAINS

SUBROUTINE TOY_ADVANCE(STATE, DT, T)
  CLASS(*),      INTENT(INOUT) :: STATE
  REAL(KIND=dp), INTENT(IN)    :: DT, T
  SELECT TYPE (STATE)
  TYPE IS (TOY_STATE_T)
    STATE%Y      = STATE%Y * (1.0_dp - TOY_K*DT)
    STATE%T_LAST = T
  END SELECT
END SUBROUTINE TOY_ADVANCE

END MODULE TEST_TOY_REGIME

PROGRAM TEST_TIMESTEPPER
USE KINDS,           ONLY: dp, i4
USE TIMESTEPPER,     ONLY: RUN
USE TEST_TOY_REGIME, ONLY: TOY_STATE_T, TOY_ADVANCE, TOY_K
IMPLICIT NONE

REAL(KIND=dp),    PARAMETER :: Y0      = 2.0_dp
REAL(KIND=dp),    PARAMETER :: DT      = 0.01_dp
INTEGER(KIND=i4), PARAMETER :: N_STEPS = 50
REAL(KIND=dp),    PARAMETER :: TOL     = 1.0E-12_dp

TYPE(TOY_STATE_T) :: STATE
REAL(KIND=dp) :: EXPECTED_Y, EXPECTED_T_LAST, ERR_Y, ERR_T
INTEGER(KIND=i4) :: N_FAIL

STATE%Y = Y0
CALL RUN(TOY_ADVANCE, STATE, DT, N_STEPS)

! Explicit Euler is exactly this recurrence -- no discretization-error
! tolerance to guess, the comparison is exact up to roundoff.
EXPECTED_Y = Y0 * (1.0_dp - TOY_K*DT)**N_STEPS
EXPECTED_T_LAST = REAL(N_STEPS-1, KIND=dp) * DT

ERR_Y = ABS(STATE%Y - EXPECTED_Y)
ERR_T = ABS(STATE%T_LAST - EXPECTED_T_LAST)

N_FAIL = 0
IF (ERR_Y > TOL) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,ES10.3)') "FAIL  toy_decay_value  err=", ERR_Y
ELSE
  WRITE(*,'(A,ES10.3)') "PASS  toy_decay_value  err=", ERR_Y
END IF

IF (ERR_T > TOL) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,ES10.3)') "FAIL  toy_decay_time_bookkeeping  err=", ERR_T
ELSE
  WRITE(*,'(A,ES10.3)') "PASS  toy_decay_time_bookkeeping  err=", ERR_T
END IF

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_timestepper"
END IF

END PROGRAM TEST_TIMESTEPPER
