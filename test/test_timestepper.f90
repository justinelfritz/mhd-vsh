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
!>
!> TOY_COMPUTE_DT/TOY_SET_DT (added 2026-08-21) exercise
!> TIMESTEPPER::RUN_ADAPTIVE's generic mechanics with a toy "CFL"
!> condition unrelated to any real stability physics: DT_MAX=TOY_CFL_C/Y
!> (an arbitrary but deterministic, hand-checkable function of state),
!> and TOY_SET_DT simply records the DT it was given in TOY_LAST_SET_DT
!> (module SAVE state, standing in for a real regime's cached
!> DT-dependent factorization) -- enough to confirm RUN_ADAPTIVE calls
!> COMPUTE_DT/SET_DT at the right cadence and actually uses the DT they
!> produce, without needing a second real regime just for this test.
USE KINDS, ONLY: dp
IMPLICIT NONE
PRIVATE
PUBLIC :: TOY_STATE_T, TOY_ADVANCE, TOY_K, TOY_COMPUTE_DT, TOY_SET_DT, &
          TOY_CFL_C, TOY_LAST_SET_DT, TOY_N_SET_DT_CALLS

REAL(KIND=dp), PARAMETER :: TOY_K = 0.3_dp
REAL(KIND=dp), PARAMETER :: TOY_CFL_C = 0.1_dp

REAL(KIND=dp)    :: TOY_LAST_SET_DT    = -1.0_dp
INTEGER          :: TOY_N_SET_DT_CALLS = 0

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

!> Matches TIMESTEPPER::REGIME_DT_I. Toy "CFL" condition, arbitrary but
!> deterministic and hand-checkable: DT_MAX = TOY_CFL_C/Y.
FUNCTION TOY_COMPUTE_DT(STATE) RESULT(DT_MAX)
  CLASS(*), INTENT(IN) :: STATE
  REAL(KIND=dp) :: DT_MAX
  DT_MAX = 0.0_dp
  SELECT TYPE (STATE)
  TYPE IS (TOY_STATE_T)
    DT_MAX = TOY_CFL_C / STATE%Y
  END SELECT
END FUNCTION TOY_COMPUTE_DT

!> Matches TIMESTEPPER::REGIME_SET_DT_I. Records what it was given
!> rather than doing anything with it (a real regime would refresh its
!> own cached DT-dependent state here -- see HALL_REGIME::HALL_SET_DT).
SUBROUTINE TOY_SET_DT(DT)
  REAL(KIND=dp), INTENT(IN) :: DT
  TOY_LAST_SET_DT = DT
  TOY_N_SET_DT_CALLS = TOY_N_SET_DT_CALLS + 1
END SUBROUTINE TOY_SET_DT

END MODULE TEST_TOY_REGIME

PROGRAM TEST_TIMESTEPPER
USE KINDS,           ONLY: dp, i4
USE TIMESTEPPER,     ONLY: RUN, RUN_ADAPTIVE
USE TEST_TOY_REGIME, ONLY: TOY_STATE_T, TOY_ADVANCE, TOY_K, TOY_COMPUTE_DT, TOY_SET_DT, &
                            TOY_CFL_C, TOY_LAST_SET_DT, TOY_N_SET_DT_CALLS
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

CALL TEST_RUN_ADAPTIVE(N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_timestepper"
END IF

CONTAINS

!> Exercises RUN_ADAPTIVE's generic mechanics (recompute cadence, DT
!> actually used, SET_DT call count) against TEST_TOY_REGIME's toy
!> DT_MAX=TOY_CFL_C/Y "CFL" condition -- an independently-coded replica
!> of the same recurrence (recompute every DT_RECOMPUTE_EVERY steps,
!> else reuse the last DT; explicit-Euler decay in between) is
!> hand-rolled below rather than calling RUN_ADAPTIVE twice, so this
!> isn't just checking RUN_ADAPTIVE against itself.
SUBROUTINE TEST_RUN_ADAPTIVE(N_FAIL)
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  INTEGER(KIND=i4), PARAMETER :: N_STEPS_A = 12
  INTEGER(KIND=i4), PARAMETER :: EVERY     = 5
  REAL(KIND=dp),    PARAMETER :: SAFETY    = 0.5_dp
  TYPE(TOY_STATE_T) :: STATE_A
  REAL(KIND=dp)     :: Y_EXP, T_EXP, T_LAST_EXP, DT_EXP, ERR
  INTEGER(KIND=i4)  :: ISTEP, N_SET_DT_EXP

  STATE_A%Y = Y0
  CALL RUN_ADAPTIVE(TOY_ADVANCE, STATE_A, TOY_COMPUTE_DT, TOY_SET_DT, N_STEPS_A, &
    DT_RECOMPUTE_EVERY=EVERY, CFL_SAFETY=SAFETY, T_START=0.0_dp)

  ! Independent replica of the recurrence. T_LAST_EXP takes T_EXP's
  ! value BEFORE each step's increment (the k-th ADVANCE call sees
  ! T=T_START+sum of the first k-1 DTs, same convention RUN itself
  ! uses -- see EXPECTED_T_LAST above), so after the loop it holds the
  ! T seen by the LAST call, matching TOY_ADVANCE's own T_LAST bookkeeping.
  Y_EXP = Y0
  T_EXP = 0.0_dp
  T_LAST_EXP = 0.0_dp
  DT_EXP = 0.0_dp
  N_SET_DT_EXP = 0
  DO ISTEP = 1, N_STEPS_A
    IF (MOD(ISTEP-1, EVERY) == 0) THEN
      DT_EXP = SAFETY * (TOY_CFL_C / Y_EXP)
      N_SET_DT_EXP = N_SET_DT_EXP + 1
    END IF
    T_LAST_EXP = T_EXP
    Y_EXP = Y_EXP * (1.0_dp - TOY_K*DT_EXP)
    T_EXP = T_EXP + DT_EXP
  END DO

  ERR = ABS(STATE_A%Y - Y_EXP)
  CALL REPORT("run_adaptive_value", ERR, TOL, N_FAIL)

  ERR = ABS(STATE_A%T_LAST - T_LAST_EXP)
  CALL REPORT("run_adaptive_time_bookkeeping", ERR, TOL, N_FAIL)

  ERR = ABS(TOY_LAST_SET_DT - DT_EXP)
  CALL REPORT("run_adaptive_last_set_dt", ERR, TOL, N_FAIL)

  IF (TOY_N_SET_DT_CALLS /= N_SET_DT_EXP) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,I0,A,I0)') "FAIL  run_adaptive_set_dt_call_count  got=", &
      TOY_N_SET_DT_CALLS, "  expected=", N_SET_DT_EXP
  ELSE
    WRITE(*,'(A,I0)') "PASS  run_adaptive_set_dt_call_count  count=", TOY_N_SET_DT_CALLS
  END IF
END SUBROUTINE TEST_RUN_ADAPTIVE

SUBROUTINE REPORT(NAME, ERR, TOL_ARG, N_FAIL)
  CHARACTER(*),     INTENT(IN)    :: NAME
  REAL(KIND=dp),    INTENT(IN)    :: ERR, TOL_ARG
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  IF (ERR > TOL_ARG) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,ES10.3)') "FAIL  ", NAME, "  err=", ERR
  ELSE
    WRITE(*,'(A,A,A,ES10.3)') "PASS  ", NAME, "  err=", ERR
  END IF
END SUBROUTINE REPORT

END PROGRAM TEST_TIMESTEPPER
