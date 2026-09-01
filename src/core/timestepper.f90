MODULE TIMESTEPPER
!> Shared fixed-step time-integration loop, written once against
!> REGIME_INTERFACE's REGIME_ADVANCE_I and reused by every regime's own
!> driver program: repeatedly calls the regime-supplied ADVANCE
!> procedure, advancing T by DT each step.
!>
!> RUN is fixed-step; RUN_ADAPTIVE (added 2026-08-21) recomputes DT
!> periodically from the regime's own state via a driver-supplied
!> callback (e.g. FIELD_DIAGNOSTICS::HALL_COURANT_TIMESTEP) -- see its
!> own docstring. Both still say nothing about output/IO beyond the
!> ON_STEP observer.
USE KINDS,            ONLY: dp, i4
USE REGIME_INTERFACE, ONLY: REGIME_ADVANCE_I
IMPLICIT NONE
PRIVATE
PUBLIC :: RUN, REGIME_ON_STEP_I, RUN_ADAPTIVE, REGIME_DT_I, REGIME_SET_DT_I

ABSTRACT INTERFACE
  !> Optional per-step observer, supplied by the driver program (not by
  !> the regime -- unlike REGIME_ADVANCE_I, which the regime itself
  !> implements). Called read-only after every ADVANCE, e.g. to log a
  !> diagnostic time series; must not mutate STATE.
  !>
  !> @param STATE The regime's state, read-only, as it stands after this
  !>   step.
  !> @param T Time after this step.
  !> @param ISTEP Step number just completed, 1..N_STEPS.
  SUBROUTINE REGIME_ON_STEP_I(STATE, T, ISTEP)
    IMPORT :: dp, i4
    CLASS(*),         INTENT(IN) :: STATE
    REAL(KIND=dp),    INTENT(IN) :: T
    INTEGER(KIND=i4), INTENT(IN) :: ISTEP
  END SUBROUTINE REGIME_ON_STEP_I

  !> Recomputes the maximum stable step size from the regime's current
  !> state (e.g. a CFL/Courant condition) -- supplied by the driver as a
  !> thin wrapper around whatever diagnostic function the regime's own
  !> physics implies; RUN_ADAPTIVE itself has no notion of what makes a
  !> step "stable" for any particular regime.
  !>
  !> @param STATE The regime's state, read-only, as it stands when this
  !>   is called (may be mid-run, not just at t=0).
  !> Returns: the raw (safety-factor-free) maximum stable step size --
  !>   RUN_ADAPTIVE itself applies CFL_SAFETY on top of this.
  FUNCTION REGIME_DT_I(STATE) RESULT(DT_MAX)
    IMPORT :: dp
    CLASS(*), INTENT(IN) :: STATE
    REAL(KIND=dp) :: DT_MAX
  END FUNCTION REGIME_DT_I

  !> Re-initializes whatever cached, DT-dependent state the regime's own
  !> ADVANCE relies on (e.g. an implicit solve's cached LU factorization)
  !> for a new step size -- supplied by the driver as a thin wrapper
  !> around the regime's own `*_INIT`. Needed because REGIME_ADVANCE_I's
  !> fixed (STATE,DT,T) signature has no room to pass RGRID/OPS/ETA/...
  !> required to rebuild such caches, so every regime instead caches
  !> those itself (module SAVE state) at `*_INIT` time -- the same
  !> reasoning that motivated that design already applies here: this
  !> callback re-invokes `*_INIT` using the regime's own cached values,
  !> with only DT actually changing.
  !>
  !> @param DT The new step size to adopt.
  SUBROUTINE REGIME_SET_DT_I(DT)
    IMPORT :: dp
    REAL(KIND=dp), INTENT(IN) :: DT
  END SUBROUTINE REGIME_SET_DT_I
END INTERFACE

CONTAINS

!> @param ADVANCE The regime's own per-step procedure, matching
!>   REGIME_ADVANCE_I.
!> @param STATE The regime's own state, in/out, advanced in place over
!>   N_STEPS calls to ADVANCE.
!> @param DT Fixed step size.
!> @param N_STEPS Number of steps to take.
!> @param T_START Optional start time (default 0). ADVANCE's k-th call
!>   (k=1..N_STEPS) sees `T = T_START + (k-1)*DT`.
!> @param ON_STEP Optional observer matching REGIME_ON_STEP_I, called
!>   after every step (e.g. to log a diagnostic time series); has no
!>   effect on the evolution itself.
SUBROUTINE RUN(ADVANCE, STATE, DT, N_STEPS, T_START, ON_STEP)
  PROCEDURE(REGIME_ADVANCE_I)             :: ADVANCE
  CLASS(*),         INTENT(INOUT)         :: STATE
  REAL(KIND=dp),    INTENT(IN)            :: DT
  INTEGER(KIND=i4), INTENT(IN)            :: N_STEPS
  REAL(KIND=dp),    INTENT(IN), OPTIONAL  :: T_START
  PROCEDURE(REGIME_ON_STEP_I), OPTIONAL   :: ON_STEP
  REAL(KIND=dp)    :: T
  INTEGER(KIND=i4) :: ISTEP

  T = 0.0_dp
  IF (PRESENT(T_START)) T = T_START

  DO ISTEP = 1, N_STEPS
    CALL ADVANCE(STATE, DT, T)
    T = T + DT
    IF (PRESENT(ON_STEP)) CALL ON_STEP(STATE, T, ISTEP)
  END DO
END SUBROUTINE RUN

!> Like RUN, but DT is recomputed from the regime's own state every
!> DT_RECOMPUTE_EVERY steps (via COMPUTE_DT) rather than held fixed for
!> the whole run, scaled by CFL_SAFETY<1 for a safety margin below the
!> raw stable limit COMPUTE_DT returns. Whenever DT is recomputed,
!> SET_DT is called so the regime can refresh any cached, DT-dependent
!> state (e.g. an implicit solve's LU factorization) before the next
!> ADVANCE call -- SET_DT is called even on the very first step, so the
!> regime never runs with whatever DT its own prior `*_INIT` call happened
!> to be given.
!>
!> DT_RECOMPUTE_EVERY is deliberately NOT defaulted to 1 (recompute-every-
!> step): recomputing DT is cheap (COMPUTE_DT is typically an O(N_r)
!> diagnostic), but SET_DT is typically NOT (e.g. HALL_REGIME's own
!> SET_DT re-factorizes an `O(N_r**3)`-per-l dense system) -- callers must
!> choose a cadence that amortizes that cost against how fast their
!> regime's own stability limit actually drifts, there's no
!> regime-agnostic "safe default" here.
!>
!> @param ADVANCE The regime's own per-step procedure, matching
!>   REGIME_ADVANCE_I.
!> @param STATE The regime's own state, in/out, advanced in place over
!>   N_STEPS calls to ADVANCE.
!> @param COMPUTE_DT Matching REGIME_DT_I -- recomputes the raw maximum
!>   stable step size from STATE.
!> @param SET_DT Matching REGIME_SET_DT_I -- re-initializes the regime's
!>   own cached DT-dependent state for a newly chosen step size.
!> @param N_STEPS Number of steps to take.
!> @param DT_RECOMPUTE_EVERY Recompute (and call SET_DT) every this many
!>   steps -- 1 means every step.
!> @param CFL_SAFETY Optional safety margin, `DT = CFL_SAFETY*COMPUTE_DT(STATE)`
!>   (default 0.5).
!> @param T_START Optional start time (default 0).
!> @param T_MAX Optional target simulated time (added 2026-08-24, for a
!>   real `tmax=1000` yr production run where the adaptive `dt`
!>   trajectory can't be predicted well enough in advance to pick
!>   N_STEPS precisely). When present, the loop also stops as soon as
!>   `T` reaches `T_MAX`, clipping the final step's `DT` so `T` lands
!>   exactly on `T_MAX` rather than overshooting it -- N_STEPS still
!>   applies as an upper safety cap (in case `DT` turns out smaller than
!>   expected and `T_MAX` would otherwise take far longer than intended
!>   to reach). Absent, behavior is exactly the original N_STEPS-only
!>   loop, unchanged code path.
!> @param ON_STEP Optional observer matching REGIME_ON_STEP_I, called
!>   after every step; has no effect on the evolution itself.
SUBROUTINE RUN_ADAPTIVE(ADVANCE, STATE, COMPUTE_DT, SET_DT, N_STEPS, DT_RECOMPUTE_EVERY, &
    CFL_SAFETY, T_START, T_MAX, ON_STEP)
  PROCEDURE(REGIME_ADVANCE_I)             :: ADVANCE
  CLASS(*),         INTENT(INOUT)         :: STATE
  PROCEDURE(REGIME_DT_I)                  :: COMPUTE_DT
  PROCEDURE(REGIME_SET_DT_I)              :: SET_DT
  INTEGER(KIND=i4), INTENT(IN)            :: N_STEPS
  INTEGER(KIND=i4), INTENT(IN)            :: DT_RECOMPUTE_EVERY
  REAL(KIND=dp),    INTENT(IN), OPTIONAL  :: CFL_SAFETY
  REAL(KIND=dp),    INTENT(IN), OPTIONAL  :: T_START
  REAL(KIND=dp),    INTENT(IN), OPTIONAL  :: T_MAX
  PROCEDURE(REGIME_ON_STEP_I), OPTIONAL   :: ON_STEP
  REAL(KIND=dp)    :: T, DT, SAFETY
  INTEGER(KIND=i4) :: ISTEP

  SAFETY = 0.5_dp
  IF (PRESENT(CFL_SAFETY)) SAFETY = CFL_SAFETY
  T = 0.0_dp
  IF (PRESENT(T_START)) T = T_START

  DO ISTEP = 1, N_STEPS
    IF (MOD(ISTEP-1, DT_RECOMPUTE_EVERY) == 0) THEN
      DT = SAFETY * COMPUTE_DT(STATE)
      CALL SET_DT(DT)
    END IF
    IF (PRESENT(T_MAX)) THEN
      IF (T+DT > T_MAX) THEN
        DT = T_MAX - T
        CALL SET_DT(DT)
      END IF
    END IF
    CALL ADVANCE(STATE, DT, T)
    T = T + DT
    IF (PRESENT(ON_STEP)) CALL ON_STEP(STATE, T, ISTEP)
    IF (PRESENT(T_MAX)) THEN
      IF (T >= T_MAX) EXIT
    END IF
  END DO
END SUBROUTINE RUN_ADAPTIVE

END MODULE TIMESTEPPER
