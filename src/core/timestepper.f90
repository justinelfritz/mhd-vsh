MODULE TIMESTEPPER
!> Shared fixed-step time-integration loop, written once against
!> REGIME_INTERFACE's REGIME_ADVANCE_I and reused by every regime's own
!> driver program: repeatedly calls the regime-supplied ADVANCE
!> procedure, advancing T by DT each step.
!>
!> Deliberately minimal -- nothing here about output, diagnostics, or
!> adaptive stepping yet; add those once there's an IO layer and a
!> second regime to design against (see REGIME_INTERFACE's header for
!> why STATE is CLASS(*)).
USE KINDS,            ONLY: dp, i4
USE REGIME_INTERFACE, ONLY: REGIME_ADVANCE_I
IMPLICIT NONE
PRIVATE
PUBLIC :: RUN, REGIME_ON_STEP_I

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
END INTERFACE

CONTAINS

!> @param ADVANCE The regime's own per-step procedure, matching
!>   REGIME_ADVANCE_I.
!> @param STATE The regime's own state, in/out, advanced in place over
!>   N_STEPS calls to ADVANCE.
!> @param DT Fixed step size.
!> @param N_STEPS Number of steps to take.
!> @param T_START Optional start time (default 0). ADVANCE's k-th call
!>   (k=1..N_STEPS) sees T = T_START + (k-1)*DT.
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

END MODULE TIMESTEPPER
