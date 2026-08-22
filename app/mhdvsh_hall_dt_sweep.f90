!> Combined resistive+Hall regime: dt-stability sweep. Per user
!> (2026-08-21), a pivot from "run to tmax=1000yr" (infeasible at the
!> validated dt_hall=0.001 -- ~1.5-2 days of compute) to first mapping
!> out the safe dt_hall range at a fixed, short tmax=10yr, using the
!> CURRENT parameterization (Lmax=30, NS-crust R_MIN/R_MAX=9/10 km,
!> ETA=1e-6, F_HALL=0.01, same single-mode Phi(1,0) seed as
!> app/mhdvsh_hall.f90).
!>
!> N_SUB is fixed at 1 here (unlike mhdvsh_hall.f90's default N_SUB=10)
!> so the command-line DT directly equals the Hall RK4 substep size --
!> this sweep is specifically about the Hall substep's own stability
!> limit, not about the outer implicit-diffusion cadence (which is
!> unconditionally stable regardless of its step size, backward Euler,
!> and at ETA=1e-6 barely matters over these timescales anyway -- see
!> ROADMAP.md 2026-08-21).
!>
!> Usage: mhdvsh_hall_dt_sweep <dt_yr> <tmax_yr> [output_data_path]
!> N_STEPS is computed as NINT(tmax_yr/dt_yr). Blow-up detection matches
!> app/mhdvsh_hall_stability_experiment.f90's convention: NaN or
!> `E > BLOWUP_FACTOR*E0` at any logged point.
PROGRAM MHDVSH_HALL_DT_SWEEP
USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,   ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS,  ONLY: TOTAL_MAGNETIC_ENERGY
USE DIFFUSION_REGIME,   ONLY: DIFFUSION_STATE_T
USE HALL_REGIME,        ONLY: HALL_INIT, HALL_ADVANCE
USE TIMESTEPPER,        ONLY: RUN
USE UNITS,              ONLY: ENERGY_UNIT_ERG
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R     = 40
REAL(KIND=dp),    PARAMETER :: R_MIN   = 9.0_dp    ! km -- core-crust boundary
REAL(KIND=dp),    PARAMETER :: R_MAX   = 10.0_dp   ! km -- stellar surface
INTEGER(KIND=i4), PARAMETER :: LMAX    = 30
REAL(KIND=dp),    PARAMETER :: ETA     = 1.0E-6_dp ! km**2/yr
REAL(KIND=dp),    PARAMETER :: F_HALL  = 0.01_dp   ! km**2/(1e12 G)/yr
INTEGER(KIND=i4), PARAMETER :: N_SUB   = 1          ! dt_hall = DT exactly
INTEGER(KIND=i4), PARAMETER :: DATA_UNIT = 25
REAL(KIND=dp),    PARAMETER :: PI = 3.14159265358979_dp
REAL(KIND=dp),    PARAMETER :: BLOWUP_FACTOR = 10.0_dp
INTEGER(KIND=i4), PARAMETER :: LOG_EVERY = 50

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(DIFFUSION_STATE_T) :: STATE
REAL(KIND=dp) :: E0, E_NOW, DT, TMAX
INTEGER(KIND=i4) :: IR, N_STEPS
CHARACTER(LEN=1024) :: DATA_PATH, ARG_STR
LOGICAL :: WRITE_DATA, BLEW_UP

IF (COMMAND_ARGUMENT_COUNT() < 2) THEN
  WRITE(*,'(A)') 'usage: mhdvsh_hall_dt_sweep <dt_yr> <tmax_yr> [output_data_path]'
  STOP 1
END IF
CALL GET_COMMAND_ARGUMENT(1, ARG_STR); READ(ARG_STR,*) DT
CALL GET_COMMAND_ARGUMENT(2, ARG_STR); READ(ARG_STR,*) TMAX
N_STEPS = NINT(TMAX/DT)

WRITE_DATA = (COMMAND_ARGUMENT_COUNT() >= 3)
IF (WRITE_DATA) CALL GET_COMMAND_ARGUMENT(3, DATA_PATH)

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)
CALL HALL_INIT(RGRID, OPS, LMAX, ETA, F_HALL, DT, N_SUB)

CALL ALLOC_SPECTRAL_SCALAR(STATE%PHI, N_R, LMAX)
CALL ALLOC_SPECTRAL_SCALAR(STATE%PSI, N_R, LMAX)

DO IR = 1, N_R
  STATE%PHI%COEF(IR, YLM_INDEX(1_i4,0_i4)) = &
    CMPLX(SIN(PI*(RGRID%R(IR)-R_MIN)/(R_MAX-R_MIN)), 0.0_dp, KIND=dp)
END DO

WRITE(*,'(A)')            'MHD-VSH Hall dt-stability sweep point'
WRITE(*,'(A,ES12.5,A)')   '  dt      = ', DT, ' yr (= dt_hall, N_SUB=1)'
WRITE(*,'(A,ES12.5,A,I0,A)') '  tmax    = ', TMAX, ' yr (', N_STEPS, ' steps)'

IF (WRITE_DATA) THEN
  OPEN(UNIT=DATA_UNIT, FILE=TRIM(DATA_PATH), STATUS='REPLACE', ACTION='WRITE')
  WRITE(DATA_UNIT,'(A)') '# step  t_yr  E_total_erg'
END IF

E0 = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)
IF (WRITE_DATA) WRITE(DATA_UNIT,'(I8,2ES16.8)') 0_i4, 0.0_dp, E0*ENERGY_UNIT_ERG
BLEW_UP = .FALSE.

CALL RUN(HALL_ADVANCE, STATE, DT, N_STEPS, ON_STEP=CHECK_STEP)

IF (.NOT. BLEW_UP) THEN
  E_NOW = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)
  WRITE(*,'(A,ES14.6,A)') '  E(t=0)   = ', E0*ENERGY_UNIT_ERG, ' erg'
  WRITE(*,'(A,ES14.6,A)') '  E(t=end) = ', E_NOW*ENERGY_UNIT_ERG, ' erg'
  WRITE(*,'(A)') '  RESULT: STABLE -- completed without blowup'
END IF
IF (WRITE_DATA) CLOSE(DATA_UNIT)

CONTAINS

!> Matches TIMESTEPPER::REGIME_ON_STEP_I. Checks TOTAL_MAGNETIC_ENERGY
!> for NaN or `> BLOWUP_FACTOR*E0` every LOG_EVERY steps (same convention
!> as app/mhdvsh_hall_stability_experiment.f90); logs to DATA_UNIT if
!> requested. STOPs the whole program on the first blowup detected --
!> TIMESTEPPER::RUN has no early-exit hook, so this is the only way to
!> avoid burning the rest of an already-diverged run's compute budget.
SUBROUTINE CHECK_STEP(STATE_ARG, T, ISTEP)
  CLASS(*),         INTENT(IN) :: STATE_ARG
  REAL(KIND=dp),    INTENT(IN) :: T
  INTEGER(KIND=i4), INTENT(IN) :: ISTEP
  IF (MOD(ISTEP, LOG_EVERY) /= 0 .AND. ISTEP /= N_STEPS) RETURN
  SELECT TYPE (STATE_ARG)
  TYPE IS (DIFFUSION_STATE_T)
    E_NOW = TOTAL_MAGNETIC_ENERGY(STATE_ARG%PHI, STATE_ARG%PSI, OPS, RGRID)
  END SELECT
  IF (WRITE_DATA) WRITE(DATA_UNIT,'(I8,2ES16.8)') ISTEP, T, E_NOW*ENERGY_UNIT_ERG
  IF (IEEE_IS_NAN(E_NOW)) THEN
    WRITE(*,'(A,ES12.5,A)') '  RESULT: BLOWUP (NaN) at t=', T, ' yr'
    BLEW_UP = .TRUE.
  ELSE IF (E_NOW > BLOWUP_FACTOR*E0) THEN
    WRITE(*,'(A,F6.2,A,ES12.5,A)') '  RESULT: BLOWUP (E > ', BLOWUP_FACTOR, &
      'x initial) at t=', T, ' yr'
    BLEW_UP = .TRUE.
  END IF
  IF (BLEW_UP) THEN
    IF (WRITE_DATA) CLOSE(DATA_UNIT)
    STOP 1
  END IF
END SUBROUTINE CHECK_STEP

END PROGRAM MHDVSH_HALL_DT_SWEEP
