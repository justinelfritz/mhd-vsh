!> Combined resistive+Hall regime driver: seeds a single poloidal mode
!> on a shell grid, evolves it under HALL_REGIME (implicit diffusion +
!> explicit Hall substeps, see hall_regime.f90's own header for the
!> splitting/BC design), and reports magnetic energy before/after. This
!> is the "primary interest" case -- resistive and Hall effects
!> together, not either studied in isolation (see mhdvsh_diffusion.f90
!> for pure diffusion, mhdvsh_hall_stability_experiment.f90 for pure
!> Hall).
!>
!> Optional command-line arguments (all positional, each requires the
!> ones before it):
!>   1: a file path to write a time-resolved energy-budget series to
!>      (step, t, E_poloidal, E_toroidal, joule_dissipation_rate,
!>      poynting_flux_rate, hall_poynting_flux_rate,
!>      energy_balance_residual) -- every term in the combined balance
!>      identity Edot_B,pol + Edot_B,tor = Edot_J + Edot_S,diffusion +
!>      Edot_S,Hall (Hall's own Joule term is exactly zero, so it
!>      doesn't appear -- see JOULE_DISSIPATION_RATE's own docstring),
!>      plus that identity's residual. scripts/plot_energy_budget.py
!>      needs no changes to plot this -- it's already written
!>      generically against the column set.
!>   2: a file path to write per-degree poloidal/toroidal energy to
!>      (see LOG_ENERGY_BY_L).
!>   3: a checkpoint file path -- if given, IO_CHECKPOINT::WRITE_CHECKPOINT
!>      is called every CHECKPOINT_EVERY steps (plus once, unconditionally,
!>      after the run completes) so a crash mid-run (see ROADMAP.md,
!>      "Recently resolved (2026-08-21, even later)" -- the reason this
!>      capability exists at all) loses at most CHECKPOINT_EVERY steps
!>      of work, not the whole run.
!>   4: a checkpoint file path to RESUME from (via
!>      IO_CHECKPOINT::READ_CHECKPOINT) instead of starting from the
!>      hardcoded seed IC at step 0 -- the loaded grid shape (N_R, LMAX,
!>      R_MIN, R_MAX) is checked against this driver's own compiled-in
!>      values before trusting the loaded state (see IO_CHECKPOINT's own
!>      @warning: a checkpoint does not itself carry ETA/F_HALL/DT/N_SUB,
!>      resuming means re-running this SAME driver, not a different
!>      configuration). N_STEPS below is the TOTAL target step count,
!>      not "how many more steps to run" -- resuming from step 12000
!>      with N_STEPS=50000 runs the remaining 38000. Arguments 1/2's
!>      files (if given) are appended to rather than overwritten, and
!>      the row already written for the resumed step is not repeated.
!>
!> Units: same convention as mhdvsh_diffusion.f90 (UNITS::ENERGY_UNIT_ERG/
!> POWER_UNIT_ERG_PER_S applied only at this reporting layer).
PROGRAM MHDVSH_HALL
USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,   ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS,  ONLY: TOTAL_MAGNETIC_ENERGY, TOTAL_POLOIDAL_MAGNETIC_ENERGY, &
                               TOTAL_TOROIDAL_MAGNETIC_ENERGY, JOULE_DISSIPATION_RATE, &
                               POYNTING_FLUX_RATE, HALL_POYNTING_FLUX_RATE, &
                               POLOIDAL_MAGNETIC_ENERGY_BY_L, TOROIDAL_MAGNETIC_ENERGY_BY_L
USE DIFFUSION_REGIME,   ONLY: DIFFUSION_STATE_T
USE HALL_REGIME,        ONLY: HALL_INIT, HALL_ADVANCE
USE TIMESTEPPER,        ONLY: RUN
USE UNITS,              ONLY: ENERGY_UNIT_ERG, POWER_UNIT_ERG_PER_S
USE IO_CHECKPOINT,      ONLY: WRITE_CHECKPOINT, READ_CHECKPOINT
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R     = 40
! Crust extent of the M=1.40 reference star (rhocgs=9.88d14 g/cm**3),
! from mhdvsh_tov (src/core/tov_solver.f90 -- ported TOV solver,
! regression-tested against ~/Desktop/EOSNS/fort.34/PL.DAT), superseding
! the earlier 9.0/10.0 km placeholder. Rerunning this driver's own
! dt-sweep/1000yr archive at these corrected radii is deferred (see the
! TOV/EOS port's own plan notes), not part of this update.
REAL(KIND=dp),    PARAMETER :: R_MIN   = 10.8033325018_dp   ! km -- core-crust boundary
REAL(KIND=dp),    PARAMETER :: R_MAX   = 11.6982211606_dp   ! km -- stellar surface
INTEGER(KIND=i4), PARAMETER :: LMAX    = 30
REAL(KIND=dp),    PARAMETER :: ETA     = 1.0E-6_dp ! km**2/yr -- realistic crustal value, per user (2026-08-21): 1e-8 to 1e-5 range
REAL(KIND=dp),    PARAMETER :: F_HALL  = 0.01_dp  ! km**2/(1e12 G)/yr
REAL(KIND=dp),    PARAMETER :: DT      = 0.01_dp  ! yr
INTEGER(KIND=i4), PARAMETER :: N_SUB   = 10        ! -> dt_hall=0.001, per Stage 2's empirical data
INTEGER(KIND=i4), PARAMETER :: N_STEPS = 200       ! TOTAL target step count (see CLI arg 4's own doc)
INTEGER(KIND=i4), PARAMETER :: CHECKPOINT_EVERY = 50  ! tune coarser for long production runs
INTEGER(KIND=i4), PARAMETER :: DATA_UNIT = 21
INTEGER(KIND=i4), PARAMETER :: ENERGY_L_UNIT = 24
REAL(KIND=dp),    PARAMETER :: PI = 3.14159265358979_dp
REAL(KIND=dp),    PARAMETER :: GRID_TOL = 1.0E-9_dp   ! resume grid-shape check

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(DIFFUSION_STATE_T) :: STATE
REAL(KIND=dp) :: E0, E1
REAL(KIND=dp) :: PREV_E_POL, PREV_E_TOR
INTEGER(KIND=i4) :: IR, L
CHARACTER(LEN=1024) :: DATA_PATH, ENERGY_L_PATH, CHECKPOINT_PATH, RESUME_PATH
CHARACTER(LEN=16) :: COL_LABEL
LOGICAL :: WRITE_DATA, WRITE_ENERGY_L, WRITE_CHECKPOINT_FLAG, RESUMING, NEED_ON_STEP
INTEGER(KIND=i4) :: ISTEP_START, N_STEPS_REMAINING
INTEGER(KIND=i4) :: N_R_CK, LMAX_CK
REAL(KIND=dp)    :: T_START, R_MIN_CK, R_MAX_CK

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)
CALL HALL_INIT(RGRID, OPS, LMAX, ETA, F_HALL, DT, N_SUB)

RESUMING = (COMMAND_ARGUMENT_COUNT() >= 4)
IF (RESUMING) THEN
  CALL GET_COMMAND_ARGUMENT(4, RESUME_PATH)
  CALL READ_CHECKPOINT(TRIM(RESUME_PATH), ISTEP_START, T_START, STATE%PHI, STATE%PSI, &
    N_R_CK, LMAX_CK, R_MIN_CK, R_MAX_CK)
  IF (N_R_CK /= N_R .OR. LMAX_CK /= LMAX .OR. &
      ABS(R_MIN_CK-R_MIN) > GRID_TOL .OR. ABS(R_MAX_CK-R_MAX) > GRID_TOL) THEN
    WRITE(*,'(A)') 'MHDVSH_HALL: checkpoint grid shape does not match this driver''s own parameters'
    WRITE(*,'(A,I0,A,I0,A,ES12.5,A,ES12.5)') '  checkpoint: N_R=', N_R_CK, ' LMAX=', LMAX_CK, &
      ' R_MIN=', R_MIN_CK, ' R_MAX=', R_MAX_CK
    WRITE(*,'(A,I0,A,I0,A,ES12.5,A,ES12.5)') '  this run:   N_R=', N_R, ' LMAX=', LMAX, &
      ' R_MIN=', R_MIN, ' R_MAX=', R_MAX
    STOP 1
  END IF
ELSE
  ISTEP_START = 0_i4
  T_START = 0.0_dp
  CALL ALLOC_SPECTRAL_SCALAR(STATE%PHI, N_R, LMAX)
  CALL ALLOC_SPECTRAL_SCALAR(STATE%PSI, N_R, LMAX)

  ! Single seed mode, same as app/mhdvsh_hall_stability_experiment.f90 --
  ! sin() profile vanishes at both ends, Psi stays zero (never set).
  DO IR = 1, N_R
    STATE%PHI%COEF(IR, YLM_INDEX(1_i4,0_i4)) = &
      CMPLX(SIN(PI*(RGRID%R(IR)-R_MIN)/(R_MAX-R_MIN)), 0.0_dp, KIND=dp)
  END DO
END IF

N_STEPS_REMAINING = N_STEPS - ISTEP_START
IF (N_STEPS_REMAINING <= 0) THEN
  WRITE(*,'(A,I0,A,I0,A)') 'MHDVSH_HALL: checkpoint step ', ISTEP_START, &
    ' already reaches or exceeds N_STEPS=', N_STEPS, ' -- nothing to do'
  STOP 1
END IF

WRITE_CHECKPOINT_FLAG = (COMMAND_ARGUMENT_COUNT() >= 3)
IF (WRITE_CHECKPOINT_FLAG) CALL GET_COMMAND_ARGUMENT(3, CHECKPOINT_PATH)

WRITE_DATA = (COMMAND_ARGUMENT_COUNT() >= 1)
IF (WRITE_DATA) THEN
  CALL GET_COMMAND_ARGUMENT(1, DATA_PATH)
  IF (RESUMING) THEN
    OPEN(UNIT=DATA_UNIT, FILE=TRIM(DATA_PATH), STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
  ELSE
    OPEN(UNIT=DATA_UNIT, FILE=TRIM(DATA_PATH), STATUS='REPLACE', ACTION='WRITE')
    WRITE(DATA_UNIT,'(A)') '# step  t_yr  E_poloidal_erg  E_toroidal_erg' // &
      '  joule_dissipation_rate_erg_per_s  poynting_flux_rate_erg_per_s' // &
      '  hall_poynting_flux_rate_erg_per_s  energy_balance_residual_erg_per_s'
  END IF
  PREV_E_POL = TOTAL_POLOIDAL_MAGNETIC_ENERGY(STATE%PHI, OPS, RGRID)
  PREV_E_TOR = TOTAL_TOROIDAL_MAGNETIC_ENERGY(STATE%PSI, RGRID)
  IF (.NOT. RESUMING) THEN
    ! No previous sample exists yet to backward-difference against, so
    ! the residual is undefined (not 0) at step 0 -- NaN leaves a
    ! visible gap in the plot rather than implying a (meaningless)
    ! perfect balance. On resume this row was already written by the
    ! original run, so it's not repeated.
    WRITE(DATA_UNIT,'(I8,7ES16.8)') 0_i4, 0.0_dp, &
      PREV_E_POL*ENERGY_UNIT_ERG, PREV_E_TOR*ENERGY_UNIT_ERG, &
      JOULE_DISSIPATION_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA)*POWER_UNIT_ERG_PER_S, &
      POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA)*POWER_UNIT_ERG_PER_S, &
      HALL_POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, F_HALL)*POWER_UNIT_ERG_PER_S, &
      IEEE_VALUE(1.0_dp, IEEE_QUIET_NAN)
  END IF

  ! Per-degree magnetic energy (poloidal from Phi, toroidal from Psi) --
  ! FIELD_DIAGNOSTICS::POLOIDAL_MAGNETIC_ENERGY_BY_L/
  ! TOROIDAL_MAGNETIC_ENERGY_BY_L, same functions/convention
  ! app/mhdvsh_hall_stability_experiment.f90 uses, now for the combined
  ! resistive+Hall regime instead of pure Hall.
  WRITE_ENERGY_L = (COMMAND_ARGUMENT_COUNT() >= 2)
  IF (WRITE_ENERGY_L) THEN
    CALL GET_COMMAND_ARGUMENT(2, ENERGY_L_PATH)
    IF (RESUMING) THEN
      OPEN(UNIT=ENERGY_L_UNIT, FILE=TRIM(ENERGY_L_PATH), STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
    ELSE
      OPEN(UNIT=ENERGY_L_UNIT, FILE=TRIM(ENERGY_L_PATH), STATUS='REPLACE', ACTION='WRITE')
      WRITE(ENERGY_L_UNIT,'(A)', ADVANCE='NO') '# step  t'
      DO L = 1, LMAX
        WRITE(COL_LABEL,'(A,I0)') '  Epol_l', L
        WRITE(ENERGY_L_UNIT,'(A)', ADVANCE='NO') TRIM(COL_LABEL)
      END DO
      DO L = 1, LMAX
        WRITE(COL_LABEL,'(A,I0)') '  Etor_l', L
        WRITE(ENERGY_L_UNIT,'(A)', ADVANCE='NO') TRIM(COL_LABEL)
      END DO
      WRITE(ENERGY_L_UNIT,'(A)') ''
      CALL LOG_ENERGY_BY_L(0_i4, 0.0_dp)
    END IF
  END IF
END IF

NEED_ON_STEP = WRITE_DATA .OR. WRITE_CHECKPOINT_FLAG

E0 = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)
IF (NEED_ON_STEP) THEN
  CALL RUN(HALL_ADVANCE, STATE, DT, N_STEPS_REMAINING, T_START=T_START, ON_STEP=LOG_ENERGY)
ELSE
  CALL RUN(HALL_ADVANCE, STATE, DT, N_STEPS_REMAINING, T_START=T_START)
END IF
E1 = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)

! Unconditional final checkpoint, regardless of CHECKPOINT_EVERY
! alignment, so a completed (or Ctrl-C'd) run always leaves a
! checkpoint reflecting its true final state.
IF (WRITE_CHECKPOINT_FLAG) THEN
  CALL WRITE_CHECKPOINT(TRIM(CHECKPOINT_PATH), N_STEPS, T_START + REAL(N_STEPS_REMAINING,KIND=dp)*DT, &
    STATE%PHI, STATE%PSI, RGRID)
END IF

IF (WRITE_DATA) CLOSE(DATA_UNIT)
IF (WRITE_ENERGY_L) CLOSE(ENERGY_L_UNIT)

WRITE(*,'(A)')            'MHD-VSH combined resistive+Hall regime'
WRITE(*,'(A,I0)')         '  N_r          = ', N_R
WRITE(*,'(A,I0)')         '  Lmax         = ', LMAX
WRITE(*,'(A,ES12.5,A)')   '  eta          = ', ETA, ' km**2/yr'
WRITE(*,'(A,ES12.5,A)')   '  F_Hall       = ', F_HALL, ' km**2/(1e12 G)/yr'
WRITE(*,'(A,ES12.5,A,I0,A,I0,A)') '  dt           = ', DT, ' yr (', N_STEPS, &
  ' steps total, ', N_SUB, ' Hall substeps each)'
IF (RESUMING) WRITE(*,'(A,I0,A,A)') '  resumed from step ', ISTEP_START, ' via ', TRIM(RESUME_PATH)
WRITE(*,'(A,ES12.5,A)')   '  energy(t=0)  = ', E0*ENERGY_UNIT_ERG, ' erg'
WRITE(*,'(A,ES12.5,A)')   '  energy(t=end)= ', E1*ENERGY_UNIT_ERG, ' erg'
IF (E1 < E0) THEN
  WRITE(*,'(A)') '  RESULT: energy decreased (diffusion dissipating, as expected)'
ELSE
  WRITE(*,'(A)') '  RESULT: UNEXPECTED -- energy did not decrease'
END IF
IF (WRITE_DATA) WRITE(*,'(A,A)') '  time series written to ', TRIM(DATA_PATH)
IF (WRITE_CHECKPOINT_FLAG) WRITE(*,'(A,A)') '  final checkpoint written to ', TRIM(CHECKPOINT_PATH)

CONTAINS

!> Matches TIMESTEPPER::REGIME_ON_STEP_I; logs every combined-regime
!> energy-balance term (when WRITE_DATA) and writes a periodic
!> checkpoint (when WRITE_CHECKPOINT_FLAG). ISTEP is RUN's own local
!> counter (1..N_STEPS_REMAINING for THIS call) -- GLOBAL_ISTEP adds
!> back ISTEP_START so logged step numbers and the checkpoint cadence
!> are correct across a resume, not restarting from 1 each time.
!> Residual = d(E_pol+E_tor)/dt (backward-differenced) minus
!> (Edot_J + Edot_S,diffusion + Edot_S,Hall); per the combined balance
!> identity this should sit at ~0. Same pattern as
!> mhdvsh_diffusion.f90's LOG_ENERGY, extended with the Hall Poynting
!> term and checkpointing.
SUBROUTINE LOG_ENERGY(STATE, T, ISTEP)
  CLASS(*),         INTENT(IN) :: STATE
  REAL(KIND=dp),    INTENT(IN) :: T
  INTEGER(KIND=i4), INTENT(IN) :: ISTEP
  REAL(KIND=dp) :: E_POL, E_TOR, EDOT_J, EDOT_S, EDOT_S_HALL, RESIDUAL
  INTEGER(KIND=i4) :: GLOBAL_ISTEP
  GLOBAL_ISTEP = ISTEP_START + ISTEP
  SELECT TYPE (STATE)
  TYPE IS (DIFFUSION_STATE_T)
    IF (WRITE_DATA) THEN
      E_POL       = TOTAL_POLOIDAL_MAGNETIC_ENERGY(STATE%PHI, OPS, RGRID)
      E_TOR       = TOTAL_TOROIDAL_MAGNETIC_ENERGY(STATE%PSI, RGRID)
      EDOT_J      = JOULE_DISSIPATION_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA)
      EDOT_S      = POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA)
      EDOT_S_HALL = HALL_POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, F_HALL)
      RESIDUAL = ((E_POL-PREV_E_POL) + (E_TOR-PREV_E_TOR))/DT - (EDOT_J+EDOT_S+EDOT_S_HALL)

      WRITE(DATA_UNIT,'(I8,7ES16.8)') GLOBAL_ISTEP, T, E_POL*ENERGY_UNIT_ERG, E_TOR*ENERGY_UNIT_ERG, &
        EDOT_J*POWER_UNIT_ERG_PER_S, EDOT_S*POWER_UNIT_ERG_PER_S, &
        EDOT_S_HALL*POWER_UNIT_ERG_PER_S, RESIDUAL*POWER_UNIT_ERG_PER_S

      PREV_E_POL = E_POL
      PREV_E_TOR = E_TOR
      IF (WRITE_ENERGY_L) CALL LOG_ENERGY_BY_L(GLOBAL_ISTEP, T)
    END IF

    IF (WRITE_CHECKPOINT_FLAG) THEN
      IF (MOD(GLOBAL_ISTEP, CHECKPOINT_EVERY) == 0) THEN
        CALL WRITE_CHECKPOINT(TRIM(CHECKPOINT_PATH), GLOBAL_ISTEP, T, STATE%PHI, STATE%PSI, RGRID)
      END IF
    END IF
  END SELECT
END SUBROUTINE LOG_ENERGY

!> Writes one row to ENERGY_L_UNIT: step, t, POLOIDAL_MAGNETIC_ENERGY_BY_L
!> (l=1..LMAX) then TOROIDAL_MAGNETIC_ENERGY_BY_L (l=1..LMAX), both in
!> erg (l=0 omitted -- always exactly zero, see HALL_INDUCTION_RHS's own
!> n=0 exclusion). Host-associates STATE/OPS/RGRID from the main program.
SUBROUTINE LOG_ENERGY_BY_L(ISTEP_ARG, T_ARG)
  INTEGER(KIND=i4), INTENT(IN) :: ISTEP_ARG
  REAL(KIND=dp),    INTENT(IN) :: T_ARG
  REAL(KIND=dp) :: E_POL_L(0:LMAX), E_TOR_L(0:LMAX)
  INTEGER(KIND=i4) :: LL
  E_POL_L = POLOIDAL_MAGNETIC_ENERGY_BY_L(STATE%PHI, OPS, RGRID)
  E_TOR_L = TOROIDAL_MAGNETIC_ENERGY_BY_L(STATE%PSI, RGRID)
  WRITE(ENERGY_L_UNIT,'(I8,ES14.6)', ADVANCE='NO') ISTEP_ARG, T_ARG
  DO LL = 1, LMAX
    WRITE(ENERGY_L_UNIT,'(ES14.6)', ADVANCE='NO') E_POL_L(LL)*ENERGY_UNIT_ERG
  END DO
  DO LL = 1, LMAX
    WRITE(ENERGY_L_UNIT,'(ES14.6)', ADVANCE='NO') E_TOR_L(LL)*ENERGY_UNIT_ERG
  END DO
  WRITE(ENERGY_L_UNIT,'(A)') ''
END SUBROUTINE LOG_ENERGY_BY_L

END PROGRAM MHDVSH_HALL
