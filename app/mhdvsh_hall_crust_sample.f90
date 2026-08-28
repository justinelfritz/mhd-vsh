!> Combined resistive+Hall regime driver -- identical grid/timestepper
!> setup to app/mhdvsh_hall.f90 (same N_R, LMAX, R_MIN, R_MAX, DT, N_SUB,
!> N_STEPS, seed IC), but with ETA/F_HALL replaced by an actual
!> physically-sourced (eta,f_H) pair instead of mhdvsh_hall.f90's own
!> round toy values (ETA=1e-6, F_HALL=0.01) -- a separate, deliberately
!> non-overwriting comparison run, not a modification of the existing
!> baseline driver (per user, 2026-08-24).
!>
!> ---------------------------------------------------------------------
!> SIMULATION INPUT PROVENANCE (documented here for later review):
!>   Source: results/eos_comparison/new_eos_profile.dat, row at
!>     r_km=1.15570341E+01, rhocgs=4.68605102E+06, n_e_cm-3=1.31075626E+30
!>     (produced by app/mhdvsh_tov.f90 from the current CONDUCT/COUL19-
!>     based crust conductivity engine, src/core/crust_conductivity.f90,
!>     see that module's own header for the underlying physics/citations).
!>   ETA    = 1.03614247E-3 km**2/yr  (that row's own eta_km2_per_yr column)
!>   F_HALL = 1.19580898E-2 km**2/(1e12 G)/yr  (that row's own f_hall column)
!>   Selected by user directly from the profile data file (IDE selection,
!>   2026-08-24) -- this radius sits in the outer-crust degenerate/
!>   non-degenerate transition band (theta=kT/E_F ~ 0.5-0.6 at this
!>   project's default T=1e9 K isothermal-crust assumption; see this
!>   session's own eta(T) scaling discussion for that characterization).
!>   All OTHER parameters (N_R, LMAX, R_MIN, R_MAX, DT, N_SUB, N_STEPS,
!>   CHECKPOINT_EVERY, seed IC) are IDENTICAL to app/mhdvsh_hall.f90's
!>   own values, unchanged, so this run is a controlled single-variable
!>   (eta, f_H) comparison against that existing baseline -- not a
!>   different grid/timestep/resolution study.
!> ---------------------------------------------------------------------
!>
!> Optional command-line arguments (all positional, each requires the
!> ones before it) -- identical meaning to mhdvsh_hall.f90's own:
!>   1: energy-budget time series output path
!>   2: per-degree poloidal/toroidal energy output path
!>   3: checkpoint output path (written every CHECKPOINT_EVERY steps,
!>      plus once unconditionally at the end)
!>   4: checkpoint path to RESUME from
!>
!> Units: same convention as mhdvsh_diffusion.f90/mhdvsh_hall.f90
!> (UNITS::ENERGY_UNIT_ERG/POWER_UNIT_ERG_PER_S applied only at this
!> reporting layer).
PROGRAM MHDVSH_HALL_CRUST_SAMPLE
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
REAL(KIND=dp),    PARAMETER :: R_MIN   = 10.3029378_dp   ! km -- core-crust boundary, matches mhdvsh_hall.f90
REAL(KIND=dp),    PARAMETER :: R_MAX   = 11.5632834_dp   ! km -- stellar surface, matches mhdvsh_hall.f90
INTEGER(KIND=i4), PARAMETER :: LMAX    = 30
! Physically-sourced (eta, f_H) pair -- see module header's PROVENANCE
! block above for exact source row/radius/density. NOT a toy value.
REAL(KIND=dp),    PARAMETER :: ETA     = 1.03614247E-3_dp ! km**2/yr
REAL(KIND=dp),    PARAMETER :: F_HALL  = 1.19580898E-2_dp ! km**2/(1e12 G)/yr
REAL(KIND=dp),    PARAMETER :: DT      = 0.01_dp  ! yr -- matches mhdvsh_hall.f90
INTEGER(KIND=i4), PARAMETER :: N_SUB   = 10        ! matches mhdvsh_hall.f90
INTEGER(KIND=i4), PARAMETER :: N_STEPS = 200       ! matches mhdvsh_hall.f90
INTEGER(KIND=i4), PARAMETER :: CHECKPOINT_EVERY = 50
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
    WRITE(*,'(A)') 'MHDVSH_HALL_CRUST_SAMPLE: checkpoint grid shape does not match this driver''s own parameters'
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

  ! Single seed mode, identical to mhdvsh_hall.f90's own IC.
  DO IR = 1, N_R
    STATE%PHI%COEF(IR, YLM_INDEX(1_i4,0_i4)) = &
      CMPLX(SIN(PI*(RGRID%R(IR)-R_MIN)/(R_MAX-R_MIN)), 0.0_dp, KIND=dp)
  END DO
END IF

N_STEPS_REMAINING = N_STEPS - ISTEP_START
IF (N_STEPS_REMAINING <= 0) THEN
  WRITE(*,'(A,I0,A,I0,A)') 'MHDVSH_HALL_CRUST_SAMPLE: checkpoint step ', ISTEP_START, &
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
    ! Full input provenance written into the data file itself (not just
    ! the source) so the record survives independent of this file's own
    ! future edits -- "document all simulation inputs for later review"
    ! (user, 2026-08-24).
    WRITE(DATA_UNIT,'(A)') '# MHDVSH_HALL_CRUST_SAMPLE run inputs:'
    WRITE(DATA_UNIT,'(A,I0)')       '#   N_R          = ', N_R
    WRITE(DATA_UNIT,'(A,I0)')       '#   LMAX         = ', LMAX
    WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   R_MIN        = ', R_MIN, ' km'
    WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   R_MAX        = ', R_MAX, ' km'
    WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   ETA          = ', ETA, ' km**2/yr'
    WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   F_HALL       = ', F_HALL, ' km**2/(1e12 G)/yr'
    WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   DT           = ', DT, ' yr'
    WRITE(DATA_UNIT,'(A,I0)')       '#   N_SUB        = ', N_SUB
    WRITE(DATA_UNIT,'(A,I0)')       '#   N_STEPS      = ', N_STEPS
    WRITE(DATA_UNIT,'(A)') '#   eta/f_H source: results/eos_comparison/new_eos_profile.dat,' // &
      ' row r_km=1.15570341E+01 rhocgs=4.68605102E+06 n_e_cm-3=1.31075626E+30'
    WRITE(DATA_UNIT,'(A)') '#   seed IC: single mode Phi(l=1,m=0)=sin(pi*(r-R_MIN)/(R_MAX-R_MIN)), Psi=0'
    WRITE(DATA_UNIT,'(A)') '# step  t_yr  E_poloidal_erg  E_toroidal_erg' // &
      '  joule_dissipation_rate_erg_per_s  poynting_flux_rate_erg_per_s' // &
      '  hall_poynting_flux_rate_erg_per_s  energy_balance_residual_erg_per_s'
  END IF
  PREV_E_POL = TOTAL_POLOIDAL_MAGNETIC_ENERGY(STATE%PHI, OPS, RGRID)
  PREV_E_TOR = TOTAL_TOROIDAL_MAGNETIC_ENERGY(STATE%PSI, RGRID)
  IF (.NOT. RESUMING) THEN
    WRITE(DATA_UNIT,'(I8,7ES16.8)') 0_i4, 0.0_dp, &
      PREV_E_POL*ENERGY_UNIT_ERG, PREV_E_TOR*ENERGY_UNIT_ERG, &
      JOULE_DISSIPATION_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA)*POWER_UNIT_ERG_PER_S, &
      POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA)*POWER_UNIT_ERG_PER_S, &
      HALL_POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, F_HALL)*POWER_UNIT_ERG_PER_S, &
      IEEE_VALUE(1.0_dp, IEEE_QUIET_NAN)
  END IF

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

IF (WRITE_CHECKPOINT_FLAG) THEN
  CALL WRITE_CHECKPOINT(TRIM(CHECKPOINT_PATH), N_STEPS, T_START + REAL(N_STEPS_REMAINING,KIND=dp)*DT, &
    STATE%PHI, STATE%PSI, RGRID)
END IF

IF (WRITE_DATA) CLOSE(DATA_UNIT)
IF (WRITE_ENERGY_L) CLOSE(ENERGY_L_UNIT)

WRITE(*,'(A)')            'MHD-VSH combined resistive+Hall regime -- crust-sampled (eta,f_H)'
WRITE(*,'(A)')            '  eta/f_H source = new_eos_profile.dat, r=11.5570341 km, rho=4.68605102E+06 g/cm**3'
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

!> Matches TIMESTEPPER::REGIME_ON_STEP_I; see mhdvsh_hall.f90's own
!> identical routine for the full explanation -- unchanged here.
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

!> Identical to mhdvsh_hall.f90's own -- see that file's docstring.
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

END PROGRAM MHDVSH_HALL_CRUST_SAMPLE
