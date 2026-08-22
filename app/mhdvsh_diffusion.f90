!> Diffusion-regime driver: seeds a smooth (Phi,Psi) profile on a shell
!> grid, evolves it under pure Ohmic diffusion via DIFFUSION_REGIME +
!> TIMESTEPPER, and reports the magnetic energy before/after as a
!> sanity check (pure diffusion must dissipate, never create, energy).
!>
!> Optional command-line argument: a file path to write a time-resolved
!> energy-budget series to (step, t, E_poloidal, E_toroidal,
!> joule_dissipation_rate, poynting_flux_rate, energy_balance_residual),
!> one row per step, for scripts/plot_energy_budget.py to read -- every
!> term in analytic_formulas/mhd-vsh-relations.tex's energy balance
!> (Edot_B,pol + Edot_B,tor = Edot_J + Edot_S), plus that identity's own
!> residual (d(E_pol+E_tor)/dt, backward-differenced across steps, minus
!> Joule+Poynting), which should sit at ~0 (conservation) except during
!> the rapid initial transient, where first-order (backward-Euler) time
!> truncation error dominates.
!>
!> Units (see UNITS, units.f90): R_MIN/R_MAX/RGRID%R are in km, ETA in
!> km**2/yr, DT/t in yr -- FIELD_DIAGNOSTICS' formulas don't care about
!> this choice (unit-agnostic), so nothing about the solve itself
!> changes; only the energy/rate columns actually written below (and
!> printed to the console) are converted, via UNITS' ENERGY_UNIT_ERG/
!> POWER_UNIT_ERG_PER_S, from code units to erg/erg-per-second. The
!> energy-balance RESIDUAL is computed in code units first (matching how
!> it's derived, as a difference/ratio of code-unit quantities) and only
!> converted to erg/s for the final printed/logged value.
PROGRAM MHDVSH_DIFFUSION
USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,   ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS,  ONLY: TOTAL_MAGNETIC_ENERGY, TOTAL_POLOIDAL_MAGNETIC_ENERGY, &
                               TOTAL_TOROIDAL_MAGNETIC_ENERGY, JOULE_DISSIPATION_RATE, &
                               POYNTING_FLUX_RATE
USE DIFFUSION_REGIME,   ONLY: DIFFUSION_STATE_T, DIFFUSION_INIT, DIFFUSION_ADVANCE
USE TIMESTEPPER,        ONLY: RUN
USE UNITS,              ONLY: ENERGY_UNIT_ERG, POWER_UNIT_ERG_PER_S
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R     = 40
REAL(KIND=dp),    PARAMETER :: R_MIN   = 0.5_dp   ! km
REAL(KIND=dp),    PARAMETER :: R_MAX   = 1.0_dp   ! km
INTEGER(KIND=i4), PARAMETER :: LMAX    = 3
REAL(KIND=dp),    PARAMETER :: ETA     = 0.05_dp  ! km**2/yr
REAL(KIND=dp),    PARAMETER :: DT      = 0.01_dp  ! yr
INTEGER(KIND=i4), PARAMETER :: N_STEPS = 200
INTEGER(KIND=i4), PARAMETER :: DATA_UNIT = 21

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(DIFFUSION_STATE_T) :: STATE
REAL(KIND=dp) :: E0, E1
REAL(KIND=dp) :: PREV_E_POL, PREV_E_TOR
INTEGER(KIND=i4) :: IR
CHARACTER(LEN=1024) :: DATA_PATH
LOGICAL :: WRITE_DATA

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)
CALL DIFFUSION_INIT(RGRID, OPS, LMAX, ETA, DT)

CALL ALLOC_SPECTRAL_SCALAR(STATE%PHI, N_R, LMAX)
CALL ALLOC_SPECTRAL_SCALAR(STATE%PSI, N_R, LMAX)

! Smooth seed profile on a couple of modes -- not a physically final
! initial condition (force-free ICs are still deferred), just enough to
! exercise the pipeline. sin(pi*(r-r_min)/(r_max-r_min)) already
! vanishes at both ends, though DIFFUSION_ADVANCE's first step would
! enforce that regardless.
DO IR = 1, N_R
  STATE%PHI%COEF(IR, YLM_INDEX(1,0)) = &
    CMPLX(SIN(3.14159265358979_dp*(RGRID%R(IR)-R_MIN)/(R_MAX-R_MIN)), 0.0_dp, KIND=dp)
  STATE%PSI%COEF(IR, YLM_INDEX(2,1)) = &
    CMPLX(0.5_dp*SIN(3.14159265358979_dp*(RGRID%R(IR)-R_MIN)/(R_MAX-R_MIN)), 0.0_dp, KIND=dp)
END DO

WRITE_DATA = (COMMAND_ARGUMENT_COUNT() >= 1)
IF (WRITE_DATA) THEN
  CALL GET_COMMAND_ARGUMENT(1, DATA_PATH)
  OPEN(UNIT=DATA_UNIT, FILE=TRIM(DATA_PATH), STATUS='REPLACE', ACTION='WRITE')
  WRITE(DATA_UNIT,'(A)') '# step  t_yr  E_poloidal_erg  E_toroidal_erg' // &
    '  joule_dissipation_rate_erg_per_s  poynting_flux_rate_erg_per_s' // &
    '  energy_balance_residual_erg_per_s'
  ! No previous sample exists yet to backward-difference against, so the
  ! residual is undefined (not 0) at step 0 -- NaN leaves a visible gap
  ! in the plot rather than implying a (meaningless) perfect balance.
  PREV_E_POL = TOTAL_POLOIDAL_MAGNETIC_ENERGY(STATE%PHI, OPS, RGRID)
  PREV_E_TOR = TOTAL_TOROIDAL_MAGNETIC_ENERGY(STATE%PSI, RGRID)
  WRITE(DATA_UNIT,'(I8,6ES16.8)') 0_i4, 0.0_dp, &
    PREV_E_POL*ENERGY_UNIT_ERG, PREV_E_TOR*ENERGY_UNIT_ERG, &
    JOULE_DISSIPATION_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA)*POWER_UNIT_ERG_PER_S, &
    POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA)*POWER_UNIT_ERG_PER_S, &
    IEEE_VALUE(1.0_dp, IEEE_QUIET_NAN)
END IF

E0 = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)
IF (WRITE_DATA) THEN
  CALL RUN(DIFFUSION_ADVANCE, STATE, DT, N_STEPS, ON_STEP=LOG_ENERGY)
ELSE
  CALL RUN(DIFFUSION_ADVANCE, STATE, DT, N_STEPS)
END IF
E1 = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)

IF (WRITE_DATA) CLOSE(DATA_UNIT)

WRITE(*,'(A)')            'MHD-VSH diffusion regime'
WRITE(*,'(A,I0)')         '  N_r          = ', N_R
WRITE(*,'(A,I0)')         '  Lmax         = ', LMAX
WRITE(*,'(A,ES12.5,A)')   '  eta          = ', ETA, ' km**2/yr'
WRITE(*,'(A,ES12.5,A,I0,A)') '  dt           = ', DT, ' yr (', N_STEPS, ' steps)'
WRITE(*,'(A,ES12.5,A)')   '  energy(t=0)  = ', E0*ENERGY_UNIT_ERG, ' erg'
WRITE(*,'(A,ES12.5,A)')   '  energy(t=end)= ', E1*ENERGY_UNIT_ERG, ' erg'
IF (E1 < E0) THEN
  WRITE(*,'(A)') '  RESULT: energy decreased, as expected for pure diffusion'
ELSE
  WRITE(*,'(A)') '  RESULT: UNEXPECTED -- energy did not decrease'
END IF
IF (WRITE_DATA) WRITE(*,'(A,A)') '  time series written to ', TRIM(DATA_PATH)

CONTAINS

!> Matches TIMESTEPPER::REGIME_ON_STEP_I; logs every energy-balance term
!> (step, t, E_poloidal, E_toroidal, joule_dissipation_rate,
!> poynting_flux_rate, energy_balance_residual) to the already-open
!> DATA_UNIT. The residual is d(E_pol+E_tor)/dt -- backward-differenced
!> against PREV_E_POL/PREV_E_TOR, the previous call's values, host
!> state updated at the end of this call -- minus Joule+Poynting; per
!> the balance identity this should sit at ~0. Host-associates
!> OPS/RGRID/PREV_E_POL/PREV_E_TOR from the main program rather than
!> needing module-level state, since this procedure is only ever used
!> as an actual argument within this one run.
SUBROUTINE LOG_ENERGY(STATE, T, ISTEP)
  CLASS(*),         INTENT(IN) :: STATE
  REAL(KIND=dp),    INTENT(IN) :: T
  INTEGER(KIND=i4), INTENT(IN) :: ISTEP
  REAL(KIND=dp) :: E_POL, E_TOR, EDOT_J, EDOT_S, RESIDUAL
  SELECT TYPE (STATE)
  TYPE IS (DIFFUSION_STATE_T)
    E_POL  = TOTAL_POLOIDAL_MAGNETIC_ENERGY(STATE%PHI, OPS, RGRID)
    E_TOR  = TOTAL_TOROIDAL_MAGNETIC_ENERGY(STATE%PSI, RGRID)
    EDOT_J = JOULE_DISSIPATION_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA)
    EDOT_S = POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA)
    RESIDUAL = ((E_POL-PREV_E_POL) + (E_TOR-PREV_E_TOR))/DT - (EDOT_J+EDOT_S)

    ! E_POL/E_TOR/EDOT_J/EDOT_S/RESIDUAL above stay in code units (so
    ! PREV_E_POL/PREV_E_TOR bookkeeping and the residual's own
    ! backward-difference are computed consistently); only the values
    ! actually written are converted to erg/erg-per-second.
    WRITE(DATA_UNIT,'(I8,6ES16.8)') ISTEP, T, E_POL*ENERGY_UNIT_ERG, E_TOR*ENERGY_UNIT_ERG, &
      EDOT_J*POWER_UNIT_ERG_PER_S, EDOT_S*POWER_UNIT_ERG_PER_S, RESIDUAL*POWER_UNIT_ERG_PER_S

    PREV_E_POL = E_POL
    PREV_E_TOR = E_TOR
  END SELECT
END SUBROUTINE LOG_ENERGY

END PROGRAM MHDVSH_DIFFUSION
