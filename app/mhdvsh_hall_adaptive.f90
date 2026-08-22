!> Combined resistive+Hall regime driver with a DYNAMIC, Hall-CFL-limited
!> outer timestep, instead of app/mhdvsh_hall.f90's fixed DT -- otherwise
!> identical setup (same grid, seed mode, ETA/F_HALL/N_SUB). Built on
!> TIMESTEPPER::RUN_ADAPTIVE + HALL_REGIME::HALL_COMPUTE_DT/HALL_SET_DT
!> (see their own docstrings for the exact tc=min(dr/(F_HALL*|current|))
!> formula, Justin Elfritz, 2026-08-21, and the RMS-vs-pointwise-max
!> caveat in how |current| is estimated).
!>
!> DT is recomputed every N_SUB outer steps (amortizing the cost of
!> re-factorizing DIFFUSION_INIT's implicit solve, which HALL_SET_DT
!> triggers on every recompute) and set to 0.5*tc (a safety margin below
!> the raw Hall-CFL limit) -- both per Justin Elfritz's own choice
!> (2026-08-21), not hardcoded defaults; see mhdvsh_hall.f90 for the
!> fixed-DT counterpart these choices are compared against.
!>
!> Optional command-line argument: a file path to write a time-resolved
!> series of (step, t, dt_yr, tc_raw_yr, E_poloidal, E_toroidal) --
!> tc_raw is the un-scaled Hall-CFL limit (dt=0.5*tc_raw), logged
!> alongside dt itself so a growing gap between them (tc_raw shrinking
!> as the Hall term drives a current cascade) is visible directly, not
!> just inferred from dt's own trend.
!>
!> Units: same convention as mhdvsh_hall.f90 (UNITS::ENERGY_UNIT_ERG
!> applied only at this reporting layer; dt/tc are already in yr, this
!> project's code time unit -- see UNITS).
PROGRAM MHDVSH_HALL_ADAPTIVE
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,   ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS,  ONLY: TOTAL_MAGNETIC_ENERGY, TOTAL_POLOIDAL_MAGNETIC_ENERGY, &
                               TOTAL_TOROIDAL_MAGNETIC_ENERGY, HALL_COURANT_TIMESTEP
USE DIFFUSION_REGIME,   ONLY: DIFFUSION_STATE_T
USE HALL_REGIME,        ONLY: HALL_INIT, HALL_ADVANCE, HALL_COMPUTE_DT, HALL_SET_DT
USE TIMESTEPPER,        ONLY: RUN_ADAPTIVE
USE UNITS,              ONLY: ENERGY_UNIT_ERG
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R     = 40
! Same crust extent as mhdvsh_hall.f90 (M=1.40 reference star, see its
! own comment for provenance).
REAL(KIND=dp),    PARAMETER :: R_MIN   = 10.8033325018_dp   ! km
REAL(KIND=dp),    PARAMETER :: R_MAX   = 11.6982211606_dp   ! km
INTEGER(KIND=i4), PARAMETER :: LMAX    = 30
REAL(KIND=dp),    PARAMETER :: ETA     = 1.0E-6_dp  ! km**2/yr
REAL(KIND=dp),    PARAMETER :: F_HALL  = 0.01_dp    ! km**2/(1e12 G)/yr
INTEGER(KIND=i4), PARAMETER :: N_SUB   = 10          ! Hall substeps per outer step
INTEGER(KIND=i4), PARAMETER :: N_STEPS = 200
REAL(KIND=dp),    PARAMETER :: CFL_SAFETY = 0.5_dp   ! dt = CFL_SAFETY*tc
INTEGER(KIND=i4), PARAMETER :: DATA_UNIT = 21
REAL(KIND=dp),    PARAMETER :: PI = 3.14159265358979_dp

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(DIFFUSION_STATE_T) :: STATE
REAL(KIND=dp) :: E0, E1, TC0, DT0
INTEGER(KIND=i4) :: IR
CHARACTER(LEN=1024) :: DATA_PATH
LOGICAL :: WRITE_DATA

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)

CALL ALLOC_SPECTRAL_SCALAR(STATE%PHI, N_R, LMAX)
CALL ALLOC_SPECTRAL_SCALAR(STATE%PSI, N_R, LMAX)

! Same single seed mode as mhdvsh_hall.f90.
DO IR = 1, N_R
  STATE%PHI%COEF(IR, YLM_INDEX(1_i4,0_i4)) = &
    CMPLX(SIN(PI*(RGRID%R(IR)-R_MIN)/(R_MAX-R_MIN)), 0.0_dp, KIND=dp)
END DO

! Initial DT from the seed field's own Hall-CFL limit, not an arbitrary
! placeholder -- RUN_ADAPTIVE's own SET_DT call on step 1 (see its
! docstring) would overwrite an arbitrary guess anyway, but there's no
! reason to pass one when the true value is this cheap to compute
! before HALL_INIT (HALL_COURANT_TIMESTEP only needs OPS/RGRID, both
! already built above).
TC0 = HALL_COURANT_TIMESTEP(STATE%PHI, STATE%PSI, OPS, RGRID, F_HALL)
DT0 = CFL_SAFETY * TC0
CALL HALL_INIT(RGRID, OPS, LMAX, ETA, F_HALL, DT0, N_SUB)

WRITE_DATA = (COMMAND_ARGUMENT_COUNT() >= 1)
IF (WRITE_DATA) THEN
  CALL GET_COMMAND_ARGUMENT(1, DATA_PATH)
  OPEN(UNIT=DATA_UNIT, FILE=TRIM(DATA_PATH), STATUS='REPLACE', ACTION='WRITE')
  WRITE(DATA_UNIT,'(A)') '# step  t_yr  dt_yr  tc_raw_yr  E_poloidal_erg  E_toroidal_erg'
  WRITE(DATA_UNIT,'(I8,5ES16.8)') 0_i4, 0.0_dp, DT0, TC0, &
    TOTAL_POLOIDAL_MAGNETIC_ENERGY(STATE%PHI, OPS, RGRID)*ENERGY_UNIT_ERG, &
    TOTAL_TOROIDAL_MAGNETIC_ENERGY(STATE%PSI, RGRID)*ENERGY_UNIT_ERG
END IF

E0 = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)
IF (WRITE_DATA) THEN
  CALL RUN_ADAPTIVE(HALL_ADVANCE, STATE, HALL_COMPUTE_DT, HALL_SET_DT, N_STEPS, &
    DT_RECOMPUTE_EVERY=N_SUB, CFL_SAFETY=CFL_SAFETY, ON_STEP=LOG_STATE)
ELSE
  CALL RUN_ADAPTIVE(HALL_ADVANCE, STATE, HALL_COMPUTE_DT, HALL_SET_DT, N_STEPS, &
    DT_RECOMPUTE_EVERY=N_SUB, CFL_SAFETY=CFL_SAFETY)
END IF
E1 = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)

IF (WRITE_DATA) CLOSE(DATA_UNIT)

WRITE(*,'(A)')            'MHD-VSH combined resistive+Hall regime -- dynamic Hall-CFL dt'
WRITE(*,'(A,I0)')         '  N_r          = ', N_R
WRITE(*,'(A,I0)')         '  Lmax         = ', LMAX
WRITE(*,'(A,ES12.5,A)')   '  eta          = ', ETA, ' km**2/yr'
WRITE(*,'(A,ES12.5,A)')   '  F_Hall       = ', F_HALL, ' km**2/(1e12 G)/yr'
WRITE(*,'(A,ES12.5,A,ES12.5,A)') '  initial dt   = ', DT0, ' yr (CFL_SAFETY*tc, tc = ', TC0, ' yr)'
WRITE(*,'(A,I0,A,I0,A)')  '  N_steps      = ', N_STEPS, ' (dt recomputed every ', N_SUB, ' steps)'
WRITE(*,'(A,ES12.5,A)')   '  energy(t=0)  = ', E0*ENERGY_UNIT_ERG, ' erg'
WRITE(*,'(A,ES12.5,A)')   '  energy(t=end)= ', E1*ENERGY_UNIT_ERG, ' erg'
IF (E1 < E0) THEN
  WRITE(*,'(A)') '  RESULT: energy decreased (diffusion dissipating, as expected)'
ELSE
  WRITE(*,'(A)') '  RESULT: UNEXPECTED -- energy did not decrease'
END IF
IF (WRITE_DATA) WRITE(*,'(A,A)') '  time series written to ', TRIM(DATA_PATH)

CONTAINS

!> Matches TIMESTEPPER::REGIME_ON_STEP_I. Logs (step, t, dt, tc_raw,
!> E_poloidal, E_toroidal). dt/tc_raw are recomputed fresh from the
!> post-step STATE via HALL_COMPUTE_DT every logged step -- on a step
!> where RUN_ADAPTIVE itself also recomputed DT (every N_SUB steps),
!> this matches the DT actually used going forward exactly; in between,
!> it's the CFL limit implied by the state as it now stands, generally
!> different from (typically smaller than, as the field evolves) the DT
!> actually in use for the next few steps until the next recompute.
SUBROUTINE LOG_STATE(STATE, T, ISTEP)
  CLASS(*),         INTENT(IN) :: STATE
  REAL(KIND=dp),    INTENT(IN) :: T
  INTEGER(KIND=i4), INTENT(IN) :: ISTEP
  REAL(KIND=dp) :: TC_NOW, DT_NOW
  SELECT TYPE (STATE)
  TYPE IS (DIFFUSION_STATE_T)
    TC_NOW = HALL_COMPUTE_DT(STATE)
    DT_NOW = CFL_SAFETY * TC_NOW
    WRITE(DATA_UNIT,'(I8,5ES16.8)') ISTEP, T, DT_NOW, TC_NOW, &
      TOTAL_POLOIDAL_MAGNETIC_ENERGY(STATE%PHI, OPS, RGRID)*ENERGY_UNIT_ERG, &
      TOTAL_TOROIDAL_MAGNETIC_ENERGY(STATE%PSI, RGRID)*ENERGY_UNIT_ERG
  END SELECT
END SUBROUTINE LOG_STATE

END PROGRAM MHDVSH_HALL_ADAPTIVE
