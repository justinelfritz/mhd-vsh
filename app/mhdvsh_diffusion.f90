!> Diffusion-regime driver: seeds a smooth (Phi,Psi) profile on a shell
!> grid, evolves it under pure Ohmic diffusion via DIFFUSION_REGIME +
!> TIMESTEPPER, and reports the magnetic energy before/after as a
!> sanity check (pure diffusion must dissipate, never create, energy).
!>
!> Optional command-line argument: a file path to write a time-resolved
!> (t, total_energy) series to, one row per step, for
!> scripts/plot_energy_budget.py to read. Only total magnetic energy is
!> logged -- the individual Joule/Poynting balance terms aren't derived
!> in the code yet (see FIELD_DIAGNOSTICS's header), so there's nothing
!> else to log until that physics exists.
PROGRAM MHDVSH_DIFFUSION
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,   ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS,  ONLY: TOTAL_MAGNETIC_ENERGY
USE DIFFUSION_REGIME,   ONLY: DIFFUSION_STATE_T, DIFFUSION_INIT, DIFFUSION_ADVANCE
USE TIMESTEPPER,        ONLY: RUN
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R     = 40
REAL(KIND=dp),    PARAMETER :: R_MIN   = 0.5_dp
REAL(KIND=dp),    PARAMETER :: R_MAX   = 1.0_dp
INTEGER(KIND=i4), PARAMETER :: LMAX    = 3
REAL(KIND=dp),    PARAMETER :: ETA     = 0.05_dp
REAL(KIND=dp),    PARAMETER :: DT      = 0.01_dp
INTEGER(KIND=i4), PARAMETER :: N_STEPS = 200
INTEGER(KIND=i4), PARAMETER :: DATA_UNIT = 21

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(DIFFUSION_STATE_T) :: STATE
REAL(KIND=dp) :: E0, E1
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
  WRITE(DATA_UNIT,'(A)') '# step  t  total_magnetic_energy'
  WRITE(DATA_UNIT,'(I8,2ES16.8)') 0, 0.0_dp, &
    TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)
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
WRITE(*,'(A,ES12.5)')     '  eta          = ', ETA
WRITE(*,'(A,ES12.5,A,I0,A)') '  dt           = ', DT, '  (', N_STEPS, ' steps)'
WRITE(*,'(A,ES12.5)')     '  energy(t=0)  = ', E0
WRITE(*,'(A,ES12.5)')     '  energy(t=end)= ', E1
IF (E1 < E0) THEN
  WRITE(*,'(A)') '  RESULT: energy decreased, as expected for pure diffusion'
ELSE
  WRITE(*,'(A)') '  RESULT: UNEXPECTED -- energy did not decrease'
END IF
IF (WRITE_DATA) WRITE(*,'(A,A)') '  time series written to ', TRIM(DATA_PATH)

CONTAINS

!> Matches TIMESTEPPER::REGIME_ON_STEP_I; logs (step, t, total energy)
!> to the already-open DATA_UNIT. Host-associates OPS/RGRID from the
!> main program rather than needing module-level state, since this
!> procedure is only ever used as an actual argument within this one run.
SUBROUTINE LOG_ENERGY(STATE, T, ISTEP)
  CLASS(*),         INTENT(IN) :: STATE
  REAL(KIND=dp),    INTENT(IN) :: T
  INTEGER(KIND=i4), INTENT(IN) :: ISTEP
  SELECT TYPE (STATE)
  TYPE IS (DIFFUSION_STATE_T)
    WRITE(DATA_UNIT,'(I8,2ES16.8)') ISTEP, T, &
      TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)
  END SELECT
END SUBROUTINE LOG_ENERGY

END PROGRAM MHDVSH_DIFFUSION
