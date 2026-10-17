!> One-off diagnostic (2026-08-29): dumps Phi(r) near the boundary from a
!> REAL, already-completed simulation's own checkpoint file (not a fresh
!> IC reconstruction -- see mhdvsh_ic_boundary_check.f90 for that) for
!> every mode with non-negligible amplitude, plus the FD-operator
!> (OPS%D1, the actual discretized derivative the code uses) evaluation
!> of the Robin-BC residual at the outer boundary for each such mode.
!>
!> Usage: mhdvsh_checkpoint_boundary_check <checkpoint_path> <output_data_path>
PROGRAM MHDVSH_CHECKPOINT_BOUNDARY_CHECK
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,   ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T
USE IO_CHECKPOINT,      ONLY: READ_CHECKPOINT
IMPLICIT NONE

REAL(KIND=dp), PARAMETER :: ACTIVE_TOL = 1.0E-6_dp  ! relative to the mode's own peak |Phi|

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(SPECTRAL_SCALAR_T) :: PHI, PSI
REAL(KIND=dp) :: R_MIN, R_MAX, T_CK, LAMBDA_L, DPHI_DR_BND, ROBIN_REQUIRED, ROBIN_RESIDUAL
REAL(KIND=dp) :: PEAK_ABS
INTEGER(KIND=i4) :: N_R, LMAX, ISTEP_CK
INTEGER(KIND=i4) :: L, M, IDX, IR, UNIT
CHARACTER(LEN=1024) :: CHECKPOINT_PATH, OUT_PATH

IF (COMMAND_ARGUMENT_COUNT() < 2) THEN
  WRITE(*,'(A)') 'Usage: mhdvsh_checkpoint_boundary_check <checkpoint_path> <output_data_path>'
  STOP 1
END IF
CALL GET_COMMAND_ARGUMENT(1, CHECKPOINT_PATH)
CALL GET_COMMAND_ARGUMENT(2, OUT_PATH)

CALL READ_CHECKPOINT(TRIM(CHECKPOINT_PATH), ISTEP_CK, T_CK, PHI, PSI, N_R, LMAX, R_MIN, R_MAX)
WRITE(*,'(A,I0,A,ES14.6,A,I0,A,I0)') 'checkpoint: step=', ISTEP_CK, ' t=', T_CK, &
  ' N_R=', N_R, ' LMAX=', LMAX

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)

OPEN(NEWUNIT=UNIT, FILE=TRIM(OUT_PATH), STATUS='REPLACE', ACTION='WRITE')
WRITE(UNIT,'(A)') '# l  m  ir  r_km  Re_Phi(r)  Im_Phi(r)  B_r(r)=l(l+1)/r**2*Re_Phi(r)'

WRITE(*,'(A)') 'Robin-BC residual (real FD operator) by mode, active modes only:'
WRITE(*,'(A)') '   l   m       Phi(Rout)         dPhi/dr|Rout    required(-l/R*Phi)       residual'

DO L = 0, LMAX
  LAMBDA_L = REAL(L*(L+1), KIND=dp)
  DO M = -L, L
    IDX = YLM_INDEX(L, M)
    PEAK_ABS = MAXVAL(ABS(PHI%COEF(:,IDX)))
    IF (PEAK_ABS < ACTIVE_TOL * MAXVAL(ABS(PHI%COEF))) CYCLE

    DO IR = 1, N_R
      WRITE(UNIT,'(2I5,I8,4ES18.10)') L, M, IR, RGRID%R(IR), &
        REAL(PHI%COEF(IR,IDX),KIND=dp), AIMAG(PHI%COEF(IR,IDX)), &
        LAMBDA_L/RGRID%R(IR)**2 * REAL(PHI%COEF(IR,IDX),KIND=dp)
    END DO

    DPHI_DR_BND = SUM(OPS%D1(N_R,:) * REAL(PHI%COEF(:,IDX), KIND=dp))
    ROBIN_REQUIRED = -(REAL(L,KIND=dp)/RGRID%R(N_R)) * REAL(PHI%COEF(N_R,IDX),KIND=dp)
    ROBIN_RESIDUAL = DPHI_DR_BND - ROBIN_REQUIRED
    WRITE(*,'(2I4,4ES18.8)') L, M, REAL(PHI%COEF(N_R,IDX),KIND=dp), DPHI_DR_BND, &
      ROBIN_REQUIRED, ROBIN_RESIDUAL
  END DO
END DO
CLOSE(UNIT)
WRITE(*,'(A,A)') 'full per-mode Phi(r) profiles written to ', TRIM(OUT_PATH)

END PROGRAM MHDVSH_CHECKPOINT_BOUNDARY_CHECK
