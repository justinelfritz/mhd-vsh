!> TOV stellar-structure driver: solves a single neutron star's radial
!> profile (SOLVE_TOV_STAR) at a fixed central density, then computes a
!> self-consistent resistivity/Hall-prefactor profile from it
!> (ETA_AND_F_HALL_AT), printing the M/R summary and (optionally) writing
!> the full radial profile to a data file.
!>
!> Central density (9.88e14 `g/cm**3`) matches nstot.f's own commented
!> "M=1.40" reference star -- the regression target for this port (see
!> test/test_tov_solver.f90, checked against the real, already-solved
!> 2015 output at ~/Desktop/EOSNS/fort.34 and PL.DAT).
!>
!> Command-line arguments (all optional):
!>   1: EOS table path (default data/eos/lowd-eos.ja.tab)
!>   2: output data file path for the radial profile (r, rho, n_e,
!>      eta, f_H, all in this project's own code units -- see UNITS)
PROGRAM MHDVSH_TOV
USE KINDS,              ONLY: dp, i4
USE TOV_SOLVER,         ONLY: TOV_PROFILE_T, SOLVE_TOV_STAR
USE CRUST_CONDUCTIVITY, ONLY: ETA_AND_F_HALL_AT
IMPLICIT NONE

REAL(KIND=dp),    PARAMETER :: RHOCGS  = 9.88E14_dp   ! M=1.40 reference star
INTEGER(KIND=i4), PARAMETER :: NPOINTS = 414
INTEGER(KIND=i4), PARAMETER :: DATA_UNIT = 22

TYPE(TOV_PROFILE_T) :: PROFILE
REAL(KIND=dp), ALLOCATABLE :: ETA_PROFILE(:), F_HALL_PROFILE(:)
REAL(KIND=dp), ALLOCATABLE :: R_CRUST(:), RHOCGS_CRUST(:), NEL_CGS_CRUST(:)
CHARACTER(LEN=1024) :: EOS_PATH, DATA_PATH
LOGICAL :: WRITE_DATA
INTEGER(KIND=i4) :: I, N_CRUST

IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
  CALL GET_COMMAND_ARGUMENT(1, EOS_PATH)
ELSE
  EOS_PATH = 'data/eos/lowd-eos.ja.tab'
END IF

CALL SOLVE_TOV_STAR(RHOCGS, TRIM(EOS_PATH), NPOINTS, PROFILE)

WRITE(*,'(A)')          'MHD-VSH TOV solver'
WRITE(*,'(A,A)')        '  EOS table    = ', TRIM(EOS_PATH)
WRITE(*,'(A,ES12.5,A)') '  rho_central  = ', RHOCGS, ' g/cm**3'
WRITE(*,'(A,ES12.5,A)') '  Radius       = ', PROFILE%RADIUS_KM, ' km'
WRITE(*,'(A,ES12.5,A)') '  Mass         = ', PROFILE%MASS_MSUN, ' Msun'

! Crust extent: rows with a genuine nucleus (A_TABLE>0) are where
! ETA_AND_F_HALL_AT's own gate accepts input -- this is the same
! R_MIN/R_MAX range app/mhdvsh_hall.f90 should be updated to use.
N_CRUST = COUNT(PROFILE%A_TABLE > 0.0_dp)
IF (N_CRUST > 0) THEN
  WRITE(*,'(A,ES12.5,A,ES12.5,A)') '  Crust extent = ', &
    MINVAL(PROFILE%R, MASK=PROFILE%A_TABLE > 0.0_dp), ' -- ', &
    MAXVAL(PROFILE%R, MASK=PROFILE%A_TABLE > 0.0_dp), ' km'
END IF

! ETA_AND_F_HALL_AT only accepts rows with a genuine nucleus (A_TABLE>0)
! -- restrict to that sub-range rather than passing the full profile.
! R_CRUST/RHOCGS_CRUST/NEL_CGS_CRUST come back index-aligned with
! ETA_PROFILE/F_HALL_PROFILE (all sized N_CRUST), unlike PROFILE's own
! arrays (sized N, center-to-surface).
CALL ETA_AND_F_HALL_AT_CRUST(PROFILE, ETA_PROFILE, F_HALL_PROFILE, &
  R_CRUST, RHOCGS_CRUST, NEL_CGS_CRUST)

WRITE_DATA = (COMMAND_ARGUMENT_COUNT() >= 2)
IF (WRITE_DATA) THEN
  CALL GET_COMMAND_ARGUMENT(2, DATA_PATH)
  OPEN(UNIT=DATA_UNIT, FILE=TRIM(DATA_PATH), STATUS='REPLACE', ACTION='WRITE')
  WRITE(DATA_UNIT,'(A)') '# r_km  rhocgs  n_e_cm-3  eta_km2_per_yr  f_hall_km2_per_1e12G_per_yr'
  DO I = 1, N_CRUST
    WRITE(DATA_UNIT,'(5ES16.8)') R_CRUST(I), RHOCGS_CRUST(I), NEL_CGS_CRUST(I), &
      ETA_PROFILE(I), F_HALL_PROFILE(I)
  END DO
  CLOSE(DATA_UNIT)
  WRITE(*,'(A,A)') '  profile written to ', TRIM(DATA_PATH)
END IF

CONTAINS

!> ETA_AND_F_HALL_AT itself STOPs on any row without a genuine nucleus
!> (by design -- see its own header), so build a crust-only sub-profile
!> first rather than passing PROFILE's full center-to-surface extent
!> (which includes homogeneous-matter rows near the center). Also
!> returns the crust-only R/RHOCGS/n_e(cgs) arrays, index-aligned with
!> ETA_OUT/F_HALL_OUT, since PROFILE's own R(1:N_CRUST) does NOT
!> correspond to the crust rows (the crust sits near the surface, not
!> at the start of the center-to-surface profile).
SUBROUTINE ETA_AND_F_HALL_AT_CRUST(FULL_PROFILE, ETA_OUT, F_HALL_OUT, R_OUT, RHOCGS_OUT, NEL_CGS_OUT)
  TYPE(TOV_PROFILE_T), INTENT(IN) :: FULL_PROFILE
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: ETA_OUT(:), F_HALL_OUT(:)
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: R_OUT(:), RHOCGS_OUT(:), NEL_CGS_OUT(:)
  TYPE(TOV_PROFILE_T) :: CRUST
  INTEGER(KIND=i4) :: J, K

  CRUST%N = COUNT(FULL_PROFILE%A_TABLE > 0.0_dp)
  ALLOCATE(CRUST%R(CRUST%N), CRUST%RHOCGS(CRUST%N), CRUST%PCGS(CRUST%N))
  ALLOCATE(CRUST%NBFM(CRUST%N), CRUST%NEL(CRUST%N))
  ALLOCATE(CRUST%AH(CRUST%N), CRUST%ZH(CRUST%N), CRUST%XH(CRUST%N))
  ALLOCATE(CRUST%YE(CRUST%N), CRUST%YN(CRUST%N), CRUST%A_TABLE(CRUST%N))

  K = 0
  DO J = 1, FULL_PROFILE%N
    IF (FULL_PROFILE%A_TABLE(J) <= 0.0_dp) CYCLE
    K = K + 1
    CRUST%R(K)      = FULL_PROFILE%R(J)
    CRUST%RHOCGS(K) = FULL_PROFILE%RHOCGS(J)
    CRUST%PCGS(K)   = FULL_PROFILE%PCGS(J)
    CRUST%NBFM(K)   = FULL_PROFILE%NBFM(J)
    CRUST%NEL(K)    = FULL_PROFILE%NEL(J)
    CRUST%AH(K)     = FULL_PROFILE%AH(J)
    CRUST%ZH(K)     = FULL_PROFILE%ZH(J)
    CRUST%XH(K)     = FULL_PROFILE%XH(J)
    CRUST%YE(K)     = FULL_PROFILE%YE(J)
    CRUST%YN(K)     = FULL_PROFILE%YN(J)
    CRUST%A_TABLE(K)= FULL_PROFILE%A_TABLE(J)
  END DO

  CALL ETA_AND_F_HALL_AT(CRUST, ETA_OUT, F_HALL_OUT)
  R_OUT = CRUST%R
  RHOCGS_OUT = CRUST%RHOCGS
  NEL_CGS_OUT = CRUST%NEL * 1.0E39_dp   ! fm**-3 -> cm**-3
END SUBROUTINE ETA_AND_F_HALL_AT_CRUST

END PROGRAM MHDVSH_TOV
