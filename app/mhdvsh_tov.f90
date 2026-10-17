!> TOV stellar-structure driver: solves a single neutron star's radial
!> profile (SOLVE_TOV_STAR) at a fixed central density, then computes a
!> self-consistent resistivity/Hall-prefactor profile from it
!> (ETA_AND_F_HALL_AT), printing the M/R summary and (optionally) writing
!> the full radial profile to a data file.
!>
!> Central baryon number density (0.5447307 `fm**-3`) matches NSCool's own
!> bundled "M=1.40" reference star -- the regression target for this
!> port (see test/test_tov_solver.f90, checked against NSCool's own
!> bundled `Prof_APR_Cat_1.4.dat`/`prod_APR_EOS_Cat.dat`).
!>
!> Command-line arguments (all optional):
!>   1: EOS table path (default data/eos/APR_EOS_Cat.dat)
!>   2: output data file path for the radial profile (r, rho, n_e,
!>      eta, f_H, all in this project's own code units -- see UNITS)
!>   3: N_R, number of nodes for an equidistributed simulation radial
!>      grid (skipped entirely if absent -- this is a genuinely free
!>      parameter, independent of NPOINTS above, which only controls
!>      TOV-integration accuracy)
!>   4: output path for the equidistributed grid (default
!>      data/grid/radial_grid.dat), only used if (3) is given
!>
!> Equidistributed grid (2026-10-16, per user request): places N_R
!> nodes so that a monitor function M(r) = |d ln(eta)/dr| + |d
!> ln(f_H)/dr| -- the two quantities that directly set the Hall-CFL
!> constraint -- integrates to equal increments between consecutive
!> nodes (cf. grid_radial.f90's own header for the motivating ~67 vs
!> ~8e5 Rmag gap that a uniform grid produces). A small floor is added
!> to M before integrating so no region (e.g. deep in the crust, where
!> eta/f_H vary slowly) is left with zero node density.
PROGRAM MHDVSH_TOV
USE KINDS,              ONLY: dp, i4
USE TOV_SOLVER,         ONLY: TOV_PROFILE_T, SOLVE_TOV_STAR
USE CRUST_CONDUCTIVITY, ONLY: ETA_AND_F_HALL_AT
USE ODE_INTEGRATOR,     ONLY: LOCATE_TABLE
IMPLICIT NONE

REAL(KIND=dp),    PARAMETER :: NBAR_CENTRAL = 0.5447307_dp   ! fm**-3, M=1.40 reference star
INTEGER(KIND=i4), PARAMETER :: NPOINTS = 414
INTEGER(KIND=i4), PARAMETER :: DATA_UNIT = 22
INTEGER(KIND=i4), PARAMETER :: GRID_UNIT = 23
REAL(KIND=dp),    PARAMETER :: RHOL_CGS = 2.2E14_dp   ! core-crust boundary, matches CRUST_CONDUCTIVITY
REAL(KIND=dp),    PARAMETER :: MONITOR_FLOOR_FRAC = 0.05_dp   ! floor = this * mean(M)

TYPE(TOV_PROFILE_T) :: PROFILE
REAL(KIND=dp), ALLOCATABLE :: ETA_PROFILE(:), F_HALL_PROFILE(:), N_E_PROFILE(:)
REAL(KIND=dp), ALLOCATABLE :: R_CRUST(:), RHOCGS_CRUST(:)
REAL(KIND=dp), ALLOCATABLE :: R_GRID(:), DR_GRID(:)
CHARACTER(LEN=1024) :: EOS_PATH, DATA_PATH, GRID_PATH, NR_ARG
LOGICAL :: WRITE_DATA, WRITE_GRID
INTEGER(KIND=i4) :: I, N_CRUST, N_R

IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
  CALL GET_COMMAND_ARGUMENT(1, EOS_PATH)
ELSE
  EOS_PATH = 'data/eos/APR_EOS_Cat.dat'
END IF

CALL SOLVE_TOV_STAR(NBAR_CENTRAL, TRIM(EOS_PATH), NPOINTS, PROFILE)

WRITE(*,'(A)')          'MHD-VSH TOV solver'
WRITE(*,'(A,A)')        '  EOS table    = ', TRIM(EOS_PATH)
WRITE(*,'(A,ES12.5,A)') '  nbar_central = ', NBAR_CENTRAL, ' fm**-3'
WRITE(*,'(A,ES12.5,A)') '  Radius       = ', PROFILE%RADIUS_KM, ' km'
WRITE(*,'(A,ES12.5,A)') '  Mass         = ', PROFILE%MASS_MSUN, ' Msun'

! Crust extent: rows below the core-crust boundary density are where
! ETA_AND_F_HALL_AT's own gate accepts input -- this is the same
! R_MIN/R_MAX range app/mhdvsh_hall.f90 should be updated to use.
N_CRUST = COUNT(PROFILE%RHOCGS <= RHOL_CGS)
IF (N_CRUST > 0) THEN
  WRITE(*,'(A,ES12.5,A,ES12.5,A)') '  Crust extent = ', &
    MINVAL(PROFILE%R, MASK=PROFILE%RHOCGS <= RHOL_CGS), ' -- ', &
    MAXVAL(PROFILE%R, MASK=PROFILE%RHOCGS <= RHOL_CGS), ' km'
END IF

! ETA_AND_F_HALL_AT only accepts rows below the core-crust boundary --
! restrict to that sub-range rather than passing the full profile.
! R_CRUST/RHOCGS_CRUST come back index-aligned with ETA_PROFILE/
! F_HALL_PROFILE/N_E_PROFILE (all sized N_CRUST), unlike PROFILE's own
! arrays (sized N, center-to-surface).
CALL ETA_AND_F_HALL_AT_CRUST(PROFILE, ETA_PROFILE, F_HALL_PROFILE, N_E_PROFILE, &
  R_CRUST, RHOCGS_CRUST)

WRITE_DATA = (COMMAND_ARGUMENT_COUNT() >= 2)
IF (WRITE_DATA) THEN
  CALL GET_COMMAND_ARGUMENT(2, DATA_PATH)
  OPEN(UNIT=DATA_UNIT, FILE=TRIM(DATA_PATH), STATUS='REPLACE', ACTION='WRITE')
  WRITE(DATA_UNIT,'(A)') '# r_km  rhocgs  n_e_cm-3  eta_km2_per_yr  f_hall_km2_per_1e12G_per_yr'
  DO I = 1, N_CRUST
    WRITE(DATA_UNIT,'(5ES16.8)') R_CRUST(I), RHOCGS_CRUST(I), N_E_PROFILE(I), &
      ETA_PROFILE(I), F_HALL_PROFILE(I)
  END DO
  CLOSE(DATA_UNIT)
  WRITE(*,'(A,A)') '  profile written to ', TRIM(DATA_PATH)
END IF

WRITE_GRID = (COMMAND_ARGUMENT_COUNT() >= 3)
IF (WRITE_GRID) THEN
  CALL GET_COMMAND_ARGUMENT(3, NR_ARG)
  READ(NR_ARG, *) N_R

  IF (COMMAND_ARGUMENT_COUNT() >= 4) THEN
    CALL GET_COMMAND_ARGUMENT(4, GRID_PATH)
  ELSE
    GRID_PATH = 'data/grid/radial_grid.dat'
  END IF

  CALL BUILD_EQUIDISTRIBUTED_GRID(R_CRUST, ETA_PROFILE, F_HALL_PROFILE, N_R, R_GRID, DR_GRID)

  OPEN(UNIT=GRID_UNIT, FILE=TRIM(GRID_PATH), STATUS='REPLACE', ACTION='WRITE')
  WRITE(GRID_UNIT,'(A)') '# i  r_km  dr_km'
  DO I = 1, N_R
    WRITE(GRID_UNIT,'(I6,2ES16.8)') I, R_GRID(I), DR_GRID(I)
  END DO
  CLOSE(GRID_UNIT)
  WRITE(*,'(A,I0,A,A)') '  equidistributed grid (N_R=', N_R, ') written to ', TRIM(GRID_PATH)
END IF

CONTAINS

!> ETA_AND_F_HALL_AT itself STOPs on any row above the core-crust
!> boundary density (by design -- see its own header), so build a
!> crust-only sub-profile first rather than passing PROFILE's full
!> center-to-surface extent (which includes core rows near the
!> center). Also returns the crust-only R/RHOCGS arrays, index-aligned
!> with ETA_OUT/F_HALL_OUT/N_E_OUT, since PROFILE's own R(1:N_CRUST)
!> does NOT correspond to the crust rows (the crust sits near the
!> surface, not at the start of the center-to-surface profile).
SUBROUTINE ETA_AND_F_HALL_AT_CRUST(FULL_PROFILE, ETA_OUT, F_HALL_OUT, N_E_OUT, R_OUT, RHOCGS_OUT)
  TYPE(TOV_PROFILE_T), INTENT(IN) :: FULL_PROFILE
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: ETA_OUT(:), F_HALL_OUT(:), N_E_OUT(:)
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: R_OUT(:), RHOCGS_OUT(:)
  TYPE(TOV_PROFILE_T) :: CRUST
  INTEGER(KIND=i4) :: J, K

  CRUST%N = COUNT(FULL_PROFILE%RHOCGS <= RHOL_CGS)
  ALLOCATE(CRUST%R(CRUST%N), CRUST%RHOCGS(CRUST%N), CRUST%PCGS(CRUST%N), CRUST%NBFM(CRUST%N))

  K = 0
  DO J = 1, FULL_PROFILE%N
    IF (FULL_PROFILE%RHOCGS(J) > RHOL_CGS) CYCLE
    K = K + 1
    CRUST%R(K)      = FULL_PROFILE%R(J)
    CRUST%RHOCGS(K) = FULL_PROFILE%RHOCGS(J)
    CRUST%PCGS(K)   = FULL_PROFILE%PCGS(J)
    CRUST%NBFM(K)   = FULL_PROFILE%NBFM(J)
  END DO

  CALL ETA_AND_F_HALL_AT(CRUST, ETA_OUT, F_HALL_OUT, N_E_OUT)
  R_OUT = CRUST%R
  RHOCGS_OUT = CRUST%RHOCGS
END SUBROUTINE ETA_AND_F_HALL_AT_CRUST

!> Places N_R nodes on [R_DENSE(1), R_DENSE(SIZE(R_DENSE))] by
!> equidistributing the monitor function M(r) = |d ln(eta)/dr| + |d
!> ln(f_H)/dr|, built from the TOV-crust profile R_DENSE
!> (ETA_DENSE/F_HALL_DENSE index-aligned with it). A floor
!> MONITOR_FLOOR_FRAC*mean(M) is added to M before integrating so flat
!> regions still get a baseline node density rather than none.
!>
!> R_DENSE is NOT differentiated directly: the TOV solver's own
!> log(pressure)-based adaptive stepping places dozens of points
!> within ~1e-12 km of each other right at the surface (where
!> pressure/density fall off fastest), so naive adjacent-point
!> differencing there divides by a near-zero spacing and produces a
!> single enormous spike that swamps MEAN_M and flattens the whole
!> grid (found by direct inspection of R_DENSE's spacing while
!> building this). Instead, ETA_DENSE/F_HALL_DENSE are first
!> log-linearly resampled (via LOCATE_TABLE bisection) onto a uniform
!> table of N_SAMPLE points spanning the same [R_MIN,R_MAX] -- this
!> decouples "where the TOV integrator's pressure-based stepping
!> happened to land" from "where the monitor function is evaluated",
!> and the monitor function is then well-behaved (uniform spacing, no
!> near-zero denominators).
!>
!> Algorithm (standard equidistribution meshing), on the uniform
!> resampled table: compute M via centered (one-sided at the ends)
!> finite differences of ln(eta), ln(f_H); integrate M by the
!> trapezoidal rule to get the cumulative weight S(r), S(R_MIN)=0;
!> normalize to S_HAT in [0,1]; place node k (k=1..N_R) at the r where
!> S_HAT(r) = (k-1)/(N_R-1), found by linear inversion against the
!> resampled table (S_HAT is strictly increasing since M>0 everywhere
!> after the floor). Endpoints are copied exactly from R_DENSE, not
!> interpolated, to avoid roundoff drift off its own boundary values.
!>
!> DR_OUT is the centered cell width: DR_OUT(i) =
!> (R_OUT(i+1)-R_OUT(i-1))/2 for interior i, and the one-sided gap to
!> the single neighbor at i=1 and i=N_R (confirmed convention).
SUBROUTINE BUILD_EQUIDISTRIBUTED_GRID(R_DENSE, ETA_DENSE, F_HALL_DENSE, N_R, R_OUT, DR_OUT)
  REAL(KIND=dp),    INTENT(IN)  :: R_DENSE(:), ETA_DENSE(:), F_HALL_DENSE(:)
  INTEGER(KIND=i4), INTENT(IN)  :: N_R
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: R_OUT(:), DR_OUT(:)

  INTEGER(KIND=i4), PARAMETER :: N_SAMPLE = 2000_i4
  INTEGER(KIND=i4) :: N_DENSE, J, K
  REAL(KIND=dp), ALLOCATABLE :: LOGETA_DENSE(:), LOGFH_DENSE(:)
  REAL(KIND=dp) :: R_UNIFORM(N_SAMPLE), LOGETA(N_SAMPLE), LOGFH(N_SAMPLE)
  REAL(KIND=dp) :: MONITOR(N_SAMPLE), S(N_SAMPLE), S_HAT(N_SAMPLE)
  REAL(KIND=dp) :: MEAN_M, S_TARGET, FRAC, DR_SAMPLE
  INTEGER(KIND=i4) :: JLO

  N_DENSE = SIZE(R_DENSE)
  ALLOCATE(LOGETA_DENSE(N_DENSE), LOGFH_DENSE(N_DENSE))
  LOGETA_DENSE = LOG(ETA_DENSE)
  LOGFH_DENSE  = LOG(F_HALL_DENSE)

  ! Uniform-in-r resampling of ln(eta), ln(f_H) via bisection + linear
  ! interpolation on the (generally non-uniform, near-duplicate-r-
  ! afflicted) native table.
  DR_SAMPLE = (R_DENSE(N_DENSE) - R_DENSE(1)) / REAL(N_SAMPLE-1, KIND=dp)
  DO J = 1, N_SAMPLE
    R_UNIFORM(J) = R_DENSE(1) + REAL(J-1, KIND=dp)*DR_SAMPLE
  END DO
  R_UNIFORM(N_SAMPLE) = R_DENSE(N_DENSE)   ! exact, not accumulated-roundoff

  DO J = 1, N_SAMPLE
    CALL LOCATE_TABLE(R_DENSE, R_UNIFORM(J), JLO)
    JLO = MAX(1, MIN(JLO, N_DENSE-1))
    FRAC = (R_UNIFORM(J) - R_DENSE(JLO)) / (R_DENSE(JLO+1) - R_DENSE(JLO))
    LOGETA(J) = LOGETA_DENSE(JLO) + FRAC*(LOGETA_DENSE(JLO+1) - LOGETA_DENSE(JLO))
    LOGFH(J)  = LOGFH_DENSE(JLO)  + FRAC*(LOGFH_DENSE(JLO+1)  - LOGFH_DENSE(JLO))
  END DO

  ! |d ln(eta)/dr| + |d ln(f_H)/dr|, centered differences on the now-
  ! uniform sample table, one-sided at the two ends.
  MONITOR(1) = ABS((LOGETA(2)-LOGETA(1))/DR_SAMPLE) + ABS((LOGFH(2)-LOGFH(1))/DR_SAMPLE)
  DO J = 2, N_SAMPLE-1
    MONITOR(J) = ABS((LOGETA(J+1)-LOGETA(J-1))/(2.0_dp*DR_SAMPLE)) &
               + ABS((LOGFH(J+1)-LOGFH(J-1))/(2.0_dp*DR_SAMPLE))
  END DO
  MONITOR(N_SAMPLE) = ABS((LOGETA(N_SAMPLE)-LOGETA(N_SAMPLE-1))/DR_SAMPLE) &
                     + ABS((LOGFH(N_SAMPLE)-LOGFH(N_SAMPLE-1))/DR_SAMPLE)

  MEAN_M = SUM(MONITOR) / REAL(N_SAMPLE, KIND=dp)
  MONITOR = MONITOR + MONITOR_FLOOR_FRAC*MEAN_M

  ! Cumulative trapezoidal integral, normalized to [0,1].
  S(1) = 0.0_dp
  DO J = 2, N_SAMPLE
    S(J) = S(J-1) + 0.5_dp*(MONITOR(J)+MONITOR(J-1))*DR_SAMPLE
  END DO
  S_HAT = S / S(N_SAMPLE)

  ALLOCATE(R_OUT(N_R), DR_OUT(N_R))
  R_OUT(1)   = R_DENSE(1)
  R_OUT(N_R) = R_DENSE(N_DENSE)
  DO K = 2, N_R-1
    S_TARGET = REAL(K-1, KIND=dp) / REAL(N_R-1, KIND=dp)
    ! Locate the bracketing sample-table interval (S_HAT strictly
    ! increasing) and linearly interpolate r.
    DO J = 1, N_SAMPLE-1
      IF (S_HAT(J) <= S_TARGET .AND. S_TARGET <= S_HAT(J+1)) THEN
        FRAC = (S_TARGET - S_HAT(J)) / (S_HAT(J+1) - S_HAT(J))
        R_OUT(K) = R_UNIFORM(J) + FRAC*(R_UNIFORM(J+1)-R_UNIFORM(J))
        EXIT
      END IF
    END DO
  END DO

  DR_OUT(1)   = R_OUT(2) - R_OUT(1)
  DR_OUT(N_R) = R_OUT(N_R) - R_OUT(N_R-1)
  DO K = 2, N_R-1
    DR_OUT(K) = 0.5_dp*(R_OUT(K+1) - R_OUT(K-1))
  END DO
END SUBROUTINE BUILD_EQUIDISTRIBUTED_GRID

END PROGRAM MHDVSH_TOV
