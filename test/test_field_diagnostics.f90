!> Checks FIELD_DIAGNOSTICS' energy/Joule/Poynting formulas implement
!> what's documented (matching analytic_formulas/mhd-vsh-relations.tex
!> Eqs.3,4,6,8) correctly -- NOT a check that those equations are the
!> final word on the physics (that's the user's own derivation; this
!> just checks the code matches it). Two hand-picked modes with known
!> closed-form profiles (Phi'=1 or 0, Phi''=0, Psi'=0 for both, chosen
!> so every derivative is exactly reproduced by the FD operators) give
!> exactly-computable expected values at every grid point; an
!> independently-coded trapezoidal sum gives the expected totals, so
!> nothing here tolerates discretization error against a continuum
!> answer -- only against another discrete evaluation of the same
!> formula.
PROGRAM TEST_FIELD_DIAGNOSTICS
USE KINDS,               ONLY: dp, i4
USE GLOBALS,             ONLY: pi
USE VSH,                 ONLY: YLM_INDEX
USE GRID_RADIAL,         ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,    ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,         ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS,   ONLY: POLOIDAL_ENERGY_DENSITY, TOROIDAL_ENERGY_DENSITY, &
                                TOTAL_POLOIDAL_MAGNETIC_ENERGY, &
                                TOTAL_TOROIDAL_MAGNETIC_ENERGY, TOTAL_MAGNETIC_ENERGY, &
                                POLOIDAL_MAGNETIC_ENERGY_BY_L, TOROIDAL_MAGNETIC_ENERGY_BY_L, &
                                JOULE_DISSIPATION_RATE, POYNTING_FLUX_RATE, &
                                HALL_POYNTING_FLUX_RATE, CURRENT_DENSITY_SQUARED_BY_R, &
                                HALL_COURANT_TIMESTEP
USE VSH,                 ONLY: GWI
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R   = 30
REAL(KIND=dp),    PARAMETER :: R_MIN = 0.5_dp
REAL(KIND=dp),    PARAMETER :: R_MAX = 1.0_dp
INTEGER(KIND=i4), PARAMETER :: LMAX  = 4
REAL(KIND=dp),    PARAMETER :: ETA   = 0.7_dp
REAL(KIND=dp),    PARAMETER :: TOL   = 1.0E-9_dp

! Mode 1: L=2,M=1, Phi(r)=r (Phi'=1, Phi''=0 exactly), Psi=const complex.
INTEGER(KIND=i4), PARAMETER :: L1 = 2, M1 = 1
COMPLEX(KIND=dp), PARAMETER :: PSI1 = CMPLX(1.2_dp, -0.5_dp, KIND=dp)
! Mode 2: L=3,M=-2, Phi=const complex (Phi'=Phi''=0 exactly), Psi=0.
INTEGER(KIND=i4), PARAMETER :: L2 = 3, M2 = -2
COMPLEX(KIND=dp), PARAMETER :: PHI2 = CMPLX(0.8_dp, 0.6_dp, KIND=dp)

! Hall Poynting check: a separate, smaller mode set (own LMAX) so the
! independent triple-mode-sum re-implementation below stays tractable.
INTEGER(KIND=i4), PARAMETER :: LMAX_H = 2
REAL(KIND=dp),    PARAMETER :: F_HALL = 0.4_dp
REAL(KIND=dp),    PARAMETER :: F_HALL_TEST = 0.25_dp

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(SPECTRAL_SCALAR_T) :: PHI, PSI
TYPE(SPECTRAL_SCALAR_T) :: PHI_H, PSI_H
REAL(KIND=dp) :: LAMBDA1, LAMBDA2
REAL(KIND=dp), ALLOCATABLE :: EXP_POL(:), EXP_TOR(:), EXP_JOULE(:)
COMPLEX(KIND=dp) :: CURV_PHI1, CURV_PHI2, FLUX_OUT, FLUX_IN
REAL(KIND=dp) :: EXP_E_POL, EXP_E_TOR, EXP_EDOT_J, EXP_EDOT_S
REAL(KIND=dp) :: GOT_E_POL, GOT_E_TOR, GOT_E_TOTAL, GOT_EDOT_J, GOT_EDOT_S
REAL(KIND=dp) :: GOT_EDOT_HS, EXP_EDOT_HS
REAL(KIND=dp), ALLOCATABLE :: EXP_POL_MODE1(:), EXP_POL_MODE2(:)
REAL(KIND=dp), ALLOCATABLE :: GOT_E_POL_BY_L(:), GOT_E_TOR_BY_L(:)
REAL(KIND=dp) :: EXP_E_POL_MODE1, EXP_E_POL_MODE2
REAL(KIND=dp) :: EXP_TC, GOT_TC, GOT_TC_PROFILE
REAL(KIND=dp) :: J_RMS_EXP(N_R), DR_I, J_CELL
REAL(KIND=dp), ALLOCATABLE :: F_HALL_PROFILE_TEST(:)
INTEGER(KIND=i4) :: IR, N_FAIL

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)
CALL ALLOC_SPECTRAL_SCALAR(PHI, N_R, LMAX)
CALL ALLOC_SPECTRAL_SCALAR(PSI, N_R, LMAX)

DO IR = 1, N_R
  PHI%COEF(IR, YLM_INDEX(L1,M1)) = CMPLX(RGRID%R(IR), 0.0_dp, KIND=dp)
END DO
PSI%COEF(:, YLM_INDEX(L1,M1)) = PSI1
PHI%COEF(:, YLM_INDEX(L2,M2)) = PHI2

LAMBDA1 = REAL(L1*(L1+1), KIND=dp)
LAMBDA2 = REAL(L2*(L2+1), KIND=dp)

! ---- Energy density/totals ----
ALLOCATE(EXP_POL(N_R), EXP_TOR(N_R), EXP_JOULE(N_R))
DO IR = 1, N_R
  EXP_POL(IR) = ( LAMBDA1*(LAMBDA1*RGRID%R(IR)**2/RGRID%R(IR)**2 + 1.0_dp) + &
                  LAMBDA2*(LAMBDA2*ABS(PHI2)**2/RGRID%R(IR)**2) ) / (8.0_dp*pi)
  EXP_TOR(IR) = ( LAMBDA1*ABS(PSI1)**2 ) / (8.0_dp*pi)

  ! Phi''-Lambda*Phi/r**2: mode1 Phi=r (Phi''=0), mode2 Phi=const (Phi''=0)
  CURV_PHI1 = -LAMBDA1*RGRID%R(IR)/RGRID%R(IR)**2
  CURV_PHI2 = -LAMBDA2*PHI2/RGRID%R(IR)**2
  EXP_JOULE(IR) = LAMBDA1*(LAMBDA1*ABS(PSI1)**2/RGRID%R(IR)**2 + ABS(CURV_PHI1)**2) + &
                  LAMBDA2*ABS(CURV_PHI2)**2
END DO

EXP_E_POL  = TRAPZ(EXP_POL, RGRID%R, N_R)
EXP_E_TOR  = TRAPZ(EXP_TOR, RGRID%R, N_R)
EXP_EDOT_J = -(ETA/(4.0_dp*pi)) * TRAPZ(EXP_JOULE, RGRID%R, N_R)

! ---- Per-degree-l poloidal energy: isolate mode1 (L1)/mode2 (L2)'s own
! contribution to EXP_POL above (they're additive, different l, so this
! is just splitting that same sum back into its two summands).
ALLOCATE(EXP_POL_MODE1(N_R), EXP_POL_MODE2(N_R))
DO IR = 1, N_R
  EXP_POL_MODE1(IR) = LAMBDA1*(LAMBDA1*RGRID%R(IR)**2/RGRID%R(IR)**2 + 1.0_dp) / (8.0_dp*pi)
  EXP_POL_MODE2(IR) = LAMBDA2*(LAMBDA2*ABS(PHI2)**2/RGRID%R(IR)**2) / (8.0_dp*pi)
END DO
EXP_E_POL_MODE1 = TRAPZ(EXP_POL_MODE1, RGRID%R, N_R)
EXP_E_POL_MODE2 = TRAPZ(EXP_POL_MODE2, RGRID%R, N_R)

! ---- Poynting: boundary-only, outer minus inner ----
! Psi'=0 for both modes (mode1 const, mode2 zero); Phi'=1 (mode1), 0 (mode2)
CURV_PHI1 = -LAMBDA1*RGRID%R(N_R)/RGRID%R(N_R)**2
CURV_PHI2 = -LAMBDA2*PHI2/RGRID%R(N_R)**2
FLUX_OUT = LAMBDA1*(CMPLX(0.0_dp,0.0_dp,KIND=dp)*CONJG(PSI1) - CONJG(CMPLX(1.0_dp,0.0_dp,KIND=dp))*CURV_PHI1) + &
           LAMBDA2*(CMPLX(0.0_dp,0.0_dp,KIND=dp) - CONJG(CMPLX(0.0_dp,0.0_dp,KIND=dp))*CURV_PHI2)
CURV_PHI1 = -LAMBDA1*RGRID%R(1)/RGRID%R(1)**2
CURV_PHI2 = -LAMBDA2*PHI2/RGRID%R(1)**2
FLUX_IN  = LAMBDA1*(CMPLX(0.0_dp,0.0_dp,KIND=dp)*CONJG(PSI1) - CONJG(CMPLX(1.0_dp,0.0_dp,KIND=dp))*CURV_PHI1) + &
           LAMBDA2*(CMPLX(0.0_dp,0.0_dp,KIND=dp) - CONJG(CMPLX(0.0_dp,0.0_dp,KIND=dp))*CURV_PHI2)
EXP_EDOT_S = -(ETA/(4.0_dp*pi)) * REAL(FLUX_OUT - FLUX_IN, KIND=dp)

! ---- Hall Poynting: own smaller mode set, own field pair ----
! All m=0 (axisymmetric), chosen so GWI(1,0,2,0,1,0) -- a triangle- and
! parity-valid triple (k=1,k'=2,n=1: |1-2|<=1<=3, 1+2+1=4 even) -- picks
! out a genuinely nonzero contribution via TERM2's Psi'_kl*Phi_k'l'*
! conj(Phi'_nm) piece: Psi(1,0) needs a nonzero derivative too (not just
! Phi), hence the linear (not constant) Psi(1,0) profile below.
CALL ALLOC_SPECTRAL_SCALAR(PHI_H, N_R, LMAX_H)
CALL ALLOC_SPECTRAL_SCALAR(PSI_H, N_R, LMAX_H)
DO IR = 1, N_R
  PHI_H%COEF(IR, YLM_INDEX(1,0)) = CMPLX(RGRID%R(IR), 0.0_dp, KIND=dp)
  PSI_H%COEF(IR, YLM_INDEX(1,0)) = CMPLX(0.9_dp*RGRID%R(IR), 0.0_dp, KIND=dp)
END DO
PHI_H%COEF(:, YLM_INDEX(2,0)) = CMPLX(0.8_dp, 0.6_dp, KIND=dp)
PSI_H%COEF(:, YLM_INDEX(2,0)) = CMPLX(0.3_dp, 0.9_dp, KIND=dp)

GOT_EDOT_HS = HALL_POYNTING_FLUX_RATE(PHI_H, PSI_H, OPS, RGRID, F_HALL)
EXP_EDOT_HS = INDEPENDENT_HALL_POYNTING(PHI_H, PSI_H, OPS, RGRID, F_HALL, LMAX_H)

! ---- Hall-CFL Courant timestep: reuses the (PHI,PSI) mode pair and
! EXP_JOULE already derived above for the Joule-rate check -- EXP_JOULE
! IS exactly CURRENT_DENSITY_SQUARED_BY_R's expected per-r output (both
! are the angle-integrated |curl(B)|**2, pre-ETA, pre-radial-integral;
! DPsi_dr's own contribution is correctly absent from EXP_JOULE since
! both modes' Psi is radially constant -- PSI1 literal, mode2 zero).
! EXP_TC below independently re-derives HALL_COURANT_TIMESTEP's own
! cell-averaged-RMS-current formula from EXP_JOULE, not by calling it.
DO IR = 1, N_R
  J_RMS_EXP(IR) = SQRT(EXP_JOULE(IR)/(4.0_dp*pi))
END DO
EXP_TC = HUGE(1.0_dp)
DO IR = 1, N_R-1
  DR_I  = RGRID%R(IR+1) - RGRID%R(IR)
  J_CELL = 0.5_dp*(J_RMS_EXP(IR)+J_RMS_EXP(IR+1))
  EXP_TC = MIN(EXP_TC, DR_I/(F_HALL_TEST*J_CELL))
END DO
GOT_TC = HALL_COURANT_TIMESTEP(PHI, PSI, OPS, RGRID, F_HALL_TEST)

ALLOCATE(F_HALL_PROFILE_TEST(N_R))
F_HALL_PROFILE_TEST = F_HALL_TEST
GOT_TC_PROFILE = HALL_COURANT_TIMESTEP(PHI, PSI, OPS, RGRID, F_HALL_TEST, &
  F_HALL_PROFILE=F_HALL_PROFILE_TEST)

GOT_E_POL   = TOTAL_POLOIDAL_MAGNETIC_ENERGY(PHI, OPS, RGRID)
GOT_E_TOR   = TOTAL_TOROIDAL_MAGNETIC_ENERGY(PSI, RGRID)
GOT_E_TOTAL = TOTAL_MAGNETIC_ENERGY(PHI, PSI, OPS, RGRID)
GOT_EDOT_J  = JOULE_DISSIPATION_RATE(PHI, PSI, OPS, RGRID, ETA)
GOT_EDOT_S  = POYNTING_FLUX_RATE(PHI, PSI, OPS, RGRID, ETA)
ALLOCATE(GOT_E_POL_BY_L(0:LMAX), GOT_E_TOR_BY_L(0:LMAX))
GOT_E_POL_BY_L = POLOIDAL_MAGNETIC_ENERGY_BY_L(PHI, OPS, RGRID)
GOT_E_TOR_BY_L = TOROIDAL_MAGNETIC_ENERGY_BY_L(PSI, RGRID)

N_FAIL = 0
CALL CHECK("poloidal_energy_density", MAXVAL(ABS(POLOIDAL_ENERGY_DENSITY(PHI,OPS,RGRID)-EXP_POL)), N_FAIL)
CALL CHECK("toroidal_energy_density", MAXVAL(ABS(TOROIDAL_ENERGY_DENSITY(PSI,RGRID)-EXP_TOR)), N_FAIL)
CALL CHECK("total_poloidal_energy",   ABS(GOT_E_POL-EXP_E_POL), N_FAIL)
CALL CHECK("total_toroidal_energy",   ABS(GOT_E_TOR-EXP_E_TOR), N_FAIL)
CALL CHECK("total_energy_is_sum",     ABS(GOT_E_TOTAL-(EXP_E_POL+EXP_E_TOR)), N_FAIL)
CALL CHECK("joule_dissipation_rate",  ABS(GOT_EDOT_J-EXP_EDOT_J), N_FAIL)
CALL CHECK("poynting_flux_rate",      ABS(GOT_EDOT_S-EXP_EDOT_S), N_FAIL)
CALL CHECK("hall_poynting_flux_rate", ABS(GOT_EDOT_HS-EXP_EDOT_HS), N_FAIL)
CALL CHECK("poloidal_energy_by_l_mode1", ABS(GOT_E_POL_BY_L(L1)-EXP_E_POL_MODE1), N_FAIL)
CALL CHECK("poloidal_energy_by_l_mode2", ABS(GOT_E_POL_BY_L(L2)-EXP_E_POL_MODE2), N_FAIL)
CALL CHECK("poloidal_energy_by_l_sums_to_total", ABS(SUM(GOT_E_POL_BY_L)-GOT_E_POL), N_FAIL)
CALL CHECK("toroidal_energy_by_l_mode1", ABS(GOT_E_TOR_BY_L(L1)-EXP_E_TOR), N_FAIL)
CALL CHECK("toroidal_energy_by_l_mode2_is_zero", ABS(GOT_E_TOR_BY_L(L2)-0.0_dp), N_FAIL)
CALL CHECK("toroidal_energy_by_l_sums_to_total", ABS(SUM(GOT_E_TOR_BY_L)-GOT_E_TOR), N_FAIL)
! Looser tolerance than every other check here: EXP_JOULE's CURV_PHI1/2
! use the exactly-continuum curvature (-Lambda*Phi/r**2, since both
! modes' Phi''=0 analytically), while CURRENT_DENSITY_SQUARED_BY_R
! computes Phi'' via the FD D2 operator -- a genuine, small (~1e-9,
! largest at the one-sided boundary-row stencils) discretization
! residual, not a bug; JOULE_DISSIPATION_RATE's own check above already
! carries this same residual, just partly smoothed out by its
! trapezoidal radial integral (hence its looser-than-machine-precision
! 1e-11 pass, not 1e-14 like the pure-energy checks).
CALL CHECK_TOL("current_density_squared_by_r", &
  MAXVAL(ABS(CURRENT_DENSITY_SQUARED_BY_R(PHI,PSI,OPS,RGRID)-EXP_JOULE)), N_FAIL, 1.0E-7_dp)
CALL CHECK("hall_courant_timestep", ABS(GOT_TC-EXP_TC), N_FAIL)
CALL CHECK("hall_courant_timestep_profile_matches_scalar", ABS(GOT_TC_PROFILE-GOT_TC), N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_field_diagnostics"
END IF

CONTAINS

SUBROUTINE CHECK(NAME, ERR, N_FAIL)
  CHARACTER(*),     INTENT(IN)    :: NAME
  REAL(KIND=dp),    INTENT(IN)    :: ERR
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  CALL CHECK_TOL(NAME, ERR, N_FAIL, TOL)
END SUBROUTINE CHECK

SUBROUTINE CHECK_TOL(NAME, ERR, N_FAIL, TOL_ARG)
  CHARACTER(*),     INTENT(IN)    :: NAME
  REAL(KIND=dp),    INTENT(IN)    :: ERR
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp),    INTENT(IN)    :: TOL_ARG
  IF (ERR > TOL_ARG) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,ES10.3)') "FAIL  ", NAME, "  err=", ERR
  ELSE
    WRITE(*,'(A,A,A,ES10.3)') "PASS  ", NAME, "  err=", ERR
  END IF
END SUBROUTINE CHECK_TOL

FUNCTION TRAPZ(F, R, N) RESULT(TOTAL)
  INTEGER(KIND=i4), INTENT(IN) :: N
  REAL(KIND=dp),    INTENT(IN) :: F(N), R(N)
  REAL(KIND=dp) :: TOTAL
  INTEGER(KIND=i4) :: I
  TOTAL = 0.0_dp
  DO I = 1, N-1
    TOTAL = TOTAL + 0.5_dp*(F(I)+F(I+1))*(R(I+1)-R(I))
  END DO
END FUNCTION TRAPZ

!> Independent re-implementation of HALL_POYNTING_FLUX_RATE's formula
!> (analytic_formulas/mhd-vsh-relations.tex, "Energy Budget in Hall
!> Limit"), deliberately structured differently from the production code
!> (two separate full boundary passes via BOUNDARY_PASS, values/
!> derivatives fetched fresh per mode via FIELD_AT/CURVATURE_AT rather
!> than precomputed into shared arrays) so this isn't just a copy of the
!> same code -- reduces, though can never eliminate, the chance of the
!> same transcription error appearing in both places.
FUNCTION INDEPENDENT_HALL_POYNTING(PHIFIELD, PSIFIELD, OPSARG, RGRIDARG, F_HALL_ARG, LMAXH) RESULT(E_DOT_S)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHIFIELD, PSIFIELD
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPSARG
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRIDARG
  REAL(KIND=dp),           INTENT(IN) :: F_HALL_ARG
  INTEGER(KIND=i4),        INTENT(IN) :: LMAXH
  REAL(KIND=dp) :: E_DOT_S
  COMPLEX(KIND=dp) :: SUM_OUT, SUM_IN
  SUM_OUT = BOUNDARY_PASS(PHIFIELD, PSIFIELD, OPSARG, RGRIDARG, LMAXH, RGRIDARG%N)
  SUM_IN  = BOUNDARY_PASS(PHIFIELD, PSIFIELD, OPSARG, RGRIDARG, LMAXH, 1_i4)
  E_DOT_S = -(F_HALL_ARG/(8.0_dp*pi)) * REAL(SUM_OUT - SUM_IN, KIND=dp)
END FUNCTION INDEPENDENT_HALL_POYNTING

!> One full triple-mode-sum pass of the Hall Poynting bracket, evaluated
!> at radial row IROW.
FUNCTION BOUNDARY_PASS(PHIFIELD, PSIFIELD, OPSARG, RGRIDARG, LMAXH, IROW) RESULT(TOTAL)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHIFIELD, PSIFIELD
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPSARG
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRIDARG
  INTEGER(KIND=i4),        INTENT(IN) :: LMAXH, IROW
  COMPLEX(KIND=dp) :: TOTAL
  COMPLEX(KIND=dp) :: PSI_KL_V, PSI_KL_D, PHI_KL_CURV
  COMPLEX(KIND=dp) :: PHI_K2L2_V, PHI_K2L2_D, PSI_K2L2_V, PSI_K2L2_D
  COMPLEX(KIND=dp) :: PHI_NM_V, PHI_NM_D, PSI_NM_V, PSI_NM_D
  COMPLEX(KIND=dp) :: COUPLING
  REAL(KIND=dp)    :: LK, LK2, LN
  INTEGER(KIND=i4) :: K, L, K2, L2, N, M

  TOTAL = (0.0_dp, 0.0_dp)
  DO K = 0, LMAXH
    LK = REAL(K*(K+1), KIND=dp)
    DO L = -K, K
      CALL FIELD_AT(PSIFIELD, OPSARG, K, L, IROW, PSI_KL_V, PSI_KL_D)
      CALL CURVATURE_AT(PHIFIELD, OPSARG, RGRIDARG, K, L, IROW, PHI_KL_CURV)
      DO K2 = 0, LMAXH
        LK2 = REAL(K2*(K2+1), KIND=dp)
        DO L2 = -K2, K2
          CALL FIELD_AT(PHIFIELD, OPSARG, K2, L2, IROW, PHI_K2L2_V, PHI_K2L2_D)
          CALL FIELD_AT(PSIFIELD, OPSARG, K2, L2, IROW, PSI_K2L2_V, PSI_K2L2_D)
          DO N = 0, LMAXH
            LN = REAL(N*(N+1), KIND=dp)
            DO M = -N, N
              COUPLING = GWI(K, L, K2, L2, N, M)
              IF (COUPLING == (0.0_dp,0.0_dp)) CYCLE
              CALL FIELD_AT(PHIFIELD, OPSARG, N, M, IROW, PHI_NM_V, PHI_NM_D)
              CALL FIELD_AT(PSIFIELD, OPSARG, N, M, IROW, PSI_NM_V, PSI_NM_D)

              TOTAL = TOTAL + COUPLING * ( &
                -LK*(LK2+LN-LK) * PSI_KL_V * &
                  (PSI_K2L2_V*CONJG(PSI_NM_V) + PHI_K2L2_D*CONJG(PHI_NM_D)) &
                + LK2*(LK+LN-LK2) * PHI_K2L2_V * &
                  (PSI_KL_D*CONJG(PHI_NM_D) - PSI_NM_V*PHI_KL_CURV) )
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
END FUNCTION BOUNDARY_PASS

!> Value and radial derivative of FIELD's (K,L) mode at radial row IROW,
!> via direct FD stencil dot products.
SUBROUTINE FIELD_AT(FIELDARG, OPSARG, K, L, IROW, VAL, DERIV)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN)  :: FIELDARG
  TYPE(RADIAL_OPERATOR_T), INTENT(IN)  :: OPSARG
  INTEGER(KIND=i4),        INTENT(IN)  :: K, L, IROW
  COMPLEX(KIND=dp),        INTENT(OUT) :: VAL, DERIV
  INTEGER(KIND=i4) :: IDX
  IDX = YLM_INDEX(K, L)
  VAL   = FIELDARG%COEF(IROW, IDX)
  DERIV = SUM(OPSARG%D1(IROW,:)*FIELDARG%COEF(:,IDX))
END SUBROUTINE FIELD_AT

!> Phi''-l(l+1)/r**2*Phi for PHIFIELD's (K,L) mode at radial row IROW.
SUBROUTINE CURVATURE_AT(PHIFIELD, OPSARG, RGRIDARG, K, L, IROW, CURV)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN)  :: PHIFIELD
  TYPE(RADIAL_OPERATOR_T), INTENT(IN)  :: OPSARG
  TYPE(RADIAL_GRID_T),     INTENT(IN)  :: RGRIDARG
  INTEGER(KIND=i4),        INTENT(IN)  :: K, L, IROW
  COMPLEX(KIND=dp),        INTENT(OUT) :: CURV
  COMPLEX(KIND=dp) :: D2VAL
  INTEGER(KIND=i4) :: IDX
  IDX = YLM_INDEX(K, L)
  D2VAL = SUM(OPSARG%D2(IROW,:)*PHIFIELD%COEF(:,IDX))
  CURV = D2VAL - REAL(K*(K+1),KIND=dp)*PHIFIELD%COEF(IROW,IDX)/RGRIDARG%R(IROW)**2
END SUBROUTINE CURVATURE_AT

END PROGRAM TEST_FIELD_DIAGNOSTICS
