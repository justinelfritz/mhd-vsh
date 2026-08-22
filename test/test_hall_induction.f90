!> Checks HALL_INDUCTION_RHS against a hand-derived closed form, not
!> against itself -- an independently-coded reference computation, same
!> methodology as test_field_diagnostics.f90's INDEPENDENT_HALL_POYNTING
!> check on HALL_POYNTING_FLUX_RATE.
!>
!> Setup: exactly two populated modes, A=(K_A,0) and B=(K_B,1), both with
!> CONSTANT (r-independent) Phi/Psi values -- chosen deliberately, not
!> just for simplicity. Constant profiles make Phi'=Psi'=0 EXACTLY for
!> both modes (any consistent FD stencil reproduces a constant's
!> derivative as exactly zero), which kills two of the formula's four
!> term-groups outright (Phi_dot's I-term, and Psi_dot's term2 both
!> multiply by a first derivative that's now identically zero) --
!> leaving only Phi_dot's J-term and Psi_dot's term1+term4 to verify,
!> all hand-derivable in closed form (term1 still needs a genuine radial
!> derivative of a non-constant r-dependent bracket, computed here via
!> the SAME OPS%D1 operator HALL_INDUCTION_RHS uses -- matching, not
!> approximating, its discrete derivative, per this codebase's existing
!> "checked against another discrete evaluation of the same formula, not
!> continuum truth" test philosophy).
!>
!> Target order M=1=0+1 (forced by the m=l+l' selection rule); choosing
!> M different from 2*0=0 and 2*1=2 means the "self-coupling" terms
!> (A,A) and (B,B) are forced to exactly zero by that same rule, so only
!> the (k,l,k',l')=(A,B) and (B,A) cross-assignments can contribute --
!> both are included in the reference sum below.
PROGRAM TEST_HALL_INDUCTION
USE KINDS,             ONLY: dp, i4
USE VSH,               ONLY: YLM_INDEX, GWI, GWJ
USE GRID_RADIAL,       ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,  ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,       ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE HALL_INDUCTION,    ONLY: HALL_INDUCTION_RHS
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R = 20
REAL(KIND=dp),    PARAMETER :: R_MIN = 0.5_dp
REAL(KIND=dp),    PARAMETER :: R_MAX = 1.0_dp
INTEGER(KIND=i4), PARAMETER :: LMAX_SEARCH = 6
REAL(KIND=dp),    PARAMETER :: F_HALL = 0.7_dp
REAL(KIND=dp),    PARAMETER :: TOL = 1.0E-9_dp

INTEGER(KIND=i4), PARAMETER :: L_A = 0, L_B = 1, M_TARGET = 1
COMPLEX(KIND=dp), PARAMETER :: PHI_A = CMPLX(1.3_dp, -0.4_dp, KIND=dp)
COMPLEX(KIND=dp), PARAMETER :: PSI_A = CMPLX(0.6_dp,  0.2_dp, KIND=dp)
COMPLEX(KIND=dp), PARAMETER :: PHI_B = CMPLX(-0.5_dp, 0.9_dp, KIND=dp)
COMPLEX(KIND=dp), PARAMETER :: PSI_B = CMPLX(0.7_dp, -0.3_dp, KIND=dp)

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(SPECTRAL_SCALAR_T) :: PHI, PSI, PHI_DOT, PSI_DOT
INTEGER(KIND=i4) :: K_A, K_B, N_TARGET
INTEGER(KIND=i4) :: N_FAIL
COMPLEX(KIND=dp), ALLOCATABLE :: REF_PHI_DOT(:), REF_PSI_DOT(:)

CALL FIND_NONTRIVIAL_TRIPLE(K_A, K_B, N_TARGET)
WRITE(*,'(A,I0,A,I0,A,I0,A)') 'Using (K_A,K_B,N) = (', K_A, ',', K_B, ',', N_TARGET, ')'

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)
CALL ALLOC_SPECTRAL_SCALAR(PHI, N_R, LMAX_SEARCH)
CALL ALLOC_SPECTRAL_SCALAR(PSI, N_R, LMAX_SEARCH)
PHI%COEF(:, YLM_INDEX(K_A,L_A)) = PHI_A
PSI%COEF(:, YLM_INDEX(K_A,L_A)) = PSI_A
PHI%COEF(:, YLM_INDEX(K_B,L_B)) = PHI_B
PSI%COEF(:, YLM_INDEX(K_B,L_B)) = PSI_B

CALL HALL_INDUCTION_RHS(PHI, PSI, OPS, RGRID, F_HALL, PHI_DOT, PSI_DOT)

ALLOCATE(REF_PHI_DOT(N_R), REF_PSI_DOT(N_R))
CALL REFERENCE_RHS(K_A, K_B, N_TARGET, REF_PHI_DOT, REF_PSI_DOT)

N_FAIL = 0
CALL CHECK("phi_dot_nm", MAXVAL(ABS(PHI_DOT%COEF(:,YLM_INDEX(N_TARGET,M_TARGET)) - REF_PHI_DOT)), N_FAIL)
CALL CHECK("psi_dot_nm", MAXVAL(ABS(PSI_DOT%COEF(:,YLM_INDEX(N_TARGET,M_TARGET)) - REF_PSI_DOT)), N_FAIL)
CALL CHECK("phi_dot_nm_nontrivial", -MAXVAL(ABS(REF_PHI_DOT)), N_FAIL, MUST_BE_NEGATIVE=.TRUE.)
CALL CHECK("psi_dot_nm_nontrivial", -MAXVAL(ABS(REF_PSI_DOT)), N_FAIL, MUST_BE_NEGATIVE=.TRUE.)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_hall_induction"
END IF

CONTAINS

!> Finds a (K_A,K_B,N) triple (K_A,K_B>=1, N in 0..LMAX_SEARCH) for which
!> at least one of the couplings feeding the surviving terms (Phi_dot's
!> J-term, Psi_dot's term1/term4) is nonzero for both cross-assignments
!> (A,B) and (B,A), so the test isn't a trivial 0=0 comparison.
SUBROUTINE FIND_NONTRIVIAL_TRIPLE(KA, KB, NT)
  INTEGER(KIND=i4), INTENT(OUT) :: KA, KB, NT
  INTEGER(KIND=i4) :: TRY_KA, TRY_KB, TRY_N
  LOGICAL :: FOUND
  FOUND = .FALSE.
  DO TRY_KA = 1, LMAX_SEARCH
    DO TRY_KB = 1, LMAX_SEARCH
      DO TRY_N = 0, LMAX_SEARCH
        IF (ABS(GWJ(TRY_KA,L_A,TRY_KB,L_B,TRY_N,M_TARGET)) > 1.0E-10_dp .OR. &
            ABS(GWJ(TRY_KB,L_B,TRY_KA,L_A,TRY_N,M_TARGET)) > 1.0E-10_dp .OR. &
            ABS(GWI(TRY_KB,L_B,TRY_KA,L_A,TRY_N,M_TARGET)) > 1.0E-10_dp .OR. &
            ABS(GWI(TRY_KA,L_A,TRY_KB,L_B,TRY_N,M_TARGET)) > 1.0E-10_dp) THEN
          KA = TRY_KA; KB = TRY_KB; NT = TRY_N
          FOUND = .TRUE.
          RETURN
        END IF
      END DO
    END DO
  END DO
  IF (.NOT. FOUND) THEN
    WRITE(*,'(A)') 'No nontrivial (K_A,K_B,N) triple found -- widen LMAX_SEARCH.'
    STOP 1
  END IF
END SUBROUTINE FIND_NONTRIVIAL_TRIPLE

!> Independently-derived reference RHS at mode (N_TARGET,M_TARGET), for
!> the two-constant-mode setup above. Sums the two surviving
!> cross-assignments (k,l,k',l')=(A,B) and (B,A); the self-assignments
!> (A,A),(B,B) are omitted entirely since the m=l+l' selection rule
!> forces them to zero for this M_TARGET by construction (not computed
!> and dropped -- never nonzero in the first place).
SUBROUTINE REFERENCE_RHS(KA, KB, NT, PHI_DOT_OUT, PSI_DOT_OUT)
  INTEGER(KIND=i4), INTENT(IN)  :: KA, KB, NT
  COMPLEX(KIND=dp), INTENT(OUT) :: PHI_DOT_OUT(N_R), PSI_DOT_OUT(N_R)
  PHI_DOT_OUT = (0.0_dp,0.0_dp)
  PSI_DOT_OUT = (0.0_dp,0.0_dp)
  CALL ADD_ASSIGNMENT(KA, L_A, PHI_A, PSI_A, KB, L_B, PHI_B, PSI_B, NT, PHI_DOT_OUT, PSI_DOT_OUT)
  CALL ADD_ASSIGNMENT(KB, L_B, PHI_B, PSI_B, KA, L_A, PHI_A, PSI_A, NT, PHI_DOT_OUT, PSI_DOT_OUT)
END SUBROUTINE REFERENCE_RHS

!> Adds one (k,l,k',l')=(KX,LX)->(KY,LY) assignment's contribution to the
!> target mode (NT,M_TARGET) -- Phi_dot's J-term (only surviving Phi_dot
!> term for constant profiles) plus Psi_dot's term1+term4 (only surviving
!> Psi_dot terms), evaluated at every radial row.
SUBROUTINE ADD_ASSIGNMENT(KX, LX, PHIX, PSIX, KY, LY, PHIY, PSIY, NT, PHI_DOT_OUT, PSI_DOT_OUT)
  INTEGER(KIND=i4), INTENT(IN)    :: KX, LX, KY, LY, NT
  COMPLEX(KIND=dp), INTENT(IN)    :: PHIX, PSIX, PHIY, PSIY
  COMPLEX(KIND=dp), INTENT(INOUT) :: PHI_DOT_OUT(N_R), PSI_DOT_OUT(N_R)
  REAL(KIND=dp)    :: LX_, LY_, LN_, R2(N_R)
  COMPLEX(KIND=dp) :: CURV_PHIX(N_R), CURV_PHIY(N_R), BRACKET1(N_R), DBRACKET1(N_R)
  COMPLEX(KIND=dp) :: COUPLING_J_STD, COUPLING_J_SWP, COUPLING_I_SWP
  INTEGER(KIND=i4) :: IR

  LX_ = REAL(KX*(KX+1), KIND=dp)
  LY_ = REAL(KY*(KY+1), KIND=dp)
  LN_ = REAL(NT*(NT+1), KIND=dp)
  R2 = RGRID%R**2
  CURV_PHIX = -LX_*PHIX/R2   ! Phi_X''=0 exactly (constant), curvature = -Lambda_X/r**2 * Phi_X
  CURV_PHIY = -LY_*PHIY/R2

  COUPLING_J_STD = GWJ(KX,LX,KY,LY,NT,M_TARGET)   ! J^{nm}_{klk'l'}, (k,l)=(X), (k',l')=(Y)
  COUPLING_J_SWP = GWJ(KY,LY,KX,LX,NT,M_TARGET)   ! J^{nm}_{k'l'kl}
  COUPLING_I_SWP = GWI(KY,LY,KX,LX,NT,M_TARGET)   ! I^{nm}_{k'l'kl}

  ! Phi_dot J-term: -(f_H/Lambda_n)*(Lambda_Y/r**2)*(Psi_X*Psi_Y + Phi_Y*Phi^(1)_X) * J^{nm}_{klk'l'}
  PHI_DOT_OUT = PHI_DOT_OUT - (F_HALL/LN_)*(LY_/R2)*(PSIX*PSIY + PHIY*CURV_PHIX) * COUPLING_J_STD

  ! Psi_dot term1: d/dr[(f_H/r**2)*(Phi^(1)_Y*Phi_X + Psi_X*Psi_Y)] * ((Ln+Ly-Lx)/(2Ln))*Lx * I^{nm}_{k'l'kl}
  BRACKET1 = (F_HALL/R2) * (CURV_PHIY*PHIX + PSIX*PSIY)
  DBRACKET1 = MATMUL(OPS%D1, BRACKET1)
  PSI_DOT_OUT = PSI_DOT_OUT + DBRACKET1 * ((LN_+LY_-LX_)/(2.0_dp*LN_)) * LX_ * COUPLING_I_SWP

  ! Psi_dot term4: (f_H/r**2)*(Psi_X*Phi^(1)_Y - Psi'_Y*Phi'_X) * J^{nm}_{k'l'kl}  [Psi'_Y=Phi'_X=0]
  PSI_DOT_OUT = PSI_DOT_OUT + (F_HALL/R2)*(PSIX*CURV_PHIY) * COUPLING_J_SWP

  ! Silence unused-arg warning; IR intentionally unused (array ops throughout)
  IR = 0
END SUBROUTINE ADD_ASSIGNMENT

SUBROUTINE CHECK(NAME, ERR, N_FAIL, MUST_BE_NEGATIVE)
  CHARACTER(LEN=*),  INTENT(IN)    :: NAME
  REAL(KIND=dp),     INTENT(IN)    :: ERR
  INTEGER(KIND=i4),  INTENT(INOUT) :: N_FAIL
  LOGICAL, OPTIONAL, INTENT(IN)    :: MUST_BE_NEGATIVE
  LOGICAL :: FAILED
  IF (PRESENT(MUST_BE_NEGATIVE)) THEN
    FAILED = ERR >= 0.0_dp   ! ERR = -MAXVAL(ABS(x)); must be < 0, i.e. x nonzero somewhere
  ELSE
    FAILED = ERR > TOL
  END IF
  IF (FAILED) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,ES10.3)') "FAIL  ", NAME, "  err=", ERR
  ELSE
    WRITE(*,'(A,A,A,ES10.3)') "PASS  ", NAME, "  err=", ERR
  END IF
END SUBROUTINE CHECK

END PROGRAM TEST_HALL_INDUCTION
