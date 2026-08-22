MODULE HALL_INDUCTION
!> The weak-Hall induction equations' RHS (Phi-dot, Psi-dot), transcribed
!> directly from analytic_formulas/mhd-vsh-relations.tex, "The Hall
!> induction equations are:" (lines 392-404). A pure, side-effect-free
!> evaluator of state -- no boundary condition is applied to its output
!> (that's the caller's job, same convention as RADIAL_OPERATORS
!> excluding BC/curvature "applied later, at system-assembly time"), and
!> it owns no state type and no time-integration decision, matching
!> FIELD_DIAGNOSTICS' layering (a regime-agnostic physics primitive, not
!> regime plumbing) -- see HALL_POYNTING_FLUX_RATE there, the direct
!> precedent this module's loop structure mirrors.
!>
!> F_HALL is a plain REAL(dp) scalar, spatially uniform, unless the
!> optional F_HALL_PROFILE argument is present (one value per radial
!> row, e.g. from CRUST_CONDUCTIVITY::ETA_AND_F_HALL_AT's self-consistent
!> f_H(r)=c/(4*pi*e*n_e(r))) -- absent, behavior is bit-identical to the
!> uniform case (test_hall_induction.f90 confirms this).
!>
!> @warning Index-order care: the four terms below use TWO different
!>   orderings of the (k,l),(k',l') pair through the I/J coupling
!>   coefficients -- I^{nm}_{klk'l'} vs I^{nm}_{k'l'kl}, J^{nm}_{klk'l'}
!>   vs J^{nm}_{k'l'kl} -- NOT a typo in the source document. GWI is
!>   exactly SYMMETRIC under that swap and GWJ exactly ANTISYMMETRIC
!>   (confirmed numerically to ~1e-14 across 1392 nonzero-coupling
!>   sextuples, not just asserted -- Justin Elfritz, 2026-08-20), so only
!>   the "standard"-order GWI/GWJ is actually called below; the swapped
!>   value is obtained algebraically (+COUPLING_I_STD, -COUPLING_J_STD)
!>   rather than a second GWI/GWJ call, halving the coupling-coefficient
!>   cost. Terms are still written out with their tex-matching swapped
!>   subscript in comments for direct correspondence to the source
!>   document.
!>
!> @warning Cost: two independent restrictions keep this tractable at
!>   large LMAX. (1) The outer (k,l),(k',l') pair only ranges over modes
!>   actually nonzero in the current PHI/PSI (an exact optimization, not
!>   an approximation -- any pair with a zero mode contributes exactly
!>   zero regardless of coupling, so skipping it changes nothing). This
!>   is what makes a sparse/growing-from-a-single-mode state (e.g. an
!>   axisymmetric cascade seeded from one mode) cheap even at LMAX=30,
!>   where the full (LMAX+1)**2=961-mode space would otherwise force an
!>   infeasible ~961**3 raw loop. (2) The target m is forced to l+l' by
!>   the coupling coefficients' own selection rule (GWI/GWJ vanish unless
!>   m=l+l'), so for each active (k,l),(k',l') pair only ONE m is ever
!>   evaluated (computed directly, not searched), and n is restricted to
!>   the triangle band |k-k'|<=n<=k+k' (also exact, not approximate --
!>   GWI/GWJ are zero outside it). What's NOT yet done: even within an
!>   active pair's own bracket precompute, the per-radial-node cost is
!>   still full (N_R each); for a genuinely dense/high-order state (most
!>   of the (LMAX+1)**2 modes populated) this reduces to the previously
!>   documented O(Nlm**3 * N_R) worst case -- deferred until a workload
!>   actually needs it, per the project's own Hall-regime design plan
!>   (see ROADMAP.md).
USE KINDS,            ONLY: dp, i4
USE GRID_RADIAL,      ONLY: RADIAL_GRID_T
USE RADIAL_OPERATORS, ONLY: RADIAL_OPERATOR_T, ADD_CURVATURE_TERM
USE FIELD_TYPES,      ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE VSH,              ONLY: YLM_INDEX, GWI, GWJ
IMPLICIT NONE
PRIVATE
PUBLIC :: HALL_INDUCTION_RHS

CONTAINS

!> @param PHI, PSI Current state (poloidal/toroidal potentials).
!> @param OPS Radial derivative operators (D1, D2) PHI/PSI were built on.
!> @param RGRID Radial grid PHI/OPS were built on.
!> @param F_HALL Hall prefactor c/(4*pi*e*n_e), spatially uniform unless
!>   F_HALL_PROFILE is present (see module header).
!> @param F_HALL_PROFILE Optional per-radial-row Hall prefactor (size
!>   RGRID%N), overriding F_HALL row-by-row when present.
!> Returns: PHI_DOT, PSI_DOT (freshly allocated SPECTRAL_SCALAR_T, same
!>   shape as PHI/PSI) -- the Hall induction equations' RHS, interior
!>   values only, no boundary condition applied.
SUBROUTINE HALL_INDUCTION_RHS(PHI, PSI, OPS, RGRID, F_HALL, PHI_DOT, PSI_DOT, F_HALL_PROFILE)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN)  :: PHI, PSI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN)  :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN)  :: RGRID
  REAL(KIND=dp),           INTENT(IN)  :: F_HALL
  TYPE(SPECTRAL_SCALAR_T), INTENT(OUT) :: PHI_DOT, PSI_DOT
  REAL(KIND=dp), OPTIONAL, INTENT(IN)  :: F_HALL_PROFILE(:)

  COMPLEX(KIND=dp), ALLOCATABLE :: PHI_ALL(:,:), PSI_ALL(:,:)
  COMPLEX(KIND=dp), ALLOCATABLE :: DPHI_ALL(:,:), DPSI_ALL(:,:), CURV_PHI_ALL(:,:)
  COMPLEX(KIND=dp), ALLOCATABLE :: BRACKET1(:), BRACKET3(:), DBRACKET1(:), DBRACKET3(:)
  REAL(KIND=dp),    ALLOCATABLE :: D2_CURV(:,:)
  REAL(KIND=dp),    ALLOCATABLE :: R2(:), F_HALL_ARR(:)
  COMPLEX(KIND=dp) :: COUPLING_I_STD, COUPLING_I_SWP, COUPLING_J_STD, COUPLING_J_SWP
  REAL(KIND=dp)    :: LAMBDA_K, LAMBDA_K2, LAMBDA_N
  INTEGER(KIND=i4) :: N_R, NLM
  INTEGER(KIND=i4) :: K, L, K2, L2, N, M, IDX_KL, IDX_K2L2, IDX_NM
  INTEGER(KIND=i4), ALLOCATABLE :: ACTIVE_K(:), ACTIVE_L(:), ACTIVE_IDX(:)
  INTEGER(KIND=i4) :: N_ACTIVE, IA, IA2, KK, LL, IDX, N_LO, N_HI
  REAL(KIND=dp), PARAMETER :: ACTIVE_TOL = 1.0E-13_dp

  N_R = RGRID%N
  NLM = PHI%NLM
  CALL ALLOC_SPECTRAL_SCALAR(PHI_DOT, N_R, PHI%LMAX)
  CALL ALLOC_SPECTRAL_SCALAR(PSI_DOT, N_R, PHI%LMAX)
  PHI_DOT%COEF = (0.0_dp, 0.0_dp)
  PSI_DOT%COEF = (0.0_dp, 0.0_dp)

  ALLOCATE(PHI_ALL(N_R,NLM), PSI_ALL(N_R,NLM))
  ALLOCATE(DPHI_ALL(N_R,NLM), DPSI_ALL(N_R,NLM), CURV_PHI_ALL(N_R,NLM))
  ALLOCATE(BRACKET1(N_R), BRACKET3(N_R), DBRACKET1(N_R), DBRACKET3(N_R))
  ALLOCATE(D2_CURV(N_R,N_R), R2(N_R), F_HALL_ARR(N_R))

  IF (PRESENT(F_HALL_PROFILE)) THEN
    F_HALL_ARR = F_HALL_PROFILE
  ELSE
    F_HALL_ARR = F_HALL
  END IF
  R2 = RGRID%R**2
  PHI_ALL = PHI%COEF
  PSI_ALL = PSI%COEF
  DPHI_ALL = MATMUL(OPS%D1, PHI%COEF)
  DPSI_ALL = MATMUL(OPS%D1, PSI%COEF)
  DO K = 0, PHI%LMAX
    D2_CURV = OPS%D2
    CALL ADD_CURVATURE_TERM(D2_CURV, K, RGRID)
    DO L = -K, K
      IDX_KL = YLM_INDEX(K, L)
      CURV_PHI_ALL(:,IDX_KL) = MATMUL(D2_CURV, PHI%COEF(:,IDX_KL))
    END DO
  END DO

  ! Active-mode list: any (k,l) with a nonzero Phi or Psi coefficient
  ! anywhere in r. Exact restriction, not an approximation -- see module
  ! header. Bounds the (k,l),(k',l') pair enumeration by N_ACTIVE**2
  ! instead of NLM**2.
  ALLOCATE(ACTIVE_K(NLM), ACTIVE_L(NLM), ACTIVE_IDX(NLM))
  N_ACTIVE = 0
  DO KK = 0, PHI%LMAX
    DO LL = -KK, KK
      IDX = YLM_INDEX(KK, LL)
      IF (ANY(ABS(PHI_ALL(:,IDX)) > ACTIVE_TOL) .OR. ANY(ABS(PSI_ALL(:,IDX)) > ACTIVE_TOL)) THEN
        N_ACTIVE = N_ACTIVE + 1
        ACTIVE_K(N_ACTIVE) = KK
        ACTIVE_L(N_ACTIVE) = LL
        ACTIVE_IDX(N_ACTIVE) = IDX
      END IF
    END DO
  END DO

  DO IA = 1, N_ACTIVE
    K = ACTIVE_K(IA); L = ACTIVE_L(IA); IDX_KL = ACTIVE_IDX(IA)
    LAMBDA_K = REAL(K*(K+1), KIND=dp)
    DO IA2 = 1, N_ACTIVE
      K2 = ACTIVE_K(IA2); L2 = ACTIVE_L(IA2); IDX_K2L2 = ACTIVE_IDX(IA2)
      LAMBDA_K2 = REAL(K2*(K2+1), KIND=dp)

          ! Psi_dot's term1/term3 brackets depend only on (k,l,k',l'),
          ! not on the target (n,m) -- precomputed once per pair and
          ! reused for every (n,m) below.
          BRACKET1 = (F_HALL_ARR/R2) * ( CURV_PHI_ALL(:,IDX_K2L2)*PHI_ALL(:,IDX_KL) + &
                                      PSI_ALL(:,IDX_KL)*PSI_ALL(:,IDX_K2L2) )
          DBRACKET1 = MATMUL(OPS%D1, BRACKET1)
          BRACKET3 = (F_HALL_ARR/R2) * ( PSI_ALL(:,IDX_K2L2)*DPHI_ALL(:,IDX_KL) - &
                                      DPSI_ALL(:,IDX_KL)*PHI_ALL(:,IDX_K2L2) )
          DBRACKET3 = MATMUL(OPS%D1, BRACKET3)

          ! m is forced to l+l' by GWI/GWJ's own selection rule -- only
          ! that one m is ever evaluated, not searched. n is restricted
          ! to the triangle band |k-k'|<=n<=k+k' (also exact: GWI/GWJ
          ! vanish outside it), intersected with |m|<=n and n<=LMAX.
          ! n=0 excluded unconditionally: both equations carry an
          ! explicit 1/Lambda_n factor (f_H/Lambda_n in Phi_dot;
          ! (Ln+Lk'-Lk)/(2Ln) and Lk'/Ln in Psi_dot), singular at
          ! Lambda_0=0 -- physically consistent, since a degree-0
          ! poloidal/toroidal potential carries no field anyway (same
          ! reason every seed/test IC in this codebase starts at l>=1).
          M = L + L2
          IF (ABS(M) > PHI%LMAX) CYCLE
          N_LO = MAX(1_i4, ABS(K-K2), ABS(M))
          N_HI = MIN(PHI%LMAX, K+K2)
          DO N = N_LO, N_HI
            LAMBDA_N = REAL(N*(N+1), KIND=dp)
              COUPLING_I_STD = GWI(K, L, K2, L2, N, M)   ! I^{nm}_{klk'l'}
              COUPLING_J_STD = GWJ(K, L, K2, L2, N, M)   ! J^{nm}_{klk'l'}
              IF (ABS(COUPLING_I_STD) < 1.0E-13_dp .AND. ABS(COUPLING_J_STD) < 1.0E-13_dp) CYCLE
              ! I^{nm}_{k'l'kl}=I^{nm}_{klk'l'} (GWI symmetric under the
              ! (k,l)<->(k',l') swap) and J^{nm}_{k'l'kl}=-J^{nm}_{klk'l'}
              ! (GWJ antisymmetric under it) -- confirmed numerically to
              ! ~1e-14 across 1392 nonzero sextuples (k,l,k',l' up to
              ! degree 4), not just asserted; halves the GWI/GWJ call
              ! count versus evaluating both orderings independently.
              COUPLING_I_SWP =  COUPLING_I_STD
              COUPLING_J_SWP = -COUPLING_J_STD
              IDX_NM = YLM_INDEX(N, M)

              ! Phi_dot_nm: (f_H/Lambda_n) * [ I-term - J-term ]
              !   I-term: (Lambda_k'/r**2)*((Ln+Lk-Lk')/2)*(Psi_k'l'*Phi'_kl - Phi_k'l'*Psi'_kl) * I^{nm}_{klk'l'}
              !   J-term: (Lambda_k'/r**2)*(Psi_kl*Psi_k'l' + Phi_k'l'*Phi^(1)_kl) * J^{nm}_{klk'l'}
              PHI_DOT%COEF(:,IDX_NM) = PHI_DOT%COEF(:,IDX_NM) + (F_HALL_ARR/LAMBDA_N) * ( &
                (LAMBDA_K2/R2)*((LAMBDA_N+LAMBDA_K-LAMBDA_K2)/2.0_dp) * &
                  (PSI_ALL(:,IDX_K2L2)*DPHI_ALL(:,IDX_KL) - PHI_ALL(:,IDX_K2L2)*DPSI_ALL(:,IDX_KL)) * &
                  COUPLING_I_STD &
                - (LAMBDA_K2/R2)*(PSI_ALL(:,IDX_KL)*PSI_ALL(:,IDX_K2L2) + &
                  PHI_ALL(:,IDX_K2L2)*CURV_PHI_ALL(:,IDX_KL)) * COUPLING_J_STD )

              ! Psi_dot_nm: term1 + term2 - term3 + term4
              !   term1: d/dr[bracket1] * ((Ln+Lk'-Lk)/(2Ln)) * Lk * I^{nm}_{k'l'kl}
              !   term2: (f_H/r**2)*(Psi'_k'l'*Psi_kl + Phi'_kl*Phi^(1)_k'l') * ((Lk+Lk'-Ln)/2) * I^{nm}_{klk'l'}
              !   term3: d/dr[bracket3] * (Lk'/Ln) * J^{nm}_{klk'l'}
              !   term4: (f_H/r**2)*(Psi_kl*Phi^(1)_k'l' - Psi'_k'l'*Phi'_kl) * J^{nm}_{k'l'kl}
              PSI_DOT%COEF(:,IDX_NM) = PSI_DOT%COEF(:,IDX_NM) &
                + DBRACKET1 * ((LAMBDA_N+LAMBDA_K2-LAMBDA_K)/(2.0_dp*LAMBDA_N)) * LAMBDA_K * COUPLING_I_SWP &
                + (F_HALL_ARR/R2)*(DPSI_ALL(:,IDX_K2L2)*PSI_ALL(:,IDX_KL) + &
                    DPHI_ALL(:,IDX_KL)*CURV_PHI_ALL(:,IDX_K2L2)) * &
                    ((LAMBDA_K+LAMBDA_K2-LAMBDA_N)/2.0_dp) * COUPLING_I_STD &
                - DBRACKET3 * (LAMBDA_K2/LAMBDA_N) * COUPLING_J_STD &
                + (F_HALL_ARR/R2)*(PSI_ALL(:,IDX_KL)*CURV_PHI_ALL(:,IDX_K2L2) - &
                    DPSI_ALL(:,IDX_K2L2)*DPHI_ALL(:,IDX_KL)) * COUPLING_J_SWP
          END DO
    END DO
  END DO

  DEALLOCATE(PHI_ALL, PSI_ALL, DPHI_ALL, DPSI_ALL, CURV_PHI_ALL)
  DEALLOCATE(BRACKET1, BRACKET3, DBRACKET1, DBRACKET3, D2_CURV, R2, F_HALL_ARR)
  DEALLOCATE(ACTIVE_K, ACTIVE_L, ACTIVE_IDX)
END SUBROUTINE HALL_INDUCTION_RHS

END MODULE HALL_INDUCTION
