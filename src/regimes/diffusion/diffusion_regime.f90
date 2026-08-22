MODULE DIFFUSION_REGIME
!> Pure linear Ohmic diffusion of the (Phi,Psi) potential pair,
!> decoupled per (l,m): d(Phi_lm)/dt = ETA*(D2 - l(l+1)/r**2)*Phi_lm,
!> same equation for Psi_lm. The simplest possible regime -- no
!> advection, no Hall, no magnetofriction -- built first to prove the
!> whole pipeline (grid -> radial operators -> boundary conditions ->
!> linear_solve -> regime_interface/timestepper) end to end before any
!> nonlinear regime is attempted.
!>
!> Backward Euler in time: per l, A_l = I - dt*ETA*(D2-l(l+1)/r**2), with
!> the inner Dirichlet (BOUNDARY_CONDITIONS::APPLY_INNER_DIRICHLET_BC)
!> and outer vacuum (APPLY_VACUUM_BC_POLOIDAL for Phi,
!> APPLY_VACUUM_BC_TOROIDAL for Psi) rows spliced in. Since A_l depends
!> only on l (not m) and DT/ETA are fixed for a run, each A_l is
!> factorized once (DIFFUSION_INIT, called by the driver before handing
!> DIFFUSION_ADVANCE to TIMESTEPPER::RUN) and its LU factors reused for
!> every m at that l and every subsequent timestep -- the actual payoff
!> of LINEAR_SOLVE's factorize-once/solve-many design.
!>
!> A_l is a purely real matrix (FD stencils, curvature term, and BC rows
!> are all real), even though Phi/Psi's coefficients are complex, so
!> each complex solve is done as two real solves (real and imaginary
!> parts as separate right-hand-side columns of the same factorization)
!> rather than needing a complex LAPACK routine.
!>
!> The factorizations are cached in module-level (SAVE'd) state rather
!> than threaded through DIFFUSION_STATE_T, since REGIME_ADVANCE_I's
!> fixed signature (STATE, DT, T) has no room for anything beyond the
!> evolving state itself.
USE KINDS,               ONLY: dp, i4
USE GRID_RADIAL,         ONLY: RADIAL_GRID_T
USE RADIAL_OPERATORS,    ONLY: RADIAL_OPERATOR_T, ADD_CURVATURE_TERM
USE BOUNDARY_CONDITIONS, ONLY: APPLY_VACUUM_BC_POLOIDAL, APPLY_VACUUM_BC_TOROIDAL, &
                                APPLY_INNER_DIRICHLET_BC
USE LINEAR_SOLVE,        ONLY: DENSE_FACTORS_T, FACTORIZE_DENSE, SOLVE_FACTORED
USE FIELD_TYPES,         ONLY: SPECTRAL_SCALAR_T
USE VSH,                 ONLY: YLM_INDEX
IMPLICIT NONE
PRIVATE
PUBLIC :: DIFFUSION_STATE_T, DIFFUSION_INIT, DIFFUSION_ADVANCE

TYPE :: DIFFUSION_STATE_T
  TYPE(SPECTRAL_SCALAR_T) :: PHI, PSI
END TYPE DIFFUSION_STATE_T

TYPE(DENSE_FACTORS_T), ALLOCATABLE, SAVE :: PHI_FACTORS(:), PSI_FACTORS(:)

CONTAINS

!> One-time setup: builds and factorizes the per-l backward-Euler
!> systems for Phi and Psi (0<=l<=LMAX), caching them for every
!> subsequent DIFFUSION_ADVANCE call. Must be called once before handing
!> DIFFUSION_ADVANCE to TIMESTEPPER::RUN, with the same LMAX every
!> STATE%PHI/PSI in that run were allocated with. DT and ETA are baked
!> into the cached factorizations -- DIFFUSION_ADVANCE's own DT argument
!> is trusted to match what was passed here (TIMESTEPPER::RUN uses a
!> single fixed DT for an entire run, so this holds as long as
!> DIFFUSION_INIT was called with that same value).
!>
!> @param RGRID Radial grid (a shell grid; DIFFUSION_REGIME assumes an
!>   inner Dirichlet + outer vacuum boundary, so a full-sphere grid isn't
!>   supported here).
!> @param OPS Radial derivative operators built on RGRID.
!> @param LMAX Maximum spherical-harmonic degree to evolve.
!> @param ETA Resistivity (diffusivity) coefficient, used uniformly
!>   across all radial rows unless ETA_PROFILE is present.
!> @param DT Fixed timestep (must match every DT passed to
!>   DIFFUSION_ADVANCE for this initialization to remain valid).
!> @param ETA_PROFILE Optional per-radial-row resistivity (size RGRID%N,
!>   e.g. from CRUST_CONDUCTIVITY::ETA_AND_F_HALL_AT), overriding ETA
!>   row-by-row when present. Absent, behavior is bit-identical to the
!>   uniform-ETA case (test_diffusion_regime.f90 confirms this).
SUBROUTINE DIFFUSION_INIT(RGRID, OPS, LMAX, ETA, DT, ETA_PROFILE)
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  INTEGER(KIND=i4),        INTENT(IN) :: LMAX
  REAL(KIND=dp),           INTENT(IN) :: ETA, DT
  REAL(KIND=dp), OPTIONAL, INTENT(IN) :: ETA_PROFILE(:)
  REAL(KIND=dp), ALLOCATABLE :: A_PHI(:,:), A_PSI(:,:)
  INTEGER(KIND=i4) :: L, N, I

  N = RGRID%N

  IF (ALLOCATED(PHI_FACTORS)) DEALLOCATE(PHI_FACTORS)
  IF (ALLOCATED(PSI_FACTORS)) DEALLOCATE(PSI_FACTORS)
  ALLOCATE(PHI_FACTORS(0:LMAX), PSI_FACTORS(0:LMAX))

  ALLOCATE(A_PHI(N,N), A_PSI(N,N))
  DO L = 0, LMAX
    A_PHI = OPS%D2
    CALL ADD_CURVATURE_TERM(A_PHI, L, RGRID)
    IF (PRESENT(ETA_PROFILE)) THEN
      DO I = 1, N
        A_PHI(I,:) = -DT*ETA_PROFILE(I)*A_PHI(I,:)
      END DO
    ELSE
      A_PHI = -DT*ETA*A_PHI
    END IF
    DO I = 1, N
      A_PHI(I,I) = A_PHI(I,I) + 1.0_dp
    END DO
    A_PSI = A_PHI   ! identical PDE part; BC rows below make them differ

    A_PHI(N,:) = OPS%D1(N,:)
    CALL APPLY_VACUUM_BC_POLOIDAL(A_PHI, L, RGRID)
    CALL APPLY_INNER_DIRICHLET_BC(A_PHI)
    CALL FACTORIZE_DENSE(PHI_FACTORS(L), A_PHI)

    CALL APPLY_VACUUM_BC_TOROIDAL(A_PSI, RGRID)
    CALL APPLY_INNER_DIRICHLET_BC(A_PSI)
    CALL FACTORIZE_DENSE(PSI_FACTORS(L), A_PSI)
  END DO
  DEALLOCATE(A_PHI, A_PSI)
END SUBROUTINE DIFFUSION_INIT

!> Advances STATE (a DIFFUSION_STATE_T) by one backward-Euler step,
!> using the factorizations cached by DIFFUSION_INIT. Matches
!> REGIME_INTERFACE::REGIME_ADVANCE_I; DT and T are both unused (this
!> regime trusts DIFFUSION_INIT's baked-in DT rather than the per-call
!> one, and has no explicit time dependence) -- required by the shared
!> interface, not by this regime's own physics.
SUBROUTINE DIFFUSION_ADVANCE(STATE, DT, T)
  CLASS(*),      INTENT(INOUT) :: STATE
  REAL(KIND=dp), INTENT(IN)    :: DT, T
  ASSOCIATE (UNUSED_DT => DT, UNUSED_T => T); END ASSOCIATE

  SELECT TYPE (STATE)
  TYPE IS (DIFFUSION_STATE_T)
    CALL SOLVE_FIELD(STATE%PHI, PHI_FACTORS)
    CALL SOLVE_FIELD(STATE%PSI, PSI_FACTORS)
  END SELECT
END SUBROUTINE DIFFUSION_ADVANCE

!> Solves FIELD's per-l backward-Euler update in place, using FACTORS(l)
!> for each l -- shared by both Phi and Psi (only which FACTORS array is
!> passed differs). A_l is real, so each l's complex right-hand side
!> (one column per m, both Re and Im parts) is solved as one real
!> multi-column solve rather than needing a complex factorization; rows
!> 1 and N of the right-hand side are zeroed first since those rows of
!> A_l are homogeneous boundary constraints, not the PDE.
SUBROUTINE SOLVE_FIELD(FIELD, FACTORS)
  TYPE(SPECTRAL_SCALAR_T), INTENT(INOUT) :: FIELD
  TYPE(DENSE_FACTORS_T),   INTENT(IN)    :: FACTORS(0:)
  REAL(KIND=dp), ALLOCATABLE :: RHS(:,:), X(:,:)
  INTEGER(KIND=i4) :: L, M, IDX, NCOL, N

  N = FIELD%N_R
  DO L = 0, FIELD%LMAX
    NCOL = 2*L+1
    ALLOCATE(RHS(N, 2*NCOL), X(N, 2*NCOL))
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      RHS(:, M+L+1)      = REAL(FIELD%COEF(:,IDX), KIND=dp)
      RHS(:, NCOL+M+L+1) = AIMAG(FIELD%COEF(:,IDX))
    END DO
    RHS(1,:) = 0.0_dp
    RHS(N,:) = 0.0_dp

    CALL SOLVE_FACTORED(X, FACTORS(L), RHS)

    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      FIELD%COEF(:,IDX) = CMPLX(X(:,M+L+1), X(:,NCOL+M+L+1), KIND=dp)
    END DO
    DEALLOCATE(RHS, X)
  END DO
END SUBROUTINE SOLVE_FIELD

END MODULE DIFFUSION_REGIME
