MODULE FIELD_DIAGNOSTICS
!> Time-resolved diagnostics computed directly from the (Phi,Psi)
!> potential representation -- never from a reconstructed physical B.
!> Currently: total magnetic energy and its radial density.
!>
!> @warning PROVISIONAL PHYSICS. The per-mode coefficients below (the
!>   R_l=r/sqrt(l(l+1)) rescaling, the split of Phi_lm/Phi_lm' between
!>   VSH_POL_DN/VSH_POL_UP, and hence the exact l(l+1) factors in
!>   MAGNETIC_ENERGY_DENSITY's formula) come from a draft manuscript
!>   (magfric/writeup_full/Jul2024/Jul2024.tex Eq.396,400) not yet
!>   independently re-derived/confirmed -- see BOUNDARY_CONDITIONS for
!>   the same source used for the outer BC. The *structure* (that
!>   FORTVSH's standard/J-coupled VSH basis VSH_TOR/VSH_POL_UP/
!>   VSH_POL_DN is complete and orthonormal, so the angle-integrated
!>   energy density is exactly a sum of |coefficient|**2 terms with no
!>   mode-coupling/Gaunt integrals needed) should still hold regardless;
!>   only the specific per-term coefficients are expected to change once
!>   re-derived. Update MAGNETIC_ENERGY_DENSITY's inner formula when that
!>   happens -- everything else in this module (the radial quadrature in
!>   TOTAL_MAGNETIC_ENERGY, the public interface) is independent of it.
!>
!> Current (provisional) derivation: B_pol(l,m,r) = (1/R_l**2)*
!> [Phi_lm*VSH_POL_DN(l,m) + R_l*Phi_lm'*VSH_POL_UP(l,m)], B_tor(l,m,r) =
!> -(i/R_l)*Psi_lm*VSH_TOR(l,m). Squaring and using orthonormality gives,
!> at each radius,
!>   integral |B|**2 dOmega = sum_lm l(l+1) * [ l(l+1)*|Phi_lm|**2/r**2
!>                                              + |Phi_lm'|**2 + |Psi_lm|**2 ]
!>
!> @warning The discrete analogue of div(B)=0 (a numerical-consistency
!>   check, since it's an identity by construction in the continuous
!>   theory) is deliberately not implemented yet -- it needs its own
!>   careful derivation in terms of Phi/Psi and their FD derivatives,
!>   not a quick add-on to this module.
USE KINDS,            ONLY: dp, i4
USE GRID_RADIAL,      ONLY: RADIAL_GRID_T
USE RADIAL_OPERATORS, ONLY: RADIAL_OPERATOR_T
USE FIELD_TYPES,      ONLY: SPECTRAL_SCALAR_T
USE VSH,              ONLY: YLM_INDEX
IMPLICIT NONE
PRIVATE
PUBLIC :: MAGNETIC_ENERGY_DENSITY, TOTAL_MAGNETIC_ENERGY

CONTAINS

!> @param PHI Poloidal potential (see module header), SPECTRAL_SCALAR_T.
!> @param PSI Toroidal potential, SPECTRAL_SCALAR_T; same N_R/LMAX as PHI.
!> @param OPS Radial derivative operator (for Phi'); OPS%D1 built on RGRID.
!> @param RGRID Radial grid PHI/PSI/OPS were built on.
!> Returns: angle-integrated |B|**2 at every radial node, size (RGRID%N).
!>   Integrate in r (times 1/2, or whatever unit convention) for total
!>   energy -- see TOTAL_MAGNETIC_ENERGY.
FUNCTION MAGNETIC_ENERGY_DENSITY(PHI, PSI, OPS, RGRID) RESULT(E_OF_R)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHI, PSI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp) :: E_OF_R(RGRID%N)
  COMPLEX(KIND=dp), ALLOCATABLE :: DPHI_DR(:,:)
  REAL(KIND=dp) :: LAMBDA_L
  INTEGER(KIND=i4) :: L, M, IDX, IR

  ALLOCATE(DPHI_DR(PHI%N_R, PHI%NLM))
  DPHI_DR = MATMUL(OPS%D1, PHI%COEF)

  E_OF_R = 0.0_dp
  DO L = 0, PHI%LMAX
    LAMBDA_L = REAL(L*(L+1), KIND=dp)
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      DO IR = 1, RGRID%N
        E_OF_R(IR) = E_OF_R(IR) + LAMBDA_L * ( &
          LAMBDA_L*ABS(PHI%COEF(IR,IDX))**2/RGRID%R(IR)**2 + &
          ABS(DPHI_DR(IR,IDX))**2 + ABS(PSI%COEF(IR,IDX))**2 )
      END DO
    END DO
  END DO
  DEALLOCATE(DPHI_DR)
END FUNCTION MAGNETIC_ENERGY_DENSITY

!> @param PHI, PSI, OPS, RGRID as MAGNETIC_ENERGY_DENSITY.
!> Returns: total magnetic energy, 1/2 * the radial (trapezoidal)
!>   integral of MAGNETIC_ENERGY_DENSITY; caller applies any further
!>   unit prefactor (4*pi, mu_0, ...) of their own convention.
FUNCTION TOTAL_MAGNETIC_ENERGY(PHI, PSI, OPS, RGRID) RESULT(E_TOTAL)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHI, PSI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp) :: E_TOTAL
  REAL(KIND=dp) :: E_OF_R(RGRID%N)
  INTEGER(KIND=i4) :: IR

  E_OF_R = MAGNETIC_ENERGY_DENSITY(PHI, PSI, OPS, RGRID)

  E_TOTAL = 0.0_dp
  DO IR = 1, RGRID%N-1
    E_TOTAL = E_TOTAL + 0.5_dp*(E_OF_R(IR)+E_OF_R(IR+1))*(RGRID%R(IR+1)-RGRID%R(IR))
  END DO
  E_TOTAL = 0.5_dp * E_TOTAL
END FUNCTION TOTAL_MAGNETIC_ENERGY

END MODULE FIELD_DIAGNOSTICS
