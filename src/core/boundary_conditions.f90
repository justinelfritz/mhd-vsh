MODULE BOUNDARY_CONDITIONS
!> Physics content of the boundary conditions applied to a magnetic
!> field carried as a poloidal/toroidal potential pair (Phi,Psi) --
!> FORTVSH's "standard" (J-coupled) VSH_POL_UP/VSH_POL_DN/VSH_TOR basis,
!> per Geppert & Wiebicke (1991) as used in the accompanying manuscript
!> (magfric/writeup_full/Jul2024/Jul2024.tex, Sec.3 "Choices of initial
!> and boundary conditions") -- each an (N_r) radial profile per (l,m),
!> the shape of a FIELD_TYPES SPECTRAL_SCALAR_T column. Deliberately
!> stops at mutating the row of a caller-owned radial operator; how that
!> operator is assembled into a full timestep system (time-discretization
!> coefficients, RHS forcing) is LINEAR_SOLVE's job, not yet written --
!> mirrors RADIAL_OPERATORS' own APPLY_ORIGIN_REGULARITY, which does the
!> same kind of in-place row replacement for the opposite (inner/origin)
!> boundary.
!>
!> **Vacuum (insulating) outer boundary, Model NS-A.** Matching the field
!> at r=r_out to a current-free exterior multipole expansion (the paper's
!> Eq.25, citing GW1) gives, independent of m:
!>   dPhi_nm/dr |_{r_out} = -(l/r_out) * Phi_nm(r_out)   (Robin)
!>   Psi_nm(r_out) = 0                                    (Dirichlet)
!> Unlike the polar (a_rad,a_pol,a_tor) basis this module's first version
!> was built against, Phi here is a single potential whose own radial
!> derivative carries independent physical information, so the vacuum
!> match is a genuine derivative (Robin) condition -- it needs the
!> existing D1 finite-difference stencil at the boundary row, not just an
!> algebraic value-to-value ratio.
USE KINDS,       ONLY: dp, i4
USE GRID_RADIAL, ONLY: RADIAL_GRID_T
IMPLICIT NONE
PRIVATE
PUBLIC :: APPLY_VACUUM_BC_POLOIDAL, APPLY_VACUUM_BC_TOROIDAL, APPLY_INNER_DIRICHLET_BC

CONTAINS

!> Overwrites row N (the outermost radial node, r=r_out) of a
!> caller-owned copy of Phi's first-derivative operator so that row
!> computes d(Phi)/dr + (l/r_out)*Phi instead of just d(Phi)/dr -- the
!> LHS of the vacuum Robin condition (see module header), meant to
!> replace whatever evolution-equation row would otherwise apply there,
!> with the corresponding system RHS forced to zero.
!>
!> D1's existing row N (built by RADIAL_OPERATORS' BUILD_RADIAL_OPERATORS)
!> is already the correct one-sided FD closure for d/dr at r_out, so no
!> stencil recomputation is needed here -- only the diagonal augmentation.
!>
!> @warning D1 is l-dependent through this augmentation (unlike the
!>   shared baseline RADIAL_OPERATOR_T%D1 it is meant to be copied from):
!>   pass a fresh per-l copy, the same convention RADIAL_OPERATORS'
!>   APPLY_ORIGIN_REGULARITY already uses for the inner boundary.
!>
!> @param D1 In/out first-derivative operator (a per-l copy of
!>   RADIAL_OPERATOR_T%D1), shape (N,N); row N is overwritten in place.
!> @param L Spherical-harmonic degree the condition applies to, l>=0.
!> @param RGRID Radial grid D1 was built on (a shell grid; r_out is
!>   RGRID%R(RGRID%N)).
SUBROUTINE APPLY_VACUUM_BC_POLOIDAL(D1, L, RGRID)
  REAL(KIND=dp),        INTENT(INOUT) :: D1(:,:)
  INTEGER(KIND=i4),     INTENT(IN)    :: L
  TYPE(RADIAL_GRID_T),  INTENT(IN)    :: RGRID
  INTEGER(KIND=i4) :: N
  N = RGRID%N
  D1(N,N) = D1(N,N) + REAL(L, KIND=dp)/RGRID%R(N)
END SUBROUTINE APPLY_VACUUM_BC_POLOIDAL

!> Overwrites row N (the outermost radial node, r=r_out) of a
!> caller-owned operator acting on Psi with the trivial Dirichlet row
!> Psi(r_out)=0 -- the vacuum condition's toroidal half (see module
!> header). Unlike APPLY_VACUUM_BC_POLOIDAL this doesn't depend on l or
!> on any existing FD content in that row, since it replaces the row
!> outright rather than augmenting it.
!>
!> @param OP In/out operator acting on Psi's radial profile (any (N,N)
!>   system matrix an evolution equation for Psi would otherwise place a
!>   row N in), overwritten in place.
!> @param RGRID Radial grid OP was built on (a shell grid; r_out is
!>   RGRID%R(RGRID%N)).
SUBROUTINE APPLY_VACUUM_BC_TOROIDAL(OP, RGRID)
  REAL(KIND=dp),       INTENT(INOUT) :: OP(:,:)
  TYPE(RADIAL_GRID_T), INTENT(IN)    :: RGRID
  INTEGER(KIND=i4) :: N
  N = RGRID%N
  OP(N,:) = 0.0_dp
  OP(N,N) = 1.0_dp
END SUBROUTINE APPLY_VACUUM_BC_TOROIDAL

!> Overwrites row 1 (the innermost radial node, r=r_min) of a
!> caller-owned operator with the trivial Dirichlet row f(r_min)=0 --
!> the "vanishing magnetic field at the core-crust boundary"
!> simplification, shared by both Phi and Psi since they get the
!> identical condition there (unlike the outer boundary, where Phi gets
!> the Robin condition and Psi gets Dirichlet -- see
!> APPLY_VACUUM_BC_POLOIDAL/TOROIDAL). Same row-replacement pattern as
!> those, just at the opposite boundary.
!>
!> @param OP In/out operator acting on Phi's or Psi's radial profile,
!>   overwritten in place.
SUBROUTINE APPLY_INNER_DIRICHLET_BC(OP)
  REAL(KIND=dp), INTENT(INOUT) :: OP(:,:)
  OP(1,:) = 0.0_dp
  OP(1,1) = 1.0_dp
END SUBROUTINE APPLY_INNER_DIRICHLET_BC

END MODULE BOUNDARY_CONDITIONS
