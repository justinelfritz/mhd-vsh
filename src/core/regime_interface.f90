MODULE REGIME_INTERFACE
!> The calling contract between a regime's own physics code and the
!> shared TIMESTEPPER: every regime provides a concrete ADVANCE
!> procedure matching REGIME_ADVANCE_I, fully opaque to the timestepper
!> -- it advances the regime's own state by one step of size DT however
!> its own physics requires (explicit forcing, any implicit solve via
!> LINEAR_SOLVE, ...), with none of that visible outside the regime.
!>
!> STATE is CLASS(*) (unlimited polymorphic) because each regime's state
!> is a genuinely different concrete type (e.g. one SPECTRAL_SCALAR_T
!> for a scalar diffusion regime, a Phi/Psi pair for a magnetic regime,
!> see FIELD_TYPES) -- this interface has to describe any of them
!> without knowing which. Each regime's own ADVANCE implementation does
!> a SELECT TYPE internally to recover its concrete state type; nothing
!> else in this codebase uses polymorphism, type-bound procedures, or
!> inheritance, and this interface doesn't require any of that either.
USE KINDS, ONLY: dp
IMPLICIT NONE
PRIVATE
PUBLIC :: REGIME_ADVANCE_I

ABSTRACT INTERFACE
  !> @param STATE The regime's own state, in/out, advanced in place.
  !> @param DT Step size.
  !> @param T Time at the start of this step (regimes with explicitly
  !>   time-dependent forcing can use this; most won't need it).
  SUBROUTINE REGIME_ADVANCE_I(STATE, DT, T)
    IMPORT :: dp
    CLASS(*),      INTENT(INOUT) :: STATE
    REAL(KIND=dp), INTENT(IN)    :: DT, T
  END SUBROUTINE REGIME_ADVANCE_I
END INTERFACE

END MODULE REGIME_INTERFACE
