# MHD-VSH Roadmap

Master to-do list. Update this alongside the code — it's the durable,
shared record of what's done, what's blocked on analytical work, and
what's scaffoldable next, so it doesn't only live in a chat session.

## Status

Working end-to-end: grid → radial operators → boundary conditions →
`linear_solve.f90` → `regime_interface.f90`/`timestepper.f90` →
`field_diagnostics.f90`, with the first real regime
(`diffusion_regime.f90`, pure linear Ohmic diffusion) running, tested,
and producing a time-resolved energy plot on every build. See
`project.md` (FORD-generated API docs, `ford project.md` →
`docs/index.html`) for the module-by-module reference.

## Analytical tasks (yours)

These block specific pieces of code from moving past "provisional" or
"placeholder." Nothing below is guessed at in the code — each is
explicitly flagged where it's used.

1. **Re-derive the `Phi`/`Psi` → `B_pol`/`B_tor` coefficients.**
   `field_diagnostics.f90`'s energy formula (and anything built on the
   same `R_l = r/sqrt(l(l+1))`, `VSH_POL_DN`/`VSH_POL_UP` split) is
   currently sourced from an unverified draft manuscript
   (`magfric/writeup_full/Jul2024/Jul2024.tex` Eq.396,400). Blocks:
   trusting energy numbers as physics, the Joule/Poynting decomposition
   below, and the `(Phi,Psi)->B` converter.
2. **Derive the Hall term** (`J x B`, fully nonlinear, small
   coefficient) for the weak electron-MHD regime. Currently an explicit
   placeholder in the (not-yet-built) Hall regime module.
3. **Derive the Joule dissipation + Poynting flux terms** (the energy
   *balance*, not just the reservoir `field_diagnostics.f90` already
   computes) -- needed for the full "every contribution to system
   energy" plot. The reference manuscript's own "Energy budget" section
   is an unfinished stub, so there's no shortcut to lean on; the
   Poynting flux specifically will need mode-coupling (Gaunt)
   integrals, a genuinely harder derivation than the energy reservoir
   was.
4. **Dynamo / ion-MHD limit analytical work.** Explicitly deferred --
   don't start the dynamo regime or resolve the kinematic-vs-full-MHD
   question until this is done.
5. *(Possibly joint)* the discrete `div(B)=0` self-consistency check --
   a numerical-discretization question (does the FD scheme preserve an
   identity that's exact in the continuous theory), deferred rather
   than guessed at.

## Programming tasks (scaffoldable now, independent of the above)

1. **Weak-Hall electron-MHD regime module** -- structure and IMEX time-split
   (implicit diffusion via `linear_solve.f90`, explicit Hall forcing)
   can be built now, reusing `diffusion_regime.f90`'s solve machinery,
   with the Hall term itself left as a placeholder per item 2 above.
2. **`(Phi,Psi) -> physical B` converter**, using `VSH_POL_UP_ALL`/
   `VSH_POL_DN_ALL` (formulas already identified, just not implemented)
   -- needed eventually for I/O/visualization and for evaluating
   nonlinear terms pseudospectrally.
3. **`io_*.f90`** -- proper checkpoint/restart and field-snapshot
   output. Only ad-hoc energy time-series logging exists today (the
   `ON_STEP` callback in `mhdvsh_diffusion.f90`).
4. **Force-free Bessel-Riccati initial conditions** -- unlike 1-4 above,
   this has a solid reference already (Igoshev, Elfritz & Popov 2016,
   MNRAS 462, 3689, Appendix A; arXiv:1608.08806), so it's an
   implementation task, not an open derivation: spherical Bessel
   functions `j_l`/`n_l` (FORTVSH provides none) + a root-finder for
   the transcendental eigenvalue equation.
5. **More artifacts as regimes come online** -- e.g. a radial
   energy-density profile plot, boundary-residual plots -- following
   the `build/artifacts/<regime>/<diagnostic>/{data,plots}/` convention
   in `CMakeLists.txt`'s "Artifacts" section.
6. **`induction`/`magnetofriction`/`mhd_full` regime scaffolding** as
   each is actually started.

## Infrastructure

- **Documentation**: `project.md` (FORD config) is set up and verified
  (`ford project.md` runs clean, zero warnings). Doc-comment style
  throughout matches FORTVSH's own `!>` convention.
- **CI / GitHub**: `.github/workflows/ci.yml` (build FORTVSH from
  source, build MHD-VSH against it, run `ctest`, upload
  `build/artifacts/**`) and `.github/workflows/docs.yml` (FORD build +
  GitHub Pages deploy, mirroring FORTVSH's own working setup almost
  exactly) are written and ready, but **inert until this becomes an
  actual git repository with a GitHub remote** -- there is no `.git`
  here yet.
