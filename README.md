# MHD-VSH

A spherical MHD solver built on [FORTVSH](https://github.com/justinelfritz/FORTVSH)'s
vector spherical harmonics. Magnetic fields are advanced as a
poloidal/toroidal potential pair (`Phi`, `Psi`) rather than raw field
components, matched to a vacuum exterior at the outer boundary.

See `ROADMAP.md` for current status and what's next, and `project.md`
for the FORD-generated API documentation.

## Building

Requires a Fortran compiler (gfortran), CMake >= 3.14, LAPACK, and a
built [FORTVSH](https://github.com/justinelfritz/FORTVSH) install.

```bash
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/fortvsh-install
cmake --build build
cd build && ctest --output-on-failure
```

`cmake --build build` also regenerates the diagnostic plots under
`build/artifacts/<regime>/<diagnostic>/plots/` as a normal build target
-- see `CMakeLists.txt`'s "Artifacts" section.

### Regime drivers

Each regime is its own executable (no runtime dispatch):

```bash
./build/mhdvsh_diffusion [output_data_file]
```

## Generating documentation

```bash
pip install ford
ford project.md
```

Writes to `docs/`; open `docs/index.html`.

## Project layout

- `src/core/` -- regime-agnostic building blocks: grids, radial
  operators, transforms, boundary conditions, the linear solver, the
  regime/timestepper contract, diagnostics.
- `src/regimes/<name>/` -- one module per physics regime, each owning
  its own PDE assembly.
- `app/` -- one driver program per regime.
- `test/` -- CTest suite, one executable per module/regime.
- `scripts/` -- plotting scripts consuming `build/artifacts/` data.
