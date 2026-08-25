# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project overview

OpticTrace.jl is a Julia package for optical ray tracing of illumination
systems (as opposed to imaging-only raytracers): sequential/non-sequential
tracing of rays through lens surfaces, aperture handling, glass/refractive
index catalogs (Edmund, Thorlabs, Zemax import), mesh-based geometry, and
3D visualization via GLMakie.

## Environment

- Julia 1.12+, standard `Pkg` workflow (`Project.toml` / `Manifest.toml`).
- Instantiate deps: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Run the test suite: `julia --project=. -e 'using Pkg; Pkg.test()'`
  - Tests use GLMakie, which needs a display. CI runs headless via
    `xvfb-run` (see `.github/workflows/CI.yml`, which uses
    `xvfb-run -s '-screen 0 1024x768x24'` as the test-run prefix); do the
    same locally if running tests in a headless environment with no
    windowing system at all (e.g.
    `xvfb-run -s '-screen 0 1024x768x24' julia --project=. -e '...'`).
  - On a machine that already has a real display (e.g. a normal
    interactive session), `xvfb-run` isn't needed -- `test/plotting.jl`
    calls `GLMakie.activate!(visible=false)` itself before creating any
    figure, which renders everything off-screen (no window flashes, no
    blocking on `display(fig)`) and persists for the rest of the test
    run even through `plotting.jl`'s own internal `activate!(...)` calls
    (e.g. inside `multipleFigures`), since GLMakie merges into one
    config dict rather than resetting it on each call.

## Project structure

- `src/OpticTrace.jl` — module entry point; the `include(...)` order here
  defines load order and is significant (later files depend on types/
  functions defined earlier).
- `src/` is organized by concern, not mirrored 1:1 by test files:
  - `constants.jl`, `lens_definitions.jl` — core types (e.g. `Ray`,
    `SurfProfile*` types) and physical constants.
  - `tracing.jl`, `surfaces.jl`, `extended_geo.jl` — ray/surface
    intersection (`deltaToSurf`, `surfNormal`, `sag`), and
    `traceGeometry`/`traceSurf` (+ `!` in-place variants).
  - `aperture.jl`, `surface_manipulation.jl` — aperture and surface
    geometry helpers.
  - `mesh_primitives.jl` — `GeometryBasics`/`MeshIO` mesh construction.
  - `lens_refractive_index.jl`, `lens_edmund.jl`, `lens_thorlabs.jl` —
    glass/catalog data and refractive index models.
  - `zemax.jl` — Zemax file import.
  - `characterization.jl` — spot diagrams, system characterization.
  - `plotting.jl` — GLMakie-based visualization.
  - `printing.jl` — `Base.show`/pretty-printing for core types.
- `test/runtests.jl` includes `test/allocations.jl` plus one file per
  phase of a now-complete, phased test-writing effort (see `TODO.md`'s
  "Test-writing plan" section — all 12 phases are done): `foundations.jl`,
  `optics.jl`, `surface_builders.jl`, `mesh_primitives.jl`,
  `trace_geometry.jl`, `surface_manipulation.jl`, `characterization.jl`,
  `refractive_index.jl`, `lens_catalogs.jl`, `zemax.jl`, `printing.jl`,
  `plotting.jl`. Every real, reachable function in `src/` has functional
  coverage (or a documented exclusion reason, see `TODO.md`); known-broken
  cases are captured as `Test.@test_broken`/`@test_throws` rather than
  left uncovered. `test/testing.jl` is currently commented out there —
  check before assuming it runs. `test/helper.jl` has shared test helpers
  (e.g. `normal_from_sag` computed via `ForwardDiff` for validating
  analytic surface normals, and a `riFunc` stub for catalog-independent
  lens-builder tests). `test/fixtures/` holds checked-in test data (e.g.
  a minimal `.zmx` file for `zemax.jl`'s tests) — prefer this over
  depending on machine-local paths when a test needs a real file to read.
  `test/scratch.jl` is scratch/manual exploration, not wired into the
  test suite. `refractive_index.jl`, `lens_catalogs.jl`, and part of
  `optics.jl` (via `lens_TLAC254_060`) depend on real glass-catalog
  `.yml` files under `OpticTrace.dirBaseRefractiveIndex`, a hardcoded
  path on the original author's machine and not part of the repo (see
  `TODO.md`) — on a checkout without that directory (e.g. CI), those
  tests are skipped rather than failed: `test/helper.jl` defines
  `HAS_GLASS_CATALOG = isdir(OpticTrace.dirBaseRefractiveIndex)`, and
  each catalog-dependent testset checks it, printing an `@info` and
  skipping instead of running when it's `false`. This means CI currently
  gets **zero coverage** of catalog-dependent code paths — vendoring the
  glass database (or at least a subset, as `test/fixtures/` does for
  Zemax) would fix that but hasn't been done yet.
- No `scripts/` directory currently exists in this repo. A `docs/`
  directory does exist, but only holds `docs/zemax_reference.md` (a
  kept Python reference for a future `.zar`-archive-reading feature,
  see `TODO.md`) — there's no generated/Documenter.jl-style
  documentation site (unlike what
  `.github/instructions/copilot-instructions.md` implies).

## Code conventions actually used in this codebase

There is a Copilot instructions file at
`.github/instructions/copilot-instructions.md`. Treat it as background,
but note where it diverges from the code as it exists today:

- **Naming**: `src/` predominantly uses `camelCase` for functions
  (`deltaToSurf`, `surfNormal`, `traceGeometry`, `modFunc`) and
  `PascalCase` for types (`Ray`, `SurfProfileConic`, `OptSurface`).
  Some test/helper code uses `snake_case` instead (`normal_from_sag`),
  but that's inconsistency, not a second accepted style. **For new code,
  use `camelCase` for variables and functions, `PascalCase` for types and
  structs** (per `copilot-instructions.md`) — don't propagate the
  `snake_case` outlier into new functions.
- **Multiple dispatch over branching**: surface-specific behavior
  (`sag`, `deltaToSurf`, `surfNormal`) is implemented as dispatch across
  `SurfProfile*`/`OptSurface`/`ModelSurface` subtypes — add new surface
  types by adding new methods, not by branching on a type tag.
- **In-place variants**: performance-sensitive tracing functions have
  `!`-suffixed in-place counterparts (`traceGeometry!`, `traceSurf!`,
  `Trace!`). When adding tracing functionality, consider whether an
  in-place variant is warranted.
- **Type parameters for AD compatibility**: numeric types are generally
  parameterized (`Ray{3,T}`, `SurfProfileConic{U}`) rather than hardcoded
  to `Float64`, so functions stay compatible with `ForwardDiff` dual
  numbers. Preserve this when editing tracing/surface code — hardcoding
  `Float64` will silently break autodiff-based normal/gradient
  computations (see `test/helper.jl`'s `normal_from_sag`).
- **Docstrings**: a full `src/`-wide audit/coverage pass (see `TODO.md`)
  brought every method definition up to having its own docstring, so
  this is no longer the sparse-in-practice exception the Copilot
  instructions' blanket "docstrings on all public functions" rule might
  suggest needs fixing -- it's now actually true. Keep it that way for
  new code: docstrings stay per-method (not shared across a whole
  dispatch family) but should describe parameters consistently with
  their dispatch siblings; document non-obvious physical/geometric
  conventions (sign conventions, coordinate frames) since those are the
  easiest things to get subtly wrong. Watch for the binding gotcha a
  later re-scan turned up: a comment line left sitting between a
  docstring and its target definition silently breaks the doc binding
  entirely (`@doc` returns nothing despite the docstring being right
  there in the source) -- verify with `@doc` after adding one, don't
  just trust that the text is in the file.

## Testing

- Framework: `Test` stdlib, `@testset`/`@test` in `test/`.
- A common pattern for verifying analytic surface normals/intersections
  is to cross-check against `ForwardDiff`-computed gradients of the `sag`
  function (see `test/helper.jl`). Prefer this approach when adding
  tests for new surface types.
- Recent history (`git log`) shows this project fixing tracing bugs by
  adding regression tests for the specific surface/case that broke
  (e.g. aspheric traces) — follow that pattern for bug fixes here: add a
  regression test alongside the fix.
- `TODO.md`'s "Code issues" section lists this project's known,
  confirmed-by-testing bugs (many surfaced by writing the test suite
  itself) — check it before assuming odd behavior you hit is a bug you
  just found; it may already be tracked, with a `@test_broken`/
  `@test_throws` regression test already covering it.
- The phased test-writing effort that built out `test/` (see `TODO.md`'s
  "Test-writing plan") is complete: every real, reachable function in
  `src/` has functional coverage, or an explicit, documented reason it's
  excluded (unused *and* unexported, or entirely non-functional).
  Memory-allocation/type-stability testing is explicitly out of scope
  for that effort and remains open future work.

## Dependencies

Key deps (see `Project.toml`): `StaticArrays`, `GeometryBasics`,
`CoordinateTransformations` (core geometry/types), `GLMakie` (plotting),
`ForwardDiff` (autodiff, used both in the library and in tests),
`Optim`/`Roots` (numerical solving), `DataInterpolations`, `YAML`/`FileIO`/
`MeshIO` (data and mesh I/O), `StatsBase` (histogram binning behind
`rayHeatmap`/`rayHeatmap!`, `src/plotting.jl`). Don't add new
dependencies without updating `Project.toml`'s `[deps]` and `[compat]`.
`test/runtests.jl` also directly `using`s several of these (`GLMakie`,
`StatsBase`, `GeometryBasics`, `ForwardDiff`, `LinearAlgebra`,
`StaticArrays`) for use inside the test files themselves, not just
transitively through `OpticTrace` — add to that `using` list there if a
new test file needs direct access to one of these packages' own
exports/types.
