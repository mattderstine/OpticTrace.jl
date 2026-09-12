# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project overview

OpticTrace.jl is a Julia package for optical ray tracing of illumination
systems (as opposed to imaging-only raytracers): tracing of rays through
lens surfaces either via the ordered `traceGeometry` walk over a fixed
`Vector{AbstractSurface}` or "non-sequentially" via manual, one-surface-
at-a-time `traceSurf` calls (**not** Zemax's `MODE NSC` sense of
non-sequential -- see the `zemax.jl` note under Project structure below;
this package's own tracer is sequential-only end to end), aperture
handling, glass/refractive index catalogs (Edmund, Thorlabs, Zemax
import), mesh-based geometry,
3D visualization via GLMakie, a Bonito.jl-based web UI
(`src/zemax_browser.jl`) for browsing a directory of `.zmx`/`.zar`/`.zmf`
Zemax files and extracting archive contents, and a generic, reusable
Bonito.jl file/directory picker component (`src/UItools/filepicker.jl`).

## Environment

- Julia 1.12+, standard `Pkg` workflow (`Project.toml` / `Manifest.toml`).
- Instantiate deps: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Run the test suite: `julia --project=. -e 'using Pkg; Pkg.test()'`
  - `./test.sh` wraps this exact invocation (full output to a temp log,
    last 150 lines printed, exit code preserved) as a fixed command so it
    doesn't need a fresh permission prompt each run (see
    `.claude/settings.json`) -- prefer it over ad hoc `julia -e '...'`
    invocations when just running the suite.
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
  - `agf.jl` — loads glasses from Zemax `.agf` glass-catalog text files
    (`loadAGFCatalog!`/`loadAGFCatalog`, paralleling
    `lens_refractive_index.jl`'s `loadRICatalog!`/`loadRICatalog`, and
    populating the same `defaultGlassCatalog` dict) by reusing the
    already-implemented `riFormula1`/`riFormula3` dispersion math rather
    than new formula code -- only AGF dispersion-formula codes 1
    (Schott) and 2 (Sellmeier1) are supported, since those are the only
    ones that are exact special cases of `riFormula3`/`riFormula1`
    respectively; other codes are skipped (see `TODO.md` #37).
  - `zemax.jl` — Zemax file import (`.zmx` parsing, `.zar`/`.zmf`
    archive reading and extraction). `readZemaxSystem` (returning an
    `OpticalSystem`, `lens_definitions.jl`) is the canonical entry
    point -- it owns the object-surface split and the infinite-object-
    distance case, rather than callers slicing `readZemax`'s raw
    `Vector{ZemaxSurf}` themselves. Supported surface `TYPE`s:
    `STANDARD`, `EVENASPH`, `TOROIDAL`, `COORDBRK`, `TILTSURF`,
    `ODDASPHE`, `XPOLYNOM`, `PARAXIAL` (see `TODO.md` bug #4 for
    what's still unsupported, e.g. `GRID_SAG`/`FZERNSAG`). This
    package's tracer is sequential-only end to end -- a Zemax `MODE
    NSC` (non-sequential) file is rejected outright by
    `readZemaxSystem`, not partially supported. A Zemax file's own
    length unit (`UNIT` -- `MM`/`CM`/`IN`/`M`) is converted to this
    package's canonical `LENGTH_UNIT` (mm, `constants.jl`) once, inside
    `readZemax` (via `convertZemaxUnitsToMM!`/`zemaxUnitToMM`), before
    any `ZemaxSurf`/`ZemaxHeader` value is used elsewhere -- so
    `OpticalSystem.geo`/`objectDistance`/`apertureValue` (when `ENPD`)/
    `fields` (when `zemaxFieldTypeIsHeight`) are always mm regardless of
    the source file's `UNIT`; `OpticalSystem.units`/`ZemaxHeader.units`
    itself is left as the *source* file's original unit string, kept as
    provenance only. Wavelengths (`WAVM`) are never unit-converted --
    always `WAVELENGTH_UNIT` (μm, `constants.jl`), independent of
    `UNIT` -- and any calculation combining a length with a wavelength
    (e.g. OPD-to-waves in `plotting.jl`) goes through the
    `LENGTH_TO_WAVELENGTH` constant rather than a hardcoded literal.
  - `zemax_browser.jl` — Bonito.jl web UI over `zemax.jl`'s functions:
    browsing `.zmx`/`.zar`/`.zmf` files via `UItools/filepicker.jl`'s
    `filePicker` (`:file` mode, filtered to those three extensions),
    content preview, and archive extraction -- including browsing for
    an extraction output directory via a "Browse..." toggle over
    another `filePicker` (`:directory` mode) in `archiveContentPane`
    (`zemaxBrowser` is the entry point). See the "Bonito/GLMakie name
    collision" note under Dependencies below before touching this
    file's imports.
  - `UItools/filepicker.jl` — a standalone, generic Bonito.jl
    file/directory picker component (`filePicker` is the entry point,
    `filePickerApp` a standalone-app wrapper for previewing/testing it),
    used by `zemax_browser.jl` (see above; this is `TODO.md`/`FIXED.md`
    item #24). Same "qualify every Bonito symbol" rule applies here too.
  - `characterization.jl` — spot diagrams, system characterization.
  - `plotting.jl` — GLMakie-based visualization.
  - `printing.jl` — `Base.show`/pretty-printing for core types.
- `test/runtests.jl` includes `test/allocations.jl` plus one file per
  phase of a now-complete, phased test-writing effort (originally
  planned as a "Test-writing plan" section in `TODO.md`, removed once
  all 12 phases were done -- see `git log -p -- TODO.md` if the
  original plan's wording is ever needed; `FIXED.md` entries mention
  individual phases in passing, e.g. "phase 8", "phase 12"):
  `foundations.jl`,
  `optics.jl`, `surface_builders.jl`, `mesh_primitives.jl`,
  `trace_geometry.jl`, `surface_manipulation.jl`, `characterization.jl`,
  `refractive_index.jl`, `lens_catalogs.jl`, `zemax.jl`,
  `zemax_browser.jl`, `filepicker.jl`, `printing.jl`, `plotting.jl`
  (`zemax_browser.jl` and `filepicker.jl` were added after the phased
  plan completed, not part of it, but follow the same "every reachable
  function gets coverage" bar). Every real,
  reachable function in `src/` has functional
  coverage (or a documented exclusion reason, see `TODO.md`); known-broken
  cases are captured as `Test.@test_broken`/`@test_throws` rather than
  left uncovered. `test/testing.jl` is currently commented out there —
  check before assuming it runs. `test/helper.jl` has shared test helpers
  (e.g. `normal_from_sag` computed via `ForwardDiff` for validating
  analytic surface normals, and a `riFunc` stub for catalog-independent
  lens-builder tests). `test/fixtures/` holds checked-in test data: a
  minimal `.zmx` file for `zemax.jl`'s tests, plus synthetic
  `test_archive.zar`/`test_catalog.zmf` (built by
  `test/fixtures/generate_synthetic_zemax.jl`, a maintenance script *not*
  wired into `runtests.jl` — rerun it by hand if the `.zar`/`.zmf` byte
  layout understanding ever changes) so `zemax_browser.jl`'s archive
  tests get real CI coverage instead of depending on the machine-local
  `HAS_ZAR_SAMPLE`/`HAS_ZMF_SAMPLE` samples `zemax.jl`'s own tests still
  use (see below) — prefer checked-in fixtures like these over depending
  on machine-local paths when a test needs a real file to read.
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
- No top-level `scripts/` directory exists in this repo (the closest
  thing, `test/fixtures/generate_synthetic_zemax.jl`, is a test-fixture
  maintenance script, not a general scripts area). `docs/` holds a
  minimal Documenter.jl site (`docs/make.jl`, `docs/src/index.md`,
  `docs/src/api.md`, plus `docs/Project.toml`; `docs/build/` and
  `docs/Manifest.toml` are gitignored local build output), which
  auto-generates its API reference page from `src/`'s own docstrings via
  `@autodocs` rather than listing functions by hand -- so it needs no
  maintenance as functions are added/removed. Build it locally with
  `julia docs/make.jl` (`Pkg.instantiate()`s the `docs/` environment
  itself); it's not wired into CI or published anywhere yet. (There used
  to also be a `docs/references/zemax_reference.md` -- a ported Python
  reference for `.zar`/`.zmf` binary-format parsing, from back before
  that parsing was implemented in `src/zemax.jl`, see `FIXED.md` #5/#23
  -- removed as no longer needed now that `src/zemax.jl` is the working
  implementation.)

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
  `snake_case` outlier into new functions. Functions should not be
  prefixed with a leading underscore (`_foo`) to signal "private"/
  internal — Julia gives leading underscores no special meaning, so this
  repo does not use that convention; keep internal helper names the same
  `camelCase` as public functions.
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
- `TODO.md`'s "Bugs" section lists this project's known,
  confirmed-by-testing defects (many surfaced by writing the test suite
  itself; "Code issues" is a separate section for missing
  features/cleanup/open design questions, not confirmed bugs) — check
  it before assuming odd behavior you hit is a bug you just found; it
  may already be tracked, with a `@test_broken`/`@test_throws`
  regression test already covering it. Both sections use a stable,
  never-renumbered `**N.**` item-number convention (see `TODO.md`'s own
  header for why, and `FIXED.md`, which archives resolved entries under
  their original number rather than deleting them) — don't renumber
  surrounding items when removing or adding one.
- The phased test-writing effort that built out `test/` (see the note
  under "Project structure" above -- its planning section has since
  been removed from `TODO.md`) is complete: every real, reachable
  function in
  `src/` has functional coverage, or an explicit, documented reason it's
  excluded (unused *and* unexported, or entirely non-functional).
  Memory-allocation/type-stability testing is explicitly out of scope
  for that effort and remains open future work.
- **Testing Bonito UI code** (`test/zemax_browser.jl`, `test/filepicker.jl`): split coverage
  into (1) plain unit tests of any logic layer that doesn't touch Bonito
  at all (dispatch, path handling, struct-building — the bulk of the
  coverage, cheapest to write and to trust), and (2) Bonito "smoke"
  tests that build the real `Bonito.App`/DOM nodes and render them via
  `Bonito.export_static` to a temp HTML file, then assert on substrings
  in the rendered output. That second layer only checks the *initial*
  render, not click-driven interactivity (button callbacks, live
  Observable updates) — deliberately out of scope for this project's
  automated suite; a real headless-browser-driven test would be the only
  way to cover that and isn't worth the added CI dependency/flakiness at
  this project's size. No `xvfb-run` is needed for either layer: Bonito
  is a plain HTTP/WebSocket server, unlike GLMakie's GLFW/OpenGL
  requirement.

## Dependencies

Key deps (see `Project.toml`): `StaticArrays`, `GeometryBasics`,
`CoordinateTransformations` (core geometry/types), `GLMakie` (plotting),
`Bonito` (the `zemax_browser.jl` web UI; Julia 1.12+ floor, see
"Environment," is partly driven by `Bonito@5`'s own `julia = "1.11"`
compat floor), `ForwardDiff` (autodiff, used both in the library and in
tests), `Optim`/`Roots` (numerical solving), `DataInterpolations`,
`YAML`/`FileIO`/`MeshIO` (data and mesh I/O), `StatsBase` (histogram
binning behind `rayHeatmap`/`rayHeatmap!`, `src/plotting.jl`),
`BenchmarkTools` (used directly by `test/runtests.jl` and
`test/scratch.jl`, not by any file under `src/`). Don't add
new dependencies without updating `Project.toml`'s `[deps]` and
`[compat]`.
`test/runtests.jl` also directly `using`s several of these (`GLMakie`,
`StatsBase`, `GeometryBasics`, `ForwardDiff`, `LinearAlgebra`,
`StaticArrays`, `BenchmarkTools`) for use inside the test files
themselves, not just transitively through `OpticTrace` — add to that
`using` list there if a new test file needs direct access to one of
these packages' own exports/types.

`Project.toml` also lists `IterTools`, `using`'d in
`src/OpticTrace.jl` -- but no file under `src/` actually calls any of
its functions (see `TODO.md` Code issues #25); treat it as unused
rather than as a dependency backing any current feature.

**Bonito/GLMakie name collision — never add `using Bonito` to a shared
`using` block.** GLMakie and Bonito both export several identical names
bound to unrelated types: `Button`, `Slider`, `Checkbox`, `Dropdown`.
Confirmed directly: after `using GLMakie, Bonito` in the same scope, a
bare `Button` throws `UndefVarError` (Julia leaves a genuinely
conflicting exported name unresolved rather than picking one) — and
`src/plotting.jl`'s `multipleFigures` already uses GLMakie's unqualified
`Button`, so a top-level `using Bonito` in `src/OpticTrace.jl` (which
already has `using GLMakie`) would break it. `src/OpticTrace.jl` instead
does `import Bonito` (not `using`) once, module-wide; `src/zemax_browser.jl`
and `src/UItools/filepicker.jl` rely on that binding (they're `include`'d
into the same module) rather than importing it themselves, and both
qualify every reference (`Bonito.App`, `Bonito.DOM`, `Bonito.Button`,
`Bonito.TextField`, `Bonito.Server`, `Bonito.on`, `Bonito.Observable`,
...). `test/zemax_browser.jl` and `test/filepicker.jl` are separate
top-level scripts, not `include`'d into the module, so each does its own
local, testset-scoped `import Bonito` for the identical reason — not
added to `test/runtests.jl`'s shared `using` block. Keep this pattern for
any future Bonito-dependent code in this package.

**Before writing or reviewing Bonito code, check Bonito.jl's own
`AGENTS.md`** (`https://github.com/SimonDanisch/Bonito.jl/blob/master/AGENTS.md`,
fetch the raw file — `.../raw/master/AGENTS.md` — for the full text) for
architecture guidance: widget/state-ownership patterns, the three-tier
Julia↔JS communication model, session-scoped `on`/`map` (deregister
listeners when a session closes, not the bare forms), `Bonito.Styles`/
`CSS` instead of inline `style="..."` strings, and its anti-pattern
checklist. `src/zemax_browser.jl` and `src/UItools/filepicker.jl` are
this repo's only Bonito-dependent code.
