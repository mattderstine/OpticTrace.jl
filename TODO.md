# TODO

Running list of known issues to come back to. Not exhaustive -- add to it
as things are found; remove entries once fixed.

## Code issues

- `src/zemax.jl`: `ZemaxGeometry` struct is defined but never constructed
  anywhere -- `readZemax` returns its parsed data as a plain
  `(zsurfs, name, units, wavelengths)` tuple instead of wrapping it in a
  `ZemaxGeometry`. Either start using `ZemaxGeometry` as the return type,
  or remove it if it's not needed.
- `src/zemax.jl`: `readZemax`'s `basept`/`dir` keyword args are accepted
  but not actually used during parsing (the `basecurrent`/`dircurrent`
  variables they seed are only read by commented-out code). Either wire
  them up (the commented-out `zemaxsurfToSurface!`/`push!` lines suggest
  the original intent) or drop the unused args.
- `src/beamlet_decomposition.jl`: `gaussBeamParams` is an unimplemented
  stub (`function gaussBeamParams() end`) with a docstring describing a
  5-argument signature and return value that don't exist yet. See the
  `# TODO` comment already left in the file.
- `src/zemax.jl` (`zemaxsurfToSurface`): only Zemax surface `type`s
  `"STANDARD"` and `"EVENASPH"` are supported; any other type (e.g.
  toroidal, coordinate breaks) throws `error("Zemax surface type ...
  not implemented yet")`.
- `src/zemax.jl` / `docs/zemax_reference.md`: reading `.zar` archives
  (zipped Zemax file bundles) isn't implemented at all -- only the
  Python reference implementation exists, kept in
  `docs/zemax_reference.md` as a starting point for a future
  `readZemaxArchive`-style function.
- `src/surfaces.jl` (`lensASinglet`, `lensEASinglet`, lines ~518 & ~548):
  both have a `#ToDo` comment -- "should check if the input is really an
  asphere. if not make the surface spherical" -- that validation isn't
  implemented; passing non-aspheric coefficients silently proceeds as-is.
- **Bug** -- `src/surfaces.jl`, `lensSinglet`'s `order = "reverse"`
  branch: passes `Base.compute_assumed_setting` (a real but completely
  unrelated Julia compiler internal -- almost certainly a stray
  autocomplete/typo) as the `coating` argument to the second
  `refractSphere` call, instead of `coating`. Confirmed by direct
  testing: throws a `MethodError` (`coating::APorString` doesn't match
  a `Function`). Reachable through the public API
  (`lensSinglet(...; order="reverse")`, or via any of its four
  `lens_edmund.jl` callers with `order="reverse"`) but not exercised by
  any test.
- **Bug** -- `src/surfaces.jl`, `reflectOAConic`: its
  `attributesSurfaces` keyword's default value is written
  `attributesSurfaces = attributeSurfaces` (missing the second `s`),
  referencing an undefined variable. Confirmed by direct testing:
  calling `reflectOAConic(...)` without explicitly passing
  `attributesSurfaces` throws `UndefVarError: attributeSurfaces not
  defined`. `reflectOAConic` is exported, so this is reachable by any
  direct external caller; the only in-repo caller (`reflectOAP`) avoids
  it by always passing `attributesSurfaces` explicitly.
- `src/surfaces.jl`: `refractAsphere`'s `asphere` parameter is typed
  `AbstractVector{Float64}` (hardcoded), unlike its sibling
  `refractEvenAsphere`'s `asphere::AbstractVector{T}` (generic). Per
  `CLAUDE.md`'s stated convention ("numeric types are generally
  parameterized... rather than hardcoded to Float64... so functions
  stay compatible with ForwardDiff"), this hardcoding will silently
  break autodiff-based normal/gradient computations through
  `refractAsphere` specifically, unlike through `refractEvenAsphere`.
- `src/surfaces.jl`: `surfNormal(r::Point3{T}, s::NoProfile)`,
  `deltaToSurf(r::Ray{T}, p::NoProfile)`, and `modFunc(ray::Ray{T},
  normal::Vec3{T}, d::NoBendIndex)` are all effectively dead code --
  each has a more specific same-named method in `src/tracing.jl`
  (`NoProfile{T}`/`NoBendIndex{T}` tied to the ray's own type `T`) that
  Julia's dispatch always prefers when both apply, confirmed via
  `@which`. Not incorrect, just redundant -- low-priority cleanup
  candidates.
- `src/lens_thorlabs.jl`, `lens_ACL12708U(base, dir, wl)`: `wl` is
  actually used as a refractive index (passed as `rinOut`/`rinIn` to
  the two surface constructors), not a wavelength despite the name --
  misleading for any caller who reasonably expects to pass a
  wavelength like its siblings (`lens_TLF220APC`, `lens_TLF357775_405`)
  do. Not a crash bug, but worth fixing the parameter name/behavior for
  consistency. Also has no `order`/`lensname` keywords, unlike every
  other lens builder in the file.
- `src/lens_edmund.jl`, `lens_EO38398`: not exported (missing from this
  file's `export` line, unlike its three siblings
  `lens_EO68001`/`lens_EO67548`/`lens_EO67652`) -- likely an oversight;
  currently only reachable as `OpticTrace.lens_EO38398(...)`.
- `src/lens_thorlabs.jl`, `lensAC127019AB`: **also not exported**
  (missing from this file's two `export` lines, unlike every other
  builder in the file) -- found while writing `test/lens_catalogs.jl`
  (phase 9); only reachable as `OpticTrace.lensAC127019AB(...)`, same
  pattern as `lens_EO38398` above. Separately, unlike its two
  structurally identical siblings (`lensAC508180AB`, `lensAC127050A`),
  doesn't validate `order` -- any value other than exactly `"forward"` is
  silently treated as `"reverse"` instead of erroring on an
  unrecognized value.
- `src/mesh_primitives.jl`: `GeometryBasics.radius`/`widths` for
  `OptSurface` dispatch to `gbRadius`/`gbWidths`, which only have a
  method for the `(SizeLens, SurfProfileConic)` aperture/profile
  combination. Any `OptSurface` using a different aperture or a
  non-conic profile (sphere, asphere, even-asphere, cylinder, toroid,
  off-axis conic) will throw a `MethodError` when its mesh bounds are
  computed (e.g. for 3D plotting).
- **Bugs** -- `src/mesh_primitives.jl`, several confirmed (independently
  flagged by the IDE's linter as "Possible method call error"):
  - `gbWidths(a::Washer, p::NoProfile)` and `gbWidths(a::Disk,
    p::NoProfile)`: both reference `a.SemiDiameter`/`a.SemiDiamater`,
    neither of which is a real field (the actual field is
    `semiDiameter`; `SemiDiamater` is also a typo for `SemiDiameter`).
    Throws a field-access error if called.
  - `GeometryBasics.radius(c::Washer)`/`widths(c::Washer)` and the
    `Disk` equivalents: all four call `gbRadius`/`gbWidths` with
    `c.semiDiameter` (a bare number) as the first argument, but those
    functions expect the whole `Washer`/`Disk` object. No matching
    method exists for a number, so these throw a `MethodError`.
    Reachable from live plotting (`plotModelSurf!`/`Makie.mesh!`)
    whenever a `ModelSurface`'s aperture is a `RoundAperture` with a
    nonzero `obscure`, if Makie's mesh conversion calls these methods.
  - `gbRadius(a::RectAperture, profile::NoProfile)`: **correction**
    (found while writing `test/mesh_primitives.jl` for the phased
    test-writing plan, phase 4) to the branch description above -- the
    bug is actually in the *finite*-aperture branch (i.e. whenever
    `wclear`/`lclear` are **not** both infinite, the ordinary/everyday
    case for a `RectAperture`), not the `wclear == ∞ && lclear == ∞`
    branch. That both-infinite branch (`sqrt(a.wo^2+a.lo^2)`) works
    correctly; the finite branch references `a.clear`, which isn't a
    field of `RectAperture` (the real fields are `wclear`/`lclear`). No
    live code currently constructs an `OptSurface`/`ModelSurface` with a
    `RectAperture` aperture at all, so this is unreachable today, but it
    would break for the common finite-aperture case (not just an edge
    case) if that changed.
  - Likely fix for all of the above: `s/SemiDiameter/semiDiameter/`,
    `s/SemiDiamater/semiDiameter/`, `gbRadius(c, c.profile)` /
    `gbWidths(c, c.profile)` instead of `c.semiDiameter`, and
    `a.clear` -> `a.lclear`.
  - Also: `GeometryBasics.normals(s::AbstractSurface, nvertices=60)`
    calls `inOrOut(s)`, which only has a method for `OptSurface` --
    calling `normals` on any other `AbstractSurface` (e.g.
    `ModelSurface`) throws a `MethodError`, despite the generic
    `s::AbstractSurface` signature implying it should work for any
    surface type.
  - Minor/non-bug: `GeometryBasics.normals(s::OptSurface, nvertices=60)`
    is byte-for-byte identical to the generic
    `GeometryBasics.normals(s::AbstractSurface, nvertices=60)` method
    above it -- redundant, not incorrect.
- `src/plotting.jl` (~line 250-261): a non-mutating `trcAndPlotRay`
  (counterpart to `trcAndPlotRay!`) is commented out with the note "see
  if this method is needed" -- open question on whether to implement it.
- **Bug** -- `src/plotting.jl`, `perimeterRays`: on any ray that fails
  to trace (`status != 0`), it stores a `NaN` `Ray` and then does a bare
  `return` -- which returns `nothing`, discarding the whole `rays`
  vector (including rays already successfully traced) instead of
  continuing to the next perimeter angle. Almost certainly a `continue`
  was intended. `plotPerimeterRays`/`plotPerimeterRays!` (both variants)
  will fail if fed a `geo` where any perimeter ray misses, since they
  iterate over the `nothing` return value.
  - **Separate, more severe bug** -- found while writing
    `test/plotting.jl` (phase 12): `perimeterRays` never gets far enough
    to reach the bug above. Its `Ray(...)` call builds the direction
    argument as a raw `[cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)]` literal --
    a plain `Vector`, not a `Vec3` -- but `Ray` only has a constructor
    for `(Point{N,T}, Vec{N,T})`. This throws `MethodError`
    unconditionally, for every call regardless of input, before the
    per-ray tracing loop ever runs. Confirmed by direct testing.
    `perimeterRays`/`plotPerimeterRays`/`plotPerimeterRays!` (both
    variants) are all exported and all fail this way currently -- none
    of the four is currently callable at all. Likely fix: wrap both the
    base-offset and direction literals as `Vec3(...)`/use `Point3`
    arithmetic consistently, matching the rest of `src/`'s convention.
- **Bug** -- `src/plotting.jl`, `plotOPD!(scene, h::Float64,
  egeo::ExtendedGeometry; ...)`: its final line returns `θr, opdx, opdy`,
  but `θr` is never assigned anywhere in this method (it loops over
  `x = LinRange(-sizeP, sizeP, points)`, not an angle range) --
  `UndefVarError: θr not defined`, thrown only after all the actual
  tracing/plotting work has completed. No call site anywhere in
  `src/`/`test/`.
  - **Separate, more severe bug, shared with `plotOPD3D!`** -- found
    while writing `test/plotting.jl` (phase 12): both this method and
    `plotOPD3D!` never get far enough to reach the `θr` bug above (or,
    for `plotOPD3D!`, to complete at all). Both contain the identical
    line `dirRef = normalize!(usedgeo[1].base.base .- r)`. `.base.base`
    is a `Point3` (an immutable `StaticArrays`-style point), so
    `.- r` produces another immutable `Point3`; `normalize!` is the
    in-place variant and tries to `setindex!` into it, which `Point3`
    doesn't support. Throws `ErrorException("setindex!(::Point{3,
    Float64}, value, ::Int) is not defined. Hint: Use MArray or
    SizedArray...")` unconditionally, for any input, near the start of
    both methods -- confirmed by direct testing. Likely fix: `normalize`
    (non-mutating) instead of `normalize!`.
- **Likely bug** -- `src/plotting.jl`, `plotSpotDiagram(fig, spts;
  title = "Spot Diagram")` (the 2-arg method): calls `Axis(fig[1,1],
  title, tellwidth = false)`, passing `title` positionally. Its sibling
  5-arg method correctly writes `Axis(fig[1,1]; title, tellwidth =
  false)` (with the `;`). `Makie.Axis` doesn't accept a title
  positionally, so this is likely a missing-semicolon typo that throws
  a `MethodError`. No call site anywhere in `src/`/`test/`.
- `src/characterization.jl` (`traceLoss`, line ~236): docstring ends
  with a bare `TBW` ("to be written") placeholder -- description is
  incomplete/unfinished.
- `test/optics.jl` (line ~74): `sag` tests for `SurfProfileCyl` and
  `SurfProfileToroid` are written but disabled inside a `#= =#` block,
  with the comment "Tests not implemented for SurfProfileCyl &
  SurfProfileToroid". Worth revisiting -- unclear if they're disabled
  because the expected values are wrong or because the feature is
  incomplete.
- `src/lens_definitions.jl` (`ExtendedGeometry.geo` field, line ~220):
  typed as `Array{AbstractSurface}` with the comment "needs to be
  changed to `AbstractOpticalObject`" -- that abstract type already
  exists (`lens_definitions.jl:18`) but nothing uses it yet.
- `src/surfaces.jl` (`planeMirror`, line ~458): docstring notes it
  "could use `NoProfile` to speed things up but then would need to
  potentially overload other functions" -- a known, deliberately
  deferred performance optimization.
- `src/lens_definitions.jl` (`SurfProfileToroid`): only a `sag` method
  exists (`src/tracing.jl`), and it's explicitly code-commented as
  "likely incorrect" -- there is no `deltaToSurf`/`surfNormal` method
  for this type at all, so toroidal surfaces are defined but not
  actually traceable yet.
- **Bug** -- `src/tracing.jl`, `deltaToSurf(r::Ray{3,T}, profile::SurfProfileCyl{T})`
  (~line 292): the function body refers to a variable `p` throughout,
  but the parameter is named `profile` -- `p` is never defined, so
  calling this method throws `UndefVarError: p not defined`. Confirmed
  by the IDE's own linter ("Possible method call error"). Currently has
  no live call site anywhere (`SurfProfileCyl` is never constructed
  outside a disabled test block, see the `test/optics.jl` entry above),
  which is presumably why this hasn't been caught. Fix is presumably
  `s/p\./profile\./g` within that method.
- `src/tracing.jl`, `deltaToSurf(r::Ray{3,T}, p::SurfProfileOAConic{T})`
  (~line 275): preceded by the original author's own code comments
  (`# logic is flawed in this one` / `#change to add offset to ray to
  put it into the coordinate system of the offset parabola`) flagging
  it as believed-incorrect. Unlike the `SurfProfileCyl` bug above, this
  one *is* reachable (`SurfProfileOAConic` is constructed by
  `reflectOAConic` in `src/surfaces.jl`), so any live use of
  `reflectOAConic` may be tracing incorrectly. Needs investigation.
- **Bugs** -- `src/surface_manipulation.jl`, `reverseProfile!` (three
  confirmed field-name typos, plus one dispatch-fallback gap; nothing
  in this file has any call site anywhere, which is presumably why none
  of these have surfaced):
  - `reverseProfile!(profile::SurfProfileSphere)`: assigns
    `profile.curve` -- not a real field (`SurfProfileSphere`'s field is
    `curv`).
  - `reverseProfile!(profile::SurfProfileToroid)`: assigns
    `profile.curveX`/`profile.curveY` -- not real fields
    (`SurfProfileToroid`'s fields are `curvX`/`curvY`).
  - `reverseProfile!(profile::T) where T<:AbstractSurfProfile` (the
    generic fallback for profile types without their own specific
    method): assumes every such type has an `a` field. That's false for
    `SurfProfileOAConic` (`curv`/`ϵ`/`offset`) and `NoProfile` (`curv`
    only), both of which fall through to this method -- calling it on
    either throws a field-access error.
  - `reverseGeo`/`reverseSurface!`: despite `reverseGeo`'s
    `Vector{T} where T<:AbstractSurface` signature, `reverseSurface!`
    only has a method for `OptSurface` -- a `geo` containing a
    `ModelSurface` (e.g. built by `roundAperture`/`rectAperture`, both
    plausible in a real lens system with an aperture stop) throws a
    `MethodError`.
  - Separately (not a confirmed bug, just an open question worth
    checking): `reverseSurface!` leaves `surf.base.dir`/`surf.base.ydir`
    unchanged when reversing a surface -- only position and
    profile/index data are flipped. Whether propagation-direction
    reversal should also flip the surface's local orientation isn't
    obvious from the code alone.
- `src/lens_definitions.jl`: four abstract types appear to be unused
  scaffolding -- `AbstractRay`/`AbstractSurfBase` each have exactly one
  subtype (`Ray`/`SurfBase`) and are never used as a dispatch target
  anywhere; `AbstractOpticalObject`/`AbstractTrace` have zero subtypes
  at all. Low priority -- candidates for removal, or for actually being
  put to use (see the `ExtendedGeometry.geo` entry above for
  `AbstractOpticalObject`).

- **Bug** -- `src/characterization.jl`: `sizeOpticSurface` is exported
  (line 6) but never defined anywhere in the repo -- any call to it
  throws immediately (`UndefVarError`/no matching method). Found during
  the initial codebase survey for the phased test-writing plan (see
  "Test-writing plan" below); either implement it or remove it from the
  export list.
- **Bug** -- `src/tracing.jl`, `traceSurf`/`traceSurf!` (all four
  methods: `OptSurface` and `ModelSurface`, mutating and non-mutating):
  each checks `if delta == NaN` to detect a ray that missed the surface
  (`deltaToSurf` returns `NaN` on a miss). In IEEE 754, `NaN == NaN` is
  always `false`, so this branch is dead code -- status 1 ("missed") is
  unreachable through this check in all four methods. Confirmed by
  direct testing (test-writing plan, phase 5, `test/trace_geometry.jl`):
  a ray that should miss instead silently propagates `NaN` through the
  rest of the surface-normal/`modFunc`/`clipAperture` math and ends up
  reporting a false "success" (status 0) with a garbage `NaN` trace,
  rather than erroring or correctly reporting a miss. Likely fix:
  `isnan(delta)` instead of `delta == NaN`.
- **Bug** -- `src/extended_geo.jl`, `updateEGeo!`: named and documented
  as mutating its `ExtendedGeometry` argument ("calls the function
  necessary to create a static geometry for tracing", per its own
  docstring), but the current implementation only calls
  `defaultSetupGeo(...)` and returns the result -- it never assigns back
  into `egeo.geo`. Confirmed by direct testing (test-writing plan, phase
  1, `test/foundations.jl`). Likely fix: `egeo.geo =
  defaultSetupGeo(...)`.
- `src/lens_refractive_index.jl`: `dirBaseRefractiveIndex` is a
  hardcoded absolute path on the original author's machine
  (`/Users/matt/Development/Projects/refractiveindex/database/data`),
  not part of the repo. `test/refractive_index.jl` (phase 8),
  `test/lens_catalogs.jl` (phase 9), and (via an otherwise-unused
  `simplesystem` variable that calls `lens_TLAC254_060`) the
  pre-existing `test/optics.jl` all read real glass `.yml` files from
  this path with no override -- **this was a CI-breaking gap**: pointing
  `dirBaseRefractiveIndex` at a nonexistent directory and re-running the
  suite confirmed `getRefractiveIndexFunc` threw an uncaught
  `SystemError` ("No such file or directory") in all three files, each
  failing its whole `@testset`/file (not just an individual `@test`),
  which `.github/workflows/CI.yml`'s plain `Pkg.test()` against a fresh
  `actions/checkout@v3` clone would have hit every time (that path never
  exists there). **Fixed**: all three now check `HAS_GLASS_CATALOG`
  (`test/helper.jl`, `= isdir(OpticTrace.dirBaseRefractiveIndex)`) and
  skip (print an `@info`, don't fail) their catalog-dependent portion
  when it's false, verified by re-running the whole suite with
  `dirBaseRefractiveIndex` pointed at a nonexistent path: 0 errors, only
  the catalog-dependent tests (123 of them) skipped. This is a
  stopgap, not a real fix -- CI still gets zero coverage of the
  catalog-dependent code paths this way. Follow-up items that would
  actually close that gap (identified while building the test plan):
  - **Option: vendor a glass-file subset into `test/fixtures/`.** Same
    pattern already used for the Zemax fixture
    (`test/fixtures/test_singlet.zmx`) -- check the specific `.yml`
    files the current tests actually read into the repo and point
    `getRefractiveIndexFunc`/the lens builders at that directory instead
    of `OpticTrace.dirBaseRefractiveIndex` when running under CI. Lower
    effort than the full database (below); gives CI real coverage of
    the catalog-dependent code paths without depending on an external
    download. The exact file list currently exercised (11 files):
    `glass/schott/{N-SF11,N-BK7,N-LAK22,N-SF6,N-SF2,N-LAK10,N-SF57}.yml`,
    `glass/cdgm/{D-ZK3,D-LAK6}.yml`, `glass/hoya/{BAF11,E-FD10}.yml`
    (per each glass file's own license/redistribution terms on
    refractiveindex.info -- check before committing).
  - **Option: vendor/download the full refractive-index database** (e.g.
    from refractiveindex.info) so it's available to both automated tests
    and package users, instead of depending on a machine-local path that
    isn't part of the repo or a clean checkout. Larger effort than the
    subset option above (whole-database size/licensing to work out), but
    also fixes this for real package users, not just CI -- currently
    `dirBaseRefractiveIndex` is unusable for anyone except the original
    author regardless of testing.
  - Add a configuration-file mechanism to store the refractive-index
    database's location/state (replacing the hardcoded path), with room
    for other future package configuration alongside it -- useful either
    way, but especially once either vendoring option above exists and
    needs a location to point at.
  - **Process note**: this entry itself was wrong twice while being
    written -- first (phase 8) it correctly flagged `test/optics.jl` as
    also depending on this path, but gave no specifics; then a later
    pass "corrected" that to say `test/optics.jl` *doesn't* depend on it,
    based on a `grep` for the literal identifiers
    `dirBaseRefractiveIndex`/`getRefractiveIndexFunc` -- which missed
    the indirect dependency through `lens_TLAC254_060` (a lens-builder
    call that itself calls `getRefractiveIndexFunc` internally). Only
    actually simulating the missing-directory case (not grepping for
    identifier names) caught it. Grepping for a dependency's own name is
    not sufficient when the dependency can be reached indirectly through
    another function call.
- **Bug** -- `src/characterization.jl`, `traceMonteCarloRays`: its very
  first executable line builds `badray = Ray((NaN, NaN, NaN), (NaN, NaN,
  NaN))` from raw tuples, but `Ray` requires an actual `Point{N,T}`/
  `Vec{N,T}` pair -- this throws a `MethodError` unconditionally, before
  the function ever reaches its own arguments or the ray-tracing loop.
  Confirmed by direct testing (test-writing plan, phase 7,
  `test/characterization.jl`): no call to this function can currently
  succeed regardless of inputs. Likely fix: `badray = Ray(Point3(NaN,
  NaN, NaN), Vec3(NaN, NaN, NaN))`.

## Test-writing plan

A phased plan to give every real, reachable function in `src/`
functional test coverage (memory allocation/type-stability testing is
separate future work). Functions that are both unused and unexported, or
entirely non-functional, are excluded from the plan rather than given
tests -- most of the individual bugs and gaps this section and "Tests
needed" below describe are exactly why. Known-broken cases inside an
otherwise-working function are captured as `Test.@test_broken` alongside
normal passing tests, rather than left uncovered. This mostly supersedes
the "Tests needed" section below (which predates these exclusion-list
decisions) for the phases already done.

- [x] Phase 1 -- Foundations (`test/foundations.jl`)
- [x] Phase 2 -- Surface profile math (extended `test/optics.jl`)
- [x] Phase 3 -- Surface & lens builders (`test/surface_builders.jl`)
- [x] Phase 4 -- Mesh geometry (`test/mesh_primitives.jl`)
- [x] Phase 5 -- Full ray tracing (`test/trace_geometry.jl`)
- [x] Phase 6 -- Surface manipulation (`test/surface_manipulation.jl`)
- [x] Phase 7 -- Characterization / system analysis
      (`test/characterization.jl`)
- [x] Phase 8 -- Refractive-index catalog (`test/refractive_index.jl`)
- [x] Phase 9 -- Lens catalogs (`test/lens_catalogs.jl`)
- [x] Phase 10 -- Zemax import (`test/zemax.jl`)
- [x] Phase 11 -- Printing / text output (`test/printing.jl`)
- [x] Phase 12 -- Plotting (graphical output) (`test/plotting.jl`)

**All 12 phases complete.** Every real, reachable function in `src/`
now has functional test coverage (or an explicit exclusion reason, see
above); known-broken cases are captured as `Test.@test_broken`/
`@test_throws` rather than left uncovered. Memory allocation/
type-stability testing remains separate future work, as noted above.

## Docstring audit

Done. The src/-wide docstring review (characterization.jl, aperture.jl,
mesh_primitives.jl, printing.jl, lens_thorlabs.jl, surfaces.jl,
extended_geo.jl, zemax.jl, tracing.jl, plotting.jl, constants.jl,
lens_definitions.jl) is complete -- all found mismatches in *existing*
docstrings fixed, and constants.jl/lens_definitions.jl (previously
undocumented) now have docstrings throughout. `beamlet_decomposition.jl`
was intentionally skipped (see Code issues above).

## Docstring coverage pass (methods with no docstring at all)

A separate pass found 123 method/function definitions across 11 files
in `src/` with no docstring at all (as opposed to the audit above,
which was about *existing* docstrings being wrong). Policy for this
work: docstrings stay per-method, but each method's docstring must
describe its parameters consistently with its dispatch siblings (see
project memory `feedback-docstring-policy`); ambiguous cases get
flagged for the user rather than guessed at. Doing this file-by-file;
once every file below is done, a final pass checks consistency across
function-name families that span multiple files (`sag`/`deltaToSurf`/
`surfNormal`/`modFunc`/`gbRadius`/`gbWidths`, split across
`tracing.jl`/`surfaces.jl`/`mesh_primitives.jl`).

Test coverage in this project is almost nonexistent, so most functions
having no call sites elsewhere in the codebase is normal here, not a
sign of dead code (see project memory `project-optictrace-test-coverage`).
When this pass turns up a function/method with no call sites, that does
**not** get written into its docstring as "unused"/"dead"/"a candidate
for removal" -- it gets a plain, descriptive docstring, and a line in
"Tests needed" below instead (see project memory
`feedback-unused-functions`).

Caveat on the "123 methods" count: it came from a regex-based scan that
missed one-liner method definitions whose argument list contains a
default value with `=` before the closing paren (e.g.
`f(x; y = false) = ...`) -- confirmed by finding one such case by hand
in `plotting.jl` (`plotSurface3D!(scene, s::OptSurface;
transparency = false) = ...`) that the scan hadn't flagged. Each file
is still being read in full as it's done, so misses like this get
caught during the file's own pass, but the original per-file counts
below may be undercounts, and files already marked done were only as
thorough as that file's own full read -- not re-verified against a
fixed scanner.

- [x] `lens_refractive_index.jl` (13 methods)
- [x] `tracing.jl` (25 methods)
- [x] `mesh_primitives.jl` (22 methods)
- [x] `plotting.jl` (18 methods, +1 found by hand during the pass)
- [x] `surface_manipulation.jl` (10 methods)
- [x] `printing.jl` (8 methods)
- [x] `surfaces.jl` (8 methods)
- [x] `characterization.jl` (7 methods)
- [x] `lens_thorlabs.jl` (5 methods)
- [x] `lens_edmund.jl` (4 methods)
- [x] `aperture.jl` (3 methods)
- [x] final cross-file function-name-family consistency pass

**Coverage pass complete.** A final full-`src/` re-scan confirmed zero
remaining undocumented method definitions. That same final pass also
turned up an important Julia semantics gotcha (see project memory
`feedback-julia-docstring-binding`): a comment line sitting between a
docstring and its target definition silently breaks the doc binding
entirely (confirmed via `@doc` returning nothing despite the docstring
being right there in the source). Two docstrings written earlier this
session had exactly this problem (`printFigure` and `plotSurface3D!(
scene, s::OptSurface; ...)`, both in `src/plotting.jl`, where a
pre-existing comment had been left between the new docstring and the
function) -- both found and fixed by a dedicated re-scan, and verified
via `@doc` afterward. A second full-`src/` re-scan after the fix found
zero remaining cases of either problem.

## Tests needed

Functions/methods found to have no call sites anywhere in `src/`/`test/`
during the docstring coverage pass above. Not necessarily exhaustive or
prioritized -- originally a running list to work from once test-writing
started; now superseded item-by-item as the phased "Test-writing plan"
above covers each one (kept here, marked, rather than deleted, since the
original reasoning is still useful context for anyone reading this
after the fact).

- ✅ **Covered** (phase 1, `test/foundations.jl`) -- ~~`src/aperture.jl`:
  `roundAperture`, `rectAperture`, `clipAperture`, and `isAperture` have
  no direct unit test~~ (only exercised indirectly, if at all, through
  higher-level code) -- no bugs found on review this time, but worth a
  straightforward direct test since it's plain numeric code with no
  Makie dependency.
- ✅ **Covered** (phase 2, extended `test/optics.jl`) -- ~~`src/tracing.jl`:
  `sag`/`deltaToSurf`/`surfNormal` for `SurfProfileOAConic`~~ -- used in
  production by `reflectOAConic` (`src/surfaces.jl`) but has no
  dedicated unit test (relevant given the `deltaToSurf` correctness
  concern noted in Code issues above). `sag`/`surfNormal` now pass;
  `deltaToSurf` is `@test_broken` per the known logic-flaw bug.
- ✅ **Covered** (phase 2, extended `test/optics.jl`) -- ~~`src/tracing.jl`:
  `sag`/`deltaToSurf`/`surfNormal` for `SurfProfileCyl`, and `sag` for
  `SurfProfileToroid`~~ -- exercised only by the disabled block in
  `test/optics.jl` (see that entry in Code issues above); no production
  code constructs either profile type currently. `SurfProfileCyl`'s
  `sag`/`surfNormal` now pass (the old disabled block used an outdated
  formula and a stale expected value, replaced rather than just
  re-enabled); `SurfProfileToroid`'s `sag` is `@test_broken`.
  `SurfProfileCyl`'s `deltaToSurf` is excluded entirely (confirmed fully
  broken, not just a wrong value -- see Code issues above) rather than
  tested.
- ✅ **Covered** (phase 5, `test/trace_geometry.jl`) -- ~~`src/tracing.jl`:
  `traceGeometryRel`/`traceGeometryRel!` are exported and used
  elsewhere, but have no direct unit test of their own~~ (as opposed to
  being exercised indirectly through higher-level code).
- ✅ **Covered** (phase 4, `test/mesh_primitives.jl`) -- ~~`src/mesh_primitives.jl`:
  none of the `GeometryBasics`/`gb*` mesh overloads for `OptSurface`,
  `Washer`, `Disk`, or `RectAperture` have any test coverage at all~~ --
  this is presumably why the several bugs listed in Code issues above
  (typo'd fields, scalar-vs-object argument mixups) have gone unnoticed.
  **Correction**: this entry originally said testing these would need a
  Makie-capable test environment -- that turned out not to be true.
  `GeometryBasics.origin`/`radius`/`widths`/`coordinates`/
  `texturecoordinates`/`faces`/`normals` and the `gb*` functions are all
  plain function calls that just return points/vectors/face-index
  arrays; only actually creating a Makie `Figure`/`Scene` (`plotting.jl`,
  phase 12) needs a display. All of it was tested directly, no
  `xvfb-run` required.
- ✅ **Covered** (phase 12, `test/plotting.jl`) -- ~~`src/plotting.jl`:
  essentially nothing in this file has test coverage~~ -- run headless
  via `GLMakie.activate!(visible=false)` (persists across the file's own
  internal `activate!` calls, e.g. inside `multipleFigures`, since
  GLMakie merges into one config dict rather than resetting it), so
  `xvfb-run` wasn't needed here either, same as the `mesh_primitives.jl`
  correction above. Confirmed the three known bugs already listed in
  Code issues (`perimeterRays`'s early return, `plotOPD!`'s `θr` typo,
  `plotSpotDiagram`'s missing `;`), and found two more in the process
  (also now in Code issues): `perimeterRays` builds its `Ray` with a raw
  `Vector` direction instead of `Vec3`, so it and both
  `plotPerimeterRays`/`plotPerimeterRays!` variants throw `MethodError`
  unconditionally, regardless of the early-return bug; and
  `plotOPD!(scene, h, egeo::ExtendedGeometry)`/`plotOPD3D!` both call
  `normalize!` on an immutable `Point3`, throwing `ErrorException`
  unconditionally, before either method reaches its own later logic.
  `opdRel` (defined in this file, but exported from `src/printing.jl`'s
  `export` line instead of this file's) turned out to be plain
  OPD-difference math with no Makie dependency, and ended up covered by
  phase 11 (`test/printing.jl`) instead of waiting on phase 12.
- ✅ **Covered** (phase 6, `test/surface_manipulation.jl`) -- ~~`src/surface_manipulation.jl`:
  nothing in this file (`reverseGeo`, `reverseSurface!`, `reverseBase!`,
  `reverseProfile!`, `reverseMod!`, `thickGeo`) has any call site
  anywhere in `src/`/`test/`~~ -- this is presumably why the field-typo
  and dispatch-fallback bugs listed in Code issues above have gone
  unnoticed. All of these are plain numeric/struct-manipulation code
  with no Makie dependency, so were straightforward to unit-test
  directly (build a small `geo`, reverse it, check the expected sign
  flips and positions).
- ✅ **Covered** (phase 3, `test/surface_builders.jl`) -- ~~`src/surfaces.jl`:
  `lensSinglet` and its four `lens_edmund.jl` callers, and
  `reflectOAConic`/`reflectOAP`, have no test coverage~~ -- presumably
  why the `Base.compute_assumed_setting` and `attributeSurfaces` bugs
  listed in Code issues above have gone unnoticed. Both bugs are now
  caught via `@test_throws MethodError`/`@test_throws UndefVarError`
  respectively. (The `lens_edmund.jl` callers themselves are covered
  separately, phase 9.)
- ✅ **Covered** (phase 9, `test/lens_catalogs.jl`) -- ~~`src/lens_edmund.jl`:
  none of its four lens builders (`lens_EO38398`, `lens_EO68001`,
  `lens_EO67548`, `lens_EO67652`) have any test coverage~~ -- same
  glass-catalog-file dependency caveat as `lens_thorlabs.jl` below. All
  four forward straight through to `lensSinglet`'s broken
  `order="reverse"` branch (see Code issues above), now caught via
  `@test_throws MethodError` for each.
- ✅ **Covered** (phase 9, `test/lens_catalogs.jl`) -- ~~`src/lens_thorlabs.jl`:
  none of its lens-builder functions (`lensAC508180AB`, `lensAC127050A`,
  `lensAC127019AB`, `lens_ACL12708U`, `lens_TLF220APC`,
  `lens_TLF357775_405`, `lens_TLAC254_060`) have any test coverage~~ --
  tested end to end against the real local glass-catalog directory
  (`OpticTrace.dirBaseRefractiveIndex`), which turned out not to be a
  blocker since all the needed `.yml` files (schott N-SF11/N-BK7/
  N-LAK22/N-SF6/N-SF2/N-LAK10/N-SF57, cdgm D-ZK3/D-LAK6, hoya
  BAF11/E-FD10) are present on this machine. Found and documented
  `lensAC127019AB`'s missing export (see Code issues above) while
  writing these tests.
- ✅ **Covered** (phase 7, `test/characterization.jl`) -- ~~`src/characterization.jl`:
  `computeRearFocalPlane` and the 3-ray `surfClosestApproach` have no
  test coverage of their own~~ (unlike `findRFP`, which exercises the
  3-ray `surfClosestApproach` indirectly, or `distClosestApproach`,
  exercised indirectly via the 2-ray `surfClosestApproach`). Also now
  covers `distClosestApproach`'s documented unit-direction assumption.
  Phase 7 additionally found and fixed the entirely-broken
  `traceMonteCarloRays` (see Code issues above) while covering
  `traceMonteCarloRays`/`traceLoss`, which weren't separately called out
  in this list originally.
- `src/lens_refractive_index.jl`: `findRefractiveIndex`,
  `findRefractiveIndexAlt` -- fixed-arity dispersion formulas (same
  math as `riFormula1`/`riFormula2`, respectively, for `findRefractiveIndex`;
  a power-series formula for `findRefractiveIndexAlt`). Straightforward
  to unit-test directly against known coefficients. **Still excluded**
  by design: both are unused and unexported, which the "Test-writing
  plan" above's exclusion criteria deliberately skip (see that section's
  linked plan) -- not forgotten, just out of scope for now. `riFormula1`/
  `riFormula2`/`riFormula3` (the formulas these duplicate) are covered
  in phase 8.
