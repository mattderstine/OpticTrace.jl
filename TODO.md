# TODO

Running list of known issues to come back to. Not exhaustive -- add to it
as things are found. Once an item is fixed, don't delete or renumber it
here -- move its full entry to `FIXED.md`, keeping its original number,
and simply drop it from the list below. Item numbers in `Bugs`/`Code
issues` are stable, permanent identifiers: never reused, never
resequenced, so a number always refers to the same thing whether it's
still open (here) or resolved (`FIXED.md`). `Bugs` and `Code issues`
are separate numbering pools (each starts at 1) -- a new entry gets
`max(existing numbers for that same section, across this file and
FIXED.md) + 1`. Entries use a
bold `**N.**` prefix rather than real Markdown ordered-list syntax
(`1.`/`2.`/...) on purpose: CommonMark renderers (GitHub, VS Code's
Markdown preview) auto-renumber ordered lists sequentially from the
first item, silently hiding any gap left by a removed item -- a bold
prefix always displays the literal number written, in any renderer or
plain-text view.

## Bugs

Confirmed defects -- something throws, returns wrong data, or silently
does the wrong thing. See "Code issues" below for missing features,
cleanup candidates, and open design questions that aren't bugs. See
`FIXED.md` for bugs already resolved.

(none currently open)

## Code issues

Missing features, cleanup candidates, and open design questions -- not
confirmed bugs (see "Bugs" above for those). See `FIXED.md` for issues
already resolved.

- **1.** `src/zemax.jl`: `ZemaxGeometry` struct is defined but never
  constructed anywhere -- `readZemax` returns its parsed data as a
  plain `(zsurfs, name, units, wavelengths)` tuple instead of wrapping
  it in a `ZemaxGeometry`. Either start using `ZemaxGeometry` as the
  return type, or remove it if it's not needed.
- **2.** `src/zemax.jl`: `readZemax`'s `basept`/`dir` keyword args are
  accepted but not actually used during parsing (the
  `basecurrent`/`dircurrent` variables they seed are only read by
  commented-out code). Either wire them up (the commented-out
  `zemaxsurfToSurface!`/`push!` lines suggest the original intent) or
  drop the unused args.
- **3.** `src/beamlet_decomposition.jl`: `gaussBeamParams` is an
  unimplemented stub (`function gaussBeamParams() end`) with a
  docstring describing a 5-argument signature and return value that
  don't exist yet. See the `# TODO` comment already left in the file.
- **4.** `src/zemax.jl` (`zemaxsurfToSurface`): only Zemax surface
  `type`s `"STANDARD"` and `"EVENASPH"` are supported; any other type
  (e.g. toroidal, coordinate breaks) throws `error("Zemax surface type
  ... not implemented yet")`.
- **6.** `src/surfaces.jl` (`lensASinglet`, `lensEASinglet`, lines ~518
  & ~548): both have a `#ToDo` comment -- "should check if the input is
  really an asphere. if not make the surface spherical" -- that
  validation isn't implemented; passing non-aspheric coefficients
  silently proceeds as-is.
- **7.** `src/surfaces.jl`: `refractAsphere`'s `asphere` parameter is
  typed `AbstractVector{Float64}` (hardcoded), unlike its sibling
  `refractEvenAsphere`'s `asphere::AbstractVector{T}` (generic). Per
  `CLAUDE.md`'s stated convention ("numeric types are generally
  parameterized... rather than hardcoded to Float64... so functions
  stay compatible with ForwardDiff"), this hardcoding will silently
  break autodiff-based normal/gradient computations through
  `refractAsphere` specifically, unlike through `refractEvenAsphere`.
- **8.** `src/surfaces.jl`: `surfNormal(r::Point3{T}, s::NoProfile)`,
  `deltaToSurf(r::Ray{T}, p::NoProfile)`, and `modFunc(ray::Ray{T},
  normal::Vec3{T}, d::NoBendIndex)` are all effectively dead code --
  each has a more specific same-named method in `src/tracing.jl`
  (`NoProfile{T}`/`NoBendIndex{T}` tied to the ray's own type `T`)
  that Julia's dispatch always prefers when both apply, confirmed via
  `@which`. Not incorrect, just redundant -- low-priority cleanup
  candidates.
- **9.** `src/lens_thorlabs.jl`, `lens_ACL12708U(base, dir, wl)`: `wl`
  is actually used as a refractive index (passed as `rinOut`/`rinIn` to
  the two surface constructors), not a wavelength despite the name --
  misleading for any caller who reasonably expects to pass a
  wavelength like its siblings (`lens_TLF220APC`, `lens_TLF357775_405`)
  do. Not a crash bug, but worth fixing the parameter name/behavior
  for consistency. Also has no `order`/`lensname` keywords, unlike
  every other lens builder in the file.
- **10.** `src/lens_edmund.jl`, `lens_EO38398`: not exported (missing
  from this file's `export` line, unlike its three siblings
  `lens_EO68001`/`lens_EO67548`/`lens_EO67652`) -- likely an
  oversight; currently only reachable as
  `OpticTrace.lens_EO38398(...)`.
- **11.** `src/lens_thorlabs.jl`, `lensAC127019AB`: **also not
  exported** (missing from this file's two `export` lines, unlike
  every other builder in the file) -- found while writing
  `test/lens_catalogs.jl` (phase 9); only reachable as
  `OpticTrace.lensAC127019AB(...)`, same pattern as `lens_EO38398`
  above. Separately, unlike its two structurally identical siblings
  (`lensAC508180AB`, `lensAC127050A`), doesn't validate `order` -- any
  value other than exactly `"forward"` is silently treated as
  `"reverse"` instead of erroring on an unrecognized value.
- **12.** `src/mesh_primitives.jl`: `GeometryBasics.radius`/`widths`
  for `OptSurface` dispatch to `gbRadius`/`gbWidths`, which only have a
  method for the `(SizeLens, SurfProfileConic)` aperture/profile
  combination. Any `OptSurface` using a different aperture or a
  non-conic profile (sphere, asphere, even-asphere, cylinder, toroid,
  off-axis conic) will throw a `MethodError` when its mesh bounds are
  computed (e.g. for 3D plotting).
- **13.** `src/plotting.jl` (~line 250-261): a non-mutating
  `trcAndPlotRay` (counterpart to `trcAndPlotRay!`) is commented out
  with the note "see if this method is needed" -- open question on
  whether to implement it.
- **14.** `src/characterization.jl` (`traceLoss`, line ~236): docstring
  ends with a bare `TBW` ("to be written") placeholder -- description
  is incomplete/unfinished.
- **15.** `test/optics.jl` (line ~74): `sag` tests for `SurfProfileCyl`
  and `SurfProfileToroid` are written but disabled inside a `#= =#`
  block, with the comment "Tests not implemented for SurfProfileCyl &
  SurfProfileToroid". Worth revisiting -- unclear if they're disabled
  because the expected values are wrong or because the feature is
  incomplete.
- **16.** `src/lens_definitions.jl` (`ExtendedGeometry.geo` field, line
  ~220): typed as `Array{AbstractSurface}` with the comment "needs to
  be changed to `AbstractOpticalObject`" -- that abstract type
  already exists (`lens_definitions.jl:18`) but nothing uses it yet.
- **17.** `src/surfaces.jl` (`planeMirror`, line ~458): docstring notes
  it "could use `NoProfile` to speed things up but then would need to
  potentially overload other functions" -- a known, deliberately
  deferred performance optimization.
- **18.** `src/lens_definitions.jl` (`SurfProfileToroid`): only a `sag`
  method exists (`src/tracing.jl`), and it's explicitly
  code-commented as "likely incorrect" -- there is no
  `deltaToSurf`/`surfNormal` method for this type at all, so toroidal
  surfaces are defined but not actually traceable yet.
- **19.** `src/tracing.jl`, `deltaToSurf(r::Ray{3,T},
  p::SurfProfileOAConic{T})` (~line 275): preceded by the original
  author's own code comments (`# logic is flawed in this one` /
  `#change to add offset to ray to put it into the coordinate system
  of the offset parabola`) flagging it as believed-incorrect. This one
  *is* reachable (`SurfProfileOAConic` is constructed by
  `reflectOAConic` in `src/surfaces.jl`), so any live use of
  `reflectOAConic` may be tracing incorrectly. Needs investigation.
- **20.** `src/lens_definitions.jl`: four abstract types appear to be
  unused scaffolding -- `AbstractRay`/`AbstractSurfBase` each have
  exactly one subtype (`Ray`/`SurfBase`) and are never used as a
  dispatch target anywhere; `AbstractOpticalObject`/`AbstractTrace`
  have zero subtypes at all. Low priority -- candidates for removal, or
  for actually being put to use (see the `ExtendedGeometry.geo` entry
  above for `AbstractOpticalObject`).
- **21.** `src/lens_refractive_index.jl`: `dirBaseRefractiveIndex` is a
  hardcoded absolute path on the original author's machine
  (`/Users/matt/Development/Projects/refractiveindex/database/data`),
  not part of the repo. `test/refractive_index.jl` (phase 8),
  `test/lens_catalogs.jl` (phase 9), and (via an otherwise-unused
  `simplesystem` variable that calls `lens_TLAC254_060`) the
  pre-existing `test/optics.jl` all read real glass `.yml` files from
  this path with no override -- **this was a CI-breaking gap**:
  pointing `dirBaseRefractiveIndex` at a nonexistent directory and
  re-running the suite confirmed `getRefractiveIndexFunc` threw an
  uncaught `SystemError` ("No such file or directory") in all three
  files, each failing its whole `@testset`/file (not just an
  individual `@test`), which `.github/workflows/CI.yml`'s plain
  `Pkg.test()` against a fresh `actions/checkout@v3` clone would have
  hit every time (that path never exists there). **Fixed**: all
  three now check `HAS_GLASS_CATALOG` (`test/helper.jl`, `=
  isdir(OpticTrace.dirBaseRefractiveIndex)`) and skip (print an
  `@info`, don't fail) their catalog-dependent portion when it's
  false, verified by re-running the whole suite with
  `dirBaseRefractiveIndex` pointed at a nonexistent path: 0 errors,
  only the catalog-dependent tests (123 of them) skipped. This is a
  stopgap, not a real fix -- CI still gets zero coverage of the
  catalog-dependent code paths this way. Follow-up items that would
  actually close that gap (identified while building the test plan):
  - **Option: vendor a glass-file subset into `test/fixtures/`.**
    Same pattern already used for the Zemax fixture
    (`test/fixtures/test_singlet.zmx`) -- check the specific `.yml`
    files the current tests actually read into the repo and point
    `getRefractiveIndexFunc`/the lens builders at that directory
    instead of `OpticTrace.dirBaseRefractiveIndex` when running under
    CI. Lower effort than the full database (below); gives CI real
    coverage of the catalog-dependent code paths without depending on
    an external download. The exact file list currently exercised
    (11 files):
    `glass/schott/{N-SF11,N-BK7,N-LAK22,N-SF6,N-SF2,N-LAK10,N-SF57}.yml`,
    `glass/cdgm/{D-ZK3,D-LAK6}.yml`, `glass/hoya/{BAF11,E-FD10}.yml`
    (per each glass file's own license/redistribution terms on
    refractiveindex.info -- check before committing).
  - **Option: vendor/download the full refractive-index database**
    (e.g. from refractiveindex.info) so it's available to both
    automated tests and package users, instead of depending on a
    machine-local path that isn't part of the repo or a clean
    checkout. Larger effort than the subset option above
    (whole-database size/licensing to work out), but also fixes this
    for real package users, not just CI -- currently
    `dirBaseRefractiveIndex` is unusable for anyone except the
    original author regardless of testing.
  - Add a configuration-file mechanism to store the refractive-index
    database's location/state (replacing the hardcoded path), with
    room for other future package configuration alongside it --
    useful either way, but especially once either vendoring option
    above exists and needs a location to point at.
  - **Process note**: this entry itself was wrong twice while being
    written -- first (phase 8) it correctly flagged `test/optics.jl`
    as also depending on this path, but gave no specifics; then a
    later pass "corrected" that to say `test/optics.jl` *doesn't*
    depend on it, based on a `grep` for the literal identifiers
    `dirBaseRefractiveIndex`/`getRefractiveIndexFunc` -- which missed
    the indirect dependency through `lens_TLAC254_060` (a
    lens-builder call that itself calls `getRefractiveIndexFunc`
    internally). Only actually simulating the missing-directory case
    (not grepping for identifier names) caught it. Grepping for a
    dependency's own name is not sufficient when the dependency can
    be reached indirectly through another function call.
- **22.** `src/surface_manipulation.jl`, `reverseProfile!(profile::SurfProfileOAConic)`:
  negates `profile.curv` but leaves `profile.offset` untouched, so the
  reversed profile keeps describing the pre-reversal off-axis geometry
  instead of the mirrored one. Added (see `FIXED.md` #6) to fix the
  field-access crash the generic fallback previously hit on this type,
  but deliberately incomplete -- the function prints a message noting
  this at call time. Needs someone who understands `offset`'s sign/
  coordinate convention relative to `curv` to finish it properly.
- **24.** `src/zemax_browser.jl`, `_archiveContentPane`: the extraction
  output-path field is a plain editable text box pre-filled from
  `defaultExtractionOutputPath` -- there's no folder-picker UI, by
  deliberate choice, not oversight. A native `<input type="file">`
  picker can't work here: browsers withhold the real filesystem path
  from that input, and extraction runs server-side (needs a real path).
  If a picker is wanted later, two options were identified: (a) a
  server-side directory-tree picker panel reusing this same file's
  `walkZemaxDirectory`/tree-rendering code (no new dependency, correct
  regardless of whether the browser and the Bonito server are on the
  same machine -- recommended if this is revisited), or (b) shelling out
  to a native OS folder dialog from the server process (only correct
  when browser and server are the same machine, needs per-OS handling, a
  new dependency, and a no-op path for headless CI -- not recommended).

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
  `SurfProfileCyl`'s `deltaToSurf` had the same undefined-variable bug,
  now fixed and covered.
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
  correction above. Confirmed the known bugs already listed in Bugs
  (`perimeterRays`'s early return -- still open; `plotOPD!`'s `θr`
  typo and `plotSpotDiagram`'s missing `;` were two more, both now
  fixed, see `FIXED.md`), and found two more in the process (also now
  fixed, see `FIXED.md`): `perimeterRays` built its `Ray` with a raw
  `Vector` direction instead of `Vec3`, so it and both
  `plotPerimeterRays`/`plotPerimeterRays!` variants threw `MethodError`
  unconditionally, regardless of the early-return bug; and
  `plotOPD!(scene, h, egeo::ExtendedGeometry)`/`plotOPD3D!` both called
  `normalize!` on an immutable `Point3`, throwing `ErrorException`
  unconditionally, before either method reached its own later logic.
  `opdRel` (defined in this file, but exported from `src/printing.jl`'s
  `export` line instead of this file's) turned out to be plain
  OPD-difference math with no Makie dependency, and ended up covered by
  phase 11 (`test/printing.jl`) instead of waiting on phase 12.
- ✅ **Covered** (phase 6, `test/surface_manipulation.jl`) -- ~~`src/surface_manipulation.jl`:
  nothing in this file (`reverseGeo`, `reverseSurface!`, `reverseBase!`,
  `reverseProfile!`, `reverseMod!`, `thickGeo`) has any call site
  anywhere in `src/`/`test/`~~ -- this is presumably why the field-typo
  and dispatch-fallback bugs have gone unnoticed (the field-typo bugs,
  `SurfProfileSphere`/`SurfProfileToroid`, and the `ModelSurface`
  `reverseSurface!` dispatch gap, are all now fixed, see `FIXED.md`;
  the generic `reverseProfile!` fallback's gap remains open for
  `SurfProfileOAConic`, see Bugs above). All of these
  are plain numeric/struct-manipulation code with no Makie dependency,
  so were straightforward to unit-test directly (build a small `geo`,
  reverse it, check the expected sign flips and positions).
- ✅ **Covered** (phase 3, `test/surface_builders.jl`) -- ~~`src/surfaces.jl`:
  `lensSinglet` and its four `lens_edmund.jl` callers, and
  `reflectOAConic`/`reflectOAP`, have no test coverage~~ -- presumably
  why the `Base.compute_assumed_setting` and `attributeSurfaces` bugs
  listed in Code issues above went unnoticed. Both bugs are now fixed,
  with the tests updated from `@test_throws MethodError`/
  `@test_throws UndefVarError` to real passing assertions. (The
  `lens_edmund.jl` callers themselves are covered separately, phase 9.)
- ✅ **Covered** (phase 9, `test/lens_catalogs.jl`) -- ~~`src/lens_edmund.jl`:
  none of its four lens builders (`lens_EO38398`, `lens_EO68001`,
  `lens_EO67548`, `lens_EO67652`) have any test coverage~~ -- same
  glass-catalog-file dependency caveat as `lens_thorlabs.jl` below. All
  four forward straight through to `lensSinglet`'s `order="reverse"`
  branch, whose bug is now fixed; their `order="reverse"` tests were
  updated from `@test_throws MethodError` to real passing assertions.
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
  Phase 7 additionally found the entirely-broken `traceMonteCarloRays`
  and added `@test_throws` coverage documenting it (since fixed, see
  `FIXED.md` #11) while covering `traceMonteCarloRays`/`traceLoss`,
  which weren't separately called out in this list originally.
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
