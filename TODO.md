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

- **14.** `src/zemax.jl`, `zemaxUnitToMM(units::String)`: only accepts
  the abbreviated Zemax `UNIT` tokens `"MM"`/`"CM"`/`"IN"`/`"M"` --
  `readZemax` (via `convertZemaxUnitsToMM!`) throws an `ArgumentError`
  on any `.zmx` file whose `UNIT` line spells the unit out as `"METER"`
  instead of `"M"`. Confirmed against real sample files: at least three
  files under the local Zemax install's `Samples/` tree (`Non-
  sequential/Miscellaneous/Multiple mirror telescope.zmx`, `Sequential/
  Image Simulation/Example 4, a diffraction limited system.ZMX`,
  `Sequential/Telescopes/Hubble.zmx`) use `UNIT METER` and fail to
  parse with the current code. Fix: add a `units == "METER" && return
  1000.0` branch (same factor as `"M"`).

## Code issues

Missing features, cleanup candidates, and open design questions -- not
confirmed bugs (see "Bugs" above for those). See `FIXED.md` for issues
already resolved.

- **3.** `src/beamlet_decomposition.jl`: `gaussBeamParams` is an
  unimplemented stub (`function gaussBeamParams() end`) with a
  docstring describing a 5-argument signature and return value that
  don't exist yet. See the `# TODO` comment already left in the file.
- **4.** `src/zemax.jl` (`zemaxsurfToSurface`): only Zemax surface
  `type`s `"STANDARD"`, `"EVENASPH"`, `"TOROIDAL"` (see `FIXED.md` #18),
  `"COORDBRK"`, `"TILTSURF"`, `"ODDASPHE"`, `"XPOLYNOM"`, and
  `"PARAXIAL"` are supported; any other type (e.g. `"GRID_SAG"`,
  `"FZERNSAG"`) throws `error("Zemax surface type ... not implemented
  yet")`.
- **6.** `src/surfaces.jl` (`lensASinglet`, `lensEASinglet`, lines ~578
  & ~608): both have a `#ToDo` comment -- "should check if the input is
  really an asphere. if not make the surface spherical" -- that
  validation isn't implemented; passing non-aspheric coefficients
  silently proceeds as-is.
- **9.** `src/lens_thorlabs.jl`, `lens_ACL12708U(base, dir, wl)`: `wl`
  is actually used as a refractive index (passed as `rinOut`/`rinIn` to
  the two surface constructors), not a wavelength despite the name --
  misleading for any caller who reasonably expects to pass a
  wavelength like its siblings (`lens_TLF220APC`, `lens_TLF357775_405`)
  do. Not a crash bug, but worth fixing the parameter name/behavior
  for consistency. Also has no `order`/`lensname` keywords, unlike
  every other lens builder in the file.
- **13.** `src/plotting.jl` (~line 325-336): a non-mutating
  `trcAndPlotRay` (counterpart to `trcAndPlotRay!`) is commented out
  with the note "see if this method is needed" -- open question on
  whether to implement it.
- **14.** `src/characterization.jl` (`traceLoss`, line ~314): docstring
  ends with a bare `TBW` ("to be written") placeholder -- description
  is incomplete/unfinished.
- **16.** `src/lens_definitions.jl` (`ExtendedGeometry.geo` field, line
  ~784): typed as `Array{AbstractSurface}` with the comment "needs to
  be changed to `AbstractOpticalObject`" -- that abstract type
  already exists (`lens_definitions.jl:93`) but nothing uses it yet.
- **17.** `src/surfaces.jl` (`planeMirror`, line ~520): docstring notes
  it "could use `NoProfile` to speed things up but then would need to
  potentially overload other functions" -- a known, deliberately
  deferred performance optimization.
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
- **25.** `src/OpticTrace.jl`: `using IterTools` (line 8) appears to be
  an unused import -- no file under `src/` calls any `IterTools`
  function (confirmed by grepping for common ones: `partition`,
  `product`, `chain`, `groupby`, `subsets`, `distinct`, etc., all with
  zero hits). `test/testing.jl` (commented out of `runtests.jl`, see
  "Project structure" in `CLAUDE.md`) also has its own `using
  IterTools`, so it isn't even exercised indirectly by the disabled
  file. Low priority -- candidate for removal from both `src/OpticTrace.jl`
  and `Project.toml`'s `[deps]`, unless it was meant to back some
  not-yet-written feature.
- **26.** `src/zemax.jl`, `readZemax`'s `FTYP` parsing: only the
  field-type code (`FTYP`'s first token) is decoded into
  `ZemaxHeader.fieldType`/`OpticalSystem.fieldType`. `FTYP`'s other
  tokens (telecentricity, a field-count value, an afocal-image-space
  flag, among others -- exact bit positions not confirmed against
  OpticStudio's own file-format reference) aren't decoded. Found while
  scanning real sample files: at least one multi-field file
  (`sc_dbga1.zmx` from the local Zemax install's `Samples/Short
  course/...` directory) writes `XFLN`/`YFLN`/`FWGN` with more entries
  than the system's true field count, apparently padding a fixed-size
  array -- so `OpticalSystem.fields`/`fieldWeight` may include spurious
  trailing entries for such a file, since nothing currently truncates
  them to `FTYP`'s real count. `readZemax`'s `XFLN`/`YFLN`/`FWGN`
  parsing takes every token on those lines as-is; needs revisiting once
  `FTYP`'s field-count token position is confirmed.
- **28.** `src/lens_refractive_index.jl`: `findRefractiveIndex`,
  `findRefractiveIndexAlt` -- fixed-arity dispersion formulas (same
  math as `riFormula1`/`riFormula2`, respectively, for
  `findRefractiveIndex`; a power-series formula for
  `findRefractiveIndexAlt`) have no test coverage. Straightforward to
  unit-test directly against known coefficients, but deliberately
  excluded from this project's completed functional-test-coverage
  effort, since both are unused and unexported. Not forgotten, just
  out of scope for now -- `riFormula1`/`riFormula2`/`riFormula3` (the
  formulas these duplicate) are covered in `test/refractive_index.jl`.
- **29.** `src/zemax.jl`: Zemax surface `TYPE`s `GRID_SAG`, `FZERNSAG`,
  `IRREGULA`, `FZERNPHA` (Zernike/measured sag-irregularity surfaces --
  a base conic/asphere plus either a Zernike polynomial correction or
  a grid of measured sag values, interpolatable via
  `DataInterpolations`, already a project dependency) are
  unimplemented -- explicitly scoped as stretch/optional ("Phase 5")
  in the "Support additional Zemax surface types" work that added
  `TOROIDAL` (`FIXED.md` #18)/`COORDBRK`/`TILTSURF`/`ODDASPHE`/
  `XPOLYNOM`/`PARAXIAL` support (see bug #4's current wording), never
  started. Low frequency in the local Zemax install's sample scan
  (1-7 hits each, out of 469 files).
- **30.** `src/zemax.jl`: Zemax surface `TYPE`s `PARAX_XY`, `FRESNELS`,
  `CHEBYSHV` turned up in the same sample-file `TYPE` frequency scan
  that drove the "Support additional Zemax surface types" work (1-2
  hits each, out of 469 sample files under the local Zemax install)
  but were never actually scoped in or out by that work -- a gap in
  the planning itself, not a deliberate exclusion like the items in
  #31 below. Nobody has looked at what these would actually need yet.
- **31.** `src/zemax.jl`: Zemax surface `TYPE`s `NONSEQCO` (non-
  sequential component -- this package's tracer is sequential-only end
  to end, a fundamentally different tracing paradigm, not just a new
  surface type), `USERSURF` (arbitrary user-compiled DLL -- no generic
  conversion possible), `GRINSUR1`/`GRINSUR8`/`GRINSUR9`/`GRINSU11`
  (gradient-index -- needs non-straight-line propagation inside a
  volume, a different tracing loop entirely), `DGRATING`/`BINARY_1`/
  `BINARY_2`/`BINARY_3` (diffraction gratings/binary-kinoform optics --
  need wavelength/diffraction-order-dependent bending, a new
  `AbstractBendType` family), `BIRE__IN`/`BIRE_OUT` (birefringent --
  needs polarization-dependent refractive index), `HOLOGRM2`: each
  needs real new physics beyond "new profile + existing Snell's-law
  bend," deliberately left out of the "Support additional Zemax
  surface types" work as a separate, much larger project per type.
  Tracked here (rather than only in that work's own, now-superseded
  planning document) so they aren't silently forgotten.
- **32.** `test/zemax.jl`: `readZemax`'s own line-by-line parsing of
  `TYPE TOROIDAL`/`COORDBRK`/`TILTSURF`/`ODDASPHE`/`XPOLYNOM`/
  `PARAXIAL` from a real `.zmx` file isn't exercised by any checked-in
  fixture -- coverage for these six types is via
  `zemaxsurfToProfileAperture`/`zemaxsurfToSurface` unit tests built
  from directly-constructed `ZemaxSurf` literals, which never go
  through `readZemax`'s actual text parser at all. (Real-sample
  round-trips through `readZemaxSystem` were done ad hoc against the
  local Zemax install during development, not captured as permanent
  regression tests.) `test/fixtures/test_singlet.zmx` only covers
  `STANDARD`/`EVENASPH`. Extending it (or adding new fixtures) per
  type would close this gap.
- **33.** `src/extended_geo.jl` (`rotationX`/`rotationY`/`rotationZ`,
  used by `zemax.jl`'s `zemaxCoordBreakFrame` for `COORDBRK`/`TILTSURF`
  tilts): hand-rolled 3×3 rotation matrices via `SMatrix`/trig rather
  than a `Rotations.jl` dependency, matching this file's existing
  `findPerpenMap` style -- flagged at the time as "a call worth
  revisiting, not a blocker," never revisited. Working and tested
  (empirically verified against real fold-mirror/toroid samples from
  the local Zemax install), so low priority.
- **34.** `src/tracing.jl`, `modFunc(::ParaxialLensT)`: two documented,
  deliberate simplifications from reusing `surfNormal`'s return slot
  for the ray's transverse offset instead of a true surface normal
  (see `ParaxialProfile`'s docstring, `lens_definitions.jl`): (1) the
  ideal-lens bend is the small-angle-consistent formula only, exact
  for near-axial rays and an approximation farther off-axis -- there's
  no way to decompose a ray's direction into axial/transverse
  components for an exact-at-any-angle formula without the true
  normal; (2) `nIn` always reports `bend.refIndexOut`, never
  `refIndexIn`, since determining which physical side a ray hit from
  (the way `DielectricT`/`MirrorR` do, via `cosI`'s sign) also needs
  the true normal. Real `PARAXIAL` samples have `refIndexIn ==
  refIndexOut` in practice, so (2) is currently harmless; (1) would
  matter for a system doing wide-angle real ray tracing (not just
  first-order layout) through a `PARAXIAL` surface.
- **35.** `src/zemax.jl`, `readZemax`'s `XDAT` parsing (`TYPE
  XPOLYNOM`'s coefficients): `s.extraData[2]` is parsed but never
  interpreted or stored -- an unidentified Zemax control flag, same
  "don't guess the bit" situation as `FTYP`'s undecoded trailing
  tokens (#26). Doesn't affect `sag`, since `SurfProfileXYPoly` never
  reads it, but its actual meaning (per Zemax's own file-format
  reference, not yet consulted) is unknown.
- **36.** `src/tracing.jl`, The surfaces just added need tests
  that test the actual ray tracing at the system level. This should
  be done by first checking the operation on the sample files and then 
  writing new .zmx files with surface tests.
