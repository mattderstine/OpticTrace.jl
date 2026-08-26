# Fixed

Archive of resolved `TODO.md` entries: what the bug/issue was, and
how/where it was fixed. Entries keep the exact number they had in
`TODO.md`'s `Bugs`/`Code issues` lists -- they are not renumbered to be
sequential in this file, since the number is a stable, permanent
identifier shared between the two files (see `TODO.md`'s header for the
full numbering convention). Entries use a bold `**N.**` prefix rather
than real Markdown ordered-list syntax, same reason as `TODO.md`: a
literal prefix survives any renderer's auto-renumbering, a real ordered
list doesn't.

## Bugs

- **1.** `src/mesh_primitives.jl`:
  `GeometryBasics.normals(s::AbstractSurface, nvertices=60)` called
  `inOrOut(s)`, which only had a method for `OptSurface` -- calling
  `normals` on any other `AbstractSurface` (e.g. `ModelSurface`) threw
  a `MethodError`, despite the generic `s::AbstractSurface` signature
  implying it should work for any surface type. (Minor/non-bug aside:
  `GeometryBasics.normals(s::OptSurface, nvertices=60)` was
  byte-for-byte identical to this generic method -- redundant, not
  incorrect.)

  **Fixed**: added `inOrOut(s::ModelSurface)` (`src/mesh_primitives.jl`,
  right after `inOrOut(s::OptSurface)`), which always returns `1` --
  `ModelSurface` has no `mod` field (no refractive-index in/out
  distinction to make, unlike `OptSurface`), so there's no other
  meaningful convention to pick. With `inOrOut` now covering every
  `AbstractSurface` subtype the package defines, the generic
  `GeometryBasics.normals(s::AbstractSurface, ...)` method works
  correctly for `ModelSurface` too, which also made the dedicated
  `GeometryBasics.normals(s::OptSurface, ...)` override (byte-for-byte
  identical body) fully redundant -- deleted, along with its docstring.
  `test/mesh_primitives.jl`'s `"normals (AbstractSurface, ModelSurface
  -- known bug)"` testset was renamed and its `@test_throws
  MethodError` replaced with real assertions (reusing the same
  `SizeLens`-apertured `ModelSurface` fixture, since `samplePoints` has
  no method for `RoundAperture`/`RectAperture`); the pre-existing
  `"normals (OptSurface)"` testset was left unchanged -- it's not a
  test of the removed method specifically, it's the only regression
  coverage `normals` has for `OptSurface`, and it keeps passing
  unchanged once dispatch falls through to the generic method.

- **2.** `src/plotting.jl`, `perimeterRays`: on any ray that fails to
  trace (`status != 0`), it stored a `NaN` `Ray` and then did a bare
  `return` -- which returned `nothing`, discarding the whole `rays`
  vector (including rays already successfully traced) instead of
  continuing to the next perimeter angle. `plotPerimeterRays`/
  `plotPerimeterRays!` (both variants) failed if fed a `geo` where any
  perimeter ray missed, since they iterate over the `nothing` return
  value.

  **Fixed**: `return` changed to `continue` (`src/plotting.jl`,
  `perimeterRays`), so a miss is now recorded (as a `NaN` `Ray`, still
  logged via `println`) and tracing continues to the next perimeter
  angle. Since a `continue` on a miss no longer advances the success
  counter `i`, `rays` (preallocated to length `points`) can have
  trailing unset slots after the loop -- the final line was changed
  from `rays` to `rays[1:i-1]`, so the function now returns only the
  rays that actually traced successfully (length `<= points`) rather
  than a fixed-length vector with `NaN`/undef placeholders for misses.
  Docstrings for `perimeterRays`/`plotPerimeterRays` (`src/plotting.jl`)
  and the corresponding note in `test/plotting.jl`'s
  `"perimeterRays / plotPerimeterRays"` testset were updated to drop
  the "known bug" language and describe the truncated-length return.

- **3.** `src/plotting.jl`, `perimeterRays` (separate, more severe bug
  than #2 -- found while writing `test/plotting.jl`, phase 12):
  `perimeterRays` never got far enough to reach bug #2. Its `Ray(...)`
  call built the direction argument as a raw
  `[cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)]` literal -- a plain `Vector`,
  not a `Vec3` -- but `Ray` only has a constructor for `(Point{N,T},
  Vec{N,T})`. This threw `MethodError` unconditionally, for every call
  regardless of input, before the per-ray tracing loop ever ran.
  Confirmed by direct testing. `perimeterRays`/`plotPerimeterRays`/
  `plotPerimeterRays!` (both variants) are all exported and all failed
  this way -- none of the four was callable at all.

  **Fixed**: wrapped the base-point offset and direction explicitly as
  `Point3`/`Vec3` in `src/plotting.jl`'s `perimeterRays`
  (`Point3(r[1] + radius*cos(ϕ), r[2] + radius*sin(ϕ), r[3])` and
  `Vec3(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ))`), matching the
  component-wise construction style used elsewhere in `src/` (e.g.
  `ORIGIN = Point3(0.,0.,0.)`). The identical defect in the same
  function's miss-handling branch (`Ray([NaN,NaN,NaN], [NaN,NaN,NaN])`)
  was fixed the same way, to `Ray(Point3(NaN,NaN,NaN),
  Vec3(NaN,NaN,NaN))` -- note that branch's separate `return`-instead-
  of-`continue` logic bug is **not** fixed by this; it was still open
  at the time (fixed later, separately, see #2 above). `test/plotting.jl`'s
  `perimeterRays`/`plotPerimeterRays` testset was updated from four
  `@test_throws MethodError` assertions to real assertions (verified
  against the existing test fixture, whose 8 perimeter rays all hit
  successfully, so bug #2's branch isn't exercised by that test).

- **4.** `src/plotting.jl`, `plotOPD!(scene, h::Float64,
  egeo::ExtendedGeometry; ...)`: its final line returned `θr, opdx,
  opdy`, but `θr` was never assigned anywhere in this method (it loops
  over `x = LinRange(-sizeP, sizeP, points)`, not an angle range) --
  `UndefVarError: θr not defined`, thrown only after all the actual
  tracing/plotting work had completed.

  **Fixed**: dropped `θr` from the return -- the method now returns
  `(opdx, opdy)` (`src/plotting.jl`, end of `plotOPD!(scene,
  h::Float64, egeo::ExtendedGeometry; ...)`); docstring updated to
  match.

- **5.** `src/plotting.jl`, `plotOPD!`/`plotOPD3D!` (separate, more
  severe bug than #4 above -- found while writing `test/plotting.jl`,
  phase 12): both methods never got far enough to reach that bug (or,
  for `plotOPD3D!`, to complete at all). Both contained the identical
  line `dirRef = normalize!(usedgeo[1].base.base .- r)`. `.base.base`
  is a `Point3` (an immutable `StaticArrays`-style point), so `.- r`
  produced another immutable `Point3`; `normalize!` is the in-place
  variant and tried to `setindex!` into it, which `Point3` doesn't
  support. Threw `ErrorException("setindex!(::Point{3, Float64},
  value, ::Int) is not defined. Hint: Use MArray or SizedArray...")`
  unconditionally, for any input, near the start of both methods.

  **Fixed**: `dirRef = normalize(Vec3(usedgeo[1].base.base .- r))` in
  both methods (`src/plotting.jl`) -- non-mutating `normalize`
  (matching the "likely fix"), plus an explicit `Vec3(...)` wrap:
  `Point3 .- Vec3` (broadcast subtraction) turned out to still produce
  a `Point3`, not a `Vec3`, which `normalize` alone would have
  preserved -- `Ray`'s direction argument needs a `Vec3`, so this
  surfaced only once `normalize!`'s crash was out of the way. See #12
  below for a second, separate bug this same fix uncovered a bit
  further into `plotOPD!`. `test/plotting.jl`'s
  `"plotOPD!(egeo) / plotOPD3D!"` testset was updated from two
  `@test_throws ErrorException` assertions to real assertions on
  `plotOPD!`'s `(opdx, opdy)` return and `plotOPD3D!`'s returned
  `scene`, with `funcGeo` in the fixture changed from a no-op
  (`p, wl -> nothing`) to one returning the real geometry array (see
  #10 below for why that was necessary too).

- **6.** `src/surface_manipulation.jl`, `reverseProfile!(profile::T)
  where T<:AbstractSurfProfile` (the generic fallback for profile types
  without their own specific method): assumed every such type has an
  `a` field. That's false for `SurfProfileOAConic` (`curv`/`ϵ`/
  `offset`), which fell through to this method -- calling it threw a
  field-access error. (`NoProfile` used to hit this same fallback too,
  but got its own dedicated `reverseProfile!(profile::NoProfile)` no-op
  method separately -- see #7 below.)

  **Fixed**: added a dedicated `reverseProfile!(profile::SurfProfileOAConic)`
  method (`src/surface_manipulation.jl`, right after the
  `SurfProfileConic` method), so `SurfProfileOAConic` no longer falls
  through to the generic fallback at all. **Incomplete on purpose**: the
  new method only negates `profile.curv` (matching its `SurfProfile{,OA}Conic`
  siblings) and leaves `profile.offset` untouched -- reversing an
  off-axis conic also needs the offset re-expressed in the reversed
  geometry's frame, which whoever added this method didn't attempt; it
  prints a message noting this at call time. Tracked as a new, separate
  Code issue (`TODO.md` #22) rather than left silently wrong. (The first
  attempt at this fix mistakenly typed the new method
  `profile::SurfProfileConic` -- identical to the existing method's
  signature just above it -- which Julia silently accepted as a
  redefinition/overwrite of that method instead of a new one, so
  `SurfProfileOAConic` still fell through to the generic fallback
  unchanged; caught by re-running `using OpticTrace` and noticing the
  "Method definition ... overwritten" warning, fixed by correcting the
  type annotation.) `test/surface_manipulation.jl`'s `"SurfProfileOAConic"`
  testset was updated from `@test_throws FieldError` to real assertions
  (`curv` negated, `ϵ` and `offset` both unchanged -- `offset` being
  unchanged is the known-incomplete behavior, not a regression check).

- **7.** `src/surface_manipulation.jl`, `reverseGeo`/`reverseSurface!`:
  despite `reverseGeo`'s `Vector{T} where T<:AbstractSurface` signature,
  `reverseSurface!` only had a method for `OptSurface` -- a `geo`
  containing a `ModelSurface` (e.g. built by
  `roundAperture`/`rectAperture`, both plausible in a real lens system
  with an aperture stop) threw a `MethodError`.

  **Fixed**: added `reverseSurface!(surf::ModelSurface, bpoint,
  epoint)` (`src/surface_manipulation.jl`, right after the `OptSurface`
  method), mirroring that method's structure (`reverseBase!` +
  `reverseProfile!` + recompute coordinate transforms) but skipping the
  `reverseMod!` call -- `ModelSurface` has no `.mod`/coating at all
  (just a fixed `.refIndex` scalar for OPD bookkeeping, no in/out pair
  to swap), so `.refIndex` is left untouched. This depended on #6's
  partial fix (above): `roundAperture`-built `ModelSurface`s use
  `NoProfile`, so exercising this new method calls `reverseProfile!` on
  a `NoProfile`, which needed its own no-op method to not throw. At the
  time, #6 stayed open in `TODO.md` since its other case,
  `SurfProfileOAConic`, was a separate, still-unfixed defect through the
  same generic fallback -- fixed later, separately (see #6 above).
  `test/surface_manipulation.jl`'s `"geo containing a ModelSurface
  (known bug, see TODO.md)"` testset was updated from `@test_throws
  MethodError` to real assertions (order swap, position reflection for
  both surfaces; curvature/index-swap checks for the `OptSurface`;
  untouched `refIndex` for the `ModelSurface`), mirroring the existing
  `"OptSurface-only geo (working)"` test's assertion style.

  Separately, TODO.md flagged an *open question*, not a confirmed bug,
  that this fix doesn't resolve: `reverseSurface!` (both the `OptSurface`
  method and this new `ModelSurface` one) leaves `surf.base.dir`/
  `surf.base.ydir` unchanged when reversing a surface -- only position
  and profile/index data are flipped. Whether propagation-direction
  reversal should also flip the surface's local orientation isn't
  obvious from the code alone; the two methods are at least now
  consistent with each other on this point.

- **8.** `src/characterization.jl`: `sizeOpticSurface` was exported
  (line 6) but never defined anywhere in the repo -- any call to it
  threw immediately (`UndefVarError`/no matching method). Found during
  the initial codebase survey for the phased test-writing plan.

  **Fixed**: dropped `sizeOpticSurface` from the `export` line
  (`src/characterization.jl:6`). Nothing else in `src/`/`test/`
  referenced it, so no other changes were needed.

- **9.** `src/tracing.jl`, `traceSurf`/`traceSurf!` (all four methods:
  `OptSurface` and `ModelSurface`, mutating and non-mutating): each
  checked `if delta == NaN` to detect a ray that missed the surface
  (`deltaToSurf` returns `NaN` on a miss). In IEEE 754, `NaN == NaN` is
  always `false`, so this branch was dead code -- status 1 ("missed")
  was unreachable through this check in all four methods. A ray that
  should miss instead silently propagated `NaN` through the rest of
  the surface-normal/`modFunc`/`clipAperture` math and ended up
  reporting a false "success" (status 0) with a garbage `NaN` trace,
  rather than erroring or correctly reporting a miss.

  **Fixed**: changed all four occurrences to `if isnan(delta)`
  (`src/tracing.jl`, `traceSurf(r, s::OptSurface)`,
  `traceSurf!(trc, r, s::OptSurface)`, `traceSurf(r, s::ModelSurface)`,
  `traceSurf!(trc, r, s::ModelSurface)`), and cleaned up the "known
  bug" language from each method's docstring. `test/trace_geometry.jl`'s
  two miss-case testsets (`OptSurface`, `ModelSurface`) had their
  `@test_broken status == 1` flipped to `@test status == 1` and were
  renamed (dropped "known bug -- see TODO.md").

- **10.** `src/extended_geo.jl`, `updateEGeo!`: named and documented as
  mutating its `ExtendedGeometry` argument ("calls the function
  necessary to create a static geometry for tracing", per its own
  docstring), but the implementation only called `defaultSetupGeo(...)`
  and returned the result -- it never assigned back into `egeo.geo`.

  **Fixed**: `egeo.geo = defaultSetupGeo(...)`
  (`src/extended_geo.jl:156`). `test/foundations.jl`'s
  `"defaultSetupGeo / updateEGeo!"` testset's `updateEGeo!` call had to
  move off its original fixture -- that fixture's `funcGeo` (`testfunc`)
  returns a plain `Tuple`, fine for exercising `defaultSetupGeo`
  directly (kept, `result1`/`result2`), but not assignable into
  `egeo.geo::Array{AbstractSurface}` -- replaced with a second
  `ExtendedGeometry`/`funcGeo` pair whose `funcGeo` returns a real
  `Array{AbstractSurface}`, and dropped the stale `@test_broken
  egeo.geo == ret`. This fix also meant every other `egeo` fixture
  built for testing (e.g. `plotOPD!`/`plotOPD3D!`'s, see #5 above) now
  needed a real `funcGeo` too, since `updateEGeo!` is called at the top
  of both and its result actually replaces `egeo.geo` now.

- **11.** `src/characterization.jl`, `traceMonteCarloRays`: its very
  first executable line built `badray = Ray((NaN, NaN, NaN), (NaN, NaN,
  NaN))` from raw tuples, but `Ray` requires an actual `Point{N,T}`/
  `Vec{N,T}` pair -- this threw a `MethodError` unconditionally, before
  the function ever reached its own arguments or the ray-tracing loop.
  No call to this function could succeed regardless of inputs.

  **Fixed**: `badray = Ray(Point3(NaN, NaN, NaN), Vec3(NaN, NaN, NaN))`
  (`src/characterization.jl:196`), matching the same fix already
  applied to `perimeterRays`'s identical bug (#3 above). Also cleaned
  up the "known bug" language from the function's docstring.
  `test/characterization.jl`'s `traceMonteCarloRays` test was updated
  from `@test_throws MethodError` to real assertions on the returned
  `(rays, cnt, missed)` tuple, reusing the exact same fixture already
  there (a flat `referencePlane` with `semiDiameter=10.0`, 50 rays
  launched within a `radius=5.0` disk straight down `+z`) -- verified
  all 50 rays hit with no clipping (`cnt == 0`, `missed` all zero, no
  `NaN`s in any returned ray).

- **12.** `src/plotting.jl`, `plotOPD!(scene, h::Float64,
  egeo::ExtendedGeometry; ...)`: not previously in `TODO.md` -- found
  while verifying the #5 fix above (this method has no call site
  anywhere in `src/`/`test/`, so it had never actually been run before).
  Further into the same method, the `x`/`y` sweep loop built each
  sample ray's base/direction as raw `SVector` literals
  (`bx = SVector(t, 0., 0.)`, `px = normalize(SVector(t, 0., z).-r)`,
  and the `by`/`py` equivalents) instead of `Point3`/`Vec3` -- same
  defect class as #3/#11 above, and equally fatal: `Ray(::SVector,
  ::SVector)` has no method, so this threw `MethodError` unconditionally
  for every call that got this far, immediately after the #5 fix
  cleared the `normalize!` crash earlier in the same method.

  **Fixed**: `bx = Point3(t, 0., 0.)`, `by = Point3(0., t, 0.)`,
  `px = normalize(Vec3(t, 0., z).-r)`, `py = normalize(Vec3(0., t,
  z).-r)` (`src/plotting.jl`, inside `plotOPD!`'s sample loop). Covered
  by the same updated `test/plotting.jl` testset as #5 above (no
  separate test needed).

## Code issues

- **5.** `src/zemax.jl` / `docs/zemax_reference.md`: reading `.zar`
  archives (Zemax file bundles) wasn't implemented at all -- only the
  Python reference implementation existed, kept in
  `docs/zemax_reference.md` as a starting point for a future
  `readZemaxArchive`-style function.

  **Fixed**: ported the reference implementation to `src/zemax.jl` as
  `ZarEntry`, `lzwDecompress`, `readZemaxArchive`, `listZemaxArchive`,
  and `extractZemaxArchive` (selected-entries and extract-all methods).
  Validated byte-for-byte against a from-scratch Python re-port of the
  same reference, run against two real sample `.zar` files covering
  both header layouts ("earlier"/`0xEA` and "latest"/`0xEC`). Covered
  by `test/zemax.jl`'s `lzwDecompress` unit test (a real, hand-verified
  compressed/decompressed byte pair) and a `HAS_ZAR_SAMPLE`-gated
  integration testset (see `test/helper.jl`) exercising listing and
  both extraction forms against a real archive.

- **23.** `src/zemax.jl`: reading `.zmf` lens-catalog files (how
  vendors like Edmund/Thorlabs publish stock lens families) wasn't
  implemented, and wasn't previously tracked anywhere in this repo --
  found and fixed in the same pass as #5 above. The format is
  undocumented by Ansys/Zemax; `docs/zemax_reference.md` now keeps the
  community reverse-engineering reference this was ported from
  (`rayopt`'s `zemax.py`), including the record layout and the
  obfuscation formula.

  **Fixed**: added `ZmfEntry`, `zmfDeobfuscate`, `readZmfCatalog`,
  `listZmfCatalog`, and `extractZmfCatalog` (selected-lenses and
  extract-all methods) to `src/zemax.jl`. `extractZmfCatalog` writes
  each decoded lens as a `<name>.zmx` file, directly readable by the
  existing `readZemax` -- no changes to that parser were needed.
  Validated byte-for-byte against a from-scratch Python re-port of
  `rayopt`'s `zmf_read`/`zmf_obfuscate` (re-implemented without
  `rayopt`'s ORM/session machinery and without its now-removed
  `numpy.fromstring`/`.tostring()` calls), run against three real
  vendor catalogs. Covered by `test/zemax.jl`'s `zmfDeobfuscate` unit
  test (a real, hand-verified obfuscated/plaintext byte pair) and a
  `HAS_ZMF_SAMPLE`-gated integration testset (see `test/helper.jl`)
  exercising listing and both extraction forms against a real catalog.
