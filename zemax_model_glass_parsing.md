# Handling MODEL Glass in a Zemax Lens Parser

> **Note (2026-09-13): checked against real files.** This doc was written
> before looking at real `.zmx` files; it's kept as-is because the local
> sample corpus is too small to rule out what it describes. In the 469
> local OpticStudio sample files (861 `GLAS` lines, `VERS` 2011–2019):
>
> - **No `SOLVETYPE` lines at all.** The solve is always in the `GLAS`
>   line itself, token 3: `GLAS <name> <solve> <pickupSurf> <Nd> <Vd>
>   <dPgF> ...` with solve `0` fixed, `1` model glass, `2` pickup, `3`
>   substitute. `readZemax` warns if it ever sees a `SOLVETYPE` line, so
>   a real example can be added as a fixture.
> - **`___BLANK` usually means a pickup, not a model glass.** 6 of the 9
>   `___BLANK` surfaces pick up from a surface with no `GLAS` line (air),
>   carrying `1.5 40` filler. The "parse trailing `GLAS` fields as Nd/Vd"
>   fallback below would turn those air gaps into glass, so the
>   implementation follows pickups instead and only reads Nd/Vd/dPgF for
>   solve 1.
>
> What was implemented: `ModelGlass`/`riFormula5`
> (`src/lens_refractive_index.jl`) and `ZEMAX_GLAS_TOKENS`/
> `ZemaxGlassSolve`/`resolveZemaxGlassPickups!`/`zemaxGlassIndexFunc`
> (`src/zemax.jl`). Open questions are tracked in `TODO.md` Code issues
> #39.

## Problem

Zemax `.ZMX` lens files reference glass on each surface's `GLAS` line. Most
surfaces name a real catalog glass (e.g. `N-BK7`), which is resolved against
an `.AGF` catalog. Some surfaces instead carry a placeholder name such as
`___BLANK` — this indicates the material is defined by a **solve** rather
than a catalog lookup, most commonly a **Model Glass** solve (idealized
dispersion via `Nd`, `Vd`, `dPgF`).

A parser that only does catalog lookup will either crash or silently
misassign a bogus/default glass whenever it hits one of these surfaces.
This doc outlines the approach for detecting and correctly handling them.

## Key facts to design around

- `GLAS` is a multi-field line: name, then several numeric fields.
- Catalog lookup only uses the **first field** (the name). The trailing
  numeric fields are meaningless leftovers *unless* a solve is active on
  that surface — in which case they (or an adjacent `SOLVETYPE` line) carry
  the real parameters.
- Model Glass parameters are always the same triple: `Nd` (index at
  587.5618 nm), `Vd` (Abbe number), `dPgF` (deviation from the normal
  line/partial dispersion).
- Field positions/count on `GLAS` have shifted across OpticStudio versions.
  The `.ZMX` text format itself is **unofficial and unsupported** —
  community-reverse-engineered, no guaranteed field layout across versions.
- Model Glass solves can also appear as a **separate `SOLVETYPE` line**
  rather than (or in addition to) trailing `GLAS` fields, especially in
  non-sequential mode (`SOLVETYPE surf, NSC_MATM_#, Nd, Vd, dPgF`) vs.
  sequential mode (`SOLVETYPE surf, GM, Nd, Vd, dPgF`).

## Proposed parsing approach

### 1. Treat glass resolution as a small state machine per surface, not a single lookup

For each surface, resolve material in this order:

1. Read the `GLAS` line's name field.
2. If the name matches an entry in the loaded `.AGF` catalog → resolve
   normally (ignore trailing numeric fields entirely).
3. If the name is a blank/placeholder token (`___BLANK`, `___BLNK`, or
   empty) → mark the surface as **solve-defined** and continue scanning
   forward for the governing solve.
4. Look for a `SOLVETYPE` line associated with the same surface index.
   - If found and its code is the Model Glass type (`GM` sequential /
     `NSC_MATM_#` non-sequential) → extract `Nd`, `Vd`, `dPgF` from the
     `SOLVETYPE` arguments.
   - If no `SOLVETYPE` line is found, fall back to parsing the trailing
     `GLAS` fields positionally as `Nd`, `Vd`, `dPgF`.
5. If neither source yields usable numbers, flag the surface as
   **unresolved glass** rather than defaulting silently (e.g. to air or
   index 1.5). Surface this in parser output/logs.

### 2. Don't hard-code column offsets — validate them

Because field layout varies by Zemax version:

- Parse `GLAS` fields by splitting on whitespace and taking them
  positionally, but keep the version/format assumption explicit and
  documented in one place (a single constant/config, not scattered
  magic indices).
- Add a sanity check on parsed values: `Nd` should be roughly in
  `[1.0, 2.5]`, `Vd` roughly in `[10, 100]`. If out of range, treat as a
  parse failure rather than propagating garbage into the ray trace.
- Keep a small set of reference `.ZMX` files (from known OpticStudio
  versions) as parser test fixtures, specifically including surfaces with
  Model Glass solves, so any future field-layout drift is caught by tests
  rather than by silently wrong refraction results.

### 3. Represent the result explicitly in the parser's data model

Rather than forcing every surface into "named catalog glass," give the
internal surface/material representation a variant type, e.g.:

- `CatalogGlass(name)` — resolved via AGF lookup at trace time.
- `ModelGlass(Nd, Vd, dPgF)` — index computed directly from the model
  formula, no catalog needed.
- `UnresolvedGlass(raw_line)` — parsing failed; carries the raw text for
  debugging, and should cause a loud warning/error rather than a silent
  default.

This keeps downstream ray-tracing code from needing to know *how* the
index was derived — it just asks the material object for `n(wavelength)`.

### 4. Implement the Model Glass dispersion formula

Once `Nd`, `Vd`, `dPgF` are extracted, compute wavelength-dependent index
using OpticStudio's standard idealized dispersion model (Conrady-style
formula built from `Nd`, `Vd`, and the partial dispersion deviation term).
This is the same math used by existing open-source implementations (e.g.
`OpticSim.jl`'s `modelglass` function) and can be validated against known
Nd/Vd/dPgF → index-at-wavelength test cases from those tools.

## Summary of the change

| Before | After |
|---|---|
| `GLAS` name → AGF lookup only | `GLAS` name → AGF lookup, or solve-defined material |
| Blank/placeholder name → crash or wrong default | Detected, routed to Model Glass (or flagged unresolved) |
| Trailing `GLAS` numeric fields ignored/misused | Used only when no catalog match and no `SOLVETYPE`, with range validation |
| Implicit assumption of fixed field layout | Explicit, version-documented, test-covered field layout |

## Open items to confirm against real files

- Exact `SOLVETYPE` code strings for Model Glass in the target OpticStudio
  version(s) you need to support (sequential vs. non-sequential mode).
- Whether trailing `GLAS` fields are ever populated *and* correct even
  when a `SOLVETYPE` line is also present (i.e. which source wins if both
  exist) — verify against a real exported file rather than assuming.
