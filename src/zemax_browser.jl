#=
    a Bonito.jl web UI for browsing a directory of .zmx/.zar/.zmf Zemax
    files and viewing/extracting their contents.

    This file qualifies every Bonito symbol (Bonito.App, Bonito.DOM,
    Bonito.Button, ...) throughout rather than `using Bonito`, relying on
    the `import Bonito` in src/OpticTrace.jl (this file is `include`'d into
    that module, so the binding is already in scope here). GLMakie and
    Bonito both export several identical names bound to unrelated types
    (Button, Slider, Checkbox, Dropdown); since src/OpticTrace.jl already
    does `using GLMakie`, a plain `using Bonito` anywhere in this module
    would make those names ambiguous/undefined module-wide, breaking
    src/plotting.jl's existing (unqualified) use of GLMakie's Button in
    multipleFigures. Confirmed by testing `using GLMakie, Bonito` together:
    a bare `Button` throws UndefVarError.
=#

export zemaxBrowser

"""
    _TABLE_CSS

Bonito's bundled `table.css` asset (plain `table`/`tr`/`th`/`td` element
selectors -- zebra striping, header shading, consistent cell padding; no
class names required), wrapped via `@path` so its bytes are embedded and
survive bundle relocation, same as Bonito's own internal uses of its
bundled CSS assets (e.g. `KaTeXCSS`, `ChoicesCSS`). Injected once as a
root child of [`zemaxBrowserApp`](@ref) so it applies to every `<table>`
this UI renders -- currently just [`renderZemaxFileSummary`](@ref)'s
per-surface table, the only `<table>` this package's Bonito UIs build.
"""
const _TABLE_CSS = Bonito.Asset(Bonito.@path(Bonito.dependency_path("table.css")))

"""
    _TABLE_DARK_MODE_CSS

Dark-mode override for [`_TABLE_CSS`](@ref)'s zebra striping. table.css's
`tr:nth-child(even)`/`tr:hover`/`th` rules (`#f2f2f2`/`#ddd`/`#e2e2e2`
background, `white` header text -- light grays close to white) have no
media query of their own, so they stay active in *every* color scheme,
dark included, unless something more specific overrides them per rule --
overriding only some of table.css's rules under `@media
(prefers-color-scheme: dark)` leaves the untouched ones still in force
(a first pass here left `nth-child(even)` unoverridden on the assumption
that meant "no styling, falls back to the page background" -- it doesn't;
it meant "table.css's own unconditional `#f2f2f2` still applies", a near-
white stripe against a dark page). Every table.css rule that sets a
color needs an explicit dark-mode counterpart here, not just the ones
that looked obviously wrong:
- `nth-child(even)`: `background-color: transparent`, canceling
  table.css's `#f2f2f2` outright so these rows actually show the page's
  own dark background (never given its own color here, since it's the
  browser/OS-dependent ambient one this file has no fixed value for).
- `nth-child(odd)`: `#4a4a4a`, a clearly-lighter-than-ambient gray so the
  zebra stripe reads against any reasonably dark ambient background.
- `tr:hover`: `#606060`, lighter again than `#4a4a4a` so hovering stays
  visually distinct from the zebra stripe itself (not just from the
  transparent rows).
- `th`: background `#1d1d1d` (kept from the first pass) plus `color:
  #999999` (medium gray, replacing table.css's `white` -- requested
  directly, distinct from the zebra-contrast fixes above).

`tr:hover`'s background carries `!important`, unlike the other three
rules here. A `table tbody tr:hover` specificity-bump selector (relying on
the `<tbody>` browsers implicitly insert around bare `<tr>` children) was
tried first to outrank `table tr:nth-child(odd)`/`table
tr:nth-child(even)` without `!important` -- on paper its specificity
(3 elements, 1 pseudo-class) already beats theirs (2 elements, 1
pseudo-class) -- but it still didn't visibly win in testing (light-mode
hover, via table.css's own unconditional `tr:hover`, worked throughout;
only the dark-mode override silently lost). Rather than keep chasing
which part of the cascade was actually deciding it, `!important` sidesteps
the question: it outranks every non-`!important` declaration regardless
of selector specificity or source order, which is exactly the guarantee
needed here against a stylesheet (table.css) whose position in the
cascade this file doesn't control.
"""
const _TABLE_DARK_MODE_CSS = Bonito.Styles(Bonito.CSS(
    "@media (prefers-color-scheme: dark)",
    Bonito.CSS("table tr:nth-child(even)", "background-color" => "transparent"),
    Bonito.CSS("table tr:nth-child(odd)", "background-color" => "#4a4a4a"),
    Bonito.CSS("table tr:hover", "background-color" => "#606060 !important"),
    Bonito.CSS("table th", "background-color" => "#1d1d1d", "color" => "#999999"),
))

"""
    zemaxFileKind(path::AbstractString) -> Union{Symbol,Nothing}

Classify `path` by extension (case-insensitive): `:zmx`, `:zar`, `:zmf`,
or `nothing` for anything else. Used to dispatch content-pane rendering
in the browser UI.
"""
function zemaxFileKind(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    ext == ".zmx" && return :zmx
    ext == ".zar" && return :zar
    ext == ".zmf" && return :zmf
    return nothing
end

"""
    ZemaxFileSummary

Display-oriented bundle of a `.zmx` file's parsed fields, as returned by
[`zemaxFileSummary`](@ref). Wraps [`readZemax`](@ref)'s `(zsurfs,
header)` return value's name/units/wavelengths fields verbatim -- no
new parsing.

Fields:
    name::String                  - system name, from readZemax
    units::String                  - length units, from readZemax
    wavelengths::Vector{Float64}   - fixed 24-element wavelength array, from readZemax
    zsurfs::Vector{ZemaxSurf}      - parsed surfaces, from readZemax
    geo::Union{Vector{AbstractSurface},Nothing}  - built geometry from
        readZemaxSystem, for the geometry preview plot; nothing if it
        couldn't be built (see geoError)
    geoError::Union{String,Nothing}  - error message when geo is
        nothing (e.g. an unsupported Zemax surface type or a MODE NSC
        file, see readZemaxSystem); nothing when geo was built
        successfully. Exactly one of geo/geoError is nothing.
"""
struct ZemaxFileSummary
    name::String
    units::String
    wavelengths::Vector{Float64}
    zsurfs::Vector{ZemaxSurf}
    geo::Union{Vector{AbstractSurface},Nothing}
    geoError::Union{String,Nothing}
end

"""
    zemaxFileSummary(path::String; glassCatalog = defaultGlassCatalog) -> ZemaxFileSummary

Read and summarize a `.zmx` file at `path` via [`readZemax`](@ref). This is
the single code path the browser UI uses to render a `.zmx` file's fields,
whether `path` is a file found directly in the directory tree or a
temporary file extracted from inside an archive (see
[`zemaxEntitySummary`](@ref)). Also builds the geometry (via
[`readZemaxSystem`](@ref), mirroring [`viewZemaxFile`](@ref)'s own
"parse, then build geometry" sequence -- `glassCatalog` is forwarded to
it unchanged, same keyword/default as `readZemaxSystem` itself) for the
geometry preview plot; a file using a still-unsupported Zemax surface
type, a material missing from `glassCatalog`, or a `MODE NSC` file
fails only at that stage (`readZemax` itself never throws on those), so
the failure is caught and stored in the returned
[`ZemaxFileSummary`](@ref)'s `geoError` rather than propagated -- the
surface table must still render even when the plot can't be built. The
browser UI itself always uses the default `glassCatalog` (no catalog
picker); the keyword exists so tests can exercise the geometry-building
success path against fixtures with fictional material names not in the
real, machine-local catalog (see `test/zemax.jl`'s own `testCatalog`
pattern).
"""
function zemaxFileSummary(path::String; glassCatalog = defaultGlassCatalog)::ZemaxFileSummary
    zsurfs, header = readZemax(path)
    geo, geoError = try
        readZemaxSystem(path; glassCatalog).geo, nothing
    catch e
        nothing, sprint(showerror, e)
    end
    return ZemaxFileSummary(header.name, header.units, header.wavelengths, zsurfs, geo, geoError)
end

"""
    ZemaxArchiveEntity

Display-oriented metadata for one entity inside a `.zar`/`.zmf` archive,
as returned by [`listZemaxArchiveEntities`](@ref).

Fields:
    name::String        - entity name, exactly as listZemaxArchive/listZmfCatalog report it
    byteSize::Int         - length of the payload extraction would write
    previewable::Bool     - true iff a .zmx detail view is possible for this entity
"""
struct ZemaxArchiveEntity
    name::String
    byteSize::Int
    previewable::Bool
end

isZmxName(name::AbstractString) = endswith(lowercase(name), ".zmx")

"""
    listZemaxArchiveEntities(archivePath::String) -> Vector{ZemaxArchiveEntity}

List the entities inside a `.zar` or `.zmf` archive (dispatched on
`archivePath`'s extension), for display in the browser UI. `.zar` entries
are previewable iff their name ends in `.zmx`; every `.zmf` entry is
previewable, since a `.zmf` lens record is always already-decoded `.zmx`
text (see [`ZmfEntry`](@ref)). Throws an error if `archivePath` isn't a
`.zar` or `.zmf` file.
"""
function listZemaxArchiveEntities(archivePath::String)::Vector{ZemaxArchiveEntity}
    ext = lowercase(splitext(archivePath)[2])
    if ext == ".zar"
        return [ZemaxArchiveEntity(e.name, length(e.data), isZmxName(e.name))
                for e in readZemaxArchive(archivePath)]
    elseif ext == ".zmf"
        return [ZemaxArchiveEntity(e.name, length(e.data), true)
                for e in readZmfCatalog(archivePath)]
    else
        error("listZemaxArchiveEntities: not a .zar or .zmf archive: $archivePath")
    end
end

"""
    zemaxEntitySummary(archivePath::String, entityName::String) -> ZemaxFileSummary

Preview a `.zmx` entity found inside a `.zar`/`.zmf` archive. Since
`ZarEntry`/`ZmfEntry` payloads are in-memory bytes rather than files on
disk, this reuses [`extractZemaxArchive`](@ref)/[`extractZmfCatalog`](@ref)
to write the single named entity to a fresh temporary directory, then
calls [`zemaxFileSummary`](@ref) on the extracted path -- no separate
byte-writing logic. Only call this for an entity whose
[`ZemaxArchiveEntity`](@ref)`.previewable` is `true`; callers must check
that flag themselves first.
"""
function zemaxEntitySummary(archivePath::String, entityName::String)::ZemaxFileSummary
    ext = lowercase(splitext(archivePath)[2])
    tmp = mktempdir()
    paths = if ext == ".zar"
        extractZemaxArchive(archivePath, [entityName]; outputPath = tmp)
    elseif ext == ".zmf"
        extractZmfCatalog(archivePath, [entityName]; outputPath = tmp)
    else
        error("zemaxEntitySummary: not a .zar or .zmf archive: $archivePath")
    end
    return zemaxFileSummary(only(paths))
end

"""
    zemaxEntityBytes(archivePath::String, entityName::String) -> Vector{UInt8}

Return the raw payload bytes for `entityName` inside a `.zar`/`.zmf`
archive at `archivePath`, without writing anything to disk -- unlike
[`zemaxEntitySummary`](@ref), which is specifically for `.zmx` entities
routed through `readZemax`. Used to preview non-`.zmx` entities (e.g. a
`.zar` archive's `.AGF` glass-catalog entries) as plain text; in
practice `.zmf` entries never need this since every `.zmf` entity is
already `previewable` (see [`listZemaxArchiveEntities`](@ref)). Throws
an error if `archivePath` isn't a `.zar`/`.zmf` archive, or if
`entityName` isn't present in it.
"""
function zemaxEntityBytes(archivePath::String, entityName::String)::Vector{UInt8}
    ext = lowercase(splitext(archivePath)[2])
    entries = if ext == ".zar"
        readZemaxArchive(archivePath)
    elseif ext == ".zmf"
        readZmfCatalog(archivePath)
    else
        error("zemaxEntityBytes: not a .zar or .zmf archive: $archivePath")
    end
    idx = findfirst(e -> e.name == entityName, entries)
    idx === nothing && error("zemaxEntityBytes: entity \"$entityName\" not found in $archivePath")
    return entries[idx].data
end

"""
    looksLikeText(bytes::AbstractVector{UInt8}) -> Bool

Heuristic check for whether `bytes` is plausibly text rather than
binary: empty counts as text; otherwise, a NUL byte anywhere in the
first 512 bytes (binary formats almost always have one early) rules it
out, and the whole buffer must decode as valid UTF-8. Not a
general-purpose file-type sniffer -- just enough to decide whether a
text preview of an unrecognized archive entity is worth attempting.
"""
function looksLikeText(bytes::AbstractVector{UInt8})
    isempty(bytes) && return true
    any(==(0x00), @view bytes[1:min(end, 512)]) && return false
    return isvalid(String, bytes)
end

"""
    extractZemaxEntities(archivePath::String, names=nothing; outputPath=nothing) -> Vector{String}

Extension-dispatched wrapper over [`extractZemaxArchive`](@ref)/
[`extractZmfCatalog`](@ref): routes to the extract-all overload when
`names === nothing`, otherwise to the named-entries overload. Inherits
their `mkpath`/error-on-missing-name behavior unchanged -- see those
functions' docstrings for the full contract.
"""
function extractZemaxEntities(archivePath::String,
                               names::Union{Nothing,AbstractVector{<:AbstractString}} = nothing;
                               outputPath = nothing)
    ext = lowercase(splitext(archivePath)[2])
    if ext == ".zar"
        return names === nothing ? extractZemaxArchive(archivePath; outputPath = outputPath) :
                                    extractZemaxArchive(archivePath, names; outputPath = outputPath)
    elseif ext == ".zmf"
        return names === nothing ? extractZmfCatalog(archivePath; outputPath = outputPath) :
                                    extractZmfCatalog(archivePath, names; outputPath = outputPath)
    else
        error("extractZemaxEntities: not a .zar or .zmf archive: $archivePath")
    end
end

"""
    defaultExtractionOutputPath(archivePath::String) -> String

The default extraction output directory for `archivePath`: its parent
directory -- the folder [`extractZemaxArchive`](@ref)/
[`extractZmfCatalog`](@ref) will create their own same-named,
extension-stripped subfolder inside of, not that subfolder itself
(showing the subfolder here would double-nest once it's created). Used
to pre-fill the browser UI's editable output-path field. Throws an
error unless `archivePath` is a `.zar` or `.zmf` file.
"""
function defaultExtractionOutputPath(archivePath::String)::String
    ext = lowercase(splitext(archivePath)[2])
    (ext == ".zar" || ext == ".zmf") ||
        error("defaultExtractionOutputPath: not a .zar or .zmf archive: $archivePath")
    return dirname(archivePath)
end

#=
    Bonito UI layer below -- every Bonito symbol is qualified (Bonito.App,
    Bonito.DOM, ...), see the module-collision note at the top of this file.
=#

"""
    renderGeometryPreview(summary::ZemaxFileSummary)

Render a live, interactive 3D preview of `summary.geo` via
[`plotGeometry3D`](@ref), sized to fill its container -- see
[`renderZemaxFileSummary`](@ref)'s docstring for the flex-column
ancestor chain (all the way up to [`zemaxBrowserApp`](@ref)'s
`height: 100vh` root) that gives this container a real, non-shrink-
wrapped height for it to fill. The `Figure` is wrapped in
`WGLMakie.WithConfig(fig; resize_to = :parent)` rather than placed bare
-- `WGLMakie` (imported, not `using`'d, in `src/OpticTrace.jl`, for the
same collision-avoidance reason as `Bonito`) defines
`Bonito.jsrender(session, ::Makie.FigureLike)` (and, for `WithConfig`,
an analogous method), so Bonito's own rendering pipeline dispatches to
it automatically, and `resize_to = :parent` makes the rendered WebGL
canvas track this returned `div`'s actual on-screen size (via a JS
`ResizeObserver` on the canvas's grandparent -- confirmed by reading
WGLMakie's `get_resize_element` in its bundled JS -- which is exactly
this `div`) rather than staying at `plotGeometry3D`'s nominal pixel
`size`, so the plot ends up matching the preview column's width and
the remaining vertical space, and keeps tracking it if the window is
resized. This still gives a mouse-manipulable (rotate/pan/zoom, via
`plotGeometry3D`'s `cam3d_cad!` camera) view over the live session, not
a static snapshot. When `summary.geo` is `nothing` (geometry
construction failed -- see [`ZemaxFileSummary`](@ref)'s `geoError`
field), renders that error message instead, with matching `flex`
sizing so the layout doesn't jump between the two cases.
"""
function renderGeometryPreview(summary::ZemaxFileSummary)
    fillStyle = Bonito.Styles("flex" => "1 1 auto", "min-height" => "0", "width" => "100%")
    summary.geo === nothing &&
        return Bonito.DOM.div("Geometry preview unavailable: ", summary.geoError; style = fillStyle)
    fig, _ax = plotGeometry3D(summary.geo; size = (800, 600))
    return Bonito.DOM.div(WGLMakie.WithConfig(fig; resize_to = :parent); style = fillStyle)
end

"""
    renderZemaxFileSummary(summary::ZemaxFileSummary)

Render a [`ZemaxFileSummary`](@ref) as a Bonito DOM node: name, units, the
populated (nonzero) wavelengths, a table with one row per surface
showing the same fields [`printZemaxSurfs`](@ref) already treats as the
interesting ones (type, curvature, distance, material, radius, stop,
conic, coating, comm), and a geometry preview plot below the table (see
[`renderGeometryPreview`](@ref)). Used both for a `.zmx` file found
directly in the tree and for a previewable `.zmx` entity found inside
an archive.

The returned `div` is itself a `flex-direction: column` flex container
with `flex: 1 1 auto`, so it stretches to fill the preview column's
full height (that column, in turn, is a flex column too -- see
[`zemaxBrowserApp`](@ref)/[`archiveContentPane`](@ref)) and lets the
geometry preview -- the one child with its own `flex: 1 1 auto` (see
[`renderGeometryPreview`](@ref)) -- grow to consume whatever's left
below `header`/the surface table, which stay their natural content
size. `min-height: 0` overrides flex's default `min-height: auto`,
which would otherwise stop this `div` (and in turn the plot) from
shrinking below its content's intrinsic height when space is tight.
"""
function renderZemaxFileSummary(summary::ZemaxFileSummary)
    usedWavelengths = filter(!=(0.0), summary.wavelengths)
    header = Bonito.DOM.div(
        Bonito.DOM.h3(summary.name),
        Bonito.DOM.div("Units: ", summary.units),
        Bonito.DOM.div("Wavelengths: ", join(string.(usedWavelengths), ", ")),
    )
    headerRow = Bonito.DOM.tr(
        Bonito.DOM.th("#"), Bonito.DOM.th("type"), Bonito.DOM.th("curvature"),
        Bonito.DOM.th("distance"), Bonito.DOM.th("material"), Bonito.DOM.th("semiDiam"),
        Bonito.DOM.th("stop"), Bonito.DOM.th("conic"), Bonito.DOM.th("coating"), Bonito.DOM.th("comm"),
    )
    rows = [Bonito.DOM.tr(
                Bonito.DOM.td(string(i)), Bonito.DOM.td(s.type), Bonito.DOM.td(string(s.curvature)),
                Bonito.DOM.td(string(s.distance)), Bonito.DOM.td(s.material), Bonito.DOM.td(string(s.radius)),
                Bonito.DOM.td(string(s.stop)), Bonito.DOM.td(string(s.conic)), Bonito.DOM.td(s.coating),
                Bonito.DOM.td(s.comm),
            ) for (i, s) in enumerate(summary.zsurfs)]
    return Bonito.DOM.div(header, Bonito.DOM.table(headerRow, rows...), renderGeometryPreview(summary);
                           style = Bonito.Styles("display" => "flex", "flex-direction" => "column",
                                                  "flex" => "1 1 auto", "min-height" => "0"))
end

"""
    renderNonPreviewableEntity(entity::ZemaxArchiveEntity)

Render a non-previewable archive entity as name + byte size + a "preview
not available" note -- no attempt to interpret its content.
"""
function renderNonPreviewableEntity(entity::ZemaxArchiveEntity)
    return Bonito.DOM.div(
        Bonito.DOM.h3(entity.name),
        Bonito.DOM.div("Size: $(entity.byteSize) bytes"),
        Bonito.DOM.div("Preview not available for this entity type."),
    )
end

"""
    renderTextPreview(name::String, bytes::AbstractVector{UInt8})

Render a plain-text preview of `bytes` (assumed to already look like
text -- see [`looksLikeText`](@ref)): `name` as a header, then the
first 10 lines, line-wrapped so long lines stay on screen instead of
forcing horizontal scroll.
"""
function renderTextPreview(name::String, bytes::AbstractVector{UInt8})
    lines = split(String(bytes), '\n')
    firstLines = join(first(lines, min(10, length(lines))), "\n")
    return Bonito.DOM.div(
        Bonito.DOM.h3(name),
        Bonito.DOM.pre(firstLines;
            style = Bonito.Styles(
                "white-space" => "pre-wrap",
                "word-break" => "break-word",
                "max-width" => "100%",
                "border" => "1px solid var(--optictrace-border, #ccc)",
                "padding" => "6px",
            )),
    )
end

"""
    archiveContentPane(path::String) -> (middle, preview)

Render the two right-hand columns for a `.zar`/`.zmf` archive at `path`.
`middle`: a title, a scrollable (bounded-height, see
[`_SCROLLABLE_LIST_STYLE`](@ref)) list of its entities (each with a
"select" and an "Extract" button; `.zmf` entity names get their implicit
`.zmx` extension appended for display -- `.zar` entries already carry
theirs), an editable output-path field pre-filled via
[`defaultExtractionOutputPath`](@ref) with a "Browse..." button that
toggles a [`filePicker`](@ref) (`:directory` mode, rooted at the output
path if it already exists, else at `dirname(path)`) open beneath it --
picking a directory (or canceling) writes the result into the field (or
leaves it unchanged) and collapses the picker again -- an "Extract All"
button (which always writes into a new subfolder named after the
archive under the output-path field's directory -- see
[`extractZemaxArchive`](@ref)/[`extractZmfCatalog`](@ref) -- unlike the
per-entity "Extract" buttons above, which write directly into it), and
a status line. `preview`: a detail pane driven by which entity is
currently selected -- a `.zmx` entity renders via
[`renderZemaxFileSummary`](@ref); a non-`.zmx` entity that looks like
text (see [`looksLikeText`](@ref)/[`zemaxEntityBytes`](@ref), e.g. a
`.zar` archive's `.AGF` glass-catalog entries) renders via
[`renderTextPreview`](@ref); anything else falls back to
[`renderNonPreviewableEntity`](@ref). Both share one `selectedEntity`
observable (reset on every re-render of this pane, i.e. whenever the
outer `selectedPath` changes to a different archive) so a click in
`middle` updates `preview`.
"""
function archiveContentPane(path::String)
    entities = listZemaxArchiveEntities(path)
    archiveKind = zemaxFileKind(path)
    selectedEntity = Bonito.Observable{Union{Nothing,String}}(nothing)
    status = Bonito.Observable("")
    outputField = Bonito.TextField(defaultExtractionOutputPath(path); style = nothing)

    browsing = Bonito.Observable(false)
    browseButton = Bonito.Button("Browse..."; style = nothing)
    Bonito.on(browseButton.value) do _clicked
        browsing[] = true
    end
    browsePane = Bonito.map!(Bonito.Observable{Any}(Bonito.DOM.div()), browsing) do isBrowsing
        isBrowsing || return Bonito.DOM.div()
        root = isdir(outputField.value[]) ? outputField.value[] : dirname(path)
        picker, pickedDir, active = filePicker(root; mode = :directory)
        Bonito.on(active) do isActive
            isActive && return
            pickedDir[] === nothing || (outputField.value[] = pickedDir[])
            browsing[] = false
        end
        return picker
    end

    entityRows = map(entities) do ent
        displayName = archiveKind === :zmf ? ent.name * ".zmx" : ent.name
        selectBtn = Bonito.Button(ent.previewable ? displayName : "$displayName (no preview)"; style = nothing)
        Bonito.on(selectBtn.value) do _clicked
            selectedEntity[] = ent.name
        end
        extractBtn = Bonito.Button("Extract"; style = nothing)
        Bonito.on(extractBtn.value) do _clicked
            try
                paths = extractZemaxEntities(path, [ent.name]; outputPath = outputField.value[])
                status[] = "Extracted to: " * join(paths, ", ")
            catch e
                status[] = "Error: $(sprint(showerror, e))"
            end
        end
        return Bonito.DOM.div(selectBtn, extractBtn, " ($(ent.byteSize) bytes)")
    end

    extractAllBtn = Bonito.Button("Extract All"; style = nothing)
    Bonito.on(extractAllBtn.value) do _clicked
        try
            paths = extractZemaxEntities(path; outputPath = outputField.value[])
            outDir = isempty(paths) ? outputField.value[] : dirname(paths[1])
            status[] = "Extracted $(length(paths)) file(s) to: $outDir"
        catch e
            status[] = "Error: $(sprint(showerror, e))"
        end
    end

    preview = Bonito.map!(Bonito.Observable{Any}(Bonito.DOM.div()), selectedEntity) do name
        if name === nothing
            return Bonito.DOM.div("Select an entity to preview it.")
        end
        idx = findfirst(e -> e.name == name, entities)
        ent = entities[idx]
        ent.previewable && return renderZemaxFileSummary(zemaxEntitySummary(path, name))
        bytes = zemaxEntityBytes(path, name)
        displayName = archiveKind === :zmf ? ent.name * ".zmx" : ent.name
        return looksLikeText(bytes) ? renderTextPreview(displayName, bytes) :
                                        renderNonPreviewableEntity(ent)
    end

    middle = Bonito.DOM.div(
        Bonito.DOM.h3("Entities in $(basename(path))"),
        Bonito.DOM.div(entityRows...; style = _SCROLLABLE_LIST_STYLE),
        Bonito.DOM.div("Output path: ", outputField, browseButton, extractAllBtn;
                        style = Bonito.Styles("margin-top" => "16px")),
        browsePane,
        Bonito.DOM.div(status),
    )
    return middle, preview
end

"""
    contentPane(path::String) -> (middle, preview)

Top-level content-pane dispatch, by [`zemaxFileKind`](@ref)`(path)`,
returning the middle and preview columns for [`zemaxBrowserApp`](@ref):
a `.zmx` file has no list of its own, so `middle` is empty and the full
summary (via [`renderZemaxFileSummary`](@ref)) goes straight to
`preview`; a `.zar`/`.zmf` archive renders both via
[`archiveContentPane`](@ref); an empty `path` (nothing selected yet) or
an unrecognized extension puts a placeholder message in `middle`,
leaving `preview` blank.
"""
function contentPane(path::String)
    isempty(path) && return Bonito.DOM.div("Select a file to view its contents."), Bonito.DOM.div()
    kind = zemaxFileKind(path)
    kind === :zmx && return Bonito.DOM.div(), renderZemaxFileSummary(zemaxFileSummary(path))
    (kind === :zar || kind === :zmf) && return archiveContentPane(path)
    return Bonito.DOM.div("Unsupported file: $path"), Bonito.DOM.div()
end

"""
    _SERVER_CLOSE_GRACE_PERIOD

Seconds [`zemaxBrowserApp`](@ref) waits, after its last live session
closes, before actually closing the server -- long enough for a page
refresh (which briefly drops and re-establishes the connection) to
reconnect without losing the server out from under it.
"""
const _SERVER_CLOSE_GRACE_PERIOD = 2.0

"""
    zemaxBrowserApp(rootDir::String; serverRef=nothing, closeOnDisconnect::Bool=false) -> Bonito.App

Build the full browser app for `rootDir` as three columns: a fixed-width
"File" column (via [`filePicker`](@ref), `:file` mode, filtered to
`.zmx`/`.zar`/`.zmf`), a fixed-width middle column, and a flexible
preview column -- the latter two from [`contentPane`](@ref), driven by
[`filePicker`](@ref)'s own `selected` observable (refreshed whenever the
user double-clicks a file or clicks "Select", not on every single click)
and, for an archive, `archiveContentPane`'s own inner `selectedEntity`
observable for the preview column alone. The middle/preview pair is
built together in one reactive step (not as two independent ones) so
both share that same inner state.

`serverRef`/`closeOnDisconnect` are [`zemaxBrowser`](@ref)'s plumbing for
closing the server when its browser tab closes -- both are inert by
default (no session tracking at all), so calling this directly (e.g. in
tests, with no real `Server` behind it) is unaffected. When
`closeOnDisconnect` is `true`, a live count of connected sessions is kept
outside the per-session closure (so a page refresh's brief
disconnect-then-reconnect only causes a net no-op, not two independent
counts); once that count drops to zero, `serverRef[]` is closed after
[`_SERVER_CLOSE_GRACE_PERIOD`](@ref) seconds, rechecked at that point in
case a refresh reconnected in the meantime.
"""
function zemaxBrowserApp(rootDir::String; serverRef = nothing, closeOnDisconnect::Bool = false)
    activeSessions = Ref(0)
    sessionsLock = ReentrantLock()
    return Bonito.App() do session
        picker, selectedPath, _active = filePicker(rootDir; mode = :file,
                                                    extensions = [".zmx", ".zar", ".zmf"])
        if closeOnDisconnect && serverRef !== nothing
            lock(sessionsLock) do
                activeSessions[] += 1
            end
            Bonito.on(session.on_close) do _closed
                remaining = lock(sessionsLock) do
                    activeSessions[] -= 1
                end
                remaining > 0 && return
                @async begin
                    sleep(_SERVER_CLOSE_GRACE_PERIOD)
                    stillIdle = lock(sessionsLock) do
                        activeSessions[] <= 0
                    end
                    s = serverRef[]
                    if stillIdle && s !== nothing
                        close(s)
                    end
                end
            end
        end
        fileColumn = Bonito.DOM.div(Bonito.DOM.h3("File"), picker;
                                     style = Bonito.Styles("width" => "300px", "flex" => "0 0 auto",
                                                            "overflow-y" => "auto"))
        detailColumns = Bonito.map!(Bonito.Observable{Any}(Bonito.DOM.div()), selectedPath) do p
            middle, preview = contentPane(p === nothing ? "" : p)
            return Bonito.DOM.div(
                Bonito.DOM.div(middle; style = Bonito.Styles("width" => "360px", "flex" => "0 0 auto",
                                                              "overflow-y" => "auto", "padding" => "0 12px")),
                Bonito.DOM.div(preview; style = Bonito.Styles("flex" => "1 1 auto", "overflow-y" => "auto",
                                                               "padding" => "0 12px", "display" => "flex",
                                                               "flex-direction" => "column"));
                style = Bonito.Styles("display" => "flex", "flex" => "1 1 auto"),
            )
        end
        # height:100vh (rather than the previous shrink-to-content sizing,
        # with align-items:flex-start) gives every column a definite height
        # to stretch to -- needed so the preview column's geometry plot (see
        # renderGeometryPreview) can flex-grow to fill the space below its
        # table down to the bottom of the window, rather than just being as
        # tall as its own content.
        return Bonito.DOM.div(_THEME_STYLES, _TABLE_CSS, _TABLE_DARK_MODE_CSS, fileColumn, detailColumns;
                               style = Bonito.Styles("display" => "flex", "height" => "100vh"))
    end
end

"""
    zemaxBrowser(rootDir::String; port::Integer=8081, openBrowser::Bool=true,
                 closeOnDisconnect::Bool=true) -> Bonito.Server

Launch a standalone Zemax file browser: starts a `Bonito.Server` on
`127.0.0.1:port` serving [`zemaxBrowserApp`](@ref)`(rootDir)`, optionally
opens a browser tab pointing at it (`openBrowser`, via `openInBrowser`
-- defined in `UItools/filepicker.jl` and shared with that file's
[`filePickerDialog`](@ref)), and returns the live `Server` immediately
without blocking. By default (`closeOnDisconnect=true`,
since this opens exactly one tab and is meant to live only as long as
that tab does) the server closes itself once its tab/window closes --
after a short grace period (see [`_SERVER_CLOSE_GRACE_PERIOD`](@ref)) so
a page refresh doesn't kill it. Pass `closeOnDisconnect=false` to keep
today's manual-lifecycle behavior instead, closing it yourself later with
`close(server)`. Throws an error if `rootDir` isn't a directory.
"""
function zemaxBrowser(rootDir::String; port::Integer = 8081, openBrowser::Bool = true,
                       closeOnDisconnect::Bool = true)
    isdir(rootDir) || error("zemaxBrowser: not a directory: $rootDir")
    serverRef = Ref{Union{Nothing,Bonito.Server}}(nothing)
    app = zemaxBrowserApp(rootDir; serverRef, closeOnDisconnect)
    server = Bonito.Server(app, "127.0.0.1", Int(port))
    serverRef[] = server
    openBrowser && openInBrowser(Bonito.online_url(server, "/"))
    return server
end
