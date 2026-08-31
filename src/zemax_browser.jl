#=
    a Bonito.jl web UI for browsing a directory of .zmx/.zar/.zmf Zemax
    files and viewing/extracting their contents.

    This file `import`s Bonito rather than `using` it, and qualifies every
    Bonito symbol (Bonito.App, Bonito.DOM, Bonito.Button, ...) throughout.
    GLMakie and Bonito both export several identical names bound to
    unrelated types (Button, Slider, Checkbox, Dropdown); since
    src/OpticTrace.jl already does `using GLMakie`, a plain `using Bonito`
    anywhere in this module would make those names ambiguous/undefined
    module-wide, breaking src/plotting.jl's existing (unqualified) use of
    GLMakie's Button in multipleFigures. Confirmed by testing `using
    GLMakie, Bonito` together: a bare `Button` throws UndefVarError.
=#
import Bonito

export zemaxBrowser

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
[`zemaxFileSummary`](@ref). Wraps [`readZemax`](@ref)'s `(zsurfs, name,
units, wavelengths)` return tuple verbatim -- no new parsing.

Fields:
    name::String                  - system name, from readZemax
    units::String                  - length units, from readZemax
    wavelengths::Vector{Float64}   - fixed 24-element wavelength array, from readZemax
    zsurfs::Vector{ZemaxSurf}      - parsed surfaces, from readZemax
"""
struct ZemaxFileSummary
    name::String
    units::String
    wavelengths::Vector{Float64}
    zsurfs::Vector{ZemaxSurf}
end

"""
    zemaxFileSummary(path::String) -> ZemaxFileSummary

Read and summarize a `.zmx` file at `path` via [`readZemax`](@ref). This is
the single code path the browser UI uses to render a `.zmx` file's fields,
whether `path` is a file found directly in the directory tree or a
temporary file extracted from inside an archive (see
[`zemaxEntitySummary`](@ref)).
"""
function zemaxFileSummary(path::String)::ZemaxFileSummary
    zsurfs, name, units, wavelengths = readZemax(path)
    return ZemaxFileSummary(name, units, wavelengths, zsurfs)
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
    renderZemaxFileSummary(summary::ZemaxFileSummary)

Render a [`ZemaxFileSummary`](@ref) as a Bonito DOM node: name, units, the
populated (nonzero) wavelengths, and a table with one row per surface
showing the same fields [`printZemaxSurfs`](@ref) already treats as the
interesting ones (type, curvature, distance, material, radius, stop,
conic, coating, comm). Used both for a `.zmx` file found directly in the
tree and for a previewable `.zmx` entity found inside an archive.
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
        Bonito.DOM.th("distance"), Bonito.DOM.th("material"), Bonito.DOM.th("radius"),
        Bonito.DOM.th("stop"), Bonito.DOM.th("conic"), Bonito.DOM.th("coating"), Bonito.DOM.th("comm"),
    )
    rows = [Bonito.DOM.tr(
                Bonito.DOM.td(string(i)), Bonito.DOM.td(s.type), Bonito.DOM.td(string(s.curvature)),
                Bonito.DOM.td(string(s.distance)), Bonito.DOM.td(s.material), Bonito.DOM.td(string(s.radius)),
                Bonito.DOM.td(string(s.stop)), Bonito.DOM.td(string(s.conic)), Bonito.DOM.td(s.coating),
                Bonito.DOM.td(s.comm),
            ) for (i, s) in enumerate(summary.zsurfs)]
    return Bonito.DOM.div(header, Bonito.DOM.table(headerRow, rows...))
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
            style = "white-space:pre-wrap; word-break:break-word; " *
                    "max-width:100%; border:1px solid #ccc; padding:6px;"),
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
                        style = "margin-top:16px;"),
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
                                     style = "width:300px; flex:0 0 auto; overflow-y:auto;")
        detailColumns = Bonito.map!(Bonito.Observable{Any}(Bonito.DOM.div()), selectedPath) do p
            middle, preview = contentPane(p === nothing ? "" : p)
            return Bonito.DOM.div(
                Bonito.DOM.div(middle; style = "width:360px; flex:0 0 auto; overflow-y:auto; padding:0 12px;"),
                Bonito.DOM.div(preview; style = "flex:1 1 auto; overflow-y:auto; padding:0 12px;");
                style = "display:flex; flex:1 1 auto;",
            )
        end
        return Bonito.DOM.div(fileColumn, detailColumns; style = "display:flex; align-items:flex-start;")
    end
end

function openInBrowser(url::String)
    try
        if Sys.isapple()
            run(`open $url`)
        elseif Sys.iswindows()
            run(`cmd /c start $url`)
        else
            run(`xdg-open $url`)
        end
    catch e
        @warn "Could not automatically open a browser tab; visit $url manually" exception = e
    end
    return nothing
end

"""
    zemaxBrowser(rootDir::String; port::Integer=8081, openBrowser::Bool=true,
                 closeOnDisconnect::Bool=true) -> Bonito.Server

Launch a standalone Zemax file browser: starts a `Bonito.Server` on
`127.0.0.1:port` serving [`zemaxBrowserApp`](@ref)`(rootDir)`, optionally
opens a browser tab pointing at it (`openBrowser`), and returns the live
`Server` immediately without blocking. By default (`closeOnDisconnect=true`,
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
