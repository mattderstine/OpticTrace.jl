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
    ZemaxFileNode

One node in a directory tree walked by [`walkZemaxDirectory`](@ref).

Fields:
    path::String                     - absolute path
    name::String                     - basename(path)
    isDir::Bool                      - whether this node is a directory
    kind::Symbol                     - :dir, :zmx, :zar, or :zmf
    children::Vector{ZemaxFileNode}  - child nodes, empty unless isDir
"""
struct ZemaxFileNode
    path::String
    name::String
    isDir::Bool
    kind::Symbol
    children::Vector{ZemaxFileNode}
end

"""
    zemaxFileKind(path::AbstractString) -> Union{Symbol,Nothing}

Classify `path` by extension (case-insensitive): `:zmx`, `:zar`, `:zmf`,
or `nothing` for anything else. Used both to filter [`walkZemaxDirectory`](@ref)'s
tree and to dispatch content-pane rendering in the browser UI.
"""
function zemaxFileKind(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    ext == ".zmx" && return :zmx
    ext == ".zar" && return :zar
    ext == ".zmf" && return :zmf
    return nothing
end

"""
    walkZemaxDirectory(rootDir::String) -> ZemaxFileNode

Recursively walk `rootDir`, returning a [`ZemaxFileNode`](@ref) tree.
Every subdirectory is included as a `:dir` node regardless of whether it
(transitively) contains a browsable file -- a plain file-tree browser, not
a filtered view. Only files whose [`zemaxFileKind`](@ref) is non-`nothing`
become leaf nodes; other files are silently skipped. Symlink loops are not
handled.
"""
function walkZemaxDirectory(rootDir::String)::ZemaxFileNode
    children = ZemaxFileNode[]
    for entry in sort(readdir(rootDir))
        path = joinpath(rootDir, entry)
        if isdir(path)
            push!(children, walkZemaxDirectory(path))
        else
            kind = zemaxFileKind(path)
            kind === nothing && continue
            push!(children, ZemaxFileNode(path, entry, false, kind, ZemaxFileNode[]))
        end
    end
    name = basename(rootDir)
    isempty(name) && (name = rootDir)
    return ZemaxFileNode(rootDir, name, true, :dir, children)
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

_isZmxName(name::AbstractString) = endswith(lowercase(name), ".zmx")

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
        return [ZemaxArchiveEntity(e.name, length(e.data), _isZmxName(e.name))
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

The default extraction output directory for `archivePath`, exactly as
[`extractZemaxArchive`](@ref)/[`extractZmfCatalog`](@ref) would compute it
themselves (via the same `_stripExtension` helper) -- used to pre-fill the
browser UI's editable output-path field with a value that matches what
extraction would default to on its own.
"""
function defaultExtractionOutputPath(archivePath::String)::String
    ext = lowercase(splitext(archivePath)[2])
    if ext == ".zar"
        return _stripExtension(archivePath, ".zar")
    elseif ext == ".zmf"
        return _stripExtension(archivePath, ".zmf")
    else
        error("defaultExtractionOutputPath: not a .zar or .zmf archive: $archivePath")
    end
end

#=
    Bonito UI layer below -- every Bonito symbol is qualified (Bonito.App,
    Bonito.DOM, ...), see the module-collision note at the top of this file.
=#

function _kindLabel(kind::Symbol)
    kind === :zmx && return "[zmx] "
    kind === :zar && return "[zar] "
    kind === :zmf && return "[zmf] "
    return ""
end

"""
    _treeNode(node::ZemaxFileNode, selectedPath::Bonito.Observable{String})

Recursively render `node` as a Bonito DOM node: directories as an always-
expanded labeled group (no collapse/expand state in v1), files as a
clickable button that sets `selectedPath[]` to the file's path.
"""
function _treeNode(node::ZemaxFileNode, selectedPath::Bonito.Observable{String})
    if node.isDir
        childNodes = [_treeNode(c, selectedPath) for c in node.children]
        return Bonito.DOM.div(
            Bonito.DOM.div(node.name; style = "font-weight:bold; margin-top:4px;"),
            Bonito.DOM.div(childNodes...; style = "margin-left:14px;"),
        )
    else
        button = Bonito.Button(_kindLabel(node.kind) * node.name; style = nothing)
        Bonito.on(button.value) do _clicked
            selectedPath[] = node.path
        end
        return Bonito.DOM.div(button; style = "margin:1px 0;")
    end
end

"""
    _renderZemaxFileSummary(summary::ZemaxFileSummary)

Render a [`ZemaxFileSummary`](@ref) as a Bonito DOM node: name, units, the
populated (nonzero) wavelengths, and a table with one row per surface
showing the same fields [`printZemaxSurfs`](@ref) already treats as the
interesting ones (type, curvature, distance, material, radius, stop,
conic, coating, comm). Used both for a `.zmx` file found directly in the
tree and for a previewable `.zmx` entity found inside an archive.
"""
function _renderZemaxFileSummary(summary::ZemaxFileSummary)
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
    _renderNonPreviewableEntity(entity::ZemaxArchiveEntity)

Render a non-previewable archive entity as name + byte size + a "preview
not available" note -- no attempt to interpret its content.
"""
function _renderNonPreviewableEntity(entity::ZemaxArchiveEntity)
    return Bonito.DOM.div(
        Bonito.DOM.h3(entity.name),
        Bonito.DOM.div("Size: $(entity.byteSize) bytes"),
        Bonito.DOM.div("Preview not available for this entity type."),
    )
end

"""
    _archiveContentPane(path::String)

Render the content pane for a `.zar`/`.zmf` archive at `path`: a list of
its entities (each with a "select" and an "Extract" button), an editable
output-path field pre-filled via [`defaultExtractionOutputPath`](@ref), an
"Extract All" button, a status line, and a detail sub-pane driven by which
entity is currently selected (reset on every re-render of this pane, i.e.
whenever the outer `selectedPath` changes to a different archive).
"""
function _archiveContentPane(path::String)
    entities = listZemaxArchiveEntities(path)
    selectedEntity = Bonito.Observable{Union{Nothing,String}}(nothing)
    status = Bonito.Observable("")
    outputField = Bonito.TextField(defaultExtractionOutputPath(path); style = nothing)

    entityRows = map(entities) do ent
        selectBtn = Bonito.Button(ent.previewable ? ent.name : "$(ent.name) (no preview)"; style = nothing)
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
            status[] = "Extracted $(length(paths)) file(s) to: $(outputField.value[])"
        catch e
            status[] = "Error: $(sprint(showerror, e))"
        end
    end

    entityDetail = Bonito.map!(Bonito.Observable{Any}(Bonito.DOM.div()), selectedEntity) do name
        if name === nothing
            return Bonito.DOM.div("Select an entity above to preview it.")
        end
        idx = findfirst(e -> e.name == name, entities)
        ent = entities[idx]
        return ent.previewable ? _renderZemaxFileSummary(zemaxEntitySummary(path, name)) :
                                  _renderNonPreviewableEntity(ent)
    end

    return Bonito.DOM.div(
        Bonito.DOM.h3("Entities in $(basename(path))"),
        Bonito.DOM.div(entityRows...),
        Bonito.DOM.div("Output path: ", outputField, extractAllBtn),
        Bonito.DOM.div(status),
        entityDetail,
    )
end

"""
    _contentPane(path::String)

Top-level content-pane dispatch, by [`zemaxFileKind`](@ref)`(path)`:
a `.zmx` file renders via [`_renderZemaxFileSummary`](@ref); a `.zar`/
`.zmf` archive renders via [`_archiveContentPane`](@ref); an empty `path`
(nothing selected yet) renders a placeholder prompt.
"""
function _contentPane(path::String)
    isempty(path) && return Bonito.DOM.div("Select a file from the tree to view its contents.")
    kind = zemaxFileKind(path)
    kind === :zmx && return _renderZemaxFileSummary(zemaxFileSummary(path))
    (kind === :zar || kind === :zmf) && return _archiveContentPane(path)
    return Bonito.DOM.div("Unsupported file: $path")
end

"""
    _zemaxBrowserApp(rootDir::String) -> Bonito.App

Build the full browser app for `rootDir`: a fixed-width tree pane (via
[`_treeNode`](@ref)/[`walkZemaxDirectory`](@ref)) laid out as a sibling of
a flexible content pane (via [`_contentPane`](@ref)), driven by one
`selectedPath` observable per session. The two-column layout is
deliberate: a future WGLMakie render pane can be added later as a third
sibling `div` fed by its own observable, without restructuring this
tree/content wiring.
"""
function _zemaxBrowserApp(rootDir::String)
    return Bonito.App() do session
        selectedPath = Bonito.Observable("")
        tree = _treeNode(walkZemaxDirectory(rootDir), selectedPath)
        content = Bonito.map!(Bonito.Observable{Any}(Bonito.DOM.div()), selectedPath) do p
            _contentPane(p)
        end
        return Bonito.DOM.div(
            Bonito.DOM.div(tree; style = "width:320px; float:left; overflow-y:auto;"),
            Bonito.DOM.div(content; style = "margin-left:340px;"),
        )
    end
end

function _openInBrowser(url::String)
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
    zemaxBrowser(rootDir::String; port::Integer=8081, openBrowser::Bool=true) -> Bonito.Server

Launch a standalone Zemax file browser: starts a `Bonito.Server` on
`127.0.0.1:port` serving [`_zemaxBrowserApp`](@ref)`(rootDir)`, optionally
opens a browser tab pointing at it (`openBrowser`), and returns the live
`Server` immediately without blocking -- close it later with `close(server)`.
Throws an error if `rootDir` isn't a directory.
"""
function zemaxBrowser(rootDir::String; port::Integer = 8081, openBrowser::Bool = true)
    isdir(rootDir) || error("zemaxBrowser: not a directory: $rootDir")
    app = _zemaxBrowserApp(rootDir)
    server = Bonito.Server(app, "127.0.0.1", Int(port))
    openBrowser && _openInBrowser(Bonito.online_url(server, "/"))
    return server
end
