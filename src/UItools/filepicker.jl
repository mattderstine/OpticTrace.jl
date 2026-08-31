#=
    a generic, reusable Bonito.jl file/directory picker component for
    GUIs built with Makie & Bonito.

    Interaction model is a classic single-pane, breadcrumb-driven file
    dialog: only the current directory's immediate children are listed
    (not a fully-expanded recursive tree -- that can't bound its own
    height), an up-arrow and a breadcrumb bar navigate around, and
    double-clicking a directory row navigates into it while
    double-clicking a file row selects-and-confirms in one step.

    Like src/zemax_browser.jl, this file `import`s Bonito rather than
    `using` it and qualifies every Bonito symbol -- see that file's
    header comment (and CLAUDE.md's "Bonito/GLMakie name collision"
    note) for why a bare `using Bonito` anywhere in this module would
    break src/plotting.jl's unqualified use of GLMakie's `Button`.
=#
import Bonito

export filePicker

"""
    FILE_PICKER_MODES

The valid `mode` values for [`filePicker`](@ref)/[`filePickerApp`](@ref):
`:file` (select one file), `:directory` (select one directory, the one
currently displayed), or `:multipleFiles` (select any number of files).
"""
const FILE_PICKER_MODES = (:file, :directory, :multipleFiles)

"""
    checkFilePickerMode(mode::Symbol)

Throw an `ArgumentError` unless `mode` is one of [`FILE_PICKER_MODES`](@ref).
"""
function checkFilePickerMode(mode::Symbol)
    mode in FILE_PICKER_MODES ||
        throw(ArgumentError("filePicker: mode must be one of $FILE_PICKER_MODES, got :$mode"))
    return nothing
end

"""
    FILE_PICKER_SORT_MODES

The valid `sortBy` values for [`listDirectory`](@ref)/[`filePicker`](@ref)/
[`filePickerApp`](@ref): `:name` (alphabetical by filename) or `:type`
(alphabetical by file-suffix/extension, ties broken by filename). Both
orderings are case-insensitive.
"""
const FILE_PICKER_SORT_MODES = (:name, :type)

"""
    checkFilePickerSortMode(sortBy::Symbol)

Throw an `ArgumentError` unless `sortBy` is one of [`FILE_PICKER_SORT_MODES`](@ref).
"""
function checkFilePickerSortMode(sortBy::Symbol)
    sortBy in FILE_PICKER_SORT_MODES ||
        throw(ArgumentError("listDirectory: sortBy must be one of $FILE_PICKER_SORT_MODES, got :$sortBy"))
    return nothing
end

"""
    listDirectory(dir::String; extensions=nothing, includeHidden::Bool=true, sortBy::Symbol=:name) -> (dirs::Vector{String}, files::Vector{String})

List the immediate children of `dir`: subdirectory names and file names,
each case-insensitively sorted. `extensions`, when given, is a collection
of lowercase extension strings (e.g. `[".zmx", ".zar"]`) used to filter
`files`; `nothing` (the default) includes every file. `includeHidden`,
when `false`, excludes any entry (directory or file) whose name starts
with `.` (e.g. `.git`, `.DS_Store`); `true` (the default) includes them,
preserving `readdir`'s own behavior. `sortBy` is one of
[`FILE_PICKER_SORT_MODES`](@ref): `:name` (the default) sorts both `dirs`
and `files` case-insensitively by filename; `:type` sorts both
case-insensitively by extension first, filename second -- directories
have no extension, so they always sort as one group by filename either
way. Only one directory level is read -- no recursion -- since the
picker UI navigates lazily, one level at a time, rather than walking the
whole tree up front. Throws an `ArgumentError` for an unrecognized
`sortBy`.
"""
function listDirectory(dir::String;
                        extensions::Union{Nothing,AbstractVector{<:AbstractString}} = nothing,
                        includeHidden::Bool = true, sortBy::Symbol = :name)
    checkFilePickerSortMode(sortBy)
    dirs = String[]
    files = String[]
    for entry in readdir(dir)
        (includeHidden || !startswith(entry, ".")) || continue
        path = joinpath(dir, entry)
        if isdir(path)
            push!(dirs, entry)
        elseif extensions === nothing || lowercase(splitext(entry)[2]) in extensions
            push!(files, entry)
        end
    end
    sortKey(name) = sortBy === :type ? (lowercase(splitext(name)[2]), lowercase(name)) : (lowercase(name),)
    return sort(dirs; by = sortKey), sort(files; by = sortKey)
end

"""
    pathBreadcrumbs(path::String) -> Vector{Pair{String,String}}

Split `path` into `label => fullPathThroughHere` pairs, one per path
segment, for breadcrumb-bar rendering. E.g. `/Users/matt/data` becomes
`["/" => "/", "Users" => "/Users", "matt" => "/Users/matt", "data" => "/Users/matt/data"]`.
`path` is normalized via `abspath` before splitting, so a relative path
still produces a well-formed, absolute breadcrumb chain.
"""
function pathBreadcrumbs(path::String)
    sep = Base.Filesystem.path_separator
    acc = sep
    crumbs = Pair{String,String}[sep => acc]
    for segment in split(abspath(path), sep; keepempty = false)
        acc = joinpath(acc, segment)
        push!(crumbs, segment => acc)
    end
    return crumbs
end

#=
    Bonito UI layer below -- every Bonito symbol is qualified (Bonito.App,
    Bonito.DOM, ...), see the module-collision note at the top of this file.
=#

const _UP_ARROW = "↑"

"""
    upButton(pathField::Bonito.TextField)

An up-arrow `Bonito.Button` that navigates `pathField` to its parent
directory (via `dirname`) on click -- a no-op at the filesystem root,
where `dirname` is idempotent.
"""
function upButton(pathField::Bonito.TextField)
    button = Bonito.Button(_UP_ARROW; style = nothing)
    Bonito.on(button.value) do _clicked
        pathField.value[] = dirname(pathField.value[])
    end
    return button
end

"""
    breadcrumbBar(pathField::Bonito.TextField)

Render `pathField`'s current directory as an up-arrow button (see
[`upButton`](@ref)) followed by [`pathBreadcrumbs`](@ref) rendered as
clickable path segments; clicking a segment navigates `pathField` to
that ancestor directory. Reactively rebuilt whenever `pathField.value`
changes.
"""
function breadcrumbBar(pathField::Bonito.TextField)
    return Bonito.map!(Bonito.Observable{Any}(Bonito.DOM.div()), pathField.value) do currentPath
        crumbButtons = map(pathBreadcrumbs(currentPath)) do crumb
            label, target = crumb
            button = Bonito.Button(label; style = nothing)
            Bonito.on(button.value) do _clicked
                pathField.value[] = target
            end
            return button
        end
        return Bonito.DOM.div(upButton(pathField), crumbButtons...;
                               style = "display:flex; flex-wrap:wrap; align-items:center;")
    end
end

"""
    directoryRow(dir::String, name::String, mode::Symbol, pathField::Bonito.TextField, pending)

Render one subdirectory row: a name-only `Bonito.Button`. Double-
clicking always navigates `pathField` into `dir/name`; in `:directory`
mode a single click additionally selects `dir/name` into `pending[]`
(navigate-only, no selection, in `:file`/`:multipleFiles` mode).
"""
function directoryRow(dir::String, name::String, mode::Symbol,
                        pathField::Bonito.TextField, pending)
    childPath = joinpath(dir, name)
    dblclick = Bonito.Observable(false)
    Bonito.on(dblclick) do _dbl
        pathField.value[] = childPath
    end
    button = Bonito.Button(name; style = nothing,
                            ondblclick = Bonito.js"event=> $(dblclick).notify(true);")
    if mode === :directory
        Bonito.on(button.value) do _clicked
            pending[] = childPath
        end
    end
    return Bonito.DOM.div(button; style = "margin:1px 0;")
end

"""
    fileRow(dir::String, name::String, mode::Symbol, pending, selected, active)

Render one file row for `:file`/`:multipleFiles` mode (`:directory`
mode renders files as plain non-interactive text instead -- see
[`directoryListing`](@ref)). In `:file` mode a click selects `dir/name`
into `pending[]`; in `:multipleFiles` mode a `Bonito.Checkbox` toggles
`dir/name`'s membership in `pending[]::Vector{String}`. Either way, a
double-click short-circuits straight to a final selection: sets
`selected[]` to just this one file and `active[] = false`, bypassing
the Select button (the same "double-click to open" shortcut a native
file dialog gives you, even in multi-select mode).
"""
function fileRow(dir::String, name::String, mode::Symbol, pending, selected, active)
    path = joinpath(dir, name)
    dblclick = Bonito.Observable(false)
    Bonito.on(dblclick) do _dbl
        selected[] = mode === :multipleFiles ? [path] : path
        active[] = false
    end
    dblAttr = Bonito.js"event=> $(dblclick).notify(true);"
    if mode === :multipleFiles
        checkbox = Bonito.Checkbox(false; ondblclick = dblAttr)
        Bonito.on(checkbox.value) do checked
            current = pending[]
            pending[] = checked ? union(current, [path]) : setdiff(current, [path])
        end
        return Bonito.DOM.div(checkbox, " ", name; style = "margin:1px 0;")
    else
        button = Bonito.Button(name; style = nothing, ondblclick = dblAttr)
        Bonito.on(button.value) do _clicked
            pending[] = path
        end
        return Bonito.DOM.div(button; style = "margin:1px 0;")
    end
end

"""
    _SCROLLABLE_LIST_STYLE

Shared CSS for a bounded-height, scrolling list box -- used by
[`directoryListing`](@ref) here and reused by `zemax_browser.jl`'s
archive entity list, so both browsing UIs in this package present long
lists the same way.
"""
const _SCROLLABLE_LIST_STYLE = "max-height:400px; overflow-y:auto; border:1px solid #ccc; padding:4px;"

"""
    hiddenCheckboxRow(default::Bool) -> (row, showHidden::Bonito.Observable{Bool})

A "Hidden" `Bonito.Checkbox` with its label rendered to its left in a
smaller font, `default`-initialized. `row` is the DOM node to embed;
`showHidden` is the checkbox's own value `Observable` -- `true` (checked)
means dotfiles/dot-directories should be included in the listing, `false`
(unchecked) means they should be filtered out (see [`listDirectory`](@ref)'s
`includeHidden` kwarg, which this drives via [`directoryListing`](@ref)).
"""
function hiddenCheckboxRow(default::Bool)
    checkbox = Bonito.Checkbox(default)
    row = Bonito.DOM.div(
        Bonito.DOM.span("Hidden"; style = "font-size:0.8em; margin-right:4px;"),
        checkbox;
        style = "display:flex; align-items:center;",
    )
    return row, checkbox.value
end

"""
    sortDropdown(default::Symbol) -> (dropdown::Bonito.Dropdown, sortLabel::Bonito.Observable)

A `Bonito.Dropdown` offering "Name"/"Type" sort-order choices (see
[`FILE_PICKER_SORT_MODES`](@ref)), initially selecting "Type" iff
`default === :type` (otherwise "Name"). `sortLabel` is the dropdown's own
value `Observable`, holding whichever of `"Name"`/`"Type"` is currently
selected -- [`directoryListing`](@ref) maps that string back onto a
`sortBy` symbol for [`listDirectory`](@ref).
"""
function sortDropdown(default::Symbol)
    labels = ["Name", "Type"]
    dropdown = Bonito.Dropdown(labels; index = default === :type ? 2 : 1)
    return dropdown, dropdown.value
end

"""
    directoryListing(pathField::Bonito.TextField, mode::Symbol, extensions, pending, selected, active, showHidden::Bonito.Observable{Bool}, sortLabel)

Render the current directory's contents (via [`listDirectory`](@ref)):
subdirectories as [`directoryRow`](@ref)s, then files -- as
[`fileRow`](@ref)s in `:file`/`:multipleFiles` mode, or as plain
non-interactive text (context only) in `:directory` mode. Wrapped in a
bounded-height scrolling container so a directory with many entries
still fits on one screen. Reactively rebuilt whenever `pathField.value`,
`showHidden` (see [`hiddenCheckboxRow`](@ref)), or `sortLabel` (see
[`sortDropdown`](@ref)) changes. If `listDirectory` throws a
`Base.IOError` (e.g. a permission-denied directory -- macOS's
TCC-protected `Documents`/`Desktop`/etc. folders being the common case),
that's caught here and rendered as an in-place error message instead of
propagating -- letting the exception escape this callback would leave
the listing pane silently stuck on its previous contents (breadcrumbs
and the path field still update, since they don't depend on this
`Observable`, making navigation look like it silently does nothing).
"""
function directoryListing(pathField::Bonito.TextField, mode::Symbol, extensions,
                            pending, selected, active,
                            showHidden::Bonito.Observable{Bool}, sortLabel)
    return Bonito.map!(Bonito.Observable{Any}(Bonito.DOM.div()), pathField.value, showHidden, sortLabel) do currentPath, show, label
        sortBy = label == "Type" ? :type : :name
        local dirs, files
        try
            dirs, files = listDirectory(currentPath; extensions = extensions, includeHidden = show, sortBy = sortBy)
        catch e
            e isa Base.IOError || rethrow()
            return Bonito.DOM.div("Cannot read this directory: $(sprint(showerror, e))";
                                   style = "color:#b00020; padding:8px;")
        end
        dirRows = [directoryRow(currentPath, name, mode, pathField, pending) for name in dirs]
        fileRows = if mode === :directory
            [Bonito.DOM.div(name; style = "margin:1px 0; color:gray;") for name in files]
        else
            [fileRow(currentPath, name, mode, pending, selected, active) for name in files]
        end
        return Bonito.DOM.div(dirRows..., fileRows...; style = _SCROLLABLE_LIST_STYLE)
    end
end

"""
    filePickerControls(mode::Symbol, pathField::Bonito.TextField, pending, selected, active, hiddenRow)

The picker's final row: "Select"/"Cancel" on the left, `hiddenRow` (see
[`hiddenCheckboxRow`](@ref)) pushed to the right so the "Hidden" checkbox
sits at the bottom-right of the whole component. Select commits a result
into `selected[]` -- `pathField.value[]` (the currently-displayed
directory) in `:directory` mode, otherwise `pending[]` -- and sets
`active[] = false`. Cancel sets `selected[]` to `nothing` (`String[]` for
`:multipleFiles`) and `active[] = false`, without committing anything
from `pending[]`.
"""
function filePickerControls(mode::Symbol, pathField::Bonito.TextField, pending, selected, active, hiddenRow)
    selectButton = Bonito.Button("Select"; style = nothing)
    Bonito.on(selectButton.value) do _clicked
        selected[] = mode === :directory ? pathField.value[] : pending[]
        active[] = false
    end
    cancelButton = Bonito.Button("Cancel"; style = nothing)
    Bonito.on(cancelButton.value) do _clicked
        selected[] = mode === :multipleFiles ? String[] : nothing
        active[] = false
    end
    return Bonito.DOM.div(
        Bonito.DOM.div(selectButton, cancelButton),
        hiddenRow;
        style = "display:flex; justify-content:space-between; align-items:center;",
    )
end

"""
    filePicker(rootDir::String; mode::Symbol=:file, extensions=nothing, showHidden::Bool=false, sortBy::Symbol=:name) -> (component, selected, active)

Build an embeddable file/directory picker component rooted at
`rootDir`. `mode` is one of [`FILE_PICKER_MODES`](@ref):
- `:file` -- select a single file; `selected[]::Union{Nothing,String}`.
- `:directory` -- select a single directory (the one currently
  displayed); `selected[]::Union{Nothing,String}`.
- `:multipleFiles` -- select any number of files;
  `selected[]::Vector{String}`.

`extensions`, when given, restricts which files are listed (a
collection of lowercase extension strings, e.g. `[".zmx"]`).

`showHidden` sets the default (checked/unchecked) state of the "Hidden"
checkbox rendered at the bottom-right of the component (see
[`hiddenCheckboxRow`](@ref)): `false` (the default) hides
dotfiles/dot-directories; `true` shows them. The checkbox lets the user
toggle this live, driving [`listDirectory`](@ref)'s `includeHidden` kwarg
directly -- checked shows hidden entries, unchecked hides them.

`sortBy` is one of [`FILE_PICKER_SORT_MODES`](@ref) and sets the initial
selection of the "Name"/"Type" sort-order dropdown (see
[`sortDropdown`](@ref)) rendered to the right of the path textbox, in
its own row -- the dropdown lets the user switch live between the two
orderings [`listDirectory`](@ref) supports.

Returns `(component, selected, active)`: `component` is the DOM node to
embed in a caller's layout (this file has no notion of a standalone
popup window -- Bonito composes components within one session, the
same way `zemax_browser.jl`'s tree/content panes compose); `selected`
is the `Bonito.Observable` the caller subscribes to for the final
result; `active` is a `Bonito.Observable{Bool}`, initially `true` and
set to `false` once the user commits a selection or cancels -- a signal
the embedding caller can use to know when to remove/hide the picker.

Navigation: the path text field at the top can be typed into directly
(Enter navigates), double-clicking a directory row navigates into it,
the up-arrow/breadcrumb bar (see [`breadcrumbBar`](@ref)) navigates to
an ancestor directory, and -- in `:file`/`:multipleFiles` mode --
double-clicking a file row selects-and-confirms that one file
immediately, bypassing the Select button.

Throws an `ArgumentError` for an unrecognized `mode` or `sortBy`, or an
`ErrorException` if `rootDir` isn't a directory.
"""
function filePicker(rootDir::String; mode::Symbol = :file,
                     extensions::Union{Nothing,AbstractVector{<:AbstractString}} = nothing,
                     showHidden::Bool = false, sortBy::Symbol = :name)
    checkFilePickerMode(mode)
    checkFilePickerSortMode(sortBy)
    isdir(rootDir) || error("filePicker: not a directory: $rootDir")

    pathField = Bonito.TextField(abspath(rootDir); style = nothing)
    if mode === :multipleFiles
        pending = Bonito.Observable{Vector{String}}(String[])
        selected = Bonito.Observable{Vector{String}}(String[])
    else
        pending = Bonito.Observable{Union{Nothing,String}}(nothing)
        selected = Bonito.Observable{Union{Nothing,String}}(nothing)
    end
    active = Bonito.Observable(true)

    breadcrumbs = breadcrumbBar(pathField)
    hiddenRow, showHiddenObs = hiddenCheckboxRow(showHidden)
    sortDropdownWidget, sortLabelObs = sortDropdown(sortBy)
    pathRow = Bonito.DOM.div(
        Bonito.DOM.div(pathField; style = "flex:1 1 auto;"),
        sortDropdownWidget;
        style = "display:flex; align-items:center; gap:6px;",
    )
    listing = directoryListing(pathField, mode, extensions, pending, selected, active,
                                 showHiddenObs, sortLabelObs)
    controls = filePickerControls(mode, pathField, pending, selected, active, hiddenRow)

    component = Bonito.DOM.div(breadcrumbs, pathRow, listing, controls)
    return component, selected, active
end

"""
    filePickerApp(rootDir::String; mode::Symbol=:file, extensions=nothing, showHidden::Bool=false, sortBy::Symbol=:name) -> Bonito.App

Thin standalone-app wrapper around [`filePicker`](@ref), analogous to
`zemax_browser.jl`'s `_zemaxBrowserApp`/`zemaxBrowser`: gives the picker
a real `Bonito.App` to render on its own, for previewing or testing.
Discards the `selected`/`active` Observables `filePicker` returns -- a
real embedding caller should call `filePicker` directly and keep them.
"""
function filePickerApp(rootDir::String; mode::Symbol = :file,
                        extensions::Union{Nothing,AbstractVector{<:AbstractString}} = nothing,
                        showHidden::Bool = false, sortBy::Symbol = :name)
    return Bonito.App() do session
        component, _selected, _active = filePicker(rootDir; mode = mode, extensions = extensions,
                                                     showHidden = showHidden, sortBy = sortBy)
        return component
    end
end
