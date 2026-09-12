#=

Functions to load glasses from Zemax .agf glass-catalog files into a
glass dictionary of the same shape as loadRICatalog! builds
(src/lens_refractive_index.jl).

=#

export loadAGFCatalog, loadAGFCatalog!, getAGFRefractiveIndexFunc, readAGFRecords

"""
    readAGFRecords(path::AbstractString)

Parse a Zemax `.agf` glass-catalog text file at `path` into one
`(name, dispform, coefficients)` named tuple per glass record. Each
record starts with an `NM` line (`NM <name> <dispform> ...`) and runs
until the next `NM` line or end of file; only the record's `CD` line
(whitespace-separated dispersion coefficients) is read -- every other
tag (`GC`, `ED`, `TD`, `OD`, `LD`, `IT`, `MW`, ...) is present in real
`.agf` files but ignored here without error, since none of it is needed
to evaluate the dispersion formulas this package supports.

Returns a `Vector` of `(name::String, dispform::Int,
coefficients::Vector{Float64})` named tuples, in file order.
"""
function readAGFRecords(path::AbstractString)
    records = NamedTuple{(:name, :dispform, :coefficients),Tuple{String,Int,Vector{Float64}}}[]
    name = ""
    dispform = 0
    coefficients = Float64[]
    havePendingRecord = false

    function flush!()
        if havePendingRecord
            push!(records, (name = name, dispform = dispform, coefficients = coefficients))
        end
    end

    for line in eachline(path)
        tokens = split(strip(line))
        isempty(tokens) && continue
        tag = tokens[1]
        if tag == "NM"
            flush!()
            name = tokens[2]
            dispform = parse(Int, tokens[3])
            coefficients = Float64[]
            havePendingRecord = true
        elseif tag == "CD" && havePendingRecord
            coefficients = map(x -> parse(Float64, x), tokens[2:end])
        end
    end
    flush!()
    return records
end

"""
    getAGFRefractiveIndexFunc(dispform::Integer, coefficients::AbstractVector{<:Real})

Build a `wavelength -> refractive index` function from one `.agf`
glass record's dispersion-formula code (`dispform`, from its `NM` line)
and coefficients (from its `CD` line), reusing this package's existing
[`riFormula1`](@ref)/[`riFormula3`](@ref) dispersion math rather than
implementing AGF's formulas from scratch:

- `dispform == 1` (Schott: `n² = a0 + a1λ² + a2λ⁻² + a3λ⁻⁴ + a4λ⁻⁶ +
  a5λ⁻⁸`) is the exact special case of `riFormula3`'s
  `n² = c[1] + Σ c[i]·λ^c[i+1]` with fixed exponents
  `(2, -2, -4, -6, -8)`, via `c = [a0, a1,2.0, a2,-2.0, a3,-4.0,
  a4,-6.0, a5,-8.0]`.
- `dispform == 2` (Sellmeier1: `n² - 1 = K1λ²/(λ²-L1) + K2λ²/(λ²-L2) +
  K3λ²/(λ²-L3)`) is the exact special case of `riFormula1`'s
  `n² = 1 + c[1] + Σ c[i]·λ²/(λ²-c[i+1]²)` with `c[1] = 0` and
  `c[i+1] = sqrt(L_i)`, via `c = [0.0, K1,sqrt(L1), K2,sqrt(L2),
  K3,sqrt(L3)]`. This assumes each `L_i` is non-negative, true for
  ordinary Sellmeier1 glasses.

Any other `dispform` code (Sellmeier2/3/4/5, Herzberger, Conrady,
Handbook of Optics 1/2, Extended1/2/3 -- see `TODO.md` #37) is not
supported and returns `nothing`, mirroring
[`getRefractiveIndexFunc`](@ref)'s "`nothing` means skip" contract.
"""
function getAGFRefractiveIndexFunc(dispform::Integer, coefficients::AbstractVector{<:Real})
    if dispform == 1
        a0, a1, a2, a3, a4, a5 = coefficients[1:6]
        c = [a0, a1, 2.0, a2, -2.0, a3, -4.0, a4, -6.0, a5, -8.0]
        return λ -> riFormula3(λ, c)
    elseif dispform == 2
        K1, L1, K2, L2, K3, L3 = coefficients[1:6]
        c = [0.0, K1, sqrt(L1), K2, sqrt(L2), K3, sqrt(L3)]
        return λ -> riFormula1(λ, c)
    else
        return nothing
    end
end

"""
    loadAGFCatalog!(dictforpath::Dict{AbstractString, Any}, pth::AbstractString)

Add a `glassname -> wavelength->index` entry to `dictforpath` for every
glass record found in the `.agf` catalog at `pth`, via
[`readAGFRecords`](@ref) and [`getAGFRefractiveIndexFunc`](@ref). `pth`
may be a single `.agf` file, or a directory, which is walked (like
[`loadRICatalog!`](@ref)) for every `.agf`/`.AGF` file found. Glasses
whose dispersion-formula code isn't supported are skipped, printing a
message naming the glass and its unsupported code.

Returns `dictforpath`, mutated in place.

Arguments:
- `dictforpath` -- the catalog dict to add entries to
- `pth` -- an `.agf` file, or a directory of `.agf` files
"""
function loadAGFCatalog!(dictforpath::Dict{AbstractString,Any}, pth::AbstractString)
    agfFiles = if isfile(pth)
        [pth]
    else
        found = String[]
        for (dir, _, files) in walkdir(pth)
            for file in files
                if lowercase(splitext(file)[2]) == ".agf"
                    push!(found, joinpath(dir, file))
                end
            end
        end
        found
    end

    for agfFile in agfFiles
        for record in readAGFRecords(agfFile)
            f = getAGFRefractiveIndexFunc(record.dispform, record.coefficients)
            if f === nothing
                println("Skipping $(record.name) in $agfFile: unsupported dispersion formula code $(record.dispform)")
                continue
            end
            dictforpath[record.name] = f
        end
    end
    return dictforpath
end

"""
    loadAGFCatalog!(pth::AbstractString)

Same as `loadAGFCatalog!(dictforpath, pth)`, but adds entries directly
to the global [`defaultGlassCatalog`](@ref) instead of a caller-supplied
dict.

Arguments:
- `pth` -- an `.agf` file, or a directory of `.agf` files
"""
function loadAGFCatalog!(pth::AbstractString)
    loadAGFCatalog!(defaultGlassCatalog, pth)
end

"""
    loadAGFCatalog(pth::AbstractString)

Build a fresh glass catalog (a `Dict{AbstractString,Any}` mapping glass
names to `wavelength -> index` functions, seeded with a `"DEFAULT"`
entry) from the `.agf` catalog at `pth` -- delegates to
[`loadAGFCatalog!`](@ref).

Arguments:
- `pth` -- an `.agf` file, or a directory of `.agf` files
"""
function loadAGFCatalog(pth::AbstractString)
    dictforpath = Dict{AbstractString,Any}()
    dictforpath["DEFAULT"] = rInDef
    loadAGFCatalog!(dictforpath, pth)
end
