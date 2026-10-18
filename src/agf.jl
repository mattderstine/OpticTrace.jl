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

The `dispform` field is read as a `Float64` and then converted to an
`Int`, since some vendor catalogs write it as a float (e.g.
`1.00000000E+00` instead of `1`); a value that isn't a whole number
(e.g. `1.5`) is rejected.

A record whose `NM` line's `dispform` field or `CD` line's coefficients
can't be parsed as the expected type is skipped with a printed message naming the source
file, line number, and offending line -- consistent with how
[`loadAGFCatalog!`](@ref) skips glasses with an unsupported dispersion
formula code -- rather than throwing and aborting the whole catalog
load.

Returns a `Vector` of `(name::String, dispform::Int,
coefficients::Vector{Float64})` named tuples, in file order.
"""
function readAGFRecords(path::AbstractString)
    records = NamedTuple{(:name, :dispform, :coefficients),Tuple{String,Int,Vector{Float64}}}[]
    name = ""
    dispform = 0
    coefficients = Float64[]
    havePendingRecord = false
    skipPendingRecord = false

    function flush!()
        if havePendingRecord && !skipPendingRecord
            push!(records, (name = name, dispform = dispform, coefficients = coefficients))
        end
    end

    for (lineNumber, line) in enumerate(eachline(path))
        tokens = split(strip(line))
        isempty(tokens) && continue
        tag = tokens[1]
        if tag == "NM"
            flush!()
            name = tokens[2]
            coefficients = Float64[]
            havePendingRecord = true
            skipPendingRecord = false
            try
                dispform = Int(parse(Float64, tokens[3]))
            catch e
                e isa Union{ArgumentError,InexactError} || rethrow()
                println("Skipping $name in $path (line $lineNumber): cannot parse dispersion-formula code as Int: $line")
                skipPendingRecord = true
            end
        elseif tag == "CD" && havePendingRecord && !skipPendingRecord
            try
                coefficients = map(x -> parse(Float64, x), tokens[2:end])
            catch e
                e isa ArgumentError || rethrow()
                println("Skipping $name in $path (line $lineNumber): cannot parse dispersion coefficients as Float64: $line")
                skipPendingRecord = true
            end
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
[`riFormula2`](@ref)/[`riFormula3`](@ref)/[`riFormula4`](@ref)/[`riFormula5`](@ref)/
[`riFormula7`](@ref) dispersion math rather than
implementing AGF's formulas from scratch:

- `dispform == 1` (Schott: `n² = a0 + a1λ² + a2λ⁻² + a3λ⁻⁴ + a4λ⁻⁶ +
  a5λ⁻⁸`) is the exact special case of `riFormula3`'s
  `n² = c[1] + Σ c[i]·λ^c[i+1]` with fixed exponents
  `(2, -2, -4, -6, -8)`, via `c = [a0, a1,2.0, a2,-2.0, a3,-4.0,
  a4,-6.0, a5,-8.0]`.
- `dispform == 2` (Sellmeier1: `n² - 1 = K1λ²/(λ²-L1) + K2λ²/(λ²-L2) +
  K3λ²/(λ²-L3)`) is the exact special case of `riFormula2`'s
  `n² = 1 + c[1] + Σ c[i]·λ²/(λ²-c[i+1])` with `c[1] = 0` and
  `c[i+1] = L_i` (unsquared), via `c = [0.0, K1,L1, K2,L2, K3,L3]`.
  No `sqrt` is taken, so a glass with a negative `L_i` is handled
  without error.
- `dispform == 5` (Conrady: `n = n0 + A/λ + B/λ^3.5`) is the exact
  special case of [`riFormula5`](@ref)'s `n = c[1] + Σ c[i]·λ^c[i+1]`
  (the Cauchy form, which gives `n` directly rather than `n²`), via
  `c = [n0, A,-1.0, B,-3.5]`.

- `dispform == 3` (Herzberger: `n = A + B·L + C·L² + D·λ² + E·λ⁴ +
  F·λ⁶`, `L = 1/(λ²-0.028)`) is [`riFormula7`](@ref) with
  `c = [A, B, C, D, E, F]`.
- `dispform == 4` (Sellmeier2: `n²-1 = A + B1λ²/(λ²-L1) + B2/(λ²-L2)`)
  is [`riFormula4`](@ref) with `c = [1+A, B1,2, L1,1, B2,0, L2,1]`.
- `dispform == 6`/`11` (Sellmeier3/Sellmeier5: Sellmeier1 with 4/5
  terms) are `riFormula2` with `c = [0, K1,L1, K2,L2, ...]`.
- `dispform == 7`/`8` (Handbook of Optics 1: `n² = A + B/(λ²-C) - Dλ²`;
  Handbook of Optics 2: `n² = A + Bλ²/(λ²-C) - Dλ²`) are `riFormula4`
  with its second pole term unused and its first power-law term
  `-D·λ²` (`c10 = -D`, `c11 = 2`).
- `dispform == 9` (Sellmeier4: `n² = A + Bλ²/(λ²-C) + Dλ²/(λ²-E)`) is
  `riFormula2` with `c = [A-1, B,C, D,E]`.
- `dispform == 10`/`12`/`13` (Extended/Extended2/Extended3) are
  `riFormula3`'s power-law sum with fixed exponent lists: Extended
  `(0, 2, -2, -4, -6, -8, -10, -12)`, Extended2
  `(0, 2, -2, -4, -6, -8, 4, 6)`, Extended3
  `(0, 2, 4, -2, -4, -6, -8, -10, -12)`.

Codes 3, 6, 7, 8, 9, 12 and 13 were checked against the catalog `nd`
of every such glass with a real `nd` in the local Zemax `Glasscat/`
(agreement ~1e-6). Codes 4, 10 and 11 had no checkable glass there, so
the formulas for 4 and 10 are written from Zemax's published
definitions and **unverified against real data**; 11 is only
spot-checked against KRS5's physical index. Treat them accordingly if
a glass using them looks wrong.

Any other `dispform` code returns `nothing`, mirroring
[`getRefractiveIndexFunc`](@ref)'s "`nothing` means skip" contract.
"""
function getAGFRefractiveIndexFunc(dispform::Integer, coefficients::AbstractVector{<:Real})
    if dispform == 1
        a0, a1, a2, a3, a4, a5 = coefficients[1:6]
        c = [a0, a1, 2.0, a2, -2.0, a3, -4.0, a4, -6.0, a5, -8.0]
        return λ -> riFormula3(λ, c)
    elseif dispform == 2
        K1, L1, K2, L2, K3, L3 = coefficients[1:6]
        c = [0.0, K1, L1, K2, L2, K3, L3]
        return λ -> riFormula2(λ, c)
    elseif dispform == 5
        n0, A, B = coefficients[1:3]
        c = [n0, A, -1.0, B, -3.5]
        return λ -> riFormula5(λ, c)
    elseif dispform == 3
        # Herzberger: n = A + B·L + C·L² + D·λ² + E·λ⁴ + F·λ⁶, L = 1/(λ²-0.028)
        c = Float64.(coefficients[1:6])
        return λ -> riFormula7(λ, c)
    elseif dispform == 4
        A, B1, L1, B2, L2 = coefficients[1:5]
        c = [1 + A, B1, 2.0, L1, 1.0, B2, 0.0, L2, 1.0]
        return λ -> riFormula4(λ, c)
    elseif dispform == 6 || dispform == 11
        nTerms = dispform == 6 ? 4 : 5
        c = [0.0; Float64.(coefficients[1:2nTerms])]
        return λ -> riFormula2(λ, c)
    elseif dispform == 7
        A, B, C, D = coefficients[1:4]
        c = [A, B, 0.0, C, 1.0, 0.0, 0.0, 0.0, 0.0, -D, 2.0]
        return λ -> riFormula4(λ, c)
    elseif dispform == 8
        A, B, C, D = coefficients[1:4]
        c = [A, B, 2.0, C, 1.0, 0.0, 0.0, 0.0, 0.0, -D, 2.0]
        return λ -> riFormula4(λ, c)
    elseif dispform == 9
        A, B, C, D, E = coefficients[1:5]
        c = [A - 1, B, C, D, E]
        return λ -> riFormula2(λ, c)
    elseif dispform == 10 || dispform == 12 || dispform == 13
        exps = dispform == 10 ? (0.0, 2.0, -2.0, -4.0, -6.0, -8.0, -10.0, -12.0) :
               dispform == 12 ? (0.0, 2.0, -2.0, -4.0, -6.0, -8.0, 4.0, 6.0) :
                                (0.0, 2.0, 4.0, -2.0, -4.0, -6.0, -8.0, -10.0, -12.0)
        c = Float64[coefficients[1]]
        for (a, e) in zip(coefficients[2:length(exps)], exps[2:end])
            push!(c, a, e)
        end
        return λ -> riFormula3(λ, c)
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
