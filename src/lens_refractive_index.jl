#=

Functions to get refractive index from the RefractiveIndex.info website

=#


export getRefractiveIndexFunc, setDirectoryBaseRefractiveIndex, loadRICatalog, getDirectoryBaseRI
export loadRICatalog!

"""
    dirBaseRefractiveIndex

Mutable global holding the local filesystem root of a
RefractiveIndex.info-style glass-data directory tree (`.yml` files
organized as `<basepath>/glass/<manufacturer>/<glass>.yml`), used as
the default `basepath` for [`getRefractiveIndexFunc`](@ref)/
[`loadRICatalog!`](@ref). Currently hardcoded to a path on the original
author's machine -- not portable to other machines as-is. Read via
[`getDirectoryBaseRI`](@ref); change it via
[`setDirectoryBaseRefractiveIndex`](@ref) rather than assigning to it
directly.
"""
global dirBaseRefractiveIndex::String ="/Users/matt/Development/Projects/refractiveindex/database/data"

"""
    setDirectoryBaseRefractiveIndex()

Reset [`dirBaseRefractiveIndex`](@ref) to its hardcoded default.
"""
function setDirectoryBaseRefractiveIndex()
    global dirBaseRefractiveIndex="/Users/matt/Development/Projects/refractiveindex/database/data"
end

"""
    setDirectoryBaseRefractiveIndex(dir::AbstractString)

Set [`dirBaseRefractiveIndex`](@ref) to `dir`.
"""
function setDirectoryBaseRefractiveIndex(dir::AbstractString)
    global dirBaseRefractiveIndex= dir
end

"""
    getDirectoryBaseRI()

Return the current value of [`dirBaseRefractiveIndex`](@ref).
"""
getDirectoryBaseRI() = dirBaseRefractiveIndex

#=

fusedSilica = (0.696166300, 0.407942600, 0.897479400, 4.67914826E-3,
        1.35120631E-2, 97.9340025)

b270 = (93.2901, 11.7559, 1.27686, 206.221, -25.1637, 0.010476)

calciumFluoride = (0.443749998, 0.444930066, 0.150133991,
        0.00178027854, 0.00788536061, 0.0124119491)

magnesiumFluoride = (0.27620, 0.60967, 0.0080, 2.14973, 0.08636^2, 18.0^2, 25.0^2)

dlak6 = (2.81430513, -0.013449512, 0.0196569671, 0.000249252546, 1.98274573e-05, -8.02641015e-07)


libRI = YAML.load_file("/Users/matt/Development/Projects/refractiveIndex/database/library.yml")
=#


"""
    findRefractiveIndex(lambda, (b1, b2, b3, c1, c2, c3))

A fixed 3-term Sellmeier dispersion formula:
`n = sqrt(1 + b1*λ²/(λ²-c1) + b2*λ²/(λ²-c2) + b3*λ²/(λ²-c3))`, with
`λ` (`lambda`) in microns. Same mathematical form as
[`riFormula1`](@ref)/[`riFormula2`](@ref), but takes its 6 coefficients
as fixed positional arguments rather than a coefficient vector. Not
exported.
"""
findRefractiveIndex(lambda, (b1, b2, b3, c1, c2, c3)) =
    sqrt(1 + (b1 * lambda^2)/(lambda^2 - c1) +
    (b2 *lambda^2)/(lambda^2 - c2) + (b3 *lambda^2)/(lambda^2 - c3))

"""
    findRefractiveIndexAlt(lambda, (c1, c2, c3, c4, c5, c6))

A fixed power-series dispersion formula:
`n = sqrt(c1 + c2*λ² + c3*λ⁻² + c4*λ⁻⁴ + c5*λ⁻⁶ + c6*λ⁻⁸)`, with `λ`
(`lambda`) in microns. Not exported.
"""
findRefractiveIndexAlt(lambda, (c1, c2, c3, c4, c5, c6))=
    sqrt(c1 + c2*lambda^2 + c3 * lambda^-2 +c4 * lambda^-4 +
        c5 * lambda^-6+ c6*lambda^-8)


"""
    riFormula2(λ::T, c::S) where {T<:Real, S<:AbstractVector{T}}

Dispersion formula: `n = sqrt(1 + c[1] + Σ c[i]*λ²/(λ²-c[i+1]))` summed
over `i = 2, 4, 6, ...`, with `λ` in microns. Corresponds to the
`"formula 2"` dispersion type in a RefractiveIndex.info-style glass
YAML record's `DATA[1]["type"]` field -- selected by
[`getRefractiveIndexFunc`](@ref) via that string; not usually called
directly by name.

Compare [`riFormula1`](@ref) (identical shape, but the resonance term
`c[i+1]` is squared) and [`riFormula3`](@ref) (a power-law sum instead
of resonance-pole terms).

Arguments:
- `λ` -- wavelength, in microns
- `c` -- dispersion coefficients, as read from the glass YAML file's
  `"coefficients"` field
"""
function riFormula2(λ::T, c::S) where {T<:Real, S<:AbstractVector{T}}
    λ2 = λ^2
    sum = 1.0 + c[1]
    for i in 2:2:length(c)
        sum += c[i] * λ2/(λ2 - c[i+1])
    end
    sqrt(sum)
end

"""
    riFormula1(λ::T, c::S) where {T<:Real, S<:AbstractVector{T}}

Dispersion formula: `n = sqrt(1 + c[1] + Σ c[i]*λ²/(λ²-c[i+1]^2))`
summed over `i = 2, 4, 6, ...`, with `λ` in microns -- the same shape
as [`riFormula2`](@ref), except the resonance term `c[i+1]` is squared
here. Corresponds to the `"formula 1"` dispersion type in a
RefractiveIndex.info-style glass YAML record; selected by
[`getRefractiveIndexFunc`](@ref) via that string.

Arguments:
- `λ` -- wavelength, in microns
- `c` -- dispersion coefficients, as read from the glass YAML file's
  `"coefficients"` field
"""
function riFormula1(λ::T, c::S) where {T<:Real, S<:AbstractVector{T}}
    λ2 = λ^2
    sum = 1.0 + c[1]
    for i in 2:2:length(c)
        sum += c[i] * λ2/(λ2 - c[i+1]^2)
    end
    sqrt(sum)
end

"""
    riFormula3(λ::T, c::S) where {T<:Real, S<:AbstractVector{T}}

Dispersion formula: `n = sqrt(c[1] + Σ c[i]*λ^c[i+1])` summed over
`i = 2, 4, 6, ...`, with `λ` in microns. Corresponds to the
`"formula 3"` dispersion type in a RefractiveIndex.info-style glass
YAML record; selected by [`getRefractiveIndexFunc`](@ref) via that
string. Compare [`riFormula1`](@ref)/[`riFormula2`](@ref) (resonance-pole
terms instead of power-law terms).

Arguments:
- `λ` -- wavelength, in microns
- `c` -- dispersion coefficients, as read from the glass YAML file's
  `"coefficients"` field
"""
function riFormula3(λ::T, c::S) where {T<:Real, S<:AbstractVector{T}}

    sum = c[1]
    for i in 2:2:length(c)
        sum += c[i] * λ ^ c[i+1]
    end
    sqrt(sum)
end

#=
riN_LAK22(l) = riFormula2(l,
        map(x->parse(Float64,x),
            split(YAML.load_file("/Users/matt/Development/Projects/refractiveIndex/N-LAK22.yml")["DATA"][1]["coefficients"]," ")
            )
        )

riN_LAK22(0.5)
=#



"""
    getRefractiveIndexFunc(basepath::AbstractString, path::AbstractString)

Build a `wavelength -> refractive index` function for one glass, by
reading its RefractiveIndex.info-style YAML record at
`joinpath(basepath, path)`, dispatching on its `DATA[1]["type"]` field
to [`riFormula1`](@ref)/[`riFormula2`](@ref)/[`riFormula3`](@ref), and
parsing its `"coefficients"` field.

Returns `nothing` if the record has no `"coefficients"` entry (e.g. a
tabulated-data-only record with no dispersion formula); otherwise
returns a 1-argument function `wavelength -> index`, suitable for use
as a [`defaultGlassCatalog`](@ref) entry.

Arguments:
- `basepath` -- root of the glass-data directory tree
- `path` -- glass YAML file path, relative to `basepath`
"""
function getRefractiveIndexFunc(basepath::AbstractString, path::AbstractString)
    funcs = Dict("formula 2"=>riFormula2, "formula 1"=> riFormula1, "formula 3"=>riFormula3)
    pth = joinpath(basepath,path)
    record = YAML.load_file(pth)
    data = record["DATA"][1]
    typeData = data["type"]
    coefentries = get(data, "coefficients", nothing)
    if coefentries == nothing
        return nothing
    end
    coefs = map(x->parse(Float64,x),split(coefentries))
    f(x) = funcs[typeData](x, coefs)
    f
end

"""
    getRefractiveIndexFunc(path::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)

Same as `getRefractiveIndexFunc(basepath, path)`, defaulting `basepath`
to the global [`dirBaseRefractiveIndex`](@ref).

Arguments:
- `path` -- glass YAML file path, relative to `basepath`
- `basepath` -- root of the glass-data directory tree (keyword,
  defaults to [`dirBaseRefractiveIndex`](@ref))
"""
function getRefractiveIndexFunc(path::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)
    getRefractiveIndexFunc(basepath, path)
end


#riN_LAK22 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-LAK22.yml")
#riN_SF6 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-SF6.yml")
#riN_SF2 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-SF2.yml")
#riN_BK7 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-BK7.yml")


"""
    loadRICatalog(pth::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)

Build a fresh glass catalog (a `Dict{AbstractString,Any}` mapping glass
names to `wavelength -> index` functions, seeded with a `"DEFAULT"`
entry) by walking `joinpath(basepath, pth)` -- delegates to
[`loadRICatalog!`](@ref).

Arguments:
- `pth` -- directory to walk, relative to `basepath` (e.g. a
  manufacturer subdirectory like `"glass/schott"`)
- `basepath` -- root of the glass-data directory tree (keyword,
  defaults to [`dirBaseRefractiveIndex`](@ref))
"""
function loadRICatalog(pth::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)
    dictforpath = Dict{AbstractString, Any}()
    dictforpath["DEFAULT"]= rInDef
    loadRICatalog!(dictforpath, pth::AbstractString; basepath)
end

"""
    loadRICatalog!(dictforpath::Dict{AbstractString, Any}, pth::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)

Walk `joinpath(basepath, pth)` and, for every `.yml` file found, add a
`glassname -> wavelength->index` entry to `dictforpath` (via
[`getRefractiveIndexFunc`](@ref)), skipping files with no dispersion
formula. Files whose name splits into exactly two dash-separated parts
(e.g. `"SCHOTT-BK7.yml"`) also get a second entry keyed by the part
after the dash (e.g. `"BK7"`), so such a glass can be looked up with or
without its manufacturer prefix. Non-`.yml` files are reported and
skipped.

Returns `dictforpath`, mutated in place.

Arguments:
- `dictforpath` -- the catalog dict to add entries to
- `pth` -- directory to walk, relative to `basepath`
- `basepath` -- root of the glass-data directory tree (keyword,
  defaults to [`dirBaseRefractiveIndex`](@ref))
"""
function loadRICatalog!(dictforpath::Dict{AbstractString, Any}, pth::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)
    path = joinpath(basepath,pth)

    for (path, dirs, files) in walkdir(path)
        #println("Files in $path")
        for file in files
            #println(file) # path to files
            glassname, ext = splitext(file)
            if ext == ".yml"
                riFunc = getRefractiveIndexFunc(path, file)
                if riFunc==nothing
                    println("Skipping $(joinpath(path,file)): no dispersion function found")
                    continue
                end
                dictforpath[glassname] = riFunc
                sfile = split(glassname, "-")
                if length(sfile)==2
                    dictforpath[sfile[2]] = riFunc #remove the manufacturer specific prefix and add the glass
                end

            else
                println("Additional file found: $(joinpath(path,file))")
            end

        end
    end
    return dictforpath
end

"""
    loadRICatalog!(pth::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)

Same as `loadRICatalog!(dictforpath, pth; basepath)`, but adds entries
directly to the global [`defaultGlassCatalog`](@ref) instead of a
caller-supplied dict.

Arguments:
- `pth` -- directory to walk, relative to `basepath`
- `basepath` -- root of the glass-data directory tree (keyword,
  defaults to [`dirBaseRefractiveIndex`](@ref))
"""
function loadRICatalog!(pth::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)
    loadRICatalog!(defaultGlassCatalog, pth; basepath)
end
