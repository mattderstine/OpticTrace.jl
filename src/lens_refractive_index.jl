#=

Functions to get refractive index from the RefractiveIndex.info website

=#


export getRefractiveIndexFunc, setDirectoryBaseRefractiveIndex, loadRICatalog, getDirectoryBaseRI
export loadRICatalog!, riFormula5, ModelGlass

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

"""
    riFormula5(λ::Real, c::AbstractVector{<:Real})

Dispersion formula: `n = c[1] + Σ c[i]*λ^c[i+1]` summed over
`i = 2, 4, 6, ...`, with `λ` in microns -- the same power-law sum as
[`riFormula3`](@ref), except it gives `n` directly rather than `n²` (no
`sqrt`). Corresponds to the `"formula 5"` (Cauchy) dispersion type in a
RefractiveIndex.info-style glass YAML record; selected by
[`getRefractiveIndexFunc`](@ref) via that string.

Also the shared math behind two Zemax dispersion models that are exact
special cases of this form: AGF dispersion-formula code 5 (Conrady,
`n = n0 + A/λ + B/λ^3.5`, see [`getAGFRefractiveIndexFunc`](@ref)) and
Zemax's Model Glass (a cubic in `1/λ`, see [`ModelGlass`](@ref)).

Unlike [`riFormula1`](@ref)/[`riFormula2`](@ref)/[`riFormula3`](@ref),
`λ` and `c` aren't required to share an element type, so either one can
carry `ForwardDiff` dual numbers independently of the other (e.g.
differentiating a `ModelGlass` index with respect to its `nd`, whose
coefficients are then duals while `λ` stays a `Float64`).

Arguments:
- `λ` -- wavelength, in microns
- `c` -- dispersion coefficients: `c[1]` is the constant term, followed
  by `(coefficient, exponent)` pairs
"""
function riFormula5(λ::Real, c::AbstractVector{<:Real})
    sum = c[1]
    for i in 2:2:length(c)
        sum += c[i] * λ ^ c[i+1]
    end
    sum
end

"""
    riCoef(c::AbstractVector, i::Integer)

Coefficient `c[i]`, or `zero(eltype(c))` when the glass record lists
fewer than `i` coefficients -- RefractiveIndex.info records may omit
trailing zero coefficients of the fixed-length formulas
([`riFormula4`](@ref), [`riFormula7`](@ref)).
"""
riCoef(c::AbstractVector, i::Integer) = i <= length(c) ? c[i] : zero(eltype(c))

"""
    riFormula4(λ::Real, c::AbstractVector{<:Real})

Dispersion formula (RefractiveIndex.info "formula 4"):
`n² = c1 + c2·λ^c3/(λ²-c4^c5) + c6·λ^c7/(λ²-c8^c9) + c10·λ^c11 +
c12·λ^c13 + c14·λ^c15 + c16·λ^c17`, with `λ` in microns. Missing
trailing coefficients are taken as zero, and a pole term whose
numerator coefficient (`c2`/`c6`) is zero is skipped, so unused terms
can't produce `0/0`.

Arguments:
- `λ` -- wavelength, in microns
- `c` -- dispersion coefficients, as read from the glass YAML file
"""
function riFormula4(λ::Real, c::AbstractVector{<:Real})
    λ2 = λ^2
    sum = riCoef(c, 1)
    for i in (2, 6)
        ci = riCoef(c, i)
        iszero(ci) && continue
        sum += ci * λ ^ riCoef(c, i+1) / (λ2 - riCoef(c, i+2) ^ riCoef(c, i+3))
    end
    for i in (10, 12, 14, 16)
        ci = riCoef(c, i)
        iszero(ci) && continue
        sum += ci * λ ^ riCoef(c, i+1)
    end
    sqrt(sum)
end

"""
    riFormula6(λ::Real, c::AbstractVector{<:Real})

Dispersion formula (RefractiveIndex.info "formula 6", gases):
`n = 1 + c[1] + Σ c[i]/(c[i+1] - λ⁻²)` summed over `i = 2, 4, 6, ...`,
with `λ` in microns.

Arguments:
- `λ` -- wavelength, in microns
- `c` -- dispersion coefficients, as read from the glass YAML file
"""
function riFormula6(λ::Real, c::AbstractVector{<:Real})
    sum = 1 + c[1]
    for i in 2:2:length(c)-1
        sum += c[i] / (c[i+1] - λ^-2)
    end
    sum
end

"""
    riFormula7(λ::Real, c::AbstractVector{<:Real})

Dispersion formula (RefractiveIndex.info "formula 7", Herzberger):
`n = c1 + c2/(λ²-0.028) + c3/(λ²-0.028)² + c4·λ² + c5·λ⁴ + c6·λ⁶`,
with `λ` in microns. A missing `c6` is taken as zero.

Arguments:
- `λ` -- wavelength, in microns
- `c` -- dispersion coefficients, as read from the glass YAML file
"""
function riFormula7(λ::Real, c::AbstractVector{<:Real})
    λ2 = λ^2
    d = λ2 - 0.028
    riCoef(c, 1) + riCoef(c, 2) / d + riCoef(c, 3) / d^2 +
        riCoef(c, 4) * λ2 + riCoef(c, 5) * λ2^2 + riCoef(c, 6) * λ2^3
end

"""
    riFormula8(λ::Real, c::AbstractVector{<:Real})

Dispersion formula (RefractiveIndex.info "formula 8", retro): with
`x = c1 + c2·λ²/(λ²-c3) + c4·λ²`, `(n²-1)/(n²+2) = x`, i.e.
`n = sqrt((1+2x)/(1-x))`, with `λ` in microns.

Arguments:
- `λ` -- wavelength, in microns
- `c` -- dispersion coefficients, as read from the glass YAML file
"""
function riFormula8(λ::Real, c::AbstractVector{<:Real})
    λ2 = λ^2
    x = c[1] + c[2] * λ2 / (λ2 - c[3]) + c[4] * λ2
    sqrt((1 + 2x) / (1 - x))
end

"""
    riFormula9(λ::Real, c::AbstractVector{<:Real})

Dispersion formula (RefractiveIndex.info "formula 9", exotic):
`n² = c1 + c2/(λ²-c3) + c4·(λ-c5)/((λ-c5)²+c6)`, with `λ` in microns.

Arguments:
- `λ` -- wavelength, in microns
- `c` -- dispersion coefficients, as read from the glass YAML file
"""
function riFormula9(λ::Real, c::AbstractVector{<:Real})
    d = λ - c[5]
    sqrt(c[1] + c[2] / (λ^2 - c[3]) + c[4] * d / (d^2 + c[6]))
end

"""
    riFormulas

Maps a RefractiveIndex.info `DATA[1]["type"]` string (`"formula 1"`
... `"formula 9"`) to its dispersion function, `riFormulaN(λ, c)`.
Used by [`getRefractiveIndexFunc`](@ref); tabulated record types have
no entry.
"""
const riFormulas = Dict{String,Function}(
    "formula 1" => riFormula1, "formula 2" => riFormula2, "formula 3" => riFormula3,
    "formula 4" => riFormula4, "formula 5" => riFormula5, "formula 6" => riFormula6,
    "formula 7" => riFormula7, "formula 8" => riFormula8, "formula 9" => riFormula9)

"""
    ModelGlass{T<:Real}

Zemax "Model Glass": an idealized glass defined only by its refractive
index at the d-line (587.5618 nm) `nd`, its Abbe number `vd`, and its
partial-dispersion deviation from the normal line `dPgF` -- the material
behind a Zemax `.zmx` surface whose `GLAS` line has solve type 1 (see
[`zemaxGlassIndexFunc`](@ref)). A `ModelGlass` is callable,
`(g::ModelGlass)(λ)`, so it can be used anywhere a glass catalog's
`wavelength -> index` function is expected.

The dispersion model follows OpticSim.jl's `modelglass` (MIT license):
the normal line from Schott's "TIE-29: Refractive Index and Dispersion"
gives the partial dispersion `a = dPgF + 0.6438 - 0.001682*vd`, and GNU
GOptical's `Material::Abbe` fit gives `C1..C4` linear in `a`, with

    n(λ) = nd + (nd - 1)/vd * (C1 + C2/λ + C3/λ² + C4/λ³)

evaluated via [`riFormula5`](@ref) with precomputed coefficients. The
fit is only meant for visible wavelengths, 0.36-0.75 μm; evaluating
outside that range warns once per session but still returns a value.
This approximates OpticStudio's own model glass and hasn't been checked
against it (see `TODO.md` #39).

Fields:
    nd::T                     - refractive index at 587.5618 nm
    vd::T                     - Abbe number at the d-line
    dPgF::T                   - deviation of partial dispersion P_g,F from the normal line
    coefficients::SVector{7,T} - precomputed `riFormula5` coefficients
"""
struct ModelGlass{T<:Real}
    nd::T
    vd::T
    dPgF::T
    coefficients::SVector{7,T}

    function ModelGlass{T}(nd, vd, dPgF) where {T<:Real}
        1.0 <= nd <= 2.5 || throw(ArgumentError("ModelGlass: nd = $nd is outside the plausible range [1.0, 2.5]"))
        10.0 <= vd <= 100.0 || throw(ArgumentError("ModelGlass: vd = $vd is outside the plausible range [10, 100]"))
        a = dPgF + 0.6438 - 0.001682 * vd
        c1 = a * -6.11873891971188577088 + 1.17752614766485175224
        c2 = a * 18.27315722388047447566 + -8.93204522498095698779
        c3 = a * -14.55275321129051135927 + 7.91015964461522003148
        c4 = a * 3.48385106908642905310 + -1.80321117937358499361
        k = (nd - 1) / vd
        coefficients = SVector{7,T}(nd + k * c1, k * c2, -1, k * c3, -2, k * c4, -3)
        new{T}(nd, vd, dPgF, coefficients)
    end
end

"""
    ModelGlass(nd::Real, vd::Real, dPgF::Real = 0.0)

Construct a [`ModelGlass`](@ref), promoting `nd`/`vd`/`dPgF` to a common
element type (so a `ForwardDiff` dual in any one of them carries
through). Throws an `ArgumentError` if `nd` is outside `[1.0, 2.5]` or
`vd` is outside `[10, 100]` -- values that far out are almost certainly
a parse error, not a real glass.

Arguments:
- `nd` -- refractive index at 587.5618 nm
- `vd` -- Abbe number at the d-line
- `dPgF` -- deviation of partial dispersion from the normal line
  (defaults to 0, i.e. a glass on the normal line)
"""
ModelGlass(nd::Real, vd::Real, dPgF::Real = 0.0) =
    ModelGlass{promote_type(typeof(nd), typeof(vd), typeof(dPgF))}(nd, vd, dPgF)

"""
    (g::ModelGlass)(λ::Real)

Refractive index of model glass `g` at wavelength `λ` (microns), via
[`riFormula5`](@ref). Logs a warning (once per session) if `λ` is
outside the model's 0.36-0.75 μm fit range, but still returns the
extrapolated value.
"""
function (g::ModelGlass)(λ::Real)
    (λ < 0.36 || λ > 0.75) &&
        @warn "ModelGlass evaluated at λ = $λ μm, outside its 0.36-0.75 μm fit range; index is extrapolated" maxlog = 1
    riFormula5(λ, g.coefficients)
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
through [`riFormulas`](@ref) to [`riFormula1`](@ref)-[`riFormula9`](@ref)
(`"formula 1"`-`"formula 9"`), and parsing its `"coefficients"` field.

Returns `nothing` if the record has no `"coefficients"` entry (e.g. a
tabulated-data-only record with no dispersion formula) or its formula
type isn't in [`riFormulas`](@ref) (skipped with a message, at load
time rather than on first evaluation); otherwise
returns a 1-argument function `wavelength -> index`, suitable for use
as a [`defaultGlassCatalog`](@ref) entry.

Arguments:
- `basepath` -- root of the glass-data directory tree
- `path` -- glass YAML file path, relative to `basepath`
"""
function getRefractiveIndexFunc(basepath::AbstractString, path::AbstractString)
    pth = joinpath(basepath,path)
    record = YAML.load_file(pth)
    data = record["DATA"][1]
    typeData = data["type"]
    coefentries = get(data, "coefficients", nothing)
    if coefentries == nothing
        return nothing
    end
    fn = get(riFormulas, typeData, nothing)
    if fn === nothing
        println("Skipping $pth: unsupported dispersion formula \"$typeData\"")
        return nothing
    end
    coefs = map(x->parse(Float64,x),split(coefentries))
    f(x) = fn(x, coefs)
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
