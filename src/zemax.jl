#=
    read a zemax file and return the data as a geometry object
=#
export readZemax, printZemaxSurfs, zemaxsurfsToGeo, viewZemaxFile
export ZemaxHeader, readZemaxSystem
export ZarEntry, lzwDecompress, readZemaxArchive, listZemaxArchive, extractZemaxArchive
export ZmfEntry, zmfDeobfuscate, readZmfCatalog, listZmfCatalog, extractZmfCatalog

const     parmlength = 20

"""
    ZemaxSurf{T}

Mutable intermediate representation of one Zemax `SURF` block, as parsed
from a `.zmx` file by [`readZemax`](@ref) before being converted into an
`OptSurface` by [`zemaxsurfToSurface`](@ref).

Fields:
    curvature::T          - surface curvature, from CURV
    distance::T            - thickness to the next surface, from DISZ
    material::String      - glass/material name, from GLAS (default "DEFAULT")
    radius::T              - surface (semi-)diameter, from DIAM
    stop::Bool             - whether this surface is the aperture stop, from STOP
    conic::T               - conic constant, from CONI
    aspherics::Vector{T}   - aspheric coefficients, indexed by Zemax PARM number
    coating::String        - coating name, from COAT
    type::String           - Zemax surface type, e.g. "STANDARD" or "EVENASPH", from TYPE
    comm::String           - surface comment, from COMM
    extraData::Vector{T}   - Zemax "Extra Data" coefficients, indexed by XDAT number
                             (e.g. TYPE XPOLYNOM's normalization radius + polynomial
                             terms) -- unlike `aspherics`, grows on demand rather than
                             being preallocated to a fixed length, since XDAT can run
                             longer than `parmlength`; empty when the surface has no
                             XDAT lines at all (the common case)
"""
mutable struct ZemaxSurf{T}
    curvature::T
    distance::T
    material::String
    radius::T
    stop::Bool
    conic::T
    aspherics::Vector{T}
    coating::String
    type::String
    comm::String
    extraData::Vector{T}
end

"""
    ZemaxSurf(curvature, distance, material, radius, stop, conic, aspherics, coating, type, comm)

Convenience constructor omitting `extraData` (defaults to an empty
vector -- the common case, since only a few Zemax surface types like
`TYPE XPOLYNOM` use `XDAT` lines at all).
"""
ZemaxSurf(curvature, distance, material, radius, stop, conic, aspherics, coating, type, comm) =
    ZemaxSurf(curvature, distance, material, radius, stop, conic, aspherics, coating, type, comm, eltype(aspherics)[])

"""
    ZemaxSurf()

Construct a `ZemaxSurf{Float64}` with default values: zero curvature,
distance, and conic; `"DEFAULT"` material; not a stop; `parmlength`
zeroed aspheric coefficients; empty coating; `"STANDARD"` type; empty
comment; empty extra data.
"""
ZemaxSurf() = ZemaxSurf(0.0, 0.0, "DEFAULT", 0.0, false, 0.0, zeros(Float64, parmlength), "", "STANDARD","")


#resetZemaxSurf!(s::ZemaxSurf) = (s.curvature = 0.0; s.distance = 0.0; s.material = "AIR"; s.radius = 0.0; s.stop = false; s.conic = 0.0; s.aspherics = zeros(Float64, parmlength); s.coating = ""; s.type = "STANDARD")

"""
    ZemaxHeader{T}

System-level metadata parsed from a Zemax `.zmx` file's header block
(everything before its first `SURF` line), as returned by
[`readZemax`](@ref) alongside the per-surface [`ZemaxSurf`](@ref)
vector. Absorbed into an [`OpticalSystem`](@ref) by
[`readZemaxSystem`](@ref).

Fields:
    name::String                   - system name, from NAME
    units::String                  - the *source* file's own length unit, from UNIT, kept verbatim as
                                       provenance -- by the time `readZemax` returns, every length-dimensioned
                                       field on this header/its accompanying `ZemaxSurf`s (`apertureValue`
                                       when ENPD, `fields` when height-typed, etc.) has already been converted
                                       to `LENGTH_UNIT` (mm) by `convertZemaxUnitsToMM!`, regardless of what
                                       this field says
    wavelengths::Vector{T}         - fixed-length (24-element) array indexed by Zemax wavelength number, from
                                       WAVM -- always in `WAVELENGTH_UNIT` (μm); unlike every other numeric
                                       field here, never affected by `units`/`UNIT`, since Zemax's own WAVM
                                       convention is unit-independent
    primaryWavelengthIndex::Int    - index into `wavelengths`, from PWAV (defaults to 1 if PWAV is absent)
    apertureType::String           - which aperture spec is active: "ENPD", "OBNA", "FNUM", or "FLOA" (no value)
    apertureValue::T                - the value for whichever `apertureType` is active (`NaN` for "FLOA");
                                       in `LENGTH_UNIT` (mm) when `apertureType == "ENPD"`, a dimensionless
                                       ratio for "OBNA"/"FNUM"
    fieldType::Int                  - field-type code, from FTYP's first field (angle/object height/image
                                       height/...) -- see `zemaxFieldTypeIsHeight`
    fields::Vector{Point2{T}}      - design field points, XFLN/YFLN zipped pairwise into (x,y); in
                                       `LENGTH_UNIT` (mm) when `zemaxFieldTypeIsHeight(fieldType)`, degrees
                                       otherwise
    fieldWeight::Vector{T}          - per-field weight, from FWGN, index-aligned with `fields` -- always
                                       dimensionless, never unit-converted
    glassCatalogs::Vector{String}  - glass catalog names, from GCAT
    mode::String                   - "SEQ" or "NSC", from MODE
    notes::String                  - freeform system notes, reconstructed from NOTE lines
"""
mutable struct ZemaxHeader{T<:Real}
    name::String
    units::String
    wavelengths::Vector{T}
    primaryWavelengthIndex::Int
    apertureType::String
    apertureValue::T
    fieldType::Int
    fields::Vector{Point2{T}}
    fieldWeight::Vector{T}
    glassCatalogs::Vector{String}
    mode::String
    notes::String
end

"""
    ZemaxHeader()

Construct a `ZemaxHeader{Float64}` with default values: `"Zemax
System"` name, `"MM"` units, 24 zeroed wavelengths, primary wavelength
index 1, no aperture spec, field type 0, no fields, no glass catalogs,
`"SEQ"` mode, and empty notes.
"""
ZemaxHeader() = ZemaxHeader("Zemax System", "MM", zeros(Float64, 24), 1,
    "", NaN, 0, Point2{Float64}[], Float64[], String[], "SEQ", "")

"""
    readZemax(filename::String)

Read a Zemax `.zmx` sequential-lens text file and parse it into a vector
of [`ZemaxSurf`](@ref) records, one per `SURF` block (plus a leading
record for the object surface), and a [`ZemaxHeader`](@ref) of
system-level metadata.

Recognized `SURF`-scoped line keywords: `CURV`, `DISZ`, `GLAS`, `DIAM`,
`STOP`, `TYPE`, `CONI`, `PARM`, `XDAT`, `COAT`, `COMM`. Recognized header
keywords: `UNIT`, `NAME`, `WAVM`, `PWAV`, `ENPD`/`OBNA`/`FNUM`/`FLOA`,
`FTYP`, `XFLN`/`YFLN`/`FWGN`, `GCAT`, `MODE`, `NOTE`. Any other keyword
is silently ignored. Surface positions/orientations are computed later,
by [`zemaxsurfToSurface`](@ref)/[`zemaxsurfsToGeo`](@ref).

Before returning, every length-dimensioned field of `zsurfs`/`header` is
converted from the file's own `UNIT` to [`LENGTH_UNIT`](@ref) (mm) via
[`convertZemaxUnitsToMM!`](@ref) -- see its docstring for exactly which
fields that covers. `header.units` itself keeps recording the file's
original, pre-conversion unit as provenance.

Returns `(zsurfs, header)`:
    zsurfs::Vector{ZemaxSurf}   - the parsed surfaces, in file order (zsurfs[1] is the object surface);
                                   length-dimensioned fields already converted to `LENGTH_UNIT` (mm)
    header::ZemaxHeader          - system-level metadata; length-dimensioned fields already converted to
                                   `LENGTH_UNIT` (mm), see `units`'s own field doc
"""
function readZemax(filename::String)

    header = ZemaxHeader()
    surfnum =0

    lines = readlines(filename)

    zsurfs = Vector{ZemaxSurf}()
    zsurf = ZemaxSurf()

    xfln = Float64[] # accumulated from XFLN, zipped with yfln into header.fields once parsing is done
    yfln = Float64[] # accumulated from YFLN

    # Parse the file to extract the necessary data
    for l in lines
        curline =  replace(string(strip(l,['\n','\r',' ', '\t', '\xff', '\xfe', '\0']) ), "\x00" => "") # Clean line endings and whitespace
        entries = split(curline)

        #=
        if length(entries) == 0
            println("Blank line: $curline")
            continue
        end
        for e in entries
            print(e, ", ")
        end
        println(" ")
        =#
        if length(entries) == 0
            continue
        end

        if startswith(entries[1], "SURF")
            #println("SURF line: $curline")
            surfnum = parse(Int, entries[2])
            if surfnum > 0
                # Finalize the previous surface before starting a new one
                # Example: geo.geo[end].aspherics = copy(parm)
                push!(zsurfs, zsurf)
                zsurf = ZemaxSurf()
                #resetZemaxSurf!(zsurf)
            end
            # Handle surface definitions
            # Example: geo.geo.push!(Spheroid(...))
        elseif startswith(entries[1], "UNIT")
            # Handle unit definitions
            header.units = entries[2]
        elseif startswith(entries[1], "NAME")
            # Handle name definitions -- a bare `NAME` line (no quoted
            # name following) is valid Zemax output for an unnamed
            # system; leave header.name empty in that case instead of
            # throwing (see TODO.md Bug #13)
            parts = split(curline, " ", limit=2)
            header.name = length(parts) > 1 ? strip(parts[2], '"') : ""
        elseif startswith(entries[1], "WAVM")
            # Handle wavelength definitions
            header.wavelengths[parse(Int32,entries[2])] = parse(Float64, entries[3])
        elseif startswith(entries[1], "PWAV")
            # Handle primary wavelength index
            header.primaryWavelengthIndex = parse(Int, entries[2])
        elseif startswith(entries[1], "ENPD") || startswith(entries[1], "OBNA") || startswith(entries[1], "FNUM")
            # Handle aperture definitions -- value is the 2nd token; a
            # trailing flag digit (e.g. afocal/telecentric) is ignored
            header.apertureType = string(entries[1])
            header.apertureValue = parse(Float64, entries[2])
        elseif startswith(entries[1], "FLOA")
            # Handle float-by-stop aperture -- bare flag, no value
            header.apertureType = "FLOA"
            header.apertureValue = NaN
        elseif startswith(entries[1], "FTYP")
            # Handle field-type definitions -- only the field-type code
            # (1st token) is parsed; other FTYP tokens (telecentricity,
            # field count, afocal-image-space flag) aren't decoded yet,
            # see TODO.md
            header.fieldType = parse(Int, entries[2])
        elseif startswith(entries[1], "XFLN")
            # Handle field X positions -- zipped with YFLN into header.fields after parsing
            append!(xfln, parse.(Float64, entries[2:end]))
        elseif startswith(entries[1], "YFLN")
            # Handle field Y positions -- zipped with XFLN into header.fields after parsing
            append!(yfln, parse.(Float64, entries[2:end]))
        elseif startswith(entries[1], "FWGN")
            # Handle field weights
            append!(header.fieldWeight, parse.(Float64, entries[2:end]))
        elseif startswith(entries[1], "GCAT")
            # Handle glass catalog list
            header.glassCatalogs = string.(entries[2:end])
        elseif startswith(entries[1], "MODE")
            # Handle system mode (SEQ/NSC)
            header.mode = entries[2]
        elseif startswith(entries[1], "NOTE")
            # Handle system notes -- only "NOTE 0 <text>" lines carry
            # text; other note-line codes (e.g. "NOTE 4") are formatting
            # markers, not text, and are skipped
            if length(entries) >= 2 && entries[2] == "0"
                noteText = length(entries) > 2 ? join(entries[3:end], " ") : ""
                header.notes = isempty(header.notes) ? noteText : header.notes * "\n" * noteText
            end
        elseif startswith(entries[1], "CURV")
            # Handle curvature definitions
            zsurf.curvature = parse(Float64, entries[2])
        elseif startswith(entries[1], "DISZ")
            # Handle distance definitions
            zsurf.distance = parse(Float64, entries[2])
        elseif startswith(entries[1], "GLAS")
            # Handle glass/material definitions
            zsurf.material = entries[2]
        elseif startswith(entries[1], "DIAM")
            # Handle diameter definitions
            zsurf.radius = parse(Float64, entries[2])
        elseif startswith(entries[1], "STOP")
            # Handle stop surface definitions
            zsurf.stop = true
        elseif startswith(entries[1], "TYPE")
            # set type of surface
            zsurf.type = string(entries[2])
        elseif startswith(entries[1], "CONI")
            # Handle conic definitions
            zsurf.conic = parse(Float64, entries[2])
        elseif startswith(entries[1], "PARM")
            # Handle aspheric parameter definitions
            idx = parse(Int, entries[2])  # 1-based index
            val = parse(Float64, entries[3])
            if 1 <= idx <= parmlength
                zsurf.aspherics[idx] = val
            else
                @warn "PARM index $idx out of bounds"
            end
        elseif startswith(entries[1], "XDAT")
            # Handle "Extra Data" definitions (e.g. TYPE XPOLYNOM's
            # normalization radius + polynomial terms) -- unlike PARM,
            # not preallocated to a fixed length: grow extraData on
            # demand so an index beyond its current length doesn't need
            # a bounds warning the way an out-of-range PARM index does
            idx = parse(Int, entries[2])  # 1-based index
            val = parse(Float64, entries[3])
            while length(zsurf.extraData) < idx
                push!(zsurf.extraData, 0.0)
            end
            zsurf.extraData[idx] = val
        elseif startswith(entries[1], "COAT")
            # Handle coating definitions
            zsurf.coating = entries[2]
        elseif startswith(entries[1], "COMM")
            # Handle comment definitions
            zsurf.comm = join(entries[2:end], " ")
        else
            # Handle other commands or ignore
        end
    end
    push!(zsurfs, zsurf)
    header.fields = Point2.(xfln, yfln)
    convertZemaxUnitsToMM!(zsurfs, header)
    return zsurfs, header
end

"""
    zemaxUnitToMM(units::String) -> Float64

Millimeters per one `units` (a Zemax `UNIT` token: `"MM"`, `"CM"`,
`"IN"`, or `"M"`). Throws an `ArgumentError` for anything else, rather
than silently treating an unrecognized token as millimeters (Zemax's own
supported `UNIT` values).
"""
function zemaxUnitToMM(units::String)::Float64
    units == "MM" && return 1.0
    units == "CM" && return 10.0
    units == "IN" && return 25.4
    units == "M" && return 1000.0
    throw(ArgumentError("zemaxUnitToMM: unrecognized Zemax UNIT \"$units\" " *
        "(expected one of \"MM\", \"CM\", \"IN\", \"M\")"))
end

"""
    zemaxFieldTypeIsHeight(fieldType::Int) -> Bool

Whether Zemax `FTYP` field-type code `fieldType` represents a
length-dimensioned field position (object height, paraxial image
height, or real image height -- codes `1`, `2`, `3`) rather than an
angle in degrees (codes `0` "Angle" and `4` "Theodolite Angle"). Used by
[`convertZemaxUnitsToMM!`](@ref) to decide whether `header.fields`
(`XFLN`/`YFLN`) needs unit conversion. Throws an `ArgumentError` for any
other code.
"""
function zemaxFieldTypeIsHeight(fieldType::Int)::Bool
    fieldType in (1, 2, 3) && return true
    fieldType in (0, 4) && return false
    throw(ArgumentError("zemaxFieldTypeIsHeight: unrecognized Zemax FTYP field-type code " *
        "$fieldType (expected 0-4)"))
end

"""
    convertZemaxUnitsToMM!(zsurfs, header::ZemaxHeader)

Convert every length-dimensioned field of `zsurfs`/`header` in place
from the source file's own unit (`header.units`, a Zemax `UNIT` token)
to [`LENGTH_UNIT`](@ref) (mm), via [`zemaxUnitToMM`](@ref). Called once
by [`readZemax`](@ref), right before it returns, so every consumer of
`readZemax`'s result -- directly, or through [`readZemaxSystem`](@ref) --
gets already-mm-canonical geometry, rather than relying on each caller
to remember a separate conversion step. A no-op when `header.units` is
already `"MM"`.

`header.units` itself is left unchanged: it keeps recording the
*source* file's original unit as provenance (see [`OpticalSystem`](@ref)'s
`units` field docs), even though every numeric field it's attached to is
now in mm regardless of what `header.units` says.

Which fields get scaled, and by what power of the mm-per-unit factor,
depends on `s.type` for the per-surface `aspherics`/`extraData` fields
-- Zemax reuses those arrays for different physical quantities
(different length dimensions, or no length dimension at all, e.g. a
tilt angle or a dimensionless flag) depending on surface type, mirroring
[`zemaxsurfToProfileAperture`](@ref)'s own per-type dispatch. Confirmed
against each type's `sag` method in `src/tracing.jl`:

- every type: `curvature` (÷ factor, it's `1/length`), `distance`,
  `radius` (× factor).
- `"EVENASPH"`: `aspherics[i]` for `i = 2:end` multiplies `r^(2i)` (see
  `sag(x, y, ::SurfProfileEvenAsphere)`), so each scales by
  `factor^(1-2i)`.
- `"ODDASPHE"`: `aspherics[i]` for `i = 1:8` multiplies `r^i` (see
  `sag(x, y, ::SurfProfileOddAsphere)`), so each scales by
  `factor^(1-i)`.
- `"TOROIDAL"`: `aspherics[1]` is `Rx`, a radius (× factor).
- `"COORDBRK"`: `aspherics[1]`/`aspherics[2]` are `dx`/`dy` decenters (×
  factor each); `aspherics[3:5]` are tilt angles in degrees and
  `aspherics[6]` is the order flag -- neither scaled.
- `"TILTSURF"`: `aspherics[1]`/`aspherics[2]` are tilt angles in
  degrees -- not scaled.
- `"PARAXIAL"`: `aspherics[1]` is a focal length (× factor);
  `aspherics[2]` (an OPD-mode flag) is not scaled.
- `"XPOLYNOM"`: `extraData[1]` is the normalization radius (× factor);
  `extraData[3:end]` are polynomial term coefficients that multiply an
  already-normalized, dimensionless `(x/normRadius)^m (y/normRadius)^n`
  (see `sag(x, y, ::SurfProfileXYPoly)`), so -- unlike `EVENASPH`/
  `ODDASPHE` -- they scale uniformly by the factor regardless of term;
  `extraData[2]` (an unidentified control flag) is not scaled.
- `"STANDARD"`: no extra length-dimensioned fields.

Also converts `header.apertureValue` (× factor), but only when
`header.apertureType == "ENPD"` (entrance pupil diameter, a length) --
`"OBNA"`/`"FNUM"` are dimensionless ratios and `"FLOA"` has no value, so
neither is scaled. `header.fields` (`XFLN`/`YFLN`) is converted (× factor)
only when [`zemaxFieldTypeIsHeight`](@ref)`(header.fieldType)`;
`header.fieldWeight` is never scaled (always dimensionless).
`header.wavelengths` is never scaled -- Zemax's `WAVM` values are always
in [`WAVELENGTH_UNIT`](@ref) regardless of `header.units`.
"""
function convertZemaxUnitsToMM!(zsurfs, header::ZemaxHeader)
    factor = zemaxUnitToMM(header.units)
    factor == 1.0 && return nothing

    for s in zsurfs
        s.curvature /= factor
        s.distance *= factor
        s.radius *= factor
        if s.type == "EVENASPH"
            for i in 2:length(s.aspherics)
                s.aspherics[i] *= factor^(1 - 2i)
            end
        elseif s.type == "ODDASPHE"
            for i in 1:min(8, length(s.aspherics))
                s.aspherics[i] *= factor^(1 - i)
            end
        elseif s.type == "TOROIDAL"
            length(s.aspherics) >= 1 && (s.aspherics[1] *= factor)
        elseif s.type == "COORDBRK"
            length(s.aspherics) >= 1 && (s.aspherics[1] *= factor)
            length(s.aspherics) >= 2 && (s.aspherics[2] *= factor)
        elseif s.type == "PARAXIAL"
            length(s.aspherics) >= 1 && (s.aspherics[1] *= factor)
        elseif s.type == "XPOLYNOM"
            length(s.extraData) >= 1 && (s.extraData[1] *= factor)
            for i in 3:length(s.extraData)
                s.extraData[i] *= factor
            end
        end
        # "STANDARD"/"TILTSURF": no length-dimensioned aspherics/extraData fields.
    end

    header.apertureType == "ENPD" && (header.apertureValue *= factor)
    if !isempty(header.fields) && zemaxFieldTypeIsHeight(header.fieldType)
        header.fields = factor .* header.fields
    end

    return nothing
end

"""
    printZemaxSurfs(zsurfs::Vector{ZemaxSurf})

Print a one-line summary (curvature, distance, material, radius, stop,
conic, coating, type) for each surface in `zsurfs` to stdout, followed
by its aspheric coefficient vector. Useful for inspecting the result of
[`readZemax`](@ref).
"""
function printZemaxSurfs(zsurfs::Vector{ZemaxSurf})
    println("surface curvature distance material radius stop conic coating type")
    for (i, s) in enumerate(zsurfs)
        println("$i  $(s.curvature) $(s.distance) $(s.material) $(s.radius) $(s.stop) $(s.conic)
        $(s.coating) $(s.type)")
        println("  Aspherics: ", s.aspherics)
    end
end

"""
    zemaxsurfToSurface(num, rinIn::Float64, rinOut::Float64, basept::Point, dir::Vec3,
        ydir::Union{Vec3,Nothing}, s::ZemaxSurf)

Convert one parsed [`ZemaxSurf`](@ref) record into an `OptSurface`
located at `basept`, oriented along `dir`/`ydir` (`ydir = nothing`
resolves to the same default guess `updateCoordChange` uses elsewhere),
with entering/exiting refractive indices `rinIn`/`rinOut`. `num` is only
used to build the surface's name (`s.comm * " \$num"`).

The profile/aperture come from [`zemaxsurfToProfileAperture`](@ref)
(see its docstring for which Zemax `type`s are supported). `TYPE
TILTSURF` is the one exception to "surface built at `(basept, dir,
ydir)` as given": its own tilt (`PARM 1`/`PARM 2`, tilt about x/y) is
applied first via [`zemaxCoordBreakFrame`](@ref), and the surface is
built at the resulting tilted frame instead.

`s.material == "MIRROR"` (Zemax's reserved name for a reflective
surface, e.g. `GLAS MIRROR`) builds a `MirrorR` bend instead of the
usual `DielectricT` -- see [`zemaxsurfsToGeo`](@ref)'s docstring for
the caller-side half of this (skipping the glass catalog for it, and
why `rinIn == rinOut` here for a mirror). `s.type == "PARAXIAL"`
builds a `ParaxialLensT` bend instead (checked before the `MIRROR`
case), with its focal length from `s.aspherics[1]` (Zemax `PARM 1`;
`PARM 2`, an "OPD Mode" flag, doesn't affect real ray tracing and is
left undecoded), and a distinct `:cyan3` display color (instead of the
usual `:lightgray`) so an ideal lens is visually identifiable in a
plot -- it has no real curvature to distinguish it otherwise.

Returns `(newbasept, newdir, newydir, rinOut, newsurf)`: `newdir`/
`newydir` are `dir`/`ydir` unchanged for every type except `TILTSURF`
(where they carry the tilt forward to whatever comes next);
`newbasept = basept + s.distance * newdir` is the starting point for
the next surface; `rinOut` is passed through unchanged so it can be
reused as the next surface's `rinIn`.
"""
function zemaxsurfToSurface(num, rinIn::Float64, rinOut::Float64, basept::Point, dir::Vec3,
        ydir::Union{Vec3,Nothing}, s::ZemaxSurf)
    color = s.type == "PARAXIAL" ? :cyan3 : :lightgray
    profile, aperture = zemaxsurfToProfileAperture(s)

    if s.type == "TILTSURF"
        basept, dir, ydir = zemaxCoordBreakFrame(basept, dir, ydir,
            zero(s.curvature), zero(s.curvature), s.aspherics[1], s.aspherics[2], zero(s.curvature), 0)
    end

    myydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
            updateCoordChange(basept, dir, ydir)
    bend = if s.type == "PARAXIAL"
        ParaxialLensT(s.aspherics[1], rinIn, rinOut)
    elseif s.material == "MIRROR"
        MirrorR(rinIn, rinOut)
    else
        DielectricT(rinIn, rinOut)
    end
    newsurf = OptSurface(s.comm * " $num",
        SurfBase(basept, dir, myydir),
        aperture,
        profile,
        bend,
        #AmpParam(coating),
        getAmpParams(s.coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
    return basept + s.distance * dir, dir, myydir, rinOut, newsurf
end

"""
    zemaxsurfToProfileAperture(s::ZemaxSurf) -> (profile, aperture)

Convert one parsed [`ZemaxSurf`](@ref) record's own geometry fields
(`type`/`curvature`/`conic`/`aspherics`/`radius`) into a
`(profile::AbstractSurfProfile, aperture::SizeLens)` pair -- the shared
step factored out of [`zemaxsurfToSurface`](@ref) (which wraps it into a
refracting `OptSurface`) and [`zemaxObjectToModelSurface`](@ref) (which
wraps it into a non-refracting `ModelSurface` for the object surface),
so both stay in sync about which Zemax `type`s are supported.

Supports Zemax `type`s `"STANDARD"` (a `SurfProfileConic`), `"EVENASPH"`
(a `SurfProfileEvenAsphere`, using `s.aspherics[2:end]`), `"TOROIDAL"`
(a `SurfProfileToroid` -- `s.curvature`/`s.conic` describe the base y-z
curve, and `s.aspherics[1]` is Zemax's own `PARM 1` "Radius of
Rotation" `Rx`, converted to the x-sweep curvature `1/Rx`; `Rx == 0`
maps to a curvature of `0`, matching Zemax's own convention that a `0`
radius of rotation means no x sweep at all rather than a literal
zero-radius one -- confirmed against real `TOROIDAL` samples under the
local Zemax install's `Samples` directory), `"TILTSURF"` (a
`SurfProfileConic`, identical to `"STANDARD"` -- Zemax's `TILTSURF` is
an otherwise-ordinary conic surface, with the tilt itself handled
separately by [`zemaxsurfToSurface`](@ref) via
[`zemaxCoordBreakFrame`](@ref)), `"ODDASPHE"` (a `SurfProfileOddAsphere`
-- Zemax's own `PARM i` multiplies `rⁱ` directly for `i = 1..8`, so
`a = s.aspherics[1:8]`, confirmed against a real `ODDASPHE` sample: an
axicon whose only nonzero term, `PARM 1`, produces the linear cone
profile an axicon is defined by), and `"XPOLYNOM"` (a
`SurfProfileXYPoly` -- built from `s.extraData` rather than
`s.aspherics`, since Zemax stores this type's coefficients via `XDAT`
lines: `s.extraData[1]` is the normalization radius, `s.extraData[2]`
an unidentified Zemax control flag not interpreted here, and
`s.extraData[3:end]` the polynomial term coefficients in Zemax's own
bivariate term order -- see [`xyPolyTermPowers`](@ref)), and
`"PARAXIAL"` (a `ParaxialProfile` -- an ideal thin lens has no real
sag; its focal length, from `s.aspherics[1]`, is handled separately by
[`zemaxsurfToSurface`](@ref), which builds a `ParaxialLensT` bend
instead of the usual `DielectricT`); any other type throws an error.
"""
function zemaxsurfToProfileAperture(s::ZemaxSurf)
    aperture = SizeLens(s.radius)
    if s.type == "STANDARD" || s.type == "TILTSURF"
        return SurfProfileConic(s.curvature, conicToϵ(s.conic)), aperture
    elseif s.type == "EVENASPH"
        return SurfProfileEvenAsphere(s.curvature, conicToϵ(s.conic), s.aspherics[2:end]), aperture
    elseif s.type == "TOROIDAL"
        rx = s.aspherics[1]
        curvX = rx == 0 ? zero(rx) : 1 / rx
        return SurfProfileToroid(s.curvature, conicToϵ(s.conic), curvX), aperture
    elseif s.type == "ODDASPHE"
        return SurfProfileOddAsphere(s.curvature, conicToϵ(s.conic), s.aspherics[1:8]), aperture
    elseif s.type == "XPOLYNOM"
        normRadius = s.extraData[1]
        return SurfProfileXYPoly(s.curvature, conicToϵ(s.conic), normRadius, s.extraData[3:end]), aperture
    elseif s.type == "PARAXIAL"
        return ParaxialProfile(zero(s.curvature)), aperture
    else
        error("Zemax surface type $(s.type) not implemented yet")
    end
end

"""
    zemaxCoordBreakFrame(basept::Point, dir::Vec3, ydir::Union{Vec3,Nothing},
        dx, dy, tiltXdeg, tiltYdeg, tiltZdeg, orderFlag) -> (newbasept, newdir, newydir)

Apply a Zemax coordinate-break decenter+tilt to the current local frame
`(basept, dir, ydir)` (`ydir = nothing` resolves to the same default
guess `updateCoordChange`/`findPerpenMap` use elsewhere), following
Zemax's own documented `COORDBRK` "Order" convention:

- `orderFlag == 0` (the common case): decenter `dx`/`dy` first, along
  the *current* (pre-tilt) local x/y axes, then tilt -- intrinsically,
  about the local axes as they evolve -- in X, then Y, then Z order.
- `orderFlag != 0`: reversed -- tilt first (intrinsically Z, then Y,
  then X), then decenter `dx`/`dy` along the *resulting* (post-tilt)
  local x/y axes. This is Zemax's own documented mechanism for undoing
  an earlier coordinate break with a second one.

Tilt angles are in degrees (as raw Zemax `PARM` values are); a positive
angle follows the standard right-hand rule about the corresponding
(current, evolving) local axis, built via [`rotationX`](@ref)/
[`rotationY`](@ref)/[`rotationZ`](@ref) -- confirmed empirically against
the local Zemax install's "Fold Mirror Using Coordinate Breaks.ZMX"
sample (a 45° `COORDBRK` tilt about local x immediately ahead of a
mirror correctly folds an on-axis beam by exactly 90°, matching a real
fold mirror's physical behavior, when traced by hand through this
formula and `modFunc(::MirrorR)`'s reflection law).

Used by [`zemaxsurfsToGeo`](@ref) for `TYPE COORDBRK` surfaces
(decenter+tilt only, no surface built) and by
[`zemaxsurfToSurface`](@ref) for `TYPE TILTSURF` surfaces (tilt only,
`dx=dy=tiltZdeg=0`, `orderFlag=0`, then a real surface *is* built at the
resulting frame).
"""
function zemaxCoordBreakFrame(basept::Point, dir::Vec3, ydir::Union{Vec3,Nothing},
        dx::T, dy::T, tiltXdeg::T, tiltYdeg::T, tiltZdeg::T, orderFlag) where T<:Real
    curydir, toGlobalDir = findPerpenMap(dir, ydir)
    curxdir = normalize(cross(curydir, dir))

    Rx = rotationX(deg2rad(tiltXdeg))
    Ry = rotationY(deg2rad(tiltYdeg))
    Rz = rotationZ(deg2rad(tiltZdeg))
    Rlocal = orderFlag == 0 ? Rx * Ry * Rz : Rz * Ry * Rx

    Mnew = toGlobalDir.linear * Rlocal
    newdir = normalize(Vec3(Mnew * ZAXIS))
    newydir = normalize(Vec3(Mnew * YAXIS))

    if orderFlag == 0
        newbasept = basept + dx * curxdir + dy * curydir
    else
        newxdir = normalize(Vec3(Mnew * XAXIS))
        newbasept = basept + dx * newxdir + dy * newydir
    end

    return newbasept, newdir, newydir
end

"""
    zemaxObjectToModelSurface(s::ZemaxSurf, basept::Point, dir::Vec3; rin = rInDef(), color = :lightgray)

Convert the object surface (Zemax `SURF 0`, `s`) into a non-refracting
`ModelSurface` carrying its own profile/aperture, via
[`zemaxsurfToProfileAperture`](@ref) -- most files' object surface is a
flat, inert placeholder (`CURV 0`), but Zemax allows real geometry on
it too (e.g. a curved source, or a reverse-traced eye model's retina),
so this preserves it rather than discarding `s` entirely the way
[`readZemaxSystem`](@ref) does for the traceable chain (`geo` never
includes the object surface).

Position: `basept - s.distance * dir` when `s.distance` is finite
(undoing the same `+distance*dir` step used to place every other
surface, since the object sits *upstream* of surface 1 by its own
distance); when `s.distance` is infinite, positioned at `basept`
itself, an arbitrary, documented anchor -- no finite global position is
meaningful for an infinite-conjugate object, only the surface's local
profile is. If `s.type == "TILTSURF"` (a real, if unusual, Zemax file --
e.g. a tilted-object test target), that tilt (`PARM 1`/`PARM 2`, no
decenter) is applied in place at that position, via
[`zemaxCoordBreakFrame`](@ref), the same way `TILTSURF` is handled for
every other surface in [`zemaxsurfToSurface`](@ref).
"""
function zemaxObjectToModelSurface(s::ZemaxSurf, basept::Point, dir::Vec3; rin = rInDef(), color = :lightgray)
    profile, aperture = zemaxsurfToProfileAperture(s)
    position = isinf(s.distance) ? basept : basept - s.distance * dir
    tiltydir = nothing

    if s.type == "TILTSURF"
        position, dir, tiltydir = zemaxCoordBreakFrame(position, dir, nothing,
            zero(s.curvature), zero(s.curvature), s.aspherics[1], s.aspherics[2], zero(s.curvature), 0)
    end

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
            updateCoordChange(position, dir, tiltydir)
    ModelSurface(isempty(s.comm) ? "Object" : s.comm,
        SurfBase(position, dir, ydir),
        aperture,
        profile,
        rin,
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end

"""
    zemaxsurfsToGeo(zemaxsurfs, base, dir, wavelength::Float64;
        ydir::Union{Vec3,Nothing} = nothing, glassCatalog::Dict{AbstractString, Any} = defaultGlassCatalog)

Convert a vector of *real* [`ZemaxSurf`](@ref) records into a traceable
geometry at the given wavelength, by repeatedly calling
[`zemaxsurfToSurface`](@ref) and threading the running local frame
(`base`/`dir`/`ydir`) and refractive index from one surface to the
next.

`TYPE COORDBRK` surfaces are handled specially: they build no
`OptSurface` at all, only update the running frame via
[`zemaxCoordBreakFrame`](@ref) (decenter + tilt, per Zemax's own
`COORDBRK` "Order" convention) and then advance `base` by their own
`DISZ` along the *resulting* (post-tilt) `dir` -- matching Zemax's rule
that a coordinate break's own thickness is always applied last,
regardless of the order its decenter/tilt were applied in. Every other
surface goes through `zemaxsurfToSurface` as before (which itself
handles `TYPE TILTSURF`'s own tilt).

`zemaxsurfs` must **not** include the object surface (`readZemax`'s
`zsurfs[1]`) -- go through [`readZemaxSystem`](@ref), which handles
that split, rather than calling this directly on a full `readZemax`
result. Any `isinf(zsurf.distance)` encountered here is treated as an
anomaly (a malformed file, or a real surface incorrectly passed as if
it were the object) and raises a clear error rather than silently
propagating `Inf`/`NaN` into later surfaces' coordinates -- the object
surface's own (possibly infinite) distance is meant to be handled
before this function ever sees `zemaxsurfs`, not by this loop.

`glassCatalog` maps each surface's `material` name to a function of
wavelength returning its refractive index (see `defaultGlassCatalog`) --
never consulted for a `COORDBRK` surface (which has no material), nor
for a surface whose material is the reserved name `"MIRROR"` (Zemax's
convention for a reflective surface, e.g. `GLAS MIRROR`; not a real
glass, so never present in a catalog): reflection doesn't change the
medium, so such a surface's `rinOut` is just its `rinIn` carried
through unchanged, and [`zemaxsurfToSurface`](@ref) builds it with a
`MirrorR` bend instead of the usual `DielectricT`.

Only a mirror bracketed by real `TYPE COORDBRK` tilts (as in every
known real sample) is handled correctly -- those already physically
rotate `dir` via rotation matrices, so a coordinate break's own
(possibly negative) `DISZ` afterward is just ordinary signed 3D
displacement along the new `dir`, no special-casing needed. A **bare**
mirror (no coordinate break) relies on a different, implicit Zemax
convention instead (the local frame doesn't itself rotate; only the
*sign* of later `DISZ` values encodes the fold) that this function has
no equivalent for and does not attempt to model. This only matters if
an OpticTrace-to-Zemax *export* path is ever written -- it would need
to always emit a mirror bracketed by a matching pair of coordinate
breaks, never a bare reflective surface.

Returns the resulting `Vector{AbstractSurface}`.
"""
function zemaxsurfsToGeo(zemaxsurfs, base, dir, wavelength::Float64;
        ydir::Union{Vec3,Nothing} = nothing, glassCatalog::Dict{AbstractString, Any} =defaultGlassCatalog)
    geo = Vector{AbstractSurface}()
    rinIn = rInDef()
    for (i,zsurf) in enumerate(zemaxsurfs)
        isinf(zsurf.distance) && error("zemaxsurfsToGeo: surface $i has an infinite distance; only the " *
            "object surface (readZemax's zsurfs[1]) may have DISZ INFINITY -- did you pass the full " *
            "readZemax result instead of using readZemaxSystem?")
        if zsurf.type == "COORDBRK"
            base, dir, ydir = zemaxCoordBreakFrame(base, dir, ydir,
                zsurf.aspherics[1], zsurf.aspherics[2], zsurf.aspherics[3], zsurf.aspherics[4],
                zsurf.aspherics[5], zsurf.aspherics[6])
            base = base + zsurf.distance * dir
        else
            rinOut = zsurf.material == "MIRROR" ? rinIn : glassCatalog[zsurf.material](wavelength)
            base, dir, ydir, rinIn, surf = zemaxsurfToSurface(i, rinIn, rinOut, base, dir, ydir, zsurf)
            push!(geo, surf)
        end
    end
    return geo
end

"""
    readZemaxSystem(filename; basept = ORIGIN, dir = ZAXIS, wavelength = 0.5, glassCatalog = defaultGlassCatalog)

Read a Zemax `.zmx` file into a materialized [`OpticalSystem`](@ref) --
the canonical entry point for the Zemax-import pipeline, replacing the
hand-glued `readZemax`+slice+`zemaxsurfsToGeo` sequence
[`viewZemaxFile`](@ref) used to do inline.

Calls [`readZemax`](@ref), then splits its result: `zsurfs[1]` (the
object surface) never becomes part of the traceable `geo` -- its own
geometry is preserved separately as `OpticalSystem.objectSurface` (via
[`zemaxObjectToModelSurface`](@ref)), and its distance (finite or
`Inf`) becomes `OpticalSystem.objectDistance`/`objectAtInfinity`. The
remaining surfaces (`zsurfs[2:end]`) are converted via
[`zemaxsurfsToGeo`](@ref) at `wavelength`.

Raises an error if the file's `MODE` is `"NSC"` (non-sequential) --
this package's tracer is sequential-only end-to-end, so a genuinely
non-sequential file can't be imported meaningfully (see `TODO.md`).

Returns an `OpticalSystem{Float64}`.
"""
function readZemaxSystem(filename; basept = ORIGIN, dir = ZAXIS, wavelength = 0.5, glassCatalog = defaultGlassCatalog)
    zsurfs, header = readZemax(filename)
    header.mode == "NSC" && error("readZemaxSystem: $filename is a non-sequential-mode (MODE NSC) file; " *
        "not supported -- this package's tracer is sequential-only, see TODO.md")

    objSurf = zsurfs[1]
    geo = zemaxsurfsToGeo(zsurfs[2:end], basept, dir, wavelength; glassCatalog)
    objectSurface = zemaxObjectToModelSurface(objSurf, basept, dir)

    OpticalSystem(geo, objectSurface, header.name, header.units, header.wavelengths,
        header.primaryWavelengthIndex, objSurf.distance, isinf(objSurf.distance),
        header.apertureType, header.apertureValue, header.fieldType, header.fields,
        header.fieldWeight, header.glassCatalogs, header.mode, header.notes)
end

"""
    viewZemaxFile(filename; glassCatalog = defaultGlassCatalog)

Read, convert, print, and plot a Zemax `.zmx` file in one call: reads
`filename` with [`readZemaxSystem`](@ref), prints the parsed surfaces
with [`printZemaxSurfs`](@ref), displays the resulting geometry with
`plotGeometry3D`, and prints it with `printGeo`.

Returns the built [`OpticalSystem`](@ref).
"""
function viewZemaxFile(filename; glassCatalog=defaultGlassCatalog)
    zgeo, header = readZemax(filename)
    println("Zemax System Name: $(header.name) Units: $(header.units)")
    printZemaxSurfs(zgeo)
    sys = readZemaxSystem(filename; basept = ORIGIN, dir = ZAXIS, wavelength = 0.5, glassCatalog)

    fig,a = plotGeometry3D(sys.geo)
    display(fig)
    printGeo(sys.geo)
    return sys
end

#=
    .zar archive reading -- a .zar bundles a .zmx file with its
    supporting data (e.g. a glass catalog) using a small custom binary
    framing (not a zip container, despite appearances). Ported from a
    Python reference implementation.
=#

"""
    ZarEntry

One packed entry from a Zemax `.zar` archive, as returned by
[`readZemaxArchive`](@ref).

Fields:
    name::String         - the entry's file name (with any `.LZW` compression suffix already stripped)
    data::Vector{UInt8}  - the entry's (decompressed, if it was compressed) raw bytes
"""
struct ZarEntry
    name::String
    data::Vector{UInt8}
end

"""
    readBits(data::Vector{UInt8}, bitIndex::Int, n::Int) -> Int

Read `n` bits starting at zero-based bit offset `bitIndex` from `data`,
treating the byte vector as one contiguous, most-significant-bit-first
bitstream. Used by [`lzwDecompress`](@ref) to read variable-width LZW
codewords.
"""
function readBits(data::Vector{UInt8}, bitIndex::Int, n::Int)::Int
    value = 0
    for i in 0:n-1
        globalBit = bitIndex + i
        byteIdx = globalBit ÷ 8
        bitOffset = 7 - (globalBit % 8)
        bit = (data[byteIdx + 1] >> bitOffset) & 0x1
        value = (value << 1) | Int(bit)
    end
    return value
end

"""
    lzwDecompress(compressed::Vector{UInt8}) -> Vector{UInt8}

Decompress a byte vector using the variable-width LZW algorithm Zemax
`.zar` archives use for `.LZW`-suffixed entries: a dictionary seeded
with the 256 single-byte sequences, codewords starting at 9 bits and
growing by 1 bit each time the dictionary size crosses a power of two,
and the classic LZW "KwKwK" fallback when a codeword isn't yet in the
dictionary.

Ported from a Python reference implementation (itself adapted from
https://gist.github.com/BertrandBordage/611a915e034c47aa5d38911fc0bc7df9),
reading bits directly from `compressed` via [`readBits`](@ref) rather
than materializing a giant binary string.
"""
function lzwDecompress(compressed::Vector{UInt8})::Vector{UInt8}
    totalBits = length(compressed) * 8
    codeWordLength = 8
    words = Vector{Vector{UInt8}}(undef, 256)
    for i in 0:255
        words[i+1] = UInt8[i]
    end

    bitIndex = 0
    previousWord = UInt8[]
    decompressed = UInt8[]

    while true
        if 2^codeWordLength <= length(words)
            codeWordLength += 1
        end
        if bitIndex + codeWordLength > totalBits
            break
        end
        code = readBits(compressed, bitIndex, codeWordLength)
        bitIndex += codeWordLength

        latestWord = code < length(words) ? words[code + 1] : vcat(previousWord, previousWord[1:1])
        append!(decompressed, latestWord)
        if length(previousWord) > 0
            push!(words, vcat(previousWord, latestWord[1:1]))
        end
        previousWord = latestWord
    end

    return decompressed
end

"""
    stripExtension(path::AbstractString, ext::AbstractString)

Return `path` with a trailing `ext` removed (case-insensitively), or
`path` unchanged if it doesn't end with `ext`. Used by
[`extractZemaxArchive`](@ref)/[`extractZmfCatalog`](@ref) to compute a
default extraction directory.
"""
function stripExtension(path::AbstractString, ext::AbstractString)
    return endswith(lowercase(path), lowercase(ext)) ? path[1:end-length(ext)] : path
end

"""
    readZemaxArchive(filename::String) -> Vector{ZarEntry}

Read a Zemax `.zar` archive (a `.zmx` file bundled with its supporting
data, e.g. a glass catalog) and return its packed entries.

Each archive entry has a small per-entry header (a 2-byte version tag,
a packed byte count, and a name) immediately followed by that many
bytes of payload; entries whose name ends in `.LZW` are LZW-compressed
(decompressed here via [`lzwDecompress`](@ref), with the `.LZW` suffix
stripped from the returned name). Two header layouts exist -- "earlier"
(tagged `0xEA`) and "latest" (tagged `0xEC`), identified by the first
byte of each entry's version tag -- with different header lengths, size
field widths, and name encodings (UTF-8 for "earlier", UTF-16LE for
"latest"; both layouts have been exercised against real sample files).
An unrecognized version byte throws an error rather than guessing.
Numeric fields are read assuming a little-endian host, as is standard
on essentially all current hardware.

Ported from a Python reference implementation.
"""
function readZemaxArchive(filename::String)::Vector{ZarEntry}
    entries = ZarEntry[]
    open(filename) do io
        while !eof(io)
            versionBytes = read(io, 2)
            length(versionBytes) < 2 && break
            versionByte = versionBytes[1]

            if versionByte == 0xEA
                headerLength = 0x14C - 2
                sizeRange = (0xC - 2 + 1):(0x10 - 2)
                nameOffset = (0x20 - 2) + 1
                latest = false
            elseif versionByte == 0xEC
                headerLength = 0x288 - 2
                sizeRange = (0x10 - 2 + 1):(0x18 - 2)
                nameOffset = (0x30 - 2) + 1
                latest = true
            else
                error("Unknown ZAR header version byte 0x$(string(versionByte, base=16))")
            end

            header = read(io, headerLength)

            if latest
                packedSize = Int(reinterpret(UInt64, header[sizeRange])[1])
                nameBytes = header[nameOffset:end]
                nameU16 = Vector{UInt16}(reinterpret(UInt16, nameBytes))
                name = transcode(String, nameU16)
                nullIdx = findfirst('\0', name)
                name = nullIdx === nothing ? name : name[1:prevind(name, nullIdx)]
            else
                packedSize = Int(reinterpret(UInt32, header[sizeRange])[1])
                nameBytes = header[nameOffset:end]
                nullPos = findfirst(==(0x00), nameBytes)
                nameBytes = nullPos === nothing ? nameBytes : nameBytes[1:nullPos-1]
                name = String(nameBytes)
            end

            payload = read(io, packedSize)
            if length(name) >= 4 && uppercase(name[end-3:end]) == ".LZW"
                payload = lzwDecompress(payload)
                name = name[1:end-4]
            end

            push!(entries, ZarEntry(name, payload))
        end
    end
    return entries
end

"""
    listZemaxArchive(filename::String) -> Vector{String}

List the entry names packed inside a Zemax `.zar` archive, without
writing anything to disk. Thin wrapper over [`readZemaxArchive`](@ref).
"""
listZemaxArchive(filename::String) = [e.name for e in readZemaxArchive(filename)]

"""
    extractZemaxArchive(filename::String, names::AbstractVector{<:AbstractString}; outputPath=nothing) -> Vector{String}

Extract only the named entries from a Zemax `.zar` archive `filename`,
writing each to `outputPath` (default: `filename` with its `.zar`
extension stripped, created via `mkpath` if it doesn't already exist).
Returns the paths written, in the order `names` was given. Throws an
error if a requested name isn't present in the archive.
"""
function extractZemaxArchive(filename::String, names::AbstractVector{<:AbstractString}; outputPath=nothing)
    entries = readZemaxArchive(filename)
    outDir = outputPath === nothing ? stripExtension(filename, ".zar") : outputPath
    mkpath(outDir)
    paths = String[]
    for name in names
        idx = findfirst(e -> e.name == name, entries)
        idx === nothing && error("Entry \"$name\" not found in archive $filename")
        path = joinpath(outDir, name)
        write(path, entries[idx].data)
        push!(paths, path)
    end
    return paths
end

"""
    extractZemaxArchive(filename::String; outputPath=nothing) -> Vector{String}

Extract every entry from a Zemax `.zar` archive `filename` into a new
subdirectory named after `filename` with its `.zar` extension stripped
(created via `mkpath` if it doesn't already exist), itself created under
`outputPath` (default: `filename`'s own directory) -- e.g.
`extractZemaxArchive("/a/b/lens.zar")` writes into `/a/b/lens/`, and
`extractZemaxArchive("/a/b/lens.zar"; outputPath="/x")` writes into
`/x/lens/`. Returns the paths written.
"""
function extractZemaxArchive(filename::String; outputPath=nothing)
    entries = readZemaxArchive(filename)
    base = outputPath === nothing ? dirname(filename) : outputPath
    outDir = joinpath(base, basename(stripExtension(filename, ".zar")))
    mkpath(outDir)
    paths = String[]
    for e in entries
        path = joinpath(outDir, e.name)
        write(path, e.data)
        push!(paths, path)
    end
    return paths
end

#=
    .zmf catalog reading -- a .zmf bundles a whole vendor's stock-lens
    catalog into one binary file: a 4-byte version header followed by
    one fixed-size metadata record per lens, each immediately followed
    by that lens's XOR-obfuscated .zmx-style description text. Format
    and obfuscation formula are an unofficial, community-reverse-
    engineered convention (not documented by Ansys/Zemax); ported from
    rayopt's zemax.py (https://github.com/quartiq/rayopt/blob/master/rayopt/zemax.py),
    which is itself the same convention several independent tools rely
    on. Empirically validated in development against real vendor
    catalogs.
=#

"""
    ZmfEntry

One decoded lens record from a Zemax `.zmf` lens-catalog file, as
returned by [`readZmfCatalog`](@ref).

Fields:
    name::String         - catalog part name/number
    efl::Float64          - effective focal length, from the record header
    enp::Float64          - entrance pupil diameter, from the record header
    elements::Int         - element count, from the record header
    data::Vector{UInt8}   - decoded lens description: `.zmx`-format text as raw bytes, directly writable to a file [`readZemax`](@ref) can read
"""
struct ZmfEntry
    name::String
    efl::Float64
    enp::Float64
    elements::Int
    data::Vector{UInt8}
end

"""
    zmfDeobfuscate(data::Vector{UInt8}, efl::Float64, enp::Float64) -> Vector{UInt8}

Reverse the XOR keystream Zemax `.zmf` catalog files use to obfuscate
each lens's description text. The keystream is derived per output byte
from a fixed trigonometric formula seeded by that lens's own `efl`/
`enp`, then reduced to a single byte by formatting the intermediate
value in `%.8e` scientific notation and taking 3 of its digit
characters -- reproduced here exactly, including that string
round-trip, so the keystream matches the reference implementation
bit-for-bit (the digit positions taken depend on whether the formatted
string has a leading `-`, which is intentional and load-bearing, not a
bug).

The same function both obfuscates and deobfuscates, since XOR is its
own inverse. Ported from `rayopt`'s `zmf_obfuscate`
(https://github.com/quartiq/rayopt/blob/master/rayopt/zemax.py).
"""
function zmfDeobfuscate(data::Vector{UInt8}, efl::Float64, enp::Float64)::Vector{UInt8}
    iv = cos(6 * efl + 3 * enp)
    iv = cos(655 * (π / 180) * iv) + iv
    out = similar(data)
    for p in 0:length(data)-1
        k = 13.2 * (iv + sin(17 * (p + 3))) * (p + 1)
        digits = (@sprintf("%.8e", k))[5:7]
        keyByte = UInt8(parse(Int, digits) % 256)
        out[p+1] = data[p+1] ⊻ keyByte
    end
    return out
end

"""
    readZmfCatalog(filename::String) -> Vector{ZmfEntry}

Read a Zemax `.zmf` lens-catalog file (as vendors like Edmund/Thorlabs
publish stock lens families) and decode each lens record.

The file is a 4-byte little-endian format-version header (only `1001`
is documented/observed; any other value throws an error) followed by
one fixed 144-byte record per lens (100-byte name + 7 `UInt32` fields +
2 `Float64` fields, all little-endian) and that record's obfuscated
description text (length given by the record, deobfuscated via
[`zmfDeobfuscate`](@ref)). Numeric fields are read assuming a
little-endian host, as is standard on essentially all current hardware.

Ported from `rayopt`'s `zmf_read`
(https://github.com/quartiq/rayopt/blob/master/rayopt/zemax.py).
"""
function readZmfCatalog(filename::String)::Vector{ZmfEntry}
    entries = ZmfEntry[]
    open(filename) do io
        version = read(io, UInt32)
        version != 1001 && error("Unsupported ZMF catalog version $version (only 1001 is supported)")
        while !eof(io)
            nameBytes = read(io, 100)
            length(nameBytes) < 100 && break
            read(io, UInt32)              # per-lens format version -- not currently surfaced on ZmfEntry
            elements = read(io, UInt32)
            read(io, UInt32)              # shape code, indexes "?EBPM" -- not currently surfaced
            read(io, UInt32)              # aspheric flag -- not currently surfaced
            read(io, UInt32)              # grin flag -- not currently surfaced
            read(io, UInt32)              # toroidal flag -- not currently surfaced
            descLen = read(io, UInt32)
            efl = read(io, Float64)
            enp = read(io, Float64)

            nullPos = findfirst(==(0x00), nameBytes)
            trimmed = nullPos === nothing ? nameBytes : nameBytes[1:nullPos-1]
            name = String(trimmed)

            description = read(io, descLen)
            length(description) != descLen && error("Truncated .zmf description for lens \"$name\"")
            data = zmfDeobfuscate(description, efl, enp)

            push!(entries, ZmfEntry(name, efl, enp, Int(elements), data))
        end
    end
    return entries
end

"""
    listZmfCatalog(filename::String) -> Vector{String}

List the lens names in a Zemax `.zmf` catalog, without decoding any
lens's obfuscated description text -- cheap even for large, many-lens
vendor catalogs, since only the always-plaintext name field of each
fixed-size record needs reading; the rest of each record (including the
description) is skipped over rather than read.
"""
function listZmfCatalog(filename::String)::Vector{String}
    names = String[]
    open(filename) do io
        version = read(io, UInt32)
        version != 1001 && error("Unsupported ZMF catalog version $version (only 1001 is supported)")
        while !eof(io)
            nameBytes = read(io, 100)
            length(nameBytes) < 100 && break
            skip(io, 4 * 6)   # recVersion, elements, shape, aspheric, grin, toroidal
            descLen = read(io, UInt32)
            skip(io, 8 * 2)   # efl, enp
            skip(io, descLen) # obfuscated description text -- not decoded for a listing

            nullPos = findfirst(==(0x00), nameBytes)
            trimmed = nullPos === nothing ? nameBytes : nameBytes[1:nullPos-1]
            push!(names, String(trimmed))
        end
    end
    return names
end

"""
    extractZmfCatalog(filename::String, names::AbstractVector{<:AbstractString}; outputPath=nothing) -> Vector{String}

Decode only the named lenses from a Zemax `.zmf` catalog `filename` and
write each as `<name>.zmx` under `outputPath` (default: `filename` with
its `.zmf`/`.ZMF` extension stripped, created via `mkpath` if it
doesn't already exist) -- each extracted file is directly readable by
[`readZemax`](@ref), since a decoded lens description is exactly
`.zmx`-format text. Returns the paths written, in the order `names` was
given. Throws an error if a requested name isn't present in the
catalog.
"""
function extractZmfCatalog(filename::String, names::AbstractVector{<:AbstractString}; outputPath=nothing)
    entries = readZmfCatalog(filename)
    outDir = outputPath === nothing ? stripExtension(filename, ".zmf") : outputPath
    mkpath(outDir)
    paths = String[]
    for name in names
        idx = findfirst(e -> e.name == name, entries)
        idx === nothing && error("Lens \"$name\" not found in catalog $filename")
        path = joinpath(outDir, name * ".zmx")
        write(path, entries[idx].data)
        push!(paths, path)
    end
    return paths
end

"""
    extractZmfCatalog(filename::String; outputPath=nothing) -> Vector{String}

Decode every lens in a Zemax `.zmf` catalog `filename` and write each as
`<name>.zmx` into a new subdirectory named after `filename` with its
`.zmf`/`.ZMF` extension stripped (created via `mkpath` if it doesn't
already exist), itself created under `outputPath` (default: `filename`'s
own directory) -- e.g. `extractZmfCatalog("/a/b/cat.zmf")` writes into
`/a/b/cat/`, and `extractZmfCatalog("/a/b/cat.zmf"; outputPath="/x")`
writes into `/x/cat/`. Returns the paths written.
"""
function extractZmfCatalog(filename::String; outputPath=nothing)
    entries = readZmfCatalog(filename)
    base = outputPath === nothing ? dirname(filename) : outputPath
    outDir = joinpath(base, basename(stripExtension(filename, ".zmf")))
    mkpath(outDir)
    paths = String[]
    for e in entries
        path = joinpath(outDir, e.name * ".zmx")
        write(path, e.data)
        push!(paths, path)
    end
    return paths
end
