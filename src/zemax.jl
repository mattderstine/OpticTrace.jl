#=
    read a zemax file and return the data as a geometry object
=#
export readZemax, printZemaxSurfs, zemaxsurfsToGeo, viewZemaxFile

"""
    ZemaxGeometry{N, T}

Container type intended to hold a fully-imported Zemax system: the traced
geometry, the base point/direction the geometry starts from, and the
system's wavelength/name/units metadata.

Fields:
    geo::Vector{AbstractSurface{N,T}}  - the imported optical surfaces
    basept::Point{N, T}                - global coordinate the geometry starts at
    dir::Vec{N, T}                     - propagation direction the geometry starts along
    wavelengths::Vector{T}             - wavelengths defined in the Zemax file
    name::String                       - system name, from the Zemax file's NAME field
    units::String                      - length units, from the Zemax file's UNIT field (e.g. "MM")

Not currently constructed anywhere in this file -- [`readZemax`](@ref)
returns its parsed data as a plain tuple rather than wrapping it in a
`ZemaxGeometry`.
"""
struct ZemaxGeometry{N, T}
    geo::Vector{AbstractSurface{N,T}}
    basept::Point{N, T}
    dir::Vec{N, T}
    wavelengths::Vector{T}
    name::String
    units::String
end


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
end

"""
    ZemaxSurf()

Construct a `ZemaxSurf{Float64}` with default values: zero curvature,
distance, and conic; `"DEFAULT"` material; not a stop; `parmlength`
zeroed aspheric coefficients; empty coating; `"STANDARD"` type; empty
comment.
"""
ZemaxSurf() = ZemaxSurf(0.0, 0.0, "DEFAULT", 0.0, false, 0.0, zeros(Float64, parmlength), "", "STANDARD","")


#resetZemaxSurf!(s::ZemaxSurf) = (s.curvature = 0.0; s.distance = 0.0; s.material = "AIR"; s.radius = 0.0; s.stop = false; s.conic = 0.0; s.aspherics = zeros(Float64, parmlength); s.coating = ""; s.type = "STANDARD")

"""
    readZemax(filename::String; basept = ORIGIN, dir = ZAXIS)

Read a Zemax `.zmx` sequential-lens text file and parse it into a vector
of [`ZemaxSurf`](@ref) records, one per `SURF` block (plus a leading
record for the object surface).

Recognized line keywords: `UNIT`, `NAME`, `WAVM`, `CURV`, `DISZ`, `GLAS`,
`DIAM`, `STOP`, `TYPE`, `CONI`, `PARM`, `COAT`, `COMM`. Any other
keyword is silently ignored. `basept`/`dir` are accepted but not
currently used during parsing itself -- surface positions are computed
later, by [`zemaxsurfToSurface`](@ref)/[`zemaxsurfsToGeo`](@ref).

Returns `(zsurfs, name, units, wavelengths)`:
    zsurfs::Vector{ZemaxSurf}    - the parsed surfaces, in file order
    name::String                 - system name, from the NAME field
    units::String                 - length units, from the UNIT field
    wavelengths::Vector{Float64}  - fixed-length (24-element) array indexed by Zemax wavelength number, from WAVM
"""
function readZemax(filename::String; basept = ORIGIN, dir = ZAXIS)


    geo = Vector{AbstractSurface}()


    units = "MM"
    name = "Zemax System"

    wavelengths = zeros(Float64, 24) #array of 24 wavelengths
    surfnum =0

    lines = readlines(filename)

    zsurfs = Vector{ZemaxSurf}()
    zsurf = ZemaxSurf()

    basecurrent = basept #this variable gets updated by zemaxsurfToSurface!
    dircurrent = dir #this variable will get updated by zemaxsurfToSurface! if coordinate breaks are implemented
    rinCur = rInDef() # initial refractive index, thhis variable gets updated by zemaxsurfToSurface!

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
                #surf = zemaxsurfToSurface!(zemaxsurf,basecurrent, dircurrent, rinCur)
                #push!(geo, surf)
                zsurf = ZemaxSurf()
                #resetZemaxSurf!(zsurf)
            end
            # Handle surface definitions
            # Example: geo.geo.push!(Spheroid(...))
        elseif startswith(entries[1], "UNIT")
            # Handle unit definitions
            units = entries[2]
        elseif startswith(entries[1], "NAME")
            # Handle name definitions
            name = strip(split(curline, " ", limit=2)[2], '"')
        elseif startswith(entries[1], "WAVM")
            # Handle wavelength definitions
            wavelengths[parse(Int32,entries[2])] = parse(Float64, entries[3])
        elseif startswith(entries[1], "CURV")
            # Handle curvature definitions
            zsurf.curvature = parse(Float64, entries[2])
        elseif startswith(entries[1], "DISZ")
            # Handle distance definitions
            zsurf.distance = parse(Float64, entries[2])
            basecurrent += zsurf.distance * dircurrent
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
    #surf = zemaxsurf_to_surface(zemaxsurf,basecurrent, dircurrent, rinCur)
    #push!(geo, surf)
    # Create the ZemaxGeometry object
    #zgeo = ZemaxGeometry(geo, basepnt, dir, wavelengths, name, units)
    push!(zsurfs, zsurf)
    return zsurfs, name, units, wavelengths
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
    zemaxsurfToSurface(num, rinIn::Float64, rinOut::Float64, basept::Point, dir::Vec3, s::ZemaxSurf)

Convert one parsed [`ZemaxSurf`](@ref) record into an `OptSurface`
located at `basept`, oriented along `dir`, with entering/exiting
refractive indices `rinIn`/`rinOut`. `num` is only used to build the
surface's name (`s.comm * " \$num"`).

Supports Zemax `type`s `"STANDARD"` (built as a `SurfProfileConic`) and
`"EVENASPH"` (built as a `SurfProfileEvenAsphere`, using
`s.aspherics[2:end]`); any other type throws an error.

Returns `(newbasept, rinOut, newsurf)`, where `newbasept = basept +
s.distance * dir` is the starting point for the next surface, and
`rinOut` is passed through unchanged so it can be reused as the next
surface's `rinIn`.
"""
function zemaxsurfToSurface(num, rinIn::Float64, rinOut::Float64, basept::Point, dir::Vec3,  s::ZemaxSurf)
    color = :lightgray

    if s.type == "STANDARD"
        ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(basept, dir, nothing)
        newsurf = OptSurface(s.comm * " $num",
            SurfBase(basept, dir, ydir),
            SizeLens(s.radius),
            SurfProfileConic( s.curvature, conicToϵ(s.conic)),
            DielectricT(rinIn, rinOut),
            #AmpParam(coating),
            getAmpParams(s.coating; attributesSurfaces),
            toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
            color
            )
    elseif s.type == "EVENASPH"
        ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(basept, dir, nothing)
        newsurf = OptSurface(s.comm * " $num",
            SurfBase(basept, dir, ydir),
            SizeLens(s.radius),
            SurfProfileEvenAsphere( s.curvature, conicToϵ(s.conic), s.aspherics[2:end]),
            DielectricT(rinIn, rinOut),
            #AmpParam(coating),
            getAmpParams(s.coating; attributesSurfaces),
            toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
            color
        )
    else
        error("Zemax surface type $(s.type) not implemented yet")
    end
    return basept + s.distance * dir, rinOut, newsurf
end

"""
    zemaxsurfsToGeo(zemaxsurfs, base, dir, wavelength::Float64; glassCatalog::Dict{AbstractString, Any} = defaultGlassCatalog)

Convert a vector of [`ZemaxSurf`](@ref) records (as returned by
[`readZemax`](@ref)) into a traceable geometry at the given wavelength,
by repeatedly calling [`zemaxsurfToSurface`](@ref) and threading the base
point and refractive index from one surface to the next.

`glassCatalog` maps each surface's `material` name to a function of
wavelength returning its refractive index (see `defaultGlassCatalog`).

Returns the resulting `Vector{AbstractSurface}`.
"""
function zemaxsurfsToGeo(zemaxsurfs, base, dir, wavelength::Float64; glassCatalog::Dict{AbstractString, Any} =defaultGlassCatalog)
    geo = Vector{AbstractSurface}()
    rinIn = rInDef()
    for (i,zsurf) in enumerate(zemaxsurfs)
        rinOut = glassCatalog[zsurf.material](wavelength)
        base, rinIn, surf = zemaxsurfToSurface(i,rinIn, rinOut, base, dir, zsurf)
        push!(geo, surf)
    end
    return geo
end

"""
    viewZemaxFile(filename; glassCatalog = defaultGlassCatalog)

Read, convert, print, and plot a Zemax `.zmx` file in one call: reads
`filename` with [`readZemax`](@ref), prints the parsed surfaces with
[`printZemaxSurfs`](@ref), builds a geometry with
[`zemaxsurfsToGeo`](@ref) (skipping the leading object surface,
`zgeo[2:end]`) at a fixed wavelength of 0.5, displays it with
`plotGeometry3D`, and prints it with `printGeo`.

Returns the built geometry (`Vector{AbstractSurface}`).
"""
function viewZemaxFile(filename; glassCatalog=defaultGlassCatalog)
    zgeo, name, units, wave = readZemax(filename; basept = ORIGIN, dir = ZAXIS)
    println("Zemax System Name: $name Units: $units")
    printZemaxSurfs(zgeo)
    geo = zemaxsurfsToGeo(zgeo[2:end], ORIGIN, ZAXIS, 0.5; glassCatalog)

    fig,a = plotGeometry3D(geo)
    display(fig)
    printGeo(geo)
    return geo
end
