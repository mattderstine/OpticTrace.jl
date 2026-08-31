#=
    read a zemax file and return the data as a geometry object
=#
export readZemax, printZemaxSurfs, zemaxsurfsToGeo, viewZemaxFile
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

#=
    .zar archive reading -- a .zar bundles a .zmx file with its
    supporting data (e.g. a glass catalog) using a small custom binary
    framing (not a zip container, despite appearances). Ported from the
    Python reference kept in docs/zemax_reference.md.
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

Ported from the Python reference kept in `docs/zemax_reference.md`
(itself adapted from
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

Ported from the Python reference kept in `docs/zemax_reference.md`.
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
    catalogs (see docs/zemax_reference.md).
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
