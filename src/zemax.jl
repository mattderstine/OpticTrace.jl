#= 
    read a zemax file and return the data as a geometry object

=#
export readZemax, printZemaxSurfs, zemaxsurfsToGeo

#=

    
def zmx_to_system(data, item=None):
    s = System()
    next_pos = 0.
    s.append(Spheroid(material=air))
    for line in data.splitlines():
        e = s[-1]
        if not line.strip():
            continue
        line = line.strip().split(" ", 1)
        cmd = line[0]
        args = len(line) == 2 and line[1] or ""
        if cmd == "UNIT":
            s.scale = {
                    "MM": 1e-3,
                    "INCH": 25.4e-3,
                    "IN": 25.4e-3,
                    }[args.split()[0]]
        elif cmd == "NAME":
            s.description = args.strip("\"")
        elif cmd == "SURF":
            s.append(Spheroid(distance=next_pos, material=air))
        elif cmd == "CURV":
            e.curvature = float(args.split()[0])
        elif cmd == "DISZ":
            next_pos = float(args)
        elif cmd == "GLAS":
            args = args.split()
            name = args[0]
            try:
                e.material = Material.make(name)
            except KeyError:
                try:
                    e.material = Material.make((float(args[3]),
                                                float(args[4])))
                except Exception as e:
                    print("material not found", name, e)
        elif cmd == "DIAM":
            e.radius = float(args.split()[0])
        elif cmd == "STOP":
            e.stop = True
        elif cmd == "WAVL":
            s.wavelengths = [float(i)*1e-6 for i in args.split() if i]
        elif cmd == "COAT":
            e.coating = args.split()[0]
        elif cmd == "CONI":
            e.conic = float(args.split()[0])
        elif cmd == "PARM":
            i, j = args.split()
            i = int(i) - 1
            j = float(j)
            if i < 0:
                if j:
                    print("aspheric 0 degree not supported", cmd, args)
                continue
            if e.aspherics is None:
                e.aspherics = []
            while len(e.aspherics) <= i:
                e.aspherics.append(0.)
            e.aspherics[i] = j
        elif cmd in ("GCAT",  # glass catalog names
                     "OPDX",  # opd
                     "RAIM",  # ray aiming
                     "CONF",  # configurations
                     "ENPD", "PUPD",  # pupil
                     "EFFL",  # focal lengths
                     "VERS",  # version
                     "MODE",  # mode
                     "NOTE",  # note
                     "TYPE",  # surface type
                     "HIDE",  # surface hide
                     "MIRR",  # surface is mirror
                     "PARM",  # aspheric parameters
                     "SQAP",  # square aperture?
                     "XDAT", "YDAT",  # xy toroidal data
                     "OBNA",  # object na
                     "PKUP",  # pickup
                     "MAZH", "CLAP", "PPAR", "VPAR", "EDGE", "VCON",
                     "UDAD", "USAP", "TOLE", "PFIL", "TCED", "FNUM",
                     "TOL", "MNUM", "MOFF", "FTYP", "SDMA", "GFAC",
                     "PUSH", "PICB", "ROPD", "PWAV", "POLS", "GLRS",
                     "BLNK", "COFN", "NSCD", "GSTD", "DMFS", "ISNA",
                     "VDSZ", "ENVD", "ZVDX", "ZVDY", "ZVCX", "ZVCY",
                     "ZVAN", "XFLN", "YFLN", "VDXN", "VDYN", "VCXN",
                     "VCYN", "VANN", "FWGT", "FWGN", "WWGT", "WWGN",
                     "WAVN", "WAVM", "XFLD", "YFLD", "MNCA", "MNEA",
                     "MNCG", "MNEG", "MXCA", "MXCG", "RGLA", "TRAC",
                     "FLAP", "TCMM", "FLOA", "PMAG", "TOTR", "SLAB",
                     "POPS", "COMM", "PZUP", "LANG", "FIMP",
                     ):
            pass
        else:
            print(cmd, "not handled", args)
            continue
    return s






=#

struct ZemaxGeometry{N, T}
    geo::Vector{AbstractSurface{N,T}}
    basept::Point{N, T}
    dir::Vec{N, T}
    wavelengths::Vector{T}
    name::String
    units::String
end


const     parmlength = 20
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

ZemaxSurf() = ZemaxSurf(0.0, 0.0, "DEFAULT", 0.0, false, 0.0, zeros(Float64, parmlength), "", "STANDARD","")


#resetZemaxSurf!(s::ZemaxSurf) = (s.curvature = 0.0; s.distance = 0.0; s.material = "AIR"; s.radius = 0.0; s.stop = false; s.conic = 0.0; s.aspherics = zeros(Float64, parmlength); s.coating = ""; s.type = "STANDARD")

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
            println("SURF line: $curline")
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

function printZemaxSurfs(zsurfs::Vector{ZemaxSurf})
    for (i, s) in enumerate(zsurfs)
        println("$i  $(s.curvature) $(s.distance) $(s.material) $(s.radius) $(s.stop) $(s.conic) 
        $(s.coating) $(s.type)")
        println("  Aspherics: ", s.aspherics)
    end
end

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

function zemaxsurfsToGeo(zemaxsurfs, base, dir, wavelength::Float64, glassCatalog::Dict{AbstractString, Any} )
    geo = Vector{AbstractSurface}()
    rinIn = rInDef()
    for (i,zsurf) in enumerate(zemaxsurfs)
        rinOut = glassCatalog[zsurf.material](wavelength)
        base, rinIn, surf = zemaxsurfToSurface(i,rinIn, rinOut, base, dir, zsurf)
        push!(geo, surf)
    end
    return geo
end

