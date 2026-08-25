#printing.jl


export printTrcCoords, trcAndPrintRay, trcAndPrintRayRel, printTrcLen
export opdRel, trcAndPrintLengthsRel, trcAndPrintLengths, tracenumFromName, surfaceFromName
export printSurfNames, printMissed, printGeo, printSurface, printTrcStatus
export numsurfFromName


"""
    tracenumFromName(surfview, geo)
        surfview -string with user defined name of the surface
        geo - geometry

    returns the index into an Array{Trace} corresponding to geo for the first surface with that name.
"""
tracenumFromName(surfview, geo) = numsurfFromName(surfview, geo)+1



"""
    numsurfFromName(surfview, geo)

        surfview -string with user defined name of the surface
        geo - geometry

    returns the surface number. tracenumFromName adds 1 to this to index properly for Trace results
"""
function numsurfFromName(surfview, geo)
    if (surfview == "end")
        surfnum = length(geo)
    else
        surfnum = findfirst(x->x.surfname==surfview, geo)
        if isnothing(surfnum)
            println("Surface: $surfview not found. Default to end")
            surfnum = length(geo)
        end
    end
    surfnum
end

"""
    surfaceFromName(name, geo)
    name -string with user defined name of the surface
    geo - geometry

    returns the surface corresponding to the name
"""
surfaceFromName(name, geo)= geo[tracenumFromName(name,geo)-1]

"""
printTrcCoords(status, trc, geo; format="normal")
    print a trace of ray on geo

"""
function printTrcCoords(status, trc, geo; format="normal")

    println("\n--- Trace ---")
    printTrcStatus(status)
#   println("length of geo = $(length(geo))")
    for i in eachindex(trc)
        a = trc[i]
    #    println("$(a.ray.base)      $(a.ray.dir)")
        name = i==1 ? "Start" : geo[i-1].surfname
        b = a.ray.base
        d = a.ray.dir
        if format == "normal"
            @printf("%24s (%10.4g, %10.4g, %10.4g)  (%10.4g, %10.4g, %10.4g)\n",name, b[1], b[2], b[3], d[1], d[2], d[3])
        else
            println("$name  $b  $d")
        end
    end
end

"""
trcAndPrintRay(ray::Ray, geo)
    trace and print a trace of ray on geo
    ray is in absolute coordinates
"""
function trcAndPrintRay(ray::Ray, geo)
    status, trc = traceGeometry(ray, geo)
    printTrcCoords(status, trc, geo)
    trc
end

"""
trcAndPrintRayRel(ray::Ray, geo)
    trace and print a trace of ray on geo
    ray is in local coordinates
    returns the trace
"""
function trcAndPrintRayRel(ray::Ray, geo)
    status, trc = traceGeometryRel(ray, geo)
    printTrcCoords(status, trc, geo)
    trc
end

"""
printTrcLen(status, trc, geo)
    print the lengths of a trace, trc, on geo

    returns the totaldelta, total OPL (optical path length) and total reduced distance
"""
function printTrcLen(status, trc, geo)
    trcStatMsg=("Normal","Missed","TIR","Clipped")
    println("\n--- Trace Lengths---")
#   println("length of geo = $(length(geo))")
    lastone = length(trc)
    println("             Surf            Length       Index       OPL       Rdc Dist")
    totaldelta = 0.
    totalOPL = 0.
    totalRD = 0.
    for i in 1:lastone
        a = trc[i]
    #    println("$(a.ray.base)      $(a.ray.dir)")
        name = i==1 ? "Start" : geo[i-1].surfname
        if status != 0 && i == lastone
            println(trcStatMsg[status+1])
        end
        delta = a.delta
        index = a.nIn
        opl = delta *index
        rd = delta/index
        totaldelta += delta
        totalOPL += opl
        totalRD += rd
        @printf("%24s   %10.4g %10.4g %10.4g %10.4g\n",name, delta, index, opl, rd)
    end
    @printf("%24s   %10.4g            %10.4g %10.4g\n","Total", totaldelta, totalOPL, totalRD)
    [totaldelta, totalOPL, totalRD]

end


"""
    trcAndPrintLengthsRel(ray::Ray, geo)

Trace `ray` (in local coordinates, relative to `geo`'s first surface --
see `traceGeometryRel`) and print its per-surface lengths/OPL/reduced
distance via `printTrcLen`. Returns whatever `printTrcLen` returns
(`[totaldelta, totalOPL, totalRD]`). See `trcAndPrintLengths` below for
the absolute-coordinates sibling.
"""
function trcAndPrintLengthsRel(ray::Ray, geo)
    status, trc = traceGeometryRel(ray, geo)
    printTrcLen(status, trc, geo)
end

"""
    trcAndPrintLengths(ray::Ray, geo)

Trace `ray` (in absolute/global coordinates -- see `traceGeometry`) and
print its per-surface lengths/OPL/reduced distance via `printTrcLen`.
Returns whatever `printTrcLen` returns (`[totaldelta, totalOPL,
totalRD]`). See `trcAndPrintLengthsRel` above for the local-coordinates
sibling.
"""
function trcAndPrintLengths(ray::Ray, geo)
    status, trc = traceGeometry(ray, geo)
    printTrcLen(status, trc, geo)
end


"""
    printSurfNames(geo; fulldir=:false)

Print a numbered, one-line-per-surface listing of `geo`: each surface's
name, base point, and local z direction, plus its local y direction too
if `fulldir` is truthy.

The default, `fulldir = :false`, looks like a `Symbol` but isn't one:
`false` is a reserved literal token, not a valid identifier, so Julia's
`:` quoting operator applied to it just evaluates to the plain boolean
`false` (`typeof(:false) === Bool`). Equivalent to writing
`fulldir = false`, just unusual style.
"""
function printSurfNames(geo; fulldir = :false)
    for (i, surf) in enumerate(geo)
        base = surf.base.base
        dir = surf.base.dir
        ydir = surf.base.ydir
        strbasedir = @sprintf("(%10.5f, %10.5f, %10.5f)  (%10.5f, %10.5f, %10.5f)",base[1], base[2], base[3], dir[1], dir[2], dir[3] )
        if fulldir
            strydir = @sprintf("  (%10.5f, %10.5f, %10.5f)",ydir[1], ydir[2], ydir[3] )
        else
            strydir =""
        end
        println("$i\t$(surf.surfname) $(" "^(20-length(surf.surfname))) $strbasedir $strydir")
    end
end

"""
printMissed(m, geo) goes with montecarlo raytraces to display the clipped rays

"""
function printMissed(m, geo)
    i = 2
    for surf in geo
        if sum(m[:,i]) > 0
            println("   $(surf.surfname) $(" "^(20-length(surf.surfname))) $(m[1, i])  $(m[2, i])  $(m[3, i])")
        end
        i += 1
    end
end

"""
    printSurface(surf)

Shorthand for `printSurface(1, surf)` -- print `surf` labeled as
surface number 1. See `printSurface(n, surf::OptSurface)`/
`printSurface(n, surf::ModelSurface)` below for what actually gets
printed.
"""
printSurface(surf) = printSurface(1, surf)

"""
    printSurface(n, surf::OptSurface)

Print one line for `surf`, labeled `n`: its name, base point, aperture,
profile, and `mod` (refract/reflect/diffuse behavior). See
`printSurface(n, surf::ModelSurface)` below for the sibling method
(same, minus `mod`, since `ModelSurface` has none).
"""
function printSurface(n, surf::OptSurface)
    base = surf.base.base
    strbase = @sprintf("(%10.5f, %10.5f, %10.5f)",base[1], base[2], base[3] )
    println("$n  $(surf.surfname) $(" "^(20-length(surf.surfname)))  $(strbase)\n$(" "^20)$(surf.aperture)  $(surf.profile)   $(surf.mod)")
end
#=
struct ModelSurface <: AbstractSurface  #use the data to overload GemoetryBasics
    surfname::String
    base::SurfBase
    aperture::AbstractSize
    profile::AbstractSurfProfile
    refIndex::Float64 #needed for OPD calculations
    toGlobalCoord::AffineMap
    toLocalCoord::AffineMap
    toGlobalDir::LinearMap
    toLocalDir::LinearMap
end
=#
"""
    printSurface(n, surf::ModelSurface)

Print one line for `surf`, labeled `n`: its name, base point, aperture,
and profile (no `mod` line, unlike `printSurface(n, surf::OptSurface)`
above -- `ModelSurface` has no `mod` field).
"""
function printSurface(n, surf::ModelSurface)
    base = surf.base.base
    strbase = @sprintf("(%10.5f, %10.5f, %10.5f)",base[1], base[2], base[3] )
    println("$n  $(surf.surfname) $(" "^(20-length(surf.surfname)))  $strbase\n$(" "^20)$(surf.aperture)  $(surf.profile)")
end

"""
    printGeo(geometry)

Print every surface in `geometry`, numbered in order, via
`printSurface` (dispatching per-surface to the `OptSurface`/
`ModelSurface` method above as appropriate).
"""
function printGeo(geometry)
    for (n, geo) in enumerate(geometry)
        printSurface(n, geo)
    end
end

"""
    printTrcStatus(status; flagNormal=false, clipmsg=true)

Print the status message for a trace result (`"Normal"`, `"Missed"`,
`"TIR"`, or `"Clipped"`, indexed by `status`), unless `clipmsg` is
false. By default (`flagNormal=false`), a normal (`status==0`) result
prints nothing; pass `flagNormal=true` to print `"Normal"` too. Despite
the name, `clipmsg` gates every status message, not just the "Clipped"
one.
"""
function printTrcStatus(status; flagNormal = false, clipmsg=true)
    trcStatMsg=("Normal","Missed","TIR","Clipped")
    if clipmsg &&(flagNormal || status != 0 )
        println(trcStatMsg[status+1])
    end
end
