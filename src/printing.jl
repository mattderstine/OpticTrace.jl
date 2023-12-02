#printing.jl


export printTrcCoords, trcAndPrintRay, trcAndPrintRayRel, printTrcLen
export opdRel, trcAndPrintLengthsRel, trcAndPrintLengths, surfnumFromName
export printSurfNames, printMissed, printGeo, printSurface, printTrcStatus


"""
    surfFromName(surfview, geo)
    surfview -string with user defined name of the surface
    geo - geometry

    returns the surface number +1. This indexes properly for Trace results
"""
function surfnumFromName(surfview, geo)
    if (surfview == "end")
        surfnum = length(geo)+1
    else
        surfnum = findfirst(x->x.surfname==surfview, geo)
        if isnothing(surfnum)
            println("Surface: $surfview not found. Default to end")
            surfnum = length(geo)+1
        else
            surfnum = surfnum +1
        end
    end
    surfnum
end


"""
printTrcCoords(ray::Ray, geo; color=:blue, clipmsg=false)
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
trcandPrintRay(ray::Ray, geo; color=:blue, clipmsg=false)
    trace and print a trace of ray on geo
    ray is in absolute coordinates    
"""
function trcAndPrintRay(ray::Ray, geo)
    status, trc = traceGeometry(ray, geo)
    printTrcCoords(status, trc, geo)
    trc
end

"""
trcandPrintRayRel(ray::Ray, geo; color=:blue, clipmsg=false)
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
    
    returns the totaldelta, total OPD and total reduced distance
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


function trcAndPrintLengthsRel(ray::Ray, geo)
    status, trc = traceGeometryRel(ray, geo)
    printTrcLen(status, trc, geo)
end

function trcAndPrintLengths(ray::Ray, geo)
    status, trc = traceGeometry(ray, geo)
    printTrcLen(status, trc, geo)
end


function printSurfNames(geo)
    for (i, surf) in enumerate(geo)
        base = surf.base.base
        dir = surf.base.dir
        strbasedir = @sprintf("(%10.5f, %10.5f, %10.5f)  (%10.5f, %10.5f, %10.5f)",base[1], base[2], base[3], dir[1], dir[2], dir[3] )
        println("$i\t$(surf.surfname) $(" "^(20-length(surf.surfname))) $strbasedir")
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
function printSurface(n, surf::ModelSurface)
    base = surf.base.base
    strbase = @sprintf("(%10.5f, %10.5f, %10.5f)",base[1], base[2], base[3] )
    println("$n  $(surf.surfname) $(" "^(20-length(surf.surfname)))  $strbase\n$(" "^20)$(surf.aperture)  $(surf.profile)")
end  

function printGeo(geometry)
    for (n, geo) in enumerate(geometry)
        printSurface(n, geo)
    end
end

function printTrcStatus(status; flagNormal = false, clipmsg=true)
    trcStatMsg=("Normal","Missed","TIR","Clipped")
    if clipmsg &&(flagNormal || status != 0 )
        println(trcStatMsg[status+1])
    end
end
