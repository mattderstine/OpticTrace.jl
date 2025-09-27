# characterization functions and traces


export sag, randomPointOnSquare, randomPointOnDisk, traceMonteCarloRays, computeRearFocalPlane, traceLoss


function randomPointOnSquare(xhalf)
    r = 2.0* rand(2) .- (1., 1.)
    a=xhalf * r
    (a[1], a[2], 0.)
end


function randomPointOnDisk(rmax)
    r = 2.0* rand(2) .- (1., 1.)
    while (norm(r) > 1.)
        r = 2.0* rand(2) .- (1., 1.)
    end
    a=rmax * r
    (a[1], a[2], 0.)
end

"""
    rays,cnt,missed= traceMonteCarloRays(radiusfunc,anglefunc, radius::Float64, θmax::Float64,
        pnts::Int64, geo::Array{AbstractSurface}; surfnum=-1)

    returns
        an array of rays at the target surface
        the number of terminated rays
        an array describing the terminated rays by surface and type

    radiusfunc is a function returns Float64 that takes one argument defining
        the radius, relative to first surface
    anglefunc is a function returning FLoat64 that takes one argument defining
        the maximum angle, relative to first surface
    radius - defining size
    θmax - defining angle
    pnts - number of rays to start
    surfnum is the number of the surface +1 to sample, -1 is last surface

"""
function traceMonteCarloRays(radiusfunc,anglefunc, radius::Float64, θmax::Float64, pnts::Int64,  geo::Array{AbstractSurface}; surfnum=-1)
    trcStatMsg=("Normal","Missed","TIR","Clipped")

    rays=Vector{Ray}(undef,pnts)
    miss=Vector{Tuple}(undef,pnts)
    badray = Ray((NaN, NaN, NaN), (NaN, NaN, NaN))
    missed = zeros(Int32,length(trcStatMsg)-1, length(geo)+1) #types of errors, length of trace one more than length of geometry
    #println("status = $(length(trcStatMsg)-1)  length = $(length(geo)+1)")
    cnt = Threads.Atomic{Int64}(0)
    Threads.@threads for i in 1:pnts
        inray = Ray(radiusfunc(radius), anglefunc(θmax))
        status, trc = traceGeometryRel(inray, geo)
        lentrc = length(trc)
        if status != 0 && (surfnum <0 || surfnum > lentrc)
            #print(trcStatMsg[status+1],"$i  ")
            rays[i] = badray  #this one is to get filtered out
            miss[i] = (status, lentrc)
            Threads.atomic_add!(cnt,1)
        else
            rays[i]= (surfnum == -1) ? trc[end].ray : trc[surfnum].ray
            miss[i]=(0,0) #done this way so can do raytrace in parallel
        end

    end

    for m in miss #combine all parallel miss data
        if m[1] !=0
            missed[m[1], m[2]]+=1
        end
    end
    #delta =  (sum(missed)) - cnt[] # cnt is an Atomic

    #println("MonteCarloRays  clipped = $(cnt[]) of $pnts")
    raysf = filter(x->!isnan(x.base[1]),rays)  #remove clipped rays
    #println("remaining rays = $(length(raysf)) == $(pnts-cnt[])")
    raysf,cnt[],missed  #cnt is an Atomic, need [] to get value
end


function computeRearFocalPlane(geo; epsilon = 0.001)
    statusb,bore = traceGeometryRel(Ray(ORIGIN, ZAXIS), geo)
    statusy,ytrace = traceGeometryRel(Ray(Point(0., epsilon, 0.), ZAXIS), geo)
    statusx,xtrace = traceGeometryRel(Ray(Point(epsilon, 0., 0.), ZAXIS), geo)

    println("computeRearFocalPlane\nBore")
    printTrcCoords(statusb, bore, geo)
    println("ytrace")
    printTrcCoords(statusy, ytrace, geo)
    println("xtrace")
    printTrcCoords(statusx, xtrace, geo)
    #check to make sure all made it through
    if statusb !=0 || statusy !=0 || statusx !=0
        return(1, NaN, NaN, ORIGIN, ORIGIN)
    end

    #compute intersection of bore and ytrace, xtrace
    #get last ray, rotate input rays to direction of bore, find when length that y, x are zero

    rayb = bore[end].ray
    rayx = xtrace[end].ray
    rayy = ytrace[end].ray

    dirb = rayb.dir
    dirx = rayx.dir
    diry = rayy.dir

    baseb = rayb.base
    basex = rayx.base
    basey = rayy.base

#=
    if dirb[3] != 1.
        println("not yet general calculations $dirb")
        return(2, NaN, NaN, ORIGIN, ORIGIN)
     end
=#

    if dirx[1] != 0.
        lenxzero = -basex[1]/dirx[1]
        lenxeps = (epsilon - basex[1])/dirx[1]
    else
        lenxzero = NaN
        lenxeps = NaN
    end

    if diry[2] != 0.
        lenyzero = -basey[2]/diry[2]
        lenyeps = (epsilon - basey[2])/diry[2]
    else
        lenyzero = NaN
        lenyeps = NaN
    end

    xplane = SVector(baseb[1], baseb[2], basex[3] + lenxzero * dirx[3])


    yplane = SVector(baseb[1], baseb[2], basey[3] + lenyzero * diry[3])

    flx = (lenxzero - lenxeps) * dirx[3]
    fly = (lenyzero - lenyeps) * diry[3]

    return (0, fly, flx, yplane, xplane)
end


"""
    traceLoss(trc, l = 1.0)
    calculate the loss of an Array{Trace}
TBW
"""
function traceLoss(trc, l = 1.0)
    for t in trc
        l *= t.pmatrix.trans[1]
    end
    return l
end

