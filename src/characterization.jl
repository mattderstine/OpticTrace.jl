# characterization functions and traces


export sag, randomPointOnSquare, randomPointOnDisk, traceMonteCarloRays, computeRearFocalPlane, traceLoss
export findRFP, surfClosestApproach, distClosestApproach

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
    surfClosestApproach(ray1, ray2; semiDiam = 5.0, ri=refIndexDefault, surfname="closest"))
    Given two rays, find the point of closest approach
    return a surface along ray1 at that point and
    the distance along ray1 to that point
    semiDiam is the semi-diameter of the returned surface

"""
function surfClosestApproach(ray1::Ray, ray2::Ray; semiDiam = 5.0, ri=refIndexDefault, surfname="closest")
    t1,t2 = distClosestApproach(ray1, ray2)
 

    closep1 =ray1.base .+ t1 * ray1.dir
    #println("closest point on ray1 = $closep1  t1 = $t1")
    refSurf = referencePlane(surfname, closep1, ray1.dir, ri, semiDiam, "none")
    return refSurf, t1
end

function surfClosestApproach(ray1::Ray, ray2::Ray, ray3::Ray; semiDiam = 5.0, ri=refIndexDefault, surfname="closest")
    t1,t2 = distClosestApproach(ray1, ray2)

    if ray1.dir != ray3.dir
        error("boresight and ytrace not parallel at start")
    end

    closep1 =ray3.base .+ t1 * ray3.dir
    #println("closest point on ray1 = $closep1  t1 = $t1")
    refSurf = referencePlane(surfname, closep1, ray3.dir, ri, semiDiam, "none")
    return refSurf, t1
end


function distClosestApproach(ray1::Ray, ray2::Ray)
    p1 = ray1.base
    d1 = ray1.dir
    p2 = ray2.base
    d2 = ray2.dir

    dp = p2 .- p1
    d1d1 = 1.0 #directions are normalized
    d2d2 = 1.0
    d1d2 = dot(d1, d2)
    dpd1 = dot(dp, d1)
    dpd2 = dot(dp, d2)

    denom = d1d1 * d2d2 - d1d2 * d1d2
    if abs(denom) < 1e-10
        println("rays are parallel")
        return Inf
    end

    num1 = dpd1 * d2d2 - dpd2 * d1d2
    num2 = dpd1 * d1d2 - dpd2 * d1d1

    t1 = num1 / denom
    t2 = num2 / denom


    return t1,t2
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
"""
    findRFP(geo; epsilon = 0.001, ydir = YAXIS)
        geo - geometry

    returns
        surfRFP - surface at rear focal plane
        fl - focal length
"""

function findRFP(geo; epsilon = 0.001, ydir = YAXIS)
    statusb,bore = traceGeometryRel(Ray(ORIGIN, ZAXIS), geo)
    statusy,ytrace = traceGeometryRel(Ray(ORIGIN+ epsilon * ydir, ZAXIS), geo)

    println("computeRearFocalPlane\nBore")
    printTrcCoords(statusb, bore, geo)
    println("ytrace")
    printTrcCoords(statusy, ytrace, geo)

    #check to make sure all made it through
    if statusb !=0 || statusy !=0
        return(1, NaN, NaN, ORIGIN, ORIGIN)
    end

    #compute intersection of bore and ytrace, xtrace
    #get last ray, rotate input rays to direction of bore, find when length that y, x are zero

    rayb = bore[end].ray
    rayy = ytrace[end].ray

    surfRFP, t1 = surfClosestApproach(rayb, rayy, semiDiam=5.0, surfname="RFP")
    if surfRFP == nothing
        error("rays are parallel")
    end

    rayyb = ytrace[begin].ray
    raybb = bore[begin].ray
    surfRPP, t2 = surfClosestApproach(rayyb, rayy, raybb, semiDiam=5.0, surfname="RPP")
    if surfRPP == nothing
        error("rays are parallel")
    end

    fl = norm(surfRFP.base.base - surfRPP.base.base) #may not have the right sign
    println("findRFP  fl = $fl")    

    return surfRFP, fl
end
