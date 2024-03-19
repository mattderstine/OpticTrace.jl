export traceGeometry, traceGeometryRel, surfAmpFunc
export traceGeometry!, traceGeometryRel!
export refractConic, refractSphere, reflectConic, cDiffuser
export reflectOAConic, reflectOAP, refractAsphere, reflectAsphere
export referencePlane, planeMirror, modFunc
export lensSinglet, lensASinglet, traceMonteCarloRays, conmputeRearFocalPlane
export sag, randomPointOnSquare, randomPointOnDisk
export attributesSurfaces, surfNormal

"""
    sag - compute the z coordinate in local coordinates for decendants of
    AbstractProfile

    returns z in local coordinates
"""
function sag(x::Float64, y::Float64, s::SurfProfileConic)
    if s.ϵ == 0. || s.curv == 0.
        z=s.curv * (x^2 + y^2) * 0.5
    else
        z=(1 - sqrt(1 - s.curv^2 * s.ϵ * (x^2 + y^2)))/(s.curv*s.ϵ)
    end
    z
end

function sag(x::Float64, y::Float64, s::SurfProfileSphere)

    z=1 - sqrt(1 - s.curv^2 * (x^2 + y^2))
end

function sag(x::Float64, y::Float64, s::SurfProfileOAConic)
    #println("at sag for SurfProfileOAConic - offset = $(s.offset)")
    a=SurfProfileConic(s.curv, s.ϵ)
    sg = sag(x-s.offset[1], y-s.offset[2], a)+s.offset[3] #add to offset for true sag
    #println("sag = $sg")
    sg
end

function sag(x::Float64, y::Float64, s::SurfProfileAsphere)
    #println("at sag for SurfProfileOAConic - offset = $(s.offset)")
    r2 = (x^2+y^2)
    r = sqrt(r2)
    sqrtarg = 1-s.ϵ * s.curv^2 * r2
    sg = s.curv * r2 /(1+ sqrt(1-s.ϵ * s.curv^2 * r2))+sum([s.a[i] * r2 * r^i for i in 1:length(s.a)])
    sg
end

function sag(x::Float64, y::Float64, s::SurfProfileCyl)
    if s.ϵ == 0. || s.curv == 0.
        z=s.curv * (y^2) * 0.5
    else
        z=(1 - sqrt(1 - s.curv^2 * s.ϵ * (y^2)))/(s.curv*s.ϵ)
    end
    z
end

function sag(x::Float64, y::Float64, s::SurfProfileToroid)
    #this is likely incorrect
    z= 1 - sqrt(1 - (s.curvY * y)^2 - (s.curvX * x)^2)
end

"""
    surfNormal - computes normal vectors for use by Makie to texture when rendering
    returns a unit vector
    returns normal in LOCAL coordinates
"""
function surfNormal(r::Point3, s::SurfProfileConic)
    #println("sag = $(sag(r[1],r[2],s))")

    if s.ϵ == 1
        denomI = 1
    else
        denomsq = 1 + s.curv * r[3] * (s.ϵ -1) * (s.curv * s.ϵ * r[3] -2)
        #next one is directly from Welford 4.30
        #denomsq = 1-2s.curv *(s.ϵ - 1) * r[3] + s.curv^2 * s.ϵ *(s.ϵ -1) * r[3]^2
        if denomsq < 0.0
            #println("r = $r   c = $(s.curv)  ϵ = $(s.ϵ)")
            #println("denomsq = $denomsq")
            denomI = NaN
        else
            denomI = 1/sqrt(denomsq)
        end
    end
    Vec3(-s.curv * r[1]*denomI, -s.curv * r[2]*denomI,
        (1.0 - s.curv * s.ϵ * r[3])*denomI )
end

function surfNormal(r::Point3, s::SurfProfileOAConic)
    #println("OA r = $r offset = $(s.offset) net = $(r .- s.offset)")
    surfNormal(r .- s.offset, SurfProfileConic(s.curv, s.ϵ))
end

function surfNormal(rr::Point3, s::SurfProfileAsphere)
    x = rr[1]
    y = rr[2]
    z = rr[3]
    r2 = (x^2+y^2)
    r = sqrt(r2)
    aspheresum = sum([i * r^(i-2) * s.a[i-2]  for i in 3:2+length(s.a)])
    normalize!([-x * (s.curv + aspheresum), -y*(s.curv +aspheresum), 1-s.curv * s.ϵ * z])
end

function surfNormal(r::Point3, s::SurfProfileCyl)
    #println("sag = $(sag(r[1],r[2],s))")

    if s.ϵ == 1
        denomI = 1
    else
        denomsq = 1 + s.curv * r[3] * (s.ϵ -1) * (s.curv * s.ϵ * r[3] -2)
        #next one is directly from Welford 4.30
        #denomsq = 1-2s.curv *(s.ϵ - 1) * r[3] + s.curv^2 * s.ϵ *(s.ϵ -1) * r[3]^2
        if denomsq < 0.0
            #println("r = $r   c = $(s.curv)  ϵ = $(s.ϵ)")
            #println("denomsq = $denomsq")
            denomI = NaN
        else
            denomI = 1/sqrt(denomsq)
        end
    end
    Vec3(0., -s.curv * r[2]*denomI, (1.0 - s.curv * s.ϵ * r[3])*denomI )
end

"""
    traceGeometry - trace a Ray through an array of OptSurfaces
    r       intial ray vector is global coordinates
    geo     the array of surfaces

    returns status, array of results from traceGeometry
"""
function traceGeometry(r::Ray, geo)
    trc = Vector{Trace}(undef, length(geo)+1)
    trc[1] = Trace(r, 1., 0., identityAmpMats()) #save the start ray etc
    curRay = r
    status = 0
    i = 1
    for surf in geo
        i += 1 #first element in trace is [2] in array
        status, trc[i] = traceSurf(curRay, surf)
        if status != 0
            break
        end
        curRay=trc[i].ray
    end
    status, trc[1:i] #only send the good ones!
end

"""
    traceGeometry! - trace a Ray through an array of OptSurfaces
    trc is a Vector{Trace} of length(geo) + 1 (the output)
    r is the intial ray vector is global coordinates
    geo is the array of surfaces

    returns status, length of trace
"""
function traceGeometry!(trc, r::Ray, geo)
    #trc = Vector{Trace}(undef, length(geo)+1)
    Trace!(trc[1], r, 1., 0., identityAmpMats()) #save the start ray etc
    curRay = r
    status = 0
    i = 1
    for surf in geo
        i += 1 #first element in trace is [2] in array
        status, - = traceSurf!(trc[i], curRay, surf)
        if status != 0
            break
        end
        curRay=trc[i].ray
    end
    status, i #only send the good ones!
end

function traceGeometryRel(rr::Ray, geo)
    trc = Vector{Trace}(undef, length(geo)+1)
    surf1=geo[1]
    curRay=Ray(surf1.toGlobalCoord(rr.base), surf1.toGlobalDir(rr.dir))
    trc[1] = Trace(curRay, 1., 0., identityAmpMats()) #save the start ray etc
    status = 0
    i = 1
    for surf in geo
        i += 1 #first element in trace is [2] in array
        status, trc[i] = traceSurf(curRay, surf)
        if status != 0
            break
        end
        curRay=trc[i].ray
    end
    status, trc[1:i] #only send the good ones!
end

function traceGeometryRel!(trc, rr::Ray, geo)
    #trc = Vector{Trace}(undef, length(geo)+1)
    surf1=geo[1]
    curRay=Ray(surf1.toGlobalCoord(rr.base), surf1.toGlobalDir(rr.dir))
    traceGeometry!(trc, curRay, geo)
end


"""
    traceSurf

    Compute the exit ray from a surface

    returns error, newRay, index, delta, polarizationMatrices

    newray is tuple of ray & direction
    polarzationMatrices is P & O matrices
    error
        0 - Normal
        1 - missed
        2 - TIR when refraction is expected - no polarization
        3 - Blocked by aperture/size of element

    newRay will contain best representation of the ray on error so it could be
    used in plotting and to continue a nonsequential raytrace
"""
function traceSurf(r::Ray,s::OptSurface)
    localRayStart=s.toLocalCoord(r.base)
    localRayDir = s.toLocalDir(r.dir)

    delta = deltaToSurf(Ray(localRayStart, localRayDir), s.profile)
    #=
    if debugFlag
        println("traceSurf...  Δ = $delta")
        println("start = $(r.base)")
    end
    =#
    if delta == NaN
        return (1, Trace(r, NaN, delta, identityAmpMats()) )
    end
    newRayBase = r.base + r.dir * delta
    newLocalBase = localRayStart + localRayDir * delta

    lnormal = surfNormal(newLocalBase, s.profile)

    normal = s.toGlobalDir(lnormal)
    #=
    if debugFlag
        println("local normal = $lnormal  global normal = $normal")
    end
    =#
    t, newRayDir, nIn = modFunc(Ray(newRayBase, r.dir), normal, s.mod)
    if !t
        return (2, Trace(Ray(newRayBase, newRayDir), nIn, delta, identityAmpMats()))
    end
    #=
    if debugFlag
        println("old dir = $(r.dir)  newdir = $newRayDir")
    end
    =#
    # can modify amplitude/polarization and direction of the ray, can't change intersection

    ampMats, newRayDir = surfAmpFunc(r.dir, newRayDir, normal, newLocalBase, s.mod, s.coating)

    return(0, Trace(Ray(newRayBase, newRayDir), nIn, delta, ampMats))
end

function Trace!(trc, ray, index, delta, ampdata)
    trc.ray = ray
    trc.nIn = index
    trc.delta = delta
    trc.pmatrix = ampdata
    return trc
end

function traceSurf!(trc, r::Ray,s::OptSurface)
    localRayStart=s.toLocalCoord(r.base)
    localRayDir = s.toLocalDir(r.dir)

    delta = deltaToSurf(Ray(localRayStart, localRayDir), s.profile)
    #=
    if debugFlag
        println("traceSurf...  Δ = $delta")
        println("start = $(r.base)")
    end
    =#
    if delta == NaN
        return (1, Trace!(trc, r, NaN, delta, identityAmpMats()) )
    end
    newRayBase = r.base + r.dir * delta
    newLocalBase = localRayStart + localRayDir * delta

    lnormal = surfNormal(newLocalBase, s.profile)

    normal = s.toGlobalDir(lnormal)
    #=
    if debugFlag
        println("local normal = $lnormal  global normal = $normal")
    end
    =#
    t, newRayDir, nIn = modFunc(Ray(newRayBase, r.dir), normal, s.mod)
    if !t
        return (2, Trace!(trc, Ray(newRayBase, newRayDir), nIn, delta, identityAmpMats()))
    end
    #=
    if debugFlag
        println("old dir = $(r.dir)  newdir = $newRayDir")
    end
    =#
    # can modify amplitude/polarization and direction of the ray, can't change intersection

    ampMats, newRayDir = surfAmpFunc(r.dir, newRayDir, normal, newLocalBase, s.mod, s.coating)

    return(0, Trace!(trc, Ray(newRayBase, newRayDir), nIn, delta, ampMats))
end

function traceSurf(r::Ray,s::ModelSurface)
    localRayStart=s.toLocalCoord(r.base)
    localRayDir = s.toLocalDir(r.dir)
    #println("\ntraceSurf - ModelSurface  name = $(s.surfname)")
    #println("aperture = $(s.aperture)")

    delta = deltaToSurf(Ray(localRayStart, localRayDir), s.profile)
    if delta == NaN
        return (1, Trace(r, s.refIndex, delta, identityAmpMats()))
    end
    newRayBase = r.base + r.dir * delta
    newLocalBase = localRayStart + localRayDir * delta


    if isAperture(s.aperture) && clipAperture(newLocalBase, s.aperture)
        stat = 3 #clipped
    else
        stat = 0 #ok to continue
    end

    return(stat, Trace(Ray(newRayBase, r.dir), s.refIndex, delta, identityAmpMats()))
end

function traceSurf!(trc, r::Ray,s::ModelSurface)
    localRayStart=s.toLocalCoord(r.base)
    localRayDir = s.toLocalDir(r.dir)
    #println("\ntraceSurf - ModelSurface  name = $(s.surfname)")
    #println("aperture = $(s.aperture)")

    delta = deltaToSurf(Ray(localRayStart, localRayDir), s.profile)
    if delta == NaN
        return (1, Trace!(trc, r, s.refIndex, delta, identityAmpMats()))
    end
    newRayBase = r.base + r.dir * delta
    newLocalBase = localRayStart + localRayDir * delta


    if isAperture(s.aperture) && clipAperture(newLocalBase, s.aperture)
        stat = 3 #clipped
    else
        stat = 0 #ok to continue
    end

    return(stat, Trace!(trc, Ray(newRayBase, r.dir), s.refIndex, delta, identityAmpMats()))
end




"""
    deltaToSurf - find distance along Ray to OptSurface
                calculation is in local coordinates
    returns distance to intersection of ray with surface

"""
function deltaToSurf(r::Ray, p::SurfProfileConic)
    x0, y0, z0 = r.base
    L, M, N = r.dir

    ffunc = p.curv * (x0^2 + y0^2 + z0^2 * p.ϵ) - 2 * z0
    gfunc = N - p.curv * (L * x0 + M * y0 + N * z0 * p.ϵ)

    if p.ϵ ==1.
        C = p.curv
    else
        C = p.curv * (L^2 + M^2 + p.ϵ * N^2)
    end
    #=
    if debugFlag
        println("N = $N  gfunc = $gfunc  ffunc = $ffunc")
        println("p.curv = $(p.curv)  C = $C")
    end
    =#

    if isapprox(C, 0., atol=1e-16)
        if isapprox(N, 0., atol=1e-16)
            qmiss = gfunc^2-p.curv*ffunc
            if qmiss < 0. || p.curv == 0.
                Δ =  NaN # root is imaginary or ray parallel to plane, miss
            else
                Δ = (gfunc - sqrt(qmiss))/p.curv
            end
        else
            Δ = ( 0.5*p.curv*(x0^2+y0^2)-z0)/N # negative sign removed 9/29/20
        end
    else
        qmiss = gfunc^2-C*ffunc
        if qmiss < 0.
            Δ =  NaN # root is imaginary, miss
        else
            Δ = (gfunc - sqrt(qmiss))/C
        end
    end
    Δ
end


# logic is flawed in this one
#change to add offset to ray to put it into the coordinate system of the offset parabola
function deltaToSurf(r::Ray, p::SurfProfileOAConic)
    #=
    if debugFlag
        println("deltaToSurf OAConic")
        println("base = $(r.base)  dir = $(r.dir)  offset = $(p.offset)  net = $(r.base .- p.offset)")
    end
    =#
    de = deltaToSurf(Ray(r.base .+ p.offset, r.dir), SurfProfileConic(p.curv, p.ϵ))
    #=
    if debugFlag
        println("delta = $de")
    end
    =#
    de
end

function deltaToSurf(r::Ray, profile::SurfProfileAsphere)
    x0, y0, z0 = r.base
    L, M, N = r.dir

    guess = deltaToSurf(r, SurfProfileConic(profile.curv, profile.ϵ))
    #=
    if debugFlag
        println("guess = $guess")
    end
    =#
    f(δ)=sag(x0 + L*δ, y0 + M * δ, profile) - z0 - N * δ
    Δl = find_zero(f, guess)
    Δl
end

function deltaToSurf(r::Ray, profile::SurfProfileCyl)
    x0, y0, z0 = r.base
    L, M, N = r.dir

    ffunc = p.curv * (y0^2 + z0^2 * p.ϵ) - 2 * z0
    gfunc = N - p.curv * (M * y0 + N * z0 * p.ϵ)

    if p.ϵ ==1.
        C = p.curv
    else
        C = p.curv * ( M^2 + p.ϵ * N^2)
    end
    #=
    if debugFlag
        println("N = $N  gfunc = $gfunc  ffunc = $ffunc")
        println("p.curv = $(p.curv)  C = $C")
    end
    =#


    if isapprox(C, 0., atol=1e-16)
        if isapprox(N, 0., atol=1e-16)
            qmiss = gfunc^2-p.curv*ffunc
            if qmiss < 0. || p.curv == 0.
                Δ =  NaN # root is imaginary or ray parallel to plane, miss
            else
                Δ = (gfunc - sqrt(qmiss))/p.curv
            end
        else
            Δ = ( 0.5*p.curv*(x0^2+y0^2)-z0)/N # negative sign removed 9/29/20
        end
    else
        qmiss = gfunc^2-C*ffunc
        if qmiss < 0.
            Δ =  NaN # root is imaginary, miss
        else
            Δ = (gfunc - sqrt(qmiss))/C
        end
    end
    Δ
end



"""
    identityAmpMats

"""
function identityAmpMats()
    #identity = [[1.0, 0.0, 0.0] [ 0.0, 1.0, 0.0] [ 0.0, 0.0, 1.0]]
    AmpData(identityPol, identityPol, [1.])
    #retirm matrices for P & O and unity transmission
end

"""
    modFunc

"""
function modFunc(ray::Ray, normal::Vec3, d::AbstractBendDielectric)
    r = ray.base
    a = ray.dir

    cosI = a ⋅ normal
    if cosI <0. #hitting surface "backwards"
        nIn = d.refIndexOut
        nOut = d.refIndexIn
        cosI = -cosI
    else
        nIn = d.refIndexIn
        nOut = d.refIndexOut
    end

    tirtest = nOut^2 - nIn^2 * (1-cosI^2)
    if tirtest < 0
        return false, -(a - (2. * (a ⋅ normal)) .* normal), nIn #TIR direction
    end

    cosIp = sqrt(tirtest)/nOut

    kparam = (nOut * cosIp - nIn * cosI)
    return true, (nIn .* a + kparam .* normal)./nOut, nIn
end



function modFunc(ray::Ray, normal::Vec3, d::AbstractBendMirror)
    r = ray.base
    a = ray.dir
    cosI = a ⋅ normal
    if cosI <0. #hitting surface "backwards"
        nIn = d.refIndexOut
    else
        nIn = d.refIndexIn
    end

    true, (a - (2. * (cosI)) .* normal), nIn  # Welford 4.46
end

function modFunc(ray::Ray, normal::Vec3, d::CDiffuser)
    r = ray.base
    a = ray.dir
    perpapprox = zeros(Float64,3)
    absa = abs.(a)
    perpapprox[findfirst(isequal(minimum(absa)), absa)]=1.  #create the vector almost perpendicular to the direction

    perp1 = normalize!(cross(a, perpapprox))
    perp2 = cross(a, perp1)
    dir = 2 .* rand(2) .- [1., 1.]
    while norm(dir)>1.
        dir = 2 .* rand(2) .- [1., 1.]
    end
    displace = d.tanθ .* (dir[1] .* perp1 + dir[2] .* perp2)

    cosI = a ⋅ normal
    if cosI <0. #hitting surface "backwards"
        nIn = d.refIndexOut
    else
        nIn = d.refIndexIn
    end

    true, normalize(a + displace), nIn
end



"""
    surfAmpFunc

"""
function surfAmpFunc(dirIn::Vec3, dirOut::Vec3, normal::Vec3, newLocalBase::Point3, sndex::DielectricT,amp::AmpParam)
    identityAmpMats(), dirOut
end

function surfAmpFunc(dirIn::Vec3, dirOut::Vec3, normal::Vec3, newLocalBase::Point3, ndex::MirrorR,amp::AmpParam)
    identityAmpMats(), dirOut
end

function surfAmpFunc(dirIn::Vec3, dirOut::Vec3, normal::Vec3, newLocalBase::Point3, ndex::CDiffuser,amp::AmpParam)
    identityAmpMats(), dirOut
end

attributesSurfaces::Dict{String,Any} = Dict{String, Any}()

function getAmpParams(s::String; attributesSurfaces=attributesSurfaces)
    amp = get!(attributesSurfaces, s, AmpParam(s))
    return amp
end

function getAmpParams(s::AbstractAmplitudeParam; attributesSurfaces=attributesSurfaces)
    return s
end


"""
    refractConic(
        surfname,
        base, dir,
        refractiveindex in, refractive index out,
        curvature,
        ϵ,
        semiDiameter,
        coating)
"""
function refractConic(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    ϵ::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :yellow, 
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3, Nothing}=nothing)

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal,ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileConic( c, ϵ),
        DielectricT(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end

"""
    refractSphere(
        surfname,
        base, dir,
        refractiveindex in, refractive index out,
        curvature,
        semiDiameter,
        coating)
"""
function refractSphere(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :aquamarine2, 
    attributesSurfaces = attributesSurfaces, ydir::Union{Vec3, Nothing}=nothing)

    refractConic(surfname,
        pointInPlane,
        planenormal,
        rinIn,
        rinOut,
        c,
        1., # ϵ
        semiDiam,
        coating
        ;color, attributesSurfaces, ydir)
end



"""
    reflectConic(
        surfname,
        base, dir,
        refractiveindex in, refractive index out,
        curvature,
        ϵ,
        semiDiameter,
        coating)
"""
function reflectConic(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    ϵ::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :aquamarine2,
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3, Nothing}=nothing)

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal,ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileConic( c, ϵ),
        MirrorR(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end

"""
    cDiffuser(
        surfname,
        base, dir,
        refractiveindex in, refractive index out,
        half deflection angle,
        semiDiameter,
        coating)
"""
function cDiffuser(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    θ::Float64,
    semiDiam::Float64,
    coating::APorString
    ; color = :brown, 
    ydir::Union{Vec3, Nothing}=nothing)

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        NoProfile(0.),
        CDiffuser(tan(θ), rinIn, rinOut),
        AmpParam(coating),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end



function reflectOAConic(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    ydir::Vec3,
    offset::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    ϵ::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :aquamarine2, attributesSurfaces = attributeSurfaces)

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileOAConic( c, ϵ, offset),
        MirrorR(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end

"""
    convention is that rotation around z axis, ray going -z direction reflects to +y dir
"""
reflectOAP(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    ydir::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :blue, attributesSurfaces = attributesSurfaces
    ) = reflectOAConic(surfname,
                        pointInPlane,
                        planenormal,
                        ydir,
                        Vec3(0, 1/c, -0.5/c),
                        rinIn,
                        rinOut,
                        c,
                        0.0,
                        semiDiam,
                        coating; color, attributesSurfaces)
#=

    refracting asphere
    reflecting asphere
    OA conic, refract & reflect - not tested  yet
    prism
    hologram/grating


=#
#=

Overloads for referencePlane

=#

"""
    refractAsphere
    

"""
function refractAsphere(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    ϵ::Float64,
    asphere::AbstractVector{Float64},
    semiDiam::Float64,
    coating::APorString
    ;color = :aquamarine,  
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3, Nothing}=nothing)

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileAsphere( c, ϵ, asphere),
        DielectricT(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end


"""
    reflectAsphere
    like refractAsphere but reflects

"""
function reflectAsphere(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    ϵ::Float64,
    asphere::AbstractVector{Float64},
    semiDiam::Float64,
    coating::APorString
    ;color = :aquamarine, 
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3, Nothing}=nothing)

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileAsphere( c, ϵ, asphere),
        MirrorR(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir, 
        color
        )
end

function sag(x::Float64, y::Float64, s::NoProfile)
    0.
end

function gbRadius(aperture::SizeLens, profile::NoProfile)
    aperture.semiDiameter
end

function gbWidths(a::SizeLens, p::NoProfile)
    diam = 2a.semiDiameter
    SVector(diam, diam, 0.)
end

function surfNormal(r::Point3, s::NoProfile)
    Vec3(0., 0., 1.)
end

function deltaToSurf(r::Ray, p::NoProfile)
    x0, y0, z0 = r.base
    L,M,N = r.dir

    if N ≈ 0.
        Δ =  NaN #ray parallel to flat surface
    else
        Δ = -z0/N
    end
    Δ
end

function modFunc(ray::Ray, normal::Vec3, d::NoBendIndex)
    true, ray.dir, d.refIndexIn
end

function surfAmpFunc(dirIn::Vec3, dirOut::Vec3, normal::Vec3, newRayBase::Point3, ndex::NoBendIndex,amp::NoAmpParam)
    identityAmpMats(), dirOut
end

function referencePlane(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :khaki3, 
    ydir::Union{Vec3, Nothing}=nothing)

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        NoProfile( 0.),
        NoBendIndex(rinIn),
        NoAmpParam(coating),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end


"""
    planeMirror is just a mirror with zero curvature and conic
    could use NoProfile to speed things up but then would need to
    potentially overload other functions
"""
function planeMirror(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :blue, 
    attributesSurfaces = attributesSurfaces,  #attributesSurfaces is a global dictionary
    ydir::Union{Vec3, Nothing}=nothing)

    reflectConic(surfname,
        pointInPlane,
        planenormal,
        rinIn,
        rinOut,
        0.,
        0.,
        semiDiam,
        coating; color, attributesSurfaces, ydir)
end

"""
helper function for singlet lens
"""
function lensSinglet(base, dir, curv1, curv2, thick, lambda, riFunc, semiDiam; order = "forward", lensname = "Singlet", coating = "default")
    ri = riFunc(lambda)
    base1 = base + thick .* dir
    if (order != "reverse" )
        lens = [
        refractSphere("$(lensname)_1", base, dir, refIndexDefault, ri,
            curv1, semiDiam, coating),
        refractSphere("$(lensname)_2", base1, dir, ri,refIndexDefault,
            curv2, semiDiam, coating),

        ]
    else
        

        lens = [
        refractSphere("$(lensname)_2", base, dir, refIndexDefault, ri, 
            -curv2, semiDiam, coating),
        refractSphere("$(lensname)_1", base1, dir, ri, refIndexDefault,
            -curv1, semiDiam, Base.compute_assumed_setting)
        ]
    end
    lens
end



"""
lensASinglet(base, dir, curv1, ϵ1, aphere1, curv2, ϵ2, aphere2, 
      thick, lambda, riFunc, semiDiam; order = "forward", lensname = "ASinglet")
"""
function lensASinglet(base, dir, curv1, ϵ1, aphere1, curv2, ϵ2, aphere2, thick, lambda, riFunc, semiDiam; order = "forward", lensname = "ASinglet", coating="default")
    ri = riFunc(lambda)
    # should check if the input is really an asphere. if not make the surface spherical #ToDo
    base1 = base + thick .* dir
    if (order != "reverse" )
        lens = [
        refractAsphere("$(lensname)_1", base, dir, refIndexDefault, ri,
            curv1, ϵ1, aphere1, semiDiam, coating),
        refractAsphere("$(lensname)_2", base1, dir, ri,refIndexDefault,
            curv2, ϵ2, aphere2, semiDiam, coating)
        ]
    else

        lens = [
        refractAsphere("$(lensname)_2", base, dir, refIndexDefault, ri, 
            -curv2, ϵ2, -aphere2, semiDiam, coating),
        refractAsphere("$(lensname)_1", base1, dir, ri, refIndexDefault, 
            -curv1, ϵ1, -aphere1, semiDiam, coating)
        ]
    end
    lens
end


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
    statusy,ytrace = traceGeometryRel(Ray(SVector(0., epsilon, 0.), ZAXIS), geo)
    statusx,xtrace = traceGeometryRel(Ray(SVector(epsilon, 0., 0.), ZAXIS), geo)

    printTrcCoords(statusb, bore, geo)

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
