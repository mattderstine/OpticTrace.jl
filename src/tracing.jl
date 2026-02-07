export traceGeometry, traceGeometryRel, surfAmpFunc
export traceGeometry!, traceGeometryRel!

export attributesSurfaces, surfNormal, modFunc

"""
    sag - compute the z coordinate in local coordinates for decendants of
    AbstractProfile

    returns z in local coordinates
"""
function sag(x::T, y::T, s::SurfProfileConic{T}) where T<:Real
    r2 = (x^2+y^2)
    sqrtarg = 1-s.ϵ * s.curv^2 * r2
    if sqrtarg < 0.
        z=NaN
    elseif s.curv ==0.
        z = 0.
    elseif s.ϵ == 0. 
        z=s.curv * r2 * 0.5
    else
       #z=(1 - sqrt(1 - s.curv^2 * s.ϵ * (x^2 + y^2)))/(s.curv*s.ϵ)
       z= s.curv * r2 /(1+ sqrt(sqrtarg))
    end
    z
end

function sag(x::T, y::T, s::SurfProfileSphere{T}) where T<:Real

    z=1 - sqrt(1 - s.curv^2 * (x^2 + y^2))
end

function sag(x::T, y::T, s::SurfProfileOAConic{T}) where T<:Real
    #println("at sag for SurfProfileOAConic - offset = $(s.offset)")
    a=SurfProfileConic(s.curv, s.ϵ)
    sg = sag(x-s.offset[1], y-s.offset[2], a)+s.offset[3] #add to offset for true sag
    #println("sag = $sg")
    sg
end


"""
    sag - compute the z coordinate in local coordinates for decendants of
    AbstractProfile

    returns z in local coordinates
    For SurfProfileAsphere, compute from 3rd order term onwards


"""
function sag(x::T, y::T, s::SurfProfileAsphere{T}) where T<:Real
    r2 = (x^2+y^2)
    r = sqrt(r2)
    sqrtarg = 1-s.ϵ * s.curv^2 * r2
    if sqrtarg < 0.
        return NaN
    end
    #sg = s.curv * r2 /(1+ sqrt(sqrtarg))+sum([ss * r2 * r^i for (i,ss) in enumerate(s.a)])

    asp = 0.
    for ss in Iterators.reverse(s.a)
        asp = (asp + ss) * r
    end
    asp *= r2 #apsehreic coefficients array starts at r3
    sg = s.curv * r2 /(1+ sqrt(sqrtarg))+asp

    sg
end

function sag(x::T, y::T, s::SurfProfileEvenAsphere{T}) where T<:Real
    r2 = (x^2+y^2)

    sqrtarg = 1-s.ϵ * s.curv^2 * r2
    if sqrtarg < 0.
        return NaN
    end
    #sg = s.curv * r2 /(1+ sqrt(sqrtarg))+sum([ss * r2 * r^i for (i,ss) in enumerate(s.a)])

    asp = 0.
    for ss in Iterators.reverse(s.a)
        asp = (asp + ss) * r2
    end
    asp *= r2 # s starts at order 4 unlike zemax which starts at order 2
    sg = s.curv * r2 /(1+ sqrt(sqrtarg))+asp

    sg
end

function sag(x::T, y::T, s::SurfProfileCyl{T}) where T<:Real
    if s.ϵ == 0. || s.curv == 0.
        z=s.curv * (y^2) * 0.5
    else
        z=(1 - sqrt(1 - s.curv^2 * s.ϵ * (y^2)))/(s.curv*s.ϵ)
    end
    z
end

function sag(x::T, y::T, s::SurfProfileToroid{T}) where T<:Real
    #this is likely incorrect
    z= 1 - sqrt(1 - (s.curvY * y)^2 - (s.curvX * x)^2)
end

"""
    surfNormal - computes normal vectors for use by Makie to texture when rendering
    returns a unit vector
    returns normal in LOCAL coordinates
"""
function surfNormal(r::Point3{T}, s::SurfProfileConic{T}) where T<:Real
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

function surfNormal(r::Point3{T}, s::SurfProfileOAConic{T}) where T<:Real
    #println("OA r = $r offset = $(s.offset) net = $(r .- s.offset)")
    surfNormal(r .- s.offset, SurfProfileConic(s.curv, s.ϵ))
end

function surfNormal(rr::Point3{T}, s::SurfProfileAsphere{T}) where T<:Real
    x = rr[1]
    y = rr[2]
    z = rr[3]
    r2 = (x^2+y^2)
    r = sqrt(r2)

    asp = 0.
    for (i,ss) in Iterators.reverse(enumerate(s.a))
        asp = (asp + (i+2) * ss) * r
    end
    asp *= r

    norm= normalize([-x * (s.curv + asp), -y*(s.curv +asp), 1-s.curv * s.ϵ * z])
    Vec3(norm[1], norm[2], norm[3])
end

function surfNormal(rr::Point3{T}, s::SurfProfileEvenAsphere{T}) where T<:Real
    x = rr[1]
    y = rr[2]
    z = rr[3]
    r2 = (x^2+y^2)

    #aspheresum = sum([i * r^(i-2) * s.a[i-2]  for i in 3:2+length(s.a)])
    #aspheresum = sum([(i+2) * r^i * s.a[i]  for i in 1:length(s.a)])
    asp = 0.

    for (i,ss) in Iterators.reverse(enumerate(s.a))
        asp = (asp + (2*(i+1)) * ss )*r2

    end
    asp *= sqrt(r2)

    norm= normalize([-x * (s.curv + asp), -y*(s.curv +asp), 1-s.curv * s.ϵ * z])
    Vec3(norm[1], norm[2], norm[3])
end

function surfNormal(r::Point3{T}, s::SurfProfileCyl{T}) where T<:Real
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
function traceGeometry(r::Ray{3,T}, geo) where T<:Real
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
function traceGeometry!(trc::Vector{Trace}, r::Ray{3,T}, geo) where T<:Real
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

function traceGeometryRel(rr::Ray{3,T}, geo) where T<:Real
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

function traceGeometryRel!(trc, rr::Ray{3,T}, geo) where T<:Real
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
function traceSurf(r::Ray{3,T},s::OptSurface{3,T}) where T<:Real
    localRayStart=s.toLocalCoord(r.base)
    localRayDir = s.toLocalDir(r.dir)

    delta = deltaToSurf(Ray(localRayStart, localRayDir), s.profile)
    #=
    if debugFlag
        println("traceSurf...  Δ = $delta")
        println("start = $(r.base)")
    end
    =#
    if delta == NaN #missed
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
    if !t #TIR or other failure
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

function Trace!(trc::Trace{T}, ray, index, delta, ampdata) where T<:Real
    trc.ray = ray
    trc.nIn = index
    trc.delta = delta
    trc.pmatrix = ampdata
    return trc
end

function traceSurf!(trc, r::Ray{3,T},s::OptSurface{3,T}) where T<:Real
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

function traceSurf(r::Ray{3,T},s::ModelSurface{3,T}) where T<:Real
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

function traceSurf!(trc, r::Ray{3,T},s::ModelSurface{3,T}) where T<:Real
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
function deltaToSurf(r::Ray{3,T}, p::SurfProfileConic{T}) where T<:Real
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
function deltaToSurf(r::Ray{3,T}, p::SurfProfileOAConic{T}) where T<:Real
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

function deltaToSurf(r::Ray{3,T}, profile::AbstractAsphericProfile{T}) where T<:Real
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

function deltaToSurf(r::Ray{3,T}, profile::SurfProfileCyl{T}) where T<:Real
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


const AMPPERFECTSURFACE = AmpData(identityPol, identityPol, [1.])
"""
    identityAmpMats

"""
identityAmpMats() = AMPPERFECTSURFACE
#=
function identityAmpMats()
    #identity = [[1.0, 0.0, 0.0] [ 0.0, 1.0, 0.0] [ 0.0, 0.0, 1.0]]
    AmpData(identityPol, identityPol, [1.])
    #retirm matrices for P & O and unity transmission
end
=#

"""
    modFunc

"""
function modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::S) where {S <: AbstractBendDielectric, T<:Real}

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



function modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::S) where {S <: AbstractBendMirror, T<:Real}
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

function modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::CDiffuser{T}) where T<:Real
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


function surfNormal(r::Point3{T}, s::NoProfile{T}) where T<:Real
    Vec3(0., 0., 1.)
end

function deltaToSurf(r::Ray{3,T}, p::NoProfile{T}) where T<:Real
    x0, y0, z0 = r.base
    L,M,N = r.dir

    if N ≈ 0.
        Δ =  NaN #ray parallel to flat surface
    else
        Δ = -z0/N
    end
    Δ
end

function modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::NoBendIndex{T}) where T<:Real
    true, ray.dir, d.refIndexIn
end



"""
    surfAmpFunc

"""
function surfAmpFunc(dirIn::Vec3{T}, dirOut::Vec3{T}, normal::Vec3{T}, newLocalBase::Point3{T}, sndex::B, amp) where {T<:Real, B <: AbstractBendType{T}}
    identityAmpMats(), dirOut
end

#=
function surfAmpFunc(dirIn::Vec3, dirOut::Vec3, normal::Vec3, newLocalBase::Point3, ndex::MirrorR,amp::AmpParam)
    identityAmpMats(), dirOut
end

function surfAmpFunc(dirIn::Vec3, dirOut::Vec3, normal::Vec3, newLocalBase::Point3, ndex::CDiffuser,amp::AmpParam)
    identityAmpMats(), dirOut
end

=#

attributesSurfaces::Dict{String,Any} = Dict{String, Any}()

function getAmpParams(s::String; attributesSurfaces=attributesSurfaces)
    amp = get!(attributesSurfaces, s, AmpParam(s))
    return amp
end

function getAmpParams(s::S; attributesSurfaces=attributesSurfaces) where S <: AbstractAmplitudeParam
    return s
end
