4+4

using OpticTrace
using StaticArrays
#=

attributesSurfaces = Dict{String, Any}()


=#

##

testGeometrytime = [
    referencePlane("start",Point3(0., 0., -50.), ZAXIS, refIndexDefault , 10., "testcoat"),
    refractSphere("Lens1i", Point3(0., 0., 0.), ZAXIS, refIndexDefault, 1.5, 1. / 50., 20., "testcoat"),
    refractSphere("Lens10", Point3(0., 0., 10.), ZAXIS, 1.5, refIndexDefault, -1. /50., 20., "testcoat"),


    referencePlane("end", Point3(0.,0, 40.), ZAXIS, refIndexDefault, 3., "testcoat")
]
figt, scenet = plotGeometry3D(testGeometrytime)
display(figt)

function testrace(n)
    sts = zeros(Int64, n)
    - , trc = traceGeometryRel(Ray(Point3(0.0, rand(), 0.0), ZAXIS), testGeometrytime)
    for i in 1:n
        sts[i], trc = traceGeometryRel(Ray(Point3(0.0, rand(), 0.0), ZAXIS), testGeometrytime)
    end
    sts
end

function testrace2(n)
    sts = zeros(Int64, n)
    - , trc = traceGeometryRel(Ray(Point3(0.0, rand(), 0.0), ZAXIS), testGeometrytime)
    for i in 1:n
        sts[i], len = traceGeometryRel!(trc,Ray(Point3(0.0, rand(), 0.0), ZAXIS), testGeometrytime)
    end
    sts
end

@time testrace(1)
@time testrace2(1)
@time st, trc = traceGeometryRel(Ray(ORIGIN, ZAXIS), testGeometrytime)
@time st, len = traceGeometryRel!(trc, Ray(ORIGIN, ZAXIS), testGeometrytime)

ottest(n) = [Point(0.0, 0.0, i) for i in 1:n]

zray = Ray(ORIGIN, ZAXIS)
@time st, trc = OpticTrace.traceSurf(zray, testGeometrytime[1])

function testtraceSurf(surf)
    st, trace =OpticTrace.traceSurf(zray, surf)
    st == 0
end

testtraceSurf(testGeometrytime[3])

testGeometrytime[3]
function testtoGlobalDir()
    for i in 0.:0.1:1.0
        testGeometrytime[2].toGlobalDir(ZAXIS)
    end
end

function testtoGlobalDir2()
    for i in 0.:0.1:1.0
        OpticTrace.matmul3(testGeometrytime[2].toGCMat, ZAXIS)
    end
end

@time testtoGlobalDir()
@time testtoGlobalDir2()

a, b = findPerpenMap(ZAXIS)

@time c = b(ZAXIS)
dds = SMatrix(hcat(XAXIS, YAXIS, ZAXIS))
d = LinearMap(dds)

function testd()
    @time a = [ZAXIS == d(ZAXIS) for i in 1:100]
end

testd()

function matmul(a::SMatrix)
    function b(x)
        c = a * x
        return SVector(c)
    end
    return b
end

f = matmul(dds)
f(ZAXIS)
f(ZAXIS) ==ZAXIS
ZAXIS
function testmatmul(dir)
     ZAXIS == f(dir)
end
@time testmatmul(ZAXIS)

function tt!(dir)
    dir =  dd*dir
end
ttt = ZAXIS
@time tt!(ttt)
t4() = [ZAXIS==f(ZAXIS) for i in 1:100]
@time ZAXIS==d(ZAXIS)

@time t4()

t5(dest, rot, start) = [ZAXIS == mul!(dest, rot, start) for i in 1:100]
function testdest()
    dest = [0.,0.,0.]
    rot = dds
    start = ZAXIS
   @showtime mul!(dest, rot, start)
   dest == ZAXIS
end

@time testdest()


a=[0.,0.,0.]
function multest(a,dir)
    #a=[0.,0.,0.]
    
    return [ZAXIS==mul!(a, dd, dir) for i in 1:100]
end

@time multest(a, XAXIS)

bb =  LinearMap(dd)
mul!(a,bb,ZAXIS)

@time c=bb(I)
@time c*a

SVector{Float64}[a...]
b = @SVector abs2

@time z3 = matmul3(dds, ZAXIS)

function multest3(dd)

    a = rand(3)
    @time [ZAXIS==OpticTrace.matmul3(dd, (1., 0., i)) for i in 1:100]
end

multest3(dds)

o = (0.0,0.0,0.0)
d = (1.0, 0.0, 0.0)

osv = Point3(o)
dsv = Vec3(d)

function rayalloctest()
    osv = Point3(rand(3))
    dsv = Vec3(rand(3))
    @timev a = Ray(osv,dsv)
    return (a.base[1])
end

rayalloctest()

##

using StaticArrays
abstract type AbstractRay{T,N} end


struct Ray{T,N} <: AbstractRay{T,N}
    base::SVector{N,T}
    dir::SVector{N,T}
end

function propray(ray::Ray, d::Float64)
    Ray(SVector(ray.base .+ d * ray.dir), ray.dir)
end

function rayalloctest()
    osv = SVector{3}(rand(3))
    dsv = SVector{3}(rand(3))
    @timev a = propray(Ray(osv,dsv), 3.)
    return (a.base)
end

rayalloctest()

##
rp=referencePlane("start", ORIGIN, ZAXIS, 1., 5., "nocoating"; ydir = nothing)

