4+4

using OpticTrace

#=

attributesSurfaces = Dict{String, Any}()


=#

using BenchmarkTools
using Test
using StaticArrays

@benchmark sag(x, y, s) setup=(x = rand(); y=rand(); s = OpticTrace.SurfProfileConic(0.01, 0.))

@allocations pp = sag(0., 1., OpticTrace.SurfProfileConic(0.01, 0.))
pp

@allocations surfprofile = OpticTrace.SurfProfileAsphere(0.01, 0., [0., 0.1, -0.2])
@allocations pp1 = sag(0. , 1., surfprofile)
@allocations sum(2*i for i in 1:20)

f(j) = sum(2*i for i in 1:j)
@allocations f(20)
pp1
k = [2., 3.]


@benchmark traceGeometry(r, geo ) setup=( r=Ray(Point3(0., 10.0, -50.), Vec3(0., -10.0/50., sqrt(1-(10.0/50.)^2))); geo=testGeometry)

trc = [Trace(r, 1., 0., OpticTrace.identityAmpMats()) for i in 1:length(testGeometry)+1]
r=Ray(Point3(0., 10.0, -50.), Vec3(0., -10.0/50., sqrt(1-(10.0/50.)^2)))
OpticTrace.Trace!(trc[1], r, 1., 0., OpticTrace.identityAmpMats()) #save the start ray etc
@benchmark traceGeometry!(trc, r, geo ) setup=(trc = [Trace(Ray(Point3(0., 10.0, -50.), Vec3(0., -10.0/50., sqrt(1-(10.0/50.)^2))), 1., 0., OpticTrace.identityAmpMats()) for i in 1:length(testGeometry)+1]; r=Ray(Point3(0., 10.0, -50.), Vec3(0., -10.0/50., sqrt(1-(10.0/50.)^2))); geo=testGeometry)

using StaticArrays

a = StaticArrays.Vec3(0, 0, 0)