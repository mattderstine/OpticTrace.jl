## Testing Set Up Figures
4+4

using OpticTrace # load all the lens design stuff
using StaticArrays
using FileIO # MeshIO should also be installed
using LinearAlgebra
using CoordinateTransformations
using IterTools
using Roots
using DataInterpolations
using StatsBase
using GeometryBasics
using GLMakie
using Printf
using MeshIO


#set up simple geometry and plot it

figs = [Figure(size = (1600,1200),backgroundcolor = RGBf(0.98, 0.98, 0.98), ) for i in 1:5]
#figs=multipleFigures()


## Test the diffuser fig1

"""test the diffuser
# 11/26/20 works as expected
"""
function testDiffuser(fig; pnts=20)

    testradius = 2.
    distance = testradius / tan(10.0 * π/180.)
    testDiffGeo = [
        referencePlane("start", ORIGIN, ZAXIS, 1., 5., "nocoating"; ydir = nothing),
        cDiffuser("diff", Point3(0., 0., 2.), ZAXIS, 1., 1., 10.0 * π/180., 5., "nocoating"),
        referencePlane("end", Point3(0., 0., distance+2), ZAXIS, 1., 5., "nocoating" ; ydir = Vec3(sin(π/2), cos(π/2), 0.))
    ]

    outer_padding = 30
    lscene = LScene(fig[1,2], scenekw = (camera = cam3d_cad!, raw = false))
    linesegments!(lscene,[Makie.Point3(0, 0,0) => Makie.Point3(0,0,1)],color=:green, linewidth=2)
    linesegments!(lscene,[Makie.Point3(0, 0,0) => Makie.Point3(1,0,0)],color=:blue, linewidth=2)
    linesegments!(lscene,[Makie.Point3(0, 0,0) => Makie.Point3(0,1,0)],color=:red, linewidth=2)

    plotGeometry3D!(lscene, testDiffGeo)
    for j in 1:pnts
        trcAndPlotRay!(lscene, Ray(ORIGIN, ZAXIS), testDiffGeo)
    end

    ax1 = fig[1,3] = Axis(fig, title = "Spots radius = $testradius")
    ax1.aspect = DataAspect()


#should produce a circle with radius "testradius".
    aa()=traceGeometryRel(Ray(Point(0.,0., 0.),ZAXIS), testDiffGeo)
    pnts = [SVector{2}(aa()[2][4].ray.base[1:2]) for i in 1:100000]
    scatter!(ax1,pnts)

    fig
end

figDiffuserTest=testDiffuser(figs[1], pnts=100)



display(figDiffuserTest)

## test 2

testgeo = [
    refractSphere("first", ORIGIN,ZAXIS,1., 1.5, 0., 5., "testcoat";ydir = nothing),
    refractSphere("second", Point3(0., 0., 3.), ZAXIS, 1.5, 1., -0.1, 5., "testcoat"),

    refractSphere("third", Point3(0., 0., 3.5), ZAXIS, 1., 1.5, 0.1, 5., "testcoat"; ydir = Vec3(sin(π/2), cos(π/2), 0.)),
    refractSphere("fourth", Point3(0., 0., 6.5), ZAXIS, 1.5, 1., 0., 5., "testcoat"; ydir = Vec3(sin(π/2), cos(π/2), 0.)),

#    refractSphere("last", SVector(0., 0., 15.), ZAXIS, 1. , 1. , 0., 10., "testcoat")
    referencePlane("RP-last",Point3(0., 0., 15.), ZAXIS, 1. , 7., "testcoat")

]

refIndexDefault = 1. 

testGeometry = [
    referencePlane("start",Point3(0., 0., -50.), ZAXIS, refIndexDefault , 10., "testcoat"),
    refractSphere("Lens1i", Point3(0., 0., 0.), ZAXIS, refIndexDefault, 1.5, 1. / 50., 20., "testcoat"),
    refractSphere("Lens10", Point3(0., 0., 10.), ZAXIS, 1.5, refIndexDefault, -1. /50., 20., "testcoat"),
    referencePlane("image", Point3(0., 0., 50.), ZAXIS, refIndexDefault, 3., "testcoat")
]


function test2(fig)
    lscene = LScene(fig[1,2], scenekw = (camera = cam3d_cad!, raw = false))
    linesegments!(lscene,[Makie.Point3f(0, 0,0) => Makie.Point3f(0,0,1)],color=:green, linewidth=2)
    linesegments!(lscene,[Makie.Point3f(0, 0,0) => Makie.Point3f(1,0,0)],color=:blue, linewidth=2)
    linesegments!(lscene,[Makie.Point3f(0, 0,0) => Makie.Point3f(0,1,0)],color=:red, linewidth=2)

    scene1 = plotGeometry3D!(lscene, testgeo)

    scene2 = trcAndPrintPlotRay!(lscene, Ray(Point3(0., 5.0, -10.), Vec3(0., 0., 1.)), testgeo )

    #refPlane = referencePlane("test reference plane",ORIGIN,YAXIS, 1., 5., "testCoating")
    lscene2 = LScene(fig[1,3], scenekw = (camera = cam3d_cad!, raw = false))
    linesegments!(lscene2,[Makie.Point3f(0, 0,0) => Makie.Point3f(0,0,1)],color=:green, linewidth=2)
    linesegments!(lscene2,[Makie.Point3f(0, 0,0) => Makie.Point3f(1,0,0)],color=:blue, linewidth=2)
    linesegments!(lscene2,[Makie.Point3f(0, 0,0) => Makie.Point3f(0,1,0)],color=:red, linewidth=2)

    scene3 = plotGeometry3D!(lscene2, testGeometry)

    scene4 = trcAndPrintPlotRay!(lscene2, Ray(Point3(0., 10.0, -50.), Vec3(0., 0., 1.)), testGeometry)

    scene5 = trcAndPrintPlotRay!(lscene2, Ray(Point3(0., 10.0, -50.), Vec3(0., -10.0/50., sqrt(1-(10.0/50.)^2))), testGeometry)
    fig
end

display(test2(figs[2]))


gbmsh = merge(GeometryBasics.mesh.(testgeo))

save("test/testgeo.stl", gbmsh)
## test 3

l2i = refractSphere("Lens2i", Point3(-131., 0., 40.), XAXIS, refIndexDefault,
    1.5, 0., 20., "coatingTransmitFresnel") # gaussian power of 1/100
l2o = refractSphere("Lens2o", Point3(-141., 0., 40.), XAXIS, 1.5,
    refIndexDefault, 1.0 /50., 20., "coatingTransmitFresnel")

testGeometry2 = [
    referencePlane("start",Point3(0., 0., -50.), ZAXIS, refIndexDefault , 10., "testcoat"),
    refractSphere("Lens1i", Point3(0., 0., 0.), ZAXIS, refIndexDefault, 1.5, 1. / 50., 20., "testcoat"),
    refractSphere("Lens10", Point3(0., 0., 10.), ZAXIS, 1.5, refIndexDefault, -1. /50., 20., "testcoat"),
    planeMirror("Mirror", Point3(0., 0., 40.), Vec3(cos(π /4), 0., sin(π/4)),
    refIndexDefault, 1.5, 8., "coatingReflectFresnel"),
    roundAperture("App", Point3(-20., 0., 40.), XAXIS, refIndexDefault, 0., 5.; ydir = Vec3(sin(π/2), cos(π/2), 0.)),
    roundAperture("Hole", Point3(-20., -2.5, 40.), XAXIS, refIndexDefault,1., ∞;ydir = nothing),
    rectAperture("RectApp", Point3(-120., -5., 40.), XAXIS,
    Vec3(0., sqrt(2.)/2, sqrt(2.)/2), refIndexDefault, 1.0, 2.0, 5., 10.),
    l2i, l2o,
    planeMirror("Mirror", Point3(-161., 0., 40.), XAXIS, refIndexDefault,
        1.5, 8., "coatingTransmitFresnel"),
    # test going backwards through refractive surfaces
    l2o, l2i,

    referencePlane("end", Point3(50., 0., 40.), XAXIS, refIndexDefault, 3., "testcoat")
]

#
fig6, scene6 = plotGeometry3D(testGeometry2)
cam = Makie.cameracontrols(scene6)
cam.eyeposition[]=SVector{3, Float64}(-80, -300., 100.)
cam.lookat[] = ORIGIN
cam.upvector[] = YAXIS
Makie.update_cam!(scene6.scene)
display(fig6)


function testrace(n)
    sts = zeros(Int64, n)
    - , trc = traceGeometryRel(Ray(Point3(0.0, rand(), 0.0), ZAXIS), testGeometry2)
    for i in 1:n
        sts[i], trc = traceGeometryRel(Ray(Point3(0.0, rand(), 0.0), ZAXIS), testGeometry2)
    end
    sts
end

function testrace2(n)
    sts = zeros(Int64, n)
    - , trc = traceGeometryRel(Ray(Point3(0.0, rand(), 0.0), ZAXIS), testGeometry2)
    for i in 1:n
        sts[i], len = traceGeometryRel!(trc,Ray(Point3(0.0, rand(), 0.0), ZAXIS), testGeometry2)
    end
    sts
end

@time testrace(10000)
@time testrace2(10000)
@time st, trc = traceGeometryRel(Ray(ORIGIN, ZAXIS), testGeometry2)
@time st, len = traceGeometryRel!(trc, Ray(ORIGIN, ZAXIS), testGeometry2)
ottest(n) = [Point(0.0, 0.0, i) for i in 1:n]
## test OAP operation test4

function testOAP(f1, f2, oapseparation, aper)

    vefl = f1 * 25.4 # 1 inch
    vefl1 = f2 * 25.4
#=
    oapseparation = 50.0
    aper = 38.1 * 0.5 #from NewportMKS catalog
=#
    sourceplane = ORIGIN
    mirror1 = Point3(sourceplane .+ (vefl1) .* ZAXIS...)
    println("mirror1 base = $mirror1")
    mirror2 = Point3(mirror1 .+ oapseparation .* YAXIS...)
    println("mirror2 base = $mirror2")
    imageplane = Point3(mirror2 .+ vefl .* ZAXIS...)
    println("image plane = $imageplane")

    c = 1.0/vefl
    c1 = 1.0/vefl1
    ϵ = 0.0 #it's a parabola

    [
        referencePlane("start",sourceplane, ZAXIS, refIndexDefault , 5.0, "testcoat"),

        #=
        reflectConic("bigmirror", mirror1 .- SVector{3, Float64}(0., vefl/2, vefl), -YAXIS,
            refIndexDefault, 1.5,-c, -ϵ, 2*aper, "dummycoating"),
        =#

        reflectOAP("OAP1", mirror1, -YAXIS,  ZAXIS,
                refIndexDefault, 1.5,
                -c, aper, "dummycoating"),

        reflectOAP("OAP2", mirror2, YAXIS, -ZAXIS,
            refIndexDefault, 1.5,
            -c1, aper, "dummycoating"),

        referencePlane("image",imageplane, ZAXIS, refIndexDefault , 5.0, "testcoat")
    ]
end



function offAxisParabolaTest(fig)

    parabolageometry3 = testOAP(1.0, 1.0, 50., 38.1/2)
    outer_padding = 30
    lscene = LScene(fig[1,2], scenekw = (camera = cam3d_cad!, raw = false))
    linesegments!(lscene,[Makie.Point3f(0, 0,0) => Makie.Point3f(0,0,1)],color=:green, linewidth=2)
    linesegments!(lscene,[Makie.Point3f(0, 0,0) => Makie.Point3f(1,0,0)],color=:blue, linewidth=2)
    linesegments!(lscene,[Makie.Point3f(0, 0,0) => Makie.Point3f(0,1,0)],color=:red, linewidth=2)

    plotGeometry3D!(lscene, parabolageometry3)

    testProfile = SurfProfileOAConic(1/25.0 ,  0.0 , Vec3(0.0, 25.0, -25.0 * 0.5))

    #θ = 0.3
    θ = 0.

    trcAndPrintPlotRay!(lscene, Ray(ORIGIN,-ZAXIS),parabolageometry3)
    
    ax1 = fig[1,3] = Axis(fig, title = "Sag")

    plotXSag!(ax1, 25.0, 0.0, testProfile)
    plotYSag!(ax1, 25.0, 0.0, testProfile)

    fig
end

offAxisParabolaTest(figs[4])

display(figs[4])

function sN(x)
    sn = surfNormal(Vec3{3, Float64}(0.0,x, sag(0.0, x, testProfile)),testProfile)
    n = norm(sn)
    println("sn = $sn   norm = $n")
end






## test 5
diffusergeo =
    [

    cDiffuser("diffuser", ORIGIN, ZAXIS, refIndexDefault, refIndexDefault, 25.0*π / 180., 12.5, "nocoating"),
    roundAperture("aperture", Point3(0., 0., 100.), ZAXIS, 1.0, 0.0, 20.0, color = :green3),
    refractSphere("Lens1i", Point3(0., 0., 100.), ZAXIS, 1.0, 101., 0., 20., "testcoat"),
    refractSphere("Lens10", Point3(0., 0., 101.), ZAXIS, 101., 1., -1. / 5000., 20., "testcoat"),
    referencePlane("image",Point3(149.969 .* ZAXIS...), ZAXIS, refIndexDefault , 50., "testcoat")
    ]

fig7, scenediff = plotGeometry3D(diffusergeo)
begin
    for i in 1:100
        trcAndPlotRay!(scenediff, Ray(ORIGIN, Vec3(normalize([0., 0., .9])...)), diffusergeo, color=:blue)
    end
end

display(fig7)

#=


perfectLensgeo = [
    refractSphere("Lens1i", SVector(0., 0., 100.), ZAXIS, 1., 1001., 0., 20., "testcoat"),
    refractSphere("Lens10", SVector(0., 0., 100.1), ZAXIS, 1001., 1., -1. / 50000., 20., "testcoat"),
    referencePlane("image",149.088 .* ZAXIS, ZAXIS, refIndexDefault , 50., "testcoat")

]

perfectLensgeoA = [
    referencePlane("object", ORIGIN, ZAXIS, refIndexDefault , 10., "testcoat")
    perfectLensImage(0., 0., 80., ZAXIS, 80., 8.0, 10.)
]


scenediff = plotGeometry3D(perfectLensgeoA)


trcAndPlotRay!(scenediff, Ray(5. .*YAXIS, ZAXIS), Array{AbstractSurface,1}(perfectLensgeoA), color=plotcolors[1])

trc=trcAndPrintRay(Ray(.01 .*YAXIS, ZAXIS), Array{AbstractSurface,1}(perfectLensgeoA))

trcAndPlotRay!(scenediff, Ray(5. .*YAXIS, ZAXIS), perfectLensgeoA, color=:blue)

plotGeometry3D(perfectLensgeoA)

=#




##
#=
sfunc(x, y) =  x^2 + y^2 < 64 ? (x^2 + y^2)/8 : NaN


function surfacetest2()
    x = -10.:.1:10
    y = -10.:.1:10
    z = [sfunc(xi,yi) for xi in x, yi in y]
    surface(x, y, z, colormap = :Spectral)
end

display(surfacetest2())
=#

## test 6
# Test a 4f system with two TLF357775_405 lenses


function funcfourf()
    bfl = 4.0 # focal length in mm
    λ = 405.0e-3 # wavelength in μm
    ffl = 1.9+0.25/1.57+ 0.25 # estimated front focal lenght
    geo = [
        [referencePlane("object",ORIGIN, ZAXIS, refIndexDefault , 1., "testcoat")]
        lens_TLF357775_405(Point3(0.,0., ffl), ZAXIS,  λ; order = "forward", lensname = "L1")
        [roundAperture("aperture", Point3(0., 0., ffl+thickTL357775_405+bfl), ZAXIS, 1.0, 0.0, 3.5/2, color = :green3)]
        lens_TLF357775_405(Point3(0.,0., ffl+thickTL357775_405+bfl+bfl), ZAXIS,  λ; order = "reverse", lensname = "L2")
        [referencePlane("image",Point3(0., 0., 2*(ffl+thickTL357775_405+bfl)), ZAXIS, refIndexDefault , 1., "testcoat")]
    ]
    return geo
end

fourfgeo = funcfourf()


fig7, scenediff = plotGeometry3D(fourfgeo)
begin
    for x in -0.4:0.1:0.4
        trcAndPlotRay!(scenediff, Ray(ORIGIN, Vec3(x, 0., sqrt(1.0-x^2))), fourfgeo, color=:blue)
    end
end


display(fig7)
##

g = lens_TLAC254_060(ORIGIN, ZAXIS,  0.45; order = "forward", lensname = "TL_AC254-060")
gp = lens_TLAC254_060(ORIGIN, ZAXIS,  0.45; order = "reverse", lensname = "TL_AC254-060")
fig8, scenediff = plotGeometry3D(g)
display(fig8)

rfp,fl = findRFP(g)
rfpp, flp = findRFP(gp)

println("RFP = $rfp   fl = $fl")
println("RPP = $rfpp   fl = $flp")
bfd = rfp.base.base[3] - g[end].base.base[3]
ffd = rfpp.base.base[3] - gp[end].base.base[3]
printGeo([rfp])
printGeo([rfpp])

