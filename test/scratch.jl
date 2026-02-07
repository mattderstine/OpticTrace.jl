4+4


using OpticTrace
using GeometryBasics
using BenchmarkTools
using Test
using StaticArrays
using GLMakie
using ForwardDiff
##
#=

attributesSurfaces = Dict{String, Any}()


=#



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


function readTest(filename::AbstractString)
    lines = readlines(filename)

    # Parse the file to extract the necessary data
    for l in lines
        curlineL =  string(strip(l,['\n','\r',' ', '\t', '\xff', '\xfe', '\0']) ) # Clean line endings and whitespace
        #println("Line: ", curline)'
        curline = replace(curlineL, "\x00" => "")  # Remove null characters
        entries = split(curline)
        #=
        if length(entries) == 0
            println("Blank line: $curline")
            continue
        else
            for e in entries
                print(e, ", ")
            end
            println(" ")
        end
        =#
        #=
        l = length(entries[1])
        if l == 8
            ent = entries[1][1]*entries[1][3]*entries[1][5]*entries[1][7]
        else
            ent = entries[1]
        end
        =#
        if length(entries) != 0

            ent = entries[1]
            if startswith(ent, "SURF")
                println("SURF line: $curline")
            end

        end


    end

 
end


filename = "/Users/matt/Desktop/Zemax lenses/AL5040G-Zemax(ZMX).zmx"
readTest(filename)



##

filename = "/Users/matt/Desktop/Zemax lenses/AL5040G-Zemax(ZMX).zmx"

zgeo, name, units, wave = readZemax(filename; basept = ORIGIN, dir = ZAXIS)

printZemaxSurfs(zgeo)

riN_LAK22 = OpticTrace.getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/schott/N-LAK22.yml")
glassCatalog = loadRICatalog( "glass/schott")
loadRICatalog!(glassCatalog, "glass/hikari")
loadRICatalog!(glassCatalog, "glass/hoya")
loadRICatalog!(glassCatalog, "glass/ohara")
fullCatalog = loadRICatalog( "glass")

riN_LAK22(0.5) == glassCatalog["N-LAK22"](0.5) == glassCatalog["LAK22"](0.5) == fullCatalog["LAK22"](0.5)

glassCatalog["DEFAULT"](0.5)==1

geo = zemaxsurfsToGeo(zgeo[2:end], ORIGIN, ZAXIS, 0.5, glassCatalog)

fig,a = plotGeometry3D(geo)
display(fig)
printGeo(geo)

c = geo[1].profile.curv
ϵ = geo[1].profile.ϵ
evenasp = geo[1].profile.a

allasp = zeros(Float64, 20)
allasp[2]= evenasp[1]
allasp[4]= evenasp[2]
allasp[6]=evenasp[3]

asurf=OpticTrace.OptSurface("old asph", geo[1].base, geo[1].aperture, OpticTrace.SurfProfileAsphere(c, ϵ, allasp), 
        geo[1].mod, geo[1].coating, 
        geo[1].toGlobalCoord,geo[1].toLocalCoord,geo[1].toGlobalDir,geo[1].toLocalDir,
        :blue)

plotGeometry3D!(a, [asurf])

function to_bytes(n::Integer; bigendian=true, len=sizeof(n))
    bytes = Array{UInt8}(undef, len)
    for byte in (bigendian ? (1:len) : reverse(1:len))
        bytes[byte] = n & 0xff
        n >>= 8
    end
    return bytes
end

ZAR = ".zar"
ZIP = ".zip"
ZAR_VERSION_LENGTH = 2  # in bytes
EARLIER_CONTENT_OFFSET = 0x14C - ZAR_VERSION_LENGTH
EARLIER_PACKED_FILE_SIZE_BEGIN = 0xC - ZAR_VERSION_LENGTH
EARLIER_PACKED_FILE_SIZE_END = 0x10 - ZAR_VERSION_LENGTH
EARLIER_PACKED_FILE_NAME_OFFSET = 0x20 - ZAR_VERSION_LENGTH
EARLIER_VERSION = to_bytes(0xEA00,len=2)
LATEST_VERSION = to_bytes(0xEC03, len = 2)
LATEST_CONTENT_OFFSET = 0x288 - ZAR_VERSION_LENGTH
LATEST_PACKED_FILE_SIZE_BEGIN = 0x10 - ZAR_VERSION_LENGTH
LATEST_PACKED_FILE_SIZE_END = 0x18 - ZAR_VERSION_LENGTH
LATEST_PACKED_FILE_NAME_OFFSET = 0x30 - ZAR_VERSION_LENGTH

##

rIN_MgF2 = getRefractiveIndexFunc("main/MgF2/Li-o.yml")
rIN_MgF2e = getRefractiveIndexFunc("main/MgF2/Li-e.yml")
rIN_MgF2d = getRefractiveIndexFunc("main/MgF2/Dodge-o.yml")
rIN_MgF2de = getRefractiveIndexFunc("main/MgF2/Dodge-e.yml")

rIN_MgF2(0.2)
rIN_MgF2(1.)

x = LinRange(0.12, 1.0, 200)
y = rIN_MgF2.(x)

f = Figure()
ax = Axis(f[1, 1],
    title = "MgF2",
    xlabel = "Wavelength (µm)",
    ylabel = "Refractive Index",
)
lines!(ax, x, rIN_MgF2.(x), label = "Li o")
#lines!(x, rIN_MgF2e.(x), color=:red)
lines!(ax, x, rIN_MgF2d.(x), color=:green, label="Dodge o")
#lines!(x, rIN_MgF2de.(x), color=:orange)
axislegend(ax, position = :rt)
rIN_MgF2(0.115) - rIN_MgF2d(0.115)

f1(n,r) = r/(n-1)
g1(lambda, r) = f1(rIN_MgF2(lambda), r)


rmirror = 10.0
radlens = 200.0

frac(.12)

g0 = g1(.12, radlens)
u = rmirror/g1(.12, radlens)

na = rmirror/g0



function filter(rmirror, radlenslist)
    x = LinRange(0.12, 1.0, 200)
    fig =Figure()
    ax = Axis(fig[1, 1],
        xlabel = "Wavelength (µm)",
        ylabel = "Focal Length",
    )

    ax2 = Axis(fig[1, 2],
        xlabel = "Wavelength (µm)",
        ylabel = "log(T)",
    )


    y1(lambda, g0, radlens) = rmirror-(rmirror/g1(lambda, radlens))* g0
    r1(lambda, g0, radlens) = y1(lambda, g0, radlens) + 0.02
    frac(lambda, g0, radlens) = (r1(0.12, g0, radlens)/r1(lambda, g0, radlens))^2
    for radlens in radlenslist
        g0 = g1(.12, radlens)
        println("Rad lens: $radlens, g0: $g0")
        na = rmirror/g0
        lines!(ax, x, g1.(x, radlens))
        lines!(ax2, x, log.(frac.(x, g0, radlens)),label = "R=$radlens")
    end
    axislegend(ax2, position = :rt)
    fig
end



filter(1.0, [10.0, 50.0, 100., 200.0, 300.0, 400.0, 500.0])


##

geosimple = [roundAperture("pupil", Point3(0., 0., 10.), ZAXIS, 1.0, 0.0, 5.0, color = :green3)]

pts, center, rmsradius, deltaZ = spotDiagramHex(geosimple, Point3(0., 0., 0.), geosimple[1],9, bestFocus = false)
f = Figure(;size=(1200, 700))
plotGeometry3D(f[1,1], geosimple)
plotSpotDiagram(f[1,2], pts, center, rmsradius, deltaZ; title="Spot Diagram")
display(f)

println("Centroid: $center RMS radius: $rmsradius")

##

spsAsphere = OpticTrace.SurfProfileAsphere(0.0, 0.0, [0.00, 0.1, 0.0, 0.01])
spsEAsphere = OpticTrace.SurfProfileEvenAsphere(0.0, 0.0, [0.1, 0.01])


"""
    define function to compute the normal vector using the gradient of the sag function
"""
function normal_from_sag(x1, y1, surfProfile)
    f(x) = sag(x[1], x[2], surfProfile)

    x0 = [x1, y1]
    # find the gradient of the sag value at the intersection point using an autodiff package
    grad_sag = ForwardDiff.gradient(f, x0)
    # the normal vector is then obtained by normalizing the gradient vector
    normal_vector = normalize(Vec3(grad_sag[1], grad_sag[2], -1.0))
    return normal_vector
end
normal_from_sag(1.0, 1.0, spsAsphere)






        sag_value = sag(1.0, 1.0, spsAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2

        sag_value = sag(1.0, 1.0, spsEAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2
        sag_value = sag(1.0, 1.0, spsAsphere)
        @test sag_value ≈ asphere_factor1 * asphere_coeff1 + asphere_factor2 * asphere_coeff2

        sag_value = sag(1.0, 1.0, spsEAsphere)
        @test sag_value ≈ asphere_factor1 * asphere_coeff1 + asphere_factor2 * asphere_coeff2

        @test normal_from_sag(1.0, 1.0, spsAsphere) ≈ OpticTrace.surfNormal(Point(1.0, 1.0, sag(1.0, 1.0, spsAsphere)), spsAsphere)



x0 = [1.0, 1.0]

sag_value = sag2(x0[1], x0[2], spsAsphere)
function g(x) 
    return sag2(x[1], x[2], spsAsphere)
end
grad_sag = ForwardDiff.gradient(g, x0)
grad = OpticTrace.normalize(Vec3(grad_sag[1], grad_sag[2], -1.0))
println(grad)

function sag2(x::T, y::T, s::OpticTrace.SurfProfileAsphere) where T<:Real
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

normal = OpticTrace.surfNormal(Point3(x0[1], x0[2], sag2(x0[1], x0[2], spsAsphere)), spsAsphere)