4+4


using OpticTrace
using GeometryBasics
using BenchmarkTools
using Test
using StaticArrays
using GLMakie
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




pts, center, rmsradius= spotDiagramHex(geosimple, Point3(0., 0., 0.), geosimple[1],9)
f,a,p =scatter(pts)
a.aspect = DataAspect()
display(f)

println("Centroid: $center RMS radius: $rmsradius")

x = [p[1] for p in pts]
y = [p[2] for p in pts]
f,a,p = scatter(x, y)
display(f)

points = Point2.(x,y)
f,a,p = scatter(points)
display(f)

fig9, scene = plotGeometry3D(geosimple)
fig9


ax = fig9[1,2] = Axis(fig9, title = "Spot Diagram")
scatter!(ax, pts, markersize=2, color=:red)
ax.aspect = DataAspect()
display(fig9)


