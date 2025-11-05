4+4

using OpticTrace
using GeometryBasics
using BenchmarkTools
using Test
using StaticArrays
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
