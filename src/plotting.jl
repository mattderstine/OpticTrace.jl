export saveFigure, printFigure, multipleFigures,plotSurface3D!,plotGeometry3D
export plotGeometry3D!
export trcAndPrintPlot!, trcAndPrintPlotRay!, trcAndPlotRay!,trcAndPlotRayRel!
export plotRayFan!,perimeterRays,plotPerimeterRays
export plotPerimeterRays!, rayHeatmap, rayHeatmap!, computeExitPupilLoc
export plotOPD!, plotOPD3D!
export plotXSag!, plotYSag!
export plotTrace!
export plotSpotDiagram




"""
    saveFigure(fileNameStub, fig; startnum = 0, directory = "")
    saves Makie figure to a file, then displays it
    returns nothing (the return value of the trailing `display(fig)` call)

    Eventually will check for the exisitence of the file and add increasing numbers to avoid overwrite
    Currently just writes to fileNameStub*".png" in the current directory, overwriting previous files

    fileNameStub    string with base name for files to be stored
    fig             Makie figure object
    startnum        number to append to the stub. If 0, nothing is appended
    directory       location of where to store the figure image

"""
function saveFigure(fileNameStub, fig; startnum = 0, directory = "")
    if startnum > 0
        fileNameStub = fileNameStub * string(startnum)
    end
    filename = joinpath(directory, fileNameStub * ".png")
    save(filename, fig)
    display(fig)
end

#try to depricate this name
"""
    printFigure(fileNameStub, fig; startnum = 0, directory = "")

Deprecated alias for `saveFigure` -- prints a deprecation notice, then
delegates to `saveFigure(fileNameStub, fig; startnum, directory)`. See
that method's docstring above for the full parameter/return contract.
"""
function printFigure(fileNameStub, fig; startnum = 0, directory = "")
    println("printFigure is depricated. Use saveFigure instead.")
    saveFigure(fileNameStub, fig; startnum, directory)
end


"""
    multipleFigures(numFigs = 5; fileNameStub = "FigurePrint", activeButtonColor = RGBf(0.8, 0.94, 0.8), activeFig = 1,  size = (1200,900))
    sets up a stack of figures with buttons to navigate between them. Includes a button to save the figure image to a file
    returns an array of Makie figures
"""
function multipleFigures(numFigs = 5; fileNameStub = "FigurePrint", activeButtonColor = RGBf(0.8, 0.94, 0.8), activeFig = 1,  size = (1200,900))
    GLMakie.activate!(title="lens Julia", inline=false)   #this only works in GLMakie
    figs = [Figure(;size,backgroundcolor = RGBf(0.98, 0.98, 0.98)) for i in 1:numFigs]
    for (i, fig) in enumerate(figs)
        #ax1 = fig[1,1] = Axis(fig, title="Figure $i")
        fig[1, 1] = buttongrid = GridLayout(tellwidth = true, tellheight=false)
        printButtonNum = numFigs + 1
        buttons = buttongrid[1:printButtonNum, 1] = [Button(fig, label = (i<=numFigs ? "$i" : "Save")) for i in 1:printButtonNum]

        printButton = buttons[printButtonNum]

        saveColor = buttons[1].buttoncolor[]


        for (j,button) in enumerate(buttons)
            if i==j
                button.buttoncolor[] = activeButtonColor
            end

            on(button.clicks) do n
                #println("$(button.label[]) was clicked $n times.")
                if button == printButton
                    printFigure(fileNameStub, fig)
                else
                    display(figs[findfirst(isequal(button), buttons)])
                end
            end
        end
    end
    display(figs[activeFig])
    figs
end








#plotSurface3D!(scene, s::OptSurface; color=:aquamarine2)=Makie.mesh!(scene, s, color=color, fxaa = true, transparency = true)#someday add options to display mesh
"""
    plotSurface3D!(scene, s::OptSurface; transparency=false)

Render `s` into `scene` as a Makie mesh, colored `s.color`
(`Makie.mesh!(scene, s, color=s.color, fxaa=true, transparency)`).
`s`'s own `GeometryBasics` overloads (`src/mesh_primitives.jl`) supply
the mesh geometry. See `plotSurface3D!(scene, s::ModelSurface)` below
for the non-`OptSurface` sibling method (it dispatches differently, via
`plotModelSurf!`, rather than a direct `Makie.mesh!` call).
"""
plotSurface3D!(scene, s::OptSurface; transparency = false)=Makie.mesh!(scene, s, color=s.color, fxaa = true, transparency = transparency)#someday add options to display mesh

"""
    plotGeometry3D!(scene, geometry)

Render every surface in `geometry` into `scene` via `plotSurface3D!`.
Returns `scene`. See `plotGeometry3D(fig, geometry)`'s docstring above
for the higher-level entry point that also sets up the `LScene`/axes
this is meant to be plotted into.
"""
function plotGeometry3D!(scene, geometry)
    for geo in geometry
        plotSurface3D!(scene, geo)
    end
    scene
end


"""
plotGeometry3D(geometry; size = (1200,700))
geometry is an array of OptSurfaces
size is the initial window size
creates a new Figure and delegates to plotGeometry3D(fig, geometry)
returns a Makie Figure and LScene located at [1,1]

"""
function plotGeometry3D(geometry; size = (1200,700))
    fig = Figure(;size)
    plotGeometry3D(fig, geometry)
end

"""
plotGeometry3D(fig, geometry)
geometry is an array of OptSurfaces
fig is an existing Makie Figure to plot into
draws coordinate axes and the geometry into an LScene at fig[1,1]
returns the Makie Figure and LScene located at [1,1]

"""
function plotGeometry3D(fig, geometry)
    ax = fig[1,1]=LScene(fig, show_axis=false, scenekw = (camera = cam3d_cad!,))
    linesegments!(ax,[Point3f(0, 0,0) => Point3f(0,0,1)],color=:green, linewidth=2)
    linesegments!(ax,[Point3f(0, 0,0) => Point3f(1,0,0)],color=:blue, linewidth=2)
    linesegments!(ax,[Point3f(0, 0,0) => Point3f(0,1,0)],color=:red, linewidth=2)
    plotGeometry3D!(ax, geometry)
    fig,ax
end

"""
    plotSurface3D!(scene, s::ModelSurface)

Render `s` into `scene` by dispatching on its aperture type to
`plotModelSurf!(scene, s.aperture, s)` (see that function's
`RectAperture`/`RoundAperture` methods below) -- `ModelSurface`s have
no mesh profile of their own to draw (see `NoProfile`), so this draws
their aperture shape instead. Compare `plotSurface3D!(scene,
s::OptSurface; transparency=false)` above, which draws the surface's
actual profile mesh directly.
"""
function plotSurface3D!(scene, s::ModelSurface)
    plotModelSurf!(scene, s.aperture, s) #dispatch on s.aperture
end

"""
    plotModelSurf!(scene, a::RectAperture, s::ModelSurface, dsize=0.3)

Render `s`'s rectangular aperture into `scene`: an obscuration
rectangle (if `a.wo`/`a.lo` are nonzero), and a "frame" of four
polygons around the clear-aperture rectangle (if `a.wclear`/`a.lclear`
are both finite; `dsize` sets the frame's width relative to the
aperture size). Draws nothing for whichever piece doesn't apply. See
`plotModelSurf!(scene, a::RoundAperture, s::ModelSurface, dsize=0.3)`
below for the round-aperture sibling method.
"""
function plotModelSurf!(scene, a::RectAperture, s::ModelSurface, dsize = 0.3)
    #this is brute force. I'm still learning...
    connect = [
    1 2 3;
    3 4 1
    ]

    if a.wo !=0f0 && a.lo !=0f0
        x0 = a.wo
        y0 = a.lo
        r1 = [-x0 -y0 0.0;
            -x0 y0 0.0;
            x0 y0 0.0;
            x0 -y0 0.0]
        r1a = [s.toGlobalCoord(r1[i,1:3]) for i in 1:4]
        Makie.poly!(scene, r1a, connect, color=s.color, fxaa = true, transparency = true)

    end

    if a.wclear != ∞ && a.lclear != ∞
        x0 = a.wclear
        y0 = a.lclear
        dx = dsize * 0.5(x0+y0) #uniform "window" frame
        dy = dx
        x1 = -x0 - dx
        x2 = -x0
        x3 = x0
        x4 = x0 + dx
        y1 = -y0 - dy
        y2 = -y0
        y3 = y0
        y4 = y0 + dy
        #could be cleaned up to one polygon, but not taking the time
        r1 = [
            x1 y1 0.;
            x1 y2 0.;
            x3 y2 0.;
            x3 y1 0.]
        r2 =[
            x1 y2 0.;
            x1 y4 0.;
            x2 y4 0.;
            x2 y2 0.]
        r3 = [
            x2 y4 0.;
            x4 y4 0.;
            x4 y3 0.;
            x2 y3 0.]
        r4 = [
            x4 y3 0.;
            x4 y1 0.;
            x3 y1 0.;
            x3 y3 0.]
        r1a = [s.toGlobalCoord(r1[i,1:3]) for i in 1:4] #would like to know how to do this without hard coding length
        r2a = [s.toGlobalCoord(r2[i,1:3]) for i in 1:4]
        r3a = [s.toGlobalCoord(r3[i,1:3]) for i in 1:4]
        r4a = [s.toGlobalCoord(r4[i,1:3]) for i in 1:4]
        Makie.poly!(scene, r1a, connect, color=s.color, fxaa = true, transparency = true)
        Makie.poly!(scene, r2a, connect, color=s.color, fxaa = true, transparency = true)
        Makie.poly!(scene, r3a, connect, color=s.color, fxaa = true, transparency = true)
        Makie.poly!(scene, r4a, connect, color=s.color, fxaa = true, transparency = true)
    end
end

#Makie.mesh!(scene, s, color = :orange,  transparency=true, camera=Makie.cam3d_cad!)#someday add options to display mesh

"""
    plotModelSurf!(scene, a::RoundAperture, s::ModelSurface, dsize=0.3)

Render `s`'s round aperture into `scene`: a `Disk` mesh for the central
obscuration (if `a.obscure != 0`), and a `Washer` mesh for the outer
clear-aperture ring (if `a.semiDiameter != ∞`). `dsize` is accepted for
signature symmetry with the `RectAperture` method above but is unused
here (`Disk`/`Washer` don't take a frame-width parameter). See
`Washer`/`Disk`'s docstrings (`src/mesh_primitives.jl`) for their known
`GeometryBasics.radius`/`widths` bugs, which this rendering path can
trigger (see `TODO.md`).
"""
function plotModelSurf!(scene, a::RoundAperture, s::ModelSurface, dsize = 0.3)
    if a.obscure != 0.
        Makie.mesh!(scene, Disk(s.base.base, s.base.dir,
            Float64(a.obscure),s.toGlobalCoord,s.toGlobalDir),
            color = s.color, fxaa = true, transparency = true)# center obscuration
    end
    if a.semiDiameter != ∞ # not just an obscuration
        Makie.mesh!(scene, Washer(s.base.base, s.base.dir,
            Float64(a.semiDiameter),s.toGlobalCoord,s.toGlobalDir),
            color = s.color, fxaa = true, transparency = true)
    end
end


"""
trcAndPrintPlot!(ray::Ray, geo; color=:blue)
    trace an absolute ray
    print the ray
    plot the ray on the last LScene
"""
function trcAndPrintPlot!(ray::Ray, geo; color=:blue)

    trc=trcAndPrintRay(ray, geo)
    Makie.lines!([a.ray.base for a in trc], color=color)
    trc
end

"""
trcAndPrintPlotRay!(scene, ray::Ray, geo; color=:blue)
    trace an absolute ray
    print the trace
    plot the ray on scene
    return the trace
"""
function trcAndPrintPlotRay!(scene, ray::Ray, geo; color=:blue)
    trc=trcAndPrintRay(ray, geo)
    Makie.lines!(scene, [a.ray.base for a in trc], color=color)
    trc
end

"""
trcAndPlotRay!(scene, ray::Ray, geo; color=:blue, clipmsg=false)
     trace an absolute ray
    plot the ray on scene
    return the status & trace
    set clipmsg=true to print a message if the ray does not reach the end of geo

"""
function trcAndPlotRay!(scene, ray::Ray, geo; color=:blue, clipmsg=false)
    status, trc = traceGeometry(ray, geo)
    printTrcStatus(status, flagNormal = false, clipmsg=clipmsg)
    Makie.lines!(scene, [a.ray.base for a in trc], color=color)
    status, trc
end

"""
    plotTrace!(scene, trace::Vector{Trace}; color=:blue)

Plot an already-computed `trace` (e.g. from `traceGeometry`/
`traceGeometryRel`) into `scene`, connecting each step's ray base
point. Returns `scene`. Unlike `trcAndPlotRay!`/`trcAndPlotRayRel!`
above/below, this doesn't trace a ray itself -- it just draws one you
already have.
"""
function plotTrace!(scene, trace::Vector{Trace}; color=:blue)
    Makie.lines!(scene, [a.ray.base for a in trace], color=color)
    scene
end

#=
see if this method is needed

function trcAndPlotRay(ray::Ray, geo; color=:blue, clipmsg=false)
    trcStatMsg=("Normal","Missed","TIR","Clipped")
    status, trc = traceGeometry(ray, geo)
    if clipmsg && status != 0
        println(trcStatMsg[status+1])
    end
    Makie.lines!(scene, [a.ray.base for a in trc], color=color)
end
=#
"""
trcAndPlotRayRel!(scene, ray::Ray, geo; color=:blue, clipmsg=false)
    trace a relative ray
    plot the ray on scene
    return the status & trace
    set clipmsg=true to print a message if the ray does not reach the end of geo

"""
function trcAndPlotRayRel!(scene, ray::Ray, geo; color=:blue, clipmsg=false)
    status, trc = traceGeometryRel(ray, geo)
    printTrcStatus(status, flagNormal = false, clipmsg=clipmsg)
    points = [Point{3,Float32}(a.ray.base...) for a in trc]
    lines!(scene, points, color=color)
    status, trc
end

"""
trcAndPlotRayRel!(ray::Ray, geo; color=:blue, clipmsg=false)
    trace a relative ray
    plot the ray on the active scene
    return the status & trace
    set clipmsg=true to print a message if the ray does not reach the end of geo

"""
function trcAndPlotRayRel!(ray::Ray, geo; color=:blue, clipmsg=false)
    status, trc = traceGeometryRel(ray, geo)
    printTrcStatus(status, flagNormal = false, clipmsg=clipmsg)
    Makie.lines!([Point{3,Float32}(a.ray.base...) for a in trc], color=color)
    status, trc
end


"""
    opdRel(ray, refTrace, geo)
    compute the OPD relative to refTrace of relative ray with optics geo

"""
function opdRel(ray, refTrace, geo)
    status, trc = traceGeometryRel(ray, geo)
    if status >0
        return NaN
    end
    opd = 0.
    for (t, r) in zip(trc, refTrace)
        opd += (t.delta * t.nIn- r.delta * r.nIn)
    end
    opd
end






#=
"""
    rayfan(point, min angle, max angle, points, geometry)
    returns a set of rays intercepts at the image plane of the geometry. scan is in y angle
"""
function rayfan(r::SVector{3, Float64}, θmin::Float64, θmax::Float64, pnts::Int64, geo;surfview = "end")
    trcStatMsg=("Normal","Missed","TIR","Clipped")
    raysy = Vector{Ray}(undef, pnts)
    raysx = Vector{Ray}(undef, pnts)
    surfnum = surfnumFromName(surfview, geo)
    i = 1
    status, trc = traceGeometry(Ray(r, [0.,0.,0.]), geo)
    if status != 0
        #println(trcStatMsg[status+1])
        println("Reference ray missed ")
        return [Ray([NaN, NaN, NaN], [NaN, NaN, NaN])],[Ray([NaN, NaN, NaN], [NaN, NaN, NaN])]
    else
        refbase = surfnum <0 ? trc[end].ray : trc[surfnum].ray
    end

    for θ in LinRange(θmin, θmax, pnts)
        status, trc = traceGeometry(Ray(r, [0.,sin(θ), cos(θ)]), geo)
        if status != 0
            #println(trcStatMsg[status+1])
            raysy[i] = Ray([NaN, NaN, NaN], [NaN, NaN, NaN])
        else
            raysy[i] = surfnum <0 ? trc[end].ray : trc[surfnum].ray
        end
        i += 1
    end
    i = 1
    for θ in LinRange(θmin, θmax, pnts)
        status, trc = traceGeometry(Ray(r, [sin(θ), 0., cos(θ)]), geo)
        if status != 0
            #println(trcStatMsg[status+1])
            raysx[i] = Ray([NaN, NaN, NaN], [NaN, NaN, NaN])
        else
            raysx[i] = surfnum <0 ? trc[end].ray : trc[surfnum].ray
        end
        i += 1
    end

    raysy, raysx

end
=#

"""
    getrefbase(point, geo, surfview)
    returns the local coodinates of the reference ray intersection at surfview

"""

function getrefbase(r::Point3, geo, surfview)
    surfnum = tracenumFromName(surfview, geo)
    s = geo[surfnum-1]  #surfnum is defined for traces

    #find mapping of local coordinates
    status, trc = traceGeometryRel(Ray(r, ZAXIS), geo) #ZAXIS is default local direction
    if status != 0 && surfnum >0 && length(trc)<surfnum
        #println(trcStatMsg[status+1])
        error("Ray did not intersect surface: $surfview")
    else
        reforigin = s.toLocalCoord(surfnum <0 ? trc[end].ray.base : trc[surfnum].ray.base)
    end
    reforigin
end


"""
    SKEWLIMIT

Threshold on `cos(Δθ)` between the local x/y reference directions used
by `plotRayFan!`: if the two directions' dot product exceeds this, a
"significant skew distortion" warning is printed, since `plotRayFan!`'s
y/x decomposition assumes they're close to orthogonal.
"""
const SKEWLIMIT = 0.1  #cos(θ) limit for skew rays in rayfan calculation

"""
    plotRayFan!(fig, point, max angle, geometry; color=:blue ,surfview=Surface Name, points=33, θmin = min angle )
    plots the y- and x-fan ray-intercept curves into fig
    returns refbase, the reference intercept point the fan is measured relative to

    current version assumes telecentric pupil/stop (i.e. reference ray θ=0)

"""
function plotRayFan!(fig, r::Point3, θmax::Float64, geo; surfview = "end", color = :blue, points=33, θmin::Float64=NaN, printTrace=false, rayscene=nothing )
    raysy = Vector{SVector{3, Float64}}(undef, points)
    raysx = Vector{SVector{3, Float64}}(undef, points)
    surfnum = tracenumFromName(surfview, geo)
    s = geo[surfnum-1]  #surfnum is defined for traces

    #find mapping of local coordinates

    reforigin = getrefbase(ORIGIN, geo, surfview)
    refxdir = normalize(getrefbase(Point3(0.1, 0.0, 0.0), geo, surfview) - reforigin)
    refydir = normalize(getrefbase(Point3(0.0, 0.1, 0.0), geo, surfview) - reforigin)
    if refxdir ⋅ refydir > SKEWLIMIT
        println("signficant skew distortion in rayfan calculation: cos(Δθ) = $(refxdir ⋅ refydir)")
    #=
    else
        println("refxdir = $refxdir   refydir = $refydir   cos(Δθ) = $(refxdir ⋅ refydir)")
        println("θx = $(s.toLocalDir(refxdir)⋅XAXIS)   θy = $(s.toLocalDir(refydir)⋅YAXIS)")
    =#
    end


    #find the reference local reference intercept coordinates
    refbase = getrefbase(r, geo, surfview)

    θm = isnan(θmin) ? -θmax : θmin
    θr =  LinRange(θm, θmax, points)
    println("max angle = $θmax  max NA = $(sin(θmax))")
    i = 1
    for θ in θr
        #println("dir = $([0.,sin(θ), cos(θ)])")
        status, trc = traceGeometryRel(Ray(r, Vec3(0.,sin(θ), cos(θ))), geo)
        if printTrace
            println("\ny $θ   -----------------------")
            printTrcCoords(status, trc, geoOpticChannelLensOnly(λ0); format="normal")
        end
        if !isnothing(rayscene)
            Makie.lines!(rayscene, [a.ray.base for a in trc], color=color)
        end
        if status != 0 && surfnum >0 && length(trc)<surfnum
            #println(trcStatMsg[status+1])
            raysy[i] = SVector(NaN, NaN, NaN)
        else
            raysy[i] = s.toLocalCoord(surfnum <0 ? trc[end].ray.base : trc[surfnum].ray.base) - refbase
        end

        status, trc = traceGeometryRel(Ray(r, Vec3(sin(θ), 0., cos(θ))), geo)
        if printTrace
            println("\nx $θ   -----------------------")
            printTrcCoords(status, trc, geoOpticChannelLensOnly(λ0); format="normal")
        end
        if !isnothing(rayscene)
            Makie.lines!(rayscene, [a.ray.base for a in trc], color=color)
        end
        if status != 0 && surfnum >0 && length(trc)<surfnum
            #println(trcStatMsg[status+1])
            raysx[i] = SVector(NaN, NaN, NaN)
        else
            raysx[i] = s.toLocalCoord(surfnum <0 ? trc[end].ray.base : trc[surfnum].ray.base) - refbase
        end
        i += 1
    end

    #=
    the right way to do this is to either:
    a) assume the final surface has the yaxis aligned properly to the input rays <- use this one!
    b) find the axis directions from the intersection points of the x&y rays.

    =#
    #find the largest of x & y at the edge and use it for the sign (this gives the wrong answer in some situations)
    #yv = abs(raysy[end][1]) > abs(raysy[end][2]) ? 1 : 2
    #xv = abs(raysx[end][1]) > abs(raysx[end][2]) ? 1 : 2
    #println("yv = $yv   xv = $xv")
    #y = [sign(t[yv])*norm(t) for t in raysy] #use norm in case the output intersections are not on an axis
    #x = [sign(t[xv])*norm(t) for t in raysx]

    y = [t[2] for t in raysy]
    x = [t[1] for t in raysx]

    a,p= Makie.lines(fig, θr, y, color=color, linestyle = :solid)
    Makie.lines!(a, θr, x , color=color, linestyle = :dash)
    refbase
end


#=
"""
    plotRayfan(point, min angle, max angle, points, geometry; color=color surfview=Surface Name )
    returns Makie scene of plot of rayfan
"""
function plotRayfan(r::SVector{3, Float64}, θmin::Float64, θmax::Float64, pnts::Int64, geo; color = :blue, surfview = "end")
    raysy, raysx = rayfan(r, θmin, θmax, pnts, geo, surfview = surfview)

#    Makie.lines(θrange, rays)
    θ = LinRange(θmin, θmax, pnts)
    y = [norm(t.base) for t in raysy]
    x = [norm(t.base) for t in raysx]
    Makie.lines(θ, y , color=color)
    Makie.lines!(θ, x, color=color, linestyle = :dash)
end
"""
    plotRayfan!(scene, point, min angle, max angle, points, geometry; color=color, surfview=Surface Name)
    adds rayfan to scene
"""
function plotRayfan!(scene, r::SVector{3, Float64}, θmin::Float64, θmax::Float64, pnts::Int64, geo::Array{AbstractSurface}; color = :blue, dim=2, surfview = "end")

    raysy, raysx = rayfan(r, θmin, θmax, pnts, geo, surfview = surfview)
    θ = LinRange(θmin, θmax, pnts)
    y = [norm(t.base) for t in raysy]
    x = [norm(t.base) for t in raysx]
    Makie.lines!(scene, θ, y , color=color)
    Makie.lines!(scene, θ, x, color=color, linestyle = :dash)

end

function plotRayfan!(r::SVector{3, Float64}, θmin::Float64, θmax::Float64, pnts::Int64, geo; color = :blue, dim=2, surfview = "end")
    raysy, raysx = rayfan(r, θmin, θmax, pnts, geo, surfview = surfview)
    θ = LinRange(θmin, θmax, pnts)
    y = [norm(t.base) for t in raysy]
    x = [norm(t.base) for t in raysx]
    Makie.lines!(θ, y , color=color)
    Makie.lines!(θ, x, color=color, linestyle = :dash)
end

=#

"""
    perimeterRays(r::SVector{3, Float64}, radius::Float64, θ::Float64, points::Int64, geo; surfview="end")

Trace `points` rays arranged around a circle of `radius` centered at
`r` (in the local x/y plane), all launched at polar angle `θ` from the
local z axis, and collect each one's ray at `surfview` (see
`tracenumFromName`). Intended to return a `Vector{Ray}` of length
`points`.

**This method has a bug**: if any ray fails to trace (`status != 0`),
it prints a status message, stores a NaN `Ray` in `rays[i]`, and then
executes a bare `return` -- which returns `nothing`, discarding the
`rays` vector entirely (including any rays already successfully traced
before this point), rather than continuing to the next perimeter angle
or returning the partially-`NaN`-filled vector. Very likely a
`continue` was intended instead of `return`. This also means
`plotPerimeterRays`/`plotPerimeterRays!` below will fail if fed a `geo`
where any perimeter ray misses. See `TODO.md`.
"""
function perimeterRays(r::SVector{3, Float64}, radius::Float64, θ::Float64, points::Int64, geo;surfview = "end")
    trcStatMsg=("Normal","Missed","TIR","Clipped")
    rays = Vector{Ray}(undef, points)
    surfnum = tracenumFromName(surfview, geo)

    i = 1
    for ϕ in LinRange(0., 2pi, points)
        status, trc = traceGeometryRel(Ray(Point3(r[1] + radius*cos(ϕ), r[2] + radius*sin(ϕ), r[3]),
            Vec3(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ))), geo)
        if status != 0
            println(trcStatMsg[status+1])
            rays[i] = Ray(Point3(NaN, NaN, NaN), Vec3(NaN, NaN, NaN))

            return
        end
        rays[i] = trc[surfnum].ray
        i += 1
    end
    rays
end

"""
    plotPerimeterRays(r::SVector{3, Float64}, radius::Float64, θ::Float64, pnts::Int64, geo; color=:blue, surfview="end")

Trace a ring of perimeter rays via `perimeterRays` and plot them as
arrows (base point + direction) in a new figure. See
`plotPerimeterRays!`/`plotPerimeterRays!(scene, ...)` below for the
mutating variants (into the active scene, or into a given `scene`).
Subject to the `perimeterRays` bug noted above if any ray misses.
"""
function plotPerimeterRays(r::SVector{3, Float64}, radius::Float64, θ::Float64, pnts::Int64, geo; color=:blue, surfview = "end")
    pr = perimeterRays(r, radius, θ, pnts, geo, surfview = surfview)
    Makie.arrows([Makie.Point3f(a.base) for a in pr], [Makie.Point3f(a.dir) for a in pr], linecolor=color, arrowcolor=color, arrowsize=0.1)
end

"""
    plotPerimeterRays!(r::SVector{3, Float64}, radius::Float64, θ::Float64, pnts::Int64, geo; color=:blue, surfview="end")

Like `plotPerimeterRays` above, but plots into the active scene
(`Makie.arrows!` with no explicit scene) instead of creating a new
figure.
"""
function plotPerimeterRays!(r::SVector{3, Float64}, radius::Float64, θ::Float64, pnts::Int64, geo; color=:blue, surfview = "end")
    pr = perimeterRays(r, radius, θ, pnts, geo, surfview = surfview)
    Makie.arrows!([Makie.Point3f(a.base) for a in pr], [Makie.Point3f(a.dir) for a in pr], linecolor=color, arrowcolor=color, arrowsize=0.1)
end

"""
    plotPerimeterRays!(scene, r::SVector{3, Float64}, radius::Float64, θ::Float64, pnts::Int64, geo; color=:blue, surfview="end")

Like `plotPerimeterRays` above, but plots into the given `scene`.
"""
function plotPerimeterRays!(scene, r::SVector{3, Float64}, radius::Float64, θ::Float64, pnts::Int64, geo; color=:blue, surfview = "end")
    pr = perimeterRays(r, radius, θ, pnts, geo, surfview = surfview)
    Makie.arrows!(scene, [Makie.Point3f(a.base) for a in pr], [Makie.Point3f(a.dir) for a in pr], linecolor=color, arrowcolor=color, arrowsize=0.1)
end


"""
    rayHeatmap(pts; mcbins::Int64=25, center=(0.,0.), width=25.)

Bin `pts` (an iterable of 2D points, e.g. ray intercepts) into a
`mcbins`x`mcbins` 2D histogram over a `width`x`width` window centered
at `center`, and plot it as a new Makie heatmap figure.

Returns `(d, fig)`: `d` is the `StatsBase.Histogram`, `fig` is the
Makie heatmap plot object. See `rayHeatmap!` below for the mutating
variant (into an existing axis).
"""
function rayHeatmap(pts; mcbins::Int64=25, center = (0., 0.), width = 25.)
    cprime=transpose(hcat(pts...))
    xmin = center[1] - width/2
    xmax = center[1] + width/2+width/10000
    ymin = center[2] - width/2
    ymax = center[2] + width/2+width/10000
    step = width/mcbins
    #println("\nrayHeatmap\nxmin = $xmin  xmax = $xmax  ymin = $ymin  ymax = $ymax  step = $step\n")

    d = fit(Histogram, (cprime[:,1], cprime[:,2]), (xmin:step:xmax, ymin:step:ymax))
    (d,heatmap(d.edges[1], d.edges[2], d.weights))
end

"""
    rayHeatmap!(ax, pts; mcbins::Int64=25, center=(0.,0.), width=25.)

Like `rayHeatmap` above, but plots into the existing axis `ax`
(`Makie.heatmap!`) instead of creating a new figure.

Returns `(d, plt)`: `d` is the `StatsBase.Histogram`, `plt` is the
Makie heatmap plot object.
"""
function rayHeatmap!(ax, pts; mcbins::Int64=25, center = (0., 0.), width = 25.)
    cprime=transpose(hcat(pts...))
    xmin = center[1] - width/2
    xmax = center[1] + width/2+width/10000
    ymin = center[2] - width/2
    ymax = center[2] + width/2+width/10000
    step = width/mcbins
    #println("\nrayHeatmap!  xmin = $xmin xmax = $xmax ymin = $ymin ymax = $ymax")
    #println("step = $step   xmax-xim = $(xmax-xmin) ymax-ymin = $(ymax-ymin)")
    d = fit(Histogram, (cprime[:,1], cprime[:,2]), (xmin:step:xmax, ymin:step:ymax))
    #println("edges = $(d.edges)\n")
    (d,heatmap!(ax,d.edges[1], d.edges[2], d.weights))
end


"""
    computeExitPupilLoc(geo; epsilon=0.001, format="quiet")

Estimate the z location (in the "stop" surface's local frame) of the
system's exit pupil, by tracing a small-`epsilon`-angle ray from the
stop through the rest of `geo` and finding where its y (and by
assumption x) coordinate crosses zero. `format` is passed to
`printTrcCoords` when not `"quiet"` (any other value prints the trace).

Returns `(status, zpupil)`:
- `status` -- `0` on success, `1` if the reference ray didn't make it
  through `geo` (in which case the second value is `ORIGIN`, a
  `Point3`, not a `Float64` -- callers must check `status` before using
  the second value, since its type differs between the two cases)
- `zpupil::Float64` -- the estimated exit pupil z location (on success)

Only handles the case where the reference ray's final direction has a
nonzero y component (`dirb[2] != 0.`); otherwise returns `zpupil = NaN`.
The commented-out block below flags this as a known simplification
("should calculate the distance in local coordinates rather than
assume the global Z direction").
"""
function computeExitPupilLoc(geo; epsilon = 0.001, format="quiet")
    surfnumStop = tracenumFromName("stop", geo)-1
    locgeo = geo[surfnumStop:end]
    status,trc = traceGeometryRel(Ray(ORIGIN, Vec3(0., sin(epsilon), cos(epsilon))), locgeo)
    #check to make sure all made it through
    if status !=0
        return(1, ORIGIN)
    end
    if format!="quiet"
        # takes "normal" or something else
        printTrcCoords(status, trc, geo, format=format)
    end
    #compute intersection of bore and ytrace, xtrace
    #get last ray, rotate input rays to direction of bore, find when length that y, x are zero

    rayb = trc[end].ray

    dirb = rayb.dir

    baseb = rayb.base
#=
    if dirb[3] != 1.
        println("not yet general calculations $dirb")
        return(2, NaN, NaN, ORIGIN, ORIGIN)
     end
=#

#=
 TBD should calculate the distance in local coordinates rather than assume the global Z direction
 something like

    localdir = geo[end].toLocalDir(dirb)
    localbase = geo[end].toLocalCoord(baseb)

=#
    if dirb[2] != 0.
        lenzero = -baseb[2]/dirb[2]
    else
        lenzero = NaN
    end


    zpupil = baseb[3] + lenzero * dirb[3]


    return (0, zpupil)
end


"""
    plotOPD!(scene, point, max angle, geometry; color=:blue ,surfview=Surface Name, points=33, θmin = min angle )
    plots the x- and y-fan OPD curves into scene
    returns θr, opdx, opdy

    current version assumes telecentric pupil/stop (i.e. reference ray θ=0)

"""
function plotOPD!(scene, r::Point3, θmax::Float64, geo; surfview = "end", color = :blue, points=33, θmin::Float64=NaN, offset = 0., λ=1.0,label="")
    #trcStatMsg=("Normal","Missed","TIR","Clipped")
    opdy = Vector{Float64}(undef, points)
    opdx = Vector{Float64}(undef, points)

    #surfview is assumed to be an image plane

    surfnum = tracenumFromName(surfview, geo)
    localgeo = geo[1:surfnum-1]  #surfnum is defined for traces

    #find the reference local reference intercept coordinates
    #println("base = $r")
    status, refTrace = traceGeometryRel(Ray(r, ZAXIS),localgeo) #ZAXIS is default local direction
    if status != 0
        #println(trcStatMsg[status+1])
        println("Reference ray did not intersect surface: $surfview")
        return
    end

    status,zExitPupilLoc=computeExitPupilLoc(geo, epsilon = 0.0001, format="normal")

    if status != 0
        #println(trcStatMsg[status+1])
        println("Can't compute ExitPupilLoc")
        return
    end

    #put exit pupil in geo

    if abs(zExitPupilLoc)> 1e5 #it's at infinity
        zExitPupilLoc = refTrace[end-1].ray.base[3]  #put it at prior Surface
    end

    baseImage = refTrace[end].ray.base
    dirImage =refTrace[end].ray.dir

    radiusRefSphere = (baseImage[3]-zExitPupilLoc+offset)/dirImage[3]
    curvRefSphere = 1.0/radiusRefSphere
    #should have x & y = 0 if exitpupil unless telecentric system
    baseExitPupil = Point3(0., 0., zExitPupilLoc)

    println("zExitPupilLoc = $zExitPupilLoc  radiusRefSphere = $radiusRefSphere")

    testgeo = [localgeo[1:end-1]
            [refractSphere("reference sphere",  baseExitPupil, dirImage, refIndexDefault, refIndexDefault, curvRefSphere, 25., "testcoat")]
        #    [localgeo[end]]
        ]
    printSurfNames(testgeo)
    status, refTrace = traceGeometryRel(Ray(r, ZAXIS),testgeo) #ZAXIS is default local direction




    θm = isnan(θmin) ? -θmax : θmin

    θr =  LinRange(θm, θmax, points)


    for (i,θ) in enumerate(θr)
        #println("dir = $([0.,sin(θ), cos(θ)])")
        opdy[i] = opdRel(Ray(r, Vec3(0.,sin(θ), cos(θ))), refTrace,testgeo)*1000.0/λ
        opdx[i] = opdRel(Ray(r, Vec3(sin(θ), 0., cos(θ))), refTrace,testgeo)*1000.0/λ
    end
    ax = Axis(scene[1,1]; title=label)
    Makie.lines!(ax, θr, opdx, color=color, linestyle = :dash)
    Makie.lines!(ax, θr, opdy , color=color, label=label)
    θr, opdx, opdy
end

"""
    plotOPD!(scene, h::Float64, egeo::ExtendedGeometry; surfstop="stop", surfview="end", color=:blue, points=33, focusOffset=0., label="")

`ExtendedGeometry` sibling of `plotOPD!(scene, r::Point3, θmax::Float64,
geo; ...)` above: rebuilds `egeo`'s geometry (`updateEGeo!`), picks a
reference ray at normalized object height `h` (`0` to `1`) through the
`egeo.surfaceObject`'s aperture, and plots x/y OPD-vs-position curves
(rather than OPD-vs-angle, as the `geo` method does) into `scene`.
Requires `egeo.geo[1]` to be named `"stop"` or `"pupil"`; prints a
message and returns `nothing` if not (or if the reference/exit-pupil
tracing fails).

Returns `(opdx, opdy)` (unlike its `θmax`-based sibling above, this
method sweeps `x`/`y` position rather than angle, so there's no `θr`
to return alongside them).
"""
function plotOPD!(scene, h::Float64, egeo::ExtendedGeometry; surfstop = "stop", surfview = "end", color = :blue, points=33, focusOffset = 0., label="")
    #trcStatMsg=("Normal","Missed","TIR","Clipped")
    opdy = Vector{Float64}(undef, points)
    opdx = Vector{Float64}(undef, points)

    updateEGeo!(egeo)

    #surfview is assumed to be an image plane
    geo = egeo.geo
    if !(geo[1].surfname == "stop" || geo[1].surfname == "pupil")
        println("first surface not stop or pupil")
        return
    end
    localgeo = geo
    surfnum = tracenumFromName(surfview, localgeo)
    usedgeo = geo[1:surfnum-1]  #surfnum is defined for traces
    #these are just for test
    sizeO = h * sizeOptic(egeo.surfaceObject.aperture) # h is 0 to 1
    sizeP = sizeOptic(usedgeo[1].aperture)

    r = Vec3(0., sizeO, egeo.surfaceObject.base.base[3])
    z = usedgeo[1].base.base[3]

    dirRef = normalize(Vec3(usedgeo[1].base.base .- r))
    #println("r = $r  dirRef = $dirRef")
    #find the reference local reference intercept coordinates
    #println("base = $r")
    status, refTrace = traceGeometryRel(Ray(ORIGIN, dirRef),usedgeo) #ZAXIS is default local direction
    if status != 0
        #println(trcStatMsg[status+1])
        println("Reference ray did not intersect surface: $surfview")
        return
    end

    status,zExitPupilLoc=computeExitPupilLoc(usedgeo, epsilon = 0.0001, format="quiet")

    if status != 0
        #println(trcStatMsg[status+1])
        println("Can't compute ExitPupilLoc")
        return
    end

    #put exit pupil in geo

    if abs(zExitPupilLoc)> 1e5 #it's at infinity
        zExitPupilLoc = refTrace[end-1].ray.base[3]  #put it at prior Surface
    end

    baseImage = refTrace[end].ray.base
    dirImage =refTrace[end].ray.dir

    radiusRefSphere = (baseImage[3]-zExitPupilLoc+focusOffset)/dirImage[3]
    curvRefSphere = 1.0/radiusRefSphere
    #should have x & y = 0 if exitpupil unless telecentric system
    baseExitPupil = Point3(0., 0., zExitPupilLoc)

    #println("zExitPupilLoc = $zExitPupilLoc  radiusRefSphere = $radiusRefSphere")

    finalgeo = [usedgeo[1:end-1]
            [refractSphere("reference sphere",  baseExitPupil, dirImage, refIndexDefault, refIndexDefault, curvRefSphere, 25., "testcoat")]
        #    [localgeo[end]]
        ]
    #printSurfNames(finalgeo)
    status, refTrace = traceGeometryRel(Ray(ORIGIN, dirRef),finalgeo) #ZAXIS is default local direction



    x =  LinRange(-sizeP, sizeP, points)
    z = usedgeo[1].base.base[3]
    wl = egeo.wavelength[1] * 1e-3

    for (i,t) in enumerate(x)
        #print("t = $t, i = $i  ")
        #println("z = $z, r = $r")
        bx = Point3(t, 0., 0.)
        by = Point3(0., t, 0.)
        px = normalize(Vec3(t, 0., z).-r)
        py = normalize(Vec3(0., t, z).-r)
        #println("bx = $bx px = $px  by = $by py = $py  ")
        opdy[i] = opdRel(Ray(by, py), refTrace, finalgeo)/wl
        opdx[i] = opdRel(Ray(bx, px), refTrace, finalgeo)/wl
    end

    Makie.lines!(scene, x, opdx, color=color, linestyle = :dash)
    Makie.lines!(scene, x, opdy , color=color, label=label)
    opdx, opdy
end

"""
    plotOPD3D!(scene, h::Float64, egeo::ExtendedGeometry; surfstop="stop", surfview="end", color=:blue, points=33, focusOffset=0., label="")

Like `plotOPD!(scene, h::Float64, egeo::ExtendedGeometry; ...)` above
(same reference-ray/exit-pupil setup, same `"stop"`/`"pupil"`
requirement), but plots a full 2D OPD *surface* over an `x`/`y` grid
(`Makie.surface!`) instead of two 1D x/y cross-section curves.

Returns `scene`.
"""
function plotOPD3D!(scene, h::Float64, egeo::ExtendedGeometry; surfstop = "stop", surfview = "end", color = :blue, points=33, focusOffset = 0., label="")

    updateEGeo!(egeo)


    geo = egeo.geo
    if !(geo[1].surfname == "stop" || geo[1].surfname == "pupil")
        println("first surface not stop or pupil")
        return
    end
    localgeo = geo
    surfnum = tracenumFromName(surfview, localgeo)
    usedgeo = geo[1:surfnum-1]  #surfnum is defined for traces
    #these are just for test
    sizeO = h * sizeOptic(egeo.surfaceObject.aperture) # h is 0 to 1
    sizeP = sizeOptic(usedgeo[1].aperture)

    r = Vec3(0., sizeO, egeo.surfaceObject.base.base[3])

    dirRef = normalize(Vec3(usedgeo[1].base.base .- r))
    #println("r = $r  dirRef = $dirRef")
    #find the reference local reference intercept coordinates
    #println("base = $r")
    status, refTrace = traceGeometryRel(Ray(ORIGIN, dirRef),usedgeo) #ZAXIS is default local direction
    if status != 0
        #println(trcStatMsg[status+1])
        println("Reference ray did not intersect surface: $surfview")
        return
    end

    status,zExitPupilLoc=computeExitPupilLoc(usedgeo, epsilon = 0.0001, format="quiet")

    if status != 0
        #println(trcStatMsg[status+1])
        println("Can't compute ExitPupilLoc")
        return
    end

    #put exit pupil in geo

    if abs(zExitPupilLoc)> 1e5 #it's at infinity
        zExitPupilLoc = refTrace[end-1].ray.base[3]  #put it at prior Surface
    end

    baseImage = refTrace[end].ray.base
    dirImage =refTrace[end].ray.dir

    radiusRefSphere = (baseImage[3]-zExitPupilLoc+focusOffset)/dirImage[3]
    curvRefSphere = 1.0/radiusRefSphere
    #should have x & y = 0 if exitpupil unless telecentric system
    baseExitPupil = Point3(0., 0., zExitPupilLoc)

    #println("zExitPupilLoc = $zExitPupilLoc  radiusRefSphere = $radiusRefSphere")

    finalgeo = [usedgeo[1:end-1]
            [refractSphere("reference sphere",  baseExitPupil, dirImage, refIndexDefault, refIndexDefault, curvRefSphere, 25., "testcoat")]
        #    [localgeo[end]]
        ]
    #printSurfNames(finalgeo)
    status, refTrace = traceGeometryRel(Ray(ORIGIN, dirRef),finalgeo) #ZAXIS is default local direction



    x =  LinRange(-sizeP, sizeP, points)
    y =  LinRange(-sizeP, sizeP, points)
    z = finalgeo[1].base.base[3]

    wl = egeo.wavelength[1] * 1e-3
    opdfunc(xi, yi) = opdRel(Ray(Point3(xi, yi, 0.), normalize(Vec3(xi, yi, z).-r)), refTrace, finalgeo)/wl
    opd = [opdfunc(xi, yi) for xi in x, yi in y]


    Makie.surface!(scene, x,y,opd, label=label)
    #scale!(scene, 1, 1, 1000)
    scene
end


"""
    plotXSag!(scene, xmax, ycut, profile)

Plot `sag(x, ycut, profile)` for `x` ranging over `(-xmax, xmax)` (160
points) into `scene` -- a cross-section of the surface's sag along the
local x direction, at fixed y = `ycut`. See `plotYSag!` below for the
analogous y-direction cross-section.
"""
function plotXSag!(scene, xmax, ycut, profile)
    sag1(x) = sag(x, Float64(ycut), profile)
    x = range(-xmax, stop = xmax, length = 160)
    z = sag1.(x)
   lines!(scene, x,z, color=:blue)
end

"""
    plotYSag!(scene, xmax, ycut, profile)

Plot `sag(ycut, x, profile)` for `x` ranging over `(-xmax, xmax)` (160
points) into `scene` -- a cross-section of the surface's sag along the
local y direction, at fixed x = `ycut` (despite the parameter's name,
which mirrors `plotXSag!`'s -- there, `ycut` fixes y; here, it fixes
x). See `plotXSag!` above for the analogous x-direction cross-section.
"""
function plotYSag!(scene, xmax, ycut, profile)
    sag1(x) = sag(Float64(ycut), x, profile)
    x = range(-xmax, stop = xmax, length = 160)
    z = sag1.(x)
    lines!(scene, x,z, color=:red)
end

"""
    plotSpotDiagram(fig, spts, center, rmsradius, deltaz; title = "Spot Diagram", showRMS=true)
    returns a figure with the spot diagram and optionally the RMS radius and ΔZ
    It uses the results from spotDiagramHex()
"""

function plotSpotDiagram(fig, spts, center, rmsradius, deltaz; title = "Spot Diagram", showRMS=true)
    ax = Axis(fig[1,1]; title, tellwidth = false)
    scatter!(ax, spts, markersize=2, color=:blue)
    ax.aspect = DataAspect()
    if showRMS
        arc!(ax, center, rmsradius, -π, π, color=:red)
        Label(fig[2,1], @sprintf("RMS Radius: %8.3f  ΔZ: %8.3g", rmsradius, deltaz), tellwidth = false)
    end
    fig
end

"""
    plotSpotDiagram(fig, spts; title = "Spot Diagram")

Like `plotSpotDiagram(fig, spts, center, rmsradius, deltaz; ...)`
above, but without the RMS-radius circle/label overlay -- just the
scatter plot of `spts`. See that method's docstring above for the
general contract shared by both.
"""
function plotSpotDiagram(fig, spts; title = "Spot Diagram")
    ax = Axis(fig[1,1]; title, tellwidth = false)
    scatter!(ax, spts, markersize=2, color=:blue)
    ax.aspect = DataAspect()
    fig
end
