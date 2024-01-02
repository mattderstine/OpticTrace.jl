export saveFigure, printFigure, multipleFigures, sizeOptic, sizeOpticSurface,plotSurface3D!,plotGeometry3D
export plotGeometry3D!
export trcAndPrintPlot!, trcAndPrintPlotRay!, trcAndPlotRay!,trcAndPlotRayRel!
export plotRayFan!,perimeterRays,plotPerimeterRays
export plotPerimeterRays!, rayHeatmap, rayHeatmap!, computeExitPupilLoc
export plotOPD!, plotOPD3D!
export plotXSag!, plotYSag!
export plotTrace!




"""
    printFigure(fileNameStub, fig; startnum = 0, directory = "")
    saves Makie figure to a file
    returns the fig object

    Eventually will check for the exisitence of the file and add increasing numbers to avoid overwrite
    Currently just writes to fileNameStub*".png" in the current directory, overwriting previous files

    fileNameStub    string with base name for files to be stored
    fig             Makie figure object
    startnum        number to append to the stub. If 0, nothing is appended unless files exist
    directory       location of where to store the figure image

"""
function saveFigure(fileNameStub, fig; startnum = 0, directory = "")
    filename = fileNameStub * ".png"
    save(filename, fig)
    display(fig)
end

#try to depricate this name
printFigure(fileNameStub, fig; startnum = 0, directory = "") = saveFigure(fileNameStub, fig; startnum=startnum, directory = directory)


"""
    multipleFigures(numFigs = 5; fileNameStub = "FigurePrint", activeButtonColor = RGBf0(0.8, 0.94, 0.8), activeFig = 1,  size = (1200,900))
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





sizeOptic(aper::SizeLens) = aper.semiDiameter

sizeOptic(aper::RoundAperture) = aper.semiDiameter

sizeOptic(aper::RectAperture) = max(aper.wclear, aper.lclear)

sizeOpticSurface(surf::AbstractSurface) = sizeOptic(surf.aperture)


#plotSurface3D!(scene, s::OptSurface; color=:aquamarine2)=Makie.mesh!(scene, s, color=color, fxaa = true, transparency = true)#someday add options to display mesh
plotSurface3D!(scene, s::OptSurface; transparency = false)=Makie.mesh!(scene, s, color=s.color, fxaa = true, transparency = transparency)#someday add options to display mesh

function plotGeometry3D!(scene, geometry)
    for geo in geometry
        plotSurface3D!(scene, geo)
    end
    scene
end


"""
plotGeometry3D(geometry; size = (1200,700))
geometry is an array of OptSurfaces
resoltion is the initial window size
returns a Makie Figure and LScene located at [1,1] 

"""
function plotGeometry3D(geometry; size = (1200,700))
    fig = Figure(;size)
    ax = fig[1,1]=LScene(fig, show_axis=false, scenekw = (camera = cam3d_cad!,))
    linesegments!(ax,[Point3f0(0, 0,0) => Point3f0(0,0,1)],color=:green, linewidth=2)
    linesegments!(ax,[Point3f0(0, 0,0) => Point3f0(1,0,0)],color=:blue, linewidth=2)
    linesegments!(ax,[Point3f0(0, 0,0) => Point3f0(0,1,0)],color=:red, linewidth=2)
    plotGeometry3D!(ax, geometry)
    fig,ax
end

function plotSurface3D!(scene, s::ModelSurface)
    plotModelSurf!(scene, s.aperture, s) #dispatch on s.aperture
end

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
trcAndPrintPlot!(scene, ray::Ray, geo; color=:blue)
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
    return the trace
    set clipmsg=true to print a message if the ray does not reach the end of geo
        
"""
function trcAndPlotRay!(scene, ray::Ray, geo; color=:blue, clipmsg=false)
    status, trc = traceGeometry(ray, geo)
    printTrcStatus(status, flagNormal = false, clipmsg=clipmsg)
    Makie.lines!(scene, [a.ray.base for a in trc], color=color)
    status, trc
end

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
    plotRayFan!(scene, point, max angle, geometry; color=:blue ,surfview=Surface Name, points=33, θmin = min angle )
    returns Makie scene of plot of rayfan

    current version assumes telecentric pupil/stop (i.e. reference ray θ=0)

"""
function plotRayFan!(scene, r::SVector{3}, θmax::Float64, geo; surfview = "end", color = :blue, points=33, θmin::Float64=NaN, printTrace=false, rayscene=nothing )
    raysy = Vector{SVector{3, Float64}}(undef, points)
    raysx = Vector{SVector{3, Float64}}(undef, points)
    surfnum = surfnumFromName(surfview, geo)
    s = geo[surfnum-1]  #surfnum is defined for traces

    #find the reference local reference intercept coordinates
    #println("base = $r")
    status, trc = traceGeometryRel(Ray(r, ZAXIS), geo) #ZAXIS is default local direction
    if status != 0 && surfnum >0 && length(trc)<surfnum
        #println(trcStatMsg[status+1])
        println("Reference ray did not intersect surface: $surfview")
        return
    else
        refbase = s.toLocalCoord(surfnum <0 ? trc[end].ray.base : trc[surfnum].ray.base)
    end
    θm = isnan(θmin) ? -θmax : θmin

    θr =  LinRange(θm, θmax, points)

    i = 1
    for θ in θr
        #println("dir = $([0.,sin(θ), cos(θ)])")
        status, trc = traceGeometryRel(Ray(r, [0.,sin(θ), cos(θ)]), geo)
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

        status, trc = traceGeometryRel(Ray(r, [sin(θ), 0., cos(θ)]), geo)
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
    #find the largest of x & y at the edge and use it for the sign
    yv = abs(raysy[end][1]) > abs(raysy[end][2]) ? 1 : 2
    xv = abs(raysx[end][1]) > abs(raysx[end][2]) ? 1 : 2
    #println("yv = $yv   xv = $xv")
    y = [sign(t[yv])*norm(t) for t in raysy] #use norm in case the output intersections are not on an axis
    x = [sign(t[xv])*norm(t) for t in raysx]


    Makie.lines!(scene, θr, x, color=color, linestyle = :dash)
    Makie.lines!(scene, θr, y , color=color)
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

function perimeterRays(r::SVector{3, Float64}, radius::Float64, θ::Float64, points::Int64, geo;surfview = "end")
    trcStatMsg=("Normal","Missed","TIR","Clipped")
    rays = Vector{Ray}(undef, points)
    surfnum = surfnumFromName(surfview, geo)

    i = 1
    for ϕ in LinRange(0., 2pi, points)
        status, trc = traceGeometryRel(Ray(r + radius .* [cos(ϕ), sin(ϕ), 0.] , [cos(ϕ)*sin(θ),sin(ϕ)*sin(θ), cos(θ)]), geo)
        if status != 0
            println(trcStatMsg[status+1])
            rays[i] = Ray([NaN, NaN, NaN], [NaN, NaN, NaN])

            return
        end
        rays[i] = trc[surfnum].ray
        i += 1
    end
    rays
end

function plotPerimeterRays(r::SVector{3, Float64}, radius::Float64, θ::Float64, pnts::Int64, geo; color=:blue, surfview = "end")
    pr = perimeterRays(r, radius, θ, pnts, geo, surfview = surfview)
    Makie.arrows([Makie.Point3f0(a.base) for a in pr], [Makie.Point3f0(a.dir) for a in pr], linecolor=color, arrowcolor=color, arrowsize=0.1)
end

function plotPerimeterRays!(r::SVector{3, Float64}, radius::Float64, θ::Float64, pnts::Int64, geo; color=:blue, surfview = "end")
    pr = perimeterRays(r, radius, θ, pnts, geo, surfview = surfview)
    Makie.arrows!([Makie.Point3f0(a.base) for a in pr], [Makie.Point3f0(a.dir) for a in pr], linecolor=color, arrowcolor=color, arrowsize=0.1)
end

function plotPerimeterRays!(scene, r::SVector{3, Float64}, radius::Float64, θ::Float64, pnts::Int64, geo; color=:blue, surfview = "end")
    pr = perimeterRays(r, radius, θ, pnts, geo, surfview = surfview)
    Makie.arrows!(scene, [Makie.Point3f0(a.base) for a in pr], [Makie.Point3f0(a.dir) for a in pr], linecolor=color, arrowcolor=color, arrowsize=0.1)
end


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


function computeExitPupilLoc(geo; epsilon = 0.001, format="quiet")
    surfnumStop = surfnumFromName("stop", geo)-1
    locgeo = geo[surfnumStop:end]
    status,trc = traceGeometryRel(Ray(ORIGIN, SVector(0., sin(epsilon), cos(epsilon))), locgeo)
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
    returns Makie scene of plot of OPD

    current version assumes telecentric pupil/stop (i.e. reference ray θ=0)

"""
function plotOPD!(scene, r::SVector, θmax::Float64, geo; surfview = "end", color = :blue, points=33, θmin::Float64=NaN, offset = 0., λ=1.0,label="")
    #trcStatMsg=("Normal","Missed","TIR","Clipped")
    opdy = Vector{Float64}(undef, points)
    opdx = Vector{Float64}(undef, points)

    #surfview is assumed to be an image plane

    surfnum = surfnumFromName(surfview, geo)
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
    baseExitPupil = SVector{3}(0., 0., zExitPupilLoc)

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
        opdy[i] = opdRel(Ray(r, [0.,sin(θ), cos(θ)]), refTrace,testgeo)*1000.0/λ 
        opdx[i] = opdRel(Ray(r, [sin(θ), 0., cos(θ)]), refTrace,testgeo)*1000.0/λ
    end

    Makie.lines!(scene, θr, opdx, color=color, linestyle = :dash)
    Makie.lines!(scene, θr, opdy , color=color, label=label)
    θr, opdx, opdy
end

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
    surfnum = surfnumFromName(surfview, localgeo)
    usedgeo = geo[1:surfnum-1]  #surfnum is defined for traces
    #these are just for test
    sizeO = h * sizeOptic(egeo.surfaceObject.aperture) # h is 0 to 1
    sizeP = sizeOptic(usedgeo[1].aperture)

    r = SVector{3}(0., sizeO, egeo.surfaceObject.base.base[3])
    z = usedgeo[1].base.base[3]

    dirRef = normalize!(usedgeo[1].base.base .- r)
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
    baseExitPupil = SVector{3}(0., 0., zExitPupilLoc)

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
        bx = SVector(t, 0., 0.)
        by = SVector(0., t, 0.)
        px = normalize(SVector(t, 0., z).-r)
        py = normalize(SVector(0., t, z).-r)
        #println("bx = $bx px = $px  by = $by py = $py  ")
        opdy[i] = opdRel(Ray(by, py), refTrace, finalgeo)/wl
        opdx[i] = opdRel(Ray(bx, px), refTrace, finalgeo)/wl
    end

    Makie.lines!(scene, x, opdx, color=color, linestyle = :dash)
    Makie.lines!(scene, x, opdy , color=color, label=label)
    θr, opdx, opdy
end

function plotOPD3D!(scene, h::Float64, egeo::ExtendedGeometry; surfstop = "stop", surfview = "end", color = :blue, points=33, focusOffset = 0., label="")

    updateEGeo!(egeo)


    geo = egeo.geo
    if !(geo[1].surfname == "stop" || geo[1].surfname == "pupil")
        println("first surface not stop or pupil")
        return
    end
    localgeo = geo
    surfnum = surfnumFromName(surfview, localgeo)
    usedgeo = geo[1:surfnum-1]  #surfnum is defined for traces
    #these are just for test
    sizeO = h * sizeOptic(egeo.surfaceObject.aperture) # h is 0 to 1
    sizeP = sizeOptic(usedgeo[1].aperture)

    r = SVector{3}(0., sizeO, egeo.surfaceObject.base.base[3])

    dirRef = normalize!(usedgeo[1].base.base .- r)
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
    baseExitPupil = SVector{3}(0., 0., zExitPupilLoc)

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
    opdfunc(xi, yi) = opdRel(Ray(SVector{3}(xi, yi, 0.), normalize(SVector{3}(xi, yi, z).-r)), refTrace, finalgeo)/wl
    opd = [opdfunc(xi, yi) for xi in x, yi in y]

    
    Makie.surface!(scene, x,y,opd, label=label)
    #scale!(scene, 1, 1, 1000)
    scene
end


function plotXSag!(scene, xmax, ycut, profile)
    sag1(x) = sag(x, Float64(ycut), profile)
    x = range(-xmax, stop = xmax, length = 160)
    z = sag1.(x)
   lines!(scene, x,z, color=:blue)
end

function plotYSag!(scene, xmax, ycut, profile)
    sag1(x) = sag(Float64(ycut), x, profile)
    x = range(-xmax, stop = xmax, length = 160)
    z = sag1.(x)
    lines!(scene, x,z, color=:red)
end

