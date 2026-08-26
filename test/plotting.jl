@testset "plotting.jl" begin

    # GLMakie needs a display -- run headless (no window flashes, no
    # interactive event loop) rather than requiring xvfb-run locally.
    # This must run before any Figure/Scene is created; the setting
    # persists as a merged default even through plotting.jl's own
    # internal `GLMakie.activate!(title=..., inline=false)` calls (e.g.
    # inside multipleFigures), since GLMakie merges into one config dict
    # rather than resetting it.
    GLMakie.activate!(visible = false)

    riFunc(λ) = 1.5

    # geo[1]=object, geo[2]="PL_1", geo[3]="PL_2", geo[4]=image
    geo = AbstractSurface[[referencePlane("object", ORIGIN, ZAXIS, 1.0, 10.0, "none")]
        lensSinglet(Point3(0.0, 0.0, 5.0), ZAXIS, 0.02, -0.015, 3.0, 0.5, riFunc, 5.0; order = "forward", lensname = "PL");
        [referencePlane("image", Point3(0.0, 0.0, 50.0), ZAXIS, 1.0, 10.0, "none")]
    ]
    ray = Ray(Point3(0.5, 0.3, -5.0), Vec3(0.0, 0.0, 1.0))

    @testset "plotSurface3D! / plotGeometry3D! / plotGeometry3D" begin
        fig = Figure()
        f, ax = plotGeometry3D(fig, geo)
        @test f === fig
        @test ax isa LScene

        fig2, ax2 = plotGeometry3D(geo)
        @test fig2 isa Figure
        @test ax2 isa LScene

        # ModelSurface siblings: round aperture (with & without a central
        # obscuration) and rect aperture
        modelRound = roundAperture("stop_round", ORIGIN, ZAXIS, 1.0, 0.0, 5.0)
        modelRoundObsc = roundAperture("stop_round_obsc", ORIGIN, ZAXIS, 1.0, 1.0, 5.0)
        modelRect = rectAperture("stop_rect", ORIGIN, ZAXIS, YAXIS, 1.0, 0.0, 0.0, 3.0, 4.0)

        for ms in (modelRound, modelRoundObsc, modelRect)
            figm, axm = plotGeometry3D([ms])
            @test axm isa LScene
        end

        scene = plotGeometry3D!(ax, geo)
        @test scene === ax
    end

    @testset "plotModelSurf!" begin
        scene = LScene(Figure()[1, 1])

        @testset "RoundAperture" begin
            # no obscuration, finite semiDiameter -> draws the Washer ring,
            # so the last-executed `if` block's value (its return) is non-nothing
            ms1 = roundAperture("r1", ORIGIN, ZAXIS, 1.0, 0.0, 3.0)
            @test !isnothing(OpticTrace.plotModelSurf!(scene, ms1.aperture, ms1))

            # obscuration only (semiDiameter = ∞) -> only the Disk branch
            # runs; the un-taken `if a.semiDiameter != ∞` block is the
            # function's last statement, so it implicitly returns nothing
            # even though the Disk mesh was drawn -- smoke test only
            ms2 = roundAperture("r2", ORIGIN, ZAXIS, 1.0, 1.0, ∞)
            OpticTrace.plotModelSurf!(scene, ms2.aperture, ms2)
            @test true
        end

        @testset "RectAperture" begin
            # no obscuration, finite clear aperture -> draws the frame,
            # so the last-executed `if` block's value is non-nothing
            ms3 = rectAperture("q1", ORIGIN, ZAXIS, YAXIS, 1.0, 0.0, 0.0, 3.0, 4.0)
            @test !isnothing(OpticTrace.plotModelSurf!(scene, ms3.aperture, ms3))

            # obscuration present, infinite clear aperture -> only the
            # obscuration-rectangle branch runs; the un-taken frame `if`
            # is the last statement, so this implicitly returns nothing
            # even though the obscuration rectangle was drawn -- smoke
            # test only
            ms4 = rectAperture("q2", ORIGIN, ZAXIS, YAXIS, 1.0, 1.0, 1.0, ∞, ∞)
            OpticTrace.plotModelSurf!(scene, ms4.aperture, ms4)
            @test true
        end
    end

    @testset "saveFigure / printFigure" begin
        fig = Figure()
        ax = Axis(fig[1, 1])
        lines!(ax, 1:10, rand(10))

        mktempdir() do dir
            saveFigure("stub", fig; directory = dir)
            @test isfile(joinpath(dir, "stub.png"))

            saveFigure("stub", fig; startnum = 2, directory = dir)
            @test isfile(joinpath(dir, "stub2.png"))

            out = mktemp() do path, io
                redirect_stdout(io) do
                    printFigure("stub3", fig; directory = dir)
                end
                close(io)
                read(path, String)
            end
            @test occursin("depricated", out)
            @test isfile(joinpath(dir, "stub3.png"))
        end
    end

    @testset "multipleFigures" begin
        figs = multipleFigures(2)
        @test figs isa Vector{Figure}
        @test length(figs) == 2
    end

    @testset "trace-and-plot family" begin
        fig = Figure()
        scene = LScene(fig[1, 1])

        out, trc = begin
            result = Ref{Any}(nothing)
            captured = mktemp() do path, io
                redirect_stdout(io) do
                    result[] = trcAndPrintPlot!(ray, geo)
                end
                close(io)
                read(path, String)
            end
            captured, result[]
        end
        @test occursin("Start", out)
        @test length(trc) == length(geo) + 1

        _, trc2 = mktemp() do path, io
            r = Ref{Any}(nothing)
            redirect_stdout(io) do
                r[] = trcAndPrintPlotRay!(scene, ray, geo)
            end
            close(io), r[]
        end
        @test trc2[end].ray.base ≈ trc[end].ray.base

        status3, trc3 = trcAndPlotRay!(scene, ray, geo)
        @test status3 == 0
        @test trc3[end].ray.base ≈ trc[end].ray.base

        localRay = Ray(geo[1].toLocalCoord(ray.base), geo[1].toLocalDir(ray.dir))

        status4, trc4 = trcAndPlotRayRel!(scene, localRay, geo)
        @test status4 == 0
        @test trc4[end].ray.base ≈ trc[end].ray.base

        status5, trc5 = trcAndPlotRayRel!(localRay, geo) # active-scene variant
        @test status5 == 0
        @test trc5[end].ray.base ≈ trc[end].ray.base

        sceneReturned = plotTrace!(scene, trc)
        @test sceneReturned === scene
    end

    @testset "perimeterRays / plotPerimeterRays" begin
        r = SVector(0.0, 0.0, -5.0)

        pr = perimeterRays(r, 1.0, 0.01, 8, geo; surfview = "end")
        @test length(pr) == 8
        @test all(ray -> !any(isnan, ray.base) && !any(isnan, ray.dir), pr)
        # surfview="end" -> the image plane, geo's last surface (z=50)
        @test all(ray -> ray.base[3] ≈ 50.0, pr)
        # None of these 8 rays miss, so the truncated-length return path
        # (FIXED.md #2) isn't exercised here.

        plt = plotPerimeterRays(r, 1.0, 0.01, 8, geo; surfview = "end")
        @test !isnothing(plt)

        fig = Figure()
        scene = Axis3(fig[1, 1])
        @test !isnothing(plotPerimeterRays!(scene, r, 1.0, 0.01, 8, geo; surfview = "end"))
        @test !isnothing(plotPerimeterRays!(r, 1.0, 0.01, 8, geo; surfview = "end"))
    end

    @testset "rayHeatmap / rayHeatmap!" begin
        pts = [SVector(0.1i, 0.1i) for i in -10:10]

        d, fig = rayHeatmap(pts; mcbins = 10, center = (0.0, 0.0), width = 10.0)
        @test d isa StatsBase.Histogram
        @test sum(d.weights) == length(pts)

        fig2 = Figure()
        ax2 = Axis(fig2[1, 1])
        d2, plt2 = rayHeatmap!(ax2, pts; mcbins = 10, center = (0.0, 0.0), width = 10.0)
        @test d2.weights == d.weights
        @test !isnothing(plt2)
    end

    @testset "computeExitPupilLoc / getrefbase / plotRayFan!" begin
        stop = roundAperture("stop", ORIGIN, ZAXIS, 1.0, 0.0, 2.0)
        lensgeo = lensSinglet(ORIGIN, ZAXIS, 0.05, -0.04, 3.0, 0.5, riFunc, 5.0; order = "forward", lensname = "RF")
        rfGeo = AbstractSurface[stop; lensgeo]

        status, zpupil = computeExitPupilLoc(rfGeo)
        @test status == 0
        @test zpupil isa Float64
        @test !isnan(zpupil)

        # no "stop" surface at all -> tracenumFromName falls back to "end"
        # (prints a "not found" message, per its own documented behavior)
        noStopGeo = AbstractSurface[lensgeo...]
        out, (statusNS, zNS) = mktemp() do path, io
            r = Ref{Any}(nothing)
            redirect_stdout(io) do
                r[] = computeExitPupilLoc(noStopGeo)
            end
            close(io)
            read(path, String), r[]
        end
        @test occursin("not found", out)

        refbase = OpticTrace.getrefbase(ORIGIN, rfGeo, "end")
        @test refbase ≈ Point3(0.0, 0.0, 0.0)

        fig = Figure()
        result = redirect_stdout(devnull) do
            plotRayFan!(fig[1, 1], ORIGIN, 0.05, rfGeo; surfview = "end", points = 5)
        end
        @test result ≈ Point3(0.0, 0.0, 0.0)
    end

    @testset "plotOPD! (geo-based)" begin
        stop = roundAperture("stop", ORIGIN, ZAXIS, 1.0, 0.0, 2.0)
        lensgeo = lensSinglet(ORIGIN, ZAXIS, 0.05, -0.04, 3.0, 0.5, riFunc, 5.0; order = "forward", lensname = "OD")
        opdGeo = AbstractSurface[stop; lensgeo]

        fig = Figure()
        θr, opdx, opdy = redirect_stdout(devnull) do
            plotOPD!(fig, ORIGIN, 0.05, opdGeo; surfview = "end", points = 5)
        end
        @test length(θr) == 5
        @test length(opdx) == 5
        @test length(opdy) == 5
        @test all(!isnan, opdx)
        @test all(!isnan, opdy)
        @test opdx[3] ≈ 0.0 atol = 1e-6 # on-axis (θ=0) ray has ~zero OPD vs itself
        @test opdy[3] ≈ 0.0 atol = 1e-6
        # symmetric fan angles about a telecentric on-axis reference -> symmetric OPD
        @test opdx[1] ≈ opdx[end] atol = 1e-6
        @test opdy[1] ≈ opdy[end] atol = 1e-6
    end

    @testset "plotOPD!(egeo) / plotOPD3D!" begin
        object = referencePlane("object", Point3(0.0, 0.0, -100.0), ZAXIS, 1.0, 5.0, "none")
        stop = roundAperture("stop", ORIGIN, ZAXIS, 1.0, 0.0, 2.0)
        lensgeo = lensSinglet(ORIGIN, ZAXIS, 0.05, -0.04, 3.0, 0.5, riFunc, 5.0; order = "forward", lensname = "OE")
        egeoGeo = AbstractSurface[stop; lensgeo]
        # updateEGeo! now actually assigns funcGeo's result into egeo.geo, so
        # funcGeo has to return a real Array{AbstractSurface} (not a no-op).
        buildGeo(p, wl) = egeoGeo
        egeo = ExtendedGeometry(egeoGeo, buildGeo, OpticTrace.defaultSetupGeo, object, [0.5], Dict{Symbol,Any}())

        fig = Figure()
        ax = Axis(fig[1, 1])
        opdx, opdy = redirect_stdout(devnull) do
            plotOPD!(ax, 0.5, egeo; points = 5)
        end
        @test length(opdx) == 5
        @test length(opdy) == 5
        @test all(!isnan, opdx)
        @test all(!isnan, opdy)

        fig2 = Figure()
        ax2 = LScene(fig2[1, 1])
        result = redirect_stdout(devnull) do
            plotOPD3D!(ax2, 0.5, egeo; points = 5)
        end
        @test result === ax2
    end

    @testset "plotXSag! / plotYSag!" begin
        surf = refractSphere("sag1", ORIGIN, ZAXIS, 1.0, 1.5, 0.02, 5.0, "none")
        fig = Figure()
        ax = Axis(fig[1, 1])

        pltX = plotXSag!(ax, 4.0, 0.0, surf.profile)
        @test !isnothing(pltX)
        expectedZ = [sag(x, 0.0, surf.profile) for x in range(-4.0, stop = 4.0, length = 160)]
        @test [p[2] for p in pltX[1][]] ≈ expectedZ

        pltY = plotYSag!(ax, 4.0, 0.0, surf.profile)
        @test !isnothing(pltY)
        expectedZY = [sag(0.0, x, surf.profile) for x in range(-4.0, stop = 4.0, length = 160)]
        @test [p[2] for p in pltY[1][]] ≈ expectedZY
    end

    @testset "plotSpotDiagram" begin
        spts = [Point2(0.1i, -0.1i) for i in -5:5]

        fig = Figure()
        result = plotSpotDiagram(fig, spts, Point2(0.0, 0.0), 1.0, 0.01)
        @test result === fig

        fig2 = Figure()
        result2 = plotSpotDiagram(fig, spts, Point2(0.0, 0.0), 1.0, 0.01; showRMS = false)
        @test result2 === fig

        fig3 = Figure()
        result3 = plotSpotDiagram(fig3, spts)
        @test result3 === fig3
    end

end
