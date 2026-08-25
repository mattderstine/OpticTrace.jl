@testset "printing.jl" begin

    # printing.jl's functions all write straight to the global stdout
    # stream (println/@printf), rather than taking an io argument, so
    # capturing what they print means redirecting the real stdout stream
    # to a temp file -- an IOBuffer isn't supported as a redirect_stdout
    # target on this Julia version. Returns (output::String, result) so
    # callers can check both what got printed and what the function
    # returned.
    function captureOutput(f)
        result = Ref{Any}(nothing)
        output = mktemp() do path, io
            redirect_stdout(io) do
                result[] = f()
            end
            close(io)
            read(path, String)
        end
        return output, result[]
    end

    # geo[1]=object, geo[2]="TG_1", geo[3]="TG_2", geo[4]=image
    geo = [[referencePlane("object", ORIGIN, ZAXIS, 1.0, 10.0, "none")]
        lensSinglet(Point3(0.0, 0.0, 5.0), ZAXIS, 0.02, -0.015, 3.0, 0.5, riFunc, 5.0; order = "forward", lensname = "TG");
        [referencePlane("image", Point3(0.0, 0.0, 50.0), ZAXIS, 1.0, 10.0, "none")]
    ]
    ray = Ray(Point3(0.5, 0.3, -5.0), Vec3(0.0, 0.0, 1.0))
    status, trc = traceGeometry(ray, geo)
    # geo[1] (object) sits at ORIGIN/ZAXIS with an identity local<->global
    # transform (see test/trace_geometry.jl), so this local ray traces
    # identically to `ray` above via the Rel variants.
    localRay = Ray(geo[1].toLocalCoord(ray.base), geo[1].toLocalDir(ray.dir))

    @testset "numsurfFromName / tracenumFromName / surfaceFromName" begin
        @test numsurfFromName("TG_1", geo) == 2
        @test numsurfFromName("image", geo) == 4
        @test numsurfFromName("end", geo) == length(geo)

        out, val = captureOutput(() -> numsurfFromName("nope", geo))
        @test occursin("not found", out)
        @test val == length(geo) # falls back to "end" behavior

        @test tracenumFromName("TG_1", geo) == 3
        @test tracenumFromName("end", geo) == length(geo) + 1

        @test surfaceFromName("TG_1", geo) === geo[2]
        @test surfaceFromName("image", geo) === geo[end]
    end

    @testset "printTrcCoords / trcAndPrintRay / trcAndPrintRayRel" begin
        out, _ = captureOutput(() -> printTrcCoords(status, trc, geo))
        @test occursin("Start", out)
        @test occursin("TG_1", out)
        @test occursin("TG_2", out)
        @test occursin("image", out)

        outRaw, _ = captureOutput(() -> printTrcCoords(status, trc, geo; format = "raw"))
        @test occursin("Start", outRaw)

        _, trcRet = captureOutput(() -> trcAndPrintRay(ray, geo))
        @test trcRet[end].ray.base ≈ trc[end].ray.base
        @test trcRet[end].ray.dir ≈ trc[end].ray.dir

        _, trcRelRet = captureOutput(() -> trcAndPrintRayRel(localRay, geo))
        @test trcRelRet[end].ray.base ≈ trc[end].ray.base
        @test trcRelRet[end].ray.dir ≈ trc[end].ray.dir
    end

    @testset "printTrcLen / trcAndPrintLengths / trcAndPrintLengthsRel" begin
        expectedDelta = sum(t.delta for t in trc)
        expectedOPL = sum(t.delta * t.nIn for t in trc)
        expectedRD = sum(t.delta / t.nIn for t in trc)
        expected = [expectedDelta, expectedOPL, expectedRD]

        out, result = captureOutput(() -> printTrcLen(status, trc, geo))
        @test result ≈ expected
        @test occursin("Total", out)

        _, result2 = captureOutput(() -> trcAndPrintLengths(ray, geo))
        @test result2 ≈ expected

        _, result3 = captureOutput(() -> trcAndPrintLengthsRel(localRay, geo))
        @test result3 ≈ expected
    end

    @testset "printSurfNames" begin
        out, _ = captureOutput(() -> printSurfNames(geo))
        @test occursin("object", out)
        @test occursin("TG_1", out)
        @test occursin("TG_2", out)
        @test occursin("image", out)

        outFull, _ = captureOutput(() -> printSurfNames(geo; fulldir = true))
        @test length(outFull) > length(out) # extra ydir column printed per line
    end

    @testset "printMissed" begin
        # column i+1 of m corresponds to geo[i] (printMissed's loop starts
        # i=2 before the first, geo[1]-referring, iteration)
        m = zeros(Int32, 3, length(geo) + 1)
        m[2, 3] = 5 # flags geo[2] == "TG_1"

        out, _ = captureOutput(() -> printMissed(m, geo))
        @test occursin("TG_1", out)
        @test !occursin("object", out)
        @test !occursin("TG_2", out)
        @test !occursin("image", out)
    end

    @testset "printSurface / printGeo" begin
        out1, _ = captureOutput(() -> printSurface(geo[2]))
        @test occursin("TG_1", out1)

        out1n, _ = captureOutput(() -> printSurface(3, geo[2]))
        @test occursin("TG_1", out1n)
        @test startswith(out1n, "3")

        modelSurf = roundAperture("stop", ORIGIN, ZAXIS, 1.0, 1.0, 5.0)
        outModel, _ = captureOutput(() -> printSurface(1, modelSurf))
        @test occursin("stop", outModel)

        outGeo, _ = captureOutput(() -> printGeo(geo))
        @test occursin("object", outGeo)
        @test occursin("TG_1", outGeo)
        @test occursin("TG_2", outGeo)
        @test occursin("image", outGeo)
    end

    @testset "printTrcStatus" begin
        out0default, _ = captureOutput(() -> printTrcStatus(0))
        @test out0default == ""

        out0flag, _ = captureOutput(() -> printTrcStatus(0; flagNormal = true))
        @test occursin("Normal", out0flag)

        out1, _ = captureOutput(() -> printTrcStatus(1))
        @test occursin("Missed", out1)

        outSuppressed, _ = captureOutput(() -> printTrcStatus(1; clipmsg = false))
        @test outSuppressed == ""
    end

    @testset "opdRel" begin
        # a ray compared against its own trace has zero OPD
        @test opdRel(localRay, trc, geo) ≈ 0.0 atol = 1e-9

        rayOffset = Ray(localRay.base + Vec3(0.1, 0.0, 0.0), localRay.dir)
        opd = opdRel(rayOffset, trc, geo)
        @test opd isa Float64
        @test !isnan(opd)

        # status > 0 (TIR here) -> NaN, per opdRel's own early-return guard
        tirSurf = refractConic("tir", ORIGIN, ZAXIS, 1.5, 1.0, 0.0, 0.0, 5.0, "none")
        tirGeo = [tirSurf]
        θ = deg2rad(60.0) # beyond the critical angle for 1.5->1.0
        tirRay = Ray(Point3(0.0, 0.0, -1.0), Vec3(sin(θ), 0.0, cos(θ)))
        statusTir, trcTir = traceGeometryRel(tirRay, tirGeo)
        @test statusTir == 2
        @test isnan(opdRel(tirRay, trcTir, tirGeo))
    end

end
