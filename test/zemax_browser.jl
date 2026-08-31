@testset "zemax_browser.jl" begin

    fixturesDir = joinpath(@__DIR__, "fixtures")
    zmxPath = joinpath(fixturesDir, "test_singlet.zmx")
    zarPath = joinpath(fixturesDir, "test_archive.zar")
    zmfPath = joinpath(fixturesDir, "test_catalog.zmf")

    @testset "zemaxFileKind" begin
        @test OpticTrace.zemaxFileKind("foo.zmx") == :zmx
        @test OpticTrace.zemaxFileKind("foo.ZMX") == :zmx
        @test OpticTrace.zemaxFileKind("foo.zar") == :zar
        @test OpticTrace.zemaxFileKind("foo.zmf") == :zmf
        @test OpticTrace.zemaxFileKind("foo.txt") === nothing
        @test OpticTrace.zemaxFileKind("foo") === nothing
    end

    @testset "zemaxFileSummary" begin
        summary = OpticTrace.zemaxFileSummary(zmxPath)
        @test summary isa OpticTrace.ZemaxFileSummary
        @test summary.name == "Test Singlet"
        @test summary.units == "MM"
        @test length(summary.wavelengths) == 24
        @test length(summary.zsurfs) == 4
    end

    @testset "listZemaxArchiveEntities" begin
        @testset ".zar" begin
            entities = OpticTrace.listZemaxArchiveEntities(zarPath)
            @test length(entities) == 2
            byName = Dict(e.name => e for e in entities)
            @test byName["SYNTH.ZMX"].previewable
            @test !byName["SYNTH.AGF"].previewable
            @test byName["SYNTH.ZMX"].byteSize > 0
            @test byName["SYNTH.AGF"].byteSize > 0
        end

        @testset ".zmf" begin
            entities = OpticTrace.listZemaxArchiveEntities(zmfPath)
            @test length(entities) == 2
            @test all(e -> e.previewable, entities)
            @test Set(e.name for e in entities) == Set(["LENS1", "LENS2"])
        end

        @test_throws ErrorException OpticTrace.listZemaxArchiveEntities(zmxPath)
    end

    @testset "zemaxEntityBytes" begin
        bytes = OpticTrace.zemaxEntityBytes(zarPath, "SYNTH.AGF")
        @test String(bytes) == "synthetic glass catalog stub, not real glass data\n"

        zmxBytes = OpticTrace.zemaxEntityBytes(zmfPath, "LENS1")
        @test !isempty(zmxBytes)
        @test startswith(String(zmxBytes), "VERS")

        @test_throws ErrorException OpticTrace.zemaxEntityBytes(zarPath, "NOT_A_REAL_ENTRY")
        @test_throws ErrorException OpticTrace.zemaxEntityBytes(zmxPath, "SYNTH.ZMX")
    end

    @testset "looksLikeText" begin
        @test OpticTrace.looksLikeText(UInt8[])
        @test OpticTrace.looksLikeText(Vector{UInt8}("plain ascii text\nsecond line"))
        @test !OpticTrace.looksLikeText(UInt8[0x68, 0x69, 0x00, 0x6a])
        @test !OpticTrace.looksLikeText(UInt8[0xff, 0xfe, 0xfd])
    end

    @testset "zemaxEntitySummary" begin
        @testset ".zar entity" begin
            summary = OpticTrace.zemaxEntitySummary(zarPath, "SYNTH.ZMX")
            @test summary isa OpticTrace.ZemaxFileSummary
            @test summary.name == "Synthetic Archive Lens"
            @test summary.units == "MM"
        end

        @testset ".zmf entity" begin
            summary = OpticTrace.zemaxEntitySummary(zmfPath, "LENS2")
            @test summary.name == "Synthetic Lens Two"
            @test summary.units == "MM"
        end
    end

    @testset "extractZemaxEntities" begin
        @testset ".zar extract all" begin
            outdir = mktempdir()
            paths = OpticTrace.extractZemaxEntities(zarPath; outputPath = outdir)
            @test length(paths) == 2
            @test all(isfile, paths)
            nestedDir = joinpath(outdir, splitext(basename(zarPath))[1])
            @test all(p -> dirname(p) == nestedDir, paths)
        end

        @testset ".zar extract named" begin
            outdir = mktempdir()
            paths = OpticTrace.extractZemaxEntities(zarPath, ["SYNTH.ZMX"]; outputPath = outdir)
            @test paths == [joinpath(outdir, "SYNTH.ZMX")]
        end

        @testset ".zmf extract all" begin
            outdir = mktempdir()
            paths = OpticTrace.extractZemaxEntities(zmfPath; outputPath = outdir)
            @test length(paths) == 2
            @test all(isfile, paths)
            nestedDir = joinpath(outdir, splitext(basename(zmfPath))[1])
            @test all(p -> dirname(p) == nestedDir, paths)
        end

        @testset ".zmf extract named" begin
            outdir = mktempdir()
            paths = OpticTrace.extractZemaxEntities(zmfPath, ["LENS1"]; outputPath = outdir)
            @test paths == [joinpath(outdir, "LENS1.zmx")]
        end

        @test_throws ErrorException OpticTrace.extractZemaxEntities(zmxPath)
    end

    @testset "defaultExtractionOutputPath" begin
        @test OpticTrace.defaultExtractionOutputPath(zarPath) == dirname(zarPath)
        @test OpticTrace.defaultExtractionOutputPath(zmfPath) == dirname(zmfPath)
        @test_throws ErrorException OpticTrace.defaultExtractionOutputPath(zmxPath)
    end

    @testset "Bonito UI smoke tests" begin
        import Bonito

        function renderToString(nodeFn)
            app = Bonito.App() do session
                nodeFn()
            end
            path = joinpath(mktempdir(), "out.html")
            Bonito.export_static(path, app)
            return read(path, String)
        end

        @testset "zemaxBrowserApp" begin
            app = OpticTrace.zemaxBrowserApp(fixturesDir)
            @test app isa Bonito.App
            path = joinpath(mktempdir(), "out.html")
            Bonito.export_static(path, app)
            html = read(path, String)
            @test occursin("test_singlet.zmx", html)
            @test occursin("test_archive.zar", html)
            @test occursin("test_catalog.zmf", html)
            @test occursin("Select a file", html)
            @test occursin("↑", html)
            @test occursin(">File<", html)
        end

        @testset "contentPane .zmx" begin
            html = renderToString(() -> Bonito.DOM.div(OpticTrace.contentPane(zmxPath)...))
            @test occursin("Test Singlet", html)
            @test occursin("EVENASPH", html)
            @test occursin("TESTGLASS", html)
        end

        @testset "contentPane .zar" begin
            html = renderToString(() -> Bonito.DOM.div(OpticTrace.contentPane(zarPath)...))
            @test occursin("SYNTH.ZMX", html)
            @test occursin("SYNTH.AGF", html)
            @test occursin("no preview", html)
            @test occursin("Output path", html)
            @test occursin("Browse", html)
            @test occursin("Extract All", html)
        end

        @testset "contentPane .zmf" begin
            html = renderToString(() -> Bonito.DOM.div(OpticTrace.contentPane(zmfPath)...))
            @test occursin("LENS1.zmx", html)
            @test occursin("LENS2.zmx", html)
            @test occursin("Extract All", html)
        end

        @testset "contentPane empty selection" begin
            html = renderToString(() -> Bonito.DOM.div(OpticTrace.contentPane("")...))
            @test occursin("Select a file", html)
        end

        @testset "renderTextPreview" begin
            bytes = OpticTrace.zemaxEntityBytes(zarPath, "SYNTH.AGF")
            html = renderToString(() -> OpticTrace.renderTextPreview("SYNTH.AGF", bytes))
            @test occursin("SYNTH.AGF", html)
            @test occursin("synthetic glass catalog stub", html)
        end
    end

end
