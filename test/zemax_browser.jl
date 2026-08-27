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

    @testset "walkZemaxDirectory" begin
        tmp = mktempdir()
        mkpath(joinpath(tmp, "sub"))
        write(joinpath(tmp, "a.zmx"), "")
        write(joinpath(tmp, "b.zar"), "")
        write(joinpath(tmp, "ignored.txt"), "")
        write(joinpath(tmp, "sub", "c.zmf"), "")

        node = OpticTrace.walkZemaxDirectory(tmp)
        @test node.isDir
        @test node.kind == :dir
        @test node.path == tmp

        names = [c.name for c in node.children]
        @test "a.zmx" in names
        @test "b.zar" in names
        @test "sub" in names
        @test !("ignored.txt" in names)

        aNode = node.children[findfirst(==("a.zmx"), names)]
        @test aNode.kind == :zmx
        @test !aNode.isDir
        @test isempty(aNode.children)

        subNode = node.children[findfirst(==("sub"), names)]
        @test subNode.isDir
        @test length(subNode.children) == 1
        @test subNode.children[1].name == "c.zmf"
        @test subNode.children[1].kind == :zmf

        @testset "on test/fixtures itself" begin
            fixturesNode = OpticTrace.walkZemaxDirectory(fixturesDir)
            fixtureNames = [c.name for c in fixturesNode.children]
            @test "test_singlet.zmx" in fixtureNames
            @test "test_archive.zar" in fixtureNames
            @test "test_catalog.zmf" in fixtureNames
            # generate_synthetic_zemax.jl has no browsable extension
            @test !("generate_synthetic_zemax.jl" in fixtureNames)
        end
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
        end

        @testset ".zmf extract named" begin
            outdir = mktempdir()
            paths = OpticTrace.extractZemaxEntities(zmfPath, ["LENS1"]; outputPath = outdir)
            @test paths == [joinpath(outdir, "LENS1.zmx")]
        end

        @test_throws ErrorException OpticTrace.extractZemaxEntities(zmxPath)
    end

    @testset "defaultExtractionOutputPath" begin
        @test OpticTrace.defaultExtractionOutputPath(zarPath) == OpticTrace._stripExtension(zarPath, ".zar")
        @test OpticTrace.defaultExtractionOutputPath(zmfPath) == OpticTrace._stripExtension(zmfPath, ".zmf")
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

        @testset "_zemaxBrowserApp" begin
            app = OpticTrace._zemaxBrowserApp(fixturesDir)
            @test app isa Bonito.App
            path = joinpath(mktempdir(), "out.html")
            Bonito.export_static(path, app)
            html = read(path, String)
            @test occursin("test_singlet.zmx", html)
            @test occursin("test_archive.zar", html)
            @test occursin("test_catalog.zmf", html)
            @test occursin("Select a file", html)
        end

        @testset "_contentPane .zmx" begin
            html = renderToString(() -> OpticTrace._contentPane(zmxPath))
            @test occursin("Test Singlet", html)
            @test occursin("EVENASPH", html)
            @test occursin("TESTGLASS", html)
        end

        @testset "_contentPane .zar" begin
            html = renderToString(() -> OpticTrace._contentPane(zarPath))
            @test occursin("SYNTH.ZMX", html)
            @test occursin("SYNTH.AGF", html)
            @test occursin("no preview", html)
            @test occursin("Output path", html)
            @test occursin("Extract All", html)
        end

        @testset "_contentPane .zmf" begin
            html = renderToString(() -> OpticTrace._contentPane(zmfPath))
            @test occursin("LENS1", html)
            @test occursin("LENS2", html)
            @test occursin("Extract All", html)
        end

        @testset "_contentPane empty selection" begin
            html = renderToString(() -> OpticTrace._contentPane(""))
            @test occursin("Select a file", html)
        end
    end

end
