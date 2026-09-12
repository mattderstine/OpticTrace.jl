@testset "agf.jl" begin

    # Reference nd values below reuse the same known-good coefficient
    # sets already verified in test/refractive_index.jl: Hoya BAF11's
    # riFormula3 coefficients (nd = 1.66672) recast as AGF Schott
    # coefficients, and N-BK7's public Sellmeier coefficients
    # (nd = 1.5168) recast as AGF Sellmeier1 coefficients.

    AGF_FIXTURE = joinpath(@__DIR__, "fixtures", "test_glasses.agf")

    @testset "getAGFRefractiveIndexFunc (pure math, no file I/O)" begin
        @testset "dispform 1 (Schott, BAF11 coefficients)" begin
            c = [2.7272955, -0.013470236, 0.014809271, 0.0019615335, -0.00019787663, 1.0851117e-05]
            f = OpticTrace.getAGFRefractiveIndexFunc(1, c)
            @test f !== nothing
            @test f(0.5875618) ≈ 1.66672 atol = 1e-4
        end

        @testset "dispform 2 (Sellmeier1, N-BK7 coefficients)" begin
            c = [1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653]
            f = OpticTrace.getAGFRefractiveIndexFunc(2, c)
            @test f !== nothing
            @test f(0.5875618) ≈ 1.5168 atol = 1e-4
        end

        @testset "unsupported dispform returns nothing" begin
            @test OpticTrace.getAGFRefractiveIndexFunc(3, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]) === nothing
        end
    end

    @testset "readAGFRecords (fixture file)" begin
        records = OpticTrace.readAGFRecords(AGF_FIXTURE)
        @test length(records) == 3

        baf = only(filter(r -> r.name == "BAF11TEST", records))
        @test baf.dispform == 1
        @test baf.coefficients ≈ [2.7272955, -0.013470236, 0.014809271, 0.0019615335, -0.00019787663, 1.0851117e-05, 0.0, 0.0, 0.0, 0.0]

        bk7 = only(filter(r -> r.name == "NBK7TEST", records))
        @test bk7.dispform == 2
        @test bk7.coefficients ≈ [1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653, 0.0, 0.0, 0.0, 0.0]

        herz = only(filter(r -> r.name == "HERZTEST", records))
        @test herz.dispform == 3
    end

    @testset "loadAGFCatalog! (into a caller-supplied dict)" begin
        testDict = Dict{AbstractString,Any}()
        testDict["DEFAULT"] = rInDef
        loadAGFCatalog!(testDict, AGF_FIXTURE)

        @test haskey(testDict, "BAF11TEST")
        @test testDict["BAF11TEST"](0.5875618) ≈ 1.66672 atol = 1e-4
        @test haskey(testDict, "NBK7TEST")
        @test testDict["NBK7TEST"](0.5875618) ≈ 1.5168 atol = 1e-4
        @test !haskey(testDict, "HERZTEST") # unsupported dispform, skipped
    end

    @testset "loadAGFCatalog (fresh dict)" begin
        cat = loadAGFCatalog(AGF_FIXTURE)
        @test haskey(cat, "DEFAULT")
        @test cat["DEFAULT"] === rInDef
        @test haskey(cat, "BAF11TEST")
        @test haskey(cat, "NBK7TEST")
        @test !haskey(cat, "HERZTEST")
    end

    @testset "loadAGFCatalog! (into the global defaultGlassCatalog)" begin
        loadAGFCatalog!(AGF_FIXTURE)
        @test haskey(OpticTrace.defaultGlassCatalog, "BAF11TEST")
        @test OpticTrace.defaultGlassCatalog["BAF11TEST"](0.5875618) ≈ 1.66672 atol = 1e-4
        @test haskey(OpticTrace.defaultGlassCatalog, "NBK7TEST")
        @test OpticTrace.defaultGlassCatalog["NBK7TEST"](0.5875618) ≈ 1.5168 atol = 1e-4
    end

    @testset "loadAGFCatalog! over a directory" begin
        testDict = Dict{AbstractString,Any}()
        testDict["DEFAULT"] = rInDef
        loadAGFCatalog!(testDict, joinpath(@__DIR__, "fixtures"))
        @test haskey(testDict, "BAF11TEST")
        @test haskey(testDict, "NBK7TEST")
    end

end
