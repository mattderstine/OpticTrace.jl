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

        @testset "dispform 2 with negative L (TODO.md #16)" begin
            f = OpticTrace.getAGFRefractiveIndexFunc(2, [0.5, -0.1, 0.0, 1.0, 0.0, 1.0])
            @test f !== nothing
            λ = 0.5875618
            @test f(λ) ≈ sqrt(1 + 0.5 * λ^2 / (λ^2 + 0.1)) atol = 1e-12
        end

        @testset "dispform 5 (Conrady, invented coefficients)" begin
            n0, A, B = 1.45, 0.0125, 7.0e-5
            f = OpticTrace.getAGFRefractiveIndexFunc(5, [n0, A, B, 0.0, 0.0, 0.0])
            @test f !== nothing
            # Conrady written out directly, independent of riFormula5
            for λ in (0.4861327, 0.5875618, 0.6562725, 1.5)
                @test f(λ) ≈ n0 + A / λ + B / λ^3.5 atol = 1e-12
            end
        end

        # Codes 3, 6-9, 12, 13 (TODO.md #37): real glasses from the local
        # Zemax Glasscat, checked against the nd on each glass's own NM
        # line (nd at 0.5875618 um), which is independent of this code.
        # A sweep of the whole local Glasscat/ agrees to ~1e-6 for these.
        @testset "dispform 3 (Herzberger, NHG IRG104)" begin
            c = [1.64425684, 0.00697248348, -0.000295666238, -0.00225530604, -8.6130679e-6, -2.0302048e-7]
            @test OpticTrace.getAGFRefractiveIndexFunc(3, c)(0.5875618) ≈ 1.662519 atol = 1e-5
        end

        @testset "dispform 6 (Sellmeier3, NIKON 5742)" begin
            c = [1.30326513, 0.00977713845, 0.127313632, 0.0463756328, 0.895143157, 107.962139, 0.0, 0.0]
            @test OpticTrace.getAGFRefractiveIndexFunc(6, c)(0.5875618) ≈ 1.576531 atol = 1e-5
        end

        @testset "dispform 7 (Handbook of Optics 1, BBO)" begin
            c = [2.7405, 0.0184, 0.0179, 0.0155]
            @test OpticTrace.getAGFRefractiveIndexFunc(7, c)(0.5875618) ≈ 1.670737 atol = 1e-5
        end

        @testset "dispform 8 (Handbook of Optics 2, ZnO-e)" begin
            c = [2.80333, 0.9447, 0.09024016, 0.00714]
            @test OpticTrace.getAGFRefractiveIndexFunc(8, c)(0.5875618) ≈ 2.019875 atol = 1e-5
        end

        @testset "dispform 9 (Sellmeier4, LiYF4)" begin
            c = [1.38757, 0.70757, 0.00931, 0.18849, 50.99741]
            @test OpticTrace.getAGFRefractiveIndexFunc(9, c)(0.5875618) ≈ 1.453776 atol = 1e-5
        end

        @testset "dispform 12 (Extended2, NIKON-HIKARI E-LAFH3)" begin
            c = [3.10091582, -0.0131944552, 0.0377063856, 0.00215246343, -0.000117361655, 1.96473108e-5, -2.2411102e-5, 0.0]
            @test OpticTrace.getAGFRefractiveIndexFunc(12, c)(0.5875618) ≈ 1.795041 atol = 1e-5
        end

        @testset "dispform 13 (Extended3, NIKON-HIKARI J-LASFH17HS)" begin
            c = [3.81071676, -0.0163737936, 0.0, 0.0583672875, 0.00412726108, -0.000552229126, 0.000141302816, -1.45517862e-5, 7.444268e-7]
            @test OpticTrace.getAGFRefractiveIndexFunc(13, c)(0.5875618) ≈ 2.00069 atol = 1e-5
        end

        # Codes 4, 10 and 11 have no (or only a placeholder-nd) glass in the
        # local Glasscat, so these write each formula out directly with
        # invented coefficients; KRS5 (Handbook of Optics) for code 11 is
        # also physically plausible (n ~ 2.4 at 1 um).
        @testset "dispform 4 (Sellmeier2, unverified against real data)" begin
            A, B1, L1, B2, L2 = 0.01, 1.2, 0.01, 0.8, 100.0
            f = OpticTrace.getAGFRefractiveIndexFunc(4, [A, B1, L1, B2, L2, 0.0])
            for λ in (0.5, 0.5875618, 2.0)
                @test f(λ) ≈ sqrt(1 + A + B1 * λ^2 / (λ^2 - L1) + B2 / (λ^2 - L2)) atol = 1e-12
            end
        end

        @testset "dispform 10 (Extended, unverified against real data)" begin
            a = [2.5, -0.01, 0.02, 0.001, -1e-4, 1e-5, -1e-6, 1e-7]
            f = OpticTrace.getAGFRefractiveIndexFunc(10, a)
            λ = 0.5875618
            n2 = a[1] + a[2] * λ^2 + a[3] * λ^-2 + a[4] * λ^-4 + a[5] * λ^-6 + a[6] * λ^-8 + a[7] * λ^-10 + a[8] * λ^-12
            @test f(λ) ≈ sqrt(n2) atol = 1e-12
        end

        @testset "dispform 11 (Sellmeier5, KRS5)" begin
            c = [1.8293958, 0.0225, 1.6675593, 0.0625, 1.1210424, 0.1225, 0.04513366, 0.2025, 12.380234, 27089.681]
            f = OpticTrace.getAGFRefractiveIndexFunc(11, c)
            λ = 1.0
            n2 = 1 + sum(c[i] * λ^2 / (λ^2 - c[i+1]) for i in 1:2:9)
            @test f(λ) ≈ sqrt(n2) atol = 1e-12
            @test 2.3 < f(λ) < 2.6
        end

        @testset "unsupported dispform returns nothing" begin
            @test OpticTrace.getAGFRefractiveIndexFunc(99, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]) === nothing
        end
    end

    @testset "readAGFRecords (fixture file)" begin
        records = OpticTrace.readAGFRecords(AGF_FIXTURE)
        @test length(records) == 6

        baf = only(filter(r -> r.name == "BAF11TEST", records))
        @test baf.dispform == 1
        @test baf.coefficients ≈ [2.7272955, -0.013470236, 0.014809271, 0.0019615335, -0.00019787663, 1.0851117e-05, 0.0, 0.0, 0.0, 0.0]

        bk7 = only(filter(r -> r.name == "NBK7TEST", records))
        @test bk7.dispform == 2
        @test bk7.coefficients ≈ [1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653, 0.0, 0.0, 0.0, 0.0]

        conrady = only(filter(r -> r.name == "CONRADYTEST", records))
        @test conrady.dispform == 5
        @test conrady.coefficients[1:3] ≈ [1.45, 0.0125, 7.0e-5]

        # dispform written as a float is read as the equivalent Int; a
        # non-integer one is skipped by readAGFRecords
        flt = only(filter(r -> r.name == "FLOATCODETEST", records))
        @test flt.dispform == 1
        @test flt.dispform isa Int
        @test isempty(filter(r -> r.name == "FRACCODETEST", records))

        herz = only(filter(r -> r.name == "HERZTEST", records))
        @test herz.dispform == 99
    end

    @testset "loadAGFCatalog! (into a caller-supplied dict)" begin
        testDict = Dict{AbstractString,Any}()
        testDict["DEFAULT"] = rInDef
        loadAGFCatalog!(testDict, AGF_FIXTURE)

        @test haskey(testDict, "BAF11TEST")
        @test testDict["BAF11TEST"](0.5875618) ≈ 1.66672 atol = 1e-4
        # negative-L Sellmeier1 glass sits before NBK7TEST in the file and
        # must not abort loading (TODO.md #16)
        @test haskey(testDict, "NEGLTEST")
        @test haskey(testDict, "NBK7TEST")
        @test testDict["NBK7TEST"](0.5875618) ≈ 1.5168 atol = 1e-4
        @test haskey(testDict, "CONRADYTEST")
        @test testDict["CONRADYTEST"](0.5875618) ≈ 1.45 + 0.0125 / 0.5875618 + 7.0e-5 / 0.5875618^3.5 atol = 1e-12
        @test !haskey(testDict, "HERZTEST") # unsupported dispform, skipped
        @test haskey(testDict, "FLOATCODETEST")
        @test testDict["FLOATCODETEST"](0.5875618) ≈ 1.66672 atol = 1e-4
        @test !haskey(testDict, "FRACCODETEST")
    end

    @testset "loadAGFCatalog (fresh dict)" begin
        cat = loadAGFCatalog(AGF_FIXTURE)
        @test haskey(cat, "DEFAULT")
        @test cat["DEFAULT"] === rInDef
        @test haskey(cat, "BAF11TEST")
        @test haskey(cat, "NBK7TEST")
        @test haskey(cat, "CONRADYTEST")
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

    if !HAS_ZEMAX_GLASSCAT
        @info "Skipping agf.jl real-catalog Conrady cross-check: $ZEMAX_GLASSCAT_DIR not found"
    else
        @testset "Conrady glasses reproduce their catalog nd (real DOWSIL_MS.AGF, read in place)" begin
            # Each NM line's 5th token is the catalog's own nd -- an
            # independent reference for the Conrady mapping onto riFormula5
            agfPath = joinpath(ZEMAX_GLASSCAT_DIR, "DOWSIL_MS.AGF")
            testDict = Dict{AbstractString,Any}()
            loadAGFCatalog!(testDict, agfPath)
            checked = 0
            for line in eachline(agfPath)
                tokens = split(filter(!=('\0'), line))
                (length(tokens) >= 5 && tokens[1] == "NM" && tryparse(Float64, tokens[3]) == 5.0) || continue
                @test testDict[tokens[2]](0.5875618) ≈ parse(Float64, tokens[5]) atol = 1e-6
                checked += 1
            end
            @test checked > 0
        end

        @testset "codes 3, 6-9, 12, 13 reproduce their catalog nd (whole local Glasscat, TODO.md #37)" begin
            # Glasses with a placeholder nd (0 or 1.0, e.g. infrared
            # materials) can't be cross-checked this way and are skipped.
            codes = (3, 6, 7, 8, 9, 12, 13)
            checked = Dict(code => 0 for code in codes)
            for agfFile in readdir(ZEMAX_GLASSCAT_DIR; join=true)
                lowercase(splitext(agfFile)[2]) == ".agf" || continue
                catalogNd = Dict{String,Float64}()
                for line in eachline(agfFile)
                    tokens = split(filter(!=('\0'), line))
                    (length(tokens) >= 5 && tokens[1] == "NM") || continue
                    nd = tryparse(Float64, tokens[5])
                    nd === nothing || (catalogNd[tokens[2]] = nd)
                end
                for record in OpticTrace.readAGFRecords(agfFile)
                    record.dispform in codes || continue
                    nd = get(catalogNd, record.name, 0.0)
                    (nd == 0.0 || nd == 1.0) && continue
                    f = OpticTrace.getAGFRefractiveIndexFunc(record.dispform, record.coefficients)
                    @test f(0.5875618) ≈ nd atol = 1e-5
                    checked[record.dispform] += 1
                end
            end
            @test all(>(0), values(checked))
        end
    end

end
