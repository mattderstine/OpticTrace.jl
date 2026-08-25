@testset "refractive_index.jl" begin

    # Expected index values below are cross-checked against each glass
    # YAML file's own SPECS.nd entry where available (an independently
    # listed reference value, sourced -- per each file's REFERENCES field
    # -- from the manufacturer's own catalog via refractiveindex.info,
    # not derived from our own dispersion-formula code), or otherwise
    # computed directly from the file's own coefficients when no such
    # independent reference value is available (documented per case
    # below). Wavelengths are in microns, matching this codebase's
    # convention throughout.

    @testset "Dispersion formulas (pure math, no file I/O)" begin

        @testset "riFormula1 (Schott infrared/schott/infrared/IRG23.yml)" begin
            # IRG23 is an infrared-only glass (valid range 1.0-12 microns)
            # with no visible-range nd in its SPECS, so there's no
            # independent reference value to cross-check against here --
            # this is a direct, precision-computed value from the file's
            # own "formula 1" coefficients.
            c = [3.74971, 3.07065, 0.498392, 0.953206, 40.7557]
            @test OpticTrace.riFormula1(2.0, c) ≈ 2.832201860963834 atol = 1e-9
            @test OpticTrace.riFormula1(5.0, c) ≈ 2.799394510278137 atol = 1e-9
        end

        @testset "riFormula2 (Schott N-BK7, schott/N-BK7.yml)" begin
            # N-BK7 at the d-line (587.5618 nm): SPECS lists nd = 1.5168,
            # the textbook value for this glass, from Schott's own catalog.
            c = [0.0, 1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653]
            @test OpticTrace.riFormula2(0.5875618, c) ≈ 1.5168 atol = 1e-4
        end

        @testset "riFormula3 (Hoya BAF11/E-FD10, glass/hoya/*.yml)" begin
            # BAF11: SPECS lists nd = 1.66672 (Hoya catalog via
            # refractiveindex.info).
            cBAF11 = [2.7272955, -0.013470236, 2.0, 0.014809271, -2.0, 0.0019615335, -4.0, -0.00019787663, -6.0, 1.0851117e-05, -8.0]
            @test OpticTrace.riFormula3(0.5875618, cBAF11) ≈ 1.66672 atol = 1e-4

            # E-FD10: SPECS lists nd = 1.72825.
            cEFD10 = [2.881518, -0.013228312, 2.0, 0.03145559, -2.0, 0.0026851666, -4.0, -0.00022577544, -6.0, 2.4693268e-05, -8.0]
            @test OpticTrace.riFormula3(0.5875618, cEFD10) ≈ 1.72825 atol = 1e-4
        end
    end

    if !HAS_GLASS_CATALOG
        @info "Skipping refractive_index.jl catalog-access tests: OpticTrace.dirBaseRefractiveIndex ($(OpticTrace.dirBaseRefractiveIndex)) not found on this machine"
    else
        @testset "Catalog access (needs the real glass directory)" begin

        @testset "getRefractiveIndexFunc" begin
            f = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/hoya/BAF11.yml")
            @test f !== nothing
            @test f(0.5875618) ≈ 1.66672 atol = 1e-4

            f2 = getRefractiveIndexFunc("glass/hoya/BAF11.yml") # basepath defaults to dirBaseRefractiveIndex
            @test f2(0.5875618) ≈ 1.66672 atol = 1e-4
        end

        @testset "loadRICatalog! (into a caller-supplied dict)" begin
            testDict = Dict{AbstractString,Any}()
            testDict["DEFAULT"] = rInDef
            loadRICatalog!(testDict, "glass/hoya"; basepath=OpticTrace.dirBaseRefractiveIndex)

            @test haskey(testDict, "BAF11")
            @test testDict["BAF11"](0.5875618) ≈ 1.66672 atol = 1e-4
            @test haskey(testDict, "E-FD10")
            @test haskey(testDict, "FD10") # dash-split alias: "E-FD10" -> also keyed "FD10"
        end

        @testset "loadRICatalog (fresh dict)" begin
            cat = loadRICatalog("glass/hoya"; basepath=OpticTrace.dirBaseRefractiveIndex)
            @test haskey(cat, "DEFAULT")
            @test cat["DEFAULT"] === rInDef
            @test haskey(cat, "BAF11")
        end

        @testset "loadRICatalog! (into the global defaultGlassCatalog)" begin
            OpticTrace.loadRICatalog!("glass/hoya"; basepath=OpticTrace.dirBaseRefractiveIndex)
            @test haskey(OpticTrace.defaultGlassCatalog, "BAF11")
            @test OpticTrace.defaultGlassCatalog["BAF11"](0.5875618) ≈ 1.66672 atol = 1e-4
        end

        @testset "setDirectoryBaseRefractiveIndex / getDirectoryBaseRI" begin
            original = getDirectoryBaseRI()
            try
                setDirectoryBaseRefractiveIndex("/some/custom/path")
                @test getDirectoryBaseRI() == "/some/custom/path"

                setDirectoryBaseRefractiveIndex() # 0-arg: resets to the hardcoded default
                # ties to the hardcoded-path TODO item -- this literal
                # string may need updating once a config-file mechanism
                # replaces it.
                @test getDirectoryBaseRI() == "/Users/matt/Development/Projects/refractiveindex/database/data"
            finally
                setDirectoryBaseRefractiveIndex(original)
            end
            @test getDirectoryBaseRI() == original
        end
    end
    end # if HAS_GLASS_CATALOG

end
