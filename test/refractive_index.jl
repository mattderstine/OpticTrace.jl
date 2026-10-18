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

        @testset "riFormula5 (Cauchy form, invented coefficients)" begin
            # n = c1 + c2*λ^c3 + c4*λ^c5, written out directly
            c = [1.45, 0.0125, -1.0, 7.0e-5, -3.5]
            for λ in (0.4, 0.5875618, 1.2)
                @test riFormula5(λ, c) ≈ 1.45 + 0.0125 / λ + 7.0e-5 / λ^3.5 atol = 1e-12
            end
            @test riFormula5(0.5, [1.3]) == 1.3 # constant term only
            # λ and c needn't share an element type (unlike riFormula1/2/3)
            @test riFormula5(0.5, [1, 2, -1]) ≈ 5.0
            @test ForwardDiff.derivative(λ -> riFormula5(λ, c), 0.6) ≈ -0.0125 / 0.36 - 3.5 * 7.0e-5 / 0.6^4.5 atol = 1e-10
        end

        @testset "ModelGlass" begin
            g = ModelGlass(1.5168, 64.17, 0.0)

            # OpticSim.jl's own documented value for its modelglass
            @test g(0.5875618) ≈ 1.5168003970108495 atol = 1e-12

            # the fit reproduces nd, and approximately vd, for glasses
            # across the glass map (vd from the F/C lines)
            for (nd, vd, dPgF) in ((1.5168, 64.17, 0.0), (1.5598, 42.33, 0.0), (1.7847, 26.10, 0.014), (1.62, 56.9, -0.0003))
                mg = ModelGlass(nd, vd, dPgF)
                @test mg(0.5875618) ≈ nd atol = 1e-4
                @test (nd - 1) / (mg(0.4861327) - mg(0.6562725)) ≈ vd rtol = 5e-3
            end

            # normal dispersion: index falls with wavelength
            @test g(0.45) > g(0.55) > g(0.65)

            @test g.nd == 1.5168 && g.vd == 64.17 && g.dPgF == 0.0
            @test ModelGlass(1.5168, 64.17) == g # dPgF defaults to 0
            @test ModelGlass(1.5f0, 60, 0.0) isa ModelGlass{Float64} # promoted

            @testset "autodiff" begin
                @test ForwardDiff.derivative(g, 0.55) < 0
                # dual nd with Float64 λ -- needs riFormula5's independent type parameters
                dndnd = ForwardDiff.derivative(nd -> ModelGlass(nd, 64.17, 0.0)(0.5875618), 1.5168)
                @test dndnd ≈ 1.0 atol = 1e-2
            end

            @testset "input range checks" begin
                @test_throws ArgumentError ModelGlass(0.9, 50.0, 0.0)
                @test_throws ArgumentError ModelGlass(2.6, 50.0, 0.0)
                @test_throws ArgumentError ModelGlass(1.5, 5.0, 0.0)
                @test_throws ArgumentError ModelGlass(1.5, 120.0, 0.0)
            end

            @testset "warns outside the 0.36-0.75 μm fit range" begin
                @test_logs (:warn, r"fit range") ModelGlass(1.5, 60.0)(1.0)
                @test_logs ModelGlass(1.5, 60.0)(0.5) # no warning inside the range
            end
        end

        # Formulas 4, 6-9 (TODO.md #15): expected values written out
        # independently from each formula's definition (not via the
        # code under test), using coefficients from real records.
        @testset "riFormula4 (KAP Moutzouris-gamma, 11 of 17 coefficients)" begin
            c = [2.196073191, 0.015025722, 0, 0.190981743, 2, 0, 0, 0, 1, -0.006100780, 2]
            λ = 0.5875618
            n2 = 2.196073191 + 0.015025722 / (λ^2 - 0.190981743^2) - 0.006100780 * λ^2
            @test OpticTrace.riFormula4(λ, c) ≈ sqrt(n2) atol = 1e-12
            @test OpticTrace.riFormula4(λ, c) ≈ 1.497542175829694 atol = 1e-9
            # explicit zero padding gives the same answer
            @test OpticTrace.riFormula4(λ, [c; zeros(6)]) == OpticTrace.riFormula4(λ, c)
        end

        @testset "riFormula6 (He, Mansfield)" begin
            c = [0, 0.01470091, 423.98]
            @test OpticTrace.riFormula6(0.5875618, c) ≈ 1.0000349121124632 atol = 1e-12
        end

        @testset "riFormula7 (Si, Edwards; c6 omitted)" begin
            c = [3.41983, 0.159906, -0.123109, 1.26878e-6, -1.95104e-9]
            @test OpticTrace.riFormula7(3.0, c) ≈ 3.436134677527718 atol = 1e-9
        end

        @testset "riFormula8 (AgBr, Schroter)" begin
            c = [0.452505, 0.09939, 0.070537, -0.000150]
            @test OpticTrace.riFormula8(0.55, c) ≈ 2.275584479878346 atol = 1e-9
        end

        @testset "riFormula9 (urea, Rosker-e)" begin
            c = [2.51527, 0.0240, 0.0300, 0.020, 1.52, 0.8771]
            @test OpticTrace.riFormula9(0.6, c) ≈ 1.605403788031452 atol = 1e-9
        end

        @testset "getRefractiveIndexFunc unsupported formula (TODO.md #15)" begin
            mktempdir() do dir
                write(joinpath(dir, "bad.yml"), "DATA:\n  - type: formula 99\n    coefficients: 1 2 3\n")
                write(joinpath(dir, "good.yml"), "DATA:\n  - type: formula 6\n    coefficients: 0 0.01470091 423.98\n")
                write(joinpath(dir, "tab.yml"), "DATA:\n  - type: tabulated n\n    data: |\n      0.5 1.5\n")
                # skipped at load time, not a KeyError on first evaluation
                @test getRefractiveIndexFunc(dir, "bad.yml") === nothing
                @test getRefractiveIndexFunc(dir, "tab.yml") === nothing
                f = getRefractiveIndexFunc(dir, "good.yml")
                @test f(0.5875618) ≈ 1.0000349121124632 atol = 1e-12
                cat = Dict{AbstractString,Any}()
                loadRICatalog!(cat, "."; basepath=dir)
                @test haskey(cat, "good")
                @test !haskey(cat, "bad")
            end
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

        @testset "getRefractiveIndexFunc, \"formula 5\" (Cauchy) record" begin
            # glass/misc/soda-lime/Rubin-lowiron.yml: "formula 5" with
            # coefficients 1.5130 -0.003169 2 0.003962 -2, written out here
            f = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/misc/soda-lime/Rubin-lowiron.yml")
            λ = 0.5875618
            @test f(λ) ≈ 1.5130 - 0.003169 * λ^2 + 0.003962 / λ^2 atol = 1e-12
        end

        @testset "formulas 4, 6-9 real records (TODO.md #15)" begin
            base = OpticTrace.dirBaseRefractiveIndex
            f = getRefractiveIndexFunc(base, "organic/C8H5KO4 - potassium hydrogen phthalate/Moutzouris-gamma.yml")
            @test f(0.5875618) ≈ 1.497542175829694 atol = 1e-9
            f = getRefractiveIndexFunc(base, "main/He/Mansfield.yml")
            @test f(0.5875618) ≈ 1.0000349121124632 atol = 1e-12
            f = getRefractiveIndexFunc(base, "main/Si/Edwards.yml")
            @test f(3.0) ≈ 3.436134677527718 atol = 1e-9
            f = getRefractiveIndexFunc(base, "main/AgBr/Schroter.yml")
            @test f(0.55) ≈ 2.275584479878346 atol = 1e-9
            f = getRefractiveIndexFunc(base, "organic/CH4N2O - urea/Rosker-e.yml")
            @test f(0.6) ≈ 1.605403788031452 atol = 1e-9
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
