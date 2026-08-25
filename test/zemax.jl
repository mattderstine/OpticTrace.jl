@testset "zemax.jl" begin

    fixture = joinpath(@__DIR__, "fixtures", "test_singlet.zmx")

    @testset "readZemax" begin
        zsurfs, name, units, wavelengths = readZemax(fixture)

        @test name == "Test Singlet"
        @test units == "MM"
        @test length(wavelengths) == 24
        @test wavelengths[1] == 0.5875618
        @test all(wavelengths[2:end] .== 0.0)

        @test length(zsurfs) == 4

        # SURF 0 -- object surface, defaults throughout (no TYPE/GLAS/DIAM/COMM)
        s0 = zsurfs[1]
        @test s0.curvature == 0.0
        @test s0.distance == Inf
        @test s0.material == "DEFAULT"
        @test s0.radius == 0.0
        @test s0.stop == false
        @test s0.conic == 0.0
        @test s0.type == "STANDARD"
        @test s0.comm == ""

        # SURF 1 -- stop, STANDARD, glass TESTGLASS
        s1 = zsurfs[2]
        @test s1.curvature == 0.05
        @test s1.distance == 5.0
        @test s1.material == "TESTGLASS"
        @test s1.radius == 10.0
        @test s1.stop == true
        @test s1.conic == 0.0
        @test s1.coating == "TESTCOAT"
        @test s1.type == "STANDARD"
        @test s1.comm == "Front"

        # SURF 2 -- EVENASPH, conic + two PARM entries
        s2 = zsurfs[3]
        @test s2.curvature == -0.04
        @test s2.distance == 20.0
        @test s2.material == "DEFAULT" # no GLAS line on this surface
        @test s2.radius == 10.0
        @test s2.stop == false
        @test s2.conic == -0.5
        @test s2.coating == "TESTCOAT"
        @test s2.type == "EVENASPH"
        @test s2.comm == "Back"
        @test s2.aspherics[1] == 0.001  # PARM 1
        @test s2.aspherics[2] == 0.0002 # PARM 2
        @test all(s2.aspherics[3:end] .== 0.0)

        # SURF 3 -- final flat surface, DISZ 0
        s3 = zsurfs[4]
        @test s3.curvature == 0.0
        @test s3.distance == 0.0
        @test s3.type == "STANDARD"

        @testset "basept/dir kwargs accepted but not used during parsing (see TODO.md)" begin
            zsurfsA, = readZemax(fixture; basept = ORIGIN, dir = ZAXIS)
            zsurfsB, = readZemax(fixture; basept = Point3(5.0, -3.0, 2.0), dir = YAXIS)
            @test zsurfsA[2].curvature == zsurfsB[2].curvature
            @test zsurfsA[3].aspherics == zsurfsB[3].aspherics
        end
    end

    @testset "printZemaxSurfs" begin
        zsurfs, = readZemax(fixture)
        # smoke test only -- just confirm it runs without erroring over a
        # real parsed surface vector
        redirect_stdout(devnull) do
            printZemaxSurfs(zsurfs)
        end
        @test true
    end

    @testset "zemaxsurfToSurface" begin

        base = ORIGIN
        dir = ZAXIS
        rinIn = 1.0
        rinOut = 1.5

        @testset "STANDARD" begin
            s = OpticTrace.ZemaxSurf(0.03, 4.0, "TESTGLASS", 6.0, true, 0.2, zeros(Float64, 20), "coat1", "STANDARD", "S")
            newbase, newrinOut, surf = OpticTrace.zemaxsurfToSurface(1, rinIn, rinOut, base, dir, s)

            @test newbase ≈ base + s.distance .* dir
            @test newrinOut == rinOut
            @test surf.profile isa OpticTrace.SurfProfileConic
            @test surf.profile.curv == s.curvature
            @test surf.profile.ϵ == OpticTrace.conicToϵ(s.conic)
            @test surf.mod isa DielectricT
            @test surf.mod.refIndexIn == rinIn
            @test surf.mod.refIndexOut == rinOut
            @test surf.aperture.semiDiameter == s.radius
            @test surf.base.base == base
            @test surf.base.dir == dir
            @test surf.surfname == "S 1"
        end

        @testset "EVENASPH" begin
            asph = zeros(Float64, 20)
            asph[2] = 0.0007
            asph[3] = -0.00002
            s = OpticTrace.ZemaxSurf(-0.02, 8.0, "DEFAULT", 5.0, false, -0.6, asph, "coat2", "EVENASPH", "A")
            newbase, newrinOut, surf = OpticTrace.zemaxsurfToSurface(2, rinIn, rinOut, base, dir, s)

            @test newbase ≈ base + s.distance .* dir
            @test surf.profile isa OpticTrace.SurfProfileEvenAsphere
            @test surf.profile.curv == s.curvature
            @test surf.profile.ϵ == OpticTrace.conicToϵ(s.conic)
            @test surf.profile.a == s.aspherics[2:end]
            @test surf.surfname == "A 2"
        end

        @testset "unsupported surface type" begin
            s = OpticTrace.ZemaxSurf(0.0, 1.0, "DEFAULT", 5.0, false, 0.0, zeros(Float64, 20), "", "TORIC", "")
            @test_throws ErrorException OpticTrace.zemaxsurfToSurface(1, rinIn, rinOut, base, dir, s)
        end
    end

    @testset "zemaxsurfsToGeo" begin
        zsurfs, = readZemax(fixture)

        testCatalog = Dict{AbstractString,Any}(
            "DEFAULT" => (λ -> 1.0),
            "TESTGLASS" => (λ -> 1.6),
        )

        geo = zemaxsurfsToGeo(zsurfs[2:end], ORIGIN, ZAXIS, 0.5875618; glassCatalog = testCatalog)

        @test length(geo) == 3

        # surface 1: object -> TESTGLASS
        @test geo[1].base.base == ORIGIN
        @test geo[1].mod.refIndexIn == refIndexDefault
        @test geo[1].mod.refIndexOut == 1.6
        @test geo[1].profile isa OpticTrace.SurfProfileConic

        # surface 2: TESTGLASS -> DEFAULT, aspheric
        @test geo[2].base.base ≈ ORIGIN + 5.0 .* ZAXIS
        @test geo[2].mod.refIndexIn == 1.6
        @test geo[2].mod.refIndexOut == 1.0
        @test geo[2].profile isa OpticTrace.SurfProfileEvenAsphere

        # surface 3: DEFAULT -> DEFAULT, flat, zero distance
        @test geo[3].base.base ≈ ORIGIN + 25.0 .* ZAXIS
        @test geo[3].mod.refIndexIn == 1.0
        @test geo[3].mod.refIndexOut == 1.0
    end

end
