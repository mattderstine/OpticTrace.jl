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

    @testset "lzwDecompress" begin
        # First 24 bytes of the real "3479-S02.ZMX.LZW" entry's payload
        # inside ZAR_SAMPLE_PATH (see test/helper.jl), decoding to a
        # clean 21-byte ZMX header line -- cross-checked against a
        # from-scratch Python port of the reference decompressor.
        compressed = UInt8[0x2b,0x11,0x4a,0x45,0x31,0x00,0xc4,0x60,0x30,0x1a,
                            0x0c,0x46,0x82,0x01,0x80,0x80,0x64,0x36,0x1c,0x8d,
                            0xc6,0x00,0xd0,0x51]
        expected = UInt8[0x56,0x45,0x52,0x53,0x20,0x31,0x30,0x30,0x34,0x31,
                          0x34,0x20,0x30,0x20,0x32,0x36,0x39,0x37,0x30,0x0d,0x0a]
        @test OpticTrace.lzwDecompress(compressed) == expected
        @test String(OpticTrace.lzwDecompress(compressed)) == "VERS 100414 0 26970\r\n"
    end

    @testset "zmfDeobfuscate" begin
        # First 16 bytes of lens "5002"'s obfuscated description inside
        # ZMF_SAMPLE_PATH (see test/helper.jl), with that lens's real
        # efl/enp -- cross-checked against a from-scratch Python port
        # of rayopt's zmf_obfuscate.
        obfuscated = UInt8[0x9a,0x0c,0xa1,0xb0,0x16,0x51,0xc8,0xe4,
                            0xec,0x89,0x40,0x0a,0xab,0xe5,0x0b,0x87]
        expected = UInt8[0x56,0x45,0x52,0x53,0x20,0x31,0x30,0x30,
                          0x35,0x30,0x33,0x0a,0x4d,0x4f,0x44,0x45]
        @test OpticTrace.zmfDeobfuscate(obfuscated, 4.485, 5.2) == expected
        @test String(OpticTrace.zmfDeobfuscate(obfuscated, 4.485, 5.2)) == "VERS 100503\nMODE"
    end

    if HAS_ZAR_SAMPLE
        @testset "Zemax .zar archive reading" begin
            names = listZemaxArchive(ZAR_SAMPLE_PATH)
            @test names == ["3479-S02.ZMX", "SCHOTT.AGF", "INFRARED.AGF",
                             "MISC.AGF", "COATINGTHOR.DAT"]

            outdir = mktempdir()
            paths = extractZemaxArchive(ZAR_SAMPLE_PATH; outputPath = outdir)
            @test length(paths) == 5
            @test all(isfile, paths)
            nestedDir = joinpath(outdir, splitext(basename(ZAR_SAMPLE_PATH))[1])
            @test all(p -> dirname(p) == nestedDir, paths)

            zmxtext = read(joinpath(nestedDir, "3479-S02.ZMX"), String)
            @test startswith(zmxtext, "VERS 100414 0 26970")
            @test occursin("NAME LF1988 - Negative Meniscus - N-BK7", zmxtext)

            outdir2 = mktempdir()
            selected = extractZemaxArchive(ZAR_SAMPLE_PATH, ["3479-S02.ZMX"]; outputPath = outdir2)
            @test selected == [joinpath(outdir2, "3479-S02.ZMX")]
            @test readdir(outdir2) == ["3479-S02.ZMX"]

            @test_throws ErrorException extractZemaxArchive(ZAR_SAMPLE_PATH, ["NOT_A_REAL_ENTRY"]; outputPath = mktempdir())
        end
    else
        @info "Skipping Zemax .zar archive reading tests: $ZAR_SAMPLE_PATH not found"
    end

    if HAS_ZMF_SAMPLE
        @testset "Zemax .zmf catalog reading" begin
            names = listZmfCatalog(ZMF_SAMPLE_PATH)
            @test names == ["5002", "8003"]

            entries = readZmfCatalog(ZMF_SAMPLE_PATH)
            @test length(entries) == 2
            @test entries[1].name == "5002"
            @test entries[1].efl == 4.485
            @test entries[1].enp == 5.2
            @test entries[1].elements == 1

            outdir = mktempdir()
            paths = extractZmfCatalog(ZMF_SAMPLE_PATH; outputPath = outdir)
            @test length(paths) == 2
            @test all(isfile, paths)
            nestedDir = joinpath(outdir, splitext(basename(ZMF_SAMPLE_PATH))[1])
            @test all(p -> dirname(p) == nestedDir, paths)

            zmxtext = read(joinpath(nestedDir, "5002.zmx"), String)
            @test startswith(zmxtext, "VERS 100503")
            @test occursin("GLAS ACRYLIC", zmxtext)
            @test occursin("CURV 3.752562061747658500E-001", zmxtext)

            outdir2 = mktempdir()
            selected = extractZmfCatalog(ZMF_SAMPLE_PATH, ["8003"]; outputPath = outdir2)
            @test selected == [joinpath(outdir2, "8003.zmx")]
            @test readdir(outdir2) == ["8003.zmx"]

            @test_throws ErrorException extractZmfCatalog(ZMF_SAMPLE_PATH, ["NOT_A_REAL_LENS"]; outputPath = mktempdir())
        end
    else
        @info "Skipping Zemax .zmf catalog reading tests: $ZMF_SAMPLE_PATH not found"
    end

end
