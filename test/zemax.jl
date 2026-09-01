@testset "zemax.jl" begin

    fixture = joinpath(@__DIR__, "fixtures", "test_singlet.zmx")

    @testset "readZemax" begin
        zsurfs, header = readZemax(fixture)

        @test header.name == "Test Singlet"
        @test header.units == "MM"
        @test length(header.wavelengths) == 24
        @test header.wavelengths[1] == 0.5875618
        @test all(header.wavelengths[2:end] .== 0.0)

        @test length(zsurfs) == 4

        @testset "ZemaxHeader -- new header keywords" begin
            @test header.mode == "SEQ"
            @test header.notes == "Test fixture for OpticTrace.jl's Zemax import.\nSecond line of notes."
            @test header.apertureType == "ENPD"
            @test header.apertureValue == 8.0
            @test header.glassCatalogs == ["SCHOTT", "MISC"]
            @test header.primaryWavelengthIndex == 1
            @test header.fieldType == 0
            @test header.fields == [Point2(0.0, 0.0), Point2(5.0, 3.0)]
            @test header.fieldWeight == [1.0, 0.5]
        end

        @testset "ZemaxHeader defaults (PWAV/aperture/GCAT absent)" begin
            defaultHeader = OpticTrace.ZemaxHeader()
            @test defaultHeader.primaryWavelengthIndex == 1
            @test defaultHeader.apertureType == ""
            @test isnan(defaultHeader.apertureValue)
            @test isempty(defaultHeader.glassCatalogs)
            @test isempty(defaultHeader.fields)
            @test defaultHeader.mode == "SEQ"
        end

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

        @testset "bare NAME line (no quoted name) doesn't crash (TODO.md Bug #13)" begin
            bareNameFixture = tempname() * ".zmx"
            write(bareNameFixture, replace(read(fixture, String),
                "NAME \"Test Singlet\"" => "NAME"))
            _, bareHeader = readZemax(bareNameFixture)
            @test bareHeader.name == ""
            rm(bareNameFixture)
        end

        @testset "XDAT parsing (TYPE XPOLYNOM support, TODO.md #4 / FIXED.md)" begin
            # XDAT (Zemax's "Extra Data" editor, used by TYPE XPOLYNOM)
            # is a separate keyword from PARM -- insert some XDAT lines
            # into SURF 2's block, out of order and with a gap, to
            # exercise the on-demand-growable extraData array.
            xdatFixture = tempname() * ".zmx"
            write(xdatFixture, replace(read(fixture, String),
                "  COAT TESTCOAT\n  DISZ 20.0" =>
                "  COAT TESTCOAT\n  XDAT 3 -0.01\n  XDAT 1 14.0\n  DISZ 20.0"))
            xdatSurfs, = readZemax(xdatFixture)
            @test xdatSurfs[3].extraData == [14.0, 0.0, -0.01]
            rm(xdatFixture)
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
            newbase, newdir, newydir, newrinOut, surf = OpticTrace.zemaxsurfToSurface(1, rinIn, rinOut, base, dir, nothing, s)

            @test newbase ≈ base + s.distance .* dir
            @test newdir == dir # unchanged: only TILTSURF/COORDBRK change dir
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
            newbase, newdir, newydir, newrinOut, surf = OpticTrace.zemaxsurfToSurface(2, rinIn, rinOut, base, dir, nothing, s)

            @test newbase ≈ base + s.distance .* dir
            @test surf.profile isa OpticTrace.SurfProfileEvenAsphere
            @test surf.profile.curv == s.curvature
            @test surf.profile.ϵ == OpticTrace.conicToϵ(s.conic)
            @test surf.profile.a == s.aspherics[2:end]
            @test surf.surfname == "A 2"
        end

        @testset "unsupported surface type" begin
            s = OpticTrace.ZemaxSurf(0.0, 1.0, "DEFAULT", 5.0, false, 0.0, zeros(Float64, 20), "", "TORIC", "")
            @test_throws ErrorException OpticTrace.zemaxsurfToSurface(1, rinIn, rinOut, base, dir, nothing, s)
        end

        @testset "TILTSURF" begin
            # PARM 1/2 = tilt about x/y (deg); no decenter, no z tilt, no
            # order flag -- confirmed against real TILTSURF samples in the
            # local Zemax install (see TODO.md #4 / FIXED.md).
            parm = zeros(Float64, 20)
            parm[1] = 30.0
            s = OpticTrace.ZemaxSurf(0.0, 5.0, "DEFAULT", 5.0, false, 0.0, parm, "", "TILTSURF", "T")
            newbase, newdir, newydir, newrinOut, surf = OpticTrace.zemaxsurfToSurface(1, rinIn, rinOut, base, dir, nothing, s)

            @test surf.profile isa OpticTrace.SurfProfileConic
            expectedbase, expecteddir, expectedydir = OpticTrace.zemaxCoordBreakFrame(base, dir, nothing,
                0.0, 0.0, parm[1], parm[2], 0.0, 0)
            @test newdir ≈ expecteddir
            @test newydir ≈ expectedydir
            @test surf.base.base ≈ base # TILTSURF tilts in place, no decenter
            @test newbase ≈ surf.base.base + s.distance * newdir
        end

        @testset "PARAXIAL (TODO.md #4 / FIXED.md)" begin
            parm = zeros(Float64, 20)
            parm[1] = 50.0
            s = OpticTrace.ZemaxSurf(0.0, 12.47, "DEFAULT", 5.0, false, 0.0, parm, "", "PARAXIAL", "")
            newbase, newdir, newydir, newrinOut, surf = OpticTrace.zemaxsurfToSurface(1, rinIn, rinOut, base, dir, nothing, s)

            @test surf.profile isa OpticTrace.ParaxialProfile
            @test surf.mod isa OpticTrace.ParaxialLensT
            @test surf.mod.focalLength == parm[1]
            @test surf.mod.refIndexIn == rinIn
            @test surf.mod.refIndexOut == rinOut
            @test surf.color == :cyan3 # distinct from the usual :lightgray
            @test newdir == dir # PARAXIAL doesn't tilt (unlike TILTSURF)
        end
    end

    @testset "zemaxsurfToProfileAperture" begin
        # the shared helper zemaxsurfToSurface/zemaxObjectToModelSurface both
        # build on -- see TODO.md Prerequisite 2 (OpticalSystem)
        @testset "STANDARD" begin
            s = OpticTrace.ZemaxSurf(0.03, 4.0, "TESTGLASS", 6.0, true, 0.2, zeros(Float64, 20), "coat1", "STANDARD", "S")
            profile, aperture = OpticTrace.zemaxsurfToProfileAperture(s)
            @test profile isa OpticTrace.SurfProfileConic
            @test profile.curv == s.curvature
            @test profile.ϵ == OpticTrace.conicToϵ(s.conic)
            @test aperture isa SizeLens
            @test aperture.semiDiameter == s.radius
        end

        @testset "EVENASPH" begin
            asph = zeros(Float64, 20)
            asph[2] = 0.0007
            s = OpticTrace.ZemaxSurf(-0.02, 8.0, "DEFAULT", 5.0, false, -0.6, asph, "coat2", "EVENASPH", "A")
            profile, aperture = OpticTrace.zemaxsurfToProfileAperture(s)
            @test profile isa OpticTrace.SurfProfileEvenAsphere
            @test profile.a == s.aspherics[2:end]
            @test aperture.semiDiameter == s.radius
        end

        @testset "TOROIDAL" begin
            # PARM 1 (aspherics[1]) is Zemax's "Radius of Rotation" Rx --
            # confirmed against real TOROIDAL samples in the local Zemax
            # install (see TODO.md #18 / FIXED.md).
            parm = zeros(Float64, 20)
            parm[1] = 1.0e5
            s = OpticTrace.ZemaxSurf(0.02, 10.0, "BK7", 11.0, false, -0.5, parm, "", "TOROIDAL", "T")
            profile, aperture = OpticTrace.zemaxsurfToProfileAperture(s)
            @test profile isa OpticTrace.SurfProfileToroid
            @test profile.curvY == s.curvature
            @test profile.ϵY == OpticTrace.conicToϵ(s.conic)
            @test profile.curvX ≈ 1 / parm[1]
            @test aperture.semiDiameter == s.radius
        end

        @testset "TOROIDAL, Rx == 0 (no x sweep, per Zemax's own convention)" begin
            s = OpticTrace.ZemaxSurf(-0.045, 40.0, "DEFAULT", 5.0, false, 0.0, zeros(Float64, 20), "", "TOROIDAL", "")
            profile, aperture = OpticTrace.zemaxsurfToProfileAperture(s)
            @test profile isa OpticTrace.SurfProfileToroid
            @test profile.curvX == 0.0
        end

        @testset "ODDASPHE" begin
            # Zemax's own Parameter i multiplies r^i directly for i=1..8
            # -- confirmed against a real axicon sample, where PARM 1
            # alone produces the linear cone an axicon is defined by
            # (see TODO.md #4 / FIXED.md).
            parm = zeros(Float64, 20)
            parm[1] = 0.01
            s = OpticTrace.ZemaxSurf(0.0, 2.0, "BK7", 2.0, false, 0.0, parm, "", "ODDASPHE", "")
            profile, aperture = OpticTrace.zemaxsurfToProfileAperture(s)
            @test profile isa OpticTrace.SurfProfileOddAsphere
            @test profile.curv == s.curvature
            @test profile.ϵ == OpticTrace.conicToϵ(s.conic)
            @test profile.a == parm[1:8]
            @test aperture.semiDiameter == s.radius
        end

        @testset "XPOLYNOM" begin
            # coefficients come from extraData (Zemax's XDAT lines), not
            # aspherics (PARM) -- extraData[1] is the normalization
            # radius, extraData[2] an unidentified control flag not
            # stored, extraData[3:end] the polynomial term coefficients
            # (see TODO.md #4 / FIXED.md).
            extraData = zeros(Float64, 7)
            extraData[1] = 14.0
            extraData[2] = 1.0
            extraData[3] = -0.01 # x term
            extraData[5] = -0.01 # x^2 term
            s = OpticTrace.ZemaxSurf(0.0, -25.0, "MIRROR", 11.8, false, 0.0, zeros(Float64, 20), "", "XPOLYNOM", "",
                extraData)
            profile, aperture = OpticTrace.zemaxsurfToProfileAperture(s)
            @test profile isa OpticTrace.SurfProfileXYPoly
            @test profile.curv == s.curvature
            @test profile.ϵ == OpticTrace.conicToϵ(s.conic)
            @test profile.normRadius == extraData[1]
            @test profile.a == extraData[3:end]
            @test aperture.semiDiameter == s.radius
        end

        @testset "PARAXIAL" begin
            # an ideal thin lens has no real sag -- PARM 1 (focal
            # length) is handled by zemaxsurfToSurface's bend
            # selection, not here (see TODO.md #4 / FIXED.md).
            parm = zeros(Float64, 20)
            parm[1] = 50.0
            s = OpticTrace.ZemaxSurf(0.0, 12.47, "DEFAULT", 5.0, false, 0.0, parm, "", "PARAXIAL", "")
            profile, aperture = OpticTrace.zemaxsurfToProfileAperture(s)
            @test profile isa OpticTrace.ParaxialProfile
            @test aperture.semiDiameter == s.radius
        end

        @testset "unsupported surface type" begin
            s = OpticTrace.ZemaxSurf(0.0, 1.0, "DEFAULT", 5.0, false, 0.0, zeros(Float64, 20), "", "TORIC", "")
            @test_throws ErrorException OpticTrace.zemaxsurfToProfileAperture(s)
        end
    end

    @testset "zemaxObjectToModelSurface" begin
        @testset "finite object distance" begin
            s = OpticTrace.ZemaxSurf(0.01, 12.0, "DEFAULT", 3.0, false, 0.0, zeros(Float64, 20), "", "STANDARD", "Obj")
            surf = OpticTrace.zemaxObjectToModelSurface(s, ORIGIN, ZAXIS)
            @test surf isa ModelSurface
            @test surf.base.base ≈ ORIGIN - 12.0 .* ZAXIS
            @test surf.profile isa OpticTrace.SurfProfileConic
            @test surf.profile.curv == s.curvature
            @test surf.aperture.semiDiameter == s.radius
            @test surf.surfname == "Obj"
        end

        @testset "infinite object distance -- positioned at basept" begin
            s = OpticTrace.ZemaxSurf(0.0, Inf, "DEFAULT", 0.0, false, 0.0, zeros(Float64, 20), "", "STANDARD", "")
            surf = OpticTrace.zemaxObjectToModelSurface(s, ORIGIN, ZAXIS)
            @test surf.base.base == ORIGIN
            @test surf.surfname == "Object" # falls back when s.comm is empty
        end

        @testset "custom refractive index" begin
            s = OpticTrace.ZemaxSurf(0.0, 1.0, "DEFAULT", 0.0, false, 0.0, zeros(Float64, 20), "", "STANDARD", "")
            surf = OpticTrace.zemaxObjectToModelSurface(s, ORIGIN, ZAXIS; rin = 1.33)
            @test surf.refIndex == 1.33
        end
    end

    @testset "zemaxCoordBreakFrame" begin
        @testset "no-op (all zeros) leaves the frame unchanged" begin
            newbase, newdir, newydir = OpticTrace.zemaxCoordBreakFrame(ORIGIN, ZAXIS, nothing,
                0.0, 0.0, 0.0, 0.0, 0.0, 0)
            @test newbase == ORIGIN
            @test newdir ≈ ZAXIS
            @test newydir ≈ YAXIS # the same default findPerpenMap(ZAXIS, nothing) guess
        end

        @testset "45 deg tilt about X folds an on-axis beam by 90 deg (Fold Mirror Using Coordinate Breaks.ZMX)" begin
            # Empirically validated against the local Zemax install's real
            # sample: a 45 deg COORDBRK tilt about local x immediately
            # ahead of a mirror should fold a beam travelling along +Z by
            # exactly 90 deg -- confirms the right-hand-rule sign
            # convention used here matches Zemax's own.
            _, newdir, newydir = OpticTrace.zemaxCoordBreakFrame(ORIGIN, ZAXIS, nothing,
                0.0, 0.0, 45.0, 0.0, 0.0, 0)
            incoming = ZAXIS
            reflected = incoming - 2 * dot(incoming, newdir) * newdir
            @test reflected ≈ YAXIS atol=1e-12
        end

        @testset "tilt about Z only: dir unchanged, ydir rotates in the x-y plane (Toroid.zmx)" begin
            _, newdir, newydir = OpticTrace.zemaxCoordBreakFrame(ORIGIN, ZAXIS, nothing,
                0.0, 0.0, 0.0, 0.0, 45.0, 0)
            @test newdir ≈ ZAXIS
            @test newydir ≈ normalize(Vec3(-sind(45.0), cosd(45.0), 0.0))
        end

        @testset "decenter only (order flag 0): moves basept along the current x/y axes" begin
            newbase, newdir, newydir = OpticTrace.zemaxCoordBreakFrame(ORIGIN, ZAXIS, nothing,
                1.0, 2.0, 0.0, 0.0, 0.0, 0)
            @test newbase ≈ ORIGIN + 1.0 * XAXIS + 2.0 * YAXIS
            @test newdir ≈ ZAXIS
        end

        @testset "order flag: decenter-then-tilt (0) vs tilt-then-decenter (1) give different basepts" begin
            base0, dir0, _ = OpticTrace.zemaxCoordBreakFrame(ORIGIN, ZAXIS, nothing, 1.0, 0.0, 0.0, 0.0, 90.0, 0)
            base1, dir1, _ = OpticTrace.zemaxCoordBreakFrame(ORIGIN, ZAXIS, nothing, 1.0, 0.0, 0.0, 0.0, 90.0, 1)
            @test dir0 ≈ dir1 ≈ ZAXIS # a pure Z tilt doesn't change dir either way
            @test base0 ≈ ORIGIN + 1.0 * XAXIS # decenter applied in the pre-tilt (original) x axis
            @test base1 ≈ ORIGIN + 1.0 * normalize(Vec3(cosd(90.0), sind(90.0), 0.0)) # decenter applied in the post-tilt x axis
            @test !isapprox(base0, base1) # the two orders genuinely disagree here
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

        @testset "isinf guard (TODO.md Prerequisite 2): rejects an infinite distance on a non-object surface" begin
            # zsurfs itself (unsliced) still has the object surface's DISZ
            # INFINITY at index 1 -- zemaxsurfsToGeo must reject it rather
            # than silently propagating Inf/NaN into later surfaces
            @test_throws ErrorException zemaxsurfsToGeo(zsurfs, ORIGIN, ZAXIS, 0.5875618; glassCatalog = testCatalog)
        end

        @testset "COORDBRK (Phase 2): tilts/decenters the frame, builds no surface" begin
            coordbrkParm = zeros(Float64, 20)
            coordbrkParm[1] = 0.5   # decenter X
            coordbrkParm[3] = 45.0  # tilt X (deg)
            coordbrk = OpticTrace.ZemaxSurf(0.0, 0.0, "DEFAULT", 0.0, false, 0.0, coordbrkParm, "", "COORDBRK", "")
            standard = OpticTrace.ZemaxSurf(0.02, 10.0, "TESTGLASS", 5.0, false, 0.0, zeros(Float64, 20), "", "STANDARD", "S")

            geo = zemaxsurfsToGeo([coordbrk, standard], ORIGIN, ZAXIS, 0.5875618; glassCatalog = testCatalog)

            @test length(geo) == 1 # the COORDBRK itself builds no surface

            expectedbase, expecteddir, expectedydir = OpticTrace.zemaxCoordBreakFrame(ORIGIN, ZAXIS, nothing,
                coordbrkParm[1], coordbrkParm[2], coordbrkParm[3], coordbrkParm[4], coordbrkParm[5], coordbrkParm[6])
            # COORDBRK's own DISZ (0 here) is applied last, along the new dir
            @test geo[1].base.base ≈ expectedbase + coordbrk.distance * expecteddir
            @test geo[1].base.dir ≈ expecteddir
            @test geo[1].base.ydir ≈ expectedydir
            @test geo[1].mod.refIndexIn == refIndexDefault # unaffected by the coordbreak
            @test geo[1].mod.refIndexOut == 1.6
        end

        @testset "GLAS MIRROR (TODO.md #27): MirrorR bend, glass catalog not consulted" begin
            standard = OpticTrace.ZemaxSurf(0.0, 5.0, "TESTGLASS", 5.0, false, 0.0, zeros(Float64, 20), "", "STANDARD", "S")
            mirror = OpticTrace.ZemaxSurf(0.0, 10.0, "MIRROR", 5.0, false, 0.0, zeros(Float64, 20), "", "STANDARD", "M")
            # testCatalog has no "MIRROR" key -- this must not throw a KeyError
            geo = zemaxsurfsToGeo([standard, mirror], ORIGIN, ZAXIS, 0.5875618; glassCatalog = testCatalog)

            @test length(geo) == 2
            @test geo[2].mod isa OpticTrace.MirrorR
            @test geo[2].mod.refIndexIn == 1.6  # TESTGLASS -- unchanged by the mirror
            @test geo[2].mod.refIndexOut == 1.6
        end
    end

    @testset "readZemaxSystem / OpticalSystem" begin
        testCatalog = Dict{AbstractString,Any}(
            "DEFAULT" => (λ -> 1.0),
            "TESTGLASS" => (λ -> 1.6),
        )

        sys = readZemaxSystem(fixture; wavelength = 0.5875618, glassCatalog = testCatalog)

        @test sys isa OpticalSystem
        @test length(sys.geo) == 3
        @test sys.geo[1].base.base == ORIGIN # object surface excluded, chain still starts at basept

        @testset "object surface preserved separately, never in geo" begin
            @test sys.objectSurface isa ModelSurface
            @test sys.objectDistance == Inf
            @test sys.objectAtInfinity == true
            @test sys.objectSurface.base.base == ORIGIN # infinite-distance placeholder anchor
            @test all(s -> !(s isa ModelSurface), sys.geo)
        end

        @testset "header metadata absorbed" begin
            @test sys.name == "Test Singlet"
            @test sys.units == "MM"
            @test sys.primaryWavelengthIndex == 1
            @test sys.apertureType == "ENPD"
            @test sys.apertureValue == 8.0
            @test sys.glassCatalogs == ["SCHOTT", "MISC"]
            @test sys.mode == "SEQ"
            @test sys.fields == [Point2(0.0, 0.0), Point2(5.0, 3.0)]
            @test sys.fieldWeight == [1.0, 0.5]
            @test occursin("Test fixture", sys.notes)
        end

        @testset "MODE NSC rejected" begin
            nscFixture = tempname() * ".zmx"
            write(nscFixture, replace(read(fixture, String), "MODE SEQ" => "MODE NSC"))
            @test_throws ErrorException readZemaxSystem(nscFixture; glassCatalog = testCatalog)
            rm(nscFixture)
        end
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
