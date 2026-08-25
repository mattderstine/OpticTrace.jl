@testset "foundations.jl" begin

    @testset "Constants & globals" begin
        @test ∞ === Inf
        @test ORIGIN == Point3(0., 0., 0.)
        @test ZAXIS == Vec3(0., 0., 1.)
        @test YAXIS == Vec3(0., 1., 0.)
        @test XAXIS == Vec3(1., 0., 0.)
        @test OpticTrace.identityPol == [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

        @test rInDef() == 1.0
        @test rInDef(0.5) == rInDef() # rInDef(wave) ignores its argument

        old = rInDef()
        try
            setRInDef(1.33)
            @test rInDef() == 1.33
            @test rInDef(0.5) == 1.33
        finally
            setRInDef(old)
        end
        @test rInDef() == old
    end

    @testset "conicToϵ" begin
        @test OpticTrace.conicToϵ(0.0) == 1.0
        @test OpticTrace.conicToϵ(-1.0) == 0.0
    end

    @testset "Core struct construction" begin
        r = Ray(Point3(1., 2., 3.), Vec3(0., 0., 1.))
        @test r.base == Point3(1., 2., 3.)
        @test r.dir == Vec3(0., 0., 1.)
        @test r isa Ray{3,Float64}

        sb = SurfBase(Point3(0., 0., 0.), ZAXIS, YAXIS)
        @test sb.base == ORIGIN
        @test sb.dir == ZAXIS
        @test sb.ydir == YAXIS

        sl = SizeLens(5.0)
        @test sl.semiDiameter == 5.0

        ra = RoundAperture(1.0, 5.0)
        @test ra.obscure == 1.0
        @test ra.semiDiameter == 5.0

        rect = RectAperture(1.0, 2.0, 5.0, 6.0)
        @test rect.wo == 1.0
        @test rect.lo == 2.0
        @test rect.wclear == 5.0
        @test rect.lclear == 6.0

        amp = AmpData(OpticTrace.identityPol, OpticTrace.identityPol, [1.0])
        @test amp.p == OpticTrace.identityPol
        @test amp.o == OpticTrace.identityPol
        @test amp.trans == [1.0]

        trc = Trace(r, 1.0, 0.0, amp)
        @test trc.ray == r
        @test trc.nIn == 1.0
        @test trc.delta == 0.0
        @test trc.pmatrix === amp

        dT = DielectricT(1.0, 1.5)
        @test dT.refIndexIn == 1.0
        @test dT.refIndexOut == 1.5

        mR = MirrorR(1.0, 1.0)
        @test mR.refIndexIn == 1.0
        @test mR.refIndexOut == 1.0

        cd = CDiffuser(tan(0.1), 1.0, 1.5)
        @test cd.tanθ ≈ tan(0.1)
        @test cd.refIndexIn == 1.0
        @test cd.refIndexOut == 1.5

        nb = NoBendIndex(1.0, 1.5)
        @test nb.refIndexIn == 1.0
        @test nb.refIndexOut == 1.5

        nb2 = NoBendIndex(1.33) # convenience constructor: same index both sides
        @test nb2.refIndexIn == 1.33
        @test nb2.refIndexOut == 1.33

        ap = OpticTrace.AmpParam("test")
        @test ap.type == "test"

        nap = NoAmpParam("none")
        @test nap.type == "none"
    end

    @testset "Coordinate transforms" begin

        @testset "findPerpenMap" begin
            # default ydir guess
            newy, toGlobalDir = findPerpenMap(ZAXIS)
            newx = toGlobalDir(Vec3(1., 0., 0.))
            @test toGlobalDir(Vec3(0., 1., 0.)) ≈ newy
            @test toGlobalDir(Vec3(0., 0., 1.)) ≈ ZAXIS
            @test norm(newx) ≈ 1.0
            @test norm(newy) ≈ 1.0
            @test dot(newx, newy) ≈ 0.0 atol = 1e-12
            @test dot(newx, ZAXIS) ≈ 0.0 atol = 1e-12
            @test dot(newy, ZAXIS) ≈ 0.0 atol = 1e-12

            # explicit ydir
            newy2, toGlobalDir2 = findPerpenMap(ZAXIS, YAXIS)
            @test newy2 == YAXIS
            newx2 = toGlobalDir2(Vec3(1., 0., 0.))
            @test newx2 ≈ XAXIS
            @test toGlobalDir2(Vec3(0., 0., 1.)) ≈ ZAXIS

            # ydir=nothing matches the default-guess (2-arg) method exactly
            newy3, toGlobalDir3 = findPerpenMap(ZAXIS, nothing)
            @test newy3 == newy
            @test toGlobalDir3(Vec3(1., 0., 0.)) == newx
        end

        @testset "updateCoordChange" begin
            pip = Point3(1., 2., 3.)
            ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir = updateCoordChange(pip, ZAXIS)
            @test toGlobalCoord(ORIGIN) ≈ pip
            @test toLocalCoord(pip) ≈ ORIGIN

            # round-trip two sample points through both coordinate and direction maps
            for testpt in (Point3(0.5, -0.3, 2.0), Point3(-1.2, 4.0, 0.1))
                @test toLocalCoord(toGlobalCoord(testpt)) ≈ testpt
                @test toGlobalCoord(toLocalCoord(testpt)) ≈ testpt
            end
            for testdir in (Vec3(0.3, 0.4, 0.5), Vec3(-0.6, 0.1, 0.2))
                @test toLocalDir(toGlobalDir(testdir)) ≈ testdir
                @test toGlobalDir(toLocalDir(testdir)) ≈ testdir
            end

            ydir2, = updateCoordChange(pip, ZAXIS, YAXIS)
            @test ydir2 == YAXIS
        end

        @testset "updateCoordChange!" begin
            surf = referencePlane("test", ORIGIN, ZAXIS, 1.0, 5.0, "none")

            # 1-arg method: recompute in place from the surface's current base fields
            surf.base.base = Point3(1., 2., 3.)
            ydir0, toGlobalCoord0, toLocalCoord0, = updateCoordChange(surf.base.base, surf.base.dir, surf.base.ydir)
            updateCoordChange!(surf)
            @test surf.base.ydir == ydir0
            @test surf.toGlobalCoord(ORIGIN) ≈ toGlobalCoord0(ORIGIN)
            @test surf.toLocalCoord(surf.base.base) ≈ ORIGIN

            # 3-arg method: replace surf.base wholesale, then recompute
            newbase = SurfBase(Point3(4., 5., 6.), YAXIS, XAXIS)
            updateCoordChange!(surf, newbase)
            @test surf.base.base == Point3(4., 5., 6.)
            @test surf.base.dir == YAXIS
            @test surf.toLocalCoord(Point3(4., 5., 6.)) ≈ ORIGIN

            updateCoordChange!(surf, newbase, XAXIS)
            @test surf.base.ydir == XAXIS
        end

        @testset "defaultSetupGeo / updateEGeo!" begin
            testfunc(parameters, wl) = (parameters[:scale] * wl, parameters[:tag])
            object = referencePlane("object", ORIGIN, ZAXIS, 1.0, 5.0, "none")
            egeo = ExtendedGeometry(AbstractSurface[], testfunc, OpticTrace.defaultSetupGeo, object,
                [0.5, 0.6], Dict(:scale => 2.0, :tag => "hi"))

            result1 = OpticTrace.defaultSetupGeo(testfunc, object, egeo.wavelength, egeo.parameters)
            @test result1 == (1.0, "hi") # default numWL=1 -> wavelength[1]=0.5, scale*0.5=1.0

            result2 = OpticTrace.defaultSetupGeo(testfunc, object, egeo.wavelength, egeo.parameters; numWL=2)
            @test result2 == (1.2, "hi") # wavelength[2]=0.6, scale*0.6=1.2

            # updateEGeo! mutates egeo.geo, so it needs a funcGeo that
            # actually returns an Array{AbstractSurface}, matching the
            # egeo.geo field's type -- testfunc above returns a plain Tuple,
            # fine for exercising defaultSetupGeo directly (result1/result2),
            # but not assignable into egeo.geo, so it's not reused here.
            geoFunc(parameters, wl) = AbstractSurface[referencePlane("dyn", ORIGIN, ZAXIS, 1.0, parameters[:scale] * wl, "none")]
            egeo2 = ExtendedGeometry(AbstractSurface[], geoFunc, OpticTrace.defaultSetupGeo, object,
                [0.5, 0.6], Dict(:scale => 2.0, :tag => "hi"))

            ret2 = OpticTrace.updateEGeo!(egeo2)
            @test egeo2.geo === ret2
            @test egeo2.geo[1].aperture.semiDiameter == 1.0 # scale*wavelength[1] = 2.0*0.5
        end
    end

    @testset "Aperture" begin

        @testset "isAperture" begin
            @test isAperture(SizeLens(5.0)) == false
            @test isAperture(RoundAperture(1.0, 5.0)) == true
            @test isAperture(RectAperture(1.0, 1.0, 5.0, 5.0)) == true
        end

        @testset "clipAperture" begin
            ra = RoundAperture(1.0, 5.0) # obscure=1.0, semiDiameter=5.0
            @test OpticTrace.clipAperture(Point3(0.0, 0.0, 0.0), ra) == true  # inside obscuration
            @test OpticTrace.clipAperture(Point3(3.0, 0.0, 0.0), ra) == false # in clear region
            @test OpticTrace.clipAperture(Point3(6.0, 0.0, 0.0), ra) == true  # beyond semiDiameter

            rect = RectAperture(1.0, 1.0, 5.0, 5.0) # wo=lo=1.0, wclear=lclear=5.0
            @test OpticTrace.clipAperture(Point3(0.5, 0.5, 0.0), rect) == true  # inside obscuration
            @test OpticTrace.clipAperture(Point3(3.0, 3.0, 0.0), rect) == false # in clear region
            @test OpticTrace.clipAperture(Point3(6.0, 0.0, 0.0), rect) == true  # beyond wclear
            @test OpticTrace.clipAperture(Point3(0.0, 6.0, 0.0), rect) == true  # beyond lclear
        end

        @testset "roundAperture / rectAperture" begin
            surf = roundAperture("stop", ORIGIN, ZAXIS, 1.0, 1.0, 5.0)
            @test surf isa ModelSurface
            @test surf.aperture isa RoundAperture
            @test surf.aperture.obscure == 1.0
            @test surf.aperture.semiDiameter == 5.0
            @test surf.refIndex == 1.0
            @test surf.profile isa NoProfile

            surf2 = rectAperture("stop2", ORIGIN, ZAXIS, YAXIS, 1.0, 0.5, 0.5, 5.0, 5.0)
            @test surf2.aperture isa RectAperture
            @test surf2.aperture.wo == 0.5
            @test surf2.aperture.lo == 0.5
            @test surf2.aperture.wclear == 5.0
            @test surf2.aperture.lclear == 5.0
        end
    end

end
