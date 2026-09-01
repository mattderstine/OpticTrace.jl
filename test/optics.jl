


@testset "optics.jl" begin
    # Write your tests here.

    # test objects

    asphere_coeff1 = 0.1
    asphere_coeff2 = 0.01
    factor1 = 4.0
    factor2 = 8.0
    curve1 = 1.0
    curve2 = sqrt(0.5)
    curve3 = 0.0
    curve = curve2
    spc1 = OpticTrace.SurfProfileConic(curve1, 0.0)
    spc2 = OpticTrace.SurfProfileConic(curve2, 1.0)
    spc3 = OpticTrace.SurfProfileConic(curve3, 1.0)
    sps = OpticTrace.SurfProfileSphere(curve)
    spsAsphere = OpticTrace.SurfProfileAsphere(0.0, 0.0, [0.00,  asphere_coeff1,0.0, asphere_coeff2])
    spsEAsphere = OpticTrace.SurfProfileEvenAsphere(0.0, 0.0, [asphere_coeff1, asphere_coeff2])

    offsetOA = Vec3(0.3, -0.2, 0.05)
    spOA = OpticTrace.SurfProfileOAConic(curve1, 0.0, offsetOA)
    baseConicOA = OpticTrace.SurfProfileConic(curve1, 0.0)

    spsCyl1 = OpticTrace.SurfProfileCyl(curve1, 0.0, Float64[])
    spsCyl2 = OpticTrace.SurfProfileCyl(curve2, 1.0, Float64[])

    # axicon-style: a[1] (r^1 term) dominates, unlike SurfProfileAsphere
    # whose a[1] is r^3 -- see TODO.md #4 / FIXED.md
    oddAsphereCoeff1 = 0.01
    spsOddAsphere = OpticTrace.SurfProfileOddAsphere(0.0, 0.0, [oddAsphereCoeff1, 0.0])

    # reflector5.ZMX-style: x^2/y^2 terms (Zemax term indices 3 and 5)
    xyPolyCoeffs = zeros(5)
    xyPolyCoeffs[3] = -0.01 # x^2
    xyPolyCoeffs[5] = -0.01 # y^2
    spsXYPoly = OpticTrace.SurfProfileXYPoly(0.0, 0.0, 14.0, xyPolyCoeffs)

    curvYToroid = 0.3
    curvXToroid = 0.2
    # ϵY=1.0 (true circular y-z cross-section) so the x=0/y=0 cross-section
    # checks below can compare directly against SurfProfileSphere, which
    # always uses the ϵ=1 (spherical) formula, not the ϵ=0 paraxial one.
    spsToroid = OpticTrace.SurfProfileToroid(curvYToroid, 1.0, curvXToroid)

    # `simplesystem` is unused elsewhere in this file, but building it
    # calls lens_TLAC254_060, which looks up real glass files from
    # OpticTrace.dirBaseRefractiveIndex -- a machine-local directory,
    # not part of the repo (see TODO.md). Guarded so this file still
    # runs (skipping only this dead assignment) on a checkout without
    # that directory, e.g. CI.
    if HAS_GLASS_CATALOG
        simplesystem = [[referencePlane("object", ORIGIN, ZAXIS, 1.0, 10.0, "test coating")]
            lens_TLAC254_060(Point(0.0, 0.0, 10.0), ZAXIS, 0.45; order="forward", lensname="TL_AC254-060");
            [referencePlane("image", Point(0.0, 0.0, 70.0), ZAXIS, 1.0, 10.0, "test coating")]
        ]
    end

    #sag tests
    @testset "sag tests" begin

        sag_value = sag(1.0, 1.0, spc1)
        @test sag_value ≈ 1.0

        sag_value = sag(0.9999999999999999, 1.0, spc2)
        @test sag_value ≈ sqrt(2)

        sag_value = sag(1.0, 1.0, spc3)
        @test sag_value == 0.0

        sag_value = sag(0.9999999999999999, 1.0, sps)
        @test sag_value ≈ sqrt(2)

        sag_value = sag(100.0, 1.0, sps)
        @test isnan(sag_value)

        sag_value = sag(100.0, 1.0, spc2)
        @test isnan(sag_value)



        #r = sqrt(2) r2=2, r4 = 4, r6 = 8

        sag_value = sag(1.0, 1.0, spsAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2

        sag_value = sag(1.0, 1.0, spsEAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2
        sag_value = sag(1.0, 1.0, spsAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2

        sag_value = sag(1.0, 1.0, spsEAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2

        # SurfProfileCyl: x does not appear in the formula (only y does)
        sag_value = sag(1.0, 1.0, spsCyl1) # ϵ=0 branch: z = curv*y^2/2
        @test sag_value ≈ curve1 * 1.0^2 * 0.5
        @test sag(5.0, 1.0, spsCyl1) == sag_value # x is ignored

        sag_value = sag(1.0, 1.0, spsCyl2) # ϵ=1 branch (circular cross-section)
        @test sag_value ≈ sag(0.0, 1.0, sps) # reduces to SurfProfileSphere's formula along y
        @test sag(-3.0, 1.0, spsCyl2) == sag_value # x is ignored

        # SurfProfileToroid (fixed, see TODO.md #18 / FIXED.md): a toroid's
        # x=0 cross-section is exactly the base y-z curve's own sag, which
        # for ϵY=0 matches SurfProfileSphere(curvY).
        expectedToroidSag = sag(0.0, 0.5, OpticTrace.SurfProfileSphere(curvYToroid))
        actualToroidSag = sag(0.0, 0.5, spsToroid)
        @test actualToroidSag ≈ expectedToroidSag

        # y=0 cross-section is purely the x-sweep term, curvature curvX
        @test sag(0.5, 0.0, spsToroid) ≈ sag(0.5, 0.0, OpticTrace.SurfProfileSphere(curvXToroid))

        # curvX == 0 disables the x sweep entirely (pure extrusion along x)
        spsToroidNoXSweep = OpticTrace.SurfProfileToroid(curvYToroid, 0.0, 0.0)
        @test sag(3.0, 0.5, spsToroidNoXSweep) == sag(0.0, 0.5, spsToroidNoXSweep)

        # SurfProfileOddAsphere (TODO.md #4 / FIXED.md): a[1] is r^1, so
        # sag is exactly linear in r along any ray through the origin --
        # this is what makes it able to represent an axicon, unlike
        # SurfProfileAsphere (a[1] is r^3).
        @test sag(1.0, 0.0, spsOddAsphere) ≈ oddAsphereCoeff1
        @test sag(2.0, 0.0, spsOddAsphere) ≈ 2 * oddAsphereCoeff1
        @test sag(0.0, 2.0, spsOddAsphere) ≈ 2 * oddAsphereCoeff1 # radially symmetric

        # SurfProfileXYPoly (TODO.md #4 / FIXED.md): x^2/y^2 terms in
        # *normalized* (x/normRadius, y/normRadius) coordinates.
        xyPolyR = 5.0
        expectedXYPolySag = xyPolyCoeffs[3] * (xyPolyR / spsXYPoly.normRadius)^2
        @test sag(xyPolyR, 0.0, spsXYPoly) ≈ expectedXYPolySag
        @test sag(0.0, xyPolyR, spsXYPoly) ≈ expectedXYPolySag # x^2/y^2 coefficients are equal here
        @test sag(xyPolyR, xyPolyR, spsXYPoly) ≈ 2 * expectedXYPolySag

        @testset "xyPolyTermPowers" begin
            @test OpticTrace.xyPolyTermPowers(1) == (1, 0) # x
            @test OpticTrace.xyPolyTermPowers(2) == (0, 1) # y
            @test OpticTrace.xyPolyTermPowers(3) == (2, 0) # x^2
            @test OpticTrace.xyPolyTermPowers(4) == (1, 1) # xy
            @test OpticTrace.xyPolyTermPowers(5) == (0, 2) # y^2
            @test OpticTrace.xyPolyTermPowers(6) == (3, 0) # x^3
            @test OpticTrace.xyPolyTermPowers(9) == (0, 3) # y^3
            @test OpticTrace.xyPolyTermPowers(10) == (4, 0) # x^4
        end
    end

    @testset "surfNormal & deltaToSurf tests" begin



        offset = 0.2
        r2 = offset^2
        ray = Ray(Point(0.0, offset, -1.0), Vec(0.0, 0.0, 1.0))

        #sphere test
        sag_value = curve * r2 / (1.0 + sqrt(1 - curve^2 * r2))
        delta = OpticTrace.deltaToSurf(ray, sps)
        @test delta ≈ 1.0+sag_value
        point = rprop(ray, delta)
        @test point ≈ Point3(0.0, offset, sag_value)
        normal = OpticTrace.surfNormal(point, sps)
        @test normal_from_sag(point, sps) ≈ normal

        #paraboloid test
        sag_value = curve1 * offset^2 * 0.5
        delta = OpticTrace.deltaToSurf(ray, spc1)
        @test delta ≈ 1.0+sag_value
        point = rprop(ray, delta)
        @test point ≈ Point3(0.0, offset, sag_value)
        normal = OpticTrace.surfNormal(point, spc1)
        @test normal_from_sag(point, spc1) ≈ normal

        #asphere test
        sag_value = offset^4 * asphere_coeff1 + offset^6 * asphere_coeff2
        delta = OpticTrace.deltaToSurf(ray, spsAsphere)
        @test delta ≈ 1.0+sag_value
        point = rprop(ray, delta)
        @test point ≈ Point3(0.0, offset, sag_value)
        normal = OpticTrace.surfNormal(point, spsAsphere)
        @test normal_from_sag(point, spsAsphere) ≈ normalize(normal)

        #even asphere test
        sag_value = offset^4 * asphere_coeff1 + offset^6 * asphere_coeff2
        delta = OpticTrace.deltaToSurf(ray, spsEAsphere)
        @test delta ≈ 1.0+sag_value
        point = rprop(ray, delta)
        @test point ≈ Point3(0.0, offset, sag_value)
        normal = OpticTrace.surfNormal(point, spsEAsphere)
        @test normal_from_sag(point, spsEAsphere) ≈ normalize(normal)

        #off-axis conic test (sag & surfNormal correct; deltaToSurf has a
        #known bug -- author's own comment flags "logic is flawed in this
        #one", see TODO.md)
        x0, y0 = 0.1, 0.15
        z0 = sag(x0, y0, spOA)
        @test z0 ≈ sag(x0 - offsetOA[1], y0 - offsetOA[2], baseConicOA) + offsetOA[3]

        normalOA = OpticTrace.surfNormal(Point3(x0, y0, z0), spOA)
        expectedNormalOA = OpticTrace.surfNormal(Point3(x0, y0, z0) .- offsetOA, baseConicOA)
        @test normalOA ≈ expectedNormalOA
        @test normal_from_sag(x0, y0, spOA) ≈ normalize(normalOA)

        oaRay = Ray(Point3(x0, y0, -1.0), Vec3(0.0, 0.0, 1.0))
        oaDelta = OpticTrace.deltaToSurf(oaRay, spOA)
        oaPoint = rprop(oaRay, oaDelta)
        @test_broken oaPoint ≈ Point3(x0, y0, z0)

        #cylinder test
        xC, yC = 0.3, 0.4
        zC = sag(xC, yC, spsCyl2)
        pointC = Point3(xC, yC, zC)
        normalC = OpticTrace.surfNormal(pointC, spsCyl2)
        @test normal_from_sag(xC, yC, spsCyl2) ≈ normalize(normalC)

        # deltaToSurf, ϵ=1 branch: reduces to the same circular
        # cross-section formula as the sphere test above, since
        # curve2 == curve and ray's x is 0
        sag_valueCyl = curve2 * r2 / (1.0 + sqrt(1 - curve2^2 * r2))
        deltaCyl = OpticTrace.deltaToSurf(ray, spsCyl2)
        @test deltaCyl ≈ 1.0 + sag_valueCyl
        pointCyl = rprop(ray, deltaCyl)
        @test pointCyl ≈ Point3(0.0, offset, sag_valueCyl)
        normalCyl = OpticTrace.surfNormal(pointCyl, spsCyl2)
        @test normal_from_sag(pointCyl, spsCyl2) ≈ normalize(normalCyl)

        #toroid test (TODO.md #18 / FIXED.md)
        xT, yT = 0.3, 0.4
        zT = sag(xT, yT, spsToroid)
        pointT = Point3(xT, yT, zT)
        normalT = OpticTrace.surfNormal(pointT, spsToroid)
        @test normal_from_sag(xT, yT, spsToroid) ≈ normalize(normalT)

        toroidRay = Ray(Point3(0.0, offset, -1.0), Vec3(0.0, 0.0, 1.0))
        deltaToroid = OpticTrace.deltaToSurf(toroidRay, spsToroid)
        pointToroid = rprop(toroidRay, deltaToroid)
        @test pointToroid ≈ Point3(0.0, offset, sag(0.0, offset, spsToroid))

        #odd asphere test (TODO.md #4 / FIXED.md) -- deltaToSurf and
        #surfNormal both come from AbstractAsphericProfile's generic
        #fallbacks (numeric root-find, ForwardDiff gradient), not
        #hand-derived closed forms
        oddAsphereRay = Ray(Point3(0.0, offset, -1.0), Vec3(0.0, 0.0, 1.0))
        deltaOddAsphere = OpticTrace.deltaToSurf(oddAsphereRay, spsOddAsphere)
        pointOddAsphere = rprop(oddAsphereRay, deltaOddAsphere)
        @test pointOddAsphere ≈ Point3(0.0, offset, sag(0.0, offset, spsOddAsphere))
        normalOddAsphere = OpticTrace.surfNormal(pointOddAsphere, spsOddAsphere)
        @test normal_from_sag(0.0, offset, spsOddAsphere) ≈ normalize(normalOddAsphere)

        #xy-polynomial test (TODO.md #4 / FIXED.md) -- same free
        #deltaToSurf/surfNormal fallbacks as odd asphere above
        xT2, yT2 = 3.0, 1.5
        zT2 = sag(xT2, yT2, spsXYPoly)
        pointXYPoly = Point3(xT2, yT2, zT2)
        normalXYPoly = OpticTrace.surfNormal(pointXYPoly, spsXYPoly)
        @test normal_from_sag(xT2, yT2, spsXYPoly) ≈ normalize(normalXYPoly)

        xyPolyRay = Ray(Point3(xT2, yT2, -1.0), Vec3(0.0, 0.0, 1.0))
        deltaXYPoly = OpticTrace.deltaToSurf(xyPolyRay, spsXYPoly)
        pointXYPolyRay = rprop(xyPolyRay, deltaXYPoly)
        @test pointXYPolyRay ≈ pointXYPoly
    end

    @testset "ParaxialProfile / ParaxialLensT (TODO.md #4 / FIXED.md)" begin
        pp = OpticTrace.ParaxialProfile(0.0)

        @test sag(3.0, -2.0, pp) == 0.0 # flat, like NoProfile

        paraxialRay = Ray(Point3(0.0, 2.0, -1.0), Vec3(0.0, 0.0, 1.0))
        paraxialDelta = OpticTrace.deltaToSurf(paraxialRay, pp)
        @test paraxialDelta ≈ 1.0 # same flat z=0 plane as NoProfile

        # surfNormal is deliberately NOT a true normal here -- it's the
        # raw local coordinates, unnormalized, fed to modFunc(::ParaxialLensT)
        @test OpticTrace.surfNormal(Point3(1.5, -0.7, 0.0), pp) == Vec3(1.5, -0.7, 0.0)

        @test OpticTrace.reverseProfile!(pp) === pp # no-op

        f = 50.0
        bend = OpticTrace.ParaxialLensT(f, 1.0, 1.5)

        @testset "classic thin-lens regression: axis-parallel ray converges to the focal point" begin
            h = 2.0
            offsetVec = Vec3(0.0, h, 0.0) # what surfNormal(::ParaxialProfile) would hand modFunc
            status, newDir, nIn = OpticTrace.modFunc(Ray(Point3(0.0, h, 0.0), ZAXIS), offsetVec, bend)
            @test status == true
            @test nIn == bend.refIndexOut
            # propagate from the lens plane (z=0, height h) to z=f and confirm height -> 0
            endpoint = Point3(0.0, h, 0.0) + newDir * (f / newDir[3])
            @test endpoint ≈ Point3(0.0, 0.0, f)
        end

        @testset "reverseMod! swaps indices, leaves focalLength (via generic AbstractBendType fallback)" begin
            bendCopy = OpticTrace.ParaxialLensT(f, 1.0, 1.5)
            OpticTrace.reverseMod!(bendCopy)
            @test bendCopy.refIndexIn == 1.5
            @test bendCopy.refIndexOut == 1.0
            @test bendCopy.focalLength == f
        end
    end



    @testset "modFunc tests" begin


        normal = Vec3(0.0, 0.0, 1.0)
        ray = Ray(Point(0.0, 0.0, 0.0), Vec(0.0, 0.4, sqrt(1.0 - 0.4^2)))
        dT = OpticTrace.DielectricT(1.0, 1.5)
        status, dir, nIn = OpticTrace.modFunc(ray, normal, dT)
        @test status == true
        @test nIn == dT.refIndexIn
        @test dT.refIndexOut * dir[2] ≈ dT.refIndexIn * ray.dir[2]

        dR = OpticTrace.MirrorR(1.0, 1.0)
        status, dir, nIn = OpticTrace.modFunc(ray, normal, dR)
        @test status == true
        @test nIn ≈ dR.refIndexIn
        @test dir ≈ Vec3(0.0, 0.4, -sqrt(1.0 - 0.4^2))

        @testset "CDiffuser" begin
            a = normalize(Vec3(0.1, 0.2, 0.97))
            θ = 0.15
            cd = OpticTrace.CDiffuser(tan(θ), 1.0, 1.5)
            diffRay = Ray(Point3(0.0, 0.0, 0.0), a)

            dirs = Vector{Vec3{Float64}}()
            for _ in 1:20
                status, dir, nIn = OpticTrace.modFunc(diffRay, normal, cd)
                @test status == true
                @test nIn == cd.refIndexIn # a ⋅ normal > 0 -> hitting "forward"
                @test norm(dir) ≈ 1.0
                @test dot(dir, a) >= cos(θ) - 1e-9 # stays within the diffuser cone
                push!(dirs, dir)
            end
            @test length(unique(dirs)) > 1 # scattering is randomized
        end

        @testset "NoBendIndex" begin
            nbRay = Ray(Point3(0.0, 0.0, 0.0), Vec3(0.1, 0.2, sqrt(1.0 - 0.1^2 - 0.2^2)))
            nb = OpticTrace.NoBendIndex(1.33)
            status, dir, nIn = OpticTrace.modFunc(nbRay, normal, nb)
            @test status == true
            @test dir == nbRay.dir
            @test nIn == nb.refIndexIn
        end

    end
    #=
        a test set for surfNormal methods
        The tests checks that the normal vector of different surfaces and different intersection points are correct.

        This is done by comparing the normal vector obtained from the surfNormal method
        with the expected normal vector found using the gradient of the sag methods found using and autodiff package.
        The first step for this is to find the intersection point of the ray with the surface, then find the sag value at that point,
             and then find the gradient of the sag value at that point. The normal vector is then obtained by normalizing the gradient vector.
    =#

end
