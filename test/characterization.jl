@testset "characterization.jl" begin

    @testset "Sizing" begin
        @test sizeOptic(SizeLens(5.0)) == 5.0
        @test sizeOptic(RoundAperture(1.0, 5.0)) == 5.0
        @test sizeOptic(RectAperture(1.0, 1.0, 3.0, 4.0)) ≈ norm([3.0, 4.0])

        surf = refractSphere("size_test", ORIGIN, ZAXIS, 1.0, 1.5, 0.02, 7.0, "none")
        @test sizeOptic(surf) == 7.0
    end

    @testset "Sampling" begin

        @testset "randomPointOnSquare / randomPointOnDisk" begin
            for _ in 1:100
                p = randomPointOnSquare(3.0)
                @test abs(p[1]) <= 3.0
                @test abs(p[2]) <= 3.0
                @test p[3] == 0.0
            end
            for _ in 1:100
                p = randomPointOnDisk(3.0)
                @test norm([p[1], p[2]]) <= 3.0 + 1e-9
                @test p[3] == 0.0
            end
        end

        @testset "toLocalRay" begin
            surf = refractSphere("tlr", Point3(1.0, 2.0, 3.0), ZAXIS, 1.0, 1.5, 0.02, 5.0, "none")
            ray = Ray(Point3(1.0, 2.0, 3.0), ZAXIS) # already at surf's own origin/direction
            localRay = toLocalRay(ray, surf)
            @test localRay.base ≈ ORIGIN
            @test localRay.dir ≈ ZAXIS
        end

        @testset "localRaysHexapolar" begin
            geo = [referencePlane("lrh_img", ORIGIN, ZAXIS, 1.0, 10.0, "none")]
            pupil = referencePlane("lrh_pupil", Point3(0.0, 0.0, 5.0), ZAXIS, 1.0, 3.0, "none")
            basept = Point3(0.0, 0.0, -5.0)

            rays = localRaysHexapolar(geo, basept, pupil, 1)
            @test length(rays) == 7 # 1 center ray + 6 around the single ring
            # geo[end] sits at ORIGIN with an identity local<->global transform
            @test all(r -> isapprox(r.base[3], 0.0; atol=1e-9), rays)

            @test_throws ErrorException localRaysHexapolar(geo, basept, pupil, 0)
        end
    end

    @testset "Closest approach & focal plane" begin

        @testset "distClosestApproach" begin
            ray1 = Ray(Point3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0))
            ray2 = Ray(Point3(1.0, 1.0, 1.0), Vec3(0.0, 0.0, 1.0))

            t1, t2 = OpticTrace.distClosestApproach(ray1, ray2)
            @test t1 ≈ 1.0
            @test t2 ≈ -1.0
            # hand-verified: closest points are (1,0,0) on ray1 and (1,1,0)
            # on ray2 -- the perpendicular-direction y-separation of 1.0
            closePt1 = ray1.base + t1 * ray1.dir
            closePt2 = ray2.base + t2 * ray2.dir
            @test norm(closePt1 - closePt2) ≈ 1.0

            ray3 = Ray(Point3(0.0, 1.0, 0.0), Vec3(1.0, 0.0, 0.0)) # parallel to ray1
            @test_throws ErrorException OpticTrace.distClosestApproach(ray1, ray3)
        end

        @testset "surfClosestApproach (2-ray)" begin
            ray1 = Ray(Point3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0))
            ray2 = Ray(Point3(1.0, 1.0, 1.0), Vec3(0.0, 0.0, 1.0))

            refSurf, t1 = surfClosestApproach(ray1, ray2)
            @test t1 ≈ 1.0
            @test refSurf.base.base ≈ Point3(1.0, 0.0, 0.0)
            @test refSurf.base.dir == ray1.dir
        end

        @testset "surfClosestApproach (3-ray)" begin
            ray1 = Ray(Point3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0))
            ray2 = Ray(Point3(1.0, 1.0, 1.0), Vec3(0.0, 0.0, 1.0))
            ray3 = Ray(Point3(5.0, 0.0, 0.0), ray1.dir) # parallel to ray1, offset

            refSurf3, t1_3 = surfClosestApproach(ray1, ray2, ray3)
            @test t1_3 ≈ 1.0
            @test refSurf3.base.base ≈ ray3.base + t1_3 * ray3.dir
            @test refSurf3.base.dir == ray3.dir

            ray3bad = Ray(Point3(5.0, 0.0, 0.0), YAXIS)
            @test_throws ErrorException surfClosestApproach(ray1, ray2, ray3bad)
        end

        @testset "computeRearFocalPlane" begin

            @testset "no-power (flat) system" begin
                geo = [referencePlane("rfp_flat", Point3(0.0, 0.0, 10.0), ZAXIS, 1.0, 10.0, "none")]
                status, fly, flx, yplane, xplane = computeRearFocalPlane(geo)
                @test status == 0
                @test isnan(fly)
                @test isnan(flx)
            end

            @testset "converging lens, cross-checked against a manual trace" begin
                geo = lensSinglet(Point3(0.0, 0.0, 5.0), ZAXIS, 0.05, -0.04, 3.0, 0.5, riFunc, 5.0; order="forward", lensname="RFP")
                epsilon = 0.001

                status, fly, flx, yplane, xplane = computeRearFocalPlane(geo; epsilon=epsilon)
                @test status == 0

                # reproduce computeRearFocalPlane's own documented algorithm
                # by hand, using the already-validated traceGeometryRel, to
                # cross-check its orchestration.
                statusb, bore = traceGeometryRel(Ray(ORIGIN, ZAXIS), geo)
                statusy, ytrace = traceGeometryRel(Ray(Point(0.0, epsilon, 0.0), ZAXIS), geo)
                statusx, xtrace = traceGeometryRel(Ray(Point(epsilon, 0.0, 0.0), ZAXIS), geo)
                @test statusb == 0 && statusy == 0 && statusx == 0

                rayb, rayx, rayy = bore[end].ray, xtrace[end].ray, ytrace[end].ray
                lenxzero = -rayx.base[1] / rayx.dir[1]
                lenxeps = (epsilon - rayx.base[1]) / rayx.dir[1]
                lenyzero = -rayy.base[2] / rayy.dir[2]
                lenyeps = (epsilon - rayy.base[2]) / rayy.dir[2]

                expectedFlx = (lenxzero - lenxeps) * rayx.dir[3]
                expectedFly = (lenyzero - lenyeps) * rayy.dir[3]
                expectedXplane = SVector(rayb.base[1], rayb.base[2], rayx.base[3] + lenxzero * rayx.dir[3])
                expectedYplane = SVector(rayb.base[1], rayb.base[2], rayy.base[3] + lenyzero * rayy.dir[3])

                @test flx ≈ expectedFlx
                @test fly ≈ expectedFly
                @test xplane ≈ expectedXplane
                @test yplane ≈ expectedYplane

                # sanity: a biconvex converging lens should have a finite,
                # positive focal length, with the focal plane beyond the lens
                @test flx > 0
                @test fly > 0
                @test xplane[3] > geo[end].base.base[3]
                @test yplane[3] > geo[end].base.base[3]
            end
        end

        @testset "findRFP" begin
            geo = lensSinglet(Point3(0.0, 0.0, 5.0), ZAXIS, 0.05, -0.04, 3.0, 0.5, riFunc, 5.0; order="forward", lensname="RFP2")
            epsilon = 0.001

            surfRFP, fl = findRFP(geo; epsilon=epsilon)

            statusb, bore = traceGeometryRel(Ray(ORIGIN, ZAXIS), geo)
            statusy, ytrace = traceGeometryRel(Ray(ORIGIN + epsilon * YAXIS, ZAXIS), geo)
            @test statusb == 0 && statusy == 0

            rayb, rayy = bore[end].ray, ytrace[end].ray
            expectedRFP, t1 = surfClosestApproach(rayb, rayy, semiDiam=5.0, surfname="RFPcheck")
            @test surfRFP.base.base ≈ expectedRFP.base.base

            rayyb, raybb = ytrace[begin].ray, bore[begin].ray
            expectedRPP, t2 = surfClosestApproach(rayyb, rayy, raybb, semiDiam=5.0, surfname="RPPcheck")
            expectedFl = norm(expectedRFP.base.base - expectedRPP.base.base)
            @test fl ≈ expectedFl
            @test fl > 0
        end
    end

    @testset "Spot diagrams" begin

        @testset "centroidofpoints / rmsradiusofpoints" begin
            pts = [Point(0.0, 0.0), Point(2.0, 0.0), Point(0.0, 2.0), Point(2.0, 2.0)]
            center = centroidofpoints(pts)
            @test center ≈ Point(1.0, 1.0)

            rms = rmsradiusofpoints(pts, center)
            @test rms ≈ sqrt(2.0) # every point sits at distance sqrt(2) from (1,1)
        end

        @testset "spotDiagramHex" begin
            geo = [referencePlane("sd_img", ORIGIN, ZAXIS, 1.0, 10.0, "none")]
            pupil = referencePlane("sd_pupil", Point3(0.0, 0.0, 5.0), ZAXIS, 1.0, 3.0, "none")
            basept = Point3(0.0, 0.0, -5.0)

            pts, center, rmsradius, deltaZ = spotDiagramHex(geo, basept, pupil, 1; bestFocus=false)
            @test center ≈ centroidofpoints(pts)
            @test rmsradius ≈ rmsradiusofpoints(pts, center)
            @test deltaZ == 0.0
            @test norm(center) < 1e-6 # on-axis, symmetric system

            ptsF, centerF, rmsradiusF, deltaZF = spotDiagramHex(geo, basept, pupil, 1; bestFocus=true)
            @test isfinite(deltaZF)
            @test centerF ≈ centroidofpoints(ptsF)
            @test rmsradiusF ≈ rmsradiusofpoints(ptsF, centerF)
        end
    end

    @testset "Monte Carlo & loss" begin

        @testset "traceMonteCarloRays (known bug: entirely non-functional, see TODO.md)" begin
            # traceMonteCarloRays's very first executable line builds
            # `badray = Ray((NaN, NaN, NaN), (NaN, NaN, NaN))` from raw
            # tuples, but Ray requires an actual Point{N,T}/Vec{N,T} pair --
            # this throws a MethodError unconditionally, before the
            # function ever reaches its own arguments or the ray-tracing
            # loop. Discovered while writing this test (phase 7); no call
            # to this function can currently succeed regardless of inputs.
            geo = AbstractSurface[referencePlane("mc_img", ORIGIN, ZAXIS, 1.0, 10.0, "none")]
            radiusfunc = r -> Point3(randomPointOnDisk(r)...)
            anglefunc = θ -> Vec3(0.0, 0.0, 1.0)

            @test_throws MethodError traceMonteCarloRays(radiusfunc, anglefunc, 5.0, 0.0, 50, geo)
        end

        @testset "traceLoss" begin
            amp1 = AmpData(OpticTrace.identityPol, OpticTrace.identityPol, [0.9])
            amp2 = AmpData(OpticTrace.identityPol, OpticTrace.identityPol, [0.8])
            trc = [Trace(Ray(ORIGIN, ZAXIS), 1.0, 0.0, amp1), Trace(Ray(ORIGIN, ZAXIS), 1.0, 0.0, amp2)]

            @test traceLoss(trc) ≈ 0.9 * 0.8
            @test traceLoss(trc, 2.0) ≈ 2.0 * 0.9 * 0.8
        end
    end

end
