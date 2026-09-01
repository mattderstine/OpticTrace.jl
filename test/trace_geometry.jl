@testset "trace_geometry.jl" begin

    @testset "traceSurf (OptSurface)" begin

        @testset "normal refraction (status 0)" begin
            surf = refractSphere("os_normal", ORIGIN, ZAXIS, 1.0, 1.5, 0.02, 5.0, "none")
            ray = Ray(Point3(0.3, 0.2, -1.0), Vec3(0.0, 0.0, 1.0))

            status, trc = OpticTrace.traceSurf(ray, surf)
            @test status == 0

            # cross-check against the individual pieces directly -- surf
            # sits at ORIGIN/ZAXIS with an identity local<->global
            # transform, so local == global coordinates here.
            expectedDelta = OpticTrace.deltaToSurf(ray, surf.profile)
            expectedBase = ray.base + ray.dir * expectedDelta
            expectedNormal = OpticTrace.surfNormal(expectedBase, surf.profile)
            _, expectedDir, expectedNIn = OpticTrace.modFunc(Ray(expectedBase, ray.dir), expectedNormal, surf.mod)

            @test trc.delta ≈ expectedDelta
            @test trc.ray.base ≈ expectedBase
            @test trc.ray.dir ≈ expectedDir
            @test trc.nIn ≈ expectedNIn
        end

        @testset "TIR (status 2)" begin
            surf = refractConic("os_tir", ORIGIN, ZAXIS, 1.5, 1.0, 0.0, 0.0, 5.0, "none")
            θ = deg2rad(60.0) # critical angle for 1.5->1.0 is asin(1/1.5) ≈ 41.8°
            ray = Ray(Point3(0.0, 0.0, -1.0), Vec3(sin(θ), 0.0, cos(θ)))

            status, trc = OpticTrace.traceSurf(ray, surf)
            @test status == 2
            @test trc.nIn == 1.5
        end

        @testset "miss (status 1)" begin
            surf = refractConic("os_miss", ORIGIN, ZAXIS, 1.0, 1.5, 0.0, 0.0, 5.0, "none")
            ray = Ray(Point3(0.0, 0.0, -1.0), Vec3(1.0, 0.0, 0.0)) # parallel to the flat surface -> deltaToSurf returns NaN

            status, trc = OpticTrace.traceSurf(ray, surf)
            @test status == 1
            @test isnan(trc.delta)
        end

        @testset "traceSurf! matches traceSurf" begin
            surf = refractSphere("os_inplace", ORIGIN, ZAXIS, 1.0, 1.5, 0.02, 5.0, "none")
            ray = Ray(Point3(0.3, 0.2, -1.0), Vec3(0.0, 0.0, 1.0))

            status1, trc1 = OpticTrace.traceSurf(ray, surf)
            trcBuf = Trace(ray, 1.0, 0.0, OpticTrace.identityAmpMats())
            status2, trc2 = OpticTrace.traceSurf!(trcBuf, ray, surf)

            @test status1 == status2
            @test trc1.ray.base ≈ trc2.ray.base
            @test trc1.ray.dir ≈ trc2.ray.dir
            @test trc1.nIn ≈ trc2.nIn
            @test trc1.delta ≈ trc2.delta
            @test trc2 === trcBuf
        end

        @testset "paraxial thin lens, full traceSurf pipeline (TODO.md #4 / FIXED.md)" begin
            # Integration test for the surfNormal(::ParaxialProfile) ->
            # s.toGlobalDir -> modFunc(::ParaxialLensT) pipeline: unlike
            # the unit-level modFunc test in test/optics.jl, this exercises
            # the real traceSurf coordinate-transform wiring end to end,
            # including a surface not sitting at the identity transform.
            focalLength = 50.0
            profile = OpticTrace.ParaxialProfile(0.0)
            bend = OpticTrace.ParaxialLensT(focalLength, 1.0, 1.0)
            basept = Point3(1.0, 2.0, 3.0)
            dirTilt = normalize(Vec3(0.1, 0.0, 1.0))
            ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                OpticTrace.updateCoordChange(basept, dirTilt, nothing)
            surf = OptSurface("paraxial", SurfBase(basept, dirTilt, ydir), SizeLens(10.0), profile,
                bend, OpticTrace.getAmpParams("none"; OpticTrace.attributesSurfaces),
                toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir, :cyan3)

            h = 2.0 # local height above the lens's own optical axis
            localOffset = toGlobalDir(Vec3(0.0, h, 0.0))
            rayBase = basept + localOffset - 10.0 * dirTilt # well upstream, parallel to the lens axis
            ray = Ray(rayBase, dirTilt)

            status, trc = OpticTrace.traceSurf(ray, surf)
            @test status == 0

            # propagate from the lens plane to the back focal plane
            # (basept + f*dirTilt) and confirm the ray lands back on-axis
            # -- travel along the OUTGOING ray direction (not dirTilt,
            # which it's no longer parallel to) far enough that its
            # projection onto dirTilt reaches the focal plane
            focalPoint = basept + focalLength * dirTilt
            t = dot(focalPoint - trc.ray.base, dirTilt) / dot(trc.ray.dir, dirTilt)
            endpoint = trc.ray.base + trc.ray.dir * t
            @test endpoint ≈ focalPoint atol=1e-9
        end
    end

    @testset "traceSurf (ModelSurface)" begin

        @testset "normal pass-through (status 0)" begin
            modelSurf = roundAperture("ms_normal", ORIGIN, ZAXIS, 1.0, 1.0, 5.0)
            ray = Ray(Point3(3.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0)) # within clear aperture, outside obscuration

            status, trc = OpticTrace.traceSurf(ray, modelSurf)
            @test status == 0
            @test trc.ray.dir == ray.dir # pass-through: direction unchanged
            @test trc.ray.base ≈ Point3(3.0, 0.0, 0.0)
            @test trc.nIn == modelSurf.refIndex
        end

        @testset "aperture-clipped (status 3)" begin
            modelSurf = roundAperture("ms_clip", ORIGIN, ZAXIS, 1.0, 1.0, 5.0)
            ray = Ray(Point3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0)) # inside the central obscuration

            status, trc = OpticTrace.traceSurf(ray, modelSurf)
            @test status == 3
        end

        @testset "miss (status 1)" begin
            modelSurf = roundAperture("ms_miss", ORIGIN, ZAXIS, 1.0, 1.0, 5.0)
            ray = Ray(Point3(3.0, 0.0, -1.0), Vec3(1.0, 0.0, 0.0)) # parallel to the flat surface -> NaN delta

            status, trc = OpticTrace.traceSurf(ray, modelSurf)
            @test status == 1
            @test isnan(trc.delta)
        end

        @testset "traceSurf! matches traceSurf" begin
            modelSurf = roundAperture("ms_inplace", ORIGIN, ZAXIS, 1.0, 1.0, 5.0)
            ray = Ray(Point3(3.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0))

            status1, trc1 = OpticTrace.traceSurf(ray, modelSurf)
            trcBuf = Trace(ray, 1.0, 0.0, OpticTrace.identityAmpMats())
            status2, trc2 = OpticTrace.traceSurf!(trcBuf, ray, modelSurf)

            @test status1 == status2
            @test trc1.ray.base ≈ trc2.ray.base
            @test trc1.ray.dir ≈ trc2.ray.dir
            @test trc1.nIn ≈ trc2.nIn
            @test trc2 === trcBuf
        end
    end

    @testset "traceGeometry / traceGeometry! / traceGeometryRel / traceGeometryRel!" begin

        geo = [[referencePlane("object", ORIGIN, ZAXIS, 1.0, 10.0, "none")]
            lensSinglet(Point3(0.0, 0.0, 5.0), ZAXIS, 0.02, -0.015, 3.0, 0.5, riFunc, 5.0; order="forward", lensname="TG");
            [referencePlane("image", Point3(0.0, 0.0, 50.0), ZAXIS, 1.0, 10.0, "none")]
        ]

        ray = Ray(Point3(0.5, 0.3, -5.0), Vec3(0.0, 0.0, 1.0))

        @testset "traceGeometry" begin
            status, trc = traceGeometry(ray, geo)
            @test status == 0
            @test length(trc) == length(geo) + 1
            @test trc[1].ray == ray
        end

        @testset "traceGeometry!" begin
            status1, trc1 = traceGeometry(ray, geo)
            # traceGeometry! mutates each trc[i] in place via Trace! -- the
            # buffer must be pre-filled with real Trace objects, not left
            # `undef` (Trace! reads trc[i] before overwriting its fields).
            # traceGeometry! requires exactly Vector{Trace} (parametric
            # containers are invariant in Julia, so a comprehension's
            # narrower Vector{Trace{Float64,AmpData{Float64}}} won't match).
            trcBuf = Trace[Trace(ray, 1.0, 0.0, OpticTrace.identityAmpMats()) for _ in 1:(length(geo)+1)]
            status2, len2 = traceGeometry!(trcBuf, ray, geo)

            @test status2 == status1
            @test len2 == length(trc1)
            for i in 1:len2
                @test trcBuf[i].ray.base ≈ trc1[i].ray.base
                @test trcBuf[i].ray.dir ≈ trc1[i].ray.dir
                @test trcBuf[i].nIn ≈ trc1[i].nIn
                @test trcBuf[i].delta ≈ trc1[i].delta
            end
        end

        @testset "traceGeometryRel converts geo[1]-local to global" begin
            # use a geo offset from the origin so the local->global
            # conversion is actually exercised (not an identity map).
            geoOffset = [[referencePlane("object2", Point3(1.0, 2.0, 3.0), ZAXIS, 1.0, 10.0, "none")]
                lensSinglet(Point3(1.0, 2.0, 8.0), ZAXIS, 0.02, -0.015, 3.0, 0.5, riFunc, 5.0; order="forward", lensname="TG2");
                [referencePlane("image2", Point3(1.0, 2.0, 53.0), ZAXIS, 1.0, 10.0, "none")]
            ]
            localRay = Ray(Point3(0.5, 0.3, -5.0), Vec3(0.0, 0.0, 1.0))
            surf1 = geoOffset[1]
            convertedRay = Ray(surf1.toGlobalCoord(localRay.base), surf1.toGlobalDir(localRay.dir))

            statusRel, trcRel = traceGeometryRel(localRay, geoOffset)
            statusGlobal, trcGlobal = traceGeometry(convertedRay, geoOffset)

            @test statusRel == statusGlobal
            @test trcRel[end].ray.base ≈ trcGlobal[end].ray.base
            @test trcRel[end].ray.dir ≈ trcGlobal[end].ray.dir

            @testset "traceGeometryRel!" begin
                trcBufRel = Trace[Trace(localRay, 1.0, 0.0, OpticTrace.identityAmpMats()) for _ in 1:(length(geoOffset)+1)]
                statusRelBang, lenRelBang = traceGeometryRel!(trcBufRel, localRay, geoOffset)

                @test statusRelBang == statusRel
                @test trcBufRel[lenRelBang].ray.base ≈ trcRel[end].ray.base
                @test trcBufRel[lenRelBang].ray.dir ≈ trcRel[end].ray.dir
            end
        end
    end

    @testset "Trace!" begin
        t = Trace(Ray(ORIGIN, ZAXIS), 1.0, 0.0, OpticTrace.identityAmpMats())
        newRay = Ray(Point3(1.0, 2.0, 3.0), YAXIS)
        result = OpticTrace.Trace!(t, newRay, 1.5, 2.0, OpticTrace.identityAmpMats())
        @test result === t
        @test t.ray == newRay
        @test t.nIn == 1.5
        @test t.delta == 2.0
    end

    @testset "surfAmpFunc" begin
        dT = DielectricT(1.0, 1.5)
        ampMats, dirOut = OpticTrace.surfAmpFunc(Vec3(0.0, 0.0, 1.0), Vec3(0.1, 0.0, 0.99),
            Vec3(0.0, 0.0, 1.0), Point3(0.0, 0.0, 0.0), dT, OpticTrace.AmpParam("x"))
        @test ampMats === OpticTrace.identityAmpMats()
        @test dirOut == Vec3(0.1, 0.0, 0.99)
    end

    @testset "getAmpParams" begin
        testDict = Dict{String,Any}()

        a1 = OpticTrace.getAmpParams("newcoat"; attributesSurfaces=testDict)
        @test a1 isa OpticTrace.AmpParam
        @test a1.type == "newcoat"
        @test haskey(testDict, "newcoat")

        a2 = OpticTrace.getAmpParams("newcoat"; attributesSurfaces=testDict)
        @test a2 === a1 # reused from the cache, not rebuilt

        nap = NoAmpParam("passthrough")
        a3 = OpticTrace.getAmpParams(nap; attributesSurfaces=testDict)
        @test a3 === nap # identity pass-through for an already-built AbstractAmplitudeParam
    end

end
