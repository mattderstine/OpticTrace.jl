@testset "mesh_primitives.jl" begin

    @testset "SizeLens / SurfProfileConic mesh (OptSurface)" begin

        sl = SizeLens(5.0)
        spc = OpticTrace.SurfProfileConic(0.02, 0.5)

        @testset "samplePoints" begin
            pts = collect(OpticTrace.samplePoints(sl, 4))
            @test length(pts) == 16 # nvertices^2, nested r/φ grid
            @test all(p -> norm(p) <= 5.0 + 1e-9, pts)
        end

        @testset "gbWidths / gbRadius" begin
            @test OpticTrace.gbRadius(sl, spc) == 5.0
            @test OpticTrace.gbWidths(sl, spc) == SVector(10.0, 10.0, spc.curv * 5.0^2)
        end

        surf = refractSphere("mesh_test", ORIGIN, ZAXIS, 1.0, 1.5, 0.02, 5.0, "none")

        @testset "origin / radius / widths" begin
            @test GeometryBasics.origin(surf) == surf.base.base
            @test GeometryBasics.radius(surf) == surf.aperture.semiDiameter
            @test GeometryBasics.widths(surf) == SVector(10.0, 10.0, surf.profile.curv * 5.0^2)
        end

        @testset "coordinates" begin
            pts = collect(GeometryBasics.coordinates(surf, 4))
            @test length(pts) == 16
            # surf sits at ORIGIN/ZAXIS with the default ydir guess, so
            # toGlobalCoord is the identity map here -- z should match sag
            # exactly for every sampled (x,y).
            for p in pts
                @test p[3] ≈ sag(p[1], p[2], surf.profile)
            end
        end

        @testset "texturecoordinates" begin
            tex = collect(GeometryBasics.texturecoordinates(surf, 4))
            @test length(tex) == 16
            @test all(t -> 0.0 <= t[1] <= 1.0 && 0.0 <= t[2] <= 1.0, tex)
        end

        @testset "faces" begin
            f = GeometryBasics.faces(surf, 4)
            @test length(f) > 0
        end

        @testset "inOrOut" begin
            # refIndexIn=1.0 < refIndexOut=1.5 -> "in"
            @test OpticTrace.inOrOut(surf) == -1
            mirrorSurf = reflectConic("mesh_mirror", ORIGIN, ZAXIS, 1.5, 1.0, 0.02, 0.5, 5.0, "none")
            @test OpticTrace.inOrOut(mirrorSurf) == 1
        end

        @testset "normals (OptSurface)" begin
            ns = collect(GeometryBasics.normals(surf, 4))
            @test length(ns) == 16
            @test all(n -> norm(n) ≈ 1.0, ns)
        end

        @testset "normals (AbstractSurface, ModelSurface -- known bug)" begin
            # A ModelSurface with a SizeLens aperture so samplePoints
            # succeeds and the failure is isolated to inOrOut, which has
            # no method outside OptSurface (see TODO.md).
            ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir = updateCoordChange(ORIGIN, ZAXIS)
            modelSurf = ModelSurface("mesh_model", SurfBase(ORIGIN, ZAXIS, ydir), SizeLens(5.0),
                NoProfile(0.), 1.0, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir, :blue)
            @test_throws MethodError GeometryBasics.normals(modelSurf, 4)
        end
    end

    @testset "NoProfile / SizeLens mesh (referencePlane surfaces.jl overloads)" begin
        refSurf = referencePlane("mesh_ref", ORIGIN, ZAXIS, 1.0, 7.0, "none")

        @test GeometryBasics.radius(refSurf) == 7.0
        @test GeometryBasics.widths(refSurf) == SVector(14.0, 14.0, 0.0)

        pts = collect(GeometryBasics.coordinates(refSurf, 4))
        @test all(p -> p[3] == 0.0, pts)

        # refIndexIn == refIndexOut for a referencePlane -> inOrOut's "else"
        # branch, dir = -1; toGlobalDir is identity here (ORIGIN/ZAXIS).
        ns = collect(GeometryBasics.normals(refSurf, 4))
        @test all(n -> n ≈ -ZAXIS, ns)
    end

    @testset "Washer mesh" begin
        _, toGlobalCoordId, _, toGlobalDirId, _ = updateCoordChange(ORIGIN, ZAXIS)
        washer = OpticTrace.Washer(ORIGIN, ZAXIS, 3.0, toGlobalCoordId, toGlobalDirId)

        @testset "samplePoints / origin / coordinates / normals (working)" begin
            pts = collect(OpticTrace.samplePoints(washer, 4))
            @test length(pts) == 16
            @test all(p -> 3.0 - 1e-9 <= norm(p) <= 3.0 * 1.3 + 1e-9, pts)

            @test GeometryBasics.origin(washer) == ORIGIN

            gpts = collect(GeometryBasics.coordinates(washer, 4))
            @test all(p -> p[3] == 0.0, gpts)

            ns = collect(GeometryBasics.normals(washer, 4))
            @test all(n -> n ≈ ZAXIS, ns)
        end

        @testset "gbRadius (working) vs gbWidths/radius/widths (broken)" begin
            @test OpticTrace.gbRadius(washer, NoProfile(0.)) == washer.semiDiameter

            # gbWidths references misspelled fields SemiDiameter/SemiDiamater
            # (real field: semiDiameter) -- see TODO.md.
            @test_throws FieldError OpticTrace.gbWidths(washer, NoProfile(0.))

            # GeometryBasics.radius/widths(::Washer) both pass a bare number
            # to gbRadius/gbWidths instead of the Washer itself.
            @test_throws MethodError GeometryBasics.radius(washer)
            @test_throws MethodError GeometryBasics.widths(washer)
        end
    end

    @testset "Disk mesh" begin
        _, toGlobalCoordId, _, toGlobalDirId, _ = updateCoordChange(ORIGIN, ZAXIS)
        disk = OpticTrace.Disk(ORIGIN, ZAXIS, 3.0, toGlobalCoordId, toGlobalDirId)

        @testset "samplePoints / origin / coordinates / normals (working)" begin
            pts = collect(OpticTrace.samplePoints(disk, 4))
            @test length(pts) == 16
            @test all(p -> norm(p) <= 3.0 + 1e-9, pts)

            @test GeometryBasics.origin(disk) == ORIGIN

            gpts = collect(GeometryBasics.coordinates(disk, 4))
            @test all(p -> p[3] == 0.0, gpts)

            ns = collect(GeometryBasics.normals(disk, 4))
            @test all(n -> n ≈ ZAXIS, ns)
        end

        @testset "gbRadius (working) vs gbWidths/radius/widths (broken)" begin
            @test OpticTrace.gbRadius(disk, NoProfile(0.)) == disk.semiDiameter

            @test_throws FieldError OpticTrace.gbWidths(disk, NoProfile(0.))
            @test_throws MethodError GeometryBasics.radius(disk)
            @test_throws MethodError GeometryBasics.widths(disk)
        end
    end

    @testset "RectAperture / NoProfile mesh" begin
        rectFinite = RectAperture(0.5, 0.5, 5.0, 6.0)
        @test OpticTrace.gbWidths(rectFinite, NoProfile(0.)) == SVector(10.0, 12.0, 0.0)

        rectInfinite = RectAperture(0.5, 0.5, ∞, ∞)
        @test OpticTrace.gbWidths(rectInfinite, NoProfile(0.)) == SVector(1.0, 1.0, 0.0) # falls back to 2wo, 2lo

        # gbRadius: the both-infinite branch works...
        @test OpticTrace.gbRadius(rectInfinite, NoProfile(0.)) ≈ sqrt(0.5^2 + 0.5^2)

        # ...but the finite (everyday) branch references a nonexistent field
        # `a.clear` (real fields are `wclear`/`lclear`). Confirmed by direct
        # reading of the source: this is the *opposite* of which branch
        # TODO.md's prose describes as broken -- the bug is in the finite
        # case, not the both-infinite case, so it affects ordinary
        # rectangular apertures, not just an edge case.
        @test_throws FieldError OpticTrace.gbRadius(rectFinite, NoProfile(0.))
    end

end
