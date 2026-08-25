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

        @testset "normals (AbstractSurface, ModelSurface)" begin
            # A ModelSurface with a SizeLens aperture, since samplePoints
            # only has methods for SizeLens/Washer/Disk apertures.
            ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir = updateCoordChange(ORIGIN, ZAXIS)
            modelSurf = ModelSurface("mesh_model", SurfBase(ORIGIN, ZAXIS, ydir), SizeLens(5.0),
                NoProfile(0.), 1.0, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir, :blue)
            ns = collect(GeometryBasics.normals(modelSurf, 4))
            @test length(ns) == 16
            @test all(n -> norm(n) ≈ 1.0, ns)
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

        @testset "gbRadius / gbWidths / radius / widths" begin
            @test OpticTrace.gbRadius(washer, NoProfile(0.)) == washer.semiDiameter
            @test OpticTrace.gbWidths(washer, NoProfile(0.)) == SVector(2washer.semiDiameter, 2washer.semiDiameter, 0.)

            @test GeometryBasics.radius(washer) == washer.semiDiameter
            @test GeometryBasics.widths(washer) == SVector(2washer.semiDiameter, 2washer.semiDiameter, 0.)
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

        @testset "gbRadius / gbWidths / radius / widths" begin
            @test OpticTrace.gbRadius(disk, NoProfile(0.)) == disk.semiDiameter
            @test OpticTrace.gbWidths(disk, NoProfile(0.)) == SVector(2disk.semiDiameter, 2disk.semiDiameter, 0.)

            @test GeometryBasics.radius(disk) == disk.semiDiameter
            @test GeometryBasics.widths(disk) == SVector(2disk.semiDiameter, 2disk.semiDiameter, 0.)
        end
    end

    @testset "RectAperture / NoProfile mesh" begin
        rectFinite = RectAperture(0.5, 0.5, 5.0, 6.0)
        @test OpticTrace.gbWidths(rectFinite, NoProfile(0.)) == SVector(10.0, 12.0, 0.0)

        rectInfinite = RectAperture(0.5, 0.5, ∞, ∞)
        @test OpticTrace.gbWidths(rectInfinite, NoProfile(0.)) == SVector(1.0, 1.0, 0.0) # falls back to 2wo, 2lo

        @test OpticTrace.gbRadius(rectInfinite, NoProfile(0.)) ≈ sqrt(0.5^2 + 0.5^2)
        @test OpticTrace.gbRadius(rectFinite, NoProfile(0.)) ≈ sqrt(rectFinite.wclear^2 + rectFinite.lclear^2)
    end

end
