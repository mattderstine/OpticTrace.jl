@testset "surface_manipulation.jl" begin

    @testset "reverseBase!" begin
        b = SurfBase(Point3(1.0, 2.0, 3.0), ZAXIS, YAXIS)
        bpoint = Point3(0.0, 0.0, 0.0)
        epoint = Point3(10.0, 0.0, 0.0)

        result = OpticTrace.reverseBase!(b, bpoint, epoint)
        @test result === b
        @test b.base ≈ bpoint + epoint - Point3(1.0, 2.0, 3.0)
        @test b.dir == ZAXIS # unchanged
        @test b.ydir == YAXIS # unchanged
    end

    @testset "reverseProfile!" begin

        @testset "SurfProfileConic (working)" begin
            p = OpticTrace.SurfProfileConic(0.02, 0.5)
            result = OpticTrace.reverseProfile!(p)
            @test result === p
            @test p.curv == -0.02
            @test p.ϵ == 0.5
        end

        @testset "SurfProfileAsphere (working, via generic fallback)" begin
            p = OpticTrace.SurfProfileAsphere(0.02, 0.5, [0.001, 0.0, 0.0001])
            OpticTrace.reverseProfile!(p)
            @test p.curv == -0.02
            @test p.ϵ == 0.5
            @test p.a == [-0.001, 0.0, -0.0001]
        end

        @testset "SurfProfileEvenAsphere (working, via generic fallback)" begin
            p = OpticTrace.SurfProfileEvenAsphere(0.02, 0.5, [0.001, 0.0001])
            OpticTrace.reverseProfile!(p)
            @test p.curv == -0.02
            @test p.ϵ == 0.5
            @test p.a == [-0.001, -0.0001]
        end

        @testset "SurfProfileCyl (working)" begin
            p = OpticTrace.SurfProfileCyl(0.02, 0.5, [0.001])
            OpticTrace.reverseProfile!(p)
            @test p.curv == -0.02
            @test p.ϵ == 0.5
            @test p.a == [-0.001]
        end

        @testset "SurfProfileSphere (working)" begin
            p = OpticTrace.SurfProfileSphere(0.02)
            result = OpticTrace.reverseProfile!(p)
            @test result === p
            @test p.curv == -0.02
        end

        @testset "SurfProfileToroid (working)" begin
            p = OpticTrace.SurfProfileToroid(0.3, 0.2)
            result = OpticTrace.reverseProfile!(p)
            @test result === p
            @test p.curvY == -0.3
            @test p.curvX == -0.2
        end

        @testset "SurfProfileOAConic (dedicated method, offset not handled -- see TODO.md)" begin
            pOA = OpticTrace.SurfProfileOAConic(0.02, 0.5, Vec3(1.0, 2.0, 3.0))
            result = OpticTrace.reverseProfile!(pOA)
            @test result === pOA
            @test pOA.curv == -0.02
            @test pOA.ϵ == 0.5
            @test pOA.offset == Vec3(1.0, 2.0, 3.0) # unchanged -- known-incomplete, see TODO.md
        end

        @testset "NoProfile (working)" begin
            pNP = NoProfile(0.0)
            result = OpticTrace.reverseProfile!(pNP)
            @test result === pNP
            @test pNP.curv == 0.0
        end
    end

    @testset "reverseMod!" begin
        d = DielectricT(1.0, 1.5)
        result = OpticTrace.reverseMod!(d)
        @test result === d
        @test d.refIndexIn == 1.5
        @test d.refIndexOut == 1.0

        m = MirrorR(1.0, 1.0)
        OpticTrace.reverseMod!(m)
        @test m.refIndexIn == 1.0
        @test m.refIndexOut == 1.0
    end

    @testset "reverseGeo / reverseSurface!" begin

        base = Point3(0.0, 0.0, 0.0)
        dir = ZAXIS
        curv1, curv2 = 0.02, -0.015
        thick = 3.0

        @testset "OptSurface-only geo (working)" begin
            lens = lensSinglet(base, dir, curv1, curv2, thick, 0.5, riFunc, 5.0; order="forward", lensname="RG")
            beginpoint = lens[1].base.base
            endpoint = lens[2].base.base

            reversed = reverseGeo(lens)
            @test length(reversed) == 2

            # order is flipped
            @test reversed[1].surfname == lens[2].surfname
            @test reversed[2].surfname == lens[1].surfname

            # positions reflected about the midpoint of beginpoint/endpoint
            @test reversed[1].base.base ≈ beginpoint + endpoint - lens[2].base.base
            @test reversed[2].base.base ≈ beginpoint + endpoint - lens[1].base.base

            # curvature negated
            @test reversed[1].profile.curv ≈ -lens[2].profile.curv
            @test reversed[2].profile.curv ≈ -lens[1].profile.curv

            # refractive indices swapped
            @test reversed[1].mod.refIndexIn ≈ lens[2].mod.refIndexOut
            @test reversed[1].mod.refIndexOut ≈ lens[2].mod.refIndexIn
            @test reversed[2].mod.refIndexIn ≈ lens[1].mod.refIndexOut
            @test reversed[2].mod.refIndexOut ≈ lens[1].mod.refIndexIn

            # reverseGeo deep-copies -- the original geo is untouched
            @test lens[1].profile.curv == curv1
            @test lens[2].profile.curv == curv2
        end

        @testset "direction mismatch errors" begin
            s1 = refractSphere("dm1", ORIGIN, ZAXIS, 1.0, 1.5, 0.02, 5.0, "none")
            s2 = refractSphere("dm2", Point3(0.0, 0.0, 3.0), YAXIS, 1.5, 1.0, -0.015, 5.0, "none")
            @test_throws ErrorException reverseGeo([s1, s2])
        end

        @testset "geo containing a ModelSurface" begin
            mixedGeo = AbstractSurface[
                refractSphere("mg1", ORIGIN, ZAXIS, 1.0, 1.5, 0.02, 5.0, "none"),
                roundAperture("mg2", Point3(0.0, 0.0, 3.0), ZAXIS, 1.5, 0.5, 5.0)
            ]
            beginpoint = mixedGeo[1].base.base
            endpoint = mixedGeo[2].base.base

            reversed = reverseGeo(mixedGeo)
            @test length(reversed) == 2

            # order is flipped
            @test reversed[1].surfname == mixedGeo[2].surfname
            @test reversed[2].surfname == mixedGeo[1].surfname

            # positions reflected about the midpoint of beginpoint/endpoint
            @test reversed[1].base.base ≈ beginpoint + endpoint - mixedGeo[2].base.base
            @test reversed[2].base.base ≈ beginpoint + endpoint - mixedGeo[1].base.base

            # ModelSurface (reversed[1]): NoProfile reversal is a no-op,
            # refIndex has no in/out pair so it's left untouched
            @test reversed[1].refIndex == mixedGeo[2].refIndex

            # OptSurface (reversed[2]): curvature negated, indices swapped
            @test reversed[2].profile.curv ≈ -mixedGeo[1].profile.curv
            @test reversed[2].mod.refIndexIn ≈ mixedGeo[1].mod.refIndexOut
            @test reversed[2].mod.refIndexOut ≈ mixedGeo[1].mod.refIndexIn
        end
    end

    @testset "thickGeo" begin
        geo = [
            referencePlane("tg_o", ORIGIN, ZAXIS, 1.0, 10.0, "none"),
            referencePlane("tg_i", Point3(0.0, 0.0, 25.0), ZAXIS, 1.0, 10.0, "none"),
        ]
        @test thickGeo(geo) ≈ 25.0

        # off-axis shifts shouldn't matter -- the result is projected onto
        # the geometry's own initial direction.
        geoShifted = [
            referencePlane("tg_o2", Point3(1.0, 2.0, 0.0), ZAXIS, 1.0, 10.0, "none"),
            referencePlane("tg_i2", Point3(-3.0, 5.0, 25.0), ZAXIS, 1.0, 10.0, "none"),
        ]
        @test thickGeo(geoShifted) ≈ 25.0

        badGeo = [
            referencePlane("tg_o3", ORIGIN, ZAXIS, 1.0, 10.0, "none"),
            refractSphere("tg_i3", Point3(0.0, 0.0, 3.0), YAXIS, 1.0, 1.5, 0.02, 5.0, "none"),
        ]
        @test_throws ErrorException thickGeo(badGeo)
    end

end
