@testset "surface_builders.jl" begin

    coating = "test_coating"

    @testset "Single-surface builders" begin

        @testset "refractConic" begin
            surf = refractConic("s1", ORIGIN, ZAXIS, 1.0, 1.5, 0.02, 0.5, 5.0, coating)
            @test surf isa OptSurface
            @test surf.profile isa OpticTrace.SurfProfileConic
            @test surf.profile.curv == 0.02
            @test surf.profile.ϵ == 0.5
            @test surf.mod isa DielectricT
            @test surf.mod.refIndexIn == 1.0
            @test surf.mod.refIndexOut == 1.5
            @test surf.aperture isa SizeLens
            @test surf.aperture.semiDiameter == 5.0
            @test surf.coating isa OpticTrace.AmpParam
            @test surf.coating.type == coating
            @test surf.base.base == ORIGIN
            @test surf.base.dir == ZAXIS
            @test surf.toGlobalCoord(ORIGIN) ≈ ORIGIN
        end

        @testset "refractSphere" begin
            surf = refractSphere("s2", ORIGIN, ZAXIS, 1.0, 1.5, 0.02, 5.0, coating)
            @test surf.profile isa OpticTrace.SurfProfileConic
            @test surf.profile.curv == 0.02
            @test surf.profile.ϵ == 1.0 # refractSphere == refractConic with ϵ=1
            @test surf.mod isa DielectricT
            @test surf.mod.refIndexIn == 1.0
            @test surf.mod.refIndexOut == 1.5
        end

        @testset "reflectConic" begin
            surf = reflectConic("s3", ORIGIN, ZAXIS, 1.0, 1.0, 0.02, 0.5, 5.0, coating)
            @test surf.profile isa OpticTrace.SurfProfileConic
            @test surf.profile.curv == 0.02
            @test surf.profile.ϵ == 0.5
            @test surf.mod isa MirrorR
            @test surf.mod.refIndexIn == 1.0
            @test surf.mod.refIndexOut == 1.0
        end

        @testset "cDiffuser" begin
            θ = 0.1
            surf = cDiffuser("s4", ORIGIN, ZAXIS, 1.0, 1.0, θ, 5.0, coating)
            @test surf.profile isa NoProfile
            @test surf.mod isa CDiffuser
            @test surf.mod.tanθ ≈ tan(θ)
            @test surf.mod.refIndexIn == 1.0
            @test surf.mod.refIndexOut == 1.0
            @test surf.coating isa OpticTrace.AmpParam
            @test surf.coating.type == coating
        end

        @testset "reflectOAConic" begin
            offset = Vec3(0.0, 0.5, -1.0)

            # calling without an explicit attributesSurfaces hits the broken
            # default keyword value (`attributeSurfaces`, missing an "s"),
            # an undefined variable -- see TODO.md.
            @test_throws UndefVarError reflectOAConic("s5", ORIGIN, ZAXIS, YAXIS, offset, 1.0, 1.0, 0.02, 0.0, 5.0, coating)

            # passing it explicitly avoids the bug
            surf = reflectOAConic("s5", ORIGIN, ZAXIS, YAXIS, offset, 1.0, 1.0, 0.02, 0.0, 5.0, coating; attributesSurfaces=attributesSurfaces)
            @test surf.profile isa SurfProfileOAConic
            @test surf.profile.curv == 0.02
            @test surf.profile.ϵ == 0.0
            @test surf.profile.offset == offset
            @test surf.mod isa MirrorR
            @test surf.mod.refIndexIn == 1.0
            @test surf.mod.refIndexOut == 1.0
        end

        @testset "reflectOAP" begin
            c = 0.05
            surf = reflectOAP("s6", ORIGIN, ZAXIS, YAXIS, 1.0, 1.0, c, 5.0, coating)
            @test surf.profile isa SurfProfileOAConic
            @test surf.profile.curv == c
            @test surf.profile.ϵ == 0.0
            @test surf.profile.offset ≈ Vec3(0.0, 1 / c, -0.5 / c)
            @test surf.mod isa MirrorR
        end

        @testset "refractAsphere" begin
            asph = [0.001, 0.0, 0.0001]
            surf = refractAsphere("s7", ORIGIN, ZAXIS, 1.0, 1.5, 0.02, 0.5, asph, 5.0, coating)
            @test surf.profile isa OpticTrace.SurfProfileAsphere
            @test surf.profile.curv == 0.02
            @test surf.profile.ϵ == 0.5
            @test surf.profile.a == asph
            @test surf.mod isa DielectricT
            @test surf.mod.refIndexIn == 1.0
            @test surf.mod.refIndexOut == 1.5
        end

        @testset "refractEvenAsphere" begin
            asph = [0.001, 0.0001]
            surf = OpticTrace.refractEvenAsphere("s8", ORIGIN, ZAXIS, 1.0, 1.5, 0.02, 0.5, asph, 5.0, coating)
            @test surf.profile isa OpticTrace.SurfProfileEvenAsphere
            @test surf.profile.curv == 0.02
            @test surf.profile.ϵ == 0.5
            @test surf.profile.a == asph
            @test surf.mod isa DielectricT
        end

        @testset "reflectAsphere" begin
            asph = [0.001, 0.0, 0.0001]
            surf = reflectAsphere("s9", ORIGIN, ZAXIS, 1.0, 1.0, 0.02, 0.5, asph, 5.0, coating)
            @test surf.profile isa OpticTrace.SurfProfileAsphere
            @test surf.profile.a == asph
            @test surf.mod isa MirrorR
        end

        @testset "referencePlane" begin
            surf = referencePlane("s10", ORIGIN, ZAXIS, 1.33, 5.0, coating)
            @test surf.profile isa NoProfile
            @test surf.mod isa NoBendIndex
            @test surf.mod.refIndexIn == 1.33
            @test surf.mod.refIndexOut == 1.33
            @test surf.coating isa NoAmpParam
            @test surf.coating.type == coating
        end

        @testset "planeMirror" begin
            surf = planeMirror("s11", ORIGIN, ZAXIS, 1.0, 1.0, 5.0, coating)
            @test surf.profile isa OpticTrace.SurfProfileConic
            @test surf.profile.curv == 0.0
            @test surf.profile.ϵ == 0.0
            @test surf.mod isa MirrorR
            @test surf.mod.refIndexIn == 1.0
            @test surf.mod.refIndexOut == 1.0
        end
    end

    @testset "Lens builders" begin

        base = ORIGIN
        dir = ZAXIS
        thick = 3.0
        semiDiam = 5.0
        wl = 0.5

        @testset "lensSinglet" begin
            curv1, curv2 = 0.02, -0.015

            lensFwd = lensSinglet(base, dir, curv1, curv2, thick, wl, riFunc, semiDiam; order="forward", lensname="L")
            @test length(lensFwd) == 2
            @test lensFwd[1].profile.curv == curv1
            @test lensFwd[2].profile.curv == curv2
            @test lensFwd[1].base.base == base
            @test lensFwd[2].base.base ≈ base + thick .* dir
            @test lensFwd[1].mod.refIndexIn == refIndexDefault
            @test lensFwd[1].mod.refIndexOut == riFunc(wl)
            @test lensFwd[2].mod.refIndexIn == riFunc(wl)
            @test lensFwd[2].mod.refIndexOut == refIndexDefault

            # order="reverse" is a known bug: passes a stray Julia compiler
            # internal (Base.compute_assumed_setting) instead of `coating`
            # to the second refractSphere call -- see TODO.md.
            @test_throws MethodError lensSinglet(base, dir, curv1, curv2, thick, wl, riFunc, semiDiam; order="reverse", lensname="L")
        end

        @testset "lensASinglet" begin
            curv1, curv2 = 0.02, -0.015
            ϵ1, ϵ2 = 0.5, 0.5
            a1 = [0.001, 0.0, 0.0001]
            a2 = [0.0005, 0.0, 0.00005]

            lensFwd = lensASinglet(base, dir, curv1, ϵ1, a1, curv2, ϵ2, a2, thick, wl, riFunc, semiDiam; order="forward")
            @test lensFwd[1].profile.curv == curv1
            @test lensFwd[1].profile.a == a1
            @test lensFwd[2].profile.curv == curv2
            @test lensFwd[2].profile.a == a2

            lensRev = lensASinglet(base, dir, curv1, ϵ1, a1, curv2, ϵ2, a2, thick, wl, riFunc, semiDiam; order="reverse")
            @test lensRev[1].profile.curv ≈ -curv2
            @test lensRev[1].profile.a ≈ -a2
            @test lensRev[2].profile.curv ≈ -curv1
            @test lensRev[2].profile.a ≈ -a1

            @test_throws ErrorException lensASinglet(base, dir, curv1, ϵ1, a1, curv2, ϵ2, a2, thick, wl, riFunc, semiDiam; order="sideways")
        end

        @testset "lensEASinglet" begin
            curv1, curv2 = 0.02, -0.015
            ϵ1, ϵ2 = 0.5, 0.5
            a1 = [0.001, 0.0001]
            a2 = [0.0005, 0.00005]

            lensFwd = OpticTrace.lensEASinglet(base, dir, curv1, ϵ1, a1, curv2, ϵ2, a2, thick, wl, riFunc, semiDiam; order="forward")
            @test lensFwd[1].profile.curv == curv1
            @test lensFwd[1].profile.a == a1
            @test lensFwd[2].profile.curv == curv2
            @test lensFwd[2].profile.a == a2

            lensRev = OpticTrace.lensEASinglet(base, dir, curv1, ϵ1, a1, curv2, ϵ2, a2, thick, wl, riFunc, semiDiam; order="reverse")
            @test lensRev[1].profile.curv ≈ -curv2
            @test lensRev[1].profile.a ≈ -a2
            @test lensRev[2].profile.curv ≈ -curv1
            @test lensRev[2].profile.a ≈ -a1

            @test_throws ErrorException OpticTrace.lensEASinglet(base, dir, curv1, ϵ1, a1, curv2, ϵ2, a2, thick, wl, riFunc, semiDiam; order="sideways")
        end
    end

end
