@testset "lens_catalogs.jl" begin

  if !HAS_GLASS_CATALOG
    @info "Skipping lens_catalogs.jl: OpticTrace.dirBaseRefractiveIndex ($(OpticTrace.dirBaseRefractiveIndex)) not found on this machine -- every builder tested here looks up a real glass file internally, so none of this file's tests can run without it"
  else

    # These tests hit the real glass-catalog directory
    # (OpticTrace.dirBaseRefractiveIndex) rather than a stub riFunc, since
    # every builder here looks up its own glass internally. Expected
    # curvature/thickness/aspheric-coefficient literals are copied
    # straight from the corresponding src/lens_*.jl definitions (several
    # of them are function-local, not module constants, so can't be
    # referenced directly); expected refractive indices are computed via
    # getRefractiveIndexFunc against the same catalog files, cross-checked
    # where possible against the independently sourced SPECS.nd values
    # already used in test/refractive_index.jl (phase 8).

    base = ORIGIN
    dir = ZAXIS
    λ = 0.5875618 # d-line

    @testset "lens_edmund.jl" begin

        @testset "lens_EO38398 (not exported -- see TODO.md)" begin
            riN_SF11 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/schott/N-SF11.yml")

            lens = OpticTrace.lens_EO38398(base, dir, λ)
            @test length(lens) == 2
            @test lens[1].profile.curv == 8.496176720475789867E-02
            @test lens[2].profile.curv == 0.0
            @test lens[1].base.base == base
            @test lens[2].base.base ≈ base + 3.1 .* dir
            @test lens[1].aperture.semiDiameter == 5.75
            @test lens[1].mod.refIndexIn == refIndexDefault
            @test lens[1].mod.refIndexOut ≈ riN_SF11(λ)
            @test lens[2].mod.refIndexIn ≈ riN_SF11(λ)
            @test lens[2].mod.refIndexOut == refIndexDefault

            lensRev = OpticTrace.lens_EO38398(base, dir, λ; order = "reverse")
            @test length(lensRev) == 2
            @test lensRev[1].profile.curv == 0.0
            @test lensRev[2].profile.curv == -8.496176720475789867E-02
            @test lensRev[1].base.base == base
            @test lensRev[2].base.base ≈ base + 3.1 .* dir
            @test lensRev[1].mod.refIndexOut ≈ riN_SF11(λ)
            @test lensRev[2].mod.refIndexIn ≈ riN_SF11(λ)
        end

        @testset "lens_EO68001" begin
            riN_BK7 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/schott/N-BK7.yml")

            lens = lens_EO68001(base, dir, λ)
            @test lens[1].profile.curv == -3.869969040247679681E-02
            @test lens[2].profile.curv == 0.0
            @test lens[2].base.base ≈ base + 3.5 .* dir
            @test lens[1].aperture.semiDiameter == 12.5
            @test lens[1].mod.refIndexOut ≈ riN_BK7(λ)
            @test lens[2].mod.refIndexIn ≈ riN_BK7(λ)

            lensRev = lens_EO68001(base, dir, λ; order = "reverse")
            @test lensRev[1].profile.curv == 0.0
            @test lensRev[2].profile.curv == 3.869969040247679681E-02
            @test lensRev[2].base.base ≈ base + 3.5 .* dir
            @test lensRev[1].mod.refIndexOut ≈ riN_BK7(λ)
            @test lensRev[2].mod.refIndexIn ≈ riN_BK7(λ)
        end

        @testset "lens_EO67548" begin
            riN_BK7 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/schott/N-BK7.yml")

            lens = lens_EO67548(base, dir, λ)
            @test lens[1].profile.curv == 2.579979360165119903E-02
            @test lens[2].base.base ≈ base + 4.5 .* dir
            @test lens[1].aperture.semiDiameter == 12.5
            @test lens[1].mod.refIndexOut ≈ riN_BK7(λ)

            lensRev = lens_EO67548(base, dir, λ; order = "reverse")
            @test lensRev[1].profile.curv == 0.0
            @test lensRev[2].profile.curv == -2.579979360165119903E-02
            @test lensRev[2].base.base ≈ base + 4.5 .* dir
            @test lensRev[1].mod.refIndexOut ≈ riN_BK7(λ)
            @test lensRev[2].mod.refIndexIn ≈ riN_BK7(λ)
        end

        @testset "lens_EO67652" begin
            riN_BK7 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/schott/N-BK7.yml")

            lens = lens_EO67652(base, dir, λ)
            @test lens[1].profile.curv == 1.304461257500649958E-02
            @test lens[2].profile.curv == -1.304461257500649958E-02
            @test lens[2].base.base ≈ base + 3.5 .* dir
            @test lens[1].mod.refIndexOut ≈ riN_BK7(λ)

            lensRev = lens_EO67652(base, dir, λ; order = "reverse")
            @test lensRev[1].profile.curv == 1.304461257500649958E-02
            @test lensRev[2].profile.curv == -1.304461257500649958E-02
            @test lensRev[2].base.base ≈ base + 3.5 .* dir
            @test lensRev[1].mod.refIndexOut ≈ riN_BK7(λ)
        end
    end

    @testset "lens_thorlabs.jl" begin

        @testset "lensAC508180AB" begin
            riN_LAK22 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/schott/N-LAK22.yml")
            riN_SF6 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/schott/N-SF6.yml")

            lensFwd = lensAC508180AB(base, dir, λ)
            @test length(lensFwd) == 3
            @test lensFwd[1].profile.curv == 6.923078165548776121E-03
            @test lensFwd[2].profile.curv == -8.663404145164059142E-03
            @test lensFwd[3].profile.curv == -3.046912099599362322E-03
            @test lensFwd[1].mod.refIndexIn == refIndexDefault
            @test lensFwd[1].mod.refIndexOut ≈ riN_LAK22(λ)
            @test lensFwd[2].mod.refIndexIn ≈ riN_LAK22(λ)
            @test lensFwd[2].mod.refIndexOut ≈ riN_SF6(λ)
            @test lensFwd[3].mod.refIndexOut == refIndexDefault
            @test lensFwd[2].base.base ≈ base + 9.5 .* dir
            @test lensFwd[3].base.base ≈ base + 13.5 .* dir
            @test lensFwd[1].coating.type == "testcoat"

            lensRev = lensAC508180AB(base, dir, λ; order = "reverse")
            @test lensRev[1].profile.curv ≈ 3.046912099599362322E-03
            @test lensRev[2].profile.curv ≈ 8.663404145164059142E-03
            @test lensRev[3].profile.curv ≈ -6.923078165548776121E-03
            @test lensRev[1].mod.refIndexOut ≈ riN_SF6(λ)
            @test lensRev[2].mod.refIndexOut ≈ riN_LAK22(λ)
            @test lensRev[3].mod.refIndexOut == refIndexDefault
            @test lensRev[2].base.base ≈ base + 4.0 .* dir
            @test lensRev[3].base.base ≈ base + 13.5 .* dir

            @test_throws ErrorException lensAC508180AB(base, dir, λ; order = "sideways")
        end

        @testset "lensAC127050A" begin
            riN_BK7 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/schott/N-BK7.yml")
            riN_SF2 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/schott/N-SF2.yml")

            lensFwd = lensAC127050A(base, dir, λ)
            @test lensFwd[1].profile.curv == OpticTrace.curvAC127050A_1
            @test lensFwd[2].profile.curv == OpticTrace.curvAC127050A_2
            @test lensFwd[3].profile.curv == OpticTrace.curvAC127050A_3
            @test lensFwd[1].mod.refIndexOut ≈ riN_BK7(λ)
            @test lensFwd[2].mod.refIndexOut ≈ riN_SF2(λ)
            @test lensFwd[2].base.base ≈ base + OpticTrace.thickAC127050A_1 .* dir

            lensRev = lensAC127050A(base, dir, λ; order = "reverse")
            @test lensRev[1].profile.curv ≈ -OpticTrace.curvAC127050A_3
            @test lensRev[3].profile.curv ≈ -OpticTrace.curvAC127050A_1

            @test_throws ErrorException lensAC127050A(base, dir, λ; order = "sideways")
        end

        @testset "lensAC127019AB (not exported, no order validation -- see TODO.md)" begin
            riN_LAK10 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/schott/N-LAK10.yml")
            riN_SF57 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/schott/N-SF57.yml")

            lensFwd = OpticTrace.lensAC127019AB(base, dir, λ)
            @test lensFwd[1].profile.curv == OpticTrace.curvAC127019AB_1
            @test lensFwd[2].profile.curv == OpticTrace.curvAC127019AB_2
            @test lensFwd[1].mod.refIndexOut ≈ riN_LAK10(λ)
            @test lensFwd[2].mod.refIndexOut ≈ riN_SF57(λ)

            lensRev = OpticTrace.lensAC127019AB(base, dir, λ; order = "reverse")
            @test lensRev[1].profile.curv ≈ -OpticTrace.curvAC127019AB_3

            # any order value other than exactly "forward" is silently
            # treated as "reverse" -- confirmed here rather than erroring
            lensGarbage = OpticTrace.lensAC127019AB(base, dir, λ; order = "sideways")
            @test lensGarbage[1].profile.curv == lensRev[1].profile.curv
            @test lensGarbage[2].profile.curv == lensRev[2].profile.curv
            @test lensGarbage[3].profile.curv == lensRev[3].profile.curv
        end

        @testset "lens_ACL12708U" begin
            # not a wavelength despite the parameter name `wl` -- used
            # directly as the refractive index between the two surfaces,
            # see TODO.md / this function's docstring
            wl = 1.5

            lens = lens_ACL12708U(base, dir, wl)
            @test length(lens) == 2
            @test lens[1].profile isa OpticTrace.SurfProfileConic
            @test lens[1].profile.curv == -OpticTrace.osurf2C_ACL12708U
            @test lens[1].mod.refIndexIn == refIndexDefault
            @test lens[1].mod.refIndexOut == wl
            @test lens[1].aperture.semiDiameter == OpticTrace.semiDiamACL12708U

            @test lens[2].profile isa OpticTrace.SurfProfileEvenAsphere
            @test lens[2].profile.curv == -OpticTrace.osurf1C_ACL12708U
            @test lens[2].profile.ϵ == OpticTrace.osurf1ϵ_ACL12708U
            @test lens[2].profile.a == .-OpticTrace.osurf1ASP_ACL12708U
            @test lens[2].mod.refIndexIn == wl
            @test lens[2].mod.refIndexOut == refIndexDefault
            @test lens[2].base.base ≈ base + OpticTrace.thickACL12708U .* dir
        end

        @testset "lens_TLF220APC" begin
            riD_ZK3 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/cdgm/D-ZK3.yml")

            lens = lens_TLF220APC(base, dir, λ)
            @test length(lens) == 2
            @test lens[1].profile.curv == 0.0
            @test lens[2].profile.curv == OpticTrace.osurf2C_F220APC
            @test lens[2].profile.ϵ == OpticTrace.osurf2ϵ_F220APC
            @test lens[2].profile.a == OpticTrace.osurf2ASP_F220APC
            @test lens[1].mod.refIndexOut ≈ riD_ZK3(λ)
            @test lens[2].base.base ≈ base + OpticTrace.thickF220APC .* dir
            @test lens[1].aperture.semiDiameter == OpticTrace.semiDiamF220APC

            @test_throws ErrorException lens_TLF220APC(base, dir, λ; order = "sideways")
        end

        @testset "lens_TLF357775_405" begin
            riD_LAK6 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/cdgm/D-LAK6.yml")

            lens = lens_TLF357775_405(base, dir, λ)
            @test length(lens) == 2
            @test lens[1].profile.curv == 0.0
            @test lens[2].profile.curv == -3.497173380949716859E-01
            @test lens[1].mod.refIndexOut ≈ riD_LAK6(λ)
            @test lens[2].base.base ≈ base + thickTL357775_405 .* dir

            @test_throws ErrorException lens_TLF357775_405(base, dir, λ; order = "sideways")
        end

        @testset "lens_TLAC254_060" begin
            riBAF11 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/hoya/BAF11.yml")
            riEFD10 = getRefractiveIndexFunc(OpticTrace.dirBaseRefractiveIndex, "glass/hoya/E-FD10.yml")

            # cross-check against the independently sourced SPECS.nd
            # values already used in test/refractive_index.jl (phase 8)
            @test riBAF11(λ) ≈ 1.66672 atol = 1e-4
            @test riEFD10(λ) ≈ 1.72825 atol = 1e-4

            lensFwd = lens_TLAC254_060(base, dir, λ)
            @test length(lensFwd) == 3
            @test lensFwd[1].profile.curv == 2.398656752218759900E-002
            @test lensFwd[2].profile.curv == -3.863987635239570000E-002
            @test lensFwd[3].profile.curv == -4.334633723450400300E-003
            @test lensFwd[1].mod.refIndexOut ≈ riBAF11(λ)
            @test lensFwd[2].mod.refIndexOut ≈ riEFD10(λ)
            @test lensFwd[2].base.base ≈ base + 8.0 .* dir
            @test lensFwd[3].base.base ≈ base + 10.5 .* dir

            lensRev = lens_TLAC254_060(base, dir, λ; order = "reverse")
            @test lensRev[1].profile.curv ≈ 4.334633723450400300E-003
            @test lensRev[2].profile.curv ≈ 3.863987635239570000E-002
            @test lensRev[3].profile.curv ≈ -2.398656752218759900E-002

            @test_throws ErrorException lens_TLAC254_060(base, dir, λ; order = "sideways")
        end
    end

  end # if HAS_GLASS_CATALOG

end
