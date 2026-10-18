# Wavefront analysis (src/wavefront.jl) -- see plans/wavefront_analysis.md
using Test, LinearAlgebra, GeometryBasics
using OpticTrace
const OT = OpticTrace

# Concave paraboloid mirror (focal length f, vertex at p0) in a medium of index n,
# collimated light travelling along +z, image plane `shift` mm beyond the focus.
# Optional rigid transform: rotation matrix Rm and translation t applied to everything.
function parabolaSystem(f; p0 = [0., 0., 100.], shift = 0.0, n = 1.0, semiDiam = 10.0,
                        Rm = Matrix{Float64}(I, 3, 3), t = zeros(3))
    place(p) = Point3(Rm * p .+ t)
    dirv(d) = Vec3(Rm * d)
    ydir = dirv([0., 1., 0.])
    mirror = reflectConic("mirror", place(p0), dirv([0., 0., 1.]), n, n, -1 / (2f), 0.0, semiDiam, "c"; ydir)
    img = refractSphere("img", place(p0 .+ [0., 0., -f + shift]), dirv([0., 0., 1.]), n, n, 0.0, 25.0, "c"; ydir)
    Vector{OT.AbstractSurface}([mirror, img]), CollimatedField(dirv([0., 0., 1.]))
end

weightedRMS(r) = (m = r.mask; o = r.opd[m]; sqrt(sum(abs2, o .- sum(o) / length(o)) / length(o)))

@testset "wavefront" begin
    mk(name, z, nin, nout, c = 0.0) = refractSphere(name, Point3(0., 0., z), Vec3(0., 0., 1.), nin, nout, c, 25., "c")

    @testset "pathLength" begin
        # start on surface 1 (local frame of geo[1]): 0 + 10 mm of glass n=1.5 + 30 mm of air
        geo = Vector{OT.AbstractSurface}([mk("s1", 10., 1., 1.5), mk("s2", 20., 1.5, 1.), mk("img", 50., 1., 1.)])
        st, trc = traceGeometryRel(OT.Ray(Point3(0., 0., 0.), Vec3(0., 0., 1.)), geo)
        @test st == 0
        @test pathLength(trc, geo, 0.5) ≈ 45.0
        @test pathLength(trc, geo, 0.5; upto = 3) ≈ 15.0
        # a sphere start behind the vertex: delta is negative and accepted
        geo3 = Vector{OT.AbstractSurface}([mk("sph", 10., 1., 1., 0.1)])
        st, trc = traceGeometryRel(OT.Ray(Point3(0., 2., 30.), Vec3(0., 0., 1.)), geo3)
        @test st == 0
        @test trc[end].delta < 0
    end

    @testset "perfect collimated system has zero OPD" begin
        geo, field = parabolaSystem(100.0)
        r = wavefrontMap(geo, field, 0.5; N = 17, geoWavelength = 0.5)
        @test r isa WavefrontResult
        @test count(r.mask) == count(uv -> uv[1]^2 + uv[2]^2 <= 1 + 1e-12, Iterators.product(r.u, r.v))
        @test maximum(abs, r.opd[r.mask]) < 1e-6
        @test isnan(r.opd[1, 1])   # corner outside the pupil disk
        @test r.opd[9, 9] == 0     # chief ray: OPD is exactly zero at the pupil center
        @test r.refRadius ≈ 100.0 rtol = 1e-3
        @test r.refSide == 1
        @test all(r.weight[r.mask] .== 1.0)
    end

    @testset "perfect finite-conjugate system has zero OPD" begin
        # point source at the center of curvature of a concave spherical mirror
        R = 100.0
        mirror = reflectConic("mirror", Point3(0., 0., R), Vec3(0., 0., 1.), 1., 1., -1 / R, 1.0, 10., "c")
        img = refractSphere("img", Point3(0., 0., 0.), Vec3(0., 0., 1.), 1., 1., 0.0, 25., "c")
        geo = Vector{OT.AbstractSurface}([mirror, img])
        r = wavefrontMap(geo, ObjectPoint(Point3(0., 0., 0.)), 0.5; N = 17, geoWavelength = 0.5)
        @test maximum(abs, r.opd[r.mask]) < 1e-6
    end

    @testset "defocus: RMS = peak/(2*sqrt(3)), sign independent" begin
        f, δ, λ = 100.0, 0.5, 0.5
        geo, field = parabolaSystem(f; shift = δ)
        r = wavefrontMap(geo, field, λ; N = 65, geoWavelength = λ)
        h = 10.0
        θ = atan(h, f - h^2 / (4f))
        peak = δ * (1 - cos(θ)) / (λ * LENGTH_TO_WAVELENGTH)   # waves
        @test weightedRMS(r) ≈ peak / (2sqrt(3)) rtol = 0.03
        # flipping the shift flips the sign of the OPD, not its size
        geo2, field2 = parabolaSystem(f; shift = -δ)
        r2 = wavefrontMap(geo2, field2, λ; N = 65, geoWavelength = λ)
        @test weightedRMS(r2) ≈ weightedRMS(r) rtol = 0.02
        @test sign(r.opd[33, 65]) == -sign(r2.opd[33, 65])
    end

    @testset "image-space index is used" begin
        f, δ, λ = 100.0, 0.5, 0.5
        r1 = wavefrontMap(parabolaSystem(f; shift = δ)[1], CollimatedField(Vec3(0., 0., 1.)), λ; N = 33, geoWavelength = λ)
        geo, field = parabolaSystem(f; shift = δ, n = 1.5)
        r15 = wavefrontMap(geo, field, λ; N = 33, geoWavelength = λ)
        @test r15.nImage == 1.5
        @test weightedRMS(r15) ≈ 1.5 * weightedRMS(r1) rtol = 0.02
    end

    @testset "rigid motion of the whole system leaves the OPD unchanged" begin
        f, δ, λ = 100.0, 0.5, 0.5
        r0 = wavefrontMap(parabolaSystem(f; shift = δ)...,  λ; N = 17, geoWavelength = λ) |> x -> x
        a, b = 0.3, -0.5
        Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)]
        Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)]
        geo, field = parabolaSystem(f; shift = δ, Rm = Ry * Rx, t = [5., -3., 20.])
        r1 = wavefrontMap(geo, field, λ; N = 17, geoWavelength = λ)
        @test r1.mask == r0.mask
        @test maximum(abs, (r1.opd .- r0.opd)[r0.mask]) < 1e-6
    end

    @testset "pupil at infinity gives a plane reference" begin
        # stop at the focal plane of a paraboloid: rays through the stop center
        # leave the mirror parallel, so the exit pupil is at infinity
        stop = OT.referencePlane("stop", Point3(0., 0., 0.), Vec3(0., 0., 1.), 1., 10., "c")
        mirror = reflectConic("mirror", Point3(0., 0., 100.), Vec3(0., 0., 1.), 1., 1., -1 / 200, 0.0, 30., "c")
        img = refractSphere("img", Point3(0., 0., -200.), Vec3(0., 0., 1.), 1., 1., 0.0, 25., "c")
        geo = Vector{OT.AbstractSurface}([stop, mirror, img])
        r = wavefrontMap(geo, ObjectPoint(Point3(0., 0., -50.)), 0.5; N = 9, geoWavelength = 0.5)
        @test r isa WavefrontResult
        @test isinf(r.refRadius)
        @test all(isnan, r.exitPupil)
        @test all(isfinite, r.opd[r.mask])
        @test r.opd[5, 5] == 0
    end

    @testset "plane reference is the large-radius limit of the sphere" begin
        C = Point3(0., 0., 0.)
        b = Point3(3., -2., 80.)
        d = Vec3(normalize(C .- b))
        for dd in (d, Vec3(normalize(Vec3(0.02, 0.01, -1.0))))
            plane = OT.oplToReference(b, dd, 5.0, 1.0, C, Inf, 1)
            sphere = OT.oplToReference(b, dd, 5.0, 1.0, C, 1e9, 1)
            # the constant sigma*R*n is dropped, so the large sphere and the plane agree
            @test sphere ≈ plane atol = 1e-6
        end
        # sphere miss
        @test isnan(OT.oplToReference(Point3(1000., 0., 0.), Vec3(0., 0., 1.), 0.0, 1.0, C, 10.0, 1))
        # plane reference point
        @test OT.referencePoint(b, d, C, Inf, 1, Point3(0., 0., 50.), Vec3(0., 0., -1.)) ≈ Point3(1.875, -1.25, 50.) atol = 1e-9
    end

    @testset "apodization weights" begin
        geo, field = parabolaSystem(100.0)
        r = wavefrontMap(geo, field, 0.5; N = 9, geoWavelength = 0.5, apodization = GaussianPupil(0.5))
        @test r.weight[5, 5] == 1.0
        @test r.weight[5, 9] ≈ exp(-1 / 0.25)
        r = wavefrontMap(geo, field, 0.5; N = 9, geoWavelength = 0.5, apodization = FunctionPupil((u, v) -> 2u + 3))
        @test r.weight[9, 5] ≈ 5.0
        r = wavefrontMap(geo, field, 0.5; N = 9, geoWavelength = 0.5, apodization = AnglePupil(d -> d[3]))
        @test r.weight[5, 5] ≈ 1.0
        @test OT.pupilAmplitude(UniformPupil(), 0.3, 0.1, nothing) == 1.0
    end

    @testset "phase hook adds optical path" begin
        # a pass-through surface with a user-defined phase: kappa mm of path per mm of local x
        struct PhaseTilt{T} <: OT.AbstractBendDielectric{T}
            refIndexIn::T
            refIndexOut::T
            kappa::T
        end
        OT.bendPhase(m::PhaseTilt, s, rayIn, trcOut, λ) = m.kappa * s.toLocalCoord(trcOut.ray.base)[1]
        OT.modFunc(ray::OT.Ray{3,T}, normal::Vec3{T}, d::PhaseTilt{T}) where T<:Real =
            OT.modFunc(ray, normal, OT.DielectricT(d.refIndexIn, d.refIndexOut))
        s = mk("stop", 90., 1., 1.)
        s.aperture = OT.SizeLens(10.0)
        s = OT.OptSurface(s.surfname, s.base, s.aperture, s.profile, PhaseTilt(1.0, 1.0, 1e-4), s.coating,
                          s.toGlobalCoord, s.toLocalCoord, s.toGlobalDir, s.toLocalDir, s.color)
        pgeo, field = parabolaSystem(100.0)
        geo = Vector{OT.AbstractSurface}([s; pgeo])
        λ = 0.5
        r = wavefrontMap(geo, field, λ; N = 9, geoWavelength = λ)
        # a perfect paraboloid after the plate: the OPD is just the added phase,
        # kappa * x(mm) = kappa * 10 * (u or v) mm of path, in waves
        peak = 1e-4 * 10.0 / (λ * LENGTH_TO_WAVELENGTH)
        @test maximum(abs, r.opd[r.mask]) ≈ peak rtol = 1e-6
        @test r.opd[5, 5] == 0
        @test all(abs.(r.opd[r.mask] .+ reverse(r.opd)[r.mask]) .< 1e-6)   # odd symmetry
        @test count(i -> abs(r.opd[i]) > 1e-9, findall(r.mask)) == count(r.mask) - 9   # zero along the axis perpendicular to the gradient
    end

    @testset "unsupported / warned cases" begin
        geo, field = parabolaSystem(100.0)
        diff = cDiffuser("diff", Point3(0., 0., 50.), Vec3(0., 0., 1.), 1., 1., 0.1, 10., "c")
        gd = Vector{OT.AbstractSurface}([geo[1], diff, geo[2]])
        @test_logs (:warn, r"CDiffuser") @test isnothing(wavefrontMap(gd, field, 0.5; N = 5, geoWavelength = 0.5))
        @test_logs (:warn, r"geoWavelength|built for") wavefrontMap(geo, field, 0.5; N = 5, geoWavelength = 0.6)
        @test_throws ArgumentError wavefrontMap(geo, field, 0.5; N = 8)
        @test_throws ArgumentError wavefrontMap(geo, field, 0.5; aimer = StopCenterAim(2))
        @test_throws ArgumentError wavefrontMap(geo, field, 0.5; aimer = StopCenterAim("nope"))
        # an image surface that is not reachable (chief ray misses): warns and returns nothing
        far = Vector{OT.AbstractSurface}([geo[1], refractSphere("img", Point3(0., 0., 500.), Vec3(1., 0., 0.), 1., 1., 0.0, 25., "c")])
        @test_logs (:warn, r"chief ray") @test isnothing(wavefrontMap(far, field, 0.5; N = 5, geoWavelength = 0.5))
    end

    @testset "StopCenterAim" begin
        geo, _ = parabolaSystem(100.0)
        aimer = StopCenterAim("mirror")
        @test OT.pupilSurfaceIndex(aimer, geo) == 1
        @test OT.pupilPoint(aimer, geo, 0.5, 0.0) ≈ Point3(5., 0., 100.)
        st, ray = aimRay(aimer, geo, ObjectPoint(Point3(0., 0., 0.)), 0.0, 0.0)
        @test st == 0 && ray.dir ≈ Vec3(0., 0., 1.)
        st, ray = aimRay(aimer, geo, CollimatedField(Vec3(0., 0., 2.)), 1.0, 0.0)
        @test ray.dir ≈ Vec3(0., 0., 1.) && ray.base ≈ Point3(10., 0., 100.)
    end
end

@testset "wavefront statistics" begin
    λ = 0.5
    # spherical mirror (f = 100, R = 200) with collimated light: pure spherical aberration
    function sphereMirror(shift)
        mirror = reflectConic("mirror", Point3(0., 0., 100.), Vec3(0., 0., 1.), 1., 1., -1 / 200, 1.0, 10., "c")
        img = refractSphere("img", Point3(0., 0., 0.0 + shift), Vec3(0., 0., 1.), 1., 1., 0.0, 25., "c")
        Vector{OT.AbstractSurface}([mirror, img]), CollimatedField(Vec3(0., 0., 1.))
    end

    @testset "perfect system" begin
        geo, field = parabolaSystem(100.0)
        r = wavefrontMap(geo, field, λ; N = 21, geoWavelength = λ)
        st = wavefrontStats(r)
        @test st.rms < 1e-6 && st.pv < 1e-5
        @test st.strehl ≈ 1 atol = 1e-9
        @test st.strehlMarechal ≈ 1 atol = 1e-9
        @test st.npoints == count(r.mask)
        @test wavefrontRMS(r; weighted = false) < 1e-6
    end

    @testset "refocusOPD and best focus" begin
        geo, field = parabolaSystem(100.0; shift = 0.5)
        r = wavefrontMap(geo, field, λ; N = 33, geoWavelength = λ)
        @test refocusOPD(r, 0.0) ≈ r.opd nans = true
        # the image plane (z = 0.5) is 0.5 mm before the focus (z = 0) along the reflected ray (-z);
        # moving the reference point 0.5 mm along the propagation direction restores a perfect focus
        @test wavefrontRMS(r; opd = refocusOPD(r, 0.5)) < 1e-5
        δ, rms, opd = bestFocus(r)
        @test δ ≈ 0.5 atol = 1e-3
        @test rms < 1e-5
        @test opd ≈ refocusOPD(r, δ) nans = true
        @test bestFocus(r; range = 2.0)[1] ≈ 0.5 atol = 1e-3
    end

    @testset "spherical aberration: best focus cuts the RMS by 4" begin
        geo, field = sphereMirror(0.0)
        r = wavefrontMap(geo, field, λ; N = 65, geoWavelength = λ)
        rms0 = wavefrontRMS(r)
        δ, rmsBest, _ = bestFocus(r)
        # W040*rho^4 has RMS 2/(3 sqrt(5)) W040 at the paraxial focus, W040/(6 sqrt(5)) at best focus
        @test rmsBest / rms0 ≈ 0.25 rtol = 0.05
        # the best-focus defocus has peak size W040 = rms0 * 3 sqrt(5)/2 (in waves)
        W040 = rms0 * 3sqrt(5) / 2 * λ * LENGTH_TO_WAVELENGTH   # mm
        θ = 2asin(10 / 200)   # marginal ray angle after reflection
        @test abs(δ) ≈ 2W040 / θ^2 rtol = 0.15
        # RMS at the best focus really is a minimum
        for d in (δ - 0.05, δ + 0.05)
            @test wavefrontRMS(r; opd = refocusOPD(r, d)) > rmsBest
        end
    end

    @testset "Strehl: Marechal vs reference-point Strehl, PV, weights" begin
        geo, field = parabolaSystem(100.0; shift = 0.02)
        r = wavefrontMap(geo, field, λ; N = 65, geoWavelength = λ)
        rms = wavefrontRMS(r)
        @test 0.03 < rms < 0.1
        @test strehlRatio(r) ≈ strehlMarechal(rms) rtol = 0.03
        @test wavefrontPV(r) > rms
        # a pure defocus a*rho^2 has RMS a/(2 sqrt(3)) and PV a
        @test wavefrontPV(r) ≈ 2sqrt(3) * rms rtol = 0.05
        # apodization changes the weighted RMS but not the unweighted one
        ra = wavefrontMap(geo, field, λ; N = 65, geoWavelength = λ, apodization = GaussianPupil(0.5))
        @test wavefrontRMS(ra; weighted = false) ≈ rms rtol = 1e-9
        @test wavefrontRMS(ra) < rms
        # leaving out the low-weight edge reduces PV
        @test wavefrontPV(ra; minWeight = 0.5) < wavefrontPV(ra)
        # empty
        empty = WavefrontResult(r.u, r.v, fill(NaN, 65, 65), fill(false, 65, 65), zeros(65, 65), r.pupilPoints, r.lastBase,
            r.lastDir, r.oplToLast, r.chiefBase, r.chiefDir, r.chiefOPL, r.imagePoint, r.exitPupil, r.refRadius, r.refSide,
            r.nImage, r.λ, r.geoWavelength, r.field)
        @test isnan(wavefrontRMS(empty)) && isnan(wavefrontPV(empty)) && isnan(strehlRatio(empty))
    end

    @testset "legacy opdRel cross-check agrees inside its envelope" begin
        # stop, plano-convex lens, image plane; on-axis point 1000 mm in front of the stop
        stop = OT.referencePlane("stop", Point3(0., 0., 0.), Vec3(0., 0., 1.), 1., 8., "c")
        l1 = refractSphere("l1", Point3(0., 0., 5.), Vec3(0., 0., 1.), 1., 1.5, 1 / 50., 25., "c")
        l2 = refractSphere("l2", Point3(0., 0., 10.), Vec3(0., 0., 1.), 1.5, 1., 0.0, 25., "c")
        img = refractSphere("img", Point3(0., 0., 120.), Vec3(0., 0., 1.), 1., 1., 0.0, 25., "c")
        geo = Vector{OT.AbstractSurface}([stop, l1, l2, img])
        r = wavefrontMap(geo, ObjectPoint(Point3(0., 0., -1000.)), λ; N = 65, geoWavelength = λ)
        rmsNew = wavefrontRMS(r)
        rmsOld, opds = rmsWavefrontLegacy(geo, Point3(0., 0., -1000.), atan(8, 1000), λ)
        @test rmsNew > 0.3     # aberrated enough for the comparison to mean something
        @test rmsOld ≈ rmsNew rtol = 0.05
        @test length(opds) == 4 * 12^2
    end
end

@testset "wavefront PSF and spot" begin
    λ = 0.5
    λp = λ * LENGTH_TO_WAVELENGTH   # mm
    function naOf(r)
        e1, e2 = OT.defaultImageAxes(r.chiefDir)
        maximum(hypot(sum(r.lastDir[i] .* e1), sum(r.lastDir[i] .* e2)) for i in eachindex(r.opd) if r.mask[i])
    end

    @testset "diffraction-limited PSF is an Airy pattern" begin
        geo, field = parabolaSystem(100.0)
        r = wavefrontMap(geo, field, λ; N = 65, geoWavelength = λ)
        NA = naOf(r)
        psf = wavefrontPSF(r)
        @test psf isa PSFResult
        @test size(psf.intensity) == (129, 129)
        @test psf.strehl ≈ 1 atol = 1e-6
        @test psf.intensity[65, 65] ≈ 1 atol = 1e-6
        @test psf.pupilDistortion < 0.01
        @test psf.pitch ≈ λp / (8NA)
        m = spotMetrics(psf; fractions = (0.838, 0.9))
        @test hypot(m.peakX, m.peakY) < psf.pitch
        @test m.fwhmX ≈ 1.029 * λp / (2NA) rtol = 0.03
        @test m.fwhmY ≈ m.fwhmX rtol = 0.02
        @test m.fwhmArea ≈ 1.029 * λp / (2NA) rtol = 0.05
        @test m.ee[1] ≈ 0.61 * λp / NA rtol = 0.04   # 83.8% of the energy is inside the first dark ring
        @test m.rmsRadius > 0 && m.semiMajor >= m.semiMinor
        @test m.angle isa Float64
    end

    @testset "defocused PSF: Strehl falls, spot matches the geometric blur" begin
        geo, field = parabolaSystem(100.0; shift = 0.5)
        r = wavefrontMap(geo, field, λ; N = 65, geoWavelength = λ)
        grms, offs = geometricSpot(r)
        @test length(offs) == count(r.mask)
        # blur radius of a uniformly filled disc: 1/sqrt(2) of its edge radius
        edge = 0.5 * tan(2asin(10 / 200)) * 0.99
        @test grms ≈ edge / sqrt(2) rtol = 0.1
        psf = wavefrontPSF(r; pixels = 257, halfWidth = 0.07)
        @test psf.strehl < 0.1
        m = spotMetrics(psf; center = :chief)
        @test m.rmsRadius ≈ grms rtol = 0.15
        # the reference-point Strehl is the on-axis PSF value
        @test psf.intensity[129, 129] ≈ strehlRatio(r) rtol = 1e-6
        # energy is (almost) all inside the window
        @test sum(psf.intensity) * psf.pitch^2 / psf.totalEnergy ≈ 1 rtol = 0.05
        # refocusing restores the diffraction-limited peak
        psf2 = wavefrontPSF(r; opd = refocusOPD(r, 0.5))
        @test psf2.strehl ≈ 1 atol = 1e-4
    end

    @testset "a wavefront tilt moves the peak, not the Strehl" begin
        struct PhaseTiltPSF{T} <: OT.AbstractBendDielectric{T}
            refIndexIn::T
            refIndexOut::T
            kappa::T
        end
        OT.bendPhase(m::PhaseTiltPSF, s, rayIn, trcOut, λ) = m.kappa * s.toLocalCoord(trcOut.ray.base)[1]
        OT.modFunc(ray::OT.Ray{3,T}, normal::Vec3{T}, d::PhaseTiltPSF{T}) where T<:Real =
            OT.modFunc(ray, normal, OT.DielectricT(d.refIndexIn, d.refIndexOut))
        s0 = refractSphere("stop", Point3(0., 0., 90.), Vec3(0., 0., 1.), 1., 1., 0.0, 25., "c")
        s0.aperture = OT.SizeLens(10.0)
        s = OT.OptSurface(s0.surfname, s0.base, s0.aperture, s0.profile, PhaseTiltPSF(1.0, 1.0, 1e-4), s0.coating,
                          s0.toGlobalCoord, s0.toLocalCoord, s0.toGlobalDir, s0.toLocalDir, s0.color)
        pgeo, field = parabolaSystem(100.0)
        r = wavefrontMap(Vector{OT.AbstractSurface}([s; pgeo]), field, λ; N = 65, geoWavelength = λ)
        @test strehlRatio(r) < 0.1     # at the reference point
        psf = wavefrontPSF(r)
        @test psf.strehl > 0.95        # the peak is displaced, not reduced
        m = spotMetrics(psf)
        @test hypot(m.peakX, m.peakY) ≈ 100 * 1e-4 rtol = 0.1   # focal length * tilt
    end

    @testset "PSF edge cases and warnings" begin
        geo, field = parabolaSystem(100.0)
        r = wavefrontMap(geo, field, λ; N = 33, geoWavelength = λ)
        @test_throws ArgumentError wavefrontPSF(r; pixels = 128)
        @test_throws ArgumentError wavefrontPSF(r; pitch = 1e-3, halfWidth = 0.1)
        @test_logs (:warn, r"Nyquist") wavefrontPSF(r; pitch = 5e-3, pixels = 9)
        @test_logs (:warn, r"aliasing") match_mode = :any wavefrontPSF(r; halfWidth = 5.0, pixels = 9)
        none = fill(NaN, size(r.opd))
        @test_throws ArgumentError wavefrontPSF(r; opd = none)
        @test_throws ArgumentError spotMetrics(wavefrontPSF(r; pixels = 9); center = :nope)
        axes = (Vec3(0., 1., 0.), Vec3(1., 0., 0.))
        p = wavefrontPSF(r; pixels = 9, axes)
        @test p.axisX == Vec3(0., 1., 0.)
        # a window too small for the energy fractions requested gives NaN
        @test isnan(spotMetrics(wavefrontPSF(r; pixels = 9); fractions = (0.99,)).ee[1])
        # the telecentric (plane reference) case: PSF axes need only the chief direction
        @test OT.defaultImageAxes(Vec3(1., 0., 0.))[1] ≈ Vec3(0., 1., 0.)
    end
end

# a collimated beam through a pass-through plate with an x-dependent phase, then a perfect paraboloid
struct PhaseTiltPlot{T} <: OT.AbstractBendDielectric{T}
    refIndexIn::T
    refIndexOut::T
    kappa::T
end
OT.bendPhase(m::PhaseTiltPlot, s, rayIn, trcOut, λ) = m.kappa * s.toLocalCoord(trcOut.ray.base)[1]
OT.modFunc(ray::OT.Ray{3,T}, normal::Vec3{T}, d::PhaseTiltPlot{T}) where T<:Real =
    OT.modFunc(ray, normal, OT.DielectricT(d.refIndexIn, d.refIndexOut))
function tiltedSystem(kappa)
    s0 = refractSphere("stop", Point3(0., 0., 90.), Vec3(0., 0., 1.), 1., 1., 0.0, 25., "c")
    s0.aperture = OT.SizeLens(10.0)
    s = OT.OptSurface(s0.surfname, s0.base, s0.aperture, s0.profile, PhaseTiltPlot(1.0, 1.0, kappa), s0.coating,
                      s0.toGlobalCoord, s0.toLocalCoord, s0.toGlobalDir, s0.toLocalDir, s0.color)
    pgeo, field = parabolaSystem(100.0)
    Vector{OT.AbstractSurface}([s; pgeo]), field
end

@testset "wavefrontRemove" begin
    λ = 0.5
    geo, field = tiltedSystem(1e-4)
    r = wavefrontMap(geo, field, λ; N = 33, geoWavelength = λ)
    @test wavefrontRMS(r) > 0.5                       # a 4-wave tilt across the pupil
    @test wavefrontRemove(r; remove = :none) == r.opd || isequal(wavefrontRemove(r; remove = :none), r.opd)
    p = wavefrontRemove(r; remove = :piston)
    @test wavefrontRMS(r; opd = p) ≈ wavefrontRMS(r)  # piston removal does not change the RMS
    t = wavefrontRemove(r; remove = :pistonTilt)
    @test wavefrontRMS(r; opd = t) < 1e-6
    @test strehlRatio(r; opd = t) ≈ 1 atol = 1e-9
    @test wavefrontRMS(r; opd = wavefrontRemove(r; remove = :pistonTiltFocus)) < 1e-6
    @test isnan(t[1, 1]) && !isnan(t[17, 17])
    # defocus is removed by the focus term but not by tilt
    gd, fd = parabolaSystem(100.0; shift = 0.5)
    rd = wavefrontMap(gd, fd, λ; N = 33, geoWavelength = λ)
    @test wavefrontRMS(rd; opd = wavefrontRemove(rd; remove = :pistonTilt)) ≈ wavefrontRMS(rd) rtol = 1e-3
    @test wavefrontRMS(rd; opd = wavefrontRemove(rd; remove = :pistonTiltFocus)) < 0.01 * wavefrontRMS(rd)   # defocus is rho^2 to < 1%
    @test_throws ArgumentError wavefrontRemove(r; remove = :nope)
    none = fill(NaN, size(r.opd))
    @test all(isnan, wavefrontRemove(r; opd = none))
end

@testset "wavefront plots" begin
    GLMakie.activate!(visible = false)
    λ = 0.5
    geo, field = tiltedSystem(5e-5)
    r = wavefrontMap(geo, field, λ; N = 33, geoWavelength = λ)
    mid = 17

    @testset "wavefront heatmap, surface, cuts" begin
        fig = Figure()
        ax = Axis(fig[1, 1])
        hm = plotWavefrontHeatmap!(ax, r)
        z = wavefrontRemove(r)
        @test isapprox(hm[3][], z; rtol = 1e-5, nans = true)    # Makie stores Float32
        lim = maximum(abs, filter(isfinite, z))
        @test collect(hm.colorrange[]) ≈ [-lim, lim] rtol = 1e-5
        @test lim > OT.MIN_WAVEFRONT_COLOR_RANGE
        hm2 = plotWavefrontHeatmap!(ax, r; remove = :pistonTilt)   # nothing but rounding noise left
        @test collect(hm2.colorrange[]) ≈ [-1e-3, 1e-3]    # floor on the range
        ax3 = Axis3(fig[1, 2])
        sf = plotWavefront3D!(ax3, r; remove = :none)
        @test isapprox(sf[3][], r.opd; rtol = 1e-5, nans = true)
        axc = Axis(fig[2, 1:2])
        xc, yc = plotWavefrontCuts!(axc, r)
        @test last.(xc[1][]) ≈ z[:, mid] nans = true
        @test first.(xc[1][]) ≈ r.u
        @test last.(yc[1][]) ≈ z[mid, :] nans = true
        saveTestFigure(fig, "wavefront plots - heatmap, surface, cuts")
        # an all-masked OPD still plots (default color range)
        none = fill(NaN, size(r.opd))
        @test collect(plotWavefrontHeatmap!(Axis(Figure()[1, 1]), r; opd = none).colorrange[]) ≈ [-1e-3, 1e-3]
    end

    @testset "wavefrontFigure" begin
        fig = wavefrontFigure(r)
        @test fig isa Figure
        @test any(b -> b isa Axis3, fig.content)
        @test count(b -> b isa Axis, fig.content) == 2
        @test any(b -> b isa Colorbar, fig.content)
        label = only(b for b in fig.content if b isa Label)
        @test occursin("RMS", label.text[]) && occursin("Strehl", label.text[])
        saveTestFigure(wavefrontFigure(r; remove = :pistonTilt), "wavefront figure - tilt removed")
        # a defocused system, tilt-free
        gd, fd = parabolaSystem(100.0; shift = 0.2)
        rd = wavefrontMap(gd, fd, λ; N = 33, geoWavelength = λ)
        saveTestFigure(wavefrontFigure(rd), "wavefront figure - defocus")
    end

    @testset "PSF plots" begin
        psf = wavefrontPSF(r; pixels = 65)
        fig = Figure()
        ax = Axis(fig[1, 1])
        hm = plotPSFHeatmap!(ax, psf)
        @test isapprox(hm[3][], psf.intensity; rtol = 1e-5)
        @test plotPSFHeatmap!(ax, psf; scale = :log)[3][] ≈ log10.(max.(psf.intensity, 1e-5)) rtol = 1e-5
        axc = Axis(fig[1, 2])
        xc, yc = plotPSFCuts!(axc, psf; scale = :log)
        @test length(xc[1][]) == 65
        @test_throws ArgumentError plotPSFHeatmap!(ax, psf; scale = :cubic)
        f = psfFigure(psf)
        @test f isa Figure && count(b -> b isa Axis, f.content) == 2
        label = only(b for b in f.content if b isa Label)
        @test occursin("FWHM", label.text[])
        saveTestFigure(f, "PSF figure - tilted wavefront")
        saveTestFigure(psfFigure(psf; scale = :log), "PSF figure - log scale")
    end
end
