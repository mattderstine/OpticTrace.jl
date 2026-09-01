export traceGeometry, traceGeometryRel, surfAmpFunc
export traceGeometry!, traceGeometryRel!

export attributesSurfaces, surfNormal, modFunc

"""
    sag - compute the z coordinate in local coordinates for decendants of
    AbstractSurfProfile

    The type of x,y and the argument to SurfProfieXXX are different so that ForwardDiff
    can be used to find the gradient of the sag function at the intersection point of the
    ray with the surface. The normal vector is then obtained by normalizing the gradient vector. THis is used as a check

    returns z in local coordinates
"""
function sag(x::T, y::T, s::SurfProfileConic{U}) where {T<:Real,U<:Real}
    r2 = (x^2 + y^2)
    sqrtarg = 1 - s.ϵ * s.curv^2 * r2
    if sqrtarg < 0.
        z = NaN
    elseif s.curv == 0.
        z = 0.
    elseif s.ϵ == 0.
        z = s.curv * r2 * 0.5
    else
        #z=(1 - sqrt(1 - s.curv^2 * s.ϵ * (x^2 + y^2)))/(s.curv*s.ϵ)
        z = s.curv * r2 / (1 + sqrt(sqrtarg))
    end
    z
end

"""
    sag(x, y, s::SurfProfileSphere)

Sag of a spherical surface (equivalent to `SurfProfileConic` with
`ϵ = 0`, computed directly rather than by delegating): `z = curv*r² /
(1 + sqrt(1 - curv²*r²))` where `r² = x²+y²`. Returns `NaN` if `(x,y)`
is beyond the sphere's radius (`1 - curv²*r² < 0`). See
`sag(x, y, s::SurfProfileConic)`'s docstring above for the general
x/y/s/return contract shared by every `sag` method.
"""
function sag(x::T, y::T, s::SurfProfileSphere{U}) where {T<:Real,U<:Real}
    r2 = (x^2 + y^2)
    sqrtarg = 1 - s.curv^2 * r2
    if sqrtarg < 0.
        z = NaN
    elseif s.curv == 0.
        z = 0.
    else
        z = s.curv * r2 / (1 + sqrt(sqrtarg))
    end
    z
end

"""
    sag(x, y, s::SurfProfileOAConic)

Sag of an off-axis conic surface: shifts `(x,y)` by `s.offset[1:2]`,
evaluates the equivalent `SurfProfileConic(s.curv, s.ϵ)`'s sag there,
and adds `s.offset[3]` to the result. See
`sag(x, y, s::SurfProfileConic)`'s docstring above for the general
x/y/s/return contract shared by every `sag` method.
"""
function sag(x::T, y::T, s::SurfProfileOAConic{U}) where {T<:Real,U<:Real}
    #println("at sag for SurfProfileOAConic - offset = $(s.offset)")
    a = SurfProfileConic(s.curv, s.ϵ)
    sg = sag(x - s.offset[1], y - s.offset[2], a) + s.offset[3] #add to offset for true sag
    #println("sag = $sg")
    sg
end


"""
    sag - compute the z coordinate in local coordinates for decendants of
    AbstractSurfProfile

    returns z in local coordinates
    For SurfProfileAsphere, compute from 3rd order term onwards


"""
function sag(x::T, y::T, s::OpticTrace.SurfProfileAsphere{U}) where {T<:Real,U<:Real}
    r2 = (x^2 + y^2)
    r = sqrt(r2)
    sqrtarg = 1 - s.ϵ * s.curv^2 * r2
    if sqrtarg < 0.
        return NaN
    end
    #sg = s.curv * r2 /(1+ sqrt(sqrtarg))+sum([ss * r2 * r^i for (i,ss) in enumerate(s.a)])

    asp = 0.
    for ss in Iterators.reverse(s.a)
        asp = (asp + ss) * r
    end
    asp *= r2 #apsehreic coefficients array starts at r3
    sg = s.curv * r2 / (1 + sqrt(sqrtarg)) + asp

    sg
end

"""
    sag(x, y, s::SurfProfileEvenAsphere)

Sag of an even-aspheric surface: a `SurfProfileConic(s.curv, s.ϵ)` base
plus a polynomial correction evaluated from `s.a`, whose coefficients
are even orders starting at 4th order (see `SurfProfileEvenAsphere`).
Returns `NaN` if `(x,y)` is beyond the base conic's domain. Compare
`sag(x, y, s::SurfProfileAsphere)` (same idea, but `s.a`'s coefficients
start at 3rd order and include odd orders). See
`sag(x, y, s::SurfProfileConic)`'s docstring above for the general
x/y/s/return contract shared by every `sag` method.
"""
function sag(x::T, y::T, s::SurfProfileEvenAsphere{U}) where {T<:Real,U<:Real}
    r2 = (x^2 + y^2)

    sqrtarg = 1 - s.ϵ * s.curv^2 * r2
    if sqrtarg < 0.
        return NaN
    end
    #sg = s.curv * r2 /(1+ sqrt(sqrtarg))+sum([ss * r2 * r^i for (i,ss) in enumerate(s.a)])

    asp = 0.
    for ss in Iterators.reverse(s.a)
        asp = (asp + ss) * r2
    end
    asp *= r2 # s starts at order 4 unlike zemax which starts at order 2
    sg = s.curv * r2 / (1 + sqrt(sqrtarg)) + asp

    sg
end

"""
    sag(x, y, s::SurfProfileOddAsphere)

Sag of an "odd asphere" surface: a `SurfProfileConic(s.curv, s.ϵ)` base
plus a polynomial correction evaluated from `s.a`, whose coefficients
start at 1st order (`s.a[i]` is the coefficient of `r^i`) -- unlike
`sag(x, y, s::SurfProfileAsphere)` (starts at 3rd order), this needs no
extra `r²` factor after the Horner evaluation, since there's no order
offset to account for. Returns `NaN` if `(x,y)` is beyond the base
conic's domain. See `sag(x, y, s::SurfProfileConic)`'s docstring above
for the general x/y/s/return contract shared by every `sag` method.
"""
function sag(x::T, y::T, s::SurfProfileOddAsphere{U}) where {T<:Real,U<:Real}
    r2 = (x^2 + y^2)
    r = sqrt(r2)

    sqrtarg = 1 - s.ϵ * s.curv^2 * r2
    if sqrtarg < 0.
        return NaN
    end

    asp = zero(promote_type(T, U))
    for ss in Iterators.reverse(s.a)
        asp = (asp + ss) * r
    end
    s.curv * r2 / (1 + sqrt(sqrtarg)) + asp
end

"""
    xyPolyTermPowers(k::Int) -> (m::Int, n::Int)

Zemax's standard bivariate polynomial term ordering, used by
`SurfProfileXYPoly` (Zemax `TYPE XPOLYNOM`): for 1-based term index `k`,
returns the `(m, n)` exponents of that term's `x^m * y^n` monomial.
Terms are ordered degree-by-degree (`m+n` increasing), with `m`
descending within each degree: `k=1 -> (1,0)` (`x`), `k=2 -> (0,1)`
(`y`), `k=3 -> (2,0)` (`x²`), `k=4 -> (1,1)` (`xy`), `k=5 -> (0,2)`
(`y²`), `k=6 -> (3,0)` (`x³`), and so on. Generic in `k` (not limited to
any fixed term count), so it keeps working regardless of how many terms
a given `XPOLYNOM` surface actually defines.
"""
function xyPolyTermPowers(k::Int)
    d = 1
    while k > d + 1
        k -= d + 1
        d += 1
    end
    m = d - (k - 1)
    n = d - m
    m, n
end

"""
    sag(x, y, s::SurfProfileXYPoly)

Sag of a general 2D polynomial surface: a `SurfProfileConic(s.curv,
s.ϵ)` base plus `Σ s.a[k] * (x/s.normRadius)^m * (y/s.normRadius)^n`,
where `(m,n) = xyPolyTermPowers(k)` for each term `k`. Zero coefficients
are skipped (both for speed and to avoid a `0^0` ambiguity for any term
whose `(m,n)` includes a zero exponent at `x=0`/`y=0`). Returns `NaN` if
`(x,y)` is beyond the base conic's domain. See `sag(x, y,
s::SurfProfileConic)`'s docstring above for the general x/y/s/return
contract shared by every `sag` method.
"""
function sag(x::T, y::T, s::SurfProfileXYPoly{U}) where {T<:Real,U<:Real}
    base = sag(x, y, SurfProfileConic(s.curv, s.ϵ))
    isnan(base) && return base

    ρx = x / s.normRadius
    ρy = y / s.normRadius
    poly = zero(promote_type(T, U))
    for (k, coeff) in enumerate(s.a)
        coeff == 0 && continue
        m, n = xyPolyTermPowers(k)
        poly += coeff * ρx^m * ρy^n
    end
    base + poly
end

"""
    sag(x, y, s::SurfProfileCyl)

Sag of a cylindrical surface: a conic cross-section in `y` only (`x`
does not appear in the formula), using the same
`SurfProfileConic`-style sag equation with `s.curv`/`s.ϵ`. See
`sag(x, y, s::SurfProfileConic)`'s docstring above for the general
x/y/s/return contract shared by every `sag` method.
"""
function sag(x::T, y::T, s::SurfProfileCyl{U}) where {T<:Real,U<:Real}
    if s.ϵ == 0. || s.curv == 0.
        z = s.curv * (y^2) * 0.5
    else
        z = (1 - sqrt(1 - s.curv^2 * s.ϵ * (y^2))) / (s.curv * s.ϵ)
    end
    z
end

"""
    sag(x, y, s::SurfProfileToroid)

Sag of a toroidal surface: the sag of the base y-z conic curve
(`s.curvY`/`s.ϵY`, evaluated at `y` alone, the same closed form as
`sag(x, y, s::SurfProfileConic)` restricted to y) plus the sag of a
circular arc of curvature `s.curvX` evaluated at `x` alone -- these two
contributions are independent and additive, matching Zemax's own
`TOROIDAL` surface definition (the y-z profile curve swept along local
x by a circular arc of radius `1/s.curvX`). At `x=0` this reduces
exactly to the base y-z conic's own sag; when `s.curvX == 0` it reduces
to a pure extrusion of the y-z curve along x (no x sweep at all, per
`SurfProfileToroid`'s own docstring). See `sag(x, y,
s::SurfProfileConic)`'s docstring above for the general x/y/s/return
contract shared by every `sag` method.
"""
function sag(x::T, y::T, s::SurfProfileToroid{U}) where {T<:Real,U<:Real}
    if s.curvY == 0.
        zyz = zero(promote_type(T, U))
    elseif s.ϵY == 0.
        zyz = s.curvY * y^2 * 0.5
    else
        sqrtargY = 1 - s.ϵY * s.curvY^2 * y^2
        zyz = sqrtargY < 0. ? oftype(sqrtargY, NaN) : s.curvY * y^2 / (1 + sqrt(sqrtargY))
    end

    if s.curvX == 0.
        zx = zero(promote_type(T, U))
    else
        sqrtargX = 1 - s.curvX^2 * x^2
        zx = sqrtargX < 0. ? oftype(sqrtargX, NaN) : s.curvX * x^2 / (1 + sqrt(sqrtargX))
    end

    zyz + zx
end


"""
    deltaToSurf - find distance along Ray to OptSurface
                calculation is in local coordinates
    returns distance to intersection of ray with surface

"""
function deltaToSurf(r::Ray{3,T}, p::SurfProfileConic{T}) where T<:Real
    x0, y0, z0 = r.base
    L, M, N = r.dir

    ffunc = p.curv * (x0^2 + y0^2 + z0^2 * p.ϵ) - 2 * z0
    gfunc = N - p.curv * (L * x0 + M * y0 + N * z0 * p.ϵ)

    if p.ϵ == 1.
        C = p.curv
    else
        C = p.curv * (L^2 + M^2 + p.ϵ * N^2)
    end
    #=
    if debugFlag
        println("N = $N  gfunc = $gfunc  ffunc = $ffunc")
        println("p.curv = $(p.curv)  C = $C")
    end
    =#

    if isapprox(C, 0., atol=1e-16)
        if isapprox(N, 0., atol=1e-16)
            qmiss = gfunc^2 - p.curv * ffunc
            if qmiss < 0. || p.curv == 0.
                Δ = NaN # root is imaginary or ray parallel to plane, miss
            else
                Δ = (gfunc - sqrt(qmiss)) / p.curv
            end
        else
            Δ = (0.5 * p.curv * (x0^2 + y0^2) - z0) / N # negative sign removed 9/29/20
        end
    else
        qmiss = gfunc^2 - C * ffunc
        if qmiss < 0.
            Δ = NaN # root is imaginary, miss
        else
            Δ = (gfunc - sqrt(qmiss)) / C
        end
    end
    Δ
end

"""
    deltaToSurf(r::Ray{3,T}, p::SurfProfileSphere{T}) where T<:Real

Distance along `r` to its intersection with a `SurfProfileSphere`, in
local coordinates -- the same quadratic-in-`Δ` solve as
`deltaToSurf(r, p::SurfProfileConic)` (this file's canonical
`deltaToSurf` docstring, above), specialized to `ϵ = 0`.
"""
function deltaToSurf(r::Ray{3,T}, p::SurfProfileSphere{T}) where T<:Real
    x0, y0, z0 = r.base
    L, M, N = r.dir

    ffunc = p.curv * (x0^2 + y0^2 + z0^2 ) - 2 * z0
    gfunc = N - p.curv * (L * x0 + M * y0 + N * z0)


    C = p.curv


    if isapprox(C, 0., atol=1e-16)
        if isapprox(N, 0., atol=1e-16)
            qmiss = gfunc^2 - p.curv * ffunc
            if qmiss < 0. || p.curv == 0.
                Δ = NaN # root is imaginary or ray parallel to plane, miss
            else
                Δ = (gfunc - sqrt(qmiss)) / p.curv
            end
        else
            Δ = (0.5 * p.curv * (x0^2 + y0^2) - z0) / N # negative sign removed 9/29/20
        end
    else
        qmiss = gfunc^2 - C * ffunc
        if qmiss < 0.
            Δ = NaN # root is imaginary, miss
        else
            Δ = (gfunc - sqrt(qmiss)) / C
        end
    end
    Δ
end

# logic is flawed in this one
#change to add offset to ray to put it into the coordinate system of the offset parabola
"""
    deltaToSurf(r::Ray{3,T}, p::SurfProfileOAConic{T}) where T<:Real

Distance along `r` to its intersection with an off-axis conic surface:
offsets `r.base` by `p.offset` and delegates to
`deltaToSurf(r, ::SurfProfileConic)` with the equivalent
`SurfProfileConic(p.curv, p.ϵ)`.

**The preceding code comments (`# logic is flawed in this one` /
`#change to add offset to ray to put it into the coordinate system of
the offset parabola`) are the original author's own note that this is
believed incorrect** -- treat this method's results with suspicion
until that's investigated (see `TODO.md`).
"""
function deltaToSurf(r::Ray{3,T}, p::SurfProfileOAConic{T}) where T<:Real
    #=
    if debugFlag
        println("deltaToSurf OAConic")
        println("base = $(r.base)  dir = $(r.dir)  offset = $(p.offset)  net = $(r.base .- p.offset)")
    end
    =#
    de = deltaToSurf(Ray(r.base .+ p.offset, r.dir), SurfProfileConic(p.curv, p.ϵ))
    #=
    if debugFlag
        println("delta = $de")
    end
    =#
    de
end

"""
    deltaToSurf(r::Ray{3,T}, profile::AbstractAsphericProfile{T}) where T<:Real

Distance along `r` to its intersection with an aspheric surface
(`SurfProfileAsphere` or `SurfProfileEvenAsphere`): finds a numeric
root of `sag(x0+Lδ, y0+Mδ, profile) - z0 - Nδ = 0` (via
`Roots.find_zero`), using the equivalent base conic's `deltaToSurf`
solution as the initial guess.
"""
function deltaToSurf(r::Ray{3,T}, profile::AbstractAsphericProfile{T}) where T<:Real
    x0, y0, z0 = r.base
    L, M, N = r.dir

    guess = deltaToSurf(r, SurfProfileConic(profile.curv, profile.ϵ))
    #=
    if debugFlag
        println("guess = $guess")
    end
    =#
    f(δ) = sag(x0 + L * δ, y0 + M * δ, profile) - z0 - N * δ
    Δl = find_zero(f, guess)
    Δl
end

"""
    deltaToSurf(r::Ray{3,T}, profile::SurfProfileCyl{T}) where T<:Real

Computes the distance along `r` to its intersection with a
`SurfProfileCyl`, following the same quadratic-in-`Δ` solve as
`deltaToSurf(r, p::SurfProfileConic)` restricted to `y`/`z`.
"""
function deltaToSurf(r::Ray{3,T}, profile::SurfProfileCyl{T}) where T<:Real
    x0, y0, z0 = r.base
    L, M, N = r.dir

    ffunc = profile.curv * (y0^2 + z0^2 * profile.ϵ) - 2 * z0
    gfunc = N - profile.curv * (M * y0 + N * z0 * profile.ϵ)

    if profile.ϵ == 1.
        C = profile.curv
    else
        C = profile.curv * (M^2 + profile.ϵ * N^2)
    end
    #=
    if debugFlag
        println("N = $N  gfunc = $gfunc  ffunc = $ffunc")
        println("profile.curv = $(profile.curv)  C = $C")
    end
    =#


    if isapprox(C, 0., atol=1e-16)
        if isapprox(N, 0., atol=1e-16)
            qmiss = gfunc^2 - profile.curv * ffunc
            if qmiss < 0. || profile.curv == 0.
                Δ = NaN # root is imaginary or ray parallel to plane, miss
            else
                Δ = (gfunc - sqrt(qmiss)) / profile.curv
            end
        else
            Δ = (0.5 * profile.curv * (x0^2 + y0^2) - z0) / N # negative sign removed 9/29/20
        end
    else
        qmiss = gfunc^2 - C * ffunc
        if qmiss < 0.
            Δ = NaN # root is imaginary, miss
        else
            Δ = (gfunc - sqrt(qmiss)) / C
        end
    end
    Δ
end

"""
    deltaToSurf(r::Ray{3,T}, p::SurfProfileToroid{T}) where T<:Real

Distance along `r` to its intersection with a `SurfProfileToroid`: no
closed-form root exists for a general toroid (unlike the conic/sphere/
cylinder cases above), so this follows the same numerical pattern used
for `AbstractAsphericProfile` (`deltaToSurf(r, ::AbstractAsphericProfile)`
above) -- a `Roots.find_zero` refinement of `sag(x,y,p) - z0 - Nδ = 0`,
seeded from the base y-z conic's own (closed-form) intersection as the
initial guess.
"""
function deltaToSurf(r::Ray{3,T}, p::SurfProfileToroid{T}) where T<:Real
    x0, y0, z0 = r.base
    L, M, N = r.dir

    guess = deltaToSurf(r, SurfProfileConic(p.curvY, p.ϵY))
    f(δ) = sag(x0 + L * δ, y0 + M * δ, p) - z0 - N * δ
    find_zero(f, guess)
end





"""
    surfNormal - computes normal vectors for use by Makie to texture when rendering
    returns a unit vector
    returns normal in LOCAL coordinates
"""
function surfNormal(r::Point3{T}, s::SurfProfileConic{T}) where T<:Real
    #println("sag = $(sag(r[1],r[2],s))")

    if s.ϵ == 1
        denomI = 1.0
    else
        denomsq = 1 + s.curv * r[3] * (s.ϵ - 1) * (s.curv * s.ϵ * r[3] - 2)
        #next one is directly from Welford 4.30
        #denomsq = 1-2s.curv *(s.ϵ - 1) * r[3] + s.curv^2 * s.ϵ *(s.ϵ -1) * r[3]^2
        if denomsq < 0.0
            #println("r = $r   c = $(s.curv)  ϵ = $(s.ϵ)")
            #println("denomsq = $denomsq")
            denomI = NaN
        else
            denomI = 1 / sqrt(denomsq)
        end
    end
    Vec3(-s.curv * r[1] * denomI, -s.curv * r[2] * denomI,
        (1.0 - s.curv * s.ϵ * r[3]) * denomI)
end

"""
    surfNormal(r::Point3{T}, s::SurfProfileSphere{T}) where T<:Real

Surface normal of a `SurfProfileSphere` at local point `r` (which must
lie on the sphere, e.g. `r[3] == sag(r[1], r[2], s)`), in local
coordinates. See `surfNormal(r, s::SurfProfileConic)`'s docstring above
for the general contract shared by every `surfNormal` method.
"""
function surfNormal(r::Point3{T}, s::SurfProfileSphere{T}) where T<:Real
    #println("sag = $(sag(r[1],r[2],s))")
    grad = Vec3(-s.curv * r[1], -s.curv * r[2], (1.0 - s.curv * r[3]))
    #normalize(Vec3(-s.curv * r[1], -s.curv * r[2], (1.0 - s.curv * r[3])))
end

"""
    surfNormal(r::Point3{T}, s::SurfProfileOAConic{T}) where T<:Real

Surface normal of an off-axis conic surface: offsets `r` by `s.offset`
and delegates to `surfNormal(r, ::SurfProfileConic)` with the
equivalent `SurfProfileConic(s.curv, s.ϵ)`.
"""
function surfNormal(r::Point3{T}, s::SurfProfileOAConic{T}) where T<:Real
    #println("OA r = $r offset = $(s.offset) net = $(r .- s.offset)")
    surfNormal(r .- s.offset, SurfProfileConic(s.curv, s.ϵ))
end

"""
    surfNormal(rr::Point3{T}, s::SurfProfileAsphere{T}) where T<:Real

Surface normal of a `SurfProfileAsphere` at local point `rr`, in local
coordinates -- the analytic gradient of
`sag(x, y, s::SurfProfileAsphere)`, normalized. See
`surfNormal(r, s::SurfProfileConic)`'s docstring above for the general
contract shared by every `surfNormal` method.
"""
function surfNormal(rr::Point3{T}, s::SurfProfileAsphere{T}) where T<:Real
    x = rr[1]
    y = rr[2]
    z = rr[3]
    r2 = (x^2 + y^2)
    r = sqrt(r2)

    asp = 0.
    for (i, ss) in Iterators.reverse(enumerate(s.a))
        asp = (asp + (i + 2) * ss) * r
        #println("i = $i  ss = $ss  asp = $asp")
    end
    #asp *= r  # dr = x/r dx + y/r dy so skip this mulitply

    norm = normalize([-x * (s.curv + asp), -y * (s.curv + asp), 1 - s.curv * s.ϵ * z])
    Vec3(norm[1], norm[2], norm[3])
end

"""
    surfNormal(rr::Point3{T}, s::SurfProfileEvenAsphere{T}) where T<:Real

Surface normal of a `SurfProfileEvenAsphere` at local point `rr`, in
local coordinates -- the analytic gradient of
`sag(x, y, s::SurfProfileEvenAsphere)`, normalized. See
`surfNormal(r, s::SurfProfileConic)`'s docstring above for the general
contract shared by every `surfNormal` method.
"""
function surfNormal(rr::Point3{T}, s::SurfProfileEvenAsphere{T}) where T<:Real
    x = rr[1]
    y = rr[2]
    z = rr[3]
    r2 = (x^2 + y^2)

    #aspheresum = sum([i * r^(i-2) * s.a[i-2]  for i in 3:2+length(s.a)])
    #aspheresum = sum([(i+2) * r^i * s.a[i]  for i in 1:length(s.a)])
    asp = 0.

    for (i, ss) in Iterators.reverse(enumerate(s.a))
        asp = (asp + (2 * (i + 1)) * ss) * r2

    end
   #asp *= sqrt(r2)

    norm = normalize([-x * (s.curv + asp), -y * (s.curv + asp), 1 - s.curv * s.ϵ * z])
    Vec3(norm[1], norm[2], norm[3])
end

"""
    surfNormal(r::Point3{T}, s::AbstractAsphericProfile{T}) where T<:Real

Generic fallback `surfNormal` for any `AbstractAsphericProfile` subtype
that doesn't have its own dedicated method (currently
`SurfProfileOddAsphere`/`SurfProfileXYPoly`; `SurfProfileAsphere`/
`SurfProfileEvenAsphere` have more specific closed-form methods above
that Julia dispatches to instead): the `ForwardDiff` gradient of `sag`
at `r`, normalized -- mirrors the `gbRadius`/`gbWidths` generic-fallback
pattern in `src/mesh_primitives.jl`, avoiding a hand-derived closed-form
gradient for every new aspheric-family type.
"""
function surfNormal(r::Point3{T}, s::AbstractAsphericProfile{T}) where T<:Real
    grad = ForwardDiff.gradient(xy -> sag(xy[1], xy[2], s), SVector(r[1], r[2]))
    normalize(Vec3(-grad[1], -grad[2], one(T)))
end

"""
    surfNormal(r::Point3{T}, s::SurfProfileCyl{T}) where T<:Real

Surface normal of a `SurfProfileCyl` at local point `r`, in local
coordinates -- the same closed-form conic-normal formula as
`surfNormal(r, s::SurfProfileConic)`, with the `x` component fixed at
`0` (no curvature along local x). See that method's docstring above for
the general contract shared by every `surfNormal` method.
"""
function surfNormal(r::Point3{T}, s::SurfProfileCyl{T}) where T<:Real
    #println("sag = $(sag(r[1],r[2],s))")

    if s.ϵ == 1
        denomI = 1
    else
        denomsq = 1 + s.curv * r[3] * (s.ϵ - 1) * (s.curv * s.ϵ * r[3] - 2)
        #next one is directly from Welford 4.30
        #denomsq = 1-2s.curv *(s.ϵ - 1) * r[3] + s.curv^2 * s.ϵ *(s.ϵ -1) * r[3]^2
        if denomsq < 0.0
            #println("r = $r   c = $(s.curv)  ϵ = $(s.ϵ)")
            #println("denomsq = $denomsq")
            denomI = NaN
        else
            denomI = 1 / sqrt(denomsq)
        end
    end
    Vec3(0., -s.curv * r[2] * denomI, (1.0 - s.curv * s.ϵ * r[3]) * denomI)
end

"""
    surfNormal(r::Point3{T}, s::SurfProfileToroid{T}) where T<:Real

Surface normal of a `SurfProfileToroid` at local point `r`, in local
coordinates -- the analytic gradient of `sag(x, y, s::SurfProfileToroid)`
(itself a sum of an independent x term and y term, so its gradient is
just the two 1D conic derivatives stacked together), normalized. See
`surfNormal(r, s::SurfProfileConic)`'s docstring above for the general
contract shared by every `surfNormal` method.
"""
function surfNormal(r::Point3{T}, s::SurfProfileToroid{T}) where T<:Real
    x, y = r[1], r[2]

    dzdx = s.curvX == 0 ? zero(T) : s.curvX * x / sqrt(1 - s.curvX^2 * x^2)
    dzdy = s.curvY == 0 ? zero(T) : s.curvY * y / sqrt(1 - s.ϵY * s.curvY^2 * y^2)

    normalize(Vec3(-dzdx, -dzdy, one(T)))
end

"""
    traceGeometry - trace a Ray through an array of OptSurfaces
    r       intial ray vector is global coordinates
    geo     the array of surfaces

    returns status, array of results from traceGeometry
"""
function traceGeometry(r::Ray{3,T}, geo) where T<:Real
    trc = Vector{Trace}(undef, length(geo) + 1)
    trc[1] = Trace(r, 1., 0., identityAmpMats()) #save the start ray etc
    curRay = r
    status = 0
    i = 1
    for surf in geo
        i += 1 #first element in trace is [2] in array
        status, trc[i] = traceSurf(curRay, surf)
        if status != 0
            break
        end
        curRay = trc[i].ray
    end
    status, trc[1:i] #only send the good ones!
end

"""
    traceGeometry! - trace a Ray through an array of OptSurfaces
    trc is a Vector{Trace} of length(geo) + 1 (the output)
    r is the intial ray vector is global coordinates
    geo is the array of surfaces

    returns status, length of trace
"""
function traceGeometry!(trc::Vector{Trace}, r::Ray{3,T}, geo) where T<:Real
    #trc = Vector{Trace}(undef, length(geo)+1)
    Trace!(trc[1], r, 1., 0., identityAmpMats()) #save the start ray etc
    curRay = r
    status = 0
    i = 1
    for surf in geo
        i += 1 #first element in trace is [2] in array
        status, - = traceSurf!(trc[i], curRay, surf)
        if status != 0
            break
        end
        curRay = trc[i].ray
    end
    status, i #only send the good ones!
end

"""
    traceGeometryRel - trace a Ray through an array of OptSurfaces
    rr      initial ray vector, in LOCAL coordinates relative to geo[1]
            (converted to global coordinates via geo[1].toGlobalCoord/
            toGlobalDir before tracing)
    geo     the array of surfaces

    returns status, array of results from traceGeometry
"""
function traceGeometryRel(rr::Ray{3,T}, geo) where T<:Real
    trc = Vector{Trace}(undef, length(geo) + 1)
    surf1 = geo[1]
    curRay = Ray(surf1.toGlobalCoord(rr.base), surf1.toGlobalDir(rr.dir))
    trc[1] = Trace(curRay, 1., 0., identityAmpMats()) #save the start ray etc
    status = 0
    i = 1
    for surf in geo
        i += 1 #first element in trace is [2] in array
        status, trc[i] = traceSurf(curRay, surf)
        if status != 0
            break
        end
        curRay = trc[i].ray
    end
    status, trc[1:i] #only send the good ones!
end

"""
    traceGeometryRel! - trace a Ray through an array of OptSurfaces
    trc is a Vector{Trace} of length(geo) + 1 (the output)
    rr is the initial ray vector, in LOCAL coordinates relative to
       geo[1] (converted to global coordinates via geo[1].toGlobalCoord/
       toGlobalDir, then delegates to traceGeometry!)
    geo is the array of surfaces

    returns status, length of trace
"""
function traceGeometryRel!(trc, rr::Ray{3,T}, geo) where T<:Real
    #trc = Vector{Trace}(undef, length(geo)+1)
    surf1 = geo[1]
    curRay = Ray(surf1.toGlobalCoord(rr.base), surf1.toGlobalDir(rr.dir))
    traceGeometry!(trc, curRay, geo)
end


"""
    traceSurf(r::Ray{3,T}, s::OptSurface{3,T}) where T<:Real

    Compute the exit ray from a refracting/reflecting OptSurface

    returns (status, trc::Trace)

    trc.ray is the new ray (base & direction); trc.pmatrix holds the P & O
    polarization matrices
    status
        0 - Normal
        1 - missed
        2 - TIR when refraction is expected - no polarization

    This method does not check `s.aperture`, so it never returns status 3
    (aperture-clipped) -- that status is only returned by the ModelSurface
    methods of traceSurf/traceSurf! (see their docstrings).

    trc will contain best representation of the ray on error so it could be
    used in plotting and to continue a nonsequential raytrace
"""
function traceSurf(r::Ray{3,T}, s::OptSurface{3,T}) where T<:Real
    localRayStart = s.toLocalCoord(r.base)
    localRayDir = s.toLocalDir(r.dir)

    delta = deltaToSurf(Ray(localRayStart, localRayDir), s.profile)
    #=
    if debugFlag
        println("traceSurf...  Δ = $delta")
        println("start = $(r.base)")
    end
    =#
    if isnan(delta) #missed
        return (1, Trace(r, NaN, delta, identityAmpMats()))
    end
    newRayBase = r.base + r.dir * delta
    newLocalBase = localRayStart + localRayDir * delta

    lnormal = surfNormal(newLocalBase, s.profile)

    normal = s.toGlobalDir(lnormal)
    #=
    if debugFlag
        println("local normal = $lnormal  global normal = $normal")
    end
    =#
    t, newRayDir, nIn = modFunc(Ray(newRayBase, r.dir), normal, s.mod)
    if !t #TIR or other failure
        return (2, Trace(Ray(newRayBase, newRayDir), nIn, delta, identityAmpMats()))
    end
    #=
    if debugFlag
        println("old dir = $(r.dir)  newdir = $newRayDir")
    end
    =#
    # can modify amplitude/polarization and direction of the ray, can't change intersection

    ampMats, newRayDir = surfAmpFunc(r.dir, newRayDir, normal, newLocalBase, s.mod, s.coating)

    return (0, Trace(Ray(newRayBase, newRayDir), nIn, delta, ampMats))
end

"""
    Trace!(trc::Trace{T}, ray, index, delta, ampdata) where T<:Real

Overwrite `trc`'s fields (`ray`, `nIn`, `delta`, `pmatrix`) in place
with the given values, and return it. The in-place counterpart to
constructing a new `Trace(ray, index, delta, ampdata)`; used by
`traceGeometry!`/`traceSurf!` to avoid allocating a new `Trace` per
step.
"""
function Trace!(trc::Trace{T}, ray, index, delta, ampdata) where T<:Real
    trc.ray = ray
    trc.nIn = index
    trc.delta = delta
    trc.pmatrix = ampdata
    return trc
end

"""
    traceSurf!(trc, r::Ray{3,T}, s::OptSurface{3,T}) where T<:Real

In-place version of `traceSurf(r, s::OptSurface{3,T})` (see that
method's docstring above for the full status-code contract) -- writes
the result into `trc` via `Trace!` instead of allocating a new `Trace`.

returns (status, trc::Trace)
"""
function traceSurf!(trc, r::Ray{3,T}, s::OptSurface{3,T}) where T<:Real
    localRayStart = s.toLocalCoord(r.base)
    localRayDir = s.toLocalDir(r.dir)

    delta = deltaToSurf(Ray(localRayStart, localRayDir), s.profile)
    #=
    if debugFlag
        println("traceSurf...  Δ = $delta")
        println("start = $(r.base)")
    end
    =#
    if isnan(delta)
        return (1, Trace!(trc, r, NaN, delta, identityAmpMats()))
    end
    newRayBase = r.base + r.dir * delta
    newLocalBase = localRayStart + localRayDir * delta

    lnormal = surfNormal(newLocalBase, s.profile)

    normal = s.toGlobalDir(lnormal)
    #=
    if debugFlag
        println("local normal = $lnormal  global normal = $normal")
    end
    =#
    t, newRayDir, nIn = modFunc(Ray(newRayBase, r.dir), normal, s.mod)
    if !t
        return (2, Trace!(trc, Ray(newRayBase, newRayDir), nIn, delta, identityAmpMats()))
    end
    #=
    if debugFlag
        println("old dir = $(r.dir)  newdir = $newRayDir")
    end
    =#
    # can modify amplitude/polarization and direction of the ray, can't change intersection

    ampMats, newRayDir = surfAmpFunc(r.dir, newRayDir, normal, newLocalBase, s.mod, s.coating)

    return (0, Trace!(trc, Ray(newRayBase, newRayDir), nIn, delta, ampMats))
end

"""
    traceSurf(r::Ray{3,T}, s::ModelSurface{3,T}) where T<:Real

    Compute the exit ray from a non-refracting ModelSurface (a
    pass-through reference/model surface used for characterization --
    the ray's direction is unchanged, only its base point advances to the
    surface and `s.aperture` is checked).

    returns (status, trc::Trace)

    trc.ray is the new ray (base advanced to the surface, direction
    unchanged)
    status
        0 - Normal
        1 - missed
        3 - Blocked by aperture/size of element (`isAperture(s.aperture) &&
            clipAperture(...)`)

    Unlike traceSurf(r, s::OptSurface{3,T}), this method has no TIR/status-2
    case, since no refraction or reflection is computed here.

    trc will contain best representation of the ray on error so it could be
    used in plotting and to continue a nonsequential raytrace
"""
function traceSurf(r::Ray{3,T}, s::ModelSurface{3,T}) where T<:Real
    localRayStart = s.toLocalCoord(r.base)
    localRayDir = s.toLocalDir(r.dir)
    #println("\ntraceSurf - ModelSurface  name = $(s.surfname)")
    #println("aperture = $(s.aperture)")

    delta = deltaToSurf(Ray(localRayStart, localRayDir), s.profile)
    if isnan(delta)
        return (1, Trace(r, s.refIndex, delta, identityAmpMats()))
    end
    newRayBase = r.base + r.dir * delta
    newLocalBase = localRayStart + localRayDir * delta


    if isAperture(s.aperture) && clipAperture(newLocalBase, s.aperture)
        stat = 3 #clipped
    else
        stat = 0 #ok to continue
    end

    return (stat, Trace(Ray(newRayBase, r.dir), s.refIndex, delta, identityAmpMats()))
end

"""
    traceSurf!(trc, r::Ray{3,T}, s::ModelSurface{3,T}) where T<:Real

    In-place version of traceSurf(r, s::ModelSurface{3,T}) (see that
    method's docstring above for the full status-code contract) --
    writes the result into `trc` (via Trace!) instead of allocating a
    new Trace.

    returns (status, trc::Trace)
    status
        0 - Normal
        1 - missed
        3 - Blocked by aperture/size of element (`isAperture(s.aperture) &&
            clipAperture(...)`)
"""
function traceSurf!(trc, r::Ray{3,T}, s::ModelSurface{3,T}) where T<:Real
    localRayStart = s.toLocalCoord(r.base)
    localRayDir = s.toLocalDir(r.dir)
    #println("\ntraceSurf - ModelSurface  name = $(s.surfname)")
    #println("aperture = $(s.aperture)")

    delta = deltaToSurf(Ray(localRayStart, localRayDir), s.profile)
    if isnan(delta)
        return (1, Trace!(trc, r, s.refIndex, delta, identityAmpMats()))
    end
    newRayBase = r.base + r.dir * delta
    newLocalBase = localRayStart + localRayDir * delta


    if isAperture(s.aperture) && clipAperture(newLocalBase, s.aperture)
        stat = 3 #clipped
    else
        stat = 0 #ok to continue
    end

    return (stat, Trace!(trc, Ray(newRayBase, r.dir), s.refIndex, delta, identityAmpMats()))
end




"""
    AMPPERFECTSURFACE

The `AmpData` value representing a surface that doesn't alter
polarization or transmission at all: `p`/`o` are both `identityPol`
and `trans` is `[1.]`. Returned by `identityAmpMats`.
"""
const AMPPERFECTSURFACE = AmpData(identityPol, identityPol, [1.])

"""
    identityAmpMats()

Return `AMPPERFECTSURFACE`, the shared identity `AmpData` value. Used
as the polarization/amplitude result wherever a trace step doesn't (or
can't yet) compute real polarization data -- e.g. missed/TIR traces,
and `ModelSurface` tracing, which doesn't track polarization at all.
"""
identityAmpMats() = AMPPERFECTSURFACE
#=
function identityAmpMats()
    #identity = [[1.0, 0.0, 0.0] [ 0.0, 1.0, 0.0] [ 0.0, 0.0, 1.0]]
    AmpData(identityPol, identityPol, [1.])
    #retirm matrices for P & O and unity transmission
end
=#

"""
    modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::S) where {S<:AbstractBendDielectric,T<:Real}

Refract `ray` at a dielectric boundary with surface normal `normal` and
refractive indices `d.refIndexIn`/`d.refIndexOut` (swapped if the ray
is hitting the surface from the "back" side, `cosI < 0`), via Snell's
law.

Returns `(ok, newDir, nIn)`:
- `ok::Bool` -- `false` on total internal reflection (in which case
  `newDir` is the reflected direction instead of a refracted one);
  `true` otherwise
- `newDir::Vec3` -- the new (unit) ray direction
- `nIn` -- the refractive index the ray was traveling in before this
  surface (needed by the caller for OPD bookkeeping)

See `modFunc(ray, normal, d::S) where S<:AbstractBendMirror`,
`modFunc(ray, normal, d::CDiffuser)`, and
`modFunc(ray, normal, d::NoBendIndex)` below for the other `mod` types
this generic function dispatches on; all share this `(ok, newDir, nIn)`
return shape.
"""
function modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::S) where {S<:AbstractBendDielectric,T<:Real}

    r = ray.base
    a = ray.dir

    cosI = a ⋅ normal
    if cosI < 0. #hitting surface "backwards"
        nIn = d.refIndexOut
        nOut = d.refIndexIn
        cosI = -cosI
    else
        nIn = d.refIndexIn
        nOut = d.refIndexOut
    end

    tirtest = nOut^2 - nIn^2 * (1 - cosI^2)
    if tirtest < 0
        return false, -(a - (2. * (a ⋅ normal)) .* normal), nIn #TIR direction
    end

    cosIp = sqrt(tirtest) / nOut

    kparam = (nOut * cosIp - nIn * cosI)
    return true, (nIn .* a + kparam .* normal) ./ nOut, nIn
end



"""
    modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::S) where {S<:AbstractBendMirror,T<:Real}

Reflect `ray` off a mirror surface with normal `normal` (Welford 4.46).
`d.refIndexIn`/`d.refIndexOut` select which side's index to report as
`nIn`, depending on which way the ray hits the surface (`cosI < 0` or
not) -- a mirror doesn't itself change refractive index.

Returns `(true, newDir, nIn)` -- see
`modFunc(ray, normal, d::S) where S<:AbstractBendDielectric`'s
docstring above for the shared `(ok, newDir, nIn)` return shape (`ok`
is always `true` here; reflection can't "fail" the way refraction can
via TIR).
"""
function modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::S) where {S<:AbstractBendMirror,T<:Real}
    r = ray.base
    a = ray.dir
    cosI = a ⋅ normal
    if cosI < 0. #hitting surface "backwards"
        nIn = d.refIndexOut
    else
        nIn = d.refIndexIn
    end

    true, (a - (2. * (cosI)) .* normal), nIn  # Welford 4.46
end

"""
    modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::CDiffuser{T}) where T<:Real

Scatter `ray` within a cone around its own direction: builds two unit
vectors perpendicular to `ray.dir`, samples a random point in the unit
disk, scales it by `d.tanθ`, and adds the result to `ray.dir` before
renormalizing. `d.refIndexIn`/`d.refIndexOut` select `nIn` the same way
as the other `modFunc` methods.

Returns `(true, newDir, nIn)` -- see
`modFunc(ray, normal, d::S) where S<:AbstractBendDielectric`'s
docstring above for the shared `(ok, newDir, nIn)` return shape.
"""
function modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::CDiffuser{T}) where T<:Real
    r = ray.base
    a = ray.dir
    perpapprox = zeros(Float64, 3)
    absa = abs.(a)
    perpapprox[findfirst(isequal(minimum(absa)), absa)] = 1.  #create the vector almost perpendicular to the direction

    perp1 = normalize!(cross(a, perpapprox))
    perp2 = cross(a, perp1)
    dir = 2 .* rand(2) .- [1., 1.]
    while norm(dir) > 1.
        dir = 2 .* rand(2) .- [1., 1.]
    end
    displace = d.tanθ .* (dir[1] .* perp1 + dir[2] .* perp2)

    cosI = a ⋅ normal
    if cosI < 0. #hitting surface "backwards"
        nIn = d.refIndexOut
    else
        nIn = d.refIndexIn
    end

    true, normalize(a + displace), nIn
end


"""
    surfNormal(r::Point3{T}, s::NoProfile{T}) where T<:Real

Surface normal of a flat `NoProfile` plane: always `(0,0,1)`,
regardless of `r`. See `surfNormal(r, s::SurfProfileConic)`'s docstring
above for the general contract shared by every `surfNormal` method.
"""
function surfNormal(r::Point3{T}, s::NoProfile{T}) where T<:Real
    Vec3(0., 0., 1.)
end

"""
    deltaToSurf(r::Ray{3,T}, p::NoProfile{T}) where T<:Real

Distance along `r` to its intersection with the local `z=0` plane
(`NoProfile`'s implicit flat surface): `Δ = -z0/N`. Returns `NaN` if
the ray is parallel to the plane (`N ≈ 0`). See
`deltaToSurf(r, p::SurfProfileConic)`'s docstring above for the general
contract shared by every `deltaToSurf` method.
"""
function deltaToSurf(r::Ray{3,T}, p::NoProfile{T}) where T<:Real
    x0, y0, z0 = r.base
    L, M, N = r.dir

    if N ≈ 0.
        Δ = NaN #ray parallel to flat surface
    else
        Δ = -z0 / N
    end
    Δ
end

"""
    deltaToSurf(r::Ray{3,T}, p::ParaxialProfile{T}) where T<:Real

Distance along `r` to its intersection with the local `z=0` plane
(`ParaxialProfile`'s implicit flat surface, same as `NoProfile`'s) --
identical formula to `deltaToSurf(r, ::NoProfile)`. See
`deltaToSurf(r, p::SurfProfileConic)`'s docstring above for the general
contract shared by every `deltaToSurf` method.
"""
function deltaToSurf(r::Ray{3,T}, p::ParaxialProfile{T}) where T<:Real
    x0, y0, z0 = r.base
    L, M, N = r.dir

    if N ≈ 0.
        Δ = NaN #ray parallel to flat surface
    else
        Δ = -z0 / N
    end
    Δ
end

"""
    surfNormal(r::Point3{T}, s::ParaxialProfile{T}) where T<:Real

**Not a true normal.** Returns the local intersection coordinates
`(x, y, 0)` unchanged (not unit length) -- see `ParaxialProfile`'s own
docstring (`src/lens_definitions.jl`) for why: an ideal thin lens's
`modFunc(::ParaxialLensT)` needs the ray's position on the lens, not a
surface normal, and this is how that position reaches it through
`traceSurf`'s existing (otherwise unmodified) `surfNormal` ->
`s.toGlobalDir` -> `modFunc` pipeline. Because `GeometryBasics.normals`
(`src/mesh_primitives.jl`) also calls this generic `surfNormal`, that
file has a dedicated, more specific method for `ParaxialProfile`-
profiled surfaces that returns the real `(0,0,1)` normal for mesh
shading instead of calling this method -- don't remove that override
without keeping shading normals correct in mind.
"""
function surfNormal(r::Point3{T}, s::ParaxialProfile{T}) where T<:Real
    Vec3(r[1], r[2], zero(T))
end

"""
    modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::NoBendIndex{T}) where T<:Real

Pass `ray` through unchanged: returns `(true, ray.dir, d.refIndexIn)`.
See `modFunc(ray, normal, d::S) where S<:AbstractBendDielectric`'s
docstring above for the shared `(ok, newDir, nIn)` return shape.
"""
function modFunc(ray::Ray{3,T}, normal::Vec3{T}, d::NoBendIndex{T}) where T<:Real
    true, ray.dir, d.refIndexIn
end

"""
    modFunc(ray::Ray{3,T}, offset::Vec3{T}, bend::ParaxialLensT{T}) where T<:Real

Ideal thin-lens ray transfer for a `ParaxialLensT` bend. `offset` is
**not a surface normal** -- per `ParaxialProfile`'s docstring
(`src/lens_definitions.jl`), it's the ray's global transverse offset
from the lens's own optical axis, computed by `traceSurf` from
`surfNormal(::ParaxialProfile)`'s repurposed return value the same way
every other bend type's real normal is computed, just carrying
different information. The bend itself is a direct vector form of the
textbook paraxial thin-lens transfer law (`slope' = slope - height/f`):

    newDir = normalize(ray.dir - offset / bend.focalLength)

This is exact for near-axial rays and an approximation farther
off-axis (reusing the `normal` argument slot for position means this
method has no way to decompose `ray.dir` into axial/transverse
components the way an exact-for-any-angle ideal-lens formula would
need) -- an intentional scope match to a surface literally named
"paraxial", not a shortfall.

Always succeeds (`true`, no TIR-style failure mode for an ideal lens),
and always reports `bend.refIndexOut` as `nIn`: unlike `DielectricT`/
`MirrorR`, this method has no true normal available to test which
physical side a ray hit from, so it can't pick between
`refIndexIn`/`refIndexOut` the way those methods do. A narrow,
documented simplification -- real Zemax `PARAXIAL` samples sit in an
unchanged medium on both sides, so `refIndexIn == refIndexOut` in
practice anyway. See `modFunc(ray, normal, d::S) where
S<:AbstractBendDielectric`'s docstring above for the shared `(ok,
newDir, nIn)` return shape.
"""
function modFunc(ray::Ray{3,T}, offset::Vec3{T}, bend::ParaxialLensT{T}) where T<:Real
    true, normalize(ray.dir - offset / bend.focalLength), bend.refIndexOut
end



"""
    surfAmpFunc(dirIn::Vec3{T}, dirOut::Vec3{T}, normal::Vec3{T}, newLocalBase::Point3{T}, sndex::B, amp) where {T<:Real,B<:AbstractBendType{T}}

Compute the amplitude/polarization transfer data and (possibly
modified) outgoing direction for a trace step. This catch-all method
does no actual polarization/amplitude calculation -- it returns
`(identityAmpMats(), dirOut)` unchanged, for every `mod`/`amp`
combination. (Two more specific methods, for `MirrorR`/`CDiffuser`, are
commented out below with the same identity behavior -- polarization
tracking isn't implemented for any surface type yet.)

Arguments:
- `dirIn`, `dirOut` -- incoming/outgoing ray direction
- `normal` -- local surface normal at the intersection
- `newLocalBase` -- local intersection point
- `sndex` -- the surface's `mod` (`B<:AbstractBendType{T}`)
- `amp` -- the surface's coating/amplitude parameter
"""
function surfAmpFunc(dirIn::Vec3{T}, dirOut::Vec3{T}, normal::Vec3{T}, newLocalBase::Point3{T}, sndex::B, amp) where {T<:Real,B<:AbstractBendType{T}}
    identityAmpMats(), dirOut
end

#=
function surfAmpFunc(dirIn::Vec3, dirOut::Vec3, normal::Vec3, newLocalBase::Point3, ndex::MirrorR,amp::AmpParam)
    identityAmpMats(), dirOut
end

function surfAmpFunc(dirIn::Vec3, dirOut::Vec3, normal::Vec3, newLocalBase::Point3, ndex::CDiffuser,amp::AmpParam)
    identityAmpMats(), dirOut
end

=#

"""
    attributesSurfaces

Global `Dict{String,Any}` cache mapping coating names to
`AmpParam`/`NoAmpParam` objects, populated lazily by `getAmpParams` the
first time a given coating name is seen. Threaded through as a keyword
argument (`attributesSurfaces=...`) by most surface constructors in
`src/surfaces.jl`.
"""
attributesSurfaces::Dict{String,Any} = Dict{String,Any}()

"""
    getAmpParams(s::String; attributesSurfaces=attributesSurfaces)

Look `s` up in `attributesSurfaces`, creating (and caching) a new
`AmpParam(s)` if it isn't already present. Called when a surface's
`coating` is given as a raw string name; see
`getAmpParams(s::S) where S<:AbstractAmplitudeParam`'s docstring below
for the other half of the `coating::APorString` dispatch.
"""
function getAmpParams(s::String; attributesSurfaces=attributesSurfaces)
    amp = get!(attributesSurfaces, s, AmpParam(s))
    return amp
end

"""
    getAmpParams(s::S; attributesSurfaces=attributesSurfaces) where S<:AbstractAmplitudeParam

Return `s` unchanged -- called when a surface's `coating` is already an
`AbstractAmplitudeParam` object (e.g. `AmpParam`, `NoAmpParam`) rather
than a string name. The `attributesSurfaces` keyword is accepted for a
uniform call signature with `getAmpParams(s::String)`'s docstring
above, but unused here.
"""
function getAmpParams(s::S; attributesSurfaces=attributesSurfaces) where S<:AbstractAmplitudeParam
    return s
end
