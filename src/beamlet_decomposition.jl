#beamlet decomposition

struct Beamlet{T,N}
    base::SurfBase{T,N}
    amp::T
    wx::T
    wy::T
    ϕ::T
    curvx::T
    curvy::T
    chief::Ray{T,N}
    rayWx::Ray{T,N}
    rayWy::Ray{T,N}
    rayDx::Ray{T,N}
    rayDy::Ray{T,N}
end

#=
1. Define input amplitude and phase funciton (initial Gaussian beam)
1. Define an array of beamlets
3. Assign curvature and direction to beamlets
4. Define sampling plane
5. Compute intrinsic complex amplitude of beamlets at sampling points (Matrix A)
6. Compute complex amplitude of input function at sampling points (Vector B)
7. Compute complex coefficients for each beamlet (compute psuedo inverse of A * B)
8. Verify fit by reconstructing input from beamlets at sampling points
9. Propagate beamlets to waist of input and reconstruct. Compare

Functions needed

Input phase, direction, curvature
beamlet sampling (positions)
Assign curvature and direction to beamlets
Compute Matrix A
Compute Vector B
Reconstruct field from beamlets
Propagate beamlets (initial translation then raytrace)

=#

"""
gaussBeamParams(base, dir, w, curv, λ)
    base is the global coordinate of the origin
    dir is the diredtion of propagation
    curv is the curvature (1/R(z)) 
    λ is the wavelength

    Returns 
"""
function gaussBeamParams(base, dir, w, curv, λ)

end

    