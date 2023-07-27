


const origin = SVector(0., 0., 0.)
const zaxis = SVector(0., 0., 1.)
const yaxis = SVector(0.,1.,0.)
const xaxis = SVector(1.,0.,0.)

const ∞ = Inf

struct Ray
    base::SVector{3, Float64}
    dir::SVector{3, Float64}
end

struct SurfBase
    base::SVector{3, Float64}
    dir::SVector{3, Float64}
    ydir::SVector{3, Float64}
end

abstract type AbstractSize end
abstract type AbstractSurfProfile end
abstract type AbstractAmplitudeParam end
abstract type AbstractBendType end #some info on how to do the bending
abstract type AbstractOpticalObject end

struct SizeLens <: AbstractSize  #size for round optics
    semiDiameter::Float32
end

abstract type TrueAperture <: AbstractSize end

struct RoundAperture <: TrueAperture
    obscure::Float32
    semiDiameter::Float32
end

struct RectAperture <: TrueAperture
    wo::Float32
    lo::Float32
    wclear::Float32
    lclear::Float32
end

#=
struct SurfProfileSphere <: AbstractSurfProfile
    curv :: Float64
end
=#

struct SurfProfileConic <: AbstractSurfProfile
    curv :: Float64
    ϵ::Float64 #see Welford for definition of ϵ
end

struct SurfProfileOAConic <: AbstractSurfProfile
    curv :: Float64
    ϵ::Float64 #see Welford for definition of ϵ
    offset::SVector{3, Float64}
end

struct SurfProfileAsphere <: AbstractSurfProfile
    curv :: Float64
    ϵ::Float64 #see Welford for definition of ϵ
    a::AbstractArray{Float64} #hope this works...
end

struct SurfProfileCyl <: AbstractSurfProfile
    curv :: Float64
    ϵ::Float64 #see Welford for definition of ϵ    ϵy::Float64 #see Welford for definition of ϵ
    a::AbstractArray{Float64} #hope this works...
end

struct NoProfile <: AbstractSurfProfile
    curv :: Float64 #but never used
end


struct AmpParam <: AbstractAmplitudeParam
    type::String
end

abstract type AbstractAmpData end

struct AmpData <: AbstractAmpData
    p::Array{Float64,2}
    o::Array{Float64,2}
    trans::Float64 #for intensity systems, the "transmission" coefficient of the surface
end


"""
    Trace(ray, nIn, delta, pmatrix)
    Element of a raytrace.
        ray     location of interesection with surface and exiting direction
        nIn     refractive index leaving surface
        delta   length of ray leaving surface
        pmatrix polarization p & o matrices
"""
struct Trace
    ray::Ray
    nIn::Float64
    delta::Float64
    pmatrix::AbstractAmpData
end


struct DielectricT <: AbstractBendType
    refIndexIn :: Float64
    refIndexOut :: Float64
end

struct MirrorR <: AbstractBendType
    refIndexIn :: Float64
    refIndexOut :: Float64
end

struct CDiffuser <: AbstractBendType
    tanθ::Float64
    refIndexIn::Float64
    refIndexOut::Float64
end

struct NoBendIndex <: AbstractBendType
    refIndexIn :: Float64
    refIndexOut::Float64
end

function NoBendIndex(n::Float64)
    #NoBendIndex(Float64(n), Float64(n))
    NoBendIndex(n, n)
end

struct NoAmpParam <: AbstractAmplitudeParam
    type::String
end

abstract type AbstractSurface <: GeometryBasics.GeometryPrimitive{3, Float64} end
#abstract type AbstractSurface <: GeometryBasics.AbstractGeometry{3, Float64} end

struct OptSurface <: AbstractSurface  #use the data to overload GemoetryBasics
    surfname::String
    base::SurfBase
    aperture::AbstractSize
    profile::AbstractSurfProfile
    mod::AbstractBendType
    coating::AbstractAmplitudeParam
    toGlobalCoord::AffineMap
    toLocalCoord::AffineMap
    toGlobalDir::LinearMap
    toLocalDir::LinearMap
end

struct ModelSurface <: AbstractSurface  #use the data to overload GemoetryBasics
    surfname::String
    base::SurfBase
    aperture::AbstractSize
    profile::AbstractSurfProfile
    refIndex::Float64 #needed for OPD calculations
    toGlobalCoord::AffineMap
    toLocalCoord::AffineMap
    toGlobalDir::LinearMap
    toLocalDir::LinearMap
end



"""
    ExtendedGeometry(geo, funcGeo, funcSetup, surfaceObject, wavelenth, parameters)

    Struct of the information to full characterize an optical system
    including pointers to functions that create the geometry and setup the values
    It includes a dictionary for definitions that can be used by the functions
        geo         array of AbstractSurfaces, the geometry. typically used for analysis
        funcGeo     the function that generates geo. enables adjustable parameters
        funcSetup   function called to set up conditions for funcGeo
        surfaceObject   define an Object to analyze
        wavelength  array of wavelengths in microns
        parameters  dictionary for parameters used by funcSetup and funcGeo

"""
struct ExtendedGeometry
    geo::Array{AbstractSurface}  #needs to be changed to AbstractOpticalObject
    funcGeo
    funcSetup
    surfaceObject::AbstractSurface
    wavelength::Array{Float64}
    parameters::Dict
end

