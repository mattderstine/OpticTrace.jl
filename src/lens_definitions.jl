


export Ray, SurfBase, Trace, OptSurface, ModelSurface, ExtendedGeometry, SurfProfileOAConic, SizeLens, RoundAperture, RectAperture
export NoProfile, NoBendIndex, NoAmpParam, updateCoordChange, AbstractSurface, AbstractAmplitudeParam, DielectricT,MirrorR,CDiffuser, NoBendIndex
export AmpData


abstract type AbstractRay{T<:Real,N} end
abstract type AbstractSurfBase{T<:Real,N} end

abstract type AbstractSize{T<:Real} end
abstract type AbstractSurfProfile{T<:Real} end
abstract type AbstractAmplitudeParam end
abstract type AbstractBendType{T<:Real} end #some info on how to do the bending
abstract type AbstractBendDielectric{T<:Real} <: AbstractBendType{T} end
abstract type AbstractBendMirror{T<:Real} <: AbstractBendType{T} end
abstract type AbstractOpticalObject{T<:Real} end
abstract type AbstractTrace{T<:Real,N} end

APorString = Union{AbstractAmplitudeParam, String}

struct Ray{T,N} <: AbstractRay{T,N}
    base::SVector{N,T}
    dir::SVector{N,T}
end

#=
function Ray(b, d)
    Ray(Point{3}(b), Vec{3}(d))
end
=#

mutable struct SurfBase{T, N} <: AbstractSurfBase{T,N}
    base::SVector{N,T}
    dir::SVector{N,T}
    ydir::SVector{N,T}
end


struct SizeLens{T} <: AbstractSize{T}  #size for round optics
    semiDiameter::T
end

abstract type TrueAperture{T} <: AbstractSize{T} end

struct RoundAperture{T} <: TrueAperture{T}
    obscure::T
    semiDiameter::T
end

struct RectAperture{T} <: TrueAperture{T}
    wo::T
    lo::T
    wclear::T
    lclear::T
end


struct SurfProfileSphere{T} <: AbstractSurfProfile{T}
    curv :: T
end


struct SurfProfileConic{T} <: AbstractSurfProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ
end

struct SurfProfileOAConic{T} <: AbstractSurfProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ
    offset::SVector{3, T}
end

struct SurfProfileAsphere{T} <: AbstractSurfProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ
    a::Array{T} #hope this works...
end

struct SurfProfileCyl{T} <: AbstractSurfProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ    
#   ϵy::T #see Welford for definition of ϵ
    a::Array{T} #hope this works...
end

struct SurfProfileToroid{T} <: AbstractSurfProfile{T}
    curvY :: T
    curvX :: T 
end

struct NoProfile{T} <: AbstractSurfProfile{T}
    curv :: T #but never used
end


struct AmpParam <: AbstractAmplitudeParam
    type::String
end

abstract type AbstractAmpData{T <:Real, N} end

mutable struct AmpData{T,N,M, L} <: AbstractAmpData{T,N}
    p::SMatrix{N,M,T, L}
    o::SMatrix{N,M,T, L}
    trans::T #for intensity systems, the "transmission" coefficient of the surface
end


"""
    Trace(ray, nIn, delta, pmatrix)
    Element of a raytrace.
        ray     location of interesection with surface and exiting direction
        nIn     refractive index leaving surface
        delta   length of ray leaving surface
        pmatrix polarization p & o matrices
"""
mutable struct Trace{T,N} <:AbstractTrace{T,N}
    ray::Ray{T,N}
    nIn::T
    delta::T
    pmatrix::AbstractAmpData{T,N}
end

#=

AbstractBendTypes

=#


struct DielectricT{T} <: AbstractBendDielectric{T}
    refIndexIn :: T
    refIndexOut :: T
end

struct MirrorR{T} <: AbstractBendMirror{T}
    refIndexIn :: T
    refIndexOut :: T
end

struct CDiffuser{T} <: AbstractBendType{T}
    tanθ::T
    refIndexIn::T
    refIndexOut::T
end

struct NoBendIndex{T} <: AbstractBendType{T}
    refIndexIn :: T
    refIndexOut::T
end
NoBendIndex(n) = NoBendIndex(n,n)

#=
function NoBendIndex(n::Float64)
    #NoBendIndex(Float64(n), Float64(n))
    NoBendIndex(n, n)
end
=#

struct NoAmpParam <: AbstractAmplitudeParam
    type::String
end

abstract type AbstractSurface{N,T} <: GeometryBasics.GeometryPrimitive{N, T} end
#abstract type AbstractSurface <: GeometryBasics.AbstractGeometry{3, Float64} end

mutable struct OptSurface{N,T,L} <: AbstractSurface{N,T} #use the data to overload GemoetryBasics
    surfname::String
    base::SurfBase{T,N}
    aperture::AbstractSize{T}
    profile::AbstractSurfProfile{T}
    mod::AbstractBendType{T}
    coating::AbstractAmplitudeParam
    toGlobalCoord::AffineMap
    toLocalCoord::AffineMap
    toGlobalDir::LinearMap
    toLocalDir::LinearMap
    toGCMat::SMatrix{N,N,T,L}
    toLCMat::SMatrix{N,N,T,L}
    color
end

mutable struct ModelSurface{N,T,L} <: AbstractSurface{N,T}  #use the data to overload GemoetryBasics
    surfname::String
    base::SurfBase{T,N}
    aperture::AbstractSize{T}
    profile::AbstractSurfProfile{T}
    refIndex::T #needed for OPD calculations
    toGlobalCoord::AffineMap
    toLocalCoord::AffineMap
    toGlobalDir::LinearMap
    toLocalDir::LinearMap
    toGCMat::SMatrix{N,N,T,L}
    toLCMat::SMatrix{N,N,T,L}
    color
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

