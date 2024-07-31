


export Ray, SurfBase, Trace, OptSurface, ModelSurface, ExtendedGeometry, SurfProfileOAConic, SizeLens, RoundAperture, RectAperture
export NoProfile, NoBendIndex, NoAmpParam, updateCoordChange, AbstractSurface, AbstractAmplitudeParam, DielectricT,MirrorR,CDiffuser, NoBendIndex
export AmpData


abstract type AbstractRay{N, T<:Real} end
abstract type AbstractSurfBase{N, T<:Real} end

abstract type AbstractSize{T<:Real} end
abstract type AbstractSurfProfile{T<:Real} end
abstract type AbstractAmplitudeParam end
abstract type AbstractBendType{T<:Real} end #some info on how to do the bending
abstract type AbstractBendDielectric{T} <: AbstractBendType{T} end
abstract type AbstractBendMirror{T} <: AbstractBendType{T} end
abstract type AbstractOpticalObject{T<:Real} end
abstract type AbstractTrace{T<:Real} end

APorString = Union{AbstractAmplitudeParam, String}

struct Ray{N, T} <: AbstractRay{N, T}
    base::Point{N, T}
    dir::Vec{N, T}

end

mutable struct SurfBase{N, T} <: AbstractSurfBase{N, T}
    base::Point{N, T}
    dir::Vec{N, T}
    ydir::Vec{N, T}
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
    offset::Vec{3, T}
end

struct SurfProfileAsphere{T} <: AbstractSurfProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ
    a::Vector{T} #hope this works...
end

struct SurfProfileCyl{T} <: AbstractSurfProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ    
#   ϵy::T #see Welford for definition of ϵ
    a::Vector{T} #hope this works...
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

abstract type AbstractAmpData{T <:Real} end

mutable struct AmpData{T} <: AbstractAmpData{T}
    p::SMatrix{3,3,T}
    o::SMatrix{3,3,T}
    trans::Array{T,1} #for intensity systems, the "transmission" coefficient of the surface
end


"""
    Trace(ray, nIn, delta, pmatrix)
    Element of a raytrace.
        ray     location of interesection with surface and exiting direction
        nIn     refractive index leaving surface
        delta   length of ray leaving surface
        pmatrix polarization p & o matrices
"""
mutable struct Trace{T <:Real, S <:AbstractAmpData{T}} 
    ray::Ray{3, T}
    nIn::T
    delta::T
    pmatrix::S
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

mutable struct OptSurface{N,T, S <:AbstractAmplitudeParam, U<:AbstractSize{T}, V<:AbstractSurfProfile{T}, W<:AbstractBendType{T}} <: AbstractSurface{N,T} #use the data to overload GemoetryBasics
    surfname::String
    base::SurfBase{N, T}
    aperture::U
    profile::V
    mod::W
    coating::S
    toGlobalCoord::AffineMap
    toLocalCoord::AffineMap
    toGlobalDir::LinearMap
    toLocalDir::LinearMap
    color
end

mutable struct ModelSurface{N,T, U<:AbstractSize{T}, V<:AbstractSurfProfile{T}, } <: AbstractSurface{N,T}  #use the data to overload GemoetryBasics
    surfname::String
    base::SurfBase{N, T}
    aperture::U
    profile::V
    refIndex::T #needed for OPD calculations
    toGlobalCoord::AffineMap
    toLocalCoord::AffineMap
    toGlobalDir::LinearMap
    toLocalDir::LinearMap
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

