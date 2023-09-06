


export Ray, SurfBase, Trace, OptSurface, ModelSurface, ExtendedGeometry, SurfProfileOAConic

abstract type AbstractRay{T<:Real,N} end
abstract type AbstractSurfBase{T<:Real,N} end

struct Ray{T,N} <: AbstractRay{T,N}
    base::SVector{N, T}
    dir::SVector{N, T}
#    function Ray(b::Tuple{T, T, T}, d::Tuple{T, T, T}) where {T}
#        return new{T,3}(SVector(b), SVector(d))
#    end
end

struct SurfBase{T,N} <: AbstractSurfBase{T,N}
    base::SVector{N, T}
    dir::SVector{N, T}
    ydir::SVector{N, T}
end

abstract type AbstractSize{T<:Real} end
abstract type AbstractSurfProfile{T<:Real} end
abstract type AbstractAmplitudeParam end
abstract type AbstractBendType{T<:Real} end #some info on how to do the bending
abstract type AbstractOpticalObject{T<:Real} end
abstract type AbstractTrace{T<:Real,N} end

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
    offset::SVector{3,T}
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

abstract type AbstractAmpData{T <:Real} end

struct AmpData{T} <: AbstractAmpData{T}
    p::Matrix{T}
    o::Matrix{T}
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
struct Trace{T,N} <:AbstractTrace{T,N}
    ray::Ray{T,N}
    nIn::T
    delta::T
    pmatrix::AbstractAmpData{T}
end


struct DielectricT{T} <: AbstractBendType{T}
    refIndexIn :: T
    refIndexOut :: T
end

struct MirrorR{T} <: AbstractBendType{T}
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
    function NoBendIndex(n::T) where {T<:Real}
        return new{T}(n,n)
    end
end

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

struct OptSurface{N,T} <: AbstractSurface{N,T}  #use the data to overload GemoetryBasics
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
end

struct ModelSurface{N,T} <: AbstractSurface{N,T}  #use the data to overload GemoetryBasics
    surfname::String
    base::SurfBase{T,N}
    aperture::AbstractSize{T}
    profile::AbstractSurfProfile{T}
    refIndex::T #needed for OPD calculations
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

