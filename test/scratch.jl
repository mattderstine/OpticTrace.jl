4+4

using OpticTrace

#=

attributesSurfaces = Dict{String, Any}()


=#

struct AmpParaInterpR <: AbstractAmplitudeParam
    type::String #label the surface coating
    funcScatter #the deflection function including parameters, returns delta vector
    funcReflect #the reflection function, returns Float64 reflectivity
end

function surfAmpFunc!(dirIn::Vec3, dirOut::Vec3, normal::Vec3, newRayBase::Point3, ndex::CDiffuser,amp::AmpParaInterpR)
    

    AmpData(identity, identity, funcReflect(newRayBase))
end

