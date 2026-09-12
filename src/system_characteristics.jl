#=

    # Functions to determine the basic characteristics of an optical system

    The parameters to be found are:
    - Paraxial(small angle ray) Focal length
    - Paraxial Entrance pupil size
    - Paraxial Entrance pupil position
    - F/# and NA
    - Paraxial Exit pupil size
    - Paraxial Exit pupil position
    - Distortion - Maybe
    - Vignetting - Maybe

    - detection of afocal and doubly telecentric systems
    - detection of afocal image spaces

=#

"""
    docstring needed

"""
mutable struct SystemCharacteritics{T<:Real}
    paraxialFocalLength::T
    fNumber::T
    numericalAperture::T
    entrancePupilSemiDiam::T
    entrancePupilPosition::Point{T, 3}
    exitPupilSemiDiam::T
    exitPupilPosition::Point{T, 3}
end

"""
    stub to be updated

    findSystemCharacteristics!(sc::SystemCharacteritics, sys::OpticalSystem)

    Find the basic characteristics of an optical system.

    # Arguments
    - `sc::SystemCharacteritics`: The structure to store the system characteristics.
    - `sys::OpticalSystem`: The optical system to analyze.

    # Returns
    - `SystemCharacteritics`: A structure containing the basic characteristics of the optical system.
"""
function findSystemCharacteristics!(sc::SystemCharacteritics, sys::OpticalSystem)
    #thse functions are placeholders. The actual implementation may be different.
    
    # Find the paraxial focal length
    sc.paraxialFocalLength = findParaxialFocalLength(sys)

    # Find the paraxial entrance pupil size and position
    sc.entrancePupilSemiDiam, sc.entrancePupilPosition = findParaxialEntrancePupil(sys)

    # Find the F/# and NA
    sc.fNumber, sc.numericalAperture = findFNumberAndNA(sys)

    # Find the paraxial exit pupil size and position
    sc.exitPupilSemiDiam, sc.exitPupilPosition = findParaxialExitPupil(sys)

    # Return a structureof the characteristics
    return sc
end