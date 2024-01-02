#extended geometry functions
export findPerpenMap, updateCoordChange!
"""
    findPerpenMap(planenormal)
Find a perpendicular to plane normal - "local x direction" by using y dir as the
guess. 
"""
function findPerpenMap(planenormal::SVector{3})
    tryy = YAXIS
    newx = cross(tryy, planenormal) # if planenormal == ZAXIS this give x axis
    if norm(newx)==0. # planenormal is along y axis
        tryy = XAXIS
        newx = cross(tryy, planenormal)
    end

    newx=normalize(newx)
    newy = SVector{3}(cross(planenormal, newx)...)
    newy, LinearMap(SMatrix{3,3}([newx; newy; planenormal]))
end

"""
    findPerpenMap(planenormal, ydir)
Find a perpendicular in the plane to planenormal and ydir
"""
function findPerpenMap(planenormal::SVector{3}, ydir::Union{Nothing, SVector{3}})
    if isnothing(ydir)
        tryy = YAXIS
        newx = cross(tryy, planenormal) # if planenormal == ZAXIS this give x axis
        if norm(newx) < 0.1 # planenormal is very close to Y axis
            tryy = XAXIS
            newx = cross(tryy, planenormal)
        end
        newx = normalize(newx)
        newy = SVector{3}(cross(planenormal, newx)...)
    else
        newx = cross(ydir, planenormal) # if planenormal == ZAXIS this give x axis
        newx = normalize(newx)
        newy = ydir
    end
    newy, LinearMap(SMatrix{3,3}([newx; newy; planenormal]))
end



"""
    updateCoordChange()
update the coordinate transformation functions when the location and/or
direction changes
"""
function updateCoordChange(pointInPlane::SVector{3}, planenormal::SVector{3})
    pip = pointInPlane
    ydir, toGlobalDir = findPerpenMap(planenormal)
    toGlobalCoord = compose(Translation(pip), toGlobalDir)
    toLocalCoord = inv(toGlobalCoord)
    toLocalDir = inv(toGlobalDir)
    id = hcat(XAXIS, YAXIS, ZAXIS)
    toLCMat =SMatrix{3,3,Float64,9}(toLocalDir(id))
    toGCMat = SMatrix{3,3,Float64,9}(toGlobalDir(id))
    return ydir,toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir#, toGCMat, toLCMat
end

function updateCoordChange(pointInPlane::SVector{3}, planenormal::SVector{3}, ydir::Union{SVector{3},Nothing})
    pip = pointInPlane
    myydir, toGlobalDir = findPerpenMap(planenormal, ydir)
    toGlobalCoord = compose(Translation(pip), toGlobalDir)
    toLocalCoord = inv(toGlobalCoord)
    toLocalDir = inv(toGlobalDir)
    id = hcat(XAXIS, YAXIS, ZAXIS)
    toLCMat =SMatrix{3,3,Float64,9}(toLocalDir(id))
    toGCMat = SMatrix{3,3,Float64,9}(toGlobalDir(id))
    return myydir,toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir#, toGCMat, toLCMat
end

function updateCoordChange!(surf)
    -, surf.toGlobalCoord, surf.toLocalCoord, surf.toGlobalDir, surf.toLocalDir = updateCoordChange(surf.base.base, surf.base.dir, surf.base.ydir)
    #, surf.toLCMat, surf.toCGMat
end

"""
    EGeo(func, baseObject, size, wavelength; dir = ZAXIS, setup = defaultSetupGeo, parameters = Dict{Any,Any}() )
    simple way to set up ExtendedGeometry

    func - function defining geometry
    baseObject - location of object
    size - size of object
    wavelength - array of wavelengths in microns

    dir - direction of object
    setup - function to set parameters for func
    parameters - dictionary for parameters used by funcSetup and funcGeo

"""
function EGeo(func, baseObject, size, wavelength; dir = ZAXIS, setup = defaultSetupGeo, parameters = Dict{Any,Any}() )
    object = referencePlane("object", baseObject, dir, refIndexDefault, size, "none")
    newgeo = setup(func, object, wavelength, parameters)
    a = ExtendedGeometry(newgeo, func, setup, object, wavelength, parameters)
end


"""
    defaultSetupGeo(func, object, wavelength, parameters; numWL = 1)
    predefined function for simplest setup: nothing but running func

"""
function defaultSetupGeo(func, object, wavelength, parameters; numWL = 1)
    # perform setup functions for the GeometryBasics
    geo = func(parameters, wavelength[numWL])
end


"""
    updateEGeo!(egeo; numWL=1) calls the function necessary to create a static geometry for tracing
        egeo    extended geometry to update
        numWL   which wavelength to use
"""
function updateEGeo!(egeo::ExtendedGeometry; numWL = 1)
    defaultSetupGeo(egeo.funcGeo, egeo.surfaceObject, egeo.wavelength, egeo.parameters; numWL = numWL)
end


#=
    ampltude

=#
