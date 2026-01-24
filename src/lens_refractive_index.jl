#=

Functions to get refractive index from the RefractiveIndex.info website

=#


export getRefractiveIndexFunc, setDirectoryBaseRefractiveIndex, loadRICatalog, getDirectoryBaseRI
export loadRICatalog!

global dirBaseRefractiveIndex::String ="/Users/matt/Development/Projects/refractiveindex/database/data"

function setDirectoryBaseRefractiveIndex()
    global dirBaseRefractiveIndex="/Users/matt/Development/Projects/refractiveindex/database/data"
end

function setDirectoryBaseRefractiveIndex(dir::AbstractString)
    global dirBaseRefractiveIndex= dir
end

getDirectoryBaseRI() = dirBaseRefractiveIndex

global defaultGlassCatalog = Dict{AbstractString, Any}()
defaultGlassCatalog["DEFAULT"]= rInDef

#=

fusedSilica = (0.696166300, 0.407942600, 0.897479400, 4.67914826E-3,
        1.35120631E-2, 97.9340025)

b270 = (93.2901, 11.7559, 1.27686, 206.221, -25.1637, 0.010476)

calciumFluoride = (0.443749998, 0.444930066, 0.150133991,
        0.00178027854, 0.00788536061, 0.0124119491)

magnesiumFluoride = (0.27620, 0.60967, 0.0080, 2.14973, 0.08636^2, 18.0^2, 25.0^2)

dlak6 = (2.81430513, -0.013449512, 0.0196569671, 0.000249252546, 1.98274573e-05, -8.02641015e-07)


libRI = YAML.load_file("/Users/matt/Development/Projects/refractiveIndex/database/library.yml")
=#


findRefractiveIndex(lambda, (b1, b2, b3, c1, c2, c3)) =
    sqrt(1 + (b1 * lambda^2)/(lambda^2 - c1) +
    (b2 *lambda^2)/(lambda^2 - c2) + (b3 *lambda^2)/(lambda^2 - c3))

findRefractiveIndexAlt(lambda, (c1, c2, c3, c4, c5, c6))=
    sqrt(c1 + c2*lambda^2 + c3 * lambda^-2 +c4 * lambda^-4 +
        c5 * lambda^-6+ c6*lambda^-8)


function riFormula2(λ, c)
    λ2 = λ^2
    sum = 1.0 + c[1]
    for i in 2:2:length(c)
        sum += c[i] * λ2/(λ2 - c[i+1])
    end
    sqrt(sum)
end

function riFormula1(λ, c)
    λ2 = λ^2
    sum = 1.0 + c[1]
    for i in 2:2:length(c)
        sum += c[i] * λ2/(λ2 - c[i+1]^2)
    end
    sqrt(sum)
end

function riFormula3(λ, c)

    sum = c[1]
    for i in 2:2:length(c)
        sum += c[i] * λ ^ c[i+1]
    end
    sqrt(sum)
end

#=
riN_LAK22(l) = riFormula2(l,
        map(x->parse(Float64,x),
            split(YAML.load_file("/Users/matt/Development/Projects/refractiveIndex/N-LAK22.yml")["DATA"][1]["coefficients"]," ")
            )
        )

riN_LAK22(0.5)
=#



function getRefractiveIndexFunc(basepath::AbstractString, path::AbstractString)
    funcs = Dict("formula 2"=>riFormula2, "formula 1"=> riFormula1, "formula 3"=>riFormula3)
    pth = joinpath(basepath,path)
    record = YAML.load_file(pth)
    data = record["DATA"][1]
    typeData = data["type"]
    coefentries = get(data, "coefficients", nothing)
    if coefentries == nothing
        return nothing
    end
    coefs = map(x->parse(Float64,x),split(coefentries))
    f(x) = funcs[typeData](x, coefs)
    f
end

function getRefractiveIndexFunc(path::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)
    getRefractiveIndexFunc(basepath, path)
end


#riN_LAK22 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-LAK22.yml")
#riN_SF6 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-SF6.yml")
#riN_SF2 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-SF2.yml")
#riN_BK7 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-BK7.yml")


function loadRICatalog(pth::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)
    dictforpath = Dict{AbstractString, Any}()
    dictforpath["DEFAULT"]= rInDef
    loadRICatalog!(dictforpath, pth::AbstractString; basepath)
end

function loadRICatalog!(dictforpath::Dict{AbstractString, Any}, pth::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)
    path = joinpath(basepath,pth)

    for (path, dirs, files) in walkdir(path)
        #println("Files in $path")
        for file in files
            #println(file) # path to files
            glassname, ext = splitext(file)
            if ext == ".yml"
                riFunc = getRefractiveIndexFunc(path, file)
                if riFunc==nothing
                    println("Skipping $(joinpath(path,file)): no dispersion function found")
                    continue
                end
                dictforpath[glassname] = riFunc
                sfile = split(glassname, "-")
                if length(sfile)==2
                    dictforpath[sfile[2]] = riFunc #remove the manufacturer specific prefix and add the glass
                end

            else
                println("Additional file found: $(joinpath(path,file))")
            end

        end
    end
    return dictforpath
end

function loadRICatalog!(pth::AbstractString; basepath::AbstractString = dirBaseRefractiveIndex)
    loadRICatalog!(defaultGlassCatalog, pth; basepath)
end
