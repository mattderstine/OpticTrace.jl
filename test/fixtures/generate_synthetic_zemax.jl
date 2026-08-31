#=
    Maintenance script (not part of the test suite -- not `include`d by
    runtests.jl) that (re)generates the small synthetic .zar/.zmf fixtures
    checked into this directory: test_archive.zar and test_catalog.zmf.

    Binary archive formats aren't hand-editable the way test_singlet.zmx
    is, so this readable/diffable generator is the auditable source of
    truth for those fixtures -- run it manually if the .zar/.zmf format
    understanding in docs/zemax_reference.md / src/zemax.jl ever changes:

        julia --project=. test/fixtures/generate_synthetic_zemax.jl

    Byte layout below matches src/zemax.jl's readZemaxArchive (the
    "earlier"/0xEA header) and readZmfCatalog exactly -- see those
    functions' docstrings for the field-by-field reference this was
    built from.
=#
using OpticTrace

function buildZarEntryBytes(name::String, payload::Vector{UInt8})
    io = IOBuffer()
    write(io, UInt8(0xEA)) # version byte (checked: "earlier" header layout)
    write(io, UInt8(0x00)) # second version-tag byte (not inspected by the reader)

    header = zeros(UInt8, 330) # 0x14C - 2
    sizeBytes = reinterpret(UInt8, [UInt32(length(payload))])
    header[11:14] .= sizeBytes # (0xC-2+1):(0x10-2), packed byte count

    nameBytes = Vector{UInt8}(name)
    header[31:30+length(nameBytes)] .= nameBytes # (0x20-2)+1 onward, null-terminated by the zero padding

    write(io, header)
    write(io, payload)
    return take!(io)
end

function generateSyntheticZar(path::String)
    zmxText = """
    VERS 100414 0 26970
    NAME Synthetic Archive Lens
    UNIT MM
    SURF 0
    SURF 1
      TYPE STANDARD
      CURV 0.05
      DISZ 5.0
      GLAS TESTGLASS 0 0 1.5
      DIAM 10.0
    """
    entry1 = buildZarEntryBytes("SYNTH.ZMX", Vector{UInt8}(zmxText))
    entry2 = buildZarEntryBytes("SYNTH.AGF", Vector{UInt8}("synthetic glass catalog stub, not real glass data\n"))

    io = IOBuffer()
    write(io, entry1)
    write(io, entry2)
    write(path, take!(io))
end

function buildZmfEntryBytes(name::String, elements::Int, efl::Float64, enp::Float64, plaintext::Vector{UInt8})
    io = IOBuffer()
    nameBytes = zeros(UInt8, 100)
    nb = Vector{UInt8}(name)
    nameBytes[1:length(nb)] .= nb
    write(io, nameBytes)

    write(io, UInt32(0))         # per-lens format version -- not surfaced/checked
    write(io, UInt32(elements))
    write(io, UInt32(0))         # shape code -- not surfaced/checked
    write(io, UInt32(0))         # aspheric flag -- not surfaced/checked
    write(io, UInt32(0))         # grin flag -- not surfaced/checked
    write(io, UInt32(0))         # toroidal flag -- not surfaced/checked
    write(io, UInt32(length(plaintext)))
    write(io, Float64(efl))
    write(io, Float64(enp))

    obfuscated = OpticTrace.zmfDeobfuscate(plaintext, efl, enp) # self-inverse: plaintext -> "obfuscated"
    write(io, obfuscated)
    return take!(io)
end

function generateSyntheticZmf(path::String)
    lens1 = Vector{UInt8}("VERS 100503\nMODE SEQ\nNAME Synthetic Lens One\nUNIT MM\n")
    lens2 = Vector{UInt8}("VERS 100503\nMODE SEQ\nNAME Synthetic Lens Two\nUNIT MM\n")

    io = IOBuffer()
    write(io, UInt32(1001))
    write(io, buildZmfEntryBytes("LENS1", 1, 4.485, 5.2, lens1))
    write(io, buildZmfEntryBytes("LENS2", 2, 10.0, 8.0, lens2))
    write(path, take!(io))
end

if abspath(PROGRAM_FILE) == @__FILE__
    generateSyntheticZar(joinpath(@__DIR__, "test_archive.zar"))
    generateSyntheticZmf(joinpath(@__DIR__, "test_catalog.zmf"))
    println("Wrote test_archive.zar and test_catalog.zmf")
end
