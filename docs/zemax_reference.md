# Zemax import reference

This file collects two things that used to live as comments inside
`src/zemax.jl`:

1. The docstrings now attached to each type/function in `src/zemax.jl`
   (reproduced here for a single-page overview).
2. The original reference Python code those docstrings/implementations
   were informed by (or, for the `.zar` reader, code that has not been
   ported to Julia at all yet).

## Docstrings

### `ZemaxGeometry{N, T}`

Container type intended to hold a fully-imported Zemax system: the traced
geometry, the base point/direction the geometry starts from, and the
system's wavelength/name/units metadata.

Fields:
- `geo::Vector{AbstractSurface{N,T}}` — the imported optical surfaces
- `basept::Point{N, T}` — global coordinate the geometry starts at
- `dir::Vec{N, T}` — propagation direction the geometry starts along
- `wavelengths::Vector{T}` — wavelengths defined in the Zemax file
- `name::String` — system name, from the Zemax file's NAME field
- `units::String` — length units, from the Zemax file's UNIT field (e.g. "MM")

Not currently constructed anywhere in `zemax.jl` — `readZemax` returns
its parsed data as a plain tuple rather than wrapping it in a
`ZemaxGeometry`.

### `ZemaxSurf{T}`

Mutable intermediate representation of one Zemax `SURF` block, as parsed
from a `.zmx` file by `readZemax` before being converted into an
`OptSurface` by `zemaxsurfToSurface`.

Fields:
- `curvature::T` — surface curvature, from CURV
- `distance::T` — thickness to the next surface, from DISZ
- `material::String` — glass/material name, from GLAS (default "DEFAULT")
- `radius::T` — surface (semi-)diameter, from DIAM
- `stop::Bool` — whether this surface is the aperture stop, from STOP
- `conic::T` — conic constant, from CONI
- `aspherics::Vector{T}` — aspheric coefficients, indexed by Zemax PARM number
- `coating::String` — coating name, from COAT
- `type::String` — Zemax surface type, e.g. "STANDARD" or "EVENASPH", from TYPE
- `comm::String` — surface comment, from COMM

### `ZemaxSurf()`

Construct a `ZemaxSurf{Float64}` with default values: zero curvature,
distance, and conic; `"DEFAULT"` material; not a stop; `parmlength`
zeroed aspheric coefficients; empty coating; `"STANDARD"` type; empty
comment.

### `readZemax(filename::String; basept = ORIGIN, dir = ZAXIS)`

Read a Zemax `.zmx` sequential-lens text file and parse it into a vector
of `ZemaxSurf` records, one per `SURF` block (plus a leading record for
the object surface).

Recognized line keywords: `UNIT`, `NAME`, `WAVM`, `CURV`, `DISZ`, `GLAS`,
`DIAM`, `STOP`, `TYPE`, `CONI`, `PARM`, `COAT`, `COMM`. Any other
keyword is silently ignored. `basept`/`dir` are accepted but not
currently used during parsing itself — surface positions are computed
later, by `zemaxsurfToSurface`/`zemaxsurfsToGeo`.

Returns `(zsurfs, name, units, wavelengths)`:
- `zsurfs::Vector{ZemaxSurf}` — the parsed surfaces, in file order
- `name::String` — system name, from the NAME field
- `units::String` — length units, from the UNIT field
- `wavelengths::Vector{Float64}` — fixed-length (24-element) array indexed by Zemax wavelength number, from WAVM

### `printZemaxSurfs(zsurfs::Vector{ZemaxSurf})`

Print a one-line summary (curvature, distance, material, radius, stop,
conic, coating, type) for each surface in `zsurfs` to stdout, followed
by its aspheric coefficient vector. Useful for inspecting the result of
`readZemax`.

### `zemaxsurfToSurface(num, rinIn::Float64, rinOut::Float64, basept::Point, dir::Vec3, s::ZemaxSurf)`

Convert one parsed `ZemaxSurf` record into an `OptSurface` located at
`basept`, oriented along `dir`, with entering/exiting refractive indices
`rinIn`/`rinOut`. `num` is only used to build the surface's name
(`s.comm * " $num"`).

Supports Zemax `type`s `"STANDARD"` (built as a `SurfProfileConic`) and
`"EVENASPH"` (built as a `SurfProfileEvenAsphere`, using
`s.aspherics[2:end]`); any other type throws an error.

Returns `(newbasept, rinOut, newsurf)`, where `newbasept = basept +
s.distance * dir` is the starting point for the next surface, and
`rinOut` is passed through unchanged so it can be reused as the next
surface's `rinIn`.

### `zemaxsurfsToGeo(zemaxsurfs, base, dir, wavelength::Float64; glassCatalog = defaultGlassCatalog)`

Convert a vector of `ZemaxSurf` records (as returned by `readZemax`)
into a traceable geometry at the given wavelength, by repeatedly calling
`zemaxsurfToSurface` and threading the base point and refractive index
from one surface to the next.

`glassCatalog` maps each surface's `material` name to a function of
wavelength returning its refractive index (see `defaultGlassCatalog`).

Returns the resulting `Vector{AbstractSurface}`.

### `viewZemaxFile(filename; glassCatalog = defaultGlassCatalog)`

Read, convert, print, and plot a Zemax `.zmx` file in one call: reads
`filename` with `readZemax`, prints the parsed surfaces with
`printZemaxSurfs`, builds a geometry with `zemaxsurfsToGeo` (skipping the
leading object surface, `zgeo[2:end]`) at a fixed wavelength of 0.5,
displays it with `plotGeometry3D`, and prints it with `printGeo`.

Returns the built geometry (`Vector{AbstractSurface}`).

## Python reference code

### `zmx_to_system` — reference `.zmx` parser

The general structure of `readZemax`'s line-by-line parsing (recognizing
`SURF`, `CURV`, `DISZ`, `GLAS`, `DIAM`, `STOP`, `WAVL`, `COAT`, `CONI`,
`PARM`, and a long list of ignored fields) follows this Python reference
implementation.

```python
def zmx_to_system(data, item=None):
    s = System()
    next_pos = 0.
    s.append(Spheroid(material=air))
    for line in data.splitlines():
        e = s[-1]
        if not line.strip():
            continue
        line = line.strip().split(" ", 1)
        cmd = line[0]
        args = len(line) == 2 and line[1] or ""
        if cmd == "UNIT":
            s.scale = {
                    "MM": 1e-3,
                    "INCH": 25.4e-3,
                    "IN": 25.4e-3,
                    }[args.split()[0]]
        elif cmd == "NAME":
            s.description = args.strip("\"")
        elif cmd == "SURF":
            s.append(Spheroid(distance=next_pos, material=air))
        elif cmd == "CURV":
            e.curvature = float(args.split()[0])
        elif cmd == "DISZ":
            next_pos = float(args)
        elif cmd == "GLAS":
            args = args.split()
            name = args[0]
            try:
                e.material = Material.make(name)
            except KeyError:
                try:
                    e.material = Material.make((float(args[3]),
                                                float(args[4])))
                except Exception as e:
                    print("material not found", name, e)
        elif cmd == "DIAM":
            e.radius = float(args.split()[0])
        elif cmd == "STOP":
            e.stop = True
        elif cmd == "WAVL":
            s.wavelengths = [float(i)*1e-6 for i in args.split() if i]
        elif cmd == "COAT":
            e.coating = args.split()[0]
        elif cmd == "CONI":
            e.conic = float(args.split()[0])
        elif cmd == "PARM":
            i, j = args.split()
            i = int(i) - 1
            j = float(j)
            if i < 0:
                if j:
                    print("aspheric 0 degree not supported", cmd, args)
                continue
            if e.aspherics is None:
                e.aspherics = []
            while len(e.aspherics) <= i:
                e.aspherics.append(0.)
            e.aspherics[i] = j
        elif cmd in ("GCAT",  # glass catalog names
                     "OPDX",  # opd
                     "RAIM",  # ray aiming
                     "CONF",  # configurations
                     "ENPD", "PUPD",  # pupil
                     "EFFL",  # focal lengths
                     "VERS",  # version
                     "MODE",  # mode
                     "NOTE",  # note
                     "TYPE",  # surface type
                     "HIDE",  # surface hide
                     "MIRR",  # surface is mirror
                     "PARM",  # aspheric parameters
                     "SQAP",  # square aperture?
                     "XDAT", "YDAT",  # xy toroidal data
                     "OBNA",  # object na
                     "PKUP",  # pickup
                     "MAZH", "CLAP", "PPAR", "VPAR", "EDGE", "VCON",
                     "UDAD", "USAP", "TOLE", "PFIL", "TCED", "FNUM",
                     "TOL", "MNUM", "MOFF", "FTYP", "SDMA", "GFAC",
                     "PUSH", "PICB", "ROPD", "PWAV", "POLS", "GLRS",
                     "BLNK", "COFN", "NSCD", "GSTD", "DMFS", "ISNA",
                     "VDSZ", "ENVD", "ZVDX", "ZVDY", "ZVCX", "ZVCY",
                     "ZVAN", "XFLN", "YFLN", "VDXN", "VDYN", "VCXN",
                     "VCYN", "VANN", "FWGT", "FWGN", "WWGT", "WWGN",
                     "WAVN", "WAVM", "XFLD", "YFLD", "MNCA", "MNEA",
                     "MNCG", "MNEG", "MXCA", "MXCG", "RGLA", "TRAC",
                     "FLAP", "TCMM", "FLOA", "PMAG", "TOTR", "SLAB",
                     "POPS", "COMM", "PZUP", "LANG", "FIMP",
                     ):
            pass
        else:
            print(cmd, "not handled", args)
            continue
    return s
```

### `.zar` archive reading — not yet ported to Julia

Zemax can bundle a `.zmx` file with its supporting data into a `.zar`
archive. Nothing in `zemax.jl` reads `.zar` files yet; this Python
reference implementation (LZW-decompressing the archive's packed
entries) is kept here as the starting point for a future
`readZemaxArchive`-style function.

```python
import logging
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Generator, List, Union

log = logging.getLogger(__name__)

__all__ = ['read', 'UnpackedData', 'extract', 'repack']

ZAR = '.zar'
ZIP = '.zip'
ZAR_VERSION_LENGTH = 2  # in bytes
EARLIER_CONTENT_OFFSET = 0x14C - ZAR_VERSION_LENGTH
EARLIER_PACKED_FILE_SIZE_BEGIN = 0xC - ZAR_VERSION_LENGTH
EARLIER_PACKED_FILE_SIZE_END = 0x10 - ZAR_VERSION_LENGTH
EARLIER_PACKED_FILE_NAME_OFFSET = 0x20 - ZAR_VERSION_LENGTH
EARLIER_VERSION = 0xEA00.to_bytes(2, 'big')
LATEST_VERSION = 0xEC03.to_bytes(2, 'big')
LATEST_CONTENT_OFFSET = 0x288 - ZAR_VERSION_LENGTH
LATEST_PACKED_FILE_SIZE_BEGIN = 0x10 - ZAR_VERSION_LENGTH
LATEST_PACKED_FILE_SIZE_END = 0x18 - ZAR_VERSION_LENGTH
LATEST_PACKED_FILE_NAME_OFFSET = 0x30 - ZAR_VERSION_LENGTH


def _decompress_lzw(compressed: bytes) -> bytes:
    """
    Decompresses bytes using the variable LZW algorithm, starting with code strings of length 9.

    This function is used internally by the read function.
    General information about LZW: https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv%E2%80%93Welch
    Adapted partly from https://gist.github.com/BertrandBordage/611a915e034c47aa5d38911fc0bc7df9

    :param compressed: The compressed bytes without header.
    :return: The decompressed bytes.
    """
    # Convert input to bits
    compressed_bits: str = bin(int.from_bytes(compressed, 'big'))[2:].zfill(len(compressed) * 8)
    # convert to binary string and pad to 8-fold length

    code_word_length = 8
    words: List[bytes] = [_.to_bytes(1, 'big') for _ in range(2**code_word_length)]
    # integer codes refer to a words in an expanding dictionary

    bit_index = 0
    previous_word: bytes = b''
    decompressed: List[bytes] = []

    while True:
        if 2**code_word_length <= len(words):  # If the dictionary is full
            code_word_length += 1              # increase the code word length
        if bit_index + code_word_length > len(compressed_bits):
            break  # stop when the bits run out
        # Get the next code word from the data bit string
        code = int(compressed_bits[bit_index:bit_index + code_word_length], 2)
        bit_index += code_word_length

        # If word in dictionary, use it; else add it as a new word
        latest_word: bytes = words[code] if code < len(words) else previous_word + previous_word[:1]
        decompressed.append(latest_word)  # Update result
        if len(previous_word) > 0:  # Skip first iteration
            words.append(previous_word + latest_word[:1])  # Add as new encoding

        previous_word = latest_word

    return b''.join(decompressed)  # convert to bytes


@dataclass
class UnpackedData(object):
    """A structure to represent the file blocks in a zar-archive.

    Parameters:
        name: A string with the name of the file contained in the archive.
        unpacked_contents: The unpacked (decompressed) bytes of this file.
    """

    file_name: str
    unpacked_contents: bytes


def read(input_full_file: Union[Path, str]) -> Generator[UnpackedData, None, None]:
    """
    Reads a zar archive file and generates a series of (unpacked file name, unpacked file contents) tuples.

    The returned Generator produces tuples in the order found in the archive.

    :param input_full_file: The archive or the path to the archive.
    :return: A Generator of name-data tuples.
    """
    # Make sure that the input arguments are both pathlib.Path-s
    if isinstance(input_full_file, str):
        input_full_file = Path(input_full_file.strip())
    with open(input_full_file, 'rb') as input_file:
        while True:
            version = input_file.read(ZAR_VERSION_LENGTH)
            if len(version) < ZAR_VERSION_LENGTH:
                break  # end of file
            if version[0] == LATEST_VERSION[0]:
                header_length = LATEST_CONTENT_OFFSET
            elif version[0] == EARLIER_VERSION[0]:
                header_length = EARLIER_CONTENT_OFFSET
            else:
                log.warning(f'Unknown ZAR header "{version.hex()}"!')
                header_length = LATEST_CONTENT_OFFSET
                version = LATEST_VERSION  # override and cross fingers

            header = input_file.read(header_length)

            if version[0] == LATEST_VERSION[0]:
                packed_file_size = int.from_bytes(
                    header[LATEST_PACKED_FILE_SIZE_BEGIN:LATEST_PACKED_FILE_SIZE_END],
                    byteorder='little',
                    signed=False,
                )
                packed_file_name = header[LATEST_PACKED_FILE_NAME_OFFSET:].decode('utf-16-le')
                packed_file_name = packed_file_name[:packed_file_name.find('\0')]  # ignore all 0's on the right
            else:
                packed_file_size = int.from_bytes(
                    header[EARLIER_PACKED_FILE_SIZE_BEGIN:EARLIER_PACKED_FILE_SIZE_END],
                    byteorder='little',
                    signed=False,
                )
                packed_file_name_bytes = header[EARLIER_PACKED_FILE_NAME_OFFSET:]
                packed_file_name_bytes = packed_file_name_bytes[:packed_file_name_bytes.find(0x0)]
                packed_file_name = packed_file_name_bytes.decode('utf-8')
            log.debug(f'Version {version.hex()}. Packed file {packed_file_name} has size {packed_file_size} bytes.')

            # Read and process data
            archive_data = input_file.read(packed_file_size)
            if packed_file_name[-4:].upper() == '.LZW':
                archive_data = _decompress_lzw(archive_data)
                packed_file_name = packed_file_name[:-4]

            # Yield a series of tuples from the Generator
            yield UnpackedData(file_name=packed_file_name, unpacked_contents=archive_data)


def extract(input_full_file: Union[Path, str], output_path: Union[Path, str, None] = None) -> None:
    """
    Imports the data from a zar archive file and writes it as a regular directory.

    :param input_full_file: The path to zar-file.
    :param output_path: The path where the files should be saved. Default: the same as the input_full_file but
        without the extension.
    """
    # Make sure that the input arguments are both pathlib.Path-s
    if isinstance(input_full_file, str):
        input_full_file = Path(input_full_file.strip())
    if output_path is None:  # By default, just drop the .zar extension for the output names
        output_path = input_full_file.parent / (
            input_full_file.stem if input_full_file.suffix.lower() == ZAR else input_full_file
        )
    elif isinstance(output_path, str):
        output_path = Path(output_path.strip())
    Path.mkdir(output_path, exist_ok=True, parents=True)
    log.debug(f'Extracting {input_full_file} to directory {output_path}/...')

    # Unpack and store the recovered data
    for unpacked_data in read(input_full_file):
        with open(output_path / unpacked_data.file_name, 'wb') as unpacked_file:
            unpacked_file.write(unpacked_data.unpacked_contents)

    log.info(f'Extracted {input_full_file} to directory {output_path}/.')


def repack(input_full_file: Union[Path, str], output_full_file: Union[Path, str, None] = None) -> None:
    """
    Imports the data from a zar archive file and writes it as a regular zip file.

    :param input_full_file: The file path, including the file name, of the zar-file.
    :param output_full_file: TThe file path, including the file name, of the destination zip-file.
        Default: the same as `input_full_file` but with the extension changed to 'zip'.
    """
    # Make sure that the input arguments are both pathlib.Path-s
    if isinstance(input_full_file, str):
        input_full_file = Path(input_full_file.strip())
    if output_full_file is None:  # By default, just change .zar to .zip for the destination archive
        if input_full_file.suffix.lower() == ZAR:
            output_full_file = input_full_file.with_suffix(ZIP)
        else:  # or tag on .zip when it hasn't the .zar extension
            output_full_file = input_full_file.parent / (input_full_file.name + ZIP)
    else:
        if isinstance(output_full_file, str):
            if not output_full_file.lower().endswith(ZIP):
                output_full_file += '/' + input_full_file.name + ZIP
            output_full_file = Path(output_full_file.strip())
        elif isinstance(output_full_file, Path) and not output_full_file.name.lower().endswith(ZIP):
            output_full_file /= input_full_file.name + ZIP
        Path.mkdir(output_full_file.parent, exist_ok=True, parents=True)
    log.debug(f'Converting {input_full_file} to zip archive {output_full_file}...')

    # Open the output archive and start storing unpacked files
    repack_directory = output_full_file.stem  # all but the extension
    with zipfile.ZipFile(
        output_full_file,
        mode='a',
        compression=zipfile.ZIP_DEFLATED,
        allowZip64=False,
        compresslevel=9,
    ) as archive_file:
        # Unpack and store the recovered data
        for unpacked_data in read(input_full_file):
            archive_file.writestr(f'{repack_directory}/{unpacked_data.file_name}', unpacked_data.unpacked_contents)

    log.info(f'Converted {input_full_file} to zip archive {output_full_file}.')
```
