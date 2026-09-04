# Zemax import reference

This file collects the original reference Python code
`src/zemax.jl`'s implementation was ported from/informed by: the
general structure of `readZemax`'s line-by-line `.zmx` parsing, plus
the separate Python references the `.zar`/`.zmf` archive-reading code
was ported from. Kept as historical context on the binary format
layouts (`.zar`/`.zmf`) and the original parsing approach -- not a
description of `src/zemax.jl`'s current API, which has grown
substantially since (`ZemaxHeader`, `OpticalSystem`, `readZemaxSystem`,
support for many more Zemax surface `TYPE`s, ...). For current API
documentation, see each type/function's own docstring in
`src/zemax.jl` (via Julia's `@doc` or the source directly) -- this
file previously also reproduced a snapshot of those docstrings, but
that snapshot was removed once it had drifted far enough from the
real API to be actively misleading rather than useful; `CLAUDE.md`'s
Project Structure section has a brief, current overview instead.

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

### `.zar` archive reading — ported to Julia

Zemax can bundle a `.zmx` file with its supporting data into a `.zar`
archive. `zemax.jl`'s `readZemaxArchive`/`listZemaxArchive`/
`extractZemaxArchive` and `lzwDecompress` are a direct port of this
Python reference implementation (LZW-decompressing the archive's packed
entries), kept here for reference. Both header layouts described below
("earlier"/`0xEA` and "latest"/`0xEC`) have been validated against real
sample `.zar` files.

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

### `.zmf` catalog reading — ported to Julia

Zemax lens vendors (Edmund Optics, Thorlabs, and many others) publish
whole families of stock lenses as a single `.zmf` lens-catalog file.
Unlike `.zar`, this format is not documented anywhere in Ansys/Zemax's
own materials; the layout below is an unofficial, community
reverse-engineering that several independent tools (this codebase
included, now) rely on. `zemax.jl`'s `readZmfCatalog`/
`listZmfCatalog`/`extractZmfCatalog` and `zmfDeobfuscate` are a direct
port of `rayopt`'s `zmf_read`/`zmf_obfuscate`
(https://github.com/quartiq/rayopt/blob/master/rayopt/zemax.py),
re-implemented without `rayopt`'s SQLAlchemy-backed `Catalog`/session
machinery (not needed here) and without the deprecated
`numpy.fromstring`/`.tostring()` calls the original uses (removed in
numpy >= 2.0).

Binary layout, all fields little-endian:

- A 4-byte `UInt32` file header: the catalog format version. Only
  `1001` has ever been observed/documented; `readZmfCatalog` throws an
  error on any other value rather than guessing.
- Then, repeated to end of file, one fixed 144-byte record per lens:
  - 100 bytes: the lens's catalog name/part number (NUL-padded, not
    obfuscated).
  - 7 × `UInt32`: per-lens format version (matches the `VERS` line in
    that lens's decoded description), element count, a shape-code
    index into `"?EBPM"`, and aspheric/GRIN/toroidal flags. None of
    these beyond element count are currently surfaced on `ZmfEntry`.
  - 2 × `Float64`: effective focal length (`efl`) and entrance pupil
    diameter (`enp`).
  - Immediately followed by that lens's *description*: `descLen` bytes
    (from one of the `UInt32` fields above) of XOR-obfuscated text
    which, once deobfuscated, is byte-for-byte the same `.zmx`-format
    grammar `readZemax` already parses (starting with a `VERS ######`
    line matching the record's version field).

The obfuscation keystream (`zmf_obfuscate` in the Python source below;
the same function both obfuscates and deobfuscates, since XOR is its
own inverse) is derived per output byte from a fixed trigonometric
formula seeded by that lens's own `efl`/`enp`, then reduced to a byte by
formatting the intermediate value in `%.8e` scientific notation and
taking 3 of its digit characters. This was validated empirically during
development (not just read off the reference source) by running a
from-scratch, dependency-free Python re-implementation against several
real vendor `.zmf` catalogs and checking that the decoded description's
own `VERS` line matched its record's version field exactly.

```python
from struct import Struct

head = Struct("<I")
lens = Struct("<100sIIIIIIIdd")
shapes = "?EBPM"


def zmf_read(file, session):
    cat = Catalog()
    cat.load(file)
    f = open(file, "rb")
    cat.version, = head.unpack(f.read(head.size))
    assert cat.version in (1001, )
    while True:
        l = Lens()
        li = f.read(lens.size)
        if len(li) != lens.size:
            break
        li = list(lens.unpack(li))
        l.name = li[0].decode("latin1").strip("\0")
        l.shape = shapes[li[3]]
        l.elements = li[2]
        l.aspheric = li[4]
        l.version = li[1]
        l.grin = li[5]
        l.toroidal = li[6]
        l.efl = li[8]
        l.enp = li[9]
        description = f.read(li[7])
        description = zmf_obfuscate(description, l.efl, l.enp)
        description = description.decode("latin1")
        assert description.startswith(f"VERS {l.version:06d}\n")
        l.data = description
        cat.lenses.append(l)
    return cat


def zmf_obfuscate(data, a, b):
    iv = np.cos(6*a + 3*b)
    iv = np.cos(655*(np.pi/180)*iv) + iv
    p = np.arange(len(data))
    k = 13.2*(iv + np.sin(17*(p + 3)))*(p + 1)
    k = (int((f"{_:.8e}")[4:7]) for _ in k)
    data = np.fromstring(data, np.uint8)
    data ^= np.fromiter(k, np.uint8, len(data))
    return data.tostring()
```
