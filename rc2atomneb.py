#!/usr/bin/env python3
"""
rc2atomneb.py

Convert recombination-line atomic data from the AtomNeb rc_data source tree
into AtomNeb-style FITS files under an atomic-data-rc directory.

This script is a Python analogue of the IDL generator gen_rc_atomneb.pro.\n\nVersion notes:\n  v2 fixes Ne III BR output, blank List string cells, and SH95 trailing-row handling.\n  v3 excludes SSB17 from default builds and requires the full OIIlines_ABC data when SSB17 is requested.

Input tree expected, for example:
    rc_data/
      rc_collection/
      rc_PPB91/
      rc_SH95/
      rc_he_ii_PFSD12/
      rc_n_iii_FSL13/
      rc_o_iii_SSB17/

Outputs:
    rc_collection.fits
    rc_PPB91.fits
    rc_SH95.fits
    rc_he_ii_PFSD12.fits
    rc_n_iii_FSL13.fits
    rc_o_iii_SSB17.fits              # optional: only when --collections includes ssb17
    rc_o_iii_SSB17_orl_case_b.fits   # optional: only when --collections includes ssb17

Requires:
    numpy
    astropy

Notes
-----
* SSB17 O II files are not built by default. To build them, request
  --collections ssb17 or include ssb17 in the collection list. The full
  OIIlines_ABC.txt or OIIlines_ABC file must be present in rc_o_iii_SSB17.
* This script intentionally writes List and References HDUs followed by data
  extensions, matching the generator layout used by gen_rc_atomneb.pro.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

def require_astropy():
    try:
        from astropy.io import fits
    except Exception as exc:
        raise SystemExit(
            "ERROR: astropy is required to write FITS files.\n"
            "Install it with: pip install astropy numpy"
        ) from exc
    return fits


def ascii_clean(s: Any, maxlen: int | None = None) -> str:
    out = "" if s is None else str(s)
    out = out.replace("\u2013", "-").replace("\u2014", "-")
    out = out.replace("\u2018", "'").replace("\u2019", "'")
    out = out.replace("\u201c", '"').replace("\u201d", '"')
    out = out.encode("ascii", "replace").decode("ascii").strip()
    if maxlen is not None:
        out = out[:maxlen]
    return out


def parse_float(x: Any, default: float = np.nan) -> float:
    try:
        s = str(x).strip().replace("D", "E").replace("d", "e")
        if s == "" or s in {"-", "..."}:
            return default
        return float(s)
    except Exception:
        return default


def parse_int(x: Any, default: int = 0) -> int:
    try:
        s = str(x).strip()
        if s == "" or s in {"-", "..."}:
            return default
        return int(float(s))
    except Exception:
        return default


def noncomment_lines(path: Path) -> list[str]:
    return [
        ln.rstrip("\n")
        for ln in path.read_text(errors="replace").splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]


def read_list_file(path: Path) -> list[str]:
    if not path.exists():
        return []
    out = []
    for ln in path.read_text(errors="replace").splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        out.append(s)
    return out


def read_references(path: Path) -> list[dict[str, str]]:
    """Read AtomNeb-style reference files: atomic_data & reference."""

    rows: list[dict[str, str]] = []
    if not path.exists():
        return rows
    for ln in path.read_text(errors="replace").splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        if "&" in s:
            key, ref = s.split("&", 1)
        else:
            key, ref = "", s
        rows.append(
            {
                "AtomicData": ascii_clean(key, 128),
                "Reference": ascii_clean(ref, 512),
            }
        )
    return rows


def make_primary_hdu(fits, comments: Sequence[str] | None = None):
    hdu = fits.PrimaryHDU()
    hdr = hdu.header
    hdr["CREATOR"] = ("rc2atomneb.py", "created for AtomNeb-style RC data")
    if comments:
        for line in comments:
            hdr.add_comment(ascii_clean(line, 72))
    return hdu


def infer_string_width(values: Iterable[Any], minimum: int = 1, maximum: int = 512) -> int:
    width = minimum
    for val in values:
        width = max(width, len(ascii_clean(val)))
    return max(minimum, min(maximum, width))


def make_table_hdu(fits, rows: Sequence[dict[str, Any]], extname: str):
    """Create a binary table HDU from a list of dictionaries.

    String columns are ASCII. Integer columns are K. Float columns are D.
    """

    if not rows:
        # Minimal empty table, useful for incomplete short SSB17 examples.
        cols = [fits.Column(name="EMPTY", format="1A", array=np.array([], dtype="S1"))]
        hdu = fits.BinTableHDU.from_columns(cols)
        hdu.header["EXTNAME"] = extname
        return hdu

    names = list(rows[0].keys())
    cols = []
    for name in names:
        vals = [row.get(name) for row in rows]
        # Determine column type.
        if all(isinstance(v, (int, np.integer)) and not isinstance(v, bool) for v in vals):
            arr = np.array(vals, dtype=np.int64)
            cols.append(fits.Column(name=name, format="K", array=arr))
        elif all(isinstance(v, (float, int, np.floating, np.integer)) and not isinstance(v, bool) for v in vals):
            arr = np.array(vals, dtype=np.float64)
            cols.append(fits.Column(name=name, format="D", array=arr))
        else:
            width = infer_string_width(vals, minimum=1, maximum=512)
            # Important for FITS viewers such as fv:
            # a zero-length byte string in a fixed-width character column can
            # appear as NULL.  Store an all-space field for logically blank
            # values so the cell displays as blank instead.
            out_vals = []
            for v in vals:
                s = ascii_clean(v, width)
                if s == "":
                    s = " " * width
                out_vals.append(s.encode("ascii", "replace"))
            arr = np.array(out_vals, dtype=f"S{width}")
            cols.append(fits.Column(name=name, format=f"{width}A", array=arr))
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.header["EXTNAME"] = extname
    return hdu


def make_image_hdu(fits, data: np.ndarray, extname: str):
    arr = np.asarray(data, dtype=np.float64)
    hdu = fits.ImageHDU(data=arr)
    hdu.header["EXTNAME"] = extname
    return hdu


def write_hdul(fits, out_path: Path, hdus: list[Any], overwrite: bool = True):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList(hdus).writeto(out_path, overwrite=overwrite, checksum=True)


def add_extensions_from_list(entries: Sequence[str], offset: int = 3) -> list[dict[str, Any]]:
    rows = []
    for i, fname in enumerate(entries):
        stem = Path(fname).stem
        rows.append({"Aeff_Data": stem, "Extension": i + offset})
    return rows


def read_numeric_table(path: Path, skip_first: int = 0) -> np.ndarray:
    rows = []
    for ln in path.read_text(errors="replace").splitlines()[skip_first:]:
        s = ln.strip()
        if not s or s.startswith("#") or s.startswith("..."):
            continue
        vals = [parse_float(tok) for tok in s.split()]
        vals = [v for v in vals if np.isfinite(v)]
        if vals:
            rows.append(vals)
    if not rows:
        return np.zeros((0, 0), dtype=float)
    n = max(len(r) for r in rows)
    out = np.full((len(rows), n), np.nan, dtype=float)
    for i, row in enumerate(rows):
        out[i, : len(row)] = row
    return out


# -----------------------------------------------------------------------------
# PPB91
# -----------------------------------------------------------------------------

def parse_ppb91_file(path: Path) -> list[dict[str, Any]]:
    rows = []
    for ln in path.read_text(errors="replace").splitlines():
        s = ln.strip()
        if not s:
            continue
        parts = [p.strip() for p in s.split("|")]
        if len(parts) < 10:
            # fallback whitespace, but PPB91 source is pipe separated
            parts = re.split(r"\s+", s)
        if len(parts) < 10:
            continue
        rows.append(
            {
                "Ion": ascii_clean(parts[0], 16),
                "Case1": ascii_clean(parts[1], 8),
                "Wavelength": parse_float(parts[2]),
                "a": parse_float(parts[3]),
                "b": parse_float(parts[4]),
                "c": parse_float(parts[5]),
                "d": parse_float(parts[6]),
                "br": parse_float(parts[7]),
                "Q": ascii_clean(parts[8], 16),
                "Y": parse_float(parts[9]),
            }
        )
    return rows


def build_rc_ppb91(fits, root: Path, out_dir: Path, overwrite: bool = True):
    folder = root / "rc_PPB91"
    if not folder.exists():
        return None
    entries = read_list_file(folder / "list_rc_PPB91.txt")
    list_rows = add_extensions_from_list(entries)
    ref_rows = read_references(folder / "rc_PPB91_references.txt")

    hdus = [
        make_primary_hdu(
            fits,
            [
                "PPB91 recombination coefficients.",
                "Pequignot, Petitjean and Boisson 1991, A&A, 251, 680.",
            ],
        ),
        make_table_hdu(fits, list_rows, "List"),
        make_table_hdu(fits, ref_rows, "References"),
    ]
    for fname in entries:
        stem = Path(fname).stem
        hdus.append(make_table_hdu(fits, parse_ppb91_file(folder / fname), stem))
    out = out_dir / "rc_PPB91.fits"
    write_hdul(fits, out, hdus, overwrite=overwrite)
    return out


# -----------------------------------------------------------------------------
# SH95
# -----------------------------------------------------------------------------

def parse_sh95_file(path: Path) -> np.ndarray:
    """Parse one Storey & Hummer 1995 coefficient grid.

    The first line gives ``temp_num dens_num``.  The original
    ``gen_rc_atomneb.pro`` reads exactly ``temp_num * dens_num`` subsequent
    rows into a 302 x N image array and ignores any trailing row(s).  Some SH95
    source files contain one extra final line after the physical grid; this
    must not be written to the FITS extension.
    """

    lines = noncomment_lines(path)
    if not lines:
        return np.zeros((0, 0), dtype=float)

    header = lines[0].split()
    temp_num = parse_int(header[0], 0) if len(header) >= 1 else 0
    dens_num = parse_int(header[1], 0) if len(header) >= 2 else 0
    ngrid = int(temp_num) * int(dens_num)

    rows = []
    for ln in lines[1:1 + ngrid]:
        vals = [parse_float(tok) for tok in ln.split()]
        vals = [v for v in vals if np.isfinite(v)]
        if vals:
            # IDL uses temp1=dblarr(302) and writes temp1[0:301].
            vals = vals[:302]
            rows.append(vals)

    if not rows:
        return np.zeros((0, 0), dtype=float)

    ncol = 302
    out = np.full((len(rows), ncol), np.nan, dtype=float)
    for i, row in enumerate(rows):
        out[i, : min(len(row), ncol)] = row[:ncol]
    return out


def build_rc_sh95(fits, root: Path, out_dir: Path, overwrite: bool = True):
    folder = root / "rc_SH95"
    if not folder.exists():
        return None
    entries = read_list_file(folder / "list_rc_SH95.txt")
    list_rows = add_extensions_from_list(entries)
    ref_rows = read_references(folder / "rc_SH95_references.txt")

    hdus = [
        make_primary_hdu(
            fits,
            [
                "SH95 hydrogenic recombination coefficients.",
                "Storey and Hummer 1995, MNRAS, 272, 41.",
            ],
        ),
        make_table_hdu(fits, list_rows, "List"),
        make_table_hdu(fits, ref_rows, "References"),
    ]
    for fname in entries:
        stem = Path(fname).stem
        hdus.append(make_image_hdu(fits, parse_sh95_file(folder / fname), stem))
    out = out_dir / "rc_SH95.fits"
    write_hdul(fits, out, hdus, overwrite=overwrite)
    return out


# -----------------------------------------------------------------------------
# PFSD12 He I
# -----------------------------------------------------------------------------

def parse_pfsd_wavelength(path: Path) -> list[dict[str, Any]]:
    rows = []
    for ln in path.read_text(errors="replace").splitlines():
        s = ln.strip()
        if not s:
            continue
        parts = [p.strip() for p in s.split("&")]
        if len(parts) < 3:
            continue
        rows.append(
            {
                "Wavelength": parse_float(parts[0]),
                "LowerTerm": ascii_clean(parts[1], 64),
                "UpperTerm": ascii_clean(parts[2], 64),
            }
        )
    return rows


def build_rc_pfsd12(fits, root: Path, out_dir: Path, overwrite: bool = True):
    folder = root / "rc_he_ii_PFSD12"
    if not folder.exists():
        return None
    entries = read_list_file(folder / "list_rc_he_ii_PFSD12.txt")
    list_rows = add_extensions_from_list(entries)
    ref_rows = read_references(folder / "rc_he_ii_PFSD12_references.txt")

    hdus = [
        make_primary_hdu(
            fits,
            [
                "He I recombination data from Porter et al. 2012/2013.",
            ],
        ),
        make_table_hdu(fits, list_rows, "List"),
        make_table_hdu(fits, ref_rows, "References"),
    ]
    for fname in entries:
        stem = Path(fname).stem
        path = folder / fname
        if stem == "he_ii_wavelength":
            hdus.append(make_table_hdu(fits, parse_pfsd_wavelength(path), stem))
        else:
            hdus.append(make_image_hdu(fits, read_numeric_table(path), stem))
    out = out_dir / "rc_he_ii_PFSD12.fits"
    write_hdul(fits, out, hdus, overwrite=overwrite)
    return out


# -----------------------------------------------------------------------------
# RC collection
# -----------------------------------------------------------------------------

def parse_collection_aeff6(path: Path) -> list[dict[str, Any]]:
    rows = []
    for vals in read_numeric_table(path):
        if len(vals) >= 6:
            rows.append(
                {
                    "Wavelength": vals[0],
                    "a": vals[1],
                    "b": vals[2],
                    "c": vals[3],
                    "d": vals[4],
                    "f": vals[5],
                }
            )
    return rows


def parse_collection_ne_iii_aeff(path: Path) -> list[dict[str, Any]]:
    """Parse Ne III RC-collection coefficients.

    ``ne_iii_aeff.dat`` has seven numeric columns.  The final column is the
    branching ratio ``BR`` and is present in AtomNeb's ``ne_iii_aeff`` table.
    """

    rows = []
    for vals in read_numeric_table(path):
        if len(vals) >= 7:
            rows.append(
                {
                    "Wavelength": vals[0],
                    "a": vals[1],
                    "b": vals[2],
                    "c": vals[3],
                    "d": vals[4],
                    "f": vals[5],
                    "BR": vals[6],
                }
            )
    return rows


def parse_collection_n_iii_aeff(path: Path) -> list[dict[str, Any]]:
    rows = []
    for vals in read_numeric_table(path):
        if len(vals) >= 3:
            rows.append({"a": vals[0], "b": vals[1], "c": vals[2]})
    return rows


def parse_collection_o_iii_aeff(path: Path) -> list[dict[str, Any]]:
    rows = []
    lines = path.read_text(errors="replace").splitlines()
    # First line is a header in the AtomNeb source.
    for ln in lines[1:]:
        s = ln.strip()
        if not s:
            continue
        parts = [p.strip() for p in s.split("|")]
        if len(parts) < 9:
            continue
        rows.append(
            {
                "Term": ascii_clean(parts[0], 64),
                "Case1": ascii_clean(parts[1], 8),
                "a2": parse_float(parts[2]),
                "a4": parse_float(parts[3]),
                "a5": parse_float(parts[4]),
                "a6": parse_float(parts[5]),
                "b": parse_float(parts[6]),
                "c": parse_float(parts[7]),
                "d": parse_float(parts[8]),
            }
        )
    return rows


def parse_br_fixed(path: Path, kind: str) -> list[dict[str, Any]]:
    """Parse N II/O II branching-ratio fixed-width rows.

    The IDL generator keeps only wavelength, g values, multiplet/terms, and BR
    values.  This parser follows that reduced AtomNeb structure.
    """

    rows = []
    for ln in path.read_text(errors="replace").splitlines():
        if not ln.strip():
            continue
        # Fixed-width-ish source; token parsing is adequate for the columns kept.
        toks = ln.split()
        if len(toks) < 12:
            continue
        wavelength = parse_float(toks[1])

        # Find the last one or three BR floats depending on file.
        nums = [parse_float(t) for t in toks]
        if kind == "n":
            br = nums[-1]
            row = {
                "Wavelength": wavelength,
                "BR": br,
                "g1": parse_int(toks[9]) if len(toks) > 9 else 0,
                "g2": parse_int(toks[14]) if len(toks) > 14 else 0,
                "Mult1": ascii_clean(toks[7] if len(toks) > 7 else "", 32),
                "LowerTerm": ascii_clean(toks[11] if len(toks) > 11 else "", 64),
                "UpperTerm": ascii_clean(toks[16] if len(toks) > 16 else "", 64),
            }
        else:
            row = {
                "Wavelength": wavelength,
                "Br_A": nums[-3] if len(nums) >= 3 else np.nan,
                "Br_B": nums[-2] if len(nums) >= 2 else np.nan,
                "Br_C": nums[-1] if len(nums) >= 1 else np.nan,
                "g1": parse_int(toks[9]) if len(toks) > 9 else 0,
                "g2": parse_int(toks[14]) if len(toks) > 14 else 0,
                "Mult1": ascii_clean(toks[7] if len(toks) > 7 else "", 32),
                "LowerTerm": ascii_clean(toks[11] if len(toks) > 11 else "", 64),
                "UpperTerm": ascii_clean(toks[16] if len(toks) > 16 else "", 64),
            }
        rows.append(row)
    return rows


def parse_collection_file(path: Path) -> list[dict[str, Any]]:
    stem = path.stem
    if stem == "c_iii_aeff":
        return parse_collection_aeff6(path)
    if stem == "ne_iii_aeff":
        return parse_collection_ne_iii_aeff(path)
    if stem == "n_iii_aeff":
        return parse_collection_n_iii_aeff(path)
    if stem == "o_iii_aeff":
        return parse_collection_o_iii_aeff(path)
    if stem == "n_iii_br":
        return parse_br_fixed(path, "n")
    if stem == "o_iii_br":
        return parse_br_fixed(path, "o")
    return []


def build_rc_collection(fits, root: Path, out_dir: Path, overwrite: bool = True):
    folder = root / "rc_collection"
    if not folder.exists():
        return None
    entries = read_list_file(folder / "list_rc_collection.txt")
    list_rows = add_extensions_from_list(entries)
    ref_rows = read_references(folder / "rc_collection_references.txt")

    hdus = [
        make_primary_hdu(
            fits,
            [
                "RC Collection recombination coefficients and branching ratios.",
            ],
        ),
        make_table_hdu(fits, list_rows, "List"),
        make_table_hdu(fits, ref_rows, "References"),
    ]
    for fname in entries:
        stem = Path(fname).stem
        hdus.append(make_table_hdu(fits, parse_collection_file(folder / fname), stem))
    out = out_dir / "rc_collection.fits"
    write_hdul(fits, out, hdus, overwrite=overwrite)
    return out


# -----------------------------------------------------------------------------
# FSL13 N II
# -----------------------------------------------------------------------------

def parse_fsl_line(line: str):
    # IDL format: A22, A33, F9, A2, then seven E8 values.
    tr = line[0:22].strip()
    trans = line[22:55].strip()
    wavelength = parse_float(line[55:64])
    tx = line[64:66].strip()
    tail = line[66:].split()
    vals = [parse_float(x) for x in tail[:7]]
    if len(vals) < 7:
        return None
    return tr, trans, wavelength, tx, vals


def build_rc_fsl13(fits, root: Path, out_dir: Path, overwrite: bool = True):
    folder = root / "rc_n_iii_FSL13"
    if not folder.exists():
        return None

    tables = [folder / f"table{i}.dat" for i in (3, 4, 5, 6)]
    if not all(p.exists() for p in tables):
        return None

    lines_by_table = [p.read_text(errors="replace").splitlines() for p in tables]
    n = min(len(x) for x in lines_by_table)
    list_rows = []
    data_arrays = []
    for i in range(n):
        parsed = [parse_fsl_line(lines_by_table[k][i]) for k in range(4)]
        if any(p is None for p in parsed):
            continue
        tr, trans, wave, tx, _vals = parsed[0]
        ind = len(list_rows) + 1
        extname = f"n_iii_aeff_{ind}"
        list_rows.append(
            {
                "Aeff_Data": extname,
                "Extension": ind + 2,
                "IND": ind,
                "Wavelength": wave,
                "Tr": ascii_clean(tr, 64),
                "Trans": ascii_clean(trans, 128),
                "T_X": ascii_clean(tx, 8),
            }
        )
        # Store as 4 density/case rows x 7 temperature coefficients, so aeff[0]
        # returns the seven values shown in AtomNeb examples.
        data_arrays.append(np.array([p[4] for p in parsed], dtype=float))

    refs_txt = ""
    ref_file = folder / "references.txt"
    if ref_file.exists():
        refs_txt = "; ".join([x.strip() for x in ref_file.read_text(errors="replace").splitlines() if x.strip()])
    if not refs_txt:
        refs_txt = "Fang X., Storey P.J., and Liu X.-W., 2011, A&A, 530, A18; 2013, A&A, 550, C2"
    ref_rows = [{"AtomicData": "n_iii_aeff", "Reference": ascii_clean(refs_txt, 512)}]

    hdus = [
        make_primary_hdu(fits, ["FSL13 N II recombination coefficients."]),
        make_table_hdu(fits, list_rows, "List"),
        make_table_hdu(fits, ref_rows, "References"),
    ]
    for row, arr in zip(list_rows, data_arrays):
        hdus.append(make_image_hdu(fits, arr, row["Aeff_Data"]))
    out = out_dir / "rc_n_iii_FSL13.fits"
    write_hdul(fits, out, hdus, overwrite=overwrite)
    return out


# -----------------------------------------------------------------------------
# SSB17 O II
# -----------------------------------------------------------------------------

SSB17_META_SLICES = [
    (0, 7, "IND"),
    (7, 28, "lower_term"),
    (28, 31, "_dash1"),
    (31, 56, "upper_term"),
    (56, 61, "IPI"),
    (61, 65, "NLI"),
    (65, 68, "JI2"),
    (68, 72, "PI"),
    (72, 75, "_dash2"),
    (75, 80, "IPF"),
    (80, 84, "NLF"),
    (84, 87, "JF2"),
    (87, 91, "PF"),
    (91, 102, "Wavelength"),
    (102, 107, "T_X"),
    (107, 119, "E_I"),
    (119, 131, "E_F"),
    (131, 139, "JPI2"),
    (139, 143, "LI"),
    (143, 146, "KI2"),
    (146, 149, "_dash3"),
    (149, 153, "JPF2"),
    (153, 157, "LF"),
    (157, 160, "KF2"),
]


def parse_ssb17_metadata_line(line: str) -> dict[str, Any] | None:
    if not re.match(r"^\s*\d+\s", line):
        return None
    # Skip coefficient rows such as "2.000 3.47E-31 ..."
    if re.match(r"^\s*\d+\.\d+", line):
        return None

    values = {}
    for a, b, name in SSB17_META_SLICES:
        values[name] = line[a:b].strip() if a < len(line) else ""
    ind = parse_int(values.get("IND"), 0)
    if ind <= 0:
        return None
    return {
        "IND": ind,
        "lower_term": ascii_clean(values.get("lower_term"), 64),
        "upper_term": ascii_clean(values.get("upper_term"), 64),
        "IPI": ascii_clean(values.get("IPI"), 16),
        "NLI": ascii_clean(values.get("NLI"), 16),
        "JI2": parse_int(values.get("JI2")),
        "PI": parse_int(values.get("PI")),
        "IPF": ascii_clean(values.get("IPF"), 16),
        "NLF": ascii_clean(values.get("NLF"), 16),
        "JF2": parse_int(values.get("JF2")),
        "PF": parse_int(values.get("PF")),
        "Wavelength": parse_float(values.get("Wavelength")),
        "T_X": ascii_clean(values.get("T_X"), 8),
        "E_I": parse_float(values.get("E_I")),
        "E_F": parse_float(values.get("E_F")),
        "JPI2": parse_int(values.get("JPI2")),
        "LI": parse_int(values.get("LI")),
        "KI2": parse_int(values.get("KI2")),
        "JPF2": parse_int(values.get("JPF2")),
        "LF": parse_int(values.get("LF")),
        "KF2": parse_int(values.get("KF2")),
    }


def parse_ssb17_coefficients(lines: list[str]) -> dict[tuple[int, str], np.ndarray]:
    """Parse coefficient blocks keyed by (line index, case)."""

    blocks: dict[tuple[int, str], np.ndarray] = {}
    i = 0
    header_re = re.compile(r"Line\s+(\d+)\s+Wavelength\s+\[A\]\s+([0-9.Ee+-]+).*Case\s+([ABC])", re.I)
    while i < len(lines):
        m = header_re.search(lines[i])
        if not m:
            i += 1
            continue
        ind = int(m.group(1))
        case = m.group(3).upper()
        i += 1
        # optional density header
        if i < len(lines) and "log10" in lines[i]:
            i += 1
        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                break
            if header_re.search(lines[i]):
                break
            if s.startswith("..."):
                i += 1
                continue
            toks = s.split()
            vals = [parse_float(t) for t in toks]
            vals = [v for v in vals if np.isfinite(v)]
            if len(vals) >= 2:
                # First value is log10(T); IDL stores only coefficient vector.
                rows.append(vals[1:])
            i += 1
        if rows:
            ncol = max(len(r) for r in rows)
            arr = np.full((len(rows), ncol), np.nan, dtype=float)
            for r, row in enumerate(rows):
                arr[r, : len(row)] = row
            blocks[(ind, case)] = arr
    return blocks


def parse_ssb17(path: Path):
    lines = path.read_text(errors="replace").splitlines()
    metadata: dict[int, dict[str, Any]] = {}
    for ln in lines:
        rec = parse_ssb17_metadata_line(ln)
        if rec is not None:
            metadata[rec["IND"]] = rec
    coeffs = parse_ssb17_coefficients(lines)
    return metadata, coeffs


def ssb17_list_row(meta: dict[str, Any], case: str, extension: int, extname: str) -> dict[str, Any]:
    row = {
        "Aeff_Data": extname,
        "Extension": extension,
        "IND": int(meta.get("IND", 0)),
        "Wavelength": float(meta.get("Wavelength", np.nan)),
        "Case1": case,
        "lower_term": meta.get("lower_term", ""),
        "upper_term": meta.get("upper_term", ""),
        "IPI": meta.get("IPI", ""),
        "NLI": meta.get("NLI", ""),
        "JI2": int(meta.get("JI2", 0)),
        "PI": int(meta.get("PI", 0)),
        "IPF": meta.get("IPF", ""),
        "NLF": meta.get("NLF", ""),
        "JF2": int(meta.get("JF2", 0)),
        "PF": int(meta.get("PF", 0)),
        "T_X": meta.get("T_X", ""),
        "E_I": float(meta.get("E_I", np.nan)),
        "E_F": float(meta.get("E_F", np.nan)),
        "JPI2": int(meta.get("JPI2", 0)),
        "LI": int(meta.get("LI", 0)),
        "KI2": int(meta.get("KI2", 0)),
        "JPF2": int(meta.get("JPF2", 0)),
        "LF": int(meta.get("LF", 0)),
        "KF2": int(meta.get("KF2", 0)),
    }
    return row


SSB17_DOWNLOAD_MESSAGE = """
SSB17 O II recombination data were requested, but the required full file was
not found.

Please download OIIlines_ABC from one of:

  https://cdsarc.cds.unistra.fr/ftp/VI/150/DataFiles/OIIlines_ABC
  https://cdsarc.cds.unistra.fr/viz-bin/cat/VI/150

Then place it in:

  rc_data/rc_o_iii_SSB17/OIIlines_ABC

or:

  rc_data/rc_o_iii_SSB17/OIIlines_ABC.txt

The shortened test copy is not enough to reproduce AtomNeb's complete
rc_o_iii_SSB17.fits and rc_o_iii_SSB17_orl_case_b.fits products.
""".strip()


def find_ssb17_source(folder: Path) -> Path:
    """Return the required SSB17 OIIlines_ABC source file.

    AtomNeb's full SSB17 products require the large CDS table.  Accept both the
    original CDS filename without extension and the common .txt name.
    """

    candidates = [folder / "OIIlines_ABC.txt", folder / "OIIlines_ABC"]
    for path in candidates:
        if path.exists() and path.is_file():
            return path
    raise FileNotFoundError(SSB17_DOWNLOAD_MESSAGE)




def build_rc_ssb17_one(
    fits,
    folder: Path,
    out_path: Path,
    mode: str,
    overwrite: bool = True,
    wavelength_min: float = 3500.0,
    wavelength_max: float = 9000.0,
):
    source = find_ssb17_source(folder)

    metadata, coeffs = parse_ssb17(source)
    rows = []
    arrays = []

    # If metadata are incomplete, still allow blocks to be written using a
    # minimal metadata record.
    for (ind, case), arr in sorted(coeffs.items()):
        if mode == "case_b_optical":
            if case != "B":
                continue
        meta = metadata.get(ind)
        if meta is None:
            meta = {"IND": ind, "Wavelength": np.nan}
        if mode == "case_b_optical":
            wave = float(meta.get("Wavelength", np.nan))
            if not (wavelength_min <= wave <= wavelength_max):
                continue
        extname = f"o_iii_aeff_{case.lower()}_{ind}"
        rows.append(ssb17_list_row(meta, case, len(rows) + 3, extname))
        arrays.append(arr)

    # For full mode, include metadata cases even when coefficient blocks are
    # absent only if the full source has no parsed blocks.  This avoids creating
    # thousands of empty image extensions accidentally.
    ref_rows = [
        {
            "AtomicData": "o_iii_aeff",
            "Reference": "Storey, P.J., Sochi, T. and Bastin, R. 2017, MNRAS, 470, 379; VizieR On-line Data Catalog: VI/150",
        }
    ]
    ref_file = folder / "references.txt"
    if ref_file.exists():
        refs = [x.strip() for x in ref_file.read_text(errors="replace").splitlines() if x.strip()]
        if refs:
            ref_rows[0]["Reference"] = "; ".join(refs)

    hdus = [
        make_primary_hdu(fits, ["SSB17 O II recombination coefficients."]),
        make_table_hdu(fits, rows, "List"),
        make_table_hdu(fits, ref_rows, "References"),
    ]
    for row, arr in zip(rows, arrays):
        hdus.append(make_image_hdu(fits, arr, row["Aeff_Data"]))
    write_hdul(fits, out_path, hdus, overwrite=overwrite)
    return out_path


def build_rc_ssb17(
    fits,
    root: Path,
    out_dir: Path,
    overwrite: bool = True,
    make_full: bool = True,
    make_case_b_optical: bool = True,
):
    folder = root / "rc_o_iii_SSB17"
    if not folder.exists():
        return []
    outs = []
    if make_full:
        out = build_rc_ssb17_one(fits, folder, out_dir / "rc_o_iii_SSB17.fits", "full", overwrite=overwrite)
        if out:
            outs.append(out)
    if make_case_b_optical:
        out = build_rc_ssb17_one(
            fits,
            folder,
            out_dir / "rc_o_iii_SSB17_orl_case_b.fits",
            "case_b_optical",
            overwrite=overwrite,
        )
        if out:
            outs.append(out)
    return outs


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def build_all(args) -> list[Path]:
    fits = require_astropy()
    root = Path(args.rc_root).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs: list[Path] = []

    builders = {
        "collection": lambda: build_rc_collection(fits, root, out_dir, overwrite=not args.no_overwrite),
        "ppb91": lambda: build_rc_ppb91(fits, root, out_dir, overwrite=not args.no_overwrite),
        "sh95": lambda: build_rc_sh95(fits, root, out_dir, overwrite=not args.no_overwrite),
        "pfsd12": lambda: build_rc_pfsd12(fits, root, out_dir, overwrite=not args.no_overwrite),
        "fsl13": lambda: build_rc_fsl13(fits, root, out_dir, overwrite=not args.no_overwrite),
    }

    requested = [x.strip().lower() for x in args.collections.split(",") if x.strip()]
    if "all" in requested:
        requested = ["collection", "ppb91", "sh95", "pfsd12", "fsl13"]

    for name in requested:
        if name == "ssb17":
            outputs.extend(
                build_rc_ssb17(
                    fits,
                    root,
                    out_dir,
                    overwrite=not args.no_overwrite,
                    make_full=not args.ssb17_no_full,
                    make_case_b_optical=not args.ssb17_no_case_b_optical,
                )
            )
            continue
        builder = builders.get(name)
        if builder is None:
            print(f"WARNING: unknown collection {name!r}; skipping", file=sys.stderr)
            continue
        out = builder()
        if out is not None:
            outputs.append(out)

    return outputs


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Convert AtomNeb recombination source data to AtomNeb-style FITS files."
    )
    ap.add_argument(
        "--rc-root",
        default="./rc_data",
        help="Root directory containing rc_collection, rc_PPB91, rc_SH95, etc.",
    )
    ap.add_argument(
        "--out-dir",
        default="./atomic-data-rc",
        help="Output directory for FITS files.",
    )
    ap.add_argument(
        "--collections",
        default="all",
        help=(
            "Comma-separated collections to build. Default/all builds collection, "
            "ppb91, sh95, pfsd12, and fsl13. Add ssb17 explicitly to build the "
            "large optional O II SSB17 products."
        ),
    )
    ap.add_argument(
        "--ssb17-no-full",
        action="store_true",
        help="Do not write rc_o_iii_SSB17.fits.",
    )
    ap.add_argument(
        "--ssb17-no-case-b-optical",
        action="store_true",
        help="Do not write rc_o_iii_SSB17_orl_case_b.fits.",
    )
    ap.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Do not overwrite existing FITS files.",
    )
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        outputs = build_all(args)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    if not outputs:
        print("No FITS files were written. Check --rc-root and --collections.", file=sys.stderr)
        return 1
    print("Wrote:")
    for out in outputs:
        print(f"  {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
