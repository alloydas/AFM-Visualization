"""Minimal Bruker Nanoscope .spm reader (Height / AFM image channels)."""

from __future__ import annotations

import re
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HEADER_MARK = b"\\*Ciao image list"
DATA_OFFSET_RE = re.compile(r"\\Data offset:\s*(\d+)")
DATA_LEN_RE = re.compile(r"\\Data length:\s*(\d+)")
BYTES_PER_PIXEL_RE = re.compile(r"\\Bytes/pixel:\s*(\d+)")
SAMPS_RE = re.compile(r"\\Samps/line:\s*(\d+)")
LINES_RE = re.compile(r"\\Number of lines:\s*(\d+)")
SCAN_SIZE_RE = re.compile(
    r"\\Scan Size:\s*([\d.]+)\s*([\d.]+)?\s*(?:~m|nm|um|µm)?", re.I
)
CHANNEL_RE = re.compile(r'\\@2:Image Data:\s*S\s*\[[^\]]+\]\s*"([^"]+)"')
ZSCALE_RE = re.compile(
    r"\\@2:Z scale:\s*V\s*\[[^\]]+\]\s*\(([\d.eE+-]+)\s*[^)]+\)\s*([\d.eE+-]+)\s*V"
)
ZMAGNIFY_RE = re.compile(r"\\@Z magnify:\s*C\s*\[[^\]]+\]\s*([\d.eE+-]+)")
LINE_DIR_RE = re.compile(r"\\Line Direction:\s*(\w+)")
PLANE_FIT_RE = re.compile(r"\\Plane fit:\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+(\d+)")


@dataclass
class SpmImage:
    name: str
    data: np.ndarray
    scan_size_x_nm: float
    scan_size_y_nm: float
    line_direction: str
    z_sensitivity_nm_per_v: float
    plane_fit: tuple[float, float, float, int] | None


@dataclass
class SpmFile:
    path: Path
    images: list[SpmImage]
    header_scan_size_x_nm: float | None
    header_scan_size_y_nm: float | None


def _read_header_text(raw: bytes) -> str:
    end = raw.find(b"\x0c\x04")
    if end == -1:
        end = min(len(raw), 2_000_000)
    return raw[:end].decode("latin-1", errors="replace")


def _parse_scan_size_nm(text: str) -> tuple[float | None, float | None]:
    m = SCAN_SIZE_RE.search(text)
    if not m:
        return None, None
    x = float(m.group(1))
    y = float(m.group(2)) if m.group(2) else x
    if x < 500:
        x *= 1000.0
        y *= 1000.0
    return x, y


def _find_zsens_nm_per_v(header: str) -> float:
    m = re.search(r"\\@Sens\.\s+ZsensSens:\s*V\s*([\d.eE+-]+)\s*nm/V", header)
    if m:
        return float(m.group(1))
    m = re.search(r"\\@Sens\.\s+ZsensSens:\s*V\s*([\d.eE+-]+)", header)
    return float(m.group(1)) if m else 1592.393


def _decode_heights(
    raw: bytes,
    offset: int,
    length: int,
    bpp: int,
    width: int,
    height: int,
    zscale_v: float,
    zsens_nm_per_v: float,
    magnify: float,
) -> np.ndarray:
    if bpp not in (2, 4):
        raise ValueError(f"Unsupported bytes/pixel: {bpp}")
    count = length // bpp
    expected = width * height
    if count != expected:
        if count % width == 0:
            height = count // width
            expected = count
        else:
            raise ValueError(
                f"Pixel count mismatch: {count} pixels from length/bpp, width={width}"
            )
    dtype = ">i2" if bpp == 2 else ">i4"
    arr = np.frombuffer(raw, dtype=dtype, count=count, offset=offset).astype(np.float64)
    if arr.size != expected:
        raise ValueError(f"Pixel count mismatch: got {arr.size}, expected {expected}")
    # Gwyddion/Nanoscope: height = raw * 2 * ZScale * ZSens * magnify / 2^(8*bpp)
    denom = float(1 << (8 * bpp))
    scale = 2.0 * zscale_v * zsens_nm_per_v * magnify / denom
    h = arr.reshape(height, width) * scale
    return h


def _split_image_blocks(header: str) -> list[str]:
    parts = header.split("\\*Ciao image list")
    return [p for p in parts[1:] if p.strip()]


def load_spm(path: str | Path) -> SpmFile:
    path = Path(path)
    raw = path.read_bytes()
    header = _read_header_text(raw)
    hx, hy = _parse_scan_size_nm(header)
    zsens = _find_zsens_nm_per_v(header)
    images: list[SpmImage] = []

    for block in _split_image_blocks(header):
        ch = CHANNEL_RE.search(block)
        if not ch:
            continue
        name = ch.group(1)
        off_m = DATA_OFFSET_RE.search(block)
        len_m = DATA_LEN_RE.search(block)
        bpp_m = BYTES_PER_PIXEL_RE.search(block)
        sw_m = SAMPS_RE.search(block)
        ln_m = LINES_RE.search(block)
        zs_m = ZSCALE_RE.search(block)
        mag_m = ZMAGNIFY_RE.search(block)
        if not all([off_m, len_m, bpp_m, sw_m, zs_m]):
            continue
        if not ln_m:
            continue
        offset = int(off_m.group(1))
        length = int(len_m.group(1))
        bpp = int(bpp_m.group(1))
        width = int(sw_m.group(1))
        height = int(ln_m.group(1))
        lsb_v = float(zs_m.group(1))
        zscale_v = float(zs_m.group(2))
        magnify = float(mag_m.group(1)) if mag_m else 1.0
        _ = lsb_v
        sx, sy = _parse_scan_size_nm(block)
        if sx is None:
            sx, sy = hx or 0.0, hy or sx or 0.0
        line_dir = LINE_DIR_RE.search(block)
        pf = PLANE_FIT_RE.search(block)
        plane = (
            tuple(float(pf.group(i)) for i in range(1, 5))  # type: ignore[misc]
            if pf
            else None
        )
        data = _decode_heights(
            raw, offset, length, bpp, width, height, zscale_v, zsens, magnify
        )
        images.append(
            SpmImage(
                name=name,
                data=data,
                scan_size_x_nm=sx,
                scan_size_y_nm=sy,
                line_direction=line_dir.group(1) if line_dir else "",
                z_sensitivity_nm_per_v=zsens,
                plane_fit=plane,
            )
        )

    return SpmFile(path=path, images=images, header_scan_size_x_nm=hx, header_scan_size_y_nm=hy)


def get_channel(spm: SpmFile, substring: str) -> SpmImage:
    key = substring.lower()
    for img in spm.images:
        if key in img.name.lower():
            return img
    names = [i.name for i in spm.images]
    raise KeyError(f"No channel matching {substring!r}; available: {names}")


def subtract_plane(z: np.ndarray) -> np.ndarray:
    gh, gw = z.shape
    yy, xx = np.mgrid[0:gh, 0:gw]
    a = np.column_stack([xx.ravel(), yy.ravel(), np.ones(gw * gh)])
    coef, _, _, _ = np.linalg.lstsq(a, z.ravel(), rcond=None)
    plane = (coef[0] * xx + coef[1] * yy + coef[2]).astype(np.float64)
    return z - plane
