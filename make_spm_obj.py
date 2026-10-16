"""Convert Bruker NanoScope .spm AFM scans into Wavefront .obj meshes.

Written for the TGXYZ02 scans in TGXYZ02.zip, but it reads any NanoScope file
with a height channel. Unlike the idealized artifact meshes here
(tgxyz02_*.obj, vgrp_*.obj), these are MEASURED surfaces: they carry the real
tip-convolution artifact, noise and scanner drift.

Scaling. NanoScope stores height as an integer whose full scale is 2^(8*bytes),
so

    height_nm = (raw / 2**(8*bytes_per_pixel)) * Z_scale_volts * sensitivity_nm_per_V

and the sensitivity to use is the one named in brackets in the Z scale line
("[Sens. ZsensSens]"), not the similarly-named Sens. Zsens. Pixels carrying the
invalid marker (|raw| at the int32 limit) are filled from their neighbours.

Each mesh is a closed solid: the scan as a height field on top, side walls, and
a flat bottom, so it has volume for CAD or printing. Units are NANOMETRES, +Y
up, with the lowest measured point at y = 0.

Usage:
    python3 make_spm_obj.py                        # the two default scans below
    python3 make_spm_obj.py path/to/scan.spm       # any scan
    python3 make_spm_obj.py scan.spm --step 2      # downsample by 2 (fewer triangles)
"""
import os
import re
import sys

import numpy as np

BASE = 1_000.0  # nm of solid below the lowest measured point

DEFAULTS = [
    # (file, downsample) -- one from each pitch area of the TGXYZ02
    ("TGXYZ02/06241807.0_00001.spm", 1),
    ("TGXYZ02/06241746.0_00001.spm", 1),
]


def read_height(path):
    """Return (heights_nm, scan_size_nm, n) for the first height channel."""
    raw = open(path, "rb").read()
    txt = raw[:120_000].decode("latin-1")
    for blk in txt.split("\\*Ciao image list")[1:]:
        def g(k):
            m = re.search(r"\\" + re.escape(k) + r":\s*(.*)", blk)
            return m.group(1).strip() if m else None
        if "Height" not in (g("@2:Image Data") or ""):
            continue
        n = int(g("Samps/line"))
        lines = int(g("Number of lines"))
        bpp = int(g("Bytes/pixel"))
        off = int(g("Data offset"))
        size_txt = g("Scan Size")
        size = float(re.match(r"([\d.]+)", size_txt).group(1))
        if "~m" in size_txt or "um" in size_txt:
            size *= 1000.0                      # header writes µm as "~m"
        zs = g("@2:Z scale")
        zvolt = float(re.search(r"\)\s*([-\d.eE]+)\s*V", zs).group(1))
        key = re.search(r"\[(Sens\.[^\]]+)\]", zs).group(1)
        sens = float(re.search(r"\\@" + re.escape(key) + r":\s*V\s+([\d.]+)", txt).group(1))
        a = np.frombuffer(raw, dtype=np.int32 if bpp == 4 else np.int16,
                          count=n * lines, offset=off).astype(np.float64)
        a[np.abs(a) > 2.1e9] = np.nan           # invalid-pixel marker
        z = (a / 2.0 ** (8 * bpp)) * zvolt * sens
        return z.reshape(lines, n), size, n
    raise ValueError("no height channel in " + path)


def clean(z):
    """Remove sample tilt and fill invalid pixels from their neighbours."""
    n = z.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    ok = np.isfinite(z)
    A = np.c_[xx[ok], yy[ok], np.ones(ok.sum())]
    c, *_ = np.linalg.lstsq(A, z[ok], rcond=None)
    z = z - (c[0] * xx + c[1] * yy + c[2])
    if not ok.all():                            # nearest-neighbour fill, a few passes
        z = z.copy()
        for _ in range(8):
            bad = ~np.isfinite(z)
            if not bad.any():
                break
            filled = z.copy()
            for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                shifted = np.roll(np.roll(z, dy, 0), dx, 1)
                take = bad & np.isfinite(shifted)
                filled[take] = shifted[take]
            z = filled
        z[~np.isfinite(z)] = np.nanmedian(z)
    return z - z.min()


def step_height(z):
    """Separation of the two dominant height levels (nm)."""
    v = z.ravel()
    hist, edges = np.histogram(v, bins=400)
    mids = 0.5 * (edges[1:] + edges[:-1])
    order = np.argsort(hist)[::-1]
    a = mids[order[0]]
    b = next((mids[i] for i in order if abs(mids[i] - a) > 0.25 * (v.max() - v.min())), None)
    return abs(b - a) if b is not None else float("nan")


def mesh(z, size_nm):
    """Closed solid: height field on top, walls, flat bottom."""
    n = z.shape[0]
    px = size_nm / (n - 1)
    verts, faces = [], []
    for i in range(n):                          # top surface
        for j in range(n):
            verts.append((j * px, z[i, j], i * px))
    top = lambda i, j: i * n + j + 1
    for i in range(n - 1):
        for j in range(n - 1):
            a, b, c, d = top(i, j), top(i, j + 1), top(i + 1, j + 1), top(i + 1, j)
            faces += [(a, b, c), (a, c, d)]
    span = (n - 1) * px
    base = -BASE
    corner = {}
    def bv(x, zc):
        key = (round(x, 3), round(zc, 3))
        if key not in corner:
            verts.append((x, base, zc))
            corner[key] = len(verts)
        return corner[key]
    for j in range(n - 1):                      # walls along z = 0 and z = span
        faces += [(top(0, j), bv(j * px, 0), bv((j + 1) * px, 0)), (top(0, j), bv((j + 1) * px, 0), top(0, j + 1))]
        faces += [(top(n - 1, j + 1), bv((j + 1) * px, span), bv(j * px, span)), (top(n - 1, j + 1), bv(j * px, span), top(n - 1, j))]
    for i in range(n - 1):                      # walls along x = 0 and x = span
        faces += [(top(i + 1, 0), bv(0, (i + 1) * px), bv(0, i * px)), (top(i + 1, 0), bv(0, i * px), top(i, 0))]
        faces += [(top(i, n - 1), bv(span, i * px), bv(span, (i + 1) * px)), (top(i, n - 1), bv(span, (i + 1) * px), top(i + 1, n - 1))]
    # bottom, split on the same grid so its edges match the walls
    for j in range(n - 1):
        for i in range(n - 1):
            a, b = bv(j * px, i * px), bv((j + 1) * px, i * px)
            c, d = bv((j + 1) * px, (i + 1) * px), bv(j * px, (i + 1) * px)
            faces += [(a, c, b), (a, d, c)]
    return verts, faces


def convert(path, down=1, out=None):
    z, size, n = read_height(path)
    z = clean(z)
    if down > 1:
        z = z[::down, ::down]
        n = z.shape[0]
    verts, faces = mesh(z, size)
    faces = [(a, c, b) for (a, b, c) in faces]      # wind CCW seen from outside
    tag = "" if down == 1 else "_%dpx" % n     # keep downsampled output separate
    out = out or ("scan_%s_%.0fum%s.obj" % (os.path.basename(path).split(".")[0], size / 1000, tag))
    with open(out, "w") as fh:
        for line in [
            "AFM scan converted from %s" % os.path.basename(path),
            "MEASURED surface: includes the tip-convolution artifact, noise and drift.",
            "Scan %.0f x %.0f um, %d x %d points (%.1f nm/px), height range %.1f nm." % (
                size / 1000, size / 1000, n, n, size / (n - 1), z.max() - z.min()),
            "Two-level step height here: %.1f nm." % step_height(z),
            "Lowest point at y = 0, %.0f um solid below. Units: nanometres, +Y up." % (BASE / 1000),
            "Generated by make_spm_obj.py",
        ]:
            fh.write("# " + line + "\n")
        fh.write("#\n# vertices: %d   triangles: %d\n\n" % (len(verts), len(faces)))
        fh.write("o %s\n" % out[:-4])
        for (x, y, zc) in verts:
            fh.write("v %.3f %.3f %.3f\n" % (x, y, zc))
        for (a, b, c) in faces:
            fh.write("f %d %d %d\n" % (a, b, c))
    print("%-34s %6d verts %7d tris   %.0f µm scan, %dx%d, step %.1f nm"
          % (out, len(verts), len(faces), size / 1000, n, n, step_height(z)))


def main():
    argv = sys.argv[1:]
    down = 1
    if "--step" in argv:
        i = argv.index("--step")
        down = int(argv[i + 1])
        argv = argv[:i] + argv[i + 2:]          # drop the flag and its value
    args = [a for a in argv if not a.startswith("--")]
    if args:
        for a in args:
            convert(a, down)
    else:
        for path, d in DEFAULTS:
            if os.path.exists(path):
                convert(path, d)
            else:
                print("missing %s (unzip TGXYZ02.zip first)" % path)


if __name__ == "__main__":
    main()
