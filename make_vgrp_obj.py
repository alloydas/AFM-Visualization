"""Generate Wavefront .obj meshes of the Bruker VGRP SPM calibration grating.

VGRP-15M (https://www.brukerafmprobes.com/product/3538/vgrp-15m) and VGRP-GS
(https://www.brukerafmprobes.com/product/3536/vgrp-gs) are the same artifact --
"Calibration Artifact: 180nm Depth, 10um Pitch, Pt Coated" -- a 2D array of
square pits. They differ only in the mount (15 mm metal disc vs glass slide),
which is not part of the scanned geometry, so one mesh serves both.

    pitch      10 µm in x and y
    pit        5 µm square (half the pitch), 180 nm deep, vertical walls
    coating    Pt (not modelled; it is conformal and thin)

Two meshes are written:

  vgrp_grating_20um.obj   a 20 x 20 µm patch (2 x 2 pitch cells, 4 pits) --
                          the same field afm-3d.html shows at its 20 µm scan.
  vgrp_cell_10um.obj      one 10 µm unit cell (1 pit), tileable.

Both are closed solids: pit floor at y = 0, plateau at y = 180 nm, and a 1 µm
support slab below, so they have volume for CAD, FEA or ray tracing rather than
being an open surface. Units are NANOMETRES, +Y up, matching the simulator.

Run:  python3 make_vgrp_obj.py
"""

PITCH = 10_000.0      # nm
DEPTH = 180.0         # nm
DUTY = 0.5            # pit side as a fraction of the pitch
BASE = 1_000.0        # support slab thickness below the pit floor (nm)
CELL = PITCH * DUTY   # 5 µm square grid: each square is all pit or all plateau


def build(cells):
    """Blocky height field on a `cells` x `cells` grid of CELL-sized squares.

    A square is a pit when both its indices are even, which puts one 5 µm pit
    in each 10 µm pitch cell. Returns merged vertices and outward triangles.
    """
    top = [[0.0 if (i % 2 == 0 and j % 2 == 0) else DEPTH for j in range(cells)]
           for i in range(cells)]
    verts, index, faces = [], {}, []

    def v(x, y, z):
        key = (round(x, 4), round(y, 4), round(z, 4))
        if key not in index:
            index[key] = len(verts) + 1          # .obj is 1-based
            verts.append(key)
        return index[key]

    def quad(a, b, c, d):
        faces.append((a, b, c))
        faces.append((a, c, d))

    for i in range(cells):
        for j in range(cells):
            x0, x1 = i * CELL, (i + 1) * CELL
            z0, z1 = j * CELL, (j + 1) * CELL
            h = top[i][j]
            # top face (up)
            quad(v(x0, h, z0), v(x0, h, z1), v(x1, h, z1), v(x1, h, z0))
            # walls to the neighbour in +x and +z where the height steps
            for (di, dj) in ((1, 0), (0, 1)):
                ni, nj = i + di, j + dj
                if ni >= cells or nj >= cells:
                    continue
                nh = top[ni][nj]
                if nh == h:
                    continue
                lo, hi = min(h, nh), max(h, nh)
                if di:                                    # wall runs along z at x1
                    a, b = v(x1, lo, z0), v(x1, lo, z1)
                    c, d = v(x1, hi, z1), v(x1, hi, z0)
                    quad(a, b, c, d) if h < nh else quad(a, d, c, b)
                else:                                     # wall runs along x at z1
                    a, b = v(x0, lo, z1), v(x1, lo, z1)
                    c, d = v(x1, hi, z1), v(x0, hi, z1)
                    quad(a, d, c, b) if h < nh else quad(a, b, c, d)

    span = cells * CELL
    # Outer skirt, in two bands so its edges line up with the pit walls: every
    # boundary square gets -BASE..0, and plateau squares get 0..DEPTH on top.
    # (One tall quad per square instead would leave T-junctions at pit corners.)
    for i in range(cells):
        for j in range(cells):
            x0, x1 = i * CELL, (i + 1) * CELL
            z0, z1 = j * CELL, (j + 1) * CELL
            bands = [(-BASE, 0.0)] + ([(0.0, DEPTH)] if top[i][j] > 0 else [])
            for (lo, hi) in bands:
                if i == 0:
                    quad(v(0, lo, z0), v(0, lo, z1), v(0, hi, z1), v(0, hi, z0))
                if i == cells - 1:
                    quad(v(span, lo, z1), v(span, lo, z0), v(span, hi, z0), v(span, hi, z1))
                if j == 0:
                    quad(v(x1, lo, 0), v(x0, lo, 0), v(x0, hi, 0), v(x1, hi, 0))
                if j == cells - 1:
                    quad(v(x0, lo, span), v(x1, lo, span), v(x1, hi, span), v(x0, hi, span))
    # slab bottom (down), split per square so its edges match the skirt's
    for i in range(cells):
        for j in range(cells):
            x0, x1 = i * CELL, (i + 1) * CELL
            z0, z1 = j * CELL, (j + 1) * CELL
            quad(v(x0, -BASE, z0), v(x1, -BASE, z0), v(x1, -BASE, z1), v(x0, -BASE, z1))
    return verts, faces


def write_obj(path, header, verts, faces):
    with open(path, "w") as fh:
        for line in header:
            fh.write("# " + line + "\n")
        fh.write("#\n# vertices: %d   triangles: %d\n\n" % (len(verts), len(faces)))
        fh.write("o %s\n" % path[:-4])
        for (x, y, z) in verts:
            fh.write("v %.4f %.4f %.4f\n" % (x, y, z))
        for (a, b, c) in faces:
            fh.write("f %d %d %d\n" % (a, b, c))


def main():
    common = [
        "Bruker VGRP SPM calibration grating (VGRP-15M and VGRP-GS are the same",
        "artifact; they differ only in the mount, 15 mm disc vs glass slide).",
        "2D array of square pits: %.0f µm pitch, %.0f µm pits, %.0f nm deep, vertical walls."
        % (PITCH / 1000, PITCH * DUTY / 1000, DEPTH),
        "Pit floor at y = 0, plateau at y = %.0f nm, %.0f µm support slab below."
        % (DEPTH, BASE / 1000),
        "Pt coating not modelled (thin and conformal). Units: nanometres, +Y up.",
        "Generated by make_vgrp_obj.py",
    ]
    for name, cells, note in [
        ("vgrp_grating_20um.obj", 4,
         "20 x 20 µm patch (2 x 2 pitch cells, 4 pits) — the field afm-3d.html shows at 20 µm."),
        ("vgrp_cell_10um.obj", 2,
         "One 10 µm unit cell (1 pit), tileable."),
    ]:
        verts, faces = build(cells)
        write_obj(name, [common[0], common[1], note] + common[2:], verts, faces)
        span = cells * CELL / 1000
        print("%-24s %5d verts %5d tris   %.0f x %.0f µm" % (name, len(verts), len(faces), span, span))


if __name__ == "__main__":
    main()
