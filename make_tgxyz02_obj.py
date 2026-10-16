"""Generate Wavefront .obj meshes of the MikroMasch TGXYZ02 calibration grating.

TGXYZ02 (https://www.spmtips.com/test-structures-TGXYZ02) is a 3D XYZ standard:
rectangular SiO2 steps on a Si wafer, 100 nm step height (±3%), with a 5 µm
pitch area (circular features and X/Y lines) and a 10 µm pitch area (square
features), pitch accurate to 0.1 µm.

Two areas are modelled, and their feature sizes come from the 27 AFM scans in
TGXYZ02.zip, since the vendor publishes pitch and step height but not feature
width (measured from line profiles across the scans, step 104.9 ± 2.6 nm):

    circles   5 µm pitch, circular holes Ø 2.54 µm, 100 nm deep
    squares   10 µm pitch, square pillars 5 µm on a side, 100 nm tall

Both meshes are a 20 x 20 µm patch — the field afm-3d.html shows at its 20 µm
scan — as closed solids on a 1 µm support slab. Units are NANOMETRES, +Y up,
the lower surface at y = 0, matching the other files here.

Run:  python3 make_tgxyz02_obj.py
"""
import math

STEP = 100.0          # nm, datasheet
BASE = 1_000.0        # support slab below the lower surface (nm)
SPAN = 20_000.0       # patch size (nm)

SQ_PITCH = 10_000.0   # nm
SQ_SIDE = 5_000.0     # nm (measured)
CIRC_PITCH = 5_000.0  # nm
CIRC_DIA = 2_540.0    # nm (measured)
N_ARC = 64            # samples around each circular hole (multiple of 8)


class Mesh:
    def __init__(self):
        self.verts, self.index, self.faces = [], {}, []

    def v(self, x, y, z):
        key = (round(x, 3), round(y, 3), round(z, 3))
        if key not in self.index:
            self.index[key] = len(self.verts) + 1
            self.verts.append(key)
        return self.index[key]

    def tri(self, a, b, c):
        if a != b and b != c and a != c:
            self.faces.append((a, b, c))

    def quad(self, a, b, c, d):
        self.tri(a, b, c)
        self.tri(a, c, d)

    def write(self, path, header):
        with open(path, "w") as fh:
            for line in header:
                fh.write("# " + line + "\n")
            fh.write("#\n# vertices: %d   triangles: %d\n\n" % (len(self.verts), len(self.faces)))
            fh.write("o %s\n" % path[:-4])
            for (x, y, z) in self.verts:
                fh.write("v %.4f %.4f %.4f\n" % (x, y, z))
            for (a, b, c) in self.faces:
                fh.write("f %d %d %d\n" % (a, b, c))


def skirt_and_base(m, cells, cell, top_of):
    """Outer walls (in bands, so edges meet) and the slab bottom, per cell."""
    span = cells * cell
    for i in range(cells):
        for j in range(cells):
            x0, x1 = i * cell, (i + 1) * cell
            z0, z1 = j * cell, (j + 1) * cell
            h = top_of(i, j)
            bands = [(-BASE, 0.0)] + ([(0.0, h)] if h > 0 else [])
            for (lo, hi) in bands:
                if i == 0:
                    m.quad(m.v(0, lo, z0), m.v(0, lo, z1), m.v(0, hi, z1), m.v(0, hi, z0))
                if i == cells - 1:
                    m.quad(m.v(span, lo, z1), m.v(span, lo, z0), m.v(span, hi, z0), m.v(span, hi, z1))
                if j == 0:
                    m.quad(m.v(x1, lo, 0), m.v(x0, lo, 0), m.v(x0, hi, 0), m.v(x1, hi, 0))
                if j == cells - 1:
                    m.quad(m.v(x0, lo, span), m.v(x1, lo, span), m.v(x1, hi, span), m.v(x0, hi, span))
            m.quad(m.v(x0, -BASE, z0), m.v(x1, -BASE, z0), m.v(x1, -BASE, z1), m.v(x0, -BASE, z1))


def squares():
    """10 µm pitch, square pillars: each pitch cell splits into 2x2 half-cells."""
    cell = SQ_PITCH / 2.0                      # 5 µm squares, pillar when both even
    cells = int(round(SPAN / cell))
    top = lambda i, j: STEP if (i % 2 == 0 and j % 2 == 0) else 0.0
    m = Mesh()
    for i in range(cells):
        for j in range(cells):
            x0, x1 = i * cell, (i + 1) * cell
            z0, z1 = j * cell, (j + 1) * cell
            h = top(i, j)
            m.quad(m.v(x0, h, z0), m.v(x0, h, z1), m.v(x1, h, z1), m.v(x1, h, z0))
            for (di, dj) in ((1, 0), (0, 1)):
                ni, nj = i + di, j + dj
                if ni >= cells or nj >= cells:
                    continue
                nh = top(ni, nj)
                if nh == h:
                    continue
                lo, hi = min(h, nh), max(h, nh)
                if di:
                    a, b = m.v(x1, lo, z0), m.v(x1, lo, z1)
                    c, d = m.v(x1, hi, z1), m.v(x1, hi, z0)
                    m.quad(a, b, c, d) if h < nh else m.quad(a, d, c, b)
                else:
                    a, b = m.v(x0, lo, z1), m.v(x1, lo, z1)
                    c, d = m.v(x1, hi, z1), m.v(x0, hi, z1)
                    m.quad(a, d, c, b) if h < nh else m.quad(a, b, c, d)
    skirt_and_base(m, cells, cell, top)
    return m, cells, cell


def circles():
    """5 µm pitch, circular holes: top plane fanned out from each hole rim.

    Every cell samples the same angles, so the boundary points two neighbouring
    cells put on their shared edge coincide and the mesh stays closed.
    """
    cell = CIRC_PITCH
    cells = int(round(SPAN / cell))
    r = CIRC_DIA / 2.0
    half = cell / 2.0
    m = Mesh()
    angles = [2 * math.pi * k / N_ARC for k in range(N_ARC)]

    def on_square(a):
        """Where the ray at angle a leaves the cell (local coords)."""
        ca, sa = math.cos(a), math.sin(a)
        t = min(half / abs(ca) if ca else 1e18, half / abs(sa) if sa else 1e18)
        return t * ca, t * sa

    span = cells * cell
    for i in range(cells):
        for j in range(cells):
            cx, cz = (i + 0.5) * cell, (j + 0.5) * cell
            for k in range(N_ARC):
                a0, a1 = angles[k], angles[(k + 1) % N_ARC]
                rim0 = m.v(cx + r * math.cos(a0), STEP, cz + r * math.sin(a0))
                rim1 = m.v(cx + r * math.cos(a1), STEP, cz + r * math.sin(a1))
                sx0, sz0 = on_square(a0)
                sx1, sz1 = on_square(a1)
                ox0, oz0 = cx + sx0, cz + sz0
                ox1, oz1 = cx + sx1, cz + sz1
                out0, out1 = m.v(ox0, STEP, oz0), m.v(ox1, STEP, oz1)
                m.quad(rim0, out0, out1, rim1)              # top plane, outward
                flo0 = m.v(cx + r * math.cos(a0), 0.0, cz + r * math.sin(a0))
                flo1 = m.v(cx + r * math.cos(a1), 0.0, cz + r * math.sin(a1))
                m.quad(rim0, rim1, flo1, flo0)              # hole wall
                m.tri(m.v(cx, 0.0, cz), flo0, flo1)         # hole floor
                # slab bottom: the same fan, so every edge has a partner
                bot0, bot1 = m.v(ox0, -BASE, oz0), m.v(ox1, -BASE, oz1)
                m.tri(m.v(cx, -BASE, cz), bot1, bot0)
                # outer wall wherever this segment lies on the patch boundary
                on_edge = (
                    (abs(ox0) < 1e-6 and abs(ox1) < 1e-6)
                    or (abs(ox0 - span) < 1e-6 and abs(ox1 - span) < 1e-6)
                    or (abs(oz0) < 1e-6 and abs(oz1) < 1e-6)
                    or (abs(oz1 - span) < 1e-6 and abs(oz0 - span) < 1e-6)
                )
                if on_edge:
                    m.quad(out0, bot0, bot1, out1)
    # the fan is wound the other way round than the blocky builder
    m.faces = [(a, c, b) for (a, b, c) in m.faces]
    return m, cells, cell


def main():
    common = [
        "MikroMasch TGXYZ02 XYZ calibration grating (SiO2 steps on Si).",
        "Step height %.0f nm (datasheet, +/-3%%); pitch accurate to 0.1 um." % STEP,
        "Feature widths measured from the 27 AFM scans in TGXYZ02.zip",
        "(step there measured 104.9 +/- 2.6 nm, pitches 5.00 and 10.00 um).",
        "Lower surface at y = 0, %.0f um support slab below. Units: nanometres, +Y up." % (BASE / 1000),
        "Generated by make_tgxyz02_obj.py",
    ]
    m, cells, cell = squares()
    m.write("tgxyz02_squares_20um.obj",
            [common[0], "10 um pitch area: square pillars %.1f um on a side, %.0f nm tall, over %.0f x %.0f um."
             % (SQ_SIDE / 1000, STEP, SPAN / 1000, SPAN / 1000)] + common[1:])
    print("tgxyz02_squares_20um.obj  %5d verts %5d tris   %d x %d cells of %.1f um"
          % (len(m.verts), len(m.faces), cells, cells, cell / 1000))
    m, cells, cell = circles()
    m.write("tgxyz02_circles_20um.obj",
            [common[0], "5 um pitch area: circular holes %.2f um across, %.0f nm deep, over %.0f x %.0f um."
             % (CIRC_DIA / 1000, STEP, SPAN / 1000, SPAN / 1000)] + common[1:])
    print("tgxyz02_circles_20um.obj  %5d verts %5d tris   %d x %d cells of %.1f um"
          % (len(m.verts), len(m.faces), cells, cells, cell / 1000))


if __name__ == "__main__":
    main()
