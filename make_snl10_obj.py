"""Generate Wavefront .obj meshes of the Bruker SNL-10 probe tip.

Datasheet (https://www.brukerafmprobes.com/p-3693-snl-10.aspx): a sharpened
silicon pyramid on a silicon-nitride cantilever,

    front angle 15 +/- 2.5 deg      (+Z, the sharp face)
    back  angle 25 +/- 2.5 deg      (-Z, the blunt face)
    side  angle 22.5 +/- 2.5 deg    (+/-X)
    tip radius 2 nm nominal, 12 nm max
    tip height 2.5 - 8.0 um

Angles are half-angles from vertical, so the sharp pyramid has

    T_pyr(x,z) = max( |x|*cot(side),  z*cot(front) if z >= 0 else |z|*cot(back) )

where T is the height of the tip surface above the apex.

APEX ROUNDING.  A radius-R apex is not "max(pyramid, sphere)" -- for a sphere
sitting at the sharp apex that term never wins, and the apex stays sharp.  The
correct shape is the morphological OPENING of the pyramid by a ball of radius R:
every ball of radius R that fits inside the pyramid, unioned together.  That
restores the four faces exactly, rounds the apex to radius R, and also fillets
the four edges with radius R, which is what a real sharpened tip looks like.
It is computed here as

    T(x,z) = min over |d| <= R of [ T'(x+dx, z+dz) - sqrt(R^2 - |d|^2) ]

with T' the pyramid whose faces are pushed inward by R (the set of allowed ball
centres), i.e. each face raised by R/sin(its angle).  Because the four faces of
an asymmetric pyramid have no common inscribed sphere, the apex is a blended
patch rather than one spherical cap -- as it is on a real probe.

Two meshes are written:

  snl10_probe.obj   whole tip, height 5000 nm (middle of the 2.5-8 um range).
  snl10_apex.obj    the first 100 nm above the apex, where the 2 nm rounding
                    and the edge fillets are actually resolvable.

Convention (matching afm-3d.html): units are NANOMETRES, +Y up, apex at the
origin, +Z the sharp front face.  Meshes are closed and manifold.

Run:  python3 make_snl10_obj.py
"""
import math
import numpy as np

# ---- datasheet geometry ----
FRONT, BACK, SIDE = 15.0, 25.0, 22.5      # half-angles from vertical (deg)
R_APEX = 2.0                              # nominal apex radius (nm); max 12
PROBE_H = 5000.0                          # tip height to model (nm)
APEX_H = 100.0                            # near-apex cutaway height (nm)

N_THETA = 256                             # directions around the axis
N_RING = 140                              # rings from apex to rim
N_DISC = 81                               # ball sampling (odd, so offset 0 is included)

cotF, cotB, cotS = (1.0 / math.tan(math.radians(a)) for a in (FRONT, BACK, SIDE))
tanF, tanB, tanS = (math.tan(math.radians(a)) for a in (FRONT, BACK, SIDE))
sinF, sinB, sinS = (math.sin(math.radians(a)) for a in (FRONT, BACK, SIDE))


def t_pyramid(x, z):
    """Sharp pyramid: height of the tip surface above its apex."""
    return np.maximum(np.abs(x) * cotS, np.where(z >= 0, z * cotF, -z * cotB))


def t_eroded(x, z, R):
    """Pyramid with every face pushed inward by R -- the allowed ball centres.

    Every face constrains everywhere, so this is the max over all four raised
    planes: unlike the sharp pyramid, the front plane still binds at z < 0
    (it has been raised by R/sin(front)) and the back plane still binds at z > 0.
    """
    return np.maximum(
        np.abs(x) * cotS + R / sinS,
        np.maximum(z * cotF + R / sinF, -z * cotB + R / sinB),
    )


def _disc_offsets(R):
    """Sample points of a disc of radius R, with their sqrt(R^2-|d|^2) lift."""
    g = np.linspace(-R, R, N_DISC)          # odd count -> includes 0 exactly
    dx, dz = np.meshgrid(g, g, indexing="ij")
    dx, dz = dx.ravel(), dz.ravel()
    keep = dx * dx + dz * dz <= R * R
    dx, dz = dx[keep], dz[keep]
    lift = np.sqrt(np.maximum(0.0, R * R - dx * dx - dz * dz))
    # include the rim exactly (lift 0) so the faces come back undistorted
    ang = np.linspace(0, 2 * math.pi, 4 * N_DISC, endpoint=False)
    dx = np.concatenate([dx, R * np.cos(ang)])
    dz = np.concatenate([dz, R * np.sin(ang)])
    lift = np.concatenate([lift, np.zeros_like(ang)])
    return dx, dz, lift


def t_rounded(x, z, R, offsets, chunk=800):
    """Opening of the pyramid by a ball of radius R, evaluated at (x, z)."""
    dx, dz, lift = offsets
    x = np.asarray(x, float).ravel()
    z = np.asarray(z, float).ravel()
    out = np.empty_like(x)
    for i in range(0, x.size, chunk):
        xs = x[i:i + chunk, None] + dx[None, :]
        zs = z[i:i + chunk, None] + dz[None, :]
        out[i:i + chunk] = np.min(t_eroded(xs, zs, R) - lift[None, :], axis=1)
    return out


def find_apex(R, offsets):
    """Lowest point of the rounded surface.

    Front and back differ, so the lowest inscribed ball sits off axis toward the
    blunt back face; x stays on the axis because both side faces share an angle.
    (Strictly the minimum is a short horizontal edge in x, since the ball can
    slide sideways before the side faces bind -- the apex of an asymmetric
    pyramid is a tiny cylindrical fillet, not a point. x = 0 is its centre.)
    """
    z0, win = 0.0, 4.0 * R
    for _ in range(6):
        g = z0 + np.linspace(-win, win, 401)
        y = t_rounded(np.zeros_like(g), g, R, offsets)
        z0 = float(g[int(np.argmin(y))])
        win /= 20.0
    return 0.0, z0, float(t_rounded([0.0], [z0], R, offsets)[0])


def build(height, R):
    """Polar mesh of the tip surface from the apex up to `height`, plus a cap."""
    offsets = _disc_offsets(R)
    x0, z0, apex_y = find_apex(R, offsets)
    surf = lambda x, z: t_rounded(x, z, R, offsets) - apex_y

    thetas = np.linspace(0, 2 * math.pi, N_THETA, endpoint=False)
    cos_t, sin_t = np.cos(thetas), np.sin(thetas)

    # rim radius per direction: bisect surf(r) = height along each ray from the apex
    hi = np.full(N_THETA, height * max(tanS, tanB, tanF) * 1.5 + 4 * R)
    lo = np.zeros(N_THETA)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        too_low = surf(x0 + mid * cos_t, z0 + mid * sin_t) < height
        lo = np.where(too_low, mid, lo)
        hi = np.where(too_low, hi, mid)
    rim = 0.5 * (lo + hi)

    # geometric ring spacing so the apex is finely resolved
    frac = np.geomspace(R / (40.0 * max(rim)), 1.0, N_RING)

    verts = [(x0, 0.0, z0)]                        # apex (lowest point, y = 0)
    ring_start = []
    for f in frac:
        r = rim * f
        x, z = x0 + r * cos_t, z0 + r * sin_t
        y = surf(x, z)
        ring_start.append(len(verts) + 1)
        verts.extend(zip(x.tolist(), y.tolist(), z.tolist()))

    faces = []
    first = ring_start[0]
    for k in range(N_THETA):
        faces.append((1, first + (k + 1) % N_THETA, first + k))
    for i in range(len(ring_start) - 1):
        a, b = ring_start[i], ring_start[i + 1]
        for k in range(N_THETA):
            k2 = (k + 1) % N_THETA
            faces.append((a + k, a + k2, b + k2))
            faces.append((a + k, b + k2, b + k))
    verts.append((x0, height, z0))                 # flat top cap
    centre = len(verts)
    last = ring_start[-1]
    for k in range(N_THETA):
        faces.append((centre, last + k, last + (k + 1) % N_THETA))
    return verts, faces, rim, (x0, z0)


def write_obj(path, header, verts, faces):
    with open(path, "w") as fh:
        for line in header:
            fh.write("# " + line + "\n")
        fh.write("#\n# vertices: %d   triangles: %d\n\n" % (len(verts), len(faces)))
        fh.write("o %s\n" % path[:-4])
        for (x, y, z) in verts:
            fh.write("v %.6f %.6f %.6f\n" % (x, y, z))
        for (a, b, c) in faces:
            fh.write("f %d %d %d\n" % (a, c, b))   # wound CCW seen from outside


def main():
    common = [
        "Bruker SNL-10 AFM probe tip (sharpened silicon pyramid).",
        "Half-angles from vertical: front (+Z) %.1f deg, back (-Z) %.1f deg, side (+/-X) %.1f deg."
        % (FRONT, BACK, SIDE),
        "Apex radius %.1f nm (datasheet: 2 nm nominal, 12 nm max), rounded as the",
        "morphological opening of the pyramid by a ball of that radius, which also",
        "fillets the four edges. Units: nanometres. +Y up, apex at the origin.",
        "Generated by make_snl10_obj.py",
    ]
    for name, height, note in [
        ("snl10_probe.obj", PROBE_H,
         "Whole tip, height %.0f nm (datasheet range 2.5-8.0 um)." % PROBE_H),
        ("snl10_apex.obj", APEX_H,
         "Near-apex cutaway, first %.0f nm, where the %.1f nm rounding is resolvable."
         % (APEX_H, R_APEX)),
    ]:
        verts, faces, rim, apex_xz = build(height, R_APEX)
        hdr = [common[0], note] + common[1:]
        hdr = [h % R_APEX if "%.1f nm (datasheet" in h else h for h in hdr]
        write_obj(name, hdr, verts, faces)
        print("%-16s %6d verts %7d tris   height %7.0f nm   rim radius %.1f-%.1f nm"
              % (name, len(verts), len(faces), height, rim.min(), rim.max()))
        print("                 apex offset from the pyramid axis: x %+.3f nm, z %+.3f nm (toward the blunt back face)" % apex_xz)


if __name__ == "__main__":
    main()
