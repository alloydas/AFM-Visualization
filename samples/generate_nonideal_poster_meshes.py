"""Generate non-ideal OBJ meshes for AFM simulator demonstrations.

Coordinates are nanometres and use the simulator's convention: X/Z form the
lateral plane and Y is height.  Run from the repository root:

    python samples/generate_nonideal_poster_meshes.py
"""

from __future__ import annotations

import math
from pathlib import Path


OUT = Path(__file__).parent


def write_obj(path: Path, name: str, vertices: list[tuple[float, float, float]], faces: list[tuple[int, int, int]], comment: str) -> None:
    lines = [
        f"# {comment}",
        "# Coordinates are nanometres; Y is height.",
        f"o {name}",
    ]
    lines.extend(f"v {x:.5f} {y:.5f} {z:.5f}" for x, y, z in vertices)
    lines.extend(f"f {a} {b} {c}" for a, b, c in faces)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_subtle_worn_tip() -> None:
    """Make the original gently worn, tilted, asymmetric tip."""
    radial_steps = 18
    angular_steps = 48
    radius_nm = 38.0
    vertices: list[tuple[float, float, float]] = []

    for ir in range(radial_steps + 1):
        u = ir / radial_steps
        radius = u * radius_nm
        for it in range(angular_steps):
            theta = 2 * math.pi * it / angular_steps
            x = radius * math.cos(theta)
            z = radius * math.sin(theta)
            elliptical_radius = math.sqrt((x / 1.10) ** 2 + (z / 0.82) ** 2)
            cone = 0.68 * elliptical_radius
            tilt = 0.13 * x - 0.05 * z
            facet = 1.25 * math.sin(3 * theta + 0.4) * u**2
            shoulder = 5.5 * math.exp(-((x - 17) ** 2 / 72 + (z + 5) ** 2 / 38)) * u
            blunt = 2.1 * (1 - math.exp(-(radius / 5.8) ** 2))
            vertices.append((x, max(0.0, cone + tilt + facet + shoulder + blunt), z))

    faces: list[tuple[int, int, int]] = []
    for ir in range(radial_steps):
        for it in range(angular_steps):
            nxt = (it + 1) % angular_steps
            a = ir * angular_steps + it + 1
            b = ir * angular_steps + nxt + 1
            c = (ir + 1) * angular_steps + it + 1
            d = (ir + 1) * angular_steps + nxt + 1
            faces.extend(((a, c, d), (a, d, b)))

    write_obj(
        OUT / "tip_subtle_worn_tilted.obj",
        "tip_subtle_worn_tilted",
        vertices,
        faces,
        "Non-ideal AFM tip: subtle asymmetric wear, tilt, facets, and worn shoulder.",
    )


def make_worn_asymmetric_tip() -> None:
    """Make a visibly damaged split-apex tip for poster demonstrations."""
    radial_steps = 18
    angular_steps = 48
    radius_nm = 38.0
    vertices: list[tuple[float, float, float]] = []

    for ir in range(radial_steps + 1):
        u = ir / radial_steps
        radius = u * radius_nm
        for it in range(angular_steps):
            theta = 2 * math.pi * it / angular_steps
            x = radius * math.cos(theta)
            z = radius * math.sin(theta)

            # A damaged probe often has more than one contacting asperity.
            # The primary apex is sharp and offset left; the secondary apex is
            # blunter, lower-resolution, and offset right.  Taking the lower
            # envelope creates a clearly visible split/double-apex tip.
            primary = 0.78 * math.sqrt((x + 13) ** 2 + (z + 2) ** 2)
            secondary = 4.0 + 0.48 * math.sqrt((x - 15) ** 2 + (z - 5) ** 2)
            twin_apex = min(primary, secondary)

            # Wear and a tilted fractured flank make the two lobes asymmetric.
            tilt = 0.10 * x - 0.07 * z
            facet = 2.8 * math.sin(3 * theta + 0.4) * u**2
            fractured_flank = 8.0 * math.exp(
                -((x - 24) ** 2 / 105 + (z + 11) ** 2 / 55)
            ) * u
            y = max(0.0, twin_apex + tilt + facet + fractured_flank)
            vertices.append((x, y, z))

    faces: list[tuple[int, int, int]] = []
    for ir in range(radial_steps):
        for it in range(angular_steps):
            nxt = (it + 1) % angular_steps
            a = ir * angular_steps + it + 1
            b = ir * angular_steps + nxt + 1
            c = (ir + 1) * angular_steps + it + 1
            d = (ir + 1) * angular_steps + nxt + 1
            faces.extend(((a, c, d), (a, d, b)))

    write_obj(
        OUT / "tip_worn_tilted_asymmetric.obj",
        "tip_worn_tilted_asymmetric",
        vertices,
        faces,
        "Damaged AFM tip: split/double apex, asymmetric wear, tilt, and fractured flank.",
    )


def make_rough_nanocluster_surface() -> None:
    """Make a high-relief synthetic surface with obvious non-ideal features."""
    n = 61
    extent_nm = 100.0
    base_nm = 18.0
    vertices: list[tuple[float, float, float]] = []

    clusters = [
        (-45, -30, 23, 13, 18),
        (-13, 31, 39, 15, 25),
        (26, -5, 52, 18, 14),
        (54, 36, 22, 10, 16),
        (59, -49, 17, 14, 10),
        (-64, 47, 15, 12, 10),
    ]

    for row in range(n):
        z = -extent_nm + 2 * extent_nm * row / (n - 1)
        for col in range(n):
            x = -extent_nm + 2 * extent_nm * col / (n - 1)
            y = base_nm

            # Multi-scale background roughness.
            y += 2.6 * math.sin(0.105 * x + 0.4) * math.cos(0.13 * z - 0.7)
            y += 1.5 * math.sin(0.31 * x - 0.18 * z)
            y += 0.9 * math.cos(0.47 * z + 0.08 * x)

            # Irregular nanoparticles / agglomerates.
            for cx, cz, amp, sx, sz in clusters:
                dx = x - cx
                dz = z - cz
                y += amp * math.exp(-(dx * dx / (2 * sx * sx) + dz * dz / (2 * sz * sz)))

            # A raised, irregular terrace in the upper-left.
            terrace_gate = 1 / (1 + math.exp(-(x + 32) / 3.3))
            terrace_gate *= 1 / (1 + math.exp((x - 4) / 3.3))
            terrace_gate *= 1 / (1 + math.exp(-(z - 25) / 3.3))
            terrace_gate *= 1 / (1 + math.exp((z - 70) / 3.3))
            y += 19 * terrace_gate

            # Two broad pits and a strong diagonal scratch.
            y -= 17 * math.exp(-((x + 18) ** 2 / 170 + (z + 43) ** 2 / 115))
            y -= 12 * math.exp(-((x - 62) ** 2 / 95 + (z - 23) ** 2 / 130))
            scratch_distance = abs(z - 0.34 * x + 8) / math.sqrt(1.0 + 0.34**2)
            y -= 8.0 * math.exp(-(scratch_distance / 3.3) ** 2)
            vertices.append((x, max(2.0, y), z))

    faces: list[tuple[int, int, int]] = []
    for row in range(n - 1):
        for col in range(n - 1):
            a = row * n + col + 1
            b = a + 1
            c = a + n
            d = c + 1
            faces.extend(((a, c, d), (a, d, b)))

    write_obj(
        OUT / "surface_rough_nanoclusters.obj",
        "surface_rough_nanoclusters",
        vertices,
        faces,
        "High-relief non-ideal sample: roughness, agglomerates, raised terrace, pits, and scratch.",
    )


if __name__ == "__main__":
    make_subtle_worn_tip()
    make_worn_asymmetric_tip()
    make_rough_nanocluster_surface()
    print(
        "Generated tip_subtle_worn_tilted.obj, "
        "tip_worn_tilted_asymmetric.obj, and surface_rough_nanoclusters.obj"
    )
