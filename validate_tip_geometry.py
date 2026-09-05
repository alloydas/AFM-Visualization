"""Validate exported AFM tip meshes against the simulator's analytical geometry.

The simulator uses X/Z as its lateral plane and Y as height above the contact
reference. This script independently parses a Wavefront OBJ, samples its lower
Y envelope, and compares that envelope with the analytical tip profile used by
``afm-3d.html``.

Examples
--------
Run equation-level checks only:

    python validate_tip_geometry.py --self-test

Validate an exported or sample cone:

    python validate_tip_geometry.py samples/tip_cone_15deg.obj \
        --tip cone --param angle=15

Reports are written to ``validation_results/`` by default. Exit status is zero
only when all checks pass.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence


BLOCKED = 1e9
DEFAULT_KERNEL_RADIUS_NM = 30.0
EPS = 1e-9

TIP_DEFAULTS: dict[str, dict[str, float]] = {
    "cone": {"angle": 15.0},
    "sphere": {"R": 20.0},
    "hyperboloid": {"a": 7.0},
    "flat": {"R": 18.0},
    "inv_sphere": {"R": 20.0, "lip": 8.0},
    "asymmetric": {"angle": 20.0, "tilt": 15.0},
    "double_tip": {"angle": 15.0, "gap": 30.0},
    "sphere_cone": {"R": 30.0, "angle": 15.0},
    "triangular_pyramid": {"angle": 15.0},
    "faceted_pyramid": {
        "R": 2.0,
        "front": 15.0,
        "back": 25.0,
        "side": 22.5,
    },
}


@dataclass(frozen=True)
class Vertex:
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class Triangle:
    a: int
    b: int
    c: int


@dataclass
class Sample:
    x_nm: float
    z_nm: float
    expected_y_nm: float
    observed_y_nm: float | None
    error_nm: float | None
    status: str


@dataclass
class Metrics:
    expected_samples: int
    intersected_samples: int
    missing_samples: int
    coverage_percent: float
    mean_error_nm: float | None
    mean_abs_error_nm: float | None
    rmse_nm: float | None
    max_abs_error_nm: float | None
    p95_abs_error_nm: float | None


def parse_params(values: Sequence[str]) -> dict[str, float]:
    parsed: dict[str, float] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"Parameter must use NAME=VALUE syntax: {value!r}")
        name, raw = value.split("=", 1)
        name = name.strip()
        if not name:
            raise ValueError(f"Parameter name is empty: {value!r}")
        parsed[name] = float(raw)
    return parsed


def resolved_params(tip: str, overrides: dict[str, float]) -> dict[str, float]:
    if tip not in TIP_DEFAULTS:
        raise ValueError(f"Unsupported tip {tip!r}")
    result = TIP_DEFAULTS[tip].copy()
    result.update(overrides)
    return result


def validate_params(tip: str, p: dict[str, float]) -> list[str]:
    problems: list[str] = []

    for key in ("R", "a", "gap"):
        if key in p and p[key] <= 0:
            problems.append(f"{key} must be greater than zero")

    for key in ("angle", "front", "back", "side"):
        if key in p and not 0 < p[key] < 90:
            problems.append(f"{key} must be between 0 and 90 degrees")

    if p.get("lip", 0) < 0:
        problems.append("lip must not be negative")

    if tip == "sphere_cone" and p["angle"] >= p["R"]:
        problems.append(
            "sphere_cone tangent half-width (angle parameter) must be smaller "
            "than R; the current UI permits invalid combinations"
        )

    if tip == "asymmetric":
        alpha = math.radians(p["angle"])
        tilt = abs(math.radians(p["tilt"]))
        if alpha + tilt >= math.pi / 2:
            problems.append("angle + abs(tilt) must be less than 90 degrees")

    return problems


def tip_height(tip: str, x_nm: float, z_nm: float, p: dict[str, float]) -> float:
    """Return analytical solid height T(x,z), matching afm-3d.html's kernel."""
    r_nm = math.hypot(x_nm, z_nm)

    if tip == "cone":
        return r_nm / math.tan(math.radians(p["angle"]))

    if tip == "sphere":
        radius = p["R"]
        return (
            radius - math.sqrt(max(0.0, radius * radius - r_nm * r_nm))
            if r_nm < radius
            else BLOCKED
        )

    if tip == "hyperboloid":
        a = p["a"]
        return math.sqrt(r_nm * r_nm + a * a) - a

    if tip == "flat":
        return 0.0 if r_nm <= p["R"] else BLOCKED

    if tip == "inv_sphere":
        radius, outer = p["R"], p["R"] + p.get("lip", 0.0)
        if r_nm > outer:
            return BLOCKED
        if r_nm <= radius:
            return math.sqrt(max(0.0, radius * radius - r_nm * r_nm))
        return 0.0

    if tip == "asymmetric":
        alpha = math.radians(p["angle"])
        tilt = math.radians(p["tilt"])
        return max(0.0, r_nm / math.tan(alpha) + x_nm * math.tan(tilt))

    if tip == "double_tip":
        alpha = math.radians(p["angle"])
        half_gap = p["gap"] / 2
        r1 = math.hypot(x_nm - half_gap, z_nm)
        r2 = math.hypot(x_nm + half_gap, z_nm)
        return min(r1, r2) / math.tan(alpha)

    if tip == "sphere_cone":
        radius, half_width = p["R"], p["angle"]
        side_vertical = math.sqrt(radius * radius - half_width * half_width)
        tangent_y = radius - side_vertical
        if r_nm <= half_width:
            return radius - math.sqrt(
                max(0.0, radius * radius - r_nm * r_nm)
            )
        return tangent_y + (r_nm - half_width) * half_width / side_vertical

    if tip == "triangular_pyramid":
        # A regular three-sided pyramid. At height y, its equilateral
        # cross-section has inradius y*tan(angle). These are the three outward
        # face-normal projections, separated by 120 degrees.
        face_distance = max(
            x_nm,
            -0.5 * x_nm + (math.sqrt(3) / 2) * z_nm,
            -0.5 * x_nm - (math.sqrt(3) / 2) * z_nm,
        )
        return max(0.0, face_distance / math.tan(math.radians(p["angle"])))

    if tip == "faceted_pyramid":
        cot_front = 1 / math.tan(math.radians(p["front"]))
        cot_back = 1 / math.tan(math.radians(p["back"]))
        cot_side = 1 / math.tan(math.radians(p["side"]))
        planar = max(
            abs(x_nm) * cot_side,
            z_nm * cot_front if z_nm >= 0 else -z_nm * cot_back,
        )
        if r_nm < p["R"]:
            spherical = p["R"] - math.sqrt(
                max(0.0, p["R"] ** 2 - r_nm**2)
            )
            return max(planar, spherical)
        return planar

    raise ValueError(f"Unsupported tip {tip!r}")


def parse_obj(path: Path) -> tuple[list[Vertex], list[Triangle]]:
    vertices: list[Vertex] = []
    triangles: list[Triangle] = []

    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8", errors="replace").splitlines(), 1
    ):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()

        if fields[0] == "v" and len(fields) >= 4:
            vertices.append(Vertex(float(fields[1]), float(fields[2]), float(fields[3])))
        elif fields[0] == "f" and len(fields) >= 4:
            face: list[int] = []
            for token in fields[1:]:
                raw_index = token.split("/", 1)[0]
                index = int(raw_index)
                if index < 0:
                    index = len(vertices) + index
                else:
                    index -= 1
                if index < 0 or index >= len(vertices):
                    raise ValueError(
                        f"{path}:{line_number}: face index {raw_index} is out of range"
                    )
                face.append(index)
            for i in range(1, len(face) - 1):
                triangles.append(Triangle(face[0], face[i], face[i + 1]))

    if not vertices:
        raise ValueError(f"{path} contains no OBJ vertices")
    if not triangles:
        raise ValueError(f"{path} contains no OBJ faces")
    return vertices, triangles


def infer_origin(vertices: Sequence[Vertex]) -> tuple[float, float, float]:
    """Infer lateral center and contact-reference Y from lowest vertices."""
    min_y = min(v.y for v in vertices)
    y_span = max(v.y for v in vertices) - min_y
    tolerance = max(1e-5, y_span * 1e-5)
    contact = [v for v in vertices if v.y <= min_y + tolerance]
    return (
        statistics.fmean(v.x for v in contact),
        min_y,
        statistics.fmean(v.z for v in contact),
    )


def normalize_vertices(
    vertices: Sequence[Vertex], origin: tuple[float, float, float]
) -> list[Vertex]:
    ox, oy, oz = origin
    return [Vertex(v.x - ox, v.y - oy, v.z - oz) for v in vertices]


def projected_radius(vertices: Sequence[Vertex]) -> float:
    return max(math.hypot(v.x, v.z) for v in vertices)


def recommended_radius(
    tip: str, p: dict[str, float], mesh_radius: float
) -> float:
    cap = min(DEFAULT_KERNEL_RADIUS_NM, mesh_radius * 0.98)
    if tip == "sphere":
        return min(cap, p["R"] * 0.95)
    if tip == "flat":
        return min(mesh_radius * 0.98, p["R"] * 0.98)
    if tip == "inv_sphere":
        return min(mesh_radius * 0.98, (p["R"] + p.get("lip", 0.0)) * 0.98)
    if tip == "double_tip":
        return min(mesh_radius * 0.98, p["gap"] / 2 + cap / 2)
    if tip == "sphere_cone":
        return min(mesh_radius * 0.98, p["angle"] * 1.9)
    if tip == "triangular_pyramid":
        # Circumradius = 2 * inradius for an equilateral triangle. Restrict
        # the circular comparison region to the interior of the triangular
        # base so every sampled vertical ray intersects the finite mesh.
        return mesh_radius * 0.475
    return cap


def triangle_height_at_xz(
    x: float, z: float, a: Vertex, b: Vertex, c: Vertex
) -> float | None:
    """Intersect a vertical line with one triangle using X/Z barycentrics."""
    denominator = (b.z - c.z) * (a.x - c.x) + (c.x - b.x) * (a.z - c.z)
    if abs(denominator) <= 1e-12:
        return None

    wa = ((b.z - c.z) * (x - c.x) + (c.x - b.x) * (z - c.z)) / denominator
    wb = ((c.z - a.z) * (x - c.x) + (a.x - c.x) * (z - c.z)) / denominator
    wc = 1.0 - wa - wb
    bary_tolerance = 1e-8
    if (
        wa < -bary_tolerance
        or wb < -bary_tolerance
        or wc < -bary_tolerance
    ):
        return None
    return wa * a.y + wb * b.y + wc * c.y


def lower_envelope_height(
    x: float,
    z: float,
    vertices: Sequence[Vertex],
    triangles: Sequence[Triangle],
) -> float | None:
    intersections: list[float] = []
    for tri in triangles:
        height = triangle_height_at_xz(
            x, z, vertices[tri.a], vertices[tri.b], vertices[tri.c]
        )
        if height is not None:
            intersections.append(height)
    return min(intersections) if intersections else None


def sample_coordinates(radius: float, grid_size: int) -> Iterable[tuple[float, float]]:
    if grid_size < 5 or grid_size % 2 == 0:
        raise ValueError("grid-size must be an odd integer of at least 5")
    step = 2 * radius / (grid_size - 1)
    for row in range(grid_size):
        z = -radius + row * step
        for col in range(grid_size):
            x = -radius + col * step
            if math.hypot(x, z) <= radius + EPS:
                yield x, z


def compare_mesh(
    tip: str,
    p: dict[str, float],
    vertices: Sequence[Vertex],
    triangles: Sequence[Triangle],
    radius: float,
    grid_size: int,
) -> list[Sample]:
    samples: list[Sample] = []
    for x, z in sample_coordinates(radius, grid_size):
        expected = tip_height(tip, x, z, p)
        if expected >= BLOCKED / 2:
            continue
        observed = lower_envelope_height(x, z, vertices, triangles)
        if observed is None:
            samples.append(Sample(x, z, expected, None, None, "missing_mesh"))
        else:
            error = observed - expected
            samples.append(Sample(x, z, expected, observed, error, "compared"))
    return samples


def percentile(values: Sequence[float], fraction: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def calculate_metrics(samples: Sequence[Sample]) -> Metrics:
    errors = [s.error_nm for s in samples if s.error_nm is not None]
    absolute = [abs(e) for e in errors]
    expected = len(samples)
    intersected = len(errors)
    missing = expected - intersected
    return Metrics(
        expected_samples=expected,
        intersected_samples=intersected,
        missing_samples=missing,
        coverage_percent=100 * intersected / expected if expected else 0.0,
        mean_error_nm=statistics.fmean(errors) if errors else None,
        mean_abs_error_nm=statistics.fmean(absolute) if absolute else None,
        rmse_nm=math.sqrt(statistics.fmean(e * e for e in errors))
        if errors
        else None,
        max_abs_error_nm=max(absolute) if absolute else None,
        p95_abs_error_nm=percentile(absolute, 0.95) if absolute else None,
    )


def write_samples_csv(path: Path, samples: Sequence[Sample]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "x_nm",
                "z_nm",
                "expected_y_nm",
                "observed_y_nm",
                "signed_error_nm",
                "status",
            ]
        )
        for sample in samples:
            writer.writerow(
                [
                    f"{sample.x_nm:.9g}",
                    f"{sample.z_nm:.9g}",
                    f"{sample.expected_y_nm:.9g}",
                    ""
                    if sample.observed_y_nm is None
                    else f"{sample.observed_y_nm:.9g}",
                    "" if sample.error_nm is None else f"{sample.error_nm:.9g}",
                    sample.status,
                ]
            )


def svg_polyline(
    points: Sequence[tuple[float, float]],
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    width: float,
    height: float,
    margin: float,
) -> str:
    def sx(x: float) -> float:
        return margin + (x - x_min) / max(EPS, x_max - x_min) * (width - 2 * margin)

    def sy(y: float) -> float:
        return height - margin - (y - y_min) / max(EPS, y_max - y_min) * (
            height - 2 * margin
        )

    return " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in points)


def write_profile_svg(path: Path, samples: Sequence[Sample], radius: float) -> None:
    center_step = min((abs(s.z_nm) for s in samples), default=0.0)
    profile = sorted(
        (s for s in samples if abs(abs(s.z_nm) - center_step) <= 1e-7),
        key=lambda s: s.x_nm,
    )
    expected = [(s.x_nm, s.expected_y_nm) for s in profile]
    observed = [
        (s.x_nm, s.observed_y_nm)
        for s in profile
        if s.observed_y_nm is not None
    ]
    all_y = [y for _, y in expected] + [y for _, y in observed]
    y_min = min(all_y, default=0.0)
    y_max = max(all_y, default=1.0)
    if math.isclose(y_min, y_max):
        y_max = y_min + 1.0

    width, height, margin = 800.0, 460.0, 62.0
    expected_points = svg_polyline(
        expected, -radius, radius, y_min, y_max, width, height, margin
    )
    observed_points = svg_polyline(
        observed, -radius, radius, y_min, y_max, width, height, margin
    )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="800" height="460" viewBox="0 0 800 460">
<rect width="800" height="460" fill="white"/>
<text x="400" y="26" text-anchor="middle" font-family="sans-serif" font-size="18">Tip centerline profile comparison</text>
<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="#444"/>
<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="#444"/>
<text x="400" y="448" text-anchor="middle" font-family="sans-serif" font-size="13">Lateral X (nm)</text>
<text x="16" y="230" text-anchor="middle" transform="rotate(-90 16 230)" font-family="sans-serif" font-size="13">Height above contact reference (nm)</text>
<polyline points="{expected_points}" fill="none" stroke="#1769aa" stroke-width="3"/>
<polyline points="{observed_points}" fill="none" stroke="#d1495b" stroke-width="2" stroke-dasharray="7 4"/>
<line x1="585" y1="43" x2="615" y2="43" stroke="#1769aa" stroke-width="3"/>
<text x="622" y="48" font-family="sans-serif" font-size="12">Analytical</text>
<line x1="585" y1="62" x2="615" y2="62" stroke="#d1495b" stroke-width="2" stroke-dasharray="7 4"/>
<text x="622" y="67" font-family="sans-serif" font-size="12">OBJ mesh</text>
<text x="{margin}" y="{height-margin+18}" font-family="sans-serif" font-size="11">{-radius:.3g}</text>
<text x="{width-margin}" y="{height-margin+18}" text-anchor="end" font-family="sans-serif" font-size="11">{radius:.3g}</text>
<text x="{margin-8}" y="{height-margin}" text-anchor="end" font-family="sans-serif" font-size="11">{y_min:.3g}</text>
<text x="{margin-8}" y="{margin+4}" text-anchor="end" font-family="sans-serif" font-size="11">{y_max:.3g}</text>
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def analytical_mesh(
    tip: str,
    p: dict[str, float],
    radius: float,
    grid_size: int = 41,
) -> tuple[list[list[float]], list[int]]:
    """Build a clipped analytical height-field mesh for the 3D report."""
    vertices: list[list[float]] = []
    lookup: dict[tuple[int, int], int] = {}
    step = 2 * radius / (grid_size - 1)

    for row in range(grid_size):
        z = -radius + row * step
        for col in range(grid_size):
            x = -radius + col * step
            if math.hypot(x, z) > radius + EPS:
                continue
            y = tip_height(tip, x, z, p)
            if y >= BLOCKED / 2:
                continue
            lookup[(row, col)] = len(vertices)
            vertices.append([x, y, z])

    indices: list[int] = []
    for row in range(grid_size - 1):
        for col in range(grid_size - 1):
            a = lookup.get((row, col))
            b = lookup.get((row, col + 1))
            c = lookup.get((row + 1, col))
            d = lookup.get((row + 1, col + 1))
            if None not in (a, b, c, d):
                indices.extend([a, c, d, a, d, b])
    return vertices, indices


def physics_kernel_points(
    tip: str,
    p: dict[str, float],
    radius: float = DEFAULT_KERNEL_RADIUS_NM,
    spacing: float = 2.5,
) -> list[list[float]]:
    """Sample the simulator's default 2.5 nm physics-kernel lattice."""
    extent = math.floor(radius / spacing)
    points: list[list[float]] = []
    for row in range(-extent, extent + 1):
        z = row * spacing
        for col in range(-extent, extent + 1):
            x = col * spacing
            if math.hypot(x, z) > radius + EPS:
                continue
            y = tip_height(tip, x, z, p)
            if y < BLOCKED / 2:
                points.append([x, y, z])
    return points


def error_color(error: float, tolerance: float) -> list[float]:
    strength = min(1.0, abs(error) / max(tolerance, EPS))
    pale = 0.9 - 0.75 * strength
    return [1.0, pale, pale] if error >= 0 else [pale, pale, 1.0]


def write_3d_html(
    path: Path,
    tip: str,
    params: dict[str, float],
    vertices: Sequence[Vertex],
    triangles: Sequence[Triangle],
    samples: Sequence[Sample],
    metrics: Metrics,
    tolerance: float,
    validation_radius: float,
) -> None:
    """Write an interactive Three.js comparison of analytical and OBJ geometry."""
    analytical_radius = max(DEFAULT_KERNEL_RADIUS_NM, validation_radius)
    expected_vertices, expected_indices = analytical_mesh(
        tip, params, analytical_radius
    )
    kernel = physics_kernel_points(tip, params)
    compared = [sample for sample in samples if sample.observed_y_nm is not None]
    all_heights = [vertex.y for vertex in vertices]
    all_heights.extend(point[1] for point in expected_vertices)
    y_max = max(all_heights, default=48.0)

    payload = {
        "title": f"{tip} geometry validation",
        "status": (
            "PASS"
            if metrics.max_abs_error_nm is not None
            and metrics.max_abs_error_nm <= tolerance
            else "FAIL"
        ),
        "parameters": params,
        "tolerance": tolerance,
        "metrics": asdict(metrics),
        "cameraScale": max(analytical_radius, y_max * 0.55, 20.0),
        "obj": {
            "positions": [[v.x, v.y, v.z] for v in vertices],
            "indices": [[t.a, t.b, t.c] for t in triangles],
        },
        "analytical": {
            "positions": expected_vertices,
            "indices": expected_indices,
        },
        "kernel": kernel,
        "errors": {
            "positions": [
                [s.x_nm, s.observed_y_nm + max(0.08, tolerance * 0.3), s.z_nm]
                for s in compared
            ],
            "colors": [error_color(s.error_nm or 0.0, tolerance) for s in compared],
        },
    }

    html = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>AFM Tip Geometry — Interactive 3D Validation</title>
<style>
  * { box-sizing: border-box; }
  html, body { width: 100%; height: 100%; margin: 0; overflow: hidden; font-family: Arial, sans-serif; background: #10151d; color: #e6edf3; }
  #view { position: absolute; inset: 0; }
  .panel { position: absolute; z-index: 2; top: 16px; left: 16px; width: min(360px, calc(100vw - 32px)); padding: 16px; background: rgba(15,21,29,.94); border: 1px solid #3a4654; border-radius: 6px; }
  h1 { margin: 0 0 5px; font-size: 18px; }
  .subtitle { color: #9caaba; font-size: 12px; margin-bottom: 12px; }
  .status { display: inline-block; padding: 3px 8px; border-radius: 10px; color: #07110d; background: #4ade80; font-size: 11px; font-weight: 700; }
  .status.fail { background: #fb7185; }
  .metric { display: grid; grid-template-columns: 1fr auto; gap: 8px; margin-top: 6px; font-size: 12px; }
  .metric span:first-child { color: #9caaba; }
  .controls { margin-top: 13px; border-top: 1px solid #303a46; padding-top: 10px; display: grid; gap: 7px; font-size: 12px; }
  label { display: flex; align-items: center; gap: 8px; cursor: pointer; }
  .swatch { width: 18px; height: 3px; display: inline-block; }
  .mesh { background: #f59e0b; }
  .analytic { background: #38bdf8; }
  .kernel { background: #c084fc; }
  .error { background: linear-gradient(90deg,#3b82f6,#fff,#ef4444); }
  .hint { position: absolute; z-index: 2; right: 14px; bottom: 12px; color: #9caaba; font-size: 11px; }
</style>
</head>
<body>
<div id="view"></div>
<section class="panel">
  <h1 id="title"></h1>
  <div class="subtitle">Rendered triangular mesh versus analytical geometry</div>
  <span class="status" id="status"></span>
  <div class="metric"><span>RMSE</span><strong id="rmse"></strong></div>
  <div class="metric"><span>Maximum absolute error</span><strong id="maxError"></strong></div>
  <div class="metric"><span>Mesh coverage</span><strong id="coverage"></strong></div>
  <div class="metric"><span>Acceptance tolerance</span><strong id="tolerance"></strong></div>
  <div class="controls">
    <label><input id="showMesh" type="checkbox" checked><span class="swatch mesh"></span>Rendered OBJ mesh</label>
    <label><input id="showAnalytical" type="checkbox" checked><span class="swatch analytic"></span>Analytical surface wireframe</label>
    <label><input id="showKernel" type="checkbox" checked><span class="swatch kernel"></span>Physics kernel, 2.5 nm lattice</label>
    <label><input id="showErrors" type="checkbox"><span class="swatch error"></span>Signed-error samples</label>
  </div>
</section>
<div class="hint">Drag to orbit · wheel to zoom · right-drag to pan</div>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
<script>
const data = __AFM_VALIDATION_DATA__;
const container = document.getElementById("view");
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x10151d);

const camera = new THREE.PerspectiveCamera(42, innerWidth / innerHeight, 0.05, 2000);
const scale = data.cameraScale;
camera.position.set(scale * 1.65, scale * 1.25, scale * 1.65);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.setSize(innerWidth, innerHeight);
container.appendChild(renderer.domElement);

const controls = new THREE.OrbitControls(camera, renderer.domElement);
controls.target.set(0, scale * 0.45, 0);
controls.enableDamping = true;

scene.add(new THREE.HemisphereLight(0xffffff, 0x263241, 1.25));
const key = new THREE.DirectionalLight(0xffffff, 0.9);
key.position.set(scale, scale * 2, scale);
scene.add(key);

function geometryFrom(dataPart) {
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.Float32BufferAttribute(dataPart.positions.flat(), 3));
  geometry.setIndex(dataPart.indices.flat());
  geometry.computeVertexNormals();
  return geometry;
}

const objMesh = new THREE.Mesh(
  geometryFrom(data.obj),
  new THREE.MeshPhongMaterial({ color: 0xf59e0b, transparent: true, opacity: 0.58, side: THREE.DoubleSide, shininess: 70 })
);
scene.add(objMesh);

const analyticalLines = new THREE.LineSegments(
  new THREE.WireframeGeometry(geometryFrom(data.analytical)),
  new THREE.LineBasicMaterial({ color: 0x38bdf8, transparent: true, opacity: 0.8 })
);
scene.add(analyticalLines);

function pointCloud(positions, color, size, colors) {
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions.flat(), 3));
  if (colors) geometry.setAttribute("color", new THREE.Float32BufferAttribute(colors.flat(), 3));
  return new THREE.Points(
    geometry,
    new THREE.PointsMaterial({ color, size, sizeAttenuation: true, vertexColors: Boolean(colors) })
  );
}

const kernelPoints = pointCloud(data.kernel, 0xc084fc, Math.max(0.28, scale * 0.012));
scene.add(kernelPoints);
const errorPoints = pointCloud(data.errors.positions, 0xffffff, Math.max(0.22, scale * 0.009), data.errors.colors);
errorPoints.visible = false;
scene.add(errorPoints);

const grid = new THREE.GridHelper(scale * 2.8, 20, 0x536273, 0x293442);
grid.position.y = -0.05;
scene.add(grid);
const axes = new THREE.AxesHelper(scale * 0.45);
scene.add(axes);

document.getElementById("title").textContent = data.title;
const status = document.getElementById("status");
status.textContent = data.status;
if (data.status !== "PASS") status.classList.add("fail");
const fmt = value => value == null ? "n/a" : `${value.toFixed(4)} nm`;
document.getElementById("rmse").textContent = fmt(data.metrics.rmse_nm);
document.getElementById("maxError").textContent = fmt(data.metrics.max_abs_error_nm);
document.getElementById("coverage").textContent = `${data.metrics.coverage_percent.toFixed(2)}%`;
document.getElementById("tolerance").textContent = `${data.tolerance.toFixed(3)} nm`;
document.getElementById("showMesh").onchange = event => objMesh.visible = event.target.checked;
document.getElementById("showAnalytical").onchange = event => analyticalLines.visible = event.target.checked;
document.getElementById("showKernel").onchange = event => kernelPoints.visible = event.target.checked;
document.getElementById("showErrors").onchange = event => errorPoints.visible = event.target.checked;

addEventListener("resize", () => {
  camera.aspect = innerWidth / innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(innerWidth, innerHeight);
});

function animate() {
  requestAnimationFrame(animate);
  controls.update();
  renderer.render(scene, camera);
}
animate();
</script>
</body>
</html>
"""
    path.write_text(
        html.replace("__AFM_VALIDATION_DATA__", json.dumps(payload)),
        encoding="utf-8",
    )


def assert_close(name: str, actual: float, expected: float, tolerance: float = 1e-8) -> None:
    if abs(actual - expected) > tolerance:
        raise AssertionError(f"{name}: expected {expected}, got {actual}")


def run_self_checks() -> list[str]:
    """Verify analytical equations at points with independently known results."""
    checks: list[str] = []

    cone = resolved_params("cone", {"angle": 45})
    assert_close("cone apex", tip_height("cone", 0, 0, cone), 0)
    assert_close("cone 45-degree slope", tip_height("cone", 5, 0, cone), 5)
    checks.append("cone apex and slope")

    sphere = resolved_params("sphere", {"R": 10})
    y = tip_height("sphere", 6, 0, sphere)
    assert_close("sphere profile", y, 2)
    assert_close("sphere implicit equation", 6**2 + (y - 10) ** 2, 10**2)
    checks.append("sphere radius invariant")

    hyperboloid = resolved_params("hyperboloid", {"a": 4})
    assert_close(
        "hyperboloid profile", tip_height("hyperboloid", 3, 0, hyperboloid), 1
    )
    checks.append("hyperboloid profile")

    flat = resolved_params("flat", {"R": 5})
    assert_close("flat contact", tip_height("flat", 4, 0, flat), 0)
    if tip_height("flat", 6, 0, flat) < BLOCKED / 2:
        raise AssertionError("flat support boundary was not blocked")
    checks.append("flat contact plane and support")

    bowl = resolved_params("inv_sphere", {"R": 10, "lip": 2})
    assert_close("concave center", tip_height("inv_sphere", 0, 0, bowl), 10)
    assert_close("concave rim", tip_height("inv_sphere", 10, 0, bowl), 0)
    assert_close("concave lip", tip_height("inv_sphere", 11, 0, bowl), 0)
    checks.append("concave sphere center, rim and lip")

    asymmetric = resolved_params("asymmetric", {"angle": 45, "tilt": 10})
    assert_close("asymmetric apex", tip_height("asymmetric", 0, 0, asymmetric), 0)
    if tip_height("asymmetric", 3, 0, asymmetric) <= tip_height(
        "asymmetric", -3, 0, asymmetric
    ):
        raise AssertionError("asymmetric directional slope is reversed")
    checks.append("asymmetric directional slope")

    double = resolved_params("double_tip", {"angle": 45, "gap": 12})
    assert_close("double left apex", tip_height("double_tip", -6, 0, double), 0)
    assert_close("double right apex", tip_height("double_tip", 6, 0, double), 0)
    assert_close("double center", tip_height("double_tip", 0, 0, double), 6)
    checks.append("double-tip apex separation")

    joined = resolved_params("sphere_cone", {"R": 10, "angle": 6})
    tangent = tip_height("sphere_cone", 6, 0, joined)
    assert_close("sphere-cone tangent position", tangent, 2)
    expected_slope = 6 / 8
    delta = 1e-5
    left_slope = (
        tip_height("sphere_cone", 6, 0, joined)
        - tip_height("sphere_cone", 6 - delta, 0, joined)
    ) / delta
    right_slope = (
        tip_height("sphere_cone", 6 + delta, 0, joined)
        - tip_height("sphere_cone", 6, 0, joined)
    ) / delta
    assert_close("sphere-cone left tangent", left_slope, expected_slope, 2e-6)
    assert_close("sphere-cone right tangent", right_slope, expected_slope, 2e-6)
    checks.append("sphere-cone C1 tangent join")

    triangular = resolved_params("triangular_pyramid", {"angle": 45})
    assert_close(
        "triangular pyramid apex",
        tip_height("triangular_pyramid", 0, 0, triangular),
        0,
    )
    assert_close(
        "triangular pyramid positive-X face",
        tip_height("triangular_pyramid", 1, 0, triangular),
        1,
    )
    assert_close(
        "triangular pyramid opposite vertex",
        tip_height("triangular_pyramid", -2, 0, triangular),
        1,
    )
    checks.append("triangular-pyramid face planes")

    pyramid = resolved_params(
        "faceted_pyramid", {"R": 2, "front": 45, "back": 30, "side": 45}
    )
    assert_close("pyramid side", tip_height("faceted_pyramid", 3, 0, pyramid), 3)
    assert_close("pyramid front", tip_height("faceted_pyramid", 0, 3, pyramid), 3)
    assert_close(
        "pyramid back",
        tip_height("faceted_pyramid", 0, -3, pyramid),
        3 / math.tan(math.radians(30)),
    )
    checks.append("faceted-pyramid directional face angles")

    return checks


def print_metrics(metrics: Metrics, tolerance: float, min_coverage: float) -> bool:
    passed = (
        metrics.max_abs_error_nm is not None
        and metrics.max_abs_error_nm <= tolerance
        and metrics.coverage_percent >= min_coverage
    )
    status = "PASS" if passed else "FAIL"
    print(f"\nGeometry validation: {status}")
    print(
        f"  coverage       {metrics.coverage_percent:.2f}% "
        f"({metrics.intersected_samples}/{metrics.expected_samples})"
    )
    if metrics.rmse_nm is not None:
        print(f"  RMSE           {metrics.rmse_nm:.6f} nm")
        print(f"  mean abs error {metrics.mean_abs_error_nm:.6f} nm")
        print(f"  max abs error  {metrics.max_abs_error_nm:.6f} nm")
        print(f"  p95 abs error  {metrics.p95_abs_error_nm:.6f} nm")
        print(f"  mean bias      {metrics.mean_error_nm:.6f} nm")
    print(f"  tolerance      {tolerance:.6f} nm")
    print(f"  min coverage   {min_coverage:.2f}%")
    return passed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare an AFM tip OBJ lower envelope with the analytical geometry "
            "used by afm-3d.html."
        )
    )
    parser.add_argument("obj", nargs="?", type=Path, help="tip OBJ to validate")
    parser.add_argument("--tip", choices=sorted(TIP_DEFAULTS), help="analytical tip type")
    parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="NAME=VALUE",
        help="override a tip parameter; may be repeated",
    )
    parser.add_argument(
        "--radius",
        type=float,
        help="lateral validation radius in nm (default: safe automatic radius)",
    )
    parser.add_argument(
        "--grid-size",
        type=int,
        default=31,
        help="odd Cartesian sample-grid size (default: 31)",
    )
    parser.add_argument(
        "--tolerance-nm",
        type=float,
        default=0.25,
        help="maximum permitted absolute vertical error (default: 0.25)",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=99.0,
        help="minimum percentage of analytical samples covered by the mesh",
    )
    parser.add_argument(
        "--origin",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        help="OBJ contact origin in nm; default infers it from lowest vertices",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_results"),
        help="report directory (default: validation_results)",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="run analytical equation checks; with no OBJ, stop afterward",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.self_test or args.obj is None:
        try:
            checks = run_self_checks()
        except AssertionError as exc:
            print(f"Analytical self-test: FAIL\n  {exc}", file=sys.stderr)
            return 1
        print(f"Analytical self-test: PASS ({len(checks)} checks)")
        for check in checks:
            print(f"  - {check}")
        if args.obj is None:
            return 0

    if args.tip is None:
        parser.error("--tip is required when an OBJ path is provided")
    if not args.obj.is_file():
        parser.error(f"OBJ file does not exist: {args.obj}")
    if args.tolerance_nm < 0:
        parser.error("--tolerance-nm must not be negative")
    if not 0 <= args.min_coverage <= 100:
        parser.error("--min-coverage must be between 0 and 100")

    try:
        params = resolved_params(args.tip, parse_params(args.param))
        problems = validate_params(args.tip, params)
        if problems:
            for problem in problems:
                print(f"Invalid geometry: {problem}", file=sys.stderr)
            return 1

        vertices, triangles = parse_obj(args.obj)
        origin = tuple(args.origin) if args.origin else infer_origin(vertices)
        normalized = normalize_vertices(vertices, origin)
        mesh_radius = projected_radius(normalized)
        radius = (
            args.radius
            if args.radius is not None
            else recommended_radius(args.tip, params, mesh_radius)
        )
        if radius <= 0:
            raise ValueError("validation radius must be greater than zero")
        if radius >= mesh_radius:
            raise ValueError(
                f"validation radius {radius:g} nm reaches beyond the mesh's "
                f"{mesh_radius:g} nm projected radius"
            )

        samples = compare_mesh(
            args.tip,
            params,
            normalized,
            triangles,
            radius,
            args.grid_size,
        )
        metrics = calculate_metrics(samples)
    except (OSError, ValueError) as exc:
        print(f"Validation setup failed: {exc}", file=sys.stderr)
        return 2

    passed = print_metrics(metrics, args.tolerance_nm, args.min_coverage)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{args.obj.stem}__{args.tip}"
    csv_path = args.output_dir / f"{stem}_samples.csv"
    json_path = args.output_dir / f"{stem}_report.json"
    svg_path = args.output_dir / f"{stem}_profile.svg"
    html_path = args.output_dir / f"{stem}_3d.html"
    write_samples_csv(csv_path, samples)
    write_profile_svg(svg_path, samples, radius)
    write_3d_html(
        html_path,
        args.tip,
        params,
        normalized,
        triangles,
        samples,
        metrics,
        args.tolerance_nm,
        radius,
    )

    report = {
        "status": "PASS" if passed else "FAIL",
        "source_obj": str(args.obj.resolve()),
        "tip": args.tip,
        "parameters": params,
        "coordinate_convention": {
            "lateral_plane": "X/Z",
            "height_axis": "Y",
            "units": "nm",
            "normalization_origin_xyz": list(origin),
        },
        "mesh": {
            "vertices": len(vertices),
            "triangles": len(triangles),
            "projected_radius_nm": mesh_radius,
        },
        "sampling": {
            "radius_nm": radius,
            "grid_size": args.grid_size,
            "tolerance_nm": args.tolerance_nm,
            "minimum_coverage_percent": args.min_coverage,
        },
        "metrics": asdict(metrics),
        "artifacts": {
            "samples_csv": str(csv_path.resolve()),
            "profile_svg": str(svg_path.resolve()),
            "interactive_3d_html": str(html_path.resolve()),
        },
    }
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print("\nEvidence written to:")
    print(f"  {json_path}")
    print(f"  {csv_path}")
    print(f"  {svg_path}")
    print(f"  {html_path}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
