"""AFM contact-mode forward model (3D grid) — mirrors afm-3d.html physics."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

GW_DEFAULT = 80
GH_DEFAULT = 80
X_NM_DEFAULT = 200.0
Y_NM_DEFAULT = 200.0
TIP_R_NM = 28.0
H_MAX = 95.0
BASE = 25.0
BLOCKED = 1e9
SNL10_ANGLES = {"front": 15.0, "back": 25.0, "side": 22.5}


@dataclass
class InvertState:
    inv_z: bool = False
    mirror_x: bool = False
    mirror_y: bool = False
    inv_tip: bool = False
    erosion: bool = False


@dataclass
class SimConfig:
    tip: str = "sphere"
    tip_params: dict[str, float] = field(default_factory=lambda: {"R": 20.0})
    surface: str = "hemisphere"
    surface_params: dict[str, float] = field(default_factory=lambda: {"R": 30.0})
    gw: int = GW_DEFAULT
    gh: int = GH_DEFAULT
    x_nm: float = X_NM_DEFAULT
    y_nm: float = Y_NM_DEFAULT
    noise_nm: float = 0.0
    noise_seed: int = 42
    inverts: InvertState = field(default_factory=InvertState)
    h_max: float = H_MAX
    h_min: float = 0.0
    tip_kernel_radius_nm: float = TIP_R_NM


def _d_nm(cfg: SimConfig) -> float:
    return cfg.x_nm / cfg.gw


def _tip_r_pixels(cfg: SimConfig) -> int:
    return max(1, int(math.ceil(cfg.tip_kernel_radius_nm / _d_nm(cfg))))


def seeded_rand(i: int, seed: int) -> float:
    x = math.sin(i * 127.1 + seed * 311.7) * 43758.5453
    return x - math.floor(x)


def build_raw_surface(cfg: SimConfig) -> np.ndarray:
    """Build raw surface heights (GH, GW) in nm."""
    gw, gh = cfg.gw, cfg.gh
    d_nm = _d_nm(cfg)
    p = cfg.surface_params
    surf = np.zeros((gh, gw), dtype=np.float64)
    cx, cy = cfg.x_nm / 2, cfg.y_nm / 2

    for r in range(gh):
        for c in range(gw):
            x, y = c * d_nm, r * d_nm
            dx, dy = x - cx, y - cy
            dist = math.hypot(dx, dy)
            match cfg.surface:
                case "sine":
                    a, t = p.get("amp", 15), p.get("period", 55)
                    surf[r, c] = (
                        BASE
                        + a * math.sin(2 * math.pi * x / t) * 0.7
                        + a * 0.3 * math.sin(2 * math.pi * y / t)
                    )
                case "hills":
                    a, s = p.get("amp", 12), p.get("scale", 60)
                    surf[r, c] = BASE + a * (
                        math.sin(2 * math.pi * x / s + 0.5)
                        * math.cos(2 * math.pi * y / (s * 0.8))
                        + 0.5 * math.sin(2 * math.pi * (x + y) / (s * 1.3) + 1.2)
                    )
                case "hemisphere":
                    rad = p.get("R", 25)
                    surf[r, c] = (
                        BASE + math.sqrt(rad * rad - dist * dist)
                        if dist < rad
                        else BASE
                    )
                case "pit":
                    rad = p.get("R", 25)
                    th = BASE + rad
                    surf[r, c] = (
                        th - math.sqrt(rad * rad - dist * dist) if dist < rad else th
                    )
                case "trench":
                    t, dep = p.get("period", 60), p.get("depth", 22)
                    surf[r, c] = BASE if (x % t) / t < 0.45 else BASE + dep
                case "rough":
                    s = p.get("scale", 10) * 0.25
                    surf[r, c] = (
                        BASE
                        + s * 4
                        + s * 3 * math.sin(2 * math.pi * x / 48 + 0.9)
                        + s * 2 * math.cos(2 * math.pi * y / 36 + 1.3)
                        + s * 1.5 * math.sin(2 * math.pi * (x + y) / 22 + 2.1)
                        + s * math.sin(2 * math.pi * x / 11 + 1.4)
                    )
                case "pyramid":
                    h, w = p.get("h", 30), p.get("w", 80) / 2
                    t = max(0.0, 1 - max(abs(dx), abs(dy)) / w)
                    surf[r, c] = BASE + h * t
                case "dna":
                    a, t = p.get("amp", 18), p.get("period", 50)
                    surf[r, c] = BASE + a * abs(math.sin(math.pi * x / t))
                case "chirp":
                    a = p.get("amp", 15)
                    cr = p.get("rate", 8) * 0.0002
                    surf[r, c] = BASE + a * math.sin(2 * math.pi * x * (0.015 + cr * x))
                case "lattice":
                    a, sp = p.get("amp", 10), p.get("spacing", 22)
                    surf[r, c] = (
                        BASE
                        + a
                        * (0.5 + 0.5 * math.cos(2 * math.pi * x / sp))
                        * (0.5 + 0.5 * math.cos(2 * math.pi * y / sp))
                    )
                case "flat":
                    surf[r, c] = BASE
                case "vgrp_15m":
                    pitch = p.get("pitch_nm", 10_000.0)
                    depth = p.get("depth_nm", 180.0)
                    fw = p.get("feature_width_nm", 6000.0)
                    margin = (pitch - fw) / 2.0
                    px = (x + p.get("phase_x_nm", 0.0)) % pitch
                    py = (y + p.get("phase_y_nm", 0.0)) % pitch
                    in_well = margin <= px <= pitch - margin and margin <= py <= pitch - margin
                    surf[r, c] = 0.0 if in_well else depth
                case _:
                    surf[r, c] = BASE
    return surf


def apply_transforms(raw: np.ndarray, cfg: SimConfig) -> np.ndarray:
    inv = cfg.inverts
    flat = raw.ravel().copy()
    n = flat.size

    if cfg.noise_nm > 0:
        for i in range(n):
            flat[i] += (seeded_rand(i, cfg.noise_seed) - 0.5) * 2 * cfg.noise_nm

    surf = np.maximum(2.0, flat).reshape(raw.shape)

    if inv.inv_z:
        mn, mx = surf.min(), surf.max()
        surf = mn + mx - surf

    if inv.mirror_x:
        surf = surf[:, ::-1].copy()

    if inv.mirror_y:
        surf = surf[::-1, :].copy()

    return surf


def tip_height(
    r_nm: float,
    dc_nm: float,
    tip: str,
    p: dict[str, float],
    dr_nm: float = 0.0,
) -> float:
    if tip == "cone":
        alpha = math.radians(p.get("angle", 15))
        return r_nm / math.tan(alpha)
    if tip == "sphere":
        r = p.get("R", 20)
        return r - math.sqrt(r * r - r_nm * r_nm) if r_nm < r else BLOCKED
    if tip == "hyperboloid":
        a = p.get("a", 7)
        return math.sqrt(r_nm * r_nm + a * a) - a
    if tip == "flat":
        r = p.get("R", 18)
        return 0.0 if r_nm <= r else BLOCKED
    if tip == "inv_sphere":
        r = p.get("R", 20)
        lip = max(0.0, p.get("lip", 8))
        outer = r + lip
        if r_nm > outer:
            return BLOCKED
        if r_nm <= r:
            return math.sqrt(max(0.0, r * r - r_nm * r_nm))
        return 0.0
    if tip == "asymmetric":
        alpha = math.radians(p.get("angle", 20))
        tilt = math.radians(p.get("tilt", 15))
        return max(0.0, r_nm / math.tan(alpha) + dc_nm * math.tan(tilt))
    if tip == "double_tip":
        alpha = math.radians(p.get("angle", 15))
        g = p.get("gap", 30) / 2
        r1 = abs(dc_nm - g)
        r2 = abs(dc_nm + g)
        return min(r1, r2) / math.tan(alpha)
    if tip == "sphere_cone":
        r, half_w = p.get("R", 30), p.get("angle", 15)
        sv = math.sqrt(max(1e-6, r * r - half_w * half_w))
        y_tang = r - sv
        if r_nm <= half_w:
            return r - math.sqrt(max(0.0, r * r - r_nm * r_nm))
        return y_tang + (r_nm - half_w) * half_w / sv
    if tip in ("faceted_pyramid", "snl10"):
        return _faceted_pyramid_tip(r_nm, dc_nm, dr_nm, p)
    return r_nm


def _faceted_pyramid_tip(
    r_nm: float, dc_nm: float, dr_nm: float, p: dict[str, float]
) -> float:
    """Asymmetric 4-sided pyramid + rounded apex (Bruker SNL-10 / faceted_pyramid)."""
    r_cap = p.get("R", 2.0)
    front = p.get("front", SNL10_ANGLES["front"])
    back = p.get("back", SNL10_ANGLES["back"])
    side = p.get("side", SNL10_ANGLES["side"])
    cot_f = 1.0 / math.tan(math.radians(front))
    cot_b = 1.0 / math.tan(math.radians(back))
    cot_s = 1.0 / math.tan(math.radians(side))
    tx = abs(dc_nm) * cot_s
    ty = dr_nm * cot_f if dr_nm >= 0 else (-dr_nm) * cot_b
    tp = max(tx, ty)
    if r_nm < r_cap:
        tp = max(tp, r_cap - math.sqrt(max(0.0, r_cap * r_cap - r_nm * r_nm)))
    return tp


def build_tip_kernel(cfg: SimConfig) -> np.ndarray:
    tip_r = _tip_r_pixels(cfg)
    kw = 2 * tip_r + 1
    d_nm = _d_nm(cfg)
    p = cfg.tip_params
    kernel = np.zeros((kw, kw), dtype=np.float64)

    for dr in range(-tip_r, tip_r + 1):
        for dc in range(-tip_r, tip_r + 1):
            r_nm = math.hypot(dc * d_nm, dr * d_nm)
            kernel[dr + tip_r, dc + tip_r] = tip_height(
                r_nm, dc * d_nm, cfg.tip, p, dr_nm=dr * d_nm
            )
    return kernel


def _do_erosion(cfg: SimConfig) -> bool:
    inv = cfg.inverts
    # Contact-mode AFM defaults to dilation for every tip (including concave).
    # User toggles can flip to erosion for teaching / morphological duals.
    return inv.erosion != inv.inv_tip


def compute_measured(surface: np.ndarray, kernel: np.ndarray, cfg: SimConfig) -> np.ndarray:
    gh, gw = surface.shape
    tip_r = _tip_r_pixels(cfg)
    erosion = _do_erosion(cfg)
    measured = np.zeros_like(surface)

    for r in range(gh):
        for c in range(gw):
            best = math.inf if erosion else -math.inf
            for dr in range(-tip_r, tip_r + 1):
                nr = r + dr
                if nr < 0 or nr >= gh:
                    continue
                for dc in range(-tip_r, tip_r + 1):
                    nc = c + dc
                    if nc < 0 or nc >= gw:
                        continue
                    d = kernel[dr + tip_r, dc + tip_r]
                    if d > 1e6:
                        continue
                    v = surface[nr, nc] + d if erosion else surface[nr, nc] - d
                    if erosion:
                        best = min(best, v)
                    else:
                        best = max(best, v)
            measured[r, c] = np.clip(best, cfg.h_min, cfg.h_max)

    return measured


def simulate(cfg: SimConfig) -> tuple[np.ndarray, np.ndarray]:
    """Return (true_surface, measured_afm) height maps in nm."""
    raw = build_raw_surface(cfg)
    surf = apply_transforms(raw, cfg)
    return surf, forward_from_surface(surf, cfg)


def forward_from_surface(surface: np.ndarray, cfg: SimConfig) -> np.ndarray:
    kernel = build_tip_kernel(cfg)
    return compute_measured(surface, kernel, cfg)


def upsample_surface(surface: np.ndarray, factor: int) -> np.ndarray:
    if factor <= 1:
        return surface.astype(np.float64, copy=False)
    from scipy.ndimage import zoom

    return zoom(surface, factor, order=1)


def grid_stats(true: np.ndarray, measured: np.ndarray) -> dict[str, float]:
    err = np.abs(measured - true)
    h_min = float(min(true.min(), measured.min()))
    h_max = float(max(true.max(), measured.max()))
    if h_min == h_max:
        h_max = h_min + 1.0
    return {
        "err_max_nm": float(err.max()),
        "err_rms_nm": float(np.sqrt(np.mean(err * err))),
        "h_min_nm": h_min,
        "h_max_nm": h_max,
    }


def height_to_uint16(h: np.ndarray, h_min: float, h_max: float) -> np.ndarray:
    t = (h - h_min) / (h_max - h_min)
    return np.clip(np.round(t * 65535), 0, 65535).astype(np.uint16)


def decode_uint16_png(arr: np.ndarray, h_min: float, h_max: float) -> np.ndarray:
    return arr.astype(np.float64) / 65535.0 * (h_max - h_min) + h_min


def params_dict(cfg: SimConfig, stats: dict[str, float]) -> dict[str, Any]:
    inv = cfg.inverts
    return {
        "simulator": "3d",
        "tip": cfg.tip,
        "tip_params": dict(cfg.tip_params),
        "surface": cfg.surface,
        "surface_params": dict(cfg.surface_params),
        "grid": {
            "GW": cfg.gw,
            "GH": cfg.gh,
            "x_nm": cfg.x_nm,
            "y_nm": cfg.y_nm,
            "d_nm": _d_nm(cfg),
        },
        "noise_nm": cfg.noise_nm,
        "inverts": {
            "invZ": inv.inv_z,
            "mirrorX": inv.mirror_x,
            "mirrorY": inv.mirror_y,
            "invTip": inv.inv_tip,
            "erosion": inv.erosion,
        },
        "operation": "erosion" if _do_erosion(cfg) else "dilation",
        "png_scale": {
            "h_min_nm": stats["h_min_nm"],
            "h_max_nm": stats["h_max_nm"],
            "decode": "h_nm = uint16 / 65535 * (h_max_nm - h_min_nm) + h_min_nm",
        },
        "stats": {
            "err_max_nm": stats["err_max_nm"],
            "err_rms_nm": stats["err_rms_nm"],
        },
    }
