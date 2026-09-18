"""Compare forward-model prediction to real Bruker .spm measurements on calibration samples."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from afm_forward_model import SimConfig, forward_from_surface, upsample_surface
from spm_io import get_channel, load_spm, subtract_plane

ROOT = Path(__file__).resolve().parent
SPECS_PATH = ROOT / "samples" / "calibration_specs.json"
DEFAULT_OUT = ROOT / "validation_results" / "vgrp15m_snl10"


def load_specs() -> dict:
    with SPECS_PATH.open(encoding="utf-8") as f:
        return json.load(f)


def estimate_grating_phase_nm(
    profile: np.ndarray,
    pitch_nm: float,
    depth_nm: float,
    d_nm: float,
    feature_width_nm: float,
) -> float:
    """Find phase offset (nm) along one axis via cross-correlation with an ideal square grating."""
    n = profile.size
    t = np.arange(n, dtype=np.float64) * d_nm
    margin = (pitch_nm - feature_width_nm) / 2.0
    best_phase, best_score = 0.0, -np.inf
    for phase in np.linspace(0, pitch_nm, num=48, endpoint=False):
        px = (t + phase) % pitch_nm
        in_well = (px >= margin) & (px <= pitch_nm - margin)
        template = np.where(in_well, 0.0, depth_nm)
        template = template - template.mean()
        prof = profile - profile.mean()
        score = float(np.dot(prof, template))
        if score > best_score:
            best_score = score
            best_phase = phase
    return best_phase


def build_vgrp_surface(
    gh: int,
    gw: int,
    x_nm: float,
    y_nm: float,
    pitch_nm: float,
    depth_nm: float,
    feature_width_nm: float,
    phase_x_nm: float,
    phase_y_nm: float,
) -> np.ndarray:
    d_x = x_nm / gw
    d_y = y_nm / gh
    margin = (pitch_nm - feature_width_nm) / 2.0
    surf = np.zeros((gh, gw), dtype=np.float64)
    for r in range(gh):
        y = r * d_y
        py = (y + phase_y_nm) % pitch_nm
        for c in range(gw):
            x = c * d_x
            px = (x + phase_x_nm) % pitch_nm
            in_well = margin <= px <= pitch_nm - margin and margin <= py <= pitch_nm - margin
            surf[r, c] = 0.0 if in_well else depth_nm
    return surf


def profile_metrics(true: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    err = pred - true
    return {
        "rmse_nm": float(np.sqrt(np.mean(err * err))),
        "mae_nm": float(np.mean(np.abs(err))),
        "max_abs_nm": float(np.max(np.abs(err))),
    }


def save_profile_png(path: Path, x_nm: np.ndarray, measured: np.ndarray, sim: np.ndarray) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(x_nm, measured, label="SPM (detrended)", lw=1.2)
    ax.plot(x_nm, sim, label="Forward model", lw=1.2, alpha=0.85)
    ax.set_xlabel("Fast axis (nm)")
    ax.set_ylabel("Height (nm)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_case(
    spm_path: Path,
    channel: str,
    upsample: int,
    out_dir: Path,
) -> dict:
    specs = load_specs()
    vgrp = specs["vgrp_15m"]
    tip_spec = specs["snl10_b"]["forward_model_defaults"]

    spm = load_spm(spm_path)
    img = get_channel(spm, channel)
    measured = subtract_plane(img.data.copy())

    gh, gw = measured.shape
    x_nm = img.scan_size_x_nm
    y_nm = img.scan_size_y_nm
    pitch_nm = float(vgrp["pitch_nm"])
    depth_nm = float(vgrp["depth_nm"])
    d_x = x_nm / gw
    d_y = y_nm / gh
    fw = float(vgrp["feature_width_nm"])

    row_mean = measured.mean(axis=1)
    phase_y_nm = estimate_grating_phase_nm(row_mean, pitch_nm, depth_nm, d_y, fw)
    col_mean = measured.mean(axis=0)
    phase_x_nm = estimate_grating_phase_nm(col_mean, pitch_nm, depth_nm, d_x, fw)

    true_s = build_vgrp_surface(
        gh,
        gw,
        x_nm,
        y_nm,
        pitch_nm,
        depth_nm,
        float(vgrp["feature_width_nm"]),
        phase_x_nm,
        phase_y_nm,
    )

    cfg = SimConfig(
        tip="snl10",
        tip_params=dict(tip_spec),
        surface="vgrp_15m",
        gw=gw,
        gh=gh,
        x_nm=x_nm,
        y_nm=y_nm,
        h_max=500.0,
        h_min=-50.0,
        tip_kernel_radius_nm=40.0,
    )

    if upsample > 1:
        true_f = upsample_surface(true_s, upsample)
        cfg_f = SimConfig(
            tip=cfg.tip,
            tip_params=cfg.tip_params,
            surface=cfg.surface,
            gw=gw * upsample,
            gh=gh * upsample,
            x_nm=x_nm,
            y_nm=y_nm,
            h_max=cfg.h_max,
            h_min=cfg.h_min,
            tip_kernel_radius_nm=cfg.tip_kernel_radius_nm,
        )
        sim_f = forward_from_surface(true_f, cfg_f)
        from scipy.ndimage import zoom

        simulated = zoom(sim_f, 1.0 / upsample, order=1)
    else:
        simulated = forward_from_surface(true_s, cfg)

    simulated = subtract_plane(simulated)
    z_off = float(np.median(measured - simulated))
    simulated_aligned = simulated + z_off

    metrics = profile_metrics(measured, simulated_aligned)
    prof = measured.ravel()
    simv = simulated_aligned.ravel()
    fit = np.linalg.lstsq(
        np.column_stack([simv, np.ones_like(simv)]), prof, rcond=None
    )[0]
    scale, offset = float(fit[0]), float(fit[1])
    scaled_rmse = float(
        np.sqrt(np.mean((prof - (scale * simv + offset)) ** 2))
    )
    metrics["linear_scale"] = scale
    metrics["linear_offset_nm"] = offset
    metrics["rmse_after_linear_fit_nm"] = scaled_rmse
    metrics["measured_step_p90_p10_nm"] = float(
        np.percentile(measured, 90) - np.percentile(measured, 10)
    )
    metrics["simulated_step_p90_p10_nm"] = float(
        np.percentile(simulated_aligned, 90) - np.percentile(simulated_aligned, 10)
    )
    mid = gh // 2
    prof_m = measured[mid, :]
    prof_s = simulated_aligned[mid, :]
    x_axis = np.linspace(0, x_nm, gw, endpoint=False)

    out_dir.mkdir(parents=True, exist_ok=True)
    save_profile_png(out_dir / "midline_profile.png", x_axis, prof_m, prof_s)

    report = {
        "spm_file": str(spm_path),
        "channel": img.name,
        "scan_size_nm": [x_nm, y_nm],
        "grid": [gh, gw],
        "upsample": upsample,
        "tip": {"type": "snl10", "params": tip_spec},
        "surface": {
            "type": "vgrp_15m",
            "pitch_nm": pitch_nm,
            "depth_nm": depth_nm,
            "feature_width_nm": vgrp["feature_width_nm"],
            "phase_x_nm": phase_x_nm,
            "phase_y_nm": phase_y_nm,
        },
        "alignment": {"z_offset_nm": z_off},
        "metrics_vs_measured": metrics,
        "caveats": [
            "Instrument: Peak Force / Height Sensor — geometric dilation only.",
            "VGRP layout uses 6 µm feature width where PDF lists pitch/depth only.",
        ],
    }
    with (out_dir / "report.json").open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    np.savetxt(out_dir / "measured_midline.csv", np.column_stack([x_axis, prof_m]), delimiter=",")
    np.savetxt(out_dir / "simulated_midline.csv", np.column_stack([x_axis, prof_s]), delimiter=",")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        default="vgrp15m_snl10",
        help="Key in samples/calibration_specs.json cases section",
    )
    parser.add_argument("--spm", type=Path, default=None, help="Override .spm path")
    parser.add_argument("--upsample", type=int, default=4, help="Grid upsample factor")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    specs = load_specs()
    case = specs["cases"].get(args.case)
    if case is None:
        raise SystemExit(f"Unknown case {args.case!r}")

    spm_path = args.spm or (ROOT / case["spm_glob"])
    if not spm_path.is_file():
        raise SystemExit(f"SPM not found: {spm_path}")

    report = run_case(
        spm_path,
        case["channel_substring"],
        max(1, args.upsample),
        args.out,
    )
    m = report["metrics_vs_measured"]
    print(f"Wrote {args.out / 'report.json'}")
    print(f"RMSE vs measured: {m['rmse_nm']:.2f} nm  MAE: {m['mae_nm']:.2f} nm")
    if "rmse_after_linear_fit_nm" in m:
        print(
            f"After linear scale fit (diagnostic): RMSE {m['rmse_after_linear_fit_nm']:.2f} nm, "
            f"scale={m['linear_scale']:.3f}"
        )
    if "measured_step_p90_p10_nm" in m:
        print(
            f"Step spread p90-p10 - measured: {m['measured_step_p90_p10_nm']:.1f} nm, "
            f"sim: {m['simulated_step_p90_p10_nm']:.1f} nm (spec depth {180} nm)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
