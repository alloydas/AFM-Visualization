#!/usr/bin/env python3
"""Make the paper's summary figures from the saved validation outputs (read-only).

    python3 paper/make_figures.py [--results validation_results/realdata] [--out paper/figures]

Writes fig_concept, fig_gratings, fig_tip, fig_cs20ng and fig_polymers as .pdf and .png.
Every number comes from rmse_report/*.csv and cs20ng_icon_saasph10um/*.csv; nothing is refitted.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]

# Reference palette (dataviz skill, light mode): categorical slots in fixed order, chrome inks.
S1, S2, S3, S4, S5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SURF_FILL = "#c3c2b7"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "legend.frameon": False,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.titlecolor": INK, "text.color": INK, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.axisbelow": True, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 1.6, "savefig.dpi": 300,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
})

FOLDERS = {  # folder key -> short label (order = Table 2)
    "vgrp15m_icon_scanasystair_probe4": "15M Icon\nSA-Air",
    "vgrp15m_icon_rtespa300_125": "15M Icon\nRTESPA-300",
    "vgrpgs_resolve_rtespa150": "GS Resolve\nRTESPA-150",
    "vgrpgs_resolve_scanasystair": "GS Resolve\nSA-Air",
    "vgrp15m_resolve_snl10b": "15M Resolve\nSNL10 B",
}


def read(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def f(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def save(fig: plt.Figure, out: Path, name: str) -> None:
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("wrote", out / f"{name}.pdf")


def panel_label(ax: plt.Axes, s: str) -> None:
    ax.text(-0.02, 1.04, s, transform=ax.transAxes, fontweight="bold", fontsize=9, ha="right", va="bottom")


# --------------------------------------------------------------------------- concept
def dilate_1d(s: np.ndarray, t: np.ndarray) -> np.ndarray:
    """I[i] = max_j S[i+j] - T[j], T centred (len odd)."""
    k = len(t) // 2
    pad = np.pad(s, k, mode="edge")
    return np.array([np.max(pad[i:i + len(t)] - t) for i in range(len(s))])


def erode_1d(img: np.ndarray, t: np.ndarray) -> np.ndarray:
    """E[i] = min_j I[i-j] + T[j]."""
    k = len(t) // 2
    pad = np.pad(img, k, mode="edge")
    tr = t[::-1]
    return np.array([np.min(pad[i:i + len(t)] + tr) for i in range(len(img))])


def fig_concept(out: Path) -> None:
    dx = 0.5
    x = np.arange(0, 600 + dx, dx)
    s = np.zeros_like(x)
    s[(x >= 90) & (x < 230)] = 30.0                       # mesa
    s[(x >= 300) & (x < 335)] = -35.0                     # narrow trench
    s += 40.0 * np.exp(-0.5 * ((x - 470) / 4.0) ** 2)     # sharp asperity
    d = np.arange(-160, 160 + dx, dx)

    def sphere_cone(r: float, half_angle_deg: float) -> np.ndarray:
        a = np.radians(half_angle_deg)
        xt = r * np.cos(a)                                  # sphere-cone tangent point
        sph = r - np.sqrt(np.clip(r * r - d * d, 0, None))
        cone = (r - r * np.sin(a)) + (np.abs(d) - xt) / np.tan(a)
        return np.where(np.abs(d) <= xt, sph, cone)

    tips = [("sharp tip ($R$ 5 nm, 20°)", sphere_cone(5.0, 20.0), S1),
            ("blunt tip ($R$ 60 nm, 35°)", sphere_cone(60.0, 35.0), S2)]

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.5), gridspec_kw={"width_ratios": [2.2, 1]})
    ax = axes[0]
    ax.fill_between(x, s, -60, color=SURF_FILL, alpha=0.55, lw=0, label="true surface $S$")
    ax.plot(x, s, color=INK2, lw=0.8)
    blunt = None
    for label, t, col in tips:
        img = dilate_1d(s, t)
        ax.plot(x, img, color=col, label=f"image $I$, {label}")
        if col == S2:
            blunt = (t, img)
    t, img = blunt
    e = erode_1d(img, t)
    ax.fill_between(x, s, e, where=e > s + 0.5, color=S2, alpha=0.18, lw=0,
                    label="unreachable by blunt tip ($E > S$)")
    for xc in (317.5,):                                   # draw the blunt tip over the trench
        i = int(round(xc / dx))
        near = np.abs(d) <= 70
        ax.plot(xc + d[near], img[i] + t[near], color=S2, lw=0.8, ls="--")
    ax.set_xlim(0, 600)
    ax.set_ylim(-45, 95)
    ax.set_xlabel("lateral position (nm)")
    ax.set_ylabel("height (nm)")
    ax.legend(loc="upper left", ncol=2, handlelength=1.6, columnspacing=1.0)
    panel_label(ax, "a")

    ax = axes[1]
    r = np.logspace(0, 4.3, 200)
    h = 20.3
    ax.plot(r, (np.sqrt(1.8) - np.sqrt(0.2)) * np.sqrt(r * h), color=S1, label=r"sphere on step, $0.894\sqrt{Rh}$")
    ax.axhspan(52, 80, color=S3, alpha=0.25, lw=0, label="measured CS-20NG, 52–80 nm")
    ax.plot([1e4], [0.894 * np.sqrt(1e4 * h)], "o", ms=6, color=S2, mec="white", mew=1.2,
            label=r"header sphere, $R$ 10 µm")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, 2e4)
    ax.set_ylim(2, 1000)
    ax.set_xlabel("sphere radius $R$ (nm)")
    ax.set_ylabel("10–90 % edge width (nm)")
    ax.set_title(f"edge width, step $h$ = {h} nm", loc="left")
    ax.legend(loc="lower right")
    panel_label(ax, "b")
    fig.tight_layout(w_pad=2.0)
    save(fig, out, "fig_concept")


# --------------------------------------------------------------------------- gratings
def fig_gratings(rep: Path, out: Path) -> None:
    rows = [r for r in read(rep / "rmse_by_scan.csv")
            if r["folder_key"] in FOLDERS and r["role"] == "predict" and r["primary_held_out"] == "True"]
    rows.sort(key=lambda r: (list(FOLDERS).index(r["folder_key"]), r["scan_file"]))
    labels, groups = [], []
    for r in rows:
        labels.append(r["scan_file"][-9:-4])
        groups.append(r["folder_key"])
    xs = np.arange(len(rows), dtype=float)
    # gap between folders
    for i in range(1, len(rows)):
        if groups[i] != groups[i - 1]:
            xs[i:] += 3.0

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.2), gridspec_kw={"width_ratios": [1, 1]})
    ax = axes[0]
    series = [("M0 (nominal preset)", "M0_rmse_nm", S2, "s"),
              ("M1 (6 µm wells)", "M1_rmse_nm", S4, "D"),
              ("M2 prior (sharp corners)", "M2_prior_rmse_nm", S5, "^"),
              ("M2 held-out", "primary_rmse_nm", S1, "o")]
    for name, col, c, m in series:
        ax.plot(xs, [f(r[col]) for r in rows], m, ms=5.5, color=c, mec="white", mew=0.8, label=name, ls="none")
    ax.set_yscale("log")
    ax.set_ylim(1.5, 80)
    ax.set_ylabel("RMSE on held-out frame (nm)")
    panel_label(ax, "a")

    ax = axes[1]
    w = 0.36
    ax.bar(xs - w / 2, [f(r["primary_flat_rmse_nm"]) for r in rows], w, color=S3, label="flat pixels")
    ax.bar(xs + w / 2, [f(r["primary_edge_rmse_nm"]) for r in rows], w, color=S1, label="wall band")
    ax.plot(xs - w / 2, [f(r["noise_nm"]) for r in rows], "_", ms=9, mew=1.6, color=INK, label="flat-level scatter")
    ax.set_ylabel("M2 held-out RMSE (nm)")
    ax.set_ylim(0, 23)
    ax.legend(loc="upper left", ncol=3, handlelength=1.4, columnspacing=1.0)
    panel_label(ax, "b")
    axes[0].set_ylim(1.5, 150)
    axes[0].legend(loc="upper left", ncol=2, handletextpad=0.3, columnspacing=1.0)

    for ax in axes:
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, rotation=90, fontsize=6.5)
        ax.grid(axis="x", visible=False)
        trans = ax.get_xaxis_transform()
        for key in FOLDERS:
            sel = [x for x, g in zip(xs, groups) if g == key]
            if sel:
                ax.text(np.mean(sel), -0.25, FOLDERS[key], transform=trans, ha="center", va="top",
                        fontsize=5.6, color=INK2, linespacing=1.0)
    fig.tight_layout(w_pad=2.0)
    fig.subplots_adjust(bottom=0.32)
    save(fig, out, "fig_gratings")


# --------------------------------------------------------------------------- tip sensitivity
def fig_tip(rep: Path, out: Path) -> None:
    full = read(rep / "rmse_full.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6), sharey=False)
    for ax, key, title in [(axes[0], "vgrp15m_icon_rtespa300_125", "VGRP-15M Icon, RTESPA-300 (3 held-out frames)"),
                           (axes[1], "vgrp15m_icon_scanasystair_probe4", "VGRP-15M Icon, ScanAsyst-Air (3 held-out frames)")]:
        predict = {r["scan_file"] for r in full if r["folder_key"] == key and r["role"] == "predict"}
        for model, col, name in [("M0", S2, "M0: nominal geometry held fixed"),
                                 ("M2_holdout", S1, "M2 held-out: geometry recalibrated per probe")]:
            acc: dict[float, list[float]] = defaultdict(list)
            base: list[float] = []
            for r in full:
                if r["folder_key"] != key or r["model"] != model or r["scan_file"] not in predict:
                    continue
                if r["probe_key"] == "no_tip":
                    base.append(f(r["rmse_nm"]))
                elif r["yaw"] == "+90":
                    acc[f(r["probe_R_nm"])].append(f(r["rmse_nm"]))
            b = float(np.mean(base))
            rr = sorted(acc)
            dv = [np.mean(acc[x]) - b for x in rr]
            ax.plot(rr, dv, "-o", ms=4.5, color=col, mec="white", mew=0.8, label=name)
        ax.axhline(0, color=AXIS, lw=0.8)
        ax.set_xscale("log")
        ax.set_xlabel("probe apex radius $R$ (nm), yaw +90°")
        ax.set_title(title, loc="left")
        ticks = sorted({x for x in (2, 5, 8, 12, 30, 60, 125, 250) if x <= max(rr)})
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(x) for x in ticks])
        ax.minorticks_off()
    axes[0].set_ylabel("RMSE change vs no-tip control (nm)")
    axes[0].legend(loc="center left", bbox_to_anchor=(0.0, 0.62))
    panel_label(axes[0], "a")
    panel_label(axes[1], "b")
    fig.tight_layout(w_pad=2.0)
    save(fig, out, "fig_tip")


# --------------------------------------------------------------------------- CS-20NG
def fig_cs20ng(res: Path, out: Path) -> None:
    d = res / "cs20ng_icon_saasph10um"
    esf = read(d / "cs20ng_esf.csv")
    prof = read(d / "cs20ng_r_profiles.csv")
    file0 = "CaptureFile.0_00000.spm"
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.5))
    for ax, fam, title in [(axes[0], "slow", "slow-axis edges, frame 00000"),
                           (axes[1], "fast", "fast-axis edges, frame 00000")]:
        sel = [r for r in esf if r["file"] == file0 and r["family"] == fam and r["variant"] == "resp_vertical"]
        for lab, col, name, style in [("R10000", S2, "simulated, header sphere $R$ 10 µm", "-"),
                                      ("best", S1, "simulated, joint fit ($R$ ≈ 95 nm)", "-")]:
            pts = sorted((f(r["d_nm"]), f(r["mean_nm"])) for r in sel if r["label"] == lab)
            ax.plot(*zip(*pts), style, color=col, label=name)
        m = sorted((f(r["d_nm"]), f(r["mean_nm"]), f(r["sem_nm"])) for r in sel if r["label"] == "measured")
        mx, my, me = map(np.array, zip(*m))
        ax.errorbar(mx, my, yerr=me, fmt="o", ms=2.6, color=INK, elinewidth=0.6, label="measured (mean ± SE)",
                    zorder=3)
        ax.set_xlim(-1000, 1000)
        ax.set_xlabel("distance to registered edge (nm)")
        ax.set_title(title, loc="left")
    axes[0].set_ylabel("height (nm)")
    axes[0].legend(loc="lower left", fontsize=6.3)
    ax = axes[2]
    for fn, col in [("CaptureFile.0_00000.spm", S1), ("CaptureFile.0_00001.spm", S3), ("CaptureFile.0_00002.spm", S4)]:
        pts = sorted((f(r["R_nm"]), f(r["t_vs_best"])) for r in prof
                     if r["file"] == fn and r["variant"] == "resp_vertical" and f(r["R_nm"]) > 0)
        ax.plot(*zip(*pts), color=col, label=f"frame {fn[-9:-4]}")
    ax.axhline(1.645, color=INK2, lw=0.9, ls="--")
    ax.text(30, 2.1, "$t$ = 1.645", color=INK2, fontsize=6.5)
    ax.axvline(1e4, color=S2, lw=0.9, ls=":")
    ax.set_xscale("log")
    ax.set_xlim(25, 2e4)
    ax.set_xlabel("sphere radius $R$ (nm)")
    ax.set_ylabel("paired $t$ vs best fit")
    ax.set_title("joint fits, fitted fast-axis response", loc="left")
    ax.legend(loc="upper left")
    for a, s in zip(axes, "abc"):
        panel_label(a, s)
    fig.tight_layout(w_pad=1.4)
    save(fig, out, "fig_cs20ng")


# --------------------------------------------------------------------------- polymers
def fig_polymers(rep: Path, out: Path) -> None:
    samples = [("fsilica_icon_scanasystair_probe2", "FSILICA", S1),
               ("psfilm_icon_scanasystair_probe2", "PSFILM", S2),
               ("psldpe_resolve_scanasystair", "PS-LDPE", S3)]
    stats = read(rep / "polymer_stats.csv")
    cross = [r for r in read(rep / "rmse_crossscale.csv") if r["is_primary"] in ("True", "true", "1")]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    ax = axes[0]
    for key, name, col in samples:
        sel = [r for r in stats if r["folder_key"] == key and r["source"] == "nominal"]
        px = [f(r["pixel_nm"]) for r in sel]
        rat = [f(r["sq_bg_ratio_noise_matched"]) for r in sel]
        ax.plot(px, rat, "o", ms=5, color=col, mec="white", mew=0.8, label=name, ls="none")
    ax.axhline(1, color=INK2, lw=0.9, ls="--")
    ax.text(60, 1.08, "preset = measured", color=INK2, fontsize=6.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("pixel size (nm)")
    ax.set_ylabel("background $S_q$, preset / measured")
    ax.set_title("nominal random-surface presets", loc="left")
    ax.legend(loc="upper right")
    panel_label(ax, "a")

    ax = axes[1]
    xpos = 0.0
    ticks, ticklabels = [], []
    for key, name, col in samples:
        sel = [r for r in cross if r["folder_key"] == key]
        for r in sel:
            ax.bar(xpos, f(r["rmse_nm"]), 0.7, color=col)
            ax.plot(xpos, f(r["parent_noise_nm"]), "_", ms=8, mew=1.6, color=INK)
            ticks.append(xpos)
            ticklabels.append(f"{r['child'][-9:-4]}→{r['parent'][-9:-4]}")
            xpos += 1
        ax.text(xpos - len(sel) / 2 - 0.5, 1.72, name, ha="center", color=col, fontsize=7, fontweight="bold")
        xpos += 0.8
    ax.plot([], [], "_", ms=8, mew=1.6, color=INK, label="parent noise level")
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels, rotation=90, fontsize=5.5)
    ax.grid(axis="x", visible=False)
    ax.set_ylim(0, 1.85)
    ax.set_ylabel("zoom → parent RMSE (nm)")
    ax.set_title("cross-scale: measured zoom through the engine", loc="left")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.9))
    panel_label(ax, "b")
    fig.tight_layout(w_pad=2.0)
    save(fig, out, "fig_polymers")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results", type=Path, default=REPO / "validation_results" / "realdata")
    ap.add_argument("--out", type=Path, default=REPO / "paper" / "figures")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    rep = a.results / "rmse_report"
    fig_concept(a.out)
    fig_gratings(rep, a.out)
    fig_tip(rep, a.out)
    fig_cs20ng(a.results, a.out)
    fig_polymers(rep, a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
