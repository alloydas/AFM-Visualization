# AFM Visualization

Interactive AFM tip–sample simulators (2D / 3D) plus mesh generators and (on the validation line) forward-model checks against real SPM data.

## Layout

| Path | Role |
|------|------|
| `index.html`, `afm-2d.html`, `afm-3d.html`, `README.html` | Deployed web app (Vercel routes these at the repo root) |
| `vercel.json` | Static hosting routes |
| `afm_forward_model.py` | Shared forward-model physics (imported by validation scripts) |
| `tools/` | Mesh generators (`make_*.py`) and, when merged, validators |
| `meshes/{tips,surfaces,scans}/` | Generated Wavefront `.obj` meshes |
| `samples/` | Calibration assets, datasheets, and experimental SPM folders (**pipeline contract** — do not rename) |
| `validation_results/` | Validation outputs (**pipeline contract**; gitignored except `.gitkeep`) |
| `data/raw/` | Large raw SPM archives (e.g. TGXYZ02); gitignored |
| `docs/` | Paper drafts, notes, presentation material |

## Quick start (web)

Open `index.html` locally, or deploy the repo with Vercel (no build step beyond static HTML).

## Mesh generators

From the repo root:

```bash
python tools/make_snl10_obj.py
python tools/make_vgrp_obj.py
python tools/make_cs20ng_obj.py
```

Outputs land under `meshes/`. SPM → mesh conversion reads `data/raw/TGXYZ02/` by default:

```bash
python tools/make_spm_obj.py
```

## Validation paths (do not break)

Case definitions and SPM globs live under `samples/` (including `samples/calibration_specs.json` on the validation branch). Suite reports are written to `validation_results/`. Keep those prefixes stable when adding or merging validation tools.
