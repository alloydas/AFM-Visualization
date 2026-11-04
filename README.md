# AFM Visualization

Interactive AFM tip–sample simulators (2D / 3D), mesh generators, and forward-model validation against real Bruker SPM data.

## Layout

| Path | Role |
|------|------|
| `index.html`, `afm-2d.html`, `afm-3d.html`, `README.html` | Deployed web app (Vercel routes these at the repo root) |
| `vercel.json` | Static hosting routes |
| `afm_forward_model.py`, `spm_io.py` | Forward-model physics + SPM I/O (**stay at root** for imports) |
| `validate_forward_experiment.py`, `validate_tip_geometry.py` | Validation runners (**stay at root** — paths resolve from repo root) |
| `tools/` | Mesh generators (`make_*.py`) only |
| `meshes/{tips,surfaces,scans}/` | Generated Wavefront `.obj` meshes |
| `samples/` | Calibration assets, datasheets, SPM folders (**do not rename** — pipeline contract) |
| `validation_results/` | Validation outputs (**do not rename** — pipeline contract; gitignored except `.gitkeep`) |
| `paper/` | Forward-model validation paper + figures |
| `data/raw/` | Large raw SPM archives (e.g. TGXYZ02); gitignored |
| `docs/` | Notes and drafts |

## Validation (reproducibility)

Paths in `samples/calibration_specs.json` and the validators are rooted at the **repository root**:

```bash
python validate_forward_experiment.py --case vgrp15m_snl10
# reads:  samples/.../*.spm  + samples/calibration_specs.json
# writes: validation_results/vgrp15m_snl10/
```

Do not move `samples/` or `validation_results/`. Keep the validators next to `afm_forward_model.py` / `spm_io.py` so `ROOT = Path(__file__).parent` stays correct.

## Mesh generators

```bash
python tools/make_snl10_obj.py
python tools/make_vgrp_obj.py
python tools/make_spm_obj.py   # reads data/raw/TGXYZ02/
```

Outputs land under `meshes/`.
