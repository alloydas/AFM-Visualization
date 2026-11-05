"""Shared paths for mesh generators under tools/."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MESHES = ROOT / "meshes"
TIPS = MESHES / "tips"
SURFACES = MESHES / "surfaces"
SCANS = MESHES / "scans"
RAW = ROOT / "data" / "raw"


def mesh_path(subdir: Path, name: str) -> str:
    subdir.mkdir(parents=True, exist_ok=True)
    return str(subdir / name)
