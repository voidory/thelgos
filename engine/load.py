"""File loading for the engine (the only I/O layer). Reads params/ and groups/ from a bundle root."""
import csv
import hashlib
import json
from pathlib import Path

from lgo_engine import ParameterSet

ROOT = Path(__file__).resolve().parents[1]


def _csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def load_params(root=ROOT):
    root = Path(root)
    return ParameterSet(_csv(root / "params/onsite.csv"), _csv(root / "params/pairs.csv"),
                        json.loads((root / "params/presets.json").read_text())["presets"],
                        json.loads((root / "params/channel_defaults.json").read_text())["defaults"])


def load_groups(root=ROOT):
    return {p.stem: json.loads(p.read_text()) for p in sorted(Path(root, "groups").glob("*.json"))}


def conventions_sha256(root=ROOT):
    return hashlib.sha256(Path(root, "CONVENTIONS.md").read_bytes()).hexdigest()


def load_extras(root=ROOT):
    """Reading notes and curated reference data for the app panel (not used by the calculation)."""
    root = Path(root)
    return {"notes": json.loads((root / "params/case_notes.json").read_text()),
            "curated": json.loads((root / "params/curated_checks.json").read_text())}
