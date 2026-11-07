"""File loading for the engine (the only I/O layer). Reads params/ and groups/ from a bundle root."""
import csv
import hashlib
import re
import json
from pathlib import Path

from lgo_engine import ParameterSet

ROOT = Path(__file__).resolve().parents[1]


def _csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def load_params(root=ROOT):
    root = Path(root)
    P = ParameterSet(_csv(root / "params/onsite.csv"), _csv(root / "params/pairs.csv"),
                     json.loads((root / "params/presets.json").read_text())["presets"],
                     json.loads((root / "params/channel_defaults.json").read_text())["defaults"])
    # Release 3: molecule packs (packs/index.json) and the names of the group files, so a pack whose group
    # is missing can fall back to C1 (F2.14). Packs are data beside the CONVENTIONS fixtures.
    P.packs = load_packs(root)
    P.group_names = {p.stem for p in Path(root, "groups").glob("*.json")}
    return P


def load_packs(root=ROOT):
    idx = Path(root, "packs/index.json")
    if not idx.exists():
        return []
    return [json.loads(Path(root, "packs", f"{i}.json").read_text()) for i in json.loads(idx.read_text())["packs"]]


def load_groups(root=ROOT):
    return {p.stem: json.loads(p.read_text()) for p in sorted(Path(root, "groups").glob("*.json"))}


def conventions_sha256(root=ROOT):
    return hashlib.sha256(Path(root, "CONVENTIONS.md").read_bytes()).hexdigest()


def conventions_version(root=ROOT):
    """The version and status line from CONVENTIONS.md's header, e.g. ("1.4", "not frozen")."""
    head = Path(root, "CONVENTIONS.md").read_text().split("\n", 5)
    for line in head:
        m = re.match(r"\*\*Version: ([0-9.]+)\s*[—-]\s*(.*?)\*\*", line)
        if m:
            return m.group(1), ("frozen" if "NOT FROZEN" not in m.group(2).upper() else "not frozen")
    return None, None


def load_extras(root=ROOT):
    """Reading notes and curated reference data for the app panel (not used by the calculation)."""
    root = Path(root)
    return {"notes": json.loads((root / "params/case_notes.json").read_text()),
            "curated": json.loads((root / "params/curated_checks.json").read_text())}
