"""`make add-molecule PACK=packs/<id>.json`: the molecule-pack gate (spec section 14, Release 3 Milestone A).

Schema check, then the engine-side checks in lgo_engine.packs.gate: formula, parameters (owner decision
R3-1), electron count, declared symmetry or C1 fallback, R4 count, ledger totals, the pairing schedule and
morph invariants. The matcher fixture check arrives with the matcher (Milestone B). Exits non-zero, naming
the failing check, if any fails. Nothing is written."""
import argparse
import json
import sys
from pathlib import Path

from jsonschema import Draft202012Validator, ValidationError

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "engine")]
from load import load_params, load_groups  # noqa: E402
from lgo_engine import packs  # noqa: E402


def main(argv=None):
    a = argparse.ArgumentParser()
    a.add_argument("pack", type=Path)
    path = a.parse_args(argv).pack
    pack = json.loads(path.read_text())
    schema = json.loads((ROOT / "packs/pack.schema.json").read_text())
    try:
        Draft202012Validator(schema).validate(pack)
        results = [("schema", True, "")]
    except ValidationError as e:
        print(f"FAIL schema: {e.message} at {'/'.join(map(str, e.absolute_path))}")
        return 1
    results += packs.gate(load_params(ROOT), load_groups(ROOT), pack)
    results.append(("matcher fixture", True, "not yet available: the matcher arrives in Release 3 Milestone B"))
    listed = pack["id"] in json.loads((ROOT / "packs/index.json").read_text())["packs"]
    for name, ok, detail in results:
        print(f"{'PASS' if ok else 'FAIL'} {name}{': ' + detail if detail else ''}")
    if not listed:
        print(f"NOTE {pack['id']} is not listed in packs/index.json, so the app will not show it")
    failed = [n for n, ok, _ in results if not ok]
    print(f"{pack['id']}: {'passes' if not failed else 'fails: ' + '; '.join(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
