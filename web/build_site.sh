#!/usr/bin/env bash
# Assemble the static site for GitHub Pages into the directory given (default: _site).
# Only what the browser app fetches is published: the web app, the engine sources, params/, groups/, packs/
# and CONVENTIONS.md (read for its version line and hash). Evidence snapshots, checker scripts and
# tests stay in the repository but are not served as part of the site.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-_site}"
rm -rf "$OUT"
mkdir -p "$OUT/web" "$OUT/engine/lgo_engine" "$OUT/params" "$OUT/groups" "$OUT/packs" "$OUT/matcher"
cp web/index.html web/app.js web/worker.js web/style.css web/proposal.css "$OUT/web/"
cp engine/load.py "$OUT/engine/"
cp engine/lgo_engine/*.py "$OUT/engine/lgo_engine/"
cp params/*.csv params/*.json "$OUT/params/"
cp groups/*.json "$OUT/groups/"
cp packs/*.json "$OUT/packs/"
cp matcher/library.json "$OUT/matcher/"
cp CONVENTIONS.md "$OUT/"
# Jekyll would drop files starting with "_" (the engine's __init__.py); serve files as they are.
touch "$OUT/.nojekyll"
cat > "$OUT/index.html" <<'HTML'
<!doctype html>
<meta charset="utf-8">
<title>LGO Method Explorer</title>
<meta http-equiv="refresh" content="0; url=web/">
<link rel="canonical" href="web/">
<p><a href="web/">LGO Method Explorer</a></p>
HTML
echo "site assembled in $OUT"
