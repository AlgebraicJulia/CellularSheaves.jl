#!/bin/bash
# Build docs/viz/tropical_demo.html: bundle src/main.ts with esbuild (a standalone Go
# binary: ESBUILD=/path/to/esbuild, or `go install github.com/evanw/esbuild/cmd/esbuild@latest`)
# and inline it, with data.json (from docs/scripts/tropical_viz_export.jl), into template.html.
set -euo pipefail
cd "$(dirname "$0")"
ESBUILD=${ESBUILD:-esbuild}
"$ESBUILD" src/main.ts --bundle --minify --target=es2020 --format=iife --outfile=app.js --log-level=warning
awk -v data=data.json -v app=app.js '
    index($0, "/*DATA*/") { sub(/\/\*DATA\*\/.*/, ""); printf "%s", $0; while ((getline l < data) > 0) printf "%s", l; print "</script>"; next }
    index($0, "/*APP*/")  { sub(/\/\*APP\*\/.*/, "");  printf "%s", $0; while ((getline l < app) > 0) print l; print "</script>"; next }
    { print }' template.html > tropical_demo.html
rm -f app.js
echo "wrote docs/viz/tropical_demo.html ($(du -h tropical_demo.html | cut -f1))"
