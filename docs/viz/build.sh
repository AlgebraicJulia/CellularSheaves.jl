#!/bin/bash
# Build docs/viz/tropical_demo.html: bundle src/main.ts with esbuild (a standalone Go
# binary: ESBUILD=/path/to/esbuild, or `go install github.com/evanw/esbuild/cmd/esbuild@latest`)
# and inline it, with data.json and the gzipped state fields (*.bin, from
# docs/scripts/tropical_viz_export.jl), into template.html.
set -euo pipefail
cd "$(dirname "$0")"
ESBUILD=${ESBUILD:-esbuild}
"$ESBUILD" src/main.ts --bundle --minify --target=es2020 --format=iife --outfile=app.js --log-level=warning
: > blobs.html
for name in vstar vup lower16; do
    printf '<script id="blob-%s" type="application/octet-stream">' "$name" >> blobs.html
    gzip -9 -c "$name.bin" | base64 -w0 >> blobs.html
    printf '</script>\n' >> blobs.html
done
awk -v data=data.json -v app=app.js -v blobs=blobs.html '
    index($0, "/*DATA*/")  { sub(/\/\*DATA\*\/.*/, ""); printf "%s", $0; while ((getline l < data) > 0) printf "%s", l; print "</script>"; next }
    index($0, "/*BLOBS*/") { while ((getline l < blobs) > 0) print l; next }
    index($0, "/*APP*/")   { sub(/\/\*APP\*\/.*/, "");  printf "%s", $0; while ((getline l < app) > 0) print l; print "</script>"; next }
    { print }' template.html > tropical_demo.html
rm -f app.js blobs.html
echo "wrote docs/viz/tropical_demo.html ($(stat -c %s tropical_demo.html) bytes)"
