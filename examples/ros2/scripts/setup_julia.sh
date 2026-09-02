#!/usr/bin/env bash
# setup_julia.sh : prepare the container's Julia environment. Idempotent.
#
# The container does NOT reuse examples/ros2/Manifest.toml. That manifest belongs
# to the host: it pins a different Julia version and it may dev-depend on sibling
# checkouts (CliqueTrees.jl, for one) that exist beside the repo on a workstation
# and nowhere inside an image. Sharing it means the host and the container take
# turns re-resolving each other's work, and the container's resolve fails outright
# on the missing sibling.
#
# So the container gets its own environment directory on a volume, holding a copy
# of Project.toml with CellularSheaves developed from the mounted repo. The repo's
# own manifest is never touched.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
REPO="$(cd "$ROOT/../.." && pwd)"
ENV_DIR="${JULIA_ENV_DIR:-/ws/julia_env}"

mkdir -p "$ENV_DIR"
cp "$ROOT/Project.toml" "$ENV_DIR/Project.toml"

echo "[setup_julia] env=$ENV_DIR  repo=$REPO  depot=${JULIA_DEPOT_PATH:-default}"
SHEAF_REPO="$REPO" julia --project="$ENV_DIR" -e '
    using Pkg
    Pkg.develop(path = ENV["SHEAF_REPO"])
    Pkg.instantiate()
    Pkg.precompile()
'
echo "[setup_julia] ready"
