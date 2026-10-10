#!/bin/bash
# Run Julia on the CellularSheaves dependency image built by sysimage/build.jl:
#   sysimage/julia.sh [julia arguments...]
# Uses the image, the environment it was built from (with its MPI and CUDA
# preferences) unless --project is given, and the CPU target it was built for.
# CellularSheaves, CliqueTrees and Mumblebee load on top as ordinary packages.
dir=${CS_SYSIMAGE_DIR:-$HOME/sysimage}
image=$dir/current.so
[[ -e $image ]] || { echo "sysimage/julia.sh: no image at $image; run sysimage/build.jl first" >&2; exit 1; }
export JULIA_CPU_TARGET=$(cat "$dir/cpu_target")
project=(--project="$dir/env")
for arg in "$@"; do [[ $arg == --project* ]] && project=(); done
exec julia --sysimage="$image" "${project[@]}" "$@"
