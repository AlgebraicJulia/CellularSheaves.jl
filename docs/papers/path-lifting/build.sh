#!/bin/bash
# Build paper.pdf (and page previews with the argument "preview"). Run on a compute node:
#   srun -A fairbanksj -p hpg-dev -t 10 -c 1 --mem=2G in-sandbox bash build.sh [preview]
set -e
source /etc/profile.d/modules.sh
cd "$(dirname "$0")"
module load julia/1.12.6
julia --startup-file=no figures/make_reeb_figure.jl
module load texlive/2023
pdflatex -interaction=nonstopmode paper.tex > /dev/null || true
bibtex paper > /dev/null || true
pdflatex -interaction=nonstopmode paper.tex > /dev/null || true
pdflatex -interaction=nonstopmode paper.tex | grep -E 'undefined|^!' || true
grep -A4 '^!' paper.log || true
if [ "${1:-}" = preview ]; then
  module load ghostscript/9.22
  rm -f preview-*.png
  gs -q -sDEVICE=png16m -r70 -o preview-%02d.png paper.pdf
fi
