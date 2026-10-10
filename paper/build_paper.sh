#!/usr/bin/env bash
# Rebuild paper/main.pdf. Figures and tables are regenerated first from the data in the repository.
set -e
cd "$(dirname "$0")/.."
python postprocessing/acoustics/post_acoustique_combine.py
python postprocessing/aero/perf_aero.py
python paper/scripts/dipole_check.py
cd paper
pdflatex -interaction=nonstopmode main.tex >/dev/null
bibtex main >/dev/null
pdflatex -interaction=nonstopmode main.tex >/dev/null
pdflatex -interaction=nonstopmode main.tex
