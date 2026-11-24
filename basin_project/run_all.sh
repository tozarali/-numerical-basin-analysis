#!/usr/bin/env bash
# Linux/macOS equivalent of run_all.cmd. Arguments go to basinlab.run_all.
set -euo pipefail
cd "$(dirname "$0")"
if ! conda env list | awk '{print $1}' | grep -qx basin-project; then
  conda env create -f environment.yml
fi
conda run --no-capture-output -n basin-project python -m basinlab.experiments.check_environment
conda run --no-capture-output -n basin-project python -m unittest discover -s tests
conda run --no-capture-output -n basin-project python -m basinlab.run_all "$@"
