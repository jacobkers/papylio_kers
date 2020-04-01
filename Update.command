#!/bin/bash
cd "${0%/*}"
git pull --ff-only
conda update conda
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda env update --file environment.yml --prune
