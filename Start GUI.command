#!/bin/bash
cd "${0%/*}"
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate papylio
python "Start GUI.py"
conda deactivate

#python tt.py
  