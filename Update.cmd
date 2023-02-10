ECHO OFF
CALL conda update conda
CALL conda activate trace_analysis
git pull --ff-only
CALL conda env update --file environment.yml --prune
CALL conda deactivate
pause
