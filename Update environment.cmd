ECHO OFF
CALL conda activate trace_analysis
CALL conda env update --file environment.yml --prune
pause
