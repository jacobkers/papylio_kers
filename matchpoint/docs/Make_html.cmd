call conda activate trace_analysis
call sphinx-build -M clean . _build
call sphinx-build -M html . _build
call conda deactivate
pause