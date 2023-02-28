call conda activate trace_analysis
call make clean
call make html
call conda deactivate
xcopy "%cd%\_build\html" "M:\tnw\bn\cmj\Shared\Code\Documentation" /w /s
pause