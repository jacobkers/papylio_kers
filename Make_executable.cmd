ECHO OFF
CALL conda activate trace_analysis
::pyinstaller GUI.py --paths trace_analysis\mapping --hidden-import trace_analysis.plugins --collect-data distributed --onedir --add-data "trace_analysis\default_configuration.yml;."
pyinstaller GUI.spec --noconfirm
call conda deactivate
::xcopy "%cd%\dist\traceAnalysisGUI" "M:\tnw\bn\cmj\Shared\Code\traceAnalysisGUI" /w /s
pause