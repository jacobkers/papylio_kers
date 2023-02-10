ECHO OFF
CALL conda activate trace_analysis
pyinstaller GUI.py --paths trace_analysis\mapping --hidden-import trace_analysis.plugins --collect-data distributed --onefolder --add-data "trace_analysis\default_configuration.yml;."
pause