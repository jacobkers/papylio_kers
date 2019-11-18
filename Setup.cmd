ECHO OFF
IF EXIST "%~dp0%.git\" (
	ECHO 'Already a git repository present'
) ELSE (
	git clone https://gitlab.tudelft.nl/ivoseverins/traceanalysis.git
	CALL conda env create -f traceAnalysis\environment.yml
)
pause
