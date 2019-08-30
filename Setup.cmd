ECHO OFF
IF EXIST "%~dp0%.git\" (
	ECHO 'Already a git repository present'
) ELSE (
	git clone --depth=1 https://gitlab.tudelft.nl/ivoseverins/traceanalysis.git
	CALL conda create --yes --name trace_analysis python=3.7
	CALL traceanalysis\Libraries.cmd
)
pause
