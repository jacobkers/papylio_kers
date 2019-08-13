IF NOT EXIST %~dp0\.git\ (
	git clone --depth=1 https://gitlab.tudelft.nl/ivoseverins/traceanalysis.git
	.\traceAnalysis\Libraries.cmd
) ELSE (
	ECHO 'Already a git repository present'
)
pause