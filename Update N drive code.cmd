call conda activate trace_analysis
rmdir "N:\tnw\BN\Shared\Ivo\traceAnalysis" /s
mkdir "N:\tnw\BN\Shared\Ivo\traceAnalysis"
git clone https://gitlab.tudelft.nl/ivoseverins/traceanalysis.git "N:\tnw\BN\Shared\Ivo\traceAnalysis" --branch develop
pause