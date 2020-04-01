ECHO OFF
git pull --ff-only
CALL conda update conda
CALL conda env update --file environment.yml --prune
pause
