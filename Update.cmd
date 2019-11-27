ECHO OFF
git pull --ff-only
CALL conda env update --file environment.yml --prune
pause
