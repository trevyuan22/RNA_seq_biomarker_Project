@echo off
cd /d "%~dp0"
where python >nul 2>nul && (set PY=python) || (set PY=py -3)
%PY% -c "import numpy,pandas,scipy,sklearn,matplotlib,joblib" >nul 2>nul || (echo Installing the required packages. This only happens the first time and takes a few minutes. & %PY% -m pip install -r requirements.txt)
%PY% -c "import pydeseq2" >nul 2>nul || (echo Installing optional extras. & %PY% -m pip install -r requirements-optional.txt)
echo Starting the app. A page will open in your web browser.
%PY% app.py
pause
