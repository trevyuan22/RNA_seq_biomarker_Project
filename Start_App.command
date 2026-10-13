#!/bin/bash
cd "$(dirname "$0")"
python3 -c "import numpy,pandas,scipy,sklearn,matplotlib,joblib" 2>/dev/null || { echo "Installing the required packages. This only happens the first time."; python3 -m pip install -r requirements.txt; }
python3 -c "import pydeseq2" 2>/dev/null || python3 -m pip install -r requirements-optional.txt || true
echo "Starting the app. A page will open in your web browser."
python3 app.py
