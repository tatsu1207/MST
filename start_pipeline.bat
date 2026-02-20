@echo off
:: MST Web — End-to-end pipeline UI (app2.py)
:: Double-click this file to launch the full QC → DADA2 → SourceTracker → Pathogen pipeline.

title MST Pipeline

:: Verify conda is available
where conda >nul 2>&1
if errorlevel 1 (
    echo [ERR] conda not found on PATH.
    echo       Install Miniconda and restart your terminal:
    echo       https://docs.conda.io/en/latest/miniconda.html
    pause
    exit /b 1
)

:: Verify the ST environment exists
conda env list | findstr /B "ST " >nul 2>&1
if errorlevel 1 (
    echo [ERR] Conda environment "ST" not found.
    echo       Run setup_wsl.sh inside WSL to create it first.
    pause
    exit /b 1
)

echo.
echo  Starting MST pipeline app ^(app2.py^)...
echo  Open your browser at:  http://localhost:8501
echo  Press Ctrl+C here to stop the server.
echo.

conda run --no-capture-output -n ST streamlit run scripts/app2.py

pause
