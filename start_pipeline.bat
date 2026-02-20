@echo off
:: MST Web — End-to-end pipeline UI (app2.py)
:: Double-click this file to launch the full QC → DADA2 → SourceTracker → Pathogen pipeline.

title MST Pipeline

:: Verify WSL is available
where wsl >nul 2>&1
if errorlevel 1 (
    echo [ERR] WSL not found. Please enable Windows Subsystem for Linux.
    pause
    exit /b 1
)

echo.
echo  Starting MST pipeline app ^(app2.py^)...
echo  Open your browser at:  http://localhost:8501
echo  Press Ctrl+C here to stop the server.
echo.

wsl bash -c "cd /mnt/d/github/MST_test && /home/unnot/miniforge3/bin/conda run --no-capture-output -n ST streamlit run scripts/app2.py"

pause
