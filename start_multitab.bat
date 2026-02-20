@echo off
:: MST Web — Multi-tab UI (app.py)
:: Double-click this file to launch all tools in a single tabbed interface.

title MST Multi-tab

:: Verify WSL is available
where wsl >nul 2>&1
if errorlevel 1 (
    echo [ERR] WSL not found. Please enable Windows Subsystem for Linux.
    pause
    exit /b 1
)

echo.
echo  Starting MST multi-tab app ^(app.py^)...
echo  Open your browser at:  http://localhost:8501
echo  Press Ctrl+C here to stop the server.
echo.

wsl bash -c "cd /mnt/d/github/MST_test && /home/unnot/miniforge3/bin/conda run --no-capture-output -n ST streamlit run scripts/app.py"

pause
