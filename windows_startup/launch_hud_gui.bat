@echo off
title mtDNA Pipeline - Graphical HUD
echo =====================================================================
echo  Launching mtDNA Pipeline HUD (WSL2 Backend)
echo =====================================================================
echo.

:: Wait 2 seconds in background for Streamlit server startup, then open browser
start /min cmd /c "timeout /t 2 /nobreak >nul & start "" http://localhost:8501"

:: Run Streamlit headless within WSL2 active environment
wsl.exe -e bash -lic "cd ~/amtdna-variant-pipeline && conda activate adna_qc && streamlit run src/gui/app.py --server.headless true --server.port 8501"

pause