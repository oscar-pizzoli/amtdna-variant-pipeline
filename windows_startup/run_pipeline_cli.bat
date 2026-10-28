@echo off
title mtDNA Pipeline - Terminal Runner
echo =====================================================================
echo  Executing mtDNA Pipeline via WSL2 Terminal
echo =====================================================================
echo.

wsl.exe -e bash -lic "cd ~/amtdna-variant-pipeline && conda activate adna_qc && python3 -u scripts/run_complete_pipeline.py"

echo.
echo Pipeline execution completed.
pause