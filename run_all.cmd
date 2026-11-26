@echo off
rem Double-click or run from a terminal; extra arguments go to basinlab.run_all
rem   run_all.cmd --quick    run_all.cmd --list    run_all.cmd --only regime_map
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_all.ps1" %*
set "result=%ERRORLEVEL%"
echo.
echo Exit code: %result%
pause
exit /b %result%
