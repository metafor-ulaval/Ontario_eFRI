@echo off
rem Launch the eFRI Shiny app with the portable R shipped in this folder.
rem Closing the browser tab stops the app.

setlocal
set "APP_DIR=%~dp0"
set "APP_DIR=%APP_DIR:~0,-1%"
set "RSCRIPT=%APP_DIR%\R-portable\bin\Rscript.exe"

if not exist "%RSCRIPT%" (
  echo Portable R was not found: "%RSCRIPT%"
  pause
  exit /b 1
)

rem Give R its home folder: otherwise R asks Windows for the Documents folder,
rem which can be returned with a wrong encoding when it contains accents
rem (e.g. OneDrive - Universite Laval) and makes the background computation fail
set "HOME=%USERPROFILE%"
set "R_USER=%USERPROFILE%"

rem Only use the packages shipped with the portable R
set "R_LIBS="
set "R_LIBS_USER=%APP_DIR%\R-portable\library"
set "R_LIBS_SITE="
set "EFRI_DESKTOP=1"

set "APP_DIR_R=%APP_DIR:\=/%"
"%RSCRIPT%" -e "shiny::runApp('%APP_DIR_R%', launch.browser = TRUE)"

if errorlevel 1 pause
