echo (Start script) - %DATE%_%TIME% >> process.log

:iteration

echo (Test) - %DATE%_%TIME% >> process.log

start /WAIT  matlab -nodesktop -nojvm -nosplash -minimize -logfile out.txt  -wait -r runOptimization;exit;
call :GetUnixTime UNIX_TIME
findstr /v "To get started, type one of these: helpwin, helpdesk, or demo." out.txt > optimizeLuminariesNResult_%UNIX_TIME%.txt
ren optimizeLuminariesNResult.mat optimizeLuminariesNResult_%UNIX_TIME%.mat

echo (Test Fin) - %DATE%_%TIME% >> process.log

GOTO iteration

echo (End script) - %DATE%_%TIME% >> process.log

REM method to get time as unix
:GetUnixTime
setlocal enableextensions
for /f %%x in ('wmic path win32_utctime get /format:list ^| findstr "="') do (
    set %%x)
set /a z=(14-100%Month%%%100)/12, y=10000%Year%%%10000-z
set /a ut=y*365+y/4-y/100+y/400+(153*(100%Month%%%100+12*z-3)+2)/5+Day-719469
set /a ut=ut*86400+100%Hour%%%100*3600+100%Minute%%%100*60+100%Second%%%100
endlocal & set "%1=%ut%" & goto :EOF